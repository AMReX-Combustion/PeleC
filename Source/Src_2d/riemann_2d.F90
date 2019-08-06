module riemann_module

  use amrex_fort_module
  use amrex_constants_module
  use riemann_util_module

  use meth_params_module, only : NQ, QVAR, NVAR, QRHO, QU, QV, QW, &
                                 QPRES, QREINT, QFS, &
                                 QFX, URHO, UMX, UMY, UEDEN, UEINT, &
                                 GDPRES, GDGAME, &
                                 NGDNV, small_dens, small_pres, small_temp, &
                                 cg_maxiter, cg_tol, cg_blend, &
                                 npassive, upass_map, qpass_map, &
                                 riemann_solver, ppm_temp_fix, hybrid_riemann, &
                                 allow_negative_energy

  implicit none

  private

  public cmpflx, shock, riemanncg

  real (amrex_real), parameter :: smallu = 1.e-12_amrex_real

contains

! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine cmpflx(qm,qp,qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
                    flx,flx_l1,flx_l2,flx_h1,flx_h2, &
                    qint, qg_l1,qg_l2,qg_h1,qg_h2, &
                    gamc,csml,c,qd_l1,qd_l2,qd_h1,qd_h2, &
                    bcMask, bcMask_l1, bcMask_l2, bcMask_h1, bcMask_h2, &
                    shk,s_l1,s_l2,s_h1,s_h2, &
                    idir,ilo,ihi,jlo,jhi,domlo,domhi)

    use eos_type_module
    use eos_module, only: eos_re, eos_rt, mine

    implicit none
    integer, intent(in) :: qpd_l1,qpd_l2,qpd_h1,qpd_h2
    integer, intent(in) :: flx_l1,flx_l2,flx_h1,flx_h2
    integer, intent(in) :: qg_l1,qg_l2,qg_h1,qg_h2
    integer, intent(in) :: qd_l1,qd_l2,qd_h1,qd_h2
    integer, intent(in) :: bcMask_l1, bcMask_l2, bcMask_h1, bcMask_h2

    integer, intent(in) :: s_l1,s_l2,s_h1,s_h2
    integer, intent(in) :: idir,ilo,ihi,jlo,jhi
    integer, intent(in) :: domlo(2),domhi(2)

    integer, intent(inout) :: bcMask(bcMask_l1:bcMask_h1,bcMask_l2:bcMask_h2)
    
    double precision, intent(inout) :: qint(qg_l1:qg_h1,qg_l2:qg_h2,NGDNV)

    double precision, intent(inout) ::  qm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
    double precision, intent(inout) ::  qp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
    double precision, intent(inout) :: flx(flx_l1:flx_h1,flx_l2:flx_h2,NVAR)

    double precision, intent(in) :: gamc(qd_l1:qd_h1,qd_l2:qd_h2)
    double precision, intent(in) ::    c(qd_l1:qd_h1,qd_l2:qd_h2)
    double precision, intent(in) :: csml(qd_l1:qd_h1,qd_l2:qd_h2)
    double precision, intent(in) ::  shk( s_l1: s_h1, s_l2: s_h2)

    ! Local variables
    integer i, j

    double precision, allocatable :: smallc(:,:), cavg(:,:)
    double precision, allocatable :: gamcm(:,:), gamcp(:,:)
    integer :: imin, imax, jmin, jmax
    integer :: is_shock
    double precision :: cl, cr
    type (eos_t) :: eos_state

    allocate ( smallc(ilo-1:ihi+1,jlo-1:jhi+1) )
    allocate (   cavg(ilo-1:ihi+1,jlo-1:jhi+1) )
    allocate (  gamcm(ilo-1:ihi+1,jlo-1:jhi+1) )
    allocate (  gamcp(ilo-1:ihi+1,jlo-1:jhi+1) )
    if (idir == 1) then
       do j = jlo, jhi
          do i = ilo, ihi+1
             smallc(i,j) = max( csml(i,j), csml(i-1,j) )
             cavg(i,j) = HALF*( c(i,j) + c(i-1,j) )
             gamcm(i,j) = gamc(i-1,j)
             gamcp(i,j) = gamc(i,j)
          enddo
       enddo

    else
       do j = jlo, jhi+1
          do i = ilo, ihi
             smallc(i,j) = max( csml(i,j), csml(i,j-1) )
             cavg(i,j) = HALF*( c(i,j) + c(i,j-1) )
             gamcm(i,j) = gamc(i,j-1)
             gamcp(i,j) = gamc(i,j)
          enddo
       enddo
    endif

    if (ppm_temp_fix == 2) then
       ! recompute the thermodynamics on the interface to make it
       ! all consistent -- THIS PROBABLY DOESN"T WORK WITH RADIATION

       ! we want to take the edge states of rho, p, and X, and get
       ! new values for gamc and (rho e) on the edges that are
       ! thermodynamically consistent.
       if (idir == 1) then
          imin = ilo
          imax = ihi+1
          jmin = jlo
          jmax = jhi
       else
          imin = ilo
          imax = ihi
          jmin = jlo
          jmax = jhi+1
       endif

       call build(eos_state)

       do j = jmin, jmax
          do i = imin, imax

             ! this is an initial guess for iterations, since we
             ! can't be certain that temp is on interfaces
             eos_state%T = 1000.0d0

             ! minus state
             eos_state % rho      = qm(i,j,QRHO)
             eos_state % p        = qm(i,j,QPRES)
             eos_state % e        = qm(i,j,QREINT)/qm(i,j,QRHO)
             eos_state % massfrac = qm(i,j,QFS:QFS-1+nspecies)
             eos_state % aux      = qm(i,j,QFX:QFX-1+naux)

             ! Protect against negative energies

             if (allow_negative_energy .eq. 0 .and. eos_state % e < mine) then
                eos_state % T = small_temp
                call eos_rt(eos_state)
             else
                call eos_re(eos_state)
             endif

             qm(i,j,QREINT) = qm(i,j,QRHO)*eos_state%e
             qm(i,j,QPRES) = eos_state%p
             gamcm(i,j) = eos_state%gam1


             ! plus state
             eos_state % rho      = qp(i,j,QRHO)
             eos_state % p        = qp(i,j,QPRES)
             eos_state % e        = qp(i,j,QREINT)/qp(i,j,QRHO)
             eos_state % massfrac = qp(i,j,QFS:QFS-1+nspecies)
             eos_state % aux      = qp(i,j,QFX:QFX-1+naux)

             ! Protect against negative energies

             if (allow_negative_energy .eq. 0 .and. eos_state % e < mine) then
                eos_state % T = small_temp
                call eos_rt(eos_state)
             else
                call eos_re(eos_state)
             endif

             qp(i,j,QREINT) = qp(i,j,QRHO)*eos_state%e
             qp(i,j,QPRES) = eos_state%p
             gamcp(i,j) = eos_state%gam1

          enddo
       enddo

       call destroy(eos_state)

    endif

    ! Solve Riemann problem (godunov state passed back, but only (u,p) saved)
    if (riemann_solver == 0) then
       ! Colella, Glaz, & Ferguson solver
       call riemannus(qm, qp, qpd_l1, qpd_l2, qpd_h1, qpd_h2, &
                      gamcm, gamcp, cavg, smallc, ilo-1, jlo-1, ihi+1, jhi+1, &
                      flx, flx_l1, flx_l2, flx_h1, flx_h2, &
                      qint, qg_l1, qg_l2, qg_h1, qg_h2, &
                      bcMask, bcMask_l1, bcMask_l2, bcMask_h1, bcMask_h2, &
                      idir, ilo, ihi, jlo, jhi, domlo, domhi)

    elseif (riemann_solver == 1) then
       ! Colella & Glaz solver
       call riemanncg(qm, qp, qpd_l1, qpd_l2, qpd_h1, qpd_h2, &
                      gamcm, gamcp, cavg, smallc, ilo-1, jlo-1, ihi+1, jhi+1, &
                      flx, flx_l1, flx_l2, flx_h1, flx_h2, &
                      qint, qg_l1, qg_l2, qg_h1, qg_h2, &
                      bcMask, bcMask_l1, bcMask_l2, bcMask_h1, bcMask_h2, &
                      idir, ilo, ihi, jlo, jhi, domlo, domhi)

    elseif (riemann_solver == 2) then
       ! HLLC
       call HLLC(qm, qp, qpd_l1, qpd_l2, qpd_h1, qpd_h2, &
                 gamcm, gamcp, cavg, smallc, ilo-1, jlo-1, ihi+1, jhi+1, &
                 flx, flx_l1, flx_l2, flx_h1, flx_h2, &
                 qint, qg_l1, qg_l2, qg_h1, qg_h2, &
                 bcMask, bcMask_l1, bcMask_l2, bcMask_h1, bcMask_h2, &
                 idir, ilo, ihi, jlo, jhi, domlo, domhi)

    elseif (riemann_solver == 3) then
       ! Marc solver (for Fuego gases)
       call riemannmd(qm, qp, qpd_l1, qpd_l2, qpd_h1, qpd_h2, &
                      gamcm, gamcp, cavg, smallc, ilo-1, jlo-1, ihi+1, jhi+1, &
                      flx, flx_l1, flx_l2, flx_h1, flx_h2, &
                      qint, qg_l1, qg_l2, qg_h1, qg_h2, &
                      bcMask, bcMask_l1, bcMask_l2, bcMask_h1, bcMask_h2, &
                      idir, ilo, ihi, jlo, jhi, domlo, domhi)

    else
       call bl_error("ERROR: invalid value of riemann_solver")
    endif

    if (hybrid_riemann == 1) then
       ! correct the fluxes using an HLL scheme if we are in a shock
       ! and doing the hybrid approach
       if (idir == 1) then
          imin = ilo
          imax = ihi+1
          jmin = jlo
          jmax = jhi
       else
          imin = ilo
          imax = ihi
          jmin = jlo
          jmax = jhi+1
       endif

       do j = jmin, jmax
          do i = imin, imax

             if (idir == 1) then
                is_shock = shk(i-1,j) + shk(i,j)
             else
                is_shock = shk(i,j-1) + shk(i,j)
             endif

             if (is_shock >= 1) then

                if (idir == 1) then
                   cl = c(i-1,j)
                   cr = c(i,j)
                else
                   cl = c(i,j-1)
                   cr = c(i,j)
                endif

                call HLL(qm(i,j,:), qp(i,j,:), cl, cr, &
                         idir, 2, flx(i,j,:))

             endif

          enddo
       enddo

    endif

    deallocate(smallc,cavg,gamcm,gamcp)
  end subroutine cmpflx


  subroutine shock(q,qd_l1,qd_l2,qd_h1,qd_h2, &
                   shk,s_l1,s_l2,s_h1,s_h2, &
                   ilo1,ilo2,ihi1,ihi2,dx,dy)

    use prob_params_module, only : coord_type

    implicit none
    integer, intent(in) :: qd_l1, qd_l2, qd_h1, qd_h2
    integer, intent(in) :: s_l1, s_l2, s_h1, s_h2
    integer, intent(in) :: ilo1, ilo2, ihi1, ihi2
    double precision, intent(in) :: dx, dy
    double precision, intent(in) :: q(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
    double precision, intent(inout) :: shk(s_l1:s_h1,s_l2:s_h2)

    integer :: i, j

    double precision :: divU
    double precision :: px_pre, px_post, py_pre, py_post
    double precision :: e_x, e_y, d
    double precision :: p_pre, p_post, pjump

    double precision :: rc, rm, rp

    double precision, parameter :: small = 1.d-10
    double precision, parameter :: eps = 0.33d0

    ! This is a basic multi-dimensional shock detection algorithm.
    ! This implementation follows Flash, which in turn follows
    ! AMRA and a Woodward (1995) (supposedly -- couldn't locate that).
    !
    ! The spirit of this follows the shock detection in Colella &
    ! Woodward (1984)

    do j = ilo2-1, ihi2+1
       do i = ilo1-1, ihi1+1

          ! construct div{U}
          if (coord_type == 0) then
             divU = HALF*(q(i+1,j,QU) - q(i-1,j,QU))/dx + &
                    HALF*(q(i,j+1,QV) - q(i,j-1,QV))/dy
          else if (coord_type == 1) then
             ! r-z
             rc = dble(i + HALF)*dx
             rm = dble(i - 1 + HALF)*dx
             rp = dble(i + 1 + HALF)*dx

             divU = HALF*(rp*q(i+1,j,QU) - rm*q(i-1,j,QU))/(rc*dx) + &
                    HALF*(q(i,j+1,QV) - q(i,j-1,QV))/dy
          else
             call bl_error("ERROR: invalid coord_type in shock")
          endif
             
          ! find the pre- and post-shock pressures in each direction
          if (q(i+1,j,QPRES) - q(i-1,j,QPRES) < ZERO) then
             px_pre  = q(i+1,j,QPRES)
             px_post = q(i-1,j,QPRES)
          else
             px_pre  = q(i-1,j,QPRES)
             px_post = q(i+1,j,QPRES)
          endif

          if (q(i,j+1,QPRES) - q(i,j-1,QPRES) < ZERO) then
             py_pre  = q(i,j+1,QPRES)
             py_post = q(i,j-1,QPRES)
          else
             py_pre  = q(i,j-1,QPRES)
             py_post = q(i,j+1,QPRES)
          endif

          ! use compression to create unit vectors for the shock direction
          e_x = (q(i+1,j,QU) - q(i-1,j,QU))**2
          e_y = (q(i,j+1,QV) - q(i,j-1,QV))**2
          d = ONE/(e_x + e_y + small)

          e_x = e_x*d
          e_y = e_y*d

          ! project the pressures onto the shock direction
          p_pre  = e_x*px_pre + e_y*py_pre
          p_post = e_x*px_post + e_y*py_post

          ! test for compression + pressure jump to flag a shock
          if (p_pre == ZERO) then
             ! this can arise if e_x = e_y = 0 (U = 0)
             pjump = ZERO
          else
             pjump = eps - (p_post - p_pre)/p_pre
          endif

          if (pjump < ZERO .and. divU < ZERO) then
             shk(i,j) = ONE
          else
             shk(i,j) = ZERO
          endif

       enddo
    enddo

  end subroutine shock


! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine riemanncg(ql,qr,qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
                       gamcl,gamcr,cav,smallc,gd_l1,gd_l2,gd_h1,gd_h2, &
                       uflx,uflx_l1,uflx_l2,uflx_h1,uflx_h2, &
                       qint, qg_l1,qg_l2,qg_h1,qg_h2, &
                       bcMask, bcMask_l1, bcMask_l2, bcMask_h1, bcMask_h2, &
                       idir,ilo1,ihi1,ilo2,ihi2,domlo,domhi)

    ! this implements the approximate Riemann solver of Colella & Glaz (1985)

    use amrex_fort_module
    use network, only : nspecies, naux
    use eos_type_module
    use eos_module
    use prob_params_module, only : coord_type

    implicit none
    double precision, parameter:: small = 1.d-8
    double precision, parameter :: small_u = 1.d-10

    integer :: qpd_l1,qpd_l2,qpd_h1,qpd_h2
    integer :: gd_l1,gd_l2,gd_h1,gd_h2
    integer :: uflx_l1,uflx_l2,uflx_h1,uflx_h2
    integer :: qg_l1,qg_l2,qg_h1,qg_h2
    integer :: idir,ilo1,ihi1,ilo2,ihi2
    integer :: qd_l1, qd_h1,qd_l2, qd_h2
    integer :: bcMask_l1, bcMask_l2, bcMask_h1, bcMask_h2
    integer :: domlo(2),domhi(2)

    double precision :: ql(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
    double precision :: qr(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
    double precision ::  gamcl(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision ::  gamcr(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision ::    cav(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision :: smallc(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision :: uflx(uflx_l1:uflx_h1,uflx_l2:uflx_h2,NVAR)
    double precision :: qint(qg_l1:qg_h1,qg_l2:qg_h2,NGDNV)
    integer :: bcMask(bcMask_l1:bcMask_h1,bcMask_l2:bcMask_h2)

    integer :: i,j,ilo,jlo,ihi,jhi, ipassive
    integer :: n, nqp

    double precision :: rgdnv,vgdnv,wgdnv,ustar,gamgdnv
    double precision :: rl, ul, vl, v2l, pl, rel
    double precision :: rr, ur, vr, v2r, pr, rer
    double precision :: wl, wr, rhoetot
    double precision :: rstar, cstar, pstar
    double precision :: ro, uo, po, co, gamco
    double precision :: sgnm, spin, spout, ushock, frac
    double precision :: wsmall, csmall,qavg

    double precision :: gcl, gcr
    double precision :: clsq, clsql, clsqr, wlsq, wosq, wrsq, wo
    double precision :: zl, zr
    double precision :: denom, dpditer, dpjmp
    double precision :: gamc_bar, game_bar
    double precision :: gamel, gamer, gameo, gamstar, gmin, gmax, gdot

    integer :: iter, iter_max
    double precision :: tol
    double precision :: err

    logical :: converged

    double precision :: pstar_old
    double precision :: taul, taur, tauo
    double precision :: ustar_r, ustar_l, ustar_r_old, ustar_l_old
    double precision :: pstar_lo, pstar_hi

    double precision, parameter :: weakwv = 1.d-3

    double precision, allocatable :: pstar_hist(:), pstar_hist_extra(:)

    type (eos_t) :: eos_state

    integer :: iu, iv1, iv2
    integer :: idx, idy

    call build(eos_state)

    if (cg_blend .eq. 2 .and. cg_maxiter < 5) then

       call bl_error("Error: need cg_maxiter >= 5 to do a bisection search on secant iteration failure.")

    endif

    tol = cg_tol
    iter_max = cg_maxiter

    !  set min/max based on normal direction
    if (idir == 1) then
       ilo = ilo1
       ihi = ihi1 + 1
       jlo = ilo2
       jhi = ihi2

       iu = QU
       iv1 = QV
       iv2 = QW
    else
       ilo = ilo1
       ihi = ihi1
       jlo = ilo2
       jhi = ihi2+1

       iu = QV
       iv1 = QU
       iv2 = QW
    endif

    allocate (pstar_hist(iter_max))
    allocate (pstar_hist_extra(iter_max))

    do j = jlo, jhi
       do i = ilo, ihi

          ! left state
          rl = max(ql(i,j,QRHO),small_dens)

          ! pick left velocities based on direction
          ul = ql(i,j,iu)
          vl = ql(i,j,iv1)
          v2l = ql(i,j,iv2)

          pl  = ql(i,j,QPRES )
          rel = ql(i,j,QREINT)
          gcl = gamcl(i,j)

          ! sometimes we come in here with negative energy or pressure
          ! note: reset both in either case, to remain thermo
          ! consistent
          if (rel <= ZERO .or. pl <= small_pres) then
             print *, "WARNING: (rho e)_l < 0 or pl < small_pres in Riemann: ", rel, pl, small_pres
             eos_state % T        = small_temp
             eos_state % rho      = rl
             eos_state % massfrac = ql(i,j,QFS:QFS-1+nspecies)
             eos_state % aux      = ql(i,j,QFX:QFX-1+naux)

             call eos_rt(eos_state)

             rel = rl*eos_state%e
             pl  = eos_state%p
             gcl = eos_state%gam1
          endif

          ! right state
          rr = max(qr(i,j,QRHO),small_dens)

          ! pick right velocities based on direction
          if (idir == 1) then
             ur = qr(i,j,QU)
             vr = qr(i,j,QV)
             v2r = qr(i,j,QW)
          else
             ur = qr(i,j,QV)
             vr = qr(i,j,QU)
             v2r = qr(i,j,QW)
          endif

          pr  = qr(i,j,QPRES)
          rer = qr(i,j,QREINT)
          gcr = gamcr(i,j)

          if (rer <= ZERO .or. pr <= small_pres) then
             print *, "WARNING: (rho e)_r < 0 or pr < small_pres in Riemann: ", rer, pr, small_pres
             eos_state % T        = small_temp
             eos_state % rho      = rr
             eos_state % massfrac = qr(i,j,QFS:QFS-1+nspecies)
             eos_state % aux      = qr(i,j,QFX:QFX-1+naux)

             call eos_rt(eos_state)

             rer = rr*eos_state%e
             pr  = eos_state%p
             gcr = eos_state%gam1
          endif

          ! common quantities
          taul = ONE/rl
          taur = ONE/rr

          ! lagrangian sound speeds
          clsql = gcl*pl*rl
          clsqr = gcr*pr*rr


          ! Note: in the original Colella & Glaz paper, they predicted
          ! gamma_e to the interfaces using a special (non-hyperbolic)
          ! evolution equation.  In PeleC, we instead bring (rho e)
          ! to the edges, so we construct the necessary gamma_e here from
          ! what we have on the interfaces.
          gamel = pl/rel + ONE
          gamer = pr/rer + ONE

          ! these should consider a wider average of the cell-centered
          ! gammas
          gmin = min(gamel, gamer, ONE, FOUR3RD)
          gmax = max(gamel, gamer, TWO, FIVE3RD)

          game_bar = HALF*(gamel + gamer)
          gamc_bar = HALF*(gcl + gcr)

          gdot = TWO*(ONE - game_bar/gamc_bar)*(game_bar - ONE)

          csmall = smallc(i,j)
          wsmall = small_dens*csmall
          wl = max(wsmall,sqrt(abs(clsql)))
          wr = max(wsmall,sqrt(abs(clsqr)))

          ! make an initial guess for pstar -- this is a two-shock
          ! approximation
          !pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))/(wl + wr)
          pstar = pl + ( (pr - pl) - wr*(ur - ul) )*wl/(wl+wr)
          pstar = max(pstar,small_pres)

          ! get the shock speeds -- this computes W_s from CG Eq. 34
          call wsqge(pl,taul,gamel,gdot,  &
                     gamstar,pstar,wlsq,clsql,gmin,gmax)

          call wsqge(pr,taur,gamer,gdot,  &
                     gamstar,pstar,wrsq,clsqr,gmin,gmax)

          pstar_old = pstar

          wl = sqrt(wlsq)
          wr = sqrt(wrsq)

          ! R-H jump conditions give ustar across each wave -- these should
          ! be equal when we are done iterating
          ustar_l = ul - (pstar-pl)/wl
          ustar_r = ur + (pstar-pr)/wr

          ! revise our pstar guess
          !pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))/(wl + wr)
          pstar = pl + ( (pr - pl) - wr*(ur - ul) )*wl/(wl+wr)
          pstar = max(pstar,small_pres)

          ! secant iteration
          converged = .false.
          iter = 1
          do while ((iter <= iter_max .and. .not. converged) .or. iter <= 2)

             call wsqge(pl,taul,gamel,gdot,  &
                        gamstar,pstar,wlsq,clsql,gmin,gmax)

             call wsqge(pr,taur,gamer,gdot,  &
                        gamstar,pstar,wrsq,clsqr,gmin,gmax)

             wl = ONE / sqrt(wlsq)
             wr = ONE / sqrt(wrsq)

             ustar_l_old = ustar_l
             ustar_r_old = ustar_r

             ! note that wl, wr here are already inverses
             ustar_l = ul - (pstar - pl)*wl
             ustar_r = ur + (pstar - pr)*wr

             dpditer = pstar - pstar_old

             zl = ustar_l - ustar_l_old
             !if (zp-weakwv*cav(i,j) <= ZERO) then
             !   zp = dpditer*wl
             !endif

             zr = ustar_r - ustar_r_old
             !if (zm-weakwv*cav(i,j) <= ZERO) then
             !   zm = dpditer*wr
             !endif
             
             ! the new pstar is found via CG Eq. 18

             !denom = zp + zm
             denom = zl - zr
             denom = (ustar_l - ustar_r) - (ustar_l_old - ustar_r_old)

             pstar_old = pstar

             if (abs(denom) > small_u .and. abs(dpditer) > small_pres) then
                pstar = pstar - (ustar_l - ustar_r)*dpditer/denom
             endif

             pstar = max(pstar, small_pres)

             err = abs(pstar - pstar_old)
             if (err < tol*pstar) converged = .true.

             pstar_hist(iter) = pstar

             iter = iter + 1

          enddo

          ! If we failed to converge using the secant iteration, we can either
          ! stop here; or, revert to the original two-shock estimate for pstar;
          ! or do a bisection root find using the bounds established by the most
          ! recent iterations.

          if (.not. converged) then

             if (cg_blend .eq. 0) then

                print *, 'pstar history: '
                do iter = 1, iter_max
                   print *, iter, pstar_hist(iter)
                enddo

                print *, ' '
                print *, 'left state  (r,u,p,re,gc): ', rl, ul, pl, rel, gcl
                print *, 'right state (r,u,p,re,gc): ', rr, ur, pr, rer, gcr
                print *, 'cav, smallc:',  cav(i,j), csmall
                call bl_error("ERROR: non-convergence in the Riemann solver")

             else if (cg_blend .eq. 1) then

                pstar = pl + ( (pr - pl) - wr*(ur - ul) )*wl/(wl+wr)

             else if (cg_blend .eq. 2) then

                ! first try to find a reasonable bounds 
                pstar_lo = minval(pstar_hist(iter_max-5:iter_max))
                pstar_hi = maxval(pstar_hist(iter_max-5:iter_max))

                call pstar_bisection(pstar_lo, pstar_hi, &
                                     ul, pl, taul, gamel, clsql, &
                                     ur, pr, taur, gamer, clsqr, &
                                     gdot, gmin, gmax, &
                                     pstar, gamstar, converged, pstar_hist_extra)

                if (.not. converged) then
                   ! abort -- doesn't seem solvable
                   print *, 'pstar history: '
                   do iter = 1, iter_max
                      print *, iter, pstar_hist(iter)
                   enddo

                   do iter = 1, iter_max
                      print *, iter+iter_max, pstar_hist_extra(iter)
                   enddo

                   print *, ' '
                   print *, 'left state  (r,u,p,re,gc): ', rl, ul, pl, rel, gcl
                   print *, 'right state (r,u,p,re,gc): ', rr, ur, pr, rer, gcr
                   print *, 'cav, smallc:',  cav(i,j), csmall
                   call bl_error("ERROR: non-convergence in the Riemann solver")

                endif

             else

                call bl_error("ERROR: unrecognized cg_blend option.")

             endif

          endif


          ! we converged!  construct the single ustar for the region
          ! between the left and right waves, using the updated wave speeds
          ustar_r = ur-(pr-pstar)*wr  ! careful -- here wl, wr are 1/W
          ustar_l = ul+(pl-pstar)*wl

          ustar = HALF*(ustar_l + ustar_r)

          ! for symmetry preservation, if ustar is really small, then we
          ! set it to zero
          if (abs(ustar) < smallu*HALF*(abs(ul) + abs(ur))) then
             ustar = ZERO
          endif

          ! sample the solution -- here we look first at the direction
          ! that the contact is moving.  This tells us if we need to
          ! worry about the L/L* states or the R*/R states.
          if (ustar .gt. ZERO) then
             ro = rl
             uo = ul
             po = pl
             tauo = taul
             !reo = rel
             gamco = gcl
             gameo = gamel

          else if (ustar .lt. ZERO) then
             ro = rr
             uo = ur
             po = pr
             tauo = taur
             !reo = rer
             gamco = gcr
             gameo = gamer
          else
             ro = HALF*(rl+rr)
             uo = HALF*(ul+ur)
             po = HALF*(pl+pr)
             tauo = HALF*(taul+taur)
             !reo = HALF*(rel+rer)
             gamco = HALF*(gcl+gcr)
             gameo = HALF*(gamel + gamer)
          endif

          ! use tau = 1/rho as the independent variable here
          ro = max(small_dens,ONE/tauo)
          tauo = ONE/ro

          co = sqrt(abs(gamco*po/ro))
          co = max(csmall,co)
          clsq = (co*ro)**2

          ! now that we know which state (left or right) we need to worry
          ! about, get the value of gamstar and wosq across the wave we
          ! are dealing with.
          call wsqge(po,tauo,gameo,gdot,   &
                     gamstar,pstar,wosq,clsq,gmin,gmax)

          sgnm = sign(ONE,ustar)

          wo = sqrt(wosq)
          dpjmp = pstar - po

          ! is this max really necessary?
          !rstar=max(ONE-ro*dpjmp/wosq, (gameo-ONE)/(gameo+ONE))
          rstar=ONE-ro*dpjmp/wosq
          rstar=ro/rstar
          rstar = max(small_dens,rstar)

          !entho = (reo/ro + po/ro)/co**2
          !estar = reo + (pstar - po)*entho

          cstar = sqrt(abs(gamco*pstar/rstar))
          cstar = max(cstar,csmall)


          spout = co - sgnm*uo
          spin = cstar - sgnm*ustar

          !ushock = HALF*(spin + spout)
          ushock = wo/ro - sgnm*uo

          if (pstar-po .ge. ZERO) then
             spin = ushock
             spout = ushock
          endif
          ! if (spout-spin .eq. ZERO) then
          !    scr = small*cav(i,j)
          ! else
          !    scr = spout-spin
          ! endif
          ! frac = (ONE + (spout + spin)/scr)*HALF
          ! frac = max(ZERO,min(ONE,frac))

          frac = HALF*(ONE + (spin + spout)/max(spout-spin,spin+spout, small*cav(i,j)))

          ! the transverse velocity states only depend on the
          ! direction that the contact moves
          if (ustar .gt. ZERO) then
             vgdnv = vl
             wgdnv = v2l
          else if (ustar .lt. ZERO) then
             vgdnv = vr
             wgdnv = v2r
          else
             vgdnv = HALF*(vl+vr)
             wgdnv = HALF*(v2l+v2r)
          endif

          ! linearly interpolate between the star and normal state -- this covers the
          ! case where we are inside the rarefaction fan.
          rgdnv = frac*rstar + (ONE - frac)*ro
          qint(i,j,iu) = frac*ustar + (ONE - frac)*uo
          qint(i,j,GDPRES) = frac*pstar + (ONE - frac)*po
          gamgdnv =  frac*gamstar + (ONE-frac)*gameo

          ! now handle the cases where instead we are fully in the
          ! star or fully in the original (l/r) state
          if (spout .lt. ZERO) then
             rgdnv = ro
             qint(i,j,iu) = uo
             qint(i,j,GDPRES) = po
             gamgdnv = gameo
          endif
          if (spin .ge. ZERO) then
             rgdnv = rstar
             qint(i,j,iu) = ustar
             qint(i,j,GDPRES) = pstar
             gamgdnv = gamstar
          endif

          qint(i,j,GDGAME) = gamgdnv

          qint(i,j,GDPRES) = max(qint(i,j,GDPRES),small_pres)

          ! Enforce that fluxes through a symmetry plane or wall are hard zero.
          ! Here the NSCBC info about if we have a wall or not is contained in the ghost-cell
          qint(i,j,iu) = bc_test(idir, i, j, &
                                 bcMask, bcMask_l1, bcMask_l2, bcMask_h1, bcMask_h2, &
                                 domlo, domhi) * qint(i,j,iu)

          ! Compute fluxes, order as conserved state (not q)
          uflx(i,j,URHO) = rgdnv*qint(i,j,iu)

          ! note: for axisymmetric geometries, we do not include the
          ! pressure in the r-direction, since div{F} + grad{p} cannot
          ! be written in a flux difference form
          if (idir == 1) then
             uflx(i,j,UMX) = uflx(i,j,URHO)*qint(i,j,iu)
             uflx(i,j,UMY) = uflx(i,j,URHO)*vgdnv
             if (coord_type == 0) then
                uflx(i,j,UMX) = uflx(i,j,UMX) + qint(i,j,GDPRES)
             endif
          else
             uflx(i,j,UMX) = uflx(i,j,URHO)*vgdnv
             uflx(i,j,UMY) = uflx(i,j,URHO)*qint(i,j,iu) + qint(i,j,GDPRES)
          endif

          ! compute the total energy from the internal, p/(gamma - 1), and the kinetic
          rhoetot = qint(i,j,GDPRES)/(gamgdnv - ONE) + &
               HALF*rgdnv*(qint(i,j,iu)**2 + vgdnv**2 + wgdnv**2)

          uflx(i,j,UEDEN) = qint(i,j,iu)*(rhoetot + qint(i,j,GDPRES))
          uflx(i,j,UEINT) = qint(i,j,iu)*qint(i,j,GDPRES)/(gamgdnv - ONE)

          ! advected quantities -- only the contact matters
          ! note: this includes the z-velocity flux
          do ipassive = 1, npassive
             n  = upass_map(ipassive)
             nqp = qpass_map(ipassive)

             if (ustar .gt. ZERO) then
                uflx(i,j,n) = uflx(i,j,URHO)*ql(i,j,nqp)
             else if (ustar .lt. ZERO) then
                uflx(i,j,n) = uflx(i,j,URHO)*qr(i,j,nqp)
             else
                qavg = HALF * (ql(i,j,nqp) + qr(i,j,nqp))
                uflx(i,j,n) = uflx(i,j,URHO)*qavg
             endif
          enddo

       enddo
    enddo

    deallocate(pstar_hist_extra)
    deallocate(pstar_hist)

  end subroutine riemanncg

! :::
! ::: ------------------------------------------------------------------
! :::


  subroutine riemannus(ql, qr, qpd_l1, qpd_l2, qpd_h1, qpd_h2, &
                       gamcl, gamcr, cav, smallc, gd_l1, gd_l2, gd_h1, gd_h2, &
                       uflx, uflx_l1, uflx_l2, uflx_h1, uflx_h2, &
                       qint, qg_l1, qg_l2, qg_h1, qg_h2, &
                       bcMask, bcMask_l1, bcMask_l2, bcMask_h1, bcMask_h2, &
                       idir, ilo1, ihi1, ilo2, ihi2, domlo, domhi)

    use prob_params_module, only : coord_type

    implicit none
    double precision, parameter:: small = 1.d-8
    integer :: qpd_l1, qpd_l2, qpd_h1, qpd_h2
    integer :: gd_l1, gd_l2, gd_h1, gd_h2
    integer :: uflx_l1, uflx_l2, uflx_h1, uflx_h2
    integer :: qg_l1, qg_l2, qg_h1, qg_h2
    integer :: idir, ilo1, ihi1, ilo2, ihi2
    integer :: qd_l1, qd_h1,qd_l2, qd_h2
    integer :: bcMask_l1, bcMask_l2, bcMask_h1, bcMask_h2
    integer :: domlo(2),domhi(2)

    double precision :: ql(qpd_l1:qpd_h1,qpd_l2:qpd_h2,NQ)
    double precision :: qr(qpd_l1:qpd_h1,qpd_l2:qpd_h2,NQ)

    double precision :: gamcl(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision :: gamcr(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision :: cav(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision :: smallc(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision :: uflx(uflx_l1:uflx_h1,uflx_l2:uflx_h2,NVAR)
    double precision :: qint(qg_l1:qg_h1,qg_l2:qg_h2,NGDNV)
    integer :: bcMask(bcMask_l1:bcMask_h1,bcMask_l2:bcMask_h2)
    integer :: ilo,ihi,jlo,jhi
    integer :: n, nqp
    integer :: i, j, ipassive

    double precision :: rgd, vgd, wgd, regd, ustar
    double precision :: rl, ul, vl, v2l, pl, rel
    double precision :: rr, ur, vr, v2r, pr, rer
    double precision :: wl, wr, rhoetot, scr
    double precision :: rstar, cstar, estar, pstar
    double precision :: ro, uo, po, reo, co, gamco, entho, drho
    double precision :: sgnm, spin, spout, ushock, frac
    double precision :: wsmall, csmall,qavg

    integer :: iu, iv1, iv2
    integer :: idx, idy

    !************************************************************
    !  set min/max based on normal direction
    if (idir == 1) then
       ilo = ilo1
       ihi = ihi1 + 1
       jlo = ilo2
       jhi = ihi2

       iu = QU
       iv1 = QV
       iv2 = QW
    else
       ilo = ilo1
       ihi = ihi1
       jlo = ilo2
       jhi = ihi2+1

       iu = QV
       iv1 = QU
       iv2 = QW
    endif

    do j = jlo, jhi
       do i = ilo, ihi

          rl = ql(i,j,QRHO)

          !  pick left velocities based on direction
          ul = ql(i,j,iu)
          vl = ql(i,j,iv1)
          v2l = ql(i,j,iv2)

          pl = ql(i,j,QPRES)
          rel = ql(i,j,QREINT)
          rr = qr(i,j,QRHO)

          !  pick right velocities based on direction
          if (idir == 1) then
             ur = qr(i,j,QU)
             vr = qr(i,j,QV)
             v2r = qr(i,j,QW)
          else
             ur = qr(i,j,QV)
             vr = qr(i,j,QU)
             v2r = qr(i,j,QW)
          endif

          pr = qr(i,j,QPRES)
          rer = qr(i,j,QREINT)

          csmall = smallc(i,j)
          wsmall = small_dens*csmall
          wl = max(wsmall,sqrt(abs(gamcl(i,j)*pl*rl)))
          wr = max(wsmall,sqrt(abs(gamcr(i,j)*pr*rr)))

          pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))/(wl + wr)
          ustar = ((wl*ul + wr*ur) + (pl - pr))/(wl + wr)

          pstar = max(pstar,small_pres)
          ! for symmetry preservation, if ustar is really small, then we
          ! set it to zero
          if (abs(ustar) < smallu*HALF*(abs(ul) + abs(ur))) then
             ustar = ZERO
          endif

          if (ustar .gt. ZERO) then
             ro = rl
             uo = ul
             po = pl

             reo = rel
             gamco = gamcl(i,j)

          else if (ustar .lt. ZERO) then
             ro = rr
             uo = ur
             po = pr

             reo = rer
             gamco = gamcr(i,j)

          else
             ro = HALF*(rl+rr)
             uo = HALF*(ul+ur)
             po = HALF*(pl+pr)

             reo = HALF*(rel+rer)
             gamco = HALF*(gamcl(i,j)+gamcr(i,j))

          endif

          ro = max(small_dens,ro)

          co = sqrt(abs(gamco*po/ro))
          co = max(csmall,co)

          drho = (pstar - po)/co**2
          rstar = ro + drho
          rstar = max(small_dens,rstar)

          entho = (reo/ro + po/ro)/co**2
          estar = reo + (pstar - po)*entho
          cstar = sqrt(abs(gamco*pstar/rstar))
          cstar = max(cstar,csmall)

          sgnm = sign(ONE,ustar)
          spout = co - sgnm*uo
          spin = cstar - sgnm*ustar

          ushock = HALF*(spin + spout)

          if (pstar-po >= ZERO) then
             spin = ushock
             spout = ushock
          endif

          if (spout-spin == ZERO) then
             scr = small*cav(i,j)
          else
             scr = spout-spin
          endif

          frac = (ONE + (spout + spin)/scr)*HALF
          frac = max(ZERO,min(ONE,frac))

          if (ustar .gt. ZERO) then
             vgd = vl
             wgd = v2l
          else if (ustar .lt. ZERO) then
             vgd = vr
             wgd = v2r
          else
             vgd = HALF*(vl+vr)
             wgd = HALF*(v2l+v2r)
          endif

          rgd = frac*rstar + (ONE - frac)*ro

          qint(i,j,iu) = frac*ustar + (ONE - frac)*uo
          qint(i,j,iv1) = vgd
          qint(i,j,iv2) = wgd

          qint(i,j,GDPRES) = frac*pstar + (ONE - frac)*po
          regd = frac*estar + (ONE - frac)*reo
          if (spout < ZERO) then
             rgd = ro
             qint(i,j,iu) = uo
             qint(i,j,GDPRES) = po
             regd = reo
          endif

          if (spin >= ZERO) then
             rgd = rstar
             qint(i,j,iu) = ustar
             qint(i,j,GDPRES) = pstar
             regd = estar
          endif

          ! not sure what this should be for radiation?
          qint(i,j,GDGAME) = qint(i,j,GDPRES)/regd + ONE
          
          ! enforce that the fluxes through a symmetry plane or wall are zero
          ! Here the NSCBC info about if we have a wall or not is contained in the ghost-cell
          !if (idir == 1) write(*,*) 'DEBUG RIEMAMNN',i,j,bcMask(0,j)
          !if (idir == 2) write(*,*) 'DEBUG RIEMAMNN',i,j,bcMask(i,domhi(2)+1)
          qint(i,j,iu) = bc_test(idir, i, j, &
                                 bcMask, bcMask_l1, bcMask_l2, bcMask_h1, bcMask_h2, &
                                 domlo, domhi) * qint(i,j,iu)

          ! Compute fluxes, order as conserved state (not q)
          uflx(i,j,URHO) = rgd*qint(i,j,iu)

          ! note: for axisymmetric geometries, we do not include the
          ! pressure in the r-direction, since div{F} + grad{p} cannot
          ! be written in a flux difference form
          if (idir == 1) then
             uflx(i,j,UMX) = uflx(i,j,URHO)*qint(i,j,iu)
             uflx(i,j,UMY) = uflx(i,j,URHO)*vgd
             if (coord_type == 0) then
                uflx(i,j,UMX) = uflx(i,j,UMX) + qint(i,j,GDPRES)
             endif
          else
             uflx(i,j,UMX) = uflx(i,j,URHO)*vgd
             uflx(i,j,UMY) = uflx(i,j,URHO)*qint(i,j,iu) + qint(i,j,GDPRES)
          endif

          rhoetot = regd + HALF*rgd*(qint(i,j,iu)**2 + vgd**2 + wgd**2)

          uflx(i,j,UEDEN) = qint(i,j,iu)*(rhoetot + qint(i,j,GDPRES))
          uflx(i,j,UEINT) = qint(i,j,iu)*regd
          do ipassive = 1, npassive
             n  = upass_map(ipassive)
             nqp = qpass_map(ipassive)

             if (ustar .gt. ZERO) then
                uflx(i,j,n) = uflx(i,j,URHO)*ql(i,j,nqp)
             else if (ustar .lt. ZERO) then
                uflx(i,j,n) = uflx(i,j,URHO)*qr(i,j,nqp)
             else
                qavg = HALF * (ql(i,j,nqp) + qr(i,j,nqp))
                uflx(i,j,n) = uflx(i,j,URHO)*qavg
             endif
          enddo

       enddo
    enddo
  end subroutine riemannus


! :::
! ::: ------------------------------------------------------------------
! :::


  subroutine riemannmd(ql, qr, qpd_l1, qpd_l2, qpd_h1, qpd_h2, &
                       gamcl, gamcr, cav, smallc, gd_l1, gd_l2, gd_h1, gd_h2, &
                       uflx, uflx_l1, uflx_l2, uflx_h1, uflx_h2, &
                       qint, qg_l1, qg_l2, qg_h1, qg_h2, &
                       bcMask, bcMask_l1, bcMask_l2, bcMask_h1, bcMask_h2, &
                       idir, ilo1, ihi1, ilo2, ihi2, domlo, domhi)

    use prob_params_module, only : coord_type
    use network, only : nspecies

    use eos_module

    implicit none
    integer :: qpd_l1, qpd_l2, qpd_h1, qpd_h2
    integer :: gd_l1, gd_l2, gd_h1, gd_h2
    integer :: uflx_l1, uflx_l2, uflx_h1, uflx_h2
    integer :: qg_l1, qg_l2, qg_h1, qg_h2
    integer :: idir, ilo1, ihi1, ilo2, ihi2
    integer :: qd_l1, qd_h1,qd_l2, qd_h2
    integer :: bcMask_l1, bcMask_l2, bcMask_h1, bcMask_h2
    integer :: domlo(2),domhi(2)

    double precision :: ql(qpd_l1:qpd_h1,qpd_l2:qpd_h2,NQ)
    double precision :: qr(qpd_l1:qpd_h1,qpd_l2:qpd_h2,NQ)

    double precision :: gamcl(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision :: gamcr(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision :: cav(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision :: smallc(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision :: uflx(uflx_l1:uflx_h1,uflx_l2:uflx_h2,NVAR)
    double precision :: qint(qg_l1:qg_h1,qg_l2:qg_h2,NGDNV)
    double precision :: qavg, uflx_w_dummy
    integer :: bcMask(bcMask_l1:bcMask_h1,bcMask_l2:bcMask_h2)
    integer :: ilo,ihi,jlo,jhi
    integer :: n, nqp
    integer :: i, j, ipassive
    integer :: idx, idy


    ! Copy these into scratch here because we may modify them
    ! for outflow hack
    double precision :: ul, vl, v2l, rel, ur, vr, v2r, rer 

    double precision :: rgd, regd, ustar
    integer :: bc_test_mask
  

    integer :: iu, iv1, iv2

    type(eos_t) :: eos_state, gdnv_state

    call build(eos_state)
    call build(gdnv_state)

    !************************************************************
    !  set min/max based on normal direction
    if (idir == 1) then
       ilo = ilo1
       ihi = ihi1 + 1
       jlo = ilo2
       jhi = ihi2

       iu = QU
       iv1 = QV
       iv2 = QW
    else
       ilo = ilo1
       ihi = ihi1
       jlo = ilo2
       jhi = ihi2+1

       iu = QV
       iv1 = QU
       iv2 = QW
    endif

    do j = jlo, jhi
       do i = ilo, ihi
       
        ! Here the NSCBC info about if we have a wall or not is contained in the ghost-cell
        bc_test_mask = bc_test(idir, i, j, &
                               bcMask, bcMask_l1, bcMask_l2, bcMask_h1, bcMask_h2, &
                               domlo, domhi)
                               
        !TODO: consider transposing ql, qr on pass into this routine so that passing the species doesn't make a
        ! strided copy into a temporary

         ul = ql(i,j,iu)
         vl = ql(i,j,iv1)
        v2l = ql(i,j,iv2)
        rel = ql(i,j,QREINT)

         ur = qr(i,j,iu)
         vr = qr(i,j,iv1)
        v2r = qr(i,j,iv2)
        rer = qr(i,j,QREINT)
        
        call outflow_hack(ul,ur,vl,vr,v2l,v2r,rel,rer,&
                          bcMask, bcMask_l1, bcMask_l2, bcMask_h1, bcMask_h2, &
                          idir, i, j, domlo, domhi)
!

        call riemann_md_singlepoint( &
          ql(i,j,QRHO), ul, vl, v2l, ql(i,j,QPRES), rel, ql(i,j,QFS:QFS+nspecies-1), gamcl(i,j), &
          qr(i,j,QRHO), ur, vr, v2r, qr(i,j,QPRES), rer, qr(i,j,QFS:QFS+nspecies-1), gamcr(i,j), &
          qint(i,j,iu), qint(i,j,iv1), qint(i,j,iv2), qint(i,j,GDPRES),qint(i,j,GDGAME), &
          regd, rgd, ustar, &
          eos_state, gdnv_state, nspecies, &
          uflx(i,j,URHO), uflx(i,j,UMX), uflx(i,j,UMY), uflx_w_dummy, uflx(i,j,UEDEN), uflx(i,j,UEINT), &
          idir, coord_type, bc_test_mask, smallc(i,j), cav(i,j) )



  
          do ipassive = 1, npassive
             n  = upass_map(ipassive)
             nqp = qpass_map(ipassive)

             if (ustar .gt. ZERO) then
                uflx(i,j,n) = uflx(i,j,URHO)*ql(i,j,nqp)
             else if (ustar .lt. ZERO) then
                uflx(i,j,n) = uflx(i,j,URHO)*qr(i,j,nqp)
             else
                qavg = HALF * (ql(i,j,nqp) + qr(i,j,nqp))
                uflx(i,j,n) = uflx(i,j,URHO)*qavg
             endif
          enddo

       enddo
    enddo

    call destroy(eos_state)
    call destroy(gdnv_state)

  end subroutine riemannmd


! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine HLLC(ql, qr, qpd_l1, qpd_l2, qpd_h1, qpd_h2, &
                  gamcl, gamcr, cav, smallc, gd_l1, gd_l2, gd_h1, gd_h2, &
                  uflx, uflx_l1, uflx_l2, uflx_h1, uflx_h2, &
                  qint, qg_l1, qg_l2, qg_h1, qg_h2, &
                  bcMask, bcMask_l1, bcMask_l2, bcMask_h1, bcMask_h2, &
                  idir, ilo1, ihi1, ilo2, ihi2, domlo, domhi)

    ! this is an implementation of the HLLC solver described in Toro's
    ! book.  it uses the simplest estimate of the wave speeds, since
    ! those should work for a general EOS.  We also initially do the
    ! CGF Riemann construction to get pstar and ustar, since we'll need
    ! to know the pressure and velocity on the interface for the grad p
    ! term in momentum and for an internal energy update

    implicit none

    double precision, parameter:: small = 1.d-8
    integer :: qpd_l1, qpd_l2, qpd_h1, qpd_h2
    integer :: gd_l1, gd_l2, gd_h1, gd_h2
    integer :: uflx_l1, uflx_l2, uflx_h1, uflx_h2
    integer :: qg_l1, qg_l2, qg_h1, qg_h2
    integer :: idir, ilo1, ihi1, ilo2, ihi2
    integer :: qd_l1, qd_h1,qd_l2, qd_h2
    integer :: bcMask_l1, bcMask_l2, bcMask_h1, bcMask_h2
    integer :: domlo(2),domhi(2)

    double precision :: ql(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
    double precision :: qr(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
    double precision :: gamcl(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision :: gamcr(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision :: cav(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision :: smallc(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision :: uflx(uflx_l1:uflx_h1,uflx_l2:uflx_h2,NVAR)
    double precision :: qint(qg_l1:qg_h1,qg_l2:qg_h2,NGDNV)
    integer :: bcMask(bcMask_l1:bcMask_h1,bcMask_l2:bcMask_h2)

    integer :: ilo,ihi,jlo,jhi
    integer :: i, j
    integer :: bnd_fac
    
    !double precision :: regd
    double precision :: ustar
    double precision :: rl, ul, pl, rel
    double precision :: rr, ur, pr, rer
    double precision :: wl, wr, scr
    double precision :: rstar, cstar, pstar
    double precision :: ro, uo, po, co, gamco
    double precision :: sgnm, spin, spout, ushock, frac
    double precision :: wsmall, csmall

    double precision :: U_hllc_state(nvar), U_state(nvar), F_state(nvar)
    double precision :: S_l, S_r, S_c

    integer :: iu, iv1, iv2
    integer :: idx, idy

    !  set min/max based on normal direction
    if (idir == 1) then
       ilo = ilo1
       ihi = ihi1 + 1
       jlo = ilo2
       jhi = ihi2

       iu = QU
       iv1 = QV
       iv2 = QW
    else
       ilo = ilo1
       ihi = ihi1
       jlo = ilo2
       jhi = ihi2+1

       iu = QV
       iv1 = QU
       iv2 = QW
    endif

    do j = jlo, jhi
       do i = ilo, ihi

          rl = ql(i,j,QRHO)

          ! pick left velocities based on direction
          ! ul is always normal to the interface
          if (idir == 1) then
             ul  = ql(i,j,QU)
          else
             ul  = ql(i,j,QV)
          endif

          pl = ql(i,j,QPRES)
          rel = ql(i,j,QREINT)

          rr = qr(i,j,QRHO)

          ! pick right velocities based on direction
          ! ur is always normal to the interface
          if (idir == 1) then
             ur  = qr(i,j,QU)
          else
             ur  = qr(i,j,QV)
          endif

          pr = qr(i,j,QPRES)
          rer = qr(i,j,QREINT)

          ! now we essentially do the CGF solver to get p and u on the
          ! interface, but we won't use these in any flux construction.
          csmall = smallc(i,j)
          wsmall = small_dens*csmall
          wl = max(wsmall,sqrt(abs(gamcl(i,j)*pl*rl)))
          wr = max(wsmall,sqrt(abs(gamcr(i,j)*pr*rr)))

          pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))/(wl + wr)
          ustar = ((wl*ul + wr*ur) + (pl - pr))/(wl + wr)

          pstar = max(pstar,small_pres)

          ! for symmetry preservation, if ustar is really small, then we
          ! set it to zero
          if (abs(ustar) < smallu*HALF*(abs(ul) + abs(ur))) then
             ustar = ZERO
          endif

          if (ustar > ZERO) then
             ro = rl
             uo = ul
             po = pl
             gamco = gamcl(i,j)

          else if (ustar < ZERO) then
             ro = rr
             uo = ur
             po = pr
             gamco = gamcr(i,j)

          else
             ro = HALF*(rl+rr)
             uo = HALF*(ul+ur)
             po = HALF*(pl+pr)
             gamco = HALF*(gamcl(i,j)+gamcr(i,j))
          endif

          ro = max(small_dens,ro)

          co = sqrt(abs(gamco*po/ro))
          co = max(csmall,co)

          rstar = ro + (pstar - po)/co**2
          rstar = max(small_dens,rstar)

          cstar = sqrt(abs(gamco*pstar/rstar))
          cstar = max(cstar,csmall)

          sgnm = sign(ONE,ustar)
          spout = co - sgnm*uo
          spin = cstar - sgnm*ustar

          ushock = HALF*(spin + spout)

          if (pstar-po > ZERO) then
             spin = ushock
             spout = ushock
          endif

          if (spout-spin == ZERO) then
             scr = small*cav(i,j)
          else
             scr = spout-spin
          endif
          frac = (ONE + (spout + spin)/scr)*HALF
          frac = max(ZERO,min(ONE,frac))

          qint(i,j,iu) = frac*ustar + (ONE - frac)*uo
          qint(i,j,GDPRES) = frac*pstar + (ONE - frac)*po

          ! TODO
          !gegdnv(i,j) = pgdnv(i,j)/regd + ONE

          ! now we do the HLLC construction

          ! Here the NSCBC info about if we have a wall or not is contained in the ghost-cell
          bnd_fac = bc_test(idir, i, j, &
                            bcMask, bcMask_l1, bcMask_l2, bcMask_h1, bcMask_h2, &
                            domlo, domhi)
          
          ! use the simplest estimates of the wave speeds
          S_l = min(ul - sqrt(gamcl(i,j)*pl/rl), ur - sqrt(gamcr(i,j)*pr/rr))
          S_r = max(ul + sqrt(gamcl(i,j)*pl/rl), ur + sqrt(gamcr(i,j)*pr/rr))

          ! estimate of the contact speed -- this is Toro Eq. 10.8
          S_c = (pr - pl + rl*ul*(S_l - ul) - rr*ur*(S_r - ur))/ &
             (rl*(S_l - ul) - rr*(S_r - ur))

          if (S_r <= ZERO) then
             ! R region
             call cons_state(qr(i,j,:), U_state)
             call compute_flux(idir, 2, bnd_fac, U_state, pr, F_state)

          else if (S_r > ZERO .and. S_c <= ZERO) then
             ! R* region
             call cons_state(qr(i,j,:), U_state)
             call compute_flux(idir, 2, bnd_fac, U_state, pr, F_state)

             call HLLC_state(idir, S_r, S_c, qr(i,j,:), U_hllc_state)

             ! correct the flux
             F_state(:) = F_state(:) + S_r*(U_hllc_state(:) - U_state(:))

          else if (S_c > ZERO .and. S_l < ZERO) then
             ! L* region
             call cons_state(ql(i,j,:), U_state)
             call compute_flux(idir, 2, bnd_fac, U_state, pl, F_state)

             call HLLC_state(idir, S_l, S_c, ql(i,j,:), U_hllc_state)

             ! correct the flux
             F_state(:) = F_state(:) + S_l*(U_hllc_state(:) - U_state(:))

          else
             ! L region
             call cons_state(ql(i,j,:), U_state)
             call compute_flux(idir, 2, bnd_fac, U_state, pl, F_state)

          endif


          ! Enforce that fluxes through a symmetry plane or wall are hard zero.
          ! and store the fluxes
          uflx(i,j,:) = F_state(:)

       enddo
    enddo
  end subroutine HLLC

end module riemann_module
