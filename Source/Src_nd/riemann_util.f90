module riemann_util_module

  use amrex_fort_module
  use amrex_constants_module

  implicit none

contains

  pure function bc_test(idir, i, j, &
                        bcMask, bcMask_l1, bcMask_l2, bcMask_h1, bcMask_h2, &
                        domlo, domhi) result (f)

    use prob_params_module, only : physbc_lo, physbc_hi, Symmetry, SlipWall, NoSlipWall
    
    implicit none
    integer, intent(in) :: idir, i, j, domlo(*), domhi(*)
    integer, intent(in) :: bcMask_l1, bcMask_l2, bcMask_h1, bcMask_h2
    integer, intent(in) :: bcMask(bcMask_l1:bcMask_h1,bcMask_l2:bcMask_h2)
    integer :: f

    ! Enforce that fluxes through a symmetry plane or wall are hard zero.
    f = 1
      
    if (idir == 1) then
      
      if (i == domlo(1) .and. &
              (bcMask(i,j) == Symmetry .or. &
              bcMask(i,j) == SlipWall .or. &
              bcMask(i,j) == NoSlipWall )) then
            f = 0
      endif
      
      if (i == domhi(1)+1 .and. &
              (bcMask(i,j) == Symmetry .or. &
              bcMask(i,j) == SlipWall .or. &
              bcMask(i,j) == NoSlipWall )) then
            f = 0
      endif
    end if
      
    if (idir == 2) then
            
      if (j == domlo(2) .and. &
              (bcMask(i,j) == Symmetry .or. &
              bcMask(i,j) == SlipWall .or. &
              bcMask(i,j) == NoSlipWall )) then
            f = 0
      endif
      
      if (j == domhi(2)+1 .and. &
              (bcMask(i,j) == Symmetry .or. &
              bcMask(i,j) == SlipWall .or. &
              bcMask(i,j) == NoSlipWall )) then
            f = 0
      end if
    endif
  
  end function bc_test
  
  pure function bc_test_3d(idir, i, j, k, &
                        bcMask, bcMask_l1, bcMask_l2, bcMask_l3, bcMask_h1, bcMask_h2, bcMask_h3, &
                        domlo, domhi) result (f)

    use prob_params_module, only : physbc_lo, physbc_hi, Symmetry, SlipWall, NoSlipWall
    
    implicit none
    integer, intent(in) :: idir, i, j, k, domlo(*), domhi(*)
    integer, intent(in) :: bcMask_l1, bcMask_l2, bcMask_l3, bcMask_h1, bcMask_h2, bcMask_h3
    integer, intent(in) :: bcMask(bcMask_l1:bcMask_h1,bcMask_l2:bcMask_h2,bcMask_l3:bcMask_h3)
    integer :: f

    ! Enforce that fluxes through a symmetry plane or wall are hard zero.
    f = 1
     
    if (idir == 1) then
      
      if (i == domlo(1) .and. &
              (bcMask(i,j,k) == SlipWall .or. &
               bcMask(i,j,k) == Symmetry .or. &
               bcMask(i,j,k) == NoSlipWall  )) then
            f = 0
      endif
      
      if (i == domhi(1)+1 .and. &
              (bcMask(i,j,k) == SlipWall .or. &
               bcMask(i,j,k) == Symmetry .or. &
               bcMask(i,j,k) == NoSlipWall  )) then
            f = 0
      endif
    end if
      
    if (idir == 2) then
            
      if (j == domlo(2) .and. &
              (bcMask(i,j,k) == SlipWall .or. &
               bcMask(i,j,k) == Symmetry .or. &
               bcMask(i,j,k) == NoSlipWall  )) then
            f = 0
      endif
      
      if (j == domhi(2)+1 .and. &
              (bcMask(i,j,k) == SlipWall .or. &
               bcMask(i,j,k) == Symmetry .or. &
               bcMask(i,j,k) == NoSlipWall  )) then
            f = 0
      end if
    endif
      
    if (idir == 3) then
            
      if (k == domlo(3) .and. &
              (bcMask(i,j,k) == SlipWall .or. &
               bcMask(i,j,k) == Symmetry .or. &
               bcMask(i,j,k) == NoSlipWall  )) then
            f = 0
      endif
      
      if (k == domhi(3)+1 .and. &
              (bcMask(i,j,k) == SlipWall .or. &
               bcMask(i,j,k) == Symmetry .or. &
               bcMask(i,j,k) == NoSlipWall  )) then
            f = 0
      end if
    endif
        
  end function bc_test_3d

  pure subroutine outflow_hack(ul,ur,vl,vr,v2l,v2r,rel,rer,&
                               bcMask, bcMask_l1, bcMask_l2, bcMask_h1, bcMask_h2, &
                               idir, i, j, domlo, domhi)

    use prob_params_module, only : physbc_lo, physbc_hi, Outflow

    implicit none
    double precision, intent(inout) :: ul,ur,vl,vr,v2l,v2r,rel,rer
    integer, intent(in) :: idir, i, j, domlo(*), domhi(*)
    integer, intent(in) :: bcMask_l1, bcMask_l2, bcMask_h1, bcMask_h2
    integer, intent(in) :: bcMask(bcMask_l1:bcMask_h1,bcMask_l2:bcMask_h2)
    integer :: idx

    if (idir == 1) then
       idx = i
    else
       idx = j
    endif

    if (idx == domlo(idir) .and. bcMask(i,j) == Outflow) then
       ul = ur
       vl = vr
       v2l = v2r
       rel = rer
    endif

    if (idx == domhi(idir)+1 .and. bcMask(i,j) == Outflow) then
       ur = ul
       vr = vl
       v2r = v2l
       rer = rel
    endif

  end subroutine outflow_hack

  pure subroutine wsqge(p,v,gam,gdot,gstar,pstar,wsq,csq,gmin,gmax)

    ! compute the lagrangian wave speeds.

    implicit none
    double precision, intent(in) :: p,v,gam,gdot,pstar,csq,gmin,gmax
    double precision, intent(out) :: wsq, gstar

    double precision, parameter :: smlp1 = 1.d-10
    double precision, parameter :: small = 1.d-7

    double precision :: alpha, beta

    ! First predict a value of game across the shock

    ! CG Eq. 31
    gstar = (pstar-p)*gdot/(pstar+p) + gam
    gstar = max(gmin, min(gmax, gstar))

    ! Now use that predicted value of game with the R-H jump conditions
    ! to compute the wave speed.

    ! this is CG Eq. 34
    alpha = pstar - (gstar-ONE)*p/(gam-ONE)
    if (alpha == ZERO) alpha = smlp1*(pstar + p)

    beta = pstar + HALF*(gstar-ONE)*(pstar+p)

    wsq = (pstar-p)*beta/(v*alpha)

    if (abs(pstar - p) < smlp1*(pstar + p)) then
       wsq = csq
    endif
    wsq = max(wsq, (HALF*(gam-ONE)/gam)*csq)

    return
  end subroutine wsqge


  pure subroutine pstar_bisection(pstar_lo, pstar_hi, &
       ul, pl, taul, gamel, clsql, &
       ur, pr, taur, gamer, clsqr, &
       gdot, gmin, gmax, &
       pstar, gamstar, converged, pstar_hist_extra)

    ! we want to zero                                                                     
    ! f(p*) = u*_l(p*) - u*_r(p*)                                                         
    ! we'll do bisection                                                                  

    use meth_params_module, only : cg_maxiter, cg_tol

    implicit none
    double precision, intent(inout) :: pstar_lo, pstar_hi
    double precision, intent(in) :: ul, pl, taul, gamel, clsql
    double precision, intent(in) :: ur, pr, taur, gamer, clsqr
    double precision, intent(in) :: gdot, gmin, gmax
    double precision, intent(out) :: pstar, gamstar
    logical, intent(out) :: converged
    double precision, intent(out) :: pstar_hist_extra(:)

    double precision :: pstar_c, ustar_l, ustar_r, f_lo, f_hi, f_c
    double precision :: wl, wr, wlsq, wrsq

    integer :: iter


    ! lo bounds
    call wsqge(pl, taul, gamel, gdot,  &
         gamstar, pstar_lo, wlsq, clsql, gmin, gmax)

    call wsqge(pr, taur, gamer, gdot,  &
         gamstar, pstar_lo, wrsq, clsqr, gmin, gmax)

    wl = ONE / sqrt(wlsq)
    wr = ONE / sqrt(wrsq)

    ustar_l = ul - (pstar_lo - pstar)*wl
    ustar_r = ur + (pstar_lo - pstar)*wr

    f_lo = ustar_l - ustar_r


    ! hi bounds
    call wsqge(pl, taul, gamel, gdot,  &
         gamstar, pstar_hi, wlsq, clsql, gmin, gmax)

    call wsqge(pr, taur, gamer, gdot,  &
         gamstar, pstar_hi, wrsq, clsqr, gmin, gmax)

    wl = ONE / sqrt(wlsq)
    wr = ONE / sqrt(wrsq)

    ustar_l = ul - (pstar_hi - pstar)*wl
    ustar_r = ur + (pstar_hi - pstar)*wr

    f_hi = ustar_l - ustar_r

    ! bisection
    iter = 1
    do while (iter <= cg_maxiter)

       pstar_c = HALF * (pstar_lo + pstar_hi)

       pstar_hist_extra(iter) = pstar_c

       call wsqge(pl, taul, gamel, gdot,  &
            gamstar, pstar_c, wlsq, clsql, gmin, gmax)

       call wsqge(pr, taur, gamer, gdot,  &
            gamstar, pstar_c, wrsq, clsqr, gmin, gmax)

       wl = ONE / sqrt(wlsq)
       wr = ONE / sqrt(wrsq)

       ustar_l = ul - (pstar_c - pl)*wl
       ustar_r = ur - (pstar_c - pr)*wr

       f_c = ustar_l - ustar_r

       if ( HALF * abs(pstar_lo - pstar_hi) < cg_tol * pstar_c ) then
          converged = .true.
          exit
       endif

       if (f_lo * f_c < ZERO) then
          ! root is in the left half
          pstar_hi = pstar_c
          f_hi = f_c
       else
          pstar_lo = pstar_c
          f_lo = f_c
       endif
    enddo

    pstar = pstar_c

  end subroutine pstar_bisection


  subroutine HLL(ql, qr, cl, cr, idir, ndim, f)

    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, QPRES, QREINT, &
         URHO, UMX, UMY, UMZ, UEDEN, UEINT, &
         npassive, upass_map, qpass_map
    use prob_params_module, only : coord_type

    implicit none
    double precision, intent(in) :: ql(QVAR), qr(QVAR), cl, cr
    double precision, intent(inout) :: f(NVAR)
    integer, intent(in) :: idir, ndim

    integer :: ivel, ivelt, iveltt, imom, imomt, imomtt
    double precision :: a1, a4, bd, bl, bm, bp, br
    double precision :: cavg, uavg
    double precision :: fl_tmp, fr_tmp
    double precision :: rhod, rhoEl, rhoEr, rhol_sqrt, rhor_sqrt
    integer :: n, nq

    integer :: ipassive

    double precision, parameter :: small = 1.d-10

    select case (idir)
    case (1)
       ivel = QU
       ivelt = QV
       iveltt = QW

       imom = UMX
       imomt = UMY
       imomtt = UMZ

    case (2)
       ivel = QV
       ivelt = QU
       iveltt = QW

       imom = UMY
       imomt = UMX
       imomtt = UMZ

    case (3)
       ivel = QW
       ivelt = QU
       iveltt = QV

       imom = UMZ
       imomt = UMX
       imomtt = UMY

    end select

    rhol_sqrt = sqrt(ql(QRHO))
    rhor_sqrt = sqrt(qr(QRHO))

    rhod = ONE/(rhol_sqrt + rhor_sqrt)



    ! compute the average sound speed. This uses an approximation from
    ! E88, eq. 5.6, 5.7 that assumes gamma falls between 1
    ! and 5/3
    cavg = sqrt( (rhol_sqrt*cl**2 + rhor_sqrt*cr**2)*rhod + &
         HALF*rhol_sqrt*rhor_sqrt*rhod**2*(qr(ivel) - ql(ivel))**2 )


    ! Roe eigenvalues (E91, eq. 5.3b)
    uavg = (rhol_sqrt*ql(ivel) + rhor_sqrt*qr(ivel))*rhod

    a1 = uavg - cavg
    a4 = uavg + cavg


    ! signal speeds (E91, eq. 4.5)
    bl = min(a1, ql(ivel) - cl)
    br = max(a4, qr(ivel) + cr)

    bm = min(ZERO, bl)
    bp = max(ZERO, br)

    bd = bp - bm

    if (abs(bd) < small*max(abs(bm),abs(bp))) return

    bd = ONE/bd


    ! compute the fluxes according to E91, eq. 4.4b -- note that the
    ! min/max above picks the correct flux if we are not in the star
    ! region

    ! density flux
    fl_tmp = ql(QRHO)*ql(ivel)
    fr_tmp = qr(QRHO)*qr(ivel)

    f(URHO) = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr(QRHO) - ql(QRHO))


    ! normal momentum flux.  Note for 1-d and 2-d non cartesian
    ! r-coordinate, we leave off the pressure term and handle that
    ! separately in the update, to accommodate different geometries
    if (ndim == 1 .or. (ndim == 2 .and. coord_type == 1 .and. idir == 1)) then
       fl_tmp = ql(QRHO)*ql(ivel)**2
       fr_tmp = qr(QRHO)*qr(ivel)**2
    else
       fl_tmp = ql(QRHO)*ql(ivel)**2 + ql(QPRES)
       fr_tmp = qr(QRHO)*qr(ivel)**2 + qr(QPRES)
    endif

    f(imom) = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr(QRHO)*qr(ivel) - ql(QRHO)*ql(ivel))


    ! transverse momentum flux
    fl_tmp = ql(QRHO)*ql(ivel)*ql(ivelt)
    fr_tmp = qr(QRHO)*qr(ivel)*qr(ivelt)

    f(imomt) = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr(QRHO)*qr(ivelt) - ql(QRHO)*ql(ivelt))


    fl_tmp = ql(QRHO)*ql(ivel)*ql(iveltt)
    fr_tmp = qr(QRHO)*qr(ivel)*qr(iveltt)

    f(imomtt) = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr(QRHO)*qr(iveltt) - ql(QRHO)*ql(iveltt))


    ! total energy flux
    rhoEl = ql(QREINT) + HALF*ql(QRHO)*(ql(ivel)**2 + ql(ivelt)**2 + ql(iveltt)**2)
    fl_tmp = ql(ivel)*(rhoEl + ql(QPRES))

    rhoEr = qr(QREINT) + HALF*qr(QRHO)*(qr(ivel)**2 + qr(ivelt)**2 + qr(iveltt)**2)
    fr_tmp = qr(ivel)*(rhoEr + qr(QPRES))

    f(UEDEN) = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(rhoEr - rhoEl)


    ! eint flux
    fl_tmp = ql(QREINT)*ql(ivel)
    fr_tmp = qr(QREINT)*qr(ivel)

    f(UEINT) = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr(QREINT) - ql(QREINT))


    ! passively-advected scalar fluxes
    do ipassive = 1, npassive
       n  = upass_map(ipassive)
       nq = qpass_map(ipassive)

       fl_tmp = ql(QRHO)*ql(nq)*ql(ivel)
       fr_tmp = qr(QRHO)*qr(nq)*qr(ivel)

       f(n) = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr(QRHO)*qr(nq) - ql(QRHO)*ql(nq))
    enddo

  end subroutine HLL


  pure subroutine cons_state(q, U)

    use meth_params_module, only: QVAR, QRHO, QU, QV, QW, QREINT, &
         NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, &
         npassive, upass_map, qpass_map

    implicit none
    real (amrex_real), intent(in)  :: q(QVAR)
    real (amrex_real), intent(out) :: U(NVAR)

    integer :: ipassive, n, nq

    U(URHO) = q(QRHO)

    ! since we advect all 3 velocity components regardless of dimension, this
    ! will be general
    U(UMX)  = q(QRHO)*q(QU)
    U(UMY)  = q(QRHO)*q(QV)
    U(UMZ)  = q(QRHO)*q(QW)

    U(UEDEN) = q(QREINT) + HALF*q(QRHO)*(q(QU)**2 + q(QV)**2 + q(QW)**2)
    U(UEINT) = q(QREINT)

    ! we don't care about T here, but initialize it to make NaN
    ! checking happy
    U(UTEMP) = ZERO

    do ipassive = 1, npassive
       n  = upass_map(ipassive)
       nq = qpass_map(ipassive)
       U(n) = q(QRHO)*q(nq)
    enddo

  end subroutine cons_state


  pure subroutine HLLC_state(idir, S_k, S_c, q, U)

    use meth_params_module, only: QVAR, QRHO, QU, QV, QW, QREINT, QPRES, &
         NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, &
         npassive, upass_map, qpass_map

    implicit none
    integer, intent(in) :: idir
    real (amrex_real), intent(in)  :: S_k, S_c
    real (amrex_real), intent(in)  :: q(QVAR)
    real (amrex_real), intent(out) :: U(NVAR)

    real (amrex_real) :: hllc_factor, u_k
    integer :: ipassive, n, nq

    if (idir == 1) then
       u_k = q(QU)
    elseif (idir == 2) then
       u_k = q(QV)
    elseif (idir == 3) then
       u_k = q(QW)
    endif

    hllc_factor = q(QRHO)*(S_k - u_k)/(S_k - S_c)
    U(URHO) = hllc_factor
    if (idir == 1) then
       U(UMX)  = hllc_factor*S_c
       U(UMY)  = hllc_factor*q(QV)
       U(UMZ)  = hllc_factor*q(QW)
    elseif (idir == 2) then
       U(UMX)  = hllc_factor*q(QU)
       U(UMY)  = hllc_factor*S_c
       U(UMZ)  = hllc_factor*q(QW)
    elseif (idir == 3) then
       U(UMX)  = hllc_factor*q(QU)
       U(UMY)  = hllc_factor*q(QV)
       U(UMZ)  = hllc_factor*S_c
    endif

    U(UEDEN) = hllc_factor*(q(QREINT)/q(QRHO) + &
         HALF*(q(QU)**2 + q(QV)**2 + q(QW)**2) + &
         (S_c - u_k)*(S_c + q(QPRES)/(q(QRHO)*(S_k - u_k))))
    U(UEINT) = hllc_factor*q(QREINT)/q(QRHO)

    U(UTEMP) = ZERO  ! we don't evolve T

    do ipassive = 1, npassive
       n  = upass_map(ipassive)
       nq = qpass_map(ipassive)
       U(n) = hllc_factor*q(nq)
    enddo

  end subroutine HLLC_state


  pure subroutine compute_flux(idir, ndim, bnd_fac, U, p, F)

    use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, &
         npassive, upass_map
    use prob_params_module, only : coord_type

    implicit none
    integer, intent(in) :: idir, ndim, bnd_fac
    real (amrex_real), intent(in) :: U(NVAR)
    real (amrex_real), intent(in) :: p
    real (amrex_real), intent(out) :: F(NVAR)

    integer :: ipassive, n
    real (amrex_real) :: u_flx

    if (idir == 1) then
       u_flx = U(UMX)/U(URHO)
    elseif (idir == 2) then
       u_flx = U(UMY)/U(URHO)
    elseif (idir == 3) then
       u_flx = U(UMZ)/U(URHO)
    endif

    if (bnd_fac == 0) then
       u_flx = ZERO
    endif

    F(URHO) = U(URHO)*u_flx

    F(UMX) = U(UMX)*u_flx
    F(UMY) = U(UMY)*u_flx
    F(UMZ) = U(UMZ)*u_flx

    if (ndim == 3 .or. (ndim == 2 .and. coord_type == 0) .or. &
         (ndim == 2 .and. coord_type == 1 .and. idir == 2)) then
       ! we do not include the pressure term in any non-Cartesian
       ! coordinate directions
       F(UMX-1+idir) = F(UMX-1+idir) + p
    endif

    F(UEINT) = U(UEINT)*u_flx
    F(UEDEN) = (U(UEDEN) + p)*u_flx

    F(UTEMP) = ZERO

    do ipassive = 1, npassive
       n = upass_map(ipassive)
       F(n) = U(n)*u_flx
    enddo

  end subroutine compute_flux

  subroutine riemann_md_singlepoint( rl, ul, vl, v2l, pl, rel, spl, gamcl, &
       rr, ur, vr, v2r, pr, rer, spr, gamcr,&
       qint_iu, qint_iv1, qint_iv2, qint_gdpres, qint_gdgame, &
       regd, rgd, ustar, eos_state, gdnv_state, nsp, &
       uflx_rho, uflx_u, uflx_v, uflx_w, uflx_eden, uflx_eint, &
       idir, coord_type, bc_test_val, csmall, cav)

    use eos_module
    use meth_params_module, only : small_dens, small_pres

    implicit none

    ! Inputs
    integer, intent(in) :: nsp
    double precision, intent(in) :: rl, ul, vl, v2l, pl, rel, spl(nsp), gamcl ! Left state
    double precision, intent(in) :: rr, ur, vr, v2r, pr, rer, spr(nsp), gamcr ! Right state
    integer, intent(in) :: idir, coord_type
    integer, intent(in) :: bc_test_val
    double precision, intent(in) :: csmall, cav

    ! Work values sent back to compute passive scalar flux
    double precision, intent(out) :: rgd, regd, ustar 
    double precision, intent(out) :: qint_iu, qint_iv1, qint_iv2, qint_gdpres, qint_gdgame

    ! Work arrays passed in becase allocation/construction is expensive, do it elsewhere
    type(eos_t), intent(inout) :: eos_state, gdnv_state


    ! Outputs

    double precision, intent(inout):: uflx_rho, uflx_u, uflx_v, uflx_w, uflx_eden, uflx_eint

    ! Local work variables
    double precision :: vgd, wgd, csr, csl, wl, wr, rhoetot, scr 
    double precision :: rstar, cstar, estar, pstar
    double precision :: ro, uo, po, reo, co, drho
    double precision :: sgnm, spin, spout, ushock, frac
    double precision :: wsmall

    ! Yuck.
    double precision, parameter:: small = 1.d-8
    real (amrex_real), parameter :: smallu = 1.e-12_amrex_real

    wsmall = small_dens*csmall

    gdnv_state%rho = rl
    gdnv_state%p = pl
    !unncessary for gamma law
    gdnv_state%massfrac = spl
    !dir$ inline recursive
    call eos_rp(gdnv_state)
    csl = gdnv_state%cs

    gdnv_state%rho = rr
    gdnv_state%p = pr
    !unncessary for gamma law
    gdnv_state%massfrac = spr
    !dir$ inline recursive
    call eos_rp(gdnv_state)
    csr = gdnv_state%cs

    wl = max(wsmall,sqrt(abs(gamcl*pl*rl)))
    wr = max(wsmall,sqrt(abs(gamcr*pr*rr)))

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
       !unncessary for gamma law
       gdnv_state%massfrac = spl
    else if (ustar .lt. ZERO) then
       ro = rr
       uo = ur
       po = pr
       !unncessary for gamma law
       gdnv_state%massfrac = spr
    else
       ro = HALF*(rl+rr)
       uo = HALF*(ul+ur)
       po = HALF*(pl+pr)
       !unncessary for gamma law
       gdnv_state%massfrac = HALF*(spl+spr)
    endif

    gdnv_state%rho = ro
    gdnv_state%p = po
    !dir$ inline recursive
    call eos_rp(gdnv_state)
    reo = gdnv_state%rho * gdnv_state%e
    co = gdnv_state%cs

    drho = (pstar - po)/co**2
    rstar = ro + drho
    rstar = max(small_dens,rstar)

    ! At star state, mass fractions are upwinded, have rho, p.  Calc c and rhoe
    gdnv_state%rho = rstar
    gdnv_state%p = pstar
    !dir$ inline recursive
    call eos_rp(gdnv_state)
    cstar = gdnv_state%cs
    estar = gdnv_state%rho * gdnv_state%e

    sgnm = sign(ONE,ustar)
    spout = co - sgnm*uo
    spin = cstar - sgnm*ustar

    ushock = HALF*(spin + spout)

    if (pstar-po >= ZERO) then
       spin = ushock
       spout = ushock
    endif

    if (spout-spin == ZERO) then
       scr = small*cav
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

    qint_iu = frac*ustar + (ONE - frac)*uo
    qint_iv1 = vgd
    qint_iv2 = wgd

    qint_gdpres = frac*pstar + (ONE - frac)*po


    gdnv_state%rho = rgd
    gdnv_state%p = qint_gdpres
    !dir$ inline recursive
    call eos_rp(gdnv_state)
    regd = gdnv_state%rho * gdnv_state%e

    if (spout < ZERO) then
       rgd = ro
       qint_iu = uo
       qint_gdpres = po
       regd = reo
    endif

    if (spin >= ZERO) then
       rgd = rstar
       qint_iu = ustar
       qint_gdpres = pstar
       regd = estar
    endif

    gdnv_state%rho = rgd
    gdnv_state%p = qint_gdpres
    !dir$ inline recursive
    call eos_rp(gdnv_state)
    regd = gdnv_state%rho * gdnv_state%e


    ! not sure what this should be for radiation?
    qint_gdgame = qint_gdpres/regd + ONE

    ! enforce that the fluxes through a symmetry plane or wall are zero
    qint_iu = bc_test_val * qint_iu

    ! Compute fluxes, order as conserved state (not q)
    uflx_rho = rgd*qint_iu

    ! note: for axisymmetric geometries, we do not include the
    ! pressure in the r-direction, since div{F} + grad{p} cannot
    ! be written in a flux difference form
    if (idir == 1) then
       uflx_u = uflx_rho*qint_iu
       uflx_v = uflx_rho*qint_iv1
       uflx_w = uflx_rho*qint_iv2
       if (coord_type == 0) then
          uflx_u = uflx_u + qint_gdpres
       endif
    else if (idir == 2) then
       uflx_u = uflx_rho*qint_iv1
       uflx_v = uflx_rho*qint_iu + qint_gdpres
       uflx_w = uflx_rho*qint_iv2
    else
       uflx_u = uflx_rho*qint_iv1
       uflx_v = uflx_rho*qint_iv2
       uflx_w = uflx_rho*qint_iu + qint_gdpres
    endif

    rhoetot = regd + HALF*rgd*(qint_iu**2 + qint_iv1**2 + qint_iv2**2)

    uflx_eden = qint_iu*(rhoetot + qint_gdpres)
    uflx_eint = qint_iu*regd


  end subroutine riemann_md_singlepoint

  !
  !> Single point wrapper for the core of the Riemann solver
  !  this is called from the ?d versions of the routines;
  !  where performance is important one needs to check that it 
  !  is (1) inlined and (2) that vector instructions are generated

  ! Single points of work array (qint_iu, qint_iv1, etc) are passed in direclty;
  ! partially so that this has no dependence on the component index variables (e.g. QPRES),
  ! and partially to avoid creation of a temporary (passing qint(i,j,:) is bad...)
  subroutine riemann_md_vec( rl, ul, vl, v2l, pl, rel, spl, gamcl, &
       rr, ur, vr, v2r, pr, rer, spr, gamcr,&
       qint_iu, vgd, wgd, qint_gdpres, qint_gdgame, &
       regd, rgd, ustar, gdnv_state, nsp, &
       uflx_rho, uflx_u, uflx_v, uflx_w, uflx_eden, uflx_eint, &
       bc_test_val, csmall, cav, VECLEN)

    use eos_module
    use meth_params_module, only : small_dens, small_pres

    implicit none

    ! Inputs
    integer, intent(in) :: nsp, VECLEN
    double precision, intent(in), dimension(VECLEN) :: rl, ul, vl, v2l, pl, rel, gamcl 
    double precision, intent(in), dimension(VECLEN) :: rr, ur, vr, v2r, pr, rer,  gamcr ! Right state
    double precision, intent(in), dimension(VECLEN,nsp) :: spl, spr
    ! always idir=1          integer, intent(in) :: idir, coord_type
    integer, intent(in) :: bc_test_val
    double precision, intent(in), dimension(VECLEN) :: csmall, cav

    ! Work values sent back to compute passive scalar flux
    double precision, intent(out), dimension(VECLEN) :: rgd, regd, ustar 
    double precision, intent(out), dimension(VECLEN) :: qint_iu, vgd, wgd, qint_gdpres, qint_gdgame

    ! Work arrays passed in becase allocation/construction is expensive, do it elsewhere
    type(eos_t), intent(inout) :: gdnv_state


    ! Outputs

    double precision, intent(inout), dimension(VECLEN) :: uflx_rho, uflx_u, uflx_v, uflx_w, uflx_eden, uflx_eint

    ! Local work variables
    double precision, dimension(VECLEN) :: csr, csl, wl, wr, rhoetot, scr 
    double precision, dimension(VECLEN) :: rstar, cstar, estar, pstar
    double precision, dimension(VECLEN) :: ro, uo, po, reo, co, drho
    double precision, dimension(VECLEN) :: sgnm, spin, spout, ushock, frac
    double precision, dimension(VECLEN) :: wsmall, psmall, rsmall
    double precision, dimension(VECLEN,nsp) :: sp

    integer :: vii, n

    ! Yuck.
    double precision, parameter:: small = 1.d-8
    real (amrex_real), parameter :: smallu = 1.e-12_amrex_real
    real (amrex_real), parameter :: Hsmallu = 0.5e-12_amrex_real

    rsmall = small_dens
    wsmall = small_dens*csmall
    psmall = small_pres

    do vii = 1, VECLEN
       gdnv_state%rho = rl(vii)
       gdnv_state%p = pl(vii)
       !unncessary for gamma law
       gdnv_state%massfrac = spl(vii,:)
       !dir$ inline recursive
       call eos_rp(gdnv_state)
       csl(vii) = gdnv_state%cs
    enddo

    do vii = 1, VECLEN
       gdnv_state%rho = rr(vii)
       gdnv_state%p = pr(vii)
       !unncessary for gamma law
       gdnv_state%massfrac = spr(vii,:)
       !dir$ inline recursive
       call eos_rp(gdnv_state)
       csr(vii) = gdnv_state%cs
    enddo

    wl = sqrt(abs(gamcl*pl*rl))
    wl = merge(wsmall, wl, wl<wsmall)

    wr = sqrt(abs(gamcr*pr*rr))
    wr = merge(wsmall, wr, wr<wsmall)

    pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))/(wl + wr)
    ustar = ((wl*ul + wr*ur) + (pl - pr))/(wl + wr)

    pstar = merge(psmall,pstar,pstar<psmall)

    ro = merge(rl, rr, ustar > ZERO)
    uo = merge(ul, ur, ustar > ZERO)
    po = merge(pl, pr, ustar > ZERO)
    do n=1,nsp
       sp(:,n) = merge(spl(:,n), spr(:,n), ustar > ZERO)
    enddo
    ! for symmetry preservation, if ustar is really small, then we
    ! set it to zero
    where(abs(ustar) < Hsmallu*(abs(ul) + abs(ur)) .or. ustar.eq.0)
       ustar = ZERO
       ro = HALF*(rl+rr)
       uo = HALF*(ul+ur)
       po = HALF*(pl+pr)
    end where
    do n=1,nsp
       where(abs(ustar) < Hsmallu*(abs(ul) + abs(ur)) .or. ustar.eq.0)
          sp(:,n) = HALF*(spl(:,n)+spr(:,n))
       end where
    enddo

    do vii = 1, VECLEN
       gdnv_state%rho = ro(vii)
       gdnv_state%p = po(vii)
       !unncessary for gamma law
       gdnv_state%massfrac = sp(vii,:)
       !dir$ inline recursive
       call eos_rp(gdnv_state)
       reo(vii) = gdnv_state%rho * gdnv_state%e
       co(vii) = gdnv_state%cs
    enddo

    drho = (pstar - po)/co**2
    rstar = ro + drho
    rstar = merge(rsmall, rstar, rstar<rsmall)

    ! At star state, mass fractions are upwinded, have rho, p.  Calc c and rhoe
    do vii = 1, VECLEN
       gdnv_state%rho = rstar(vii)
       gdnv_state%p = pstar(vii)
       gdnv_state%massfrac = sp(vii,:)
       !dir$ inline recursive
       call eos_rp(gdnv_state)
       cstar(vii) = gdnv_state%cs
       estar(vii) = gdnv_state%rho * gdnv_state%e
    enddo

    ! No is there a vector version of sign? dunno.
    do vii = 1, VECLEN
       sgnm(vii) = sign(ONE,ustar(vii))
    enddo

    spout = co - sgnm*uo
    spin = cstar - sgnm*ustar

    ushock = HALF*(spin + spout)

    spout = merge(spout, ushock, pstar<po)
    spin = merge(spin, ushock, pstar < po)

    scr = merge(small*cav, spout-spin, spout .eq. spin)

    frac = (ONE + (spout + spin)/scr)*HALF
    frac = max(ZERO,min(ONE,frac))

    vgd = merge(vl, vr, ustar > ZERO)
    wgd = merge(v2l, v2r, ustar > ZERO)

    where (ustar .eq. ZERO)
       vgd = HALF*(vl+vr)
       wgd = HALF*(v2l+v2r)
    end where

    rgd = frac*rstar + (ONE - frac)*ro

    qint_iu = frac*ustar + (ONE - frac)*uo
    qint_gdpres = frac*pstar + (ONE - frac)*po

    do vii = 1, VECLEN
       gdnv_state%rho = rgd(vii)
       gdnv_state%p = qint_gdpres(vii)
       gdnv_state%massfrac = sp(vii,:)
       !dir$ inline recursive
       call eos_rp(gdnv_state)
       regd(vii) = gdnv_state%rho * gdnv_state%e
    enddo

    where (spout < ZERO)
       rgd = ro
       qint_iu = uo
       qint_gdpres = po
       regd = reo
    end where

    where (spin >= ZERO) 
       rgd = rstar
       qint_iu = ustar
       qint_gdpres = pstar
       regd = estar
    end where

    do vii = 1, VECLEN
       gdnv_state%rho = rgd(vii)
       gdnv_state%p = qint_gdpres(vii)
       gdnv_state%massfrac = sp(vii,:)
       !dir$ inline recursive
       call eos_rp(gdnv_state)
       regd(vii) = gdnv_state%rho * gdnv_state%e
    enddo

    ! not sure what this should be for radiation?
    qint_gdgame = qint_gdpres/regd + ONE

    ! enforce that the fluxes through a symmetry plane or wall are zero
    qint_iu = bc_test_val * qint_iu

    ! Compute fluxes, order as conserved state (not q)
    uflx_rho = rgd*qint_iu

    ! note: for axisymmetric geometries, we should not include the
    ! pressure in the r-direction, since div{F} + grad{p} cannot
    ! be written in a flux difference form
    uflx_u = uflx_rho*qint_iu + qint_gdpres
    uflx_v = uflx_rho*vgd
    uflx_w = uflx_rho*wgd

    rhoetot = regd + HALF*rgd*(qint_iu**2 + vgd**2 + wgd**2)

    uflx_eden = qint_iu*(rhoetot + qint_gdpres)
    uflx_eint = qint_iu*regd


  end subroutine riemann_md_vec

  subroutine riemann_cg_singlepoint( rl, ul, vl, v2l, pl, rel, spl, gamcl, &
       rr, ur, vr, v2r, pr, rer, spr, gamcr,&
       ugd, v1gd, v2gd, pgd, gamegd, regd, rgd, ustar, nsp, &
       uflx_rho, uflx_u, uflx_v, uflx_w, uflx_eden, uflx_eint, &
       idir, coord_type, bc_test_val, csmall, cav)

    use meth_params_module, only : small_dens, small_pres, cg_maxiter, cg_tol, cg_blend

    implicit none

    ! Inputs
    integer, intent(in) :: nsp
    double precision, intent(in) :: rl, ul, vl, v2l, pl, rel, spl(nsp), gamcl ! Left state
    double precision, intent(in) :: rr, ur, vr, v2r, pr, rer, spr(nsp), gamcr ! Right state
    integer, intent(in) :: idir, coord_type
    integer, intent(in) :: bc_test_val
    double precision, intent(in) :: csmall, cav

    ! Work values sent back to compute passive scalar flux
    double precision, intent(in) :: regd !currently not used, and only set to intent(in) to avoid warning
    double precision, intent(out) :: rgd, ustar
    double precision, intent(out) :: ugd, v1gd, v2gd, pgd, gamegd

    ! Outputs
    double precision, intent(inout):: uflx_rho, uflx_u, uflx_v, uflx_w, uflx_eden, uflx_eint

    double precision, parameter:: small = 1.d-8
    double precision, parameter :: small_u = 1.d-10

    double precision :: wl, wr, rstar, cstar, pstar
    double precision :: ro, uo, po, co, gamco
    double precision :: sgnm, spin, spout, ushock, frac

    double precision :: clsq, clsql, clsqr, wlsq, wosq, wrsq, wo
    double precision :: zl, zr
    double precision :: denom, dpditer, dpjmp
    double precision :: gamc_bar, game_bar
    double precision :: gamel, gamer, gameo, gamstar, gmin, gmax, gdot

    integer :: iter
    double precision :: err

    logical :: converged

    double precision :: pstar_old,wsmall
    double precision :: taul, taur, tauo
    double precision :: ustar_r, ustar_l, ustar_r_old, ustar_l_old
    double precision :: pstar_lo, pstar_hi

    double precision, parameter :: weakwv = 1.d-3

    double precision :: pstar_hist(cg_maxiter), pstar_hist_extra(cg_maxiter)

    ! common quantities
    taul = ONE/rl
    taur = ONE/rr

    ! lagrangian sound speeds
    clsql = gamcl*pl*rl
    clsqr = gamcr*pr*rr


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
    gamc_bar = HALF*(gamcl + gamcr)

    gdot = TWO*(ONE - game_bar/gamc_bar)*(game_bar - ONE)

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
    do while ((iter <= cg_maxiter .and. .not. converged) .or. iter <= 2)

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
       if (err < cg_tol*pstar) converged = .true.

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
          do iter = 1, cg_maxiter
             print *, iter, pstar_hist(iter)
          enddo

          print *, ' '
          print *, 'left state  (r,u,p,re,gc): ', rl, ul, pl, rel, gamcl
          print *, 'right state (r,u,p,re,gc): ', rr, ur, pr, rer, gamcr
          print *, 'cav, smallc:',  cav, csmall
          call bl_error("ERROR: non-convergence in the Riemann solver")

       else if (cg_blend .eq. 1) then

          pstar = pl + ( (pr - pl) - wr*(ur - ul) )*wl/(wl+wr)

       else if (cg_blend .eq. 2) then

          ! first try to find a reasonable bounds 
          pstar_lo = minval(pstar_hist(cg_maxiter-5:cg_maxiter))
          pstar_hi = maxval(pstar_hist(cg_maxiter-5:cg_maxiter))

          call pstar_bisection(pstar_lo, pstar_hi, &
               ul, pl, taul, gamel, clsql, &
               ur, pr, taur, gamer, clsqr, &
               gdot, gmin, gmax, &
               pstar, gamstar, converged, pstar_hist_extra)

          if (.not. converged) then
             ! abort -- doesn't seem solvable
             print *, 'pstar history: '
             do iter = 1, cg_maxiter
                print *, iter, pstar_hist(iter)
             enddo

             do iter = 1, cg_maxiter
                print *, iter+cg_maxiter, pstar_hist_extra(iter)
             enddo

             print *, ' '
             print *, 'left state  (r,u,p,re,gc): ', rl, ul, pl, rel, gamcl
             print *, 'right state (r,u,p,re,gc): ', rr, ur, pr, rer, gamcr
             print *, 'cav, smallc:',  cav, csmall
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
    if (abs(ustar) < csmall*HALF*(abs(ul) + abs(ur))) then
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
       gamco = gamcl
       gameo = gamel
    else if (ustar .lt. ZERO) then
       ro = rr
       uo = ur
       po = pr
       tauo = taur
       !reo = rer
       gamco = gamcr
       gameo = gamer
    else
       ro = HALF*(rl+rr)
       uo = HALF*(ul+ur)
       po = HALF*(pl+pr)
       tauo = HALF*(taul+taur)
       !reo = HALF*(rel+rer)
       gamco = HALF*(gamcl + gamcr)
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

    frac = HALF*(ONE + (spin + spout)/max(spout-spin,spin+spout, small*cav))

    ! the transverse velocity states only depend on the
    ! direction that the contact moves
    if (ustar .gt. ZERO) then
       v1gd = vl
       v2gd = v2l
    else if (ustar .lt. ZERO) then
       v1gd = vr
       v2gd = v2r
    else
       v1gd = HALF*(vl+vr)
       v2gd = HALF*(v2l+v2r)
    endif

    ! linearly interpolate between the star and normal state -- this covers the
    ! case where we are inside the rarefaction fan.
    rgd = frac*rstar + (ONE - frac)*ro
    ugd = frac*ustar + (ONE - frac)*uo
    pgd = frac*pstar + (ONE - frac)*po
    gamegd =  frac*gamstar + (ONE-frac)*gameo

    ! now handle the cases where instead we are fully in the
    ! star or fully in the original (l/r) state
    if (spout .lt. ZERO) then
       rgd = ro
       ugd = uo
       pgd = po
       gamegd = gameo
    endif
    if (spin .ge. ZERO) then
       rgd = rstar
       ugd = ustar
       pgd = pstar
       gamegd = gamstar
    endif

    pgd = max(pgd,small_pres)

    ! Enforce that fluxes through a symmetry plane or wall are hard zero.
    ! Here the NSCBC info about if we have a wall or not is contained in the ghost-cell
    ugd = bc_test_val * ugd

  end subroutine riemann_cg_singlepoint

end module riemann_util_module
