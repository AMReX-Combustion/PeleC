module hyp_advection_module 

  use amrex_ebcellflag_module, only : get_neighbor_cells
  use pelec_eb_stencil_types_module, only : eb_bndry_geom
  use riemann_util_module, only : riemann_md_singlepoint, riemann_cg_singlepoint
  use prob_params_module, only: dim

  implicit none 
  private 
  public pc_hyp_mol_flux
  contains 

  !> Computes fluxes for hyperbolic conservative update.
  !> @brief 
  !> Uses MOL formulation
  !! @param[inout] flux1  flux in X direction on X edges
  !> @param[in] q        (const)  input state, primitives
  !> @param[in] flatn    (const)  flattening parameter
  !> @param[in] src      (const)  source
  !> @param[in] nx       (const)  number of cells in X direction
  !> @param[in] ny       (const)  number of cells in Y direction
  !> @param[in] nz       (const)  number of cells in Z direction
  !> @param[in] dx       (const)  grid spacing in X direction
  !> @param[in] dy       (const)  grid spacing in Y direction
  !> @param[in] dz       (const)  grid spacing in Z direction
  !> @param[in] dt       (const)  time stepsize
  !> @param[inout] flux1    (modify) flux in X direction on X edges
  !> @param[inout] flux2    (modify) flux in Y direction on Y edges
  !> @param[inout] flux3    (modify) flux in Z direction on Z edges
  subroutine pc_hyp_mol_flux(lo, hi, &
                     domlo, domhi, &
                     q, qd_lo, qd_hi, &
                     qaux, qa_lo, qa_hi, &
                     Ax,  Axlo,  Axhi,&
                     flux1, fd1_lo, fd1_hi, &
                     Ay,  Aylo,  Ayhi,&
                     flux2, fd2_lo, fd2_hi, &
                     flatn, fltd_lo, fltd_hi, &
                     V, Vlo, Vhi, &
                     D, Dlo, Dhi,&
#ifdef PELEC_USE_EB
                     vfrac, vflo, vfhi, &
                     flag, fglo, fghi, &
                     ebg, Nebg, ebflux, nebflux, &
#endif
                     h) bind(C,name="pc_hyp_mol_flux")

    use amrex_mempool_module, only : bl_allocate, bl_deallocate
    use meth_params_module, only : QVAR, NVAR, QPRES, QRHO, QU, QV, QFS, QC, QCSML, NQAUX, &
                                   URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, UTEMP, UFA, UFX, nadv, &
                                   riemann_solver

    use slope_module, only : slopex, slopey
    use network, only : nspecies, naux
    use eos_type_module
    use eos_module, only : eos_t, eos_rp
    use riemann_module, only: cmpflx, shock
    use amrex_constants_module
    use amrex_fort_module, only : amrex_real

    implicit none

    integer, parameter :: VECLEN = 16
    integer :: is, ie, ic, vi, vii

    integer, intent(in) ::      qd_lo(2),   qd_hi(2)
    integer, intent(in) ::      qa_lo(2),   qa_hi(2)
    integer, intent(in) ::         lo(2),      hi(2)
    integer, intent(in) ::      domlo(2),   domhi(2)
    integer, intent(in) ::       Axlo(2),    Axhi(2)
    integer, intent(in) ::     fd1_lo(2),  fd1_hi(2)
    integer, intent(in) ::       Aylo(2),    Ayhi(2)
    integer, intent(in) ::     fd2_lo(2),  fd2_hi(2)
    integer, intent(in) ::    fltd_lo(2), fltd_hi(2)
    integer, intent(in) ::        Vlo(2),     Vhi(2)
    integer, intent(in) ::        Dlo(2),     Dhi(2)
    double precision, intent(in) :: h(2)

#ifdef PELEC_USE_EB
    integer, intent(in) ::  vflo(2),    vfhi(2)
    integer, intent(in) ::  fglo(2),    fghi(2)
    integer, intent(in) :: flag(fglo(1):fghi(1),fglo(2):fghi(2))
    real(amrex_real), intent(in) :: vfrac(vflo(1):vfhi(1),vflo(2):vfhi(2))

    integer, intent(in) :: nebflux
    real(amrex_real), intent(inout) ::   ebflux(0:nebflux-1,1:NVAR)
    integer,            intent(in   ) :: Nebg
    type(eb_bndry_geom),intent(in   ) :: ebg(0:Nebg-1)    
    real(amrex_real) :: eb_norm(2), tflux(3), full_area
#endif
    double precision, intent(in) ::     q(  qd_lo(1):  qd_hi(1),  qd_lo(2):  qd_hi(2),QVAR)  !> State
    double precision, intent(in) ::  qaux(  qa_lo(1):  qa_hi(1),  qa_lo(2):  qa_hi(2),NQAUX) !> Auxiliary state
    double precision, intent(in) :: flatn(fltd_lo(1):fltd_hi(1),fltd_lo(2):fltd_hi(2))

    double precision, intent(in   ) ::    Ax(  Axlo(1):  Axhi(1),  Axlo(2):  Axhi(2))
    double precision, intent(inout) :: flux1(fd1_lo(1):fd1_hi(1),fd1_lo(2):fd1_hi(2),NVAR)
    double precision, intent(in   ) ::    Ay(  Aylo(1):  Ayhi(1),  Aylo(2):  Ayhi(2))
    double precision, intent(inout) :: flux2(fd2_lo(1):fd2_hi(1),fd2_lo(2):fd2_hi(2),NVAR)
    double precision, intent(inout) ::     V(   Vlo(1):   Vhi(1),   Vlo(2):   Vhi(2))
    double precision, intent(inout) ::     D(   Dlo(1):   Dhi(1),   Dlo(2):   Dhi(2),NVAR)

    integer :: i, j, nsp, L, ivar
    integer :: qt_lo(2), qt_hi(2)

    ! Left and right state arrays (edge centered, cell centered)
    double precision, pointer :: dqx(:,:,:), dqy(:,:,:)

    ! Other left and right state arrays
    double precision :: qtempl(VECLEN,1:5+nspecies)
    double precision :: qtempr(VECLEN,1:5+nspecies)
    double precision :: rhoe_l(VECLEN)
    double precision :: rhoe_r(VECLEN)
    double precision :: cspeed(VECLEN)
    double precision :: gamc_l(VECLEN)
    double precision :: gamc_r(VECLEN)
    double precision :: cav(VECLEN)
    double precision :: csmall(VECLEN)

    ! Riemann solve work arrays
    double precision, dimension(VECLEN) :: u_gd, v_gd, w_gd, &
                                           p_gd, game_gd, re_gd, &
                                           r_gd, ustar
    double precision :: flux_tmp(VECLEN, NVAR)
    integer, parameter :: idir = 1
    integer :: nextra
    integer, parameter :: coord_type = 0
    integer, parameter :: bc_test_val = 1

    type (eos_t) :: eos_state, gdnv_state

    integer, parameter :: R_RHO = 1
    integer, parameter :: R_UN  = 2
    integer, parameter :: R_UT1 = 3
    integer, parameter :: R_UT2 = 4
    integer, parameter :: R_P   = 5
    integer, parameter :: R_Y   = 6

!   concept is to advance cells lo to hi
!   need fluxes on the boundary
!   if tile is eb need to expand by 2 cells in each directions
!   would like to do this tile by tile
#ifdef PELEC_USE_EB
    nextra = 3
#else
    nextra = 0
#endif
    do L=1,dim
       qt_lo(L) = lo(L) - nextra
       qt_hi(L) = hi(L) + nextra
       if (qt_lo(L)-1 .lt. qd_lo(L) .or. qt_hi(L)+1 .gt. qd_hi(L)) then
          call bl_pd_abort()
       endif
    enddo

    ! allocate space for eos calls
    call build(eos_state)
    call build(gdnv_state)

    !  allocate space for slopes
    call bl_allocate ( dqx, qt_lo(1), qt_hi(1), qt_lo(2), qt_hi(2), 1, QVAR)
    call bl_allocate ( dqy, qt_lo(1), qt_hi(1), qt_lo(2), qt_hi(2), 1, QVAR)

#ifdef PELEC_USE_EB
    call slopex(q,flatn,qd_lo,qd_hi, &
                   dqx,qt_lo,qt_hi, &
                   lo(1)-nextra,lo(2)-nextra,  &
                   hi(1)+nextra,hi(2)+nextra,QVAR,NQAUX, &
                   domlo,domhi, &
                   qaux, qa_lo, qa_hi, &
                   flag, fglo, fghi)

#else
    call slopex(q,flatn,qd_lo,qd_hi, &
                   dqx,qt_lo,qt_hi, &
                   lo(1),lo(2),hi(1),hi(2),QVAR &
                   qaux, qa_lo, qa_hi)
#endif

    do j = lo(2)-nextra, hi(2)+nextra
       do i = lo(1)-nextra+1, hi(1)+nextra, VECLEN
          is = i
          ie = min(is+VECLEN-1, hi(1)+nextra)
          ic = ie - is + 1

          cspeed(1:ic) = qaux(is-1:ie-1,j,QC)
          qtempl(1:ic,R_UN ) = q(is-1:ie-1,j,QU   ) + 0.5d0*((dqx(is-1:ie-1,j,2)-dqx(is-1:ie-1,j,1))/q(is-1:ie-1,j,QRHO))
          qtempl(1:ic,R_P  ) = q(is-1:ie-1,j,QPRES) + 0.5d0*( dqx(is-1:ie-1,j,1)+dqx(is-1:ie-1,j,2) )*cspeed(1:ic)
          qtempl(1:ic,R_UT1) = q(is-1:ie-1,j,QV   ) + 0.5d0 * dqx(is-1:ie-1,j,3)
          qtempl(1:ic,R_UT2) = 0.d0

          qtempl(1:ic,R_RHO) = 0.d0
          do nsp = 1,nspecies
            qtempl(1:ic,R_Y-1+nsp) = q(is-1:ie-1,j,QFS-1+nsp)*q(is-1:ie-1,j,QRHO) + 0.5d0*(dqx(is-1:ie-1,j,4+nsp) &
                 + q(is-1:ie-1,j,QFS-1+nsp)*(dqx(is-1:ie-1,j,1)+dqx(is-1:ie-1,j,2))/cspeed(1:ic) )
            qtempl(1:ic,R_RHO) = qtempl(1:ic,R_RHO) + qtempl(1:ic,R_Y-1+nsp)
          enddo

          do nsp = 1,nspecies
            qtempl(1:ic,R_Y-1+nsp) = qtempl(1:ic,R_Y-1+nsp)/qtempl(1:ic,R_RHO)
          enddo

          cspeed(1:ic) = qaux(is  :ie  ,j,QC)
          qtempr(1:ic,R_UN ) = q(is  :ie  ,j,QU   ) - 0.5d0*((dqx(is  :ie  ,j,2)-dqx(is  :ie  ,j,1))/q(is  :ie  ,j,QRHO))
          qtempr(1:ic,R_P  ) = q(is  :ie  ,j,QPRES) - 0.5d0*( dqx(is  :ie  ,j,1)+dqx(is  :ie  ,j,2))*cspeed(1:ic)
          qtempr(1:ic,R_UT1) = q(is  :ie  ,j,QV   ) - 0.5d0 * dqx(is  :ie  ,j,3)
          qtempr(1:ic,R_UT2) = 0.d0

          qtempr(1:ic,R_RHO) = 0.d0
          do nsp = 1,nspecies
            qtempr(1:ic,R_Y-1+nsp) = q(is  :ie  ,j,QFS-1+nsp)*q(is  :ie  ,j,QRHO) - 0.5d0*(dqx(is  :ie  ,j,4+nsp) &
                 + q(is  :ie  ,j,QFS-1+nsp)*(dqx(is  :ie  ,j,1)+dqx(is  :ie  ,j,2))/cspeed(1:ic) )
            qtempr(1:ic,R_RHO) = qtempr(1:ic,R_RHO) + qtempr(1:ic,R_Y-1+nsp)
          enddo

          do nsp = 1,nspecies
            qtempr(1:ic,5+nsp) = qtempr(1:ic,5+nsp)/qtempr(1:ic,1)
          enddo
 
          ! Small and avg c
          cav(1:ic) = HALF * ( qaux(is:ie,j,QC) + qaux(is-1:ie-1,j,QC) )
          csmall(1:ic) = min( qaux(is:ie,j,QCSML), qaux(is-1:ie-1,j,QCSML) )

         ! TODO: Make this loop a call to a vector EOS routine
          do vii = 1, ic
            ! Given rho, p, Y, evaluate T, e
            eos_state%rho      = qtempl(vii,R_RHO)
            eos_state%p        = qtempl(vii,R_P)
            eos_state%massfrac = qtempl(vii,R_Y:R_Y-1+nspecies)
            call eos_rp(eos_state)
            rhoe_l(vii) = eos_state%rho * eos_state%e
            gamc_l(vii) = eos_state%gam1
    
            eos_state%rho      = qtempr(vii,R_RHO)
            eos_state%p        = qtempr(vii,R_P)
            eos_state%massfrac = qtempr(vii,R_Y:R_Y-1+nspecies)
            call eos_rp(eos_state)
            rhoe_r(vii) = eos_state%rho * eos_state%e
            gamc_r(vii) = eos_state%gam1
          enddo
 
          ! Single point version of multi-component Riemann solve
          ! Argument order:
          ! Input Left state: rho, vel, trans vel 1, trans vel 2, pressure, reint, 1:nspecies species, gamc
          ! Input Right states: rho, vel, trans vel 1, trans vel 2, pressure, reint, 1:nspecies species, gamc,
          ! Output Godunov states: iu, trans vel 1, trans vel 2, pressure, game, regd, rgd, ustar
          ! Work arrays: eos_state, gdnv_state
          ! Array sizes: nspecies
          ! Fluxes: rho, umx, umy, ueden, ueint
          ! Tests: idir=1, coord_type=0, bc_test = 1.0, 
          ! smallc = max( csml(i), csml(i-1)
          ! cav = HALF*( c(i,j), c(i-1,j))
          ! TODO: Make this loop a call to a vector Riemann solver
          vii = 0
          do vi = is, ie ! source/dest array index space   
             vii = vii + 1 ! work array index

             if (riemann_solver .eq. 3) then
                call riemann_md_singlepoint( &
                     qtempl(vii,R_RHO), qtempl(vii,R_UN), qtempl(vii,R_UT1), qtempl(vii,R_UT2), qtempl(vii,R_P), rhoe_l(vii), qtempl(vii,R_Y:R_Y-1+nspecies), gamc_l(vii),&
                     qtempr(vii,R_RHO), qtempr(vii,R_UN), qtempr(vii,R_UT1), qtempr(vii,R_UT2), qtempr(vii,R_P), rhoe_r(vii), qtempr(vii,R_Y:R_Y-1+nspecies), gamc_r(vii),&
                     u_gd(vii), v_gd(vii), w_gd(vii), p_gd(vii), game_gd(vii), re_gd(vii), r_gd(vii), ustar(vii), eos_state, gdnv_state, nspecies,&
                     flux_tmp(vii,URHO), flux_tmp(vii,UMX), flux_tmp(vii,UMY), flux_tmp(vii,UMZ), flux_tmp(vii,UEDEN), flux_tmp(vii,UEINT), &
                     idir, coord_type, bc_test_val, csmall(vii), cav(vii) )
             else if (riemann_solver .eq. 1) then
                call riemann_cg_singlepoint( &
                     qtempl(vii,R_RHO), qtempl(vii,R_UN), qtempl(vii,R_UT1), qtempl(vii,R_UT2), qtempl(vii,R_P), rhoe_l(vii), qtempl(vii,R_Y:R_Y-1+nspecies), gamc_l(vii),&
                     qtempr(vii,R_RHO), qtempr(vii,R_UN), qtempr(vii,R_UT1), qtempr(vii,R_UT2), qtempr(vii,R_P), rhoe_r(vii), qtempr(vii,R_Y:R_Y-1+nspecies), gamc_r(vii),&
                     u_gd(vii), v_gd(vii), w_gd(vii), p_gd(vii), game_gd(vii), re_gd(vii), r_gd(vii), ustar(vii),nspecies,&
                     flux_tmp(vii,URHO), flux_tmp(vii,UMX), flux_tmp(vii,UMY), flux_tmp(vii,UMZ), flux_tmp(vii,UEDEN), flux_tmp(vii,UEINT), &
                     idir, coord_type, bc_test_val, csmall(vii), cav(vii) )
             else
                write(*,*) "Aborting, no valid Riemann sovler", riemann_solver
                call abort()
             endif

             ! Get upwinded massfractions
             if (ustar(vii) .gt. ZERO) then
                eos_state%massfrac = qtempl(vii,R_Y:R_Y+nspecies-1)
             else if (ustar(vii) .lt. ZERO) then
                eos_state%massfrac = qtempr(vii,R_Y:R_Y+nspecies-1)
             else
                eos_state%massfrac = HALF*(qtempl(vii,R_Y:R_Y+nspecies-1) + qtempr(vii,R_Y:R_Y+nspecies-1))
             endif
             eos_state%rho = r_gd(vii)
             eos_state%p   = p_gd(vii)
             call eos_rp(eos_state)
             re_gd(vii) = r_gd(vii) * eos_state%e

             flux_tmp(vii,URHO)            = r_gd(vii) * u_gd(vii)
             flux_tmp(vii,UFS:UFS+nspecies-1) = flux_tmp(vii,URHO) * eos_state%massfrac
             flux_tmp(vii,UMX)             = flux_tmp(vii,URHO) * u_gd(vii) + p_gd(vii)
             flux_tmp(vii,UMY)             = flux_tmp(vii,URHO) * v_gd(vii)
             flux_tmp(vii,UMZ)             = ZERO
             flux_tmp(vii,UEINT)           = u_gd(vii) * re_gd(vii)
             flux_tmp(vii,UEDEN)           = u_gd(vii) * ( re_gd(vii) &
                  +  HALF * r_gd(vii) * (u_gd(vii)**2 + v_gd(vii)**2 + w_gd(vii)**2) + p_gd(vii) )
             
             ! Clear unused flux slots
             flux_tmp(vii, UTEMP) = 0.0
             if (naux .gt. 0) then
                flux_tmp(vii, UFX:UFX+naux) = 0.0
             endif
             if (nadv .gt. 0) then
                flux_tmp(vii, UFA:UFA+nadv) = 0.0
             endif

          enddo ! end of loop over vi

          do ivar = 1, NVAR
             flux1(is:ie,j,ivar) = flux1(is:ie,j,ivar) + flux_tmp(1:ic,ivar) * Ax(is:ie,j)
          enddo
       enddo
    enddo

#ifdef PELEC_USE_EB
    call slopey(q,flatn,qd_lo,qd_hi, &
         dqy,qt_lo,qt_hi, &
         lo(1)-nextra,lo(2)-nextra,  &
         hi(1)+nextra,hi(2)+nextra,QVAR,NQAUX,&
         domlo,domhi, &
         qaux, qa_lo, qa_hi, &
         flag, fglo, fghi)
#else
    call slopey(q,flatn,qd_lo,qd_hi, &
         dqy,qt_lo,qt_hi, &
         lo(1),lo(2),hi(1),hi(2),QVAR)
#endif

    do j = lo(2)-nextra+1, hi(2)+nextra
       do i = lo(1)-nextra, hi(1)+nextra, VECLEN
          is = i
          ie = min(is+VECLEN-1, hi(1)+nextra)
          ic = ie - is + 1

          cspeed(1:ic) = qaux(is:ie,j-1,QC)
          qtempl(1:ic,R_UN ) = q(is:ie,j-1,QV   ) + 0.5d0*((dqy(is:ie,j-1,2)-dqy(is:ie,j-1,1))/q(is:ie,j-1,QRHO))
          qtempl(1:ic,R_P  ) = q(is:ie,j-1,QPRES) + 0.5d0*( dqy(is:ie,j-1,1)+dqy(is:ie,j-1,2))*cspeed(1:ic)
          qtempl(1:ic,R_UT1) = q(is:ie,j-1,QU   ) + 0.5d0 * dqy(is:ie,j-1,3)
          qtempl(1:ic,R_UT2) = 0.d0

          qtempl(1:ic,R_RHO) = 0.d0
          do nsp = 1,nspecies
             qtempl(1:ic,R_Y-1+nsp) = q(is:ie,j-1,QFS-1+nsp)*q(is:ie,j-1,QRHO) + 0.5d0*(dqy(is:ie,j-1,4+nsp) &
                  + q(is:ie,j-1,QFS-1+nsp)*(dqy(is:ie,j-1,1)+dqy(is:ie,j-1,2))/cspeed(1:ic) )
             qtempl(1:ic,R_RHO) = qtempl(1:ic,R_RHO) + qtempl(1:ic,R_Y-1+nsp)
          enddo

          do nsp = 1,nspecies
             qtempl(1:ic,R_Y-1+nsp) = qtempl(1:ic,R_Y-1+nsp)/qtempl(1:ic,R_RHO)
          enddo

          cspeed(1:ic) = qaux(is:ie,j  ,QC)
          qtempr(1:ic,R_UN ) = q(is:ie,j  ,QV   ) - 0.5d0*((dqy(is:ie,j  ,2)-dqy(is:ie,j  ,1))/q(is:ie,j  ,QRHO))
          qtempr(1:ic,R_P  ) = q(is:ie,j  ,QPRES) - 0.5d0*( dqy(is:ie,j  ,1)+dqy(is:ie,j  ,2))*cspeed(1:ic)
          qtempr(1:ic,R_UT1) = q(is:ie,j  ,QU   ) - 0.5d0 * dqy(is:ie,j,3)
          qtempr(1:ic,R_UT2) = 0.d0

          qtempr(1:ic,R_RHO) = 0.d0
          do nsp = 1,nspecies
             qtempr(1:ic,R_Y-1+nsp) = q(is:ie,j  ,QFS-1+nsp)*q(is:ie,j  ,QRHO) - 0.5d0*(dqy(is:ie,j  ,4+nsp) &
                  + q(is:ie,j  ,QFS-1+nsp)*(dqy(is:ie,j  ,1)+dqy(is:ie,j  ,2))/cspeed(1:ic) )
             qtempr(1:ic,R_RHO) = qtempr(1:ic,R_RHO) + qtempr(1:ic,R_Y-1+nsp)
          enddo
              
          do nsp = 1,nspecies
             qtempr(1:ic,R_Y-1+nsp) = qtempr(1:ic,R_Y-1+nsp)/qtempr(1:ic,R_RHO)
          enddo

          ! Small and avg c
          cav(1:ic) = HALF * ( qaux(is:ie,j,QC) + qaux(is:ie,j-1,QC) )
          csmall(1:ic) = min( qaux(is:ie,j,QCSML), qaux(is:ie,j-1,QCSML) )

          ! TODO: Make this loop a call to a vector EOS routine
          do vii = 1, ic
             ! Given rho, p, Y, evaluate T, e
             eos_state%rho      = qtempl(vii,R_RHO)
             eos_state%p        = qtempl(vii,R_P)
             eos_state%massfrac = qtempl(vii,R_Y:R_Y-1+nspecies)
             call eos_rp(eos_state)
             rhoe_l(vii) = eos_state%rho * eos_state%e
             gamc_l(vii) = eos_state%gam1
        
             eos_state%rho      = qtempr(vii,R_RHO)
             eos_state%p        = qtempr(vii,R_P)
             eos_state%massfrac = qtempr(vii,R_Y:R_Y-1+nspecies)
             call eos_rp(eos_state)
             rhoe_r(vii) = eos_state%rho * eos_state%e
             gamc_r(vii) = eos_state%gam1
          enddo

          ! Single point version of multi-component Riemann solve
          ! Argument order:
          ! Input Left state: rho, vel, trans vel 1, trans vel 2, pressure, reint, 1:nspecies species, gamc
          ! Input Right states: rho, vel, trans vel 1, trans vel 2, pressure, reint, 1:nspecies species, gamc,
          ! Output Godunov states: iu, trans vel 1, trans vel 2, pressure, game, regd, rgd, ustar
          ! Work arrays: eos_state, gdnv_state
          ! Array sizes: nspecies
          ! Fluxes: rho, umx, umy, ueden, ueint
          ! Tests: idir=1, coord_type=0, bc_test = 1.0, 
          ! smallc = max( csml(i), csml(i-1)
          ! cav = HALF*( c(i,j), c(i-1,j))
    
          ! TODO: Make this loop a call to a vector Riemann solver
          vii = 0
          do vi = is, ie ! source/dest array index space
             vii = vii + 1 ! work array index

             if (riemann_solver .eq. 3) then
                call riemann_md_singlepoint( &
                     qtempl(vii,R_RHO), qtempl(vii,R_UN), qtempl(vii,R_UT1), qtempl(vii,R_UT2), qtempl(vii,R_P), rhoe_l(vii), qtempl(vii,R_Y:R_Y-1+nspecies), gamc_l(vii),&
                     qtempr(vii,R_RHO), qtempr(vii,R_UN), qtempr(vii,R_UT1), qtempr(vii,R_UT2), qtempr(vii,R_P), rhoe_r(vii), qtempr(vii,R_Y:R_Y-1+nspecies), gamc_r(vii),&
                     v_gd(vii), u_gd(vii), w_gd(vii), p_gd(vii), game_gd(vii), re_gd(vii), r_gd(vii), ustar(vii), eos_state, gdnv_state, nspecies,&
                     flux_tmp(vii,URHO), flux_tmp(vii,UMY), flux_tmp(vii,UMX), flux_tmp(vii,UMZ), flux_tmp(vii,UEDEN), flux_tmp(vii,UEINT), &
                     idir, coord_type, bc_test_val, csmall(vii), cav(vii) )
             else if (riemann_solver .eq. 1) then
                call riemann_cg_singlepoint( &
                     qtempl(vii,R_RHO), qtempl(vii,R_UN), qtempl(vii,R_UT1), qtempl(vii,R_UT2), qtempl(vii,R_P), rhoe_l(vii), qtempl(vii,R_Y:R_Y-1+nspecies), gamc_l(vii),&
                     qtempr(vii,R_RHO), qtempr(vii,R_UN), qtempr(vii,R_UT1), qtempr(vii,R_UT2), qtempr(vii,R_P), rhoe_r(vii), qtempr(vii,R_Y:R_Y-1+nspecies), gamc_r(vii),&
                     v_gd(vii), u_gd(vii), w_gd(vii), p_gd(vii), game_gd(vii), re_gd(vii), r_gd(vii), ustar(vii), nspecies,&
                     flux_tmp(vii,URHO), flux_tmp(vii,UMY), flux_tmp(vii,UMX), flux_tmp(vii,UMZ), flux_tmp(vii,UEDEN), flux_tmp(vii,UEINT), &
                     idir, coord_type, bc_test_val, csmall(vii), cav(vii) )
             else
                call abort()
             endif

             ! Get upwinded massfractions
             if (ustar(vii) .gt. ZERO) then
                eos_state%massfrac = qtempl(vii,R_Y:R_Y+nspecies-1)
             else if (ustar(vii) .lt. ZERO) then
                eos_state%massfrac = qtempr(vii,R_Y:R_Y+nspecies-1)
             else
                eos_state%massfrac = HALF*(qtempl(vii,R_Y:R_Y+nspecies-1) + qtempr(vii,R_Y:R_Y+nspecies-1))
             endif
             eos_state%rho = r_gd(vii)
             eos_state%p   = p_gd(vii)
             call eos_rp(eos_state)
             re_gd(vii) = r_gd(vii) * eos_state%e

             flux_tmp(vii,URHO)            = r_gd(vii) * v_gd(vii)
             flux_tmp(vii,UFS:UFS+nspecies-1) = flux_tmp(vii,URHO) * eos_state%massfrac
             flux_tmp(vii,UMX)             = flux_tmp(vii,URHO) * u_gd(vii)
             flux_tmp(vii,UMY)             = flux_tmp(vii,URHO) * v_gd(vii) + p_gd(vii)
             flux_tmp(vii,UMZ)             = ZERO
             flux_tmp(vii,UEINT)           = v_gd(vii) * re_gd(vii)
             flux_tmp(vii,UEDEN)           = v_gd(vii) * ( re_gd(vii) &
                  +  HALF * r_gd(vii) * (u_gd(vii)**2 + v_gd(vii)**2 + w_gd(vii)**2) + p_gd(vii) )

             ! Clear unused flux slots
             flux_tmp(vii, UTEMP) = 0.0
             if (naux .gt. 0) then
                flux_tmp(vii, UFX:UFX+naux) = 0.0
             endif
             if (nadv .gt. 0) then
                flux_tmp(vii, UFA:UFA+nadv) = 0.0
             endif
          enddo

          do ivar = 1, NVAR
             flux2(is:ie,j, ivar) = flux2(is:ie,j,ivar) + flux_tmp(1:ic,ivar) * Ay(is:ie,j)
          enddo
       enddo
    enddo

    ! Done computing flux through regular faces.

    ! TODO: Flux through wall face here...
#ifdef PELEC_USE_EB

    full_area = h(1)**(dim - 1)

    ! Loop over cut cells only - need to pass in a list of these (lift the list Marc built for diffusion)
    do L = 0, nebflux-1, VECLEN
       is = L
       ie = min(is+VECLEN-1, nebflux-1)
       ic = ie - is +1

       do vii = 1, ic
          i = ebg(is+vii-1) % iv(0)
          j = ebg(is+vii-1) % iv(1)
          if (i.ge.lo(1)-nextra .and. i.le.hi(1)+nextra &
               .and. j.ge.lo(2)-nextra .and. j.le.hi(2)+nextra ) then

             ! Get normal - checkto make sure this is facing outward!
             eb_norm = ebg(is+vii-1)%eb_normal
             eb_norm = eb_norm / sqrt(eb_norm(1)**2 + eb_norm(2)**2)

             ! Assume left state is the cell centered state - normal veclocity
             cspeed(vii) = qaux(i,j,QC)

             qtempl(vii,R_UN ) =  - q(i,j,QU)*eb_norm(1) &
                  -                 q(i,j,QV)*eb_norm(2)
             qtempl(vii,R_UT1) = 0.d0
             qtempl(vii,R_UT2) = 0.d0
             qtempl(vii,R_P  ) = q(i,j,QPRES)

             qtempl(vii,R_RHO) = q(i,j,QRHO)

             do nsp = 1,nspecies
                qtempl(vii,R_Y-1+nsp) = q(i,j,QFS-1+nsp)
             enddo

             ! Flip the velocity about the normal for the right state - will use left
             ! state for remainder of right state
             qtempr(vii,R_UN) = -1.0*qtempl(vii,R_UN)
          
             ! Small and avg c
             cav(vii) =  qaux(i,j,QC) 
             csmall(vii) = qaux(i,j,QCSML)

          endif
       enddo

       !     TODO: Make this loop a call to a vector EOS routine
       do vii = 1, ic
          ! Given rho, p, Y, evaluate T, e
          i = ebg(is+vii-1) % iv(0)
          j = ebg(is+vii-1) % iv(1)
          if (i.ge.lo(1)-nextra .and. i.le.hi(1)+nextra &
               .and. j.ge.lo(2)-nextra .and. j.le.hi(2)+nextra ) then
                eos_state%rho      = qtempl(vii,R_RHO)
                eos_state%p        = qtempl(vii,R_P)
                eos_state%massfrac = qtempl(vii,R_Y:R_Y-1+nspecies)
                call eos_rp(eos_state)
                rhoe_l(vii) = eos_state%rho * eos_state%e
                gamc_l(vii) = eos_state%gam1
          endif

       enddo
       
       ! Solve Riemann problem; store flux in flux4 - cell centered data structure
       vii = 0
       do vi = is, ie ! source/dest array index space
          vii = vii + 1 ! work array index

          i = ebg(is+vii-1) % iv(0)
          j = ebg(is+vii-1) % iv(1)
          ! this is going to prevent vectorization, but until vector riemann solver
          ! is in place might as well avoid copying into a temporary 
          if (i.ge.lo(1)-nextra .and. i.le.hi(1)+nextra &
               .and. j.ge.lo(2)-nextra .and. j.le.hi(2)+nextra ) then

             if (riemann_solver .eq. 3) then
                call riemann_md_singlepoint( &
                     qtempl(vii,R_RHO), qtempl(vii,R_UN), qtempl(vii,R_UT1), qtempl(vii,R_UT2), qtempl(vii,R_P), rhoe_l(vii), qtempl(vii,R_Y:R_Y-1+nspecies), gamc_l(vii),&
                     qtempl(vii,R_RHO), qtempr(vii,R_UN), qtempl(vii,R_UT1), qtempl(vii,R_UT2), qtempl(vii,R_P), rhoe_l(vii), qtempl(vii,R_Y:R_Y-1+nspecies), gamc_l(vii),&
                     u_gd(vii), v_gd(vii), w_gd(vii), p_gd(vii), game_gd(vii), re_gd(vii), r_gd(vii), ustar(vii), eos_state, gdnv_state, nspecies,&
                     flux_tmp(vii,URHO), tflux(1), tflux(2), tflux(3),  flux_tmp(vii,UEDEN), flux_tmp(vii,UEINT), &
                     idir, coord_type, bc_test_val, csmall(vii), cav(vii) )
             else if (riemann_solver .eq. 1) then
                call riemann_cg_singlepoint( &
                     qtempl(vii,R_RHO), qtempl(vii,R_UN), qtempl(vii,R_UT1), qtempl(vii,R_UT2), qtempl(vii,R_P), rhoe_l(vii), qtempl(vii,R_Y:R_Y-1+nspecies), gamc_l(vii),&
                     qtempl(vii,R_RHO), qtempr(vii,R_UN), qtempl(vii,R_UT1), qtempl(vii,R_UT2), qtempl(vii,R_P), rhoe_l(vii), qtempl(vii,R_Y:R_Y-1+nspecies), gamc_l(vii),&
                     u_gd(vii), v_gd(vii), w_gd(vii), p_gd(vii), game_gd(vii), re_gd(vii), r_gd(vii), ustar(vii), nspecies,&
                     flux_tmp(vii,URHO), tflux(1), tflux(2), tflux(3), flux_tmp(vii,UEDEN), flux_tmp(vii,UEINT), &
                     idir, coord_type, bc_test_val, csmall(vii), cav(vii) )
             else
                call abort()
             endif

             ! Get upwinded massfractions
             eos_state%massfrac = qtempl(vii,R_Y:R_Y+nspecies-1)
             eos_state%rho = r_gd(vii)
             eos_state%p   = p_gd(vii)
             call eos_rp(eos_state)
             re_gd(vii) = r_gd(vii) * eos_state%e

             flux_tmp(vii,URHO)            = r_gd(vii) * u_gd(vii)
             flux_tmp(vii,UFS:UFS+nspecies-1) = flux_tmp(vii,URHO) * eos_state%massfrac

             eb_norm = ebg(is+vii-1)%eb_normal
             eb_norm = eb_norm / sqrt(eb_norm(1)**2 + eb_norm(2)**2)
             flux_tmp(vii,UMX) = -tflux(1) * eb_norm(1)
             flux_tmp(vii,UMY) = -tflux(1) * eb_norm(2)
             flux_tmp(vii,UMZ)             = ZERO
             flux_tmp(vii,UEINT)           = u_gd(vii) * re_gd(vii)
             flux_tmp(vii,UEDEN)           = u_gd(vii) * ( re_gd(vii) &
                  +  HALF * r_gd(vii) * (u_gd(vii)**2 + v_gd(vii)**2 + w_gd(vii)**2) + p_gd(vii) )

             !   Compute species flux like passive scalar from intermediate state
             do nsp = 0, nspecies-1
                flux_tmp(vii,UFS+nsp) = flux_tmp(vii,URHO)*qtempl(vii,R_Y+nsp)
             enddo
          endif

       enddo ! End future vector loop

       ! Copy result into ebflux vector (should this be joined with above loop?)
       vii = 0
       do vi = is, ie ! source/dest array index space
          vii = vii + 1 ! work array index

          i = ebg(is+vii-1) % iv(0)
          j = ebg(is+vii-1) % iv(1)
          ! this is going to prevent vectorization, but until vector riemann solver
          ! is in place might as well avoid copying into a temporary 
          if (i.ge.lo(1)-nextra .and. i.le.hi(1)+nextra &
               .and. j.ge.lo(2)-nextra .and. j.le.hi(2)+nextra ) then

             do ivar = 1, NVAR
                ebflux(vi,ivar) = ebflux(vi,ivar) + flux_tmp(vii,ivar) * ebg(is+vii-1)%eb_area * full_area
             enddo
             
          endif
       enddo ! End future vector loop
    enddo ! End loop over cut cells
#endif

    ! Deallocate arrays

    call bl_deallocate (dqx)
    call bl_deallocate (dqy)
    call destroy(eos_state)
    call destroy(gdnv_state)

    ! Flux1(i,j,N) has edge based flux of rho,u,v,w,eden,eint,species indexed by
    ! URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS:UFS+nspecies-1, and similar for flux2, flux3
    ! EB face flux assuming wall BC is loaded into ebflux with same ordering

    ! After this, must call routine to multiply fluxes by aperatures and interpolate to face centers 
    do ivar=1,NVAR
       do j = lo(2)-nextra, hi(2)+nextra
          do i = lo(1)-nextra, hi(1)+nextra
             D(i,j,ivar) = - (flux1(i+1,j,ivar) - flux1(i,j,ivar) &
                  +           flux2(i,j+1,ivar) - flux2(i,j,ivar) )/V(i,j)
          enddo
       enddo
    enddo

  end subroutine pc_hyp_mol_flux
end module hyp_advection_module 
