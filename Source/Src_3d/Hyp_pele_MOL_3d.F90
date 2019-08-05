module hyp_advection_module 

  use amrex_ebcellflag_module, only : get_neighbor_cells
  use pelec_eb_stencil_types_module, only : eb_bndry_geom
  use riemann_util_module, only : riemann_md_singlepoint, riemann_md_vec
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
                     Az,  Azlo,  Azhi,&
                     flux3, fd3_lo, fd3_hi, &
                     flatn, fltd_lo, fltd_hi, &
                     V, Vlo, Vhi, &
                     D, Dlo, Dhi,&
#ifdef PELEC_USE_EB
                     vfrac, vflo, vfhi, &
                     flag, fglo, fghi, &
                     ebg, Nebg, ebflux, nebflux, &
#endif
                     h) &
                     bind(C,name="pc_hyp_mol_flux")



    use amrex_mempool_module, only : bl_allocate, bl_deallocate
    use meth_params_module, only : QVAR, NVAR, QPRES, QRHO, QU, QV, QW, &
                                   QFS,  &
                                   QC, QCSML, NQAUX, nadv, &
                                   URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, UTEMP, UFX, UFA, &
                                   eb_small_vfrac
    use slope_module, only : slopex, slopey, slopez
    use network, only : nspecies, naux
    use eos_type_module
    use eos_module, only : eos_t, eos_rp
    use riemann_module, only: cmpflx, shock
    use amrex_constants_module
    use amrex_fort_module, only : amrex_real

    implicit none

    integer, parameter :: VECLEN = 16
    integer :: vis, vie, vic ! Loop bounds for vector blocking
    integer :: vi, vii ! Loop indicies for unrolled loops over 

    integer, intent(in) ::      qd_lo(3),   qd_hi(3)
    integer, intent(in) ::      qa_lo(3),   qa_hi(3)
    integer, intent(in) ::         lo(3),      hi(3)
    integer, intent(in) ::      domlo(3),   domhi(3)
    integer, intent(in) ::       Axlo(3),    Axhi(3)
    integer, intent(in) ::     fd1_lo(3),  fd1_hi(3)
    integer, intent(in) ::       Aylo(3),    Ayhi(3)
    integer, intent(in) ::     fd2_lo(3),  fd2_hi(3)
    integer, intent(in) ::       Azlo(3),    Azhi(3)
    integer, intent(in) ::     fd3_lo(3),  fd3_hi(3)
    integer, intent(in) ::    fltd_lo(3), fltd_hi(3)
    integer, intent(in) ::        Vlo(3),     Vhi(3)
    integer, intent(in) ::        Dlo(3),     Dhi(3)
    double precision, intent(in) :: h(3)

#ifdef PELEC_USE_EB
    integer, intent(in) ::  fglo(3),    fghi(3)
    integer, intent(in) ::  vflo(3),    vfhi(3)
    integer, intent(in) :: flag(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3))
    real(amrex_real), intent(in) :: vfrac(vflo(1):vfhi(1),vflo(2):vfhi(2),vflo(3):vfhi(3))

    integer, intent(in) :: nebflux
    real(amrex_real), intent(inout) ::   ebflux(0:nebflux-1,1:NVAR)
    integer,            intent(in   ) :: Nebg
    type(eb_bndry_geom),intent(in   ) :: ebg(0:Nebg-1)    
    real(amrex_real) :: eb_norm(3), full_area
    real(amrex_real) :: sum_kappa, sum_nbrs
#endif
    double precision, intent(in) ::     q(  qd_lo(1):  qd_hi(1),  qd_lo(2):  qd_hi(2),  qd_lo(3):  qd_hi(3),QVAR)  !> State
    double precision, intent(in) ::  qaux(  qa_lo(1):  qa_hi(1),  qa_lo(2):  qa_hi(2),  qa_lo(3):  qa_hi(3),NQAUX) !> Auxiliary state
    double precision, intent(in) :: flatn(fltd_lo(1):fltd_hi(1),fltd_lo(2):fltd_hi(2),fltd_lo(3):fltd_hi(3))

    double precision, intent(in   ) ::    Ax(  Axlo(1):  Axhi(1),  Axlo(2):  Axhi(2),  Axlo(3):  Axhi(3))
    double precision, intent(inout) :: flux1(fd1_lo(1):fd1_hi(1),fd1_lo(2):fd1_hi(2),fd1_lo(3):fd1_hi(3),NVAR)
    double precision, intent(in   ) ::    Ay(  Aylo(1):  Ayhi(1),  Aylo(2):  Ayhi(2),  Aylo(3):  Ayhi(3))
    double precision, intent(inout) :: flux2(fd2_lo(1):fd2_hi(1),fd2_lo(2):fd2_hi(2),fd2_lo(3):fd2_hi(3),NVAR)
    double precision, intent(in   ) ::    Az(  Azlo(1):  Azhi(1),  Azlo(2):  Azhi(2),  Azlo(3):  Azhi(3))
    double precision, intent(inout) :: flux3(fd3_lo(1):fd3_hi(1),fd3_lo(2):fd3_hi(2),fd3_lo(3):fd3_hi(3),NVAR)
    double precision, intent(inout) ::     V(   Vlo(1):   Vhi(1),   Vlo(2):   Vhi(2),   Vlo(3):   Vhi(3))
    double precision, intent(inout) ::     D(   Dlo(1):   Dhi(1),   Dlo(2):   Dhi(2),   Dlo(3):   Dhi(3),NVAR)

    integer :: i, j, k, nsp, L, ivar
    integer :: qt_lo(3), qt_hi(3)

    ! Left and right state arrays (edge centered, cell centered)
    double precision, pointer :: dqx(:,:,:,:), dqy(:,:,:,:), dqz(:,:,:,:)

    ! Other left and right state arrays
    double precision :: qtempl(VECLEN,1:5+nspecies)
    double precision :: qtempr(VECLEN,1:5+nspecies)
    double precision :: rhoe_l(VECLEN)
    double precision :: rhoe_r(VECLEN)
    double precision :: cspeed(VECLEN)
    double precision :: gamc_l(VECLEN)
    double precision :: gamc_r(VECLEN)
    double precision :: cavg(VECLEN)
    double precision :: csmall(VECLEN)

    ! Scratch for neighborhood of cut cells
    integer :: nbr(-1:1,-1:1,-1:1)

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

   !initialize flux_tmp to 0
   !don't want fortran to fill it with wrong values
    flux_tmp = 0.d0

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

    !  allocate spaece for slopes
    call bl_allocate ( dqx, qt_lo(1), qt_hi(1), qt_lo(2), qt_hi(2), qt_lo(3), qt_hi(3), 1, QVAR)
    call bl_allocate ( dqy, qt_lo(1), qt_hi(1), qt_lo(2), qt_hi(2), qt_lo(3), qt_hi(3), 1, QVAR)
    call bl_allocate ( dqz, qt_lo(1), qt_hi(1), qt_lo(2), qt_hi(2), qt_lo(3), qt_hi(3), 1, QVAR)

    call bl_proffortfuncstart_int(0)
    ! Compute all slopes at kc (k3d)
#ifdef PELEC_USE_EB
    call slopex(q,flatn,qd_lo,qd_hi, &
                   dqx,qt_lo,qt_hi, &
                   lo(1)-nextra,lo(2)-nextra,lo(3)-nextra,  &
                   hi(1)+nextra,hi(2)+nextra,hi(3)+nextra,QVAR,NQAUX, &
                   domlo,domhi, &
                   qaux, qa_lo, qa_hi, &
                   flag, fglo, fghi)
#else
    call slopex(q,flatn,qd_lo,qd_hi, &
                   dqx,qt_lo,qt_hi, &
                   lo(1),lo(2),lo(3),hi(1),hi(2),hi(3),QVAR,NQAUX &
                   qaux, qa_lo, qa_hi)
#endif
    call bl_proffortfuncstop_int(0)
      !  alphas   1,2,3 correspond to u-c,  u+c  repeated u  repsectively, 
      !  3 and 4 are transverse velocities
      !  species alphas are for rho Y_k
      ! right eigenvectors are rho, u, p, v, w
      ! i-1 => vis-1:vie-1
      ! i => vis:vie
      ! noidx => 1:vic
      !   in qtemp, 1 is rho, 2 is u, 3 is v , 4 is w and 5 is p,
      !   species after  . . note predicted rho is sum of predicted rho Y_k

    call bl_proffortfuncstart_int(1)
    do k = lo(3)-nextra, hi(3)+nextra
       do j = lo(2)-nextra, hi(2)+nextra
          do i = lo(1)-nextra+1, hi(1)+nextra, VECLEN
             vis = i
             vie = min(vis+VECLEN-1, hi(1)+nextra)
             vic = vie - vis + 1

             cspeed(1:vic) = qaux(vis-1:vie-1,j,k,QC)

             ! Left u
             qtempl(1:vic,R_UN ) = q(vis-1:vie-1,j,k,QU   ) + 0.5d0*((dqx(vis-1:vie-1,j,k,2)-dqx(vis-1:vie-1,j,k,1))/q(vis-1:vie-1,j,k,QRHO))

             ! Left p
             qtempl(1:vic,R_P  ) = q(vis-1:vie-1,j,k,QPRES) + 0.5d0*( dqx(vis-1:vie-1,j,k,1)+dqx(vis-1:vie-1,j,k,2))*cspeed(1:vic)

             ! Left v
             qtempl(1:vic,R_UT1) = q(vis-1:vie-1,j,k,QV   ) + 0.5d0*  dqx(vis-1:vie-1,j,k,3)

             ! Left w
             qtempl(1:vic,R_UT2) = q(vis-1:vie-1,j,k,QW   ) + 0.5d0*  dqx(vis-1:vie-1,j,k,4)

             ! Left rho - computed as sum(rhoY_k) below after species
             qtempl(1:vic,R_RHO) = 0.d0
             do nsp = 1,nspecies

                qtempl(1:vic,R_Y - 1 +nsp) = q(vis-1:vie-1,j,k,QFS-1+nsp) * q(vis-1:vie-1,j,k,QRHO) + 0.5d0*(dqx(vis-1:vie-1,j,k,4+nsp) &
                     + q(vis-1:vie-1,j,k,QFS-1+nsp) * (dqx(vis-1:vie-1,j,k,1) + dqx(vis-1:vie-1,j,k,2))/cspeed(1:vic) )
                qtempl(1:vic,R_RHO) = qtempl(1:vic,R_RHO) + qtempl(1:vic,R_Y - 1 +nsp)

             enddo

             do nsp = 1,nspecies
               qtempl(1:vic,R_Y -1 +nsp) = qtempl(1:vic,R_Y - 1 +nsp)/qtempl(1:vic,R_RHO)
             enddo

             cspeed(1:vic) = qaux(vis:vie,j,k,QC  )

             ! Right u
             qtempr(1:vic,R_UN ) = q(vis:vie,j,k,QU   ) - 0.5d0 * ((dqx(vis:vie,j,k,2)-dqx(vis:vie,j,k,1))/q(vis:vie,j,k,QRHO))

             ! Right p
             qtempr(1:vic,R_P  ) = q(vis:vie,j,k,QPRES) - 0.5d0 * ( dqx(vis:vie,j,k,1)+dqx(vis:vie,j,k,2))*cspeed(1:vic)

             ! Right v
             qtempr(1:vic,R_UT1) = q(vis:vie,j,k,QV   ) - 0.5d0 *   dqx(vis:vie,j,k,3)

             ! Right w
             qtempr(1:vic,R_UT2) = q(vis:vie,j,k,QW   ) - 0.5d0 *   dqx(vis:vie,j,k,4)

             ! Right rho - computed as sum(rhoY_k) below after species
             qtempr(1:vic,R_RHO) = 0.d0

             do nsp = 1,nspecies

                qtempr(1:vic,R_Y - 1 +nsp) = q(vis:vie,j,k,QFS-1+nsp)*q(vis:vie,j,k,QRHO) - 0.5d0*(dqx(vis:vie,j,k,4+nsp) &
                     + q(vis:vie,j,k,QFS-1+nsp)*(dqx(vis:vie,j,k,1)+dqx(vis:vie,j,k,2))/cspeed(1:vic) )
                qtempr(1:vic,R_RHO) = qtempr(1:vic,R_RHO) + qtempr(1:vic,R_Y - 1 + nsp)

             enddo

             do nsp = 1,nspecies
                qtempr(1:vic,R_Y - 1 +nsp) = qtempr(1:vic,R_Y - 1 +nsp)/qtempr(1:vic,R_RHO)
             enddo

             ! Small and avg c
             cavg(1:vic) = HALF * ( qaux(vis:vie,j,k,QC) + qaux(vis-1:vie-1,j,k,QC) )
             csmall(1:vic) = min( qaux(vis:vie,j,k,QCSML), qaux(vis-1:vie-1,j,k,QCSML) )

             ! TODO: Make this loop a call to a vector EOS routine
             do vii = 1, vic
                ! Have p, rhoY (composition is rhoY), rho 
                !  - evaluate T, use that to evaluate internal energy
                eos_state%rho = qtempl(vii,R_RHO)
                eos_state%p = qtempl(vii,R_P)
                eos_state%massfrac = qtempl(vii,R_Y:R_Y-1+nspecies)
                !dir$ inline recursive
                call eos_rp(eos_state)
                rhoe_l(vii) = eos_state%rho * eos_state%e
                gamc_l(vii) = eos_state%gam1

                eos_state%rho = qtempr(vii,R_RHO)
                eos_state%p = qtempr(vii,R_P)
                eos_state%massfrac = qtempr(vii,R_Y:R_Y-1+nspecies)
                !dir$ inline recursive
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
             call riemann_md_vec( &
                  qtempl(1:vic,R_RHO), qtempl(1:vic,R_UN), qtempl(1:vic,R_UT1), qtempl(1:vic,R_UT2), qtempl(1:vic,R_P), rhoe_l(1:vic), qtempl(1:vic,R_Y:R_Y-1+nspecies), gamc_l(1:vic),&
                  qtempr(1:vic,R_RHO), qtempr(1:vic,R_UN), qtempr(1:vic,R_UT1), qtempr(1:vic,R_UT2), qtempr(1:vic,R_P), rhoe_r(1:vic), qtempr(1:vic,R_Y:R_Y-1+nspecies), gamc_r(1:vic),&
                  u_gd(1:vic), v_gd(1:vic), w_gd(1:vic), p_gd(1:vic), game_gd(1:vic), re_gd(1:vic), r_gd(1:vic), ustar(1:vic),&
                  eos_state, nspecies,&
                  flux_tmp(1:vic,URHO), flux_tmp(1:vic,UMX), flux_tmp(1:vic,UMY), flux_tmp(1:vic,UMZ), flux_tmp(1:vic,UEDEN), flux_tmp(1:vic,UEINT), &
                  bc_test_val, csmall(1:vic), cavg(1:vic), vic )

             vii = 0
             do vi = vis, vie ! source/dest array index space   
                vii = vii + 1 ! work array index

                ! Compute species flux like passive scalar from intermediate state
                do nsp = 0, nspecies-1
                   if (ustar(vii) .gt. ZERO) then
                      flux_tmp(vii,UFS+nsp) = flux_tmp(vii,URHO)*qtempl(vii,R_Y+nsp)
                   else if (ustar(vii) .lt. ZERO) then
                      flux_tmp(vii,UFS+nsp) = flux_tmp(vii,URHO)*qtempr(vii,R_Y+nsp)
                   else
                      flux_tmp(vii,UFS+nsp) = flux_tmp(vii,URHO)&
                           *HALF*(qtempl(vii,R_Y+nsp) + qtempr(vii,R_Y+nsp) )
                   endif
                enddo

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
                flux1(vis:vie,j,k, ivar) = flux1(vis:vie,j,k, ivar ) + flux_tmp(1:vic,ivar) * Ax(vis:vie,j,k)

             enddo
          enddo
       enddo
    enddo
    call bl_proffortfuncstop_int(1)

    call bl_proffortfuncstart_int(2)
#ifdef PELEC_USE_EB
    call slopey(q,flatn,qd_lo,qd_hi, &
         dqy,qt_lo,qt_hi, &
         lo(1)-nextra,lo(2)-nextra,lo(3)-nextra,  &
         hi(1)+nextra,hi(2)+nextra,hi(3)+nextra,QVAR,NQAUX,&
         domlo,domhi, &
         qaux, qa_lo, qa_hi, &
         flag, fglo, fghi)
#else
    call slopey(q,flatn,qd_lo,qd_hi, &
         dqy,qt_lo,qt_hi, &
         lo(1),lo(2),lo(3),hi(1),hi(2),hi(3),QVAR,NQAUX, &
         qaux, qa_lo, qa_hi)
#endif
    call bl_proffortfuncstop_int(2)

    call bl_proffortfuncstart_int(3)
    do k = lo(3)-nextra, hi(3)+nextra
       do j = lo(2)-nextra+1, hi(2)+nextra

          do i = lo(1)-nextra, hi(1)+nextra, VECLEN
             vis = i
             vie = min(vis+VECLEN-1, hi(1)+nextra)
             vic = vie - vis + 1


             !     1,2,3 correspond to u-c, u, u+c repsectively
             ! right eigenvectors are rho, v, p, u, w
             !   in qtemp, 1 is rho, 2 is v, 3 is q , 4 is u and 5 is w


              cspeed(1:vic      ) = qaux(vis:vie,j-1,k,QC)

              ! Left v
              qtempl(1:vic,R_UN ) = q(vis:vie,j-1,k,QV   ) + 0.5d0*((dqy(vis:vie,j-1,k,2)-dqy(vis:vie,j-1,k,1))/q(vis:vie,j-1,k,QRHO))
     
              ! Left p
              qtempl(1:vic,R_P  ) = q(vis:vie,j-1,k,QPRES) + 0.5d0* (dqy(vis:vie,j-1,k,1)+dqy(vis:vie,j-1,k,2))*cspeed(1:vic)

              ! Left u
              qtempl(1:vic,R_UT1) = q(vis:vie,j-1,k,QU   ) + 0.5d0*  dqy(vis:vie,j-1,k,3)

              ! Left w
              qtempl(1:vic,R_UT2) = q(vis:vie,j-1,k,QW   ) + 0.5d0*  dqy(vis:vie,j-1,k,4)

              ! Left rho - computed as sum(rhoY_k) below after species
              qtempl(1:vic,R_RHO) = 0.d0
              do nsp = 1,nspecies

                qtempl(1:vic,R_Y-1+nsp) = q(vis:vie,j-1,k,QFS-1+nsp) * q(vis:vie,j-1,k,QRHO) + 0.5d0*(dqy(vis:vie,j-1,k,4+nsp) &
                                          + q(vis:vie,j-1,k,QFS-1+nsp) * (dqy(vis:vie,j-1,k,1) + dqy(vis:vie,j-1,k,2))/cspeed(1:vic) )
                qtempl(1:vic,R_RHO) = qtempl(1:vic,R_RHO) + qtempl(1:vic,R_Y-1+nsp)

               enddo

              do nsp = 1,nspecies
                qtempl(1:vic,R_Y-1+nsp) = qtempl(1:vic,R_Y-1+nsp)/qtempl(1:vic,R_RHO)
              enddo

              cspeed(1:vic) = qaux(vis:vie,j,k,QC)

              ! Right v
              qtempr(1:vic,R_UN ) = q(vis:vie,j,k,QV   ) - 0.5d0 * ((dqy(vis:vie,j,k,2)-dqy(vis:vie,j,k,1))/q(vis:vie,j,k,QRHO))

              ! Right p
              qtempr(1:vic,R_P  ) = q(vis:vie,j,k,QPRES) - 0.5d0 * ( dqy(vis:vie,j,k,1)+dqy(vis:vie,j,k,2))*cspeed(1:vic)

              ! Right u
              qtempr(1:vic,R_UT1) = q(vis:vie,j,k,QU   ) - 0.5d0 *   dqy(vis:vie,j,k,3)

              ! Right w
              qtempr(1:vic,R_UT2) = q(vis:vie,j,k,QW) - 0.5d0 * dqy(vis:vie,j,k,4)

              ! Right rho - computed as sum(rhoY_k) below after species
              qtempr(1:vic,R_RHO) = 0.d0
              do nsp = 1,nspecies

                qtempr(1:vic,R_Y-1+nsp) = q(vis:vie,j,k,QFS-1+nsp) &
                                          * q(vis:vie,j,k,QRHO) - 0.5d0*(dqy(vis:vie,j,k,4+nsp) &
                                          + q(vis:vie,j,k,QFS-1+nsp) &
                                             * (dqy(vis:vie,j,k,1) + dqy(vis:vie,j,k,2)) &
                                             /cspeed(1:vic) )
                qtempr(1:vic,R_RHO) = qtempr(1:vic,R_RHO) + qtempr(1:vic,R_Y-1+nsp)

               enddo
              do nsp = 1,nspecies
                qtempr(1:vic,R_Y-1+nsp) = qtempr(1:vic,R_Y-1+nsp)/qtempr(1:vic,R_RHO)
              enddo

              ! Small and avg c
             cavg(1:vic) = HALF * ( qaux(vis:vie,j,k,QC) + qaux(vis:vie,j-1,k,QC) )
             csmall(1:vic) = min( qaux(vis:vie,j,k,QCSML), qaux(vis:vie,j-1,k,QCSML) )

            ! TODO: Make this loop a call to a vector EOS routine
             do vii = 1, vic
               ! Have p, rhoY (composition is rhoY), rho 
               !  - evaluate T, use that to evaluate internal energy

               eos_state%rho = qtempl(vii,R_RHO)
               eos_state%p = qtempl(vii,R_P)
               eos_state%massfrac = qtempl(vii,R_Y:R_Y-1+nspecies)
               !dir$ inline recursive
               call eos_rp(eos_state)
               rhoe_l(vii) = eos_state%rho * eos_state%e
               gamc_l(vii) = eos_state%gam1

               eos_state%rho = qtempr(vii,R_RHO)
               eos_state%p = qtempr(vii,R_P)
               eos_state%massfrac = qtempr(vii,R_Y:R_Y-1+nspecies)
               !dir$ inline recursive
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
        call riemann_md_vec( &
             qtempl(1:vic,R_RHO), qtempl(1:vic,R_UN), qtempl(1:vic,R_UT1), qtempl(1:vic,R_UT2), qtempl(1:vic,R_P), rhoe_l(1:vic), qtempl(1:vic,R_Y:R_Y-1+nspecies), gamc_l(1:vic),&
             qtempr(1:vic,R_RHO), qtempr(1:vic,R_UN), qtempr(1:vic,R_UT1), qtempr(1:vic,R_UT2), qtempr(1:vic,R_P), rhoe_r(1:vic), qtempr(1:vic,R_Y:R_Y-1+nspecies), gamc_r(1:vic),&
             v_gd(1:vic), u_gd(1:vic), w_gd(1:vic), p_gd(1:vic), game_gd(1:vic), re_gd(1:vic), r_gd(1:vic), ustar(1:vic),&
             eos_state, nspecies,&
             flux_tmp(1:vic,URHO), flux_tmp(1:vic,UMY), flux_tmp(1:vic,UMX), flux_tmp(1:vic,UMZ), flux_tmp(1:vic,UEDEN), flux_tmp(1:vic,UEINT), &
             bc_test_val, csmall(1:vic), cavg(1:vic), vic )

        vii = 0
        do vi = vis, vie ! source/dest array index space
          vii = vii + 1 ! work array index


                ! Compute species flux like passive scalar from intermediate state
                do nsp = 0, nspecies-1
                   if (ustar(vii) .gt. ZERO) then
                      flux_tmp(vii,UFS+nsp) = flux_tmp(vii,URHO)*qtempl(vii,R_Y+nsp)
                   else if (ustar(vii) .lt. ZERO) then
                      flux_tmp(vii,UFS+nsp) = flux_tmp(vii,URHO)*qtempr(vii,R_Y+nsp)
                   else
                      flux_tmp(vii,UFS+nsp) = flux_tmp(vii,URHO)&
                           *HALF*(qtempl(vii,R_Y+nsp) + qtempr(vii,R_Y+nsp) )
                   endif
                enddo

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
                 flux2(vis:vie,j,k, ivar) = flux2(vis:vie,j,k, ivar ) + flux_tmp(1:vic,ivar) * Ay(vis:vie,j,k)

             enddo
          enddo
       enddo
    enddo
    call bl_proffortfuncstop_int(3)

    call bl_proffortfuncstart_int(4)
          ! Compute all slopes at kc (k3d)
#ifdef PELEC_USE_EB
    call slopez(q,flatn,qd_lo,qd_hi, &
         dqz,qt_lo,qt_hi, &
         lo(1)-nextra,lo(2)-nextra,lo(3)-nextra,  &
         hi(1)+nextra,hi(2)+nextra,hi(3)+nextra,QVAR,NQAUX, &
         domlo,domhi, &
         qaux, qa_lo, qa_hi, &
         flag, fglo, fghi)
#else
    call slopez(q,flatn,qd_lo,qd_hi, &
         dqz,qt_lo,qt_hi, &
         qaux, qa_lo, qa_hi, &
         lo(1),lo(2),lo(3),hi(1),hi(2),hi(3),QVAR,NQAUX, &
         qaux, qa_lo, qa_hi)
#endif
    call bl_proffortfuncstop_int(4)

    call bl_proffortfuncstart_int(5)
    do k = lo(3)-nextra+1, hi(3)+nextra
       do j = lo(2)-nextra, hi(2)+nextra

          do i = lo(1)-nextra, hi(1)+nextra, VECLEN
             vis = i
             vie = min(vis+VECLEN-1, hi(1)+nextra)
             vic = vie - vis + 1


             cspeed(1:vic) = qaux(vis:vie,j,k-1,QC)

             ! Left w
             qtempl(1:vic,R_UN ) = q(vis:vie,j,k-1,QW   ) + 0.5d0 *((dqz(vis:vie,j,k-1,2)-dqz(vis:vie,j,k-1,1))/q(vis:vie,j,k-1,QRHO))

             ! Left p
             qtempl(1:vic,R_P  )=  q(vis:vie,j,k-1,QPRES) + 0.5d0 * (dqz(vis:vie,j,k-1,1)+dqz(vis:vie,j,k-1,2))*cspeed(1:vic)

             ! Left u
             qtempl(1:vic,R_UT1) = q(vis:vie,j,k-1,QU   ) + 0.5d0 *  dqz(vis:vie,j,k-1,3)

             ! Left v
             qtempl(1:vic,R_UT2) = q(vis:vie,j,k-1,QV   ) + 0.5d0 *  dqz(vis:vie,j,k-1,4)


             ! Left rho - computed as sum(rhoY_k) below after species
             qtempl(1:vic,R_RHO) = 0.d0
             do nsp = 1,nspecies

               qtempl(1:vic,R_Y-1+nsp) = q(vis:vie,j,k-1,QFS-1+nsp) &
                                         * q(vis:vie,j,k-1,QRHO) + 0.5d0*(dqz(vis:vie,j,k-1,4+nsp) &
                                         + q(vis:vie,j,k-1,QFS-1+nsp) &
                                            * (dqz(vis:vie,j,k-1,1) + dqz(vis:vie,j,k-1,2)) &
                                            /cspeed(1:vic) )
               qtempl(1:vic,R_RHO) = qtempl(1:vic,R_RHO) + qtempl(1:vic,R_Y-1+nsp)
              enddo

             do nsp = 1,nspecies
               qtempl(1:vic,R_Y-1+nsp) = qtempl(1:vic,R_Y-1+nsp)/qtempl(1:vic,R_RHO)
             enddo

             cspeed(1:vic) = qaux(vis:vie,j,k,QC)

             ! Right w
             qtempr(1:vic,R_UN ) = q(vis:vie,j,k,QW   ) - 0.5d0 * ((dqz(vis:vie,j,k,2)-dqz(vis:vie,j,k,1))/q(vis:vie,j,k,QRHO))

             ! Right p
             qtempr(1:vic,R_P  ) = q(vis:vie,j,k,QPRES) - 0.5d0 *  (dqz(vis:vie,j,k,1)+dqz(vis:vie,j,k,2))*cspeed(1:vic)

             ! Right u
             qtempr(1:vic,R_UT1) = q(vis:vie,j,k,QU   ) - 0.5d0 *   dqz(vis:vie,j,k,3)

             ! Right v
             qtempr(1:vic,R_UT2) = q(vis:vie,j,k,QV   ) - 0.5d0 *   dqz(vis:vie,j,k,4)

             ! Right rho - computed as sum(rhoY_k) below after species
             qtempr(1:vic,R_RHO) = 0.d0
             do nsp = 1,nspecies

                qtempr(1:vic,R_Y-1+nsp) = q(vis:vie,j,k,QFS-1+nsp) &
                                          * q(vis:vie,j,k,QRHO) - 0.5d0*(dqz(vis:vie,j,k,4+nsp) &
                                          + q(vis:vie,j,k,QFS-1+nsp) &
                                             * (dqz(vis:vie,j,k,1) + dqz(vis:vie,j,k,2)) &
                                             /cspeed(1:vic) )
                qtempr(1:vic,R_RHO) = qtempr(1:vic,R_RHO) + qtempr(1:vic,R_Y-1+nsp)


             enddo
             do nsp = 1,nspecies
               qtempr(1:vic,R_Y-1+nsp) = qtempr(1:vic,R_Y-1+nsp)/qtempr(1:vic,R_RHO)
             enddo

             ! Small and avg c
             cavg(1:vic) = HALF * ( qaux(vis:vie,j,k,QC) + qaux(vis:vie,j,k-1,QC) )
             csmall(1:vic) = min( qaux(vis:vie,j,k,QCSML), qaux(vis:vie,j,k-1,QCSML) )
             ! TODO: Make this loop a call to a vector EOS routine
             do vii = 1, vic
                ! Have p, rhoY (composition is rhoY), rho 
                !  - evaluate T, use that to evaluate internal energy
                eos_state%rho = qtempl(vii,R_RHO)
                eos_state%p = qtempl(vii,R_P)
                eos_state%massfrac = qtempl(vii,R_Y:R_Y-1+nspecies)
                !dir$ inline recursive
                call eos_rp(eos_state)
                rhoe_l(vii) = eos_state%rho * eos_state%e
                gamc_l(vii) = eos_state%gam1

                eos_state%rho = qtempr(vii,R_RHO)
                eos_state%p = qtempr(vii,R_P)
                eos_state%massfrac = qtempr(vii,R_Y:R_Y-1+nspecies)
                !dir$ inline recursive
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
             call riemann_md_vec( &
                  qtempl(1:vic,R_RHO), qtempl(1:vic,R_UN), qtempl(1:vic,R_UT1), qtempl(1:vic,R_UT2), qtempl(1:vic,R_P), rhoe_l(1:vic), qtempl(1:vic,R_Y:R_Y-1+nspecies), gamc_l(1:vic),&
                  qtempr(1:vic,R_RHO), qtempr(1:vic,R_UN), qtempr(1:vic,R_UT1), qtempr(1:vic,R_UT2), qtempr(1:vic,R_P), rhoe_r(1:vic), qtempr(1:vic,R_Y:R_Y-1+nspecies), gamc_r(1:vic),&
                  w_gd(1:vic), u_gd(1:vic), v_gd(1:vic), p_gd(1:vic), game_gd(1:vic), re_gd(1:vic), r_gd(1:vic), ustar(1:vic),&
                  eos_state, nspecies,&
                  flux_tmp(1:vic,URHO), flux_tmp(1:vic,UMZ), flux_tmp(1:vic,UMX), flux_tmp(1:vic,UMY), flux_tmp(1:vic,UEDEN), flux_tmp(1:vic,UEINT), &
                  bc_test_val, csmall(1:vic), cavg(1:vic), vic )

             vii = 0
             do vi = vis, vie ! source/dest array index space
                vii = vii + 1 ! work array index


                do nsp = 0, nspecies-1
                   if (ustar(vii) .gt. ZERO) then
                      flux_tmp(vii,UFS+nsp) = flux_tmp(vii,URHO)*qtempl(vii,R_Y+nsp)
                   else if (ustar(vii) .lt. ZERO) then
                      flux_tmp(vii,UFS+nsp) = flux_tmp(vii,URHO)*qtempr(vii,R_Y+nsp)
                   else
                      flux_tmp(vii,UFS+nsp) = flux_tmp(vii,URHO)&
                           *HALF*(qtempl(vii,R_Y+nsp) + qtempr(vii,R_Y+nsp) )
                   endif
                enddo

                
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
                flux3(vis:vie,j,k, ivar) = flux3(vis:vie,j,k,ivar ) + flux_tmp(1:vic,ivar) * Az(vis:vie,j,k)
             enddo
          enddo
       enddo
    enddo
    call bl_proffortfuncstop_int(5)

    ! Done computing flux through regular faces.

    ! TODO: Flux through wall face here...
#ifdef PELEC_USE_EB

    call bl_proffortfuncstart_int(6)
    full_area = h(1)**(dim - 1)

    ! Loop over cut cells only - need to pass in a list of these (lift the list Marc built for diffusion)
    do L = 0, nebflux-1, VECLEN
       vis = L
       vie = min(vis+VECLEN-1, nebflux-1)
       vic = vie - vis +1

       do vii = 1, vic
          i = ebg(vis+vii-1) % iv(0)
          j = ebg(vis+vii-1) % iv(1)
          k = ebg(vis+vii-1) % iv(2)
          if (i.ge.lo(1)-nextra .and. i.le.hi(1)+nextra &
               .and. j.ge.lo(2)-nextra .and. j.le.hi(2)+nextra &
               .and. k.ge.lo(3)-nextra .and. k.le.hi(3)+nextra ) then
           ! Get normal - checkto make sure this is facing outward!
             eb_norm = ebg(vis+vii-1)%eb_normal

    !        eb_norm(1) = ebg(vis+vii-1)%eb_apertureX(2) - ebg(vis+vii-1)%eb_apertureX(1)
    !        eb_norm(2) = ebg(vis+vii-1)%eb_apertureY(2) - ebg(vis+vii-1)%eb_apertureY(1)
    !        eb_norm(3) = ebg(vis+vii-1)%eb_apertureZ(2) - ebg(vis+vii-1)%eb_apertureZ(1)

             eb_norm = eb_norm / sqrt(eb_norm(1)**2 + eb_norm(2)**2 + eb_norm(3)**2)


             ! Replace q for cut cell with average of neighborhood so that ebflux is consistent with cell merging
             call get_neighbor_cells(flag(i,j,k),nbr)

             if( vfrac(i,j,k) < eb_small_vfrac ) then
                nbr(0,0,0) = 0
                sum_kappa = sum(nbr(-1:1,-1:1,-1:1) * vfrac(i-1:i+1,j-1:j+1,k-1:k+1))

                ! Construct left state from volume weighted average of neighbourhood
                sum_nbrs =   sum(nbr(-1:1,-1:1,-1:1) * vfrac(i-1:i+1,j-1:j+1,k-1:k+1) * qaux(i-1:i+1,j-1:j+1,k-1:k+1,QC))
                cspeed(vii) = sum_nbrs/sum_kappa

                qtempl(vii,R_UN) = 0.0d0
                sum_nbrs =   sum(nbr(-1:1,-1:1,-1:1) * vfrac(i-1:i+1,j-1:j+1,k-1:k+1) * q(i-1:i+1,j-1:j+1,k-1:k+1,QU))
                qtempl(vii,R_UN) = qtempl(vii,R_UN) -sum_nbrs/sum_kappa*eb_norm(1)

                sum_nbrs =   sum(nbr(-1:1,-1:1,-1:1) * vfrac(i-1:i+1,j-1:j+1,k-1:k+1) * q(i-1:i+1,j-1:j+1,k-1:k+1,QV))
                qtempl(vii,R_UN) = qtempl(vii,R_UN) -sum_nbrs/sum_kappa*eb_norm(2)

                sum_nbrs =   sum(nbr(-1:1,-1:1,-1:1) * vfrac(i-1:i+1,j-1:j+1,k-1:k+1) * q(i-1:i+1,j-1:j+1,k-1:k+1,QW))
                qtempl(vii,R_UN) = qtempl(vii,R_UN) -sum_nbrs/sum_kappa*eb_norm(3)

                qtempl(vii,R_UT1) = 0.0
                qtempl(vii,R_UT2) = 0.0

                sum_nbrs =   sum(nbr(-1:1,-1:1,-1:1) * vfrac(i-1:i+1,j-1:j+1,k-1:k+1) * q(i-1:i+1,j-1:j+1,k-1:k+1,QPRES))
                qtempl(vii,R_P  ) = sum_nbrs/sum_kappa

                sum_nbrs =   sum(nbr(-1:1,-1:1,-1:1) * vfrac(i-1:i+1,j-1:j+1,k-1:k+1) * q(i-1:i+1,j-1:j+1,k-1:k+1,QRHO))

                qtempl(vii,R_RHO) = sum_nbrs/sum_kappa

                do nsp = 1,nspecies
                   sum_nbrs =   sum(nbr(-1:1,-1:1,-1:1) * vfrac(i-1:i+1,j-1:j+1,k-1:k+1) * q(i-1:i+1,j-1:j+1,k-1:k+1,QFS-1+nsp))
                   qtempl(vii,R_Y-1+nsp) = sum_nbrs/sum_kappa
                enddo

                ! Flip the velocity about the normal for the right state - will use left
                ! state for remainder of right state
                qtempr(vii,R_UN  ) = -1.0*qtempl(vii,R_UN)

                ! Small and avg c
                sum_nbrs =   sum(nbr(-1:1,-1:1,-1:1) * vfrac(i-1:i+1,j-1:j+1,k-1:k+1) * qaux(i-1:i+1,j-1:j+1,k-1:k+1,QC))
                cavg(vii)  = sum_nbrs/sum_kappa

                sum_nbrs =   sum(nbr(-1:1,-1:1,-1:1) * vfrac(i-1:i+1,j-1:j+1,k-1:k+1) * qaux(i-1:i+1,j-1:j+1,k-1:k+1,QCSML))
                csmall(vii) = sum_nbrs/sum_kappa


             else
                ! Assume left state is the cell centered state - normal veclocity
                cspeed(vii) = qaux(i,j,k,QC)

                qtempl(vii,R_UN ) = - q(i,j,k,QU)*eb_norm(1) &
                     -                q(i,j,k,QV)*eb_norm(2) &
                     -                q(i,j,k,QW)*eb_norm(3)

                qtempl(vii,R_UT1) = 0.0
                qtempl(vii,R_UT2) = 0.0
                qtempl(vii,R_P  ) = q(i,j,k,QPRES)
                qtempl(vii,R_RHO) = q(i,j,k,QRHO)

                do nsp = 1,nspecies
                   qtempl(vii,R_Y-1+nsp) = q(i,j,k,QFS-1+nsp)
                enddo

                ! Flip the velocity about the normal for the right state - will use left
                ! state for remainder of right state
                qtempr(vii,R_UN  ) = -1.0*qtempl(vii,R_UN)

                ! Small and avg c
                cavg(vii)  =  qaux(i,j,k,QC)
                csmall(vii) = qaux(i,j,k,QCSML)

             endif




          endif
       enddo


    !     TODO: Make this loop a call to a vector EOS routine
       do vii = 1, vic

          ! Have p, rhoY (composition is rhoY), rho 
          i = ebg(vis+vii-1) % iv(0)
          j = ebg(vis+vii-1) % iv(1)
          k = ebg(vis+vii-1) % iv(2)
          if (i.ge.lo(1)-nextra .and. i.le.hi(1)+nextra &
               .and. j.ge.lo(2)-nextra .and. j.le.hi(2)+nextra &
               .and. k.ge.lo(3)-nextra .and. k.le.hi(3)+nextra ) then
          !  - evaluate T, use that to evaluate internal energy
        
                eos_state%rho =      qtempl(vii,R_RHO)
                eos_state%p   =      qtempl(vii,R_P)
                eos_state%massfrac = qtempl(vii,R_Y:R_Y-1+nspecies)
                !dir$ inline recursive
                call eos_rp(eos_state)
                rhoe_l(vii) = eos_state%rho * eos_state%e
                gamc_l(vii) = eos_state%gam1
         endif

       enddo

       ! Solve Riemann problem; store flux in flux4 - cell centered data structure
       vii = 0
       do vi = vis, vie ! source/dest array index space
          vii = vii + 1 ! work array index

          i = ebg(vis+vii-1) % iv(0)
          j = ebg(vis+vii-1) % iv(1)
          k = ebg(vis+vii-1) % iv(2)
          ! this is going to prevent vectorization, but until vector riemann solver
          ! is in place might as well avoid copying into a temporary 
          if (       i.ge.lo(1)-nextra+1 .and. i.le.hi(1)+nextra-1 &
               .and. j.ge.lo(2)-nextra+1 .and. j.le.hi(2)+nextra-1 &
               .and. k.ge.lo(3)-nextra+1 .and. k.le.hi(3)+nextra-1 ) then

             call riemann_md_singlepoint( &
                  qtempl(vii,R_RHO), qtempl(vii,R_UN), qtempl(vii,R_UT1), qtempl(vii,R_UT2), qtempl(vii,R_P), rhoe_l(vii), qtempl(vii,R_Y:R_Y-1+nspecies), gamc_l(vii),&
                  qtempl(vii,R_RHO), qtempr(vii,R_UN), qtempl(vii,R_UT1), qtempl(vii,R_UT2), qtempl(vii,R_P), rhoe_l(vii), qtempl(vii,R_Y:R_Y-1+nspecies), gamc_l(vii),&
                  u_gd(vii), v_gd(vii), w_gd(vii), p_gd(vii), game_gd(vii), re_gd(vii), r_gd(vii), ustar(vii),&
                  eos_state, gdnv_state, nspecies,&
                  flux_tmp(vii,URHO), flux_tmp(vii,UMX), flux_tmp(vii,UMY), flux_tmp(vii,UMZ), flux_tmp(vii,UEDEN), flux_tmp(vii,UEINT), &
                  idir, coord_type, bc_test_val, csmall(vii), cavg(vii) )
          
  
             eb_norm = ebg(vis+vii-1)%eb_normal

             eb_norm = eb_norm / sqrt(eb_norm(1)**2 + eb_norm(2)**2 + eb_norm(3)**2)

             flux_tmp(vii,UMY) = -flux_tmp(vii,UMX) * eb_norm(2)
             flux_tmp(vii,UMZ) = -flux_tmp(vii,UMX) * eb_norm(3)
             flux_tmp(vii,UMX) = -flux_tmp(vii,UMX) * eb_norm(1)
             !   Compute species flux like passive scalar from intermediate state
             do nsp = 0, nspecies-1
                flux_tmp(vii,UFS+nsp) = flux_tmp(vii,URHO)*qtempl(vii,R_Y+nsp)
             enddo
          endif

       enddo ! End future vector loop

       ! Copy result into ebflux vector. Being a bit chicken here and only copy values where ebg % iv is within box
       vii = 0
       do vi = vis, vie ! source/dest array index space
          vii = vii + 1 ! work array index

          i = ebg(vis+vii-1) % iv(0)
          j = ebg(vis+vii-1) % iv(1)
          k = ebg(vis+vii-1) % iv(2)
          ! this is going to prevent vectorization, but until vector riemann solver
          ! is in place might as well avoid copying into a temporary 
          if (       i.ge.lo(1)-nextra+1 .and. i.le.hi(1)+nextra-1 &
               .and. j.ge.lo(2)-nextra+1 .and. j.le.hi(2)+nextra-1 &
               .and. k.ge.lo(3)-nextra+1 .and. k.le.hi(3)+nextra-1 ) then

             do ivar = 1, NVAR
                ebflux(vi,ivar) = ebflux(vi,ivar) + flux_tmp(vii,ivar) * ebg(vi)%eb_area * full_area
             enddo

          endif
       enddo ! End future vector loop
    enddo ! End loop over cut cells

    call bl_proffortfuncstop_int(6)
#endif

    ! Deallocate arrays
    call bl_deallocate (dqx)
    call bl_deallocate (dqy)
    call bl_deallocate (dqz)

    call destroy(eos_state)
    call destroy(gdnv_state)

    ! Now, all faces flux done - now ready for Marc's magic; flux in x-direction 
    ! is loaded into flux1
    ! where flux1(i,j,k,N) has edge based flux of rho,u,v,w,eden,eint,species indexed by
    ! URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS:UFS+nspecies-1, and similar for flux2, flux3
    ! EB face flux assuming wall BC is loaded into ebflux with same ordering

    ! Call Marc's routine to multiply fluxes by aperatures and interpolate to face centers 
    ! (includes multiplying covered faces by 0.0) - or just return the fluxes and add 'em
    ! up with the diffusion fluxes before calling Marc's routine to compute hybrid divergence

    call bl_proffortfuncstart_int(7)

    do ivar=1,NVAR
       do k = lo(3)-nextra+1, hi(3)+nextra-1
          do j = lo(2)-nextra+1, hi(2)+nextra-1
             do i = lo(1)-nextra+1, hi(1)+nextra-1

                D(i,j,k,ivar) = - (flux1(i+1,j,k,ivar) - flux1(i,j,k,ivar) &
                     +             flux2(i,j+1,k,ivar) - flux2(i,j,k,ivar) &
                     +             flux3(i,j,k+1,ivar) - flux3(i,j,k,ivar) )/V(i,j,k)

             enddo
          enddo
       enddo
    enddo

    call bl_proffortfuncstop_int(7)

  end subroutine pc_hyp_mol_flux
end module hyp_advection_module 
