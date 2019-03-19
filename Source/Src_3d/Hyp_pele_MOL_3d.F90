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

    use meth_params_module, only : plm_iorder, QVAR, NVAR, QPRES, QRHO, QU, QV, QW, &
                                   QFS, QC, QCSML, NQAUX, nadv, &
                                   URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, UTEMP, UFX, UFA, &
                                   small_dens, small_pres
    use slope_module, only : slopex, slopey, slopez
    use actual_network, only : naux
    use eos_module, only : eos_rp1
    use chemistry_module, only: Ru

    implicit none

    integer, parameter  :: nspec_2=9
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
    double precision, intent(inout) ::   ebflux(0:nebflux-1,1:NVAR)
    integer,            intent(in   ) :: Nebg
    type(eb_bndry_geom),intent(in   ) :: ebg(0:Nebg-1)    
    double precision :: eb_norm(3), full_area
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

    integer :: i, j, k, n, nsp, L, ivar
    integer :: qt_lo(3), qt_hi(3)
    integer :: ilo1, ilo2, ihi1, ihi2, ilo3, ihi3

    ! Left and right state arrays (edge centered, cell centered)
    double precision, pointer :: dqx(:,:,:,:), dqy(:,:,:,:), dqz(:,:,:,:)

    ! Other left and right state arrays
    double precision :: qtempl(1:5+nspec_2)
    double precision :: qtempr(1:5+nspec_2)
    double precision :: rhoe_l
    double precision :: rhoe_r
    double precision :: cspeed
    double precision :: gamc_l
    double precision :: gamc_r
    double precision :: cavg
    double precision :: csmall

    ! Scratch for neighborhood of cut cells
    integer :: nbr(-1:1,-1:1,-1:1)

    ! Riemann solve work arrays
    double precision:: u_gd, v_gd, w_gd, &
         p_gd, game_gd, re_gd, &
         r_gd, ustar
    double precision :: flux_tmp(NVAR)
    integer, parameter :: idir = 1
    integer :: nextra
    integer, parameter :: coord_type = 0
    integer, parameter :: bc_test_val = 1
    
    double precision :: eos_state_rho
    double precision :: eos_state_p
    double precision :: eos_state_massfrac(nspec_2)
    double precision :: eos_state_gam1
    double precision :: eos_state_e
    double precision :: eos_state_cs

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

    flux_tmp = 0.d0
    ilo1=lo(1)-nextra
    ilo2=lo(2)-nextra
    ilo3=lo(3)-nextra
    ihi1=hi(1)+nextra
    ihi2=hi(2)+nextra
    ihi3=hi(3)+nextra

    do L=1,3
       qt_lo(L) = lo(L) - nextra
       qt_hi(L) = hi(L) + nextra
       !if (qt_lo(L)-1 .lt. qd_lo(L) .or. qt_hi(L)+1 .gt. qd_hi(L)) then
       !   stop 1
       !endif
    enddo

    allocate(dqx(qt_lo(1):qt_hi(1), qt_lo(2):qt_hi(2), qt_lo(3):qt_hi(3), 1:QVAR))
    allocate(dqy(qt_lo(1):qt_hi(1), qt_lo(2):qt_hi(2), qt_lo(3):qt_hi(3), 1:QVAR))
    allocate(dqz(qt_lo(1):qt_hi(1), qt_lo(2):qt_hi(2), qt_lo(3):qt_hi(3), 1:QVAR))

    !$acc update device(ru,naux,plm_iorder,qvar,nvar,qpres,qrho,qu,qv,qw,qfs,qc,qcsml,nqaux,nadv,urho,umx,umy,umz,ueden,ueint,ufs,utemp,ufx,ufa)
    !$acc enter data create(dqx,dqy,dqz) copyin(flux1,flux2,flux3,d) copyin(v,ax,ay,az,q,flatn,qd_lo,qd_hi,qt_lo,qt_hi,domlo,domhi,qaux,qa_lo,qa_hi,flag,fglo,fghi,lo,hi)

    !$acc parallel loop gang vector collapse(4) present(qt_lo,qt_hi,dqx)
    do n = 1, qvar
       do k = qt_lo(3), qt_hi(3)
          do j = qt_lo(2), qt_hi(2)
             do i = qt_lo(1), qt_hi(1)
                dqx(i,j,k,n) = 0.d0
             enddo
          enddo
       enddo
    enddo
    !$acc end parallel
    if(plm_iorder.ne.1) then
       call slopex(q,flatn,qd_lo,qd_hi, &
                  dqx,qt_lo,qt_hi, &
                  ilo1,ilo2,ilo3, &
                  ihi1,ihi2,ihi3,qvar,nqaux, &
                  domlo,domhi, &
                  qaux,qa_lo,qa_hi, &
                  flag,fglo,fghi)
    end if

    !$acc kernels present(q,qaux,dqx,ax,flux1)
    !$acc loop gang vector collapse(3) private(n,qtempl,qtempr,gamc_l,rhoe_l,rhoe_r,gamc_r,u_gd, v_gd, w_gd, p_gd, game_gd, re_gd, r_gd, ustar, eos_state_rho, eos_state_p, eos_state_massfrac, eos_state_e, eos_state_gam1, eos_state_cs, flux_tmp, csmall, cavg, vic, ivar) 
    do k = ilo3, ihi3
       do j = ilo2, ihi2
          do i = ilo1+1, ihi1
             qtempl(R_UN) = q(i-1,j,k,QU) + 0.5d0 * ((dqx(i-1,j,k,2) - dqx(i-1,j,k,1)) / q(i-1,j,k,QRHO))
             qtempl(R_P) = q(i-1,j,k,QPRES) + 0.5d0 * (dqx(i-1,j,k,1) + dqx(i-1,j,k,2)) * qaux(i-1,j,k,QC)
             qtempl(R_UT1) = q(i-1,j,k,QV) + 0.5d0 * dqx(i-1,j,k,3)
             qtempl(R_UT2) = q(i-1,j,k,QW) + 0.5d0 * dqx(i-1,j,k,4)
             qtempl(R_RHO) = 0.d0
             do n = 1,nspec_2
                qtempl(R_Y-1+n) = q(i-1,j,k,QFS-1+n) * q(i-1,j,k,QRHO) + 0.5d0 * (dqx(i-1,j,k,4+n) &
                                  + q(i-1,j,k,QFS-1+n) * (dqx(i-1,j,k,1) + dqx(i-1,j,k,2)) / qaux(i-1,j,k,QC))
                qtempl(R_RHO) = qtempl(R_RHO) + qtempl(R_Y-1+n)
             enddo

             do n = 1,nspec_2
               qtempl(R_Y-1+n) = qtempl(R_Y-1+n) / qtempl(R_RHO)
             enddo

             qtempr(R_UN) = q(i,j,k,QU) - 0.5d0 * ((dqx(i,j,k,2) - dqx(i,j,k,1)) / q(i,j,k,QRHO))
             qtempr(R_P) = q(i,j,k,QPRES) - 0.5d0 * (dqx(i,j,k,1) + dqx(i,j,k,2)) * qaux(i,j,k,QC)
             qtempr(R_UT1) = q(i,j,k,QV) - 0.5d0 * dqx(i,j,k,3)
             qtempr(R_UT2) = q(i,j,k,QW) - 0.5d0 * dqx(i,j,k,4)
             qtempr(R_RHO) = 0.d0

             do n = 1,nspec_2
                qtempr(R_Y-1+n) = q(i,j,k,QFS-1+n) * q(i,j,k,QRHO) - 0.5d0 * &
                                    (dqx(i,j,k,4+n) + q(i,j,k,QFS-1+n) * &
                                    (dqx(i,j,k,1) + dqx(i,j,k,2)) / qaux(i,j,k,QC))
                qtempr(R_RHO) = qtempr(R_RHO) + qtempr(R_Y-1+n)
             enddo

             do n = 1,nspec_2
                qtempr(R_Y-1+n) = qtempr(R_Y-1+n)/qtempr(R_RHO)
             enddo

             cavg = 0.5d0 * (qaux(i,j,k,QC) + qaux(i-1,j,k,QC))
             csmall = min(qaux(i,j,k,QCSML), qaux(i-1,j,k,QCSML))

             eos_state_rho = qtempl(R_RHO)
             eos_state_p = qtempl(R_P)
             eos_state_massfrac = qtempl(R_Y:R_Y-1+nspec_2)
             call eos_rp1(eos_state_rho, eos_state_p, eos_state_massfrac, eos_state_e, eos_state_gam1, eos_state_cs, nspec_2)
             rhoe_l = eos_state_rho * eos_state_e
             gamc_l = eos_state_gam1

             eos_state_rho = qtempr(R_RHO)
             eos_state_p = qtempr(R_P)
             eos_state_massfrac = qtempr(R_Y:R_Y-1+nspec_2)
             call eos_rp1(eos_state_rho, eos_state_p, eos_state_massfrac, eos_state_e, eos_state_gam1, eos_state_cs, nspec_2)
             rhoe_r = eos_state_rho * eos_state_e
             gamc_r = eos_state_gam1

             call riemann_md_vec(qtempl(R_RHO), qtempl(R_UN), qtempl(R_UT1), qtempl(R_UT2), &
                                 qtempl(R_P), rhoe_l, qtempl(R_Y:R_Y-1+nspec_2), gamc_l, &
                                 qtempr(R_RHO), qtempr(R_UN), qtempr(R_UT1), qtempr(R_UT2), &
                                 qtempr(R_P), rhoe_r, qtempr(R_Y:R_Y-1+nspec_2), gamc_r,&
                                 u_gd, v_gd, w_gd, p_gd, game_gd, re_gd, r_gd, ustar, &
                                 eos_state_rho, eos_state_p, eos_state_massfrac, &
                                 eos_state_e, eos_state_gam1, eos_state_cs, nspec_2, &
                                 flux_tmp(URHO), flux_tmp(UMX), flux_tmp(UMY), flux_tmp(UMZ), &
                                 flux_tmp(UEDEN), flux_tmp(UEINT), bc_test_val, csmall, cavg, vic)

             do n = 0, nspec_2-1
                flux_tmp(UFS+n) = merge(flux_tmp(URHO)*qtempl(R_Y+n), flux_tmp(URHO)*qtempr(R_Y+n), ustar .ge. 0.d0)
                flux_tmp(UFS+n) = merge(flux_tmp(URHO)*0.5d0*(qtempl(R_Y+n) + qtempr(R_Y+n)), flux_tmp(UFS+n), ustar .eq. 0.d0)
             enddo

             flux_tmp(UTEMP) = 0.0
             do n = UFX, UFX+naux
                flux_tmp(n) = merge(0.d0, flux_tmp(n), naux .gt. 0)
             enddo
             do n = UFA, UFA+nadv
                flux_tmp(n) = merge(0.d0, flux_tmp(n), nadv .gt. 0)
             enddo

             do ivar = 1, NVAR
                flux1(i,j,k,ivar) = flux1(i,j,k,ivar) + flux_tmp(ivar) * ax(i,j,k)
             enddo
          enddo
       enddo
    enddo
    !$acc end kernels

    !$acc parallel loop gang vector collapse(4) present(qt_lo,qt_hi,dqy)
    do n = 1, qvar
       do k = qt_lo(3), qt_hi(3)
          do j = qt_lo(2), qt_hi(2)
             do i = qt_lo(1), qt_hi(1)
                dqy(i,j,k,n) = 0.d0
             enddo
          enddo
       enddo
    enddo
    !$acc end parallel
    if(plm_iorder.ne.1) then
       call slopey(q,flatn,qd_lo,qd_hi, &
                  dqy,qt_lo,qt_hi, &
                  ilo1,ilo2,ilo3, &
                  ihi1,ihi2,ihi3,qvar,nqaux, &
                  domlo,domhi, &
                  qaux,qa_lo,qa_hi, &
                  flag,fglo,fghi)
    end if

    !$acc kernels present(dqy,qaux,q,ay,flux2)
    !$acc loop gang vector collapse(3) private(n,qtempl,qtempr,gamc_l,rhoe_l,rhoe_r,gamc_r,u_gd, v_gd, w_gd, p_gd, game_gd, re_gd, r_gd, ustar, eos_state_rho, eos_state_p, eos_state_massfrac, eos_state_e, eos_state_gam1, eos_state_cs, flux_tmp, csmall, cavg, vic, ivar) 
    do k = ilo3, ihi3
       do j = ilo2+1, ihi2
          do i = ilo1, ihi1
             qtempl(R_UN) = q(i,j-1,k,QV) + 0.5d0 * ((dqy(i,j-1,k,2) - dqy(i,j-1,k,1)) / q(i,j-1,k,QRHO))
             qtempl(R_P) = q(i,j-1,k,QPRES) + 0.5d0 * (dqy(i,j-1,k,1) + dqy(i,j-1,k,2)) * qaux(i,j-1,k,QC)
             qtempl(R_UT1) = q(i,j-1,k,QU) + 0.5d0 * dqy(i,j-1,k,3)
             qtempl(R_UT2) = q(i,j-1,k,QW) + 0.5d0 * dqy(i,j-1,k,4)
             qtempl(R_RHO) = 0.d0
             do n = 1,nspec_2
               qtempl(R_Y-1+n) = q(i,j-1,k,QFS-1+n) * q(i,j-1,k,QRHO) + 0.5d0 * (dqy(i,j-1,k,4+n) &
                                 + q(i,j-1,k,QFS-1+n) * (dqy(i,j-1,k,1) + dqy(i,j-1,k,2)) / qaux(i,j-1,k,QC))
               qtempl(R_RHO) = qtempl(R_RHO) + qtempl(R_Y-1+n)
             enddo

             do n = 1,nspec_2
               qtempl(R_Y-1+n) = qtempl(R_Y-1+n)/qtempl(R_RHO)
             enddo

             qtempr(R_UN) = q(i,j,k,QV) - 0.5d0 * ((dqy(i,j,k,2) - dqy(i,j,k,1)) / q(i,j,k,QRHO))
             qtempr(R_P) = q(i,j,k,QPRES) - 0.5d0 * (dqy(i,j,k,1) + dqy(i,j,k,2)) * qaux(i,j,k,QC)
             qtempr(R_UT1) = q(i,j,k,QU) - 0.5d0 * dqy(i,j,k,3)
             qtempr(R_UT2) = q(i,j,k,QW) - 0.5d0 * dqy(i,j,k,4)
             qtempr(R_RHO) = 0.d0

             do n = 1,nspec_2
               qtempr(R_Y-1+n) = q(i,j,k,QFS-1+n) &
                                   * q(i,j,k,QRHO) - 0.5d0*(dqy(i,j,k,4+n) &
                                   + q(i,j,k,QFS-1+n) &
                                   * (dqy(i,j,k,1) + dqy(i,j,k,2)) &
                                   / qaux(i,j,k,QC))
               qtempr(R_RHO) = qtempr(R_RHO) + qtempr(R_Y-1+n)
             enddo

             do n = 1,nspec_2
               qtempr(R_Y-1+n) = qtempr(R_Y-1+n)/qtempr(R_RHO)
             enddo

             cavg = 0.5d0 * (qaux(i,j,k,QC) + qaux(i,j-1,k,QC))
             csmall = min(qaux(i,j,k,QCSML), qaux(i,j-1,k,QCSML))

             eos_state_rho = qtempl(R_RHO)
             eos_state_p = qtempl(R_P)
             eos_state_massfrac = qtempl(R_Y:R_Y-1+nspec_2)
             call eos_rp1(eos_state_rho, eos_state_p, eos_state_massfrac, eos_state_e, eos_state_gam1, eos_state_cs, nspec_2)
             rhoe_l = eos_state_rho * eos_state_e
             gamc_l = eos_state_gam1

             eos_state_rho = qtempr(R_RHO)
             eos_state_p = qtempr(R_P)
             eos_state_massfrac = qtempr(R_Y:R_Y-1+nspec_2)
             call eos_rp1(eos_state_rho, eos_state_p, eos_state_massfrac, eos_state_e, eos_state_gam1, eos_state_cs, nspec_2)
             rhoe_r = eos_state_rho * eos_state_e
             gamc_r = eos_state_gam1

             call riemann_md_vec( &
                  qtempl(R_RHO), qtempl(R_UN), qtempl(R_UT1), qtempl(R_UT2), qtempl(R_P), rhoe_l, qtempl(R_Y:R_Y-1+nspec_2), gamc_l,&
                  qtempr(R_RHO), qtempr(R_UN), qtempr(R_UT1), qtempr(R_UT2), qtempr(R_P), rhoe_r, qtempr(R_Y:R_Y-1+nspec_2), gamc_r,&
                  v_gd, u_gd, w_gd, p_gd, game_gd, re_gd, r_gd, ustar,&
                  eos_state_rho, eos_state_p, eos_state_massfrac, &
                  eos_state_e, eos_state_gam1, eos_state_cs, nspec_2,&
                  flux_tmp(URHO), flux_tmp(UMY), flux_tmp(UMX), flux_tmp(UMZ), flux_tmp(UEDEN), flux_tmp(UEINT), &
                  bc_test_val, csmall, cavg, vic)

             do n = 0, nspec_2-1
                flux_tmp(UFS+n) = merge(flux_tmp(URHO)*qtempl(R_Y+n), flux_tmp(URHO)*qtempr(R_Y+n), ustar .ge. 0.d0)
                flux_tmp(UFS+n) = merge(flux_tmp(URHO)*0.5d0*(qtempl(R_Y+n) + qtempr(R_Y+n)), flux_tmp(UFS+n), ustar .eq. 0.d0)
             enddo

             flux_tmp(UTEMP) = 0.0
             do n = UFX, UFX+naux
                flux_tmp(n) = merge(0.d0, flux_tmp(n), naux .gt. 0)
             enddo
             do n = UFA, UFA+nadv
                flux_tmp(n) = merge(0.d0, flux_tmp(n), nadv .gt. 0)
             enddo

             do ivar = 1, NVAR
                flux2(i,j,k,ivar) = flux2(i,j,k,ivar) + flux_tmp(ivar) * ay(i,j,k)
             enddo
          enddo
       enddo
    enddo
    !$acc end kernels

    !$acc parallel loop gang vector collapse(4) present(qt_lo,qt_hi,dqz)
    do n = 1, qvar
       do k = qt_lo(3), qt_hi(3)
          do j = qt_lo(2), qt_hi(2)
             do i = qt_lo(1), qt_hi(1)
                dqz(i,j,k,n) = 0.d0
             enddo
          enddo
       enddo
    enddo
    !$acc end parallel
    if(plm_iorder.ne.1) then
       call slopez(q,flatn,qd_lo,qd_hi, &
                  dqz,qt_lo,qt_hi, &
                  ilo1,ilo2,ilo3, &
                  ihi1,ihi2,ihi3,qvar,nqaux, &
                  domlo,domhi, &
                  qaux,qa_lo,qa_hi, &
                  flag,fglo,fghi)
    end if

    !$acc kernels present(q,qaux,dqz,az,flux3)
    !$acc loop gang vector collapse(3) private(n,qtempl,qtempr,gamc_l,rhoe_l,rhoe_r,gamc_r,u_gd, v_gd, w_gd, p_gd, game_gd, re_gd, r_gd, ustar, eos_state_rho, eos_state_p, eos_state_massfrac, eos_state_e, eos_state_gam1, eos_state_cs, flux_tmp, csmall, cavg, vic, ivar)
    do k = ilo3+1, ihi3
       do j = ilo2, ihi2
          do i = ilo1, ihi1
             qtempl(R_UN) = q(i,j,k-1,QW) + 0.5d0 * ((dqz(i,j,k-1,2) - dqz(i,j,k-1,1)) / q(i,j,k-1,QRHO))
             qtempl(R_P) = q(i,j,k-1,QPRES) + 0.5d0 * (dqz(i,j,k-1,1) + dqz(i,j,k-1,2)) * qaux(i,j,k-1,QC)
             qtempl(R_UT1) = q(i,j,k-1,QU) + 0.5d0 * dqz(i,j,k-1,3)
             qtempl(R_UT2) = q(i,j,k-1,QV) + 0.5d0 * dqz(i,j,k-1,4)
             qtempl(R_RHO) = 0.d0
             do n = 1,nspec_2
               qtempl(R_Y-1+n) = q(i,j,k-1,QFS-1+n) &
                                 * q(i,j,k-1,QRHO) + 0.5d0*(dqz(i,j,k-1,4+n) &
                                 + q(i,j,k-1,QFS-1+n) &
                                 * (dqz(i,j,k-1,1) + dqz(i,j,k-1,2)) &
                                 / qaux(i,j,k-1,QC))
               qtempl(R_RHO) = qtempl(R_RHO) + qtempl(R_Y-1+n)
             enddo

             do n = 1,nspec_2
               qtempl(R_Y-1+n) = qtempl(R_Y-1+n)/qtempl(R_RHO)
             enddo

             qtempr(R_UN) = q(i,j,k,QW) - 0.5d0 * ((dqz(i,j,k,2) - dqz(i,j,k,1)) / q(i,j,k,QRHO))
             qtempr(R_P) = q(i,j,k,QPRES) - 0.5d0 * (dqz(i,j,k,1) + dqz(i,j,k,2)) * qaux(i,j,k,QC)
             qtempr(R_UT1) = q(i,j,k,QU) - 0.5d0 * dqz(i,j,k,3)
             qtempr(R_UT2) = q(i,j,k,QV) - 0.5d0 * dqz(i,j,k,4)
             qtempr(R_RHO) = 0.d0

             do n = 1,nspec_2
                qtempr(R_Y-1+n) = q(i,j,k,QFS-1+n) &
                                  * q(i,j,k,QRHO) - 0.5d0*(dqz(i,j,k,4+n) &
                                  + q(i,j,k,QFS-1+n) &
                                  * (dqz(i,j,k,1) + dqz(i,j,k,2)) &
                                  / qaux(i,j,k,QC))
                qtempr(R_RHO) = qtempr(R_RHO) + qtempr(R_Y-1+n)
             enddo

             do n = 1,nspec_2
               qtempr(R_Y-1+n) = qtempr(R_Y-1+n)/qtempr(R_RHO)
             enddo

             cavg = 0.5d0 * (qaux(i,j,k,QC) + qaux(i,j,k-1,QC))
             csmall = min(qaux(i,j,k,QCSML), qaux(i,j,k-1,QCSML))

             eos_state_rho = qtempl(R_RHO)
             eos_state_p = qtempl(R_P)
             eos_state_massfrac = qtempl(R_Y:R_Y-1+nspec_2)
             call eos_rp1(eos_state_rho, eos_state_p, eos_state_massfrac, eos_state_e, eos_state_gam1, eos_state_cs, nspec_2)
             rhoe_l = eos_state_rho * eos_state_e
             gamc_l = eos_state_gam1

             eos_state_rho = qtempr(R_RHO)
             eos_state_p = qtempr(R_P)
             eos_state_massfrac = qtempr(R_Y:R_Y-1+nspec_2)
             call eos_rp1(eos_state_rho, eos_state_p, eos_state_massfrac, eos_state_e, eos_state_gam1, eos_state_cs, nspec_2)
             rhoe_r = eos_state_rho * eos_state_e
             gamc_r = eos_state_gam1

             call riemann_md_vec( &
                  qtempl(R_RHO), qtempl(R_UN), qtempl(R_UT1), qtempl(R_UT2), qtempl(R_P), rhoe_l, qtempl(R_Y:R_Y-1+nspec_2), gamc_l,&
                  qtempr(R_RHO), qtempr(R_UN), qtempr(R_UT1), qtempr(R_UT2), qtempr(R_P), rhoe_r, qtempr(R_Y:R_Y-1+nspec_2), gamc_r,&
                  w_gd, u_gd, v_gd, p_gd, game_gd, re_gd, r_gd, ustar,&
                  eos_state_rho, eos_state_p, eos_state_massfrac, &
                  eos_state_e, eos_state_gam1, eos_state_cs, nspec_2,&
                  flux_tmp(URHO), flux_tmp(UMZ), flux_tmp(UMX), flux_tmp(UMY), flux_tmp(UEDEN), flux_tmp(UEINT), &
                  bc_test_val, csmall, cavg, vic)

             do n = 0, nspec_2-1
                flux_tmp(UFS+n) = merge(flux_tmp(URHO)*qtempl(R_Y+n), flux_tmp(URHO)*qtempr(R_Y+n), ustar .ge. 0.d0)
                flux_tmp(UFS+n) = merge(flux_tmp(URHO)*0.5d0*(qtempl(R_Y+n) + qtempr(R_Y+n)), flux_tmp(UFS+n), ustar .eq. 0.d0)
             enddo

             flux_tmp(UTEMP) = 0.0
             do n = UFX, UFX+naux
                flux_tmp(n) = merge(0.d0, flux_tmp(n), naux .gt. 0)
             enddo
             do n = UFA, UFA+nadv
                flux_tmp(n) = merge(0.d0, flux_tmp(n), nadv .gt. 0)
             enddo

             do ivar = 1, NVAR
                flux3(i,j,k,ivar) = flux3(i,j,k,ivar) + flux_tmp(ivar) * az(i,j,k)
             enddo
          enddo
       enddo
    enddo
    !$acc end kernels

    !$acc parallel loop gang vector collapse(4) present(lo,hi,d,flux1,flux2,flux3,v)
    do ivar=1,NVAR
       do k = lo(3)-nextra+1, hi(3)+nextra-1
          do j = lo(2)-nextra+1, hi(2)+nextra-1
             do i = lo(1)-nextra+1, hi(1)+nextra-1
                d(i,j,k,ivar) = - (flux1(i+1,j,k,ivar) - flux1(i,j,k,ivar) &
                                +  flux2(i,j+1,k,ivar) - flux2(i,j,k,ivar) &
                                +  flux3(i,j,k+1,ivar) - flux3(i,j,k,ivar)) / v(i,j,k)
             enddo
          enddo
       enddo
    enddo
    !$acc end parallel
    !$acc exit data delete(dqx,dqy,dqz) copyout(flux1,flux2,flux3,d) copyout(v,ax,ay,az,q,flatn,qd_lo,qd_hi,qt_lo,qt_hi,domlo,domhi,qaux,qa_lo,qa_hi,flag,fglo,fghi,lo,hi)

    deallocate(dqx)
    deallocate(dqy)
    deallocate(dqz)

    !do ivar=1,NVAR
    !   do k = lo(3)-nextra+1, hi(3)+nextra-1
    !      do j = lo(2)-nextra+1, hi(2)+nextra-1
    !         do i = lo(1)-nextra+1, hi(1)+nextra-1
    !            print *, "FLUX3: ", flux3(i,j,k,ivar)
    !         enddo
    !      enddo
    !   enddo
    !enddo

  end subroutine pc_hyp_mol_flux
end module hyp_advection_module 
