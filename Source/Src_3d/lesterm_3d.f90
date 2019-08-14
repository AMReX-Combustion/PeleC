! compute the LES source terms and fluxes for all the
! conservative equations.
!
! AT THE MOMENT:
! - no multifluid component of the equations are solved
!
module lesterm_module

  implicit none

  private :: get_sfs_stresses_xdir,&
       get_sfs_stresses_ydir,&
       get_sfs_stresses_zdir

  public :: pc_smagorinsky_sfs_term,&
       pc_dynamic_smagorinsky_sfs_term,&
       pc_dynamic_smagorinsky_quantities,&
       pc_dynamic_smagorinsky_coeffs

contains

  subroutine pc_smagorinsky_sfs_term(lo,  hi,&
       dmnlo, dmnhi,&
       Q,   Qlo,   Qhi,&
       tx,  txlo,  txhi,&
       Ax,  Axlo,  Axhi,&
       fx,  fxlo,  fxhi,&
       ty,  tylo,  tyhi,&
       Ay,  Aylo,  Ayhi,&
       fy,  fylo,  fyhi,&
       tz,  tzlo,  tzhi,&
       Az,  Azlo,  Azhi,&
       fz,  fzlo,  fzhi,&
       V,   Vlo,   Vhi,&
       L,   Llo,   Lhi,&
       deltax) bind(C, name = "pc_smagorinsky_sfs_term")

    use network, only   : nspecies
    use meth_params_module, only : NVAR, UMX, UMY, UMZ, UEDEN, QVAR, QRHO, QU, QW, QTEMP, QFS, Cs, CI, PrT
    use amrex_constants_module
    use eos_type_module
    use eos_module
    use prob_params_module, only : physbc_lo, physbc_hi, Inflow

    implicit none

    integer, intent(in) ::     lo(3),    hi(3)
    integer, intent(in) ::  dmnlo(3), dmnhi(3)
    integer, intent(in) ::    Qlo(3),   Qhi(3)
    integer, intent(in) ::   txlo(3),  txhi(3)
    integer, intent(in) ::   Axlo(3),  Axhi(3)
    integer, intent(in) ::   fxlo(3),  fxhi(3)
    integer, intent(in) ::   tylo(3),  tyhi(3)
    integer, intent(in) ::   Aylo(3),  Ayhi(3)
    integer, intent(in) ::   fylo(3),  fyhi(3)
    integer, intent(in) ::   tzlo(3),  tzhi(3)
    integer, intent(in) ::   Azlo(3),  Azhi(3)
    integer, intent(in) ::   fzlo(3),  fzhi(3)
    integer, intent(in) ::    Vlo(3),   Vhi(3)
    integer, intent(in) ::    Llo(3),   Lhi(3)

    double precision, intent(in   ) ::    Q(   Qlo(1):   Qhi(1),   Qlo(2):   Qhi(2),   Qlo(3):   Qhi(3), QVAR)
    double precision, intent(in   ) ::   tx(  txlo(1):  txhi(1),  txlo(2):  txhi(2),  txlo(3):  txhi(3), 6)
    double precision, intent(in   ) ::   Ax(  Axlo(1):  Axhi(1),  Axlo(2):  Axhi(2),  Axlo(3):  Axhi(3))
    double precision, intent(inout) ::   fx(  fxlo(1):  fxhi(1),  fxlo(2):  fxhi(2),  fxlo(3):  fxhi(3), NVAR)
    double precision, intent(in   ) ::   ty(  tylo(1):  tyhi(1),  tylo(2):  tyhi(2),  tylo(3):  tyhi(3), 6)
    double precision, intent(in   ) ::   Ay(  Aylo(1):  Ayhi(1),  Aylo(2):  Ayhi(2),  Aylo(3):  Ayhi(3))
    double precision, intent(inout) ::   fy(  fylo(1):  fyhi(1),  fylo(2):  fyhi(2),  fylo(3):  fyhi(3), NVAR)
    double precision, intent(in   ) ::   tz(  tzlo(1):  tzhi(1),  tzlo(2):  tzhi(2),  tzlo(3):  tzhi(3), 6)
    double precision, intent(in   ) ::   Az(  Azlo(1):  Azhi(1),  Azlo(2):  Azhi(2),  Azlo(3):  Azhi(3))
    double precision, intent(inout) ::   fz(  fzlo(1):  fzhi(1),  fzlo(2):  fzhi(2),  fzlo(3):  fzhi(3), NVAR)
    double precision, intent(inout) ::    L(   Llo(1):   Lhi(1),   Llo(2):   Lhi(2),   Llo(3):   Lhi(3), NVAR)
    double precision, intent(in   ) ::    V(   Vlo(1):   Vhi(1),   Vlo(2):   Vhi(2),   Vlo(3):   Vhi(3))
    double precision, intent(in   ) :: deltax(3)

    integer :: i, j, k, n
    double precision :: sigmaxx, sigmaxy, sigmaxz, sigmayx, sigmayy, sigmayz, sigmazx, sigmazy, sigmazz
    double precision :: alphaij(3,3), alpha(3), flux_T(3)
    double precision :: Uface(3)
    double precision :: gfaci(lo(1):hi(1)+1)
    double precision :: gfacj(lo(2):hi(2)+1)
    double precision :: gfack(lo(3):hi(3)+1)
    double precision :: dxinv(3)
    double precision :: deltabar
    double precision :: Cs2
    type(eos_t) :: sfs_eos_state

    dxinv = 1.d0/deltax
    Cs2 = Cs**2
    call build(sfs_eos_state)

    gfaci = dxinv(1)
    if (lo(1).le.dmnlo(1) .and. physbc_lo(1).eq.Inflow) gfaci(dmnlo(1)) = gfaci(dmnlo(1)) * TWO
    if (hi(1).gt.dmnhi(1) .and. physbc_hi(1).eq.Inflow) gfaci(dmnhi(1)+1) = gfaci(dmnhi(1)+1) * TWO
    gfacj = dxinv(2)
    if (lo(2).le.dmnlo(2) .and. physbc_lo(2).eq.Inflow) gfacj(dmnlo(2)) = gfacj(dmnlo(2)) * TWO
    if (hi(2).gt.dmnhi(2) .and. physbc_hi(2).eq.Inflow) gfacj(dmnhi(2)+1) = gfacj(dmnhi(2)+1) * TWO
    gfack = dxinv(3)
    if (lo(3).le.dmnlo(3) .and. physbc_lo(3).eq.Inflow) gfack(dmnlo(3)) = gfack(dmnlo(3)) * TWO
    if (hi(3).gt.dmnhi(3) .and. physbc_hi(3).eq.Inflow) gfack(dmnhi(3)+1) = gfack(dmnhi(3)+1) * TWO

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)+1
             ! SFS stress
             deltabar = deltax(1)
             call get_sfs_stresses_xdir(i, j, k, Q, Qlo, Qhi, tx, txlo, txhi, alphaij(1,1), alphaij(1,2), alphaij(1,3), alpha(1), flux_T(1), gfaci(i), deltabar)
             sigmaxx = Cs2 * alphaij(1,1) - THIRD * CI * alpha(1)
             sigmaxy = Cs2 * alphaij(1,2)
             sigmaxz = Cs2 * alphaij(1,3)
             fx(i,j,k,UMX)   = - sigmaxx
             fx(i,j,k,UMY)   = - sigmaxy
             fx(i,j,k,UMZ)   = - sigmaxz

             ! SFS turbulent diffusion
             Uface(:) = HALF*(Q(i,j,k,QU:QW) + Q(i-1,j,k,QU:QW))
             fx(i,j,k,UEDEN) = - sigmaxx*Uface(1) - sigmaxy*Uface(2) - sigmaxz*Uface(3)

             ! SFS heat flux
             sfs_eos_state % massfrac(:) = Q(i,j,k,QFS:QFS+nspecies-1)
             sfs_eos_state % T           = Q(i,j,k,QTEMP)
             call eos_cv(sfs_eos_state)
             fx(i,j,k,UEDEN) = fx(i,j,k,UEDEN) - sfs_eos_state%gam1 * sfs_eos_state%cv * Cs2 / PrT * flux_T(1)
          end do
       end do
    end do

    ! Scale fluxes by area
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)+1
             fx(i,j,k,UMX)   = fx(i,j,k,UMX)   * Ax(i,j,k)
             fx(i,j,k,UMY)   = fx(i,j,k,UMY)   * Ax(i,j,k)
             fx(i,j,k,UMZ)   = fx(i,j,k,UMZ)   * Ax(i,j,k)
             fx(i,j,k,UEDEN) = fx(i,j,k,UEDEN) * Ax(i,j,k)
          enddo
       enddo
    enddo

    do k=lo(3),hi(3)
       do i=lo(1),hi(1)

          do j=lo(2),hi(2)+1
             ! SFS stress
             deltabar = deltax(2)
             call get_sfs_stresses_ydir(i, j, k, Q, Qlo, Qhi, ty, tylo, tyhi, alphaij(2,1), alphaij(2,2), alphaij(2,3), alpha(2), flux_T(2), gfacj(j), deltabar)
             sigmayx = Cs2 * alphaij(2,1)
             sigmayy = Cs2 * alphaij(2,2) - THIRD * CI * alpha(2)
             sigmayz = Cs2 * alphaij(2,3)
             fy(i,j,k,UMX)   = - sigmayx
             fy(i,j,k,UMY)   = - sigmayy
             fy(i,j,k,UMZ)   = - sigmayz

             ! SFS turbulent diffusion
             Uface(:) = HALF*(Q(i,j,k,QU:QW) + Q(i,j-1,k,QU:QW))
             fy(i,j,k,UEDEN) = - sigmayx*Uface(1) - sigmayy*Uface(2) - sigmayz*Uface(3)

             ! SFS heat flux
             sfs_eos_state % massfrac(:) = Q(i,j,k,QFS:QFS+nspecies-1)
             sfs_eos_state % T           = Q(i,j,k,QTEMP)
             call eos_cv(sfs_eos_state)
             fy(i,j,k,UEDEN) = fy(i,j,k,UEDEN) - sfs_eos_state%gam1 * sfs_eos_state%cv * Cs2 / PrT * flux_T(2)
          end do
       end do
    end do

    ! Scale fluxes by area
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)
             fy(i,j,k,UMX)   = fy(i,j,k,UMX)   * Ay(i,j,k)
             fy(i,j,k,UMY)   = fy(i,j,k,UMY)   * Ay(i,j,k)
             fy(i,j,k,UMZ)   = fy(i,j,k,UMZ)   * Ay(i,j,k)
             fy(i,j,k,UEDEN) = fy(i,j,k,UEDEN) * Ay(i,j,k)
          end do
       end do
    end do

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)

          do k=lo(3),hi(3)+1
             ! SFS stress
             deltabar = deltax(3)
             call get_sfs_stresses_zdir(i, j, k, Q, Qlo, Qhi, tz, tzlo, tzhi, alphaij(3,1), alphaij(3,2), alphaij(3,3), alpha(3), flux_T(3), gfack(k), deltabar)
             sigmazx = Cs2 * alphaij(3,1)
             sigmazy = Cs2 * alphaij(3,2)
             sigmazz = Cs2 * alphaij(3,3) - THIRD * CI * alpha(3)
             fz(i,j,k,UMX)   = - sigmazx
             fz(i,j,k,UMY)   = - sigmazy
             fz(i,j,k,UMZ)   = - sigmazz

             ! SFS turbulent diffusion
             Uface(:) = HALF*(Q(i,j,k,QU:QW) + Q(i,j,k-1,QU:QW))
             fz(i,j,k,UEDEN) = - sigmazx*Uface(1) - sigmazy*Uface(2) - sigmazz*Uface(3)

             ! SFS heat flux
             sfs_eos_state % massfrac(:) = Q(i,j,k,QFS:QFS+nspecies-1)
             sfs_eos_state % T           = Q(i,j,k,QTEMP)
             call eos_cv(sfs_eos_state)
             fz(i,j,k,UEDEN) = fz(i,j,k,UEDEN) - sfs_eos_state%gam1 * sfs_eos_state%cv * Cs2 / PrT * flux_T(3)
          end do
       end do
    end do

    ! Scale fluxes by area
    do k=lo(3),hi(3)+1
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             fz(i,j,k,UMX)   = fz(i,j,k,UMX)   * Az(i,j,k)
             fz(i,j,k,UMY)   = fz(i,j,k,UMY)   * Az(i,j,k)
             fz(i,j,k,UMZ)   = fz(i,j,k,UMZ)   * Az(i,j,k)
             fz(i,j,k,UEDEN) = fz(i,j,k,UEDEN) * Az(i,j,k)
          end do
       end do
    end do

    call destroy(sfs_eos_state)

    do n=1,NVAR
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                L(i,j,k,n) = - (fx(i+1,j,k,n)-fx(i,j,k,n) &
                               +fy(i,j+1,k,n)-fy(i,j,k,n) &
                               +fz(i,j,k+1,n)-fz(i,j,k,n) ) / V(i,j,k)
             end do
          end do
       end do
    end do

  end subroutine pc_smagorinsky_sfs_term


  subroutine pc_dynamic_smagorinsky_sfs_term(lo,  hi,&
                                             Q,   Qlo,   Qhi,&
                                             alphaij, alphaijlo, alphaijhi,&
                                             alpha, alphalo, alphahi,&
                                             flux_T, flux_Tlo, flux_Thi,&
                                             Cs2x, Cs2xlo, Cs2xhi, &
                                             CIx, CIxlo, CIxhi, &
                                             PrTx, PrTxlo, PrTxhi, &
                                             Ax,  Axlo,  Axhi,&
                                             fx,  fxlo,  fxhi,&
                                             Cs2y, Cs2ylo, Cs2yhi, &
                                             CIy, CIylo, CIyhi, &
                                             PrTy, PrTylo, PrTyhi, &
                                             Ay,  Aylo,  Ayhi,&
                                             fy,  fylo,  fyhi,&
                                             Cs2z, Cs2zlo, Cs2zhi, &
                                             CIz, CIzlo, CIzhi, &
                                             PrTz, PrTzlo, PrTzhi, &
                                             Az,  Azlo,  Azhi,&
                                             fz,  fzlo,  fzhi,&
                                             V,   Vlo,   Vhi,&
                                             L,   Llo,   Lhi,&
                                             deltax) bind(C, name = "pc_dynamic_smagorinsky_sfs_term")

    use network, only : nspecies
    use meth_params_module, only : NVAR, UMX, UMY, UMZ, UEDEN, QVAR, QRHO, QU, QW, QTEMP, QFS
    use amrex_constants_module
    use eos_type_module
    use eos_module

    implicit none

    integer, intent(in) ::        lo(3),         hi(3)
    integer, intent(in) ::       Qlo(3),        Qhi(3)
    integer, intent(in) :: alphaijlo(3),  alphaijhi(3)
    integer, intent(in) ::   alphalo(3),    alphahi(3)
    integer, intent(in) ::  flux_Tlo(3),   flux_Thi(3)
    integer, intent(in) ::    Cs2xlo(3),     Cs2xhi(3)
    integer, intent(in) ::     CIxlo(3),      CIxhi(3)
    integer, intent(in) ::    PrTxlo(3),     PrTxhi(3)
    integer, intent(in) ::      Axlo(3),       Axhi(3)
    integer, intent(in) ::      fxlo(3),       fxhi(3)
    integer, intent(in) ::    Cs2ylo(3),     Cs2yhi(3)
    integer, intent(in) ::     CIylo(3),      CIyhi(3)
    integer, intent(in) ::    PrTylo(3),     PrTyhi(3)
    integer, intent(in) ::      Aylo(3),       Ayhi(3)
    integer, intent(in) ::      fylo(3),       fyhi(3)
    integer, intent(in) ::    Cs2zlo(3),     Cs2zhi(3)
    integer, intent(in) ::     CIzlo(3),      CIzhi(3)
    integer, intent(in) ::    PrTzlo(3),     PrTzhi(3)
    integer, intent(in) ::      Azlo(3),       Azhi(3)
    integer, intent(in) ::      fzlo(3),       fzhi(3)
    integer, intent(in) ::       Vlo(3),        Vhi(3)
    integer, intent(in) ::       Llo(3),        Lhi(3)

    double precision, intent(in   ) :: Q       (  Qlo       (1):Qhi       (1),  Qlo       (2):Qhi       (2),  Qlo       (3):Qhi       (3), QVAR)
    double precision, intent(in   ) :: alphaij (  alphaijlo (1):alphaijhi (1),  alphaijlo (2):alphaijhi (2),  alphaijlo (3):alphaijhi (3), 3*3)
    double precision, intent(in   ) :: alpha   (  alphalo   (1):alphahi   (1),  alphalo   (2):alphahi   (2),  alphalo   (3):alphahi   (3), 3)
    double precision, intent(in   ) :: flux_T  (  flux_Tlo  (1):flux_Thi  (1),  flux_Tlo  (2):flux_Thi  (2),  flux_Tlo  (3):flux_Thi  (3), 3)
    double precision, intent(in   ) :: Cs2x    (  Cs2xlo    (1):Cs2xhi    (1),  Cs2xlo    (2):Cs2xhi    (2),  Cs2xlo    (3):Cs2xhi    (3))
    double precision, intent(in   ) :: CIx     (  CIxlo     (1):CIxhi     (1),  CIxlo     (2):CIxhi     (2),  CIxlo     (3):CIxhi     (3))
    double precision, intent(in   ) :: PrTx    (  PrTxlo    (1):PrTxhi    (1),  PrTxlo    (2):PrTxhi    (2),  PrTxlo    (3):PrTxhi    (3))
    double precision, intent(in   ) :: Ax      (  Axlo      (1):Axhi      (1),  Axlo      (2):Axhi      (2),  Axlo      (3):Axhi      (3))
    double precision, intent(inout) :: fx      (  fxlo      (1):fxhi      (1),  fxlo      (2):fxhi      (2),  fxlo      (3):fxhi      (3), NVAR)
    double precision, intent(in   ) :: Cs2y    (  Cs2ylo    (1):Cs2yhi    (1),  Cs2ylo    (2):Cs2yhi    (2),  Cs2ylo    (3):Cs2yhi    (3))
    double precision, intent(in   ) :: CIy     (  CIylo     (1):CIyhi     (1),  CIylo     (2):CIyhi     (2),  CIylo     (3):CIyhi     (3))
    double precision, intent(in   ) :: PrTy    (  PrTylo    (1):PrTyhi    (1),  PrTylo    (2):PrTyhi    (2),  PrTylo    (3):PrTyhi    (3))
    double precision, intent(in   ) :: Ay      (  Aylo      (1):Ayhi      (1),  Aylo      (2):Ayhi      (2),  Aylo      (3):Ayhi      (3))
    double precision, intent(inout) :: fy      (  fylo      (1):fyhi      (1),  fylo      (2):fyhi      (2),  fylo      (3):fyhi      (3), NVAR)
    double precision, intent(in   ) :: Cs2z    (  Cs2zlo    (1):Cs2zhi    (1),  Cs2zlo    (2):Cs2zhi    (2),  Cs2zlo    (3):Cs2zhi    (3))
    double precision, intent(in   ) :: CIz     (  CIzlo     (1):CIzhi     (1),  CIzlo     (2):CIzhi     (2),  CIzlo     (3):CIzhi     (3))
    double precision, intent(in   ) :: PrTz    (  PrTzlo    (1):PrTzhi    (1),  PrTzlo    (2):PrTzhi    (2),  PrTzlo    (3):PrTzhi    (3))
    double precision, intent(in   ) :: Az      (  Azlo      (1):Azhi      (1),  Azlo      (2):Azhi      (2),  Azlo      (3):Azhi      (3))
    double precision, intent(inout) :: fz      (  fzlo      (1):fzhi      (1),  fzlo      (2):fzhi      (2),  fzlo      (3):fzhi      (3), NVAR)
    double precision, intent(inout) :: L       (  Llo       (1):Lhi       (1),  Llo       (2):Lhi       (2),  Llo       (3):Lhi       (3), NVAR)
    double precision, intent(in   ) :: V       (  Vlo       (1):Vhi       (1),  Vlo       (2):Vhi       (2),  Vlo       (3):Vhi       (3))
    double precision, intent(in   ) :: deltax(3)


    integer :: i, j, k, n
    double precision :: sigmaxx, sigmaxy, sigmaxz, sigmayx, sigmayy, sigmayz, sigmazx, sigmazy, sigmazz
    double precision :: Uface(3)
    double precision :: dxinv(3)
    type(eos_t) :: sfs_eos_state

    integer :: i11, i12, i13, i21, i22, i23, i31, i32, i33
    i11 = (1-1) * 3 + 1 ! index for entry (1,1)
    i12 = (2-1) * 3 + 1 ! index for entry (1,2)
    i13 = (3-1) * 3 + 1 ! index for entry (1,3)
    i21 = (1-1) * 3 + 2 ! index for entry (2,1)
    i22 = (2-1) * 3 + 2 ! index for entry (2,2)
    i23 = (3-1) * 3 + 2 ! index for entry (2,3)
    i31 = (1-1) * 3 + 3 ! index for entry (3,1)
    i32 = (2-1) * 3 + 3 ! index for entry (3,2)
    i33 = (3-1) * 3 + 3 ! index for entry (3,3)

    dxinv = 1.d0/deltax
    call build(sfs_eos_state)

    do k= lo(3), hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)+1
             ! SFS stress
             sigmaxx = Cs2x(i,j,k) * alphaij(i,j,k,i11) - THIRD * CIx(i,j,k) * alpha(i,j,k,1)
             sigmaxy = Cs2x(i,j,k) * alphaij(i,j,k,i12)
             sigmaxz = Cs2x(i,j,k) * alphaij(i,j,k,i13)
             fx(i,j,k,UMX)   = - sigmaxx
             fx(i,j,k,UMY)   = - sigmaxy
             fx(i,j,k,UMZ)   = - sigmaxz

             ! SFS turbulent diffusion
             Uface(:) = HALF*(Q(i,j,k,QU:QW) + Q(i-1,j,k,QU:QW))
             fx(i,j,k,UEDEN) = - sigmaxx*Uface(1) - sigmaxy*Uface(2) - sigmaxz*Uface(3)

             ! SFS heat flux
             sfs_eos_state % massfrac(:) = Q(i,j,k,QFS:QFS+nspecies-1)
             sfs_eos_state % T           = Q(i,j,k,QTEMP)
             call eos_cv(sfs_eos_state)
             fx(i,j,k,UEDEN) = fx(i,j,k,UEDEN) - sfs_eos_state%gam1 * sfs_eos_state%cv * Cs2x(i,j,k) / PrTx(i,j,k) * flux_T(i,j,k,1)
          end do
       end do
    end do

    ! Scale fluxes by area
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)+1
             fx(i,j,k,UMX)   = fx(i,j,k,UMX)   * Ax(i,j,k)
             fx(i,j,k,UMY)   = fx(i,j,k,UMY)   * Ax(i,j,k)
             fx(i,j,k,UMZ)   = fx(i,j,k,UMZ)   * Ax(i,j,k)
             fx(i,j,k,UEDEN) = fx(i,j,k,UEDEN) * Ax(i,j,k)
          end do
       end do
    end do

    do k=lo(3),hi(3)
       do i=lo(1),hi(1)
          do j=lo(2),hi(2)+1
             ! SFS stress
             sigmayx = Cs2y(i,j,k) * alphaij(i,j,k,i21)
             sigmayy = Cs2y(i,j,k) * alphaij(i,j,k,i22) - THIRD * CIy(i,j,k) * alpha(i,j,k,2)
             sigmayz = Cs2y(i,j,k) * alphaij(i,j,k,i23)
             fy(i,j,k,UMX)   = - sigmayx
             fy(i,j,k,UMY)   = - sigmayy
             fy(i,j,k,UMZ)   = - sigmayz

             ! SFS turbulent diffusion
             Uface(:) = HALF*(Q(i,j,k,QU:QW) + Q(i,j-1,k,QU:QW))
             fy(i,j,k,UEDEN) = - sigmayx*Uface(1) - sigmayy*Uface(2) - sigmayz*Uface(3)

             ! SFS heat flux
             sfs_eos_state % massfrac(:) = Q(i,j,k,QFS:QFS+nspecies-1)
             sfs_eos_state % T           = Q(i,j,k,QTEMP)
             call eos_cv(sfs_eos_state)
             fy(i,j,k,UEDEN) = fy(i,j,k,UEDEN) - sfs_eos_state%gam1 * sfs_eos_state%cv * Cs2y(i,j,k) / PrTy(i,j,k) * flux_T(i,j,k,2)
          end do
       end do
    end do

    ! Scale fluxes by area
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)
             fy(i,j,k,UMX)   = fy(i,j,k,UMX)   * Ay(i,j,k)
             fy(i,j,k,UMY)   = fy(i,j,k,UMY)   * Ay(i,j,k)
             fy(i,j,k,UMZ)   = fy(i,j,k,UMZ)   * Ay(i,j,k)
             fy(i,j,k,UEDEN) = fy(i,j,k,UEDEN) * Ay(i,j,k)
          end do
       end do
    end do

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          do k=lo(3),hi(3)+1
             ! SFS stress
             sigmazx = Cs2z(i,j,k) * alphaij(i,j,k,i31)
             sigmazy = Cs2z(i,j,k) * alphaij(i,j,k,i32)
             sigmazz = Cs2z(i,j,k) * alphaij(i,j,k,i33) - THIRD * CIz(i,j,k) * alpha(i,j,k,3)
             fz(i,j,k,UMX)   = - sigmazx
             fz(i,j,k,UMY)   = - sigmazy
             fz(i,j,k,UMZ)   = - sigmazz

             ! SFS turbulent diffusion
             Uface(:) = HALF*(Q(i,j,k,QU:QW) + Q(i,j,k-1,QU:QW))
             fz(i,j,k,UEDEN) = - sigmazx*Uface(1) - sigmazy*Uface(2) - sigmazz*Uface(3)

             ! SFS heat flux
             sfs_eos_state % massfrac(:) = Q(i,j,k,QFS:QFS+nspecies-1)
             sfs_eos_state % T           = Q(i,j,k,QTEMP)
             call eos_cv(sfs_eos_state)
             fz(i,j,k,UEDEN) = fz(i,j,k,UEDEN) - sfs_eos_state%gam1 * sfs_eos_state%cv * Cs2z(i,j,k) / PrTz(i,j,k) * flux_T(i,j,k,3)
          end do
       end do
    end do

    ! Scale fluxes by area
    do k=lo(3),hi(3)+1
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             fz(i,j,k,UMX)   = fz(i,j,k,UMX)   * Az(i,j,k)
             fz(i,j,k,UMY)   = fz(i,j,k,UMY)   * Az(i,j,k)
             fz(i,j,k,UMZ)   = fz(i,j,k,UMZ)   * Az(i,j,k)
             fz(i,j,k,UEDEN) = fz(i,j,k,UEDEN) * Az(i,j,k)
          end do
       end do
    end do

    call destroy(sfs_eos_state)

    do n=1,NVAR
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                L(i,j,k,n) = - (fx(i+1,j,k,n)-fx(i,j,k,n) &
                               +fy(i,j+1,k,n)-fy(i,j,k,n) &
                               +fz(i,j,k+1,n)-fz(i,j,k,n) ) / V(i,j,k)
             end do
          end do
       end do
    end do

  end subroutine pc_dynamic_smagorinsky_sfs_term


  subroutine pc_dynamic_smagorinsky_quantities(lo, hi,&
                                               dmnlo, dmnhi,&
                                               Q, Qlo, Qhi,&
                                               tx,  txlo,  txhi,&
                                               ty,  tylo,  tyhi,&
                                               tz,  tzlo,  tzhi,&
                                               Kij, Kijlo, Kijhi,&
                                               RUT, RUTlo, RUThi,&
                                               alphaij, alphaijlo, alphaijhi,&
                                               alpha, alphalo, alphahi,&
                                               flux_T, flux_Tlo, flux_Thi,&
                                               fgr,&
                                               deltax) bind(C, name = "pc_dynamic_smagorinsky_quantities")

    use meth_params_module, only : QVAR, QRHO, QU, QV, QW, QTEMP
    use amrex_constants_module
    use prob_params_module, only : physbc_lo, physbc_hi, Inflow

    implicit none

    integer, intent(in) :: lo       (3), hi       (3)
    integer, intent(in) :: dmnlo    (3), dmnhi    (3)
    integer, intent(in) :: Qlo      (3), Qhi      (3)
    integer, intent(in) :: txlo     (3), txhi     (3)
    integer, intent(in) :: tylo     (3), tyhi     (3)
    integer, intent(in) :: tzlo     (3), tzhi     (3)
    integer, intent(in) :: Kijlo    (3), Kijhi    (3)
    integer, intent(in) :: RUTlo    (3), RUThi    (3)
    integer, intent(in) :: alphaijlo(3), alphaijhi(3)
    integer, intent(in) :: alphalo  (3), alphahi  (3)
    integer, intent(in) :: flux_Tlo (3), flux_Thi (3)
    integer, intent(in) :: fgr

    double precision, intent(in   ) :: Q      (  Qlo      (1):Qhi      (1),  Qlo     (2): Qhi     (2),  Qlo       (3):Qhi      (3), QVAR)
    double precision, intent(in   ) :: tx     (  txlo     (1):txhi     (1),  txlo    (2): txhi    (2),  txlo      (3):txhi     (3), 6)
    double precision, intent(in   ) :: ty     (  tylo     (1):tyhi     (1),  tylo    (2): tyhi    (2),  tylo      (3):tyhi     (3), 6)
    double precision, intent(in   ) :: tz     (  tzlo     (1):tzhi     (1),  tzlo    (2): tzhi    (2),  tzlo      (3):tzhi     (3), 6)
    double precision, intent(inout) :: Kij    (  Kijlo    (1):Kijhi    (1),  Kijlo    (2):Kijhi    (2),  Kijlo    (3):Kijhi    (3), 6)
    double precision, intent(inout) :: RUT    (  RUTlo    (1):RUThi    (1),  RUTlo    (2):RUThi    (2),  RUTlo    (3):RUThi    (3), 3)
    double precision, intent(inout) :: alphaij(  alphaijlo(1):alphaijhi(1),  alphaijlo(2):alphaijhi(2),  alphaijlo(3):alphaijhi(3), 3*3)
    double precision, intent(inout) :: alpha  (  alphalo  (1):alphahi  (1),  alphalo  (2):alphahi  (2),  alphalo  (3):alphahi  (3), 3)
    double precision, intent(inout) :: flux_T (  flux_Tlo (1):flux_Thi (1),  flux_Tlo (2):flux_Thi (2),  flux_Tlo (3):flux_Thi (3), 3)
    double precision, intent(in   ) :: deltax(3)

    integer :: i, j, k
    double precision :: dxinv(3)
    double precision :: deltabar
    double precision :: gfaci(lo(1):hi(1)+1)
    double precision :: gfacj(lo(2):hi(2)+1)
    double precision :: gfack(lo(3):hi(3)+1)

    integer :: i11, i12, i13, i21, i22, i23, i31, i32, i33
    i11 = (1-1) * 3 + 1 ! index for entry (1,1)
    i12 = (2-1) * 3 + 1 ! index for entry (1,2)
    i13 = (3-1) * 3 + 1 ! index for entry (1,3)
    i21 = (1-1) * 3 + 2 ! index for entry (2,1)
    i22 = (2-1) * 3 + 2 ! index for entry (2,2)
    i23 = (3-1) * 3 + 2 ! index for entry (2,3)
    i31 = (1-1) * 3 + 3 ! index for entry (3,1)
    i32 = (2-1) * 3 + 3 ! index for entry (3,2)
    i33 = (3-1) * 3 + 3 ! index for entry (3,3)

    dxinv = 1.d0/deltax

    gfaci = dxinv(1)
    if (lo(1).le.dmnlo(1) .and. physbc_lo(1).eq.Inflow) gfaci(dmnlo(1)) = gfaci(dmnlo(1)) * TWO
    if (hi(1).gt.dmnhi(1) .and. physbc_hi(1).eq.Inflow) gfaci(dmnhi(1)+1) = gfaci(dmnhi(1)+1) * TWO
    gfacj = dxinv(2)
    if (lo(2).le.dmnlo(2) .and. physbc_lo(2).eq.Inflow) gfacj(dmnlo(2)) = gfacj(dmnlo(2)) * TWO
    if (hi(2).gt.dmnhi(2) .and. physbc_hi(2).eq.Inflow) gfacj(dmnhi(2)+1) = gfacj(dmnhi(2)+1) * TWO
    gfack = dxinv(3)
    if (lo(3).le.dmnlo(3) .and. physbc_lo(3).eq.Inflow) gfack(dmnlo(3)) = gfack(dmnlo(3)) * TWO
    if (hi(3).gt.dmnhi(3) .and. physbc_hi(3).eq.Inflow) gfack(dmnhi(3)+1) = gfack(dmnhi(3)+1) * TWO

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             Kij(i,j,k, 1) = Q(i,j,k,QRHO) * Q(i,j,k,QU) * Q(i,j,k,QU)
             Kij(i,j,k, 2) = Q(i,j,k,QRHO) * Q(i,j,k,QU) * Q(i,j,k,QV)
             Kij(i,j,k, 3) = Q(i,j,k,QRHO) * Q(i,j,k,QU) * Q(i,j,k,QW)
             Kij(i,j,k, 4) = Q(i,j,k,QRHO) * Q(i,j,k,QV) * Q(i,j,k,QV)
             Kij(i,j,k, 5) = Q(i,j,k,QRHO) * Q(i,j,k,QV) * Q(i,j,k,QW)
             Kij(i,j,k, 6) = Q(i,j,k,QRHO) * Q(i,j,k,QW) * Q(i,j,k,QW)

             RUT(i,j,k,:) = Q(i,j,k,QRHO) * Q(i,j,k,QU:QW) * Q(i,j,k,QTEMP)

             deltabar = fgr*deltax(1)
             call get_sfs_stresses_xdir(i, j, k, Q, Qlo, Qhi, tx, txlo, txhi, alphaij(i,j,k,i11), alphaij(i,j,k,i12), alphaij(i,j,k,i13), alpha(i,j,k,1), flux_T(i,j,k,1), gfaci(i), deltabar)

             deltabar = fgr*deltax(2)
             call get_sfs_stresses_ydir(i, j, k, Q, Qlo, Qhi, ty, tylo, tyhi, alphaij(i,j,k,i21), alphaij(i,j,k,i22), alphaij(i,j,k,i23), alpha(i,j,k,2), flux_T(i,j,k,2), gfacj(j), deltabar)

             deltabar = fgr*deltax(3)
             call get_sfs_stresses_zdir(i, j, k, Q, Qlo, Qhi, tz, tzlo, tzhi, alphaij(i,j,k,i31), alphaij(i,j,k,i32), alphaij(i,j,k,i33), alpha(i,j,k,3), flux_T(i,j,k,3), gfack(k), deltabar)
          end do
       end do
    end do

  end subroutine pc_dynamic_smagorinsky_quantities


  subroutine pc_dynamic_smagorinsky_coeffs(lo,  hi,&
                                           dmnlo, dmnhi,&
                                           Q,   Qlo,   Qhi,&
                                           tx,  txlo,  txhi,&
                                           ty,  tylo,  tyhi,&
                                           tz,  tzlo,  tzhi,&
                                           Kij, Kijlo, Kijhi,&
                                           RUT, RUTlo, RUThi,&
                                           alphaij, alphaijlo, alphaijhi,&
                                           alpha, alphalo, alphahi,&
                                           flux_T, flux_Tlo, flux_Thi,&
                                           Cs2, Cs2lo, Cs2hi, &
                                           CI, CIlo, CIhi, &
                                           PrT, PrTlo, PrThi, &
                                           fgr,&
                                           deltax) bind(C, name = "pc_dynamic_smagorinsky_coeffs")

    use meth_params_module, only : QVAR, QRHO, QU, QV, QW, QTEMP
    use amrex_constants_module
    use prob_params_module, only : physbc_lo, physbc_hi, Inflow

    implicit none

    integer, intent(in) :: lo       (3), hi       (3)
    integer, intent(in) :: dmnlo    (3), dmnhi    (3)
    integer, intent(in) :: Qlo      (3), Qhi      (3)
    integer, intent(in) :: txlo     (3), txhi     (3)
    integer, intent(in) :: tylo     (3), tyhi     (3)
    integer, intent(in) :: tzlo     (3), tzhi     (3)
    integer, intent(in) :: Kijlo    (3), Kijhi    (3)
    integer, intent(in) :: RUTlo    (3), RUThi    (3)
    integer, intent(in) :: alphaijlo(3), alphaijhi(3)
    integer, intent(in) :: alphalo  (3), alphahi  (3)
    integer, intent(in) :: flux_Tlo (3), flux_Thi (3)
    integer, intent(in) :: Cs2lo    (3), Cs2hi    (3)
    integer, intent(in) :: CIlo     (3), CIhi     (3)
    integer, intent(in) :: PrTlo    (3), PrThi    (3)
    integer, intent(in) :: fgr

    double precision, intent(in   ) :: Q      (  Qlo      (1):Qhi      (1),   Qlo     (2): Qhi     (2),   Qlo     (3): Qhi     (3), QVAR)
    double precision, intent(in   ) :: tx     (  txlo     (1):txhi     (1),  txlo    (2): txhi    (2),  txlo      (3):txhi     (3), 6)
    double precision, intent(in   ) :: ty     (  tylo     (1):tyhi     (1),  tylo    (2): tyhi    (2),  tylo      (3):tyhi     (3), 6)
    double precision, intent(in   ) :: tz     (  tzlo     (1):tzhi     (1),  tzlo    (2): tzhi    (2),  tzlo      (3):tzhi     (3), 6)
    double precision, intent(in   ) :: Kij    (  Kijlo    (1):Kijhi    (1),  Kijlo    (2):Kijhi    (2),  Kijlo    (3):Kijhi    (3), 6)
    double precision, intent(in   ) :: RUT    (  RUTlo    (1):RUThi    (1),  RUTlo    (2):RUThi    (2),  RUTlo    (3):RUThi    (3), 3)
    double precision, intent(in   ) :: alphaij(  alphaijlo(1):alphaijhi(1),  alphaijlo(2):alphaijhi(2),  alphaijlo(3):alphaijhi(3), 3*3)
    double precision, intent(in   ) :: alpha  (  alphalo  (1):alphahi  (1),  alphalo  (2):alphahi  (2),  alphalo  (3):alphahi  (3), 3)
    double precision, intent(in   ) :: flux_T (  flux_Tlo (1):flux_Thi (1),  flux_Tlo (2):flux_Thi (2),  flux_Tlo (3):flux_Thi (3), 3)
    double precision, intent(inout) :: Cs2    (  Cs2lo    (1):Cs2hi    (1),  Cs2lo    (2):Cs2hi    (2),  Cs2lo    (3):Cs2hi    (3))
    double precision, intent(inout) :: CI     (  CIlo     (1):CIhi     (1),  CIlo     (2):CIhi     (2),  CIlo     (3):CIhi     (3))
    double precision, intent(inout) :: PrT    (  PrTlo    (1):PrThi    (1),  PrTlo    (2):PrThi    (2),  PrTlo    (3):PrThi    (3))
    double precision, intent(in   ) :: deltax(3)

    double precision :: L(3*3), M(3*3), betaij(3*3), beta(3), T(3), KE(3)
    double precision :: LM
    double precision :: MM
    double precision :: Lkk
    double precision :: bma
    double precision :: TT
    double precision :: KT

    integer :: i, j, k
    double precision :: dxinv(3)
    double precision :: deltahat
    double precision, parameter :: small_num = 1d-8
    double precision :: gfaci(lo(1):hi(1)+1)
    double precision :: gfacj(lo(2):hi(2)+1)
    double precision :: gfack(lo(3):hi(3)+1)
    integer :: i11, i12, i13, i21, i22, i23, i31, i32, i33
    i11 = (1-1) * 3 + 1 ! index for entry (1,1)
    i12 = (2-1) * 3 + 1 ! index for entry (1,2)
    i13 = (3-1) * 3 + 1 ! index for entry (1,3)
    i21 = (1-1) * 3 + 2 ! index for entry (2,1)
    i22 = (2-1) * 3 + 2 ! index for entry (2,2)
    i23 = (3-1) * 3 + 2 ! index for entry (2,3)
    i31 = (1-1) * 3 + 3 ! index for entry (3,1)
    i32 = (2-1) * 3 + 3 ! index for entry (3,2)
    i33 = (3-1) * 3 + 3 ! index for entry (3,3)

    dxinv = 1.d0/deltax

    gfaci = dxinv(1)
    if (lo(1).le.dmnlo(1) .and. physbc_lo(1).eq.Inflow) gfaci(dmnlo(1)) = gfaci(dmnlo(1)) * TWO
    if (hi(1).gt.dmnhi(1) .and. physbc_hi(1).eq.Inflow) gfaci(dmnhi(1)+1) = gfaci(dmnhi(1)+1) * TWO
    gfacj = dxinv(2)
    if (lo(2).le.dmnlo(2) .and. physbc_lo(2).eq.Inflow) gfacj(dmnlo(2)) = gfacj(dmnlo(2)) * TWO
    if (hi(2).gt.dmnhi(2) .and. physbc_hi(2).eq.Inflow) gfacj(dmnhi(2)+1) = gfacj(dmnhi(2)+1) * TWO
    gfack = dxinv(3)
    if (lo(3).le.dmnlo(3) .and. physbc_lo(3).eq.Inflow) gfack(dmnlo(3)) = gfack(dmnlo(3)) * TWO
    if (hi(3).gt.dmnhi(3) .and. physbc_hi(3).eq.Inflow) gfack(dmnhi(3)+1) = gfack(dmnhi(3)+1) * TWO

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             ! Make betaij and beta
             deltahat = fgr*deltax(1)
             call get_sfs_stresses_xdir(i, j, k, Q, Qlo, Qhi, tx, txlo, txhi, betaij(i11), betaij(i12), betaij(i13), beta(1), T(1), gfaci(i), deltahat)
             T(1) = flux_T(i,j,k,1) - T(1)

             deltahat = fgr*deltax(2)
             call get_sfs_stresses_ydir(i, j, k, Q, Qlo, Qhi, ty, tylo, tyhi, betaij(i21), betaij(i22), betaij(i23), beta(2), T(2), gfacj(j), deltahat)
             T(2) = flux_T(i,j,k,2) - T(2)

             deltahat = fgr*deltax(3)
             call get_sfs_stresses_zdir(i, j, k, Q, Qlo, Qhi, tz, tzlo, tzhi, betaij(i31), betaij(i32), betaij(i33), beta(3), T(3), gfack(k), deltahat)
             T(3) = flux_T(i,j,k,3) - T(3)

             ! "resolved turbulent stresses" and others
             M(:) = betaij(:) - alphaij(i,j,k,:)
             L(i11) = Kij(i,j,k, 1) - Q(i,j,k,QRHO) * Q(i,j,k,QU) * Q(i,j,k,QU)
             L(i12) = Kij(i,j,k, 2) - Q(i,j,k,QRHO) * Q(i,j,k,QU) * Q(i,j,k,QV)
             L(i13) = Kij(i,j,k, 3) - Q(i,j,k,QRHO) * Q(i,j,k,QU) * Q(i,j,k,QW)
             L(i21) = L(i12)
             L(i22) = Kij(i,j,k, 4) - Q(i,j,k,QRHO) * Q(i,j,k,QV) * Q(i,j,k,QV)
             L(i23) = Kij(i,j,k, 5) - Q(i,j,k,QRHO) * Q(i,j,k,QV) * Q(i,j,k,QW)
             L(i31) = L(i13)
             L(i32) = L(i23)
             L(i33) = Kij(i,j,k, 6) - Q(i,j,k,QRHO) * Q(i,j,k,QW) * Q(i,j,k,QW)
             KE(:) = RUT(i,j,k,:) - Q(i,j,k,QRHO) * Q(i,j,k,QU:QW) * Q(i,j,k,QTEMP)

             ! Contractions
             LM = sum(L(:) * M(:)) + small_num
             MM = sum(M(:) * M(:)) + small_num
             Lkk = (L(i11) + L(i22) + L(i33)) + small_num
             bma = THIRD*sum(beta(:) - alpha(i,j,k,:)) + small_num
             TT = sum(T(:)*T(:)) + small_num
             KT = sum(KE(:)*T(:)) + small_num

             ! Coefficients
             Cs2(i,j,k) = max(LM / MM, small_num)
             CI(i,j,k) = max(Lkk / bma, small_num)
             PrT(i,j,k) = max(TT / KT, small_num)
          end do
       end do
    end do

    ! scale Prandtl with Cs2
    PrT(:,:,:) = Cs2(:,:,:) * PrT(:,:,:)

  end subroutine pc_dynamic_smagorinsky_coeffs


  subroutine get_sfs_stresses_xdir(i, j, k, Q, Qlo, Qhi, tx, txlo, txhi, alphaij_xx, alphaij_xy, alphaij_xz, alpha, flux_T, gfac, deltabar)

    use meth_params_module, only : QVAR, QRHO, QU, QV, QW, QTEMP
    use amrex_constants_module

    implicit none

    integer, intent(in) :: i, j, k
    integer, intent(in) :: Qlo(3), Qhi(3)
    integer, intent(in) :: txlo(3), txhi(3)
    double precision, intent(in   ) :: Q( Qlo(1):Qhi(1),  Qlo(2): Qhi(2),  Qlo(3):Qhi(3), QVAR)
    double precision, intent(in   ) :: tx( txlo(1):txhi(1), txlo(2):txhi(2), txlo(3):txhi(3), 6)
    double precision, intent(inout) :: alphaij_xx, alphaij_xy, alphaij_xz
    double precision, intent(inout) :: alpha
    double precision, intent(inout) :: flux_T
    double precision, intent(in   ) :: gfac
    double precision, intent(in   ) :: deltabar

    ! local
    double precision :: Sijmag, Skk, mut
    double precision :: dUdx(3,3), S(3,3)
    double precision :: dTdx

    dUdx(1,1) = gfac * (Q(i,j,k,QU)    - Q(i-1,j,k,QU))
    dUdx(1,2) = tx(i,j,k,1)
    dUdx(1,3) = tx(i,j,k,4)
    dUdx(2,1) = gfac * (Q(i,j,k,QV)    - Q(i-1,j,k,QV))
    dUdx(2,2) = tx(i,j,k,2)
    dUdx(2,3) = tx(i,j,k,5)
    dUdx(3,1) = gfac * (Q(i,j,k,QW)    - Q(i-1,j,k,QW))
    dUdx(3,2) = tx(i,j,k,3)
    dUdx(3,3) = tx(i,j,k,6)

    S(:,:) = HALF * (dUdx(:,:) + transpose(dUdx(:,:)))
    Skk = S(1,1) + S(2,2) + S(3,3)
    Sijmag = sqrt(TWO * sum(S(:,:)**2))
    mut = Q(i,j,k,QRHO) * deltabar**2 * Sijmag

    alphaij_xx = TWO * mut * ( S(1,1) - THIRD * Skk )
    alphaij_xy = TWO * mut * S(1,2)
    alphaij_xz = TWO * mut * S(1,3)
    alpha      = TWO * mut * Sijmag

    dTdx = gfac * (Q(i,j,k,QTEMP) - Q(i-1,j,k,QTEMP))
    flux_T = mut * dTdx

  end subroutine get_sfs_stresses_xdir


  subroutine get_sfs_stresses_ydir(i, j, k, Q, Qlo, Qhi, ty, tylo, tyhi, alphaij_yx, alphaij_yy, alphaij_yz, alpha, flux_T, gfac, deltabar)

    use meth_params_module, only : QVAR, QRHO, QU, QV, QW, QTEMP
    use amrex_constants_module

    implicit none

    integer, intent(in) :: i, j, k
    integer, intent(in) :: Qlo(3), Qhi(3)
    integer, intent(in) :: tylo(3), tyhi(3)
    double precision, intent(in   ) :: Q( Qlo(1):Qhi(1),  Qlo(2): Qhi(2),  Qlo(3):Qhi(3), QVAR)
    double precision, intent(in   ) :: ty( tylo(1):tyhi(1), tylo(2):tyhi(2), tylo(3):tyhi(3), 6)
    double precision, intent(inout) :: alphaij_yx, alphaij_yy, alphaij_yz
    double precision, intent(inout) :: alpha
    double precision, intent(inout) :: flux_T
    double precision, intent(in   ) :: gfac
    double precision, intent(in   ) :: deltabar

    ! local
    double precision :: Sijmag, Skk, mut
    double precision :: dUdx(3,3), S(3,3)
    double precision :: dTdy

    dUdx(1,1) = ty(i,j,k,1)
    dUdx(1,2) = gfac * (Q(i,j,k,QU)    - Q(i,j-1,k,QU))
    dUdx(1,3) = ty(i,j,k,4)
    dUdx(2,1)= ty(i,j,k,2)
    dudx(2,2) = gfac * (Q(i,j,k,QV)    - Q(i,j-1,k,QV))
    dUdx(2,3)= ty(i,j,k,5)
    dUdx(3,1) = ty(i,j,k,3)
    dUdx(3,2) = gfac * (Q(i,j,k,QW)    - Q(i,j-1,k,QW))
    dUdx(3,3) = ty(i,j,k,6)

    S(:,:) = HALF * (dUdx(:,:) + transpose(dUdx(:,:)))
    Skk = S(1,1) + S(2,2) + S(3,3)
    Sijmag = sqrt(TWO * sum(S(:,:)**2))
    mut = Q(i,j,k,QRHO) * deltabar**2 * Sijmag

    alphaij_yx = TWO * mut * S(2,1)
    alphaij_yy = TWO * mut * ( S(2,2) - THIRD * Skk)
    alphaij_yz = TWO * mut * S(2,3)
    alpha      = TWO * mut * Sijmag

    dTdy = gfac * (Q(i,j,k,QTEMP) - Q(i,j-1,k,QTEMP))
    flux_T = mut * dTdy

  end subroutine get_sfs_stresses_ydir


  subroutine get_sfs_stresses_zdir(i, j, k, Q, Qlo, Qhi, tz, tzlo, tzhi, alphaij_zx, alphaij_zy, alphaij_zz, alpha, flux_T, gfac, deltabar)

    use meth_params_module, only : QVAR, QRHO, QU, QV, QW, QTEMP
    use amrex_constants_module

    implicit none

    integer, intent(in) :: i, j, k
    integer, intent(in) :: Qlo(3), Qhi(3)
    integer, intent(in) :: tzlo(3), tzhi(3)
    double precision, intent(in   ) :: Q( Qlo(1):Qhi(1),  Qlo(2): Qhi(2),  Qlo(3):Qhi(3), QVAR)
    double precision, intent(in   ) :: tz( tzlo(1):tzhi(1), tzlo(2):tzhi(2), tzlo(3):tzhi(3), 6)
    double precision, intent(inout) :: alphaij_zx, alphaij_zy, alphaij_zz
    double precision, intent(inout) :: alpha
    double precision, intent(inout) :: flux_T
    double precision, intent(in   ) :: gfac
    double precision, intent(in   ) :: deltabar

    ! local
    double precision :: Sijmag, Skk, mut
    double precision :: dUdx(3,3), S(3,3)
    double precision :: dTdz

    dUdx(1,1) = tz(i,j,k,1)
    dUdx(1,2) = tz(i,j,k,4)
    dUdx(1,3) = gfac * (Q(i,j,k,QU)    - Q(i,j,k-1,QU))
    dUdx(2,1) = tz(i,j,k,2)
    dUdx(2,2) = tz(i,j,k,5)
    dUdx(2,3) = gfac * (Q(i,j,k,QV)    - Q(i,j,k-1,QV))
    dUdx(3,1) = tz(i,j,k,3)
    dUdx(3,2) = tz(i,j,k,6)
    dUdx(3,3) = gfac * (Q(i,j,k,QW)    - Q(i,j,k-1,QW))

    S(:,:) = HALF * (dUdx(:,:) + transpose(dUdx(:,:)))
    Skk = S(1,1) + S(2,2) + S(3,3)
    Sijmag = sqrt(TWO * sum(S(:,:)**2))
    mut = Q(i,j,k,QRHO) * deltabar**2 * Sijmag

    alphaij_zx = TWO * mut * S(3,1)
    alphaij_zy = TWO * mut * S(3,2)
    alphaij_zz = TWO * mut * ( S(3,3) - THIRD * Skk)
    alpha      = TWO * mut * Sijmag

    dTdz = gfac * (Q(i,j,k,QTEMP) - Q(i,j,k-1,QTEMP))
    flux_T = mut * dTdz

  end subroutine get_sfs_stresses_zdir

end module lesterm_module
