! compute the LES source terms and fluxes for all the
! conservative equations.
!
! AT THE MOMENT:
! - no multifluid component of the equations are solved
!
module lesterm_module

  implicit none

  private :: get_sfs_stresses_xdir,&
       get_sfs_stresses_ydir

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
                                     V,   Vlo,   Vhi,&
                                     L,   Llo,   Lhi,&
                                     deltax) bind(C, name = "pc_smagorinsky_sfs_term")

    use network, only : nspecies
    use meth_params_module, only : NVAR, UMX, UMY, UEDEN, QVAR, QRHO, QU, QV, QTEMP, QFS, Cs, CI, PrT
    use amrex_constants_module
    use eos_type_module
    use eos_module
    use prob_params_module, only : physbc_lo, physbc_hi, Inflow

    implicit none

    integer, intent(in) ::     lo(2),    hi(2)
    integer, intent(in) ::  dmnlo(2), dmnhi(2)
    integer, intent(in) ::    Qlo(2),   Qhi(2)
    integer, intent(in) ::   txlo(2),  txhi(2)
    integer, intent(in) ::   Axlo(2),  Axhi(2)
    integer, intent(in) ::   fxlo(2),  fxhi(2)
    integer, intent(in) ::   tylo(2),  tyhi(2)
    integer, intent(in) ::   Aylo(2),  Ayhi(2)
    integer, intent(in) ::   fylo(2),  fyhi(2)
    integer, intent(in) ::    Vlo(2),   Vhi(2)
    integer, intent(in) ::    Llo(2),   Lhi(2)

    double precision, intent(in   ) ::    Q(   Qlo(1):   Qhi(1),   Qlo(2):   Qhi(2), QVAR)
    double precision, intent(in   ) ::   tx(  txlo(1):  txhi(1),  txlo(2):  txhi(2), 2)
    double precision, intent(in   ) ::   Ax(  Axlo(1):  Axhi(1),  Axlo(2):  Axhi(2))
    double precision, intent(inout) ::   fx(  fxlo(1):  fxhi(1),  fxlo(2):  fxhi(2), NVAR)
    double precision, intent(in   ) ::   ty(  tylo(1):  tyhi(1),  tylo(2):  tyhi(2), 2)
    double precision, intent(in   ) ::   Ay(  Aylo(1):  Ayhi(1),  Aylo(2):  Ayhi(2))
    double precision, intent(inout) ::   fy(  fylo(1):  fyhi(1),  fylo(2):  fyhi(2), NVAR)
    double precision, intent(inout) ::    L(   Llo(1):   Lhi(1),   Llo(2):   Lhi(2), NVAR)
    double precision, intent(in   ) ::    V(   Vlo(1):   Vhi(1),   Vlo(2):   Vhi(2))
    double precision, intent(in   ) :: deltax(2)

    integer :: i, j, n
    double precision :: sigmaxx, sigmaxy, sigmayx, sigmayy
    double precision :: alphaij(2,2), alpha(2), flux_T(2)
    double precision :: Uface(2)
    double precision :: gfaci(lo(1):hi(1)+1)
    double precision :: gfacj(lo(2):hi(2)+1)
    double precision :: dxinv(2)
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

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)+1
          ! SFS stress
          deltabar = deltax(1)
          call get_sfs_stresses_xdir(i, j, Q, Qlo, Qhi, tx, txlo, txhi, alphaij(1,1), alphaij(1,2), alpha(1), flux_T(1), gfaci(i), deltabar)
          sigmaxx = Cs2 * alphaij(1,1) - THIRD * CI * alpha(1)
          sigmaxy = Cs2 * alphaij(1,2)
          fx(i,j,UMX)   = - sigmaxx
          fx(i,j,UMY)   = - sigmaxy

          ! SFS turbulent diffusion
          Uface(:) = HALF*(Q(i,j,QU:QV) + Q(i-1,j,QU:QV))
          fx(i,j,UEDEN) = - sigmaxx*Uface(1) - sigmaxy*Uface(2)

          ! SFS heat flux
          sfs_eos_state % massfrac(:) = Q(i,j,QFS:QFS+nspecies-1)
          sfs_eos_state % T           = Q(i,j,QTEMP)
          call eos_cv(sfs_eos_state)
          fx(i,j,UEDEN) = fx(i,j,UEDEN) - sfs_eos_state%gam1 * sfs_eos_state%cv * Cs2 / PrT * flux_T(1)
       end do
    end do

    ! Scale fluxes by area
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)+1
          fx(i,j,UMX)   = fx(i,j,UMX)   * Ax(i,j)
          fx(i,j,UMY)   = fx(i,j,UMY)   * Ax(i,j)
          fx(i,j,UEDEN) = fx(i,j,UEDEN) * Ax(i,j)
       end do
    end do

    do i=lo(1),hi(1)
       do j=lo(2),hi(2)+1
          ! SFS stress
          deltabar = deltax(2)
          call get_sfs_stresses_ydir(i, j, Q, Qlo, Qhi, ty, tylo, tyhi, alphaij(2,1), alphaij(2,2), alpha(2), flux_T(2), gfacj(j), deltabar)
          sigmayx = Cs2 * alphaij(2,1)
          sigmayy = Cs2 * alphaij(2,2) - THIRD * CI * alpha(2)
          fy(i,j,UMX)   = - sigmayx
          fy(i,j,UMY)   = - sigmayy

          ! SFS turbulent diffusion
          Uface(:) = HALF*(Q(i,j,QU:QV) + Q(i,j-1,QU:QV))
          fy(i,j,UEDEN) = - sigmayx*Uface(1) - sigmayy*Uface(2)

          ! SFS heat flux
          sfs_eos_state % massfrac(:) = Q(i,j,QFS:QFS+nspecies-1)
          sfs_eos_state % T           = Q(i,j,QTEMP)
          call eos_cv(sfs_eos_state)
          fy(i,j,UEDEN) = fy(i,j,UEDEN) - sfs_eos_state%gam1 * sfs_eos_state%cv * Cs2 / PrT * flux_T(2)
       end do
    end do

    ! Scale fluxes by area
    do j=lo(2),hi(2)+1
       do i=lo(1),hi(1)
          fy(i,j,UMX)   = fy(i,j,UMX)   * Ay(i,j)
          fy(i,j,UMY)   = fy(i,j,UMY)   * Ay(i,j)
          fy(i,j,UEDEN) = fy(i,j,UEDEN) * Ay(i,j)
       end do
    end do

    call destroy(sfs_eos_state)

    do n=1,NVAR
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             L(i,j,n) = - (fx(i+1,j,n) - fx(i,j,n) &
                         + fy(i,j+1,n) - fy(i,j,n) ) / V(i,j)
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
                                             V,   Vlo,   Vhi,&
                                             L,   Llo,   Lhi,&
                                             deltax) bind(C, name = "pc_dynamic_smagorinsky_sfs_term")

    use network, only : nspecies
    use meth_params_module, only : NVAR, UMX, UMY, UEDEN, QVAR, QRHO, QU, QV, QTEMP, QFS
    use amrex_constants_module
    use eos_type_module
    use eos_module

    implicit none

    integer, intent(in) ::        lo(2),        hi(2)
    integer, intent(in) ::       Qlo(2),       Qhi(2)
    integer, intent(in) :: alphaijlo(2), alphaijhi(2)
    integer, intent(in) ::   alphalo(2),   alphahi(2)
    integer, intent(in) :: flux_Tlo(2),   flux_Thi(2)
    integer, intent(in) ::   Cs2xlo(2),     Cs2xhi(2)
    integer, intent(in) ::    CIxlo(2),      CIxhi(2)
    integer, intent(in) ::   PrTxlo(2),     PrTxhi(2)
    integer, intent(in) ::     Axlo(2),       Axhi(2)
    integer, intent(in) ::     fxlo(2),       fxhi(2)
    integer, intent(in) ::   Cs2ylo(2),     Cs2yhi(2)
    integer, intent(in) ::    CIylo(2),      CIyhi(2)
    integer, intent(in) ::   PrTylo(2),     PrTyhi(2)
    integer, intent(in) ::     Aylo(2),       Ayhi(2)
    integer, intent(in) ::     fylo(2),       fyhi(2)
    integer, intent(in) ::      Vlo(2),        Vhi(2)
    integer, intent(in) ::      Llo(2),        Lhi(2)

    double precision, intent(in   ) :: Q       (  Qlo       (1):Qhi       (1),  Qlo       (2):Qhi       (2), QVAR)
    double precision, intent(in   ) :: alphaij (  alphaijlo (1):alphaijhi (1),  alphaijlo (2):alphaijhi (2), 2*2)
    double precision, intent(in   ) :: alpha   (  alphalo   (1):alphahi   (1),  alphalo   (2):alphahi   (2), 2)
    double precision, intent(in   ) :: flux_T  (  flux_Tlo  (1):flux_Thi  (1),  flux_Tlo  (2):flux_Thi  (2), 2)
    double precision, intent(in   ) :: Cs2x    (  Cs2xlo    (1):Cs2xhi    (1),  Cs2xlo    (2):Cs2xhi    (2))
    double precision, intent(in   ) :: CIx     (  CIxlo     (1):CIxhi     (1),  CIxlo     (2):CIxhi     (2))
    double precision, intent(in   ) :: PrTx    (  PrTxlo    (1):PrTxhi    (1),  PrTxlo    (2):PrTxhi    (2))
    double precision, intent(in   ) :: Ax      (  Axlo      (1):Axhi      (1),  Axlo      (2):Axhi      (2))
    double precision, intent(inout) :: fx      (  fxlo      (1):fxhi      (1),  fxlo      (2):fxhi      (2), NVAR)
    double precision, intent(in   ) :: Cs2y    (  Cs2ylo    (1):Cs2yhi    (1),  Cs2ylo    (2):Cs2yhi    (2))
    double precision, intent(in   ) :: CIy     (  CIylo     (1):CIyhi     (1),  CIylo     (2):CIyhi     (2))
    double precision, intent(in   ) :: PrTy    (  PrTylo    (1):PrTyhi    (1),  PrTylo    (2):PrTyhi    (2))
    double precision, intent(in   ) :: Ay      (  Aylo      (1):Ayhi      (1),  Aylo      (2):Ayhi      (2))
    double precision, intent(inout) :: fy      (  fylo      (1):fyhi      (1),  fylo      (2):fyhi      (2), NVAR)
    double precision, intent(inout) :: L       (  Llo       (1):Lhi       (1),  Llo       (2):Lhi       (2), NVAR)
    double precision, intent(in   ) :: V       (  Vlo       (1):Vhi       (1),  Vlo       (2):Vhi       (2))
    double precision, intent(in   ) :: deltax(2)


    integer :: i, j, n
    double precision :: sigmaxx, sigmaxy, sigmayx, sigmayy
    double precision :: Uface(2)
    double precision :: dxinv(2)
    type(eos_t) :: sfs_eos_state

    integer :: i11, i12, i21, i22
    i11 = (1-1) * 2 + 1 ! index for entry (1,1)
    i12 = (2-1) * 2 + 1 ! index for entry (1,2)
    i21 = (1-1) * 2 + 2 ! index for entry (2,1)
    i22 = (2-1) * 2 + 2 ! index for entry (2,2)

    dxinv = 1.d0/deltax
    call build(sfs_eos_state)

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)+1
          ! SFS stress
          sigmaxx = Cs2x(i,j) * alphaij(i,j,i11) - THIRD * CIx(i,j) * alpha(i,j,1)
          sigmaxy = Cs2x(i,j) * alphaij(i,j,i12)
          fx(i,j,UMX)   = - sigmaxx
          fx(i,j,UMY)   = - sigmaxy

          ! SFS turbulent diffusion
          Uface(:) = HALF*(Q(i,j,QU:QV) + Q(i-1,j,QU:QV))
          fx(i,j,UEDEN) = - sigmaxx*Uface(1) - sigmaxy*Uface(2)

          ! SFS heat flux
          sfs_eos_state % massfrac(:) = Q(i,j,QFS:QFS+nspecies-1)
          sfs_eos_state % T           = Q(i,j,QTEMP)
          call eos_cv(sfs_eos_state)
          fx(i,j,UEDEN) = fx(i,j,UEDEN) - sfs_eos_state%gam1 * sfs_eos_state%cv * Cs2x(i,j) / PrTx(i,j) * flux_T(i,j,1)
       end do
    end do

    ! Scale fluxes by area
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)+1
          fx(i,j,UMX)   = fx(i,j,UMX)   * Ax(i,j)
          fx(i,j,UMY)   = fx(i,j,UMY)   * Ax(i,j)
          fx(i,j,UEDEN) = fx(i,j,UEDEN) * Ax(i,j)
       end do
    end do

    do i=lo(1),hi(1)
       do j=lo(2),hi(2)+1
          ! SFS stres
          sigmayx = Cs2y(i,j) * alphaij(i,j,i21)
          sigmayy = Cs2y(i,j) * alphaij(i,j,i22) - THIRD * CIy(i,j) * alpha(i,j,2)
          fy(i,j,UMX)   = - sigmayx
          fy(i,j,UMY)   = - sigmayy

          ! SFS turbulent diffusion
          Uface(:) = HALF*(Q(i,j,QU:QV) + Q(i,j-1,QU:QV))
          fy(i,j,UEDEN) = - sigmayx*Uface(1) - sigmayy*Uface(2)

          ! SFS heat flux
          sfs_eos_state % massfrac(:) = Q(i,j,QFS:QFS+nspecies-1)
          sfs_eos_state % T           = Q(i,j,QTEMP)
          call eos_cv(sfs_eos_state)
          fy(i,j,UEDEN) = fy(i,j,UEDEN) - sfs_eos_state%gam1 * sfs_eos_state%cv * Cs2y(i,j) / PrTy(i,j) * flux_T(i,j,2)
       end do
    end do

    ! Scale fluxes by area
    do j=lo(2),hi(2)+1
       do i=lo(1),hi(1)
          fy(i,j,UMX)   = fy(i,j,UMX)   * Ay(i,j)
          fy(i,j,UMY)   = fy(i,j,UMY)   * Ay(i,j)
          fy(i,j,UEDEN) = fy(i,j,UEDEN) * Ay(i,j)
       end do
    end do

    call destroy(sfs_eos_state)

    do n=1,NVAR
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             L(i,j,n) = - (fx(i+1,j,n) - fx(i,j,n) &
                         + fy(i,j+1,n) - fy(i,j,n) ) / V(i,j)
          end do
       end do
    end do

  end subroutine pc_dynamic_smagorinsky_sfs_term


  subroutine pc_dynamic_smagorinsky_quantities(lo, hi,&
                                               dmnlo, dmnhi,&
                                               Q, Qlo, Qhi,&
                                               tx,  txlo,  txhi,&
                                               ty,  tylo,  tyhi,&
                                               Kij, Kijlo, Kijhi,&
                                               RUT, RUTlo, RUThi,&
                                               alphaij, alphaijlo, alphaijhi,&
                                               alpha, alphalo, alphahi,&
                                               flux_T, flux_Tlo, flux_Thi,&
                                               fgr,&
                                               deltax) bind(C, name = "pc_dynamic_smagorinsky_quantities")

    use meth_params_module, only : QVAR, QRHO, QU, QV, QTEMP
    use amrex_constants_module
    use prob_params_module, only : physbc_lo, physbc_hi, Inflow

    implicit none

    integer, intent(in) :: lo       (2), hi       (2)
    integer, intent(in) :: dmnlo    (2), dmnhi    (2)
    integer, intent(in) :: Qlo      (2), Qhi      (2)
    integer, intent(in) :: txlo     (2), txhi     (2)
    integer, intent(in) :: tylo     (2), tyhi     (2)
    integer, intent(in) :: Kijlo    (2), Kijhi    (2)
    integer, intent(in) :: RUTlo    (2), RUThi    (2)
    integer, intent(in) :: alphaijlo(2), alphaijhi(2)
    integer, intent(in) :: alphalo  (2), alphahi  (2)
    integer, intent(in) :: flux_Tlo (2), flux_Thi (2)
    integer, intent(in) :: fgr

    double precision, intent(in   ) :: Q      (  Qlo      (1):Qhi      (1),  Qlo      (2):Qhi      (2), QVAR)
    double precision, intent(in   ) :: tx     (  txlo     (1):txhi     (1),  txlo     (2):txhi     (2), 2)
    double precision, intent(in   ) :: ty     (  tylo     (1):tyhi     (1),  tylo     (2):tyhi     (2), 2)
    double precision, intent(inout) :: Kij    (  Kijlo    (1):Kijhi    (1),  Kijlo    (2):Kijhi    (2), 3)
    double precision, intent(inout) :: RUT    (  RUTlo    (1):RUThi    (1),  RUTlo    (2):RUThi    (2), 2)
    double precision, intent(inout) :: alphaij(  alphaijlo(1):alphaijhi(1),  alphaijlo(2):alphaijhi(2), 2*2)
    double precision, intent(inout) :: alpha  (  alphalo  (1):alphahi  (1),  alphalo  (2):alphahi  (2), 2)
    double precision, intent(inout) :: flux_T (  flux_Tlo (1):flux_Thi (1),  flux_Tlo (2):flux_Thi (2), 2)
    double precision, intent(in   ) :: deltax(2)

    integer :: i, j
    double precision :: dxinv(2)
    double precision :: deltabar
    double precision :: gfaci(lo(1):hi(1)+1)
    double precision :: gfacj(lo(2):hi(2)+1)

    integer :: i11, i12, i21, i22
    i11 = (1-1) * 2 + 1 ! index for entry (1,1)
    i12 = (2-1) * 2 + 1 ! index for entry (1,2)
    i21 = (1-1) * 2 + 2 ! index for entry (2,1)
    i22 = (2-1) * 2 + 2 ! index for entry (2,2)

    dxinv = 1.d0/deltax

    gfaci = dxinv(1)
    if (lo(1).le.dmnlo(1) .and. physbc_lo(1).eq.Inflow) gfaci(dmnlo(1)) = gfaci(dmnlo(1)) * TWO
    if (hi(1).gt.dmnhi(1) .and. physbc_hi(1).eq.Inflow) gfaci(dmnhi(1)+1) = gfaci(dmnhi(1)+1) * TWO
    gfacj = dxinv(2)
    if (lo(2).le.dmnlo(2) .and. physbc_lo(2).eq.Inflow) gfacj(dmnlo(2)) = gfacj(dmnlo(2)) * TWO
    if (hi(2).gt.dmnhi(2) .and. physbc_hi(2).eq.Inflow) gfacj(dmnhi(2)+1) = gfacj(dmnhi(2)+1) * TWO

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)

          Kij(i,j, 1) = Q(i,j,QRHO) * Q(i,j,QU) * Q(i,j,QU)
          Kij(i,j, 2) = Q(i,j,QRHO) * Q(i,j,QU) * Q(i,j,QV)
          Kij(i,j, 3) = Q(i,j,QRHO) * Q(i,j,QV) * Q(i,j,QV)

          RUT(i,j,:) = Q(i,j,QRHO) * Q(i,j,QU:QV) * Q(i,j,QTEMP)

          deltabar = fgr*deltax(1)
          call get_sfs_stresses_xdir(i, j, Q, Qlo, Qhi, tx, txlo, txhi, alphaij(i,j,i11), alphaij(i,j,i12), alpha(i,j,1), flux_T(i,j,1), gfaci(i), deltabar)

          deltabar = fgr*deltax(2)
          call get_sfs_stresses_ydir(i, j, Q, Qlo, Qhi, ty, tylo, tyhi, alphaij(i,j,i21), alphaij(i,j,i22), alpha(i,j,2), flux_T(i,j,2), gfacj(j), deltabar)

       end do
    end do

  end subroutine pc_dynamic_smagorinsky_quantities


  subroutine pc_dynamic_smagorinsky_coeffs(lo,  hi,&
                                           dmnlo, dmnhi,&
                                           Q,   Qlo,   Qhi,&
                                           tx,  txlo,  txhi,&
                                           ty,  tylo,  tyhi,&
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

    use meth_params_module, only : QVAR, QRHO, QU, QV, QTEMP
    use amrex_constants_module
    use prob_params_module, only : physbc_lo, physbc_hi, Inflow

    implicit none

    integer, intent(in) :: lo       (2), hi       (2)
    integer, intent(in) :: dmnlo    (2), dmnhi    (2)
    integer, intent(in) :: Qlo      (2), Qhi      (2)
    integer, intent(in) :: txlo     (2), txhi     (2)
    integer, intent(in) :: tylo     (2), tyhi     (2)
    integer, intent(in) :: Kijlo    (2), Kijhi    (2)
    integer, intent(in) :: RUTlo    (2), RUThi    (2)
    integer, intent(in) :: alphaijlo(2), alphaijhi(2)
    integer, intent(in) :: alphalo  (2), alphahi  (2)
    integer, intent(in) :: flux_Tlo (2), flux_Thi (2)
    integer, intent(in) :: Cs2lo    (2), Cs2hi    (2)
    integer, intent(in) :: CIlo     (2), CIhi     (2)
    integer, intent(in) :: PrTlo    (2), PrThi    (2)
    integer, intent(in) :: fgr

    double precision, intent(in   ) :: Q      (  Qlo      (1):Qhi      (1),   Qlo     (2): Qhi     (2), QVAR)
    double precision, intent(in   ) :: tx     (  txlo     (1):txhi     (1),  txlo     (2):txhi     (2), 2)
    double precision, intent(in   ) :: ty     (  tylo     (1):tyhi     (1),  tylo     (2):tyhi     (2), 2)
    double precision, intent(in   ) :: Kij    (  Kijlo    (1):Kijhi    (1),  Kijlo    (2):Kijhi    (2), 3)
    double precision, intent(in   ) :: RUT    (  RUTlo    (1):RUThi    (1),  RUTlo    (2):RUThi    (2), 2)
    double precision, intent(in   ) :: alphaij(  alphaijlo(1):alphaijhi(1),  alphaijlo(2):alphaijhi(2), 2*2)
    double precision, intent(in   ) :: alpha  (  alphalo  (1):alphahi  (1),  alphalo  (2):alphahi  (2), 2)
    double precision, intent(in   ) :: flux_T (  flux_Tlo (1):flux_Thi (1),  flux_Tlo (2):flux_Thi (2), 2)
    double precision, intent(inout) :: Cs2    (  Cs2lo    (1):Cs2hi    (1),  Cs2lo    (2):Cs2hi    (2))
    double precision, intent(inout) :: CI     (  CIlo     (1):CIhi     (1),  CIlo     (2):CIhi     (2))
    double precision, intent(inout) :: PrT    (  PrTlo    (1):PrThi    (1),  PrTlo    (2):PrThi    (2))
    double precision, intent(in   ) :: deltax(2)

    double precision :: L(2*2), M(2*2), betaij(2*2), beta(2), T(2), KE(2)
    double precision :: LM
    double precision :: MM
    double precision :: Lkk
    double precision :: bma
    double precision :: TT
    double precision :: KT

    integer :: i, j
    double precision :: dxinv(2)
    double precision :: deltahat
    double precision, parameter :: small_num = 1d-8
    double precision :: gfaci(lo(1):hi(1)+1)
    double precision :: gfacj(lo(2):hi(2)+1)
    integer :: i11, i12, i21, i22
    i11 = (1-1) * 2 + 1 ! index for entry (1,1)
    i12 = (2-1) * 2 + 1 ! index for entry (1,2)
    i21 = (1-1) * 2 + 2 ! index for entry (2,1)
    i22 = (2-1) * 2 + 2 ! index for entry (2,2)

    dxinv = 1.d0/deltax

    gfaci = dxinv(1)
    if (lo(1).le.dmnlo(1) .and. physbc_lo(1).eq.Inflow) gfaci(dmnlo(1)) = gfaci(dmnlo(1)) * TWO
    if (hi(1).gt.dmnhi(1) .and. physbc_hi(1).eq.Inflow) gfaci(dmnhi(1)+1) = gfaci(dmnhi(1)+1) * TWO
    gfacj = dxinv(2)
    if (lo(2).le.dmnlo(2) .and. physbc_lo(2).eq.Inflow) gfacj(dmnlo(2)) = gfacj(dmnlo(2)) * TWO
    if (hi(2).gt.dmnhi(2) .and. physbc_hi(2).eq.Inflow) gfacj(dmnhi(2)+1) = gfacj(dmnhi(2)+1) * TWO

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          ! Make betaij and beta
          deltahat = fgr*deltax(1)
          call get_sfs_stresses_xdir(i, j, Q, Qlo, Qhi, tx, txlo, txhi, betaij(i11), betaij(i12), beta(1), T(1), gfaci(i), deltahat)
          T(1) = flux_T(i,j,1) - T(1)

          deltahat = fgr*deltax(2)
          call get_sfs_stresses_ydir(i, j, Q, Qlo, Qhi, ty, tylo, tyhi, betaij(i21), betaij(i22), beta(2), T(2), gfacj(j), deltahat)
          T(2) = flux_T(i,j,2) - T(2)

          ! "resolved turbulent stresses" and others
          M(:) = betaij(:) - alphaij(i,j,:)
          L(i11) = Kij(i,j,1) - Q(i,j,QRHO) * Q(i,j,QU) * Q(i,j,QU)
          L(i12) = Kij(i,j,2) - Q(i,j,QRHO) * Q(i,j,QU) * Q(i,j,QV)
          L(i22) = Kij(i,j,3) - Q(i,j,QRHO) * Q(i,j,QV) * Q(i,j,QV)
          L(i21) = L(i12)
          KE(:) = RUT(i,j,:) - Q(i,j,QRHO) * Q(i,j,QU:QV) * Q(i,j,QTEMP)

          ! Contractions
          LM = sum(L(:) * M(:)) + small_num
          MM = sum(M(:) * M(:)) + small_num
          Lkk = (L(i11) + L(i22)) + small_num
          bma = HALF*sum(beta(:) - alpha(i,j,:)) + small_num
          TT = sum(T(:)*T(:)) + small_num
          KT = sum(KE(:)*T(:)) + small_num

          ! Coefficients
          Cs2(i,j) =  max(LM / MM, small_num)
          CI(i,j) =  max(Lkk / bma, small_num)
          PrT(i,j) =  max(TT / KT, small_num)
       end do
    end do

    ! scale Prandtl with Cs2
    PrT(:,:) = Cs2(:,:) * PrT(:,:)

  end subroutine pc_dynamic_smagorinsky_coeffs


  subroutine get_sfs_stresses_xdir(i, j, Q, Qlo, Qhi, tx, txlo, txhi, alphaij_xx, alphaij_xy, alpha, flux_T, gfac, deltabar)

    use meth_params_module, only : QVAR, QRHO, QU, QV, QTEMP
    use amrex_constants_module

    implicit none

    integer, intent(in) :: i, j
    integer, intent(in) :: Qlo(2), Qhi(2)
    integer, intent(in) :: txlo(2), txhi(2)
    double precision, intent(in   ) :: Q( Qlo(1):Qhi(1),  Qlo(2): Qhi(2), QVAR)
    double precision, intent(in   ) :: tx( txlo(1):txhi(1),  txlo(2): txhi(2), 2)
    double precision, intent(inout) :: alphaij_xx, alphaij_xy
    double precision, intent(inout) :: alpha
    double precision, intent(inout) :: flux_T
    double precision, intent(in   ) :: gfac
    double precision, intent(in   ) :: deltabar

    ! local
    double precision :: Sijmag, Skk, mut
    double precision :: dUdx(2,2), S(2,2)
    double precision :: dTdx

    ! dUdx(:,1) = dxinv(1)*(Q(i,j,QU:QV) - Q(i-1,j,QU:QV))
    ! dUdx(:,2) = FOURTH*dxinv(2)*(Q(i,j+1,QU:QV)+Q(i-1,j+1,QU:QV)-Q(i,j-1,QU:QV)-Q(i-1,j-1,QU:QV))
    dUdx(1,1) = gfac * (Q(i,j,QU)    - Q(i-1,j,QU))
    dUdx(1,2) = tx(i,j,1)
    dUdx(2,1) = gfac * (Q(i,j,QV)    - Q(i-1,j,QV))
    dUdx(2,2) = tx(i,j,2)

    S(:,:) = HALF * (dUdx(:,:) + transpose(dUdx(:,:)))
    Skk = S(1,1) + S(2,2)
    Sijmag = sqrt(TWO * sum(S(:,:)**2))
    mut = Q(i,j,QRHO) * deltabar**2 * Sijmag

    alphaij_xx = TWO * mut * ( S(1,1) - THIRD * Skk )
    alphaij_xy = TWO * mut * S(1,2)
    alpha      = TWO * mut * Sijmag

    dTdx = gfac * (Q(i,j,QTEMP) - Q(i-1,j,QTEMP))
    flux_T = mut * dTdx

  end subroutine get_sfs_stresses_xdir


  subroutine get_sfs_stresses_ydir(i, j, Q, Qlo, Qhi, ty, tylo, tyhi, alphaij_yx, alphaij_yy, alpha, flux_T, gfac, deltabar)

    use meth_params_module, only : QVAR, QRHO, QU, QV, QTEMP
    use amrex_constants_module

    implicit none

    integer, intent(in) :: i, j
    integer, intent(in) :: Qlo(2), Qhi(2)
    integer, intent(in) :: tylo(2), tyhi(2)
    double precision, intent(in   ) :: Q( Qlo(1):Qhi(1),  Qlo(2): Qhi(2), QVAR)
    double precision, intent(in   ) :: ty( tylo(1):tyhi(1),  tylo(2): tyhi(2), 2)
    double precision, intent(inout) :: alphaij_yx, alphaij_yy
    double precision, intent(inout) :: alpha
    double precision, intent(inout) :: flux_T
    double precision, intent(in   ) :: gfac
    double precision, intent(in   ) :: deltabar

    ! local
    double precision :: Sijmag, Skk, mut
    double precision :: dUdx(2,2), S(2,2)
    double precision :: dTdy

    ! dUdx(:,1) = FOURTH*dxinv(1)*(Q(i+1,j,QU:QV)+Q(i+1,j-1,QU:QV)-Q(i-1,j,QU:QV)-Q(i-1,j-1,QU:QV))
    ! dUdx(:,2) = dxinv(2)*(Q(i,j,QU:QV) - Q(i,j-1,QU:QV))
    dUdx(1,1) = ty(i,j,1)
    dUdx(1,2) = gfac * (Q(i,j,QU)    - Q(i,j-1,QU))
    dUdx(2,1) = ty(i,j,2)
    dUdx(2,2) = gfac * (Q(i,j,QV)    - Q(i,j-1,QV))

    S(:,:) = HALF * (dUdx(:,:) + transpose(dUdx(:,:)))
    Skk = S(1,1) + S(2,2)
    Sijmag = sqrt(TWO * sum(S(:,:)**2))
    mut = Q(i,j,QRHO) * deltabar**2 * Sijmag

    alphaij_yx = TWO * mut * S(2,1)
    alphaij_yy = TWO * mut * ( S(2,2) - THIRD * Skk)
    alpha      = TWO * mut * Sijmag

    dTdy = gfac * (Q(i,j,QTEMP) - Q(i,j-1,QTEMP))
    flux_T = mut * dTdy

  end subroutine get_sfs_stresses_ydir

end module lesterm_module
