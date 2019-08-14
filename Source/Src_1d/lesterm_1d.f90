! compute the LES source terms and fluxes for all the
! conservative equations.
!
! AT THE MOMENT:
! - no multifluid component of the equations are solved
!
module lesterm_module

  implicit none

  private :: get_sfs_stresses_xdir

  public :: pc_smagorinsky_sfs_term,&
       pc_dynamic_smagorinsky_sfs_term,&
       pc_dynamic_smagorinsky_quantities,&
       pc_dynamic_smagorinsky_coeffs

contains

  subroutine pc_smagorinsky_sfs_term(lo,  hi,&
                         dmnlo, dmnhi,&
                         Q,   Qlo,   Qhi,&
                         Ax,  Axlo,  Axhi,&
                         fx,  fxlo,  fxhi,&
                         V,   Vlo,   Vhi,&
                         L,   Llo,   Lhi,&
                         deltax) bind(C, name = "pc_smagorinsky_sfs_term")

    use network, only : nspecies
    use meth_params_module, only : NVAR, UMX, UEDEN, QVAR, QRHO, QU, QW, QTEMP, QFS, Cs, CI, PrT
    use amrex_constants_module
    use eos_type_module
    use eos_module
    use prob_params_module, only : physbc_lo, physbc_hi, Inflow

    implicit none

    integer, intent(in) ::     lo(1),    hi(1)
    integer, intent(in) ::  dmnlo(1), dmnhi(1)
    integer, intent(in) ::    Qlo(1),   Qhi(1)
    integer, intent(in) ::   Axlo(1),  Axhi(1)
    integer, intent(in) ::   fxlo(1),  fxhi(1)
    integer, intent(in) ::    Vlo(1),   Vhi(1)
    integer, intent(in) ::    Llo(1),   Lhi(1)

    double precision, intent(in   ) ::    Q(   Qlo(1):   Qhi(1), QVAR)
    double precision, intent(in   ) ::   Ax(  Axlo(1):  Axhi(1))
    double precision, intent(inout) ::   fx(  fxlo(1):  fxhi(1), NVAR)
    double precision, intent(inout) ::    L(   Llo(1):   Lhi(1), NVAR)
    double precision, intent(in   ) ::    V(   Vlo(1):   Vhi(1))
    double precision, intent(in   ) :: deltax(1)

    integer :: i, n
    double precision :: sigmaxx
    double precision :: alphaij, alpha, flux_T
    double precision :: uface
    double precision :: gfaci(lo(1):hi(1)+1)
    double precision :: dxinv(1)
    double precision :: deltabar
    double precision :: Cs2
    type(eos_t) :: sfs_eos_state

    dxinv = 1.d0/deltax
    Cs2 = Cs**2
    call build(sfs_eos_state)

    gfaci = dxinv(1)
    if (lo(1).le.dmnlo(1) .and. physbc_lo(1).eq.Inflow) gfaci(dmnlo(1)) = gfaci(dmnlo(1)) * TWO
    if (hi(1).gt.dmnhi(1) .and. physbc_hi(1).eq.Inflow) gfaci(dmnhi(1)+1) = gfaci(dmnhi(1)+1) * TWO

    do i=lo(1),hi(1)+1
       ! SFS stress
       deltabar = deltax(1)
       call get_sfs_stresses_xdir(i, Q, Qlo, Qhi, alphaij, alpha, flux_T, gfaci(i), deltabar)
       sigmaxx = Cs2 * alphaij - THIRD * CI * alpha
       fx(i,UMX)   = - sigmaxx

       ! SFS turbulent diffusion
       uface = HALF*(Q(i,QU)    + Q(i-1,QU))
       fx(i,UEDEN) = - sigmaxx*uface

       ! SFS heat flux
       sfs_eos_state % massfrac(:) = Q(i,QFS:QFS+nspecies-1)
       sfs_eos_state % T           = Q(i,QTEMP)
       call eos_cv(sfs_eos_state)
       fx(i,UEDEN) = fx(i,UEDEN) - sfs_eos_state%gam1 * sfs_eos_state%cv * Cs2 / PrT * flux_T
    end do

    ! Scale fluxes by area
    do i=lo(1),hi(1)+1
       fx(i,UMX)   = fx(i,UMX)   * Ax(i)
       fx(i,UEDEN) = fx(i,UEDEN) * Ax(i)
    enddo

    call destroy(sfs_eos_state)

    do n=1,NVAR
       do i = lo(1), hi(1)
          L(i,n) = - (fx(i+1,n)-fx(i,n)) / V(i)
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
                                             V,   Vlo,   Vhi,&
                                             L,   Llo,   Lhi,&
                                             deltax) bind(C, name = "pc_dynamic_smagorinsky_sfs_term")

    use network, only : nspecies
    use meth_params_module, only : NVAR, UMX, UEDEN, QVAR, QRHO, QU, QTEMP, QFS
    use amrex_constants_module
    use eos_type_module
    use eos_module

    implicit none

    integer, intent(in) ::        lo(1),        hi(1)
    integer, intent(in) ::       Qlo(1),       Qhi(1)
    integer, intent(in) :: alphaijlo(1), alphaijhi(1)
    integer, intent(in) ::   alphalo(1),   alphahi(1)
    integer, intent(in) :: flux_Tlo(1),   flux_Thi(1)
    integer, intent(in) ::    Cs2xlo(1),      Cs2xhi(1)
    integer, intent(in) ::    CIxlo(1),      CIxhi(1)
    integer, intent(in) ::   PrTxlo(1),     PrTxhi(1)
    integer, intent(in) ::     Axlo(1),       Axhi(1)
    integer, intent(in) ::     fxlo(1),       fxhi(1)
    integer, intent(in) ::      Vlo(1),        Vhi(1)
    integer, intent(in) ::      Llo(1),        Lhi(1)

    double precision, intent(in   ) :: Q       (  Qlo       (1):Qhi       (1), QVAR)
    double precision, intent(in   ) :: alphaij (  alphaijlo (1):alphaijhi (1), 1)
    double precision, intent(in   ) :: alpha   (  alphalo   (1):alphahi   (1), 1)
    double precision, intent(in   ) :: flux_T  (  flux_Tlo  (1):flux_Thi  (1), 1)
    double precision, intent(in   ) :: Cs2x    (  Cs2xlo    (1):Cs2xhi    (1))
    double precision, intent(in   ) :: CIx     (  CIxlo     (1):CIxhi     (1))
    double precision, intent(in   ) :: PrTx    (  PrTxlo    (1):PrTxhi    (1))
    double precision, intent(in   ) :: Ax      (  Axlo      (1):Axhi      (1))
    double precision, intent(inout) :: fx      (  fxlo      (1):fxhi      (1), NVAR)
    double precision, intent(inout) :: L       (  Llo       (1):Lhi       (1), NVAR)
    double precision, intent(in   ) :: V       (  Vlo       (1):Vhi       (1))
    double precision, intent(in   ) :: deltax(1)

    integer :: i, n
    double precision :: sigmaxx
    double precision :: uface
    double precision :: dxinv(1)
    type(eos_t) :: sfs_eos_state

    integer :: i11
    i11 = (1-1) * 1 + 1 ! index for entry (1,1)

    dxinv = 1.d0/deltax
    call build(sfs_eos_state)

    do i=lo(1),hi(1)+1
       ! SFS stress
       sigmaxx = Cs2x(i) * alphaij(i,i11) - THIRD * CIx(i) * alpha(i,1)
       fx(i,UMX)   = - sigmaxx

       ! SFS turbulent diffusion
       uface = HALF*(Q(i,QU)    + Q(i-1,QU))
       fx(i,UEDEN) = - sigmaxx*uface

       ! SFS heat flux
       sfs_eos_state % massfrac(:) = Q(i,QFS:QFS+nspecies-1)
       sfs_eos_state % T           = Q(i,QTEMP)
       call eos_cv(sfs_eos_state)
       fx(i,UEDEN) = fx(i,UEDEN) - sfs_eos_state%gam1 * sfs_eos_state%cv * Cs2x(i) / PrTx(i) * flux_T(i,1)
    end do

    ! Scale fluxes by area
    do i=lo(1),hi(1)+1
       fx(i,UMX)   = fx(i,UMX)   * Ax(i)
       fx(i,UEDEN) = fx(i,UEDEN) * Ax(i)
    enddo

    call destroy(sfs_eos_state)

    do n=1,NVAR
       do i = lo(1), hi(1)
          L(i,n) = - (fx(i+1,n)-fx(i,n)) / V(i)
       end do
    end do

  end subroutine pc_dynamic_smagorinsky_sfs_term


  subroutine pc_dynamic_smagorinsky_quantities(lo, hi,&
                                               dmnlo, dmnhi,&
                                               Q, Qlo, Qhi,&
                                               Kij, Kijlo, Kijhi,&
                                               RUT, RUTlo, RUThi,&
                                               alphaij, alphaijlo, alphaijhi,&
                                               alpha, alphalo, alphahi,&
                                               flux_T, flux_Tlo, flux_Thi,&
                                               fgr,&
                                               deltax) bind(C, name = "pc_dynamic_smagorinsky_quantities")

    use meth_params_module, only : QVAR, QRHO, QU, QTEMP
    use amrex_constants_module
    use prob_params_module, only : physbc_lo, physbc_hi, Inflow

    implicit none

    integer, intent(in) :: lo       (1), hi       (1)
    integer, intent(in) :: dmnlo    (1), dmnhi    (1)
    integer, intent(in) :: Qlo      (1), Qhi      (1)
    integer, intent(in) :: Kijlo    (1), Kijhi    (1)
    integer, intent(in) :: RUTlo    (1), RUThi    (1)
    integer, intent(in) :: alphaijlo(1), alphaijhi(1)
    integer, intent(in) :: alphalo  (1), alphahi  (1)
    integer, intent(in) :: flux_Tlo (1), flux_Thi (1)
    integer, intent(in) :: fgr

    double precision, intent(in   ) :: Q      (  Qlo      (1):Qhi      (1), QVAR)
    double precision, intent(inout) :: Kij    (  Kijlo    (1):Kijhi    (1), 1)
    double precision, intent(inout) :: RUT    (  RUTlo    (1):RUThi    (1), 1)
    double precision, intent(inout) :: alphaij(  alphaijlo(1):alphaijhi(1), 1)
    double precision, intent(inout) :: alpha  (  alphalo  (1):alphahi  (1), 1)
    double precision, intent(inout) :: flux_T (  flux_Tlo (1):flux_Thi (1), 1)
    double precision, intent(in   ) :: deltax(1)

    integer :: i
    double precision :: dxinv(1)
    double precision :: deltabar
    double precision :: gfaci(lo(1):hi(1)+1)

    integer :: i11
    i11 = (1-1) * 2 + 1 ! index for entry (1,1)

    dxinv = 1.d0/deltax

    gfaci = dxinv(1)
    if (lo(1).le.dmnlo(1) .and. physbc_lo(1).eq.Inflow) gfaci(dmnlo(1)) = gfaci(dmnlo(1)) * TWO
    if (hi(1).gt.dmnhi(1) .and. physbc_hi(1).eq.Inflow) gfaci(dmnhi(1)+1) = gfaci(dmnhi(1)+1) * TWO

    do i=lo(1),hi(1)

       Kij(i,1) = Q(i,QRHO) * Q(i,QU) * Q(i,QU)
       RUT(i,:) = Q(i,QRHO) * Q(i,QU) * Q(i,QTEMP)

       deltabar = fgr*deltax(1)
       call get_sfs_stresses_xdir(i, Q, Qlo, Qhi, alphaij(i,i11), alpha(i,1), flux_T(i,1), gfaci(i), deltabar)

    end do

  end subroutine pc_dynamic_smagorinsky_quantities


  subroutine pc_dynamic_smagorinsky_coeffs(lo,  hi,&
                                           dmnlo, dmnhi,&
                                           Q,   Qlo,   Qhi,&
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

    use meth_params_module, only : QVAR, QRHO, QU, QTEMP
    use amrex_constants_module
    use prob_params_module, only : physbc_lo, physbc_hi, Inflow

    implicit none

    integer, intent(in) :: lo       (1), hi       (1)
    integer, intent(in) :: dmnlo    (1), dmnhi    (1)
    integer, intent(in) :: Qlo      (1), Qhi      (1)
    integer, intent(in) :: Kijlo    (1), Kijhi    (1)
    integer, intent(in) :: RUTlo    (1), RUThi    (1)
    integer, intent(in) :: alphaijlo(1), alphaijhi(1)
    integer, intent(in) :: alphalo  (1), alphahi  (1)
    integer, intent(in) :: flux_Tlo (1), flux_Thi (1)
    integer, intent(in) :: Cs2lo    (1), Cs2hi    (1)
    integer, intent(in) :: CIlo     (1), CIhi     (1)
    integer, intent(in) :: PrTlo    (1), PrThi    (1)
    integer, intent(in) :: fgr

    double precision, intent(in   ) :: Q      (  Qlo      (1):Qhi      (1), QVAR)
    double precision, intent(in   ) :: Kij    (  Kijlo    (1):Kijhi    (1), 1)
    double precision, intent(in   ) :: RUT    (  RUTlo    (1):RUThi    (1), 1)
    double precision, intent(in   ) :: alphaij(  alphaijlo(1):alphaijhi(1), 1*1)
    double precision, intent(in   ) :: alpha  (  alphalo  (1):alphahi  (1), 1)
    double precision, intent(in   ) :: flux_T (  flux_Tlo (1):flux_Thi (1), 1)
    double precision, intent(inout) :: Cs2    (  Cs2lo    (1):Cs2hi    (1))
    double precision, intent(inout) :: CI     (  CIlo     (1):CIhi     (1))
    double precision, intent(inout) :: PrT    (  PrTlo    (1):PrThi    (1))
    double precision, intent(in   ) :: deltax(1)

    double precision :: L(1), M(1), betaij(1), beta(1), T(1), KE(1)
    double precision :: LM
    double precision :: MM
    double precision :: Lkk
    double precision :: bma
    double precision :: TT
    double precision :: KT

    integer :: i
    double precision :: dxinv(1)
    double precision :: deltahat
    double precision, parameter :: small_num = 1d-8
    double precision :: gfaci(lo(1):hi(1)+1)
    integer :: i11
    i11 = (1-1) * 2 + 1 ! index for entry (1,1)

    dxinv = 1.d0/deltax

    gfaci = dxinv(1)
    if (lo(1).le.dmnlo(1) .and. physbc_lo(1).eq.Inflow) gfaci(dmnlo(1)) = gfaci(dmnlo(1)) * TWO
    if (hi(1).gt.dmnhi(1) .and. physbc_hi(1).eq.Inflow) gfaci(dmnhi(1)+1) = gfaci(dmnhi(1)+1) * TWO

    do i=lo(1),hi(1)

       ! Make betaij and beta
       deltahat = fgr*deltax(1)
       call get_sfs_stresses_xdir(i, Q, Qlo, Qhi, betaij(i11), beta(1), T(1), gfaci(i), deltahat)
       T(1) = flux_T(i,1) - T(1)

       ! "resolved turbulent stresses" and others
       M(:) = betaij(:) - alphaij(i,:)
       L(i11) = Kij(i,1) - Q(i,QRHO) * Q(i,QU) * Q(i,QU)
       KE(:) = RUT(i,:) - Q(i,QRHO) * Q(i,QU) * Q(i,QTEMP)

       ! Contractions
       LM = sum(L(:) * M(:)) + small_num
       MM = sum(M(:) * M(:)) + small_num
       Lkk = L(i11) + small_num
       bma = sum(beta(:) - alpha(i,:)) + small_num
       TT = sum(T(:)*T(:)) + small_num
       KT = sum(KE(:)*T(:)) + small_num

       ! Coefficients
       Cs2(i) =  max(LM / MM, small_num)
       CI(i) = max(Lkk / bma , small_num)
       PrT(i) = max(TT / KT , small_num)
    end do

    ! scale Prandtl with Cs2
    PrT(:) = Cs2(:) * PrT(:)

  end subroutine pc_dynamic_smagorinsky_coeffs


  subroutine get_sfs_stresses_xdir(i, Q, Qlo, Qhi, alphaij_xx, alpha, flux_T, gfac, deltabar)

    use meth_params_module, only : QVAR, QRHO, QU, QTEMP
    use amrex_constants_module

    implicit none

    integer, intent(in) :: i
    integer, intent(in) :: Qlo(1), Qhi(1)
    double precision, intent(in   ) :: Q( Qlo(1):Qhi(1), QVAR)
    double precision, intent(inout) :: alphaij_xx
    double precision, intent(inout) :: alpha
    double precision, intent(inout) :: flux_T
    double precision, intent(in   ) :: gfac
    double precision, intent(in   ) :: deltabar

    ! local
    double precision :: Sijmag, Skk, mut
    double precision :: dudx, S
    double precision :: dTdx

    dudx = gfac * (Q(i,QU)    - Q(i-1,QU))
    S = dudx
    Skk = dudx
    Sijmag = sqrt(TWO * S**2)
    mut = Q(i,QRHO) * deltabar**2 * Sijmag

    alphaij_xx = TWO * mut * ( S - THIRD * Skk )
    alpha      = TWO * mut * Sijmag

    dTdx = gfac * (Q(i,QTEMP) - Q(i-1,QTEMP))
    flux_T = mut * dTdx

  end subroutine get_sfs_stresses_xdir

end module lesterm_module
