! compute the LES source terms and fluxes for all the
! conservative equations.
!
! Constant and Dynamic Smagorinsky models are implemented
!
! AT THE MOMENT:
! - no multifluid component of the equations are solved
!
! - Strictly, the formulation assumes perfect gas / constant specific heats / GammaLaw EoS.
!     This assumption (e + p/rho = h = Cp*T) is used to relate the turbulent transport of internal
!     energy + pressure term to the turbulent temperature transport and thus the temperature gradient:
!
!     with <> = filter, {} = density weighted filter, state = [rho, rho*E, rho*u_j]
!     <rho>{u_j h} - <rho>{u_j}{h}  = Smagorinsky Model = <rho> nu_T/Pr_T * d{h}/dx_j
!                                     -> modeled as ->  <rho> nu_T/Pr_T * Cp(<state>) * dT(<state>)/dx_j
!
!     For non-perfect gasses there are two issues with the formulation: h =/= Cp*T, and T is a nonlinear
!     function of the state variables so T(<state>) =/= {T}. For ideal gases where
!     the specific heats are relatively weak functions of temperaure, the error induced from both
!     of these issues is expected to be small relative to the overall error in the Smagorinsky
!     models because the error in the approximation delta_h = Cp*delta_T is ~ dCp/dT*delta_T.
!     But the present formulation is unsuited for real gas equations and superctritical fluids
!     where properties such as Cp are very strong nonlinear functions of the state variables.
!     Converting this term to an enthalpy-based Smagorinsky-like model would alleviate the former issue,
!     but the latter would still remain: h(<state>) =/= {h}.
!
! - The formulation also assumes that the pressure term in the velocity equation can be solved
!     using p(<state>), rather than <p>, which are not equivalent for non-perfect gasses (even for ideal gasses
!     where <p> = <rho>*R*{T} (assuming R is constant) {T} is a non-linear function of the state, so this makes
!     p(<state>) =/= <p>. Again, for weak nonlinearity this is probably not the leading order error in the
!     modeling approach. An alternative explanation is that any discrepancy is part of what is modeled by the
!     Smagorinsky approach.
!
! - The last major modeling assumption is that molecular transport may be closed using, for example,
!     <q> = < -lambda dT/dx > = -lambda(<state>) d{T}/dx
!     This type of closure is standard in the literature, or at least similar closures using <lambda> instead
!     of lambda(<state>) are standard. Although perhaps it has not been verified for real gas EoS. In any event,
!     for coarse LES grids the molecular transport coeffs are much smaller than the turbulent one, so this is
!     standard modeling approach is probably justified. 

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
    gfacj = dxinv(2)
    gfack = dxinv(3)

    ! Fluxes through x-faces
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

             ! SFS heat flux - move state from cell centers to faces to compute cp for flux at face
             sfs_eos_state % massfrac(:) = HALF*( Q(i,j,k,QFS:QFS+nspecies-1) + Q(i-1,j,k,QFS:QFS+nspecies-1) )
             sfs_eos_state % T           = HALF*( Q(i,j,k,QTEMP) + Q(i-1,j,k,QTEMP) )
             call eos_cp(sfs_eos_state)
             fx(i,j,k,UEDEN) = fx(i,j,k,UEDEN) - sfs_eos_state%cp * Cs2 / PrT * flux_T(1)
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

    ! Fluxes through y-faces
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

             ! SFS heat flux - move state to faces for cp calculation
             sfs_eos_state % massfrac(:) = HALF*( Q(i,j,k,QFS:QFS+nspecies-1) + Q(i,j-1,k,QFS:QFS+nspecies-1) )
             sfs_eos_state % T           = HALF*( Q(i,j,k,QTEMP) + Q(i,j-1,k,QTEMP) )
             call eos_cp(sfs_eos_state)
             fy(i,j,k,UEDEN) = fy(i,j,k,UEDEN) - sfs_eos_state%cp * Cs2 / PrT * flux_T(2)
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

    ! Fluxes through z-faces
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
             ! SFS heat flux - move state to faces to calculate cp
             sfs_eos_state % massfrac(:) = HALF* ( Q(i,j,k,QFS:QFS+nspecies-1) + Q(i,j,k-1,QFS:QFS+nspecies-1) )
             sfs_eos_state % T           = HALF* ( Q(i,j,k,QTEMP) + Q(i,j,k-1,QTEMP))
             call eos_cp(sfs_eos_state)
             fz(i,j,k,UEDEN) = fz(i,j,k,UEDEN) - sfs_eos_state%cp * Cs2 / PrT * flux_T(3)
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
                                             alphaijx, alphaijxlo, alphaijxhi,&
                                             alphax, alphaxlo, alphaxhi,&
                                             flux_Tx, flux_Txlo, flux_Txhi,&
                                             Cs2x, Cs2xlo, Cs2xhi, &
                                             CIx, CIxlo, CIxhi, &
                                             Cs2ovPrTx, Cs2ovPrTxlo, Cs2ovPrTxhi, &
                                             Ax,  Axlo,  Axhi,&
                                             fx,  fxlo,  fxhi,&
                                             alphaijy, alphaijylo, alphaijyhi,&
                                             alphay, alphaylo, alphayhi,&
                                             flux_Ty, flux_Tylo, flux_Tyhi,&
                                             Cs2y, Cs2ylo, Cs2yhi, &
                                             CIy, CIylo, CIyhi, &
                                             Cs2ovPrTy, Cs2ovPrTylo, Cs2ovPrTyhi, &
                                             Ay,  Aylo,  Ayhi,&
                                             fy,  fylo,  fyhi,&
                                             alphaijz, alphaijzlo, alphaijzhi,&
                                             alphaz, alphazlo, alphazhi,&
                                             flux_Tz, flux_Tzlo, flux_Tzhi,&
                                             Cs2z, Cs2zlo, Cs2zhi, &
                                             CIz, CIzlo, CIzhi, &
                                             Cs2ovPrTz, Cs2ovPrTzlo, Cs2ovPrTzhi, &
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

    integer, intent(in) ::         lo(3),         hi(3)
    integer, intent(in) ::        Qlo(3),        Qhi(3)
    integer, intent(in) :: alphaijxlo(3), alphaijxhi(3)
    integer, intent(in) ::   alphaxlo(3),   alphaxhi(3)
    integer, intent(in) ::  flux_Txlo(3),  flux_Txhi(3)
    integer, intent(in) ::     Cs2xlo(3),     Cs2xhi(3)
    integer, intent(in) ::      CIxlo(3),      CIxhi(3)
    integer, intent(in) ::Cs2ovPrTxlo(3),Cs2ovPrTxhi(3)
    integer, intent(in) ::       Axlo(3),       Axhi(3)
    integer, intent(in) ::       fxlo(3),       fxhi(3)
    integer, intent(in) :: alphaijylo(3), alphaijyhi(3)
    integer, intent(in) ::   alphaylo(3),   alphayhi(3)
    integer, intent(in) ::  flux_Tylo(3),  flux_Tyhi(3)
    integer, intent(in) ::     Cs2ylo(3),     Cs2yhi(3)
    integer, intent(in) ::      CIylo(3),      CIyhi(3)
    integer, intent(in) ::Cs2ovPrTylo(3),Cs2ovPrTyhi(3)
    integer, intent(in) ::       Aylo(3),       Ayhi(3)
    integer, intent(in) ::       fylo(3),       fyhi(3)
    integer, intent(in) :: alphaijzlo(3), alphaijzhi(3)
    integer, intent(in) ::   alphazlo(3),   alphazhi(3)
    integer, intent(in) ::  flux_Tzlo(3),  flux_Tzhi(3)
    integer, intent(in) ::     Cs2zlo(3),     Cs2zhi(3)
    integer, intent(in) ::      CIzlo(3),      CIzhi(3)
    integer, intent(in) ::Cs2ovPrTzlo(3),Cs2ovPrTzhi(3)
    integer, intent(in) ::       Azlo(3),       Azhi(3)
    integer, intent(in) ::       fzlo(3),       fzhi(3)
    integer, intent(in) ::        Vlo(3),        Vhi(3)
    integer, intent(in) ::        Llo(3),        Lhi(3)

    double precision, intent(in   ) :: Q        (  Qlo       (1):Qhi        (1),  Qlo       (2):Qhi        (2),  Qlo       (3):Qhi        (3), QVAR)
    double precision, intent(in   ) :: alphaijx (  alphaijxlo(1):alphaijxhi (1),  alphaijxlo(2):alphaijxhi (2),  alphaijxlo(3):alphaijxhi (3), 3)
    double precision, intent(in   ) :: alphax   (  alphaxlo  (1):alphaxhi   (1),  alphaxlo  (2):alphaxhi   (2),  alphaxlo  (3):alphaxhi   (3))
    double precision, intent(in   ) :: flux_Tx  (  flux_Txlo (1):flux_Txhi  (1),  flux_Txlo (2):flux_Txhi  (2),  flux_Txlo (3):flux_Txhi  (3))
    double precision, intent(in   ) :: Cs2x     (  Cs2xlo    (1):Cs2xhi     (1),  Cs2xlo    (2):Cs2xhi     (2),  Cs2xlo    (3):Cs2xhi     (3))
    double precision, intent(in   ) :: CIx      (  CIxlo     (1):CIxhi      (1),  CIxlo     (2):CIxhi      (2),  CIxlo     (3):CIxhi      (3))
    double precision, intent(in   ) :: Cs2ovPrTx( Cs2ovPrTxlo(1):Cs2ovPrTxhi(1), Cs2ovPrTxlo(2):Cs2ovPrTxhi(2), Cs2ovPrTxlo(3):Cs2ovPrTxhi(3))
    double precision, intent(in   ) :: Ax       (  Axlo      (1):Axhi       (1),  Axlo      (2):Axhi       (2),  Axlo      (3):Axhi       (3))
    double precision, intent(inout) :: fx       (  fxlo      (1):fxhi       (1),  fxlo      (2):fxhi       (2),  fxlo      (3):fxhi       (3), NVAR)
    double precision, intent(in   ) :: alphaijy (  alphaijylo(1):alphaijyhi (1),  alphaijylo(2):alphaijyhi (2),  alphaijylo(3):alphaijyhi (3), 3)
    double precision, intent(in   ) :: alphay   (  alphaylo  (1):alphayhi   (1),  alphaylo  (2):alphayhi   (2),  alphaylo  (3):alphayhi   (3))
    double precision, intent(in   ) :: flux_Ty  (  flux_Tylo (1):flux_Tyhi  (1),  flux_Tylo (2):flux_Tyhi  (2),  flux_Tylo (3):flux_Tyhi  (3))
    double precision, intent(in   ) :: Cs2y     (  Cs2ylo    (1):Cs2yhi     (1),  Cs2ylo    (2):Cs2yhi     (2),  Cs2ylo    (3):Cs2yhi     (3))
    double precision, intent(in   ) :: CIy      (  CIylo     (1):CIyhi      (1),  CIylo     (2):CIyhi      (2),  CIylo     (3):CIyhi      (3))
    double precision, intent(in   ) :: Cs2ovPrTy( Cs2ovPrTylo(1):Cs2ovPrTyhi(1), Cs2ovPrTylo(2):Cs2ovPrTyhi(2), Cs2ovPrTylo(3):Cs2ovPrTyhi(3))
    double precision, intent(in   ) :: Ay       (  Aylo      (1):Ayhi       (1),  Aylo      (2):Ayhi       (2),  Aylo      (3):Ayhi       (3))
    double precision, intent(inout) :: fy       (  fylo      (1):fyhi       (1),  fylo      (2):fyhi       (2),  fylo      (3):fyhi       (3), NVAR)
    double precision, intent(in   ) :: alphaijz (  alphaijzlo(1):alphaijzhi (1),  alphaijzlo(2):alphaijzhi (2),  alphaijzlo(3):alphaijzhi (3), 3)
    double precision, intent(in   ) :: alphaz   (  alphazlo  (1):alphazhi   (1),  alphazlo  (2):alphazhi   (2),  alphazlo  (3):alphazhi   (3))
    double precision, intent(in   ) :: flux_Tz  (  flux_Tzlo (1):flux_Tzhi  (1),  flux_Tzlo (2):flux_Tzhi  (2),  flux_Tzlo (3):flux_Tzhi  (3))
    double precision, intent(in   ) :: Cs2z     (  Cs2zlo    (1):Cs2zhi     (1),  Cs2zlo    (2):Cs2zhi     (2),  Cs2zlo    (3):Cs2zhi     (3))
    double precision, intent(in   ) :: CIz      (  CIzlo     (1):CIzhi      (1),  CIzlo     (2):CIzhi      (2),  CIzlo     (3):CIzhi      (3))
    double precision, intent(in   ) :: Cs2ovPrTz( Cs2ovPrTzlo(1):Cs2ovPrTzhi(1), Cs2ovPrTzlo(2):Cs2ovPrTzhi(2), Cs2ovPrTzlo(3):Cs2ovPrTzhi(3))
    double precision, intent(in   ) :: Az       (  Azlo      (1):Azhi       (1),  Azlo      (2):Azhi       (2),  Azlo      (3):Azhi       (3))
    double precision, intent(inout) :: fz       (  fzlo      (1):fzhi       (1),  fzlo      (2):fzhi       (2),  fzlo      (3):fzhi       (3), NVAR)
    double precision, intent(inout) :: L        (  Llo       (1):Lhi        (1),  Llo       (2):Lhi        (2),  Llo       (3):Lhi        (3), NVAR)
    double precision, intent(in   ) :: V        (  Vlo       (1):Vhi        (1),  Vlo       (2):Vhi        (2),  Vlo       (3):Vhi        (3))
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

    ! Fluxes through x-faces
    do k= lo(3), hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)+1
             ! SFS stress
             sigmaxx = Cs2x(i,j,k) * alphaijx(i,j,k,1) - THIRD * CIx(i,j,k) * alphax(i,j,k)
             sigmaxy = Cs2x(i,j,k) * alphaijx(i,j,k,2)
             sigmaxz = Cs2x(i,j,k) * alphaijx(i,j,k,3)
             fx(i,j,k,UMX)   = - sigmaxx
             fx(i,j,k,UMY)   = - sigmaxy
             fx(i,j,k,UMZ)   = - sigmaxz

             ! SFS turbulent diffusion
             Uface(:) = HALF*(Q(i,j,k,QU:QW) + Q(i-1,j,k,QU:QW))
             fx(i,j,k,UEDEN) = - sigmaxx*Uface(1) - sigmaxy*Uface(2) - sigmaxz*Uface(3)

             ! SFS heat flux - move state from cell centers to faces to compute cp for flux at face
             sfs_eos_state % massfrac(:) = HALF*( Q(i,j,k,QFS:QFS+nspecies-1) + Q(i-1,j,k,QFS:QFS+nspecies-1) )
             sfs_eos_state % T           = HALF*( Q(i,j,k,QTEMP) + Q(i-1,j,k,QTEMP) )
             call eos_cp(sfs_eos_state)
             fx(i,j,k,UEDEN) = fx(i,j,k,UEDEN) - sfs_eos_state%cp * Cs2ovPrTx(i,j,k) * flux_Tx(i,j,k)
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

    ! Fluxes through y-faces
    do k=lo(3),hi(3)
       do i=lo(1),hi(1)
          do j=lo(2),hi(2)+1
             ! SFS stress
             sigmayx = Cs2y(i,j,k) * alphaijy(i,j,k,1)
             sigmayy = Cs2y(i,j,k) * alphaijy(i,j,k,2) - THIRD * CIy(i,j,k) * alphay(i,j,k)
             sigmayz = Cs2y(i,j,k) * alphaijy(i,j,k,3)
             fy(i,j,k,UMX)   = - sigmayx
             fy(i,j,k,UMY)   = - sigmayy
             fy(i,j,k,UMZ)   = - sigmayz

             ! SFS turbulent diffusion
             Uface(:) = HALF*(Q(i,j,k,QU:QW) + Q(i,j-1,k,QU:QW))
             fy(i,j,k,UEDEN) = - sigmayx*Uface(1) - sigmayy*Uface(2) - sigmayz*Uface(3)

             ! SFS heat flux - move state to faces for cp calculation
             sfs_eos_state % massfrac(:) = HALF*( Q(i,j,k,QFS:QFS+nspecies-1) + Q(i,j-1,k,QFS:QFS+nspecies-1) )
             sfs_eos_state % T           = HALF*( Q(i,j,k,QTEMP) + Q(i,j-1,k,QTEMP) )
             call eos_cp(sfs_eos_state)
             fy(i,j,k,UEDEN) = fy(i,j,k,UEDEN) -  sfs_eos_state%cp * Cs2ovPrTy(i,j,k) * flux_Ty(i,j,k)
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

    ! Fluxes through z-faces
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          do k=lo(3),hi(3)+1
             ! SFS stress
             sigmazx = Cs2z(i,j,k) * alphaijz(i,j,k,1)
             sigmazy = Cs2z(i,j,k) * alphaijz(i,j,k,2)
             sigmazz = Cs2z(i,j,k) * alphaijz(i,j,k,3) - THIRD * CIz(i,j,k) * alphaz(i,j,k)
             fz(i,j,k,UMX)   = - sigmazx
             fz(i,j,k,UMY)   = - sigmazy
             fz(i,j,k,UMZ)   = - sigmazz

             ! SFS turbulent diffusion
             Uface(:) = HALF*(Q(i,j,k,QU:QW) + Q(i,j,k-1,QU:QW))
             fz(i,j,k,UEDEN) = - sigmazx*Uface(1) - sigmazy*Uface(2) - sigmazz*Uface(3)

             ! SFS heat flux - move state to faces to calculate cp
             sfs_eos_state % massfrac(:) = HALF* ( Q(i,j,k,QFS:QFS+nspecies-1) + Q(i,j,k-1,QFS:QFS+nspecies-1) )
             sfs_eos_state % T           = HALF* ( Q(i,j,k,QTEMP) + Q(i,j,k-1,QTEMP))
             call eos_cp(sfs_eos_state)
             fz(i,j,k,UEDEN) = fz(i,j,k,UEDEN) - sfs_eos_state%cp * Cs2ovPrTz(i,j,k) * flux_Tz(i,j,k)
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

  ! Quantities
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

    use meth_params_module, only : QVAR, QRHO, QU, QV, QW, QTEMP
    use amrex_constants_module
    use prob_params_module, only : physbc_lo, physbc_hi, Inflow

    implicit none

    integer, intent(in) :: lo       (3), hi       (3)
    integer, intent(in) :: dmnlo    (3), dmnhi    (3)
    integer, intent(in) :: Qlo      (3), Qhi      (3)
    integer, intent(in) :: Kijlo    (3), Kijhi    (3)
    integer, intent(in) :: RUTlo    (3), RUThi    (3)
    integer, intent(in) :: alphaijlo(3), alphaijhi(3)
    integer, intent(in) :: alphalo  (3), alphahi  (3)
    integer, intent(in) :: flux_Tlo (3), flux_Thi (3)
    integer, intent(in) :: fgr

    double precision, intent(in   ) :: Q      (  Qlo      (1):Qhi      (1),  Qlo     (2): Qhi     (2),  Qlo       (3):Qhi      (3), QVAR)
    double precision, intent(inout) :: Kij    (  Kijlo    (1):Kijhi    (1),  Kijlo    (2):Kijhi    (2),  Kijlo    (3):Kijhi    (3), 6)
    double precision, intent(inout) :: RUT    (  RUTlo    (1):RUThi    (1),  RUTlo    (2):RUThi    (2),  RUTlo    (3):RUThi    (3), 3)
    double precision, intent(inout) :: alphaij(  alphaijlo(1):alphaijhi(1),  alphaijlo(2):alphaijhi(2),  alphaijlo(3):alphaijhi(3), 3*3)
    double precision, intent(inout) :: alpha  (  alphalo  (1):alphahi  (1),  alphalo  (2):alphahi  (2),  alphalo  (3):alphahi  (3))
    double precision, intent(inout) :: flux_T (  flux_Tlo (1):flux_Thi (1),  flux_Tlo (2):flux_Thi (2),  flux_Tlo (3):flux_Thi (3), 3)
    double precision, intent(in   ) :: deltax(3)

    integer :: i, j, k
    double precision :: dxinv(3)
    double precision :: deltabar
    double precision :: gfaci(lo(1):hi(1)+1)
    double precision :: gfacj(lo(2):hi(2)+1)
    double precision :: gfack(lo(3):hi(3)+1)
    double precision :: gfacijk(3)

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
    gfacj = dxinv(2)
    gfack = dxinv(3)

    ! Maybe revisit how to determine deltahat in the case where dx1 != dx2 != dx3 (irrelevant for now as this is presently not allowed)
    deltabar = THIRD*fgr*(deltax(1)+deltax(2)+deltax(3))
    
    do k=lo(3),hi(3)
       gfacijk(3) = gfack(k)
       do j=lo(2),hi(2)
          gfacijk(2) = gfacj(j)
          do i=lo(1),hi(1)
             gfacijk(1) = gfaci(i)
             call get_sfs_stresses_cc(i, j, k, Q, Qlo, Qhi, &
                                      alphaij(i,j,k,i11), alphaij(i,j,k,i12), alphaij(i,j,k,i13), alphaij(i,j,k,i22), alphaij(i,j,k,i23), alphaij(i,j,k,i33), &
                                      alpha(i,j,k), flux_T(i,j,k,:), gfacijk, deltabar)
             alphaij(i,j,k,i21) = alphaij(i,j,k,i12)
             alphaij(i,j,k,i31) = alphaij(i,j,k,i13)
             alphaij(i,j,k,i32) = alphaij(i,j,k,i23)

             Kij(i,j,k, 1) = Q(i,j,k,QRHO) * Q(i,j,k,QU) * Q(i,j,k,QU)
             Kij(i,j,k, 2) = Q(i,j,k,QRHO) * Q(i,j,k,QU) * Q(i,j,k,QV)
             Kij(i,j,k, 3) = Q(i,j,k,QRHO) * Q(i,j,k,QU) * Q(i,j,k,QW)
             Kij(i,j,k, 4) = Q(i,j,k,QRHO) * Q(i,j,k,QV) * Q(i,j,k,QV)
             Kij(i,j,k, 5) = Q(i,j,k,QRHO) * Q(i,j,k,QV) * Q(i,j,k,QW)
             Kij(i,j,k, 6) = Q(i,j,k,QRHO) * Q(i,j,k,QW) * Q(i,j,k,QW)

             RUT(i,j,k,:) = Q(i,j,k,QRHO) * Q(i,j,k,QU:QW) * Q(i,j,k,QTEMP)
             
          end do
       end do
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
                                           Cs2ovPrT, Cs2ovPrTlo, Cs2ovPrThi, &
                                           fgr,&
                                           deltax) bind(C, name = "pc_dynamic_smagorinsky_coeffs")

    use meth_params_module, only : QVAR, QRHO, QU, QV, QW, QTEMP
    use amrex_constants_module
    use prob_params_module, only : physbc_lo, physbc_hi, Inflow

    implicit none

    integer, intent(in) :: lo        (3), hi        (3)
    integer, intent(in) :: dmnlo     (3), dmnhi     (3)
    integer, intent(in) :: Qlo       (3), Qhi       (3)
    integer, intent(in) :: Kijlo     (3), Kijhi     (3)
    integer, intent(in) :: RUTlo     (3), RUThi     (3)
    integer, intent(in) :: alphaijlo (3), alphaijhi (3)
    integer, intent(in) :: alphalo   (3), alphahi   (3)
    integer, intent(in) :: flux_Tlo  (3), flux_Thi  (3)
    integer, intent(in) :: Cs2lo     (3), Cs2hi     (3)
    integer, intent(in) :: CIlo      (3), CIhi      (3)
    integer, intent(in) :: Cs2ovPrTlo(3), Cs2ovPrThi(3)
    integer, intent(in) :: fgr

    double precision, intent(in   ) :: Q       (  Qlo      (1):Qhi       (1),   Qlo     (2): Qhi      (2),   Qlo     (3): Qhi      (3), QVAR)
    double precision, intent(in   ) :: Kij     (  Kijlo    (1):Kijhi     (1),  Kijlo    (2):Kijhi     (2),  Kijlo    (3):Kijhi     (3), 6)
    double precision, intent(in   ) :: RUT     (  RUTlo    (1):RUThi     (1),  RUTlo    (2):RUThi     (2),  RUTlo    (3):RUThi     (3), 3)
    double precision, intent(in   ) :: alphaij (  alphaijlo(1):alphaijhi (1),  alphaijlo(2):alphaijhi (2),  alphaijlo(3):alphaijhi (3), 3*3)
    double precision, intent(in   ) :: alpha   (  alphalo  (1):alphahi   (1),  alphalo  (2):alphahi   (2),  alphalo  (3):alphahi   (3))
    double precision, intent(in   ) :: flux_T  (  flux_Tlo (1):flux_Thi  (1),  flux_Tlo (2):flux_Thi  (2),  flux_Tlo (3):flux_Thi  (3), 3)
    double precision, intent(inout) :: Cs2     (  Cs2lo    (1):Cs2hi     (1),  Cs2lo    (2):Cs2hi     (2),  Cs2lo    (3):Cs2hi     (3))
    double precision, intent(inout) :: CI      (  CIlo     (1):CIhi      (1),  CIlo     (2):CIhi      (2),  CIlo     (3):CIhi      (3))
    double precision, intent(inout) :: Cs2ovPrT( Cs2ovPrTlo(1):Cs2ovPrThi(1), Cs2ovPrTlo(2):Cs2ovPrThi(2), Cs2ovPrTlo(3):Cs2ovPrThi(3))
    double precision, intent(in   ) :: deltax(3)

    double precision :: L(3*3), M(3*3), betaij(3*3), beta, T(3), KE(3)
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
    double precision :: gfacijk(3)
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
    gfacj = dxinv(2)
    gfack = dxinv(3)

    ! Calculate terms based on fluxes (located at edges), and move the values to cell centers
    do k=lo(3),hi(3)
       gfacijk(3) = gfack(k)
       do j=lo(2),hi(2)
          gfacijk(2) = gfacj(j)
          do i=lo(1),hi(1)
             ! Maybe revisit how to determine deltahat in the case where dx1 != dx2 != dx3 (irrelevant for now as this is presently not allowed)
             deltahat = THIRD*fgr*(deltax(1)+deltax(2)+deltax(3))
             gfacijk(1) = gfaci(i)
             call get_sfs_stresses_cc(i, j, k, Q, Qlo, Qhi, betaij(i11), betaij(i12), betaij(i13), betaij(i22), betaij(i23), betaij(i33), beta, T(:), gfacijk, deltahat)
             betaij(i21) = betaij(i12)
             betaij(i31) = betaij(i13)
             betaij(i32) = betaij(i23)
             
             T(:) = flux_T(i,j,k,:) - T(:)
             
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
             KE(:)  = RUT(i,j,k,:) - Q(i,j,k,QRHO) * Q(i,j,k,QU:QW) * Q(i,j,k,QTEMP)

             ! Contractions
             LM = sum(L(:) * M(:))
             MM = sum(M(:) * M(:))
             Lkk = (L(i11) + L(i22) + L(i33))
             bma = (beta - alpha(i,j,k))
             TT = sum(T(:)*T(:))
             KT = sum(KE(:)*T(:))

             ! Coefficients
             Cs2(i,j,k)      = max(LM  / (MM  + small_num), small_num)
             CI (i,j,k)      = max(Lkk / (bma + small_num), small_num)
             Cs2ovPrT(i,j,k) = max(KT  / (TT  + small_num), small_num)

             ! Limit CI to reasonable bounds
             CI(i,j,k) = min(CI(i,j,k), 1.0d+0)
             
          end do
       end do
    end do

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
    mut = HALF*(Q(i,j,k,QRHO) + Q(i-1,j,k,QRHO)) * deltabar**2 * Sijmag

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
    mut = HALF*(Q(i,j,k,QRHO) + Q(i,j-1,k,QRHO)) * deltabar**2 * Sijmag

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
    mut = HALF*(Q(i,j,k,QRHO) + Q(i,j,k-1,QRHO)) * deltabar**2 * Sijmag

    alphaij_zx = TWO * mut * S(3,1)
    alphaij_zy = TWO * mut * S(3,2)
    alphaij_zz = TWO * mut * ( S(3,3) - THIRD * Skk)
    alpha      = TWO * mut * Sijmag

    dTdz = gfac * (Q(i,j,k,QTEMP) - Q(i,j,k-1,QTEMP))
    flux_T = mut * dTdz

  end subroutine get_sfs_stresses_zdir

  ! Above subroutines get edge values (required for calculating actual momentum fluxes)
  ! Below subroutine gets node values (required for dynamic procedure to determine coeffs)
  
  subroutine get_sfs_stresses_cc(i, j, k, Q, Qlo, Qhi, alpha11, alpha12, alpha13, alpha22, alpha23, alpha33, alpha, flux_T, gfac, deltabar)
    
    use meth_params_module, only : QVAR, QRHO, QU, QV, QW, QTEMP
    use amrex_constants_module

    implicit none

    integer, intent(in) :: i, j, k
    integer, intent(in) :: Qlo(3), Qhi(3)
    double precision, intent(in   ) :: Q( Qlo(1):Qhi(1),  Qlo(2): Qhi(2),  Qlo(3):Qhi(3), QVAR)
    double precision, intent(inout) :: alpha11, alpha12, alpha13, alpha22, alpha23, alpha33
    double precision, intent(inout) :: alpha
    double precision, intent(inout) :: flux_T(3)
    double precision, intent(in   ) :: gfac(3)
    double precision, intent(in   ) :: deltabar

    ! local
    double precision :: Sijmag, Skk, mut
    double precision :: dUdx(3,3), S(3,3)
    double precision :: dTdx(3)

    ! Calculate derivatives at cell centers, second order central difference
    dUdx(1,1) = gfac(1) * (Q(i+1,j,k,QU)    - Q(i-1,j,k,QU))
    dUdx(1,2) = gfac(2) * (Q(i,j+1,k,QU)    - Q(i,j-1,k,QU))
    dUdx(1,3) = gfac(3) * (Q(i,j,k+1,QU)    - Q(i,j,k-1,QU))
    dUdx(2,1) = gfac(1) * (Q(i+1,j,k,QV)    - Q(i-1,j,k,QV))
    dUdx(2,2) = gfac(2) * (Q(i,j+1,k,QV)    - Q(i,j-1,k,QV))
    dUdx(2,3) = gfac(3) * (Q(i,j,k+1,QV)    - Q(i,j,k-1,QV))
    dUdx(3,1) = gfac(1) * (Q(i+1,j,k,QW)    - Q(i-1,j,k,QW))
    dUdx(3,2) = gfac(2) * (Q(i,j+1,k,QW)    - Q(i,j-1,k,QW))
    dUdx(3,3) = gfac(3) * (Q(i,j,k+1,QW)    - Q(i,j,k-1,QW))
    dUdx = HALF * dUdx

    S(:,:) = HALF * (dUdx(:,:) + transpose(dUdx(:,:)))
    Skk = S(1,1) + S(2,2) + S(3,3)
    Sijmag = sqrt(TWO * sum(S(:,:)**2))
    mut = Q(i,j,k,QRHO) * deltabar**2 * Sijmag

    alpha11    = TWO * mut * (S(1,1) - THIRD * Skk)
    alpha12    = TWO * mut * (S(1,2))
    alpha13    = TWO * mut * (S(1,3))
    alpha22    = TWO * mut * (S(2,2) - THIRD * Skk)
    alpha23    = TWO * mut * (S(2,3))
    alpha33    = TWO * mut * (S(3,3) - THIRD * Skk)
    alpha      = TWO * mut * Sijmag
    
    dTdx(1) = gfac(1) * (Q(i+1,j,k,QTEMP) - Q(i-1,j,k,QTEMP))
    dTdx(2) = gfac(2) * (Q(i,j+1,k,QTEMP) - Q(i,j-1,k,QTEMP))
    dTdx(3) = gfac(3) * (Q(i,j,k+1,QTEMP) - Q(i,j,k-1,QTEMP))
    dTdx = HALF * dTdx
    flux_T = mut * dTdx
    
  end subroutine get_sfs_stresses_cc
  
end module lesterm_module
