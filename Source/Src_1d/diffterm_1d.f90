  ! compute the diffusion source terms and fluxes for all the
  ! conservative equations.
  ! 
module diffterm_module

  implicit none

  private

  public :: pc_diffterm

contains

  subroutine pc_diffterm(lo,  hi,&
                         dmnlo, dmnhi,&
                         Q,   Qlo,   Qhi,&
                         Dx,  Dxlo,  Dxhi,&
                         mux, muxlo, muxhi,&
                         xix, xixlo, xixhi,&
                         lamx,lamxlo,lamxhi,&
                         Ax,  Axlo,  Axhi,&
                         fx,  fxlo,  fxhi,&
                         V,   Vlo,   Vhi,&
                         D,   Dlo,   Dhi,&
                         deltax) bind(C, name = "pc_diffterm")

    use network, only : nspecies
    use meth_params_module, only : NVAR, UMX, UMY, UMZ, UEDEN, UFS, QVAR, QU, QV, QPRES, QTEMP, QFS, QRHO
    use amrex_constants_module
    use eos_type_module
    use eos_module
    use prob_params_module, only : physbc_lo, physbc_hi

    implicit none

    integer, intent(in) ::     lo(1),    hi(1)
    integer, intent(in) ::  dmnlo(1), dmnhi(1)
    integer, intent(in) ::    Qlo(1),   Qhi(1)

    integer, intent(in) ::   Dxlo(1),  Dxhi(1)
    integer, intent(in) ::  muxlo(1), muxhi(1)
    integer, intent(in) ::  xixlo(1), xixhi(1)
    integer, intent(in) :: lamxlo(1),lamxhi(1)
    integer, intent(in) ::   Axlo(1),  Axhi(1)
    integer, intent(in) ::   fxlo(1),  fxhi(1)

    integer, intent(in) ::    Dlo(1),   Dhi(1)
    integer, intent(in) ::    Vlo(1),   Vhi(1)

    double precision, intent(in   ) ::    Q(   Qlo(1):   Qhi(1), QVAR)
    double precision, intent(in   ) ::   Dx(  Dxlo(1):  Dxhi(1), nspecies)
    double precision, intent(in   ) ::  mux( muxlo(1): muxhi(1) )
    double precision, intent(in   ) ::  xix( xixlo(1): xixhi(1) )
    double precision, intent(in   ) :: lamx(lamxlo(1):lamxhi(1) )
    double precision, intent(in   ) ::   Ax(  Axlo(1):  Axhi(1) )
    double precision, intent(inout) ::   fx(  fxlo(1):  fxhi(1), NVAR)
    double precision, intent(inout) ::    D(   Dlo(1):   Dhi(1), NVAR)
    double precision, intent(in   ) ::    V(   Vlo(1):   Vhi(1) )
    double precision, intent(in   ) :: deltax(1)

    integer :: i, n
    double precision :: tauxx, divu
    double precision :: Uface, dudx
    double precision :: pface, hface, Xface, Yface
    double precision :: dTdx, dXdx, Vd
    double precision :: Vci(lo(1):hi(1)+1), dlnpi(lo(1):hi(1)+1), gfaci(lo(1):hi(1)+1)
    type(eos_t) :: eosi(lo(1)-1:hi(1)+1)

    double precision, parameter :: twoThirds = 2.d0/3.d0
    double precision :: dxinv(1)

    dxinv = 1.d0/deltax
    do i=lo(1)-1,hi(1)+1
       call build(eosi(i))
    enddo

    gfaci = dxinv(1)
 
    do i=lo(1),hi(1)+1
       dTdx = gfaci(i) * (Q(i,QTEMP) - Q(i-1,QTEMP))
       dudx = gfaci(i) * (Q(i,QU)    - Q(i-1,QU))

       divu = dudx
       tauxx = mux(i)*(2.d0*dudx-twoThirds*divu) + xix(i)*divu
       Uface    = HALF*(Q(i,QU) + Q(i-1,QU))
       pface    = HALF*(Q(i,QPRES) + Q(i-1,QPRES))

       fx(i,UMX)   = - tauxx
       fx(i,UMY)   = 0.d0
       fx(i,UMZ)   = 0.d0
       fx(i,UEDEN) = - tauxx*Uface - lamx(i)*dTdx

       dlnpi(i) = gfaci(i) * (Q(i,QPRES) - Q(i-1,QPRES)) / pface
    end do

    do i=lo(1)-1,hi(1)+1
       eosi(i) % massfrac(:) = Q(i,QFS:QFS+nspecies-1)
       eosi(i) % T           = Q(i,QTEMP)
       eosi(i) % rho         = Q(i,QRHO)
       eosi(i) % p           = Q(i,QPRES)
       call eos_ytx(eosi(i))
       call eos_hi(eosi(i))
    end do

    ! Get species/enthalpy diffusion, compute correction velocity
    Vci = 0.d0
    do n=1,nspecies
       do i = lo(1), hi(1)+1
          Xface = HALF*(eosi(i)%molefrac(n) + eosi(i-1)%molefrac(n))
          Yface = HALF*(eosi(i)%massfrac(n) + eosi(i-1)%massfrac(n))
          hface = HALF*(eosi(i)%hi(n)       + eosi(i-1)%hi(n))

          dXdx = gfaci(i) * (eosi(i)%molefrac(n) - eosi(i-1)%molefrac(n))
          Vd = -Dx(i,n)*(dXdx + (Xface - Yface) * dlnpi(i))

          fx(i,UFS+n-1) = Vd
          Vci(i) = Vci(i) + Vd
          fx(i,UEDEN) = fx(i,UEDEN) + Vd*hface
       end do
    end do

    ! Add correction velocity
    do n=1,nspecies
       do i = lo(1), hi(1)+1
          Yface = HALF*(eosi(i)%massfrac(n) + eosi(i-1)%massfrac(n))
          hface = HALF*(eosi(i)%hi(n)       + eosi(i-1)%hi(n))

          fx(i,UFS+n-1) = fx(i,UFS+n-1) - Yface*Vci(i)
          fx(i,UEDEN)   = fx(i,UEDEN)   - Yface*Vci(i)*hface
       end do
    end do
    
    ! Sscale fluxes by area
    do i=lo(1),hi(1)+1
       fx(i,UMX)   = fx(i,UMX)   * Ax(i)
       fx(i,UEDEN) = fx(i,UEDEN) * Ax(i)
       fx(i,UFS:UFS+nspecies-1) = fx(i,UFS:UFS+nspecies-1) * Ax(i)
    enddo


    do i=lo(1)+1,hi(1)-1
       call destroy(eosi(i))
    enddo

    do n=1,NVAR
       do i = lo(1), hi(1)
          D(i,n) = - (fx(i+1,n)-fx(i,n))/V(i)
       end do
    end do

  end subroutine pc_diffterm

end module diffterm_module
