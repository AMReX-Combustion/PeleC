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
    use meth_params_module, only : NVAR, UMX, UEDEN, UFS, QVAR, QU, QPRES, QTEMP, QFS, QRHO
    use amrex_constants_module
    use eos_type_module
    use eos_module

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

    integer, intent(in) ::    Vlo(1),   Vhi(1)
    integer, intent(in) ::    Dlo(1),   Dhi(1)

    double precision, intent(in   ) ::    Q(   Qlo(1):   Qhi(1), QVAR)
    double precision, intent(in   ) ::   Dx(  Dxlo(1):  Dxhi(1), nspecies)
    double precision, intent(in   ) ::  mux( muxlo(1): muxhi(1)  )
    double precision, intent(in   ) ::  xix( xixlo(1): xixhi(1)  )
    double precision, intent(in   ) :: lamx(lamxlo(1):lamxhi(1)  )
    double precision, intent(in   ) ::   Ax(  Axlo(1):  Axhi(1)  )
    double precision, intent(inout) ::   fx(  fxlo(1):  fxhi(1), NVAR)
    double precision, intent(inout) ::    D(   Dlo(1):   Dhi(1), NVAR)
    double precision, intent(in   ) ::    V(   Vlo(1):   Vhi(1)  )
    double precision, intent(in   ) :: deltax(1)

    integer :: i, n, nn
    double precision :: tauxx, dudx, divu
    double precision :: uface, pface, hface, Yface
    double precision :: dTdx, Vd
    double precision :: Vc(lo(1):hi(1)+1)
    double precision :: ddrive(lo(1):hi(1)+1,nspecies), gradY(lo(1):hi(1)+1,nspecies)
    double precision :: gradP(lo(1):hi(1)+1), dsum(lo(1):hi(1)+1)
    double precision, parameter :: twoThirds = 2.d0/3.d0
    double precision :: dxinv(1)
    type(eos_t) :: eos_state(lo(1)-1:hi(1)+1)

    dxinv = 1.d0/deltax
    do i=lo(1)-1,hi(1)+1
       call build(eos_state(i))
    enddo

    do i=lo(1),hi(1)+1
       ! viscous stress
       dudx = dxinv(1)*(Q(i,QU) - Q(i-1,QU))
       divu = dudx
       tauxx = mux(i)*(2.d0*dudx-twoThirds*divu) + xix(i)*divu
       uface = HALF*(Q(i,QU)    + Q(i-1,QU))
       pface = HALF*(Q(i,QPRES) + Q(i-1,QPRES))

       fx(i,UMX)   = -tauxx
       fx(i,UEDEN) = -tauxx*uface

       ! thermal conduction
       dTdx = dxinv(1) * (Q(i,QTEMP) - Q(i-1,QTEMP))
       fx(i,UEDEN) = fx(i,UEDEN) - lamx(i)*dTdx

       ! (1/p)(dp/dx)
   !    dlnp(i) = dxinv(1) * (Q(i,QPRES) - Q(i-1,QPRES)) / pface
       gradP(i) = dxinv(1) * (Q(i,QPRES) - Q(i-1,QPRES))
       Vc(i) = 0.d0
    end do

    do i=lo(1)-1,hi(1)+1
       eos_state(i) % massfrac(:) = Q(i,QFS:QFS+nspecies-1)
       eos_state(i) % T           = Q(i,QTEMP)
       eos_state(i) % rho         = Q(i,QRHO)
       call eos_ytx(eos_state(i))
       call eos_hi(eos_state(i))
       call eos_get_transport(eos_state(i))
    end do

      do n=1,nspecies
          do i = lo(1), hi(1)+1

             gradY(i,n) = dxinv(1) * (eos_state(i)%massfrac(n) - eos_state(i-1)%massfrac(n))

!    put in P term
             ddrive(i,n) = 0.5d0*(eos_state(i)% diP(n) + eos_state(i-1)% diP(n)) * gradP(i)

          enddo
       enddo

       do n=1,nspecies
       do nn=1,nspecies
          do i = lo(1), hi(1)+1

            ddrive(i,n) = ddrive(i,n)+ 0.5d0* (eos_state(i) % dijY(n,nn) &
                           + eos_state(i-1) % dijY(n,nn)) * gradY(i,nn)

          enddo
       enddo
       enddo

       dsum = 0.d0
       do n=1,nspecies
          do i = lo(1), hi(1)+1

            dsum(i) = dsum(i) + ddrive(i,n)

          enddo
       enddo
       do n=1,nspecies
          do i = lo(1), hi(1)+1

            ddrive(i,n) =  ddrive(i,n) - eos_state(i)%massfrac(n) * dsum(i)

          enddo
       enddo


    ! Get species/enthalpy diffusion, compute correction velocity
    do n=1,nspecies
       do i = lo(1), hi(1)+1
          hface = HALF*(eos_state(i)%hi(n)       + eos_state(i-1)%hi(n))

          Vd = -Dx(i,n)*ddrive(i,n)

          fx(i,UFS+n-1) = Vd
          Vc(i) = Vc(i) + Vd
          fx(i,UEDEN) = fx(i,UEDEN) + Vd*hface
       end do
    end do

    ! Add correction velocity
    do n=1,nspecies
       do i = lo(1), hi(1)+1
          Yface = HALF*(eos_state(i)%massfrac(n) + eos_state(i-1)%massfrac(n))
          hface = HALF*(eos_state(i)%hi(n)       + eos_state(i-1)%hi(n))

          fx(i,UFS+n-1) = fx(i,UFS+n-1) - Yface*Vc(i)
          fx(i,UEDEN)   = fx(i,UEDEN)   - Yface*Vc(i)*hface
       end do
    end do

    ! Scale fluxes by area
    do i=lo(1),hi(1)+1
       fx(i,UMX)   = fx(i,UMX)   * Ax(i)
       fx(i,UEDEN) = fx(i,UEDEN) * Ax(i)
       fx(i,UFS:UFS+nspecies-1) = fx(i,UFS:UFS+nspecies-1) * Ax(i)
    enddo

    do i=lo(1)+1,hi(1)-1
       call destroy(eos_state(i))
    enddo

    do n=1,NVAR
       do i = lo(1), hi(1)
          D(i,n) = - (fx(i+1,n)-fx(i,n)) / V(i)
       end do
    end do

  end subroutine pc_diffterm

end module diffterm_module
