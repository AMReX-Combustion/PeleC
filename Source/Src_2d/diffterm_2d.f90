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
                         tx,  txlo,  txhi,&
                         Ax,  Axlo,  Axhi,&
                         fx,  fxlo,  fxhi,&
                         Dy,  Dylo,  Dyhi,&  
                         muy, muylo, muyhi,& 
                         xiy, xiylo, xiyhi,& 
                         lamy,lamylo,lamyhi,&
                         ty,  tylo,  tyhi,&
                         Ay,  Aylo,  Ayhi,&
                         fy,  fylo,  fyhi,&
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

    integer, intent(in) ::     lo(2),    hi(2)
    integer, intent(in) ::  dmnlo(2), dmnhi(2)
    integer, intent(in) ::    Qlo(2),   Qhi(2)

    integer, intent(in) ::   Dxlo(2),  Dxhi(2)
    integer, intent(in) ::  muxlo(2), muxhi(2)
    integer, intent(in) ::  xixlo(2), xixhi(2)
    integer, intent(in) :: lamxlo(2),lamxhi(2)
    integer, intent(in) ::   txlo(2),  txhi(2)
    integer, intent(in) ::   Axlo(2),  Axhi(2)
    integer, intent(in) ::   fxlo(2),  fxhi(2)

    integer, intent(in) ::   Dylo(2),  Dyhi(2)
    integer, intent(in) ::  muylo(2), muyhi(2)
    integer, intent(in) ::  xiylo(2), xiyhi(2)
    integer, intent(in) :: lamylo(2),lamyhi(2)
    integer, intent(in) ::   tylo(2),  tyhi(2)
    integer, intent(in) ::   Aylo(2),  Ayhi(2)
    integer, intent(in) ::   fylo(2),  fyhi(2)

    integer, intent(in) ::    Dlo(2),   Dhi(2)
    integer, intent(in) ::    Vlo(2),   Vhi(2)
    
    double precision, intent(in   ) ::    Q(   Qlo(1):   Qhi(1),   Qlo(2):   Qhi(2), QVAR)
    double precision, intent(in   ) ::   Dx(  Dxlo(1):  Dxhi(1),  Dxlo(2):  Dxhi(2), nspecies)
    double precision, intent(in   ) ::  mux( muxlo(1): muxhi(1), muxlo(2): muxhi(2) )
    double precision, intent(in   ) ::  xix( xixlo(1): xixhi(1), xixlo(2): xixhi(2) )
    double precision, intent(in   ) :: lamx(lamxlo(1):lamxhi(1),lamxlo(2):lamxhi(2) )
    double precision, intent(in   ) ::   tx(  txlo(1):  txhi(1),  txlo(2):  txhi(2), 2)
    double precision, intent(in   ) ::   Ax(  Axlo(1):  Axhi(1),  Axlo(2):  Axhi(2) )
    double precision, intent(inout) ::   fx(  fxlo(1):  fxhi(1),  fxlo(2):  fxhi(2), NVAR)
    double precision, intent(in   ) ::   Dy(  Dylo(1):  Dyhi(1),  Dylo(2):  Dyhi(2), nspecies)
    double precision, intent(in   ) ::  muy( muylo(1): muyhi(1), muylo(2): muyhi(2) )
    double precision, intent(in   ) ::  xiy( xiylo(1): xiyhi(1), xiylo(2): xiyhi(2) )
    double precision, intent(in   ) :: lamy(lamylo(1):lamyhi(1),lamylo(2):lamyhi(2) )
    double precision, intent(in   ) ::   ty(  tylo(1):  tyhi(1),  tylo(2):  tyhi(2), 2)
    double precision, intent(in   ) ::   Ay(  Aylo(1):  Ayhi(1),  Aylo(2):  Ayhi(2) )
    double precision, intent(inout) ::   fy(  fylo(1):  fyhi(1),  fylo(2):  fyhi(2), NVAR)
    double precision, intent(inout) ::    D(   Dlo(1):   Dhi(1),   Dlo(2):   Dhi(2), NVAR)
    double precision, intent(in   ) ::    V(   Vlo(1):   Vhi(1),   Vlo(2):   Vhi(2) )
    double precision, intent(in   ) :: deltax(2)

    integer :: i, j, k, n
    double precision :: tauxx, tauxy, tauyx, tauyy, divu
    double precision :: Uface(2), dudx,dvdx,dudy,dvdy
    double precision :: pface, hface, Xface, Yface
    double precision :: dTdx, dTdy, dXdx, dXdy, Vd
    double precision :: Vci(lo(1):hi(1)+1), dlnpi(lo(1):hi(1)+1), gfaci(lo(1):hi(1)+1)
    double precision :: Vcj(lo(2):hi(2)+1), dlnpj(lo(2):hi(2)+1), gfacj(lo(2):hi(2)+1)
    type(eos_t) :: eosi(lo(1)-1:hi(1)+1)
    type(eos_t) :: eosj(lo(2)-1:hi(2)+1)

    double precision, parameter :: twoThirds = 2.d0/3.d0
    double precision :: dxinv(2)

    dxinv = 1.d0/deltax
    do i=lo(1)-1,hi(1)+1
       call build(eosi(i))
    enddo

    do j=lo(2),hi(2)

       gfaci = dxinv(1)
       
       do i=lo(1),hi(1)+1
          dTdx = gfaci(i) * (Q(i,j,QTEMP) - Q(i-1,j,QTEMP))
          dudx = gfaci(i) * (Q(i,j,QU)    - Q(i-1,j,QU))
          dvdx = gfaci(i) * (Q(i,j,QV)    - Q(i-1,j,QV))
          dudy = tx(i,j,1)
          dvdy = tx(i,j,2)

          divu = dudx + dvdy
          tauxx = mux(i,j)*(2.d0*dudx-twoThirds*divu) + xix(i,j)*divu
          tauxy = mux(i,j)*(dudy+dvdx)
          Uface(1) = HALF*(Q(i,j,QU) + Q(i-1,j,QU))
          Uface(2) = HALF*(Q(i,j,QV) + Q(i-1,j,QV))
          pface    = HALF*(Q(i,j,QPRES) + Q(i-1,j,QPRES))

          fx(i,j,UMX)   = - tauxx
          fx(i,j,UMY)   = - tauxy
          fx(i,j,UMZ)   = 0.d0
          fx(i,j,UEDEN) = - tauxx*Uface(1) - tauxy*Uface(2) - lamx(i,j)*dTdx

          dlnpi(i) = gfaci(i) * (Q(i,j,QPRES) - Q(i-1,j,QPRES)) / pface
       end do

       do i=lo(1)-1,hi(1)+1
          eosi(i) % massfrac(:) = Q(i,j,QFS:QFS+nspecies-1)
          eosi(i) % T           = Q(i,j,QTEMP)
          eosi(i) % rho         = Q(i,j,QRHO)
          eosi(i) % p           = Q(i,j,QPRES)
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
             Vd = -Dx(i,j,n)*(dXdx + (Xface - Yface) * dlnpi(i))

             fx(i,j,UFS+n-1) = Vd
             Vci(i) = Vci(i) + Vd
             fx(i,j,UEDEN) = fx(i,j,UEDEN) + Vd*hface
          end do
       end do

       ! Add correction velocity
       do n=1,nspecies
          do i = lo(1), hi(1)+1
             Yface = HALF*(eosi(i)%massfrac(n) + eosi(i-1)%massfrac(n))
             hface = HALF*(eosi(i)%hi(n)       + eosi(i-1)%hi(n))

             fx(i,j,UFS+n-1) = fx(i,j,UFS+n-1) - Yface*Vci(i)
             fx(i,j,UEDEN)   = fx(i,j,UEDEN)   - Yface*Vci(i)*hface
          end do
       end do
    end do
    
    ! Sscale fluxes by area
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)+1
          fx(i,j,UMX)   = fx(i,j,UMX)   * Ax(i,j)
          fx(i,j,UMY)   = fx(i,j,UMY)   * Ax(i,j)
          fx(i,j,UEDEN) = fx(i,j,UEDEN) * Ax(i,j)
          fx(i,j,UFS:UFS+nspecies-1) = fx(i,j,UFS:UFS+nspecies-1) * Ax(i,j)
       enddo
    enddo


    do i=lo(1)+1,hi(1)-1
       call destroy(eosi(i))
    enddo
    do j=lo(2)-1,hi(2)+1
       call build(eosj(j))
    enddo

    do i=lo(1),hi(1)

       gfacj = dxinv(2)
       
       do j=lo(2),hi(2)+1
              
          dTdy = gfacj(j) * (Q(i,j,QTEMP) - Q(i,j-1,QTEMP))
          dudy = gfacj(j) * (Q(i,j,QU)    - Q(i,j-1,QU))
          dvdy = gfacj(j) * (Q(i,j,QV)    - Q(i,j-1,QV))
          dudx = ty(i,j,1)
          dvdx = ty(i,j,2)

          divu = dudx + dvdy
          tauyx = muy(i,j)*(dudy+dvdx)
          tauyy = muy(i,j)*(2.d0*dvdy-twoThirds*divu) + xiy(i,j)*divu
          Uface(1) = HALF*(Q(i,j,QU) + Q(i,j-1,QU))
          Uface(2) = HALF*(Q(i,j,QV) + Q(i,j-1,QV))
          pface    = HALF*(Q(i,j,QPRES) + Q(i,j-1,QPRES))

          fy(i,j,UMX)   = - tauyx
          fy(i,j,UMY)   = - tauyy
          fy(i,j,UMZ)   = 0.d0
          fy(i,j,UEDEN) = - tauyx*Uface(1) - tauyy*Uface(2) - lamy(i,j)*dTdy

          dlnpj(j) = gfacj(j) * (Q(i,j,QPRES) - Q(i,j-1,QPRES)) / pface
       end do

       do j=lo(2)-1,hi(2)+1
          eosj(j) % massfrac(:) = Q(i,j,QFS:QFS+nspecies-1)
          eosj(j) % T           = Q(i,j,QTEMP)
          eosj(j) % rho         = Q(i,j,QRHO)
          eosj(j) % p           = Q(i,j,QPRES)
          call eos_ytx(eosj(j))
          call eos_hi(eosj(j))
       end do

       ! Get species/enthalpy diffusion, compute correction velocity
       Vcj = 0.d0
       do n=1,nspecies
          do j = lo(2), hi(2)+1
             Xface = HALF*(eosj(j)%molefrac(n) + eosj(j-1)%molefrac(n))
             Yface = HALF*(eosj(j)%massfrac(n) + eosj(j-1)%massfrac(n))
             hface = HALF*(eosj(j)%hi(n)       + eosj(j-1)%hi(n))

             dXdy = gfacj(j) * (eosj(j)%molefrac(n) - eosj(j-1)%molefrac(n))
             Vd = -Dy(i,j,n)*(dXdy + (Xface - Yface) * dlnpj(j))

             fy(i,j,UFS+n-1) = Vd
             Vcj(j) = Vcj(j) + Vd
             fy(i,j,UEDEN) = fy(i,j,UEDEN) + Vd*hface
          end do
       end do

       ! Add correction velocity
       do n=1,nspecies
          do j = lo(2), hi(2)+1
             Yface = HALF*(eosj(j)%massfrac(n) + eosj(j-1)%massfrac(n))
             hface = HALF*(eosj(j)%hi(n)       + eosj(j-1)%hi(n))

             fy(i,j,UFS+n-1) = fy(i,j,UFS+n-1) - Yface*Vcj(j)
             fy(i,j,UEDEN)   = fy(i,j,UEDEN)   - Yface*Vcj(j)*hface
          end do
       end do
    end do

    ! Sscale fluxes by area
    do j=lo(2),hi(2)+1
       do i=lo(1),hi(1)
          fy(i,j,UMX)   = fy(i,j,UMX)   * Ay(i,j)
          fy(i,j,UMY)   = fy(i,j,UMY)   * Ay(i,j)
          fy(i,j,UEDEN) = fy(i,j,UEDEN) * Ay(i,j)
          fy(i,j,UFS:UFS+nspecies-1) = fy(i,j,UFS:UFS+nspecies-1) * Ay(i,j)
       enddo
    enddo

    do j=lo(2)+1,hi(2)-1
       call destroy(eosj(j))
    enddo

    do n=1,NVAR
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             D(i,j,n) = - (fx(i+1,j,n)-fx(i,j,n) &
                  +        fy(i,j+1,n)-fy(i,j,n))/V(i,j)
          end do
       end do
    end do

  end subroutine pc_diffterm

end module diffterm_module
