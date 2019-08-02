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
                         Dz,  Dzlo,  Dzhi,&  
                         muz, muzlo, muzhi,& 
                         xiz, xizlo, xizhi,& 
                         lamz,lamzlo,lamzhi,&
                         tz,  tzlo,  tzhi,&
                         Az,  Azlo,  Azhi,&
                         fz,  fzlo,  fzhi,&
                         V,   Vlo,   Vhi,&
                         D,   Dlo,   Dhi,&
                         deltax) bind(C, name = "pc_diffterm")

    use network, only : nspecies
    use meth_params_module, only : NVAR, UMX, UMY, UMZ, UEDEN, UFS, QVAR, QU, QV, QW, QPRES, QTEMP, QFS, QRHO
    use amrex_constants_module
    use eos_type_module
    use eos_module
    use prob_params_module, only : physbc_lo, physbc_hi

    implicit none

    integer, intent(in) ::     lo(3),    hi(3)
    integer, intent(in) ::  dmnlo(3), dmnhi(3)
    integer, intent(in) ::    Qlo(3),   Qhi(3)

    integer, intent(in) ::   Dxlo(3),  Dxhi(3)
    integer, intent(in) ::  muxlo(3), muxhi(3)
    integer, intent(in) ::  xixlo(3), xixhi(3)
    integer, intent(in) :: lamxlo(3),lamxhi(3)
    integer, intent(in) ::   txlo(3),  txhi(3)
    integer, intent(in) ::   Axlo(3),  Axhi(3)
    integer, intent(in) ::   fxlo(3),  fxhi(3)

    integer, intent(in) ::   Dylo(3),  Dyhi(3)
    integer, intent(in) ::  muylo(3), muyhi(3)
    integer, intent(in) ::  xiylo(3), xiyhi(3)
    integer, intent(in) :: lamylo(3),lamyhi(3)
    integer, intent(in) ::   tylo(3),  tyhi(3)
    integer, intent(in) ::   Aylo(3),  Ayhi(3)
    integer, intent(in) ::   fylo(3),  fyhi(3)

    integer, intent(in) ::   Dzlo(3),  Dzhi(3)
    integer, intent(in) ::  muzlo(3), muzhi(3)
    integer, intent(in) ::  xizlo(3), xizhi(3)
    integer, intent(in) :: lamzlo(3),lamzhi(3)
    integer, intent(in) ::   tzlo(3),  tzhi(3)
    integer, intent(in) ::   Azlo(3),  Azhi(3)
    integer, intent(in) ::   fzlo(3),  fzhi(3)

    integer, intent(in) ::    Dlo(3),   Dhi(3)
    integer, intent(in) ::    Vlo(3),   Vhi(3)

    double precision, intent(in   ) ::    Q(   Qlo(1):   Qhi(1),   Qlo(2):   Qhi(2),   Qlo(3):   Qhi(3), QVAR)
    double precision, intent(in   ) ::   Dx(  Dxlo(1):  Dxhi(1),  Dxlo(2):  Dxhi(2),  Dxlo(3):  Dxhi(3), nspecies)
    double precision, intent(in   ) ::  mux( muxlo(1): muxhi(1), muxlo(2): muxhi(2), muxlo(3): muxhi(3)  )
    double precision, intent(in   ) ::  xix( xixlo(1): xixhi(1), xixlo(2): xixhi(2), xixlo(3): xixhi(3)  )
    double precision, intent(in   ) :: lamx(lamxlo(1):lamxhi(1),lamxlo(2):lamxhi(2),lamxlo(3):lamxhi(3)  )
    double precision, intent(in   ) ::   tx(  txlo(1):  txhi(1),  txlo(2):  txhi(2),  txlo(3):  txhi(3), 4)
    double precision, intent(in   ) ::   Ax(  Axlo(1):  Axhi(1),  Axlo(2):  Axhi(2),  Axlo(3):  Axhi(3)  )
    double precision, intent(inout) ::   fx(  fxlo(1):  fxhi(1),  fxlo(2):  fxhi(2),  fxlo(3):  fxhi(3), NVAR)
    double precision, intent(in   ) ::   Dy(  Dylo(1):  Dyhi(1),  Dylo(2):  Dyhi(2),  Dylo(3):  Dyhi(3), nspecies)
    double precision, intent(in   ) ::  muy( muylo(1): muyhi(1), muylo(2): muyhi(2), muylo(3): muyhi(3)  )
    double precision, intent(in   ) ::  xiy( xiylo(1): xiyhi(1), xiylo(2): xiyhi(2), xiylo(3): xiyhi(3)  )
    double precision, intent(in   ) :: lamy(lamylo(1):lamyhi(1),lamylo(2):lamyhi(2),lamylo(3):lamyhi(3)  )
    double precision, intent(in   ) ::   ty(  tylo(1):  tyhi(1),  tylo(2):  tyhi(2),  tylo(3):  tyhi(3), 4)
    double precision, intent(in   ) ::   Ay(  Aylo(1):  Ayhi(1),  Aylo(2):  Ayhi(2),  Aylo(3):  Ayhi(3)  )
    double precision, intent(inout) ::   fy(  fylo(1):  fyhi(1),  fylo(2):  fyhi(2),  fylo(3):  fyhi(3), NVAR)
    double precision, intent(in   ) ::   Dz(  Dzlo(1):  Dzhi(1),  Dzlo(2):  Dzhi(2),  Dzlo(3):  Dzhi(3), nspecies)
    double precision, intent(in   ) ::  muz( muzlo(1): muzhi(1), muzlo(2): muzhi(2), muzlo(3): muzhi(3)  )
    double precision, intent(in   ) ::  xiz( xizlo(1): xizhi(1), xizlo(2): xizhi(2), xizlo(3): xizhi(3)  )
    double precision, intent(in   ) :: lamz(lamzlo(1):lamzhi(1),lamzlo(2):lamzhi(2),lamzlo(3):lamzhi(3)  )
    double precision, intent(in   ) ::   tz(  tzlo(1):  tzhi(1),  tzlo(2):  tzhi(2),  tzlo(3):  tzhi(3), 4)
    double precision, intent(in   ) ::   Az(  Azlo(1):  Azhi(1),  Azlo(2):  Azhi(2),  Azlo(3):  Azhi(3)  )
    double precision, intent(inout) ::   fz(  fzlo(1):  fzhi(1),  fzlo(2):  fzhi(2),  fzlo(3):  fzhi(3), NVAR)
    double precision, intent(inout) ::    D(   Dlo(1):   Dhi(1),   Dlo(2):   Dhi(2),   Dlo(3):   Dhi(3), NVAR)
    double precision, intent(in   ) ::    V(   Vlo(1):   Vhi(1),   Vlo(2):   Vhi(2),   Vlo(3):   Vhi(3)  )
    double precision, intent(in   ) :: deltax(3)

    integer :: i, j, k, n, nn
    double precision :: tauxx, tauxy, tauxz, tauyx, tauyy, tauyz, tauzx, tauzy, tauzz, divu
    double precision :: Uface(3), dudx,dvdx,dwdx,dudy,dvdy,dwdy,dudz,dvdz,dwdz
    double precision :: hface, Yface
    double precision :: dTdx, dTdy, dTdz, Vd

    double precision :: gradPi(lo(1):hi(1)+1), dsumi(lo(1):hi(1)+1)
    double precision :: Vci(lo(1):hi(1)+1), gfaci(lo(1):hi(1)+1)
    double precision :: ddrivei(lo(1):hi(1)+1,nspecies), gradYi(lo(1):hi(1)+1,nspecies)

    double precision :: Vcj(lo(2):hi(2)+1), gfacj(lo(2):hi(2)+1)
    double precision :: ddrivej(lo(2):hi(2)+1,nspecies), gradYj(lo(2):hi(2)+1,nspecies)
    double precision :: gradPj(lo(2):hi(2)+1), dsumj(lo(2):hi(2)+1)

    double precision :: Vck(lo(3):hi(3)+1), gfack(lo(3):hi(3)+1)
    double precision :: ddrivek(lo(3):hi(3)+1,nspecies), gradYk(lo(3):hi(3)+1,nspecies)
    double precision :: gradPk(lo(3):hi(3)+1), dsumk(lo(3):hi(3)+1)

    double precision, parameter :: twoThirds = 2.d0/3.d0
    double precision :: dxinv(3)
    type(eos_t) :: eosi(lo(1)-1:hi(1)+1)
    type(eos_t) :: eosj(lo(2)-1:hi(2)+1)
    type(eos_t) :: eosk(lo(3)-1:hi(3)+1)

    dxinv = 1.d0/deltax

    gfaci = dxinv(1)
    gfacj = dxinv(2)
    gfack = dxinv(3)
   
    do i=lo(1)-1,hi(1)+1
       call build(eosi(i))
    enddo

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)+1
             dTdx = gfaci(i) * (Q(i,j,k,QTEMP) - Q(i-1,j,k,QTEMP))
             dudx = gfaci(i) * (Q(i,j,k,QU)    - Q(i-1,j,k,QU))
             dvdx = gfaci(i) * (Q(i,j,k,QV)    - Q(i-1,j,k,QV))
             dwdx = gfaci(i) * (Q(i,j,k,QW)    - Q(i-1,j,k,QW))
             dudy = tx(i,j,k,1)
             dvdy = tx(i,j,k,2)
             dudz = tx(i,j,k,3)
             dwdz = tx(i,j,k,4)

             divu = dudx + dvdy + dwdz
             tauxx = mux(i,j,k)*(2.d0*dudx-twoThirds*divu) + xix(i,j,k)*divu
             tauxy = mux(i,j,k)*(dudy+dvdx)
             tauxz = mux(i,j,k)*(dudz+dwdx)

             Uface(1) = HALF*(Q(i,j,k,QU) + Q(i-1,j,k,QU))
             Uface(2) = HALF*(Q(i,j,k,QV) + Q(i-1,j,k,QV))
             Uface(3) = HALF*(Q(i,j,k,QW) + Q(i-1,j,k,QW))

             fx(i,j,k,UMX)   = - tauxx
             fx(i,j,k,UMY)   = - tauxy
             fx(i,j,k,UMZ)   = - tauxz
             fx(i,j,k,UEDEN) = - tauxx*Uface(1) - tauxy*Uface(2) - tauxz*Uface(3)

             fx(i,j,k,UEDEN) = fx(i,j,k,UEDEN) - lamx(i,j,k)*dTdx
             gradPi(i) = gfaci(i) * (Q(i,j,k,QPRES) - Q(i-1,j,k,QPRES)) 
             Vci(i) = 0.d0
          end do
          do i=lo(1)-1,hi(1)+1
             eosi(i) % massfrac(:) = Q(i,j,k,QFS:QFS+nspecies-1)
             eosi(i) % T           = Q(i,j,k,QTEMP)
             eosi(i) % rho         = Q(i,j,k,QRHO)
             call eos_ytx(eosi(i))
             call eos_hi(eosi(i))
             call eos_get_transport(eosi(i))
          end do
          do n=1,nspecies
             do i = lo(1), hi(1)+1
                gradYi(i,n) = gfaci(i) * (eosi(i)%massfrac(n) - eosi(i-1)%massfrac(n))
!    put in P term
                ddrivei(i,n) = 0.5d0*(eosi(i)% diP(n) + eosi(i-1)% diP(n)) * gradPi(i)
             enddo
          enddo
          do n=1,nspecies
            do nn=1,nspecies
              do i = lo(1), hi(1)+1
                ddrivei(i,n) = ddrivei(i,n)+ 0.5d0* (eosi(i) % dijY(n,nn) &
                             + eosi(i-1) % dijY(n,nn)) * gradYi(i,nn)
              enddo
            enddo
          enddo
          dsumi = 0.d0
          do n=1,nspecies
             do i = lo(1), hi(1)+1
               dsumi(i) = dsumi(i) + ddrivei(i,n)
             enddo
          enddo
          do n=1,nspecies
             do i = lo(1), hi(1)+1
               ddrivei(i,n) =  ddrivei(i,n) - eosi(i)%massfrac(n) * dsumi(i)
             enddo
          enddo
          ! Get species/enthalpy diffusion, compute correction velocity
          do n=1,nspecies
             do i = lo(1), hi(1)+1
                hface = HALF*(eosi(i)%hi(n) + eosi(i-1)%hi(n))
                Vd = -Dx(i,j,k,n)*ddrivei(i,n)
                fx(i,j,k,UFS+n-1) = Vd
                Vci(i) = Vci(i) + Vd
                fx(i,j,k,UEDEN) = fx(i,j,k,UEDEN) + Vd*hface
             end do
          end do
          ! Add correction velocity
          do n=1,nspecies
             do i = lo(1), hi(1)+1
                Yface = HALF*(eosi(i)%massfrac(n) + eosi(i-1)%massfrac(n))
                hface = HALF*(eosi(i)%hi(n)       + eosi(i-1)%hi(n))

                fx(i,j,k,UFS+n-1) = fx(i,j,k,UFS+n-1) - Yface*Vci(i)
                fx(i,j,k,UEDEN)   = fx(i,j,k,UEDEN)   - Yface*Vci(i)*hface
             end do
          end do
       end do
    end do
    
    ! Sscale fluxes by area
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)+1
             fx(i,j,k,UMX)   = fx(i,j,k,UMX)   * Ax(i,j,k)
             fx(i,j,k,UMY)   = fx(i,j,k,UMY)   * Ax(i,j,k)
             fx(i,j,k,UMZ)   = fx(i,j,k,UMZ)   * Ax(i,j,k)
             fx(i,j,k,UEDEN) = fx(i,j,k,UEDEN) * Ax(i,j,k)
             fx(i,j,k,UFS:UFS+nspecies-1) = fx(i,j,k,UFS:UFS+nspecies-1) * Ax(i,j,k)
          enddo
       enddo
    enddo

    do i=lo(1)+1,hi(1)-1
       call destroy(eosi(i))
    enddo

    do j=lo(2)-1,hi(2)+1
       call build(eosj(j))
    enddo

    do k=lo(3),hi(3)
       do i=lo(1),hi(1)
          do j=lo(2),hi(2)+1
             dTdy = gfacj(j) * (Q(i,j,k,QTEMP) - Q(i,j-1,k,QTEMP))
             dudy = gfacj(j) * (Q(i,j,k,QU)    - Q(i,j-1,k,QU))
             dvdy = gfacj(j) * (Q(i,j,k,QV)    - Q(i,j-1,k,QV))
             dwdy = gfacj(j) * (Q(i,j,k,QW)    - Q(i,j-1,k,QW))
             dudx = ty(i,j,k,1)
             dvdx = ty(i,j,k,2)
             dvdz = ty(i,j,k,3)
             dwdz = ty(i,j,k,4)

             divu = dudx + dvdy + dwdz
             tauyx = muy(i,j,k)*(dudy+dvdx)
             tauyy = muy(i,j,k)*(2.d0*dvdy-twoThirds*divu) + xiy(i,j,k)*divu
             tauyz = muy(i,j,k)*(dwdy+dvdz)
             Uface(1) = HALF*(Q(i,j,k,QU) + Q(i,j-1,k,QU))
             Uface(2) = HALF*(Q(i,j,k,QV) + Q(i,j-1,k,QV))
             Uface(3) = HALF*(Q(i,j,k,QW) + Q(i,j-1,k,QW))

             fy(i,j,k,UMX)   = - tauyx
             fy(i,j,k,UMY)   = - tauyy
             fy(i,j,k,UMZ)   = - tauyz
             fy(i,j,k,UEDEN) = - tauyx*Uface(1) - tauyy*Uface(2) - tauyz*Uface(3)

             fy(i,j,k,UEDEN) = fy(i,j,k,UEDEN) - lamy(i,j,k)*dTdy
             gradPj(j) = gfacj(j) * (Q(i,j,k,QPRES) - Q(i,j-1,k,QPRES)) 
             Vcj(j) = 0.d0
          end do
          do j=lo(2)-1,hi(2)+1
             eosj(j) % massfrac(:) = Q(i,j,k,QFS:QFS+nspecies-1)
             eosj(j) % T           = Q(i,j,k,QTEMP)
             eosj(j) % rho         = Q(i,j,k,QRHO)
             call eos_ytx(eosj(j))
             call eos_hi(eosj(j))
             call eos_get_transport(eosj(j))
          end do
          do n=1,nspecies
            do j = lo(2), hi(2)+1
             gradYj(j,n) = gfacj(j) * (eosj(j)%massfrac(n) - eosj(j-1)%massfrac(n))
!    put in P term
             ddrivej(j,n) = 0.5d0*(eosj(j)% diP(n) + eosj(j-1)% diP(n)) * gradPj(j)
!            if(i.eq.1 .and. j.eq.207)then
!                write(6,*)" in species pterm ",n,ddrivej(j,n)
!            endif
!            ddrivej(j,n) = ddrivej(j,n)- 0.5d0*( eosj(j)%massfrac(n)* eosj(j)%wbar/(Ru*eosj(j) % T * eosj(j) % rho ) &
!                        + eosj(j-1)%massfrac(n)* eosj(j-1)%wbar/(Ru*eosj(j-1) % T * eosj(j-1) % rho )) * gradP(j)
            enddo
          enddo
          do n=1,nspecies
            do nn=1,nspecies
              do j = lo(2), hi(2)+1
                ddrivej(j,n) = ddrivej(j,n)+ 0.5d0* (eosj(j) % dijY(n,nn) + eosj(j-1) % dijY(n,nn)) * gradYj(j,nn)
              enddo
            enddo
          enddo
          dsumj = 0.d0
          do n=1,nspecies
             do j = lo(2), hi(2)+1
               dsumj(j) = dsumj(j) + ddrivej(j,n)
             enddo
          enddo
          do n=1,nspecies
             do j = lo(2), hi(2)+1
               ddrivej(j,n) =  ddrivej(j,n) - eosj(j)%massfrac(n) * dsumj(j)
!               if(i.eq.1 .and. j.eq.207)then
!                   write(6,*)" species total ",n,ddrivej(j,n), dsumj(j)
!               endif
             enddo
          enddo
          ! Get species/enthalpy diffusion, compute correction velocity
          do n=1,nspecies
             do j = lo(2), hi(2)+1
                hface = HALF*(eosj(j)%hi(n)       + eosj(j-1)%hi(n))
                Vd = -Dy(i,j,k,n)*ddrivej(j,n)
                fy(i,j,k,UFS+n-1) = Vd
                Vcj(j) = Vcj(j) + Vd
                fy(i,j,k,UEDEN) = fy(i,j,k,UEDEN) + Vd*hface
             end do
          end do
          ! Add correction velocity
          do n=1,nspecies
             do j = lo(2), hi(2)+1
                Yface = HALF*(eosj(j)%massfrac(n) + eosj(j-1)%massfrac(n))
                hface = HALF*(eosj(j)%hi(n)       + eosj(j-1)%hi(n))

                fy(i,j,k,UFS+n-1) = fy(i,j,k,UFS+n-1) - Yface*Vcj(j)
                fy(i,j,k,UEDEN)   = fy(i,j,k,UEDEN)   - Yface*Vcj(j)*hface
             end do
          end do
       end do
    end do
    ! Sscale fluxes by area
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)
             fy(i,j,k,UMX)   = fy(i,j,k,UMX)   * Ay(i,j,k)
             fy(i,j,k,UMY)   = fy(i,j,k,UMY)   * Ay(i,j,k)
             fy(i,j,k,UMZ)   = fy(i,j,k,UMZ)   * Ay(i,j,k)
             fy(i,j,k,UEDEN) = fy(i,j,k,UEDEN) * Ay(i,j,k)
             fy(i,j,k,UFS:UFS+nspecies-1) = fy(i,j,k,UFS:UFS+nspecies-1) * Ay(i,j,k)
          enddo
       enddo
    enddo

    do j=lo(2)+1,hi(2)-1
       call destroy(eosj(j))
    enddo

    do k=lo(3)-1,hi(3)+1
       call build(eosk(k))
    enddo

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          do k=lo(3),hi(3)+1
             dTdz = gfack(k) * (Q(i,j,k,QTEMP) - Q(i,j,k-1,QTEMP))
             dudz = gfack(k) * (Q(i,j,k,QU)    - Q(i,j,k-1,QU))
             dvdz = gfack(k) * (Q(i,j,k,QV)    - Q(i,j,k-1,QV))
             dwdz = gfack(k) * (Q(i,j,k,QW)    - Q(i,j,k-1,QW))
             dudx = tz(i,j,k,1)
             dwdx = tz(i,j,k,2)
             dvdy = tz(i,j,k,3)
             dwdy = tz(i,j,k,4)

             divu = dudx + dvdy + dwdz
             tauzx = muz(i,j,k)*(dudz+dwdx)
             tauzy = muz(i,j,k)*(dvdz+dwdy)
             tauzz = muz(i,j,k)*(2.d0*dwdz-twoThirds*divu) + xiz(i,j,k)*divu
             Uface(1) = HALF*(Q(i,j,k,QU) + Q(i,j,k-1,QU))
             Uface(2) = HALF*(Q(i,j,k,QV) + Q(i,j,k-1,QV))
             Uface(3) = HALF*(Q(i,j,k,QW) + Q(i,j,k-1,QW))

             fz(i,j,k,UMX)   = - tauzx
             fz(i,j,k,UMY)   = - tauzy
             fz(i,j,k,UMZ)   = - tauzz
             fz(i,j,k,UEDEN) = - tauzx*Uface(1) - tauzy*Uface(2) - tauzz*Uface(3)

             fz(i,j,k,UEDEN) = fz(i,j,k,UEDEN) - lamz(i,j,k)*dTdz
             gradPk(k) = gfack(k) * (Q(i,j,k,QPRES) - Q(i,j,k-1,QPRES)) 
             Vck(k) = 0.d0
          end do
          do k=lo(3)-1,hi(3)+1
             eosk(k) % massfrac(:) = Q(i,j,k,QFS:QFS+nspecies-1)
             eosk(k) % T           = Q(i,j,k,QTEMP)
             eosk(k) % rho         = Q(i,j,k,QRHO)
             call eos_ytx(eosk(k))
             call eos_hi(eosk(k))
             call eos_get_transport(eosk(k))
          end do
          do n=1,nspecies
            do k = lo(3), hi(3)+1
              gradYk(k,n) = gfack(k) * (eosk(k)%massfrac(n) - eosk(k-1)%massfrac(n))
!    put in P term
              ddrivek(k,n) = 0.5d0*(eosk(k)% diP(n) + eosk(k-1)% diP(n)) * gradPk(k)
            enddo
          enddo
          do n=1,nspecies
            do nn=1,nspecies
              do k = lo(3), hi(3)+1
                ddrivek(k,n) = ddrivek(k,n)+ 0.5d0* (eosk(k) % dijY(n,nn) + eosk(k-1) % dijY(n,nn)) * gradYk(k,nn)
              enddo
            enddo
          enddo
          dsumk = 0.d0
          do n=1,nspecies
             do k = lo(3), hi(3)+1
               dsumk(k) = dsumk(k) + ddrivek(k,n)
             enddo
          enddo
          do n=1,nspecies
             do k = lo(3), hi(3)+1
               ddrivek(k,n) =  ddrivek(k,n) - eosk(k)%massfrac(n) * dsumk(k)
             enddo
          enddo
          ! Get species/enthalpy diffusion, compute correction velocity
          do n=1,nspecies
             do k = lo(3), hi(3)+1
                hface = HALF*(eosk(k)%hi(n)       + eosk(k-1)%hi(n))
                Vd = -Dz(i,j,k,n)*ddrivek(k,n)
                fz(i,j,k,UFS+n-1) = Vd
                Vck(k) = Vck(k) + Vd
                fz(i,j,k,UEDEN) = fz(i,j,k,UEDEN) + Vd*hface
             end do
          end do
          ! Add correction velocity
          do n=1,nspecies
             do k = lo(3), hi(3)+1
                Yface = HALF*(eosk(k)%massfrac(n) + eosk(k-1)%massfrac(n))
                hface = HALF*(eosk(k)%hi(n)       + eosk(k-1)%hi(n))
                fz(i,j,k,UFS+n-1) = fz(i,j,k,UFS+n-1) - Yface*Vck(k)
                fz(i,j,k,UEDEN)   = fz(i,j,k,UEDEN)   - Yface*Vck(k)*hface
             end do
          end do
       end do
    end do

    ! Sscale fluxes by area
    do k=lo(3),hi(3)+1
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             fz(i,j,k,UMX)   = fz(i,j,k,UMX)   * Az(i,j,k)
             fz(i,j,k,UMY)   = fz(i,j,k,UMY)   * Az(i,j,k)
             fz(i,j,k,UMZ)   = fz(i,j,k,UMZ)   * Az(i,j,k)
             fz(i,j,k,UEDEN) = fz(i,j,k,UEDEN) * Az(i,j,k)
             fz(i,j,k,UFS:UFS+nspecies-1) = fz(i,j,k,UFS:UFS+nspecies-1) * Az(i,j,k)
          enddo
       enddo
    enddo

    do k=lo(3)+1,hi(3)-1
       call destroy(eosk(k))
    enddo

    do n=1,NVAR
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                D(i,j,k,n) = - (fx(i+1,j,k,n)-fx(i,j,k,n) &
                               +fy(i,j+1,k,n)-fy(i,j,k,n) &
                               +fz(i,j,k+1,n)-fz(i,j,k,n) )/V(i,j,k)
             end do
          end do
       end do
    end do

  end subroutine pc_diffterm

end module diffterm_module
