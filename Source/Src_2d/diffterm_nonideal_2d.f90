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
    use meth_params_module, only : NVAR, UMX, UMY, UEDEN, UFS, QVAR, QU, QV, QPRES, QTEMP, QFS, QRHO
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

    integer :: i, j, n, nn
    double precision :: tauxx, tauxy, tauyx, tauyy, divu
    double precision :: Uface(2), dudx,dvdx,dudy,dvdy
    double precision :: pface, hface, Yface
    double precision :: dTdx, dTdy, Vd
    double precision :: Vci(lo(1):hi(1)+1)
    double precision :: ddrivei(lo(1):hi(1)+1,nspecies), gradYi(lo(1):hi(1)+1,nspecies)
    double precision :: gradPi(lo(1):hi(1)+1), dsumi(lo(1):hi(1)+1), gfaci(lo(1):hi(1)+1)
    double precision :: Vcj(lo(2):hi(2)+1)
    double precision :: ddrivej(lo(2):hi(2)+1,nspecies), gradYj(lo(2):hi(2)+1,nspecies)
    double precision :: gradPj(lo(2):hi(2)+1), dsumj(lo(2):hi(2)+1), gfacj(lo(2):hi(2)+1)
    double precision, parameter :: twoThirds = 2.d0/3.d0
    double precision :: dxinv(2)
    type(eos_t) :: eos_statei(lo(1)-1:hi(1)+1)
    type(eos_t) :: eos_statej(lo(2)-1:hi(2)+1)

    dxinv = 1.d0/deltax
    do i=lo(1)-1,hi(1)+1
       call build(eos_statei(i))
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

          Uface(:) = HALF*(Q(i,j,QU:QV) + Q(i-1,j,QU:QV))
          pface    = HALF*(Q(i,j,QPRES) + Q(i-1,j,QPRES))

          fx(i,j,UMX)   = - tauxx

      !   write(6,*)" in diff x ",i,j,tauxx, mux(i,j),xix(i,j),divu, dxinv(1)
      !   write(6,*)" in diff x vels ",Q(i,j,QU),Q(i-1,j,QU),Q(i,j,QV),Q(i-1,j,QV)
      !   stop
          fx(i,j,UMY)   = - tauxy
          fx(i,j,UEDEN) = - tauxx*Uface(1) - tauxy*Uface(2)

          ! thermal conduction
          fx(i,j,UEDEN) = fx(i,j,UEDEN) - lamx(i,j)*dTdx

          ! (1/p)(dp/dx)
          !   dlnp(i) = dxinv(1) * (Q(i,j,QPRES) - Q(i-1,j,QPRES)) / pface
          gradPi(i) = gfaci(i) * (Q(i,j,QPRES) - Q(i-1,j,QPRES))
          Vci(i) = 0.d0
       end do

       do i=lo(1)-1,hi(1)+1
          eos_statei(i) % massfrac(:) = Q(i,j,QFS:QFS+nspecies-1)
          eos_statei(i) % T           = Q(i,j,QTEMP)
          eos_statei(i) % rho         = Q(i,j,QRHO)
          call eos_ytx(eos_statei(i))
          call eos_hi(eos_statei(i))
          call eos_get_transport(eos_statei(i))
       end do

       do n=1,nspecies
          do i = lo(1), hi(1)+1

             gradYi(i,n) = gfaci(i) * (eos_statei(i)%massfrac(n) - eos_statei(i-1)%massfrac(n))

!    put in P term
             ddrivei(i,n) = 0.5d0*(eos_statei(i)% diP(n) + eos_statei(i-1)% diP(n)) * gradPi(i)

          enddo
       enddo

       do n=1,nspecies
       do nn=1,nspecies
          do i = lo(1), hi(1)+1

            ddrivei(i,n) = ddrivei(i,n)+ 0.5d0* (eos_statei(i) % dijY(n,nn) &
                           + eos_statei(i-1) % dijY(n,nn)) * gradYi(i,nn)

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

            ddrivei(i,n) =  ddrivei(i,n) - eos_statei(i)%massfrac(n) * dsumi(i)

          enddo
       enddo


       ! Get species/enthalpy diffusion, compute correction velocity
       do n=1,nspecies
          do i = lo(1), hi(1)+1
             hface = HALF*(eos_statei(i)%hi(n) + eos_statei(i-1)%hi(n))

             Vd = -Dx(i,j,n)*ddrivei(i,n)
             
             fx(i,j,UFS+n-1) = Vd
             Vci(i) = Vci(i) + Vd
             fx(i,j,UEDEN) = fx(i,j,UEDEN) + Vd*hface
          end do
       end do

       ! Add correction velocity
       do n=1,nspecies
          do i = lo(1), hi(1)+1
             Yface = HALF*(eos_statei(i)%massfrac(n) + eos_statei(i-1)%massfrac(n))
             hface = HALF*(eos_statei(i)%hi(n)       + eos_statei(i-1)%hi(n))

             fx(i,j,UFS+n-1) = fx(i,j,UFS+n-1) - Yface*Vci(i)
             fx(i,j,UEDEN)   = fx(i,j,UEDEN)   - Yface*Vci(i)*hface
          end do
       end do
    end do

    ! Scale fluxes by area
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)+1
          fx(i,j,UMX)   = fx(i,j,UMX)   * Ax(i,j)
          fx(i,j,UMY)   = fx(i,j,UMY)   * Ax(i,j)
          fx(i,j,UEDEN) = fx(i,j,UEDEN) * Ax(i,j)
          fx(i,j,UFS:UFS+nspecies-1) = fx(i,j,UFS:UFS+nspecies-1) * Ax(i,j)
       enddo
    enddo

    do i=lo(1)+1,hi(1)-1
       call destroy(eos_statei(i))
    enddo


    do j=lo(2)-1,hi(2)+1
       call build(eos_statej(j))
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

          Uface(:) = HALF*(Q(i,j,QU:QV) + Q(i,j-1,QU:QV))
          pface    = HALF*(Q(i,j,QPRES) + Q(i,j-1,QPRES))

          fy(i,j,UMX)   = - tauyx
          fy(i,j,UMY)   = - tauyy
          fy(i,j,UEDEN) = - tauyx*Uface(1) - tauyy*Uface(2)

          ! thermal conduction
          fy(i,j,UEDEN) = fy(i,j,UEDEN) - lamy(i,j)*dTdy

          ! (1/p)(dp/dy)
          !  dlnp(j) = dxinv(2) * (Q(i,j,QPRES) - Q(i,j-1,QPRES)) / pface
          gradPj(j) = gfacj(j) * (Q(i,j,QPRES) - Q(i,j-1,QPRES))
          Vcj(j) = 0.d0
!           if(i.eq.1.and. j.eq.207)then
!                write(6,*)" in diffterm ",j,Q(i,j,QTEMP), muy(i,j), xiy(i,j), lamy(i,j)
!           endif
       end do

       do j=lo(2)-1,hi(2)+1
          eos_statej(j) % massfrac(:) = Q(i,j,QFS:QFS+nspecies-1)
          eos_statej(j) % T           = Q(i,j,QTEMP)
          eos_statej(j) % rho         = Q(i,j,QRHO)
          call eos_ytx(eos_statej(j))
          call eos_hi(eos_statej(j))
          call eos_get_transport(eos_statej(j))
       end do

       do n=1,nspecies
          do j = lo(2), hi(2)+1

             gradYj(j,n) = gfacj(j) * (eos_statej(j)%massfrac(n) - eos_statej(j-1)%massfrac(n))

!    put in P term
             ddrivej(j,n) = 0.5d0*(eos_statej(j)% diP(n) + eos_statej(j-1)% diP(n)) * gradPj(j)

!            if(i.eq.1 .and. j.eq.207)then
!                write(6,*)" in species pterm ",n,ddrive(j,n)
!            endif
!            ddrive(j,n) = ddrive(j,n)- 0.5d0*( eos_state(j)%massfrac(n)* eos_state(j)%wbar/(Ru*eos_state(j) % T * eos_state(j) % rho ) &
!                        + eos_state(j-1)%massfrac(n)* eos_state(j-1)%wbar/(Ru*eos_state(j-1) % T * eos_state(j-1) % rho )) * gradP(j)


          enddo
       enddo

       do n=1,nspecies
       do nn=1,nspecies
          do j = lo(2), hi(2)+1

            ddrivej(j,n) = ddrivej(j,n)+ 0.5d0* (eos_statej(j) % dijY(n,nn) + eos_statej(j-1) % dijY(n,nn)) * gradYj(j,nn)

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

            ddrivej(j,n) =  ddrivej(j,n) - eos_statej(j)%massfrac(n) * dsumj(j)
!            if(i.eq.1 .and. j.eq.207)then
!                write(6,*)" species total ",n,ddrive(j,n), dsum(j)
!            endif

          enddo
       enddo

       ! Get species/enthalpy diffusion, compute correction velocity
       do n=1,nspecies
          do j = lo(2), hi(2)+1
             hface = HALF*(eos_statej(j)%hi(n)       + eos_statej(j-1)%hi(n))

             Vd = -Dy(i,j,n)*ddrivej(j,n)
             
             fy(i,j,UFS+n-1) = Vd
             Vcj(j) = Vcj(j) + Vd
             fy(i,j,UEDEN) = fy(i,j,UEDEN) + Vd*hface
   !         if(i.eq.1 .and. j.eq.207)then
   !             write(6,*)" species diff coeff, vell ",n,Dy(i,j,n),Vd
   !         endif
          end do
       end do

       ! Add correction velocity
       do n=1,nspecies
          do j = lo(2), hi(2)+1
             Yface = HALF*(eos_statej(j)%massfrac(n) + eos_statej(j-1)%massfrac(n))
             hface = HALF*(eos_statej(j)%hi(n)       + eos_statej(j-1)%hi(n))

             fy(i,j,UFS+n-1) = fy(i,j,UFS+n-1) - Yface*Vcj(j)
             fy(i,j,UEDEN)   = fy(i,j,UEDEN)   - Yface*Vcj(j)*hface
!          if(i.eq.1.and.j.eq.207)then
!                write(6,*)" in species flux ",n, fy(i,j,UFS+n-1)
!            endif
          end do
       end do
    end do

    ! Scale fluxes by area
    do j=lo(2),hi(2)+1
       do i=lo(1),hi(1)
          fy(i,j,UMX)   = fy(i,j,UMX)   * Ay(i,j)
          fy(i,j,UMY)   = fy(i,j,UMY)   * Ay(i,j)
          fy(i,j,UEDEN) = fy(i,j,UEDEN) * Ay(i,j)
          fy(i,j,UFS:UFS+nspecies-1) = fy(i,j,UFS:UFS+nspecies-1) * Ay(i,j)
       enddo
    enddo

    do j=lo(2)+1,hi(2)-1
       call destroy(eos_statej(j))
    enddo

    do n=1,NVAR
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             D(i,j,n) = - (fx(i+1,j,n) - fx(i,j,n) &
                          +fy(i,j+1,n) - fy(i,j,n) )/V(i,j)
          end do
       end do
    end do

  end subroutine pc_diffterm

end module diffterm_module
