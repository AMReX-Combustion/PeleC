module forcing_src_module

#include <AMReX_CONSTANTS.H>

  implicit none

  private

  public :: pc_forcing_src

contains

  subroutine pc_forcing_src(lo,hi, &
                            old_state,os_lo,os_hi, &
                            new_state,ns_lo,ns_hi, &
                            src,src_lo,src_hi,problo,dx,xlo,xhi, &
                            time,dt) bind(C, name = "pc_forcing_src")

    use bl_constants_module, only: ZERO, HALF
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ
    use prob_params_module, only: dim
    use probdata_module

    implicit none

    integer,          intent(in   ) :: lo(3),hi(3)
    integer,          intent(in   ) :: os_lo(3),os_hi(3)
    integer,          intent(in   ) :: ns_lo(3),ns_hi(3)
    integer,          intent(in   ) :: src_lo(3),src_hi(3)
    double precision, intent(in   ) :: xlo(3), xhi(3)
    double precision, intent(in   ) :: old_state(os_lo(1):os_hi(1),os_lo(2):os_hi(2),os_lo(3):os_hi(3),NVAR)
    double precision, intent(in   ) :: new_state(ns_lo(1):ns_hi(1),ns_lo(2):ns_hi(2),ns_lo(3):ns_hi(3),NVAR)
    double precision, intent(inout) :: src(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),NVAR)
    double precision, intent(in   ) :: problo(3),dx(3),time,dt

    ! local
    integer :: i,j,k
    integer :: kx, ky, kz
    integer :: modx, mody, modz, xstep, ystep, zstep
    integer :: nxmodes, nymodes, nzmodes
    double precision :: f1, f2, f3
    double precision :: HLx, HLy, HLz, hx, hy, hz 
    double precision :: infl_time, kappa, kappaMax
    double precision :: kxd, kyd, kzd
    double precision :: Lmin, twicePi, xT, x, y, z


!write(*,*) 'DEBUG ENTERING FORCING'



    src(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = ZERO

!     Homogeneous Isotropic Turbulence
    twicePi=two*Pi

!     Adjust z offset for probtype 15
    if (time_offset.gt.zero) then
      infl_time = time + time_offset
    else
      infl_time = time
    endif

    Lmin = min(Lx,Ly,Lz)

    kappaMax = dfloat(nmodes)/Lmin + 1.0d-8
    nxmodes = nmodes*int(0.5+Lx/Lmin)
    nymodes = nmodes*int(0.5+Ly/Lmin)
    nzmodes = nmodes*int(0.5+Lz/Lmin)

    xstep = int(Lx/Lmin+0.5)
    ystep = int(Ly/Lmin+0.5)
    zstep = int(Lz/Lmin+0.5)

    HLx = Lx
    HLy = Ly
    HLz = Lz

!    write(*,*) 'DEBUG FORCING',lo,hi,xlo
!    write(*,*) 'DEBUG FORCING 2',infl_time,div_free_force,mode_start,nmodes,nxmodes,nymodes,nzmodes,xstep,ystep,zstep,HLx,HLy,HLz 

!               do kz = mode_start*zstep, nzmodes, zstep
!                  kzd = dfloat(kz)
!                  do ky = mode_start*ystep, nymodes, ystep
!                     kyd = dfloat(ky)
!                     do kx = mode_start*xstep, nxmodes, xstep
!                        kxd = dfloat(kx)
!                        kappa = sqrt( (kxd*kxd)/(Lx*Lx) + (kyd*kyd)/(Ly*Ly) + (kzd*kzd)/(Lz*Lz) )
!                        if (kappa.le.kappaMax) then
!                           write (*,*) "Mode"
!                           write (*,*) "kappa = ",kx,ky,kz,kappa
!                           write (*,*) "Amplitudes - C"
!                           write (*,*) FAX(kx,ky,kz), FAY(kx,ky,kz), FAZ(kx,ky,kz)
!                           write (*,*) "Frequencies"
!                           write (*,*) FTX(kx,ky,kz), FTY(kx,ky,kz), FTZ(kx,ky,kz)
!                           write (*,*) "Phases"
!                           write (*,*) FPX(kx,ky,kz), FPY(kx,ky,kz), FPZ(kx,ky,kz)
!                           write (*,*) "TAT"
!                           write (*,*) TAT(kx,ky,kz)
!                        endif
!                     enddo
!                  enddo
!               enddo


!write(*,*) ""
!write(*,*) "OK UNTIL HERE",infl_time
!write(*,*) ""

!write(*,*) 'DEBUG lo,hi,xlo',lo,hi,xlo

    do k = lo(3), hi(3)
       z = xlo(3) + dx(3)*(dble(k-lo(3)) + HALF)
       do j = lo(2), hi(2)
          y = xlo(2) + dx(2)*(dble(j-lo(2)) + HALF)
          do i = lo(1), hi(1)
             x = xlo(1) + dx(1)*(dble(i-lo(1)) + HALF)
             f1 = zero
             f2 = zero
             f3 = zero

             do kz = mode_start*zstep, nzmodes, zstep
                kzd = dfloat(kz)
                do ky = mode_start*ystep, nymodes, ystep
                   kyd = dfloat(ky)
                   do kx = mode_start*xstep, nxmodes, xstep
                      kxd = dfloat(kx)
                      kappa = sqrt( (kxd*kxd)/(Lx*Lx) + (kyd*kyd)/(Ly*Ly) + (kzd*kzd)/(Lz*Lz) )
                      if (kappa.le.kappaMax) then
                         xT = cos(FTX(kx,ky,kz)*infl_time+TAT(kx,ky,kz))
                         if (div_free_force.eq.1) then
                             f1 = f1 + xT * ( FAZ(kx,ky,kz)*twicePi*(kyd/HLy) &
                                          * sin(twicePi*kxd*x/HLx+FPZX(kx,ky,kz)) &
                                          * cos(twicePi*kyd*y/HLy+FPZY(kx,ky,kz)) &
                                          * sin(twicePi*kzd*z/HLz+FPZZ(kx,ky,kz)) &
                                         -    FAY(kx,ky,kz)*twicePi*(kzd/HLz) &
                                          * sin(twicePi*kxd*x/HLx+FPYX(kx,ky,kz)) &
                                          * sin(twicePi*kyd*y/HLy+FPYY(kx,ky,kz)) &
                                          * cos(twicePi*kzd*z/HLz+FPYZ(kx,ky,kz)) )

                             f2 = f2 + xT * ( FAX(kx,ky,kz)*twicePi*(kzd/HLz) &
                                          * sin(twicePi*kxd*x/HLx+FPXX(kx,ky,kz)) &
                                          * sin(twicePi*kyd*y/HLy+FPXY(kx,ky,kz)) &
                                          * cos(twicePi*kzd*z/HLz+FPXZ(kx,ky,kz)) &
                                         -    FAZ(kx,ky,kz)*twicePi*(kxd/HLx) &
                                          * cos(twicePi*kxd*x/HLx+FPZX(kx,ky,kz)) &
                                          * sin(twicePi*kyd*y/HLy+FPZY(kx,ky,kz)) &
                                          * sin(twicePi*kzd*z/HLz+FPZZ(kx,ky,kz)) )

                             f3 = f3 + xT * ( FAY(kx,ky,kz)*twicePi*(kxd/HLx) &
                                          * cos(twicePi*kxd*x/HLx+FPYX(kx,ky,kz)) &
                                          * sin(twicePi*kyd*y/HLy+FPYY(kx,ky,kz)) &
                                          * sin(twicePi*kzd*z/HLz+FPYZ(kx,ky,kz)) &
                                         -    FAX(kx,ky,kz)*twicePi*(kyd/HLy) &
                                          * sin(twicePi*kxd*x/HLx+FPXX(kx,ky,kz)) &
                                          * cos(twicePi*kyd*y/HLy+FPXY(kx,ky,kz)) &
                                          * sin(twicePi*kzd*z/HLz+FPXZ(kx,ky,kz)) )

                         else
                             f1 = f1 + xT*FAX(kx,ky,kz)* cos(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) &
                                                       * sin(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) &
                                                       * sin(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
                             f2 = f2 + xT*FAY(kx,ky,kz)* sin(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) &
                                                       * cos(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) &
                                                       * sin(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))
                             f3 = f3 + xT*FAZ(kx,ky,kz)* sin(twicePi*kxd*x/HLx+FPX(kx,ky,kz)) &
                                                       * sin(twicePi*kyd*y/HLy+FPY(kx,ky,kz)) &
                                                       * cos(twicePi*kzd*z/HLz+FPZ(kx,ky,kz))

             !          write(*,*) 'DEBUG F1, F2, F3',kx,ky,kz,kappa,xT,f1,f2,f3
             !          write(*,*) 'DEBUG x,y,z',x,y,z,kxd,kyd,kzd
                         endif
                      endif
                   enddo
                enddo
             enddo


             do kz = 1, zstep - 1
                kzd = dfloat(kz)
                do ky = mode_start, nymodes
                   kyd = dfloat(ky)
                   do kx = mode_start, nxmodes
                      kxd = dfloat(kx)
                      kappa = sqrt( (kxd*kxd)/(Lx*Lx) + (kyd*kyd)/(Ly*Ly) + (kzd*kzd)/(Lz*Lz) )
                      if (kappa.le.kappaMax) then
                       write(*,*) 'COUCOU'
                         xT = cos(FTX(kx,ky,kz)*infl_time+TAT(kx,ky,kz))
                         if (div_free_force.eq.1) then
                             f1 = f1 + xT * ( FAZ(kx,ky,kz)*twicePi*(kyd/HLy) &
                                          * sin(twicePi*kxd*x/HLx+FPZX(kx,ky,kz)) &
                                          * cos(twicePi*kyd*y/HLy+FPZY(kx,ky,kz)) &
                                          * sin(twicePi*kzd*z/HLz+FPZZ(kx,ky,kz)) &
                                         -    FAY(kx,ky,kz)*twicePi*(kzd/HLz) &
                                          * sin(twicePi*kxd*x/HLx+FPYX(kx,ky,kz)) &
                                          * sin(twicePi*kyd*y/HLy+FPYY(kx,ky,kz)) &
                                          * cos(twicePi*kzd*z/HLz+FPYZ(kx,ky,kz)) )

                             f2 = f2 + xT * ( FAX(kx,ky,kz)*twicePi*(kzd/HLz) &
                                          * sin(twicePi*kxd*x/HLx+FPXX(kx,ky,kz)) &
                                          * sin(twicePi*kyd*y/HLy+FPXY(kx,ky,kz)) &
                                          * cos(twicePi*kzd*z/HLz+FPXZ(kx,ky,kz)) &
                                         -    FAZ(kx,ky,kz)*twicePi*(kxd/HLx) &
                                          * cos(twicePi*kxd*x/HLx+FPZX(kx,ky,kz)) &
                                          * sin(twicePi*kyd*y/HLy+FPZY(kx,ky,kz)) &
                                          * sin(twicePi*kzd*z/HLz+FPZZ(kx,ky,kz)) )

                             f3 = f3 + xT * ( FAY(kx,ky,kz)*twicePi*(kxd/HLx) &
                                          * cos(twicePi*kxd*x/HLx+FPYX(kx,ky,kz)) &
                                          * sin(twicePi*kyd*y/HLy+FPYY(kx,ky,kz)) &
                                          * sin(twicePi*kzd*z/HLz+FPYZ(kx,ky,kz)) &
                                         -    FAX(kx,ky,kz)*twicePi*(kyd/HLy) & 
                                          * sin(twicePi*kxd*x/HLx+FPXX(kx,ky,kz)) &
                                          * cos(twicePi*kyd*y/HLy+FPXY(kx,ky,kz)) &
                                          * sin(twicePi*kzd*z/HLz+FPXZ(kx,ky,kz)) )

                         else
                             f1 = f1 + xT*FAX(kx,ky,kz)* cos(twicePi*kxd*x/Lx+FPX(kx,ky,kz)) &
                                                       * sin(twicePi*kyd*y/Ly+FPY(kx,ky,kz)) &
                                                       * sin(twicePi*kzd*z/Lz+FPZ(kx,ky,kz))
                             f2 = f2 + xT*FAY(kx,ky,kz)* sin(twicePi*kxd*x/Lx+FPX(kx,ky,kz)) &
                                                       * cos(twicePi*kyd*y/Ly+FPY(kx,ky,kz)) &
                                                       * sin(twicePi*kzd*z/Lz+FPZ(kx,ky,kz))
                             f3 = f3 + xT*FAZ(kx,ky,kz)* sin(twicePi*kxd*x/Lx+FPX(kx,ky,kz)) &
                                                       * sin(twicePi*kyd*y/Ly+FPY(kx,ky,kz)) &
                                                       * cos(twicePi*kzd*z/Lz+FPZ(kx,ky,kz))
!write(*,*) 'DEBUG F1, F2, F3',kx,ky,kz,kappa,xT,f1,f2,f3
                          endif
                      endif
                   enddo
                enddo
             enddo

             src(i,j,k,UMX) = f1*new_state(i,j,k,URHO)
             src(i,j,k,UMY) = f2*new_state(i,j,k,URHO)
             src(i,j,k,UMZ) = f3*new_state(i,j,k,URHO)

!write(*,*) "DEBUG",i,j,k,src(i,j,k,UMX),src(i,j,k,UMY),src(i,j,k,UMZ)

          end do
       end do
    end do



  end subroutine pc_forcing_src

end module forcing_src_module
