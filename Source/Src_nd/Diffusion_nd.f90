module diffusion_module

  implicit none

  public

contains

  subroutine pc_diffextrap(lo, hi, vlo, vhi, dif, d_lo, d_hi, nc) &
       bind(C, name="pc_diffextrap")

    ! this routine extrapolates the diffusion term into the
    ! ghostcells

    implicit none

    integer :: lo(3), hi(3)
    integer :: vlo(3), vhi(3)
    integer :: d_lo(3), d_hi(3)
    integer :: nc
    double precision :: dif(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)

    ! Local variables

    integer :: i, j, k, n

    do n=1,nc
       ! left side
       if (lo(1) .eq. vlo(1)) then
          do i = d_lo(1), lo(1)-1
             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   dif(i,j,k,n) = dif(lo(1),j,k,n)
                end do
             end do
          end do
       endif

       ! right side
       if (hi(1) .eq. vhi(1)) then
          do i = hi(1)+1, d_hi(1)
             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   dif(i,j,k,n) = dif(hi(1),j,k,n)
                end do
             end do
          end do
       endif

       ! bottom side
       if (lo(2) .eq. vlo(2)) then
          do j = d_lo(2), lo(2)-1
             do k = lo(3), hi(3)
                do i = lo(1), hi(1)
                   dif(i,j,k,n) = dif(i,lo(2),k,n)
                end do
             end do
          end do
       endif

       ! top side
       if (hi(2) .eq. vhi(2)) then
          do j = hi(2)+1, d_hi(2)
             do k = lo(3), hi(3)
                do i = lo(1), hi(1)
                   dif(i,j,k,n) = dif(i,hi(2),k,n)
                end do
             end do
          end do
       endif

       ! down side
       if (lo(3) .eq. vlo(3)) then
          do k = d_lo(3), lo(3)-1
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   dif(i,j,k,n) = dif(i,j,lo(3),n)
                end do
             end do
          end do
       endif

       ! up side
       if (hi(3) .eq. vhi(3)) then
          do k = hi(3)+1, d_hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   dif(i,j,k,n) = dif(i,j,hi(3),n)
                end do
             end do
          end do
       endif

       ! k-edges
       if (lo(1) .eq. vlo(1) .and. lo(2) .eq. vlo(2)) then
          do j = d_lo(2), lo(2)-1
             do i = d_lo(1), lo(1)-1
                do k = lo(3), hi(3)
                   dif(i,j,k,n) = dif(lo(1),lo(2),k,n)
                end do
             end do
          end do
       endif

       if (lo(1) .eq. vlo(1) .and. hi(2) .eq. vhi(2)) then
          do j = hi(2)+1, d_hi(2)
             do i = d_lo(1), d_lo(1)-1
                do k = lo(3), hi(3)
                   dif(i,j,k,n) = dif(lo(1),hi(2),k,n)
                end do
             end do
          end do
       endif

       if (hi(1) .eq. vhi(1) .and. lo(2) .eq. vlo(2)) then
          do j = d_lo(2), lo(2)-1
             do i = hi(1)+1, d_hi(1)
                do k = lo(3), hi(3)
                   dif(i,j,k,n) = dif(hi(1),lo(2),k,n)
                end do
             end do
          end do
       endif

       if (hi(1) .eq. vhi(1) .and. hi(2) .eq. vhi(2)) then
          do j = hi(2)+1, d_hi(2)
             do i = hi(1)+1, d_hi(1)
                do k = lo(3), hi(3)
                   dif(i,j,k,n) = dif(hi(1),hi(2),k,n)
                end do
             end do
          end do
       endif

       ! j-edges
       if (lo(1) .eq. vlo(1) .and. lo(3) .eq. vlo(3)) then
          do k = d_lo(3), lo(3)-1
             do i = d_lo(1), lo(1)-1
                do j = lo(2), hi(2)
                   dif(i,j,k,n) = dif(lo(1),j,lo(3),n)
                end do
             end do
          end do
       endif

       if (lo(1) .eq. vlo(1) .and. hi(3) .eq. vhi(3)) then
          do k = hi(3)+1, d_hi(3)
             do i = d_lo(1), lo(1)-1
                do j = lo(2), hi(2)
                   dif(i,j,k,n) = dif(lo(1),j,hi(3),n)
                end do
             end do
          end do
       endif

       if (hi(1) .eq. vhi(1) .and. lo(3) .eq. vlo(3)) then
          do k = d_lo(3), lo(3)-1
             do i = hi(1)+1, d_hi(1)
                do j = lo(2), hi(2)
                   dif(i,j,k,n) = dif(hi(1),j,lo(3),n)
                end do
             end do
          end do
       endif

       if (hi(1) .eq. vhi(1) .and. hi(3) .eq. vhi(3)) then
          do k = hi(3)+1, d_hi(3)
             do i = hi(1)+1, d_hi(1)
                do j = lo(2), hi(2)
                   dif(i,j,k,n) = dif(hi(1),j,hi(3),n)
                end do
             end do
          end do
       endif

       ! i-edges
       if (lo(2) .eq. vlo(2) .and. lo(3) .eq. vlo(3)) then
          do k = d_lo(3), lo(3)-1
             do j = d_lo(2), lo(2)-1
                do i = lo(1), hi(1)
                   dif(i,j,k,n) = dif(i,lo(2),lo(3),n)
                end do
             end do
          end do
       endif

       if (lo(2) .eq. vlo(2) .and. hi(3) .eq. vhi(3)) then
          do k = hi(3)+1, d_hi(3)
             do j = d_lo(2), lo(2)-1
                do i = lo(1), hi(1)
                   dif(i,j,k,n) = dif(i,lo(2),hi(3),n)
                end do
             end do
          end do
       endif

       if (hi(2) .eq. vhi(2) .and. lo(3) .eq. vlo(3)) then
          do k = d_lo(3), lo(3)-1
             do j = hi(2)+1, d_hi(2)
                do i = lo(1), hi(1)
                   dif(i,j,k,n) = dif(i,hi(2),lo(3),n)
                end do
             end do
          end do
       endif

       if (hi(2) .eq. vhi(2) .and. hi(3) .eq. vhi(3)) then
          do k = hi(3)+1, d_hi(3)
             do j = hi(2)+1, d_hi(2)
                do i = lo(1), hi(1)
                   dif(i,j,k,n) = dif(i,hi(2),hi(3),n)
                end do
             end do
          end do
       endif

       ! corners
       if (lo(1) .eq. vlo(1) .and. lo(2) .eq. vlo(2) .and. lo(3) .eq. vlo(3)) then
          do k = d_lo(3), lo(3)-1
             do j = d_lo(2), lo(2)-1
                do i = d_lo(1), lo(1)-1
                   dif(i,j,k,n) = dif(lo(1),lo(2),lo(3),n)
                end do
             end do
          end do
       endif

       if (lo(1) .eq. vlo(1) .and. hi(2) .eq. vhi(2) .and. lo(3) .eq. vlo(3)) then
          do k = d_lo(3), lo(3)-1
             do j = hi(2)+1, d_hi(2)
                do i = d_lo(1), lo(1)-1
                   dif(i,j,k,n) = dif(lo(1),hi(2),lo(3),n)
                end do
             end do
          end do
       endif

       if (hi(1) .eq. vhi(1) .and. lo(2) .eq. vlo(2) .and. lo(3) .eq. vlo(3)) then
          do k = d_lo(3), lo(3)-1
             do j = d_lo(2), lo(2)-1
                do i = hi(1)+1, d_hi(1)
                   dif(i,j,k,n) = dif(hi(1),lo(2),lo(3),n)
                end do
             end do
          end do
       endif

       if (hi(1) .eq. vhi(1) .and. hi(2) .eq. vhi(2) .and. lo(3) .eq. vlo(3)) then
          do k = d_lo(3), lo(3)-1
             do j = hi(2)+1, d_hi(2)
                do i = hi(1)+1, d_hi(1)
                   dif(i,j,k,n) = dif(hi(1),hi(2),lo(3),n)
                end do
             end do
          end do
       endif
       
       if (lo(1) .eq. vlo(1) .and. lo(2) .eq. vlo(2) .and. hi(3) .eq. vhi(3)) then
          do k = hi(3)+1, d_hi(3)
             do j = d_lo(2), lo(2)-1
                do i = d_lo(1), lo(1)-1
                   dif(i,j,k,n) = dif(lo(1),lo(2),hi(3),n)
                end do
             end do
          end do
       endif

       if (lo(1) .eq. vlo(1) .and. hi(2) .eq. vhi(2) .and. hi(3) .eq. vhi(3)) then
          do k = d_hi(3), hi(3)-1
             do j = hi(2)+1, d_hi(2)
                do i = d_lo(1), lo(1)-1
                   dif(i,j,k,n) = dif(lo(1),hi(2),hi(3),n)
                end do
             end do
          end do
       endif

       if (hi(1) .eq. vhi(1) .and. lo(2) .eq. vlo(2) .and. hi(3) .eq. vhi(3)) then
          do k = d_hi(3), hi(3)-1
             do j = d_lo(2), lo(2)-1
                do i = hi(1)+1, d_hi(1)
                   dif(i,j,k,n) = dif(hi(1),lo(2),hi(3),n)
                end do
             end do
          end do
       endif

       if (hi(1) .eq. vhi(1) .and. hi(2) .eq. vhi(2) .and. hi(3) .eq. vhi(3)) then
          do k = d_hi(3), hi(3)-1
             do j = hi(2)+1, d_hi(2)
                do i = hi(1)+1, d_hi(1)
                   dif(i,j,k,n) = dif(hi(1),hi(2),hi(3),n)
                end do
             end do
          end do
       endif

    enddo

  end subroutine pc_diffextrap

  subroutine pc_move_transport_coeffs_to_ec(lo,hi,dlo,dhi, &
       cfab,c_lo,c_hi, &
       efab,e_lo,e_hi, dir, nc, do_harmonic) &
       bind(C, name="pc_move_transport_coeffs_to_ec")

    use prob_params_module, only : physbc_lo, physbc_hi
    use amrex_constants_module

    implicit none

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: dlo(3), dhi(3)
    integer         , intent(in   ) :: c_lo(3), c_hi(3)
    integer         , intent(in   ) :: e_lo(3), e_hi(3)
    integer         , intent(in   ) :: dir, nc, do_harmonic
    real (amrex_real), intent(in   ) :: cfab(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3),nc)
    real (amrex_real), intent(inout) :: efab(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),nc)

    ! local variables
    integer          :: i, j, k, n

    if (do_harmonic .eq. 0) then
       if (dir .EQ. 0) then
          do n = 1,nc
             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = lo(1), hi(1)+1
                      efab(i,j,k,n) = half*(cfab(i,j,k,n) + cfab(i-1,j,k,n))
                   end do
                end do
             end do
          end do
       else if (dir .EQ. 1) then
          do n = 1,nc
             do k = lo(3), hi(3)
                do j = lo(2), hi(2)+1
                   do i = lo(1), hi(1)
                      efab(i,j,k,n) = half*(cfab(i,j,k,n) + cfab(i,j-1,k,n))
                   end do
                end do
             end do
          end do
       else if (dir .EQ. 2) then
          do n = 1,nc
             do k = lo(3), hi(3)+1
                do j = lo(2), hi(2)
                   do i = lo(1), hi(1)
                      efab(i,j,k,n) = half*(cfab(i,j,k,n) + cfab(i,j,k-1,n))
                   end do
                end do
             end do
          end do
       end if
    else
       if (dir .EQ. 0) then
          do n = 1,nc
             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = lo(1), hi(1)+1
                      if ((cfab(i,j,k,n) * cfab(i-1,j,k,n)) .gt.zero) then
                         efab(i,j,k,n) =&
                              2*(cfab(i,j,k,n) * cfab(i-1,j,k,n))&
                              /(cfab(i,j,k,n) + cfab(i-1,j,k,n))
                      else
                         efab(i,j,k,n) = zero
                      endif
                   end do
                end do
             end do
          end do
       else if (dir .EQ. 1) then
          do n = 1,nc
             do k = lo(3), hi(3)
                do j = lo(2), hi(2)+1
                   do i = lo(1), hi(1)
                      if((cfab(i,j,k,n) * cfab(i,j-1,k,n)).gt.zero) then
                         efab(i,j,k,n) =&
                              2*(cfab(i,j,k,n) * cfab(i,j-1,k,n))&
                              /(cfab(i,j,k,n) + cfab(i,j-1,k,n))
                      else
                         efab(i,j,k,n) = zero
                      endif
                   end do
                end do
             end do
          end do
       else if (dir .EQ. 2) then
          do n = 1,nc
             do k = lo(3), hi(3)+1
                do j = lo(2), hi(2)
                   do i = lo(1), hi(1)
                      if((cfab(i,j,k,n) * cfab(i,j,k-1,n)).gt.zero) then
                         efab(i,j,k,n) =&
                              2*(cfab(i,j,k,n) * cfab(i,j,k-1,n))&
                              /(cfab(i,j,k,n) + cfab(i,j,k-1,n))
                      else
                         efab(i,j,k,n) = zero
                      endif
                   end do
                end do
             end do
          end do
       end if
    end if

  end subroutine pc_move_transport_coeffs_to_ec


    ! One function for all directions
    subroutine pc_diffterm_aux(lo,  hi,&
                         dmnlo, dmnhi,&
                         Q,   Qlo,   Qhi,&
                         Daux,Dauxlo,Dauxhi,&
                         A,  Alo,  Ahi,&
                         f,  flo,  fhi,&
                         V,   Vlo,   Vhi,&
                         D,   Dlo,   Dhi,&
                         deltax, dir) bind(C, name = "pc_diffterm_aux")

    use network, only : naux
    use meth_params_module, only : NVAR, UFX, QVAR, QFX
    use amrex_constants_module
    use prob_params_module, only : physbc_lo, physbc_hi

    implicit none

    integer, intent(in) ::     lo(3),    hi(3)
    integer, intent(in) ::  dmnlo(3), dmnhi(3)
    integer, intent(in) ::    Qlo(3),   Qhi(3)
    integer, intent(in) :: Dauxlo(3),Dauxhi(3)
    integer, intent(in) ::   Alo(3),  Ahi(3)
    integer, intent(in) ::   flo(3),  fhi(3)
    integer, intent(in) ::    Dlo(3),   Dhi(3)
    integer, intent(in) ::    Vlo(3),   Vhi(3)
    integer, intent(in) ::   dir
    
    double precision, intent(in   ) ::    Q(   Qlo(1):   Qhi(1),   Qlo(2):   Qhi(2),   Qlo(3):   Qhi(3), QVAR)
    double precision, intent(in   ) :: Daux(Dauxlo(1):Dauxhi(1),Dauxlo(2):Dauxhi(2),Dauxlo(3):Dauxhi(3), naux)
    double precision, intent(in   ) ::   A(  Alo(1):  Ahi(1),  Alo(2):  Ahi(2),  Alo(3):  Ahi(3)  )
    double precision, intent(inout) ::   f(  flo(1):  fhi(1),  flo(2):  fhi(2),  flo(3):  fhi(3), NVAR)
    double precision, intent(inout) ::    D(   Dlo(1):   Dhi(1),   Dlo(2):   Dhi(2),   Dlo(3):   Dhi(3), NVAR)
    double precision, intent(in   ) ::    V(   Vlo(1):   Vhi(1),   Vlo(2):   Vhi(2),   Vlo(3):   Vhi(3)  )
    double precision, intent(in   ) :: deltax(3)

    integer :: i, j, k, n, direction
    integer, dimension(3) :: isdir
    double precision :: dAddir
    double precision :: gfac(3)
    double precision :: dxinv(3)

    direction = dir + 1 ! Convert to x=1, y=2, z=3
    isdir = 0
    isdir(direction) = 1
    
    dxinv = 1.d0/deltax
    gfac = dxinv 
    
    ! calculate flux
    do n=1,naux
       do k=lo(3),hi(3)+isdir(3)
          do j=lo(2),hi(2)+isdir(2)
             do i=lo(1),hi(1)+isdir(1)
                dAddir = gfac(direction) * (Q(i,j,k,QFX+n-1) - Q(i-isdir(1),j-isdir(2),k-isdir(3),QFX+n-1))
                f(i,j,k,UFX+n-1) = - dAddir * Daux(i,j,k,n)
             end do
          enddo
       enddo
    enddo
    
    ! calculate mass flow rate across surface
    do n=UFX,UFX+naux-1
       do k=lo(3),hi(3)+isdir(3)
          do j=lo(2),hi(2)+isdir(2)
             do i=lo(1),hi(1)+isdir(1)
                f(i,j,k,n) = f(i,j,k,n) * A(i,j,k)
             enddo
          enddo
       enddo
    enddo

    ! Calculate net flow
    do n=UFX,UFX+naux-1
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                D(i,j,k,n) = D(i,j,k,n) - (f(i+isdir(1),j+isdir(2),k+isdir(3),n)-f(i,j,k,n) )/V(i,j,k)
             end do
          end do
       end do
    end do

  end subroutine pc_diffterm_aux

end module diffusion_module
