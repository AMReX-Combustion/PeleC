module diffusion_module

  implicit none

  public

contains

  subroutine pc_diffextrap(lo, hi, dif, d_lo, d_hi, nc) &
       bind(C, name="pc_diffextrap")

    ! this routine extrapolates the diffusion term into the
    ! ghostcells

    use prob_params_module, only: dg

    implicit none

    integer :: lo(3), hi(3)
    integer :: d_lo(3), d_hi(3)
    integer :: nc
    double precision :: dif(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)

    ! Local variables

    integer :: i, j, k, n

    do n=1,nc
       ! left side
       if (d_lo(1) .lt. lo(1)) then
          i = lo(1)-1*dg(1)
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                dif(i,j,k,n) = dif(i+1*dg(1),j,k,n)
             end do
          end do
       endif

       ! right side
       if (d_hi(1) .gt. hi(1)) then
          i = hi(1)+1*dg(1)
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                dif(i,j,k,n) = dif(i-1*dg(1),j,k,n)
             end do
          end do
       endif

       ! bottom side
       if (d_lo(2) .lt. lo(2)) then
          j = lo(2)-1*dg(2)
          do k = lo(3), hi(3)
             do i = lo(1), hi(1)
                dif(i,j,k,n) = dif(i,j+1*dg(2),k,n)
             end do
          end do
       endif

       ! top side
       if (d_hi(2) .gt. hi(2)) then
          j = hi(2)+1*dg(2)
          do k = lo(3), hi(3)
             do i = lo(1), hi(1)
                dif(i,j,k,n) = dif(i,j-1*dg(2),k,n)
             end do
          end do
       endif

       ! down side
       if (d_lo(3) .lt. lo(3)) then
          k = lo(3)-1*dg(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                dif(i,j,k,n) = dif(i,j,k+1*dg(3),n)
             end do
          end do
       endif

       ! up side
       if (d_hi(3) .gt. hi(3)) then
          k = hi(3)+1*dg(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                dif(i,j,k,n) = dif(i,j,k-1*dg(3),n)
             end do
          end do
       endif

       ! k-edges
       if (d_lo(1) .lt. lo(1) .and. d_lo(2) .lt. lo(2)) then
          i = lo(1)-1*dg(1)
          j = lo(2)-1*dg(2)
          do k = lo(3), hi(3)
             dif(i,j,k,n) = dif(i+1*dg(1),j+1*dg(2),k,n)
          end do
       endif

       if (d_lo(1) .lt. lo(1) .and. d_hi(2) .gt. hi(2)) then
          i = lo(1)-1*dg(1)
          j = hi(2)+1*dg(2)
          do k = lo(3), hi(3)
             dif(i,j,k,n) = dif(i+1*dg(1),j-1*dg(2),k,n)
          end do
       endif

       if (d_hi(1) .gt. hi(1) .and. d_lo(2) .lt. lo(2)) then
          i = hi(1)+1*dg(1)
          j = lo(2)-1*dg(2)
          do k = lo(3), hi(3)
             dif(i,j,k,n) = dif(i-1*dg(1),j+1*dg(2),k,n)
          end do
       endif

       if (d_hi(1) .gt. hi(1) .and. d_hi(2) .gt. hi(2)) then
          i = hi(1)+1*dg(1)
          j = hi(2)+1*dg(2)
          do k = lo(3), hi(3)
             dif(i,j,k,n) = dif(i-1*dg(1),j-1*dg(2),k,n)
          end do
       endif

       ! j-edges
       if (d_lo(1) .lt. lo(1) .and. d_lo(3) .lt. lo(3)) then
          i = lo(1)-1*dg(1)
          k = lo(3)-1*dg(3)
          do j = lo(2), hi(2)
             dif(i,j,k,n) = dif(i+1*dg(1),j,k+1*dg(3),n)
          end do
       endif

       if (d_lo(1) .lt. lo(1) .and. d_hi(3) .gt. hi(3)) then
          i = lo(1)-1*dg(1)
          k = hi(3)+1*dg(3)
          do j = lo(2), hi(2)
             dif(i,j,k,n) = dif(i+1*dg(1),j,k-1*dg(3),n)
          end do
       endif

       if (d_hi(1) .gt. hi(1) .and. d_lo(3) .lt. lo(3)) then
          i = hi(1)+1*dg(1)
          k = lo(3)-1*dg(3)
          do j = lo(2), hi(2)
             dif(i,j,k,n) = dif(i-1*dg(1),j,k+1*dg(3),n)
          end do
       endif

       if (d_hi(1) .gt. lo(1) .and. d_hi(3) .gt. lo(3)) then
          i = hi(1)+1*dg(1)
          k = hi(3)+1*dg(3)
          do j = lo(2), hi(2)
             dif(i,j,k,n) = dif(i-1*dg(1),j,k-1*dg(3),n)
          end do
       endif

       ! i-edges
       if (d_lo(2) .lt. lo(2) .and. d_lo(3) .lt. lo(3)) then
          j = lo(2)-1*dg(2)
          k = lo(3)-1*dg(3)
          do i = lo(1), hi(1)
             dif(i,j,k,n) = dif(i,j+1*dg(2),k+1*dg(3),n)
          end do
       endif

       if (d_lo(2) .lt. lo(2) .and. d_hi(3) .gt. hi(3)) then
          j = lo(2)-1*dg(2)
          k = hi(3)+1*dg(3)
          do i = lo(1), hi(1)
             dif(i,j,k,n) = dif(i,j+1*dg(2),k-1*dg(3),n)
          end do
       endif

       if (d_hi(2) .gt. hi(2) .and. d_lo(3) .lt. lo(3)) then
          j = hi(2)+1*dg(2)
          k = lo(3)-1*dg(3)
          do i = lo(1), hi(1)
             dif(i,j,k,n) = dif(i,j-1*dg(2),k+1*dg(3),n)
          end do
       endif

       if (d_hi(2) .gt. hi(2) .and. d_hi(3) .gt. hi(3)) then
          j = hi(2)+1*dg(2)
          k = hi(3)+1*dg(3)
          do i = lo(1), hi(1)
             dif(i,j,k,n) = dif(i,j-1*dg(2),k-1*dg(3),n)
          end do
       endif

       ! corners
       if (d_lo(1) .lt. lo(1) .and. d_lo(2) .lt. lo(2) .and. d_lo(3) .lt. lo(3)) then
          i = lo(1)-1*dg(1)
          j = lo(2)-1*dg(2)
          k = lo(3)-1*dg(3)
          dif(i,j,k,n) = dif(i+1*dg(1),j+1*dg(2),k+1*dg(3),n)
       endif

       if (d_lo(1) .lt. lo(1) .and. d_hi(2) .gt. hi(2) .and. d_lo(3) .lt. lo(3)) then
          i = lo(1)-1*dg(1)
          j = hi(2)+1*dg(2)
          k = lo(3)-1*dg(3)
          dif(i,j,k,n) = dif(i+1*dg(1),j-1*dg(2),k+1*dg(3),n)
       endif

       if (d_hi(1) .gt. hi(1) .and. d_lo(2) .lt. lo(2) .and. d_lo(3) .lt. lo(3)) then
          i = hi(1)+1*dg(1)
          j = hi(2)+1*dg(2)
          k = lo(3)-1*dg(3)
          dif(i,j,k,n) = dif(i-1*dg(1),j-1*dg(2),k+1*dg(3),n)
       endif

       if (d_hi(1) .gt. hi(1) .and. d_hi(2) .gt. hi(2) .and. d_lo(3) .lt. lo(3)) then
          i = hi(1)+1*dg(1)
          j = lo(2)-1*dg(2)
          k = lo(3)-1*dg(3)
          dif(i,j,k,n) = dif(i-1*dg(1),j+1*dg(2),k+1*dg(3),n)
       endif

       if (d_lo(1) .lt. lo(1) .and. d_lo(2) .lt. lo(2) .and. d_hi(3) .gt. hi(3)) then
          i = lo(1)-1*dg(1)
          j = lo(2)-1*dg(2)
          k = hi(3)+1*dg(3)
          dif(i,j,k,n) = dif(i+1*dg(1),j+1*dg(2),k-1*dg(3),n)
       endif

       if (d_lo(1) .lt. lo(1) .and. d_hi(2) .gt. hi(2) .and. d_hi(3) .gt. hi(3)) then
          i = lo(1)-1*dg(1)
          j = hi(2)+1*dg(2)
          k = hi(3)+1*dg(3)
          dif(i,j,k,n) = dif(i+1*dg(1),j-1*dg(2),k-1*dg(3),n)
       endif

       if (d_hi(1) .gt. hi(1) .and. d_lo(2) .lt. lo(2) .and. d_hi(3) .gt. hi(3)) then
          i = hi(1)+1*dg(1)
          j = lo(2)-1*dg(2)
          k = hi(3)+1*dg(3)
          dif(i,j,k,n) = dif(i-1*dg(1),j+1*dg(2),k-1*dg(3),n)
       endif

       if (d_hi(1) .gt. hi(1) .and. d_hi(2) .gt. hi(2) .and. d_hi(3) .gt. hi(3)) then
          i = hi(1)+1*dg(1)
          j = hi(2)+1*dg(2)
          k = hi(3)+1*dg(3)
          dif(i,j,k,n) = dif(i-1*dg(1),j-1*dg(2),k-1*dg(3),n)
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

end module diffusion_module
