! Filtering functions
module filter_module

  implicit none

  public :: filter

contains

  ! ------------------------------------------------------------------------------
  !> Multidimensional filter kernel given filtering weights
  ! ------------------------------------------------------------------------------
  subroutine filter(lo,  hi,&
                    Q,   Qlo,   Qhi,&
                    Qh,  Qhlo,  Qhhi,&
                    ng,&
                    w,&
                    nstart,&
                    ncnt,&
                    ncomp) bind(C, name = "filter")

    use amrex_constants_module

    implicit none

    integer, intent(in) :: lo  (3), hi  (3)
    integer, intent(in) :: Qlo (3), Qhi (3)
    integer, intent(in) :: Qhlo(3), Qhhi(3)
    integer, intent(in) :: ng ! (2*ng+1) filter points
    integer, intent(in) :: nstart
    integer, intent(in) :: ncnt
    integer, intent(in) :: ncomp

    double precision, intent (in   ) :: Q (  Qlo  (1):Qhi  (1),   Qlo  (2): Qhi  (2), Qlo  (3): Qhi  (3), ncomp)
    double precision, intent (inout) :: Qh(  Qhlo (1):Qhhi (1),   Qhlo (2): Qhhi (2), Qhlo (3): Qhhi (3), ncomp)
    double precision, intent (in   ) :: w(-ng:ng)

    integer :: i, j, k, nc, l, m, n

    do nc=nstart, nstart+ncnt-1

       do n = -ng, ng
          do k=lo(3),hi(3)
             do m = -ng, ng
                do j=lo(2),hi(2)
                   do l = -ng, ng
                      do i=lo(1),hi(1)
                         Qh(i,j,k,nc) = Qh(i,j,k,nc) + w(l) * w(m) * w(n) * Q(i+l,j+m,k+n,nc)
                      end do
                   end do
                end do
             end do
          end do
       end do

    end do

  end subroutine filter

end module filter_module
