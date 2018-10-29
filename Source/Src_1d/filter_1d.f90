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

    integer, intent(in) :: lo  (1), hi  (1)
    integer, intent(in) :: Qlo (1), Qhi (1)
    integer, intent(in) :: Qhlo(1), Qhhi(1)
    integer, intent(in) :: ng ! (2*ng+1) filter points
    integer, intent(in) :: nstart
    integer, intent(in) :: ncnt
    integer, intent(in) :: ncomp

    double precision, intent (in   ) :: Q ( Qlo  (1):Qhi  (1), ncomp)
    double precision, intent (inout) :: Qh( Qhlo (1):Qhhi (1), ncomp)
    double precision, intent (in   ) :: w(-ng:ng)

    integer :: i, nc, l

    do nc=nstart, nstart+ncnt-1

       do l = -ng, ng
          do i=lo(1),hi(1)
             Qh(i,nc) = Qh(i,nc) + w(l) * Q(i+l,nc)
          end do
       end do

    end do

  end subroutine filter

end module filter_module
