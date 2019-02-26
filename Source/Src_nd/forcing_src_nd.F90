  ! Compute the linear forcing term for turbulence (Rosales PoF 2015)
module forcing_src_module

  implicit none

  private

  public :: pc_forcing_src
  double precision, save, public :: u0, v0, w0, forcing

contains

  subroutine pc_forcing_src(lo,hi, &
                            old_state,os_lo,os_hi, &
                            new_state,ns_lo,ns_hi, &
                            src,src_lo,src_hi,problo,dx,xlo,xhi, &
                            time,dt) bind(C, name = "pc_forcing_src")

    use amrex_constants_module, only: ZERO
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ
    use prob_params_module, only: dim

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

    src(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = ZERO
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! Linear forcing term
             src(i,j,k,UMX) = forcing * new_state(i,j,k,URHO) * (new_state(i,j,k,UMX) - u0)
             if (dim .ge. 2) then
                src(i,j,k,UMY) = forcing * new_state(i,j,k,URHO) * (new_state(i,j,k,UMY) - v0)
                if (dim .ge. 3) then
                   src(i,j,k,UMZ) = forcing * new_state(i,j,k,URHO) * (new_state(i,j,k,UMZ) - w0)
                endif
             endif

          end do
       end do
    end do

  end subroutine pc_forcing_src

end module forcing_src_module
