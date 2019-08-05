  ! Compute the MMS source for the conservative equations.
module mms_src_module

  implicit none

  private

  public :: pc_mms_src

contains

  subroutine pc_mms_src(lo,hi,&
                        S,S_lo,S_hi,&
                        src,src_lo,src_hi,&
                        problo,dx,time) bind(C, name = "pc_mms_src")


    use meth_params_module, only : NVAR

#ifdef USE_MASA
    use network, only : nspecies
    use amrex_constants_module, only: HALF
    use meth_params_module, only : URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS
    use prob_params_module, only: dim
    use masa
#endif

    implicit none

    integer          :: lo(3),hi(3)
    integer          :: S_lo(3),S_hi(3),src_lo(3),src_hi(3)
    double precision :: S(S_lo(1):S_hi(1),S_lo(2):S_hi(2),S_lo(3):S_hi(3),NVAR)
    double precision :: src(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),NVAR)
    double precision :: problo(3),dx(3),time

#ifdef USE_MASA

    ! local
    integer :: i,j,k
    double precision :: x,y,z

    do k = lo(3), hi(3)
       z = problo(3) + dx(3)*(dble(k) + HALF)

       do j = lo(2), hi(2)
          y = problo(2) + dx(2)*(dble(j) + HALF)

          do i = lo(1), hi(1)
             x = problo(1) + dx(1)*(dble(i) + HALF)

             ! Evaluate the source using MASA
             src(i,j,k,URHO) = masa_eval_3d_source_rho(x,y,z)
             src(i,j,k,UMX) = masa_eval_3d_source_rho_u(x,y,z)
             if (dim .ge. 2) then
                src(i,j,k,UMY) = masa_eval_3d_source_rho_v(x,y,z)
                if (dim .ge. 3) then
                   src(i,j,k,UMZ) = masa_eval_3d_source_rho_w(x,y,z)
                endif
             endif
             src(i,j,k,UEDEN) = masa_eval_3d_source_rho_e(x,y,z)
             src(i,j,k,UFS:UFS+nspecies-1) = src(i,j,k,URHO) * S(i,j,k,UFS:UFS+nspecies-1) / masa_eval_3d_exact_rho(x,y,z)

          end do
       end do
    end do

#else
    call bl_error('MASA is not turned on. Turn on with USE_MASA=TRUE.')
#endif

  end subroutine pc_mms_src

end module mms_src_module
