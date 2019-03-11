module weno_module

 use amrex_fort_module, only : amrex_real

  implicit none

  private

  public weno5js_face, weno5z_face
      
  real(amrex_real), parameter :: b1=13.d0/12.d0, oneSixth=1.d0/6.d0, oneHalf=1.d0/2.d0
  real(amrex_real), dimension(0:2), parameter :: weno5_face_wghts_1 = &
      (/ 0.3d0, 0.6d0, 0.1d0 /)

  integer         , save :: wenop = 1
  real(amrex_real), save :: eps = 1.d-40
  
contains

  ! :::
  ! ::: ----------------------------------------------------------------
  ! :::  Routines below are for WENO interpolation
  ! ::: ----------------------------------------------------------------
  ! :::
  !
  ! WENO Variants
  !     0: WENO-JS 5th order
  !     1: WENO-Z  5th order
        
 subroutine weno5z_face(v, vl, vr)
    real(amrex_real), intent(in)  :: v(-2:3)
    real(amrex_real), intent(out) :: vl , vr   ! left and right at i+1/2

    real(amrex_real) :: vr_2, vr_1, vr_0, vl_2, vl_1, vl_0
    real(amrex_real) :: beta(0:2)
    real(amrex_real) :: alpha(0:2)
    real(amrex_real) :: alpha1
    real(amrex_real) :: tau

    beta(2) = b1*(v(-2)-2.d0*v(-1)+v(0))**2 + 0.25d0*(v(-2)-4.d0*v(-1)+3.d0*v(0))**2
    beta(1) = b1*(v(-1)-2.d0*v( 0)+v(1))**2 + 0.25d0*(v(-1)-v(1))**2
    beta(0) = b1*(v( 0)-2.d0*v( 1)+v(2))**2 + 0.25d0*(3.d0*v(0)-4.d0*v(1)+v(2))**2

    tau = abs(beta(2) - beta(0))

    beta(2) = 1 + (tau / (eps+beta(2)))**wenop
    beta(1) = 1 + (tau / (eps+beta(1)))**wenop
    beta(0) = 1 + (tau / (eps+beta(0)))**wenop

    alpha(2) =      beta(2)
    alpha(1) = 6.d0*beta(1)
    alpha(0) = 3.d0*beta(0)
    alpha1 = 1.d0/(alpha(2) + alpha(1) + alpha(0))

    vl_2 = 2.d0*v(-2) - 7.d0*v(-1) + 11.d0*v(0)
    vl_1 =     -v(-1) + 5.d0*v( 0) +  2.d0*v(1)
    vl_0 = 2.d0*v( 0) + 5.d0*v( 1) -       v(2)

    vl = oneSixth*alpha1*(alpha(2)*vl_2 + alpha(1)*vl_1 + alpha(0)*vl_0)

    !-----------------------------------------------------------

    beta(2) = b1*(v(3)-2.d0*v( 2)+v( 1))**2 + 0.25d0*(v(3)-4.d0*v(2)+3.d0*v(1))**2
    beta(1) = b1*(v(2)-2.d0*v( 1)+v( 0))**2 + 0.25d0*(v(2)-v(0))**2
    beta(0) = b1*(v(1)-2.d0*v( 0)+v(-1))**2 + 0.25d0*(3.d0*v(1)-4.d0*v(0)+v(-1))**2

    tau = abs(beta(2) - beta(0))

    beta(2) = 1 + (tau / (eps+beta(2)))**wenop
    beta(1) = 1 + (tau / (eps+beta(1)))**wenop
    beta(0) = 1 + (tau / (eps+beta(0)))**wenop

    alpha(2) =      beta(2)
    alpha(1) = 6.d0*beta(1)
    alpha(0) = 3.d0*beta(0)
    alpha1 = 1.d0/(alpha(2) + alpha(1) + alpha(0))

    vr_2 = 11.d0*v( 1) - 7.d0*v( 2) +  2.d0*v(3)
    vr_1 = vl_0
    vr_0 = vl_1

    vr = oneSixth*alpha1*(alpha(2)*vr_2 + alpha(1)*vr_1 + alpha(0)*vr_0)

    return
  end subroutine weno5z_face

  ! :::
  ! ::: ----------------------------------------------------------------
  ! ::: ----------------------------------------------------------------
  ! :::
  
  subroutine weno5js_face(v, vl, vr)
    
    real(amrex_real), intent(in)  :: v(-2:3)
    real(amrex_real), intent(out) :: vl , vr   ! left and right at i+1/2

    real(amrex_real) :: vr_2, vr_1, vr_0, vl_2, vl_1, vl_0
    real(amrex_real) :: beta(0:2)
    real(amrex_real) :: alpha(0:2)
    real(amrex_real) :: alpha1
    
    beta(2) = b1*(v(-2)-2.d0*v(-1)+v(0))**2 + 0.25d0*(v(-2)-4.d0*v(-1)+3.d0*v(0))**2
    beta(1) = b1*(v(-1)-2.d0*v( 0)+v(1))**2 + 0.25d0*(v(-1)-v(1))**2
    beta(0) = b1*(v( 0)-2.d0*v( 1)+v(2))**2 + 0.25d0*(3.d0*v(0)-4.d0*v(1)+v(2))**2

    beta(2) = 1.d0/(eps+beta(2))**wenop
    beta(1) = 1.d0/(eps+beta(1))**wenop
    beta(0) = 1.d0/(eps+beta(0))**wenop

    alpha(2) =      beta(2)
    alpha(1) = 6.d0*beta(1)
    alpha(0) = 3.d0*beta(0)
    alpha1 = 1.d0/(alpha(2) + alpha(1) + alpha(0))
    
    vl_2 = 2.d0*v(-2) - 7.d0*v(-1) + 11.d0*v(0)
    vl_1 =     -v(-1) + 5.d0*v( 0) +  2.d0*v(1)
    vl_0 = 2.d0*v( 0) + 5.d0*v( 1) -       v(2)
    
    vl = oneSixth*alpha1*(alpha(2)*vl_2 + alpha(1)*vl_1 + alpha(0)*vl_0)

    !-----------------------------------------------------------

    beta(2) = b1*(v(3)-2.d0*v( 2)+v( 1))**2 + 0.25d0*(v(3)-4.d0*v(2)+3.d0*v(1))**2
    beta(1) = b1*(v(2)-2.d0*v( 1)+v( 0))**2 + 0.25d0*(v(2)-v(0))**2
    beta(0) = b1*(v(1)-2.d0*v( 0)+v(-1))**2 + 0.25d0*(3.d0*v(1)-4.d0*v(0)+v(-1))**2

    beta(2) = 1.d0/(eps+beta(2))**wenop
    beta(1) = 1.d0/(eps+beta(1))**wenop
    beta(0) = 1.d0/(eps+beta(0))**wenop

    alpha(2) =      beta(2)
    alpha(1) = 6.d0*beta(1)
    alpha(0) = 3.d0*beta(0)
    alpha1 = 1.d0/(alpha(2) + alpha(1) + alpha(0))
    
    vr_2 = 11.d0*v( 1) - 7.d0*v( 2) +  2.d0*v(3)
    vr_1 = vl_0
    vr_0 = vl_1
    
    vr = oneSixth*alpha1*(alpha(2)*vr_2 + alpha(1)*vr_1 + alpha(0)*vr_0)

    return
  end subroutine weno5js_face
  
  
  ! :::
  ! ::: ----------------------------------------------------------------
  ! ::: ----------------------------------------------------------------
  ! :::
  
  ! Normalizes ``alphas'' into weights.
  subroutine getWghts(alpha, len)
      integer         , intent(in)    :: len
      real(amrex_real), intent(inout) :: alpha(1:len)

      real(amrex_real) :: alpha1
      integer          :: i

      alpha1 = 1.d0 / sum(alpha)
      alpha  = alpha1 * alpha
  end subroutine
  
end module weno_module
