module weno_module

 use amrex_fort_module, only : amrex_real

  implicit none

  private

  public weno5js_face, weno5z_face, weno3z_face, weno7z_face
      
  real(amrex_real), parameter :: b1=13.d0/12.d0, oneSixth=1.d0/6.d0, oneHalf=1.d0/2.d0, oneTwelve=1.0d0/12.0d0
  real(amrex_real), dimension(0:2), parameter :: weno5_face_wghts_1 = &
      (/ 0.3d0, 0.6d0, 0.1d0 /)
  real(amrex_real), dimension(0:3), parameter :: weno7_face_wghts_1 = &
      (/ 1.0d0/35.0d0, 12.0d0/35.0d0, 18.0d0/35.0d0, 4.0d0/35.0d0 /)

  integer         , save :: wenop = 2
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
  !     2: WENO-Z  3th order
  !     3: WENO-Z  7th order
        
 subroutine weno5z_face(v, vl, vr)
    implicit none
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
  
  subroutine weno7z_face(v, vl, vr)
    implicit none
    real(amrex_real), intent(in)  :: v(-3:4)
    real(amrex_real), intent(out) :: vl , vr   ! left and right at i+1/2

    real(amrex_real) :: vr_3, vr_2, vr_1, vr_0, vl_3, vl_2, vl_1, vl_0
    real(amrex_real) :: beta(0:3)
    real(amrex_real) :: alpha(0:3)
    real(amrex_real) :: alpha1
    real(amrex_real) :: tau

    beta(3) = v(-3)*(547.0d0*v(-3) - 3882.0d0*v(-2) + 4642.0d0*v(-1) - 1854.0d0*v(0)) + &
              v(-2)*(7043.0d0*v(-2) - 17246.0d0*v(-1) + 7042.0d0*v(0)) + &
              v(-1)*(11003.0d0*v(-1) - 9402.0d0*v(0)) + 2107.0d0*v(0)**2
    
    beta(2) = v(-2)*(267.0d0*v(-2) - 1642.0d0*v(-1) + 1602.0d0*v(0) - 494.0d0*v(1)) + &
              v(-1)*(2843.0d0*v(-1) - 5966.0d0*v(0) + 1922.0d0*v(1)) + &
              v( 0)*(3443.0d0*v(0) - 2522.0d0*v(1)) + 547.0d0*v(1)**2

    beta(1) = v(-1)*(547.0d0*v(-1) - 2522.0d0*v(0) + 1922.0d0*v(1) - 494.0d0*v(2)) + &
              v( 0)*(3443.0d0*v(0) - 5966.0d0*v(1) + 1602.0d0*v(2)) + &
              v( 1)*(2843.0d0*v(1) - 1642.0d0*v(2)) + 267.0d0*v(2)**2
              
    beta(0) = v( 0)*(2107.0d0*v(0) - 9402.0d0*v(1) + 7042.0d0*v(2) - 1854.0d0*v(3)) + &
              v( 1)*(11003.0d0*v(1) - 17246.0d0*v(2) + 4642.0d0*v(3)) + &
              v( 2)*(7043.0d0*v(2) - 3882.0d0*v(3)) + 547.0d0*v(3)**2          

    tau = abs(beta(3) - beta(0))

    beta(3) = 1 + (tau / (eps+beta(3)))**wenop
    beta(2) = 1 + (tau / (eps+beta(2)))**wenop
    beta(1) = 1 + (tau / (eps+beta(1)))**wenop
    beta(0) = 1 + (tau / (eps+beta(0)))**wenop

    alpha(3) = weno7_face_wghts_1(0) * beta(3)
    alpha(2) = weno7_face_wghts_1(1) * beta(2)
    alpha(1) = weno7_face_wghts_1(2) * beta(1)
    alpha(0) = weno7_face_wghts_1(3) * beta(0)
    alpha1 = 1.d0/(alpha(3) + alpha(2) + alpha(1) + alpha(0))

    vl_3 = -3.0d0*v(-3) + 13.0d0*v(-2) - 23.0d0*v(-1) + 25.0d0*v(0)
    vl_2 = 1.0d0*v(-2) - 5.0d0*v(-1) + 13.0d0*v(0) + 3.0d0*v(1)
    vl_1 = -1.0d0*v(-1) + 7.0d0*v(0) + 7.0d0*v(1) - 1.0d0*v(2)
    vl_0 = 3.0d0*v(0) + 13.0d0*v(1) - 5.0d0*v(2) + 1.0d0*v(3)

    vl = oneTwelve*alpha1*(alpha(3)*vl_3 + alpha(2)*vl_2 + alpha(1)*vl_1 + alpha(0)*vl_0)

    !-----------------------------------------------------------

    beta(3) = v(4)*(547.0d0*v(4) - 3882.0d0*v(3) + 4642.0d0*v(2) - 1854.0d0*v(1)) + &
              v(3)*(7043.0d0*v(3) - 17246.0d0*v(2) + 7042.0d0*v(1)) + &
              v(2)*(11003.0d0*v(2) - 9402.0d0*v(1)) + 2107.0d0*v(1)**2
    
    beta(2) = v(3)*(267.0d0*v(3) - 1642.0d0*v(2) + 1602.0d0*v(1) - 494.0d0*v(0)) + &
              v(2)*(2843.0d0*v(2) - 5966.0d0*v(1) + 1922.0d0*v(0)) + &
              v( 1)*(3443.0d0*v(1) - 2522.0d0*v(0)) + 547.0d0*v(0)**2

    beta(1) = v(2)*(547.0d0*v(2) - 2522.0d0*v(1) + 1922.0d0*v(0) - 494.0d0*v(-1)) + &
              v( 1)*(3443.0d0*v(1) - 5966.0d0*v(0) + 1602.0d0*v(-1)) + &
              v( 0)*(2843.0d0*v(0) - 1642.0d0*v(-1)) + 267.0d0*v(-1)**2
              
    beta(0) = v( 1)*(2107.0d0*v(1) - 9402.0d0*v(0) + 7042.0d0*v(-1) - 1854.0d0*v(-2)) + &
              v( 0)*(11003.0d0*v(0) - 17246.0d0*v(-1) + 4642.0d0*v(-2)) + &
              v( -1)*(7043.0d0*v(-1) - 3882.0d0*v(-2)) + 547.0d0*v(-2)**2

    tau = abs(beta(3) - beta(0))

    beta(3) = 1 + (tau / (eps+beta(3)))**wenop
    beta(2) = 1 + (tau / (eps+beta(2)))**wenop
    beta(1) = 1 + (tau / (eps+beta(1)))**wenop
    beta(0) = 1 + (tau / (eps+beta(0)))**wenop

    alpha(3) = weno7_face_wghts_1(0) * beta(3)
    alpha(2) = weno7_face_wghts_1(1) * beta(2)
    alpha(1) = weno7_face_wghts_1(2) * beta(1)
    alpha(0) = weno7_face_wghts_1(3) * beta(0)
    alpha1 = 1.d0/(alpha(3) + alpha(2) + alpha(1) + alpha(0))

    vr_3 = 25.0d0*v(1) - 23.0d0*v(2) + 13.0d0*v(3) - 3.0d0*v(4)
    vr_2 = vl_0
    vr_1 = vl_1
    vr_0 = vl_2

    vr = oneTwelve*alpha1*(alpha(3)*vr_3 + alpha(2)*vr_2 + alpha(1)*vr_1 + alpha(0)*vr_0)

    return
  end subroutine weno7z_face

  ! :::
  ! ::: ----------------------------------------------------------------
  ! ::: ----------------------------------------------------------------
  ! :::
  
  subroutine weno5js_face(v, vl, vr)
    implicit none
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
   
  
    subroutine weno3z_face(v, vl, vr)
      implicit none
      real(amrex_real), intent(in)  :: v(-2:3)
      real(amrex_real), intent(out) :: vl, vr  ! left and right at i + 1/2

      real(amrex_real) :: vr_1, vr_0, vl_1, vl_0
      real(amrex_real) :: beta(0:1)
      real(amrex_real) :: alpha(0:1)
      real(amrex_real) :: alpha1
      real(amrex_real) :: tau

      beta(1) = (v(-1) - v( 0))**2
      beta(0) = (v( 0) - v( 1))**2

      tau = abs(beta(1) - beta(0))

      beta(1) = 1 + (tau / (eps+beta(1)))**wenop
      beta(0) = 1 + (tau / (eps+beta(0)))**wenop

      alpha(1) = 2.d0*beta(1)
      alpha(0) = 1.d0*beta(0)
      alpha1 = 1.d0/(alpha(1) + alpha(0))

      vl_1 = -1.d0*v(-1) +  3.d0*v( 0)
      vl_0 =  1.d0*v( 0) +  1.d0*v( 1)

      vl   = oneHalf*alpha1 * (alpha(1)*vl_1 + alpha(0)*vl_0)

      beta(1) = (v( 2) - v( 1))**2
      beta(0) = (v( 1) - v( 0))**2

      tau = abs(beta(1) - beta(0))

      beta(1) = 1 + (tau / (eps+beta(1)))**wenop
      beta(0) = 1 + (tau / (eps+beta(0)))**wenop

      alpha(1) = 2.d0*beta(1)
      alpha(0) = 1.d0*beta(0)
      alpha1 = 1.d0/(alpha(1) + alpha(0))

      vr_1 =  3.d0*v( 1) -  1.d0*v( 2)
      vr_0 = vl_0

      vr   = oneHalf*alpha1 * (alpha(1)*vr_1 + alpha(0)*vr_0)
      
  end subroutine weno3z_face

  
  ! :::
  ! ::: ----------------------------------------------------------------
  ! ::: ----------------------------------------------------------------
  ! :::
  
  ! Normalizes ``alphas'' into weights.
  subroutine getWghts(alpha, len)
      implicit none
      integer         , intent(in)    :: len
      real(amrex_real), intent(inout) :: alpha(1:len)

      real(amrex_real) :: alpha1
      integer          :: i

      alpha1 = 1.d0 / sum(alpha)
      alpha  = alpha1 * alpha
  end subroutine
  
end module weno_module
