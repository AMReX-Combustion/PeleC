module weno_module

 use amrex_fort_module, only : amrex_real

  implicit none

  private

  public weno5js_face, weno5z_face
      
  real(amrex_real), parameter :: b1=13.d0/12.d0, oneSixth=1.d0/6.d0, oneHalf=1.d0/2.d0
  real(amrex_real), dimension(0:2), parameter :: weno5_face_wghts_1 = &
      (/ 0.3d0, 0.6d0, 0.1d0 /)
  real(amrex_real), dimension(0:2), parameter :: weno5_face_wghts_teno = &
      (/ 0.6d0, 0.3d0, 0.1d0 /)
  real(amrex_real), dimension(0:3), parameter :: weno6_face_wghts_1 = &
      (/ 9.0d0/20.d0, 6.0d0/20.0d0, 1.0d0/20.0d0, 4.0d0/20.0d0 /)
  real(amrex_real), dimension(0:4), parameter :: weno7_face_wghts_1 = &
      (/ 18.d0/35.d0, 9.d0/35.d0, 3.d0/35.d0, 4.d0/35.d0, 1.d0/35.d0 /)
  !real(amrex_real), dimension(0:5), parameter :: weno8_face_wghts_1 = &
  !    (/ 30.d0/70.d0, 18.0d0/70.0d0, 4.0d0/70.0d0, 12.0d0/70.0d0, 1.0d0/70.0d0, 5.0d0/70.0d0 /)
  !    
  real(amrex_real), dimension(0:5), parameter :: weno8_face_wghts_1 = &
      (/ 0.4130855804023061d0, 0.2193140179474727d0, 0.06459052081827904d0, &
      0.1439236310125986d0, 0.02746675592227204d0, 0.1316194638970745d0 /)
  
  real(amrex_real), parameter, private :: dsp = 0.0463783d0
  real(amrex_real), parameter, private :: dss = 0.01d0
  real(amrex_real), dimension(0:3), parameter, private :: Ck = &
       (/ 1.5d0*dsp+1.5d0*dss, 0.5d0-1.5d0*dsp+4.5d0*dss, 0.5d0-1.5d0*dsp-4.5d0*dss, 1.5d0*dsp-1.5d0*dss /)

      
  integer         , save :: wenop = 1
  real(amrex_real), save :: eps = 1.d-40
  ! WENO Variants
  !     0: WENO-JS 5th order
  !     1: WENO-Z  5th order

  real(amrex_real), parameter :: wenot_cutoff = 1.0d-7

contains

  ! :::
  ! ::: ----------------------------------------------------------------
  ! :::  Routines below are for WENO interpolation
  ! ::: ----------------------------------------------------------------
  ! :::
        
    subroutine weno_mdcd_face(v, vl, vr)
    real(amrex_real), intent(in) :: v(-2:3)
    real(amrex_real), intent(out) :: vl, vr

    real(amrex_real) :: vk(0:3), wk(0:3), IS(0:3)

    IS(0) = eps + b1*(v(-2)-2.d0*v(-1)+v(0))**2 + 0.25d0*(v(-2)-4.d0*v(-1)+3.d0*v(0))**2
    IS(1) = eps + b1*(v(-1)-2.d0*v( 0)+v(1))**2 + 0.25d0*(v(-1)-v(1))**2
    IS(2) = eps + b1*(v( 0)-2.d0*v( 1)+v(2))**2 + 0.25d0*(3.d0*v(0)-4.d0*v(1)+v(2))**2
    IS(3) = eps + b1*(v( 1)-2.d0*v( 2)+v(3))**2 + 0.25d0*(5.d0*v(1)-8.d0*v(2)+3.d0*v(3))**2
    IS(3) = maxval(IS)

    wk = Ck / IS**wenop

    vk(0) = ( 2.d0*v(-2) - 7.d0*v(-1) + 11.d0*v(0))*oneSixth
    vk(1) = (     -v(-1) + 5.d0*v( 0) +  2.d0*v(1))*oneSixth
    vk(2) = ( 2.d0*v( 0) + 5.d0*v( 1) -       v(2))*oneSixth
    vk(3) = (11.d0*v( 1) - 7.d0*v( 2) +  2.d0*v(3))*oneSixth

    vl = sum(wk*vk)/sum(wk)

    !-----------------------------------------------------------

    IS(0) = eps + b1*(v(3)-2.d0*v( 2)+v( 1))**2 + 0.25d0*(v(3)-4.d0*v(2)+3.d0*v(1))**2
    IS(1) = eps + b1*(v(2)-2.d0*v( 1)+v( 0))**2 + 0.25d0*(v(2)-v(0))**2
    IS(2) = eps + b1*(v(1)-2.d0*v( 0)+v(-1))**2 + 0.25d0*(3.d0*v(1)-4.d0*v(0)+v(-1))**2
    IS(3) = eps + b1*(v(0)-2.d0*v(-1)+v(-2))**2 + 0.25d0*(5.d0*v(0)-8.d0*v(-1)+3.d0*v(-2))**2
    IS(3) = maxval(IS)

    wk = Ck / IS**wenop

    vk = vk(3:0:-1)

    vr = sum(wk*vk)/sum(wk)

  end subroutine weno_mdcd_face


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
  
  subroutine weno5m_face(v, vl, vr)
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

    call getWghts(alpha, 3)
    call henrick_func(alpha, weno5_face_wghts_1, 3)
    call getWghts(alpha, 3)
    
    vl_2 = 2.d0*v(-2) - 7.d0*v(-1) + 11.d0*v(0)
    vl_1 =     -v(-1) + 5.d0*v( 0) +  2.d0*v(1)
    vl_0 = 2.d0*v( 0) + 5.d0*v( 1) -       v(2)
    
    vl = oneSixth*(alpha(2)*vl_2 + alpha(1)*vl_1 + alpha(0)*vl_0)

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

    call getWghts(alpha, 3)
    call henrick_func(alpha, weno5_face_wghts_1, 3)
    call getWghts(alpha, 3)
    
    vr_2 = 11.d0*v( 1) - 7.d0*v( 2) +  2.d0*v(3)
    vr_1 = vl_0
    vr_0 = vl_1
    
    vr = oneSixth*(alpha(2)*vr_2 + alpha(1)*vr_1 + alpha(0)*vr_0)

    return
  end subroutine weno5m_face

  ! :::
  ! ::: ----------------------------------------------------------------
  ! ::: ----------------------------------------------------------------
  ! :::
    
  subroutine weno5teno_face(v, vl, vr)
          
    real(amrex_real), intent(in)  :: v(-2:3)
    real(amrex_real), intent(out) :: vl , vr   ! left and right at i+1/2

    real(amrex_real) :: vr_2, vr_1, vr_0, vl_2, vl_1, vl_0
    real(amrex_real) :: beta(0:2)
    real(amrex_real) :: alpha(0:2)
    real(amrex_real) :: tau
    real(amrex_real) :: dpi_l, dpip1_l, dpim1_l, dpmag, rmag, dpi_r, dpip1_r, dpim1_r
    real(amrex_real) :: Ct_param, delta, r_cut_off, amp_min, amp_max

    beta(:)=0.0d0
    beta(2) = b1*(v(-2)-2.d0*v(-1)+v(0))**2.0d0 + 0.25d0*(v(-2)-4.d0*v(-1)+3.d0*v(0))**2.0d0
    beta(1) = b1*(v(-1)-2.d0*v( 0)+v(1))**2.0d0 + 0.25d0*(v(-1)-v(1))**2.0d0
    beta(0) = b1*(v( 0)-2.d0*v( 1)+v(2))**2.0d0 + 0.25d0*(3.d0*v(0)-4.d0*v(1)+v(2))**2.0d0

    tau = abs(beta(0) - beta(2))
    
    alpha(2) = (1.0d0 + tau / (beta(2) + eps))**6
    alpha(1) = (1.0d0 + tau / (beta(1) + eps))**6
    alpha(0) = (1.0d0 + tau / (beta(0) + eps))**6
    
    dpi_l = 0.25d0*(-v(1)+2.0d0*v(0)-v(-1))
    dpip1_l = 0.25d0*(-v(2)+2.0d0*v(1)-v(0))
    dpim1_l = 0.25d0*(-v(0)+2.0d0*v(-1)-v(-2))
    dpmag = 0.5d0*((dpi_l-dpip1_l)**2.0d0 + (dpi_l-dpim1_l)**2.0d0)
    rmag = (dpmag/(v(0)**2.0d0))+eps
    if (dpmag == 0.0d0) rmag=eps
    !write(*,*) 'DEBUG 0.10',rmag,dpmag,eps
    !write(*,*) 'DEBUG 0',v
    
    amp_max = 7.5d0
    amp_min = 2.5d0
    r_cut_off = 1.0d-5
    delta = 1.0d-5
    Ct_param = amp_min+(0.5d0*(amp_max-amp_min)*(1.0d0+tanh((r_cut_off-rmag)/delta)))
    Ct_param = 10.0d0**(-floor(Ct_param))
     


    !Ct_param = wenot_cutoff
    
    call getWghts(alpha, 3)
    call thresh(alpha, weno5_face_wghts_teno, Ct_param, 3)
    call getWghts(alpha, 3)
    
    vl_2 =  2.d0*v(-2) - 7.d0*v(-1) + 11.d0*v(0)
    vl_1 =  2.d0*v( 0) + 5.d0*v( 1) - 1.0d0*v(2)
    vl_0 = -1.d0*v(-1) + 5.d0*v( 0) + 2.0d0*v(1)
    
    ! Below was for WENO, not TENO 
    !vl_1 =     -v(-1) + 5.d0*v( 0) +  2.d0*v(1)
    !vl_0 = 2.d0*v( 0) + 5.d0*v( 1) -       v(2)
    
    vl = oneSixth*(alpha(2)*vl_2 + alpha(1)*vl_1 + alpha(0)*vl_0)

    !-----------------------------------------------------------

    beta(2) = b1*(v(3)-2.d0*v( 2)+v( 1))**2 + 0.25d0*(v(3)-4.d0*v(2)+3.d0*v(1))**2
    beta(1) = b1*(v(2)-2.d0*v( 1)+v( 0))**2 + 0.25d0*(v(2)-v(0))**2
    beta(0) = b1*(v(1)-2.d0*v( 0)+v(-1))**2 + 0.25d0*(3.d0*v(1)-4.d0*v(0)+v(-1))**2

    tau = abs(beta(0) - beta(2))

    alpha(2) = (1.0d0 + tau / (beta(2) + eps))**6
    alpha(1) = (1.0d0 + tau / (beta(1) + eps))**6
    alpha(0) = (1.0d0 + tau / (beta(0) + eps))**6
    
    
    dpi_r = dpip1_l !(-v(0)+2.0d0*v(1)-v(2))/4.0d0
    dpip1_r = dpi_l !(-v(-1)+2.0d0*v(0)-v(1))/4.0d0
    dpim1_r = 0.25d0*(-v(1)+2.0d0*v(2)-v(3))
    dpmag = 0.5d0*((dpi_r-dpip1_r)**2.0d0 + (dpi_r-dpim1_r)**2.0d0)
    rmag = (dpmag/(v(1)**2.0d0))+eps
    if (dpmag == 0.0d0) rmag=eps
    
    

    Ct_param = amp_min+(0.5d0*(amp_max-amp_min)*(1.0d0+tanh((r_cut_off-rmag)/delta)))
    Ct_param = 10.0d0**(-floor(Ct_param))
        !write(*,*) 'DEBUG 2',Ct_param,rmag
    !write(*,*) ''   

    !Ct_param = wenot_cutoff
    
    call getWghts(alpha, 3)
    call thresh(alpha, weno5_face_wghts_teno, Ct_param, 3)
    call getWghts(alpha, 3)
    
    vr_2 = 11.d0*v( 1) - 7.d0*v( 2) +  2.d0*v(3)
    vr_1 = vl_0
    vr_0 = vl_1
    
    vr = oneSixth*(alpha(2)*vr_2 + alpha(1)*vr_1 + alpha(0)*vr_0)

    return
  end subroutine weno5teno_face

  ! :::
  ! ::: ----------------------------------------------------------------
  ! ::: ----------------------------------------------------------------
  ! :::
  
  subroutine weno7teno_face(v, vl, vr)

    real(amrex_real), intent(in)  :: v(-3:4)
    real(amrex_real), intent(out) :: vl , vr   ! left and right at i+1/2

    real(amrex_real) :: vr_4, vr_3, vr_2, vr_1, vr_0, vl_4, vl_3, vl_2, vl_1, vl_0
    real(amrex_real) :: beta(0:4)
    real(amrex_real) :: alpha(0:4)
    real(amrex_real) :: tau
    real(amrex_real) :: coeff_1_36, coeff_13_12, coeff_781_720
    
    real(amrex_real) :: weno7t_cutoff

    coeff_1_36 = 1.0d0/36.0d0
    coeff_13_12 = 13.0d0/12.0d0
    coeff_781_720 = 781.0d0/720.0d0
    
    beta(4) =  coeff_1_36*(11.0d0*v(0) - 18.0d0*v(-1) + 9.d0*v(-2) - 2.d0*v(-3))**2 &
             + coeff_13_12*(2.d0*v(0) - 5.d0*v(-1) + 4.d0*v(-2) - v(-3))**2  &
             + coeff_781_720*(v(0) - 3.d0*v(-1) + 3.d0*v(-2) - v(-3))**2
    beta(3) =  coeff_1_36*(-11.0d0*v(0) + 18.0d0*v(1) - 9.d0*v(2) + 2.d0*v(3))**2 &
             + coeff_13_12*(2.d0*v(0) - 5.d0*v(1) + 4.d0*v(2) - v(3))**2  &
             + coeff_781_720*(-v(0) + 3.d0*v(1) - 3.d0*v(2) + v(3))**2
    beta(2) = b1*(v(-2)-2.d0*v(-1)+v(0))**2 + 0.25d0*(v(-2)-4.d0*v(-1)+3.d0*v(0))**2
    beta(1) = b1*(v(-1)-2.d0*v( 0)+v(1))**2 + 0.25d0*(v(-1)-v(1))**2
    beta(0) = b1*(v( 0)-2.d0*v( 1)+v(2))**2 + 0.25d0*(3.d0*v(0)-4.d0*v(1)+v(2))**2

    tau = abs(beta(4) - 1.d0/6.d0*(beta(1)+beta(2)+4.0d0*(beta(0))))

    alpha(4) = (1.0d0 + tau / (beta(4) + eps))**6
    alpha(3) = (1.0d0 + tau / (beta(3) + eps))**6
    alpha(2) = (1.0d0 + tau / (beta(2) + eps))**6
    alpha(1) = (1.0d0 + tau / (beta(1) + eps))**6
    alpha(0) = (1.0d0 + tau / (beta(0) + eps))**6

    call compute_Ct_cutoff(v(-3:3),weno7t_cutoff,1)

    call getWghts(alpha, 5)
    call thresh(alpha, weno7_face_wghts_1, weno7t_cutoff, 5)
    call getWghts(alpha, 5)

    vl_4 =  25.d0*v(0) - 23.d0*v(-1) + 13.d0*v(-2) - 3.d0*v(-3)
    vl_3 =  3.d0*v(0) + 13.d0*v(1) - 5.d0*v(2) + v(3) 
    vl_2 =  2.d0*v(-2) - 7.d0*v(-1) + 11.d0*v(0)
    vl_1 =  2.d0*v( 0) + 5.d0*v( 1) - 1.0d0*v(2)
    vl_0 = -1.d0*v(-1) + 5.d0*v( 0) + 2.0d0*v(1)

    vl =  oneSixth*(alpha(2)*vl_2 + alpha(1)*vl_1 + alpha(0)*vl_0) &
        + (1.0d0/12.0d0)*(alpha(3)*vl_3 + alpha(4)*vl_4)

    !-----------------------------------------------------------

    beta(4) =  coeff_1_36*(11.0d0*v(1) - 18.0d0*v(2) + 9.d0*v(3) - 2.d0*v(4))**2 &
             + coeff_13_12*(2.d0*v(1) - 5.d0*v(2) + 4.d0*v(3) - v(4))**2  &
             + coeff_781_720*(v(1) - 3.d0*v(2) + 3.d0*v(3) - v(4))**2
    beta(3) =  coeff_1_36*(-11.0d0*v(1) + 18.0d0*v(0) - 9.d0*v(-1) + 2.d0*v(-2))**2 &
             + coeff_13_12*(2.d0*v(1) - 5.d0*v(0) + 4.d0*v(-1) - v(-2))**2  &
             + coeff_781_720*(-v(1) + 3.d0*v(0) - 3.d0*v(-1) + v(-2))**2
    beta(2) = b1*(v(3)-2.d0*v( 2)+v( 1))**2 + 0.25d0*(v(3)-4.d0*v(2)+3.d0*v(1))**2
    beta(1) = b1*(v(2)-2.d0*v( 1)+v( 0))**2 + 0.25d0*(v(2)-v(0))**2
    beta(0) = b1*(v(1)-2.d0*v( 0)+v(-1))**2 + 0.25d0*(3.d0*v(1)-4.d0*v(0)+v(-1))**2

    tau = abs(beta(4) - 1.d0/6.d0*(beta(1)+beta(2)+4.0d0*(beta(0))))

    alpha(4) = (1.0d0 + tau / (beta(4) + eps))**6
    alpha(3) = (1.0d0 + tau / (beta(3) + eps))**6
    alpha(2) = (1.0d0 + tau / (beta(2) + eps))**6
    alpha(1) = (1.0d0 + tau / (beta(1) + eps))**6
    alpha(0) = (1.0d0 + tau / (beta(0) + eps))**6

    call compute_Ct_cutoff(v(-2:4),weno7t_cutoff,2)

    call getWghts(alpha, 5)
    call thresh(alpha, weno7_face_wghts_1, weno7t_cutoff, 5)
    call getWghts(alpha, 5)

    vr_4 =  25.d0*v(1) - 23.d0*v(2) + 13.d0*v(3) - 3.d0*v(4)
    vr_3 =  3.d0*v(1) + 13.d0*v(0) - 5.d0*v(-1) + v(-2)
    vr_2 = 11.d0*v( 1) - 7.d0*v( 2) +  2.d0*v(3)
    vr_1 = vl_0
    vr_0 = vl_1

    vr = oneSixth*(alpha(2)*vr_2 + alpha(1)*vr_1 + alpha(0)*vr_0)  &
        + (1.0d0/12.0d0)*(alpha(3)*vr_3 + alpha(4)*vr_4)


    return
  end subroutine weno7teno_face

  ! :::
  ! ::: ----------------------------------------------------------------
  ! ::: ----------------------------------------------------------------
  ! :::
  
  subroutine weno8teno_face(v, vl, vr)

    real(amrex_real), intent(in)  :: v(-3:5)
    real(amrex_real), intent(out) :: vl , vr   ! left and right at i+1/2

    real(amrex_real) :: vr_5, vr_4, vr_3, vr_2, vr_1, vr_0, vl_5, vl_4, vl_3, vl_2, vl_1, vl_0
    real(amrex_real) :: beta(0:5)
    real(amrex_real) :: alpha(0:5)
    real(amrex_real) :: tau
    real(amrex_real) :: coeff_1_36, coeff_13_12, coeff_781_720
    
    real(amrex_real) :: weno8t_cutoff

    coeff_1_36 = 1.0d0/36.0d0
    coeff_13_12 = 13.0d0/12.0d0
    coeff_781_720 = 781.0d0/720.0d0
    
    beta(5) = (1.d0/144.d0)*(-25.d0*v(0) + 48.d0*v(1) - 36.d0*v(2) + 16.d0*v(3) - 3.d0*v(4))**2 &
            + (13.d0/1728.d0)*(35.d0*v(0) - 104.d0*v(1) + 114.d0*v(2) - 56.d0*v(3) + 11.d0*v(4))**2 &
            + (781.d0/2880.d0)*(-5.d0*v(0) + 18.d0*v(1) - 24.d0*v(2) + 14.d0*v(3) - 3.d0*v(4))**2 &
            - (1.d0/4320.d0)*(35.d0*v(0) - 104.d0*v(1) + 114.d0*v(2) - 56.d0*v(3) + 11.d0*v(4))* &
              (v(0) - 4.d0*v(1) + 6.d0*v(2) - 4.d0*v(3) + v(4)) &
            + (32803.d0/30240.d0)*(v(0) - 4.0d0*v(1) + 6.0d0*v(2) - 4.0d0*v(3) + v(4))**2
            
    beta(4) =  coeff_1_36*(11.0d0*v(0) - 18.0d0*v(-1) + 9.d0*v(-2) - 2.d0*v(-3))**2 &
             + coeff_13_12*(2.d0*v(0) - 5.d0*v(-1) + 4.d0*v(-2) - v(-3))**2  &
             + coeff_781_720*(v(0) - 3.d0*v(-1) + 3.d0*v(-2) - v(-3))**2
    beta(3) =  coeff_1_36*(-11.0d0*v(0) + 18.0d0*v(1) - 9.d0*v(2) + 2.d0*v(3))**2 &
             + coeff_13_12*(2.d0*v(0) - 5.d0*v(1) + 4.d0*v(2) - v(3))**2  &
             + coeff_781_720*(-v(0) + 3.d0*v(1) - 3.d0*v(2) + v(3))**2
    beta(2) = b1*(v(-2)-2.d0*v(-1)+v(0))**2 + 0.25d0*(v(-2)-4.d0*v(-1)+3.d0*v(0))**2
    beta(1) = b1*(v(-1)-2.d0*v( 0)+v(1))**2 + 0.25d0*(v(-1)-v(1))**2
    beta(0) = b1*(v( 0)-2.d0*v( 1)+v(2))**2 + 0.25d0*(3.d0*v(0)-4.d0*v(1)+v(2))**2

    tau = abs(beta(5) - 1.d0/6.d0*(beta(1)+beta(2)+4.0d0*(beta(0))))

    alpha(5) = (1.0d0 + tau / (beta(5) + eps))**6
    alpha(4) = (1.0d0 + tau / (beta(4) + eps))**6
    alpha(3) = (1.0d0 + tau / (beta(3) + eps))**6
    alpha(2) = (1.0d0 + tau / (beta(2) + eps))**6
    alpha(1) = (1.0d0 + tau / (beta(1) + eps))**6
    alpha(0) = (1.0d0 + tau / (beta(0) + eps))**6

    call compute_Ct_cutoff(v(-3:3),weno8t_cutoff,1)

    call getWghts(alpha, 6)
    call thresh(alpha, weno8_face_wghts_1, weno8t_cutoff, 6)
    call getWghts(alpha, 6)

    vl_5 =  12.d0*v(0) + 77.d0*v(1) - 43.d0*v(2) + 17.d0*v(3) - 3.d0*v(4)
    vl_4 =  25.d0*v(0) - 23.d0*v(-1) + 13.d0*v(-2) - 3.d0*v(-3)
    vl_3 =  3.d0*v(0) + 13.d0*v(1) - 5.d0*v(2) + v(3) 
    vl_2 =  2.d0*v(-2) - 7.d0*v(-1) + 11.d0*v(0)
    vl_1 =  2.d0*v( 0) + 5.d0*v( 1) - 1.0d0*v(2)
    vl_0 = -1.d0*v(-1) + 5.d0*v( 0) + 2.0d0*v(1)

    vl =  oneSixth*(alpha(2)*vl_2 + alpha(1)*vl_1 + alpha(0)*vl_0) &
        + (1.0d0/12.0d0)*(alpha(3)*vl_3 + alpha(4)*vl_4) &
        + (1.0d0/60.0d0)*alpha(5)*vl_5

    !-----------------------------------------------------------
    beta(5) = (1.d0/144.d0)*(-25.d0*v(1) + 48.d0*v(0) - 36.d0*v(-1) + 16.d0*v(-2) - 3.d0*v(-3))**2 &
            + (13.d0/1728.d0)*(35.d0*v(1) - 104.d0*v(0) + 114.d0*v(-1) - 56.d0*v(-2) + 11.d0*v(-3))**2 &
            + (781.d0/2880.d0)*(-5.d0*v(1) + 18.d0*v(0) - 24.d0*v(-1) + 14.d0*v(-2) - 3.d0*v(-3))**2 &
            - (1.d0/4320.d0)*(35.d0*v(1) - 104.d0*v(0) + 114.d0*v(-1) - 56.d0*v(-2) + 11.d0*v(-3))* &
              (v(1) - 4.d0*v(0) + 6.d0*v(-1) - 4.d0*v(-2) + v(-3)) &
            + (32803.d0/30240.d0)*(v(1) - 4.0d0*v(0) + 6.0d0*v(-1) - 4.0d0*v(-2) + v(-3))**2
    
    beta(4) =  coeff_1_36*(11.0d0*v(1) - 18.0d0*v(2) + 9.d0*v(3) - 2.d0*v(4))**2 &
             + coeff_13_12*(2.d0*v(1) - 5.d0*v(2) + 4.d0*v(3) - v(4))**2  &
             + coeff_781_720*(v(1) - 3.d0*v(2) + 3.d0*v(3) - v(4))**2
    beta(3) =  coeff_1_36*(-11.0d0*v(1) + 18.0d0*v(0) - 9.d0*v(-1) + 2.d0*v(-2))**2 &
             + coeff_13_12*(2.d0*v(1) - 5.d0*v(0) + 4.d0*v(-1) - v(-2))**2  &
             + coeff_781_720*(-v(1) + 3.d0*v(0) - 3.d0*v(-1) + v(-2))**2
    beta(2) = b1*(v(3)-2.d0*v( 2)+v( 1))**2 + 0.25d0*(v(3)-4.d0*v(2)+3.d0*v(1))**2
    beta(1) = b1*(v(2)-2.d0*v( 1)+v( 0))**2 + 0.25d0*(v(2)-v(0))**2
    beta(0) = b1*(v(1)-2.d0*v( 0)+v(-1))**2 + 0.25d0*(3.d0*v(1)-4.d0*v(0)+v(-1))**2

    tau = abs(beta(5) - 1.d0/6.d0*(beta(1)+beta(2)+4.0d0*(beta(0))))

    alpha(5) = (1.0d0 + tau / (beta(5) + eps))**6
    alpha(4) = (1.0d0 + tau / (beta(4) + eps))**6
    alpha(3) = (1.0d0 + tau / (beta(3) + eps))**6
    alpha(2) = (1.0d0 + tau / (beta(2) + eps))**6
    alpha(1) = (1.0d0 + tau / (beta(1) + eps))**6
    alpha(0) = (1.0d0 + tau / (beta(0) + eps))**6

    call compute_Ct_cutoff(v(-2:4),weno8t_cutoff,2)

    call getWghts(alpha, 6)
    call thresh(alpha, weno8_face_wghts_1, weno8t_cutoff, 6)
    call getWghts(alpha, 6)

    vr_5 =  12.d0*v(1) + 77.d0*v(0) - 43.d0*v(-1) + 17.d0*v(-2) - 3.d0*v(-3)
    vr_4 =  25.d0*v(1) - 23.d0*v(2) + 13.d0*v(3) - 3.d0*v(4)
    vr_3 =  3.d0*v(1) + 13.d0*v(0) - 5.d0*v(-1) + v(-2)
    vr_2 = 11.d0*v( 1) - 7.d0*v( 2) +  2.d0*v(3)
    vr_1 = vl_0
    vr_0 = vl_1

    vr = oneSixth*(alpha(2)*vr_2 + alpha(1)*vr_1 + alpha(0)*vr_0)  &
        + (1.0d0/12.0d0)*(alpha(3)*vr_3 + alpha(4)*vr_4) &
        + (1.0d0/60.0d0)*alpha(5)*vl_5


    return
  end subroutine weno8teno_face
  
  ! :::
  ! ::: ----------------------------------------------------------------
  ! ::: ----------------------------------------------------------------
  ! :::
  
  subroutine weno6tenoN_face(v, vl, vr)

    real(amrex_real), intent(in)  :: v(-2:3)
    real(amrex_real), intent(out) :: vl , vr   ! left and right at i+1/2

    real(amrex_real) :: beta(0:3)
    real(amrex_real) :: alpha(0:3)
    real(amrex_real) :: tau
    real(amrex_real) :: coeff_1_36, coeff_13_12, coeff_781_720
    
    real(amrex_real) :: weno7t_cutoff

    coeff_1_36 = 1.0d0/36.0d0
    coeff_13_12 = 13.0d0/12.0d0
    coeff_781_720 = 781.0d0/720.0d0
    
    beta(3) =  coeff_1_36*(-11.0d0*v(0) + 18.0d0*v(1) - 9.d0*v(2) + 2.d0*v(3))**2 &
             + coeff_13_12*(2.d0*v(0) - 5.d0*v(1) + 4.d0*v(2) - v(3))**2  &
             + coeff_781_720*(-v(0) + 3.d0*v(1) - 3.d0*v(2) + v(3))**2
    beta(2) = b1*(v(-2)-2.d0*v(-1)+v(0))**2 + 0.25d0*(v(-2)-4.d0*v(-1)+3.d0*v(0))**2
    beta(1) = b1*(v(-1)-2.d0*v( 0)+v(1))**2 + 0.25d0*(v(-1)-v(1))**2
    beta(0) = b1*(v( 0)-2.d0*v( 1)+v(2))**2 + 0.25d0*(3.d0*v(0)-4.d0*v(1)+v(2))**2

    tau = abs(beta(3) - 1.d0/6.d0*(beta(1)+beta(2)+4.0d0*(beta(0))))

    alpha(3) = (1.0d0 + tau / (beta(3) + eps))**6
    alpha(2) = (1.0d0 + tau / (beta(2) + eps))**6
    alpha(1) = (1.0d0 + tau / (beta(1) + eps))**6
    alpha(0) = (1.0d0 + tau / (beta(0) + eps))**6

    !call compute_Ct_cutoff(v(-3:3),weno7t_cutoff,1)

    call getWghts(alpha, 4)
    call thresh(alpha, weno6_face_wghts_1, wenot_cutoff, 4)
    if ((alpha(0) > 0.0d0) .and. (alpha(1) > 0.0d0) .and. (alpha(2) > 0.0d0) .and. (alpha(3) > 0.0d0)) then
      vl = 0.57d0*(1.0d0/60.0d0)*(2.0d0*v(-2) - 13.0d0*v(-1) + 47.0d0*v(0) + 27.0d0*v(1) - 3.0d0*v(2)) &
       + (1.0d0-0.57d0)*(1.0d0/60.0d0)*(-3.0d0*v(-1) + 27.0d0*v(0) + 47.0d0*v(1) - 13.0d0*v(2) + 2.0d0*v(3))
    else if ((alpha(0) > 0.0d0) .and. (alpha(1) > 0.0d0) .and. (alpha(2) == 0.0d0) .and. (alpha(3) > 0.0d0)) then
      vl = (1.0d0/60.0d0)*(-3.0d0*v(-1) + 27.0d0*v(0) + 47.0d0*v(1) - 13.0d0*v(2) + 2.0d0*v(3))
    else if ((alpha(0) == 0.0d0) .and. (alpha(1) > 0.0d0) .and. (alpha(2) == 0.0d0) .and. (alpha(3) > 0.0d0)) then
      vl = (1.0d0/12.0d0)*(3.0d0*v(0) + 13.0d0*v(1) - 5.0d0*v(2) + v(3))
    else if ((alpha(0) > 0.0d0) .and. (alpha(1) == 0.0d0) .and. (alpha(2) > 0.0d0) .and. (alpha(3) == 0.0d0)) then
      vl = (1.0d0/12.0d0)*(v(-2) - 5.0d0*v(-1) + 13.0d0*v(0) + 3.0d0*v(1))
    else if ((alpha(0) > 0.0d0) .and. (alpha(1) > 0.0d0) .and. (alpha(2) > 0.0d0) .and. (alpha(3) == 0.0d0)) then
      vl = (1.0d0/60.0d0)*(2.0d0*v(-2) - 13.0d0*v(-1) + 47.0d0*v(0) + 27.0d0*v(1) - 3.0d0*v(2))
    else if ((alpha(0) > 0.0d0) .and. (alpha(1) > 0.0d0) .and. (alpha(2) == 0.0d0) .and. (alpha(3) == 0.0d0)) then
      vl = (1.0d0/12.0d0)*(-1.0d0*v(-1) + 7.0d0*v(0) + 7.0d0*v(1) - 1.0d0*v(2))
    else if ((alpha(0) > 0.0d0) .and. (alpha(1) == 0.0d0) .and. (alpha(2) == 0.0d0) .and. (alpha(3) == 0.0d0)) then
      vl = (1.0d0/6.0d0)*(-1.0d0*v(-1) + 5.0d0*v(0) + 2.0d0*v(1))
    else if ((alpha(0) == 0.0d0) .and. (alpha(1) > 0.0d0) .and. (alpha(2) == 0.0d0) .and. (alpha(3) == 0.0d0)) then
      vl = (1.0d0/6.0d0)*(2.0d0*v(0) + 5.0d0*v(1) - 1.0d0*v(2))
    else if ((alpha(0) == 0.0d0) .and. (alpha(1) == 0.0d0) .and. (alpha(2) > 0.0d0) .and. (alpha(3) == 0.0d0)) then
      vl = (1.0d0/6.0d0)*(2.0d0*v(-2) - 7.0d0*v(-1) + 11.0d0*v(0))
    end if
    
    !-----------------------------------------------------------


    beta(3) =  coeff_1_36*(-11.0d0*v(1) + 18.0d0*v(0) - 9.d0*v(-1) + 2.d0*v(-2))**2 &
             + coeff_13_12*(2.d0*v(1) - 5.d0*v(0) + 4.d0*v(-1) - v(-2))**2  &
             + coeff_781_720*(-v(1) + 3.d0*v(0) - 3.d0*v(-1) + v(-2))**2
    beta(2) = b1*(v(3)-2.d0*v( 2)+v( 1))**2 + 0.25d0*(v(3)-4.d0*v(2)+3.d0*v(1))**2
    beta(1) = b1*(v(2)-2.d0*v( 1)+v( 0))**2 + 0.25d0*(v(2)-v(0))**2
    beta(0) = b1*(v(1)-2.d0*v( 0)+v(-1))**2 + 0.25d0*(3.d0*v(1)-4.d0*v(0)+v(-1))**2

    tau = abs(beta(3) - 1.d0/6.d0*(beta(1)+beta(2)+4.0d0*(beta(0))))

    alpha(3) = (1.0d0 + tau / (beta(3) + eps))**6
    alpha(2) = (1.0d0 + tau / (beta(2) + eps))**6
    alpha(1) = (1.0d0 + tau / (beta(1) + eps))**6
    alpha(0) = (1.0d0 + tau / (beta(0) + eps))**6

    !call compute_Ct_cutoff(v(-2:4),weno7t_cutoff,2)

    call getWghts(alpha, 4)
    call thresh(alpha, weno6_face_wghts_1, wenot_cutoff, 4)

    if ((alpha(0) > 0.0d0) .and. (alpha(1) > 0.0d0) .and. (alpha(2) > 0.0d0) .and. (alpha(3) > 0.0d0)) then
      vr = 0.57d0*(1.0d0/60.0d0)*(2.0d0*v(3) - 13.0d0*v(2) + 47.0d0*v(1) + 27.0d0*v(0) - 3.0d0*v(-1)) &
       + (1.0d0-0.57d0)*(1.0d0/60.0d0)*(-3.0d0*v(2) + 27.0d0*v(1) + 47.0d0*v(0) - 13.0d0*v(-1) + 2.0d0*v(-2))
    else if ((alpha(0) > 0.0d0) .and. (alpha(1) > 0.0d0) .and. (alpha(2) == 0.0d0) .and. (alpha(3) > 0.0d0)) then
      vr = (1.0d0/60.0d0)*(-3.0d0*v(2) + 27.0d0*v(1) + 47.0d0*v(0) - 13.0d0*v(-1) + 2.0d0*v(-2))
    else if ((alpha(0) == 0.0d0) .and. (alpha(1) > 0.0d0) .and. (alpha(2) == 0.0d0) .and. (alpha(3) > 0.0d0)) then
      vr = (1.0d0/12.0d0)*(3.0d0*v(1) + 13.0d0*v(0) - 5.0d0*v(-1) + v(-2))
    else if ((alpha(0) > 0.0d0) .and. (alpha(1) == 0.0d0) .and. (alpha(2) > 0.0d0) .and. (alpha(3) == 0.0d0)) then
      vr = (1.0d0/12.0d0)*(v(3) - 5.0d0*v(2) + 13.0d0*v(1) + 3.0d0*v(0))
    else if ((alpha(0) > 0.0d0) .and. (alpha(1) > 0.0d0) .and. (alpha(2) > 0.0d0) .and. (alpha(3) == 0.0d0)) then
      vr = (1.0d0/60.0d0)*(2.0d0*v(3) - 13.0d0*v(2) + 47.0d0*v(1) + 27.0d0*v(0) - 3.0d0*v(-1))
    else if ((alpha(0) > 0.0d0) .and. (alpha(1) > 0.0d0) .and. (alpha(2) == 0.0d0) .and. (alpha(3) == 0.0d0)) then
      vr = (1.0d0/12.0d0)*(-1.0d0*v(2) + 7.0d0*v(1) + 7.0d0*v(0) - 1.0d0*v(-1))
    else if ((alpha(0) > 0.0d0) .and. (alpha(1) == 0.0d0) .and. (alpha(2) == 0.0d0) .and. (alpha(3) == 0.0d0)) then
      vr = (1.0d0/6.0d0)*(-1.0d0*v(2) + 5.0d0*v(1) + 2.0d0*v(0))
    else if ((alpha(0) == 0.0d0) .and. (alpha(1) > 0.0d0) .and. (alpha(2) == 0.0d0) .and. (alpha(3) == 0.0d0)) then
      vr = (1.0d0/6.0d0)*(2.0d0*v(1) + 5.0d0*v(0) - 1.0d0*v(-1))
    else if ((alpha(0) == 0.0d0) .and. (alpha(1) == 0.0d0) .and. (alpha(2) > 0.0d0) .and. (alpha(3) == 0.0d0)) then
      vr = (1.0d0/6.0d0)*(2.0d0*v(3) - 7.0d0*v(2) + 11.0d0*v(1))
    end if
    



    return
  end subroutine weno6tenoN_face  
  
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

  ! Thresholding for weights, returns zero if weight is below given threshold.
  subroutine thresh(wghts, altwghts, cutoff, len)
      integer         , intent(in)     :: len
      real(amrex_real), intent(inout)  :: wghts   (1:len)
      real(amrex_real), intent(in)     :: altwghts(1:len)
      real(amrex_real), intent(in)     :: cutoff

      integer :: i

      do i = 1, len
        if (wghts(i) .LT. cutoff) then
            wghts(i) = 0.d0
        else
            wghts(i) = altwghts(i)
        endif
      end do
  end subroutine thresh
  
  ! The WENO-M(apped) function of Henrick et al, 2005
  subroutine henrick_func(w, wbar, len)
      integer, intent(in)             :: len
      real(amrex_real), intent(inout) :: w   (1:len)
      real(amrex_real), intent(in)    :: wbar(1:len)

      integer :: i

      do i = 1, len
        w(i) = w(i) * (wbar(i) + wbar(i)*wbar(i) - 3.d0 * wbar(i) * w(i) + w(i)*w(i)) / (wbar(i)*wbar(i) + w(i) * (1 - 2.d0 * wbar(i)))
      end do
  end subroutine
    
  subroutine compute_Ct_cutoff(v, cutoff, flag_lr)
      integer,          intent(in)      :: flag_lr
      real(amrex_real), intent(in)      :: v(-3:3)
      real(amrex_real), intent(out)     :: cutoff

      integer :: i, ind
      real(amrex_real) :: eta_half, m_fact, Cr, tseta, g_func, beta_param, coeff_1, coeff_2, weno7t_cutoff
      real(amrex_real) :: delta_p, delta_m, eps_eta
      real(amrex_real) :: eta(0:4)

      Cr = 0.25d0
      tseta = 1.0d-3
      coeff_1 = 10.5d0  !10
      coeff_2 = 3.5d0     !8
      eps_eta = ((0.9*Cr)/(1.0d0-0.9d0*Cr))*tseta**2
    
      do i = 0,4
        if (flag_lr == 1) then
          delta_m = v(i-2)-v(i-3)  
          delta_p = v(i-1)-v(i-2)
        else if (flag_lr == 2) then
          ind = -i
          delta_m = v(ind+2)-v(ind+3)  
          delta_p = v(ind+1)-v(ind+2)
        end if
        eta(i) = (abs(2.0d0*delta_m*delta_p) + eps_eta) &
                /(delta_m**2 + delta_p**2 + eps_eta)
      end do

      eta_half = min(eta(0),eta(1),eta(2),eta(3),eta(4))
      m_fact = 1.0d0 - min(1.0d0,(eta_half/Cr))
      g_func = (1.0d0+4.0d0*m_fact)*((1.0d0-m_fact)**4)
      beta_param = coeff_1 - coeff_2*(1.0d0-g_func)
      cutoff = 10.d0**(-floor(beta_param))


  end subroutine compute_Ct_cutoff

end module weno_module
