module probdata_module
  use eos_module
  implicit none

  double precision, save :: initP, Tstart, Tend, delta
  integer, save :: TestNum, numPoints
  character(len=128) :: OutputFile

contains

  subroutine RunEOStests

    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEINT, UEDEN, UTEMP, UFS
    use chemistry_module, only : nspecies, get_species_index,Ru
    use eos_type_module
    use transport_module

    integer :: iN2, iC2H6, iCO2
    type(eos_t) :: eos_state
    type(eos_t) :: eos_rhoT
    type(eos_t) :: eos_rhoP
    type(eos_t) :: eos_rhoE
    double precision :: tauDelta
    double precision :: dP_dT, dP_dtau, dH_dtau
    double precision :: PatT, PatTplus
    double precision :: PatTau, PatTauplus
    double precision :: HatTau, HatTauplus
    double precision :: amT, amTplus, amTminus
    double precision :: dam_dT, d2am_dT2, dam_dT_analy, d2am_dT2_analy
    double precision :: hmix,sumhyk
    double precision :: invdelta
    double precision :: dtaudY_analy(nspec),dtaudY_chain(nspec)
    double precision :: fYkplus,fYk,muK,sumFyk,muatYk
    type(trv_t)  :: TestTrans
    type(wtr_t)  :: whichTest
    character(len=128) :: FileName
    integer :: iunit,i
    double precision :: errorP, errorT, errorRHO, errorE
    double precision :: MaXerrorP_rhoT,MaXerrorP_rhoE
    double precision :: MaXerrorT_rhoP,MaXerrorT_rhoE
    double precision :: MaXerrorE_rhoP,MaXerrorE_rhoT

    call build(eos_state)
    call build(eos_rhoT)
    call build(eos_rhoP)
    call build(eos_rhoE)

    call transport_init()

    whichTest % wtr_get_xi    = .false.
    whichTest % wtr_get_mu    = .true.
    whichTest % wtr_get_lam   = .true.
    whichTest % wtr_get_Ddiag = .false.

    call trv_build(TestTrans,1)    

    iN2 = get_species_index("N2")
    iC2H6 = get_species_index("C2H6")
    iCO2 = get_species_index("CO2")

!!$      Test the species thermodynamic quantities one species at a time
!!$      rho, Cp, Cv, Speed of Sound, Internal energy and enthalpy

    print*,"=== Testing Thermodynamic and transport property calculation for single species, N2 ====="

    ! N2
    FileName = trim(OutputFile)//'_vsNIST_N2.txt'
    iunit = 10
    OPEN(UNIT=iunit,FILE=trim(FileName),FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")

    eos_state % molefrac(1:nspec) = 0.0d0
    eos_state % molefrac(iN2) = 1.0d0

    ! get mass fractions from mole fractions
    call eos_xty(eos_state) 
    eos_state % P = initP

    TestTrans % eos_state(1) % massfrac(1:nspec) = eos_state% massfrac(1:nspec)
    TestTrans % eos_state(1) % P = eos_state% P

    do i = 1,numPoints

       ! Generate these quantities over a range of temperature
       eos_state % T = Tstart + (Tend-Tstart)*dble(i-1)/(numPoints-1)

       ! call to compute rho and e
       call eos_tp(eos_state)

       ! Call to test transport - viscosity and thermal conductivity 
       TestTrans % eos_state(1) % T = eos_state% T
       TestTrans % eos_state(1) % rho = eos_state% rho

       ! call to compute mu and lambda
       call transport(whichTest, TestTrans)

       write(iunit,'(8(E20.8,4x))') eos_state%P, eos_state % T, eos_state % rho, eos_state%cp, eos_state%cv, eos_state%cs, TestTrans%mu(1),TestTrans%lam(1)

    end do
    CLOSE(UNIT=iunit)

!!$       Test mixture rules
    print*,"=========Testing mixture rules==================="

    print*,"======= Running CO2 + C2H6 mixture eq.ratio test = 1 ========="

    FileName = trim(OutputFile)//'_CO2_C2H6_Xco2_0_0.txt'
    iunit = 10
    OPEN(UNIT=iunit,FILE=trim(FileName),FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")

    eos_state % molefrac(1:nspec) = 0.0d0

    eos_state % molefrac(iC2H6) = 1.0d0
    eos_state % molefrac(iCO2)  = 0.0000d0

    ! get mass fractions from mole fractions
    call eos_xty(eos_state) 

    eos_state % T = 350.0d0

    do i = 1,numPoints

       ! Generate these quantities over a range of temperature
       eos_state % P = 1.0d0*1e6*10.0d0 + 40.0d0*1e6*10.0*dble(i-1)/(numPoints-1)

       ! call to compute rho and e
       call eos_tp(eos_state)

       write(iunit,'(6(E20.8,4x))') eos_state%P, eos_state % T, eos_state % rho, eos_state%cp, eos_state%cv, eos_state%cs

    end do
    CLOSE(UNIT=iunit)

    print*,"======= Running CO2 + C2H6 mixture eq. ratio test = 2 ========="

    FileName = trim(OutputFile)//'_CO2_C2H6_Xco2_0_251.txt'
    iunit = 10
    OPEN(UNIT=iunit,FILE=trim(FileName),FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")

    eos_state % molefrac(1:nspec) = 0.0d0

    eos_state % molefrac(iCO2)  = 0.25166d0
    eos_state % molefrac(iC2H6) = 1.0d0 - eos_state % molefrac(iCO2)

    ! get mass fractions from mole fractions
    call eos_xty(eos_state) 

    eos_state % T = 350.0d0

    do i = 1,numPoints

       ! Generate these quantities over a range of temperature
       eos_state % P = 1.0d0*1e6*10.0d0 + 40.0d0*1e6*10.0*dble(i-1)/(numPoints-1)

       ! call to compute rho and e
       call eos_tp(eos_state)

       write(iunit,'(6(E20.8,4x))') eos_state%P, eos_state % T, eos_state % rho, eos_state%cp, eos_state%cv, eos_state%cs

    end do
    CLOSE(UNIT=iunit)

    print*,"======= Running CO2 + C2H6 mixture eq. ratio test = 3 ========="

    FileName = trim(OutputFile)//'_CO2_C2H6_Xco2_0_49245.txt'
    iunit = 10
    OPEN(UNIT=iunit,FILE=trim(FileName),FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")

    eos_state % molefrac(1:nspec) = 0.0d0

    eos_state % molefrac(iCO2)  = 0.49245d0
    eos_state % molefrac(iC2H6) = 1.0d0 - eos_state % molefrac(iCO2)

    ! get mass fractions from mole fractions
    call eos_xty(eos_state) 

    eos_state % T = 350.0d0

    do i = 1,numPoints

       ! Generate these quantities over a range of temperature
       eos_state % P = 1.0d0*1e6*10.0d0 + 40.0d0*1e6*10.0*dble(i-1)/(numPoints-1)

       ! call to compute rho and e
       call eos_tp(eos_state)

       write(iunit,'(6(E20.8,4x))') eos_state%P, eos_state % T, eos_state % rho, eos_state%cp, eos_state%cv, eos_state%cs

    end do
    CLOSE(UNIT=iunit)

    print*,"======= Running CO2 + C2H6 mixture eq. ratio test = 4 ========="

    FileName = trim(OutputFile)//'_CO2_C2H6_Xco2_0_73978.txt'
    iunit = 10
    OPEN(UNIT=iunit,FILE=trim(FileName),FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")

    eos_state % molefrac(1:nspec) = 0.0d0

    eos_state % molefrac(iCO2)  = 0.73978d0
    eos_state % molefrac(iC2H6) = 1.0d0 - eos_state % molefrac(iCO2)

    ! get mass fractions from mole fractions
    call eos_xty(eos_state) 

    eos_state % T = 350.0d0

    do i = 1,numPoints

       ! Generate these quantities over a range of temperature
       eos_state % P = 1.0d0*1e6*10.0d0 + 40.0d0*1e6*10.0*dble(i-1)/(numPoints-1)

       ! call to compute rho and e
       call eos_tp(eos_state)

       write(iunit,'(6(E20.8,4x))') eos_state%P, eos_state % T, eos_state % rho, eos_state%cp, eos_state%cv, eos_state%cs

    end do
    CLOSE(UNIT=iunit)

    print*,"======= Running CO2 + C2H6 mixture eq. ratio test = 5 ========="

    FileName = trim(OutputFile)//'_CO2_C2H6_Xco2_0_90367.txt'
    iunit = 10
    OPEN(UNIT=iunit,FILE=trim(FileName),FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")

    eos_state % molefrac(1:nspec) = 0.0d0

    eos_state % molefrac(iCO2)  = 0.90367d0
    eos_state % molefrac(iC2H6) = 1.0d0 - eos_state % molefrac(iCO2)

    ! get mass fractions from mole fractions
    call eos_xty(eos_state) 

    eos_state % T = 350.0d0

    do i = 1,numPoints

       ! Generate these quantities over a range of temperature
       eos_state % P = 1.0d0*1e6*10.0d0 + 40.0d0*1e6*10.0*dble(i-1)/(numPoints-1)

       ! call to compute rho and e
       call eos_tp(eos_state)

       write(iunit,'(6(E20.8,4x))') eos_state%P, eos_state % T, eos_state % rho, eos_state%cp, eos_state%cv, eos_state%cs

    end do
    CLOSE(UNIT=iunit)

    print*,"======= Running CO2 + C2H6 mixture eq ratio test 6 ========="

    FileName = trim(OutputFile)//'_CO2_C2H6_Xco2_1_0.txt'
    iunit = 10
    OPEN(UNIT=iunit,FILE=trim(FileName),FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")

    eos_state % molefrac(1:nspec) = 0.0d0

    eos_state % molefrac(iCO2)  = 1.0d0
    eos_state % molefrac(iC2H6) = 0.0d0 

    ! get mass fractions from mole fractions
    call eos_xty(eos_state) 

    eos_state % T = 350.0d0

    do i = 1,numPoints

       ! Generate these quantities over a range of temperature
       eos_state % P = 1.0d0*1e6*10.0d0 + 40.0d0*1e6*10.0*dble(i-1)/(numPoints-1)

       ! call to compute rho and e
       call eos_tp(eos_state)

       write(iunit,'(6(E20.8,4x))') eos_state%P, eos_state % T, eos_state % rho, eos_state%cp, eos_state%cv, eos_state%cs

    end do
    CLOSE(UNIT=iunit)
    

!!$       Conversion between various state variables
    print*,"=====Starting testing of state variable conversion=========="
    print*,"          === Testing RHO & T as input -> Compute P and e "
    print*,"          === Testing RHO & P as input -> Compute T and e "
    print*,"          === Testing RHO & E as input -> Compute T and P "
    
    FileName = trim(OutputFile)//'_Convert_StateVar.txt'
    iunit = 10
    OPEN(UNIT=iunit,FILE=trim(FileName),FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")

    ! Default to N2
    eos_state % molefrac(1:nspec) = 0.0d0
    eos_state % molefrac(iN2) = 1.0d0

    eos_rhoT % molefrac(1:nspec) = 0.0d0
    eos_rhoT % molefrac(iN2) = 1.0d0

    eos_rhoP % molefrac(1:nspec) = 0.0d0
    eos_rhoP % molefrac(iN2) = 1.0d0

    eos_rhoE % molefrac(1:nspec) = 0.0d0
    eos_rhoE % molefrac(iN2) = 1.0d0

    ! get mass fractions from mole fractions
    call eos_xty(eos_state)
    call eos_xty(eos_rhoT)
    call eos_xty(eos_rhoP)
    call eos_xty(eos_rhoE)

    eos_state % P = initP

    MaXerrorP_rhoT = -1e12;MaXerrorP_rhoE = -1e12
    MaXerrorT_rhoP = -1e12;MaXerrorT_rhoE = -1e12
    MaXerrorE_rhoP = -1e12;MaXerrorE_rhoT = -1e12

    do i = 1,numPoints

       ! Generate these quantities over a range of temperature
       eos_state % T = Tstart + (Tend-Tstart)*dble(i-1)/(numPoints-1)

       ! call to compute rho and e
       call eos_tp(eos_state)

       ! Use this rho and T as inputs and compute P and e
       eos_rhoT % rho = eos_state % rho
       eos_rhoT % T   = eos_state % T
       
       call eos_rt(eos_rhoT)

       errorP = abs(eos_state%P - eos_rhoT%P)*100.0d0/eos_state%P
       errorE = abs(eos_state%e - eos_rhoT%e)*100.0d0/abs(eos_state%e)

       MaXerrorE_rhoT = max(MaXerrorE_rhoT,errorE)
       MaXerrorP_rhoT = max(MaXerrorP_rhoT,errorP)

       ! Use this rho and P as inputs and compute T and e
       eos_rhoP % rho = eos_state % rho
       eos_rhoP  % P  = eos_state % P
       call eos_rp(eos_rhoP)

       errorE = abs(eos_state%e - eos_rhoP%e)*100.0d0/abs(eos_state%e)
       errorT = abs(eos_state%T - eos_rhoP%T)*100.0d0/eos_state%T

       MaXerrorE_rhoP = max(MaXerrorE_rhoP,errorE)
       MaXerrorT_rhoP = max(MaXerrorT_rhoP,errorT)

       ! Use this rho and e as inputs and compute T and P
       eos_rhoE % rho = eos_state % rho
       eos_rhoE % e =   eos_state % e
       eos_rhoE % P =   eos_state % P
       eos_rhoE % T =   eos_state % T !! Used to generate a guess for Newton solver for Temperature
       call eos_re(eos_rhoE)

       errorP = abs(eos_state%P - eos_rhoE%P)*100.0d0/eos_state%P
       errorT = abs(eos_state%T - eos_rhoE%T)*100.0d0/eos_state%T

       MaXerrorP_rhoE = max(MaXerrorP_rhoE,errorP)
       MaXerrorT_rhoE = max(MaXerrorT_rhoE,errorT)

       write(iunit,'(9(E20.8,4x))') eos_state%P, eos_rhoT%P, eos_rhoE%P, eos_state % T,eos_rhoP % T, eos_rhoE % T , eos_state % e, eos_rhoT % e , eos_rhoP % e

    end do

    CLOSE(UNIT=iunit)
    
    print*,"===== RHO & T as input -> Compute P and e ==="
    print*,"Error in P computation =", MaXerrorP_rhoT, "  ","Error in E computation=", MaXerrorE_rhoT
     
    print*,"===== RHO & P as input -> Compute T and e ==="
    print*,"Error in T computation =", MaXerrorT_rhoP, "  ","Error in E computation=", MaXerrorE_rhoP
     
    print*,"===== RHO & E as input -> Compute P and T  ==="
    print*,"Error in P computation =", MaXerrorP_rhoE, "  ","Error in T computation=", MaXerrorT_rhoE
    
    print*,"=====Finished testing of state variable conversion=========="

!!$    Test derivatives dP/dT, dP/dtau, dhm/dtau and thereafter test species Enthalpy
    print*,"=== Testing various derivatives of EOS  ====="
    
    FileName = trim(OutputFile)//'_Derivatives.txt'
    iunit = 10
    OPEN(UNIT=iunit,FILE=trim(FileName),FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")

    invdelta = 1.0d0/delta

    ! Default to N2
    eos_state % molefrac(1:nspec) = 0.0d0
    eos_state % molefrac(iN2) = 1.0d0

    ! get mass fractions from mole fractions
    call eos_xty(eos_state) 

    eos_state % P = initP

    do i = 1,numPoints

       ! Generate these quantities over a range of temperature
       eos_state % T = Tstart + (Tend-Tstart)*dble(i-1)/(numPoints-1)

       ! call to compute rho and e
       call eos_tp(eos_state)

       ! Analytical expressions for the derivatives
       call eos_deriv(eos_state)

       dtaudY_analy(:) = eos_state%taui(:)

       ! Second derivative of am analytical expression
       d2am_dT2_analy = eos_state%d2AmdT2

       ! Mixture enthalpy
       call eos_h(eos_state)
       hmix = eos_state%h

       ! call to compute P at T and tau,Y
       call eos_rt(eos_state)
       PatT = eos_state%p
       dam_dT_analy = eos_state%damdT

       ! Store am at T
       amT = eos_state%am

       ! call to compute P at T+epsilon and tau,Y
       eos_state % T = eos_state % T + delta
       call eos_rt(eos_state)
       PatTplus = eos_state%p
       amTplus  = eos_state%am

       ! dp/dT at constant tau (or rho)
       dP_dT = (PatTplus - PatT)/delta

       ! Test derivative dam/dT and d2amdT2
       dam_dT   = (amTplus-amT)/delta

       ! call to compute P at T+epsilon and tau,Y
       eos_state % T = eos_state % T - 2.d0* delta
       call eos_rt(eos_state)
       amTminus = eos_state%am

       ! Second derivative of am w.r.t. to Temperature, T
       d2am_dT2 = (amTminus - 2.0d0*amT + amTplus)*invdelta*invdelta

       ! Compute derivative w.r.t to tau (specific volume) constant T 
       eos_state % T = Tstart + (Tend-Tstart)*dble(i-1)/(numPoints-1)

       ! call to compute P at tau and T 
       call eos_rt(eos_state)
       PatTau = eos_state%p

       ! call to compute H at tau and T 
       call eos_h(eos_state)
       HatTau = eos_state%h

       ! call to compute P at tau+Delta and T, Y
       tauDelta = 1.0d0/eos_state % rho + delta
       eos_state % rho = 1.0d0/tauDelta

       call eos_rt(eos_state)
       PatTauplus = eos_state%p

       ! call to compute H at tau+Delta and T,Y
       call eos_h(eos_state)
       HatTauplus = eos_state%h

       ! dp/dTau and dHdTau at constant T
       dP_dtau = (PatTauplus - PatTau)/delta
       dH_dtau = (HatTauplus - HatTau)/delta

       ! Reset rho
       eos_state % P = initP
       call eos_tp(eos_state)

       ! Get species h_i
       call eos_hi(eos_state)
       sumhyk = sum(eos_state%massFrac(:)*eos_state%hi(:))

       ! Derivatives w.r.t to species Y_k

       ! Evaluate Free energy at Yk
       eos_state % massfrac(1:nspec) = 0.0d0
       eos_state % massfrac(iN2) = 1.0d0
       call eos_mui(eos_state)
       fYk = eos_state%f
       sumFyk = sum(eos_state%massFrac(:)*eos_state%mui(:))
       muatYk = eos_state%mui(iN2)

       ! Evaluate Free energy at Yk
       eos_state % massfrac(1:nspec) = 0.0d0
       eos_state % massfrac(iN2) = 1.0d0 + delta
       call eos_mui(eos_state)
       fYkplus = eos_state%f

       muK = (fYkplus - fYk)*invdelta

       ! Output in the file - filename provided by the user
       write(iunit,'(19(E20.8,4x))') initP, eos_state % T, eos_state % rho,eos_state%dpdT,dP_dT, eos_state%dpdtau,dP_dtau,eos_state%dhdr,dH_dtau,hmix,sumhyk,dam_dT_analy, dam_dT,d2am_dT2_analy,d2am_dT2,muK,muatYk,fYk,sumFyk

    end do

    CLOSE(UNIT=iunit)

    call destroy(eos_state)
    call destroy(eos_rhoT)
    call destroy(eos_rhoP)
    call destroy(eos_rhoE)
    call trv_destroy(TestTrans)
    
    stop '====== Finished runnning EOS tests============'

  end subroutine RunEOStests


end module probdata_module
