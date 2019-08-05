module pc_prob_module

#include <AMReX_CONSTANTS.H>

  implicit none

  private

  public :: amrex_probinit, pc_initdata, pc_prob_close

contains

  subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(C, name = "amrex_probinit")

    use parallel
    use probdata_module
    use bl_error_module
    use bl_constants_module, only: HALF
    use network, only: nspecies, naux, molec_wt
    use extern_probin_module, only: const_viscosity, const_bulk_viscosity, const_conductivity, const_diffusivity
    use eos_module

    implicit none

    integer :: init, namlen
    integer :: name(namlen)
    double precision :: problo(3), probhi(3)
    integer untin,i

    namelist /fortin/ restart, lambda0, reynolds_lambda0, mach_t0, prandtl, &
                      spectrum_type, mode_start, nmodes, &
                      force_scale, forcing_time_scale_min, forcing_time_scale_max, &
                      forcing_time_scale_min, forcing_time_scale_max, &
                      time_offset, div_free_force

    ! Build "probin" filename -- the name of file containing fortin namelist.
    integer, parameter :: maxlen = 256
    character probin*(maxlen)

    ! Local
    double precision :: twicePi, kxd, kyd, kzd, rn
    double precision :: thetaTmp, phiTmp
    double precision :: cosThetaTmp, cosPhiTmp
    double precision :: sinThetaTmp, sinPhiTmp
    double precision :: px, py, pz, mp2, Ekh
    integer :: kx, ky, kz, mode_count, reduced_mode_count
    integer :: modx, mody, modz, xstep, ystep, zstep
    integer :: iseed
    integer :: nxmodes, nymodes, nzmodes

    double precision :: Lmin
    double precision :: kappa, kappaMax, freqMin, freqMax, freqDiff, pdk


    integer, parameter :: out_unit=20
    type(eos_t) :: eos_state
    integer(kind=8) :: n

    if (namlen .gt. maxlen) then
       call bl_error('probin file name too long')
    end if

    do i = 1, namlen
       probin(i:i) = char(name(i))
    end do

    ! set namelist defaults here
    restart = .false.
    lambda0 = 0.5_dp_t
    reynolds_lambda0 = 100.0_dp_t
    mach_t0 = 0.1_dp_t
    prandtl = 0.71_dp_t

    div_free_force = 0 
    spectrum_type          = 0
    mode_start             = 1
    nmodes                 = 4
    force_scale            = 1.0d10
    forcing_time_scale_min = 0.0
    forcing_time_scale_max = 0.0
  

    ! Initial pressure and temperature
    p0 = 1013250.0d0 ! [erg cm^-3]
    T0 = 1200.d0

    ! Read namelists
    untin = 9
    open(untin,file=probin(1:namlen),form='formatted',status='old')
    read(untin,fortin)
    close(unit=untin)

    ! set local variable defaults
    Lx = probhi(1) - problo(1)
    Ly = probhi(2) - problo(2)
    Lz = probhi(3) - problo(3)

    ! Define the molecular weight for air
    molec_wt = 28.97

    ! Wavelength associated to Taylor length scale
    k0 = 2.d0/lambda0

    ! Set the equation of state variables
    call build(eos_state)
    eos_state % p = p0
    eos_state % T = T0
    eos_state % massfrac    = 0.d0
    eos_state % massfrac(1) = 1.0d0
    call eos_tp(eos_state)

    ! Initial density, velocity, and material properties
    rho0  = eos_state % rho
    eint0 = eos_state % e
    urms0 = mach_t0 * eos_state % cs / sqrt(3.d0)
    tau  = lambda0 / urms0
    const_bulk_viscosity = 0.d0
    const_diffusivity = 0.d0
    const_viscosity = rho0 * urms0 * lambda0 / reynolds_lambda0
    const_conductivity = const_viscosity * eos_state % cp / prandtl

    ! Write this out to file (might be useful for postprocessing)
     if ( parallel_ioprocessor() ) then
       open(unit=out_unit,file="ic.txt",action="write",status="replace")
       write(out_unit,*)"lambda0, k0, rho0, urms0, tau, p0, T0, mu, k, c_s0, Reynolds, Mach, Prandtl"
       write(out_unit,*) lambda0, "," , k0, "," , eos_state % rho, "," , urms0, "," , tau, "," , &
            eos_state % p, "," , eos_state % T, "," , const_viscosity, "," , &
            const_conductivity, "," , eos_state % cs, "," , &
            reynolds_lambda0, "," , mach_t0, "," , prandtl
       close(out_unit)
    endif


    if ( parallel_ioprocessor() ) then
       write (*,*) "Initialising random number generator..."
    end if

    twicePi = two*Pi

    if (blrandseed.gt.0) then
       call blutilinitrand(blrandseed)
       call blutilrand(rn)
       call blutilinitrand(blrandseed)
       if ( parallel_ioprocessor() ) then
          write (*,*) "blrandseed = ",blrandseed
          write (*,*) "first random number = ",rn
       endif
    else
       call blutilinitrand(111397)
    endif

    if ( parallel_ioprocessor() ) then
       write(*,*) "Lx = ",Lx
       write(*,*) "Ly = ",Ly
       write(*,*) "Lz = ",Lz
    endif

    Lmin = min(Lx,Ly,Lz)
    kappaMax = dfloat(nmodes)/Lmin + 1.0d-8
    nxmodes = nmodes*int(0.5+Lx/Lmin)
    nymodes = nmodes*int(0.5+Ly/Lmin)
    nzmodes = nmodes*int(0.5+Lz/Lmin)
    if ( parallel_ioprocessor() ) then 
      write(*,*) "Lmin = ",Lmin
      write(*,*) "kappaMax = ",kappaMax
      write(*,*) "nxmodes = ",nxmodes
      write(*,*) "nymodes = ",nymodes
      write(*,*) "nzmodes = ",nzmodes
    endif

    allocate(FTX(1:nxmodes, 1:nymodes, 1:nzmodes))
    allocate(FTY(1:nxmodes, 1:nymodes, 1:nzmodes))
    allocate(FTZ(1:nxmodes, 1:nymodes, 1:nzmodes))
    allocate(TAT(1:nxmodes, 1:nymodes, 1:nzmodes))
    allocate(TAP(1:nxmodes, 1:nymodes, 1:nzmodes))

    allocate(FPX(1:nxmodes, 1:nymodes, 1:nzmodes))
    allocate(FPY(1:nxmodes, 1:nymodes, 1:nzmodes))
    allocate(FPZ(1:nxmodes, 1:nymodes, 1:nzmodes))
    allocate(FAX(1:nxmodes, 1:nymodes, 1:nzmodes))
    allocate(FAY(1:nxmodes, 1:nymodes, 1:nzmodes))
    allocate(FAZ(1:nxmodes, 1:nymodes, 1:nzmodes))

    allocate(FPXX(1:nxmodes, 1:nymodes, 1:nzmodes))
    allocate(FPYX(1:nxmodes, 1:nymodes, 1:nzmodes))
    allocate(FPZX(1:nxmodes, 1:nymodes, 1:nzmodes))
    allocate(FPXY(1:nxmodes, 1:nymodes, 1:nzmodes))
    allocate(FPYY(1:nxmodes, 1:nymodes, 1:nzmodes))
    allocate(FPZY(1:nxmodes, 1:nymodes, 1:nzmodes))
    allocate(FPXZ(1:nxmodes, 1:nymodes, 1:nzmodes))
    allocate(FPYZ(1:nxmodes, 1:nymodes, 1:nzmodes))
    allocate(FPZZ(1:nxmodes, 1:nymodes, 1:nzmodes))


    if (forcing_time_scale_min.eq.zero) then
      forcing_time_scale_min = half
    endif
    if (forcing_time_scale_max.eq.zero) then
       forcing_time_scale_max = one
    endif

    freqMin = one/forcing_time_scale_max
    freqMax = one/forcing_time_scale_min
    freqDiff= freqMax-freqMin

    if ( parallel_ioprocessor() ) then
       write(*,*) "forcing_time_scale_min = ",forcing_time_scale_min
       write(*,*) "forcing_time_scale_max = ",forcing_time_scale_max
       write(*,*) "freqMin = ",freqMin
       write(*,*) "freqMax = ",freqMax
       write(*,*) "freqDiff = ",freqDiff
    endif

    mode_count = 0

    xstep = int(Lx/Lmin+0.5)
    ystep = int(Ly/Lmin+0.5)
    zstep = int(Lz/Lmin+0.5)
    if ( parallel_ioprocessor() ) then
      write (*,*) "Mode step ",xstep, ystep, zstep
    endif

    do kz = mode_start*zstep, nzmodes, zstep
       kzd = dfloat(kz)
       do ky = mode_start*ystep, nymodes, ystep
          kyd = dfloat(ky)
          do kx = mode_start*xstep, nxmodes, xstep
             kxd = dfloat(kx)

             kappa = sqrt( (kxd*kxd)/(Lx*Lx) + (kyd*kyd)/(Ly*Ly) + (kzd*kzd)/(Lz*Lz) )

             if (kappa.le.kappaMax) then
                call blutilrand(rn)
                FTX(kx,ky,kz) = (freqMin + freqDiff*rn)*twicePi
                call blutilrand(rn)
                FTY(kx,ky,kz) = (freqMin + freqDiff*rn)*twicePi
                call blutilrand(rn)
                FTZ(kx,ky,kz) = (freqMin + freqDiff*rn)*twicePi
!     Translation angles, theta=0..2Pi and phi=0..Pi
                call blutilrand(rn)
                TAT(kx,ky,kz) = rn*twicePi
                call blutilrand(rn)
                TAP(kx,ky,kz) = rn*Pi
!     Phases
                call blutilrand(rn)
                FPX(kx,ky,kz) = rn*twicePi
                call blutilrand(rn)
                FPY(kx,ky,kz) = rn*twicePi
                call blutilrand(rn)
                FPZ(kx,ky,kz) = rn*twicePi
                if (div_free_force.eq.1) then
                  call blutilrand(rn)
                  FPXX(kx,ky,kz) = rn*twicePi
                  call blutilrand(rn)
                  FPYX(kx,ky,kz) = rn*twicePi
                  call blutilrand(rn)
                  FPZX(kx,ky,kz) = rn*twicePi
                  call blutilrand(rn)
                  FPXY(kx,ky,kz) = rn*twicePi
                  call blutilrand(rn)
                  FPYY(kx,ky,kz) = rn*twicePi
                  call blutilrand(rn)
                  FPZY(kx,ky,kz) = rn*twicePi
                  call blutilrand(rn)
                  FPXZ(kx,ky,kz) = rn*twicePi
                  call blutilrand(rn)
                  FPYZ(kx,ky,kz) = rn*twicePi
                  call blutilrand(rn)
                  FPZZ(kx,ky,kz) = rn*twicePi
                endif
!     Amplitudes (alpha)
                call blutilrand(rn)
                thetaTmp      = rn*twicePi
                cosThetaTmp   = cos(thetaTmp)
                sinThetaTmp   = sin(thetaTmp)
                call blutilrand(rn)
                phiTmp        = rn*Pi
                cosPhiTmp     = cos(phiTmp)
                sinPhiTmp     = sin(phiTmp)

                px = cosThetaTmp * sinPhiTmp
                py = sinThetaTmp * sinPhiTmp
                pz =               cosPhiTmp

                mp2           = px*px + py*py + pz*pz
                if (kappa .lt. 0.000001) then
                  if ( parallel_ioprocessor() ) then
                     write(*,*) "ZERO AMPLITUDE MODE ",kx,ky,kz
                  endif
                  FAX(kx,ky,kz) = zero
                  FAY(kx,ky,kz) = zero
                  FAZ(kx,ky,kz) = zero
                else
!     Count modes that contribute
                   mode_count = mode_count + 1
!     Set amplitudes
                   if (spectrum_type.eq.1) then
                       Ekh        = one / kappa
                   else if (spectrum_type.eq.2) then
                       Ekh        = one / (kappa*kappa)
                   else
                       Ekh        = one
                   endif

                   if (div_free_force.eq.1) then
                       Ekh = Ekh / kappa
                   endif

                   if (moderate_zero_modes.eq.1) then
                           if (kx.eq.0) Ekh = Ekh / two
                           if (ky.eq.0) Ekh = Ekh / two
                           if (kz.eq.0) Ekh = Ekh / two
                   endif

                   if (force_scale.gt.zero) then
                           FAX(kx,ky,kz) = force_scale * px * Ekh / mp2
                           FAY(kx,ky,kz) = force_scale * py * Ekh / mp2
                           FAZ(kx,ky,kz) = force_scale * pz * Ekh / mp2
                   else
                           FAX(kx,ky,kz) = px * Ekh / mp2
                           FAY(kx,ky,kz) = py * Ekh / mp2
                           FAZ(kx,ky,kz) = pz * Ekh / mp2
                   endif

                   if ( parallel_ioprocessor() ) then 
                           write (*,*) "Mode"
                           write (*,*) "kappa = ",kx,ky,kz,kappa, sqrt(FAX(kx,ky,kz)**2+FAY(kx,ky,kz)**2+FAZ(kx,ky,kz)**2)
                           write (*,*) "Amplitudes - A"
                           write (*,*) FAX(kx,ky,kz), FAY(kx,ky,kz), FAZ(kx,ky,kz)
                           write (*,*) "Frequencies"
                           write (*,*) FTX(kx,ky,kz), FTY(kx,ky,kz), FTZ(kx,ky,kz)
                    endif

                endif
             endif


          enddo
       enddo
    enddo


!    Now let's break symmetry, have to assume high aspect ratio in z for now
    reduced_mode_count = 0
    do kz = 1, zstep - 1
      kzd = dfloat(kz)
      do ky = mode_start, nymodes
         kyd = dfloat(ky)
         do kx = mode_start, nxmodes
            kxd = dfloat(kx)

            kappa = sqrt( (kxd*kxd)/(Lx*Lx) + (kyd*kyd)/(Ly*Ly) + (kzd*kzd)/(Lz*Lz) )
            if (kappa.le.kappaMax) then
                call blutilrand(rn)
                FTX(kx,ky,kz) = (freqMin + freqDiff*rn)*twicePi
                call blutilrand(rn)
                FTY(kx,ky,kz) = (freqMin + freqDiff*rn)*twicePi
                call blutilrand(rn)
                FTZ(kx,ky,kz) = (freqMin + freqDiff*rn)*twicePi
!     Translation angles, theta=0..2Pi and phi=0..Pi
                call blutilrand(rn)
                TAT(kx,ky,kz) = rn*twicePi
                call blutilrand(rn)
                TAP(kx,ky,kz) = rn*Pi
!     Phases
                call blutilrand(rn)
                FPX(kx,ky,kz) = rn*twicePi
                call blutilrand(rn)
                FPY(kx,ky,kz) = rn*twicePi
                call blutilrand(rn)
                FPZ(kx,ky,kz) = rn*twicePi
                if (div_free_force.eq.1) then
                   call blutilrand(rn)
                   FPXX(kx,ky,kz) = rn*twicePi
                   call blutilrand(rn)
                   FPYX(kx,ky,kz) = rn*twicePi
                   call blutilrand(rn)
                   FPZX(kx,ky,kz) = rn*twicePi
                   call blutilrand(rn)
                   FPXY(kx,ky,kz) = rn*twicePi
                   call blutilrand(rn)
                   FPYY(kx,ky,kz) = rn*twicePi
                   call blutilrand(rn)
                   FPZY(kx,ky,kz) = rn*twicePi
                   call blutilrand(rn)
                   FPXZ(kx,ky,kz) = rn*twicePi
                   call blutilrand(rn)
                   FPYZ(kx,ky,kz) = rn*twicePi
                   call blutilrand(rn)
                   FPZZ(kx,ky,kz) = rn*twicePi
                endif
!     Amplitudes (alpha)
                call blutilrand(rn)
                thetaTmp      = rn*twicePi
                cosThetaTmp   = cos(thetaTmp)
                sinThetaTmp   = sin(thetaTmp)
                call blutilrand(rn)
                phiTmp        = rn*Pi
                cosPhiTmp     = cos(phiTmp)
                sinPhiTmp     = sin(phiTmp)

                px = cosThetaTmp * sinPhiTmp
                py = sinThetaTmp * sinPhiTmp
                pz =               cosPhiTmp

                mp2           = px*px + py*py + pz*pz
                if (kappa .lt. 0.000001) then
                   if ( parallel_ioprocessor() ) then  
                      write(*,*) "ZERO AMPLITUDE MODE ",kx,ky,kz
                   endif
                   FAX(kx,ky,kz) = zero
                   FAY(kx,ky,kz) = zero
                   FAZ(kx,ky,kz) = zero
                else
!     Count modes that contribute
                  reduced_mode_count = reduced_mode_count + 1
!     Set amplitudes
                  if (spectrum_type.eq.1) then
                    Ekh        = one / kappa
                  else if (spectrum_type.eq.2) then
                    Ekh        = one / (kappa*kappa)
                  else
                    Ekh        = one
                  endif
         
                  if (div_free_force.eq.1) then
                    Ekh = Ekh / kappa
                  endif

                  if (moderate_zero_modes.eq.1) then
                    if (kx.eq.0) Ekh = Ekh / two
                    if (ky.eq.0) Ekh = Ekh / two
                    if (kz.eq.0) Ekh = Ekh / two
                  endif
                  if (force_scale.gt.zero) then
                    FAX(kx,ky,kz) = forcing_epsilon * force_scale * px * Ekh / mp2
                    FAY(kx,ky,kz) = forcing_epsilon * force_scale * py * Ekh / mp2
                    FAZ(kx,ky,kz) = forcing_epsilon * force_scale * pz * Ekh / mp2
                  else
                    FAX(kx,ky,kz) = forcing_epsilon * px * Ekh / mp2
                    FAY(kx,ky,kz) = forcing_epsilon * py * Ekh / mp2
                    FAZ(kx,ky,kz) = forcing_epsilon * pz * Ekh / mp2
                  endif

                  if ( parallel_ioprocessor() ) then       
                     write (*,*) "Mode"
                     write (*,*) "kappa = ",kx,ky,kz,kappa, sqrt(FAX(kx,ky,kz)**2+FAY(kx,ky,kz)**2+FAZ(kx,ky,kz)**2)
                     write (*,*) "Amplitudes - B"
                     write (*,*) FAX(kx,ky,kz), FAY(kx,ky,kz), FAZ(kx,ky,kz)
                     write (*,*) "Frequencies"
                     write (*,*) FTX(kx,ky,kz), FTY(kx,ky,kz), FTZ(kx,ky,kz)
                  endif
                endif
            endif


          enddo
       enddo
    enddo

    if ( parallel_ioprocessor() ) then
       write(*,*) "mode_count = ",mode_count
       write(*,*) "reduced_mode_count = ",reduced_mode_count
       if (spectrum_type.eq.1) then
          write (*,*) "Spectrum type 1"
       else if (spectrum_type.eq.2) then
          write (*,*) "Spectrum type 2"
       else
          write (*,*) "Spectrum type OTHER"
       endif
    endif


  end subroutine amrex_probinit


  ! ::: -----------------------------------------------------------
  ! ::: This routine is called at problem setup time and is used
  ! ::: to initialize data on each grid.
  ! :::
  ! ::: NOTE:  all arrays have one cell of ghost zones surrounding
  ! :::        the grid interior.  Values in these cells need not
  ! :::        be set here.
  ! :::
  ! ::: INPUTS/OUTPUTS:
  ! :::
  ! ::: level     => amr level of grid
  ! ::: time      => time at which to init data
  ! ::: lo,hi     => index limits of grid interior (cell centered)
  ! ::: nstate    => number of state components.  You should know
  ! :::		   this already!
  ! ::: state     <=  Scalar array
  ! ::: delta     => cell size
  ! ::: xlo,xhi   => physical locations of lower left and upper
  ! :::              right hand corner of grid.  (does not include
  ! :::		   ghost region).
  ! ::: -----------------------------------------------------------
  subroutine pc_initdata(level,time,lo,hi,nvar, &
       state,state_lo,state_hi, &
       delta,xlo,xhi) bind(C, name = "pc_initdata")

    use probdata_module
    use network, only: nspecies, naux, molec_wt
    use eos_type_module
    use meth_params_module, only : URHO, UMX, UMY, UMZ, &
         UEDEN, UEINT, UFS, UTEMP
    use bl_constants_module, only: ZERO, HALF, M_PI
    use eos_module

    implicit none

    integer :: level, nvar
    integer :: lo(3), hi(3)
    integer :: state_lo(3),state_hi(3)
    double precision :: xlo(3), xhi(3), time, delta(3)
    double precision :: state(state_lo(1):state_hi(1), &
         state_lo(2):state_hi(2), &
         state_lo(3):state_hi(3),nvar)

    type(eos_t) :: eos_state
    integer :: i, j, k
    double precision :: x, y, z
    double precision :: uinterp, vinterp, winterp

    ! Define the molecular weight for air
    molec_wt = 28.97

    ! Set the equation of state variables
    call build(eos_state)
    eos_state % p = p0
    eos_state % T = T0
    eos_state % massfrac    = 0.d0
    eos_state % massfrac(1) = 1.0d0

    call eos_tp(eos_state)

    ! Uniform density, temperature, pressure (internal energy)
    state(:,:,:,URHO)            = rho0
    do i = 1,nspecies
       state(:,:,:,UFS+i-1)      = rho0 * eos_state % massfrac(i)
    enddo
    state(:,:,:,UTEMP)           = T0
    state(:,:,:,UEINT)           = rho0 * eint0

    ! Fill in the velocities and energy.
    do k = lo(3), hi(3)
       z = xlo(3) + delta(3)*(dble(k-lo(3)) + HALF)

       do j = lo(2), hi(2)
          y = xlo(2) + delta(2)*(dble(j-lo(2)) + HALF)

          do i = lo(1), hi(1)
             x = xlo(1) + delta(1)*(dble(i-lo(1)) + HALF)

             uinterp = 0.0d0
             vinterp = 0.0d0
             winterp = 0.0d0

             state(i,j,k,UMX)   = state(i,j,k,URHO) * uinterp
             state(i,j,k,UMY)   = state(i,j,k,URHO) * vinterp
             state(i,j,k,UMZ)   = state(i,j,k,URHO) * winterp
             state(i,j,k,UEDEN) = state(i,j,k,UEINT) + HALF * state(i,j,k,URHO) * (uinterp**2 + vinterp**2 + winterp**2)
          enddo
       enddo
    enddo

  end subroutine pc_initdata

  subroutine pc_prob_close() &
       bind(C, name="pc_prob_close")

  end subroutine pc_prob_close



  ! ::: -----------------------------------------------------------
  ! ::: Calculate the Taylor length scale in x-direction.
  ! :::    lambda = <u^2>/<(du/dx)^2>
  ! :::
  ! ::: INPUTS/OUTPUTS:
  ! :::
  ! ::: u(n,n,n) => input velocity array
  ! ::: x(n,n,n) => input coordinate array
  ! ::: n        => size of array
  ! ::: lambda  <=> lambda
  ! ::: -----------------------------------------------------------
  subroutine calculate_taylor_microscale(u, x, n, lambda)

    implicit none

    double precision, intent(in) :: u(n,n,n)
    double precision, intent(in) :: x(n,n,n)
    integer, intent(in) :: n
    double precision, intent(out) :: lambda

    ! Local variables
    integer :: i
    double precision :: dudx2 = 0, u2 = 0

    ! Square of velocity
    u2 = sum(u(:,:,:)*u(:,:,:))

    ! Calculate the gradients. Assume periodicity for the edges.
    do i = 2, n-1
       dudx2 = dudx2 + sum((u(i+1,:,:) - u(i-1,:,:))**2 / (x(i+1,:,:) - x(i-1,:,:))**2)
    enddo
    dudx2 = dudx2 + sum((u(2,:,:) - u(n,:,:))**2 / (x(2,:,:) - x(n,:,:))**2)
    dudx2 = dudx2 + sum((u(1,:,:) - u(n-1,:,:))**2 / (x(1,:,:) - x(n-1,:,:))**2)

    lambda = sqrt(u2/dudx2)

    return
  end subroutine calculate_taylor_microscale

end module pc_prob_module
