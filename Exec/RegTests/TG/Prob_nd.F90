module pc_prob_module

  implicit none

  private

  public :: amrex_probinit, pc_initdata, pc_prob_close

contains

  subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(C, name = "amrex_probinit")

    use probdata_module
    use amrex_fort_module
    implicit none

    integer :: init, namlen
    integer :: name(namlen)
    double precision :: problo(3), probhi(3)

    integer untin,i

    namelist /fortin/ reynolds, mach, prandtl, convecting, omega_x, omega_y, omega_z

    ! Build "probin" filename -- the name of file containing fortin namelist.
    integer, parameter :: maxlen = 256
    character probin*(maxlen)

    if (namlen .gt. maxlen) then
       call bl_error('probin file name too long')
    end if

    do i = 1, namlen
       probin(i:i) = char(name(i))
    end do

    ! set namelist defaults here
    reynolds = 1600.0_amrex_real
    mach = 0.1_amrex_real
    prandtl = 0.71_amrex_real
    convecting = .false.
    omega_x = 1.0_amrex_real
    omega_y = 1.0_amrex_real
    omega_z = 1.0_amrex_real

    ! Read namelists
    untin = 9
    open(untin,file=probin(1:namlen),form='formatted',status='old')
    read(untin,fortin)
    close(unit=untin)

    ! set local variable defaults
    L_x = probhi(1) - problo(1)
    L_y = probhi(2) - problo(2)
    L_z = probhi(3) - problo(3)

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
       delta,xlo,xhi) bind(C, name="pc_initdata")

    use amrex_paralleldescriptor_module, only: amrex_pd_ioprocessor
    use probdata_module
    use network, only: nspecies, naux, molec_wt
    use eos_type_module
    use meth_params_module, only : URHO, UMX, UMY, UMZ, &
         UEDEN, UEINT, UFS, UTEMP, small_temp
    use amrex_constants_module, only: ZERO, HALF, M_PI
    use extern_probin_module, only: const_viscosity, const_bulk_viscosity, const_conductivity, const_diffusivity
    use eos_module

    implicit none

    integer :: level, nvar
    integer :: lo(3), hi(3)
    integer :: state_lo(3),state_hi(3)
    double precision :: xlo(3), xhi(3), time, delta(3)
    double precision :: state(state_lo(1):state_hi(1), &
         state_lo(2):state_hi(2), &
         state_lo(3):state_hi(3),nvar)

    integer :: i,j,k
    double precision :: x,y,z,rho,u,v,w,p,eint
    double precision :: L, rho0, v0, p0, T0
    integer, parameter :: out_unit=20

    type(eos_t) :: eos_state

    call build(eos_state)

    ! Define the molecular weight for air
    molec_wt = 28.97

    ! Define the length scale
    L = 1.d0/M_PI

    ! Initial pressure and temperature
    p0 = 1.013d6 ! [erg cm^-3]
    T0 = 300.d0

    ! Set the equation of state variables
    eos_state % p = p0
    eos_state % T = T0
    eos_state % massfrac    = 0.d0
    eos_state % massfrac(1) = 1.d0
    call eos_tp(eos_state)

    ! Initial density, velocity, and material properties
    rho0 = eos_state % rho
    v0   = mach * eos_state % cs
    const_bulk_viscosity = 0.d0
    const_diffusivity = 0.d0
    const_viscosity = rho0 * v0 * L / reynolds
    const_conductivity = const_viscosity * eos_state % cp / prandtl
    state(:,:,:,UTEMP) = T0

    ! Write this out to file (might be useful for postprocessing)
    if ( amrex_pd_ioprocessor() ) then
       open(unit=out_unit,file="ic.txt",action="write",status="replace")
       write(out_unit,*)"L, rho0, v0, p0, T0, gamma, mu, k, c_s0, Reynolds, Mach, Prandtl, omega_x, omega_y, omega_z"
       write(out_unit,*) L, "," , eos_state % rho, "," , v0, "," , eos_state % p, "," , &
            eos_state % T, "," , gamma_const, "," , const_viscosity, "," , &
            const_conductivity, "," , eos_state % cs, "," , &
            reynolds, "," , mach, "," , prandtl, ",", omega_x, ",", omega_y, ",", omega_z
       close(out_unit)
    endif

    do k = lo(3), hi(3)
       z = xlo(3) + delta(3)*(dble(k-lo(3)) + HALF)

       do j = lo(2), hi(2)
          y = xlo(2) + delta(2)*(dble(j-lo(2)) + HALF)

          do i = lo(1), hi(1)
             x = xlo(1) + delta(1)*(dble(i-lo(1)) + HALF)

             u =  v0 * sin(omega_x * x/L) * cos(omega_y * y/L) * cos(omega_z * z/L)
             v = -v0 * cos(omega_x * x/L) * sin(omega_y * y/L) * cos(omega_z * z/L)
             if (convecting) then
                u = u + v0
                v = v + v0
             endif
             w = 0.d0

             eos_state % p   =  p0 + rho0*v0*v0/16.d0 * &
                  ( cos(2.d0*omega_x * x/L) + cos(2.d0*omega_y * y/L) ) * (cos(2.d0*omega_z * z/L) + 2.d0)

             ! Call EOS by specifying the temperature and pressure
             call eos_tp(eos_state)
             rho  = eos_state % rho
             eint = eos_state % e

             ! Fill the states
             state(i,j,k,URHO)            = rho
             state(i,j,k,UFS:UFS+nspecies-1) = rho * eos_state % massfrac(1:nspecies)
             state(i,j,k,UMX)             = rho * u
             state(i,j,k,UMY)             = rho * v
             state(i,j,k,UMZ)             = rho * w
             state(i,j,k,UEINT)           = rho * eint
             state(i,j,k,UEDEN)           = rho * (eint + HALF * (u**2 + v**2 + w**2))

          enddo
       enddo
    enddo

  end subroutine pc_initdata

  subroutine pc_prob_close() &
       bind(C, name="pc_prob_close")
  end subroutine pc_prob_close

end module pc_prob_module
