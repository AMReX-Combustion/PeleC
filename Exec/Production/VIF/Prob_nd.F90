module pc_prob_module

  implicit none

  private

  public :: amrex_probinit, pc_initdata, pc_prob_close

contains
  
  subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(C, name = "amrex_probinit")
    
    use parallel
    use probdata_module
    use bl_error_module
    use bl_constants_module, only: HALF
    use network, only: nspec, naux, molec_wt
    use extern_probin_module, only: const_viscosity, const_bulk_viscosity, const_conductivity, const_diffusivity
    use eos_module

    implicit none

    integer :: init, namlen
    integer :: name(namlen)
    double precision :: problo(3), probhi(3)

    integer untin,i

    namelist /fortin/ p_ref, T_ref, reynolds, mach, prandtl, a_u, omega_u, omega_t

    ! Build "probin" filename -- the name of file containing fortin namelist.
    integer, parameter :: maxlen = 256
    character probin*(maxlen)

    ! Local
    type(eos_t) :: eos_state
    integer, parameter :: out_unit=20

    if (namlen .gt. maxlen) then
       call bl_error('probin file name too long')
    end if

    do i = 1, namlen
       probin(i:i) = char(name(i))
    end do

    ! set namelist defaults here
    p_ref = 1.013d6 ! [erg cm^-3]
    T_ref = 300.d0
    reynolds = 10000.0_dp_t
    mach = 0.1_dp_t
    prandtl = 0.71_dp_t
    a_u = 0.5_dp_t
    omega_u = 5.0_dp_t
    omega_t = 1.0_dp_t

    ! Read namelists
    untin = 9
    open(untin,file=probin(1:namlen),form='formatted',status='old')
    read(untin,fortin)
    close(unit=untin)

    ! set local variable defaults
    L_x = probhi(1) - problo(1)
    L_y = probhi(2) - problo(2)
    L = HALF * L_y

    ! Define the molecular weight for air
    molec_wt = 28.97

    ! Set the equation of state variables
    call build(eos_state)
    eos_state % p = p_ref
    eos_state % T = T_ref
    eos_state % massfrac = 0.d0
    eos_state % massfrac(1) = 1.d0
    call eos_tp(eos_state)

    ! Initial density, velocity, and material properties
    rho0 = eos_state % rho
    eint0 = eos_state % e
    v0 = mach * eos_state % cs
    tau = L / v0
    const_bulk_viscosity = 0.d0
    const_diffusivity = 0.d0
    const_viscosity = rho0 * v0 * L / reynolds
    const_conductivity = const_viscosity * eos_state % cp / prandtl

    ! Write this out to file (might be useful for postprocessing)
    if ( parallel_ioprocessor() ) then
       open(unit=out_unit,file="ic.txt",action="write",status="replace")
       write(out_unit,*)"L, rho0, v0, tau, p_ref, T_ref, gamma, mu, k, c_s0, Reynolds, Mach, Prandtl"
       write(out_unit,*) L, "," , eos_state % rho, "," , v0, "," , tau, "," , eos_state % p, "," , &
            eos_state % T, "," , gamma_const, "," , const_viscosity, "," , &
            const_conductivity, "," , eos_state % cs, "," , &
            reynolds, "," , mach, "," , prandtl
       close(out_unit)
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
       delta,xlo,xhi) bind(C, name="pc_initdata")

    use probdata_module
    use network, only: nspec, naux, molec_wt
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

    integer :: n
    double precision :: u,v,w

    type(eos_t) :: eos_state

    ! Set the equation of state variables
    call build(eos_state)
    eos_state % p = p_ref
    eos_state % T = T_ref
    eos_state % massfrac = 0.d0
    eos_state % massfrac(1) = 1.d0
    call eos_tp(eos_state)

    ! Fill the states
    u = v0
    v = 0.d0
    w = 0.d0
    state(:,:,:,URHO)       = rho0
    do n = 1,nspec
       state(:,:,:,UFS+n-1) = rho0 * eos_state % massfrac(n)
    enddo
    state(:,:,:,UMX)        = rho0 * u
    state(:,:,:,UMY)        = rho0 * v
    state(:,:,:,UMZ)        = rho0 * w
    state(:,:,:,UEINT)      = rho0 * eint0
    state(:,:,:,UEDEN)      = rho0 * (eint0 + HALF * (u**2 + v**2 + w**2))
    state(:,:,:,UTEMP)      = T_ref

  end subroutine pc_initdata

  subroutine pc_prob_close() &
       bind(C, name="pc_prob_close")
  end subroutine pc_prob_close

end module pc_prob_module
