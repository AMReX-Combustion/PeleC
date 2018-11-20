module pc_prob_module

  implicit none

  private

  public :: amrex_probinit, pc_initdata, pc_prob_close

contains

  subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(C, name = "amrex_probinit")

    use probdata_module
    use bl_error_module
    implicit none

    integer :: init, namlen
    integer :: name(namlen)
    double precision :: problo(3), probhi(3)

    integer untin,i

    namelist /fortin/ mach, beta, radius, xc, yc

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
    mach = 0.05
    beta = 0.02
    radius = 0.5
    xc = 5.0
    yc = 5.0

    ! Read namelists
    untin = 9
    open(untin,file=probin(1:namlen),form='formatted',status='old')
    read(untin,fortin)
    close(unit=untin)

    ! set local variable defaults
    L_x = probhi(1) - problo(1)
    L_y = probhi(2) - problo(2)

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
    use network, only: nspec, naux, molec_wt
    use eos_type_module
    use meth_params_module, only : URHO, UMX, UMY, &
         UEDEN, UEINT, UFS, UTEMP
    use bl_constants_module, only: ZERO, ONE, HALF, TWO
    use extern_probin_module, only: const_viscosity, const_bulk_viscosity, const_conductivity, const_diffusivity
    use eos_module

    implicit none

    integer :: level, nvar
    integer :: lo(2), hi(2)
    integer :: state_lo(2),state_hi(2)
    double precision :: xlo(2), xhi(2), time, delta(2)
    double precision :: state(state_lo(1):state_hi(1), &
         state_lo(2):state_hi(2), nvar)

    integer :: i, j
    double precision :: x, y, r2
    double precision :: rho0, u0, p0, T0, tau
    double precision :: u, v, eint
    integer, parameter :: out_unit=20

    type(eos_t) :: eos_state

    ! Define the molecular weight for air
    molec_wt = 28.97

    ! Initial pressure and temperature
    p0 = 1.0d6 ! [erg cm^-3]
    T0 = 300.d0

    ! Set the equation of state variables
    call build(eos_state)
    eos_state % p = p0
    eos_state % T = T0
    eos_state % massfrac    = 0.d0
    eos_state % massfrac(1) = 1.d0
    call eos_tp(eos_state)

    ! Initial density, velocity, and material properties
    rho0  = eos_state % rho
    u0 = mach * eos_state % cs
    tau  = L_x / u0
    const_bulk_viscosity = 0.d0
    const_diffusivity = 0.d0
    const_viscosity = 0.d0
    const_conductivity = 0.d0

    ! Write this out to file (might be useful for postprocessing)
    open(unit=out_unit,file="ic.txt",action="write",status="replace")
    write(out_unit,*)"rho0, u0, tau, p0, T0, gamma, c_s0, Mach, beta, radius, xc, yc"
    write(out_unit,*) eos_state % rho, "," , u0, "," , tau, "," , eos_state % p, "," , &
         eos_state % T, "," , gamma_const, ",", eos_state % cs, "," , mach, "," , &
         beta, "," , radius, "," ,  xc, "," , yc
    close(out_unit)

    ! Fill in the velocities and energy.
    do j = lo(2), hi(2)
       y = xlo(2) + delta(2)*(dble(j-lo(2)) + HALF)

       do i = lo(1), hi(1)
          x = xlo(1) + delta(1)*(dble(i-lo(1)) + HALF)

          r2 = ((x-xc)**2 + (y-yc)**2) / (radius**2)
          u =  u0 * ( ONE - beta * (y-yc)/radius * exp(-HALF*r2))
          v =  u0 * beta * (x-xc)/radius * exp(-HALF*r2)

          ! Call EOS by specifying the temperature and density
          eos_state % T = T0 - (u0*beta)**2 * exp(-HALF*r2) / (TWO * eos_state % cp)
          eos_state % rho = rho0 * (eos_state % T /T0)**(ONE/(gamma_const - ONE))
          call eos_rt(eos_state)
          eint = eos_state % e

          ! Fill the states
          state(i,j,URHO)            = eos_state % rho
          state(i,j,UFS:UFS+nspec-1) = eos_state % rho * eos_state % massfrac(1:nspec)
          state(i,j,UMX)             = eos_state % rho * u
          state(i,j,UMY)             = eos_state % rho * v
          state(i,j,UTEMP)           = eos_state % T
          state(i,j,UEINT)           = eos_state % rho * eint
          state(i,j,UEDEN)           = eos_state % rho * (eint + HALF * (u**2 + v**2))
       enddo
    enddo

  end subroutine pc_initdata

  subroutine pc_prob_close() &
       bind(C, name="pc_prob_close")

  end subroutine pc_prob_close

end module pc_prob_module
