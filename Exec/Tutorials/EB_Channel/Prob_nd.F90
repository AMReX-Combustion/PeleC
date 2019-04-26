module pc_prob_module

  implicit none

  private

  public :: amrex_probinit, pc_initdata, pc_prob_close

contains
  
  subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(C, name = "amrex_probinit")
    
    use probdata_module
    use amrex_fort_module
    use eos_module

    implicit none
    integer :: init, namlen
    integer :: name(namlen)
    double precision :: problo(3), probhi(3)

    integer :: untin,i

    type(eos_t) :: eos_state
    namelist /fortin/ reynolds, mach, prandtl, angle

    ! Build "probin" filename -- the name of file containing fortin namelist.
    integer, parameter :: maxlen = 256
    character :: probin*(maxlen)

    if (namlen .gt. maxlen) then
       call bl_error('probin file name too long')
    end if

    do i = 1, namlen
       probin(i:i) = char(name(i))
    end do

    ! set namelist defaults
    reynolds = 100.0_amrex_real
    mach = 0.1_amrex_real
    prandtl = 0.71_amrex_real
    half_channel = 1.0_amrex_real
    angle = 0.0_amrex_real
   
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
  ! ::: nvar      => number of state components.
  ! ::: state     <= scalar array
  ! ::: delta     => cell size
  ! ::: xlo, xhi  => physical locations of lower left and upper
  ! :::              right hand corner of grid.  (does not include
  ! :::		   ghost region).
  ! ::: -----------------------------------------------------------
  subroutine pc_initdata(level,time,lo,hi,nvar, &
       state,state_lo,state_hi, &
       delta,xlo,xhi) bind(C, name="pc_initdata")
    use amrex_paralleldescriptor_module, only: amrex_pd_ioprocessor
    use network, only: nspec, molec_wt
    use probdata_module
    use meth_params_module, only : URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, UFS
    use amrex_constants_module, only: M_PI, HALF
    use extern_probin_module, only: const_viscosity, const_bulk_viscosity, const_conductivity, const_diffusivity
    use eos_type_module
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
    integer, parameter :: out_unit=20
    double precision :: x, y, z, cost, sint, yp, u, v, w
    type(eos_t) :: eos_state

    call build(eos_state)

    ! Define the molecular weight for air
    molec_wt = 28.97

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
    cost = cos(M_PI/180.d0 * angle)
    sint = sin(M_PI/180.d0 * angle)
    u0 = mach * eos_state % cs * cost
    v0 = mach * eos_state % cs * sint
    w0 = 0.d0
    umag0 = sqrt(u0**2 + v0**2 + w0**2)
    eint0 = eos_state % e
    const_bulk_viscosity = 0.d0
    const_diffusivity = 0.d0
    const_viscosity = rho0 * umag0 * half_channel / reynolds
    const_conductivity = const_viscosity * eos_state % cp / prandtl
    dpdx = - 2.d0 * const_viscosity * umag0 / (half_channel ** 2) 

    ! Write this out to file (might be useful for postprocessing)
    if ( amrex_pd_ioprocessor() ) then
       open(unit=out_unit,file="ic.txt",action="write",status="replace")
       write(out_unit,*)"half_channel, rho0, u0, v0, w0, p0, T0, gamma, mu, k, c_s0, dpdx, Reynolds, Mach, Prandtl, angle"
       write(out_unit,*) half_channel, "," , eos_state % rho, "," , u0, "," , v0, ",", w0, "," , eos_state % p, "," , &
            eos_state % T, "," , gamma_const, "," , const_viscosity, "," , &
            const_conductivity, "," , eos_state % cs, "," , &
            dpdx, "," , reynolds, "," , mach, "," , prandtl, ",", angle
       close(out_unit)
    endif

    do k = lo(3), hi(3)
       z = xlo(3) + delta(3)*(dble(k-lo(3)) + HALF)
       do j = lo(2), hi(2)
          y = xlo(2) + delta(2)*(dble(j-lo(2)) + HALF)
          do i = lo(1), hi(1)
             x = xlo(1) + delta(1)*(dble(i-lo(1)) + HALF)

             eos_state % p = dpdx * (x * cost + y * sint) + p0
             eos_state % T = T0
             eos_state % massfrac    = 0.d0
             eos_state % massfrac(1) = 1.d0
             call eos_tp(eos_state)

             ! Use a parabolic profile (perturbed by 10% relative to
             ! the exact) to initialize the velocities
             yp = -x * sint + y * cost
             u = 1.1 * u0 * (1.d0 - (yp/half_channel)**2)
             v = 1.1 * v0 * (1.d0 - (yp/half_channel)**2)
             w = w0
             
             state(i,j,k,URHO)            = eos_state % rho
             state(i,j,k,UFS:UFS+nspec-1) = eos_state % rho * eos_state % massfrac(:)
             state(i,j,k,UMX)             = eos_state % rho * u
             state(i,j,k,UMY)             = eos_state % rho * v
             state(i,j,k,UMZ)             = eos_state % rho * w
             state(i,j,k,UEINT)           = eos_state % rho * eint0
             state(i,j,k,UEDEN)           = eos_state % rho * (eint0 + HALF * (u**2 + v**2 + w**2))
             state(i,j,k,UTEMP)           = T0
          enddo
       enddo
    enddo
     
  end subroutine pc_initdata

  subroutine pc_prob_close() &
       bind(C, name="pc_prob_close")
  end subroutine pc_prob_close

end module pc_prob_module
