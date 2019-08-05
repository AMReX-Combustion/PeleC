subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use probdata_module
  use amrex_error_module
  implicit none

  integer :: init, namlen
  integer :: name(namlen)
  double precision :: problo(3), probhi(3)

  integer untin,i

  namelist /fortin/ reynolds, mach, prandtl
    
  ! Build "probin" filename -- the name of file containing fortin namelist.
  integer, parameter :: maxlen = 256
  character probin*(maxlen)

  if (namlen .gt. maxlen) then
     call amrex_error('probin file name too long')
  end if

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! set namelist defaults here
  reynolds = 1600.0
  mach = 0.1
  prandtl = 0.71

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

  use probdata_module
  use network, only: nspecies, naux
  use chemistry_module, only : nspecies, get_species_index
  use eos_type_module
  use meth_params_module, only : URHO, UMX, UMY, UMZ, &
       UEDEN, UEINT, UFS, UTEMP, small_temp
  use amrex_constants_module, only: ZERO, HALF, M_PI
  !use extern_probin_module, only: const_viscosity, const_bulk_viscosity, const_conductivity, const_diffusivity
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
  double precision :: L, rho0, v0, p0, T0, a0
  integer, parameter :: out_unit=20
  
  type(eos_t) :: eos_state
  
  integer :: iN2

  if ( nspecies .ne. nspecies ) then
      write(*,*) 'Something very wrong'
      write(*,*) ' nspecies =', nspecies
      write(*,*) ' nspecies =', nspecies
      stop
  endif

  
  iN2 = get_species_index("N2")
  call build(eos_state)

  ! Define the length scale
  L = 1.d0/M_PI

  ! Initial pressure and temperature
  p0 = 1.013d6 ! [erg cm^-3]
  T0 = 400.d0

  ! Set the equation of state variables
  eos_state % p = 1.0e6
  eos_state % T = 400.0
  eos_state % massfrac    = 0.d0
  eos_state % massfrac(iN2) = 1.d0
  call eos_tp(eos_state)
  write(*,*) 'Initial T', eos_state % T
  write(*,*) 'Initial p ', eos_state % p
  write(*,*) 'Initial rho ', eos_state % rho
  write(*,*) 'Initial e ',  eos_state % e

  ! Initial density, velocity, and material properties
  rho0 = eos_state % rho
  a0   = 1d0 
  !const_bulk_viscosity = 0.d0
  !const_diffusivity = 0.d0
  !const_viscosity = rho0 * v0 * L / reynolds
  !const_conductivity = const_viscosity * eos_state % cp / prandtl
  state(:,:,:,UTEMP) = T0

  ! Write this out to file (might be useful for postprocessing)
  open(unit=out_unit,file="ic.txt",action="write",status="replace")
  write(out_unit,*)"L, rho0, v0, p0, T0, gamma, mu, k, c_s0, Reynolds, Mach, Prandtl"
  !write(out_unit,*) L, "," , eos_state % rho, "," , v0, "," , eos_state % p, "," , &
  !     eos_state % T, ",", const_viscosity, "," , &
  !     const_conductivity, "," , eos_state % cs, "," , &
  !     reynolds, "," , mach, "," , prandtl
  close(out_unit)
  
  do k = lo(3), hi(3)
     z = (k+HALF)*delta(3)
     w = 0.d0

     do j = lo(2), hi(2)
        y = xlo(2) + delta(2)*(dble(j-lo(2)) + HALF)
!       y = (j+HALF)*delta(2)

        do i = lo(1), hi(1)
           x = xlo(1) + delta(1)*(dble(i-lo(1)) + HALF)
!          x = (i+HALF)*delta(1)

           u = a0 * y/L
           v = a0 * x/L

           eos_state % p   =  p0 

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
