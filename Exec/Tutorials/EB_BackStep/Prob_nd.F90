module pc_prob_module

  implicit none

  private

  public :: amrex_probinit, pc_initdata, pc_prob_close

contains
  
  subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(C, name = "amrex_probinit")
    
    use probdata_module
    use amrex_fort_module
    use extern_probin_module, only: const_viscosity, const_bulk_viscosity, const_conductivity, const_diffusivity
    use eos_module

    implicit none
    integer :: init, namlen
    integer :: name(namlen)
    double precision :: problo(3), probhi(3)

    integer :: untin,i

    type(eos_t) :: eos_state
    namelist /fortin/ p_domain, dens_domain, vx_in, vy_in, Re_L, Pr

    ! Build "probin" filename -- the name of file containing fortin namelist.
    integer, parameter :: maxlen = 256
    character :: probin*(maxlen)
    double precision :: L

    if (namlen .gt. maxlen) then
       call bl_error('probin file name too long')
    end if

    do i = 1, namlen
       probin(i:i) = char(name(i))
    end do

    ! set namelist defaults

    p_domain    = 1013250.0        ! initial conditions - pressure
    dens_domain = 0.00116      ! initial conditions - density
    vx_in       = 5000.0     ! vx at inflow
    vy_in       = 0.d0        ! vy at inflow
    Re_L        = 10000.d0
    Pr          = 0.7
   
    !two third of domain 
    L = (probhi(1) - problo(1))*0.66 

    ! Read namelists
    untin = 9
    open(untin,file=probin(1:namlen),form='formatted',status='old')
    read(untin,fortin)
    close(unit=untin)

    !default values
    call build(eos_state)
    eos_state % p = p_domain
    eos_state % rho = dens_domain
    eos_state % massfrac = 0.d0
    eos_state % massfrac(1) = 1.d0
    call eos_rp(eos_state)

    !compute transport from Re and Pr
    const_viscosity    = dens_domain*vx_in*L/Re_L
    const_conductivity = const_viscosity*eos_state%cp/Pr


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
    use network, only: nspecies
    use probdata_module
    use meth_params_module, only : URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, UFS
    use amrex_constants_module, only: M_PI, FOUR3RD, ONE, HALF, ZERO
    use eos_type_module
    use eos_module
    implicit none

    integer :: level, nvar
    integer :: lo(3), hi(3)
    integer :: state_lo(3), state_hi(3)
    double precision :: xlo(3), xhi(3), time, delta(3)
    double precision :: state(state_lo(1):state_hi(1), &
         state_lo(2):state_hi(2), &
         state_lo(3):state_hi(3),nvar)


    double precision :: xmin,ymin,zmin

    integer :: i,j,k

    type (eos_t) :: eos_state

    call build(eos_state)

    eos_state % rho = dens_domain
    eos_state % p   = p_domain
    eos_state % T   = 300.0
    eos_state % massfrac = ZERO
    eos_state % massfrac(nspecies) = ONE
    call eos_rp(eos_state)


    do k = lo(3), hi(3)
       zmin = xlo(3) + delta(3)*dble(k-lo(3)) 

        do j = lo(2), hi(2)
           ymin = xlo(2) + delta(2)*dble(j-lo(2))

           do i = lo(1), hi(1)
              xmin = xlo(1) + delta(1)*dble(i-lo(1))

                state(i,j,k,URHO) = eos_state % rho
                state(i,j,k,UMX) = 0.d0
                state(i,j,k,UMY) = 0.d0
                state(i,j,k,UMZ) = 0.d0
                state(i,j,k,UEINT) = eos_state % rho * eos_state % e
                state(i,j,k,UEDEN) = state(i,j,k,UEINT)
                state(i,j,k,UFS:UFS+nspecies-1) = eos_state % massfrac(:) * state(i,j,k,URHO)
                state(i,j,k,UTEMP) = eos_state % T

             enddo
          enddo
       enddo

    call destroy(eos_state)

  end subroutine pc_initdata

  subroutine pc_prob_close() &
       bind(C, name="pc_prob_close")
  end subroutine pc_prob_close

end module pc_prob_module
