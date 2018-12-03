module pc_prob_module

  implicit none

  private

  public :: amrex_probinit, pc_initdata, pc_prob_close

contains
  
  subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(C, name = "amrex_probinit")
    
    use probdata_module
    use bl_error_module
    use chemistry_module, only : nspecies, get_species_index

    implicit none
    integer :: init, namlen
    integer :: name(namlen)
    double precision :: problo(3), probhi(3)

    integer :: untin,i

    namelist /fortin/ Pres_domain, Temp_domain, Yfuel_domain, Yox_domain, &
    YN2_domain, dens_jet, vel_jet, Yox_jet, Yfuel_jet, YN2_jet, centx, centz, &
    r_circ, r_hole, nholes, cone_angle 

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
    Pres_domain    = 10132500.0
    Temp_domain    = 500.0
    Yfuel_domain   = 0.d0
    Yox_domain     = 0.25d0
    YN2_domain     = 0.75d0
    dens_jet       = 0.0001d0
    vel_jet        = 5000.d0
    Yfuel_jet      = 1.d0
    Yox_jet        = 0.d0
    YN2_jet        = 0.d0
    centx          = 2.5d0
    centz          = 2.5d0
    r_circ         = 3.d0
    r_hole         = 0.4d0
    nholes         = 4
    cone_angle     = 20.d0

    ! Read namelists
    untin = 9
    open(untin,file=probin(1:namlen),form='formatted',status='old')
    read(untin,fortin)
    close(unit=untin)

    fuel_id = get_species_index("CH4")
    ox_id   = get_species_index("O2")
    N2_id   = get_species_index("N2") 


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
    use network, only: nspec
    use probdata_module
    use meth_params_module, only : URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, UFS
    use bl_constants_module, only: M_PI, FOUR3RD, ONE, HALF
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

    type (eos_t) :: eos_state
    integer :: i,j,k

    call build(eos_state)

    eos_state % p   = Pres_domain
    eos_state % T   = Temp_domain
    eos_state % massfrac = ZERO
    eos_state % massfrac(fuel_id) = Yfuel_domain
    eos_state % massfrac(ox_id)   = Yox_domain
    eos_state % massfrac(N2_id)   = YN2_domain

    call eos_tp(eos_state)

    do k = lo(3), hi(3)

        do j = lo(2), hi(2)

           do i = lo(1), hi(1)

                state(i,j,k,URHO) = eos_state % rho
                state(i,j,k,UMX) = 0.d0
                state(i,j,k,UMY) = 0.d0
                state(i,j,k,UMZ) = 0.d0
                state(i,j,k,UEINT) = eos_state % rho * eos_state % e
                state(i,j,k,UEDEN) = state(i,j,k,UEINT)
                state(i,j,k,UFS:UFS+nspec-1) = eos_state % massfrac(:) * state(i,j,k,URHO)
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
