module pc_prob_module

  implicit none

  private

  public :: amrex_probinit, pc_initdata, pc_prob_close

contains
  
  subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(C, name = "amrex_probinit")
    
    use probdata_module
    use prob_params_module, only: dim
    use amrex_error_module
    implicit none

    integer :: init, namlen
    integer :: name(namlen)
    double precision :: problo(dim), probhi(dim)

    integer untin,i

    namelist /fortin/ p_init, Y_init_H2, Y_init_O2, Y_init_N2, T_init

    !
    !     Build "probin" filename -- the name of file containing fortin namelist.
    !     
    integer, parameter :: maxlen = 256
    character probin*(maxlen)

    if (namlen .gt. maxlen) then
       write(6,*) 'probin file name too long'
       stop
    end if

    do i = 1, namlen
       probin(i:i) = char(name(i))
    end do

    p_init    = 1013250.d0  ! 1 atm
    Y_init_H2 = 0.06d0
    Y_init_O2 = 0.5d0
    Y_init_N2 = 0.44d0
    T_init    = 940.d0

    !     Read namelists
    untin = 9
    open(untin,file=probin(1:namlen),form='formatted',status='old')
    read(untin,fortin)
    close(unit=untin)

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
  ! ::: dx        => cell size
  ! ::: xlo, xhi  => physical locations of lower left and upper
  ! :::              right hand corner of grid.  (does not include
  ! :::		   ghost region).
  ! ::: -----------------------------------------------------------

  subroutine pc_initdata(level,time,lo,hi,nvar, &
       state,state_lo,state_hi, &
       dx,xlo,xhi) bind(C, name="pc_initdata")

    use eos_type_module
    use probdata_module
    use meth_params_module, only : URHO, UMX, UMY, UMZ, UEINT, UEDEN, UTEMP, UFS
    use eos_module
    use network, only: nspecies
    use chemistry_module, only : get_species_index
    use amrex_constants_module, only: M_PI, HALF

    implicit none

    integer :: level, nvar
    integer :: lo(3), hi(3)
    integer :: state_lo(3), state_hi(3)
    double precision :: xlo(3), xhi(3), time, dx(3)
    double precision :: state(state_lo(1):state_hi(1), &
         state_lo(2):state_hi(2), &
         state_lo(3):state_hi(3),nvar)

    ! local variables
    integer :: iH2,iO2,iN2
    integer :: i,j,k
    
    type (eos_t) :: eos_state
    call build(eos_state)
    
    iH2   = get_species_index("H2")
    iO2   = get_species_index("O2")
    iN2   = get_species_index("N2")

    eos_state % p = p_init
    eos_state % T = T_init
    eos_state % massfrac(1:nspecies) = 0.d0
    eos_state % massfrac(iH2)     = Y_init_H2
    eos_state % massfrac(iO2)     = Y_init_O2
    eos_state % massfrac(iN2)     = Y_init_N2
    call eos_tp(eos_state)

    do k=lo(3),hi(3)
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)

                state(i,j,k,URHO ) = eos_state % rho
                state(i,j,k,UMX  ) = 0.d0
                state(i,j,k,UMY  ) = 0.d0
                state(i,j,k,UMZ  ) = 0.d0
                state(i,j,k,UEINT) = eos_state % rho  *  eos_state % e
                state(i,j,k,UEDEN) = eos_state % rho  *  eos_state % e
                state(i,j,k,UTEMP) = eos_state % T
                state(i,j,k,UFS:UFS+nspecies-1) = eos_state % rho  *  eos_state % massfrac(1:nspecies)

           enddo
        enddo
     enddo

    call destroy(eos_state)

  end subroutine pc_initdata

  subroutine pc_prob_close() &
       bind(C, name="pc_prob_close")

    use probdata_module

  end subroutine pc_prob_close

end module pc_prob_module
