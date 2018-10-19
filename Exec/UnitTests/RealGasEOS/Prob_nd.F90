module pc_prob_module

  implicit none

  private

  public :: amrex_probinit, pc_initdata, pc_prob_close

contains

  subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(C, name = "amrex_probinit")

    use probdata_module
    use prob_params_module, only: dim
    use bl_error_module
    implicit none

    integer :: init, namlen
    integer :: name(namlen)
    double precision :: problo(dim), probhi(dim)
    integer untin,i

    namelist /fortin/ initP,Tstart,Tend,TestNum,numPoints, delta, OutputFile

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

    initP = 1013250.d0  ! 1 atm

    Tstart = 300.0d0   ! Start temperature variation
    Tend   = 700.0d0   ! End temperature variation

    ! Run all tests
    TestNum = 1

    ! Number of points between Tstart and Tend
    numPoints = 128

    ! For finite difference testing of derivative
    delta = 1.0e-10

    ! Output file name
    OutputFile = 'TestDerivatives_SRK.txt'

    !     Read namelists
    untin = 9
    open(untin,file=probin(1:namlen),form='formatted',status='old')
    read(untin,fortin)
    close(unit=untin)

    call RunEOStests

  end subroutine amrex_probinit

  subroutine pc_prob_close() &
       bind(C, name="pc_prob_close")
  end subroutine pc_prob_close

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
       dx,xlo,xhi) bind(C, name = "pc_initdata")
    use eos_type_module
    use probdata_module
    use meth_params_module, only : URHO, UMX, UMY, UMZ, UEINT, UEDEN, UTEMP, UFS
    use prob_params_module, only : problo, dim
    use eos_module
    use network, only: nspec
    use bl_constants_module, only: M_PI

    implicit none

    integer :: level, nvar
    integer :: lo(3), hi(3)
    integer :: state_lo(3), state_hi(3)
    double precision :: xlo(3), xhi(3), time, dx(3)
    double precision :: state(state_lo(1):state_hi(1), &
         state_lo(2):state_hi(2), &
         state_lo(3):state_hi(3),nvar)

  end subroutine pc_initdata

end module pc_prob_module
