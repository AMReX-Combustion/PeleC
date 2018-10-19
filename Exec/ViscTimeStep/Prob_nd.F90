module pc_prob_module

  implicit none

  private

  public :: amrex_probinit, pc_initdata, pc_prob_close, test_pc_estdt_vel_diffusion

contains
  
  subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(C, name = "amrex_probinit")
    
    use probdata_module
    use bl_error_module
    implicit none

    integer :: init, namlen
    integer :: name(namlen)
    double precision :: problo(3), probhi(3)

    integer untin,i

    namelist /fortin/ tolerance

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
    tolerance = 1.d-13

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
       delta,xlo,xhi) bind(C, name = "pc_initdata")

    use probdata_module

    implicit none

    integer :: level, nvar
    integer :: lo(3), hi(3)
    integer :: state_lo(3),state_hi(3)
    double precision :: xlo(3), xhi(3), time, delta(3)
    double precision :: state(state_lo(1):state_hi(1), &
         state_lo(2):state_hi(2), &
         state_lo(3):state_hi(3),nvar)

    integer num, numtests, count_fail
    logical, allocatable :: test_flags(:)

    numtests = 1
    allocate(test_flags(numtests))
    test_flags(:) = .true.
    count_fail = 0

    write(*,'("========================================================")')
    write(*,'("Start of test suite")')

    do num = 1, numtests

       write(*,'(/,"Running test",I2)')num

       ! First test
       if (num == 1) then
          call test_pc_estdt_vel_diffusion(tolerance,test_flags(num))
       endif

       ! Output test result
       if (test_flags(num)) then
          write(*,'("....Test passed")')
       else
          write(*,'("....Test failed")')
          count_fail = count_fail + 1
       endif

    enddo

    ! Summary
    write(*,'(/,"Test suite summary:")')
    if (count_fail > 0) then
       write(*,'(I2," /",I2, " tests failed")')count_fail, numtests
       write(*,'("========================================================")')
       stop 9
    else
       write(*,'("All tests passed")')
       write(*,'("========================================================")')
       stop 0
    endif

  end subroutine pc_initdata

  subroutine pc_prob_close() &
       bind(C, name="pc_prob_close")
  end subroutine pc_prob_close

  subroutine test_pc_estdt_vel_diffusion(tolerance,flag) bind(C, name = "test_pc_estdt_vel_diffusion")

    use timestep_module
    use meth_params_module, only: NVAR, URHO
    use extern_probin_module, only: const_viscosity

    implicit none

    double precision, intent(in) :: tolerance
    logical, intent(inout) :: flag

    integer :: i, j, k
    integer          :: lo(3), hi(3)
    integer          :: s_lo(3), s_hi(3)
    double precision, allocatable :: state(:,:,:,:)
    double precision :: dx(3), dt, rho0

    ! Initialize data
    lo(1) = 1
    lo(2) = 1
    lo(3) = 1
    hi(1) = 2
    hi(2) = 3
    hi(3) = 3
    s_lo(1) = 1
    s_lo(2) = 1
    s_lo(3) = 1
    s_hi(1) = 2
    s_hi(2) = 3
    s_hi(3) = 3
    allocate(state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR))
    dx(:) = 1.d-1
    dt = 1.d200

    ! Dummy data
    const_viscosity = 3.4d0
    rho0 = 1.d0

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          state(1,j,k,URHO) = 2 * rho0
          state(2,j,k,URHO) = 1.d0
          state(3,j,k,URHO) = 3 * rho0
       enddo
    enddo

    ! Call the function to test
    call pc_estdt_vel_diffusion(lo,hi,state,s_lo,s_hi,dx,dt)

    ! Comparison
    if (abs(dt - 0.5*(dx(1)**2)*rho0/const_viscosity) > tolerance) then
       flag = .false.
    endif

  end subroutine test_pc_estdt_vel_diffusion

end module pc_prob_module
