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
    double precision :: problo(3), probhi(3)

    integer untin,i

    namelist /fortin/ frac, &
         rho_1, rho_2, p0_base

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
    frac = 0.5d0
    rho_1 = 1.0d0
    rho_2 = 2.0d0
    p0_base = 5.0d0

    ! Read namelists
    untin = 9
    open(untin,file=probin(1:namlen),form='formatted',status='old')
    read(untin,fortin)
    close(unit=untin)


    ! set local variable defaults
    do i = 1,dim
       split(i) = frac*(problo(i)+probhi(i))
    end do

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

    use eos_type_module
    use eos_module
    use probdata_module
    use prob_params_module, only: dim
    use meth_params_module, only : URHO, UMX, UMY, UMZ, &
         UEDEN, UEINT, UFS, UTEMP, small_temp
    use bl_constants_module, only: ZERO, HALF, M_PI

    implicit none

    integer :: level, nvar
    integer :: lo(3), hi(3)
    integer :: state_lo(3),state_hi(3)
    double precision :: xlo(3), xhi(3), time, delta(3)
    double precision :: state(state_lo(1):state_hi(1), &
         state_lo(2):state_hi(2), &
         state_lo(3):state_hi(3),nvar)

    integer :: i,j,k
    double precision :: x(3),r2d,pres,presmid,pertheight

    type (eos_t) :: eos_state

    call build(eos_state)

    presmid  = p0_base - rho_1*split(dim)


    do k = lo(3), hi(3)
       x(3) = (k+HALF)*delta(3)

       do j = lo(2), hi(2)
          x(2) = (j+HALF)*delta(2)

          do i = lo(1), hi(1)
             x(1) = (i+HALF)*delta(1)

             if (dim .eq.3) then
                r2d = min(sqrt((x(1)-split(1))**2+(x(2)-split(2))**2), 0.5d0*L_x)
                pertheight = 0.5d0 - 0.01d0*cos(2.0d0*M_PI*r2d/L_x)
             else
                pertheight = 0.01d0*HALF*(cos(2.0d0*M_PI*x(1)/L_x) + &
                     cos(2.0d0*M_PI*(L_x-x(1))/L_x)) + 0.5d0
             end if
             state(i,j,k,URHO) = rho_1 + ((rho_2-rho_1)/2.0d0)* &
                  (1+tanh((x(dim)-pertheight)/0.005d0))
             state(i,j,k,UFS) = state(i,j,k,URHO)
             state(i,j,k,UFS+1:UFS+nspec-1) = 0.d0

             eos_state%rho = state(i,j,k,URHO)
             eos_state%T = 100000.d0  ! initial guess
             eos_state%massfrac(:) = state(i,j,k,UFS:UFS+nspec-1)

             if (x(dim) .lt. split(dim)) then
                eos_state%p = p0_base - rho_1*x(dim)
             else
                eos_state%p = presmid - rho_2*(x(dim)-split(dim))
             end if

             call eos_rp(eos_state)

             state(i,j,k,UEINT) = eos_state % e * eos_state % rho
             state(i,j,k,UEDEN) = state(i,j,k,UEINT) &
                  + HALF * (state(i,j,k,UMX)**2 + state(i,j,k,UMY)**2 + state(i,j,k,UMZ)**2)/eos_state%rho

          enddo
       enddo
    enddo

    call destroy(eos_state)

  end subroutine pc_initdata

  subroutine pc_prob_close() &
       bind(C, name="pc_prob_close")
  end subroutine pc_prob_close

end module pc_prob_module
