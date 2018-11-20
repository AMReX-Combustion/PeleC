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
    double precision :: problo(2), probhi(2)

    integer untin,i

    namelist /fortin/ prob_type, turbfile, pamb, phi_in, T_in, vn_in, T_co, vn_co, &
         splitx, xfrontw, Tfrontw, blobr, blobx, bloby, blobT, inflow_period, inflow_vnmag, &
         splity, yfrontw, turb_boost_factor, &
         max_tracerr_lev, tracerr, max_vorterr_lev, vorterr, max_tempgrad_lev, tempgrad

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

    prob_type = 0

    pamb = 4.053d7  ! 40 Patm

    phi_in = 0.2d0   ! mole fraction of CH3OCH3 in fuel
    T_in   = 400.d0  ! temperature of fuel
    vn_in  = 5.12d3  ! fuel injection velocity

    T_co  = 1525.d0   ! temperature of air
    vn_co = 5.12d2    ! air injection velocity

    splitx  = 0.00569d0  ! where fuel and air split
    xfrontw = .2d0       ! controls the width of split
    Tfrontw = 0.0075d0

    inflow_period = 1.d-5 ! period of sinusoidal variation of inflow velocity
    inflow_vnmag  = 1.d3  ! magnitude of the variation

    blobr = -1.d0
    blobx = 0.d0
    bloby = 0.027d0
    blobT = 1500.d0

    splity  = 0.001d0
    yfrontw = 0.0004d0

    turb_boost_factor = 1.d0

    max_tracerr_lev = -1
    tracerr = 1.d-8

    max_vorterr_lev = -1
    vorterr = 5.d4

    max_tempgrad_lev = -1
    tempgrad = 10.d0

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
    use prob_params_module, only : problo
    use eos_module

    implicit none

    integer :: level, nvar
    integer :: lo(3), hi(3)
    integer :: state_lo(3), state_hi(3)
    double precision :: xlo(3), xhi(3), time, dx(3)
    double precision :: state(state_lo(1):state_hi(1), &
         state_lo(2):state_hi(2), &
         state_lo(3):state_hi(3),nvar)

    ! local variables
    integer :: i, j, k, n
    double precision :: x(3), u(3)
    type (eos_t) :: eos_state

    call build(eos_state)

    if (.not. jet_initialized) then
       call init_jet()
    end if

    do k = lo(3), hi(3)
       x(3) = problo(3) + dx(3)*(dble(k) + HALF)

       do j = lo(2), hi(2)
          x(2) = problo(2) + dx(2)*(dble(j) + HALF)

          do i = lo(1), hi(1)
             x(1) = problo(1) + dx(1)*(dble(i) + HALF)

             call inflow_boundary(x, eos_state, u)

             state(i,j,k,URHO ) = eos_state % rho
             state(i,j,k,UMX  ) = eos_state % rho  *  u(1)
             state(i,j,k,UMY  ) = eos_state % rho  *  u(2)
             state(i,j,k,UMZ  ) = eos_state % rho  *  u(3)
             state(i,j,k,UEINT) = eos_state % rho  *  eos_state % e
             state(i,j,k,UEDEN) = eos_state % rho * (eos_state % e + 0.5d0 * (u(1)**2 + u(2)**2 + u(3)**2))
             state(i,j,k,UTEMP) = eos_state % T
             do n=1, nspec
                state(i,j,k,UFS+n-1) = eos_state % rho  *  eos_state % massfrac(n)
             end do

          end do
       end do
    end do

  end subroutine pc_initdata

  subroutine pc_prob_close() &
       bind(C, name="pc_prob_close")
  end subroutine pc_prob_close

end module pc_prob_module
