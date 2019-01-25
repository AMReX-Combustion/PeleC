module pc_prob_module

  implicit none

  private

  public :: amrex_probinit, pc_initdata, pc_prob_close

contains

  subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(C, name = "amrex_probinit")

    use bl_types
    use prob_params_module, only: center, dim
    use probdata_module
    use bl_error_module
    use eos_type_module
    use eos_module
    use network, only : nspec
    use extern_probin_module, only: const_conductivity

    implicit none

    integer init, namlen
    integer name(namlen)
    double precision problo(dim), probhi(dim)

    integer untin,i

    namelist /fortin/ diff_coeff, T1, T2, t_0

    ! Build "probin" filename -- the name of file containing fortin namelist.

    integer, parameter :: maxlen = 256
    character probin*(maxlen)
    real (kind=dp_t) :: Y(nspec)

    type (eos_t) :: eos_state
    call build(eos_state)

    if (namlen .gt. maxlen) call bl_error("probin file name too long")

    do i = 1, namlen
       probin(i:i) = char(name(i))
    end do

    ! Set namelist defaults
    T1 = 1.0_dp_t
    T2 = 2.0_dp_t
    t_0 = 0.001_dp_t
    diff_coeff = 1.0_dp_t

    ! set center, domain extrema
    center(1:dim) = (problo(1:dim)+probhi(1:dim))/2.d0

    ! Read namelists
    untin = 9
    open(untin,file=probin(1:namlen),form='formatted',status='old')
    read(untin,fortin)
    close(unit=untin)

    ! the conductivity is the physical quantity that appears in the
    ! diffusion term of the energy equation.  It is set via
    ! diffusion.conductivity in the inputs file.  For this test problem,
    ! we want to set the diffusion coefficient, D = k/(rho c_v), so the
    ! free parameter we have to play with is rho.  Note that for an
    ! ideal gas, c_v does not depend on rho, so we can call it the EOS
    ! with any density.
    Y(:) = 0.d0
    Y(1) = 1.d0

    eos_state%T = T1
    eos_state%rho = 1.0
    eos_state%massfrac(:) = Y(:)

    call eos_rt(eos_state)

    ! diffusion coefficient is D = k/(rho c_v). we are doing an ideal
    ! gas, so c_v is constant, so find the rho that combines with
    ! the conductivity
    rho0 = const_conductivity/(diff_coeff*eos_state%cv)

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

    use probdata_module, only : T1, T2, diff_coeff, t_0, rho0
    use eos_module
    use network, only: nspec
    use meth_params_module, only : URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, UTEMP
    use prob_params_module, only : problo, center, dim
    use bl_constants_module, only: ZERO, HALF
    
    implicit none

    integer :: level, nvar
    integer :: lo(3), hi(3)
    integer :: state_lo(3),state_hi(3)
    double precision :: xlo(3), xhi(3), time, delta(3)
    double precision :: state(state_lo(1):state_hi(1), &
         state_lo(2):state_hi(2), &
         state_lo(3):state_hi(3),nvar)

    double precision :: xc(3)
    double precision :: Y(nspec), temp
    double precision :: dist2
    integer :: i,j,k,n

    type (eos_t) :: eos_state
    call build(eos_state)

    ! set the composition
    Y(:) = 0.d0
    Y(1) = 1.d0

    do k = lo(3), hi(3)
       xc(3) = problo(3) + delta(3)*(dble(k) + HALF)

       do j = lo(2), hi(2)
          xc(2) = problo(2) + delta(2)*(dble(j) + HALF)

          do i = lo(1), hi(1)
             xc(1) = problo(1) + delta(1)*(dble(i) + HALF)

             state(i,j,k,URHO) = rho0

             dist2 = (xc(1) - center(1))**2
             do n = 2,dim
                dist2 = dist2 + (xc(n) - center(n))**2
             end do

             temp = (T2 - T1)*exp(-0.25*dist2/(diff_coeff*t_0) ) + T1
             state(i,j,k,UTEMP) = temp

             ! compute the internal energy and temperature
             eos_state%T = temp
             eos_state%rho = state(i,j,k,URHO)
             eos_state%massfrac(:) = Y

             call eos_rt(eos_state)

             state(i,j,k,UMX) = ZERO
             state(i,j,k,UMY) = ZERO
             state(i,j,k,UMZ) = ZERO

             state(i,j,k,UEDEN) = rho0*eos_state%e +  &
                  0.5d0*(state(i,j,k,UMX)**2 + &
                  state(i,j,k,UMY)**2 + &
                  state(i,j,k,UMZ)**2)/state(i,j,k,URHO)

             state(i,j,k,UEINT) = rho0*eos_state%e

             state(i,j,k,UFS:UFS-1+nspec) = rho0*Y(:)
          enddo
       enddo
    enddo

  end subroutine pc_initdata

  subroutine pc_prob_close() &
       bind(C, name="pc_prob_close")
  end subroutine pc_prob_close

end module pc_prob_module
