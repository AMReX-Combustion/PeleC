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
    double precision :: problo(3), probhi(3)

    integer untin,i

    namelist /fortin/ p0, u0, phi0, T0_bg, T0_pk, sigma

    integer, parameter :: maxlen = 256
    character probin*(maxlen)

    ! Defaults
    p0     = 1013250.d0
    u0     = 0.d0
    phi0   = 1.d0
    T0_bg  = 300.d0
    T0_pk  = 1200.d0
    sigma  = 2.d0

    !
    ! Build "probin" filename -- the name of file containing fortin namelist.
    !     
    if (namlen .gt. maxlen) then
       write(6,*) 'probin file name too long'
       stop
    end if

    do i = 1, namlen
       probin(i:i) = char(name(i))
    end do

    ! Read namelists
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
    use chemistry_module, only : get_species_index, nspecies
    use prob_params_module, only : problo, probhi, dim
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

    ! local variables
    integer :: i, j, k, iN2, iO2, iH2, d
    double precision :: a, x, y, z, r2, arg, gauss, mid(3), L(3)
    type (eos_t) :: eos_state, eos_state_bg, eos_state_pk

    call build(eos_state)
    call build(eos_state_bg)
    call build(eos_state_pk)

    iN2 = get_species_index("N2")
    iO2 = get_species_index("O2")
    iH2 = get_species_index("H2")

    eos_state_bg % molefrac = 0.d0
    eos_state_bg % T        = T0_bg
    eos_state_bg % p        = p0

    ! hc + a.O2 -> b.CO2 + c.H2O
    ! for hc = H2, a=0.5, b=0, c=1
    a = 0.5d0
    eos_state_bg % molefrac(iO2) = 1.d0/(1.d0 + phi0/a  + 0.79d0/0.21d0)
    eos_state_bg % molefrac(iH2) = phi0 * eos_state_bg % molefrac(iO2) / a
    eos_state_bg % molefrac(iN2) = (0.79d0/0.21d0) * eos_state_bg % molefrac(iO2)

    !eos_state_pk % molefrac = eos_state_bg % molefrac
!!! hack
    eos_state_pk % molefrac = 0.d0
    eos_state_pk % molefrac(iN2) = 0.79d0
    eos_state_pk % molefrac(iO2) = 0.21d0
!!!

    eos_state_pk % T        = T0_pk
    eos_state_pk % p        = p0

    do d=1,dim
       mid(d) = 0.5d0 * (problo(d) + probhi(d))
       L(d) = probhi(d) - problo(d)
    enddo

    do k = lo(3), hi(3)
       z = problo(3) + (k+0.5d0)*dx(3)
       do j = lo(2), hi(2)
          y = problo(2) + (j+0.5d0)*dx(2)
          do i = lo(1), hi(1)
             x = problo(1) + (i+0.5d0)*dx(1)

             ! r2 = (x-mid(1))**2 + (y-mid(2))**2
             ! arg = -r2 / (2.0*sigma**2)
             !if (ABS(arg).gt.50.0) arg=-50.d0
             ! gauss = 1.d0/(sigma * sqrt(2*M_PI)) * exp(arg)
             ! gauss = exp(arg)

             eos_state%molefrac = eos_state_pk%molefrac
             eos_state%p        = eos_state_pk%p
             !eos_state%T        = gauss*eos_state_pk%T        + (1.d0 - gauss)*eos_state_bg%T

!!! hack
             r2 = (x-mid(1))**2 + (y-mid(2))**2
             !r2 = (y-mid(2))**2
             if (sqrt(r2) < L(2)*0.25) then
                eos_state%T = T0_pk
                eos_state = eos_state_pk
             else
                eos_state%T = T0_bg
                eos_state = eos_state_bg
             endif
!!!

             call eos_xty(eos_state)
             call eos_tp(eos_state)

             state(i,j,k,URHO )   = eos_state % rho
             state(i,j,k,UMX:UMZ) = 0.d0
             state(i,j,k,UEINT)   = eos_state % rho  * eos_state % e
             state(i,j,k,UEDEN)   = eos_state % rho  * eos_state % e
             state(i,j,k,UTEMP)   = eos_state % T
             state(i,j,k,UFS:UFS+nspec-1) = eos_state % rho  *  eos_state % massfrac(1:nspec)

          end do
       end do
    end do

    call destroy(eos_state_pk)
    call destroy(eos_state_bg)
    call destroy(eos_state)

  end subroutine pc_initdata

  subroutine pc_prob_close() &
       bind(C, name="pc_prob_close")

    use probdata_module

    ! do nothing

  end subroutine pc_prob_close

end module pc_prob_module
