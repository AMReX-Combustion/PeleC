module pc_prob_module

  implicit none

  private

  public :: amrex_probinit, pc_initdata, pc_prob_close

contains
  
  subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(C, name = "amrex_probinit")
    
    use probdata_module
    use amrex_fort_module

    implicit none
    integer :: init, namlen
    integer :: name(namlen)
    double precision :: problo(3), probhi(3)

    integer :: untin,i

    namelist /fortin/ p_domain, dens_domain, vx_in, vy_in

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

    p_domain    = 1.d0        ! initial conditions - pressure
    dens_domain = 1.d0        ! initial conditions - density
    vx_in       = 3.16227     ! vx at inflow
    vy_in       = 0.d0        ! vy at inflow

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
    double precision :: xx, yy, zz, xl, xr
    double precision :: dist
    double precision :: vctr, p_exp, vol_ambient, vol_pert, dx_sub

    integer :: i,j,k, ii, jj, kk
    integer :: npert, nambient

    type (eos_t) :: eos_state

    call build(eos_state)

    eos_state % rho = dens_domain
    eos_state % p   = p_domain
    eos_state % T   = 100000.0
    eos_state % massfrac = ZERO
    eos_state % massfrac(nspecies) = ONE
    call eos_rp(eos_state)

    print *,"rhoe:",eos_state%rho * eos_state%e

    do k = lo(3), hi(3)
       zmin = xlo(3) + delta(3)*dble(k-lo(3)) 

        do j = lo(2), hi(2)
           ymin = xlo(2) + delta(2)*dble(j-lo(2))

           do i = lo(1), hi(1)
              xmin = xlo(1) + delta(1)*dble(i-lo(1))

                state(i,j,k,URHO) = eos_state % rho
                state(i,j,k,UMX) = dens_domain * vx_in
                state(i,j,k,UMY) = dens_domain * vy_in
                state(i,j,k,UMZ) = 0.d0
                state(i,j,k,UEINT) = eos_state % rho * eos_state % e
                state(i,j,k,UEDEN) = state(i,j,k,UEINT) +  &
                     0.5d0*dens_domain*(vx_in**2 + vy_in**2)
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
