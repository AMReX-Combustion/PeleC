module pc_prob_module

  implicit none

  private

  public :: amrex_probinit, pc_initdata, pc_prob_close

contains
  
  subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(C, name = "amrex_probinit")
    
    use eos_module
    use eos_type_module
    use bl_error_module 
    use network
    use probdata_module

    implicit none

    integer init, namlen
    integer name(namlen)
    double precision problo(3), probhi(3)
    double precision massfrac(nspec)

    integer untin,i

    type (eos_t) :: eos_state

    namelist /fortin/ p_l, u_l, rho_l, p_r, u_r, rho_r, frac, idir

    !
    !     Build "probin" filename -- the name of file containing fortin namelist.
    !     
    integer maxlen
    parameter (maxlen=256)
    character probin*(maxlen)

    call build(eos_state)

    if (namlen .gt. maxlen) then
       call bl_error("probin file name too long")
    end if

    do i = 1, namlen
       probin(i:i) = char(name(i))
    end do

    ! set namelist defaults

    p_l = 10132500.d0       ! left pressure (erg/cc)
    u_l = 0.0               ! left velocity (cm/s)
    rho_l = 1.0             ! left density (g/cc)

    p_r = 0.1*p_l           ! right pressure (erg/cc)
    u_r = 0.0               ! right velocity (cm/s)
    rho_r = 0.1*rho_l       ! right density (g/cc)

    idir = 1                ! direction across which to jump
    frac = 0.5              ! fraction of the domain for the interface

    !     Read namelists
    untin = 9
    open(untin,file=probin(1:namlen),form='formatted',status='old')
    read(untin,fortin)
    close(unit=untin)

    split(1) = frac*(problo(1)+probhi(1))
    split(2) = frac*(problo(2)+probhi(2))
    split(3) = frac*(problo(3)+probhi(3))

    ! compute the internal energy (erg/cc) for the left and right state
    massfrac(:) = 0.0d0
    massfrac(1) = 1.0d0

    eos_state%rho = rho_l
    eos_state%p = p_l
    eos_state%T = 100000.d0  ! initial guess
    eos_state%massfrac(:) = massfrac(:)

    call eos_rp(eos_state)

    rhoe_l = rho_l*eos_state%e
    T_l = eos_state%T

    eos_state%rho = rho_r
    eos_state%p = p_r
    eos_state%T = 100000.d0  ! initial guess
    eos_state%massfrac(:) = massfrac(:)

    call eos_rp(eos_state)

    rhoe_r = rho_r*eos_state%e
    T_r = eos_state%T

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
    use prob_params_module, only: dim
    implicit none

    integer :: level, nvar
    integer :: lo(3), hi(3)
    integer :: state_lo(3), state_hi(3)
    double precision :: xlo(3), xhi(3), time, delta(3)
    double precision :: state(state_lo(1):state_hi(1), &
         state_lo(2):state_hi(2), &
         state_lo(3):state_hi(3),nvar)

    double precision xcen,ycen,zcen
    integer i,j,k

    if (idir .gt. dim .or. idir .lt.1) then
       call bl_abort('invalid idir')
    end if

    do k = lo(3), hi(3)
       zcen = xlo(3) + delta(3)*(float(k-lo(3)) + 0.5d0)

       do j = lo(2), hi(2)
          ycen = xlo(2) + delta(2)*(float(j-lo(2)) + 0.5d0)

          do i = lo(1), hi(1)
             xcen = xlo(1) + delta(1)*(float(i-lo(1)) + 0.5d0)

             if (idir == 1) then
                if (xcen <= split(1)) then
                   state(i,j,k,URHO) = rho_l
                   state(i,j,k,UMX) = rho_l*u_l
                   state(i,j,k,UMY) = 0.d0
                   state(i,j,k,UMZ) = 0.d0
                   state(i,j,k,UEDEN) = rhoe_l + 0.5*rho_l*u_l*u_l
                   state(i,j,k,UEINT) = rhoe_l
                   state(i,j,k,UTEMP) = T_l
                else
                   state(i,j,k,URHO) = rho_r
                   state(i,j,k,UMX) = rho_r*u_r
                   state(i,j,k,UMY) = 0.d0
                   state(i,j,k,UMZ) = 0.d0
                   state(i,j,k,UEDEN) = rhoe_r + 0.5*rho_r*u_r*u_r
                   state(i,j,k,UEINT) = rhoe_r
                   state(i,j,k,UTEMP) = T_r
                endif

             else if (idir == 2) then
                if (ycen <= split(2)) then
                   state(i,j,k,URHO) = rho_l
                   state(i,j,k,UMX) = 0.d0
                   state(i,j,k,UMY) = rho_l*u_l
                   state(i,j,k,UMZ) = 0.d0
                   state(i,j,k,UEDEN) = rhoe_l + 0.5*rho_l*u_l*u_l
                   state(i,j,k,UEINT) = rhoe_l
                   state(i,j,k,UTEMP) = T_l
                else
                   state(i,j,k,URHO) = rho_r
                   state(i,j,k,UMX) = 0.d0
                   state(i,j,k,UMY) = rho_r*u_r
                   state(i,j,k,UMZ) = 0.d0
                   state(i,j,k,UEDEN) = rhoe_r + 0.5*rho_r*u_r*u_r
                   state(i,j,k,UEINT) = rhoe_r
                   state(i,j,k,UTEMP) = T_r
                endif

             else if (idir == 3) then
                if (zcen <= split(3)) then
                   state(i,j,k,URHO) = rho_l
                   state(i,j,k,UMX) = 0.d0
                   state(i,j,k,UMY) = 0.d0
                   state(i,j,k,UMZ) = rho_l*u_l
                   state(i,j,k,UEDEN) = rhoe_l + 0.5*rho_l*u_l*u_l
                   state(i,j,k,UEINT) = rhoe_l
                   state(i,j,k,UTEMP) = T_l
                else
                   state(i,j,k,URHO) = rho_r
                   state(i,j,k,UMX) = 0.d0
                   state(i,j,k,UMY) = 0.d0
                   state(i,j,k,UMZ) = rho_r*u_r
                   state(i,j,k,UEDEN) = rhoe_r + 0.5*rho_r*u_r*u_r
                   state(i,j,k,UEINT) = rhoe_r
                   state(i,j,k,UTEMP) = T_r
                endif
             endif

             state(i,j,k,UFS:UFS-1+nspec) = 0.0d0
             state(i,j,k,UFS  ) = state(i,j,k,URHO)


          enddo
       enddo
    enddo

  end subroutine pc_initdata

  subroutine pc_prob_close() &
       bind(C, name="pc_prob_close")
  end subroutine pc_prob_close

end module pc_prob_module
