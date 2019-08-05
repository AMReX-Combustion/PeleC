module pc_prob_module

  implicit none

  private

  public :: amrex_probinit, pc_initdata, pc_prob_close

contains
  
  subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(C, name = "amrex_probinit")
    
    use probdata_module
    use amrex_fort_module
    use amrex_constants_module, only: M_PI, HALF
    use eos_module

    implicit none
    integer :: init, namlen
    integer :: name(namlen)
    double precision :: problo(3), probhi(3)

    integer :: untin,i

    type(eos_t) :: eos_state
    namelist /fortin/ pl,rhol,pr,rhor,angle

    ! Build "probin" filename -- the name of file containing fortin namelist.
    integer, parameter :: maxlen = 256
    character :: probin*(maxlen)
    double precision :: cost, sint, lo_xd,hi_xd

    if (namlen .gt. maxlen) then
       call bl_error('probin file name too long')
    end if

    do i = 1, namlen
       probin(i:i) = char(name(i))
    end do

    pl=1.d0
    rhol=1.d0
    pr=0.1d0
    rhor=0.125d0
    angle=0.d0

    ! Read namelists
    untin = 9
    open(untin,file=probin(1:namlen),form='formatted',status='old')
    read(untin,fortin)
    close(unit=untin)
    
    cost = cos(M_PI/180.d0 * angle)
    sint = sin(M_PI/180.d0 * angle)


    lo_xd=problo(1)*cost+problo(2)*sint
    hi_xd=probhi(1)*cost+probhi(2)*sint

    midx=0.5d0*(lo_xd+hi_xd)


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
    use network, only: nspecies, molec_wt
    use probdata_module
    use meth_params_module, only : URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, UFS
    use amrex_constants_module, only: M_PI, HALF
    use eos_type_module
    use eos_module
    implicit none

    integer :: level, nvar
    integer :: lo(3), hi(3)
    integer :: state_lo(3),state_hi(3)
    double precision :: xlo(3), xhi(3), time, delta(3)
    double precision :: state(state_lo(1):state_hi(1), &
         state_lo(2):state_hi(2), &
         state_lo(3):state_hi(3),nvar)

    integer :: i,j,k
    double precision :: x, y, z, cost, sint, xp
    type(eos_t) :: eos_state_l,eos_state_r

    call build(eos_state_l)
    call build(eos_state_r)

    ! Define the molecular weight for air
    molec_wt = 28.97

    ! Set the equation of state variables
    eos_state_l % p    = pl
    eos_state_l % rho = rhol
    eos_state_l % massfrac    = 0.d0
    eos_state_l % massfrac(1) = 1.d0
    call eos_rp(eos_state_l)
    
    ! Set the equation of state variables
    eos_state_r % p    = pr
    eos_state_r % rho = rhor
    eos_state_r % massfrac    = 0.d0
    eos_state_r % massfrac(1) = 1.d0
    call eos_rp(eos_state_r)

    ! Initial density, velocity, and material properties
    cost = cos(M_PI/180.d0 * angle)
    sint = sin(M_PI/180.d0 * angle)

    do k = lo(3), hi(3)
       z = xlo(3) + delta(3)*(dble(k-lo(3)) + HALF)
       do j = lo(2), hi(2)
          y = xlo(2) + delta(2)*(dble(j-lo(2)) + HALF)
          do i = lo(1), hi(1)
             x = xlo(1) + delta(1)*(dble(i-lo(1)) + HALF)

             xp = x * cost + y * sint

             if(xp .lt. midx) then 
                state(i,j,k,URHO)            = eos_state_l % rho
                state(i,j,k,UFS:UFS+nspecies-1) = eos_state_l % rho * eos_state_l % massfrac(:)
                state(i,j,k,UMX)             = 0.d0
                state(i,j,k,UMY)             = 0.d0
                state(i,j,k,UMZ)             = 0.d0
                state(i,j,k,UEINT)           = eos_state_l % rho * eos_state_l % e
                state(i,j,k,UEDEN)           = eos_state_l % rho * eos_state_l % e
                state(i,j,k,UTEMP)           = eos_state_l % T
             else
                state(i,j,k,URHO)            = eos_state_r % rho
                state(i,j,k,UFS:UFS+nspecies-1) = eos_state_r % rho * eos_state_r % massfrac(:)
                state(i,j,k,UMX)             = 0.d0
                state(i,j,k,UMY)             = 0.d0
                state(i,j,k,UMZ)             = 0.d0
                state(i,j,k,UEINT)           = eos_state_r % rho * eos_state_r % e
                state(i,j,k,UEDEN)           = eos_state_r % rho * eos_state_r % e
                state(i,j,k,UTEMP)           = eos_state_r % T
             endif

          enddo
       enddo
    enddo

    call destroy(eos_state_l)
    call destroy(eos_state_r)
     
  end subroutine pc_initdata

  subroutine pc_prob_close() &
       bind(C, name="pc_prob_close")
  end subroutine pc_prob_close

end module pc_prob_module
