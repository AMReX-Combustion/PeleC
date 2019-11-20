module pc_prob_module

  implicit none

  private

  public :: amrex_probinit, pc_initdata, pc_prob_close

contains
  
  subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(C, name = "amrex_probinit")
    
    use probdata_module
    use amrex_error_module
    implicit none

    integer :: init, namlen
    integer :: name(namlen)
    double precision :: problo(3), probhi(3)

    integer untin,i

    namelist /fortin/ p_ref, u_ref, T_ref

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
    p_ref = 101300.0d0
    u_ref = 100.0d0
    T_ref = 300.0d0

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
    use network, only: nspecies, naux, molec_wt
    use meth_params_module, only : URHO, UMX, UMY, UMZ, UEINT, UEDEN, UTEMP, UFS
    use prob_params_module, only : problo, probhi
    use eos_module
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
    integer :: i, j, k, n
    double precision :: x, y ,z, xysqr, xcen, ycen, delta
    double precision :: p0, T0, c_0, RC, RCsqr, C
    type (eos_t) :: eos_state

    call build(eos_state)

    xcen = (probhi(1)/2.0d0)
    ycen = (probhi(2)/2.0d0)

    eos_state % p = p_ref
    eos_state % T = T_ref
    eos_state % massfrac(1) = 1.d0
    molec_wt = 289703.374673617 

    call eos_tp(eos_state)
    c_0 = eos_state % cs

    !  c_0 = sqrt(gamma*287.84956*300.0d0)
    delta = 0.001
    RC = 0.1 * probhi(1)
    RCsqr = RC**2.0
    C =  0.11


    do k = lo(3), hi(3)
       z = problo(3) + dx(3)*(dble(k) + HALF)

       do j = lo(2), hi(2)
          y = problo(2) + dx(2)*(dble(j) + HALF)

          do i = lo(1), hi(1)
             x = problo(1) + dx(1)*(dble(i) + HALF)

             xysqr = (((x-xcen)**2.0d0) + ((y-ycen)**2.0d0))

             if (sqrt(xysqr) <= (4.0d0*RC)) then
                eos_state % p = p_ref*exp(-0.5d0*gamma_const*((C/(c_0*RC))**2.0d0)*exp(-xysqr/RCsqr))
                eos_state % T = T_ref            
                call eos_tp(eos_state)
                state(i,j,k,URHO )  = eos_state % rho

                state(i,j,k,UMX  )  = state(i,j,k,URHO ) *(u_ref +   ((y-ycen)*C/(RCsqr))*exp(-xysqr/(2.0d0*RCsqr)) )
                state(i,j,k,UMY  )  = -state(i,j,k,URHO )*( ((x-xcen)*C/(RCsqr))*exp(-xysqr/(2.0d0*RCsqr)))

             else
                eos_state % p = p_ref
                eos_state % T = T_ref           
                call eos_tp(eos_state)
                state(i,j,k,URHO )  = eos_state % rho

                state(i,j,k,UMX  )  = state(i,j,k,URHO ) * u_ref
                state(i,j,k,UMY  )  = 0.0d0

             end if


             state(i,j,k,UMZ  ) = 0.0 
             state(i,j,k,UEINT) = eos_state % rho  *  eos_state % e
             state(i,j,k,UEDEN) = eos_state % rho * (eos_state % e + 0.5d0 * (state(i,j,k,UMX)**2 + state(i,j,k,UMY)**2))
             state(i,j,k,UTEMP) = eos_state % T
             do n=1, nspecies
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
