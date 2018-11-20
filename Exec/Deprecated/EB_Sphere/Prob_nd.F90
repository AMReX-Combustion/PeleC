subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

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

  namelist /fortin/ rpulse,rho0,drho0,p0

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
! ::: delta     => cell size
! ::: xlo, xhi  => physical locations of lower left and upper
! :::              right hand corner of grid.  (does not include
! :::		   ghost region).
! ::: -----------------------------------------------------------

subroutine pc_initdata(level,time,lo,hi,nvar, &
                       state,state_lo,state_hi, &
                       delta,xlo,xhi)
  use eos_module
  use eos_type_module
  use network, only: nspec
  use probdata_module
  use meth_params_module, only : URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, UFS
  implicit none

  integer :: level, nvar
  integer :: lo(3), hi(3)
  integer :: state_lo(3), state_hi(3)
  double precision :: xlo(3), xhi(3), time, delta(3)
  double precision :: state(state_lo(1):state_hi(1), &
                            state_lo(2):state_hi(2), &
                            state_lo(3):state_hi(3),nvar)

  double precision :: r,xcen,ycen,zcen,rhoe,T,Pt
  integer i,j,k
  double precision massfrac(nspec)
  double precision, parameter ::  Pi = 4.d0*atan(1.d0)

  type (eos_t) :: eos_state

  massfrac(:) = 0.0d0
  massfrac(1) = 1.0d0

  call build(eos_state)


  do k = lo(3), hi(3)
     zcen = xlo(3) + delta(3)*(float(k-lo(3)) + 0.5d0)
     
     do j = lo(2), hi(2)
        ycen = xlo(2) + delta(2)*(float(j-lo(2)) + 0.5d0)

        do i = lo(1), hi(1)
           xcen = xlo(1) + delta(1)*(float(i-lo(1)) + 0.5d0)

           r = sqrt(xcen*xcen + ycen*ycen + zcen*zcen)

           if (r .gt. rpulse) then
              state(i,j,k,urho) = rho0
           else
              state(i,j,k,urho) = rho0 + drho0*exp(-16.d0*r*r)*(cos(Pi*r))**6
           end if

           eos_state%rho = state(i,j,k,urho)
           eos_state%p = Pt
           eos_state%T = 100000.d0  ! initial guess
           eos_state%massfrac(:) = massfrac(:)

           call eos_rp(eos_state)
 
           rhoe = state(i,j,k,urho)*eos_state%e
           T = eos_state%T

           state(i,j,k,ueint) = rhoe

           state(i,j,k,umx:umz) = 0.d0
           state(i,j,k,ueden) = rhoe
           state(i,j,k,utemp) = T

           
        enddo
     enddo
  enddo

end subroutine pc_initdata

subroutine pc_prob_close() &
     bind(C, name="pc_prob_close")

end subroutine pc_prob_close
