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
  double precision :: problo(dim), probhi(dim)
  character(len=72) :: pmf_datafile = ''

  integer untin,i

  namelist /fortin/ pamb, pmf_datafile, prob_axis, flame_offset

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

  pamb = 1013250.d0*100.0  ! 1 atm
  
  flame_offset = 0.d0
  
  !     Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

  call initialize_pmf(pmf_datafile)
  call init_bc()

  L = 1.d0
  L(1:dim) = probhi(1:dim) - problo(1:dim)
  
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
  use prob_params_module, only : problo, dim
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
  integer :: i, j, k, n, nPMF
  double precision :: yl, yr, x, y, z, u(3), pert
  double precision :: pmf_vals(nspec+4)

  type (eos_t) :: eos_state

  if (.not. bc_initialized) then
     call bl_error('bc_init not called prior to pc_initdata')
  endif

  call build(eos_state)

  do k = lo(3), hi(3)
     z = problo(3) + dx(3)*(dble(k) + HALF)
     do j = lo(2), hi(2)
        y = problo(2) + dx(2)*(dble(j) + HALF)
        do i = lo(1), hi(1)
           x = problo(1) + dx(1)*(dble(i) + HALF)

              pert = 0.d0
              if (prob_axis .eq. 1) then

                 call pmf(x+flame_offset,x+flame_offset,pmf_vals,nPMF)

              else if (prob_axis .eq. 2) then

                 call pmf(y+flame_offset,y+flame_offset,pmf_vals,nPMF)

              else

                 call pmf(z+flame_offset,z+flame_offset,pmf_vals,nPMF)

              endif

              eos_state % molefrac(1:nspec) = pmf_vals(4:3+nspec)
              eos_state % T = pmf_vals(1)

              eos_state % p = pmf_vals(nspec+4) !pamb
              u = 0
              u(prob_axis) = pmf_vals(2)


              call eos_xty(eos_state) ! get mass fractions from mole fractions
              call eos_tp(eos_state)

              state(i,j,k,URHO ) = eos_state % rho
              state(i,j,k,UMX  ) = eos_state % rho  *  u(1)
              state(i,j,k,UMY  ) = eos_state % rho  *  u(2)
              state(i,j,k,UMZ  ) = eos_state % rho  *  u(3)
              state(i,j,k,UEINT) = eos_state % rho  *  eos_state % e
              state(i,j,k,UEDEN) = eos_state % rho  * (eos_state % e + 0.5d0 * (u(1)**2 + u(2)**2 + u(3)**2))
              state(i,j,k,UTEMP) = eos_state % T
              state(i,j,k,UFS:UFS+nspec-1) = eos_state % rho  *  eos_state % massfrac(1:nspec)

        end do
     end do
  end do

  call destroy(eos_state)

end subroutine pc_initdata

subroutine pc_prob_close() &
     bind(C, name="pc_prob_close")

  use probdata_module

  call clear_bc()

end subroutine pc_prob_close

end module pc_prob_module
