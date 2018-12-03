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

    !namelist /fortin/ pamb, phi_in, T_in, vn_in, pertmag, pmf_datafile

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

    !     Read namelists
    untin = 9
    open(untin,file=probin(1:namlen),form='formatted',status='old')
    !read(untin,fortin)
    close(unit=untin)

    call init_bc()

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
       dx,xlo,xhi) bind(C, name = "pc_initdata")
    use eos_type_module
    use probdata_module
    use meth_params_module, only : URHO, UMX, UMY, UMZ, UEINT, UEDEN, UTEMP, UFS
    use prob_params_module, only : problo, probhi, dim
    use eos_module
    use chemistry_module, only : nspecies, get_species_index
    use bl_constants_module, only: M_PI

    implicit none

    integer :: level, nvar
    integer :: lo(3), hi(3)
    integer :: state_lo(3), state_hi(3)
    double precision :: xlo(3), xhi(3), time, dx(3)
    double precision :: state(state_lo(1):state_hi(1), &
         state_lo(2):state_hi(2), &
         state_lo(3):state_hi(3),nvar)


    integer :: iN2, iO2, iH2, nPMF, i, j, k, nimages(3), ii, jj, kk
    double precision :: phi, a, cellshift, r, Tt, Pt, Xt(nspecies), rfire, u(3), x, y, z
    double precision T0, T1, pamb

    type(eos_t) :: eos_state

    call build(eos_state)

    iN2 = get_species_index("N2")
    iO2 = get_species_index("O2")
    iH2 = get_species_index("H2")

    phi = 1.d0
    T0 = 300.d0
    T1 = 1400.d0
    a = 0.5d0 ! for H2-air
    pamb = 1013250.d0 * 10
    nimages = 0
    nimages(1:dim) = 3
    rfire = .2d0 * (probhi(1) - problo(1))
    cellshift = 0.5d0

    call build(eos_state)

    do k = lo(3), hi(3)
       z = problo(3) + dx(3)*(dble(k) + HALF)
       do j = lo(2), hi(2)
          y = problo(2) + dx(2)*(dble(j) + HALF)
          do i = lo(1), hi(1)
             x = problo(1) + dx(1)*(dble(i) + HALF)

             Xt = 0.d0
             Xt(iO2) = 1.d0/(1.d0 + phi/a  + 0.79d0/0.21d0)
             Xt(iH2) = phi * Xt(iO2) / a

             Tt = T0
             Pt = pamb

             do kk = -nimages(3), nimages(3)
                do jj = -nimages(2), nimages(2)
                   do ii = -nimages(1), nimages(1)

                      if (dim .gt. 1) then
                         y = problo(2) + dx(2)*(j+cellshift) + jj * (probhi(2) - problo(2))
                      else
                         y = 0.d0
                      end if

                      if (dim .gt. 2) then
                         z = problo(3) + dx(3)*(k+cellshift) + kk * (probhi(3) - problo(3))
                      else
                         z = 0.d0
                      end if

                      x =    problo(1) + dx(1)*(i+cellshift) + ii * (probhi(1) - problo(1))

                      r = sqrt(x**2+y**2+z**2)

                      Pt = Pt        + 0.1d0*pamb * exp(-(r / rfire)**2)
                      Tt = Tt           + (T1-T0) * exp(-(r / rfire)**2)
                      Xt(iH2) = Xt(iH2) + 0.025d0 * exp(-(r / rfire)**2)
                      Xt(iO2) = Xt(iO2) - 0.050d0 * exp(-(r / rfire)**2)

                   end do

                end do

             end do

             Xt(iN2) = 1.d0 - Xt(iH2) - Xt(iO2)

             eos_state % molefrac = Xt
             eos_state % T        = Tt
             eos_state % p        = Pt
             u = 0.d0

             call eos_xty(eos_state)
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
