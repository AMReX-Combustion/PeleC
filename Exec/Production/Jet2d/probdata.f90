module probdata_module

  implicit none

  ! parameters for jet
  integer, save :: prob_type
  character (len=128), save :: turbfile = ""  
  double precision, save :: turb_boost_factor
  double precision, save :: pamb, phi_in, T_in, vn_in, T_co, vn_co
  double precision, save :: splitx, xfrontw, Tfrontw
  double precision, save :: blobr, blobx, bloby, blobT
  double precision, save :: inflow_period, inflow_vnmag
  double precision, save :: splity, yfrontw
  double precision, allocatable, save :: fuel_Y(:), air_Y(:)
  double precision, allocatable, save :: fuel_state(:), air_state(:)

  logical, save :: jet_initialized = .false.

  double precision, save :: center(3)

  ! These determine the refinement criteria
  double precision, save :: denerr,   dengrad
  double precision, save :: velerr,   velgrad
  double precision, save :: presserr, pressgrad
  double precision, save :: temperr,  tempgrad
  double precision, save :: vorterr,  vortgrad
  double precision, save :: tracerr
  integer         , save :: max_denerr_lev    = -1
  integer         , save :: max_dengrad_lev   = -1
  integer         , save :: max_velerr_lev    = -1
  integer         , save :: max_velgrad_lev   = -1
  integer         , save :: max_presserr_lev  = -1 
  integer         , save :: max_pressgrad_lev = -1
  integer         , save :: max_temperr_lev   = -1
  integer         , save :: max_tempgrad_lev  = -1
  integer         , save :: max_vorterr_lev   = -1
  integer         , save :: max_vortgrad_lev  = -1
  integer         , save :: max_tracerr_lev   = -1

contains

  subroutine init_jet

    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UTEMP, UFS
    use chemistry_module, only : nspecies, get_species_index
    use eos_module

    integer :: iN2, iO2, iH2
    double precision :: vt, ek

    type(eos_t) :: eos_state

    call build(eos_state)

    allocate(fuel_state(NVAR))
    allocate( air_state(NVAR))
    allocate(fuel_Y(nspecies))
    allocate( air_Y(nspecies))

    iN2 = get_species_index("N2")
    iO2 = get_species_index("O2")
    iH2 = get_species_index("H2")

    ! ----- Fuel -----
    eos_state % molefrac      = 0.d0
    eos_state % molefrac(iH2) = phi_in
    eos_state % molefrac(iN2) = 1.d0 - eos_state % molefrac(iH2)

    eos_state % p = pamb
    eos_state % T = T_in

    call eos_xty(eos_state) ! get mass fractions from mole fractions
    call eos_tp(eos_state)

    vt = vn_in
    ek = 0.5d0*vt**2

    fuel_state(URHO ) = eos_state % rho
    fuel_state(UMX  ) = 0.d0
    fuel_state(UMY  ) = eos_state % rho * vt
    fuel_state(UMZ  ) = 0.d0
    fuel_state(UEDEN) = eos_state % rho * (eos_state % e + ek)
    fuel_state(UTEMP) = eos_state % T
    fuel_state(UFS:UFS+nspecies-1) = eos_state % rho * eos_state % massfrac(1:nspecies)
    fuel_Y = eos_state % massfrac

    ! ----- Air -----
    eos_state % molefrac      = 0.d0
    eos_state % molefrac(iN2) = 0.79d0
    eos_state % molefrac(iO2) = 1.d0 - eos_state % molefrac(iN2)

    eos_state % p = pamb
    eos_state % T = T_co

    call eos_xty(eos_state) ! get mass fractions from mole fractions
    call eos_tp(eos_state)

    vt = vn_co
    ek = 0.5d0*vt**2

    air_state(URHO ) = eos_state % rho
    air_state(UMX  ) = 0.d0
    air_state(UMY  ) = eos_state % rho * vt
    air_state(UMZ  ) = 0.d0
    air_state(UEDEN) = eos_state % rho * (eos_state % e + ek)
    air_state(UTEMP) = eos_state % T
    air_state(UFS:UFS+nspecies-1) = eos_state % rho * eos_state % massfrac(1:nspecies)
    air_Y = eos_state % massfrac

    jet_initialized = .true.

  end subroutine init_jet

  subroutine inflow_boundary(xx, eos_state, uu)

    use eos_module
    use eos_type_module
    use chemistry_module, only : nspecies

    implicit none

    double precision, intent(in)    :: xx(3)
    type (eos_t),     intent(inout) :: eos_state
    double precision, intent(inout) :: uu(3)

    double precision :: x, y, z, eta, eta1, r, u1, u2, u3, sigma
    integer :: n

    x = xx(1)
    y = xx(2)
    z = xx(3)
    sigma = 2.5d0*xfrontw*splitx

    if (prob_type .eq. 0) then

       eta = 0.5d0 * (tanh((x + splitx)/sigma)   &
            &       - tanh((x - splitx)/sigma))

       do n=1,nspecies
          eos_state % massfrac(n) = eta*fuel_Y(n) + (1.d0-eta)*air_Y(n)
       end do
       eos_state % T  = eta * T_in + (1.d0-eta) * T_co
       u1 = 0.d0
       u2 = eta  * vn_in + (1.d0-eta ) * vn_co
       u3 = 0.d0

    else if (prob_type .eq. 1) then

       eta = 0.5d0 * (tanh((x + splitx)/Tfrontw)  &
            &       - tanh((x - splitx)/Tfrontw))
       eta1 = 0.5d0 * (tanh((x + blobr)/xfrontw)  &
            &        - tanh((x - blobr)/xfrontw))

       do n=1,nspecies
          eos_state % massfrac(n) = eta*fuel_Y(n) + (1.d0-eta)*air_Y(n)
       end do
       eos_state % T  = eta * T_in + (1.d0-eta) * T_co
       u1 = 0.d0
       u2 = eta1 * vn_in + (1.d0-eta1) * vn_co
       u3 = 0.d0

       if (blobr .gt. 0.d0) then
          eta = 0.5d0*(1.d0 - TANH(-2.d0*(y-bloby)/Tfrontw))
          eos_state % T  = eta * T_co + (1.d0-eta) * eos_state % T
          do n=1,nspecies
             eos_state % massfrac(n) = eta*air_Y(n) + (1.d0-eta)*eos_state % massfrac(n)
          end do
                    
          ! Superimpose blob of hot air
          r = SQRT((x-blobx)**2 + (y-bloby)**2)
          eta = 0.5d0*(1.d0 - TANH(2.d0*(r-blobr)/Tfrontw))
          do n=1,nspecies
             eos_state % massfrac(n) = eta*air_Y(n) + (1.d0-eta)*eos_state % massfrac(n)
          enddo
          eos_state % T  = eta * blobT + (1.d0-eta) * eos_state % T
       end if

    else if (prob_type .eq. 2) then

       eta = 0.5d0 * (tanh((x + splitx)/xfrontw)  &
            &       - tanh((x - splitx)/xfrontw))
       if ((y-splity) < 5.d0*yfrontw) then
          eta = eta * 0.5d0*(1.d0-tanh((y -splity)/yfrontw))
       else
          eta = 0.d0
       end if

       do n=1,nspecies
          eos_state % massfrac(n) = eta*fuel_Y(n) + (1.d0-eta)*air_Y(n)
       end do
       eos_state % T  = eta * T_in + (1.d0-eta) * T_co
       u1 = 0.d0
       u2 = eta * vn_in + (1.d0-eta) * vn_co
       u3 = 0.d0

    end if

    eos_state % p = pamb

    call eos_tp(eos_state)

    uu(1) = u1
    uu(2) = u2
    uu(3) = u3

  end subroutine inflow_boundary

end module probdata_module
