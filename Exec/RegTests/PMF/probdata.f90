module probdata_module

  use eos_module
  use pmf_module

  implicit none

  double precision, save :: pamb, phi_in, T_in, vn_in, L(3), pertmag
  double precision, save, allocatable :: fuel_state(:)
  logical, save :: bc_initialized = .false.

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

  subroutine init_bc

    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEINT, UEDEN, UTEMP, UFS
    use chemistry_module, only : nspecies, get_species_index
    use eos_type_module

    integer :: iN2, iO2, iH2, nPMF
    double precision :: vt, ek, a, yl, yr, sumY
    double precision, allocatable :: pmf_vals(:)

    type(eos_t) :: eos_state

    call build(eos_state)
    allocate(pmf_vals(nspec+3))

    iN2 = get_species_index("N2")
    iO2 = get_species_index("O2")
    iH2 = get_species_index("H2")

    ! ----- Fuel -----
    if (phi_in .lt. 0) then
       yl = 0.d0
       yr = 0.d0 ! FIXME: get plo into saved data somehow
       call pmf(yl,yr,pmf_vals,nPMF)
       eos_state % molefrac(1:nspec) = MAX(0.d0,pmf_vals(4:3+nspec))
       eos_state % molefrac(iN2) = 1.d0 - (sum(eos_state % molefrac(1:nspec)) - eos_state % molefrac(iN2))
       eos_state % T = pmf_vals(1)
       vn_in = pmf_vals(2)
    else
       a = 0.5d0 ! for H2-air
       eos_state % molefrac = 0.d0
       eos_state % molefrac(iO2) = 1.d0/(1.d0 + phi_in/a  + 0.79d0/0.21d0)
       eos_state % molefrac(iH2) = phi_in * eos_state % molefrac(iO2) / a
       eos_state % molefrac(iN2) = 1.d0 - eos_state % molefrac(iH2) - eos_state % molefrac(iO2)
       eos_state % T = T_in
    endif
    eos_state % p = pamb

    call eos_xty(eos_state) ! get mass fractions from mole fractions
    call eos_tp(eos_state)

    vt = vn_in
    ek = 0.5d0*vt**2

    if (allocated(fuel_state)) then
       call bl_error('fuel_state already allocated')
    endif
    allocate(fuel_state(NVAR))

    fuel_state(URHO ) = eos_state % rho
    fuel_state(UMX  ) = 0.d0
    fuel_state(UMY  ) = eos_state % rho * vt
    fuel_state(UMZ  ) = 0.d0
    fuel_state(UEINT) = eos_state % rho * eos_state % e
    fuel_state(UEDEN) = eos_state % rho * (eos_state % e + ek)
    fuel_state(UTEMP) = eos_state % T
    fuel_state(UFS:UFS+nspecies-1) = eos_state % rho * eos_state % massfrac(1:nspecies)

    bc_initialized = .true.

    call destroy(eos_state)
    deallocate(pmf_vals)

  end subroutine init_bc


  subroutine clear_bc()

    deallocate(fuel_state)

  end subroutine clear_bc

end module probdata_module
