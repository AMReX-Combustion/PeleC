
! This file is automatically created by parse_pelec_params.py.  To update
! or add runtime parameters, please edit _cpp_parameters and then run
! mk_params.sh

! This module stores the runtime parameters and integer names for 
! indexing arrays.
!
! The Fortran-specific parameters are initialized in set_method_params(),
! and the ones that we are mirroring from C++ and obtaining through the
! ParmParse module are initialized in set_pelec_method_params().

module meth_params_module

  use amrex_fort_module

  implicit none

  ! number of ghost cells for the hyperbolic solver
  integer, parameter     :: NHYP    = 4

  ! NTHERM: number of thermodynamic variables
  integer, save :: NTHERM, NVAR
  integer, save :: URHO, UMX, UMY, UMZ, UMR, UML, UMP, UEDEN, UEINT, UTEMP, UFA, UFS, UFX
  integer, save :: USHK

  ! QTHERM: number of primitive variables
  integer, save :: QTHERM, QVAR
  integer, parameter :: QRHO=1, QU=2, QV=3, QW=4, QPRES=6, QREINT=7, QTEMP=8, QGAME=5
  integer, save :: QFS=9
  integer, save :: NQAUX, QGAMC, QC, QCSML, QDPDR, QDPDE, QRSPEC
  integer, save :: QFA, QFX
  !integer, save :: QRHO, QU, QV, QW, QPRES, QREINT, QTEMP, QGAME
  !integer, save :: QFA, QFS, QFX

  integer, save :: nadv

  ! NQ will be the total number of primitive variables, hydro + radiation
  integer, save :: NQ         

  integer, save :: npassive
  integer, save, allocatable :: qpass_map(:), upass_map(:)

  ! These are used for the Godunov state
  ! Note that the velocity indices here are picked to be the same value
  ! as in the primitive variable array
  integer, save :: NGDNV, GDRHO, GDU, GDV, GDW, GDPRES, GDGAME


  ! This for keeping track of particles states, and 
  integer, save :: PLOC, PVEL, PTEMP, PDIA, PRHO, PSPC
  integer, save :: PFVEL, PFRHO, PFTEMP, PFP, PFSPC


  integer         , save :: numpts_1d

  double precision, save, allocatable :: outflow_data_old(:,:)
  double precision, save, allocatable :: outflow_data_new(:,:)
  double precision, save :: outflow_data_old_time
  double precision, save :: outflow_data_new_time
  logical         , save :: outflow_data_allocated
  double precision, save :: max_dist

  double precision, save :: diffuse_cutoff_density

  ! these flags are for interpreting the EXT_DIR BCs
  integer, parameter :: EXT_UNDEFINED = -1
  integer, parameter :: EXT_HSE = 1
  
  integer, save :: xl_ext, yl_ext, zl_ext, xr_ext, yr_ext, zr_ext

  ! Create versions of these variables on the GPU
  ! the device update is then done in PeleC_nd.f90

  !$acc declare &
  !$acc create(NTHERM, NVAR) &
  !$acc create(URHO, UMX, UMY, UMZ, UMR, UML, UMP, UEDEN, UEINT, UTEMP, UFA, UFS,UFX) &
  !$acc create(USHK) &
  !$acc create(QTHERM, QVAR) &
  !$acc create(QRHO, QU, QV, QW, QPRES, QREINT, QTEMP) &
  !$acc create(QGAMC, QGAME) &
  !$acc create(NQ) &
  !$acc create(QFA, QFS, QFX) &
  !$acc create(xl_ext, yl_ext, zl_ext, xr_ext, yr_ext, zr_ext)

  ! Begin the declarations of the ParmParse parameters

  integer         , save :: levmsk_interior
  integer         , save :: levmsk_covered
  integer         , save :: levmsk_notcovered
  integer         , save :: levmsk_physbnd
  double precision, save :: difmag
  double precision, save :: small_dens
  double precision, save :: small_massfrac
  double precision, save :: small_temp
  double precision, save :: small_pres
  double precision, save :: small_ener
  integer         , save :: do_hydro
  integer         , save :: do_mol_AD
  integer         , save :: nscbc_adv
  integer         , save :: nscbc_diff
  integer         , save :: hybrid_hydro
  integer         , save :: ppm_type
  integer         , save :: weno_variant
  integer         , save :: ppm_trace_sources
  integer         , save :: ppm_temp_fix
  integer         , save :: ppm_predict_gammae
  integer         , save :: ppm_reference_eigenvectors
  integer         , save :: plm_iorder
  integer         , save :: hybrid_riemann
  integer         , save :: riemann_solver
  integer         , save :: cg_maxiter
  double precision, save :: cg_tol
  integer         , save :: cg_blend
  integer         , save :: use_flattening
  integer         , save :: transverse_use_eos
  integer         , save :: transverse_reset_density
  integer         , save :: transverse_reset_rhoe
  integer         , save :: dual_energy_update_E_from_e
  double precision, save :: dual_energy_eta1
  double precision, save :: dual_energy_eta2
  double precision, save :: dual_energy_eta3
  integer         , save :: use_pslope
  integer         , save :: fix_mass_flux
  integer         , save :: limit_fluxes_on_small_dens
  integer         , save :: density_reset_method
  integer         , save :: allow_negative_energy
  integer         , save :: allow_small_energy
  integer         , save :: first_order_hydro
  character (len=128), save :: xl_ext_bc_type
  character (len=128), save :: xr_ext_bc_type
  character (len=128), save :: yl_ext_bc_type
  character (len=128), save :: yr_ext_bc_type
  character (len=128), save :: zl_ext_bc_type
  character (len=128), save :: zr_ext_bc_type
  double precision, save :: eb_small_vfrac
  integer         , save :: do_mms
  double precision, save :: cfl
  double precision, save :: dtnuc_e
  double precision, save :: dtnuc_X
  integer         , save :: dtnuc_mode
  double precision, save :: dxnuc
  integer         , save :: do_react
  double precision, save :: react_T_min
  double precision, save :: react_T_max
  double precision, save :: react_rho_min
  double precision, save :: react_rho_max
  integer         , save :: disable_shock_burning
  integer         , save :: do_acc
  integer         , save :: track_grid_losses

  !$acc declare &
  !$acc create(levmsk_interior, levmsk_covered, levmsk_notcovered) &
  !$acc create(levmsk_physbnd, difmag, small_dens) &
  !$acc create(small_massfrac, small_temp, small_pres) &
  !$acc create(small_ener, do_hydro, do_mol_AD) &
  !$acc create(nscbc_adv, nscbc_diff, hybrid_hydro) &
  !$acc create(ppm_type, weno_variant, ppm_trace_sources) &
  !$acc create(ppm_temp_fix, ppm_predict_gammae, ppm_reference_eigenvectors) &
  !$acc create(plm_iorder, hybrid_riemann, riemann_solver) &
  !$acc create(cg_maxiter, cg_tol, cg_blend) &
  !$acc create(use_flattening, transverse_use_eos, transverse_reset_density) &
  !$acc create(transverse_reset_rhoe, dual_energy_update_E_from_e, dual_energy_eta1) &
  !$acc create(dual_energy_eta2, dual_energy_eta3, use_pslope) &
  !$acc create(fix_mass_flux, limit_fluxes_on_small_dens, density_reset_method) &
  !$acc create(allow_negative_energy, allow_small_energy, first_order_hydro) &
  !$acc create(eb_small_vfrac, do_mms, cfl) &
  !$acc create(dtnuc_e, dtnuc_X, dtnuc_mode) &
  !$acc create(dxnuc, do_react, react_T_min) &
  !$acc create(react_T_max, react_rho_min, react_rho_max) &
  !$acc create(disable_shock_burning, do_acc, track_grid_losses)

  ! End the declarations of the ParmParse parameters

  double precision, save :: rot_vec(3)

contains

  subroutine set_pelec_method_params() bind(C,name="set_pelec_method_params")

    use parmparse_module, only: parmparse_build, parmparse_destroy, ParmParse

    implicit none

    type (ParmParse) :: pp

    call parmparse_build(pp, "pelec")

    levmsk_interior = 0;
    levmsk_covered = 1;
    levmsk_notcovered = 2;
    levmsk_physbnd = 3;
    difmag = 0.1d0;
    small_dens = 1.d-200;
    small_massfrac = 1.d-200;
    small_temp = 1.d-200;
    small_pres = 1.d-200;
    small_ener = -1.d200;
    do_hydro = -1;
    do_mol_AD = 0;
    nscbc_adv = 1;
    nscbc_diff = 0;
    hybrid_hydro = 0;
    ppm_type = 1;
    weno_variant = 1;
    ppm_trace_sources = 1;
    ppm_temp_fix = 0;
    ppm_predict_gammae = 0;
    ppm_reference_eigenvectors = 0;
    plm_iorder = 2;
    hybrid_riemann = 0;
    riemann_solver = 0;
    cg_maxiter = 12;
    cg_tol = 1.0d-5;
    cg_blend = 2;
    use_flattening = 1;
    transverse_use_eos = 0;
    transverse_reset_density = 1;
    transverse_reset_rhoe = 0;
    dual_energy_update_E_from_e = 1;
    dual_energy_eta1 = 1.0d0;
    dual_energy_eta2 = 1.0d-4;
    dual_energy_eta3 = 1.0d0;
    use_pslope = 1;
    fix_mass_flux = 0;
    limit_fluxes_on_small_dens = 0;
    density_reset_method = 1;
    allow_negative_energy = 1;
    allow_small_energy = 1;
    first_order_hydro = 0;
    xl_ext_bc_type = "";
    xr_ext_bc_type = "";
    yl_ext_bc_type = "";
    yr_ext_bc_type = "";
    zl_ext_bc_type = "";
    zr_ext_bc_type = "";
    eb_small_vfrac = 1.0d-2;
    do_mms = 0;
    cfl = 0.8d0;
    dtnuc_e = 1.d200;
    dtnuc_X = 1.d200;
    dtnuc_mode = 1;
    dxnuc = 1.d200;
    do_react = 0;
    react_T_min = 0.0d0;
    react_T_max = 1.d200;
    react_rho_min = 0.0d0;
    react_rho_max = 1.d200;
    disable_shock_burning = 0;
    do_acc = -1;
    track_grid_losses = 0;

    call pp%query("levmsk_interior", levmsk_interior)
    call pp%query("levmsk_covered", levmsk_covered)
    call pp%query("levmsk_notcovered", levmsk_notcovered)
    call pp%query("levmsk_physbnd", levmsk_physbnd)
    call pp%query("difmag", difmag)
    call pp%query("small_dens", small_dens)
    call pp%query("small_massfrac", small_massfrac)
    call pp%query("small_temp", small_temp)
    call pp%query("small_pres", small_pres)
    call pp%query("small_ener", small_ener)
    call pp%query("do_hydro", do_hydro)
    call pp%query("do_mol_AD", do_mol_AD)
    call pp%query("nscbc_adv", nscbc_adv)
    call pp%query("nscbc_diff", nscbc_diff)
    call pp%query("hybrid_hydro", hybrid_hydro)
    call pp%query("ppm_type", ppm_type)
    call pp%query("weno_variant", weno_variant)
    call pp%query("ppm_trace_sources", ppm_trace_sources)
    call pp%query("ppm_temp_fix", ppm_temp_fix)
    call pp%query("ppm_predict_gammae", ppm_predict_gammae)
    call pp%query("ppm_reference_eigenvectors", ppm_reference_eigenvectors)
    call pp%query("plm_iorder", plm_iorder)
    call pp%query("hybrid_riemann", hybrid_riemann)
    call pp%query("riemann_solver", riemann_solver)
    call pp%query("cg_maxiter", cg_maxiter)
    call pp%query("cg_tol", cg_tol)
    call pp%query("cg_blend", cg_blend)
    call pp%query("use_flattening", use_flattening)
    call pp%query("transverse_use_eos", transverse_use_eos)
    call pp%query("transverse_reset_density", transverse_reset_density)
    call pp%query("transverse_reset_rhoe", transverse_reset_rhoe)
    call pp%query("dual_energy_update_E_from_e", dual_energy_update_E_from_e)
    call pp%query("dual_energy_eta1", dual_energy_eta1)
    call pp%query("dual_energy_eta2", dual_energy_eta2)
    call pp%query("dual_energy_eta3", dual_energy_eta3)
    call pp%query("use_pslope", use_pslope)
    call pp%query("fix_mass_flux", fix_mass_flux)
    call pp%query("limit_fluxes_on_small_dens", limit_fluxes_on_small_dens)
    call pp%query("density_reset_method", density_reset_method)
    call pp%query("allow_negative_energy", allow_negative_energy)
    call pp%query("allow_small_energy", allow_small_energy)
    call pp%query("first_order_hydro", first_order_hydro)
    call pp%query("xl_ext_bc_type", xl_ext_bc_type)
    call pp%query("xr_ext_bc_type", xr_ext_bc_type)
    call pp%query("yl_ext_bc_type", yl_ext_bc_type)
    call pp%query("yr_ext_bc_type", yr_ext_bc_type)
    call pp%query("zl_ext_bc_type", zl_ext_bc_type)
    call pp%query("zr_ext_bc_type", zr_ext_bc_type)
    call pp%query("eb_small_vfrac", eb_small_vfrac)
    call pp%query("do_mms", do_mms)
    call pp%query("cfl", cfl)
    call pp%query("dtnuc_e", dtnuc_e)
    call pp%query("dtnuc_X", dtnuc_X)
    call pp%query("dtnuc_mode", dtnuc_mode)
    call pp%query("dxnuc", dxnuc)
    call pp%query("do_react", do_react)
    call pp%query("react_T_min", react_T_min)
    call pp%query("react_T_max", react_T_max)
    call pp%query("react_rho_min", react_rho_min)
    call pp%query("react_rho_max", react_rho_max)
    call pp%query("disable_shock_burning", disable_shock_burning)
    call pp%query("do_acc", do_acc)
    call pp%query("track_grid_losses", track_grid_losses)

    !$acc update &
    !$acc device(levmsk_interior, levmsk_covered, levmsk_notcovered) &
    !$acc device(levmsk_physbnd, difmag, small_dens) &
    !$acc device(small_massfrac, small_temp, small_pres) &
    !$acc device(small_ener, do_hydro, do_mol_AD) &
    !$acc device(nscbc_adv, nscbc_diff, hybrid_hydro) &
    !$acc device(ppm_type, weno_variant, ppm_trace_sources) &
    !$acc device(ppm_temp_fix, ppm_predict_gammae, ppm_reference_eigenvectors) &
    !$acc device(plm_iorder, hybrid_riemann, riemann_solver) &
    !$acc device(cg_maxiter, cg_tol, cg_blend) &
    !$acc device(use_flattening, transverse_use_eos, transverse_reset_density) &
    !$acc device(transverse_reset_rhoe, dual_energy_update_E_from_e, dual_energy_eta1) &
    !$acc device(dual_energy_eta2, dual_energy_eta3, use_pslope) &
    !$acc device(fix_mass_flux, limit_fluxes_on_small_dens, density_reset_method) &
    !$acc device(allow_negative_energy, allow_small_energy, first_order_hydro) &
    !$acc device(eb_small_vfrac, do_mms, cfl) &
    !$acc device(dtnuc_e, dtnuc_X, dtnuc_mode) &
    !$acc device(dxnuc, do_react, react_T_min) &
    !$acc device(react_T_max, react_rho_min, react_rho_max) &
    !$acc device(disable_shock_burning, do_acc, track_grid_losses)


    ! now set the external BC flags
    select case (xl_ext_bc_type)
    case ("hse", "HSE")
       xl_ext = EXT_HSE
    case default
       xl_ext = EXT_UNDEFINED
    end select

    select case (yl_ext_bc_type)
    case ("hse", "HSE")
       yl_ext = EXT_HSE
    case default
       yl_ext = EXT_UNDEFINED
    end select

    select case (zl_ext_bc_type)
    case ("hse", "HSE")
       zl_ext = EXT_HSE
    case default
       zl_ext = EXT_UNDEFINED
    end select

    select case (xr_ext_bc_type)
    case ("hse", "HSE")
       xr_ext = EXT_HSE
    case default
       xr_ext = EXT_UNDEFINED
    end select

    select case (yr_ext_bc_type)
    case ("hse", "HSE")
       yr_ext = EXT_HSE
    case default
       yr_ext = EXT_UNDEFINED
    end select

    select case (zr_ext_bc_type)
    case ("hse", "HSE")
       zr_ext = EXT_HSE
    case default
       zr_ext = EXT_UNDEFINED
    end select

    !$acc update device(xl_ext, yl_ext, zl_ext, xr_ext, yr_ext, zr_ext)

    call parmparse_destroy(pp)

  end subroutine set_pelec_method_params

end module meth_params_module
