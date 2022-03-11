
// This file is automatically created by parse_pelec_params.py.  To update
// or add runtime parameters, please edit _cpp_parameters and then run
// mk_params.sh

int PeleC::levmsk_interior = 0;
int PeleC::levmsk_covered = 1;
int PeleC::levmsk_notcovered = 2;
int PeleC::levmsk_physbnd = 3;
int PeleC::state_interp_order = 1;
int PeleC::lin_limit_state_interp = 1;
int PeleC::state_nghost = 0;
int PeleC::do_reflux = 1;
int PeleC::do_avg_down = 1;
amrex::Real PeleC::difmag = 0.1;
amrex::Real PeleC::small_pres = 1.e-200;
int PeleC::do_hydro = -1;
int PeleC::do_mol = 0;
int PeleC::nscbc_adv = 1;
int PeleC::nscbc_diff = 0;
int PeleC::add_ext_src = 0;
int PeleC::add_forcing_src = 0;
int PeleC::ppm_type = 0;
int PeleC::ppm_trace_sources = 0;
int PeleC::plm_iorder = 2;
int PeleC::use_laxf_flux = 0;
int PeleC::riemann_solver = 0;
int PeleC::use_flattening = 1;
int PeleC::dual_energy_update_E_from_e = 1;
amrex::Real PeleC::dual_energy_eta2 = 1.0e-4;
int PeleC::allow_negative_energy = 1;
int PeleC::allow_small_energy = 1;
std::string PeleC::xl_ext_bc_type = "";
std::string PeleC::xr_ext_bc_type = "";
std::string PeleC::yl_ext_bc_type = "";
std::string PeleC::yr_ext_bc_type = "";
std::string PeleC::zl_ext_bc_type = "";
std::string PeleC::zr_ext_bc_type = "";
int PeleC::do_les = 0;
int PeleC::use_explicit_filter = 0;
amrex::Real PeleC::Cs = 0.0;
amrex::Real PeleC::CI = 0.0;
amrex::Real PeleC::PrT = 1.0;
amrex::Real PeleC::eb_boundary_T = 1.0;
int PeleC::eb_isothermal = 1;
int PeleC::eb_noslip = 1;
std::string PeleC::redistribution_type = "StateRedist";
int PeleC::eb_clean_massfrac = 1;
amrex::Real PeleC::eb_clean_massfrac_threshold = 0.0;
int PeleC::eb_srd_max_order = 2;
int PeleC::do_mms = 0;
std::string PeleC::masa_solution_name = "ad_cns_3d_les";
amrex::Real PeleC::fixed_dt = -1.0;
amrex::Real PeleC::initial_dt = -1.0;
amrex::Real PeleC::dt_cutoff = 0.0;
amrex::Real PeleC::max_dt = 1.e200;
amrex::Real PeleC::cfl = 0.8;
amrex::Real PeleC::init_shrink = 1.0;
amrex::Real PeleC::change_max = 1.1;
int PeleC::sdc_iters = 1;
int PeleC::mol_iters = 1;
int PeleC::do_react = 0;
std::string PeleC::chem_integrator = "ReactorNull";
int PeleC::bndry_func_thread_safe = 1;
#ifdef AMREX_DEBUG
int PeleC::print_energy_diagnostics = 1;
#else
int PeleC::print_energy_diagnostics = 0;
#endif
int PeleC::track_grid_losses = 0;
int PeleC::sum_interval = -1;
int PeleC::track_extrema = 1;
std::string PeleC::extrema_spec_name = "";
amrex::Real PeleC::sum_per = -1.0e0;
int PeleC::hard_cfl_limit = 1;
std::string PeleC::job_name = "";
std::string PeleC::flame_trac_name = "";
std::string PeleC::fuel_name = "";
