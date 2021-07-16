#include <memory>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_Vector.H>
#include <AMReX_TagBox.H>

#ifdef PELEC_USE_EB
#include <AMReX_EBMultiFabUtil.H>
#include "hydro_redistribution.H"
#endif

#ifdef AMREX_PARTICLES
#include <AMReX_Particles.H>
#endif

#ifdef PELEC_USE_MASA
#include <masa.h>
using namespace MASA;
#endif

#include "PeleC.H"
#include "Derive.H"
#include "prob.H"
#include "PelePhysics.H"
#include "Timestep.H"
#include "Utilities.H"
#include "Tagging.H"
#include "IndexDefines.H"
#if defined(PELEC_USE_REACTIONS) && defined(USE_SUNDIALS_PP)
#include "reactor.H"
#endif

#ifdef PELEC_ENABLE_FPE_TRAP
#if defined(__linux__)
#include <cfenv>
#elif defined(__APPLE__)
#include <fenv.h>
#endif
#endif

bool PeleC::signalStopJob = false;
bool PeleC::dump_old = false;
int PeleC::verbose = 0;
int PeleC::radius_grow = 1;
amrex::BCRec PeleC::phys_bc;
amrex::Real PeleC::frac_change = std::numeric_limits<amrex::Real>::max();
int PeleC::Density = -1;
int PeleC::Eden = -1;
int PeleC::Eint = -1;
int PeleC::Temp = -1;
int PeleC::Xmom = -1;
int PeleC::Ymom = -1;
int PeleC::Zmom = -1;
int PeleC::FirstSpec = -1;
int PeleC::FirstAux = -1;
int PeleC::NumAdv = 0;
int PeleC::FirstAdv = -1;
int PeleC::pstateVel = -1;
int PeleC::pstateT = -1;
int PeleC::pstateDia = -1;
int PeleC::pstateRho = -1;
int PeleC::pstateY = -1;
int PeleC::pstateNum = 0;

#include "pelec_defaults.H"

int PeleC::diffuse_temp = 0;
int PeleC::diffuse_enth = 0;
int PeleC::diffuse_spec = 0;
int PeleC::diffuse_vel = 0;
amrex::Real PeleC::diffuse_cutoff_density =
  -std::numeric_limits<amrex::Real>::max();
bool PeleC::do_diffuse = false;

#ifdef PELEC_USE_MASA
bool PeleC::mms_initialized = false;
#endif

int PeleC::use_hybrid_weno = 0;
int PeleC::weno_scheme = 1;

int PeleC::les_model = 0;
int PeleC::les_filter_type = no_filter;
int PeleC::les_filter_fgr = 1;
int PeleC::les_test_filter_type = box_3pt_optimized_approx;
int PeleC::les_test_filter_fgr = 2;

bool PeleC::eb_in_domain = false;
#ifdef PELEC_USE_EB
bool PeleC::eb_initialized = false;
bool PeleC::body_state_set = false;
amrex::GpuArray<amrex::Real, NVAR> PeleC::body_state;
#endif

bool PeleC::do_react_load_balance = false;
bool PeleC::do_mol_load_balance = false;

amrex::Vector<std::string> PeleC::spec_names;

amrex::Vector<int> PeleC::src_list;

// this will be reset upon restart
amrex::Real PeleC::previousCPUTimeUsed = 0.0;
amrex::Real PeleC::startCPUTime = 0.0;
int PeleC::num_state_type = 0;

#ifdef PELEC_USE_EB
static bool eb_initialized = false;

bool
ebInitialized()
{
  return eb_initialized;
}

void
ebInitialized(bool eb_init_val)
{
  eb_initialized = eb_init_val;
}
#endif

void
PeleC::read_params()
{
  static bool read_params_done = false;

  if (read_params_done) {
    return;
  }

  read_params_done = true;

  amrex::ParmParse pp("pelec");

#include "pelec_queries.H"

  pp.query("v", verbose);
  pp.query("sum_interval", sum_interval);
  pp.query("dump_old", dump_old);

  // Get boundary conditions
  amrex::Vector<std::string> lo_bc_char(AMREX_SPACEDIM);
  amrex::Vector<std::string> hi_bc_char(AMREX_SPACEDIM);
  pp.getarr("lo_bc", lo_bc_char, 0, AMREX_SPACEDIM);
  pp.getarr("hi_bc", hi_bc_char, 0, AMREX_SPACEDIM);

  amrex::Vector<int> lo_bc(AMREX_SPACEDIM);
  amrex::Vector<int> hi_bc(AMREX_SPACEDIM);
  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    if (lo_bc_char[dir] == "Interior") {
      lo_bc[dir] = 0;
    } else if (lo_bc_char[dir] == "Hard") {
      lo_bc[dir] = 1;
    } else if (lo_bc_char[dir] == "FOExtrap") {
      lo_bc[dir] = 2;
    } else if (lo_bc_char[dir] == "Symmetry") {
      lo_bc[dir] = 3;
    } else if (lo_bc_char[dir] == "SlipWall") {
      lo_bc[dir] = 4;
    } else if (lo_bc_char[dir] == "NoSlipWall") {
      lo_bc[dir] = 5;
    } else if (lo_bc_char[dir] == "UserBC") {
      lo_bc[dir] = 6;
    } else {
      amrex::Abort("Wrong boundary condition word in lo_bc, please use: "
                   "Interior, UserBC, Symmetry, SlipWall, NoSlipWall");
    }

    if (hi_bc_char[dir] == "Interior") {
      hi_bc[dir] = 0;
    } else if (hi_bc_char[dir] == "Hard") {
      hi_bc[dir] = 1;
    } else if (hi_bc_char[dir] == "FOExtrap") {
      hi_bc[dir] = 2;
    } else if (hi_bc_char[dir] == "Symmetry") {
      hi_bc[dir] = 3;
    } else if (hi_bc_char[dir] == "SlipWall") {
      hi_bc[dir] = 4;
    } else if (hi_bc_char[dir] == "NoSlipWall") {
      hi_bc[dir] = 5;
    } else if (hi_bc_char[dir] == "UserBC") {
      hi_bc[dir] = 6;
    } else {
      amrex::Abort("Wrong boundary condition word in hi_bc, please use: "
                   "Interior, UserBC, Symmetry, SlipWall, NoSlipWall");
    }
  }

  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    phys_bc.setLo(dir, lo_bc[dir]);
    phys_bc.setHi(dir, hi_bc[dir]);
  }

  // Check phys_bc against possible periodic geometry
  // if periodic, must have internal BC marked.
  // Check, periodic means interior in those directions.
  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    if (amrex::DefaultGeometry().isPeriodic(dir)) {
      if (lo_bc[dir] != Interior && amrex::ParallelDescriptor::IOProcessor()) {
        std::cerr << "PeleC::read_params:periodic in direction " << dir
                  << " but low BC is not Interior\n";
        amrex::Error();
      }
      if (hi_bc[dir] != Interior && amrex::ParallelDescriptor::IOProcessor()) {
        std::cerr << "PeleC::read_params:periodic in direction " << dir
                  << " but high BC is not Interior\n";
        amrex::Error();
      }
    } else {
      // If not periodic, should not be interior.
      if (lo_bc[dir] == Interior && amrex::ParallelDescriptor::IOProcessor()) {
        std::cerr << "PeleC::read_params:interior bc in direction " << dir
                  << " but not periodic\n";
        amrex::Error();
      }
      if (hi_bc[dir] == Interior && amrex::ParallelDescriptor::IOProcessor()) {
        std::cerr << "PeleC::read_params:interior bc in direction " << dir
                  << " but not periodic\n";
        amrex::Error();
      }
    }
  }

  if (amrex::DefaultGeometry().IsRZ() && (lo_bc[0] != Symmetry)) {
    amrex::Error("PeleC::read_params: must set r=0 boundary condition to "
                 "Symmetry for r-z");
  }

  // TODO: Any reason to support spherical in PeleC?
  if (amrex::DefaultGeometry().IsRZ()) {
    amrex::Abort("We don't support cylindrical coordinate systems in 3D");
  } else if (amrex::DefaultGeometry().IsSPHERICAL()) {
    amrex::Abort("We don't support spherical coordinate systems in 3D");
  }

  pp.query("diffuse_temp", diffuse_temp);
  pp.query("diffuse_enth", diffuse_enth);
  pp.query("diffuse_spec", diffuse_spec);
  pp.query("diffuse_vel", diffuse_vel);
  pp.query("diffuse_cutoff_density", diffuse_cutoff_density);

  do_diffuse = diffuse_temp || diffuse_enth || diffuse_spec || diffuse_vel;

  // sanity checks
  if (cfl <= 0.0 || cfl > 1.0) {
    amrex::Error("Invalid CFL factor; must be between zero and one.");
  }
  if ((do_hydro == 1) && (do_mol == 1) && (cfl > 0.3)) {
    amrex::Print() << "WARNING -- CFL should be <= 0.3 when using MOL hydro."
                   << std::endl;
  }

  if ((do_les || use_explicit_filter) && (AMREX_SPACEDIM != 3)) {
    amrex::Abort("Using LES/filtering currently requires 3d.");
  }

  if (do_les) {
    pp.query("les_model", les_model);
    pp.query("les_test_filter_type", les_test_filter_type);
    pp.query("les_test_filter_fgr", les_test_filter_fgr);
  }

  if (use_explicit_filter) {
    pp.query("les_filter_type", les_filter_type);
    pp.query("les_filter_fgr", les_filter_fgr);
  }

  // Check on PPM type
  if ((do_hydro == 1) && (do_mol == 0)) {
    if (ppm_type != 0 && ppm_type != 1) {
      amrex::Error("PeleC::ppm_type must be 0 (PLM) or 1 (PPM)");
    }
  }

  pp.query("use_hybrid_weno", use_hybrid_weno);
  pp.query("weno_scheme", weno_scheme);
  if (ppm_type != 1 && use_hybrid_weno == 1) {
    amrex::Error("PeleC::ppm_type must be 1 (PPM) to use WENO method");
  }

  // for the moment, ppm_type = 0 does not support ppm_trace_sources --
  // we need to add the momentum sources to the states (and not
  // add it in trans_3d
  if (ppm_type == 0 && ppm_trace_sources == 1) {
    amrex::Print()
      << "WARNING: ppm_trace_sources = 1 not implemented for ppm_type = 0"
      << std::endl;
    ppm_trace_sources = 0;
    pp.add("ppm_trace_sources", ppm_trace_sources);
  }
  /*
    if (ppm_temp_fix > 0 && AMREX_SPACEDIM == 1) {
      amrex::Error("ppm_temp_fix > 0 not implemented in 1-d");
    }

    if (hybrid_riemann == 1 && AMREX_SPACEDIM == 1) {
      amrex::Error("hybrid_riemann only implemented in 2- and 3-d");
    }
  */
  if (
    hybrid_riemann == 1 && (amrex::DefaultGeometry().IsSPHERICAL() ||
                            amrex::DefaultGeometry().IsRZ())) {
    amrex::Error(
      "hybrid_riemann should only be used for Cartesian coordinates");
  }

  if (use_colglaz >= 0) {
    amrex::Error("use_colglaz is deprecated. Use riemann_solver instead");
  }

  if (max_dt < fixed_dt) {
    amrex::Error("Cannot have max_dt < fixed_dt");
  }

#ifdef AMREX_PARTICLES
  readParticleParams();
#endif

#ifdef PELEC_USE_EB
  if ((do_mol == 0) && (eb_in_domain)) {
    amrex::Abort("Must do_mol = 1 when using EB\n");
  }
#endif

  // Read tagging parameters
  read_tagging_params();

  // TODO: What is this?
  amrex::StateDescriptor::setBndryFuncThreadSafety(bndry_func_thread_safe);

  // Get some useful amr inputs
  amrex::ParmParse ppa("amr");

  // This turns on the lb stuff inside Amr, but we use our own flag to signal
  // whether to gather data
  ppa.query("loadbalance_with_workestimates", do_mol_load_balance);
  ppa.query("loadbalance_with_workestimates", do_react_load_balance);
}

PeleC::PeleC()
  : old_sources(num_src),
    new_sources(num_src)
#ifdef PELEC_USE_MASA
    ,
    mms_src_evaluated(false)
#endif
{
  nGrowF = 0;
  // Is this relevant for PeleC?
  for (int i = 0; i < n_lost; i++) {
    material_lost_through_boundary_cumulative[i] = 0.0;
    material_lost_through_boundary_temp[i] = 0.0;
  }
}

PeleC::PeleC(
  amrex::Amr& papa,
  int lev,
  const amrex::Geometry& level_geom,
  const amrex::BoxArray& bl,
  const amrex::DistributionMapping& dm,
  amrex::Real time)
  : AmrLevel(papa, lev, level_geom, bl, dm, time),
    old_sources(num_src),
    new_sources(num_src)
#ifdef PELEC_USE_MASA
    ,
    mms_src_evaluated(false)
#endif
{
  buildMetrics();

#ifdef PELEC_USE_EB
  init_eb(level_geom, bl, dm);
#endif

  amrex::MultiFab& S_new = get_new_data(State_Type);

  for (int n = 0; n < src_list.size(); ++n) {
    int oldGrow = numGrow();
    int newGrow = S_new.nGrow();
#ifdef AMREX_PARTICLES
    if (src_list[n] == spray_src) {
      oldGrow = 1;
      newGrow = amrex::max<amrex::Real>(1, newGrow);
    }
#endif
    old_sources[src_list[n]] = std::make_unique<amrex::MultiFab>(
      grids, dmap, NVAR, oldGrow, amrex::MFInfo(), Factory());
    new_sources[src_list[n]] = std::make_unique<amrex::MultiFab>(
      grids, dmap, NVAR, newGrow, amrex::MFInfo(), Factory());
  }

  if (do_hydro) {
    Sborder.define(grids, dmap, NVAR, numGrow(), amrex::MFInfo(), Factory());
  } else if (do_diffuse) {
    Sborder.define(grids, dmap, NVAR, numGrow(), amrex::MFInfo(), Factory());
  }
#ifdef AMREX_PARTICLES
  else if (do_spray_particles) {
    Sborder.define(grids, dmap, NVAR, numGrow(), amrex::MFInfo(), Factory());
  }
#endif

  if (!do_mol) {
    if (do_hydro) {
      hydro_source.define(
        grids, dmap, NVAR, numGrow(), amrex::MFInfo(), Factory());

      // This array holds the sum of all source terms that affect the
      // hydrodynamics. If we are doing the source term predictor, we'll also
      // use this after the hydro update to store the sum of the new-time
      // sources, so that we can compute the time derivative of the source
      // terms.
      sources_for_hydro.define(
        grids, dmap, NVAR, numGrow(), amrex::MFInfo(), Factory());
    }
  } else {
    Sborder.define(grids, dmap, NVAR, numGrow(), amrex::MFInfo(), Factory());
  }

  // Is this relevant for PeleC?
  for (int i = 0; i < n_lost; i++) {
    material_lost_through_boundary_cumulative[i] = 0.0;
    material_lost_through_boundary_temp[i] = 0.0;
  }

  if (do_reflux && level > 0) {
    flux_reg.define(
      bl, papa.boxArray(level - 1), dm, papa.DistributionMap(level - 1),
      level_geom, papa.Geom(level - 1), papa.refRatio(level - 1), level, NVAR);

    if (!amrex::DefaultGeometry().IsCartesian()) {
      pres_reg.define(
        bl, papa.boxArray(level - 1), dm, papa.DistributionMap(level - 1),
        level_geom, papa.Geom(level - 1), papa.refRatio(level - 1), level, 1);
    }
  }

#ifdef PELEC_USE_REACTIONS
  get_new_data(Reactions_Type).setVal(0.0);
#endif

  // Don't need this in pure C++?
  // initialize the Godunov state array used in hydro -- we wait
  // until here so that ngroups is defined (if needed) in
  // rad_params_module
  // if (do_hydro)
  //{
  //  init_godunov_indices();
  //}

  // initialize LES variables
  if (do_les) {
    init_les();
  }

  // initialize filters and variables
  nGrowF = 0;
  if (use_explicit_filter) {
    init_filters();
  }
}

PeleC::~PeleC() = default;

void
PeleC::buildMetrics()
{
  const int ngrd = grids.size();

  radius.resize(ngrd);

  const amrex::Real* dx = geom.CellSize();

  for (int i = 0; i < ngrd; i++) {
    const amrex::Box& b = grids[i];
    int ilo = b.smallEnd(0) - radius_grow;
    int ihi = b.bigEnd(0) + radius_grow;
    int len = ihi - ilo + 1;

    radius[i].resize(len);

    amrex::Real* rad = radius[i].dataPtr();

    if (amrex::DefaultGeometry().IsCartesian()) {
      for (int j = 0; j < len; j++) {
        rad[j] = 1.0;
      }
    } else {
      amrex::RealBox gridloc =
        amrex::RealBox(grids[i], geom.CellSize(), geom.ProbLo());

      const amrex::Real xlo = gridloc.lo(0) + (0.5 - radius_grow) * dx[0];

      for (int j = 0; j < len; j++) {
        rad[j] = xlo + j * dx[0];
      }
    }
  }

  volume.clear();
  volume.define(
    grids, dmap, 1, numGrow(), amrex::MFInfo(), amrex::FArrayBoxFactory());
  geom.GetVolume(volume);

  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    area[dir].clear();
    area[dir].define(
      getEdgeBoxArray(dir), dmap, 1, numGrow(), amrex::MFInfo(),
      amrex::FArrayBoxFactory());
    geom.GetFaceArea(area[dir], dir);
  }

#ifdef PELEC_USE_EB
  vfrac.clear();
  vfrac.define(grids, dmap, 1, numGrow(), amrex::MFInfo(), Factory());
  const auto& ebfactory =
    dynamic_cast<amrex::EBFArrayBoxFactory const&>(Factory());
  amrex::MultiFab::Copy(vfrac, ebfactory.getVolFrac(), 0, 0, 1, numGrow());
  areafrac = ebfactory.getAreaFrac();
#endif

  level_mask.clear();
  level_mask.define(grids, dmap, 1, 1);
  level_mask.BuildMask(
    geom.Domain(), geom.periodicity(), levmsk_covered, levmsk_notcovered,
    levmsk_physbnd, levmsk_interior);

  if (level == 0) {
    setGridInfo();
  }
}

void
PeleC::setTimeLevel(amrex::Real time, amrex::Real dt_old, amrex::Real dt_new)
{
  AmrLevel::setTimeLevel(time, dt_old, dt_new);
}

void
PeleC::setGridInfo()
{
  // Send refinement data to Fortran. We do it here
  // because now the grids have been initialized and
  // we need this data for setting up the problem.
  // Note that this routine will always get called
  // on level 0, even if we are doing a restart,
  // so it is safe to put this here.
  if (level == 0) {
    const int max_level = parent->maxLevel();
    const int nlevs = max_level + 1;
    const int size = 3 * nlevs;

    amrex::Vector<amrex::Real> dx_level(size);
    amrex::Vector<int> domlo_level(size);
    amrex::Vector<int> domhi_level(size);

    const amrex::Real* dx_coarse = geom.CellSize();

    const int* domlo_coarse = geom.Domain().loVect();
    const int* domhi_coarse = geom.Domain().hiVect();

    for (int dir = 0; dir < 3; dir++) {
      dx_level[dir] = (ZFILL(dx_coarse))[dir];

      domlo_level[dir] = (ARLIM_3D(domlo_coarse))[dir];
      domhi_level[dir] = (ARLIM_3D(domhi_coarse))[dir];
    }

    for (int lev = 1; lev <= max_level; lev++) {
      amrex::IntVect ref_ratio = parent->refRatio(lev - 1);

      // Note that we are explicitly calculating here what the
      // data would be on refined levels rather than getting the
      // data directly from those levels, because some potential
      // refined levels may not exist at the beginning of the simulation.

      for (int dir = 0; dir < 3; dir++) {
        if (dir < AMREX_SPACEDIM) {
          dx_level[3 * lev + dir] =
            dx_level[3 * (lev - 1) + dir] / ref_ratio[dir];
          int ncell = (domhi_level[3 * (lev - 1) + dir] -
                       domlo_level[3 * (lev - 1) + dir] + 1) *
                      ref_ratio[dir];
          domlo_level[3 * lev + dir] = domlo_level[dir];
          domhi_level[3 * lev + dir] = domlo_level[3 * lev + dir] + ncell - 1;
        } else {
          dx_level[3 * lev + dir] = 0.0;
          domlo_level[3 * lev + dir] = 0;
          domhi_level[3 * lev + dir] = 0;
        }
      }
    }

    // Old fortran call
    // set_grid_info(max_level, dx_level, domlo_level, domhi_level);
  }
}

/*
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_prob_close()
{
}
*/

void
PeleC::initData()
{
  BL_PROFILE("PeleC::initData()");

  // Copy problem parameter structs to device
  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, PeleC::h_prob_parm_device,
    PeleC::h_prob_parm_device + 1, PeleC::d_prob_parm_device);

  // int ns = NVAR;
  amrex::MultiFab& S_new = get_new_data(State_Type);
  // amrex::Real cur_time = state[State_Type].curTime();

  S_new.setVal(0.0);

#if AMREX_SPACEDIM > 1
  // make sure dx = dy = dz -- that's all we guarantee to support
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();
  const amrex::Real small = 1.e-13;
  if (
    amrex::max<amrex::Real>(AMREX_D_DECL(
      static_cast<amrex::Real>(0.0),
      static_cast<amrex::Real>(amrex::Math::abs(dx[0] - dx[1])),
      static_cast<amrex::Real>(amrex::Math::abs(dx[0] - dx[2])))) >
    small * dx[0]) {
    amrex::Abort("dx != dy != dz not supported");
  }
#endif

  if (verbose) {
    amrex::Print() << "Initializing the data at level " << level << std::endl;
  }

#ifdef PELEC_USE_REACTIONS
  get_new_data(Reactions_Type).setVal(0.0);
#endif

  if (do_mol_load_balance || do_react_load_balance) {
    get_new_data(Work_Estimate_Type).setVal(1.0);
  }

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  for (amrex::MFIter mfi(S_new, amrex::TilingIfNotGPU()); mfi.isValid();
       ++mfi) {
    const amrex::Box& box = mfi.tilebox();
    auto sfab = S_new.array(mfi);
    const auto geomdata = geom.data();

    const ProbParmDevice* lprobparm = d_prob_parm_device;

    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      pc_initdata(i, j, k, sfab, geomdata, *lprobparm);
      // Verify that the sum of (rho Y)_i = rho at every cell
      pc_check_initial_species(i, j, k, sfab);
    });
  }

  enforce_consistent_e(S_new);

  // computeTemp(S_new,0);

#ifdef PELEC_USE_EB
  set_body_state(S_new);
  InitialRedistribution();
#endif

#ifdef AMREX_PARTICLES
  if (level == 0) {
    initParticles();
  } else {
    // TODO: Determine how many ghost cells to use here
    int nGrow = 0;
    particleRedistribute(level - 1, nGrow, 0, false);
  }
#endif

  if (verbose) {
    amrex::Print() << "Done initializing level " << level << " data "
                   << std::endl;
  }
}

void
PeleC::init(AmrLevel& old)
{
  BL_PROFILE("PeleC::init(old)");

  auto* oldlev = (PeleC*)&old;

  // Create new grid data by fillpatching from old.
  amrex::Real dt_new = parent->dtLevel(level);
  amrex::Real cur_time = oldlev->state[State_Type].curTime();
  amrex::Real prev_time = oldlev->state[State_Type].prevTime();
  amrex::Real dt_old = cur_time - prev_time;
  setTimeLevel(cur_time, dt_old, dt_new);

  amrex::MultiFab& S_new = get_new_data(State_Type);
  FillPatch(old, S_new, 0, cur_time, State_Type, 0, NVAR);

#ifdef PELEC_USE_REACTIONS
  amrex::MultiFab& React_new = get_new_data(Reactions_Type);

  if (do_react) {
    FillPatch(
      old, React_new, 0, cur_time, Reactions_Type, 0, React_new.nComp());
  } else {
    React_new.setVal(0);
  }
#endif

  if (do_mol_load_balance || do_react_load_balance) {
    amrex::MultiFab& work_estimate_new = get_new_data(Work_Estimate_Type);
    FillPatch(
      old, work_estimate_new, 0, cur_time, Work_Estimate_Type, 0,
      work_estimate_new.nComp());
  }
}

void
PeleC::init()
{
  // This version inits the data on a new level that did not
  // exist before regridding.
  BL_PROFILE("PeleC::init()");

  amrex::Real dt = parent->dtLevel(level);
  amrex::Real cur_time = getLevel(level - 1).state[State_Type].curTime();
  amrex::Real prev_time = getLevel(level - 1).state[State_Type].prevTime();

  amrex::Real dt_old =
    (cur_time - prev_time) / (amrex::Real)parent->MaxRefRatio(level - 1);

  setTimeLevel(cur_time, dt_old, dt);
  amrex::MultiFab& S_new = get_new_data(State_Type);
  FillCoarsePatch(S_new, 0, cur_time, State_Type, 0, NVAR);

  if (do_mol_load_balance || do_react_load_balance) {
    amrex::MultiFab& work_estimate_new = get_new_data(Work_Estimate_Type);
    int ncomp = work_estimate_new.nComp();
    FillCoarsePatch(
      work_estimate_new, 0, cur_time, Work_Estimate_Type, 0, ncomp);
  }
}

amrex::Real
PeleC::initialTimeStep()
{
  BL_PROFILE("PeleC::initialTimeStep()");

  amrex::Real init_dt = 0.0;

  if (initial_dt > 0.0) {
    init_dt = initial_dt;
  } else {
    const amrex::Real dummy_dt = 0.0;
    init_dt = init_shrink * estTimeStep(dummy_dt);
  }

  return init_dt;
}

amrex::Real PeleC::estTimeStep(amrex::Real /*dt_old*/)
{
  BL_PROFILE("PeleC::estTimeStep()");

  if (fixed_dt > 0.0) {
    return fixed_dt;
  }

  // set_amr_info(level, -1, -1, -1.0, -1.0);

  amrex::Real estdt = max_dt;

  const amrex::MultiFab& stateMF = get_new_data(State_Type);

  const amrex::Real* dx = geom.CellSize();

  std::string limiter = "pelec.max_dt";

  // Start the hydro with the max_dt value, but divide by CFL
  // to account for the fact that we multiply by it at the end.
  // This ensures that if max_dt is more restrictive than the hydro
  // criterion, we will get exactly max_dt for a timestep.

  const amrex::Real max_dt_over_cfl = max_dt / cfl;
  amrex::Real estdt_hydro = max_dt_over_cfl;
  amrex::Real estdt_vdif = max_dt_over_cfl;
  amrex::Real estdt_tdif = max_dt_over_cfl;
  amrex::Real estdt_edif = max_dt_over_cfl;
  if (do_hydro || do_mol || diffuse_vel || diffuse_temp || diffuse_enth) {

#ifdef PELEC_USE_EB
    auto const& fact =
      dynamic_cast<amrex::EBFArrayBoxFactory const&>(stateMF.Factory());
    auto const& flags = fact.getMultiEBCellFlagFab();
#endif

    prefetchToDevice(stateMF); // This should accelerate the below operations.
    amrex::Real AMREX_D_DECL(dx1 = dx[0], dx2 = dx[1], dx3 = dx[2]);

    if (do_hydro) {
      amrex::Real dt = amrex::ReduceMin(
        stateMF,
#ifdef PELEC_USE_EB
        flags,
#endif
        0,
        [=] AMREX_GPU_HOST_DEVICE(
          amrex::Box const& bx, const amrex::Array4<const amrex::Real>& fab_arr
#ifdef PELEC_USE_EB
          ,
          const amrex::Array4<const amrex::EBCellFlag>& flag_arr
#endif
          ) noexcept -> amrex::Real {
          return pc_estdt_hydro(
            bx, fab_arr,
#ifdef PELEC_USE_EB
            flag_arr,
#endif
            AMREX_D_DECL(dx1, dx2, dx3));
        });
      estdt_hydro = amrex::min<amrex::Real>(estdt_hydro, dt);
    }

    if (diffuse_vel) {
      pele::physics::transport::TransParm const* ltransparm =
        pele::physics::transport::trans_parm_g;
      amrex::Real dt = amrex::ReduceMin(
        stateMF,
#ifdef PELEC_USE_EB
        flags,
#endif
        0,
        [=] AMREX_GPU_HOST_DEVICE(
          amrex::Box const& bx, const amrex::Array4<const amrex::Real>& fab_arr
#ifdef PELEC_USE_EB
          ,
          const amrex::Array4<const amrex::EBCellFlag>& flag_arr
#endif
          ) noexcept -> amrex::Real {
          return pc_estdt_veldif(
            bx, fab_arr,
#ifdef PELEC_USE_EB
            flag_arr,
#endif
            AMREX_D_DECL(dx1, dx2, dx3), ltransparm);
        });
      estdt_vdif = amrex::min<amrex::Real>(estdt_vdif, dt);
    }

    if (diffuse_temp) {
      pele::physics::transport::TransParm const* ltransparm =
        pele::physics::transport::trans_parm_g;
      amrex::Real dt = amrex::ReduceMin(
        stateMF,
#ifdef PELEC_USE_EB
        flags,
#endif
        0,
        [=] AMREX_GPU_HOST_DEVICE(
          amrex::Box const& bx, const amrex::Array4<const amrex::Real>& fab_arr
#ifdef PELEC_USE_EB
          ,
          const amrex::Array4<const amrex::EBCellFlag>& flag_arr
#endif
          ) noexcept -> amrex::Real {
          return pc_estdt_tempdif(
            bx, fab_arr,
#ifdef PELEC_USE_EB
            flag_arr,
#endif
            AMREX_D_DECL(dx1, dx2, dx3), ltransparm);
        });
      estdt_tdif = amrex::min<amrex::Real>(estdt_tdif, dt);
    }

    if (diffuse_enth) {
      pele::physics::transport::TransParm const* ltransparm =
        pele::physics::transport::trans_parm_g;
      amrex::Real dt = amrex::ReduceMin(
        stateMF,
#ifdef PELEC_USE_EB
        flags,
#endif
        0,
        [=] AMREX_GPU_HOST_DEVICE(
          amrex::Box const& bx, const amrex::Array4<const amrex::Real>& fab_arr
#ifdef PELEC_USE_EB
          ,
          const amrex::Array4<const amrex::EBCellFlag>& flag_arr
#endif
          ) noexcept -> amrex::Real {
          return pc_estdt_enthdif(
            bx, fab_arr,
#ifdef PELEC_USE_EB
            flag_arr,
#endif
            AMREX_D_DECL(dx1, dx2, dx3), ltransparm);
        });
      estdt_edif = amrex::min<amrex::Real>(estdt_edif, dt);
    }

    estdt_hydro = amrex::min<amrex::Real>(
      estdt_hydro,
      amrex::min<amrex::Real>(
        estdt_vdif, amrex::min<amrex::Real>(estdt_tdif, estdt_edif)));

    amrex::ParallelDescriptor::ReduceRealMin(estdt_hydro);
    estdt_hydro *= cfl;

    if (verbose) {
      amrex::Print() << "...estimated hydro-limited timestep at level " << level
                     << ": " << estdt_hydro << std::endl;
    }

    // Determine if this is more restrictive than the maximum timestep limiting
    if (estdt_hydro < estdt) {
      limiter = "hydro";
      estdt = estdt_hydro;
    }
  }

#ifdef AMREX_PARTICLES
  amrex::Real estdt_particle = max_dt;
  if (do_spray_particles) {
    particleEstTimeStep(estdt_particle);
    if (estdt_particle < estdt) {
      limiter = "particles";
      estdt = estdt_particle;
    }
  }
#endif

  if (verbose) {
    amrex::Print() << "PeleC::estTimeStep (" << limiter << "-limited) at level "
                   << level << ":  estdt = " << estdt << '\n';
  }

  return estdt;
}

void
PeleC::computeNewDt(
  int finest_level,
  int /*sub_cycle*/,
  amrex::Vector<int>& n_cycle,
  const amrex::Vector<amrex::IntVect>& /*ref_ratio*/,
  amrex::Vector<amrex::Real>& dt_min,
  amrex::Vector<amrex::Real>& dt_level,
  amrex::Real stop_time,
  int post_regrid_flag)
{
  BL_PROFILE("PeleC::computeNewDt()");

  // We are at the start of a coarse grid timecycle.
  // Compute the timesteps for the next iteration.
  if (level > 0) {
    return;
  }

  amrex::Real dt_0 = 1.0e+100;
  int n_factor = 1;
  for (int i = 0; i <= finest_level; i++) {
    PeleC& adv_level = getLevel(i);
    dt_min[i] = adv_level.estTimeStep(dt_level[i]);
  }

  if (fixed_dt <= 0.0) {
    if (post_regrid_flag == 1) {
      // Limit dt's by pre-regrid dt
      for (int i = 0; i <= finest_level; i++) {
        dt_min[i] = amrex::min<amrex::Real>(dt_min[i], dt_level[i]);
      }
    } else {
      // Limit dt's by change_max * old dt
      for (int i = 0; i <= finest_level; i++) {
        if (verbose && amrex::ParallelDescriptor::IOProcessor()) {
          if (dt_min[i] > change_max * dt_level[i]) {
            amrex::Print() << "PeleC::compute_new_dt : limiting dt at level "
                           << i << '\n';
            amrex::Print() << " ... new dt computed: " << dt_min[i] << '\n';
            amrex::Print() << " ... but limiting to: "
                           << change_max * dt_level[i] << " = " << change_max
                           << " * " << dt_level[i] << '\n';
          }
        }
        dt_min[i] =
          amrex::min<amrex::Real>(dt_min[i], change_max * dt_level[i]);
      }
    }
  }

  // Find the minimum over all levels
  for (int i = 0; i <= finest_level; i++) {
    n_factor *= n_cycle[i];
    dt_0 = amrex::min<amrex::Real>(dt_0, n_factor * dt_min[i]);
  }

  // Limit dt's by the value of stop_time.
  const amrex::Real dt_eps = 0.001 * dt_0;
  amrex::Real cur_time = state[State_Type].curTime();
  if (stop_time >= 0.0) {
    if ((cur_time + dt_0) > (stop_time - dt_eps)) {
      dt_0 = stop_time - cur_time;
    }
  }

  n_factor = 1;
  for (int i = 0; i <= finest_level; i++) {
    n_factor *= n_cycle[i];
    dt_level[i] = dt_0 / n_factor;
  }
}

void
PeleC::computeInitialDt(
  int finest_level,
  int /*sub_cycle*/,
  amrex::Vector<int>& n_cycle,
  const amrex::Vector<amrex::IntVect>& /*ref_ratio*/,
  amrex::Vector<amrex::Real>& dt_level,
  amrex::Real stop_time)
{
  BL_PROFILE("PeleC::computeInitialDt()");

  // Grids have been constructed, compute dt for all levels.
  if (level > 0) {
    return;
  }

  amrex::Real dt_0 = 1.0e+100;
  int n_factor = 1;
  // TODO: This will need to change for optimal subcycling.
  for (int i = 0; i <= finest_level; i++) {
    dt_level[i] = getLevel(i).initialTimeStep();
    n_factor *= n_cycle[i];
    dt_0 = amrex::min<amrex::Real>(dt_0, n_factor * dt_level[i]);
  }

  // Limit dt's by the value of stop_time.
  const amrex::Real dt_eps = 0.001 * dt_0;
  amrex::Real cur_time = state[State_Type].curTime();
  if (stop_time >= 0.0) {
    if ((cur_time + dt_0) > (stop_time - dt_eps)) {
      dt_0 = stop_time - cur_time;
    }
  }

  n_factor = 1;
  for (int i = 0; i <= finest_level; i++) {
    n_factor *= n_cycle[i];
    dt_level[i] = dt_0 / n_factor;
  }
}

void
PeleC::post_timestep(int
#ifdef AMREX_PARTICLES
                       iteration
#endif
)
{
  BL_PROFILE("PeleC::post_timestep()");

  const int finest_level = parent->finestLevel();

#ifdef AMREX_PARTICLES
  const int ncycle = parent->nCycle(level);
  if (do_spray_particles) {
    // Remove virtual particles at this level if we have any.
    if (theVirtPC() != 0)
      removeVirtualParticles();

    // Remove Ghost particles on the final iteration
    if (iteration == ncycle)
      removeGhostParticles();

    // Sync up if we're level 0 or if we have particles that may have moved
    // off the next finest level and need to be added to our own level.
    if ((iteration < ncycle && level < finest_level) || level == 0) {
      // TODO: Determine how many ghost cells to use here
      int nGrow = iteration;
      theSprayPC()->Redistribute(level, theSprayPC()->finestLevel(), nGrow);
    }
  }
#endif

  if (do_reflux && level < finest_level) {
    reflux();

    // We need to do this before anything else because refluxing changes the
    // values of coarse cells underneath fine grids with the assumption they'll
    // be over-written by averaging down
    if (level < finest_level) {
      avgDown();
    }

    // Clean up any aberrant state data generated by the reflux.
    // amrex::MultiFab& S_new_crse = get_new_data(State_Type);
    // clean_state(S_new_crse);
  }

  // Re-compute temperature after all the other updates.
  amrex::MultiFab& S_new = get_new_data(State_Type);
  int ng_pts = 0;
  computeTemp(S_new, ng_pts);

  problem_post_timestep();

  if (level == 0) {
    int nstep = parent->levelSteps(0);
    amrex::Real dtlev = parent->dtLevel(0);
    amrex::Real cumtime = parent->cumTime() + dtlev;

    bool sum_int_test = (sum_interval > 0 && nstep % sum_interval == 0);

    bool sum_per_test = false;

    if (sum_per > 0.0) {
      const int num_per_old =
        static_cast<int>(amrex::Math::floor((cumtime - dtlev) / sum_per));
      const int num_per_new =
        static_cast<int>(amrex::Math::floor((cumtime) / sum_per));

      if (num_per_old != num_per_new) {
        sum_per_test = true;
      }
    }

    if (sum_int_test || sum_per_test) {
      sum_integrated_quantities();
    }
  }
}

void
PeleC::post_restart()
{
  BL_PROFILE("PeleC::post_restart()");

  // Copy problem parameter structs to device
  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, PeleC::h_prob_parm_device,
    PeleC::h_prob_parm_device + 1, PeleC::d_prob_parm_device);

  // amrex::Real cur_time = state[State_Type].curTime();

#ifdef AMREX_PARTICLES
  if (do_spray_particles) {
    particlePostRestart(parent->theRestartFile());
  }
#endif

  // Don't need this in pure C++?
  // initialize the Godunov state array used in hydro -- we wait
  // until here so that ngroups is defined (if needed) in
  // rad_params_module
  // if (do_hydro) {
  //  init_godunov_indices();
  //}

  // initialize LES variables
  if (do_les) {
    init_les();
  }

  // initialize filters and variables
  nGrowF = 0;
  if (use_explicit_filter) {
    init_filters();
  }

  problem_post_restart();
}

void
PeleC::postCoarseTimeStep(amrex::Real cumtime)
{
  BL_PROFILE("PeleC::postCoarseTimeStep()");
  AmrLevel::postCoarseTimeStep(cumtime);
}

void
PeleC::post_regrid(
  int
#ifdef AMREX_PARTICLES
    lbase
#endif
  ,
  int /*new_finest*/)
{
  BL_PROFILE("PeleC::post_regrid()");
  fine_mask.clear();

#ifdef AMREX_PARTICLES
  if (do_spray_particles && theSprayPC() != 0 && level == lbase) {
    // TODO: Determine how many ghost cells to use here
    int nGrow = 0;
    particleRedistribute(lbase);
  }
#endif
}

void PeleC::post_init(amrex::Real /*stop_time*/)
{
  BL_PROFILE("PeleC::post_init()");

  amrex::Real dtlev = parent->dtLevel(level);
  amrex::Real cumtime = parent->cumTime();

  // Fill Reactions_Type data based on initial dt
#ifdef PELEC_USE_REACTIONS
  if (do_react == 1) {
    bool react_init = true;
    react_state(cumtime, dtlev, react_init);
  }
#endif

  if (level > 0) {
    return;
  }

  // Average data down from finer levels
  // so that conserved data is consistent between levels.
  if (do_avg_down != 0) {
    int finest_level = parent->finestLevel();
    for (int k = finest_level - 1; k >= 0; k--) {
      getLevel(k).avgDown();
    }
  }

  // Allow the user to define their own post_init functions.
  problem_post_init();

  int nstep = parent->levelSteps(0);
  if (cumtime != 0.0) {
    cumtime += dtlev;
  }

  bool sum_int_test = false;

  if (sum_interval > 0) {
    if (nstep % sum_interval == 0) {
      sum_int_test = true;
    }
  }

  bool sum_per_test = false;

  if (sum_per > 0.0) {
    const int num_per_old =
      static_cast<int>(amrex::Math::floor((cumtime - dtlev) / sum_per));
    const int num_per_new =
      static_cast<int>(amrex::Math::floor((cumtime) / sum_per));

    if (num_per_old != num_per_new) {
      sum_per_test = true;
    }
  }

  if (sum_int_test || sum_per_test) {
    sum_integrated_quantities();
  }
}

int
PeleC::okToContinue()
{
  if (level > 0) {
    return 1;
  }

  int test = 1;

  if (signalStopJob) {
    test = 0;

    amrex::Print()
      << " Signalling a stop of the run due to signalStopJob = true."
      << std::endl;
  } else if (parent->dtLevel(0) < dt_cutoff) {
    test = 0;

    amrex::Print() << " Signalling a stop of the run because dt < dt_cutoff."
                   << std::endl;
  }

  return test;
}

void
PeleC::reflux()
{
  BL_PROFILE("PeleC::reflux()");

  AMREX_ASSERT(level < parent->finestLevel());

  const amrex::Real strt = amrex::ParallelDescriptor::second();

  PeleC& fine_level = getLevel(level + 1);
  amrex::MultiFab& S_crse = get_new_data(State_Type);

#ifdef PELEC_USE_EB

  amrex::MultiFab& S_fine = fine_level.get_new_data(State_Type);

  fine_level.flux_reg.Reflux(S_crse, vfrac, S_fine, fine_level.vfrac);

  if (!amrex::DefaultGeometry().IsCartesian()) {
    amrex::Abort("rz not yet compatible with EB");
  }

#else

  fine_level.flux_reg.Reflux(S_crse);

  if (!amrex::DefaultGeometry().IsCartesian()) {
    amrex::MultiFab dr(
      volume.boxArray(), volume.DistributionMap(), 1, volume.nGrow(),
      amrex::MFInfo(), amrex::FArrayBoxFactory());
    dr.setVal(geom.CellSize(0));
    amrex::Abort("PeleC reflux not yet ready for r-z");
    // fine_level.pres_reg.Reflux(S_crse,dr,1.0,0,Xmom,1,geom);
  }
#endif

  if (verbose) {
    const int IOProc = amrex::ParallelDescriptor::IOProcessorNumber();
    amrex::Real end = amrex::ParallelDescriptor::second() - strt;

#ifdef AMREX_LAZY
    Lazy::QueueReduction([=]() mutable {
#endif
      amrex::ParallelDescriptor::ReduceRealMax(end, IOProc);

      amrex::Print() << "PeleC::reflux() at level " << level
                     << " : time = " << end << std::endl;
#ifdef AMREX_LAZY
    });
#endif
  }
}

void
PeleC::avgDown()
{
  BL_PROFILE("PeleC::avgDown()");

  if (level == parent->finestLevel()) {
    return;
  }

  avgDown(State_Type);

#ifdef PELEC_USE_REACTIONS
  avgDown(Reactions_Type);
#endif
}

void
PeleC::normalize_species(amrex::MultiFab& /*S*/)
{
  amrex::Abort("We don't normalize species!");
}

void
PeleC::enforce_consistent_e(amrex::MultiFab& S)
{
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  for (amrex::MFIter mfi(S, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    const amrex::Box& tbox = mfi.tilebox();
    const auto Sfab = S.array(mfi);
    amrex::ParallelFor(
      tbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const amrex::Real rhoInv = 1.0 / Sfab(i, j, k, URHO);
        const amrex::Real u = Sfab(i, j, k, UMX) * rhoInv;
        const amrex::Real v = Sfab(i, j, k, UMY) * rhoInv;
        const amrex::Real w = Sfab(i, j, k, UMZ) * rhoInv;
        Sfab(i, j, k, UEDEN) = Sfab(i, j, k, UEINT) + 0.5 *
                                                        Sfab(i, j, k, URHO) *
                                                        (u * u + v * v + w * w);
      });
  }
}

amrex::Real
PeleC::enforce_min_density(
  amrex::MultiFab& /*S_old*/, const amrex::MultiFab& S_new)
{
  // This routine sets the density in S_new to be larger than the density
  // floor. Note that it will operate everywhere on S_new, including ghost
  // zones. S_old is present so that, after the hydro call, we know what the old
  // density was so that we have a reference for comparison. If you are calling
  // it elsewhere and there's no meaningful reference state, just pass in the
  // same amrex::MultiFab twice.
  //  @return  The return value is the the negative fractional change in the
  // state that has the largest magnitude. If there is no reference state, this
  // is meaningless.

  amrex::Real dens_change = 1.0;
  amrex::Real mass_added = 0.0; // cppcheck-suppress variableScope
  amrex::Real eint_added = 0.0; // cppcheck-suppress variableScope
  amrex::Real eden_added = 0.0; // cppcheck-suppress variableScope

#ifdef PELEC_USE_EB
  auto const& fact =
    dynamic_cast<amrex::EBFArrayBoxFactory const&>(S_new.Factory());
  auto const& flags = fact.getMultiEBCellFlagFab();
#endif

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion()) \
  reduction(min                                           \
            : dens_change)
#endif
  for (amrex::MFIter mfi(S_new, amrex::TilingIfNotGPU()); mfi.isValid();
       ++mfi) {

#ifdef PELEC_USE_EB
    const amrex::Box& bx = mfi.growntilebox();
    const auto& flag_fab = flags[mfi];
    amrex::FabType typ = flag_fab.getType(bx);
    if (typ == amrex::FabType::covered) {
      continue;
    }
#endif

    // Not enabled on the GPU
    // const auto& stateold = S_old[mfi];
    // auto& statenew = S_new[mfi];
    // const auto& vol = volume[mfi];
    // enforce_minimum_density(stateold.dataPtr(), ARLIM_3D(stateold.loVect()),
    // ARLIM_3D(stateold.hiVect()),
    //                        statenew.dataPtr(), ARLIM_3D(statenew.loVect()),
    //                        ARLIM_3D(statenew.hiVect()), vol.dataPtr(),
    //                        ARLIM_3D(vol.loVect()), ARLIM_3D(vol.hiVect()),
    //                        ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
    //                        &mass_added, &eint_added, &eden_added,
    //                        &dens_change, &verbose);
  }

  if (print_energy_diagnostics) {
    amrex::Real foo[3] = {mass_added, eint_added, eden_added};

#ifdef AMREX_LAZY
    Lazy::QueueReduction([=]() mutable {
#endif
      amrex::ParallelDescriptor::ReduceRealSum(
        foo, 3, amrex::ParallelDescriptor::IOProcessorNumber());

      if (amrex::ParallelDescriptor::IOProcessor()) {
        mass_added = foo[0];
        eint_added = foo[1];
        eden_added = foo[2];

        if (amrex::Math::abs(mass_added) != 0.0) {
          amrex::Print() << "   Mass added from negative density correction : "
                         << mass_added << std::endl;
          amrex::Print() << "(rho e) added from negative density correction : "
                         << eint_added << std::endl;
          amrex::Print() << "(rho E) added from negative density correction : "
                         << eden_added << std::endl;
        }
      }
#ifdef AMREX_LAZY
    });
#endif
  }

  return dens_change;
}

void
PeleC::avgDown(int state_indx)
{
  BL_PROFILE("PeleC::avgDown(state_indx)");

  if (level == parent->finestLevel()) {
    return;
  }

  amrex::MultiFab& S_crse = get_new_data(state_indx);
  amrex::MultiFab& S_fine = getLevel(level + 1).get_new_data(state_indx);

#ifdef PELEC_USE_EB
  PeleC& fine_lev = getLevel(level + 1);

  amrex::EB_average_down(
    S_fine, S_crse, fine_lev.Volume(), fine_lev.volFrac(), 0, S_fine.nComp(),
    fine_ratio);

  if (state_indx == State_Type) {
    set_body_state(S_crse); // TODO: Is this necessary?
  }
#else

  const amrex::Geometry& fgeom = getLevel(level + 1).geom;
  const amrex::Geometry& cgeom = geom;

  amrex::average_down(
    S_fine, S_crse, fgeom, cgeom, 0, S_fine.nComp(), fine_ratio);

#endif
}

void
PeleC::allocOldData()
{
  for (int k = 0; k < num_state_type; k++) {
    state[k].allocOldData();
  }
}

void
PeleC::removeOldData()
{
  AmrLevel::removeOldData();
}

void
PeleC::errorEst(
  amrex::TagBoxArray& tags,
  int /*clearval*/,
  int /*tagval*/,
  amrex::Real time,
  int /*n_error_buf*/,
  int /*ngrow*/)
{
  BL_PROFILE("PeleC::errorEst()");

  amrex::MultiFab S_data(
    get_new_data(State_Type).boxArray(),
    get_new_data(State_Type).DistributionMap(), NVAR, 1);
  const amrex::Real cur_time = state[State_Type].curTime();
  FillPatch(
    *this, S_data, S_data.nGrow(), cur_time, State_Type, Density, NVAR, 0);

  amrex::Vector<amrex::BCRec> bcs(NVAR);
  const char tagval = amrex::TagBox::SET;
  // const char clearval = amrex::TagBox::CLEAR;

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  {
    for (amrex::MFIter mfi(S_data, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
      const amrex::Box& tilebox = mfi.tilebox();
      const auto Sfab = S_data.array(mfi);
      auto tag_arr = tags.array(mfi);
      const auto datbox = S_data[mfi].box();

#ifdef PELEC_USE_EB
      const auto vfrac_arr = vfrac.array(mfi);
#endif

      amrex::FArrayBox S_derData(datbox, 1, amrex::The_Async_Arena);
      auto S_derarr = S_derData.array();
      const int ncp = S_derData.nComp();
      const int* bc = bcs[0].data();

      // Tagging density
      if (level < tagging_parm->max_denerr_lev) {
        const amrex::Real captured_denerr = tagging_parm->denerr;
        amrex::ParallelFor(
          tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            tag_error(i, j, k, tag_arr, Sfab, captured_denerr, tagval);
          });
      }
      if (level < tagging_parm->max_dengrad_lev) {
        const amrex::Real captured_dengrad = tagging_parm->dengrad;
        amrex::ParallelFor(
          tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            tag_graderror(i, j, k, tag_arr, Sfab, captured_dengrad, tagval);
          });
      }

      // Tagging pressure
      S_derData.setVal<amrex::RunOn::Device>(0.0);
      pc_derpres(
        datbox, S_derData, ncp, Sfab.nComp(), S_data[mfi], geom, time, bc,
        level);
      if (level < tagging_parm->max_presserr_lev) {
        const amrex::Real captured_presserr = tagging_parm->presserr;
        amrex::ParallelFor(
          tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            tag_error(i, j, k, tag_arr, S_derarr, captured_presserr, tagval);
          });
      }
      if (level < tagging_parm->max_pressgrad_lev) {
        const amrex::Real captured_pressgrad = tagging_parm->pressgrad;
        amrex::ParallelFor(
          tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            tag_graderror(
              i, j, k, tag_arr, S_derarr, captured_pressgrad, tagval);
          });
      }

      // Tagging vel_x
      S_derData.setVal<amrex::RunOn::Device>(0.0);
      pc_dervelx(
        datbox, S_derData, ncp, Sfab.nComp(), S_data[mfi], geom, time, bc,
        level);
      if (level < tagging_parm->max_velerr_lev) {
        const amrex::Real captured_velerr = tagging_parm->velerr;
        amrex::ParallelFor(
          tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            tag_abserror(i, j, k, tag_arr, S_derarr, captured_velerr, tagval);
          });
      }
      if (level < tagging_parm->max_velgrad_lev) {
        const amrex::Real captured_velgrad = tagging_parm->velgrad;
        amrex::ParallelFor(
          tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            tag_graderror(i, j, k, tag_arr, S_derarr, captured_velgrad, tagval);
          });
      }

      // Tagging vel_y
      S_derData.setVal<amrex::RunOn::Device>(0.0);
      pc_dervely(
        datbox, S_derData, ncp, Sfab.nComp(), S_data[mfi], geom, time, bc,
        level);
      if (level < tagging_parm->max_velerr_lev) {
        const amrex::Real captured_velerr = tagging_parm->velerr;
        amrex::ParallelFor(
          tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            tag_abserror(i, j, k, tag_arr, S_derarr, captured_velerr, tagval);
          });
      }
      if (level < tagging_parm->max_velgrad_lev) {
        const amrex::Real captured_velgrad = tagging_parm->velgrad;
        amrex::ParallelFor(
          tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            tag_graderror(i, j, k, tag_arr, S_derarr, captured_velgrad, tagval);
          });
      }

      // Tagging vel_z
      S_derData.setVal<amrex::RunOn::Device>(0.0);
      pc_dervelz(
        datbox, S_derData, ncp, Sfab.nComp(), S_data[mfi], geom, time, bc,
        level);
      if (level < tagging_parm->max_velerr_lev) {
        const amrex::Real captured_velerr = tagging_parm->velerr;
        amrex::ParallelFor(
          tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            tag_abserror(i, j, k, tag_arr, S_derarr, captured_velerr, tagval);
          });
      }
      if (level < tagging_parm->max_velgrad_lev) {
        const amrex::Real captured_velgrad = tagging_parm->velgrad;
        amrex::ParallelFor(
          tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            tag_graderror(i, j, k, tag_arr, S_derarr, captured_velgrad, tagval);
          });
      }

      // Tagging magnitude of vorticity
      S_derData.setVal<amrex::RunOn::Device>(0.0);
      pc_dermagvort(
        tilebox, S_derData, ncp, Sfab.nComp(), S_data[mfi], geom, time, bc,
        level);
      if (level < tagging_parm->max_vorterr_lev) {
        const amrex::Real vorterr =
          tagging_parm->vorterr * std::pow(2.0, level);
        amrex::ParallelFor(
          tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            tag_abserror(i, j, k, tag_arr, S_derarr, vorterr, tagval);
          });
      }

      // Tagging temperature
      S_derData.setVal<amrex::RunOn::Device>(0.0);
      pc_dertemp(
        datbox, S_derData, ncp, Sfab.nComp(), S_data[mfi], geom, time, bc,
        level);
      if (level < tagging_parm->max_temperr_lev) {
        const amrex::Real captured_temperr = tagging_parm->temperr;
        amrex::ParallelFor(
          tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            tag_error(i, j, k, tag_arr, S_derarr, captured_temperr, tagval);
          });
      }
      if (level < tagging_parm->max_tempgrad_lev) {
        const amrex::Real captured_tempgrad = tagging_parm->tempgrad;
        amrex::ParallelFor(
          tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            tag_graderror(
              i, j, k, tag_arr, S_derarr, captured_tempgrad, tagval);
          });
      }

      // Tagging flame tracer
      if (!flame_trac_name.empty()) {
        int idx = -1;
        for (int i = 0; i < spec_names.size(); ++i) {
          if (flame_trac_name == spec_names[i]) {
            idx = i;
          }
        }

        if (idx >= 0) {
          // const std::string name = "Y("+flame_trac_name+")";
          // if (amrex::ParallelDescriptor::IOProcessor())
          // amrex::Print() << " Flame tracer will be " << name << '\n';

          S_derData.setVal<amrex::RunOn::Device>(0.0);
          pc_derspectrac(
            datbox, S_derData, ncp, Sfab.nComp(), S_data[mfi], geom, time, bc,
            level, idx);

          if (level < tagging_parm->max_ftracerr_lev) {
            const amrex::Real captured_ftracerr = tagging_parm->ftracerr;
            amrex::ParallelFor(
              tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                tag_error(
                  i, j, k, tag_arr, S_derarr, captured_ftracerr, tagval);
              });
          }
          if (level < tagging_parm->max_ftracgrad_lev) {
            const amrex::Real captured_ftracgrad = tagging_parm->ftracgrad;
            amrex::ParallelFor(
              tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                tag_graderror(
                  i, j, k, tag_arr, S_derarr, captured_ftracgrad, tagval);
              });
          }

        } else {
          amrex::Abort("Unknown species identified as flame_trac_name");
        }
      }

#ifdef PELEC_USE_EB
      // Tagging volume fraction
      if (level < tagging_parm->max_vfracerr_lev) {
        amrex::ParallelFor(
          tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            tag_error_bounds(i, j, k, tag_arr, vfrac_arr, 0.0, 1.0, tagval);
          });
      }
#endif

      // Problem specific tagging
      const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx =
        geom.CellSizeArray();
      const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_lo =
        geom.ProbLoArray();
      const auto captured_level = level;
      amrex::ParallelFor(
        tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          set_problem_tags<ProblemTags>(
            i, j, k, tag_arr, Sfab, tagval, dx, prob_lo, time, captured_level);
        });

      // Now update the tags in the TagBox.
      // tag_arr.tags(itags, tilebox);
    }
  }
}

std::unique_ptr<amrex::MultiFab>
PeleC::derive(const std::string& name, amrex::Real time, int ngrow)
{

  if ((do_les) && (name == "C_s2")) {
    std::unique_ptr<amrex::MultiFab> derive_dat(
      new amrex::MultiFab(grids, dmap, 1, 0));
    amrex::MultiFab::Copy(*derive_dat, LES_Coeffs, comp_Cs2, 0, 1, 0);
    return derive_dat;
  }
  if ((do_les) && (name == "C_I")) {
    std::unique_ptr<amrex::MultiFab> derive_dat(
      new amrex::MultiFab(grids, dmap, 1, 0));
    amrex::MultiFab::Copy(*derive_dat, LES_Coeffs, comp_CI, 0, 1, 0);
    return derive_dat;
  }
  if ((do_les) && (les_model != 1) && (name == "Pr_T")) {
    std::unique_ptr<amrex::MultiFab> derive_dat(
      new amrex::MultiFab(grids, dmap, 1, 0));
    amrex::MultiFab::Copy(*derive_dat, LES_Coeffs, comp_PrT, 0, 1, 0);
    return derive_dat;
  }
  if ((do_les) && (les_model == 1) && (name == "Pr_T")) {
    std::unique_ptr<amrex::MultiFab> derive_dat(
      new amrex::MultiFab(grids, dmap, 1, 0));
    amrex::MultiFab::Copy(*derive_dat, LES_Coeffs, comp_Cs2ovPrT, 0, 1, 0);
    amrex::MultiFab::Divide(*derive_dat, LES_Coeffs, comp_Cs2, 0, 1, 0);
    return derive_dat;
  }

#ifdef PELEC_USE_EB
  if (name == "vfrac") {
    std::unique_ptr<amrex::MultiFab> mf(
      new amrex::MultiFab(grids, dmap, 1, ngrow, amrex::MFInfo()));
    if (ngrow > 0) {
      mf->setBndry(0);
    }
    amrex::MultiFab::Copy(*mf, vfrac, 0, 0, 1, 0);
    return mf;
  }
#endif

#ifdef AMREX_PARTICLES
  return particleDerive(name, time, ngrow);
#else
  return AmrLevel::derive(name, time, ngrow);
#endif
}

void
PeleC::derive(
  const std::string& name, amrex::Real time, amrex::MultiFab& mf, int dcomp)
{
#ifdef PELEC_USE_EB
  if (name == "vfrac") {
    amrex::MultiFab::Copy(mf, vfrac, 0, dcomp, 1, 0);
  } else
#endif
  {
    AmrLevel::derive(name, time, mf, dcomp);
  }
}

void
PeleC::clear_prob()
{
  pc_prob_close();
}

#ifdef PELEC_USE_REACTIONS
void
PeleC::init_reactor()
{
#ifdef USE_SUNDIALS_PP
#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    reactor_init(1, 1);
  }
#endif
}

void
PeleC::close_reactor()
{
}
#endif

void
PeleC::init_les()
{
  // Fill with default coefficient values
  LES_Coeffs.define(grids, dmap, nCompC, 1);
  LES_Coeffs.setVal(0.0);
  LES_Coeffs.setVal(Cs * Cs, comp_Cs2, 1, LES_Coeffs.nGrow());
  LES_Coeffs.setVal(CI, comp_CI, 1, LES_Coeffs.nGrow());
  if (les_model == 1) {
    LES_Coeffs.setVal(Cs * Cs * PrT, comp_Cs2ovPrT, 1, LES_Coeffs.nGrow());
  } else {
    LES_Coeffs.setVal(PrT, comp_PrT, 1, LES_Coeffs.nGrow());
  }

  amrex::Print() << "WARNING: LES with Fuego assumes Cp is a weak function of T"
                 << std::endl;
  if (NUM_SPECIES > 2) {
    amrex::Abort("LES is not supported for multi-component systems");
  } else if (NUM_SPECIES == 2) {
    amrex::Print()
      << "WARNING: LES is not supported for multi-component systems"
      << std::endl;
  }
#ifdef PELEC_USE_SRK
  amrex::Abort("LES is not supported for non-ideal equations of state");
#endif
}

void
PeleC::init_filters()
{
  if (level > 0) {
    amrex::IntVect ref_ratio = parent->refRatio(level - 1);
    les_filter = Filter(
      les_filter_type,
      static_cast<int>(les_filter_fgr * std::pow(ref_ratio[0], level)));
  } else {
    les_filter = Filter(les_filter_type, les_filter_fgr);
  }

  nGrowF = les_filter.get_filter_ngrow();

  // Add grow cells necessary for explicit filtering of source terms
  if (do_hydro) {
    Sborder.define(
      grids, dmap, NVAR, Sborder.nGrow() + nGrowF, amrex::MFInfo(), Factory());
    hydro_source.define(
      grids, dmap, NVAR, hydro_source.nGrow() + nGrowF, amrex::MFInfo(),
      Factory());
    sources_for_hydro.define(
      grids, dmap, NVAR, sources_for_hydro.nGrow() + nGrowF, amrex::MFInfo(),
      Factory());
  }

  volume.clear();
  volume.define(
    grids, dmap, 1, numGrow() + nGrowF, amrex::MFInfo(),
    amrex::FArrayBoxFactory());
  geom.GetVolume(volume);

  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    area[dir].clear();
    area[dir].define(
      getEdgeBoxArray(dir), dmap, 1, numGrow() + nGrowF, amrex::MFInfo(),
      amrex::FArrayBoxFactory());
    geom.GetFaceArea(area[dir], dir);
  }
}
#ifdef PELEC_USE_MASA
void
PeleC::init_mms()
{
  if (!mms_initialized) {
    if (verbose && amrex::ParallelDescriptor::IOProcessor()) {
      amrex::Print() << "Initializing MMS" << std::endl;
    }
// Shut of FPE for MASA initialization because it has FPEs
#ifdef PELEC_ENABLE_FPE_TRAP
#if defined(__linux__)
    unsigned int prev_fpe_excepts = fegetexcept();
    fedisableexcept(prev_fpe_excepts);
#elif defined(__APPLE__)
    static fenv_t prev_fpe_excepts;
    fegetenv(&prev_fpe_excepts);
    static fenv_t new_fpe_excepts;
    new_fpe_excepts.__control |= FE_ALL_EXCEPT;
    new_fpe_excepts.__mxcsr |= FE_ALL_EXCEPT << 7;
    fesetenv(&new_fpe_excepts);
#endif
#endif
    masa_init("mms", masa_solution_name.c_str());
    masa_set_param("Cs", PeleC::Cs);
    masa_set_param("CI", PeleC::CI);
    masa_set_param("PrT", PeleC::PrT);
    mms_initialized = true;
#ifdef PELEC_ENABLE_FPE_TRAP
#if defined(__linux__)
    if (prev_fpe_excepts != 0) {
      feenableexcept(prev_fpe_excepts);
    }
#elif defined(__APPLE__)
    fesetenv(&prev_fpe_excepts);
#endif
#endif
  }
}
#endif

void
PeleC::reset_internal_energy(amrex::MultiFab& S_new, int ng)
{
#ifndef AMREX_USE_GPU
  amrex::Real sum0 = 0.0;
  if (parent->finestLevel() == 0 && print_energy_diagnostics) {
    // Pass in the multifab and the component
    sum0 = volWgtSumMF(S_new, Eden, true);
  }
#endif
  // Ensure (rho e) isn't too small or negative
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  const auto captured_allow_small_energy = allow_small_energy;
  const auto captured_allow_negative_energy = allow_negative_energy;
  const auto captured_dual_energy_update_E_from_e = dual_energy_update_E_from_e;
  const auto captured_verbose = verbose;
  const auto captured_dual_energy_eta2 = dual_energy_eta2;
  for (amrex::MFIter mfi(S_new, amrex::TilingIfNotGPU()); mfi.isValid();
       ++mfi) {
    const amrex::Box& bx = mfi.growntilebox(ng);
    const auto& sarr = S_new.array(mfi);
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      pc_rst_int_e(
        i, j, k, sarr, captured_allow_small_energy,
        captured_allow_negative_energy, captured_dual_energy_update_E_from_e,
        captured_dual_energy_eta2, captured_verbose);
    });
  }

#ifndef AMREX_USE_GPU
  if (parent->finestLevel() == 0 && print_energy_diagnostics) {
    // Pass in the multifab and the component
    amrex::Real sum = volWgtSumMF(S_new, Eden, true);
#ifdef AMREX_LAZY
    Lazy::QueueReduction([=]() mutable {
#endif
      amrex::ParallelDescriptor::ReduceRealSum(sum0);
      amrex::ParallelDescriptor::ReduceRealSum(sum);
      if (
        amrex::ParallelDescriptor::IOProcessor() &&
        amrex::Math::abs(sum - sum0) > 0) {
        amrex::Print() << "(rho E) added from reset terms                 : "
                       << sum - sum0 << " out of " << sum0 << std::endl;
      }
#ifdef AMREX_LAZY
    });
#endif
  }
#endif
}

void
PeleC::computeTemp(amrex::MultiFab& S, int ng)
{
  reset_internal_energy(S, ng);

#ifdef PELEC_USE_EB
  auto const& fact =
    dynamic_cast<amrex::EBFArrayBoxFactory const&>(S.Factory());
  auto const& flags = fact.getMultiEBCellFlagFab();
#endif

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  for (amrex::MFIter mfi(S, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    const amrex::Box& bx = mfi.growntilebox(ng);

#ifdef PELEC_USE_EB
    const auto& flag_fab = flags[mfi];
    amrex::FabType typ = flag_fab.getType(bx);
    if (typ == amrex::FabType::covered) {
      continue;
    }
#endif

    const auto& sarr = S.array(mfi);
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      pc_cmpTemp(i, j, k, sarr);
    });
  }
}

amrex::Real
PeleC::getCPUTime()
{
  int numCores = amrex::ParallelDescriptor::NProcs();
#ifdef _OPENMP
  numCores = numCores * omp_get_max_threads();
#endif

  amrex::Real T =
    numCores * (amrex::ParallelDescriptor::second() - startCPUTime) +
    previousCPUTimeUsed;

  return T;
}

amrex::MultiFab&
PeleC::build_fine_mask()
{
  // Mask for zeroing covered cells
  AMREX_ASSERT(level > 0);

  if (!fine_mask.empty()) {
    return fine_mask;
  }

  const amrex::BoxArray& cba = parent->boxArray(level - 1);
  const amrex::DistributionMapping& cdm = parent->DistributionMap(level - 1);

  fine_mask.define(cba, cdm, 1, 0, amrex::MFInfo(), amrex::FArrayBoxFactory());
  fine_mask.setVal(1.0);

  amrex::BoxArray fba = parent->boxArray(level);
  amrex::iMultiFab ifine_mask = makeFineMask(cba, cdm, fba, crse_ratio, 1, 0);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  for (amrex::MFIter mfi(fine_mask, amrex::TilingIfNotGPU()); mfi.isValid();
       ++mfi) {
    auto& fab = fine_mask[mfi];
    auto& ifab = ifine_mask[mfi];
    const auto arr = fab.array();
    const auto iarr = ifab.array();
    amrex::ParallelFor(
      fab.box(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
#ifdef _OPENMP
#pragma omp atomic write
#endif
        arr(i, j, k) = iarr(i, j, k);
      });
  }

  return fine_mask;
}

const amrex::iMultiFab*
PeleC::build_interior_boundary_mask(int ng)
{
  for (int i = 0; i < ib_mask.size(); ++i) {
    if (ib_mask[i]->nGrow() == ng) {
      return ib_mask[i].get();
    }
  }

  // If we got here, we need to build a new one
  if (ib_mask.empty()) {
    ib_mask.resize(0);
  }

  ib_mask.push_back(std::make_unique<amrex::iMultiFab>(
    grids, dmap, 1, ng, amrex::MFInfo(),
    amrex::DefaultFabFactory<amrex::IArrayBox>()));

  amrex::iMultiFab* imf = ib_mask.back().get();
  int ghost_covered_by_valid = 0;
  int other_cells =
    1; // uncovered ghost, valid, and outside domain cells are set to 1

  imf->BuildMask(
    geom.Domain(), geom.periodicity(), ghost_covered_by_valid, other_cells,
    other_cells, other_cells);

  return imf;
}

amrex::Real
PeleC::clean_state(amrex::MultiFab& S)
{
  // Enforce a minimum density.

  amrex::MultiFab temp_state(
    S.boxArray(), S.DistributionMap(), S.nComp(), S.nGrow(), amrex::MFInfo(),
    Factory());

  amrex::MultiFab::Copy(temp_state, S, 0, 0, S.nComp(), S.nGrow());

  amrex::Real frac_change_t = enforce_min_density(temp_state, S);

  // normalize_species(S);

  return frac_change_t;
}

amrex::Real
PeleC::clean_state(const amrex::MultiFab& S, amrex::MultiFab& S_old)
{
  // Enforce a minimum density.

  amrex::Real frac_change_t = enforce_min_density(S_old, S);

  // normalize_species(S);

  return frac_change_t;
}

#ifdef PELEC_USE_EB
void
PeleC::InitialRedistribution()
{
  BL_PROFILE("PeleC::InitialRedistribution()");

  // Next we must redistribute the initial solution if we are going to use
  // MergeRedist or StateRedist redistribution schemes
  if (
    (eb_in_domain) && ((redistribution_type != "StateRedist") &&
                       (redistribution_type != "MergeRedist"))) {
    return;
  }

  if (redistribution_type == "MergeRedist") {
    amrex::Abort("MergeRedist is unsupported. Check with AMReX-Hydro if that "
                 "has been fixed");
  }

  if (verbose) {
    amrex::Print() << "Doing initial redistribution... " << std::endl;
  }

  // Initial data are set at new time step
  amrex::MultiFab& S_new = get_new_data(State_Type);
  amrex::MultiFab tmp(
    grids, dmap, S_new.nComp(), numGrow(), amrex::MFInfo(), Factory());

  amrex::MultiFab::Copy(tmp, S_new, 0, 0, S_new.nComp(), S_new.nGrow());
  const amrex::Real time = state[State_Type].curTime();
  FillPatch(*this, tmp, numGrow(), time, State_Type, 0, S_new.nComp());
  EB_set_covered(tmp, 0.0);

  const amrex::StateDescriptor* desc = state[State_Type].descriptor();
  const auto& bcs = desc->getBCs();
  amrex::Gpu::DeviceVector<amrex::BCRec> d_bcs(desc->nComp());
  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, bcs.begin(), bcs.end(), d_bcs.begin());

  for (amrex::MFIter mfi(S_new, amrex::TilingIfNotGPU()); mfi.isValid();
       ++mfi) {
    const amrex::Box& bx = mfi.validbox();

    auto const& fact =
      dynamic_cast<amrex::EBFArrayBoxFactory const&>(S_new.Factory());

    auto const& flags = fact.getMultiEBCellFlagFab()[mfi];
    amrex::Array4<const amrex::EBCellFlag> const& flag_arr =
      flags.const_array();

    if (
      (flags.getType(amrex::grow(bx, 1)) != amrex::FabType::covered) &&
      (flags.getType(amrex::grow(bx, 1)) != amrex::FabType::regular)) {
      amrex::Array4<const amrex::Real> AMREX_D_DECL(fcx, fcy, fcz), ccc,
        AMREX_D_DECL(apx, apy, apz);

      AMREX_D_TERM(fcx = fact.getFaceCent()[0]->const_array(mfi);
                   , fcy = fact.getFaceCent()[1]->const_array(mfi);
                   , fcz = fact.getFaceCent()[2]->const_array(mfi););

      ccc = fact.getCentroid().const_array(mfi);

      AMREX_D_TERM(apx = fact.getAreaFrac()[0]->const_array(mfi);
                   , apy = fact.getAreaFrac()[1]->const_array(mfi);
                   , apz = fact.getAreaFrac()[2]->const_array(mfi););

      Redistribution::ApplyToInitialData(
        bx, NVAR, S_new.array(mfi), tmp.array(mfi), flag_arr,
        AMREX_D_DECL(apx, apy, apz), vfrac.const_array(mfi),
        AMREX_D_DECL(fcx, fcy, fcz), ccc, d_bcs.dataPtr(), geom,
        redistribution_type);
    }
  }
  set_body_state(S_new);
}
#endif
