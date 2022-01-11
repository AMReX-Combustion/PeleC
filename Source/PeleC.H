#ifndef _PELEC_H_
#define _PELEC_H_

#include <AMReX_BC_TYPES.H>
#include <AMReX_AmrLevel.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_ParmParse.H>

#ifdef PELEC_USE_EB
#include <AMReX_EBFArrayBox.H>
#include "EBStencilTypes.H"
#include "SparseData.H"
#include <AMReX_EBFluxRegister.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_MultiCutFab.H>
#else
#include <AMReX_YAFluxRegister.H>
#endif

#ifdef AMREX_LAZY
#include <AMReX_Lazy.H>
#endif

#ifdef AMREX_PARTICLES
#include <AMReX_Particles.H>
#include <SprayParticles.H>
#endif

#include "Filter.H"
#include "Tagging.H"
#include "IndexDefines.H"
#include "prob_parm.H"
#include "PelePhysics.H"
#include "ReactorBase.H"

#define UserBC 6

enum StateType { State_Type = 0, Reactions_Type, Work_Estimate_Type };

// Create storage for all source terms.

enum sources {
  ext_src = 0,
  forcing_src,
#ifdef AMREX_PARTICLES
  spray_src,
#endif
  diff_src,
  les_src,
#ifdef PELEC_USE_MASA
  mms_src,
#endif
  num_src
};

/*
static amrex::Box
the_same_box(const amrex::Box& b)
{
  return b;
}

static amrex::Box
grow_box_by_one(const amrex::Box& b)
{
  return amrex::grow(b, 1);
}
*/

// AmrLevel-derived class for hyperbolic conservation equations for stellar
// media

class PeleC : public amrex::AmrLevel
{
public:
  // Default constructor.  Builds invalid object.
  PeleC();

  PeleC(
    amrex::Amr& papa,
    int lev,
    const amrex::Geometry& level_geom,
    const amrex::BoxArray& bl,
    const amrex::DistributionMapping& dm,
    amrex::Real time);

  virtual ~PeleC();

  // Restart from a checkpoint file.
  virtual void restart(
    amrex::Amr& papa, std::istream& is, bool bReadSpecial = false) override;

  // This is called only when we restart from an old checkpoint.
  virtual void
  set_state_in_checkpoint(amrex::Vector<int>& state_in_checkpoint) override;
  bool check_state_in_checkpoint(const StateType state_type);

  // Call AmrLevel::checkPoint and then add radiation info
  virtual void checkPoint(
    const std::string& dir,
    std::ostream& os,
    amrex::VisMF::How how,
    bool dump_old) override;

  virtual void setPlotVariables() override;

  // Write a plotfile to specified directory.
  virtual void writePlotFile(
    const std::string& dir, std::ostream& os, amrex::VisMF::How how) override;
  virtual void writeSmallPlotFile(
    const std::string& dir, std::ostream& os, amrex::VisMF::How how) override;
  void writeJobInfo(const std::string& dir);
  static void writeBuildInfo(std::ostream& os);

  // Define data descriptors.
  static void variableSetUp();

  // Cleanup data descriptors at end of run.
  static void variableCleanUp();

  // Initialize grid data at problem start-up.
  virtual void initData() override;

#ifdef AMREX_PARTICLES

  // Read particle-related inputs
  static void readParticleParams();

  // Define particle related variables
  static void defineParticles();

  // Initialize particle locations and velocities (and strengths if relevant)
  void initParticles();

  // Timestamp particles
  void particleTimestamp(int ngrow);

  // Default verbosity of Particle class
  static int particle_verbose;

  // How to initialize at restart
  void particlePostRestart(
    const std::string& restart_file, bool is_checkpoint = true);

  // Redistribute
  virtual void particleRedistribute(
    int lbase = 0, int nGrow = 0, int local = 0, bool init_part = false);

  // Setup virtual particles if necessary
  void setupVirtualParticles();

  // Remove virtual particles if necessary
  void removeVirtualParticles();

  // Setup ghost particles (for finer levels) if necessary
  void setupGhostParticles(int ngrow);

  // Remove ghost particles (for this level) if necessary
  void removeGhostParticles();

  // Time step control based on particles
  void particleEstTimeStep(amrex::Real& est_dt);

  // Derived quantities associated with particles
  std::unique_ptr<amrex::MultiFab>
  particleDerive(const std::string& name, amrex::Real time, int ngrow);

  // Default cfl of particles in Particle class
  static amrex::Real particle_cfl;

  // Should we write particles into plotfiles?
  static int write_particle_plotfiles;

  // Should we write particle ascii files?
  static int write_spray_ascii_files;

  void setSprayGridInfo(
    const int amr_iteration,
    const int amr_ncycle,
    int& ghost_width,
    int& where_width,
    int& spray_n_grow,
    int& tmp_src_width);

  static int particle_mass_tran;
  static int particle_heat_tran;
  static int particle_mom_tran;
  static amrex::Vector<std::string> sprayFuelNames;
  static amrex::Real sprayRefT;
#endif

  // Get problem metrics.
  amrex::MultiFab* Area();
  amrex::MultiFab& Area(int dir);

  // The volume of the multifab.
  amrex::MultiFab& Volume();

  // Set time levels of state data.
  virtual void setTimeLevel(
    amrex::Real time, amrex::Real dt_old, amrex::Real dt_new) override;

  // Initialize data on this level from another PeleC (during regrid).
  virtual void init(amrex::AmrLevel& old) override;

  // Initialize data on this level after regridding if old level did not
  // previously exist
  virtual void init() override;

  // Initialize EB geometry for finest_level and level grids for
  // other levels for the Amr class to do timed load balances.
  virtual int WorkEstType() override { return Work_Estimate_Type; }

#ifdef PELEC_USE_EB
  static bool DoMOLLoadBalance() { return do_mol_load_balance; }

  const amrex::MultiFab& volFrac() const { return vfrac; }

  void init_eb(
    const amrex::Geometry& level_geom,
    const amrex::BoxArray& ba,
    const amrex::DistributionMapping& dm);

  void initialize_eb2_structs();

  void define_body_state();

  void set_body_state(amrex::MultiFab& S);

  void zero_in_body(amrex::MultiFab& S) const;
#endif

  void set_typical_values_chem();

  // Proceed with next timestep?
  virtual int okToContinue() override;

  // Advance grids at this level in time.
  virtual amrex::Real
  advance(amrex::Real time, amrex::Real dt, int iteration, int ncycle) override;

  amrex::Real do_mol_advance(
    amrex::Real time, amrex::Real dt, int amr_iteration, int amr_ncycle);

  amrex::Real do_sdc_advance(
    amrex::Real time, amrex::Real dt, int amr_iteration, int amr_ncycle);

  void initialize_sdc_advance(
    amrex::Real time, amrex::Real dt, int amr_iteration, int amr_ncycle);

  void finalize_sdc_advance(
    amrex::Real time, amrex::Real dt, int amr_iteration, int amr_ncycle);

  amrex::Real do_sdc_iteration(
    amrex::Real time,
    amrex::Real dt,
    int amr_iteration,
    int amr_ncycle,
    int sub_iteration,
    int sub_ncycle);

  void initialize_sdc_iteration(
    amrex::Real time,
    amrex::Real dt,
    int amr_iteration,
    int amr_ncycle,
    int sdc_iteration,
    int sdc_ncycle);

  static void finalize_sdc_iteration(
    amrex::Real time,
    amrex::Real dt,
    int amr_iteration,
    int amr_ncycle,
    int sdc_iteration,
    int sdc_ncycle);

  void construct_Snew(
    amrex::MultiFab& S_new, const amrex::MultiFab& S_old, amrex::Real dt);

  void construct_hydro_source(
    const amrex::MultiFab& S,
    amrex::Real time,
    amrex::Real dt,
    int amr_iteration,
    int amr_ncycle,
    int sub_iteration,
    int sub_ncycle);

  void construct_old_source(
    int src,
    amrex::Real time,
    amrex::Real dt,
    int amr_iteration,
    int amr_ncycle,
    int sub_iteration,
    int sub_ncycle);

  void construct_new_source(
    int src,
    amrex::Real time,
    amrex::Real dt,
    int amr_iteration,
    int amr_ncycle,
    int sub_iteration,
    int sub_ncycle);

  void sum_of_sources(amrex::MultiFab& source);

  void construct_old_ext_source(amrex::Real time, amrex::Real dt);

  void construct_new_ext_source(amrex::Real time, amrex::Real dt);

  static void fill_ext_source(
    amrex::Real time,
    amrex::Real dt,
    const amrex::MultiFab& state_old,
    const amrex::MultiFab& S_new,
    amrex::MultiFab& ext_src,
    int ng);

  void construct_old_forcing_source(amrex::Real time, amrex::Real dt);

  void construct_new_forcing_source(amrex::Real time, amrex::Real dt);

  static void fill_forcing_source(
    const amrex::MultiFab& state_old,
    const amrex::MultiFab& state_new,
    amrex::MultiFab& forcing_src,
    int ng);

#ifdef PELEC_USE_MASA
  void construct_old_mms_source(amrex::Real time);

  void construct_new_mms_source(amrex::Real time);

  void fill_mms_source(
    amrex::Real time,
    const amrex::MultiFab& S,
    amrex::MultiFab& mms_src,
    int ng);
#endif

  void getLESTerm(
    amrex::Real time,
    amrex::Real dt,
    amrex::MultiFab& LESTerm,
    amrex::Real flux_factor);

  void getSmagorinskyLESTerm(
    amrex::Real time,
    amrex::Real dt,
    amrex::MultiFab& LESTerm,
    amrex::Real flux_factor);

  void getDynamicSmagorinskyLESTerm(
    amrex::Real time,
    amrex::Real dt,
    amrex::MultiFab& LESTerm,
    amrex::Real flux_factor);

  void construct_old_les_source(
    amrex::Real time, amrex::Real dt, int sub_iteration, int sub_ncycle);

  void construct_new_les_source(
    amrex::Real time, amrex::Real dt, int sub_iteration, int sub_ncycle);

  static void set_active_sources();

  // Estimate time step.
  amrex::Real estTimeStep(amrex::Real dt_old);

  // Compute initial time step.
  amrex::Real initialTimeStep();

  // Compute initial `dt'.
  virtual void computeInitialDt(
    int finest_level,
    int sub_cycle,
    amrex::Vector<int>& n_cycle,
    const amrex::Vector<amrex::IntVect>& ref_ratio,
    amrex::Vector<amrex::Real>& dt_level,
    amrex::Real stop_time) override;

  // Compute new `dt'.
  virtual void computeNewDt(
    int finest_level,
    int sub_cycle,
    amrex::Vector<int>& n_cycle,
    const amrex::Vector<amrex::IntVect>& ref_ratio,
    amrex::Vector<amrex::Real>& dt_min,
    amrex::Vector<amrex::Real>& dt_level,
    amrex::Real stop_time,
    int post_regrid_flag) override;

  // Allocate data at old time.
  virtual void allocOldData() override;

  // Remove data at old time.
  virtual void removeOldData() override;

  // Passes some data about the grid to a Fortran module.
  void setGridInfo();

  // Print information about energy budget.
  void do_energy_diagnostics();

  // Do work after timestep().
  virtual void post_timestep(int iteration) override;

  // Contains operations to be done only after a full coarse timestep.
  virtual void postCoarseTimeStep(amrex::Real cumtime) override;

  // Do work after regrid().
  virtual void post_regrid(int lbase, int new_finest) override;

  // Do work after a restart().
  virtual void post_restart() override;

  // Do work after init().
  virtual void post_init(amrex::Real stop_time) override;

  // Error estimation for regridding.
  virtual void errorEst(
    amrex::TagBoxArray& tags,
    int clearval,
    int tagval,
    amrex::Real time,
    int n_error_buf = 0,
    int ngrow = 0) override;

  // Returns a MultiFab containing the derived data for this level.
  // The user is responsible for deleting this pointer when done
  // with it.  If ngrow>0 the MultiFab is built on the appropriately
  // grown BoxArray.
  virtual std::unique_ptr<amrex::MultiFab>
  derive(const std::string& name, amrex::Real time, int ngrow) override;

  // This version of derive() fills the dcomp'th component of mf with the
  // derived quantity.
  virtual void derive(
    const std::string& name,
    amrex::Real time,
    amrex::MultiFab& mf,
    int dcomp) override;

  static int numGrow();

  void react_state(
    amrex::Real time,
    amrex::Real dt,
    bool init = false,
    amrex::MultiFab* aux_src = nullptr);

  void reset_internal_energy(amrex::MultiFab& S_new, int ng);

  void computeTemp(amrex::MultiFab& State, int ng);

  void getMOLSrcTerm(
    const amrex::MultiFab& S,
    amrex::MultiFab& MOLSrcTerm,
    amrex::Real time,
    amrex::Real dt,
    amrex::Real flux_factor);

  static void enforce_consistent_e(amrex::MultiFab& S);

  amrex::Real volWgtSum(
    const std::string& name,
    amrex::Real time,
    bool local = false,
    bool finemask = true);
  amrex::Real volWgtSquaredSum(
    const std::string& name, amrex::Real time, bool local = false);
  amrex::Real volWgtSumMF(
    const amrex::MultiFab& mf,
    int comp,
    bool local = false,
    bool finemask = false);
  amrex::Real
  volWgtSquaredSumDiff(int comp, amrex::Real time, bool local = false);
  amrex::Real
  maxDerive(const std::string& name, amrex::Real time, bool local = false);
  amrex::Real
  minDerive(const std::string& name, amrex::Real time, bool local = false);

  // derives that need variables part of this class (e.g. trans_parm)
  static void pc_derviscosity(
    const amrex::Box& bx,
    amrex::FArrayBox& derfab,
    int dcomp,
    int ncomp,
    const amrex::FArrayBox& datfab,
    const amrex::Geometry& geomdata,
    amrex::Real time,
    const int* bcrec,
    const int level);
  static void pc_derbulkviscosity(
    const amrex::Box& bx,
    amrex::FArrayBox& derfab,
    int dcomp,
    int ncomp,
    const amrex::FArrayBox& datfab,
    const amrex::Geometry& geomdata,
    amrex::Real time,
    const int* bcrec,
    const int level);
  static void pc_derconductivity(
    const amrex::Box& bx,
    amrex::FArrayBox& derfab,
    int dcomp,
    int ncomp,
    const amrex::FArrayBox& datfab,
    const amrex::Geometry& geomdata,
    amrex::Real time,
    const int* bcrec,
    const int level);
  static void pc_derdiffusivity(
    const amrex::Box& bx,
    amrex::FArrayBox& derfab,
    int dcomp,
    int ncomp,
    const amrex::FArrayBox& datfab,
    const amrex::Geometry& geomdata,
    amrex::Real time,
    const int* bcrec,
    const int level);

  // static int NVAR;
  static int Density, Xmom, Ymom, Zmom, Eden, Eint, Temp;

  static int FirstAdv, NumAdv;
  static int FirstSpec;
  static int FirstAux;

  /// Index locations into particle state
  static int pstateVel, pstateT, pstateDia, pstateRho, pstateY, pstateNum;

  // This MultiFab is on the coarser level.  This is useful for the coarser
  // level to mask out the finer level.  We only build this when it is needed.
  // This coarse MultiFab has to live on the fine level because it must be
  // updated even when only the fine level changes.
  amrex::MultiFab fine_mask;
  amrex::MultiFab& build_fine_mask();

  static bool eb_in_domain;
  AMREX_FORCE_INLINE static bool ebInDomain()
  {
    amrex::ParmParse pp("eb2");
    std::string geom_type("all_regular");
    pp.query("geom_type", geom_type);
    return (geom_type != "all_regular");
  }

#ifdef PELEC_USE_EB

  amrex::EBFluxRegister flux_reg;

  amrex::EBFluxRegister& getFluxReg();
  amrex::EBFluxRegister& getFluxReg(int lev);

  amrex::EBFluxRegister pres_reg;

  amrex::EBFluxRegister& getPresReg();
  amrex::EBFluxRegister& getPresReg(int lev);

  // static EB helper functions called from outside class
  static bool ebInitialized();
#else
  amrex::YAFluxRegister flux_reg;

  amrex::YAFluxRegister& getFluxReg();
  amrex::YAFluxRegister& getFluxReg(int lev);

  amrex::YAFluxRegister pres_reg;

  amrex::YAFluxRegister& getPresReg();
  amrex::YAFluxRegister& getPresReg(int lev);
#endif
  void avgDown();
  void avgDown(int state_indx);

  // Parameters for using WENO in the PPM method
  static int use_hybrid_weno;
  static int weno_scheme;

  static ProbParmDevice* h_prob_parm_device;
  static ProbParmDevice* d_prob_parm_device;
  static ProbParmHost* prob_parm_host;
  static TaggingParm* tagging_parm;
  static PassMap* h_pass_map;
  static PassMap* d_pass_map;
  static pele::physics::transport::TransportParams<
    pele::physics::PhysicsType::transport_type>
    trans_parms;

protected:
  amrex::iMultiFab level_mask;
#ifdef PELEC_USE_EB
  std::array<const amrex::MultiCutFab*, AMREX_SPACEDIM> areafrac;
#endif

  // Build a mask that ghost cells overlapping with interior cells in the same
  // multifab are set to 0, whereas others are set to 1.
  amrex::Vector<std::unique_ptr<amrex::iMultiFab>> ib_mask;
  const amrex::iMultiFab* build_interior_boundary_mask(int ng);

  // A state array with ghost zones.
  amrex::MultiFab Sborder;

  // Source terms to the hydrodynamics solve.
  amrex::MultiFab sources_for_hydro;

  // Source term representing hydrodynamics update.
  amrex::MultiFab hydro_source;

  // Non-hydro source terms.
  amrex::Vector<std::unique_ptr<amrex::MultiFab>> old_sources;
  amrex::Vector<std::unique_ptr<amrex::MultiFab>> new_sources;

  std::unique_ptr<pele::physics::reactions::ReactorBase> reactor;
  void init_reactor();
  void close_reactor();

  void init_les();
  void init_filters();

#ifdef PELEC_USE_MASA
  static void init_mms();
#endif

  static void clear_prob();

  static void read_params();

  static void read_tagging_params();

  PeleC& getLevel(int lev);

  void reflux();

  static void normalize_species(amrex::MultiFab& S_new);

  static amrex::Real
  enforce_min_density(amrex::MultiFab& S_old, const amrex::MultiFab& S_new);

  amrex::Real clean_state(const amrex::MultiFab& S);

  static amrex::Real
  clean_state(const amrex::MultiFab& S, amrex::MultiFab& S_old);

  void InitialRedistribution();

  void buildMetrics();

  // integrate derived quantities over domain

  amrex::Real
  sumDerive(const std::string& name, amrex::Real time, bool local = false);

  int find_datalog_index(const std::string& logname);

  void sum_integrated_quantities();

  void monitor_extrema();

  void write_info();

  void stopJob();

  // The data.
  amrex::MultiFab volume;
  amrex::MultiFab area[AMREX_SPACEDIM];
  amrex::MultiFab dLogArea[1];
  amrex::Vector<amrex::Vector<amrex::Real>> radius;
  // int m_ngrow;

  // Static data members.
#include "pelec_params.H"

  static bool signalStopJob;
  static bool dump_old;
  static int radius_grow;
  static int verbose;
  static amrex::BCRec phys_bc;

  // for keeping track of mass changes from negative density resets
  static amrex::Real frac_change;

  // For keeping track of fluid quantities lost at physical grid boundaries.
  // This should persist through restarts, but right now only on level 0.
  static const int n_lost = 8;
  amrex::Real material_lost_through_boundary_cumulative[n_lost];
  amrex::Real material_lost_through_boundary_temp[n_lost];

  // for keeping track of the amount of CPU time used -- this will persist
  // after restarts
  static amrex::Real previousCPUTimeUsed;
  static amrex::Real startCPUTime;

  static amrex::Real getCPUTime();

  static int diffuse_temp;
  static int diffuse_enth;
  static int diffuse_spec;
  static int diffuse_vel;
  static amrex::Real diffuse_cutoff_density;
  static bool do_diffuse;

  static int les_model;
  static int les_filter_type;
  static int les_filter_fgr;
  Filter les_filter;
  int nGrowF;
  static int les_test_filter_type;
  static int les_test_filter_fgr;
  amrex::MultiFab LES_Coeffs;
  amrex::MultiFab filtered_les_source;

#ifdef PELEC_USE_MASA
  static bool mms_initialized;
  bool mms_src_evaluated;
  amrex::MultiFab mms_source;
#endif

  // bool FillPatchedOldState_ok;

  // There can be only one Diffusion object, it covers all levels:
  static class Diffusion* diffusion;

#ifdef AMREX_PARTICLES
  static int do_spray_particles;
  static SprayParticleContainer* theSprayPC();
  static SprayParticleContainer* theVirtPC();
  static SprayParticleContainer* theGhostPC();
#endif

  static int num_state_type;

  static amrex::Vector<std::string> spec_names;

  static amrex::Vector<int> src_list;

  static bool use_typical_vals_chem;
  static bool use_typical_vals_chem_usr;
  static amrex::Real typical_rhoY_val_min;
  static int reset_typical_vals_int;
  static amrex::Vector<amrex::Real> typical_values_chem_usr;

  // problem-specific includes
#include <Problem.H>

#ifdef PELEC_USE_EB

  static amrex::GpuArray<amrex::Real, NVAR> body_state;
  static bool body_state_set;
  static bool eb_initialized;

  amrex::MultiFab vfrac;

  amrex::Vector<amrex::Gpu::DeviceVector<EBBndryGeom>> sv_eb_bndry_geom;
  amrex::Vector<amrex::Gpu::DeviceVector<EBBndrySten>> sv_eb_bndry_grad_stencil;
  amrex::
    GpuArray<amrex::Vector<amrex::Gpu::DeviceVector<FaceSten>>, AMREX_SPACEDIM>
      flux_interp_stencil;

  amrex::Vector<SparseData<amrex::Real, EBBndrySten>> sv_eb_flux;
  amrex::Vector<SparseData<amrex::Real, EBBndrySten>> sv_eb_bcval;
#endif
  static bool do_react_load_balance;
  static bool do_mol_load_balance;
};

void pc_bcfill_hyp(
  amrex::Box const& bx,
  amrex::FArrayBox& data,
  const int dcomp,
  const int numcomp,
  amrex::Geometry const& geom,
  const amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcr,
  const int bcomp,
  const int scomp);

void pc_reactfill_hyp(
  amrex::Box const& bx,
  amrex::FArrayBox& data,
  const int dcomp,
  const int numcomp,
  amrex::Geometry const& geom,
  const amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcr,
  const int bcomp,
  const int scomp);

void pc_nullfill(
  amrex::Box const& bx,
  amrex::FArrayBox& data,
  const int dcomp,
  const int numcomp,
  amrex::Geometry const& geom,
  const amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcr,
  const int bcomp,
  const int scomp);

// Inlines

AMREX_FORCE_INLINE
int
PeleC::numGrow()
{
  int ng = 4;
  if (
    (eb_in_domain) && ((redistribution_type == "StateRedist") ||
                       (redistribution_type == "NewStateRedist"))) {
    ng += 1;
  }
  return ng;
}

AMREX_FORCE_INLINE
amrex::MultiFab*
PeleC::Area()
{
  return area;
}

AMREX_FORCE_INLINE
amrex::MultiFab&
PeleC::Area(int dir)
{
  return area[dir];
}

AMREX_FORCE_INLINE
amrex::MultiFab&
PeleC::Volume()
{
  return volume;
}

AMREX_FORCE_INLINE
PeleC&
PeleC::getLevel(int lev)
{
  return *(PeleC*)&parent->getLevel(lev);
}

AMREX_FORCE_INLINE
void
PeleC::stopJob()
{
  signalStopJob = true;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_check_initial_species(
  const int i, const int j, const int k, amrex::Array4<amrex::Real> const& Sfab)
{
  amrex::Real spec_sum = 0.0;
  for (int n = 0; n < NUM_SPECIES; n++) {
    spec_sum = spec_sum + Sfab(i, j, k, UFS + n);
  }
  if (
    amrex::Math::abs(Sfab(i, j, k, URHO) - spec_sum) >
    1.e-8 * Sfab(i, j, k, URHO)) {
    // print *,'Sum of (rho X)_i vs rho at (i,j,k):
    // ',i,j,k,spec_sum,state(i,j,k,URHO)
    amrex::Abort("Error:: Failed check of initial species summing to 1");
  }
}

#ifdef PELEC_USE_EB
AMREX_FORCE_INLINE
amrex::EBFluxRegister&
PeleC::getFluxReg()
{
  return flux_reg;
}

AMREX_FORCE_INLINE
amrex::EBFluxRegister&
PeleC::getFluxReg(int lev)
{
  return getLevel(lev).getFluxReg();
}

AMREX_FORCE_INLINE
amrex::EBFluxRegister&
PeleC::getPresReg()
{
  return pres_reg;
}

AMREX_FORCE_INLINE
amrex::EBFluxRegister&
PeleC::getPresReg(int lev)
{
  return getLevel(lev).getPresReg();
}

#else

AMREX_FORCE_INLINE
amrex::YAFluxRegister&
PeleC::getFluxReg()
{
  return flux_reg;
}

AMREX_FORCE_INLINE
amrex::YAFluxRegister&
PeleC::getFluxReg(int lev)
{
  return getLevel(lev).getFluxReg();
}

AMREX_FORCE_INLINE
amrex::YAFluxRegister&
PeleC::getPresReg()
{
  return pres_reg;
}

AMREX_FORCE_INLINE
amrex::YAFluxRegister&
PeleC::getPresReg(int lev)
{
  return getLevel(lev).getPresReg();
}
#endif

#endif
