
#ifndef WIN32
#include <unistd.h>
#endif

#include <iomanip>

#include <algorithm>
#include <cstdio>
#include <vector>
#include <iostream>
#include <string>
#include <ctime>

using std::cout;
using std::cerr;
using std::endl;
using std::istream;
using std::ostream;
using std::pair;
using std::string;

#include <AMReX_Utility.H>
#include <AMReX_CONSTANTS.H>
#include <PeleC.H>
#include <PeleC_F.H>
#include <Derive_F.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_ParmParse.H>

#include <PeleC_error_F.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBMultiFabUtil.H>
#endif

#ifdef AMREX_PARTICLES
#include <AMReX_Particles.H>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace amrex;

#ifdef USE_MASA
#include <masa.h>
using namespace MASA;
#endif

bool         PeleC::signalStopJob = false;

bool         PeleC::dump_old      = false;

int          PeleC::verbose       = 0;
int          PeleC::radius_grow   = 1;
BCRec        PeleC::phys_bc;
int          PeleC::NUM_STATE     = -1;
int          PeleC::NUM_GROW      = -1;
int          PeleC::QTHERM        = -1;
int          PeleC::QVAR          = -1;
int          PeleC::cQRHO         = -1;
int          PeleC::cQU           = -1;
int          PeleC::cQV           = -1;
int          PeleC::cQW           = -1;
int          PeleC::cQGAME        = -1;
int          PeleC::cQPRES        = -1;
int          PeleC::cQREINT       = -1;
int          PeleC::cQTEMP        = -1;
int          PeleC::cQFA          = -1;
int          PeleC::cQFS          = -1;
int          PeleC::cQFX          = -1;
int          PeleC::NQAUX         = -1;
int          PeleC::cQGAMC        = -1;
int          PeleC::cQC           = -1;
int          PeleC::cQCSML        = -1;
int          PeleC::cQDPDR        = -1;
int          PeleC::cQDPDE        = -1;
int          PeleC::cQRSPEC        = -1;

Real         PeleC::frac_change   = 1.e200;

int          PeleC::Density       = -1;
int          PeleC::Eden          = -1;
int          PeleC::Eint          = -1;
int          PeleC::Temp          = -1;
int          PeleC::Xmom          = -1;
int          PeleC::Ymom          = -1;
int          PeleC::Zmom          = -1;
int          PeleC::NumSpec       = 0;
int          PeleC::FirstSpec     = -1;

int          PeleC::NumAux        = 0;
int          PeleC::FirstAux      = -1;

int          PeleC::NumAdv        = 0;
int          PeleC::FirstAdv      = -1;

int          PeleC::pstate_loc = -1;
int          PeleC::pstate_vel = -1;
int          PeleC::pstate_T   = -1;
int          PeleC::pstate_dia = -1;
int          PeleC::pstate_rho = -1;
int          PeleC::pstate_spc = -1;
int          PeleC::n_pstate   =  0;
int          PeleC::pfld_vel   = -1;
int          PeleC::pfld_rho   = -1;
int          PeleC::pfld_T     = -1;
int          PeleC::pfld_p     = -1;
int          PeleC::pfld_spc   = -1;
int          PeleC::n_pfld     = 0;

#include <pelec_defaults.H>

int          PeleC::nGrowTr      = 4;
int          PeleC::diffuse_temp = 0;
int          PeleC::diffuse_enth = 0;
int          PeleC::diffuse_spec = 0;
int          PeleC::diffuse_vel  = 0;
Real         PeleC::diffuse_cutoff_density = -1.e200;
bool         PeleC::do_diffuse   = false;

#ifdef USE_MASA
bool         PeleC::mms_initialized = false;
#endif

int           PeleC::les_model = 0;
int           PeleC::les_filter_type = no_filter;
int           PeleC::les_filter_fgr = 1;
int           PeleC::les_test_filter_type = box_3pt_optimized_approx;
int           PeleC::les_test_filter_fgr = 2;
int           PeleC::comp_Cs2 = 0;
int           PeleC::comp_CI = PeleC::comp_Cs2 + 1;
int           PeleC::comp_PrT = PeleC::comp_CI + 1;
int           PeleC::nCompC = PeleC::comp_PrT + 1;

#ifdef PELE_USE_EB
bool         PeleC::eb_initialized      = false;
bool         PeleC::no_eb_in_domain     = true;
bool         PeleC::body_state_set      = false;
std::vector<Real> PeleC::body_state;
#endif
bool         PeleC::do_react_load_balance = false;
bool         PeleC::do_mol_load_balance = false;


#ifdef AMREX_PARTICLES
SprayParticleContainer* PeleC::SprayPC = nullptr;
#endif

std::string  PeleC::probin_file = "probin";
std::vector<std::string> PeleC::spec_names;

std::vector<int> PeleC::src_list;

#if BL_SPACEDIM == 1
IntVect      PeleC::hydro_tile_size(1024);
#elif BL_SPACEDIM == 2
//IntVect      PeleC::hydro_tile_size(1024,16);
IntVect      PeleC::hydro_tile_size(1024,1024);
#else
IntVect      PeleC::hydro_tile_size(1024,16,16);
#endif

// this will be reset upon restart
Real         PeleC::previousCPUTimeUsed = 0.0;

Real         PeleC::startCPUTime = 0.0;

int          PeleC::num_state_type = 0;

#ifdef PELE_USE_EB
static bool eb_initialized = false;

bool ebInitialized()
{
  return eb_initialized;
}

void ebInitialized(bool eb_init_val)
{
  eb_initialized = eb_init_val;
}
#endif

void
PeleC::variableCleanUp ()
{
#ifdef AMREX_PARTICLES
  delete SprayPC;
  SprayPC = nullptr;;
#endif

  desc_lst.clear();

  clear_method_params();

  close_transport();

#ifdef REACTIONS
  if (do_react == 1)
  {
    close_reactor();
  }
#endif

  close_network();

  clear_grid_info();

  clear_prob();

#ifdef PELE_USE_EB
  eb_initialized = false;
#endif
}

void
PeleC::read_params ()
{
  static bool read_params_done = false;

  if (read_params_done) return;

  read_params_done = true;

  ParmParse pp("pelec");

#include <pelec_queries.H>
  
  pp.query("v",verbose);
  pp.query("sum_interval",sum_interval);
  pp.query("dump_old",dump_old);

  // Get boundary conditions
  Vector<string> lo_bc_char(BL_SPACEDIM);
  Vector<string> hi_bc_char(BL_SPACEDIM);
  pp.getarr("lo_bc",lo_bc_char,0,BL_SPACEDIM);
  pp.getarr("hi_bc",hi_bc_char,0,BL_SPACEDIM); 
  
  Vector<int> lo_bc(BL_SPACEDIM), hi_bc(BL_SPACEDIM);
  for (int dir = 0; dir<BL_SPACEDIM; dir++){
    if (!lo_bc_char[dir].compare("Interior")){
      lo_bc[dir] = 0;
    } else if (!lo_bc_char[dir].compare("Hard")){
      lo_bc[dir] = 1;
    } else if (!lo_bc_char[dir].compare("FOExtrap")){
      lo_bc[dir] = 2;
    } else if (!lo_bc_char[dir].compare("Symmetry")){
      lo_bc[dir] = 3;
    } else if (!lo_bc_char[dir].compare("SlipWall")){
      lo_bc[dir] = 4;
    } else if (!lo_bc_char[dir].compare("NoSlipWall")){
      lo_bc[dir] = 5;
    } else if (!lo_bc_char[dir].compare("UserBC")){
      lo_bc[dir] = 6;
    } else {
      amrex::Abort("Wrong boundary condition word in lo_bc, please use: Interior, UserBC, Symmetry, SlipWall, NoSlipWall");
    }
    
    if (!hi_bc_char[dir].compare("Interior")){
      hi_bc[dir] = 0;
    } else if (!hi_bc_char[dir].compare("Hard")){
      hi_bc[dir] = 1;
    } else if (!hi_bc_char[dir].compare("FOExtrap")){
      hi_bc[dir] = 2;
    } else if (!hi_bc_char[dir].compare("Symmetry")){
      hi_bc[dir] = 3;
    } else if (!hi_bc_char[dir].compare("SlipWall")){
      hi_bc[dir] = 4;
    } else if (!hi_bc_char[dir].compare("NoSlipWall")){
      hi_bc[dir] = 5;
    } else if (!hi_bc_char[dir].compare("UserBC")){
      hi_bc[dir] = 6;
    } else {
      amrex::Abort("Wrong boundary condition word in hi_bc, please use: Interior, UserBC, Symmetry, SlipWall, NoSlipWall");
    }
  }

  
  for (int i = 0; i < BL_SPACEDIM; i++)
  {
    phys_bc.setLo(i,lo_bc[i]);
    phys_bc.setHi(i,hi_bc[i]);
  }
    
  //
  // Check phys_bc against possible periodic geometry
  // if periodic, must have internal BC marked.
  //
  if (DefaultGeometry().isAnyPeriodic())
  {
    //
    // Do idiot check.  Periodic means interior in those directions.
    //
    for (int dir = 0; dir<BL_SPACEDIM; dir++)
    {
      if (DefaultGeometry().isPeriodic(dir))
      {
        if (lo_bc[dir] != Interior)
        {
          std::cerr << "PeleC::read_params:periodic in direction "
                    << dir
                    << " but low BC is not Interior\n";
          amrex::Error();
        }
        if (hi_bc[dir] != Interior)
        {
          std::cerr << "PeleC::read_params:periodic in direction "
                    << dir
                    << " but high BC is not Interior\n";
          amrex::Error();
        }
      }
    }
  }
  else
  {
    //
    // Do idiot check.  If not periodic, should be no interior.
    //
    for (int dir=0; dir<BL_SPACEDIM; dir++)
    {
      if (lo_bc[dir] == Interior)
      {
        std::cerr << "PeleC::read_params:interior bc in direction "
                  << dir
                  << " but not periodic\n";
        amrex::Error();
      }
      if (hi_bc[dir] == Interior)
      {
        std::cerr << "PeleC::read_params:interior bc in direction "
                  << dir
                  << " but not periodic\n";
        amrex::Error();
      }
    }
  }

  if ( DefaultGeometry().IsRZ() && (lo_bc[0] != Symmetry) ) {
    std::cerr << "ERROR:PeleC::read_params: must set r=0 boundary condition to Symmetry for r-z\n";
    amrex::Error();
  }

  // TODO: Any reason to support spherical in PeleC?
#if (BL_SPACEDIM == 1)
  if ( DefaultGeometry().IsSPHERICAL() )
  {
    if ( (lo_bc[0] != Symmetry) && (DefaultGeometry().ProbLo(0) == 0.0) )
    {
      std::cerr << "ERROR:PeleC::read_params: must set r=0 boundary condition to Symmetry for spherical\n";
      amrex::Error();
    }
  }
#elif (BL_SPACEDIM == 2)
  if ( DefaultGeometry().IsSPHERICAL() )
  {
    amrex::Abort("We don't support spherical coordinate systems in 2D");
  }
#elif (BL_SPACEDIM == 3)
  if ( DefaultGeometry().IsRZ() )
  {
    amrex::Abort("We don't support cylindrical coordinate systems in 3D");
  }
  else if ( DefaultGeometry().IsSPHERICAL() )
  {
    amrex::Abort("We don't support spherical coordinate systems in 3D");
  }
#endif

  pp.query("diffuse_temp",diffuse_temp);
  pp.query("diffuse_enth",diffuse_enth);
  pp.query("diffuse_spec",diffuse_spec);
  pp.query("diffuse_vel",diffuse_vel);
  pp.query("diffuse_cutoff_density",diffuse_cutoff_density);

  do_diffuse = diffuse_temp || diffuse_enth || diffuse_spec || diffuse_vel;
  
  // sanity checks
  if (cfl <= 0.0 || cfl > 1.0) {
    amrex::Error("Invalid CFL factor; must be between zero and one.");
  }

#ifdef PELEC_USE_MOL
  if (!do_mol_AD)
  {
    amrex::Abort("Must do_mol_AD=1 when compiled with HYP_TYPE = MOL\n");
  }
#else
  if (do_mol_AD)
  {
    amrex::Abort("Cannot do_mol_AD when compiled with HYP_TYPE != MOL\n");
  }
#endif

  if (do_les){
    pp.query("les_model",les_model);
    pp.query("les_test_filter_type",les_test_filter_type);
    pp.query("les_test_filter_fgr",les_test_filter_fgr);
  }

  if (use_explicit_filter){
    pp.query("les_filter_type",les_filter_type);
    pp.query("les_filter_fgr",les_filter_fgr);
  }

  // for the moment, ppm_type = 0 does not support ppm_trace_sources --
  // we need to add the momentum sources to the states (and not
  // add it in trans_3d
  if (ppm_type == 0 && ppm_trace_sources == 1)
  {
    amrex::Print() << "WARNING: ppm_trace_sources = 1 not implemented for ppm_type = 0" << std::endl;
    ppm_trace_sources = 0;
    pp.add("ppm_trace_sources",ppm_trace_sources);
  }

  if (ppm_temp_fix > 0 && BL_SPACEDIM == 1)
  {
    std::cerr << "ppm_temp_fix > 0 not implemented in 1-d \n";
    amrex::Error();
  }

  if (hybrid_riemann == 1 && BL_SPACEDIM == 1)
  {
    std::cerr << "hybrid_riemann only implemented in 2- and 3-d\n";
    amrex::Error();
  }

  if (hybrid_riemann == 1 && (DefaultGeometry().IsSPHERICAL() || DefaultGeometry().IsRZ() ))
  {
    std::cerr << "hybrid_riemann should only be used for Cartesian coordinates\n";
    amrex::Error();
  }

  if (use_colglaz >= 0)
  {
    std::cerr << "ERROR:: use_colglaz is deprecated.  Use riemann_solver instead\n";
    amrex::Error();
  }

  if (max_dt < fixed_dt)
  {
    std::cerr << "cannot have max_dt < fixed_dt\n";
    amrex::Error();
  }

#ifdef AMREX_PARTICLES
  read_particle_params();
#endif

  // TODO: What is this?
  StateDescriptor::setBndryFuncThreadSafety(bndry_func_thread_safe);

  // Get some useful amr inputs
  ParmParse ppa("amr");
  ppa.query("probin_file",probin_file);

  Vector<int> tilesize(BL_SPACEDIM);
  if (ppa.queryarr("hydro_tile_size", tilesize, 0, BL_SPACEDIM))
  {
    for (int i=0; i<BL_SPACEDIM; i++) hydro_tile_size[i] = tilesize[i];
  }

  // This turns on the lb stuff inside Amr, but we use our own flag to signal whether to gather data
  ppa.query("loadbalance_with_workestimates",do_mol_load_balance);
  ppa.query("loadbalance_with_workestimates",do_react_load_balance);
}

PeleC::PeleC ()
  : old_sources(num_src),
    new_sources(num_src)
#ifdef USE_MASA
  ,mms_src_evaluated(false)
#endif
{
}

PeleC::PeleC (Amr&            papa,
	      int             lev,
	      const Geometry& level_geom,
	      const BoxArray& bl,
	      const DistributionMapping& dm,
	      Real            time)
  : AmrLevel(papa,lev,level_geom,bl,dm,time),
    old_sources(num_src),
    new_sources(num_src)
#ifdef USE_MASA
  ,mms_src_evaluated(false)
#endif
{
  buildMetrics();
    
#ifdef PELE_USE_EB
  init_eb(level_geom, bl, dm);
  // #define PELE_UNIT_TEST_DN
#ifdef PELE_UNIT_TEST_DN
  test_dn();
#endif
#endif

  MultiFab& S_new = get_new_data(State_Type);

  for (int n = 0; n < src_list.size(); ++n)
  {
    old_sources[src_list[n]] = std::unique_ptr<MultiFab>(new MultiFab(grids, dmap, NUM_STATE, NUM_GROW,
                                                                      MFInfo(), Factory()));
    new_sources[src_list[n]] = std::unique_ptr<MultiFab>(new MultiFab(grids, dmap, NUM_STATE, S_new.nGrow(),
                                                                      MFInfo(), Factory()));
  }

#ifdef AMREX_PARTICLES
    Sborder.define(grids,dmap,NUM_STATE,NUM_GROW,MFInfo(),Factory());
#endif

  if (do_hydro)
  {
    Sborder.define(grids,dmap,NUM_STATE,NUM_GROW,MFInfo(),Factory());
  }
  else if (do_diffuse)
  {
    Sborder.define(grids,dmap,NUM_STATE,nGrowTr,MFInfo(),Factory());
  }

  if (!do_mol_AD)
  {
    if (do_hydro)
    {
      hydro_source.define(grids,dmap, NUM_STATE,NUM_GROW,MFInfo(),Factory());

      // This array holds the sum of all source terms that affect the hydrodynamics.
      // If we are doing the source term predictor, we'll also use this after the
      // hydro update to store the sum of the new-time sources, so that we can
      // compute the time derivative of the source terms.
      sources_for_hydro.define(grids,dmap,NUM_STATE,NUM_GROW,MFInfo(),Factory());
    }
  }
  else
  {
      Sborder.define(grids,dmap,NUM_STATE,nGrowTr,MFInfo(),Factory());
  }

  // Is this relevant for PeleC?
  for (int i = 0; i < n_lost; i++) {
    material_lost_through_boundary_cumulative[i] = 0.0;
    material_lost_through_boundary_temp[i] = 0.;
  }

  if (do_reflux && level > 0) {
    flux_reg.define(bl, papa.boxArray(level-1),
		    dm, papa.DistributionMap(level-1),
		    level_geom, papa.Geom(level-1),
		    papa.refRatio(level-1), level, NUM_STATE);
    
    if (!DefaultGeometry().IsCartesian())
    {
      pres_reg.define(bl, papa.boxArray(level-1),
		      dm, papa.DistributionMap(level-1),
		      level_geom, papa.Geom(level-1),
		      papa.refRatio(level-1), level, 1);
    }
  }
  
#ifdef REACTIONS
  get_new_data(Reactions_Type).setVal(0.0);
#endif

  // initialize the Godunov state array used in hydro -- we wait
  // until here so that ngroups is defined (if needed) in
  // rad_params_module
  if (do_hydro)
  {
    init_godunov_indices();
  }

  // initialize LES variables
  if (do_les){
    init_les();
  }
  
  // initialize filters and variables
  nGrowF = 0;
  if (use_explicit_filter){
    init_filters();
  }

}

PeleC::~PeleC ()
{
}


void
PeleC::buildMetrics ()
{
  const int ngrd = grids.size();

  radius.resize(ngrd);

  const Real* dx = geom.CellSize();

  for (int i = 0; i < ngrd; i++)
  {
    const Box& b = grids[i];
    int ilo      = b.smallEnd(0)-radius_grow;
    int ihi      = b.bigEnd(0)+radius_grow;
    int len      = ihi - ilo + 1;

    radius[i].resize(len);

    Real* rad = radius[i].dataPtr();

    if (DefaultGeometry().IsCartesian())
    {
      for (int j = 0; j < len; j++)
      {
        rad[j] = 1.0;
      }
    }
    else
    {
      RealBox gridloc = RealBox(grids[i],geom.CellSize(),geom.ProbLo());

      const Real xlo = gridloc.lo(0) + (0.5 - radius_grow)*dx[0];

      for (int j = 0; j < len; j++)
      {
        rad[j] = xlo + j*dx[0];
      }
    }
  }

  volume.clear();
  volume.define(grids,dmap,1,NUM_GROW,MFInfo(),FArrayBoxFactory());
  geom.GetVolume(volume);

  for (int dir = 0; dir < BL_SPACEDIM; dir++)
  {
    area[dir].clear();
    area[dir].define(getEdgeBoxArray(dir),dmap,1,NUM_GROW,MFInfo(),FArrayBoxFactory());
    geom.GetFaceArea(area[dir],dir);
  }

  dLogArea[0].clear();
#if (BL_SPACEDIM <= 2)
  geom.GetDLogA(dLogArea[0],grids,dmap,0,NUM_GROW);
#endif

#ifdef AMREX_USE_EB
  vfrac.clear();
  vfrac.define(grids,dmap,1,NUM_GROW,MFInfo(),Factory());
  const auto& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(Factory());
  MultiFab::Copy(vfrac, ebfactory.getVolFrac(), 0, 0, 1, NUM_GROW);
  areafrac = ebfactory.getAreaFrac();
#endif    

  level_mask.clear();
  level_mask.define(grids,dmap,1,1);
  level_mask.BuildMask(geom.Domain(), geom.periodicity(), 
          levmsk_covered,
          levmsk_notcovered,
          levmsk_physbnd,
          levmsk_interior);

  if (level == 0) setGridInfo();
}

    void
PeleC::setTimeLevel (Real time,
        Real dt_old,
        Real dt_new)
{
    AmrLevel::setTimeLevel(time,dt_old,dt_new);
}

    void
PeleC::setGridInfo ()
{
    /** Send refinement data to Fortran. We do it here
      because now the grids have been initialized and
      we need this data for setting up the problem.
      Note that this routine will always get called
      on level 0, even if we are doing a restart,
      so it is safe to put this here.
      */

    if (level == 0)
    {
        int max_level = parent->maxLevel();
        int nlevs = max_level + 1;

        Real dx_level[3*nlevs];
        int domlo_level[3*nlevs];
        int domhi_level[3*nlevs];

        const Real* dx_coarse = geom.CellSize();

        const int* domlo_coarse = geom.Domain().loVect();
        const int* domhi_coarse = geom.Domain().hiVect();

        for (int dir = 0; dir < 3; dir++)
        {
            dx_level[dir] = (ZFILL(dx_coarse))[dir];

            domlo_level[dir] = (ARLIM_3D(domlo_coarse))[dir];
            domhi_level[dir] = (ARLIM_3D(domhi_coarse))[dir];
        }

        for (int lev = 1; lev <= max_level; lev++)
        {
            IntVect ref_ratio = parent->refRatio(lev-1);

            // Note that we are explicitly calculating here what the
            // data would be on refined levels rather than getting the
            // data directly from those levels, because some potential
            // refined levels may not exist at the beginning of the simulation.

            for (int dir = 0; dir < 3; dir++)
            {
                if (dir < BL_SPACEDIM) {
                    dx_level[3 * lev + dir] = dx_level[3 * (lev - 1) + dir] / ref_ratio[dir];
                    int ncell = (domhi_level[3 * (lev - 1) + dir] - domlo_level[3 * (lev - 1) + dir] + 1) * ref_ratio[dir];
                    domlo_level[3 * lev + dir] = domlo_level[dir];
                    domhi_level[3 * lev + dir] = domlo_level[3 * lev + dir] + ncell - 1;
                } else {
                    dx_level[3 * lev + dir] = 0.0;
                    domlo_level[3 * lev + dir] = 0;
                    domhi_level[3 * lev + dir] = 0;
                }
            }
        }

        set_grid_info(max_level, dx_level, domlo_level, domhi_level);
    }
}

    void
PeleC::initData ()
{
    BL_PROFILE("PeleC::initData()");

    int ns          = NUM_STATE;
    const Real* dx  = geom.CellSize();
    MultiFab& S_new = get_new_data(State_Type);
    Real cur_time   = state[State_Type].curTime();

    S_new.setVal(0.);

    // make sure dx = dy = dz -- that's all we guarantee to support
#if (BL_SPACEDIM == 2)
    const Real SMALL = 1.e-13;
    if (fabs(dx[0] - dx[1]) > SMALL*dx[0])
    {
        amrex::Abort("dx != dy not supported");
    }
#elif (BL_SPACEDIM == 3)
    const Real SMALL = 1.e-13;
  if ( (fabs(dx[0] - dx[1]) > SMALL*dx[0]) || (fabs(dx[0] - dx[2]) > SMALL*dx[0]) )
  {
    amrex::Abort("dx != dy != dz not supported");
  }
#endif

  set_amr_info(level, -1, -1, -1.0, -1.0);

  if (verbose)
  {
    amrex::Print() << "Initializing the data at level " << level << std::endl;
  }

#ifdef REACTIONS
  get_new_data(Reactions_Type).setVal(0.);
#endif

  if (do_mol_load_balance || do_react_load_balance)
  {
    get_new_data(Work_Estimate_Type).setVal(1.0);
  }

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
  {
    RealBox gridloc = RealBox(grids[mfi.index()],geom.CellSize(),geom.ProbLo());
    const Box& box     = mfi.validbox();
    const int* lo      = box.loVect();
    const int* hi      = box.hiVect();

    pc_initdata(level, cur_time, ARLIM_3D(lo), ARLIM_3D(hi), ns,
                BL_TO_FORTRAN_3D(S_new[mfi]), ZFILL(dx),
                ZFILL(gridloc.lo()), ZFILL(gridloc.hi()));
    
    // Verify that the sum of (rho Y)_i = rho at every cell
    pc_check_initial_species(ARLIM_3D(lo), ARLIM_3D(hi), BL_TO_FORTRAN_3D(S_new[mfi]));
  }

  enforce_consistent_e(S_new);

  //computeTemp(S_new,0);

#ifdef PELE_USE_EB
  set_body_state(S_new);
#endif

  set_special_tagging_flag(cur_time);

#ifdef AMREX_PARTICLES
  if (level == 0) {
    init_particles();
  } else {
    particle_redistribute(level-1, true);
  }
#endif

  if (verbose)
  {
    amrex::Print() << "Done initializing level " << level << " data " << std::endl;
  }
}


void
PeleC::init (AmrLevel &old)
{
  BL_PROFILE("PeleC::init(old)");

  PeleC* oldlev = (PeleC*) &old;

  //
  // Create new grid data by fillpatching from old.
  //
  Real dt_new    = parent->dtLevel(level);
  Real cur_time  = oldlev->state[State_Type].curTime();
  Real prev_time = oldlev->state[State_Type].prevTime();
  Real dt_old    = cur_time - prev_time;
  setTimeLevel(cur_time,dt_old,dt_new);

  MultiFab& S_new = get_new_data(State_Type);
  FillPatch(old,S_new,0,cur_time,State_Type,0,NUM_STATE);

#ifdef REACTIONS
  MultiFab& React_new = get_new_data(Reactions_Type);

  if (do_react)
  {
    FillPatch(old,React_new,0,cur_time,Reactions_Type,0,React_new.nComp());
  }
  else
  {
    React_new.setVal(0);
  }
#endif

  if (do_mol_load_balance || do_react_load_balance)
  {
    MultiFab& work_estimate_new = get_new_data(Work_Estimate_Type);
    FillPatch(old,work_estimate_new,0,cur_time,Work_Estimate_Type,
              0,work_estimate_new.nComp());
  }
}

void
PeleC::init ()
{
  /**
     This version inits the data on a new level that did not
     exist before regridding.
  */
  BL_PROFILE("PeleC::init()");

  Real dt        = parent->dtLevel(level);
  Real cur_time  = getLevel(level-1).state[State_Type].curTime();
  Real prev_time = getLevel(level-1).state[State_Type].prevTime();

  Real dt_old = (cur_time - prev_time)/(Real)parent->MaxRefRatio(level-1);

  setTimeLevel(cur_time,dt_old,dt);
  MultiFab& S_new = get_new_data(State_Type);
  FillCoarsePatch(S_new, 0, cur_time, State_Type, 0, NUM_STATE);

  if (do_mol_load_balance || do_react_load_balance)
  {
    MultiFab& work_estimate_new = get_new_data(Work_Estimate_Type);
    int ncomp = work_estimate_new.nComp();
    FillCoarsePatch(work_estimate_new, 0, cur_time, Work_Estimate_Type, 0, ncomp);
  }
}

Real
PeleC::initialTimeStep ()
{
  BL_PROFILE("PeleC::initialTimeStep()");

  Real dummy_dt = 0.0;
  Real init_dt  = 0.0;

  if (initial_dt > 0.0)
  {
    init_dt = initial_dt;
  }
  else
  {
    init_dt = init_shrink*estTimeStep(dummy_dt);
  }

  return init_dt;
}

Real
PeleC::estTimeStep (Real dt_old)
{
  BL_PROFILE("PeleC::estTimeStep()");

  if (fixed_dt > 0.0) return fixed_dt;

  set_amr_info(level, -1, -1, -1.0, -1.0);

  Real estdt = max_dt;

  const MultiFab& stateMF = get_new_data(State_Type);

  const Real* dx = geom.CellSize();

  std::string limiter = "pelec.max_dt";

  // Start the hydro with the max_dt value, but divide by CFL
  // to account for the fact that we multiply by it at the end.
  // This ensures that if max_dt is more restrictive than the hydro
  // criterion, we will get exactly max_dt for a timestep.

  Real estdt_hydro = max_dt / cfl;

  if (do_hydro || do_mol_AD || diffuse_vel || diffuse_temp || diffuse_enth)
  {


#ifdef PELE_USE_EB
    auto const& fact = dynamic_cast<EBFArrayBoxFactory const&>(stateMF.Factory());
    auto const& flags = fact.getMultiEBCellFlagFab();
#endif

#ifdef _OPENMP
#pragma omp parallel reduction(min:estdt_hydro)
#endif
    {
      for (MFIter mfi(stateMF,true); mfi.isValid(); ++mfi)
      {
        const Box& box = mfi.tilebox();

#ifdef PELE_USE_EB
        const auto& flag_fab = flags[mfi];
        FabType typ = flag_fab.getType(box);
        if (typ == FabType::covered) {
          continue;
        }
#endif

        const auto& Sfab = stateMF[mfi];

        if (do_hydro)
        {
          Real dt = max_dt / cfl;
          pc_estdt(ARLIM_3D(box.loVect()), ARLIM_3D(box.hiVect()),
                   BL_TO_FORTRAN_3D(Sfab),
                   ZFILL(dx),&dt);
          estdt_hydro = std::min(estdt_hydro,dt);
        }

        if (diffuse_vel)
        {
          Real dt = max_dt / cfl;
          pc_estdt_vel_diffusion(ARLIM_3D(box.loVect()), ARLIM_3D(box.hiVect()),
                                 BL_TO_FORTRAN_3D(Sfab),ZFILL(dx),&dt);
          estdt_hydro = std::min(estdt_hydro,dt);
        }

        if (diffuse_temp)
        {
          Real dt = max_dt / cfl;
          pc_estdt_temp_diffusion(ARLIM_3D(box.loVect()), ARLIM_3D(box.hiVect()),
                                  BL_TO_FORTRAN_3D(Sfab),ZFILL(dx),&dt);
          estdt_hydro = std::min(estdt_hydro,dt);
        }

        if (diffuse_enth)
        {
          Real dt = max_dt / cfl;
          pc_estdt_enth_diffusion(ARLIM_3D(box.loVect()), ARLIM_3D(box.hiVect()),
                                  BL_TO_FORTRAN_3D(Sfab),ZFILL(dx),&dt);
          estdt_hydro = std::min(estdt_hydro,dt);
        }
      }
    }

    ParallelDescriptor::ReduceRealMin(estdt_hydro);
    estdt_hydro *= cfl;

    if (verbose)
    {
      amrex::Print() << "...estimated hydro-limited timestep at level "
                     << level << ": " << estdt_hydro << std::endl;
    }

    // Determine if this is more restrictive than the maximum timestep limiting
    if (estdt_hydro < estdt) {
      limiter = "hydro";
      estdt = estdt_hydro;
    }
  }

#ifdef AMREX_PARTICLES
  Real estdt_particle = max_dt;
  if (do_spray_particles)
  {
    particle_est_time_step(estdt_particle);
    if (estdt_particle < estdt) {
      limiter = "particles";
      estdt = estdt_particle;
    }
  }
#endif

  if (verbose)
  {
    amrex::Print() << "PeleC::estTimeStep (" << limiter << "-limited) at level "
                   << level << ":  estdt = " << estdt << '\n';
  }

  return estdt;
}

void
PeleC::computeNewDt (int                   finest_level,
		     int                   sub_cycle,
		     Vector<int>&           n_cycle,
		     const Vector<IntVect>& ref_ratio,
		     Vector<Real>&          dt_min,
		     Vector<Real>&          dt_level,
		     Real                  stop_time,
		     int                   post_regrid_flag)
{
  BL_PROFILE("PeleC::computeNewDt()");

  //
  // We are at the start of a coarse grid timecycle.
  // Compute the timesteps for the next iteration.
  //
  if (level > 0) return;

  Real dt_0 = 1.0e+100;
  int n_factor = 1;
  for (int i = 0; i <= finest_level; i++)
  {
    PeleC& adv_level = getLevel(i);
    dt_min[i] = adv_level.estTimeStep(dt_level[i]);
  }

  if (fixed_dt <= 0.0)
  {
    if (post_regrid_flag == 1)
    {
      // Limit dt's by pre-regrid dt
      for (int i = 0; i <= finest_level; i++)
      {
        dt_min[i] = std::min(dt_min[i],dt_level[i]);
      }
    }
    else
    {
      // Limit dt's by change_max * old dt
      for (int i = 0; i <= finest_level; i++)
      {
        if (verbose && ParallelDescriptor::IOProcessor())
        {
          if (dt_min[i] > change_max*dt_level[i])
          {
            cout << "PeleC::compute_new_dt : limiting dt at level "
                 << i << '\n';
            cout << " ... new dt computed: " << dt_min[i]
                 << '\n';
            cout << " ... but limiting to: "
                 << change_max * dt_level[i] << " = " << change_max
                 << " * " << dt_level[i] << '\n';
          }
        }
        dt_min[i] = std::min(dt_min[i],change_max*dt_level[i]);
      }
    }
  }

  // Find the minimum over all levels
  for (int i = 0; i <= finest_level; i++)
  {
    n_factor *= n_cycle[i];
    dt_0 = std::min(dt_0,n_factor*dt_min[i]);
  }

  // Limit dt's by the value of stop_time.
  const Real eps = 0.001*dt_0;
  Real cur_time  = state[State_Type].curTime();
  if (stop_time >= 0.0)
  {
    if ((cur_time + dt_0) > (stop_time - eps))
    {
      dt_0 = stop_time - cur_time;
    }
  }

  n_factor = 1;
  for (int i = 0; i <= finest_level; i++)
  {
    n_factor *= n_cycle[i];
    dt_level[i] = dt_0/n_factor;
  }
}

void
PeleC::computeInitialDt (int                   finest_level,
			 int                   sub_cycle,
			 Vector<int>&           n_cycle,
			 const Vector<IntVect>& ref_ratio,
			 Vector<Real>&          dt_level,
			 Real                  stop_time)
{
  BL_PROFILE("PeleC::computeInitialDt()");
  
  // Grids have been constructed, compute dt for all levels.
  if (level > 0) return;

  Real dt_0 = 1.0e+100;
  int n_factor = 1;
  ///TODO/DEBUG: This will need to change for optimal subcycling.
  for (int i = 0; i <= finest_level; i++)
  {
    dt_level[i] = getLevel(i).initialTimeStep();
    n_factor   *= n_cycle[i];
    dt_0 = std::min(dt_0,n_factor*dt_level[i]);
  }

  // Limit dt's by the value of stop_time.
  const Real eps = 0.001*dt_0;
  Real cur_time  = state[State_Type].curTime();
  if (stop_time >= 0.0)
  {
    if ((cur_time + dt_0) > (stop_time - eps))
    {
      dt_0 = stop_time - cur_time;
    }
  }

  n_factor = 1;
  for (int i = 0; i <= finest_level; i++)
  {
    n_factor *= n_cycle[i];
    dt_level[i] = dt_0/n_factor;
  }
}

void
PeleC::post_timestep (int iteration)
{
  BL_PROFILE("PeleC::post_timestep()");

  const int finest_level = parent->finestLevel();
  const int ncycle       = parent->nCycle(level);

#ifdef AMREX_PARTICLES
  if (do_spray_particles)
  {
    //
    // Remove virtual particles at this level if we have any.
    //
    remove_virtual_particles();

    //
    // Remove Ghost particles on the final iteration
    //
    if (iteration == parent->nCycle(level))
      remove_ghost_particles();

    //
    // Sync up if we're level 0 or if we have particles that may have moved
    // off the next finest level and need to be added to our own level.
    //
    if ((iteration < ncycle and level < finest_level) || level == 0)
    {
      theSprayPC()->Redistribute(level,theSprayPC()->finestLevel(),iteration);
    }
  }
#endif

  if (do_reflux && level < finest_level)
  {
    reflux();

    // We need to do this before anything else because refluxing changes the values of coarse cells
    //    underneath fine grids with the assumption they'll be over-written by averaging down
    if (level < finest_level)
    {
      avgDown();
    }

    // Clean up any aberrant state data generated by the reflux.
    MultiFab& S_new_crse = get_new_data(State_Type);
    // clean_state(S_new_crse);

    // Flush Fortran output
    if (verbose)
    {
      flush_output();
    }
  }

  // Re-compute temperature after all the other updates.
  MultiFab& S_new = get_new_data(State_Type);
  int ng_pts = 0;
  computeTemp(S_new, ng_pts);

#ifdef DO_PROBLEM_POST_TIMESTEP

  problem_post_timestep();

#endif

  if (level == 0)
  {
    int nstep = parent->levelSteps(0);
    Real dtlev = parent->dtLevel(0);
    Real cumtime = parent->cumTime() + dtlev;

    bool sum_int_test = (sum_interval > 0  &&  nstep%sum_interval == 0);

    bool sum_per_test = false;

    if (sum_per > 0.0)
    {
      const int num_per_old = floor((cumtime - dtlev) / sum_per);
      const int num_per_new = floor((cumtime        ) / sum_per);

      if (num_per_old != num_per_new)
      {
        sum_per_test = true;
      }
    }

    if (sum_int_test || sum_per_test)
    {
      sum_integrated_quantities();
    }
  }
}

void
PeleC::post_restart ()
{
  BL_PROFILE("PeleC::post_restart()");

  Real cur_time = state[State_Type].curTime();

#ifdef AMREX_PARTICLES
  if (do_spray_particles)
  {
    particle_post_restart(parent->theRestartFile());
  }
#endif

  set_special_tagging_flag(cur_time);

  // initialize the Godunov state array used in hydro -- we wait
  // until here so that ngroups is defined (if needed) in
  // rad_params_module
  if (do_hydro)
  {
    init_godunov_indices();
  }

  // initialize LES variables
  if (do_les){
    init_les();
  }

  // initialize filters and variables
  nGrowF = 0;
  if (use_explicit_filter){
    init_filters();
  }

#ifdef DO_PROBLEM_POST_RESTART
  problem_post_restart();
#endif
}

void
PeleC::postCoarseTimeStep (Real cumtime)
{
  BL_PROFILE("PeleC::postCoarseTimeStep()");
  AmrLevel::postCoarseTimeStep(cumtime);
}

void
PeleC::post_regrid (int lbase,
		    int new_finest)
{
  BL_PROFILE("PeleC::post_regrid()");
  fine_mask.clear();

#ifdef AMREX_PARTICLES
  if (do_spray_particles && SprayPC && level == lbase)
  {
    SprayPC->Redistribute(false, false, lbase);
  }
#endif
}

void
PeleC::post_init (Real stop_time)
{
  BL_PROFILE("PeleC::post_init()");

  Real dtlev = parent->dtLevel(level);
  Real cumtime = parent->cumTime();

  // Fill Reactions_Type data based on initial dt
#ifdef REACTIONS
  if (do_react == 1)
  {
    bool react_init = true;
    react_state(cumtime,dtlev,react_init);
  }
#endif

  if (level > 0)
    return;

  //
  // Average data down from finer levels
  // so that conserved data is consistent between levels.
  //
  if (do_avg_down != 0) {
    int finest_level = parent->finestLevel();
    for (int k = finest_level-1; k>= 0; k--)
    {
      getLevel(k).avgDown();
    }
  }

#ifdef DO_PROBLEM_POST_INIT
  //
  // Allow the user to define their own post_init functions.
  //
  problem_post_init();
#endif

  int nstep = parent->levelSteps(0);
  if (cumtime != 0.0) cumtime += dtlev;

  bool sum_int_test = false;

  if (sum_interval > 0)
  {
    if (nstep%sum_interval == 0)
    {
      sum_int_test = true;
    }
  }

  bool sum_per_test = false;

  if (sum_per > 0.0)
  {
    const int num_per_old = floor((cumtime - dtlev) / sum_per);
    const int num_per_new = floor((cumtime        ) / sum_per);

    if (num_per_old != num_per_new)
    {
      sum_per_test = true;
    }
  }

  if (sum_int_test || sum_per_test)
  {
    sum_integrated_quantities();
  }
}

int
PeleC::okToContinue ()
{
  if (level > 0)
  {
    return 1;
  }

  int test = 1;

  if (signalStopJob)
  {
    test = 0;

    amrex::Print() << " Signalling a stop of the run due to signalStopJob = true." << std::endl;
  }
  else if (parent->dtLevel(0) < dt_cutoff)
  {
    test = 0;

    amrex::Print() << " Signalling a stop of the run because dt < dt_cutoff." << std::endl;
  }

  return test;
}

void
PeleC::reflux ()
{
  BL_PROFILE("PeleC::reflux()");

  BL_ASSERT(level<parent->finestLevel());

  const Real strt = ParallelDescriptor::second();

  PeleC& fine_level = getLevel(level+1);
  MultiFab& S_crse = get_new_data(State_Type);

#ifdef PELE_USE_EB

  MultiFab& S_fine = fine_level.get_new_data(State_Type);

  fine_level.flux_reg.Reflux(S_crse, vfrac, S_fine, fine_level.vfrac);

  if (!DefaultGeometry().IsCartesian())
  {
    amrex::Abort("rz not yet compatible with EB");
  }

#else

  fine_level.flux_reg.Reflux(S_crse);

  if (!DefaultGeometry().IsCartesian())
  {
    MultiFab dr(volume.boxArray(),volume.DistributionMap(),1,volume.nGrow(),
                MFInfo(),FArrayBoxFactory());
    dr.setVal(geom.CellSize(0));
    amrex::Abort("PeleC reflux not yet ready for r-z");
    //fine_level.pres_reg.Reflux(S_crse,dr,1.0,0,Xmom,1,geom);
  }
#endif

  if (verbose)
  {
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    Real      end    = ParallelDescriptor::second() - strt;

#ifdef BL_LAZY
    Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(end,IOProc);

        amrex::Print() << "PeleC::reflux() at level " << level << " : time = " << end << std::endl;
#ifdef BL_LAZY
      });
#endif
  }
}

void
PeleC::avgDown ()
{
  BL_PROFILE("PeleC::avgDown()");

  if (level == parent->finestLevel()) return;

  avgDown(State_Type);

#ifdef REACTIONS
  avgDown(Reactions_Type);
#endif

}

void
PeleC::normalize_species (MultiFab& S)
{
    amrex::Abort("We don't normalize species!");
}

void
PeleC::enforce_consistent_e (MultiFab& S)
{

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S,true); mfi.isValid(); ++mfi)
    {
        const Box& box     = mfi.tilebox();
        const int* lo      = box.loVect();
        const int* hi      = box.hiVect();

        pc_enforce_consistent_e(ARLIM_3D(lo), ARLIM_3D(hi), BL_TO_FORTRAN_3D(S[mfi]));
    }
}

Real
PeleC::enforce_min_density (MultiFab& S_old, MultiFab& S_new)
{

  /** This routine sets the density in S_new to be larger than the density floor.
      Note that it will operate everywhere on S_new, including ghost zones.
      S_old is present so that, after the hydro call, we know what the old density
      was so that we have a reference for comparison. If you are calling it elsewhere
      and there's no meaningful reference state, just pass in the same MultiFab twice.
      @return  The return value is the the negative fractional change in the state that has the
      largest magnitude. If there is no reference state, this is meaningless.
  */


  Real dens_change = 1.e0;

  Real mass_added = 0.0;
  Real eint_added = 0.0;
  Real eden_added = 0.0;

#ifdef PELE_USE_EB
  auto const& fact = dynamic_cast<EBFArrayBoxFactory const&>(S_new.Factory());
  auto const& flags = fact.getMultiEBCellFlagFab();
#endif

#ifdef _OPENMP
#pragma omp parallel reduction(min:dens_change)
#endif
  for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi) {

    const Box& bx = mfi.growntilebox();

#ifdef PELE_USE_EB
    const auto& flag_fab = flags[mfi];
    FabType typ = flag_fab.getType(bx);
    if (typ == FabType::covered)
    {
      continue;
    }
#endif

    const auto& stateold = S_old[mfi];
    auto& statenew = S_new[mfi];
    const auto& vol = volume[mfi];

    enforce_minimum_density(stateold.dataPtr(), ARLIM_3D(stateold.loVect()), ARLIM_3D(stateold.hiVect()),
                            statenew.dataPtr(), ARLIM_3D(statenew.loVect()), ARLIM_3D(statenew.hiVect()),
                            vol.dataPtr(), ARLIM_3D(vol.loVect()), ARLIM_3D(vol.hiVect()),
                            ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
                            &mass_added, &eint_added, &eden_added, &dens_change,
                            &verbose);

  }

  if (print_energy_diagnostics)
  {
    Real foo[3] = {mass_added, eint_added, eden_added};

#ifdef BL_LAZY
    Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealSum(foo, 3, ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
        {
          mass_added = foo[0];
          eint_added = foo[1];
          eden_added = foo[2];

          if (std::abs(mass_added) != 0.0)
          {
            std::cout << "   Mass added from negative density correction : " <<
              mass_added << std::endl;
            std::cout << "(rho e) added from negative density correction : " <<
              eint_added << std::endl;
            std::cout << "(rho E) added from negative density correction : " <<
              eden_added << std::endl;
          }
        }
#ifdef BL_LAZY
      });
#endif

  }

  return dens_change;
}

void
PeleC::avgDown (int state_indx)
{
  BL_PROFILE("PeleC::avgDown(state_indx)");

  if (level == parent->finestLevel()) return;

  MultiFab&  S_crse   = get_new_data(state_indx);
  MultiFab&  S_fine   = getLevel(level+1).get_new_data(state_indx);

#ifdef PELE_USE_EB

  PeleC& fine_lev = getLevel(level+1);

  amrex::EB_average_down(S_fine, S_crse, fine_lev.Volume(), fine_lev.volFrac(),
                         0, S_fine.nComp(), fine_ratio);

  if (state_indx == State_Type) {
    set_body_state(S_crse); // TODO: Is this necessary?
  }

#else

  const Geometry& fgeom = getLevel(level+1).geom;
  const Geometry& cgeom =                   geom;

  amrex::average_down(S_fine, S_crse,
		      fgeom, cgeom,
		      0, S_fine.nComp(), fine_ratio);

#endif
}

void
PeleC::allocOldData ()
{
  for (int k = 0; k < num_state_type; k++)
  {
    state[k].allocOldData();
  }
}

void
PeleC::removeOldData()
{
  AmrLevel::removeOldData();
}

void
PeleC::errorEst (TagBoxArray& tags,
		 int          clearval,
		 int          tagval,
		 Real         time,
		 int          n_error_buf,
		 int          ngrow)
{
  BL_PROFILE("PeleC::errorEst()");

  const Real cur_time = state[State_Type].curTime();
  MultiFab S_data(get_new_data(State_Type).boxArray(), get_new_data(State_Type).DistributionMap(), NUM_STATE, 1);
  FillPatch(*this, S_data, S_data.nGrow(), cur_time, State_Type, Density, NUM_STATE, 0);
 
  const int*  domlo = geom.Domain().loVect();
  const int*  domhi = geom.Domain().hiVect();
  const Real* dx        = geom.CellSize();
  Real dt = parent->dtLevel(level);
  const Real* prob_lo   = geom.ProbLo();

  Vector<BCRec>       bcs(NUM_STATE);
  
#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    FArrayBox S_derData;
    Vector<int>  itags;

    for (MFIter mfi(S_data,false); mfi.isValid(); ++mfi)
    {
      FArrayBox   &datfab = S_data[mfi];
      auto&       tagfab  = tags[mfi];
      const Box&  tilebx  = mfi.tilebox();
      const int&  grid_no  = mfi.index();
      
      const RealBox& pbx  = RealBox(tilebx,geom.CellSize(),geom.ProbLo());
      const Box&  datbox  = datfab.box();

      // We cannot pass tagfab to Fortran becuase it is BaseFab<char>.
      // So we are going to get a temporary integer array.
      tagfab.get_itags(itags, tilebx);

      int*        tptr    = itags.dataPtr();
      const int*  tlo     = tilebx.loVect();
      const int*  thi     = tilebx.hiVect();
      const int*  lo      = tlo;
      const int*  hi      = thi;
      const Real* xlo     = pbx.lo();
      const int*  dlo     = datbox.loVect();
      const int*  dhi     = datbox.hiVect();
      const int   ncomp   = datfab.nComp();

#ifdef PELEC_USE_EB
      FArrayBox   &vfracfab = vfrac[mfi];
#endif

      S_derData.resize(datbox, 1);
      const int   ncp   = S_derData.nComp();
      const int* bc =  bcs[0].data();
      
      // Tagging Density
      pc_denerror(tptr,ARLIM_3D(tlo), ARLIM_3D(thi),
                  &tagval, &clearval,
                  BL_TO_FORTRAN_3D(S_data[mfi]),
                  ARLIM_3D(lo),ARLIM_3D(hi), &ncomp, domlo,domhi, 
                  ZFILL(dx), ZFILL(xlo), ZFILL(prob_lo), &time, &level);

      //----------------------
      // Recasting pressure
      // Warning: bcs are dummy values, and level is passed as grid_no in the fortran routine
      //          they are not used, but one may want the correct values 
      S_derData.setVal(0.0);
      pc_derpres(S_derData.dataPtr(), ARLIM_3D(S_derData.loVect()), ARLIM_3D(S_derData.hiVect()),&ncp,
                 BL_TO_FORTRAN_3D(S_data[mfi]),&ncomp,
                 ARLIM_3D(dlo),ARLIM_3D(dhi),domlo,domhi,
                 ZFILL(dx), ZFILL(xlo),&time,&dt,bc,&level,&grid_no);
      
      // Tagging Pressure
      pc_presserror(tptr,ARLIM_3D(tlo), ARLIM_3D(thi),
                  &tagval, &clearval,
                  S_derData.dataPtr(), ARLIM_3D(S_derData.loVect()), ARLIM_3D(S_derData.hiVect()),
                  ARLIM_3D(lo),ARLIM_3D(hi), &ncomp, domlo,domhi, 
                  ZFILL(dx), ZFILL(xlo), ZFILL(prob_lo), &time, &level);

      //----------------------
      // Recasting vel_x
      S_derData.setVal(0.0);
      pc_dervelx(S_derData.dataPtr(), ARLIM_3D(S_derData.loVect()), ARLIM_3D(S_derData.hiVect()),&ncp,
                 BL_TO_FORTRAN_3D(S_data[mfi]),&ncomp,
                 ARLIM_3D(dlo),ARLIM_3D(dhi),domlo,domhi,
                 ZFILL(dx), ZFILL(xlo),&time,&dt,bc,&level,&grid_no);

      // Tagging vel_x
      pc_velerror(tptr,ARLIM_3D(tlo), ARLIM_3D(thi),
                  &tagval, &clearval,
                  S_derData.dataPtr(), ARLIM_3D(S_derData.loVect()), ARLIM_3D(S_derData.hiVect()),
                  ARLIM_3D(lo),ARLIM_3D(hi), &ncomp, domlo,domhi, 
                  ZFILL(dx), ZFILL(xlo), ZFILL(prob_lo), &time, &level);
      
#if (BL_SPACEDIM >= 2) 
      //----------------------
      // Recasting vel_y
      S_derData.setVal(0.0);
      pc_dervely(S_derData.dataPtr(), ARLIM_3D(S_derData.loVect()), ARLIM_3D(S_derData.hiVect()),&ncp,
                 BL_TO_FORTRAN_3D(S_data[mfi]),&ncomp,
                 ARLIM_3D(dlo),ARLIM_3D(dhi),domlo,domhi,
                 ZFILL(dx), ZFILL(xlo),&time,&dt,bc,&level,&grid_no);

      // Tagging vel_y
      pc_velerror(tptr,ARLIM_3D(tlo), ARLIM_3D(thi),
                  &tagval, &clearval,
                  S_derData.dataPtr(), ARLIM_3D(S_derData.loVect()), ARLIM_3D(S_derData.hiVect()),
                  ARLIM_3D(lo),ARLIM_3D(hi), &ncomp, domlo,domhi, 
                  ZFILL(dx), ZFILL(xlo), ZFILL(prob_lo), &time, &level);
#endif

#if (BL_SPACEDIM == 3)
      //----------------------
      // Recasting vel_z
      S_derData.setVal(0.0);
      pc_dervelz(S_derData.dataPtr(), ARLIM_3D(S_derData.loVect()), ARLIM_3D(S_derData.hiVect()),&ncp,
                 BL_TO_FORTRAN_3D(S_data[mfi]),&ncomp,
                 ARLIM_3D(dlo),ARLIM_3D(dhi),domlo,domhi,
                 ZFILL(dx), ZFILL(xlo),&time,&dt,bc,&level,&grid_no);

      // Tagging vel_z
      pc_velerror(tptr,ARLIM_3D(tlo), ARLIM_3D(thi),
                  &tagval, &clearval,
                  S_derData.dataPtr(), ARLIM_3D(S_derData.loVect()), ARLIM_3D(S_derData.hiVect()),
                  ARLIM_3D(lo),ARLIM_3D(hi), &ncomp, domlo,domhi, 
                  ZFILL(dx), ZFILL(xlo), ZFILL(prob_lo), &time, &level);
#endif

      //----------------------
      // Recasting magVort
      S_derData.setVal(0.0);
      pc_dermagvort(S_derData.dataPtr(), ARLIM_3D(S_derData.loVect()), ARLIM_3D(S_derData.hiVect()),&ncp,
                 BL_TO_FORTRAN_3D(S_data[mfi]),&ncomp,
                 ARLIM_3D(lo),ARLIM_3D(hi),domlo,domhi,
                 ZFILL(dx), ZFILL(xlo),&time,&dt,bc,&level,&grid_no);
      
      // Tagging magVorticity
      pc_vorterror(tptr,ARLIM_3D(tlo), ARLIM_3D(thi),
                  &tagval, &clearval,
                  S_derData.dataPtr(), ARLIM_3D(S_derData.loVect()), ARLIM_3D(S_derData.hiVect()),
                  ARLIM_3D(lo),ARLIM_3D(hi), &ncomp, domlo,domhi, 
                  ZFILL(dx), ZFILL(xlo), ZFILL(prob_lo), &time, &level);
      
      //----------------------
      // Recasting Temperature
      S_derData.setVal(0.0);
      pc_dertemp(S_derData.dataPtr(), ARLIM_3D(S_derData.loVect()), ARLIM_3D(S_derData.hiVect()),&ncp,
                 BL_TO_FORTRAN_3D(S_data[mfi]),&ncomp,
                 ARLIM_3D(dlo),ARLIM_3D(dhi),domlo,domhi,
                 ZFILL(dx), ZFILL(xlo),&time,&dt,bc,&level,&grid_no);
      
      // Tagging Temperature
      pc_temperror(tptr,ARLIM_3D(tlo), ARLIM_3D(thi),
                  &tagval, &clearval,
                  S_derData.dataPtr(), ARLIM_3D(S_derData.loVect()), ARLIM_3D(S_derData.hiVect()),
                  ARLIM_3D(lo),ARLIM_3D(hi), &ncomp, domlo,domhi, 
                  ZFILL(dx), ZFILL(xlo), ZFILL(prob_lo), &time, &level);
      
      //----------------------
      // Recasting Flame Tracer
      if (!flame_trac_name.empty())
      {
        int idx = -1;
        for (int i=0; i<spec_names.size(); ++i)
        {
          if (flame_trac_name == spec_names[i])
          {
            idx = i;
          }
        }
    
        if (idx >= 0)
        {
          //const std::string name = "Y("+flame_trac_name+")";
          //if (ParallelDescriptor::IOProcessor())
          //  std::cout << " Flame tracer will be " << name << '\n';
        
        S_derData.setVal(0.0);
        pc_derspectrac(S_derData.dataPtr(), ARLIM_3D(S_derData.loVect()), ARLIM_3D(S_derData.hiVect()),&ncp,
                 BL_TO_FORTRAN_3D(S_data[mfi]),&ncomp,
                 ARLIM_3D(dlo),ARLIM_3D(dhi),domlo,domhi,
                 ZFILL(dx), ZFILL(xlo),&time,&dt,bc,&level,&grid_no,&idx);
        
        // Tagging Flame Tracer
        pc_ftracerror(tptr,ARLIM_3D(tlo), ARLIM_3D(thi),
                  &tagval, &clearval,
                  S_derData.dataPtr(), ARLIM_3D(S_derData.loVect()), ARLIM_3D(S_derData.hiVect()),
                  ARLIM_3D(lo),ARLIM_3D(hi), &ncomp, domlo,domhi, 
                  ZFILL(dx), ZFILL(xlo), ZFILL(prob_lo), &time, &level);
      
        }
        else
        {
          amrex::Abort("Unknown species identified as flame_trac_name");
        }
      } 
          
      //----------------------

#ifdef PELEC_USE_EB
      pc_vfracerror(tptr,ARLIM_3D(tlo), ARLIM_3D(thi),
                  &tagval, &clearval,
                  vfracfab.dataPtr(), ARLIM_3D(vfracfab.loVect()), ARLIM_3D(vfracfab.hiVect()),
                  ARLIM_3D(lo),ARLIM_3D(hi), &ncomp, domlo,domhi, 
                  ZFILL(dx), ZFILL(xlo), ZFILL(prob_lo), &time, &level);
#endif

      //----------------------
      // Problem specific tagging
      set_problem_tags(tptr,ARLIM_3D(tlo), ARLIM_3D(thi),
                       BL_TO_FORTRAN_3D(S_data[mfi]),
                       &tagval, &clearval,
                       ARLIM_3D(lo),ARLIM_3D(hi),
                       ZFILL(dx), ZFILL(prob_lo), &time, &level);

      // Now update the tags in the TagBox.
      tagfab.tags(itags, tilebx);
      
    }
  } 
}

std::unique_ptr<MultiFab>
PeleC::derive (const std::string& name,
	       Real           time,
	       int            ngrow)
{

  if ((do_les) && (name == "C_s2")) {
    std::unique_ptr<MultiFab> derive_dat(new MultiFab(grids,dmap,1,0));
    MultiFab::Copy(*derive_dat,LES_Coeffs,comp_Cs2,0,1,0);
    return derive_dat;
  }
  else if ((do_les) && (name == "C_I")) {
    std::unique_ptr<MultiFab> derive_dat(new MultiFab(grids,dmap,1,0));
    MultiFab::Copy(*derive_dat,LES_Coeffs,comp_CI,0,1,0);
    return derive_dat;
  }
  else if ((do_les) && (name == "Pr_T")) {
    std::unique_ptr<MultiFab> derive_dat(new MultiFab(grids,dmap,1,0));
    MultiFab::Copy(*derive_dat,LES_Coeffs,comp_PrT,0,1,0);
    return derive_dat;
  }

#ifdef PELE_USE_EB
  if (name == "vfrac") {
    std::unique_ptr<MultiFab> mf(new MultiFab(grids, dmap, 1, ngrow, MFInfo()));
    if (ngrow > 0) {
      mf->setBndry(0);
    }
    MultiFab::Copy(*mf,vfrac,0,0,1,0);
    return mf;
  }
#endif

#ifdef AMREX_PARTICLES
  return particle_derive(name,time,ngrow);
#else
  return AmrLevel::derive(name,time,ngrow);
#endif
}

void
PeleC::derive (const std::string& name,
	       Real           time,
	       MultiFab&      mf,
	       int            dcomp)
{
#ifdef PELE_USE_EB
  if (name == "vfrac") {
    MultiFab::Copy(mf,vfrac,0,dcomp,1,0);
  } else
#endif
    if (1) {
      AmrLevel::derive(name,time,mf,dcomp);
    }
}

void
PeleC::init_network ()
{
  pc_network_init();
}

void
PeleC::close_network ()
{
  pc_network_close();
}

void
PeleC::clear_prob ()
{
  pc_prob_close();
}

void
PeleC::init_extern ()
{
  // initialize the external runtime parameters -- these will
  // live in the probin

  amrex::Print() << "reading extern runtime parameters ..." << std::endl;

  int probin_file_length = probin_file.length();
  Vector<int> probin_file_name(probin_file_length);

  for (int i = 0; i < probin_file_length; i++)
  {
    probin_file_name[i] = probin_file[i];
  }

  pc_extern_init(probin_file_name.dataPtr(),&probin_file_length);
}

#ifdef REACTIONS
void
PeleC::init_reactor ()
{
  pc_reactor_init();
}

void
PeleC::close_reactor ()
{
  pc_reactor_close();
}
#endif

void
PeleC::init_transport ()
{
  pc_transport_init();
}

void
PeleC::close_transport ()
{
  pc_transport_close();
}

void
PeleC::init_les ()
{
    // Fill with default coefficient values
    LES_Coeffs.define(grids, dmap, nCompC, 1);
    LES_Coeffs.setVal(0.0);
    LES_Coeffs.setVal(Cs*Cs,comp_Cs2,1,LES_Coeffs.nGrow());
    LES_Coeffs.setVal(CI,comp_CI,1,LES_Coeffs.nGrow());
    LES_Coeffs.setVal(PrT,comp_PrT,1,LES_Coeffs.nGrow());
}

void
PeleC::init_filters ()
{
  if (level > 0){
    IntVect ref_ratio = parent->refRatio(level-1);
    les_filter = Filter(les_filter_type, les_filter_fgr * std::pow(ref_ratio[0], level));
  }
  else{
    les_filter = Filter(les_filter_type, les_filter_fgr);
  }

  nGrowF = les_filter.get_filter_ngrow();

  // Add grow cells necessary for explicit filtering of source terms
  if (do_hydro)
  {
    Sborder.define(grids,dmap,NUM_STATE,NUM_GROW+nGrowF,MFInfo(),Factory());
    hydro_source.define(grids,dmap, NUM_STATE, hydro_source.nGrow()+nGrowF, MFInfo(),Factory());
    sources_for_hydro.define(grids,dmap,NUM_STATE,sources_for_hydro.nGrow()+nGrowF,MFInfo(),Factory());
  }

  volume.clear();
  volume.define(grids,dmap,1,NUM_GROW+nGrowF,MFInfo(),FArrayBoxFactory());
  geom.GetVolume(volume);

  for (int dir = 0; dir < BL_SPACEDIM; dir++)
  {
    area[dir].clear();
    area[dir].define(getEdgeBoxArray(dir),dmap,1,NUM_GROW+nGrowF,MFInfo(),FArrayBoxFactory());
    geom.GetFaceArea(area[dir],dir);
  }

  dLogArea[0].clear();
#if (BL_SPACEDIM <= 2)
  geom.GetDLogA(dLogArea[0],grids,dmap,0,NUM_GROW+nGrowF);
#endif
}

#ifdef USE_MASA
void
PeleC::init_mms()
{
  if (!mms_initialized)
  {
    if (verbose && ParallelDescriptor::IOProcessor())
    {
      std::cout << "Initializing MMS" << std::endl;
    }
#ifdef USE_MASA
    masa_init("mms", masa_solution_name.c_str());
#endif
    mms_initialized = true;
  }
}
#endif

void
PeleC::reset_internal_energy(MultiFab& S_new, int ng)
{
    Real sum  = 0.;
    Real sum0 = 0.;

    if (parent->finestLevel() == 0 && print_energy_diagnostics)
    {
        // Pass in the multifab and the component
        sum0 = volWgtSumMF(&S_new,Eden,true);
    }

    // Ensure (rho e) isn't too small or negative
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(ng);

        reset_internal_e(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
                         BL_TO_FORTRAN_3D(S_new[mfi]),
			 print_fortran_warnings);
    }

    // Flush Fortran output

    if (verbose)
	flush_output();

    if (parent->finestLevel() == 0 && print_energy_diagnostics)
    {
        // Pass in the multifab and the component
        sum = volWgtSumMF(&S_new,Eden,true);
#ifdef BL_LAZY
	Lazy::QueueReduction( [=] () mutable {
#endif
		ParallelDescriptor::ReduceRealSum(sum0);
		ParallelDescriptor::ReduceRealSum(sum);
		if (ParallelDescriptor::IOProcessor() && std::abs(sum-sum0) > 0)
		    std::cout << "(rho E) added from reset terms                 : " << sum-sum0 << " out of " << sum0 << std::endl;
#ifdef BL_LAZY
	    });
#endif
    }
}

void
PeleC::computeTemp(MultiFab& S, int ng)
{
  reset_internal_energy(S, ng);

#ifdef PELE_USE_EB
  auto const& fact = dynamic_cast<EBFArrayBoxFactory const&>(S.Factory());
  auto const& flags = fact.getMultiEBCellFlagFab();
#endif

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(S,true); mfi.isValid(); ++mfi)
  {
    const Box& bx = mfi.growntilebox(ng);

#ifdef PELE_USE_EB
    const auto& flag_fab = flags[mfi];
    FabType typ = flag_fab.getType(bx);
    if (typ == FabType::covered) {
      continue;
    }
#endif

    auto& Sfab = S[mfi];
    compute_temp(ARLIM_3D(bx.loVect()),ARLIM_3D(bx.hiVect()),BL_TO_FORTRAN_3D(Sfab));
  }
}

void
PeleC::set_special_tagging_flag(Real time)
{
  if (!do_special_tagging) return;

  MultiFab& S_new = get_new_data(State_Type);
  Real max_den = S_new.norm0(Density);

  int flag_was_changed = 0;
  pc_set_special_tagging_flag(max_den,&flag_was_changed);
  if (ParallelDescriptor::IOProcessor())
  {
    if (flag_was_changed == 1)
    {
      std::ofstream os("Bounce_time",std::ios::out);
      os << "T_Bounce " << time << std::endl;
      os.close();
    }
  }
}

Real
PeleC::getCPUTime()
{
  int numCores = ParallelDescriptor::NProcs();
#ifdef _OPENMP
  numCores = numCores*omp_get_max_threads();
#endif

  Real T = numCores*(ParallelDescriptor::second() - startCPUTime) +
    previousCPUTimeUsed;

  return T;
}


MultiFab&
PeleC::build_fine_mask()
{
  // Mask for zeroing covered cells
  BL_ASSERT(level > 0);

  if (!fine_mask.empty()) return fine_mask;

  BoxArray baf = parent->boxArray(level);
  baf.coarsen(crse_ratio);

  const BoxArray& bac = parent->boxArray(level-1);
  const DistributionMapping& dmc = parent->DistributionMap(level-1);
  fine_mask.define(bac,dmc,1,0,MFInfo(),FArrayBoxFactory());
  fine_mask.setVal(1.0);

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(fine_mask); mfi.isValid(); ++mfi)
  {
    auto& fab = fine_mask[mfi];

    const std::vector< std::pair<int,Box> >& isects = baf.intersections(fab.box());

    for (int ii = 0; ii < isects.size(); ++ii)
    {
      fab.setVal(0.0,isects[ii].second,0);
    }
  }

  return fine_mask;
}

const iMultiFab*
PeleC::build_interior_boundary_mask (int ng)
{
  for (int i = 0; i < ib_mask.size(); ++i)
  {
    if (ib_mask[i]->nGrow() == ng) {
      return ib_mask[i].get();
    }
  }

    //  If we got here, we need to build a new one
  if (ib_mask.size() == 0) {
    ib_mask.resize(0);
  }

  ib_mask.push_back(std::unique_ptr<iMultiFab>(new iMultiFab(grids, dmap, 1, ng,
                                                             MFInfo(),
                                                             DefaultFabFactory<IArrayBox>())));

  iMultiFab* imf = ib_mask.back().get();
  int ghost_covered_by_valid = 0;
  int other_cells = 1; // uncovered ghost, valid, and outside domain cells are set to 1

  imf->BuildMask(geom.Domain(), geom.periodicity(),
                 ghost_covered_by_valid, other_cells, other_cells, other_cells);

  return imf;
}

Real
PeleC::clean_state(MultiFab& S)
{
    // Enforce a minimum density.
  
  MultiFab temp_state(S.boxArray(), S.DistributionMap(), S.nComp(), S.nGrow(),MFInfo(),Factory());

  MultiFab::Copy(temp_state, S, 0, 0, S.nComp(), S.nGrow());

  Real frac_change = enforce_min_density(temp_state, S);

  //normalize_species(S);

  return frac_change;
}

Real
PeleC::clean_state(MultiFab& S,
                   MultiFab& S_old)
{
  // Enforce a minimum density.

  Real frac_change = enforce_min_density(S_old, S);

  //normalize_species(S);

  return frac_change;
}
