#ifndef _PeleC_H_
#define _PeleC_H_

#include <AMReX_BC_TYPES.H>
#include <AMReX_AmrLevel.H>
#include <AMReX_iMultiFab.H>
#ifdef PELE_USE_EB
#include <AMReX_EBFArrayBox.H>
#include <EBStencilTypes.H>
#include <SparseData.H>
#include <AMReX_EBFluxRegister.H>
#if BL_SPACEDIM > 1
#endif
#else
#include <AMReX_YAFluxRegister.H>
#endif

#ifdef BL_LAZY
#include <AMReX_Lazy.H>
#endif

#ifdef AMREX_PARTICLES
#include <AMReX_Particles.H>
#include <SprayParticles.H>
#endif

#include <Filter.H>

#include <iostream>

using std::istream;
using std::ostream;

#define UserBC 6


enum StateType { State_Type = 0
#ifdef REACTIONS
                 ,Reactions_Type
#endif
                 , Work_Estimate_Type
};

// Create storage for all source terms.

enum sources { ext_src = 0,
               forcing_src,
#ifdef AMREX_PARTICLES
               spray_src,
#endif
               diff_src,
               les_src,
#ifdef USE_MASA
               mms_src,
#endif
               num_src };

//
// AmrLevel-derived class for hyperbolic conservation equations for stellar media
//

class PeleC
    :
    public amrex::AmrLevel
{
public:
    //
    //Default constructor.  Builds invalid object.
    //
    PeleC ();
    //
    //The basic constructor.
    //
    PeleC (amrex::Amr&            papa,
           int                    lev,
           const amrex::Geometry& level_geom,
           const amrex::BoxArray& bl,
           const amrex::DistributionMapping& dm,
           amrex::Real                   time);
    //
    //The destructor.
    //
    virtual ~PeleC ();
    //
    //Restart from a checkpoint file.
    //
    virtual void restart (amrex::Amr&     papa,
                          istream&        is,
			  bool            bReadSpecial = false) override;
    //
    //This is called only when we restart from an old checkpoint.
    //
    virtual void set_state_in_checkpoint (amrex::Vector<int>& 
                                          state_in_checkpoint) override;
    //
    //Call AmrLevel::checkPoint and then add radiation info
    //
    virtual void checkPoint(const std::string& dir,
                            std::ostream&      os,
                            amrex::VisMF::How         how,
                            bool               dump_old) override;

    /*A string written as the first item in writePlotFile() at
      level zero. It is so we can distinguish between different
      types of plot files. For PeleC it has the form: PeleC-Vnnn
    */
    virtual std::string thePlotFileType () const override;

    virtual void setPlotVariables () override;
    //
    //Write a plotfile to specified directory.
    //
    virtual void writePlotFile (const std::string& dir,
                                ostream&       os,
                                amrex::VisMF::How     how) override;
    virtual void writeSmallPlotFile (const std::string& dir,
				     ostream&       os,
				     amrex::VisMF::How     how) override;
    void writeJobInfo (const std::string& dir);
    static void writeBuildInfo (std::ostream& os);

    //
    // Define data descriptors.
    //
    static void variableSetUp ();
    //
    // Cleanup data descriptors at end of run.
    //
    static void variableCleanUp ();
    //
    // Initialize grid data at problem start-up.
    //
    virtual void initData () override;

#ifdef AMREX_PARTICLES
    //
    // Read particle-related inputs
    //
    static void read_particle_params ();
    //
    // Initialize particle locations and velocities (and strengths if relevant)
    //
    void init_particles ();
    //
    // Write particles in checkpoint directories
    //
    void ParticleCheckPoint (const std::string& dir);
    //
    // Timestamp particles
    //
    void TimestampParticles (int ngrow); 
    //
    // Advance the particles by dt
    //
    void advance_particles (int iteration, amrex::Real time, amrex::Real dt);
    //
    // Default verbosity of Particle class
    //
    static int particle_verbose;

    //
    // Write particles in plotfile directories
    //
    void particle_plot_file(const std::string& dir);

    //
    // How to initialize at restart
    //
    void particle_post_restart(const std::string& restart_file, bool is_checkpoint = true);

    //
    // Redistribute
    //
    virtual void particle_redistribute(int lbase = 0, bool init_part = false);

    //
    // Setup virtual particles if necessary
    //
    void setup_virtual_particles();

    //
    // Remove virtual particles if necessary
    //
    void remove_virtual_particles();
    //
    // Setup ghost particles (for finer levels) if necessary
    //
    void setup_ghost_particles(int ngrow);

    //
    // Remove ghost particles (for this level) if necessary
    //
    void remove_ghost_particles();

    //
    // Time step control based on particles
    //
    void particle_est_time_step(amrex::Real& est_dt);

    //
    //
    // Derived quantities associated with particles
    //
    std::unique_ptr<amrex::MultiFab> particle_derive (const std::string& name, amrex::Real time, int ngrow);

    // Default cfl of particles in Particle class
    //
    static amrex::Real particle_cfl;

    //
    // Should we write particles into plotfiles?
    //
    static int write_particles_in_plotfile;

    void set_spray_grid_info(int amr_iteration,
                             int amr_ncycle,
                             int ghost_width,
                             int where_width,
                             int spray_n_grow,
                             int tmp_src_width);
#endif

    //
    // Get problem metrics.
    //
    amrex::MultiFab* Area ();
    amrex::MultiFab& Area (int dir);
    //
    // The volume of the multifab.
    //
    amrex::MultiFab& Volume ();

    // Set time levels of state data.
    //
    virtual void setTimeLevel (amrex::Real time,
                               amrex::Real dt_old,
                               amrex::Real dt_new) override;
    //
    // Initialize data on this level from another PeleC (during regrid).
    //
    virtual void init (amrex::AmrLevel& old) override;
    //
    // Initialize data on this level after regridding if old level did not previously exist
    //
    virtual void init () override;
    
    /**
     * Initialize EB geometry for finest_level and level grids for
     * other levels
     */
#ifdef PELE_USE_EB
    static bool DoMOLLoadBalance () {return do_mol_load_balance;}

    //for the Amr class to do timed load balances. 
    virtual int WorkEstType () override { return Work_Estimate_Type; }

    const amrex::MultiFab& volFrac () const { return vfrac; }

    void init_eb (const amrex::Geometry& level_geom, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm);

    void initialize_eb2_structs();

#ifdef PELE_UNIT_TEST_DN
    void test_dn();
#endif

    void define_body_state();

    void set_body_state(amrex::MultiFab& S);

    void zero_in_body(amrex::MultiFab& S) const;
#endif
    //
    // Proceed with next timestep?
    //
    virtual int okToContinue () override;
    //
    // Advance grids at this level in time.
    //
    virtual amrex::Real advance (amrex::Real time,
                          amrex::Real dt,
                          int  iteration,
                          int  ncycle) override;

    amrex::Real do_mol_advance(amrex::Real time,
                               amrex::Real dt,
                               int  amr_iteration,
                               int  amr_ncycle);

    amrex::Real do_sdc_advance(amrex::Real time,
                               amrex::Real dt,
                               int  amr_iteration,
                               int  amr_ncycle);

    void initialize_sdc_advance(amrex::Real time, amrex::Real dt, int amr_iteration, int amr_ncycle);

    void finalize_sdc_advance(amrex::Real time, amrex::Real dt, int amr_iteration, int amr_ncycle);

    amrex::Real do_sdc_iteration (amrex::Real time,
                           amrex::Real dt,
                           int  amr_iteration,
                           int  amr_ncycle,
                           int  sub_iteration,
                           int  sub_ncycle);

    void initialize_sdc_iteration(amrex::Real time, amrex::Real dt,
                                  int amr_iteration, int amr_ncycle, int sdc_iteration, int sdc_ncycle);

    void finalize_sdc_iteration(amrex::Real time, amrex::Real dt,
                                int amr_iteration, int amr_ncycle, int sdc_iteration, int sdc_ncycle);

    void construct_Snew(amrex::MultiFab& S_new, const amrex::MultiFab& S_old, amrex::Real dt);

    void construct_hydro_source(const amrex::MultiFab& S, amrex::Real time, amrex::Real dt, int amr_iteration, int amr_ncycle, int sub_iteration, int sub_ncycle);
 
    void construct_old_source(int src, amrex::Real time, amrex::Real dt, int amr_iteration, int amr_ncycle, int sub_iteration, int sub_ncycle);

    void construct_new_source(int src, amrex::Real time, amrex::Real dt, int amr_iteration, int amr_ncycle, int sub_iteration, int sub_ncycle);

    void sum_of_sources(amrex::MultiFab& source);

    void construct_old_ext_source(amrex::Real time, amrex::Real dt);

    void construct_new_ext_source(amrex::Real time, amrex::Real dt);

    void fill_ext_source(amrex::Real time, amrex::Real dt,
                         const amrex::MultiFab& S_old,
                         const amrex::MultiFab& S_new,
                         amrex::MultiFab& ext_src, int ng);

    void construct_old_forcing_source(amrex::Real time, amrex::Real dt);

    void construct_new_forcing_source(amrex::Real time, amrex::Real dt);

    void fill_forcing_source(amrex::Real time, amrex::Real dt,
                             const amrex::MultiFab& S_old,
                             const amrex::MultiFab& S_new,
                             amrex::MultiFab& forcing_src, int ng);

#ifdef USE_MASA
    void construct_old_mms_source(amrex::Real time);

    void construct_new_mms_source(amrex::Real time);

    void fill_mms_source (amrex::Real time, const amrex::MultiFab& S, amrex::MultiFab& mms_src, int ng);
#endif

    void getLESTerm (amrex::Real time, amrex::Real dt, amrex::MultiFab& LESTerm, amrex::Real flux_factor);

    void getSmagorinskyLESTerm (amrex::Real time, amrex::Real dt, amrex::MultiFab& LESTerm, amrex::Real flux_factor);

    void getDynamicSmagorinskyLESTerm (amrex::Real time, amrex::Real dt, amrex::MultiFab& LESTerm, amrex::Real flux_factor);

    void construct_old_les_source(amrex::Real time, amrex::Real dt, int sub_iteration, int sub_ncycle);

    void construct_new_les_source(amrex::Real time, amrex::Real dt, int sub_iteration, int sub_ncycle);

    static void set_active_sources();

    //
    // Estimate time step.
    //
    amrex::Real estTimeStep (amrex::Real dt_old);
    //
    // Compute initial time step.
    //
    amrex::Real initialTimeStep ();
    //
    // Compute initial `dt'.
    //
    virtual void computeInitialDt (int                                 finest_level,
                                   int                                 sub_cycle,
                                   amrex::Vector<int>&                  n_cycle,
                                   const amrex::Vector<amrex::IntVect>& ref_ratio,
                                   amrex::Vector<amrex::Real>&                 dt_level,
                                   amrex::Real                                stop_time) override;
    //
    // Compute new `dt'.
    //
    virtual void computeNewDt (int                                 finest_level,
                               int                                 sub_cycle,
                               amrex::Vector<int>&                  n_cycle,
                               const amrex::Vector<amrex::IntVect>& ref_ratio,
                               amrex::Vector<amrex::Real>&                 dt_min,
                               amrex::Vector<amrex::Real>&                 dt_level,
                               amrex::Real                                stop_time,
                               int                                 post_regrid_flag) override;
    //
    // Allocate data at old time.
    //
    virtual void allocOldData () override;
    //
    // Remove data at old time.
    //
    virtual void removeOldData () override;
    //
    // Passes some data about the grid to a Fortran module.
    //
    void setGridInfo();
    //
    // Print information about energy budget.
    //
    void do_energy_diagnostics();
    //
    // Do work after timestep().
    //
    virtual void post_timestep (int iteration) override;
    //
    // Contains operations to be done only after a full coarse timestep.
    //
    virtual void postCoarseTimeStep (amrex::Real cumtime) override;
    //
    // Do work after regrid().
    //
    virtual void post_regrid (int lbase,
                              int new_finest) override;
    //
    // Do work after a restart().
    //
    virtual void post_restart () override;
    //
    // Do work after init().
    //
    virtual void post_init (amrex::Real stop_time) override;
    //
    // Error estimation for regridding.
    //
    virtual void errorEst (amrex::TagBoxArray& tb,
                           int                 clearval,
                           int                 tagval,
                           amrex::Real                time,
                           int                 n_error_buf = 0,
                           int                 ngrow = 0) override;

    // Returns a MultiFab containing the derived data for this level.
    // The user is responsible for deleting this pointer when done
    // with it.  If ngrow>0 the MultiFab is built on the appropriately
    // grown BoxArray.
    virtual std::unique_ptr<amrex::MultiFab> derive (const std::string& name,
                                     amrex::Real               time,
                                     int                ngrow) override;

    // This version of derive() fills the dcomp'th component of mf with the derived quantity.
    virtual void derive (const std::string& name,
                         amrex::Real               time,
                         amrex::MultiFab&   mf,
                         int                dcomp) override;

    static int numGrow();

#ifdef REACTIONS
    void react_state(amrex::Real time, amrex::Real dt, bool init=false, amrex::MultiFab* A_aux = nullptr);
#endif

    void reset_internal_energy (amrex::MultiFab& State, int ng);

    void computeTemp (amrex::MultiFab& State, int ng);

    void getMOLSrcTerm (const amrex::MultiFab& S,
                        amrex::MultiFab&       MOLSrcTerm,
                        amrex::Real time,
                        amrex::Real                   dt,
                        amrex::Real                   flux_factor);

    amrex::Real volWgtSum (const std::string& name, amrex::Real time, bool local=false, bool finemask=true);
    amrex::Real volWgtSquaredSum (const std::string& name, amrex::Real time, bool local=false);
    amrex::Real volWgtSumMF (amrex::MultiFab* mf, int comp, bool local=false, bool finemask=false);


#ifdef USE_MASA
    amrex::Real mms_volWgtSquaredSumDiff (int comp, amrex::Real time, bool local=false);
#endif

    void set_special_tagging_flag (amrex::Real time);

    static int       NUM_STATE;
    static int       Density, Xmom, Ymom, Zmom, Eden, Eint, Temp;

    static int       FirstAdv,  NumAdv;
    static int       FirstSpec, NumSpec;
    static int       FirstAux,  NumAux;

    /// Index locations into particle state
    static int       pstate_loc, pstate_vel, pstate_T, pstate_dia, pstate_rho, pstate_spc, n_pstate;

    /// Component locations for field variables that get interpolated to particle position
    static int       pfld_vel, pfld_rho, pfld_T, pfld_p, pfld_spc, n_pfld;
    //
    // This MultiFab is on the coarser level.  This is useful for the coarser level
    // to mask out the finer level.  We only build this when it is needed.
    // This coarse MultiFab has to live on the fine level because it must be updated
    // even when only the fine level changes.
    //
    amrex::MultiFab fine_mask;
    amrex::MultiFab& build_fine_mask();

#ifdef PELE_USE_EB
    static bool      no_eb_in_domain;

    amrex::EBFluxRegister     flux_reg;

    amrex::EBFluxRegister& getFluxReg ();
    amrex::EBFluxRegister& getFluxReg (int lev);

    amrex::EBFluxRegister     pres_reg;

    amrex::EBFluxRegister& getPresReg ();
    amrex::EBFluxRegister& getPresReg (int lev);

    // static EB helper functions called from outside class
    static std::string convertIntGG(int number);
    static bool ebInitialized();
#else
    amrex::YAFluxRegister     flux_reg;

    amrex::YAFluxRegister& getFluxReg ();
    amrex::YAFluxRegister& getFluxReg (int lev);

    amrex::YAFluxRegister     pres_reg;

    amrex::YAFluxRegister& getPresReg ();
    amrex::YAFluxRegister& getPresReg (int lev);
#endif
    void avgDown ();
    void avgDown (int state_indx);

protected:

    amrex::iMultiFab level_mask;
#ifdef PELEC_USE_EB
    std::array<const amrex::MultiCutFab*,AMREX_SPACEDIM> areafrac;
#endif
    //
    // Build a mask that ghost cells overlapping with interior cells in the same multifab
    // are set to 0, whereas others are set to 1.
    //
    amrex::Vector<std::unique_ptr<amrex::iMultiFab>> ib_mask;
    const amrex::iMultiFab* build_interior_boundary_mask (int ng);

    ///
    /// A state array with ghost zones.
    ///
    amrex::MultiFab Sborder;
    ///
    /// Source terms to the hydrodynamics solve.
    ///
    amrex::MultiFab sources_for_hydro;

    ///
    /// Source term representing hydrodynamics update.
    ///
    amrex::MultiFab hydro_source;

    ///
    /// Non-hydro source terms.
    ///
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> old_sources;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> new_sources;

    ///
    ///  Call extern/networks/.../network.f90::network_init()
    ///
    static void init_extern ();

    static void init_network ();
    static void close_network ();

#ifdef REACTIONS
    static void init_reactor ();
    static void close_reactor ();
#endif

    static void init_transport ();
    static void close_transport ();
  
    void init_les ();
    void init_filters ();

#ifdef USE_MASA
    static void init_mms();
#endif

    static void clear_prob ();

    static void read_params ();

    PeleC& getLevel (int lev);

    void reflux ();

    void normalize_species (amrex::MultiFab& S_new);

    void enforce_consistent_e (amrex::MultiFab& S);

    amrex::Real enforce_min_density (amrex::MultiFab& S_old, amrex::MultiFab& S_new);

    amrex::Real clean_state (amrex::MultiFab& state);

    amrex::Real clean_state (amrex::MultiFab& state, amrex::MultiFab& state_old);

    void buildMetrics ();


    /// integrate derived quantities over domain

    amrex::Real sumDerive (const std::string& name, amrex::Real time, bool local=false);

    void sum_integrated_quantities ();

    void write_info ();

    void stopJob ();

    //
    // The data.
    //
    amrex::MultiFab             volume;
    amrex::MultiFab             area[BL_SPACEDIM];
    amrex::MultiFab             dLogArea[1];
    amrex::Vector< amrex::Vector<amrex::Real> > radius;
    int m_ngrow;
    //
    // Static data members.
    //
#include <pelec_params.H>

    static bool      signalStopJob;
    static bool      dump_old;
    static int       radius_grow;
    static int       verbose;
    static amrex::BCRec     phys_bc;
    static int       NUM_GROW, QTHERM, QVAR, cQRHO, cQU, cQV, cQW, cQGAME, cQPRES, cQREINT, cQTEMP, cQFA, cQFS, cQFX, NQAUX, cQGAMC, cQC, cQCSML, cQDPDR, cQDPDE, cQRSPEC;

    /// for keeping track of mass changes from negative density resets
    static amrex::Real      frac_change;

    /// For keeping track of fluid quantities lost at physical grid boundaries.
    /// This should persist through restarts, but right now only on level 0.
    static const int n_lost = 8;
    amrex::Real             material_lost_through_boundary_cumulative[n_lost];
    amrex::Real             material_lost_through_boundary_temp[n_lost];

    /// for keeping track of the amount of CPU time used -- this will persist
    /// after restarts
    static amrex::Real      previousCPUTimeUsed;
    static amrex::Real      startCPUTime;

    static amrex::Real getCPUTime();

    static int       diffuse_temp;
    static int       diffuse_enth;
    static int       diffuse_spec;
    static int       diffuse_vel;
    static amrex::Real      diffuse_cutoff_density;
    static bool      do_diffuse;

    static int       les_model;
    static int       les_filter_type;
    static int       les_filter_fgr;
    Filter           les_filter;
    int              nGrowF;
    static int       les_test_filter_type;
    static int       les_test_filter_fgr;
    amrex::MultiFab  LES_Coeffs;
    static int       comp_Cs2;
    static int       comp_CI;
    static int       comp_PrT;
    static int       nCompC;
    amrex::MultiFab  filtered_les_source;

#ifdef USE_MASA
    static bool      mms_initialized;
    bool             mms_src_evaluated;
    amrex::MultiFab  mms_source;
#endif

#ifdef AMREX_PARTICLES
    static int       do_spray_particles;
#endif

    bool             FillPatchedOldState_ok;

    // There can be only one Diffusion object, it covers all levels:
    static class Diffusion *diffusion;

#ifdef AMREX_PARTICLES
    static SprayParticleContainer* SprayPC;
    static SprayParticleContainer* theSprayPC ();
    static SprayParticleContainer* theVirtPC ();
    static SprayParticleContainer* theGhostPC ();
#endif

    // Name of the probin file and its length.
    static std::string probin_file;

    static amrex::IntVect hydro_tile_size;

    static int num_state_type;

    static std::vector<std::string> spec_names;

    static std::vector<int> src_list;

/* problem-specific includes */
#include <Problem.H>

    static int nGrowTr;

#ifdef PELE_USE_EB

    amrex::FabArray<amrex::BaseFab<int>> ebmask;
    static std::vector<amrex::Real> body_state;
    static bool body_state_set;
    static bool eb_initialized;


    amrex::MultiFab vfrac;

    std::vector<std::vector<EBBndryGeom>> sv_eb_bndry_geom;
    std::vector<std::vector<EBBndrySten>> sv_eb_bndry_grad_stencil;
    std::vector<std::vector<FaceSten>> flux_interp_stencil[BL_SPACEDIM];

    std::vector<SparseData<amrex::Real,EBBndrySten>> sv_eb_flux;
    std::vector<SparseData<amrex::Real,EBBndrySten>> sv_eb_bcval;
#endif
  static bool do_react_load_balance;
  static bool do_mol_load_balance;

};

    void pc_bcfill_hyp (amrex::Box const& bx, amrex::FArrayBox& data,
                 const int dcomp, const int numcomp,
                 amrex::Geometry const& geom, const amrex::Real time,
                 const amrex::Vector<amrex::BCRec>& bcr, const int bcomp,
                 const int scomp);
    
    void pc_reactfill_hyp (amrex::Box const& bx, amrex::FArrayBox& data,
                 const int dcomp, const int numcomp,
                 amrex::Geometry const& geom, const amrex::Real time,
                 const amrex::Vector<amrex::BCRec>& bcr, const int bcomp,
                 const int scomp);


//
// Inlines.
//

inline
int
PeleC::numGrow()
{
    return NUM_GROW;
}

inline
amrex::MultiFab*
PeleC::Area ()
{
    return area;
}

inline
amrex::MultiFab&
PeleC::Area (int dir)
{
    return area[dir];
}

inline
amrex::MultiFab&
PeleC::Volume ()
{
    return volume;
}

inline
PeleC&
PeleC::getLevel (int lev)
{
    return *(PeleC *) &parent->getLevel(lev);
}

inline
void
PeleC::stopJob()
{
    signalStopJob = true;
}

#ifdef PELE_USE_EB
inline
amrex::EBFluxRegister&
PeleC::getFluxReg ()
{
    return flux_reg;
}

inline
amrex::EBFluxRegister&
PeleC::getFluxReg (int lev)
{
    return getLevel(lev).getFluxReg();
}

inline
amrex::EBFluxRegister&
PeleC::getPresReg ()
{
    return pres_reg;
}

inline
amrex::EBFluxRegister&
PeleC::getPresReg (int lev)
{
    return getLevel(lev).getPresReg();
}

#else

inline
amrex::YAFluxRegister&
PeleC::getFluxReg ()
{
    return flux_reg;
}

inline
amrex::YAFluxRegister&
PeleC::getFluxReg (int lev)
{
    return getLevel(lev).getFluxReg();
}

inline
amrex::YAFluxRegister&
PeleC::getPresReg ()
{
    return pres_reg;
}

inline
amrex::YAFluxRegister&
PeleC::getPresReg (int lev)
{
    return getLevel(lev).getPresReg();
}
#endif


#endif /*_PeleC_H_*/
