

#ifdef PELEC_USE_SPRAY
#include "PeleC.H"
#include "SprayParticles.H"

namespace {
bool virtual_particles_set = false;

// Containers for the real "active" Particles
SprayParticleContainer* SprayPC = nullptr;

// Container for temporary, virtual Particles
SprayParticleContainer* VirtPC = nullptr;

// Container for temporary, ghost Particles
SprayParticleContainer* GhostPC = nullptr;

SprayData sprayData;
// Indices for spray source MultiFab
int sprayRhoSrcIndx = 0;
int sprayMomSrcIndx = 1;
int sprayEngSrcIndx = 1 + AMREX_SPACEDIM;
int spraySpecSrcIndx = 2 + AMREX_SPACEDIM;
SprayComps scomps;

void
RemoveParticlesOnExit()
{
  delete SprayPC;
  SprayPC = nullptr;
  delete GhostPC;
  GhostPC = nullptr;
  delete VirtPC;
  VirtPC = nullptr;
}

std::string init_file;
int init_function = 1;
int particle_verbose = 0;
amrex::Real particle_cfl = 0.5;
amrex::Real wall_temp = 300.;
int mass_trans = 1;
int mom_trans = 1;
int plot_spray_src = 0;
} // namespace

int PeleC::write_spray_ascii_files = 0;
// momentum + density + fuel species + energy
int PeleC::num_spray_src = AMREX_SPACEDIM + 2 + SPRAY_FUEL_NUM;
std::string PeleC::sprayFuelNames[SPRAY_FUEL_NUM];

SprayParticleContainer*
PeleC::theSprayPC()
{
  return SprayPC;
}

SprayParticleContainer*
PeleC::theVirtPC()
{
  return VirtPC;
}

SprayParticleContainer*
PeleC::theGhostPC()
{
  return GhostPC;
}

void
PeleC::estTimeStepParticles(amrex::Real& est_dt)
{
  if (!do_spray_particles) {
    return;
  }
  BL_PROFILE("PeleC::particleEstTimeStep()");
  amrex::Real est_dt_particle = theSprayPC()->estTimestep(level, particle_cfl);

  if (est_dt_particle > 0) {
    est_dt = amrex::min<amrex::Real>(est_dt, est_dt_particle);
  }

  if (verbose && amrex::ParallelDescriptor::IOProcessor()) {
    if (est_dt_particle > 0) {
      amrex::Print() << "...estdt from particles at level " << level << ": "
                     << est_dt_particle << '\n';
    } else {
      amrex::Print() << "...there are no particles at level " << level << '\n';
    }
  }
}

void
PeleC::readSprayParams()
{
  amrex::ParmParse pp("pelec");

  pp.query("do_spray_particles", do_spray_particles);
  SprayParticleContainer::readSprayParams(
    particle_verbose, particle_cfl, wall_temp, mass_trans, mom_trans,
    write_spray_ascii_files, plot_spray_src, init_function, init_file,
    sprayData, sprayFuelNames);
}

void
PeleC::defineParticles()
{
  if (!do_spray_particles) {
    return;
  }
  // There must be at least as many fuel species in the spray as
  // there are species in the fluid
  if (SPRAY_FUEL_NUM > NUM_SPECIES) {
    amrex::Abort("Cannot have more spray fuel species than fluid species");
  }
  if (
    (std::is_same<
      pele::physics::PhysicsType::eos_type, pele::physics::eos::SRK>::value) ||
    (std::is_same<
      pele::physics::PhysicsType::eos_type,
      pele::physics::eos::Fuego>::value)) {
    for (int i = 0; i < SPRAY_FUEL_NUM; ++i) {
      for (int ns = 0; ns < NUM_SPECIES; ++ns) {
        std::string gas_spec = spec_names[ns];
        if (gas_spec == sprayFuelNames[i]) {
          sprayData.indx[i] = ns;
        }
      }
      if (sprayData.indx[i] < 0) {
        amrex::Print() << "Fuel " << sprayFuelNames[i]
                       << " not found in species list" << std::endl;
        amrex::Abort();
      }
    }
  } else {
    for (int ns = 0; ns < SPRAY_FUEL_NUM; ++ns) {
      sprayData.indx[ns] = 0;
    }
  }
  amrex::Vector<amrex::Real> fuelEnth(NUM_SPECIES);
  auto eos = pele::physics::PhysicsType::eos();
  eos.T2Hi(sprayData.ref_T, fuelEnth.data());
  for (int ns = 0; ns < SPRAY_FUEL_NUM; ++ns) {
    const int fspec = sprayData.indx[ns];
    sprayData.latent[ns] -= fuelEnth[fspec];
  }
  scomps.mass_tran = mass_trans;
  scomps.mom_tran = mom_trans;
  scomps.rhoIndx = PeleC::Density;
  scomps.momIndx = PeleC::Xmom;
  scomps.engIndx = PeleC::Eden;
  scomps.utempIndx = PeleC::Temp;
  scomps.specIndx = PeleC::FirstSpec;
  scomps.rhoSrcIndx = sprayRhoSrcIndx;
  scomps.momSrcIndx = sprayMomSrcIndx;
  scomps.specSrcIndx = spraySpecSrcIndx;
  scomps.engSrcIndx = sprayEngSrcIndx;
}

int
PeleC::sprayStateGhosts(const int amr_ncycle)
{
  if (!do_spray_particles) {
    return 0;
  }
  int finest_level = parent->finestLevel();
  return SprayParticleContainer::getStateGhostCells(
    level, finest_level, amr_ncycle);
}

void
PeleC::defineSpraySource(const int amr_ncycle)
{
  int finest_level = parent->finestLevel();
  int spray_source_ghosts = SprayParticleContainer::getSourceGhostCells(
    level, finest_level, amr_ncycle);
  // We must make a temporary spray source term to ensure number of ghost
  // cells are correct
  tmp_spray_source.define(
    grids, dmap, num_spray_src, spray_source_ghosts, amrex::MFInfo(),
    Factory());
}

void
PeleC::setupVirtualParticles()
{
  BL_PROFILE("PeleC::setupVirtualParticles()");
  if (theSprayPC() != nullptr && !virtual_particles_set) {
    if (level < parent->finestLevel()) {
      SprayParticleContainer::AoS virts;
      ((PeleC*)&parent->getLevel(level + 1))->setupVirtualParticles();
      theVirtPC()->CreateVirtualParticles(level + 1, virts);
      theVirtPC()->AddParticlesAtLevel(virts, level);

      theSprayPC()->CreateVirtualParticles(level + 1, virts);
      theVirtPC()->AddParticlesAtLevel(virts, level);
    }
    virtual_particles_set = true;
  }
}

void
PeleC::removeVirtualParticles()
{
  if (theVirtPC() != nullptr) {
    theVirtPC()->RemoveParticlesAtLevel(level);
  }
  virtual_particles_set = false;
}

void
PeleC::setupGhostParticles(int ngrow)
{
  BL_PROFILE("PeleC::setupGhostParticles()");
  AMREX_ASSERT(level < parent->finestLevel());
  if (PeleC::theSprayPC() != nullptr) {
    SprayParticleContainer::AoS ghosts;
    theSprayPC()->CreateGhostParticles(level, ngrow, ghosts);
    theGhostPC()->AddParticlesAtLevel(ghosts, level + 1, ngrow);
  }
}

void
PeleC::removeGhostParticles()
{
  if (theGhostPC() != nullptr) {
    theGhostPC()->RemoveParticlesAtLevel(level);
  }
}

// Create new particle data
void
PeleC::createDataParticles()
{
  SprayPC =
    new SprayParticleContainer(parent, &phys_bc, sprayData, scomps, wall_temp);
  theSprayPC()->SetVerbose(particle_verbose);
  VirtPC =
    new SprayParticleContainer(parent, &phys_bc, sprayData, scomps, wall_temp);
  GhostPC =
    new SprayParticleContainer(parent, &phys_bc, sprayData, scomps, wall_temp);
}

// Initialize the particles on the grid at level 0
void
PeleC::initParticles()
{
  BL_PROFILE("PeleC::initParticles()");

  if (level > 0) {
    return;
  }

  // Make sure to call RemoveParticlesOnExit() on exit.
  amrex::ExecOnFinalize(RemoveParticlesOnExit);

  if (do_spray_particles) {
    AMREX_ASSERT(theSprayPC() == 0);
    createDataParticles();

    if (!init_file.empty()) {
      theSprayPC()->InitFromAsciiFile(init_file, NSR_SPR + NAR_SPR);
    } else if (init_function > 0) {
      const ProbParmHost* lprobparm = prob_parm_host;
      const ProbParmDevice* lprobparm_d = h_prob_parm_device;
      theSprayPC()->InitSprayParticles(*lprobparm, *lprobparm_d);
    } else {
      amrex::Abort(
        "Must initialize spray particles with particles.init_function or "
        "particles.init_file");
    }
  }
}

void
PeleC::postRestartParticles(bool is_checkpoint)
{
  if (level > 0) {
    defineSpraySource(parent->MaxRefRatio(level - 1));
    return;
  }

  if (do_spray_particles) {
    defineSpraySource(1);
    AMREX_ASSERT(SprayPC == 0);
    createDataParticles();

    // Make sure to call RemoveParticlesOnExit() on exit.
    amrex::ExecOnFinalize(RemoveParticlesOnExit);
    {
      theSprayPC()->Restart(
        parent->theRestartFile(), "particles", is_checkpoint);
      amrex::Gpu::Device::streamSynchronize();
    }
  }
}

void
PeleC::particleMKD(
  amrex::Real time,
  amrex::Real dt,
  int sub_iteration,
  int sub_ncycle,
  int amr_ncycle)
{
  if (sub_iteration != 0) {
    return;
  }
  old_sources[spray_src]->setVal(0.);
  tmp_spray_source.setVal(0.);
  // Setup ghost particles for use in finer levels. Note that ghost
  // particles that will be used by this level have already been created,
  // the particles being set here are only used by finer levels.
  int finest_level = parent->finestLevel();

  if (level < finest_level) {
    // Setup the virtual particles that represent particles on finer levels
    setupVirtualParticles();
    // Make a copy of the particles on this level into ghost particles
    // for the finer level
    int finer_ref = parent->MaxRefRatio(level);
    // Determine the number of ghost cells on the next level we need for making
    // ghost particles
    int ghost_width = SprayParticleContainer::getGhostPartCells(
      level + 1, finest_level, finer_ref);
    setupGhostParticles(ghost_width);
  }

  // Advance the particle velocities to the half-time and the positions to
  // the new time
  if (particle_verbose) {
    amrex::Print()
      << "moveKickDrift ... updating particle positions and velocity\n";
  }
  int spray_state_ghosts = sprayStateGhosts(amr_ncycle);
  int spray_source_ghosts = tmp_spray_source.nGrow();
  auto const* ltransparm = PeleC::trans_parms.device_trans_parm();
  // We will make a temporary copy of the source term array inside
  // moveKickDrift and we are only going to use the spray force out to one ghost
  // cell so we need only define spray_force_old with one ghost cell

  AMREX_ASSERT(old_sources[spray_src]->nGrow() >= 1);

  // Do the valid particles themselves
  bool isVirt = false;
  bool isGhost = false;
  bool doMove = true;
  theSprayPC()->moveKickDrift(
    Sborder, tmp_spray_source, level, dt, time, isVirt, isGhost,
    spray_state_ghosts, spray_source_ghosts, doMove, ltransparm);

  // Only need the coarsest virtual particles here.
  if (level < finest_level) {
    isVirt = true;
    isGhost = false;
    theVirtPC()->moveKickDrift(
      Sborder, tmp_spray_source, level, dt, time, isVirt, isGhost,
      spray_state_ghosts, spray_source_ghosts, doMove, ltransparm);
  }

  // Miiiight need all Ghosts
  if (theGhostPC() != nullptr && level != 0) {
    isVirt = false;
    isGhost = true;
    theGhostPC()->moveKickDrift(
      Sborder, tmp_spray_source, level, dt, time, isVirt, isGhost,
      spray_state_ghosts, spray_source_ghosts, doMove, ltransparm);
  }
  // Must call transfer source after moveKick and moveKickDrift
  // on all particle types
  theSprayPC()->transferSource(
    spray_source_ghosts, level, tmp_spray_source, *old_sources[spray_src]);
}

void
PeleC::particleMK(
  amrex::Real time,
  amrex::Real dt,
  int sub_iteration,
  int sub_ncycle,
  int amr_ncycle)
{
  if (sub_iteration != sub_ncycle - 1 && sub_ncycle != 0) {
    return;
  }
  int spray_state_ghosts = sprayStateGhosts(amr_ncycle);
  int spray_source_ghosts = tmp_spray_source.nGrow();
  new_sources[spray_src]->setVal(0.);
  auto const* ltransparm = PeleC::trans_parms.device_trans_parm();
  if (particle_verbose) {
    amrex::Print() << "moveKick ... updating velocity only\n";
  }
  bool isVirt = false;
  bool isGhost = false;
  theSprayPC()->moveKick(
    Sborder, tmp_spray_source, level, dt, time, isVirt, isGhost,
    spray_state_ghosts, spray_source_ghosts, ltransparm);

  if (level < parent->finestLevel()) {
    isVirt = true;
    isGhost = false;
    theVirtPC()->moveKick(
      Sborder, tmp_spray_source, level, dt, time, isVirt, isGhost,
      spray_state_ghosts, spray_source_ghosts, ltransparm);
  }

  if (theGhostPC() != nullptr && level != 0) {
    isVirt = false;
    isGhost = true;
    theGhostPC()->moveKick(
      Sborder, tmp_spray_source, level, dt, time, isVirt, isGhost,
      spray_state_ghosts, spray_source_ghosts, ltransparm);
  }
  theSprayPC()->transferSource(
    spray_source_ghosts, level, tmp_spray_source, *new_sources[spray_src]);
}

void
PeleC::postTimeStepParticles(int iteration)
{
  const int finest_level = parent->finestLevel();
  const int ncycle = parent->nCycle(level);
  if (do_spray_particles) {
    // Remove virtual particles at this level if we have any.
    if (theVirtPC() != nullptr) {
      removeVirtualParticles();
    }

    // Remove Ghost particles on the final iteration
    if (iteration == ncycle) {
      removeGhostParticles();
    }

    // Do particle injection
    int nstep = parent->levelSteps(0);
    amrex::Real dtlev = parent->dtLevel(0);
    amrex::Real cumtime = parent->cumTime() + dtlev;
    const ProbParmHost* lprobparm = prob_parm_host;
    const ProbParmDevice* lprobparm_d = h_prob_parm_device;
    BL_PROFILE_VAR("SprayParticles::injectParticles()", INJECT_SPRAY);
    bool injectParts = theSprayPC()->injectParticles(
      cumtime, dtlev, nstep, level, finest_level, *lprobparm, *lprobparm_d);
    BL_PROFILE_VAR_STOP(INJECT_SPRAY);
    // Sync up if we're level 0 or if we have particles that may have moved
    // off the next finest level and need to be added to our own level, or
    // if we injected particles
    if (
      (iteration < ncycle && level < finest_level) || level == 0 ||
      injectParts) {
      // TODO: Determine how many ghost cells to use here
      int nGrow = iteration;
      theSprayPC()->Redistribute(level, theSprayPC()->finestLevel(), nGrow);
    }
  }
}

void
PeleC::postInitParticles()
{
  amrex::Real dtlev = parent->dtLevel(level);
  amrex::Real cumtime = parent->cumTime();
  if (do_spray_particles) {
    const ProbParmHost* lprobparm = prob_parm_host;
    const ProbParmDevice* lprobparm_d = h_prob_parm_device;
    BL_PROFILE_VAR("SprayParticles::injectParticles()", INJECT_SPRAY);
    bool injectParts = theSprayPC()->injectParticles(
      cumtime, dtlev, 0, level, parent->finestLevel(), *lprobparm,
      *lprobparm_d);
    BL_PROFILE_VAR_STOP(INJECT_SPRAY);
    if (injectParts) {
      theSprayPC()->Redistribute(level, theSprayPC()->finestLevel(), 0);
    }
  }
}

void
PeleC::particle_redistribute(int lbase, bool init_part)
{
  if (!do_spray_particles) {
    return;
  }
  BL_PROFILE("PeleC::particle_redistribute()");
  int flev = parent->finestLevel();
  if (theSprayPC()) {
    // If we are calling with init_part = true, then we want to force the
    // redistribute without checking whether the grids have changed.
    if (init_part) {
      theSprayPC()->Redistribute(lbase);
      return;
    }

    // These are usually the BoxArray and DMap from the last regridding.
    static amrex::Vector<amrex::BoxArray> ba;
    static amrex::Vector<amrex::DistributionMapping> dm;

    bool changed = false;

    while (parent->getAmrLevels()[flev] == nullptr) {
      flev--;
    }

    if (ba.size() != flev + 1) {
      ba.resize(flev + 1);
      dm.resize(flev + 1);
      changed = true;
    } else {
      for (int i = 0; i <= flev && !changed; i++) {
        // Check if BoxArrays have changed during regridding
        if (ba[i] != parent->boxArray(i)) {
          changed = true;
        }

        if (!changed) {
          // Check DistributionMaps have changed during regridding
          if (dm[i] != parent->getLevel(i).get_new_data(0).DistributionMap()) {
            changed = true;
          }
        }
      }
    }

    if (changed) {
      // We only need to call Redistribute if the BoxArrays or DistMaps have
      // changed.
      // We also only call it for particles >= lbase. This is
      // because if we called redistribute during a subcycle, there may be
      // particles not in the proper position on coarser levels.
      if (verbose && amrex::ParallelDescriptor::IOProcessor()) {
        amrex::Print() << "Calling redistribute because grid has changed "
                       << '\n';
      }
      if (flev == 0) {
        // Do a local redistribute
        theSprayPC()->Redistribute(lbase, theSprayPC()->finestLevel(), 0, 1);
      } else {
        theSprayPC()->Redistribute(lbase, theSprayPC()->finestLevel(), 0);
      }

      // Use the new BoxArray and DistMap to define ba and dm for next time.
      for (int i = 0; i <= flev; i++) {
        ba[i] = parent->boxArray(i);
        dm[i] = parent->getLevel(i).get_new_data(0).DistributionMap();
      }
    } else {
      if (verbose && amrex::ParallelDescriptor::IOProcessor()) {
        amrex::Print()
          << "NOT calling redistribute because grid has NOT changed " << '\n';
      }
    }
  }
}

#endif
