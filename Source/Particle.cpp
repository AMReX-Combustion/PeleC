#include <iomanip>
#include <vector>
#include <algorithm>
#include <string>

#include "PeleC.H"

using namespace amrex;

#ifdef PELEC_SPRAY

namespace {
bool virtual_particles_set = false;

// Containers for the real "active" Particles
SprayParticleContainer* SprayPC = 0;

// Container for temporary, virtual Particles
SprayParticleContainer* VirtPC = 0;

// Container for temporary, ghost Particles
SprayParticleContainer* GhostPC = 0;

Gpu::HostVector<Real> sprayCritT;
Gpu::HostVector<Real> sprayBoilT;
Gpu::HostVector<Real> sprayCp;
Gpu::HostVector<Real> sprayLatent;
Gpu::HostVector<int> sprayIndxMap;

void
RemoveParticlesOnExit()
{
  delete SprayPC;
  SprayPC = 0;
  delete GhostPC;
  GhostPC = 0;
  delete VirtPC;
  VirtPC = 0;
}
} // namespace

int PeleC::particle_verbose = 0;
Real PeleC::particle_cfl = 0.4;

int PeleC::write_particle_plotfiles = 1;
int PeleC::write_spray_ascii_files = 1;
int PeleC::particle_mass_tran = 0;
int PeleC::particle_heat_tran = 0;
int PeleC::particle_mom_tran = 0;
Vector<std::string> PeleC::sprayFuelNames;
Real PeleC::sprayRefT;

namespace {
std::string particle_init_file;
int particle_init_uniform = 0;
std::string timestamp_dir;
std::vector<int> timestamp_indices;
} // namespace

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
PeleC::particleEstTimeStep(Real& est_dt)
{
  BL_PROFILE("PeleC::particleEstTimeStep()");
  Real est_dt_particle = theSprayPC()->estTimestep(level, particle_cfl);

  if (est_dt_particle > 0)
    est_dt = amrex::min<amrex::Real>(est_dt, est_dt_particle);

  if (verbose && ParallelDescriptor::IOProcessor()) {
    if (est_dt_particle > 0) {
      amrex::Print() << "...estdt from particles at level " << level << ": "
                     << est_dt_particle << '\n';
    } else {
      amrex::Print() << "...there are no particles at level " << level << '\n';
    }
  }
}

void
PeleC::readParticleParams()
{
  ParmParse pp("pelec");

  pp.query("do_spray_particles", do_spray_particles);

  ParmParse ppp("particles");

  // Control the verbosity of the Particle class
  ppp.query("v", particle_verbose);

  ppp.get("mass_transfer", particle_mass_tran);
  ppp.get("heat_transfer", particle_heat_tran);
  ppp.get("mom_transfer", particle_mom_tran);
  ppp.query("particle_cfl", particle_cfl);
  if (particle_cfl > 0.5)
    Abort("particle_cfl must be <= 0.5");
  // Number of fuel species in spray droplets
  // Must match the number specified at compile time
  const int nfuel = ppp.countval("fuel_species");
  if (nfuel != SPRAY_FUEL_NUM) {
    Abort("Number of fuel species in input file must match SPRAY_FUEL_NUM");
  }

  sprayFuelNames.assign(nfuel, "");
  sprayCritT.resize(nfuel);
  sprayBoilT.resize(nfuel);
  sprayLatent.resize(nfuel);
  sprayCp.resize(nfuel);
  sprayIndxMap.resize(nfuel);
  std::vector<std::string> fuel_names;
  std::vector<Real> crit_T;
  std::vector<Real> boil_T;
  std::vector<Real> latent;
  std::vector<Real> spraycp;
  ppp.getarr("fuel_species", fuel_names);
  ppp.getarr("fuel_crit_temp", crit_T);
  ppp.getarr("fuel_boil_temp", boil_T);
  ppp.getarr("fuel_latent", latent);
  ppp.getarr("fuel_cp", spraycp);
  for (int i = 0; i < nfuel; ++i) {
    sprayFuelNames[i] = fuel_names[i];
    sprayCritT[i] = crit_T[i];
    sprayBoilT[i] = boil_T[i];
    sprayLatent[i] = latent[i];
    sprayCp[i] = spraycp[i];
  }
  // Must use same reference temperature for all fuels
  // TODO: This means the reference temperature must be the same for all fuel
  // species
  ppp.get("fuel_ref_temp", sprayRefT);

  // Set if particle plot files should be written
  ppp.query("write_particle_plotfiles", write_particle_plotfiles);

  // Set if spray ascii files should be written
  ppp.query("write_spray_ascii_files", write_spray_ascii_files);

  // Used in initData() on startup to read in a file of particles.
  ppp.query("particle_init_file", particle_init_file);

  // Used in initData() on startup to set a uniform particle field
  ppp.query("particle_init_uniform", particle_init_uniform);

  // Used in post_restart() to read in a file of particles.

  // This must be true the first time you try to restart from a checkpoint
  // that was written with USE_PARTICLES=FALSE; i.e. one that doesn't have
  // the particle checkpoint stuff (even if there are no active particles).
  // Otherwise the code will fail when trying to read the checkpointed
  // particles.

  // ppp.query("restart_from_nonparticle_chkfile",
  // restart_from_nonparticle_chkfile);

  // The directory in which to store timestamp files.
  ppp.query("timestamp_dir", timestamp_dir);

  // Only the I/O processor makes the directory if it doesn't already exist.
  if (ParallelDescriptor::IOProcessor())
    if (!amrex::UtilCreateDirectory(timestamp_dir, 0755))
      amrex::CreateDirectoryFailed(timestamp_dir);

  // Force other processors to wait till directory is built.
  ParallelDescriptor::Barrier();
}

void
PeleC::defineParticles()
{
  // There must be at least as many fuel species in the spray as
  // there are species in the fluid
  if (SPRAY_FUEL_NUM > NUM_SPECIES) {
    Abort("Cannot have more spray fuel species than fluid species");
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
          sprayIndxMap[i] = ns;
        }
      }
      if (sprayIndxMap[i] < 0) {
        Print() << "Fuel " << sprayFuelNames[i] << " not found in species list"
                << std::endl;
        Abort();
      }
    }
  } else {
    for (int ns = 0; ns < SPRAY_FUEL_NUM; ++ns) {
      sprayIndxMap[ns] = 0;
    }
  }
}

void
PeleC::setupVirtualParticles()
{
  BL_PROFILE("PeleC::setupVirtualParticles()");
  amrex::Gpu::LaunchSafeGuard lsg(true);
  if (PeleC::theSprayPC() != 0 && !virtual_particles_set) {
    SprayParticleContainer::AoS virts;
    if (level < parent->finestLevel()) {
      ((PeleC*)&parent->getLevel(level + 1))->setupVirtualParticles();
      PeleC::theVirtPC()->CreateVirtualParticles(level + 1, virts);
      PeleC::theVirtPC()->AddParticlesAtLevel(virts, level);

      PeleC::theSprayPC()->CreateVirtualParticles(level + 1, virts);
      PeleC::theVirtPC()->AddParticlesAtLevel(virts, level);
    }
    virtual_particles_set = true;
  }
}

void
PeleC::removeVirtualParticles()
{
  BL_PROFILE("PeleC::removeVirtualParticles()");
  amrex::Gpu::LaunchSafeGuard lsg(true);
  if (VirtPC != 0)
    VirtPC->RemoveParticlesAtLevel(level);
  virtual_particles_set = false;
}

void
PeleC::setupGhostParticles(int ngrow)
{
  BL_PROFILE("PeleC::setupGhostParticles()");
  AMREX_ASSERT(level < parent->finestLevel());
  amrex::Gpu::LaunchSafeGuard lsg(true);
  if (PeleC::theSprayPC() != 0) {
    SprayParticleContainer::AoS ghosts;
    PeleC::theSprayPC()->CreateGhostParticles(level, ngrow, ghosts);
    PeleC::theGhostPC()->AddParticlesAtLevel(ghosts, level + 1, ngrow);
  }
}

void
PeleC::removeGhostParticles()
{
  BL_PROFILE("PeleC::removeGhostParticles()");
  amrex::Gpu::LaunchSafeGuard lsg(true);
  if (GhostPC != 0)
    GhostPC->RemoveParticlesAtLevel(level);
}

// Initialize the particles on the grid at level 0
void
PeleC::initParticles()
{
  BL_PROFILE("PeleC::initParticles()");

  if (level > 0)
    return;

  // Make sure to call RemoveParticlesOnExit() on exit.
  amrex::ExecOnFinalize(RemoveParticlesOnExit);

  if (do_spray_particles) {
    AMREX_ASSERT(theSprayPC() == 0);
    // Whether we need to use ghost and virtual particles
    bool gvParticles = false;
    if (parent->subCycle()) {
      gvParticles = true;
    }

    SprayPC = new SprayParticleContainer(parent, &phys_bc);
    theSprayPC()->SetVerbose(particle_verbose);

    if (gvParticles) {
      VirtPC = new SprayParticleContainer(parent, &phys_bc);
      GhostPC = new SprayParticleContainer(parent, &phys_bc);
    }
    // Pass constant reference data and memory allocations to GPU
    theSprayPC()->buildFuelData(
      sprayCritT, sprayBoilT, sprayCp, sprayLatent, sprayIndxMap, sprayRefT);
    if (gvParticles) {
      theGhostPC()->buildFuelData(
        sprayCritT, sprayBoilT, sprayCp, sprayLatent, sprayIndxMap, sprayRefT);
      theVirtPC()->buildFuelData(
        sprayCritT, sprayBoilT, sprayCp, sprayLatent, sprayIndxMap, sprayRefT);
    }

    if (!particle_init_file.empty()) {
      theSprayPC()->InitFromAsciiFile(particle_init_file, NSR_SPR + NAR_SPR);
    } else if (particle_init_uniform > 0) {
      theSprayPC()->InitParticlesUniform(this, level, particle_init_uniform);
    }
  }
}

void
PeleC::particlePostRestart(const std::string& restart_file, bool is_checkpoint)
{
  amrex::Gpu::setLaunchRegion(false);
  if (level > 0)
    return;

  if (do_spray_particles) {
    AMREX_ASSERT(SprayPC == 0);

    SprayPC = new SprayParticleContainer(parent, &phys_bc);
    theSprayPC()->SetVerbose(particle_verbose);
    theSprayPC()->buildFuelData(
      sprayCritT, sprayBoilT, sprayCp, sprayLatent, sprayIndxMap, sprayRefT);

    if (parent->subCycle()) {
      VirtPC = new SprayParticleContainer(parent, &phys_bc);
      GhostPC = new SprayParticleContainer(parent, &phys_bc);
      theGhostPC()->buildFuelData(
        sprayCritT, sprayBoilT, sprayCp, sprayLatent, sprayIndxMap, sprayRefT);
      theVirtPC()->buildFuelData(
        sprayCritT, sprayBoilT, sprayCp, sprayLatent, sprayIndxMap, sprayRefT);
    }

    // Make sure to call RemoveParticlesOnExit() on exit.
    amrex::ExecOnFinalize(RemoveParticlesOnExit);
    {
      amrex::Gpu::LaunchSafeGuard lsg(true);
      theSprayPC()->Restart(
        parent->theRestartFile(), "particles", is_checkpoint);
      amrex::Gpu::Device::streamSynchronize();
    }
  }
}

// TODO: This has not been checked or updated, use with caution
std::unique_ptr<MultiFab>
PeleC::particleDerive(const std::string& name, Real time, int ngrow)
{
  BL_PROFILE("PeleC::particleDerive()");

  if (theSprayPC() && name == "particle_count") {
    Abort("Should not be called until it is updated");
    std::unique_ptr<MultiFab> derive_dat(new MultiFab(grids, dmap, 1, 0));
    MultiFab temp_dat(grids, dmap, 1, 0);
    temp_dat.setVal(0);
    theSprayPC()->Increment(temp_dat, level);
    MultiFab::Copy(*derive_dat, temp_dat, 0, 0, 1, 0);
    return derive_dat;
  } else if (theSprayPC() && name == "total_particle_count") {
    Abort("Should not be called until it is updated");

    // We want the total particle count at this level or higher.
    std::unique_ptr<MultiFab> derive_dat =
      particleDerive("particle_count", time, ngrow);

    IntVect trr(AMREX_D_DECL(1, 1, 1));

    for (int lev = level + 1; lev <= parent->finestLevel(); lev++) {
      auto ba = parent->boxArray(lev);
      const auto& dm = parent->DistributionMap(lev);
      MultiFab temp_dat(ba, dm, 1, 0);

      trr *= parent->refRatio(lev - 1);

      ba.coarsen(trr);
      MultiFab ctemp_dat(ba, dm, 1, 0);

      temp_dat.setVal(0);
      ctemp_dat.setVal(0);

      theSprayPC()->Increment(temp_dat, lev);

      for (MFIter mfi(temp_dat); mfi.isValid(); ++mfi) {
        const FArrayBox& ffab = temp_dat[mfi];
        FArrayBox& cfab = ctemp_dat[mfi];
        const Box& fbx = ffab.box();

        AMREX_ASSERT(cfab.box() == amrex::coarsen(fbx, trr));

        for (IntVect p = fbx.smallEnd(); p <= fbx.bigEnd(); fbx.next(p)) {
          const Real val = ffab(p);
          if (val > 0)
            cfab(amrex::coarsen(p, trr)) += val;
        }
      }

      temp_dat.clear();

      MultiFab dat(grids, dmap, 1, 0);
      dat.setVal(0);
      dat.copy(ctemp_dat);

      MultiFab::Add(*derive_dat, dat, 0, 0, 1, 0);
    }

    return derive_dat;
  } else {
    return AmrLevel::derive(name, time, ngrow);
  }
}

void
PeleC::particleRedistribute(int lbase, int nGrow, int local, bool init_part)
{
  BL_PROFILE("PeleC::particleRedistribute()");
  int flev = parent->finestLevel();
  if (theSprayPC()) {
    amrex::Gpu::LaunchSafeGuard lsg(true);

    // If we are calling with init_part = true, then we want to force the
    // redistribute without checking whether the grids have changed.
    if (init_part) {
      theSprayPC()->Redistribute(lbase);
      return;
    }

    // These are usually the BoxArray and DMap from the last regridding.
    static Vector<BoxArray> ba;
    static Vector<DistributionMapping> dm;

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
        if (ba[i] != parent->boxArray(i))
          changed = true;

        if (!changed) {
          // Check DistributionMaps have changed during regridding
          if (dm[i] != parent->getLevel(i).get_new_data(0).DistributionMap())
            changed = true;
        }
      }
    }

    if (changed) {
      // We only need to call Redistribute if the BoxArrays or DistMaps have
      // changed.
      // We also only call it for particles >= lbase. This is
      // because if we called redistribute during a subcycle, there may be
      // particles not in the proper position on coarser levels.
      if (verbose && ParallelDescriptor::IOProcessor())
        amrex::Print() << "Calling redistribute because grid has changed "
                       << '\n';
      if (flev == 0) {
        // Do a local redistribute
        theSprayPC()->Redistribute(lbase, -1, nGrow, true);
      } else {
        theSprayPC()->Redistribute(lbase, -1, nGrow, false);
      }

      // Use the new BoxArray and DistMap to define ba and dm for next time.
      for (int i = 0; i <= flev; i++) {
        ba[i] = parent->boxArray(i);
        dm[i] = parent->getLevel(i).get_new_data(0).DistributionMap();
      }
    } else {
      if (verbose && ParallelDescriptor::IOProcessor())
        amrex::Print()
          << "NOT calling redistribute because grid has NOT changed " << '\n';
    }
  }
}

void
PeleC::particleTimestamp(int ngrow)
{
#if 1
  return;
  amrex::Abort("Need to fill in PeleC::TimestampParticles");
#else
  static bool first = true;
  static int imax = -1;
  if (first) {
    first = false;

    ParmParse ppp("particles");

    // have to do it here, not in read_particle_params, because Density, ...,
    // are set after read_particle_params is called.

    int timestamp_density = 1;
    ppp.query("timestamp_density", timestamp_density);
    if (timestamp_density) {
      timestamp_indices.push_back(Density);
      amrex::Print() << "Density = " << Density << std::endl;
    }

    int timestamp_temperature = 0;
    ppp.query("timestamp_temperature", timestamp_temperature);
    if (timestamp_temperature) {
      timestamp_indices.push_back(Temp);
      amrex::Print() << "Temp = " << Temp << std::endl;
    }

    if (!timestamp_indices.empty()) {
      imax = *(amrex::max<amrex::Real> _element(
        timestamp_indices.begin(), timestamp_indices.end()));
    }
  }

  if (theSprayPC() && !timestamp_dir.empty()) {
    std::string basename = timestamp_dir;

    if (basename[basename.length() - 1] != '/')
      basename += '/';

    basename += "Timestamp";

    int finest_level = parent->finestLevel();
    Real time = state[State_Type].curTime();

    for (int lev = level; lev <= finest_level; lev++) {
      if (theSprayPC()->NumberOfParticlesAtLevel(lev) <= 0)
        continue;

      MultiFab& S_new = parent->getLevel(lev).get_new_data(State_Type);

      if (imax >= 0) { // FillPatchIterator will fail otherwise
        int ng = (lev == level) ? ngrow : 1;
        FillPatchIterator fpi(
          parent->getLevel(lev), S_new, ng, time, State_Type, 0, imax + 1);
        const MultiFab& S = fpi.get_mf();
        theSprayPC()->Timestamp(basename, S, lev, time, timestamp_indices);
      } else {
        theSprayPC()->Timestamp(basename, S_new, lev, time, timestamp_indices);
      }
    }
  }
#endif
}

#endif
