#include <iomanip>
#include <vector>
#include <algorithm>
#include <string>

#include "PeleC.H"

using namespace amrex;

#ifdef AMREX_PARTICLES

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
Gpu::HostVector<Real> sprayRho;
Gpu::HostVector<Real> sprayMu;
Gpu::HostVector<int> sprayIndxMap;
amrex::Real sprayRefT = 300.;
amrex::Real parcelSize = 1.;
amrex::Real spraySigma = -1.; // Surface tension
amrex::Real T_wall = -1.;
SprayComps scomps;
bool splash_model = true;

std::string particle_init_file;
int particle_init_function = 1;

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

int PeleC::do_spray_particles = 0;
int PeleC::particle_verbose = 0;
Real PeleC::particle_cfl = 0.5;

int PeleC::write_particle_plotfiles = 1;
int PeleC::write_spray_ascii_files = 1;
int PeleC::particle_mass_tran = 0;
int PeleC::particle_heat_tran = 0;
int PeleC::particle_mom_tran = 0;
Vector<std::string> PeleC::sprayFuelNames;

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
  if (!do_spray_particles) return;
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
    amrex::Abort("particle_cfl must be <= 0.5");
  // Number of fuel species in spray droplets
  // Must match the number specified at compile time
  const int nfuel = ppp.countval("fuel_species");
  if (nfuel != SPRAY_FUEL_NUM)
    amrex::Abort("Number of fuel species in input file must match SPRAY_FUEL_NUM");

  sprayFuelNames.assign(nfuel, "");
  sprayCritT.resize(nfuel);
  sprayBoilT.resize(nfuel);
  sprayCp.resize(nfuel);
  sprayLatent.resize(nfuel);
  sprayRho.resize(nfuel);
  sprayMu.resize(nfuel);
  sprayIndxMap.resize(nfuel);
  std::vector<std::string> fuel_names;
  std::vector<Real> crit_T;
  std::vector<Real> boil_T;
  std::vector<Real> spraycp;
  std::vector<Real> latent;
  std::vector<Real> sprayrho;
  std::vector<Real> mu(nfuel, 0.);
  ppp.getarr("fuel_species", fuel_names);
  ppp.getarr("fuel_crit_temp", crit_T);
  ppp.getarr("fuel_boil_temp", boil_T);
  ppp.getarr("fuel_cp", spraycp);
  ppp.getarr("fuel_latent", latent);
  ppp.getarr("fuel_rho", sprayrho);
  ppp.queryarr("fuel_mu", mu);
  for (int i = 0; i != nfuel; ++i) {
    sprayFuelNames[i] = fuel_names[i];
    sprayCritT[i] = crit_T[i];
    sprayBoilT[i] = boil_T[i];
    sprayCp[i] = spraycp[i];
    sprayLatent[i] = latent[i];
    sprayRho[i] = sprayrho[i];
    sprayMu[i] = mu[i];
  }

  // Set the number of particles per parcel
  ppp.query("parcel_size", parcelSize);
  ppp.query("use_splash_model", splash_model);
  if (splash_model) {
    if (!ppp.contains("fuel_sigma") || !ppp.contains("wall_temp")) {
      Print() << "fuel_sigma or wall_temp must be set for splash model. "
              << "Set use_splash_model = false to turn off splash model" << std::endl;
      Abort();
    }
    // Set the fuel surface tension and contact angle
    ppp.get("fuel_sigma", spraySigma);
    // TODO: Have this retrieved from proper boundary data during runtime
    ppp.get("wall_temp", T_wall);
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

  // Used in initData() on startup to set the particle field using the
  // SprayParticlesInitInsert.cpp problem specific function
  ppp.query("particle_init_function", particle_init_function);

  if (verbose && ParallelDescriptor::IOProcessor()) {
    amrex::Print() << "Spray fuel species " << sprayFuelNames[0];
    for (int i = 1; i != SPRAY_FUEL_NUM; ++i)
      amrex::Print() << ", " << sprayFuelNames[i];
    amrex::Print() << std::endl;
    amrex::Print() << "Number of particles per parcel " << parcelSize << std::endl;
  }
  // Force other processors to wait till directory is built.
  ParallelDescriptor::Barrier();
}

void
PeleC::defineParticles()
{
  // There must be at least as many fuel species in the spray as
  // there are species in the fluid
  if (SPRAY_FUEL_NUM > NUM_SPECIES) {
    amrex::Abort("Cannot have more spray fuel species than fluid species");
  }
#ifdef PELEC_EOS_FUEGO
  for (int i = 0; i != SPRAY_FUEL_NUM; ++i) {
    for (int ns = 0; ns != NUM_SPECIES; ++ns) {
      std::string gas_spec = spec_names[ns];
      if (gas_spec == sprayFuelNames[i]) {
        sprayIndxMap[i] = ns;
      }
    }
    if (sprayIndxMap[i] < 0) {
      amrex::Print() << "Fuel " << sprayFuelNames[i] << " not found in species list"
                     << std::endl;
      amrex::Abort();
    }
  }
#else
  for (int ns = 0; ns != SPRAY_FUEL_NUM; ++ns) {
    sprayIndxMap[ns] = 0;
  }
#endif
  scomps.heat_tran = PeleC::particle_heat_tran;
  scomps.mass_tran = PeleC::particle_mass_tran;
  scomps.mom_tran = PeleC::particle_mom_tran;
  scomps.rhoIndx = PeleC::Density;
  scomps.momIndx = PeleC::Xmom;
  scomps.engIndx = PeleC::Eden;
  scomps.utempIndx = PeleC::Temp;
  scomps.specIndx = PeleC::FirstSpec;
}

void
PeleC::setupVirtualParticles()
{
  BL_PROFILE("PeleC::setupVirtualParticles()");
  amrex::Gpu::LaunchSafeGuard lsg(true);
  if (PeleC::theSprayPC() != 0 && !virtual_particles_set) {
    if (level < parent->finestLevel()) {
#ifdef USE_SPRAY_SOA
      SprayParticleContainer::ParticleTileType virts;
#else
      SprayParticleContainer::AoS virts;
#endif
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
#ifdef USE_SPRAY_SOA
    SprayParticleContainer::ParticleTileType ghosts;
#else
    SprayParticleContainer::AoS ghosts;
#endif
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

// Create new particle data
void
PeleC::createParticleData()
{
  SprayPC = new SprayParticleContainer(parent, &phys_bc);
  theSprayPC()->SetVerbose(particle_verbose);
  VirtPC = new SprayParticleContainer(parent, &phys_bc);
  GhostPC = new SprayParticleContainer(parent, &phys_bc);
  // Pass constant reference data and memory allocations to GPU
  theSprayPC()->buildFuelData(sprayCritT, sprayBoilT, sprayCp, sprayLatent,
                              sprayRho, sprayMu, sprayIndxMap, scomps,
                              parcelSize, sprayRefT, spraySigma, T_wall);
  theGhostPC()->buildFuelData(sprayCritT, sprayBoilT, sprayCp, sprayLatent,
                              sprayRho, sprayMu, sprayIndxMap, scomps,
                              parcelSize, sprayRefT, spraySigma, T_wall);
  theVirtPC()->buildFuelData(sprayCritT, sprayBoilT, sprayCp, sprayLatent,
                             sprayRho, sprayMu, sprayIndxMap, scomps,
                             parcelSize, sprayRefT, spraySigma, T_wall);
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
    createParticleData();

    if (!particle_init_file.empty()) {
      theSprayPC()->InitFromAsciiFile(particle_init_file, NSR_SPR + NAR_SPR);
    } else if (particle_init_function > 0) {
      theSprayPC()->InitSprayParticles();
    }
  }
}

void
PeleC::particlePostRestart(const std::string& restart_file, bool is_checkpoint)
{
  if (level > 0)
    return;

  if (do_spray_particles) {
    AMREX_ASSERT(SprayPC == 0);
    createParticleData();

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

void
PeleC::particleMKD (const Real       time,
                    const Real       dt,
                    const int        ghost_width,
                    const int        spray_n_grow,
                    const int        tmp_src_width,
                    const int        where_width,
                    amrex::MultiFab& tmp_spray_source)
{
  amrex::Gpu::LaunchSafeGuard lsg(true);
  // Setup ghost particles for use in finer levels. Note that ghost
  // particles that will be used by this level have already been created,
  // the particles being set here are only used by finer levels.
  int finest_level = parent->finestLevel();

  // Check if I need to insert new particles
  int nstep = parent->levelSteps(0);

  BL_PROFILE_VAR("SprayParticles::injectParticles()", INJECT_SPRAY);
  bool injectParts = theSprayPC()->
    injectParticles(time, dt, nstep, level, finest_level);
  bool insertParts = theSprayPC()->
    insertParticles(time, dt, nstep, level, finest_level);
  // Only redistribute if we injected or inserted particles
  if (injectParts || insertParts)
    theSprayPC()->Redistribute(level);
  BL_PROFILE_VAR_STOP(INJECT_SPRAY);

  // Setup the virtual particles that represent particles on finer levels
  if (level < finest_level)
    setupVirtualParticles();

  // Make a copy of the particles on this level into ghost particles
  // for the finer level
  if (level < finest_level)
    setupGhostParticles(ghost_width);

  // Advance the particle velocities to the half-time and the positions to
  // the new time
  if (particle_verbose)
    amrex::Print()
      << "moveKickDrift ... updating particle positions and velocity\n";

  // We will make a temporary copy of the source term array inside
  // moveKickDrift and we are only going to use the spray force out to one ghost cell 
  // so we need only define spray_force_old with one ghost cell

  AMREX_ASSERT(old_sources[spray_src]->nGrow() >= 1);

  // Do the valid particles themselves
  theSprayPC()->moveKickDrift(
    Sborder, tmp_spray_source,
    level, dt, time,
    false, // not virtual particles
    false, // not ghost particles
    spray_n_grow,
    tmp_src_width,
    true, // Move the particles
    where_width);

  // Only need the coarsest virtual particles here.
  if (level < finest_level)
    theVirtPC()->moveKickDrift(
      Sborder, tmp_spray_source,
      level, dt, time, true, false,
      spray_n_grow, tmp_src_width, true, where_width);

  // Miiiight need all Ghosts
  if (theGhostPC() != 0 && level != 0)
    theGhostPC()->moveKickDrift(
      Sborder, tmp_spray_source,
      level, dt, time, false, true,
      spray_n_grow, tmp_src_width, true, where_width);
  // Must call transfer source after moveKick and moveKickDrift
  // on all particle types
  theSprayPC()->transferSource(
    tmp_src_width, level, tmp_spray_source, *old_sources[spray_src]);
}

void
PeleC::particleMK (const Real       time,
                   const Real       dt,
                   const int        spray_n_grow,
                   const int        tmp_src_width,
                   const int        amr_iteration,
                   const int        amr_ncycle,
                   amrex::MultiFab& tmp_spray_source)
{
  theSprayPC()->moveKick(
    Sborder, tmp_spray_source,
    level, dt, time, false, false,
    spray_n_grow, tmp_src_width);

  if (level < parent->finestLevel())
    theVirtPC()->moveKick(
      Sborder, tmp_spray_source,
      level, dt, time, true, false,
      spray_n_grow, tmp_src_width);

  if (theGhostPC() != 0 && level != 0)
    theGhostPC()->moveKick(
      Sborder, tmp_spray_source,
      level, dt, time, false, true,
      spray_n_grow, tmp_src_width);
  theSprayPC()->transferSource(
    tmp_src_width, level, tmp_spray_source, *new_sources[spray_src]);
}

// TODO: This has not been checked or updated, use with caution
std::unique_ptr<MultiFab>
PeleC::particleDerive(const std::string& name, Real time, int ngrow)
{
  BL_PROFILE("PeleC::particleDerive()");

  if (theSprayPC() && name == "particle_count") {
    amrex::Abort("Should not be called until it is updated");
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
PeleC::particle_redistribute(int lbase, bool init_part)
{
  if (!do_spray_particles) return;
  BL_PROFILE("PeleC::particle_redistribute()");
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
