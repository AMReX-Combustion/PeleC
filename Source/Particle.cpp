
#include "PeleC.H"

using namespace amrex;

#ifdef AMREX_PARTICLES

namespace {
bool virtual_particles_set = false;

// Containers for the real "active" Particles
SprayParticleContainer* SprayPC = nullptr;

// Container for temporary, virtual Particles
SprayParticleContainer* VirtPC = nullptr;

// Container for temporary, ghost Particles
SprayParticleContainer* GhostPC = nullptr;

SprayData sprayData;
amrex::Real spray_ref_temp = 300.;
amrex::Real parcel_size = 1.;
amrex::Real spray_sigma = -1.; // Surface tension
amrex::Real wall_temp = -1.;
// Indices for spray source MultiFab
int sprayRhoSrcIndx = 0;
int sprayMomSrcIndx = 1;
int sprayEngSrcIndx = 1 + AMREX_SPACEDIM;
int spraySpecSrcIndx = 2 + AMREX_SPACEDIM;
SprayComps scomps;
bool splash_model = true;

std::string particle_init_file;
int particle_init_function = 1;

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
} // namespace

int PeleC::do_spray_particles = 1;
int PeleC::particle_verbose = 0;
Real PeleC::particle_cfl = 0.5;

int PeleC::write_spray_ascii_files = 0;
// momentum + density + fuel species + energy
int PeleC::num_spray_src = AMREX_SPACEDIM + 2 + SPRAY_FUEL_NUM;
int PeleC::particle_mass_tran = 1;
int PeleC::particle_mom_tran = 1;
Vector<std::string> PeleC::spray_fuel_names;

void
getPSatCoef(
  Real* psat_coef, ParmParse& ppp, std::string fuel_name, const int spf)
{
  std::string psat_read = fuel_name + "_psat";
  std::vector<Real> inp_coef(4);
  ppp.getarr(psat_read.c_str(), inp_coef);
  for (int i = 0; i < 4; ++i) {
    psat_coef[4 * spf + i] = inp_coef[i];
  }
}

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
  if (do_spray_particles == 0)
    return;
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
  ppp.get("mom_transfer", particle_mom_tran);
  ppp.query("cfl", particle_cfl);
  if (particle_cfl > 0.5)
    amrex::Abort("particles.cfl must be <= 0.5");
  // Number of fuel species in spray droplets
  // Must match the number specified at compile time
  const int nfuel = ppp.countval("fuel_species");
  if (nfuel != SPRAY_FUEL_NUM)
    amrex::Abort(
      "Number of fuel species in input file must match SPRAY_FUEL_NUM");

  std::vector<std::string> fuel_names;
  std::vector<Real> crit_T;
  std::vector<Real> boil_T;
  std::vector<Real> spraycp;
  std::vector<Real> latent;
  std::vector<Real> sprayrho;
  std::vector<Real> mu(nfuel, -1.);
  std::vector<Real> lambda(nfuel, -1.);
  if (particle_mass_tran) {
    spray_fuel_names.assign(nfuel, "");
    ppp.getarr("fuel_species", fuel_names);
    ppp.getarr("fuel_crit_temp", crit_T);
    ppp.getarr("fuel_boil_temp", boil_T);
    ppp.getarr("fuel_cp", spraycp);
    ppp.getarr("fuel_latent", latent);
    ppp.getarr("fuel_rho", sprayrho);
    ppp.queryarr("fuel_mu", mu);
    ppp.queryarr("fuel_lambda", lambda);
    for (int i = 0; i < nfuel; ++i) {
      spray_fuel_names[i] = fuel_names[i];
      sprayData.critT[i] = crit_T[i];
      sprayData.boilT[i] = boil_T[i];
      sprayData.cp[i] = spraycp[i];
      sprayData.latent[i] = latent[i];
      sprayData.ref_latent[i] = latent[i];
      sprayData.rho[i] = sprayrho[i];
      sprayData.mu[i] = mu[i];
      sprayData.lambda[i] = lambda[i];
      getPSatCoef(sprayData.psat_coef.data(), ppp, fuel_names[i], i);
    }
  }

  // Set the number of particles per parcel
  ppp.query("parcel_size", parcel_size);
  ppp.query("use_splash_model", splash_model);
  if (splash_model) {
    if (
      !ppp.contains("fuel_sigma") || !ppp.contains("wall_temp") || mu[0] < 0. ||
      lambda[0] < 0.) {
      Print() << "particles.fuel_sigma, wall_temp, fuel_mu, and fuel_lambda "
              << "must be set for splash model. Set use_splash_model = false "
              << "to turn off splash model" << std::endl;
      Abort();
    }
    // Set the fuel surface tension and contact angle
    ppp.get("fuel_sigma", spray_sigma);
    // TODO: Have this retrieved from proper boundary data during runtime
    ppp.get("wall_temp", wall_temp);
  }

  // Must use same reference temperature for all fuels
  // TODO: This means the reference temperature must be the same for all fuel
  // species
  ppp.get("fuel_ref_temp", spray_ref_temp);

  // Set if spray ascii files should be written
  ppp.query("write_spray_ascii_files", write_spray_ascii_files);

  // Used in initData() on startup to read in a file of particles.
  ppp.query("init_file", particle_init_file);

  // Used in initData() on startup to set the particle field using the
  // SprayParticlesInitInsert.cpp problem specific function
  ppp.query("init_function", particle_init_function);

  // Set the data for the liquid fuel and spray droplets
  sprayData.num_ppp = parcel_size;
  sprayData.ref_T = spray_ref_temp;
  sprayData.sigma = spray_sigma;

  if (verbose && ParallelDescriptor::IOProcessor()) {
    amrex::Print() << "Spray fuel species " << spray_fuel_names[0];
    for (int i = 1; i < SPRAY_FUEL_NUM; ++i)
      amrex::Print() << ", " << spray_fuel_names[i];
    amrex::Print() << std::endl;
    amrex::Print() << "Number of particles per parcel " << parcel_size
                   << std::endl;
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
  for (int i = 0; i < SPRAY_FUEL_NUM; ++i) {
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
      std::string gas_spec = spec_names[ns];
      if (gas_spec == spray_fuel_names[i]) {
        sprayData.indx[i] = ns;
      }
    }
    if (sprayData.indx[i] < 0) {
      amrex::Print() << "Fuel " << spray_fuel_names[i]
                     << " not found in species list" << std::endl;
      amrex::Abort();
    }
  }
#else
  for (int ns = 0; ns < SPRAY_FUEL_NUM; ++ns)
    sprayData.indx[ns] = 0;
#endif
  amrex::Vector<Real> fuelEnth(NUM_SPECIES);
  auto eos = pele::physics::PhysicsType::eos();
  eos.T2Hi(sprayData.ref_T, fuelEnth.data());
  for (int ns = 0; ns < SPRAY_FUEL_NUM; ++ns) {
    const int fspec = sprayData.indx[ns];
    sprayData.latent[ns] -= fuelEnth[fspec];
  }
  scomps.mass_tran = PeleC::particle_mass_tran;
  scomps.mom_tran = PeleC::particle_mom_tran;
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

void
PeleC::setupVirtualParticles()
{
  BL_PROFILE("PeleC::setupVirtualParticles()");
  if (PeleC::theSprayPC() != nullptr && !virtual_particles_set) {
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
  if (VirtPC != nullptr)
    VirtPC->RemoveParticlesAtLevel(level);
  virtual_particles_set = false;
}

void
PeleC::setupGhostParticles(int ngrow)
{
  BL_PROFILE("PeleC::setupGhostParticles()");
  AMREX_ASSERT(level < parent->finestLevel());
  if (PeleC::theSprayPC() != nullptr) {
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
  if (GhostPC != nullptr)
    GhostPC->RemoveParticlesAtLevel(level);
}

// Create new particle data
void
PeleC::createParticleData()
{
  SprayPC = new SprayParticleContainer(
    parent, &phys_bc, sprayData, scomps, parcel_size, wall_temp);
  theSprayPC()->SetVerbose(particle_verbose);
  VirtPC = new SprayParticleContainer(
    parent, &phys_bc, sprayData, scomps, parcel_size, wall_temp);
  GhostPC = new SprayParticleContainer(
    parent, &phys_bc, sprayData, scomps, parcel_size, wall_temp);
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

  if (do_spray_particles == 1) {
    AMREX_ASSERT(theSprayPC() == 0);
    createParticleData();

    if (!particle_init_file.empty()) {
      theSprayPC()->InitFromAsciiFile(particle_init_file, NSR_SPR + NAR_SPR);
    } else if (particle_init_function > 0) {
      const ProbParmHost* lprobparm = prob_parm_host;
      const ProbParmDevice* lprobparm_d = h_prob_parm_device;
      theSprayPC()->InitSprayParticles(*lprobparm, *lprobparm_d);
    } else {
      Abort("Must initialize spray particles with particles.init_function or "
            "particles.init_file");
    }
  }
}

void
PeleC::particlePostRestart(const std::string& restart_file, bool is_checkpoint)
{
  if (level > 0)
    return;

  if (do_spray_particles == 1) {
    AMREX_ASSERT(SprayPC == 0);
    createParticleData();

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
  const Real time,
  const Real dt,
  const int ghost_width,
  const int spray_n_grow,
  const int tmp_src_width,
  const int where_width,
  amrex::MultiFab& tmp_spray_source)
{
  // Setup ghost particles for use in finer levels. Note that ghost
  // particles that will be used by this level have already been created,
  // the particles being set here are only used by finer levels.
  int finest_level = parent->finestLevel();

  // Check if I need to insert new particles
  int nstep = parent->levelSteps(0);

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
  // moveKickDrift and we are only going to use the spray force out to one ghost
  // cell so we need only define spray_force_old with one ghost cell

  AMREX_ASSERT(old_sources[spray_src]->nGrow() >= 1);

  // Do the valid particles themselves
  theSprayPC()->moveKickDrift(
    Sborder, tmp_spray_source, level, dt, time,
    false, // not virtual particles
    false, // not ghost particles
    spray_n_grow, tmp_src_width,
    true, // Move the particles
    where_width);

  // Only need the coarsest virtual particles here.
  if (level < finest_level)
    theVirtPC()->moveKickDrift(
      Sborder, tmp_spray_source, level, dt, time, true, false, spray_n_grow,
      tmp_src_width, true, where_width);

  // Miiiight need all Ghosts
  if (theGhostPC() != nullptr && level != 0)
    theGhostPC()->moveKickDrift(
      Sborder, tmp_spray_source, level, dt, time, false, true, spray_n_grow,
      tmp_src_width, true, where_width);
  // Must call transfer source after moveKick and moveKickDrift
  // on all particle types
  theSprayPC()->transferSource(
    tmp_src_width, level, tmp_spray_source, *old_sources[spray_src]);
}

void
PeleC::particleMK(
  const Real time,
  const Real dt,
  const int spray_n_grow,
  const int tmp_src_width,
  const int amr_iteration,
  const int amr_ncycle,
  amrex::MultiFab& tmp_spray_source)
{
  theSprayPC()->moveKick(
    Sborder, tmp_spray_source, level, dt, time, false, false, spray_n_grow,
    tmp_src_width);

  if (level < parent->finestLevel())
    theVirtPC()->moveKick(
      Sborder, tmp_spray_source, level, dt, time, true, false, spray_n_grow,
      tmp_src_width);

  if (theGhostPC() != nullptr && level != 0)
    theGhostPC()->moveKick(
      Sborder, tmp_spray_source, level, dt, time, false, true, spray_n_grow,
      tmp_src_width);
  theSprayPC()->transferSource(
    tmp_src_width, level, tmp_spray_source, *new_sources[spray_src]);
}

void
PeleC::particle_redistribute(int lbase, bool init_part)
{
  if (do_spray_particles == 0)
    return;
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

#endif
