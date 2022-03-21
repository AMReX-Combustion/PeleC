#include "mechanism.H"

#include "PeleC.H"
#include "IndexDefines.H"

amrex::Real
PeleC::advance(
  amrex::Real time, amrex::Real dt, int amr_iteration, int amr_ncycle)
{
  // The main driver for a single level implementing the time advance.
  //        @param time the current simulation time
  //        @param dt the timestep to advance (e.g., go from time to time + dt)
  //        @param amr_iteration where we are in the current AMR subcycle.  Each
  //                        level will take a number of steps to reach the
  //                        final time of the coarser level below it.  This
  //                        counter starts at 1
  //        @param amr_ncycle  the number of subcycles at this level

  BL_PROFILE("PeleC::advance()");

  int finest_level = parent->finestLevel();

  if (level < finest_level && do_reflux) {
    getFluxReg(level + 1).reset();

    if (!amrex::DefaultGeometry().IsCartesian()) {
      amrex::Abort("Flux registers not r-z compatible yet");
    }
  }

  amrex::Real dt_new;
  if (do_mol) {
    dt_new = do_mol_advance(time, dt, amr_iteration, amr_ncycle);
  } else {
    dt_new = do_sdc_advance(time, dt, amr_iteration, amr_ncycle);
  }

  return dt_new;
}

amrex::Real
PeleC::do_mol_advance(
  amrex::Real time, amrex::Real dt, int amr_iteration, int amr_ncycle)
{
  BL_PROFILE("PeleC::do_mol_advance()");

  // Check that we are not asking to advance stuff we don't know to
  // if (src_list.size() > 0) amrex::Abort("Have not integrated other sources
  // into MOL advance yet");

  for (int i = 0; i < num_state_type; ++i) {
    if (!(i == Reactions_Type && do_react)) {
      state[i].allocOldData();
      state[i].swapTimeLevels(dt);
    }
  }

  if (do_mol_load_balance || do_react_load_balance) {
    get_new_data(Work_Estimate_Type).setVal(0.0);
  }

  amrex::MultiFab& S_old = get_old_data(State_Type);
  amrex::MultiFab& S_new = get_new_data(State_Type);

  amrex::MultiFab molSrc(grids, dmap, NVAR, 0, amrex::MFInfo(), Factory());

  amrex::MultiFab molSrc_old;
  amrex::MultiFab molSrc_new;
  if (mol_iters > 1) {
    molSrc_old.define(grids, dmap, NVAR, 0, amrex::MFInfo(), Factory());
    molSrc_new.define(grids, dmap, NVAR, 0, amrex::MFInfo(), Factory());
  }

  if (!do_react) {
    get_new_data(Reactions_Type).setVal(0.0);
  }
  const amrex::MultiFab& I_R = get_new_data(Reactions_Type);

  set_body_state(S_old);
  set_body_state(S_new);

  // Compute S^{n} = MOLRhs(U^{n})
  if (verbose != 0) {
    amrex::Print() << "... Computing MOL source term at t^{n} " << std::endl;
  }
  FillPatch(*this, Sborder, numGrow() + nGrowF, time, State_Type, 0, NVAR);
  amrex::Real flux_factor = 0;
  getMOLSrcTerm(Sborder, molSrc, time, dt, flux_factor);

  // Build other (neither spray nor diffusion) sources at t_old
  for (int n = 0; n < src_list.size(); ++n) {
    if (
      src_list[n] != diff_src
#ifdef AMREX_PARTICLES
      && src_list[n] != spray_src
#endif
    ) {
      construct_old_source(
        src_list[n], time, dt, amr_iteration, amr_ncycle, 0, 0);

      // add sources to molsrc
      amrex::MultiFab::Saxpy(
        molSrc, 1.0, *old_sources[src_list[n]], 0, 0, NVAR, 0);
    }
  }

  if (mol_iters > 1) {
    amrex::MultiFab::Copy(molSrc_old, molSrc, 0, 0, NVAR, 0);
  }

  // U^* = U^n + dt*S^n
  amrex::MultiFab::LinComb(S_new, 1.0, Sborder, 0, dt, molSrc, 0, 0, NVAR, 0);

  // U^{n+1,*} = U^n + dt*S^n + dt*I_R
  if (do_react) {
    amrex::MultiFab::Saxpy(S_new, dt, I_R, 0, FirstSpec, NUM_SPECIES, 0);
    amrex::MultiFab::Saxpy(S_new, dt, I_R, NUM_SPECIES, Eden, 1, 0);
  }

  computeTemp(S_new, 0);

  // Compute S^{n+1} = MOLRhs(U^{n+1,*})
  if (verbose != 0) {
    amrex::Print() << "... Computing MOL source term at t^{n+1} " << std::endl;
  }
  FillPatch(*this, Sborder, numGrow() + nGrowF, time + dt, State_Type, 0, NVAR);
  flux_factor = mol_iters > 1 ? 0 : 1;
  getMOLSrcTerm(Sborder, molSrc, time, dt, flux_factor);

  // Build other (neither spray nor diffusion) sources at t_new
  for (int n = 0; n < src_list.size(); ++n) {
    if (
      src_list[n] != diff_src
#ifdef AMREX_PARTICLES
      && src_list[n] != spray_src
#endif
    ) {
      construct_new_source(
        src_list[n], time + dt, dt, amr_iteration, amr_ncycle, 0, 0);

      // add sources to molsrc
      amrex::MultiFab::Saxpy(
        molSrc, 1.0, *new_sources[src_list[n]], 0, 0, NVAR, 0);
    }
  }

  // U^{n+1.**} = 0.5*(U^n + U^{n+1,*}) + 0.5*dt*S^{n+1} = U^n + 0.5*dt*S^n +
  // 0.5*dt*S^{n+1} + 0.5*dt*I_R
  amrex::MultiFab::LinComb(S_new, 0.5, Sborder, 0, 0.5, S_old, 0, 0, NVAR, 0);
  amrex::MultiFab::Saxpy(
    S_new, 0.5 * dt, molSrc, 0, 0, NVAR,
    0); //  NOTE: If I_R=0, we are done and U_new is the final new-time state

  if (do_react) {
    amrex::MultiFab::Saxpy(S_new, 0.5 * dt, I_R, 0, FirstSpec, NUM_SPECIES, 0);
    amrex::MultiFab::Saxpy(S_new, 0.5 * dt, I_R, NUM_SPECIES, Eden, 1, 0);

    // F_{AD} = (1/dt)(U^{n+1,**} - U^n) - I_R = 0.5*(S^{n}+S^{n+1}(which is a
    // guess!))
    amrex::MultiFab::LinComb(
      molSrc, 1.0 / dt, S_new, 0, -1.0 / dt, S_old, 0, 0, NVAR, 0);
    amrex::MultiFab::Subtract(molSrc, I_R, 0, FirstSpec, NUM_SPECIES, 0);
    amrex::MultiFab::Subtract(molSrc, I_R, NUM_SPECIES, Eden, 1, 0);

    // Compute I_R and U^{n+1} = U^n + dt*(F_{AD} + I_R)
    react_state(time, dt, false, &molSrc);
  }

  computeTemp(S_new, 0);

  if (do_react) {
    for (int mol_iter = 2; mol_iter <= mol_iters; ++mol_iter) {
      if (verbose != 0) {
        amrex::Print() << "... Re-computing MOL source term at t^{n+1} (iter = "
                       << mol_iter << " of " << mol_iters << ")" << std::endl;
      }
      FillPatch(
        *this, Sborder, numGrow() + nGrowF, time + dt, State_Type, 0, NVAR);
      flux_factor = mol_iter == mol_iters ? 1 : 0;
      getMOLSrcTerm(Sborder, molSrc_new, time, dt, flux_factor);

      // F_{AD} = (1/2)(molSrc_old + molSrc_new)
      amrex::MultiFab::LinComb(
        molSrc, 0.5, molSrc_old, 0, 0.5, molSrc_new, 0, 0, NVAR, 0);

      // Compute I_R and U^{n+1} = U^n + dt*(F_{AD} + I_R)
      react_state(time, dt, false, &molSrc);

      computeTemp(S_new, 0);
    }
  }

  set_body_state(S_new);

  return dt;
}

#ifdef AMREX_PARTICLES
void
PeleC::setSprayGridInfo(
  const int amr_iteration,
  const int amr_ncycle,
  int& ghost_width,
  int& where_width,
  int& spray_n_grow,
  int& tmp_src_width)
{
  // TODO: Re-evaluate these numbers and include the particle cfl into the
  // calcuation A particle in cell (i) can affect cell values in (i-1) to (i+1)
  int stencil_deposition_width = 1;

  // A particle in cell (i) may need information from cell values in (i-1) to
  // (i+1)
  //   to update its position (typically via interpolation of the acceleration
  //   from the grid)
  int stencil_interpolation_width = 1;

  // A particle that starts in cell (i + amr_ncycle) can reach
  //   cell (i) in amr_ncycle number of steps .. after "amr_iteration" steps
  //   the particle has to be within (i + amr_ncycle+1-amr_iteration) to reach
  //   cell (i) in the remaining (amr_ncycle-amr_iteration) steps

  // *** ghost_width ***  is used
  //   *) to set how many cells are used to hold ghost particles i.e copies of
  //   particles that live on (level-1) can affect the grid over all of the
  //   amr_ncycle steps. We define ghost cells at the coarser level to cover all
  //   iterations so we can't reduce this number as amr_iteration increases.

  ghost_width = 0;
  if (parent->subCycle() && parent->maxLevel() > 0)
    ghost_width += amr_ncycle + stencil_deposition_width;

  // *** where_width ***  is used
  //   *) to set how many cells the Where call in moveKickDrift tests = max of
  //     {ghost_width + (1-amr_iteration) - 1}:
  //      the minus 1 arises because this occurs *after* the move} and
  //     {amr_iteration}:
  //     the number of cells out that a cell initially in the fine grid may
  //     have moved and we don't want to just lose it (we will redistribute it
  //     when we're done}

  where_width = amrex::max<amrex::Real>(
    ghost_width + (1 - amr_iteration) - 1, amr_iteration);

  // *** spray_n_grow *** is used
  //   *) to determine how many ghost cells we need to fill in the MultiFab from
  //      which the particle interpolates its acceleration

  spray_n_grow = ghost_width + stencil_interpolation_width;

  // *** tmp_src_width ***  is used
  //   *) to set how many ghost cells are needed in the tmp_src_ptr MultiFab
  //   that we
  //      define inside moveKickDrift and moveKick.   This needs to be big
  //      enough to hold the contribution from all the particles within
  //      ghost_width so that we don't have to test on whether the particles are
  //      trying to write out of bounds

  tmp_src_width = ghost_width + stencil_deposition_width;
}
#endif

amrex::Real
PeleC::do_sdc_advance(
  amrex::Real time, amrex::Real dt, int amr_iteration, int amr_ncycle)
{
  BL_PROFILE("PeleC::do_sdc_advance()");

  amrex::Real dt_new = dt;

  // This routine will advance the old state data (called S_old here)
  // to the new time, for a single level.  The new data is called
  // S_new here.  The update includes reactions (if we are not doing
  // SDC), hydro, and the source terms.

  initialize_sdc_advance(time, dt, amr_iteration, amr_ncycle);

  if (do_react_load_balance) {
    get_new_data(Work_Estimate_Type).setVal(0.0);
  }

  for (int sdc_iter = 0; sdc_iter < sdc_iters; ++sdc_iter) {
    if (sdc_iters > 1) {
      amrex::Print() << "SDC iteration " << sdc_iter + 1 << " of " << sdc_iters
                     << ".\n";
    }

    dt_new = do_sdc_iteration(
      time, dt, amr_iteration, amr_ncycle, sdc_iter, sdc_iters);
  }

  finalize_sdc_advance(time, dt, amr_iteration, amr_ncycle);

  return dt_new;
}

amrex::Real
PeleC::do_sdc_iteration(
  amrex::Real time,
  amrex::Real dt,
  int amr_iteration,
  int amr_ncycle,
  int sub_iteration,
  int sub_ncycle)
{
  // This routine will advance the old state data (called S_old here)
  // to the new time, for a single level.  The new data is called
  // S_new here.  The update includes hydro, and the source terms.

  BL_PROFILE("PeleC::do_sdc_iteration()");

  const amrex::MultiFab& S_old = get_old_data(State_Type);
  amrex::MultiFab& S_new = get_new_data(State_Type);

  initialize_sdc_iteration(
    time, dt, amr_iteration, amr_ncycle, sub_iteration, sub_ncycle);

  // Create Sborder if hydro or diffuse, with the appropriate number of grow
  // cells
  int nGrow_Sborder = 0;
  bool fill_Sborder = false;

  if (do_hydro) {
    fill_Sborder = true;
    nGrow_Sborder = numGrow() + nGrowF;
  } else if (do_diffuse) {
    fill_Sborder = true;
    nGrow_Sborder = numGrow();
  }
#ifdef AMREX_PARTICLES
  bool use_ghost_parts = false; // Use ghost particles
  bool use_virt_parts = false;  // Use virtual particles
  if (parent->finestLevel() > 0 && level < parent->finestLevel()) {
    use_ghost_parts = true;
    use_virt_parts = true;
  }
  int ghost_width = 0;
  int where_width = 0;
  int spray_n_grow = 0;
  int tmp_src_width = 0;

  if (do_spray_particles) {
    setSprayGridInfo(
      amr_iteration, amr_ncycle, ghost_width, where_width, spray_n_grow,
      tmp_src_width);
    fill_Sborder = true;
    nGrow_Sborder = std::max(nGrow_Sborder, spray_n_grow);
  }
  if (fill_Sborder && Sborder.nGrow() < nGrow_Sborder) {
    Print() << "PeleC::do_sdc_iteration thinks Sborder needs " << nGrow_Sborder
            << " grow cells, but Sborder defined with only " << Sborder.nGrow()
            << std::endl;
    Abort();
  }

#endif

  if (fill_Sborder) {
    FillPatch(*this, Sborder, nGrow_Sborder, time, State_Type, 0, NVAR);
  }

  if (sub_iteration == 0) {
#ifdef AMREX_PARTICLES

    // Compute drag terms from particles at old positions, move particles to new
    // positions  based on old-time velocity field
    // TODO: Maybe move this mess into construct_old_source?
    if (do_spray_particles) {
      amrex::Gpu::LaunchSafeGuard lsg(true);

      // Setup ghost particles for use in finer levels. Note that ghost
      // particles that will be used by this level have already been created,
      // the particles being set here are only used by finer levels.
      int finest_level = parent->finestLevel();

      // Setup the virtual particles that particles on finer levels
      if (level < finest_level && use_virt_parts)
        setupVirtualParticles();

      // Check if I need to insert new particles
      amrex::Real cur_time = state[State_Type].curTime();
      int nstep = parent->levelSteps(0);

      bool injectParts = false;
      bool insertParts = false;
      if (level == finest_level)
        injectParts = theSprayPC()->injectParticles(cur_time, nstep, level);
      if (level == finest_level)
        insertParts = theSprayPC()->insertParticles(cur_time, nstep, level);

      // Only redistribute if we injected or inserted particles
      if (injectParts || insertParts) {
        // TODO: Determine the number of ghost cells needed here
        int nGrow = 1;
        particleRedistribute(level, nGrow, 0);
      }

      // Make a copy of the particles on this level into ghost particles
      // for the finer level
      if (use_ghost_parts)
        setupGhostParticles(ghost_width);

      // Advance the particle velocities to the half-time and the positions to
      // the new time
      if (particle_verbose)
        amrex::Print()
          << "moveKickDrift ... updating particle positions and velocity\n";

      // We will make a temporary copy of the source term array inside
      // moveKickDrift and we are only going to use the spray force out to one
      // ghost cell so we need only define spray_force_old with one ghost cell

      AMREX_ASSERT(old_sources[spray_src]->nGrow() >= 1);
      old_sources[spray_src]->setVal(0.);

      // Do the valid particles themselves
      theSprayPC()->moveKickDrift(
        Sborder, *old_sources[spray_src], level, dt, cur_time,
        false, // not virtual particles
        false, // not ghost particles
        tmp_src_width,
        true, // Move the particles
        where_width);

      // Only need the coarsest virtual particles here.
      if (level < finest_level && use_virt_parts)
        theVirtPC()->moveKickDrift(
          Sborder, *old_sources[spray_src], level, dt, cur_time, true, false,
          tmp_src_width, true, where_width);

      // Might need all Ghosts
      if (use_ghost_parts && theGhostPC() != 0)
        theGhostPC()->moveKickDrift(
          Sborder, *old_sources[spray_src], level, dt, cur_time, false, true,
          tmp_src_width, true, where_width);
    }
#endif

    // Build other (neither spray nor diffusion) sources at t_old
    for (int n = 0; n < src_list.size(); ++n) {
      if (
        src_list[n] != diff_src
#ifdef AMREX_PARTICLES
        && src_list[n] != spray_src
#endif
      ) {
        construct_old_source(
          src_list[n], time, dt, amr_iteration, amr_ncycle, sub_iteration,
          sub_ncycle);
      }
    }

    // Get diffusion source separate from other sources, since it requires grow
    // cells, and we may want to reuse what we fill-patched for hydro
    if (do_diffuse) {
      if (verbose != 0) {
        amrex::Print() << "... Computing diffusion terms at t^(n)" << std::endl;
      }
      AMREX_ASSERT(
        !do_mol); // Currently this combo only managed through MOL integrator
      amrex::Real flux_factor_old = 0.5;

      getMOLSrcTerm(Sborder, *old_sources[diff_src], time, dt, flux_factor_old);
    }

    // Initialize sources at t_new by copying from t_old
    for (int n = 0; n < src_list.size(); ++n) {
      amrex::MultiFab::Copy(
        *new_sources[src_list[n]], *old_sources[src_list[n]], 0, 0, NVAR, 0);
    }
  }

  // Construct hydro source, will use old and current iterate of new sources.
  if (do_hydro) {
    construct_hydro_source(
      Sborder, time, dt, amr_iteration, amr_ncycle, sub_iteration, sub_ncycle);
  }

  // Construct S_new with current iterate of all sources
  construct_Snew(S_new, S_old, dt);

  int ng_src = 0;
  computeTemp(S_new, ng_src);

  // Now update t_new sources (diffusion separate because it requires a fill
  // patch)
  if (do_diffuse) {
    if (verbose != 0) {
      amrex::Print() << "... Computing diffusion terms at t^(n+1,"
                     << sub_iteration + 1 << ")" << std::endl;
    }
    FillPatch(*this, Sborder, numGrow(), time + dt, State_Type, 0, NVAR);
    amrex::Real flux_factor_new = sub_iteration == sub_ncycle - 1 ? 0.5 : 0;
    getMOLSrcTerm(Sborder, *new_sources[diff_src], time, dt, flux_factor_new);
  }

  // Build other (neither spray nor diffusion) sources at t_new
  for (int n = 0; n < src_list.size(); ++n) {
    if (
      src_list[n] != diff_src
#ifdef AMREX_PARTICLES
      && src_list[n] != spray_src
#endif
    ) {
      construct_new_source(
        src_list[n], time + dt, dt, amr_iteration, amr_ncycle, sub_iteration,
        sub_ncycle);
    }
  }

#ifdef AMREX_PARTICLES
  if (do_spray_particles) {
    // Advance the particle velocities by dt/2 to the new time.
    if (particle_verbose)
      amrex::Print() << "moveKick ... updating velocity only\n";

    if (!do_diffuse) { // Else, this was already done above.  No need to redo
      FillPatch(*this, Sborder, nGrow_Sborder, time + dt, State_Type, 0, NVAR);
    }

    new_sources[spray_src]->setVal(0.);

    theSprayPC()->moveKick(
      Sborder, *new_sources[spray_src], level, dt, time + dt, false, false,
      tmp_src_width);

    // Virtual particles will be recreated, so we need not kick them.
    // TODO: Is this true with SDC iterations??

    // Ghost particles need to be kicked except during the final iteration.
    if (amr_iteration != amr_ncycle && use_ghost_parts && theGhostPC() != 0)
      theGhostPC()->moveKick(
        Sborder, *new_sources[spray_src], level, dt, time + dt, false, true,
        tmp_src_width);
  }
#endif

  // Update I_R and rebuild S_new accordingly
  if (do_react) {
    react_state(time, dt);
  } else {
    construct_Snew(S_new, S_old, dt);
    get_new_data(Reactions_Type).setVal(0);
  }

  computeTemp(S_new, ng_src);

  finalize_sdc_iteration(
    time, dt, amr_iteration, amr_ncycle, sub_iteration, sub_ncycle);

  return dt;
}

void
PeleC::construct_Snew(
  amrex::MultiFab& S_new, const amrex::MultiFab& S_old, amrex::Real dt)
{
  int ng = 0;

  amrex::MultiFab::Copy(S_new, S_old, 0, 0, NVAR, ng);
  for (int n = 0; n < src_list.size(); ++n) {
    amrex::MultiFab::Saxpy(
      S_new, 0.5 * dt, *new_sources[src_list[n]], 0, 0, NVAR, ng);
    amrex::MultiFab::Saxpy(
      S_new, 0.5 * dt, *old_sources[src_list[n]], 0, 0, NVAR, ng);
  }
  if (do_hydro) {
    amrex::MultiFab::Saxpy(S_new, dt, hydro_source, 0, 0, NVAR, ng);
  }

  if (do_react) {
    amrex::MultiFab& I_R = get_new_data(Reactions_Type);
    amrex::MultiFab::Saxpy(S_new, dt, I_R, 0, FirstSpec, NUM_SPECIES, 0);
    amrex::MultiFab::Saxpy(S_new, dt, I_R, NUM_SPECIES, Eden, 1, 0);
  }
}

void
PeleC::initialize_sdc_iteration(
  amrex::Real /*time*/,
  amrex::Real /*dt*/,
  int /*amr_iteration*/,
  int /*amr_ncycle*/,
  int /*sdc_iteration*/,
  int /*sdc_ncycle*/)
{
  BL_PROFILE("PeleC::initialize_sdc_iteration()");

  // Reset the change from density resets
  frac_change = 1;

  // Reset the grid loss tracking.
  if (track_grid_losses) {
    for (amrex::Real& i : material_lost_through_boundary_temp) {
      i = 0.0;
    }
  }
}

void
PeleC::finalize_sdc_iteration(
  amrex::Real /*time*/,
  amrex::Real /*dt*/,
  int /*amr_iteration*/,
  int /*amr_ncycle*/,
  int /*sdc_iteration*/,
  int /*sdc_ncycle*/)
{
  BL_PROFILE("PeleC::finalize_sdc_iteration()");
}

void
PeleC::initialize_sdc_advance(
  amrex::Real /*time*/,
  amrex::Real dt,
  int /*amr_iteration*/,
  int /*amr_ncycle*/)
{
  BL_PROFILE("PeleC::initialize_sdc_advance()");

  for (int i = 0; i < num_state_type; ++i) {
    state[i].allocOldData();
    state[i].swapTimeLevels(dt);
  }

  if (do_react) {
    // Initialize I_R with value from previous time step
    amrex::MultiFab::Copy(
      get_new_data(Reactions_Type), get_old_data(Reactions_Type), 0, 0,
      get_new_data(Reactions_Type).nComp(),
      get_new_data(Reactions_Type).nGrow());
  }
}

void
PeleC::finalize_sdc_advance(
  amrex::Real /*time*/,
  amrex::Real /*dt*/,
  int /*amr_iteration*/,
  int /*amr_ncycle*/)
{
  BL_PROFILE("PeleC::finalize_sdc_advance()");

  // Add the material lost in this timestep to the cumulative losses.
  if (track_grid_losses) {
    amrex::ParallelDescriptor::ReduceRealSum(
      material_lost_through_boundary_temp, n_lost);

    for (int i = 0; i < n_lost; i++) {
      material_lost_through_boundary_cumulative[i] +=
        material_lost_through_boundary_temp[i];
    }
  }

  // amrex::Real cur_time = state[State_Type].curTime();
}
