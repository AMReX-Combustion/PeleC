#include <cmath>

#include "mechanism.h"

#include "PeleC.H"
#include "IndexDefines.H"

amrex::Real
PeleC::advance(
  amrex::Real time, amrex::Real dt, int amr_iteration, int amr_ncycle)
{
  /** the main driver for a single level implementing the time advance.

         @param time the current simulation time
         @param dt the timestep to advance (e.g., go from time to time + dt)
         @param amr_iteration where we are in the current AMR subcycle.  Each
                         level will take a number of steps to reach the
                         final time of the coarser level below it.  This
                         counter starts at 1
         @param amr_ncycle  the number of subcycles at this level
  */

  BL_PROFILE("PeleC::advance()");

  int finest_level = parent->finestLevel();

  if (level < finest_level && do_reflux) {
    getFluxReg(level + 1).reset();

    if (!amrex::DefaultGeometry().IsCartesian()) {
      amrex::Abort("Flux registers not r-z compatible yet");
      getPresReg(level + 1).reset();
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
#ifdef PELEC_USE_REACTIONS
    if (!(i == Reactions_Type && do_react)) {
#endif
      state[i].allocOldData();
      state[i].swapTimeLevels(dt);
#ifdef PELEC_USE_REACTIONS
    }
#endif
  }

  if (do_mol_load_balance || do_react_load_balance) {
    get_new_data(Work_Estimate_Type).setVal(0.0);
  }

  // get old and new state
  amrex::MultiFab& S_old = get_old_data(State_Type);
  amrex::MultiFab& S_new = get_new_data(State_Type);

  // define sourceterm
  amrex::MultiFab molSrc(grids, dmap, NVAR, 0, amrex::MFInfo(), Factory());

  amrex::MultiFab molSrc_old, molSrc_new;
  if (mol_iters > 1) {
    molSrc_old.define(grids, dmap, NVAR, 0, amrex::MFInfo(), Factory());
    molSrc_new.define(grids, dmap, NVAR, 0, amrex::MFInfo(), Factory());
  }

#ifdef PELEC_USE_REACTIONS
  amrex::MultiFab& I_R = get_new_data(Reactions_Type);
#endif

#ifdef PELEC_USE_EB
  set_body_state(S_old);
  set_body_state(S_new);
#endif

  // Compute S^{n} = MOLRhs(U^{n})
  if (verbose) {
    amrex::Print() << "... Computing MOL source term at t^{n} " << std::endl;
  }
  int nGrow_Sborder = nGrowTr;
#ifdef AMREX_PARTICLES
  int ghost_width = 0;
  int where_width = 0;
  int spray_n_grow = 0;
  int tmp_src_width = 0;

  if (do_spray_particles) {
    int finest_level = parent->finestLevel();
    int finer_ref = 0;
    if (level < finest_level) finer_ref = parent->MaxRefRatio(level);
    theSprayPC()->setSprayGridInfo(
      level, finest_level, amr_ncycle, amr_iteration, finer_ref,
      ghost_width, where_width, spray_n_grow, tmp_src_width);
    nGrow_Sborder = std::max(nGrow_Sborder, spray_n_grow);
  }
  if (Sborder.nGrow() < nGrow_Sborder) {
    std::string abortStr =
      "Sborder has " + std::to_string(Sborder.nGrow()) +
      " ghost cells, but needs " + std::to_string(nGrow_Sborder);
    amrex::Abort(abortStr);
  }
#endif
  FillPatch(*this, Sborder, nGrow_Sborder, time, State_Type, 0, NVAR);
  amrex::Real flux_factor = 0;
  getMOLSrcTerm(Sborder, molSrc, time, dt, flux_factor);

#ifdef AMREX_PARTICLES
  // We must make a temporary spray source term to ensure number of ghost
  // cells are correct
  amrex::MultiFab tmp_spray_source(
    grids, dmap, NVAR, tmp_src_width, amrex::MFInfo(), Factory());
  tmp_spray_source.setVal(0.);
  if (do_spray_particles) {
    old_sources[spray_src]->setVal(0.);
    particleMKD(
      time, dt, ghost_width, spray_n_grow,
      tmp_src_width, where_width, tmp_spray_source);
    amrex::MultiFab::Saxpy(molSrc, 1.0, *old_sources[spray_src], 0, 0, NVAR, 0);
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
        src_list[n], time, dt, amr_iteration, amr_ncycle, 0, 0);

      // add sources to molsrc
      amrex::MultiFab::Saxpy(
        molSrc, 1.0, *old_sources[src_list[n]], 0, 0, NVAR, 0);
    }
  }

  if (mol_iters > 1)
    amrex::MultiFab::Copy(molSrc_old, molSrc, 0, 0, NVAR, 0);

  // U^* = U^n + dt*S^n
  amrex::MultiFab::LinComb(S_new, 1.0, Sborder, 0, dt, molSrc, 0, 0, NVAR, 0);

#ifdef PELEC_USE_REACTIONS
  // U^{n+1,*} = U^n + dt*S^n + dt*I_R
  if (do_react == 1) {
    amrex::MultiFab::Saxpy(S_new, dt, I_R, 0, FirstSpec, NUM_SPECIES, 0);
    amrex::MultiFab::Saxpy(S_new, dt, I_R, NUM_SPECIES, Eden, 1, 0);
  }
#endif

  computeTemp(S_new, 0);

  // Compute S^{n+1} = MOLRhs(U^{n+1,*})
  if (verbose) {
    amrex::Print() << "... Computing MOL source term at t^{n+1} " << std::endl;
  }
  FillPatch(*this, Sborder, nGrow_Sborder, time + dt, State_Type, 0, NVAR);
  flux_factor = mol_iters > 1 ? 0 : 1;
  getMOLSrcTerm(Sborder, molSrc, time, dt, flux_factor);

#ifdef AMREX_PARTICLES
  if (do_spray_particles) {
    new_sources[spray_src]->setVal(0.);
    particleMK(
      time + dt, dt, spray_n_grow, tmp_src_width, amr_iteration, amr_ncycle, tmp_spray_source);
    amrex::MultiFab::Saxpy(molSrc, 1.0, *new_sources[spray_src], 0, 0, NVAR, 0);
  }
#endif

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

#ifdef PELEC_USE_REACTIONS
  if (do_react == 1) {
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
#endif

  computeTemp(S_new, 0);

#ifdef PELEC_USE_REACTIONS
  if (do_react == 1) {
    for (int mol_iter = 2; mol_iter <= mol_iters; ++mol_iter) {
      if (verbose) {
        amrex::Print() << "... Re-computing MOL source term at t^{n+1} (iter = "
                       << mol_iter << " of " << mol_iters << ")" << std::endl;
      }
      FillPatch(*this, Sborder, nGrow_Sborder, time + dt, State_Type, 0, NVAR);
      flux_factor = mol_iter == mol_iters ? 1 : 0;
      getMOLSrcTerm(Sborder, molSrc_new, time, dt, flux_factor);

      // F_{AD} = (1/2)(S_old + S_new)
      amrex::MultiFab::LinComb(
        molSrc, 0.5, molSrc_old, 0, 0.5, molSrc_new, 0, 0, NVAR, 0);

      // Compute I_R and U^{n+1} = U^n + dt*(F_{AD} + I_R)
      react_state(time, dt, false, &molSrc);

      computeTemp(S_new, 0);
    }
  }
#endif

#ifdef PELEC_USE_EB
  set_body_state(S_new);
#endif

  return dt;
}

amrex::Real
PeleC::do_sdc_advance(
  amrex::Real time, amrex::Real dt, int amr_iteration, int amr_ncycle)
{
  BL_PROFILE("PeleC::do_sdc_advance()");

  amrex::Real dt_new = dt;

  /** This routine will advance the old state data (called S_old here)
      to the new time, for a single level.  The new data is called
      S_new here.  The update includes reactions (if we are not doing
      SDC), hydro, and the source terms.
  */

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
  /** This routine will advance the old state data (called S_old here)
      to the new time, for a single level.  The new data is called
      S_new here.  The update includes hydro, and the source terms.
  */

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
    nGrow_Sborder = NUM_GROW + nGrowF;
  } else if (do_diffuse) {
    fill_Sborder = true;
    nGrow_Sborder = NUM_GROW;
  }
#ifdef AMREX_PARTICLES
  int ghost_width = 0;
  int where_width = 0;
  int spray_n_grow = 0;
  int tmp_src_width = 0;

  if (do_spray_particles) {
    int finest_level = parent->finestLevel();
    int finer_ref = 0;
    if (level < finest_level) finer_ref = parent->MaxRefRatio(level);
    theSprayPC()->setSprayGridInfo(
      level, finest_level, amr_ncycle, amr_iteration, finer_ref,
      ghost_width, where_width, spray_n_grow, tmp_src_width);
    fill_Sborder = true;
    nGrow_Sborder = std::max(nGrow_Sborder, spray_n_grow);
  }
  if (fill_Sborder && Sborder.nGrow() < nGrow_Sborder) {
    std::string abortStr =
      "Sborder has " + std::to_string(Sborder.nGrow()) +
      " ghost cells, but needs " + std::to_string(nGrow_Sborder);
    amrex::Abort(abortStr);
  }

#endif

  if (fill_Sborder) {
    FillPatch(*this, Sborder, nGrow_Sborder, time, State_Type, 0, NVAR);
  }

#ifdef AMREX_PARTICLES
  // We must make a temporary spray source term to ensure number of ghost
  // cells are correct
  amrex::MultiFab tmp_spray_source(
    grids, dmap, NVAR, tmp_src_width, amrex::MFInfo(), Factory());
  tmp_spray_source.setVal(0.);
#endif
  if (sub_iteration == 0) {
#ifdef AMREX_PARTICLES
    //
    // Compute drag terms from particles at old positions, move particles to new
    // positions
    //  based on old-time velocity field
    //
    // TODO: Maybe move this mess into construct_old_source?
    if (do_spray_particles) {
      old_sources[spray_src]->setVal(0.);
      particleMKD(
        time, dt, ghost_width, spray_n_grow,
        tmp_src_width, where_width, tmp_spray_source);
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
    // cells, and we
    //  may want to reuse what we fill-patched for hydro
    if (do_diffuse) {
      if (verbose) {
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
    if (verbose) {
      amrex::Print() << "... Computing diffusion terms at t^(n+1,"
                     << sub_iteration + 1 << ")" << std::endl;
    }
    int nGrowDiff = nGrowTr;
#ifdef AMREX_PARTICLES
    if (do_spray_particles && level > 0) {
      int maxref = parent->MaxRefRatio(level - 1);
      if (maxref > 2) nGrowDiff = nGrow_Sborder;
    }
#endif
    FillPatch(*this, Sborder, nGrowDiff, time + dt, State_Type, 0, NVAR);
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
  if (do_spray_particles && sub_iteration == sub_ncycle - 1) {
    new_sources[spray_src]->setVal(0.);
    // Advance the particle velocities by dt/2 to the new time.
    if (particle_verbose)
      amrex::Print() << "moveKick ... updating velocity only\n";

    if (!do_diffuse) { // Else, this was already done above.  No need to redo
      FillPatch(*this, Sborder, nGrow_Sborder, time + dt, State_Type, 0, NVAR);
    }
    particleMK(
      time + dt, dt, spray_n_grow, tmp_src_width, amr_iteration, amr_ncycle, tmp_spray_source);
  }
#endif

#ifdef PELEC_USE_REACTIONS
  // Update I_R and rebuild S_new accordingly
  if (do_react == 1) {
    react_state(time, dt);
  } else {
    construct_Snew(S_new, S_old, dt);
    get_new_data(Reactions_Type).setVal(0);
  }
#else
  construct_Snew(S_new, S_old, dt);
#endif

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

#ifdef PELEC_USE_REACTIONS
  if (do_react == 1) {
    amrex::MultiFab& I_R = get_new_data(Reactions_Type);
    amrex::MultiFab::Saxpy(S_new, dt, I_R, 0, FirstSpec, NUM_SPECIES, 0);
    amrex::MultiFab::Saxpy(S_new, dt, I_R, NUM_SPECIES, Eden, 1, 0);
  }
#endif
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
    for (int i = 0; i < n_lost; i++) {
      material_lost_through_boundary_temp[i] = 0.0;
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

#ifdef PELEC_USE_REACTIONS
  if (do_react == 1) {
    // Initialize I_R with value from previous time step
    amrex::MultiFab::Copy(
      get_new_data(Reactions_Type), get_old_data(Reactions_Type), 0, 0,
      get_new_data(Reactions_Type).nComp(),
      get_new_data(Reactions_Type).nGrow());
  }
#endif
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
