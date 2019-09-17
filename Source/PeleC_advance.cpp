#include <PeleC.H>
#include <PeleC_F.H>

#include <cmath>

using namespace amrex;

using std::string;

Real
PeleC::advance (Real time,
		Real dt,
		int  amr_iteration,
		int  amr_ncycle)
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

  if (level < finest_level && do_reflux)
  {
    getFluxReg(level+1).reset();

    if (!DefaultGeometry().IsCartesian()) {
      amrex::Abort("Flux registers not r-z compatible yet");
      getPresReg(level+1).reset();
    }
  }

  Real dt_new = dt;
  if (do_mol_AD)
  {
    dt_new = do_mol_advance(time, dt, amr_iteration, amr_ncycle);
  }
  else
  {
    dt_new = do_sdc_advance(time, dt, amr_iteration, amr_ncycle);
  }

  return dt_new;
}

Real
PeleC::do_mol_advance(Real time,
                      Real dt,
                      int  amr_iteration,
                      int  amr_ncycle)
{
  BL_PROFILE("PeleC::do_mol_advance()");

  // Check that we are not asking to advance stuff we don't know to
  //if (src_list.size() > 0) amrex::Abort("Have not integrated other sources into MOL advance yet");

  for (int i = 0; i < num_state_type; ++i) {
    bool skip = false;
#ifdef REACTIONS
    skip = i == Reactions_Type && do_react;
#endif
    if (!skip) {
      state[i].allocOldData();
      state[i].swapTimeLevels(dt);
    }
  }

  if (do_mol_load_balance || do_react_load_balance)
  {
    get_new_data(Work_Estimate_Type).setVal(0.0);
  }

  MultiFab& U_old = get_old_data(State_Type);
  MultiFab& U_new = get_new_data(State_Type);
  MultiFab S(grids,dmap,NUM_STATE,0,MFInfo(),Factory());

  MultiFab S_old, S_new;
  if (mol_iters > 1) {
    S_old.define(grids,dmap,NUM_STATE,0,MFInfo(),Factory());
    S_new.define(grids,dmap,NUM_STATE,0,MFInfo(),Factory());
  }

#ifdef REACTIONS
  MultiFab& I_R = get_new_data(Reactions_Type);
#endif

#ifdef PELE_USE_EB
  set_body_state(U_old);
  set_body_state(U_new);
#endif


  // Compute S^{n} = MOLRhs(U^{n})
  if (verbose) { amrex::Print() << "... Computing MOL source term at t^{n} " << std::endl; }
  FillPatch(*this, Sborder, nGrowTr, time, State_Type, 0, NUM_STATE);
  Real flux_factor = 0;
  getMOLSrcTerm(Sborder, S, time, dt, flux_factor);

    // Build other (neither spray nor diffusion) sources at t_old
    for (int n = 0; n < src_list.size(); ++n)
    {
      if (src_list[n] != diff_src
#ifdef AMREX_PARTICLES
          && src_list[n] != spray_src
#endif
        )
      {
	construct_old_source(src_list[n], time, dt, amr_iteration, amr_ncycle, 0, 0);
        MultiFab::Saxpy(S, 1.0, *old_sources[src_list[n]], 0, 0, NUM_STATE, 0);
      }
    }

  if (mol_iters > 1) MultiFab::Copy(S_old,S,0,0,NUM_STATE,0);

  // U^* = U^n + dt*S^n
  MultiFab::LinComb(U_new, 1.0, Sborder, 0, dt, S, 0, 0, NUM_STATE, 0);

#ifdef REACTIONS
    // U^{n+1,*} = U^n + dt*S^n + dt*I_R)
  if (do_react == 1) {
    MultiFab::Saxpy(U_new, dt, I_R, 0,      FirstSpec, NumSpec, 0);
    MultiFab::Saxpy(U_new, dt, I_R, NumSpec,Eden,      1,       0);
  }
#endif

  computeTemp(U_new,0);


  // Compute S^{n+1} = MOLRhs(U^{n+1,*})
  if (verbose) { amrex::Print() << "... Computing MOL source term at t^{n+1} " << std::endl; }
  FillPatch(*this, Sborder, nGrowTr, time+dt, State_Type, 0, NUM_STATE);
  flux_factor = mol_iters > 1 ?  0 : 1;
  getMOLSrcTerm(Sborder, S, time, dt, flux_factor);

  // Build other (neither spray nor diffusion) sources at t_new
  for (int n = 0; n < src_list.size(); ++n)
  {
    if (src_list[n] != diff_src
#ifdef AMREX_PARTICLES
      && src_list[n] != spray_src
#endif
      )
    {
      construct_new_source(src_list[n], time + dt, dt, amr_iteration, amr_ncycle, 0, 0);
      MultiFab::Saxpy(S, 1.0, *new_sources[src_list[n]], 0, 0, NUM_STATE, 0);
    }
  }

  // U^{n+1.**} = 0.5*(U^n + U^{n+1,*}) + 0.5*dt*S^{n+1} = U^n + 0.5*dt*S^n + 0.5*dt*S^{n+1} + 0.5*dt*I_R
  MultiFab::LinComb(U_new, 0.5, Sborder, 0, 0.5, U_old, 0, 0, NUM_STATE, 0);
  MultiFab::Saxpy(U_new, 0.5*dt, S, 0, 0, NUM_STATE, 0);   //  NOTE: If I_R=0, we are done and U_new is the final new-time state

#ifdef REACTIONS
  if (do_react == 1) {
    MultiFab::Saxpy(U_new, 0.5*dt, I_R, 0,       FirstSpec, NumSpec, 0);
    MultiFab::Saxpy(U_new, 0.5*dt, I_R, NumSpec, Eden,      1,       0);

    // F_{AD} = (1/dt)(U^{n+1,**} - U^n) - I_R
    MultiFab::LinComb(S, 1.0/dt, U_new, 0, -1.0/dt, U_old, 0, 0, NUM_STATE, 0);
    MultiFab::Subtract(S, I_R, 0,      FirstSpec, NumSpec, 0);
    MultiFab::Subtract(S, I_R, NumSpec,Eden,      1,       0);

    // Compute I_R and U^{n+1} = U^n + dt*(F_{AD} + I_R)
    react_state(time, dt, false, &S);  // false = not react_init
  }
#endif

  computeTemp(U_new,0);



#ifdef REACTIONS
  if (do_react == 1)
  {
    for (int mol_iter = 2; mol_iter<=mol_iters; ++mol_iter)
    {
      if (verbose) { amrex::Print() << "... Re-computing MOL source term at t^{n+1} (iter = " << mol_iter << " of " << mol_iters << ")" << std::endl; }
      FillPatch(*this, Sborder, nGrowTr, time + dt, State_Type, 0, NUM_STATE);
      flux_factor = mol_iter==mol_iters  ?  1  : 0;
      getMOLSrcTerm(Sborder, S_new, time, dt, flux_factor);

      // F_{AD} = (1/2)(S_old + S_new)
      MultiFab::LinComb(S, 0.5, S_old, 0, 0.5, S_new, 0, 0, NUM_STATE, 0);

      // Compute I_R and U^{n+1} = U^n + dt*(F_{AD} + I_R)
      react_state(time, dt, false, &S);  // false = not react_init
      
      computeTemp(U_new,0);
    }
  }
#endif

#ifdef PELE_USE_EB
  set_body_state(U_new);
#endif

  return dt;
}

#ifdef AMREX_PARTICLES
void
PeleC::set_spray_grid_info(int amr_iteration,
                           int amr_ncycle,
                           int ghost_width,
                           int where_width,
                           int spray_n_grow,
                           int tmp_src_width)
{
  // A particle in cell (i) can affect cell values in (i-1) to (i+1)
  int stencil_deposition_width = 1;

  // A particle in cell (i) may need information from cell values in (i-1) to (i+1)
  //   to update its position (typically via interpolation of the acceleration from the grid)
  int stencil_interpolation_width = 1;

  // A particle that starts in cell (i + amr_ncycle) can reach
  //   cell (i) in amr_ncycle number of steps .. after "amr_iteration" steps
  //   the particle has to be within (i + amr_ncycle+1-amr_iteration) to reach cell (i)
  //   in the remaining (amr_ncycle-amr_iteration) steps

  // *** ghost_width ***  is used
  //   *) to set how many cells are used to hold ghost particles i.e copies of particles
  //      that live on (level-1) can affect the grid over all of the amr_ncycle steps.
  //      We define ghost cells at the coarser level to cover all iterations so
  //      we can't reduce this number as amr_iteration increases.

  ghost_width = amr_ncycle + stencil_deposition_width;

  // *** where_width ***  is used
  //   *) to set how many cells the Where call in moveKickDrift tests = max of
  //     {ghost_width + (1-amr_iteration) - 1}:
  //      the minus 1 arises because this occurs *after* the move} and
  //     {amr_iteration}:
  //     the number of cells out that a cell initially in the fine grid may
  //     have moved and we don't want to just lose it (we will redistribute it when we're done}

  where_width =  std::max(ghost_width + (1-amr_iteration) - 1, amr_iteration);

  // *** spray_n_grow *** is used
  //   *) to determine how many ghost cells we need to fill in the MultiFab from
  //      which the particle interpolates its acceleration

  //spray_n_grow = ghost_width + (1-amr_iteration) + stencil_interpolation_width;
  spray_n_grow = ghost_width + (1-amr_iteration) + (amr_iteration-1) +
                     stencil_interpolation_width ;


  // *** tmp_src_width ***  is used
  //   *) to set how many ghost cells are needed in the tmp_src_ptr MultiFab that we
  //      define inside moveKickDrift and moveKick.   This needs to be big enough
  //      to hold the contribution from all the particles within ghost_width so that
  //      we don't have to test on whether the particles are trying to write out of bounds

  tmp_src_width = ghost_width + stencil_deposition_width;

//   std::cout << "src_width="<<tmp_src_width<< "spray_n_grow="<<spray_n_grow 
//             << "where_width="<<where_width<< "ghost_width="<<ghost_width
//             << "stencil_deposition_width="<<stencil_deposition_width
//             << "stencil_interpolation_width="<<stencil_interpolation_width; 
}
#endif

Real
PeleC::do_sdc_advance(Real time,
                      Real dt,
                      int  amr_iteration,
                      int  amr_ncycle)
{
  BL_PROFILE("PeleC::do_sdc_advance()");

  Real dt_new = dt;

  /** This routine will advance the old state data (called S_old here)
      to the new time, for a single level.  The new data is called
      S_new here.  The update includes reactions (if we are not doing
      SDC), hydro, and the source terms.
  */

  initialize_sdc_advance(time, dt, amr_iteration, amr_ncycle);

  if (do_react_load_balance) {
    get_new_data(Work_Estimate_Type).setVal(0.0);
  }

  for (int sdc_iter = 0; sdc_iter < sdc_iters; ++sdc_iter)
  {
    if (sdc_iters > 1)
    {
      amrex::Print() << "SDC iteration " << sdc_iter + 1 << " of " << sdc_iters << ".\n";
    }

    dt_new = do_sdc_iteration(time, dt, amr_iteration, amr_ncycle, sdc_iter, sdc_iters);
  }

  finalize_sdc_advance(time, dt, amr_iteration, amr_ncycle);

  return dt_new;
}

Real
PeleC::do_sdc_iteration (Real time,
                         Real dt,
                         int  amr_iteration,
                         int  amr_ncycle,
                         int  sub_iteration,
                         int  sub_ncycle)
{
  /** This routine will advance the old state data (called S_old here)
      to the new time, for a single level.  The new data is called
      S_new here.  The update includes hydro, and the source terms.
  */

  BL_PROFILE("PeleC::do_sdc_iteration()");

  MultiFab& S_old = get_old_data(State_Type);
  MultiFab& S_new = get_new_data(State_Type);

  initialize_sdc_iteration(time, dt, amr_iteration, amr_ncycle, 
                           sub_iteration, sub_ncycle);

  // Create Sborder if hydro or diffuse, with the appropriate number of grow cells
  int nGrow_Sborder = 0;
  bool fill_Sborder = false;

  if (do_hydro)
  {
    fill_Sborder = true;
    nGrow_Sborder = NUM_GROW+nGrowF;
  }
  else if (do_diffuse)
  {
    fill_Sborder = true;
    nGrow_Sborder = NUM_GROW;
  }
#ifdef AMREX_PARTICLES
  fill_Sborder = true;
  int ghost_width = 0;
  int where_width = 0;
  int spray_n_grow = 0;
  int tmp_src_width = 0;
  set_spray_grid_info(amr_iteration,amr_ncycle,ghost_width,where_width,spray_n_grow,tmp_src_width);

  if (do_spray_particles)
  {
    nGrow_Sborder = std::max(nGrow_Sborder, spray_n_grow);
  }
  if (fill_Sborder &&  Sborder.nGrow() < nGrow_Sborder) {
    Print() << "PeleC::do_sdc_iteration thinks Sborder needs " << nGrow_Sborder
            << " grow cells, but Sborder defined with only " << Sborder.nGrow() << std::endl;
    Abort();
  }

#endif

  if (fill_Sborder)
  {
    FillPatch(*this, Sborder, nGrow_Sborder, time, State_Type, 0, NUM_STATE);
  }
  
  if (sub_iteration==0)
  {
#ifdef AMREX_PARTICLES
    //
    // Compute drag terms from particles at old positions, move particles to new positions
    //  based on old-time velocity field
    //
    // TODO: Maybe move this mess into construct_old_source?
    if (do_spray_particles)
    {
      //
      // Setup the virtual particles that represent finer level particles
      //
      setup_virtual_particles();

      //
      // Setup ghost particles for use in finer levels. Note that ghost particles
      // that will be used by this level have already been created, the
      // particles being set here are only used by finer levels.
      //
      int finest_level = parent->finestLevel();

      //  
      // Check if I need to insert new particles
      //
        Real cur_time = state[State_Type].curTime();
        int nstep = parent->levelSteps(0);

//      if (level == finest_level) 
//        theSprayPC()->injectParticles(cur_time,nstep,level);
        if (level == finest_level) 
          theSprayPC()->insertParticles(cur_time,nstep,level);

      particle_redistribute(level,false);

      if (level < finest_level)
        setup_ghost_particles(ghost_width);

      // Advance the particle velocities to the half-time and the positions to the new time
      if (particle_verbose)
        amrex::Print() << "moveKickDrift ... updating particle positions and velocity\n";

      // We will make a temporary copy of the source term array inside moveKickDrift
      //    and we are only going to use the spray force out to one ghost cell
      //    so we need only define spray_force_old with one ghost cell

      BL_ASSERT(old_sources[spray_src]->nGrow() >= 1);
      old_sources[spray_src]->setVal(0.);

      // Do the valid particles themselves
      theSprayPC()->moveKickDrift(Sborder,*old_sources[spray_src], level, dt, tmp_src_width, where_width);

      // Only need the coarsest virtual particles here.
      if (level < finest_level)
        theVirtPC()->moveKickDrift(Sborder,*old_sources[spray_src], level, dt, tmp_src_width, where_width);

      // Miiiight need all Ghosts
      if (theGhostPC() != 0)
         theGhostPC()->moveKickDrift(Sborder,*old_sources[spray_src], level, dt, tmp_src_width, where_width);
    }
#endif

    // Build other (neither spray nor diffusion) sources at t_old
    for (int n = 0; n < src_list.size(); ++n)
    {
      if (src_list[n] != diff_src
#ifdef AMREX_PARTICLES
          && src_list[n] != spray_src
#endif
        )
      {
	construct_old_source(src_list[n], time, dt, amr_iteration, amr_ncycle,
			     sub_iteration, sub_ncycle);
      }
    }

    // Get diffusion source separate from other sources, since it requires grow cells, and we
    //  may want to reuse what we fill-patched for hydro
    if (do_diffuse)
    {
      if (verbose) {
	amrex::Print() << "... Computing diffusion terms at t^(n)" << std::endl;
      }
      BL_ASSERT(!do_mol_AD); // Currently this combo only managed through MOL integrator
      Real flux_factor_old = 0.5;
      getMOLSrcTerm(Sborder,*old_sources[diff_src],time,dt,flux_factor_old);
    }

    // Initialize sources at t_new by copying from t_old
    for (int n = 0; n < src_list.size(); ++n)
    {
      MultiFab::Copy(*new_sources[src_list[n]],
		     *old_sources[src_list[n]],0,0,NUM_STATE,0);
    }
  }

  // Construct hydro source, will use old and current iterate of new sources.
  if (do_hydro)
  {
    construct_hydro_source(Sborder, time, dt, amr_iteration, amr_ncycle, sub_iteration, sub_ncycle);
  }

  // Construct S_new with current iterate of all sources
  construct_Snew(S_new, S_old, dt);

  int ng_src = 0;
  computeTemp(S_new, ng_src);

  // Now update t_new sources (diffusion separate because it requires a fill patch)
  if (do_diffuse)
  {
    if (verbose) {
      amrex::Print() << "... Computing diffusion terms at t^(n+1," << sub_iteration+1 << ")" << std::endl;
    }
    FillPatch(*this, Sborder, nGrowTr, time + dt, State_Type, 0, NUM_STATE);
    Real flux_factor_new = sub_iteration==sub_ncycle-1 ? 0.5 : 0;
    getMOLSrcTerm(Sborder,*new_sources[diff_src],time,dt,flux_factor_new);
  }

  // Build other (neither spray nor diffusion) sources at t_new
  for (int n = 0; n < src_list.size(); ++n)
  {
    if (src_list[n] != diff_src
#ifdef AMREX_PARTICLES
      && src_list[n] != spray_src
#endif
      )
    {
      construct_new_source(src_list[n], time + dt, dt, amr_iteration, amr_ncycle, sub_iteration, sub_ncycle);
    }
  }

#ifdef AMREX_PARTICLES
  if (do_spray_particles)
  {
    // Advance the particle velocities by dt/2 to the new time.
    if (particle_verbose)
      amrex::Print() << "moveKick ... updating velocity only\n";

    if (! fill_Sborder) {
      Print() << "Sborder already defined, has nGrow = " << Sborder.nGrow()
              << ", needs "<< nGrow_Sborder << std::endl;
      Abort("Sborder should already be defined, why are we defining here again?");
      Sborder.define(grids, dmap, NUM_STATE, nGrow_Sborder);
    }

    if (!do_diffuse) { // Else, this was already done above.  No need to redo
      FillPatch(*this, Sborder, nGrow_Sborder, time + dt, State_Type, 0, NUM_STATE);
    }

    new_sources[spray_src]->setVal(0.);

    theSprayPC()->moveKick(Sborder, *new_sources[spray_src], level, dt, tmp_src_width);

    // Virtual particles will be recreated, so we need not kick them.
    // TODO: Is this true with SDC iterations??

    // Ghost particles need to be kicked except during the final iteration.
    if (amr_iteration != amr_ncycle)
      theGhostPC()->moveKick(Sborder, *new_sources[spray_src], level, dt, tmp_src_width);
  }
#endif

#ifdef REACTIONS
  // Update I_R and rebuild S_new accordingly
  if (do_react == 1)
  {
    react_state(time, dt);
  }
  else
  {
    construct_Snew(S_new, S_old, dt);
    get_new_data(Reactions_Type).setVal(0);
  }
#else
  construct_Snew(S_new, S_old, dt);
#endif


  computeTemp(S_new, ng_src);

  finalize_sdc_iteration(time, dt, amr_iteration, amr_ncycle, sub_iteration, sub_ncycle);

  return dt;
}

void
PeleC::construct_Snew(MultiFab& S_new, const MultiFab& S_old, Real dt)
{
  int ng = 0;

  MultiFab::Copy(S_new,S_old,0,0,NUM_STATE,ng);
  for (int n = 0; n < src_list.size(); ++n)
  {
    MultiFab::Saxpy(S_new,0.5*dt,*new_sources[src_list[n]],0,0,NUM_STATE,ng);
    MultiFab::Saxpy(S_new,0.5*dt,*old_sources[src_list[n]],0,0,NUM_STATE,ng);
  }
  if (do_hydro)
  {
    MultiFab::Saxpy(S_new,dt,hydro_source,0,0,NUM_STATE,ng);
  }
  
#ifdef REACTIONS
  if (do_react == 1)
  {
    MultiFab& I_R = get_new_data(Reactions_Type);
    MultiFab::Saxpy(S_new,dt,I_R,0,FirstSpec,NumSpec,0);
    MultiFab::Saxpy(S_new,dt,I_R,NumSpec,Eden,1,0);
  }
#endif
}

void
PeleC::initialize_sdc_iteration(Real time,
                                Real dt,
                                int  amr_iteration,
                                int  amr_ncycle,
                                int  sdc_iteration,
                                int  sdc_ncycle)
{
  BL_PROFILE("PeleC::initialize_sdc_iteration()");

  // Reset the change from density resets
  frac_change = 1;

  // Reset the grid loss tracking.
  if (track_grid_losses)
  {
    for (int i = 0; i < n_lost; i++) {
      material_lost_through_boundary_temp[i] = 0.0;
    }
  }
}

void
PeleC::finalize_sdc_iteration(Real time,
                              Real dt,
                              int  amr_iteration,
                              int  amr_ncycle,
                              int  sdc_iteration,
                              int  sdc_ncycle)
{
  BL_PROFILE("PeleC::finalize_sdc_iteration()");
}

void
PeleC::initialize_sdc_advance(Real time,
                              Real dt,
                              int  amr_iteration,
                              int  amr_ncycle)
{
  BL_PROFILE("PeleC::initialize_sdc_advance()");

  // Pass some information about the state of the simulation to a Fortran module.
  set_amr_info(level, amr_iteration, amr_ncycle, time, dt);

  for (int i = 0; i < num_state_type; ++i)
  {
    state[i].allocOldData();
    state[i].swapTimeLevels(dt);
  }

#ifdef REACTIONS
  if (do_react == 1)
  {
    // Initialize I_R with value from previous time step
    MultiFab::Copy(get_new_data(Reactions_Type),get_old_data(Reactions_Type),
                   0,0,get_new_data(Reactions_Type).nComp(),get_new_data(Reactions_Type).nGrow());
  }
#endif
}


void
PeleC::finalize_sdc_advance(Real time, Real dt, int amr_iteration, int amr_ncycle)
{
  BL_PROFILE("PeleC::finalize_sdc_advance()");

  // Add the material lost in this timestep to the cumulative losses.
  if (track_grid_losses)
  {
    ParallelDescriptor::ReduceRealSum(material_lost_through_boundary_temp, n_lost);

    for (int i = 0; i < n_lost; i++) {
      material_lost_through_boundary_cumulative[i] += material_lost_through_boundary_temp[i];
    }
  }

  Real cur_time = state[State_Type].curTime();
  set_special_tagging_flag(cur_time);
}
