#include "React.H"

void
PeleC::react_state(
  amrex::Real /*time*/,
  amrex::Real dt,
  bool react_init,
  amrex::MultiFab* aux_src)
{
  // Update I_R, and recompute S_new
  BL_PROFILE("PeleC::react_state()");

  const amrex::Real strt_time = amrex::ParallelDescriptor::second();

  AMREX_ASSERT(do_react == 1);

  if (verbose && amrex::ParallelDescriptor::IOProcessor()) {
    if (react_init) {
      amrex::Print() << "... Initializing reactions, using interval dt = " << dt
                     << std::endl;
    } else {
      amrex::Print() << "... Computing reactions for dt = " << dt << std::endl;
    }
  }

  amrex::MultiFab& S_new = get_new_data(State_Type);
  const int ng = S_new.nGrow();
  prefetchToDevice(S_new);

  // Create a MultiFab with all of the non-reacting source terms.
  amrex::MultiFab non_react_src_tmp;
  amrex::MultiFab* non_react_src = nullptr;

  if (react_init) {
    non_react_src_tmp.define(grids, dmap, NVAR, ng, amrex::MFInfo(), Factory());
    non_react_src_tmp.setVal(0);
    non_react_src = &non_react_src_tmp;
  } else {
    // Only do this if we are not at the first step
    // Build non-reacting source term, and an S_new that does not include
    // reactions
    if (aux_src == nullptr) {
      non_react_src_tmp.define(
        grids, dmap, NVAR, ng, amrex::MFInfo(), Factory());
      non_react_src_tmp.setVal(0);
      non_react_src = &non_react_src_tmp;

      for (int n = 0; n < src_list.size(); ++n) {
        amrex::MultiFab::Saxpy(
          non_react_src_tmp, 0.5, *new_sources[src_list[n]], 0, 0, NVAR, ng);
        amrex::MultiFab::Saxpy(
          non_react_src_tmp, 0.5, *old_sources[src_list[n]], 0, 0, NVAR, ng);
      }

      if (do_hydro && !do_mol) {
        amrex::MultiFab::Add(non_react_src_tmp, hydro_source, 0, 0, NVAR, ng);
      }
    } else {
      // in MOL update all non-reacting sources
      // are passed into auxillary sources
      non_react_src = aux_src;
    }

    // S_new = S_old + dt*(non reacting source terms)
    amrex::MultiFab& S_old = get_old_data(State_Type);
    amrex::MultiFab::Copy(S_new, S_old, 0, 0, NVAR, ng);
    amrex::MultiFab::Saxpy(S_new, dt, *non_react_src, 0, 0, NVAR, ng);
  }

  amrex::MultiFab& react_src = get_new_data(Reactions_Type);
  react_src.setVal(0.0);
  prefetchToDevice(react_src);

#ifdef USE_SUNDIALS_PP
  // for sundials box integration
  amrex::MultiFab STemp(grids, dmap, NUM_SPECIES + 2, 0);
  amrex::MultiFab extsrc_rY(grids, dmap, NUM_SPECIES, 0);
  amrex::MultiFab extsrc_rE(grids, dmap, 1, 0);
  amrex::iMultiFab dummyMask(grids, dmap, 1, 0);
  amrex::MultiFab fctCount(grids, dmap, 1, 0);
  dummyMask.setVal(1);

  if (!react_init) {
    amrex::MultiFab& S_old = get_old_data(State_Type);
    amrex::MultiFab::Copy(STemp, S_old, UFS, 0, NUM_SPECIES, STemp.nGrow());
    amrex::MultiFab::Copy(STemp, S_old, UTEMP, NUM_SPECIES, 1, STemp.nGrow());
    amrex::MultiFab::Copy(
      STemp, S_old, UEINT, NUM_SPECIES + 1, 1, STemp.nGrow());
  } else {
    amrex::MultiFab::Copy(STemp, S_new, UFS, 0, NUM_SPECIES, STemp.nGrow());
    amrex::MultiFab::Copy(STemp, S_new, UTEMP, NUM_SPECIES, 1, STemp.nGrow());
    amrex::MultiFab::Copy(
      STemp, S_new, UEINT, NUM_SPECIES + 1, 1, STemp.nGrow());
  }
  amrex::MultiFab::Copy(
    extsrc_rY, *non_react_src, UFS, 0, NUM_SPECIES, STemp.nGrow());
#endif

#ifdef PELEC_USE_EB
  auto const& fact =
    dynamic_cast<amrex::EBFArrayBoxFactory const&>(S_new.Factory());
  auto const& flags = fact.getMultiEBCellFlagFab();
#endif

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  {
    for (amrex::MFIter mfi(S_new, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {

      const amrex::Box& bx = mfi.growntilebox(ng);

      // old state or the state at t=0
      auto const& sold_arr =
        react_init ? S_new.array(mfi) : get_old_data(State_Type).array(mfi);

      // new state
      auto const& snew_arr = S_new.array(mfi);
      auto const& nonrs_arr = non_react_src->array(mfi);
      auto const& I_R = react_src.array(mfi);

      // only update beyond first step
      // TODO: Update here? Or just get reaction source?
      const int do_update = react_init ? 0 : 1;

      const int captured_clean_massfrac = clean_massfrac;

#ifdef PELEC_USE_EB
      const auto& flag_fab = flags[mfi];
      amrex::FabType typ = flag_fab.getType(bx);
      if (typ == amrex::FabType::covered) {
        if (do_react_load_balance) {
          const amrex::Box vbox = mfi.tilebox();
          get_new_data(Work_Estimate_Type)[mfi].plus<amrex::RunOn::Device>(
            0.0, vbox);
        }
        continue;
      }
      if (typ == amrex::FabType::singlevalued || typ == amrex::FabType::regular)
#endif
      {
        if (chem_integrator == 1) {
          // for rk64 we set minimum, maximum and guess
          // number of sub-iterations
          const int nsubsteps_min = adaptrk_nsubsteps_min;
          const int nsubsteps_max = adaptrk_nsubsteps_max;
          const int nsubsteps_guess = adaptrk_nsubsteps_guess;

          // for rk64 we set the error tolerance
          const amrex::Real errtol = adaptrk_errtol;

          amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              pc_expl_reactions(
                i, j, k, sold_arr, snew_arr, nonrs_arr, I_R, dt, nsubsteps_min,
                nsubsteps_max, nsubsteps_guess, errtol, do_update,
                captured_clean_massfrac);
            });
        } else if (chem_integrator == 2) {
#ifdef USE_SUNDIALS_PP
          amrex::Real wt =
            amrex::ParallelDescriptor::second(); // timing for each fab

          amrex::Real current_time = 0.0;

          auto const& rhoY = STemp.array(mfi);
          auto const& T = STemp.array(mfi, NUM_SPECIES);
          auto const& rhoE = STemp.array(mfi, NUM_SPECIES + 1);
          auto const& frcExt = extsrc_rY.array(mfi);
          auto const& frcEExt = extsrc_rE.array(mfi);
          auto const& mask = dummyMask.array(mfi);
          auto const& fc = fctCount.array(mfi);

          amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              // work on old state
              amrex::Real rhou = sold_arr(i, j, k, UMX);
              amrex::Real rhov = sold_arr(i, j, k, UMY);
              amrex::Real rhow = sold_arr(i, j, k, UMZ);
              amrex::Real rho_old = sold_arr(i, j, k, URHO);
              amrex::Real rhoInv = 1.0 / rho_old;

              amrex::Real e_old =
                (sold_arr(i, j, k, UEDEN) // total energy
                 -
                 0.5 * (rhou * rhou + rhov * rhov + rhow * rhow) * rhoInv) // KE
                * rhoInv;

              // work on new state
              rhou = snew_arr(i, j, k, UMX);
              rhov = snew_arr(i, j, k, UMY);
              rhow = snew_arr(i, j, k, UMZ);
              rhoInv = 1.0 / snew_arr(i, j, k, URHO);

              amrex::Real rhoedot_ext =
                (snew_arr(i, j, k, UEDEN) // new total energy
                 - 0.5 * (rhou * rhou + rhov * rhov + rhow * rhow) *
                     rhoInv // new KE
                 - rho_old * e_old) /
                dt;

              frcEExt(i, j, k) = rhoedot_ext;
              if (captured_clean_massfrac == 1) {
                clip_normalize_rYarr(i, j, k, sold_arr, rhoY);
              }
            });

          const int reactor_type = 1;
          react(
            bx, rhoY, frcExt, T, rhoE, frcEExt, fc, mask, dt, current_time,
            reactor_type
#ifdef AMREX_USE_GPU
            ,
            amrex::Gpu::gpuStream()
#endif
          );

          amrex::Gpu::Device::streamSynchronize();

          // unpack data
          amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              // work on old state
              amrex::Real rhou = sold_arr(i, j, k, UMX);
              amrex::Real rhov = sold_arr(i, j, k, UMY);
              amrex::Real rhow = sold_arr(i, j, k, UMZ);
              amrex::Real rho_old = sold_arr(i, j, k, URHO);
              amrex::Real rhoInv = 1.0 / rho_old;

              amrex::Real e_old =
                (sold_arr(i, j, k, UEDEN) // old total energy
                 -
                 0.5 * (rhou * rhou + rhov * rhov + rhow * rhow) * rhoInv) // KE
                * rhoInv;

              rhou = snew_arr(i, j, k, UMX);
              rhov = snew_arr(i, j, k, UMY);
              rhow = snew_arr(i, j, k, UMZ);
              rhoInv = 1.0 / snew_arr(i, j, k, URHO);

              amrex::Real rhoedot_ext =
                (snew_arr(i, j, k, UEDEN) // new total energy
                 -
                 0.5 * (rhou * rhou + rhov * rhov + rhow * rhow) * rhoInv // KE
                 - rho_old * e_old) // old internal energy
                / dt;

              amrex::Real umnew =
                sold_arr(i, j, k, UMX) + dt * nonrs_arr(i, j, k, UMX);
              amrex::Real vmnew =
                sold_arr(i, j, k, UMY) + dt * nonrs_arr(i, j, k, UMY);
              amrex::Real wmnew =
                sold_arr(i, j, k, UMZ) + dt * nonrs_arr(i, j, k, UMZ);

              // get new rho
              amrex::Real rhonew = 0.0;

              for (int nsp = 0; nsp < NUM_SPECIES; nsp++) {
                rhonew += rhoY(i, j, k, nsp);
              }

              if (do_update) {
                snew_arr(i, j, k, URHO) = rhonew;
                snew_arr(i, j, k, UMX) = umnew;
                snew_arr(i, j, k, UMY) = vmnew;
                snew_arr(i, j, k, UMZ) = wmnew;

                for (int nsp = 0; nsp < NUM_SPECIES; nsp++) {
                  snew_arr(i, j, k, UFS + nsp) = rhoY(i, j, k, nsp);
                }
                snew_arr(i, j, k, UTEMP) = T(i, j, k);

                snew_arr(i, j, k, UEINT) = rho_old * e_old + dt * rhoedot_ext;
                snew_arr(i, j, k, UEDEN) =
                  snew_arr(i, j, k, UEINT) +
                  0.5 * (umnew * umnew + vmnew * vmnew + wmnew * wmnew) /
                    rhonew;
              }

              for (int nsp = 0; nsp < NUM_SPECIES; nsp++) {
                I_R(i, j, k, nsp) = (rhoY(i, j, k, nsp)              // new rhoy
                                     - sold_arr(i, j, k, UFS + nsp)) // old rhoy
                                      / dt -
                                    nonrs_arr(i, j, k, UFS + nsp);
              }

              I_R(i, j, k, NUM_SPECIES) =
                (rho_old * e_old + dt * rhoedot_ext // new internal energy
                 + 0.5 * (umnew * umnew + vmnew * vmnew + wmnew * wmnew) /
                     rhonew                  // new KE
                 - sold_arr(i, j, k, UEDEN)) // old total energy
                  / dt -
                nonrs_arr(i, j, k, UEDEN);
            });

          wt = (amrex::ParallelDescriptor::second() - wt) / bx.d_numPts();

          if (do_react_load_balance) {
            const amrex::Box vbox = mfi.tilebox();
            get_new_data(Work_Estimate_Type)[mfi].plus<amrex::RunOn::Device>(
              wt, vbox);
          }
#else
          amrex::Abort(
            "chem_integrator=2 which requires Sundials to be enabled");
#endif
        } else {
          amrex::Abort("chem_integrator must be equal to 1 or 2");
        }

        // update heat release
        amrex::ParallelFor(
          bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            I_R(i, j, k, NUM_SPECIES + 1) = 0.0;
            auto eos = pele::physics::PhysicsType::eos();

            amrex::Real hi[NUM_SPECIES] = {0.0};

#ifndef PELEC_USE_SRK
            eos.T2Hi(snew_arr(i, j, k, UTEMP), hi);
#else
            amrex::Real Yspec[NUM_SPECIES] = {0.0};
            for(int nsp=0;nsp<NUM_SPECIES;nsp++)
            {
                Yspec[nsp]=snew_arr(i, j, k, UFS+nsp)/snew_arr(i, j, k, URHO);
            }
            eos.RTY2Hi(snew_arr(i,j,k,URHO), snew_arr(i,j,k,UTEMP), Yspec, hi);
#endif

            for (int nsp = 0; nsp < NUM_SPECIES; nsp++) {
              I_R(i, j, k, NUM_SPECIES + 1) -= hi[nsp] * I_R(i, j, k, nsp);
            }
          });
      }
    }
  }

  if (ng > 0) {
    S_new.FillBoundary(geom.periodicity());
  }

  if (verbose > 1) {
    const int IOProc = amrex::ParallelDescriptor::IOProcessorNumber();
    amrex::Real run_time = amrex::ParallelDescriptor::second() - strt_time;

#ifdef AMREX_LAZY
    Lazy::QueueReduction([=]() mutable {
#endif
      amrex::ParallelDescriptor::ReduceRealMax(run_time, IOProc);

      if (amrex::ParallelDescriptor::IOProcessor()) {
        amrex::Print() << "PeleC::react_state() time = " << run_time << "\n";
      }
#ifdef AMREX_LAZY
    });
#endif
  }
}
