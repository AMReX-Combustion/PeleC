#include <AMReX_DistributionMapping.H>

#include "PeleC.H"
#include "React.H"

#ifdef PELEC_USE_EXPLICIT_REACT
void
PeleC::react_state_explicit(
  amrex::Real time, amrex::Real dt, bool react_init, amrex::MultiFab* A_aux)
{
  /*
    Update I_R, and recompute S_new
   */
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

  amrex::MultiFab Atmp, *Ap;

  if (A_aux == nullptr || react_init) {
    Atmp.define(grids, dmap, NVAR, ng, amrex::MFInfo(), Factory());
    Atmp.setVal(0);
    Ap = &Atmp;
  }

  if (!react_init) {
    // Build non-reacting source term, and an S_new that does not include
    // reactions

    if (A_aux == nullptr) {
      for (int n = 0; n < src_list.size(); ++n) {
        amrex::MultiFab::Saxpy(
          Atmp, 0.5, *new_sources[src_list[n]], 0, 0, NVAR, ng);
        amrex::MultiFab::Saxpy(
          Atmp, 0.5, *old_sources[src_list[n]], 0, 0, NVAR, ng);
      }
      if (do_hydro && !do_mol) {
        amrex::MultiFab::Add(Atmp, hydro_source, 0, 0, NVAR, ng);
      }
    } else {
      Ap = A_aux;
    }

    amrex::MultiFab& S_old = get_old_data(State_Type);
    amrex::MultiFab::Copy(S_new, S_old, 0, 0, NVAR, ng);
    amrex::MultiFab::Saxpy(S_new, dt, *Ap, 0, 0, NVAR, ng);
  }

  amrex::MultiFab& reactions = get_new_data(Reactions_Type);
  reactions.setVal(0.0);
  prefetchToDevice(reactions);
  if (use_reactions_work_estimate) {
    amrex::Abort("Need to implement redistribution of chemistry work");
  }

#ifdef AMREX_USE_EB
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

      auto const& uold =
        react_init ? S_new.array(mfi) : get_old_data(State_Type).array(mfi);
      auto const& unew = S_new.array(mfi);
      auto const& a = Ap->array(mfi);
      amrex::FArrayBox w(bx, 1);
      amrex::Elixir w_eli = w.elixir();
      auto const& w_arr = w.array();
      auto const& I_R = reactions.array(mfi);
      const int do_update =
        react_init ? 0 : 1; // TODO: Update here? Or just get reaction source?

#ifdef AMREX_USE_EB
      const auto& flag_fab = flags[mfi];
      amrex::FabType typ = flag_fab.getType(bx);
      if (typ == amrex::FabType::covered) {
        continue;
      } else if (
        typ == amrex::FabType::singlevalued || typ == amrex::FabType::regular)
#endif
      {
        if (chem_integrator == 1) {
          amrex::Abort("Implicit Chemistry is not implemented yet on GPU,  "
                       "only explicit (use pelec.chem_integrator=2).");
          /*                pc_react_state(ARLIM_3D(bx.loVect()),
             ARLIM_3D(bx.hiVect()), uold.dataPtr(),  ARLIM_3D(uold.loVect()),
             ARLIM_3D(uold.hiVect()), unew.dataPtr(),  ARLIM_3D(unew.loVect()),
             ARLIM_3D(unew.hiVect()), a.dataPtr(),     ARLIM_3D(a.loVect()),
             ARLIM_3D(a.hiVect()), m.dataPtr(),     ARLIM_3D(m.loVect()),
             ARLIM_3D(m.hiVect()), w.dataPtr(),     ARLIM_3D(w.loVect()),
             ARLIM_3D(w.hiVect()), I_R.dataPtr(),   ARLIM_3D(I_R.loVect()),
             ARLIM_3D(I_R.hiVect()), time, dt, do_update); */
        } else {
          const int nsubsteps_min = adaptrk_nsubsteps_min;
          const int nsubsteps_max = adaptrk_nsubsteps_max;
          const int nsubsteps_guess = adaptrk_nsubsteps_guess;
          const amrex::Real errtol = adaptrk_errtol;

          amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              pc_expl_reactions(
                i, j, k, uold, unew, a, w_arr, I_R, dt, nsubsteps_min, nsubsteps_max,
                nsubsteps_guess, errtol, do_update);
            });
        }

        if (do_react_load_balance || do_mol_load_balance) {
          get_new_data(Work_Estimate_Type)[mfi].plus<amrex::RunOn::Device>(w);
        }
      }
    }
  }

  if (ng > 0)
    S_new.FillBoundary(geom.periodicity());

  if (verbose > 1) {

    const int IOProc = amrex::ParallelDescriptor::IOProcessorNumber();
    amrex::Real run_time = amrex::ParallelDescriptor::second() - strt_time;

#ifdef AMREX_LAZY
    Lazy::QueueReduction([=]() mutable {
#endif
      amrex::ParallelDescriptor::ReduceRealMax(run_time, IOProc);

      if (amrex::ParallelDescriptor::IOProcessor())
        amrex::Print() << "PeleC::react_state() time = " << run_time << "\n";
#ifdef AMREX_LAZY
    });
#endif
  }
}
#endif
