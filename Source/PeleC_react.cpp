
#include "PeleC.H"
#include "PeleC_F.H"

#include "AMReX_DistributionMapping.H"

#ifdef PELE_USE_EB
#include <AMReX_MultiCutFab.H>
#endif

using std::string;
using namespace amrex;

void
PeleC::react_state(Real time, Real dt, bool react_init, MultiFab* A_aux)
{
  /*
    Update I_R, and recompute S_new
   */
    BL_PROFILE("PeleC::react_state()");

    const Real strt_time = ParallelDescriptor::second();

    BL_ASSERT(do_react == 1);

    if (verbose && ParallelDescriptor::IOProcessor()) {
        if (react_init) {
            std::cout << "... Initializing reactions, using interval dt = " << dt << std::endl;
        }
        else {
            std::cout << "... Computing reactions for dt = " << dt << std::endl;
        }
    }

    MultiFab& S_new = get_new_data(State_Type);

    // Build the burning mask, in case the state has ghost zones.

    const int ng = S_new.nGrow();
    auto interior_mask = build_interior_boundary_mask(ng);

    // Create a MultiFab with all of the non-reacting source terms.

    MultiFab Atmp, *Ap;

    if (A_aux == nullptr || react_init)
    {
      Atmp.define(grids, dmap, NUM_STATE, ng, MFInfo(), Factory());
      Atmp.setVal(0);
      Ap = &Atmp;
    }

    if (!react_init)
    {
      // Build non-reacting source term, and an S_new that does not include reactions

      if (A_aux == nullptr)
      {
        for (int n = 0; n < src_list.size(); ++n)
        {
          MultiFab::Saxpy(Atmp,0.5,*new_sources[src_list[n]],0,0,NUM_STATE,ng);
          MultiFab::Saxpy(Atmp,0.5,*old_sources[src_list[n]],0,0,NUM_STATE,ng);
        }
        if (do_hydro && !do_mol_AD)
        {
          MultiFab::Add(Atmp,hydro_source,0,0,NUM_STATE,ng);
        }
      }
      else
      {
        Ap = A_aux;
      }

      MultiFab& S_old = get_old_data(State_Type);
      MultiFab::Copy(S_new,S_old,0,0,NUM_STATE,ng);
      MultiFab::Saxpy(S_new,dt,*Ap,0,0,NUM_STATE,ng);
    }

    MultiFab& reactions = get_new_data(Reactions_Type);
    reactions.setVal(0.0);

    if (use_reactions_work_estimate) {
	amrex::Abort("Need to implement redistribution of chemistry work");
    }

#ifdef _OPENMP
#pragma omp parallel
#endif
    {

        FArrayBox w;
        for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
        {

            const Box& bx = mfi.growntilebox(ng);

            const FArrayBox& uold = react_init ? S_new[mfi] : get_old_data(State_Type)[mfi];
            FArrayBox& unew       = S_new[mfi];
            FArrayBox& a          = (*Ap)[mfi];
            const IArrayBox& m    = (*interior_mask)[mfi];
            w.resize(bx,1);
            FArrayBox& I_R        = reactions[mfi];
            int do_update         = react_init ? 0 : 1;  // TODO: Update here? Or just get reaction source?

#ifdef PELE_USE_EB
            const EBFArrayBox& ufab = static_cast<const EBFArrayBox&>(unew);
            const auto& flag_fab = ufab.getEBCellFlagFab();
            FabType typ = flag_fab.getType(bx);
            if (typ == FabType::singlevalued || typ == FabType::regular) {
#else
            {
#endif

            if(chem_integrator==1)
            {
                pc_react_state(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
                        uold.dataPtr(),  ARLIM_3D(uold.loVect()),  ARLIM_3D(uold.hiVect()),
                        unew.dataPtr(),  ARLIM_3D(unew.loVect()),  ARLIM_3D(unew.hiVect()),
                        a.dataPtr(),     ARLIM_3D(a.loVect()),     ARLIM_3D(a.hiVect()),
                        m.dataPtr(),     ARLIM_3D(m.loVect()),     ARLIM_3D(m.hiVect()),
                        w.dataPtr(),     ARLIM_3D(w.loVect()),     ARLIM_3D(w.hiVect()),
                        I_R.dataPtr(),   ARLIM_3D(I_R.loVect()),   ARLIM_3D(I_R.hiVect()),
#ifdef PELE_USE_EB
                        BL_TO_FORTRAN_ANYD(flag_fab),
#endif
                        time, dt, do_update);
            }
            else
            {
                //right now hard-coding number of substeps when doing explicit 
                //reaction update. This will be changed later to an adaptive RK
                //approach and will no longer be used.
                int nsubsteps_min=20;
                int nsubsteps_max=200;
                int nsubsteps_guess=100;

                pc_react_state_expl(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
                        uold.dataPtr(),  ARLIM_3D(uold.loVect()),  ARLIM_3D(uold.hiVect()),
                        unew.dataPtr(),  ARLIM_3D(unew.loVect()),  ARLIM_3D(unew.hiVect()),
                        a.dataPtr(),     ARLIM_3D(a.loVect()),     ARLIM_3D(a.hiVect()),
                        m.dataPtr(),     ARLIM_3D(m.loVect()),     ARLIM_3D(m.hiVect()),
                        w.dataPtr(),     ARLIM_3D(w.loVect()),     ARLIM_3D(w.hiVect()),
                        I_R.dataPtr(),   ARLIM_3D(I_R.loVect()),   ARLIM_3D(I_R.hiVect()),
                        time, dt, do_update,adaptrk_nsubsteps_min,adaptrk_nsubsteps_max,adaptrk_nsubsteps_guess,adaptrk_errtol);
            }


            if (do_react_load_balance || do_mol_load_balance)
            {
                get_new_data(Work_Estimate_Type)[mfi].plus(w);
            }
            }
        }
    }
    if (ng > 0)
        S_new.FillBoundary(geom.periodicity());

    if (verbose > 1) {

        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
                ParallelDescriptor::ReduceRealMax(run_time, IOProc);

                if (ParallelDescriptor::IOProcessor())
                std::cout << "PeleC::react_state() time = " << run_time << "\n";
#ifdef BL_LAZY
                });
#endif

    }
}
