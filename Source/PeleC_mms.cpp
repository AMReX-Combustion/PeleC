
#include <PeleC.H>
#include <PeleC_F.H>

using std::string;
using namespace amrex;

#ifdef _OPENMP
#include <omp.h>
#endif

#include <masa.h>
using namespace MASA;

void
PeleC::construct_old_mms_source (Real time)
{
  auto& S_old = get_old_data(State_Type);
  
  int ng = 0; // None filled

  old_sources[mms_src]->setVal(0.0);

  fill_mms_source(time, S_old, *old_sources[mms_src], ng);

  old_sources[mms_src]->FillBoundary(geom.periodicity());

}

void
PeleC::construct_new_mms_source(Real time)
{
    auto& S_old = get_old_data(State_Type);

    int ng = 0;

    new_sources[mms_src]->setVal(0.0);

    fill_mms_source(time, S_old, *new_sources[mms_src], ng);

}

// **********************************************************************************************
/**
 * Calculate the MMS term by calling the MASA library
 * Across all conserved state components, compute the MMS "source term"
 **/
void
PeleC::fill_mms_source (Real time, const MultiFab& S, MultiFab& mms_src, int ng)
{
    BL_PROFILE("PeleC::fill_mms_source()");

    if (do_mms == 0){
      return;
    }

    else if (mms_src_evaluated){
      if (verbose && ParallelDescriptor::IOProcessor()) {
        std::cout << "... Reusing MMS source at time " << time << std::endl;
      }
      MultiFab::Copy(mms_src, mms_source, 0, 0, NUM_STATE, ng);
      return;
    }

    else{
      if (verbose && ParallelDescriptor::IOProcessor()) {
        std::cout << "... Computing MMS source at time " << time << std::endl;
      }

#ifdef USE_MASA

      // Store the source for later reuse
      mms_source.define(grids,dmap, NUM_STATE, ng);
      mms_source.setVal(0.0);

      const Real* dx = geom.CellSize();
      const Real* prob_lo = geom.ProbLo();

#ifdef PELE_USE_EB
      auto const& fact = dynamic_cast<EBFArrayBoxFactory const&>(S.Factory());
      auto const& flags = fact.getMultiEBCellFlagFab();
#endif

      // FIXME: Reuse fillpatched data used for adv and diff...
      FillPatchIterator fpi(*this, mms_source, ng, time, State_Type, 0, NUM_STATE);

#ifdef _OPENMP
#pragma omp parallel
#endif
      {
        for (MFIter mfi(mms_source, MFItInfo().EnableTiling(hydro_tile_size).SetDynamic(true)); mfi.isValid(); ++mfi)
        {
          const Box& bx = mfi.growntilebox(ng);

#ifdef PELE_USE_EB
          const auto& flag_fab = flags[mfi];
          FabType typ = flag_fab.getType(bx);
          if (typ == FabType::covered) {
            continue;
          }
#endif

          const auto& Sfab = S[mfi];

          pc_mms_src(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
                     BL_TO_FORTRAN_3D(Sfab),
                     BL_TO_FORTRAN_3D(mms_source[mfi]),
                     ZFILL(prob_lo),ZFILL(dx),&time);
        }
      }

      MultiFab::Copy(mms_src, mms_source, 0, 0, NUM_STATE, ng);
      mms_src_evaluated = true;

#else
      amrex::Error("MASA is not turned on. Turn on with USE_MASA=TRUE.");
#endif
    }
}
