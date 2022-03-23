#ifdef _OPENMP
#include <omp.h>
#endif

#include <masa.h>

#include "PeleC.H"
#include "IndexDefines.H"

using std::string;
using namespace MASA;

void
PeleC::construct_old_mms_source(amrex::Real time)
{
  const auto& S_old = get_old_data(State_Type);

  int ng = 0; // None filled

  old_sources[mms_src]->setVal(0.0);

  fill_mms_source(time, S_old, *old_sources[mms_src], ng);

  old_sources[mms_src]->FillBoundary(geom.periodicity());
}

void
PeleC::construct_new_mms_source(amrex::Real time)
{
  const auto& S_old = get_old_data(State_Type);

  int ng = 0;

  new_sources[mms_src]->setVal(0.0);

  fill_mms_source(time, S_old, *new_sources[mms_src], ng);
}

void
PeleC::fill_mms_source(
  amrex::Real time, const amrex::MultiFab& S, amrex::MultiFab& mms_src, int ng)
{
  BL_PROFILE("PeleC::fill_mms_source()");

  if (do_mms == 0) {
    return;
  }

  if (mms_src_evaluated) {
    if (verbose && amrex::ParallelDescriptor::IOProcessor()) {
      amrex::Print() << "... Reusing MMS source at time " << time << std::endl;
    }
    amrex::MultiFab::Copy(mms_src, mms_source, 0, 0, NVAR, ng);
    return;
  }

  if (verbose && amrex::ParallelDescriptor::IOProcessor()) {
    amrex::Print() << "... Computing MMS source at time " << time << std::endl;
  }

#ifdef PELEC_USE_MASA

  // Store the source for later reuse
  mms_source.define(grids, dmap, NVAR, ng);
  mms_source.setVal(0.0);

  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_lo =
    geom.ProbLoArray();

  auto const& fact =
    dynamic_cast<amrex::EBFArrayBoxFactory const&>(S.Factory());
  auto const& flags = fact.getMultiEBCellFlagFab();

  // FIXME: Reuse fillpatched data used for adv and diff...
  amrex::FillPatchIterator fpi(
    *this, mms_source, ng, time, State_Type, 0, NVAR);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  {
    for (amrex::MFIter mfi(mms_source, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
      const amrex::Box& bx = mfi.growntilebox(ng);

      const auto& flag_fab = flags[mfi];
      amrex::FabType typ = flag_fab.getType(bx);
      if (typ == amrex::FabType::covered) {
        continue;
      }

      auto const& s = S.array(mfi);
      auto const& src = mms_source.array(mfi);

      // Evaluate the source using MASA
      amrex::ParallelFor(
        bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
          const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
          const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

          src(i, j, k, URHO) = masa_eval_3d_source_rho(x, y, z);
          src(i, j, k, UMX) = masa_eval_3d_source_rho_u(x, y, z);
          src(i, j, k, UMY) = masa_eval_3d_source_rho_v(x, y, z);
          src(i, j, k, UMZ) = masa_eval_3d_source_rho_w(x, y, z);
          src(i, j, k, UEDEN) = masa_eval_3d_source_rho_e(x, y, z);
          for (int n = 0; n < NUM_SPECIES; n++) {
            src(i, j, k, UFS + n) = src(i, j, k, URHO) * s(i, j, k, UFS + n) /
                                    masa_eval_3d_exact_rho(x, y, z);
          }
        });
    }
  }

  amrex::MultiFab::Copy(mms_src, mms_source, 0, 0, NVAR, ng);
  mms_src_evaluated = true;

#else
  amrex::Error("MASA is not turned on. Turn on with PELEC_USE_MASA=TRUE.");
#endif
}
