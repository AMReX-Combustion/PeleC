#include <memory> // for allocator

#include "AMReX_Array4.H"           // for Array4
#include "AMReX_EBCellFlag.H"       // for EBCellFlagFab
#include "AMReX_EBFabFactory.H"     // for EBFArrayBoxFactory
#include "AMReX_FabArray.H"         // for FabArray
#include "AMReX_FabFactory.H"       // for FabType, FabFactory, FabType::co...
#include "AMReX_Geometry.H"         // for Geometry
#include "AMReX_GpuLaunchFunctsC.H" // for ParallelFor
#include "AMReX_GpuQualifiers.H"    // for AMREX_GPU_DEVICE
#include "AMReX_MFIter.H"           // for MFIter, TilingIfNotGPU
#include "AMReX_MultiFab.H"         // for MultiFab
#include "AMReX_REAL.H"             // for Real
#include "AMReX_Vector.H"           // for Vector

#include "IndexDefines.H" // for NVAR
#include "PeleC.H"        // for PeleC, ext_src, State_Type

namespace amrex {
class Box;
} // namespace amrex

void
PeleC::construct_old_ext_source(amrex::Real time, amrex::Real dt)
{
  const amrex::MultiFab& S_old = get_old_data(State_Type);

  int ng = 0; // None filled

  old_sources[ext_src]->setVal(0.0);

  if (!add_ext_src) {
    return;
  }

  fill_ext_source(time, dt, S_old, S_old, *old_sources[ext_src], ng);

  old_sources[ext_src]->FillBoundary(geom.periodicity());
}

void
PeleC::construct_new_ext_source(amrex::Real time, amrex::Real dt)
{
  const amrex::MultiFab& S_old = get_old_data(State_Type);
  const amrex::MultiFab& S_new = get_new_data(State_Type);

  int ng = 0;

  new_sources[ext_src]->setVal(0.0);

  if (!add_ext_src) {
    return;
  }

  fill_ext_source(time, dt, S_old, S_new, *new_sources[ext_src], ng);
}

void
PeleC::fill_ext_source(
  amrex::Real /*time*/,
  amrex::Real /*dt*/,
  const amrex::MultiFab& state_old
  /*unused*/,
  const amrex::MultiFab& /*state_new*/,
  amrex::MultiFab& ext_src,
  int ng)
{
  auto const& fact =
    dynamic_cast<amrex::EBFArrayBoxFactory const&>(state_old.Factory());
  auto const& flags = fact.getMultiEBCellFlagFab();

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  for (amrex::MFIter mfi(ext_src, amrex::TilingIfNotGPU()); mfi.isValid();
       ++mfi) {
    const amrex::Box& bx = mfi.growntilebox(ng);

    const auto& flag_fab = flags[mfi];
    amrex::FabType typ = flag_fab.getType(bx);
    if (typ == amrex::FabType::covered) {
      continue;
    }

    // auto const& So = state_old.array(mfi);
    // auto const& Sn = state_new.array(mfi);
    auto const& Farr = ext_src.array(mfi);

    // Evaluate the external source
    amrex::ParallelFor(
      bx, NVAR, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
        Farr(i, j, k, n) = 0.0;
      });
  }
}
