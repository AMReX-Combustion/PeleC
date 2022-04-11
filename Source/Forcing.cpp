#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

#include "PeleC.H"
#include "IndexDefines.H"

void
PeleC::construct_old_forcing_source(amrex::Real /*time*/, amrex::Real /*dt*/)
{
  const amrex::MultiFab& S_old = get_old_data(State_Type);

  int ng = 0;

  old_sources[forcing_src]->setVal(0.0);

  if (!add_forcing_src) {
    return;
  }

  fill_forcing_source(S_old, S_old, *old_sources[forcing_src], ng);

  old_sources[forcing_src]->FillBoundary(geom.periodicity());
}

void
PeleC::construct_new_forcing_source(amrex::Real /*time*/, amrex::Real /*dt*/)
{
  const amrex::MultiFab& S_old = get_old_data(State_Type);
  const amrex::MultiFab& S_new = get_new_data(State_Type);

  int ng = 0;

  new_sources[forcing_src]->setVal(0.0);

  if (!add_forcing_src) {
    return;
  }

  fill_forcing_source(S_old, S_new, *new_sources[forcing_src], ng);
}

void
PeleC::fill_forcing_source(
  const amrex::MultiFab& state_old
  /*unused*/,
  const amrex::MultiFab& state_new,
  amrex::MultiFab& forcing_src,
  int ng)
{
  auto const& fact =
    dynamic_cast<amrex::EBFArrayBoxFactory const&>(state_old.Factory());
  auto const& flags = fact.getMultiEBCellFlagFab();

  amrex::Real u0 = forcing_u0;
  amrex::Real v0 = forcing_v0;
  amrex::Real w0 = forcing_w0;
  amrex::Real force = forcing_force;

  auto const& sarrs = state_new.const_arrays();
  auto const& srcs = forcing_src.arrays();
  auto const& flagarrs = flags.const_arrays();
  const amrex::IntVect ngs(ng);
  amrex::ParallelFor(
    forcing_src, ngs,
    [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
      if (!flagarrs[nbx](i, j, k).isCovered()) {
        const auto& sarr = sarrs[nbx];
        const auto& src = srcs[nbx];
        src(i, j, k, UMX) =
          force * sarr(i, j, k, URHO) * (sarr(i, j, k, UMX) - u0);
        src(i, j, k, UMY) =
          force * sarr(i, j, k, URHO) * (sarr(i, j, k, UMY) - v0);
        src(i, j, k, UMZ) =
          force * sarr(i, j, k, URHO) * (sarr(i, j, k, UMZ) - w0);
      }
    });
  amrex::Gpu::synchronize();
}
