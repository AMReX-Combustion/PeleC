#include "PeleC.H"
#include "Forcing.H"

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

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  for (amrex::MFIter mfi(forcing_src, amrex::TilingIfNotGPU()); mfi.isValid();
       ++mfi) {
    const amrex::Box& bx = mfi.growntilebox(ng);

    const auto& flag_fab = flags[mfi];
    amrex::FabType typ = flag_fab.getType(bx);
    if (typ == amrex::FabType::covered) {
      continue;
    }

    auto const& sarr = state_new.array(mfi);
    auto const& src = forcing_src.array(mfi);

    amrex::Real u0 = forcing_u0;
    amrex::Real v0 = forcing_v0;
    amrex::Real w0 = forcing_w0;
    amrex::Real force = forcing_force;

    // Evaluate the linear forcing term
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      src(i, j, k, UMX) =
        force * sarr(i, j, k, URHO) * (sarr(i, j, k, UMX) - u0);
      src(i, j, k, UMY) =
        force * sarr(i, j, k, URHO) * (sarr(i, j, k, UMY) - v0);
      src(i, j, k, UMZ) =
        force * sarr(i, j, k, URHO) * (sarr(i, j, k, UMZ) - w0);
    });
  }
}
