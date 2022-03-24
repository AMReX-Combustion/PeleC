#include "PeleC.H"
#include "IndexDefines.H"

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

  // auto const& Sos = state_old.const_arrays();
  // auto const& Sns = state_new.const_arrays();
  auto const& Farrs = ext_src.arrays();
  auto const& flagarrs = flags.const_arrays();
  const amrex::IntVect ngs(ng);
  amrex::ParallelFor(
    ext_src, ngs, NVAR,
    [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) noexcept {
      if (!flagarrs[nbx](i, j, k).isCovered()) {
        Farrs[nbx](i, j, k, n) = 0.0;
      }
    });
}
