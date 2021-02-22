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
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

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
  const amrex::MultiFab&
#ifdef PELEC_USE_EB
    state_old
#endif
  ,
  const amrex::MultiFab& /*state_new*/,
  amrex::MultiFab& ext_src,
  int ng)
{
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_lo =
    geom.ProbLoArray();

#ifdef PELEC_USE_EB
  auto const& fact =
    dynamic_cast<amrex::EBFArrayBoxFactory const&>(state_old.Factory());
  auto const& flags = fact.getMultiEBCellFlagFab();
#endif

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  for (amrex::MFIter mfi(ext_src, amrex::TilingIfNotGPU()); mfi.isValid();
       ++mfi) {
    const amrex::Box& bx = mfi.growntilebox(ng);

#ifdef PELEC_USE_EB
    const auto& flag_fab = flags[mfi];
    amrex::FabType typ = flag_fab.getType(bx);
    if (typ == amrex::FabType::covered) {
      continue;
    }
#endif

    // auto const& So = state_old.array(mfi);
    // auto const& Sn = state_new.array(mfi);
    auto const& Farr = ext_src.array(mfi);

    ProbParmDevice const* pp = PeleC::prob_parm_device.get();
    // Evaluate the external source
    amrex::ParallelFor(
      bx, NVAR, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
        const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
        const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
        const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];
        const amrex::Real rad = std::sqrt((x - 10) * (x - 10) + y * y + z * z);
        if (rad < 0.5) {
          const amrex::Real rho = 0.0013;
          const amrex::Real uin = 1000.0;
          const amrex::Real area = dx[1] * dx[2];
          const amrex::Real vol = dx[0] * dx[1] * dx[2];
          const amrex::Real mdot = rho * uin * area;
          Farr(i, j, k, URHO) = mdot / vol;
          Farr(i, j, k, UMX) = mdot * uin / vol;
          Farr(i, j, k, UEDEN) = mdot * uin * uin / vol;
        }
      });
  }
}
