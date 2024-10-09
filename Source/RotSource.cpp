#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

#include "PeleC.H"
#include "IndexDefines.H"

void
PeleC::construct_old_rot_source(amrex::Real /*time*/, amrex::Real /*dt*/)
{
  const amrex::MultiFab& S_old = get_old_data(State_Type);

  int ng = 0;

  old_sources[rot_src]->setVal(0.0);

  if (do_rf) {
    fill_rot_source(S_old, S_old, *old_sources[rot_src], ng);
    old_sources[rot_src]->FillBoundary(geom.periodicity());
  }
}

void
PeleC::construct_new_rot_source(amrex::Real /*time*/, amrex::Real /*dt*/)
{
  const amrex::MultiFab& S_old = get_old_data(State_Type);
  const amrex::MultiFab& S_new = get_new_data(State_Type);

  int ng = 0;

  new_sources[rot_src]->setVal(0.0);

  if (do_rf) {
    fill_rot_source(S_old, S_new, *new_sources[rot_src], ng);
  }
}

#if AMREX_SPACEDIM > 2
void
PeleC::fill_rot_source(
  const amrex::MultiFab& state_old,
  const amrex::MultiFab& state_new,
  amrex::MultiFab& rot_src,
  int ng)
{
  auto const& fact =
    dynamic_cast<amrex::EBFArrayBoxFactory const&>(state_old.Factory());
  auto const& flags = fact.getMultiEBCellFlagFab();

  const auto geomdata = geom.data();
  const bool rf_on = do_rf;
  amrex::Real omega = rf_omega;
  int axis = rf_axis;

  // need a check to prevent mrf in 2 and 1d
  // read a vector instead for mrf_axis_x and y
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> axis_loc = {
    rf_axis_x, rf_axis_y, rf_axis_z};

  auto const& sarrs = state_new.const_arrays();
  auto const& srcs = rot_src.arrays();
  auto const& flagarrs = flags.const_arrays();
  const amrex::IntVect ngs(ng);
  amrex::ParallelFor(
    rot_src, ngs, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
      if (!flagarrs[nbx](i, j, k).isCovered()) {
        const auto& sarr = sarrs[nbx];
        const auto& src = srcs[nbx];

        amrex::RealVect r(AMREX_D_DECL(0.0, 0.0, 0.0));
        amrex::RealVect w(AMREX_D_DECL(0.0, 0.0, 0.0));
        amrex::RealVect v(AMREX_D_DECL(0.0, 0.0, 0.0));

        amrex::IntVect iv(AMREX_D_DECL(i, j, k));
        r = get_rotaxis_vec(iv, axis_loc, geomdata);
        w[axis] = omega;

        if (rf_on) {
          v[0] = sarr(i, j, k, UMX) / sarr(i, j, k, URHO);
          v[1] = sarr(i, j, k, UMY) / sarr(i, j, k, URHO);
          v[2] = sarr(i, j, k, UMZ) / sarr(i, j, k, URHO);

          // non-inertial frame case
          amrex::RealVect w_cross_v = w.crossProduct(v);
          amrex::RealVect w_cross_r = w.crossProduct(r);
          amrex::RealVect w_cross_w_cross_r = w.crossProduct(w_cross_r);

          src(i, j, k, UMX) =
            -sarr(i, j, k, URHO) * (2.0 * w_cross_v[0] + w_cross_w_cross_r[0]);
          src(i, j, k, UMY) =
            -sarr(i, j, k, URHO) * (2.0 * w_cross_v[1] + w_cross_w_cross_r[1]);

          src(i, j, k, UMZ) =
            -sarr(i, j, k, URHO) * (2.0 * w_cross_v[2] + w_cross_w_cross_r[2]);
        }
      }
    });

  amrex::Gpu::synchronize();
}
#else
void
PeleC::fill_rot_source(
  const amrex::MultiFab& /*state_old*/,
  const amrex::MultiFab& /*state_new*/,
  amrex::MultiFab& /*rot_src*/,
  int /*ng*/)
{
}
#endif
