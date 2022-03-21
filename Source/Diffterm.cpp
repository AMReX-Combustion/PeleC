#include "Diffterm.H"

// This file contains the driver for generating the diffusion fluxes, which are
// then used to generate the diffusion flux divergence.
// pc_compute_diffusion_flux utilizes functions from the pc_diffterm_3D
// header: pc_move_transcoefs_to_ec -> Moves Cell Centered Transport
// Coefficients to Edge Centers pc_compute_tangential_vel_derivs -> Computes
// the Tangential Velocity Derivatives pc_diffusion_flux -> Computes the
// diffusion flux per direction with the coefficients and velocity derivatives.

void
pc_compute_diffusion_flux(
  const amrex::Box& box,
  const amrex::Array4<const amrex::Real>& q,
  const amrex::Array4<const amrex::Real>& coef,
  const amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM> flx,
  const amrex::GpuArray<const amrex::Array4<const amrex::Real>, AMREX_SPACEDIM>
    area,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> del,
  const int do_harmonic,
  const amrex::FabType typ,
  const int Ncut,
  const EBBndryGeom* ebg,
  const amrex::Array4<amrex::EBCellFlag const>& flags)
{
  {
    // Compute Extensive diffusion fluxes for X, Y, Z
    BL_PROFILE("PeleC::diffusion_flux()");
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
      const amrex::Real delta = del[dir];
      amrex::Real d1 = 0.0;
      amrex::Real d2 = 0.0;
      amrex::Box ebox = amrex::surroundingNodes(box, dir);
      if (dir == 0) {
        AMREX_D_TERM(d2 = 1.;, d1 = del[1];, d2 = del[2];);
      } else if (dir == 1) {
        AMREX_D_TERM(d2 = 1.;, d1 = del[0];, d2 = del[2];);
      } else if (dir == 2) {
        d1 = del[0];
        d2 = del[1];
      }

      amrex::FArrayBox tander_ec(ebox, GradUtils::nCompTan);
      amrex::Elixir tander_eli = tander_ec.elixir();
      auto const& tander = tander_ec.array();
      amrex::ParallelFor(
        ebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          pc_compute_tangential_vel_derivs(i, j, k, q, dir, d1, d2, tander);
        });

      // Reset tangential derivatives to avoid using covered (invalid) data
      if (typ == amrex::FabType::singlevalued) {
        if (Ncut > 0) {
          BL_PROFILE("PeleC::pc_compute_tangential_vel_derivs_eb()");
          pc_compute_tangential_vel_derivs_eb(
            ebox, dir, d1, d2, ebg, Ncut, q, flags, tander);
        }
      } else if (typ == amrex::FabType::multivalued) {
        amrex::Abort(
          "multi-valued eb tangential derivatives to be implemented");
      }

      amrex::ParallelFor(
        ebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          amrex::Real c[dComp_lambda + 1];
          for (int n = 0; n < dComp_lambda + 1; n++) {
            pc_move_transcoefs_to_ec(i, j, k, n, coef, c, dir, do_harmonic);
          }
          pc_diffusion_flux(
            i, j, k, q, c, tander, area[dir], flx[dir], delta, dir);
        });
    }
  }
}
