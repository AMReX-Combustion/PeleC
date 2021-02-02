#ifndef _DIFFTERM_H_
#define _DIFFTERM_H_

#include <AMReX_FArrayBox.H>
#ifdef PELEC_USE_EB
#include <AMReX_EBCellFlag.H>
#endif

#include "IndexDefines.H"
#include "EBStencilTypes.H"
#include "EOS.H"
#include "Utilities.H"
#include "GradUtil.H"
#include "Diffusion.H"

// This header file contains functions and declarations for diffterm in 3D for
// PeleC GPU. As per the convention of AMReX, inlined device functions are
// defined here. Where as non-inline host functions are declared here.

// The original Diffterm function passes too many arguments, which contribute to
// kernel launch overhead. In order to combat this, we split some of the
// computation up into separate Device functions which will then be launched by
// a CPU function called pc_diffterm. pc_diffusion_flux calculates the diffusion
// flux per diction.

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_diffusion_flux(
  const int i,
  const int j,
  const int k,
  const amrex::Array4<const amrex::Real>& q,
  const amrex::Real coef[],
  const amrex::Array4<const amrex::Real>& td,
  const amrex::Array4<const amrex::Real>& a,
  const amrex::Array4<amrex::Real>& flx,
  const amrex::Real delta,
  const int dir)
{
  const int bdim[3] = {dir == 0, dir == 1, dir == 2};
  const int im = i - bdim[0];
  const int jm = j - bdim[1];
  const int km = k - bdim[2];
  const amrex::Real dxinv = 1.0 / delta;
  amrex::Real taux = 0.0, tauy = 0.0, tauz = 0.0;
  if (dir == 0) {
    AMREX_D_TERM(
      const amrex::Real dudx = dxinv * (q(i, j, k, QU) - q(im, jm, km, QU));
      , const amrex::Real dvdx = dxinv * (q(i, j, k, QV) - q(im, jm, km, QV));
      const amrex::Real dudy = td(i, j, k, 0);
      const amrex::Real dvdy = td(i, j, k, 1);
      , const amrex::Real dwdx = dxinv * (q(i, j, k, QW) - q(im, jm, km, QW));
      // const amrex::Real dwdy = td(i, j, k, 2);
      const amrex::Real dudz = td(i, j, k, 3);
      // const amrex::Real dvdz = td(i, j, k, 4);
      const amrex::Real dwdz = td(i, j, k, 5););
    const amrex::Real divu = AMREX_D_TERM(dudx, +dvdy, +dwdz);

    taux =
      coef[dComp_mu] * (2.0 * dudx - 2.0 / 3.0 * divu) + coef[dComp_xi] * divu;
    AMREX_D_TERM(, tauy = coef[dComp_mu] * (dudy + dvdx);
                 , tauz = coef[dComp_mu] * (dudz + dwdx););
  } else if (dir == 1) {
    const amrex::Real dudx = td(i, j, k, 0);
    const amrex::Real dvdx = td(i, j, k, 1);
    const amrex::Real dudy = dxinv * (q(i, j, k, QU) - q(im, jm, km, QU));
    const amrex::Real dvdy = dxinv * (q(i, j, k, QV) - q(im, jm, km, QV));
#if AMREX_SPACEDIM == 3
    const amrex::Real dwdy = dxinv * (q(i, j, k, QW) - q(im, jm, km, QW));
    // const amrex::Real dwdx = td(i, j, k, 2);
    // const amrex::Real dudz = td(i, j, k, 3);
    const amrex::Real dvdz = td(i, j, k, 4);
    const amrex::Real dwdz = td(i, j, k, 5);
    tauz = coef[dComp_mu] * (dwdy + dvdz);
#endif
    const amrex::Real divu = AMREX_D_TERM(dudx, +dvdy, +dwdz);

    taux = coef[dComp_mu] * (dudy + dvdx);
    tauy =
      coef[dComp_mu] * (2.0 * dvdy - 2.0 / 3.0 * divu) + coef[dComp_xi] * divu;
  } else if (dir == 2) {
    const amrex::Real dudx = td(i, j, k, 0);
    // const amrex::Real dvdx = td(i, j, k, 1);
    const amrex::Real dwdx = td(i, j, k, 2);
    // const amrex::Real dudy = td(i, j, k, 3);
    const amrex::Real dvdy = td(i, j, k, 4);
    const amrex::Real dwdy = td(i, j, k, 5);
    const amrex::Real dudz = dxinv * (q(i, j, k, QU) - q(im, jm, km, QU));
    const amrex::Real dvdz = dxinv * (q(i, j, k, QV) - q(im, jm, km, QV));
    const amrex::Real dwdz = dxinv * (q(i, j, k, QW) - q(im, jm, km, QW));
    const amrex::Real divu = dudx + dvdy + dwdz;

    taux = coef[dComp_mu] * (dudz + dwdx);
    tauy = coef[dComp_mu] * (dvdz + dwdy);
    tauz =
      coef[dComp_mu] * (2.0 * dwdz - 2.0 / 3.0 * divu) + coef[dComp_xi] * divu;
  }
  flx(i, j, k, UMX) = -taux;
  flx(i, j, k, UMY) = -tauy;
  flx(i, j, k, UMZ) = -tauz;
  flx(i, j, k, UEDEN) =
    0.5 * (AMREX_D_TERM(
            -taux * (q(i, j, k, QU) + q(im, jm, km, QU)),
            -tauy * (q(i, j, k, QV) + q(im, jm, km, QV)),
            -tauz * (q(i, j, k, QW) + q(im, jm, km, QW)))) -
    coef[dComp_lambda] * (dxinv * (q(i, j, k, QTEMP) - q(im, jm, km, QTEMP)));
  const amrex::Real dlnp = dxinv * (q(i, j, k, QPRES) - q(im, jm, km, QPRES)) /
                           (0.5 * (q(i, j, k, QPRES) + q(im, jm, km, QPRES)));

  // Get massfrac, molefrac, enthalpy
  amrex::Real mass1[NUM_SPECIES], mass2[NUM_SPECIES];
  amrex::Real mole1[NUM_SPECIES], mole2[NUM_SPECIES];
  amrex::Real hi1[NUM_SPECIES], hi2[NUM_SPECIES];
  for (int ns = 0; ns < NUM_SPECIES; ++ns) {
    mass1[ns] = q(i, j, k, ns + QFS);
  }
  EOS::Y2X(mass1, mole1);
  amrex::Real T = q(i, j, k, QTEMP);
  EOS::T2Hi(T, hi1);
  for (int ns = 0; ns < NUM_SPECIES; ++ns) {
    mass2[ns] = q(im, jm, km, ns + QFS);
  }
  T = q(im, jm, km, QTEMP);
  EOS::Y2X(mass2, mole2);
  EOS::T2Hi(T, hi2);

  // Get species/enthalpy diffusion, compute correction vel
  amrex::Real Vc = 0.0;
  for (int ns = 0; ns < NUM_SPECIES; ++ns) {
    const amrex::Real Xface = 0.5 * (mole1[ns] + mole2[ns]);
    const amrex::Real Yface = 0.5 * (mass1[ns] + mass2[ns]);
    const amrex::Real hface = 0.5 * (hi1[ns] + hi2[ns]);
    const amrex::Real dXdx = dxinv * (mole1[ns] - mole2[ns]);
    const amrex::Real Vd =
      -coef[dComp_rhoD + ns] * (dXdx + (Xface - Yface) * dlnp);
    flx(i, j, k, UFS + ns) = Vd;
    Vc += Vd;
    flx(i, j, k, UEDEN) += Vd * hface;
  }
  // Add correction velocity to fluxes
  for (int ns = 0; ns < NUM_SPECIES; ++ns) {
    const amrex::Real Yface = 0.5 * (mass1[ns] + mass2[ns]);
    const amrex::Real hface = 0.5 * (hi1[ns] + hi2[ns]);
    flx(i, j, k, UFS + ns) -= Yface * Vc;
    flx(i, j, k, UEDEN) -= Yface * hface * Vc;
  }
  // Scale by area
  AMREX_D_TERM(flx(i, j, k, UMX) *= a(i, j, k);
               , flx(i, j, k, UMY) *= a(i, j, k);
               , flx(i, j, k, UMZ) *= a(i, j, k););
  flx(i, j, k, UEDEN) *= a(i, j, k);
  for (int ns = 0; ns < NUM_SPECIES; ++ns)
    flx(i, j, k, UFS + ns) *= a(i, j, k);
}

// This function computes the flux divergence.
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_flux_div(
  const int i,
  const int j,
  const int k,
  const int n,
  AMREX_D_DECL(
    const amrex::Array4<const amrex::Real>& fx,
    const amrex::Array4<const amrex::Real>& fy,
    const amrex::Array4<const amrex::Real>& fz),
  const amrex::Array4<const amrex::Real>& V,
  const amrex::Array4<amrex::Real>& D)
{
  D(i, j, k, n) =
    -(AMREX_D_TERM(
      fx(i + 1, j, k, n) - fx(i, j, k, n), +fy(i, j + 1, k, n) - fy(i, j, k, n),
      +fz(i, j, k + 1, n) - fz(i, j, k, n))) /
    V(i, j, k);
}

void pc_compute_diffusion_flux(
  const amrex::Box& box,
  const amrex::Array4<const amrex::Real>& q,
  const amrex::Array4<const amrex::Real>& coef,
  const amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM> flx,
  const amrex::GpuArray<const amrex::Array4<const amrex::Real>, AMREX_SPACEDIM>
    a,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> del,
  const int do_harmonic
#ifdef PELEC_USE_EB
  ,
  const amrex::FabType typ,
  const int Ncut,
  const EBBndryGeom* ebg,
  const amrex::Array4<amrex::EBCellFlag const>& flags
#endif
);

#endif
