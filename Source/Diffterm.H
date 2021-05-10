#ifndef _DIFFTERM_H_
#define _DIFFTERM_H_

#include <AMReX_FArrayBox.H>
#ifdef PELEC_USE_EB
#include <AMReX_EBCellFlag.H>
#endif

#include "IndexDefines.H"
#include "EBStencilTypes.H"
#include "PelePhysics.H"
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

template <typename EOSType>
struct SpeciesEnergyFlux
{
  AMREX_GPU_HOST_DEVICE
  void operator()(
    const int i,
    const int j,
    const int k,
    const int im,
    const int jm,
    const int km,
    const amrex::Real dxinv,
    const amrex::Real coef[],
    const amrex::Array4<const amrex::Real>& q,
    const amrex::Array4<amrex::Real>& flx)
  {
    auto eos = pele::physics::PhysicsType::eos();

    // Get massfrac, molefrac, enthalpy
    amrex::Real mass1[NUM_SPECIES], mass2[NUM_SPECIES];
    amrex::Real mole1[NUM_SPECIES], mole2[NUM_SPECIES];
    amrex::Real hi1[NUM_SPECIES], hi2[NUM_SPECIES];
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
      mass1[ns] = q(i, j, k, ns + QFS);
    }
    eos.Y2X(mass1, mole1);
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
      mass2[ns] = q(im, jm, km, ns + QFS);
    }
    eos.Y2X(mass2, mole2);

    // Compute species and enthalpy fluxes for ideal EOS
    // Get species/enthalpy diffusion, compute correction vel
    amrex::Real T = q(i, j, k, QTEMP);
    eos.T2Hi(T, hi1);
    T = q(im, jm, km, QTEMP);
    eos.T2Hi(T, hi2);
    amrex::Real Vc = 0.0;
    const amrex::Real dpdx = dxinv * (q(i, j, k, QPRES) - q(im, jm, km, QPRES));
    const amrex::Real dlnp =
      dpdx / (0.5 * (q(i, j, k, QPRES) + q(im, jm, km, QPRES)));
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
  }
};

// FIXME remove once we compile in all the eos
#ifdef PELEC_USE_SRK
template <>
struct SpeciesEnergyFlux<pele::physics::eos::SRK>
{
  AMREX_GPU_HOST_DEVICE
  void operator()(
    const int i,
    const int j,
    const int k,
    const int im,
    const int jm,
    const int km,
    const amrex::Real dxinv,
    const amrex::Real coef[],
    const amrex::Array4<const amrex::Real>& q,
    const amrex::Array4<amrex::Real>& flx)
  {
    pele::physics::eos::SRK eos;

    // Get massfrac, molefrac, enthalpy
    amrex::Real mass1[NUM_SPECIES], mass2[NUM_SPECIES];
    amrex::Real mole1[NUM_SPECIES], mole2[NUM_SPECIES];
    amrex::Real hi1[NUM_SPECIES], hi2[NUM_SPECIES];
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
      mass1[ns] = q(i, j, k, ns + QFS);
    }
    eos.Y2X(mass1, mole1);
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
      mass2[ns] = q(im, jm, km, ns + QFS);
    }
    eos.Y2X(mass2, mole2);

    // Compute species and enthalpy fluxes accounting for nonideal EOS
    // Implementation note: nonideal EOS coeffs are evaluated at cell centers,
    // then inteprolated to faces. Could be slightly more efficient to
    // interpolate state to face, then do the EOS calls. Current implementation
    // matches old Fortran code.

    const amrex::Real rho1 = q(i, j, k, QRHO);
    const amrex::Real rho2 = q(im, jm, km, QRHO);
    const amrex::Real T1 = q(i, j, k, QTEMP);
    const amrex::Real T2 = q(im, jm, km, QTEMP);
    const amrex::Real dpdx = dxinv * (q(i, j, k, QPRES) - q(im, jm, km, QPRES));

    amrex::Real Vc = 0.0;
    amrex::Real diP1[NUM_SPECIES], dijY1[NUM_SPECIES][NUM_SPECIES];
    eos.RTY2transport(rho1, T1, mass1, diP1, dijY1);
    eos.RTY2Hi(rho1, T1, mass1, hi1);
    amrex::Real diP2[NUM_SPECIES], dijY2[NUM_SPECIES][NUM_SPECIES];
    eos.RTY2transport(rho2, T2, mass2, diP2, dijY2);
    eos.RTY2Hi(rho2, T2, mass2, hi2);
    amrex::Real dYdx[NUM_SPECIES], ddrive[NUM_SPECIES];
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
      dYdx[ns] = dxinv * (mass1[ns] - mass2[ns]);
      ddrive[ns] = 0.5 * (diP1[ns] + diP2[ns]) * dpdx;
    }
    for (int ns1 = 0; ns1 < NUM_SPECIES; ++ns1) {
      for (int ns2 = 0; ns2 < NUM_SPECIES; ++ns2) {
        ddrive[ns1] += 0.5 * (dijY1[ns1][ns2] + dijY2[ns1][ns2]) * dYdx[ns2];
      }
    }
    amrex::Real dsum = 0.0;
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
      dsum += ddrive[ns];
    }
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
      const amrex::Real Yface = 0.5 * (mass1[ns] + mass2[ns]);
      const amrex::Real hface = 0.5 * (hi1[ns] + hi2[ns]);
      ddrive[ns] -= Yface * dsum;
      const amrex::Real Vd = -coef[dComp_rhoD + ns] * ddrive[ns];
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
  }
};
#endif

struct FluxTypes
{
  using SpeciesEnergyFluxType = SpeciesEnergyFlux<pele::physics::EosType>;
};

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
  const amrex::Array4<const amrex::Real>& area,
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

  FluxTypes::SpeciesEnergyFluxType()(i, j, k, im, jm, km, dxinv, coef, q, flx);

  // Scale by area
  AMREX_D_TERM(flx(i, j, k, UMX) *= area(i, j, k);
               , flx(i, j, k, UMY) *= area(i, j, k);
               , flx(i, j, k, UMZ) *= area(i, j, k););
  flx(i, j, k, UEDEN) *= area(i, j, k);
  for (int ns = 0; ns < NUM_SPECIES; ++ns)
    flx(i, j, k, UFS + ns) *= area(i, j, k);
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
    area,
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
