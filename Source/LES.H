#ifndef LES_H
#define LES_H

#include "PeleC.H"
#include "Utilities.H"
#include "Diffusion.H"
#include "Diffterm.H"
#include "GradUtil.H"
#include "PelePhysics.H"
#include "Filter.H"

#ifdef AMREX_USE_OMP
#include <omp.h>
#endif

#if AMREX_SPACEDIM == 3
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
get_sfs_stresses_xdir(
  const int i,
  const int j,
  const int k,
  const amrex::Array4<const amrex::Real>& q,
  const amrex::Array4<const amrex::Real>& td,
  const amrex::Real dxinv,
  const amrex::Real deltabar,
  amrex::Real& alphaij_xx,
  amrex::Real& alphaij_xy,
  amrex::Real& alphaij_xz,
  amrex::Real& alpha,
  amrex::Real& flux_T)
{

  amrex::Real dUdx[3][3] = {{0.0}};
  dUdx[0][0] = dxinv * (q(i, j, k, QU) - q(i - 1, j, k, QU)); // dudx
  dUdx[0][1] = td(i, j, k, 0);                                // dudy
  dUdx[0][2] = td(i, j, k, 3);                                // dudz
  dUdx[1][0] = dxinv * (q(i, j, k, QV) - q(i - 1, j, k, QV)); // dvdx
  dUdx[1][1] = td(i, j, k, 1);                                // dvdy
  dUdx[1][2] = td(i, j, k, 4);                                // dvdz
  dUdx[2][0] = dxinv * (q(i, j, k, QW) - q(i - 1, j, k, QW)); // dwdx
  dUdx[2][1] = td(i, j, k, 2);                                // dwdy
  dUdx[2][2] = td(i, j, k, 5);                                // dwdz

  amrex::Real S[3][3] = {{0.0}};
  amrex::Real Sijmag = 0;
  for (int m = 0; m < 3; m++) {
    for (int n = 0; n < 3; n++) {
      S[m][n] = 0.5 * (dUdx[m][n] + dUdx[n][m]);
      Sijmag += S[m][n] * S[m][n];
    }
  }
  Sijmag = std::sqrt(2.0 * Sijmag);
  const amrex::Real Skk = S[0][0] + S[1][1] + S[2][2];

  const amrex::Real mut = 0.5 * (q(i, j, k, QRHO) + q(i - 1, j, k, QRHO)) *
                          deltabar * deltabar * Sijmag;

  alphaij_xx = 2.0 * mut * (S[0][0] - Skk / 3.0);
  alphaij_xy = 2.0 * mut * S[0][1];
  alphaij_xz = 2.0 * mut * S[0][2];
  alpha = 2.0 * mut * Sijmag;

  const amrex::Real dTdx = dxinv * (q(i, j, k, QTEMP) - q(i - 1, j, k, QTEMP));
  flux_T = mut * dTdx;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
get_sfs_stresses_ydir(
  const int i,
  const int j,
  const int k,
  const amrex::Array4<const amrex::Real>& q,
  const amrex::Array4<const amrex::Real>& td,
  const amrex::Real dxinv,
  const amrex::Real deltabar,
  amrex::Real& alphaij_yx,
  amrex::Real& alphaij_yy,
  amrex::Real& alphaij_yz,
  amrex::Real& alpha,
  amrex::Real& flux_T)
{

  amrex::Real dUdx[3][3] = {{0.0}};
  dUdx[0][0] = td(i, j, k, 0);                                // dudx
  dUdx[0][1] = dxinv * (q(i, j, k, QU) - q(i, j - 1, k, QU)); // dudy
  dUdx[0][2] = td(i, j, k, 3);                                // dudz
  dUdx[1][0] = td(i, j, k, 1);                                // dvdx
  dUdx[1][1] = dxinv * (q(i, j, k, QV) - q(i, j - 1, k, QV)); // dvdy
  dUdx[1][2] = td(i, j, k, 4);                                // dvdz
  dUdx[2][0] = td(i, j, k, 2);                                // dwdx
  dUdx[2][1] = dxinv * (q(i, j, k, QW) - q(i, j - 1, k, QW)); // dwdy
  dUdx[2][2] = td(i, j, k, 5);                                // dwdz

  amrex::Real S[3][3] = {{0.0}};
  amrex::Real Sijmag = 0;
  for (int m = 0; m < 3; m++) {
    for (int n = 0; n < 3; n++) {
      S[m][n] = 0.5 * (dUdx[m][n] + dUdx[n][m]);
      Sijmag += S[m][n] * S[m][n];
    }
  }
  Sijmag = std::sqrt(2.0 * Sijmag);
  const amrex::Real Skk = S[0][0] + S[1][1] + S[2][2];

  const amrex::Real mut = 0.5 * (q(i, j, k, QRHO) + q(i, j - 1, k, QRHO)) *
                          deltabar * deltabar * Sijmag;

  alphaij_yx = 2.0 * mut * S[1][0];
  alphaij_yy = 2.0 * mut * (S[1][1] - Skk / 3.0);
  alphaij_yz = 2.0 * mut * S[1][2];
  alpha = 2.0 * mut * Sijmag;

  const amrex::Real dTdy = dxinv * (q(i, j, k, QTEMP) - q(i, j - 1, k, QTEMP));
  flux_T = mut * dTdy;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
get_sfs_stresses_zdir(
  const int i,
  const int j,
  const int k,
  const amrex::Array4<const amrex::Real>& q,
  const amrex::Array4<const amrex::Real>& td,
  const amrex::Real dxinv,
  const amrex::Real deltabar,
  amrex::Real& alphaij_zx,
  amrex::Real& alphaij_zy,
  amrex::Real& alphaij_zz,
  amrex::Real& alpha,
  amrex::Real& flux_T)
{

  amrex::Real dUdx[3][3] = {{0.0}};
  dUdx[0][0] = td(i, j, k, 0);                                // dudx
  dUdx[0][1] = td(i, j, k, 3);                                // dudy
  dUdx[0][2] = dxinv * (q(i, j, k, QU) - q(i, j, k - 1, QU)); // dudz
  dUdx[1][0] = td(i, j, k, 1);                                // dvdx
  dUdx[1][1] = td(i, j, k, 4);                                // dvdy
  dUdx[1][2] = dxinv * (q(i, j, k, QV) - q(i, j, k - 1, QV)); // dvdz
  dUdx[2][0] = td(i, j, k, 2);                                // dwdx
  dUdx[2][1] = td(i, j, k, 5);                                // dwdy
  dUdx[2][2] = dxinv * (q(i, j, k, QW) - q(i, j, k - 1, QW)); // dwdz

  amrex::Real S[3][3] = {{0.0}};
  amrex::Real Sijmag = 0;
  for (int m = 0; m < 3; m++) {
    for (int n = 0; n < 3; n++) {
      S[m][n] = 0.5 * (dUdx[m][n] + dUdx[n][m]);
      Sijmag += S[m][n] * S[m][n];
    }
  }
  Sijmag = std::sqrt(2.0 * Sijmag);
  const amrex::Real Skk = S[0][0] + S[1][1] + S[2][2];

  const amrex::Real mut = 0.5 * (q(i, j, k, QRHO) + q(i, j, k - 1, QRHO)) *
                          deltabar * deltabar * Sijmag;

  alphaij_zx = 2.0 * mut * S[2][0];
  alphaij_zy = 2.0 * mut * S[2][1];
  alphaij_zz = 2.0 * mut * (S[2][2] - Skk / 3.0);
  alpha = 2.0 * mut * Sijmag;

  const amrex::Real dTdz = dxinv * (q(i, j, k, QTEMP) - q(i, j, k - 1, QTEMP));
  flux_T = mut * dTdz;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_smagorinsky_sfs_term(
  const int i,
  const int j,
  const int k,
  const amrex::Array4<const amrex::Real>& q,
  const amrex::Array4<const amrex::Real>& td,
  const amrex::Array4<const amrex::Real>& a,
  const amrex::Real dx,
  const int dir,
  const amrex::Real Cs,
  const amrex::Real CI,
  const amrex::Real PrT,
  const amrex::Array4<amrex::Real>& flx)
{

  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  const amrex::IntVect ivm(iv - amrex::IntVect::TheDimensionVector(dir));

  // SFS stress
  const amrex::Real dxinv = 1.0 / dx;
  const amrex::Real deltabar = dx;
  const amrex::Real Cs2 = Cs * Cs;
  const amrex::Real third = 1.0 / 3.0;
  amrex::Real alphaij[AMREX_SPACEDIM] = {0.0}, alpha, flux_T;
  if (dir == 0) {
    get_sfs_stresses_xdir(
      i, j, k, q, td, dxinv, deltabar, alphaij[0], alphaij[1], alphaij[2],
      alpha, flux_T);
  } else if (dir == 1) {
    get_sfs_stresses_ydir(
      i, j, k, q, td, dxinv, deltabar, alphaij[0], alphaij[1], alphaij[2],
      alpha, flux_T);
  } else {
    get_sfs_stresses_zdir(
      i, j, k, q, td, dxinv, deltabar, alphaij[0], alphaij[1], alphaij[2],
      alpha, flux_T);
  }
  const amrex::Real sigmadx =
    Cs2 * alphaij[0] - third * CI * alpha * static_cast<amrex::Real>(dir == 0);
  const amrex::Real sigmady =
    Cs2 * alphaij[1] - third * CI * alpha * static_cast<amrex::Real>(dir == 1);
  const amrex::Real sigmadz =
    Cs2 * alphaij[2] - third * CI * alpha * static_cast<amrex::Real>(dir == 2);
  flx(iv, UMX) = -sigmadx;
  flx(iv, UMY) = -sigmady;
  flx(iv, UMZ) = -sigmadz;

  // SFS turbulent diffusion
  const amrex::Real Uface[AMREX_SPACEDIM] = {
    static_cast<amrex::Real>(0.5 * (q(iv, QU) + q(ivm, QU))),
    static_cast<amrex::Real>(0.5 * (q(iv, QV) + q(ivm, QV))),
    static_cast<amrex::Real>(0.5 * (q(iv, QW) + q(ivm, QW)))};
  flx(iv, UEDEN) =
    -sigmadx * Uface[0] - sigmady * Uface[1] - sigmadz * Uface[2];

  // SFS heat flux - move state from cell centers to faces to compute cp for
  // flux at face
  amrex::Real T = 0.5 * (q(iv, QTEMP) + q(ivm, QTEMP));
  amrex::Real massfrac[NUM_SPECIES];
  for (int n = 0; n < NUM_SPECIES; n++) {
    massfrac[n] = 0.5 * (q(iv, QFS + n) + q(ivm, QFS + n));
  }
  amrex::Real cp = 0.0;
  auto eos = pele::physics::PhysicsType::eos();
  eos.TY2Cp(T, massfrac, cp);
  flx(iv, UEDEN) = flx(iv, UEDEN) - cp * Cs2 / PrT * flux_T;

  // Scale fluxes by area
  flx(iv, UMX) = flx(iv, UMX) * a(iv);
  flx(iv, UMY) = flx(iv, UMY) * a(iv);
  flx(iv, UMZ) = flx(iv, UMZ) * a(iv);
  flx(iv, UEDEN) = flx(iv, UEDEN) * a(iv);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
get_sfs_stresses_cc(
  const int i,
  const int j,
  const int k,
  const amrex::Array4<const amrex::Real>& q,
  const amrex::Real dxinv[AMREX_SPACEDIM],
  const amrex::Real deltabar,
  amrex::Real& alpha00,
  amrex::Real& alpha01,
  amrex::Real& alpha02,
  amrex::Real& alpha11,
  amrex::Real& alpha12,
  amrex::Real& alpha22,
  amrex::Real& alpha,
  amrex::Real& flux_T0,
  amrex::Real& flux_T1,
  amrex::Real& flux_T2)
{
  // Calculate derivatives at cell centers, second order central difference
  amrex::Real dUdx[3][3] = {{0.0}};
  dUdx[0][0] = 0.5 * dxinv[0] * (q(i + 1, j, k, QU) - q(i - 1, j, k, QU));
  dUdx[0][1] = 0.5 * dxinv[1] * (q(i, j + 1, k, QU) - q(i, j - 1, k, QU));
  dUdx[0][2] = 0.5 * dxinv[2] * (q(i, j, k + 1, QU) - q(i, j, k - 1, QU));
  dUdx[1][0] = 0.5 * dxinv[0] * (q(i + 1, j, k, QV) - q(i - 1, j, k, QV));
  dUdx[1][1] = 0.5 * dxinv[1] * (q(i, j + 1, k, QV) - q(i, j - 1, k, QV));
  dUdx[1][2] = 0.5 * dxinv[2] * (q(i, j, k + 1, QV) - q(i, j, k - 1, QV));
  dUdx[2][0] = 0.5 * dxinv[0] * (q(i + 1, j, k, QW) - q(i - 1, j, k, QW));
  dUdx[2][1] = 0.5 * dxinv[1] * (q(i, j + 1, k, QW) - q(i, j - 1, k, QW));
  dUdx[2][2] = 0.5 * dxinv[2] * (q(i, j, k + 1, QW) - q(i, j, k - 1, QW));

  amrex::Real S[3][3] = {{0.0}};
  amrex::Real Sijmag = 0;
  for (int m = 0; m < 3; m++) {
    for (int n = 0; n < 3; n++) {
      S[m][n] = 0.5 * (dUdx[m][n] + dUdx[n][m]);
      Sijmag += S[m][n] * S[m][n];
    }
  }
  Sijmag = std::sqrt(2.0 * Sijmag);
  const amrex::Real Skk = S[0][0] + S[1][1] + S[2][2];
  const amrex::Real mut = q(i, j, k, QRHO) * deltabar * deltabar * Sijmag;

  alpha00 = 2.0 * mut * (S[0][0] - Skk / 3.0);
  alpha01 = 2.0 * mut * (S[0][1]);
  alpha02 = 2.0 * mut * (S[0][2]);
  alpha11 = 2.0 * mut * (S[1][1] - Skk / 3.0);
  alpha12 = 2.0 * mut * (S[1][2]);
  alpha22 = 2.0 * mut * (S[2][2] - Skk / 3.0);
  alpha = 2.0 * mut * Sijmag;

  flux_T0 =
    mut * 0.5 * dxinv[0] * (q(i + 1, j, k, QTEMP) - q(i - 1, j, k, QTEMP));
  flux_T1 =
    mut * 0.5 * dxinv[1] * (q(i, j + 1, k, QTEMP) - q(i, j - 1, k, QTEMP));
  flux_T2 =
    mut * 0.5 * dxinv[2] * (q(i, j, k + 1, QTEMP) - q(i, j, k - 1, QTEMP));
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_dynamic_smagorinsky_quantities(
  const int i,
  const int j,
  const int k,
  const amrex::Array4<const amrex::Real>& q,
  const int fgr,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx,
  const amrex::Array4<amrex::Real>& Kij,
  const amrex::Array4<amrex::Real>& RUT,
  const amrex::Array4<amrex::Real>& alphaij,
  const amrex::Array4<amrex::Real>& alpha,
  const amrex::Array4<amrex::Real>& flux_T)
{
  // indexes
  const int i00 = 0 * AMREX_SPACEDIM + 0;
  const int i01 = 1 * AMREX_SPACEDIM + 0;
  const int i02 = 2 * AMREX_SPACEDIM + 0;
  const int i10 = 0 * AMREX_SPACEDIM + 1;
  const int i11 = 1 * AMREX_SPACEDIM + 1;
  const int i12 = 2 * AMREX_SPACEDIM + 1;
  const int i20 = 0 * AMREX_SPACEDIM + 2;
  const int i21 = 1 * AMREX_SPACEDIM + 2;
  const int i22 = 2 * AMREX_SPACEDIM + 2;

  // Maybe revisit how to determine deltahat in the case where dx1 != dx2 != dx3
  // (irrelevant for now as this is presently not allowed)
  const amrex::Real deltabar = fgr * (dx[0] + dx[1] + dx[2]) / 3.0;
  const amrex::Real dxinv[AMREX_SPACEDIM] = {
    static_cast<amrex::Real>(1.0 / dx[0]),
    static_cast<amrex::Real>(1.0 / dx[1]),
    static_cast<amrex::Real>(1.0 / dx[2])};

  get_sfs_stresses_cc(
    i, j, k, q, dxinv, deltabar, alphaij(i, j, k, i00), alphaij(i, j, k, i01),
    alphaij(i, j, k, i02), alphaij(i, j, k, i11), alphaij(i, j, k, i12),
    alphaij(i, j, k, i22), alpha(i, j, k), flux_T(i, j, k, 0),
    flux_T(i, j, k, 1), flux_T(i, j, k, 2));
  alphaij(i, j, k, i10) = alphaij(i, j, k, i01);
  alphaij(i, j, k, i20) = alphaij(i, j, k, i02);
  alphaij(i, j, k, i21) = alphaij(i, j, k, i12);

  Kij(i, j, k, 0) = q(i, j, k, QRHO) * q(i, j, k, QU) * q(i, j, k, QU);
  Kij(i, j, k, 1) = q(i, j, k, QRHO) * q(i, j, k, QU) * q(i, j, k, QV);
  Kij(i, j, k, 2) = q(i, j, k, QRHO) * q(i, j, k, QU) * q(i, j, k, QW);
  Kij(i, j, k, 3) = q(i, j, k, QRHO) * q(i, j, k, QV) * q(i, j, k, QV);
  Kij(i, j, k, 4) = q(i, j, k, QRHO) * q(i, j, k, QV) * q(i, j, k, QW);
  Kij(i, j, k, 5) = q(i, j, k, QRHO) * q(i, j, k, QW) * q(i, j, k, QW);

  RUT(i, j, k, 0) = q(i, j, k, QRHO) * q(i, j, k, QU) * q(i, j, k, QTEMP);
  RUT(i, j, k, 1) = q(i, j, k, QRHO) * q(i, j, k, QV) * q(i, j, k, QTEMP);
  RUT(i, j, k, 2) = q(i, j, k, QRHO) * q(i, j, k, QW) * q(i, j, k, QTEMP);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_dynamic_smagorinsky_coeffs(
  const int i,
  const int j,
  const int k,
  const amrex::Array4<const amrex::Real>& q,
  const int fgr,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx,
  const amrex::Array4<const amrex::Real>& Kij,
  const amrex::Array4<const amrex::Real>& RUT,
  const amrex::Array4<const amrex::Real>& alphaij,
  const amrex::Array4<const amrex::Real>& alpha,
  const amrex::Array4<const amrex::Real>& flux_T,
  const amrex::Array4<amrex::Real>& coeff)
{
  // indexes
  const int i00 = 0 * AMREX_SPACEDIM + 0;
  const int i01 = 1 * AMREX_SPACEDIM + 0;
  const int i02 = 2 * AMREX_SPACEDIM + 0;
  const int i10 = 0 * AMREX_SPACEDIM + 1;
  const int i11 = 1 * AMREX_SPACEDIM + 1;
  const int i12 = 2 * AMREX_SPACEDIM + 1;
  const int i20 = 0 * AMREX_SPACEDIM + 2;
  const int i21 = 1 * AMREX_SPACEDIM + 2;
  const int i22 = 2 * AMREX_SPACEDIM + 2;

  // Maybe revisit how to determine deltahat in the case where dx1 != dx2 != dx3
  // (irrelevant for now as this is presently not allowed)
  const amrex::Real deltahat = fgr * (dx[0] + dx[1] + dx[2]) / 3.0;
  const amrex::Real dxinv[AMREX_SPACEDIM] = {
    static_cast<amrex::Real>(1.0 / dx[0]),
    static_cast<amrex::Real>(1.0 / dx[1]),
    static_cast<amrex::Real>(1.0 / dx[2])};

  amrex::Real L[AMREX_SPACEDIM * AMREX_SPACEDIM],
    M[AMREX_SPACEDIM * AMREX_SPACEDIM], betaij[AMREX_SPACEDIM * AMREX_SPACEDIM],
    beta, T[AMREX_SPACEDIM], KE[AMREX_SPACEDIM];

  get_sfs_stresses_cc(
    i, j, k, q, dxinv, deltahat, betaij[i00], betaij[i01], betaij[i02],
    betaij[i11], betaij[i12], betaij[i22], beta, T[0], T[1], T[2]);
  betaij[i10] = betaij[i01];
  betaij[i20] = betaij[i02];
  betaij[i21] = betaij[i12];

  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    T[dir] = flux_T(i, j, k, dir) - T[dir];
  }

  // "resolved turbulent stresses" and others
  for (int dir = 0; dir < AMREX_SPACEDIM * AMREX_SPACEDIM; dir++) {
    M[dir] = betaij[dir] - alphaij(i, j, k, dir);
  }

  L[i00] = Kij(i, j, k, 0) - q(i, j, k, QRHO) * q(i, j, k, QU) * q(i, j, k, QU);
  L[i01] = Kij(i, j, k, 1) - q(i, j, k, QRHO) * q(i, j, k, QU) * q(i, j, k, QV);
  L[i02] = Kij(i, j, k, 2) - q(i, j, k, QRHO) * q(i, j, k, QU) * q(i, j, k, QW);
  L[i10] = L[i01];
  L[i11] = Kij(i, j, k, 3) - q(i, j, k, QRHO) * q(i, j, k, QV) * q(i, j, k, QV);
  L[i12] = Kij(i, j, k, 4) - q(i, j, k, QRHO) * q(i, j, k, QV) * q(i, j, k, QW);
  L[i20] = L[i02];
  L[i21] = L[i12];
  L[i22] = Kij(i, j, k, 5) - q(i, j, k, QRHO) * q(i, j, k, QW) * q(i, j, k, QW);

  KE[0] =
    RUT(i, j, k, 0) - q(i, j, k, QRHO) * q(i, j, k, QU) * q(i, j, k, QTEMP);
  KE[1] =
    RUT(i, j, k, 1) - q(i, j, k, QRHO) * q(i, j, k, QV) * q(i, j, k, QTEMP);
  KE[2] =
    RUT(i, j, k, 2) - q(i, j, k, QRHO) * q(i, j, k, QW) * q(i, j, k, QTEMP);

  // Contractions
  amrex::Real TT = 0.0, KT = 0.0, LM = 0.0, MM = 0.0;
  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    TT += T[dir] * T[dir];
    KT += KE[dir] * T[dir];
    for (int n = 0; n < AMREX_SPACEDIM; n++) {
      LM += L[dir * AMREX_SPACEDIM + n] * M[dir * AMREX_SPACEDIM + n];
      MM += M[dir * AMREX_SPACEDIM + n] * M[dir * AMREX_SPACEDIM + n];
    }
  }
  const amrex::Real Lkk = L[i00] + L[i11] + L[i22];
  const amrex::Real bma = beta - alpha(i, j, k);

  // Coefficients
  coeff(i, j, k, comp_Cs2) = amrex::max<amrex::Real>(
    LM / (MM + constants::small_num()), constants::small_num());
  coeff(i, j, k, comp_CI) = amrex::max<amrex::Real>(
    Lkk / (bma + constants::small_num()), constants::small_num());
  coeff(i, j, k, comp_Cs2ovPrT) = amrex::max<amrex::Real>(
    KT / (TT + constants::small_num()), constants::small_num());

  // Limit CI to reasonable bounds
  coeff(i, j, k, comp_CI) =
    amrex::min<amrex::Real>(coeff(i, j, k, comp_CI), 1.0);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_dynamic_smagorinsky_sfs_term(
  const int i,
  const int j,
  const int k,
  const amrex::Array4<const amrex::Real>& q,
  const amrex::Array4<const amrex::Real>& alphaij,
  const amrex::Array4<const amrex::Real>& alpha,
  const amrex::Array4<const amrex::Real>& flux_T,
  const amrex::Array4<const amrex::Real>& coeff,
  const amrex::Array4<const amrex::Real>& a,
  const int dir,
  const amrex::Array4<amrex::Real>& flx)
{
  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  const amrex::IntVect ivm(iv - amrex::IntVect::TheDimensionVector(dir));

  // SFS stress
  const amrex::Real sigmadx =
    coeff(iv, comp_Cs2) * alphaij(iv, 0) -
    coeff(iv, comp_CI) * alpha(iv) / 3.0 * static_cast<amrex::Real>(dir == 0);
  const amrex::Real sigmady =
    coeff(iv, comp_Cs2) * alphaij(iv, 1) -
    coeff(iv, comp_CI) * alpha(iv) / 3.0 * static_cast<amrex::Real>(dir == 1);
  const amrex::Real sigmadz =
    coeff(iv, comp_Cs2) * alphaij(iv, 2) -
    coeff(iv, comp_CI) * alpha(iv) / 3.0 * static_cast<amrex::Real>(dir == 2);
  flx(iv, UMX) = -sigmadx;
  flx(iv, UMY) = -sigmady;
  flx(iv, UMZ) = -sigmadz;

  // SFS turbulent diffusion
  const amrex::Real Uface[AMREX_SPACEDIM] = {
    static_cast<amrex::Real>(0.5 * (q(iv, QU) + q(ivm, QU))),
    static_cast<amrex::Real>(0.5 * (q(iv, QV) + q(ivm, QV))),
    static_cast<amrex::Real>(0.5 * (q(iv, QW) + q(ivm, QW)))};
  flx(iv, UEDEN) =
    -sigmadx * Uface[0] - sigmady * Uface[1] - sigmadz * Uface[2];

  // SFS heat flux - move state from cell centers to faces to compute cp for
  // flux at face
  amrex::Real T = 0.5 * (q(iv, QTEMP) + q(ivm, QTEMP));
  amrex::Real massfrac[NUM_SPECIES];
  for (int n = 0; n < NUM_SPECIES; n++) {
    massfrac[n] = 0.5 * (q(iv, QFS + n) + q(ivm, QFS + n));
  }
  amrex::Real cp = 0.0;
  auto eos = pele::physics::PhysicsType::eos();
  eos.TY2Cp(T, massfrac, cp);
  flx(iv, UEDEN) = flx(iv, UEDEN) - cp * coeff(iv, comp_Cs2ovPrT) * flux_T(iv);

  // Scale fluxes by area
  flx(iv, UMX) = flx(iv, UMX) * a(iv);
  flx(iv, UMY) = flx(iv, UMY) * a(iv);
  flx(iv, UMZ) = flx(iv, UMZ) * a(iv);
  flx(iv, UEDEN) = flx(iv, UEDEN) * a(iv);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
get_wale_sfs_stresses_xdir(
  const int i,
  const int j,
  const int k,
  const amrex::Array4<const amrex::Real>& q,
  const amrex::Array4<const amrex::Real>& td,
  const amrex::Real dxinv,
  const amrex::Real deltabar,
  amrex::Real& alphaij_xx,
  amrex::Real& alphaij_xy,
  amrex::Real& alphaij_xz,
  amrex::Real& alpha,
  amrex::Real& flux_T)
{

  amrex::Real g[3][3] = {{0.0}};
  g[0][0] = dxinv * (q(i, j, k, QU) - q(i - 1, j, k, QU)); // dudx
  g[0][1] = td(i, j, k, 0);                                // dudy
  g[0][2] = td(i, j, k, 3);                                // dudz
  g[1][0] = dxinv * (q(i, j, k, QV) - q(i - 1, j, k, QV)); // dvdx
  g[1][1] = td(i, j, k, 1);                                // dvdy
  g[1][2] = td(i, j, k, 4);                                // dvdz
  g[2][0] = dxinv * (q(i, j, k, QW) - q(i - 1, j, k, QW)); // dwdx
  g[2][1] = td(i, j, k, 2);                                // dwdy
  g[2][2] = td(i, j, k, 5);                                // dwdz

  amrex::Real g2[3][3] = {{0.0}};
  g2[0][0] = g[0][0] * g[0][0] + g[0][1] * g[1][0] + g[0][2] * g[2][0];
  g2[0][1] = g[0][0] * g[0][1] + g[0][1] * g[1][1] + g[0][2] * g[2][1];
  g2[0][2] = g[0][0] * g[0][2] + g[0][1] * g[1][2] + g[0][2] * g[2][2];
  g2[1][0] = g[1][0] * g[0][0] + g[1][1] * g[1][0] + g[1][2] * g[2][0];
  g2[1][1] = g[1][0] * g[0][1] + g[1][1] * g[1][1] + g[1][2] * g[2][1];
  g2[1][2] = g[1][0] * g[0][2] + g[1][1] * g[1][2] + g[1][2] * g[2][2];
  g2[2][0] = g[2][0] * g[0][0] + g[2][1] * g[1][0] + g[2][2] * g[2][0];
  g2[2][1] = g[2][0] * g[0][1] + g[2][1] * g[1][1] + g[2][2] * g[2][1];
  g2[2][2] = g[2][0] * g[0][2] + g[2][1] * g[1][2] + g[2][2] * g[2][2];

  amrex::Real g2kk = g2[0][0] + g2[1][1] + g2[2][2];

  amrex::Real S[3][3] = {{0.0}};
  amrex::Real Sijmag = 0;
  amrex::Real SijSij = 0;
  amrex::Real D[3][3] = {{0.0}};
  amrex::Real DijDij = 0;
  for (int m = 0; m < 3; m++) {
    for (int n = 0; n < 3; n++) {
      S[m][n] = 0.5 * (g[m][n] + g[n][m]);
      SijSij += S[m][n] * S[m][n];
      D[m][n] = 0.5 * (g2[m][n] + g2[n][m]);
      if (m == n) {
        D[m][n] = D[m][n] - g2kk / 3.0;
      }
      DijDij += D[m][n] * D[m][n];
    }
  }
  Sijmag = std::sqrt(2.0 * SijSij);
  const amrex::Real Skk = S[0][0] + S[1][1] + S[2][2];

  amrex::Real DijDij3 = DijDij * DijDij * DijDij;
  amrex::Real DijDij5 = DijDij3 * DijDij * DijDij;
  amrex::Real SijSij5 = SijSij * SijSij * SijSij * SijSij * SijSij;

  const amrex::Real mut =
    0.5 * (q(i, j, k, QRHO) + q(i - 1, j, k, QRHO)) * deltabar * deltabar *
    (std::sqrt(DijDij3) / (std::sqrt(SijSij5) + std::sqrt(std::sqrt(DijDij5))));

  alphaij_xx = 2.0 * mut * (S[0][0] - Skk / 3.0);
  alphaij_xy = 2.0 * mut * S[0][1];
  alphaij_xz = 2.0 * mut * S[0][2];
  alpha = 2.0 * 0.5 * (q(i, j, k, QRHO) + q(i - 1, j, k, QRHO)) * deltabar *
          deltabar * Sijmag * Sijmag;

  const amrex::Real dTdx = dxinv * (q(i, j, k, QTEMP) - q(i - 1, j, k, QTEMP));
  flux_T = mut * dTdx;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
get_wale_sfs_stresses_ydir(
  const int i,
  const int j,
  const int k,
  const amrex::Array4<const amrex::Real>& q,
  const amrex::Array4<const amrex::Real>& td,
  const amrex::Real dxinv,
  const amrex::Real deltabar,
  amrex::Real& alphaij_yx,
  amrex::Real& alphaij_yy,
  amrex::Real& alphaij_yz,
  amrex::Real& alpha,
  amrex::Real& flux_T)
{

  amrex::Real g[3][3] = {{0.0}};
  g[0][0] = td(i, j, k, 0);                                // dudx
  g[0][1] = dxinv * (q(i, j, k, QU) - q(i, j - 1, k, QU)); // dudy
  g[0][2] = td(i, j, k, 3);                                // dudz
  g[1][0] = td(i, j, k, 1);                                // dvdx
  g[1][1] = dxinv * (q(i, j, k, QV) - q(i, j - 1, k, QV)); // dvdy
  g[1][2] = td(i, j, k, 4);                                // dvdz
  g[2][0] = td(i, j, k, 2);                                // dwdx
  g[2][1] = dxinv * (q(i, j, k, QW) - q(i, j - 1, k, QW)); // dwdy
  g[2][2] = td(i, j, k, 5);                                // dwdz

  amrex::Real g2[3][3] = {{0.0}};
  g2[0][0] = g[0][0] * g[0][0] + g[0][1] * g[1][0] + g[0][2] * g[2][0];
  g2[0][1] = g[0][0] * g[0][1] + g[0][1] * g[1][1] + g[0][2] * g[2][1];
  g2[0][2] = g[0][0] * g[0][2] + g[0][1] * g[1][2] + g[0][2] * g[2][2];
  g2[1][0] = g[1][0] * g[0][0] + g[1][1] * g[1][0] + g[1][2] * g[2][0];
  g2[1][1] = g[1][0] * g[0][1] + g[1][1] * g[1][1] + g[1][2] * g[2][1];
  g2[1][2] = g[1][0] * g[0][2] + g[1][1] * g[1][2] + g[1][2] * g[2][2];
  g2[2][0] = g[2][0] * g[0][0] + g[2][1] * g[1][0] + g[2][2] * g[2][0];
  g2[2][1] = g[2][0] * g[0][1] + g[2][1] * g[1][1] + g[2][2] * g[2][1];
  g2[2][2] = g[2][0] * g[0][2] + g[2][1] * g[1][2] + g[2][2] * g[2][2];

  amrex::Real g2kk = g2[0][0] + g2[1][1] + g2[2][2];

  amrex::Real S[3][3] = {{0.0}};
  amrex::Real Sijmag = 0;
  amrex::Real SijSij = 0;
  amrex::Real D[3][3] = {{0.0}};
  amrex::Real DijDij = 0;
  for (int m = 0; m < 3; m++) {
    for (int n = 0; n < 3; n++) {
      S[m][n] = 0.5 * (g[m][n] + g[n][m]);
      SijSij += S[m][n] * S[m][n];
      D[m][n] = 0.5 * (g2[m][n] + g2[n][m]);
      if (m == n) {
        D[m][n] = D[m][n] - g2kk / 3.0;
      }
      DijDij += D[m][n] * D[m][n];
    }
  }
  Sijmag = std::sqrt(2.0 * SijSij);
  const amrex::Real Skk = S[0][0] + S[1][1] + S[2][2];

  amrex::Real DijDij3 = DijDij * DijDij * DijDij;
  amrex::Real DijDij5 = DijDij3 * DijDij * DijDij;
  amrex::Real SijSij5 = SijSij * SijSij * SijSij * SijSij * SijSij;

  const amrex::Real mut =
    0.5 * (q(i, j, k, QRHO) + q(i, j - 1, k, QRHO)) * deltabar * deltabar *
    (std::sqrt(DijDij3) / (std::sqrt(SijSij5) + std::sqrt(std::sqrt(DijDij5))));

  alphaij_yx = 2.0 * mut * S[1][0];
  alphaij_yy = 2.0 * mut * (S[1][1] - Skk / 3.0);
  alphaij_yz = 2.0 * mut * S[1][2];
  alpha = 2.0 * 0.5 * (q(i, j, k, QRHO) + q(i, j - 1, k, QRHO)) * deltabar *
          deltabar * Sijmag * Sijmag;

  const amrex::Real dTdy = dxinv * (q(i, j, k, QTEMP) - q(i, j - 1, k, QTEMP));
  flux_T = mut * dTdy;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
get_wale_sfs_stresses_zdir(
  const int i,
  const int j,
  const int k,
  const amrex::Array4<const amrex::Real>& q,
  const amrex::Array4<const amrex::Real>& td,
  const amrex::Real dxinv,
  const amrex::Real deltabar,
  amrex::Real& alphaij_zx,
  amrex::Real& alphaij_zy,
  amrex::Real& alphaij_zz,
  amrex::Real& alpha,
  amrex::Real& flux_T)
{

  amrex::Real g[3][3] = {{0.0}};
  g[0][0] = td(i, j, k, 0);                                // dudx
  g[0][1] = td(i, j, k, 3);                                // dudy
  g[0][2] = dxinv * (q(i, j, k, QU) - q(i, j, k - 1, QU)); // dudz
  g[1][0] = td(i, j, k, 1);                                // dvdx
  g[1][1] = td(i, j, k, 4);                                // dvdy
  g[1][2] = dxinv * (q(i, j, k, QV) - q(i, j, k - 1, QV)); // dvdz
  g[2][0] = td(i, j, k, 2);                                // dwdx
  g[2][1] = td(i, j, k, 5);                                // dwdy
  g[2][2] = dxinv * (q(i, j, k, QW) - q(i, j, k - 1, QW)); // dwdz

  amrex::Real g2[3][3] = {{0.0}};
  g2[0][0] = g[0][0] * g[0][0] + g[0][1] * g[1][0] + g[0][2] * g[2][0];
  g2[0][1] = g[0][0] * g[0][1] + g[0][1] * g[1][1] + g[0][2] * g[2][1];
  g2[0][2] = g[0][0] * g[0][2] + g[0][1] * g[1][2] + g[0][2] * g[2][2];
  g2[1][0] = g[1][0] * g[0][0] + g[1][1] * g[1][0] + g[1][2] * g[2][0];
  g2[1][1] = g[1][0] * g[0][1] + g[1][1] * g[1][1] + g[1][2] * g[2][1];
  g2[1][2] = g[1][0] * g[0][2] + g[1][1] * g[1][2] + g[1][2] * g[2][2];
  g2[2][0] = g[2][0] * g[0][0] + g[2][1] * g[1][0] + g[2][2] * g[2][0];
  g2[2][1] = g[2][0] * g[0][1] + g[2][1] * g[1][1] + g[2][2] * g[2][1];
  g2[2][2] = g[2][0] * g[0][2] + g[2][1] * g[1][2] + g[2][2] * g[2][2];

  amrex::Real g2kk = g2[0][0] + g2[1][1] + g2[2][2];

  amrex::Real S[3][3] = {{0.0}};
  amrex::Real Sijmag = 0;
  amrex::Real SijSij = 0;
  amrex::Real D[3][3] = {{0.0}};
  amrex::Real DijDij = 0;
  for (int m = 0; m < 3; m++) {
    for (int n = 0; n < 3; n++) {
      S[m][n] = 0.5 * (g[m][n] + g[n][m]);
      SijSij += S[m][n] * S[m][n];
      D[m][n] = 0.5 * (g2[m][n] + g2[n][m]);
      if (m == n) {
        D[m][n] = D[m][n] - g2kk / 3.0;
      }
      DijDij += D[m][n] * D[m][n];
    }
  }
  Sijmag = std::sqrt(2.0 * SijSij);
  const amrex::Real Skk = S[0][0] + S[1][1] + S[2][2];

  amrex::Real DijDij3 = DijDij * DijDij * DijDij;
  amrex::Real DijDij5 = DijDij3 * DijDij * DijDij;
  amrex::Real SijSij5 = SijSij * SijSij * SijSij * SijSij * SijSij;

  const amrex::Real mut =
    0.5 * (q(i, j, k, QRHO) + q(i, j, k - 1, QRHO)) * deltabar * deltabar *
    (std::sqrt(DijDij3) / (std::sqrt(SijSij5) + std::sqrt(std::sqrt(DijDij5))));

  alphaij_zx = 2.0 * mut * S[2][0];
  alphaij_zy = 2.0 * mut * S[2][1];
  alphaij_zz = 2.0 * mut * (S[2][2] - Skk / 3.0);
  alpha = 2.0 * 0.5 * (q(i, j, k, QRHO) + q(i, j, k - 1, QRHO)) * deltabar *
          deltabar * Sijmag * Sijmag;

  const amrex::Real dTdz = dxinv * (q(i, j, k, QTEMP) - q(i, j, k - 1, QTEMP));
  flux_T = mut * dTdz;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_wale_sfs_term(
  const int i,
  const int j,
  const int k,
  const amrex::Array4<const amrex::Real>& q,
  const amrex::Array4<const amrex::Real>& td,
  const amrex::Array4<const amrex::Real>& a,
  const amrex::Real dx,
  const int dir,
  const amrex::Real Cw,
  const amrex::Real CI,
  const amrex::Real PrT,
  const amrex::Array4<amrex::Real>& flx)
{

  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  const amrex::IntVect ivm(iv - amrex::IntVect::TheDimensionVector(dir));

  // SFS stress
  const amrex::Real dxinv = 1.0 / dx;
  const amrex::Real deltabar = dx;
  const amrex::Real Cw2 = Cw * Cw;
  const amrex::Real third = 1.0 / 3.0;
  amrex::Real alphaij[AMREX_SPACEDIM] = {0.0}, alpha, flux_T;
  if (dir == 0) {
    get_wale_sfs_stresses_xdir(
      i, j, k, q, td, dxinv, deltabar, alphaij[0], alphaij[1], alphaij[2],
      alpha, flux_T);
  } else if (dir == 1) {
    get_wale_sfs_stresses_ydir(
      i, j, k, q, td, dxinv, deltabar, alphaij[0], alphaij[1], alphaij[2],
      alpha, flux_T);
  } else {
    get_wale_sfs_stresses_zdir(
      i, j, k, q, td, dxinv, deltabar, alphaij[0], alphaij[1], alphaij[2],
      alpha, flux_T);
  }
  const amrex::Real sigmadx =
    Cw2 * alphaij[0] - third * CI * alpha * static_cast<amrex::Real>(dir == 0);
  const amrex::Real sigmady =
    Cw2 * alphaij[1] - third * CI * alpha * static_cast<amrex::Real>(dir == 1);
  const amrex::Real sigmadz =
    Cw2 * alphaij[2] - third * CI * alpha * static_cast<amrex::Real>(dir == 2);
  flx(iv, UMX) = -sigmadx;
  flx(iv, UMY) = -sigmady;
  flx(iv, UMZ) = -sigmadz;

  // SFS turbulent diffusion
  const amrex::Real Uface[AMREX_SPACEDIM] = {
    static_cast<amrex::Real>(0.5 * (q(iv, QU) + q(ivm, QU))),
    static_cast<amrex::Real>(0.5 * (q(iv, QV) + q(ivm, QV))),
    static_cast<amrex::Real>(0.5 * (q(iv, QW) + q(ivm, QW)))};
  flx(iv, UEDEN) =
    -sigmadx * Uface[0] - sigmady * Uface[1] - sigmadz * Uface[2];

  // SFS heat flux - move state from cell centers to faces to compute cp for
  // flux at face
  amrex::Real T = 0.5 * (q(iv, QTEMP) + q(ivm, QTEMP));
  amrex::Real massfrac[NUM_SPECIES];
  for (int n = 0; n < NUM_SPECIES; n++) {
    massfrac[n] = 0.5 * (q(iv, QFS + n) + q(ivm, QFS + n));
  }
  amrex::Real cp = 0.0;
  auto eos = pele::physics::PhysicsType::eos();
  eos.TY2Cp(T, massfrac, cp);
  flx(iv, UEDEN) = flx(iv, UEDEN) - cp * Cw2 / PrT * flux_T;

  // Scale fluxes by area
  flx(iv, UMX) = flx(iv, UMX) * a(iv);
  flx(iv, UMY) = flx(iv, UMY) * a(iv);
  flx(iv, UMZ) = flx(iv, UMZ) * a(iv);
  flx(iv, UEDEN) = flx(iv, UEDEN) * a(iv);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
get_vreman_sfs_stresses_xdir(
  const int i,
  const int j,
  const int k,
  const amrex::Array4<const amrex::Real>& q,
  const amrex::Array4<const amrex::Real>& td,
  const amrex::Real dxinv,
  const amrex::Real deltabar,
  amrex::Real& alphaij_xx,
  amrex::Real& alphaij_xy,
  amrex::Real& alphaij_xz,
  amrex::Real& alpha,
  amrex::Real& flux_T)
{

  amrex::Real alfa[3][3] = {{0.0}};
  alfa[0][0] = dxinv * (q(i, j, k, QU) - q(i - 1, j, k, QU)); // dudx  d1v1
  alfa[0][1] = td(i, j, k, 0);                                // dudy  d2v1
  alfa[0][2] = td(i, j, k, 3);                                // dudz  d3v1
  alfa[1][0] = dxinv * (q(i, j, k, QV) - q(i - 1, j, k, QV)); // dvdx  d1v2
  alfa[1][1] = td(i, j, k, 1);                                // dvdy  d2v2
  alfa[1][2] = td(i, j, k, 4);                                // dvdz  d3v2
  alfa[2][0] = dxinv * (q(i, j, k, QW) - q(i - 1, j, k, QW)); // dwdx  d1v3
  alfa[2][1] = td(i, j, k, 2);                                // dwdy  d2v3
  alfa[2][2] = td(i, j, k, 5);                                // dwdz  d3v3

  amrex::Real S[3][3] = {{0.0}};
  amrex::Real beta[3][3] = {{0.0}};
  amrex::Real Sijmag = 0;
  amrex::Real SijSij = 0;
  amrex::Real aijaij = 0;
  for (int m = 0; m < 3; m++) {
    for (int n = 0; n < 3; n++) {
      S[m][n] = 0.5 * (alfa[m][n] + alfa[n][m]);
      SijSij += S[m][n] * S[m][n];
      aijaij += alfa[m][n] * alfa[m][n];
      beta[m][n] = deltabar * deltabar *
                   (alfa[0][m] * alfa[0][n] + alfa[1][m] * alfa[1][n] +
                    alfa[2][m] * alfa[2][n]);
    }
  }
  amrex::Real bbeta = 0;
  bbeta = beta[0][0] * beta[1][1] - beta[0][1] * beta[0][1] +
          beta[0][0] * beta[2][2] - beta[0][2] * beta[0][2] +
          beta[1][1] * beta[2][2] - beta[1][2] * beta[1][2];

  Sijmag = std::sqrt(2.0 * SijSij);
  const amrex::Real Skk = S[0][0] + S[1][1] + S[2][2];

  const amrex::Real mut = 0.5 * (q(i, j, k, QRHO) + q(i - 1, j, k, QRHO)) *
                          (std::sqrt(bbeta / aijaij));

  alphaij_xx = 2.0 * mut * (S[0][0] - Skk / 3.0);
  alphaij_xy = 2.0 * mut * S[0][1];
  alphaij_xz = 2.0 * mut * S[0][2];
  alpha = 2.0 * 0.5 * (q(i, j, k, QRHO) + q(i - 1, j, k, QRHO)) * deltabar *
          deltabar * Sijmag * Sijmag;

  const amrex::Real dTdx = dxinv * (q(i, j, k, QTEMP) - q(i - 1, j, k, QTEMP));
  flux_T = mut * dTdx;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
get_vreman_sfs_stresses_ydir(
  const int i,
  const int j,
  const int k,
  const amrex::Array4<const amrex::Real>& q,
  const amrex::Array4<const amrex::Real>& td,
  const amrex::Real dxinv,
  const amrex::Real deltabar,
  amrex::Real& alphaij_yx,
  amrex::Real& alphaij_yy,
  amrex::Real& alphaij_yz,
  amrex::Real& alpha,
  amrex::Real& flux_T)
{

  amrex::Real alfa[3][3] = {{0.0}};
  alfa[0][0] = td(i, j, k, 0);                                // dudx
  alfa[0][1] = dxinv * (q(i, j, k, QU) - q(i, j - 1, k, QU)); // dudy
  alfa[0][2] = td(i, j, k, 3);                                // dudz
  alfa[1][0] = td(i, j, k, 1);                                // dvdx
  alfa[1][1] = dxinv * (q(i, j, k, QV) - q(i, j - 1, k, QV)); // dvdy
  alfa[1][2] = td(i, j, k, 4);                                // dvdz
  alfa[2][0] = td(i, j, k, 2);                                // dwdx
  alfa[2][1] = dxinv * (q(i, j, k, QW) - q(i, j - 1, k, QW)); // dwdy
  alfa[2][2] = td(i, j, k, 5);                                // dwdz

  amrex::Real S[3][3] = {{0.0}};
  amrex::Real beta[3][3] = {{0.0}};
  amrex::Real Sijmag = 0;
  amrex::Real SijSij = 0;
  amrex::Real aijaij = 0;
  for (int m = 0; m < 3; m++) {
    for (int n = 0; n < 3; n++) {
      S[m][n] = 0.5 * (alfa[m][n] + alfa[n][m]);
      SijSij += S[m][n] * S[m][n];
      aijaij += alfa[m][n] * alfa[m][n];
      beta[m][n] = deltabar * deltabar *
                   (alfa[0][m] * alfa[0][n] + alfa[1][m] * alfa[1][n] +
                    alfa[2][m] * alfa[2][n]);
    }
  }
  amrex::Real bbeta = 0;
  bbeta = beta[0][0] * beta[1][1] - beta[0][1] * beta[0][1] +
          beta[0][0] * beta[2][2] - beta[0][2] * beta[0][2] +
          beta[1][1] * beta[2][2] - beta[1][2] * beta[1][2];

  Sijmag = std::sqrt(2.0 * SijSij);
  const amrex::Real Skk = S[0][0] + S[1][1] + S[2][2];

  const amrex::Real mut = 0.5 * (q(i, j, k, QRHO) + q(i - 1, j, k, QRHO)) *
                          (std::sqrt(bbeta / aijaij));

  alphaij_yx = 2.0 * mut * S[1][0];
  alphaij_yy = 2.0 * mut * (S[1][1] - Skk / 3.0);
  alphaij_yz = 2.0 * mut * S[1][2];
  alpha = 2.0 * 0.5 * (q(i, j, k, QRHO) + q(i, j - 1, k, QRHO)) * deltabar *
          deltabar * Sijmag * Sijmag;

  const amrex::Real dTdy = dxinv * (q(i, j, k, QTEMP) - q(i, j - 1, k, QTEMP));
  flux_T = mut * dTdy;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
get_vreman_sfs_stresses_zdir(
  const int i,
  const int j,
  const int k,
  const amrex::Array4<const amrex::Real>& q,
  const amrex::Array4<const amrex::Real>& td,
  const amrex::Real dxinv,
  const amrex::Real deltabar,
  amrex::Real& alphaij_zx,
  amrex::Real& alphaij_zy,
  amrex::Real& alphaij_zz,
  amrex::Real& alpha,
  amrex::Real& flux_T)
{

  amrex::Real alfa[3][3] = {{0.0}};
  alfa[0][0] = td(i, j, k, 0);                                // dudx
  alfa[0][1] = td(i, j, k, 3);                                // dudy
  alfa[0][2] = dxinv * (q(i, j, k, QU) - q(i, j, k - 1, QU)); // dudz
  alfa[1][0] = td(i, j, k, 1);                                // dvdx
  alfa[1][1] = td(i, j, k, 4);                                // dvdy
  alfa[1][2] = dxinv * (q(i, j, k, QV) - q(i, j, k - 1, QV)); // dvdz
  alfa[2][0] = td(i, j, k, 2);                                // dwdx
  alfa[2][1] = td(i, j, k, 5);                                // dwdy
  alfa[2][2] = dxinv * (q(i, j, k, QW) - q(i, j, k - 1, QW)); // dwdz

  amrex::Real S[3][3] = {{0.0}};
  amrex::Real beta[3][3] = {{0.0}};
  amrex::Real Sijmag = 0;
  amrex::Real SijSij = 0;
  amrex::Real aijaij = 0;
  for (int m = 0; m < 3; m++) {
    for (int n = 0; n < 3; n++) {
      S[m][n] = 0.5 * (alfa[m][n] + alfa[n][m]);
      SijSij += S[m][n] * S[m][n];
      aijaij += alfa[m][n] * alfa[m][n];
      beta[m][n] = deltabar * deltabar *
                   (alfa[0][m] * alfa[0][n] + alfa[1][m] * alfa[1][n] +
                    alfa[2][m] * alfa[2][n]);
    }
  }
  amrex::Real bbeta = 0;
  bbeta = beta[0][0] * beta[1][1] - beta[0][1] * beta[0][1] +
          beta[0][0] * beta[2][2] - beta[0][2] * beta[0][2] +
          beta[1][1] * beta[2][2] - beta[1][2] * beta[1][2];

  Sijmag = std::sqrt(2.0 * SijSij);
  const amrex::Real Skk = S[0][0] + S[1][1] + S[2][2];

  const amrex::Real mut = 0.5 * (q(i, j, k, QRHO) + q(i - 1, j, k, QRHO)) *
                          (std::sqrt(bbeta / aijaij));

  alphaij_zx = 2.0 * mut * S[2][0];
  alphaij_zy = 2.0 * mut * S[2][1];
  alphaij_zz = 2.0 * mut * (S[2][2] - Skk / 3.0);
  alpha = 2.0 * 0.5 * (q(i, j, k, QRHO) + q(i, j, k - 1, QRHO)) * deltabar *
          deltabar * Sijmag * Sijmag;

  const amrex::Real dTdz = dxinv * (q(i, j, k, QTEMP) - q(i, j, k - 1, QTEMP));
  flux_T = mut * dTdz;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_vreman_sfs_term(
  const int i,
  const int j,
  const int k,
  const amrex::Array4<const amrex::Real>& q,
  const amrex::Array4<const amrex::Real>& td,
  const amrex::Array4<const amrex::Real>& a,
  const amrex::Real dx,
  const int dir,
  const amrex::Real Cw,
  const amrex::Real CI,
  const amrex::Real PrT,
  const amrex::Array4<amrex::Real>& flx)
{

  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  const amrex::IntVect ivm(iv - amrex::IntVect::TheDimensionVector(dir));

  // SFS stress
  const amrex::Real dxinv = 1.0 / dx;
  const amrex::Real deltabar = dx;
  const amrex::Real Cw2 = Cw * Cw;
  const amrex::Real third = 1.0 / 3.0;
  amrex::Real alphaij[AMREX_SPACEDIM] = {0.0}, alpha, flux_T;
  if (dir == 0) {
    get_vreman_sfs_stresses_xdir(
      i, j, k, q, td, dxinv, deltabar, alphaij[0], alphaij[1], alphaij[2],
      alpha, flux_T);
  } else if (dir == 1) {
    get_vreman_sfs_stresses_ydir(
      i, j, k, q, td, dxinv, deltabar, alphaij[0], alphaij[1], alphaij[2],
      alpha, flux_T);
  } else {
    get_vreman_sfs_stresses_zdir(
      i, j, k, q, td, dxinv, deltabar, alphaij[0], alphaij[1], alphaij[2],
      alpha, flux_T);
  }
  const amrex::Real sigmadx =
    Cw2 * alphaij[0] - third * CI * alpha * static_cast<amrex::Real>(dir == 0);
  const amrex::Real sigmady =
    Cw2 * alphaij[1] - third * CI * alpha * static_cast<amrex::Real>(dir == 1);
  const amrex::Real sigmadz =
    Cw2 * alphaij[2] - third * CI * alpha * static_cast<amrex::Real>(dir == 2);
  flx(iv, UMX) = -sigmadx;
  flx(iv, UMY) = -sigmady;
  flx(iv, UMZ) = -sigmadz;

  // SFS turbulent diffusion
  const amrex::Real Uface[AMREX_SPACEDIM] = {
    static_cast<amrex::Real>(0.5 * (q(iv, QU) + q(ivm, QU))),
    static_cast<amrex::Real>(0.5 * (q(iv, QV) + q(ivm, QV))),
    static_cast<amrex::Real>(0.5 * (q(iv, QW) + q(ivm, QW)))};
  flx(iv, UEDEN) =
    -sigmadx * Uface[0] - sigmady * Uface[1] - sigmadz * Uface[2];

  // SFS heat flux - move state from cell centers to faces to compute cp for
  // flux at face
  amrex::Real T = 0.5 * (q(iv, QTEMP) + q(ivm, QTEMP));
  amrex::Real massfrac[NUM_SPECIES];
  for (int n = 0; n < NUM_SPECIES; n++) {
    massfrac[n] = 0.5 * (q(iv, QFS + n) + q(ivm, QFS + n));
  }
  amrex::Real cp = 0.0;
  auto eos = pele::physics::PhysicsType::eos();
  eos.TY2Cp(T, massfrac, cp);
  flx(iv, UEDEN) = flx(iv, UEDEN) - cp * Cw2 / PrT * flux_T;

  // Scale fluxes by area
  flx(iv, UMX) = flx(iv, UMX) * a(iv);
  flx(iv, UMY) = flx(iv, UMY) * a(iv);
  flx(iv, UMZ) = flx(iv, UMZ) * a(iv);
  flx(iv, UEDEN) = flx(iv, UEDEN) * a(iv);
}

#endif // AMREX_SPACEDIM == 3
#endif
