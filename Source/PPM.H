#ifndef PPM_H
#define PPM_H

#include <AMReX_FArrayBox.H>
#include <AMReX_IArrayBox.H>

#include "Constants.H"
#include "IndexDefines.H"
#include "PelePhysics.H"
#include "Riemann.H"

constexpr int im2 = 0;
constexpr int im1 = 1;
constexpr int i0 = 2;
constexpr int ip1 = 3;
constexpr int ip2 = 4;

// Compute the coefficients of a parabolic reconstruction of the data in a
// zone. This uses the standard PPM limiters described in Colella & Woodward
// (1984)
//
// @param s      Real[5] the state to be reconstructed in zones i-2, i-1, i,
// i+1, i+2
// @param flatn  The flattening coefficient
// @param sm     The value of the parabola on the left edge of the zone
// @param sp     The value of the parabola on the right edge of the zone
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
ppm_reconstruct(
  const amrex::Real* s,
  const amrex::Real flatn,
  amrex::Real& sm,
  amrex::Real& sp)
{
  // Compute van Leer slopes

  amrex::Real dsl = 2.0 * (s[im1] - s[im2]);
  amrex::Real dsr = 2.0 * (s[i0] - s[im1]);

  amrex::Real dsvl_l = 0.0;
  if (dsl * dsr > 0.0) {
    amrex::Real dsc = 0.5 * (s[i0] - s[im2]);
    dsvl_l =
      amrex::Math::copysign(1.0, dsc) *
      amrex::min<amrex::Real>(
        amrex::Math::abs(dsc),
        amrex::min<amrex::Real>(amrex::Math::abs(dsl), amrex::Math::abs(dsr)));
  }

  dsl = 2.0 * (s[i0] - s[im1]);
  dsr = 2.0 * (s[ip1] - s[i0]);

  amrex::Real dsvl_r = 0.0;
  if (dsl * dsr > 0.0) {
    amrex::Real dsc = 0.5 * (s[ip1] - s[im1]);
    dsvl_r =
      amrex::Math::copysign(1.0, dsc) *
      amrex::min<amrex::Real>(
        amrex::Math::abs(dsc),
        amrex::min<amrex::Real>(amrex::Math::abs(dsl), amrex::Math::abs(dsr)));
  }

  // Interpolate s to edges
  sm = 0.5 * (s[i0] + s[im1]) - (1.0 / 6.0) * (dsvl_r - dsvl_l);

  // Make sure sedge lies in between adjacent cell-centered values
  sm = amrex::max<amrex::Real>(sm, amrex::min<amrex::Real>(s[i0], s[im1]));
  sm = amrex::min<amrex::Real>(sm, amrex::max<amrex::Real>(s[i0], s[im1]));

  // Compute van Leer slopes
  dsl = 2.0 * (s[i0] - s[im1]);
  dsr = 2.0 * (s[ip1] - s[i0]);

  dsvl_l = 0.0;
  if (dsl * dsr > 0.0) {
    amrex::Real dsc = 0.5 * (s[ip1] - s[im1]);
    dsvl_l =
      amrex::Math::copysign(1.0, dsc) *
      amrex::min<amrex::Real>(
        amrex::Math::abs(dsc),
        amrex::min<amrex::Real>(amrex::Math::abs(dsl), amrex::Math::abs(dsr)));
  }

  dsl = 2.0 * (s[ip1] - s[i0]);
  dsr = 2.0 * (s[ip2] - s[ip1]);

  dsvl_r = 0.0;
  if (dsl * dsr > 0.0) {
    amrex::Real dsc = 0.5 * (s[ip2] - s[i0]);
    dsvl_r =
      amrex::Math::copysign(1.0, dsc) *
      amrex::min<amrex::Real>(
        amrex::Math::abs(dsc),
        amrex::min<amrex::Real>(amrex::Math::abs(dsl), amrex::Math::abs(dsr)));
  }

  // Interpolate s to edges
  sp = 0.5 * (s[ip1] + s[i0]) - (1.0 / 6.0) * (dsvl_r - dsvl_l);

  // Make sure sedge lies in between adjacent cell-centered values
  sp = amrex::max<amrex::Real>(sp, amrex::min<amrex::Real>(s[ip1], s[i0]));
  sp = amrex::min<amrex::Real>(sp, amrex::max<amrex::Real>(s[ip1], s[i0]));

  // Flatten the parabola
  sm = flatn * sm + (1.0 - flatn) * s[i0];
  sp = flatn * sp + (1.0 - flatn) * s[i0];

  // Modify using quadratic limiters -- note this version of the limiting comes
  // from Colella and Sekora (2008), not the original PPM paper.
  if ((sp - s[i0]) * (s[i0] - sm) <= 0.0) {
    sp = s[i0];
    sm = s[i0];

  } else if (
    amrex::Math::abs(sp - s[i0]) >= 2.0 * amrex::Math::abs(sm - s[i0])) {
    sp = 3.0 * s[i0] - 2.0 * sm;

  } else if (
    amrex::Math::abs(sm - s[i0]) >= 2.0 * amrex::Math::abs(sp - s[i0])) {
    sm = 3.0 * s[i0] - 2.0 * sp;
  }
}

// Integrate under the parabola using from the left and right edges
// with the wave speeds u-c, u, u+c
//
// @param sm    The parabola value on the left edge of the zone
// @param sp    The parabola value on the right edge of the zone
// @param sc    The state value at the zone center
// @param u     The fluid velocity
// @param c     The sound speed
// @param dtdx  dt/dx (timestep / cell width)
// @param Ip    amrex::Real[3], the integrals under the parabola over
//              the 3 waves from the right zone edge
// @param Im    amrex::Real[3], the integrals under the parabola over
//              the 3 waves from the left zone edge
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
ppm_int_profile(
  const amrex::Real sm,
  const amrex::Real sp,
  const amrex::Real sc,
  const amrex::Real u,
  const amrex::Real c,
  const amrex::Real dtdx,
  amrex::Real* Ip,
  amrex::Real* Im)
{

  // Integrate the parabolic profile to the edge of the cell.

  // compute x-component of Ip and Im
  amrex::Real s6 = 6.0 * sc - 3.0 * (sm + sp);

  // Ip/m is the integral under the parabola for the extent
  // that a wave can travel over a timestep

  // Ip integrates to the right edge of a cell
  // Im integrates to the left edge of a cell

  // u-c wave
  amrex::Real speed = u - c;
  amrex::Real sigma = amrex::Math::abs(speed) * dtdx;

  // if speed == ZERO, then either branch is the same
  if (speed <= 0.0) {
    Ip[0] = sp;
    Im[0] = sm + 0.5 * sigma * (sp - sm + (1.0 - (2.0 / 3.0) * sigma) * s6);
  } else {
    Ip[0] = sp - 0.5 * sigma * (sp - sm - (1.0 - (2.0 / 3.0) * sigma) * s6);
    Im[0] = sm;
  }

  // u wave
  speed = u;
  sigma = amrex::Math::abs(speed) * dtdx;

  if (speed <= 0.0) {
    Ip[1] = sp;
    Im[1] = sm + 0.5 * sigma * (sp - sm + (1.0 - (2.0 / 3.0) * sigma) * s6);
  } else {
    Ip[1] = sp - 0.5 * sigma * (sp - sm - (1.0 - (2.0 / 3.0) * sigma) * s6);
    Im[1] = sm;
  }

  // u+c wave
  speed = u + c;
  sigma = amrex::Math::abs(speed) * dtdx;

  if (speed <= 0.0) {
    Ip[2] = sp;
    Im[2] = sm + 0.5 * sigma * (sp - sm + (1.0 - (2.0 / 3.0) * sigma) * s6);
  } else {
    Ip[2] = sp - 0.5 * sigma * (sp - sm - (1.0 - (2.0 / 3.0) * sigma) * s6);
    Im[2] = sm;
  }
}

// Integrate under the parabola from the left and right edges for a
// single wave speed

// @param sm    The parabola value on the left edge of the zone
// @param sp    The parabola value on the right edge of the zone
// @param sc    The state value at the zone center
// @param lam   The wave speed
// @param dtdx  dt/dx (timestep / cell width)
// @param Ip    integral under the parabola from the right zone edge
// @param Im    integral under the parabola from the left zone edge
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
ppm_int_profile_single(
  const amrex::Real sm,
  const amrex::Real sp,
  const amrex::Real sc,
  const amrex::Real lam,
  const amrex::Real dtdx,
  amrex::Real& Ip,
  amrex::Real& Im)
{

  // Integrate the parabolic profile to the edge of the cell.
  // This is the MHD version.  We come in already with the eigenvalues.

  // compute x-component of Ip and Im
  amrex::Real s6 = 6.0 * sc - 3.0 * (sm + sp);

  // Ip/m is the integral under the parabola for the extent
  // that a wave can travel over a timestep

  // Ip integrates to the right edge of a cell
  // Im integrates to the left edge of a cell

  amrex::Real sigma = amrex::Math::abs(lam) * dtdx;

  if (lam <= 0.0) {
    Ip = sp;
    Im = sm + 0.5 * sigma * (sp - sm + (1.0 - (2.0 / 3.0) * sigma) * s6);
  } else {
    Ip = sp - 0.5 * sigma * (sp - sm - (1.0 - (2.0 / 3.0) * sigma) * s6);
    Im = sm;
  }
}

// Host Functions
void trace_ppm(
  const amrex::Box& bx,
  const int idir,
  amrex::Array4<amrex::Real const> const& q_arr,
  amrex::Array4<amrex::Real const> const& srcQ,
  amrex::Array4<amrex::Real> const& qm,
  amrex::Array4<amrex::Real> const& qp,
  const amrex::Box& vbx,
  const amrex::Real dt,
  const amrex::Real* dx,
  const bool use_flattening,
  const bool use_hybrid_weno,
  const int weno_scheme);

#endif
