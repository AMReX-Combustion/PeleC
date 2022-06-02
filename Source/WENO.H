#ifndef WENO_H
#define WENO_H

#include <AMReX_FArrayBox.H>
#include <AMReX_IArrayBox.H>

// WENO 5 JS (Original Jiang GS, Shu CW. Efficient Implementation of Weighted
// ENO Schemes.
//             J. Comput. Phys. 1996; 126(1): 202–228.
//             doi: 10.1006/jcph.1996.0130)
// @param s      Stencils i-2, i-1, i, i+1, i+2
// @param sm     The value of the interpolation on the left edge of the i cell
// @param sp     The value of the interpolation on the right edge of the i cell
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
weno_reconstruct_5js(const amrex::Real* s, amrex::Real& sm, amrex::Real& sp)
{
  constexpr amrex::Real eps = 1e-6;

  amrex::Real vl[3];
  amrex::Real vr[3];
  amrex::Real beta[3];
  amrex::Real alpha[3];

  beta[2] = (13.0 / 12.0) * pow(s[0] - 2.0 * s[1] + s[2], 2) +
            0.25 * pow(s[0] - 4.0 * s[1] + 3.0 * s[2], 2);
  beta[1] = (13.0 / 12.0) * pow(s[1] - 2.0 * s[2] + s[3], 2) +
            0.25 * (s[1] - s[3]) * (s[1] - s[3]);
  beta[0] = (13.0 / 12.0) * pow(s[2] - 2.0 * s[3] + s[4], 2) +
            0.25 * pow(3.0 * s[2] - 4.0 * s[3] + s[4], 2);

  beta[2] = 1.0 / ((eps + beta[2]) * (eps + beta[2]));
  beta[1] = 1.0 / ((eps + beta[1]) * (eps + beta[1]));
  beta[0] = 1.0 / ((eps + beta[0]) * (eps + beta[0]));

  alpha[2] = beta[2];
  alpha[1] = 6.0 * beta[1];
  alpha[0] = 3.0 * beta[0];
  amrex::Real alpha1 = 1.0 / (alpha[2] + alpha[1] + alpha[0]);

  vl[2] = 2.0 * s[0] - 7.0 * s[1] + 11.0 * s[2];
  vl[1] = -s[1] + 5.0 * s[2] + 2.0 * s[3];
  vl[0] = 2.0 * s[2] + 5.0 * s[3] - s[4];

  sp = (1.0 / 6.0) * alpha1 *
       (alpha[2] * vl[2] + alpha[1] * vl[1] + alpha[0] * vl[0]);

  beta[2] = (13.0 / 12.0) * pow(s[4] - 2.0 * s[3] + s[2], 2) +
            0.25 * pow(s[4] - 4.0 * s[3] + 3.0 * s[2], 2);
  beta[1] = (13.0 / 12.0) * pow(s[3] - 2.0 * s[2] + s[1], 2) +
            0.25 * (s[3] - s[1]) * (s[3] - s[1]);
  beta[0] = (13.0 / 12.0) * pow(s[2] - 2.0 * s[1] + s[0], 2) +
            0.25 * pow(3.0 * s[2] - 4.0 * s[1] + s[0], 2);

  beta[2] = 1.0 / ((eps + beta[2]) * (eps + beta[2]));
  beta[1] = 1.0 / ((eps + beta[1]) * (eps + beta[1]));
  beta[0] = 1.0 / ((eps + beta[0]) * (eps + beta[0]));

  alpha[2] = beta[2];
  alpha[1] = 6.0 * beta[1];
  alpha[0] = 3.0 * beta[0];
  alpha1 = 1.0 / (alpha[2] + alpha[1] + alpha[0]);

  vr[2] = 11.0 * s[2] - 7.0 * s[3] + 2.0 * s[4];
  vr[1] = -s[3] + 5.0 * s[2] + 2.0 * s[1];
  vr[0] = 2.0 * s[2] + 5.0 * s[1] - s[0];

  sm = (1.0 / 6.0) * alpha1 *
       (alpha[2] * vr[2] + alpha[1] * vr[1] + alpha[0] * vr[0]);
}

// WENO 5 Z (Borges R, Carmona M, Costa B, Don WS. An improved weighted
//           essentially non-oscillatory scheme for hyperbolic conservation
//           laws. J. Comput. Phys. 2008; 227(6): 3191–3211.
//           doi: 10.1016/j.jcp.2007.11.038.)
// @param s      Stencils i-2, i-1, i, i+1, i+2
// @param sm     The value of the interpolation on the left edge of the i cell
// @param sp     The value of the interpolation on the right edge of the i cell
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
weno_reconstruct_5z(const amrex::Real* s, amrex::Real& sm, amrex::Real& sp)
{
  constexpr amrex::Real eps = 1e-6;

  amrex::Real vl[3];
  amrex::Real vr[3];
  amrex::Real beta[3];
  amrex::Real alpha[3];

  beta[2] = (13.0 / 12.0) * pow(s[0] - 2.0 * s[1] + s[2], 2) +
            0.25 * pow(s[0] - 4.0 * s[1] + 3.0 * s[2], 2);
  beta[1] = (13.0 / 12.0) * pow(s[1] - 2.0 * s[2] + s[3], 2) +
            0.25 * (s[1] - s[3]) * (s[1] - s[3]);
  beta[0] = (13.0 / 12.0) * pow(s[2] - 2.0 * s[3] + s[4], 2) +
            0.25 * pow(3.0 * s[2] - 4.0 * s[3] + s[4], 2);

  amrex::Real tau = amrex::Math::abs(beta[2] - beta[0]);

  beta[2] = 1.0 + (tau / (eps + beta[2])) * (tau / (eps + beta[2]));
  beta[1] = 1.0 + (tau / (eps + beta[1])) * (tau / (eps + beta[1]));
  beta[0] = 1.0 + (tau / (eps + beta[0])) * (tau / (eps + beta[0]));

  alpha[2] = beta[2];
  alpha[1] = 6.0 * beta[1];
  alpha[0] = 3.0 * beta[0];
  amrex::Real alpha1 = 1.0 / (alpha[2] + alpha[1] + alpha[0]);

  vl[2] = 2.0 * s[0] - 7.0 * s[1] + 11.0 * s[2];
  vl[1] = -s[1] + 5.0 * s[2] + 2.0 * s[3];
  vl[0] = 2.0 * s[2] + 5.0 * s[3] - s[4];

  sp = (1.0 / 6.0) * alpha1 *
       (alpha[2] * vl[2] + alpha[1] * vl[1] + alpha[0] * vl[0]);

  beta[2] = (13.0 / 12.0) * pow(s[4] - 2.0 * s[3] + s[2], 2) +
            0.25 * pow(s[4] - 4.0 * s[3] + 3.0 * s[2], 2);
  beta[1] = (13.0 / 12.0) * pow(s[3] - 2.0 * s[2] + s[1], 2) +
            0.25 * (s[3] - s[1]) * (s[3] - s[1]);
  beta[0] = (13.0 / 12.0) * pow(s[2] - 2.0 * s[1] + s[0], 2) +
            0.25 * pow(3.0 * s[2] - 4.0 * s[1] + s[0], 2);

  tau = amrex::Math::abs(beta[2] - beta[0]);

  beta[2] = 1.0 + (tau / (eps + beta[2])) * (tau / (eps + beta[2]));
  beta[1] = 1.0 + (tau / (eps + beta[1])) * (tau / (eps + beta[1]));
  beta[0] = 1.0 + (tau / (eps + beta[0])) * (tau / (eps + beta[0]));

  alpha[2] = beta[2];
  alpha[1] = 6.0 * beta[1];
  alpha[0] = 3.0 * beta[0];
  alpha1 = 1.0 / (alpha[2] + alpha[1] + alpha[0]);

  vr[2] = 11.0 * s[2] - 7.0 * s[3] + 2.0 * s[4];
  vr[1] = -s[3] + 5.0 * s[2] + 2.0 * s[1];
  vr[0] = 2.0 * s[2] + 5.0 * s[1] - s[0];

  sm = (1.0 / 6.0) * alpha1 *
       (alpha[2] * vr[2] + alpha[1] * vr[1] + alpha[0] * vr[0]);
}

// WENO 7 Z (Balsara DS, ShuCW. Monotonicity Preserving Weighted Essentially
//           Non-oscillatory Schemes with Increasingly High Order of Accuracy.
//           J. Comput. Phys. 2000; 160(2): 405 - 452.
//           doi: https://doi.org/10.1006/jcph.2000.6443)
// @param s      Stencils i-3, i-2, i-1, i, i+1, i+2, i+3
// @param sm     The value of the interpolation on the left edge of the i cell
// @param sp     The value of the interpolation on the right edge of the i cell
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
weno_reconstruct_7z(const amrex::Real* s, amrex::Real& sm, amrex::Real& sp)
{
  constexpr amrex::Real eps = 1e-6;

  amrex::Real vl[4];
  amrex::Real vr[4];
  amrex::Real beta[4];
  amrex::Real alpha[4];
  const amrex::Real weno7_face_wghts_1[4] = {
    1.0 / 35.0, 12.0 / 35.0, 18.0 / 35.0, 4.0 / 35.0};

  beta[3] =
    s[0] * (547.0 * s[0] - 3882.0 * s[1] + 4642.0 * s[2] - 1854.0 * s[3]) +
    s[1] * (7043.0 * s[1] - 17246.0 * s[2] + 7042.0 * s[3]) +
    s[2] * (11003.0 * s[2] - 9402.0 * s[3]) + 2107.0 * (s[3] * s[3]);

  beta[2] =
    s[1] * (267.0 * s[1] - 1642.0 * s[2] + 1602.0 * s[3] - 494.0 * s[4]) +
    s[2] * (2843.0 * s[2] - 5966.0 * s[3] + 1922.0 * s[4]) +
    s[3] * (3443.0 * s[3] - 2522.0 * s[4]) + 547.0 * (s[4] * s[4]);

  beta[1] =
    s[2] * (547.0 * s[2] - 2522.0 * s[3] + 1922.0 * s[4] - 494.0 * s[5]) +
    s[3] * (3443.0 * s[3] - 5966.0 * s[4] + 1602.0 * s[5]) +
    s[4] * (2843.0 * s[4] - 1642.0 * s[5]) + 267.0 * (s[5] * s[5]);

  beta[0] =
    s[3] * (2107.0 * s[3] - 9402.0 * s[4] + 7042.0 * s[5] - 1854.0 * s[6]) +
    s[4] * (11003.0 * s[4] - 17246.0 * s[5] + 4642.0 * s[6]) +
    s[5] * (7043.0 * s[5] - 3882.0 * s[6]) + 547.0 * (s[6] * s[6]);

  amrex::Real tau = amrex::Math::abs(beta[3] - beta[0]);

  beta[3] = 1.0 + (tau / (eps + beta[3])) * (tau / (eps + beta[3]));
  beta[2] = 1.0 + (tau / (eps + beta[2])) * (tau / (eps + beta[2]));
  beta[1] = 1.0 + (tau / (eps + beta[1])) * (tau / (eps + beta[1]));
  beta[0] = 1.0 + (tau / (eps + beta[0])) * (tau / (eps + beta[0]));

  alpha[3] = weno7_face_wghts_1[0] * beta[3];
  alpha[2] = weno7_face_wghts_1[1] * beta[2];
  alpha[1] = weno7_face_wghts_1[2] * beta[1];
  alpha[0] = weno7_face_wghts_1[3] * beta[0];
  amrex::Real alpha1 = 1.0 / (alpha[3] + alpha[2] + alpha[1] + alpha[0]);

  vl[3] = -3.0 * s[0] + 13.0 * s[1] - 23.0 * s[2] + 25.0 * s[3];
  vl[2] = 1.0 * s[1] - 5.0 * s[2] + 13.0 * s[3] + 3.0 * s[4];
  vl[1] = -1.0 * s[2] + 7.0 * s[3] + 7.0 * s[4] - 1.0 * s[5];
  vl[0] = 3.0 * s[3] + 13.0 * s[4] - 5.0 * s[5] + 1.0 * s[6];

  sp =
    (1.0 / 12.0) * alpha1 *
    (alpha[3] * vl[3] + alpha[2] * vl[2] + alpha[1] * vl[1] + alpha[0] * vl[0]);

  beta[3] =
    s[6] * (547.0 * s[6] - 3882.0 * s[5] + 4642.0 * s[4] - 1854.0 * s[3]) +
    s[5] * (7043.0 * s[5] - 17246.0 * s[4] + 7042.0 * s[3]) +
    s[4] * (11003.0 * s[4] - 9402.0 * s[3]) + 2107.0 * (s[3] * s[3]);

  beta[2] =
    s[5] * (267.0 * s[5] - 1642.0 * s[4] + 1602.0 * s[3] - 494.0 * s[2]) +
    s[4] * (2843.0 * s[4] - 5966.0 * s[3] + 1922.0 * s[2]) +
    s[3] * (3443.0 * s[3] - 2522.0 * s[2]) + 547.0 * (s[2] * s[2]);

  beta[1] =
    s[4] * (547.0 * s[4] - 2522.0 * s[3] + 1922.0 * s[2] - 494.0 * s[1]) +
    s[3] * (3443.0 * s[3] - 5966.0 * s[2] + 1602.0 * s[1]) +
    s[2] * (2843.0 * s[2] - 1642.0 * s[1]) + 267.0 * (s[1] * s[1]);

  beta[0] =
    s[3] * (2107.0 * s[3] - 9402.0 * s[2] + 7042.0 * s[1] - 1854.0 * s[0]) +
    s[2] * (11003.0 * s[2] - 17246.0 * s[1] + 4642.0 * s[0]) +
    s[1] * (7043.0 * s[1] - 3882.0 * s[0]) + 547.0 * (s[0] * s[0]);

  tau = amrex::Math::abs(beta[3] - beta[0]);

  beta[3] = 1.0 + (tau / (eps + beta[3])) * (tau / (eps + beta[3]));
  beta[2] = 1.0 + (tau / (eps + beta[2])) * (tau / (eps + beta[2]));
  beta[1] = 1.0 + (tau / (eps + beta[1])) * (tau / (eps + beta[1]));
  beta[0] = 1.0 + (tau / (eps + beta[0])) * (tau / (eps + beta[0]));

  alpha[3] = weno7_face_wghts_1[0] * beta[3];
  alpha[2] = weno7_face_wghts_1[1] * beta[2];
  alpha[1] = weno7_face_wghts_1[2] * beta[1];
  alpha[0] = weno7_face_wghts_1[3] * beta[0];
  alpha1 = 1.0 / (alpha[3] + alpha[2] + alpha[1] + alpha[0]);

  vr[3] = -3.0 * s[6] + 13.0 * s[5] - 23.0 * s[4] + 25.0 * s[3];
  vr[2] = 1.0 * s[5] - 5.0 * s[4] + 13.0 * s[3] + 3.0 * s[2];
  vr[1] = -1.0 * s[4] + 7.0 * s[3] + 7.0 * s[2] - 1.0 * s[1];
  vr[0] = 3.0 * s[3] + 13.0 * s[2] - 5.0 * s[1] + 1.0 * s[0];

  sm =
    (1.0 / 12.0) * alpha1 *
    (alpha[3] * vr[3] + alpha[2] * vr[2] + alpha[1] * vr[1] + alpha[0] * vr[0]);
}

// WENO 3 Z
// @param s      Stencils i-1, i, i+1
// @param sm     The value of the interpolation on the left edge of the i cell
// @param sp     The value of the interpolation on the right edge of the i cell
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
weno_reconstruct_3z(const amrex::Real* s, amrex::Real& sm, amrex::Real& sp)
{
  constexpr amrex::Real eps = 1e-6;

  amrex::Real vl[2];
  amrex::Real vr[2];
  amrex::Real beta[2];
  amrex::Real alpha[2];

  beta[1] = (s[0] - s[1]) * (s[0] - s[1]);
  beta[0] = (s[1] - s[2]) * (s[1] - s[2]);

  amrex::Real tau = amrex::Math::abs(beta[1] - beta[0]);

  beta[1] = 1.0 + (tau / (eps + beta[1])) * (tau / (eps + beta[1]));
  beta[0] = 1.0 + (tau / (eps + beta[0])) * (tau / (eps + beta[0]));

  alpha[1] = 2.0 * beta[1];
  alpha[0] = 1.0 * beta[0];
  amrex::Real alpha1 = 1.0 / (alpha[1] + alpha[0]);

  vl[1] = -s[0] + 3.0 * s[1];
  vl[0] = s[1] + s[2];

  sp = 0.5 * alpha1 * (alpha[1] * vl[1] + alpha[0] * vl[0]);

  beta[1] = (s[2] - s[1]) * (s[2] - s[1]);
  beta[0] = (s[1] - s[0]) * (s[1] - s[0]);

  tau = amrex::Math::abs(beta[1] - beta[0]);

  beta[1] = 1.0 + (tau / (eps + beta[1])) * (tau / (eps + beta[1]));
  beta[0] = 1.0 + (tau / (eps + beta[0])) * (tau / (eps + beta[0]));

  alpha[1] = 2.0 * beta[1];
  alpha[0] = 1.0 * beta[0];
  alpha1 = 1.0 / (alpha[1] + alpha[0]);

  vr[1] = 3.0 * s[1] - s[2];
  vr[0] = s[0] + s[1];

  sm = 0.5 * alpha1 * (alpha[1] * vr[1] + alpha[0] * vr[0]);
}

#endif
