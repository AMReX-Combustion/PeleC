#ifndef _WENO_H_
#define _WENO_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_IArrayBox.H>

constexpr amrex::Real b1=13.0/12.0;
constexpr amrex::Real oneSixth=1.0/6.0;
constexpr amrex::Real oneHalf=1.0/2.0;
constexpr amrex::Real oneTwelve=1.0/12.0;


// Compute the coefficients of a parabolic reconstruction of the data in a
// zone. This uses the standard PPM limiters described in Colella & Woodward
// (1984)
//
// @param s      Real[5] the state to be reconstructed in zones i-2, i-1, i,
// i+1, i+2
// @param flatn  The flattening coefficient
// @param sm     The value of the parabola on the left edge of the zone
// @param sp     The value of the parabola on the right edge of the zone
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
weno_reconstruct(
  const amrex::Real* s,
  amrex::Real& sm,
  amrex::Real& sp)
{

  int wenop = 2;
  amrex::Real eps = 1.0e-40;

  amrex::Real vl[2];
  amrex::Real vr[2];
  amrex::Real beta[2];
  amrex::Real alpha[2];
  amrex::Real alpha1;

  // all index shieft by 2  -2 -> 0
  beta[2] = b1*pow(s[0]-2.0*s[1]+s[2],2) + 0.20*pow(s[0]-4.0*s[1]+3.0*s[2],2);
  beta[1] = b1*pow(s[1]-2.0*s[2]+s[3],2) + 0.25*pow(s[1]-s[3],2);
  beta[0] = b1*pow(s[2]-2.0*s[3]+s[4],2) + 0.25*pow(3.0*s[2]-4.0*s[3]+s[4],2);
//   amrex::Print() << "\n DEBUG DEBUG beta[0] = " << beta[0] << "\n"; 
//   amrex::Print() << "\n DEBUG DEBUG beta[1] = " << beta[1] << "\n"; 
//   amrex::Print() << "\n DEBUG DEBUG beta[2] = " << beta[2] << "\n";

  beta[2] = 1.0/pow(eps+beta[2],wenop);
  beta[1] = 1.0/pow(eps+beta[1],wenop);
  beta[0] = 1.0/pow(eps+beta[0],wenop);

  alpha[2] =     beta[2];
  alpha[1] = 6.0*beta[1];
  alpha[0] = 3.0*beta[0];
  alpha1 = 1.0/(alpha[2] + alpha[1] + alpha[0]);

  vl[2] = 2.0*s[0] - 7.0*s[1] + 11.0*s[2];
  vl[1] =    -s[1] + 5.0*s[2] +  2.0*s[3];
  vl[0] = 2.0*s[2] + 5.0*s[3] -      s[4];

  sp = oneSixth*alpha1*(alpha[2]*vl[2] + alpha[1]*vl[1] + alpha[0]*vl[0]);

  //-----------------------------------------------------------

  beta[2] = b1*pow(s[6]-2.0*s[5]+s[4],2) + 0.25*pow(s[6]-4.0*s[5]+3.0*s[4],2);
  beta[1] = b1*pow(s[5]-2.0*s[4]+s[3],2) + 0.25*pow(s[5]-s[3],2);
  beta[0] = b1*pow(s[4]-2.0*s[3]+s[2],2) + 0.25*pow(3.0*s[4]-4.0*s[3]+s[2],2);

  beta[2] = 1.0/pow(eps+beta[2],wenop);
  beta[1] = 1.0/pow(eps+beta[1],wenop);
  beta[0] = 1.0/pow(eps+beta[0],wenop);

  alpha[2] =     beta[2];
  alpha[1] = 6.0*beta[1];
  alpha[0] = 3.0*beta[0];
  alpha1 = 1.0/(alpha[2] + alpha[1] + alpha[0]);

  vr[2] = 11.0*s[4] - 7.0*s[5] +  2.0*s[6];
  vr[1] = vl[0];
  vr[0] = vl[1];

  sm = oneSixth*alpha1*(alpha[2]*vr[2] + alpha[1]*vr[1] + alpha[0]*vr[0]);



/* 


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
*/

}


#endif
