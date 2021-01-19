#ifndef _WENO_H_
#define _WENO_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_IArrayBox.H>

constexpr amrex::Real b1=13.0/12.0;
constexpr amrex::Real oneSixth=1.0/6.0;
constexpr amrex::Real oneHalf=1.0/2.0;
constexpr amrex::Real oneTwelve=1.0/12.0;


// @param s      Stencils i-2, i-1, i, i+1, i+2
// @param sm     The value of the interpolation on the left edge of the i cell
// @param sp     The value of the interpolation on the right edge of the i cell
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

  beta[2] = b1*pow(s[0]-2.0*s[1]+s[2],2) + 0.20*pow(s[0]-4.0*s[1]+3.0*s[2],2);
  beta[1] = b1*pow(s[1]-2.0*s[2]+s[3],2) + 0.25*pow(s[1]-s[3],2);
  beta[0] = b1*pow(s[2]-2.0*s[3]+s[4],2) + 0.25*pow(3.0*s[2]-4.0*s[3]+s[4],2);

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

  beta[2] = b1*pow(s[4]-2.0*s[3]+s[2],2) + 0.25*pow(s[4]-4.0*s[3]+3.0*s[2],2);
  beta[1] = b1*pow(s[3]-2.0*s[2]+s[1],2) + 0.25*pow(s[3]-s[1],2);
  beta[0] = b1*pow(s[2]-2.0*s[1]+s[0],2) + 0.25*pow(3.0*s[2]-4.0*s[1]+s[0],2);

  beta[2] = 1.0/pow(eps+beta[2],wenop);
  beta[1] = 1.0/pow(eps+beta[1],wenop);
  beta[0] = 1.0/pow(eps+beta[0],wenop);

  alpha[2] =     beta[2];
  alpha[1] = 6.0*beta[1];
  alpha[0] = 3.0*beta[0];
  alpha1 = 1.0/(alpha[2] + alpha[1] + alpha[0]);

  vr[2] = 11.0*s[2] - 7.0*s[3] +  2.0*s[4];
  vr[1] =     -s[3] + 5.0*s[2] +  2.0*s[1];
  vr[0] =  2.0*s[2] + 5.0*s[1] -      s[0];

  sm = oneSixth*alpha1*(alpha[2]*vr[2] + alpha[1]*vr[1] + alpha[0]*vr[0]);


}


#endif
