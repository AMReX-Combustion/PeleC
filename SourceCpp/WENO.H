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
weno_reconstruct_5js(
  const amrex::Real* s,
  amrex::Real& sm,
  amrex::Real& sp)
{

  int wenop = 2;
  amrex::Real eps = 1.0e-40;

  amrex::Real vl[3];
  amrex::Real vr[3];
  amrex::Real beta[3];
  amrex::Real alpha[3];
  amrex::Real alpha1;

  beta[2] = b1*pow(s[0]-2.0*s[1]+s[2],2) + 0.25*pow(s[0]-4.0*s[1]+3.0*s[2],2);
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


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
weno_reconstruct_5z(
  const amrex::Real* s,
  amrex::Real& sm,
  amrex::Real& sp)
{

  int wenop = 2;
  amrex::Real eps = 1.0e-40;

  amrex::Real vl[3];
  amrex::Real vr[3];
  amrex::Real beta[3];
  amrex::Real alpha[3];
  amrex::Real alpha1;
  amrex::Real tau;

  beta[2] = b1*pow(s[0]-2.0*s[1]+s[2],2) + 0.25*pow(s[0]-4.0*s[1]+3.0*s[2],2);
  beta[1] = b1*pow(s[1]-2.0*s[2]+s[3],2) + 0.25*pow(s[1]-s[3],2);
  beta[0] = b1*pow(s[2]-2.0*s[3]+s[4],2) + 0.25*pow(3.0*s[2]-4.0*s[3]+s[4],2);

  tau = abs(beta[2] - beta[0]);

  beta[2] = 1.0 + pow(tau / (eps+beta[2]),wenop);
  beta[1] = 1.0 + pow(tau / (eps+beta[1]),wenop);
  beta[0] = 1.0 + pow(tau / (eps+beta[0]),wenop);

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

  tau = abs(beta[2] - beta[0]);

  beta[2] = 1.0 + pow(tau / (eps+beta[2]),wenop);
  beta[1] = 1.0 + pow(tau / (eps+beta[1]),wenop);
  beta[0] = 1.0 + pow(tau / (eps+beta[0]),wenop);

  alpha[2] =     beta[2];
  alpha[1] = 6.0*beta[1];
  alpha[0] = 3.0*beta[0];
  alpha1 = 1.0/(alpha[2] + alpha[1] + alpha[0]);

  vr[2] = 11.0*s[2] - 7.0*s[3] +  2.0*s[4];
  vr[1] =     -s[3] + 5.0*s[2] +  2.0*s[1];
  vr[0] =  2.0*s[2] + 5.0*s[1] -      s[0];

  sm = oneSixth*alpha1*(alpha[2]*vr[2] + alpha[1]*vr[1] + alpha[0]*vr[0]);

}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
weno_reconstruct_7z(
  const amrex::Real* s,
  amrex::Real& sm,
  amrex::Real& sp)
{

  int wenop = 2;
  amrex::Real eps = 1.0e-40;

  amrex::Real vl[4];
  amrex::Real vr[4];
  amrex::Real beta[4];
  amrex::Real alpha[4];
  amrex::Real alpha1;
  amrex::Real tau;
  amrex::Real weno7_face_wghts_1[4];

  weno7_face_wghts_1[0]  =  1./35.;
  weno7_face_wghts_1[1]  = 12./35.;
  weno7_face_wghts_1[2]  = 18./35.;
  weno7_face_wghts_1[3]  =  4./35.; 

  beta[3] = s[0]*(547.0*s[0] - 3882.0*s[1] + 4642.0*s[2] - 1854.0*s[3]) + 
            s[1]*(7043.0*s[1] - 17246.0*s[2] + 7042.0*s[3]) + 
            s[2]*(11003.0*s[2] - 9402.0*s[3]) + 2107.0*pow(s[3],2);

  beta[2] = s[1]*(267.0*s[1] - 1642.0*s[2] + 1602.0*s[3] - 494.0*s[4]) + 
            s[2]*(2843.0*s[2] - 5966.0*s[3] + 1922.0*s[4]) + 
            s[3]*(3443.0*s[3] - 2522.0*s[4]) + 547.0*pow(s[4],2);

  beta[1] = s[2]*(547.0*s[2] - 2522.0*s[3] + 1922.0*s[4] - 494.0*s[5]) + 
            s[3]*(3443.0*s[3] - 5966.0*s[4] + 1602.0*s[5]) + 
            s[4]*(2843.0*s[4] - 1642.0*s[5]) + 267.0*pow(s[5],2);

  beta[0] = s[3]*(2107.0*s[3] - 9402.0*s[4] + 7042.0*s[5] - 1854.0*s[6]) + 
            s[4]*(11003.0*s[4] - 17246.0*s[5] + 4642.0*s[6]) + 
            s[5]*(7043.0*s[5] - 3882.0*s[6]) + 547.0*pow(s[6],2);

  tau = abs(beta[3] - beta[0]);

  beta[3] = 1.0 + pow(tau / (eps+beta[3]),wenop);
  beta[2] = 1.0 + pow(tau / (eps+beta[2]),wenop);
  beta[1] = 1.0 + pow(tau / (eps+beta[1]),wenop);
  beta[0] = 1.0 + pow(tau / (eps+beta[0]),wenop);

  alpha[3] = weno7_face_wghts_1[0] *  beta[3];
  alpha[2] = weno7_face_wghts_1[1] *  beta[2];
  alpha[1] = weno7_face_wghts_1[2] *  beta[1];
  alpha[0] = weno7_face_wghts_1[3] *  beta[0];
  alpha1 = 1.0/(alpha[3] + alpha[2] + alpha[1] + alpha[0]);

  vl[3] = -3.0*s[0] + 13.0*s[1] - 23.0*s[2] + 25.0*s[3];
  vl[2] = 1.0*s[1] - 5.0*s[2] + 13.0*s[3] + 3.0*s[4];
  vl[1] = -1.0*s[2] + 7.0*s[3] + 7.0*s[4] - 1.0*s[5];
  vl[0] = 3.0*s[3] + 13.0*s[4] - 5.0*s[5] + 1.0*s[6];

  sp = oneTwelve*alpha1*(alpha[3]*vl[3] + alpha[2]*vl[2] 
                       + alpha[1]*vl[1] + alpha[0]*vl[0]);

  //-----------------------------------------------------------

  beta[3] = s[6]*(547.0*s[6] - 3882.0*s[5] + 4642.0*s[4] - 1854.0*s[3]) +
            s[5]*(7043.0*s[5] - 17246.0*s[4] + 7042.0*s[3]) +
            s[4]*(11003.0*s[4] - 9402.0*s[3]) + 2107.0*pow(s[3],2);

  beta[2] = s[5]*(267.0*s[5] - 1642.0*s[4] + 1602.0*s[3] - 494.0*s[2]) +
            s[4]*(2843.0*s[4] - 5966.0*s[3] + 1922.0*s[2]) +
            s[3]*(3443.0*s[3] - 2522.0*s[2]) + 547.0*pow(s[2],2);

  beta[1] = s[4]*(547.0*s[4] - 2522.0*s[3] + 1922.0*s[2] - 494.0*s[1]) +
            s[3]*(3443.0*s[3] - 5966.0*s[2] + 1602.0*s[1]) +
            s[2]*(2843.0*s[2] - 1642.0*s[1]) + 267.0*pow(s[1],2);

  beta[0] = s[3]*(2107.0*s[3] - 9402.0*s[2] + 7042.0*s[1] - 1854.0*s[0]) +
            s[2]*(11003.0*s[2] - 17246.0*s[1] + 4642.0*s[0]) +
            s[1]*(7043.0*s[1] - 3882.0*s[0]) + 547.0*pow(s[0],2);

  tau = abs(beta[3] - beta[0]);

  beta[3] = 1.0 + pow(tau / (eps+beta[3]),wenop);
  beta[2] = 1.0 + pow(tau / (eps+beta[2]),wenop);
  beta[1] = 1.0 + pow(tau / (eps+beta[1]),wenop);
  beta[0] = 1.0 + pow(tau / (eps+beta[0]),wenop);

  alpha[3] = weno7_face_wghts_1[0] *  beta[3];
  alpha[2] = weno7_face_wghts_1[1] *  beta[2];
  alpha[1] = weno7_face_wghts_1[2] *  beta[1];
  alpha[0] = weno7_face_wghts_1[3] *  beta[0];
  alpha1 = 1.0/(alpha[3] + alpha[2] + alpha[1] + alpha[0]);

  vr[3] = -3.0*s[6] + 13.0*s[5] - 23.0*s[4] + 25.0*s[3];
  vr[2] = 1.0*s[5] - 5.0*s[4] + 13.0*s[3] + 3.0*s[2];
  vr[1] = -1.0*s[4] + 7.0*s[3] + 7.0*s[2] - 1.0*s[1];
  vr[0] = 3.0*s[3] + 13.0*s[2] - 5.0*s[1] + 1.0*s[0];  

  sm = oneTwelve*alpha1*(alpha[3]*vr[3] + alpha[2]*vr[2] 
                       + alpha[1]*vr[1] + alpha[0]*vr[0]);


}

#endif
