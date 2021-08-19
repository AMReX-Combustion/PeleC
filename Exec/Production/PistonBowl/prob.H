#ifndef _PROB_H_
#define _PROB_H_

#include <AMReX_EB2.H>
#include <AMReX_EB2_GeometryShop.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Lathe.H>
#include <AMReX_EB2_IF_Spline.H>
#include <AMReX_EB2_IF_Union.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include "mechanism.H"

#include "Constants.H"
#include "IndexDefines.H"
#include "PelePhysics.H"
#include "ProblemDerive.H"
#include "Tagging.H"
#include "Transport.H"
#include "Utilities.H"
#include "prob_parm.H"

AMREX_FORCE_INLINE
int
read_binary_int(std::ifstream& infile)
{
  int tmp = 0;
  infile.read(reinterpret_cast<char*>(&tmp), sizeof(tmp));
  return tmp;
}

AMREX_FORCE_INLINE
double
read_binary_double(std::ifstream& infile)
{
  double tmp = 0;
  infile.read(reinterpret_cast<char*>(&tmp), sizeof(tmp));
  return tmp;
}

void ReadPipeInflow(const std::string iname);

void EBLinePistonCylinder(
  const amrex::Geometry& geom, const int required_level, const int max_level);

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
InterpolateVelInflow(
  const amrex::Real x[3],
  const amrex::Real radius,
  amrex::Real& theta,
  amrex::Real& u,
  amrex::Real& v,
  amrex::Real& w,
  const amrex::Real timeInp,
  const ProbParmDevice& prob_parm)
{
  // Convert from cm to m
  amrex::Real radSI = radius * 0.01;
  // Get radius location for interpolation
  int idR = 0;
  locate(prob_parm.d_rM, prob_parm.nr, radSI, idR);
  int idRp1 = idR + 1;

  theta = std::fmod(theta, prob_parm.thetaMax);
  int idT = 0;
  locate(prob_parm.d_thetaM, prob_parm.nt, theta, idT);
  const int idTp1 = idT + 1;

  amrex::Real timeInflow = std::fmod(timeInp, prob_parm.timeInflowMax);
  int indxTime = 0;
  locate(prob_parm.d_timeInput, prob_parm.inflowNtime, timeInflow, indxTime);
  const int indxTimeP1 = indxTime + 1;

  // Interpolate in space
  amrex::Real fR1 = 0.0, fR0 = 0.0;
  if (radSI <= prob_parm.d_rM[0]) {
    fR1 = 0.0;
    fR0 = 1.0 - fR1;
  } else if (radSI >= prob_parm.d_rM[prob_parm.nr - 1]) {
    fR1 = 1.0;
    fR0 = 0.0;
    idRp1 = prob_parm.nr - 1;
    idR = prob_parm.nr - 2;
  } else {
    fR1 = (radSI - prob_parm.d_rM[idR]) /
          (prob_parm.d_rM[idRp1] - prob_parm.d_rM[idR]);
    fR0 = 1.0 - fR1;
  }

  amrex::Real fT1 = 0.0, fT0 = 0.0;
  if (theta <= prob_parm.d_thetaM[0]) {
    fT1 = 0.0;
    fT0 = 1.0 - fT1;
  } else {
    fT1 = (theta - prob_parm.d_thetaM[idT]) /
          (prob_parm.d_thetaM[idTp1] - prob_parm.d_thetaM[idT]);
    fT0 = 1.0 - fT1;
  }

  // Interpolate in time
  amrex::Real fS1 = 0.0, fS0 = 0.0;
  if (timeInflow <= prob_parm.d_timeInput[0]) {
    fS1 = 0.0;
    fS0 = 1.0 - fS1;
  } else {
    fS1 = (timeInflow - prob_parm.d_timeInput[indxTime]) /
          (prob_parm.d_timeInput[indxTimeP1] - prob_parm.d_timeInput[indxTime]);
    fS0 = 1.0 - fS1;
  }

  // Now Interpolate
  const amrex::Real ur_dns =
    fT0 *
      (fR0 *
         prob_parm.d_Ur[(indxTime * prob_parm.nr + idR) * prob_parm.nt + idT] +
       fR1 * prob_parm
               .d_Ur[(indxTime * prob_parm.nr + idRp1) * prob_parm.nt + idT]) +
    fT1 *
      (fR0 * prob_parm
               .d_Ur[(indxTime * prob_parm.nr + idR) * prob_parm.nt + idTp1] +
       fR1 * prob_parm
               .d_Ur[(indxTime * prob_parm.nr + idRp1) * prob_parm.nt + idTp1]);
  const amrex::Real ut_dns =
    fT0 *
      (fR0 *
         prob_parm.d_Ut[(indxTime * prob_parm.nr + idR) * prob_parm.nt + idT] +
       fR1 * prob_parm
               .d_Ut[(indxTime * prob_parm.nr + idRp1) * prob_parm.nt + idT]) +
    fT1 *
      (fR0 * prob_parm
               .d_Ut[(indxTime * prob_parm.nr + idR) * prob_parm.nt + idTp1] +
       fR1 * prob_parm
               .d_Ut[(indxTime * prob_parm.nr + idRp1) * prob_parm.nt + idTp1]);
  const amrex::Real uz_dns =
    fT0 *
      (fR0 *
         prob_parm.d_Uz[(indxTime * prob_parm.nr + idR) * prob_parm.nt + idT] +
       fR1 * prob_parm
               .d_Uz[(indxTime * prob_parm.nr + idRp1) * prob_parm.nt + idT]) +
    fT1 *
      (fR0 * prob_parm
               .d_Uz[(indxTime * prob_parm.nr + idR) * prob_parm.nt + idTp1] +
       fR1 * prob_parm
               .d_Uz[(indxTime * prob_parm.nr + idRp1) * prob_parm.nt + idTp1]);

  const amrex::Real ur_dns_1 =
    fT0 *
      (fR0 * prob_parm
               .d_Ur[(indxTimeP1 * prob_parm.nr + idR) * prob_parm.nt + idT] +
       fR1 *
         prob_parm
           .d_Ur[(indxTimeP1 * prob_parm.nr + idRp1) * prob_parm.nt + idT]) +
    fT1 *
      (fR0 * prob_parm
               .d_Ur[(indxTimeP1 * prob_parm.nr + idR) * prob_parm.nt + idTp1] +
       fR1 *
         prob_parm
           .d_Ur[(indxTimeP1 * prob_parm.nr + idRp1) * prob_parm.nt + idTp1]);
  const amrex::Real ut_dns_1 =
    fT0 *
      (fR0 * prob_parm
               .d_Ut[(indxTimeP1 * prob_parm.nr + idR) * prob_parm.nt + idT] +
       fR1 *
         prob_parm
           .d_Ut[(indxTimeP1 * prob_parm.nr + idRp1) * prob_parm.nt + idT]) +
    fT1 *
      (fR0 * prob_parm
               .d_Ut[(indxTimeP1 * prob_parm.nr + idR) * prob_parm.nt + idTp1] +
       fR1 *
         prob_parm
           .d_Ut[(indxTimeP1 * prob_parm.nr + idRp1) * prob_parm.nt + idTp1]);
  const amrex::Real uz_dns_1 =
    fT0 *
      (fR0 * prob_parm
               .d_Uz[(indxTimeP1 * prob_parm.nr + idR) * prob_parm.nt + idT] +
       fR1 *
         prob_parm
           .d_Uz[(indxTimeP1 * prob_parm.nr + idRp1) * prob_parm.nt + idT]) +
    fT1 *
      (fR0 * prob_parm
               .d_Uz[(indxTimeP1 * prob_parm.nr + idR) * prob_parm.nt + idTp1] +
       fR1 *
         prob_parm
           .d_Uz[(indxTimeP1 * prob_parm.nr + idRp1) * prob_parm.nt + idTp1]);

  const amrex::Real uzInterp = fS0 * uz_dns + fS1 * uz_dns_1;
  const amrex::Real urInterp = fS0 * ur_dns + fS1 * ur_dns_1;
  const amrex::Real utInterp = fS0 * ut_dns + fS1 * ut_dns_1;

  // Transform to cartesian coordinates & scale back to cm/s
  u = (urInterp * cos(theta) - radSI * utInterp * sin(theta)) * 100.0;
  v = (urInterp * sin(theta) + radSI * utInterp * cos(theta)) * 100.0;
  w = uzInterp * 100.0;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_initdata(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& state,
  amrex::GeometryData const& geomdata,
  ProbParmDevice const& prob_parm)
{
  // Geometry
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();
  const amrex::Real* dx = geomdata.CellSize();
  amrex::Real x[3] = {
    prob_lo[0] + (i + 0.5) * dx[0], prob_lo[1] + (j + 0.5) * dx[1],
    prob_lo[2] + (k + 0.5) * dx[2]};

  // Set the default state
  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real rho, eint;
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  massfrac[CH4_ID] = prob_parm.Yfuel_domain;
  massfrac[O2_ID] = prob_parm.Yox_domain;
  massfrac[N2_ID] = prob_parm.YN2_domain;
  eos.PYT2RE(prob_parm.Pres_domain, massfrac, prob_parm.Temp_domain, rho, eint);
  state(i, j, k, URHO) = rho;
  state(i, j, k, UMX) = 0.0;
  state(i, j, k, UMY) = 0.0;
  state(i, j, k, UMZ) = 0.0;
  state(i, j, k, UEINT) = rho * eint;
  state(i, j, k, UEDEN) = state(i, j, k, UEINT);
  state(i, j, k, UTEMP) = prob_parm.Temp_domain;
  for (int n = 0; n < NUM_SPECIES; n++)
    state(i, j, k, UFS + n) = rho * massfrac[n];

  if (prob_parm.hitIC) {
    amrex::Real u[3] = {0.0};
    amrex::Real mod[3] = {0.0};
    int idx[3] = {0};
    int idxp1[3] = {0};
    amrex::Real slp[3] = {0.0};
    for (int cnt = 0; cnt < 3; cnt++) {
      mod[cnt] = std::fmod(x[cnt], prob_parm.Linput);
      locate(prob_parm.d_xarray, prob_parm.inres, mod[cnt], idx[cnt]);
      idxp1[cnt] = (idx[cnt] + 1) % prob_parm.inres;
      slp[cnt] =
        (mod[cnt] - prob_parm.d_xarray[idx[cnt]]) / prob_parm.d_xdiff[idx[cnt]];
    }

    const amrex::Real f0 = (1 - slp[0]) * (1 - slp[1]) * (1 - slp[2]);
    const amrex::Real f1 = slp[0] * (1 - slp[1]) * (1 - slp[2]);
    const amrex::Real f2 = (1 - slp[0]) * slp[1] * (1 - slp[2]);
    const amrex::Real f3 = (1 - slp[0]) * (1 - slp[1]) * slp[2];
    const amrex::Real f4 = slp[0] * (1 - slp[1]) * slp[2];
    const amrex::Real f5 = (1 - slp[0]) * slp[1] * slp[2];
    const amrex::Real f6 = slp[0] * slp[1] * (1 - slp[2]);
    const amrex::Real f7 = slp[0] * slp[1] * slp[2];
    u[0] =
      prob_parm.d_uinput
          [idx[0] + prob_parm.inres * (idx[1] + prob_parm.inres * idx[2])] *
        f0 +
      prob_parm.d_uinput
          [idxp1[0] + prob_parm.inres * (idx[1] + prob_parm.inres * idx[2])] *
        f1 +
      prob_parm.d_uinput
          [idx[0] + prob_parm.inres * (idxp1[1] + prob_parm.inres * idx[2])] *
        f2 +
      prob_parm.d_uinput
          [idx[0] + prob_parm.inres * (idx[1] + prob_parm.inres * idxp1[2])] *
        f3 +
      prob_parm.d_uinput
          [idxp1[0] + prob_parm.inres * (idx[1] + prob_parm.inres * idxp1[2])] *
        f4 +
      prob_parm.d_uinput
          [idx[0] + prob_parm.inres * (idxp1[1] + prob_parm.inres * idxp1[2])] *
        f5 +
      prob_parm.d_uinput
          [idxp1[0] + prob_parm.inres * (idxp1[1] + prob_parm.inres * idx[2])] *
        f6 +
      prob_parm.d_uinput
          [idxp1[0] +
           prob_parm.inres * (idxp1[1] + prob_parm.inres * idxp1[2])] *
        f7;
    u[1] =
      prob_parm.d_vinput
          [idx[0] + prob_parm.inres * (idx[1] + prob_parm.inres * idx[2])] *
        f0 +
      prob_parm.d_vinput
          [idxp1[0] + prob_parm.inres * (idx[1] + prob_parm.inres * idx[2])] *
        f1 +
      prob_parm.d_vinput
          [idx[0] + prob_parm.inres * (idxp1[1] + prob_parm.inres * idx[2])] *
        f2 +
      prob_parm.d_vinput
          [idx[0] + prob_parm.inres * (idx[1] + prob_parm.inres * idxp1[2])] *
        f3 +
      prob_parm.d_vinput
          [idxp1[0] + prob_parm.inres * (idx[1] + prob_parm.inres * idxp1[2])] *
        f4 +
      prob_parm.d_vinput
          [idx[0] + prob_parm.inres * (idxp1[1] + prob_parm.inres * idxp1[2])] *
        f5 +
      prob_parm.d_vinput
          [idxp1[0] + prob_parm.inres * (idxp1[1] + prob_parm.inres * idx[2])] *
        f6 +
      prob_parm.d_vinput
          [idxp1[0] +
           prob_parm.inres * (idxp1[1] + prob_parm.inres * idxp1[2])] *
        f7;
    u[2] =
      prob_parm.d_winput
          [idx[0] + prob_parm.inres * (idx[1] + prob_parm.inres * idx[2])] *
        f0 +
      prob_parm.d_winput
          [idxp1[0] + prob_parm.inres * (idx[1] + prob_parm.inres * idx[2])] *
        f1 +
      prob_parm.d_winput
          [idx[0] + prob_parm.inres * (idxp1[1] + prob_parm.inres * idx[2])] *
        f2 +
      prob_parm.d_winput
          [idx[0] + prob_parm.inres * (idx[1] + prob_parm.inres * idxp1[2])] *
        f3 +
      prob_parm.d_winput
          [idxp1[0] + prob_parm.inres * (idx[1] + prob_parm.inres * idxp1[2])] *
        f4 +
      prob_parm.d_winput
          [idx[0] + prob_parm.inres * (idxp1[1] + prob_parm.inres * idxp1[2])] *
        f5 +
      prob_parm.d_winput
          [idxp1[0] + prob_parm.inres * (idxp1[1] + prob_parm.inres * idx[2])] *
        f6 +
      prob_parm.d_winput
          [idxp1[0] +
           prob_parm.inres * (idxp1[1] + prob_parm.inres * idxp1[2])] *
        f7;

    state(i, j, k, UMX) = state(i, j, k, URHO) * u[0];
    state(i, j, k, UMY) = state(i, j, k, URHO) * u[1];
    state(i, j, k, UMZ) = state(i, j, k, URHO) * u[2];
    state(i, j, k, UEDEN) =
      state(i, j, k, UEINT) +
      0.5 * state(i, j, k, URHO) * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real x[AMREX_SPACEDIM],
  const amrex::Real s_int[NVAR],
  amrex::Real s_ext[NVAR],
  const int idir,
  const int sgn,
  const amrex::Real time,
  amrex::GeometryData const& geomdata,
  ProbParmDevice const& prob_parm)
{
  const amrex::Real aYf = -33.33333333;
  const amrex::Real bYf = -5.0 - aYf;

  // default wall condition
  s_ext[URHO] = s_int[URHO];
  s_ext[UMX] = -s_int[UMX];
  s_ext[UMY] = -s_int[UMY];
  s_ext[UMZ] = -s_int[UMZ];
  s_ext[UEINT] = s_int[UEINT];
  s_ext[UEDEN] = s_int[UEDEN];
  s_ext[UTEMP] = s_int[UTEMP];
  for (int n = 0; n < NUM_SPECIES; n++)
    s_ext[UFS + n] = s_int[UFS + n];

  if (time > prob_parm.inj_time) {
    return;
  }

  // Assuming this boundary is the top XY plane
  // and in-plane angle is measured from X axis
  const amrex::Real dphi_angle = 2.0 * constants::PI() / prob_parm.nholes;

  auto eos = pele::physics::PhysicsType::eos();

  for (int nh = 1; nh <= prob_parm.nholes; nh++) {
    const amrex::Real phi_angle = nh * dphi_angle;
    const amrex::Real hole_cx =
      prob_parm.centx + prob_parm.r_circ * cos(phi_angle);
    const amrex::Real hole_cy =
      prob_parm.centz + prob_parm.r_circ * sin(phi_angle);
    if (
      ((x[0] - hole_cx) * (x[0] - hole_cx) +
       (x[1] - hole_cy) * (x[1] - hole_cy)) <
      (prob_parm.r_hole * prob_parm.r_hole)) {

      const amrex::Real rSq = (x[0] - hole_cx) * (x[0] - hole_cx) +
                              (x[1] - hole_cy) * (x[1] - hole_cy);
      const amrex::Real rad = sqrt(rSq);

      amrex::Real p_int;
      amrex::Real massfrac[NUM_SPECIES] = {0.0};
      for (int n = 0; n < NUM_SPECIES; n++)
        massfrac[n] = s_int[UFS + n] / s_int[URHO];
      eos.RTY2P(s_int[URHO], s_int[UTEMP], massfrac, p_int);

      // subsonic velocity inlet
      amrex::Real e_ext, T_ext;
      amrex::Real jet_massfrac[NUM_SPECIES] = {0.0};
      jet_massfrac[CH4_ID] =
        prob_parm.Yfuel_jet *
        (0.5 + 0.5 * tanh(aYf * (rad / prob_parm.r_hole) + bYf));
      jet_massfrac[O2_ID] = (1.0 - jet_massfrac[CH4_ID]) * 0.23;
      jet_massfrac[N2_ID] = (1.0 - jet_massfrac[CH4_ID]) * 0.77;
      eos.RYP2E(prob_parm.dens_jet, jet_massfrac, p_int, e_ext);
      eos.EY2T(e_ext, jet_massfrac, T_ext);

      amrex::Real vx_in = 0.0, vy_in = 0.0, vz_in = 0.0;
      if (prob_parm.turb_inflow_type == 1) {
        amrex::Real angle =
          constants::PI() - atan2(x[1] - hole_cy, -(x[0] - hole_cx));
        amrex::Real vxPipe = 0.0, vyPipe = 0.0, vzPipe = 0.0;
        InterpolateVelInflow(
          x, rad, angle, vxPipe, vyPipe, vzPipe, time, prob_parm);
        vy_in = vxPipe * cos(phi_angle) -
                vyPipe * sin(prob_parm.cone_angle * constants::PI() / 180.0) *
                  sin(phi_angle) +
                vzPipe * cos(prob_parm.cone_angle * constants::PI() / 180.0) *
                  sin(phi_angle);
        vx_in = -vxPipe * sin(phi_angle) -
                vyPipe * sin(prob_parm.cone_angle * constants::PI() / 180.0) *
                  cos(phi_angle) +
                vzPipe * cos(prob_parm.cone_angle * constants::PI() / 180.0) *
                  cos(phi_angle);
        vz_in =
          -vyPipe * cos(prob_parm.cone_angle * constants::PI() / 180.0) -
          vzPipe *
            sin(prob_parm.cone_angle * constants::PI() / 180.0); // top XY plane
      } else {
        const amrex::Real vjet_pipeflow =
          prob_parm.vel_jet *
          (1.0 - rSq / (prob_parm.r_hole * prob_parm.r_hole));

        // find velocity vector
        vz_in =
          -vjet_pipeflow *
          cos(prob_parm.cone_angle * constants::PI() / 180.0); // top XY plane
        const amrex::Real vr_in =
          vjet_pipeflow * sin(prob_parm.cone_angle * constants::PI() / 180.0);
        vx_in = vr_in * cos(phi_angle);
        vy_in = vr_in * sin(phi_angle);
      }

      s_ext[URHO] = prob_parm.dens_jet;
      s_ext[UMX] = prob_parm.dens_jet * vx_in;
      s_ext[UMY] = prob_parm.dens_jet * vy_in;
      s_ext[UMZ] = prob_parm.dens_jet * vz_in;
      s_ext[UEINT] = prob_parm.dens_jet * e_ext;
      s_ext[UEDEN] =
        s_ext[UEINT] + 0.5 * prob_parm.dens_jet *
                         (vx_in * vx_in + vy_in * vy_in + vz_in * vz_in);
      s_ext[UTEMP] = T_ext;
      for (int n = 0; n < NUM_SPECIES; n++)
        s_ext[UFS + n] = prob_parm.dens_jet * jet_massfrac[n];
    }
  }
}

void pc_prob_close();

using ProblemTags = EmptyProbTagStruct;
using ProblemDerives = EmptyProbDeriveStruct;

#endif
