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
#include "Geometry.H"

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

class EBLinePistonCylinder
  : public pele::pelec::Geometry::Register<EBLinePistonCylinder>
{
public:
  static std::string identifier() { return "Line-Piston-Cylinder"; }

  void build(const amrex::Geometry& geom, const int max_coarsening_level);
};

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
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  massfrac[CH4_ID] = prob_parm.Y_CH4_chamber;
  massfrac[O2_ID] = prob_parm.Y_O2_chamber;
  massfrac[N2_ID] = 1.0 - massfrac[O2_ID] - massfrac[CH4_ID];
  amrex::Real rho, eint;
  eos.PYT2RE(prob_parm.P_mean, massfrac, prob_parm.T_mean, rho, eint);
  state(i, j, k, URHO) = rho;
  state(i, j, k, UMX) = 0.0;
  state(i, j, k, UMY) = 0.0;
  state(i, j, k, UMZ) = 0.0;
  state(i, j, k, UEINT) = rho * eint;
  state(i, j, k, UEDEN) = state(i, j, k, UEINT);
  state(i, j, k, UTEMP) = prob_parm.T_mean;
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

    const amrex::Real decayx =
      (0.5 *
       (std::tanh((x[0] - prob_parm.win_lo[0]) * prob_parm.win_slope) + 1)) *
      (0.5 *
       (std::tanh(-(x[0] - prob_parm.win_hi[0]) * prob_parm.win_slope) + 1));
    const amrex::Real decayy =
      (0.5 *
       (std::tanh((x[1] - prob_parm.win_lo[1]) * prob_parm.win_slope) + 1)) *
      (0.5 *
       (std::tanh(-(x[1] - prob_parm.win_hi[1]) * prob_parm.win_slope) + 1));
    const amrex::Real decayz =
      (0.5 *
       (std::tanh((x[2] - prob_parm.win_lo[2]) * prob_parm.win_slope) + 1)) *
      (0.5 *
       (std::tanh(-(x[2] - prob_parm.win_hi[2]) * prob_parm.win_slope) + 1));
    u[0] *= decayx * decayy * decayz;
    u[1] *= decayx * decayy * decayz;
    u[2] *= decayx * decayy * decayz;

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

  if (
    (time < prob_parm.inj_start) ||
    (prob_parm.inj_start + prob_parm.inj_dur < time)) {
    return;
  }

  // Chamber composition
  amrex::Real massfrac_ch[NUM_SPECIES] = {0.0};
  massfrac_ch[CH4_ID] = prob_parm.Y_CH4_chamber;
  massfrac_ch[O2_ID] = prob_parm.Y_O2_chamber;
  massfrac_ch[N2_ID] = 1.0 - massfrac_ch[O2_ID] - massfrac_ch[CH4_ID];

  // Pure fuel composition
  amrex::Real massfrac_fuel[NUM_SPECIES] = {0.0};
  massfrac_fuel[NC12H26_ID] = 1.0;

  // Jet composition
  amrex::Real massfrac_jet[NUM_SPECIES] = {0.0};
  for (int n = 0; n < NUM_SPECIES; n++) {
    massfrac_jet[n] =
      prob_parm.Z * massfrac_fuel[n] + (1.0 - prob_parm.Z) * massfrac_ch[n];
  }

  // Assuming this boundary is the top XY plane
  // and in-plane angle is measured from X axis
  const amrex::Real dphi_angle = 2.0 * constants::PI() / prob_parm.nholes;
  auto eos = pele::physics::PhysicsType::eos();

  for (int nh = 1; nh <= prob_parm.nholes; nh++) {
    // Current hole center
    const amrex::Real phi_angle = nh * dphi_angle;
    const amrex::Real cx = prob_parm.centx + prob_parm.r_circ * cos(phi_angle);
    const amrex::Real cy = prob_parm.centy + prob_parm.r_circ * sin(phi_angle);

    // Injecting on an ellipse with R_minor = prob_parm.r_hole and R_major =
    // prob_parm.r_hole / cos(prob_parm.cone_angle)
    // Local coordinates:
    const amrex::Real n_loc =
      (x[0] - cx) * cos(phi_angle) - (x[1] - cy) * sin(phi_angle);
    const amrex::Real t_loc =
      (x[0] - cx) * sin(phi_angle) + (x[1] - cy) * cos(phi_angle);

    // Ellipse: major in normal, minor is tangential
    const amrex::Real elrSq =
      n_loc * n_loc /
        (prob_parm.r_hole / cos(prob_parm.cone_angle * M_PI / 180.0) *
         prob_parm.r_hole / cos(prob_parm.cone_angle * M_PI / 180.0)) +
      t_loc * t_loc / (prob_parm.r_hole * prob_parm.r_hole);

    if (elrSq < 1.2) {
      // Interior pressure
      amrex::Real p_int;
      amrex::Real massfrac_int[NUM_SPECIES] = {0.0};
      for (int n = 0; n < NUM_SPECIES; n++)
        massfrac_int[n] = s_int[UFS + n] / s_int[URHO];
      eos.RTY2P(s_int[URHO], s_int[UTEMP], massfrac_int, p_int);

      // Tanh profile, assuming 10% tanh width
      const amrex::Real eta = 0.5 * (1.0 - std::tanh((elrSq - 1.0) / 0.1));

      // Composition
      amrex::Real massfrac[NUM_SPECIES] = {0.0};
      for (int n = 0; n < NUM_SPECIES; n++) {
        massfrac[n] = eta * massfrac_jet[n] + (1.0 - eta) * massfrac_ch[n];
      }

      // subsonic inlet
      amrex::Real rho, eint;
      eos.PYT2RE(p_int, massfrac, prob_parm.T_jet, rho, eint);

      // velocity
      amrex::Real vx_in = 0.0, vy_in = 0.0, vz_in = 0.0;
      const amrex::Real vjet_pipeflow = prob_parm.vel_jet * eta;
      vz_in = -vjet_pipeflow *
              cos(prob_parm.cone_angle * M_PI / 180.0); // top XY plane
      const amrex::Real vr_in =
        vjet_pipeflow * sin(prob_parm.cone_angle * M_PI / 180.0);
      vx_in = vr_in * cos(phi_angle);
      vy_in = vr_in * sin(phi_angle);

      s_ext[URHO] = rho;
      s_ext[UMX] = rho * vx_in;
      s_ext[UMY] = rho * vy_in;
      s_ext[UMZ] = rho * vz_in;
      s_ext[UEINT] = rho * eint;
      s_ext[UEDEN] =
        rho * (eint + 0.5 * (vx_in * vx_in + vy_in * vy_in + vz_in * vz_in));
      s_ext[UTEMP] = prob_parm.T_jet;
      for (int n = 0; n < NUM_SPECIES; n++)
        s_ext[UFS + n] = rho * massfrac[n];
    }
  }
}

void pc_prob_close();

using ProblemTags = EmptyProbTagStruct;
using ProblemDerives = EmptyProbDeriveStruct;

#endif