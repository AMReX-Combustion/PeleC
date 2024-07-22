#ifndef PROB_H
#define PROB_H

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

#include "mechanism.H"

#include "PeleC.H"
#include "IndexDefines.H"
#include "Constants.H"
#include "PelePhysics.H"
#include "Tagging.H"
#include "ProblemSpecificFunctions.H"
#include "prob_parm.H"

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_initdata(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& state,
  amrex::GeometryData const& /*geomdata*/,
  ProbParmDevice const& prob_parm)
{

  const amrex::Real p = prob_parm.p;
  amrex::Real rho = 0.0, eint = 0.0;
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  for (int n = 0; n < NUM_SPECIES; n++) {
    massfrac[n] = prob_parm.massfrac[n];
  }
  auto eos = pele::physics::PhysicsType::eos();
  eos.PYT2RE(p, massfrac, prob_parm.T, rho, eint);

  state(i, j, k, URHO) = rho;
  state(i, j, k, UMX) = rho * prob_parm.uavg;
  state(i, j, k, UMY) = 0.0;
  state(i, j, k, UMZ) = 0.0;
  state(i, j, k, UEINT) = rho * eint;
  state(i, j, k, UEDEN) =
    rho * (eint + 0.5 * (prob_parm.uavg * prob_parm.uavg));
  state(i, j, k, UTEMP) = prob_parm.T;
  for (int n = 0; n < NUM_SPECIES; n++) {
    state(i, j, k, UFS + n) = rho * prob_parm.massfrac[n];
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real x[AMREX_SPACEDIM],
  const amrex::Real s_inter[NVAR],
  amrex::Real s_ext[NVAR],
  const int idir,
  const int sgn,
  const amrex::Real /*time*/,
  amrex::GeometryData const& geomdata,
  ProbParmDevice const& prob_parm,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& /*turb_fluc*/)
{
  if (idir == 0) {

    amrex::Real rho = 0.0, u = 0.0, v = 0.0, w = 0.0, eint = 0.0, T = 0.0;
    auto eos = pele::physics::PhysicsType::eos();

    if (sgn == 1) {

      T = prob_parm.T;
      amrex::Real p_inter = 0.0;
      amrex::Real massfrac_inter[NUM_SPECIES] = {0.0};
      for (int n = 0; n < NUM_SPECIES; n++) {
        massfrac_inter[n] = s_inter[UFS + n] / s_inter[URHO];
      }

      // Ghost state p_d = p_inter (dp/dx = 0)
      eos.RTY2P(s_inter[URHO], s_inter[UTEMP], massfrac_inter, p_inter);

      // Ghost state rho and eint (constant T)
      eos.PYT2RE(p_inter, massfrac_inter, T, rho, eint);

      // Ghost state for velocity
      u = prob_parm.umax;
      v = 0.0;
      w = 0.0;

    } else if (sgn == -1) {

      // Following Blazek p 279, eq. 8.23

      // Interior state (point d)
      const amrex::Real* prob_hi = geomdata.ProbHi();
      const amrex::Real* dx = geomdata.CellSize();
      const amrex::Real xd = prob_hi[0] - 0.5 * dx[0];
      const amrex::Real rho_inter = s_inter[URHO];
      const amrex::Real u_inter = s_inter[UMX] / rho_inter;
      const amrex::Real v_inter = s_inter[UMY] / rho_inter;
      const amrex::Real w_inter = s_inter[UMZ] / rho_inter;
      const amrex::Real T_inter = s_inter[UTEMP];
      amrex::Real p_inter = 0.0, cs_inter = 0.0;
      amrex::Real massfrac[NUM_SPECIES] = {0.0};
      for (int n = 0; n < NUM_SPECIES; n++) {
        massfrac[n] = s_inter[UFS + n] / s_inter[URHO];
      }
      eos.RTY2P(rho_inter, T_inter, massfrac, p_inter);
      eos.RTY2Cs(rho_inter, T_inter, massfrac, cs_inter);

      // Boundary state (point b)
      const amrex::Real xb = prob_hi[0];
      const amrex::Real pb = prob_parm.p;
      const amrex::Real rhob =
        s_inter[URHO] + (pb - p_inter) / (cs_inter * cs_inter);
      const amrex::Real ub = u_inter + (p_inter - pb) / (rho_inter * cs_inter);
      const amrex::Real vb = v_inter;
      const amrex::Real wb = w_inter;

      // Ghost state (point a). Linear extrapolation from d and b
      rho = (rhob - rho_inter) / (xb - xd) * (x[0] - xd) + rho_inter;
      const amrex::Real p = (pb - p_inter) / (xb - xd) * (x[0] - xd) + p_inter;

      eos.RYP2E(rho, massfrac, p, eint);
      eos.EY2T(eint, massfrac, T);

      u = (ub - u_inter) / (xb - xd) * (x[0] - xd) + u_inter;
      v = (vb - v_inter) / (xb - xd) * (x[0] - xd) + v_inter;
      w = (wb - w_inter) / (xb - xd) * (x[0] - xd) + w_inter;
    }

    s_ext[URHO] = rho;
    s_ext[UMX] = rho * u;
    s_ext[UMY] = rho * v;
    s_ext[UMZ] = rho * w;
    s_ext[UEINT] = rho * eint;
    s_ext[UEDEN] = rho * (eint + 0.5 * (u * u + v * v + w * w));
    s_ext[UTEMP] = T;
    for (int n = 0; n < NUM_SPECIES; n++) {
      s_ext[UFS + n] = rho * prob_parm.massfrac[n];
    }
  }
}

void pc_prob_close();

using ProblemSpecificFunctions = DefaultProblemSpecificFunctions;

#endif
