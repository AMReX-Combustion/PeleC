#ifndef _PROB_H_
#define _PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>

#include "mechanism.h"

#include "PeleC.H"
#include "IndexDefines.H"
#include "Constants.H"
#include "EOS.H"
#include "Tagging.H"
#include "ProblemDerive.H"
#include "prob_parm.H"
#include "Transport.H"

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_initdata(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& state,
  amrex::GeometryData const& geomdata)
{
  // Geometry
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* dx = geomdata.CellSize();

  const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
  // const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
  // const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

  const amrex::Real p = ProbParm::dpdx * x + ProbParm::p;
  amrex::Real rho = 0.0, eint = 0.0;
  EOS::PYT2RE(p, ProbParm::massfrac.begin(), ProbParm::T, rho, eint);

  state(i, j, k, URHO) = rho;
  state(i, j, k, UMX) = rho * ProbParm::uavg;
  state(i, j, k, UMY) = 0.0;
  state(i, j, k, UMZ) = 0.0;
  state(i, j, k, UEINT) = rho * eint;
  state(i, j, k, UEDEN) =
    rho * (eint + 0.5 * (ProbParm::uavg * ProbParm::uavg));
  state(i, j, k, UTEMP) = ProbParm::T;
  for (int n = 0; n < NUM_SPECIES; n++)
    state(i, j, k, UFS + n) = rho * ProbParm::massfrac[n];
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
  amrex::GeometryData const& geomdata)
{
  if (idir == 0) {
    amrex::Real rho = 0.0, u = 0.0, v = 0.0, w = 0.0, eint = 0.0, T = 0.0;

    if (sgn == 1) {
      const amrex::Real p = ProbParm::dpdx * x[0] + ProbParm::p;
      T = ProbParm::T;
      EOS::PYT2RE(p, ProbParm::massfrac.begin(), T, rho, eint);

      u = s_inter[UMX] / s_inter[URHO];
      v = s_inter[UMY] / s_inter[URHO];
      w = s_inter[UMZ] / s_inter[URHO];

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
      EOS::RTY2P(rho_inter, T_inter, ProbParm::massfrac.begin(), p_inter);
      EOS::RTY2Cs(rho_inter, T_inter, ProbParm::massfrac.begin(), cs_inter);

      // Boundary state (point b)
      const amrex::Real xb = prob_hi[0];
      const amrex::Real pb = ProbParm::dpdx * xb + ProbParm::p;
      const amrex::Real rhob =
        s_inter[URHO] + (pb - p_inter) / (cs_inter * cs_inter);
      const amrex::Real ub = u_inter + (p_inter - pb) / (rho_inter * cs_inter);
      const amrex::Real vb = v_inter;
      const amrex::Real wb = w_inter;

      // Ghost state (point a). Linear extrapolation from d and b
      rho = (rhob - rho_inter) / (xb - xd) * (x[0] - xd) + rho_inter;
      const amrex::Real p = (pb - p_inter) / (xb - xd) * (x[0] - xd) + p_inter;

      EOS::RYP2E(rho, ProbParm::massfrac.begin(), p, eint);
      EOS::EY2T(eint, ProbParm::massfrac.begin(), T);

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
    for (int n = 0; n < NUM_SPECIES; n++)
      s_ext[UFS + n] = rho * ProbParm::massfrac[n];
  }
}

void pc_prob_close();

using ProblemTags = EmptyProbTagStruct;
using ProblemDerives = EmptyProbDeriveStruct;

#endif
