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
  const amrex::Real s_int[NVAR],
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

      u = s_int[UMX] / s_int[URHO];
      v = s_int[UMY] / s_int[URHO];
      w = s_int[UMZ] / s_int[URHO];

    } else if (sgn == -1) {

      // Following Blazek p 279, eq. 8.23

      // Interior state (point d)
      const amrex::Real* prob_hi = geomdata.ProbHi();
      const amrex::Real* dx = geomdata.CellSize();
      const amrex::Real xd = prob_hi[0] - 0.5 * dx[0];
      const amrex::Real rho_int = s_int[URHO];
      const amrex::Real u_int = s_int[UMX] / rho_int;
      const amrex::Real v_int = s_int[UMY] / rho_int;
      const amrex::Real w_int = s_int[UMZ] / rho_int;
      const amrex::Real T_int = s_int[UTEMP];
      amrex::Real p_int = 0.0, cs_int = 0.0;
      EOS::RTY2P(rho_int, T_int, ProbParm::massfrac.begin(), p_int);
      EOS::RTY2Cs(rho_int, T_int, ProbParm::massfrac.begin(), cs_int);

      // Boundary state (point b)
      const amrex::Real xb = prob_hi[0];
      const amrex::Real pb = ProbParm::dpdx * xb + ProbParm::p;
      const amrex::Real rhob = s_int[URHO] + (pb - p_int) / (cs_int * cs_int);
      const amrex::Real ub = u_int + (p_int - pb) / (rho_int * cs_int);
      const amrex::Real vb = v_int;
      const amrex::Real wb = w_int;

      // Ghost state (point a). Linear extrapolation from d and b
      rho = (rhob - rho_int) / (xb - xd) * (x[0] - xd) + rho_int;
      const amrex::Real p = (pb - p_int) / (xb - xd) * (x[0] - xd) + p_int;

      EOS::RYP2E(rho, ProbParm::massfrac.begin(), p, eint);
      EOS::EY2T(eint, ProbParm::massfrac.begin(), T);

      u = (ub - u_int) / (xb - xd) * (x[0] - xd) + u_int;
      v = (vb - v_int) / (xb - xd) * (x[0] - xd) + v_int;
      w = (wb - w_int) / (xb - xd) * (x[0] - xd) + w_int;
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
