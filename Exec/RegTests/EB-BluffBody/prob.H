#ifndef _PROB_H_
#define _PROB_H_

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
#include "ProblemDerive.H"
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
  // Geometry
  // const amrex::Real* prob_lo = geomdata.ProbLo();
  // const amrex::Real* dx = geomdata.CellSize();

  // const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
  // const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
  // const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

  state(i, j, k, URHO) = prob_parm.rho;
  state(i, j, k, UMX) = prob_parm.rho * prob_parm.vx_in;
  state(i, j, k, UMY) = prob_parm.rho * prob_parm.vy_in;
  state(i, j, k, UMZ) = 0.0;
  state(i, j, k, UEINT) = prob_parm.rho * prob_parm.eint;
  state(i, j, k, UEDEN) =
    prob_parm.rho *
    (prob_parm.eint + 0.5 * (prob_parm.vx_in * prob_parm.vx_in +
                             prob_parm.vy_in * prob_parm.vy_in));
  state(i, j, k, UTEMP) = prob_parm.T;
  for (int n = 0; n < NUM_SPECIES; n++) {
    state(i, j, k, UFS + n) = prob_parm.rho * prob_parm.massfrac[n];
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real* /*x[AMREX_SPACEDIM]*/,
  const amrex::Real* /*s_int[NVAR]*/,
  amrex::Real s_ext[NVAR],
  const int /*idir*/,
  const int /*sgn*/,
  const amrex::Real /*time*/,
  amrex::GeometryData const& /*geomdata*/,
  ProbParmDevice const& prob_parm)
{
  s_ext[URHO] = prob_parm.rho;
  s_ext[UMX] = prob_parm.rho * prob_parm.vx_in;
  s_ext[UMY] = prob_parm.rho * prob_parm.vy_in;
  s_ext[UMZ] = 0.0;
  s_ext[UEINT] = prob_parm.rho * prob_parm.eint;
  s_ext[UEDEN] = prob_parm.rho *
                 (prob_parm.eint + 0.5 * (prob_parm.vx_in * prob_parm.vx_in +
                                          prob_parm.vy_in * prob_parm.vy_in));
  s_ext[UTEMP] = prob_parm.T;
  for (int n = 0; n < NUM_SPECIES; n++) {
    s_ext[UFS + n] = prob_parm.rho * prob_parm.massfrac[n];
  }
}

void pc_prob_close();

using ProblemTags = EmptyProbTagStruct;
using ProblemDerives = EmptyProbDeriveStruct;

#endif
