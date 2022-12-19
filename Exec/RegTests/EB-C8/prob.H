#ifndef PROB_H
#define PROB_H

#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
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
  amrex::GeometryData const& geomdata,
  ProbParmDevice const& prob_parm)
{
  // Geometry
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* dx = geomdata.CellSize();

  const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
  const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];

  const amrex::Real cost = cos(M_PI / 180.0 * prob_parm.angle);
  const amrex::Real sint = sin(M_PI / 180.0 * prob_parm.angle);
  const amrex::Real xp = x * cost + y * sint;

  for (int n = 0; n < NUM_SPECIES; n++) {
    state(i, j, k, UFS + n) = 0.0;
  }

  // Set the states
  if (xp < 0.5 * prob_parm.L) {
    state(i, j, k, URHO) = prob_parm.rho_l;
    state(i, j, k, UMX) = prob_parm.rho_l * prob_parm.u_l * cost;
    state(i, j, k, UMY) = prob_parm.rho_l * prob_parm.u_l * sint;
    state(i, j, k, UMZ) = 0.0;
    state(i, j, k, UEINT) = prob_parm.rhoe_l;
    state(i, j, k, UEDEN) =
      state(i, j, k, UEINT) + 0.5 *
                                (state(i, j, k, UMX) * state(i, j, k, UMX) +
                                 state(i, j, k, UMY) * state(i, j, k, UMY) +
                                 state(i, j, k, UMZ) * state(i, j, k, UMZ)) /
                                state(i, j, k, URHO);
    state(i, j, k, UTEMP) = prob_parm.T_l;
    state(i, j, k, UFS + prob_parm.left_gas_id) = state(i, j, k, URHO);
  } else {
    state(i, j, k, URHO) = prob_parm.rho_r;
    state(i, j, k, UMX) = prob_parm.rho_r * prob_parm.u_r * cost;
    state(i, j, k, UMY) = prob_parm.rho_r * prob_parm.u_r * sint;
    state(i, j, k, UMZ) = 0.0;
    state(i, j, k, UEINT) = prob_parm.rhoe_r;
    state(i, j, k, UEDEN) =
      state(i, j, k, UEINT) + 0.5 *
                                (state(i, j, k, UMX) * state(i, j, k, UMX) +
                                 state(i, j, k, UMY) * state(i, j, k, UMY) +
                                 state(i, j, k, UMZ) * state(i, j, k, UMZ)) /
                                state(i, j, k, URHO);
    state(i, j, k, UTEMP) = prob_parm.T_r;
    state(i, j, k, UFS + prob_parm.right_gas_id) = state(i, j, k, URHO);
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real* /*x[AMREX_SPACEDIM]*/,
  const amrex::Real* /*s_int[NVAR]*/,
  amrex::Real* /*s_ext[NVAR]*/,
  const int /*idir*/,
  const int /*sgn*/,
  const amrex::Real /*time*/,
  amrex::GeometryData const& /*geomdata*/,
  ProbParmDevice const& /*prob_parm*/)
{
}

void pc_prob_close();

using ProblemTags = EmptyProbTagStruct;
using ProblemDerives = EmptyProbDeriveStruct;

#endif
