#ifndef _PROB_H_
#define _PROB_H_

// AMReX
#include "AMReX.H"               // for Abort
#include "AMReX_Array.H"         // for GpuArray
#include "AMReX_Array4.H"        // for Array4
#include "AMReX_Config.H"        // for AMREX_SPACEDIM
#include "AMReX_Extension.H"     // for AMREX_FORCE_INLINE
#include "AMReX_Geometry.H"      // for GeometryData
#include "AMReX_GpuQualifiers.H" // for AMREX_GPU_DEVICE
#include "AMReX_REAL.H"          // for Real

// PelePhysics
#include "mechanism.H" // for NUM_SPECIES

// PeleC
#include "IndexDefines.H"  // for URHO, UEDEN, UEINT, UMX, UMY, UMZ
#include "ProblemDerive.H" // for EmptyProbDeriveStruct
#include "prob_parm.H"     // for ProbParmDevice
#include "Tagging.H"       // for EmptyProbTagStruct

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
  // const amrex::Real* prob_hi = geomdata.ProbHi();
  const amrex::Real* dx = geomdata.CellSize();
  const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
  const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
  amrex::Real z = 0;
#if AMREX_SPACEDIM == 3
  z = prob_lo[2] + (k + 0.5) * dx[2];
#endif

  // Set the states
  if (prob_parm.idir == 1) {
    if (x <= prob_parm.split[0]) {
      state(i, j, k, URHO) = prob_parm.rho_l;
      state(i, j, k, UMX) = prob_parm.rho_l * prob_parm.u_l;
      state(i, j, k, UMY) = 0.0;
      state(i, j, k, UMZ) = 0.0;
      state(i, j, k, UEDEN) = prob_parm.rhoe_l + 0.5 * prob_parm.rho_l *
                                                   prob_parm.u_l *
                                                   prob_parm.u_l;
      state(i, j, k, UEINT) = prob_parm.rhoe_l;
      state(i, j, k, UTEMP) = prob_parm.T_l;
    } else {
      state(i, j, k, URHO) = prob_parm.rho_r;
      state(i, j, k, UMX) = prob_parm.rho_r * prob_parm.u_r;
      state(i, j, k, UMY) = 0.0;
      state(i, j, k, UMZ) = 0.0;
      state(i, j, k, UEDEN) = prob_parm.rhoe_r + 0.5 * prob_parm.rho_r *
                                                   prob_parm.u_r *
                                                   prob_parm.u_r;
      state(i, j, k, UEINT) = prob_parm.rhoe_r;
      state(i, j, k, UTEMP) = prob_parm.T_r;
    }
  } else if (prob_parm.idir == 2) {
    if (y <= prob_parm.split[1]) {
      state(i, j, k, URHO) = prob_parm.rho_l;
      state(i, j, k, UMX) = 0.0;
      state(i, j, k, UMY) = prob_parm.rho_l * prob_parm.u_l;
      state(i, j, k, UMZ) = 0.0;
      state(i, j, k, UEDEN) = prob_parm.rhoe_l + 0.5 * prob_parm.rho_l *
                                                   prob_parm.u_l *
                                                   prob_parm.u_l;
      state(i, j, k, UEINT) = prob_parm.rhoe_l;
      state(i, j, k, UTEMP) = prob_parm.T_l;
    } else {
      state(i, j, k, URHO) = prob_parm.rho_r;
      state(i, j, k, UMX) = 0.0;
      state(i, j, k, UMY) = prob_parm.rho_r * prob_parm.u_r;
      state(i, j, k, UMZ) = 0.0;
      state(i, j, k, UEDEN) = prob_parm.rhoe_r + 0.5 * prob_parm.rho_r *
                                                   prob_parm.u_r *
                                                   prob_parm.u_r;
      state(i, j, k, UEINT) = prob_parm.rhoe_r;
      state(i, j, k, UTEMP) = prob_parm.T_r;
    }
  } else if (prob_parm.idir == 3) {
    if (z <= prob_parm.split[2]) {
      state(i, j, k, URHO) = prob_parm.rho_l;
      state(i, j, k, UMX) = 0.0;
      state(i, j, k, UMY) = 0.0;
      state(i, j, k, UMZ) = prob_parm.rho_l * prob_parm.u_l;
      state(i, j, k, UEDEN) = prob_parm.rhoe_l + 0.5 * prob_parm.rho_l *
                                                   prob_parm.u_l *
                                                   prob_parm.u_l;
      state(i, j, k, UEINT) = prob_parm.rhoe_l;
      state(i, j, k, UTEMP) = prob_parm.T_l;
    } else {
      state(i, j, k, URHO) = prob_parm.rho_r;
      state(i, j, k, UMX) = 0.0;
      state(i, j, k, UMY) = 0.0;
      state(i, j, k, UMZ) = prob_parm.rho_r * prob_parm.u_r;
      state(i, j, k, UEDEN) = prob_parm.rhoe_r + 0.5 * prob_parm.rho_r *
                                                   prob_parm.u_r *
                                                   prob_parm.u_r;
      state(i, j, k, UEINT) = prob_parm.rhoe_r;
      state(i, j, k, UTEMP) = prob_parm.T_r;
    }
  } else {
    amrex::Abort("invalid idir");
  }
  for (int n = 0; n < NUM_SPECIES; n++) {
    state(i, j, k, UFS + n) = 0;
  }
  state(i, j, k, UFS) = state(i, j, k, URHO);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real* /*x[AMREX_SPACEDIM]*/,
  const amrex::Real s_int[NVAR],
  amrex::Real s_ext[NVAR],
  const int /*idir*/,
  const int /*sgn*/,
  const amrex::Real /*time*/,
  amrex::GeometryData const& /*geomdata*/,
  ProbParmDevice const& /*prob_parm*/)
{
  for (int n = 0; n < NVAR; n++) {
    s_ext[n] = s_int[n];
  }
}

void pc_prob_close();

using ProblemTags = EmptyProbTagStruct;
using ProblemDerives = EmptyProbDeriveStruct;

#endif
