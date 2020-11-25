#ifndef _PROB_H_
#define _PROB_H_

#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>

#include "mechanism.h"

#include "PeleC.H"
#include "IndexDefines.H"
#include "Constants.H"
#include "EOS.H"
#include "Tagging.H"
#include "Transport.H"
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
  amrex::GeometryData const& geomdata)
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
  for (int n = 0; n < NUM_SPECIES; n++)
    state(i, j, k, UFS + n) = 0;

  // Set the states
  if (ProbParm::idir == 1) {
    if (x <= ProbParm::split[0]) {
      state(i, j, k, URHO) = ProbParm::rho_l;
      state(i, j, k, UMX) = ProbParm::rho_l * ProbParm::u_l;
      state(i, j, k, UMY) = 0.0;
      state(i, j, k, UMZ) = 0.0;
      state(i, j, k, UEDEN) = ProbParm::rhoe_l + 0.5 * ProbParm::rho_l *
                                                   ProbParm::u_l *
                                                   ProbParm::u_l;
      state(i, j, k, UEINT) = ProbParm::rhoe_l;
      state(i, j, k, UTEMP) = ProbParm::T_l;
      state(i, j, k, UFS + ProbParm::left_gas_id) = state(i, j, k, URHO);
    } else {
      state(i, j, k, URHO) = ProbParm::rho_r;
      state(i, j, k, UMX) = ProbParm::rho_r * ProbParm::u_r;
      state(i, j, k, UMY) = 0.0;
      state(i, j, k, UMZ) = 0.0;
      state(i, j, k, UEDEN) = ProbParm::rhoe_r + 0.5 * ProbParm::rho_r *
                                                   ProbParm::u_r *
                                                   ProbParm::u_r;
      state(i, j, k, UEINT) = ProbParm::rhoe_r;
      state(i, j, k, UTEMP) = ProbParm::T_r;
      state(i, j, k, UFS + ProbParm::right_gas_id) = state(i, j, k, URHO);
    }
  } else if (ProbParm::idir == 2) {
    if (y <= ProbParm::split[1]) {
      state(i, j, k, URHO) = ProbParm::rho_l;
      state(i, j, k, UMX) = 0.0;
      state(i, j, k, UMY) = ProbParm::rho_l * ProbParm::u_l;
      state(i, j, k, UMZ) = 0.0;
      state(i, j, k, UEDEN) = ProbParm::rhoe_l + 0.5 * ProbParm::rho_l *
                                                   ProbParm::u_l *
                                                   ProbParm::u_l;
      state(i, j, k, UEINT) = ProbParm::rhoe_l;
      state(i, j, k, UTEMP) = ProbParm::T_l;
      state(i, j, k, UFS + ProbParm::left_gas_id) = state(i, j, k, URHO);
    } else {
      state(i, j, k, URHO) = ProbParm::rho_r;
      state(i, j, k, UMX) = 0.0;
      state(i, j, k, UMY) = ProbParm::rho_r * ProbParm::u_r;
      state(i, j, k, UMZ) = 0.0;
      state(i, j, k, UEDEN) = ProbParm::rhoe_r + 0.5 * ProbParm::rho_r *
                                                   ProbParm::u_r *
                                                   ProbParm::u_r;
      state(i, j, k, UEINT) = ProbParm::rhoe_r;
      state(i, j, k, UTEMP) = ProbParm::T_r;
      state(i, j, k, UFS + ProbParm::right_gas_id) = state(i, j, k, URHO);
    }
  } else if (ProbParm::idir == 3) {
    if (z <= ProbParm::split[2]) {
      state(i, j, k, URHO) = ProbParm::rho_l;
      state(i, j, k, UMX) = 0.0;
      state(i, j, k, UMY) = 0.0;
      state(i, j, k, UMZ) = ProbParm::rho_l * ProbParm::u_l;
      state(i, j, k, UEDEN) = ProbParm::rhoe_l + 0.5 * ProbParm::rho_l *
                                                   ProbParm::u_l *
                                                   ProbParm::u_l;
      state(i, j, k, UEINT) = ProbParm::rhoe_l;
      state(i, j, k, UTEMP) = ProbParm::T_l;
      state(i, j, k, UFS + ProbParm::left_gas_id) = state(i, j, k, URHO);
    } else {
      state(i, j, k, URHO) = ProbParm::rho_r;
      state(i, j, k, UMX) = 0.0;
      state(i, j, k, UMY) = 0.0;
      state(i, j, k, UMZ) = ProbParm::rho_r * ProbParm::u_r;
      state(i, j, k, UEDEN) = ProbParm::rhoe_r + 0.5 * ProbParm::rho_r *
                                                   ProbParm::u_r *
                                                   ProbParm::u_r;
      state(i, j, k, UEINT) = ProbParm::rhoe_r;
      state(i, j, k, UTEMP) = ProbParm::T_r;
      state(i, j, k, UFS + ProbParm::right_gas_id) = state(i, j, k, URHO);
    }
  } else {
    amrex::Abort("invalid idir");
  }
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
  amrex::GeometryData const& /*geomdata*/)
{
  for (int n = 0; n < NVAR; n++) {
    s_ext[n] = s_int[n];
  }
}

void pc_prob_close();

using ProblemTags = EmptyProbTagStruct;
using ProblemDerives = EmptyProbDeriveStruct;

#endif
