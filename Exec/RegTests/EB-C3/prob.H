#ifndef _PROB_H_
#define _PROB_H_

#include "AMReX_Array.H"         // for GpuArray
#include "AMReX_Array4.H"        // for Array4
#include "AMReX_Extension.H"     // for AMREX_FORCE_INLINE
#include "AMReX_GpuQualifiers.H" // for AMREX_GPU_DEVICE
#include "AMReX_REAL.H"          // for Real

#include "IndexDefines.H"  // for UEDEN, UEINT, UFS, UMX, UMY, UMZ, URHO
#include "mechanism.H"     // for NUM_SPECIES
#include "ProblemDerive.H" // for EmptyProbDeriveStruct
#include "prob_parm.H"     // for ProbParmDevice
#include "Tagging.H"       // for EmptyProbTagStruct

namespace amrex {
struct GeometryData;
} // namespace amrex

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
  // Set the state
  state(i, j, k, URHO) = prob_parm.rho_init;
  state(i, j, k, UMX) = 0.0;
  state(i, j, k, UMY) = 0.0;
  state(i, j, k, UMZ) = 0.0;
  state(i, j, k, UEINT) = prob_parm.rho_init * prob_parm.e_init;
  state(i, j, k, UEDEN) = prob_parm.rho_init * prob_parm.e_init;
  state(i, j, k, UTEMP) = prob_parm.T_init;
  for (int n = 0; n < NUM_SPECIES; n++) {
    state(i, j, k, UFS + n) = prob_parm.rho_init * prob_parm.massfrac[n];
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
