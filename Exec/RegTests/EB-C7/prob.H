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
  // const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

  const amrex::Real cost = cos(M_PI / 180.0 * prob_parm.angle);
  const amrex::Real sint = sin(M_PI / 180.0 * prob_parm.angle);
  const amrex::Real xp = x * cost + y * sint;

  if (xp < 0.5 * prob_parm.L) {
    state(i, j, k, URHO) = prob_parm.rhol;
    state(i, j, k, UMX) = 0.0;
    state(i, j, k, UMY) = 0.0;
    state(i, j, k, UMZ) = 0.0;
    state(i, j, k, UEINT) = prob_parm.rhol * prob_parm.eintl;
    state(i, j, k, UEDEN) = prob_parm.rhol * prob_parm.eintl;
    state(i, j, k, UTEMP) = prob_parm.Tl;
    for (int n = 0; n < NUM_SPECIES; n++) {
      state(i, j, k, UFS + n) = prob_parm.rhol * prob_parm.massfrac[n];
    }
  } else {
    state(i, j, k, URHO) = prob_parm.rhor;
    state(i, j, k, UMX) = 0.0;
    state(i, j, k, UMY) = 0.0;
    state(i, j, k, UMZ) = 0.0;
    state(i, j, k, UEINT) = prob_parm.rhor * prob_parm.eintr;
    state(i, j, k, UEDEN) = prob_parm.rhor * prob_parm.eintr;
    state(i, j, k, UTEMP) = prob_parm.Tr;
    for (int n = 0; n < NUM_SPECIES; n++) {
      state(i, j, k, UFS + n) = prob_parm.rhor * prob_parm.massfrac[n];
    }
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
