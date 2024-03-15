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
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  massfrac[0] = 1.0;

  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real rho, eint;
  eos.PYT2RE(prob_parm.P_mean, massfrac, prob_parm.T_mean, rho, eint);

  const amrex::Real u = prob_parm.u0;
  const amrex::Real v = prob_parm.v0;
#if AMREX_SPACEDIM == 3
  const amrex::Real w = prob_parm.w0;
#endif
  state(i, j, k, URHO) = rho;
  AMREX_D_TERM(state(i, j, k, UMX) = rho * u;, state(i, j, k, UMY) = rho * v;
               , state(i, j, k, UMZ) = rho * w;)
  state(i, j, k, UEDEN) =
    rho * eint + 0.5 * rho * (AMREX_D_TERM(u * u, +v * v, +w * w));
  state(i, j, k, UEINT) = rho * eint;
  state(i, j, k, UTEMP) = prob_parm.T_mean;
  for (int n = 0; n < NUM_SPECIES; n++) {
    state(i, j, k, UFS + n) = rho * massfrac[n];
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real* /*x[AMREX_SPACEDIM]*/,
  const amrex::Real* /*s_int[NVAR]*/,
  amrex::Real s_ext[NVAR],
  const int idir,
  const int /*sgn*/,
  const amrex::Real /*time*/,
  amrex::GeometryData const& /*geomdata*/,
  ProbParmDevice const& prob_parm,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& /*turb_fluc*/)
{
  if (prob_parm.wall_type[idir] == 1) {
    // SlipWall: flip tangential velocities so they are reflected even
    for (int i = 0; i < AMREX_SPACEDIM; i++) {
      amrex::Real flip = (i == idir) ? 1.0 : -1.0;
      s_ext[UMX + i] *= flip;
    }
  } else if (prob_parm.wall_type[idir] == 0) {
    // NoSlipWall: do nothing, s_ext comes pre-filled for NoSlip wall
  } else {
    amrex::Abort("prob.wall_type mus be 0 (NoSlip) or 1 (Slip/Symmetry)");
  }
}

void pc_prob_close();

using ProblemSpecificFunctions = DefaultProblemSpecificFunctions;

#endif
