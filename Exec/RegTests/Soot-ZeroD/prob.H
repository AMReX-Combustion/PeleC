#ifndef PROB_H
#define PROB_H

#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>

#include "mechanism.H"

#include "IndexDefines.H"
#include "Constants.H"
#include "PelePhysics.H"
#include "Tagging.H"
#include "Transport.H"
#include "ProblemSpecificFunctions.H"
#include "prob_parm.H"
#include "SootModel.H"

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
  // const amrex::Real* prob_hi = geomdata.ProbHi();
  // const amrex::Real* dx = geomdata.CellSize();
  // const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
  // const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
  // const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

  // TG functions
  amrex::Real u[3] = {0.0};
  const amrex::Real p = prob_parm.p0;
  amrex::Real rho, eint;
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  massfrac[O2_ID] = prob_parm.Y_O2;
  massfrac[N2_ID] = prob_parm.Y_N2;
  massfrac[prob_parm.fuelIndx] = prob_parm.Y_Fuel;
  amrex::Real Tval = prob_parm.T0;
  auto eos = pele::physics::PhysicsType::eos();
  eos.PYT2RE(p, massfrac, Tval, rho, eint);

  // Set the state
  state(i, j, k, URHO) = rho;
  state(i, j, k, UMX) = rho * u[0];
  state(i, j, k, UMY) = rho * u[1];
  state(i, j, k, UMZ) = rho * u[2];
  state(i, j, k, UEINT) = rho * eint;
  state(i, j, k, UEDEN) =
    rho * (eint + 0.5 * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]));
  state(i, j, k, UTEMP) = Tval;
  for (int n = 0; n < NUM_SPECIES; n++) {
    state(i, j, k, UFS + n) = rho * massfrac[n];
  }
  for (int n = 0; n < NUM_SOOT_MOMENTS + 1; n++) {
    state(i, j, k, UFSOOT + n) = prob_parm.soot_vals[n];
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real* /*x*/,
  const amrex::Real s_int[NVAR],
  amrex::Real s_ext[NVAR],
  const int /*idir*/,
  const int /*sgn*/,
  const amrex::Real /*time*/,
  amrex::GeometryData const& /*geomdata*/,
  ProbParmDevice const& /*prob_parm*/,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& /*turb_fluc*/)
{
  // This sets the values at the boundary equal to the
  // adjacent interior values
  for (int n = 0; n < NVAR; n++) {
    s_ext[n] = s_int[n];
  }
}

void pc_prob_close();

using ProblemSpecificFunctions = DefaultProblemSpecificFunctions;
#endif
