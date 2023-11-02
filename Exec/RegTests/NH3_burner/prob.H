#ifndef PROB_H
#define PROB_H

#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

#include "PeleC.H"
#include "IndexDefines.H"
#include "PelePhysics.H"
#include "Tagging.H"
#include "ProblemDerive.H"
#include "prob_parm.H"
#include "Constants.H"

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
  amrex::Real x[3] = {
    prob_lo[0] + static_cast<amrex::Real>(i + 0.5) * dx[0],
    prob_lo[1] + static_cast<amrex::Real>(j + 0.5) * dx[1],
    prob_lo[2] + static_cast<amrex::Real>(k + 0.5) * dx[2]};
  
  amrex::Real massfrac[NUM_SPECIES] = {{0.0}};
  amrex::Real u[3] = {0.0};
  amrex::Real T    = prob_parm.T_in;

  massfrac[H2_ID]  = 0.00;
  massfrac[NH3_ID] = 0.1418434;
  massfrac[O2_ID]  = 0.19987257;
  // // massfrac[N2_ID]  = 0.65828399;
  massfrac[N2_ID]  = 1.0 - massfrac[NH3_ID]-massfrac[O2_ID];

  auto eos = pele::physics::PhysicsType::eos();

  const amrex::Real pres = prob_parm.pamb*10.; //[dyne/cm3]

  amrex::Real rho    = 0.0;
  amrex::Real energy = 0.0;
  eos.PYT2RE(pres, massfrac, T, rho, energy);

  u[0] = 0.0;
  u[1] = 0.0;
  u[2] = 0.0;

  state(i, j, k, URHO)  = rho;
  state(i, j, k, UEINT) = rho * energy;
  state(i, j, k, UTEMP) = T;
  state(i, j, k, UMX)   = rho * u[0];
  state(i, j, k, UMY)   = rho * u[1];
  state(i, j, k, UMZ)   = rho * u[2];
  state(i, j, k, UEDEN) =
    rho *
    (energy + 0.5 * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]));
  
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
  const int sgn,
  const amrex::Real /*time*/,
  amrex::GeometryData const& geomdata,
  ProbParmDevice const& prob_parm,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& turb_fluc)
{

  amrex::Real massfrac[NUM_SPECIES] = {{0.0}};
  amrex::Real T          = prob_parm.T_in;
  const amrex::Real pres = prob_parm.pamb*10.; //[dyne/cm3]
  amrex::Real u[3] = {0.0};

  massfrac[H2_ID]  = 0.00;
  massfrac[NH3_ID] = 0.1418434;
  massfrac[O2_ID]  = 0.19987257;
  // massfrac[N2_ID]  = 0.65828399;

  massfrac[N2_ID]  = 1.0 - massfrac[NH3_ID]-massfrac[O2_ID];

  amrex::Real rho    = 0.0;
  amrex::Real energy = 0.0;

  auto eos = pele::physics::PhysicsType::eos();
  eos.PYT2RE(pres, massfrac, T, rho, energy);

  u[0] = turb_fluc[0];
  u[1] = turb_fluc[1];
  u[2] = turb_fluc[2];

  s_ext[URHO]  = rho;
  s_ext[UMX]   = rho * u[0];
  s_ext[UMY]   = rho * u[1];
  s_ext[UMZ]   = rho * u[2];
  s_ext[UEINT] = rho * energy;
  s_ext[UEDEN] =
    rho * (energy + 0.5 * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]));
  s_ext[UTEMP] = T;
  for (int n = 0; n < NUM_SPECIES; n++) {
    s_ext[UFS + n] = rho * massfrac[n];
  }

}

void pc_prob_close();

using ProblemTags = EmptyProbTagStruct;
using ProblemDerives = EmptyProbDeriveStruct;

#endif
