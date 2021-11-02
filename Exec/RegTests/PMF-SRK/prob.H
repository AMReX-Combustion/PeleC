#ifndef _PROB_H_
#define _PROB_H_

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
#include "PMF.H"
#include "PMFData.H"

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
  if (prob_parm.phi_in < 0.0) {
    const amrex::Real* prob_lo = geomdata.ProbLo();
    // const amrex::Real* prob_hi = geomdata.ProbHi();
    const amrex::Real* dx = geomdata.CellSize();
    const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
    const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
    const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];
    amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {{0.0}};
    amrex::Real pert = 0.0;
    if (prob_parm.pertmag > 0.0) {
      pert =
        prob_parm.pertmag *
        (1.0 * std::sin(2 * constants::PI() * 4 * x / prob_parm.L[0]) *
           std::sin(2 * constants::PI() * 5 * y / prob_parm.L[1]) +
         1.023 *
           std::sin(2 * constants::PI() * 2 * (x - 0.4598) / prob_parm.L[0]) *
           std::sin(2 * constants::PI() * 4 * (y - 0.53765) / prob_parm.L[1]) +
         0.945 *
           std::sin(2 * constants::PI() * 3 * (x - 0.712435) / prob_parm.L[0]) *
           std::sin(2 * constants::PI() * 3 * (y - 2.137) / prob_parm.L[1]) +
         1.017 *
           std::sin(2 * constants::PI() * 5 * (x - 0.33) / prob_parm.L[0]) *
           std::sin(2 * constants::PI() * 6 * (y - 1.8) / prob_parm.L[1]) +
         0.982 *
           std::sin(2 * constants::PI() * 5 * (x - 1.4234) / prob_parm.L[0]));
    }
    pele::physics::PMF::pmf(
      prob_parm.pmf_container, z + pert, z + pert, pmf_vals);
    amrex::Real molefrac[NUM_SPECIES] = {0.0};
    for (int n = 0; n < NUM_SPECIES; n++)
      molefrac[n] = pmf_vals[3 + n];
    const amrex::Real T = pmf_vals[0];
    const amrex::Real pres = prob_parm.pamb;
    amrex::Real u[3] = {0.0};
    u[AMREX_SPACEDIM - 1] = pmf_vals[1];
    amrex::Real massfrac[NUM_SPECIES] = {0.0};
    auto eos = pele::physics::PhysicsType::eos();
    eos.X2Y(molefrac, massfrac);
    amrex::Real rho = 0.0;
    amrex::Real energy = 0.0;
    eos.PYT2RE(pres, massfrac, T, rho, energy);
    state(i, j, k, URHO) = rho;
    state(i, j, k, UMX) = rho * u[0];
    state(i, j, k, UMY) = rho * u[1];
    state(i, j, k, UMZ) = rho * u[2];
    state(i, j, k, UEINT) = rho * energy;
    state(i, j, k, UEDEN) =
      rho * (energy + 0.5 * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]));
    state(i, j, k, UTEMP) = T;
    for (int n = 0; n < NUM_SPECIES; n++)
      state(i, j, k, UFS + n) = rho * massfrac[n];
  } else {
    for (int n = 0; n < NVAR; n++)
      state(i, j, k, n) = prob_parm.fuel_state[n];
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
  ProbParmDevice const& prob_parm)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();
  amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {{0.0}};
  amrex::Real u[3] = {0.0};
  amrex::Real molefrac[NUM_SPECIES] = {0.0};
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  amrex::Real /* pert, */ rho, energy, T, pres;

  auto eos = pele::physics::PhysicsType::eos();
  if (sgn == -1) {
    pele::physics::PMF::pmf(
      prob_parm.pmf_container, prob_hi[idir], prob_hi[idir], pmf_vals);
    for (int n = 0; n < NUM_SPECIES; n++)
      molefrac[n] = pmf_vals[3 + n];
    T = pmf_vals[0];
    pres = prob_parm.pamb;
    u[0] = 0.0;
    u[1] = 0.0;
    u[2] = 0.0;
    u[AMREX_SPACEDIM - 1] = pmf_vals[1];
    eos.X2Y(molefrac, massfrac);
    eos.PYT2RE(pres, massfrac, T, rho, energy);
    s_ext[URHO] = rho;
    s_ext[UMX] = rho * u[0];
    s_ext[UMY] = rho * u[1];
    s_ext[UMZ] = rho * u[2];
    s_ext[UEINT] = rho * energy;
    s_ext[UEDEN] =
      rho * (energy + 0.5 * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]));
    s_ext[UTEMP] = T;
    for (int n = 0; n < NUM_SPECIES; n++)
      s_ext[UFS + n] = rho * massfrac[n];

  } else {
    pele::physics::PMF::pmf(
      prob_parm.pmf_container, prob_lo[idir], prob_lo[idir], pmf_vals);
    for (int n = 0; n < NUM_SPECIES; n++)
      molefrac[n] = pmf_vals[3 + n];
    T = pmf_vals[0];
    pres = prob_parm.pamb;
    u[0] = 0.0;
    u[1] = 0.0;
    u[2] = 0.0;
    u[AMREX_SPACEDIM - 1] = pmf_vals[1];
    eos.X2Y(molefrac, massfrac);
    eos.PYT2RE(pres, massfrac, T, rho, energy);
    s_ext[URHO] = rho;
    s_ext[UMX] = rho * u[0];
    s_ext[UMY] = rho * u[1];
    s_ext[UMZ] = rho * u[2];
    s_ext[UEINT] = rho * energy;
    s_ext[UEDEN] =
      rho * (energy + 0.5 * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]));
    s_ext[UTEMP] = T;
    for (int n = 0; n < NUM_SPECIES; n++)
      s_ext[UFS + n] = rho * massfrac[n];
  }
}

void pc_prob_close();

using ProblemTags = EmptyProbTagStruct;
using ProblemDerives = EmptyProbDeriveStruct;

#endif
