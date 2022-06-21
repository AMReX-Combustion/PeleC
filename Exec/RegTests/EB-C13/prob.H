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
  const amrex::Real rad = std::sqrt(x * x + y * y);
  auto eos = pele::physics::PhysicsType::eos();
  const amrex::Real theta = std::atan2(y, x);

  amrex::Real rho = prob_parm.rho;
  if ((prob_parm.r_inner - dx[0] < rad) && (rad < prob_parm.r_outer + dx[0])) {
    rho = prob_parm.rho *
          std::pow(
            1 + 0.5 * (eos.gamma - 1) * prob_parm.mach * prob_parm.mach *
                  (1 - prob_parm.r_inner * prob_parm.r_inner / (rad * rad)),
            1.0 / (eos.gamma - 1));
  }

  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  for (int n = 0; n < NUM_SPECIES; n++) {
    massfrac[n] = prob_parm.massfrac[n];
  }

  const amrex::Real p = std::pow(rho, eos.gamma) / eos.gamma;
  amrex::Real T = 0.0, eint = 0.0;
  eos.RYP2T(rho, massfrac, p, T);
  eos.RTY2E(rho, T, massfrac, eint);

  const amrex::Real u =
    prob_parm.c_s * prob_parm.mach * prob_parm.r_inner / rad * std::sin(theta);
  const amrex::Real v =
    -prob_parm.c_s * prob_parm.mach * prob_parm.r_inner / rad * std::cos(theta);
  const amrex::Real w = 0.0;

  state(i, j, k, URHO) = rho;
  state(i, j, k, UMX) = rho * u;
  state(i, j, k, UMY) = rho * v;
  state(i, j, k, UMZ) = rho * w;
  state(i, j, k, UEINT) = rho * eint;
  state(i, j, k, UEDEN) = rho * (eint + 0.5 * (u * u + v * v + w * w));
  state(i, j, k, UTEMP) = T;
  for (int n = 0; n < NUM_SPECIES; n++)
    state(i, j, k, UFS + n) = rho * prob_parm.massfrac[n];
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real x[AMREX_SPACEDIM],
  const amrex::Real s_inter[NVAR],
  amrex::Real s_ext[NVAR],
  const int idir,
  const int sgn,
  const amrex::Real /*time*/,
  amrex::GeometryData const& geomdata,
  ProbParmDevice const& prob_parm)
{
  if (((idir == 0) && (sgn == 1)) || ((idir == 1) && (sgn == 1))) {

    auto eos = pele::physics::PhysicsType::eos();
    const amrex::Real* dx = geomdata.CellSize();
    const amrex::Real rad = std::sqrt(x[0] * x[0] + x[1] * x[1]);
    const amrex::Real theta = std::atan2(x[1], x[0]);
    amrex::Real rho = prob_parm.rho;
    if (
      (prob_parm.r_inner - dx[0] < rad) && (rad < prob_parm.r_outer + dx[0])) {
      rho = prob_parm.rho *
            std::pow(
              1 + 0.5 * (eos.gamma - 1) * prob_parm.mach * prob_parm.mach *
                    (1 - prob_parm.r_inner * prob_parm.r_inner / (rad * rad)),
              1.0 / (eos.gamma - 1));
    }

    amrex::Real massfrac[NUM_SPECIES] = {0.0};
    for (int n = 0; n < NUM_SPECIES; n++) {
      massfrac[n] = prob_parm.massfrac[n];
    }

    const amrex::Real p = std::pow(rho, eos.gamma) / eos.gamma;
    amrex::Real T = 0.0, eint = 0.0;
    eos.RYP2T(rho, massfrac, p, T);
    eos.RTY2E(rho, T, massfrac, eint);

    const amrex::Real u = prob_parm.c_s * prob_parm.mach * prob_parm.r_inner /
                          rad * std::sin(theta);
    const amrex::Real v = -prob_parm.c_s * prob_parm.mach * prob_parm.r_inner /
                          rad * std::cos(theta);
    const amrex::Real w = 0.0;

    s_ext[URHO] = rho;
    s_ext[UMX] = rho * u;
    s_ext[UMY] = rho * v;
    s_ext[UMZ] = rho * w;
    s_ext[UEINT] = rho * eint;
    s_ext[UEDEN] = rho * (eint + 0.5 * (u * u + v * v + w * w));
    s_ext[UTEMP] = T;
    for (int n = 0; n < NUM_SPECIES; n++)
      s_ext[UFS + n] = rho * prob_parm.massfrac[n];
  }
}

void pc_prob_close();

using ProblemTags = EmptyProbTagStruct;
using ProblemDerives = EmptyProbDeriveStruct;

#endif
