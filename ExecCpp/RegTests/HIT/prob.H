#ifndef _PROB_H_
#define _PROB_H_

#include <AMReX_Print.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>

#include "mechanism.h"

#include "PeleC.H"
#include "IndexDefines.H"
#include "EOS.H"
#include "Tagging.H"
#include "ProblemDerive.H"
#include "prob_parm.H"
#include "Forcing.H"
#include "Utilities.H"
#include "Transport.H"

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
  amrex::Real x[3] = {
    prob_lo[0] + (i + 0.5) * dx[0], prob_lo[1] + (j + 0.5) * dx[1],
    prob_lo[2] + (k + 0.5) * dx[2]};

  // Uniform density, temperature, pressure (internal energy)
  state(i, j, k, URHO) = ProbParm::rho0;
  state(i, j, k, UEINT) = ProbParm::rho0 * ProbParm::eint0;
  state(i, j, k, UTEMP) = ProbParm::T0;
  amrex::Real massfrac[NUM_SPECIES] = {1.0};
  for (int n = 0; n < NUM_SPECIES; n++)
    state(i, j, k, UFS + n) = ProbParm::rho0 * massfrac[n];

  // Fill in the velocities and energy.
  amrex::Real u[3] = {0.0};
  amrex::Real uinterp[3] = {0.0};

  // Interpolation factors
  amrex::Real mod[3] = {0.0};
  int idx[3] = {0};
  int idxp1[3] = {0};
  amrex::Real slp[3] = {0.0};
  for (int cnt = 0; cnt < 3; cnt++) {
    mod[cnt] = std::fmod(x[cnt], ProbParm::Linput);
    locate(ProbParm::xarray, ProbParm::inres, mod[cnt], idx[cnt]);
    idxp1[cnt] = (idx[cnt] + 1) % ProbParm::inres;
    slp[cnt] =
      (mod[cnt] - ProbParm::xarray[idx[cnt]]) / ProbParm::xdiff[idx[cnt]];
  }

  const amrex::Real f0 = (1 - slp[0]) * (1 - slp[1]) * (1 - slp[2]);
  const amrex::Real f1 = slp[0] * (1 - slp[1]) * (1 - slp[2]);
  const amrex::Real f2 = (1 - slp[0]) * slp[1] * (1 - slp[2]);
  const amrex::Real f3 = (1 - slp[0]) * (1 - slp[1]) * slp[2];
  const amrex::Real f4 = slp[0] * (1 - slp[1]) * slp[2];
  const amrex::Real f5 = (1 - slp[0]) * slp[1] * slp[2];
  const amrex::Real f6 = slp[0] * slp[1] * (1 - slp[2]);
  const amrex::Real f7 = slp[0] * slp[1] * slp[2];

  uinterp[0] =
    ProbParm::uinput
        [idx[0] + ProbParm::inres * (idx[1] + ProbParm::inres * idx[2])] *
      f0 +
    ProbParm::uinput
        [idxp1[0] + ProbParm::inres * (idx[1] + ProbParm::inres * idx[2])] *
      f1 +
    ProbParm::uinput
        [idx[0] + ProbParm::inres * (idxp1[1] + ProbParm::inres * idx[2])] *
      f2 +
    ProbParm::uinput
        [idx[0] + ProbParm::inres * (idx[1] + ProbParm::inres * idxp1[2])] *
      f3 +
    ProbParm::uinput
        [idxp1[0] + ProbParm::inres * (idx[1] + ProbParm::inres * idxp1[2])] *
      f4 +
    ProbParm::uinput
        [idx[0] + ProbParm::inres * (idxp1[1] + ProbParm::inres * idxp1[2])] *
      f5 +
    ProbParm::uinput
        [idxp1[0] + ProbParm::inres * (idxp1[1] + ProbParm::inres * idx[2])] *
      f6 +
    ProbParm::uinput
        [idxp1[0] + ProbParm::inres * (idxp1[1] + ProbParm::inres * idxp1[2])] *
      f7;
  uinterp[1] =
    ProbParm::vinput
        [idx[0] + ProbParm::inres * (idx[1] + ProbParm::inres * idx[2])] *
      f0 +
    ProbParm::vinput
        [idxp1[0] + ProbParm::inres * (idx[1] + ProbParm::inres * idx[2])] *
      f1 +
    ProbParm::vinput
        [idx[0] + ProbParm::inres * (idxp1[1] + ProbParm::inres * idx[2])] *
      f2 +
    ProbParm::vinput
        [idx[0] + ProbParm::inres * (idx[1] + ProbParm::inres * idxp1[2])] *
      f3 +
    ProbParm::vinput
        [idxp1[0] + ProbParm::inres * (idx[1] + ProbParm::inres * idxp1[2])] *
      f4 +
    ProbParm::vinput
        [idx[0] + ProbParm::inres * (idxp1[1] + ProbParm::inres * idxp1[2])] *
      f5 +
    ProbParm::vinput
        [idxp1[0] + ProbParm::inres * (idxp1[1] + ProbParm::inres * idx[2])] *
      f6 +
    ProbParm::vinput
        [idxp1[0] + ProbParm::inres * (idxp1[1] + ProbParm::inres * idxp1[2])] *
      f7;
  uinterp[2] =
    ProbParm::winput
        [idx[0] + ProbParm::inres * (idx[1] + ProbParm::inres * idx[2])] *
      f0 +
    ProbParm::winput
        [idxp1[0] + ProbParm::inres * (idx[1] + ProbParm::inres * idx[2])] *
      f1 +
    ProbParm::winput
        [idx[0] + ProbParm::inres * (idxp1[1] + ProbParm::inres * idx[2])] *
      f2 +
    ProbParm::winput
        [idx[0] + ProbParm::inres * (idx[1] + ProbParm::inres * idxp1[2])] *
      f3 +
    ProbParm::winput
        [idxp1[0] + ProbParm::inres * (idx[1] + ProbParm::inres * idxp1[2])] *
      f4 +
    ProbParm::winput
        [idx[0] + ProbParm::inres * (idxp1[1] + ProbParm::inres * idxp1[2])] *
      f5 +
    ProbParm::winput
        [idxp1[0] + ProbParm::inres * (idxp1[1] + ProbParm::inres * idx[2])] *
      f6 +
    ProbParm::winput
        [idxp1[0] + ProbParm::inres * (idxp1[1] + ProbParm::inres * idxp1[2])] *
      f7;

  u[0] = uinterp[0] + forcing_params::u0;
  u[1] = uinterp[1] + forcing_params::v0;
  u[2] = uinterp[2] + forcing_params::w0;
  state(i, j, k, UMX) = ProbParm::rho0 * u[0];
  state(i, j, k, UMY) = ProbParm::rho0 * u[1];
  state(i, j, k, UMZ) = ProbParm::rho0 * u[2];
  state(i, j, k, UEDEN) =
    ProbParm::rho0 *
    (ProbParm::eint0 + 0.5 * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]));
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
  amrex::GeometryData const& /*geomdata*/)
{
}

void pc_prob_close();

using ProblemTags = EmptyProbTagStruct;
using ProblemDerives = EmptyProbDeriveStruct;

#endif
