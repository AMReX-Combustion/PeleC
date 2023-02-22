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
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();
  const amrex::Real* dx = geomdata.CellSize();

  AMREX_D_TERM(const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
               , const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
               , const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];);

  AMREX_D_TERM(const amrex::Real Lx = prob_hi[0] - prob_lo[0];
               , const amrex::Real Ly = prob_hi[1] - prob_lo[1];
               , const amrex::Real Lz = prob_hi[2] - prob_lo[2]);

#if AMREX_SPACEDIM == 2
  const amrex::Real z = 0.0;
  const amrex::Real Lz = 1.0;
#endif

  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  massfrac[O2_ID] = prob_parm.massfrac_O2;
  massfrac[N2_ID] = prob_parm.massfrac_N2;

  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real rho, eint;
  eos.PYT2RE(prob_parm.P_mean, massfrac, prob_parm.T_mean, rho, eint);

  const amrex::Real c_f = 0.0062418;
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> vel = {0.0};
  vel[0] = std::sqrt(prob_parm.dpdx * Ly / (1e-3 * 0.4688 * c_f)) *
           (1.0 - std::pow(std::abs(2.0 * y / Ly), 10.0));

  AMREX_D_TERM(const amrex::Real u =
                 vel[0] *
                 (1.0 +
                  std::sin(2 * constants::PI() * 5 * y / Ly) *
                    std::sin(2 * constants::PI() * 14 * z / Lz) +
                  std::sin(2 * constants::PI() * 7 * y / Ly) *
                    std::sin(2 * constants::PI() * 6 * (z - 0.541) / Lz));
               , const amrex::Real v =
                   0.2 * vel[0] *
                   (std::sin(2 * constants::PI() * 5 * y / Ly) *
                      std::sin(2 * constants::PI() * 7 * (x - 0.13) / Lx) +
                    std::sin(2 * constants::PI() * 6 * y / Ly) *
                      std::sin(2 * constants::PI() * 13 * x / Lx));
               , const amrex::Real w =
                   0.2 * vel[0] *
                   (std::sin(2 * constants::PI() * 11 * z / Lz) *
                      std::sin(2 * constants::PI() * 4 * x / Lx) +
                    std::sin(2 * constants::PI() * 3 * z / Lz) *
                      std::sin(2 * constants::PI() * 9 * (x - 0.779) / Lx)));

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
