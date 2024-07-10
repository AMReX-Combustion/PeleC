#ifndef PROB_H
#define PROB_H

#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_IntVect.H>

#include "mechanism.H"

#include "IndexDefines.H"
#include "Constants.H"
#include "PelePhysics.H"
#include "Tagging.H"
#include "Transport.H"
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
  amrex::GeometryData const& geomdata,
  ProbParmDevice const& prob_parm)
{

  // Geometry
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* dx = geomdata.CellSize();
  amrex::Real L = prob_parm.L;
  amrex::Real x = 0.;
  amrex::Real y = 0.;
  amrex::Real z = 0.;
  AMREX_D_TERM(x = prob_lo[0] + (i + 0.5) * dx[0];
               , y = prob_lo[1] + (j + 0.5) * dx[1];
               , z = prob_lo[2] + (k + 0.5) * dx[2];);
  const amrex::Real coef = 2. * M_PI;

  // TG functions
  amrex::Real u[3] = {0.0};
  u[0] = prob_parm.v0 * std::sin(coef * x / L) * std::cos(coef * y / L) *
         std::cos(coef * z / L);
  u[1] = -prob_parm.v0 * std::cos(coef * x / L) * std::sin(coef * y / L) *
         std::cos(coef * z / L);
  if (prob_parm.convecting) {
    u[0] += prob_parm.v0;
    u[1] += prob_parm.v0;
  }
  const amrex::Real p =
    prob_parm.p0 +
    prob_parm.rho0 * prob_parm.v0 * prob_parm.v0 *
      (std::cos(2.0 * coef * x / L) + std::cos(2.0 * coef * y / L)) *
      (std::cos(2.0 * coef * z / L) + 2.0) / 16.;
  amrex::Real rho, eint;
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  massfrac[N2_ID] = 1.;
  auto eos = pele::physics::PhysicsType::eos();
  eos.PYT2RE(p, massfrac, prob_parm.T0, rho, eint);

  // Set the state
  state(i, j, k, URHO) = rho;
  state(i, j, k, UMX) = rho * u[0];
  state(i, j, k, UMY) = rho * u[1];
  state(i, j, k, UMZ) = rho * u[2];
  state(i, j, k, UEINT) = rho * eint;
  state(i, j, k, UEDEN) =
    rho * (eint + 0.5 * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]));
  state(i, j, k, UTEMP) = prob_parm.T0;
  for (int n = 0; n < NUM_SPECIES; n++) {
    state(i, j, k, UFS + n) = rho * massfrac[n];
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real* /*x*/,
  const amrex::Real* /*s_int[NVAR]*/,
  amrex::Real* /*s_ext[NVAR]*/,
  const int /*idir*/,
  const int /*sgn*/,
  const amrex::Real /*time*/,
  amrex::GeometryData const& /*geomdata*/,
  ProbParmDevice const& /*prob_parm*/,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& /*turb_fluc*/)
{
}

void pc_prob_close();

using ProblemSpecificFunctions = DefaultProblemSpecificFunctions;

void diameter_solver(
  const amrex::Real rho0,
  const amrex::Real mu,
  const amrex::Real tau_d,
  const amrex::Real partRho,
  const amrex::Real refRho,
  const amrex::Real refU,
  amrex::Real& dia,
  amrex::Real& Re_d);

void prob_parser(
  const amrex_real* problo,
  const amrex_real* probhi,
  amrex::Real& reynolds,
  amrex::Real& mach,
  bool& convecting,
  amrex::Real& p0,
  amrex::Real& T0,
  amrex::Real& rho0,
  amrex::Real& Stmod,
  amrex::Real& rhoRatio,
  amrex::Real& L,
  amrex::Real& v0,
  amrex::Real& mu,
  amrex::Real& cs,
  amrex::Real& cp,
  amrex::Real& St_num);

#endif
