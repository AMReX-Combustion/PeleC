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
  // Geometry
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* dx = geomdata.CellSize();
  const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
  const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
  amrex::Real z = 0;
#if AMREX_SPACEDIM == 3
  z = prob_lo[2] + (k + 0.5) * dx[2];
#endif

  // TG functions
  amrex::Real u[3] = {0.0};
  u[0] = prob_parm.v0 * sin(prob_parm.omega_x * x / prob_parm.L) *
         cos(prob_parm.omega_y * y / prob_parm.L) *
         cos(prob_parm.omega_z * z / prob_parm.L);
  u[1] = -prob_parm.v0 * cos(prob_parm.omega_x * x / prob_parm.L) *
         sin(prob_parm.omega_y * y / prob_parm.L) *
         cos(prob_parm.omega_z * z / prob_parm.L);
  if (prob_parm.convecting) {
    u[0] += prob_parm.v0;
    u[1] += prob_parm.v0;
  }
  const amrex::Real p =
    prob_parm.p0 + prob_parm.rho0 * prob_parm.v0 * prob_parm.v0 / 16.0 *
                     (cos(2.0 * prob_parm.omega_x * x / prob_parm.L) +
                      cos(2.0 * prob_parm.omega_y * y / prob_parm.L)) *
                     (cos(2.0 * prob_parm.omega_z * z / prob_parm.L) + 2.0);
  amrex::Real rho, eint;
  amrex::Real massfrac[NUM_SPECIES] = {1.0};
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

struct MyProbTagStruct
{
  AMREX_GPU_DEVICE
  AMREX_FORCE_INLINE
  static void set_problem_tags(
    const int /*i*/,
    const int /*j*/,
    const int /*k*/,
    amrex::Array4<char> const& /*tag*/,
    amrex::Array4<amrex::Real const> const& /*field*/,
    char /*tagval*/,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> /*dx*/,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> /*prob_lo*/,
    const amrex::Real /*time*/,
    const int /*level*/,
    ProbParmDevice const& /*d_prob_parm_device*/) noexcept
  {
    // could do problem specific tagging here
  }
};

using ProblemTags = MyProbTagStruct;

struct MyProbDeriveStruct
{
  static void
  add(amrex::DeriveList& /*derive_lst*/, amrex::DescriptorList& /*desc_lst*/)
  {
    // Add derives as follows and define the derive function below:
    // derive_lst.add(
    //  "varname", amrex::IndexType::TheCellType(), 1, pc_varname,
    //  the_same_box);
    // derive_lst.addComponent("varname", desc_lst, State_Type, 0, NVAR);
  }

  static void pc_varname(
    const amrex::Box& /*bx*/,
    amrex::FArrayBox& /*derfab*/,
    int /*dcomp*/,
    int /*ncomp*/,
    const amrex::FArrayBox& /*datfab*/,
    const amrex::Geometry& /*geomdata*/,
    amrex::Real /*time*/,
    const int* /*bcrec*/,
    int /*level*/)
  {
  }
};

void pc_prob_close();

using ProblemDerives = MyProbDeriveStruct;

#endif
