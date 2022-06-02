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
  // const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
  // const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

  const amrex::Real rho_p =
    prob_parm.alpha * std::exp(-(x / prob_parm.sigma) * (x / prob_parm.sigma));
  const amrex::Real rho = prob_parm.rho + rho_p;
  const amrex::Real p = prob_parm.p + rho_p * prob_parm.cs * prob_parm.cs;
  const amrex::Real u_p = prob_parm.cs * rho_p / prob_parm.rho;

  amrex::Real T, eint = 0.0;
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  for (int n = 0; n < NUM_SPECIES; n++) {
    massfrac[n] = prob_parm.massfrac[n];
  }
  auto eos = pele::physics::PhysicsType::eos();
  eos.RYP2E(rho, massfrac, p, eint);
  eos.RYP2T(rho, massfrac, p, T);

  state(i, j, k, URHO) = rho;
  state(i, j, k, UMX) = rho * u_p;
  state(i, j, k, UMY) = 0.0;
  state(i, j, k, UMZ) = 0.0;
  state(i, j, k, UEINT) = rho * eint;
  state(i, j, k, UEDEN) = rho * (eint + 0.5 * u_p * u_p);
  state(i, j, k, UTEMP) = T;
  for (int n = 0; n < NUM_SPECIES; n++) {
    state(i, j, k, UFS + n) = rho * prob_parm.massfrac[n];
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

struct MyProbDeriveStruct
{
  static amrex::Box the_same_box(const amrex::Box& b) { return b; }

  static void
  add(amrex::DeriveList& derive_lst, const amrex::DescriptorList& desc_lst)
  {
    // Add derives as follows and define the derive function below:
    derive_lst.add(
      "rhoerror", amrex::IndexType::TheCellType(), 1, pc_rhoerror,
      the_same_box);
    derive_lst.addComponent("rhoerror", desc_lst, State_Type, 0, NVAR);
  }

  static void pc_rhoerror(
    const amrex::Box& bx,
    amrex::FArrayBox& derfab,
    int /*dcomp*/,
    int /*ncomp*/,
    const amrex::FArrayBox& datfab,
    const amrex::Geometry& geomdata,
    amrex::Real time,
    const int* /*bcrec*/,
    int /*level*/)
  {
    auto const dat = datfab.array();
    auto arr = derfab.array();
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_lo =
      geomdata.ProbLoArray();
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx =
      geomdata.CellSizeArray();

    ProbParmDevice const* prob_parm = PeleC::d_prob_parm_device;

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
      const amrex::Real rho_exact =
        prob_parm->rho +
        prob_parm->alpha * std::exp(
                             -((x - prob_parm->cs * time) / prob_parm->sigma) *
                             ((x - prob_parm->cs * time) / prob_parm->sigma));
      arr(i, j, k, 0) = dat(i, j, k, URHO) - rho_exact;
    });
  }
};

using ProblemDerives = MyProbDeriveStruct;

#endif
