#ifndef _PROB_H_
#define _PROB_H_

#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>

#include "mechanism.h"

#include "PeleC.H"
#include "IndexDefines.H"
#include "Constants.H"
#include "EOS.H"
#include "Transport.H"
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
  amrex::GeometryData const& geomdata)
{
  // Geometry
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* dx = geomdata.CellSize();

  const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
  // const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
  // const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

  const amrex::Real rho_p =
    ProbParm::alpha * std::exp(-(x / ProbParm::sigma) * (x / ProbParm::sigma));
  const amrex::Real rho = ProbParm::rho + rho_p;
  const amrex::Real p = ProbParm::p + rho_p * ProbParm::cs * ProbParm::cs;
  const amrex::Real u_p = ProbParm::cs * rho_p / ProbParm::rho;

  amrex::Real T, eint = 0.0;
  EOS::RYP2E(rho, ProbParm::massfrac.begin(), p, eint);
  EOS::RYP2T(rho, ProbParm::massfrac.begin(), p, T);

  state(i, j, k, URHO) = rho;
  state(i, j, k, UMX) = rho * u_p;
  state(i, j, k, UMY) = 0.0;
  state(i, j, k, UMZ) = 0.0;
  state(i, j, k, UEINT) = rho * eint;
  state(i, j, k, UEDEN) = rho * (eint + 0.5 * u_p * u_p);
  state(i, j, k, UTEMP) = T;
  for (int n = 0; n < NUM_SPECIES; n++)
    state(i, j, k, UFS + n) = rho * ProbParm::massfrac[n];
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

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
      const amrex::Real rho_exact =
        ProbParm::rho +
        ProbParm::alpha * std::exp(
                            -((x - ProbParm::cs * time) / ProbParm::sigma) *
                            ((x - ProbParm::cs * time) / ProbParm::sigma));
      arr(i, j, k, 0) = dat(i, j, k, URHO) - rho_exact;
    });
  }
};

using ProblemDerives = MyProbDeriveStruct;

#endif
