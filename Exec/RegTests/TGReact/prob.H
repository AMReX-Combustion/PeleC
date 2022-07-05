#ifndef _PROB_H_
#define _PROB_H_

#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>

#include "mechanism.h"

#include "IndexDefines.H"
#include "Constants.H"
#include "EOS.H"
#include "Tagging.H"
#include "Transport.H"
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
  const amrex::Real* prob_hi = geomdata.ProbHi();
  const amrex::Real* dx = geomdata.CellSize();
  const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
  const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
  const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];
  const amrex::Real p_si2cgs = 10.0;

  // TG functions
  amrex::Real u[3] = {0.0};
  u[0] = ProbParm::u0 * sin(ProbParm::omega_x * x / ProbParm::L_x) *
         cos(ProbParm::omega_y * y / ProbParm::L_y) *
         cos(ProbParm::omega_z * z / ProbParm::L_z);
  u[1] = -ProbParm::u0 * cos(ProbParm::omega_x * x / ProbParm::L_x) *
         sin(ProbParm::omega_y * y / ProbParm::L_y) *
         cos(ProbParm::omega_z * z / ProbParm::L_z);
  amrex::Real rho, eint, p, T;
  amrex::Real massfrac[NUM_SPECIES] = {0.0};

  if ((ProbParm::case_type_int == 1) || (ProbParm::case_type_int == 2)) {
    amrex::Real xmod = std::fmod(x, ProbParm::L_x);
    int m = 0;
    locate(ProbParm::xarray, ProbParm::nx, xmod, m);
    const int mp1 = (m + 1) % ProbParm::nx;
    const amrex::Real fact =
      (x - ProbParm::input[ProbParm::input_x + m * ProbParm::nvars]) /
      ProbParm::dxinput[m];
    const amrex::Real f0 = (1 - fact);
    const amrex::Real f1 = fact;

    p = (ProbParm::input[ProbParm::input_p + m * ProbParm::nvars] * f0 +
         ProbParm::input[ProbParm::input_p + mp1 * ProbParm::nvars] * f1) *
        p_si2cgs;
    T = ProbParm::input[ProbParm::input_Tad + m * ProbParm::nvars] * f0 +
        ProbParm::input[ProbParm::input_Tad + mp1 * ProbParm::nvars] * f1;
#ifdef PELEC_USE_REACTIONS
    if (ProbParm::case_type_int == 1) {
      massfrac[H2_ID] =
        ProbParm::input[ProbParm::input_H2_init + m * ProbParm::nvars] * f0 +
        ProbParm::input[ProbParm::input_H2_init + mp1 * ProbParm::nvars] * f1;
      massfrac[O2_ID] =
        ProbParm::input[ProbParm::input_O2_init + m * ProbParm::nvars] * f0 +
        ProbParm::input[ProbParm::input_O2_init + mp1 * ProbParm::nvars] * f1;
      massfrac[N2_ID] =
        ProbParm::input[ProbParm::input_N2_init + m * ProbParm::nvars] * f0 +
        ProbParm::input[ProbParm::input_N2_init + mp1 * ProbParm::nvars] * f1;
    } else if (ProbParm::case_type_int == 2) {
      massfrac[H2_ID] =
        ProbParm::input[ProbParm::input_H2 + m * ProbParm::nvars] * f0 +
        ProbParm::input[ProbParm::input_H2 + mp1 * ProbParm::nvars] * f1;
      massfrac[O2_ID] =
        ProbParm::input[ProbParm::input_O2 + m * ProbParm::nvars] * f0 +
        ProbParm::input[ProbParm::input_O2 + mp1 * ProbParm::nvars] * f1;
      massfrac[N2_ID] =
        ProbParm::input[ProbParm::input_N2 + m * ProbParm::nvars] * f0 +
        ProbParm::input[ProbParm::input_N2 + mp1 * ProbParm::nvars] * f1;
      massfrac[H2O_ID] =
        ProbParm::input[ProbParm::input_H2O + m * ProbParm::nvars] * f0 +
        ProbParm::input[ProbParm::input_H2O + mp1 * ProbParm::nvars] * f1;
      massfrac[H_ID] =
        ProbParm::input[ProbParm::input_H + m * ProbParm::nvars] * f0 +
        ProbParm::input[ProbParm::input_H + mp1 * ProbParm::nvars] * f1;
      massfrac[O_ID] =
        ProbParm::input[ProbParm::input_O + m * ProbParm::nvars] * f0 +
        ProbParm::input[ProbParm::input_O + mp1 * ProbParm::nvars] * f1;
      massfrac[OH_ID] =
        ProbParm::input[ProbParm::input_OH + m * ProbParm::nvars] * f0 +
        ProbParm::input[ProbParm::input_OH + mp1 * ProbParm::nvars] * f1;
      massfrac[HO2_ID] =
        ProbParm::input[ProbParm::input_HO2 + m * ProbParm::nvars] * f0 +
        ProbParm::input[ProbParm::input_HO2 + mp1 * ProbParm::nvars] * f1;
      massfrac[H2O2_ID] =
        ProbParm::input[ProbParm::input_H2O2 + m * ProbParm::nvars] * f0 +
        ProbParm::input[ProbParm::input_H2O2 + mp1 * ProbParm::nvars] * f1;
    }
#endif
    EOS::PYT2RE(p, massfrac, T, rho, eint);
  } else {
    massfrac[0] = 1.0;
    rho = ProbParm::rho0;
    T = ProbParm::T0;
    p = ProbParm::p0;
    EOS::RYP2E(rho, massfrac, p, eint);
  }

  // Set the state
  state(i, j, k, URHO) = rho;
  state(i, j, k, UMX) = rho * u[0];
  state(i, j, k, UMY) = rho * u[1];
  state(i, j, k, UMZ) = rho * u[2];
  state(i, j, k, UEINT) = rho * eint;
  state(i, j, k, UEDEN) =
    rho * (eint + 0.5 * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]));
  state(i, j, k, UTEMP) = T;
  for (int n = 0; n < NUM_SPECIES; n++)
    state(i, j, k, UFS + n) = rho * massfrac[n];
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real x[AMREX_SPACEDIM],
  const amrex::Real s_int[NVAR],
  amrex::Real s_ext[NVAR],
  const int idir,
  const int sgn,
  const amrex::Real time,
  amrex::GeometryData const& geomdata)
{
}

using ProblemTags = EmptyProbTagStruct;

struct MyProbDeriveStruct
{
  static void
  add(amrex::DeriveList& derive_lst, amrex::DescriptorList& desc_lst)
  {
    derive_lst.add(
      "z_vorticity", amrex::IndexType::TheCellType(), 1, pc_dervortz,
      grow_box_by_one);
    derive_lst.addComponent("z_vorticity", desc_lst, State_Type, 0, NVAR);

#ifdef PELEC_USE_REACTIONS
    derive_lst.add(
      "heat_release", amrex::IndexType::TheCellType(), 1, pc_derheatrelease,
      the_same_box);
    derive_lst.addComponent("heat_release", desc_lst, State_Type, 0, NVAR);
#endif
  }

  static void pc_dervortz(
    const amrex::Box& bx,
    amrex::FArrayBox& derfab,
    int dcomp,
    int /*ncomp*/,
    const amrex::FArrayBox& datfab,
    const amrex::Geometry& geomdata,
    amrex::Real /*time*/,
    const int* /*bcrec*/,
    int level)
  {
    auto const dat = datfab.array();
    auto vort = derfab.array();

    const amrex::Box& gbx = amrex::grow(bx, 1);

    amrex::FArrayBox local(gbx, 2);
    amrex::Elixir local_eli = local.elixir();
    auto larr = local.array();

    // Convert momentum to velocity.
    amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const amrex::Real rhoInv = 1.0 / dat(i, j, k, URHO);
      larr(i, j, k, 0) = dat(i, j, k, UMX) * rhoInv;
      larr(i, j, k, 1) = dat(i, j, k, UMY) * rhoInv;
    });

    const amrex::Real dx = geomdata.CellSize(0);
    const amrex::Real dy = geomdata.CellSize(1);
    const amrex::Real dz = geomdata.CellSize(2);

    // Calculate vorticity.
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const amrex::Real vx =
        0.5 * (larr(i + 1, j, k, 1) - larr(i - 1, j, k, 1)) / dx;
      const amrex::Real uy =
        0.5 * (larr(i, j + 1, k, 0) - larr(i, j - 1, k, 0)) / dy;
      vort(i, j, k) = vx - uy;
    });
  }

#ifdef PELEC_USE_REACTIONS
  static void pc_derheatrelease(
    const amrex::Box& bx,
    amrex::FArrayBox& derfab,
    int dcomp,
    int /*ncomp*/,
    const amrex::FArrayBox& datfab,
    const amrex::Geometry& geomdata,
    amrex::Real /*time*/,
    const int* /*bcrec*/,
    int level)
  {
    auto const dat = datfab.array();
    auto hr = derfab.array();

    // Calculate heat release.
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      amrex::Real rho = dat(i, j, k, URHO);
      amrex::Real T = dat(i, j, k, UTEMP);
      amrex::Real massfrac[NUM_SPECIES] = {0.0};
      const amrex::Real rhoInv = 1.0 / rho;
      for (int n = 0; n < NUM_SPECIES; n++)
        massfrac[n] = dat(i, j, k, UFS + n) * rhoInv;
      amrex::Real wdot[NUM_SPECIES] = {0.0};
      amrex::Real hi[NUM_SPECIES] = {0.0};
      EOS::RTY2WDOT(rho, T, massfrac, wdot);
      EOS::T2Hi(T, hi);

      hr(i, j, k) = 0.0;
      for (int n = 0; n < NUM_SPECIES; n++)
        hr(i, j, k) += wdot[n] * hi[n];
    });
  }
#endif
};

void pc_prob_close();

using ProblemDerives = MyProbDeriveStruct;

#endif
