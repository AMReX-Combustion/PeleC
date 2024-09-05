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
  // Standard state variable initialization
  const amrex::Real p = prob_parm.p;
  amrex::Real rho = 0.0, eint = 0.0;
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  for (int n = 0; n < NUM_SPECIES; n++) {
    massfrac[n] = prob_parm.massfrac[n];
  }
  auto eos = pele::physics::PhysicsType::eos();
  eos.PYT2RE(p, massfrac, prob_parm.T, rho, eint);

  state(i, j, k, URHO) = rho;
  state(i, j, k, UMX) = rho * prob_parm.U;
  state(i, j, k, UMY) = 0.0;
  state(i, j, k, UMZ) = 0.0;
  state(i, j, k, UEINT) = rho * eint;
  state(i, j, k, UEDEN) = rho * (eint + 0.5 * (prob_parm.U * prob_parm.U));
  state(i, j, k, UTEMP) = prob_parm.T;
  for (int n = 0; n < NUM_SPECIES; n++) {
    state(i, j, k, UFS + n) = rho * prob_parm.massfrac[n];
  }

  // Initial AUX quantity
  {
#if NUM_AUX > 0
    amrex::Real aux[NUM_AUX] = {0.0};

    // Initial state for aux quantities
    amrex::Real aux_0[NUM_AUX] = {0.0};
    for (int n = 0; n < NUM_AUX; n++) {
      aux_0[n] = (n + 1) * prob_parm.aux_IC;
    }

    // Current x,y,z locations
    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* dx = geomdata.CellSize();
    const amrex::Real x[AMREX_SPACEDIM] = {AMREX_D_DECL(
      prob_lo[0] + static_cast<amrex::Real>(i + 0.5) * dx[0],
      prob_lo[1] + static_cast<amrex::Real>(j + 0.5) * dx[1],
      prob_lo[2] + static_cast<amrex::Real>(k + 0.5) * dx[2])};

    // Create NUM_AUX boxes from (x_strt,y_strt) to (x_end,y_end)
    for (int n = 0; n < NUM_AUX; n++) {
      amrex::Real x_strt = prob_parm.aux_xy_lo + n * 4. * prob_parm.aux_length;
      amrex::Real x_end = x_strt + prob_parm.aux_length;
      amrex::Real y_strt = prob_parm.aux_xy_lo;
      amrex::Real y_end = y_strt + prob_parm.aux_height;
      if ((x[0] > x_strt) && (x[1] > y_strt)) {
        if ((x[0] < x_end) && (x[1] < y_end)) {
          aux[n] = aux_0[n];
        }
      }
    }
    // Add to state variable
    for (int n = 0; n < NUM_AUX; n++) {
      state(i, j, k, UFX + n) = aux[n];
    }
#endif
  }

  // Initial ADV quantity
  {
#if NUM_ADV > 0
    amrex::Real adv[NUM_ADV] = {0.0};

    // Initial state for adv quantities
    amrex::Real adv_0[NUM_ADV] = {0.0};
    for (int n = 0; n < NUM_ADV; n++) {
      adv_0[n] = (n + 1) * prob_parm.aux_IC;
    }

    // Current locations
    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* dx = geomdata.CellSize();
    const amrex::Real x[AMREX_SPACEDIM] = {AMREX_D_DECL(
      prob_lo[0] + static_cast<amrex::Real>(i + 0.5) * dx[0],
      prob_lo[1] + static_cast<amrex::Real>(j + 0.5) * dx[1],
      prob_lo[2] + static_cast<amrex::Real>(k + 0.5) * dx[2])};

    amrex::Real adv_xy_lo = prob_parm.aux_xy_lo + prob_parm.aux_length;
    for (int n = 0; n < NUM_ADV; n++) {
      // Create NUM_ADV boxes from (x_strt,y_strt) to (x_end,y_end)
      amrex::Real x_strt = adv_xy_lo + n * 4. * prob_parm.aux_length;
      amrex::Real x_end = x_strt + prob_parm.aux_length;
      amrex::Real y_strt = prob_parm.aux_xy_lo;
      amrex::Real y_end = y_strt + prob_parm.aux_height;
      if ((x[0] > x_strt) && (x[1] > y_strt)) {
        if ((x[0] < x_end) && (x[1] < y_end)) {
          adv[n] = adv_0[n];
        }
      }
    }

    // Add to state variable
    for (int n = 0; n < NUM_ADV; n++) {
      state(i, j, k, UFA + n) = rho * adv[n];
    }
#endif
  }
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
  ProbParmDevice const& prob_parm,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& /*turb_fluc*/)
{
  if (idir == 0) {

    amrex::Real rho = 0.0, u = 0.0, v = 0.0, w = 0.0, eint = 0.0, T = 0.0;
    auto eos = pele::physics::PhysicsType::eos();

    if (sgn == 1) {

      T = prob_parm.T;
      amrex::Real p_inter = 0.0;
      amrex::Real massfrac_inter[NUM_SPECIES] = {0.0};
      for (int n = 0; n < NUM_SPECIES; n++) {
        massfrac_inter[n] = s_inter[UFS + n] / s_inter[URHO];
      }

      // Ghost state p_d = p_inter (dp/dx = 0)
      eos.RTY2P(s_inter[URHO], s_inter[UTEMP], massfrac_inter, p_inter);

      // Ghost state rho and eint (constant T)
      eos.PYT2RE(p_inter, massfrac_inter, T, rho, eint);

      // Ghost state for velocity
      u = prob_parm.U;
      v = 0.0;
      w = 0.0;

    } else if (sgn == -1) {

      // Following Blazek p 279, eq. 8.23

      // Interior state (point d)
      const amrex::Real* prob_hi = geomdata.ProbHi();
      const amrex::Real* dx = geomdata.CellSize();
      const amrex::Real xd = prob_hi[0] - 0.5 * dx[0];
      const amrex::Real rho_inter = s_inter[URHO];
      const amrex::Real u_inter = s_inter[UMX] / rho_inter;
      const amrex::Real v_inter = s_inter[UMY] / rho_inter;
      const amrex::Real w_inter = s_inter[UMZ] / rho_inter;
      const amrex::Real T_inter = s_inter[UTEMP];
      amrex::Real p_inter = 0.0, cs_inter = 0.0;
      amrex::Real massfrac[NUM_SPECIES] = {0.0};
      for (int n = 0; n < NUM_SPECIES; n++) {
        massfrac[n] = s_inter[UFS + n] / s_inter[URHO];
      }
      eos.RTY2P(rho_inter, T_inter, massfrac, p_inter);
      eos.RTY2Cs(rho_inter, T_inter, massfrac, cs_inter);

      // Boundary state (point b)
      const amrex::Real xb = prob_hi[0];
      const amrex::Real pb = prob_parm.p;
      const amrex::Real rhob =
        s_inter[URHO] + (pb - p_inter) / (cs_inter * cs_inter);
      const amrex::Real ub = u_inter + (p_inter - pb) / (rho_inter * cs_inter);
      const amrex::Real vb = v_inter;
      const amrex::Real wb = w_inter;

      // Ghost state (point a). Linear extrapolation from d and b
      rho = (rhob - rho_inter) / (xb - xd) * (x[0] - xd) + rho_inter;
      const amrex::Real p = (pb - p_inter) / (xb - xd) * (x[0] - xd) + p_inter;

      eos.RYP2E(rho, massfrac, p, eint);
      eos.EY2T(eint, massfrac, T);

      u = (ub - u_inter) / (xb - xd) * (x[0] - xd) + u_inter;
      v = (vb - v_inter) / (xb - xd) * (x[0] - xd) + v_inter;
      w = (wb - w_inter) / (xb - xd) * (x[0] - xd) + w_inter;
    }

    s_ext[URHO] = rho;
    s_ext[UMX] = rho * u;
    s_ext[UMY] = rho * v;
    s_ext[UMZ] = rho * w;
    s_ext[UEINT] = rho * eint;
    s_ext[UEDEN] = rho * (eint + 0.5 * (u * u + v * v + w * w));
    s_ext[UTEMP] = T;
    for (int n = 0; n < NUM_SPECIES; n++) {
      s_ext[UFS + n] = rho * prob_parm.massfrac[n];
    }
  }
}

void pc_prob_close();

struct MyProblemSpecificFunctions : public DefaultProblemSpecificFunctions
{
  static void set_aux_names(amrex::Vector<std::string>& a_aux_names)
  {
    a_aux_names.resize(NUM_AUX);
#if NUM_AUX > 0
    for (int n = 0; n < NUM_AUX; ++n) {
      a_aux_names[n] = "B_" + std::to_string(n);
    }
#endif
  }

  static void set_adv_names(amrex::Vector<std::string>& a_adv_names)
  {
    a_adv_names.resize(NUM_ADV);
#if NUM_ADV > 0
    for (int n = 0; n < NUM_ADV; ++n) {
      a_adv_names[n] = "A_" + std::to_string(n);
    }
#endif
  }

  static void problem_modify_ext_sources(
    amrex::Real /*time*/,
    amrex::Real /*dt*/,
    const amrex::MultiFab& state_old,
    const amrex::MultiFab& state_new,
    amrex::MultiFab& ext_src,
    int /*ng*/,
    amrex::GeometryData const& /*geomdata*/,
    ProbParmDevice const& prob_parm)
  {
    /* Notes: ext_src contains sources from velocity forcing coming in
       This function should add to rather than overwrite ext_src.
    */

    const amrex::IntVect ngs(0);
    auto const& Sns = state_new.const_arrays();
    auto const& Farrs = ext_src.arrays();

    amrex::ParallelFor(
      state_old, ngs,
      [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
        // Source terms for the auxiliarry quantities
        amrex::Real Bn = 0.0;
        for (int n = 0; n < NUM_AUX; n++) {
          Bn = Sns[nbx](i, j, k, UFX + n);
          Farrs[nbx](i, j, k, UFX + n) += prob_parm.aux_srcstrength * Bn;
        }
      });
    amrex::Gpu::synchronize();
  }
};
using ProblemSpecificFunctions = MyProblemSpecificFunctions;
#endif
