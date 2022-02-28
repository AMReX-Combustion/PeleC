#ifndef _PROB_H_
#define _PROB_H_

#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>
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
#include "Utilities.H"
#include "turbinflow.H"

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
  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real massfrac[NUM_SPECIES] = {0.0};

  for (int n = 0; n < NUM_SPECIES; n++) {
    massfrac[n] = prob_parm.Y_ox[n];
  }

  amrex::Real x_velocity = 0.0;
  amrex::Real y_velocity = 0.0;
  amrex::Real z_velocity = 0.0;

  amrex::Real rho = 0.0, eint = 0;
  eos.PYT2R(prob_parm.P_mean, massfrac, prob_parm.T_ox, rho);
  eos.RTY2E(rho, prob_parm.T_ox, massfrac, eint);

  // Set the state
  state(i, j, k, URHO) = rho;
  state(i, j, k, UMX) = rho * x_velocity;
  state(i, j, k, UMY) = rho * y_velocity;
  state(i, j, k, UMZ) = rho * z_velocity;
  state(i, j, k, UEINT) = rho * eint;
  state(i, j, k, UEDEN) =
    rho * (eint + 0.5 * (x_velocity * x_velocity + y_velocity * y_velocity +
                         z_velocity * z_velocity));
  state(i, j, k, UTEMP) = prob_parm.T_ox;
  for (int n = 0; n < NUM_SPECIES; n++) {
    state(i, j, k, UFS + n) = rho * massfrac[n];
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
time_func(
  amrex::Real time, amrex::Real start, amrex::Real duration, amrex::Real tau)
{
  amrex::Real eta_time;
  eta_time = 0.5 * (std::tanh((time - start) / (duration * tau)) *
                      std::tanh((-time + start + duration) / (duration * tau)) +
                    1.0);
  eta_time = std::min(std::max(eta_time, 0.0), 1.0);
  return eta_time;
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
  amrex::GeometryData const& geomdata,
  ProbParmDevice const& prob_parm)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();
  amrex::Real velx, vely, velz;
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  auto eos = pele::physics::PhysicsType::eos();

  amrex::Real radius = {0.0};
  amrex::Real int_energy;

  if (sgn == 1 && idir == 2) {

    if (prob_parm.turbulence) {
      velx = s_ext[UMX];
      vely = s_ext[UMY];
      velz = s_ext[UMZ] + prob_parm.v_in;
    } else {
      velx = 0.0;
      vely = 0.0;
      velz = 0.0 + prob_parm.v_in;
    }

    radius = std::sqrt(
      std::pow(x[0] - prob_parm.center_xy[0], 2) +
      std::pow(x[1] - prob_parm.center_xy[1], 2));
    amrex::Real eta = {0};
    eta =
      0.5 *
      (1.0 - std::tanh(
               2.0 * (radius - prob_parm.D * 0.5) / (0.1 * prob_parm.D * 0.5)));
    eta *=
      time_func(time, prob_parm.inj_start, prob_parm.inj_dur, prob_parm.tau);
    for (int n = 0; n < NUM_SPECIES; n++) {
      massfrac[n] = eta * prob_parm.Y_fuel[n] + (1.0 - eta) * prob_parm.Y_ox[n];
    }
    int_energy = eta * prob_parm.U_fuel + (1.0 - eta) * prob_parm.U_ox;

    s_ext[UTEMP] = 500.0; // initial guess to help EY2T converge
    eos.EY2T(int_energy, massfrac, s_ext[UTEMP]);
    //    GET_T_GIVEN_HY(&enthalpy,massfrac,&s_ext[DEF_Temp],&ierr);

    amrex::Real rho, eint;
    eos.PYT2R(prob_parm.P_mean, massfrac, s_ext[UTEMP], rho);
    s_ext[URHO] = rho;

    eos.RTY2E(rho, s_ext[UTEMP], massfrac, eint);
    s_ext[UEINT] = rho * eint;

    for (int n = 0; n < NUM_SPECIES; n++) {
      s_ext[UFS + n] = rho * massfrac[n];
    }
    // Velocity
    s_ext[UMX] = rho * velx * eta;
    s_ext[UMY] = rho * vely * eta;
    s_ext[UMZ] = rho * velz * eta;
    s_ext[UEDEN] = s_ext[UEINT] + rho * eta * eta * 0.5 *
                                    (velx * velx + vely * vely + velz * velz);
  }
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
    ProbParmDevice /*const& d_prob_parm_device*/) noexcept
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
    // auto const dat = datfab.array();
    // auto arr = derfab.array();
    // amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    // { do something with arr
    // });
  }
};

void pc_prob_close();

using ProblemDerives = MyProbDeriveStruct;

#endif
