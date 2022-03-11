#ifndef _PROB_H_
#define _PROB_H_

#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>
#include <AMReX_ParmParse.H>

#include "IndexDefines.H"
#include "Constants.H"
#include "EOS.H"
#include "Tagging.H"
#include "Transport.H"
#include "ProblemDerive.H"
#include "prob_parm.H"
#include "EB.H"
#include <AMReX_EB2_IF_Rotation.H>
#include "Utilities.H"
#include "Geometry.H"

void parse_params(std::unique_ptr<ProbParmDevice>& prob_parm_device);

class EBConvergingNozzle
  : public pele::pelec::Geometry::Register<EBConvergingNozzle>
{
public:
  static std::string identifier() { return "converging-nozzle"; }

  void
  build(const amrex::Geometry& geom, const int max_coarsening_level) override;
};

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_initdata(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& state,
  amrex::GeometryData const& /*geomdata*/,
  ProbParmDevice const& prob_parm)
{

  const amrex::Real u = 0.0;
  const amrex::Real v = 0.0;
  const amrex::Real w = 0.0;
  const amrex::Real p = prob_parm.p0;
  const amrex::Real T = prob_parm.T0;
  amrex::Real massfrac[NUM_SPECIES] = {1.0};
  amrex::Real rho, eint;
  auto eos = pele::physics::PhysicsType::eos();
  eos.PYT2RE(p, massfrac, T, rho, eint);

  state(i, j, k, URHO) = rho;
  state(i, j, k, UMX) = rho * u;
  state(i, j, k, UMY) = rho * v;
  state(i, j, k, UMZ) = rho * w;
  state(i, j, k, UEINT) = rho * eint;
  state(i, j, k, UEDEN) = rho * (eint + 0.5 * (u * u + v * v + w * w));
  state(i, j, k, UTEMP) = T;
  for (int n = 0; n < NUM_SPECIES; n++) {
    state(i, j, k, UFS + n) = rho * massfrac[n];
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real* /*x[AMREX_SPACEDIM]*/,
  const amrex::Real s_int[NVAR],
  amrex::Real s_ext[NVAR],
  const int idir,
  const int sgn,
  const amrex::Real /*time*/,
  amrex::GeometryData const& /*geomdata*/,
  ProbParmDevice const& prob_parm)
{

  if ((idir == 0) && (sgn == 1)) {

    s_ext[URHO] = s_int[URHO];
    s_ext[UMX] = -s_int[UMX];
    s_ext[UMY] = -s_int[UMY];
    s_ext[UMZ] = -s_int[UMZ];
    s_ext[UEINT] = s_int[UEINT];
    s_ext[UEDEN] = s_int[UEDEN];
    s_ext[UTEMP] = s_int[UTEMP];
    for (int n = 0; n < NUM_SPECIES; n++) {
      s_ext[UFS + n] = s_int[UFS + n];
    }

    // Swirler conditions
    const amrex::Real u = prob_parm.u_swirl_ax;
    const amrex::Real v = 0.0;
    const amrex::Real w = 0.0;

    amrex::Real massfrac[NUM_SPECIES] = {1.0};
    auto eos = pele::physics::PhysicsType::eos();
    amrex::Real p_int = 0.0;
    eos.RTY2P(s_int[URHO], s_int[UTEMP], massfrac, p_int);

    // subsonic velocity inlet
    const amrex::Real T = prob_parm.T_swirl;
    amrex::Real rho, eint;
    amrex::Real swirl_massfrac[NUM_SPECIES] = {1.0};
    eos.PYT2RE(p_int, swirl_massfrac, T, rho, eint);

    s_ext[URHO] = rho;
    s_ext[UMX] = rho * u;
    s_ext[UMY] = rho * v;
    s_ext[UMZ] = rho * w;
    s_ext[UEINT] = rho * eint;
    s_ext[UEDEN] = rho * (eint + 0.5 * (u * u + v * v + w * w));
    s_ext[UTEMP] = T;
    for (int n = 0; n < NUM_SPECIES; n++) {
      s_ext[UFS + n] = rho * swirl_massfrac[n];
    }
  }
}

void pc_prob_close();

using ProblemTags = EmptyProbTagStruct;
using ProblemDerives = EmptyProbDeriveStruct;

#endif
