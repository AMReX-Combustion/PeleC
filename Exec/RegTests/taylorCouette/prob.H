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
#include "EOS.H"

#include "Tagging.H"
#include "ProblemDerive.H"
#include "prob_parm.H"
#include "EB.H"
#include <AMReX_EB2_IF_Rotation.H>
#include "Utilities.H"
#include "Geometry.H"

class EBConcCylinders : public pele::pelec::Geometry::Register<EBConcCylinders>
{
public:
  static std::string identifier() { return "conc_cylinders"; }

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
  amrex::GeometryData const& geomdata,
  ProbParmDevice const& prob_parm)
{
  // Geometry
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* dx = geomdata.CellSize();

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

  if (prob_parm.do_rf) {
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> axis_loc = {
      prob_parm.rf_axis_x, prob_parm.rf_axis_y, prob_parm.rf_axis_z};
    amrex::Real rotenrg = get_rot_energy(
      i, j, k, prob_parm.rf_omega, prob_parm.rf_axis, axis_loc, geomdata);
    state(i, j, k, UEDEN) -= rho * rotenrg;
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
  ProbParmDevice const& /*prob_parm*/,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& /*turb_fluc*/)
{
}

void pc_prob_close();

using ProblemTags = EmptyProbTagStruct;
using ProblemDerives = EmptyProbDeriveStruct;

/*using ProblemTags = EmptyProbTagStruct;

  struct MyProbDeriveStruct
  {
  static amrex::Box the_same_box(const amrex::Box& b) { return b; }

  static void
  add(amrex::DeriveList& derive_lst, const amrex::DescriptorList& desc_lst)
  {
// Add derives as follows and define the derive function below:
derive_lst.add(
"velx_if", amrex::IndexType::TheCellType(), 1, pc_velx_if,
the_same_box);
derive_lst.addComponent("velx_if", desc_lst, State_Type, 0, NVAR);

derive_lst.add(
"vely_if", amrex::IndexType::TheCellType(), 1, pc_vely_if,
the_same_box);
derive_lst.addComponent("vely_if", desc_lst, State_Type, 0, NVAR);

derive_lst.add(
"velz_if", amrex::IndexType::TheCellType(), 1, pc_velz_if,
the_same_box);
derive_lst.addComponent("velz_if", desc_lst, State_Type, 0, NVAR);
}

static void pc_rhoerror(
const amrex::Box& bx,
amrex::FArrayBox& derfab,
int dcomp,
int ncomp,
const amrex::FArrayBox& datfab,
const amrex::Geometry& geomdata,
amrex::Real time,
const int* bcrec,
int level)
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

using ProblemDerives = MyProbDeriveStruct;*/

#endif
