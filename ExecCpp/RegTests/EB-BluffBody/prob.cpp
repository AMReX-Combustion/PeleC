#include "prob.H"

namespace ProbParm {
AMREX_GPU_DEVICE_MANAGED amrex::Real p = 1013250.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real T = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real rho = 0.00116;
AMREX_GPU_DEVICE_MANAGED amrex::Real eint = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real vx_in = 9000.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real vy_in = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real Re_L = 2500.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real Pr = 0.7;
AMREX_GPU_DEVICE_MANAGED amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {
  1.0};
} // namespace ProbParm

void
pc_prob_close()
{
}

extern "C" {
void
amrex_probinit(
  const int* /*init*/,
  const int* /*name*/,
  const int* /*namelen*/,
  const amrex_real* problo,
  const amrex_real* probhi)
{
  // Parse params
  {
    amrex::ParmParse pp("prob");
    pp.query("p", ProbParm::p);
    pp.query("rho", ProbParm::rho);
    pp.query("vx_in", ProbParm::vx_in);
    pp.query("vy_in", ProbParm::vy_in);
    pp.query("Re_L", ProbParm::Re_L);
    pp.query("Pr", ProbParm::Pr);
  }

  amrex::Real L = (probhi[0] - problo[0]) * 0.2;

  amrex::Real cp = 0.0;
  ProbParm::massfrac[0] = 1.0;
  EOS::RYP2E(
    ProbParm::rho, ProbParm::massfrac.begin(), ProbParm::p, ProbParm::eint);
  EOS::EY2T(ProbParm::eint, ProbParm::massfrac.begin(), ProbParm::T);
  EOS::TY2Cp(ProbParm::T, ProbParm::massfrac.begin(), cp);

  TransParm trans_parm;

  // Default
  trans_parm.const_viscosity = 0.0;
  trans_parm.const_bulk_viscosity = 0.0;
  trans_parm.const_conductivity = 0.0;
  trans_parm.const_diffusivity = 0.0;

  // User-specified
  {
    amrex::ParmParse pp("transport");
    pp.query("const_viscosity", trans_parm.const_viscosity);
    pp.query("const_bulk_viscosity", trans_parm.const_bulk_viscosity);
    pp.query("const_conductivity", trans_parm.const_conductivity);
    pp.query("const_diffusivity", trans_parm.const_diffusivity);
  }

  trans_parm.const_bulk_viscosity = 0.0;
  trans_parm.const_diffusivity = 0.0;
  trans_parm.const_viscosity =
    ProbParm::rho * ProbParm::vx_in * L / ProbParm::Re_L;
  trans_parm.const_conductivity =
    trans_parm.const_viscosity * cp / ProbParm::Pr;

#ifdef AMREX_USE_GPU
  amrex::Gpu::htod_memcpy(trans_parm_g, &trans_parm, sizeof(trans_parm));
#else
  std::memcpy(trans_parm_g, &trans_parm, sizeof(trans_parm));
#endif
}
}

void
PeleC::problem_post_timestep()
{
}

void
PeleC::problem_post_init()
{
}

void
PeleC::problem_post_restart()
{
}
