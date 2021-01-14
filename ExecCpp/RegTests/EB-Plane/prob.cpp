#include "prob.H"

namespace ProbParm {
AMREX_GPU_DEVICE_MANAGED amrex::Real p = 1013250.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real T = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real rho = 0.00116;
AMREX_GPU_DEVICE_MANAGED amrex::Real eint = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real vx_in = 0.0;
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
  amrex::ParmParse pp("prob");
  pp.query("p", ProbParm::p);
  pp.query("rho", ProbParm::rho);
  pp.query("vx_in", ProbParm::vx_in);
  pp.query("vy_in", ProbParm::vy_in);
  pp.query("Re_L", ProbParm::Re_L);
  pp.query("Pr", ProbParm::Pr);

  amrex::Real L = (probhi[0] - problo[0]) * 0.2;

  amrex::Real cp = 0.0;
  ProbParm::massfrac[0] = 1.0;
  EOS::RYP2E(
    ProbParm::rho, ProbParm::massfrac.begin(), ProbParm::p, ProbParm::eint);
  EOS::EY2T(ProbParm::eint, ProbParm::massfrac.begin(), ProbParm::T);
  EOS::TY2Cp(ProbParm::T, ProbParm::massfrac.begin(), cp);

  transport_params::const_bulk_viscosity = 0.0;
  transport_params::const_diffusivity = 0.0;
  transport_params::const_viscosity =
    ProbParm::rho * ProbParm::vx_in * L / ProbParm::Re_L;
  transport_params::const_conductivity =
    transport_params::const_viscosity * cp / ProbParm::Pr;
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
