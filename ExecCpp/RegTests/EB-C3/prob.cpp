#include "prob.H"

namespace ProbParm {
AMREX_GPU_DEVICE_MANAGED amrex::Real p_init = 1013250.0; // 1 atm
AMREX_GPU_DEVICE_MANAGED amrex::Real Y_init_H2 = 0.06;
AMREX_GPU_DEVICE_MANAGED amrex::Real Y_init_O2 = 0.5;
AMREX_GPU_DEVICE_MANAGED amrex::Real Y_init_N2 = 0.44;
AMREX_GPU_DEVICE_MANAGED amrex::Real T_init = 940.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real rho_init = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real e_init = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {
  0.0};
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
  const amrex_real* /*problo*/,
  const amrex_real* /*probhi*/)
{
  // Parse params
  amrex::ParmParse pp("prob");
  pp.query("p_init", ProbParm::p_init);
  pp.query("Y_init_H2", ProbParm::Y_init_H2);
  pp.query("Y_init_O2", ProbParm::Y_init_O2);
  pp.query("Y_init_N2", ProbParm::Y_init_N2);
  pp.query("T_init", ProbParm::T_init);

  // Initial values
  ProbParm::massfrac[H2_ID] = ProbParm::Y_init_H2;
  ProbParm::massfrac[O2_ID] = ProbParm::Y_init_O2;
  ProbParm::massfrac[N2_ID] = ProbParm::Y_init_N2;
  EOS::PYT2RE(
    ProbParm::p_init, ProbParm::massfrac.begin(), ProbParm::T_init,
    ProbParm::rho_init, ProbParm::e_init);
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
