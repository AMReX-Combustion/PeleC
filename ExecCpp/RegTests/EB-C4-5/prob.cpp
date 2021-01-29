#include "prob.H"

namespace ProbParm {
AMREX_GPU_DEVICE_MANAGED amrex::Real p_init = 1013250.0; // 1 atm
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
  pp.query("p_init", PeleC::prob_parm_device->p_init);
  pp.query("T_init", PeleC::prob_parm_device->T_init);

  // Initial values
  PeleC::prob_parm_device->massfrac[0] = 1.0;
  EOS::PYT2RE(
    PeleC::prob_parm_device->p_init, PeleC::prob_parm_device->massfrac.begin(),
    PeleC::prob_parm_device->T_init, PeleC::prob_parm_device->rho_init,
    PeleC::prob_parm_device->e_init);
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
