#include "prob.H"

namespace ProbParm {
AMREX_GPU_DEVICE_MANAGED amrex::Real p0 = 1.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real rho0 = 1.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real u0 = 1.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real A = 0.2;
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
  const amrex_real* /*problo*/,
  const amrex_real* /*probhi*/)
{
  ProbParm::massfrac[0] = 1.0;
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
