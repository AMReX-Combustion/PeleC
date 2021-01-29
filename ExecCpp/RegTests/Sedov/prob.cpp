#include "prob.H"

namespace ProbParm {
AMREX_GPU_DEVICE_MANAGED amrex::Real p_ambient;
AMREX_GPU_DEVICE_MANAGED amrex::Real dens_ambient;
AMREX_GPU_DEVICE_MANAGED amrex::Real exp_energy;
AMREX_GPU_DEVICE_MANAGED amrex::Real r_init;
AMREX_GPU_DEVICE_MANAGED int nsub;
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
  pp.query("p_ambient", PeleC::prob_parm_device->p_ambient);
  pp.query("dens_ambient", PeleC::prob_parm_device->dens_ambient);
  pp.query("exp_energy", PeleC::prob_parm_device->exp_energy);
  pp.query("r_init", PeleC::prob_parm_device->r_init);
  pp.query("nsub", PeleC::prob_parm_device->nsub);
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
