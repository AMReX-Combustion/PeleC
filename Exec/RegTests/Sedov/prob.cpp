#include "prob.H"

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
  const amrex::Real* /*problo*/,
  const amrex::Real* /*probhi*/)
{
  // Parse params
  amrex::ParmParse pp("prob");
  pp.query("p_ambient", PeleC::h_prob_parm_device->p_ambient);
  pp.query("dens_ambient", PeleC::h_prob_parm_device->dens_ambient);
  pp.query("exp_energy", PeleC::h_prob_parm_device->exp_energy);
  pp.query("r_init", PeleC::h_prob_parm_device->r_init);
  pp.query("nsub", PeleC::h_prob_parm_device->nsub);
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
