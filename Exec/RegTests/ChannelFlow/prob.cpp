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
  pp.query("T_mean", PeleC::h_prob_parm_device->T_mean);
  pp.query("P_mean", PeleC::h_prob_parm_device->P_mean);
  pp.query("dpdx", PeleC::h_prob_parm_device->dpdx);
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
