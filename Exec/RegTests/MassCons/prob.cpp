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
  pp.query("u0", PeleC::h_prob_parm_device->u0);
  pp.query("v0", PeleC::h_prob_parm_device->v0);
  pp.query("w0", PeleC::h_prob_parm_device->w0);
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
