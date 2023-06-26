#include "prob.H"

void
pc_prob_close()
{
}

void
parse_params(ProbParmDevice* prob_parm_device)
{
  // Parse params
  amrex::ParmParse pp("prob");
  pp.query("P_mean", prob_parm_device->P_mean);
  pp.query("T_mean", prob_parm_device->T_mean);
  pp.query("rvort", prob_parm_device->rvort);
  pp.query("xvort", prob_parm_device->xvort);
  pp.query("yvort", prob_parm_device->yvort);
  pp.query("forcevort", prob_parm_device->forcevort);
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
  parse_params(PeleC::h_prob_parm_device);
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

void
EBConvergingNozzle::build(
  const amrex::Geometry& geom, const int max_coarsening_level)
{
}
