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
  const amrex::Real* problo,
  const amrex::Real* probhi)
{
  amrex::ParmParse pp("prob");
  pp.query("inject_fuel", PeleC::h_prob_parm_device->inject_fuel);
  pp.query("centx", PeleC::h_prob_parm_device->centx);
  pp.query("centz", PeleC::h_prob_parm_device->centz);
  pp.query("r_hole", PeleC::h_prob_parm_device->r_hole);
  pp.query("init_type", PeleC::h_prob_parm_device->init_type);
  pp.query("cavity_depth", PeleC::h_prob_parm_device->cavity_depth);
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
}
