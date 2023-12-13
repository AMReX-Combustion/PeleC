#include "prob.H"

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
  pp.query("init_v", PeleC::h_prob_parm_device->v0);
  pp.get("ref_p", PeleC::h_prob_parm_device->p0);
  pp.get("ref_T", PeleC::h_prob_parm_device->T0);
  if (pp.contains("X_N2")) {
    pp.query("X_N2", PeleC::h_prob_parm_device->YX_N2);
    pp.query("X_O2", PeleC::h_prob_parm_device->YX_O2);
    pp.query("X_H2O", PeleC::h_prob_parm_device->YX_H2O);
    pp.query("X_CO2", PeleC::h_prob_parm_device->YX_CO2);
    PeleC::h_prob_parm_device->mol_fracs = true;
  } else {
    pp.query("Y_N2", PeleC::h_prob_parm_device->YX_N2);
    pp.query("Y_O2", PeleC::h_prob_parm_device->YX_O2);
    pp.query("Y_H2O", PeleC::h_prob_parm_device->YX_H2O);
    pp.query("Y_CO2", PeleC::h_prob_parm_device->YX_CO2);
  }

}
}

void
pc_prob_close()
{
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
