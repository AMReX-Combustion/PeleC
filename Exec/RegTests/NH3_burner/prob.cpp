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
  // Parse params
  {
    amrex::ParmParse pp("prob");
    pp.query("P_mean",   PeleC::h_prob_parm_device->pamb);
    pp.query("v_in",  PeleC::h_prob_parm_device->v_in);
    pp.query("T_in",  PeleC::h_prob_parm_device->T_in);

    pp.query("T_wall1", PeleC::h_prob_parm_device->T_wall1);
    pp.query("z_wall1", PeleC::h_prob_parm_device->z_wall1);
    pp.query("T_wall2", PeleC::h_prob_parm_device->T_wall2);
    pp.query("z_wall2", PeleC::h_prob_parm_device->z_wall2);
    pp.query("T_wall3", PeleC::h_prob_parm_device->T_wall3);
    pp.query("z_wall3", PeleC::h_prob_parm_device->z_wall3);
    pp.query("T_wall4", PeleC::h_prob_parm_device->T_wall4);
    pp.query("z_wall4", PeleC::h_prob_parm_device->z_wall4);
    pp.query("T_wall5", PeleC::h_prob_parm_device->T_wall5);
    pp.query("z_wall5", PeleC::h_prob_parm_device->z_wall5);
    
    amrex::ParmParse ppEB("EB");
    ppEB.query("in_diam",PeleC::h_prob_parm_device->inflow_diam);

  }
   

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
