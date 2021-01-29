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
  const amrex_real* /*problo*/,
  const amrex_real* /*probhi*/)
{
  // Parse params
  amrex::ParmParse pp("prob");
  pp.query("p_init", PeleC::prob_parm_device->p_init);
  pp.query("Y_init_H2", PeleC::prob_parm_device->Y_init_H2);
  pp.query("Y_init_O2", PeleC::prob_parm_device->Y_init_O2);
  pp.query("Y_init_N2", PeleC::prob_parm_device->Y_init_N2);
  pp.query("T_init", PeleC::prob_parm_device->T_init);

  // Initial values
  PeleC::prob_parm_device->massfrac[H2_ID] = PeleC::prob_parm_device->Y_init_H2;
  PeleC::prob_parm_device->massfrac[O2_ID] = PeleC::prob_parm_device->Y_init_O2;
  PeleC::prob_parm_device->massfrac[N2_ID] = PeleC::prob_parm_device->Y_init_N2;
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
