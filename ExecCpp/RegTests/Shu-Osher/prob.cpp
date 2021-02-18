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
  const amrex_real* problo,
  const amrex_real* probhi)
{
  // Parse params
  amrex::ParmParse pp("prob");
  pp.query("p_l", PeleC::prob_parm_device->p_l);
  pp.query("u_l", PeleC::prob_parm_device->u_l);
  pp.query("rho_l", PeleC::prob_parm_device->rho_l);
  pp.query("T_l", PeleC::prob_parm_device->T_l);
  pp.query("p_r", PeleC::prob_parm_device->p_r);
  pp.query("u_r", PeleC::prob_parm_device->u_r);
  pp.query("rho_r_base", PeleC::prob_parm_device->rho_r_base);
  pp.query("rho_r_amp", PeleC::prob_parm_device->rho_r_amp);
  pp.query("rho_r_osc", PeleC::prob_parm_device->rho_r_osc);
  pp.query("T_r", PeleC::prob_parm_device->T_r);
  pp.query("frac", PeleC::prob_parm_device->frac);
  pp.query("idir", PeleC::prob_parm_device->idir);

  for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
    PeleC::prob_parm_device->split[idir] =
      PeleC::prob_parm_device->frac * (problo[idir] + probhi[idir]);
  }

  amrex::Real e_l;
  amrex::Real e_r /*, cs, cp */;
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  massfrac[0] = 1.0;

  EOS::RYP2E(
    PeleC::prob_parm_device->rho_l, massfrac, PeleC::prob_parm_device->p_l,
    e_l);
  EOS::EY2T(e_l, massfrac, PeleC::prob_parm_device->T_l);
  PeleC::prob_parm_device->rhoe_l = PeleC::prob_parm_device->rho_l * e_l;
  EOS::RYP2E(
    PeleC::prob_parm_device->rho_r_base, massfrac, PeleC::prob_parm_device->p_r,
    e_r);
  EOS::EY2T(e_r, massfrac, PeleC::prob_parm_device->T_r);
  PeleC::prob_parm_device->rhoe_r = PeleC::prob_parm_device->rho_r_base * e_r;
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
