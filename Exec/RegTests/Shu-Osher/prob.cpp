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
  amrex::ParmParse pp("prob");
  pp.query("p_l", PeleC::h_prob_parm_device->p_l);
  pp.query("u_l", PeleC::h_prob_parm_device->u_l);
  pp.query("rho_l", PeleC::h_prob_parm_device->rho_l);
  pp.query("T_l", PeleC::h_prob_parm_device->T_l);
  pp.query("p_r", PeleC::h_prob_parm_device->p_r);
  pp.query("u_r", PeleC::h_prob_parm_device->u_r);
  pp.query("rho_r_base", PeleC::h_prob_parm_device->rho_r_base);
  pp.query("rho_r_amp", PeleC::h_prob_parm_device->rho_r_amp);
  pp.query("rho_r_osc", PeleC::h_prob_parm_device->rho_r_osc);
  pp.query("T_r", PeleC::h_prob_parm_device->T_r);
  pp.query("frac", PeleC::h_prob_parm_device->frac);
  pp.query("idir", PeleC::h_prob_parm_device->idir);

  for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
    PeleC::h_prob_parm_device->split[idir] =
      PeleC::h_prob_parm_device->frac * (problo[idir] + probhi[idir]);
  }

  amrex::Real e_l;
  amrex::Real e_r /*, cs, cp */;
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  massfrac[0] = 1.0;

  auto eos = pele::physics::PhysicsType::eos();
  eos.RYP2E(
    PeleC::h_prob_parm_device->rho_l, massfrac, PeleC::h_prob_parm_device->p_l,
    e_l);
  eos.EY2T(e_l, massfrac, PeleC::h_prob_parm_device->T_l);
  PeleC::h_prob_parm_device->rhoe_l = PeleC::h_prob_parm_device->rho_l * e_l;
  eos.RYP2E(
    PeleC::h_prob_parm_device->rho_r_base, massfrac,
    PeleC::h_prob_parm_device->p_r, e_r);
  eos.EY2T(e_r, massfrac, PeleC::h_prob_parm_device->T_r);
  PeleC::h_prob_parm_device->rhoe_r =
    PeleC::h_prob_parm_device->rho_r_base * e_r;
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
