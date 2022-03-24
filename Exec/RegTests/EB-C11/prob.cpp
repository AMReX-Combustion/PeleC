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
  pp.query("p_l", PeleC::h_prob_parm_device->p_l);
  pp.query("u_l", PeleC::h_prob_parm_device->u_l);
  pp.query("rho_l", PeleC::h_prob_parm_device->rho_l);
  pp.query("p_r", PeleC::h_prob_parm_device->p_r);
  pp.query("u_r", PeleC::h_prob_parm_device->u_r);
  pp.query("rho_r", PeleC::h_prob_parm_device->rho_r);
  pp.query("frac", PeleC::h_prob_parm_device->frac);
  pp.query("idir", PeleC::h_prob_parm_device->idir);
  pp.get("left_gas", PeleC::prob_parm_host->gasL);
  pp.get("right_gas", PeleC::prob_parm_host->gasR);

  for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
    PeleC::h_prob_parm_device->split[idir] = 0.5;
  }

  if (PeleC::prob_parm_host->gasL == "N2") {
    PeleC::h_prob_parm_device->left_gas_id = N2_ID;
    PeleC::h_prob_parm_device->right_gas_id = HE_ID;
  } else {
    PeleC::h_prob_parm_device->left_gas_id = HE_ID;
    PeleC::h_prob_parm_device->right_gas_id = N2_ID;
  }

  amrex::Real e_l;
  amrex::Real e_r /*, cs, cp*/;
  amrex::Real massfrac_l[NUM_SPECIES] = {0.0};
  amrex::Real massfrac_r[NUM_SPECIES] = {0.0};
  massfrac_l[PeleC::h_prob_parm_device->left_gas_id] = 1.0;
  massfrac_r[PeleC::h_prob_parm_device->right_gas_id] = 1.0;

  auto eos = pele::physics::PhysicsType::eos();
  eos.RYP2E(
    PeleC::h_prob_parm_device->rho_l, massfrac_l,
    PeleC::h_prob_parm_device->p_l, e_l);
  eos.EY2T(e_l, massfrac_l, PeleC::h_prob_parm_device->T_l);
  PeleC::h_prob_parm_device->rhoe_l = PeleC::h_prob_parm_device->rho_l * e_l;

  eos.RYP2E(
    PeleC::h_prob_parm_device->rho_r, massfrac_r,
    PeleC::h_prob_parm_device->p_r, e_r);
  eos.EY2T(e_r, massfrac_r, PeleC::h_prob_parm_device->T_r);
  PeleC::h_prob_parm_device->rhoe_r = PeleC::h_prob_parm_device->rho_r * e_r;
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
