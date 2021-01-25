#include "prob.H"

void
pc_prob_close()
{
}

extern "C" {
void
amrex_probinit(
  const int* init,
  const int* name,
  const int* namelen,
  const amrex_real* problo,
  const amrex_real* probhi)
{
  // Parse params
  amrex::ParmParse pp("prob");
  pp.query("p_l", PeleC::h_prob_parm_device->p_l);
  pp.query("T_l", PeleC::h_prob_parm_device->T_l);
  pp.query("p_r", PeleC::h_prob_parm_device->p_r);
  pp.query("T_r", PeleC::h_prob_parm_device->T_r);
  pp.query("U_r", PeleC::h_prob_parm_device->U_r);
  pp.get("left_gas", PeleC::prob_parm_host->gasL);
  pp.get("right_gas", PeleC::prob_parm_host->gasR);

  if (PeleC::prob_parm_host->gasL == "N2") {
    PeleC::h_prob_parm_device->left_gas_id = N2_ID;
  } else {
    PeleC::h_prob_parm_device->left_gas_id = O2_ID;
  }
  if (PeleC::prob_parm_host->gasR == "N2") {
    PeleC::h_prob_parm_device->right_gas_id = N2_ID;
  } else {
    PeleC::h_prob_parm_device->right_gas_id = O2_ID;
  }
  amrex::Real e_l, e_r, cs, cp;
  amrex::Real massfrac_l[NUM_SPECIES] = {0.0};
  amrex::Real massfrac_r[NUM_SPECIES] = {0.0};
  massfrac_l[PeleC::h_prob_parm_device->left_gas_id] = 1.0;
  massfrac_r[PeleC::h_prob_parm_device->right_gas_id] = 1.0;
  std::cout << "left " << PeleC::h_prob_parm_device->left_gas_id << std::endl;
  std::cout << "right " << PeleC::h_prob_parm_device->right_gas_id << std::endl;

  auto eos = pele::physics::PhysicsType::eos();
  eos.PYT2RE(
    PeleC::h_prob_parm_device->p_l, massfrac_l, PeleC::h_prob_parm_device->T_l,
    PeleC::h_prob_parm_device->rho_l, e_l);
  eos.PYT2RE(
    PeleC::h_prob_parm_device->p_r, massfrac_r, PeleC::h_prob_parm_device->T_r,
    PeleC::h_prob_parm_device->rho_r, e_r);
  PeleC::h_prob_parm_device->rhoe_l = PeleC::h_prob_parm_device->rho_l * e_l;
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
