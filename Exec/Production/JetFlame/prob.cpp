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

  pp.query("P_mean", PeleC::h_prob_parm_device->P_mean);
  pp.query("inj_start", PeleC::h_prob_parm_device->inj_start);
  pp.query("inj_dur", PeleC::h_prob_parm_device->inj_dur);
  pp.query("v_in", PeleC::h_prob_parm_device->v_in);
  pp.query("D", PeleC::h_prob_parm_device->D);
  pp.query("Z", PeleC::h_prob_parm_device->Z);
  pp.query("T_fu", PeleC::h_prob_parm_device->T_fu);
  pp.query("T_ox", PeleC::h_prob_parm_device->T_ox);
  pp.query("tau", PeleC::h_prob_parm_device->tau);
  pp.query("turbulence", PeleC::h_prob_parm_device->turbulence);

  std::string fu_spec = "";
  std::string fu_ox_spec = "";
  pp.query("fu_spec", fu_spec);
  pp.query("fu_ox_spec", fu_ox_spec);
  amrex::Real Y_O2_ox = {0.};
  amrex::Real Y_fu_ox = {0.};
  pp.query("Y_O2_ox", Y_O2_ox);
  pp.query("Y_fu_ox", Y_fu_ox);

  amrex::Vector<std::string> sname;
  pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(sname);
  amrex::Real Y_pure_fuel[NUM_SPECIES] = {0.0};
  int fu_indx = -1;
  int o2_indx = -1;
  int n2_indx = -1;
  int fu_ox_indx = -1;
  for (int n = 0; n < sname.size(); n++) {
    if (sname[n] == fu_spec)
      fu_indx = n;
    if (sname[n] == "O2")
      o2_indx = n;
    if (sname[n] == "N2")
      n2_indx = n;
    if (sname[n] == fu_ox_spec)
      fu_ox_indx = n;
  }

  if (fu_indx < 0)
    amrex::Abort("Fuel species not found.");

  Y_pure_fuel[fu_indx] = 1.0;

  PeleC::h_prob_parm_device->Y_ox[o2_indx] = Y_O2_ox;
  PeleC::h_prob_parm_device->Y_ox[fu_ox_indx] = Y_fu_ox;
  PeleC::h_prob_parm_device->Y_ox[n2_indx] = 1.0 - Y_fu_ox - Y_O2_ox;

  for (int n = 0; n < NUM_SPECIES; n++) {
    PeleC::h_prob_parm_device->Y_fuel[n] =
      PeleC::h_prob_parm_device->Z * Y_pure_fuel[n] +
      (1. - PeleC::h_prob_parm_device->Z) * PeleC::h_prob_parm_device->Y_ox[n];
  }

  CKUBMS(
    &PeleC::h_prob_parm_device->T_fu, PeleC::h_prob_parm_device->Y_fuel,
    &PeleC::h_prob_parm_device->U_fuel);
  CKUBMS(
    &PeleC::h_prob_parm_device->T_ox, PeleC::h_prob_parm_device->Y_ox,
    &PeleC::h_prob_parm_device->U_ox);

  PeleC::h_prob_parm_device->center_xy[0] = 0.5 * (probhi[0] + problo[0]);
  PeleC::h_prob_parm_device->center_xy[1] = 0.5 * (probhi[1] + problo[1]);
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
