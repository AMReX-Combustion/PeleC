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
  pp.query("init_p", PeleC::h_prob_parm_device->p0);
  pp.query("init_T", PeleC::h_prob_parm_device->T0);
  pp.query("init_N2", PeleC::h_prob_parm_device->Y_N2);
  pp.query("init_O2", PeleC::h_prob_parm_device->Y_O2);
  pp.query("init_fuel", PeleC::h_prob_parm_device->Y_Fuel);
  std::string fuel_name = "C2H4";
  pp.query("fuel_name", fuel_name);
  amrex::Vector<std::string> spec_names;
  pele::physics::eos::speciesNames<pele::physics::EosType>(spec_names);
  for (int sp = 0; sp < NUM_SPECIES; ++sp) {
    std::string spec_name = spec_names[sp];
    if (spec_name == fuel_name)
      PeleC::h_prob_parm_device->fuelIndx = sp;
  }
  if (PeleC::h_prob_parm_device->fuelIndx < 0)
    amrex::Abort("Fuel not found in chemistry mechanism");

  // Initial density, velocity, and material properties
  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real eint, cs, cp;
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  massfrac[N2_ID] = PeleC::h_prob_parm_device->Y_N2;
  massfrac[O2_ID] = PeleC::h_prob_parm_device->Y_O2;
  massfrac[PeleC::h_prob_parm_device->fuelIndx] =
    PeleC::h_prob_parm_device->Y_Fuel;
  eos.PYT2RE(
    PeleC::h_prob_parm_device->p0, massfrac, PeleC::h_prob_parm_device->T0,
    PeleC::h_prob_parm_device->rho0, eint);
  eos.RTY2Cs(
    PeleC::h_prob_parm_device->rho0, PeleC::h_prob_parm_device->T0, massfrac,
    cs);
  eos.TY2Cp(PeleC::h_prob_parm_device->T0, massfrac, cp);
  amrex::Real moments[NUM_SOOT_MOMENTS + 1];
  SootData* const sd = PeleC::soot_model.getSootData();
  sd->initialSmallMomVals(moments);
  for (int n = 0; n < NUM_SOOT_MOMENTS + 1; ++n) {
    PeleC::h_prob_parm_device->soot_vals[n] = moments[n];
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
