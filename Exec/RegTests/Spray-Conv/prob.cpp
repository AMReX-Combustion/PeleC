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
  pp.query("mach", PeleC::h_prob_parm_device->mach);
  pp.query("ref_p", PeleC::h_prob_parm_device->p0);
  pp.query("ref_T", PeleC::h_prob_parm_device->T0);
  pp.query("init_O2", PeleC::h_prob_parm_device->Y_O2);
  pp.query("init_N2", PeleC::h_prob_parm_device->Y_N2);

  // Initial density, velocity, and material properties
  amrex::Real cs;
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  massfrac[O2_ID] = PeleC::h_prob_parm_device->Y_O2;
  massfrac[N2_ID] = PeleC::h_prob_parm_device->Y_N2;
  amrex::Real wbar, gamma0;
  auto eos = pele::physics::PhysicsType::eos();
  eos.TY2G(PeleC::h_prob_parm_device->T0, massfrac, gamma0);
  eos.Y2WBAR(massfrac, wbar);
  // Compute the speed of sound
  cs = std::sqrt(
    PeleC::h_prob_parm_device->T0 * gamma0 * pele::physics::Constants::RU /
    wbar);
  PeleC::h_prob_parm_device->v0 = PeleC::h_prob_parm_device->mach * cs;
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
