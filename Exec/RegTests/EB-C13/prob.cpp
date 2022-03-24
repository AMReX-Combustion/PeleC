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
    pp.query("mach", PeleC::h_prob_parm_device->mach);
    pp.query("rho", PeleC::h_prob_parm_device->rho);
    pp.query("p", PeleC::h_prob_parm_device->p);
  }
  {
    amrex::ParmParse pp("eb2");
    pp.query("r_inner", PeleC::h_prob_parm_device->r_inner);
    pp.query("r_outer", PeleC::h_prob_parm_device->r_outer);
  }
  auto eos = pele::physics::PhysicsType::eos();
  eos.RPY2Cs(
    PeleC::h_prob_parm_device->rho, PeleC::h_prob_parm_device->p,
    PeleC::h_prob_parm_device->massfrac.data(), PeleC::h_prob_parm_device->c_s);
  eos.RYP2T(
    PeleC::h_prob_parm_device->rho, PeleC::h_prob_parm_device->massfrac.data(),
    PeleC::h_prob_parm_device->p, PeleC::h_prob_parm_device->T);
  eos.RTY2E(
    PeleC::h_prob_parm_device->rho, PeleC::h_prob_parm_device->T,
    PeleC::h_prob_parm_device->massfrac.data(),
    PeleC::h_prob_parm_device->eint);

  // Output IC
  std::ofstream ofs("ic.txt", std::ofstream::out);
  amrex::Print(ofs) << "mach, rho, p, c_s, T" << std::endl;
  amrex::Print(ofs).SetPrecision(17)
    << PeleC::h_prob_parm_device->mach << "," << PeleC::h_prob_parm_device->rho
    << "," << PeleC::h_prob_parm_device->p << ","
    << PeleC::h_prob_parm_device->c_s << "," << PeleC::h_prob_parm_device->T
    << std::endl;
  ofs.close();
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
