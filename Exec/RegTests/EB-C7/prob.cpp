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
  pp.query("pl", PeleC::h_prob_parm_device->pl);
  pp.query("rhol", PeleC::h_prob_parm_device->rhol);
  pp.query("pr", PeleC::h_prob_parm_device->pr);
  pp.query("rhor", PeleC::h_prob_parm_device->rhor);
  pp.query("angle", PeleC::h_prob_parm_device->angle);

  PeleC::h_prob_parm_device->L = (probhi[0] - problo[0]);

  PeleC::h_prob_parm_device->massfrac[0] = 1.0;

  auto eos = pele::physics::PhysicsType::eos();
  eos.RYP2T(
    PeleC::h_prob_parm_device->rhol,
    PeleC::h_prob_parm_device->massfrac.begin(), PeleC::h_prob_parm_device->pl,
    PeleC::h_prob_parm_device->Tl);
  eos.RYP2E(
    PeleC::h_prob_parm_device->rhol,
    PeleC::h_prob_parm_device->massfrac.begin(), PeleC::h_prob_parm_device->pl,
    PeleC::h_prob_parm_device->eintl);

  eos.RYP2T(
    PeleC::h_prob_parm_device->rhor,
    PeleC::h_prob_parm_device->massfrac.begin(), PeleC::h_prob_parm_device->pr,
    PeleC::h_prob_parm_device->Tr);
  eos.RYP2E(
    PeleC::h_prob_parm_device->rhor,
    PeleC::h_prob_parm_device->massfrac.begin(), PeleC::h_prob_parm_device->pr,
    PeleC::h_prob_parm_device->eintr);

  // Output IC
  std::ofstream ofs("ic.txt", std::ofstream::out);
  amrex::Print(ofs) << "L, rhol, pl, Tl, rhor, pr, Tr, gamma, angle"
                    << std::endl;
  amrex::Print(ofs).SetPrecision(17)
    << PeleC::h_prob_parm_device->L << "," << PeleC::h_prob_parm_device->rhol
    << "," << PeleC::h_prob_parm_device->pl << ","
    << PeleC::h_prob_parm_device->Tl << "," << PeleC::h_prob_parm_device->rhor
    << "," << PeleC::h_prob_parm_device->pr << ","
    << PeleC::h_prob_parm_device->Tr << "," << eos.gamma << ","
    << PeleC::h_prob_parm_device->angle << std::endl;
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
