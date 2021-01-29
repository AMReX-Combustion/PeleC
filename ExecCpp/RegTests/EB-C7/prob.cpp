#include "prob.H"

namespace ProbParm {
AMREX_GPU_DEVICE_MANAGED amrex::Real pl = 1.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real rhol = 1.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real Tl = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real eintl = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real pr = 0.1;
AMREX_GPU_DEVICE_MANAGED amrex::Real rhor = 0.125;
AMREX_GPU_DEVICE_MANAGED amrex::Real eintr = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real Tr = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real angle = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real L = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {
  1.0};
} // namespace ProbParm

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
  pp.query("pl", PeleC::prob_parm_device->pl);
  pp.query("rhol", PeleC::prob_parm_device->rhol);
  pp.query("pr", PeleC::prob_parm_device->pr);
  pp.query("rhor", PeleC::prob_parm_device->rhor);
  pp.query("angle", PeleC::prob_parm_device->angle);

  PeleC::prob_parm_device->L = (probhi[0] - problo[0]);

  PeleC::prob_parm_device->massfrac[0] = 1.0;

  EOS::RYP2T(
    PeleC::prob_parm_device->rhol, PeleC::prob_parm_device->massfrac.begin(),
    PeleC::prob_parm_device->pl, PeleC::prob_parm_device->Tl);
  EOS::RYP2E(
    PeleC::prob_parm_device->rhol, PeleC::prob_parm_device->massfrac.begin(),
    PeleC::prob_parm_device->pl, PeleC::prob_parm_device->eintl);

  EOS::RYP2T(
    PeleC::prob_parm_device->rhor, PeleC::prob_parm_device->massfrac.begin(),
    PeleC::prob_parm_device->pr, PeleC::prob_parm_device->Tr);
  EOS::RYP2E(
    PeleC::prob_parm_device->rhor, PeleC::prob_parm_device->massfrac.begin(),
    PeleC::prob_parm_device->pr, PeleC::prob_parm_device->eintr);

  // Output IC
  std::ofstream ofs("ic.txt", std::ofstream::out);
  amrex::Print(ofs) << "L, rhol, pl, Tl, rhor, pr, Tr, gamma, angle"
                    << std::endl;
  amrex::Print(ofs).SetPrecision(17)
    << PeleC::prob_parm_device->L << "," << PeleC::prob_parm_device->rhol << ","
    << PeleC::prob_parm_device->pl << "," << PeleC::prob_parm_device->Tl << ","
    << PeleC::prob_parm_device->rhor << "," << PeleC::prob_parm_device->pr
    << "," << PeleC::prob_parm_device->Tr << "," << EOS::gamma << ","
    << PeleC::prob_parm_device->angle << std::endl;
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
