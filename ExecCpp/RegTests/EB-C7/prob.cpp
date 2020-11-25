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
  pp.query("pl", ProbParm::pl);
  pp.query("rhol", ProbParm::rhol);
  pp.query("pr", ProbParm::pr);
  pp.query("rhor", ProbParm::rhor);
  pp.query("angle", ProbParm::angle);

  ProbParm::L = (probhi[0] - problo[0]);

  ProbParm::massfrac[0] = 1.0;

  EOS::RYP2T(
    ProbParm::rhol, ProbParm::massfrac.begin(), ProbParm::pl, ProbParm::Tl);
  EOS::RYP2E(
    ProbParm::rhol, ProbParm::massfrac.begin(), ProbParm::pl, ProbParm::eintl);

  EOS::RYP2T(
    ProbParm::rhor, ProbParm::massfrac.begin(), ProbParm::pr, ProbParm::Tr);
  EOS::RYP2E(
    ProbParm::rhor, ProbParm::massfrac.begin(), ProbParm::pr, ProbParm::eintr);

  // Output IC
  std::ofstream ofs("ic.txt", std::ofstream::out);
  amrex::Print(ofs) << "L, rhol, pl, Tl, rhor, pr, Tr, gamma, angle"
                    << std::endl;
  amrex::Print(ofs).SetPrecision(17)
    << ProbParm::L << "," << ProbParm::rhol << "," << ProbParm::pl << ","
    << ProbParm::Tl << "," << ProbParm::rhor << "," << ProbParm::pr << ","
    << ProbParm::Tr << "," << EOS::gamma << "," << ProbParm::angle << std::endl;
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
