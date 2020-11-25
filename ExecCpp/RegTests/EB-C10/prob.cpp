#include "prob.H"

namespace ProbParm {
AMREX_GPU_DEVICE_MANAGED amrex::Real p = 1013250.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real dpdx = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real T = 300.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real rho = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real eint = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real umax = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real uavg = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real Re = 100.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real Ma = 0.1;
AMREX_GPU_DEVICE_MANAGED amrex::Real Pr = 0.7;
AMREX_GPU_DEVICE_MANAGED amrex::Real radius = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real G = 0.0;
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
  pp.query("p", ProbParm::p);
  pp.query("T", ProbParm::T);
  pp.query("Re", ProbParm::Re);
  pp.query("Ma", ProbParm::Ma);
  pp.query("Pr", ProbParm::Pr);

  amrex::ParmParse ppeb("eb2");
  ppeb.query("cylinder_radius", ProbParm::radius);

  amrex::Real L = (probhi[0] - problo[0]);

  amrex::Real cp = 0.0, cs = 0.0;
  ProbParm::massfrac[0] = 1.0;

  EOS::PYT2RE(
    ProbParm::p, ProbParm::massfrac.begin(), ProbParm::T, ProbParm::rho,
    ProbParm::eint);
  EOS::RTY2Cs(ProbParm::rho, ProbParm::T, ProbParm::massfrac.begin(), cs);
  EOS::TY2Cp(ProbParm::T, ProbParm::massfrac.begin(), cp);

  ProbParm::umax = ProbParm::Ma * cs;
  ProbParm::uavg = 0.5 * ProbParm::umax;

  transport_params::const_bulk_viscosity = 0.0;
  transport_params::const_diffusivity = 0.0;
  transport_params::const_viscosity =
    ProbParm::rho * ProbParm::umax * L / ProbParm::Re;
  transport_params::const_conductivity =
    transport_params::const_viscosity * cp / ProbParm::Pr;

  ProbParm::G = ProbParm::umax * 4 * transport_params::const_viscosity /
                (ProbParm::radius * ProbParm::radius);
  ProbParm::dpdx = -ProbParm::G;

  // Output IC
  std::ofstream ofs("ic.txt", std::ofstream::out);
  amrex::Print(ofs)
    << "L, rho, umax, p, T, gamma, mu, k, Re, Ma, Pr, dpdx, G, radius"
    << std::endl;
  amrex::Print(ofs).SetPrecision(17)
    << L << "," << ProbParm::rho << "," << ProbParm::umax << "," << ProbParm::p
    << "," << ProbParm::T << "," << EOS::gamma << ","
    << transport_params::const_viscosity << ","
    << transport_params::const_conductivity << "," << ProbParm::Re << ","
    << ProbParm::Ma << "," << ProbParm::Pr << "," << ProbParm::dpdx << ","
    << ProbParm::G << "," << ProbParm::radius << std::endl;
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
