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
  {
    amrex::ParmParse pp("prob");
    pp.query("p", ProbParm::p);
    pp.query("T", ProbParm::T);
    pp.query("Re", ProbParm::Re);
    pp.query("Ma", ProbParm::Ma);
    pp.query("Pr", ProbParm::Pr);
  }

  {
    amrex::ParmParse pp("eb2");
    pp.query("cylinder_radius", ProbParm::radius);
  }

  amrex::Real L = (probhi[0] - problo[0]);

  amrex::Real cp = 0.0;
  amrex::Real cs = 0.0;
  ProbParm::massfrac[0] = 1.0;

  EOS::PYT2RE(
    ProbParm::p, ProbParm::massfrac.begin(), ProbParm::T, ProbParm::rho,
    ProbParm::eint);
  EOS::RTY2Cs(ProbParm::rho, ProbParm::T, ProbParm::massfrac.begin(), cs);
  EOS::TY2Cp(ProbParm::T, ProbParm::massfrac.begin(), cp);

  ProbParm::umax = ProbParm::Ma * cs;
  ProbParm::uavg = 0.5 * ProbParm::umax;

  TransParm trans_parm;

  // Default
  trans_parm.const_viscosity = 0.0;
  trans_parm.const_bulk_viscosity = 0.0;
  trans_parm.const_conductivity = 0.0;
  trans_parm.const_diffusivity = 0.0;

  // User-specified
  {
    amrex::ParmParse pp("transport");
    pp.query("const_viscosity", trans_parm.const_viscosity);
    pp.query("const_bulk_viscosity", trans_parm.const_bulk_viscosity);
    pp.query("const_conductivity", trans_parm.const_conductivity);
    pp.query("const_diffusivity", trans_parm.const_diffusivity);
  }

  trans_parm.const_bulk_viscosity = 0.0;
  trans_parm.const_diffusivity = 0.0;
  trans_parm.const_viscosity =
    ProbParm::rho * ProbParm::umax * L / ProbParm::Re;
  trans_parm.const_conductivity =
    trans_parm.const_viscosity * cp / ProbParm::Pr;

#ifdef AMREX_USE_GPU
  amrex::Gpu::htod_memcpy(trans_parm_g, &trans_parm, sizeof(trans_parm));
#else
  std::memcpy(trans_parm_g, &trans_parm, sizeof(trans_parm));
#endif

  ProbParm::G = ProbParm::umax * 4 * trans_parm.const_viscosity /
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
    << trans_parm.const_viscosity << "," << trans_parm.const_conductivity << ","
    << ProbParm::Re << "," << ProbParm::Ma << "," << ProbParm::Pr << ","
    << ProbParm::dpdx << "," << ProbParm::G << "," << ProbParm::radius
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
