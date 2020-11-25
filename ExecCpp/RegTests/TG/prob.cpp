#include "prob.H"

namespace ProbParm {
AMREX_GPU_DEVICE_MANAGED amrex::Real reynolds = 1600.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real mach = 0.1;
AMREX_GPU_DEVICE_MANAGED amrex::Real prandtl = 0.71;
AMREX_GPU_DEVICE_MANAGED bool convecting = false;
AMREX_GPU_DEVICE_MANAGED amrex::Real omega_x = 1.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real omega_y = 1.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real omega_z = 1.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real L_x = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real L_y = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real L_z = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real L = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real p0 = 1.013e6; // [erg cm^-3]
AMREX_GPU_DEVICE_MANAGED amrex::Real T0 = 300.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real rho0 = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real v0 = 0.0;
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
  pp.query("reynolds", ProbParm::reynolds);
  pp.query("mach", ProbParm::mach);
  pp.query("prandtl", ProbParm::prandtl);
  pp.query("convecting", ProbParm::convecting);
  pp.query("omega_x", ProbParm::omega_x);
  pp.query("omega_y", ProbParm::omega_y);
  pp.query("omega_z", ProbParm::omega_z);

  // Define the length scale
  ProbParm::L = 1.0 / PI;
  ProbParm::L_x = probhi[0] - problo[0];
  ProbParm::L_y = probhi[1] - problo[1];
  ProbParm::L_z = probhi[2] - problo[2];

  // Initial density, velocity, and material properties
  amrex::Real eint, cs, cp;
  amrex::Real massfrac[NUM_SPECIES] = {1.0};
  EOS::PYT2RE(ProbParm::p0, massfrac, ProbParm::T0, ProbParm::rho0, eint);
  EOS::RTY2Cs(ProbParm::rho0, ProbParm::T0, massfrac, cs);
  EOS::TY2Cp(ProbParm::T0, massfrac, cp);

  ProbParm::v0 = ProbParm::mach * cs;
  transport_params::const_bulk_viscosity = 0.0;
  transport_params::const_diffusivity = 0.0;
  transport_params::const_viscosity =
    ProbParm::rho0 * ProbParm::v0 * ProbParm::L / ProbParm::reynolds;
  transport_params::const_conductivity =
    transport_params::const_viscosity * cp / ProbParm::prandtl;

  // Output IC
  std::ofstream ofs("ic.txt", std::ofstream::out);
  amrex::Print(ofs) << "L, rho0, v0, p0, T0, gamma, mu, k, c_s0, Reynolds, "
                       "Mach, Prandtl, omega_x, omega_y, omega_z"
                    << std::endl;
  amrex::Print(ofs).SetPrecision(17)
    << ProbParm::L << "," << ProbParm::rho0 << "," << ProbParm::v0 << ","
    << ProbParm::p0 << "," << ProbParm::T0 << "," << EOS::gamma << ","
    << transport_params::const_viscosity << ","
    << transport_params::const_conductivity << "," << cs << ","
    << ProbParm::reynolds << "," << ProbParm::mach << "," << ProbParm::prandtl
    << "," << ProbParm::omega_x << "," << ProbParm::omega_y << ","
    << ProbParm::omega_z << std::endl;
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
