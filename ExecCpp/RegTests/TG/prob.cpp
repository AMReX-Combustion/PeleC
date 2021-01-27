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
  const amrex_real* problo,
  const amrex_real* probhi)
{
  // Parse params
  amrex::ParmParse pp("prob");
  pp.query("reynolds", PeleC::prob_parm_device->reynolds);
  pp.query("mach", PeleC::prob_parm_device->mach);
  pp.query("prandtl", PeleC::prob_parm_device->prandtl);
  pp.query("convecting", PeleC::prob_parm_device->convecting);
  pp.query("omega_x", PeleC::prob_parm_device->omega_x);
  pp.query("omega_y", PeleC::prob_parm_device->omega_y);
  pp.query("omega_z", PeleC::prob_parm_device->omega_z);

  // Define the length scale
  PeleC::prob_parm_device->L = 1.0 / PI;
  PeleC::prob_parm_device->L_x = probhi[0] - problo[0];
  PeleC::prob_parm_device->L_y = probhi[1] - problo[1];
  PeleC::prob_parm_device->L_z = probhi[2] - problo[2];

  // Initial density, velocity, and material properties
  amrex::Real eint;
  amrex::Real cs;
  amrex::Real cp;
  amrex::Real massfrac[NUM_SPECIES] = {1.0};
  EOS::PYT2RE(
    PeleC::prob_parm_device->p0, massfrac, PeleC::prob_parm_device->T0,
    PeleC::prob_parm_device->rho0, eint);
  EOS::RTY2Cs(
    PeleC::prob_parm_device->rho0, PeleC::prob_parm_device->T0, massfrac, cs);
  EOS::TY2Cp(PeleC::prob_parm_device->T0, massfrac, cp);

  PeleC::prob_parm_device->v0 = PeleC::prob_parm_device->mach * cs;
  transport_params::const_bulk_viscosity = 0.0;
  transport_params::const_diffusivity = 0.0;
  transport_params::const_viscosity =
    PeleC::prob_parm_device->rho0 * PeleC::prob_parm_device->v0 *
    PeleC::prob_parm_device->L / PeleC::prob_parm_device->reynolds;
  transport_params::const_conductivity =
    transport_params::const_viscosity * cp / PeleC::prob_parm_device->prandtl;

  // Output IC
  std::ofstream ofs("ic.txt", std::ofstream::out);
  amrex::Print(ofs) << "L, rho0, v0, p0, T0, gamma, mu, k, c_s0, Reynolds, "
                       "Mach, Prandtl, omega_x, omega_y, omega_z"
                    << std::endl;
  amrex::Print(ofs).SetPrecision(17)
    << PeleC::prob_parm_device->L << "," << PeleC::prob_parm_device->rho0 << ","
    << PeleC::prob_parm_device->v0 << "," << PeleC::prob_parm_device->p0 << ","
    << PeleC::prob_parm_device->T0 << "," << EOS::gamma << ","
    << transport_params::const_viscosity << ","
    << transport_params::const_conductivity << "," << cs << ","
    << PeleC::prob_parm_device->reynolds << "," << PeleC::prob_parm_device->mach
    << "," << PeleC::prob_parm_device->prandtl << ","
    << PeleC::prob_parm_device->omega_x << ","
    << PeleC::prob_parm_device->omega_y << ","
    << PeleC::prob_parm_device->omega_z << std::endl;
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
