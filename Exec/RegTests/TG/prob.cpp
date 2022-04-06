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
    pp.query("reynolds", PeleC::h_prob_parm_device->reynolds);
    pp.query("mach", PeleC::h_prob_parm_device->mach);
    pp.query("prandtl", PeleC::h_prob_parm_device->prandtl);
    pp.query("convecting", PeleC::h_prob_parm_device->convecting);
    pp.query("omega_x", PeleC::h_prob_parm_device->omega_x);
    pp.query("omega_y", PeleC::h_prob_parm_device->omega_y);
    pp.query("omega_z", PeleC::h_prob_parm_device->omega_z);
  }

  // Define the length scale
  PeleC::h_prob_parm_device->L = 1.0 / constants::PI();
  PeleC::h_prob_parm_device->L_x = probhi[0] - problo[0];
  PeleC::h_prob_parm_device->L_y = probhi[1] - problo[1];
  PeleC::h_prob_parm_device->L_z = probhi[2] - problo[2];

  // Initial density, velocity, and material properties
  amrex::Real eint;
  amrex::Real cs;
  amrex::Real cp;
  amrex::Real massfrac[NUM_SPECIES] = {1.0};
  auto eos = pele::physics::PhysicsType::eos();
  eos.PYT2RE(
    PeleC::h_prob_parm_device->p0, massfrac, PeleC::h_prob_parm_device->T0,
    PeleC::h_prob_parm_device->rho0, eint);
  eos.RTY2Cs(
    PeleC::h_prob_parm_device->rho0, PeleC::h_prob_parm_device->T0, massfrac,
    cs);
  eos.TY2Cp(PeleC::h_prob_parm_device->T0, massfrac, cp);
  PeleC::h_prob_parm_device->v0 = PeleC::h_prob_parm_device->mach * cs;

  auto& trans_parm = PeleC::trans_parms.host_trans_parm();
  trans_parm.const_bulk_viscosity = 0.0;
  trans_parm.const_diffusivity = 0.0;
  trans_parm.const_viscosity =
    PeleC::h_prob_parm_device->rho0 * PeleC::h_prob_parm_device->v0 *
    PeleC::h_prob_parm_device->L / PeleC::h_prob_parm_device->reynolds;
  trans_parm.const_conductivity =
    trans_parm.const_viscosity * cp / PeleC::h_prob_parm_device->prandtl;
  PeleC::trans_parms.sync_to_device();

  // Output IC
  std::ofstream ofs("ic.txt", std::ofstream::out);
  amrex::Print(ofs) << "L, rho0, v0, p0, T0, gamma, mu, k, c_s0, Reynolds, "
                       "Mach, Prandtl, omega_x, omega_y, omega_z"
                    << std::endl;
  amrex::Print(ofs).SetPrecision(17)
    << PeleC::h_prob_parm_device->L << "," << PeleC::h_prob_parm_device->rho0
    << "," << PeleC::h_prob_parm_device->v0 << ","
    << PeleC::h_prob_parm_device->p0 << "," << PeleC::h_prob_parm_device->T0
    << "," << eos.gamma << "," << trans_parm.const_viscosity << ","
    << trans_parm.const_conductivity << "," << cs << ","
    << PeleC::h_prob_parm_device->reynolds << ","
    << PeleC::h_prob_parm_device->mach << ","
    << PeleC::h_prob_parm_device->prandtl << ","
    << PeleC::h_prob_parm_device->omega_x << ","
    << PeleC::h_prob_parm_device->omega_y << ","
    << PeleC::h_prob_parm_device->omega_z << std::endl;
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
