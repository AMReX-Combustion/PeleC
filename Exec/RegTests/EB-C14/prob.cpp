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
  {
    amrex::ParmParse pp("prob");
    pp.query("p0", PeleC::h_prob_parm_device->p0);
    pp.query("p1", PeleC::h_prob_parm_device->p1);
    pp.query("rho0", PeleC::h_prob_parm_device->rho0);
    pp.query("rho1", PeleC::h_prob_parm_device->rho1);
    pp.query("x1", PeleC::h_prob_parm_device->x1);
  }

  PeleC::h_prob_parm_device->massfrac[0] = 1.0;

  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real cp = 0.0;
  eos.RYP2T(
    PeleC::h_prob_parm_device->rho0,
    PeleC::h_prob_parm_device->massfrac.begin(), PeleC::h_prob_parm_device->p0,
    PeleC::h_prob_parm_device->T0);
  eos.RYP2E(
    PeleC::h_prob_parm_device->rho0,
    PeleC::h_prob_parm_device->massfrac.begin(), PeleC::h_prob_parm_device->p0,
    PeleC::h_prob_parm_device->e0);

  eos.RYP2T(
    PeleC::h_prob_parm_device->rho1,
    PeleC::h_prob_parm_device->massfrac.begin(), PeleC::h_prob_parm_device->p1,
    PeleC::h_prob_parm_device->T1);
  eos.RYP2E(
    PeleC::h_prob_parm_device->rho1,
    PeleC::h_prob_parm_device->massfrac.begin(), PeleC::h_prob_parm_device->p1,
    PeleC::h_prob_parm_device->e1);

  eos.TY2Cp(
    PeleC::h_prob_parm_device->T1, PeleC::h_prob_parm_device->massfrac.begin(),
    cp);

  pele::physics::transport::TransParm trans_parm;

  // Default
  trans_parm.const_viscosity = 1.21e-2;
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

  amrex::Real Pr = 0.7;
  trans_parm.const_conductivity = trans_parm.const_viscosity * cp / Pr;

#ifdef AMREX_USE_GPU
  amrex::Gpu::htod_memcpy(
    pele::physics::transport::trans_parm_g, &trans_parm, sizeof(trans_parm));
#else
  std::memcpy(
    pele::physics::transport::trans_parm_g, &trans_parm, sizeof(trans_parm));
#endif
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
