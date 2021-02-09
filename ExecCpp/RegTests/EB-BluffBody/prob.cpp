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
    pp.query("p", PeleC::prob_parm_device->p);
    pp.query("rho", PeleC::prob_parm_device->rho);
    pp.query("vx_in", PeleC::prob_parm_device->vx_in);
    pp.query("vy_in", PeleC::prob_parm_device->vy_in);
    pp.query("Re_L", PeleC::prob_parm_device->Re_L);
    pp.query("Pr", PeleC::prob_parm_device->Pr);
  }

  amrex::Real L = (probhi[0] - problo[0]) * 0.2;

  amrex::Real cp = 0.0;
  PeleC::prob_parm_device->massfrac[0] = 1.0;
  EOS::RYP2E(
    PeleC::prob_parm_device->rho, PeleC::prob_parm_device->massfrac.begin(),
    PeleC::prob_parm_device->p, PeleC::prob_parm_device->eint);
  EOS::EY2T(
    PeleC::prob_parm_device->eint, PeleC::prob_parm_device->massfrac.begin(),
    PeleC::prob_parm_device->T);
  EOS::TY2Cp(
    PeleC::prob_parm_device->T, PeleC::prob_parm_device->massfrac.begin(), cp);

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
  trans_parm.const_viscosity = PeleC::prob_parm_device->rho *
                               PeleC::prob_parm_device->vx_in * L /
                               PeleC::prob_parm_device->Re_L;
  trans_parm.const_conductivity =
    trans_parm.const_viscosity * cp / PeleC::prob_parm_device->Pr;

#ifdef AMREX_USE_GPU
  amrex::Gpu::htod_memcpy(trans_parm_g, &trans_parm, sizeof(trans_parm));
#else
  std::memcpy(trans_parm_g, &trans_parm, sizeof(trans_parm));
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
