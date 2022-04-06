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
    pp.query("p", PeleC::h_prob_parm_device->p);
    pp.query("rho", PeleC::h_prob_parm_device->rho);
    pp.query("vx_in", PeleC::h_prob_parm_device->vx_in);
    pp.query("vy_in", PeleC::h_prob_parm_device->vy_in);
    pp.query("Re_L", PeleC::h_prob_parm_device->Re_L);
    pp.query("Pr", PeleC::h_prob_parm_device->Pr);
  }

  amrex::Real L = (probhi[0] - problo[0]) * 0.2;

  amrex::Real cp = 0.0;
  PeleC::h_prob_parm_device->massfrac[0] = 1.0;
  auto eos = pele::physics::PhysicsType::eos();
  eos.RYP2E(
    PeleC::h_prob_parm_device->rho, PeleC::h_prob_parm_device->massfrac.begin(),
    PeleC::h_prob_parm_device->p, PeleC::h_prob_parm_device->eint);
  eos.EY2T(
    PeleC::h_prob_parm_device->eint,
    PeleC::h_prob_parm_device->massfrac.begin(), PeleC::h_prob_parm_device->T);
  eos.TY2Cp(
    PeleC::h_prob_parm_device->T, PeleC::h_prob_parm_device->massfrac.begin(),
    cp);

  auto& trans_parm = PeleC::trans_parms.host_trans_parm();
  trans_parm.const_bulk_viscosity = 0.0;
  trans_parm.const_diffusivity = 0.0;
  trans_parm.const_viscosity = PeleC::h_prob_parm_device->rho *
                               PeleC::h_prob_parm_device->vx_in * L /
                               PeleC::h_prob_parm_device->Re_L;
  trans_parm.const_conductivity =
    trans_parm.const_viscosity * cp / PeleC::h_prob_parm_device->Pr;
  PeleC::trans_parms.sync_to_device();
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
