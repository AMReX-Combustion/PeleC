// AMReX
#include "AMReX_ParmParse.H" // for ParmParse
#include "AMReX_REAL.H"      // for Real, amrex_real

// PelePhysics
#include "Fuego.H"       // for Fuego
#include "mechanism.H"   // for H2_ID, N2_ID, O2_ID
#include "PelePhysics.H" // for PhysicsType

// PeleC
#include "PeleC.H" // for PeleC, PeleC::h_prob_parm_device
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
  const amrex::Real* /*problo*/,
  const amrex::Real* /*probhi*/)
{
  // Parse params
  {
    amrex::ParmParse pp("prob");
    pp.query("p_init", PeleC::h_prob_parm_device->p_init);
    pp.query("Y_init_H2", PeleC::h_prob_parm_device->Y_init_H2);
    pp.query("Y_init_O2", PeleC::h_prob_parm_device->Y_init_O2);
    pp.query("Y_init_N2", PeleC::h_prob_parm_device->Y_init_N2);
    pp.query("T_init", PeleC::h_prob_parm_device->T_init);
  }

  // Initial values
  PeleC::h_prob_parm_device->massfrac[H2_ID] =
    PeleC::h_prob_parm_device->Y_init_H2;
  PeleC::h_prob_parm_device->massfrac[O2_ID] =
    PeleC::h_prob_parm_device->Y_init_O2;
  PeleC::h_prob_parm_device->massfrac[N2_ID] =
    PeleC::h_prob_parm_device->Y_init_N2;
  auto eos = pele::physics::PhysicsType::eos();
  eos.PYT2RE(
    PeleC::h_prob_parm_device->p_init,
    PeleC::h_prob_parm_device->massfrac.begin(),
    PeleC::h_prob_parm_device->T_init, PeleC::h_prob_parm_device->rho_init,
    PeleC::h_prob_parm_device->e_init);
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
