
#include <SprayParticles.H>
#include <AMReX_Particles.H>
#include <PeleC.H>
#include "prob.H"

bool
SprayParticleContainer::injectParticles(
  amrex::Real /*time*/,
  amrex::Real /*dt*/,
  int /*nstep*/,
  int /*lev*/,
  int /*finest_level*/,
  ProbParmHost const& /*prob_parm*/,
  ProbParmDevice const& /*prob_parm_d*/)
{
  return false;
}

void
SprayParticleContainer::InitSprayParticles(
  const bool /*init_parts*/,
  ProbParmHost const& /*prob_parm*/,
  ProbParmDevice const& /*prob_parm_d*/)
{
}
