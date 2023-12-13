
#include "SprayInjection.H"
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
  const bool init_parts,
  ProbParmHost const& prob_parm,
  ProbParmDevice const& /*prob_parm_d*/)
{
  if (!init_parts) {
    return;
  }
  const int level = 0;
  amrex::ParmParse pp("prob");
  amrex::IntVect partNum(AMREX_D_DECL(100, 100, 100));
  pp.query("num_particles", partNum);
  amrex::RealVect partVel = amrex::RealVect::TheZeroVector();
  amrex::Real partDia = prob_parm.partDia;
  amrex::Real partTemp = prob_parm.partTemp;
  int numRedist = 1;
  std::array<amrex::Real, SPRAY_FUEL_NUM> partY = {0.0};
  partY[0] = 1.;
  uniformSprayInit(
    partNum, partVel, partDia, partTemp, partY.begin(), level, numRedist);
}
