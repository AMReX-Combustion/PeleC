
#include "SprayInjection.H"
#include "PeleC.H"
#include "prob.H"

bool
SprayParticleContainer::injectParticles( // NOLINT
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
  ProbParmHost const& /*prob_parm*/,
  ProbParmDevice const& /*prob_parm_d*/)
{
  if (!init_parts) {
    return;
  }
  const int level = 0;
  amrex::ParmParse pp("prob");
  int numRedist = -1;
  amrex::IntVect partNum(AMREX_D_DECL(100, 100, 100));
  // Find the number of redistributions during particle initialization
  pp.query("init_redist", numRedist);
  pp.query("num_particles", partNum);
  amrex::RealVect partVel;
  std::array<amrex::Real, AMREX_SPACEDIM> pvel = {0.0};
  pp.query<amrex::Real>("part_vel", pvel);
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    partVel[dir] = pvel[dir];
  }
  amrex::Real partDia, partTemp;
  pp.get("part_dia", partDia);
  pp.get("part_temp", partTemp);
  std::array<amrex::Real, SPRAY_FUEL_NUM> partY = {0.0};
  partY[0] = 1.;
  uniformSprayInit(
    partNum, partVel, partDia, partTemp, partY.begin(), level, numRedist);
}
