
#include "SprayInjection.H"
#include <PeleC.H>
#include "prob.H"

bool
SprayParticleContainer::injectParticles(
  amrex::Real /*time*/,
  amrex::Real /*dt*/,
  int /*nstep*/,
  int /*lev*/,
  int /*finest_level*/)
{
  return false;
}

void
SprayParticleContainer::InitSprayParticles(const bool init_parts)
{
  if (!init_parts) {
    return;
  }
  const int level = 0;

  amrex::Real partRho;
  amrex::ParmParse ppp("particles");
  ppp.get("fuel_rho", partRho);
  amrex::ParmParse pp("prob");
  amrex::IntVect partNum(AMREX_D_DECL(100, 100, 100));
  pp.query("num_particles", partNum);
  amrex::RealVect partVel = amrex::RealVect::TheZeroVector();

  amrex::Real Stmod = 5.;
  amrex::Real rhoRatio = 1000.;
  amrex::Real mach = 0.1;
  // const amrex::Real Pr = 0.71;
  amrex::Real reynolds, p0, T0, rho0, L, v0;
  bool convecting;
  amrex::Real mu, cs, cp, St_num;
  // Parse params
  const amrex::Real* problo = Geom(level).ProbLo();
  const amrex::Real* probhi = Geom(level).ProbHi();
  prob_parser(
    problo, probhi, reynolds, mach, convecting, p0, T0, rho0, Stmod, rhoRatio,
    L, v0, mu, cs, cp, St_num);

  amrex::Real refRho = rho0;
  amrex::Real refL = L;
  amrex::Real refU = v0;
  amrex::Real tau_g = refL / refU;
  amrex::Real tau_d = St_num * tau_g;
  amrex::Real partDia = std::sqrt(18. * mu * tau_d / partRho);
  amrex::Real Re_d = partDia * refU * rho0 / mu;
  diameter_solver(rho0, mu, tau_d, partRho, refRho, refU, partDia, Re_d);
  amrex::Real partTemp;
  pp.query("ref_T", partTemp);
  int numRedist = 1;
  std::array<amrex::Real, SPRAY_FUEL_NUM> partY = {0.0};
  partY[0] = 1.;
  uniformSprayInit(
    partNum, partVel, partDia, partTemp, partY.begin(), level, numRedist);
}
