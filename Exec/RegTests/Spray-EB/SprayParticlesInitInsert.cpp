
#include "SprayParticles.H"
#include "SprayInjection.H"
#include "prob.H"

bool
SprayParticleContainer::injectParticles(
  amrex::Real time,
  amrex::Real dt,
  int /*nstep*/,
  int lev,
  int /*finest_level*/)
{
  if (lev != 0) {
    return false;
  }

  bool do_inject = false;
  amrex::ParmParse ps("spray.jet1");
  ps.query("do_inject", do_inject);
  if (do_inject) {
    SprayJet* js = m_sprayJets[0].get();
    sprayInjection(time, js, dt, lev);
  }

  // Redistribute is done outside of this function
  return true;
}

void
SprayParticleContainer::InitSprayParticles(const bool /*init_parts*/)
{
  std::string jet_name = "jet1";
  bool do_inject = false;
  amrex::ParmParse ps("spray.jet1");
  ps.query("do_inject", do_inject);
  if (do_inject) {
    m_sprayJets.push_back(std::make_unique<SprayJet>(jet_name, Geom(0)));
  }
  return;
}
