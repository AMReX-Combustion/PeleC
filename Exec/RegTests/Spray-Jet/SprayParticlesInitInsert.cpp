
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
  SprayJet* js = m_sprayJets[0].get();
  if (!js->jet_active(time)) {
    return false;
  }

  sprayInjection(time, js, dt, lev);

  // Redistribute is done outside of this function
  return true;
}

void
SprayParticleContainer::InitSprayParticles(const bool /*init_parts*/)
{
  m_sprayJets.resize(1);
  std::string jet_name = "jet1";
  m_sprayJets[0] = std::make_unique<SprayJet>(jet_name, Geom(0));
  return;
}
