
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
  bool inject = false;
  for (int jindx = 0; jindx < m_sprayJets.size(); ++jindx) {
    SprayJet* js = m_sprayJets[jindx].get();
    if (js->jet_active(time)) {
      sprayInjection(time, js, dt, lev);
      inject = true;
    }
  }

  // Redistribute is done outside of this function
  return inject;
}

void
SprayParticleContainer::InitSprayParticles(const bool /*init_parts*/)
{
  amrex::ParmParse ps("spray");
  amrex::Real jet_vel = 6.E4;
  amrex::Real jet_dia = 1.E-3;
  amrex::Real mass_flow_rate = 2.3E-3;
  amrex::Real part_temp = 300.;
  amrex::Real jet_start_time = 0.;
  amrex::Real jet_end_time = 10000.;
  amrex::Real spread_angle = 20.;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> Y_jet = {{0.0}};
  ps.get("jet_vel", jet_vel);
  ps.query("jet_start", jet_start_time);
  ps.query("jet_end", jet_end_time);
  // The cells are divided by this value when prescribing the jet inlet
  ps.get("jet_dia", jet_dia);
  ps.get("T", part_temp);
  ps.get("mass_flow_rate", mass_flow_rate);
  ps.get("spread_angle", spread_angle);
  std::vector<int> jets_per_dir(AMREX_SPACEDIM);
  ps.getarr("jets_per_dir", jets_per_dir);
  std::vector<amrex::Real> in_Y_jet(SPRAY_FUEL_NUM, 0.);
  in_Y_jet[0] = 1.;
  ps.queryarr("jet_mass_fracs", in_Y_jet);
  amrex::Real sumY = 0.;
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    Y_jet[spf] = in_Y_jet[spf];
    sumY += in_Y_jet[spf];
  }
  if (std::abs(sumY - 1.) > 1.E-8) {
    amrex::Abort("'jet_mass_fracs' must sum to 1");
  }
  // Total number of jets
  int num_jets = AMREX_D_TERM(jets_per_dir[0], *1, *jets_per_dir[2]);
  m_sprayJets.resize(num_jets);
  const auto plo = this->m_gdb->Geom(0).ProbLoArray();
  const auto phi = this->m_gdb->Geom(0).ProbHiArray();
  amrex::Real div_lenx =
    (phi[0] - plo[0]) / (static_cast<amrex::Real>(jets_per_dir[0]));
  int jetz = 1;
#if AMREX_SPACEDIM == 3
  amrex::Real div_lenz = 0.;
  amrex::Real zlo = 0.;
  div_lenz = (phi[2] - plo[2]) / (static_cast<amrex::Real>(jets_per_dir[2]));
  zlo = plo[2];
  jetz = jets_per_dir[2];
#endif
  amrex::Real yloc = plo[1];
  int jindx = 0;
  amrex::RealVect jet_norm(AMREX_D_DECL(0., 1., 0.));
  std::string dist_type = "Uniform";
  // Create each jet struct
  for (int i = 0; i < jets_per_dir[0]; ++i) {
    amrex::Real xloc = plo[0] + div_lenx * (static_cast<amrex::Real>(i) + 0.5);
    for (int k = 0; k < jetz; ++k) {
      amrex::RealVect jet_cent(AMREX_D_DECL(
        xloc, yloc, zlo + div_lenz * (static_cast<amrex::Real>(k) + 0.5)));
      std::string jet_name = "jet" + std::to_string(jindx);
      m_sprayJets[jindx] = std::make_unique<SprayJet>(
        jet_name, Geom(0), jet_cent, jet_norm, spread_angle, jet_dia, jet_vel,
        mass_flow_rate, part_temp, Y_jet, dist_type, jet_start_time,
        jet_end_time);
      jindx++;
    }
  }
  return;
}
