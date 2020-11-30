#include "prob.H"

namespace ProbParm {
AMREX_GPU_DEVICE_MANAGED amrex::Real p_l = 1.0;   // left pressure (erg/cc)
AMREX_GPU_DEVICE_MANAGED amrex::Real u_l = 0.0;   // left velocity (cm/s)
AMREX_GPU_DEVICE_MANAGED amrex::Real rho_l = 0.0; // left density (g/cc)
AMREX_GPU_DEVICE_MANAGED amrex::Real rhoe_l;
AMREX_GPU_DEVICE_MANAGED amrex::Real T_l = 1.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real p_r = 0.1;     // right pressure (erg/cc)
AMREX_GPU_DEVICE_MANAGED amrex::Real u_r = 0.0;     // right velocity (cm/s)
AMREX_GPU_DEVICE_MANAGED amrex::Real rho_r = 0.125; // right density (g/cc)
AMREX_GPU_DEVICE_MANAGED amrex::Real rhoe_r;
AMREX_GPU_DEVICE_MANAGED amrex::Real T_r = 1.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real frac =
  0.5; // fraction of the domain for the interface
AMREX_GPU_DEVICE_MANAGED bool use_Tinit =
  false; // optionally use T_l/r instead of p_l/r for initialization
AMREX_GPU_DEVICE_MANAGED int idir = 1; // direction across which to jump
AMREX_GPU_DEVICE_MANAGED amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> split;
} // namespace ProbParm

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
  amrex::ParmParse pp("prob");
  pp.query("p_l", ProbParm::p_l);
  pp.query("u_l", ProbParm::u_l);
  pp.query("rho_l", ProbParm::rho_l);
  pp.query("T_l", ProbParm::T_l);
  pp.query("p_r", ProbParm::p_r);
  pp.query("u_r", ProbParm::u_r);
  pp.query("rho_r", ProbParm::rho_r);
  pp.query("T_r", ProbParm::T_r);
  pp.query("frac", ProbParm::frac);
  pp.query("idir", ProbParm::idir);
  pp.query("use_Tinit", ProbParm::use_Tinit);

  for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
    ProbParm::split[idir] = ProbParm::frac * (problo[idir] + probhi[idir]);
  }

  amrex::Real e_l, e_r /*, cs, cp */;
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  massfrac[0] = 1.0;

  if (ProbParm::use_Tinit) {
    EOS::RTY2P(ProbParm::rho_l, ProbParm::T_l, massfrac, ProbParm::p_l);
    EOS::RYP2E(ProbParm::rho_l, massfrac, ProbParm::p_l, e_l);
    ProbParm::rhoe_l = ProbParm::rho_l * e_l;
    EOS::RTY2P(ProbParm::rho_r, ProbParm::T_r, massfrac, ProbParm::p_r);
    EOS::RYP2E(ProbParm::rho_r, massfrac, ProbParm::p_r, e_r);
    ProbParm::rhoe_r = ProbParm::rho_r * e_r;
  } else {
    EOS::RYP2E(ProbParm::rho_l, massfrac, ProbParm::p_l, e_l);
    EOS::EY2T(e_l, massfrac, ProbParm::T_l);
    ProbParm::rhoe_l = ProbParm::rho_l * e_l;
    EOS::RYP2E(ProbParm::rho_r, massfrac, ProbParm::p_r, e_r);
    EOS::EY2T(e_r, massfrac, ProbParm::T_r);
    ProbParm::rhoe_r = ProbParm::rho_r * e_r;
  }
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
