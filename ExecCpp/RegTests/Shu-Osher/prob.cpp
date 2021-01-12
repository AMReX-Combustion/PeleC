#include "prob.H"

namespace ProbParm {
AMREX_GPU_DEVICE_MANAGED amrex::Real p_l = 10.33333;   // left pressure (erg/cc)
AMREX_GPU_DEVICE_MANAGED amrex::Real u_l = 2.629369;   // left velocity (cm/s)
AMREX_GPU_DEVICE_MANAGED amrex::Real rho_l = 3.857143; // left density (g/cc)
AMREX_GPU_DEVICE_MANAGED amrex::Real rhoe_l;
AMREX_GPU_DEVICE_MANAGED amrex::Real T_l = 1.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real p_r = 1.0;     // right pressure (erg/cc)
AMREX_GPU_DEVICE_MANAGED amrex::Real u_r = 0.0;     // right velocity (cm/s)
AMREX_GPU_DEVICE_MANAGED amrex::Real rho_r_base = 1.0; // right density (g/cc)
AMREX_GPU_DEVICE_MANAGED amrex::Real rho_r_amp = 0.2;
AMREX_GPU_DEVICE_MANAGED amrex::Real rho_r_osc = 5.0;  
AMREX_GPU_DEVICE_MANAGED amrex::Real rhoe_r;
AMREX_GPU_DEVICE_MANAGED amrex::Real T_r = 1.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real frac =
  0.1; // fraction of the domain for the interface
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
  pp.query("rho_r_base", ProbParm::rho_r_base);
  pp.query("rho_r_amp", ProbParm::rho_r_amp);
  pp.query("rho_r_osc", ProbParm::rho_r_osc);
  pp.query("T_r", ProbParm::T_r);
  pp.query("frac", ProbParm::frac);
  pp.query("idir", ProbParm::idir);

  for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
    ProbParm::split[idir] = ProbParm::frac * (problo[idir] + probhi[idir]);
  }

  amrex::Real e_l;
  amrex::Real e_r /*, cs, cp */;
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  massfrac[0] = 1.0;

  EOS::RYP2E(ProbParm::rho_l, massfrac, ProbParm::p_l, e_l);
  EOS::EY2T(e_l, massfrac, ProbParm::T_l);
  ProbParm::rhoe_l = ProbParm::rho_l * e_l;
  EOS::RYP2E(ProbParm::rho_r_base, massfrac, ProbParm::p_r, e_r);
  EOS::EY2T(e_r, massfrac, ProbParm::T_r);
  ProbParm::rhoe_r = ProbParm::rho_r_base * e_r;
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
