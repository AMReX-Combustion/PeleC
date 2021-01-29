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
AMREX_GPU_DEVICE_MANAGED int idir = 1; // direction across which to jump
AMREX_GPU_DEVICE_MANAGED amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> split;
AMREX_GPU_DEVICE_MANAGED int left_gas_id = N2_ID;
AMREX_GPU_DEVICE_MANAGED int right_gas_id = HE_ID;
std::string gasL = "N2";
std::string gasR = "HE";
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
  const amrex_real* /*problo*/,
  const amrex_real* /*probhi*/)
{
  // Parse params
  amrex::ParmParse pp("prob");
  pp.query("p_l", PeleC::prob_parm_device->p_l);
  pp.query("u_l", PeleC::prob_parm_device->u_l);
  pp.query("rho_l", PeleC::prob_parm_device->rho_l);
  pp.query("p_r", PeleC::prob_parm_device->p_r);
  pp.query("u_r", PeleC::prob_parm_device->u_r);
  pp.query("rho_r", PeleC::prob_parm_device->rho_r);
  pp.query("frac", PeleC::prob_parm_device->frac);
  pp.query("idir", PeleC::prob_parm_device->idir);
  pp.get("left_gas", PeleC::prob_parm_device->gasL);
  pp.get("right_gas", PeleC::prob_parm_device->gasR);

  for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
    PeleC::prob_parm_device->split[idir] = 0.5;
  }

  if (PeleC::prob_parm_device->gasL == "N2") {
    PeleC::prob_parm_device->left_gas_id = N2_ID;
    PeleC::prob_parm_device->right_gas_id = HE_ID;
  } else {
    PeleC::prob_parm_device->left_gas_id = HE_ID;
    PeleC::prob_parm_device->right_gas_id = N2_ID;
  }

  amrex::Real e_l;
  amrex::Real e_r /*, cs, cp*/;
  amrex::Real massfrac_l[NUM_SPECIES] = {0.0};
  amrex::Real massfrac_r[NUM_SPECIES] = {0.0};
  massfrac_l[PeleC::prob_parm_device->left_gas_id] = 1.0;
  massfrac_r[PeleC::prob_parm_device->right_gas_id] = 1.0;

  EOS::RYP2E(PeleC::prob_parm_device->rho_l, massfrac_l, PeleC::prob_parm_device->p_l, e_l);
  EOS::EY2T(e_l, massfrac_l, PeleC::prob_parm_device->T_l);
  PeleC::prob_parm_device->rhoe_l = PeleC::prob_parm_device->rho_l * e_l;

  EOS::RYP2E(PeleC::prob_parm_device->rho_r, massfrac_r, PeleC::prob_parm_device->p_r, e_r);
  EOS::EY2T(e_r, massfrac_r, PeleC::prob_parm_device->T_r);
  PeleC::prob_parm_device->rhoe_r = PeleC::prob_parm_device->rho_r * e_r;
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
