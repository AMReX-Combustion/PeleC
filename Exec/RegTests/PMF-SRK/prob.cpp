#include "prob.H"

void
init_bc()
{
  amrex::Real vt;
  amrex::Real ek;
  amrex::Real T;
  amrex::Real rho;
  amrex::Real e;
  amrex::Real molefrac[NUM_SPECIES];
  amrex::Real massfrac[NUM_SPECIES];
  amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {{0.0}};

  if (PeleC::h_prob_parm_device->phi_in < 0) {
    const amrex::Real yl = 0.0;
    const amrex::Real yr = 0.0;
    pele::physics::PMF::pmf(
      PeleC::h_prob_parm_device->pmf_container, yl, yr, pmf_vals);
    amrex::Real mysum = 0.0;
    for (int n = 0; n < NUM_SPECIES; n++) {
      molefrac[n] = amrex::max<amrex::Real>(0.0, pmf_vals[3 + n]);
      mysum += molefrac[n];
    }
    molefrac[N2_ID] = 1.0 - (mysum - molefrac[N2_ID]);
    T = pmf_vals[0];
    PeleC::h_prob_parm_device->vn_in = pmf_vals[1];
  } else {
    const amrex::Real a = 0.5;
    for (amrex::Real& n : molefrac) {
      n = 0.0;
    }
    molefrac[O2_ID] =
      1.0 / (1.0 + PeleC::h_prob_parm_device->phi_in / a + 0.79 / 0.21);
    molefrac[H2_ID] = PeleC::h_prob_parm_device->phi_in * molefrac[O2_ID] / a;
    molefrac[N2_ID] = 1.0 - molefrac[H2_ID] - molefrac[O2_ID];
    T = PeleC::h_prob_parm_device->T_in;
  }
  const amrex::Real p = PeleC::h_prob_parm_device->pamb;

  auto eos = pele::physics::PhysicsType::eos();
  eos.X2Y(molefrac, massfrac);
  eos.PYT2RE(p, massfrac, T, rho, e);

  vt = PeleC::h_prob_parm_device->vn_in;
  ek = 0.5 * (vt * vt);

  PeleC::h_prob_parm_device->fuel_state[URHO] = rho;
  PeleC::h_prob_parm_device->fuel_state[UMX] = 0.0;
  PeleC::h_prob_parm_device->fuel_state[UMY] = rho * vt;
  PeleC::h_prob_parm_device->fuel_state[UMZ] = 0.0;
  PeleC::h_prob_parm_device->fuel_state[UEINT] = rho * e;
  PeleC::h_prob_parm_device->fuel_state[UEDEN] = rho * (e + ek);
  PeleC::h_prob_parm_device->fuel_state[UTEMP] = T;
  for (int n = 0; n < NUM_SPECIES; n++) {
    PeleC::h_prob_parm_device->fuel_state[UFS + n - 1] = rho * massfrac[n];
  }
}

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
  amrex::ParmParse pp("prob");
  pp.query("pamb", PeleC::h_prob_parm_device->pamb);
  pp.query("phi_in", PeleC::h_prob_parm_device->phi_in);
  pp.query("T_in", PeleC::h_prob_parm_device->T_in);
  pp.query("vn_in", PeleC::h_prob_parm_device->vn_in);
  pp.query("pertmag", PeleC::h_prob_parm_device->pertmag);

  PeleC::h_prob_parm_device->L[0] = probhi[0] - problo[0];
  PeleC::h_prob_parm_device->L[1] = probhi[1] - problo[1];
  PeleC::h_prob_parm_device->L[2] = probhi[2] - problo[2];

  PeleC::prob_parm_host->pmf_data.initialize();
  PeleC::h_prob_parm_device->pmf_container =
    PeleC::prob_parm_host->pmf_data.getDeviceData();

  init_bc();
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
