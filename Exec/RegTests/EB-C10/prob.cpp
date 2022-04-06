#include "prob.H"

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
  const amrex::Real* problo,
  const amrex::Real* probhi)
{
  // Parse params
  {
    amrex::ParmParse pp("prob");
    pp.query("p", PeleC::h_prob_parm_device->p);
    pp.query("T", PeleC::h_prob_parm_device->T);
    pp.query("Re", PeleC::h_prob_parm_device->Re);
    pp.query("Ma", PeleC::h_prob_parm_device->Ma);
    pp.query("Pr", PeleC::h_prob_parm_device->Pr);
  }

  {
    amrex::ParmParse pp("eb2");
    pp.query("cylinder_radius", PeleC::h_prob_parm_device->radius);
  }

  amrex::Real L = (probhi[0] - problo[0]);

  amrex::Real cp = 0.0;
  amrex::Real cs = 0.0;
  PeleC::h_prob_parm_device->massfrac[0] = 1.0;

  auto eos = pele::physics::PhysicsType::eos();
  eos.PYT2RE(
    PeleC::h_prob_parm_device->p, PeleC::h_prob_parm_device->massfrac.begin(),
    PeleC::h_prob_parm_device->T, PeleC::h_prob_parm_device->rho,
    PeleC::h_prob_parm_device->eint);
  eos.RTY2Cs(
    PeleC::h_prob_parm_device->rho, PeleC::h_prob_parm_device->T,
    PeleC::h_prob_parm_device->massfrac.begin(), cs);
  eos.TY2Cp(
    PeleC::h_prob_parm_device->T, PeleC::h_prob_parm_device->massfrac.begin(),
    cp);

  PeleC::h_prob_parm_device->umax = PeleC::h_prob_parm_device->Ma * cs;
  PeleC::h_prob_parm_device->uavg = 0.5 * PeleC::h_prob_parm_device->umax;

  auto& trans_parm = PeleC::trans_parms.host_trans_parm();
  trans_parm.const_bulk_viscosity = 0.0;
  trans_parm.const_diffusivity = 0.0;
  trans_parm.const_viscosity = PeleC::h_prob_parm_device->rho *
                               PeleC::h_prob_parm_device->umax * L /
                               PeleC::h_prob_parm_device->Re;
  trans_parm.const_conductivity =
    trans_parm.const_viscosity * cp / PeleC::h_prob_parm_device->Pr;
  PeleC::trans_parms.sync_to_device();

  PeleC::h_prob_parm_device->G =
    PeleC::h_prob_parm_device->umax * 4 * trans_parm.const_viscosity /
    (PeleC::h_prob_parm_device->radius * PeleC::h_prob_parm_device->radius);
  PeleC::h_prob_parm_device->dpdx = -PeleC::h_prob_parm_device->G;

  // Output IC
  std::ofstream ofs("ic.txt", std::ofstream::out);
  amrex::Print(ofs)
    << "L, rho, umax, p, T, gamma, mu, k, Re, Ma, Pr, dpdx, G, radius"
    << std::endl;
  amrex::Print(ofs).SetPrecision(17)
    << L << "," << PeleC::h_prob_parm_device->rho << ","
    << PeleC::h_prob_parm_device->umax << "," << PeleC::h_prob_parm_device->p
    << "," << PeleC::h_prob_parm_device->T << "," << eos.gamma << ","
    << trans_parm.const_viscosity << "," << trans_parm.const_conductivity << ","
    << PeleC::h_prob_parm_device->Re << "," << PeleC::h_prob_parm_device->Ma
    << "," << PeleC::h_prob_parm_device->Pr << ","
    << PeleC::h_prob_parm_device->dpdx << "," << PeleC::h_prob_parm_device->G
    << "," << PeleC::h_prob_parm_device->radius << std::endl;
  ofs.close();
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
