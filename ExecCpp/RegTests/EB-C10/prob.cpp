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
  const amrex_real* problo,
  const amrex_real* probhi)
{
  // Parse params
  {
    amrex::ParmParse pp("prob");
    pp.query("p", PeleC::prob_parm_device->p);
    pp.query("T", PeleC::prob_parm_device->T);
    pp.query("Re", PeleC::prob_parm_device->Re);
    pp.query("Ma", PeleC::prob_parm_device->Ma);
    pp.query("Pr", PeleC::prob_parm_device->Pr);
  }

  {
    amrex::ParmParse pp("eb2");
    pp.query("cylinder_radius", PeleC::prob_parm_device->radius);
  }

  amrex::Real L = (probhi[0] - problo[0]);

  amrex::Real cp = 0.0;
  amrex::Real cs = 0.0;
  PeleC::prob_parm_device->massfrac[0] = 1.0;

  EOS::PYT2RE(
    PeleC::prob_parm_device->p, PeleC::prob_parm_device->massfrac.begin(),
    PeleC::prob_parm_device->T, PeleC::prob_parm_device->rho,
    PeleC::prob_parm_device->eint);
  EOS::RTY2Cs(
    PeleC::prob_parm_device->rho, PeleC::prob_parm_device->T,
    PeleC::prob_parm_device->massfrac.begin(), cs);
  EOS::TY2Cp(
    PeleC::prob_parm_device->T, PeleC::prob_parm_device->massfrac.begin(), cp);

  PeleC::prob_parm_device->umax = PeleC::prob_parm_device->Ma * cs;
  PeleC::prob_parm_device->uavg = 0.5 * PeleC::prob_parm_device->umax;

  TransParm trans_parm;

  // Default
  trans_parm.const_viscosity = 0.0;
  trans_parm.const_bulk_viscosity = 0.0;
  trans_parm.const_conductivity = 0.0;
  trans_parm.const_diffusivity = 0.0;

  // User-specified
  {
    amrex::ParmParse pp("transport");
    pp.query("const_viscosity", trans_parm.const_viscosity);
    pp.query("const_bulk_viscosity", trans_parm.const_bulk_viscosity);
    pp.query("const_conductivity", trans_parm.const_conductivity);
    pp.query("const_diffusivity", trans_parm.const_diffusivity);
  }

  trans_parm.const_bulk_viscosity = 0.0;
  trans_parm.const_diffusivity = 0.0;
  trans_parm.const_viscosity = PeleC::prob_parm_device->rho *
                               PeleC::prob_parm_device->umax * L /
                               PeleC::prob_parm_device->Re;
  trans_parm.const_conductivity =
    trans_parm.const_viscosity * cp / PeleC::prob_parm_device->Pr;

#ifdef AMREX_USE_GPU
  amrex::Gpu::htod_memcpy(trans_parm_g, &trans_parm, sizeof(trans_parm));
#else
  std::memcpy(trans_parm_g, &trans_parm, sizeof(trans_parm));
#endif

  PeleC::prob_parm_device->G =
    PeleC::prob_parm_device->umax * 4 * trans_parm.const_viscosity /
    (PeleC::prob_parm_device->radius * PeleC::prob_parm_device->radius);
  PeleC::prob_parm_device->dpdx = -PeleC::prob_parm_device->G;

  // Output IC
  std::ofstream ofs("ic.txt", std::ofstream::out);
  amrex::Print(ofs)
    << "L, rho, umax, p, T, gamma, mu, k, Re, Ma, Pr, dpdx, G, radius"
    << std::endl;
  amrex::Print(ofs).SetPrecision(17)
    << L << "," << PeleC::prob_parm_device->rho << ","
    << PeleC::prob_parm_device->umax << "," << PeleC::prob_parm_device->p << ","
    << PeleC::prob_parm_device->T << "," << EOS::gamma << ","
    << trans_parm.const_viscosity << "," << trans_parm.const_conductivity << ","
    << PeleC::prob_parm_device->Re << "," << PeleC::prob_parm_device->Ma << ","
    << PeleC::prob_parm_device->Pr << "," << PeleC::prob_parm_device->dpdx
    << "," << PeleC::prob_parm_device->G << ","
    << PeleC::prob_parm_device->radius << std::endl;
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
