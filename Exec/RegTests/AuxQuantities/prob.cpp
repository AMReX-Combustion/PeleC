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
    pp.query("U", PeleC::h_prob_parm_device->U);
    pp.query("Pr", PeleC::h_prob_parm_device->Pr);
    pp.query("aux_xy_lo", PeleC::h_prob_parm_device->aux_xy_lo);
    pp.query("aux_length", PeleC::h_prob_parm_device->aux_length);
    pp.query("aux_height", PeleC::h_prob_parm_device->aux_height);
    pp.query("aux_srcstrength", PeleC::h_prob_parm_device->aux_srcstrength);
  }

  // Characteristic lengths
  amrex::Real L = (probhi[0] - problo[0]);
  amrex::Real H = probhi[1] - problo[1];

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

  PeleC::h_prob_parm_device->Ma = PeleC::h_prob_parm_device->U / cs;

  // Transport properties
  auto& trans_parm = PeleC::trans_parms.host_parm();
  trans_parm.const_bulk_viscosity = 0.0;
  trans_parm.const_diffusivity = 0.0;
  trans_parm.const_viscosity = PeleC::h_prob_parm_device->rho *
                               PeleC::h_prob_parm_device->U * H /
                               PeleC::h_prob_parm_device->Re;
  trans_parm.const_conductivity =
    trans_parm.const_viscosity * cp / PeleC::h_prob_parm_device->Pr;
  PeleC::trans_parms.sync_to_device();

  // Output IC
  amrex::Print() << "\nGenerate ic.txt" << std::endl;
  std::ofstream ofs("ic.txt", std::ofstream::out);
  amrex::Print(ofs)
    << "L = " << L << "\nH = " << H << "\np = " << PeleC::h_prob_parm_device->p
    << "\nT = " << PeleC::h_prob_parm_device->T << "\ngamma = " << eos.gamma
    << "\ncs = " << cs << "\nU = " << PeleC::h_prob_parm_device->U
    << "\nrho = " << PeleC::h_prob_parm_device->rho
    << "\nviscosity = " << trans_parm.const_viscosity
    << "\nconductivity = " << trans_parm.const_conductivity
    << "\nRe = " << PeleC::h_prob_parm_device->Re
    << "\nMa = " << PeleC::h_prob_parm_device->Ma
    << "\nPr = " << PeleC::h_prob_parm_device->Pr << "\nNUM_AUX = " << NUM_AUX
    << "\naux_xy_lo = " << PeleC::h_prob_parm_device->aux_xy_lo
    << "\naux_length = " << PeleC::h_prob_parm_device->aux_length
    << "\naux_height = " << PeleC::h_prob_parm_device->aux_height << std::endl;
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
