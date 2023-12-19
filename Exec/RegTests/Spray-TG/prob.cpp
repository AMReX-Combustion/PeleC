#include "prob.H"

void
pc_prob_close()
{
}

void
diameter_solver(
  const amrex::Real rho0,
  const amrex::Real mu,
  const amrex::Real tau_d,
  const amrex::Real partRho,
  const amrex::Real refRho,
  const amrex::Real refU,
  amrex::Real& dia,
  amrex::Real& Re_d)
{

  amrex::Real error = 1000.;
  const amrex::Real tol = 1.E-6;
  int k = 0;
  const int maxIter = 500;
  while (error > tol) {
    amrex::Real oldDia = dia;
    Re_d = dia * refU * rho0 / mu;
    amrex::Real C_D = 24. / Re_d;
    if (Re_d > 1.)
      C_D *= (1. + std::pow(Re_d, 2. / 3.) / 6.);
    dia = 0.75 * refRho * C_D * refU * tau_d / partRho;
    error = std::abs(oldDia - dia) / dia;
    k++;
    if (k > maxIter) {
      amrex::Abort("Failed to converge a particle diameter");
    }
  }
}

void
prob_parser(
  const amrex_real* problo,
  const amrex_real* probhi,
  amrex::Real& reynolds,
  amrex::Real& mach,
  bool& convecting,
  amrex::Real& p0,
  amrex::Real& T0,
  amrex::Real& rho0,
  amrex::Real& Stmod,
  amrex::Real& rhoRatio,
  amrex::Real& L,
  amrex::Real& v0,
  amrex::Real& mu,
  amrex::Real& cs,
  amrex::Real& cp,
  amrex::Real& St_num)
{

  // Parse params
  amrex::ParmParse pp("prob");
  pp.query("reynolds", reynolds);
  pp.query("mach", mach);
  pp.query("convecting", convecting);
  pp.query("ref_p", p0);
  pp.query("ref_T", T0);
  pp.query("st_mod", Stmod);
  pp.query("density_ratio", rhoRatio);

  // Define the length scale
  L = probhi[0] - problo[0];

  // Initial density, velocity, and material properties
  amrex::Real eint;
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  massfrac[N2_ID] = 1.;
  auto eos = pele::physics::PhysicsType::eos();
  eos.PYT2RE(p0, massfrac, T0, rho0, eint);
  eos.RTY2Cs(rho0, T0, massfrac, cs);
  eos.TY2Cp(T0, massfrac, cp);

  v0 = mach * cs;
  mu = rho0 * v0 * L / reynolds;

  St_num = Stmod / (8. * M_PI);
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
  amrex::Real Stmod = 5.;
  amrex::Real rhoRatio = 1000.;
  amrex::Real mach = 0.1;
  const amrex::Real Pr = 0.71;
  amrex::Real mu, cs, cp, St_num;
  // Parse params
  prob_parser(
    problo, probhi, PeleC::h_prob_parm_device->reynolds, mach,
    PeleC::h_prob_parm_device->convecting, PeleC::h_prob_parm_device->p0,
    PeleC::h_prob_parm_device->T0, PeleC::h_prob_parm_device->rho0, Stmod,
    rhoRatio, PeleC::h_prob_parm_device->L, PeleC::h_prob_parm_device->v0, mu,
    cs, cp, St_num);

  auto& trans_parm = PeleC::trans_parms.host_trans_parm();
  trans_parm.const_bulk_viscosity = 0.0;
  trans_parm.const_diffusivity = 0.0;
  trans_parm.const_viscosity = mu;
  trans_parm.const_conductivity = mu * cp / Pr;
  PeleC::trans_parms.sync_to_device();

  amrex::Real partRho;
  amrex::ParmParse ppp("particles");
  ppp.get("fuel_rho", partRho);
  amrex::Real refRho = PeleC::h_prob_parm_device->rho0;
  if (std::abs(rhoRatio - partRho / refRho) > 10.) {
    amrex::Print() << "Restart solution with particles.fuel_rho = "
                   << refRho * rhoRatio << std::endl;
    amrex::Abort();
  }
  // Time scale for Eulerian phase
  amrex::Real refL = PeleC::h_prob_parm_device->L;
  amrex::Real refU = PeleC::h_prob_parm_device->v0;
  amrex::Real tau_g = refL / refU;
  amrex::Real tau_d = St_num * tau_g;
  amrex::Real dia = std::sqrt(18. * mu * tau_d / partRho);
  amrex::Real Re_d = dia * refU * PeleC::h_prob_parm_device->rho0 / mu;
  // Iteratively solve for the diameter of the spray droplets
  diameter_solver(
    PeleC::h_prob_parm_device->rho0, mu, tau_d, partRho, refRho, refU, dia,
    Re_d);
  PeleC::prob_parm_host->partDia = dia;
  PeleC::prob_parm_host->partTemp = PeleC::h_prob_parm_device->T0;
  // Output IC
  if (amrex::ParallelDescriptor::IOProcessor()) {
    std::ofstream ofs("ic.txt", std::ofstream::out);
    amrex::Print(ofs) << "rho0: " << PeleC::h_prob_parm_device->rho0
                      << std::endl;
    amrex::Print(ofs) << "cs: " << cs << std::endl;
    amrex::Print(ofs) << "U: " << PeleC::h_prob_parm_device->v0 << std::endl;
    amrex::Print(ofs) << "mu: " << trans_parm.const_viscosity << std::endl;
    amrex::Print(ofs) << "Re: " << PeleC::h_prob_parm_device->reynolds
                      << std::endl;
    amrex::Print(ofs) << "Stokes number: " << Stmod << "*Stc" << std::endl;
    amrex::Print(ofs) << "particle diameter: " << PeleC::prob_parm_host->partDia
                      << std::endl;
    amrex::Print(ofs) << "tau_d: " << tau_d << std::endl;
    amrex::Print(ofs) << "Re_d: " << Re_d << std::endl;
    amrex::Print(ofs) << "tau_g: " << tau_g << std::endl;
    ofs.close();
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
