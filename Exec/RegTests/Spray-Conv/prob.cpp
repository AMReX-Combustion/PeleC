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
  const amrex_real* /*problo*/,
  const amrex_real* /*probhi*/)
{
  // Parse params
  amrex::ParmParse pp("prob");
  pp.query("mach", PeleC::h_prob_parm_device->mach);
  pp.query("ref_p", PeleC::h_prob_parm_device->p0);
  pp.query("ref_T", PeleC::h_prob_parm_device->T0);
  pp.query("init_O2", PeleC::h_prob_parm_device->Y_O2);
  pp.query("init_N2", PeleC::h_prob_parm_device->Y_N2);
  pp.query("init_redist", PeleC::prob_parm_host->numRedist);
  pp.query("num_particles", PeleC::prob_parm_host->partNum);
  std::array<amrex::Real, AMREX_SPACEDIM> pvel;
  pp.query<amrex::Real>("part_vel", pvel);
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    PeleC::prob_parm_host->partVel[dir] = pvel[dir];
  }
  pp.get("part_dia", PeleC::prob_parm_host->partDia);
  pp.get("part_temp", PeleC::prob_parm_host->partTemp);

  // Initial density, velocity, and material properties
  amrex::Real cs;
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  massfrac[O2_ID] = PeleC::h_prob_parm_device->Y_O2;
  massfrac[N2_ID] = PeleC::h_prob_parm_device->Y_N2;
  amrex::Real wbar, gamma0;
  auto eos = pele::physics::PhysicsType::eos();
  eos.TY2G(PeleC::h_prob_parm_device->T0, massfrac, gamma0);
  eos.Y2WBAR(massfrac, wbar);
  // Compute the speed of sound
  cs = std::sqrt(
    PeleC::h_prob_parm_device->T0 * gamma0 * pele::physics::Constants::RU /
    wbar);
  PeleC::h_prob_parm_device->v0 = PeleC::h_prob_parm_device->mach * cs;

  // Output IC
  if (amrex::ParallelDescriptor::IOProcessor()) {
    std::ofstream ofs("ic.txt", std::ofstream::out);
    amrex::Print(ofs) << "number of particles: "
                      << PeleC::prob_parm_host->partNum[0];
    for (int dir = 1; dir < AMREX_SPACEDIM; ++dir) {
      amrex::Print(ofs) << ", " << PeleC::prob_parm_host->partNum[dir];
    }
    amrex::Print(ofs) << std::endl;
    amrex::Print(ofs) << "particle velocity: "
                      << PeleC::prob_parm_host->partVel[0];
    for (int dir = 1; dir < AMREX_SPACEDIM; ++dir) {
      amrex::Print(ofs) << ", " << PeleC::prob_parm_host->partVel[dir];
    }
    amrex::Print(ofs) << std::endl;
    amrex::Print(ofs) << "p0: " << PeleC::h_prob_parm_device->p0 << std::endl;
    amrex::Print(ofs) << "cs: " << cs << std::endl;
    amrex::Print(ofs) << "U: " << PeleC::h_prob_parm_device->v0 << std::endl;
    amrex::Print(ofs) << "particle diameter: " << PeleC::prob_parm_host->partDia
                      << std::endl;
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
