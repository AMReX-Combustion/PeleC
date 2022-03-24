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
  amrex::ParmParse pp("prob");
  pp.query("rho", PeleC::h_prob_parm_device->rho);
  pp.query("p", PeleC::h_prob_parm_device->p);
  pp.query("alpha", PeleC::h_prob_parm_device->alpha);
  pp.query("sigma", PeleC::h_prob_parm_device->sigma);

  amrex::ParmParse ppeb("eb2");
  ppeb.query("cylinder_radius", PeleC::h_prob_parm_device->radius);

  PeleC::h_prob_parm_device->L = (probhi[0] - problo[0]);

  PeleC::h_prob_parm_device->massfrac[0] = 1.0;

  auto eos = pele::physics::PhysicsType::eos();
  eos.RYP2T(
    PeleC::h_prob_parm_device->rho, PeleC::h_prob_parm_device->massfrac.begin(),
    PeleC::h_prob_parm_device->p, PeleC::h_prob_parm_device->T);
  eos.RPY2Cs(
    PeleC::h_prob_parm_device->rho, PeleC::h_prob_parm_device->p,
    PeleC::h_prob_parm_device->massfrac.begin(), PeleC::h_prob_parm_device->cs);

  // Output IC
  std::ofstream ofs("ic.txt", std::ofstream::out);
  amrex::Print(ofs) << "L, rho, p, T, gamma, cs, radius, alpha, sigma"
                    << std::endl;
  amrex::Print(ofs).SetPrecision(17)
    << PeleC::h_prob_parm_device->L << "," << PeleC::h_prob_parm_device->rho
    << "," << PeleC::h_prob_parm_device->p << ","
    << PeleC::h_prob_parm_device->T << "," << eos.gamma << ","
    << PeleC::h_prob_parm_device->cs << "," << PeleC::h_prob_parm_device->radius
    << "," << PeleC::h_prob_parm_device->alpha << ","
    << PeleC::h_prob_parm_device->sigma << std::endl;
  ofs.close();
}
}

void
PeleC::problem_post_timestep()
{
  if (verbose <= 0) {
    return;
  }

  int finest_level = parent->finestLevel();
  amrex::Real time = state[State_Type].curTime();
  amrex::Real rho_err = 0.0;

  if (level == 0) {
    if (amrex::ParallelDescriptor::IOProcessor()) {
      amrex::Print() << "... Compute error problem post timestep" << std::endl;
    }

    // Calculate the errors
    for (int lev = 0; lev <= finest_level; lev++) {
      PeleC& pc_lev = getLevel(lev);

      bool local_flag = true;
      rho_err += pc_lev.volWgtSquaredSum("rhoerror", time, local_flag);
    }

    // Reductions
    amrex::ParallelDescriptor::ReduceRealSum(
      &rho_err, 1, amrex::ParallelDescriptor::IOProcessorNumber());

    // Get the norm and normalize it
    amrex::Real V = volume.sum(0, false);
    rho_err = std::sqrt(rho_err / V);

    if (amrex::ParallelDescriptor::IOProcessor()) {
      amrex::Print() << "TIME= " << time << " RHO ERROR  = " << rho_err << '\n';

      const int log_index = find_datalog_index("mmslog");
      if (log_index >= 0) {
        std::ostream& data_log2 = parent->DataLog(log_index);

        // Write the quantities at this time
        const int datwidth = 14;
        const int datprecision = 6;
        data_log2 << std::setw(datwidth) << time;
        data_log2 << std::setw(datwidth) << std::setprecision(datprecision)
                  << rho_err;
        data_log2 << std::endl;
      }
    }
  }
}

void
PeleC::problem_post_init()
{
  if (verbose <= 0) {
    return;
  }

  amrex::Real time = state[State_Type].curTime();

  if (level == 0) {
    if (amrex::ParallelDescriptor::IOProcessor()) {
      const int log_index = find_datalog_index("mmslog");
      if (log_index >= 0) {
        std::ostream& data_log2 = parent->DataLog(log_index);
        if (time == 0.0) {
          const int datwidth = 14;
          data_log2 << std::setw(datwidth) << "          time";
          data_log2 << std::setw(datwidth) << "       rho_err";
          data_log2 << std::endl;
        }
      }
    }
  }
}

void
PeleC::problem_post_restart()
{
}
