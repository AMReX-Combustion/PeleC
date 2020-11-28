#include "prob.H"

namespace ProbParm {
AMREX_GPU_DEVICE_MANAGED amrex::Real alpha = 1e-6;
AMREX_GPU_DEVICE_MANAGED amrex::Real sigma = 10.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real p = 100000.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real rho = 0.0014;
AMREX_GPU_DEVICE_MANAGED amrex::Real T = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real cs = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real radius = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real L = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {
  1.0};
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
  pp.query("rho", ProbParm::rho);
  pp.query("p", ProbParm::p);
  pp.query("alpha", ProbParm::alpha);
  pp.query("sigma", ProbParm::sigma);

  amrex::ParmParse ppeb("eb2");
  ppeb.query("cylinder_radius", ProbParm::radius);

  ProbParm::L = (probhi[0] - problo[0]);

  ProbParm::massfrac[0] = 1.0;

  EOS::RYP2T(
    ProbParm::rho, ProbParm::massfrac.begin(), ProbParm::p, ProbParm::T);
  EOS::RPY2Cs(
    ProbParm::rho, ProbParm::p, ProbParm::massfrac.begin(), ProbParm::cs);

  // Output IC
  std::ofstream ofs("ic.txt", std::ofstream::out);
  amrex::Print(ofs) << "L, rho, p, T, gamma, cs, radius, alpha, sigma"
                    << std::endl;
  amrex::Print(ofs).SetPrecision(17)
    << ProbParm::L << "," << ProbParm::rho << "," << ProbParm::p << ","
    << ProbParm::T << "," << EOS::gamma << "," << ProbParm::cs << ","
    << ProbParm::radius << "," << ProbParm::alpha << "," << ProbParm::sigma
    << std::endl;
  ofs.close();
}
}

void
PeleC::problem_post_timestep()
{
  if (verbose <= 0)
    return;

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
      if (parent->NumDataLogs() > 1) {
        std::ostream& data_log2 = parent->DataLog(1);

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
  if (verbose <= 0)
    return;

  amrex::Real time = state[State_Type].curTime();

  if (level == 0) {
    if (amrex::ParallelDescriptor::IOProcessor()) {
      if (parent->NumDataLogs() > 1) {
        std::ostream& data_log2 = parent->DataLog(1);
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
