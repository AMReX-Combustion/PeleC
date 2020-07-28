#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>

#include "mechanism.h"

#include "EOS.H"
#include "prob_parm.H"
#include "Transport.H"

namespace ProbParm {
AMREX_GPU_DEVICE_MANAGED amrex::Real alpha = 1e-4;
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
  const int* init,
  const int* name,
  const int* namelen,
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
