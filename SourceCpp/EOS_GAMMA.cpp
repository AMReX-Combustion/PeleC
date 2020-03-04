#include <AMReX_Arena.H>
#include <AMReX_ParmParse.H>

#include "EOS_GAMMA.H"

namespace EOS {

AMREX_GPU_DEVICE_MANAGED amrex::Real gamma = 1.4;

void
init()
{
  amrex::ParmParse pp("eos");
  pp.query("gamma", gamma);
}

} // namespace EOS
