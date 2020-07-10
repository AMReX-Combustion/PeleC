#include <AMReX_ParmParse.H>

#include "PeleC.H"
#include "Tagging.H"

namespace TaggingParm {
AMREX_GPU_DEVICE_MANAGED amrex::Real denerr = 1.0e10;
AMREX_GPU_DEVICE_MANAGED int max_denerr_lev = 10;
AMREX_GPU_DEVICE_MANAGED amrex::Real dengrad = 1.0e10;
AMREX_GPU_DEVICE_MANAGED int max_dengrad_lev = 10;

AMREX_GPU_DEVICE_MANAGED amrex::Real presserr = 1.0e10;
AMREX_GPU_DEVICE_MANAGED int max_presserr_lev = 10;
AMREX_GPU_DEVICE_MANAGED amrex::Real pressgrad = 1.0e10;
AMREX_GPU_DEVICE_MANAGED int max_pressgrad_lev = 10;

AMREX_GPU_DEVICE_MANAGED amrex::Real velerr = 1.0e10;
AMREX_GPU_DEVICE_MANAGED int max_velerr_lev = 10;
AMREX_GPU_DEVICE_MANAGED amrex::Real velgrad = 1.0e10;
AMREX_GPU_DEVICE_MANAGED int max_velgrad_lev = 10;

AMREX_GPU_DEVICE_MANAGED amrex::Real vorterr = 1.0e10;
AMREX_GPU_DEVICE_MANAGED int max_vorterr_lev = 10;

AMREX_GPU_DEVICE_MANAGED amrex::Real temperr = 1.0e10;
AMREX_GPU_DEVICE_MANAGED int max_temperr_lev = 10;
AMREX_GPU_DEVICE_MANAGED amrex::Real tempgrad = 1.0e10;
AMREX_GPU_DEVICE_MANAGED int max_tempgrad_lev = 10;

AMREX_GPU_DEVICE_MANAGED amrex::Real ftracerr = 1.0e10;
AMREX_GPU_DEVICE_MANAGED int max_ftracerr_lev = 10;
AMREX_GPU_DEVICE_MANAGED amrex::Real ftracgrad = 1.0e10;
AMREX_GPU_DEVICE_MANAGED int max_ftracgrad_lev = 10;

AMREX_GPU_DEVICE_MANAGED amrex::Real vfracerr = 1.0e10;
AMREX_GPU_DEVICE_MANAGED int max_vfracerr_lev = 10;
} // namespace TaggingParm

void
PeleC::read_tagging_params()
{
  amrex::ParmParse pp("tagging");

  pp.query("denerr", TaggingParm::denerr);
  pp.query("max_denerr_lev", TaggingParm::max_denerr_lev);
  pp.query("dengrad", TaggingParm::dengrad);
  pp.query("max_dengrad_lev", TaggingParm::max_dengrad_lev);

  pp.query("presserr", TaggingParm::presserr);
  pp.query("max_presserr_lev", TaggingParm::max_presserr_lev);
  pp.query("pressgrad", TaggingParm::pressgrad);
  pp.query("max_pressgrad_lev", TaggingParm::max_pressgrad_lev);

  pp.query("velerr", TaggingParm::velerr);
  pp.query("max_velerr_lev", TaggingParm::max_velerr_lev);
  pp.query("velgrad", TaggingParm::velgrad);
  pp.query("max_velgrad_lev", TaggingParm::max_velgrad_lev);

  pp.query("vorterr", TaggingParm::vorterr);
  pp.query("max_vorterr_lev", TaggingParm::max_vorterr_lev);

  pp.query("temperr", TaggingParm::temperr);
  pp.query("max_temperr_lev", TaggingParm::max_temperr_lev);
  pp.query("tempgrad", TaggingParm::tempgrad);
  pp.query("max_tempgrad_lev", TaggingParm::max_tempgrad_lev);

  pp.query("ftracerr", TaggingParm::ftracerr);
  pp.query("max_ftracerr_lev", TaggingParm::max_ftracerr_lev);
  pp.query("ftracgrad", TaggingParm::ftracgrad);
  pp.query("max_ftracgrad_lev", TaggingParm::max_ftracgrad_lev);

  pp.query("vfracerr", TaggingParm::vfracerr);
  pp.query("max_vfracerr_lev", TaggingParm::max_vfracerr_lev);
}
