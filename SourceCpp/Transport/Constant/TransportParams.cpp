#include "TransportParams.H"

namespace transport_params {

AMREX_GPU_DEVICE_MANAGED amrex::Real const_viscosity = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real const_bulk_viscosity = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real const_diffusivity = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real const_conductivity = 0.0;

void
init()
{
  amrex::ParmParse pp("transport");
  pp.query("const_viscosity", const_viscosity);
  pp.query("const_bulk_viscosity", const_bulk_viscosity);
  pp.query("const_conductivity", const_conductivity);
  pp.query("const_diffusivity", const_diffusivity);
}

void
finalize()
{
}

} // namespace transport_params
