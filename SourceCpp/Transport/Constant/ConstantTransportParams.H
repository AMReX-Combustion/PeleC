#ifndef _CONSTANT_TRANSPORT_PARAMS_H_
#define _CONSTANT_TRANSPORT_PARAMS_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

namespace transport_params {

extern AMREX_GPU_DEVICE_MANAGED amrex::Real const_viscosity;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real const_bulk_viscosity;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real const_diffusivity;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real const_conductivity;

void init();
void finalize();

} // namespace transport_params

#endif
