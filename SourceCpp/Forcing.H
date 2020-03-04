#ifndef _FORCING_H_
#define _FORCING_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include "IndexDefines.H"

namespace forcing_params {
extern AMREX_GPU_DEVICE_MANAGED amrex::Real u0;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real v0;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real w0;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real forcing;
} // namespace forcing_params

#endif
