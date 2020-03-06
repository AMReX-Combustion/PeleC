#ifndef _TRANSPORT_PARAMS_H_
#define _TRANSPORT_PARAMS_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

namespace transport_params {

void init();
void finalize();

// AMREX_GPU_DEVICE void test ();

extern AMREX_GPU_DEVICE_MANAGED amrex::Real* wt;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real* iwt;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real* eps;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real* sig;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real* dip;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real* pol;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real* zrot;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real* fitmu;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real* fitlam;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real* fitdbin;
extern AMREX_GPU_DEVICE_MANAGED int* nlin;
extern AMREX_GPU_DEVICE_MANAGED int array_size;
extern AMREX_GPU_DEVICE_MANAGED int fit_length;

} // namespace transport_params

#endif
