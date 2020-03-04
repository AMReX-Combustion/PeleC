#ifndef _PROB_PARM_H_
#define _PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

namespace ProbParm {
extern AMREX_GPU_DEVICE_MANAGED amrex::Real p;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real T;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real rho;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real eint;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real vx_in;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real vy_in;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real Re_L;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real Pr;
extern AMREX_GPU_DEVICE_MANAGED amrex::GpuArray<amrex::Real, NUM_SPECIES>
  massfrac;
} // namespace ProbParm

#endif
