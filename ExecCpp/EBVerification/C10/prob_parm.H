#ifndef _PROB_PARM_H_
#define _PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

namespace ProbParm {
extern AMREX_GPU_DEVICE_MANAGED amrex::Real p;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real dpdx;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real T;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real rho;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real eint;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real umax;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real uavg;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real Ma;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real Pr;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real radius;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real G;
extern AMREX_GPU_DEVICE_MANAGED amrex::GpuArray<amrex::Real, NUM_SPECIES>
  massfrac;
} // namespace ProbParm

#endif
