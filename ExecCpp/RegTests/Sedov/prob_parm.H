#ifndef _PROB_PARM_H_
#define _PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

namespace ProbParm {
extern AMREX_GPU_DEVICE_MANAGED amrex::Real p_ambient;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real dens_ambient;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real exp_energy;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real r_init;
extern AMREX_GPU_DEVICE_MANAGED int nsub;
extern AMREX_GPU_DEVICE_MANAGED amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>
  split;
} // namespace ProbParm

#endif
