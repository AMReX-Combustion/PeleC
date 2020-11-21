#ifndef _PROB_PARM_H_
#define _PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

namespace ProbParm {
extern AMREX_GPU_DEVICE_MANAGED amrex::Real p_init;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real T_init;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real rho_init;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real e_init;
extern AMREX_GPU_DEVICE_MANAGED amrex::GpuArray<amrex::Real, NUM_SPECIES>
  massfrac;
} // namespace ProbParm

#endif
