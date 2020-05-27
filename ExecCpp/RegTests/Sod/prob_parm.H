#ifndef _PROB_PARM_H_
#define _PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

namespace ProbParm {
extern AMREX_GPU_DEVICE_MANAGED amrex::Real p_l;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real u_l;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real rho_l;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real rhoe_l;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real T_l;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real p_r;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real u_r;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real rho_r;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real rhoe_r;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real T_r;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real frac;
extern AMREX_GPU_DEVICE_MANAGED bool use_Tinit;
extern AMREX_GPU_DEVICE_MANAGED int idir;
extern AMREX_GPU_DEVICE_MANAGED amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>
  split;
} // namespace ProbParm

#endif
