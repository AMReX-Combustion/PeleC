#ifndef _PROB_PARM_H_
#define _PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

namespace ProbParm {
extern AMREX_GPU_DEVICE_MANAGED amrex::Real reynolds;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real mach;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real prandtl;
extern AMREX_GPU_DEVICE_MANAGED bool convecting;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real omega_x;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real omega_y;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real omega_z;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real L_x;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real L_y;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real L_z;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real L;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real p0;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real T0;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real rho0;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real v0;
} // namespace ProbParm

#endif
