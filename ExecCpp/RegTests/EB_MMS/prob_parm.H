#ifndef _PROB_PARM_H_
#define _PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

namespace ProbParm {
extern AMREX_GPU_DEVICE_MANAGED amrex::Real reynolds;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real mach;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real prandtl;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real rho_x_fact;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real rho_y_fact;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real rho_z_fact;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real u_0_fact;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real v_0_fact;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real w_0_fact;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real u_r_fact;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real v_r_fact;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real w_r_fact;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real p_r_fact;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real a_rhox;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real a_rhoy;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real a_rhoz;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real a_ur;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real a_vr;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real a_wr;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real a_pr;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real L_x;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real L_y;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real L_z;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real p0;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real T0;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real rho0;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real u0;
} // namespace ProbParm

#endif
