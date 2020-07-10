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
extern AMREX_GPU_DEVICE_MANAGED amrex::Real u_x_fact;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real u_y_fact;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real u_z_fact;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real v_0_fact;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real v_x_fact;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real v_y_fact;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real v_z_fact;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real w_0_fact;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real w_x_fact;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real w_y_fact;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real w_z_fact;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real p_x_fact;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real p_y_fact;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real p_z_fact;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real a_rhox;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real a_rhoy;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real a_rhoz;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real a_ux;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real a_uy;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real a_uz;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real a_vx;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real a_vy;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real a_vz;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real a_wx;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real a_wy;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real a_wz;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real a_px;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real a_py;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real a_pz;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real L_x;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real L_y;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real L_z;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real p0;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real T0;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real rho0;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real u0;
} // namespace ProbParm

#endif
