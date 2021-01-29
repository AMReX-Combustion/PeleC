#ifndef _PROB_PARM_H_
#define _PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice : amrex::Gpu::Managed
{
  AMREX_GPU_DEVICE_MANAGED amrex::Real reynolds = 1.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real mach = 1.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real prandtl = 1.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real rho_x_fact = 0.1;
  AMREX_GPU_DEVICE_MANAGED amrex::Real rho_y_fact = 1.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real rho_z_fact = 1.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real u_0_fact = 1.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real u_x_fact = 1.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real u_y_fact = 1.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real u_z_fact = 1.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real v_0_fact = 1.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real v_x_fact = 1.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real v_y_fact = 1.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real v_z_fact = 1.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real w_0_fact = 1.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real w_x_fact = 1.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real w_y_fact = 1.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real w_z_fact = 1.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real p_x_fact = 0.2;
  AMREX_GPU_DEVICE_MANAGED amrex::Real p_y_fact = 1.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real p_z_fact = 1.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real a_rhox = 2.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real a_rhoy = 4.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real a_rhoz = 6.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real a_ux = 2.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real a_uy = 4.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real a_uz = 6.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real a_vx = 4.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real a_vy = 2.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real a_vz = 6.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real a_wx = 6.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real a_wy = 4.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real a_wz = 2.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real a_px = 6.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real a_py = 2.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real a_pz = 4.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real L_x = 0.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real L_y = 0.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real L_z = 0.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real p0 = 1.013e6; // [erg cm^-3]
  AMREX_GPU_DEVICE_MANAGED amrex::Real T0 = 300.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real rho0 = 0.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real u0 = 0.0;
};

struct ProbParmHost : amrex::Gpu::Managed
{
};

#endif
