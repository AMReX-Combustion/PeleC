#ifndef _PROB_PARM_H_
#define _PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice : amrex::Gpu::Managed
{
  AMREX_GPU_DEVICE_MANAGED amrex::Real p_l = 1.0;   // left pressure (erg/cc)
  AMREX_GPU_DEVICE_MANAGED amrex::Real u_l = 0.0;   // left velocity (cm/s)
  AMREX_GPU_DEVICE_MANAGED amrex::Real rho_l = 0.0; // left density (g/cc)
  AMREX_GPU_DEVICE_MANAGED amrex::Real rhoe_l;
  AMREX_GPU_DEVICE_MANAGED amrex::Real T_l = 1.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real p_r = 0.1;     // right pressure (erg/cc)
  AMREX_GPU_DEVICE_MANAGED amrex::Real u_r = 0.0;     // right velocity (cm/s)
  AMREX_GPU_DEVICE_MANAGED amrex::Real rho_r = 0.125; // right density (g/cc)
  AMREX_GPU_DEVICE_MANAGED amrex::Real rhoe_r;
  AMREX_GPU_DEVICE_MANAGED amrex::Real T_r = 1.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real frac =
    0.5; // fraction of the domain for the interface
  AMREX_GPU_DEVICE_MANAGED int idir = 1; // direction across which to jump
  AMREX_GPU_DEVICE_MANAGED amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> split;
  AMREX_GPU_DEVICE_MANAGED int left_gas_id = N2_ID;
  AMREX_GPU_DEVICE_MANAGED int right_gas_id = HE_ID;
  std::string gasL = "N2";
  std::string gasR = "HE";
};

struct ProbParmHost : amrex::Gpu::Managed
{
};
#endif
