#ifndef _PROB_PARM_H_
#define _PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice : amrex::Gpu::Managed
{
  std::string iname;
  AMREX_GPU_DEVICE_MANAGED bool binfmt = false;
  AMREX_GPU_DEVICE_MANAGED bool restart = false;
  AMREX_GPU_DEVICE_MANAGED amrex::Real lambda0 = 0.5;
  AMREX_GPU_DEVICE_MANAGED amrex::Real reynolds_lambda0 = 100.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real mach_t0 = 0.1;
  AMREX_GPU_DEVICE_MANAGED amrex::Real prandtl = 0.71;
  AMREX_GPU_DEVICE_MANAGED int inres = 0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real uin_norm = 1.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real L_x = 0.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real L_y = 0.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real L_z = 0.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real Linput = 0.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real k0 = 0.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real rho0 = 0.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real urms0 = 0.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real tau = 0.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real p0 = 1.013e6; // [erg cm^-3]
  AMREX_GPU_DEVICE_MANAGED amrex::Real T0 = 300.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real eint0 = 0.0;
  amrex::Gpu::ManagedVector<amrex::Real>* v_xinput = nullptr;
  amrex::Gpu::ManagedVector<amrex::Real>* v_uinput = nullptr;
  amrex::Gpu::ManagedVector<amrex::Real>* v_vinput = nullptr;
  amrex::Gpu::ManagedVector<amrex::Real>* v_winput = nullptr;
  amrex::Gpu::ManagedVector<amrex::Real>* v_xarray = nullptr;
  amrex::Gpu::ManagedVector<amrex::Real>* v_xdiff = nullptr;

  AMREX_GPU_DEVICE_MANAGED amrex::Real* xinput = nullptr;
  AMREX_GPU_DEVICE_MANAGED amrex::Real* uinput = nullptr;
  AMREX_GPU_DEVICE_MANAGED amrex::Real* vinput = nullptr;
  AMREX_GPU_DEVICE_MANAGED amrex::Real* winput = nullptr;
  AMREX_GPU_DEVICE_MANAGED amrex::Real* xarray = nullptr;
  AMREX_GPU_DEVICE_MANAGED amrex::Real* xdiff = nullptr;
};

struct ProbParmHost : amrex::Gpu::Managed
{
};

#endif
