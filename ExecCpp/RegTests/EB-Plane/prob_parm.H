#ifndef _PROB_PARM_H_
#define _PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice : amrex::Gpu::Managed
{
  AMREX_GPU_DEVICE_MANAGED amrex::Real p = 1013250.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real T = 0.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real rho = 0.00116;
  AMREX_GPU_DEVICE_MANAGED amrex::Real eint = 0.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real vx_in = 0.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real vy_in = 0.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real Re_L = 2500.0;
  AMREX_GPU_DEVICE_MANAGED amrex::Real Pr = 0.7;
  AMREX_GPU_DEVICE_MANAGED amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac =
    {1.0};
};

struct ProbParmHost : amrex::Gpu::Managed
{
};

#endif
