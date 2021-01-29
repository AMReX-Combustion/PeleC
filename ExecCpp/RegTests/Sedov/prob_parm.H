#ifndef _PROB_PARM_H_
#define _PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice : amrex::Gpu::Managed
{
  AMREX_GPU_DEVICE_MANAGED amrex::Real p_ambient;
  AMREX_GPU_DEVICE_MANAGED amrex::Real dens_ambient;
  AMREX_GPU_DEVICE_MANAGED amrex::Real exp_energy;
  AMREX_GPU_DEVICE_MANAGED amrex::Real r_init;
  AMREX_GPU_DEVICE_MANAGED int nsub;
  AMREX_GPU_DEVICE_MANAGED amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> split;
};

struct ProbParmHost : amrex::Gpu::Managed
{
};
#endif
