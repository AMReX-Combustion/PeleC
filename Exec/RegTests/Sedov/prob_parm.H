#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice
{
  amrex::Real p_ambient;
  amrex::Real dens_ambient;
  amrex::Real exp_energy;
  amrex::Real r_init;
  int nsub;
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> split;
};

struct ProbParmHost
{
  ProbParmHost() = default;
};
#endif
