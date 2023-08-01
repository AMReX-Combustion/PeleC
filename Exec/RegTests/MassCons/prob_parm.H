#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice
{
  amrex::Real T_mean = 298.0;
  amrex::Real P_mean = 1013250.0;
  amrex::Real u0 = 1000.0;
  amrex::Real v0 = 1000.0;
  amrex::Real w0 = 1000.0;
  int wall_type[AMREX_SPACEDIM] = {0};
};

struct ProbParmHost
{
  ProbParmHost() = default;
};

#endif
