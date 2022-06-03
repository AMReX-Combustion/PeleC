#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice
{
  amrex::Real P_mean = 1013250.0;
  amrex::Real T_mean = 298.0;
  amrex::Real rvort = 0.1;
  amrex::Real xvort = 0.0;
  amrex::Real yvort = 0.0;
  amrex::Real forcevort = 3.0;
};

struct ProbParmHost
{
  ProbParmHost() {}
};

#endif
