#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice
{
  amrex::Real p0 = 1.0;
  amrex::Real rho0 = 1.0;
  amrex::Real u0 = 1.0;
  amrex::Real A = 0.2;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {0.0};
};

struct ProbParmHost
{
  ProbParmHost() = default;
};

#endif
