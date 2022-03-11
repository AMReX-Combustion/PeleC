#ifndef _PROB_PARM_H_
#define _PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice
{
  amrex::Real alpha = 1e-6;
  amrex::Real sigma = 10.0;
  amrex::Real p = 100000.0;
  amrex::Real rho = 0.0014;
  amrex::Real T = 0.0;
  amrex::Real cs = 0.0;
  amrex::Real radius = 0.0;
  amrex::Real L = 0.0;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {1.0};
};

struct ProbParmHost
{
  ProbParmHost() = default;
};

#endif
