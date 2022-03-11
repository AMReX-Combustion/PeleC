#ifndef _PROB_PARM_H_
#define _PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice
{
  amrex::Real p = 1013250.0;
  amrex::Real T = 0.0;
  amrex::Real rho = 0.00116;
  amrex::Real eint = 0.0;
  amrex::Real vx_in = 9000.0;
  amrex::Real vy_in = 0.0;
  amrex::Real Re_L = 2500.0;
  amrex::Real Pr = 0.7;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {1.0};
};

struct ProbParmHost
{
  ProbParmHost() = default;
};

#endif
