#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice
{
  amrex::Real rho = 1.0;
  amrex::Real p = 1.0;
  amrex::Real mach = 1.0;
  amrex::Real T = 300.0;
  amrex::Real c_s = 0.0;
  amrex::Real eint = 0.0;
  amrex::Real r_inner = 1.0;
  amrex::Real r_outer = 2.0;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {1.0};
};

struct ProbParmHost
{
  ProbParmHost() {}
};

#endif
