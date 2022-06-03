#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice
{
  amrex::Real p0 = 1.95e4;
  amrex::Real p1 = 7.42e6;
  amrex::Real rho0 = 3.29e-5;
  amrex::Real rho1 = 3.61e-4;
  amrex::Real T0 = 300.0;
  amrex::Real T1 = 300.0;
  amrex::Real e0 = 1.0;
  amrex::Real e1 = 1.0;
  amrex::Real x1 = 4.2;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {1.0};
};

struct ProbParmHost
{
  ProbParmHost() {}
};

#endif
