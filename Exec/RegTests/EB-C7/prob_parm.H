#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice
{
  amrex::Real pl = 1.0;
  amrex::Real rhol = 1.0;
  amrex::Real Tl = 0.0;
  amrex::Real eintl = 0.0;
  amrex::Real pr = 0.1;
  amrex::Real rhor = 0.125;
  amrex::Real eintr = 0.0;
  amrex::Real Tr = 0.0;
  amrex::Real angle = 0.0;
  amrex::Real L = 0.0;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {1.0};
};

struct ProbParmHost
{
  ProbParmHost() = default;
};

#endif
