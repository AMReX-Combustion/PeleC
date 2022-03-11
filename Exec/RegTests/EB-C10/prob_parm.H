#ifndef _PROB_PARM_H_
#define _PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice
{
  amrex::Real p = 1013250.0;
  amrex::Real dpdx = 0.0;
  amrex::Real T = 300.0;
  amrex::Real rho = 0.0;
  amrex::Real eint = 0.0;
  amrex::Real umax = 0.0;
  amrex::Real uavg = 0.0;
  amrex::Real Re = 100.0;
  amrex::Real Ma = 0.1;
  amrex::Real Pr = 0.7;
  amrex::Real radius = 0.0;
  amrex::Real G = 0.0;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {1.0};
};

struct ProbParmHost
{
  ProbParmHost() = default;
};

#endif
