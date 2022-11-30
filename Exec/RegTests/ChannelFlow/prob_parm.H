#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice
{
  amrex::Real T_mean = 298.0;
  amrex::Real P_mean = 1013250.0;
  amrex::Real massfrac_O2 = 0.233;
  amrex::Real massfrac_N2 = 0.767;
  amrex::Real dpdx = 70.979; // = 709.79 Pa/m
};

struct ProbParmHost
{
  ProbParmHost() = default;
};

#endif
