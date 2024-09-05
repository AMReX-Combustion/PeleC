#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice
{
  amrex::Real p = 1013250.0;
  amrex::Real T = 300.0;
  amrex::Real rho = 0.0;
  amrex::Real eint = 0.0;
  amrex::Real U = 0.0;
  amrex::Real Re = 100.0;
  amrex::Real Ma = 0.1;
  amrex::Real Pr = 0.7;
  amrex::Real aux_IC = 10.0;
  amrex::Real aux_xy_lo = -1.5;
  amrex::Real aux_length = 0.5;
  amrex::Real aux_height = 1.0;
  amrex::Real aux_srcstrength = 0.0;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {1.0};
};

struct ProbParmHost
{
  ProbParmHost() = default;
};

#endif
