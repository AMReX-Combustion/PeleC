#ifndef _PROB_PARM_H_
#define _PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice
{
  amrex::Real p_init = 1013250.0; // 1 atm
  amrex::Real T_init = 940.0;
  amrex::Real T_inner = 1200.0;
  amrex::Real T_outer = 1000.0;
  amrex::Real rho_init = 0.0;
  amrex::Real e_init = 0.0;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {0.0};
  amrex::Real rad_temp = 0.025;
  amrex::Real rad_jet = 0.01;
  amrex::Real u_jet = 0.0;
  amrex::Real T_jet = 300.0;
};

struct ProbParmHost
{
  ProbParmHost() {}
};

#endif
