#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice
{
  amrex::Real p_init = 1013250.0; // 1 atm
  amrex::Real Y_init_H2 = 0.06;
  amrex::Real Y_init_O2 = 0.5;
  amrex::Real Y_init_N2 = 0.44;
  amrex::Real T_init = 940.0;
  amrex::Real rho_init = 0.0;
  amrex::Real e_init = 0.0;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {0.0};
};

struct ProbParmHost
{
  ProbParmHost() = default;
};
#endif
