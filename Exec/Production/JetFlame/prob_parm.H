#ifndef _PROB_PARM_H_
#define _PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

// using namespace amrex::literals;

struct ProbParmDevice
{
  amrex::Real P_mean = 1013250.0;
  amrex::Real inj_start = 1e-6;
  amrex::Real inj_dur = 5e-4;
  amrex::Real v_in = 2800;
  amrex::Real D = 1.7e-4;
  amrex::Real Z = 0.45;
  amrex::Real T_ox = 900;
  amrex::Real T_fu = 470;
  amrex::Real tau = 0.01;
  bool turbulence = false;

  amrex::Real Y_fuel[NUM_SPECIES] = {0.0};
  amrex::Real Y_ox[NUM_SPECIES] = {0.0};
  amrex::Real center_xy[2] = {0.0};
  amrex::Real U_fuel = {0.0};
  amrex::Real U_ox = {0.0};
};

struct ProbParmHost
{
  ProbParmHost() {}
};

#endif
