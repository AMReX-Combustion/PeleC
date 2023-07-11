#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

// using namespace amrex::literals;

struct ProbParmDevice
{
  amrex::Real Pres_amb = 1013250.0;
  amrex::Real Temp_amb = 1000.0;
  amrex::Real Yox_amb = 0.23;
  amrex::Real temp_jet = 300.0;
  amrex::Real vel_jet = 5000.0;
  amrex::Real Yox_channel = 0.23;
  amrex::Real Yfuel_jet = 1.0;

  // C1 case
  amrex::Real centx = 1.0;
  amrex::Real centz = 0.225;

  // C2 case
  // amrex::Real centx = 1.5;
  // amrex::Real centz = 0.0;

  amrex::Real r_hole = 0.05;
  int inject_fuel = 0;
  amrex::Real Mach_channel = 1.5;
  int init_type = 0;
  amrex::Real cavity_depth = 0.305;
};

struct ProbParmHost
{
  ProbParmHost() = default;
};

#endif
