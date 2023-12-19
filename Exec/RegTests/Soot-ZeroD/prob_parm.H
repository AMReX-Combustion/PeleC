#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

struct ProbParmDevice
{
  amrex::Real p0 = 1.013e6; // [erg cm^-3]
  amrex::Real T0 = 300.0;
  amrex::Real rho0 = 0.0;
  amrex::Real Y_N2 = 0.6693174609128006;
  amrex::Real Y_O2 = 0.20333443771796805;
  amrex::Real Y_Fuel = 0.12734810136923122;
  amrex::GpuArray<amrex::Real, NUM_SOOT_MOMENTS + 1> soot_vals;
  int fuelIndx = -1;
};

struct ProbParmHost
{
};

#endif
