#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice
{
  amrex::Real reynolds = 1600.0;
  amrex::Real mach = 0.1;
  amrex::Real prandtl = 0.71;
  bool convecting = false;
  amrex::Real omega_x = 1.0;
  amrex::Real omega_y = 1.0;
  amrex::Real omega_z = 1.0;
  amrex::Real L_x = 0.0;
  amrex::Real L_y = 0.0;
  amrex::Real L_z = 0.0;
  amrex::Real L = 0.0;
  amrex::Real p0 = 1.013e6; // [erg cm^-3]
  amrex::Real T0 = 300.0;
  amrex::Real rho0 = 0.0;
  amrex::Real v0 = 0.0;
};

struct ProbParmHost
{
  ProbParmHost() = default;
};

#endif
