#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IntVect.H>

struct ProbParmDevice
{
  amrex::Real reynolds = 1600.0;
  amrex::Real mach = 0.1;
  bool convecting = false;
  amrex::Real L = 1.0;
  amrex::Real p0 = 1.013e6; // [erg cm^-3]
  amrex::Real T0 = 1000.0;
  amrex::Real rho0 = 0.0;
  amrex::Real v0 = 0.0;
  amrex::Real Y_O2 = 0.233;
  amrex::Real Y_N2 = 0.767;
};

struct ProbParmHost
{
  amrex::Real partTemp = 300.;
  amrex::Real partDia = 1.E-3;
};

#endif
