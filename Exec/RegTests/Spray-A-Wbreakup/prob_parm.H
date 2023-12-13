#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice
{
  amrex::Real p0 = 1.013e6; // [erg cm^-3]
  amrex::Real T0 = 300.0;
  amrex::Real v0 = 0.;
  amrex::Real YX_O2 = 0.233;
  amrex::Real YX_N2 = 0.767;
  amrex::Real YX_H2O = 0.;
  amrex::Real YX_CO2 = 0.;
  bool mol_fracs = false;
};

struct ProbParmHost
{
};

#endif
