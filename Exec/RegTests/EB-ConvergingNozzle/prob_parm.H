#ifndef _PROB_PARM_H_
#define _PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice
{
  amrex::Real d_cyl = 7.542;
  amrex::Real l_nozzle = 4.789;
  amrex::Real d_nozzle = 1.780;
  amrex::Real l_combustor = 23.711;
  amrex::Real l_exit = 1.0;
  amrex::Real d_swirler = 2.8;
  amrex::Real u_swirl_ax = 0.0;
  amrex::Real T_swirl = 300.0;
  amrex::Real p0 = 1.013e6; // [erg cm^-3]
  amrex::Real T0 = 300.0;
};

struct ProbParmHost
{
  ProbParmHost() = default;
};

#endif
