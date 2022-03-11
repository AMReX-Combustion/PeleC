#ifndef _PROB_PARM_H_
#define _PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice
{
  amrex::Real p_l = 1.0;   // left pressure (erg/cc)
  amrex::Real u_l = 0.0;   // left velocity (cm/s)
  amrex::Real rho_l = 0.0; // left density (g/cc)
  amrex::Real rhoe_l;
  amrex::Real T_l = 1.0;
  amrex::Real p_r = 0.1;     // right pressure (erg/cc)
  amrex::Real u_r = 0.0;     // right velocity (cm/s)
  amrex::Real rho_r = 0.125; // right density (g/cc)
  amrex::Real rhoe_r;
  amrex::Real T_r = 1.0;
  amrex::Real frac = 0.5; // fraction of the domain for the interface
  bool use_Tinit =
    false;      // optionally use T_l/r instead of p_l/r for initialization
  int idir = 1; // direction across which to jump
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> split;
};

struct ProbParmHost
{
  ProbParmHost() = default;
};
#endif
