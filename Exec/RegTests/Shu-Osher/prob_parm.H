#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice
{
  amrex::Real p_l = 10.33333;   // left pressure (erg/cc)
  amrex::Real u_l = 2.629369;   // left velocity (cm/s)
  amrex::Real rho_l = 3.857143; // left density (g/cc)
  amrex::Real rhoe_l;
  amrex::Real T_l = 1.0;
  amrex::Real p_r = 1.0;        // right pressure (erg/cc)
  amrex::Real u_r = 0.0;        // right velocity (cm/s)
  amrex::Real rho_r_base = 1.0; // right density (g/cc)
  amrex::Real rho_r_amp = 0.2;
  amrex::Real rho_r_osc = 5.0;
  amrex::Real rhoe_r;
  amrex::Real T_r = 1.0;
  amrex::Real frac = 0.1; // fraction of the domain for the interface
  int idir = 1;           // direction across which to jump
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> split;
};

struct ProbParmHost
{
  ProbParmHost() = default;
};
#endif
