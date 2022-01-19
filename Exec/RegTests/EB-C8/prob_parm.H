#ifndef _PROB_PARM_H_
#define _PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice
{
  amrex::Real p_l = 1.0;   // left pressure (erg/cc)
  amrex::Real rho_l = 1.0; // left density (g/cc)
  amrex::Real rhoe_l = 0.0;
  amrex::Real T_l = 1.0;
  amrex::Real u_l = 0.0;
  amrex::Real p_r = 0.1;     // right pressure (erg/cc)
  amrex::Real rho_r = 0.125; // right density (g/cc)
  amrex::Real rhoe_r = 0.0;
  amrex::Real T_r = 1.0;
  amrex::Real u_r = 0.0;
  amrex::Real angle = 0.0;
  amrex::Real L = 0.0;
  int left_gas_id = N2_ID;
  int right_gas_id = HE_ID;
};

struct ProbParmHost
{
  std::string gasL = "N2";
  std::string gasR = "HE";
  ProbParmHost() {}
};

#endif
