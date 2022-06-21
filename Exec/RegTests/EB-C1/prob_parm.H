#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice
{
  amrex::Real reynolds = 1.0;
  amrex::Real mach = 1.0;
  amrex::Real prandtl = 1.0;
  amrex::Real rho_x_fact = 0.1;
  amrex::Real rho_y_fact = 1.0;
  amrex::Real rho_z_fact = 1.0;
  amrex::Real u_0_fact = 0.0;
  amrex::Real v_0_fact = 0.0;
  amrex::Real w_0_fact = 0.0;
  amrex::Real u_r_fact = 0.1;
  amrex::Real v_r_fact = 0.1;
  amrex::Real w_r_fact = 0.1;
  amrex::Real p_r_fact = 0.0;
  amrex::Real a_rhox = 2.0;
  amrex::Real a_rhoy = 2.0;
  amrex::Real a_rhoz = 2.0;
  amrex::Real a_ur = 2.0;
  amrex::Real a_vr = 2.0;
  amrex::Real a_wr = 2.0;
  amrex::Real a_pr = 0.0;
  amrex::Real L_x = 0.0;
  amrex::Real L_y = 0.0;
  amrex::Real L_z = 0.0;
  amrex::Real p0 = 1.013e6; // [erg cm^-3]
  amrex::Real T0 = 300.0;
  amrex::Real rho0 = 0.0;
  amrex::Real u0 = 0.0;
};

struct ProbParmHost
{
  ProbParmHost() {}
};

#endif
