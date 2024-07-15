#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice
{
  amrex::Real d_inner = 25.0;
  amrex::Real d_outer = 75.0;
  amrex::Real p0 = 1e6;
  amrex::Real T0 = 300.0;
  amrex::Real mu = 0.1;
  amrex::Real kappa = 1.0;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {1.0};
  int do_rf = 0;
  amrex::Real rf_omega = 0;
  amrex::Real rf_axis = 2;
  amrex::Real rf_axis_x = 0.0;
  amrex::Real rf_axis_y = 0.0;
  amrex::Real rf_axis_z = 0.0;
};

struct ProbParmHost
{
  ProbParmHost() = default;
};

#endif
