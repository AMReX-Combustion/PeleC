#ifndef _PROB_PARM_H_
#define _PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>
#include "PMFData.H"

struct ProbParmDevice
{
  amrex::Real pamb = 1013250.0 * 100.0;
  amrex::Real phi_in = -0.2;
  amrex::Real T_in = 298.0;
  amrex::Real vn_in = 0.2;
  amrex::Real pertmag = 0.0;
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> L = {{1.0}};
  amrex::GpuArray<amrex::Real, NVAR> fuel_state = {{0.0}};
  pele::physics::PMF::PmfData::DataContainer* pmf_container;
};

struct ProbParmHost
{
  pele::physics::PMF::PmfData pmf_data;
  ProbParmHost() {}
};

#endif
