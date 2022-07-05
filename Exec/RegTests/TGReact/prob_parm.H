#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

struct ProbParmDevice
{
  amrex::Real rho0 = 1.0e-3;
  amrex::Real u0 = 100.0;
  amrex::Real T0 = 300.0;
  amrex::Real p0 = 1.013e6; // [erg cm^-3]
  amrex::Real mu = 0.00625;
  amrex::Real omega_x = 2.0 * constants::PI();
  amrex::Real omega_y = 2.0 * constants::PI();
  amrex::Real omega_z = 2.0 * constants::PI();
  amrex::Real L_x = 0.0;
  amrex::Real L_y = 0.0;
  amrex::Real L_z = 0.0;
  int case_type_int = 0;
  int nx = 0;
  int nvars = 15;
  amrex::Real* d_input = nullptr;
  amrex::Real* d_xarray = nullptr;
  amrex::Real* d_dxinput = nullptr;
  int input_x = 0;
  int input_H2_init = 1;
  int input_O2_init = 2;
  int input_N2_init = 3;
  int input_Tad = 4;
  int input_p = 5;
  int input_H2 = 6;
  int input_O2 = 7;
  int input_N2 = 8;
  int input_H2O = 9;
  int input_H = 10;
  int input_O = 11;
  int input_OH = 12;
  int input_HO2 = 13;
  int input_H2O2 = 14;
};

struct ProbParmHost
{
  std::string iname = "";
  std::string case_type = "cold";
  amrex::Vector<amrex::Real> h_input;
  amrex::Vector<amrex::Real> h_xarray;
  amrex::Vector<amrex::Real> h_dxinput;
  amrex::Gpu::DeviceVector<amrex::Real> input;
  amrex::Gpu::DeviceVector<amrex::Real> xarray;
  amrex::Gpu::DeviceVector<amrex::Real> dxinput;
  ProbParmHost() : input(0), xarray(0), dxinput(0) {}
};

#endif
