#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

struct ProbParmDevice
{
  // Chamber conditions
  amrex::Real P_mean = 6.0795e7;
  amrex::Real T_mean = 900.0;
  amrex::Real Y_CH4_chamber = 0.020594;
  amrex::Real Y_O2_chamber = 0.164305;

  // Injection parameters
  int nholes = 4;
  amrex::Real cone_angle = 45.0;
  amrex::Real centx = 0.0;
  amrex::Real centy = 0.0;
  amrex::Real r_circ = 0.15;
  amrex::Real r_hole = 0.0085;
  amrex::Real T_jet = 470.0;
  amrex::Real vel_jet = 2800.0;
  amrex::Real inj_start = 1e-6;
  amrex::Real inj_dur = .5e-3;
  amrex::Real tau = 0.01;
  amrex::Real Z = 0.45;
  bool turbulent_inflow = false;

  bool binfmt = false;
  bool hitIC = false;
  int inres = 0;
  amrex::Real uin_norm = 1.0;
  amrex::Real lscale = 1.0;
  amrex::Real offset = 0.0;
  amrex::Real urms0 = 0.0;
  amrex::Real Linput = 0.0;
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> win_lo = {0.0, 0.0, 0.0};
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> win_hi = {0.0, 0.0, 0.0};
  amrex::Real win_slope = 1.0;

  amrex::Real* d_xinput = nullptr;
  amrex::Real* d_uinput = nullptr;
  amrex::Real* d_vinput = nullptr;
  amrex::Real* d_winput = nullptr;
  amrex::Real* d_xarray = nullptr;
  amrex::Real* d_xdiff = nullptr;
};

struct ProbParmHost
{
  amrex::Vector<amrex::Real> h_xinput;
  amrex::Vector<amrex::Real> h_uinput;
  amrex::Vector<amrex::Real> h_vinput;
  amrex::Vector<amrex::Real> h_winput;
  amrex::Vector<amrex::Real> h_xarray;
  amrex::Vector<amrex::Real> h_xdiff;
  amrex::Gpu::DeviceVector<amrex::Real> xinput;
  amrex::Gpu::DeviceVector<amrex::Real> uinput;
  amrex::Gpu::DeviceVector<amrex::Real> vinput;
  amrex::Gpu::DeviceVector<amrex::Real> winput;
  amrex::Gpu::DeviceVector<amrex::Real> xarray;
  amrex::Gpu::DeviceVector<amrex::Real> xdiff;
  std::string iname;

  ProbParmHost()
    : xinput(0), uinput(0), vinput(0), winput(0), xarray(0), xdiff(0)
  {
  }
};

#endif
