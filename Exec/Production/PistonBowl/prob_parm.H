#ifndef _PROB_PARM_H_
#define _PROB_PARM_H_

#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

struct ProbParmDevice
{
  amrex::Real Pres_domain = 10132500.0;
  amrex::Real Temp_domain = 500.0;
  amrex::Real Yfuel_domain = 0.;
  amrex::Real Yox_domain = 0.25;
  amrex::Real YN2_domain = 0.75;
  amrex::Real dens_jet = 0.0001;
  amrex::Real vel_jet = 5000.;
  amrex::Real Yox_jet = 1.0;
  amrex::Real Yfuel_jet = 0.0;
  amrex::Real YN2_jet = 0.0;
  amrex::Real centx = 2.5;
  amrex::Real centz = 2.5;
  amrex::Real r_circ = 3.0;
  amrex::Real r_hole = 0.4;
  int nholes = 4;
  amrex::Real cone_angle = 20.0;
  amrex::Real inj_time = 0.00014309;
  bool binfmt = false;
  bool hitIC = false;
  int inres = 0;
  amrex::Real uin_norm = 1.0;
  amrex::Real lscale = 1.0;
  amrex::Real offset = 0.0;
  amrex::Real urms0 = 0.0;
  amrex::Real Linput = 0.0;
  int turb_inflow_type = 0;
  int inflowNtime = 0;
  int nr = 0;
  int nt = 0;
  amrex::Real inflowFreq = 0.0;
  amrex::Real timeFromInflow = 0.0;
  amrex::Real timeInflowMax = 0.0;
  amrex::Real thetaMax = 0.0;
  amrex::Real* d_xinput = nullptr;
  amrex::Real* d_uinput = nullptr;
  amrex::Real* d_vinput = nullptr;
  amrex::Real* d_winput = nullptr;
  amrex::Real* d_xarray = nullptr;
  amrex::Real* d_xdiff = nullptr;
  amrex::Real* d_timeInput = nullptr;
  amrex::Real* d_rM = nullptr;
  amrex::Real* d_thetaM = nullptr;
  amrex::Real* d_Uz = nullptr;
  amrex::Real* d_Ur = nullptr;
  amrex::Real* d_Ut = nullptr;
};

struct ProbParmHost
{
  amrex::Vector<amrex::Real> h_xinput;
  amrex::Vector<amrex::Real> h_uinput;
  amrex::Vector<amrex::Real> h_vinput;
  amrex::Vector<amrex::Real> h_winput;
  amrex::Vector<amrex::Real> h_xarray;
  amrex::Vector<amrex::Real> h_xdiff;
  amrex::Vector<amrex::Real> h_timeInput;
  amrex::Vector<amrex::Real> h_rM;
  amrex::Vector<amrex::Real> h_thetaM;
  amrex::Vector<amrex::Real> h_Uz;
  amrex::Vector<amrex::Real> h_Ur;
  amrex::Vector<amrex::Real> h_Ut;
  amrex::Gpu::DeviceVector<amrex::Real> xinput;
  amrex::Gpu::DeviceVector<amrex::Real> uinput;
  amrex::Gpu::DeviceVector<amrex::Real> vinput;
  amrex::Gpu::DeviceVector<amrex::Real> winput;
  amrex::Gpu::DeviceVector<amrex::Real> xarray;
  amrex::Gpu::DeviceVector<amrex::Real> xdiff;
  amrex::Gpu::DeviceVector<amrex::Real> timeInput;
  amrex::Gpu::DeviceVector<amrex::Real> rM;
  amrex::Gpu::DeviceVector<amrex::Real> thetaM;
  amrex::Gpu::DeviceVector<amrex::Real> Uz;
  amrex::Gpu::DeviceVector<amrex::Real> Ur;
  amrex::Gpu::DeviceVector<amrex::Real> Ut;
  std::string iname = "";
  std::string inflowProfileFile = "";
  amrex::Vector<amrex::Real> rInput;
  amrex::Vector<amrex::Real> thetaInput;

  ProbParmHost()
    : xinput(0),
      uinput(0),
      vinput(0),
      winput(0),
      xarray(0),
      xdiff(0),
      timeInput(0),
      rM(0),
      thetaM(0),
      Uz(0),
      Ur(0),
      Ut(0)
  {
  }
};

#endif
