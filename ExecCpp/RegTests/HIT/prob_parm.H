#ifndef _PROB_PARM_H_
#define _PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice : amrex::Gpu::Managed
{
  std::string iname;
   bool binfmt = false;
   bool restart = false;
   amrex::Real lambda0 = 0.5;
   amrex::Real reynolds_lambda0 = 100.0;
   amrex::Real mach_t0 = 0.1;
   amrex::Real prandtl = 0.71;
   int inres = 0;
   amrex::Real uin_norm = 1.0;
   amrex::Real L_x = 0.0;
   amrex::Real L_y = 0.0;
   amrex::Real L_z = 0.0;
   amrex::Real Linput = 0.0;
   amrex::Real k0 = 0.0;
   amrex::Real rho0 = 0.0;
   amrex::Real urms0 = 0.0;
   amrex::Real tau = 0.0;
   amrex::Real p0 = 1.013e6; // [erg cm^-3]
   amrex::Real T0 = 300.0;
   amrex::Real eint0 = 0.0;
  amrex::Gpu::ManagedVector<amrex::Real>* v_xinput = nullptr;
  amrex::Gpu::ManagedVector<amrex::Real>* v_uinput = nullptr;
  amrex::Gpu::ManagedVector<amrex::Real>* v_vinput = nullptr;
  amrex::Gpu::ManagedVector<amrex::Real>* v_winput = nullptr;
  amrex::Gpu::ManagedVector<amrex::Real>* v_xarray = nullptr;
  amrex::Gpu::ManagedVector<amrex::Real>* v_xdiff = nullptr;

   amrex::Real* xinput = nullptr;
   amrex::Real* uinput = nullptr;
   amrex::Real* vinput = nullptr;
   amrex::Real* winput = nullptr;
   amrex::Real* xarray = nullptr;
   amrex::Real* xdiff = nullptr;
};

struct ProbParmHost : amrex::Gpu::Managed
{
};

#endif
