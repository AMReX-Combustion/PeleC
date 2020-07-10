#ifndef _PROB_PARM_H_
#define _PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

namespace ProbParm {
extern std::string iname;
extern AMREX_GPU_DEVICE_MANAGED bool binfmt;
extern AMREX_GPU_DEVICE_MANAGED bool restart;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real lambda0;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real reynolds_lambda0;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real mach_t0;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real prandtl;
extern AMREX_GPU_DEVICE_MANAGED int inres;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real uin_norm;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real L_x;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real L_y;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real L_z;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real Linput;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real k0;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real rho0;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real urms0;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real tau;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real p0;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real T0;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real eint0;

extern amrex::Gpu::ManagedVector<amrex::Real>* v_xinput;
extern amrex::Gpu::ManagedVector<amrex::Real>* v_uinput;
extern amrex::Gpu::ManagedVector<amrex::Real>* v_vinput;
extern amrex::Gpu::ManagedVector<amrex::Real>* v_winput;
extern amrex::Gpu::ManagedVector<amrex::Real>* v_xarray;
extern amrex::Gpu::ManagedVector<amrex::Real>* v_xdiff;

extern AMREX_GPU_DEVICE_MANAGED amrex::Real* xinput;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real* uinput;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real* vinput;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real* winput;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real* xarray;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real* xdiff;
} // namespace ProbParm

#endif
