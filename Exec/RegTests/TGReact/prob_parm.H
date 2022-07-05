#ifndef _PROB_PARM_H_
#define _PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

namespace ProbParm {
extern AMREX_GPU_DEVICE_MANAGED amrex::Real rho0;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real u0;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real T0;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real p0;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real mu;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real omega_x;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real omega_y;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real omega_z;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real L_x;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real L_y;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real L_z;
extern std::string case_type;
extern AMREX_GPU_DEVICE_MANAGED int case_type_int;
extern std::string iname;
extern AMREX_GPU_DEVICE_MANAGED int nx;
extern AMREX_GPU_DEVICE_MANAGED int nvars;
extern amrex::Gpu::ManagedVector<amrex::Real>* v_input;
extern amrex::Gpu::ManagedVector<amrex::Real>* v_dxinput;
extern amrex::Gpu::ManagedVector<amrex::Real>* v_xarray;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real* input;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real* dxinput;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real* xarray;
extern AMREX_GPU_DEVICE_MANAGED int input_x;
extern AMREX_GPU_DEVICE_MANAGED int input_H2_init;
extern AMREX_GPU_DEVICE_MANAGED int input_O2_init;
extern AMREX_GPU_DEVICE_MANAGED int input_N2_init;
extern AMREX_GPU_DEVICE_MANAGED int input_Tad;
extern AMREX_GPU_DEVICE_MANAGED int input_p;
extern AMREX_GPU_DEVICE_MANAGED int input_H2;
extern AMREX_GPU_DEVICE_MANAGED int input_O2;
extern AMREX_GPU_DEVICE_MANAGED int input_N2;
extern AMREX_GPU_DEVICE_MANAGED int input_H2O;
extern AMREX_GPU_DEVICE_MANAGED int input_H;
extern AMREX_GPU_DEVICE_MANAGED int input_O;
extern AMREX_GPU_DEVICE_MANAGED int input_OH;
extern AMREX_GPU_DEVICE_MANAGED int input_HO2;
extern AMREX_GPU_DEVICE_MANAGED int input_H2O2;
} // namespace ProbParm

#endif
