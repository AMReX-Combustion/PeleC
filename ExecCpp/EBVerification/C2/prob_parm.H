#ifndef _PROB_PARM_H_
#define _PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

namespace ProbParm {
extern AMREX_GPU_DEVICE_MANAGED amrex::Real p_l;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real rho_l;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real rhoe_l;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real T_l;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real p_r;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real rho_r;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real rhoe_r;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real T_r;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real angle;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real L;
extern AMREX_GPU_DEVICE_MANAGED std::string gasL;
extern AMREX_GPU_DEVICE_MANAGED std::string gasR;
extern AMREX_GPU_DEVICE_MANAGED int left_gas_id;
extern AMREX_GPU_DEVICE_MANAGED int right_gas_id;
} // namespace ProbParm

#endif
