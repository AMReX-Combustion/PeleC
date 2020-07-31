#ifndef _PROB_PARM_H_
#define _PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

namespace ProbParm {
extern AMREX_GPU_DEVICE_MANAGED amrex::Real pl;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real rhol;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real Tl;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real eintl;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real pr;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real rhor;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real Tr;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real eintr;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real angle;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real L;
extern AMREX_GPU_DEVICE_MANAGED amrex::GpuArray<amrex::Real, NUM_SPECIES>
  massfrac;
} // namespace ProbParm

#endif
