#ifndef _PROB_PARM_H_
#define _PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice : amrex::Gpu::Managed
{
};

struct ProbParmHost : amrex::Gpu::Managed
{
};

namespace ProbParm {
} // namespace ProbParm

#endif
