#ifndef CONSTANTS_H
#define CONSTANTS_H

namespace constants {
AMREX_GPU_HOST_DEVICE constexpr amrex::Real
PI()
{
  return 3.14159265358979323846264338327950288;
}
AMREX_GPU_HOST_DEVICE constexpr amrex::Real
smallu()
{
  return 1.0e-12;
}
AMREX_GPU_HOST_DEVICE constexpr amrex::Real
small_num()
{
  return 1.0e-8;
}
AMREX_GPU_HOST_DEVICE constexpr amrex::Real
very_small_num()
{
  return std::numeric_limits<amrex::Real>::epsilon() * 1e-100;
}
} // namespace constants
#endif
