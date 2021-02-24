#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_

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
} // namespace constants
#endif
