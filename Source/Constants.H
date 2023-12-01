#ifndef CONSTANTS_H
#define CONSTANTS_H

namespace constants {
AMREX_GPU_HOST_DEVICE constexpr amrex::Real
PI()
{
  return amrex::Math::pi<amrex::Real>();
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

AMREX_GPU_HOST_DEVICE constexpr int
level_mask_interior()
{
  return 0;
}

AMREX_GPU_HOST_DEVICE constexpr int
level_mask_covered()
{
  return 1;
}

AMREX_GPU_HOST_DEVICE constexpr int
level_mask_notcovered()
{
  return 2;
}

AMREX_GPU_HOST_DEVICE constexpr int
level_mask_physbnd()
{
  return 3;
}

} // namespace constants
#endif
