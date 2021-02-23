#ifndef _EBStencilTypes_H_
#define _EBStencilTypes_H_

#include <AMReX_REAL.H>
#include <AMReX_IntVect.H>

static amrex::Box stencil_volume_box(
  amrex::IntVect(AMREX_D_DECL(-1, -1, -1)),
  amrex::IntVect(AMREX_D_DECL(+1, +1, +1)));

struct VolSten
{
  amrex::Real val[3][3][3];
  amrex::IntVect iv;
};

struct FaceSten
{
  amrex::Real val[3][3];
  amrex::IntVect iv;
};

struct EBBndrySten
{
  amrex::Real val[3][3][3];
  amrex::Real bcval_sten;
  amrex::IntVect iv;
  amrex::IntVect iv_base;
};

struct EBBndryGeom
{
  amrex::Real eb_normal[AMREX_SPACEDIM];
  amrex::Real eb_centroid[AMREX_SPACEDIM];
  amrex::Real eb_area;
  amrex::Real eb_vfrac;
  amrex::IntVect iv;

  AMREX_GPU_DEVICE
  bool operator<(const EBBndryGeom& rhs) const { return iv < rhs.iv; }
};

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
// Comparison operator for thrust sort
struct EBBndryGeomCmp
{
  AMREX_GPU_DEVICE
  bool operator()(const EBBndryGeom& a, const EBBndryGeom& b)
  {
    return a.iv < b.iv;
  }
};
#endif

#endif
