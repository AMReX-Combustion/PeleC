#ifndef EBSTENCILTYPES_H
#define EBSTENCILTYPES_H

#include <AMReX_REAL.H>
#include <AMReX_IntVect.H>

static amrex::Box stencil_volume_box(
  amrex::IntVect(AMREX_D_DECL(-1, -1, -1)),
  amrex::IntVect(AMREX_D_DECL(+1, +1, +1)));

struct VolSten
{
  amrex::Real val AMREX_D_TERM([3], [3], [3]);
  amrex::IntVect iv;
};

struct FaceSten
{
  amrex::Real val AMREX_D_TERM(, [3], [3]);
  amrex::IntVect iv;
};

struct EBBndrySten
{
  amrex::Real val AMREX_D_TERM([3], [3], [3]);
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

#ifdef AMREX_USE_GPU
// Comparison operator for device sort
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
