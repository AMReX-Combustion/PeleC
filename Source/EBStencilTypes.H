#ifndef _EBStencilTypes_H
#define _EBStencilTypes_H

#include <AMReX_REAL.H>
#include <AMReX_IntVect.H>

#ifndef BL_LANG_FORT

static amrex::Box stencil_volume_box(amrex::IntVect(D_DECL(-1,-1,-1)),
                                     amrex::IntVect(D_DECL(+1,+1,+1)));

#if BL_SPACEDIM == 2

struct VolSten
{
    amrex::Real val[3][3];
    amrex::IntVect iv;
};

struct FaceSten
{
    amrex::Real val[3];
    amrex::IntVect iv;
};

struct EBBndrySten
{
    amrex::Real val[3][3];
    amrex::Real bcval_sten;
    amrex::IntVect iv;
    amrex::IntVect iv_base;
};

#elif BL_SPACEDIM == 3

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

#endif

struct EBBndryGeom
{
  amrex::Real eb_normal[BL_SPACEDIM];
  amrex::Real eb_centroid[BL_SPACEDIM];
  amrex::Real eb_area;
  amrex::Real eb_vfrac;
  amrex::IntVect iv;

  bool operator<(const EBBndryGeom& rhs) const {
    return iv <rhs.iv;
  }

};

#endif

#endif
