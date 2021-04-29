#ifndef _EB_H_
#define _EB_H_

#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_BoxArray.H>
#include <AMReX_VisMF.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Union.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Complement.H>
#include <AMReX_EB2_IF_Scale.H>
#include <AMReX_EB2_IF_Translation.H>
#include <AMReX_EB2_IF_Lathe.H>
#include <AMReX_EB2_IF_Box.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Ellipsoid.H>
#include <AMReX_EB2_IF_Sphere.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Extrusion.H>
#include <AMReX_EB2_GeometryShop.H>
#include <AMReX_EB2_IF_Difference.H>
#include <AMReX_EBFluxRegister_C.H>

#include "PeleC.H"
#include "IndexDefines.H"

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
bool
is_inside(
  const int i,
  const int j,
  const int k,
  const amrex::Dim3 lo,
  const amrex::Dim3 hi,
  const int ng)
{
  return (i >= lo.x - ng && i <= hi.x + ng) &&
         (j >= lo.y - ng && j <= hi.y + ng) &&
         (k >= lo.z - ng && k <= hi.z + ng);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
bool
is_inside(
  const int i,
  const int j,
  const int k,
  const amrex::Dim3 lo,
  const amrex::Dim3 hi)
{
  return is_inside(i, j, k, lo, hi, 0);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
int
argmax(const amrex::Real array[AMREX_SPACEDIM])
{
  const amrex::Real max_num =
    amrex::max<amrex::Real>(AMREX_D_DECL(array[0], array[1], array[2]));
  for (int i = 0; i < AMREX_SPACEDIM; i++) {
    if (max_num == array[i])
      return i;
  }
  return 0;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
idxsort(const amrex::Real array[AMREX_SPACEDIM], int idx[AMREX_SPACEDIM])
{
  amrex::Real absarray[AMREX_SPACEDIM] = {0.0};
  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    absarray[dir] = amrex::Math::abs(array[dir]);
  }

  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    const int i = argmax(absarray);
    idx[dir] = i;
    absarray[i] = -1.0;
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_set_body_state(
  const int i,
  const int j,
  const int k,
  const int n,
  const amrex::Array4<const int>& mask,
  const amrex::GpuArray<amrex::Real, NVAR> body_state,
  const int covered_val,
  const amrex::Array4<amrex::Real>& S)
{
  if (mask(i, j, k) == covered_val) {
    S(i, j, k, n) = body_state[n];
  }
}

void pc_fill_sv_ebg(
  const amrex::Box&,
  const int,
  const amrex::Array4<const amrex::Real>&,
  const amrex::Array4<const amrex::Real>&,
  const amrex::Array4<const amrex::Real>&,
  const amrex::Array4<const amrex::Real>&,
  const amrex::Array4<const amrex::Real>&,
  EBBndryGeom*);

void pc_fill_bndry_grad_stencil(
  const amrex::Box&,
  const amrex::Real,
  const int,
  const EBBndryGeom*,
  const int,
  EBBndrySten*);

void pc_fill_flux_interp_stencil(
  const amrex::Box&,
  const amrex::Box,
  const int,
  const amrex::Array4<const amrex::Real>&,
  const amrex::Array4<const amrex::Real>&,
  FaceSten*);

void pc_apply_face_stencil(
  const amrex::Box&,
  const amrex::Box,
  const FaceSten*,
  const int,
  const int,
  const int,
  const amrex::Array4<amrex::Real>&);

void pc_eb_div(
  const amrex::Box&,
  const amrex::Real,
  const int,
  const EBBndryGeom*,
  const int,
  const amrex::Array4<const amrex::Real>&,
  const amrex::Array4<const amrex::Real>&,
  const amrex::Array4<const amrex::Real>&,
  const amrex::Real*,
  const amrex::Array4<const amrex::Real>&,
  const amrex::Array4<amrex::Real>&);

void pc_apply_eb_boundry_visc_flux_stencil(
  const amrex::Box&,
  const EBBndrySten*,
  const int,
  const EBBndryGeom*,
  const int,
  amrex::Array4<const amrex::Real> const&,
  amrex::Array4<const amrex::Real> const&,
  const amrex::Real*,
  const int,
  amrex::Real*,
  const int);

void pc_apply_eb_boundry_flux_stencil(
  const amrex::Box&,
  const EBBndrySten*,
  const int,
  amrex::Array4<const amrex::Real> const&,
  const int,
  amrex::Array4<const amrex::Real> const&,
  const int,
  const amrex::Real*,
  const int,
  amrex::Real*,
  const int,
  const int);

#endif
