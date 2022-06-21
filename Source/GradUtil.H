#ifndef GRADUTIL_H
#define GRADUTIL_H

#include <AMReX_FArrayBox.H>
#include "IndexDefines.H"
#include "EB.H"

// This header file contains functions and declarations for Gradient Utilities
// in 3D for PeleC GPU. As per the convention of AMReX, inlined functions are
// defined here. Where as non-inline functions are declared here.

namespace GradUtils {
const int nCompTan = AMREX_D_PICK(0, 2, 6);
} // namespace GradUtils

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_compute_tangential_vel_derivs(
  const int i,
  const int j,
  const int k,
  const amrex::Array4<const amrex::Real>& q,
  const int dir,
  const amrex::Real dx1,
  const amrex::Real dx2,
  const amrex::Array4<amrex::Real>& td)
{
  // dx1 and dx2 will be the trangential grid spacing
  const amrex::Real dx1inv = 1.0 / dx1;
  const amrex::Real dx2inv = 1.0 / dx2;
  // This is because it is a tangential velocity derivative.
  if (dir == 0) {
    // dx1 = dy, dx2 = dz
#if AMREX_SPACEDIM > 1
    td(i, j, k, 0) = 0.25 * dx1inv *
                     (q(i, j + 1, k, QU) + q(i - 1, j + 1, k, QU) -
                      q(i, j - 1, k, QU) - q(i - 1, j - 1, k, QU)); // dudy
    td(i, j, k, 1) = 0.25 * dx1inv *
                     (q(i, j + 1, k, QV) + q(i - 1, j + 1, k, QV) -
                      q(i, j - 1, k, QV) - q(i - 1, j - 1, k, QV)); // dvdy
#endif
#if AMREX_SPACEDIM > 2
    td(i, j, k, 2) = 0.25 * dx1inv *
                     (q(i, j + 1, k, QW) + q(i - 1, j + 1, k, QW) -
                      q(i, j - 1, k, QW) - q(i - 1, j - 1, k, QW)); // dwdy
    td(i, j, k, 3) = 0.25 * dx2inv *
                     (q(i, j, k + 1, QU) + q(i - 1, j, k + 1, QU) -
                      q(i, j, k - 1, QU) - q(i - 1, j, k - 1, QU)); // dudz
    td(i, j, k, 4) = 0.25 * dx2inv *
                     (q(i, j, k + 1, QV) + q(i - 1, j, k + 1, QV) -
                      q(i, j, k - 1, QV) - q(i - 1, j, k - 1, QV)); // dvdz
    td(i, j, k, 5) = 0.25 * dx2inv *
                     (q(i, j, k + 1, QW) + q(i - 1, j, k + 1, QW) -
                      q(i, j, k - 1, QW) - q(i - 1, j, k - 1, QW)); // dwdz
#endif
  } else if (dir == 1) {
    // dx1 = dx, dx2 = dz
    td(i, j, k, 0) = 0.25 * dx1inv *
                     (q(i + 1, j, k, QU) + q(i + 1, j - 1, k, QU) -
                      q(i - 1, j, k, QU) - q(i - 1, j - 1, k, QU)); // dudx
    td(i, j, k, 1) = 0.25 * dx1inv *
                     (q(i + 1, j, k, QV) + q(i + 1, j - 1, k, QV) -
                      q(i - 1, j, k, QV) - q(i - 1, j - 1, k, QV)); // dvdx
#if AMREX_SPACEDIM > 2
    td(i, j, k, 2) = 0.25 * dx1inv *
                     (q(i + 1, j, k, QW) + q(i + 1, j - 1, k, QW) -
                      q(i - 1, j, k, QW) - q(i - 1, j - 1, k, QW)); // dwdx
    td(i, j, k, 3) = 0.25 * dx2inv *
                     (q(i, j, k + 1, QU) + q(i, j - 1, k + 1, QU) -
                      q(i, j, k - 1, QU) - q(i, j - 1, k - 1, QU)); // dudz
    td(i, j, k, 4) = 0.25 * dx2inv *
                     (q(i, j, k + 1, QV) + q(i, j - 1, k + 1, QV) -
                      q(i, j, k - 1, QV) - q(i, j - 1, k - 1, QV)); // dvdz
    td(i, j, k, 5) = 0.25 * dx2inv *
                     (q(i, j, k + 1, QW) + q(i, j - 1, k + 1, QW) -
                      q(i, j, k - 1, QW) - q(i, j - 1, k - 1, QW)); // dwdz
#endif
  } else if (dir == 2) {
    // dx1 = dx, dx2 = dy
    td(i, j, k, 0) = 0.25 * dx1inv *
                     (q(i + 1, j, k, QU) + q(i + 1, j, k - 1, QU) -
                      q(i - 1, j, k, QU) - q(i - 1, j, k - 1, QU)); // dudx
    td(i, j, k, 1) = 0.25 * dx1inv *
                     (q(i + 1, j, k, QV) + q(i + 1, j, k - 1, QV) -
                      q(i - 1, j, k, QV) - q(i - 1, j, k - 1, QV)); // dvdx
    td(i, j, k, 2) = 0.25 * dx1inv *
                     (q(i + 1, j, k, QW) + q(i + 1, j, k - 1, QW) -
                      q(i - 1, j, k, QW) - q(i - 1, j, k - 1, QW)); // dwdx
    td(i, j, k, 3) = 0.25 * dx2inv *
                     (q(i, j + 1, k, QU) + q(i, j + 1, k - 1, QU) -
                      q(i, j - 1, k, QU) - q(i, j - 1, k - 1, QU)); // dudy
    td(i, j, k, 4) = 0.25 * dx2inv *
                     (q(i, j + 1, k, QV) + q(i, j + 1, k - 1, QV) -
                      q(i, j - 1, k, QV) - q(i, j - 1, k - 1, QV)); // dvdy
    td(i, j, k, 5) = 0.25 * dx2inv *
                     (q(i, j + 1, k, QW) + q(i, j + 1, k - 1, QW) -
                      q(i, j - 1, k, QW) - q(i, j - 1, k - 1, QW)); // dwdy
  }
}

void pc_compute_tangential_vel_derivs_eb(
  const amrex::Box& /*bx*/,
  const int /*dir*/,
  const amrex::Real /*dx1*/,
  const amrex::Real /*dx2*/,
  const EBBndryGeom* /*sv_ebg*/,
  const int /*Ncut*/,
  const amrex::Array4<const amrex::Real>& /*q*/,
  const amrex::Array4<amrex::EBCellFlag const>& /*flags*/,
  const amrex::Array4<amrex::Real>& /*td*/);
#endif
