#ifndef _MOL_H_
#define _MOL_H_

#ifdef PELEC_USE_EB
#include <AMReX_EBCellFlag.H>
#endif

#include "Diffterm.H"
#include "IndexDefines.H"
#include "PeleC.H"
#include "Riemann.H"
#include "PelePhysics.H"

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
mol_slope(
  const int i,
  const int j,
  const int k,
  const amrex::GpuArray<const int, 3> bdim,
  const amrex::GpuArray<const int, 3> q_idx,
  const amrex::Array4<const amrex::Real>& q,
  const amrex::Array4<const amrex::Real>& qaux,
  const amrex::Array4<amrex::Real>& dq
#ifdef PELEC_USE_EB
  ,
  const amrex::Array4<amrex::EBCellFlag const>& flags
#endif
)
{
#ifdef PELEC_USE_EB
  const bool flagArrayL =
    flags(i, j, k).isConnected(-bdim[0], -bdim[1], -bdim[2]) &&
    !flags(i, j, k).isCovered();
  const bool flagArrayR =
    flags(i, j, k).isConnected(+bdim[0], +bdim[1], +bdim[2]) &&
    !flags(i, j, k).isCovered();
#else
  const bool flagArrayL = true;
  const bool flagArrayR = true;
#endif

  amrex::Real dlft[QVAR] = {0.0};
  amrex::Real drgt[QVAR] = {0.0};
  const int im = i - bdim[0];
  const int jm = j - bdim[1];
  const int km = k - bdim[2];
  const int ip = i + bdim[0];
  const int jp = j + bdim[1];
  const int kp = k + bdim[2];

  dlft[0] = flagArrayL ? 0.5 * (q(i, j, k, QPRES) - q(im, jm, km, QPRES)) /
                             qaux(i, j, k, QC) -
                           0.5 * q(i, j, k, QRHO) *
                             (q(i, j, k, q_idx[0]) - q(im, jm, km, q_idx[0]))
                       : 0.0;
  dlft[1] = flagArrayL ? 0.5 * (q(i, j, k, QPRES) - q(im, jm, km, QPRES)) /
                             qaux(i, j, k, QC) +
                           0.5 * q(i, j, k, QRHO) *
                             (q(i, j, k, q_idx[0]) - q(im, jm, km, q_idx[0]))
                       : 0.0;
  dlft[2] = flagArrayL ? q(i, j, k, q_idx[1]) - q(im, jm, km, q_idx[1]) : 0.0;
  dlft[3] = flagArrayL ? q(i, j, k, q_idx[2]) - q(im, jm, km, q_idx[2]) : 0.0;

  drgt[0] = flagArrayR ? 0.5 * (q(ip, jp, kp, QPRES) - q(i, j, k, QPRES)) /
                             qaux(i, j, k, QC) -
                           0.5 * q(i, j, k, QRHO) *
                             (q(ip, jp, kp, q_idx[0]) - q(i, j, k, q_idx[0]))
                       : 0.0;
  drgt[1] = flagArrayR ? 0.5 * (q(ip, jp, kp, QPRES) - q(i, j, k, QPRES)) /
                             qaux(i, j, k, QC) +
                           0.5 * q(i, j, k, QRHO) *
                             (q(ip, jp, kp, q_idx[0]) - q(i, j, k, q_idx[0]))
                       : 0.0;
  drgt[2] = flagArrayR ? q(ip, jp, kp, q_idx[1]) - q(i, j, k, q_idx[1]) : 0.0;
  drgt[3] = flagArrayR ? q(ip, jp, kp, q_idx[2]) - q(i, j, k, q_idx[2]) : 0.0;

  for (int n = 0; n < NUM_SPECIES; n++) {
    dlft[4 + n] = flagArrayL ? q(i, j, k, QRHO) * q(i, j, k, QFS + n) -
                                 q(im, jm, km, QRHO) * q(im, jm, km, QFS + n) -
                                 q(i, j, k, QFS + n) *
                                   (q(i, j, k, QPRES) - q(im, jm, km, QPRES)) /
                                   (qaux(i, j, k, QC) * qaux(i, j, k, QC))
                             : 0.0;
    drgt[4 + n] = flagArrayR ? q(ip, jp, kp, QRHO) * q(ip, jp, kp, QFS + n) -
                                 q(i, j, k, QRHO) * q(i, j, k, QFS + n) -
                                 q(i, j, k, QFS + n) *
                                   (q(ip, jp, kp, QPRES) - q(i, j, k, QPRES)) /
                                   (qaux(i, j, k, QC) * qaux(i, j, k, QC))
                             : 0.0;
  }

  for (int n = 0; n < QVAR; n++) {
    const amrex::Real dcen = 0.5 * (dlft[n] + drgt[n]);
    const amrex::Real dlim =
      dlft[n] * drgt[n] >= 0.0
        ? 2.0 * amrex::min<amrex::Real>(
                  amrex::Math::abs(dlft[n]), amrex::Math::abs(drgt[n]))
        : 0.0;
    dq(i, j, k, n) = amrex::Math::copysign(1.0, dcen) *
                     amrex::min<amrex::Real>(dlim, amrex::Math::abs(dcen));
  }
}

void pc_compute_hyp_mol_flux(
  const amrex::Box& cbox,
  const amrex::Array4<const amrex::Real>& q,
  const amrex::Array4<const amrex::Real>& qaux,
  const amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM> flx,
  const amrex::GpuArray<const amrex::Array4<const amrex::Real>, AMREX_SPACEDIM>
    area,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> del,
  const int plm_iorder
#ifdef PELEC_USE_EB
  ,
  const amrex::Real eb_small_vfrac,
  const amrex::Array4<const amrex::Real>& vfrac,
  const amrex::Array4<amrex::EBCellFlag const>& flags,
  const EBBndryGeom* ebg,
  const int Nebg,
  amrex::Real* ebflux,
  const int nebflux
#endif
);

#endif
