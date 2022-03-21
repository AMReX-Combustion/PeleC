#ifndef _MOL_H_
#define _MOL_H_

#include <AMReX_EBCellFlag.H>

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
  const int dir,
  const amrex::GpuArray<const int, 3> q_idx,
  const amrex::Array4<const amrex::Real>& q,
  const amrex::Array4<const amrex::Real>& qaux,
  const amrex::Array4<amrex::Real>& dq,
  const amrex::Array4<amrex::EBCellFlag const>& flags)
{
  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  const amrex::IntVect dm(-amrex::IntVect::TheDimensionVector(dir));
  const amrex::IntVect dp(amrex::IntVect::TheDimensionVector(dir));
  const amrex::IntVect ivm(iv + dm);
  const amrex::IntVect ivp(iv + dp);

  const bool flagArrayL =
    flags(i, j, k).isConnected(dm) && !flags(i, j, k).isCovered();
  const bool flagArrayR =
    flags(i, j, k).isConnected(dp) && !flags(i, j, k).isCovered();

  amrex::Real dlft[QVAR] = {0.0};
  amrex::Real drgt[QVAR] = {0.0};

  dlft[0] = flagArrayL
              ? 0.5 * (q(iv, QPRES) - q(ivm, QPRES)) / qaux(iv, QC) -
                  0.5 * q(iv, QRHO) * (q(iv, q_idx[0]) - q(ivm, q_idx[0]))
              : 0.0;
  dlft[1] = flagArrayL
              ? 0.5 * (q(iv, QPRES) - q(ivm, QPRES)) / qaux(iv, QC) +
                  0.5 * q(iv, QRHO) * (q(iv, q_idx[0]) - q(ivm, q_idx[0]))
              : 0.0;
  dlft[2] = flagArrayL ? q(iv, q_idx[1]) - q(ivm, q_idx[1]) : 0.0;
  dlft[3] = flagArrayL ? q(iv, q_idx[2]) - q(ivm, q_idx[2]) : 0.0;

  drgt[0] = flagArrayR
              ? 0.5 * (q(ivp, QPRES) - q(iv, QPRES)) / qaux(iv, QC) -
                  0.5 * q(iv, QRHO) * (q(ivp, q_idx[0]) - q(iv, q_idx[0]))
              : 0.0;
  drgt[1] = flagArrayR
              ? 0.5 * (q(ivp, QPRES) - q(iv, QPRES)) / qaux(iv, QC) +
                  0.5 * q(iv, QRHO) * (q(ivp, q_idx[0]) - q(iv, q_idx[0]))
              : 0.0;
  drgt[2] = flagArrayR ? q(ivp, q_idx[1]) - q(iv, q_idx[1]) : 0.0;
  drgt[3] = flagArrayR ? q(ivp, q_idx[2]) - q(iv, q_idx[2]) : 0.0;

  for (int n = 0; n < NUM_SPECIES; n++) {
    dlft[4 + n] = flagArrayL
                    ? q(iv, QRHO) * q(iv, QFS + n) -
                        q(ivm, QRHO) * q(ivm, QFS + n) -
                        q(iv, QFS + n) * (q(iv, QPRES) - q(ivm, QPRES)) /
                          (qaux(iv, QC) * qaux(iv, QC))
                    : 0.0;
    drgt[4 + n] = flagArrayR
                    ? q(ivp, QRHO) * q(ivp, QFS + n) -
                        q(iv, QRHO) * q(iv, QFS + n) -
                        q(iv, QFS + n) * (q(ivp, QPRES) - q(iv, QPRES)) /
                          (qaux(iv, QC) * qaux(iv, QC))
                    : 0.0;
  }

  for (int n = 0; n < QVAR; n++) {
    const amrex::Real dcen = 0.5 * (dlft[n] + drgt[n]);
    const amrex::Real dlim =
      dlft[n] * drgt[n] >= 0.0
        ? 2.0 * amrex::min<amrex::Real>(
                  amrex::Math::abs(dlft[n]), amrex::Math::abs(drgt[n]))
        : 0.0;
    dq(iv, n) = amrex::Math::copysign(1.0, dcen) *
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
  const int plm_iorder,
  const bool use_laxf_flux,
  const amrex::Array4<amrex::EBCellFlag const>& flags,
  const EBBndryGeom* ebg,
  const int Nebg,
  amrex::Real* ebflux,
  const int nebflux);

#endif
