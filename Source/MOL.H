#ifndef MOL_H
#define MOL_H

#include <AMReX_EBCellFlag.H>

#include "IndexDefines.H"
#include "PeleC.H"
#include "Riemann.H"
#include "PelePhysics.H"

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_lax_cmpflx_passive(
  const amrex::Real ul,
  const amrex::Real ur,
  const amrex::Real rl,
  const amrex::Real rr,
  const amrex::Real ql,
  const amrex::Real qr,
  const amrex::Real maxeigval,
  amrex::Real& flux_tmp)
{
  flux_tmp =
    0.5 * ((rl * ul * ql + rr * ur * qr) - maxeigval * (rr * qr - rl * ql));
}

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
#if NUM_ADV > 0
  for (int n = 0; n < NUM_ADV; n++) {
    dlft[QFA + n] = flagArrayL ? q(iv, QFA + n) - q(ivm, QFA + n) : 0.0;
    drgt[QFA + n] = flagArrayR ? q(ivp, QFA + n) - q(iv, QFA + n) : 0.0;
  }
#endif
  for (int n = 0; n < NUM_SPECIES; n++) {
    dlft[QFS + n] = flagArrayL
                      ? q(iv, QRHO) * q(iv, QFS + n) -
                          q(ivm, QRHO) * q(ivm, QFS + n) -
                          q(iv, QFS + n) * (q(iv, QPRES) - q(ivm, QPRES)) /
                            (qaux(iv, QC) * qaux(iv, QC))
                      : 0.0;
    drgt[QFS + n] = flagArrayR
                      ? q(ivp, QRHO) * q(ivp, QFS + n) -
                          q(iv, QRHO) * q(iv, QFS + n) -
                          q(iv, QFS + n) * (q(ivp, QPRES) - q(iv, QPRES)) /
                            (qaux(iv, QC) * qaux(iv, QC))
                      : 0.0;
  }
#if NUM_AUX > 0
  for (int n = 0; n < NUM_AUX; n++) {
    dlft[QFX + n] = flagArrayL ? q(iv, QFX + n) - q(ivm, QFX + n) : 0.0;
    drgt[QFX + n] = flagArrayR ? q(ivp, QFX + n) - q(iv, QFX + n) : 0.0;
  }
#endif
#if NUM_LIN > 0
  for (int n = 0; n < NUM_LIN; n++) {
    dlft[QLIN + n] = flagArrayL ? q(iv, QLIN + n) - q(ivm, QLIN + n) : 0.0;
    drgt[QLIN + n] = flagArrayR ? q(ivp, QLIN + n) - q(iv, QLIN + n) : 0.0;
  }
#endif

  for (int n = 0; n < QVAR; n++) {
    const amrex::Real dcen = 0.5 * (dlft[n] + drgt[n]);
    const amrex::Real dlim =
      dlft[n] * drgt[n] >= 0.0
        ? 2.0 * amrex::min<amrex::Real>(std::abs(dlft[n]), std::abs(drgt[n]))
        : 0.0;
    dq(iv, n) =
      std::copysign(1.0, dcen) * amrex::min<amrex::Real>(dlim, std::abs(dcen));
  }
}

void pc_compute_hyp_mol_flux(
  const amrex::Box& cbox,
  const amrex::Array4<const amrex::Real>& q,
  const amrex::Array4<const amrex::Real>& qaux,
  const amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx,
  const amrex::GpuArray<const amrex::Array4<const amrex::Real>, AMREX_SPACEDIM>&
    area,
  const int mol_iorder,
  const bool use_laxf_flux,
  const amrex::Array4<amrex::EBCellFlag const>& flags);

void pc_compute_hyp_mol_flux_eb(
  amrex::Geometry const& geom,
  const amrex::Box& cbox,
  const amrex::Array4<const amrex::Real>& q,
  const amrex::Array4<const amrex::Real>& qaux,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx,
  const bool use_laxf_flux,
  const amrex::Array4<const amrex::Real>& vfrac,
  const EBBndryGeom* ebg,
  const int Nebg,
  amrex::Real* ebflux,
  const int nebflux);

#endif
