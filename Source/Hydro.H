#ifndef _HYDRO_H_
#define _HYDRO_H_

#include "PeleC.H"
#include "IndexDefines.H"
#include "PelePhysics.H"
#include "Utilities.H"
#include "Godunov.H"

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_srctoprim(
  const int i,
  const int j,
  const int k,
  amrex::Array4<const amrex::Real> const& q,
  amrex::Array4<const amrex::Real> const& qa,
  amrex::Array4<const amrex::Real> const& src,
  amrex::Array4<amrex::Real> const& srcq,
  PassMap const& pmap)
{
  // Local loading Vars ==> Load once from global memory use multiple times
  const amrex::Real srcrho = src(i, j, k, URHO);
  const amrex::Real rhoinv = 1.0 / q(i, j, k, QRHO);
  const amrex::Real vx = q(i, j, k, QU);
  const amrex::Real vy = q(i, j, k, QV);
  const amrex::Real vz = q(i, j, k, QW);
  const amrex::Real umx = src(i, j, k, UMX);
  const amrex::Real umy = src(i, j, k, UMY);
  const amrex::Real umz = src(i, j, k, UMZ);
  // SrcQ Calc
  srcq(i, j, k, QRHO) = srcrho;
  srcq(i, j, k, QU) = (umx - vx * srcrho) * rhoinv;
  srcq(i, j, k, QV) = (umy - vy * srcrho) * rhoinv;
  srcq(i, j, k, QW) = (umz - vz * srcrho) * rhoinv;
  srcq(i, j, k, QREINT) = src(i, j, k, UEDEN) - vx * umx - vy * umy - vz * umz +
                          srcrho * (vx * vx + vy * vy + vz * vz);
  srcq(i, j, k, QPRES) =
    qa(i, j, k, QDPDE) *
      (srcq(i, j, k, QREINT) - q(i, j, k, QREINT) * srcrho * rhoinv) * rhoinv +
    qa(i, j, k, QDPDR) * srcrho;

  for (int ip = 0; ip < NPASSIVE; ++ip) {
    const int n = pmap.upassMap[ip];
    const int nq = pmap.qpassMap[ip];
    srcq(i, j, k, nq) = (src(i, j, k, n) - q(i, j, k, nq) * srcrho) * rhoinv;
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_divu(
  const int i,
  const int j,
  const int k,
  amrex::Array4<const amrex::Real> const& q,
  AMREX_D_DECL(
    const amrex::Real dx, const amrex::Real dy, const amrex::Real dz),
  amrex::Array4<amrex::Real> const& divu)
{
#if AMREX_SPACEDIM == 1
  amrex::Real ux = (q(i, j, k, QU) - q(i - 1, j, k, QU)) / dx;
#elif AMREX_SPACEDIM == 2
  amrex::Real ux = 0.5 *
                   (q(i, j, k, QU) - q(i - 1, j, k, QU) + q(i, j - 1, k, QU) -
                    q(i - 1, j - 1, k, QU)) /
                   dx;
  amrex::Real vy = 0.5 *
                   (q(i, j, k, QV) - q(i, j - 1, k, QV) + q(i - 1, j, k, QV) -
                    q(i - 1, j - 1, k, QV)) /
                   dy;
#elif AMREX_SPACEDIM == 3
  amrex::Real ux =
    0.25 *
    (q(i, j, k, QU) - q(i - 1, j, k, QU) + q(i, j, k - 1, QU) -
     q(i - 1, j, k - 1, QU) + q(i, j - 1, k, QU) - q(i - 1, j - 1, k, QU) +
     q(i, j - 1, k - 1, QU) - q(i - 1, j - 1, k - 1, QU)) /
    dx;
  amrex::Real vy =
    0.25 *
    (q(i, j, k, QV) - q(i, j - 1, k, QV) + q(i, j, k - 1, QV) -
     q(i, j - 1, k - 1, QV) + q(i - 1, j, k, QV) - q(i - 1, j - 1, k, QV) +
     q(i - 1, j, k - 1, QV) - q(i - 1, j - 1, k - 1, QV)) /
    dy;
  amrex::Real wz =
    0.25 *
    (q(i, j, k, QW) - q(i, j, k - 1, QW) + q(i - 1, j, k, QW) -
     q(i - 1, j, k - 1, QW) + q(i, j - 1, k, QW) - q(i, j - 1, k - 1, QW) +
     q(i - 1, j - 1, k, QW) - q(i - 1, j - 1, k - 1, QW)) /
    dz;
#endif
  divu(i, j, k) = AMREX_D_TERM(ux, +vy, +wz);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_norm_spec_flx(
  const int i, const int j, const int k, amrex::Array4<amrex::Real> const& flx)
{
  const int num = NUM_SPECIES + UFS;
  amrex::Real summ, fac;
  summ = 0.0;
  for (int n = UFS; n < num; ++n)
    summ += flx(i, j, k, n);

  fac = (summ != 0.0) ? flx(i, j, k, URHO) / summ : 1.0;
  for (int n = UFS; n < num; ++n)
    flx(i, j, k, n) *= fac;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_ext_flx(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::Real> const& flx,
  amrex::Array4<const amrex::Real> const& a)
{
  amrex::Real area = a(i, j, k); // avoid reading from global memory NVAR times
  for (int n = 0; n < NVAR; ++n)
    flx(i, j, k, n) *= area;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_update(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::Real> const& update,
  const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM> flx,
  amrex::Array4<const amrex::Real> const& vol,
  amrex::Array4<const amrex::Real> const& pdivu)
{
  amrex::Real v = vol(i, j, k);
  for (int n = 0; n < NVAR; ++n) {
    update(i, j, k, n) += (AMREX_D_TERM(
                            flx[0](i, j, k, n) - flx[0](i + 1, j, k, n),
                            +flx[1](i, j, k, n) - flx[1](i, j + 1, k, n),
                            +flx[2](i, j, k, n) - flx[2](i, j, k + 1, n))) /
                          v;
  }

  update(i, j, k, UEINT) -= pdivu(i, j, k);
}

// Host functions
void pc_umdrv(
  const int is_finest_level,
  const amrex::Real time,
  amrex::Box const& bx,
  const int* domlo,
  const int* domhi,
  const int* bclo,
  const int* bchi,
  amrex::Array4<const amrex::Real> const& uin,
  amrex::Array4<amrex::Real> const& uout,
  amrex::Array4<const amrex::Real> const& q,
  amrex::Array4<const amrex::Real> const& qaux,
  amrex::Array4<const amrex::Real> const&
    src_q, // amrex::IArrayBox const& bcMask,
  const amrex::Real* dx,
  const amrex::Real dt,
  const int ppm_type,
  const int use_flattening,
  const amrex::Real difmag,
  const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM> flx,
  const amrex::GpuArray<const amrex::Array4<const amrex::Real>, AMREX_SPACEDIM>
    a,
  amrex::Array4<amrex::Real> const& vol,
  amrex::Real cflLoc);

void pc_consup(
  amrex::Box const& bx,
  amrex::Array4<const amrex::Real> const& u,
  amrex::Array4<amrex::Real> const& update,
  const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM> flx,
  const amrex::GpuArray<const amrex::Array4<const amrex::Real>, AMREX_SPACEDIM>
    a,
  amrex::Array4<const amrex::Real> const& vol,
  amrex::Array4<const amrex::Real> const& div,
  amrex::Array4<const amrex::Real> const& pdivu,
  amrex::Real const* del,
  amrex::Real const difmag);
#endif
