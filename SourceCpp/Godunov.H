#ifndef _GODUNOV_H_
#define _GODUNOV_H_

#include <cmath>

#include <AMReX_FArrayBox.H>
#include <AMReX_IArrayBox.H>

#include "Constants.H"
#include "IndexDefines.H"
#include "Riemann.H"

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
flatten(
  const int i,
  const int j,
  const int k,
  const int dir,
  amrex::Array4<const amrex::Real> const& q)
{
  const int bdim[3] = {dir == 0, dir == 1, dir == 2};
  const int n = QPRES;
  // Parameters from uflatten
  const amrex::Real small_pres = 1.e-200;
  const amrex::Real shktst = 0.33e0;
  const amrex::Real zcut1 = 0.75e0;
  const amrex::Real zcut2 = 0.85e0;
  const amrex::Real dzcut = 1.0 / (zcut2 - zcut1);

  // amrex::Real dlft, drgt, dcen, dfm, dfp, dlim, dsgn, dtemp, drgtm1, dlftp1,
  // slop; const amrex::Real qc = q(i, j, k, n);
  const amrex::Real qm =
    q(i - 1 * bdim[0], j - 1 * bdim[1], k - 1 * bdim[2], n);
  const amrex::Real qp =
    q(i + 1 * bdim[0], j + 1 * bdim[1], k + 1 * bdim[2], n);
  // assumes QU,QV,QW ordering
  const amrex::Real qmun =
    q(i - 1 * bdim[0], j - 1 * bdim[1], k - 1 * bdim[2], QU + dir);
  const amrex::Real qpun =
    q(i + 1 * bdim[0], j + 1 * bdim[1], k + 1 * bdim[2], QU + dir);
  const amrex::Real qmm =
    q(i - 2 * bdim[0], j - 2 * bdim[1], k - 2 * bdim[2], n);
  const amrex::Real qpp =
    q(i + 2 * bdim[0], j + 2 * bdim[1], k + 2 * bdim[2], n);

  amrex::Real dp = qp - qm;

  const int ishft = (dp > 0.0) ? 1 : -1;

  amrex::Real denom = amrex::max(small_pres, amrex::Math::abs(qpp - qmm));
  amrex::Real zeta = amrex::Math::abs(dp) / denom;
  const amrex::Real z =
    amrex::min(1.0, amrex::max(0.0, dzcut * (zeta - zcut1)));

  amrex::Real tst = (qmun - qpun >= 0.0) ? 1.0 : 0.0;

  amrex::Real tmp = amrex::min(qp, qm);

  const amrex::Real chi = (amrex::Math::abs(dp) / tmp) > shktst ? tst : 0.0;

  const amrex::Real qpshft =
    q(i + (1 - ishft) * bdim[0], j + (1 - ishft) * bdim[1],
      k + (1 - ishft) * bdim[2], n);
  const amrex::Real qmshft =
    q(i - (1 + ishft) * bdim[0], j - (1 + ishft) * bdim[1],
      k - (1 + ishft) * bdim[2], n);
  const amrex::Real qpunshft =
    q(i + (1 - ishft) * bdim[0], j + (1 - ishft) * bdim[1],
      k + (1 - ishft) * bdim[2], QU + dir);
  const amrex::Real qmunshft =
    q(i - (1 + ishft) * bdim[0], j - (1 + ishft) * bdim[1],
      k - (1 + ishft) * bdim[2], QU + dir);
  const amrex::Real qppshft =
    q(i + (2 - ishft) * bdim[0], j + (2 - ishft) * bdim[1],
      k + (2 - ishft) * bdim[2], n);
  const amrex::Real qmmshft =
    q(i - (2 + ishft) * bdim[0], j - (2 + ishft) * bdim[1],
      k - (2 + ishft) * bdim[2], n);

  dp = qpshft - qmshft;

  denom = amrex::max(small_pres, amrex::Math::abs(qppshft - qmmshft));
  zeta = amrex::Math::abs(dp) / denom;
  const amrex::Real z2 =
    amrex::min(1.0, amrex::max(0.0, dzcut * (zeta - zcut1)));

  tst = (qmunshft - qpunshft >= 0.0) ? 1.0 : 0.0;

  tmp = amrex::min(qpshft, qmshft);

  const amrex::Real chi2 = (amrex::Math::abs(dp) / tmp) > shktst ? tst : 0.0;

  return 1.0 - amrex::max(chi2 * z2, chi * z);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_cmpflx(
  const int i,
  const int j,
  const int k,
  const int bclo,
  const int bchi,
  const int domlo,
  const int domhi,
  amrex::Array4<const amrex::Real> const& ql,
  amrex::Array4<const amrex::Real> const& qr,
  amrex::Array4<amrex::Real> const& flx,
  amrex::Array4<amrex::Real> const& q,
  amrex::Array4<const amrex::Real> const& qa,
  // amrex::Array4<const int> const& bcMask,
  const int dir)
{
  amrex::Real cav, ustar;
  amrex::Real spl[NUM_SPECIES];
  amrex::Real spr[NUM_SPECIES];
  amrex::Real ul, ur, vl, vr, v2l, v2r, rel, rer, gamcl, gamcr;
  int idx;
  int IU, IV, IV2;
  int GU, GV, GV2;
  int f_idx[3];
  if (dir == 0) {
    IU = QU;
    IV = QV;
    IV2 = QW;
    GU = GDU;
    GV = GDV;
    GV2 = GDW;
    gamcl = qa(i - 1, j, k, QGAMC);
    gamcr = qa(i, j, k, QGAMC);
    cav = 0.5 * (qa(i, j, k, QC) + qa(i - 1, j, k, QC));
    f_idx[0] = UMX;
    f_idx[1] = UMY;
    f_idx[2] = UMZ;
  } else if (dir == 1) {
    IU = QV;
    IV = QU;
    IV2 = QW;
    GU = GDV;
    GV = GDU;
    GV2 = GDW;
    gamcl = qa(i, j - 1, k, QGAMC);
    gamcr = qa(i, j, k, QGAMC);
    cav = 0.5 * (qa(i, j, k, QC) + qa(i, j - 1, k, QC));
    f_idx[0] = UMY;
    f_idx[1] = UMX;
    f_idx[2] = UMZ;
  } else {
    IU = QW;
    IV = QU;
    IV2 = QV;
    GU = GDW;
    GV = GDU;
    GV2 = GDV;
    gamcl = qa(i, j, k - 1, QGAMC);
    gamcr = qa(i, j, k, QGAMC);
    cav = 0.5 * (qa(i, j, k, QC) + qa(i, j, k - 1, QC));
    f_idx[0] = UMZ;
    f_idx[1] = UMX;
    f_idx[2] = UMY;
  }

  for (int sp = 0; sp < NUM_SPECIES; ++sp) {
    spl[sp] = ql(i, j, k, QFS + sp);
    spr[sp] = qr(i, j, k, QFS + sp);
  }
  ul = ql(i, j, k, IU);
  vl = ql(i, j, k, IV);
  v2l = ql(i, j, k, IV2);
  rel = ql(i, j, k, QREINT);

  ur = qr(i, j, k, IU);
  vr = qr(i, j, k, IV);
  v2r = qr(i, j, k, IV2);
  rer = qr(i, j, k, QREINT);

  // Outflow Hack
  if (dir == 2)
    idx = k;
  else
    idx = (dir == 0) ? i : j;
  if (bclo == Outflow && idx == domlo) {
    ul = ur;
    vl = vr;
    v2l = v2r;
    rel = rer;
  }
  if (bchi == Outflow && idx == domhi + 1) {
    ur = ul;
    vr = vl;
    v2r = v2l;
    rer = rel;
  }

  const int bc_test_val = 1;
  riemann(
    ql(i, j, k, QRHO), ul, vl, v2l, ql(i, j, k, QPRES), rel, spl, gamcl,
    qr(i, j, k, QRHO), ur, vr, v2r, qr(i, j, k, QPRES), rer, spr, gamcr,
    bc_test_val, qa(i, j, k, QCSML), cav, ustar, flx(i, j, k, URHO),
    flx(i, j, k, f_idx[0]), flx(i, j, k, f_idx[1]), flx(i, j, k, f_idx[2]),
    flx(i, j, k, UEDEN), flx(i, j, k, UEINT), q(i, j, k, GU), q(i, j, k, GV),
    q(i, j, k, GV2), q(i, j, k, GDPRES), q(i, j, k, GDGAME));

  amrex::Real flxrho = flx(i, j, k, URHO);
  for (int ipass = 0; ipass < NPASSIVE; ++ipass) {
    int n = indxmap::upass_map(ipass);
    int nqp = indxmap::qpass_map(ipass);
    if (ustar > 0.0)
      flx(i, j, k, n) = flxrho * ql(i, j, k, nqp);
    else if (ustar < 0.0)
      flx(i, j, k, n) = flxrho * qr(i, j, k, nqp);
    else
      flx(i, j, k, n) = flxrho * 0.50 * (ql(i, j, k, nqp) + qr(i, j, k, nqp));
  }
}

// First Transverse Correction for Predicted Y-states, using X-Flux
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_transx1(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::Real> const& qm,
  amrex::Array4<amrex::Real> const& qp,
  amrex::Array4<const amrex::Real> const& qym,
  amrex::Array4<const amrex::Real> const& qyp,
  amrex::Array4<const amrex::Real> const& flxx,
  amrex::Array4<const amrex::Real> const& qa,
  amrex::Array4<const amrex::Real> const& qint,
  const amrex::Real cdtdx)
{
  // Local variables to avoid reading from global memory too much
  amrex::Real compn, rr, rrnew, compo, rwr, c;
  amrex::Real rrr, rur, rvr, ekinr, rer, rrnewr;
  amrex::Real runewr, rvnewr, rwnewr, renewr, rhoekinr;
  amrex::Real rrl, rul, rvl, ekinl, rel, rrnewl, runewl, rvnewl, rwnewl, renewl;
  amrex::Real rhoekinl, rwl;
  amrex::Real pggp, pggm, ugp, ugm;
  amrex::Real dAup, pav, dAu;
  amrex::Real flxrho, flxu, flxv, flxw, flxe;

  flxrho = cdtdx * (flxx(i + 1, j, k, URHO) - flxx(i, j, k, URHO));
  flxu = cdtdx * (flxx(i + 1, j, k, UMX) - flxx(i, j, k, UMX));
  flxv = cdtdx * (flxx(i + 1, j, k, UMY) - flxx(i, j, k, UMY));
  flxw = cdtdx * (flxx(i + 1, j, k, UMZ) - flxx(i, j, k, UMZ));
  flxe = cdtdx * (flxx(i + 1, j, k, UEDEN) - flxx(i, j, k, UEDEN));
  c = qa(i, j, k, QGAMC);
  // Update passive variables
  for (int ipass = 0; ipass < NPASSIVE; ++ipass) {
    int n = indxmap::upass_map(ipass);
    int nqp = indxmap::qpass_map(ipass);
    compn = cdtdx * (flxx(i + 1, j, k, n) - flxx(i, j, k, n));
    rr = qyp(i, j, k, QRHO);
    rrnew = rr - flxrho;
    compo = rr * qyp(i, j, k, nqp) - compn;
    qp(i, j, k, nqp) = compo / rrnew;

    rr = qym(i, j + 1, k, QRHO);
    rrnew = rr - flxrho;
    compo = rr * qym(i, j + 1, k, nqp) - compn;
    qm(i, j + 1, k, nqp) = compo / rrnew;
  }
  // Update hydro vars

  pggp = qint(i + 1, j, k, GDPRES);
  pggm = qint(i, j, k, GDPRES);
  ugp = qint(i + 1, j, k, GDU);
  ugm = qint(i, j, k, GDU);

  dAup = pggp * ugp - pggm * ugm;
  pav = 0.5 * (pggp + pggm);
  dAu = ugp - ugm;

  // QP

  // Convert to conservative
  rrr = qyp(i, j, k, QRHO);
  rur = qyp(i, j, k, QU);
  rvr = qyp(i, j, k, QV);
  rwr = qyp(i, j, k, QW);
  ekinr = 0.5 * rrr * (rur * rur + rvr * rvr + rwr * rwr);
  rur *= rrr;
  rvr *= rrr;
  rwr *= rrr;

  rer = qyp(i, j, k, QREINT) + ekinr;
  // Add transverse predictor
  rrnewr = rrr - flxrho;
  runewr = rur - flxu;
  rvnewr = rvr - flxv;
  rwnewr = rwr - flxw;
  renewr = rer - flxe;

  // Convert back to primitive
  qp(i, j, k, QRHO) = rrnewr;
  qp(i, j, k, QU) = runewr / rrnewr;
  qp(i, j, k, QV) = rvnewr / rrnewr;
  qp(i, j, k, QW) = rwnewr / rrnewr;
  rhoekinr =
    0.5 * (runewr * runewr + rvnewr * rvnewr + rwnewr * rwnewr) / rrnewr;
  amrex::Real pnewr =
    qyp(i, j, k, QPRES) - cdtdx * (dAup + pav * dAu * (c - 1.0));
  qp(i, j, k, QPRES) = amrex::max(pnewr, SMALL_PRES);
  qp(i, j, k, QREINT) = renewr - rhoekinr;
  // QM

  // Conversion to Conservative
  rrl = qym(i, j + 1, k, QRHO);
  rul = qym(i, j + 1, k, QU);
  rvl = qym(i, j + 1, k, QV);
  rwl = qym(i, j + 1, k, QW);
  ekinl = 0.5 * rrl * (rul * rul + rvl * rvl + rwl * rwl);
  rul *= rrl;
  rvl *= rrl;
  rwl *= rrl;
  rel = qym(i, j + 1, k, QREINT) + ekinl;

  // Transverse fluxes
  rrnewl = rrl - flxrho;
  runewl = rul - flxu;
  rvnewl = rvl - flxv;
  rwnewl = rwl - flxw;
  renewl = rel - flxe;

  qm(i, j + 1, k, QRHO) = rrnewl;
  qm(i, j + 1, k, QU) = runewl / rrnewl;
  qm(i, j + 1, k, QV) = rvnewl / rrnewl;
  qm(i, j + 1, k, QW) = rwnewl / rrnewl;
  rhoekinl =
    0.5 * (runewl * runewl + rvnewl * rvnewl + rwnewl * rwnewl) / rrnewl;

  qm(i, j + 1, k, QREINT) = renewl - rhoekinl;
  amrex::Real pnewl =
    qym(i, j + 1, k, QPRES) - cdtdx * (dAup + pav * dAu * (c - 1.0));
  qm(i, j + 1, k, QPRES) = amrex::max(pnewl, SMALL_PRES);
}

// Corrects the left and right x states with y flux
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_transy1(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::Real> const& qm,
  amrex::Array4<amrex::Real> const& qp,
  amrex::Array4<const amrex::Real> const& qxm,
  amrex::Array4<const amrex::Real> const& qxp,
  amrex::Array4<const amrex::Real> const& flxy,
  amrex::Array4<const amrex::Real> const& qa,
  amrex::Array4<const amrex::Real> const& qint,
  const amrex::Real cdtdy)
{

  amrex::Real compn, rr, rrnew, compo, rwr, c;
  amrex::Real rrr, rur, rvr, ekinr, rer, rrnewr;
  amrex::Real runewr, rvnewr, rwnewr, renewr, rhoekinr;
  amrex::Real rrl, rul, rvl, ekinl, rel, rrnewl, runewl, rvnewl, rwnewl, renewl;
  amrex::Real rhoekinl, rwl;
  amrex::Real pggp, pggm, ugp, ugm;
  // Local variables to avoid reading from global memory too much
  amrex::Real dup, pav, du, pnewr, pnewl;
  amrex::Real flxrho, flxu, flxv, flxw, flxe;

  flxrho = cdtdy * (flxy(i, j + 1, k, URHO) - flxy(i, j, k, URHO));
  flxu = cdtdy * (flxy(i, j + 1, k, UMX) - flxy(i, j, k, UMX));
  flxv = cdtdy * (flxy(i, j + 1, k, UMY) - flxy(i, j, k, UMY));
  flxw = cdtdy * (flxy(i, j + 1, k, UMZ) - flxy(i, j, k, UMZ));
  flxe = cdtdy * (flxy(i, j + 1, k, UEDEN) - flxy(i, j, k, UEDEN));
  c = qa(i, j, k, QGAMC);
  // -  Update Passive Variables
  for (int ipass = 0; ipass < NPASSIVE; ++ipass) {
    int n = indxmap::upass_map(ipass);
    int nqp = indxmap::qpass_map(ipass);
    compn = cdtdy * (flxy(i, j + 1, k, n) - flxy(i, j, k, n));
    rr = qxp(i, j, k, QRHO);
    rrnew = rr - flxrho;
    compo = rr * qxp(i, j, k, nqp) - compn;
    qp(i, j, k, nqp) = compo / rrnew;

    rr = qxm(i + 1, j, k, QRHO);
    rrnew = rr - flxrho;
    compo = rr * qxm(i + 1, j, k, nqp) - compn;
    qm(i + 1, j, k, nqp) = compo / rrnew;
  }

  //- Update hydro vars -
  // pressure divergence and such.
  pggp = qint(i, j + 1, k, GDPRES);
  pggm = qint(i, j, k, GDPRES);
  ugp = qint(i, j + 1, k, GDV);
  ugm = qint(i, j, k, GDV);

  dup = pggp * ugp - pggm * ugm;
  pav = 0.5 * (pggp + pggm);
  du = ugp - ugm;

  // QP

  // Convert to conservative
  rrr = qxp(i, j, k, QRHO);
  rur = qxp(i, j, k, QU);
  rvr = qxp(i, j, k, QV);
  rwr = qxp(i, j, k, QW);
  ekinr = 0.5 * rrr * (rur * rur + rvr * rvr + rwr * rwr);
  rur *= rrr;
  rvr *= rrr;
  rwr *= rrr;

  rer = qxp(i, j, k, QREINT) + ekinr;
  // Add transverse predictor
  rrnewr = rrr - flxrho;
  runewr = rur - flxu;
  rvnewr = rvr - flxv;
  rwnewr = rwr - flxw;
  renewr = rer - flxe;

  // Convert back to primitive
  qp(i, j, k, QRHO) = rrnewr; // + hdt*srcr;
  qp(i, j, k, QU) = runewr / rrnewr;
  qp(i, j, k, QV) = rvnewr / rrnewr;
  qp(i, j, k, QW) = rwnewr / rrnewr;
  rhoekinr =
    0.5 * (runewr * runewr + rvnewr * rvnewr + rwnewr * rwnewr) / rrnewr;

  pnewr = qxp(i, j, k, QPRES) - cdtdy * (dup + pav * du * (c - 1.0));
  qp(i, j, k, QPRES) = amrex::max(pnewr, SMALL_PRES);
  qp(i, j, k, QREINT) = renewr - rhoekinr;

  // QM

  // Conversion to Conservative
  rrl = qxm(i + 1, j, k, QRHO);
  rul = qxm(i + 1, j, k, QU);
  rvl = qxm(i + 1, j, k, QV);
  rwl = qxm(i + 1, j, k, QW);
  ekinl = 0.5 * rrl * (rul * rul + rvl * rvl + rwl * rwl);
  rul *= rrl;
  rvl *= rrl;
  rwl *= rrl;
  rel = qxm(i + 1, j, k, QREINT) + ekinl;

  // Transverse fluxes
  rrnewl = rrl - flxrho;
  runewl = rul - flxu;
  rvnewl = rvl - flxv;
  rwnewl = rwl - flxw;
  renewl = rel - flxe;

  qm(i + 1, j, k, QRHO) = rrnewl;
  qm(i + 1, j, k, QU) = runewl / rrnewl;
  qm(i + 1, j, k, QV) = rvnewl / rrnewl;
  qm(i + 1, j, k, QW) = rwnewl / rrnewl;
  rhoekinl =
    0.5 * (runewl * runewl + rvnewl * rvnewl + rwnewl * rwnewl) / rrnewl;

  pnewl = qxm(i + 1, j, k, QPRES) - cdtdy * (dup + pav * du * (c - 1.0));
  qm(i + 1, j, k, QPRES) = amrex::max(pnewl, SMALL_PRES);
  qm(i + 1, j, k, QREINT) = renewl - rhoekinl;
}

/* Transverse Correction for Predicted Z-states, using X-Flux */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_transx2(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::Real> const& qm,
  amrex::Array4<amrex::Real> const& qp,
  amrex::Array4<const amrex::Real> const& qzm,
  amrex::Array4<const amrex::Real> const& qzp,
  amrex::Array4<const amrex::Real> const& flxx,
  amrex::Array4<const amrex::Real> const& qa,
  amrex::Array4<const amrex::Real> const& qint,
  const amrex::Real cdtdx)
{
  // Local variables to avoid reading from global memory too much
  amrex::Real compn, rr, rrnew, compo, rwr, c;
  amrex::Real rrr, rur, rvr, ekinr, rer, rrnewr;
  amrex::Real runewr, rvnewr, rwnewr, renewr, rhoekinr;
  amrex::Real rrl, rul, rvl, ekinl, rel, rrnewl, runewl, rvnewl, rwnewl, renewl;
  amrex::Real rhoekinl, rwl;
  amrex::Real pggp, pggm, ugp, ugm;
  amrex::Real dAup, pav, dAu;
  amrex::Real flxrho, flxu, flxv, flxw, flxe;

  flxrho = cdtdx * (flxx(i + 1, j, k, URHO) - flxx(i, j, k, URHO));
  flxu = cdtdx * (flxx(i + 1, j, k, UMX) - flxx(i, j, k, UMX));
  flxv = cdtdx * (flxx(i + 1, j, k, UMY) - flxx(i, j, k, UMY));
  flxw = cdtdx * (flxx(i + 1, j, k, UMZ) - flxx(i, j, k, UMZ));
  flxe = cdtdx * (flxx(i + 1, j, k, UEDEN) - flxx(i, j, k, UEDEN));
  c = qa(i, j, k, QGAMC);
  // -  Update Passive Variables
  for (int ipass = 0; ipass < NPASSIVE; ++ipass) {
    int n = indxmap::upass_map(ipass);
    int nqp = indxmap::qpass_map(ipass);
    compn = cdtdx * (flxx(i + 1, j, k, n) - flxx(i, j, k, n));
    rr = qzp(i, j, k, QRHO);
    rrnew = rr - flxrho;
    compo = rr * qzp(i, j, k, nqp) - compn;
    qp(i, j, k, nqp) = compo / rrnew;

    rr = qzm(i, j, k + 1, QRHO);
    rrnew = rr - flxrho;
    compo = rr * qzm(i, j, k + 1, nqp) - compn;
    qm(i, j, k + 1, nqp) = compo / rrnew;
  }
  //- Update hydro vars -

  pggp = qint(i + 1, j, k, GDPRES);
  pggm = qint(i, j, k, GDPRES);
  ugp = qint(i + 1, j, k, GDU);
  ugm = qint(i, j, k, GDU);

  dAup = pggp * ugp - pggm * ugm;
  pav = 0.5 * (pggp + pggm);
  dAu = ugp - ugm;

  // QP

  // Convert to conservative
  rrr = qzp(i, j, k, QRHO);
  rur = qzp(i, j, k, QU);
  rvr = qzp(i, j, k, QV);
  rwr = qzp(i, j, k, QW);
  ekinr = 0.5 * rrr * (rur * rur + rvr * rvr + rwr * rwr);
  rur *= rrr;
  rvr *= rrr;
  rwr *= rrr;

  rer = qzp(i, j, k, QREINT) + ekinr;
  // Add transverse predictor
  rrnewr = rrr - flxrho;
  runewr = rur - flxu;
  rvnewr = rvr - flxv;
  rwnewr = rwr - flxw;
  renewr = rer - flxe;

  // Convert back to primitive
  qp(i, j, k, QRHO) = rrnewr;
  qp(i, j, k, QU) = runewr / rrnewr;
  qp(i, j, k, QV) = rvnewr / rrnewr;
  qp(i, j, k, QW) = rwnewr / rrnewr;
  rhoekinr =
    0.5 * (runewr * runewr + rvnewr * rvnewr + rwnewr * rwnewr) / rrnewr;
  amrex::Real pnewr =
    qzp(i, j, k, QPRES) - cdtdx * (dAup + pav * dAu * (c - 1.0));
  qp(i, j, k, QPRES) = amrex::max(pnewr, SMALL_PRES);
  qp(i, j, k, QREINT) = renewr - rhoekinr;
  // QM

  // Conversion to Conservative
  rrl = qzm(i, j, k + 1, QRHO);
  rul = qzm(i, j, k + 1, QU);
  rvl = qzm(i, j, k + 1, QV);
  rwl = qzm(i, j, k + 1, QW);
  ekinl = 0.5 * rrl * (rul * rul + rvl * rvl + rwl * rwl);
  rul *= rrl;
  rvl *= rrl;
  rwl *= rrl;
  rel = qzm(i, j, k + 1, QREINT) + ekinl;

  // Transverse fluxes
  rrnewl = rrl - flxrho;
  runewl = rul - flxu;
  rvnewl = rvl - flxv;
  rwnewl = rwl - flxw;
  renewl = rel - flxe;

  qm(i, j, k + 1, QRHO) = rrnewl;
  qm(i, j, k + 1, QU) = runewl / rrnewl;
  qm(i, j, k + 1, QV) = rvnewl / rrnewl;
  qm(i, j, k + 1, QW) = rwnewl / rrnewl;
  rhoekinl =
    0.5 * (runewl * runewl + rvnewl * rvnewl + rwnewl * rwnewl) / rrnewl;

  qm(i, j, k + 1, QREINT) = renewl - rhoekinl;
  amrex::Real pnewl =
    qzm(i, j, k + 1, QPRES) - cdtdx * (dAup + pav * dAu * (c - 1.0));
  qm(i, j, k + 1, QPRES) = amrex::max(pnewl, SMALL_PRES);
}

// Corrects the left and right z states with y flux
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_transy2(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::Real> const& qm,
  amrex::Array4<amrex::Real> const& qp,
  amrex::Array4<const amrex::Real> const& qzm,
  amrex::Array4<const amrex::Real> const& qzp,
  amrex::Array4<const amrex::Real> const& flxy,
  amrex::Array4<const amrex::Real> const& qa,
  amrex::Array4<const amrex::Real> const& qint,
  const amrex::Real cdtdy)
{
  // Local variables to avoid reading from global memory too much
  amrex::Real compn, rr, rrnew, compo, rwr, c;
  amrex::Real rrr, rur, rvr, ekinr, rer, rrnewr;
  amrex::Real runewr, rvnewr, rwnewr, renewr, rhoekinr;
  amrex::Real rrl, rul, rvl, ekinl, rel, rrnewl, runewl, rvnewl, rwnewl, renewl;
  amrex::Real rhoekinl, rwl;
  amrex::Real pggp, pggm, ugp, ugm;
  amrex::Real dAup, pav, dAu;
  amrex::Real flxrho, flxu, flxv, flxw, flxe;

  flxrho = cdtdy * (flxy(i, j + 1, k, URHO) - flxy(i, j, k, URHO));
  flxu = cdtdy * (flxy(i, j + 1, k, UMX) - flxy(i, j, k, UMX));
  flxv = cdtdy * (flxy(i, j + 1, k, UMY) - flxy(i, j, k, UMY));
  flxw = cdtdy * (flxy(i, j + 1, k, UMZ) - flxy(i, j, k, UMZ));
  flxe = cdtdy * (flxy(i, j + 1, k, UEDEN) - flxy(i, j, k, UEDEN));
  c = qa(i, j, k, QGAMC);
  // -  Update Passive Variables
  for (int ipass = 0; ipass < NPASSIVE; ++ipass) {
    int n = indxmap::upass_map(ipass);
    int nqp = indxmap::qpass_map(ipass);
    compn = cdtdy * (flxy(i, j + 1, k, n) - flxy(i, j, k, n));
    rr = qzp(i, j, k, QRHO);
    rrnew = rr - flxrho;
    compo = rr * qzp(i, j, k, nqp) - compn;
    qp(i, j, k, nqp) = compo / rrnew;

    rr = qzm(i, j, k + 1, QRHO);
    rrnew = rr - flxrho;
    compo = rr * qzm(i, j, k + 1, nqp) - compn;
    qm(i, j, k + 1, nqp) = compo / rrnew;
  }
  //- Update hydro vars -

  pggp = qint(i, j + 1, k, GDPRES);
  pggm = qint(i, j, k, GDPRES);
  ugp = qint(i, j + 1, k, GDV);
  ugm = qint(i, j, k, GDV);

  dAup = pggp * ugp - pggm * ugm;
  pav = 0.5 * (pggp + pggm);
  dAu = ugp - ugm;

  // QP

  // Convert to conservative
  rrr = qzp(i, j, k, QRHO);
  rur = qzp(i, j, k, QU);
  rvr = qzp(i, j, k, QV);
  rwr = qzp(i, j, k, QW);
  ekinr = 0.5 * rrr * (rur * rur + rvr * rvr + rwr * rwr);
  rur *= rrr;
  rvr *= rrr;
  rwr *= rrr;

  rer = qzp(i, j, k, QREINT) + ekinr;
  // Add transverse predictor
  rrnewr = rrr - flxrho;
  runewr = rur - flxu;
  rvnewr = rvr - flxv;
  rwnewr = rwr - flxw;
  renewr = rer - flxe;

  // Convert back to primitive
  qp(i, j, k, QRHO) = rrnewr;
  qp(i, j, k, QU) = runewr / rrnewr;
  qp(i, j, k, QV) = rvnewr / rrnewr;
  qp(i, j, k, QW) = rwnewr / rrnewr;
  rhoekinr =
    0.5 * (runewr * runewr + rvnewr * rvnewr + rwnewr * rwnewr) / rrnewr;
  amrex::Real pnewr =
    qzp(i, j, k, QPRES) - cdtdy * (dAup + pav * dAu * (c - 1.0));
  qp(i, j, k, QPRES) = amrex::max(pnewr, SMALL_PRES);
  qp(i, j, k, QREINT) = renewr - rhoekinr;
  // QM

  // Conversion to Conservative
  rrl = qzm(i, j, k + 1, QRHO);
  rul = qzm(i, j, k + 1, QU);
  rvl = qzm(i, j, k + 1, QV);
  rwl = qzm(i, j, k + 1, QW);
  ekinl = 0.5 * rrl * (rul * rul + rvl * rvl + rwl * rwl);
  rul *= rrl;
  rvl *= rrl;
  rwl *= rrl;
  rel = qzm(i, j, k + 1, QREINT) + ekinl;

  // Transverse fluxes
  rrnewl = rrl - flxrho;
  runewl = rul - flxu;
  rvnewl = rvl - flxv;
  rwnewl = rwl - flxw;
  renewl = rel - flxe;

  qm(i, j, k + 1, QRHO) = rrnewl;
  qm(i, j, k + 1, QU) = runewl / rrnewl;
  qm(i, j, k + 1, QV) = rvnewl / rrnewl;
  qm(i, j, k + 1, QW) = rwnewl / rrnewl;
  rhoekinl =
    0.5 * (runewl * runewl + rvnewl * rvnewl + rwnewl * rwnewl) / rrnewl;

  qm(i, j, k + 1, QREINT) = renewl - rhoekinl;
  amrex::Real pnewl =
    qzm(i, j, k + 1, QPRES) - cdtdy * (dAup + pav * dAu * (c - 1.0));
  qm(i, j, k + 1, QPRES) = amrex::max(pnewl, SMALL_PRES);
}

// Corrects the left and right x states with z flux
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_transz1(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::Real> const& qm,
  amrex::Array4<amrex::Real> const& qp,
  amrex::Array4<const amrex::Real> const& qxm,
  amrex::Array4<const amrex::Real> const& qxp,
  amrex::Array4<const amrex::Real> const& flxz,
  amrex::Array4<const amrex::Real> const& qa,
  amrex::Array4<const amrex::Real> const& qint,
  const amrex::Real cdtdz)
{
  amrex::Real compn, rr, rrnew, compo, rwr, c;
  amrex::Real rrr, rur, rvr, ekinr, rer, rrnewr;
  amrex::Real runewr, rvnewr, rwnewr, renewr, rhoekinr;
  amrex::Real rrl, rul, rvl, ekinl, rel, rrnewl, runewl, rvnewl, rwnewl, renewl;
  amrex::Real rhoekinl, rwl;
  amrex::Real pggp, pggm, ugp, ugm;
  // Local variables to avoid reading from global memory too much
  amrex::Real dup, pav, du, pnewr, pnewl;
  amrex::Real flxrho, flxu, flxv, flxw, flxe;

  flxrho = cdtdz * (flxz(i, j, k + 1, URHO) - flxz(i, j, k, URHO));
  flxu = cdtdz * (flxz(i, j, k + 1, UMX) - flxz(i, j, k, UMX));
  flxv = cdtdz * (flxz(i, j, k + 1, UMY) - flxz(i, j, k, UMY));
  flxw = cdtdz * (flxz(i, j, k + 1, UMZ) - flxz(i, j, k, UMZ));
  flxe = cdtdz * (flxz(i, j, k + 1, UEDEN) - flxz(i, j, k, UEDEN));
  c = qa(i, j, k, QGAMC);
  // -  Update Passive Variables
  for (int ipass = 0; ipass < NPASSIVE; ++ipass) {
    int n = indxmap::upass_map(ipass);
    int nqp = indxmap::qpass_map(ipass);
    compn = cdtdz * (flxz(i, j, k + 1, n) - flxz(i, j, k, n));
    rr = qxp(i, j, k, QRHO);
    rrnew = rr - flxrho;
    compo = rr * qxp(i, j, k, nqp) - compn;
    qp(i, j, k, nqp) = compo / rrnew;

    rr = qxm(i + 1, j, k, QRHO);
    rrnew = rr - flxrho;
    compo = rr * qxm(i + 1, j, k, nqp) - compn;
    qm(i + 1, j, k, nqp) = compo / rrnew;
  }

  //- Update hydro vars -
  // pressure divergence and such.
  pggp = qint(i, j, k + 1, GDPRES);
  pggm = qint(i, j, k, GDPRES);
  ugp = qint(i, j, k + 1, GDW);
  ugm = qint(i, j, k, GDW);

  dup = pggp * ugp - pggm * ugm;
  pav = 0.5 * (pggp + pggm);
  du = ugp - ugm;

  // QP

  // Convert to conservative
  rrr = qxp(i, j, k, QRHO);
  rur = qxp(i, j, k, QU);
  rvr = qxp(i, j, k, QV);
  rwr = qxp(i, j, k, QW);
  ekinr = 0.5 * rrr * (rur * rur + rvr * rvr + rwr * rwr);
  rur *= rrr;
  rvr *= rrr;
  rwr *= rrr;

  rer = qxp(i, j, k, QREINT) + ekinr;
  // Add transverse predictor
  rrnewr = rrr - flxrho;
  runewr = rur - flxu;
  rvnewr = rvr - flxv;
  rwnewr = rwr - flxw;
  renewr = rer - flxe;

  // Convert back to primitive
  qp(i, j, k, QRHO) = rrnewr; // + hdt*srcr;
  qp(i, j, k, QU) = runewr / rrnewr;
  qp(i, j, k, QV) = rvnewr / rrnewr;
  qp(i, j, k, QW) = rwnewr / rrnewr;
  rhoekinr =
    0.5 * (runewr * runewr + rvnewr * rvnewr + rwnewr * rwnewr) / rrnewr;

  pnewr = qxp(i, j, k, QPRES) - cdtdz * (dup + pav * du * (c - 1.0));
  qp(i, j, k, QPRES) = amrex::max(pnewr, SMALL_PRES);
  qp(i, j, k, QREINT) = renewr - rhoekinr;

  // QM

  // Conversion to Conservative
  rrl = qxm(i + 1, j, k, QRHO);
  rul = qxm(i + 1, j, k, QU);
  rvl = qxm(i + 1, j, k, QV);
  rwl = qxm(i + 1, j, k, QW);
  ekinl = 0.5 * rrl * (rul * rul + rvl * rvl + rwl * rwl);
  rul *= rrl;
  rvl *= rrl;
  rwl *= rrl;
  rel = qxm(i + 1, j, k, QREINT) + ekinl;

  // Transverse fluxes
  rrnewl = rrl - flxrho;
  runewl = rul - flxu;
  rvnewl = rvl - flxv;
  rwnewl = rwl - flxw;
  renewl = rel - flxe;

  qm(i + 1, j, k, QRHO) = rrnewl;
  qm(i + 1, j, k, QU) = runewl / rrnewl;
  qm(i + 1, j, k, QV) = rvnewl / rrnewl;
  qm(i + 1, j, k, QW) = rwnewl / rrnewl;
  rhoekinl =
    0.5 * (runewl * runewl + rvnewl * rvnewl + rwnewl * rwnewl) / rrnewl;

  pnewl = qxm(i + 1, j, k, QPRES) - cdtdz * (dup + pav * du * (c - 1.0));
  qm(i + 1, j, k, QPRES) = amrex::max(pnewl, SMALL_PRES);
  qm(i + 1, j, k, QREINT) = renewl - rhoekinl;
}

// Transverse Correction for Predicted Y-states, using Z-Flux
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_transz2(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::Real> const& qm,
  amrex::Array4<amrex::Real> const& qp,
  amrex::Array4<const amrex::Real> const& qym,
  amrex::Array4<const amrex::Real> const& qyp,
  amrex::Array4<const amrex::Real> const& flxz,
  amrex::Array4<const amrex::Real> const& qa,
  amrex::Array4<const amrex::Real> const& qint,
  const amrex::Real cdtdz)
{
  // Local variables to avoid reading from global memory too much
  amrex::Real compn, rr, rrnew, compo, rwr, c;
  amrex::Real rrr, rur, rvr, ekinr, rer, rrnewr;
  amrex::Real runewr, rvnewr, rwnewr, renewr, rhoekinr;
  amrex::Real rrl, rul, rvl, ekinl, rel, rrnewl, runewl, rvnewl, rwnewl, renewl;
  amrex::Real rhoekinl, rwl;
  amrex::Real pggp, pggm, ugp, ugm;
  amrex::Real dAup, pav, dAu;
  amrex::Real flxrho, flxu, flxv, flxw, flxe;

  flxrho = cdtdz * (flxz(i, j, k + 1, URHO) - flxz(i, j, k, URHO));
  flxu = cdtdz * (flxz(i, j, k + 1, UMX) - flxz(i, j, k, UMX));
  flxv = cdtdz * (flxz(i, j, k + 1, UMY) - flxz(i, j, k, UMY));
  flxw = cdtdz * (flxz(i, j, k + 1, UMZ) - flxz(i, j, k, UMZ));
  flxe = cdtdz * (flxz(i, j, k + 1, UEDEN) - flxz(i, j, k, UEDEN));
  c = qa(i, j, k, QGAMC);
  // Update passive variables
  for (int ipass = 0; ipass < NPASSIVE; ++ipass) {
    int n = indxmap::upass_map(ipass);
    int nqp = indxmap::qpass_map(ipass);
    compn = cdtdz * (flxz(i, j, k + 1, n) - flxz(i, j, k, n));
    rr = qyp(i, j, k, QRHO);
    rrnew = rr - flxrho;
    compo = rr * qyp(i, j, k, nqp) - compn;
    qp(i, j, k, nqp) = compo / rrnew;

    rr = qym(i, j + 1, k, QRHO);
    rrnew = rr - flxrho;
    compo = rr * qym(i, j + 1, k, nqp) - compn;
    qm(i, j + 1, k, nqp) = compo / rrnew;
  }

  // Update hydro vars
  pggp = qint(i, j, k + 1, GDPRES);
  pggm = qint(i, j, k, GDPRES);
  ugp = qint(i, j, k + 1, GDW);
  ugm = qint(i, j, k, GDW);

  dAup = pggp * ugp - pggm * ugm;
  pav = 0.5 * (pggp + pggm);
  dAu = ugp - ugm;

  // QP

  // Convert to conservative
  rrr = qyp(i, j, k, QRHO);
  rur = qyp(i, j, k, QU);
  rvr = qyp(i, j, k, QV);
  rwr = qyp(i, j, k, QW);
  ekinr = 0.5 * rrr * (rur * rur + rvr * rvr + rwr * rwr);
  rur *= rrr;
  rvr *= rrr;
  rwr *= rrr;
  rer = qyp(i, j, k, QREINT) + ekinr;
  // Add transverse predictor
  rrnewr = rrr - flxrho;
  runewr = rur - flxu;
  rvnewr = rvr - flxv;
  rwnewr = rwr - flxw;
  renewr = rer - flxe;

  // Convert back to primitive
  qp(i, j, k, QRHO) = rrnewr;
  qp(i, j, k, QU) = runewr / rrnewr;
  qp(i, j, k, QV) = rvnewr / rrnewr;
  qp(i, j, k, QW) = rwnewr / rrnewr;
  rhoekinr =
    0.5 * (runewr * runewr + rvnewr * rvnewr + rwnewr * rwnewr) / rrnewr;
  amrex::Real pnewr =
    qyp(i, j, k, QPRES) - cdtdz * (dAup + pav * dAu * (c - 1.0));
  qp(i, j, k, QPRES) = amrex::max(pnewr, SMALL_PRES);
  qp(i, j, k, QREINT) = renewr - rhoekinr;

  // QM

  // Conversion to Conservative
  rrl = qym(i, j + 1, k, QRHO);
  rul = qym(i, j + 1, k, QU);
  rvl = qym(i, j + 1, k, QV);
  rwl = qym(i, j + 1, k, QW);
  ekinl = 0.5 * rrl * (rul * rul + rvl * rvl + rwl * rwl);
  rul *= rrl;
  rvl *= rrl;
  rwl *= rrl;
  rel = qym(i, j + 1, k, QREINT) + ekinl;

  // Transverse fluxes
  rrnewl = rrl - flxrho;
  runewl = rul - flxu;
  rvnewl = rvl - flxv;
  rwnewl = rwl - flxw;
  renewl = rel - flxe;

  qm(i, j + 1, k, QRHO) = rrnewl;
  qm(i, j + 1, k, QU) = runewl / rrnewl;
  qm(i, j + 1, k, QV) = rvnewl / rrnewl;
  qm(i, j + 1, k, QW) = rwnewl / rrnewl;
  rhoekinl =
    0.5 * (runewl * runewl + rvnewl * rvnewl + rwnewl * rwnewl) / rrnewl;

  qm(i, j + 1, k, QREINT) = renewl - rhoekinl;
  amrex::Real pnewl =
    qym(i, j + 1, k, QPRES) - cdtdz * (dAup + pav * dAu * (c - 1.0));
  qm(i, j + 1, k, QPRES) = amrex::max(pnewl, SMALL_PRES);
}

// Z corrected from xy
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_transxy(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::Real> const& qm,
  amrex::Array4<amrex::Real> const& qp,
  amrex::Array4<const amrex::Real> const& qzm,
  amrex::Array4<const amrex::Real> const& qzp,
  amrex::Array4<const amrex::Real> const& flxx,
  amrex::Array4<const amrex::Real> const& flxy,
  amrex::Array4<const amrex::Real> const& qx,
  amrex::Array4<const amrex::Real> const& qy,
  amrex::Array4<const amrex::Real> const& qa,
  amrex::Array4<const amrex::Real> const& srcq,
  const amrex::Real hdt,
  const amrex::Real cdtdx,
  const amrex::Real cdtdy)
{
  amrex::Real compn, compo, rwr, c;
  amrex::Real rrr, rur, rvr, ekinr, rer, rrnewr;
  amrex::Real runewr, rvnewr, rwnewr, renewr, rhoekinr;
  amrex::Real rrl, rul, rvl, ekinl, rel, rrnewl, runewl, rvnewl, rwnewl, renewl;
  amrex::Real rhoekinl, rwl;
  amrex::Real pggpx, pggmx, ugpx, ugmx;
  amrex::Real dAupx, pavx, dAux;
  amrex::Real pggpy, pggmy, ugpy, ugmy;
  amrex::Real dAupy, pavy, dAuy;
  amrex::Real flxrho, flxu, flxv, flxw, flxe;
  amrex::Real srcpass, srcrho, srcu, srcv, srcw, srce, srcp;

  flxrho = cdtdx * (flxx(i + 1, j, k, URHO) - flxx(i, j, k, URHO)) +
           cdtdy * (flxy(i, j + 1, k, URHO) - flxy(i, j, k, URHO));
  flxu = cdtdx * (flxx(i + 1, j, k, UMX) - flxx(i, j, k, UMX)) +
         cdtdy * (flxy(i, j + 1, k, UMX) - flxy(i, j, k, UMX));
  flxv = cdtdx * (flxx(i + 1, j, k, UMY) - flxx(i, j, k, UMY)) +
         cdtdy * (flxy(i, j + 1, k, UMY) - flxy(i, j, k, UMY));
  flxw = cdtdx * (flxx(i + 1, j, k, UMZ) - flxx(i, j, k, UMZ)) +
         cdtdy * (flxy(i, j + 1, k, UMZ) - flxy(i, j, k, UMZ));
  flxe = cdtdx * (flxx(i + 1, j, k, UEDEN) - flxx(i, j, k, UEDEN)) +
         cdtdy * (flxy(i, j + 1, k, UEDEN) - flxy(i, j, k, UEDEN));
  c = qa(i, j, k, QGAMC);
  srcrho = srcq(i, j, k, QRHO);
  srcu = srcq(i, j, k, QU);
  srcv = srcq(i, j, k, QV);
  srcw = srcq(i, j, k, QW);
  srce = srcq(i, j, k, QREINT);
  srcp = srcq(i, j, k, QPRES);

  // Update passive variables
  rrr = qzp(i, j, k, QRHO);
  rrl = qzm(i, j, k + 1, QRHO);
  rrnewl = rrl - flxrho;
  rrnewr = rrr - flxrho;
  for (int ipass = 0; ipass < NPASSIVE; ++ipass) {
    int n = indxmap::upass_map(ipass);
    int nqp = indxmap::qpass_map(ipass);
    srcpass = srcq(i, j, k, nqp);
    compn = cdtdx * (flxx(i + 1, j, k, n) - flxx(i, j, k, n)) +
            cdtdy * (flxy(i, j + 1, k, n) - flxy(i, j, k, n));
    // qp
    compo = rrr * qzp(i, j, k, nqp) - compn;
    qp(i, j, k, nqp) = compo / rrnewr + hdt * srcpass;
    // qm
    compo = rrl * qzm(i, j, k + 1, nqp) - compn;
    qm(i, j, k + 1, nqp) = compo / rrnewl + hdt * srcpass;
  }

  // Update hydro vars
  pggpx = qx(i + 1, j, k, GDPRES);
  pggmx = qx(i, j, k, GDPRES);
  ugpx = qx(i + 1, j, k, GDU);
  ugmx = qx(i, j, k, GDU);

  dAupx = pggpx * ugpx - pggmx * ugmx;
  pavx = 0.5 * (pggpx + pggmx);
  dAux = ugpx - ugmx;

  pggpy = qy(i, j + 1, k, GDPRES);
  pggmy = qy(i, j, k, GDPRES);
  ugpy = qy(i, j + 1, k, GDV);
  ugmy = qy(i, j, k, GDV);

  dAupy = pggpy * ugpy - pggmy * ugmy;
  pavy = 0.5 * (pggpy + pggmy);
  dAuy = ugpy - ugmy;
  amrex::Real pxnew = cdtdx * (dAupx + pavx * dAux * (c - 1.0));
  amrex::Real pynew = cdtdy * (dAupy + pavy * dAuy * (c - 1.0));

  // qp state
  rur = rrr * qzp(i, j, k, QU);
  rvr = rrr * qzp(i, j, k, QV);
  rwr = rrr * qzp(i, j, k, QW);
  ekinr = 0.5 * (rur * rur + rvr * rvr + rwr * rwr) / rrr;
  rer = qzp(i, j, k, QREINT) + ekinr;

  runewr = rur - flxu;
  rvnewr = rvr - flxv;
  rwnewr = rwr - flxw;
  renewr = rer - flxe;

  qp(i, j, k, QRHO) = rrnewr + hdt * srcrho;
  qp(i, j, k, QU) = runewr / rrnewr + hdt * srcu;
  qp(i, j, k, QV) = rvnewr / rrnewr + hdt * srcv;
  qp(i, j, k, QW) = rwnewr / rrnewr + hdt * srcw;
  rhoekinr =
    0.5 * (runewr * runewr + rvnewr * rvnewr + rwnewr * rwnewr) / rrnewr;
  qp(i, j, k, QREINT) = renewr - rhoekinr + hdt * srce;
  amrex::Real temppres = qzp(i, j, k, QPRES) - pxnew - pynew + hdt * srcp;
  qp(i, j, k, QPRES) = amrex::max(temppres, SMALL_PRES);

  // qm state
  rul = rrl * qzm(i, j, k + 1, QU);
  rvl = rrl * qzm(i, j, k + 1, QV);
  rwl = rrl * qzm(i, j, k + 1, QW);
  ekinl = 0.5 * (rul * rul + rvl * rvl + rwl * rwl) / rrl;
  rel = qzm(i, j, k + 1, QREINT) + ekinl;

  runewl = rul - flxu;
  rvnewl = rvl - flxv;
  rwnewl = rwl - flxw;
  renewl = rel - flxe;

  qm(i, j, k + 1, QRHO) = rrnewl + hdt * srcrho;
  qm(i, j, k + 1, QU) = runewl / rrnewl + hdt * srcu;
  qm(i, j, k + 1, QV) = rvnewl / rrnewl + hdt * srcv;
  qm(i, j, k + 1, QW) = rwnewl / rrnewl + hdt * srcw;
  rhoekinl =
    0.5 * (runewl * runewl + rvnewl * rvnewl + rwnewl * rwnewl) / rrnewl;
  qm(i, j, k + 1, QREINT) = renewl - rhoekinl + hdt * srce;
  temppres = qzm(i, j, k + 1, QPRES) - pxnew - pynew + hdt * srcp;
  qm(i, j, k + 1, QPRES) = amrex::max(temppres, SMALL_PRES);
}

// Y corrected by x and z
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_transxz(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::Real> const& qm,
  amrex::Array4<amrex::Real> const& qp,
  amrex::Array4<const amrex::Real> const& qym,
  amrex::Array4<const amrex::Real> const& qyp,
  amrex::Array4<const amrex::Real> const& flxx,
  amrex::Array4<const amrex::Real> const& flxz,
  amrex::Array4<const amrex::Real> const& qx,
  amrex::Array4<const amrex::Real> const& qz,
  amrex::Array4<const amrex::Real> const& qa,
  amrex::Array4<const amrex::Real> const& srcq,
  const amrex::Real hdt,
  const amrex::Real cdtdx,
  const amrex::Real cdtdz)
{
  amrex::Real compn, compo, rwr, c;
  amrex::Real rrr, rur, rvr, ekinr, rer, rrnewr;
  amrex::Real runewr, rvnewr, rwnewr, renewr, rhoekinr;
  amrex::Real rrl, rul, rvl, ekinl, rel, rrnewl, runewl, rvnewl, rwnewl, renewl;
  amrex::Real rhoekinl, rwl;
  amrex::Real pggpx, pggmx, ugpx, ugmx;
  amrex::Real dAupx, pavx, dAux;
  amrex::Real pggpz, pggmz, ugpz, ugmz;
  amrex::Real dAupz, pavz, dAuz;
  amrex::Real flxrho, flxu, flxv, flxw, flxe;
  amrex::Real srcpass, srcrho, srcu, srcv, srcw, srce, srcp;

  flxrho = cdtdx * (flxx(i + 1, j, k, URHO) - flxx(i, j, k, URHO)) +
           cdtdz * (flxz(i, j, k + 1, URHO) - flxz(i, j, k, URHO));
  flxu = cdtdx * (flxx(i + 1, j, k, UMX) - flxx(i, j, k, UMX)) +
         cdtdz * (flxz(i, j, k + 1, UMX) - flxz(i, j, k, UMX));
  flxv = cdtdx * (flxx(i + 1, j, k, UMY) - flxx(i, j, k, UMY)) +
         cdtdz * (flxz(i, j, k + 1, UMY) - flxz(i, j, k, UMY));
  flxw = cdtdx * (flxx(i + 1, j, k, UMZ) - flxx(i, j, k, UMZ)) +
         cdtdz * (flxz(i, j, k + 1, UMZ) - flxz(i, j, k, UMZ));
  flxe = cdtdx * (flxx(i + 1, j, k, UEDEN) - flxx(i, j, k, UEDEN)) +
         cdtdz * (flxz(i, j, k + 1, UEDEN) - flxz(i, j, k, UEDEN));
  c = qa(i, j, k, QGAMC);
  srcrho = srcq(i, j, k, QRHO);
  srcu = srcq(i, j, k, QU);
  srcv = srcq(i, j, k, QV);
  srcw = srcq(i, j, k, QW);
  srce = srcq(i, j, k, QREINT);
  srcp = srcq(i, j, k, QPRES);

  // Update passive variables
  rrr = qyp(i, j, k, QRHO);
  rrl = qym(i, j + 1, k, QRHO);
  rrnewl = rrl - flxrho;
  rrnewr = rrr - flxrho;
  for (int ipass = 0; ipass < NPASSIVE; ++ipass) {
    int n = indxmap::upass_map(ipass);
    int nqp = indxmap::qpass_map(ipass);
    srcpass = srcq(i, j, k, nqp);
    compn = cdtdx * (flxx(i + 1, j, k, n) - flxx(i, j, k, n)) +
            cdtdz * (flxz(i, j, k + 1, n) - flxz(i, j, k, n));
    // qp
    compo = rrr * qyp(i, j, k, nqp) - compn;
    qp(i, j, k, nqp) = compo / rrnewr + hdt * srcpass;
    // qm
    compo = rrl * qym(i, j + 1, k, nqp) - compn;
    qm(i, j + 1, k, nqp) = compo / rrnewl + hdt * srcpass;
  }

  // Update hydro vars
  pggpx = qx(i + 1, j, k, GDPRES);
  pggmx = qx(i, j, k, GDPRES);
  ugpx = qx(i + 1, j, k, GDU);
  ugmx = qx(i, j, k, GDU);

  dAupx = pggpx * ugpx - pggmx * ugmx;
  pavx = 0.5 * (pggpx + pggmx);
  dAux = ugpx - ugmx;

  pggpz = qz(i, j, k + 1, GDPRES);
  pggmz = qz(i, j, k, GDPRES);
  ugpz = qz(i, j, k + 1, GDW);
  ugmz = qz(i, j, k, GDW);

  dAupz = pggpz * ugpz - pggmz * ugmz;
  pavz = 0.5 * (pggpz + pggmz);
  dAuz = ugpz - ugmz;
  amrex::Real pxnew = cdtdx * (dAupx + pavx * dAux * (c - 1.0));
  amrex::Real pznew = cdtdz * (dAupz + pavz * dAuz * (c - 1.0));

  // qp state
  rur = rrr * qyp(i, j, k, QU);
  rvr = rrr * qyp(i, j, k, QV);
  rwr = rrr * qyp(i, j, k, QW);
  ekinr = 0.5 * (rur * rur + rvr * rvr + rwr * rwr) / rrr;
  rer = qyp(i, j, k, QREINT) + ekinr;

  runewr = rur - flxu;
  rvnewr = rvr - flxv;
  rwnewr = rwr - flxw;
  renewr = rer - flxe;

  qp(i, j, k, QRHO) = rrnewr + hdt * srcrho;
  qp(i, j, k, QU) = runewr / rrnewr + hdt * srcu;
  qp(i, j, k, QV) = rvnewr / rrnewr + hdt * srcv;
  qp(i, j, k, QW) = rwnewr / rrnewr + hdt * srcw;
  rhoekinr =
    0.5 * (runewr * runewr + rvnewr * rvnewr + rwnewr * rwnewr) / rrnewr;
  qp(i, j, k, QREINT) = renewr - rhoekinr + hdt * srce;
  amrex::Real temppres = qyp(i, j, k, QPRES) - pxnew - pznew + hdt * srcp;
  qp(i, j, k, QPRES) = amrex::max(temppres, SMALL_PRES);

  // qm state
  rul = rrl * qym(i, j + 1, k, QU);
  rvl = rrl * qym(i, j + 1, k, QV);
  rwl = rrl * qym(i, j + 1, k, QW);
  ekinl = 0.5 * (rul * rul + rvl * rvl + rwl * rwl) / rrl;
  rel = qym(i, j + 1, k, QREINT) + ekinl;

  runewl = rul - flxu;
  rvnewl = rvl - flxv;
  rwnewl = rwl - flxw;
  renewl = rel - flxe;

  qm(i, j + 1, k, QRHO) = rrnewl + hdt * srcrho;
  qm(i, j + 1, k, QU) = runewl / rrnewl + hdt * srcu;
  qm(i, j + 1, k, QV) = rvnewl / rrnewl + hdt * srcv;
  qm(i, j + 1, k, QW) = rwnewl / rrnewl + hdt * srcw;
  rhoekinl =
    0.5 * (runewl * runewl + rvnewl * rvnewl + rwnewl * rwnewl) / rrnewl;
  qm(i, j + 1, k, QREINT) = renewl - rhoekinl + hdt * srce;
  temppres = qym(i, j + 1, k, QPRES) - pxnew - pznew + hdt * srcp;
  qm(i, j + 1, k, QPRES) = amrex::max(temppres, SMALL_PRES);
}

// x corrected by y and z
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_transyz(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::Real> const& qm,
  amrex::Array4<amrex::Real> const& qp,
  amrex::Array4<const amrex::Real> const& qxm,
  amrex::Array4<const amrex::Real> const& qxp,
  amrex::Array4<const amrex::Real> const& flxy,
  amrex::Array4<const amrex::Real> const& flxz,
  amrex::Array4<const amrex::Real> const& qy,
  amrex::Array4<const amrex::Real> const& qz,
  amrex::Array4<const amrex::Real> const& qa,
  amrex::Array4<const amrex::Real> const& srcq,
  const amrex::Real hdt,
  const amrex::Real cdtdy,
  const amrex::Real cdtdz)
{
  amrex::Real compn, compo, rwr, c;
  amrex::Real rrr, rur, rvr, ekinr, rer, rrnewr;
  amrex::Real runewr, rvnewr, rwnewr, renewr, rhoekinr;
  amrex::Real rrl, rul, rvl, ekinl, rel, rrnewl, runewl, rvnewl, rwnewl, renewl;
  amrex::Real rhoekinl, rwl;
  amrex::Real pggpz, pggmz, ugpz, ugmz;
  amrex::Real dAupz, pavz, dAuz;
  amrex::Real pggpy, pggmy, ugpy, ugmy;
  amrex::Real dAupy, pavy, dAuy;
  amrex::Real flxrho, flxu, flxv, flxw, flxe;
  amrex::Real srcpass, srcrho, srcu, srcv, srcw, srce, srcp;

  flxrho = cdtdy * (flxy(i, j + 1, k, URHO) - flxy(i, j, k, URHO)) +
           cdtdz * (flxz(i, j, k + 1, URHO) - flxz(i, j, k, URHO));
  flxu = cdtdy * (flxy(i, j + 1, k, UMX) - flxy(i, j, k, UMX)) +
         cdtdz * (flxz(i, j, k + 1, UMX) - flxz(i, j, k, UMX));
  flxv = cdtdy * (flxy(i, j + 1, k, UMY) - flxy(i, j, k, UMY)) +
         cdtdz * (flxz(i, j, k + 1, UMY) - flxz(i, j, k, UMY));
  flxw = cdtdy * (flxy(i, j + 1, k, UMZ) - flxy(i, j, k, UMZ)) +
         cdtdz * (flxz(i, j, k + 1, UMZ) - flxz(i, j, k, UMZ));
  flxe = cdtdy * (flxy(i, j + 1, k, UEDEN) - flxy(i, j, k, UEDEN)) +
         cdtdz * (flxz(i, j, k + 1, UEDEN) - flxz(i, j, k, UEDEN));
  c = qa(i, j, k, QGAMC);
  srcrho = srcq(i, j, k, QRHO);
  srcu = srcq(i, j, k, QU);
  srcv = srcq(i, j, k, QV);
  srcw = srcq(i, j, k, QW);
  srce = srcq(i, j, k, QREINT);
  srcp = srcq(i, j, k, QPRES);

  // Update passive variables
  rrr = qxp(i, j, k, QRHO);
  rrl = qxm(i + 1, j, k, QRHO);
  rrnewl = rrl - flxrho;
  rrnewr = rrr - flxrho;

  for (int ipass = 0; ipass < NPASSIVE; ++ipass) {
    int n = indxmap::upass_map(ipass);
    int nqp = indxmap::qpass_map(ipass);
    srcpass = srcq(i, j, k, nqp);
    compn = cdtdy * (flxy(i, j + 1, k, n) - flxy(i, j, k, n)) +
            cdtdz * (flxz(i, j, k + 1, n) - flxz(i, j, k, n));
    // qp
    compo = rrr * qxp(i, j, k, nqp) - compn;
    qp(i, j, k, nqp) = compo / rrnewr + hdt * srcpass;
    // qm
    compo = rrl * qxm(i + 1, j, k, nqp) - compn;
    qm(i + 1, j, k, nqp) = compo / rrnewl + hdt * srcpass;
  }

  // Update hydro vars
  pggpz = qz(i, j, k + 1, GDPRES);
  pggmz = qz(i, j, k, GDPRES);
  ugpz = qz(i, j, k + 1, GDW);
  ugmz = qz(i, j, k, GDW);

  dAupz = pggpz * ugpz - pggmz * ugmz;
  pavz = 0.5 * (pggpz + pggmz);
  dAuz = ugpz - ugmz;

  pggpy = qy(i, j + 1, k, GDPRES);
  pggmy = qy(i, j, k, GDPRES);
  ugpy = qy(i, j + 1, k, GDV);
  ugmy = qy(i, j, k, GDV);

  dAupy = pggpy * ugpy - pggmy * ugmy;
  pavy = 0.5 * (pggpy + pggmy);
  dAuy = ugpy - ugmy;
  amrex::Real pznew = cdtdz * (dAupz + pavz * dAuz * (c - 1.0));
  amrex::Real pynew = cdtdy * (dAupy + pavy * dAuy * (c - 1.0));

  // qp state
  rur = rrr * qxp(i, j, k, QU);
  rvr = rrr * qxp(i, j, k, QV);
  rwr = rrr * qxp(i, j, k, QW);
  ekinr = 0.5 * (rur * rur + rvr * rvr + rwr * rwr) / rrr;
  rer = qxp(i, j, k, QREINT) + ekinr;

  runewr = rur - flxu;
  rvnewr = rvr - flxv;
  rwnewr = rwr - flxw;
  renewr = rer - flxe;

  qp(i, j, k, QRHO) = rrnewr + hdt * srcrho;
  qp(i, j, k, QU) = runewr / rrnewr + hdt * srcu;
  qp(i, j, k, QV) = rvnewr / rrnewr + hdt * srcv;
  qp(i, j, k, QW) = rwnewr / rrnewr + hdt * srcw;
  rhoekinr =
    0.5 * (runewr * runewr + rvnewr * rvnewr + rwnewr * rwnewr) / rrnewr;
  qp(i, j, k, QREINT) = renewr - rhoekinr + hdt * srce;
  amrex::Real temppres = qxp(i, j, k, QPRES) - pznew - pynew + hdt * srcp;
  qp(i, j, k, QPRES) = amrex::max(temppres, SMALL_PRES);

  // qm state
  rul = rrl * qxm(i + 1, j, k, QU);
  rvl = rrl * qxm(i + 1, j, k, QV);
  rwl = rrl * qxm(i + 1, j, k, QW);
  ekinl = 0.5 * (rul * rul + rvl * rvl + rwl * rwl) / rrl;
  rel = qxm(i + 1, j, k, QREINT) + ekinl;

  runewl = rul - flxu;
  rvnewl = rvl - flxv;
  rwnewl = rwl - flxw;
  renewl = rel - flxe;

  qm(i + 1, j, k, QRHO) = rrnewl + hdt * srcrho;
  qm(i + 1, j, k, QU) = runewl / rrnewl + hdt * srcu;
  qm(i + 1, j, k, QV) = rvnewl / rrnewl + hdt * srcv;
  qm(i + 1, j, k, QW) = rwnewl / rrnewl + hdt * srcw;
  rhoekinl =
    0.5 * (runewl * runewl + rvnewl * rvnewl + rwnewl * rwnewl) / rrnewl;
  qm(i + 1, j, k, QREINT) = renewl - rhoekinl + hdt * srce;
  temppres = qxm(i + 1, j, k, QPRES) - pznew - pynew + hdt * srcp;
  qm(i + 1, j, k, QPRES) = amrex::max(temppres, SMALL_PRES);
}

// Use interface states from Riemann solver for pdivu.
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_pdivu(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::Real> const& pdivu,
  AMREX_D_DECL(
    amrex::Array4<const amrex::Real> const& q1,
    amrex::Array4<const amrex::Real> const& q2,
    amrex::Array4<const amrex::Real> const& q3),
  AMREX_D_DECL(
    amrex::Array4<const amrex::Real> const& a1,
    amrex::Array4<const amrex::Real> const& a2,
    amrex::Array4<const amrex::Real> const& a3),
  amrex::Array4<const amrex::Real> const& vol)
{
  pdivu(i, j, k) = 0.5 *
                   (AMREX_D_TERM(
                     (q1(i + 1, j, k, GDPRES) + q1(i, j, k, GDPRES)) *
                       (q1(i + 1, j, k, GDU) * a1(i + 1, j, k) -
                        q1(i, j, k, GDU) * a1(i, j, k)),
                     +(q2(i, j + 1, k, GDPRES) + q2(i, j, k, GDPRES)) *
                       (q2(i, j + 1, k, GDV) * a2(i, j + 1, k) -
                        q2(i, j, k, GDV) * a2(i, j, k)),
                     +(q3(i, j, k + 1, GDPRES) + q3(i, j, k, GDPRES)) *
                       (q3(i, j, k + 1, GDW) * a3(i, j, k + 1) -
                        q3(i, j, k, GDW) * a3(i, j, k)))) /
                   vol(i, j, k);
}

// Consup Support Functions
// Designed for CPU or GPU
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_artif_visc(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::Real> const& flx,
  amrex::Array4<const amrex::Real> const& div,
  amrex::Array4<const amrex::Real> const& u,
  amrex::Real const dx,
  amrex::Real const difmag,
  const int dir)
{
  amrex::Real div1;
  if (dir == 0) { // X direction
    AMREX_D_PICK(div1 = div(i, j, k);
                 , div1 = 0.5 * (div(i, j, k) + div(i, j + 1, k));
                 , div1 = 0.25 * (div(i, j, k) + div(i, j + 1, k) +
                                  div(i, j, k + 1) + div(i, j + 1, k + 1)););
    div1 = difmag * amrex::min(0.0, div1);
    for (int n = 0; n < NVAR; ++n) {
      if (n != UTEMP) {
        // Passing in difmag TODO not pass in.
        flx(i, j, k, n) += dx * div1 * (u(i, j, k, n) - u(i - 1, j, k, n));
      }
    }
  } else if (dir == 1) { // Y direction
    AMREX_D_PICK(, div1 = 0.5 * (div(i, j, k) + div(i + 1, j, k));
                 , div1 = 0.25 * (div(i, j, k) + div(i + 1, j, k) +
                                  div(i, j, k + 1) + div(i + 1, j, k + 1)););
    div1 = difmag * amrex::min(0.0, div1);
    for (int n = 0; n < NVAR; ++n) {
      if (n != UTEMP) {
        flx(i, j, k, n) +=
          dx * div1 * (u(i, j, k, n) - u(i, j - 1, k, n)); // Here dx is dy
      }
    }
  } else { // Z direction
    div1 = 0.25 * (div(i, j, k) + div(i + 1, j, k) + div(i, j + 1, k) +
                   div(i + 1, j + 1, k));
    div1 = difmag * amrex::min(0.0, div1);
    for (int n = 0; n < NVAR; ++n) {
      if (n != UTEMP) {
        flx(i, j, k, n) +=
          dx * div1 * (u(i, j, k, n) - u(i, j, k - 1, n)); // Here dx is dz
      }
    }
  }
  flx(i, j, k, UTEMP) = 0.0;
}

// Host Functions
void pc_umeth_3D(
  amrex::Box const& bx,
  const int* bclo,
  const int* bchi,
  const int* domlo,
  const int* domhi,
  amrex::Array4<const amrex::Real> const& q,
  amrex::Array4<const amrex::Real> const& qaux,
  amrex::Array4<const amrex::Real> const& srcQ,
  amrex::Array4<amrex::Real> const& flx1,
  amrex::Array4<amrex::Real> const& flx2,
  amrex::Array4<amrex::Real> const& flx3,
  amrex::Array4<amrex::Real> const& q1,
  amrex::Array4<amrex::Real> const& q2,
  amrex::Array4<amrex::Real> const& q3,
  amrex::Array4<const amrex::Real> const& a1,
  amrex::Array4<const amrex::Real> const& a2,
  amrex::Array4<const amrex::Real> const& a3,
  amrex::Array4<amrex::Real> const& pdivu,
  amrex::Array4<const amrex::Real> const& vol,
  const amrex::Real* del,
  const amrex::Real dt,
  const int ppm_type,
  const int use_flattening);

void pc_umeth_2D(
  amrex::Box const& bx,
  const int* bclo,
  const int* bchi,
  const int* domlo,
  const int* domhi,
  amrex::Array4<const amrex::Real> const& q,
  amrex::Array4<const amrex::Real> const& qaux,
  amrex::Array4<const amrex::Real> const& srcQ,
  amrex::Array4<amrex::Real> const& flx1,
  amrex::Array4<amrex::Real> const& flx2,
  // amrex::Array4<const amrex::Real> const& dloga,
  amrex::Array4<amrex::Real> const& q1,
  amrex::Array4<amrex::Real> const& q2,
  amrex::Array4<const amrex::Real> const& a1,
  amrex::Array4<const amrex::Real> const& a2,
  amrex::Array4<amrex::Real> const& pdivu,
  amrex::Array4<const amrex::Real> const& vol,
  const amrex::Real* del,
  const amrex::Real dt,
  const int ppm_type,
  const int use_flattening);

#endif
