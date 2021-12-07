#ifndef _GODUNOV_H_
#define _GODUNOV_H_

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

  if ((!q.contains(i + 3 * bdim[0], j + 3 * bdim[1], k + 3 * bdim[2]) ||
       !q.contains(i - 3 * bdim[0], j - 3 * bdim[1], k - 3 * bdim[2]))) {
    return 1.0;
  }
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

  amrex::Real denom =
    amrex::max<amrex::Real>(small_pres, amrex::Math::abs(qpp - qmm));
  amrex::Real zeta = amrex::Math::abs(dp) / denom;
  const amrex::Real z = amrex::min<amrex::Real>(
    1.0, amrex::max<amrex::Real>(0.0, dzcut * (zeta - zcut1)));

  amrex::Real tst = (qmun - qpun >= 0.0) ? 1.0 : 0.0;

  amrex::Real tmp = amrex::min<amrex::Real>(qp, qm);

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

  denom =
    amrex::max<amrex::Real>(small_pres, amrex::Math::abs(qppshft - qmmshft));
  zeta = amrex::Math::abs(dp) / denom;
  const amrex::Real z2 = amrex::min<amrex::Real>(
    1.0, amrex::max<amrex::Real>(0.0, dzcut * (zeta - zcut1)));

  tst = (qmunshft - qpunshft >= 0.0) ? 1.0 : 0.0;

  tmp = amrex::min<amrex::Real>(qpshft, qmshft);

  const amrex::Real chi2 = (amrex::Math::abs(dp) / tmp) > shktst ? tst : 0.0;

  return 1.0 - amrex::max<amrex::Real>(chi2 * z2, chi * z);
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
  const int dir,
  PassMap const& pmap)
{
  amrex::Real cav, ustar;
  amrex::Real spl[NUM_SPECIES];
  amrex::Real spr[NUM_SPECIES];
  amrex::Real ul, ur, vl, vr, v2l, v2r;
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

  ur = qr(i, j, k, IU);
  vr = qr(i, j, k, IV);
  v2r = qr(i, j, k, IV2);

  // Outflow Hack
  if (dir == 2)
    idx = k;
  else
    idx = (dir == 0) ? i : j;
  if (bclo == Outflow && idx == domlo) {
    ul = ur;
    vl = vr;
    v2l = v2r;
  }
  if (bchi == Outflow && idx == domhi + 1) {
    ur = ul;
    vr = vl;
    v2r = v2l;
  }

  const int bc_test_val = 1;
  riemann(
    ql(i, j, k, QRHO), ul, vl, v2l, ql(i, j, k, QPRES), spl, qr(i, j, k, QRHO),
    ur, vr, v2r, qr(i, j, k, QPRES), spr, bc_test_val, cav, ustar,
    flx(i, j, k, URHO), flx(i, j, k, f_idx[0]), flx(i, j, k, f_idx[1]),
    flx(i, j, k, f_idx[2]), flx(i, j, k, UEDEN), flx(i, j, k, UEINT),
    q(i, j, k, GU), q(i, j, k, GV), q(i, j, k, GV2), q(i, j, k, GDPRES),
    q(i, j, k, GDGAME));

  amrex::Real flxrho = flx(i, j, k, URHO);
  for (int ipass = 0; ipass < NPASSIVE; ++ipass) {
    int n = pmap.upassMap[ipass];
    int nqp = pmap.qpassMap[ipass];
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
  const amrex::Real cdtdx,
  PassMap const& pmap)
{
  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  const amrex::IntVect ivp0{AMREX_D_DECL(i + 1, j, k)};
  const amrex::IntVect ivp1{AMREX_D_DECL(i, j + 1, k)};
  const amrex::Real flxrho = cdtdx * (flxx(ivp0, URHO) - flxx(iv, URHO));
  const amrex::Real flxu = cdtdx * (flxx(ivp0, UMX) - flxx(iv, UMX));
  const amrex::Real flxv = cdtdx * (flxx(ivp0, UMY) - flxx(iv, UMY));
  const amrex::Real flxw = cdtdx * (flxx(ivp0, UMZ) - flxx(iv, UMZ));
  const amrex::Real flxe = cdtdx * (flxx(ivp0, UEDEN) - flxx(iv, UEDEN));
  amrex::Real c = qa(iv, QGAMC);
  // Update passive variables
  for (int ipass = 0; ipass < NPASSIVE; ++ipass) {
    const int n = pmap.upassMap[ipass];
    const int nqp = pmap.qpassMap[ipass];
    amrex::Real compn = cdtdx * (flxx(ivp0, n) - flxx(iv, n));
    amrex::Real rr = qyp(iv, QRHO);
    amrex::Real rrnew = rr - flxrho;
    amrex::Real compo = rr * qyp(iv, nqp) - compn;
    qp(iv, nqp) = compo / rrnew;

    rr = qym(ivp1, QRHO);
    rrnew = rr - flxrho;
    compo = rr * qym(ivp1, nqp) - compn;
    qm(ivp1, nqp) = compo / rrnew;
  }
  // Update hydro vars

  const amrex::Real pggp = qint(ivp0, GDPRES);
  const amrex::Real pggm = qint(iv, GDPRES);
  const amrex::Real ugp = qint(ivp0, GDU);
  const amrex::Real ugm = qint(iv, GDU);

  const amrex::Real dAup = pggp * ugp - pggm * ugm;
  const amrex::Real pav = 0.5 * (pggp + pggm);
  const amrex::Real dAu = ugp - ugm;

  // QP

  // Convert to conservative
  const amrex::Real rrr = qyp(iv, QRHO);
  amrex::Real rur = qyp(iv, QU);
  amrex::Real rvr = qyp(iv, QV);
  amrex::Real rwr = qyp(iv, QW);
  const amrex::Real ekinr = 0.5 * rrr * (rur * rur + rvr * rvr + rwr * rwr);
  rur *= rrr;
  rvr *= rrr;
  rwr *= rrr;

  const amrex::Real rer = qyp(iv, QREINT) + ekinr;
  // Add transverse predictor
  const amrex::Real rrnewr = rrr - flxrho;
  const amrex::Real runewr = rur - flxu;
  const amrex::Real rvnewr = rvr - flxv;
  const amrex::Real rwnewr = rwr - flxw;
  const amrex::Real renewr = rer - flxe;

  // Convert back to primitive
  qp(iv, QRHO) = rrnewr;
  qp(iv, QU) = runewr / rrnewr;
  qp(iv, QV) = rvnewr / rrnewr;
  qp(iv, QW) = rwnewr / rrnewr;
  const amrex::Real rhoekinr =
    0.5 * (runewr * runewr + rvnewr * rvnewr + rwnewr * rwnewr) / rrnewr;
  amrex::Real pnewr = qyp(iv, QPRES) - cdtdx * (dAup + pav * dAu * (c - 1.0));
  qp(iv, QPRES) =
    amrex::max<amrex::Real>(pnewr, std::numeric_limits<amrex::Real>::min());
  qp(iv, QREINT) = renewr - rhoekinr;
  // QM

  // Conversion to Conservative
  const amrex::Real rrl = qym(ivp1, QRHO);
  amrex::Real rul = qym(ivp1, QU);
  amrex::Real rvl = qym(ivp1, QV);
  amrex::Real rwl = qym(ivp1, QW);
  const amrex::Real ekinl = 0.5 * rrl * (rul * rul + rvl * rvl + rwl * rwl);
  rul *= rrl;
  rvl *= rrl;
  rwl *= rrl;
  const amrex::Real rel = qym(ivp1, QREINT) + ekinl;

  // Transverse fluxes
  const amrex::Real rrnewl = rrl - flxrho;
  const amrex::Real runewl = rul - flxu;
  const amrex::Real rvnewl = rvl - flxv;
  const amrex::Real rwnewl = rwl - flxw;
  const amrex::Real renewl = rel - flxe;

  qm(ivp1, QRHO) = rrnewl;
  qm(ivp1, QU) = runewl / rrnewl;
  qm(ivp1, QV) = rvnewl / rrnewl;
  qm(ivp1, QW) = rwnewl / rrnewl;
  const amrex::Real rhoekinl =
    0.5 * (runewl * runewl + rvnewl * rvnewl + rwnewl * rwnewl) / rrnewl;

  qm(ivp1, QREINT) = renewl - rhoekinl;
  amrex::Real pnewl = qym(ivp1, QPRES) - cdtdx * (dAup + pav * dAu * (c - 1.0));
  qm(ivp1, QPRES) =
    amrex::max<amrex::Real>(pnewl, std::numeric_limits<amrex::Real>::min());
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
  const amrex::Real cdtdy,
  PassMap const& pmap)
{
  const int dir = 1;
  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  const amrex::IntVect ivp0{AMREX_D_DECL(i, j + 1, k)};
  const amrex::IntVect ivp1{AMREX_D_DECL(i + 1, j, k)};
  const amrex::Real flxrho = cdtdy * (flxy(ivp0, URHO) - flxy(iv, URHO));
  const amrex::Real flxu = cdtdy * (flxy(ivp0, UMX) - flxy(iv, UMX));
  const amrex::Real flxv = cdtdy * (flxy(ivp0, UMY) - flxy(iv, UMY));
  const amrex::Real flxw = cdtdy * (flxy(ivp0, UMZ) - flxy(iv, UMZ));
  const amrex::Real flxe = cdtdy * (flxy(ivp0, UEDEN) - flxy(iv, UEDEN));
  const amrex::Real c = qa(iv, QGAMC);
  // -  Update Passive Variables
  for (int ipass = 0; ipass < NPASSIVE; ++ipass) {
    const int n = pmap.upassMap[ipass];
    const int nqp = pmap.qpassMap[ipass];
    const amrex::Real compn = cdtdy * (flxy(ivp0, n) - flxy(iv, n));
    amrex::Real rr = qxp(iv, QRHO);
    amrex::Real rrnew = rr - flxrho;
    amrex::Real compo = rr * qxp(iv, nqp) - compn;
    qp(iv, nqp) = compo / rrnew;

    rr = qxm(ivp1, QRHO);
    rrnew = rr - flxrho;
    compo = rr * qxm(ivp1, nqp) - compn;
    qm(ivp1, nqp) = compo / rrnew;
  }

  //- Update hydro vars -
  // pressure divergence and such.
  const amrex::Real pggp = qint(ivp0, GDPRES);
  const amrex::Real pggm = qint(iv, GDPRES);
  const amrex::Real ugp = qint(ivp0, GDV);
  const amrex::Real ugm = qint(iv, GDV);

  const amrex::Real dup = pggp * ugp - pggm * ugm;
  const amrex::Real pav = 0.5 * (pggp + pggm);
  const amrex::Real du = ugp - ugm;

  // QP

  // Convert to conservative
  const amrex::Real rrr = qxp(iv, QRHO);
  amrex::Real rur = qxp(iv, QU);
  amrex::Real rvr = qxp(iv, QV);
  amrex::Real rwr = qxp(iv, QW);
  amrex::Real ekinr = 0.5 * rrr * (rur * rur + rvr * rvr + rwr * rwr);
  rur *= rrr;
  rvr *= rrr;
  rwr *= rrr;

  const amrex::Real rer = qxp(iv, QREINT) + ekinr;
  // Add transverse predictor
  const amrex::Real rrnewr = rrr - flxrho;
  const amrex::Real runewr = rur - flxu;
  const amrex::Real rvnewr = rvr - flxv;
  const amrex::Real rwnewr = rwr - flxw;
  const amrex::Real renewr = rer - flxe;

  // Convert back to primitive
  qp(iv, QRHO) = rrnewr; // + hdt*srcr;
  qp(iv, QU) = runewr / rrnewr;
  qp(iv, QV) = rvnewr / rrnewr;
  qp(iv, QW) = rwnewr / rrnewr;
  const amrex::Real rhoekinr =
    0.5 * (runewr * runewr + rvnewr * rvnewr + rwnewr * rwnewr) / rrnewr;

  const amrex::Real pnewr =
    qxp(iv, QPRES) - cdtdy * (dup + pav * du * (c - 1.0));
  qp(iv, QPRES) =
    amrex::max<amrex::Real>(pnewr, std::numeric_limits<amrex::Real>::min());
  qp(iv, QREINT) = renewr - rhoekinr;

  // QM

  // Conversion to Conservative
  const amrex::Real rrl = qxm(ivp1, QRHO);
  amrex::Real rul = qxm(ivp1, QU);
  amrex::Real rvl = qxm(ivp1, QV);
  amrex::Real rwl = qxm(ivp1, QW);
  amrex::Real ekinl = 0.5 * rrl * (rul * rul + rvl * rvl + rwl * rwl);
  rul *= rrl;
  rvl *= rrl;
  rwl *= rrl;
  const amrex::Real rel = qxm(ivp1, QREINT) + ekinl;

  // Transverse fluxes
  const amrex::Real rrnewl = rrl - flxrho;
  const amrex::Real runewl = rul - flxu;
  const amrex::Real rvnewl = rvl - flxv;
  const amrex::Real rwnewl = rwl - flxw;
  const amrex::Real renewl = rel - flxe;

  qm(ivp1, QRHO) = rrnewl;
  qm(ivp1, QU) = runewl / rrnewl;
  qm(ivp1, QV) = rvnewl / rrnewl;
  qm(ivp1, QW) = rwnewl / rrnewl;
  const amrex::Real rhoekinl =
    0.5 * (runewl * runewl + rvnewl * rvnewl + rwnewl * rwnewl) / rrnewl;

  const amrex::Real pnewl =
    qxm(ivp1, QPRES) - cdtdy * (dup + pav * du * (c - 1.0));
  qm(ivp1, QPRES) =
    amrex::max<amrex::Real>(pnewl, std::numeric_limits<amrex::Real>::min());
  qm(ivp1, QREINT) = renewl - rhoekinl;
}

// Transverse Correction for Predicted Z-states, using X-Flux
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
  const amrex::Real cdtdx,
  PassMap const& pmap)
{
  const int dir = 0;
  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  const amrex::Real flxrho = cdtdx * (flxx(i + 1, j, k, URHO) - flxx(iv, URHO));
  const amrex::Real flxu = cdtdx * (flxx(i + 1, j, k, UMX) - flxx(iv, UMX));
  const amrex::Real flxv = cdtdx * (flxx(i + 1, j, k, UMY) - flxx(iv, UMY));
  const amrex::Real flxw = cdtdx * (flxx(i + 1, j, k, UMZ) - flxx(iv, UMZ));
  const amrex::Real flxe = cdtdx * (flxx(i + 1, j, k, UEDEN) - flxx(iv, UEDEN));
  const amrex::Real c = qa(iv, QGAMC);
  // -  Update Passive Variables
  for (int ipass = 0; ipass < NPASSIVE; ++ipass) {
    const int n = pmap.upassMap[ipass];
    const int nqp = pmap.qpassMap[ipass];
    const amrex::Real compn = cdtdx * (flxx(i + 1, j, k, n) - flxx(iv, n));
    amrex::Real rr = qzp(iv, QRHO);
    amrex::Real rrnew = rr - flxrho;
    amrex::Real compo = rr * qzp(iv, nqp) - compn;
    qp(iv, nqp) = compo / rrnew;

    rr = qzm(i, j, k + 1, QRHO);
    rrnew = rr - flxrho;
    compo = rr * qzm(i, j, k + 1, nqp) - compn;
    qm(i, j, k + 1, nqp) = compo / rrnew;
  }
  //- Update hydro vars -

  const amrex::Real pggp = qint(i + 1, j, k, GDPRES);
  const amrex::Real pggm = qint(iv, GDPRES);
  const amrex::Real ugp = qint(i + 1, j, k, GDU);
  const amrex::Real ugm = qint(iv, GDU);

  const amrex::Real dAup = pggp * ugp - pggm * ugm;
  const amrex::Real pav = 0.5 * (pggp + pggm);
  const amrex::Real dAu = ugp - ugm;

  // QP

  // Convert to conservative
  const amrex::Real rrr = qzp(iv, QRHO);
  amrex::Real rur = qzp(iv, QU);
  amrex::Real rvr = qzp(iv, QV);
  amrex::Real rwr = qzp(iv, QW);
  const amrex::Real ekinr = 0.5 * rrr * (rur * rur + rvr * rvr + rwr * rwr);
  rur *= rrr;
  rvr *= rrr;
  rwr *= rrr;

  const amrex::Real rer = qzp(iv, QREINT) + ekinr;
  // Add transverse predictor
  const amrex::Real rrnewr = rrr - flxrho;
  const amrex::Real runewr = rur - flxu;
  const amrex::Real rvnewr = rvr - flxv;
  const amrex::Real rwnewr = rwr - flxw;
  const amrex::Real renewr = rer - flxe;

  // Convert back to primitive
  qp(iv, QRHO) = rrnewr;
  qp(iv, QU) = runewr / rrnewr;
  qp(iv, QV) = rvnewr / rrnewr;
  qp(iv, QW) = rwnewr / rrnewr;
  const amrex::Real rhoekinr =
    0.5 * (runewr * runewr + rvnewr * rvnewr + rwnewr * rwnewr) / rrnewr;
  const amrex::Real pnewr =
    qzp(iv, QPRES) - cdtdx * (dAup + pav * dAu * (c - 1.0));
  qp(iv, QPRES) =
    amrex::max<amrex::Real>(pnewr, std::numeric_limits<amrex::Real>::min());
  qp(iv, QREINT) = renewr - rhoekinr;
  // QM

  // Conversion to Conservative
  const amrex::Real rrl = qzm(i, j, k + 1, QRHO);
  amrex::Real rul = qzm(i, j, k + 1, QU);
  amrex::Real rvl = qzm(i, j, k + 1, QV);
  amrex::Real rwl = qzm(i, j, k + 1, QW);
  const amrex::Real ekinl = 0.5 * rrl * (rul * rul + rvl * rvl + rwl * rwl);
  rul *= rrl;
  rvl *= rrl;
  rwl *= rrl;
  const amrex::Real rel = qzm(i, j, k + 1, QREINT) + ekinl;

  // Transverse fluxes
  const amrex::Real rrnewl = rrl - flxrho;
  const amrex::Real runewl = rul - flxu;
  const amrex::Real rvnewl = rvl - flxv;
  const amrex::Real rwnewl = rwl - flxw;
  const amrex::Real renewl = rel - flxe;

  qm(i, j, k + 1, QRHO) = rrnewl;
  qm(i, j, k + 1, QU) = runewl / rrnewl;
  qm(i, j, k + 1, QV) = rvnewl / rrnewl;
  qm(i, j, k + 1, QW) = rwnewl / rrnewl;
  const amrex::Real rhoekinl =
    0.5 * (runewl * runewl + rvnewl * rvnewl + rwnewl * rwnewl) / rrnewl;

  qm(i, j, k + 1, QREINT) = renewl - rhoekinl;
  const amrex::Real pnewl =
    qzm(i, j, k + 1, QPRES) - cdtdx * (dAup + pav * dAu * (c - 1.0));
  qm(i, j, k + 1, QPRES) =
    amrex::max<amrex::Real>(pnewl, std::numeric_limits<amrex::Real>::min());
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
  const amrex::Real cdtdy,
  PassMap const& pmap)
{
  const int dir = 1;
  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  const amrex::Real flxrho = cdtdy * (flxy(i, j + 1, k, URHO) - flxy(iv, URHO));
  const amrex::Real flxu = cdtdy * (flxy(i, j + 1, k, UMX) - flxy(iv, UMX));
  const amrex::Real flxv = cdtdy * (flxy(i, j + 1, k, UMY) - flxy(iv, UMY));
  const amrex::Real flxw = cdtdy * (flxy(i, j + 1, k, UMZ) - flxy(iv, UMZ));
  const amrex::Real flxe = cdtdy * (flxy(i, j + 1, k, UEDEN) - flxy(iv, UEDEN));
  const amrex::Real c = qa(iv, QGAMC);
  // -  Update Passive Variables
  for (int ipass = 0; ipass < NPASSIVE; ++ipass) {
    const int n = pmap.upassMap[ipass];
    const int nqp = pmap.qpassMap[ipass];
    const amrex::Real compn = cdtdy * (flxy(i, j + 1, k, n) - flxy(iv, n));
    amrex::Real rr = qzp(iv, QRHO);
    amrex::Real rrnew = rr - flxrho;
    amrex::Real compo = rr * qzp(iv, nqp) - compn;
    qp(iv, nqp) = compo / rrnew;

    rr = qzm(i, j, k + 1, QRHO);
    rrnew = rr - flxrho;
    compo = rr * qzm(i, j, k + 1, nqp) - compn;
    qm(i, j, k + 1, nqp) = compo / rrnew;
  }
  //- Update hydro vars -

  const amrex::Real pggp = qint(i, j + 1, k, GDPRES);
  const amrex::Real pggm = qint(iv, GDPRES);
  const amrex::Real ugp = qint(i, j + 1, k, GDV);
  const amrex::Real ugm = qint(iv, GDV);

  const amrex::Real dAup = pggp * ugp - pggm * ugm;
  const amrex::Real pav = 0.5 * (pggp + pggm);
  const amrex::Real dAu = ugp - ugm;

  // QP

  // Convert to conservative
  const amrex::Real rrr = qzp(iv, QRHO);
  amrex::Real rur = qzp(iv, QU);
  amrex::Real rvr = qzp(iv, QV);
  amrex::Real rwr = qzp(iv, QW);
  const amrex::Real ekinr = 0.5 * rrr * (rur * rur + rvr * rvr + rwr * rwr);
  rur *= rrr;
  rvr *= rrr;
  rwr *= rrr;

  const amrex::Real rer = qzp(iv, QREINT) + ekinr;
  // Add transverse predictor
  const amrex::Real rrnewr = rrr - flxrho;
  const amrex::Real runewr = rur - flxu;
  const amrex::Real rvnewr = rvr - flxv;
  const amrex::Real rwnewr = rwr - flxw;
  const amrex::Real renewr = rer - flxe;

  // Convert back to primitive
  qp(iv, QRHO) = rrnewr;
  qp(iv, QU) = runewr / rrnewr;
  qp(iv, QV) = rvnewr / rrnewr;
  qp(iv, QW) = rwnewr / rrnewr;
  const amrex::Real rhoekinr =
    0.5 * (runewr * runewr + rvnewr * rvnewr + rwnewr * rwnewr) / rrnewr;
  const amrex::Real pnewr =
    qzp(iv, QPRES) - cdtdy * (dAup + pav * dAu * (c - 1.0));
  qp(iv, QPRES) =
    amrex::max<amrex::Real>(pnewr, std::numeric_limits<amrex::Real>::min());
  qp(iv, QREINT) = renewr - rhoekinr;
  // QM

  // Conversion to Conservative
  const amrex::Real rrl = qzm(i, j, k + 1, QRHO);
  amrex::Real rul = qzm(i, j, k + 1, QU);
  amrex::Real rvl = qzm(i, j, k + 1, QV);
  amrex::Real rwl = qzm(i, j, k + 1, QW);
  const amrex::Real ekinl = 0.5 * rrl * (rul * rul + rvl * rvl + rwl * rwl);
  rul *= rrl;
  rvl *= rrl;
  rwl *= rrl;
  const amrex::Real rel = qzm(i, j, k + 1, QREINT) + ekinl;

  // Transverse fluxes
  const amrex::Real rrnewl = rrl - flxrho;
  const amrex::Real runewl = rul - flxu;
  const amrex::Real rvnewl = rvl - flxv;
  const amrex::Real rwnewl = rwl - flxw;
  const amrex::Real renewl = rel - flxe;

  qm(i, j, k + 1, QRHO) = rrnewl;
  qm(i, j, k + 1, QU) = runewl / rrnewl;
  qm(i, j, k + 1, QV) = rvnewl / rrnewl;
  qm(i, j, k + 1, QW) = rwnewl / rrnewl;
  const amrex::Real rhoekinl =
    0.5 * (runewl * runewl + rvnewl * rvnewl + rwnewl * rwnewl) / rrnewl;

  qm(i, j, k + 1, QREINT) = renewl - rhoekinl;
  const amrex::Real pnewl =
    qzm(i, j, k + 1, QPRES) - cdtdy * (dAup + pav * dAu * (c - 1.0));
  qm(i, j, k + 1, QPRES) =
    amrex::max<amrex::Real>(pnewl, std::numeric_limits<amrex::Real>::min());
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
  const amrex::Real cdtdz,
  PassMap const& pmap)
{
  const int dir = 1;
  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  const amrex::IntVect ivp0{AMREX_D_DECL(i, j, k + 1)};
  const amrex::IntVect ivp1{AMREX_D_DECL(i + 1, j, k)};
  const amrex::Real flxrho = cdtdz * (flxz(ivp0, URHO) - flxz(iv, URHO));
  const amrex::Real flxu = cdtdz * (flxz(ivp0, UMX) - flxz(iv, UMX));
  const amrex::Real flxv = cdtdz * (flxz(ivp0, UMY) - flxz(iv, UMY));
  const amrex::Real flxw = cdtdz * (flxz(ivp0, UMZ) - flxz(iv, UMZ));
  const amrex::Real flxe = cdtdz * (flxz(ivp0, UEDEN) - flxz(iv, UEDEN));
  const amrex::Real c = qa(iv, QGAMC);
  // -  Update Passive Variables
  for (int ipass = 0; ipass < NPASSIVE; ++ipass) {
    const int n = pmap.upassMap[ipass];
    const int nqp = pmap.qpassMap[ipass];
    const amrex::Real compn = cdtdz * (flxz(ivp0, n) - flxz(iv, n));
    amrex::Real rr = qxp(iv, QRHO);
    amrex::Real rrnew = rr - flxrho;
    amrex::Real compo = rr * qxp(iv, nqp) - compn;
    qp(iv, nqp) = compo / rrnew;

    rr = qxm(ivp1, QRHO);
    rrnew = rr - flxrho;
    compo = rr * qxm(ivp1, nqp) - compn;
    qm(ivp1, nqp) = compo / rrnew;
  }

  //- Update hydro vars -
  // pressure divergence and such.
  const amrex::Real pggp = qint(ivp0, GDPRES);
  const amrex::Real pggm = qint(iv, GDPRES);
  const amrex::Real ugp = qint(ivp0, GDW);
  const amrex::Real ugm = qint(iv, GDW);

  const amrex::Real dup = pggp * ugp - pggm * ugm;
  const amrex::Real pav = 0.5 * (pggp + pggm);
  const amrex::Real du = ugp - ugm;

  // QP

  // Convert to conservative
  const amrex::Real rrr = qxp(iv, QRHO);
  amrex::Real rur = qxp(iv, QU);
  amrex::Real rvr = qxp(iv, QV);
  amrex::Real rwr = qxp(iv, QW);
  const amrex::Real ekinr = 0.5 * rrr * (rur * rur + rvr * rvr + rwr * rwr);
  rur *= rrr;
  rvr *= rrr;
  rwr *= rrr;

  const amrex::Real rer = qxp(iv, QREINT) + ekinr;
  // Add transverse predictor
  const amrex::Real rrnewr = rrr - flxrho;
  const amrex::Real runewr = rur - flxu;
  const amrex::Real rvnewr = rvr - flxv;
  const amrex::Real rwnewr = rwr - flxw;
  const amrex::Real renewr = rer - flxe;

  // Convert back to primitive
  qp(iv, QRHO) = rrnewr; // + hdt*srcr;
  qp(iv, QU) = runewr / rrnewr;
  qp(iv, QV) = rvnewr / rrnewr;
  qp(iv, QW) = rwnewr / rrnewr;
  const amrex::Real rhoekinr =
    0.5 * (runewr * runewr + rvnewr * rvnewr + rwnewr * rwnewr) / rrnewr;

  const amrex::Real pnewr =
    qxp(iv, QPRES) - cdtdz * (dup + pav * du * (c - 1.0));
  qp(iv, QPRES) =
    amrex::max<amrex::Real>(pnewr, std::numeric_limits<amrex::Real>::min());
  qp(iv, QREINT) = renewr - rhoekinr;

  // QM

  // Conversion to Conservative
  const amrex::Real rrl = qxm(ivp1, QRHO);
  amrex::Real rul = qxm(ivp1, QU);
  amrex::Real rvl = qxm(ivp1, QV);
  amrex::Real rwl = qxm(ivp1, QW);
  const amrex::Real ekinl = 0.5 * rrl * (rul * rul + rvl * rvl + rwl * rwl);
  rul *= rrl;
  rvl *= rrl;
  rwl *= rrl;
  const amrex::Real rel = qxm(ivp1, QREINT) + ekinl;

  // Transverse fluxes
  const amrex::Real rrnewl = rrl - flxrho;
  const amrex::Real runewl = rul - flxu;
  const amrex::Real rvnewl = rvl - flxv;
  const amrex::Real rwnewl = rwl - flxw;
  const amrex::Real renewl = rel - flxe;

  qm(ivp1, QRHO) = rrnewl;
  qm(ivp1, QU) = runewl / rrnewl;
  qm(ivp1, QV) = rvnewl / rrnewl;
  qm(ivp1, QW) = rwnewl / rrnewl;
  const amrex::Real rhoekinl =
    0.5 * (runewl * runewl + rvnewl * rvnewl + rwnewl * rwnewl) / rrnewl;

  const amrex::Real pnewl =
    qxm(ivp1, QPRES) - cdtdz * (dup + pav * du * (c - 1.0));
  qm(ivp1, QPRES) =
    amrex::max<amrex::Real>(pnewl, std::numeric_limits<amrex::Real>::min());
  qm(ivp1, QREINT) = renewl - rhoekinl;
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
  const amrex::Real cdtdz,
  PassMap const& pmap)
{
  const int dir = 2;
  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  const amrex::Real flxrho = cdtdz * (flxz(i, j, k + 1, URHO) - flxz(iv, URHO));
  const amrex::Real flxu = cdtdz * (flxz(i, j, k + 1, UMX) - flxz(iv, UMX));
  const amrex::Real flxv = cdtdz * (flxz(i, j, k + 1, UMY) - flxz(iv, UMY));
  const amrex::Real flxw = cdtdz * (flxz(i, j, k + 1, UMZ) - flxz(iv, UMZ));
  const amrex::Real flxe = cdtdz * (flxz(i, j, k + 1, UEDEN) - flxz(iv, UEDEN));
  const amrex::Real c = qa(iv, QGAMC);
  // Update passive variables
  for (int ipass = 0; ipass < NPASSIVE; ++ipass) {
    const int n = pmap.upassMap[ipass];
    const int nqp = pmap.qpassMap[ipass];
    const amrex::Real compn = cdtdz * (flxz(i, j, k + 1, n) - flxz(iv, n));
    amrex::Real rr = qyp(iv, QRHO);
    amrex::Real rrnew = rr - flxrho;
    amrex::Real compo = rr * qyp(iv, nqp) - compn;
    qp(iv, nqp) = compo / rrnew;

    rr = qym(i, j + 1, k, QRHO);
    rrnew = rr - flxrho;
    compo = rr * qym(i, j + 1, k, nqp) - compn;
    qm(i, j + 1, k, nqp) = compo / rrnew;
  }

  // Update hydro vars
  const amrex::Real pggp = qint(i, j, k + 1, GDPRES);
  const amrex::Real pggm = qint(iv, GDPRES);
  const amrex::Real ugp = qint(i, j, k + 1, GDW);
  const amrex::Real ugm = qint(iv, GDW);

  const amrex::Real dAup = pggp * ugp - pggm * ugm;
  const amrex::Real pav = 0.5 * (pggp + pggm);
  const amrex::Real dAu = ugp - ugm;

  // QP

  // Convert to conservative
  const amrex::Real rrr = qyp(iv, QRHO);
  amrex::Real rur = qyp(iv, QU);
  amrex::Real rvr = qyp(iv, QV);
  amrex::Real rwr = qyp(iv, QW);
  const amrex::Real ekinr = 0.5 * rrr * (rur * rur + rvr * rvr + rwr * rwr);
  rur *= rrr;
  rvr *= rrr;
  rwr *= rrr;
  const amrex::Real rer = qyp(iv, QREINT) + ekinr;
  // Add transverse predictor
  const amrex::Real rrnewr = rrr - flxrho;
  const amrex::Real runewr = rur - flxu;
  const amrex::Real rvnewr = rvr - flxv;
  const amrex::Real rwnewr = rwr - flxw;
  const amrex::Real renewr = rer - flxe;

  // Convert back to primitive
  qp(iv, QRHO) = rrnewr;
  qp(iv, QU) = runewr / rrnewr;
  qp(iv, QV) = rvnewr / rrnewr;
  qp(iv, QW) = rwnewr / rrnewr;
  const amrex::Real rhoekinr =
    0.5 * (runewr * runewr + rvnewr * rvnewr + rwnewr * rwnewr) / rrnewr;
  const amrex::Real pnewr =
    qyp(iv, QPRES) - cdtdz * (dAup + pav * dAu * (c - 1.0));
  qp(iv, QPRES) =
    amrex::max<amrex::Real>(pnewr, std::numeric_limits<amrex::Real>::min());
  qp(iv, QREINT) = renewr - rhoekinr;

  // QM

  // Conversion to Conservative
  const amrex::Real rrl = qym(i, j + 1, k, QRHO);
  amrex::Real rul = qym(i, j + 1, k, QU);
  amrex::Real rvl = qym(i, j + 1, k, QV);
  amrex::Real rwl = qym(i, j + 1, k, QW);
  const amrex::Real ekinl = 0.5 * rrl * (rul * rul + rvl * rvl + rwl * rwl);
  rul *= rrl;
  rvl *= rrl;
  rwl *= rrl;
  const amrex::Real rel = qym(i, j + 1, k, QREINT) + ekinl;

  // Transverse fluxes
  const amrex::Real rrnewl = rrl - flxrho;
  const amrex::Real runewl = rul - flxu;
  const amrex::Real rvnewl = rvl - flxv;
  const amrex::Real rwnewl = rwl - flxw;
  const amrex::Real renewl = rel - flxe;

  qm(i, j + 1, k, QRHO) = rrnewl;
  qm(i, j + 1, k, QU) = runewl / rrnewl;
  qm(i, j + 1, k, QV) = rvnewl / rrnewl;
  qm(i, j + 1, k, QW) = rwnewl / rrnewl;
  const amrex::Real rhoekinl =
    0.5 * (runewl * runewl + rvnewl * rvnewl + rwnewl * rwnewl) / rrnewl;

  qm(i, j + 1, k, QREINT) = renewl - rhoekinl;
  const amrex::Real pnewl =
    qym(i, j + 1, k, QPRES) - cdtdz * (dAup + pav * dAu * (c - 1.0));
  qm(i, j + 1, k, QPRES) =
    amrex::max<amrex::Real>(pnewl, std::numeric_limits<amrex::Real>::min());
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
  const amrex::Real cdtdy,
  PassMap const& pmap)
{
  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  const amrex::IntVect ivp0{AMREX_D_DECL(i + 1, j, k)};
  const amrex::IntVect ivp1{AMREX_D_DECL(i, j + 1, k)};
  const amrex::IntVect ivp2{AMREX_D_DECL(i, j, k + 1)};
  const amrex::Real flxrho = cdtdx * (flxx(ivp0, URHO) - flxx(iv, URHO)) +
                             cdtdy * (flxy(ivp1, URHO) - flxy(iv, URHO));
  const amrex::Real flxu = cdtdx * (flxx(ivp0, UMX) - flxx(iv, UMX)) +
                           cdtdy * (flxy(ivp1, UMX) - flxy(iv, UMX));
  const amrex::Real flxv = cdtdx * (flxx(ivp0, UMY) - flxx(iv, UMY)) +
                           cdtdy * (flxy(ivp1, UMY) - flxy(iv, UMY));
  const amrex::Real flxw = cdtdx * (flxx(ivp0, UMZ) - flxx(iv, UMZ)) +
                           cdtdy * (flxy(ivp1, UMZ) - flxy(iv, UMZ));
  const amrex::Real flxe = cdtdx * (flxx(ivp0, UEDEN) - flxx(iv, UEDEN)) +
                           cdtdy * (flxy(ivp1, UEDEN) - flxy(iv, UEDEN));
  const amrex::Real c = qa(iv, QGAMC);
  const amrex::Real srcrho = srcq(iv, QRHO);
  const amrex::Real srcu = srcq(iv, QU);
  const amrex::Real srcv = srcq(iv, QV);
  const amrex::Real srcw = srcq(iv, QW);
  const amrex::Real srce = srcq(iv, QREINT);
  const amrex::Real srcp = srcq(iv, QPRES);

  // Update passive variables
  const amrex::Real rrr = qzp(iv, QRHO);
  const amrex::Real rrl = qzm(ivp2, QRHO);
  const amrex::Real rrnewl = rrl - flxrho;
  const amrex::Real rrnewr = rrr - flxrho;
  for (int ipass = 0; ipass < NPASSIVE; ++ipass) {
    const int n = pmap.upassMap[ipass];
    const int nqp = pmap.qpassMap[ipass];
    const amrex::Real srcpass = srcq(iv, nqp);
    const amrex::Real compn = cdtdx * (flxx(ivp0, n) - flxx(iv, n)) +
                              cdtdy * (flxy(ivp1, n) - flxy(iv, n));
    // qp
    amrex::Real compo = rrr * qzp(iv, nqp) - compn;
    qp(iv, nqp) = compo / rrnewr + hdt * srcpass;
    // qm
    compo = rrl * qzm(ivp2, nqp) - compn;
    qm(ivp2, nqp) = compo / rrnewl + hdt * srcpass;
  }

  // Update hydro vars
  const amrex::Real pggpx = qx(ivp0, GDPRES);
  const amrex::Real pggmx = qx(iv, GDPRES);
  const amrex::Real ugpx = qx(ivp0, GDU);
  const amrex::Real ugmx = qx(iv, GDU);

  const amrex::Real dAupx = pggpx * ugpx - pggmx * ugmx;
  const amrex::Real pavx = 0.5 * (pggpx + pggmx);
  const amrex::Real dAux = ugpx - ugmx;

  const amrex::Real pggpy = qy(ivp1, GDPRES);
  const amrex::Real pggmy = qy(iv, GDPRES);
  const amrex::Real ugpy = qy(ivp1, GDV);
  const amrex::Real ugmy = qy(iv, GDV);

  const amrex::Real dAupy = pggpy * ugpy - pggmy * ugmy;
  const amrex::Real pavy = 0.5 * (pggpy + pggmy);
  const amrex::Real dAuy = ugpy - ugmy;
  const amrex::Real pxnew = cdtdx * (dAupx + pavx * dAux * (c - 1.0));
  const amrex::Real pynew = cdtdy * (dAupy + pavy * dAuy * (c - 1.0));

  // qp state
  const amrex::Real rur = rrr * qzp(iv, QU);
  const amrex::Real rvr = rrr * qzp(iv, QV);
  const amrex::Real rwr = rrr * qzp(iv, QW);
  const amrex::Real ekinr = 0.5 * (rur * rur + rvr * rvr + rwr * rwr) / rrr;
  const amrex::Real rer = qzp(iv, QREINT) + ekinr;

  const amrex::Real runewr = rur - flxu;
  const amrex::Real rvnewr = rvr - flxv;
  const amrex::Real rwnewr = rwr - flxw;
  const amrex::Real renewr = rer - flxe;

  qp(iv, QRHO) = rrnewr + hdt * srcrho;
  qp(iv, QU) = runewr / rrnewr + hdt * srcu;
  qp(iv, QV) = rvnewr / rrnewr + hdt * srcv;
  qp(iv, QW) = rwnewr / rrnewr + hdt * srcw;
  const amrex::Real rhoekinr =
    0.5 * (runewr * runewr + rvnewr * rvnewr + rwnewr * rwnewr) / rrnewr;
  qp(iv, QREINT) = renewr - rhoekinr + hdt * srce;
  amrex::Real temppres = qzp(iv, QPRES) - pxnew - pynew + hdt * srcp;
  qp(iv, QPRES) =
    amrex::max<amrex::Real>(temppres, std::numeric_limits<amrex::Real>::min());

  // qm state
  const amrex::Real rul = rrl * qzm(ivp2, QU);
  const amrex::Real rvl = rrl * qzm(ivp2, QV);
  const amrex::Real rwl = rrl * qzm(ivp2, QW);
  const amrex::Real ekinl = 0.5 * (rul * rul + rvl * rvl + rwl * rwl) / rrl;
  const amrex::Real rel = qzm(ivp2, QREINT) + ekinl;

  const amrex::Real runewl = rul - flxu;
  const amrex::Real rvnewl = rvl - flxv;
  const amrex::Real rwnewl = rwl - flxw;
  const amrex::Real renewl = rel - flxe;

  qm(ivp2, QRHO) = rrnewl + hdt * srcrho;
  qm(ivp2, QU) = runewl / rrnewl + hdt * srcu;
  qm(ivp2, QV) = rvnewl / rrnewl + hdt * srcv;
  qm(ivp2, QW) = rwnewl / rrnewl + hdt * srcw;
  const amrex::Real rhoekinl =
    0.5 * (runewl * runewl + rvnewl * rvnewl + rwnewl * rwnewl) / rrnewl;
  qm(ivp2, QREINT) = renewl - rhoekinl + hdt * srce;
  temppres = qzm(ivp2, QPRES) - pxnew - pynew + hdt * srcp;
  qm(ivp2, QPRES) =
    amrex::max<amrex::Real>(temppres, std::numeric_limits<amrex::Real>::min());
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
  const amrex::Real cdtdz,
  PassMap const& pmap)
{
  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  const amrex::IntVect ivp0{AMREX_D_DECL(i + 1, j, k)};
  const amrex::IntVect ivp1{AMREX_D_DECL(i, j, k + 1)};
  const amrex::IntVect ivp2{AMREX_D_DECL(i, j + 1, k)};
  const amrex::Real flxrho = cdtdx * (flxx(ivp0, URHO) - flxx(iv, URHO)) +
                             cdtdz * (flxz(ivp1, URHO) - flxz(iv, URHO));
  const amrex::Real flxu = cdtdx * (flxx(ivp0, UMX) - flxx(iv, UMX)) +
                           cdtdz * (flxz(ivp1, UMX) - flxz(iv, UMX));
  const amrex::Real flxv = cdtdx * (flxx(ivp0, UMY) - flxx(iv, UMY)) +
                           cdtdz * (flxz(ivp1, UMY) - flxz(iv, UMY));
  const amrex::Real flxw = cdtdx * (flxx(ivp0, UMZ) - flxx(iv, UMZ)) +
                           cdtdz * (flxz(ivp1, UMZ) - flxz(iv, UMZ));
  const amrex::Real flxe = cdtdx * (flxx(ivp0, UEDEN) - flxx(iv, UEDEN)) +
                           cdtdz * (flxz(ivp1, UEDEN) - flxz(iv, UEDEN));
  const amrex::Real c = qa(iv, QGAMC);
  const amrex::Real srcrho = srcq(iv, QRHO);
  const amrex::Real srcu = srcq(iv, QU);
  const amrex::Real srcv = srcq(iv, QV);
  const amrex::Real srcw = srcq(iv, QW);
  const amrex::Real srce = srcq(iv, QREINT);
  const amrex::Real srcp = srcq(iv, QPRES);

  // Update passive variables
  const amrex::Real rrr = qyp(iv, QRHO);
  const amrex::Real rrl = qym(ivp2, QRHO);
  const amrex::Real rrnewl = rrl - flxrho;
  const amrex::Real rrnewr = rrr - flxrho;
  for (int ipass = 0; ipass < NPASSIVE; ++ipass) {
    const int n = pmap.upassMap[ipass];
    const int nqp = pmap.qpassMap[ipass];
    const amrex::Real srcpass = srcq(iv, nqp);
    const amrex::Real compn = cdtdx * (flxx(ivp0, n) - flxx(iv, n)) +
                              cdtdz * (flxz(ivp1, n) - flxz(iv, n));
    // qp
    amrex::Real compo = rrr * qyp(iv, nqp) - compn;
    qp(iv, nqp) = compo / rrnewr + hdt * srcpass;
    // qm
    compo = rrl * qym(ivp2, nqp) - compn;
    qm(ivp2, nqp) = compo / rrnewl + hdt * srcpass;
  }

  // Update hydro vars
  const amrex::Real pggpx = qx(ivp0, GDPRES);
  const amrex::Real pggmx = qx(iv, GDPRES);
  const amrex::Real ugpx = qx(ivp0, GDU);
  const amrex::Real ugmx = qx(iv, GDU);

  const amrex::Real dAupx = pggpx * ugpx - pggmx * ugmx;
  const amrex::Real pavx = 0.5 * (pggpx + pggmx);
  const amrex::Real dAux = ugpx - ugmx;

  const amrex::Real pggpz = qz(ivp1, GDPRES);
  const amrex::Real pggmz = qz(iv, GDPRES);
  const amrex::Real ugpz = qz(ivp1, GDW);
  const amrex::Real ugmz = qz(iv, GDW);

  const amrex::Real dAupz = pggpz * ugpz - pggmz * ugmz;
  const amrex::Real pavz = 0.5 * (pggpz + pggmz);
  const amrex::Real dAuz = ugpz - ugmz;
  const amrex::Real pxnew = cdtdx * (dAupx + pavx * dAux * (c - 1.0));
  const amrex::Real pznew = cdtdz * (dAupz + pavz * dAuz * (c - 1.0));

  // qp state
  const amrex::Real rur = rrr * qyp(iv, QU);
  const amrex::Real rvr = rrr * qyp(iv, QV);
  const amrex::Real rwr = rrr * qyp(iv, QW);
  const amrex::Real ekinr = 0.5 * (rur * rur + rvr * rvr + rwr * rwr) / rrr;
  const amrex::Real rer = qyp(iv, QREINT) + ekinr;

  const amrex::Real runewr = rur - flxu;
  const amrex::Real rvnewr = rvr - flxv;
  const amrex::Real rwnewr = rwr - flxw;
  const amrex::Real renewr = rer - flxe;

  qp(iv, QRHO) = rrnewr + hdt * srcrho;
  qp(iv, QU) = runewr / rrnewr + hdt * srcu;
  qp(iv, QV) = rvnewr / rrnewr + hdt * srcv;
  qp(iv, QW) = rwnewr / rrnewr + hdt * srcw;
  const amrex::Real rhoekinr =
    0.5 * (runewr * runewr + rvnewr * rvnewr + rwnewr * rwnewr) / rrnewr;
  qp(iv, QREINT) = renewr - rhoekinr + hdt * srce;
  amrex::Real temppres = qyp(iv, QPRES) - pxnew - pznew + hdt * srcp;
  qp(iv, QPRES) =
    amrex::max<amrex::Real>(temppres, std::numeric_limits<amrex::Real>::min());

  // qm state
  const amrex::Real rul = rrl * qym(ivp2, QU);
  const amrex::Real rvl = rrl * qym(ivp2, QV);
  const amrex::Real rwl = rrl * qym(ivp2, QW);
  const amrex::Real ekinl = 0.5 * (rul * rul + rvl * rvl + rwl * rwl) / rrl;
  const amrex::Real rel = qym(ivp2, QREINT) + ekinl;

  const amrex::Real runewl = rul - flxu;
  const amrex::Real rvnewl = rvl - flxv;
  const amrex::Real rwnewl = rwl - flxw;
  const amrex::Real renewl = rel - flxe;

  qm(ivp2, QRHO) = rrnewl + hdt * srcrho;
  qm(ivp2, QU) = runewl / rrnewl + hdt * srcu;
  qm(ivp2, QV) = rvnewl / rrnewl + hdt * srcv;
  qm(ivp2, QW) = rwnewl / rrnewl + hdt * srcw;
  const amrex::Real rhoekinl =
    0.5 * (runewl * runewl + rvnewl * rvnewl + rwnewl * rwnewl) / rrnewl;
  qm(ivp2, QREINT) = renewl - rhoekinl + hdt * srce;
  temppres = qym(ivp2, QPRES) - pxnew - pznew + hdt * srcp;
  qm(ivp2, QPRES) =
    amrex::max<amrex::Real>(temppres, std::numeric_limits<amrex::Real>::min());
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
  const amrex::Real cdtdz,
  PassMap const& pmap)
{
  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  const amrex::IntVect ivp0{AMREX_D_DECL(i, j + 1, k)};
  const amrex::IntVect ivp1{AMREX_D_DECL(i, j, k + 1)};
  const amrex::IntVect ivp2{AMREX_D_DECL(i + 1, j, k)};
  const amrex::Real flxrho = cdtdy * (flxy(ivp0, URHO) - flxy(iv, URHO)) +
                             cdtdz * (flxz(ivp1, URHO) - flxz(iv, URHO));
  const amrex::Real flxu = cdtdy * (flxy(ivp0, UMX) - flxy(iv, UMX)) +
                           cdtdz * (flxz(ivp1, UMX) - flxz(iv, UMX));
  const amrex::Real flxv = cdtdy * (flxy(ivp0, UMY) - flxy(iv, UMY)) +
                           cdtdz * (flxz(ivp1, UMY) - flxz(iv, UMY));
  const amrex::Real flxw = cdtdy * (flxy(ivp0, UMZ) - flxy(iv, UMZ)) +
                           cdtdz * (flxz(ivp1, UMZ) - flxz(iv, UMZ));
  const amrex::Real flxe = cdtdy * (flxy(ivp0, UEDEN) - flxy(iv, UEDEN)) +
                           cdtdz * (flxz(ivp1, UEDEN) - flxz(iv, UEDEN));
  const amrex::Real c = qa(iv, QGAMC);
  const amrex::Real srcrho = srcq(iv, QRHO);
  const amrex::Real srcu = srcq(iv, QU);
  const amrex::Real srcv = srcq(iv, QV);
  const amrex::Real srcw = srcq(iv, QW);
  const amrex::Real srce = srcq(iv, QREINT);
  const amrex::Real srcp = srcq(iv, QPRES);

  // Update passive variables
  const amrex::Real rrr = qxp(iv, QRHO);
  const amrex::Real rrl = qxm(ivp2, QRHO);
  const amrex::Real rrnewl = rrl - flxrho;
  const amrex::Real rrnewr = rrr - flxrho;

  for (int ipass = 0; ipass < NPASSIVE; ++ipass) {
    const int n = pmap.upassMap[ipass];
    const int nqp = pmap.qpassMap[ipass];
    const amrex::Real srcpass = srcq(iv, nqp);
    const amrex::Real compn = cdtdy * (flxy(ivp0, n) - flxy(iv, n)) +
                              cdtdz * (flxz(ivp1, n) - flxz(iv, n));
    // qp
    amrex::Real compo = rrr * qxp(iv, nqp) - compn;
    qp(iv, nqp) = compo / rrnewr + hdt * srcpass;
    // qm
    compo = rrl * qxm(ivp2, nqp) - compn;
    qm(ivp2, nqp) = compo / rrnewl + hdt * srcpass;
  }

  // Update hydro vars
  const amrex::Real pggpz = qz(ivp1, GDPRES);
  const amrex::Real pggmz = qz(iv, GDPRES);
  const amrex::Real ugpz = qz(ivp1, GDW);
  const amrex::Real ugmz = qz(iv, GDW);

  const amrex::Real dAupz = pggpz * ugpz - pggmz * ugmz;
  const amrex::Real pavz = 0.5 * (pggpz + pggmz);
  const amrex::Real dAuz = ugpz - ugmz;

  const amrex::Real pggpy = qy(ivp0, GDPRES);
  const amrex::Real pggmy = qy(iv, GDPRES);
  const amrex::Real ugpy = qy(ivp0, GDV);
  const amrex::Real ugmy = qy(iv, GDV);

  const amrex::Real dAupy = pggpy * ugpy - pggmy * ugmy;
  const amrex::Real pavy = 0.5 * (pggpy + pggmy);
  const amrex::Real dAuy = ugpy - ugmy;
  const amrex::Real pznew = cdtdz * (dAupz + pavz * dAuz * (c - 1.0));
  const amrex::Real pynew = cdtdy * (dAupy + pavy * dAuy * (c - 1.0));

  // qp state
  const amrex::Real rur = rrr * qxp(iv, QU);
  const amrex::Real rvr = rrr * qxp(iv, QV);
  const amrex::Real rwr = rrr * qxp(iv, QW);
  const amrex::Real ekinr = 0.5 * (rur * rur + rvr * rvr + rwr * rwr) / rrr;
  const amrex::Real rer = qxp(iv, QREINT) + ekinr;

  const amrex::Real runewr = rur - flxu;
  const amrex::Real rvnewr = rvr - flxv;
  const amrex::Real rwnewr = rwr - flxw;
  const amrex::Real renewr = rer - flxe;

  qp(iv, QRHO) = rrnewr + hdt * srcrho;
  qp(iv, QU) = runewr / rrnewr + hdt * srcu;
  qp(iv, QV) = rvnewr / rrnewr + hdt * srcv;
  qp(iv, QW) = rwnewr / rrnewr + hdt * srcw;
  const amrex::Real rhoekinr =
    0.5 * (runewr * runewr + rvnewr * rvnewr + rwnewr * rwnewr) / rrnewr;
  qp(iv, QREINT) = renewr - rhoekinr + hdt * srce;
  amrex::Real temppres = qxp(iv, QPRES) - pznew - pynew + hdt * srcp;
  qp(iv, QPRES) =
    amrex::max<amrex::Real>(temppres, std::numeric_limits<amrex::Real>::min());

  // qm state
  const amrex::Real rul = rrl * qxm(ivp2, QU);
  const amrex::Real rvl = rrl * qxm(ivp2, QV);
  const amrex::Real rwl = rrl * qxm(ivp2, QW);
  const amrex::Real ekinl = 0.5 * (rul * rul + rvl * rvl + rwl * rwl) / rrl;
  const amrex::Real rel = qxm(ivp2, QREINT) + ekinl;

  const amrex::Real runewl = rul - flxu;
  const amrex::Real rvnewl = rvl - flxv;
  const amrex::Real rwnewl = rwl - flxw;
  const amrex::Real renewl = rel - flxe;

  qm(ivp2, QRHO) = rrnewl + hdt * srcrho;
  qm(ivp2, QU) = runewl / rrnewl + hdt * srcu;
  qm(ivp2, QV) = rvnewl / rrnewl + hdt * srcv;
  qm(ivp2, QW) = rwnewl / rrnewl + hdt * srcw;
  const amrex::Real rhoekinl =
    0.5 * (runewl * runewl + rvnewl * rvnewl + rwnewl * rwnewl) / rrnewl;
  qm(ivp2, QREINT) = renewl - rhoekinl + hdt * srce;
  temppres = qxm(ivp2, QPRES) - pznew - pynew + hdt * srcp;
  qm(ivp2, QPRES) =
    amrex::max<amrex::Real>(temppres, std::numeric_limits<amrex::Real>::min());
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
    div1 = difmag * amrex::min<amrex::Real>(0.0, div1);
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
    div1 = difmag * amrex::min<amrex::Real>(0.0, div1);
    for (int n = 0; n < NVAR; ++n) {
      if (n != UTEMP) {
        flx(i, j, k, n) +=
          dx * div1 * (u(i, j, k, n) - u(i, j - 1, k, n)); // Here dx is dy
      }
    }
  } else { // Z direction
    div1 = 0.25 * (div(i, j, k) + div(i + 1, j, k) + div(i, j + 1, k) +
                   div(i + 1, j + 1, k));
    div1 = difmag * amrex::min<amrex::Real>(0.0, div1);
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
#if AMREX_SPACEDIM == 3
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

#elif AMREX_SPACEDIM == 2

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

#endif
