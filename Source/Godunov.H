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

  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  const amrex::IntVect dvec(amrex::IntVect::TheDimensionVector(dir));
  const amrex::IntVect ivm3(iv - 3 * dvec);
  const amrex::IntVect ivm2(iv - 2 * dvec);
  const amrex::IntVect ivm(iv - dvec);
  const amrex::IntVect ivp(iv + dvec);
  const amrex::IntVect ivp2(iv + 2 * dvec);
  const amrex::IntVect ivp3(iv + 3 * dvec);

  if ((!q.contains(ivp3[0], ivp3[1], ivp3[2]) ||
       !q.contains(ivm3[0], ivm3[1], ivm3[2]))) {
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
  const amrex::Real qm = q(ivm, n);
  const amrex::Real qp = q(ivp, n);
  // assumes QU,QV,QW ordering
  const amrex::Real qmun = q(ivm, QU + dir);
  const amrex::Real qpun = q(ivp, QU + dir);
  const amrex::Real qmm = q(ivm2, n);
  const amrex::Real qpp = q(ivp2, n);

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

  const amrex::Real qpshft = q(iv + (1 - ishft) * dvec, n);
  const amrex::Real qmshft = q(iv - (1 + ishft) * dvec, n);
  const amrex::Real qpunshft = q(iv + (1 - ishft) * dvec, QU + dir);
  const amrex::Real qmunshft = q(iv - (1 + ishft) * dvec, QU + dir);
  const amrex::Real qppshft = q(iv + (2 - ishft) * dvec, n);
  const amrex::Real qmmshft = q(iv - (2 + ishft) * dvec, n);

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
  if (dir == 2) {
    idx = k;
  } else {
    idx = (dir == 0) ? i : j;
  }
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
  amrex::Real dummy_flx[NUM_SPECIES] = {0.0}; // FIXME later
  riemann(
    ql(i, j, k, QRHO), ul, vl, v2l, ql(i, j, k, QPRES), spl, qr(i, j, k, QRHO),
    ur, vr, v2r, qr(i, j, k, QPRES), spr, bc_test_val, cav, ustar,
    flx(i, j, k, URHO), dummy_flx, flx(i, j, k, f_idx[0]),
    flx(i, j, k, f_idx[1]), flx(i, j, k, f_idx[2]), flx(i, j, k, UEDEN),
    flx(i, j, k, UEINT), q(i, j, k, GU), q(i, j, k, GV), q(i, j, k, GV2),
    q(i, j, k, GDPRES), q(i, j, k, GDGAME));

  amrex::Real flxrho = flx(i, j, k, URHO);
  for (int ipass = 0; ipass < NPASSIVE; ++ipass) {
    int n = pmap.upassMap[ipass];
    int nqp = pmap.qpassMap[ipass];
    if (ustar > 0.0) {
      flx(i, j, k, n) = flxrho * ql(i, j, k, nqp);
    } else if (ustar < 0.0) {
      flx(i, j, k, n) = flxrho * qr(i, j, k, nqp);
    } else {
      flx(i, j, k, n) = flxrho * 0.50 * (ql(i, j, k, nqp) + qr(i, j, k, nqp));
    }
  }
}

// Transverse Correction for Predicted dir-states, using other_dir-Flux
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_transdo(
  const int i,
  const int j,
  const int k,
  const int dir,
  const int other_dir,
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
  const amrex::IntVect ivp0(iv + amrex::IntVect::TheDimensionVector(other_dir));
  const amrex::IntVect ivp1(iv + amrex::IntVect::TheDimensionVector(dir));
  const int qvidx = (other_dir == 0) ? GDU : ((other_dir == 1) ? GDV : GDW);

  const amrex::Real flxrho = cdtdx * (flxx(ivp0, URHO) - flxx(iv, URHO));
  const amrex::Real flxu = cdtdx * (flxx(ivp0, UMX) - flxx(iv, UMX));
  const amrex::Real flxv = cdtdx * (flxx(ivp0, UMY) - flxx(iv, UMY));
  const amrex::Real flxw = cdtdx * (flxx(ivp0, UMZ) - flxx(iv, UMZ));
  const amrex::Real flxe = cdtdx * (flxx(ivp0, UEDEN) - flxx(iv, UEDEN));
  const amrex::Real c = qa(iv, QGAMC);
  // Update passive variables
  for (int ipass = 0; ipass < NPASSIVE; ++ipass) {
    const int n = pmap.upassMap[ipass];
    const int nqp = pmap.qpassMap[ipass];
    const amrex::Real compn = cdtdx * (flxx(ivp0, n) - flxx(iv, n));
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
  const amrex::Real ugp = qint(ivp0, qvidx);
  const amrex::Real ugm = qint(iv, qvidx);

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
  const amrex::Real pnewl =
    qym(ivp1, QPRES) - cdtdx * (dAup + pav * dAu * (c - 1.0));
  qm(ivp1, QPRES) =
    amrex::max<amrex::Real>(pnewl, std::numeric_limits<amrex::Real>::min());
}

// dir corrected from other two dirs
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_transdd(
  const int i,
  const int j,
  const int k,
  const int dir,
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
  const int qvidx0 = (dir == 0) ? GDV : GDU;
  const int qvidx1 = (dir == 2) ? GDV : GDW;
  const amrex::GpuArray<const int, 3> bdim{
    {static_cast<int>(dir == 0), static_cast<int>(dir == 1),
     static_cast<int>(dir == 2)}};
  const amrex::GpuArray<const int, 3> l_idx{
    {bdim[0] * 0 + bdim[1] * 1 + bdim[2] * 2,
     bdim[0] * 1 + bdim[1] * 0 + bdim[2] * 0,
     bdim[0] * 2 + bdim[1] * 2 + bdim[2] * 1}};
  const amrex::IntVect ivp0(iv + amrex::IntVect::TheDimensionVector(l_idx[1]));
  const amrex::IntVect ivp1(iv + amrex::IntVect::TheDimensionVector(l_idx[2]));
  const amrex::IntVect ivp2(iv + amrex::IntVect::TheDimensionVector(dir));

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
  const amrex::Real ugpx = qx(ivp0, qvidx0);
  const amrex::Real ugmx = qx(iv, qvidx0);

  const amrex::Real dAupx = pggpx * ugpx - pggmx * ugmx;
  const amrex::Real pavx = 0.5 * (pggpx + pggmx);
  const amrex::Real dAux = ugpx - ugmx;

  const amrex::Real pggpy = qy(ivp1, GDPRES);
  const amrex::Real pggmy = qy(iv, GDPRES);
  const amrex::Real ugpy = qy(ivp1, qvidx1);
  const amrex::Real ugmy = qy(iv, qvidx1);

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

// 2D version of transdd and transdo
// Corrects the left and right dir-states with the other-Flux
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_transd(
  const int i,
  const int j,
  const int k,
  const int dir,
  amrex::Array4<amrex::Real> const& qm,
  amrex::Array4<amrex::Real> const& qp,
  amrex::Array4<const amrex::Real> const& qym,
  amrex::Array4<const amrex::Real> const& qyp,
  amrex::Array4<const amrex::Real> const& flxx,
  amrex::Array4<const amrex::Real> const& srcQ,
  amrex::Array4<const amrex::Real> const& qa,
  amrex::Array4<const amrex::Real> const& qint,
  const amrex::Real hdt,
  const amrex::Real cdtdx,
  PassMap const& pmap)
{
  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  const amrex::IntVect ivp0(
    iv + amrex::IntVect::TheDimensionVector(dir == 0 ? 1 : 0));
  const amrex::IntVect ivp1(
    iv + amrex::IntVect::TheDimensionVector(dir == 0 ? 0 : 1));
  const int qvidx = (dir == 0) ? GDV : GDU;

  const amrex::Real flxrho = cdtdx * (flxx(ivp0, URHO) - flxx(iv, URHO));
  const amrex::Real flxu = cdtdx * (flxx(ivp0, UMX) - flxx(iv, UMX));
  const amrex::Real flxv = cdtdx * (flxx(ivp0, UMY) - flxx(iv, UMY));
  const amrex::Real flxe = cdtdx * (flxx(ivp0, UEDEN) - flxx(iv, UEDEN));
  const amrex::Real srcr = srcQ(iv, QRHO);
  const amrex::Real srce = srcQ(iv, QREINT);
  const amrex::Real srcp = srcQ(iv, QPRES);
  const amrex::Real c = qa(iv, QGAMC);
  // Update passive variables
  for (int ipass = 0; ipass < NPASSIVE; ++ipass) {
    const int n = pmap.upassMap[ipass];
    const int nqp = pmap.qpassMap[ipass];
    const amrex::Real srcpass = srcQ(iv, nqp);
    const amrex::Real compn = cdtdx * (flxx(ivp0, n) - flxx(iv, n));
    amrex::Real rr = qyp(iv, QRHO);
    amrex::Real rrnew = rr - flxrho;
    amrex::Real compo = rr * qyp(iv, nqp) - compn;
    qp(iv, nqp) = compo / rrnew + hdt * srcpass;

    rr = qym(ivp1, QRHO);
    rrnew = rr - flxrho;
    compo = rr * qym(ivp1, nqp) - compn;
    qm(ivp1, nqp) = compo / rrnew + hdt * srcpass;
  }

  const amrex::Real pggp = qint(ivp0, GDPRES);
  const amrex::Real pggm = qint(iv, GDPRES);
  const amrex::Real ugp = qint(ivp0, qvidx);
  const amrex::Real ugm = qint(iv, qvidx);

  const amrex::Real dAup = pggp * ugp - pggm * ugm;
  const amrex::Real pav = 0.5 * (pggp + pggm);
  const amrex::Real dAu = ugp - ugm;

  // QP

  // Convert to conservative
  const amrex::Real rrr = qyp(iv, QRHO);
  amrex::Real rur = qyp(iv, QU);
  amrex::Real rvr = qyp(iv, QV);
  const amrex::Real ekinr = 0.5 * rrr * (rur * rur + rvr * rvr);
  rur *= rrr;
  rvr *= rrr;

  const amrex::Real rer = qyp(iv, QREINT) + ekinr;
  // Add transverse predictor
  const amrex::Real rrnewr = rrr - flxrho;
  const amrex::Real runewr = rur - flxu;
  const amrex::Real rvnewr = rvr - flxv;
  const amrex::Real renewr = rer - flxe;

  // Convert back to primitive
  qp(iv, QRHO) = rrnewr + hdt * srcr;
  qp(iv, QU) = runewr / rrnewr + hdt * srcQ(iv, QU);
  qp(iv, QV) = rvnewr / rrnewr + hdt * srcQ(iv, QV);
  const amrex::Real rhoekinr =
    0.5 * (runewr * runewr + rvnewr * rvnewr) / rrnewr;
  amrex::Real pnewr = qyp(iv, QPRES) - cdtdx * (dAup + pav * dAu * (c - 1.));
  pnewr += hdt * srcp;
  qp(iv, QPRES) =
    amrex::max<amrex::Real>(pnewr, std::numeric_limits<amrex::Real>::min());
  qp(iv, QREINT) = renewr - rhoekinr + hdt * srce;

  // QM

  // Conversion to Conservative
  const amrex::Real rrl = qym(ivp1, QRHO);
  amrex::Real rul = qym(ivp1, QU);
  amrex::Real rvl = qym(ivp1, QV);
  const amrex::Real ekinl = 0.5 * rrl * (rul * rul + rvl * rvl);
  rul *= rrl;
  rvl *= rrl;
  const amrex::Real rel = qym(ivp1, QREINT) + ekinl;

  // Transverse fluxes
  const amrex::Real rrnewl = rrl - flxrho;
  const amrex::Real runewl = rul - flxu;
  const amrex::Real rvnewl = rvl - flxv;
  const amrex::Real renewl = rel - flxe;

  qm(ivp1, QRHO) = rrnewl + hdt * srcr;
  qm(ivp1, QU) = runewl / rrnewl + hdt * srcQ(iv, QU);
  qm(ivp1, QV) = rvnewl / rrnewl + hdt * srcQ(iv, QV);
  const amrex::Real rhoekinl =
    0.5 * (runewl * runewl + rvnewl * rvnewl) / rrnewl;

  qm(ivp1, QREINT) = renewl - rhoekinl + hdt * srce;
  const amrex::Real pnewl =
    qym(ivp1, QPRES) - cdtdx * (dAup + pav * dAu * (c - 1.));
  qm(ivp1, QPRES) = amrex::max<amrex::Real>(
    pnewl + hdt * srcp, std::numeric_limits<amrex::Real>::min());
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
  amrex::Array4<const amrex::Real> const& divu,
  amrex::Array4<const amrex::Real> const& u,
  amrex::Real const dx,
  amrex::Real const difmag,
  const int dir)
{
  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  const amrex::IntVect ivm(iv - amrex::IntVect::TheDimensionVector(dir));
  const amrex::GpuArray<const int, 3> bdim{
    {static_cast<int>(dir == 0), static_cast<int>(dir == 1),
     static_cast<int>(dir == 2)}};
  const amrex::GpuArray<const int, 3> l_idx{
    {bdim[0] * 0 + bdim[1] * 1 + bdim[2] * 2,
     bdim[0] * 1 + bdim[1] * 0 + bdim[2] * 0,
     bdim[0] * 2 + bdim[1] * 2 + bdim[2] * 1}};

  const amrex::IntVect ivpj(iv + amrex::IntVect::TheDimensionVector(l_idx[1]));
  const amrex::IntVect ivpk(iv + amrex::IntVect::TheDimensionVector(l_idx[2]));
  const amrex::IntVect ivpp(
    iv + amrex::IntVect::TheDimensionVector(l_idx[1]) +
    amrex::IntVect::TheDimensionVector(l_idx[2]));
  const amrex::Real div =
    difmag *
    amrex::min<amrex::Real>(
      0.0, AMREX_D_PICK(
             divu(iv), 0.5 * (divu(iv) + divu(ivpj)),
             0.25 * (divu(iv) + divu(ivpj) + divu(ivpk) + divu(ivpp))));

  for (int n = 0; n < NVAR; ++n) {
    if (n != UTEMP) {
      flx(iv, n) += dx * div * (u(iv, n) - u(ivm, n));
    }
  }
  flx(iv, UTEMP) = 0.0;
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
  const bool use_flattening,
  const bool use_hybrid_weno,
  const int weno_scheme);

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
  const bool use_flattening,
  const bool use_hybrid_weno,
  const int weno_scheme);
#endif

#endif
