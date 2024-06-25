#ifndef GODUNOV_H
#define GODUNOV_H

#include <AMReX_FArrayBox.H>
#include <AMReX_IArrayBox.H>

#include "Constants.H"
#include "IndexDefines.H"
#include "Riemann.H"
#include "EB.H"

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_cmpflx_passive(
  const amrex::Real ustar,
  const amrex::Real flxrho,
  const amrex::Real& ql,
  const amrex::Real& qr,
  amrex::Real& flx)
{
  flx = (ustar > 0.0)   ? flxrho * ql
        : (ustar < 0.0) ? flxrho * qr
                        : flxrho * 0.5 * (ql + qr);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_transdo_passive(
  const amrex::IntVect iv,
  const amrex::IntVect ivpn,
  const amrex::IntVect ivpt,
  const int n,
  const int u_offset,
  const int q_offset,
  const amrex::Real cdtdx,
  const amrex::Real flxrho,
  const bool divide_by_rho,
  amrex::Array4<const amrex::Real> const& qnormp,
  amrex::Array4<const amrex::Real> const& qnormm,
  amrex::Array4<const amrex::Real> const& flxx,
  amrex::Array4<amrex::Real> const& qp,
  amrex::Array4<amrex::Real> const& qm,
  const bool no_cov_face,
  const bool lo_face_not_covered,
  const bool hi_face_not_covered)
{
  const int uc = u_offset + n;
  const int qc = q_offset + n;
  const amrex::Real compn =
    no_cov_face ? cdtdx * (flxx(ivpt, uc) - flxx(iv, uc)) : amrex::Real(0.0);

  amrex::Real rr, rrnew, compo;

  if (lo_face_not_covered) {
    rr = qnormp(iv, QRHO);
    rrnew = divide_by_rho ? rr - flxrho : 1.0;
    compo = rr * qnormp(iv, qc) - compn;
    qp(iv, qc) = compo / rrnew;
  } else {
    qp(iv, qc) = qnormp(iv, qc);
  }

  if (hi_face_not_covered) {
    rr = qnormm(ivpn, QRHO);
    rrnew = divide_by_rho ? rr - flxrho : 1.0;
    compo = rr * qnormm(ivpn, qc) - compn;
    qm(ivpn, qc) = compo / rrnew;
  } else {
    qm(ivpn, qc) = qnormm(ivpn, qc);
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_transdd_passive(
  const amrex::IntVect iv,
  const amrex::IntVect ivpn,
  const amrex::IntVect ivpt0,
  const amrex::IntVect ivpt1,
  const int n,
  const int u_offset,
  const int q_offset,
  const amrex::Real cdtdx0,
  const amrex::Real cdtdx1,
  const amrex::Real hdt,
  const amrex::Real rrr,
  const amrex::Real rrl,
  const amrex::Real rrnewr,
  const amrex::Real rrnewl,
  amrex::Array4<const amrex::Real> const& srcq,
  amrex::Array4<const amrex::Real> const& qnormp,
  amrex::Array4<const amrex::Real> const& qnormm,
  amrex::Array4<const amrex::Real> const& flxx,
  amrex::Array4<const amrex::Real> const& flxy,
  amrex::Array4<amrex::Real> const& qp,
  amrex::Array4<amrex::Real> const& qm,
  const bool no_cov_face,
  const bool lo_face_not_covered,
  const bool hi_face_not_covered)
{

  const int uc = u_offset + n;
  const int qc = q_offset + n;
  const amrex::Real srcpass = srcq(iv, qc);

  const amrex::Real compn = no_cov_face
                              ? cdtdx0 * (flxx(ivpt0, uc) - flxx(iv, uc)) +
                                  cdtdx1 * (flxy(ivpt1, uc) - flxy(iv, uc))
                              : amrex::Real(0.0);

  amrex::Real compo;

  // qp
  if (lo_face_not_covered) {
    compo = rrr * qnormp(iv, qc) - compn;
    qp(iv, qc) = compo / rrnewr + hdt * srcpass;
  } else {
    qp(iv, qc) = qnormp(iv, qc);
  }

  // qm

  if (hi_face_not_covered) {
    compo = rrl * qnormm(ivpn, qc) - compn;
    qm(ivpn, qc) = compo / rrnewl + hdt * srcpass;
  } else {
    qm(ivpn, qc) = qnormm(ivpn, qc);
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_transd_passive(
  const amrex::IntVect iv,
  const amrex::IntVect ivpn,
  const amrex::IntVect ivpt,
  const int n,
  const int u_offset,
  const int q_offset,
  const amrex::Real cdtdx,
  const amrex::Real hdt,
  const amrex::Real flxrho,
  const bool divide_by_rho,
  amrex::Array4<const amrex::Real> const& srcQ,
  amrex::Array4<const amrex::Real> const& qnormp,
  amrex::Array4<const amrex::Real> const& qnormm,
  amrex::Array4<const amrex::Real> const& flxx,
  amrex::Array4<amrex::Real> const& qp,
  amrex::Array4<amrex::Real> const& qm,
  const bool no_cov_face,
  const bool lo_face_not_covered,
  const bool hi_face_not_covered)
{
  const int uc = u_offset + n;
  const int qc = q_offset + n;
  const amrex::Real srcpass = srcQ(iv, qc);
  const amrex::Real compn =
    no_cov_face ? cdtdx * (flxx(ivpt, uc) - flxx(iv, uc)) : 0.0;

  amrex::Real rr, rrnew, compo;

  if (lo_face_not_covered) {
    rr = qnormp(iv, QRHO);
    rrnew = divide_by_rho ? rr - flxrho : 1.0;
    compo = rr * qnormp(iv, qc) - compn;
    qp(iv, qc) = compo / rrnew + hdt * srcpass;
  } else {
    qp(iv, qc) = qnormp(iv, qc);
  }

  if (hi_face_not_covered) {
    rr = qnormm(ivpn, QRHO);
    rrnew = divide_by_rho ? rr - flxrho : 1.0;
    compo = rr * qnormm(ivpn, qc) - compn;
    qm(ivpn, qc) = compo / rrnew + hdt * srcpass;
  } else {
    qm(ivpn, qc) = qnormm(ivpn, qc);
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
flatten(
  AMREX_D_DECL(const int i, const int j, const int k),
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

  if ((!q.contains(ivp3[0], ivp3[1], AMREX_D_PICK(0, 0, ivp3[2])) ||
       !q.contains(ivm3[0], ivm3[1], AMREX_D_PICK(0, 0, ivm3[2])))) {
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

  amrex::Real denom = amrex::max<amrex::Real>(small_pres, std::abs(qpp - qmm));
  amrex::Real zeta = std::abs(dp) / denom;
  const amrex::Real z = amrex::min<amrex::Real>(
    1.0, amrex::max<amrex::Real>(0.0, dzcut * (zeta - zcut1)));

  amrex::Real tst = (qmun - qpun >= 0.0) ? 1.0 : 0.0;

  amrex::Real tmp = amrex::min<amrex::Real>(qp, qm);

  const amrex::Real chi = (std::abs(dp) / tmp) > shktst ? tst : 0.0;

  const amrex::Real qpshft = q(iv + (1 - ishft) * dvec, n);
  const amrex::Real qmshft = q(iv - (1 + ishft) * dvec, n);
  const amrex::Real qpunshft = q(iv + (1 - ishft) * dvec, QU + dir);
  const amrex::Real qmunshft = q(iv - (1 + ishft) * dvec, QU + dir);
  const amrex::Real qppshft = q(iv + (2 - ishft) * dvec, n);
  const amrex::Real qmmshft = q(iv - (2 + ishft) * dvec, n);

  dp = qpshft - qmshft;

  denom = amrex::max<amrex::Real>(small_pres, std::abs(qppshft - qmmshft));
  zeta = std::abs(dp) / denom;
  const amrex::Real z2 = amrex::min<amrex::Real>(
    1.0, amrex::max<amrex::Real>(0.0, dzcut * (zeta - zcut1)));

  tst = (qmunshft - qpunshft >= 0.0) ? 1.0 : 0.0;

  tmp = amrex::min<amrex::Real>(qpshft, qmshft);

  const amrex::Real chi2 = (std::abs(dp) / tmp) > shktst ? tst : 0.0;

  return 1.0 - amrex::max<amrex::Real>(chi2 * z2, chi * z);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
flatten_eb(
  const int i,
  const int j,
  const int k,
  const int dir,
  amrex::Array4<amrex::EBCellFlag const> const& flags,
  amrex::Array4<const amrex::Real> const& q)
{

  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  const amrex::IntVect dvec(amrex::IntVect::TheDimensionVector(dir));

  const amrex::IntVect ivm3(iv - 3 * dvec);
  const amrex::IntVect ivp3(iv + 3 * dvec);

#if (AMREX_SPACEDIM == 2)
  if (
    (!q.contains(ivp3[0], ivp3[1], k)) || (!q.contains(ivm3[0], ivm3[1], k))) {
    return 1.0;
  }
#elif (AMREX_SPACEDIM == 3)
  if (
    !q.contains(ivp3[0], ivp3[1], ivp3[2]) ||
    !q.contains(ivm3[0], ivm3[1], ivm3[2])) {
    return 1.0;
  }
#endif

  const amrex::IntVect ng(3 * dvec);
  if (is_cut_neighborhood(iv, flags, ng)) {
    return 1.0;
  }
  return flatten(AMREX_D_DECL(i, j, k), dir, q);
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

  amrex::Real ul = ql(i, j, k, IU);
  amrex::Real vl = ql(i, j, k, IV);
  amrex::Real v2l = ql(i, j, k, IV2);
  amrex::Real pl = ql(i, j, k, QPRES);
  amrex::Real rhol = ql(i, j, k, QRHO);

  amrex::Real ur = qr(i, j, k, IU);
  amrex::Real vr = qr(i, j, k, IV);
  amrex::Real v2r = qr(i, j, k, IV2);
  amrex::Real pr = qr(i, j, k, QPRES);
  amrex::Real rhor = qr(i, j, k, QRHO);

  // Boundary condition corrections
  if (dir == 2) {
    idx = k;
  } else {
    idx = (dir == 0) ? i : j;
  }
  if (idx == domlo) {
    if (
      bclo == PCPhysBCType::no_slip_wall || bclo == PCPhysBCType::slip_wall ||
      bclo == PCPhysBCType::symmetry) {
      ul = -ur;
      vl = vr;
      v2l =
        v2r; // NoSlip: this is fine because Godunov velocity normal will be 0
      pl = pr;
      rhol = rhor;
    } else if (bclo == PCPhysBCType::outflow) {
      ul = ur;
      vl = vr;
      v2l = v2r;
      pl = pr;
      rhol = rhor;
    }
  } else if (idx == domhi + 1) {
    if (
      bchi == PCPhysBCType::no_slip_wall || bchi == PCPhysBCType::slip_wall ||
      bchi == PCPhysBCType::symmetry) {
      ur = -ul;
      vr = vl;
      v2r =
        v2l; // NoSlip: this is fine because Godunov velocity normal will be 0
      pr = pl;
      rhor = rhol;
    } else if (bchi == PCPhysBCType::outflow) {
      ur = ul;
      vr = vl;
      v2r = v2l;
      pr = pl;
      rhor = rhol;
    }
  }

  const int bc_test_val = 1;
  amrex::Real dummy_flx[NUM_SPECIES] = {0.0};

  // WIP: hard coding rotating frame stuff for now
  // soon to be fixed
  riemann(
    rhol, ul, vl, v2l, pl, spl, rhor, ur, vr, v2r, pr, spr, bc_test_val, cav,
    ustar, flx(i, j, k, URHO), dummy_flx, flx(i, j, k, f_idx[0]),
    flx(i, j, k, f_idx[1]), flx(i, j, k, f_idx[2]), flx(i, j, k, UEDEN),
    flx(i, j, k, UEINT), q(i, j, k, GU), q(i, j, k, GV), q(i, j, k, GV2),
    q(i, j, k, GDPRES), q(i, j, k, GDGAME));

  amrex::Real flxrho = flx(i, j, k, URHO);
  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
#if NUM_ADV > 0
  for (int n = 0; n < NUM_ADV; n++) {
    const int qc = QFA + n;
    pc_cmpflx_passive(ustar, flxrho, ql(iv, qc), qr(iv, qc), flx(iv, UFA + n));
  }
#endif
  for (int n = 0; n < NUM_SPECIES; n++) {
    const int qc = QFS + n;
    pc_cmpflx_passive(ustar, flxrho, ql(iv, qc), qr(iv, qc), flx(iv, UFS + n));
  }
#if NUM_AUX > 0
  for (int n = 0; n < NUM_AUX; n++) {
    const int qc = QFX + n;
    pc_cmpflx_passive(ustar, flxrho, ql(iv, qc), qr(iv, qc), flx(iv, UFX + n));
  }
#endif
#if NUM_LIN > 0
  for (int n = 0; n < NUM_LIN; n++) {
    const int qc = QLIN + n;
    pc_cmpflx_passive(
      ustar, q(i, j, k, GU), ql(iv, qc), qr(iv, qc), flx(iv, ULIN + n));
  }
#endif
}

// Transverse Correction for Predicted dir-states, using other_dir-Flux
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_transdo(
  AMREX_D_DECL(const int i, const int j, const int k),
  const int dir,
  const int other_dir,
  amrex::Array4<amrex::Real> const& qm,
  amrex::Array4<amrex::Real> const& qp,
  amrex::Array4<const amrex::Real> const& qnormm,
  amrex::Array4<const amrex::Real> const& qnormp,
  amrex::Array4<const amrex::Real> const& flxx,
  amrex::Array4<const amrex::Real> const& qa,
  amrex::Array4<const amrex::Real> const& qint,
  const amrex::Real cdtdx,
  amrex::Array4<const amrex::Real> const& norm_area = {},
  amrex::Array4<const amrex::Real> const& trans_area = {})
{
  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  const amrex::IntVect ivpn(iv + amrex::IntVect::TheDimensionVector(dir));
  const amrex::IntVect ivpt(iv + amrex::IntVect::TheDimensionVector(other_dir));
  const int qvidx = (other_dir == 0) ? GDU : ((other_dir == 1) ? GDV : GDW);

  // clang-format off
  // ivpn is an offset by 1 in the     normal direction
  // ivpt is an offset by 1 in the transverse direction
  // If (dir == 0 and other_dir == 1),  ivpt = (i,j+1,k) and ivpn = (i+1,j,k) and qvidx = GDV
  // If (dir == 0 and other_dir == 2),  ivpt = (i,j,k+1) and ivpn = (i+1,j,k) and qvidx = GDW
  // If (dir == 1 and other_dir == 0),  ivpt = (i+1,j,k) and ivpn = (i,j+1,k) and qvidx = GDU
  // If (dir == 1 and other_dir == 2),  ivpt = (i,j,k+1) and ivpn = (i,j+1,k) and qvidx = GDW
  // If (dir == 2 and other_dir == 0),  ivpt = (i+1,j,k) and ivpn = (i,j,k+1) and qvidx = GDU
  // If (dir == 2 and other_dir == 1),  ivpt = (i,j+1,k) and ivpn = (i,j,k+1) and qvidx = GDV
  // clang-format on

  bool no_cov_face = true;
  bool lo_face_not_covered = true;
  bool hi_face_not_covered = true;

  if (trans_area) {
    no_cov_face = (trans_area(ivpt) > 0.0) && (trans_area(iv) > 0.0);
  }
  if (norm_area) {
    lo_face_not_covered = (norm_area(iv) > amrex::Real(0.0));
    hi_face_not_covered = (norm_area(ivpn) > amrex::Real(0.0));
  }

  const amrex::Real flxrho = no_cov_face
                               ? cdtdx * (flxx(ivpt, URHO) - flxx(iv, URHO))
                               : amrex::Real(0.0);
  const amrex::Real flxu =
    no_cov_face ? cdtdx * (flxx(ivpt, UMX) - flxx(iv, UMX)) : amrex::Real(0.0);
  const amrex::Real flxv =
    no_cov_face ? cdtdx * (flxx(ivpt, UMY) - flxx(iv, UMY)) : amrex::Real(0.0);
  const amrex::Real flxw =
    no_cov_face ? cdtdx * (flxx(ivpt, UMZ) - flxx(iv, UMZ)) : amrex::Real(0.0);
  const amrex::Real flxe = no_cov_face
                             ? cdtdx * (flxx(ivpt, UEDEN) - flxx(iv, UEDEN))
                             : amrex::Real(0.0);
  const amrex::Real c = qa(iv, QGAMC);

  // Update passive variables
#if NUM_ADV > 0
  for (int n = 0; n < NUM_ADV; n++) {
    pc_transdo_passive(
      iv, ivpn, ivpt, n, UFA, QFA, cdtdx, flxrho, true, qnormp, qnormm, flxx,
      qp, qm, no_cov_face, lo_face_not_covered, hi_face_not_covered);
  }
#endif
  for (int n = 0; n < NUM_SPECIES; n++) {
    pc_transdo_passive(
      iv, ivpn, ivpt, n, UFS, QFS, cdtdx, flxrho, true, qnormp, qnormm, flxx,
      qp, qm, no_cov_face, lo_face_not_covered, hi_face_not_covered);
  }
#if NUM_AUX > 0
  for (int n = 0; n < NUM_AUX; n++) {
    pc_transdo_passive(
      iv, ivpn, ivpt, n, UFX, QFX, cdtdx, flxrho, true, qnormp, qnormm, flxx,
      qp, qm, no_cov_face, lo_face_not_covered, hi_face_not_covered);
  }
#endif
#if NUM_LIN > 0
  for (int n = 0; n < NUM_LIN; n++) {
    pc_transdo_passive(
      iv, ivpn, ivpt, n, ULIN, QLIN, cdtdx, 0., false, qnormp, qnormm, flxx, qp,
      qm, no_cov_face, lo_face_not_covered, hi_face_not_covered);
  }
#endif

  // Update hydro vars
  amrex::Real dAup = 0.0, pav = 0.0, dAu = 0.0;
  if (no_cov_face) {
    const amrex::Real pggp = qint(ivpt, GDPRES);
    const amrex::Real pggm = qint(iv, GDPRES);
    const amrex::Real ugp = qint(ivpt, qvidx);
    const amrex::Real ugm = qint(iv, qvidx);

    dAup = pggp * ugp - pggm * ugm;
    pav = 0.5 * (pggp + pggm);
    dAu = ugp - ugm;
  }

  // QP

  if (lo_face_not_covered) {
    // Convert to conservative
    const amrex::Real rrr = qnormp(iv, QRHO);
    amrex::Real rur = qnormp(iv, QU);
    amrex::Real rvr = qnormp(iv, QV);
    amrex::Real rwr = qnormp(iv, QW);
    const amrex::Real ekinr = 0.5 * rrr * (rur * rur + rvr * rvr + rwr * rwr);
    rur *= rrr;
    rvr *= rrr;
    rwr *= rrr;

    const amrex::Real rer = qnormp(iv, QREINT) + ekinr;
    // Add transverse predictor
    amrex::Real rrnewr = rrr - flxrho;
    amrex::Real runewr = rur - flxu;
    amrex::Real rvnewr = rvr - flxv;
    amrex::Real rwnewr = rwr - flxw;
    amrex::Real renewr = rer - flxe;

    // Convert back to primitive
    qp(iv, QRHO) = rrnewr;
    qp(iv, QU) = runewr / rrnewr;
    qp(iv, QV) = rvnewr / rrnewr;
    qp(iv, QW) = rwnewr / rrnewr;

    const amrex::Real rhoekinr =
      0.5 * (runewr * runewr + rvnewr * rvnewr + rwnewr * rwnewr) / rrnewr;

    qp(iv, QREINT) = renewr - rhoekinr;

    amrex::Real pnewr =
      qnormp(iv, QPRES) - cdtdx * (dAup + pav * dAu * (c - 1.0));
    ;
    qp(iv, QPRES) = amrex::max<amrex::Real>(pnewr, constants::very_small_num());
  } else {
    qp(iv, QRHO) = qnormp(iv, QRHO);
    qp(iv, QU) = qnormp(iv, QU);
    qp(iv, QV) = qnormp(iv, QV);
    qp(iv, QREINT) = qnormp(iv, QREINT);
    qp(iv, QPRES) = qnormp(iv, QPRES);
  }

  // QM

  if (hi_face_not_covered) {
    // Conversion to Conservative
    const amrex::Real rrl = qnormm(ivpn, QRHO);
    amrex::Real rul = qnormm(ivpn, QU);
    amrex::Real rvl = qnormm(ivpn, QV);
    amrex::Real rwl = qnormm(ivpn, QW);
    const amrex::Real ekinl = 0.5 * rrl * (rul * rul + rvl * rvl + rwl * rwl);
    rul *= rrl;
    rvl *= rrl;
    rwl *= rrl;
    const amrex::Real rel = qnormm(ivpn, QREINT) + ekinl;

    // Transverse fluxes
    amrex::Real rrnewl = rrl - flxrho;
    amrex::Real runewl = rul - flxu;
    amrex::Real rvnewl = rvl - flxv;
    amrex::Real rwnewl = rwl - flxw;
    amrex::Real renewl = rel - flxe;

    qm(ivpn, QRHO) = rrnewl;
    qm(ivpn, QU) = runewl / rrnewl;
    qm(ivpn, QV) = rvnewl / rrnewl;
    qm(ivpn, QW) = rwnewl / rrnewl;
    const amrex::Real rhoekinl =
      0.5 * (runewl * runewl + rvnewl * rvnewl + rwnewl * rwnewl) / rrnewl;

    qm(ivpn, QREINT) = renewl - rhoekinl;

    amrex::Real pnewl =
      qnormm(ivpn, QPRES) - cdtdx * (dAup + pav * dAu * (c - 1.0));
    qm(ivpn, QPRES) =
      amrex::max<amrex::Real>(pnewl, constants::very_small_num());

  } else {
    qm(ivpn, QRHO) = qnormm(ivpn, QRHO);
    qm(ivpn, QU) = qnormm(ivpn, QU);
    qm(ivpn, QV) = qnormm(ivpn, QV);
    qm(ivpn, QW) = qnormm(ivpn, QW);
    qm(ivpn, QREINT) = qnormm(ivpn, QREINT);
    qm(ivpn, QPRES) = qnormm(ivpn, QPRES);
  }
}

// dir corrected from other two dirs
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_transdd(
  AMREX_D_DECL(const int i, const int j, const int k),
  const int dir,
  amrex::Array4<amrex::Real> const& qm,
  amrex::Array4<amrex::Real> const& qp,
  amrex::Array4<const amrex::Real> const& qnormm,
  amrex::Array4<const amrex::Real> const& qnormp,
  amrex::Array4<const amrex::Real> const& flxx,
  amrex::Array4<const amrex::Real> const& flxy,
  amrex::Array4<const amrex::Real> const& qx,
  amrex::Array4<const amrex::Real> const& qy,
  amrex::Array4<const amrex::Real> const& qa,
  amrex::Array4<const amrex::Real> const& srcq,
  const amrex::Real hdt,
  const amrex::Real cdtdx0,
  const amrex::Real cdtdx1,
  amrex::Array4<const amrex::Real> const& norm_area = {},
  amrex::Array4<const amrex::Real> const& trans0_area = {},
  amrex::Array4<const amrex::Real> const& trans1_area = {})
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
  const amrex::IntVect ivpt0(iv + amrex::IntVect::TheDimensionVector(l_idx[1]));
  const amrex::IntVect ivpt1(iv + amrex::IntVect::TheDimensionVector(l_idx[2]));
  const amrex::IntVect ivpn(iv + amrex::IntVect::TheDimensionVector(dir));

  // clang-format off
  // ivpn is an offset by 1 in the normal direction
  // ivpt0 is an offset by 1 in the first transverse direction
  // ivpt1 is an offset by 1 in the second transverse direction
  // If (dir == 0),  ivpt0 = (i,j+1,k) and ivpt1 = (i,j,k+1) and ivpn = (i+1,j,k) and qvidx0 = GDV and qvidx1 = GVW
  // If (dir == 1),  ivpt0 = (i+1,j,k) and ivpt1 = (i,j,k+1) and ivpn = (i,j+1,k) and qvidx0 = GDU and qvidx1 = GVW
  // If (dir == 2),  ivpt0 = (i+1,j,k) and ivpt1 = (i,j+1,k) and ivpn = (i,j,k+1) and qvidx0 = GDU and qvidx1 = GVV
  // clang-format on

  bool no_cov_face = true;
  bool lo_face_not_covered = true;
  bool hi_face_not_covered = true;
#ifdef AMREX_USE_EB
  if (trans0_area && trans1_area) {
    no_cov_face = (trans0_area(ivpt0) > 0.0) && (trans1_area(ivpt1) > 0.0) &&
                  (trans0_area(iv) > 0.0) && (trans1_area(iv) > 0.0);
  }
  if (norm_area) {
    lo_face_not_covered = (norm_area(iv) > amrex::Real(0.0));
    hi_face_not_covered = (norm_area(ivpn) > amrex::Real(0.0));
  }
#endif

  const amrex::Real flxrho = no_cov_face
                               ? cdtdx0 * (flxx(ivpt0, URHO) - flxx(iv, URHO)) +
                                   cdtdx1 * (flxy(ivpt1, URHO) - flxy(iv, URHO))
                               : amrex::Real(0.0);
  const amrex::Real flxu = no_cov_face
                             ? cdtdx0 * (flxx(ivpt0, UMX) - flxx(iv, UMX)) +
                                 cdtdx1 * (flxy(ivpt1, UMX) - flxy(iv, UMX))
                             : amrex::Real(0.0);
  const amrex::Real flxv = no_cov_face
                             ? cdtdx0 * (flxx(ivpt0, UMY) - flxx(iv, UMY)) +
                                 cdtdx1 * (flxy(ivpt1, UMY) - flxy(iv, UMY))
                             : amrex::Real(0.0);
  const amrex::Real flxw = no_cov_face
                             ? cdtdx0 * (flxx(ivpt0, UMZ) - flxx(iv, UMZ)) +
                                 cdtdx1 * (flxy(ivpt1, UMZ) - flxy(iv, UMZ))
                             : amrex::Real(0.0);
  const amrex::Real flxe = no_cov_face
                             ? cdtdx0 * (flxx(ivpt0, UEDEN) - flxx(iv, UEDEN)) +
                                 cdtdx1 * (flxy(ivpt1, UEDEN) - flxy(iv, UEDEN))
                             : amrex::Real(0.0);
  const amrex::Real c = qa(iv, QGAMC);
  const amrex::Real srcrho = srcq(iv, QRHO);
  const amrex::Real srcu = srcq(iv, QU);
  const amrex::Real srcv = srcq(iv, QV);
  const amrex::Real srcw = srcq(iv, QW);
  const amrex::Real srce = srcq(iv, QREINT);
  const amrex::Real srcp = srcq(iv, QPRES);

  // Update passive variables
  const amrex::Real rrr = qnormp(iv, QRHO);
  const amrex::Real rrl = qnormm(ivpn, QRHO);
  const amrex::Real rrnewl = rrl - flxrho;
  const amrex::Real rrnewr = rrr - flxrho;
#if NUM_ADV > 0
  for (int n = 0; n < NUM_ADV; n++) {
    pc_transdd_passive(
      iv, ivpn, ivpt0, ivpt1, n, UFA, QFA, cdtdx0, cdtdx1, hdt, rrr, rrl,
      rrnewr, rrnewl, srcq, qnormp, qnormm, flxx, flxy, qp, qm, no_cov_face,
      lo_face_not_covered, hi_face_not_covered);
  }
#endif
  for (int n = 0; n < NUM_SPECIES; n++) {
    pc_transdd_passive(
      iv, ivpn, ivpt0, ivpt1, n, UFS, QFS, cdtdx0, cdtdx1, hdt, rrr, rrl,
      rrnewr, rrnewl, srcq, qnormp, qnormm, flxx, flxy, qp, qm, no_cov_face,
      lo_face_not_covered, hi_face_not_covered);
  }
#if NUM_AUX > 0
  for (int n = 0; n < NUM_AUX; n++) {
    pc_transdd_passive(
      iv, ivpn, ivpt0, ivpt1, n, UFX, QFX, cdtdx0, cdtdx1, hdt, rrr, rrl,
      rrnewr, rrnewl, srcq, qnormp, qnormm, flxx, flxy, qp, qm, no_cov_face,
      lo_face_not_covered, hi_face_not_covered);
  }
#endif
#if NUM_LIN > 0
  for (int n = 0; n < NUM_LIN; n++) {
    pc_transdd_passive(
      iv, ivpn, ivpt0, ivpt1, n, ULIN, QLIN, cdtdx0, cdtdx1, hdt, 1., 1., 1.,
      1., srcq, qnormp, qnormm, flxx, flxy, qp, qm, no_cov_face,
      lo_face_not_covered, hi_face_not_covered);
  }
#endif

  // Update hydro vars
  // Update hydro vars
  amrex::Real dAupx = 0.0, pavx = 0.0, dAux = 0.0;
  amrex::Real dAupy = 0.0, pavy = 0.0, dAuy = 0.0;
  if (no_cov_face) {
    const amrex::Real pggpx = qx(ivpt0, GDPRES);
    const amrex::Real pggmx = qx(iv, GDPRES);
    const amrex::Real ugpx = qx(ivpt0, qvidx0);
    const amrex::Real ugmx = qx(iv, qvidx0);

    dAupx = pggpx * ugpx - pggmx * ugmx;
    pavx = 0.5 * (pggpx + pggmx);
    dAux = ugpx - ugmx;

    const amrex::Real pggpy = qy(ivpt1, GDPRES);
    const amrex::Real pggmy = qy(iv, GDPRES);
    const amrex::Real ugpy = qy(ivpt1, qvidx1);
    const amrex::Real ugmy = qy(iv, qvidx1);

    dAupy = pggpy * ugpy - pggmy * ugmy;
    pavy = 0.5 * (pggpy + pggmy);
    dAuy = ugpy - ugmy;
  }
  const amrex::Real pxnew = cdtdx0 * (dAupx + pavx * dAux * (c - 1.0));
  const amrex::Real pynew = cdtdx1 * (dAupy + pavy * dAuy * (c - 1.0));

  // qp state
  if (lo_face_not_covered) {
    const amrex::Real rur = rrr * qnormp(iv, QU);
    const amrex::Real rvr = rrr * qnormp(iv, QV);
    const amrex::Real rwr = rrr * qnormp(iv, QW);
    const amrex::Real ekinr = 0.5 * (rur * rur + rvr * rvr + rwr * rwr) / rrr;
    const amrex::Real rer = qnormp(iv, QREINT) + ekinr;

    amrex::Real runewr = rur - flxu;
    amrex::Real rvnewr = rvr - flxv;
    amrex::Real rwnewr = rwr - flxw;
    amrex::Real renewr = rer - flxe;

    qp(iv, QRHO) = rrnewr + hdt * srcrho;
    qp(iv, QU) = runewr / rrnewr + hdt * srcu;
    qp(iv, QV) = rvnewr / rrnewr + hdt * srcv;
    qp(iv, QW) = rwnewr / rrnewr + hdt * srcw;

    const amrex::Real rhoekinr =
      0.5 * (runewr * runewr + rvnewr * rvnewr + rwnewr * rwnewr) / rrnewr;
    qp(iv, QREINT) = renewr - rhoekinr + hdt * srce;

    qp(iv, QPRES) = qnormp(iv, QPRES) - pxnew - pynew + hdt * srcp;
    qp(iv, QPRES) = amrex::max<amrex::Real>(
      qp(iv, QPRES), std::numeric_limits<amrex::Real>::min());
  } else {
    qp(iv, QRHO) = qnormp(iv, QRHO);
    qp(iv, QU) = qnormp(iv, QU);
    qp(iv, QV) = qnormp(iv, QV);
    qp(iv, QREINT) = qnormp(iv, QREINT);
    qp(iv, QPRES) = qnormp(iv, QPRES);
  }

  // qm state
  if (hi_face_not_covered) {
    const amrex::Real rul = rrl * qnormm(ivpn, QU);
    const amrex::Real rvl = rrl * qnormm(ivpn, QV);
    const amrex::Real rwl = rrl * qnormm(ivpn, QW);
    const amrex::Real ekinl = 0.5 * (rul * rul + rvl * rvl + rwl * rwl) / rrl;
    const amrex::Real rel = qnormm(ivpn, QREINT) + ekinl;

    amrex::Real runewl = rul - flxu;
    amrex::Real rvnewl = rvl - flxv;
    amrex::Real rwnewl = rwl - flxw;
    amrex::Real renewl = rel - flxe;

    qm(ivpn, QRHO) = rrnewl + hdt * srcrho;
    qm(ivpn, QU) = runewl / rrnewl + hdt * srcu;
    qm(ivpn, QV) = rvnewl / rrnewl + hdt * srcv;
    qm(ivpn, QW) = rwnewl / rrnewl + hdt * srcw;

    const amrex::Real rhoekinl =
      0.5 * (runewl * runewl + rvnewl * rvnewl + rwnewl * rwnewl) / rrnewl;
    qm(ivpn, QREINT) = renewl - rhoekinl + hdt * srce;

    qm(ivpn, QPRES) = qnormm(ivpn, QPRES) - pxnew - pynew + hdt * srcp;
    qm(ivpn, QPRES) = amrex::max<amrex::Real>(
      qm(ivpn, QPRES), std::numeric_limits<amrex::Real>::min());

  } else {
    qm(ivpn, QRHO) = qnormm(ivpn, QRHO);
    qm(ivpn, QU) = qnormm(ivpn, QU);
    qm(ivpn, QV) = qnormm(ivpn, QV);
    qm(ivpn, QW) = qnormm(ivpn, QW);
    qm(ivpn, QREINT) = qnormm(ivpn, QREINT);
    qm(ivpn, QPRES) = qnormm(ivpn, QPRES);
  }
}

// 2D version of transdd and transdo
// Corrects the left and right dir-states with the other-Flux
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_transd(
  AMREX_D_DECL(const int i, const int j, const int k),
  const int dir,
  amrex::Array4<amrex::Real> const& qm,
  amrex::Array4<amrex::Real> const& qp,
  amrex::Array4<const amrex::Real> const& qnormm,
  amrex::Array4<const amrex::Real> const& qnormp,
  amrex::Array4<const amrex::Real> const& flxx,
  amrex::Array4<const amrex::Real> const& srcQ,
  amrex::Array4<const amrex::Real> const& qa,
  amrex::Array4<const amrex::Real> const& qint,
  const amrex::Real hdt,
  const amrex::Real cdtdx,
  amrex::Array4<const amrex::Real> const& norm_area = {},
  amrex::Array4<const amrex::Real> const& trans_area = {})
{
  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  // ivpt is an offset by 1 in the transverse direction
  // ivpn is an offset by 1 in the     normal direction
  // If (dir == 0),  ivpt = (i,j+1) and ivpn = (i+1,j) and qvidx = GDV
  // If (dir == 1),  ivpt = (i+1,j) and ivpn = (i,j+1) and qvidx = GDU
  const amrex::IntVect ivpt(
    iv + amrex::IntVect::TheDimensionVector(dir == 0 ? 1 : 0));
  const amrex::IntVect ivpn(
    iv + amrex::IntVect::TheDimensionVector(dir == 0 ? 0 : 1));
  const int qvidx = (dir == 0) ? GDV : GDU;

  bool no_cov_face = true;
  bool lo_face_not_covered = true;
  bool hi_face_not_covered = true;
#ifdef AMREX_USE_EB
  if (trans_area) {
    no_cov_face = (trans_area(ivpt) > 0.0) && (trans_area(iv) > 0.0);
  }
  if (norm_area) {
    lo_face_not_covered = (norm_area(iv) > amrex::Real(0.0));
    hi_face_not_covered = (norm_area(ivpn) > amrex::Real(0.0));
  }
#endif

  const amrex::Real flxrho = no_cov_face
                               ? cdtdx * (flxx(ivpt, URHO) - flxx(iv, URHO))
                               : amrex::Real(0.0);
  const amrex::Real flxu =
    no_cov_face ? cdtdx * (flxx(ivpt, UMX) - flxx(iv, UMX)) : amrex::Real(0.0);
  const amrex::Real flxv =
    no_cov_face ? cdtdx * (flxx(ivpt, UMY) - flxx(iv, UMY)) : amrex::Real(0.0);
  const amrex::Real flxe = no_cov_face
                             ? cdtdx * (flxx(ivpt, UEDEN) - flxx(iv, UEDEN))
                             : amrex::Real(0.0);
  const amrex::Real srcr = srcQ(iv, QRHO);
  const amrex::Real srce = srcQ(iv, QREINT);
  const amrex::Real srcp = srcQ(iv, QPRES);
  const amrex::Real c = qa(iv, QGAMC);

  // Update passive variables
#if AMREX_SPACEDIM == 2
  pc_transd_passive(
    iv, ivpn, ivpt, 0, UMZ, QW, cdtdx, hdt, flxrho, true, srcQ, qnormp, qnormm,
    flxx, qp, qm, no_cov_face, lo_face_not_covered, hi_face_not_covered);
#endif
  for (int n = 0; n < NUM_ADV; n++) {
    pc_transd_passive(
      iv, ivpn, ivpt, n, UFA, QFA, cdtdx, hdt, flxrho, true, srcQ, qnormp,
      qnormm, flxx, qp, qm, no_cov_face, lo_face_not_covered,
      hi_face_not_covered);
  }
  for (int n = 0; n < NUM_SPECIES; n++) {
    pc_transd_passive(
      iv, ivpn, ivpt, n, UFS, QFS, cdtdx, hdt, flxrho, true, srcQ, qnormp,
      qnormm, flxx, qp, qm, no_cov_face, lo_face_not_covered,
      hi_face_not_covered);
  }
  for (int n = 0; n < NUM_AUX; n++) {
    pc_transd_passive(
      iv, ivpn, ivpt, n, UFX, QFX, cdtdx, hdt, flxrho, true, srcQ, qnormp,
      qnormm, flxx, qp, qm, no_cov_face, lo_face_not_covered,
      hi_face_not_covered);
  }
  for (int n = 0; n < NUM_LIN; n++) {
    pc_transd_passive(
      iv, ivpn, ivpt, n, ULIN, QLIN, cdtdx, hdt, 0., false, srcQ, qnormp,
      qnormm, flxx, qp, qm, no_cov_face, lo_face_not_covered,
      hi_face_not_covered);
  }

  amrex::Real dAup = 0.0, pav = 0.0, dAu = 0.0;
  if (no_cov_face) {
    const amrex::Real pggp = qint(ivpt, GDPRES);
    const amrex::Real pggm = qint(iv, GDPRES);
    const amrex::Real ugp = qint(ivpt, qvidx);
    const amrex::Real ugm = qint(iv, qvidx);

    dAup = pggp * ugp - pggm * ugm;
    pav = 0.5 * (pggp + pggm);
    dAu = ugp - ugm;
  }

  // QP

  if (lo_face_not_covered) {
    // Convert to conservative
    const amrex::Real rrr = qnormp(iv, QRHO);
    amrex::Real rur = qnormp(iv, QU);
    amrex::Real rvr = qnormp(iv, QV);
    const amrex::Real ekinr = 0.5 * rrr * (rur * rur + rvr * rvr);
    rur *= rrr;
    rvr *= rrr;

    const amrex::Real rer = qnormp(iv, QREINT) + ekinr;
    // Add transverse predictor
    amrex::Real rrnewr = rrr - flxrho;
    amrex::Real runewr = rur - flxu;
    amrex::Real rvnewr = rvr - flxv;
    amrex::Real renewr = rer - flxe;

    // Convert back to primitive
    qp(iv, QRHO) = rrnewr + hdt * srcr;
    qp(iv, QU) = runewr / rrnewr + hdt * srcQ(iv, QU);
    qp(iv, QV) = rvnewr / rrnewr + hdt * srcQ(iv, QV);

    const amrex::Real rhoekinr =
      0.5 * (runewr * runewr + rvnewr * rvnewr) / rrnewr;
    qp(iv, QREINT) = renewr - rhoekinr + hdt * srce;

    amrex::Real pnewr =
      qnormp(iv, QPRES) - cdtdx * (dAup + pav * dAu * (c - 1.));
    pnewr += hdt * srcp;
    qp(iv, QPRES) = pnewr;

    qp(iv, QPRES) = amrex::max<amrex::Real>(
      qp(iv, QPRES), std::numeric_limits<amrex::Real>::min());

  } else {
    qp(iv, QRHO) = qnormp(iv, QRHO);
    qp(iv, QU) = qnormp(iv, QU);
    qp(iv, QV) = qnormp(iv, QV);
    qp(iv, QREINT) = qnormp(iv, QREINT);
    qp(iv, QPRES) = qnormp(iv, QPRES);
  }

  // QM

  if (hi_face_not_covered) {
    // Conversion to Conservative
    const amrex::Real rrl = qnormm(ivpn, QRHO);
    amrex::Real rul = qnormm(ivpn, QU);
    amrex::Real rvl = qnormm(ivpn, QV);
    const amrex::Real ekinl = 0.5 * rrl * (rul * rul + rvl * rvl);
    rul *= rrl;
    rvl *= rrl;
    const amrex::Real rel = qnormm(ivpn, QREINT) + ekinl;

    // Transverse fluxes
    amrex::Real rrnewl = rrl - flxrho;
    amrex::Real runewl = rul - flxu;
    amrex::Real rvnewl = rvl - flxv;
    amrex::Real renewl = rel - flxe;

    qm(ivpn, QRHO) = rrnewl + hdt * srcr;
    qm(ivpn, QU) = runewl / rrnewl + hdt * srcQ(iv, QU);
    qm(ivpn, QV) = rvnewl / rrnewl + hdt * srcQ(iv, QV);

    const amrex::Real rhoekinl =
      0.5 * (runewl * runewl + rvnewl * rvnewl) / rrnewl;
    qm(ivpn, QREINT) = renewl - rhoekinl + hdt * srce;

    amrex::Real pnewl =
      qnormm(ivpn, QPRES) - cdtdx * (dAup + pav * dAu * (c - 1.));
    qm(ivpn, QPRES) = pnewl + hdt * srcp;

    qm(ivpn, QPRES) = amrex::max<amrex::Real>(
      qm(ivpn, QPRES), std::numeric_limits<amrex::Real>::min());

  } else {
    qm(ivpn, QRHO) = qnormm(ivpn, QRHO);
    qm(ivpn, QU) = qnormm(ivpn, QU);
    qm(ivpn, QV) = qnormm(ivpn, QV);
    qm(ivpn, QREINT) = qnormm(ivpn, QREINT);
    qm(ivpn, QPRES) = qnormm(ivpn, QPRES);
  }
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
  AMREX_D_DECL(const int i, const int j, const int k),
  amrex::Array4<amrex::Real> const& flx,
  amrex::Array4<const amrex::Real> const& divu,
  amrex::Array4<const amrex::Real> const& u,
  amrex::Real const dx,
  amrex::Real const difmag,
  const int dir,
  const int domlo,
  const int domhi,
  const int bclo,
  const int bchi)
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
  AMREX_D_TERM(, ,
               const amrex::IntVect ivpk(
                 iv + amrex::IntVect::TheDimensionVector(l_idx[2]));
               const amrex::IntVect ivpp(
                 iv + amrex::IntVect::TheDimensionVector(l_idx[1]) +
                 amrex::IntVect::TheDimensionVector(l_idx[2])););
  const amrex::Real div =
    difmag *
    amrex::min<amrex::Real>(
      0.0, AMREX_D_PICK(
             divu(iv), 0.5 * (divu(iv) + divu(ivpj)),
             0.25 * (divu(iv) + divu(ivpj) + divu(ivpk) + divu(ivpp))));

  bool at_bndry =
    ((dir == 0) && (i == domlo) && (bclo == amrex::PhysBCType::noslipwall)) ||
    ((dir == 0) && (i == domhi + 1) &&
     (bchi == amrex::PhysBCType::noslipwall)) ||
    ((dir == 1) && (j == domlo) && (bclo == amrex::PhysBCType::noslipwall)) ||
    ((dir == 1) && (j == domhi + 1) && (bchi == amrex::PhysBCType::noslipwall));
#if (AMREX_SPACEDIM == 3)
  at_bndry =
    at_bndry ||
    ((dir == 2) && (k == domlo) && (bclo == amrex::PhysBCType::noslipwall)) ||
    ((dir == 2) && (k == domhi + 1) && (bchi == amrex::PhysBCType::noslipwall));
#endif

  for (int n = 0; n < NVAR; ++n) {
    if (n != UTEMP) {
#if (AMREX_SPACEDIM == 2)
      bool is_tang_vel = (dir == 0 && n == UMY) || (dir == 1 && n == UMX);
#elif (AMREX_SPACEDIM == 3)
      bool is_tang_vel = (dir == 0 && (n == UMY || n == UMZ)) ||
                         (dir == 1 && (n == UMX || n == UMZ)) ||
                         (dir == 2 && (n == UMX || n == UMY));
#endif
      if (!(at_bndry && is_tang_vel)) {
        flx(iv, n) += dx * div * (u(iv, n) - u(ivm, n));
      }
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
  const amrex::GpuArray<const amrex::Array4<amrex::Real>, 3>& flx,
  const amrex::GpuArray<const amrex::Array4<amrex::Real>, 3>& qec,
  const amrex::GpuArray<const amrex::Array4<const amrex::Real>, 3>& a,
  amrex::Array4<amrex::Real> const& pdivu,
  amrex::Array4<const amrex::Real> const& vol,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& del,
  const amrex::Real dt,
  const int ppm_type,
  const int plm_iorder,
  const bool use_flattening,
  const bool use_hybrid_weno,
  const int weno_scheme);

void pc_umeth_eb_3D(
  amrex::Box const& bx_to_fill,
  const int* bclo,
  const int* bchi,
  const int* domlo,
  const int* domhi,
  amrex::Array4<const amrex::Real> const& q,
  amrex::Array4<const amrex::Real> const& qaux,
  amrex::Array4<const amrex::Real> const& srcQ,
  const amrex::GpuArray<const amrex::Array4<amrex::Real>, 3>& flx,
  const amrex::GpuArray<const amrex::Array4<amrex::Real>, 3>& qec,
  const amrex::GpuArray<const amrex::Array4<const amrex::Real>, 3>& ap,
  amrex::Array4<amrex::EBCellFlag const> const& flag_arr,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> del,
  const amrex::Real dt,
  const int ppm_type,
  const bool use_flattening,
  const int plm_iorder);

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
  const amrex::GpuArray<const amrex::Array4<amrex::Real>, 2>& flx,
  const amrex::GpuArray<const amrex::Array4<amrex::Real>, 2>& qec,
  const amrex::GpuArray<const amrex::Array4<const amrex::Real>, 2>& a,
  amrex::Array4<amrex::Real> const& pdivu,
  amrex::Array4<const amrex::Real> const& vol,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& del,
  const amrex::Real dt,
  const int ppm_type,
  const int plm_iorder,
  const bool use_flattening,
  const bool use_hybrid_weno,
  const int weno_scheme);

void pc_umeth_eb_2D(
  amrex::Box const& bx_to_fill,
  const int* bclo,
  const int* bchi,
  const int* domlo,
  const int* domhi,
  amrex::Array4<const amrex::Real> const& q,
  amrex::Array4<const amrex::Real> const& qaux,
  amrex::Array4<const amrex::Real> const& srcQ,
  const amrex::GpuArray<const amrex::Array4<amrex::Real>, 2>& flx,
  const amrex::GpuArray<const amrex::Array4<amrex::Real>, 2>& qec,
  const amrex::GpuArray<const amrex::Array4<const amrex::Real>, 2>& ap,
  amrex::Array4<amrex::EBCellFlag const> const& flag_arr,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> del,
  const amrex::Real dt,
  const int ppm_type,
  const bool use_flattening,
  const int plm_iorder);
#endif

#endif
