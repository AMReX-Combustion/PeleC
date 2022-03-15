#ifndef _RIEMANN_H_
#define _RIEMANN_H_
#include "PeleC.H"
#include "PelePhysics.H"

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
riemann(
  const amrex::Real rl,
  const amrex::Real ul,
  const amrex::Real vl,
  const amrex::Real v2l,
  const amrex::Real pl,
  const amrex::Real spl[NUM_SPECIES],
  const amrex::Real rr,
  const amrex::Real ur,
  const amrex::Real vr,
  const amrex::Real v2r,
  const amrex::Real pr,
  const amrex::Real spr[NUM_SPECIES],
  const int bc_test_val,
  const amrex::Real cav,
  amrex::Real& ustar,
  amrex::Real& uflx_rho,
  amrex::Real uflx_rhoY[NUM_SPECIES],
  amrex::Real& uflx_u,
  amrex::Real& uflx_v,
  amrex::Real& uflx_w,
  amrex::Real& uflx_eden,
  amrex::Real& uflx_eint,
  amrex::Real& qint_iu,
  amrex::Real& qint_iv1,
  amrex::Real& qint_iv2,
  amrex::Real& qint_gdpres,
  amrex::Real& qint_gdgame)
{
  const amrex::Real wsmall = std::numeric_limits<amrex::Real>::min();

  auto eos = pele::physics::PhysicsType::eos();

  amrex::Real gdnv_state_massfrac[NUM_SPECIES];
  for (int n = 0; n < NUM_SPECIES; n++) {
    gdnv_state_massfrac[n] = spl[n];
  }
  amrex::Real cl = 0.0;
  eos.RPY2Cs(rl, pl, gdnv_state_massfrac, cl);

  for (int n = 0; n < NUM_SPECIES; n++) {
    gdnv_state_massfrac[n] = spr[n];
  }
  amrex::Real cr = 0.0;
  eos.RPY2Cs(rr, pr, gdnv_state_massfrac, cr);

  const amrex::Real wl = amrex::max<amrex::Real>(wsmall, cl * rl);
  const amrex::Real wr = amrex::max<amrex::Real>(wsmall, cr * rr);
  const amrex::Real pstar = amrex::max<amrex::Real>(
    std::numeric_limits<amrex::Real>::min(),
    ((wr * pl + wl * pr) + wl * wr * (ul - ur)) / (wl + wr));
  ustar = ((wl * ul + wr * ur) + (pl - pr)) / (wl + wr);

  bool mask = ustar > 0.0;
  amrex::Real ro = 0.0;
  amrex::Real rspo[NUM_SPECIES];
  for (int n = 0; n < NUM_SPECIES; n++) {
    rspo[n] = mask ? rl * spl[n] : rr * spr[n];
    ro += rspo[n];
  }
  amrex::Real uo = mask ? ul : ur;
  amrex::Real po = mask ? pl : pr;

  mask =
    amrex::Math::abs(ustar) < constants::smallu() * 0.5 *
                                (amrex::Math::abs(ul) + amrex::Math::abs(ur)) ||
    ustar == 0.0;
  ustar = mask ? 0.0 : ustar;
  ro = 0.0;
  for (int n = 0; n < NUM_SPECIES; n++) {
    rspo[n] = mask ? 0.5 * (rl * spl[n] + rr * spr[n]) : rspo[n];
    ro += rspo[n];
  }
  uo = mask ? 0.5 * (ul + ur) : uo;
  po = mask ? 0.5 * (pl + pr) : po;

  amrex::Real gdnv_state_rho = ro;
  amrex::Real gdnv_state_p = po;
  for (int n = 0; n < NUM_SPECIES; n++) {
    gdnv_state_massfrac[n] = rspo[n] / ro;
  }
  amrex::Real gdnv_state_e;
  eos.RYP2E(gdnv_state_rho, gdnv_state_massfrac, gdnv_state_p, gdnv_state_e);
  amrex::Real co;
  eos.RPY2Cs(gdnv_state_rho, gdnv_state_p, gdnv_state_massfrac, co);

  const amrex::Real drho = (pstar - po) / (co * co);
  amrex::Real rstar = 0.0;
  amrex::Real rspstar[NUM_SPECIES];
  for (int n = 0; n < NUM_SPECIES; n++) {
    const amrex::Real spon = rspo[n] / ro;
    rspstar[n] = amrex::max<amrex::Real>(0.0, rspo[n] + drho * spon);
    rstar += rspstar[n];
  }
  gdnv_state_rho = rstar;
  gdnv_state_p = pstar;
  for (int n = 0; n < NUM_SPECIES; n++) {
    gdnv_state_massfrac[n] = rspstar[n] / rstar;
  }
  eos.RYP2E(gdnv_state_rho, gdnv_state_massfrac, gdnv_state_p, gdnv_state_e);
  amrex::Real cstar;
  eos.RPY2Cs(gdnv_state_rho, gdnv_state_p, gdnv_state_massfrac, cstar);

  const amrex::Real sgnm = amrex::Math::copysign(1.0, ustar);

  amrex::Real spout = co - sgnm * uo;
  amrex::Real spin = cstar - sgnm * ustar;
  const amrex::Real ushock = 0.5 * (spin + spout);

  mask = pstar < po;
  spout = mask ? spout : ushock;
  spin = mask ? spin : ushock;

  const amrex::Real scr =
    (spout == spin) ? constants::small_num() * cav : spout - spin;
  const amrex::Real frac = amrex::max<amrex::Real>(
    0.0, amrex::min<amrex::Real>(1.0, (1.0 + (spout + spin) / scr) * 0.5));

  mask = ustar > 0.0;
  qint_iv1 = mask ? vl : vr;
  qint_iv2 = mask ? v2l : v2r;

  mask = (ustar == 0.0);
  qint_iv1 = mask ? 0.5 * (vl + vr) : qint_iv1;
  qint_iv2 = mask ? 0.5 * (v2l + v2r) : qint_iv2;
  amrex::Real rgd = 0.0;
  amrex::Real rspgd[NUM_SPECIES];
  for (int n = 0; n < NUM_SPECIES; n++) {
    rspgd[n] = frac * rspstar[n] + (1.0 - frac) * rspo[n];
    rgd += rspgd[n];
  }
  qint_iu = frac * ustar + (1.0 - frac) * uo;
  qint_gdpres = frac * pstar + (1.0 - frac) * po;
  gdnv_state_rho = rgd;
  gdnv_state_p = qint_gdpres;
  for (int n = 0; n < NUM_SPECIES; n++) {
    gdnv_state_massfrac[n] = rspgd[n] / rgd;
  }
  eos.RYP2E(gdnv_state_rho, gdnv_state_massfrac, gdnv_state_p, gdnv_state_e);

  mask = (spout < 0.0);
  rgd = 0.0;
  for (int n = 0; n < NUM_SPECIES; n++) {
    rspgd[n] = mask ? rspo[n] : rspgd[n];
    rgd += rspgd[n];
  }
  qint_iu = mask ? uo : qint_iu;
  qint_gdpres = mask ? po : qint_gdpres;

  mask = (spin >= 0.0);
  rgd = 0.0;
  for (int n = 0; n < NUM_SPECIES; n++) {
    rspgd[n] = mask ? rspstar[n] : rspgd[n];
    rgd += rspgd[n];
  }
  qint_iu = mask ? ustar : qint_iu;
  qint_gdpres = mask ? pstar : qint_gdpres;

  gdnv_state_rho = rgd;
  gdnv_state_p = qint_gdpres;
  for (int n = 0; n < NUM_SPECIES; n++) {
    gdnv_state_massfrac[n] = rspgd[n] / rgd;
  }
  eos.RYP2E(gdnv_state_rho, gdnv_state_massfrac, gdnv_state_p, gdnv_state_e);
  amrex::Real regd = gdnv_state_rho * gdnv_state_e;

  qint_gdgame = qint_gdpres / regd + 1.0;
  qint_iu = bc_test_val * qint_iu;
  uflx_rho = rgd * qint_iu;
  for (int n = 0; n < NUM_SPECIES; n++) {
    uflx_rhoY[n] = rspgd[n] * qint_iu;
  }
  uflx_u = uflx_rho * qint_iu + qint_gdpres;
  uflx_v = uflx_rho * qint_iv1;
  uflx_w = uflx_rho * qint_iv2;
  const amrex::Real rhoetot =
    regd +
    0.5 * rgd * (qint_iu * qint_iu + qint_iv1 * qint_iv1 + qint_iv2 * qint_iv2);
  uflx_eden = qint_iu * (rhoetot + qint_gdpres);
  uflx_eint = qint_iu * regd;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
laxfriedrich_flux(
  const amrex::Real rl,
  const amrex::Real ul,
  const amrex::Real vl,
  const amrex::Real v2l,
  const amrex::Real pl,
  const amrex::Real spl[NUM_SPECIES],
  const amrex::Real rr,
  const amrex::Real ur,
  const amrex::Real vr,
  const amrex::Real v2r,
  const amrex::Real pr,
  const amrex::Real spr[NUM_SPECIES],
  const int /*unused*/,
  const amrex::Real /*unused*/,
  amrex::Real& ustar,
  amrex::Real& maxeigval,
  amrex::Real& uflx_rho,
  amrex::Real uflx_rhoY[NUM_SPECIES],
  amrex::Real& uflx_u,
  amrex::Real& uflx_v,
  amrex::Real& uflx_w,
  amrex::Real& uflx_eden,
  amrex::Real& uflx_eint)
{

  const amrex::Real wsmall = std::numeric_limits<amrex::Real>::min();

  auto eos = pele::physics::PhysicsType::eos();

  amrex::Real cl = 0.0;
  amrex::Real mfrac_l[NUM_SPECIES];
  amrex::Real mfrac_r[NUM_SPECIES];
  for (int n = 0; n < NUM_SPECIES; n++) {
    mfrac_l[n] = spl[n];
    mfrac_r[n] = spr[n];
  }
  eos.RPY2Cs(rl, pl, mfrac_l, cl);

  amrex::Real cr = 0.0;
  eos.RPY2Cs(rr, pr, mfrac_r, cr);

  amrex::Real el = 0.0;
  eos.RYP2E(rl, mfrac_l, pl, el);

  amrex::Real er = 0.0;
  eos.RYP2E(rr, mfrac_r, pr, er);

  // rho e total
  amrex::Real r_elt = rl * (el + 0.5 * (ul * ul + vl * vl + v2l * v2l));
  amrex::Real r_ert = rr * (er + 0.5 * (ur * ur + vr * vr + v2r * v2r));

  const amrex::Real wl = amrex::max<amrex::Real>(wsmall, cl * rl);
  const amrex::Real wr = amrex::max<amrex::Real>(wsmall, cr * rr);
  ustar = ((wl * ul + wr * ur) + (pl - pr)) / (wl + wr);

  bool mask =
    amrex::Math::abs(ustar) < constants::smallu() * 0.5 *
                                (amrex::Math::abs(ul) + amrex::Math::abs(ur)) ||
    ustar == 0.0;
  ustar = mask ? 0.0 : ustar;

  amrex::Real chalf = 0.5 * (cr + cl);
  amrex::Real max_wavespd = amrex::Math::abs(ustar) + chalf;

  maxeigval = max_wavespd;

  // density and mass frac eqn
  uflx_rho = 0.0;
  for (int n = 0; n < NUM_SPECIES; n++) {
    uflx_rhoY[n] = 0.5 * ((rl * ul * spl[n] + rr * ur * spr[n]) -
                          max_wavespd * (rr * spr[n] - rl * spl[n]));
    uflx_rho += uflx_rhoY[n];
  }

  // three momentum eqns
  uflx_u = 0.5 * ((rl * ul * ul + pl + rr * ur * ur + pr) -
                  max_wavespd * (rr * ur - rl * ul));
  uflx_v =
    0.5 * ((rl * vl * ul + rr * vr * ur) - max_wavespd * (rr * vr - rl * vl));
  uflx_w = 0.5 * ((rl * v2l * ul + rr * v2r * ur) -
                  max_wavespd * (rr * v2r - rl * v2l));

  // total energy eqn
  uflx_eden = 0.5 * (((r_elt + pl) * ul + (r_ert + pr) * ur) -
                     max_wavespd * (r_ert - r_elt));

  // not really used
  uflx_eint =
    0.5 * ((rl * el * ul + rr * er * ur) - max_wavespd * (rr * er - rl * el));
}
#endif
