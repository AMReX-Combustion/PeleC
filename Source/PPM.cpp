#include "Godunov.H"
#include "PPM.H"
#include "WENO.H"

void
trace_ppm(
  const amrex::Box& bx,
  const int idir,
  amrex::Array4<amrex::Real const> const& q_arr,
  amrex::Array4<amrex::Real const> const& /*srcQ*/,
  amrex::Array4<amrex::Real> const& qm,
  amrex::Array4<amrex::Real> const& qp,
  const amrex::Box& vbx,
  const amrex::Real dt,
  const amrex::Real* dx,
  const bool use_flattening,
  const bool use_hybrid_weno,
  const int weno_scheme)
{
  // here, lo and hi are the range we loop over -- this can include ghost cells
  // vlo and vhi are the bounds of the valid box (no ghost cells)

  //
  // rho : mass density
  // u, v, w : velocities
  // p : gas (hydro) pressure
  // ptot : total pressure (note for pure hydro, this is
  //        just the gas pressure)
  // rhoe_g : gas specific internal energy
  // cgas : sound speed for just the gas contribution
  // cc : total sound speed
  // h_g : gas specific enthalpy / cc**2
  //
  // for pure hydro, we will only consider:
  //    rho, u, v, w, ptot, rhoe_g, cc, h_g
  // amrex::Real hdt = 0.5 * dt;
  amrex::Real dtdx = dt / dx[idir];

  // auto lo = bx.loVect3d();
  // auto hi = bx.hiVect3d();

  auto vlo = vbx.loVect3d();
  auto vhi = vbx.hiVect3d();

  // This does the characteristic tracing to build the interface
  // states using the normal predictor only (no transverse terms).
  //
  // For each zone, we construct Im and Ip arrays -- these are the averages
  // of the various primitive state variables under the parabolic
  // interpolant over the region swept out by one of the 3 different
  // characteristic waves.
  //
  // Im is integrating to the left interface of the current zone
  // (which will be used to build the right ("p") state at that interface)
  // and Ip is integrating to the right interface of the current zone
  // (which will be used to build the left ("m") state at that interface).
  //
  //
  // The choice of reference state is designed to minimize the
  // effects of the characteristic projection.  We subtract the I's
  // off of the reference state, project the quantity such that it is
  // in terms of the characteristic varaibles, and then add all the
  // jumps that are moving toward the interface to the reference
  // state to get the full state on that interface.

  int QUN = 0;
  int QUT = 0;
  int QUTT = 0;

  if (idir == 0) {
    QUN = QU;
    QUT = QV;
    QUTT = QW;
  } else if (idir == 1) {
    QUN = QV;
    QUT = QW;
    QUTT = QU;
  } else if (idir == 2) {
    QUN = QW;
    QUT = QU;
    QUTT = QV;
  }

  // Trace to left and right edges using upwind PPM
  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    // amrex::Real rho = q_arr(i, j, k, QRHO);

    const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
    const amrex::IntVect ivm2(
      iv - 2 * amrex::IntVect::TheDimensionVector(idir));
    const amrex::IntVect ivm1(
      iv - 1 * amrex::IntVect::TheDimensionVector(idir));
    const amrex::IntVect ivp1(
      iv + 1 * amrex::IntVect::TheDimensionVector(idir));
    const amrex::IntVect ivp2(
      iv + 2 * amrex::IntVect::TheDimensionVector(idir));

    auto eos = pele::physics::PhysicsType::eos();

    amrex::Real massfrac[NUM_SPECIES];
    for (int species = 0; species < NUM_SPECIES; ++species) {
      massfrac[species] = q_arr(iv, species + QFS);
    }

    amrex::Real cc = 0;
    eos.RPY2Cs(q_arr(iv, QRHO), q_arr(iv, QPRES), massfrac, cc);

    amrex::Real un = q_arr(iv, QUN);

    // do the parabolic reconstruction and compute the
    // integrals under the characteristic waves

    amrex::Real flat = 1.0;
    // Calculate flattening in-place
    if (use_flattening) {
      for (int dir_flat = 0; dir_flat < AMREX_SPACEDIM; dir_flat++) {
        flat = amrex::min<amrex::Real>(flat, flatten(i, j, k, dir_flat, q_arr));
      }
    }

    amrex::Real Ip[QVAR][3];
    amrex::Real Im[QVAR][3];

    for (int n = 0; n < QVAR; n++) {
      if (use_hybrid_weno && ((weno_scheme == 0) || (weno_scheme == 1))) {

        amrex::Real s_weno5[5];
        s_weno5[0] = q_arr(ivm2, n);
        s_weno5[1] = q_arr(ivm1, n);
        s_weno5[2] = q_arr(iv, n);
        s_weno5[3] = q_arr(ivp1, n);
        s_weno5[4] = q_arr(ivp2, n);

        amrex::Real sm = 0.0;
        amrex::Real sp = 0.0;
        if (weno_scheme == 0) {
          weno_reconstruct_5js(s_weno5, sm, sp);
        } else if (weno_scheme == 1) {
          weno_reconstruct_5z(s_weno5, sm, sp);
        }
        ppm_int_profile(sm, sp, s_weno5[2], un, cc, dtdx, Ip[n], Im[n]);

      } else if (use_hybrid_weno && weno_scheme == 2) {

        amrex::Real s_weno7[7];
        const amrex::IntVect ivm3(
          iv - 3 * amrex::IntVect::TheDimensionVector(idir));
        const amrex::IntVect ivp3(
          iv + 3 * amrex::IntVect::TheDimensionVector(idir));
        s_weno7[0] = q_arr(ivm3, n);
        s_weno7[1] = q_arr(ivm2, n);
        s_weno7[2] = q_arr(ivm1, n);
        s_weno7[3] = q_arr(iv, n);
        s_weno7[4] = q_arr(ivp1, n);
        s_weno7[5] = q_arr(ivp2, n);
        s_weno7[6] = q_arr(ivp3, n);

        amrex::Real sm;
        amrex::Real sp;
        weno_reconstruct_7z(s_weno7, sm, sp);
        ppm_int_profile(sm, sp, s_weno7[3], un, cc, dtdx, Ip[n], Im[n]);

      } else if (use_hybrid_weno && weno_scheme == 3) {

        amrex::Real s_weno3[3];
        if (idir == 0) {
          s_weno3[0] = q_arr(ivm1, n);
          s_weno3[1] = q_arr(iv, n);
          s_weno3[2] = q_arr(ivp1, n);
        } else if (idir == 1) {
          s_weno3[0] = q_arr(i, j - 1, k, n);
          s_weno3[1] = q_arr(iv, n);
          s_weno3[2] = q_arr(i, j + 1, k, n);
        } else {
          s_weno3[0] = q_arr(i, j, k - 1, n);
          s_weno3[1] = q_arr(iv, n);
          s_weno3[2] = q_arr(i, j, k + 1, n);
        }

        amrex::Real sm;
        amrex::Real sp;
        weno_reconstruct_3z(s_weno3, sm, sp);
        ppm_int_profile(sm, sp, s_weno3[1], un, cc, dtdx, Ip[n], Im[n]);

        // ORIGINAL PPM
      } else {

        amrex::Real s[5];
        s[im2] = q_arr(ivm2, n);
        s[im1] = q_arr(ivm1, n);
        s[i0] = q_arr(iv, n);
        s[ip1] = q_arr(ivp1, n);
        s[ip2] = q_arr(ivp2, n);
        amrex::Real sm;
        amrex::Real sp;
        ppm_reconstruct(s, flat, sm, sp);
        ppm_int_profile(sm, sp, s[2], un, cc, dtdx, Ip[n], Im[n]);
      }
    }

    // PeleC does source term tracing in pc_transx, pc_transy, and
    // pc_transz. So to be consistent we remove the source term
    // tracing here. However, Nyx and Castro do the source term
    // tracing here instead of in the trans routines. If we wanted
    // to do the tracing here, we would have to 1) remove the tracing
    // in the trans routines AND 2) add the tracing in the pc_plm_x,
    // pc_plm_y, pc_plm_z routines.
    //
    // To do the tracing here: Uncomment the chunk below and
    // anything that uses Im_src and Ip_src.

    // // source terms
    // amrex::Real Ip_src[QVAR][3];
    // amrex::Real Im_src[QVAR][3];

    // for (int n = 0; n < QVAR; n++) {
    //   s[im2] = srcQ(ivm2, n);
    //   s[im1] = srcQ(ivm1, n);
    //   s[i0] = srcQ(iv, n);
    //   s[ip1] = srcQ(ivp1, n);
    //   s[ip2] = srcQ(ivp2, n);
    //   ppm_reconstruct(s, flat, sm, sp);
    //   ppm_int_profile(sm, sp, s[i0], un, cc, dtdx, Ip_src[n], Im_src[n]);
    // }

    for (int n = QFS; n < NUM_SPECIES + QFS; n++) {
      // Plus state on face i
      if (
        (idir == 0 && i >= vlo[0]) || (idir == 1 && j >= vlo[1]) ||
        (idir == 2 && k >= vlo[2])) {

        // We have
        //
        // q_l = q_ref - Proj{(q_ref - I)}
        //
        // and Proj{} represents the characteristic projection.
        // But for these, there is only 1-wave that matters, the u
        // wave, so no projection is needed.  Since we are not
        // projecting, the reference state doesn't matter

        qp(iv, n) = Im[n][1];
      }

      // Minus state on face i+1
      if (
        (idir == 0 && i <= vhi[0]) || (idir == 1 && j <= vhi[1]) ||
        (idir == 2 && k <= vhi[2])) {
        qm(ivp1, n) = Ip[n][1];
      }
    }

    // plus state on face i
    if (
      (idir == 0 && i >= vlo[0]) || (idir == 1 && j >= vlo[1]) ||
      (idir == 2 && k >= vlo[2])) {

      // Set the reference state
      // This will be the fastest moving state to the left --
      // this is the method that Miller & Colella and Colella &
      // Woodward use
      amrex::Real rho_ref = Im[QRHO][0];
      amrex::Real un_ref = Im[QUN][0];

      amrex::Real p_ref = Im[QPRES][0];
      amrex::Real rhoe_g_ref = Im[QREINT][0];

      rho_ref = amrex::max<amrex::Real>(
        rho_ref, std::numeric_limits<amrex::Real>::min());
      amrex::Real rho_ref_inv = 1.0 / rho_ref;
      p_ref =
        amrex::max<amrex::Real>(p_ref, std::numeric_limits<amrex::Real>::min());

      amrex::Real massfrac_ref[NUM_SPECIES];
      for (int species = 0; species < NUM_SPECIES; ++species) {
        massfrac_ref[species] = Im[species + QFS][0];
      }

      // For tracing
      amrex::Real cc_ref = 0;
      eos.RPY2Cs(rho_ref, p_ref, massfrac_ref, cc_ref);
      amrex::Real csq_ref = cc_ref * cc_ref;
      amrex::Real cc_ref_inv = 1.0 / cc_ref;
      amrex::Real h_g_ref = (p_ref + rhoe_g_ref) * rho_ref_inv;

      // *m are the jumps carried by un-c
      // *p are the jumps carried by un+c

      // Note: for the transverse velocities, the jump is carried
      //       only by the u wave (the contact)

      // we also add the sources here so they participate in the tracing
      amrex::Real dum = un_ref - Im[QUN][0] /*- hdt * Im_src[QUN][0]*/;
      amrex::Real dptotm = p_ref - Im[QPRES][0] /*- hdt * Im_src[QPRES][0]*/;

      amrex::Real drho = rho_ref - Im[QRHO][1] /*- hdt * Im_src[QRHO][1]*/;
      amrex::Real dptot = p_ref - Im[QPRES][1] /*- hdt * Im_src[QPRES][1]*/;
      amrex::Real drhoe_g =
        rhoe_g_ref - Im[QREINT][1] /*- hdt * Im_src[QREINT][1]*/;

      amrex::Real dup = un_ref - Im[QUN][2] /*- hdt * Im_src[QUN][2]*/;
      amrex::Real dptotp = p_ref - Im[QPRES][2] /*- hdt * Im_src[QPRES][2]*/;

      // {rho, u, p, (rho e)} eigensystem

      // These are analogous to the beta's from the original PPM
      // paper (except we work with rho instead of tau).  This is
      // simply (l . dq), where dq = qref - I(q)

      amrex::Real alpham =
        0.5 * (dptotm * rho_ref_inv * cc_ref_inv - dum) * rho_ref * cc_ref_inv;
      amrex::Real alphap =
        0.5 * (dptotp * rho_ref_inv * cc_ref_inv + dup) * rho_ref * cc_ref_inv;
      amrex::Real alpha0r = drho - dptot / csq_ref;
      amrex::Real alpha0e_g = drhoe_g - dptot * h_g_ref / csq_ref;

      alpham = un - cc > 0.0 ? 0.0 : -alpham;
      alphap = un + cc > 0.0 ? 0.0 : -alphap;
      alpha0r = un > 0.0 ? 0.0 : -alpha0r;
      alpha0e_g = un > 0.0 ? 0.0 : -alpha0e_g;

      // The final interface states are just
      // q_s = q_ref - sum(l . dq) r
      // note that the a{mpz}right as defined above have the minus already
      qp(iv, QRHO) = amrex::max<amrex::Real>(
        std::numeric_limits<amrex::Real>::min(),
        rho_ref + alphap + alpham + alpha0r);
      qp(iv, QUN) = un_ref + (alphap - alpham) * cc_ref * rho_ref_inv;
      // qp(i,j,k,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g_ref +
      // alpha0e_g;
      qp(iv, QPRES) = amrex::max<amrex::Real>(
        std::numeric_limits<amrex::Real>::min(),
        p_ref + (alphap + alpham) * csq_ref);

      // Transverse velocities -- there's no projection here, so we
      // don't need a reference state.  We only care about the state
      // traced under the middle wave

      // Recall that I already takes the limit of the parabola
      // in the event that the wave is not moving toward the
      // interface
      qp(iv, QUT) = Im[QUT][1] /*+ hdt * Im_src[QUT][1]*/;
      qp(iv, QUTT) = Im[QUTT][1] /*+ hdt * Im_src[QUTT][1]*/;

      // This allows the (rho e) to take advantage of (pressure > small_pres)
      amrex::Real eint = 0;
      amrex::Real massfrac_p[NUM_SPECIES];
      for (int species = 0; species < NUM_SPECIES; ++species) {
        massfrac_p[species] = qp(iv, species + QFS);
      }
      eos.RYP2E(qp(iv, QRHO), massfrac_p, qp(iv, QPRES), eint);
      qp(iv, QREINT) = qp(iv, QRHO) * eint;
    }

    // minus state on face i + 1
    if (
      (idir == 0 && i <= vhi[0]) || (idir == 1 && j <= vhi[1]) ||
      (idir == 2 && k <= vhi[2])) {

      // Set the reference state
      // This will be the fastest moving state to the right
      amrex::Real rho_ref = Ip[QRHO][2];
      amrex::Real un_ref = Ip[QUN][2];

      amrex::Real p_ref = Ip[QPRES][2];
      amrex::Real rhoe_g_ref = Ip[QREINT][2];

      rho_ref = amrex::max<amrex::Real>(
        rho_ref, std::numeric_limits<amrex::Real>::min());
      amrex::Real rho_ref_inv = 1.0 / rho_ref;
      p_ref =
        amrex::max<amrex::Real>(p_ref, std::numeric_limits<amrex::Real>::min());

      amrex::Real massfrac_ref[NUM_SPECIES];
      for (int species = 0; species < NUM_SPECIES; ++species) {
        massfrac_ref[species] = Ip[species + QFS][2];
      }

      // For tracing
      amrex::Real cc_ref = 0;
      eos.RPY2Cs(rho_ref, p_ref, massfrac_ref, cc_ref);
      amrex::Real csq_ref = cc_ref * cc_ref;
      amrex::Real cc_ref_inv = 1.0 / cc_ref;
      amrex::Real h_g_ref = (p_ref + rhoe_g_ref) * rho_ref_inv;

      // *m are the jumps carried by u-c
      // *p are the jumps carried by u+c

      amrex::Real dum = un_ref - Ip[QUN][0] /*- hdt * Ip_src[QUN][0]*/;
      amrex::Real dptotm = p_ref - Ip[QPRES][0] /*- hdt * Ip_src[QPRES][0]*/;

      amrex::Real drho = rho_ref - Ip[QRHO][1] /*- hdt * Ip_src[QRHO][1]*/;
      amrex::Real dptot = p_ref - Ip[QPRES][1] /*- hdt * Ip_src[QPRES][1]*/;
      amrex::Real drhoe_g =
        rhoe_g_ref - Ip[QREINT][1] /*- hdt * Ip_src[QREINT][1]*/;

      amrex::Real dup = un_ref - Ip[QUN][2] /*- hdt * Ip_src[QUN][2]*/;
      amrex::Real dptotp = p_ref - Ip[QPRES][2] /*- hdt * Ip_src[QPRES][2]*/;

      // {rho, u, p, (rho e)} eigensystem

      // These are analogous to the beta's from the original PPM
      // paper (except we work with rho instead of tau).  This is
      // simply (l . dq), where dq = qref - I(q)

      amrex::Real alpham =
        0.5 * (dptotm * rho_ref_inv * cc_ref_inv - dum) * rho_ref * cc_ref_inv;
      amrex::Real alphap =
        0.5 * (dptotp * rho_ref_inv * cc_ref_inv + dup) * rho_ref * cc_ref_inv;
      amrex::Real alpha0r = drho - dptot / csq_ref;
      amrex::Real alpha0e_g = drhoe_g - dptot * h_g_ref / csq_ref;

      alpham = un - cc > 0.0 ? -alpham : 0.0;
      alphap = un + cc > 0.0 ? -alphap : 0.0;
      alpha0r = un > 0.0 ? -alpha0r : 0.0;
      alpha0e_g = un > 0.0 ? -alpha0e_g : 0.0;

      // The final interface states are just
      // q_s = q_ref - sum (l . dq) r
      // note that the a{mpz}left as defined above have the minus already
      qm(ivp1, QRHO) = amrex::max<amrex::Real>(
        std::numeric_limits<amrex::Real>::min(),
        rho_ref + alphap + alpham + alpha0r);
      qm(ivp1, QUN) = un_ref + (alphap - alpham) * cc_ref * rho_ref_inv;
      // qm(i+1,j,k,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g_ref +
      // alpha0e_g;
      qm(ivp1, QPRES) = amrex::max<amrex::Real>(
        std::numeric_limits<amrex::Real>::min(),
        p_ref + (alphap + alpham) * csq_ref);

      // transverse velocities
      qm(ivp1, QUT) = Ip[QUT][1] /*+ hdt * Ip_src[QUT][1]*/;
      qm(ivp1, QUTT) = Ip[QUTT][1] /*+ hdt * Ip_src[QUTT][1]*/;

      // This allows the (rho e) to take advantage of (pressure >
      // small_pres)
      amrex::Real eint = 0;
      amrex::Real massfrac_m[NUM_SPECIES];
      for (int species = 0; species < NUM_SPECIES; ++species) {
        massfrac_m[species] = qm(ivp1, species + QFS);
      }
      eos.RYP2E(qm(ivp1, QRHO), massfrac_m, qm(ivp1, QPRES), eint);
      qm(ivp1, QREINT) = qm(ivp1, QRHO) * eint;
    }
  });
}
