#include "EB.H"

void
pc_fill_sv_ebg(
  const amrex::Box& bx,
  const int Nebg,
  const amrex::Array4<const amrex::Real>& vfrac,
  const amrex::Array4<const amrex::Real>& bcent,
  AMREX_D_DECL(
    const amrex::Array4<const amrex::Real>& apx,
    const amrex::Array4<const amrex::Real>& apy,
    const amrex::Array4<const amrex::Real>& apz),
  EBBndryGeom* ebg)
{
  const auto lo = amrex::lbound(bx);
  const auto hi = amrex::ubound(bx);
  amrex::ParallelFor(Nebg, [=] AMREX_GPU_DEVICE(int L) {
    const auto& iv = ebg[L].iv;
    if (is_inside(iv, lo, hi)) {
      const amrex::Real axm = apx(iv);
      const amrex::Real axp = apx(iv + amrex::IntVect::TheDimensionVector(0));
      const amrex::Real aym = apy(iv);
      const amrex::Real ayp = apy(iv + amrex::IntVect::TheDimensionVector(1));
#if AMREX_SPACEDIM > 2
      const amrex::Real azm = apz(iv);
      const amrex::Real azp = apz(iv + amrex::IntVect::TheDimensionVector(2));
#endif
      const amrex::Real apnorm = std::sqrt(AMREX_D_TERM(
        (axm - axp) * (axm - axp), +(aym - ayp) * (aym - ayp),
        +(azm - azp) * (azm - azp)));
      if (apnorm == 0.0) {
        amrex::Abort("pc_fill_sv_ebg: zero apnorm");
      }
      const amrex::Real apnorminv = -1.0 / apnorm;
      // pointing to the wall
      AMREX_D_TERM(ebg[L].eb_normal[0] = (axm - axp) * apnorminv;
                   , ebg[L].eb_normal[1] = (aym - ayp) * apnorminv;
                   , ebg[L].eb_normal[2] = (azm - azp) * apnorminv;)
      ebg[L].eb_area = apnorm;
      AMREX_D_TERM(ebg[L].eb_centroid[0] = bcent(iv, 0);
                   , ebg[L].eb_centroid[1] = bcent(iv, 1);
                   , ebg[L].eb_centroid[2] = bcent(iv, 2);)
      ebg[L].eb_vfrac = vfrac(iv);
    }
  });
}

// see Johansen and Collela paper
// Johansen, H., & Colella, P. (1998). A Cartesian grid embedded boundary method
// for Poisson's equation on irregular domains.
// Journal of Computational Physics, 147(1), 60-85
//
// Johansen and Collela also say
//"By constructing the gradients in this fashion, we impose one more constraint
// on the discretization of the domain: the interpolation stencil must not reach
// into cells with zero volume. For the quadratic gradient stencil, this may
// imply certain constraints on the discretization of the domain. However, the
// fact that a zero-volume cell is within two rows of another partial cell would
// indicate that the local boundary is substantially under-resolved. Such
// domains are more appropriately treated with adaptive mesh refinement"
//
//"which means that there is a chance that this stencil may
// dip into covered cells":Hari S

void
pc_fill_bndry_grad_stencil_quadratic(
  const amrex::Box& bx,
  const amrex::Real dx,
  const int /*Nebg*/,
  const EBBndryGeom* ebg,
  const int Nsten,
  EBBndrySten* grad_stencil)
{
  const auto lo = amrex::lbound(bx);
  const auto hi = amrex::ubound(bx);
  const amrex::Real area = std::pow(dx, AMREX_SPACEDIM - 1);
  const amrex::Real fac = area / dx;

  amrex::ParallelFor(Nsten, [=] AMREX_GPU_DEVICE(int L) {
    if (is_inside(ebg[L].iv, lo, hi)) {
      const amrex::Real n[AMREX_SPACEDIM] = {AMREX_D_DECL(
        ebg[L].eb_normal[0], ebg[L].eb_normal[1], ebg[L].eb_normal[2])};

      int c[AMREX_SPACEDIM] = {0};
      idxsort(n, c);
      const int ivs[AMREX_SPACEDIM] = {
        AMREX_D_DECL(ebg[L].iv[0], ebg[L].iv[1], ebg[L].iv[2])};
      const int s[AMREX_SPACEDIM] = {AMREX_D_DECL(
        (int)amrex::Math::copysign(1.0, n[c[0]]),
        (int)amrex::Math::copysign(1.0, n[c[1]]),
        (int)amrex::Math::copysign(1.0, n[c[2]]))};
      amrex::Real b[AMREX_SPACEDIM] = {AMREX_D_DECL(
        ebg[L].eb_centroid[c[0]] * s[0], ebg[L].eb_centroid[c[1]] * s[1],
        ebg[L].eb_centroid[c[2]] * s[2])};

      // From ivs, move to center of stencil, then move to lower-left of that
      const int baseiv[AMREX_SPACEDIM] = {AMREX_D_DECL(
        ivs[0] + (int)amrex::Math::copysign(1.0, n[0]) - 1,
        ivs[1] + (int)amrex::Math::copysign(1.0, n[1]) - 1,
        ivs[2] + (int)amrex::Math::copysign(1.0, n[2]) - 1)};

      const amrex::Real x[2] = {1.0, 2.0};
      amrex::Real y[2] = {
        b[1] + (x[0] - b[0]) * amrex::Math::abs(n[c[1]] / n[c[0]]),
        b[1] + (x[1] - b[0]) * amrex::Math::abs(n[c[1]] / n[c[0]])};
#if AMREX_SPACEDIM > 2
      amrex::Real z[2] = {
        b[2] + (x[0] - b[0]) * amrex::Math::abs(n[c[2]] / n[c[0]]),
        b[2] + (x[1] - b[0]) * amrex::Math::abs(n[c[2]] / n[c[0]])};
#endif

      int sh[AMREX_SPACEDIM] = {0};
      if (y[0] < 0.0 || y[1] < 0.0) {
        sh[c[1]] = -s[1]; // Slide stencil down to avoid extrapolating, push
                          // up eb, shift down base later
        b[1] += 1;
        y[0] = b[1] + (x[0] - b[0]) * amrex::Math::abs(n[c[1]] / n[c[0]]);
        y[1] = b[1] + (x[1] - b[0]) * amrex::Math::abs(n[c[1]] / n[c[0]]);
      }
#if AMREX_SPACEDIM > 2
      if (z[0] < 0.0 || z[1] < 0.0) {
        sh[c[2]] = -s[2]; // Slide stencil down to avoid extrapolating, push
                          // up eb, shift down base later
        b[2] += 1;
        z[0] = b[2] + (x[0] - b[0]) * amrex::Math::abs(n[c[2]] / n[c[0]]);
        z[1] = b[2] + (x[1] - b[0]) * amrex::Math::abs(n[c[2]] / n[c[0]]);
      }
#endif
      const amrex::Real d[2] = {
        std::sqrt(AMREX_D_TERM(
          (x[0] - b[0]) * (x[0] - b[0]), +(y[0] - b[1]) * (y[0] - b[1]),
          +(z[0] - b[2]) * (z[0] - b[2]))),
        std::sqrt(AMREX_D_TERM(
          (x[1] - b[0]) * (x[1] - b[0]), +(y[1] - b[1]) * (y[1] - b[1]),
          +(z[1] - b[2]) * (z[1] - b[2])))};

      amrex::Real sten AMREX_D_TERM([3],[3],[3]) = AMREX_D_TERM({,{,{) 0.0 AMREX_D_TERM(},},});
      // The two intersections, that are d(1) and d(2) away from the eb
      // centroid, are both in y-z planes, bounded in (0:2)x(0:2) in normalized
      // coordinates For point m, we interpolate z=0,1,2 lines, to (y(m),0),
      // (y(m),1) and (y(m),2), and then interpolate along y=y(m) to (y(m),z(m))
      amrex::Real cy[3];
#if AMREX_SPACEDIM > 2
      amrex::Real cz[3];
#endif
      for (int m = 0; m < 2; m++) {
        cy[0] = 0.5 * (y[m] - 1.0) * (y[m] - 2.0);
        cy[1] = -y[m] * (y[m] - 2.0);
        cy[2] = 0.5 * y[m] * (y[m] - 1.0);
#if AMREX_SPACEDIM > 2
        cz[0] = 0.5 * (z[m] - 1.0) * (z[m] - 2.0);
        cz[1] = -z[m] * (z[m] - 2.0);
        cz[2] = 0.5 * z[m] * (z[m] - 1.0);
#endif

#if AMREX_SPACEDIM > 2
        for (int kk = 0; kk < 3; kk++) {
#endif
          for (int jj = 0; jj < 3; jj++) {
            sten AMREX_D_TERM([m + 1], [jj], [kk]) =
              AMREX_D_TERM(1.0, *cy[jj], *cz[kk]);
          }
#if AMREX_SPACEDIM > 2
        }
#endif
      }

      for (int jj = 0; jj < 3; jj++) {
#if AMREX_SPACEDIM > 2
        for (int kk = 0; kk < 3; kk++) {
#endif
          sten AMREX_D_TERM([1], [jj], [kk]) *= d[1] / (d[0] * (d[1] - d[0]));
          sten AMREX_D_TERM([2], [jj], [kk]) *= d[0] / (d[1] * (d[0] - d[1]));
#if AMREX_SPACEDIM > 2
        }
#endif
      }
      const amrex::Real bcs = -(d[0] + d[1]) / (d[0] * d[1]);

      // Transform stencil into regular stencil structure
      amrex::Real tsten AMREX_D_TERM([3],[3],[3]) = AMREX_D_TERM({,{,{) 0.0 AMREX_D_TERM(},},});
      int ivl[AMREX_SPACEDIM] = {0};
      for (int ii = 0; ii < 3; ii++) {
        for (int jj = 0; jj < 3; jj++) {
#if AMREX_SPACEDIM > 2
          for (int kk = 0; kk < 3; kk++) {
#endif
            AMREX_D_TERM(ivl[c[0]] = ii * s[0] + ivs[c[0]] - baseiv[c[0]];
                         , ivl[c[1]] = jj * s[1] + ivs[c[1]] - baseiv[c[1]];
                         , ivl[c[2]] = kk * s[2] + ivs[c[2]] - baseiv[c[2]];)
            tsten AMREX_D_TERM([ivl[0]], [ivl[1]], [ivl[2]]) =
              sten AMREX_D_TERM([ii], [jj], [kk]);
#if AMREX_SPACEDIM > 2
          }
#endif
        }
      }

      for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
        grad_stencil[L].iv[dir] = ebg[L].iv[dir];
        // Shift base down, if required;
        grad_stencil[L].iv_base[dir] = baseiv[dir] + sh[dir];
      }

      for (int ii = 0; ii < 3; ii++) {
        for (int jj = 0; jj < 3; jj++) {
#if AMREX_SPACEDIM > 2
          for (int kk = 0; kk < 3; kk++) {
#endif
            grad_stencil[L].val PELEC_D_TERM_REVERSE([kk], [jj], [ii]) =
              fac * ebg[L].eb_area * tsten AMREX_D_TERM([ii], [jj], [kk]);
#if AMREX_SPACEDIM > 2
          }
#endif
        }
      }
      grad_stencil[L].bcval_sten = fac * ebg[L].eb_area * bcs;
    }
  });
}

// This least-squares procedure is adapted from what is
// done in unstructured grid codes
// see pages 162-164 in chapter
//"Spatial discretization: unstructured finite volume schemes"
// Blazek's text book
//(Computational Fluid Dynamics: Principles and Applications)
// The QR decomposition procedure is also described in
// Anderson, W. K., & Bonhaus, D. L. (1994). An implicit upwind algorithm for
// computing turbulent flows on unstructured grids.
// Computers & Fluids, 23(1), 1-21.
//

// Anderson and Bonhaus also go on to say in their paper that

//"Numerical experiments have been conducted which indicate that for
// reconstructing nonlinear data on highly stretched meshes using equation 18,
// the unweighted formulation is far superior to either inverse distance
// weighting or the use of gradients calculated with Green’s theorem. It was
// found, however, that for computing the actual values of gradients, inverse
// distance weighting and Greens theorem give very similar results, both of
// which are more accurate than unweighted least squares. Their failure to
// accurately reproduce surrounding data via equation 18 is attributable to the
// flawed assumption of linearly varying data. Therefore, for reconstructing
// data on boundaries of control volumes, the unweighted least squares procedure
// is used. When actual gradients are required, as in the production terms for
// the turbulence models, Green’s theorem"

//"I guess the upshot from Anderson and Bonhaus's statement is that
// least-squares is robust and more accurate on stretched weird grids but may
// not be necassarily more accurate than Green-Gauss or even quadratic Johansen
// on nicer grids": Hari S

void
pc_fill_bndry_grad_stencil_ls(
  const amrex::Box& bx,
  const amrex::Real dx,
  const int /*Nebg*/,
  const EBBndryGeom* ebg,
  const int Nsten,
  const amrex::Array4<amrex::EBCellFlag const>& flags,
  EBBndrySten* grad_stencil)
{

  AMREX_ASSERT(AMREX_SPACEDIM > 1);

  const auto lo = amrex::lbound(bx);
  const auto hi = amrex::ubound(bx);
  const amrex::Real area = std::pow(dx, AMREX_SPACEDIM - 1);
  const amrex::Real fac = area / dx;

  amrex::ParallelFor(Nsten, [=] AMREX_GPU_DEVICE(int L) {
    if (is_inside(ebg[L].iv, lo, hi)) {
      const amrex::Real n[AMREX_SPACEDIM] = {AMREX_D_DECL(
        ebg[L].eb_normal[0], ebg[L].eb_normal[1], ebg[L].eb_normal[2])};

      amrex::Real centcoord[AMREX_SPACEDIM] = {AMREX_D_DECL(
        ebg[L].eb_centroid[0], ebg[L].eb_centroid[1], ebg[L].eb_centroid[2])};

      const int ivs[AMREX_SPACEDIM] = {
        AMREX_D_DECL(ebg[L].iv[0], ebg[L].iv[1], ebg[L].iv[2])};

      // From ivs, move to center of stencil, then move to lower-left of that
      const int baseiv[AMREX_SPACEDIM] = {AMREX_D_DECL(
        ivs[0] + (int)amrex::Math::copysign(1.0, n[0]) - 1,
        ivs[1] + (int)amrex::Math::copysign(1.0, n[1]) - 1,
        ivs[2] + (int)amrex::Math::copysign(1.0, n[2]) - 1)};

      // set iv and iv_base
      for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
        grad_stencil[L].iv[dir] = ebg[L].iv[dir];
        grad_stencil[L].iv_base[dir] = baseiv[dir];
      }

      // do not include self cell
      amrex::Real AMREX_D_DECL(
        delta_x_i[NLSPTS], delta_y_i[NLSPTS], delta_z_i[NLSPTS]);
      amrex::Real qmat[NEL_TRIMAT], wvec[NLSPTS][AMREX_SPACEDIM];
      amrex::Real xi[AMREX_SPACEDIM];

      int iter = 0;
      amrex::IntVect sten_iv;

      for (int ii = 0; ii < 3; ii++) {
        for (int jj = 0; jj < 3; jj++) {
#if AMREX_SPACEDIM > 2
          for (int kk = 0; kk < 3; kk++) {
#endif
            AMREX_D_TERM(sten_iv[0] = baseiv[0] + ii;
                         , sten_iv[1] = baseiv[1] + jj;
                         , sten_iv[2] = baseiv[2] + kk;)

            if (is_inside(sten_iv, lo, hi)) {

              if (
                !(AMREX_D_TERM(
                  sten_iv[0] == ivs[0], &&sten_iv[1] == ivs[1],
                  &&sten_iv[2] == ivs[2]))

                && !flags(sten_iv).isCovered()) {
                AMREX_D_TERM(xi[0] = sten_iv[0] - ivs[0];
                             , xi[1] = sten_iv[1] - ivs[1];
                             , xi[2] = sten_iv[2] - ivs[2];)

                AMREX_D_TERM(delta_x_i[iter] = (xi[0] - centcoord[0]);
                             , delta_y_i[iter] = (xi[1] - centcoord[1]);
                             , delta_z_i[iter] = (xi[2] - centcoord[2]);)

                iter++;
              }
            }

#if AMREX_SPACEDIM > 2
          }
#endif
        }
      }

      if (iter > AMREX_SPACEDIM) {

        get_qmat(AMREX_D_DECL(delta_x_i, delta_y_i, delta_z_i), iter, qmat);
        get_weightvec(
          AMREX_D_DECL(delta_x_i, delta_y_i, delta_z_i), iter, qmat, wvec);

        iter = 0;
        amrex::Real selfweight = 0.0;
        for (int ii = 0; ii < 3; ii++) {
          for (int jj = 0; jj < 3; jj++) {
#if AMREX_SPACEDIM > 2
            for (int kk = 0; kk < 3; kk++) {
#endif
              AMREX_D_TERM(sten_iv[0] = baseiv[0] + ii;
                           , sten_iv[1] = baseiv[1] + jj;
                           , sten_iv[2] = baseiv[2] + kk;)

              if (is_inside(sten_iv, lo, hi)) {

                if (
                  !(AMREX_D_TERM(
                    sten_iv[0] == ivs[0], &&sten_iv[1] == ivs[1],
                    &&sten_iv[2] == ivs[2]))

                  && !flags(sten_iv).isCovered()) {
                  grad_stencil[L].val PELEC_D_TERM_REVERSE([kk], [jj], [ii]) =
                    fac * ebg[L].eb_area *
                    (AMREX_D_TERM(
                      wvec[iter][0] * n[0], +wvec[iter][1] * n[1],
                      +wvec[iter][2] * n[2]));
                  selfweight +=
                    grad_stencil[L].val PELEC_D_TERM_REVERSE([kk], [jj], [ii]);
                  iter++;
                } else {
                  grad_stencil[L].val PELEC_D_TERM_REVERSE([kk], [jj], [ii]) =
                    0.0;
                }
              } else {
                grad_stencil[L].val PELEC_D_TERM_REVERSE([kk], [jj], [ii]) =
                  0.0;
              }
#if AMREX_SPACEDIM > 2
            }
#endif
          }
        }
        grad_stencil[L].bcval_sten = -selfweight;
      } else {
        grad_stencil[L].bcval_sten = 0.0;

        for (int ii = 0; ii < 3; ii++) {

          for (int jj = 0; jj < 3; jj++) {

#if AMREX_SPACEDIM > 2
            for (int kk = 0; kk < 3; kk++) { //NOLINT
#endif
              grad_stencil[L].val PELEC_D_TERM_REVERSE([kk], [jj], [ii]) = 0.0;

#if AMREX_SPACEDIM > 2
            }
#endif
          }
        }
      }
    }
  });
}

void
pc_fill_flux_interp_stencil(
  const amrex::Box& bx,
  const amrex::Box /*fbx*/,
  const int Nsten,
  const amrex::Array4<const amrex::Real>& fc,
  const amrex::Array4<const amrex::Real>& fa,
  FaceSten* sten)
{
  const auto lo = amrex::lbound(bx);
  const auto hi = amrex::ubound(bx);
  amrex::ParallelFor(Nsten, [=] AMREX_GPU_DEVICE(int L) {
    const auto& iv = sten[L].iv;
    if (is_inside(iv, lo, hi)) {
#if AMREX_SPACEDIM == 2
      for (amrex::Real& jj : sten[L].val) {
        jj = 0.0;
      }
      const amrex::Real ct = fc(iv, 0);
      const int tn = (int)amrex::Math::copysign(1.0, ct);
      const amrex::Real act = amrex::Math::abs(ct);
      sten[L].val[1] = fa(iv) * (1.0 - act);
      sten[L].val[tn + 1] = fa(iv) * act;
#elif AMREX_SPACEDIM == 3
      for (auto& ii : sten[L].val) {
        for (amrex::Real& jj : ii) {
          jj = 0.0;
        }
      }
      const amrex::Real ct0 = fc(iv, 0);
      const amrex::Real ct1 = fc(iv, 1);
      const int t0n = (int)amrex::Math::copysign(1.0, ct0);
      const int t1n = (int)amrex::Math::copysign(1.0, ct1);
      const amrex::Real act0 = amrex::Math::abs(ct0);
      const amrex::Real act1 = amrex::Math::abs(ct1);
      sten[L].val[1][1] = fa(iv) * (1.0 - act0) * (1.0 - act1);
      sten[L].val[1][t0n + 1] = fa(iv) * act0 * (1.0 - act1);
      sten[L].val[t1n + 1][1] = fa(iv) * (1.0 - act0) * act1;
      sten[L].val[t1n + 1][t0n + 1] = fa(iv) * act0 * act1;
#endif
    }
  });
}

void
pc_apply_face_stencil(
  const amrex::Box& bx,
  const amrex::Box /*sbx*/,
  const FaceSten* sten,
  const int Nsten,
  const int dir,
  const int nc,
  const amrex::Array4<amrex::Real>& vout)
{
  const auto lo = amrex::lbound(bx);
  const auto hi = amrex::ubound(bx);

  for (int n = 0; n < nc; n++) {
    amrex::Gpu::DeviceVector<amrex::Real> newval(Nsten);
    amrex::Real* d_newval = newval.data();

    amrex::ParallelFor(Nsten, [=] AMREX_GPU_DEVICE(int L) {
      const auto& iv = sten[L].iv;
      d_newval[L] = 0.0;
      if (is_inside(iv, lo, hi)) {
        if (dir == 0) {
          for (int t0 = 0; t0 < 3; t0++) {
#if AMREX_SPACEDIM > 2
            for (int t1 = 0; t1 < 3; t1++) {
#endif
              const amrex::IntVect ivd = amrex::IntVect{
                AMREX_D_DECL(iv[0], iv[1] - 1 + t0, iv[2] - 1 + t1)};
              d_newval[L] +=
                sten[L].val PELEC_D_TERM_REVERSE([t1], [t0], ) * vout(ivd, n);
#if AMREX_SPACEDIM > 2
            }
#endif
          }
        } else if (dir == 1) {
          for (int t0 = 0; t0 < 3; t0++) {
#if AMREX_SPACEDIM > 2
            for (int t1 = 0; t1 < 3; t1++) {
#endif
              const amrex::IntVect ivd = amrex::IntVect{
                AMREX_D_DECL(iv[0] - 1 + t0, iv[1], iv[2] - 1 + t1)};
              d_newval[L] +=
                sten[L].val PELEC_D_TERM_REVERSE([t1], [t0], ) * vout(ivd, n);
#if AMREX_SPACEDIM > 2
            }
#endif
          }
        } else if (dir == 2) {
#if AMREX_SPACEDIM > 2
          for (int t0 = 0; t0 < 3; t0++) {
            for (int t1 = 0; t1 < 3; t1++) {
              const amrex::IntVect ivd = amrex::IntVect{
                AMREX_D_DECL(iv[0] - 1 + t0, iv[1] - 1 + t1, iv[2])};
              d_newval[L] +=
                sten[L].val PELEC_D_TERM_REVERSE([t1], [t0], ) * vout(ivd, n);
            }
          }
#endif
        }
      }
    });

    amrex::ParallelFor(Nsten, [=] AMREX_GPU_DEVICE(int L) {
      if (is_inside(sten[L].iv, lo, hi)) {
        vout(sten[L].iv, n) = d_newval[L];
      }
    });
  }
}

void
pc_eb_div(
  const amrex::Box& bx,
  const amrex::Real vol,
  const int nc,
  const EBBndryGeom* sv_ebg,
  const int Ncut,
  AMREX_D_DECL(
    const amrex::Array4<const amrex::Real>& f0,
    const amrex::Array4<const amrex::Real>& f1,
    const amrex::Array4<const amrex::Real>& f2),
  const amrex::Real* ebflux,
  const amrex::Array4<const amrex::Real>& vf,
  const amrex::Array4<amrex::Real>& DC)
{
  const auto lo = amrex::lbound(bx);
  const auto hi = amrex::ubound(bx);
  const amrex::Real volinv = 1.0 / vol;

  for (int n = 0; n < nc; n++) {
    // Recompute conservative divergence, DC, on cut cells...need DC in 2 grow
    // cells for final result
    amrex::ParallelFor(Ncut, [=] AMREX_GPU_DEVICE(int L) {
      const auto& iv = sv_ebg[L].iv;
      if (is_inside(iv, lo, hi, 2)) {
        const amrex::Real kappa_inv =
          1.0 / amrex::max<amrex::Real>(vf(iv), 1.0e-12);
        amrex::Real tmp;
#ifdef _OPENMP
#pragma omp atomic read
#endif
        tmp = ebflux[n * Ncut + L];
        DC(iv, n) =
          -(AMREX_D_TERM(
              f0(iv + amrex::IntVect::TheDimensionVector(0), n) - f0(iv, n),
              +f1(iv + amrex::IntVect::TheDimensionVector(1), n) - f1(iv, n),
              +f2(iv + amrex::IntVect::TheDimensionVector(2), n) - f2(iv, n)) +
            tmp) *
          volinv * kappa_inv;
      }
    });
  }
}

void
pc_apply_eb_boundry_visc_flux_stencil(
  const amrex::Box& bx,
  const EBBndrySten* sten,
  const int Nsten,
  const EBBndryGeom* ebg,
  const int /*Nebg*/,
  amrex::Array4<const amrex::Real> const& q,
  amrex::Array4<const amrex::Real> const& coeff,
  const amrex::Real* bcval,
  const int /*Nvals*/,
  amrex::Real* bcflux,
  const int Nflux)
{
  const auto lo = amrex::lbound(bx);
  const auto hi = amrex::ubound(bx);

  amrex::ParallelFor(Nsten, [=] AMREX_GPU_DEVICE(int L) {
    const auto& iv = sten[L].iv;
    if (is_inside(iv, lo, hi)) {
      const amrex::Real Nmag = std::sqrt(AMREX_D_TERM(
        ebg[L].eb_normal[0] * ebg[L].eb_normal[0],
        +ebg[L].eb_normal[1] * ebg[L].eb_normal[1],
        +ebg[L].eb_normal[2] * ebg[L].eb_normal[2]));
      const amrex::Real norm[AMREX_SPACEDIM] = {AMREX_D_DECL(
        ebg[L].eb_normal[0] / Nmag, ebg[L].eb_normal[1] / Nmag,
        ebg[L].eb_normal[2] / Nmag)};

#if AMREX_SPACEDIM == 2
      const amrex::Real t1[AMREX_SPACEDIM] = {-norm[1], norm[0]};
#elif AMREX_SPACEDIM == 3
      amrex::Real alpha[AMREX_SPACEDIM] = {0.0};
      int c[AMREX_SPACEDIM] = {0};
      idxsort(norm, c);
      alpha[c[AMREX_D_PICK(0, 1, 2)]] = 1.0;
      const amrex::Real ndota = AMREX_D_TERM(
        norm[0] * alpha[0], +norm[1] * alpha[1], +norm[2] * alpha[2]);
      amrex::Real t1[AMREX_SPACEDIM];
      for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
        t1[idir] = alpha[idir] - ndota * norm[idir];
      }

      const amrex::Real denom =
        1.0 /
        std::sqrt(AMREX_D_TERM(t1[0] * t1[0], +t1[1] * t1[1], +t1[2] * t1[2]));
      for (amrex::Real& idir : t1) {
        idir *= denom;
      }

      const amrex::Real t2[AMREX_SPACEDIM] = {AMREX_D_DECL(
        norm[1] * t1[2] - norm[2] * t1[1], norm[2] * t1[0] - norm[0] * t1[2],
        norm[0] * t1[1] - norm[1] * t1[0])};
#endif

      amrex::Real Qt[AMREX_SPACEDIM][AMREX_SPACEDIM];
      for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
        AMREX_D_TERM(Qt[0][idir] = norm[idir];, Qt[1][idir] = t1[idir];
                     , Qt[2][idir] = t2[idir];)
      }

      // Transform velocities at stencil points to coordinates aligned with EB
      amrex::Real Uo AMREX_D_TERM([3], [3], [3])[AMREX_SPACEDIM];
      for (int ii = 0; ii < 3; ii++) {
        for (int jj = 0; jj < 3; jj++) {
#if AMREX_SPACEDIM > 2
          for (int kk = 0; kk < 3; kk++) {
#endif
            for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
              const amrex::IntVect ivp =
                sten[L].iv_base + amrex::IntVect{AMREX_D_DECL(ii, jj, kk)};
              Uo AMREX_D_TERM([ii], [jj], [kk])[idir] = q(ivp, QU + idir);
            }
#if AMREX_SPACEDIM > 2
          }
#endif
        }
      }
      amrex::Real Ut AMREX_D_TERM([3], [3], [3])[AMREX_SPACEDIM];
      for (int ii = 0; ii < 3; ii++) {
        for (int jj = 0; jj < 3; jj++) {
#if AMREX_SPACEDIM > 2
          for (int kk = 0; kk < 3; kk++) {
#endif
            for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
              Ut AMREX_D_TERM([ii], [jj], [kk])[idir] = AMREX_D_TERM(
                Qt[idir][0] * Uo AMREX_D_TERM([ii], [jj], [kk])[0],
                +Qt[idir][1] * Uo AMREX_D_TERM([ii], [jj], [kk])[1],
                +Qt[idir][2] * Uo AMREX_D_TERM([ii], [jj], [kk])[2]);
            }
#if AMREX_SPACEDIM > 2
          }
#endif
        }
      }

      // Transform eb boundary velocities to coordinates aligned with EB
      amrex::Real bco[AMREX_SPACEDIM];
      for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
        bco[idir] = bcval[idir * Nsten + L];
      }

      amrex::Real bct[AMREX_SPACEDIM];
      for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
        bct[idir] = AMREX_D_TERM(
          Qt[idir][0] * bco[0], +Qt[idir][1] * bco[1], +Qt[idir][2] * bco[2]);
      }

      // Compute normal derivative (times eb area) using precomputed stencil
      amrex::Real sum[AMREX_SPACEDIM] = {0.0};
      for (int ii = 0; ii < 3; ii++) {
        for (int jj = 0; jj < 3; jj++) {
#if AMREX_SPACEDIM > 2
          for (int kk = 0; kk < 3; kk++) {
#endif
            for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
              sum[idir] += sten[L].val PELEC_D_TERM_REVERSE([kk], [jj], [ii]) *
                           Ut AMREX_D_TERM([ii], [jj], [kk])[idir];
            }
#if AMREX_SPACEDIM > 2
          }
#endif
        }
      }
      amrex::Real dUtdn[AMREX_SPACEDIM];
      for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
        dUtdn[idir] = sum[idir] + bct[idir] * sten[L].bcval_sten;
      }

      const amrex::Real tauDotN[AMREX_SPACEDIM] = {AMREX_D_DECL(
        ((4.0 / 3.0) * coeff(iv, dComp_mu) + coeff(iv, dComp_xi)) * dUtdn[0],
        coeff(iv, dComp_mu) * dUtdn[1], coeff(iv, dComp_mu) * dUtdn[2])};

      for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
        bcflux[idir * Nflux + L] = AMREX_D_TERM(
          Qt[0][idir] * tauDotN[0], +Qt[1][idir] * tauDotN[1],
          +Qt[2][idir] * tauDotN[2]);
      }
    }
  });
}

void
pc_apply_eb_boundry_flux_stencil(
  const amrex::Box& bx,
  const EBBndrySten* sten,
  const int Nsten,
  amrex::Array4<const amrex::Real> const& s,
  const int scomp,
  amrex::Array4<const amrex::Real> const& D,
  const int Dcomp,
  const amrex::Real* bcval,
  const int /*Nvals*/,
  amrex::Real* bcflux,
  const int Nflux,
  const int nc)
{
  const auto lo = amrex::lbound(bx);
  const auto hi = amrex::ubound(bx);

  amrex::ParallelFor(Nsten, [=] AMREX_GPU_DEVICE(int L) {
    const amrex::IntVect iv = sten[L].iv;
    if (is_inside(iv, lo, hi)) {
      for (int n = 0; n < nc; n++) {
        amrex::Real sum = 0.0;
        for (int ii = 0; ii < 3; ii++) {
          for (int jj = 0; jj < 3; jj++) {
#if AMREX_SPACEDIM > 2
            for (int kk = 0; kk < 3; kk++) {
#endif
              const amrex::IntVect ivp = amrex::IntVect(AMREX_D_DECL(
                sten[L].iv_base[0] + ii, sten[L].iv_base[1] + jj,
                sten[L].iv_base[2] + kk));
              sum += sten[L].val PELEC_D_TERM_REVERSE([kk], [jj], [ii]) *
                     s(ivp, scomp + n);
#if AMREX_SPACEDIM > 2
            }
#endif
          }
        }
        bcflux[n * Nflux + L] =
          D(iv, Dcomp + n) * (bcval[n * Nsten + L] * sten[L].bcval_sten + sum);
      }
    }
  });
}
