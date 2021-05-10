#include "EB.H"

void
pc_fill_sv_ebg(
  const amrex::Box& bx,
  const int Nebg,
  const amrex::Array4<const amrex::Real>& vfrac,
  const amrex::Array4<const amrex::Real>& bcent,
  const amrex::Array4<const amrex::Real>& apx,
  const amrex::Array4<const amrex::Real>& apy,
  const amrex::Array4<const amrex::Real>& apz,
  EBBndryGeom* ebg)
{
  const auto lo = amrex::lbound(bx);
  const auto hi = amrex::ubound(bx);
  amrex::ParallelFor(Nebg, [=] AMREX_GPU_DEVICE(int L) {
    const int i = ebg[L].iv[0];
    const int j = ebg[L].iv[1];
    const int k = ebg[L].iv[2];
    if (is_inside(i, j, k, lo, hi)) {
      const amrex::Real axm = apx(i, j, k);
      const amrex::Real axp = apx(i + 1, j, k);
      const amrex::Real aym = apy(i, j, k);
      const amrex::Real ayp = apy(i, j + 1, k);
      const amrex::Real azm = apz(i, j, k);
      const amrex::Real azp = apz(i, j, k + 1);
      const amrex::Real apnorm = std::sqrt(
        (axm - axp) * (axm - axp) + (aym - ayp) * (aym - ayp) +
        (azm - azp) * (azm - azp));
      if (apnorm == 0.0) {
        amrex::Abort("pc_fill_sv_ebg: zero apnorm");
      }
      const amrex::Real apnorminv = -1.0 / apnorm;
      ebg[L].eb_normal[0] = (axm - axp) * apnorminv; // pointing to the wall
      ebg[L].eb_normal[1] = (aym - ayp) * apnorminv;
      ebg[L].eb_normal[2] = (azm - azp) * apnorminv;
      ebg[L].eb_area = apnorm;
      ebg[L].eb_centroid[0] = bcent(i, j, k, 0);
      ebg[L].eb_centroid[1] = bcent(i, j, k, 1);
      ebg[L].eb_centroid[2] = bcent(i, j, k, 2);
      ebg[L].eb_vfrac = vfrac(i, j, k);
    }
  });
}

void
pc_fill_bndry_grad_stencil(
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
    const int i = ebg[L].iv[0];
    const int j = ebg[L].iv[1];
    const int k = ebg[L].iv[2];
    if (is_inside(i, j, k, lo, hi)) {
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
      amrex::Real z[2] = {
        b[2] + (x[0] - b[0]) * amrex::Math::abs(n[c[2]] / n[c[0]]),
        b[2] + (x[1] - b[0]) * amrex::Math::abs(n[c[2]] / n[c[0]])};

      int sh[AMREX_SPACEDIM] = {0};
      if (y[0] < 0.0 || y[1] < 0.0) {
        sh[c[1]] = -s[1]; // Slide stencil down to avoid extrapolating, push
                          // up eb, shift down base later
        b[1] += 1;
        y[0] = b[1] + (x[0] - b[0]) * amrex::Math::abs(n[c[1]] / n[c[0]]);
        y[1] = b[1] + (x[1] - b[0]) * amrex::Math::abs(n[c[1]] / n[c[0]]);
      }
      if (z[0] < 0.0 || z[1] < 0.0) {
        sh[c[2]] = -s[2]; // Slide stencil down to avoid extrapolating, push
                          // up eb, shift down base later
        b[2] += 1;
        z[0] = b[2] + (x[0] - b[0]) * amrex::Math::abs(n[c[2]] / n[c[0]]);
        z[1] = b[2] + (x[1] - b[0]) * amrex::Math::abs(n[c[2]] / n[c[0]]);
      }
      const amrex::Real d[2] = {
        std::sqrt(
          (x[0] - b[0]) * (x[0] - b[0]) + (y[0] - b[1]) * (y[0] - b[1]) +
          (z[0] - b[2]) * (z[0] - b[2])),
        std::sqrt(
          (x[1] - b[0]) * (x[1] - b[0]) + (y[1] - b[1]) * (y[1] - b[1]) +
          (z[1] - b[2]) * (z[1] - b[2]))};

      amrex::Real sten[3][3][3] = {{{0.0}}};
      // The two intersections, that are d(1) and d(2) away from the eb
      // centroid, are both in y-z planes, bounded in (0:2)x(0:2) in normalized
      // coordinates For point m, we interpolate z=0,1,2 lines, to (y(m),0),
      // (y(m),1) and (y(m),2), and then interpolate along y=y(m) to (y(m),z(m))
      amrex::Real cy[3];
      amrex::Real cz[3];
      for (int m = 0; m < 2; m++) {
        cy[0] = 0.5 * (y[m] - 1.0) * (y[m] - 2.0);
        cy[1] = -y[m] * (y[m] - 2.0);
        cy[2] = 0.5 * y[m] * (y[m] - 1.0);
        cz[0] = 0.5 * (z[m] - 1.0) * (z[m] - 2.0);
        cz[1] = -z[m] * (z[m] - 2.0);
        cz[2] = 0.5 * z[m] * (z[m] - 1.0);

        for (int kk = 0; kk < 3; kk++) {
          for (int jj = 0; jj < 3; jj++) {
            sten[m + 1][jj][kk] = cy[jj] * cz[kk];
          }
        }
      }

      for (int jj = 0; jj < 3; jj++) {
        for (int kk = 0; kk < 3; kk++) {
          sten[1][jj][kk] *= d[1] / (d[0] * (d[1] - d[0]));
          sten[2][jj][kk] *= d[0] / (d[1] * (d[0] - d[1]));
        }
      }
      const amrex::Real bcs = -(d[0] + d[1]) / (d[0] * d[1]);

      // Transform stencil into regular stencil structure
      amrex::Real tsten[3][3][3] = {{{0.0}}};
      int iv[3] = {0};
      for (int ii = 0; ii < 3; ii++) {
        for (int jj = 0; jj < 3; jj++) {
          for (int kk = 0; kk < 3; kk++) {
            iv[c[0]] = ii * s[0] + ivs[c[0]] - baseiv[c[0]];
            iv[c[1]] = jj * s[1] + ivs[c[1]] - baseiv[c[1]];
            iv[c[2]] = kk * s[2] + ivs[c[2]] - baseiv[c[2]];
            tsten[iv[0]][iv[1]][iv[2]] = sten[ii][jj][kk];
          }
        }
      }

      for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
        grad_stencil[L].iv[dir] = ebg[L].iv[dir];
        // Shift base down, if required;
        grad_stencil[L].iv_base[dir] = baseiv[dir] + sh[dir];
        for (int jj = 0; jj < 3; jj++) {
          for (int kk = 0; kk < 3; kk++) {
            grad_stencil[L].val[kk][jj][dir] =
              fac * ebg[L].eb_area * tsten[dir][jj][kk];
          }
        }
      }
      grad_stencil[L].bcval_sten = fac * ebg[L].eb_area * bcs;
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
    const int i = sten[L].iv[0];
    const int j = sten[L].iv[1];
    const int k = sten[L].iv[2];
    if (is_inside(i, j, k, lo, hi)) {
      for (auto& ii : sten[L].val) {
        for (amrex::Real& jj : ii) {
          jj = 0.0;
        }
      }
      const amrex::Real ct0 = fc(i, j, k, 0);
      const amrex::Real ct1 = fc(i, j, k, 1);
      const int t0n = (int)amrex::Math::copysign(1.0, ct0);
      const int t1n = (int)amrex::Math::copysign(1.0, ct1);
      const amrex::Real act0 = amrex::Math::abs(ct0);
      const amrex::Real act1 = amrex::Math::abs(ct1);
      sten[L].val[1][1] = fa(i, j, k) * (1.0 - act0) * (1.0 - act1);
      sten[L].val[1][t0n + 1] = fa(i, j, k) * act0 * (1.0 - act1);
      sten[L].val[t1n + 1][1] = fa(i, j, k) * (1.0 - act0) * act1;
      sten[L].val[t1n + 1][t0n + 1] = fa(i, j, k) * act0 * act1;
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
      const int i = sten[L].iv[0];
      const int j = sten[L].iv[1];
      const int k = sten[L].iv[2];
      d_newval[L] = 0.0;
      if (is_inside(i, j, k, lo, hi)) {
        if (dir == 0) {
          for (int t0 = 0; t0 < 3; t0++) {
            for (int t1 = 0; t1 < 3; t1++) {
              d_newval[L] +=
                sten[L].val[t1][t0] * vout(i, j - 1 + t0, k - 1 + t1, n);
            }
          }
        } else if (dir == 1) {
          for (int t0 = 0; t0 < 3; t0++) {
            for (int t1 = 0; t1 < 3; t1++) {
              d_newval[L] +=
                sten[L].val[t1][t0] * vout(i - 1 + t0, j, k - 1 + t1, n);
            }
          }
        } else if (dir == 2) {
          for (int t0 = 0; t0 < 3; t0++) {
            for (int t1 = 0; t1 < 3; t1++) {
              d_newval[L] +=
                sten[L].val[t1][t0] * vout(i - 1 + t0, j - 1 + t1, k, n);
            }
          }
        }
      }
    });

    amrex::ParallelFor(Nsten, [=] AMREX_GPU_DEVICE(int L) {
      const int i = sten[L].iv[0];
      const int j = sten[L].iv[1];
      const int k = sten[L].iv[2];
      if (is_inside(i, j, k, lo, hi)) {
        vout(i, j, k, n) = d_newval[L];
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
  const amrex::Array4<const amrex::Real>& f0,
  const amrex::Array4<const amrex::Real>& f1,
  const amrex::Array4<const amrex::Real>& f2,
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
      const int i = sv_ebg[L].iv[0];
      const int j = sv_ebg[L].iv[1];
      const int k = sv_ebg[L].iv[2];
      if (is_inside(i, j, k, lo, hi, 2)) {
        const amrex::Real kappa_inv =
          1.0 / amrex::max<amrex::Real>(vf(i, j, k), 1.0e-12);
        amrex::Real tmp;
#ifdef _OPENMP
#pragma omp atomic read
#endif
        tmp = ebflux[n * Ncut + L];
        DC(i, j, k, n) =
          -(f0(i + 1, j, k, n) - f0(i, j, k, n) + f1(i, j + 1, k, n) -
            f1(i, j, k, n) + f2(i, j, k + 1, n) - f2(i, j, k, n) + tmp) *
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
    const int i = sten[L].iv[0];
    const int j = sten[L].iv[1];
    const int k = sten[L].iv[2];
    if (is_inside(i, j, k, lo, hi)) {
      const amrex::Real Nmag = std::sqrt(
        ebg[L].eb_normal[0] * ebg[L].eb_normal[0] +
        ebg[L].eb_normal[1] * ebg[L].eb_normal[1] +
        ebg[L].eb_normal[2] * ebg[L].eb_normal[2]);
      const amrex::Real norm[AMREX_SPACEDIM] = {
        ebg[L].eb_normal[0] / Nmag, ebg[L].eb_normal[1] / Nmag,
        ebg[L].eb_normal[2] / Nmag};
      amrex::Real alpha[AMREX_SPACEDIM] = {0.0};
      int c[AMREX_SPACEDIM] = {0};
      idxsort(norm, c);
      alpha[c[2]] = 1.0;
      const amrex::Real ndota =
        norm[0] * alpha[0] + norm[1] * alpha[1] + norm[2] * alpha[2];
      amrex::Real t1[AMREX_SPACEDIM];
      for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
        t1[idir] = alpha[idir] - ndota * norm[idir];
      }

      const amrex::Real denom =
        1.0 / std::sqrt(t1[0] * t1[0] + t1[1] * t1[1] + t1[2] * t1[2]);
      for (amrex::Real& idir : t1) {
        idir *= denom;
      }

      amrex::Real t2[AMREX_SPACEDIM];
      t2[0] = norm[1] * t1[2] - norm[2] * t1[1];
      t2[1] = norm[2] * t1[0] - norm[0] * t1[2];
      t2[2] = norm[0] * t1[1] - norm[1] * t1[0];

      amrex::Real Qt[AMREX_SPACEDIM][AMREX_SPACEDIM];
      Qt[0][0] = norm[0];
      Qt[0][1] = norm[1];
      Qt[0][2] = norm[2];
      Qt[1][0] = t1[0];
      Qt[1][1] = t1[1];
      Qt[1][2] = t1[2];
      Qt[2][0] = t2[0];
      Qt[2][1] = t2[1];
      Qt[2][2] = t2[2];

      // Transform velocities at stencil points to coordinates aligned with EB
      amrex::Real Uo[AMREX_SPACEDIM][AMREX_SPACEDIM][AMREX_SPACEDIM]
                    [AMREX_SPACEDIM];
      for (int ii = 0; ii < 3; ii++) {
        for (int jj = 0; jj < 3; jj++) {
          for (int kk = 0; kk < 3; kk++) {
            for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
              Uo[ii][jj][kk][idir] =
                q(sten[L].iv_base[0] + ii, sten[L].iv_base[1] + jj,
                  sten[L].iv_base[2] + kk, QU + idir);
            }
          }
        }
      }
      amrex::Real Ut[AMREX_SPACEDIM][AMREX_SPACEDIM][AMREX_SPACEDIM]
                    [AMREX_SPACEDIM];
      for (int ii = 0; ii < 3; ii++) {
        for (int jj = 0; jj < 3; jj++) {
          for (int kk = 0; kk < 3; kk++) {
            for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
              Ut[ii][jj][kk][idir] = Qt[idir][0] * Uo[ii][jj][kk][0] +
                                     Qt[idir][1] * Uo[ii][jj][kk][1] +
                                     Qt[idir][2] * Uo[ii][jj][kk][2];
            }
          }
        }
      }

      // Transform eb boundary velocities to coordinates aligned with EB
      amrex::Real bco[AMREX_SPACEDIM];
      for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
        bco[idir] = bcval[idir * Nsten + L];
      }

      amrex::Real bct[AMREX_SPACEDIM];
      for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
        bct[idir] =
          Qt[idir][0] * bco[0] + Qt[idir][1] * bco[1] + Qt[idir][2] * bco[2];
      }

      // Compute normal derivative (times eb area) using precomputed stencil
      amrex::Real sum[AMREX_SPACEDIM] = {0.0};
      for (int ii = 0; ii < 3; ii++) {
        for (int jj = 0; jj < 3; jj++) {
          for (int kk = 0; kk < 3; kk++) {
            for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
              sum[idir] += sten[L].val[kk][jj][ii] * Ut[ii][jj][kk][idir];
            }
          }
        }
      }
      amrex::Real dUtdn[AMREX_SPACEDIM];
      for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
        dUtdn[idir] = sum[idir] + bct[idir] * sten[L].bcval_sten;
      }

      amrex::Real tauDotN[AMREX_SPACEDIM];
      tauDotN[0] =
        ((4.0 / 3.0) * coeff(i, j, k, dComp_mu) + coeff(i, j, k, dComp_xi)) *
        dUtdn[0];
      tauDotN[1] = coeff(i, j, k, dComp_mu) * dUtdn[1];
      tauDotN[2] = coeff(i, j, k, dComp_mu) * dUtdn[2];

      for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
        bcflux[idir * Nflux + L] = Qt[0][idir] * tauDotN[0] +
                                   Qt[1][idir] * tauDotN[1] +
                                   Qt[2][idir] * tauDotN[2];
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
    const int i = sten[L].iv[0];
    const int j = sten[L].iv[1];
    const int k = sten[L].iv[2];
    if (is_inside(i, j, k, lo, hi)) {
      for (int n = 0; n < nc; n++) {
        amrex::Real sum = 0.0;
        for (int ii = 0; ii < 3; ii++) {
          for (int jj = 0; jj < 3; jj++) {
            for (int kk = 0; kk < 3; kk++) {
              sum += sten[L].val[kk][jj][ii] *
                     s(sten[L].iv_base[0] + ii, sten[L].iv_base[1] + jj,
                       sten[L].iv_base[2] + kk, scomp + n);
            }
          }
        }
        bcflux[n * Nflux + L] =
          D(i, j, k, Dcomp + n) *
          (bcval[n * Nsten + L] * sten[L].bcval_sten + sum);
      }
    }
  });
}
