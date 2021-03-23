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
pc_fix_div_and_redistribute(
  const amrex::Box& bx,
  const amrex::Real vol,
  const amrex::Real dt,
  const int nc,
  const amrex::Real eb_small_vfrac,
  const bool levmsk_notcovered,
  const EBBndryGeom* sv_ebg,
  const int Ncut,
  const amrex::Array4<amrex::EBCellFlag const>& flags,
  const amrex::Array4<const amrex::Real>& f0,
  const amrex::Array4<const amrex::Real>& f1,
  const amrex::Array4<const amrex::Real>& f2,
  const amrex::Real* ebflux,
  const int /*nebflux*/,
  const amrex::Array4<const amrex::Real>& vf,
  const amrex::Array4<const amrex::Real>& W,
  const bool as_crse,
  const bool as_fine,
  const amrex::Array4<const int>& levmsk,
  const amrex::Array4<const int>& rr_flag_crse,
  const amrex::Array4<amrex::Real>& DC,
  const amrex::Array4<amrex::Real>& rr_drho_crse,
  const amrex::Array4<amrex::Real>& dm_as_fine)
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

    // Compute non-conservative and hybrid divergence, DNC and HD, and
    // redistribution mass dM in cut cells. Will need in 1 grow cells (see
    // below), so it depends on having a conservative div in 2 grow cells
    amrex::Gpu::DeviceVector<amrex::Real> dM_vec(Ncut);
    amrex::Gpu::DeviceVector<amrex::Real> HD_vec(Ncut);
    amrex::Real* dM = dM_vec.data();
    amrex::Real* HD = HD_vec.data();
    amrex::Gpu::synchronize();
    amrex::ParallelFor(Ncut, [=] AMREX_GPU_DEVICE(int L) {
      const int i = sv_ebg[L].iv[0];
      const int j = sv_ebg[L].iv[1];
      const int k = sv_ebg[L].iv[2];
      if (is_inside(i, j, k, lo, hi, 1)) {
        amrex::Real sum_kappa = 0.0;
        amrex::Real sum_div = 0.0;
        for (int ii = -1; ii <= 1; ii++) {
          for (int jj = -1; jj <= 1; jj++) {
            for (int kk = -1; kk <= 1; kk++) {
              sum_kappa += flags(i, j, k).isConnected(ii, jj, kk) *
                           vf(i + ii, j + jj, k + kk);
              sum_div += flags(i, j, k).isConnected(ii, jj, kk) *
                         vf(i + ii, j + jj, k + kk) *
                         DC(i + ii, j + jj, k + kk, n);
            }
          }
        }
        const amrex::Real DNC = sum_div / sum_kappa;
        if (sv_ebg[L].eb_vfrac < eb_small_vfrac) {
          dM[L] = vf(i, j, k) * DC(i, j, k, n);
          HD[L] = 0.0;
        } else {
          dM[L] = vf(i, j, k) * (1.0 - vf(i, j, k)) * (DC(i, j, k, n) - DNC);
          HD[L] = vf(i, j, k) * DC(i, j, k, n) + (1.0 - vf(i, j, k)) * DNC;
        }
      }
    });

    // Now that we finished computing HD and dM everywhere, it is safe to
    // increment DC to hold HD
    amrex::ParallelFor(Ncut, [=] AMREX_GPU_DEVICE(int L) {
      const int i = sv_ebg[L].iv[0];
      const int j = sv_ebg[L].iv[1];
      const int k = sv_ebg[L].iv[2];
      if (is_inside(i, j, k, lo, hi, 1)) {
        DC(i, j, k, n) = HD[L];
      }
    });

    // Redistribute dM - THIS REQUIRES THAT DC BE GOOD IN 1 GROW CELL
    const amrex::Real reredistribution_threshold =
      amrex_eb_get_reredistribution_threshold();
    amrex::ParallelFor(Ncut, [=] AMREX_GPU_DEVICE(int L) {
      const int i = sv_ebg[L].iv[0];
      const int j = sv_ebg[L].iv[1];
      const int k = sv_ebg[L].iv[2];
      if (is_inside(i, j, k, lo, hi, 1)) {
        amrex::Real sum_kappa = 0.0;
        for (int ii = -1; ii <= 1; ii++) {
          for (int jj = -1; jj <= 1; jj++) {
            for (int kk = -1; kk <= 1; kk++) {
              int nbr = flags(i, j, k).isConnected(ii, jj, kk);
              if ((ii == 0) && (jj == 0) && (kk == 0)) {
                nbr = 0;
              }
              if (vf(i + ii, j + jj, k + kk) < eb_small_vfrac) {
                nbr = 0;
              }
              sum_kappa +=
                nbr * vf(i + ii, j + jj, k + kk) * W(i + ii, j + jj, k + kk);
            }
          }
        }
        const amrex::Real sum_kappa_inv = 1.0 / sum_kappa;
        for (int ii = -1; ii <= 1; ii++) {
          for (int jj = -1; jj <= 1; jj++) {
            for (int kk = -1; kk <= 1; kk++) {
              int nbr = flags(i, j, k).isConnected(ii, jj, kk);
              if ((ii == 0) && (jj == 0) && (kk == 0)) {
                nbr = 0;
              }
              if (vf(i + ii, j + jj, k + kk) < eb_small_vfrac) {
                nbr = 0;
              }
              amrex::Gpu::Atomic::Add(
                &DC(i + ii, j + jj, k + kk, n),
                dM[L] * nbr * W(i + ii, j + jj, k + kk) * sum_kappa_inv);
              // DC(i + ii, j + jj, k + kk, n) += dM[L] * nbr * W(i + ii, j +
              // jj, k + kk) * sum_kappa_inv;
            }
          }
        }

        // re redistribution book keeping
        bool as_crse_crse_cell = false;
        bool as_crse_covered_cell = false;
        if (as_crse) {
          as_crse_crse_cell =
            is_inside(i, j, k, lo, hi) &&
            (rr_flag_crse(i, j, k) == amrex_yafluxreg_crse_fine_boundary_cell);
          as_crse_covered_cell =
            rr_flag_crse(i, j, k) == amrex_yafluxreg_fine_cell;
        }

        bool as_fine_valid_cell = false; // valid cells near box boundary
        bool as_fine_ghost_cell =
          false; // ghost cells just outside valid region
        if (as_fine) {
          as_fine_valid_cell = is_inside(i, j, k, lo, hi);
          as_fine_ghost_cell =
            (levmsk(i, j, k) ==
             levmsk_notcovered); // not covered by other grids
        }

        for (int ii = -1; ii <= 1; ii++) {
          for (int jj = -1; jj <= 1; jj++) {
            for (int kk = -1; kk <= 1; kk++) {
              if (
                ((ii != 0) || (jj != 0) || (kk != 0)) &&
                flags(i, j, k).isConnected(ii, jj, kk)) {

                const int iii = i + ii;
                const int jjj = j + jj;
                const int kkk = k + kk;

                const amrex::Real drho =
                  dM[L] * sum_kappa_inv * W(iii, jjj, kkk);
                const bool valid_dst_cell = is_inside(iii, jjj, kkk, lo, hi);

                if (
                  (as_crse_crse_cell) &&
                  (rr_flag_crse(iii, jjj, kkk) == amrex_yafluxreg_fine_cell) &&
                  (vf(i, j, k) > reredistribution_threshold)) {
                  rr_drho_crse(i, j, k, n) +=
                    dt * drho * (vf(iii, jjj, kkk) / vf(i, j, k));
                }

                if (
                  (as_crse_covered_cell) && (valid_dst_cell) &&
                  (rr_flag_crse(iii, jjj, kkk) ==
                   amrex_yafluxreg_crse_fine_boundary_cell) &&
                  (vf(iii, jjj, kkk) > reredistribution_threshold)) {
                  // the recipient is a crse/fine boundary cell
                  rr_drho_crse(iii, jjj, kkk, n) -= dt * drho;
                }

                if ((as_fine_valid_cell) && (!valid_dst_cell)) {
                  dm_as_fine(iii, jjj, kkk, n) += dt * drho * vf(iii, jjj, kkk);
                }

                if ((as_fine_ghost_cell) && (valid_dst_cell)) {
                  dm_as_fine(i, j, k, n) -= dt * drho * vf(iii, jjj, kkk);
                }
              }
            }
          }
        }
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
