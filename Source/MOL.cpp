#include "MOL.H"

void
pc_compute_hyp_mol_flux(
  const amrex::Box& cbox,
  const amrex::Array4<const amrex::Real>& q,
  const amrex::Array4<const amrex::Real>& qaux,
  const amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM> flx,
  const amrex::GpuArray<const amrex::Array4<const amrex::Real>, AMREX_SPACEDIM>
    area,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>
#ifdef PELEC_USE_EB
    del
#endif
  ,
  const int plm_iorder,
  const int use_laxf_flux
#ifdef PELEC_USE_EB
  ,
  const amrex::Array4<amrex::EBCellFlag const>& flags,
  const EBBndryGeom* ebg,
  const int /*Nebg*/,
  amrex::Real* ebflux,
  const int nebflux
#endif
)
{
  const int R_RHO = 0;
  const int R_UN = 1;
  const int R_UT1 = 2;
  const int R_UT2 = 3;
  const int R_P = 4;
  const int R_Y = 5;
  const int bc_test_val = 1;

  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    amrex::FArrayBox dq_fab(cbox, QVAR);
    amrex::Elixir dq_fab_eli = dq_fab.elixir();
    auto const& dq = dq_fab.array();
    setV(cbox, QVAR, dq, 0.0);

    // dimensional indexing
    const amrex::GpuArray<const int, 3> bdim{{dir == 0, dir == 1, dir == 2}};
    const amrex::GpuArray<const int, 3> q_idx{
      {bdim[0] * QU + bdim[1] * QV + bdim[2] * QW,
       bdim[0] * QV + bdim[1] * QU + bdim[2] * QU,
       bdim[0] * QW + bdim[1] * QW + bdim[2] * QV}};
    const amrex::GpuArray<const int, 3> f_idx{
      {bdim[0] * UMX + bdim[1] * UMY + bdim[2] * UMZ,
       bdim[0] * UMY + bdim[1] * UMX + bdim[2] * UMX,
       bdim[0] * UMZ + bdim[1] * UMZ + bdim[2] * UMY}};

    if (plm_iorder != 1) {
      amrex::ParallelFor(
        cbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          mol_slope(
            i, j, k, dir, q_idx, q, qaux, dq
#ifdef PELEC_USE_EB
            ,
            flags
#endif
          );
        });
    }
    const amrex::Box tbox = amrex::grow(cbox, dir, -1);
    const amrex::Box ebox = amrex::surroundingNodes(tbox, dir);
    amrex::ParallelFor(
      ebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
        const amrex::IntVect ivm(iv - amrex::IntVect::TheDimensionVector(dir));

        amrex::Real qtempl[5 + NUM_SPECIES] = {0.0};
        qtempl[R_UN] =
          q(ivm, q_idx[0]) + 0.5 * ((dq(ivm, 1) - dq(ivm, 0)) / q(ivm, QRHO));
        qtempl[R_P] =
          q(ivm, QPRES) + 0.5 * (dq(ivm, 0) + dq(ivm, 1)) * qaux(ivm, QC);
        qtempl[R_UT1] = q(ivm, q_idx[1]) + 0.5 * dq(ivm, 2);
        qtempl[R_UT2] =
          AMREX_D_PICK(0.0, 0.0, q(ivm, q_idx[2]) + 0.5 * dq(ivm, 3));
        qtempl[R_RHO] = 0.0;
        for (int n = 0; n < NUM_SPECIES; n++) {
          qtempl[R_Y + n] =
            q(ivm, QFS + n) * q(ivm, QRHO) +
            0.5 * (dq(ivm, 4 + n) +
                   q(ivm, QFS + n) * (dq(ivm, 0) + dq(ivm, 1)) / qaux(ivm, QC));
          qtempl[R_RHO] += qtempl[R_Y + n];
        }

        for (int n = 0; n < NUM_SPECIES; n++) {
          qtempl[R_Y + n] = qtempl[R_Y + n] / qtempl[R_RHO];
        }

        amrex::Real qtempr[5 + NUM_SPECIES] = {0.0};
        qtempr[R_UN] =
          q(iv, q_idx[0]) - 0.5 * ((dq(iv, 1) - dq(iv, 0)) / q(iv, QRHO));
        qtempr[R_P] =
          q(iv, QPRES) - 0.5 * (dq(iv, 0) + dq(iv, 1)) * qaux(iv, QC);
        qtempr[R_UT1] = q(iv, q_idx[1]) - 0.5 * dq(iv, 2);
        qtempr[R_UT2] =
          AMREX_D_PICK(0.0, 0.0, q(iv, q_idx[2]) - 0.5 * dq(iv, 3));
        qtempr[R_RHO] = 0.0;
        for (int n = 0; n < NUM_SPECIES; n++) {
          qtempr[R_Y + n] =
            q(iv, QFS + n) * q(iv, QRHO) -
            0.5 * (dq(iv, 4 + n) +
                   q(iv, QFS + n) * (dq(iv, 0) + dq(iv, 1)) / qaux(iv, QC));
          qtempr[R_RHO] += qtempr[R_Y + n];
        }
        for (int n = 0; n < NUM_SPECIES; n++) {
          qtempr[R_Y + n] = qtempr[R_Y + n] / qtempr[R_RHO];
        }

        const amrex::Real cavg = 0.5 * (qaux(iv, QC) + qaux(ivm, QC));

        amrex::Real spl[NUM_SPECIES];
        for (int n = 0; n < NUM_SPECIES; n++) {
          spl[n] = qtempl[R_Y + n];
        }

        amrex::Real spr[NUM_SPECIES];
        for (int n = 0; n < NUM_SPECIES; n++) {
          spr[n] = qtempr[R_Y + n];
        }

        amrex::Real flux_tmp[NVAR] = {0.0};
        amrex::Real ustar = 0.0;
        amrex::Real maxeigval = 0.0;

        if (!use_laxf_flux) {
          amrex::Real tmp0 = 0.0, tmp1 = 0.0, tmp2 = 0.0, tmp3 = 0.0,
                      tmp4 = 0.0;
          riemann(
            qtempl[R_RHO], qtempl[R_UN], qtempl[R_UT1], qtempl[R_UT2],
            qtempl[R_P], spl, qtempr[R_RHO], qtempr[R_UN], qtempr[R_UT1],
            qtempr[R_UT2], qtempr[R_P], spr, bc_test_val, cavg, ustar,
            flux_tmp[URHO], flux_tmp[f_idx[0]], flux_tmp[f_idx[1]],
            flux_tmp[f_idx[2]], flux_tmp[UEDEN], flux_tmp[UEINT], tmp0, tmp1,
            tmp2, tmp3, tmp4);
        } else {
          laxfriedrich_flux(
            qtempl[R_RHO], qtempl[R_UN], qtempl[R_UT1], qtempl[R_UT2],
            qtempl[R_P], spl, qtempr[R_RHO], qtempr[R_UN], qtempr[R_UT1],
            qtempr[R_UT2], qtempr[R_P], spr, bc_test_val, cavg, ustar,
            maxeigval, flux_tmp[URHO], flux_tmp[f_idx[0]], flux_tmp[f_idx[1]],
            flux_tmp[f_idx[2]], flux_tmp[UEDEN], flux_tmp[UEINT]);
        }

        if (!use_laxf_flux) {
          for (int n = 0; n < NUM_SPECIES; n++) {
            flux_tmp[UFS + n] = (ustar > 0.0)
                                  ? flux_tmp[URHO] * qtempl[R_Y + n]
                                  : flux_tmp[URHO] * qtempr[R_Y + n];
            flux_tmp[UFS + n] =
              (ustar == 0.0)
                ? flux_tmp[URHO] * 0.5 * (qtempl[R_Y + n] + qtempr[R_Y + n])
                : flux_tmp[UFS + n];
          }
        } else {
          for (int n = 0; n < NUM_SPECIES; n++) {

            // central
            flux_tmp[UFS + n] =
              0.5 * (qtempl[R_RHO] * qtempl[R_UN] * qtempl[R_Y + n] +
                     qtempr[R_RHO] * qtempr[R_UN] * qtempr[R_Y + n]);

            // dissipation
            flux_tmp[UFS + n] += -0.5 * maxeigval *
                                 (qtempr[R_RHO] * qtempr[R_Y + n] -
                                  qtempl[R_RHO] * qtempl[R_Y + n]);
          }
        }

        flux_tmp[UTEMP] = 0.0;
        for (int n = UFX; n < UFX + NUM_AUX; n++) {
          flux_tmp[n] = (NUM_AUX > 0) ? 0.0 : flux_tmp[n];
        }
        for (int n = UFA; n < UFA + NUM_ADV; n++) {
          flux_tmp[n] = (NUM_ADV > 0) ? 0.0 : flux_tmp[n];
        }

        for (int ivar = 0; ivar < NVAR; ivar++) {
          flx[dir](iv, ivar) += flux_tmp[ivar] * area[dir](i, j, k);
        }
      });
  }

#ifdef PELEC_USE_EB
  // nextra was 3 for EB in PeleC but we are operating on a different
  // box here, so this should be zero.
  const int nextra = 0;

  const amrex::Real full_area = std::pow(del[0], AMREX_SPACEDIM - 1);
  const amrex::Box bxg = amrex::grow(cbox, nextra - 1);

  amrex::ParallelFor(nebflux, [=] AMREX_GPU_DEVICE(int L) {
    const amrex::IntVect& iv = ebg[L].iv;
    if (bxg.contains(iv)) {
      amrex::Real ebnorm[AMREX_SPACEDIM] = {AMREX_D_DECL(
        ebg[L].eb_normal[0], ebg[L].eb_normal[1], ebg[L].eb_normal[2])};
      const amrex::Real ebnorm_mag = std::sqrt(AMREX_D_TERM(
        ebnorm[0] * ebnorm[0], +ebnorm[1] * ebnorm[1], +ebnorm[2] * ebnorm[2]));
      for (amrex::Real& dir : ebnorm) {
        dir /= ebnorm_mag;
      }

      amrex::Real flux_tmp[NVAR] = {0.0};
      AMREX_D_TERM(flux_tmp[UMX] = -q(iv, QPRES) * ebnorm[0];
                   , flux_tmp[UMY] = -q(iv, QPRES) * ebnorm[1];
                   , flux_tmp[UMZ] = -q(iv, QPRES) * ebnorm[2];)

      // Copy result into ebflux vector. Being a bit chicken here and only
      // copy values where ebg % iv is within box
      for (int n = 0; n < NVAR; n++) {
        ebflux[n * nebflux + L] += flux_tmp[n] * ebg[L].eb_area * full_area;
      }
    }
  });
#endif
}
