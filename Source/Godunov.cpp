#include "Godunov.H"
#include "PLM.H"
#include "PPM.H"

// Host function that overwrites fluxes with lower-order approximation
void
pc_low_order_boundary(
  amrex::Box const& bfbx,
  const int bclo,
  const int bchi,
  const int domlo,
  const int domhi,
  const int plm_iorder,
  const bool use_flattening,
  const int idir,
  const amrex::Real dx,
  const amrex::Real dt,
  amrex::Array4<const amrex::Real> const& q,
  amrex::Array4<const amrex::Real> const& qa,
  amrex::Array4<amrex::Real> const& flx,
  amrex::Array4<amrex::Real> const& qdir)
{
  // Compute left and right states
  amrex::Box bdbx = enclosedCells(bfbx).grow(idir, 1);
  amrex::Box bdbx2 = growHi(bdbx, idir, 1);
  amrex::FArrayBox qbm(bdbx2, QVAR, amrex::The_Async_Arena());
  amrex::FArrayBox qbp(bdbx, QVAR, amrex::The_Async_Arena());
  auto const& qbmarr = qbm.array();
  auto const& qbparr = qbp.array();
  amrex::ParallelFor(bdbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    amrex::Real slope[QVAR];

    amrex::Real flat = 1.0;
    // Calculate flattening in-place
    if (use_flattening) {
      for (int dir_flat = 0; dir_flat < AMREX_SPACEDIM; dir_flat++) {
        flat = amrex::min<amrex::Real>(
          flat, flatten(AMREX_D_DECL(i, j, k), dir_flat, q));
      }
    }

    for (int n = 0; n < QVAR; ++n) {
      slope[n] = plm_slope(AMREX_D_DECL(i, j, k), n, idir, q, flat, plm_iorder);
    }
    pc_plm_d(
      AMREX_D_DECL(i, j, k), idir, qbmarr, qbparr, slope, q, qa(i, j, k, QC),
      dx, dt);
  });

  // Recompute fluxes
  amrex::ParallelFor(bfbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    pc_cmpflx(
      i, j, k, bclo, bchi, domlo, domhi, qbmarr, qbparr, flx, qdir, qa, idir);
  });
}

// Host function to call gpu hydro functions
#if AMREX_SPACEDIM == 3
void
pc_umeth_3D(
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
  const int weno_scheme)
{
  amrex::Real const dx = del[0];
  amrex::Real const dy = del[1];
  amrex::Real const dz = del[2];
  amrex::Real const hdtdx = 0.5 * dt / dx;
  amrex::Real const hdtdy = 0.5 * dt / dy;
  amrex::Real const hdtdz = 0.5 * dt / dz;
  amrex::Real const cdtdx = 1.0 / 3.0 * dt / dx;
  amrex::Real const cdtdy = 1.0 / 3.0 * dt / dy;
  amrex::Real const cdtdz = 1.0 / 3.0 * dt / dz;
  amrex::Real const hdt = 0.5 * dt;

  const int bclx = bclo[0];
  const int bcly = bclo[1];
  const int bclz = bclo[2];
  const int bchx = bchi[0];
  const int bchy = bchi[1];
  const int bchz = bchi[2];
  const int dlx = domlo[0];
  const int dly = domlo[1];
  const int dlz = domlo[2];
  const int dhx = domhi[0];
  const int dhy = domhi[1];
  const int dhz = domhi[2];

  const amrex::Box& bxg1 = grow(bx, 1);
  const amrex::Box& bxg2 = grow(bx, 2);

  // X data
  int cdir = 0;
  const amrex::Box& xmbx = growHi(bxg2, cdir, 1);
  const amrex::Box& xflxbx = surroundingNodes(grow(bxg2, cdir, -1), cdir);
  amrex::FArrayBox qxm(xmbx, QVAR, amrex::The_Async_Arena());
  amrex::FArrayBox qxp(bxg2, QVAR, amrex::The_Async_Arena());
  auto const& qxmarr = qxm.array();
  auto const& qxparr = qxp.array();

  // Y data
  cdir = 1;
  const amrex::Box& ymbx = growHi(bxg2, cdir, 1);
  const amrex::Box& yflxbx = surroundingNodes(grow(bxg2, cdir, -1), cdir);
  amrex::FArrayBox qym(ymbx, QVAR, amrex::The_Async_Arena());
  amrex::FArrayBox qyp(bxg2, QVAR, amrex::The_Async_Arena());
  auto const& qymarr = qym.array();
  auto const& qyparr = qyp.array();

  // Z data
  cdir = 2;
  const amrex::Box& zmbx = growHi(bxg2, cdir, 1);
  const amrex::Box& zflxbx = surroundingNodes(grow(bxg2, cdir, -1), cdir);
  amrex::FArrayBox qzm(zmbx, QVAR, amrex::The_Async_Arena());
  amrex::FArrayBox qzp(bxg2, QVAR, amrex::The_Async_Arena());
  auto const& qzmarr = qzm.array();
  auto const& qzparr = qzp.array();

  // Put the PLM and slopes in the same kernel launch to avoid unnecessary
  // launch
  if (ppm_type == 0) {
    amrex::ParallelFor(
      bxg2, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        amrex::Real slope[QVAR];

        amrex::Real flat = 1.0;
        // Calculate flattening in-place
        if (use_flattening) {
          for (int dir_flat = 0; dir_flat < AMREX_SPACEDIM; dir_flat++) {
            flat = amrex::min<amrex::Real>(flat, flatten(i, j, k, dir_flat, q));
          }
        }

        // X slopes and interp
        int idir = 0;
        for (int n = 0; n < QVAR; ++n) {
          slope[n] =
            plm_slope(AMREX_D_DECL(i, j, k), n, 0, q, flat, plm_iorder);
        }
        pc_plm_d(
          AMREX_D_DECL(i, j, k), idir, qxmarr, qxparr, slope, q,
          qaux(i, j, k, QC), dx, dt);

        // Y slopes and interp
        idir = 1;
        for (int n = 0; n < QVAR; n++) {
          slope[n] =
            plm_slope(AMREX_D_DECL(i, j, k), n, 1, q, flat, plm_iorder);
        }
        pc_plm_d(
          AMREX_D_DECL(i, j, k), idir, qymarr, qyparr, slope, q,
          qaux(i, j, k, QC), dy, dt);

        // Z slopes and interp
        idir = 2;
        for (int n = 0; n < QVAR; ++n) {
          slope[n] =
            plm_slope(AMREX_D_DECL(i, j, k), n, 2, q, flat, plm_iorder);
        }
        pc_plm_d(
          AMREX_D_DECL(i, j, k), idir, qzmarr, qzparr, slope, q,
          qaux(i, j, k, QC), dz, dt);
      });
  } else if (ppm_type == 1) {
    // Compute the normal interface states by reconstructing
    // the primitive variables using the piecewise parabolic method
    // and doing characteristic tracing.  We do not apply the
    // transverse terms here.

    int idir = 0;
    trace_ppm(
      bxg2, idir, q, srcQ, qxmarr, qxparr, bxg2, dt, del, use_flattening,
      use_hybrid_weno, weno_scheme);

    idir = 1;
    trace_ppm(
      bxg2, idir, q, srcQ, qymarr, qyparr, bxg2, dt, del, use_flattening,
      use_hybrid_weno, weno_scheme);

    idir = 2;
    trace_ppm(
      bxg2, idir, q, srcQ, qzmarr, qzparr, bxg2, dt, del, use_flattening,
      use_hybrid_weno, weno_scheme);

  } else {
    amrex::Error("PeleC::ppm_type must be 0 (PLM) or 1 (PPM)");
  }

  // These are the first flux estimates as per the corner-transport-upwind
  // method X initial fluxes
  cdir = 0;
  amrex::FArrayBox fx(xflxbx, NVAR, amrex::The_Async_Arena());
  auto const& fxarr = fx.array();
  amrex::FArrayBox qgdx(xflxbx, NGDNV, amrex::The_Async_Arena());
  auto const& gdtempx = qgdx.array();
  amrex::ParallelFor(
    xflxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      pc_cmpflx(
        i, j, k, bclx, bchx, dlx, dhx, qxmarr, qxparr, fxarr, gdtempx, qaux,
        cdir);
    });

  // Y initial fluxes
  cdir = 1;
  amrex::FArrayBox fy(yflxbx, NVAR, amrex::The_Async_Arena());
  auto const& fyarr = fy.array();
  amrex::FArrayBox qgdy(yflxbx, NGDNV, amrex::The_Async_Arena());
  auto const& gdtempy = qgdy.array();
  amrex::ParallelFor(
    yflxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      pc_cmpflx(
        i, j, k, bcly, bchy, dly, dhy, qymarr, qyparr, fyarr, gdtempy, qaux,
        cdir);
    });

  // Z initial fluxes
  cdir = 2;
  amrex::FArrayBox fz(zflxbx, NVAR, amrex::The_Async_Arena());
  auto const& fzarr = fz.array();
  amrex::FArrayBox qgdz(zflxbx, NGDNV, amrex::The_Async_Arena());
  auto const& gdtempz = qgdz.array();
  amrex::ParallelFor(
    zflxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      pc_cmpflx(
        i, j, k, bclz, bchz, dlz, dhz, qzmarr, qzparr, fzarr, gdtempz, qaux,
        cdir);
    });

  // X interface corrections
  cdir = 0;
  const amrex::Box& txbx = grow(bxg1, cdir, 1);
  const amrex::Box& txbxm = growHi(txbx, cdir, 1);
  amrex::FArrayBox qxym(txbxm, QVAR, amrex::The_Async_Arena());
  amrex::FArrayBox qxyp(txbx, QVAR, amrex::The_Async_Arena());
  auto const& qmxy = qxym.array();
  auto const& qpxy = qxyp.array();

  amrex::FArrayBox qxzm(txbxm, QVAR, amrex::The_Async_Arena());
  amrex::FArrayBox qxzp(txbx, QVAR, amrex::The_Async_Arena());
  auto const& qmxz = qxzm.array();
  auto const& qpxz = qxzp.array();

  amrex::ParallelFor(txbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    // X|Y
    pc_transdo(
      AMREX_D_DECL(i, j, k), cdir, 1, qmxy, qpxy, qxmarr, qxparr, fyarr, qaux,
      gdtempy, cdtdy);
    // X|Z
    pc_transdo(
      AMREX_D_DECL(i, j, k), cdir, 2, qmxz, qpxz, qxmarr, qxparr, fzarr, qaux,
      gdtempz, cdtdz);
  });

  const amrex::Box& txfxbx = surroundingNodes(bxg1, cdir);
  amrex::FArrayBox fluxxy(txfxbx, NVAR, amrex::The_Async_Arena());
  amrex::FArrayBox fluxxz(txfxbx, NVAR, amrex::The_Async_Arena());
  amrex::FArrayBox gdvxyfab(txfxbx, NGDNV, amrex::The_Async_Arena());
  amrex::FArrayBox gdvxzfab(txfxbx, NGDNV, amrex::The_Async_Arena());

  auto const& flxy = fluxxy.array();
  auto const& flxz = fluxxz.array();
  auto const& qxy = gdvxyfab.array();
  auto const& qxz = gdvxzfab.array();

  // Riemann problem X|Y X|Z
  amrex::ParallelFor(
    txfxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      // X|Y
      pc_cmpflx(
        i, j, k, bclx, bchx, dlx, dhx, qmxy, qpxy, flxy, qxy, qaux, cdir);
      // X|Z
      pc_cmpflx(
        i, j, k, bclx, bchx, dlx, dhx, qmxz, qpxz, flxz, qxz, qaux, cdir);
    });

  // Y interface corrections
  cdir = 1;
  const amrex::Box& tybx = grow(bxg1, cdir, 1);
  const amrex::Box& tybxm = growHi(tybx, cdir, 1);
  amrex::FArrayBox qyxm(tybxm, QVAR, amrex::The_Async_Arena());
  amrex::FArrayBox qyxp(tybx, QVAR, amrex::The_Async_Arena());
  amrex::FArrayBox qyzm(tybxm, QVAR, amrex::The_Async_Arena());
  amrex::FArrayBox qyzp(tybx, QVAR, amrex::The_Async_Arena());
  auto const& qmyx = qyxm.array();
  auto const& qpyx = qyxp.array();
  auto const& qmyz = qyzm.array();
  auto const& qpyz = qyzp.array();

  amrex::ParallelFor(tybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    // Y|X
    pc_transdo(
      AMREX_D_DECL(i, j, k), cdir, 0, qmyx, qpyx, qymarr, qyparr, fxarr, qaux,
      gdtempx, cdtdx);
    // Y|Z
    pc_transdo(
      AMREX_D_DECL(i, j, k), cdir, 2, qmyz, qpyz, qymarr, qyparr, fzarr, qaux,
      gdtempz, cdtdz);
  });

  // Riemann problem Y|X Y|Z
  const amrex::Box& tyfxbx = surroundingNodes(bxg1, cdir);
  amrex::FArrayBox fluxyx(tyfxbx, NVAR, amrex::The_Async_Arena());
  amrex::FArrayBox fluxyz(tyfxbx, NVAR, amrex::The_Async_Arena());
  amrex::FArrayBox gdvyxfab(tyfxbx, NGDNV, amrex::The_Async_Arena());
  amrex::FArrayBox gdvyzfab(tyfxbx, NGDNV, amrex::The_Async_Arena());

  auto const& flyx = fluxyx.array();
  auto const& flyz = fluxyz.array();
  auto const& qyx = gdvyxfab.array();
  auto const& qyz = gdvyzfab.array();

  amrex::ParallelFor(
    tyfxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      // Y|X
      pc_cmpflx(
        i, j, k, bcly, bchy, dly, dhy, qmyx, qpyx, flyx, qyx, qaux, cdir);
      // Y|Z
      pc_cmpflx(
        i, j, k, bcly, bchy, dly, dhy, qmyz, qpyz, flyz, qyz, qaux, cdir);
    });

  // Z interface corrections
  cdir = 2;
  const amrex::Box& tzbx = grow(bxg1, cdir, 1);
  const amrex::Box& tzbxm = growHi(tzbx, cdir, 1);
  amrex::FArrayBox qzxm(tzbxm, QVAR, amrex::The_Async_Arena());
  amrex::FArrayBox qzxp(tzbx, QVAR, amrex::The_Async_Arena());
  amrex::FArrayBox qzym(tzbxm, QVAR, amrex::The_Async_Arena());
  amrex::FArrayBox qzyp(tzbx, QVAR, amrex::The_Async_Arena());

  auto const& qmzx = qzxm.array();
  auto const& qpzx = qzxp.array();
  auto const& qmzy = qzym.array();
  auto const& qpzy = qzyp.array();

  amrex::ParallelFor(tzbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    // Z|X
    pc_transdo(
      AMREX_D_DECL(i, j, k), cdir, 0, qmzx, qpzx, qzmarr, qzparr, fxarr, qaux,
      gdtempx, cdtdx);
    // Z|Y
    pc_transdo(
      AMREX_D_DECL(i, j, k), cdir, 1, qmzy, qpzy, qzmarr, qzparr, fyarr, qaux,
      gdtempy, cdtdy);
  });

  // Riemann problem Z|X Z|Y
  const amrex::Box& tzfxbx = surroundingNodes(bxg1, cdir);
  amrex::FArrayBox fluxzx(tzfxbx, NVAR, amrex::The_Async_Arena());
  amrex::FArrayBox fluxzy(tzfxbx, NVAR, amrex::The_Async_Arena());
  amrex::FArrayBox gdvzxfab(tzfxbx, NGDNV, amrex::The_Async_Arena());
  amrex::FArrayBox gdvzyfab(tzfxbx, NGDNV, amrex::The_Async_Arena());

  auto const& flzx = fluxzx.array();
  auto const& flzy = fluxzy.array();
  auto const& qzx = gdvzxfab.array();
  auto const& qzy = gdvzyfab.array();

  amrex::ParallelFor(
    tzfxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      // Z|X
      pc_cmpflx(
        i, j, k, bclz, bchz, dlz, dhz, qmzx, qpzx, flzx, qzx, qaux, cdir);
      // Z|Y
      pc_cmpflx(
        i, j, k, bclz, bchz, dlz, dhz, qmzy, qpzy, flzy, qzy, qaux, cdir);
    });

  // Temp Fabs for Final Fluxes
  amrex::FArrayBox qmfab(bxg2, QVAR, amrex::The_Async_Arena());
  amrex::FArrayBox qpfab(bxg1, QVAR, amrex::The_Async_Arena());
  auto const& qm = qmfab.array();
  auto const& qp = qpfab.array();

  // X | Y&Z
  cdir = 0;
  const amrex::Box& xfxbx = surroundingNodes(bx, cdir);
  const amrex::Box& tyzbx = grow(bx, cdir, 1);
  amrex::ParallelFor(tyzbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    pc_transdd(
      AMREX_D_DECL(i, j, k), cdir, qm, qp, qxmarr, qxparr, flyz, flzy, qyz, qzy,
      qaux, srcQ, hdt, hdtdy, hdtdz);
  });

  // Final X flux
  amrex::ParallelFor(xfxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    pc_cmpflx(
      i, j, k, bclx, bchx, dlx, dhx, qm, qp, flx[0], qec[0], qaux, cdir);
  });

  // Y | X&Z
  cdir = 1;
  const amrex::Box& yfxbx = surroundingNodes(bx, cdir);
  const amrex::Box& txzbx = grow(bx, cdir, 1);
  amrex::ParallelFor(txzbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    pc_transdd(
      AMREX_D_DECL(i, j, k), cdir, qm, qp, qymarr, qyparr, flxz, flzx, qxz, qzx,
      qaux, srcQ, hdt, hdtdx, hdtdz);
  });

  // Final Y flux
  amrex::ParallelFor(yfxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    pc_cmpflx(
      i, j, k, bcly, bchy, dly, dhy, qm, qp, flx[1], qec[1], qaux, cdir);
  });

  // Z | X&Y
  cdir = 2;
  const amrex::Box& zfxbx = surroundingNodes(bx, cdir);
  const amrex::Box& txybx = grow(bx, cdir, 1);
  amrex::ParallelFor(txybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    pc_transdd(
      AMREX_D_DECL(i, j, k), cdir, qm, qp, qzmarr, qzparr, flxy, flyx, qxy, qyx,
      qaux, srcQ, hdt, hdtdx, hdtdy);
  });

  // Final Z flux
  amrex::ParallelFor(zfxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    pc_cmpflx(
      i, j, k, bclz, bchz, dlz, dhz, qm, qp, flx[2], qec[2], qaux, cdir);
  });

  // Fix bcnormal boundaries - always use PLM and don't do N+1/2 predictor
  // because the user specifies conditions at N
  // bcnormal used for "Hard" (inflow) and "UserBC" (user_bc)
  for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
    if (
      bclo[idir] == PCPhysBCType::inflow ||
      bclo[idir] == PCPhysBCType::user_bc) {
      // Box for fluxes at this boundary
      amrex::Box bfbx = surroundingNodes(bx, idir);
      bfbx.setBig(idir, domlo[idir]);
      if (bfbx.ok()) {
        pc_low_order_boundary(
          bfbx, bclo[idir], bchi[idir], domlo[idir], domhi[idir], plm_iorder,
          use_flattening, idir, del[idir], dt, q, qaux, flx[idir], qec[idir]);
      }
    }
    if (
      bchi[idir] == PCPhysBCType::inflow ||
      bchi[idir] == PCPhysBCType::user_bc) {
      // Box for fluxes at this boundary
      amrex::Box bfbx = surroundingNodes(bx, idir);
      bfbx.setSmall(idir, domhi[idir] + 1);
      if (bfbx.ok()) {
        pc_low_order_boundary(
          bfbx, bclo[idir], bchi[idir], domlo[idir], domhi[idir], plm_iorder,
          use_flattening, idir, del[idir], dt, q, qaux, flx[idir], qec[idir]);
      }
    }
  }

  // Construct p div{U}
  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    pc_pdivu(
      i, j, k, pdivu, AMREX_D_DECL(qec[0], qec[1], qec[2]),
      AMREX_D_DECL(a[0], a[1], a[2]), vol);
  });
}

void
pc_umeth_eb_3D(
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
  const int plm_iorder)
{
  int cdir;

  amrex::Real const dx = del[0];
  amrex::Real const dy = del[1];
  amrex::Real const dz = del[2];

  amrex::Real const hdt = amrex::Real(0.5) * dt;
  amrex::Real const hdtdx = hdt / dx;
  amrex::Real const hdtdy = hdt / dy;
  amrex::Real const hdtdz = hdt / dz;
  amrex::Real const cdtdx = amrex::Real(1.0 / 3.0) * dt / dx;
  amrex::Real const cdtdy = amrex::Real(1.0 / 3.0) * dt / dy;
  amrex::Real const cdtdz = amrex::Real(1.0 / 3.0) * dt / dz;

  const int bclx = bclo[0];
  const int bcly = bclo[1];
  const int bclz = bclo[2];
  const int bchx = bchi[0];
  const int bchy = bchi[1];
  const int bchz = bchi[2];
  const int dlx = domlo[0];
  const int dly = domlo[1];
  const int dlz = domlo[2];
  const int dhx = domhi[0];
  const int dhy = domhi[1];
  const int dhz = domhi[2];

  const amrex::Box& bxg1 = grow(bx_to_fill, 1);
  const amrex::Box& bxg2 = grow(bx_to_fill, 2);

  // X data
  cdir = 0;
  const amrex::Box& xmbx = growHi(bxg2, cdir, 1);
  amrex::FArrayBox qxm(xmbx, QVAR, amrex::The_Async_Arena());
  auto const& qxmarr = qxm.array();
  amrex::FArrayBox qxp(bxg2, QVAR, amrex::The_Async_Arena());
  auto const& qxparr = qxp.array();

  // Y data
  cdir = 1;
  const amrex::Box& ymbx = growHi(bxg2, cdir, 1);
  amrex::FArrayBox qym(ymbx, QVAR, amrex::The_Async_Arena());
  auto const& qymarr = qym.array();
  amrex::FArrayBox qyp(bxg2, QVAR, amrex::The_Async_Arena());
  auto const& qyparr = qyp.array();

  // Z data
  cdir = 2;
  const amrex::Box& zmbx = growHi(bxg2, cdir, 1);
  amrex::FArrayBox qzm(zmbx, QVAR, amrex::The_Async_Arena());
  auto const& qzmarr = qzm.array();
  amrex::FArrayBox qzp(bxg2, QVAR, amrex::The_Async_Arena());
  auto const& qzparr = qzp.array();
  //
  // Put the PLM and slopes in the same kernel launch to avoid unnecessary
  // launch overhead Note that we compute the qm and qp values on bxg2 -5,-5,-5
  //
  if (ppm_type == 0) {
    amrex::ParallelFor(
      bxg2, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        if (flag_arr(i, j, k).isCovered()) {
          return;
        }

        amrex::Real slope[QVAR];

        amrex::Real flat = 1.0;
        // Calculate flattening in-place
        if (use_flattening) {
          for (int dir_flat = 0; dir_flat < AMREX_SPACEDIM; dir_flat++) {
            flat = amrex::min<amrex::Real>(
              flat, flatten_eb(i, j, k, dir_flat, flag_arr, q));
          }
        }

        // X slopes and interp
        int idir = 0;
        for (int n = 0; n < QVAR; ++n) {
          slope[n] = plm_slope_eb(
            AMREX_D_DECL(i, j, k), n, 0, flag_arr, q, flat, plm_iorder);
        }
        // cells -5,-5,-5
        pc_plm_d_eb(
          i, j, k, idir, qxmarr, qxparr, slope, q, qaux(i, j, k, QC), dx, dt,
          ap[idir]);

        // Y slopes and interp
        idir = 1;
        for (int n = 0; n < QVAR; n++) {
          slope[n] = plm_slope_eb(
            AMREX_D_DECL(i, j, k), n, 1, flag_arr, q, flat, plm_iorder);
        }
        // cells -5,-5,-5
        pc_plm_d_eb(
          i, j, k, idir, qymarr, qyparr, slope, q, qaux(i, j, k, QC), dy, dt,
          ap[idir]);

        // Z slopes and interp
        idir = 2;
        for (int n = 0; n < QVAR; ++n) {
          slope[n] = plm_slope_eb(
            AMREX_D_DECL(i, j, k), n, 2, flag_arr, q, flat, plm_iorder);
        }
        // cells -5,-5,-5
        pc_plm_d_eb(
          i, j, k, idir, qzmarr, qzparr, slope, q, qaux(i, j, k, QC), dz, dt,
          ap[idir]);
      });
  } else {
    amrex::Error("ppm_type must be 0 (PLM) when using EB");
  }

  // These are the first flux estimates as per the corner-transport-upwind
  // method

  // *************************************************************************************
  // X initial fluxes
  // *************************************************************************************
  cdir = 0;
  const amrex::Box& xflxbx = surroundingNodes(grow(bxg2, cdir, -1), cdir);

  amrex::FArrayBox fx(xflxbx, NVAR, amrex::The_Async_Arena());
  auto const& fxarr = fx.array();

  amrex::FArrayBox qgdx(xflxbx, NGDNV, amrex::The_Async_Arena());
  auto const& gdtempx = qgdx.array();

  // -4,-5,-5
  amrex::ParallelFor(
    xflxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      if (ap[cdir](i, j, k) > 0.) {
        pc_cmpflx(
          i, j, k, bclx, bchx, dlx, dhx, qxmarr, qxparr, fxarr, gdtempx, qaux,
          cdir);
      }
    });

  // *************************************************************************************
  // Y initial fluxes
  // *************************************************************************************
  cdir = 1;
  const amrex::Box& yflxbx = surroundingNodes(grow(bxg2, cdir, -1), cdir);

  amrex::FArrayBox fy(yflxbx, NVAR, amrex::The_Async_Arena());
  auto const& fyarr = fy.array();

  amrex::FArrayBox qgdy(yflxbx, NGDNV, amrex::The_Async_Arena());
  auto const& gdtempy = qgdy.array();

  // -5,-4,-5
  amrex::ParallelFor(
    yflxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      if (ap[cdir](i, j, k) > 0.) {
        pc_cmpflx(
          i, j, k, bcly, bchy, dly, dhy, qymarr, qyparr, fyarr, gdtempy, qaux,
          cdir);
      }
    });

  // *************************************************************************************
  // Z initial fluxes
  // *************************************************************************************
  cdir = 2;
  const amrex::Box& zflxbx = surroundingNodes(grow(bxg2, cdir, -1), cdir);

  amrex::FArrayBox fz(zflxbx, NVAR, amrex::The_Async_Arena());
  auto const& fzarr = fz.array();

  amrex::FArrayBox qgdz(zflxbx, NGDNV, amrex::The_Async_Arena());
  auto const& gdtempz = qgdz.array();

  // -5,-5,-4
  amrex::ParallelFor(
    zflxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      if (ap[cdir](i, j, k) > 0.) {
        pc_cmpflx(
          i, j, k, bclz, bchz, dlz, dhz, qzmarr, qzparr, fzarr, gdtempz, qaux,
          cdir);
      }
    });

  // *************************************************************************************
  // X interface corrections
  // *************************************************************************************
  cdir = 0;
  amrex::FArrayBox qxym(xmbx, QVAR, amrex::The_Async_Arena());
  amrex::FArrayBox qxyp(xmbx, QVAR, amrex::The_Async_Arena());
  auto const& qmxy = qxym.array();
  auto const& qpxy = qxyp.array();

  amrex::FArrayBox qxzm(xmbx, QVAR, amrex::The_Async_Arena());
  amrex::FArrayBox qxzp(xmbx, QVAR, amrex::The_Async_Arena());
  auto const& qmxz = qxzm.array();
  auto const& qpxz = qxzp.array();

  //
  amrex::Box xybx(bxg2);
  xybx.grow(1, -1);
  amrex::ParallelFor(xybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    // fyarr was made on y-faces to -5,-4,-5
    // this loop is over cells to   -5,-4,-5
    if (!flag_arr(i, j, k).isCovered()) {
      // X|Y
      pc_transdo(
        i, j, k, cdir, 1, qmxy, qpxy, qxmarr, qxparr, fyarr, qaux, gdtempy,
        cdtdy, ap[cdir], ap[1]);
    }
  });

  amrex::Box xzbx(bxg2);
  xzbx.grow(2, -1);
  amrex::ParallelFor(xzbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    // fzarr was made on z-faces to -5,-5,-4
    // this loop is over cells to   -5,-5,-4
    if (!flag_arr(i, j, k).isCovered()) {
      // X|Z
      pc_transdo(
        i, j, k, cdir, 2, qmxz, qpxz, qxmarr, qxparr, fzarr, qaux, gdtempz,
        cdtdz, ap[cdir], ap[2]);
    }
  });

  const amrex::Box& txfxbx = surroundingNodes(bxg2, cdir);
  amrex::FArrayBox fluxxy(txfxbx, NVAR, amrex::The_Async_Arena());
  amrex::FArrayBox fluxxz(txfxbx, NVAR, amrex::The_Async_Arena());
  amrex::FArrayBox gdvxyfab(txfxbx, NGDNV, amrex::The_Async_Arena());
  amrex::FArrayBox gdvxzfab(txfxbx, NGDNV, amrex::The_Async_Arena());

  auto const& flxy = fluxxy.array();
  auto const& flxz = fluxxz.array();
  auto const& qxy = gdvxyfab.array();
  auto const& qxz = gdvxzfab.array();

  // Riemann problem X|Y
  amrex::Box xycmpbx(surroundingNodes(bxg1, 0).grow(2, 1));
  // this loop is over x-faces to   -4,-4,-5
  amrex::ParallelFor(
    xycmpbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      // X|Y
      if (ap[cdir](i, j, k) > 0.) {
        pc_cmpflx(
          i, j, k, bclx, bchx, dlx, dhx, qmxy, qpxy, flxy, qxy, qaux, cdir);
      }
    });

  // Riemann problem X|Z
  amrex::Box xzcmpbx(surroundingNodes(bxg1, 0).grow(1, 1));
  // this loop is over x-faces to   -4,-5,-4
  amrex::ParallelFor(
    xzcmpbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      // X|Z
      if (ap[cdir](i, j, k) > 0.) {
        pc_cmpflx(
          i, j, k, bclx, bchx, dlx, dhx, qmxz, qpxz, flxz, qxz, qaux, cdir);
      }
    });

  // *************************************************************************************
  // Y interface corrections
  // *************************************************************************************

  cdir = 1;
  amrex::FArrayBox qyxm(ymbx, QVAR, amrex::The_Async_Arena());
  amrex::FArrayBox qyxp(ymbx, QVAR, amrex::The_Async_Arena());
  auto const& qmyx = qyxm.array();
  auto const& qpyx = qyxp.array();

  amrex::FArrayBox qyzm(ymbx, QVAR, amrex::The_Async_Arena());
  amrex::FArrayBox qyzp(ymbx, QVAR, amrex::The_Async_Arena());
  auto const& qmyz = qyzm.array();
  auto const& qpyz = qyzp.array();

  amrex::Box yxbx(bxg2);
  yxbx.grow(0, -1);
  amrex::ParallelFor(yxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    // Y|X
    if (!flag_arr(i, j, k).isCovered()) {
      pc_transdo(
        i, j, k, cdir, 0, qmyx, qpyx, qymarr, qyparr, fxarr, qaux, gdtempx,
        cdtdx, ap[cdir], ap[0]);
    }
  });

  amrex::Box yzbx(bxg2);
  yzbx.grow(2, -1);
  amrex::ParallelFor(yzbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    // Y|Z
    if (!flag_arr(i, j, k).isCovered()) {
      pc_transdo(
        i, j, k, cdir, 2, qmyz, qpyz, qymarr, qyparr, fzarr, qaux, gdtempz,
        cdtdz, ap[cdir], ap[2]);
    }
  });

  // Riemann problem Y|X Y|Z
  const amrex::Box& tyfxbx = surroundingNodes(bxg2, cdir);
  amrex::FArrayBox fluxyx(tyfxbx, NVAR, amrex::The_Async_Arena());
  amrex::FArrayBox fluxyz(tyfxbx, NVAR, amrex::The_Async_Arena());
  amrex::FArrayBox gdvyxfab(tyfxbx, NGDNV, amrex::The_Async_Arena());
  amrex::FArrayBox gdvyzfab(tyfxbx, NGDNV, amrex::The_Async_Arena());

  auto const& flyx = fluxyx.array();
  auto const& flyz = fluxyz.array();
  auto const& qyx = gdvyxfab.array();
  auto const& qyz = gdvyzfab.array();

  // Riemann problem Y|X
  amrex::Box yxcmpbx(surroundingNodes(bxg1, 1).grow(2, 1));
  amrex::ParallelFor(
    yxcmpbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      // Y|X
      if (ap[cdir](i, j, k) > 0.) {
        pc_cmpflx(
          i, j, k, bcly, bchy, dly, dhy, qmyx, qpyx, flyx, qyx, qaux, cdir);
      }
    });

  // Riemann problem Y|Z
  amrex::Box yzcmpbx(surroundingNodes(bxg1, 1).grow(0, 1));
  amrex::ParallelFor(
    yzcmpbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      // Y|Z
      if (ap[cdir](i, j, k) > 0.) {
        pc_cmpflx(
          i, j, k, bcly, bchy, dly, dhy, qmyz, qpyz, flyz, qyz, qaux, cdir);
      }
    });

  // *************************************************************************************
  // Z interface corrections
  // *************************************************************************************
  cdir = 2;

  amrex::FArrayBox qzxm(zmbx, QVAR, amrex::The_Async_Arena());
  amrex::FArrayBox qzxp(zmbx, QVAR, amrex::The_Async_Arena());
  auto const& qmzx = qzxm.array();
  auto const& qpzx = qzxp.array();

  amrex::FArrayBox qzym(zmbx, QVAR, amrex::The_Async_Arena());
  amrex::FArrayBox qzyp(zmbx, QVAR, amrex::The_Async_Arena());
  auto const& qmzy = qzym.array();
  auto const& qpzy = qzyp.array();

  amrex::Box zxbx(bxg2);
  zxbx.grow(0, -1);
  amrex::ParallelFor(zxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    // Z|X
    if (!flag_arr(i, j, k).isCovered()) {
      pc_transdo(
        i, j, k, cdir, 0, qmzx, qpzx, qzmarr, qzparr, fxarr, qaux, gdtempx,
        cdtdx, ap[cdir], ap[0]);
    }
  });

  amrex::Box zybx(bxg2);
  zybx.grow(1, -1);
  amrex::ParallelFor(zybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    // Z|Y
    if (!flag_arr(i, j, k).isCovered()) {
      pc_transdo(
        i, j, k, cdir, 1, qmzy, qpzy, qzmarr, qzparr, fyarr, qaux, gdtempy,
        cdtdy, ap[cdir], ap[1]);
    }
  });

  // Riemann problem Z|X Z|Y
  const amrex::Box& tzfxbx = surroundingNodes(bxg2, cdir);
  amrex::FArrayBox fluxzx(tzfxbx, NVAR, amrex::The_Async_Arena());
  amrex::FArrayBox fluxzy(tzfxbx, NVAR, amrex::The_Async_Arena());
  amrex::FArrayBox gdvzxfab(tzfxbx, NGDNV, amrex::The_Async_Arena());
  amrex::FArrayBox gdvzyfab(tzfxbx, NGDNV, amrex::The_Async_Arena());

  auto const& flzx = fluxzx.array();
  auto const& flzy = fluxzy.array();
  auto const& qzx = gdvzxfab.array();
  auto const& qzy = gdvzyfab.array();

  // Riemann problem Z|X
  amrex::Box zxcmpbx(surroundingNodes(bxg1, 2).grow(1, 1));
  amrex::ParallelFor(
    zxcmpbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      // Z|X
      if (ap[cdir](i, j, k) > 0.) {
        pc_cmpflx(
          i, j, k, bclz, bchz, dlz, dhz, qmzx, qpzx, flzx, qzx, qaux, cdir);
      }
    });

  // Riemann problem Z|Y
  amrex::Box zycmpbx(surroundingNodes(bxg1, 2).grow(0, 1));
  amrex::ParallelFor(
    zycmpbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      // Z|Y
      if (ap[cdir](i, j, k) > 0.) {
        pc_cmpflx(
          i, j, k, bclz, bchz, dlz, dhz, qmzy, qpzy, flzy, qzy, qaux, cdir);
      }
    });

  amrex::FArrayBox qmfab(bxg2, QVAR, amrex::The_Async_Arena());
  amrex::FArrayBox qpfab(bxg1, QVAR, amrex::The_Async_Arena());
  auto const& qm = qmfab.array();
  auto const& qp = qpfab.array();

  // *************************************************************************************
  // Final X steps
  // *************************************************************************************
  cdir = 0;

  // this loop is over cells to   -4,-4,-4
  amrex::ParallelFor(bxg1, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    if (!flag_arr(i, j, k).isCovered()) {
      pc_transdd(
        i, j, k, cdir, qm, qp, qxmarr, qxparr, flyz, flzy, qyz, qzy, qaux, srcQ,
        hdt, hdtdy, hdtdz, ap[0], ap[1], ap[2]);
    }
  });

  // This box must be grown by one in transverse direction so we can
  // do tangential interpolation when taking divergence later
  amrex::Box xfbx(surroundingNodes(bx_to_fill, 0));
  xfbx.grow(1, 1);
  xfbx.grow(2, 1);
  amrex::ParallelFor(xfbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    if (ap[cdir](i, j, k) > 0.) {
      pc_cmpflx(
        i, j, k, bclx, bchx, dlx, dhx, qm, qp, flx[cdir], qec[cdir], qaux,
        cdir);
    }
  });

  // *************************************************************************************
  // Final Y steps
  // *************************************************************************************
  cdir = 1;
  // this loop is over cells to   -4,-4,-4
  amrex::ParallelFor(bxg1, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    if (!flag_arr(i, j, k).isCovered()) {
      pc_transdd(
        i, j, k, cdir, qm, qp, qymarr, qyparr, flxz, flzx, qxz, qzx, qaux, srcQ,
        hdt, hdtdx, hdtdz, ap[1], ap[0], ap[2]);
    }
  });

  // This box must be grown by one in transverse direction so we can
  // do tangential interpolation when taking divergence later
  amrex::Box yfbx(surroundingNodes(bx_to_fill, 1));
  yfbx.grow(0, 1);
  yfbx.grow(2, 1);
  amrex::ParallelFor(yfbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    if (ap[cdir](i, j, k) > 0.) {
      pc_cmpflx(
        i, j, k, bcly, bchy, dly, dhy, qm, qp, flx[cdir], qec[cdir], qaux,
        cdir);
    }
  });

  // *************************************************************************************
  // Final Z steps
  // *************************************************************************************
  cdir = 2;
  amrex::ParallelFor(bxg1, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    if (!flag_arr(i, j, k).isCovered()) {
      pc_transdd(
        i, j, k, cdir, qm, qp, qzmarr, qzparr, flxy, flyx, qxy, qyx, qaux, srcQ,
        hdt, hdtdx, hdtdy, ap[2], ap[0], ap[1]);
    }
  });

  // This box must be grown by one in transverse direction so we can
  // do tangential interpolation when taking divergence later
  amrex::Box zfbx(surroundingNodes(bx_to_fill, 2));
  zfbx.grow(0, 1);
  zfbx.grow(1, 1);
  amrex::ParallelFor(zfbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    if (ap[cdir](i, j, k) > 0.) {
      pc_cmpflx(
        i, j, k, bclz, bchz, dlz, dhz, qm, qp, flx[cdir], qec[cdir], qaux,
        cdir);
    }
  });

  // Fix bcnormal boundaries - always use PLM and don't do N+1/2 predictor
  // because the user specifies conditions at N
  // bcnormal used for "Hard" (inflow) and "UserBC" (user_bc)
  for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
    if (
      bclo[idir] == PCPhysBCType::inflow ||
      bclo[idir] == PCPhysBCType::user_bc) {
      // Box for fluxes at this boundary
      amrex::Box bfbx = surroundingNodes(bx_to_fill, idir);
      bfbx.setBig(idir, domlo[idir]);
      if (bfbx.ok()) {
        pc_low_order_boundary(
          bfbx, bclo[idir], bchi[idir], domlo[idir], domhi[idir], plm_iorder,
          use_flattening, idir, del[idir], dt, q, qaux, flx[idir], qec[idir]);
      }
    }
    if (
      bchi[idir] == PCPhysBCType::inflow ||
      bchi[idir] == PCPhysBCType::user_bc) {
      // Box for fluxes at this boundary
      amrex::Box bfbx = surroundingNodes(bx_to_fill, idir);
      bfbx.setSmall(idir, domhi[idir] + 1);
      if (bfbx.ok()) {
        pc_low_order_boundary(
          bfbx, bclo[idir], bchi[idir], domlo[idir], domhi[idir], plm_iorder,
          use_flattening, idir, del[idir], dt, q, qaux, flx[idir], qec[idir]);
      }
    }
  }
}

#elif AMREX_SPACEDIM == 2

void
pc_umeth_2D(
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
  const int weno_scheme)
{
  amrex::Real const dx = del[0];
  amrex::Real const dy = del[1];
  amrex::Real const hdt = 0.5 * dt;
  amrex::Real const hdtdy = 0.5 * dt / dy;
  amrex::Real const hdtdx = 0.5 * dt / dx;

  const int bclx = bclo[0];
  const int bcly = bclo[1];
  const int bchx = bchi[0];
  const int bchy = bchi[1];
  const int dlx = domlo[0];
  const int dly = domlo[1];
  const int dhx = domhi[0];
  const int dhy = domhi[1];

  const amrex::Box& bxg1 = grow(bx, 1);
  const amrex::Box& bxg2 = grow(bx, 2);

  // X data
  int cdir = 0;
  const amrex::Box& xmbx = growHi(bxg2, cdir, 1);
  const amrex::Box& xflxbx = surroundingNodes(grow(bxg2, cdir, -1), cdir);
  amrex::FArrayBox qxm(xmbx, QVAR, amrex::The_Async_Arena());
  amrex::FArrayBox qxp(bxg2, QVAR, amrex::The_Async_Arena());
  auto const& qxmarr = qxm.array();
  auto const& qxparr = qxp.array();

  // Y data
  cdir = 1;
  const amrex::Box& ymbx = growHi(bxg2, cdir, 1);
  const amrex::Box& yflxbx = surroundingNodes(grow(bxg2, cdir, -1), cdir);
  amrex::FArrayBox qym(ymbx, QVAR, amrex::The_Async_Arena());
  amrex::FArrayBox qyp(bxg2, QVAR, amrex::The_Async_Arena());
  auto const& qymarr = qym.array();
  auto const& qyparr = qyp.array();

  if (ppm_type == 0) {
    amrex::ParallelFor(
      bxg2, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        amrex::Real slope[QVAR];

        amrex::Real flat = 1.0;
        // Calculate flattening in-place
        if (use_flattening) {
          for (int dir_flat = 0; dir_flat < AMREX_SPACEDIM; dir_flat++) {
            flat = amrex::min<amrex::Real>(
              flat, flatten(AMREX_D_DECL(i, j, k), dir_flat, q));
          }
        }

        // X slopes and interp
        for (int n = 0; n < QVAR; ++n)
          slope[n] =
            plm_slope(AMREX_D_DECL(i, j, k), n, 0, q, flat, plm_iorder);
        pc_plm_d(
          AMREX_D_DECL(i, j, k), 0, qxmarr, qxparr, slope, q, qaux(i, j, k, QC),
          dx, dt);

        // Y slopes and interp
        for (int n = 0; n < QVAR; n++)
          slope[n] =
            plm_slope(AMREX_D_DECL(i, j, k), n, 1, q, flat, plm_iorder);
        pc_plm_d(
          AMREX_D_DECL(i, j, k), 1, qymarr, qyparr, slope, q, qaux(i, j, k, QC),
          dy, dt);
      });
  } else if (ppm_type == 1) {
    // Compute the normal interface states by reconstructing
    // the primitive variables using the piecewise parabolic method
    // and doing characteristic tracing.  We do not apply the
    // transverse terms here.

    int idir = 0;
    trace_ppm(
      bxg2, idir, q, srcQ, qxmarr, qxparr, bxg2, dt, del, use_flattening,
      use_hybrid_weno, weno_scheme);

    idir = 1;
    trace_ppm(
      bxg2, idir, q, srcQ, qymarr, qyparr, bxg2, dt, del, use_flattening,
      use_hybrid_weno, weno_scheme);

  } else {
    amrex::Error("PeleC::ppm_type must be 0 (PLM) or 1 (PPM)");
  }
  // These are the first flux estimates as per the corner-transport-upwind
  // method X initial fluxes
  cdir = 0;
  amrex::FArrayBox fx(xflxbx, NVAR, amrex::The_Async_Arena());
  auto const& fxarr = fx.array();
  amrex::FArrayBox qgdx(bxg2, NGDNV, amrex::The_Async_Arena());
  auto const& gdtemp = qgdx.array();
  amrex::ParallelFor(
    xflxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      pc_cmpflx(
        i, j, k, bclx, bchx, dlx, dhx, qxmarr, qxparr, fxarr, gdtemp, qaux,
        cdir);
    });

  // Y initial fluxes
  cdir = 1;
  amrex::FArrayBox fy(yflxbx, NVAR, amrex::The_Async_Arena());
  auto const& fyarr = fy.array();
  amrex::ParallelFor(
    yflxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      pc_cmpflx(
        i, j, k, bcly, bchy, dly, dhy, qymarr, qyparr, fyarr, qec[1], qaux,
        cdir);
    });

  // X interface corrections
  cdir = 0;
  const amrex::Box& tybx = grow(bx, cdir, 1);
  amrex::FArrayBox qm(bxg2, QVAR, amrex::The_Async_Arena());
  amrex::FArrayBox qp(bxg1, QVAR, amrex::The_Async_Arena());
  auto const& qmarr = qm.array();
  auto const& qparr = qp.array();

  amrex::ParallelFor(
    tybx, [=] AMREX_GPU_DEVICE(
            int i, int j, AMREX_D_PICK(int /*k*/, int /*k*/, int k)) noexcept {
      pc_transd(
        AMREX_D_DECL(i, j, k), cdir, qmarr, qparr, qxmarr, qxparr, fyarr, srcQ,
        qaux, qec[1], hdt, hdtdy);
    });

  const amrex::Box& xfxbx = surroundingNodes(bx, cdir);

  // Final Riemann problem X
  amrex::ParallelFor(xfxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    pc_cmpflx(
      i, j, k, bclx, bchx, dlx, dhx, qmarr, qparr, flx[0], qec[0], qaux, cdir);
  });

  // Y interface corrections
  cdir = 1;
  const amrex::Box& txbx = grow(bx, cdir, 1);

  amrex::ParallelFor(
    txbx, [=] AMREX_GPU_DEVICE(
            int i, int j, AMREX_D_PICK(int /*k*/, int /*k*/, int k)) noexcept {
      pc_transd(
        AMREX_D_DECL(i, j, k), cdir, qmarr, qparr, qymarr, qyparr, fxarr, srcQ,
        qaux, gdtemp, hdt, hdtdx);
    });

  // Final Riemann problem Y
  const amrex::Box& yfxbx = surroundingNodes(bx, cdir);
  amrex::ParallelFor(yfxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    pc_cmpflx(
      i, j, k, bcly, bchy, dly, dhy, qmarr, qparr, flx[1], qec[1], qaux, cdir);
  });

  // Fix bcnormal boundaries - always use PLM and don't do N+1/2 predictor
  // because the user specifies conditions at N
  // bcnormal used for "Hard" (inflow) and "UserBC" (user_bc)
  for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
    if (
      bclo[idir] == PCPhysBCType::inflow ||
      bclo[idir] == PCPhysBCType::user_bc) {
      // Box for fluxes at this boundary
      amrex::Box bfbx = surroundingNodes(bx, idir);
      bfbx.setBig(idir, domlo[idir]);
      if (bfbx.ok()) {
        pc_low_order_boundary(
          bfbx, bclo[idir], bchi[idir], domlo[idir], domhi[idir], plm_iorder,
          use_flattening, idir, del[idir], dt, q, qaux, flx[idir], qec[idir]);
      }
    }
    if (
      bchi[idir] == PCPhysBCType::inflow ||
      bchi[idir] == PCPhysBCType::user_bc) {
      // Box for fluxes at this boundary
      amrex::Box bfbx = surroundingNodes(bx, idir);
      bfbx.setSmall(idir, domhi[idir] + 1);
      if (bfbx.ok()) {
        pc_low_order_boundary(
          bfbx, bclo[idir], bchi[idir], domlo[idir], domhi[idir], plm_iorder,
          use_flattening, idir, del[idir], dt, q, qaux, flx[idir], qec[idir]);
      }
    }
  }

  // Construct p div{U}
  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    pc_pdivu(i, j, k, pdivu, qec[0], qec[1], a[0], a[1], vol);
  });
}

void
pc_umeth_eb_2D(
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
  const int plm_iorder)
{
  BL_PROFILE("Godunov_umeth_2D_eb()");

  amrex::Real const dx = del[0];
  amrex::Real const dy = del[1];
  amrex::Real const hdt = 0.5 * dt;
  amrex::Real const hdtdy = 0.5 * dt / dy;
  amrex::Real const hdtdx = 0.5 * dt / dx;

  const int bclx = bclo[0];
  const int bcly = bclo[1];
  const int bchx = bchi[0];
  const int bchy = bchi[1];
  const int dlx = domlo[0];
  const int dly = domlo[1];
  const int dhx = domhi[0];
  const int dhy = domhi[1];

  const amrex::Box& bxg1 = grow(bx_to_fill, 1);
  const amrex::Box& bxg2 = grow(bx_to_fill, 2);

  // X data
  int cdir = 0;
  const amrex::Box& xmbx = growHi(bxg2, cdir, 1);
  const amrex::Box& xflxbx = surroundingNodes(bxg1, cdir);
  amrex::FArrayBox qxm(xmbx, QVAR, amrex::The_Async_Arena());
  amrex::FArrayBox qxp(bxg2, QVAR, amrex::The_Async_Arena());
  auto const& qxmarr = qxm.array();
  auto const& qxparr = qxp.array();

  // Y data
  cdir = 1;
  const amrex::Box& ymbx = growHi(bxg2, cdir, 1);
  const amrex::Box& yflxbx = surroundingNodes(bxg1, cdir);
  amrex::FArrayBox qym(ymbx, QVAR, amrex::The_Async_Arena());
  amrex::FArrayBox qyp(bxg2, QVAR, amrex::The_Async_Arena());
  auto const& qymarr = qym.array();
  auto const& qyparr = qyp.array();

  if (ppm_type == 0) {
    amrex::ParallelFor(
      bxg2, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        if (flag_arr(i, j, k).isCovered()) {
          return;
        }

        amrex::Real slope[QVAR];

        amrex::Real flat = 1.0;
        // Calculate flattening in-place
        if (use_flattening) {
          for (int dir_flat = 0; dir_flat < AMREX_SPACEDIM; dir_flat++) {
            flat = amrex::min<amrex::Real>(
              flat, flatten_eb(i, j, k, dir_flat, flag_arr, q));
          }
        }

        //
        // X slopes and interp
        //
        for (int n = 0; n < QVAR; ++n) {
          slope[n] = plm_slope_eb(
            AMREX_D_DECL(i, j, k), n, 0, flag_arr, q, flat, plm_iorder);
        }
        pc_plm_d_eb(
          i, j, k, 0, qxmarr, qxparr, slope, q, qaux(i, j, k, QC), dx, dt,
          ap[0]);

        //
        // Y slopes and interp
        //
        for (int n = 0; n < QVAR; n++) {
          slope[n] = plm_slope_eb(
            AMREX_D_DECL(i, j, k), n, 1, flag_arr, q, flat, plm_iorder);
        }
        pc_plm_d_eb(
          i, j, k, 1, qymarr, qyparr, slope, q, qaux(i, j, k, QC), dy, dt,
          ap[1]);
      });
  } else {
    amrex::Error("ppm_type must be 0 (PLM) when using EB");
  }

  // *******************************************************************************
  // These are the first flux estimates as per the corner-transport-upwind
  // method X initial fluxes
  // *******************************************************************************
  cdir = 0;
  amrex::FArrayBox fx(xflxbx, NVAR, amrex::The_Async_Arena());
  auto const& fxarr = fx.array();
  amrex::FArrayBox qgdx(bxg2, NGDNV, amrex::The_Async_Arena());
  auto const& gdtemp = qgdx.array();
  amrex::ParallelFor(
    xflxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      if (ap[cdir](i, j, k) > 0.) {
        pc_cmpflx(
          i, j, k, bclx, bchx, dlx, dhx, qxmarr, qxparr, fxarr, gdtemp, qaux,
          cdir);
      }
    });

  // *******************************************************************************
  // These are the first flux estimates as per the corner-transport-upwind
  // method Y initial fluxes
  // *******************************************************************************
  cdir = 1;
  amrex::FArrayBox fy(yflxbx, NVAR, amrex::The_Async_Arena());
  auto const& fyarr = fy.array();
  amrex::ParallelFor(
    yflxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      if (ap[cdir](i, j, k) > 0.) {
        pc_cmpflx(
          i, j, k, bcly, bchy, dly, dhy, qymarr, qyparr, fyarr, qec[cdir], qaux,
          cdir);
      }
    });

  // *******************************************************************************
  // X interface corrections
  // *******************************************************************************
  cdir = 0;
  const amrex::Box& tybx = bxg1;
  amrex::FArrayBox qm(bxg2, QVAR, amrex::The_Async_Arena());
  amrex::FArrayBox qp(bxg1, QVAR, amrex::The_Async_Arena());
  auto const& qmarr = qm.array();
  auto const& qparr = qp.array();

  amrex::ParallelFor(tybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    if (!flag_arr(i, j, k).isCovered()) {
      pc_transd(
        AMREX_D_DECL(i, j, k), cdir, qmarr, qparr, qxmarr, qxparr, fyarr, srcQ,
        qaux, qec[1], hdt, hdtdy, ap[cdir], ap[1]);
    }
  });

  // This box must be grown by one in transverse direction so we can
  // do tangential interpolation when taking divergence later
  const amrex::Box& xfxbx =
    surroundingNodes(grow(bx_to_fill, 1, 1 - cdir), cdir);

  // Final Riemann problem X
  amrex::ParallelFor(xfxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    if (ap[cdir](i, j, k) > 0.) {
      pc_cmpflx(
        i, j, k, bclx, bchx, dlx, dhx, qmarr, qparr, flx[cdir], qec[cdir], qaux,
        cdir);
    }
  });

  // *******************************************************************************
  // Y interface corrections
  // *******************************************************************************
  cdir = 1;
  const amrex::Box& txbx = bxg1;

  amrex::ParallelFor(txbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    if (!flag_arr(i, j, k).isCovered()) {
      pc_transd(
        AMREX_D_DECL(i, j, k), cdir, qmarr, qparr, qymarr, qyparr, fxarr, srcQ,
        qaux, gdtemp, hdt, hdtdx, ap[cdir], ap[0]);
    }
  });

  // This box must be grown by one in transverse direction so we can
  // do tangential interpolation when taking divergence later
  const amrex::Box& yfxbx =
    surroundingNodes(grow(bx_to_fill, 1, 1 - cdir), cdir);

  // Final Riemann problem Y
  amrex::ParallelFor(yfxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    if (ap[cdir](i, j, k) > 0.) {
      pc_cmpflx(
        i, j, k, bcly, bchy, dly, dhy, qmarr, qparr, flx[cdir], qec[cdir], qaux,
        cdir);
    }
  });

  // Fix bcnormal boundaries - always use PLM and don't do N+1/2 predictor
  // because the user specifies conditions at N
  // bcnormal used for "Hard" (inflow) and "UserBC" (user_bc)
  for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
    if (
      bclo[idir] == PCPhysBCType::inflow ||
      bclo[idir] == PCPhysBCType::user_bc) {
      // Box for fluxes at this boundary
      amrex::Box bfbx = surroundingNodes(bx_to_fill, idir);
      bfbx.setBig(idir, domlo[idir]);
      if (bfbx.ok()) {
        pc_low_order_boundary(
          bfbx, bclo[idir], bchi[idir], domlo[idir], domhi[idir], plm_iorder,
          use_flattening, idir, del[idir], dt, q, qaux, flx[idir], qec[idir]);
      }
    }
    if (
      bchi[idir] == PCPhysBCType::inflow ||
      bchi[idir] == PCPhysBCType::user_bc) {
      // Box for fluxes at this boundary
      amrex::Box bfbx = surroundingNodes(bx_to_fill, idir);
      bfbx.setSmall(idir, domhi[idir] + 1);
      if (bfbx.ok()) {
        pc_low_order_boundary(
          bfbx, bclo[idir], bchi[idir], domlo[idir], domhi[idir], plm_iorder,
          use_flattening, idir, del[idir], dt, q, qaux, flx[idir], qec[idir]);
      }
    }
  }
}

#endif
