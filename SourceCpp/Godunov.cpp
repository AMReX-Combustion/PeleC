#include "Godunov.H"
#include "PLM.H"
#include "PPM.H"

// Host function to call gpu hydro functions
void
pc_umeth_3D(
  amrex::Box const& bx,
  const int* bclo,
  const int* bchi,
  const int* domlo,
  const int* domhi,
  amrex::Array4<const amrex::Real> const& q,
  amrex::Array4<const amrex::Real> const& qaux,
  amrex::Array4<const amrex::Real> const&
    srcQ, // amrex::IArrayBox const& bcMask,
  amrex::Array4<amrex::Real> const& flx1,
  amrex::Array4<amrex::Real> const& flx2,
  amrex::Array4<amrex::Real> const&
    flx3, // amrex::Array4<const amrex::Real> const& dloga,
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
  const int use_flattening)
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

  // auto const& bcMaskarr = bcMask.array();
  const amrex::Box& bxg1 = grow(bx, 1);
  const amrex::Box& bxg2 = grow(bx, 2);

  // X data
  int cdir = 0;
  const amrex::Box& xmbx = growHi(bxg2, cdir, 1);
  amrex::FArrayBox qxm(xmbx, QVAR);
  amrex::FArrayBox qxp(bxg2, QVAR);
  amrex::Elixir qxmeli = qxm.elixir();
  amrex::Elixir qxpeli = qxp.elixir();
  auto const& qxmarr = qxm.array();
  auto const& qxparr = qxp.array();

  // Y data
  cdir = 1;
  const amrex::Box& yflxbx = surroundingNodes(grow(bxg2, cdir, -1), cdir);
  const amrex::Box& ymbx = growHi(bxg2, cdir, 1);
  amrex::FArrayBox qym(ymbx, QVAR);
  amrex::FArrayBox qyp(bxg2, QVAR);
  amrex::Elixir qymeli = qym.elixir();
  amrex::Elixir qypeli = qyp.elixir();
  auto const& qymarr = qym.array();
  auto const& qyparr = qyp.array();

  // Z data
  cdir = 2;
  const amrex::Box& zmbx = growHi(bxg2, cdir, 1);
  const amrex::Box& zflxbx = surroundingNodes(grow(bxg2, cdir, -1), cdir);
  amrex::FArrayBox qzm(zmbx, QVAR);
  amrex::FArrayBox qzp(bxg2, QVAR);
  amrex::Elixir qzmeli = qzm.elixir();
  amrex::Elixir qzpeli = qzp.elixir();
  auto const& qzmarr = qzm.array();
  auto const& qzparr = qzp.array();

  // Put the PLM and slopes in the same kernel launch to avoid unnecessary
  // launch overhead Pelec_Slope_* are SIMD as well as PeleC_plm_* which loop
  // over the same box
  if (ppm_type == 0) {
    amrex::ParallelFor(
      bxg2, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        amrex::Real slope[QVAR];
        // X slopes and interp
        for (int n = 0; n < QVAR; ++n)
          slope[n] = plm_slope(i, j, k, n, 0, q);
        pc_plm_x(i, j, k, qxmarr, qxparr, slope, q, qaux(i, j, k, QC), dx, dt);

        // Y slopes and interp
        for (int n = 0; n < QVAR; n++)
          slope[n] = plm_slope(i, j, k, n, 1, q);
        pc_plm_y(i, j, k, qymarr, qyparr, slope, q, qaux(i, j, k, QC), dy, dt);

        // Z slopes and interp
        for (int n = 0; n < QVAR; ++n)
          slope[n] = plm_slope(i, j, k, n, 2, q);
        pc_plm_z(i, j, k, qzmarr, qzparr, slope, q, qaux(i, j, k, QC), dz, dt);
      });
  } else if (ppm_type == 1) {
    // Compute the normal interface states by reconstructing
    // the primitive variables using the piecewise parabolic method
    // and doing characteristic tracing.  We do not apply the
    // transverse terms here.

    int idir = 0;
    trace_ppm(bxg1, idir, q, srcQ, qxmarr, qxparr, bx, dt, del, use_flattening);

    idir = 1;
    trace_ppm(bxg1, idir, q, srcQ, qymarr, qyparr, bx, dt, del, use_flattening);

    idir = 2;
    trace_ppm(bxg1, idir, q, srcQ, qzmarr, qzparr, bx, dt, del, use_flattening);

  } else {
    amrex::Error("PeleC::ppm_type must be 0 (PLM) or 1 (PPM)");
  }

  // These are the first flux estimates as per the corner-transport-upwind
  // method X initial fluxes
  cdir = 0;
  const amrex::Box& xflxbx = surroundingNodes(grow(bxg2, cdir, -1), cdir);
  amrex::FArrayBox fx(xflxbx, NVAR);
  amrex::Elixir fxeli = fx.elixir();
  auto const& fxarr = fx.array();
  amrex::FArrayBox qgdx(xflxbx, NGDNV);
  amrex::Elixir qgdxeli = qgdx.elixir();
  auto const& gdtempx = qgdx.array();
  amrex::ParallelFor(
    xflxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      pc_cmpflx(
        i, j, k, bclx, bchx, dlx, dhx, qxmarr, qxparr, fxarr, gdtempx, qaux,
        cdir);
    });

  // Y initial fluxes
  cdir = 1;
  amrex::FArrayBox fy(yflxbx, NVAR);
  amrex::Elixir fyeli = fy.elixir();
  auto const& fyarr = fy.array();
  amrex::FArrayBox qgdy(yflxbx, NGDNV);
  amrex::Elixir qgdyeli = qgdy.elixir();
  auto const& gdtempy = qgdy.array();
  amrex::ParallelFor(
    yflxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      pc_cmpflx(
        i, j, k, bcly, bchy, dly, dhy, qymarr, qyparr, fyarr, gdtempy, qaux,
        cdir);
    });

  // Z initial fluxes
  cdir = 2;
  amrex::FArrayBox fz(zflxbx, NVAR);
  amrex::Elixir fzeli = fz.elixir();
  auto const& fzarr = fz.array();
  amrex::FArrayBox qgdz(zflxbx, NGDNV);
  amrex::Elixir qgdzeli = qgdz.elixir();
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
  amrex::FArrayBox qxym(txbxm, QVAR);
  amrex::Elixir qxymeli = qxym.elixir();
  amrex::FArrayBox qxyp(txbx, QVAR);
  amrex::Elixir qxypeli = qxyp.elixir();
  auto const& qmxy = qxym.array();
  auto const& qpxy = qxyp.array();

  amrex::FArrayBox qxzm(txbxm, QVAR);
  amrex::Elixir qxzmeli = qxzm.elixir();
  amrex::FArrayBox qxzp(txbx, QVAR);
  amrex::Elixir qxzpeli = qxzp.elixir();
  auto const& qmxz = qxzm.array();
  auto const& qpxz = qxzp.array();

  amrex::ParallelFor(txbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    // X|Y
    pc_transy1(
      i, j, k, qmxy, qpxy, qxmarr, qxparr, fyarr, qaux, gdtempy, cdtdy);
    // X|Z
    pc_transz1(
      i, j, k, qmxz, qpxz, qxmarr, qxparr, fzarr, qaux, gdtempz, cdtdz);
  });

  const amrex::Box& txfxbx = surroundingNodes(bxg1, cdir);
  amrex::FArrayBox fluxxy(txfxbx, NVAR);
  amrex::FArrayBox fluxxz(txfxbx, NVAR);
  amrex::FArrayBox gdvxyfab(txfxbx, NGDNV);
  amrex::FArrayBox gdvxzfab(txfxbx, NGDNV);
  amrex::Elixir fluxxyeli = fluxxy.elixir(), gdvxyeli = gdvxyfab.elixir();
  amrex::Elixir fluxxzeli = fluxxz.elixir(), gdvxzeli = gdvxzfab.elixir();

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

  qxymeli.clear();
  qxypeli.clear();
  qxzmeli.clear();
  qxzpeli.clear();

  // Y interface corrections
  cdir = 1;
  const amrex::Box& tybx = grow(bxg1, cdir, 1);
  const amrex::Box& tybxm = growHi(tybx, cdir, 1);
  amrex::FArrayBox qyxm(tybxm, QVAR);
  amrex::FArrayBox qyxp(tybx, QVAR);
  amrex::FArrayBox qyzm(tybxm, QVAR);
  amrex::FArrayBox qyzp(tybx, QVAR);
  amrex::Elixir qyxmeli = qyxm.elixir(), qyxpeli = qyxp.elixir();
  amrex::Elixir qyzmeli = qyzm.elixir(), qyzpeli = qyzp.elixir();
  auto const& qmyx = qyxm.array();
  auto const& qpyx = qyxp.array();
  auto const& qmyz = qyzm.array();
  auto const& qpyz = qyzp.array();

  amrex::ParallelFor(tybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    // Y|X
    pc_transx1(
      i, j, k, qmyx, qpyx, qymarr, qyparr, fxarr, qaux, gdtempx, cdtdx);
    // Y|Z
    pc_transz2(
      i, j, k, qmyz, qpyz, qymarr, qyparr, fzarr, qaux, gdtempz, cdtdz);
  });

  fzeli.clear();
  qgdzeli.clear();

  // Riemann problem Y|X Y|Z
  const amrex::Box& tyfxbx = surroundingNodes(bxg1, cdir);
  amrex::FArrayBox fluxyx(tyfxbx, NVAR);
  amrex::FArrayBox fluxyz(tyfxbx, NVAR);
  amrex::FArrayBox gdvyxfab(tyfxbx, NGDNV);
  amrex::FArrayBox gdvyzfab(tyfxbx, NGDNV);
  amrex::Elixir fluxyxeli = fluxyx.elixir(), gdvyxeli = gdvyxfab.elixir();
  amrex::Elixir fluxyzeli = fluxyz.elixir(), gdvyzeli = gdvyzfab.elixir();

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

  qyxmeli.clear();
  qyxpeli.clear();
  qyzmeli.clear();
  qyzpeli.clear();

  // Z interface corrections
  cdir = 2;
  const amrex::Box& tzbx = grow(bxg1, cdir, 1);
  const amrex::Box& tzbxm = growHi(tzbx, cdir, 1);
  amrex::FArrayBox qzxm(tzbxm, QVAR);
  amrex::FArrayBox qzxp(tzbx, QVAR);
  amrex::FArrayBox qzym(tzbxm, QVAR);
  amrex::FArrayBox qzyp(tzbx, QVAR);
  amrex::Elixir qzxmeli = qzxm.elixir(), qzxpeli = qzxp.elixir();
  amrex::Elixir qzymeli = qzym.elixir(), qzypeli = qzyp.elixir();

  auto const& qmzx = qzxm.array();
  auto const& qpzx = qzxp.array();
  auto const& qmzy = qzym.array();
  auto const& qpzy = qzyp.array();

  amrex::ParallelFor(tzbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    // Z|X
    pc_transx2(
      i, j, k, qmzx, qpzx, qzmarr, qzparr, fxarr, qaux, gdtempx, cdtdx);
    // Z|Y
    pc_transy2(
      i, j, k, qmzy, qpzy, qzmarr, qzparr, fyarr, qaux, gdtempy, cdtdy);
  });

  fxeli.clear();
  fyeli.clear();
  qgdxeli.clear();
  qgdyeli.clear();

  // Riemann problem Z|X Z|Y
  const amrex::Box& tzfxbx = surroundingNodes(bxg1, cdir);
  amrex::FArrayBox fluxzx(tzfxbx, NVAR);
  amrex::FArrayBox fluxzy(tzfxbx, NVAR);
  amrex::FArrayBox gdvzxfab(tzfxbx, NGDNV);
  amrex::FArrayBox gdvzyfab(tzfxbx, NGDNV);
  amrex::Elixir fluxzxeli = fluxzx.elixir(), gdvzxeli = gdvzxfab.elixir();
  amrex::Elixir fluxzyeli = fluxzy.elixir(), gdvzyeli = gdvzyfab.elixir();

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

  qzxmeli.clear();
  qzxpeli.clear();
  qzymeli.clear();
  qzypeli.clear();

  // Temp Fabs for Final Fluxes
  amrex::FArrayBox qmfab(bxg2, QVAR);
  amrex::FArrayBox qpfab(bxg1, QVAR);
  amrex::Elixir qmeli = qmfab.elixir();
  amrex::Elixir qpeli = qpfab.elixir();
  auto const& qm = qmfab.array();
  auto const& qp = qpfab.array();

  // X | Y&Z
  cdir = 0;
  const amrex::Box& xfxbx = surroundingNodes(bx, cdir);
  const amrex::Box& tyzbx = grow(bx, cdir, 1);
  amrex::ParallelFor(tyzbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    pc_transyz(
      i, j, k, qm, qp, qxmarr, qxparr, flyz, flzy, qyz, qzy, qaux, srcQ, hdt,
      hdtdy, hdtdz);
  });

  fluxzyeli.clear();
  gdvzyeli.clear();
  gdvyzeli.clear();
  fluxyzeli.clear();
  qxmeli.clear();
  qxpeli.clear();
  // Final X flux
  amrex::ParallelFor(xfxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    pc_cmpflx(i, j, k, bclx, bchx, dlx, dhx, qm, qp, flx1, q1, qaux, cdir);
  });

  // Y | X&Z
  cdir = 1;
  const amrex::Box& yfxbx = surroundingNodes(bx, cdir);
  const amrex::Box& txzbx = grow(bx, cdir, 1);
  amrex::ParallelFor(txzbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    pc_transxz(
      i, j, k, qm, qp, qymarr, qyparr, flxz, flzx, qxz, qzx, qaux, srcQ, hdt,
      hdtdx, hdtdz);
  });

  fluxzxeli.clear();
  gdvzxeli.clear();
  gdvxzeli.clear();
  fluxxzeli.clear();
  qymeli.clear();
  qypeli.clear();
  // Final Y flux
  amrex::ParallelFor(yfxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    pc_cmpflx(i, j, k, bcly, bchy, dly, dhy, qm, qp, flx2, q2, qaux, cdir);
  });

  // Z | X&Y
  cdir = 2;
  const amrex::Box& zfxbx = surroundingNodes(bx, cdir);
  const amrex::Box& txybx = grow(bx, cdir, 1);
  amrex::ParallelFor(txybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    pc_transxy(
      i, j, k, qm, qp, qzmarr, qzparr, flxy, flyx, qxy, qyx, qaux, srcQ, hdt,
      hdtdx, hdtdy);
  });

  gdvyxeli.clear();
  fluxyxeli.clear();
  gdvxyeli.clear();
  fluxxyeli.clear();
  qzmeli.clear();
  qzpeli.clear();
  // Final Z flux
  amrex::ParallelFor(zfxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    pc_cmpflx(i, j, k, bclz, bchz, dlz, dhz, qm, qp, flx3, q3, qaux, cdir);
  });

  qmeli.clear();
  qpeli.clear();
  // Construct p div{U}
  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    pc_pdivu(
      i, j, k, pdivu, AMREX_D_DECL(q1, q2, q3), AMREX_D_DECL(a1, a2, a3), vol);
  });
}

void
pc_umeth_2D(
  amrex::Box const& bx,
  const int* bclo,
  const int* bchi,
  const int* domlo,
  const int* domhi,
  amrex::Array4<const amrex::Real> const& q,
  amrex::Array4<const amrex::Real> const& qaux,
  amrex::Array4<const amrex::Real> const&
    srcQ, // amrex::IArrayBox const& bcMask,
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
  const int use_flattening)
{
#if AMREX_SPACEDIM == 2
  {
    const int use_flattening_loc = use_flattening;
    amrex::Real const dx = del[0];
    amrex::Real const dy = del[1];
    amrex::Real const hdtdy = 0.5 * dt / dy;
    amrex::Real const hdt = 0.5 * dt;

    const int bclx = bclo[0];
    const int bcly = bclo[1];
    const int bchx = bchi[0];
    const int bchy = bchi[1];
    const int dlx = domlo[0];
    const int dly = domlo[1];
    const int dhx = domhi[0];
    const int dhy = domhi[1];

    // auto const& bcMaskarr = bcMask.array();
    const amrex::Box& bxg1 = grow(bx, 1);
    const amrex::Box& bxg2 = grow(bx, 2);

    // X data
    int cdir = 0;
    const amrex::Box& xslpbx = grow(bxg1, cdir, 1);
    const amrex::Box& xmbx = growHi(xslpbx, cdir, 1);
    amrex::FArrayBox qxm(xmbx, QVAR);
    amrex::FArrayBox qxp(xslpbx, QVAR);
    amrex::Elixir qxmeli = qxm.elixir();
    amrex::Elixir qxpeli = qxp.elixir();
    auto const& qxmarr = qxm.array();
    auto const& qxparr = qxp.array();

    // Y data
    cdir = 1;
    const amrex::Box& yflxbx = surroundingNodes(bxg1, cdir);
    const amrex::Box& yslpbx = grow(bxg1, cdir, 1);
    const amrex::Box& ymbx = growHi(yslpbx, cdir, 1);
    amrex::FArrayBox qym(ymbx, QVAR);
    amrex::FArrayBox qyp(yslpbx, QVAR);
    amrex::Elixir qymeli = qym.elixir();
    amrex::Elixir qypeli = qyp.elixir();
    auto const& qymarr = qym.array();
    auto const& qyparr = qyp.array();

    amrex::ParallelFor(
      xslpbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        amrex::Real slope[QVAR];
        // X slopes and interp
        for (int n = 0; n < QVAR; ++n)
          slope[n] = plm_slope(i, j, k, n, 0, q);
        pc_plm_x(i, j, k, qxmarr, qxparr, slope, q, qaux(i, j, k, QC), dx, dt);
      });

    amrex::ParallelFor(
      yslpbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        amrex::Real slope[QVAR];
        // Y slopes and interp
        for (int n = 0; n < QVAR; n++)
          slope[n] = plm_slope(i, j, k, n, 1, q);
        pc_plm_y(i, j, k, qymarr, qyparr, slope, q, qaux(i, j, k, QC), dy, dt);
      });

    // These are the first flux estimates as per the corner-transport-upwind
    // method X initial fluxes
    cdir = 0;
    const amrex::Box& xflxbx = surroundingNodes(bxg1, cdir);
    amrex::FArrayBox fx(xflxbx, NVAR);
    amrex::Elixir fxeli = fx.elixir();
    auto const& fxarr = fx.array();
    amrex::FArrayBox qgdx(bxg2, NGDNV);
    amrex::Elixir qgdxeli = qgdx.elixir();
    auto const& gdtemp = qgdx.array();
    amrex::ParallelFor(
      xflxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        pc_cmpflx(
          i, j, k, bclx, bchx, dlx, dhx, qxmarr, qxparr, fxarr, gdtemp, qaux,
          cdir);
      });

    // Y initial fluxes
    cdir = 1;
    amrex::FArrayBox fy(yflxbx, NVAR);
    amrex::Elixir fyeli = fy.elixir();
    auto const& fyarr = fy.array();
    amrex::ParallelFor(
      yflxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        pc_cmpflx(
          i, j, k, bcly, bchy, dly, dhy, qymarr, qyparr, fyarr, q2, qaux, cdir);
      });

    // X interface corrections
    cdir = 0;
    const amrex::Box& tybx = grow(bx, cdir, 1);
    amrex::FArrayBox qm(bxg2, QVAR);
    amrex::Elixir qmeli = qm.elixir();
    amrex::FArrayBox qp(bxg1, QVAR);
    amrex::Elixir qpeli = qp.elixir();
    auto const& qmarr = qm.array();
    auto const& qparr = qp.array();

    amrex::ParallelFor(
      tybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        pc_transy(
          i, j, k, qmarr, qparr, qxmarr, qxparr, fyarr, srcQ, qaux, q2, hdt,
          hdtdy);
      });

    fyeli.clear();
    qxmeli.clear();
    qxpeli.clear();
    const amrex::Box& xfxbx = surroundingNodes(bx, cdir);

    // Final Riemann problem X
    amrex::ParallelFor(
      xfxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        pc_cmpflx(
          i, j, k, bclx, bchx, dlx, dhx, qmarr, qparr, flx1, q1, qaux, cdir);
      });

    // Y interface corrections
    cdir = 1;
    const amrex::Box& txbx = grow(bx, cdir, 1);

    amrex::ParallelFor(
      txbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        pc_transx(
          i, j, k, qmarr, qparr, qymarr, qyparr, fxarr, srcQ, qaux, gdtemp, a1,
          vol, hdt);
      });
    fxeli.clear();
    qymeli.clear();
    qypeli.clear();

    // Final Riemann problem Y
    const amrex::Box& yfxbx = surroundingNodes(bx, cdir);
    amrex::ParallelFor(
      yfxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        pc_cmpflx(
          i, j, k, bcly, bchy, dly, dhy, qmarr, qparr, flx2, q2, qaux, cdir);
      });

    // Construct p div{U}
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      pc_pdivu(i, j, k, pdivu, q1, q2, a1, a2, vol);
    });
#endif // AMREX_SPACEDIM == 2
  }
