#include "mechanism.H"

#include "PelePhysics.H"
#include "Derive.H"
#include "IndexDefines.H"

void
pc_dervelx(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  auto const dat = datfab.array();
  auto velx = derfab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    velx(i, j, k) = dat(i, j, k, UMX) / dat(i, j, k, URHO);
  });
}

void
pc_dervely(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  auto const dat = datfab.array();
  auto vely = derfab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    vely(i, j, k) = dat(i, j, k, UMY) / dat(i, j, k, URHO);
  });
}

void
pc_dervelz(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  auto const dat = datfab.array();
  auto velz = derfab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    velz(i, j, k) = dat(i, j, k, UMZ) / dat(i, j, k, URHO);
  });
}

void
pc_dermagvel(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  auto const dat = datfab.array();
  auto magvel = derfab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real datinv = 1.0 / dat(i, j, k, URHO);
    const amrex::Real dat1 = (dat(i, j, k, UMX) * datinv);
    const amrex::Real dat2 = (dat(i, j, k, UMY) * datinv);
    const amrex::Real dat3 = (dat(i, j, k, UMZ) * datinv);
    magvel(i, j, k) = sqrt((dat1 * dat1) + (dat2 * dat2) + (dat3 * dat3));
  });
}

void
pc_dermagmom(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  auto const dat = datfab.array();
  auto magmom = derfab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    magmom(i, j, k) = sqrt(
      dat(i, j, k, UMX) * dat(i, j, k, UMX) +
      dat(i, j, k, UMY) * dat(i, j, k, UMY) +
      dat(i, j, k, UMZ) * dat(i, j, k, UMZ));
  });
}

void
pc_derkineng(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  auto const dat = datfab.array();
  auto kineng = derfab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real datxsq = dat(i, j, k, UMX) * dat(i, j, k, UMX);
    const amrex::Real datysq = dat(i, j, k, UMY) * dat(i, j, k, UMY);
    const amrex::Real datzsq = dat(i, j, k, UMZ) * dat(i, j, k, UMZ);
    kineng(i, j, k) = 0.5 / dat(i, j, k, URHO) * (datxsq + datysq + datzsq);
  });
}

void
pc_dereint1(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  // Compute internal energy from (rho E).
  auto const dat = datfab.array();
  auto e = derfab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real rhoInv = 1.0 / dat(i, j, k, URHO);
    const amrex::Real ux = dat(i, j, k, UMX) * rhoInv;
    const amrex::Real uy = dat(i, j, k, UMY) * rhoInv;
    const amrex::Real uz = dat(i, j, k, UMZ) * rhoInv;
    e(i, j, k) =
      dat(i, j, k, UEDEN) * rhoInv - 0.5 * (ux * ux + uy * uy + uz * uz);
  });
}

void
pc_dereint2(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  // Compute internal energy from (rho e).
  auto const dat = datfab.array();
  auto e = derfab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    e(i, j, k) = dat(i, j, k, UEINT) / dat(i, j, k, URHO);
  });
}

void
pc_derlogden(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  auto const dat = datfab.array();
  auto logden = derfab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    logden(i, j, k) = log10(dat(i, j, k));
  });
}

void
pc_derspec(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  auto const dat = datfab.array();
  auto spec = derfab.array();

  amrex::ParallelFor(
    bx, NUM_SPECIES, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
      spec(i, j, k, n) = dat(i, j, k, UFS + n) / dat(i, j, k, URHO);
    });
}

void
pc_dermagvort(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  int /*level*/)
{
  auto const dat = datfab.array();
  auto vort = derfab.array();

  const amrex::Box& gbx = amrex::grow(bx, 1);

  amrex::FArrayBox local(gbx, 3, amrex::The_Async_Arena());
  auto larr = local.array();

  // Convert momentum to velocity.
  amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real rhoInv = 1.0 / dat(i, j, k, URHO);
    AMREX_D_TERM(larr(i, j, k, 0) = dat(i, j, k, UMX) * rhoInv;
                 , larr(i, j, k, 1) = dat(i, j, k, UMY) * rhoInv;
                 , larr(i, j, k, 2) = dat(i, j, k, UMZ) * rhoInv;)
  });

  AMREX_D_TERM(const amrex::Real dx = geomdata.CellSize(0);
               , const amrex::Real dy = geomdata.CellSize(1);
               , const amrex::Real dz = geomdata.CellSize(2););

  // Calculate vorticity.
  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    AMREX_D_TERM(vort(i, j, k) = 0.0 * dx;
                 , const amrex::Real vx =
                     0.5 * (larr(i + 1, j, k, 1) - larr(i - 1, j, k, 1)) / dx;
                 const amrex::Real uy =
                   0.5 * (larr(i, j + 1, k, 0) - larr(i, j - 1, k, 0)) / dy;
                 const amrex::Real v3 = vx - uy;
                 , const amrex::Real wx =
                     0.5 * (larr(i + 1, j, k, 2) - larr(i - 1, j, k, 2)) / dx;

                 const amrex::Real wy =
                   0.5 * (larr(i, j + 1, k, 2) - larr(i, j - 1, k, 2)) / dy;

                 const amrex::Real uz =
                   0.5 * (larr(i, j, k + 1, 0) - larr(i, j, k - 1, 0)) / dz;
                 const amrex::Real vz =
                   0.5 * (larr(i, j, k + 1, 1) - larr(i, j, k - 1, 1)) / dz;

                 const amrex::Real v1 = wy - vz;
                 const amrex::Real v2 = uz - wx;);
    vort(i, j, k) = sqrt(AMREX_D_TERM(0., +v3 * v3, +v1 * v1 + v2 * v2));
  });
}

void
pc_derdivu(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  int /*level*/)
{
  auto const dat = datfab.array();
  auto divu = derfab.array();

  AMREX_D_TERM(const amrex::Real dx = geomdata.CellSize(0);
               , const amrex::Real dy = geomdata.CellSize(1);
               , const amrex::Real dz = geomdata.CellSize(2););

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    AMREX_D_TERM(
      const amrex::Real uhi = dat(i + 1, j, k, UMX) / dat(i + 1, j, k, URHO);
      const amrex::Real ulo = dat(i - 1, j, k, UMX) / dat(i - 1, j, k, URHO);
      , const amrex::Real vhi = dat(i, j + 1, k, UMY) / dat(i, j + 1, k, URHO);
      const amrex::Real vlo = dat(i, j - 1, k, UMY) / dat(i, j - 1, k, URHO);
      , const amrex::Real whi = dat(i, j, k + 1, UMZ) / dat(i, j, k + 1, URHO);
      const amrex::Real wlo = dat(i, j, k - 1, UMZ) / dat(i, j, k - 1, URHO););
    divu(i, j, k) =
      0.5 *
      (AMREX_D_TERM((uhi - ulo) / dx, +(vhi - vlo) / dy, +(whi - wlo) / dz));
  });
}

void
pc_derenstrophy(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  int /*level*/)
{
  // This routine will derive enstrophy  = 1/2 rho (x_vorticity^2 +
  // y_vorticity^2 + z_vorticity^2)
  auto const dat = datfab.array();
  auto enstrophy = derfab.array();

  const amrex::Box& gbx = amrex::grow(bx, 1);

  amrex::FArrayBox local(gbx, 3, amrex::The_Async_Arena());
  auto larr = local.array();

  // Convert momentum to velocity.
  amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real rhoInv = 1.0 / dat(i, j, k, URHO);
    larr(i, j, k, 0) = dat(i, j, k, UMX) * rhoInv;
    larr(i, j, k, 1) = dat(i, j, k, UMY) * rhoInv;
    larr(i, j, k, 2) = dat(i, j, k, UMZ) * rhoInv;
  });

  AMREX_D_TERM(const amrex::Real dx = geomdata.CellSize(0);
               , const amrex::Real dy = geomdata.CellSize(1);
               , const amrex::Real dz = geomdata.CellSize(2););

  // Calculate enstrophy.
  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    AMREX_D_TERM(enstrophy(i, j, k) = 0.0 * dx;
                 , const amrex::Real vx =
                     0.5 * (larr(i + 1, j, k, 1) - larr(i - 1, j, k, 1)) / dx;
                 const amrex::Real uy =
                   0.5 * (larr(i, j + 1, k, 0) - larr(i, j - 1, k, 0)) / dy;
                 const amrex::Real v3 = vx - uy;
                 , const amrex::Real wx =
                     0.5 * (larr(i + 1, j, k, 2) - larr(i - 1, j, k, 2)) / dx;

                 const amrex::Real wy =
                   0.5 * (larr(i, j + 1, k, 2) - larr(i, j - 1, k, 2)) / dy;

                 const amrex::Real uz =
                   0.5 * (larr(i, j, k + 1, 0) - larr(i, j, k - 1, 0)) / dz;
                 const amrex::Real vz =
                   0.5 * (larr(i, j, k + 1, 1) - larr(i, j, k - 1, 1)) / dz;

                 const amrex::Real v1 = wy - vz;
                 const amrex::Real v2 = uz - wx;);
    enstrophy(i, j, k) = 0.5 * dat(i, j, k, URHO) *
                         (AMREX_D_TERM(0., +v3 * v3, +v1 * v1 + v2 * v2));
  });
}

void
pc_dernull(
  const amrex::Box& /*bx*/,
  amrex::FArrayBox& /*derfab*/,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& /*datfab*/,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  // This routine is used by particle_count.  Yes it does nothing.
}

void
pc_dermolefrac(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  // Derive the mole fractions of the species
  auto const dat = datfab.array();
  auto spec = derfab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    amrex::Real mass[NUM_SPECIES];
    amrex::Real mole[NUM_SPECIES];
    const amrex::Real rhoInv = 1.0 / dat(i, j, k, URHO);

    for (int n = 0; n < NUM_SPECIES; n++) {
      mass[n] = dat(i, j, k, UFS + n) * rhoInv;
    }
    auto eos = pele::physics::PhysicsType::eos();
    eos.Y2X(mass, mole);
    for (int n = 0; n < NUM_SPECIES; n++) {
      spec(i, j, k, UFS + n) = mole[n];
    }
  });
}

void
pc_dersoundspeed(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  auto const dat = datfab.array();
  auto cfab = derfab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real rho = dat(i, j, k, URHO);
    const amrex::Real rhoInv = 1.0 / rho;
    const amrex::Real T = dat(i, j, k, UTEMP);
    amrex::Real massfrac[NUM_SPECIES];
    amrex::Real c;
    for (int n = 0; n < NUM_SPECIES; ++n) {
      massfrac[n] = dat(i, j, k, UFS + n) * rhoInv;
    }
    auto eos = pele::physics::PhysicsType::eos();
    eos.RTY2Cs(rho, T, massfrac, c);
    cfab(i, j, k) = c;
  });
}

void
pc_derentropy(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& /*datfab*/,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  // auto const dat = datfab.array();
  auto sfab = derfab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    amrex::Real s;
    auto eos = pele::physics::PhysicsType::eos();
    eos.S(s);
    sfab(i, j, k) = s;
  });
}

void
pc_dermachnumber(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  auto const dat = datfab.array();
  auto mach = derfab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real rho = dat(i, j, k, URHO);
    const amrex::Real rhoInv = 1.0 / rho;
    const amrex::Real T = dat(i, j, k, UTEMP);
    amrex::Real massfrac[NUM_SPECIES];
    amrex::Real c;
    for (int n = 0; n < NUM_SPECIES; ++n) {
      massfrac[n] = dat(i, j, k, UFS + n) * rhoInv;
    }
    auto eos = pele::physics::PhysicsType::eos();
    eos.RTY2Cs(rho, T, massfrac, c);
    const amrex::Real datxsq = dat(i, j, k, UMX) * dat(i, j, k, UMX);
    const amrex::Real datysq = dat(i, j, k, UMY) * dat(i, j, k, UMY);
    const amrex::Real datzsq = dat(i, j, k, UMZ) * dat(i, j, k, UMZ);
    mach(i, j, k) = sqrt(datxsq + datysq + datzsq) / dat(i, j, k, URHO) / c;
  });
}

void
pc_derpres(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  auto const dat = datfab.array();
  auto pfab = derfab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real rho = dat(i, j, k, URHO);
    const amrex::Real rhoInv = 1.0 / rho;
    amrex::Real T = dat(i, j, k, UTEMP);
    // amrex::Real e = dat(i, j, k, UEINT) * rhoInv;
    amrex::Real p;
    amrex::Real massfrac[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; ++n) {
      massfrac[n] = dat(i, j, k, UFS + n) * rhoInv;
    }
    auto eos = pele::physics::PhysicsType::eos();
    eos.RTY2P(rho, T, massfrac, p);
    pfab(i, j, k) = p;
  });
}

void
pc_dertemp(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  auto const dat = datfab.array();
  auto tfab = derfab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    tfab(i, j, k) = dat(i, j, k, UTEMP);
  });
}

void
pc_derspectrac(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/,
  const int idx)
{
  auto const dat = datfab.array();
  auto spectrac = derfab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    spectrac(i, j, k) = dat(i, j, k, UFS + idx) / dat(i, j, k, URHO);
  });
}

void
pc_derradialvel(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  int /*level*/)
{
  auto const dat = datfab.array();
  auto rvel = derfab.array();

  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_lo =
    geomdata.ProbLoArray();
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_hi =
    geomdata.ProbHiArray();
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx =
    geomdata.CellSizeArray();
  AMREX_D_TERM(const amrex::Real centerx = 0.5 * (prob_lo[0] + prob_hi[0]);
               , const amrex::Real centery = 0.5 * (prob_lo[1] + prob_hi[1]);
               , const amrex::Real centerz = 0.5 * (prob_lo[2] + prob_hi[2]));

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    AMREX_D_TERM(
      const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0] - centerx;
      , const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1] - centery;
      , const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2] - centerz;)
    const amrex::Real r = sqrt(AMREX_D_TERM(x * x, +y * y, +z * z));
    rvel(i, j, k) = (AMREX_D_TERM(
                      dat(i, j, k, UMX) * x, +dat(i, j, k, UMY) * y,
                      +dat(i, j, k, UMZ) * z)) /
                    (dat(i, j, k, URHO) * r);
  });
}

void
pc_dercp(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  int /*level*/)
{
  auto const dat = datfab.array();
  auto cp_arr = derfab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    amrex::Real mass[NUM_SPECIES];
    const amrex::Real rhoInv = 1.0 / dat(i, j, k, URHO);

    for (int n = 0; n < NUM_SPECIES; n++) {
      mass[n] = dat(i, j, k, UFS + n) * rhoInv;
    }
    auto eos = pele::physics::PhysicsType::eos();
    amrex::Real cp = 0.0;
    eos.RTY2Cp(dat(i, j, k, URHO), dat(i, j, k, UTEMP), mass, cp);
    cp_arr(i, j, k) = cp;
  });
}

void
pc_dercv(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  int /*level*/)
{
  auto const dat = datfab.array();
  auto cv_arr = derfab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    amrex::Real mass[NUM_SPECIES];
    const amrex::Real rhoInv = 1.0 / dat(i, j, k, URHO);

    for (int n = 0; n < NUM_SPECIES; n++) {
      mass[n] = dat(i, j, k, UFS + n) * rhoInv;
    }
    auto eos = pele::physics::PhysicsType::eos();
    amrex::Real cv = 0.0;
    eos.RTY2Cv(dat(i, j, k, URHO), dat(i, j, k, UTEMP), mass, cv);
    cv_arr(i, j, k) = cv;
  });
}

#ifdef PELEC_USE_MASA
void
pc_derrhommserror(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  int /*level*/)
{
  auto const dat = datfab.array();
  auto rhommserror = derfab.array();

  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_lo =
    geomdata.ProbLoArray();
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx =
    geomdata.CellSizeArray();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
    const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
    const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

    const amrex::Real rho = masa_eval_3d_exact_rho(x, y, z);
    rhommserror(i, j, k) = dat(i, j, k, URHO) - rho;
  });
}

void
pc_derummserror(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  int /*level*/)
{
  auto const dat = datfab.array();
  auto ummserror = derfab.array();

  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_lo =
    geomdata.ProbLoArray();
  // const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_hi =
  // geomdata.ProbHiArray();
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx =
    geomdata.CellSizeArray();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
    const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
    const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

    const amrex::Real u = masa_eval_3d_exact_u(x, y, z);
    ummserror(i, j, k) = dat(i, j, k, UMX) / dat(i, j, k, URHO) - u;
  });
}

void
pc_dervmmserror(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  int /*level*/)
{
  auto const dat = datfab.array();
  auto vmmserror = derfab.array();

  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_lo =
    geomdata.ProbLoArray();
  // const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_hi =
  // geomdata.ProbHiArray();
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx =
    geomdata.CellSizeArray();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
    const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
    const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

    const amrex::Real v = masa_eval_3d_exact_v(x, y, z);
    vmmserror(i, j, k) = dat(i, j, k, UMY) / dat(i, j, k, URHO) - v;
  });
}

void
pc_derwmmserror(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  int /*level*/)
{
  auto const dat = datfab.array();
  auto wmmserror = derfab.array();

  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_lo =
    geomdata.ProbLoArray();
  // const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_hi =
  // geomdata.ProbHiArray();
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx =
    geomdata.CellSizeArray();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
    const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
    const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

    const amrex::Real w = masa_eval_3d_exact_w(x, y, z);
    wmmserror(i, j, k) = dat(i, j, k, UMZ) / dat(i, j, k, URHO) - w;
  });
}

void
pc_derpmmserror(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  int /*level*/)
{
  auto const dat = datfab.array();
  auto pmmserror = derfab.array();

  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_lo =
    geomdata.ProbLoArray();
  // const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_hi =
  // geomdata.ProbHiArray();
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx =
    geomdata.CellSizeArray();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
    const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
    const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

    auto eos = pele::physics::PhysicsType::eos();
    const amrex::Real rho = dat(i, j, k, URHO);
    const amrex::Real rhoinv = 1.0 / rho;
    amrex::Real eint = dat(i, j, k, UEINT) * rhoinv;
    amrex::Real T = dat(i, j, k, UTEMP);
    amrex::Real massfrac[NUM_SPECIES] = {0.0};
    for (int n = 0; n < NUM_SPECIES; n++) {
      massfrac[n] = dat(i, j, k, UFS + n) * rhoinv;
    }

    amrex::Real pdat;
    eos.RYET2P(rho, massfrac, eint, T, pdat);
    const amrex::Real p = masa_eval_3d_exact_p(x, y, z);
    pmmserror(i, j, k) = pdat - p;
  });
}
#endif
