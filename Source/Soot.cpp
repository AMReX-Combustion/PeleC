#include "PeleC.H"
#include "SootModel.H"
#include "SootModel_derive.H"
#include "Transport.H"

void
PeleC::setSootIndx()
{
  SootComps sc;
  sc.qRhoIndx = QRHO;
  sc.qTempIndx = QTEMP;
  sc.qSpecIndx = QFS;
  sc.qSootIndx = QFSOOT;
  sc.rhoIndx = URHO;
  sc.engIndx = UEDEN;
  sc.specIndx = UFS;
  sc.sootIndx = UFSOOT;
  soot_model.setIndices(sc);
}

void
PeleC::addSootDerivePlotVars(
  amrex::DeriveList& soot_derive_lst,
  const amrex::DescriptorList& soot_desc_lst)
{
  // Add in soot variables
  amrex::Vector<std::string> sootNames = {"rho_soot", "sum_rho_soot"};
  soot_derive_lst.add(
    "soot_vars", amrex::IndexType::TheCellType(),
    static_cast<int>(sootNames.size()), sootNames, soot_genvars,
    amrex::DeriveRec::TheSameBox);
  soot_derive_lst.addComponent("soot_vars", soot_desc_lst, State_Type, URHO, 1);
  soot_derive_lst.addComponent(
    "soot_vars", soot_desc_lst, State_Type, UFSOOT, NUM_SOOT_MOMENTS + 1);

  // Variables associated with the second mode (large particles)
  amrex::Vector<std::string> large_part_names = {"NL", "soot_V_L", "soot_S_L"};
  soot_derive_lst.add(
    "soot_large_particles", amrex::IndexType::TheCellType(),
    static_cast<int>(large_part_names.size()), large_part_names,
    soot_largeparticledata, amrex::DeriveRec::TheSameBox);
  soot_derive_lst.addComponent(
    "soot_large_particles", soot_desc_lst, State_Type, UFSOOT,
    NUM_SOOT_MOMENTS + 1);
}

void
PeleC::construct_old_soot_source(amrex::Real time, amrex::Real dt)
{
  old_sources[soot_src]->setVal(0.0);
  if (!add_soot_src) {
    return;
  }
  amrex::MultiFab& S_old = get_old_data(State_Type);

  int ng = 0; // None filled

  PeleC::fill_soot_source(time, dt, S_old, *old_sources[soot_src], ng);

  old_sources[soot_src]->FillBoundary(geom.periodicity());
}

void
PeleC::construct_new_soot_source(amrex::Real time, amrex::Real dt)
{
  new_sources[soot_src]->setVal(0.0);
  if (!add_soot_src) {
    return;
  }
  amrex::MultiFab& S_new = get_new_data(State_Type);

  int ng = 0;

  PeleC::fill_soot_source(time, dt, S_new, *new_sources[soot_src], ng);
}

void
PeleC::fill_soot_source(
  amrex::Real time,
  amrex::Real dt,
  amrex::MultiFab& state,
  amrex::MultiFab& soot_src,
  int ng)
{
  BL_PROFILE("PeleC::fill_soot_source()");

  auto const& fact =
    dynamic_cast<amrex::EBFArrayBoxFactory const&>(state.Factory());
  auto const& flags = fact.getMultiEBCellFlagFab();

  // TODO: Change to use new ParallelFor type
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
  for (amrex::MFIter mfi(soot_src, amrex::TilingIfNotGPU()); mfi.isValid();
       ++mfi) {
    const amrex::Box& bx = mfi.growntilebox(ng);
    const auto& flag_fab = flags[mfi];
    amrex::FabType typ = flag_fab.getType(bx);
    if (typ == amrex::FabType::covered) {
      continue;
    }
    amrex::FArrayBox& Sfab = state[mfi];
    amrex::FArrayBox& soot_fab = soot_src[mfi];
    auto const& s_arr = Sfab.array();
    const int nqaux = NQAUX > 0 ? NQAUX : 1;
    amrex::FArrayBox mu_cc(bx, 1, amrex::The_Async_Arena());
    amrex::FArrayBox q(bx, QVAR, amrex::The_Async_Arena());
    amrex::FArrayBox qaux(bx, nqaux, amrex::The_Async_Arena());
    auto const& q_arr = q.array();
    auto const& qaux_arr = qaux.array();
    auto const& mu_arr = mu_cc.array();

    // Get primitives, Q, including (Y, T, p, rho) from conserved state
    // required for D term
    {
      BL_PROFILE("PeleC::ctoprim()");
      amrex::ParallelFor(
        bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          pc_ctoprim(i, j, k, s_arr, q_arr, qaux_arr);
        });
    }

    // Compute transport coefficients, coincident with Q
    {
      auto const& qar_yin = q.array(QFS);
      auto const& qar_Tin = q.array(QTEMP);
      auto const& qar_rhoin = q.array(QRHO);
      bool get_xi = false;
      bool get_mu = true;
      bool get_lam = false;
      bool get_diag = false;
      bool get_chi = false;
      BL_PROFILE("PeleC::get_transport_coeffs()");
      // Get Transport coefs on GPU.
      auto const* ltransparm = trans_parms.device_trans_parm();
      amrex::ParallelFor(
        bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          auto trans = pele::physics::PhysicsType::transport();
          amrex::Real T = qar_Tin(i, j, k);
          amrex::Real rho = qar_rhoin(i, j, k);
          amrex::GpuArray<amrex::Real, NUM_SPECIES> Y = {{0.0}};
          for (int n = 0; n < NUM_SPECIES; ++n) {
            Y[n] = qar_yin(i, j, k, n);
          }
          amrex::Real* diag = nullptr;
          amrex::Real mu = 0.;
          amrex::Real xi, lam;
          trans.transport(
            get_xi, get_mu, get_lam, get_diag, get_chi, T, rho, Y.data(), diag,
            nullptr, mu, xi, lam, ltransparm);
          mu_arr(i, j, k) = mu;
        });
    }
    auto const& soot_arr = soot_fab.array();
    soot_model.computeSootSourceTerm(bx, q_arr, mu_arr, soot_arr, time, dt);
  }
}

void
PeleC::clipSootMoments(amrex::MultiFab& S_new, const int ng)
{
  SootData* sd = soot_model.getSootData_d();
  auto const& s_arr = S_new.arrays();
  const amrex::IntVect ngs(ng);
  amrex::ParallelFor(
    S_new, ngs, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
      amrex::GpuArray<amrex::Real, NUM_SOOT_MOMENTS + 1> moments = {{0.0}};
      for (int mom = 0; mom < NUM_SOOT_MOMENTS + 1; ++mom) {
        moments[mom] = s_arr[nbx](i, j, k, UFSOOT + mom);
      }
      sd->momConvClipConv(moments.data());
      for (int mom = 0; mom < NUM_SOOT_MOMENTS + 1; ++mom) {
        s_arr[nbx](i, j, k, UFSOOT + mom) = moments[mom];
      }
    });
  amrex::Gpu::synchronize();
}
