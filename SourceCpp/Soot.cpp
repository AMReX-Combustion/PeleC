#include "PeleC.H"

#include <Transport.H>

void
PeleC::construct_old_soot_source(amrex::Real time,
                                 amrex::Real dt)
{
  old_sources[soot_src]->setVal(0.0);
  if (!add_soot_src)
    return;
  amrex::MultiFab& S_old = get_old_data(State_Type);

  int ng = 0; // None filled

  fill_soot_source(time, dt, S_old, *old_sources[soot_src], ng);

  old_sources[soot_src]->FillBoundary(geom.periodicity());
}

void
PeleC::construct_new_soot_source(amrex::Real time,
                                 amrex::Real dt)
{
  new_sources[soot_src]->setVal(0.0);
  if (!add_soot_src)
    return;
  amrex::MultiFab& S_new = get_new_data(State_Type);

  int ng = 0;

  fill_soot_source(time, dt, S_new, *new_sources[soot_src], ng);
}

void
PeleC::fill_soot_source (amrex::Real      time,
                         amrex::Real      dt,
                         amrex::MultiFab& state,
                         amrex::MultiFab& soot_src,
                         int              ng)
{
  const int nCompTr = dComp_lambda + 1;
  const amrex::Real* dx = geom.CellSize();
  const amrex::Real* prob_lo = geom.ProbLo();

#ifdef PELE_USE_EB
  auto const& fact = dynamic_cast<EBFArrayBoxFactory const&>(state.Factory());
  auto const& flags = fact.getMultiEBCellFlagFab();
#endif

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(soot_src, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    const amrex::Box& bx = mfi.growntilebox(ng);
    amrex::RealBox gridloc = amrex::RealBox(grids[mfi.index()],geom.CellSize(),geom.ProbLo());
#ifdef PELE_USE_EB
    const auto& flag_fab = flags[mfi];
    FabType typ = flag_fab.getType(bx);
    if (typ == FabType::covered) {
      continue;
    }
#endif
    amrex::FArrayBox& Sfab = state[mfi];
    amrex::FArrayBox& soot_fab = soot_src[mfi];
    auto const& s_arr = Sfab.array();
    const int nqaux = NQAUX > 0 ? NQAUX : 1;
    amrex::FArrayBox coeff_cc(bx, nCompTr), q(bx,QVAR), qaux(bx, nqaux);
    amrex::Elixir qeli = q.elixir();
    amrex::Elixir qauxeli = qaux.elixir();
    amrex::Elixir coeff_eli = coeff_cc.elixir();
    auto const& q_arr = q.array();
    auto const& qaux_arr = qaux.array();
    auto const& coeff_arr = coeff_cc.array();

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
      auto const& coe_rhoD = coeff_cc.array(dComp_rhoD);
      auto const& coe_mu = coeff_cc.array(dComp_mu);
      auto const& coe_xi = coeff_cc.array(dComp_xi);
      auto const& coe_lambda = coeff_cc.array(dComp_lambda);
      BL_PROFILE("PeleC::get_transport_coeffs()");
      amrex::launch(bx, [=] AMREX_GPU_DEVICE(amrex::Box const& tbx) {
        get_transport_coeffs(
          tbx, qar_yin, qar_Tin, qar_rhoin, coe_rhoD, coe_mu, coe_xi,
          coe_lambda);
        });
    }
    auto const& soot_arr = soot_fab.array();
    soot_model->addSootSourceTerm(bx, q_arr, coeff_arr, soot_arr, time, dt);
  }
  amrex::ReduceOps<amrex::ReduceOpMin> reduce_op;
  amrex::ReduceData<amrex::Real> reduce_data(reduce_op);
  using ReduceTuple = typename decltype(reduce_data)::Type;
  amrex::Real soot_dt = std::numeric_limits<amrex::Real>::max();
  SootData sd = soot_model->getSootData();
  for (MFIter mfi(soot_src, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    const amrex::Box& bx = mfi.growntilebox(ng);
    amrex::FArrayBox& Sfab = state[mfi];
    const amrex::FArrayBox& soot_fab = soot_src[mfi];
    auto const& s_arr = Sfab.array();
    auto const& soot_arr = soot_fab.array();
    {
      BL_PROFILE("PeleC::retrieveSootTimeStep()");
      reduce_op.eval(bx, reduce_data,
                     [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
                     {
                       amrex::Real rho_rate = soot_arr(i, j, k, URHO)/s_arr(i, j, k, URHO);
                       amrex::Real eng_rate = soot_arr(i ,j, k, UEDEN)/s_arr(i, j, k, UEDEN);
                       amrex::Real max_rate = amrex::max(std::abs(rho_rate), std::abs(eng_rate));
                       for (int n = 0; n != NUM_SPECIES; ++n) {
                         int indx = UFS + n;
                         amrex::Real cur_rate =
                           soot_arr(i, j, k, indx)/(s_arr(i, j, k, indx) + 1.E-12);
                         max_rate = amrex::max(max_rate, std::abs(cur_rate));
                       }
                       amrex::Real moments[NUM_SOOT_VARS];
                       for (int n = 0; n != NUM_SOOT_VARS; ++n)
                         moments[n] = s_arr(i, j, k, UFSOOT + n);
                       // Convert moments to mol, clip moments, and convert back to CGS
                       sd.momConvClipConv(moments);
                       // Limit time step based only on positive moment sources
                       for (int n = 0; n != NUM_SOOT_VARS; ++n) {
                         int indx = UFSOOT + n;
                         s_arr(i, j, k, indx) = moments[n];
                         max_rate =
                           amrex::max(max_rate, -soot_arr(i, j, k, indx)/moments[n]);
                       }
                       return 1./max_rate;
                     });
      ReduceTuple hv = reduce_data.value();
      Real ldt_cpu = amrex::get<0>(hv);
      soot_dt = amrex::min(soot_dt, ldt_cpu);
    }
  }
  soot_model->setTimeStep(soot_dt);
}
