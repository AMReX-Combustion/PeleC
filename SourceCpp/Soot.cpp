#include "PeleC.H"

#include <Transport.H>

void
PeleC::construct_old_soot_source(amrex::Real time,
                                 amrex::Real dt)
{
  amrex::MultiFab& S_old = get_old_data(State_Type);

  int ng = 0; // None filled

  old_sources[soot_src]->setVal(0.0);
  fill_soot_source(time, dt, S_old, S_old, *old_sources[soot_src], ng);

  old_sources[soot_src]->FillBoundary(geom.periodicity());
}

void
PeleC::construct_new_soot_source(amrex::Real time,
                                 amrex::Real dt)
{
  amrex::MultiFab& S_old = get_old_data(State_Type);
  amrex::MultiFab& S_new = get_new_data(State_Type);

  int ng = 0;

  new_sources[soot_src]->setVal(0.0);

  fill_soot_source(time, dt, S_old, S_new, *new_sources[soot_src], ng);

}

void
PeleC::fill_soot_source (amrex::Real            time,
                         amrex::Real            dt,
                         const amrex::MultiFab& state_old,
                         const amrex::MultiFab& state_new,
                         amrex::MultiFab&       soot_src,
                         int                    ng)
{
  const int nCompTr = dComp_lambda + 1;
  const amrex::Real* dx = geom.CellSize();
  const amrex::Real* prob_lo = geom.ProbLo();

#ifdef PELE_USE_EB
  auto const& fact = dynamic_cast<EBFArrayBoxFactory const&>(state_old.Factory());
  auto const& flags = fact.getMultiEBCellFlagFab();
#endif

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(soot_src, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
  {
    const amrex::Box& bx = mfi.growntilebox(ng);
    amrex::RealBox gridloc = amrex::RealBox(grids[mfi.index()],geom.CellSize(),geom.ProbLo());

#ifdef PELE_USE_EB
    const auto& flag_fab = flags[mfi];
    FabType typ = flag_fab.getType(bx);
    if (typ == FabType::covered) {
      continue;
    }
#endif

    auto const& s_old_arr = state_old.array(mfi);
    auto const& s_new_arr = state_new.array(mfi);
    auto& soot_arr = soot_src.array(mfi);
    const int nqaux = NQAUX > 0 ? NQAUX : 1;
    amrex::FArrayBox coeff_cc(bx, nCompTr), q(bx,QVAR), qaux(bx, nqaux);
    amrex::Elixir qeli = q.elixir();
    amrex::Elixir qauxeli = qaux.elixir();
    amrex::Elixir coef_eli = coef_cc.elixir();
    auto const& q_arr = q.array();
    auto const& qaux_arr = qaux.array();
    auto const& coef_arr = coef_cc.array();

    // Get primitives, Q, including (Y, T, p, rho) from conserved state
    // required for D term
    {
      BL_PROFILE("PeleC::ctoprim()");
      amrex::ParallelFor(
        bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          pc_ctoprim(i, j, k, s_new_arr, q_arr, qaux_arr);
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
    soot_model->addSootSourceTerm(bx, Snfab, Qfab, coeff_cc, Ffab, time, dt);
  }
}
