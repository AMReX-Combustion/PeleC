#include "PeleC.H"
#include "PeleC_F.H"

#include <Transport_F.H>
using namespace amrex;

void
PeleC::construct_old_soot_source(Real time,
				 Real dt)
{
  MultiFab& S_old = get_old_data(State_Type);

  int ng = 0; // None filled

  old_sources[soot_src]->setVal(0.0);
  fill_soot_source(time, dt, S_old, S_old, *old_sources[soot_src], ng);

  old_sources[soot_src]->FillBoundary(geom.periodicity());
}

void
PeleC::construct_new_soot_source(Real time,
				 Real dt)
{
  MultiFab& S_old = get_old_data(State_Type);
  MultiFab& S_new = get_new_data(State_Type);

  int ng = 0;

  new_sources[soot_src]->setVal(0.0);

  fill_soot_source(time, dt, S_old, S_new, *new_sources[soot_src], ng);

}

void
PeleC::fill_soot_source (Real time, Real dt,
			 const MultiFab& state_old,
			 const MultiFab& state_new,
			 MultiFab& soot_src, int ng)
{
  int dComp_rhoD = 0;
  int dComp_mu = dComp_rhoD + NumSpec;
  int dComp_xi = dComp_mu + 1;
  int dComp_lambda = dComp_xi + 1;
  int nCompTr = dComp_lambda + 1;
  const Real* dx = geom.CellSize();
  const Real* prob_lo = geom.ProbLo();

#ifdef PELE_USE_EB
  auto const& fact = dynamic_cast<EBFArrayBoxFactory const&>(state_old.Factory());
  auto const& flags = fact.getMultiEBCellFlagFab();
#endif

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(soot_src,true); mfi.isValid(); ++mfi)
  {
    const Box& bx = mfi.growntilebox(ng);
    RealBox gridloc = RealBox(grids[mfi.index()],geom.CellSize(),geom.ProbLo());

#ifdef PELE_USE_EB
    const auto& flag_fab = flags[mfi];
    FabType typ = flag_fab.getType(bx);
    if (typ == FabType::covered) {
      continue;
    }
#endif

    const auto& Sofab = state_old[mfi];
    const auto& Snfab = state_new[mfi];
    auto& Ffab = soot_src[mfi];
    FArrayBox coeff_cc, Qfab, Qaux;
    Qfab.resize(bx, QVAR);
    int nqaux = NQAUX > 0 ? NQAUX : 1;
    Qaux.resize(bx, nqaux);
    coeff_cc.resize(bx, nCompTr);
    // Get primitives, Q, including (Y, T, p, rho) from conserved state
    // required for D term
    {
      BL_PROFILE("PeleC::ctoprim call");
      ctoprim(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
	      Snfab.dataPtr(), ARLIM_3D(Snfab.loVect()), ARLIM_3D(Snfab.hiVect()),
	      Qfab.dataPtr(), ARLIM_3D(Qfab.loVect()), ARLIM_3D(Qfab.hiVect()),
	      Qaux.dataPtr(), ARLIM_3D(Qaux.loVect()), ARLIM_3D(Qaux.hiVect()));
    }
    // Compute transport coefficients, coincident with Q
    {
      BL_PROFILE("PeleC::get_transport_coeffs call");
      get_transport_coeffs(ARLIM_3D(bx.loVect()),
			   ARLIM_3D(bx.hiVect()),
			   BL_TO_FORTRAN_N_3D(Qfab, cQFS),
			   BL_TO_FORTRAN_N_3D(Qfab, cQTEMP),
			   BL_TO_FORTRAN_N_3D(Qfab, cQRHO),
			   BL_TO_FORTRAN_N_3D(coeff_cc, dComp_rhoD),
			   BL_TO_FORTRAN_N_3D(coeff_cc, dComp_mu),
			   BL_TO_FORTRAN_N_3D(coeff_cc, dComp_xi),
			   BL_TO_FORTRAN_N_3D(coeff_cc, dComp_lambda));
    }
    soot_model->addSootSourceTerm(bx, Snfab, Qfab, coeff_cc, Ffab, time, dt);
  }
}
