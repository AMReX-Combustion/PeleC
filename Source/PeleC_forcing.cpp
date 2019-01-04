#include "PeleC.H"
#include "PeleC_F.H"

using namespace amrex;

void
PeleC::construct_old_forcing_source(Real time,
                                    Real dt)
{
  MultiFab& S_old = get_old_data(State_Type);

  int ng = 0; // None filled

  old_sources[forcing_src]->setVal(0.0);

  if (!add_forcing_src) return;

  fill_forcing_source(time, dt, S_old, S_old, *old_sources[forcing_src], ng);

  old_sources[forcing_src]->FillBoundary(geom.periodicity());
}

void
PeleC::construct_new_forcing_source(Real time,
                                    Real dt)
{
  MultiFab& S_old = get_old_data(State_Type);
  MultiFab& S_new = get_new_data(State_Type);

  int ng = 0;

  new_sources[forcing_src]->setVal(0.0);

  if (!add_forcing_src) return;

  fill_forcing_source(time, dt, S_old, S_new, *new_sources[forcing_src], ng);

}



void
PeleC::fill_forcing_source (Real time, Real dt,
                            const MultiFab& state_old,
                            const MultiFab& state_new,
                            MultiFab& forcing_src, int ng)
{
  const Real* dx = geom.CellSize();
  const Real* prob_lo = geom.ProbLo();

#ifdef PELE_USE_EB
  auto const& fact = dynamic_cast<EBFArrayBoxFactory const&>(state_old.Factory());
  auto const& flags = fact.getMultiEBCellFlagFab();
#endif

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(forcing_src,true); mfi.isValid(); ++mfi)
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
    auto& Ffab = forcing_src[mfi];

    pc_forcing_src(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
                   BL_TO_FORTRAN_3D(Sofab),
                   BL_TO_FORTRAN_3D(Snfab),
                   BL_TO_FORTRAN_3D(Ffab),
                   ZFILL(prob_lo),ZFILL(dx),
                   ZFILL(gridloc.lo()), ZFILL(gridloc.hi()),
                   &time,&dt);
  }
}
