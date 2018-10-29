#include <PeleC.H>
#include <PeleC_F.H>

using std::string;
using namespace amrex;

#include <Transport_F.H>
#include <Filter.H>

#ifdef _OPENMP
#include <omp.h>
#endif

void
PeleC::construct_old_les_source(amrex::Real time, amrex::Real dt, int sub_iteration, int sub_ncycle)
{
  // Add grow cells necessary for explicit filtering of source terms
  if (use_explicit_filter)
  {
    filtered_les_source.define(grids,dmap, NUM_STATE,old_sources[les_src]->nGrow(), MFInfo(), Factory());
    old_sources[les_src]->define(grids,dmap, NUM_STATE,old_sources[les_src]->nGrow()+nGrowF, MFInfo(), Factory());
  }

  old_sources[les_src]->setVal(0.0);

  Real flux_factor_old = 0.5;
  getLESTerm (time, dt, *old_sources[les_src], flux_factor_old);

  old_sources[les_src]->FillBoundary(geom.periodicity());

}

void
PeleC::construct_new_les_source(amrex::Real time, amrex::Real dt, int sub_iteration, int sub_ncycle)
{
  // Add grow cells necessary for explicit filtering of source terms
  if (use_explicit_filter)
  {
    filtered_les_source.define(grids,dmap, NUM_STATE,new_sources[les_src]->nGrow(), MFInfo(), Factory());
    new_sources[les_src]->define(grids,dmap, NUM_STATE,new_sources[les_src]->nGrow()+nGrowF, MFInfo(), Factory());
  }

  new_sources[les_src]->setVal(0.0);

  Real flux_factor_new = sub_iteration==sub_ncycle-1 ? 0.5 : 0;
  getLESTerm (time, dt, *new_sources[les_src], flux_factor_new);

}


/**
 * Calculate the LES term by calling an SFS model
 * Across all conserved state components, compute the LES "source term"
 *    = -Div(LESFlux).
 **/
void
PeleC::getLESTerm (amrex::Real time, amrex::Real dt, amrex::MultiFab& LESTerm, amrex::Real flux_factor)
{
  BL_PROFILE("PeleC::getLESTerm()");

  if (do_les == 0){
    LESTerm.setVal(0,0,NUM_STATE,LESTerm.nGrow());
    return;
  }

  if (verbose) {
    amrex::Print() << "... Computing LES term at time " << time << std::endl;
  }

  switch(les_model){

  case 0:
    getSmagorinskyLESTerm(time, dt, LESTerm, flux_factor);
    break;

  case 1:
    getDynamicSmagorinskyLESTerm(time, dt, LESTerm, flux_factor);
    break;

  default:
    amrex::Error("Invalid les_model number.");
    break;

  }

  // Extrapolate to ghost cells
  if (LESTerm.nGrow() > 0)
  {
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(LESTerm,hydro_tile_size); mfi.isValid(); ++mfi)
    {
      BL_PROFILE("PeleC::diffextrap calls");

      const Box& vbx = mfi.validbox();
      const Box& tbx = mfi.tilebox();
      pc_diffextrap(BL_TO_FORTRAN_BOX(tbx),
                    BL_TO_FORTRAN_BOX(vbx),
                    BL_TO_FORTRAN_N_3D(LESTerm[mfi], Xmom), &amrex::SpaceDim);

      pc_diffextrap(BL_TO_FORTRAN_BOX(tbx),
                    BL_TO_FORTRAN_BOX(vbx),
                    BL_TO_FORTRAN_N_3D(LESTerm[mfi], FirstSpec), &NumSpec);

      const int one = 1;
      pc_diffextrap(BL_TO_FORTRAN_BOX(tbx),
                    BL_TO_FORTRAN_BOX(vbx),
                    BL_TO_FORTRAN_N_3D(LESTerm[mfi], Eden), &one);
    }
  }

  // Filter the SGS source term
  if (use_explicit_filter)
  {
    les_filter.apply_filter(hydro_tile_size, LESTerm, filtered_les_source);
    LESTerm.define(grids,dmap, NUM_STATE,filtered_les_source.nGrow(), MFInfo(), Factory());
    MultiFab::Copy(LESTerm,
                   filtered_les_source,0,0,NUM_STATE,filtered_les_source.nGrow());
  }
}



/**
 * Calculate the LES term using the Smagorinsky SFS model
 **/
void
PeleC::getSmagorinskyLESTerm (amrex::Real time, amrex::Real dt, amrex::MultiFab& LESTerm, amrex::Real flux_factor)
{

  int dComp_rhoD = 0;
  int dComp_mu = dComp_rhoD + NumSpec;
  int dComp_xi = dComp_mu + 1;
  int dComp_lambda = dComp_xi + 1;
  int nCompTr = dComp_lambda + 1;
  int ngrow = 1;
  const Real* dx = geom.CellSize();

  Real dx1 = dx[0];
  for (int d=1; d<BL_SPACEDIM; ++d) {
    dx1 *= dx[d];
  }
  std::array<Real,BL_SPACEDIM> dxD = {D_DECL(dx1, dx1, dx1)};
  const Real *dxDp = &(dxD[0]);

  MultiFab S(grids, dmap, NUM_STATE, ngrow);
  FillPatch(*this, S, ngrow, time, State_Type, 0, NUM_STATE); // FIXME: time+dt?

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    FArrayBox Qfab, Qaux, Lterm;
    FArrayBox flux_ec[BL_SPACEDIM], tander_ec[BL_SPACEDIM];
    IArrayBox bcMask;


    for (MFIter mfi(S, MFItInfo().EnableTiling(hydro_tile_size).SetDynamic(true)); mfi.isValid(); ++mfi)
    {
      const Box  vbox = mfi.tilebox();
      const Box  gbox = amrex::grow(vbox,ngrow);
      const Box  cbox = amrex::grow(vbox,ngrow-1);
      const Box& dbox = geom.Domain();

#ifdef PELEC_USE_EB
      const EBFArrayBox& Sfab = static_cast<const EBFArrayBox&>(S[mfi]);
      const auto& flag_fab = Sfab.getEBCellFlagFab();
      FabType typ = flag_fab.getType(cbox);
      if (typ != FabType::regular) {
        amrex::Error("LES on a non-regular EB Fab is not available.");
      }
#else
      const FArrayBox& Sfab = S[mfi];
#endif


      Qfab.resize(gbox,QVAR);
      int nqaux = NQAUX > 0 ? NQAUX : 1;
      Qaux.resize(gbox,nqaux);

      { // Get primitives, Q, including (Y, T, p, rho) from conserved state, required for L term
        BL_PROFILE("PeleC::ctoprim call");
        ctoprim(ARLIM_3D(gbox.loVect()), ARLIM_3D(gbox.hiVect()),
                Sfab.dataPtr(), ARLIM_3D(Sfab.loVect()), ARLIM_3D(Sfab.hiVect()),
                Qfab.dataPtr(), ARLIM_3D(Qfab.loVect()), ARLIM_3D(Qfab.hiVect()),
                Qaux.dataPtr(), ARLIM_3D(Qaux.loVect()), ARLIM_3D(Qaux.hiVect()));
      }

      // Container on grown region, required to support hybrid divergence and redistribution
      Lterm.resize(cbox,NUM_STATE);

      // Get the tangential derivatives
      for (int d=0; d<BL_SPACEDIM; ++d)
      {
        Box ebox = amrex::surroundingNodes(cbox,d);
        flux_ec[d].resize(ebox,NUM_STATE);  flux_ec[d].setVal(0);

#if (BL_SPACEDIM > 1)
        int nCompTan = AMREX_D_PICK(1, 2, 6);
        tander_ec[d].resize(ebox,nCompTan); tander_ec[d].setVal(0);
        {
          BL_PROFILE("PeleC::pc_compute_tangential_vel_derivs call");
          pc_compute_tangential_vel_derivs(cbox.loVect(), cbox.hiVect(),
                                           dbox.loVect(), dbox.hiVect(),
                                           BL_TO_FORTRAN_ANYD(Qfab),
                                           BL_TO_FORTRAN_ANYD(tander_ec[d]),
                                           geom.CellSize(),&d);
        }
#endif
      } // loop over dimension

      { // Compute extensive LES fluxes, F.A and (1/Vol).Div(F.A)
        BL_PROFILE("PeleC::pc_smagorinsky_sfs_term()");
        pc_smagorinsky_sfs_term(cbox.loVect(), cbox.hiVect(),
                                dbox.loVect(), dbox.hiVect(),
                                BL_TO_FORTRAN_ANYD(Qfab),
#if (BL_SPACEDIM > 1)
                                BL_TO_FORTRAN_ANYD(tander_ec[0]),
#endif
                                BL_TO_FORTRAN_ANYD(area[0][mfi]),
                                BL_TO_FORTRAN_ANYD(flux_ec[0]),
#if (BL_SPACEDIM > 1)
                                BL_TO_FORTRAN_ANYD(tander_ec[1]),
                                BL_TO_FORTRAN_ANYD(area[1][mfi]),
                                BL_TO_FORTRAN_ANYD(flux_ec[1]),
#if (BL_SPACEDIM > 2)
                                BL_TO_FORTRAN_ANYD(tander_ec[2]),
                                BL_TO_FORTRAN_ANYD(area[2][mfi]),
                                BL_TO_FORTRAN_ANYD(flux_ec[2]),
#endif
#endif
                                BL_TO_FORTRAN_ANYD(volume[mfi]),
                                BL_TO_FORTRAN_ANYD(Lterm),
                                geom.CellSize());
      }

      LESTerm[mfi].setVal(0,vbox,0, NUM_STATE);
      LESTerm[mfi].copy(Lterm,vbox,0,vbox,0,NUM_STATE);

      if (do_reflux && flux_factor != 0)
      {
        for (int d = 0; d < BL_SPACEDIM ; d++)
        {
          flux_ec[d].mult(flux_factor);
        }

        if (level < parent->finestLevel())
        {
          getFluxReg(level+1).CrseAdd(mfi,{D_DECL(&flux_ec[0],&flux_ec[1],&flux_ec[2])}, dxDp, dt, RunOn::Cpu);
        }

        if (level > 0)
        {
          getFluxReg(level).FineAdd(mfi, {D_DECL(&flux_ec[0],&flux_ec[1],&flux_ec[2])}, dxDp,dt, RunOn::Cpu);
        }
      }

    }  // End of MFIter scope
  } // End of OMP scope
}


/**
 * Calculate the LES term using the dynamic Smagorinsky SFS model
 **/
void
PeleC::getDynamicSmagorinskyLESTerm (amrex::Real time, amrex::Real dt, amrex::MultiFab& LESTerm, amrex::Real flux_factor)
{

  /*
    Note on the grow cells:
    N                                 (cbox) |----->         LESTerm
    N                                 (cbox) |----->         coeff_ec
    N + 1                            (g4box) |------>        filtered_coeff_cc
    N + 1 + nGrowC                   (g3box) |-------->      coeff_cc, filtered(K, RUT, alphaij, alpha, flux_T)
    N + 1 + nGrowC + nGrowD          (g2box) |---------->    filtered(S, Q, Qaux)
    N + 1 + nGrowC + nGrowT          (g1box) |----------->   K, RUT, alphaij, alpha, flux_T
    N + 1 + nGrowC + nGrowT + nGrowD (g0box) |-------------> S, Q, Qaux

    where
    nGrowD = number of grow cells necessary for the diffusion operator
    nGrowC = number of grow cells necessary for filtering the Smagorinsky coefficients
    nGrowT = number of grow cells necessary for filtering the derived quantities (test level)
  */
  Filter test_filter = Filter(les_test_filter_type, les_test_filter_fgr);
  Filter coeff_filter = Filter(box,6);

  const int nGrowD = 1;
  const int nGrowC = coeff_filter.get_filter_ngrow();
  const int nGrowT = test_filter.get_filter_ngrow();

  int dComp_rhoD = 0;
  int dComp_mu = dComp_rhoD + NumSpec;
  int dComp_xi = dComp_mu + 1;
  int dComp_lambda = dComp_xi + 1;
  int nCompTr = dComp_lambda + 1;
  const Real* dx = geom.CellSize();

  Real dx1 = dx[0];
  for (int d=1; d<BL_SPACEDIM; ++d) {
    dx1 *= dx[d];
  }
  std::array<Real,BL_SPACEDIM> dxD = {D_DECL(dx1, dx1, dx1)};
  const Real *dxDp = &(dxD[0]);

  // 1. Get state variable data
  MultiFab S(grids, dmap, NUM_STATE, nGrowD+nGrowC+nGrowT+1);
  FillPatch(*this, S, nGrowD+nGrowC+nGrowT+1, time, State_Type, 0, NUM_STATE); // FIXME: time+dt?
  LES_Coeffs.setVal(0.0);

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    FArrayBox Qfab, Qaux, Lterm;
    FArrayBox coeff_ec[BL_SPACEDIM], flux_ec[BL_SPACEDIM], tander_ec[BL_SPACEDIM];
    IArrayBox bcMask;
    for (MFIter mfi(S, MFItInfo().EnableTiling(hydro_tile_size).SetDynamic(true)); mfi.isValid(); ++mfi)
    {

      const Box  vbox = mfi.tilebox();
      const Box  g0box = amrex::grow(vbox,nGrowD+nGrowC+nGrowT+1);
      const Box  g1box = amrex::grow(vbox,nGrowC+nGrowT+1);
      const Box  g2box = amrex::grow(vbox,nGrowD+nGrowC+1);
      const Box  g3box = amrex::grow(vbox,nGrowC+1);
      const Box  g4box = amrex::grow(vbox,1);
      const Box  cbox = amrex::grow(vbox,0);
      const Box& dbox = geom.Domain();

#ifdef PELEC_USE_EB
      const EBFArrayBox& Sfab = static_cast<const EBFArrayBox&>(S[mfi]);
      const auto& flag_fab = Sfab.getEBCellFlagFab();
      FabType typ = flag_fab.getType(cbox);
      if (typ != FabType::regular) {
        amrex::Error("LES on a non-regular EB Fab is not available.");
      }
#else
      const FArrayBox& Sfab = S[mfi];
#endif

      Qfab.resize(g0box,QVAR);
      int nqaux = NQAUX > 0 ? NQAUX : 1;
      Qaux.resize(g0box,nqaux);

      { // Get primitives, Q, including (Y, T, p, rho) from conserved state, required for L term
        BL_PROFILE("PeleC::ctoprim call");
        ctoprim(ARLIM_3D(g0box.loVect()), ARLIM_3D(g0box.hiVect()),
                Sfab.dataPtr(), ARLIM_3D(Sfab.loVect()), ARLIM_3D(Sfab.hiVect()),
                Qfab.dataPtr(), ARLIM_3D(Qfab.loVect()), ARLIM_3D(Qfab.hiVect()),
                Qaux.dataPtr(), ARLIM_3D(Qaux.loVect()), ARLIM_3D(Qaux.hiVect()));
      }

      // Container on grown region, required to support hybrid divergence and redistribution
      Lterm.resize(cbox,NUM_STATE);

      // Get the tangential derivatives
      for (int d=0; d<BL_SPACEDIM; ++d)
      {
        Box ebox = amrex::surroundingNodes(g1box,d);

#if (BL_SPACEDIM > 1)
        int nCompTan = AMREX_D_PICK(1, 2, 6);
        tander_ec[d].resize(ebox,nCompTan); tander_ec[d].setVal(0);
        {
          BL_PROFILE("PeleC::pc_compute_tangential_vel_derivs call");
          pc_compute_tangential_vel_derivs(cbox.loVect(), cbox.hiVect(),
                                           dbox.loVect(), dbox.hiVect(),
                                           BL_TO_FORTRAN_ANYD(Qfab),
                                           BL_TO_FORTRAN_ANYD(tander_ec[d]),
                                           geom.CellSize(),&d);
        }
#endif
      } // loop over dimension


      // 2. Get dynamic Smagorinsky derived quantities after setting the
      // BC. These quantities need to be stored because we need to filter
      // them at the test filter level
      const int upper_triangle_n = static_cast<int>(0.5*BL_SPACEDIM*(BL_SPACEDIM+1));
      FArrayBox K, RUT, alphaij, alpha, flux_T;
      K.resize(g1box, upper_triangle_n);
      RUT.resize(g1box, BL_SPACEDIM);
      alphaij.resize(g1box, BL_SPACEDIM*BL_SPACEDIM);
      alpha.resize(g1box, BL_SPACEDIM);
      flux_T.resize(g1box, BL_SPACEDIM);

      {
        BL_PROFILE("PeleC::pc_smagorinsky_sfs_term()");
        pc_dynamic_smagorinsky_quantities(g1box.loVect(), g1box.hiVect(),
                                          dbox.loVect(), dbox.hiVect(),
                                          BL_TO_FORTRAN_ANYD(Qfab),
#if (BL_SPACEDIM > 1)
                                          BL_TO_FORTRAN_ANYD(tander_ec[0]),
                                          BL_TO_FORTRAN_ANYD(tander_ec[1]),
#endif
#if (BL_SPACEDIM > 2)
                                          BL_TO_FORTRAN_ANYD(tander_ec[2]),
#endif
                                          BL_TO_FORTRAN_ANYD(K),
                                          BL_TO_FORTRAN_ANYD(RUT),
                                          BL_TO_FORTRAN_ANYD(alphaij),
                                          BL_TO_FORTRAN_ANYD(alpha),
                                          BL_TO_FORTRAN_ANYD(flux_T),
                                          &les_filter_fgr,
                                          geom.CellSize());
      }



      // 3. Filter the state variables and the derived quantities at the
      // test filter level
      FArrayBox filtered_S, filtered_Q, filtered_Qaux, filtered_K, filtered_RUT, filtered_alphaij, filtered_alpha, filtered_flux_T;
      filtered_S.resize(g2box,NUM_STATE);
      filtered_Q.resize(g2box,QVAR);
      filtered_Qaux.resize(g2box,NQAUX>0?NQAUX:1);
      filtered_K.resize(g3box,upper_triangle_n);
      filtered_RUT.resize(g3box,BL_SPACEDIM);
      filtered_alphaij.resize(g3box,BL_SPACEDIM*BL_SPACEDIM);
      filtered_alpha.resize(g3box,BL_SPACEDIM);
      filtered_flux_T.resize(g3box,BL_SPACEDIM);

      test_filter.apply_filter(g2box, Sfab, filtered_S);
      ctoprim(ARLIM_3D(g2box.loVect()), ARLIM_3D(g2box.hiVect()),
              filtered_S.dataPtr(), ARLIM_3D(filtered_S.loVect()), ARLIM_3D(filtered_S.hiVect()),
              filtered_Q.dataPtr(), ARLIM_3D(filtered_Q.loVect()), ARLIM_3D(filtered_Q.hiVect()),
              filtered_Qaux.dataPtr(), ARLIM_3D(filtered_Qaux.loVect()), ARLIM_3D(filtered_Qaux.hiVect()));
      test_filter.apply_filter(g3box, K, filtered_K);
      test_filter.apply_filter(g3box, RUT, filtered_RUT);
      test_filter.apply_filter(g3box, alphaij, filtered_alphaij);
      test_filter.apply_filter(g3box, alpha, filtered_alpha);
      test_filter.apply_filter(g3box, flux_T, filtered_flux_T);

      // 4. Calculate the dynamic Smagorinsky coefficients
      int do_harmonic = 1;
      FArrayBox coeff_cc;
      coeff_cc.resize(g3box, nCompC);

      {
        BL_PROFILE("PeleC::pc_dynamic_smagorinsky_coeffs()");
        pc_dynamic_smagorinsky_coeffs(g3box.loVect(), g3box.hiVect(),
                                      dbox.loVect(), dbox.hiVect(),
                                      BL_TO_FORTRAN_ANYD(filtered_Q),
#if (BL_SPACEDIM > 1)
                                      BL_TO_FORTRAN_ANYD(tander_ec[0]),
                                      BL_TO_FORTRAN_ANYD(tander_ec[1]),
#endif
#if (BL_SPACEDIM > 2)
                                      BL_TO_FORTRAN_ANYD(tander_ec[2]),
#endif
                                      BL_TO_FORTRAN_ANYD(filtered_K),
                                      BL_TO_FORTRAN_ANYD(filtered_RUT),
                                      BL_TO_FORTRAN_ANYD(filtered_alphaij),
                                      BL_TO_FORTRAN_ANYD(filtered_alpha),
                                      BL_TO_FORTRAN_ANYD(filtered_flux_T),
                                      BL_TO_FORTRAN_N_ANYD(coeff_cc,comp_Cs2),
                                      BL_TO_FORTRAN_N_ANYD(coeff_cc,comp_CI),
                                      BL_TO_FORTRAN_N_ANYD(coeff_cc,comp_PrT),
                                      &les_test_filter_fgr,
                                      geom.CellSize());
      }

      // FIXME: REMOVE THIS LATER
      coeff_cc.setVal(0.16*0.16,g3box,comp_Cs2,1);
      coeff_cc.setVal(0.09,g3box,comp_CI,1);
      coeff_cc.setVal(0.7,g3box,comp_PrT,1);

      // 5. Filter to smooth the dynamic coefficients
      coeff_filter.apply_filter(g4box, coeff_cc, LES_Coeffs[mfi]);

      // 6. Get the SFS term
      for (int d=0; d<BL_SPACEDIM; ++d) {
        Box ebox = amrex::surroundingNodes(cbox,d);
        coeff_ec[d].resize(ebox,nCompC);
        flux_ec[d].resize(ebox,NUM_STATE); flux_ec[d].setVal(0);
        pc_move_transport_coeffs_to_ec(ARLIM_3D(cbox.loVect()), ARLIM_3D(cbox.hiVect()),
                                       ARLIM_3D(dbox.loVect()), ARLIM_3D(dbox.hiVect()),
                                       BL_TO_FORTRAN_ANYD(LES_Coeffs[mfi]),
                                       BL_TO_FORTRAN_ANYD(coeff_ec[d]),
                                       &d, &nCompC, &do_harmonic);
      }

      {
        BL_PROFILE("PeleC::pc_dynamic_smagorinsky_sfs_term()");
        pc_dynamic_smagorinsky_sfs_term(cbox.loVect(), cbox.hiVect(),
                                        BL_TO_FORTRAN_ANYD(Qfab),
                                        BL_TO_FORTRAN_ANYD(alphaij),
                                        BL_TO_FORTRAN_ANYD(alpha),
                                        BL_TO_FORTRAN_ANYD(flux_T),
                                        BL_TO_FORTRAN_N_ANYD(coeff_ec[0],comp_Cs2),
                                        BL_TO_FORTRAN_N_ANYD(coeff_ec[0],comp_CI),
                                        BL_TO_FORTRAN_N_ANYD(coeff_ec[0],comp_PrT),
                                        BL_TO_FORTRAN_ANYD(area[0][mfi]),
                                        BL_TO_FORTRAN_ANYD(flux_ec[0]),
#if (BL_SPACEDIM > 1)
                                        BL_TO_FORTRAN_N_ANYD(coeff_ec[1],comp_Cs2),
                                        BL_TO_FORTRAN_N_ANYD(coeff_ec[1],comp_CI),
                                        BL_TO_FORTRAN_N_ANYD(coeff_ec[1],comp_PrT),
                                        BL_TO_FORTRAN_ANYD(area[1][mfi]),
                                        BL_TO_FORTRAN_ANYD(flux_ec[1]),
#if (BL_SPACEDIM > 2)
                                        BL_TO_FORTRAN_N_ANYD(coeff_ec[2],comp_Cs2),
                                        BL_TO_FORTRAN_N_ANYD(coeff_ec[2],comp_CI),
                                        BL_TO_FORTRAN_N_ANYD(coeff_ec[2],comp_PrT),
                                        BL_TO_FORTRAN_ANYD(area[2][mfi]),
                                        BL_TO_FORTRAN_ANYD(flux_ec[2]),
#endif
#endif
                                        BL_TO_FORTRAN_ANYD(volume[mfi]),
                                        BL_TO_FORTRAN_ANYD(Lterm),
                                        geom.CellSize());
      }

      LESTerm[mfi].setVal(0,vbox,0, NUM_STATE);
      LESTerm[mfi].copy(Lterm,vbox,0,vbox,0,NUM_STATE);

      if (do_reflux && flux_factor != 0)
      {
        for (int d = 0; d < BL_SPACEDIM ; d++)
        {
          flux_ec[d].mult(flux_factor);
        }

        if (level < parent->finestLevel())
        {
          getFluxReg(level+1).CrseAdd(mfi,{D_DECL(&flux_ec[0],&flux_ec[1],&flux_ec[2])}, dxDp, dt, RunOn::Cpu);
        }

        if (level > 0)
        {
          getFluxReg(level).FineAdd(mfi, {D_DECL(&flux_ec[0],&flux_ec[1],&flux_ec[2])}, dxDp, dt, RunOn::Cpu);
        }
      }

    }  // End of MFIter scope
  } // End of OMP scope
}
