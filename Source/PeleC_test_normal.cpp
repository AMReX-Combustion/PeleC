#include <PeleC.H>
#include <PeleC_F.H>
#include <AMReX_ParallelDescriptor.H>

using namespace amrex;


#ifdef PELE_USE_EB
#include <PeleC_init_eb_F.H>
#include <AMReX_MultiCutFab.H>
#endif


void
PeleC::test_dn() {

  amrex::Print() << "Testing normal derivative calculation for cut cells..." << std::endl;

  // Assume we have sv_eb_bndry_geom and sv_eb_bndry_grad_stencil filled. Apply them to field and dump out a list of dphi/dn

 MultiFab S(grids,dmap,NUM_STATE,NUM_GROW,MFInfo(),Factory());


  FArrayBox Qfab, coeff_cc;


  std::ofstream ofs("EB_flux.txt", std::ofstream::out);



  for (MFIter mfi(S, MFItInfo().EnableTiling(hydro_tile_size).SetDynamic(true));
       mfi.isValid(); ++mfi) {

    const Box  vbox = mfi.tilebox();
    int ng = S.nGrow();
    const Box  gbox = amrex::grow(vbox,ng);
    const Box  cbox = amrex::grow(vbox,ng-1);
    const Box& dbox = geom.Domain();

    const int* lo = vbox.loVect();
    const int* hi = vbox.hiVect();

    const EBFArrayBox& Sfab = static_cast<const EBFArrayBox&>(S[mfi]);
    const auto& flag_fab = Sfab.getEBCellFlagFab();
    FabType typ = flag_fab.getType(cbox);


    int local_i = mfi.LocalIndex();
    int Ncut = no_eb_in_domain ? 0 : sv_eb_bndry_grad_stencil[local_i].size();

    int myproc = ParallelDescriptor::MyProc();

    if (typ == FabType::singlevalued && Ncut > 0) {
      sv_eb_flux[local_i].setVal(0);  // Default to Neumann for all fields

      int Nvals = sv_eb_bcval[local_i].numPts();
      int Nflux = sv_eb_flux[local_i].numPts();
      int nComp = 1;


      sv_eb_bcval[local_i].setVal(1.0, 0);

      Box box_to_apply = mfi.growntilebox(2);
      amrex::Print() << "gbox: " << gbox << std::endl;
      amrex::Print() << "box to apply: " << box_to_apply << std::endl;

      Qfab.resize(gbox, 1);
      coeff_cc.resize(gbox, 1);

      Qfab.setVal(1.0);
      coeff_cc.setVal(1.0);

      const Real* dx = geom.CellSize();
      const Real* prob_lo = geom.ProbLo();

      pc_set_synthetic_data(BL_TO_FORTRAN_BOX(box_to_apply), BL_TO_FORTRAN_N_ANYD(Qfab,0),prob_lo, dx);

      pc_apply_eb_boundry_flux_stencil(BL_TO_FORTRAN_BOX(box_to_apply),
                                       sv_eb_bndry_grad_stencil[local_i].data(),
                                       &Ncut,
                                       BL_TO_FORTRAN_N_ANYD(Qfab, 0),
                                       BL_TO_FORTRAN_N_ANYD(coeff_cc, 0),
                                       sv_eb_bcval[local_i].dataPtr(0),
                                       &Nvals,
                                       sv_eb_flux[local_i].dataPtr(0),
                                       &Nflux, &nComp);

      for (int L=0; L<Nflux; ++L) {
        amrex::Print(ofs) << local_i << " " << L << " " << sv_eb_bndry_geom[local_i][L].eb_vfrac << sv_eb_bndry_geom[local_i][L].iv << " " << " "<< sv_eb_flux[local_i](L,0) << std::endl;
      }
    }
  }
  ofs.close();
  amrex::Print() << "Finished testing normal derivative calculation for cut cells..." << std::endl << std::flush;
  amrex::Abort("Done unit test.");

  return;
}
