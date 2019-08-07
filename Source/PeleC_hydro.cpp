#include "PeleC.H"
#include "PeleC_F.H"

using namespace amrex;

/**
 *  Set up the source terms to go into the hydro.
 */
void
PeleC::construct_hydro_source(const MultiFab& S, Real time, Real dt, int amr_iteration, int amr_ncycle, int sub_iteration, int sub_ncycle)
{
#ifdef PELEC_USE_MOL

    if (verbose && ParallelDescriptor::IOProcessor()) {
        std::cout << "... Zeroing Godunov-based hydro advance" << std::endl;
    }
    hydro_source.setVal(0);

#else 

    if (verbose && ParallelDescriptor::IOProcessor()) {
        std::cout << "... Computing hydro advance" << std::endl;
    }

    BL_ASSERT(S.nGrow() == NUM_GROW+nGrowF);

    sources_for_hydro.setVal(0.0);

    int ng = 0; // TODO: This is currently the largest ngrow of the source data...maybe this needs fixing?

    for (int n = 0; n < src_list.size(); ++n)
    {
      MultiFab::Saxpy(sources_for_hydro,0.5,*new_sources[src_list[n]],0,0,NUM_STATE,ng);
      MultiFab::Saxpy(sources_for_hydro,0.5,*old_sources[src_list[n]],0,0,NUM_STATE,ng);
    }

#ifdef REACTIONS
    // Add I_R terms to advective forcing
    if (do_react == 1) {

      MultiFab::Add(sources_for_hydro, get_new_data(Reactions_Type), 0, FirstSpec, NumSpec, ng);
      MultiFab::Add(sources_for_hydro, get_new_data(Reactions_Type), NumSpec, Eden, 1, ng);

    }
#endif
    
    sources_for_hydro.FillBoundary(geom.periodicity());

    hydro_source.setVal(0);

    int finest_level = parent->finestLevel();

    const Real *dx = geom.CellSize();

    Real dx1 = dx[0];
    for (int d=1; d<BL_SPACEDIM; ++d) {
        dx1 *= dx[d];
    }

    std::array<Real,BL_SPACEDIM> dxD = {D_DECL(dx1, dx1, dx1)};
    const Real *dxDp = &(dxD[0]);

    Real courno    = -1.0e+200;

    MultiFab& S_new = get_new_data(State_Type);

    // note: the radiation consup currently does not fill these
    Real E_added_flux    = 0.;
    Real mass_added_flux = 0.;
    Real xmom_added_flux = 0.;
    Real ymom_added_flux = 0.;
    Real zmom_added_flux = 0.;
    Real mass_lost       = 0.;
    Real xmom_lost       = 0.;
    Real ymom_lost       = 0.;
    Real zmom_lost       = 0.;
    Real eden_lost       = 0.;
    Real xang_lost       = 0.;
    Real yang_lost       = 0.;
    Real zang_lost       = 0.;

    BL_PROFILE_VAR("PeleC::advance_hydro_pc_umdrv()", PC_UMDRV);

#ifdef _OPENMP
#pragma omp parallel reduction(+:E_added_flux,mass_added_flux)		\
    reduction(+:xmom_added_flux,ymom_added_flux,zmom_added_flux)	\
    reduction(+:mass_lost,xmom_lost,ymom_lost,zmom_lost)		\
    reduction(+:eden_lost,xang_lost,yang_lost,zang_lost) 		\
    reduction(max:courno)
#endif
    {
	FArrayBox flux[BL_SPACEDIM];

        FArrayBox filtered_flux[BL_SPACEDIM];
	FArrayBox filtered_source_out;

	FArrayBox pradial(Box::TheUnitBox(),1);
	FArrayBox q, qaux, src_q;
	IArrayBox bcMask[BL_SPACEDIM];

	Real cflLoc = -1.0e+200;
	int is_finest_level = (level == finest_level) ? 1 : 0;
  int flag_nscbc_isAnyPerio = (geom.isAnyPeriodic()) ? 1 : 0; 
  int flag_nscbc_perio[BL_SPACEDIM]; // For 3D, we will know which corners have a periodicity
  for (int d=0; d<BL_SPACEDIM; ++d) {
        flag_nscbc_perio[d] = (DefaultGeometry().isPeriodic(d)) ? 1 : 0;
    }
	const int*  domain_lo = geom.Domain().loVect();
	const int*  domain_hi = geom.Domain().hiVect();
  
	for (MFIter mfi(S_new,hydro_tile_size); mfi.isValid(); ++mfi)
	{
	    const Box& bx    = mfi.tilebox();
	    const Box& qbx = amrex::grow(bx, NUM_GROW+nGrowF);

	    const Box& fbx = amrex::grow(bx, nGrowF);

	    const int* lo = fbx.loVect();
	    const int* hi = fbx.hiVect();

	    const FArrayBox &statein  = S[mfi];
	    FArrayBox &stateout = S_new[mfi];

	    FArrayBox &source_in  = sources_for_hydro[mfi];
	    FArrayBox &source_out = hydro_source[mfi];

	    q.resize(qbx, QVAR);
	    qaux.resize(qbx, NQAUX);
	    src_q.resize(qbx, QVAR);
      
	    ctoprim(ARLIM_3D(qbx.loVect()), ARLIM_3D(qbx.hiVect()),
		    statein.dataPtr(), ARLIM_3D(statein.loVect()), ARLIM_3D(statein.hiVect()),
		    q.dataPtr(), ARLIM_3D(q.loVect()), ARLIM_3D(q.hiVect()),
		    qaux.dataPtr(), ARLIM_3D(qaux.loVect()), ARLIM_3D(qaux.hiVect()));
      
      // Imposing Ghost-Cells Navier-Stokes Characteristic BCs if "UserBC" are used
      // For the theory, see Motheau et al. AIAA J. Vol. 55, No. 10 : pp. 3399-3408, 2017. 
      //
      // The user should provide a bcnormal routine in bc_fill_module with additional optional arguments
      // to temporary fill ghost-cells for EXT_DIR and to provide target BC values.
      // See the examples.
      
      // Allocate fabs for bcMask. Note that we grow in the opposite direction
      // because the Riemann solver wants a face value in a ghost-cell
      for (int i = 0; i < BL_SPACEDIM ; i++)  {
        const Box& bxtmp = amrex::surroundingNodes(fbx,i);
        Box TestBox(bxtmp);
        for(int d=0; d<BL_SPACEDIM; ++d) {
          if (i!=d) TestBox.grow(d,1);
        }
        bcMask[i].resize(TestBox,1);
        bcMask[i].setVal(0);
      }
      
      // Becase bcMask is read in the Riemann solver in any case,
      // here we put physbc values in the appropriate faces for the non-nscbc case
      set_bc_mask(lo, hi, domain_lo, domain_hi,
                  D_DECL(BL_TO_FORTRAN(bcMask[0]),
	                       BL_TO_FORTRAN(bcMask[1]),
                         BL_TO_FORTRAN(bcMask[2])));

      if (nscbc_adv == 1)
      {
        impose_NSCBC(lo, hi, domain_lo, domain_hi,
                     BL_TO_FORTRAN(statein),
                     BL_TO_FORTRAN(q),
                     BL_TO_FORTRAN(qaux),
                     D_DECL(BL_TO_FORTRAN(bcMask[0]),
	                          BL_TO_FORTRAN(bcMask[1]),
                            BL_TO_FORTRAN(bcMask[2])),
                     &flag_nscbc_isAnyPerio, flag_nscbc_perio, 
                     &time, dx, &dt);
      }

	    srctoprim(ARLIM_3D(qbx.loVect()), ARLIM_3D(qbx.hiVect()),
		      q.dataPtr(), ARLIM_3D(q.loVect()), ARLIM_3D(q.hiVect()),
		      qaux.dataPtr(), ARLIM_3D(qaux.loVect()), ARLIM_3D(qaux.hiVect()),
		      source_in.dataPtr(), ARLIM_3D(source_in.loVect()), ARLIM_3D(source_in.hiVect()),
		      src_q.dataPtr(), ARLIM_3D(src_q.loVect()), ARLIM_3D(src_q.hiVect()));

      // Allocate fabs for fluxes
	    for (int i = 0; i < BL_SPACEDIM ; i++)  {
		const Box& bxtmp = amrex::surroundingNodes(fbx,i);
		flux[i].resize(bxtmp,NUM_STATE);
	    }

	    if (!DefaultGeometry().IsCartesian()) {
		pradial.resize(amrex::surroundingNodes(bx,0),1);
	    }

	    pc_umdrv
		(&is_finest_level, &time,
		 lo, hi, domain_lo, domain_hi,
		 BL_TO_FORTRAN(statein), 
		 BL_TO_FORTRAN(stateout),
		 BL_TO_FORTRAN(q),
		 BL_TO_FORTRAN(qaux),
		 BL_TO_FORTRAN(src_q),
		 BL_TO_FORTRAN(source_out),
		 D_DECL(BL_TO_FORTRAN(bcMask[0]),
	          BL_TO_FORTRAN(bcMask[1]),
            BL_TO_FORTRAN(bcMask[2])),
		 dx, &dt,
		 D_DECL(BL_TO_FORTRAN(flux[0]),
			BL_TO_FORTRAN(flux[1]),
			BL_TO_FORTRAN(flux[2])),
#if (BL_SPACEDIM < 3)
		 BL_TO_FORTRAN(pradial),
#endif
		 D_DECL(BL_TO_FORTRAN(area[0][mfi]),
			BL_TO_FORTRAN(area[1][mfi]),
			BL_TO_FORTRAN(area[2][mfi])),
#if (BL_SPACEDIM < 3)
		 BL_TO_FORTRAN(dLogArea[0][mfi]),
#endif
		 BL_TO_FORTRAN(volume[mfi]),
		 &cflLoc, verbose,
		 mass_added_flux,
		 xmom_added_flux,
		 ymom_added_flux,
		 zmom_added_flux,
		 E_added_flux,
		 mass_lost, xmom_lost, ymom_lost, zmom_lost,
		 eden_lost, xang_lost, yang_lost, zang_lost);

	    courno = std::max(courno,cflLoc);

            // Filter hydro source and fluxes here
            if (use_explicit_filter)
            {
              for (int i = 0; i < BL_SPACEDIM ; i++)  {
	        const Box& bxtmp = amrex::surroundingNodes(bx,i);
	        filtered_flux[i].resize(bxtmp,NUM_STATE);
                les_filter.apply_filter(bxtmp, flux[i], filtered_flux[i], Density, NUM_STATE);

                flux[i].setVal(0);
                flux[i].copy(filtered_flux[i], Density, Density, NUM_STATE);
              }

              filtered_source_out.resize(bx, NUM_STATE);
              les_filter.apply_filter(bx, source_out, filtered_source_out, Density, NUM_STATE);

              source_out.setVal(0);
              source_out.copy(filtered_source_out, Density, Density, NUM_STATE);
            }

            if (do_reflux  && sub_iteration == sub_ncycle-1 )
            {
              if (level < finest_level)
              {
                getFluxReg(level+1).CrseAdd(mfi,{D_DECL(&flux[0],&flux[1],&flux[2])}, dxDp,dt,RunOn::Cpu);

                if (!DefaultGeometry().IsCartesian()) {
                    amrex::Abort("Flux registers not r-z compatible yet");
                    //getPresReg(level+1).CrseAdd(mfi,pradial, dx,dt);
                }
              }

              if (level > 0)
              {
                getFluxReg(level).FineAdd(mfi, {D_DECL(&flux[0],&flux[1],&flux[2])}, dxDp,dt, RunOn::Cpu);

                if (!DefaultGeometry().IsCartesian()) {
                    amrex::Abort("Flux registers not r-z compatible yet");
                    //getPresReg(level).FineAdd(mfi,pradial, dx,dt);
                }
              }
            }
	} // MFIter loop

    } // end of OMP parallel region

    hydro_source.FillBoundary(geom.periodicity());
  
    BL_PROFILE_VAR_STOP(PC_UMDRV);

    // Flush Fortran output

    if (verbose)
	flush_output();

    if (track_grid_losses)
    {
	material_lost_through_boundary_temp[0] += mass_lost;
	material_lost_through_boundary_temp[1] += xmom_lost;
	material_lost_through_boundary_temp[2] += ymom_lost;
	material_lost_through_boundary_temp[3] += zmom_lost;
	material_lost_through_boundary_temp[4] += eden_lost;
	material_lost_through_boundary_temp[5] += xang_lost;
	material_lost_through_boundary_temp[6] += yang_lost;
	material_lost_through_boundary_temp[7] += zang_lost;
    }

    if (print_energy_diagnostics)
    {
	Real foo[5] = {E_added_flux, xmom_added_flux, ymom_added_flux,
		       zmom_added_flux, mass_added_flux};

#ifdef BL_LAZY
	Lazy::QueueReduction( [=] () mutable {
#endif
		ParallelDescriptor::ReduceRealSum(foo, 5, ParallelDescriptor::IOProcessorNumber());

		if (ParallelDescriptor::IOProcessor())
		{
		    E_added_flux    = foo[0];
		    xmom_added_flux = foo[1];
		    ymom_added_flux = foo[2];
		    zmom_added_flux = foo[3];
		    mass_added_flux = foo[4];
#ifdef AMREX_DEBUG
		    std::cout << "mass added from fluxes                      : " <<
			mass_added_flux << std::endl;
		    std::cout << "xmom added from fluxes                      : " <<
			xmom_added_flux << std::endl;
		    std::cout << "ymom added from fluxes                      : " <<
			ymom_added_flux << std::endl;
		    std::cout << "zmom added from fluxes                      : " <<
			zmom_added_flux << std::endl;
		    std::cout << "(rho E) added from fluxes                   : " <<
			E_added_flux << std::endl;
#endif
		}
#ifdef BL_LAZY
	    });
#endif
    }

    if (courno > 1.0) {
	std::cout << "WARNING -- EFFECTIVE CFL AT THIS LEVEL " << level << " IS " << courno << '\n';
	if (hard_cfl_limit == 1)
	    amrex::Abort("CFL is too high at this level -- go back to a checkpoint and restart with lower cfl number");
    }
#endif // PELEC_USE_MOL
}
