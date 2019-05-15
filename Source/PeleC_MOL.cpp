#include <PeleC.H>
#include <PeleC_F.H>

using std::string;
using namespace amrex;

#include <Transport_F.H>
#ifdef PELE_USE_EB
#include <PeleC_init_eb_F.H>
#include <AMReX_MultiCutFab.H>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

// **********************************************************************************************
void
PeleC::getMOLSrcTerm(const amrex::MultiFab& S,
                     amrex::MultiFab&       MOLSrcTerm,
                     amrex::Real            time,
                     amrex::Real            dt,
                     amrex::Real            flux_factor) {
  BL_PROFILE("PeleC::getMOLSrcTerm()");
  BL_PROFILE_VAR_NS("diffusion_stuff", diff);
  if (diffuse_temp == 0
      && diffuse_enth == 0
      && diffuse_spec == 0
      && diffuse_vel  == 0
      && do_hydro == 0)
  {
    MOLSrcTerm.setVal(0,0,NUM_STATE,MOLSrcTerm.nGrow());
    return;
  }
  /**
     Across all conserved state components, compute the method of lines rhs
     = -Div(Flux).  The input state, S, contained the conserved variables, and
     is "fill patched" in the usual AMReX way, where values at Dirichlet boundaries
     actually are assumed to live on the inflow face.

     1. Convert S to Q, primitive variables (since the transport coefficients typically
     depend on mass fractions and temperature). Q then will be face-centered on
     Dirichlet faces.

     2. Evaluate transport coefficients (these also will be face-centered, if Q is).

     3. Evaluate the diffusion operator over all components

     a. Evaluate tangential derivatives for strain terms, over all cells
     b. Replace these with versions that avoid covered cells, if present
     c. Evaluate face-centered diffusion fluxes, and their divergence
     d. Zero fluxes and divergence after the fact if subset of components have diffuse shut off
     (this allows that T be diffused alone, in order to support simple tests)
     e. Compute flux into EB, for temperature (heat flux, when Dirichlet) and momentum
     (when no-slip wall)

     4. (optionally) Evaluate face-centered and EB hyperbolic fluxes, and the divergence
     of the face-centered fluxes.  Increment the diffusion fluxes and divergences with these.

     5. Interpolate combined fluxes to cut face centroids

     6. Replace divergence of face-centered fluxes with hybrid divergence operator, and
     weigthed redistribution.

     Extra notes:

     A. The face-based transport coefficients that are computed with face-based
     Fill-Patched data at Dirichlet boundaries.

     B. Within the routine that computes diffusion fluxes, there is a need for computing
     species enthalpies at cell faces.  In the EGLib model, species enthalpies
     are a function of temperature.  At the moment, in diffterm, we evaluate
     enthalpies at cell centers and then take the face values to be the arithmetic
     average of cell values on either side.  Similarly, mass fractions at cell faces
     are needed to compute the barodiffusion and correction velocity expressions.
     Arithmetic averages are used there as well.  Thus, these face values are thermodynamically
     inconsistent.  Note sure what are the consequences of that.
  */
  int dComp_rhoD = 0;
  int dComp_mu = dComp_rhoD + NumSpec;
  int dComp_xi = dComp_mu + 1;
  int dComp_lambda = dComp_xi + 1;
  int nCompTr = dComp_lambda + 1;
  int do_harmonic = 1;  // TODO: parmparse this
  const Real* dx = geom.CellSize();

  Real dx1 = dx[0];
  for (int d=1; d < BL_SPACEDIM; ++d) {
    dx1 *= dx[d];
  }
  std::array<Real, BL_SPACEDIM> dxD = {D_DECL(dx1, dx1, dx1)};
  const Real *dxDp = &(dxD[0]);

#ifdef PELE_USE_EB
  MultiFab* cost = nullptr;

  if (do_mol_load_balance) cost = &(get_new_data(Work_Estimate_Type));

  EBFluxRegister* fr_as_crse = nullptr;
  if (do_reflux && level < parent->finestLevel()) {
    fr_as_crse = &getFluxReg(level+1);
  }

  EBFluxRegister* fr_as_fine = nullptr;
  if (do_reflux && level > 0) {
    fr_as_fine = &getFluxReg(level);
  }

  int as_crse = (fr_as_crse != nullptr);
  int as_fine = (fr_as_fine != nullptr);
#endif

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    FArrayBox Qfab, Qaux, coeff_cc, Dterm;
    IArrayBox bcMask[BL_SPACEDIM];
    FArrayBox coeff_ec[BL_SPACEDIM], flux_ec[BL_SPACEDIM],
      tander_ec[BL_SPACEDIM], flatn;
    FArrayBox dm_as_fine(Box::TheUnitBox(), NUM_STATE);
    FArrayBox fab_drho_as_crse(Box::TheUnitBox(), NUM_STATE);
    IArrayBox fab_rrflag_as_crse(Box::TheUnitBox());
    

    int flag_nscbc_isAnyPerio = (geom.isAnyPeriodic()) ? 1 : 0; 
    int flag_nscbc_perio[BL_SPACEDIM]; // For 3D, we will know which corners have a periodicity
    for (int d=0; d<BL_SPACEDIM; ++d) {
        flag_nscbc_perio[d] = (Geometry::isPeriodic(d)) ? 1 : 0;
    }
	  const int*  domain_lo = geom.Domain().loVect();
	  const int*  domain_hi = geom.Domain().hiVect();

    for (MFIter mfi(S, MFItInfo().EnableTiling(hydro_tile_size).SetDynamic(true));
         mfi.isValid(); ++mfi) {
#ifdef PELE_USE_EB
      Real wt = ParallelDescriptor::second();
#endif

      const Box  vbox = mfi.tilebox();
      int ng = S.nGrow();
      const Box  gbox = amrex::grow(vbox,ng);
      const Box  cbox = amrex::grow(vbox,ng-1);
      const Box& dbox = geom.Domain();
      
      const int* lo = vbox.loVect();
	  const int* hi = vbox.hiVect();

#ifdef PELE_USE_EB
      const EBFArrayBox& Sfab = static_cast<const EBFArrayBox&>(S[mfi]);

      const auto& flag_fab = Sfab.getEBCellFlagFab();
      FabType typ = flag_fab.getType(cbox);
      if (typ == FabType::covered) {
        MOLSrcTerm[mfi].setVal(0, vbox, 0, NUM_STATE);

        if (do_mol_load_balance) {
          wt = (ParallelDescriptor::second() - wt) / vbox.d_numPts();

          (*cost)[mfi].plus(wt, vbox);
        }
        continue;
      }

      int local_i = mfi.LocalIndex();
      int Ncut = no_eb_in_domain ? 0 : sv_eb_bndry_grad_stencil[local_i].size();
#else
      const FArrayBox& Sfab = S[mfi];
#endif

      BL_PROFILE_VAR_START(diff);
      Qfab.resize(gbox, QVAR);
      int nqaux = NQAUX > 0 ? NQAUX : 1;
      Qaux.resize(gbox, nqaux);
      // Get primitives, Q, including (Y, T, p, rho) from conserved state
      // required for D term
      {
        BL_PROFILE("PeleC::ctoprim call");
        ctoprim(ARLIM_3D(gbox.loVect()), ARLIM_3D(gbox.hiVect()),
                Sfab.dataPtr(), ARLIM_3D(Sfab.loVect()), ARLIM_3D(Sfab.hiVect()),
                Qfab.dataPtr(), ARLIM_3D(Qfab.loVect()), ARLIM_3D(Qfab.hiVect()),
                Qaux.dataPtr(), ARLIM_3D(Qaux.loVect()), ARLIM_3D(Qaux.hiVect()));
      }
      
      
      
      for (int i = 0; i < BL_SPACEDIM ; i++)  {
		    const Box& bxtmp = amrex::surroundingNodes(vbox,i);
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

      if (nscbc_diff == 1)
      {
        impose_NSCBC(lo, hi, domain_lo, domain_hi,
                     BL_TO_FORTRAN(Sfab),
                     BL_TO_FORTRAN(Qfab),
                     BL_TO_FORTRAN(Qaux),
                     D_DECL(BL_TO_FORTRAN(bcMask[0]),
	                          BL_TO_FORTRAN(bcMask[1]),
                            BL_TO_FORTRAN(bcMask[2])),
                     &flag_nscbc_isAnyPerio, flag_nscbc_perio, 
                     &time, dx, &dt);
      }
      
      // Compute transport coefficients, coincident with Q
      {
        BL_PROFILE("PeleC::get_transport_coeffs call");
        coeff_cc.resize(gbox, nCompTr);
        get_transport_coeffs(ARLIM_3D(gbox.loVect()),
                             ARLIM_3D(gbox.hiVect()),
                             BL_TO_FORTRAN_N_3D(Qfab, cQFS),
                             BL_TO_FORTRAN_N_3D(Qfab, cQTEMP),
                             BL_TO_FORTRAN_N_3D(Qfab, cQRHO),
                             BL_TO_FORTRAN_N_3D(coeff_cc, dComp_rhoD),
                             BL_TO_FORTRAN_N_3D(coeff_cc, dComp_mu),
                             BL_TO_FORTRAN_N_3D(coeff_cc, dComp_xi),
                             BL_TO_FORTRAN_N_3D(coeff_cc, dComp_lambda));
      }

      // Container on grown region, for hybrid divergence & redistribution
      Dterm.resize(cbox, NUM_STATE);

      for (int d=0; d<BL_SPACEDIM; ++d) {
        Box ebox = amrex::surroundingNodes(cbox,d);
        coeff_ec[d].resize(ebox,nCompTr);
        flux_ec[d].resize(ebox,NUM_STATE);
        flux_ec[d].setVal(0);
        // Get face-centered transport coefficients
        {
          BL_PROFILE("PeleC::pc_move_transport_coeffs_to_ec call");
          pc_move_transport_coeffs_to_ec(ARLIM_3D(cbox.loVect()),
                                         ARLIM_3D(cbox.hiVect()),
                                         ARLIM_3D(dbox.loVect()),
                                         ARLIM_3D(dbox.hiVect()),
                                         BL_TO_FORTRAN_3D(coeff_cc),
                                         BL_TO_FORTRAN_3D(coeff_ec[d]),
                                         &d, &nCompTr, &do_harmonic);
        }
#if (BL_SPACEDIM > 1)
        int nCompTan = AMREX_D_PICK(1, 2, 6);
        tander_ec[d].resize(ebox, nCompTan); tander_ec[d].setVal(0);
        // Tangential derivatives on faces only needed for velocity diffusion
        if (diffuse_vel == 0) {
          tander_ec[d].setVal(0);
        } else {
          {
            BL_PROFILE("PeleC::pc_compute_tangential_vel_derivs call");
            pc_compute_tangential_vel_derivs(cbox.loVect(),
                                             cbox.hiVect(),
                                             dbox.loVect(),
                                             dbox.hiVect(),
                                             BL_TO_FORTRAN_ANYD(Qfab),
                                             BL_TO_FORTRAN_ANYD(tander_ec[d]),
                                             geom.CellSize(), &d);
          }
#ifdef PELE_USE_EB
          if (typ == FabType::singlevalued) {
            // Reset tangential derivatives to avoid using covered (invalid) data
            if (Ncut > 0) {
              BL_PROFILE("PeleC::pc_compute_tangential_vel_derivs_eb call");
              pc_compute_tangential_vel_derivs_eb(cbox.loVect(),
                                                  cbox.hiVect(),
                                                  dbox.loVect(),
                                                  dbox.hiVect(),
                                                  sv_eb_bndry_geom[mfi.LocalIndex()].data(),
                                                  &Ncut,
                                                  BL_TO_FORTRAN_ANYD(Qfab),
                                                  BL_TO_FORTRAN_ANYD(tander_ec[d]),
                                                  BL_TO_FORTRAN_ANYD(flag_fab),
                                                  geom.CellSize(), &d);
            }
          } else if (typ == FabType::multivalued) {
            amrex::Abort("multi-valued eb tangential derivatives to be implemented");
          }
#endif
        }  // diffuse_vel
#endif
      }  // loop over dimension

      // Compute extensive diffusion fluxes, F.A and (1/Vol).Div(F.A)
      {
        BL_PROFILE("PeleC::pc_diffterm()");
        pc_diffterm(cbox.loVect(),
                    cbox.hiVect(),
                    dbox.loVect(),
                    dbox.hiVect(),
                    BL_TO_FORTRAN_ANYD(Qfab),
                    BL_TO_FORTRAN_N_ANYD(coeff_ec[0], dComp_rhoD),
                    BL_TO_FORTRAN_N_ANYD(coeff_ec[0], dComp_mu),
                    BL_TO_FORTRAN_N_ANYD(coeff_ec[0], dComp_xi),
                    BL_TO_FORTRAN_N_ANYD(coeff_ec[0], dComp_lambda),
#if (BL_SPACEDIM > 1)
                    BL_TO_FORTRAN_ANYD(tander_ec[0]),
#endif
                    BL_TO_FORTRAN_ANYD(area[0][mfi]),
                    BL_TO_FORTRAN_ANYD(flux_ec[0]),
#if (BL_SPACEDIM > 1)
                    BL_TO_FORTRAN_N_ANYD(coeff_ec[1], dComp_rhoD),
                    BL_TO_FORTRAN_N_ANYD(coeff_ec[1], dComp_mu),
                    BL_TO_FORTRAN_N_ANYD(coeff_ec[1], dComp_xi),
                    BL_TO_FORTRAN_N_ANYD(coeff_ec[1], dComp_lambda),
                    BL_TO_FORTRAN_ANYD(tander_ec[1]),
                    BL_TO_FORTRAN_ANYD(area[1][mfi]),
                    BL_TO_FORTRAN_ANYD(flux_ec[1]),
#if (BL_SPACEDIM > 2)
                    BL_TO_FORTRAN_N_ANYD(coeff_ec[2], dComp_rhoD),
                    BL_TO_FORTRAN_N_ANYD(coeff_ec[2], dComp_mu),
                    BL_TO_FORTRAN_N_ANYD(coeff_ec[2], dComp_xi),
                    BL_TO_FORTRAN_N_ANYD(coeff_ec[2], dComp_lambda),
                    BL_TO_FORTRAN_ANYD(tander_ec[2]),
                    BL_TO_FORTRAN_ANYD(area[2][mfi]),
                    BL_TO_FORTRAN_ANYD(flux_ec[2]),
#endif
#endif
                    BL_TO_FORTRAN_ANYD(volume[mfi]),
                    BL_TO_FORTRAN_ANYD(Dterm),
                    geom.CellSize());
      }

      // Shut off unwanted diffusion after the fact
      //    ick! Under normal conditions, you either have diffusion on all or
      //      none, so this shouldn't be done this way.  However, the regression
      //      test for diffusion works by diffusing only temperature through
      //      this process.  Ideally, we'd redo that test to diffuse a passive
      //      scalar instead....
          
      if (diffuse_temp == 0 && diffuse_enth == 0) {
        Dterm.setVal(0, Eden);
        Dterm.setVal(0, Eint);
        for (int d = 0; d < BL_SPACEDIM; d++) {
          flux_ec[d].setVal(0, Eden);
          flux_ec[d].setVal(0, Eint);
        }
      }
      if (diffuse_spec == 0) {
        Dterm.setVal(0, Dterm.box(), FirstSpec, NumSpec);
        for (int d = 0; d < BL_SPACEDIM ; d++) {
          flux_ec[d].setVal(0, flux_ec[d].box(), FirstSpec, NumSpec);
        }
      }

      if (diffuse_vel  == 0) {
        Dterm.setVal(0, Dterm.box(), Xmom, 3);
        for (int d = 0; d < BL_SPACEDIM; d++) {
          flux_ec[d].setVal(0, flux_ec[d].box(), Xmom, 3);
        }
      }

#ifdef PELE_USE_EB
      //  Set extensive flux at embedded boundary, potentially
      //  non-zero only for heat flux on isothermal boundaries,
      //  and momentum fluxes at no-slip walls
      if (typ == FabType::singlevalued && Ncut > 0) {
        sv_eb_flux[local_i].setVal(0);  // Default to Neumann for all fields

        int Nvals = sv_eb_bcval[local_i].numPts();
        int Nflux = sv_eb_flux[local_i].numPts();
        BL_ASSERT(Nvals == Ncut);
        BL_ASSERT(Nflux == Ncut);

        if (eb_isothermal && (diffuse_temp != 0 || diffuse_enth != 0)) {
          // Compute heat flux at EB wall
          int nComp = 1;
          sv_eb_bcval[local_i].setVal(eb_boundary_T, cQTEMP);
              
          Box box_to_apply = mfi.growntilebox(2);
          {
            BL_PROFILE("PeleC::pc_apply_eb_boundry_flux_stencil call");
            pc_apply_eb_boundry_flux_stencil(BL_TO_FORTRAN_BOX(box_to_apply),
                                             sv_eb_bndry_grad_stencil[local_i].data(),
                                             &Ncut,
                                             BL_TO_FORTRAN_N_ANYD(Qfab, cQTEMP),
                                             BL_TO_FORTRAN_N_ANYD(coeff_cc, dComp_lambda),
                                             sv_eb_bcval[local_i].dataPtr(cQTEMP),
                                             &Nvals,
                                             sv_eb_flux[local_i].dataPtr(Eden),
                                             &Nflux, &nComp);
          }
        }
        // Compute momentum transfer at no-slip EB wall
        if (eb_noslip && diffuse_vel == 1) {
          int nComp = BL_SPACEDIM;
          sv_eb_bcval[local_i].setVal(0, cQU, BL_SPACEDIM);

          Box box_to_apply = mfi.growntilebox(2);
          {
            BL_PROFILE("PeleC::pc_apply_eb_boundry_visc_flux_stencil call");
            pc_apply_eb_boundry_visc_flux_stencil(BL_TO_FORTRAN_BOX(box_to_apply),
                                                  sv_eb_bndry_grad_stencil[local_i].data(),
                                                  &Ncut,
                                                  sv_eb_bndry_geom[local_i].data(), &Ncut,
                                                  BL_TO_FORTRAN_N_ANYD(Qfab, cQU),
                                                  BL_TO_FORTRAN_N_ANYD(coeff_cc, dComp_mu),
                                                  BL_TO_FORTRAN_N_ANYD(coeff_cc, dComp_xi),
                                                  sv_eb_bcval[local_i].dataPtr(cQU), &Nvals,
                                                  sv_eb_flux[local_i].dataPtr(Xmom), &Nflux,
                                                  &nComp);
          }
        }
      }
#endif

      BL_PROFILE_VAR_STOP(diff);

#ifdef PELEC_USE_MOL
      /* At this point flux_ec contains the diffusive fluxes in each direction
         at face centers for the (potentially partially covered) grid-aligned 
         faces and sv_eb_flux contains the flux for the cut faces. Before 
         computing hybrid divergence, comptue and add in the hydro fluxes. 
         Also, Dterm currently contains the divergence of the face-centered
         diffusion fluxes.  Increment this with the divergence of the 
         face-centered hyperbloic fluxes.
      */
      if (do_hydro && do_mol_AD) 
      {
        flatn.resize(cbox,1);
        flatn.setVal(1.0);  // Set flattening to 1.0
#ifdef PELEC_USE_EB
        int nFlux = sv_eb_flux.size()==0 ? 0 : sv_eb_flux[local_i].numPts();
        const EBBndryGeom* sv_ebbg_ptr = (Ncut>0 ? sv_eb_bndry_geom[local_i].data() : 0);
        Real* sv_eb_flux_ptr = (nFlux>0 ? sv_eb_flux[local_i].dataPtr() : 0);
#endif

        { // Get face-centered hyperbolic fluxes and their divergences.
          // Get hyp flux at EB wall
          BL_PROFILE("PeleC::pc_hyp_mol_flux call");
          pc_hyp_mol_flux(vbox.loVect(), vbox.hiVect(),
                          geom.Domain().loVect(), geom.Domain().hiVect(),
                          BL_TO_FORTRAN_3D(Qfab),
                          BL_TO_FORTRAN_3D(Qaux),
                          BL_TO_FORTRAN_ANYD(area[0][mfi]),
                          BL_TO_FORTRAN_3D(flux_ec[0]),
#if (BL_SPACEDIM > 1)
                          BL_TO_FORTRAN_ANYD(area[1][mfi]),
                          BL_TO_FORTRAN_3D(flux_ec[1]),
#if (BL_SPACEDIM > 2)
                          BL_TO_FORTRAN_ANYD(area[2][mfi]),
                          BL_TO_FORTRAN_3D(flux_ec[2]),
#endif
#endif
                          BL_TO_FORTRAN_3D(flatn),
                          BL_TO_FORTRAN_ANYD(volume[mfi]),
                          BL_TO_FORTRAN_3D(Dterm),
#ifdef PELEC_USE_EB
                          BL_TO_FORTRAN_ANYD(vfrac[mfi]),
                          BL_TO_FORTRAN_ANYD(flag_fab),
                          sv_ebbg_ptr, &Ncut,
                          sv_eb_flux_ptr, &nFlux,
#endif
                          geom.CellSize());
        }
      }
#endif

#ifdef PELEC_USE_EB
      if (typ == FabType::singlevalued) {
        /* Interpolate fluxes from face centers to face centroids
         * Note that hybrid divergence and redistribution algorithms require that we
         *   be able to compute the conservative divergence on 2 grow cells, so we
         *   need interpolated fluxes on 2 grow cells, and therefore we need face
         *   centered fluxes on 3.
         */

        for (int idir=0; idir < BL_SPACEDIM; ++idir) {
          int Nsten = flux_interp_stencil[idir][local_i].size();
          int in_place = 1;
          const Box valid_interped_flux_box =
            Box(amrex::grow(vbox, 2)).surroundingNodes(idir);
          {
            BL_PROFILE("PeleC::pc_apply_face_stencil call");
            pc_apply_face_stencil(BL_TO_FORTRAN_BOX(valid_interped_flux_box),
                                  BL_TO_FORTRAN_BOX(stencil_volume_box),
                                  flux_interp_stencil[idir][local_i].data(),
                                  &Nsten, &idir,
                                  BL_TO_FORTRAN_ANYD(flux_ec[idir]),
                                  BL_TO_FORTRAN_ANYD(flux_ec[idir]),
                                  &NUM_STATE, &in_place);
          }
        }

        // Get "hybrid flux divergence" and redistribute
        //
        // This operation takes as input centroid-centered fluxes and a corresponding
        //  divergence on three grid cells.  Actually, we assume that div=(1/VOL)Div(flux)
        //  (VOL = volume of the full cells), and that flux is EXTENSIVE, weighted with the
        //  full face areas.
        //
        // Upon return:
        // div = kappa.(1/Vol) Div(FluxC.Area)  Vol = kappa.VOL, Area=aperture.AREA,
        //                  defined over the valid box

        // TODO: Rework this for r-z, if applicable
        Real VOL = 1;
        for (int idir = 0; idir < BL_SPACEDIM; ++idir)
          VOL *= geom.CellSize()[idir];

        // Set weighting for redistribution
        const FArrayBox& W = vfrac[mfi];
        int wComp = 0;

        int Nflux = sv_eb_flux[local_i].numPts();
        {
          FArrayBox* p_drho_as_crse = (fr_as_crse) ?
            fr_as_crse->getCrseData(mfi) : &fab_drho_as_crse;
          const IArrayBox* p_rrflag_as_crse = (fr_as_crse) ?
            fr_as_crse->getCrseFlag(mfi) : &fab_rrflag_as_crse;

          if (fr_as_fine) {
            dm_as_fine.resize(amrex::grow(vbox, 1), NUM_STATE);
          }
          BL_PROFILE("PeleC::pc_fix_div_and_redistribute call");
          pc_fix_div_and_redistribute(BL_TO_FORTRAN_BOX(vbox),
                                      sv_eb_bndry_geom[local_i].data(), &Ncut,
                                      BL_TO_FORTRAN_ANYD(flag_fab),
                                      D_DECL(BL_TO_FORTRAN_ANYD(flux_ec[0]),
                                             BL_TO_FORTRAN_ANYD(flux_ec[1]),
                                             BL_TO_FORTRAN_ANYD(flux_ec[2])),
                                      sv_eb_flux[local_i].dataPtr(), &Nflux,
                                      BL_TO_FORTRAN_ANYD(Dterm),
                                      BL_TO_FORTRAN_N_ANYD(W, wComp),
                                      BL_TO_FORTRAN_ANYD(vfrac[mfi]),
                                      &VOL, &NUM_STATE,
                                      &as_crse,
                                      BL_TO_FORTRAN_ANYD(*p_drho_as_crse),
                                      BL_TO_FORTRAN_ANYD(*p_rrflag_as_crse),
                                      &as_fine,
                                      BL_TO_FORTRAN_ANYD(dm_as_fine),
                                      BL_TO_FORTRAN_ANYD(level_mask[mfi]), &dt);
        }

        if (do_reflux && flux_factor != 0) {
          for (int d = 0; d < BL_SPACEDIM; d++) {
            flux_ec[d].mult(flux_factor);
          }

          if (fr_as_crse) {
#if (BL_SPACEDIM > 2)
            fr_as_crse->CrseAdd(mfi,
                                {D_DECL(&flux_ec[0], &flux_ec[1], &flux_ec[2])},
                                dxDp, dt, vfrac[mfi],
                                {&((*areafrac[0])[mfi]),
                                    &((*areafrac[1])[mfi]),
                                    &((*areafrac[2])[mfi])});
#else
            // TODO: EBfluxregisters are designed only for 3D, need for 2D
            Print() << "WARNING:Re redistribution crseadd for EB not implemented\n";
#endif
          }

          if (fr_as_fine) {
#if (BL_SPACEDIM > 2)
            fr_as_fine->FineAdd(mfi,
                                {D_DECL(&flux_ec[0], &flux_ec[1], &flux_ec[2])},
                                dxDp, dt,
                                vfrac[mfi],
                                {&((*areafrac[0])[mfi]),
                                    &((*areafrac[1])[mfi]),
                                    &((*areafrac[2])[mfi])},
                                dm_as_fine);
#else
            // TODO: EBfluxregisters are designed only for 3D, need for 2D
            Print() << "WARNING:Re redistribution fineadd for EB not implemented in 2D\n";
#endif
          }
        }
      } else if (typ != FabType::regular) {  // Single valued if loop
        amrex::Abort("multi-valued eb boundary fluxes to be implemented");
      }
#endif  //  PELEC_USE_EB ifdef

      MOLSrcTerm[mfi].setVal(0, vbox, 0, NUM_STATE);
      MOLSrcTerm[mfi].copy(Dterm, vbox, 0, vbox, 0, NUM_STATE);

#ifdef PELEC_USE_EB
      // do regular flux reg ops
      if (do_reflux && flux_factor != 0 && typ == FabType::regular) 
#else
        if (do_reflux && flux_factor != 0)  // no eb in problem
#endif
        {
          for (int d = 0; d < BL_SPACEDIM ; d++) {
            flux_ec[d].mult(flux_factor);
          }

          if (level < parent->finestLevel()) {
            getFluxReg(level+1).CrseAdd(mfi,
                                        {D_DECL(&flux_ec[0], &flux_ec[1], &flux_ec[2])},
                                        dxDp, dt);
          }

          if (level > 0) {
            getFluxReg(level).FineAdd(mfi,
                                      {D_DECL(&flux_ec[0], &flux_ec[1], &flux_ec[2])},
                                      dxDp, dt);
          }
        }

#ifdef PELEC_USE_EB
      if (do_mol_load_balance) {
        wt = (ParallelDescriptor::second() - wt) / vbox.d_numPts();
        (*cost)[mfi].plus(wt, vbox);
      }
#endif
    }  // End of MFIter scope
  }  // End of OMP scope

  // Extrapolate to ghost cells
  if (MOLSrcTerm.nGrow() > 0) {
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(MOLSrcTerm, hydro_tile_size); mfi.isValid(); ++mfi) {
      BL_PROFILE("PeleC::diffextrap calls");

      const Box& bx = mfi.validbox();
      pc_diffextrap(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
                    BL_TO_FORTRAN_N_3D(MOLSrcTerm[mfi], Xmom), &amrex::SpaceDim);

      pc_diffextrap(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
                    BL_TO_FORTRAN_N_3D(MOLSrcTerm[mfi], FirstSpec), &NumSpec);

      const int one = 1;
      pc_diffextrap(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
                    BL_TO_FORTRAN_N_3D(MOLSrcTerm[mfi], Eden), &one);
    }
  }
}
