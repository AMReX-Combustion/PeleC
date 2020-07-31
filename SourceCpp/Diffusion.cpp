#include "Diffusion.H"

void
PeleC::getMOLSrcTerm(
  const amrex::MultiFab& S,
  amrex::MultiFab& MOLSrcTerm,
  amrex::Real time,
  amrex::Real dt,
  amrex::Real flux_factor)
{
  BL_PROFILE("PeleC::getMOLSrcTerm()");
  BL_PROFILE_VAR_NS("diffusion_stuff", diff);
  if (
    diffuse_temp == 0 && diffuse_enth == 0 && diffuse_spec == 0 &&
    diffuse_vel == 0 && do_hydro == 0) {
    MOLSrcTerm.setVal(0, 0, NVAR, MOLSrcTerm.nGrow());
    return;
  }
  /**
     Across all conserved state components, compute the method of lines rhs
     = -Div(Flux).  The input state, S, contained the conserved variables, and
     is "fill patched" in the usual AMReX way, where values at Dirichlet
     boundaries actually are assumed to live on the inflow face.

     1. Convert S to Q, primitive variables (since the transport coefficients
     typically depend on mass fractions and temperature). Q then will be
     face-centered on Dirichlet faces.

     2. Evaluate transport coefficients (these also will be face-centered, if Q
     is).

     3. Evaluate the diffusion operator over all components

     a. Evaluate tangential derivatives for strain terms, over all cells
     b. Replace these with versions that avoid covered cells, if present
     c. Evaluate face-centered diffusion fluxes, and their divergence
     d. Zero fluxes and divergence after the fact if subset of components have
     diffuse shut off (this allows that T be diffused alone, in order to support
     simple tests)

     4. Replace divergence of face-centered fluxes with hybrid divergence
     operator, and weigthed redistribution.

     Extra notes:

     A. The face-based transport coefficients that are computed with face-based
     Fill-Patched data at Dirichlet boundaries.

     B. Within the routine that computes diffusion fluxes, there is a need for
     computing species enthalpies at cell faces.  In the EGLib model, species
     enthalpies are a function of temperature.  At the moment, in diffterm, we
     evaluate enthalpies at cell centers and then take the face values to be the
     arithmetic average of cell values on either side.  Similarly, mass
     fractions at cell faces are needed to compute the barodiffusion and
     correction velocity expressions. Arithmetic averages are used there as
     well.  Thus, these face values are thermodynamically inconsistent.  Note
     sure what are the consequences of that.

     The non-GPU version has EB, this does not at the moment. The Techniques for
     EB will need to be carefully considered before implementation for
     accelerators.
  */
  const int nCompTr = dComp_lambda + 1;
  const int do_harmonic = 1; // TODO: parmparse this
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

  amrex::Real dx1 = dx[0];
  for (int dir = 1; dir < AMREX_SPACEDIM; ++dir) {
    dx1 *= dx[dir];
  }
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dxD = {
    AMREX_D_DECL(dx1, dx1, dx1)};

  // Fetch some gpu arrays
  prefetchToDevice(S);
  prefetchToDevice(MOLSrcTerm);

#ifdef PELEC_USE_EB
  auto const& fact =
    dynamic_cast<amrex::EBFArrayBoxFactory const&>(S.Factory());
  auto const& flags = fact.getMultiEBCellFlagFab();
  // amrex::Elixir flags_eli = flags.elixir();
  amrex::MultiFab* cost = nullptr;

  if (do_mol_load_balance)
    cost = &(get_new_data(Work_Estimate_Type));

  amrex::EBFluxRegister* fr_as_crse = nullptr;
  if (do_reflux && level < parent->finestLevel()) {
    fr_as_crse = &getFluxReg(level + 1);
  }

  amrex::EBFluxRegister* fr_as_fine = nullptr;
  if (do_reflux && level > 0) {
    fr_as_fine = &getFluxReg(level);
  }

  const bool as_crse = (fr_as_crse != nullptr);
  const bool as_fine = (fr_as_fine != nullptr);
#endif

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  {
    // amrex::IArrayBox bcMask[AMREX_SPACEDIM];

    int flag_nscbc_isAnyPerio = (geom.isAnyPeriodic()) ? 1 : 0;
    int flag_nscbc_perio[AMREX_SPACEDIM]; // For 3D, we will know which corners
                                          // have a periodicity
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      flag_nscbc_perio[dir] =
        (amrex::DefaultGeometry().isPeriodic(dir)) ? 1 : 0;
    }
    const int* domain_lo = geom.Domain().loVect();
    const int* domain_hi = geom.Domain().hiVect();

    for (amrex::MFIter mfi(MOLSrcTerm, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
      const amrex::Box vbox = mfi.tilebox();
      int ng = S.nGrow();
      const amrex::Box gbox = amrex::grow(vbox, ng);
      const amrex::Box cbox = amrex::grow(vbox, ng - 1);
      auto const& MOLSrc = MOLSrcTerm.array(mfi);

#ifdef PELEC_USE_EB
      amrex::Real wt = amrex::ParallelDescriptor::second();
      const auto& flag_fab = flags[mfi];
      // amrex::Elixir flag_fab_eli = flag_fab.elixir();
      amrex::FabType typ = flag_fab.getType(vbox);
      if (typ == amrex::FabType::covered) {
        setV(vbox, NVAR, MOLSrc, 0);
        if (do_mol_load_balance) {
          wt = (amrex::ParallelDescriptor::second() - wt) / vbox.d_numPts();
          (*cost)[mfi].plus<amrex::RunOn::Device>(wt, vbox);
        }
        continue;
      }

      // TODO: Add check that this is nextra-1
      //       (better: fix bounds on ebflux computation in hyperbolic routine
      //                to be a constant, and make sure this matches it)
      const amrex::Box ebfluxbox = amrex::grow(vbox, 2);

      int local_i = mfi.LocalIndex();
      int Ncut = (!eb_in_domain) ? 0 : sv_eb_bndry_grad_stencil[local_i].size();
      SparseData<amrex::Real, EBBndrySten> eb_flux_thdlocal;
      eb_flux_thdlocal.define(sv_eb_bndry_grad_stencil[local_i], NVAR);
      auto* d_sv_eb_bndry_geom =
        (Ncut > 0 ? sv_eb_bndry_geom[local_i].data() : 0);
#endif

      const int* lo = vbox.loVect();
      const int* hi = vbox.hiVect();

      BL_PROFILE_VAR_START(diff);
      int nqaux = NQAUX > 0 ? NQAUX : 1;
      amrex::FArrayBox q(gbox, QVAR), qaux(gbox, nqaux),
        coeff_cc(gbox, nCompTr);
      amrex::Elixir qeli = q.elixir();
      amrex::Elixir qauxeli = qaux.elixir();
      amrex::Elixir coefeli = coeff_cc.elixir();
      auto const& s = S.array(mfi);
      auto const& qar = q.array();
      auto const& qauxar = qaux.array();

      // Get primitives, Q, including (Y, T, p, rho) from conserved state
      // required for D term
      {
        BL_PROFILE("PeleC::ctoprim()");
        amrex::ParallelFor(
          gbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            pc_ctoprim(i, j, k, s, qar, qauxar);
          });
      }
// TODO deal with NSCBC
#if 0       
      for (int dir = 0; dir < AMREX_SPACEDIM ; dir++)  {
        const amrex::Box& bxtmp = amrex::surroundingNodes(vbox,dir);
        amrex::Box TestBox(bxtmp);
        for(int d=0; d<AMREX_SPACEDIM; ++d) {
          if (dir!=d) TestBox.grow(d,1);
        }
        
        bcMask[dir].resize(TestBox,1);
        bcMask[dir].setVal(0);
	    }
      
      // Becase bcMask is read in the Riemann solver in any case,
      // here we put physbc values in the appropriate faces for the non-nscbc case
      set_bc_mask(lo, hi, domain_lo, domain_hi,
                  AMREX_D_DECL(AMREX_TO_FORTRAN(bcMask[0]),
	                       AMREX_TO_FORTRAN(bcMask[1]),
                         AMREX_TO_FORTRAN(bcMask[2])));

      if (nscbc_diff == 1)
      {
        impose_NSCBC(lo, hi, domain_lo, domain_hi,
                     AMREX_TO_FORTRAN(Sfab),
                     AMREX_TO_FORTRAN(q.fab()),
                     AMREX_TO_FORTRAN(qaux.fab()),
                     AMREX_D_DECL(AMREX_TO_FORTRAN(bcMask[0]),
                            AMREX_TO_FORTRAN(bcMask[1]),
                            AMREX_TO_FORTRAN(bcMask[2])),
                     &flag_nscbc_isAnyPerio, flag_nscbc_perio, 
                     &time, dx, &dt);
      }
#endif

      // Compute transport coefficients, coincident with Q
      auto const& coe_cc = coeff_cc.array();
      {
        auto const& qar_yin = q.array(QFS);
        auto const& qar_Tin = q.array(QTEMP);
        auto const& qar_rhoin = q.array(QRHO);
        auto const& coe_rhoD = coeff_cc.array(dComp_rhoD);
        auto const& coe_mu = coeff_cc.array(dComp_mu);
        auto const& coe_xi = coeff_cc.array(dComp_xi);
        auto const& coe_lambda = coeff_cc.array(dComp_lambda);
        BL_PROFILE("PeleC::get_transport_coeffs()");
        // Get Transport coefs on GPU.
        amrex::launch(gbox, [=] AMREX_GPU_DEVICE(amrex::Box const& tbx) {
          get_transport_coeffs(
            tbx, qar_yin, qar_Tin, qar_rhoin, coe_rhoD, coe_mu, coe_xi,
            coe_lambda);
        });
      }

      amrex::FArrayBox flux_ec[AMREX_SPACEDIM];
      amrex::Elixir flux_eli[AMREX_SPACEDIM];
      const amrex::Box eboxes[AMREX_SPACEDIM] = {AMREX_D_DECL(
        amrex::surroundingNodes(cbox, 0), amrex::surroundingNodes(cbox, 1),
        amrex::surroundingNodes(cbox, 2))};
      amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM> flx;
      const amrex::GpuArray<
        const amrex::Array4<const amrex::Real>, AMREX_SPACEDIM>
        a{AMREX_D_DECL(
          area[0].array(mfi), area[1].array(mfi), area[2].array(mfi))};
      for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
        flux_ec[dir].resize(eboxes[dir], NVAR);
        flux_eli[dir] = flux_ec[dir].elixir();
        flx[dir] = flux_ec[dir].array();
        setV(eboxes[dir], NVAR, flx[dir], 0);
      }

      amrex::FArrayBox Dfab(cbox, NVAR);
      amrex::Elixir Dfab_eli = Dfab.elixir();
      auto const& Dterm = Dfab.array();
      setV(cbox, NVAR, Dterm, 0.0);

      pc_compute_diffusion_flux(
        cbox, qar, coe_cc, flx, a, dx, do_harmonic
#ifdef PELEC_USE_EB
        ,
        typ, Ncut, d_sv_eb_bndry_geom, flags.array(mfi)
#endif
      );

      // Compute flux divergence (1/Vol).Div(F.A)
      {
        BL_PROFILE("PeleC::pc_flux_div()");
        auto const& vol = volume.array(mfi);
        amrex::ParallelFor(
          cbox, NVAR,
          [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
            pc_flux_div(
              i, j, k, n, AMREX_D_DECL(flx[0], flx[1], flx[2]), vol, Dterm);
          });
      }

      // Shut off unwanted diffusion after the fact
      //    ick! Under normal conditions, you either have diffusion on all or
      //      none, so this shouldn't be done this way.  However, the regression
      //      test for diffusion works by diffusing only temperature through
      //      this process.  Ideally, we'd redo that test to diffuse a passive
      //      scalar instead....

      if (diffuse_temp == 0 && diffuse_enth == 0) {
        setC(cbox, Eden, Eint, Dterm, 0.0);
        for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
          setC(eboxes[dir], Eden, Eint, flx[dir], 0.0);
        }
      }
      if (diffuse_spec == 0) {
        setC(cbox, FirstSpec, FirstSpec + NUM_SPECIES, Dterm, 0.0);
        for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
          setC(eboxes[dir], FirstSpec, FirstSpec + NUM_SPECIES, flx[dir], 0.0);
        }
      }

      if (diffuse_vel == 0) {
        setC(cbox, Xmom, Xmom + 3, Dterm, 0.0);
        for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
          setC(eboxes[dir], Xmom, Xmom + 3, flx[dir], 0.0);
        }
      }

#ifdef PELEC_USE_EB
      //  Set extensive flux at embedded boundary, potentially
      //  non-zero only for heat flux on isothermal boundaries,
      //  and momentum fluxes at no-slip walls
      const int nFlux =
        sv_eb_flux.size() == 0 ? 0 : sv_eb_flux[local_i].numPts();
      if (typ == amrex::FabType::singlevalued && Ncut > 0) {
        eb_flux_thdlocal.setVal(0); // Default to Neumann for all fields

        int Nvals = sv_eb_bcval[local_i].numPts();

        AMREX_ASSERT(Nvals == Ncut);
        AMREX_ASSERT(nFlux == Ncut);

        if (eb_isothermal && (diffuse_temp != 0 || diffuse_enth != 0)) {
          {
            BL_PROFILE("PeleC::pc_apply_eb_boundry_flux_stencil()");
            pc_apply_eb_boundry_flux_stencil(
              ebfluxbox, sv_eb_bndry_grad_stencil[local_i].data(), Ncut, qar,
              QTEMP, coe_cc, dComp_lambda, sv_eb_bcval[local_i].dataPtr(QTEMP),
              Nvals, eb_flux_thdlocal.dataPtr(Eden), nFlux, 1);
          }
        }
        // Compute momentum transfer at no-slip EB wall
        if (eb_noslip && diffuse_vel == 1) {
          {
            BL_PROFILE("PeleC::pc_apply_eb_boundry_visc_flux_stencil()");
            pc_apply_eb_boundry_visc_flux_stencil(
              ebfluxbox, sv_eb_bndry_grad_stencil[local_i].data(), Ncut,
              d_sv_eb_bndry_geom, Ncut, qar, coe_cc,
              sv_eb_bcval[local_i].dataPtr(QU), Nvals,
              eb_flux_thdlocal.dataPtr(Xmom), nFlux);
          }
        }
      }
#endif

      BL_PROFILE_VAR_STOP(diff);

      /* At this point flux_ec contains the diffusive fluxes in each direction
         at face centers for the (potentially partially covered) grid-aligned
         faces and eb_flux_thdlocal contains the flux for the cut faces. Before
         computing hybrid divergence, comptue and add in the hydro fluxes.
         Also, Dterm currently contains the divergence of the face-centered
         diffusion fluxes.  Increment this with the divergence of the
         face-centered hyperbloic fluxes.
      */
      if (do_hydro && do_mol) {
        // amrex::FArrayBox flatn(cbox, 1);
        // amrex::Elixir flatn_eli;
        // flatn_eli = flatn.elixir();
        // flatn.setVal(1.0); // Set flattening to 1.0

        // save off the diffusion source term and fluxes (don't want to filter
        // these)
        amrex::FArrayBox diffusion_flux[AMREX_SPACEDIM];
        amrex::Elixir diffusion_flux_eli[AMREX_SPACEDIM];
        amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM>
          diffusion_flux_arr;
        if (use_explicit_filter) {
          for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
            diffusion_flux[dir].resize(flux_ec[dir].box(), NVAR);
            diffusion_flux_eli[dir] = diffusion_flux[dir].elixir();
            diffusion_flux_arr[dir] = diffusion_flux[dir].array();
            copy_array4(
              flux_ec[dir].box(), flux_ec[dir].nComp(), flx[dir],
              diffusion_flux_arr[dir]);
          }
        }

        { // Get face-centered hyperbolic fluxes and their divergences.
          // Get hyp flux at EB wall
          BL_PROFILE("PeleC::pc_hyp_mol_flux()");
#ifdef PELEC_USE_EB
          amrex::Real* d_eb_flux_thdlocal =
            (nFlux > 0 ? eb_flux_thdlocal.dataPtr() : 0);
#endif
          auto const& vol = volume.array(mfi);
          pc_compute_hyp_mol_flux(
            cbox, qar, qauxar, flx, a, dx, plm_iorder
#ifdef PELEC_USE_EB
            ,
            eb_small_vfrac, vfrac.array(mfi), flags.array(mfi),
            d_sv_eb_bndry_geom, Ncut, d_eb_flux_thdlocal, nFlux
#endif
          );
        }

        // Filter hydro source term and fluxes here
        if (use_explicit_filter) {
          // Get the hydro term
          amrex::FArrayBox hydro_flux[AMREX_SPACEDIM];
          amrex::Elixir hydro_flux_eli[AMREX_SPACEDIM];
          amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM>
            hydro_flux_arr;
          for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
            hydro_flux[dir].resize(flux_ec[dir].box(), NVAR);
            hydro_flux_eli[dir] = hydro_flux[dir].elixir();
            hydro_flux_arr[dir] = hydro_flux[dir].array();
            lincomb_array4(
              flux_ec[dir].box(), Density, NVAR, flx[dir],
              diffusion_flux_arr[dir], 1.0, -1.0, hydro_flux_arr[dir]);
          }

          // Filter
          const amrex::Box fbox = amrex::grow(cbox, -nGrowF);
          for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
            const amrex::Box& bxtmp = amrex::surroundingNodes(fbox, dir);
            amrex::FArrayBox filtered_hydro_flux;
            filtered_hydro_flux.resize(bxtmp, NVAR);
            amrex::Elixir filtered_hydro_flux_eli =
              filtered_hydro_flux.elixir();
            les_filter.apply_filter(
              bxtmp, hydro_flux[dir], filtered_hydro_flux, Density, NVAR);

            setV(bxtmp, hydro_flux[dir].nComp(), hydro_flux_arr[dir], 0.0);
            copy_array4(
              bxtmp, hydro_flux[dir].nComp(), filtered_hydro_flux.array(),
              hydro_flux_arr[dir]);
          }

          // Combine with diffusion
          for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
            lincomb_array4(
              diffusion_flux[dir].box(), Density, NVAR, diffusion_flux_arr[dir],
              hydro_flux_arr[dir], 1.0, 1.0, flx[dir]);
          }
        }

        // Compute flux divergence (1/Vol).Div(F.A)
        {
          BL_PROFILE("PeleC::pc_flux_div()");
          auto const& vol = volume.array(mfi);
          amrex::ParallelFor(
            cbox, NVAR,
            [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
              pc_flux_div(
                i, j, k, n, AMREX_D_DECL(flx[0], flx[1], flx[2]), vol, Dterm);
            });
        }
      }

#ifdef AMREX_USE_GPU
      auto device = amrex::RunOn::Gpu;
#else
      auto device = amrex::RunOn::Cpu;
#endif

#ifdef PELEC_USE_EB
      amrex::Gpu::DeviceVector<int> v_eb_tile_mask(Ncut, 0);
      int* eb_tile_mask = v_eb_tile_mask.dataPtr();
      amrex::ParallelFor(Ncut, [=] AMREX_GPU_DEVICE(int icut) {
        if (ebfluxbox.contains(d_sv_eb_bndry_geom[icut].iv)) {
          eb_tile_mask[icut] = 1;
        }
      });
      if (typ == amrex::FabType::singlevalued) {
        sv_eb_flux[local_i].merge(eb_flux_thdlocal, 0, NVAR, v_eb_tile_mask);
      }

      amrex::FArrayBox dm_as_fine, fab_drho_as_crse;
      amrex::IArrayBox fab_rrflag_as_crse;
      amrex::Elixir dm_as_fine_eli, fab_drho_as_crse_eli,
        fab_rrflag_as_crse_eli;
      if (typ == amrex::FabType::singlevalued) {
        /* Interpolate fluxes from face centers to face centroids
         * Note that hybrid divergence and redistribution algorithms require
         * that we be able to compute the conservative divergence on 2 grow
         * cells, so we need interpolated fluxes on 2 grow cells, and therefore
         * we need face centered fluxes on 3.
         */

        {
          BL_PROFILE("PeleC::pc_apply_face_stencil()");
          for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            int Nsten = flux_interp_stencil[dir][local_i].size();
            int in_place = 1;
            const amrex::Box valid_interped_flux_box =
              amrex::Box(ebfluxbox).surroundingNodes(dir);
            pc_apply_face_stencil(
              valid_interped_flux_box, stencil_volume_box,
              flux_interp_stencil[dir][local_i].data(), Nsten, dir, NVAR,
              flx[dir]);
          }
        }

        // Get "hybrid flux divergence" and redistribute
        //
        // This operation takes as input centroid-centered fluxes and a
        // corresponding
        //  divergence on three grid cells.  Actually, we assume that
        //  div=(1/VOL)Div(flux) (VOL = volume of the full cells), and that flux
        //  is EXTENSIVE, weighted with the full face areas.
        //
        // Upon return:
        // div = kappa.(1/Vol) Div(FluxC.Area)  Vol = kappa.VOL,
        // Area=aperture.AREA,
        //                  defined over the valid box

        // TODO: Rework this for r-z, if applicable
        amrex::Real vol = 1;
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
          vol *= geom.CellSize()[dir];

        // Set weighting for redistribution
        auto const& W = vfrac.array(mfi);
        int wComp = 0;

        dm_as_fine.resize(amrex::Box::TheUnitBox(), NVAR);
        dm_as_fine_eli = dm_as_fine.elixir();
        fab_drho_as_crse.resize(amrex::Box::TheUnitBox(), NVAR);
        fab_drho_as_crse_eli = fab_drho_as_crse.elixir();
        fab_rrflag_as_crse.resize(amrex::Box::TheUnitBox());
        fab_rrflag_as_crse_eli = fab_rrflag_as_crse.elixir();
        {
          amrex::FArrayBox* p_drho_as_crse =
            (fr_as_crse) ? fr_as_crse->getCrseData(mfi) : &fab_drho_as_crse;
          const amrex::IArrayBox* p_rrflag_as_crse =
            (fr_as_crse) ? fr_as_crse->getCrseFlag(mfi) : &fab_rrflag_as_crse;

          if (fr_as_fine) {
            dm_as_fine.resize(amrex::grow(vbox, 1), NVAR);
            dm_as_fine_eli = dm_as_fine.elixir();
            dm_as_fine.setVal<amrex::RunOn::Device>(0.0);
          }
          BL_PROFILE("PeleC::pc_fix_div_and_redistribute()");
          pc_fix_div_and_redistribute(
            vbox, vol, dt, NVAR, eb_small_vfrac, levmsk_notcovered,
            d_sv_eb_bndry_geom, Ncut, flags.array(mfi),
            AMREX_D_DECL(flx[0], flx[1], flx[2]), sv_eb_flux[local_i].dataPtr(),
            nFlux, vfrac.array(mfi), W, as_crse, as_fine, level_mask.array(mfi),
            (*p_rrflag_as_crse).array(), Dterm, (*p_drho_as_crse).array(),
            dm_as_fine.array());
        }

        if (do_reflux && flux_factor != 0) {
          for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
            flux_ec[dir].mult<amrex::RunOn::Device>(flux_factor);
          }

          if (fr_as_crse) {
            fr_as_crse->CrseAdd(
              mfi, {D_DECL(&flux_ec[0], &flux_ec[1], &flux_ec[2])}, dxD.data(),
              dt, vfrac[mfi],
              {D_DECL(
                &((*areafrac[0])[mfi]), &((*areafrac[1])[mfi]),
                &((*areafrac[2])[mfi]))},
              device);
            if (AMREX_SPACEDIM <= 2) {
              amrex::Print()
                << "WARNING:Re redistribution crseadd for EB not tested "
                   "in 2D\n";
            }
          }

          if (fr_as_fine) {
            fr_as_fine->FineAdd(
              mfi, {D_DECL(&flux_ec[0], &flux_ec[1], &flux_ec[2])}, dxD.data(),
              dt, vfrac[mfi],
              {D_DECL(
                &((*areafrac[0])[mfi]), &((*areafrac[1])[mfi]),
                &((*areafrac[2])[mfi]))},
              dm_as_fine, device);

            if (AMREX_SPACEDIM <= 2) {
              amrex::Print()
                << "WARNING:Re redistribution fineadd for EB not tested "
                   "in 2D\n";
            }
          }
        }
      } else if (typ != amrex::FabType::regular) { // Single valued if loop
        amrex::Abort("multi-valued eb boundary fluxes to be implemented");
      }
#endif

      copy_array4(vbox, NVAR, Dterm, MOLSrc);

#ifdef PELEC_USE_EB
      // do regular flux reg ops
      if (do_reflux && flux_factor != 0 && typ == amrex::FabType::regular)
#else
      if (do_reflux && flux_factor != 0) // no eb in problem
#endif
      {
        for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
          amrex::ParallelFor(
            eboxes[dir], NVAR,
            [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
              flx[dir](i, j, k, n) *= flux_factor;
            });
        }

        if (level < parent->finestLevel()) {
          getFluxReg(level + 1).CrseAdd(
            mfi, {AMREX_D_DECL(&flux_ec[0], &flux_ec[1], &flux_ec[2])},
            dxD.data(), dt, device);
        }

        if (level > 0) {
          getFluxReg(level).FineAdd(
            mfi, {AMREX_D_DECL(&flux_ec[0], &flux_ec[1], &flux_ec[2])},
            dxD.data(), dt, device);
        }
      }

#ifdef PELEC_USE_EB
      if (do_mol_load_balance) {
        amrex::Gpu::streamSynchronize();
        wt = (amrex::ParallelDescriptor::second() - wt) / vbox.d_numPts();
        (*cost)[mfi].plus<amrex::RunOn::Device>(wt, vbox);
      }
#endif

      // Extrapolate to GhostCells
      if (MOLSrcTerm.nGrow() > 0) {
        BL_PROFILE("PeleC::diffextrap()");
        const int mg = MOLSrcTerm.nGrow();
        const amrex::Box bx = mfi.tilebox();
        auto low = bx.loVect();
        auto high = bx.hiVect();
        auto dlo = Dterm.begin;
        auto dhi = Dterm.end;
        const int AMREX_D_DECL(lx = low[0], ly = low[1], lz = low[2]);
        const int AMREX_D_DECL(hx = high[0], hy = high[1], hz = high[2]);
        amrex::ParallelFor(
          bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            pc_diffextrap(
              i, j, k, Dterm, mg, UMX, UMZ + 1, AMREX_D_DECL(lx, ly, lz),
              AMREX_D_DECL(hx, hy, hz), dlo, dhi);
            pc_diffextrap(
              i, j, k, Dterm, mg, UFS, UFS + NUM_SPECIES,
              AMREX_D_DECL(lx, ly, lz), AMREX_D_DECL(hx, hy, hz), dlo, dhi);
            pc_diffextrap(
              i, j, k, Dterm, mg, UEDEN, UEDEN + 1, AMREX_D_DECL(lx, ly, lz),
              AMREX_D_DECL(hx, hy, hz), dlo, dhi);
          });
      }
    } // End of MFIter scope
  }   // End of OMP scope
} // End of Function
