#include "Hydro.H"

// Set up the source terms to go into the hydro.
void
PeleC::construct_hydro_source(
  const amrex::MultiFab& S,
  amrex::Real time,
  amrex::Real dt,
  int /*amr_iteration*/,
  int /*amr_ncycle*/,
  int sub_iteration,
  int sub_ncycle)
{
  if (do_mol) {
    if ((verbose != 0) && amrex::ParallelDescriptor::IOProcessor()) {
      amrex::Print() << "... Zeroing Godunov-based hydro advance" << std::endl;
    }
    hydro_source.setVal(0);
  } else {
    BL_PROFILE("PeleC::advance_hydro_pc_umdrv()");

    if ((verbose != 0) && amrex::ParallelDescriptor::IOProcessor()) {
      amrex::Print() << "... Computing hydro advance" << std::endl;
    }

    AMREX_ASSERT(S.nGrow() == numGrow() + nGrowF);
    sources_for_hydro.setVal(0.0);
    int ng = 0; // TODO: This is currently the largest ngrow of the source
                // data...maybe this needs fixing?
    for (int n = 0; n < src_list.size(); ++n) {
      amrex::MultiFab::Saxpy(
        sources_for_hydro, 0.5, *new_sources[src_list[n]], 0, 0, NVAR, ng);
      amrex::MultiFab::Saxpy(
        sources_for_hydro, 0.5, *old_sources[src_list[n]], 0, 0, NVAR, ng);
    }
    // Add I_R terms to advective forcing
    if (do_react) {
      amrex::MultiFab::Add(
        sources_for_hydro, get_new_data(Reactions_Type), 0, FirstSpec,
        NUM_SPECIES, ng);
      amrex::MultiFab::Add(
        sources_for_hydro, get_new_data(Reactions_Type), NUM_SPECIES, Eden, 1,
        ng);
    }
    sources_for_hydro.FillBoundary(geom.periodicity());
    hydro_source.setVal(0);

    int finest_level = parent->finestLevel();

    const amrex::Real* dx = geom.CellSize();

    amrex::Real dx1 = dx[0];
    for (int dir = 1; dir < AMREX_SPACEDIM; ++dir) {
      dx1 *= dx[dir];
    }

    std::array<amrex::Real, AMREX_SPACEDIM> dxD = {
      {AMREX_D_DECL(dx1, dx1, dx1)}};
    const amrex::Real* dxDp = &(dxD[0]);

    amrex::Real courno = std::numeric_limits<amrex::Real>::lowest();

    amrex::MultiFab& S_new = get_new_data(State_Type);

    // note: the radiation consup currently does not fill these
    amrex::Real E_added_flux = 0.;
    amrex::Real mass_added_flux = 0.;
    amrex::Real xmom_added_flux = 0.;
    amrex::Real ymom_added_flux = 0.;
    amrex::Real zmom_added_flux = 0.;
    amrex::Real mass_lost = 0.;
    amrex::Real xmom_lost = 0.;
    amrex::Real ymom_lost = 0.;
    amrex::Real zmom_lost = 0.;
    amrex::Real eden_lost = 0.;
    amrex::Real xang_lost = 0.;
    amrex::Real yang_lost = 0.;
    amrex::Real zang_lost = 0.;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())               \
    reduction(+:E_added_flux,mass_added_flux)                           \
    reduction(+:xmom_added_flux,ymom_added_flux,zmom_added_flux)	\
    reduction(+:mass_lost,xmom_lost,ymom_lost,zmom_lost)		\
    reduction(+:eden_lost,xang_lost,yang_lost,zang_lost) 		\
    reduction(max:courno)
#endif
    {
      // amrex::IArrayBox bcMask[AMREX_SPACEDIM];
      amrex::Real cflLoc = std::numeric_limits<amrex::Real>::lowest();
      int is_finest_level = (level == finest_level) ? 1 : 0;
      // int flag_nscbc_isAnyPerio = (geom.isAnyPeriodic()) ? 1 : 0;
      // int flag_nscbc_perio[AMREX_SPACEDIM] = {0}; // For 3D, we will know
      // which corners have a periodicity
      // for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      //  flag_nscbc_perio[dir] =
      //    (amrex::DefaultGeometry().isPeriodic(dir)) ? 1 : 0;
      // }

      const int* domain_lo = geom.Domain().loVect();
      const int* domain_hi = geom.Domain().hiVect();

      // Temporary Fabs needed for Hydro Computation
      for (amrex::MFIter mfi(S_new, amrex::TilingIfNotGPU()); mfi.isValid();
           ++mfi) {

        const amrex::Box& bx = mfi.tilebox();
        const amrex::Box& qbx = amrex::grow(bx, numGrow() + nGrowF);
        const amrex::Box& fbx = amrex::grow(bx, nGrowF);
        // const int* lo = bx.loVect();
        // const int* hi = bx.hiVect();

        amrex::GpuArray<amrex::FArrayBox, AMREX_SPACEDIM> flux;
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
          const amrex::Box& efbx = surroundingNodes(fbx, dir);
          flux[dir].resize(efbx, NVAR, amrex::The_Async_Arena());
        }

        auto const& s = S.array(mfi);
        auto const& hyd_src = hydro_source.array(mfi);

        // Resize Temporary Fabs
        amrex::FArrayBox q(qbx, QVAR, amrex::The_Async_Arena());
        amrex::FArrayBox qaux(qbx, NQAUX, amrex::The_Async_Arena());
        amrex::FArrayBox src_q(qbx, QVAR, amrex::The_Async_Arena());
        // Get Arrays to pass to the gpu.
        auto const& qarr = q.array();
        auto const& qauxar = qaux.array();
        auto const& srcqarr = src_q.array();

        {
          BL_PROFILE("PeleC::ctoprim()");
          amrex::ParallelFor(
            qbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              pc_ctoprim(i, j, k, s, qarr, qauxar);
            });
        }

        // TODO GPUize NCSCBC
        // Imposing Ghost-Cells Navier-Stokes Characteristic BCs if "UserBC" are
        // used For the theory, see Motheau et al. AIAA J. Vol. 55, No. 10 : pp.
        // 3399-3408, 2017.
        //
        // The user should provide a bcnormal routine in bc_fill_module with
        // additional optional arguments to temporary fill ghost-cells for
        // EXT_DIR and to provide target BC values. See the examples.

        // Allocate fabs for bcMask. Note that we grow in the opposite direction
        // because the Riemann solver wants a face value in a ghost-cell
        /*
              for (int dir = 0; dir < AMREX_SPACEDIM ; dir++)  {
                const Box& bxtmp = amrex::surroundingNodes(fbx,dir);
                Box TestBox(bxtmp);
                for(int d=0; d<AMREX_SPACEDIM; ++d) {
                  if (dir!=d) TestBox.grow(d,1);
                }
                bcMask[dir].resize(TestBox,1, amrex::The_Async_Arena());
                bcMask[dir].setVal(0);
              }

              // Becase bcMask is read in the Riemann solver in any case,
              // here we put physbc values in the appropriate faces for the
           non-nscbc case set_bc_mask(lo, hi, domain_lo, domain_hi,
                          AMREX_D_DECL(AMREX_TO_FORTRAN(bcMask[0]),
                                 AMREX_TO_FORTRAN(bcMask[1]),
                                 AMREX_TO_FORTRAN(bcMask[2])));

              if (nscbc_adv == 1)
              {
                impose_NSCBC(lo, hi, domain_lo, domain_hi,
                             AMREX_TO_FORTRAN(*statein),
                             AMREX_TO_FORTRAN(q.fab()),
                             AMREX_TO_FORTRAN(qaux.fab()),
                             AMREX_D_DECL(AMREX_TO_FORTRAN(bcMask[0]),
                                    AMREX_TO_FORTRAN(bcMask[1]),
                                    AMREX_TO_FORTRAN(bcMask[2])),
                             &flag_nscbc_isAnyPerio, flag_nscbc_perio,
                             &time, dx, &dt);
              }
        */
        {
          BL_PROFILE("PeleC::srctoprim()");
          const auto& src_in = sources_for_hydro.array(mfi);
          amrex::ParallelFor(
            qbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              pc_srctoprim(i, j, k, qarr, qauxar, src_in, srcqarr);
            });
        }

        amrex::FArrayBox pradial(
          amrex::Box::TheUnitBox(), 1, amrex::The_Async_Arena());
        if (!amrex::DefaultGeometry().IsCartesian()) {
          pradial.resize(
            amrex::surroundingNodes(bx, 0), 1, amrex::The_Async_Arena());
        }

        const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM>
          flx_arr{
            {AMREX_D_DECL(flux[0].array(), flux[1].array(), flux[2].array())}};
        const amrex::GpuArray<
          const amrex::Array4<const amrex::Real>, AMREX_SPACEDIM>
          a{{AMREX_D_DECL(
            area[0].array(mfi), area[1].array(mfi), area[2].array(mfi))}};
        {
          BL_PROFILE("PeleC::umdrv()");
          pc_umdrv(
            is_finest_level, time, fbx, domain_lo, domain_hi, phys_bc.lo(),
            phys_bc.hi(), s, hyd_src, qarr, qauxar, srcqarr, dx, dt, ppm_type,
            use_flattening, use_hybrid_weno, weno_scheme, difmag, flx_arr, a,
            volume.array(mfi), cflLoc);
        }

        courno = amrex::max<amrex::Real>(courno, cflLoc);

        // Filter hydro source and fluxes here
        if (use_explicit_filter) {
          BL_PROFILE("PeleC::apply_filter()");
          for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
            const amrex::Box& bxtmp = amrex::surroundingNodes(bx, dir);
            amrex::FArrayBox filtered_flux(
              bxtmp, NVAR, amrex::The_Async_Arena());
            les_filter.apply_filter(
              bxtmp, flux[dir], filtered_flux, Density, NVAR);

            setV(bxtmp, flux[dir].nComp(), flx_arr[dir], 0.0);
            copy_array4(
              bxtmp, flux[dir].nComp(), filtered_flux.array(), flx_arr[dir]);
          }

          amrex::FArrayBox filtered_source_out(
            bx, NVAR, amrex::The_Async_Arena());
          les_filter.apply_filter(
            bx, hydro_source[mfi], filtered_source_out, Density, NVAR);

          setV(bx, hyd_src.nComp(), hyd_src, 0.0);
          copy_array4(
            bx, hyd_src.nComp(), filtered_source_out.array(), hyd_src);
        }

        if (do_reflux && sub_iteration == sub_ncycle - 1) {
          BL_PROFILE("PeleC::reflux()");
          if (level < finest_level) {
            getFluxReg(level + 1).CrseAdd(
              mfi, {{AMREX_D_DECL(&(flux[0]), &(flux[1]), &(flux[2]))}}, dxDp,
              dt, amrex::RunOn::Device);

            if (!amrex::DefaultGeometry().IsCartesian()) {
              amrex::Abort("Flux registers not r-z compatible yet");
            }
          }

          if (level > 0) {
            getFluxReg(level).FineAdd(
              mfi, {{AMREX_D_DECL(&(flux[0]), &(flux[1]), &(flux[2]))}}, dxDp,
              dt, amrex::RunOn::Device);

            if (!amrex::DefaultGeometry().IsCartesian()) {
              amrex::Abort("Flux registers not r-z compatible yet");
            }
          }
        }
      }
    }

    if (track_grid_losses) {
      material_lost_through_boundary_temp[0] += mass_lost;
      material_lost_through_boundary_temp[1] += xmom_lost;
      material_lost_through_boundary_temp[2] += ymom_lost;
      material_lost_through_boundary_temp[3] += zmom_lost;
      material_lost_through_boundary_temp[4] += eden_lost;
      material_lost_through_boundary_temp[5] += xang_lost;
      material_lost_through_boundary_temp[6] += yang_lost;
      material_lost_through_boundary_temp[7] += zang_lost;
    }

    if (print_energy_diagnostics) {
      amrex::Real foo[5] = {
        E_added_flux, xmom_added_flux, ymom_added_flux, zmom_added_flux,
        mass_added_flux};

      amrex::ParallelDescriptor::ReduceRealSum(
        foo, 5, amrex::ParallelDescriptor::IOProcessorNumber());

#ifdef AMREX_DEBUG
      if (amrex::ParallelDescriptor::IOProcessor()) {
        E_added_flux = foo[0];
        xmom_added_flux = foo[1];
        ymom_added_flux = foo[2];
        zmom_added_flux = foo[3];
        mass_added_flux = foo[4];
        amrex::Print() << "mass added from fluxes                      : "
                       << mass_added_flux << std::endl;
        amrex::Print() << "xmom added from fluxes                      : "
                       << xmom_added_flux << std::endl;
        amrex::Print() << "ymom added from fluxes                      : "
                       << ymom_added_flux << std::endl;
        amrex::Print() << "zmom added from fluxes                      : "
                       << zmom_added_flux << std::endl;
        amrex::Print() << "(rho E) added from fluxes                   : "
                       << E_added_flux << std::endl;
      }
#endif
    }

    if (courno > 1.0) {
      amrex::Print() << "WARNING -- EFFECTIVE CFL AT THIS LEVEL " << level
                     << " IS " << courno << '\n';
      if (hard_cfl_limit) {
        amrex::Abort("CFL is too high at this level -- go back to a checkpoint "
                     "and restart with lower cfl number");
      }
    }
  }
}

void
pc_umdrv(
  const int /*is_finest_level*/,
  const amrex::Real /*time*/,
  amrex::Box const& bx,
  const int* domlo,
  const int* domhi,
  const int* bclo,
  const int* bchi,
  amrex::Array4<const amrex::Real> const& uin,
  amrex::Array4<amrex::Real> const& uout,
  amrex::Array4<const amrex::Real> const& q,
  amrex::Array4<const amrex::Real> const& qaux,
  amrex::Array4<const amrex::Real> const&
    src_q, // amrex::IArrayBox const& bcMask,
  const amrex::Real* dx,
  const amrex::Real dt,
  const int ppm_type,
  const bool use_flattening,
  const bool use_hybrid_weno,
  const int weno_scheme,
  const amrex::Real difmag,
  const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM> flx,
  const amrex::GpuArray<const amrex::Array4<const amrex::Real>, AMREX_SPACEDIM>
    a,
  amrex::Array4<amrex::Real> const& vol,
  amrex::Real /*cflLoc*/)
{
  // Set Up for Hydro Flux Calculations
  auto const& bxg2 = grow(bx, 2);
  amrex::FArrayBox qec[AMREX_SPACEDIM];
  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    const amrex::Box eboxes = amrex::surroundingNodes(bxg2, dir);
    qec[dir].resize(eboxes, NGDNV, amrex::The_Async_Arena());
  }
  amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM> qec_arr{
    {AMREX_D_DECL(qec[0].array(), qec[1].array(), qec[2].array())}};

  // Temporary FArrayBoxes
  amrex::FArrayBox divu(bxg2, 1, amrex::The_Async_Arena());
  amrex::FArrayBox pdivu(bx, 1, amrex::The_Async_Arena());
  auto const& divuarr = divu.array();
  auto const& pdivuarr = pdivu.array();

  {
    BL_PROFILE("PeleC::umeth()");
#if AMREX_SPACEDIM == 1
    amrex::Abort("PLM isn't implemented in 1D.");
    // pc_umeth_1D(
    //  bx, bclo, bchi, domlo, domhi, q, qaux, src_q, // bcMask,
    //  flx[0], qec_arr[0], a[0], pdivuarr, vol, dx, dt);
#elif AMREX_SPACEDIM == 2
    pc_umeth_2D(
      bx, bclo, bchi, domlo, domhi, q, qaux, src_q, // bcMask,
      flx[0], flx[1], qec_arr[0], qec_arr[1], a[0], a[1], pdivuarr, vol, dx, dt,
      ppm_type, use_flattening, use_hybrid_weno, weno_scheme);
#elif AMREX_SPACEDIM == 3
    pc_umeth_3D(
      bx, bclo, bchi, domlo, domhi, q, qaux, src_q, // bcMask,
      flx[0], flx[1], flx[2], qec_arr[0], qec_arr[1], qec_arr[2], a[0], a[1],
      a[2], pdivuarr, vol, dx, dt, ppm_type, use_flattening, use_hybrid_weno,
      weno_scheme);
#endif
  }

  // divu
  AMREX_D_TERM(const amrex::Real dx0 = dx[0];, const amrex::Real dx1 = dx[1];
               , const amrex::Real dx2 = dx[2];);
  amrex::ParallelFor(bxg2, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    pc_divu(i, j, k, q, AMREX_D_DECL(dx0, dx1, dx2), divuarr);
  });

  // consup
  pc_consup(bx, uin, uout, flx, a, vol, divuarr, pdivuarr, dx, difmag);
}

void
pc_consup(
  amrex::Box const& bx,
  amrex::Array4<const amrex::Real> const& u,
  amrex::Array4<amrex::Real> const& update,
  const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM> flx,
  const amrex::GpuArray<const amrex::Array4<const amrex::Real>, AMREX_SPACEDIM>
    a,
  amrex::Array4<const amrex::Real> const& vol,
  amrex::Array4<const amrex::Real> const& divu,
  amrex::Array4<const amrex::Real> const& pdivu,
  amrex::Real const* del,
  amrex::Real const difmag)
{
  BL_PROFILE("PeleC::consup()");
  // Flux alterations
  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    amrex::Box const& fbx = surroundingNodes(bx, dir);
    const amrex::Real dx = del[dir];
    amrex::ParallelFor(fbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      pc_artif_visc(AMREX_D_DECL(i, j, k), flx[dir], divu, u, dx, difmag, dir);
      // Normalize Species Flux
      pc_norm_spec_flx(i, j, k, flx[dir]);
      // Make flux extensive
      pc_ext_flx(i, j, k, flx[dir], a[dir]);
    });
  }

  // Combine for Hydro Sources
  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    pc_update(i, j, k, update, flx, vol, pdivu);
  });
}
