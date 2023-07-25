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
    BL_PROFILE("PeleC::construct_hydro_source()");

    if ((verbose != 0) && amrex::ParallelDescriptor::IOProcessor()) {
      amrex::Print() << "... Computing hydro advance" << std::endl;
    }

    AMREX_ASSERT(S.nGrow() >= numGrow() + nGrowF);

    const int ng = 0;
    sources_for_hydro.setVal(0.0);
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

    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx =
      geom.CellSizeArray();

    amrex::Real courno = std::numeric_limits<amrex::Real>::lowest();

    const amrex::MultiFab& S_new = get_new_data(State_Type);

    auto const& fact =
      dynamic_cast<amrex::EBFArrayBoxFactory const&>(S.Factory());
    auto const& flags = fact.getMultiEBCellFlagFab();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion()) reduction(max \
                                                                    : courno)
#endif
    {
      amrex::Real cflLoc = std::numeric_limits<amrex::Real>::lowest();

      const int* domain_lo = geom.Domain().loVect();
      const int* domain_hi = geom.Domain().hiVect();

      for (amrex::MFIter mfi(S_new, amrex::TilingIfNotGPU()); mfi.isValid();
           ++mfi) {

        const amrex::Box& bx = mfi.tilebox();
        const amrex::Box& qbx = amrex::grow(bx, numGrow() + nGrowF);
        const amrex::Box& fbx = amrex::grow(bx, nGrowF);

        const auto& flag_fab = flags[mfi];
        amrex::FabType typ = flag_fab.getType(bx);
        const amrex::Array4<amrex::EBCellFlag const>& flag_arr =
          flag_fab.const_array();

        if (typ == amrex::FabType::covered) {
          continue;
        }

        amrex::GpuArray<amrex::FArrayBox, AMREX_SPACEDIM> flux;
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
          const amrex::Box& efbx = surroundingNodes(fbx, dir);
          flux[dir].resize(efbx, NVAR, amrex::The_Async_Arena());
          flux[dir].setVal<amrex::RunOn::Device>(0.0);
        }

        auto const& sarr = S.array(mfi);
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
              if (!flag_arr(i, j, k).isCovered()) {
                pc_ctoprim(i, j, k, sarr, qarr, qauxar);
              } else {
                for (int n = 0; n < QVAR; n++) {
                  qarr(i, j, k, n) = 0.0;
                }
              }
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
              if (!flag_arr(i, j, k).isCovered()) {
                pc_srctoprim(i, j, k, qarr, qauxar, src_in, srcqarr);
              } else {
                for (int n = 0; n < QVAR; n++) {
                  srcqarr(i, j, k, n) = 0.0;
                }
              }
            });
        }
        const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM>
          flx_arr{
            {AMREX_D_DECL(flux[0].array(), flux[1].array(), flux[2].array())}};
        const amrex::GpuArray<
          const amrex::Array4<const amrex::Real>, AMREX_SPACEDIM>
          a{{AMREX_D_DECL(
            area[0].const_array(mfi), area[1].const_array(mfi),
            area[2].const_array(mfi))}};

        if (typ == amrex::FabType::singlevalued) {
          amrex::Abort("EB Godunov not implemented yet");
        } else if (typ == amrex::FabType::regular) {
          BL_PROFILE("PeleC::umdrv()");
          pc_umdrv(
            time, fbx, domain_lo, domain_hi, phys_bc.lo(), phys_bc.hi(), sarr,
            hyd_src, qarr, qauxar, srcqarr, dx, dt, ppm_type, use_flattening,
            use_hybrid_weno, weno_scheme, difmag, flx_arr, a, volume.array(mfi),
            cflLoc);
        } else if (typ == amrex::FabType::multivalued) {
          amrex::Abort("multi-valued cells are not supported");
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

        // Reflux
        if (do_reflux && sub_iteration == sub_ncycle - 1) {
          const amrex::FabType gtyp = flag_fab.getType(amrex::grow(bx, 1));
          update_flux_registers(
            dt, bx, mfi, gtyp,
            {{AMREX_D_DECL(flux.data(), &(flux[1]), &(flux[2]))}});
        }
      }
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
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx,
  const amrex::Real dt,
  const int ppm_type,
  const bool use_flattening,
  const bool use_hybrid_weno,
  const int weno_scheme,
  const amrex::Real difmag,
  const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx,
  const amrex::GpuArray<const amrex::Array4<const amrex::Real>, AMREX_SPACEDIM>&
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
  amrex::ParallelFor(bxg2, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    pc_divu(i, j, k, q, AMREX_D_DECL(dx[0], dx[1], dx[2]), divuarr);
  });

  pc_adjust_fluxes(bx, uin, flx, a, divuarr, dx, difmag);

  // consup
  pc_consup(bx, uout, flx, vol, pdivuarr);
}

void
pc_adjust_fluxes(
  const amrex::Box& bx,
  const amrex::Array4<const amrex::Real>& u,
  const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx,
  const amrex::GpuArray<const amrex::Array4<const amrex::Real>, AMREX_SPACEDIM>&
    a,
  const amrex::Array4<const amrex::Real>& divu,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& del,
  const amrex::Real difmag)
{
  BL_PROFILE("PeleC::pc_adjust_fluxes()");
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
}

void
pc_consup(
  const amrex::Box& bx,
  const amrex::Array4<amrex::Real>& update,
  const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx,
  const amrex::Array4<const amrex::Real>& vol,
  const amrex::Array4<const amrex::Real>& pdivu)
{
  BL_PROFILE("PeleC::pc_consup()");
  // Combine for Hydro Sources
  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    pc_update(i, j, k, update, flx, vol, pdivu);
  });
}
