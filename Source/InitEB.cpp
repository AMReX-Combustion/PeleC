#include <memory>

#include "EB.H"
#include "prob.H"
#include "Utilities.H"
#include "Geometry.H"

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
#include <thrust/unique.h>
#include <thrust/sort.h>
#include <thrust/execution_policy.h>
#endif

inline bool
PeleC::ebInitialized()
{
  return eb_initialized;
}

void
PeleC::init_eb(
  const amrex::Geometry& /*level_geom*/,
  const amrex::BoxArray& /*ba*/,
  const amrex::DistributionMapping& /*dm*/)
{
  // Build the geometry information; this is done for each new set of grids
  initialize_eb2_structs();
}

// Set up PeleC EB Datastructures from AMReX EB2 constructs
// At the end of this routine, the following structures are populated:
//  - MultiFAB vfrac
//  - sv_eb_bndry_geom

void
PeleC::initialize_eb2_structs()
{
  BL_PROFILE("PeleC::initialize_eb2_structs()");
  amrex::Print() << "Initializing EB2 structs" << std::endl;

  static_assert(
    std::is_standard_layout<EBBndryGeom>::value,
    "EBBndryGeom is not standard layout");

  const amrex::MultiCutFab* bndrycent;
  std::array<const amrex::MultiCutFab*, AMREX_SPACEDIM> eb2areafrac;
  std::array<const amrex::MultiCutFab*, AMREX_SPACEDIM> facecent;

  const auto& ebfactory =
    dynamic_cast<amrex::EBFArrayBoxFactory const&>(Factory());

  // These are the data sources
  vfrac.clear();
  vfrac.define(grids, dmap, 1, numGrow(), amrex::MFInfo(), Factory());
  amrex::MultiFab::Copy(vfrac, ebfactory.getVolFrac(), 0, 0, 1, numGrow());
  bndrycent = &(ebfactory.getBndryCent());
  eb2areafrac = ebfactory.getAreaFrac();
  facecent = ebfactory.getFaceCent();

  // First pass over fabs to fill sparse per cut-cell ebg structures
  sv_eb_bndry_geom.resize(vfrac.local_size());
  sv_eb_bndry_grad_stencil.resize(vfrac.local_size());
  sv_eb_flux.resize(vfrac.local_size());
  sv_eb_bcval.resize(vfrac.local_size());

  auto const& flags = ebfactory.getMultiEBCellFlagFab();

  // Boundary stencil option: 0 = original, 1 = amrex way, 2 = least squares
  amrex::ParmParse pp("ebd");

  int bgs = 0;
  pp.query("boundary_grad_stencil_type", bgs);

  if (bgs == 0) {
    amrex::Print() << "Using quadratic stencil for the EB gradient\n";
  } else if (bgs == 1) {
    amrex::Print() << "Using least-squares stencil for the EB gradient\n";
  } else {
    amrex::Print() << "Unknown or unspecified EB gradient stencil type:" << bgs
                   << std::endl;
    amrex::Abort();
  }

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  for (amrex::MFIter mfi(vfrac, false); mfi.isValid(); ++mfi) {
    const amrex::Box tbox = mfi.growntilebox();
    const amrex::EBCellFlagFab& flagfab = flags[mfi];

    amrex::FabType typ = flagfab.getType(tbox);
    int iLocal = mfi.LocalIndex();

    if ((typ == amrex::FabType::regular) || (typ == amrex::FabType::covered)) {
      // do nothing
    } else if (typ == amrex::FabType::singlevalued) {
      const int Ncut = flagfab.getNumCutCells(tbox);
      sv_eb_bndry_geom[iLocal].resize(Ncut);
      auto const& flag_arr = flags.const_array(mfi);
      EBBndryGeom* d_sv_eb_bndry_geom = sv_eb_bndry_geom[iLocal].data();
      const auto lo = amrex::lbound(tbox);
      const auto hi = amrex::ubound(tbox);
      amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE(int /*dummy*/) noexcept {
        int ivec = 0;
        for (int i = lo.x; i <= hi.x; ++i) {
          for (int j = lo.y; j <= hi.y; ++j) {
            for (int k = lo.z; k <= hi.z; ++k) {
              const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
              const amrex::EBCellFlag& flag = flag_arr(iv);
              if (!(flag.isRegular() || flag.isCovered())) {
                d_sv_eb_bndry_geom[ivec].iv = iv;
                ivec++;
              }
            }
          }
        }
      });

      // int Nebg = sv_eb_bndry_geom[iLocal].size();

      // Now fill the sv_eb_bndry_geom
      auto const& vfrac_arr = vfrac.array(mfi);
      auto const& bndrycent_arr = bndrycent->array(mfi);
      AMREX_D_TERM(auto const& eb2areafrac_arr_0 = eb2areafrac[0]->array(mfi);
                   , auto const& eb2areafrac_arr_1 = eb2areafrac[1]->array(mfi);
                   ,
                   auto const& eb2areafrac_arr_2 = eb2areafrac[2]->array(mfi);)
      pc_fill_sv_ebg(
        tbox, Ncut, vfrac_arr, bndrycent_arr,
        AMREX_D_DECL(eb2areafrac_arr_0, eb2areafrac_arr_1, eb2areafrac_arr_2),
        sv_eb_bndry_geom[iLocal].data());

      sv_eb_bndry_grad_stencil[iLocal].resize(Ncut);

      // Fill in boundary gradient for cut cells in this grown tile
      const amrex::Real dx = geom.CellSize()[0];
#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
      const int sv_eb_bndry_geom_size = sv_eb_bndry_geom[iLocal].size();
      thrust::sort(
        thrust::device, sv_eb_bndry_geom[iLocal].data(),
        sv_eb_bndry_geom[iLocal].data() + sv_eb_bndry_geom_size,
        EBBndryGeomCmp());
#else
      sort<amrex::Gpu::DeviceVector<EBBndryGeom>>(sv_eb_bndry_geom[iLocal]);
#endif

      if (bgs == 0) {
        pc_fill_bndry_grad_stencil_quadratic(
          tbox, dx, Ncut, sv_eb_bndry_geom[iLocal].data(), Ncut,
          sv_eb_bndry_grad_stencil[iLocal].data());
      } else if (bgs == 1) {
        pc_fill_bndry_grad_stencil_ls(
          tbox, dx, Ncut, sv_eb_bndry_geom[iLocal].data(), Ncut,
          flags.array(mfi), sv_eb_bndry_grad_stencil[iLocal].data());
      } else {
        amrex::Print()
          << "Unknown or unspecified boundary gradient stencil type:" << bgs
          << std::endl;
        amrex::Abort();
      }

      sv_eb_flux[iLocal].define(sv_eb_bndry_grad_stencil[iLocal], NVAR);
      sv_eb_bcval[iLocal].define(sv_eb_bndry_grad_stencil[iLocal], QVAR);

      if (eb_isothermal && (diffuse_temp != 0 || diffuse_enth != 0)) {
        sv_eb_bcval[iLocal].setVal(eb_boundary_T, QTEMP);
      }
      if (eb_noslip && diffuse_vel == 1) {
        sv_eb_bcval[iLocal].setVal(0, QU, AMREX_SPACEDIM);
      }

    } else {
      amrex::Print() << "unknown (or multivalued) fab type" << std::endl;
      amrex::Abort();
    }
  }

  // Second pass over dirs and fabs to fill flux interpolation stencils
  amrex::Box fbox[AMREX_SPACEDIM];

  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    flux_interp_stencil[dir].resize(vfrac.local_size());

    fbox[dir] = amrex::bdryLo(
      amrex::Box(
        amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
        amrex::IntVect(AMREX_D_DECL(0, 0, 0))),
      dir, 1);

    for (int dir1 = 0; dir1 < AMREX_SPACEDIM; ++dir1) {
      if (dir1 != dir) {
        fbox[dir].grow(dir1, 1);
      }
    }

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(vfrac, false); mfi.isValid(); ++mfi) {
      const amrex::Box tbox = mfi.growntilebox(numGrow());
      const auto& flagfab = flags[mfi];
      amrex::FabType typ = flagfab.getType(tbox);
      int iLocal = mfi.LocalIndex();

      if (typ == amrex::FabType::regular || typ == amrex::FabType::covered) {
      } else if (typ == amrex::FabType::singlevalued) {
        const auto afrac_arr = (*eb2areafrac[dir])[mfi].array();
        const auto facecent_arr = (*facecent[dir])[mfi].array();

        // This used to be an std::set for cut_faces (it ensured
        // sorting and uniqueness)
        EBBndryGeom* d_sv_eb_bndry_geom = sv_eb_bndry_geom[iLocal].data();
        const int Nall_cut_faces = amrex::Reduce::Sum<int>(
          sv_eb_bndry_geom[iLocal].size(),
          [=] AMREX_GPU_DEVICE(int i) noexcept -> int {
            int r = 0;
            const amrex::IntVect& iv = d_sv_eb_bndry_geom[i].iv;
            for (int iside = 0; iside <= 1; iside++) {
              const amrex::IntVect iv_face = iv + iside * amrex::BASISV(dir);
              if (afrac_arr(iv_face) < 1.0) {
                r++;
              }
            }
            return r;
          });

        amrex::Gpu::DeviceVector<amrex::IntVect> v_all_cut_faces(
          Nall_cut_faces);
        amrex::IntVect* all_cut_faces = v_all_cut_faces.data();

        const int sv_eb_bndry_geom_size = sv_eb_bndry_geom[iLocal].size();
        // Serial loop on the GPU
        amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE(int /*dummy*/) {
          int cnt = 0;
          for (int i = 0; i < sv_eb_bndry_geom_size; i++) {
            const amrex::IntVect& iv = d_sv_eb_bndry_geom[i].iv;
            for (int iside = 0; iside <= 1; iside++) {
              const amrex::IntVect iv_face = iv + iside * amrex::BASISV(dir);
              if (afrac_arr(iv_face) < 1.0) {
                all_cut_faces[cnt] = iv_face;
                cnt++;
              }
            }
          }
        });

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
        const int v_all_cut_faces_size = v_all_cut_faces.size();
        thrust::sort(
          thrust::device, v_all_cut_faces.data(),
          v_all_cut_faces.data() + v_all_cut_faces_size);
        amrex::IntVect* unique_result_end = thrust::unique(
          thrust::device, v_all_cut_faces.data(),
          v_all_cut_faces.data() + v_all_cut_faces_size,
          thrust::equal_to<amrex::IntVect>());
        const int count_result =
          thrust::distance(v_all_cut_faces.data(), unique_result_end);
        amrex::Gpu::DeviceVector<amrex::IntVect> v_cut_faces(count_result);
        amrex::IntVect* d_all_cut_faces = v_all_cut_faces.data();
        amrex::IntVect* d_cut_faces = v_cut_faces.data();
        amrex::ParallelFor(
          v_cut_faces.size(), [=] AMREX_GPU_DEVICE(int i) noexcept {
            d_cut_faces[i] = d_all_cut_faces[i];
          });
#else
        sort<amrex::Gpu::DeviceVector<amrex::IntVect>>(v_all_cut_faces);
        amrex::Gpu::DeviceVector<amrex::IntVect> v_cut_faces =
          unique<amrex::Gpu::DeviceVector<amrex::IntVect>>(v_all_cut_faces);
#endif

        const int Nsten = v_cut_faces.size();
        if (Nsten > 0) {
          flux_interp_stencil[dir][iLocal].resize(Nsten);

          amrex::IntVect* cut_faces = v_cut_faces.data();
          auto* d_flux_interp_stencil = flux_interp_stencil[dir][iLocal].data();
          amrex::ParallelFor(
            v_cut_faces.size(), [=] AMREX_GPU_DEVICE(int i) noexcept {
              d_flux_interp_stencil[i].iv = cut_faces[i];
            });

          pc_fill_flux_interp_stencil(
            tbox, fbox[dir], Nsten, facecent_arr, afrac_arr,
            flux_interp_stencil[dir][iLocal].data());
        }
      } else {
        amrex::Abort("multi-valued flux interp stencil to be implemented");
      }
    }
  }
}

void
PeleC::define_body_state()
{
  BL_PROFILE("PeleC::define_body_state()");

  if (!eb_in_domain) {
    return;
  }

  // Scan over data and find a point in the fluid to use to
  // set computable values in cells outside the domain
  if (!body_state_set) {
    bool foundPt = false;
    const amrex::MultiFab& S = get_new_data(State_Type);
    auto const& fact =
      dynamic_cast<amrex::EBFArrayBoxFactory const&>(Factory());
    auto const& flags = fact.getMultiEBCellFlagFab();

    for (amrex::MFIter mfi(S, false); mfi.isValid() && !foundPt; ++mfi) {
      const amrex::Box vbox = mfi.validbox();
      auto const& farr = S.const_array(mfi);
      auto const& flag_arr = flags.const_array(mfi);

      const auto lo = amrex::lbound(vbox);
      const auto hi = amrex::ubound(vbox);
      amrex::Gpu::DeviceVector<amrex::Real> local_body_state(NVAR, -1);
      amrex::Real* p_body_state = local_body_state.begin();
      amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE(int /*dummy*/) noexcept {
        bool found = false;
        for (int i = lo.x; i <= hi.x && !found; ++i) {
          for (int j = lo.y; j <= hi.y && !found; ++j) {
            for (int k = lo.z; k <= hi.z && !found; ++k) {
              const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
              if (flag_arr(iv).isRegular()) {
                found = true;
                for (int n = 0; n < NVAR; ++n) {
                  p_body_state[n] = farr(iv, n);
                }
              }
            }
          }
        }
      });
      amrex::Gpu::copy(
        amrex::Gpu::deviceToHost, local_body_state.begin(),
        local_body_state.end(), body_state.begin());
      foundPt = body_state[0] != -1;
    }

    // Find proc with lowest rank to find valid point, use that for all
    amrex::Vector<int> found(amrex::ParallelDescriptor::NProcs(), 0);
    found[amrex::ParallelDescriptor::MyProc()] = (int)foundPt;
    amrex::ParallelDescriptor::ReduceIntSum(&(found[0]), found.size());
    int body_rank = -1;
    for (int i = 0; i < found.size(); ++i) {
      if (found[i] == 1) {
        body_rank = i;
      }
    }
    AMREX_ASSERT(body_rank >= 0);
    amrex::ParallelDescriptor::Bcast(
      &(body_state[0]), body_state.size(), body_rank); // NOLINT
    body_state_set = true;
  }
}

void
PeleC::set_body_state(amrex::MultiFab& S)
{
  BL_PROFILE("PeleC::set_body_state()");

  if (!eb_in_domain) {
    return;
  }

  if (!body_state_set) {
    define_body_state();
  }

  auto const& fact = dynamic_cast<amrex::EBFArrayBoxFactory const&>(Factory());
  auto const& flags = fact.getMultiEBCellFlagFab();

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  for (amrex::MFIter mfi(S, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    const amrex::Box& vbox = mfi.tilebox();
    auto const& Sar = S.array(mfi);
    auto const& flag_arr = flags.const_array(mfi);
    auto const captured_body_state = body_state;
    amrex::ParallelFor(
      vbox, NVAR, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
        pc_set_body_state(i, j, k, n, flag_arr, captured_body_state, Sar);
      });
  }
}

void
PeleC::zero_in_body(amrex::MultiFab& S) const
{
  BL_PROFILE("PeleC::zero_in_body()");

  if (!eb_in_domain) {
    return;
  }

  amrex::GpuArray<amrex::Real, NVAR> zeros = {0.0};
  auto const& fact = dynamic_cast<amrex::EBFArrayBoxFactory const&>(Factory());
  auto const& flags = fact.getMultiEBCellFlagFab();

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  for (amrex::MFIter mfi(S, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    const amrex::Box& vbox = mfi.tilebox();
    auto const& Sar = S.array(mfi);
    auto const& flag_arr = flags.const_array(mfi);
    amrex::ParallelFor(
      vbox, NVAR, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
        pc_set_body_state(i, j, k, n, flag_arr, zeros, Sar);
      });
  }
}

// Sets up implicit function using EB2 infrastructure
void
initialize_EB2(
  const amrex::Geometry& geom,
  const int
#ifdef LinePistonCylinder
    required_level
#endif
  ,
  const int max_level)
{
  BL_PROFILE("PeleC::initialize_EB2()");

  amrex::Print() << "Initializing EB2" << std::endl;
  amrex::ParmParse ppeb2("eb2");

  std::string geom_type("all_regular");
  ppeb2.query("geom_type", geom_type);

  int max_coarsening_level = 0;
  amrex::ParmParse ppamr("amr");
  amrex::Vector<int> ref_ratio(max_level, 2);
  ppamr.queryarr("ref_ratio", ref_ratio, 0, max_level);
  for (int lev = 0; lev < max_level; ++lev) {
    max_coarsening_level +=
      (ref_ratio[lev] == 2 ? 1
                           : 2); // Since EB always coarsening by factor of 2
  }

  // Custom types defined here - all_regular, plane, sphere, etc, will get
  // picked up by default (see AMReX_EB2.cpp around L100 )
  amrex::Vector<std::string> amrex_defaults(
    {"all_regular", "box", "cylinder", "plane", "sphere", "torus", "parser"});
  if (!(std::find(amrex_defaults.begin(), amrex_defaults.end(), geom_type) !=
        amrex_defaults.end())) {
    std::unique_ptr<pele::pelec::Geometry> geometry(
      pele::pelec::Geometry::create(geom_type));
    geometry->build(geom, max_coarsening_level);
  } else {
    amrex::EB2::Build(geom, max_level, max_level);
  }
}
