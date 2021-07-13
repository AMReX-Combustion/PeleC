#include <memory>

#include "EB.H"
#include "prob.H"
#include "Utilities.H"

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
//   - FabArray ebmask
//  - MultiFAB vfrac
//  - sv_eb_bndry_geom

void
PeleC::initialize_eb2_structs()
{
  BL_PROFILE("PeleC::initialize_eb2_structs()");
  amrex::Print() << "Initializing EB2 structs" << std::endl;

  // NOTE: THIS NEEDS TO BE REPLACED WITH A FLAGFAB

  // 1->regular, 0->irregular, -1->covered, 2->outside
  ebmask.define(grids, dmap, 1, 0);

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
    amrex::Print() << "Using quadratic stencil\n";
  } else if (bgs == 1) {
    amrex::Print() << "Using least-squares stencil\n";
  } else {
    amrex::Print() << "Unknown or unspecified boundary gradient stencil type:"
                   << bgs << std::endl;
    amrex::Abort();
  }

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  for (amrex::MFIter mfi(vfrac, false); mfi.isValid(); ++mfi) {
    amrex::BaseFab<int>& mfab = ebmask[mfi];
    const amrex::Box tbox = mfi.growntilebox();
    const amrex::EBCellFlagFab& flagfab = flags[mfi];

    amrex::FabType typ = flagfab.getType(tbox);
    int iLocal = mfi.LocalIndex();

    if (typ == amrex::FabType::regular) {
      mfab.setVal<amrex::RunOn::Device>(1);
    } else if (typ == amrex::FabType::covered) {
      mfab.setVal<amrex::RunOn::Device>(-1);
    } else if (typ == amrex::FabType::singlevalued) {
      int Ncut = 0;
      for (amrex::BoxIterator bit(tbox); bit.ok(); ++bit) {
        const amrex::EBCellFlag& flag = flagfab(bit(), 0);

        if (!(flag.isRegular() || flag.isCovered())) {
          Ncut++;
        }
      }

      sv_eb_bndry_geom[iLocal].resize(Ncut);
      int ivec = 0;
      for (amrex::BoxIterator bit(tbox); bit.ok(); ++bit) {
        const amrex::EBCellFlag& flag = flagfab(bit(), 0);

        if (!(flag.isRegular() || flag.isCovered())) {
          EBBndryGeom* d_sv_eb_bndry_geom = sv_eb_bndry_geom[iLocal].data();
          amrex::IntVect captured_bit = bit();
          // Serial loop on the GPU
          amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE(int /*dummy*/) {
            d_sv_eb_bndry_geom[ivec].iv = captured_bit;
          });
          ivec++;

          if (mfab.box().contains(bit())) {
            mfab(bit()) = 0;
          }
        } else {
          if (flag.isRegular()) {
            if (mfab.box().contains(bit())) {
              mfab(bit()) = 1;
            }
          } else if (flag.isCovered()) {
            if (mfab.box().contains(bit())) {
              mfab(bit()) = -1;
            }
          } else {
            if (mfab.box().contains(bit())) {
              mfab(bit()) = 2;
            }
          }
        }
      }

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
          v_all_cut_faces.data(), v_all_cut_faces.data() + v_all_cut_faces_size,
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
    AMREX_ASSERT(S.boxArray() == ebmask.boxArray());
    AMREX_ASSERT(S.DistributionMap() == ebmask.DistributionMap());

    for (amrex::MFIter mfi(S, false); mfi.isValid() && !foundPt; ++mfi) {
      const amrex::Box vbox = mfi.validbox();
      const amrex::BaseFab<int>& m = ebmask[mfi];
      const amrex::FArrayBox& fab = S[mfi];
      AMREX_ASSERT(m.box().contains(vbox));

      for (amrex::BoxIterator bit(vbox); bit.ok() && !foundPt; ++bit) {
        const amrex::IntVect& iv = bit();
        if (m(iv, 0) == 1) {
          foundPt = true;
          for (int n = 0; n < S.nComp(); ++n) {
            body_state[n] = fab(iv, n);
          }
        }
      }
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
      &(body_state[0]), body_state.size(), body_rank);
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

  int covered_val = -1;

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  for (amrex::MFIter mfi(S, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    const amrex::Box& vbox = mfi.tilebox();
    auto const& Sar = S.array(mfi);
    auto const& mask = ebmask.array(mfi);
    auto const captured_body_state = body_state;
    amrex::ParallelFor(
      vbox, NVAR, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
        pc_set_body_state(
          i, j, k, n, mask, captured_body_state, covered_val, Sar);
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
  int covered_val = -1;

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  for (amrex::MFIter mfi(S, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    const amrex::Box& vbox = mfi.tilebox();
    auto const& Sar = S.array(mfi);
    auto const& mask = ebmask.array(mfi);
    amrex::ParallelFor(
      vbox, NVAR, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
        pc_set_body_state(i, j, k, n, mask, zeros, covered_val, Sar);
      });
  }
}

std::string
convertIntGG(int number)
{
  std::stringstream ss; // create a stringstream
  ss << number;         // add number to the stream
  return ss.str();      // return a string with the contents of the stream
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

  int max_coarsening_level = max_level; // Because there are no mg solvers here

  // Custom types defined here - all_regular, plane, sphere, etc, will get
  // picked up by default (see AMReX_BE2.cpp around L100 )
  if (geom_type == "flat_plate") {
    amrex::Print() << "flat plate  geometry not currently supported. \n";
    amrex::Abort();
  } else if (geom_type == "ramp") {
    amrex::Print() << "ramp geometry\n";
    int upDir;
    int indepVar;
    amrex::Real startPt;
    amrex::Real slope;
    ppeb2.get("up_dir", upDir);
    ppeb2.get("indep_var", indepVar);
    ppeb2.get("start_pt", startPt);
    ppeb2.get("ramp_slope", slope);

    amrex::RealArray normal;
    normal.fill(0.0);
    normal[upDir] = 1.0;
    normal[indepVar] = -slope;

    amrex::RealArray point;
    point.fill(0.0);
    point[upDir] = -slope * startPt;

    // bool normalInside = true;

    amrex::EB2::PlaneIF ramp(point, normal);
    auto gshop = amrex::EB2::makeShop(ramp);
    amrex::EB2::Build(gshop, geom, max_level, max_level);
  } else if (geom_type == "combustor") {
    amrex::ParmParse pp("combustor");

    amrex::Real fwl;
    pp.get("far_wall_loc", fwl);

    amrex::EB2::PlaneIF farwall(
      {AMREX_D_DECL(fwl, 0., 0.)}, {AMREX_D_DECL(1., 0., 0.)});

    amrex::Vector<amrex::Real> pl1pt;
    amrex::Vector<amrex::Real> pl2pt;
    amrex::Vector<amrex::Real> pl2nm;
    amrex::Vector<amrex::Real> pl3pt;
    pp.getarr("ramp_plane1_point", pl1pt);
    pp.getarr("ramp_plane2_point", pl2pt);
    pp.getarr("ramp_plane2_normal", pl2nm);
    pp.getarr("ramp_plane3_point", pl3pt);

    amrex::EB2::PlaneIF r0(
      {AMREX_D_DECL(pl1pt[0], pl1pt[1], 0.)}, {AMREX_D_DECL(0., -1., 0.)});
    amrex::EB2::PlaneIF r1(
      {AMREX_D_DECL(pl2pt[0], pl2pt[1], 0.)},
      {AMREX_D_DECL(pl2nm[0], pl2nm[1], 0.)});
    amrex::EB2::PlaneIF r2(
      {AMREX_D_DECL(pl3pt[0], pl3pt[1], 0.)}, {AMREX_D_DECL(1., 0., 0.)});
    auto ramp = amrex::EB2::makeIntersection(r0, r1, r2);

    amrex::Vector<amrex::Real> pipelo;
    amrex::Vector<amrex::Real> pipehi;
    pp.getarr("pipe_lo", pipelo);
    pp.getarr("pipe_hi", pipehi);

    amrex::EB2::BoxIF pipe(
      {AMREX_D_DECL(pipelo[0], pipelo[1], -1.)},
      {AMREX_D_DECL(pipehi[0], pipehi[1], 1.)}, false);

    // where does plane 1 and plane 2 intersect?
    amrex::Real k2 = amrex::Math::abs(pl2nm[0] / pl2nm[1]);
    amrex::Real secty = pl2pt[1] + k2 * (pl3pt[0] - pl2pt[0]);
    // How much do we cut?
    amrex::Real dx = geom.CellSize(0);
    amrex::Real dycut =
      4. * (1. + max_coarsening_level) * std::min(dx, k2 * dx);
    amrex::EB2::BoxIF flat_corner(
      {AMREX_D_DECL(pl3pt[0], 0., -1.)},
      {AMREX_D_DECL(1.e10, secty + dycut, 1.)}, false);

    auto polys = amrex::EB2::makeUnion(farwall, ramp, pipe, flat_corner);

    amrex::Real lenx = amrex::DefaultGeometry().ProbLength(0);
    amrex::Real leny = amrex::DefaultGeometry().ProbLength(1);
    auto pr = amrex::EB2::translate(
      amrex::EB2::lathe(polys), {AMREX_D_DECL(lenx * 0.5, leny * 0.5, 0.)});

    auto gshop = amrex::EB2::makeShop(pr);
    amrex::EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level);
  } else if (geom_type == "ICE_PistonBowl") {
    // amrex::RealArray point;
    // amrex::RealArray normal;

    amrex::RealArray center;
    center[0] = 0.04 - 0.0125 - 0.02;
    center[1] = 0.0;
    center[2] = 0.0;

    amrex::Real radius;
    radius = 0.02;

    bool has_fluid_inside = false;
    amrex::EB2::SphereIF sf(radius, center, has_fluid_inside);

    amrex::EB2::CylinderIF cf1(
      0.04, 0.09, 0, {AMREX_D_DECL(0.045, 0.0, 0.0)}, true);
    amrex::EB2::CylinderIF cf2(
      0.04, 0.125, 0, {AMREX_D_DECL(-0.0125 - 0.02, 0.0, 0.0)}, false);
    amrex::EB2::CylinderIF cf3(
      0.03, 0.125, 0, {AMREX_D_DECL(-0.0125 - 0.02, 0.0, 0.0)}, false);
    auto pipe = amrex::EB2::makeDifference(cf2, cf3);
    amrex::EB2::CylinderIF cf4(
      0.03, 0.10, 0, {AMREX_D_DECL(-0.0125 - 0.02, 0.0, 0.0)}, false);

    amrex::RealArray center2;
    center2[0] = 0.0;
    center2[1] = 0.0;
    center2[2] = 0.0;

    amrex::Real radius2;
    radius2 = 0.09;

    bool has_fluid_inside2 = true;
    amrex::EB2::SphereIF sf2(radius2, center2, has_fluid_inside2);

    auto polys = amrex::EB2::makeUnion(cf1, pipe, cf4, sf, sf2);

    auto gshop = amrex::EB2::makeShop(polys);
    amrex::EB2::Build(
      gshop, geom, max_coarsening_level, max_coarsening_level, 4);
  } else if (geom_type == "extruded_triangles") {
    // setting some constants
    // the polygon is triangle
    // we can only do a maximum of 5 triangles (change if needed)
    const int npts_in_tri = 3;
    const int max_tri = 5;

    // number of user defined triangles
    int num_tri;

    amrex::ParmParse pp("extruded_triangles");
    amrex::Vector<amrex::Array<amrex::Real, AMREX_SPACEDIM>> alltri(
      npts_in_tri * max_tri);

    // initalize all triangles with some dummy values
    // that fall outside of the domain
    const amrex::Real* problo;
    const amrex::Real* probhi;
    amrex::Real maxlen;

    problo = geom.ProbLo();
    probhi = geom.ProbHi();

    maxlen = std::max(
      std::max(geom.ProbLength(0), geom.ProbLength(1)), geom.ProbLength(2));

    // setting all triangles to be waaay outside the domain initially
    for (int itri = 0; itri < max_tri; itri++) {
      alltri[npts_in_tri * itri + 0][0] = problo[0] + 100.0 * maxlen;
      alltri[npts_in_tri * itri + 0][1] = problo[1] - 100.0 * maxlen;
      alltri[npts_in_tri * itri + 0][2] = 0.0;

      alltri[npts_in_tri * itri + 1][0] = probhi[0] + 100.0 * maxlen;
      alltri[npts_in_tri * itri + 1][1] = problo[1] - 100.0 * maxlen;
      alltri[npts_in_tri * itri + 1][2] = 0.0;

      alltri[npts_in_tri * itri + 2][0] = probhi[0] + 100.0 * maxlen;
      alltri[npts_in_tri * itri + 2][1] = problo[1] + 100.0 * maxlen;
      alltri[npts_in_tri * itri + 2][2] = 0.0;
    }

    // get user defined number of triangles
    pp.get("num_tri", num_tri);

    for (int itri = 0; itri < num_tri; itri++) {
      amrex::Array<amrex::Real, AMREX_SPACEDIM> point{
        AMREX_D_DECL(0.0, 0.0, 0.0)};

      for (int ipt = 0; ipt < npts_in_tri; ipt++) {
        std::string pointstr =
          "tri_" + convertIntGG(itri) + "_point_" + convertIntGG(ipt);
        amrex::Vector<amrex::Real> vecpt;
        pp.getarr(pointstr.c_str(), vecpt, 0, AMREX_SPACEDIM);
        for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
          point[dir] = vecpt[dir];
        }
        alltri[npts_in_tri * itri + ipt] = point;
      }
    }

    // intersection of the 3 planes in a triangle for all triangles
    amrex::Vector<std::unique_ptr<amrex::EB2::IntersectionIF<
      amrex::EB2::PlaneIF, amrex::EB2::PlaneIF, amrex::EB2::PlaneIF>>>
      impfunc_triangles(max_tri);

    for (int itri = 0; itri < max_tri; itri++) {
      // make sure points are in anti clockwise direction to set the inside of
      // the triangle as solid phase correctly
      amrex::Array<amrex::Real, AMREX_SPACEDIM> norm0;
      amrex::Array<amrex::Real, AMREX_SPACEDIM> norm1;
      amrex::Array<amrex::Real, AMREX_SPACEDIM> norm2;

      amrex::Array<amrex::Real, AMREX_SPACEDIM> point0;
      amrex::Array<amrex::Real, AMREX_SPACEDIM> point1;
      amrex::Array<amrex::Real, AMREX_SPACEDIM> point2;

      point0 = alltri[npts_in_tri * itri + 0];
      point1 = alltri[npts_in_tri * itri + 1];
      point2 = alltri[npts_in_tri * itri + 2];

      norm0[0] = -(point1[1] - point0[1]);
      norm0[1] = (point1[0] - point0[0]);
      norm0[2] = 0.0;

      norm1[0] = -(point2[1] - point1[1]);
      norm1[1] = (point2[0] - point1[0]);
      norm1[2] = 0.0;

      norm2[0] = -(point0[1] - point2[1]);
      norm2[1] = (point0[0] - point2[0]);
      norm2[2] = 0.0;

      // normalize so that magnitude is 1
      amrex::Real norm = sqrt(norm0[0] * norm0[0] + norm0[1] * norm0[1]);
      norm0[0] = norm0[0] / norm;
      norm0[1] = norm0[1] / norm;

      // normalize so that magnitude is 1
      norm = sqrt(norm1[0] * norm1[0] + norm1[1] * norm1[1]);
      norm1[0] = norm1[0] / norm;
      norm1[1] = norm1[1] / norm;

      // normalize so that magnitude is 1
      norm = sqrt(norm2[0] * norm2[0] + norm2[1] * norm2[1]);
      norm2[0] = norm2[0] / norm;
      norm2[1] = norm2[1] / norm;

      amrex::EB2::PlaneIF plane0(point0, norm0);
      amrex::EB2::PlaneIF plane1(point1, norm1);
      amrex::EB2::PlaneIF plane2(point2, norm2);

      impfunc_triangles[itri] = std::make_unique<amrex::EB2::IntersectionIF<
        amrex::EB2::PlaneIF, amrex::EB2::PlaneIF, amrex::EB2::PlaneIF>>(

        plane0, plane1, plane2);
    }

    auto alltri_IF = amrex::EB2::makeUnion(
      *impfunc_triangles[0], *impfunc_triangles[1], *impfunc_triangles[2],
      *impfunc_triangles[3], *impfunc_triangles[4]);

    auto alltri_extrude_IF = amrex::EB2::extrude(alltri_IF, 2); // along z

    auto gshop = amrex::EB2::makeShop(alltri_extrude_IF);
    amrex::EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level);
  } else if (geom_type == "Line-Piston-Cylinder") {
#ifdef LinePistonCylinder
    EBLinePistonCylinder(geom, required_level, max_level);
#else
    amrex::Abort("Line-Piston-Cylinder geom_type not supported");
#endif
  } else if (geom_type == "polygon_revolution") {
    amrex::Print() << "polygon_revolution  geometry not currently supported. "
                      " combustor?\n";
    amrex::Abort();
    amrex::Print() << "creating geometry from polygon surfaces of revolution"
                   << std::endl;
    // bool insideRegular = false;

    // Data for polygons making up nozzle
    amrex::Vector<amrex::Vector<amrex::RealArray>> polygons;
    // For building each polygon - unlike original PeleEB, don't scale by
    // domain size
    int num_poly;
    ppeb2.get("num_poly", num_poly);
    polygons.resize(num_poly);
    for (int ipoly = 0; ipoly < num_poly; ipoly++) {
      std::string nptsstr = "poly_" + convertIntGG(ipoly) + "_num_pts";
      int num_pts;
      ppeb2.get(nptsstr.c_str(), num_pts);
      amrex::Vector<amrex::RealArray> polygon(num_pts);
      for (int ipt = 0; ipt < num_pts; ipt++) {
        amrex::RealArray point;
        point.fill(0.0);
        std::string pointstr =
          "poly_" + convertIntGG(ipoly) + "_point_" + convertIntGG(ipt);
        amrex::Vector<amrex::Real> vecpt;
        ppeb2.getarr(pointstr.c_str(), vecpt, 0, AMREX_SPACEDIM);
        for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
          point[dir] = vecpt[dir];
        }
        polygon[ipt] = point;
      }
      polygons[ipoly] = polygon;
    }
    amrex::Vector<amrex::EB2::PlaneIF> planes;
    for (int ipoly = 0; ipoly < num_poly; ipoly++) {
      const amrex::Vector<amrex::RealArray>& polygon = polygons[ipoly];
      int numPts = polygon.size(); // Number of pts in this polygon
      for (int n = 0; n < numPts; n++) {
        // The normal and point is space used to specify each half plane/space
        amrex::RealArray normal;
        normal.fill(0.0);
        amrex::RealArray point;

        // Set the normal remembering that the last point connects to the
        // first point
        normal[0] = -(polygon[(n + 1) % numPts][1] - polygon[n][1]);
        normal[1] = (polygon[(n + 1) % numPts][0] - polygon[n][0]);

        point = polygon[n];
        amrex::EB2::PlaneIF plane(point, normal);
        planes.push_back(plane);
      }
    }

    // PolygonIF pf(planes);
    // auto gshop = EB2::makeShop(pf);
    // EB2::Build(gshop, geom, max_level, max_level);
  } else if (geom_type == "moving_plane") {
    amrex::RealArray point;
    point[0] = 0.5;
    point[1] = 0.0;
    point[2] = 0.0;

    amrex::RealArray normal;
    normal[0] = -1.0;
    normal[1] = 0.0;
    normal[2] = 0.0;

    amrex::EB2::PlaneIF pf(point, normal);

    // amrex::EB2::GeometryShop<amrex::EB2::PlaneIF> gshop(pf);

    amrex::EB2::BoxIF pipe(
      {AMREX_D_DECL(-1.0, 0.25, -1.)}, {AMREX_D_DECL(1.5, 0.5, 1.)}, false);
    auto gshop = amrex::EB2::makeShop(pipe);

    amrex::EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level);
  } else if (geom_type == "quarter-circle") {

    amrex::Real r_inner = 1.0;
    amrex::Real r_outer = 2.0;
    ppeb2.query("r_inner", r_inner);
    ppeb2.query("r_outer", r_outer);

    amrex::EB2::CylinderIF inner(
      r_inner, 10, 2, {AMREX_D_DECL(0, 0, 0)}, false);
    amrex::EB2::CylinderIF outer(r_outer, 10, 2, {AMREX_D_DECL(0, 0, 0)}, true);

    auto polys = amrex::EB2::makeUnion(inner, outer);
    auto gshop = amrex::EB2::makeShop(polys);
    amrex::EB2::Build(
      gshop, geom, max_coarsening_level, max_coarsening_level, 4, false);
  } else if (geom_type == "converging-nozzle") {
#ifdef ConvergingNozzle
    EBConvergingNozzle(geom, max_level);
#else
    amrex::Abort("converging-nozzle geom_type not supported");
#endif
  } else {
    amrex::EB2::Build(geom, max_level, max_level);
  }
}
