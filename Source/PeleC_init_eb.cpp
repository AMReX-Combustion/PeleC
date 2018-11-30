#include "PeleC.H"
#include "PeleC_F.H"

using namespace amrex;
#ifdef PELE_USE_EB

#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_BoxArray.H>
#include "AMReX_VisMF.H"
#include "AMReX_PlotFileUtil.H"

#if BL_SPACEDIM > 1
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Union.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Complement.H>
#include <AMReX_EB2_IF_Scale.H>
#include <AMReX_EB2_IF_Translation.H>
#include <AMReX_EB2_IF_Lathe.H>
#include <AMReX_EB2_IF_Box.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Ellipsoid.H>
#include <AMReX_EB2_IF_Sphere.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_GeometryShop.H>
#endif

#include "PeleC_init_eb_F.H"

inline
bool PeleC::ebInitialized()
{
    return eb_initialized;
}

void
PeleC::init_eb (const Geometry& level_geom, const BoxArray& ba, const DistributionMapping& dm)
{
  // Build the geometry information; this is done for each new set of grids
  initialize_eb2_structs();

}

#if BL_SPACEDIM > 1

/**
 * Set up PeleC EB Datastructures from AMReX EB2 constructs
 *
 * At the end of this routine, the following structures are populated:
 *   - FabArray ebmask
 *  - MultiFAB vfrac
 *  - sv_eb_bndry_geom
 */

void
PeleC::initialize_eb2_structs() {
  BL_PROFILE("PeleC::initialize_eb2_structs()");
  amrex::Print() << "Initializing EB2 structs" << std::endl;

  // n.b., could set this to 1 if geometry is all_regular as an optimization
  no_eb_in_domain = 0;

  //  1->regular, 0->irregular, -1->covered, 2->outside
  ebmask.define(grids, dmap,  1, 0);

  static_assert(std::is_standard_layout<EBBndryGeom>::value,
                "EBBndryGeom is not standard layout");

  const amrex::MultiFab* volfrac;
  const amrex::MultiCutFab* bndrycent;
  std::array<const amrex::MultiCutFab*, AMREX_SPACEDIM> eb2areafrac;
  std::array<const amrex::MultiCutFab*, AMREX_SPACEDIM> facecent;

  const auto& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(Factory());

  // These are the data sources
  volfrac = &(ebfactory.getVolFrac());
  bndrycent = &(ebfactory.getBndryCent());
  eb2areafrac = ebfactory.getAreaFrac();
  facecent = ebfactory.getFaceCent();

  vfrac.copy(*volfrac);

  // First pass over fabs to fill sparse per cut-cell ebg structures
  sv_eb_bndry_geom.resize(vfrac.local_size());
  sv_eb_bndry_grad_stencil.resize(vfrac.local_size());
  sv_eb_flux.resize(vfrac.local_size());
  sv_eb_bcval.resize(vfrac.local_size());

  auto const& flags = ebfactory.getMultiEBCellFlagFab();

  for (MFIter mfi(vfrac, false); mfi.isValid(); ++mfi) {
    BaseFab<int>& mfab = ebmask[mfi];
    const Box tbox = mfi.growntilebox();
    const FArrayBox& vfab = vfrac[mfi];
    const EBCellFlagFab& flagfab = flags[mfi];

    FabType typ = flagfab.getType(tbox);
    int iLocal = mfi.LocalIndex();

    if (typ == FabType::regular) {
      mfab.setVal(1);
    } else if (typ == FabType::covered) {
      mfab.setVal(-1);
    } else if (typ == FabType::singlevalued) {
      int Ncut = 0;
      for (BoxIterator bit(tbox); bit.ok(); ++bit) {
        const EBCellFlag& flag = flagfab(bit(), 0);
        if (!(flag.isRegular() || flag.isCovered())) {
          Ncut++;
        }
      }

      sv_eb_bndry_geom[iLocal].resize(Ncut);
      int ivec = 0;
      for (BoxIterator bit(tbox); bit.ok(); ++bit) {
        const EBCellFlag& flag = flagfab(bit(), 0);
        if (!(flag.isRegular() || flag.isCovered())) {
          EBBndryGeom& sv_ebg = sv_eb_bndry_geom[iLocal][ivec];
          ivec++;
          sv_ebg.iv = bit();
          if (mfab.box().contains(bit())) mfab(bit()) = 0;
        } else {
          if (flag.isRegular()) {
            if (mfab.box().contains(bit())) mfab(bit()) = 1;
          } else if (flag.isCovered()) {
            if (mfab.box().contains(bit())) mfab(bit()) = -1;
          } else {
            if (mfab.box().contains(bit())) mfab(bit()) = 2;
          }
        }
      }

      int Nebg = sv_eb_bndry_geom[iLocal].size();

      // Now call fortran to fill the ebg
      pc_fill_sv_ebg(BL_TO_FORTRAN_BOX(tbox),
                     sv_eb_bndry_geom[iLocal].data(), &Ncut,
                     BL_TO_FORTRAN_ANYD((*volfrac)[mfi]),
                     BL_TO_FORTRAN_ANYD((*bndrycent)[mfi]),
                     D_DECL(BL_TO_FORTRAN_ANYD((*eb2areafrac[0])[mfi]),
                            BL_TO_FORTRAN_ANYD((*eb2areafrac[1])[mfi]),
                            BL_TO_FORTRAN_ANYD((*eb2areafrac[2])[mfi])));

      sv_eb_bndry_grad_stencil[iLocal].resize(Ncut);

      // Fill in boundary gradient for cut cells in this grown tile
      const Real dx = geom.CellSize()[0];
      auto& vec = sv_eb_bndry_geom[iLocal];
      std::sort(vec.begin(), vec.end());
      pc_fill_bndry_grad_stencil(BL_TO_FORTRAN_BOX(tbox),
                                 sv_eb_bndry_geom[iLocal].data(), &Ncut,
                                 sv_eb_bndry_grad_stencil[iLocal].data(),
                                 &Ncut, &dx);

      sv_eb_flux[iLocal].define(sv_eb_bndry_grad_stencil[iLocal], NUM_STATE);
      sv_eb_bcval[iLocal].define(sv_eb_bndry_grad_stencil[iLocal], QVAR);

    } else {
      amrex::Print() << "unknown (or multivalued) fab type" << std::endl;
      amrex::Abort();
    }
  }

  // Second pass over dirs and fabs to fill flux interpolation stencils
  Box fbox[BL_SPACEDIM];

  for (int idir=0; idir < BL_SPACEDIM; ++idir) {
    flux_interp_stencil[idir].resize(vfrac.local_size());


    fbox[idir] = amrex::bdryLo(Box(IntVect(D_DECL(0, 0, 0)),
                                   IntVect(D_DECL(0, 0, 0))), idir, 1);

    for (int idir1=0; idir1 < BL_SPACEDIM; ++idir1) {
      if (idir1 != idir) fbox[idir].grow(idir1, 1);
    }

    for (MFIter  mfi(vfrac, false); mfi.isValid(); ++mfi) {
      const Box tbox = mfi.growntilebox(nGrowTr);
      const FArrayBox& vfab = vfrac[mfi];
      const auto& flagfab = flags[mfi];
      FabType typ = flagfab.getType(tbox);
      int iLocal = mfi.LocalIndex();

      if (typ == FabType::regular || typ == FabType::covered) {
      } else if (typ == FabType::singlevalued) {
        const Box ebox = Box(tbox).surroundingNodes(idir);
        const CutFab&  afrac_fab = (*eb2areafrac[idir])[mfi];
        const CutFab&  facecent_fab = (*facecent[idir])[mfi];

        std::set<IntVect> cut_faces;

        for (auto& sv_ebg : sv_eb_bndry_geom[iLocal]) {
          const IntVect& iv = sv_ebg.iv;
          for (int iside=0; iside <= 1; iside++) {
            const IntVect iv_face = iv + iside*BASISV(idir);
            if (afrac_fab(iv_face) < 1.0) {
              cut_faces.insert(iv_face);
            }
          }
        }

        int Nsten = cut_faces.size();
        if (Nsten > 0) {
          flux_interp_stencil[idir][iLocal].resize(Nsten);
          int ivec = 0;
          for (std::set<IntVect>::const_iterator it = cut_faces.begin();
                it != cut_faces.end(); ++it, ++ivec) {
            flux_interp_stencil[idir][iLocal][ivec].iv = *it;
          }

          pc_fill_flux_interp_stencil(BL_TO_FORTRAN_BOX(tbox),
                                      BL_TO_FORTRAN_BOX(fbox[idir]),
                                      flux_interp_stencil[idir][iLocal].data(),
                                      &Nsten, &idir,
                                      BL_TO_FORTRAN_ANYD(facecent_fab),
                                      BL_TO_FORTRAN_ANYD(afrac_fab));
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

  if (no_eb_in_domain) return;

  // Scan over data and find a point in the fluid to use to 
  // set computable values in cells outside the domain
  if (!body_state_set)
  {
    bool foundPt = false;
    const MultiFab& S = get_new_data(State_Type);
    BL_ASSERT(S.boxArray() == ebmask.boxArray());
    BL_ASSERT(S.DistributionMap() == ebmask.DistributionMap());
	
    body_state.resize(S.nComp(),0);
    for (MFIter mfi(S,false); mfi.isValid() && !foundPt; ++mfi)
    {
      const Box vbox = mfi.validbox();
      const BaseFab<int>& m = ebmask[mfi];
      const FArrayBox& fab = S[mfi];
      BL_ASSERT(m.box().contains(vbox));

      // TODO: Remove this dog and do this work in fortran 
      for (BoxIterator bit(vbox); bit.ok() && !foundPt; ++bit)
      {
        const IntVect& iv = bit();
        if (m(iv,0) == 1) {
          foundPt = true;
          for (int n=0; n<S.nComp(); ++n)
          {
            body_state[n] = fab(iv,n);
          }
        }
      }
    }

    // Find proc with lowest rank to find valid point, use that for all
    std::vector<int> found(ParallelDescriptor::NProcs(),0);
    found[ParallelDescriptor::MyProc()] = (int)foundPt;
    ParallelDescriptor::ReduceIntSum(&(found[0]),found.size());
    int body_rank = -1;
    for (int i=0; i<found.size(); ++i) {
      if (found[i]==1) {
        body_rank = i;
      }
    }
    BL_ASSERT(body_rank>=0);
    ParallelDescriptor::Bcast(&(body_state[0]),body_state.size(),body_rank);
    body_state_set = true;
  }
}

void
PeleC::set_body_state(MultiFab& S)
{
  BL_PROFILE("PeleC::set_body_state()");

  if (no_eb_in_domain) return;

  if (!body_state_set)
  {
    define_body_state();
  }

  BL_ASSERT(S.nComp() == body_state.size());
  int nc = S.nComp();
  int covered_val = -1;

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(S,true); mfi.isValid(); ++mfi)
  {
    const Box& vbox = mfi.validbox();
    pc_set_body_state(vbox.loVect(), vbox.hiVect(),
                      BL_TO_FORTRAN_ANYD(S[mfi]),
                      BL_TO_FORTRAN_ANYD(ebmask[mfi]),
                      &(body_state[0]),&nc,&covered_val);
  }
}

void
PeleC::zero_in_body(MultiFab& S) const
{
  BL_PROFILE("PeleC::zero_in_body()");

  if (no_eb_in_domain) return;

  int nc=S.nComp();
  Vector<Real> zeros(nc,0);
  int covered_val = -1;

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(S,true); mfi.isValid(); ++mfi)
  {
    const Box& vbox = mfi.validbox();
    pc_set_body_state(vbox.loVect(), vbox.hiVect(),
                      BL_TO_FORTRAN_ANYD(S[mfi]),
                      BL_TO_FORTRAN_ANYD(ebmask[mfi]),
                      &(body_state[0]),&nc,&covered_val);
  }
}

std::string
convertIntGG(int number)
{
  std::stringstream ss;//create a stringstream
  ss << number;//add number to the stream
  return ss.str();//return a string with the contents of the stream
}

std::string
PeleC::convertIntGG(int number)
{
  std::stringstream ss;//create a stringstream
  ss << number;//add number to the stream
  return ss.str();//return a string with the contents of the stream
}

/**
 *
 * Sets up implicit function using EB2 infrastructure
 *
 **/
void
initialize_EB2 (const Geometry& geom, const int required_level, const int max_level) {
  BL_PROFILE("initializeEB2");

  amrex::Print() << "Initializing EB2" << std::endl;
  ParmParse ppeb2("eb2");

  std::string geom_type("all_regular");
  ppeb2.query("geom_type", geom_type);

  int max_coarsening_level = max_level; // Because there are no mg solvers here

  // Custom types defined here - all_regular, plane, sphere, etc, will get picked up by default
  // (see AMReX_BE2.cpp around L100 )
  if (geom_type == "flat_plate") {
    amrex::Print() << "flat plate  geometry not currently supported. \n";
    amrex::Abort();
  } else if (geom_type == "ramp") {
    amrex::Print() << "ramp geometry\n";
    int upDir;
    int indepVar;
    Real startPt;
    Real slope;
    ppeb2.get("up_dir",upDir);
    ppeb2.get("indep_var",indepVar);
    ppeb2.get("start_pt", startPt);
    ppeb2.get("ramp_slope", slope);

    RealArray normal; normal.fill(0.0);
    normal[upDir] = 1.0;
    normal[indepVar] = -slope;

    RealArray point; point.fill(0.0);
    point[upDir] = -slope*startPt;

    bool normalInside = true;

    EB2::PlaneIF ramp(point, normal);
    auto gshop = EB2::makeShop(ramp);
    EB2::Build(gshop, geom, max_level, max_level);
  }
#if BL_SPACEDIM > 2
  else if (geom_type == "combustor")
  {
    ParmParse pp("combustor");
        
    Real fwl; 
    pp.get("far_wall_loc",fwl);

    EB2::PlaneIF farwall({AMREX_D_DECL(fwl,0.,0.)},
                         {AMREX_D_DECL(1. ,0.,0.)});
        
    Vector<Real> pl1pt, pl2pt, pl2nm, pl3pt; 
    pp.getarr("ramp_plane1_point", pl1pt);
    pp.getarr("ramp_plane2_point", pl2pt);
    pp.getarr("ramp_plane2_normal", pl2nm);
    pp.getarr("ramp_plane3_point", pl3pt);

    auto ramp = EB2::makeIntersection(EB2::PlaneIF({pl1pt[0], pl1pt[1], 0.},
                                                   {      0.,      -1., 0.}),
                                      EB2::PlaneIF({pl2pt[0], pl2pt[1], 0.},
                                                   {pl2nm[0], pl2nm[1], 0.}),
                                      EB2::PlaneIF({pl3pt[0], pl3pt[1], 0.},
                                                   {      1.,       0., 0.}));

    Vector<Real> pipelo, pipehi; 
    pp.getarr("pipe_lo", pipelo);
    pp.getarr("pipe_hi", pipehi);
        
    EB2::BoxIF pipe({pipelo[0], pipelo[1], -1.}, {pipehi[0], pipehi[1], 1.}, false);

    // where does plane 1 and plane 2 intersect?
    Real k2 = std::abs(pl2nm[0]/pl2nm[1]);
    Real secty = pl2pt[1] + k2*(pl3pt[0]-pl2pt[0]);
    // How much do we cut?
    Real dx = geom.CellSize(0);
    Real dycut = 4.*(1.+max_coarsening_level)*std::min(dx, k2*dx);
    EB2::BoxIF flat_corner({pl3pt[0], 0., -1.}, {1.e10, secty+dycut, 1.}, false);
        
    auto polys = EB2::makeUnion(farwall, ramp, pipe, flat_corner);

    Real lenx = Geometry::ProbLength(0);
    Real leny = Geometry::ProbLength(1);
    auto pr = EB2::translate(EB2::lathe(polys), {lenx*0.5, leny*0.5, 0.});
        
    auto gshop = EB2::makeShop(pr);
    EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level);
  }
  else if (geom_type == "triangles")
  {
      int num_pts;
      int max_tri=5;
      ParmParse pp("triangles");
      pp.get("num_tri", num_tri);
      Vector<Array<Real,3>> alltri(max_tri);

      //initalize all triangles with some dummy values
      //that fall outside of the domain
      

      for(int itri=0;itri<max_tri;itri++)
      {
         
         alltri[itri][0] =  
      }


      for(int itri = 0; itri < num_tri; itri++)
      {
          Array<Real,3> point{0.0,0.0,0.0};
          std::string    pointstr = "poly_point_" + convertIntGG(ipt);
          Vector<Real> vecpt;
          pp.getarr(pointstr.c_str(), vecpt,  0, SpaceDim);
          for(int idir = 0; idir < SpaceDim; idir++)
          {
              point[idir] = vecpt[idir] ;
          }
          polygon[ipt] = point;
      }

      Vector<EB2::PlaneIF> planes;
      planes.resize(0);

      int numPts = polygon.size();

      for (int n = 0; n < numPts; n++)
      {
          Array<Real,3> normal{0.0,0.0,0.0};
          Array<Real,3> point{0.0,0.0,0.0};

          normal[0] = -(polygon[(n+1) % numPts][1] - polygon[n][1]);
          normal[1] =  (polygon[(n+1) % numPts][0] - polygon[n][0]);

          point = polygon[n];

          EB2::PlaneIF plane(point,normal);

          planes.push_back(plane);
      }

      // Intersect all the half planes/spaces to create an implicit function
      // that represents the polygon
      auto polygonIF = EB2::makeIntersection(planes[0],planes[1]);

      for (int n=2; n < numPts; n++)
      {
          polygonIF = EB2::makeIntersection(polygonIF,planes[n])
      }
    
      auto gshop = EB2::makeShop(polygonIF);
      EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level);
  }

#endif
  else if (geom_type == "polygon_revolution")
  {
      amrex::Print() << "polygon_revolution  geometry not currently supported. combustor?\n";
      amrex::Abort();
      amrex::Print() << "creating geometry from polygon surfaces of revolution" << std::endl;
      bool insideRegular = false;

      // Data for polygons making up nozzle
      Vector<Vector<RealArray> > polygons;
      // For building each polygon - unlike original PeleEB, don't scale by domain size
      int num_poly;
      ppeb2.get("num_poly", num_poly);
      polygons.resize(num_poly);
      for(int ipoly = 0; ipoly < num_poly; ipoly++) {
          std::string nptsstr = "poly_" + convertIntGG(ipoly) + "_num_pts";
          int num_pts;
          ppeb2.get(nptsstr.c_str(), num_pts);
          Vector<RealArray> polygon(num_pts);
          for(int ipt = 0; ipt < num_pts; ipt++)
          {
              RealArray point; point.fill(0.0);
              std::string pointstr = "poly_" + convertIntGG(ipoly) + "_point_" + convertIntGG(ipt);
              Vector<Real> vecpt;
              ppeb2.getarr(pointstr.c_str(), vecpt,  0, SpaceDim);
              for(int idir = 0; idir < SpaceDim; idir++)
              {
                  point[idir] = vecpt[idir] ;
              }
              polygon[ipt] = point;
          }
          polygons[ipoly] = polygon;
      }
      Vector<EB2::PlaneIF> planes;
      for (int ipoly = 0; ipoly < num_poly; ipoly++) {
          const Vector<RealArray>& polygon = polygons[ipoly];
          int numPts = polygon.size(); // Number of pts in this polygon
          for (int n = 0; n < numPts; n++) {
              // The normal and point is space used to specify each half plane/space
              RealArray normal; normal.fill(0.0);
              RealArray point;

              // Set the normal remembering that the last point connects to the first
              // point.
              normal[0] = -(polygon[(n+1) % numPts][1] - polygon[n][1]);
              normal[1] =  (polygon[(n+1) % numPts][0] - polygon[n][0]);

              point = polygon[n];
              EB2::PlaneIF plane(point, normal);
              planes.push_back(plane);
          }

      }

      //PolygonIF pf(planes);
      //auto gshop = EB2::makeShop(pf);
      //EB2::Build(gshop, geom, max_level, max_level);


  } else {
      EB2::Build(geom, max_level, max_level);
      // Need to somehow set no_eb_in_domain for this level?

  }

}
#endif
#endif

