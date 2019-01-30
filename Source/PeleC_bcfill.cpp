
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_PhysBCFunct.H>
#include "PeleC.H"
#include "PeleC_F.H"

using namespace amrex;

// bx                  : Cells outside physical domain and inside bx are filled.
// data, dcomp, numcomp: Fill numcomp components of data starting from dcomp.
// bcr, bcomp          : bcr[bcomp] specifies BC for component dcomp and so on.
// scomp               : component index for dcomp as in the desciptor set up in CNS::variableSetUp.

void pc_bcfill_hyp (Box const& bx, FArrayBox& data,
                 const int dcomp, const int numcomp,
                 Geometry const& geom, const Real time,
                 const Vector<BCRec>& bcr, const int bcomp,
                 const int scomp)
{
//AMREX_ALWAYS_ASSERT();

const int* domlo      = geom.Domain().loVect();
const int* domhi      = geom.Domain().hiVect();
const Real* dx  = geom.CellSize();

const RealBox pbx  = RealBox(bx,geom.CellSize(),geom.ProbLo());
const Real* xlo     = pbx.lo();

const int* bc =  bcr[bcomp].data();

const amrex::Real& time_f=time;

pc_hypfill(BL_TO_FORTRAN_3D(data),ARLIM_3D(domlo),ARLIM_3D(domhi),
           ZFILL(dx),ZFILL(xlo),&time_f,bc);



}


void pc_reactfill_hyp (Box const& bx, FArrayBox& data,
                 const int dcomp, const int numcomp,
                 Geometry const& geom, const Real time,
                 const Vector<BCRec>& bcr, const int bcomp,
                 const int scomp)
{
//AMREX_ALWAYS_ASSERT();

const int* domlo      = geom.Domain().loVect();
const int* domhi      = geom.Domain().hiVect();
const Real* dx  = geom.CellSize();

const RealBox pbx  = RealBox(bx,geom.CellSize(),geom.ProbLo());
const Real* xlo     = pbx.lo();

const int* bc =  bcr[bcomp].data();

const amrex::Real& time_f=time;

#ifdef REACTIONS
pc_reactfill(BL_TO_FORTRAN_3D(data),ARLIM_3D(domlo),ARLIM_3D(domhi),
           ZFILL(dx),ZFILL(xlo),&time_f,bc);
#endif


}

