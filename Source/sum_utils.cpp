#include <iomanip>

#include <PeleC.H>
#include <PeleC_F.H>

using namespace std;
using namespace amrex;

Real
PeleC::sumDerive (const std::string& name,
		  Real               time,
		  bool               local)
{
#ifdef PELE_USE_EB
  amrex::Abort("sumDerive undefined for EB");
#endif

    Real sum     = 0.0;
    auto mf = derive(name, time, 0);

    BL_ASSERT(!(mf == 0));

    if (level < parent->finestLevel())
    {
	const MultiFab& mask = getLevel(level+1).build_fine_mask();
	MultiFab::Multiply(*mf, mask, 0, 0, 1, 0);
    }

#ifdef _OPENMP
#pragma omp parallel reduction(+:sum)
#endif
    {
	for (MFIter mfi(*mf,true); mfi.isValid(); ++mfi)
	{
	    sum += (*mf)[mfi].sum(mfi.tilebox(),0);
	}
    }


    if (!local)
	ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

Real
PeleC::volWgtSum (const std::string& name,
		  Real               time,
		  bool               local,
		  bool               finemask)
{
  BL_PROFILE("PeleC::volWgtSum()");

  Real        sum = 0.0;
  const Real* dx  = geom.CellSize();
  auto   mf       = derive(name,time,0);

  BL_ASSERT(mf != 0);

  if (level < parent->finestLevel() && finemask)
  {
    const MultiFab& mask = getLevel(level+1).build_fine_mask();
    MultiFab::Multiply(*mf, mask, 0, 0, 1, 0);
  }

#ifdef PELE_USE_EB
  auto const& fact = dynamic_cast<EBFArrayBoxFactory const&>(mf->Factory());
  auto const& flags = fact.getMultiEBCellFlagFab();
#endif

#ifdef _OPENMP
#pragma omp parallel reduction(+:sum)
#endif    
  for (MFIter mfi(*mf,true); mfi.isValid(); ++mfi)
  {
    Real s = 0.0;
    const Box& box  = mfi.tilebox();

    FArrayBox vol(box);
    vol.setVal(1.0);
#ifdef PELE_USE_EB
    const auto& flag_fab = flags[mfi];
    FabType typ = flag_fab.getType(box);
    if (typ == FabType::covered) {
      continue;
    }
    vol.copy(vfrac[mfi], box, 0, box, 0, 1);
#endif

    auto& fab = (*mf)[mfi];

    const int* lo   = box.loVect();
    const int* hi   = box.hiVect();

    vol.mult(volume[mfi]);

    //
    // Note that this routine will do a volume weighted sum of
    // whatever quantity is passed in, not strictly the "mass".
    //
    pc_summass(ARLIM_3D(lo),ARLIM_3D(hi),BL_TO_FORTRAN_3D(fab),
               ZFILL(dx),BL_TO_FORTRAN_3D(vol),&s);

    sum += s;
  }


  if (!local)
    ParallelDescriptor::ReduceRealSum(sum);

  return sum;
}

Real
PeleC::volWgtSquaredSum (const std::string& name,
			 Real               time,
			 bool               local)
{
    BL_PROFILE("PeleC::volWgtSquaredSum()");

    Real        sum     = 0.0;
    const Real* dx      = geom.CellSize();
    auto   mf      = derive(name,time,0);

    BL_ASSERT(mf != 0);

    if (level < parent->finestLevel())
    {
	const MultiFab& mask = getLevel(level+1).build_fine_mask();
	MultiFab::Multiply(*mf, mask, 0, 0, 1, 0);
    }

#ifdef PELE_USE_EB
    MultiFab::Multiply(*mf,vfrac,0,0,1,0);

    auto const& fact = dynamic_cast<EBFArrayBoxFactory const&>(mf->Factory());
    auto const& flags = fact.getMultiEBCellFlagFab();
#endif

#ifdef _OPENMP
#pragma omp parallel reduction(+:sum)
#endif
    for (MFIter mfi(*mf,true); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = (*mf)[mfi];

        Real s = 0.0;
        const Box& box  = mfi.tilebox();
        const int* lo   = box.loVect();
        const int* hi   = box.hiVect();

        FArrayBox vol(box);
        vol.setVal(1.0);
#ifdef PELE_USE_EB
        const auto& flag_fab = flags[mfi];
        FabType typ = flag_fab.getType(box);
        if (typ == FabType::covered) {
          continue;
        }
        vol.copy(vfrac[mfi], box, 0, box, 0, 1);
#endif
        vol.mult(volume[mfi]);
        //
        // Note that this routine will do a volume weighted sum of
        // whatever quantity is passed in, not strictly the "mass".
        //

	pc_sumsquared(ARLIM_3D(lo),ARLIM_3D(hi),BL_TO_FORTRAN_3D(fab),
		      ZFILL(dx),BL_TO_FORTRAN_3D(vol),&s);

        sum += s;
    }


    if (!local)
	ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

Real
PeleC::volWgtSumMF (MultiFab* mf, int comp, bool local, bool finemask) {
    BL_PROFILE("PeleC::volWgtSumMF()");

    Real        sum     = 0.0;
    const Real* dx      = geom.CellSize();

    BL_ASSERT(mf != 0);
    MultiFab* mask;
    if (level < parent->finestLevel() && finemask) {
        mask = &(getLevel(level+1).build_fine_mask());
    }
#ifdef PELE_USE_EB
    auto const& fact = dynamic_cast<EBFArrayBoxFactory const&>(mf->Factory());
    auto const& flags = fact.getMultiEBCellFlagFab();
#endif

#ifdef _OPENMP
#pragma omp parallel reduction(+:sum)
#endif
    for (MFIter mfi(*mf,true); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = (*mf)[mfi];

        Real s = 0.0;
        const Box& box  = mfi.tilebox();
        const int* lo   = box.loVect();
        const int* hi   = box.hiVect();
        FArrayBox vol(box);
        vol.setVal(1.0);

#ifdef PELE_USE_EB
        const auto& flag_fab = flags[mfi];
        FabType typ = flag_fab.getType(box);
        if (typ == FabType::covered) {
          continue;
        }
        vol.copy(vfrac[mfi], box, 0, box, 0, 1);
#endif
        vol.mult(volume[mfi]);
        if (level < parent->finestLevel() && finemask) {
            vol.mult((*mask)[mfi]);
        }

        //
        // Note that this routine will do a volume weighted sum of
        // whatever quantity is passed in, not strictly the "mass".
        //

	pc_summass(ARLIM_3D(lo),ARLIM_3D(hi),BL_TO_FORTRAN_N_3D(fab,comp),
		   ZFILL(dx),BL_TO_FORTRAN_3D(volume[mfi]),&s);

        sum += s;
    }

    if (!local)
	ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}
