
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <string>

#include <REAL.H>
#include <Utility.H>
#include <FArrayBox.H>
#include <ParallelDescriptor.H>
#include <PArray.H>

#ifdef _OPENM
#include <omp.h>
#endif

#  if defined(BL_FORT_USE_UPPERCASE)
#    define FORT_GETPLANE    GETPLANE
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_GETPLANE    getplane
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_GETPLANE    getplane_
#  endif

extern "C" void FORT_GETPLANE(int* filename, int* len, Real* data, int* plane, int* ncomp, int* isswirltype);

void
FORT_GETPLANE (int* filename, int* len, Real* data, int* plane, int* ncomp, int* isswirltype)
{
#ifdef _OPENMP
    int nthreads = omp_get_max_threads();
    int tid = omp_get_thread_num();
#else
    int nthreads = 1;
    int tid = 0;
#endif

    static int         kmax;
    static bool        first = true;
    static PArray< Array<long> > offset(nthreads, PArrayManage);
    std::string        flctfile;

#ifdef _OPENMP
#pragma omp threadprivate (kmax, first)
#endif


    for (int i = 0; i < *len; i++)
    {
        char c = filename[i];

        flctfile += c;
    }

    if (first)
    {
        //
        // Read and save all the seekp() offsets in the inflow header file.
        //
        first = false;

        std::string hdr = flctfile; hdr += "/HDR";

        std::ifstream ifs;

        ifs.open(hdr.c_str(), std::ios::in);

        if (!ifs.good())
            BoxLib::FileOpenFailed(hdr);

        int  idummy;
        Real rdummy;
        //
        // Hardwire loop max to 3 regardless of spacedim.
        //
        for (int i = 0; i < 3; i++)
            ifs >> kmax;

        ifs >> rdummy >> rdummy >> rdummy;
        ifs >> idummy >> idummy >> idummy;

        if (*isswirltype)
        {
            //
            // Skip over fluct_times array.
            //
            for (int i = 0; i < kmax; i++)
                ifs >> rdummy;
        }

	offset.set(tid, new Array<long>(kmax*BL_SPACEDIM, 0));

        for (int i = 0; i < offset[tid].size(); i++)
            ifs >> offset[tid][i];
    }

    std::string dat = flctfile; dat += "/DAT";

    std::ifstream ifs;

    ifs.open(dat.c_str(), std::ios::in);

    if (!ifs.good())
        BoxLib::FileOpenFailed(dat);
    //
    // There are BL_SPACEDIM * kmax planes of FABs.
    // The first component are in the first kmax planes,
    // the second component in the next kmax planes, ....
    // Note also that both (*plane) and (*ncomp) start from
    // 1 not 0 since they're passed from Fortran.
    //
    const long start = offset[tid][((*plane) - 1) + (((*ncomp) - 1) * kmax)] ;

    ifs.seekg(start, std::ios::beg);

    if (!ifs.good())
        BoxLib::Abort("getplane(): seekg() failed");

    FArrayBox fab;

    fab.readFrom(ifs);

    memcpy(data, fab.dataPtr(), fab.box().numPts()*sizeof(Real));
}
