
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <string>

#include <AMReX_REAL.H>
#include <AMReX_Utility.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_ParallelDescriptor.H>

#ifdef _OPENM
#include <omp.h>
#endif

using namespace std;
using namespace amrex;


extern "C" void getplane(int* filename, int* len, Real* data, int* plane, int* ncomp, int* isswirltype);

void
getplane (int* filename, int* len, Real* data, int* plane, int* ncomp, int* isswirltype)
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
    static Vector< Vector<long> > offset(nthreads);
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
            amrex::FileOpenFailed(hdr);

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

        offset[tid].resize(kmax*AMREX_SPACEDIM,0);
           
        for (int i = 0; i < offset[tid].size(); i++)
            ifs >> offset[tid][i];
    }

    std::string dat = flctfile; dat += "/DAT";

    std::ifstream ifs;

    ifs.open(dat.c_str(), std::ios::in);

    if (!ifs.good())
        amrex::FileOpenFailed(dat);
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
        amrex::Abort("getplane(): seekg() failed");

    FArrayBox fab;

    fab.readFrom(ifs);

    memcpy(data, fab.dataPtr(), fab.box().numPts()*sizeof(Real));
}
