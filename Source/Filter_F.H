
#ifndef _FILTER_F_H_
#define _FILTER_F_H_
#include <AMReX_BLFort.H>

#include <AMReX.H>

extern "C"
{
  void filter(const int* lo, const int* hi,
              const BL_FORT_FAB_ARG_3D(Q),
              BL_FORT_FAB_ARG_3D(Qh),
              const int* ngrow,
              const amrex::Real* weights,
              const int* nstart,
              const int* ncnt,
              const int* ncomp);
}

#endif /*_FILTER_F_H_*/
