#include <cstdio>

#include <AMReX_Arena.H>

#include "mechanism.h"
#include "chemistry_file.H"
#include "TransportParams.H"

extern "C" {
void egtransetWT(amrex::Real* wt);
void egtransetEPS(amrex::Real* eps);
void egtransetSIG(amrex::Real* sig);
void egtransetDIP(amrex::Real* dip);
void egtransetPOL(amrex::Real* pol);
void egtransetZROT(amrex::Real* zrot);
void egtransetNLIN(int* nlin);
void egtransetCOFETA(amrex::Real* fitmu);
void egtransetCOFLAM(amrex::Real* fitlam);
void egtransetCOFD(amrex::Real* fitdbin);
}

namespace transport_params {

AMREX_GPU_DEVICE_MANAGED amrex::Real* wt;
AMREX_GPU_DEVICE_MANAGED amrex::Real* iwt;
AMREX_GPU_DEVICE_MANAGED amrex::Real* eps;
AMREX_GPU_DEVICE_MANAGED amrex::Real* sig;
AMREX_GPU_DEVICE_MANAGED amrex::Real* dip;
AMREX_GPU_DEVICE_MANAGED amrex::Real* pol;
AMREX_GPU_DEVICE_MANAGED amrex::Real* zrot;
AMREX_GPU_DEVICE_MANAGED amrex::Real* fitmu;
AMREX_GPU_DEVICE_MANAGED amrex::Real* fitlam;
AMREX_GPU_DEVICE_MANAGED amrex::Real* fitdbin;
AMREX_GPU_DEVICE_MANAGED int* nlin;
AMREX_GPU_DEVICE_MANAGED int array_size;
AMREX_GPU_DEVICE_MANAGED int fit_length;

void
init()
{
  array_size = NUM_SPECIES;
  fit_length = NUM_FIT;
  //    std::cout << " array_size " << array_size << std::endl;
  //    std::cout << " fit_length " << fit_length << std::endl;
  wt = static_cast<amrex::Real*>(
    amrex::The_Managed_Arena()->alloc(sizeof(amrex::Real) * array_size));
  iwt = static_cast<amrex::Real*>(
    amrex::The_Managed_Arena()->alloc(sizeof(amrex::Real) * array_size));
  eps = static_cast<amrex::Real*>(
    amrex::The_Managed_Arena()->alloc(sizeof(amrex::Real) * array_size));
  sig = static_cast<amrex::Real*>(
    amrex::The_Managed_Arena()->alloc(sizeof(amrex::Real) * array_size));
  dip = static_cast<amrex::Real*>(
    amrex::The_Managed_Arena()->alloc(sizeof(amrex::Real) * array_size));
  pol = static_cast<amrex::Real*>(
    amrex::The_Managed_Arena()->alloc(sizeof(amrex::Real) * array_size));
  zrot = static_cast<amrex::Real*>(
    amrex::The_Managed_Arena()->alloc(sizeof(amrex::Real) * array_size));

  fitmu = static_cast<amrex::Real*>(amrex::The_Managed_Arena()->alloc(
    sizeof(amrex::Real) * array_size * fit_length));
  fitlam = static_cast<amrex::Real*>(amrex::The_Managed_Arena()->alloc(
    sizeof(amrex::Real) * array_size * fit_length));
  fitdbin = static_cast<amrex::Real*>(amrex::The_Managed_Arena()->alloc(
    sizeof(amrex::Real) * array_size * array_size * fit_length));

  nlin = static_cast<int*>(
    amrex::The_Managed_Arena()->alloc(sizeof(int) * array_size));

  egtransetWT(wt);
  egtransetEPS(eps);
  egtransetSIG(sig);
  egtransetDIP(dip);
  egtransetPOL(pol);
  egtransetZROT(zrot);
  egtransetNLIN(nlin);
  egtransetCOFETA(fitmu);
  egtransetCOFLAM(fitlam);
  egtransetCOFD(fitdbin);

  for (int i = 0; i < array_size; ++i) {
    iwt[i] = 1. / wt[i];
  }
}

void
finalize()
{
  amrex::The_Managed_Arena()->free(wt);
  amrex::The_Managed_Arena()->free(iwt);
  amrex::The_Managed_Arena()->free(eps);
  amrex::The_Managed_Arena()->free(sig);
  amrex::The_Managed_Arena()->free(dip);
  amrex::The_Managed_Arena()->free(pol);
  amrex::The_Managed_Arena()->free(zrot);
  amrex::The_Managed_Arena()->free(fitmu);
  amrex::The_Managed_Arena()->free(fitlam);
  amrex::The_Managed_Arena()->free(fitdbin);
  amrex::The_Managed_Arena()->free(nlin);
}

} // namespace transport_params
