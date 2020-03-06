#ifndef _EOS_H_
#define _EOS_H_

//#include <cmath>

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>
#include <AMReX_Arena.H>
#include <AMReX_ParmParse.H>

#include "mechanism.h"
#include "chemistry_file.H"

#include "IndexDefines.H"
#include "Constants.H"

namespace EOS {

extern AMREX_GPU_DEVICE_MANAGED amrex::Real gamma;

void init();

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
get_wb(amrex::Real mass[], amrex::Real& wbar)
{
  const amrex::Real molec_wt[NUM_SPECIES] = {28.97};
  amrex::Real summ = 0.0;
  for (int i = 0; i < NUM_SPECIES; ++i)
    summ += mass[i] * (1.0 / molec_wt[i]);
  wbar = 1.0 / summ;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
get_cv(amrex::Real mass[], amrex::Real T, amrex::Real& cv)
{
  amrex::Real wbar;
  get_wb(mass, wbar);
  cv = RU / (wbar * (gamma - 1.0));
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
get_cp(amrex::Real T, amrex::Real mass[], amrex::Real& cp)
{
  amrex::Real cv;
  get_cv(mass, T, cv);
  cp = gamma * cv;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
get_p(
  amrex::Real rho,
  amrex::Real massfrac[],
  amrex::Real& e,
  amrex::Real& T,
  amrex::Real& p)
{
  p = (gamma - 1.0) * rho * e;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
get_cs_no_T(amrex::Real rho, amrex::Real p, amrex::Real mass[], amrex::Real& cs)
{
  cs = std::sqrt(gamma * p / rho);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
get_cs(amrex::Real rho, amrex::Real T, amrex::Real mass[], amrex::Real& cs)
{
  amrex::Real wbar, p;
  get_wb(mass, wbar);
  amrex::Real cv = RU / (wbar * (gamma - 1.0));
  amrex::Real e = cv * T;
  get_p(rho, mass, e, T, p);
  cs = std::sqrt(gamma * p / rho);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
cmpT(amrex::Real e, amrex::Real mass[], amrex::Real& T)
{
  amrex::Real poverrho, wbar;
  poverrho = (gamma - 1.0) * e;
  get_wb(mass, wbar);
  T = poverrho * wbar / RU;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
ctop(
  amrex::Real massfrac[],
  amrex::Real rho,
  amrex::Real& e,
  amrex::Real& T,
  amrex::Real& p,
  amrex::Real& dpdr_e,
  amrex::Real& dpde,
  amrex::Real& gam,
  amrex::Real& cs,
  amrex::Real& wbar)
{
  get_wb(massfrac, wbar);
  cmpT(e, massfrac, T);
  get_p(rho, massfrac, e, T, p);
  gam = gamma;
  dpdr_e = (gam - 1.0) * e;
  dpde = (gam - 1.0) * rho;
  get_cs(rho, T, massfrac, cs);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
TrhoY2p(amrex::Real rho, amrex::Real T, amrex::Real mass[], amrex::Real& p)
{
  amrex::Real cv, e;
  get_cv(mass, T, cv);
  e = cv * T;
  p = (gamma - 1.0) * rho * e;
}

/*
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void mpr2wdot(amrex::Real dens, amrex::Real T, amrex::Real mass[],
amrex::Real wdot[])
{
    CKWYR(&dens, &T, mass, wdot);
    amrex::Real mw[NUM_SPECIES];
    get_mw(mw);
    for(int n = 0; n < NUM_SPECIES; n++) wdot[n] *= mw[n];
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void get_ei(amrex::Real T, amrex::Real enrgi[])
{
    CKUMS(&T, enrgi);
}
*/

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
ytx(amrex::Real mass[], amrex::Real mole[])
{
  for (int n = 0; n < NUM_SPECIES; n++)
    mole[n] = mass[n];
}
/*
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void xty(amrex::Real mole[], amrex::Real mass[])
{
    CKXTY(mole, mass);
}
*/
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
tp(
  amrex::Real pres,
  amrex::Real massfrac[],
  amrex::Real T,
  amrex::Real& rho,
  amrex::Real& energy)
{
  amrex::Real wbar;
  get_wb(massfrac, wbar);
  rho = pres * wbar / (RU * T);
  energy = pres / (rho * (gamma - 1.0));
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
rhopY2e(
  amrex::Real rho, amrex::Real massfrac[], amrex::Real p, amrex::Real& energy)
{
  amrex::Real wbar;
  get_wb(massfrac, wbar);
  energy = p / (rho * (gamma - 1.0));
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
get_hi(amrex::Real mass[], amrex::Real T, amrex::Real enthalpy[])
{
  amrex::Real wbar;
  get_wb(mass, wbar);
  const amrex::Real cv = RU / (wbar * (gamma - 1.0));
  for (int n = 0; n < NUM_SPECIES; n++)
    enthalpy[n] = cv * T * gamma;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
get_s(amrex::Real& s)
{
  s = 1.0;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
get_gamma(amrex::Real e, amrex::Real& gam)
{
  gam = gamma;
}

// Hydro -> Advected -> Species -> Aux
// If num_adv == 0 -> QFA = QFS and UFA = UFS, see IndexDefines.H
// For explicit definitions.
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
int
upass_map(const int i)
{
  return i + UFA;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
int
qpass_map(const int i)
{
  return i + QFA;
}

}; // namespace EOS

void pc_eos_init();

#endif
