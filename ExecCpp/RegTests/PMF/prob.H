#ifndef _PROB_H_
#define _PROB_H_

#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>

#include "IndexDefines.H"
#include "EOS.H"
#include "Tagging.H"
#include "ProblemDerive.H"
#include "prob_parm.H"
#include "Constants.H"

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
pmf(
  amrex::Real xlo,
  amrex::Real xhi,
  amrex::GpuArray<amrex::Real, NUM_SPECIES + 4>& y_vector)
{
  amrex::Real sum = 0.0, xmid = 0.0;
  int lo_loside = 0, lo_hiside = 0;
  int hi_loside = 0, hi_hiside = 0;
  int loside = 0, hiside = 0;
  amrex::Real ylo = 0.0, yhi = 0.0, x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0,
              dydx = 0.0;

  if (ProbParm::pmf_do_average) {
    lo_loside = 0;
    lo_hiside = 0;
    hi_loside = 0;
    hi_hiside = 0;
    if (xlo < ProbParm::d_pmf_X[0]) {
      lo_loside = 0;
      lo_hiside = 0;
    }
    if (xhi < ProbParm::d_pmf_X[0]) {
      hi_loside = 0;
      hi_hiside = 0;
    }
    if (xlo > ProbParm::d_pmf_X[ProbParm::pmf_N - 1]) {
      lo_loside = ProbParm::pmf_N - 1;
      lo_hiside = ProbParm::pmf_N - 1;
    }
    if (xhi > ProbParm::d_pmf_X[ProbParm::pmf_N - 1]) {
      hi_loside = ProbParm::pmf_N - 1;
      hi_hiside = ProbParm::pmf_N - 1;
    }
    if (lo_loside == 0) {
      for (int i = 0; i < ProbParm::pmf_N - 1; i++) {
        if ((xlo > ProbParm::d_pmf_X[i]) && (xlo < ProbParm::d_pmf_X[i + 1])) {
          lo_loside = i;
          lo_hiside = i + 1;
        }
      }
    }
    if (hi_loside == 0) {
      for (int i = 0; i < ProbParm::pmf_N - 1; i++) {
        if ((xhi > ProbParm::d_pmf_X[i]) && (xhi < ProbParm::d_pmf_X[i + 1])) {
          hi_loside = i;
          hi_hiside = i + 1;
        }
      }
    }
    for (int j = 0; j < ProbParm::pmf_M; j++) {
      x1 = ProbParm::d_pmf_X[lo_loside];
      y1 = ProbParm::d_pmf_Y[ProbParm::pmf_N * j + lo_loside];
      x2 = ProbParm::d_pmf_X[lo_hiside];
      y2 = ProbParm::d_pmf_Y[ProbParm::pmf_N * j + lo_hiside];
      if (lo_loside == lo_hiside)
        dydx = 0.0;
      else
        dydx = (y2 - y1) / (x2 - x1);
      ylo = y1 + dydx * (xlo - x1);
      if (lo_loside == hi_loside) {
        yhi = y1 + dydx * (xhi - x1);
        y_vector[j] = 0.5 * (ylo + yhi);
      } else {
        sum = (x2 - xlo) * 0.5 * (ylo + y2);
        x1 = ProbParm::d_pmf_X[hi_loside];
        y1 = ProbParm::d_pmf_Y[ProbParm::pmf_N * j + hi_loside];
        x2 = ProbParm::d_pmf_X[hi_hiside];
        y2 = ProbParm::d_pmf_Y[ProbParm::pmf_N * j + hi_hiside];
        if (hi_loside == hi_hiside)
          dydx = 0.0;
        else
          dydx = (y2 - y1) / (x2 - x1);
        yhi = y1 + dydx * (xhi - x1);
        sum = sum + (xhi - x1) * 0.5 * (yhi + y1);
        for (int k = lo_hiside; k < hi_loside - 1; k++) {
          sum = sum + (ProbParm::d_pmf_X[k + 1] - ProbParm::d_pmf_X[k]) * 0.5 *
                        (ProbParm::d_pmf_Y[ProbParm::pmf_N * j + k] +
                         ProbParm::d_pmf_Y[ProbParm::pmf_N * j + k + 1]);
        }
        y_vector[j] = sum / (xhi - xlo);
      }
    }
  } else {
    xmid = 0.5 * (xlo + xhi);
    loside = -1;
    hiside = -1;
    if (xmid < ProbParm::d_pmf_X[0]) {
      loside = 0;
      hiside = 0;
    }
    if (xmid > ProbParm::d_pmf_X[ProbParm::pmf_N - 1]) {
      loside = ProbParm::pmf_N - 1;
      hiside = ProbParm::pmf_N - 1;
    }
    if (loside == -1) {
      for (int i = 0; i < ProbParm::pmf_N - 1; i++) {
        if (
          (xmid >= ProbParm::d_pmf_X[i]) &&
          (xmid <= ProbParm::d_pmf_X[i + 1])) {
          loside = i;
          hiside = i + 1;
        }
      }
    }
    for (int j = 0; j < ProbParm::pmf_M; j++) {
      x1 = ProbParm::d_pmf_X[loside];
      y1 = ProbParm::d_pmf_Y[ProbParm::pmf_N * j + loside];
      x2 = ProbParm::d_pmf_X[hiside];
      y2 = ProbParm::d_pmf_Y[ProbParm::pmf_N * j + hiside];
      if (loside == hiside)
        dydx = 0.0;
      else
        dydx = (y2 - y1) / (x2 - x1);
      y_vector[j] = y1 + dydx * (xlo - x1);
    }
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_initdata(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& state,
  amrex::GeometryData const& geomdata)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();
  const amrex::Real* dx = geomdata.CellSize();
  const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
  const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
  const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];
  amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {0.0};
  amrex::Real u[3] = {0.0};
  amrex::Real molefrac[NUM_SPECIES] = {0.0};
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  amrex::Real pert, rho, energy, T, pres;

  if (ProbParm::phi_in < 0.0) {
    pert = 0.0;
    if (ProbParm::pertmag > 0.0) {
      pert = ProbParm::pertmag *
             (1.0 * std::sin(2 * PI * 4 * x / ProbParm::L[0]) *
                std::sin(2 * PI * 5 * y / ProbParm::L[1]) +
              1.023 * std::sin(2 * PI * 2 * (x - 0.4598) / ProbParm::L[0]) *
                std::sin(2 * PI * 4 * (y - 0.53765) / ProbParm::L[1]) +
              0.945 * std::sin(2 * PI * 3 * (x - 0.712435) / ProbParm::L[0]) *
                std::sin(2 * PI * 3 * (y - 2.137) / ProbParm::L[1]) +
              1.017 * std::sin(2 * PI * 5 * (x - 0.33) / ProbParm::L[0]) *
                std::sin(2 * PI * 6 * (y - 1.8) / ProbParm::L[1]) +
              0.982 * std::sin(2 * PI * 5 * (x - 1.4234) / ProbParm::L[0]));
    }
    pmf(z + pert, z + pert, pmf_vals);
    for (int n = 0; n < NUM_SPECIES; n++)
      molefrac[n] = pmf_vals[3 + n];
    T = pmf_vals[0];
    pres = ProbParm::pamb;
    u[AMREX_SPACEDIM - 1] = pmf_vals[1];
    EOS::X2Y(molefrac, massfrac);
    EOS::PYT2RE(pres, massfrac, T, rho, energy);
    state(i, j, k, URHO) = rho;
    state(i, j, k, UMX) = rho * u[0];
    state(i, j, k, UMY) = rho * u[1];
    state(i, j, k, UMZ) = rho * u[2];
    state(i, j, k, UEINT) = rho * energy;
    state(i, j, k, UEDEN) =
      rho * (energy + 0.5 * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]));
    state(i, j, k, UTEMP) = T;
    for (int n = 0; n < NUM_SPECIES; n++)
      state(i, j, k, UFS + n) = rho * massfrac[n];
  } else {
    for (int n = 0; n < NVAR; n++)
      state(i, j, k, n) = ProbParm::d_fuel_state[n];
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real x[AMREX_SPACEDIM],
  const amrex::Real s_int[NVAR],
  amrex::Real s_ext[NVAR],
  const int idir,
  const int sgn,
  const amrex::Real time,
  amrex::GeometryData const& geomdata)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();
  amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {0.0};
  amrex::Real u[3] = {0.0};
  amrex::Real molefrac[NUM_SPECIES] = {0.0};
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  amrex::Real pert, rho, energy, T, pres;

  if (sgn == -1) {
    pmf(prob_hi[idir], prob_hi[idir], pmf_vals);
    for (int n = 0; n < NUM_SPECIES; n++)
      molefrac[n] = pmf_vals[3 + n];
    T = pmf_vals[0];
    pres = ProbParm::pamb;
    u[0] = 0.0;
    u[1] = 0.0;
    u[2] = 0.0;
    u[AMREX_SPACEDIM - 1] = pmf_vals[1];
    EOS::X2Y(molefrac, massfrac);
    EOS::PYT2RE(pres, massfrac, T, rho, energy);
    s_ext[URHO] = rho;
    s_ext[UMX] = rho * u[0];
    s_ext[UMY] = rho * u[1];
    s_ext[UMZ] = rho * u[2];
    s_ext[UEINT] = rho * energy;
    s_ext[UEDEN] =
      rho * (energy + 0.5 * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]));
    s_ext[UTEMP] = T;
    for (int n = 0; n < NUM_SPECIES; n++)
      s_ext[UFS + n] = rho * massfrac[n];

  } else {
    pmf(prob_lo[idir], prob_lo[idir], pmf_vals);
    for (int n = 0; n < NUM_SPECIES; n++)
      molefrac[n] = pmf_vals[3 + n];
    T = pmf_vals[0];
    pres = ProbParm::pamb;
    u[0] = 0.0;
    u[1] = 0.0;
    u[2] = 0.0;
    u[AMREX_SPACEDIM - 1] = pmf_vals[1];
    EOS::X2Y(molefrac, massfrac);
    EOS::PYT2RE(pres, massfrac, T, rho, energy);
    s_ext[URHO] = rho;
    s_ext[UMX] = rho * u[0];
    s_ext[UMY] = rho * u[1];
    s_ext[UMZ] = rho * u[2];
    s_ext[UEINT] = rho * energy;
    s_ext[UEDEN] =
      rho * (energy + 0.5 * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]));
    s_ext[UTEMP] = T;
    for (int n = 0; n < NUM_SPECIES; n++)
      s_ext[UFS + n] = rho * massfrac[n];
  }
}

void pc_prob_close();

using ProblemTags = EmptyProbTagStruct;
using ProblemDerives = EmptyProbDeriveStruct;

#endif
