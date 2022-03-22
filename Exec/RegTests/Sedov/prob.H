#ifndef _PROB_H_
#define _PROB_H_

// AMReX
#include "AMReX_Array4.H"        // for Array4
#include "AMReX_Config.H"        // for AMREX_SPACEDIM
#include "AMReX_Extension.H"     // for AMREX_FORCE_INLINE
#include "AMReX_Geometry.H"      // for GeometryData
#include "AMReX_GpuQualifiers.H" // for AMREX_GPU_DEVICE
#include "AMReX_REAL.H"          // for Real
#include "AMReX_SPACE.H"         // for AMREX_D_TERM

// PelePhysics
#include "GammaLaw.H"    // for GammaLaw
#include "mechanism.H"   // for NUM_SPECIES
#include "PelePhysics.H" // for PhysicsType

// PeleC
#include "Constants.H"     // for PI
#include "IndexDefines.H"  // for URHO, NVAR, UMX, UMY, UMZ, UFS, UEDEN
#include "ProblemDerive.H" // for EmptyProbDeriveStruct
#include "prob_parm.H"     // for ProbParmDevice
#include "Tagging.H"       // for EmptyProbTagStruct

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_initdata(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& state,
  amrex::GeometryData const& geomdata,
  ProbParmDevice const& prob_parm)
{
  // Set explosion pressure -- we will convert the point-explosion energy into
  // a corresponding pressure distributed throughout the perturbed volume
  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real vctr = 4.0 / 3.0 * constants::PI() *
                     (prob_parm.r_init * prob_parm.r_init * prob_parm.r_init);
  amrex::Real p_exp = (eos.gamma - 1.0) * prob_parm.exp_energy / vctr;

  // Geometry
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();
  const amrex::Real* dx = geomdata.CellSize();
  AMREX_D_TERM(const amrex::Real xmin = prob_lo[0] + (i + 0.5) * dx[0];
               , const amrex::Real ymin = prob_lo[1] + (j + 0.5) * dx[1];
               , const amrex::Real zmin = prob_lo[2] + (k + 0.5) * dx[2];)

  AMREX_D_TERM(
    amrex::Real dx_sub = dx[0] / static_cast<amrex::Real>(prob_parm.nsub);
    amrex::Real ctr_x = 0.5 * (prob_lo[0] + prob_hi[0]);
    , amrex::Real dy_sub = dx[1] / static_cast<amrex::Real>(prob_parm.nsub);
    amrex::Real ctr_y = 0.5 * (prob_lo[1] + prob_hi[1]);
    , amrex::Real dz_sub = dx[2] / static_cast<amrex::Real>(prob_parm.nsub);
    amrex::Real ctr_z = 0.5 * (prob_lo[2] + prob_hi[2]);)

  // We initialize by summing over subvolumes of each cell
  int npert = 0;
  int n_amb = 0;

#if AMREX_SPACEDIM > 2
  for (int kk = 0; kk < prob_parm.nsub; kk++) {
    const amrex::Real z = zmin + (kk + 0.5) * dz_sub - ctr_z;
#endif
#if AMREX_SPACEDIM > 1
    for (int jj = 0; jj < prob_parm.nsub; jj++) {
      const amrex::Real y = ymin + (jj + 0.5) * dy_sub - ctr_y;
#endif
      for (int ii = 0; ii < prob_parm.nsub; ii++) {
        const amrex::Real x = xmin + (ii + 0.5) * dx_sub - ctr_x;

        amrex::Real dist = AMREX_D_TERM(x * x, +y * y, +z * z);

        if (dist <= prob_parm.r_init * prob_parm.r_init) {
          npert += 1;
        } else {
          n_amb += 1;
        }
      }
#if AMREX_SPACEDIM > 1
    }
#endif
#if AMREX_SPACEDIM > 2
  }
#endif

  amrex::Real p_zone =
    (static_cast<amrex::Real>(npert) * p_exp +
     static_cast<amrex::Real>(n_amb) * prob_parm.p_ambient) /
    static_cast<amrex::Real>(prob_parm.nsub * prob_parm.nsub * prob_parm.nsub);

  amrex::Real eint = p_zone / (eos.gamma - 1.0);

  state(i, j, k, URHO) = prob_parm.dens_ambient;
  state(i, j, k, UMX) = 0.0;
  state(i, j, k, UMY) = 0.0;
  state(i, j, k, UMZ) = 0.0;

  state(i, j, k, UEDEN) = eint + 0.5 *
                                   (state(i, j, k, UMX) * state(i, j, k, UMX) +
                                    state(i, j, k, UMY) * state(i, j, k, UMY) +
                                    state(i, j, k, UMZ) * state(i, j, k, UMZ)) /
                                   state(i, j, k, URHO);

  state(i, j, k, UEINT) = eint;

  for (int n = 0; n < NUM_SPECIES; n++) {
    state(i, j, k, UFS + n) = 0;
  }

  state(i, j, k, UFS) = state(i, j, k, URHO);

  amrex::Real T;
  amrex::Real e = eint / state(i, j, k, URHO);
  amrex::Real massfrac[NUM_SPECIES] = {1.0};
  eos.EY2T(e, massfrac, T);
  state(i, j, k, UTEMP) = T;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real* /*x[AMREX_SPACEDIM]*/,
  const amrex::Real s_int[NVAR],
  amrex::Real s_ext[NVAR],
  const int /*idir*/,
  const int /*sgn*/,
  const amrex::Real /*time*/,
  amrex::GeometryData const& /*geomdata*/,
  ProbParmDevice const& /*prob_parm*/)
{
  for (int n = 0; n < NVAR; n++) {
    s_ext[n] = s_int[n];
  }
}

void pc_prob_close();

using ProblemTags = EmptyProbTagStruct;
using ProblemDerives = EmptyProbDeriveStruct;

#endif
