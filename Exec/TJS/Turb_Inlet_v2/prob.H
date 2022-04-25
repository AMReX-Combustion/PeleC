#ifndef _PROB_H_
#define _PROB_H_

//#include <AMReX_Print.H>
//#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>
#include <AMReX_Random.H>

#include "mechanism.H"

#include "PeleC.H"
#include "IndexDefines.H"
#include "Constants.H"
#include "PelePhysics.H"
#include "Tagging.H"
#include "ProblemDerive.H"
#include "prob_parm.H"

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
  // Get random val for initial x-momentum.
  // Consider changing seed to be non-deterministic
  amrex::Real pert = 1.0 + ((amrex::Random() - 0.5) * PeleC::h_prob_parm_device->init_pert_mag);
  amrex::Real u_pert = PeleC::h_prob_parm_device->u_in * pert;

  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();
  const amrex::Real* dx = geomdata.CellSize();
  amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
  //amrex::Real u_profile = (8.0 / (prob_hi[1] * prob_hi[1])) * ProbParm::u_in *
  //  (prob_hi[1] - y) * (y - prob_lo[1]);

  // Set the state
  state(i, j, k, URHO) = PeleC::h_prob_parm_device->rho_in;
  state(i, j, k, UMX) = u_pert * PeleC::h_prob_parm_device->rho_in;
  state(i, j, k, UMY) = 0.0;
  state(i, j, k, UMZ) = 0.0;
  state(i, j, k, UEINT) = PeleC::h_prob_parm_device->rho_in * PeleC::h_prob_parm_device->e_in;
  state(i, j, k, UEDEN) = PeleC::h_prob_parm_device->rho_in * 
    (PeleC::h_prob_parm_device->e_in + 0.5 * (u_pert * u_pert));
  state(i, j, k, UTEMP) = PeleC::h_prob_parm_device->T_in;
  for (int n = 0; n < NUM_SPECIES; n++)
    state(i, j, k, UFS + n) = PeleC::h_prob_parm_device->rho_in * PeleC::h_prob_parm_device->massfrac[n];
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
  amrex::GeometryData const& geomdata,
  ProbParmDevice const& prob_parm)
{
  // Coefficient to optionally ramp up inflow velocity
  amrex::Real t_coef = amrex::min(1.0, time / 0.00005);

  if (1) {
    amrex::Real rho = 0.0, u = 0.0, v = 0.0, w = 0.0, eint = 0.0, T = 0.0;
    amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac_int = {0.0};

    if (sgn == 1) {
      // inflow
      // amrex::Print() << "starting bcnormal. x = " << x[0] << " " << x[1] << " " << x[2] << " " << '\n';
      
      const amrex::Real* prob_lo = geomdata.ProbLo();
      const amrex::Real* prob_hi = geomdata.ProbHi();

      // turbulence inflow

      // initialize variables
      // const amrex::Real tau = PeleC::h_prob_parm_device->turb_length_scale / PeleC::h_prob_parm_device->turb_velocity;
      // const amrex::Real t_hat = 2.0 * constants::PI() * time / tau;
      amrex::Real x_hat[AMREX_SPACEDIM]  = {0.0};
      amrex::Real turb_u[AMREX_SPACEDIM] = {0.0};

      // convert to normalized turbulence terms
      for(int i = 0; i < AMREX_SPACEDIM; i++) {
      	// removing factor of 2.0 * constants::PI() for Huang et al. method (eq. 20)
      	x_hat[i] = x[i] / (PeleC::h_prob_parm_device->Ls_prefactor * PeleC::h_prob_parm_device->turb_length_scale);
      }
      // assume inflow is purely x-dimensional. Correct for inflow speed.
      x_hat[0] -= (2.0 * constants::PI() * PeleC::h_prob_parm_device->u_in * time)
                   / (PeleC::h_prob_parm_device->Ls_prefactor * PeleC::h_prob_parm_device->turb_length_scale);
      
      //amrex::Print() << "bcnormal check #1 (before loop) " << '\n';
      
      for(int i = 0; i < 3; i++) {

        for(int m = 0; m < PeleC::h_prob_parm_device->turb_num_modes; m++) {

          for(int n = 0; n < PeleC::h_prob_parm_device->sampling_number; n++) {

            amrex::Real arg = PeleC::h_prob_parm_device->k[0][n][m] * x_hat[0]
                              + PeleC::h_prob_parm_device->k[1][n][m] * x_hat[1]
                              + PeleC::h_prob_parm_device->k[2][n][m] * x_hat[2]
                              + PeleC::h_prob_parm_device->omega[n][m] * time;
	  
            // replace t_hat with time for the Huang et al. method.
            // arg += PeleC::h_prob_parm_device->omega[n][m] * time;

            // Compute turbulence and convert the value from mks to cgs
            turb_u[i] += PeleC::h_prob_parm_device->p[i][n][m] * cos(arg)
	               + PeleC::h_prob_parm_device->q[i][n][m] * sin(arg);
          }
        }
      }

      //scale up the velocities
      //for (int i=0; i < 3; i++) {
      //    turb_u[i] += PeleC::h_prob_parm_device->u_in;
      //}

      //amrex::Print() << "bcnormal check #2 (after loop) " << turb_u[0] << " " << turb_u[1] << " " << turb_u[2] << '\n';
      
      rho = PeleC::h_prob_parm_device->rho_in;
      
      //ignore everything for the sake of tracking down bug
      u = turb_u[0]+PeleC::h_prob_parm_device->u_in;
      v = turb_u[1];
      w = turb_u[2];
      // amrex::Print() << "Turb Inlet Vel (u, v, w) = (" << u << ", " << v << ", " << w << ")\n";
      eint =  PeleC::h_prob_parm_device->e_in;
      T = PeleC::h_prob_parm_device->T_in;
      for (int n = 0; n < NUM_SPECIES; n++)
	       s_ext[UFS + n] = PeleC::h_prob_parm_device->rho_in * PeleC::h_prob_parm_device->massfrac[n];
      // amrex::Print() << "bcnormal check #3 (finished) " << '\n';

    } else if (sgn == -1) {
      // outflow
      // Following Blazek p 279, eq. 8.23

      // Interior state (point d)
      const amrex::Real* prob_hi = geomdata.ProbHi();
      const amrex::Real* dx = geomdata.CellSize();
      const amrex::Real xd = prob_hi[0] - 0.5 * dx[0];
      const amrex::Real rho_int = s_int[URHO];
      const amrex::Real u_int = s_int[UMX] / rho_int;
      const amrex::Real v_int = s_int[UMY] / rho_int;
      const amrex::Real w_int = s_int[UMZ] / rho_int;
      const amrex::Real T_int = s_int[UTEMP];
      amrex::Real p_int = 0.0, cs_int = 0.0;

      // get species massfracs from interior points.
      for (int n = 0; n < NUM_SPECIES; n++)
        massfrac_int[n] = s_int[UFS + n] / rho_int;

      auto eos = pele::physics::PhysicsType::eos();
      eos.RTY2P(rho_int, T_int, massfrac_int.begin(), p_int);
      eos.RTY2Cs(rho_int, T_int, massfrac_int.begin(), cs_int);

      // Boundary state (point b)
      const amrex::Real xb = prob_hi[0];
      const amrex::Real pb = PeleC::h_prob_parm_device->p_in;
      const amrex::Real rhob = s_int[URHO] + (pb - p_int) / (cs_int * cs_int);
      const amrex::Real ub = u_int + (p_int - pb) / (rho_int * cs_int);
      const amrex::Real vb = v_int;
      const amrex::Real wb = w_int;

      // Ghost state (point a). Linear extrapolation from d and b
      rho = (rhob - rho_int) / (xb - xd) * (x[0] - xd) + rho_int;
      const amrex::Real p = (pb - p_int) / (xb - xd) * (x[0] - xd) + p_int;

      eos.RYP2E(rho, massfrac_int.begin(), p, eint);
      eos.EY2T(eint, massfrac_int.begin(), T);

      u = (ub - u_int) / (xb - xd) * (x[0] - xd) + u_int;
      v = (vb - v_int) / (xb - xd) * (x[0] - xd) + v_int;
      w = (wb - w_int) / (xb - xd) * (x[0] - xd) + w_int;

      // Outflow-specific
      for (int n = 0; n < NUM_SPECIES; n++)
	s_ext[UFS + n] = s_int[UFS+n];
      
      // amrex::Print() << "outflow bcnormal check (finished) " << '\n';
    }

    s_ext[URHO] = rho;
    s_ext[UMX] = rho * u;
    s_ext[UMY] = rho * v;
    s_ext[UMZ] = rho * w;
    s_ext[UEINT] = rho * eint;
    s_ext[UEDEN] = rho * (eint + 0.5 * (u * u + v * v + w * w));
    s_ext[UTEMP] = T;

  }
}

void pc_prob_close();

using ProblemTags = EmptyProbTagStruct;
using ProblemDerives = EmptyProbDeriveStruct;

#endif
