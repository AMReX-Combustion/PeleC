#ifndef TRANSCOEFF_H
#define TRANSCOEFF_H

#include "prob.H"

// This header file contains functions and declarations for diffterm.
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
pc_transcoeff(
  const bool get_xi,
  const bool get_mu,
  const bool get_lam,
  const bool get_Ddiag,
  const bool get_chi,
  const amrex::Real Tloc,
  const amrex::Real rholoc,
  amrex::Real* Yloc,
  amrex::Real* Ddiag,
  amrex::Real* chi_mix,
  amrex::Real& mu,
  amrex::Real& xi,
  amrex::Real& lam,
  pele::physics::transport::TransParm<
    pele::physics::PhysicsType::eos_type,
    pele::physics::PhysicsType::transport_type> const* tparm,
  ProbParmDevice const& prob_parm,
  const amrex::RealVect& x)
{
  auto trans = pele::physics::PhysicsType::transport();
  trans.transport(
    get_xi, get_mu, get_lam, get_Ddiag, get_chi, Tloc, rholoc, Yloc, Ddiag,
    chi_mix, mu, xi, lam, tparm);
  ProblemSpecificFunctions::problem_modify_transport_coeffs(
    get_xi, get_mu, get_lam, get_Ddiag, get_chi, Tloc, rholoc, Yloc, Ddiag,
    chi_mix, mu, xi, lam, tparm, prob_parm, x);
}

#endif
