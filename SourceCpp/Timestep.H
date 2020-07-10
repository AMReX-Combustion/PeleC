#ifndef _TIMESTEP_H_
#define _TIMESTEP_H_

#include <AMReX_FArrayBox.H>
#ifdef PELEC_USE_EB
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBCellFlag.H>
#endif

#include "EOS.H"
#include "IndexDefines.H"
#include "Constants.H"
#include "Transport.H"

// EstDt routines

namespace TimeStep {
extern AMREX_GPU_DEVICE_MANAGED amrex::Real max_dt;
} // namespace TimeStep

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
pc_trans4dt(
  const int which_trans,
  amrex::Real T,
  amrex::Real rho,
  amrex::Real massfrac[],
  amrex::Real& D)
{
  bool get_xi = false, get_mu = false, get_lam = false, get_Ddiag = false;
  amrex::Real dum1 = 0., dum2 = 0.;

  if (which_trans == 0) {
    get_mu = true;
    transport(
      get_xi, get_mu, get_lam, get_Ddiag, T, rho, massfrac, nullptr, D, dum1,
      dum2);
  } else if (which_trans == 1) {
    get_lam = true;
    transport(
      get_xi, get_mu, get_lam, get_Ddiag, T, rho, massfrac, nullptr, dum1, dum2,
      D);
  }
}

AMREX_GPU_HOST_DEVICE
amrex::Real pc_estdt_hydro(
  amrex::Box const& bx,
  const amrex::Array4<const amrex::Real>& u,
#ifdef PELEC_USE_EB
  const amrex::Array4<const amrex::EBCellFlag>& flags,
#endif
  AMREX_D_DECL(
    const amrex::Real& dx,
    const amrex::Real& dy,
    const amrex::Real& dz)) noexcept;

// Diffusion Velocity
AMREX_GPU_HOST_DEVICE
amrex::Real pc_estdt_veldif(
  amrex::Box const bx,
  const amrex::Array4<const amrex::Real>& u,
#ifdef PELEC_USE_EB
  const amrex::Array4<const amrex::EBCellFlag>& flags,
#endif
  AMREX_D_DECL(
    const amrex::Real& dx,
    const amrex::Real& dy,
    const amrex::Real& dz)) noexcept;

// Diffusion Temperature
AMREX_GPU_HOST_DEVICE
amrex::Real pc_estdt_tempdif(
  amrex::Box const bx,
  const amrex::Array4<const amrex::Real>& u,
#ifdef PELEC_USE_EB
  const amrex::Array4<const amrex::EBCellFlag>& flags,
#endif
  AMREX_D_DECL(
    const amrex::Real& dx,
    const amrex::Real& dy,
    const amrex::Real& dz)) noexcept;

// Diffusion Enthalpy
AMREX_GPU_HOST_DEVICE
amrex::Real pc_estdt_enthdif(
  amrex::Box const bx,
  const amrex::Array4<const amrex::Real>& u,
#ifdef PELEC_USE_EB
  const amrex::Array4<const amrex::EBCellFlag>& flags,
#endif
  AMREX_D_DECL(
    const amrex::Real& dx,
    const amrex::Real& dy,
    const amrex::Real& dz)) noexcept;

#endif
