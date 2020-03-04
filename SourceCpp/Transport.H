#ifndef _TRANSPORT_H_
#define _TRANSPORT_H_

#include <AMReX.H>
#include <AMReX_Gpu.H>
#include <AMReX_FArrayBox.H>

#ifdef SIMPLE_TRANSPORT
#include "SimpleTransportParams.H"
#elif CONSTANT_TRANSPORT
#include "ConstantTransportParams.H"
#endif

void pc_transport_init();
void pc_transport_close();

AMREX_GPU_DEVICE
void pc_get_transport_coeffs(
  amrex::Box const& bx,
  amrex::Array4<const amrex::Real> const& q,
  amrex::Array4<amrex::Real> const& D);

#endif
