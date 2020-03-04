
#ifndef _external_parameters_H_
#define _external_parameters_H_
#include <AMReX_BLFort.H>

  void init_extern_parameters();

  extern AMREX_GPU_MANAGED amrex::Real small_massfrac;

  extern AMREX_GPU_MANAGED amrex::Real react_T_min;

  extern AMREX_GPU_MANAGED amrex::Real react_T_max;

  extern AMREX_GPU_MANAGED amrex::Real react_rho_min;

  extern AMREX_GPU_MANAGED amrex::Real react_rho_max;

  extern AMREX_GPU_MANAGED int new_Jacobian_each_cell;

  extern AMREX_GPU_MANAGED int numaux;

  extern std::string auxnamesin;

  extern AMREX_GPU_MANAGED int use_bulk_viscosity;


#endif
