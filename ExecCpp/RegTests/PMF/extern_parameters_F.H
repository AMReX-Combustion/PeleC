
#ifndef _external_parameters_F_H_
#define _external_parameters_F_H_
#include <AMReX.H>
#include <AMReX_BLFort.H>

#ifdef __cplusplus
#include <AMReX.H>
extern "C"
{
#endif
  void get_f90_small_massfrac(amrex::Real* small_massfrac);

  void get_f90_react_T_min(amrex::Real* react_T_min);

  void get_f90_react_T_max(amrex::Real* react_T_max);

  void get_f90_react_rho_min(amrex::Real* react_rho_min);

  void get_f90_react_rho_max(amrex::Real* react_rho_max);

  void get_f90_new_Jacobian_each_cell(int* new_Jacobian_each_cell);

  void get_f90_numaux(int* numaux);

  void get_f90_auxnamesin(char* auxnamesin);

  void get_f90_auxnamesin_len(int& slen);

  void get_f90_use_bulk_viscosity(int* use_bulk_viscosity);


#ifdef __cplusplus
}
#endif

#endif
