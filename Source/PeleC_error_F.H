#ifndef _PeleC_error_F_H_
#define _PeleC_error_F_H_
#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

  void pc_laplac_error
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG_3D(var),
     const int* lo, const int* hi,
     const int* ncomp,
     const int* domlo, const int* domhi,
     const amrex::Real* dx, const amrex::Real* xlo, const amrex::Real* problo,
     const amrex::Real* time, const int* level);

  void pc_denerror
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG_3D(den),
     const int* lo, const int* hi,
     const int* ncomp,
     const int* domlo, const int* domhi,
     const amrex::Real* dx, const amrex::Real* xlo, const amrex::Real* problo,
     const amrex::Real* time, const int* level);

  void pc_velerror
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG_3D(vel),
     const int* lo, const int* hi,
     const int* ncomp,
     const int* domlo, const int* domhi,
     const amrex::Real* dx, const amrex::Real* xlo, const amrex::Real* problo,
     const amrex::Real* time, const int* level);

  void pc_vorterror
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG_3D(magvort),
     const int* lo, const int* hi,
     const int* ncomp,
     const int* domlo, const int* domhi,
     const amrex::Real* dx, const amrex::Real* xlo, const amrex::Real* problo,
     const amrex::Real* time, const int* level);

  void pc_temperror
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG_3D(temp),
     const int* lo, const int* hi,
     const int* ncomp,
     const int* domlo, const int* domhi,
     const amrex::Real* dx, const amrex::Real* xlo, const amrex::Real* problo,
     const amrex::Real* time, const int* level);

  void pc_presserror
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG_3D(press),
     const int* lo, const int* hi,
     const int* ncomp,
     const int* domlo, const int* domhi,
     const amrex::Real* dx, const amrex::Real* xlo, const amrex::Real* problo,
     const amrex::Real* time, const int* level);

  void pc_nucerror
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG_3D(state),
     const int* lo, const int* hi,
     const int* ncomp,
     const int* domlo, const int* domhi,
     const amrex::Real* dx, const amrex::Real* xlo, const amrex::Real* problo,
     const amrex::Real* time, const int* level);

  void pc_ftracerror
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG_3D(press),
     const int* lo, const int* hi,
     const int* ncomp,
     const int* domlo, const int* domhi,
     const amrex::Real* dx, const amrex::Real* xlo, const amrex::Real* problo,
     const amrex::Real* time, const int* level);

  void pc_enterror
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG_3D(den),
     const int* lo, const int* hi,
     const int* ncomp,
     const int* domlo, const int* domhi,
     const amrex::Real* dx, const amrex::Real* xlo, const amrex::Real* problo,
     const amrex::Real* time, const int* level);

  void pc_state_error
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG_3D(state),
     const int* lo, const int* hi,
     const int* ncomp,
     const int* domlo, const int* domhi,
     const amrex::Real* dx, const amrex::Real* xlo, const amrex::Real* problo,
     const amrex::Real* time, const int* level);

  void set_problem_tags
    (int* tag, const int* tag_lo, const int* tag_hi,
     BL_FORT_FAB_ARG_3D(state),
     const int* tagval, const int* clearval,
     const int* lo, const int* hi,
     const amrex::Real* dx, const amrex::Real* problo,
     const amrex::Real* time, const int* level);

  void pc_vfracerror
  (int* tag, const int* tag_lo, const int* tag_hi,
   const int* tagval, const int* clearval,
   BL_FORT_FAB_ARG_3D(vfrac),
   const int* lo, const int* hi,
   const int* ncomp,
   const int* domlo, const int* domhi,
   const amrex::Real* dx, const amrex::Real* xlo, const amrex::Real* problo,
   const amrex::Real* time, const int* level);


#ifdef __cplusplus
}
#endif

#endif
