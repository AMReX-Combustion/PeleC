#ifndef _Derive_F_H_
#define _Derive_F_H_
#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

  void pc_derspectrac
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
     const int* level, const int* grid_no, const int* idx);

  void pc_dertemp
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
     const int* level, const int* grid_no);

  void pc_derpres
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
     const int* level, const int* grid_no);

  void pc_dermolefrac
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
     const int* level, const int* grid_no);
     
  void pc_dereint1
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
     const int* level, const int* grid_no);

  void pc_dereint2
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
     const int* level, const int* grid_no);

  void pc_derlogden
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
     const int* level, const int* grid_no);

  void pc_deruplusc
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
     const int* level, const int* grid_no);

  void pc_deruminusc
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
     const int* level, const int* grid_no);

  void pc_dersoundspeed
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
     const int* level, const int* grid_no);

  void pc_dermachnumber
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
     const int* level, const int* grid_no);

  void pc_derentropy
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
     const int* level, const int* grid_no);

  void pc_derenuctimespcle
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec,
     const int* level, const int* grid_no);

  void pc_dervelx
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
     const int* level, const int* grid_no);

  void pc_dervely
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
     const int* level, const int* grid_no);

  void pc_dervelz
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
     const int* level, const int* grid_no);

  void pc_dermagvel
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
     const int* level, const int* grid_no);

  void pc_derradialvel
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
     const int* level, const int* grid_no);

  void pc_dermagmom
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
     const int* level, const int* grid_no);

  void pc_derangmomx
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
     const int* level, const int* grid_no);

  void pc_derangmomy
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
     const int* level, const int* grid_no);

  void pc_derangmomz
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
     const int* level, const int* grid_no);

  void pc_derkineng
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec,
     const int* level, const int* grid_no);

  void pc_derenstrophy
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec,
     const int* level, const int* grid_no);

  void pc_dernull
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec,
     const int* level, const int* grid_no);

  void pc_derspec
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
     const int* level, const int* grid_no);

  void pc_dermagvort
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
     const int* level, const int* grid_no);

  void pc_derdivu
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
     const int* level, const int* grid_no);

#ifdef DO_HIT_FORCE
  void pc_derforcing
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
		 const BL_FORT_FAB_ARG_3D(data), const int* ncomp,
		 const int* lo, const int* hi,
		 const int* domain_lo, const int* domain_hi,
		 const amrex::Real* delta, const amrex::Real* xlo,
		 const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
     const int* level, const int* grid_no);
  
  void pc_derforcex
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
		 const BL_FORT_FAB_ARG_3D(data), const int* ncomp,
		 const int* lo, const int* hi,
		 const int* domain_lo, const int* domain_hi,
		 const amrex::Real* delta, const amrex::Real* xlo,
		 const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
     const int* level, const int* grid_no);
    
  void pc_derforcey
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
		 const BL_FORT_FAB_ARG_3D(data), const int* ncomp,
		 const int* lo, const int* hi,
		 const int* domain_lo, const int* domain_hi,
		 const amrex::Real* delta, const amrex::Real* xlo,
		 const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
     const int* level, const int* grid_no);
    
  void pc_derforcez
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
		 const BL_FORT_FAB_ARG_3D(data), const int* ncomp,
		 const int* lo, const int* hi,
		 const int* domain_lo, const int* domain_hi,
		 const amrex::Real* delta, const amrex::Real* xlo,
		 const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
     const int* level, const int* grid_no);
#endif

#ifdef __cplusplus
}
#endif
  
/* problem-specific includes */
#include <Problem_Derive_F.H>

#endif
