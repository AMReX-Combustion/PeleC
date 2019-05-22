#ifndef PeleC_init_eb_F_H
#define PeleC_init_eb_F_H

#include "AMReX_ArrayLim.H"

extern "C"
{
  void pc_fill_sv_ebg(const int*  lo, const int*  hi,
                      EBBndryGeom* ebg, const int* Nebg,
                      const amrex_real* vfrac, const int* vflo, const int* vfhi,
                      const amrex_real* bcent, const int* blo, const int* bhi,
                      const amrex_real* apx, const int* axlo, const int* axhi
#if BL_SPACEDIM > 1
                      ,const amrex_real* apy, const int* aylo, const int* ayhi
#if BL_SPACEDIM > 2
                      ,const amrex_real* apz, const int* azlo, const int* azhi
#endif
#endif
    );


  void pc_fill_flux_interp_stencil(const int*  lo, const int*  hi,
                                     const int* slo, const int* shi,
                                     FaceSten* faceSten, const int* Nsten, const int* idir,
                                     const amrex_real* fcent, const int* fclo, const int* fchi,
                                     const amrex_real* fap,   const int* falo, const int* fahi);

    void pc_fill_bndry_grad_stencil(const int*  lo, const int*  hi,
                                    EBBndryGeom* ebg, const int* Ngeom,
                                    EBBndrySten* bndry_grad_stencil, const int* Ngrad,
                                    const Real* dx);

  void pc_fill_bndry_grad_stencil_amrex(const int*  lo, const int*  hi,
                                        EBBndryGeom* ebg, const int* Ngeom,
                                        EBBndrySten* bndry_grad_stencil, const int* Ngrad,
                                        const Real* dx);

  void pc_fill_bndry_grad_stencil_ls(const int*  lo, const int*  hi,
                                     EBBndryGeom* ebg, const int* Ngeom,
                                     EBBndrySten* bndry_grad_stencil, const int* Ngrad,
                                     const Real* dx);

    void pc_apply_eb_boundry_flux_stencil(const int*  lo, const int*  hi,
                                          EBBndrySten* sten, const int* Nsten,
                                          const amrex_real* s,  const int* slo, const int* shi,
                                          const amrex_real* D, const int* Dlo, const int* Dhi,
                                          const Real* bcval, const int* Nvals,
                                          Real* bcflux, const int* Nflux, const int* nc);

    void pc_apply_eb_boundry_visc_flux_stencil(const int*  lo, const int*  hi,
                                               EBBndrySten* sten, const int* Nsten,
                                               EBBndryGeom* ebg, const int* Ngeom,
                                               const amrex_real* s,  const int*  slo, const int* shi,
                                               const amrex_real* mu, const int* mulo, const int* muhi,
                                               const amrex_real* xi, const int* xilo, const int* xihi,
                                               const Real* bcval, const int* Nvals,
                                               Real* bcflux, const int* Nflux, const int* nc);

    void pc_apply_face_stencil(const int*  lo, const int*  hi,
                               const int* slo, const int* shi,
                               const FaceSten* redistSten, const int* Nsten, const int* idir,
                               const amrex_real* vin, const int* vinlo, const int* vinhi,
                               const amrex_real* vout, const int* voutlo, const int* vouthi,
                               const int* nComp, const int* in_place);

    void pc_fix_div_and_redistribute(const int*  lo, const int*  hi,
                                     const EBBndryGeom* sv_ebg, const int* Ncut,
                                     const void* flag, const int* fglo, const int* fghi,
                                     const amrex_real* f0in, const int* f0lo, const int* f0hi,
#if BL_SPACEDIM > 1
                                     const amrex_real* f1in, const int* f1lo, const int* f1hi,
#if BL_SPACEDIM > 2
                                     const amrex_real* f2in, const int* f2lo, const int* f2hi,
#endif
#endif
                                     const amrex_real* ebflux, const int* Nebflux,
                                     const amrex_real* DC,   const int* DClo,   const int* DChi,
                                     const amrex_real* vf,   const int* vflo,   const int* vfhi,
                                     const amrex_real* W,    const int* Wlo,    const int* Whi,
                                     const Real* VOL, const int* nComp,
                                     const int* as_crse,
                                     amrex_real* rr_drho_crse, const int* rdclo, const int* rdchi,
                                     const int* rr_flag_crse, const int* rfclo, const int* rfchi,
                                     const int* as_fine,
                                     amrex_real* dm_as_fine, const int* dflo, const int* dfhi,
                                     const int* levmsk, const int* lmlo, const int* lmhi,const amrex_real* dt);

}

#endif
