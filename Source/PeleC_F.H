#ifndef _PeleC_F_H_
#define _PeleC_F_H_
#include <AMReX_BLFort.H>

#ifdef __cplusplus
#include <EBStencilTypes.H>

extern "C"
{
#endif

  void pc_network_init();

  void pc_network_close();

  void pc_reactor_init();

  void pc_reactor_close();

  void pc_transport_init();

  void pc_transport_close();

  void pc_extern_init(const int* name, const int* namlen);

  void pc_prob_close();

  void get_num_spec(int* nspec);

  void get_num_aux(int* naux);

  void get_spec_names(int* spec_names, int* ispec, int* len);

  void get_aux_names(int* aux_names, int* iaux, int* len);

  void pc_set_special_tagging_flag(const amrex::Real& max_den, const int* flag);

  void set_amr_info(const int& level, const int& iteration, const int& ncycle,
		    const amrex::Real& time, const amrex::Real& dt);

  void get_method_params(int* NUM_GROW,int* QTHERM,int* QVAR,int* cQRHO,int* cQU,
                         int* cQV,int* cQW,int* cQGAME,int* cQPRES,int* cQREINT,
                         int* cQTEMP,int* cQFA,int* cQFS,int* cQFX,int* NQAUX,
                         int* cQGAMC,int* cQC,int* cQCSML,int* cQDPDR,int* cQDPDE,int* cQRSPEC);

  void clear_method_params();

  void set_method_params
    (const int& dm,
     const int& Density, const int& Xmom,
     const int& Eden,    const int& Eint,
     const int& Temp     , const int& FirstAdv,
     const int& FirstSpec, const int& FirstAux,
     const int& NumAdv,
     const amrex::Real& diffuse_cutoff_density,
     const int& pstate_loc, const int& pstate_vel,
     const int& pstate_T, const int& pstate_dia,
     const int& pstate_rho, const int& pstate_spc,
     const int& pfld_vel, const int& pfld_rho,
     const int& pfld_T, const int& pfld_p,
     const int& pfld_spc
     );

  void set_pelec_method_params();

  void init_godunov_indices();
  void set_problem_params
    (const int& dm,
     const int* physbc_lo, const int* physbc_hi,
     const int& Interior_value, const int& userBC_value, const int& Inflow_value,
     const int& Outflow_value, const int& Symmetry_value,
     const int& SlipWall_value, const int& NoSlipWall_value,
     const int& coord_type,
     const amrex::Real* problo, const amrex::Real* probhi,
     const amrex::Real* center);

  void set_grid_info
    (const int& max_level, const amrex::Real* dx_level,
     const int* domlo_level, const int* domhi_level);

  void clear_grid_info();

  void get_tagging_params(const int* name, const int* namlen);

  void ctoprim(const int* lo, const int* hi,
	       const amrex::Real* u, const int* u_lo, const int* u_hi,
	       const amrex::Real* q, const int* q_lo, const int* q_hi,
               const amrex::Real* qaux, const int* qa_lo, const int* qa_hi);

  void srctoprim(const int* lo,   const int* hi,
		 const amrex::Real* q,   const int* q_lo,   const int* q_hi,
		 const amrex::Real* qaux, const int* qa_lo, const int* qa_hi,
		 const amrex::Real* src, const int* src_lo, const int* src_hi,
		 const amrex::Real* srQ, const int* srQ_lo, const int* srQ_hi);

#if (BL_SPACEDIM > 1)
  void pc_compute_tangential_vel_derivs(const int* lo, const int* hi,
                                        const int* dlo, const int* dhi,
                                        const amrex::Real* Q, const int* Qlo, const int* Qhi,
                                        amrex::Real* td, const int* tdlo, const int* tdhi,
                                        const amrex::Real* dx, const int* idir);
#endif

#ifdef PELE_USE_EB
#if (BL_SPACEDIM > 1)
  void pc_compute_tangential_vel_derivs_eb(const int* lo,  const int* hi,
                                           const int* dlo,  const int* dhi,
                                           const void* ebg, const int* Ngeom,
                                           const amrex::Real* Q, const int* Qlo, const int* Qhi,
                                           amrex::Real* td, const int* tdlo, const int* tdhi,
                                           const void* flag, const int* fglo, const int* fghi,
                                           const amrex::Real* dx, const int* idir);
#endif

  void pc_set_body_state(const int* lo, const int* hi,
                         amrex::Real* S, const int* Slo, const int* Shi,
                         const int* mask, const int* mlo, const int* mhi,
                         const amrex::Real* b, const int* nc, const int* bval);

#endif

  void pc_diffterm(const int* lo,  const int* hi,
		   const int* dmnlo,  const int* dmnhi,
		   const BL_FORT_FAB_ARG_3D(Q),
		   const BL_FORT_FAB_ARG_3D(Dx),
		   const BL_FORT_FAB_ARG_3D(mux),
		   const BL_FORT_FAB_ARG_3D(xix),
		   const BL_FORT_FAB_ARG_3D(lamx),
#if (BL_SPACEDIM > 1)
		   const BL_FORT_FAB_ARG_3D(tdx),
#endif
		   BL_FORT_FAB_ARG_3D(Ax),
		   BL_FORT_FAB_ARG_3D(fx),
#if (BL_SPACEDIM > 1)
		   const BL_FORT_FAB_ARG_3D(Dy),
		   const BL_FORT_FAB_ARG_3D(muy),
		   const BL_FORT_FAB_ARG_3D(xiy),
		   const BL_FORT_FAB_ARG_3D(lamy),
		   const BL_FORT_FAB_ARG_3D(tdy),
		   BL_FORT_FAB_ARG_3D(Ay),
		   BL_FORT_FAB_ARG_3D(fy),
#if (BL_SPACEDIM > 2)
		   const BL_FORT_FAB_ARG_3D(Dz),
		   const BL_FORT_FAB_ARG_3D(muz),
		   const BL_FORT_FAB_ARG_3D(xiz),
		   const BL_FORT_FAB_ARG_3D(lamz),
		   const BL_FORT_FAB_ARG_3D(tdz),
		   BL_FORT_FAB_ARG_3D(Az),
		   BL_FORT_FAB_ARG_3D(fz),
#endif
#endif
		   BL_FORT_FAB_ARG_3D(Vol),
		   BL_FORT_FAB_ARG_3D(D),
		   const amrex::Real* dx);

  void pc_check_initial_species
    (const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(state));

  void enforce_minimum_density
    (const amrex::Real* S_old, const int* s_old_lo, const int* s_old_hi,
     amrex::Real* S_new, const int* s_new_lo, const int* s_new_hi,
     const amrex::Real* vol, const int* vol_lo, const int* vol_hi,
     const int* lo, const int* hi,
     amrex::Real* mass_added, amrex::Real* e_added,
     amrex::Real* E_added, amrex::Real* frac_change,
     const int* verbose);

  void compute_temp
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state));

  void pc_enforce_consistent_e
    (const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(state));

  void reset_internal_e
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(S_new),
     const int& verbose);

  void pc_hypfill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo, 
     const amrex::Real* time, const int* bc);

  void pc_xmomfill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo,
     const amrex::Real* time, const int* bc);

  void pc_ymomfill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo,
     const amrex::Real* time, const int* bc);

  void pc_zmomfill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo,
     const amrex::Real* time, const int* bc);

  void pc_zerofill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo,
     const amrex::Real* time, const int* bc);

  void pc_default_fill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo,
     const amrex::Real* time, const int* bc);

  void pc_reactfill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo, 
     const amrex::Real* time, const int* bc);

  inline void pc_nullfill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo, 
     const amrex::Real* time, const int* bc)
#ifdef __cplusplus
  { amrex::Abort("Not supposed to call pc_nullfill"); }
#else
    ;
#endif

  void pc_estdt
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state),
     const amrex::Real* dx, amrex::Real* dt);

  void pc_estdt_vel_diffusion
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state),
     const amrex::Real* dx, amrex::Real* dt);
  void pc_estdt_temp_diffusion
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state),
     const amrex::Real* dx, amrex::Real* dt);
  void pc_estdt_enth_diffusion
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state),
     const amrex::Real* dx, amrex::Real* dt);

  void pc_check_timestep
    (const BL_FORT_FAB_ARG_3D(state_old),
     const BL_FORT_FAB_ARG_3D(state_new),
     const int* lo, const int* hi,
     const amrex::Real* dx, const amrex::Real* dt_old, amrex::Real* dt);

  void impose_NSCBC
    (const int    lo[], const int    hi[],
     const int domlo[], const int domhi[],
     const BL_FORT_FAB_ARG(state),
     BL_FORT_FAB_ARG(q),
     BL_FORT_FAB_ARG(qaux),
     D_DECL(BL_FORT_IFAB_ARG(x_bcMask),
            BL_FORT_IFAB_ARG(y_bcMask),
            BL_FORT_IFAB_ARG(z_bcMask)),
     const int* flag_nscbc_isAnyPerio, const int flag_nscbc_perio[],
     const amrex::Real* time, const amrex::Real dx[], const amrex::Real* dt);

  void set_bc_mask
    (const int    lo[], const int    hi[],
     const int domlo[], const int domhi[],
      D_DECL(BL_FORT_IFAB_ARG(x_bcMask),
            BL_FORT_IFAB_ARG(y_bcMask),
            BL_FORT_IFAB_ARG(z_bcMask)));
  
  void pc_umdrv
    (const int* is_finest_level,
     const amrex::Real* time,
     const int    lo[], const int    hi[],
     const int domlo[], const int domhi[],
     const BL_FORT_FAB_ARG(state),
     BL_FORT_FAB_ARG(stateout),
     BL_FORT_FAB_ARG(q),
     BL_FORT_FAB_ARG(qaux),
     BL_FORT_FAB_ARG(src_q),
     BL_FORT_FAB_ARG(update),
     D_DECL(BL_FORT_IFAB_ARG(x_bcMask),
            BL_FORT_IFAB_ARG(y_bcMask),
            BL_FORT_IFAB_ARG(z_bcMask)),
     const amrex::Real dx[], const amrex::Real* dt,
     D_DECL(BL_FORT_FAB_ARG(xflux),
            BL_FORT_FAB_ARG(yflux),
            BL_FORT_FAB_ARG(zflux)),
#if (BL_SPACEDIM < 3)
     BL_FORT_FAB_ARG(pradial),
#endif
     D_DECL(BL_FORT_FAB_ARG(xarea),
            BL_FORT_FAB_ARG(yarea),
            BL_FORT_FAB_ARG(zarea)),
#if (BL_SPACEDIM < 3)
     const BL_FORT_FAB_ARG(dloga),
#endif
     const BL_FORT_FAB_ARG(volume),
     amrex::Real* cflLoc,
     const int&  verbose,
     amrex::Real& mass_added_flux,
     amrex::Real& xmom_added_flux,
     amrex::Real& ymom_added_flux,
     amrex::Real& zmom_added_flux,
     amrex::Real& E_added_flux,
     amrex::Real& mass_lost,
     amrex::Real& xmom_lost,
     amrex::Real& ymom_lost,
     amrex::Real& zmom_lost,
     amrex::Real& eden_lost,
     amrex::Real& xang_lost,
     amrex::Real& yang_lost,
     amrex::Real& zang_lost);

    void pc_hyp_mol_flux
    (
        const int* lo, const int* hi,
        const int* domlo, const int* domhi,
        const BL_FORT_FAB_ARG_3D(q),
        const BL_FORT_FAB_ARG_3D(qaux),
        BL_FORT_FAB_ARG_3D(area_ec_x),
        BL_FORT_FAB_ARG_3D(flux_ec_x),
#if (BL_SPACEDIM > 1)
        BL_FORT_FAB_ARG_3D(area_ec_y),
        BL_FORT_FAB_ARG_3D(flux_ec_y),
#if (BL_SPACEDIM > 2)
        BL_FORT_FAB_ARG_3D(area_ec_z),
        BL_FORT_FAB_ARG_3D(flux_ec_z),
#endif
#endif
        const BL_FORT_FAB_ARG_3D(flatn),
        const BL_FORT_FAB_ARG_3D(volume),
        const BL_FORT_FAB_ARG_3D(D),
#ifdef PELEC_USE_EB
        const BL_FORT_FAB_ARG_3D(vfrac),
        const void* flag, const int* fglo, const int* fghi,
        const void* ebg, const int* nebg,
        const void* ebflux, const int* nebf,
#endif
        const amrex::Real* h);


    void pc_hyp_mol_flux_vec
    (
        const int* lo, const int* hi,
        const int* domlo, const int* domhi,
        const BL_FORT_FAB_ARG_3D(q),
        const BL_FORT_FAB_ARG_3D(qaux),
        BL_FORT_FAB_ARG_3D(area_ec_x),
        BL_FORT_FAB_ARG_3D(flux_ec_x),
#if (BL_SPACEDIM > 1)
        BL_FORT_FAB_ARG_3D(area_ec_y),
        BL_FORT_FAB_ARG_3D(flux_ec_y),
#if (BL_SPACEDIM > 2)
        BL_FORT_FAB_ARG_3D(area_ec_z),
        BL_FORT_FAB_ARG_3D(flux_ec_z),
#endif
#endif
        const BL_FORT_FAB_ARG_3D(flatn),
        const BL_FORT_FAB_ARG_3D(volume),
        const BL_FORT_FAB_ARG_3D(D),
#ifdef PELEC_USE_EB
        const void* flag, const int* fglo, const int* fghi,
        const void* ebg, const int* nebg,
        const void* ebflux, const int* nebf,
#endif
        const amrex::Real* h);

  void flush_output();

  void pc_summass
    (const int* lo, const int* hi, const BL_FORT_FAB_ARG_3D(rho),
     const amrex::Real* dx, const BL_FORT_FAB_ARG_3D(vol), amrex::Real* mass);

  void pc_sumsquared
    (const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(rho),
     const amrex::Real* dx, BL_FORT_FAB_ARG_3D(vol), amrex::Real* mass);

#ifdef REACTIONS
  void pc_react_state
    (const int* lo, const int* hi,
     const amrex::Real*  uold, const int* uo_lo, const int* uo_hi,
     amrex::Real*        unew, const int* un_lo, const int* un_hi,
     const amrex::Real*  asrc, const int* as_lo, const int* as_hi,
     const int*   mask, const int*  m_lo, const int*  m_hi,
     amrex::Real*        cost, const int*  c_lo, const int*  c_hi,
     amrex::Real*       rYdot, const int* rY_lo, const int* rY_hi,
     const amrex::Real& time,  const amrex::Real& dt_react, const int& do_react);
  
  void pc_react_state_expl
    (const int* lo, const int* hi,
     const amrex::Real*  uold, const int* uo_lo, const int* uo_hi,
     amrex::Real*        unew, const int* un_lo, const int* un_hi,
     const amrex::Real*  asrc, const int* as_lo, const int* as_hi,
     const int*   mask, const int*  m_lo, const int*  m_hi,
     amrex::Real*        cost, const int*  c_lo, const int*  c_hi,
     amrex::Real*       rYdot, const int* rY_lo, const int* rY_hi,
     const amrex::Real& time,  const amrex::Real& dt_react, const int& do_react,
     const int& nsubsteps_min,const int &nsubsteps_max,const int &nsubsteps_guess);
#endif

  void pc_diffextrap
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(TempDiffTerm),
     const int* nc);

  void pc_tempdiffextrap
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(TempDiffTerm));

  void pc_move_transport_coeffs_to_ec
    (const int* lo, const int* hi,
     const int* dlo, const int* dhi,
     const BL_FORT_FAB_ARG_3D(cdata),
     const BL_FORT_FAB_ARG_3D(edata),
     const int* dir, const int* nc,
     const int* do_harmonic);

#ifdef USE_MASA
  void pc_mms_src(const int* lo, const int* hi,
		  const BL_FORT_FAB_ARG_3D(S),
		  BL_FORT_FAB_ARG_3D(mms_src),
		  const amrex::Real* problo, const amrex::Real* dx,
		  const amrex::Real* time);
#endif

  void pc_initdata(const int& level, const amrex::Real& time,
                   const int* lo, const int* hi,
                   const int& num_state,
                   BL_FORT_FAB_ARG_3D(state),
                   const amrex::Real* dx,
                   const amrex::Real* xlo,
                   const amrex::Real* xhi);

  void pc_ext_src(const int* lo, const int* hi,
                  const BL_FORT_FAB_ARG_3D(old_state),
                  const BL_FORT_FAB_ARG_3D(new_state),
                  BL_FORT_FAB_ARG_3D(ext_src),
                  const amrex::Real* prob_lo, const amrex::Real* dx,
                  const amrex::Real* time, const amrex::Real* dt);

  void pc_forcing_src(const int* lo, const int* hi,
                      const BL_FORT_FAB_ARG_3D(old_state),
                      const BL_FORT_FAB_ARG_3D(new_state),
                      BL_FORT_FAB_ARG_3D(forcing_src),
                      const amrex::Real* prob_lo, const amrex::Real* dx,
                      const amrex::Real* xlo, const amrex::Real* xhi,
                      const amrex::Real* time, const amrex::Real* dt);

  void pc_set_synthetic_data(const int* lo, const int*hi,
                             BL_FORT_FAB_ARG_3D(M),
                             const amrex::Real* prob_lo, const amrex::Real* dx);
#ifdef __cplusplus
}
#endif

/* problem-specific includes */
#include <Problem_F.H>

#endif
