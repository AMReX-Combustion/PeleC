#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <AMReX_Gpu.H>

extern "C"
{
AMREX_GPU_HOST_DEVICE void get_imw(double imw_new[]);
AMREX_GPU_HOST_DEVICE void get_mw(double mw_new[]);
void egtransetEPS(double *  EPS);
void egtransetSIG(double* SIG);
void atomicWeight(double *  awt);
void molecularWeight(double *  wt);
AMREX_GPU_HOST_DEVICE void gibbs(double *  species, double *  tc);
AMREX_GPU_HOST_DEVICE void helmholtz(double *  species, double *  tc);
AMREX_GPU_HOST_DEVICE void speciesInternalEnergy(double *  species, double *  tc);
AMREX_GPU_HOST_DEVICE void speciesEnthalpy(double *  species, double *  tc);
AMREX_GPU_HOST_DEVICE void speciesEntropy(double *  species, double *  tc);
AMREX_GPU_HOST_DEVICE void cp_R(double *  species, double *  tc);
AMREX_GPU_HOST_DEVICE void cv_R(double *  species, double *  tc);
void equilibriumConstants(double *  kc, double *  g_RT, double T);
AMREX_GPU_HOST_DEVICE void productionRate(double *  wdot, double *  sc, double T);
AMREX_GPU_HOST_DEVICE void comp_qfqr(double *  q_f, double *  q_r, double *  sc, double *  tc, double invT);
#ifndef AMREX_USE_CUDA
void comp_k_f(double *  tc, double invT, double *  k_f);
void comp_Kc(double *  tc, double invT, double *  Kc);
#endif
AMREX_GPU_HOST_DEVICE void progressRate(double *  qdot, double *  speciesConc, double T);
AMREX_GPU_HOST_DEVICE void progressRateFR(double *  q_f, double *  q_r, double *  speciesConc, double T);
AMREX_GPU_HOST_DEVICE void CKINIT();
AMREX_GPU_HOST_DEVICE void CKFINALIZE();
#ifndef AMREX_USE_CUDA
void GET_REACTION_MAP(int *  rmap);
void SetAllDefaults();
#endif
void CKINDX(int * mm, int * kk, int * ii, int * nfit );
void CKXNUM(char * line, int * nexp, int * lout, int * nval, double *  rval, int * kerr, int lenline);
void CKSNUM(char * line, int * nexp, int * lout, char * kray, int * nn, int * knum, int * nval, double *  rval, int * kerr, int lenline, int lenkray);
void CKSYME(int * kname, int * lenkname);
void CKSYMS(int * kname, int * lenkname);
void CKRP(double *  ru, double *  ruc, double *  pa);
void CKPX(double *  rho, double *  T, double *  x, double *  P);
AMREX_GPU_HOST_DEVICE void CKPY(double *  rho, double *  T, double *  y, double *  P);
void CKPC(double *  rho, double *  T, double *  c, double *  P);
void CKRHOX(double *  P, double *  T, double *  x, double *  rho);
AMREX_GPU_HOST_DEVICE void CKRHOY(double *  P, double *  T, double *  y, double *  rho);
void CKRHOC(double *  P, double *  T, double *  c, double *  rho);
void CKWT(double *  wt);
void CKAWT(double *  awt);
void CKMMWY(double *  y, double *  wtm);
void CKMMWX(double *  x, double *  wtm);
void CKMMWC(double *  c, double *  wtm);
AMREX_GPU_HOST_DEVICE void CKYTX(double *  y, double *  x);
void CKYTCP(double *  P, double *  T, double *  y, double *  c);
AMREX_GPU_HOST_DEVICE void CKYTCR(double *  rho, double *  T, double *  y, double *  c);
AMREX_GPU_HOST_DEVICE void CKXTY(double *  x, double *  y);
void CKXTCP(double *  P, double *  T, double *  x, double *  c);
void CKXTCR(double *  rho, double *  T, double *  x, double *  c);
void CKCTX(double *  c, double *  x);
void CKCTY(double *  c, double *  y);
void CKCPOR(double *  T, double *  cpor);
void CKHORT(double *  T, double *  hort);
void CKSOR(double *  T, double *  sor);
void CKCVML(double *  T, double *  cvml);
void CKCPML(double *  T, double *  cvml);
void CKUML(double *  T, double *  uml);
void CKHML(double *  T, double *  uml);
void CKGML(double *  T, double *  gml);
void CKAML(double *  T, double *  aml);
void CKSML(double *  T, double *  sml);
AMREX_GPU_HOST_DEVICE void CKCVMS(double *  T, double *  cvms);
AMREX_GPU_HOST_DEVICE void CKCPMS(double *  T, double *  cvms);
AMREX_GPU_HOST_DEVICE void CKUMS(double *  T, double *  ums);
AMREX_GPU_HOST_DEVICE void CKHMS(double *  T, double *  ums);
void CKGMS(double *  T, double *  gms);
void CKAMS(double *  T, double *  ams);
void CKSMS(double *  T, double *  sms);
void CKCPBL(double *  T, double *  x, double *  cpbl);
AMREX_GPU_HOST_DEVICE void CKCPBS(double *  T, double *  y, double *  cpbs);
void CKCVBL(double *  T, double *  x, double *  cpbl);
AMREX_GPU_HOST_DEVICE void CKCVBS(double *  T, double *  y, double *  cpbs);
void CKHBML(double *  T, double *  x, double *  hbml);
AMREX_GPU_HOST_DEVICE void CKHBMS(double *  T, double *  y, double *  hbms);
void CKUBML(double *  T, double *  x, double *  ubml);
AMREX_GPU_HOST_DEVICE void CKUBMS(double *  T, double *  y, double *  ubms);
void CKSBML(double *  P, double *  T, double *  x, double *  sbml);
void CKSBMS(double *  P, double *  T, double *  y, double *  sbms);
void CKGBML(double *  P, double *  T, double *  x, double *  gbml);
void CKGBMS(double *  P, double *  T, double *  y, double *  gbms);
void CKABML(double *  P, double *  T, double *  x, double *  abml);
void CKABMS(double *  P, double *  T, double *  y, double *  abms);
AMREX_GPU_HOST_DEVICE void CKWC(double *  T, double *  C, double *  wdot);
void CKWYP(double *  P, double *  T, double *  y, double *  wdot);
void CKWXP(double *  P, double *  T, double *  x, double *  wdot);
AMREX_GPU_HOST_DEVICE void CKWYR(double *  rho, double *  T, double *  y, double *  wdot);
void CKWXR(double *  rho, double *  T, double *  x, double *  wdot);
void CKQC(double *  T, double *  C, double *  qdot);
void CKKFKR(double *  P, double *  T, double *  x, double *  q_f, double *  q_r);
void CKQYP(double *  P, double *  T, double *  y, double *  qdot);
void CKQXP(double *  P, double *  T, double *  x, double *  qdot);
void CKQYR(double *  rho, double *  T, double *  y, double *  qdot);
void CKQXR(double *  rho, double *  T, double *  x, double *  qdot);
void CKNU(int * kdim, int * nuki);
#ifndef AMREX_USE_CUDA
void CKINU(int * i, int * nspec, int * ki, int * nu);
#endif
void CKNCF(int * mdim, int * ncf);
void CKABE(double *  a, double *  b, double *  e );
void CKEQC(double *  T, double *  C , double *  eqcon );
void CKEQYP(double *  P, double *  T, double *  y, double *  eqcon);
void CKEQXP(double *  P, double *  T, double *  x, double *  eqcon);
void CKEQYR(double *  rho, double *  T, double *  y, double *  eqcon);
void CKEQXR(double *  rho, double *  T, double *  x, double *  eqcon);
AMREX_GPU_HOST_DEVICE void DWDOT(double *  J, double *  sc, double *  T, int * consP);
AMREX_GPU_HOST_DEVICE void DWDOT_PRECOND(double *  J, double *  sc, double *  Tp, int * HP);
AMREX_GPU_HOST_DEVICE void SLJ_PRECOND_CSC(double *  Jsps, int * indx, int * len, double * sc, double * Tp, int * HP, double * gamma);
AMREX_GPU_HOST_DEVICE void SPARSITY_INFO(int * nJdata, int * consP, int NCELLS);
AMREX_GPU_HOST_DEVICE void SPARSITY_INFO_PRECOND(int * nJdata, int * consP);
AMREX_GPU_HOST_DEVICE void SPARSITY_PREPROC(int * rowVals, int * colPtrs, int * consP, int NCELLS);
#ifndef AMREX_USE_CUDA
void SPARSITY_PREPROC_PRECOND(int * rowVals, int * colPtrs, int * indx, int * consP);
#else
AMREX_GPU_HOST_DEVICE void SPARSITY_PREPROC_PRECOND(int * rowPtr, int * colIndx, int * consP);
#endif
AMREX_GPU_HOST_DEVICE void aJacobian(double *  J, double *  sc, double T, int consP);
AMREX_GPU_HOST_DEVICE void aJacobian_precond(double *  J, double *  sc, double T, int HP);
AMREX_GPU_HOST_DEVICE void dcvpRdT(double *  species, double *  tc);
AMREX_GPU_HOST_DEVICE void GET_T_GIVEN_EY(double *  e, double *  y, double *  t, int *ierr);
AMREX_GPU_HOST_DEVICE void GET_T_GIVEN_HY(double *  h, double *  y, double *  t, int *ierr);
void GET_CRITPARAMS(double *  Tci, double *  ai, double *  bi, double *  acentric_i);
/*vector version */
void VCKYTX(int *  np, double *  y, double *  x);
void VCKHMS(int *  np, double *  T, double *  ums);
void VCKWYR(int *  np, double *  rho, double *  T,
            double *  y,
            double *  wdot);
#ifndef AMREX_USE_CUDA
void vproductionRate(int npt, double *  wdot, double *  c, double *  T);
void VCKPY(int *  np, double *  rho, double *  T, double *  y, double *  P);
void vcomp_k_f(int npt, double *  k_f_s, double *  tc, double *  invT);
void vcomp_gibbs(int npt, double *  g_RT, double *  tc);
void vcomp_Kc(int npt, double *  Kc_s, double *  g_RT, double *  invT);
void vcomp_wdot(int npt, double *  wdot, double *  mixture, double *  sc,
                double *  k_f_s, double *  Kc_s,
                double *  tc, double *  invT, double *  T);
#endif
/*Transport function declarations */
void egtransetLENIMC(int* LENIMC);
void egtransetLENRMC(int* LENRMC);
void egtransetNO(int* NO);
void egtransetKK(int* KK);
void egtransetNLITE(int* NLITE);
void egtransetPATM(double* PATM);
void egtransetWT(double* WT);
void egtransetEPS(double* EPS);
void egtransetSIG(double* SIG);
void egtransetDIP(double* DIP);
void egtransetPOL(double* POL);
void egtransetZROT(double* ZROT);
void egtransetNLIN(int* NLIN);
void egtransetCOFETA(double* COFETA);
void egtransetCOFLAM(double* COFLAM);
void egtransetCOFD(double* COFD);
void egtransetKTDIF(int* KTDIF);
}

#ifndef AMREX_USE_CUDA
namespace thermo
{

    extern double fwd_A[0], fwd_beta[0], fwd_Ea[0];
    extern double low_A[0], low_beta[0], low_Ea[0];
    extern double rev_A[0], rev_beta[0], rev_Ea[0];
    extern double troe_a[0],troe_Ts[0], troe_Tss[0], troe_Tsss[0];
    extern double sri_a[0], sri_b[0], sri_c[0], sri_d[0], sri_e[0];
    extern double activation_units[0], prefactor_units[0], phase_units[0];
    extern int is_PD[0], troe_len[0], sri_len[0], nTB[0], *TBid[0];
    extern double *TB[0];
    extern std::vector<std::vector<double>> kiv; 
    extern std::vector<std::vector<double>> nuv; 

    extern double fwd_A_DEF[0], fwd_beta_DEF[0], fwd_Ea_DEF[0];
    extern double low_A_DEF[0], low_beta_DEF[0], low_Ea_DEF[0];
    extern double rev_A_DEF[0], rev_beta_DEF[0], rev_Ea_DEF[0];
    extern double troe_a_DEF[0],troe_Ts_DEF[0], troe_Tss_DEF[0], troe_Tsss_DEF[0];
    extern double sri_a_DEF[0], sri_b_DEF[0], sri_c_DEF[0], sri_d_DEF[0], sri_e_DEF[0];
    extern double activation_units_DEF[0], prefactor_units_DEF[0], phase_units_DEF[0];
    extern int is_PD_DEF[0], troe_len_DEF[0], sri_len_DEF[0], nTB_DEF[0], *TBid_DEF[0];
    extern double *TB_DEF[0];
    extern std::vector<int> rxn_map;
}
#endif
