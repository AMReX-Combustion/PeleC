#ifndef INDEX_DEFINES_H
#define INDEX_DEFINES_H

#include "mechanism.H"

static constexpr int URHO = 0;
static constexpr int UMX = 1;
static constexpr int UMY = 2;
static constexpr int UMZ = 3;
static constexpr int UEDEN = 4;
static constexpr int UEINT = 5;
static constexpr int UTEMP = 6;
static constexpr int UFA = 7;

static constexpr int QRHO = 0;
static constexpr int QU = 1;
static constexpr int QV = 2;
static constexpr int QW = 3;
static constexpr int QGAME = 4;
static constexpr int QPRES = 5;
static constexpr int QREINT = 6;
static constexpr int QTEMP = 7;
static constexpr int QFA = 8;

static constexpr int QGAMC = 0;
static constexpr int QC = 1;
static constexpr int QCSML = 2;
static constexpr int QDPDR = 3;
static constexpr int QDPDE = 4;
static constexpr int QRSPEC = 5;

static constexpr int GDRHO = 0;
static constexpr int GDU = 1;
static constexpr int GDV = 2;
static constexpr int GDW = 3;
static constexpr int GDPRES = 4;
static constexpr int GDGAME = 5;

static constexpr int dComp_rhoD = 0;
static constexpr int dComp_mu = (dComp_rhoD + NUM_SPECIES);
static constexpr int dComp_xi = (dComp_mu + 1);
static constexpr int dComp_lambda = (dComp_xi + 1);

static constexpr int comp_Cs2 = 0;
static constexpr int comp_CI = (comp_Cs2 + 1);
static constexpr int comp_PrT = (comp_CI + 1);
static constexpr int comp_Cs2ovPrT = comp_PrT;
static constexpr int nCompC = (comp_PrT + 1);

static constexpr int NTHERM = 7;
static constexpr int QTHERM = 8;

static constexpr int NUM_ADV = 0;
static constexpr int NUM_AUX = 0;

static constexpr int UFS = (UFA + NUM_ADV);
static constexpr int UFX = (UFS + NUM_SPECIES);
static constexpr int ULIN = (UFX + NUM_AUX);
#ifdef SOOT_MODEL
static constexpr int NUM_LIN = (NUM_SOOT_MOMENTS + 1);
static constexpr int UFSOOT = ULIN;
#else
static constexpr int NUM_LIN = 0;
#endif
static constexpr int NVAR =
  (NTHERM + NUM_ADV + NUM_SPECIES + NUM_AUX + NUM_LIN);

static constexpr int QFS = (QFA + NUM_ADV);
static constexpr int QFX = (QFS + NUM_SPECIES);
static constexpr int QLIN = (QFX + NUM_AUX);
#ifdef SOOT_MODEL
static constexpr int QFSOOT = QLIN;
#endif
static constexpr int QVAR =
  (QTHERM + NUM_ADV + NUM_SPECIES + NUM_AUX + NUM_LIN);
static constexpr int NQAUX = 6;
static constexpr int NGDNV = 6;

#if AMREX_SPACEDIM == 1
// hopefully we wont do EB in 1D
static constexpr int NLSPTS = 1;
static constexpr int NEL_TRIMAT = 1;

#elif AMREX_SPACEDIM == 2

static constexpr int NLSPTS = 8;
static constexpr int NEL_TRIMAT = 3;

#elif AMREX_SPACEDIM == 3

static constexpr int NLSPTS = 26;
static constexpr int NEL_TRIMAT = 6;

#endif

#endif
