#ifndef INDEX_DEFINES_H
#define INDEX_DEFINES_H

#include "mechanism.H"

#define URHO 0
#define UMX 1
#define UMY 2
#define UMZ 3
#define UEDEN 4
#define UEINT 5
#define UTEMP 6
#define UFA 7

#define QRHO 0
#define QU 1
#define QV 2
#define QW 3
#define QGAME 4
#define QPRES 5
#define QREINT 6
#define QTEMP 7
#define QFA 8

#define QGAMC 0
#define QC 1
#define QCSML 2
#define QDPDR 3
#define QDPDE 4
#define QRSPEC 5

#define GDRHO 0
#define GDU 1
#define GDV 2
#define GDW 3
#define GDPRES 4
#define GDGAME 5

#define dComp_rhoD 0
#define dComp_mu (dComp_rhoD + NUM_SPECIES)
#define dComp_xi (dComp_mu + 1)
#define dComp_lambda (dComp_xi + 1)

#define comp_Cs2 0
#define comp_CI (comp_Cs2 + 1)
#define comp_PrT (comp_CI + 1)
#define comp_Cs2ovPrT comp_PrT
#define nCompC (comp_PrT + 1)

#define SlipWall 4
#define NoSlipWall 5

#define NTHERM 7
#define QTHERM 8

#define NUM_ADV 0
#define NUM_AUX 0

#define UFS (UFA + NUM_ADV)
#define UFX (UFS + NUM_SPECIES)
#define ULIN (UFX + NUM_AUX)
#ifdef SOOT_MODEL
#define NUM_LIN (NUM_SOOT_MOMENTS + 1)
#define UFSOOT ULIN
#else
#define NUM_LIN 0
#endif
#define NVAR (NTHERM + NUM_ADV + NUM_SPECIES + NUM_AUX + NUM_LIN)

#define QFS (QFA + NUM_ADV)
#define QFX (QFS + NUM_SPECIES)
#define QLIN (QFX + NUM_AUX)
#ifdef SOOT_MODEL
#define QFSOOT QLIN
#endif
#define QVAR (QTHERM + NUM_ADV + NUM_SPECIES + NUM_AUX + NUM_LIN)
#define NQAUX 6
#define NGDNV 6

#if AMREX_SPACEDIM == 1
// hopefully we wont do EB in 1D
#define NLSPTS 1
#define NEL_TRIMAT 1

#elif AMREX_SPACEDIM == 2

#define NLSPTS 8
#define NEL_TRIMAT 3

#elif AMREX_SPACEDIM == 3

#define NLSPTS 26
#define NEL_TRIMAT 6

#endif

#endif
