
#include <new>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <iomanip>

#ifndef WIN32
#include <unistd.h>
#endif

#include <AMReX_CArena.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>
#include <AMReX_IntVect.H>
#include <AMReX_Box.H>
#include <AMReX_Amr.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_AmrLevel.H>
#include <AMReX_VisMF.H>

#include <time.h>

#include <PeleC.H>

using namespace amrex;

std::string inputs_name = "";

int
main (int   argc,
      char* argv[])
{

    //
    // Make sure to catch new failures.
    //
    amrex::Initialize(argc,argv);

    // Refuse to continue if we did not provide an inputs file.

    if (argc <= 1) {
	amrex::Abort("Error: no inputs file provided on command line.");
    }

    // Save the inputs file name for later.

    if (!strchr(argv[1], '=')) {
	inputs_name = argv[1];
    }

    BL_PROFILE_VAR("main()", pmain);

    std::cout << std::setprecision(10);

    Amr* amrptr = new Amr;

    Real strt_time = 0;
    Real stop_time = 1;
    amrptr->init(strt_time,stop_time);

    PeleC& p0 = dynamic_cast<PeleC&>(amrptr->getLevel(0));
    PeleC& p1 = dynamic_cast<PeleC&>(amrptr->getLevel(1));
    PeleC& p2 = dynamic_cast<PeleC&>(amrptr->getLevel(2));

    MultiFab& S0 = p0.get_new_data(State_Type);
    MultiFab& S1 = p1.get_new_data(State_Type);
    MultiFab& S2 = p2.get_new_data(State_Type);

    VisMF::Write(S0,"S0");
    VisMF::Write(S1,"S1");
    VisMF::Write(S2,"S2");

    int ns=S0.nComp();

    MultiFab D0(S0.boxArray(),S0.DistributionMap(),ns,1);
    MultiFab D1(S1.boxArray(),S1.DistributionMap(),ns,1);
    MultiFab D2(S2.boxArray(),S2.DistributionMap(),ns,1);
    int is_old = 0;

    p0.getDiffusionTerm(strt_time,D0,is_old,1,1,1,1);
    p1.getDiffusionTerm(strt_time,D1,is_old,1,1,1,1);
    p2.getDiffusionTerm(strt_time,D2,is_old,1,1,1,1);

    VisMF::Write(D0,"D0");
    VisMF::Write(D1,"D1");
    VisMF::Write(D2,"D2");

    // Use PeleCs avgDown operator to average down D2 to D1 and get diff, T1
    //   and to average down D1 to D0 and get diff, T0

    MultiFab T0(D0.boxArray(),D0.DistributionMap(),ns,0);
    MultiFab::Copy(T0,D0,0,0,ns,0);
    MultiFab::Copy(S1,D1,0,0,ns,0);
    p0.avgDown(State_Type);
    MultiFab::Subtract(T0,S0,0,0,ns,0);

    MultiFab T1(S1.boxArray(),S1.DistributionMap(),ns,0);
    MultiFab::Copy(T1,D1,0,0,ns,0);
    MultiFab::Copy(S2,D2,0,0,ns,0);
    p1.avgDown(State_Type);
    MultiFab::Subtract(T1,S1,0,0,ns,0);

    Real vol0 = 1;
    Real vol1 = 1;
    Real vol2 = 1;
    for (int d=0; d<BL_SPACEDIM; ++d)
    {
      vol0 *= p0.Geom().CellSize()[d];
      vol1 *= p1.Geom().CellSize()[d];
      vol2 *= p2.Geom().CellSize()[d];
    }

    int nc=PeleC::Eden;
    std::cout << "Eden: \n";
    std::cout << "   order0: " << std::log2(T0.norm0(nc)/T1.norm0(nc)) << std::endl;
    std::cout << "   order1: " << std::log2((vol0*T0.norm1(nc))/(vol1*T1.norm1(nc))) << std::endl;
    std::cout << "   order2: " << std::log2(std::sqrt(vol0/vol1)*T0.norm2(nc)/T1.norm2(nc)) << std::endl;

    nc=PeleC::Xmom;
    std::cout << "Xmom: \n";
    std::cout << "   order0: " << std::log2(T0.norm0(nc)/T1.norm0(nc)) << std::endl;
    std::cout << "   order1: " << std::log2((vol0*T0.norm1(nc))/(vol1*T1.norm1(nc))) << std::endl;
    std::cout << "   order2: " << std::log2(std::sqrt(vol0/vol1)*T0.norm2(nc)/T1.norm2(nc)) << std::endl;

    if (AMREX_SPACEDIM > 1)
    {
      nc=PeleC::Ymom;
      std::cout << "Ymom: \n";
      std::cout << "   order0: " << std::log2(T0.norm0(nc)/T1.norm0(nc)) << std::endl;
      std::cout << "   order1: " << std::log2((vol0*T0.norm1(nc))/(vol1*T1.norm1(nc))) << std::endl;
      std::cout << "   order2: " << std::log2(std::sqrt(vol0/vol1)*T0.norm2(nc)/T1.norm2(nc)) << std::endl;
    }

    if (AMREX_SPACEDIM > 2)
    {
      nc=PeleC::Zmom;
      std::cout << "Zmom: \n";
      std::cout << "   order0: " << std::log2(T0.norm0(nc)/T1.norm0(nc)) << std::endl;
      std::cout << "   order1: " << std::log2((vol0*T0.norm1(nc))/(vol1*T1.norm1(nc))) << std::endl;
      std::cout << "   order2: " << std::log2(std::sqrt(vol0/vol1)*T0.norm2(nc)/T1.norm2(nc)) << std::endl;
    }

    VisMF::Write(T0,"T0");
    VisMF::Write(T1,"T1");

    amrex::Finalize();

    return 0;
}
