// System
#include <stdio.h>  // for sprintf
#include <string.h> // for strchr
#include <time.h>   // for tm, gmtime, time, time_t
#include <iomanip>  // for setw, __iom_t6, setfill, setpr...
#include <iostream> // for endl, string, operator>>, cin
#include <string>   // for operator==

// AMReX
#include "AMReX.H"                    // for Abort, Finalize, Initialize
#include "AMReX_Amr.H"                // for Amr
#include "AMReX_AmrLevel.H"           // for AmrLevel
#include "AMReX_Arena.H"              // for The_Arena, Arena
#include "AMReX_BLProfiler.H"         // for BL_PROFILE_SET_RUN_TIME, BL_PR...
#include "AMReX_EBSupport.H"          // for EBSupport, EBSupport::full
#include "AMReX_CArena.H"             // for CArena
#include "AMReX_OpenMP.H"             // for get_max_threads
#include "AMReX_ParmParse.H"          // for ParmParse
#include "AMReX_ParallelDescriptor.H" // for IOProcessor, second, Barrier
#include "AMReX_Print.H"              // for Print
#include "AMReX_REAL.H"               // for Real

// Defined and initialized when in gnumake, but not defined in cmake and
// initialization done manually
#ifndef AMREX_USE_SUNDIALS
#include "AMReX_Sundials.H" // for Finalize, Initialize
#endif

// PeleC
#include "PeleC.H" // for PeleC

namespace amrex {
class Geometry;
class LevelBld;
} // namespace amrex

std::string inputs_name;

void
initialize_EB2(const amrex::Geometry& geom, int required_level, int max_level);

amrex::LevelBld* getLevelBld();

int
main(int argc, char* argv[])
{
  // Use this to trap NaNs in C++
  // feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);

  if (argc <= 1) {
    amrex::Abort("Error: no inputs file provided on command line.");
  }

  // check to see if it contains --describe
  if (argc >= 2) {
    for (auto i = 1; i < argc; i++) {
      if (std::string(argv[i]) == "--describe") {
        PeleC::writeBuildInfo(std::cout);
        return 0;
      }
    }
  }

  // Make sure to catch new failures.
  amrex::Initialize(argc, argv);
// Defined and initialized when in gnumake, but not defined in cmake and
// initialization done manually
#ifndef AMREX_USE_SUNDIALS
  amrex::sundials::Initialize(amrex::OpenMP::get_max_threads());
#endif

  // Save the inputs file name for later.
  if (strchr(argv[1], '=') == nullptr) {
    inputs_name = argv[1];
  }

  BL_PROFILE_VAR("main()", pmain);

  amrex::Real dRunTime1 = amrex::ParallelDescriptor::second();

  amrex::Print() << std::setprecision(10);

  int max_step;
  amrex::Real strt_time;
  amrex::Real stop_time;
  amrex::ParmParse pp;

  bool pause_for_debug = false;
  pp.query("pause_for_debug", pause_for_debug);
  if (pause_for_debug) {
    if (amrex::ParallelDescriptor::IOProcessor()) {
      amrex::Print() << "Enter any string to continue" << std::endl;
      std::string text;
      std::cin >> text;
    }
    amrex::ParallelDescriptor::Barrier();
  }

  max_step = -1;
  strt_time = 0.0;
  stop_time = -1.0;

  pp.query("max_step", max_step);
  pp.query("strt_time", strt_time);
  pp.query("stop_time", stop_time);

  if (strt_time < 0.0) {
    amrex::Abort("MUST SPECIFY a non-negative strt_time");
  }

  if (max_step < 0 && stop_time < 0.0) {
    amrex::Abort(
      "Exiting because neither max_step nor stop_time is non-negative.");
  }

  // Print the current date and time
  time_t time_type;
  struct tm* time_pointer;
  time(&time_type);
  time_pointer = gmtime(&time_type);

  if (amrex::ParallelDescriptor::IOProcessor()) {
    amrex::Print() << std::setfill('0') << "\nStarting run at " << std::setw(2)
                   << time_pointer->tm_hour << ":" << std::setw(2)
                   << time_pointer->tm_min << ":" << std::setw(2)
                   << time_pointer->tm_sec << " UTC on "
                   << time_pointer->tm_year + 1900 << "-" << std::setw(2)
                   << time_pointer->tm_mon + 1 << "-" << std::setw(2)
                   << time_pointer->tm_mday << "." << std::endl;
  }

  // Initialize random seed after we're running in parallel.
  auto* amrptr = new amrex::Amr(getLevelBld());

  amrex::AmrLevel::SetEBSupportLevel(
    amrex::EBSupport::full); // need both area and volume fractions
  amrex::AmrLevel::SetEBMaxGrowCells(
    5, 5,
    5); // 5 focdr ebcellflags, 4 for vfrac, 2 is not used for EBSupport::volume
  initialize_EB2(
    amrptr->Geom(amrptr->maxLevel()), amrptr->maxLevel(), amrptr->maxLevel());

  amrptr->init(strt_time, stop_time);

  // If we set the regrid_on_restart flag and if we are *not* going to take
  // a time step then we want to go ahead and regrid here.
  if (
    amrptr->RegridOnRestart() &&
    ((amrptr->levelSteps(0) >= max_step) || (amrptr->cumTime() >= stop_time))) {
    // Regrid only!
    amrptr->RegridOnly(amrptr->cumTime());
  }

  amrex::Real dRunTime2 = amrex::ParallelDescriptor::second();

  while ((amrptr->okToContinue() != 0) &&
         (amrptr->levelSteps(0) < max_step || max_step < 0) &&
         (amrptr->cumTime() < stop_time || stop_time < 0.0)) {
    // Do a timestep
    amrptr->coarseTimeStep(stop_time);
  }

  // Write final checkpoint
  if (amrptr->stepOfLastCheckPoint() < amrptr->levelSteps(0)) {
    amrptr->checkPoint();
  }

  // Write final plotfile
  if (amrptr->stepOfLastPlotFile() < amrptr->levelSteps(0)) {
    amrptr->writePlotFile();
  }

  time(&time_type);
  time_pointer = gmtime(&time_type);

  if (amrex::ParallelDescriptor::IOProcessor()) {
    amrex::Print() << std::setfill('0') << "\nEnding run at " << std::setw(2)
                   << time_pointer->tm_hour << ":" << std::setw(2)
                   << time_pointer->tm_min << ":" << std::setw(2)
                   << time_pointer->tm_sec << " UTC on "
                   << time_pointer->tm_year + 1900 << "-" << std::setw(2)
                   << time_pointer->tm_mon + 1 << "-" << std::setw(2)
                   << time_pointer->tm_mday << "." << std::endl;
  }

  delete amrptr;

  // This MUST follow the above delete as ~Amr() may dump files to disk
  const int IOProc = amrex::ParallelDescriptor::IOProcessorNumber();

  amrex::Real dRunTime3 = amrex::ParallelDescriptor::second();

  amrex::Real runtime_total = dRunTime3 - dRunTime1;
  amrex::Real runtime_timestep = dRunTime3 - dRunTime2;

  amrex::ParallelDescriptor::ReduceRealMax(runtime_total, IOProc);
  amrex::ParallelDescriptor::ReduceRealMax(runtime_timestep, IOProc);

  if (amrex::ParallelDescriptor::IOProcessor()) {
    amrex::Print() << "Run time = " << runtime_total << std::endl;
    amrex::Print() << "Run time w/o init = " << runtime_timestep << std::endl;
  }

  if (auto* arena = dynamic_cast<amrex::CArena*>(amrex::The_Arena())) {
    // A barrier to make sure our output follows that of RunStats.
    amrex::ParallelDescriptor::Barrier();
    // We're using a CArena -- output some FAB memory stats.
    // This'll output total # of bytes of heap space in the Arena.
    // It's actually the high water mark of heap space required by FABs.
    char buf[256];

    sprintf(
      buf, "CPU(%d): Heap Space (bytes) used by Coalescing FAB Arena: %zu",
      amrex::ParallelDescriptor::MyProc(), arena->heap_space_used());

    amrex::Print() << buf << std::endl;
  }

  BL_PROFILE_VAR_STOP(pmain);
  BL_PROFILE_SET_RUN_TIME(dRunTime2);

// Defined and finalized when in gnumake, but not defined in cmake and
// finalization done manually
#ifndef AMREX_USE_SUNDIALS
  amrex::sundials::Finalize();
#endif
  amrex::Finalize();

  return 0;
}
