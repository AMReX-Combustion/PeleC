#include <AMReX_CArena.H>
#include <AMReX_REAL.H>
#include <AMReX_Amr.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_AmrLevel.H>
#include <AMReX_EB2.H>

// Defined and initialized when in gnumake, but not defined in cmake and
// initialization done manually
#ifndef AMREX_USE_SUNDIALS
#include <AMReX_Sundials.H>
#endif

#include "PeleC.H"

std::string inputs_name;

void
initialize_EB2(const amrex::Geometry& geom, int required_level, int max_level);

amrex::LevelBld* getLevelBld();

void
override_default_parameters()
{
  amrex::ParmParse pp("eb2");
  if (not pp.contains("geom_type")) {
    std::string geom_type("all_regular");
    pp.add("geom_type", geom_type);
  }
}

// FIXMEs:
// this is a stupid name
// put it somewhere else?
// need to do writeSmallPlotFile as well
// consolidate
// make it work with HDF5
// parmparse the HDF5 flag
// rename the bl_profile region
class myamr : public amrex::Amr
{
  using amrex::Amr::Amr;

public:
  void writePlotFile() override
  {
    const bool hdf5 = true;
    if (hdf5) {
      if (!Plot_Files_Output()) {
        return;
      }

      BL_PROFILE_REGION_START("Amr::writePlotFile()");
      BL_PROFILE("Amr::writePlotFile()");

      if (first_plotfile) {
        first_plotfile = false;
        amr_level[0]->setPlotVariables();
      }

      // Don't continue if we have no variables to plot.

      if (statePlotVars().size() == 0) {
        return;
      }

      const std::string& pltfile =
        amrex::Concatenate(plot_file_root, level_steps[0], file_name_digits);

      if (verbose > 0) {
        amrex::Print() << "PLOTFILE: file = " << pltfile << '\n';
      }

      if (record_run_info && amrex::ParallelDescriptor::IOProcessor()) {
        runlog << "PLOTFILE: file = " << pltfile << '\n';
      }

      // FIXME: make this common with the other
      const auto& desc_lst = amrex::AmrLevel::get_desc_lst();
      auto& derive_lst = amrex::AmrLevel::get_derive_lst();

      amrex::Vector<std::pair<int, int>> plot_var_map;
      for (int typ = 0; typ < desc_lst.size(); typ++) {
        for (int comp = 0; comp < desc_lst[typ].nComp(); comp++) {
          if (
            amrex::Amr::isStatePlotVar(desc_lst[typ].name(comp)) &&
            desc_lst[typ].getType() == amrex::IndexType::TheCellType()) {
            plot_var_map.push_back(std::pair<int, int>(typ, comp));
          }
        }
      }

      int num_derive = 0;
      std::list<std::string> derive_names;
      const std::list<amrex::DeriveRec>& dlist = derive_lst.dlist();

      for (const auto& it : dlist) {
        if (amrex::Amr::isDerivePlotVar(it.name())) {
          derive_names.push_back(it.name());
          num_derive += it.numDerive();
        }
      }

      const auto n_data_items = plot_var_map.size() + num_derive;

      int cnt = 0;
      const int nGrow = 0;
      const amrex::Real cur_time =
        (amr_level[0]->get_state_data(State_Type)).curTime();

      // // amrex::Vector<const amrex::MultiFab*> plotMFs(finestLevel() + 1);
      // for (int lev = 0; lev < finestLevel() + 1; ++lev) {
      //   amrex::MultiFab plotMF(
      //     boxArray(lev), DistributionMap(lev), n_data_items, nGrow,
      //     amrex::MFInfo(), amr_level[lev]->Factory());

      //   // Cull data from state variables -- use no ghost cells.
      //   for (int i = 0; i < plot_var_map.size(); i++) {
      //     int typ = plot_var_map[i].first;
      //     int comp = plot_var_map[i].second;
      //     amrex::MultiFab& this_dat = amr_level[lev]->get_new_data(typ);
      //     amrex::MultiFab::Copy(plotMF, this_dat, comp, cnt, 1, nGrow);
      //     cnt++;
      //   }

      //   // Cull data from derived variables.
      //   if (!derive_names.empty()) {
      //     for (const auto& derive_name : derive_names) {
      //       const amrex::DeriveRec* rec = derive_lst.get(derive_name);
      //       int ncomp = rec->numDerive();

      //       auto derive_dat =
      //         amr_level[lev]->derive(derive_name, cur_time, nGrow);
      //       amrex::MultiFab::Copy(plotMF, *derive_dat, 0, cnt, ncomp, nGrow);
      //       cnt += ncomp;
      //     }
      //   }
      // }

      amrex::Vector<std::string> plt_var_names;
      for (int i = 0; i < plot_var_map.size(); i++) {
        int typ = plot_var_map[i].first;
        int comp = plot_var_map[i].second;
        plt_var_names.push_back(desc_lst[typ].name(comp));
      }

      for (const auto& derive_name : derive_names) {
        const amrex::DeriveRec* rec = derive_lst.get(derive_name);
        for (int i = 0; i < rec->numDerive(); i++) {
          plt_var_names.push_back(rec->variableName(i));
        }
      }

      // amrex::WriteMultiLevelPlotfile(
      //   pltfile, finestLevel() + 1, plotMFs, plt_var_names,
      //   Geom(), cur_time, istep, refRatio());

      // if (level == 0 && amrex::ParallelDescriptor::IOProcessor()) {
      //   amr_level[0]->writeJobInfo(pltfile);
      // }
      
      BL_PROFILE_REGION_STOP("Amr::writePlotFile()");
    } else {
      amrex::Amr::writePlotFile();
    }
  }
};

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
  amrex::Initialize(
    argc, argv, true, MPI_COMM_WORLD, override_default_parameters);
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
  auto* amrptr = new myamr(getLevelBld());

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
