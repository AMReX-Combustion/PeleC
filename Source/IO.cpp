#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <ctime>

#ifdef AMREX_USE_OMP
#include <omp.h>
#endif

#include <AMReX_Utility.H>
#include <AMReX_buildInfo.H>
#include <AMReX_ParmParse.H>
#include <AMReX_EBMultiFabUtil.H>

#include "mechanism.H"
#include "PltFileManager.H"

#include "PeleC.H"
#include "IO.H"
#include "IndexDefines.H"

#ifdef PELEC_USE_SPRAY
#include "SprayParticles.H"
#endif

// PeleC maintains an internal checkpoint version numbering system.
// This allows us to maintain backwards compatibility with checkpoints
// generated by old versions of the code, so that new versions can
// restart from them. The version number is stored in the PeleCHeader
// file inside a checkpoint. The following were the changes that were made
// in updating version numbers:
// 0: all checkpoints as of 11/21/16
// 1: add body state

namespace {
int input_version = -1;
int current_version = 1;
std::string body_state_filename = "body_state.fab";
} // namespace

// I/O routines for PeleC

bool
PeleC::check_state_in_checkpoint(const StateType state_type)
{
  const std::string state_pfx =
    "Level_" + std::to_string(level) + "/SD_" + std::to_string(state_type);

  const std::string filename = parent->theRestartFile();
  const std::string faHeaderFilesName(filename + "/FabArrayHeaders.txt");
  amrex::Vector<char> faHeaderFileChars;
  bool bExitOnError(false); // ---- dont exit if this file does not exist
  amrex::ParallelDescriptor::ReadAndBcastFile(
    faHeaderFilesName, faHeaderFileChars, bExitOnError);
  if (!faHeaderFileChars.empty()) { // ---- headers were read
    std::string faFileCharPtrString(faHeaderFileChars.dataPtr());
    std::istringstream fais(faFileCharPtrString, std::istringstream::in);
    while (!fais.eof()) {
      std::string faHeaderName;
      fais >> faHeaderName;
      if (!fais.eof()) {
        if (faHeaderName.rfind(state_pfx, 0) == 0) {
          return true;
        }
      }
    }
  }
  return false;
}

void
PeleC::restart(amrex::Amr& papa, std::istream& is, bool bReadSpecial)
{
  // Let's check PeleC checkpoint version first;
  // trying to read from checkpoint; if nonexisting, set it to 0.
  if (input_version == -1) {
    if (amrex::ParallelDescriptor::IOProcessor()) {
      std::ifstream PeleCHeaderFile;
      std::string FullPathPeleCHeaderFile = papa.theRestartFile();
      FullPathPeleCHeaderFile += "/PeleCHeader";
      PeleCHeaderFile.open(FullPathPeleCHeaderFile.c_str(), std::ios::in);
      if (PeleCHeaderFile.good()) {
        char foo[256];
        // first line: Checkpoint version: ?
        PeleCHeaderFile.getline(foo, 256, ':');
        PeleCHeaderFile >> input_version;
        PeleCHeaderFile.close();
      } else {
        input_version = 0;
      }
    }
    amrex::ParallelDescriptor::Bcast(
      &input_version, 1, amrex::ParallelDescriptor::IOProcessorNumber());
  }

  AMREX_ASSERT(input_version >= 0);

  // Also need to mod checkPoint function to store the new version in a text
  // file
  AmrLevel::restart(papa, is, bReadSpecial);

  // Deal here with new state descriptor types added, with corresponding
  // input_version > 0, if applicable
  amrex::Vector<int> state_in_checkpoint(desc_lst.size(), 1);
  set_state_in_checkpoint(state_in_checkpoint);
  for (int i = 0; i < desc_lst.size(); ++i) {
    if (state_in_checkpoint[i] == 0) {
      const amrex::Real ctime = state[i - 1].curTime();
      state[i].define(
        geom.Domain(), grids, dmap, desc_lst[i], ctime, parent->dtLevel(level),
        Factory());
      get_new_data(i).setVal(0.0);
    }
  }
  buildMetrics();

  init_eb();

  const amrex::MultiFab& S_new = get_new_data(State_Type);

  for (int n = 0; n < src_list.size(); ++n) {
    int oldGrow = numGrow();
    int newGrow = S_new.nGrow();
    old_sources[src_list[n]] = std::make_unique<amrex::MultiFab>(
      grids, dmap, NVAR, oldGrow, amrex::MFInfo(), Factory());
    new_sources[src_list[n]] = std::make_unique<amrex::MultiFab>(
      grids, dmap, NVAR, newGrow, amrex::MFInfo(), Factory());
  }

  if (do_hydro || do_diffuse) {
    Sborder.define(grids, dmap, NVAR, numGrow(), amrex::MFInfo(), Factory());
  }

  if (!do_mol) {
    if (do_hydro) {
      hydro_source.define(grids, dmap, NVAR, 0, amrex::MFInfo(), Factory());

      // This array holds the sum of all source terms that affect the
      // hydrodynamics. If we are doing the source term predictor, we'll also
      // use this after the hydro update to store the sum of the new-time
      // sources, so that we can compute the time derivative of the source
      // terms.
      sources_for_hydro.define(
        grids, dmap, NVAR, numGrow(), amrex::MFInfo(), Factory());
    }
  } else {
    Sborder.define(grids, dmap, NVAR, numGrow(), amrex::MFInfo(), Factory());
  }

  // get the elapsed CPU time to now
  if (level == 0 && amrex::ParallelDescriptor::IOProcessor()) {
    std::ifstream CPUFile;
    std::string FullPathCPUFile = parent->theRestartFile();
    FullPathCPUFile += "/CPUtime";
    CPUFile.open(FullPathCPUFile.c_str(), std::ios::in);

    CPUFile >> previousCPUTimeUsed;
    CPUFile.close();

    amrex::Print() << "read CPU time: " << previousCPUTimeUsed << "\n";
  }

  if (track_grid_losses && level == 0) {

    // get the current value of the diagnostic quantities
    std::ifstream DiagFile;
    std::string FullPathDiagFile = parent->theRestartFile();
    FullPathDiagFile += "/Diagnostics";
    DiagFile.open(FullPathDiagFile.c_str(), std::ios::in);

    for (amrex::Real& i : material_lost_through_boundary_cumulative) {
      DiagFile >> i;
    }

    DiagFile.close();
  }

  if (level > 0 && do_reflux) {
    flux_reg = std::make_unique<amrex::EBFluxRegister>(
      grids, papa.boxArray(level - 1), dmap, papa.DistributionMap(level - 1),
      geom, papa.Geom(level - 1), papa.refRatio(level - 1), level, NVAR);

    if (!amrex::DefaultGeometry().IsCartesian()) {
      // pres_reg.define(
      // grids, papa.boxArray(level - 1), dmap, papa.DistributionMap(level - 1),
      // geom, papa.Geom(level - 1), papa.refRatio(level - 1), level, 1);
      amrex::Abort("We don't do rz.");
    }
  }

  if (input_version > 0 && level == 0 && eb_in_domain) {
    if (amrex::ParallelDescriptor::IOProcessor()) {
      std::ifstream BodyFile;
      std::string FullPathBodyFile = parent->theRestartFile();
      FullPathBodyFile += "/" + body_state_filename;
      BodyFile.open(FullPathBodyFile.c_str(), std::ios::in);
      amrex::FArrayBox bstate_fab;
      bstate_fab.readFrom(BodyFile);
      BodyFile.close();
      if (bstate_fab.nComp() != NVAR) {
        amrex::Abort("Body state incompatible with checkpointed version");
      }
      const auto bx = bstate_fab.box();
      auto const& bstate_arr = bstate_fab.array();
      amrex::Gpu::DeviceVector<amrex::Real> local_body_state(NVAR, -1);
      amrex::Real* p_body_state = local_body_state.begin();
      amrex::ParallelFor(
        bx, NVAR, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
          p_body_state[n] = bstate_arr(i, j, k, n);
        });

      amrex::Gpu::copy(
        amrex::Gpu::deviceToHost, local_body_state.begin(),
        local_body_state.end(), body_state.begin());
    }
    amrex::ParallelDescriptor::Bcast(
      &(body_state[0]), body_state.size(), // NOLINT
      amrex::ParallelDescriptor::IOProcessorNumber());
    body_state_set = true;
  }
}

void
PeleC::set_state_in_checkpoint(amrex::Vector<int>& state_in_checkpoint)
{
  for (int i = 0; i < num_state_type; ++i) {
    const bool is_present =
      check_state_in_checkpoint(static_cast<StateType>(i));

    if (i == State_Type) {
      state_in_checkpoint[i] = 1;
      if ((!is_present) && (level == 0)) {
        amrex::Abort("State_Type is not present in the checkpoint file");
      }
    } else if (i == Reactions_Type) {
      if (!do_react) {
        state_in_checkpoint[i] = 0;
      } else {
        state_in_checkpoint[i] = is_present ? 1 : 0;
      }
    } else if (i == Work_Estimate_Type) {
      // Never use work estimate checkpoint
      state_in_checkpoint[i] = 0;
    } else {
      amrex::Abort("Unknown StateType");
    }
  }
}

void
PeleC::checkPoint(
  const std::string& dir,
  std::ostream& os,
  amrex::VisMF::How how,
  bool /*dump_old_default*/)
{
  amrex::AmrLevel::checkPoint(dir, os, how, dump_old);

#ifdef PELEC_USE_SPRAY
  if (SprayPC != nullptr) {
    SprayPC->SprayParticleIO(level, true, 0, dir);
  }
#endif

  if (level == 0 && amrex::ParallelDescriptor::IOProcessor()) {
    {
      std::ofstream PeleCHeaderFile;
      std::string FullPathPeleCHeaderFile = dir;
      FullPathPeleCHeaderFile += "/PeleCHeader";
      PeleCHeaderFile.open(FullPathPeleCHeaderFile.c_str(), std::ios::out);

      PeleCHeaderFile << "Checkpoint version: " << current_version << std::endl;
      PeleCHeaderFile.close();
    }

    {
      // store elapsed CPU time
      std::ofstream CPUFile;
      std::string FullPathCPUFile = dir;
      FullPathCPUFile += "/CPUtime";
      CPUFile.open(FullPathCPUFile.c_str(), std::ios::out);

      CPUFile << std::setprecision(15) << getCPUTime();
      CPUFile.close();
    }

    {
      // Store the level at which EB was generated
      std::ofstream EBLevelFile;
      std::string FullPathEBLevelFile = dir;
      FullPathEBLevelFile += "/EBMaxLevel";
      EBLevelFile.open(FullPathEBLevelFile.c_str(), std::ios::out);
      EBLevelFile << eb_max_lvl_gen;
      EBLevelFile.close();
    }

    if (track_grid_losses) {
      // store diagnostic quantities
      std::ofstream DiagFile;
      std::string FullPathDiagFile = dir;
      FullPathDiagFile += "/Diagnostics";
      DiagFile.open(FullPathDiagFile.c_str(), std::ios::out);

      for (amrex::Real i : material_lost_through_boundary_cumulative) {
        DiagFile << std::setprecision(15) << i << std::endl;
      }

      DiagFile.close();
    }
  }

  if (current_version > 0) {
    if (amrex::ParallelDescriptor::IOProcessor() && eb_in_domain) {
      amrex::IntVect iv(AMREX_D_DECL(0, 0, 0));
      const amrex::Box bx(iv, iv);
      amrex::FArrayBox bstate_fab(bx, NVAR);
      auto const& bstate_arr = bstate_fab.array();
      auto const captured_body_state = body_state;
      amrex::ParallelFor(
        bx, NVAR, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
          bstate_arr(i, j, k, n) = captured_body_state[n];
        });

      std::ofstream BodyFile;
      std::string FullPathBodyFile = dir;
      FullPathBodyFile += "/" + body_state_filename;
      BodyFile.open(FullPathBodyFile.c_str(), std::ios::out);
      bstate_fab.writeOn(BodyFile);
      BodyFile.close();
    }
  }
}

void
PeleC::setPlotVariables()
{
  amrex::AmrLevel::setPlotVariables();

  amrex::ParmParse pp("pelec");

  bool plot_cost = do_react_load_balance || do_mol_load_balance;
  pp.query("plot_cost", plot_cost);
  if (plot_cost) {
    for (int i = 0; i < desc_lst[Work_Estimate_Type].nComp(); i++) {
      amrex::Amr::addStatePlotVar(desc_lst[Work_Estimate_Type].name(i));
    }
  } else {
    for (int i = 0; i < desc_lst[Work_Estimate_Type].nComp(); i++) {
      amrex::Amr::deleteStatePlotVar(desc_lst[Work_Estimate_Type].name(i));
    }
  }

  if (!do_react) {
    for (int i = 0; i < desc_lst[Reactions_Type].nComp(); i++) {
      amrex::Amr::deleteStatePlotVar(desc_lst[Reactions_Type].name(i));
    }
  }

  bool plot_rhoy = true;
  pp.query("plot_rhoy", plot_rhoy);
  if (plot_rhoy) {
    for (int i = 0; i < NUM_SPECIES; i++) {
      amrex::Amr::addStatePlotVar(desc_lst[State_Type].name(FirstSpec + i));
    }
  } else {
    for (int i = 0; i < NUM_SPECIES; i++) {
      amrex::Amr::deleteStatePlotVar(desc_lst[State_Type].name(FirstSpec + i));
    }
  }

#ifdef PELEC_USE_SOOT
  if (plot_soot && add_soot_src) {
    for (int i = 0; i < NumSootVars; i++) {
      amrex::Amr::addStatePlotVar(desc_lst[State_Type].name(FirstSootVar + i));
    }
  } else {
    for (int i = 0; i < NumSootVars; i++) {
      amrex::Amr::deleteStatePlotVar(
        desc_lst[State_Type].name(FirstSootVar + i));
    }
  }
#endif
  int plot_reactions = 1;
  pp.query("plot_reactions", plot_reactions);
  if (plot_reactions == 0) {
    for (int i = 0; i < NUM_SPECIES + 1; i++) {
      amrex::Amr::deleteStatePlotVar(desc_lst[Reactions_Type].name(i));
    }
  }

  bool plot_massfrac = false;
  pp.query("plot_massfrac", plot_massfrac);
  if (plot_massfrac) {
    amrex::Amr::addDerivePlotVar("massfrac");
  } else {
    amrex::Amr::deleteDerivePlotVar("massfrac");
  }

#if NUM_ADV > 0
  bool plot_rho_adv = true;
  pp.query("plot_rho_adv", plot_rho_adv);
  if (plot_rho_adv) {
    for (int i = 0; i < NUM_ADV; i++) {
      amrex::Amr::addStatePlotVar(desc_lst[State_Type].name(FirstAdv + i));
    }
  } else {
    for (int i = 0; i < NUM_ADV; i++) {
      amrex::Amr::deleteStatePlotVar(desc_lst[State_Type].name(FirstAdv + i));
    }
  }

  bool plot_adv = false;
  pp.query("plot_adv", plot_adv);
  if (plot_adv) {
    amrex::Amr::addDerivePlotVar("adv");
  } else {
    amrex::Amr::deleteDerivePlotVar("adv");
  }
#endif

  bool plot_moleFrac = false;
  pp.query("plot_molefrac", plot_moleFrac);
  if (plot_moleFrac) {
    amrex::Amr::addDerivePlotVar("molefrac");
  } else {
    amrex::Amr::deleteDerivePlotVar("molefrac");
  }
}

void
PeleC::writePlotFilePost(const std::string& dir, std::ostream& /*os*/)
{
  if (level == 0 && amrex::ParallelDescriptor::IOProcessor()) {
    writeJobInfo(dir);
  }
}

void
PeleC::writeJobInfo(const std::string& dir)
{
  // job_info file with details about the run
  std::ofstream jobInfoFile;
  std::string FullPathJobInfoFile = dir;
  FullPathJobInfoFile += "/job_info";
  jobInfoFile.open(FullPathJobInfoFile.c_str(), std::ios::out);

  std::string PrettyLine = "==================================================="
                           "============================\n";
  std::string OtherLine = "----------------------------------------------------"
                          "----------------------------\n";
  std::string SkipSpace = "        ";

  // job information
  jobInfoFile << PrettyLine;
  jobInfoFile << " PeleC Job Information\n";
  jobInfoFile << PrettyLine;

  jobInfoFile << "job name: " << job_name << "\n\n";
  jobInfoFile << "inputs file: " << inputs_name << "\n\n";

  jobInfoFile << "number of MPI processes: "
              << amrex::ParallelDescriptor::NProcs() << "\n";
#ifdef AMREX_USE_OMP
  jobInfoFile << "number of threads:       " << omp_get_max_threads() << "\n";
#endif

  jobInfoFile << "\n";
  jobInfoFile << "CPU time used since start of simulation (CPU-hours): "
              << getCPUTime() / 3600.0;

  jobInfoFile << "\n\n";

  // plotfile information
  jobInfoFile << PrettyLine;
  jobInfoFile << " Plotfile Information\n";
  jobInfoFile << PrettyLine;

  // Convert now to tm struct for local timezone
  time_t now = time(nullptr);
  char time_buffer[128];
  struct tm localtm;
  localtime_r(&now, &localtm);
  strftime(time_buffer, sizeof(time_buffer), "%b %d %Y %H:%M:%S", &localtm);
  jobInfoFile << "output data / time: " << time_buffer << std::endl;

  std::string currentDir = amrex::FileSystem::CurrentPath();
  jobInfoFile << "output dir:         " << currentDir << "\n";

  jobInfoFile << "\n\n";

  // build information
  jobInfoFile << PrettyLine;
  jobInfoFile << " Build Information\n";
  jobInfoFile << PrettyLine;

  jobInfoFile << "build date:    " << amrex::buildInfoGetBuildDate() << "\n";
  jobInfoFile << "build machine: " << amrex::buildInfoGetBuildMachine() << "\n";
  jobInfoFile << "build dir:     " << amrex::buildInfoGetBuildDir() << "\n";
  jobInfoFile << "AMReX dir:     " << amrex::buildInfoGetAMReXDir() << "\n";

  jobInfoFile << "\n";

  jobInfoFile << "COMP:          " << amrex::buildInfoGetComp() << "\n";
  jobInfoFile << "COMP version:  " << amrex::buildInfoGetCompVersion() << "\n";
  jobInfoFile << "FCOMP:         " << amrex::buildInfoGetFcomp() << "\n";
  jobInfoFile << "FCOMP version: " << amrex::buildInfoGetFcompVersion() << "\n";

  jobInfoFile << "\n";

  for (int n = 1; n <= amrex::buildInfoGetNumModules(); n++) {
    jobInfoFile << amrex::buildInfoGetModuleName(n) << ": "
                << amrex::buildInfoGetModuleVal(n) << "\n";
  }

  jobInfoFile << "\n";

  const char* githash1 = amrex::buildInfoGetGitHash(1);
  const char* githash2 = amrex::buildInfoGetGitHash(2);
  const char* githash3 = amrex::buildInfoGetGitHash(3);
  const char* githash4 = amrex::buildInfoGetGitHash(4);
  const char* githash5 = amrex::buildInfoGetGitHash(5);
  if (strlen(githash1) > 0) {
    jobInfoFile << "PeleC       git hash: " << githash1 << "\n";
  }
  if (strlen(githash2) > 0) {
    jobInfoFile << "AMReX       git hash: " << githash2 << "\n";
  }
  if (strlen(githash3) > 0) {
    jobInfoFile << "PelePhysics git hash: " << githash3 << "\n";
  }
  if (strlen(githash4) > 0) {
    jobInfoFile << "AMReX-Hydro git hash: " << githash4 << "\n";
  }
  if (strlen(githash5) > 0) {
    jobInfoFile << "SUNDIALS    git hash: " << githash5 << "\n";
  }

  const char* buildgithash = amrex::buildInfoGetBuildGitHash();
  const char* buildgitname = amrex::buildInfoGetBuildGitName();
  if (strlen(buildgithash) > 0) {
    jobInfoFile << buildgitname << " git hash: " << buildgithash << "\n";
  }

  jobInfoFile << "\n\n";

  // grid information
  jobInfoFile << PrettyLine;
  jobInfoFile << " Grid Information\n";
  jobInfoFile << PrettyLine;

  int f_lev = parent->finestLevel();

  for (int i = 0; i <= f_lev; i++) {
    jobInfoFile << " level: " << i << "\n";
    jobInfoFile << "   number of boxes = " << parent->numGrids(i) << "\n";
    jobInfoFile << "   maximum zones   = ";
    for (int n = 0; n < AMREX_SPACEDIM; n++) {
      jobInfoFile << parent->Geom(i).Domain().length(n) << " ";
    }
    jobInfoFile << "\n\n";
  }

  jobInfoFile << " Boundary conditions\n";
  amrex::Vector<std::string> lo_bc_out(AMREX_SPACEDIM);
  amrex::Vector<std::string> hi_bc_out(AMREX_SPACEDIM);
  amrex::ParmParse pp("pelec");
  pp.getarr("lo_bc", lo_bc_out, 0, AMREX_SPACEDIM);
  pp.getarr("hi_bc", hi_bc_out, 0, AMREX_SPACEDIM);

  // these names correspond to the integer flags setup in the
  // Setup.cpp

  jobInfoFile << "   -x: " << lo_bc_out[0] << "\n";
  jobInfoFile << "   +x: " << hi_bc_out[0] << "\n";
  if (AMREX_SPACEDIM >= 2) {
    jobInfoFile << "   -y: " << lo_bc_out[1] << "\n";
    jobInfoFile << "   +y: " << hi_bc_out[1] << "\n";
  }
  if (AMREX_SPACEDIM == 3) {
    jobInfoFile << "   -z: " << lo_bc_out[2] << "\n";
    jobInfoFile << "   +z: " << hi_bc_out[2] << "\n";
  }

  jobInfoFile << "\n\n";

  const int mlen = 20;

  jobInfoFile << PrettyLine;
  jobInfoFile << " Species Information\n";
  jobInfoFile << PrettyLine;

  jobInfoFile << std::setw(6) << "index" << SkipSpace << std::setw(mlen + 1)
              << "name" << SkipSpace << std::setw(7) << "A" << SkipSpace
              << std::setw(7) << "Z"
              << "\n";
  jobInfoFile << OtherLine;
  jobInfoFile << "\n\n";

  // runtime parameters
  jobInfoFile << PrettyLine;
  jobInfoFile << " Inputs File Parameters\n";
  jobInfoFile << PrettyLine;

  amrex::ParmParse::dumpTable(jobInfoFile, true);
  jobInfoFile.close();
}

// PeleC::writeBuildInfo
// Similar to writeJobInfo, but the subset of information that makes sense
// without an input file to enable --describe in format similar to CASTRO
void
PeleC::writeBuildInfo(std::ostream& os)
{
  std::string PrettyLine = std::string(78, '=') + "\n";

  // build information
  os << PrettyLine;
  os << " PeleC Build Information\n";
  os << PrettyLine;

  os << "build date:    " << amrex::buildInfoGetBuildDate() << "\n";
  os << "build machine: " << amrex::buildInfoGetBuildMachine() << "\n";
  os << "build dir:     " << amrex::buildInfoGetBuildDir() << "\n";
  os << "AMReX dir:     " << amrex::buildInfoGetAMReXDir() << "\n";

  os << "\n";

  os << "COMP:          " << amrex::buildInfoGetComp() << "\n";
  os << "COMP version:  " << amrex::buildInfoGetCompVersion() << "\n";

  amrex::Print() << "C++ compiler:  " << amrex::buildInfoGetCXXName() << "\n";
  amrex::Print() << "C++ flags:     " << amrex::buildInfoGetCXXFlags() << "\n";

  os << "\n";

  os << "FCOMP:         " << amrex::buildInfoGetFcomp() << "\n";
  os << "FCOMP version: " << amrex::buildInfoGetFcompVersion() << "\n";

  os << "\n";

  amrex::Print() << "Link flags:    " << amrex::buildInfoGetLinkFlags() << "\n";
  amrex::Print() << "Libraries:     " << amrex::buildInfoGetLibraries() << "\n";

  os << "\n";

  for (int n = 1; n <= amrex::buildInfoGetNumModules(); n++) {
    os << amrex::buildInfoGetModuleName(n) << ": "
       << amrex::buildInfoGetModuleVal(n) << "\n";
  }

  os << "\n";
  const char* githash1 = amrex::buildInfoGetGitHash(1);
  const char* githash2 = amrex::buildInfoGetGitHash(2);
  const char* githash3 = amrex::buildInfoGetGitHash(3);
  const char* githash4 = amrex::buildInfoGetGitHash(4);
  const char* githash5 = amrex::buildInfoGetGitHash(5);
  if (strlen(githash1) > 0) {
    os << "PeleC       git hash: " << githash1 << "\n";
  }
  if (strlen(githash2) > 0) {
    os << "AMReX       git hash: " << githash2 << "\n";
  }
  if (strlen(githash3) > 0) {
    os << "PelePhysics git hash: " << githash3 << "\n";
  }
  if (strlen(githash4) > 0) {
    os << "AMReX-Hydro git hash: " << githash4 << "\n";
  }
  if (strlen(githash5) > 0) {
    os << "SUNDIALS    git hash: " << githash5 << "\n";
  }

  const char* buildgithash = amrex::buildInfoGetBuildGitHash();
  const char* buildgitname = amrex::buildInfoGetBuildGitName();
  if (strlen(buildgithash) > 0) {
    os << buildgitname << " git hash: " << buildgithash << "\n";
  }

  os << "\n";
  os << " PeleC Compile time variables: \n";

  int mm;
  int kk;
  int ii;
  int nfit;
  CKINDX(mm, kk, ii, nfit);
  os << std::setw(40) << std::left << "Number elements from chem cpp : " << mm
     << std::endl;
  os << std::setw(40) << std::left << "Number species from chem cpp : " << kk
     << std::endl;
  os << std::setw(40) << std::left << "Number reactions from chem cpp : " << ii
     << std::endl;

  os << "\n";
  os << " PeleC Defines: \n";
#ifdef AMREX_USE_OMP
  os << std::setw(35) << std::left << "AMREX_USE_OMP " << std::setw(6) << "ON"
     << std::endl;
#else
  os << std::setw(35) << std::left << "AMREX_USE_OMP " << std::setw(6) << "OFF"
     << std::endl;
#endif

#ifdef MPI_VERSION
  os << std::setw(35) << std::left << "MPI_VERSION " << std::setw(6)
     << MPI_VERSION << std::endl;
#else
  os << std::setw(35) << std::left << "MPI_VERSION " << std::setw(6)
     << "UNDEFINED" << std::endl;
#endif

#ifdef MPI_SUBVERSION
  os << std::setw(35) << std::left << "MPI_SUBVERSION " << std::setw(6)
     << MPI_SUBVERSION << std::endl;
#else
  os << std::setw(35) << std::left << "MPI_SUBVERSION " << std::setw(6)
     << "UNDEFINED" << std::endl;
#endif

  os << std::setw(35) << std::left << "NUM_ADV=" << NUM_ADV << std::endl;

  os << std::setw(35) << std::left << "NUM_AUX=" << NUM_AUX << std::endl;

  os << std::setw(35) << std::left << "NUM_LIN=" << NUM_LIN << std::endl;

#ifdef PELEC_USE_MASA
  os << std::setw(35) << std::left << "PELEC_USE_MASA " << std::setw(6) << "ON"
     << std::endl;
#else
  os << std::setw(35) << std::left << "PELEC_USE_MASA " << std::setw(6) << "OFF"
     << std::endl;
#endif

#ifdef PELEC_USE_SPRAY
  os << std::setw(35) << std::left << "PELEC_USE_SPRAY " << std::setw(6) << "ON"
     << std::endl;
#else
  os << std::setw(35) << std::left << "PELEC_USE_SPRAY " << std::setw(6)
     << "OFF" << std::endl;
#endif
#ifdef PELEC_USE_SOOT
  os << std::setw(35) << std::left << "PELEC_USE_SOOT " << std::setw(6) << "ON"
     << std::endl;
#else
  os << std::setw(35) << std::left << "PELEC_USE_SOOT " << std::setw(6) << "OFF"
     << std::endl;
#endif

  os << "\n\n";
}

void
PeleC::initLevelDataFromPlt(
  const int lev, const std::string& dataPltFile, amrex::MultiFab& S_new)
{
  amrex::Print() << "Using data (rho, u, T, Y) from pltfile " << dataPltFile
                 << std::endl;
  pele::physics::pltfilemanager::PltFileManager pltData(dataPltFile);
  const auto plt_vars = pltData.getVariableList();

  // Read rho, u, temperature (required)
  std::map<std::string, int> vars{
    {"density", -1}, {"x_velocity", -1}, {"Temp", -1}};
  for (auto& var : vars) {
    var.second = find_position(plt_vars, var.first);
    if (var.second == -1) {
      amrex::Abort("Unable to find variable in plot file: " + var.first);
    }
  }
  pltData.fillPatchFromPlt(lev, geom, vars["density"], URHO, 1, S_new);
  pltData.fillPatchFromPlt(
    lev, geom, vars["x_velocity"], UMX, AMREX_SPACEDIM, S_new);
  pltData.fillPatchFromPlt(lev, geom, vars["Temp"], UTEMP, 1, S_new);

  // Copy species from the plot file
  for (int n = 0; n < spec_names.size(); n++) {
    const auto& spec = spec_names.at(n);
    const int pos = find_position(plt_vars, "Y(" + spec + ")");
    if (pos != -1) {
      pltData.fillPatchFromPlt(lev, geom, pos, UFS + n, 1, S_new);
    }
  }

  // Sanity check the species, clean them up if they aren't too bad
  auto sarrs = S_new.arrays();
  const auto tol = init_pltfile_massfrac_tol;
  amrex::ParallelFor(
    S_new, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
      auto sarr = sarrs[nbx];
      amrex::Real sumY = 0.0;

      for (int n = 0; n < NUM_SPECIES; n++) {
        // if the species is not too far out of bounds, clip it
        const auto mf = sarr(i, j, k, UFS + n);
        if ((mf < 0.0) || (1.0 < mf)) {
          if (((-tol < mf) && (mf < 0.0)) || ((1.0 < mf) && (mf < 1 + tol))) {
            sarr(i, j, k, UFS + n) =
              amrex::min<amrex::Real>(1.0, amrex::max<amrex::Real>(0.0, mf));
          } else {
#ifdef AMREX_USE_GPU
            AMREX_DEVICE_PRINTF(
              "Species mass fraction is out of bounds (spec, value): (%d, %g)",
              n, mf);
            amrex::Abort();
#else
            amrex::Abort(
              "Species mass fraction is out of bounds (spec, value): ( " +
              std::to_string(n) + ", " + std::to_string(mf) + ")");
#endif
          }
        }

        sumY += sarr(i, j, k, UFS + n);
      }

      // If the sumY isn't too far from 1, renormalize
      if (std::abs(1.0 - sumY) < tol) {
        for (int n = 0; n < NUM_SPECIES; n++) {
          sarr(i, j, k, UFS + n) /= sumY;
        }
      } else {
#ifdef AMREX_USE_GPU
        AMREX_DEVICE_PRINTF(
          "Species mass fraction don't sum to 1. The sum is: %g", sumY);
        amrex::Abort();
#else
        amrex::Abort(
          "Species mass fraction don't sum to 1. The sum is: " +
          std::to_string(sumY));
#endif
      }
    });
  amrex::Gpu::synchronize();

  // Convert to conserved variables
  amrex::ParallelFor(
    S_new, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
      auto sarr = sarrs[nbx];
      const amrex::Real rho = sarr(i, j, k, URHO);
      const amrex::Real temp = sarr(i, j, k, UTEMP);
      amrex::Real massfrac[NUM_SPECIES] = {0.0};
      for (int n = 0; n < NUM_SPECIES; n++) {
        massfrac[n] = sarr(i, j, k, UFS + n);
        sarr(i, j, k, UFS + n) *= rho;
      }
      AMREX_D_TERM(sarr(i, j, k, UMX) *= rho;, sarr(i, j, k, UMY) *= rho;
                   , sarr(i, j, k, UMZ) *= rho;)

      auto eos = pele::physics::PhysicsType::eos();
      amrex::Real eint = 0.0;
      eos.RTY2E(rho, temp, massfrac, eint);

      sarr(i, j, k, UEINT) = rho * eint;
      sarr(i, j, k, UEDEN) =
        rho * eint + 0.5 *
                       (AMREX_D_TERM(
                         sarr(i, j, k, UMX) * sarr(i, j, k, UMX),
                         +sarr(i, j, k, UMY) * sarr(i, j, k, UMY),
                         +sarr(i, j, k, UMZ) * sarr(i, j, k, UMZ))) /
                       rho;
    });
  amrex::Gpu::synchronize();
}
