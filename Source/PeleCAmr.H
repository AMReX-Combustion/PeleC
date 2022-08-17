#ifndef PELECAMR_H
#define PELECAMR_H

#include <AMReX_ParmParse.H>
#include <AMReX_EB2.H>
#include <AMReX_PlotFileUtil.H>
#ifdef AMREX_USE_HDF5
#include <AMReX_PlotFileUtilHDF5.H>
#endif
#ifdef AMREX_USE_ASCENT
#include <AMReX_Conduit_Blueprint.H>
#include <ascent.hpp>
#include "PeleAscent.H"
#endif

#include "PeleC.H"

class PeleCAmr : public amrex::Amr
{
  using amrex::Amr::Amr;

public:
  void writePlotFile() override;
  void writeSmallPlotFile() override;
  void writePlotFileDoit(
    const std::string& pltfile,
    const bool regular,
    // Optional arguments
    const bool write_hdf5_plots = false,
    const std::string& hdf5_compression = "None@0");
#ifdef AMREX_USE_ASCENT
  void doInSituViz(const int step);
  pele::PeleAscent pele_ascent;
#endif
private:
  void constructPlotMF(
    const bool regular,
    amrex::Vector<std::unique_ptr<amrex::MultiFab>>& plotMFs,
    amrex::Vector<std::string>& plt_var_names);
};

#endif
