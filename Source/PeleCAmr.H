#ifndef PELECAMR_H
#define PELECAMR_H

#include <AMReX_ParmParse.H>
#include <AMReX_EB2.H>
#include <AMReX_PlotFileUtil.H>
#ifdef PELEC_USE_HDF5
#include <AMReX_PlotFileUtilHDF5.H>
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
    const bool write_hdf5_plots);
};

#endif
