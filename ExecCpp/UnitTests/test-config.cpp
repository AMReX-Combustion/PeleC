/** \file test_config.cpp
 *
 *  Tests various configurations for GPU builds
 */

#include "gtest/gtest.h"
#include "AMReX_ccse-mpi.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_Print.H"
#include "AMReX_Gpu.H"

namespace amrex {
const char* buildInfoGetGitHash(int i);
}

namespace pelec_tests {

TEST(Configuration, Build)
{
  const char* pc_git = amrex::buildInfoGetGitHash(1);
  const char* amrex_git = amrex::buildInfoGetGitHash(2);
  amrex::Print() << "PeleC SHA = " << pc_git << std::endl
                 << "AMReX SHA = " << amrex_git << std::endl;
}

TEST(Configuration, MPI)
{
#ifdef AMREX_USE_MPI
  int nprocs = amrex::ParallelDescriptor::NProcs();

  amrex::Print() << "MPI configuration: " << nprocs << " ranks" << std::endl;
  char mpi_lib_ver[MPI_MAX_LIBRARY_VERSION_STRING];
  int len;
  MPI_Get_library_version(mpi_lib_ver, &len);
  amrex::Print() << mpi_lib_ver << std::endl;
#else
  amrex::Print() << "PeleC not built with MPI support." << std::endl;
  GTEST_SKIP();
#endif
}

TEST(Configuration, CUDA)
{
#ifdef AMREX_USE_CUDA
#if defined(CUDA_VERSION)
  amrex::Print() << "CUDA configuration: "
                 << "CUDA_VERSION: " << CUDA_VERSION << " "
                 << CUDA_VERSION / 1000 << "." << (CUDA_VERSION % 1000) / 10
                 << std::endl;
#endif
  cudaError_t error;
  int ndevices;
  cudaDeviceProp dev;

  error = cudaGetDeviceCount(&ndevices);
  if (error != cudaSuccess) {
    ADD_FAILURE() << cudaGetErrorString(error);
    return;
  } else {
    std::cout << "Num. devices = " << ndevices << std::endl;
  }

  const int myrank = amrex::ParallelDescriptor::MyProc();
  const int rankDevice = amrex::Gpu::Device::deviceId();
  error = cudaGetDeviceProperties(&dev, rankDevice);
  if (error != cudaSuccess) {
    ADD_FAILURE() << cudaGetErrorString(error);
    return;
  }
  char busid[512];
  cudaDeviceGetPCIBusId(busid, 512, rankDevice);
  std::cout << "[" << myrank << "] " << rankDevice << ": " << dev.name
            << " CC: " << dev.major << "." << dev.minor << " ID: " << busid
            << " GM: " << (static_cast<double>(dev.totalGlobalMem) / (1 << 30))
            << "GB"
            << " ShMem/Blk: " << (dev.sharedMemPerBlock / (1 << 10)) << "KB"
            << std::endl;
#else
  amrex::Print() << "PeleC not built with CUDA support" << std::endl;
  GTEST_SKIP();
#endif
}

} // namespace pelec_tests
