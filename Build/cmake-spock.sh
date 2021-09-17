#!/bin/bash -l

module unload PrgEnv-cray
module load PrgEnv-amd
module unload rocm rocm-compiler
module load rocm/4.3.0
module load cmake
module list

cmake -DCMAKE_INSTALL_PREFIX:PATH=./install \
      -DCMAKE_CXX_COMPILER:STRING=$(which hipcc) \
      -DCMAKE_C_COMPILER:STRING=$(which hipcc) \
      -DCMAKE_BUILD_TYPE:STRING=Release \
      -DPELEC_DIM:STRING=3 \
      -DPELEC_ENABLE_AMREX_EB:BOOL=ON \
      -DPELEC_ENABLE_MPI:BOOL=OFF \
      -DPELEC_ENABLE_TESTS:BOOL=ON \
      -DPELEC_ENABLE_FCOMPARE:BOOL=OFF \
      -DPELEC_ENABLE_FCOMPARE_FOR_TESTS:BOOL=OFF \
      -DPELEC_ENABLE_MASA:BOOL=OFF \
      -DPELEC_ENABLE_ALL_WARNINGS:BOOL=ON \
      -DPELEC_ENABLE_HIP:BOOL=ON \
      -DAMReX_AMD_ARCH=gfx908 \
      -DAMReX_GPU_RDC=OFF \
      -DCMAKE_CXX_COMPILER_ID="Clang" \
      -DCMAKE_CXX_COMPILER_VERSION=12.0 \
      -DCMAKE_CXX_STANDARD=17 \
      -DCMAKE_CXX_STANDARD_COMPUTED_DEFAULT=17 \
      -DPYTHON_EXECUTABLE=$(which python3) \
      -DPELEC_PRECISION:STRING=DOUBLE \
      .. 
nice cmake --build . --parallel $(nproc)
#ctest -j $(nproc)
