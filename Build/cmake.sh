#!/bin/bash

# Example CMake config script for running test suite on a MacOS laptop with OpenMPI:

cmake -DCMAKE_INSTALL_PREFIX:PATH=./install \
      -DCMAKE_CXX_COMPILER:STRING=mpicxx \
      -DCMAKE_C_COMPILER:STRING=mpicc \
      -DMPIEXEC_PREFLAGS:STRING=--oversubscribe \
      -DCMAKE_BUILD_TYPE:STRING=Release \
      -DPELEC_DIM:STRING=3 \
      -DPELEC_ENABLE_AMREX_EB:BOOL=ON \
      -DPELEC_ENABLE_MPI:BOOL=ON \
      -DPELEC_ENABLE_TESTS:BOOL=ON \
      -DPELEC_ENABLE_FCOMPARE:BOOL=OFF \
      -DPELEC_ENABLE_FCOMPARE_FOR_TESTS:BOOL=OFF \
      -DPELEC_ENABLE_SUNDIALS:BOOL=OFF \
      -DPELEC_ENABLE_MASA:BOOL=ON \
      -DMASA_DIR:STRING=$(spack location -i masa) \
      -DPELEC_ENABLE_ALL_WARNINGS:BOOL=ON \
      -DPELEC_ENABLE_CPPCHECK:BOOL=OFF \
      -DPELEC_ENABLE_CLANG_TIDY:BOOL=OFF \
      -DPELEC_ENABLE_CUDA:BOOL=OFF \
      -DAMReX_CUDA_ARCH=Volta \
      -DPYTHON_EXECUTABLE=$(which python3) \
      -DPELEC_PRECISION:STRING=DOUBLE \
      .. 
#cmake --build . --parallel $(sysctl -n hw.ncpu)
#ctest -j $(sysctl -n hw.ncpu)
