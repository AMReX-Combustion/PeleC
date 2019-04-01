#!/bin/bash

# Example CMake config script for and OSX laptop with OpenMPI
cmake -DCMAKE_INSTALL_PREFIX:PATH=./install \
      -DCMAKE_BUILD_TYPE:STRING=Release \
      -DENABLE_TESTS:BOOL=OFF \
      -DENABLE_FCOMPARE:BOOL=OFF \
      -DENABLE_DOCUMENTATION:BOOL=OFF \
      -DPELEC_ENABLE_MPI:BOOL=ON \
      -DCMAKE_CXX_COMPILER:STRING=mpicxx \
      -DCMAKE_C_COMPILER:STRING=mpicc \
      -DCMAKE_Fortran_COMPILER:STRING=mpifort \
      -DMPI_CXX_COMPILER:STRING=mpicxx \
      -DMPI_C_COMPILER:STRING=mpicc \
      -DMPI_Fortran_COMPILER:STRING=mpifort \
      -DMPIEXEC_PREFLAGS:STRING=--oversubscribe \
      ..
