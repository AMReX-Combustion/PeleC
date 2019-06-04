#!/bin/bash

# Example CMake config script for and OSX laptop with OpenMPI

cmake -DCMAKE_INSTALL_PREFIX:PATH=./install \
      -DCMAKE_BUILD_TYPE:STRING=Release \
      -DPELEC_ENABLE_MPI:BOOL=ON \
      -DCMAKE_CXX_COMPILER:STRING=mpicxx \
      -DCMAKE_C_COMPILER:STRING=mpicc \
      -DCMAKE_Fortran_COMPILER:STRING=mpifort \
      ..

# Extra options
      #-DENABLE_FCOMPARE:BOOL=ON \
      #-DENABLE_FEXTREMA:BOOL=ON \
      #-DENABLE_TESTS:BOOL=ON \
      #-DENABLE_VERIFICATION:BOOL=ON \
      #-DTEST_WITH_FCOMPARE:BOOL=OFF \
      #-DTEST_WITH_FEXTREMA:BOOL=OFF \
      #-DPELEC_ENABLE_MASA:BOOL=ON \
      #-DMASA_DIR:STRING=$(spack location -i masa) \
      #-DMPIEXEC_PREFLAGS:STRING=--oversubscribe \
      #-DENABLE_DOCUMENTATION:BOOL=ON \
