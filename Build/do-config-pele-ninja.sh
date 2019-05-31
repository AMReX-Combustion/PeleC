#!/bin/bash

# Example CMake config script for and OSX laptop with OpenMPI

#export PATH=$PATH:$(spack location -i ninja)/bin # or whatever to get the kitware version of ninja with fortran in your path

cmake -G Ninja -DCMAKE_INSTALL_PREFIX:PATH=./install \
               -DCMAKE_BUILD_TYPE:STRING=Release \
               -DPELEC_ENABLE_MPI:BOOL=ON \
               -DCMAKE_CXX_COMPILER:STRING=mpicxx \
               -DCMAKE_C_COMPILER:STRING=mpicc \
               -DCMAKE_Fortran_COMPILER:STRING=mpifort \
               .. && ninja -j8

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
      #-DENABLE_DOCUMENTATION:BOOL=OFF \
