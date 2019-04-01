#!/bin/bash

# Example CMake config script for and OSX laptop with OpenMPI
#export PATH=$PATH:$(spack location -i ninja)/bin # or whatever to git the kitware version of ninja with fortran in your path
cmake -G Ninja -DCMAKE_INSTALL_PREFIX:PATH=./install \
               -DCMAKE_BUILD_TYPE:STRING=Release \
               -DPELEC_ENABLE_MPI:BOOL=ON \
               -DCMAKE_CXX_COMPILER:STRING=mpicxx \
               -DCMAKE_C_COMPILER:STRING=mpicc \
               -DCMAKE_Fortran_COMPILER:STRING=mpifort \
               -DMPI_CXX_COMPILER:STRING=mpicxx \
               -DMPI_C_COMPILER:STRING=mpicc \
               -DMPI_Fortran_COMPILER:STRING=mpifort \
               .. && ninja -j8

# Extra options
#      -DPELEC_ENABLE_MASA:BOOL=ON \
#      -DMASA_DIR:STRING=$(spack location -i masa) \
#      -DMPIEXEC_PREFLAGS:STRING=--oversubscribe \
#      -DENABLE_TESTS:BOOL=ON \
#      -DENABLE_DOCUMENTATION:BOOL=ON \
#      -DENABLE_FCOMPARE:BOOL=ON \
