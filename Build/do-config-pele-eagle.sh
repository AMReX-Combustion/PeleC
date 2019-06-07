#!/bin/bash -l

# Example CMake config script for PeleC on Eagle

MODULES=modules
COMPILER=gcc-7.4.0

cmd() {
  echo "+ $@"
  eval "$@"
}

set -e

# Load environment
cmd "module purge"
cmd "module unuse ${MODULEPATH}"
cmd "module use /nopt/nrel/ecom/hpacf/binaries/${MODULES}"
cmd "module use /nopt/nrel/ecom/hpacf/compilers/${MODULES}"
cmd "module use /nopt/nrel/ecom/hpacf/utilities/${MODULES}"
cmd "module use /nopt/nrel/ecom/hpacf/software/${MODULES}/${COMPILER}"
cmd "module load gcc"
cmd "module load git"
cmd "module load python/3.7.3"
cmd "module load mpich"
cmd "module load masa"
cmd "module load cmake"

# Clean before cmake configure
set +e
cmd "rm -rf CMakeFiles"
cmd "rm -f CMakeCache.txt"
set -e

(set -x; cmake -DCMAKE_INSTALL_PREFIX:PATH=./install \
               -DCMAKE_BUILD_TYPE:STRING=Release \
               -DPELEC_ENABLE_MPI:BOOL=ON \
               -DCMAKE_CXX_COMPILER:STRING=mpicxx \
               -DCMAKE_C_COMPILER:STRING=mpicc \
               -DCMAKE_Fortran_COMPILER:STRING=mpifort \
               -DENABLE_TESTS:BOOL=ON \
               -DENABLE_VERIFICATION:BOOL=ON \
               -DTEST_WITH_FCOMPARE:BOOL=OFF \
               -DTEST_WITH_FEXTREMA:BOOL=ON \
               -DPELEC_ENABLE_MASA:BOOL=ON \
               -DMASA_DIR:STRING=${MASA_ROOT_DIR} \
               .. && nice make -j24)

