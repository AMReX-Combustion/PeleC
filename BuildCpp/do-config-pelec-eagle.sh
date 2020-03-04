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
               -DPELEC_ENABLE_TESTS:BOOL=ON \
               -DPELEC_ENABLE_VERIFICATION:BOOL=ON \
               -DPELEC_TEST_WITH_FCOMPARE:BOOL=OFF \
               -DPELEC_TEST_WITH_FEXTREMA:BOOL=ON \
               -DPELEC_ENABLE_MASA:BOOL=ON \
               -DMASA_DIR:STRING=${MASA_ROOT_DIR} \
               -DPELEC_ENABLE_EB:BOOL=OFF \
               -DPELEC_ENABLE_MASA:BOOL=OFF \
               -DPELEC_ENABLE_REACTIONS:BOOL=OFF \
               -DPELEC_ENABLE_EXPLICIT_REACT:BOOL=ON \
               -DPELEC_EOS_MODEL:STRING=GammaLaw \
               -DPELEC_REACTIONS_MODEL:STRING=Null \
               -DPELEC_CHEMISTRY_MODEL:STRING=Null \
               -DPELEC_TRANSPORT_MODEL:STRING=Constant \
               -DPELEC_ENABLE_PARTICLES:BOOL=OFF \
               .. && nice make -j24)

