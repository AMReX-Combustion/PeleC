#!/bin/bash -l

cmd() {
  echo "+ $@"
  eval "$@"
}

OS=$(uname -s)

if [ "${LMOD_SYSTEM_NAME}" == 'summit' ]; then
  GNUMAKE_ARGS="USE_CUDA=TRUE PELE_USE_MAGMA=TRUE"
  SPACK_COMPILER="%gcc"
  SPACK_VARIANTS=""
  MACHINE="summit"
  cmd "module unload xl"
  cmd "module load gcc/10.2.0"
  cmd "module load cuda/11.4.2"
  cmd "module load cmake"
  cmd "module load binutils"
  cmd "module load netlib-lapack/3.9.1"
  cmd "export SPACK_MANAGER=${WORLDWORK}/cmb138/software/spack-manager-${LMOD_SYSTEM_NAME}"
  cmd "source ${SPACK_MANAGER}/start.sh && spack-start"
  cmd "spack env activate -d ${SPACK_MANAGER}/environments/exawind-${MACHINE}"
elif [ "${LMOD_SYSTEM_NAME}" == 'crusher' ]; then
  GNUMAKE_ARGS="USE_HIP=TRUE PELE_USE_MAGMA=TRUE"
  SPACK_COMPILER="%clang"
  SPACK_VARIANTS=""
  MACHINE="crusher"
  cmd "module load rocm/5.4.0"
  cmd "module load cray-python"
  cmd "module unload cray-libsci"
  cmd "module load cray-libsci/22.11.1.2"
  cmd "module load cmake"
  cmd "export SPACK_MANAGER=${WORLDWORK}/cmb138/software/spack-manager-${LMOD_SYSTEM_NAME}"
  cmd "source ${SPACK_MANAGER}/start.sh && spack-start"
  cmd "spack env activate -d ${SPACK_MANAGER}/environments/exawind-${MACHINE}"
elif [ "${NREL_CLUSTER}" == 'eagle' ]; then
  GNUMAKE_ARGS="COMP=intel"
  SPACK_COMPILER="%intel"
  SPACK_VARIANTS="~cuda"
  MACHINE="eagle"
  cmd "source /nopt/nrel/ecom/hpacf/env.sh"
  cmd "module load gcc/9.3.0"
  cmd "module load intel-parallel-studio/cluster.2020.2"
  cmd "module load binutils"
  cmd "module load python"
  cmd "module load mpt"
  cmd "module load cmake"
  cmd "module load git"
  cmd "export SPACK_MANAGER=/nopt/nrel/ecom/hpacf/spack-manager/2022-10-14/spack-manager"
  cmd "source ${SPACK_MANAGER}/start.sh && spack-start"
  cmd "spack env activate -d ${SPACK_MANAGER}/environments/exawind-${MACHINE}"
elif [ "${OS}" == 'Darwin' ]; then
  GNUMAKE_ARGS="COMP=llvm"
  SPACK_COMPILER="%apple-clang"
  SPACK_VARIANTS=""
  MACHINE="darwin"
  # Not the best assumption, but change this to wherever spack-manager is
  cmd "export SPACK_MANAGER=${HOME}/exawind/spack-manager"
  cmd "source ${SPACK_MANAGER}/start.sh && spack-start"
  cmd "spack env activate -d ${SPACK_MANAGER}/environments/exawind-${MACHINE}"
  cmd "spack load cmake"
  cmd "spack load mpich"
fi

#HDF5_ROOT=$(spack location -i hdf5 ${SPACK_COMPILER})
#H5Z_ROOT=$(spack location -i h5z-zfp ${SPACK_COMPILER})
#ZFP_ROOT=$(spack location -i zfp ${SPACK_COMPILER})
CONDUIT_ROOT=$(spack location -i conduit ${SPACK_COMPILER})
ASCENT_ROOT=$(spack location -i ascent ${SPACK_VARIANTS} ${SPACK_COMPILER})

#GNUMAKE_ARGS="USE_HDF5=FALSE USE_HDF5_ZFP=FALSE USE_ASCENT=TRUE USE_CONDUIT=TRUE HDF5_HOME=${HDF5_ROOT} H5Z_HOME=${H5Z_ROOT} ZFP_HOME=${ZFP_ROOT} ASCENT_DIR=${ASCENT_ROOT} CONDUIT_DIR=${CONDUIT_ROOT} ${GNUMAKE_ARGS}"
GNUMAKE_ARGS="USE_ASCENT=TRUE USE_CONDUIT=TRUE ASCENT_DIR=${ASCENT_ROOT} CONDUIT_DIR=${CONDUIT_ROOT} ${GNUMAKE_ARGS}"

cmd "make ${GNUMAKE_ARGS} TPLrealclean"
cmd "make ${GNUMAKE_ARGS} realclean"
cmd "make ${GNUMAKE_ARGS} TPL"
if [ "${LMOD_SYSTEM_NAME}" == 'summit' ]; then
  cmd "make ${GNUMAKE_ARGS} -j4"
else
  cmd "make ${GNUMAKE_ARGS} -j16"
fi
