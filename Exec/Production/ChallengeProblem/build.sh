#!/bin/bash -l

cmd() {
  echo "+ $@"
  eval "$@"
}

OS=$(uname -s)

if [ "${LMOD_SYSTEM_NAME}" == 'summit' ]; then
  cmd "module unload xl"
  cmd "module load gcc/10.2.0"
  cmd "module load cuda/11.4.2"
  cmd "module load cmake"
  cmd "module load binutils"
  cmd "module load netlib-lapack/3.9.1"
  cmd "export SPACK_MANAGER=${WORLDWORK}/cmb138/software/spack-manager-${LMOD_SYSTEM_NAME}"
  ARGS="USE_CUDA=TRUE PELE_USE_MAGMA=TRUE"
elif [ "${LMOD_SYSTEM_NAME}" == 'crusher' ]; then
  cmd "module unload PrgEnv-cray"
  cmd "module load PrgEnv-amd"
  cmd "module load rocm/5.1.0"
  cmd "module load cmake"
  cmd "export SPACK_MANAGER=${WORLDWORK}/cmb138/software/spack-manager-${LMOD_SYSTEM_NAME}"
  ARGS="USE_HIP=TRUE PELE_USE_MAGMA=TRUE"
elif [ "${NREL_CLUSTER}" == 'eagle' ]; then
  cmd "source /nopt/nrel/ecom/hpacf/env.sh"
  cmd "module load gcc/9.3.0"
  cmd "module load intel-parallel-studio/cluster.2020.2"
  cmd "module load binutils"
  cmd "module load python"
  cmd "module load mpt"
  cmd "module load cmake"
  cmd "export SPACK_MANAGER=/scratch/jrood/spack-manager-${NREL_CLUSTER}"
  ARGS="COMP=intel"
elif [ "${OS}" == 'Darwin' ]; then
  # Not the best assumption, but change this to wherever spack-manager is
  cmd "export SPACK_MANAGER=${HOME}/exawind/spack-manager"
  cmd "source ${SPACK_MANAGER}/start.sh && spack-start"
  cmd "spack load cmake"
  cmd "spack load mpich"
  ARGS="COMP=llvm"
fi

cmd "source ${SPACK_MANAGER}/start.sh && spack-start"
HDF5_ROOT=$(spack location -i hdf5)
H5Z_ROOT=$(spack location -i h5z-zfp)
ZFP_ROOT=$(spack location -i zfp)
ASCENT_ROOT=$(spack location -i ascent)
CONDUIT_ROOT=$(spack location -i conduit)

ARGS="USE_HDF5=TRUE USE_HDF5_ZFP=TRUE USE_ASCENT=TRUE USE_CONDUIT=TRUE HDF5_HOME=${HDF5_ROOT} H5Z_HOME=${H5Z_ROOT} ZFP_HOME=${ZFP_ROOT} ASCENT_DIR=${ASCENT_ROOT} CONDUIT_DIR=${CONDUIT_ROOT} ${ARGS}"

cmd "make ${ARGS} TPLrealclean"
cmd "make ${ARGS} realclean"
cmd "make ${ARGS} TPL"
cmd "make ${ARGS} -j16"
