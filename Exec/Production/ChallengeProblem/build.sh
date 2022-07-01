#!/bin/bash -l

cmd() {
  echo "+ $@"
  eval "$@"
}

if [ "${LMOD_SYSTEM_NAME}" == 'summit' ]; then
  cmd "module unload xl"
  cmd "module load gcc/9.3.0"
  cmd "module load cuda/11.4.2"
  cmd "module load cmake"
  cmd "module load binutils"
  cmd "export SPACK_MANAGER=${PROJWORK}/cmb138/jrood/spack-manager-${LMOD_SYSTEM_NAME}"
  ARGS="USE_CUDA=TRUE"
elif [ "${LMOD_SYSTEM_NAME}" == 'crusher' ]; then
  cmd "module unload PrgEnv-cray"
  cmd "module load PrgEnv-amd"
  cmd "module load rocm/5.1.0"
  cmd "module load cmake"
  cmd "export SPACK_MANAGER=${PROJWORK}/cmb138/jrood/spack-manager-${LMOD_SYSTEM_NAME}"
  ARGS="USE_HIP=TRUE"
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
fi

cmd "source ${SPACK_MANAGER}/start.sh && spack-start"
HDF5_HOME=$(spack location -i hdf5)
H5Z_HOME=$(spack location -i h5z-zfp)
ZFP_HOME=$(spack location -i zfp)
ASCENT_DIR=$(spack location -i ascent)
CONDUIT_DIR=$(spack location -i conduit)

ARGS="USE_MPI=TRUE TINY_PROFILE=TRUE USE_HDF5=TRUE USE_HDF5_ZFP=TRUE USE_ASCENT=TRUE USE_CONDUIT=TRUE HDF5_HOME=${HDF5_HOME} H5Z_HOME=${H5Z_HOME} ZFP_HOME=${ZFP_HOME} ASCENT_DIR=${ASCENT_DIR} CONDUIT_DIR=${CONDUIT_DIR} ${ARGS}"

cmd "make ${ARGS} TPLrealclean"
cmd "make ${ARGS} realclean"
cmd "make ${ARGS} TPL"
cmd "make ${ARGS} -j16"
