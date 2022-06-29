#!/bin/bash -l

cmd() {
  echo "+ $@"
  eval "$@"
}

if [ "${LMOD_SYSTEM_NAME}" == 'summit' ]; then
  cmd "module unload xl"
  cmd "module load gcc/9.3.0"
  cmd "export SPACK_MANAGER=${PROJWORK}/cmb138/jrood/spack-manager-${LMOD_SYSTEM_NAME}"
  cmd "source ${SPACK_MANAGER}/start.sh && spack-start"
  cmd "spack env activate -d ${SPACK_MANAGER}/environments/exawind-${LMOD_SYSTEM_NAME}"
  cmd "spack load exawind+amr_wind_gpu~nalu_wind_gpu"
  ARGS="USE_CUDA=TRUE"
elif [ "${LMOD_SYSTEM_NAME}" == 'crusher' ]; then
  cmd "module unload PrgEnv-cray"
  cmd "module load PrgEnv-amd"
  cmd "export SPACK_MANAGER=${PROJWORK}/cmb138/jrood/spack-manager-${LMOD_SYSTEM_NAME}"
  cmd "source ${SPACK_MANAGER}/start.sh && spack-start"
  cmd "spack env activate -d ${SPACK_MANAGER}/environments/exawind-${LMOD_SYSTEM_NAME}"
  cmd "spack load exawind+amr_wind_gpu~nalu_wind_gpu"
  ARGS="USE_HIP=TRUE"
elif [ "${NREL_CLUSTER}" == 'eagle' ]; then
  cmd "source /nopt/nrel/ecom/hpacf/env.sh"
  cmd "module load intel-parallel-studio/cluster.2020.2"
  cmd "module load binutils"
  cmd "module load python"
  cmd "module load mpt"
  cmd "export SPACK_MANAGER=/scratch/jrood/spack-manager-${NREL_CLUSTER}"
  cmd "source ${SPACK_MANAGER}/start.sh && spack-start"
  cmd "spack env activate -d ${SPACK_MANAGER}/environments/exawind-${NREL_CLUSTER}"
  cmd "spack load exawind~amr_wind_gpu~nalu_wind_gpu"
  ARGS="COMP=intel"
fi

HDF5_HOME=$(spack location -i hdf5)
H5Z_HOME=$(spack location -i h5z-zfp)
ZFP_HOME=$(spack location -i zfp)

ARGS="USE_MPI=TRUE TINY_PROFILE=TRUE USE_HDF5=TRUE USE_HDF5_ZFP=TRUE HDF5_HOME=${HDF5_HOME} H5Z_HOME=${H5Z_HOME} ZFP_HOME=${ZFP_HOME} ${ARGS}"

cmd "make ${ARGS} TPLrealclean"
cmd "make ${ARGS} realclean"
cmd "make ${ARGS} TPL"
cmd "make ${ARGS} -j16"
