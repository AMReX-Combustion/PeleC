#!/bin/bash -l

cmd() {
  echo "+ $@"
  eval "$@"
}

if [ "${LMOD_SYSTEM_NAME}" == 'summit' ]; then
  cmd "module unload xl"
  cmd "module load gcc/9.3.0"
  cmd "module load python/3.7.7"
  cmd "module load cuda/11.4.2"
  cmd "module load hdf5/1.10.7"
  cmd "module load cmake"
  BASE_DIR="${PROJWORK}/cfd116/jrood/spack-manager-summit/spack/opt/spack/linux-rhel8-ppc64le/gcc-9.3.0"
  H5Z_HOME="${BASE_DIR}/h5z-zfp-develop-czinnk7mbsrkmf7o74mxzpzv7hqp34sl"
  ZFP_HOME="${BASE_DIR}/zfp-0.5.5-dsefgrt3iinobwohmbuspnpar3yasm4l"
elif [ "${LMOD_SYSTEM_NAME}" == 'crusher' ]; then
  cmd "module unload PrgEnv-cray"
  cmd "module load PrgEnv-amd"
  cmd "module load rocm/5.1.0"
  cmd "module load hdf5/1.12.1"
  cmd "module load cmake"
  BASE_DIR="${PROJWORK}/cfd116/jrood/spack-manager-crusher/spack/opt/spack/cray-sles15-zen3/clang-14.0.0"
  H5Z_HOME="${BASE_DIR}/h5z-zfp-develop-y7qstijn2lvrhdvhswgxfnzcyh2nctfg"
  ZFP_HOME="${BASE_DIR}/zfp-0.5.5-dxggbklt74oeev6xtnici2mjq4zh2oqx"
fi

ARGS="USE_HDF5=TRUE USE_HDF5_ZFP=TRUE HDF5_HOME=${OLCF_HDF5_ROOT} H5Z_HOME=${H5Z_HOME} ZFP_HOME=${ZFP_HOME}"

cmd "make ${ARGS} TPLrealclean"
cmd "make ${ARGS} realclean"
cmd "make ${ARGS} TPL"
cmd "make ${ARGS} -j16"
