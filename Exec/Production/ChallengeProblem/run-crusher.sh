#!/bin/bash -l

#SBATCH -J pelec-challenge-crusher
#SBATCH -o %x.o%j
#SBATCH -A CMB138_crusher
#SBATCH -t 2:00:00
#SBATCH -N 32

set -e

cmd() {
  echo "+ $@"
  eval "$@"
}

cmd "module load rocm/5.1.0"
cmd "module load hdf5/1.12.1"
cmd "module load cmake"
cmd "export SPACK_MANAGER=${PROJWORK}/cfd116/jrood/spack-manager-${LMOD_SYSTEM_NAME}"
cmd "source ${SPACK_MANAGER}/start.sh && spack-start"
cmd "spack env activate -d ${SPACK_MANAGER}/environments/exawind-${LMOD_SYSTEM_NAME}"
cmd "spack load exawind+amr_wind_gpu~nalu_wind_gpu"
cmd "srun -N32 -n256 -c1 --gpus-per-node=8 --gpu-bind=closest ./PeleC3d.hip.x86-trento.TPROF.MPI.HIP.ex challenge.inp"
