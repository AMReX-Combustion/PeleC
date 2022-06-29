#!/bin/bash -l

#SBATCH -J pelec-challenge-eagle
#SBATCH -o %x.o%j
#SBATCH -A hpacf
#SBATCH -t 2:00:00
#SBATCH -N 64

set -e

cmd() {
  echo "+ $@"
  eval "$@"
}

cmd "source /nopt/nrel/ecom/hpacf/env.sh"
cmd "module load binutils"
cmd "module load python"
cmd "module load mpt"

cmd "export SPACK_MANAGER=/scratch/jrood/spack-manager-${NREL_CLUSTER}"
cmd "source ${SPACK_MANAGER}/start.sh && spack-start"
cmd "spack env activate -d ${SPACK_MANAGER}/environments/exawind-${NREL_CLUSTER}"
cmd "spack load exawind~amr_wind_gpu~nalu_wind_gpu"

cmd "mpirun -np 2304 ./PeleC3d.hip.x86-trento.TPROF.MPI.HIP.ex challenge.inp"
