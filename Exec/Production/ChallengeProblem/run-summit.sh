#!/bin/bash -l

#BSUB -J pelec-challenge-summit
#BSUB -o pelec-challenge-summit.o%J
#BSUB -P CMB138
#BSUB -W 60
#BSUB -nnodes 128
#BSUB -N

cmd() {
  echo "+ $@"
  eval "$@"
}

cmd "module unload xl"
cmd "module load gcc/9.3.0"

cmd "export SPACK_MANAGER=${PROJWORK}/cmb138/jrood/spack-manager-${LMOD_SYSTEM_NAME}"
cmd "source ${SPACK_MANAGER}/start.sh && spack-start"
cmd "spack env activate -d ${SPACK_MANAGER}/environments/exawind-${LMOD_SYSTEM_NAME}"
cmd "spack load exawind+amr_wind_gpu~nalu_wind_gpu"

cmd "jsrun -n 768 -r 6 -a 1 -c 1 -g 1 ./PeleC3d.gnu.TPROF.MPI.CUDA.ex challenge.inp max_step=500"
