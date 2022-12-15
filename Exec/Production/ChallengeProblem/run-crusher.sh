#!/bin/bash -l

#SBATCH -J pelec-challenge-crusher
#SBATCH -o %x.o%j
#SBATCH -A CMB138_crusher
#SBATCH -t 0:30:00
#SBATCH -N 32

set -e

cmd() {
  echo "+ $@"
  eval "$@"
}

cmd "module load rocm/5.4.0"
cmd "srun -N32 -n256 -c1 --gpus-per-node=8 --gpu-bind=closest ./PeleC3d.hip.x86-trento.TPROF.MPI.HIP.ex challenge.inp max_step=53210"
