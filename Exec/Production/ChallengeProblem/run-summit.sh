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
cmd "module load gcc/10.2.0"
cmd "module load cuda/11.4.2"
cmd "module load cmake"
cmd "module load binutils"
cmd "module load netlib-lapack/3.9.1"
cmd "jsrun -n 768 -r 6 -a 1 -c 1 -g 1 ./PeleC3d.gnu.TPROF.MPI.CUDA.ex challenge.inp max_step=500"
