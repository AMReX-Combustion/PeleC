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
cmd "module load gcc/9.3.0"
cmd "module load intel-parallel-studio/cluster.2020.2"
cmd "module load binutils"
cmd "module load python"
cmd "module load mpt"
cmd "module load cmake"
cmd "mpirun -np 2304 ./PeleC3d.intel.TPROF.MPI.ex challenge.inp"
