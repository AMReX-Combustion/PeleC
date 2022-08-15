#!/bin/bash -l

cmd() {
  echo "+ $@"
  eval "$@"
}

# Not the best assumption, but change this to wherever spack-manager is
cmd "export SPACK_MANAGER=${HOME}/exawind/spack-manager"
cmd "source ${SPACK_MANAGER}/start.sh && spack-start"
cmd "spack load mpich"
cmd "mpirun -np 8 ./PeleC3d.llvm.TPROF.MPI.ex challenge.inp amr.max_level=1"
