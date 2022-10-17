#!/bin/bash -l

cmd() {
  echo "+ $@"
  eval "$@"
}

# Not the best assumption, but change this to wherever spack-manager is
cmd "export SPACK_MANAGER=${HOME}/exawind/spack-manager"
cmd "source ${SPACK_MANAGER}/start.sh && spack-start"
cmd "spack env activate -d ${SPACK_MANAGER}/environments/exawind-darwin"
cmd "spack load mpich"
cmd "mpirun -np 8 ./PeleC3d.llvm.TPROF.MPI.ex challenge.inp amr.n_cell=128 128 32 amr.max_level=2 amr.checkpoint_files_output=0 cvode.solve_type=GMRES max_step=10"
