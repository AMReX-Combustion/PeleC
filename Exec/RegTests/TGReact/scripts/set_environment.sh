module purge
MODULES=modules
COMPILER=gcc-7.4.0
module unuse ${MODULEPATH}
module use /nopt/nrel/ecom/hpacf/binaries/${MODULES}
module use /nopt/nrel/ecom/hpacf/compilers/${MODULES}
module use /nopt/nrel/ecom/hpacf/utilities/${MODULES}
module use /nopt/nrel/ecom/hpacf/software/${MODULES}/${COMPILER}

module load intel-parallel-studio
module load gcc
module load git
module load python/3.7.4
module load visit/2.13.3-mesa

export dirname=${HOME}/combustion/Pele/PeleC/Exec/Production/TGReact
export paren=$(pwd)
export ranks_per_node=36
export OMP_NUM_THREADS=1  # Max hardware threads = 4
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

echo "Job name       = $SLURM_JOB_NAME"
echo "Num. nodes     = $SLURM_JOB_NUM_NODES"
echo "Num. threads   = $OMP_NUM_THREADS"
echo "Working dir    = $PWD"
