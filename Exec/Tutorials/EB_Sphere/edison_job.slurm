#!/bin/bash -l
#SBATCH -p debug
#SBATCH -N 10
#SBATCH -t 0:30:00
#SBATCH -J PeleC
#SBATCH -o PeleC.o%j
#SBATCH -L SCRATCH
#SBATCH -A m2860

srun -n 240 ./PeleC3d.intel.ivybridge.MPI.ex inputs-ebdemo 
#srun -n 240 ./PeleC3d.gnu.ivybridge.mpi.ex inputs-ebdemo 

