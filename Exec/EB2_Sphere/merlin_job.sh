#!/bin/bash -l
 
#PBS -j oe
#PBS -o pele.out
#PBS -N PeleC
#PBS -l walltime=12:00:00
###PBS -l select=10:quad-cache=true:ncpus=256
#PBS -l select=10:ncpus=24
###PBS -q __PBS_QUEUE__
###PBS -A __PBS_ACCOUNT__
###PBS -W __PBS_UMASK__
###PBS -Wsandbox=HOME
#PBS -V
 
date
printenv
 
cd $PBS_O_WORKDIR
cat $PBS_NODEFILE | sort | uniq > job_nodes
 
 
mpirun -np 240 ./PeleC3d.intel.MPI.ex inputs-ebdemo   > job_output 2>&1
