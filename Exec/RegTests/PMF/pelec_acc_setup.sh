#!/bin/bash -l

# Example script for setting up PeleC acc_prototype branch on NREL's Eagle machine

cmd() {
  echo "+ $@"
  eval "$@"
}

# cmd "mkdir ${HOME}/combustion && cd ${HOME}/combustion" or whatever
cmd "module unload gcc && module load pgi/18.10"
cmd "git clone -b acc_prototype --recurse-submodules -j 4 https://github.com/AMReX-Combustion/PeleC.git"
cmd "cd PeleC/Submodules/AMReX && git apply ../../Exec/RegTests/PMF/amrex_acc_prototype.patch && cd -"
cmd "cd PeleC/Exec/RegTests/PMF && nice make -j16"

# Notes:
# Edit PeleC/Exec/RegTests/PMF/GNUmakefile to change USE_ACC=TRUE/FALSE
# "make realclean" is the clean command
