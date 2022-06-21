# AMReX
DIM = 3
COMP = gnu
PRECISION = DOUBLE

# Profiling
PROFILE = FALSE
TINY_PROFILE = TRUE
COMM_PROFILE = FALSE
TRACE_PROFILE = FALSE
MEM_PROFILE = FALSE
USE_GPROF = FALSE

# Performance
USE_MPI = TRUE
USE_OMP = FALSE
USE_CUDA = FALSE
USE_HIP = FALSE

# Debugging
DEBUG = FALSE
FSANITIZER = FALSE
THREAD_SANITIZER = FALSE

# PeleC
Eos_Model := Fuego
Chemistry_Model := dodmethair_4sp
Transport_Model := Simple
USE_MASA = FALSE

# Libraries
USE_HDF5  = FALSE
HDF5_HOME = $(HDF5_ROOT)
# Or on Summit:
# HDF5_HOME = /gpfs/alpine/csc300/world-shared/gnu_build/hdf5/develop_build/hdf5
# To use HDF5 with ZFP compression, follow the instructions at https://h5z-zfp.readthedocs.io
# to compile ZFP and H5Z_ZFP plugin, or use the following gcc build on Summit:
USE_HDF5_ZFP = FALSE
H5Z_HOME  = /gpfs/alpine/csc300/world-shared/gnu_build/H5Z-ZFP/install
ZFP_HOME  = /gpfs/alpine/csc300/world-shared/gnu_build/zfp

# GNU Make
Bpack := ./Make.package
Blocs := .
PELEC_HOME := ../../..
include $(PELEC_HOME)/Exec/Make.PeleC
