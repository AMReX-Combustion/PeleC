# AMReX
DIM = 3
COMP = llvm
PRECISION = DOUBLE

# Profiling
PROFILE = FALSE
TINY_PROFILE = FALSE
COMM_PROFILE = FALSE
TRACE_PROFILE = FALSE
MEM_PROFILE = FALSE
USE_GPROF = FALSE

# Performance
USE_MPI = FALSE
USE_OMP = FALSE
USE_CUDA = FALSE
USE_HIP = FALSE
USE_DPCPP = FALSE

# Debugging
DEBUG = FALSE
FSANITIZER = FALSE
THREAD_SANITIZER = FALSE

# PeleC
USE_REACT = TRUE
USE_SUNDIALS_PP = FALSE
USE_EB = FALSE
USE_MASA = FALSE
Eos_dir := Fuego
Chemistry_Model := LiDryer
Transport_dir := Simple

#USE_HDF5  = FALSE
#HDF5_HOME = /Users/jrood/spack/opt/spack/darwin-catalina-x86_64/apple-clang-11.0.3/hdf5-1.10.4-nb5dk5var5kf7kbeekh6rt3dfytbdf2f

# GNU Make
Bpack := ./Make.package
Blocs := .
PELEC_HOME := ../../..
include $(PELEC_HOME)/Exec/Make.PeleC
