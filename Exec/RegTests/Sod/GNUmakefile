# AMReX
DIM = 3
COMP = gnu
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
Eos_Model := GammaLaw
Transport_Model := Constant
Chemistry_Model := Null

# GNU Make
Bpack := ./Make.package
Blocs := .
PELEC_HOME := ../../..
include $(PELEC_HOME)/Exec/Make.PeleC
