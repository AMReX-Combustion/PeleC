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
USE_SYCL = FALSE

# Debugging
DEBUG = FALSE
FSANITIZER = FALSE
THREAD_SANITIZER = FALSE

# PeleC
USE_EXPLICIT_REACT = FALSE
ifeq ($(USE_EXPLICIT_REACT), TRUE)
Eos_Model := Fuego
Chemistry_Model := LiDryer
Transport_Model := Simple
else
Eos_Model := GammaLaw
Chemistry_Model := Null
Transport_Model := Constant
endif

# GNU Make
Bpack := ./Make.package
Blocs := .
PELEC_HOME := ../../..
include $(PELEC_HOME)/Exec/Make.PeleC
