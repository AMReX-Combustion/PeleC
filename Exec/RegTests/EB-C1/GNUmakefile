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
Eos_Model := GammaLaw
Chemistry_Model := Null
Transport_Model := Constant

# MASA
USE_MASA = TRUE
ifeq ($(USE_MASA), TRUE)
  INCLUDE_LOCATIONS += ${MASA_HOME}/include ${MASA_HOME}/lib
  LIBRARY_LOCATIONS += ${MASA_HOME}/lib
  LIBRARIES += -Wl,-rpath,${MASA_HOME}/lib -lmasa
  DEFINES += -DPELEC_USE_MASA
endif

# GNU Make
Bpack := ./Make.package
Blocs := .
PELEC_HOME := ../../..
include $(PELEC_HOME)/Exec/Make.PeleC
