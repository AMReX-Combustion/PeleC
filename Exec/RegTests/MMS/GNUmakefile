PRECISION  = DOUBLE
PROFILE    = FALSE
USE_EB     = FALSE
DEBUG      = FALSE

DIM        = 3

COMP	   = gcc

USE_MPI    = FALSE
USE_OMP    = FALSE

Eos_dir     := GammaLaw
Transport_dir := Constant
Chemistry_Model := Null

# MASA
USE_MASA   = TRUE
ifeq ($(USE_MASA), TRUE)
INCLUDE_LOCATIONS += ${MASA_HOME}/include ${MASA_HOME}/lib
LIBRARY_LOCATIONS += ${MASA_HOME}/lib
LIBRARIES += -Wl,-rpath,${MASA_HOME}/lib -lmasa -lfmasa
DEFINES += -DUSE_MASA -DDO_PROBLEM_POST_TIMESTEP -DDO_PROBLEM_POST_INIT
endif

Bpack   := ./Make.package

# define the location of the PELE top directory
PELEC_HOME    := ../../..
include $(PELEC_HOME)/Exec/Make.PeleC

ifeq ($(DEBUG), TRUE)
DEFINES += -DDEBUG
endif
