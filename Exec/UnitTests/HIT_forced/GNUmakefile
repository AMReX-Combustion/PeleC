PRECISION  = DOUBLE
PROFILE    = FALSE
USE_EB     = FALSE
DEBUG      = FALSE

DIM        = 3 

COMP	   = gcc

USE_MPI    = TRUE
USE_OMP    = FALSE

#HYP_TYPE = MOL

Eos_dir     := GammaLaw
Transport_dir := Constant
Chemistry_Model := Null

Bpack   := ./Make.package

# define the location of the PELE top directory
PELEC_HOME    := ../../..
include $(PELEC_HOME)/Exec/Make.PeleC

DEFINES += -DO_HIT_FORCE

ifeq ($(DEBUG), TRUE)
DEFINES += -DDEBUG
endif
