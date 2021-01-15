PRECISION  = DOUBLE
PROFILE    = FALSE
USE_EB     = FALSE
DEBUG      = FALSE
DIM        = 3 
COMP	   = gnu
USE_MPI    = FALSE
USE_OMP    = FALSE

Eos_dir     := GammaLaw
Transport_dir := Constant
Chemistry_Model = Null

Bpack   := ./Make.package

PELEC_HOME    := ../../..
include $(PELEC_HOME)/Exec/Make.PeleC

ifeq ($(DEBUG), TRUE)
DEFINES += -DDEBUG
endif
