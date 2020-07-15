PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 1

COMP	   = gcc

USE_MPI    = TRUE 
USE_OMP    = FALSE

Eos_dir     := GammaLaw
Transport_dir := Constant
Chemistry_Model := Null

Bpack   := ./Make.package
Blocs   := .

# define the location of the PELE top directory
PELEC_HOME    := ../../../..
include $(PELEC_HOME)/Exec/Make.PeleC

ifeq ($(DEBUG), TRUE)
DEFINES += -DDEBUG
endif
