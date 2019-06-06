PRECISION  = DOUBLE
USE_EB     = FALSE
PROFILE    = FALSE

DEBUG      = TRUE
DEBUG      = FALSE

DIM        = 1
#DIM        = 2
#DIM        = 3

COMP	   = gcc

USE_MPI    = TRUE
USE_OMP    = FALSE

#HYP_TYPE = MOL

# This sets the EOS directory in $(PELE_PHYSICS_HOME)/Eos
Eos_dir     := GammaLaw

# This sets the network directory in $(PELE_PHYSICS_HOME)/Reactions
Reactions_dir := Null

# This sets the transport directory in $(PELE_PHYSICS_HOME)/Transport
Transport_dir := Constant

Bpack   := ./Make.package
Blocs   := .

# define the location of the PELE top directory
PELEC_HOME    := ../../..
include $(PELEC_HOME)/Exec/Make.PeleC
