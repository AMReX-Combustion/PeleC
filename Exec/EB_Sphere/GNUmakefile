PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = TRUE


#DIM        = 2
DIM        = 3

COMP	   = gcc

USE_MPI    = FALSE
USE_OMP    = FALSE

# Flags to use MOL EB formulation
DEFINES += -DPELEC_USE_EB -DPELE_USE_EB -DPELEC_USE_MOL

HYP_TYPE = MOL
TRANSPORT_TYPE = IDEAL_GAS
USE_EB = TRUE

# define the location of the PELE top directory
PELE_HOME    := ../../..

# This sets the EOS directory in $(PELE_PHYSICS_HOME)/Eos
Eos_dir     := GammaLaw

# This sets the network directory in $(PELE_PHYSICS_HOME)/Reactions
Reactions_dir := Null

# This sets the transport directory in $(PELE_PHYSICS_HOME)/Transport
Transport_dir := Constant

Bpack   := ./Make.package
Blocs   := .

include ../Make.PeleC
