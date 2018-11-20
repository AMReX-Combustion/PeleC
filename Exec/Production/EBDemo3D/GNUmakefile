PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = TRUE
DEBUG      = FALSE

#COMP       = Intel

#DIM        = 1
#DIM        = 2
DIM        = 3

USE_MPI    = TRUE
USE_OMP    = FALSE

USE_EB = TRUE

HYP_TYPE = MOL
#TINY_PROFILE = TRUE
#PROFILE = TRUE

# define the location of the PELE top directory
PELE_HOME    := ../../..

# This sets the EOS directory in $(PELE_PHYSICS_HOME)/Eos
Eos_dir     := GammaLaw

# This sets the network directory in $(PELE_PHYSICS_HOME)/Reactions
Reactions_dir := Null_air

# This sets the transport directory in $(PELE_PHYSICS_HOME)/Transport
Transport_dir := Constant

Bpack   := ./Make.package
Blocs   := .


include ../../Make.PeleC
