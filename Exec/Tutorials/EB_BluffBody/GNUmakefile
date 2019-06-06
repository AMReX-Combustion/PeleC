PRECISION  = DOUBLE
PROFILE    = FALSE
DEBUG      = FALSE

DIM        = 3

COMP	   = gnu

USE_MPI    = TRUE
USE_OMP    = FALSE

USE_EB = TRUE

HYP_TYPE = MOL

# This sets the EOS directory in $(PELE_PHYSICS_HOME)/Eos
Eos_dir     := GammaLaw

# This sets the network directory in $(PELE_PHYSICS_HOME)/Reactions
Reactions_dir := Null_air

# This sets the transport directory in $(PELE_PHYSICS_HOME)/Transport
Transport_dir := Constant

Bpack   := ./Make.package
Blocs   := .

# define the location of the PELE top directory
PELEC_HOME := ../../..
include $(PELEC_HOME)/Exec/Make.PeleC
