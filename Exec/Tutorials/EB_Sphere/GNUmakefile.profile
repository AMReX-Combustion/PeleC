PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = TRUE
DEBUG      = FALSE

#COMP       = Intel

#DIM        = 1
#DIM        = 2
DIM        = 3

USE_MPI    = TRUE
USE_OMP    = TRUE
USE_OMP    = FALSE


TINY_PROFILE = FALSE
PROFILE=TRUE
MEM_PROFILE=TRUE
TRACE_PROFILE=TRUE
COMM_PROFILE=TRUE
USE_PROFPARSER=TRUE
COMP = intel

DEBUG = FALSE


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
