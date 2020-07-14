PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

COMP       = gnu

DIM        = 3

USE_MPI    = FALSE
USE_OMP    = FALSE

USE_EB = TRUE
USE_REACT  = TRUE

HYP_TYPE = MOL
#TINY_PROFILE = TRUE
#PROFILE = TRUE

Eos_dir     := Fuego
Reactions_dir := Fuego
Chemistry_Model = air

# This sets the transport directory in $(PELE_PHYSICS_HOME)/Transport
Transport_dir := EGLib
#Transport_dir := Constant

Bpack   := ./Make.package
Blocs   := .

# define the location of the PELE top directory
PELEC_HOME := ../../..
include $(PELEC_HOME)/Exec/Make.PeleC
