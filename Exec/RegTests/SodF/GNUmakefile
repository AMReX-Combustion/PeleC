PRECISION  = DOUBLE
PROFILE    = FALSE
USE_EB     = FALSE

DEBUG      = TRUE
DEBUG      = FALSE


#DIM        = 1
DIM        = 2
#DIM        = 3

COMP	   = gcc

USE_MPI    = TRUE
USE_OMP    = FALSE

# define the location of the PELE top directory
PELE_HOME    := ../../..

# This sets the EOS directory in $(PELE_PHYSICS_HOME)/Eos
Eos_dir     := Fuego

# This sets the network directory in $(PELE_PHYSICS_HOME)/Reactions
Reactions_dir := Fuego

# This sets the model directory in $(PELE_PHYSICS_HOME)/Support/Fuego/Mechanism/Models
Chemistry_Model := chem-H

# This sets the transport directory in $(PELE_PHYSICS_HOME)/Transport
Transport_dir := EGLib

Bpack   := ./Make.package
Blocs   := .

include ../../Make.PeleC
