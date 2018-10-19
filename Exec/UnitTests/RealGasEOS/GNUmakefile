PRECISION  = DOUBLE

PROFILE    = FALSE

DEBUG      = TRUE

DIM        = 1

COMP	   = gcc

USE_MPI    = TRUE

USE_REACT  = TRUE

# define the location of the PELE top directory
PELE_HOME  := ../../../..

# This sets the EOS directory in $(PELE_PHYSICS_HOME)/Eos
Eos_dir     := Soave-Redlich-Kwong

# This sets the network directory in $(PELE_PHYSICS_HOME)/Reactions
Reactions_dir := Fuego

# This sets the model directory in $(PELE_PHYSICS_HOME)/Support/Fuego/Mechanism/Models
Chemistry_Model := grimech30

# This sets the transport directory in $(PELE_PHYSICS_HOME)/Transport
Transport_dir := Simple

Bpack   := ./Make.package
Blocs   := .

include ../../Make.PeleC
