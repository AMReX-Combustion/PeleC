PRECISION  = DOUBLE
PROFILE    = FALSE
USE_EB = FALSE
DEBUG      = TRUE
DEBUG      = FALSE

DIM        = 3
DIM        = 2

COMP	   = gcc

USE_MPI    = FALSE
USE_MPI    = TRUE

USE_REACT  = TRUE

# define the location of the PELE top directory
PELE_HOME  := ../../..

# This sets the EOS directory in $(PELE_PHYSICS_HOME)/Eos
Eos_dir     := Soave-Redlich-Kwong 

# This sets the network directory in $(PELE_PHYSICS_HOME)/Reactions
Reactions_dir := Fuego

# This sets the model directory in $(PELE_PHYSICS_HOME)/Support/Fuego/Mechanism/Models
Chemistry_Model := LiDryer
#Chemistry_Model := BurkeDryer_mod

# This sets the transport directory in $(PELE_PHYSICS_HOME)/Transport
Transport_dir := Simple

Bpack   := ./Make.package
Blocs   := .

include ../Make.PeleC
