USE_EB = FALSE

PRECISION  = DOUBLE

PROFILE    = FALSE

DEBUG      = FALSE
#DEBUG      = TRUE

#DIM        = 1
DIM        = 2

COMP	   = gcc

#USE_MPI    = TRUE

USE_REACT  = TRUE

#HYP_TYPE   = MOL

# define the location of the PELE top directory
PELE_HOME  := ../../..

# This sets the EOS directory in $(PELE_PHYSICS_HOME)/Eos
#Eos_dir     := Soave-Redlich-Kwong
Eos_dir     := Fuego

# This sets the network directory in $(PELE_PHYSICS_HOME)/Reactions
Reactions_dir := Fuego

# This sets the model directory in $(PELE_PHYSICS_HOME)/Support/Fuego/Mechanism/Models
#Chemistry_Model := Davis
Chemistry_Model := LiDryer

# This sets the transport directory in $(PELE_PHYSICS_HOME)/Transport
Transport_dir := EGLib
#Transport_dir := Simple

Bpack   := ./Make.package
Blocs   := .

include ${PELEC_HOME}/Exec/Make.PeleC
