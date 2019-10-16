PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 2

COMP	   = gcc

USE_MPI    = FALSE
USE_OMP    = FALSE
USE_GRAV   = FALSE

# Include support for participating or tracer particles
USE_PARTICLES = TRUE

SPRAY_COMPONENTS = 8

DEFINES += -DSPRAY_COMPONENTS=$(SPRAY_COMPONENTS)

# define the location of the PELE top directory
#PELEC_HOME    := ../../..

# This sets the EOS directory in $(PELE_PHYSICS_HOME)/Eos
Eos_dir := Fuego
#Eos_dir := GammaLaw

# This sets the network directory in $(PELE_PHYSICS_HOME)/Reactions
Reactions_dir := Fuego
#Reactions_dir := Null
#
# This sets the model directory in $(PELE_PHYSICS_HOME)/Support/Fuego/Mechanism/Models
#Chemistry_Model := chem-H
#Chemistry_Model := grimech12
Chemistry_Model := heptane_3sp
#Chemistry_Model := decane_3sp

# This sets the transport directory in $(PELE_PHYSICS_HOME)/Transport
Transport_dir := EGLib
#Transport_dir := Constant

Bpack   := ./Make.package 
Blocs   := .

# define the location of the PELE top directory
PELEC_HOME    := ../../../..
include $(PELEC_HOME)/Exec/Make.PeleC

ifeq ($(DEBUG), TRUE)
DEFINES += -DDEBUG
endif
