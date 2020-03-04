#/bin/bash

# Example script for generating Fuego code with Fortran
# (note the default pickler must be changed to "f" in 
# Support/Fuego/Pythia/pythia-0.4/packages/fuego/fuego/fmc/FMC.py

CHEMINP=LiDryer.mec
THERMINP=LiDryer.therm
FINALFILE=LiDryer.F90

FMC=${PELE_PHYSICS_HOME}/Support/Fuego/Pythia/products/bin/fmc.py
HEADERDIR=${PELE_PHYSICS_HOME}/Support/Fuego/Mechanism/Models/header

CHEMLK=chem.asc
LOG=chem.log
CHEMC=chem.F90

${FUEGO_PYTHON} ${FMC} -mechanism=${CHEMINP} -thermo=${THERMINP}  -name=${CHEMC}
echo Compiling ${FINALFILE}...
#cat ${CHEMC} \
#          ${HEADERDIR}/header.start\
#          ${HEADERDIR}/header.mec   ${CHEMINP}\
#          ${HEADERDIR}/header.therm ${THERMINP}\
#          ${HEADERDIR}/header.end > ${FINALFILE}
cat ${CHEMC} > ${FINALFILE}
rm -f ${CHEMC} ${CHEMLK} ${LOG} 
sed -i.bak '/\/\*/d' ${FINALFILE}
rm ${FINALFILE}.bak
