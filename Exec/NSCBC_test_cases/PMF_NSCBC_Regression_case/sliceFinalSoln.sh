#!/bin/bash

export pltfile=$1
export outfile=extract_${pltfile}.txt

export vars=`grep X ${pltfile}/Header | tr '\n' ' ' | sed -e 's/^/Temp x_velocity density /g' -e 's/$/pressure/g'`

# Extract variables
./fextract.Linux.gfortran.exe -p ${pltfile} -v "${vars}" -s ${outfile}

# Remove header stuff, replace with TECPLOT header (but need number of points...later)
sed -i -e '3s/^#/VARIABLES = /g' -e 's/  */ /g' -e '1,2d' ${outfile}

# Remove remaining comments and lines with only a space
sed -i -e '/^#/ d' -e '/^ $/ d' ${outfile}

# Count remaining lines, subtract 1 for TECPLOT header line
export nlines=`awk 'END{print NR}' ${outfile}`
let nlines--

# Insert TECPLOT zone info before data
sed -i '2i\'$'\n'"ZONE I=$nlines"$'\n' ${outfile}
