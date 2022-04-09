#!/bin/bash

mpi_ranks=36

for i_in in 0 1; do
    for i_out in 0 1; do
	dirname=run_$i_in$i_out 
	echo $dirname
	mkdir $dirname
	cp example.inp $dirname/run.inp
	cd $dirname
	sed -i "s/prob.inlet_type = 0/prob.inlet_type = $i_in/g" run.inp
	sed -i "s/prob.outlet_type = 0/prob.outlet_type = $i_out/g" run.inp
	srun -n ${mpi_ranks} ../PeleC3d.gnu.MPI.ex run.inp
	cd ..
    done
done

# change path to AMReX fextract utility executable as appropriate
../../../Submodules/AMReX/Tools/Plotfile/fextract.gnu.ex run_01/plt24000

# Requires: Python3 Numpy Scipy Pandas Matplotlib
python plot.py
