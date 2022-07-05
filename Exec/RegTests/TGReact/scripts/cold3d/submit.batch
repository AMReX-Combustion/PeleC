#!/bin/bash

#SBATCH --job-name=pelec_tgreact
#SBATCH --account=exact
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH -o %x.o%j

source ../set_environment.sh

mpi_ranks=$(expr $SLURM_JOB_NUM_NODES \* $ranks_per_node)
echo "Num. MPI Ranks = $mpi_ranks"

# Location of input files and executable
pele_exec=${dirname}/PeleC3d.intel.MPI.ex
iname="inputs_3d.inp"
input="${dirname}/${iname}"
plot_step=10000

# Compile executable
cd "${dirname}" || exit
make realclean
make -j ${ranks_per_node} DIM=3 USE_MPI=TRUE COMP=intel
cd "${paren}" || exit

# Loop on resolutions
resolutions=(32 64 128 256)
for res in "${resolutions[@]}"
do
    # Setup directory and files
    workdir="${paren}/${res}"
    mkdir -p "${workdir}"
    cd "${workdir}" || exit
    cp "${input}" "${workdir}"
    cp "${pele_exec}" "${workdir}/PeleC"
    rm -rf plt* chk* core* Backtrace.* datlog extralog ic.txt data

    sed -i "/amr.n_cell/c\amr.n_cell=${res} ${res} ${res}" "${iname}"
    sed -i "/amr.check_int/c\amr.check_int=${plot_step}" "${iname}"
    sed -i "/amr.plot_int/c\amr.plot_int=${plot_step}" "${iname}"

    # Run Pele
    srun -n "${mpi_ranks}" -c 1 --cpu_bind=cores "${workdir}/PeleC" "${iname}" > run.out

    # Post process
    visit -np ${ranks_per_node} -nowin -cli -s "${dirname}/visit_pp_aux_vars.py" -i "${iname}"

    # Clean
    cd "${paren}" || exit
done

# Final plots
# conda run -n tgreact python3 "${dirname}/plotter.py" -f ${resolutions[*]} -c cold3d
