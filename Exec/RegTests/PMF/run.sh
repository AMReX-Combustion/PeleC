#!/bin/bash -l

#SBATCH -J pelec-chem-integrators-gpu
#SBATCH -o %x.o%j
#SBATCH -t 01:00:00
#SBATCH -N 1
#SBATCH -p debug
#SBATCH -A exact
#SBATCH --gres=gpu:2

module load mpt cuda/10.2.89
export MPI_IGNORE_PBS=on

set -x

mpirun -np 1 ./PeleC3d.gnu.TPROF.MPI.CUDA.ex inputs_ex amr.n_cell=64 64 128 geometry.prob_lo=-1.25 -1.25 -1.5 geometry.prob_hi=1.25 1.25 3.5 amr.max_level=1 max_step=10 amr.plot_files_output=0 amr.plot_int=20 amr.checkpoint_files_output=0 pelec.do_mol=0 pelec.ppm_type=1 prob.pamb=1013250.0 prob.phi_in=-0.5 prob.pertmag=0.005 prob.pmf_datafile="PMF_CH4_1bar_300K_DRM_MixAvg.dat" tagging.max_ftracerr_lev=4 tagging.ftracerr=1.0e-4 ode.rtol=1.0e-6 ode.atol=1.0e-10 pelec.cfl=0.1 pelec.init_shrink=1.0 pelec.change_max=1.0 pelec.chem_integrator=1 ode.use_erkstep=0
mpirun -np 1 ./PeleC3d.gnu.TPROF.MPI.CUDA.ex inputs_ex amr.n_cell=64 64 128 geometry.prob_lo=-1.25 -1.25 -1.5 geometry.prob_hi=1.25 1.25 3.5 amr.max_level=1 max_step=10 amr.plot_files_output=0 amr.plot_int=20 amr.checkpoint_files_output=0 pelec.do_mol=0 pelec.ppm_type=1 prob.pamb=1013250.0 prob.phi_in=-0.5 prob.pertmag=0.005 prob.pmf_datafile="PMF_CH4_1bar_300K_DRM_MixAvg.dat" tagging.max_ftracerr_lev=4 tagging.ftracerr=1.0e-4 ode.rtol=1.0e-6 ode.atol=1.0e-10 pelec.cfl=0.1 pelec.init_shrink=1.0 pelec.change_max=1.0 pelec.chem_integrator=2 ode.use_erkstep=0
mpirun -np 1 ./PeleC3d.gnu.TPROF.MPI.CUDA.ex inputs_ex amr.n_cell=64 64 128 geometry.prob_lo=-1.25 -1.25 -1.5 geometry.prob_hi=1.25 1.25 3.5 amr.max_level=1 max_step=10 amr.plot_files_output=0 amr.plot_int=20 amr.checkpoint_files_output=0 pelec.do_mol=0 pelec.ppm_type=1 prob.pamb=1013250.0 prob.phi_in=-0.5 prob.pertmag=0.005 prob.pmf_datafile="PMF_CH4_1bar_300K_DRM_MixAvg.dat" tagging.max_ftracerr_lev=4 tagging.ftracerr=1.0e-4 ode.rtol=1.0e-6 ode.atol=1.0e-10 pelec.cfl=0.1 pelec.init_shrink=1.0 pelec.change_max=1.0 pelec.chem_integrator=3 ode.use_erkstep=0
mpirun -np 1 ./PeleC3d.gnu.TPROF.MPI.CUDA.ex inputs_ex amr.n_cell=64 64 128 geometry.prob_lo=-1.25 -1.25 -1.5 geometry.prob_hi=1.25 1.25 3.5 amr.max_level=1 max_step=10 amr.plot_files_output=0 amr.plot_int=20 amr.checkpoint_files_output=0 pelec.do_mol=0 pelec.ppm_type=1 prob.pamb=1013250.0 prob.phi_in=-0.5 prob.pertmag=0.005 prob.pmf_datafile="PMF_CH4_1bar_300K_DRM_MixAvg.dat" tagging.max_ftracerr_lev=4 tagging.ftracerr=1.0e-4 ode.rtol=1.0e-6 ode.atol=1.0e-10 pelec.cfl=0.1 pelec.init_shrink=1.0 pelec.change_max=1.0 pelec.chem_integrator=2 ode.use_erkstep=1
mpirun -np 1 ./PeleC3d.gnu.TPROF.MPI.CUDA.ex inputs_ex amr.n_cell=64 64 128 geometry.prob_lo=-1.25 -1.25 -1.5 geometry.prob_hi=1.25 1.25 3.5 amr.max_level=1 max_step=10 amr.plot_files_output=0 amr.plot_int=20 amr.checkpoint_files_output=0 pelec.do_mol=0 pelec.ppm_type=1 prob.pamb=1013250.0 prob.phi_in=-0.5 prob.pertmag=0.005 prob.pmf_datafile="PMF_CH4_1bar_300K_DRM_MixAvg.dat" tagging.max_ftracerr_lev=4 tagging.ftracerr=1.0e-4 ode.rtol=1.0e-6 ode.atol=1.0e-10 pelec.cfl=0.1 pelec.init_shrink=1.0 pelec.change_max=1.0 pelec.chem_integrator=3 ode.use_erkstep=1
