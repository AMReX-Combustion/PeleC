#!/bin/bash

# ./PeleC1d.gnu.MPI.ex inputs-1d amr.n_cell= 512 max_step=160 pelec.fixed_dt=0.625e-9
# ./PeleC1d.gnu.MPI.ex inputs-1d amr.n_cell= 256 max_step= 80 pelec.fixed_dt=1.250e-9
# ./PeleC1d.gnu.MPI.ex inputs-1d amr.n_cell= 128 max_step= 40 pelec.fixed_dt=2.500e-9

# ./PeleC1d.gnu.MPI.ex inputs-1d amr.n_cell= 512 max_step=160 pelec.fixed_dt=1.250e-9
# ./PeleC1d.gnu.MPI.ex inputs-1d amr.n_cell= 256 max_step= 80 pelec.fixed_dt=2.500e-9
# ./PeleC1d.gnu.MPI.ex inputs-1d amr.n_cell= 128 max_step= 40 pelec.fixed_dt=5.000e-9

# ../../../amrex/Tools/C_util/Convergence/RichardsonConvergenceTest1d.gnu.DEBUG.ex coarFile=plt00040 mediFile=plt00080 fineFile=plt00160 mediError=medErr coarError=coarErr


# ./PeleC1d.gnu.MPI.ex inputs-1d amr.n_cell= 512 max_step= 80 pelec.fixed_dt=0.625e-9
# ./PeleC1d.gnu.MPI.ex inputs-1d amr.n_cell= 256 max_step= 40 pelec.fixed_dt=1.250e-9
# ./PeleC1d.gnu.MPI.ex inputs-1d amr.n_cell= 128 max_step= 20 pelec.fixed_dt=2.500e-9

# ./PeleC1d.gnu.MPI.ex inputs-1d amr.n_cell= 512 max_step= 80 pelec.fixed_dt=1.250e-9
# ./PeleC1d.gnu.MPI.ex inputs-1d amr.n_cell= 256 max_step= 40 pelec.fixed_dt=2.500e-9
# ./PeleC1d.gnu.MPI.ex inputs-1d amr.n_cell= 128 max_step= 20 pelec.fixed_dt=5.000e-9

# ../../../amrex/Tools/C_util/Convergence/RichardsonConvergenceTest1d.gnu.DEBUG.ex coarFile=plt00020 mediFile=plt00040 fineFile=plt00080 mediError=medErr coarError=coarErr


# ./PeleC1d.gnu.MPI.ex inputs-1d amr.n_cell= 512 max_step=320 pelec.fixed_dt=1.250e-9
# ./PeleC1d.gnu.MPI.ex inputs-1d amr.n_cell= 256 max_step=160 pelec.fixed_dt=2.500e-9
# ./PeleC1d.gnu.MPI.ex inputs-1d amr.n_cell= 128 max_step= 80 pelec.fixed_dt=5.000e-9

# ../../../amrex/Tools/C_util/Convergence/RichardsonConvergenceTest1d.gnu.DEBUG.ex coarFile=plt00080 mediFile=plt00160 fineFile=plt00320 mediError=medErr coarError=coarErr

./PeleC1d.gnu.MPI.ex inputs-1d amr.n_cell= 512 max_step=640 pelec.fixed_dt=1.250e-9
./PeleC1d.gnu.MPI.ex inputs-1d amr.n_cell= 256 max_step=320 pelec.fixed_dt=2.500e-9
./PeleC1d.gnu.MPI.ex inputs-1d amr.n_cell= 128 max_step=160 pelec.fixed_dt=5.000e-9

../../../amrex/Tools/C_util/Convergence/RichardsonConvergenceTest1d.gnu.DEBUG.ex coarFile=plt00160 mediFile=plt00320 fineFile=plt00640 mediError=medErr coarError=coarErr

