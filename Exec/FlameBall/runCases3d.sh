#!/bin/bash

#mpiexec -np 8 ./PeleC3d.gnu.MPI.ex inputs-3d amr.n_cell=128 128 128 max_step=80 pelec.fixed_dt=1.250e-9
#mpiexec -np 8 ./PeleC3d.gnu.MPI.ex inputs-3d amr.n_cell= 64  64  64 max_step=40 pelec.fixed_dt=2.500e-9
#mpiexec -np 8 ./PeleC3d.gnu.MPI.ex inputs-3d amr.n_cell= 32  32  32 max_step=20 pelec.fixed_dt=5.000e-9

#../../../amrex/Tools/C_util/Convergence/RichardsonConvergenceTest3d.gnu.ex coarFile=plt00020 mediFile=plt00040 fineFile=plt00080 mediError=medErr coarError=coarErr

mpiexec -np 8 ./PeleC3d.gnu.MPI.ex inputs-3d amr.n_cell=128 128 128 max_step=80 pelec.fixed_dt=1.0e-9
mpiexec -np 8 ./PeleC3d.gnu.MPI.ex inputs-3d amr.n_cell= 64  64  64 max_step=40 pelec.fixed_dt=2.00e-9
mpiexec -np 8 ./PeleC3d.gnu.MPI.ex inputs-3d amr.n_cell= 32  32  32 max_step=20 pelec.fixed_dt=4.000e-9

../../../../amrex/Tools/C_util/Convergence/RichardsonConvergenceTest3d.gnu.ex coarFile=plt00020 mediFile=plt00040 fineFile=plt00080 mediError=medErr coarError=coarErr

