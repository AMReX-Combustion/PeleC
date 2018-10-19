#!/bin/bash

mpiexec -np 8 ./PeleC2d.gnu.MPI.ex inputs-2d amr.n_cell= 512 512 max_step=160 pelec.fixed_dt=1.250e-9
mpiexec -np 8 ./PeleC2d.gnu.MPI.ex inputs-2d amr.n_cell= 256 256 max_step= 80 pelec.fixed_dt=2.500e-9
mpiexec -np 8 ./PeleC2d.gnu.MPI.ex inputs-2d amr.n_cell= 128 128 max_step= 40 pelec.fixed_dt=5.000e-9

~/src/CCSE/amrex/Tools/C_util/Convergence/RichardsonConvergenceTest2d.gnu.ex coarFile=plt00040 mediFile=plt00080 fineFile=plt00160 mediError=medErr coarError=coarErr

#mpiexec -np 8 ./PeleC2d.gnu.MPI.ex inputs-2d amr.n_cell= 1024 1024 max_step=320 pelec.fixed_dt=0.25e-9
#mpiexec -np 8 ./PeleC2d.gnu.MPI.ex inputs-2d amr.n_cell=  512  512 max_step=160 pelec.fixed_dt=0.50e-9
#mpiexec -np 8 ./PeleC2d.gnu.MPI.ex inputs-2d amr.n_cell=  256  256 max_step= 80 pelec.fixed_dt=1.00e-9

#../../../../amrex/Tools/C_util/Convergence/RichardsonConvergenceTest2d.gnu.ex coarFile=plt00080 mediFile=plt00160 fineFile=plt00320 mediError=medErr coarError=coarErr

#mpiexec -np 8 ./PeleC2d.gnu.MPI.ex inputs-2d amr.n_cell=  256  256 max_step=160 pelec.fixed_dt=0.50e-9
#mpiexec -np 8 ./PeleC2d.gnu.MPI.ex inputs-2d amr.n_cell=  128  128 max_step= 80 pelec.fixed_dt=1.00e-9
#mpiexec -np 8 ./PeleC2d.gnu.MPI.ex inputs-2d amr.n_cell=   64   64 max_step= 40 pelec.fixed_dt=2.00e-9

#../../../../amrex/Tools/C_util/Convergence/RichardsonConvergenceTest2d.gnu.ex coarFile=plt00040 mediFile=plt00080 fineFile=plt00160 mediError=medErr coarError=coarErr


#mpiexec -np 8 ./PeleC2d.gnu.MPI.ex inputs-2d amr.n_cell=  256  256 max_step=320 pelec.fixed_dt=0.50e-9
#mpiexec -np 8 ./PeleC2d.gnu.MPI.ex inputs-2d amr.n_cell=  128  128 max_step=160 pelec.fixed_dt=1.00e-9
#mpiexec -np 8 ./PeleC2d.gnu.MPI.ex inputs-2d amr.n_cell=   64   64 max_step= 80 pelec.fixed_dt=2.00e-9

#../../../../amrex/Tools/C_util/Convergence/RichardsonConvergenceTest2d.gnu.ex coarFile=plt00080 mediFile=plt00160 fineFile=plt00320 mediError=medErr coarError=coarErr


