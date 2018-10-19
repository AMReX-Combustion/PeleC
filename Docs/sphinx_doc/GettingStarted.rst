 .. role:: cpp(code)
    :language: c++
 
 .. role:: fortran(code)
    :language: fortran

.. _GettingStarted:

Getting Started
===============

Setting up a problem to run with PeleC involves writing an input file and problem specific code in the run directory. PeleC is built using the AMReX build system which supports out-of-source builds but as configured in Pele requires a specific directory structure. 


In the Pele directory structure, cases are included in the Exec directory, where the directories included in the release include a suite of regression tests and various capability demonstrations:


* Docs    

  * sphinx_doc

* Exec                        

  * 2D_Kernel_PMF_NonIdeal  
  * ColdShock               
  * Combustor_fork          
  * EBDemo3D

    * Cold flow into a 3D dimensional geometry inspired by a gas turbine combustor.    

  * EB_Sphere
  * EB_OblqShock               
  * HIT

    * Box of homogeneous isotropic turbulence    

  * Jet2d                   
  * MMS                     
  * NSCBC_test_cases        
  * PMF          

    * 2D Premixed flame calculation (with inflow/outflow boundary conditions). Flow is vertically upward through the domain with no gravity. The lateral boundaries are periodic. A detailed hydrogen model is used. The solution is initialized with a wrinkled (perturbed) 1D steady flame solution computed using the PREMIX code. Two levels of solution-adaptive refinement are automatically triggered by the presence of the flame intermediate, HO2. 

  * RT                      
  * Sedov                   
  * Sod                     
  * SodF                    
  * Spray                   
  * TG 
    * Taylor-Green vortex.                      
  * UnitTests               
  * ViscTimeStep            

* Source     

  * Src_1d                      
  * Src_2d                      
  * Src_3d                      
  * Src_nd                      
  * Spray                       

* Util                        

  * BLAS                    
  * LAPACK                  
  * VODE                    
  * plot1d                

* constants               

Within each directory are the source files that specify the setup of that particular case. The build system will build source contained in the case directory preferentially.  A few key files that need to be supplied for (most) cases are:

* inputs --- a text file containing parameters that are ready by the ParmParse capability in AMReX. These include things like number of time steps, grid size, output file frequency, which physics to include, etc. A list of available data in the pele group can be found in PeleC/Source/param_includes/pelec_params.H
* probin --- a text file used to include namelists to be read at problem initialization to set values of parameters only used in fortran routines
* Prob_nd.F90 -- A fortran routine that contains routines called at:

  * Initialization (`amrex_probinit`) 
  * To set initial values on the grid (`pc_initdata`)
  * Problem teardown (`pc_prob_close`)

* probdata.f90 -- defines a fortran module (probdata_module) used to store data used only in PeleC fortran routines. Data that needs to accessed from both the fortran and c++ layers should not be stored in here, but rather added to the ParmParse input file and copied to the meth_params_module. 

* GNUMakefile --- in addition to setting options to build profiling, debugging, MPI, OpenMP, Compiler toolchain options, the chemical mechanism, transport model, equation of state model, and use of EB are set here for compile time selection. The GNUMakefile includes the ``Make.PeleC`` file from the `Exec` directory that contains build configuration common across the examples. 

* bc_fill_nd.f90 --- used to set values for boundary conditions on inflow and outflow. 





