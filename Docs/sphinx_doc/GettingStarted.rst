 .. role:: cpp(code)
    :language: c++
 
 .. role:: fortran(code)
    :language: fortran

.. _GettingStarted:

Getting Started
===============

Navigation
----------

The PeleC directory structure is as shown below:

* **Source** - C++ and fortran source code

  * Src_1d                      
  * Src_2d                      
  * Src_3d                      
  * Src_nd                      
  * Spray                       

* **Util** - third party utilities

  * BLAS                    
  * LAPACK                  
  * VODE                    
  * plot1d                

* **constants** - fundamental constants in CGS units


* **Docs**   - PeleC documentation 

  * sphinx_doc

* **Exec** - regression tests and various capability demonstrations
  
  * :ref:`Regression Tests<VandV>`

    * HIT  - Box of homogeneous isotropic turbulence 
    * MMS   - Method of manufactured solutions
    * PMF   - 2D Premixed H2-air flame calculation 
    * Sedov - Sedov blast wave test
    * Sod   - Sod shock tube test
    * TG    - Taylor-Green vortex test

  * UnitTests

    * NSCBC_test_cases - This directory contains many 1D, 2D and 3D cases to test the implementation of the Ghost-Cells Navier-Stokes Boundary Conditions (GC-NSCBC) on various different configurations.
    * HIT_forced - This test case is similar to the homogeneous isotropic turbulence present in the RegTests directory, at the exception that we are not starting from a turbulent initial solution, but from a flow at rest where we superimpose forcing sources to generate turbulence.
  
  * :ref:`Tutorials<Tutorials>`

    * EB_Sphere - Reacting flow around a sphere
    * EB_OblqShock - Supersonic flow over a wedge resulting in a steady attached oblique shock
    * EB_BluffBody -  Non-reacting flow around a diamond-shaped body
    * Spray


Setting up a problem to run with PeleC involves writing an input file and problem specific code in the run directory. 
PeleC is built using the AMReX build system which supports out-of-source builds but as configured in Pele requires a specific directory structure. 
Within each case directory in Exec, are the source files that specify the setup of that particular case. 
The user has to build each case by compiling source files using a GNUMakefile which also compiles and links together AMReX and PeleC sources.
The source files contained in the case directory are treated preferentially and can override PeleC/AMReX source files.  
A few key files that need to be supplied for (most) cases are:

**inputs** -- a text file containing parameters that are ready by the ParmParse capability in AMReX. These include things like number of time steps, grid size, output file frequency, which physics to include, etc. 
A list of available data in the Pele group can be found in PeleC/Source/param_includes/pelec_params.H

**probin** -- a text file used to include namelists to be read at problem initialization to set values of parameters only used in fortran routines

**Prob_nd.F90** -- A fortran routine that contains routines called at:

  * Initialization (`amrex_probinit`) 
  * To set initial values on the grid (`pc_initdata`)
  * Problem teardown (`pc_prob_close`)

**probdata.f90** -- defines a fortran module (probdata_module) used to store data used only in PeleC fortran routines. Data that needs to accessed from both the 
fortran and c++ layers should not be stored in here, but rather added to the ParmParse input file and copied to the meth_params_module. 

**GNUMakefile** -- in addition to setting options to build profiling, debugging, MPI, OpenMP, Compiler toolchain options, the chemical mechanism, transport model, equation of state model, 
and use of EB are set here for compile time selection. The GNUMakefile includes the ``Make.PeleC`` file from the `Exec` directory that contains build configuration common across the examples. 

**bc_fill_nd.f90** -- used to set values for user defined boundary conditions that use GC-NSCBC. 


.. include:: building.rst


.. include:: InputFiles.rst


.. include:: tutorials.rst


.. include:: testing.rst
