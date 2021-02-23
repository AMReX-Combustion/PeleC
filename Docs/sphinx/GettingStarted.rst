 .. role:: cpp(code)
    :language: c++
 
.. _GettingStarted:

Getting Started
===============

Navigation
----------

The PeleC directory structure is as shown below:

* **Source** - C++ source code

* **Util** - third party utilities

  * BLAS
  * LAPACK
  * VODE
  * plot1d

* **Docs**   - PeleC documentation 

  * sphinx

* **Exec** - regression tests and various capability demonstrations
  
  * :ref:`Regression Tests<VandV>`

    * EB_MMS   - Method of manufactured solutions
    * HIT      - Box of homogeneous isotropic turbulence 
    * MMS      - Method of manufactured solutions
    * PMF      - 2D Premixed H2-air flame calculation 
    * Sedov    - Sedov blast wave test
    * Sod      - Sod shock tube test
    * TG       - Taylor-Green vortex test
    * zeroD    - Single cell test for reaction model

  * UnitTests

    * NSCBC_test_cases - This directory contains many 1D, 2D and 3D cases to test the implementation of the Ghost-Cells Navier-Stokes Boundary Conditions (GC-NSCBC) on various different configurations.
    * HIT_forced - This test case is similar to the homogeneous isotropic turbulence present in the RegTests directory, at the exception that we are not starting from a turbulent initial solution, but from a flow at rest where we superimpose forcing sources to generate turbulence.
  
  * :ref:`Tutorials`

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

**prob.cpp** -- Routines called at:

  * Initialization (`amrex_probinit`) 
  * To set initial values on the grid (`pc_initdata`)
  * Problem teardown (`pc_prob_close`)

**prob.H** -- Something about prob.H

**prob_parm.H** -- Something about prob_parm.H

**GNUMakefile** -- In addition to setting options to build profiling, debugging, MPI, OpenMP, Compiler toolchain options, the chemical mechanism, transport model, equation of state model, and use of EB are set here for compile time selection. The GNUMakefile includes the ``Make.PeleC`` file from the `Exec` directory that contains build configuration common across the examples.


.. toctree::
   :maxdepth: 1

   Building
   InputFiles
   Tutorials
   Testing
