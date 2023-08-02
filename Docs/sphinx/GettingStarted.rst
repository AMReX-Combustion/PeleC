 .. role:: cpp(code)
    :language: c++
 
.. _GettingStarted:

Getting Started
===============

Navigation
----------

The PeleC directory structure is as shown below:

* **Source** - C++ source code

* **Docs**   - PeleC documentation

* **Exec** - `regression tests<VandV>` and various capability demonstrations

Setting up a problem to run with PeleC involves writing an input file and problem specific code in the run directory. 
PeleC is built using the AMReX build system which supports out-of-source builds but as configured in Pele requires a specific directory structure. 
Within each case directory in Exec, are the source files that specify the setup of that particular case. 
The user has to build each case by compiling source files using a GNUMakefile which also compiles and links together AMReX and PeleC sources.
The source files contained in the case directory are treated preferentially and can override PeleC/AMReX source files. A few key files that need to be supplied for (most) cases are:

* **inputs,inp** -- a text file containing parameters that are ready by the ParmParse capability in AMReX. These include things like number of time steps, grid size, output file frequency, which physics to include, etc. A list of available data in the Pele group can be found in PeleC/Source/param_includes/pelec_params.H.

* **prob.cpp** -- Routines called at:

  * Initialization (`amrex_probinit`) 
  * To set initial values on the grid (`pc_initdata`)
  * Problem teardown (`pc_prob_close`)

* **prob.H** -- prob.H header file is normally used to define the user-defined embedded boundary class, if any. It is also used to define the solution initialisation and boundary condition implementation functions.

* **prob_parm.H** -- Used to define problem specific parameters and their values.

* **GNUMakefile** -- In addition to setting options to build profiling, debugging, MPI, OpenMP, Compiler toolchain options, the chemical mechanism, transport model, equation of state model, and use of EB are set here for compile time selection. The GNUMakefile includes the ``Make.PeleC`` file from the `Exec` directory that contains build configuration common across the examples.


.. toctree::
   :maxdepth: 1

   Building
   InputFiles
   Tutorials
   Testing
