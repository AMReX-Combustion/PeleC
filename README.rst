PeleC 
-----
*A compressible AMR combustion code*

`PeleC` is an adaptive-mesh compressible hydrodynamics code for reacting
flows.

`Documentation <https://amrex-combustion.github.io/PeleC/>`_ | `Nightly Test Results <https://my.cdash.org/index.php?project=PeleC>`_

A Note on the Current Status of PeleC Development
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Please be aware that kernels previously written in Fortran in PeleC have been recently ported to C++ to allow for execution on GPUs. These new kernels are located in directories marked ``*Cpp``. The Fortran kernels are currently preserved in this repo using the original directory organization, as all functionality has not yet been ported to C++ and GPUs. However, development in Fortran will not continue and any future development should be done entirely in C++ and with the ability to execute the code within AMReX's GPU framework. The Fortran kernels will be phased out at a future date. This is done in order to run portably using CPUS, CUDA, HIP, and DPC++ programming model backends in the future.

Getting Started
~~~~~~~~~~~~~~~

* To compile and run the `Pele` suite of codes, one needs a C++ compiler that supports the C++14 standard.  A hierarchical strategy for parallelism is supported, based MPI + OpenMP, or MPI + CUDA.  The codes work with all major MPI and OpenMP implementations.  The codes should build and run with no modifications to the `make` system if using a Linux system with the GNU compilers, version 4.9.4 and above.

To build `PeleC` and run a sample 3D flame problem:

1. One can have PeleC use the default submodules for AMReX and PelePhysics in its own repo by simply performing: ::

    git clone --recursive git@github.com:AMReX-Combustion/PeleC.git
    cd PeleC/ExecCpp/RegTests/PMF
    make
    ./Pele3d.xxx.yyy.ex inputs_ex

Alternatively, one can set environment variables to use AMReX and PelePhysics repos from external locations:

1. Set the environment variable, AMREX_HOME, and clone a copy of `AMReX` there: ::

    export AMREX_HOME=<location for AMReX>    
    git clone git@github.com:AMReX-Codes/amrex.git ${AMREX_HOME}

2. Set the environment variable, PELE_PHYSICS_HOME, and clone a copy of `PelePhysics` there. You should be placed in the `development` branch: ::

    export PELE_PHYSICS_HOME=<location for PelePhysics>
    git clone git@github.com:AMReX-Combustion/PelePhysics.git ${PELE_PHYSICS_HOME}

3. Set the environment variable, PELEC_HOME, and clone a copy of `PeleC` there. You should be placed in the `development` branch: ::

    export PELEC_HOME=<location for PeleC>
    git clone git@github.com:AMReX-Combustion/PeleC.git ${PELEC_HOME}

4. Move to an example build folder, build an executable, run a test case: ::

    cd ${PELEC_HOME}/ExecCpp/RegTests/PMF
    make
    ./Pele3d.xxx.yyy.ex inputs_ex

* Notes

   A. In the exec line above, xxx.yyy is a tag identifying your compiler and various build options, and will vary across pltaform.  (Note that GNU compilers must be at least 4.8.4, and MPI should be at least version 3).
   B. The example is 3D premixed flame, flowing vertically upward through the domain with no gravity. The lateral boundaries are periodic.  A detailed hydrogen model is used.  The solution is initialized with a wrinkled (perturbed) 2D steady flame solution computed using the PREMIX code.  Two levels of solution-adaptive refinement are automatically triggered by the presence of the flame intermediate, HO2.
   C. In addition to informative output to the terminal, periodic plotfiles are written in the run folder.  These may be viewed with CCSE's Amrvis (<https://ccse.lbl.gov/Downloads/downloadAmrvis.html>) or Vis-It (<http://vis.lbl.gov/NERSC/Software/visit/>):

      1. In VisIt, direct the File->Open dialogue to select the file named "Header" that is inside each plotfile folder..
      2. With Amrvis, "amrvis3d plt00030", for example.


Dependencies
~~~~~~~~~~~~

`PeleC` was created as a renamed, stripped down version of `Maui`, and is built on the `AMReX` library.  In the process, the Microphysics folder was extracted, and reorganized into a separate repository, `PelePhysics`.  


Development model
~~~~~~~~~~~~~~~~~

To add a new feature to PeleC, the procedure is:

1. Create a branch for the new feature (locally): ::

    git checkout -b AmazingNewFeature

2. Develop the feature, merging changes often from the development branch into your AmazingNewFeature branch: ::
   
    git commit -m "Developed AmazingNewFeature"
    git checkout development
    git pull                     [fix any identified conflicts between local and remote branches of "development"]
    git checkout AmazingNewFeature
    git merge development        [fix any identified conflicts between "development" and "AmazingNewFeature"]

3a. Build and run the full test suite using CMake and CTest (See the `Build` directory for an example script). Please do not introduce warnings. PeleC is checked against `clang-tidy` and `cppcheck` in the CI. To use `cppcheck` and `clang-tidy` locally use these CMake options: ::

   -DPELEC_ENABLE_CLANG_TIDY:BOOL=ON
   -DPELEC_ENABLE_CPPCHECK:BOOL=ON

3b. Run `clang-tidy` by using an LLVM compiler and making sure `clang-tidy` is found during configure. Then `make` will run `clang-tidy` along with compilation. Once verifying `cppcheck` was found during configure, using the `make cppcheck` target should run its checks on the `compile_commands.json` database generated by CMake. More information on these checks can be seen in the CI files used for GitHub Actions in the `.github/workflows` directory.

3c. To easily format all source files before commit, use the following command: ::

    find SourceCpp ExecCpp \( -name "*.cpp" -o -name "*.H" \) -exec clang-format -i {} +

4. Push feature branch to PeleC repository: ::

    git push -u origin AmazingNewFeature [Note: -u option required only for the first push of new branch]

5. Submit a pull request through git@github.com:AMReX-Combustion/PeleC.git, and make sure you are requesting a merge against the development branch

6. Check the CI status on Github and make sure the tests passed for merge request

.. note::

   Github CI uses the CMake build system and CTest to test the core source files of PeleC. If you are adding source files, you will need to add them to the list of source files in the ``CMake`` directory for the tests to pass. Make sure to add them to the GNU make makefiles as well.


Test Status
~~~~~~~~~~~

Nightly test results for PeleC against multiple compilers and machines can be seen on its `CDash page <https://my.cdash.org/index.php?project=PeleC>`_.

Documentation
~~~~~~~~~~~~~

The full documentation for Pele exists in the Docs directory; at present this is maintained inline using
Sphinx  `Sphinx <http://www.sphinx-doc.org>`_. With 
Sphinx, documentation is written in *Restructured Text*. reST is a markup language
similar to Markdown, but with somewhat greater capabilities (and idiosyncrasies). There
are several `primers <http://thomas-cokelaer.info/tutorials/sphinx/rest_syntax.html>`_
available to get started. One gotcha is that indentation matters.

    cd Docs && mkdir build && cd build && cmake .. && make

Acknowledgment
~~~~~~~~~~~~~~

This research was supported by the Exascale Computing Project (ECP), Project
Number: 17-SC-20-SC, a collaborative effort of two DOE organizations -- the
Office of Science and the National Nuclear Security Administration --
responsible for the planning and preparation of a capable exascale ecosystem --
including software, applications, hardware, advanced system engineering, and
early testbed platforms -- to support the nation's exascale computing
imperative.
