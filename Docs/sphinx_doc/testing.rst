.. _Testing:

Testing and Verification
------------------------

Testing and verfication of PeleC can be performed using CTest, which is included in the CMake build system. If one builds PeleC with CMake, the testing suite, and the verification suite, can be enabled during the CMake configure step.

An example ``cmake`` configure command performed in the ``Build`` directory in PeleC is shown below with options relevant to the testing suite:

::

  cmake -DCMAKE_INSTALL_PREFIX:PATH=./install \
        -DCMAKE_BUILD_TYPE:STRING=Release \
        -DPELEC_ENABLE_MPI:BOOL=ON \
        -DCMAKE_CXX_COMPILER:STRING=mpicxx \
        -DCMAKE_C_COMPILER:STRING=mpicc \
        -DCMAKE_Fortran_COMPILER:STRING=mpifort \
        -DENABLE_FCOMPARE:BOOL=ON \
        -DENABLE_TESTS:BOOL=ON \
        -DTEST_WITH_FCOMPARE:BOOL=OFF \
        -DTEST_WITH_FEXTREMA:BOOL=OFF \
        -DENABLE_VERIFICATION:BOOL=ON \
        -DPELEC_ENABLE_MASA:BOOL=ON \
        -DMASA_DIR:STRING=/path/to/masa/dir \
        ..

While performing a ``cmake -LAH ..`` command will give descriptions of every option for the CMake project. Descriptions of particular options regarding the testing suite are listed below:

**ENABLE_FCOMPARE** -- builds the ``fcompare`` utility from AMReX as well as the executable(s), to allow for testing differences between plot files

**ENABLE_FEXTREMA** -- builds the ``fextrema`` utility from AMReX as well as the executable(s), to allow for testing differences between plot files

**ENABLE_TESTS** -- enables the base level regression test suite that will check whether each test will run its executable to completion successfully

**TEST_WITH_FCOMPARE** -- enables an additional step in the regression tests where the ``fcompare`` program from AMReX will also test for differences in the plots generated from the tests against "gold" files which contain previously verified results to machine precision

**TEST_WITH_FEXTREMA** -- enables an additional step in the regression tests where the ``fextrema`` program from AMReX will also test for differences in the maxima and minima of the plots generated from the tests against "gold" files which contain previously verified results to machine precision (this is the default method in which Travis CI detects diffs in results because it's more portable than storing plot files)

**ENABLE_VERIFICATION** -- enables the verification suite which checks that PeleC is second order accurate using several additional tests, but note that the verification tests can take a significant amount of time to run and certain Python modules are expected to exist on the user's system to generate PNG plot files

**PELEC_ENABLE_MASA** and **MASA_DIR** -- are required when the verification suite is enabled to perform the method of manufactured solutions


Building the Tests
~~~~~~~~~~~~~~~~~~

Once the user has performed the CMake configure step, the ``make`` command will build every executable required for each test. In this step, it is highly beneficial for the user to use the ``-j`` option for ``make`` to build source files in parallel. It is also beneficial if the user has access to the Ninja build system with Fortran support which is described in the section on building PeleC.

Running the Tests
~~~~~~~~~~~~~~~~~

Once the test executables are built, CTest also creates working directories for each test within the ``Build`` directory where plot files will be output, etc. This directory is analogous to the source location of the tests in ``Testing/test_files``.

To run the test suite, run ``ctest`` in the ``Build`` directory. CTest will run the tests and report their exit status. Useful options for CTest are ``-VV`` which runs in a verbose mode where the output of each test can be seen. ``-R`` where a regex string can be used to run specific sets of tests. ``-j`` where CTest will bin pack and run tests in parallel based on how many processes each test is specified to use and fit them into the amount of cores available on the machine. ``-L`` where the subset of tests containing a particular label will be run. Output for the last set of tests run is available in the ``Build`` directory in ``Testing/Temporary/LastTest.log``.

Adding Tests
~~~~~~~~~~~~

Developers are encouraged to add tests to PeleC and in this section we describe how the tests are organized in the CTest framework. The locations of the tests are in ``PeleC/Testing``. To add a test, first create a test directory with a name in ``PeleC/Testing/test_files/<test_name>``. Place an ``exe_options.cmake`` file in the directory, using others tests as an example. Place the input file for the test as ``PeleC/Testing/test_files/<test_name>.i`` and a "probin" file as ``PeleC/Testing/test_files/<test_name>.probin``, along with any other files necessary for the test. Any file in the test directory will be copied during CMake configure to the test's working directory. Next, edit the ``PeleC/Testing/CTestList.cmake`` file, add the test to the list, and specify the number of MPI ranks the test should be run with. Note there are different categories of tests and if your test falls outside of these categories, a new function to add the test will need to be created. After these steps, your test will be automatically added to the test suite database when doing the CMake configure with the testing suite enabled.
