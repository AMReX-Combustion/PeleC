.. highlight:: rst

Using KGen
----------

`KGen <https://github.com/NCAR/KGen>`__ is a tool developed by NCAR for extracting sections of Fortran code into a standalone program. This is useful for narrowing the scope of code involved in developing and porting computational kernels for new computer architectures or any other task that might benefit from this agility.

The way KGen works is by the developer first adding directives to the original source code files for whatever section or routine the they would like to extract. Next the developer can choose a set of invocations for the code section in which it will record input and output data to disk. Then by running KGen, it will modify the original code to capture the input and output variables (or "state") for the code section, then build and run the full application and write the state data to disk. Next it will extract the kernel into a standalone program that is hardcoded to read its input from the state files. It uses this state data to both run and validate the resulting standalone kernel output against the original application output. It also measures the elasped time for each call to the kernel since typically the goal is to reduce runtime of the kernel.

To describe this process further, we give the following example where we will use KGen to extract the 3D diffterm Fortran kernel in PeleC. KGen generally requires Linux due to a dependency on the strace tool which doesn't exist on other operating systems such as macOS. We will assume you have all the necessary repos for PeleC cloned in ${HOME}/combustion, for example PeleC cloned as ${HOME}/combustion/PeleC, etc.

We then perform the following steps:

1. Run something like the following commands: `cd ${HOME} && mkdir diffterm && cd diffterm && git clone https://github.com/NCAR/KGen.git && mkdir diffterm_kernel`
2. Create a Makefile for KGen in ${HOME}/diffterm/diffterm_kernel/ which looks similar to this example:

::

    KGEN_HOME := ${HOME}/diffterm/KGen
    KGEN := ${KGEN_HOME}/bin/kgen
    SRC_DIR := ${HOME}/combustion/PeleC/Source
    EXEC_DIR := ${HOME}/combustion/PeleC/Exec/PMF
    SRC := ${SRC_DIR}/Src_3d/diffterm_3d.f90
    
    test:
    	${KGEN} \
    		--cmd-clean "cd ${EXEC_DIR}; make realclean" \
    		--cmd-build "cd ${EXEC_DIR}; make -j8" \
    		--cmd-run "cd ${EXEC_DIR}; ./PeleC3d.gnu.ex inputs-3d-regt" \
    		--invocation 0:0:1,0:0:2,0:0:3,0:0:4,0:0:5,0:0:6,0:0:7,0:0:119,0:0:120,0:0:121,0:0:346,0:0:347,0:0:348,0:0:349,0:0:350,0:0:697,0:0:698,0:0:699,0:0:700,0:0:701 \
    		${SRC}
    
    clean:
    	${MAKE} clean -C src
    	rm -rf kernel state kgen.log strace.log include.ini _kgen_compflag_cmdwrapper.sh model model.ini elapsedtime coverage papi

3. Next we edit the `${HOME}/combustion/PeleC/Source/Src_3d/diffterm_3d.f90` file and add KGen directives around the kernel of interest which should look as such:

::

    diff --git a/Source/Src_3d/diffterm_3d.f90 b/Source/Src_3d/diffterm_3d.f90
    index e3a76c6..b29000b 100644
    --- a/Source/Src_3d/diffterm_3d.f90
    +++ b/Source/Src_3d/diffterm_3d.f90
    @@ -120,6 +120,8 @@ contains
         double precision, parameter :: twoThirds = 2.d0/3.d0
         double precision :: dxinv(3)
     
    +    !$kgen begin_callsite diffterm
    +
         dxinv = 1.d0/deltax
     
         call eos_ytx_vec(Q(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,QFS:QFS+nspecies-1),lo,hi,X,lo,hi,lo,hi,nspecies)
    @@ -392,6 +394,7 @@ contains
               end do
            end do
         end do
    +    !$kgen end_callsite
     
       end subroutine pc_diffterm

4. Now we should be able to just type the `make` command and KGen will instrument the original application, build the application, run the application to collect state data, and then extract the kernel into a `kernel` directory.

5. Note, for the diffterm kernel, there are calls to C routines which KGen is not able to account for so KGen will complain about missing symbols for the standalone kernel. To fix this, for our example we can copy the `${HOME}/combustion/PelePhysics/Support/Fuego/Mechanism/Models/LiDryer/LiDryer.c` and `${HOME}/combustion/PelePhysics/Support/Fuego/Evaluation/ReactionData.H` files from the PelePhysics repo into the `kernel` directory and include rules and the necessary dependencies for compiling into an object file and linking into the kernel executable. An example of the kernel makefile additions is listed below:

::

    $ diff Makefile.orig Makefile
    3a4
    > CC_0 := /nopt/nrel/apps/gcc/6.2.0/bin/gcc
    4a6
    > C_FLAGS_SET_0 := -g -O3 -DBL_FORT_USE_UNDERSCORE
    6c8
    < ALL_OBJS := diffterm_3d.o eos.o AMReX_constants_mod.o AMReX_fort_mod.o meth_params.o actual_network.o prob_params.o kernel_driver.o kgen_utils.o tprof_mod.o
    ---
    > ALL_OBJS := diffterm_3d.o eos.o AMReX_constants_mod.o AMReX_fort_mod.o meth_params.o actual_network.o prob_params.o kernel_driver.o kgen_utils.o tprof_mod.o LiDryer.o
    17c19
    < eos.o: eos.F90 kgen_utils.o tprof_mod.o AMReX_fort_mod.o AMReX_constants_mod.o
    ---
    > eos.o: eos.F90 kgen_utils.o tprof_mod.o AMReX_fort_mod.o AMReX_constants_mod.o LiDryer.o
    42a45,47
    > 
    > LiDryer.o: LiDryer.c
    > 	${CC_0} ${C_FLAGS_SET_0} -I . -c -o $@ $<

6. To run the kernel we `cd` to the kernel directory and type `make` which will build the kernel and run the kernel. Now we can refactor the kernel code and run it to validate the results see the effects on the runtime. When we are satisfied, it should be trivial to copy and paste the refactored kernel into the original application.
