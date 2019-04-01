#!/bin/bash -l

# Script for running nightly regression tests for PeleC on a particular set 
# of machines with a list of configurations for each machine using Spack
# to satisfy dependencies and submitting results to CDash

# Control over printing and executing commands
print_cmds=true
execute_cmds=true

# Function for printing and executing commands
cmd() {
  if ${print_cmds}; then echo "+ $@"; fi
  if ${execute_cmds}; then eval "$@"; fi
}

# Function for testing a single configuration
test_configuration() {
  COMPILER_ID="${COMPILER_NAME}@${COMPILER_VERSION}"
  printf "************************************************************\n"
  printf "Testing PeleC with:\n"
  printf "${COMPILER_ID}\n"
  printf "MPI_ENABLED: ${MPI_ENABLED}\n"
  printf "OPENMP_ENABLED: ${OPENMP_ENABLED}\n"
  printf "LIST_OF_TPLS: ${LIST_OF_TPLS}\n"
  printf "at $(date)\n"
  printf "************************************************************\n"
  printf "\n"

  # Logic for building up some constraints for use on Spack commands
  GENERAL_CONSTRAINTS=''
  MPI_ID=''
  MPI_CONSTRAINTS=''
  BLAS_ID=''
  BLAS_CONSTRAINTS=''
  if [ "${COMPILER_NAME}" == 'gcc' ] || [ "${COMPILER_NAME}" == 'clang' ]; then
    # OpenMPI 3.1.3 hangs at run time unless it was built with GCC > 7.3.0
    # so we use an older OpenMPI for GCC 4.9.4.
    MPI_ID="openmpi"
    if [ "${COMPILER_VERSION}" == '4.9.4' ]; then
      MPI_ID="openmpi@1.10.7"
    fi
    if [ "${MACHINE_NAME}" == 'eagle' ]; then
      MPI_ID="openmpi@3.1.3"
    fi
  elif [ "${COMPILER_NAME}" == 'intel' ]; then
    # For intel, we want to build against intel-mpi and intel-mkl
    MPI_ID="intel-mpi"
    BLAS_ID="intel-mkl"
  fi
  #if [ ! -z "${MPI_ID}" ]; then
  #  # Avoid listing plain openmpi without a version number
  #  if [ "${MPI_ID}" == 'openmpi' ]; then
  #    MPI_CONSTRAINTS=''
  #  else
  #    MPI_CONSTRAINTS="^${MPI_ID}"
  #  fi
  #fi
  #if [ ! -z "${BLAS_ID}" ]; then
  #  BLAS_CONSTRAINTS=" ^${BLAS_ID}"
  #fi
  #GENERAL_CONSTRAINTS="${MPI_CONSTRAINTS}${BLAS_CONSTRAINTS}"
  #printf "Using constraints: ${GENERAL_CONSTRAINTS}\n\n"

  cmd "cd ${PELEC_TESTING_ROOT_DIR}"

  printf "\nLoading modules...\n"
  if [ "${MACHINE_NAME}" == 'rhodes' ]; then
    cmd "module purge"
    cmd "module unuse ${MODULEPATH}"
    cmd "module use /opt/compilers/modules"
    cmd "module use /opt/utilities/modules"
    cmd "module load unzip"
    cmd "module load patch"
    cmd "module load bzip2"
    cmd "module load git"
    cmd "module load flex"
    cmd "module load bison"
    cmd "module load wget"
    cmd "module load bc"
    cmd "module load python/2.7.15"
    cmd "module load cppcheck"
    cmd "module load binutils"
    if [ "${COMPILER_NAME}" == 'gcc' ]; then
      cmd "module load ${COMPILER_NAME}/${COMPILER_VERSION}"
    elif [ "${COMPILER_NAME}" == 'clang' ]; then
      cmd "module load llvm/${COMPILER_VERSION}"
    elif [ "${COMPILER_NAME}" == 'intel' ]; then
      cmd "module load ${INTEL_COMPILER_MODULE}"
    fi
  elif [ "${MACHINE_NAME}" == 'peregrine' ] || [ "${MACHINE_NAME}" == 'eagle' ]; then
    cmd "module purge"
    cmd "module unuse ${MODULEPATH}"
    cmd "module use /nopt/nrel/ecom/hpacf/compilers/modules"
    cmd "module use /nopt/nrel/ecom/hpacf/utilities/modules"
    cmd "module load python/2.7.15"
    cmd "module load git"
    cmd "module load cppcheck"
    cmd "module load binutils"
    if [ "${COMPILER_NAME}" == 'gcc' ]; then
      cmd "module load ${COMPILER_NAME}/${COMPILER_VERSION}"
    elif [ "${COMPILER_NAME}" == 'intel' ]; then
      cmd "module load ${INTEL_COMPILER_MODULE}"
    fi
  fi

  # Set the TMPDIR to disk so it doesn't run out of space
  if [ "${MACHINE_NAME}" == 'peregrine' ] || [ "${MACHINE_NAME}" == 'eagle' ]; then
    printf "\nMaking and setting TMPDIR to disk...\n"
    cmd "mkdir -p /scratch/${USER}/.tmp"
    cmd "export TMPDIR=/scratch/${USER}/.tmp"
  fi

  # Uninstall packages we want to track; it's an error if they don't exist yet, but a soft error
  #printf "\nUninstalling MASA (this is fine to error when tests are first run or building MASA has previously failed)...\n"
  #cmd "spack uninstall -a -y masa %${COMPILER_ID} || true"

  # Update packages we want to track; it's an error if they don't exist yet, but a soft error
  #printf "\nUpdating MASA (this is fine to error when tests are first run)...\n"
  #cmd "spack cd masa %${COMPILER_ID} ${GENERAL_CONSTRAINTS} && pwd && git fetch --all && git reset --hard origin/master && git clean -df && git status -uno || true"

  cmd "cd ${PELEC_TESTING_ROOT_DIR}" # Change directories to avoid any stale file handles

  TPL_VARIANTS=''
  TPLS=(${LIST_OF_TPLS//;/ })
  for TPL in ${TPLS[*]}; do
    TPL_VARIANTS+="+${TPL}"
  done

  if [ "${MACHINE_NAME}" != 'mac' ]; then
    cmd "module list"
  fi

  printf "\nInstalling PeleC dependencies using ${COMPILER_ID}...\n"
  #cmd "spack install --only dependencies pelec ${TPL_VARIANTS} %${COMPILER_ID} ${GENERAL_CONSTRAINTS}"
  cmd "spack install masa %${COMPILER_ID} ${GENERAL_CONSTRAINTS}"
  cmd "spack install cmake %${COMPILER_ID} ${GENERAL_CONSTRAINTS}"
  if [ "${COMPILER_ID}" == 'gcc@4.9.4' ]; then
    cmd "spack install openmpi@1.10.7 %${COMPILER_ID} ${GENERAL_CONSTRAINTS}"
  else
    cmd "spack install openmpi %${COMPILER_ID} ${GENERAL_CONSTRAINTS}"
  fi

  #STAGE_DIR=$(spack location -S)
  #if [ ! -z "${STAGE_DIR}" ]; then
  #  #Haven't been able to find another robust way to rm with exclude
  #  printf "\nRemoving all staged directories except Trilinos...\n"
  #  cmd "cd ${STAGE_DIR} && rm -rf a* b* c* d* e* f* g* h* i* j* k* l* m* n* o* p* q* r* s* tar* ti* u* v* w* x* y* z*"
  #  #printf "\nRemoving all staged directories except Trilinos and OpenFAST...\n"
  #  #cmd "cd ${STAGE_DIR} && rm -rf a* b* c* d* e* f* g* h* i* j* k* l* m* n* openmpi* p* q* r* s* tar* u* v* w* x* y* z*"
  #  #find ${STAGE_DIR}/ -maxdepth 0 -type d -not -name "trilinos*" -exec rm -r {} \;
  #fi

  # Refresh available modules (this is only really necessary on the first run of this script
  # because cmake and openmpi will already have been built and module files registered in subsequent runs)
  cmd "source ${SPACK_ROOT}/share/spack/setup-env.sh"

  printf "\nLoading Spack modules into environment for CMake and MPI to use during CTest...\n"
  if [ "${MACHINE_NAME}" == 'mac' ]; then
    cmd "export PATH=$(spack location -i cmake %${COMPILER_ID})/bin:${PATH}"
    cmd "export PATH=$(spack location -i ${MPI_ID} %${COMPILER_ID})/bin:${PATH}"
  else
    cmd "spack load cmake %${COMPILER_ID}"
    cmd "spack load ${MPI_ID} %${COMPILER_ID}"
  fi

  printf "\nSetting variables to pass to CTest...\n"
  CMAKE_CONFIGURE_ARGS=''
  for TPL in ${TPLS[*]}; do
    if [ "${TPL}" == 'masa' ]; then
      MASA_DIR=$(spack location -i masa %${COMPILER_ID})
      CMAKE_CONFIGURE_ARGS="-DPELEC_ENABLE_MASA:BOOL=ON -DMASA_DIR:PATH=${MASA_DIR} ${CMAKE_CONFIGURE_ARGS}"
      printf "MASA_DIR=${MASA_DIR}\n"
    fi
  done

  # Set the extra identifiers for CDash build description
  EXTRA_BUILD_NAME="-${COMPILER_NAME}-${COMPILER_VERSION}"

  if [ ! -z "${PELEC_DIR}" ]; then
    printf "\nCleaning PeleC directory...\n"
    cmd "cd ${PELEC_DIR} && git reset --hard origin/ctest && git clean -df && git status -uno"
    cmd "cd ${PELEC_DIR}/build && rm -rf ${PELEC_DIR}/build/*"
    # Update all the submodules recursively in case the previous ctest update failed because of submodule updates
    cmd "cd ${PELEC_DIR} && git submodule update --init --recursive"
    cmd "ln -s ${HOME}/combustion/PeleCGoldFiles ${PELEC_DIR}/Testing/PeleCGoldFiles"
  fi

  #if [ "${OPENMP_ENABLED}" == 'true' ]; then
  #  printf "\nSetting OpenMP stuff...\n"
  #  cmd "export OMP_NUM_THREADS=1"
  #  cmd "export OMP_PROC_BIND=false"
  #fi

  # Run static analysis and let ctest know we have static analysis output
  #if [ "${MACHINE_NAME}" == 'peregrine' ] || \
  #   [ "${MACHINE_NAME}" == 'mac' ] || \
  #   [ "${MACHINE_NAME}" == 'rhodes' ]; then
  #  printf "\nRunning cppcheck static analysis (PeleC not updated until after this step)...\n"
  #  cmd "rm ${LOGS_DIR}/pelec-static-analysis.txt"
  #  cmd "cppcheck --enable=all --quiet -j 8 --output-file=${LOGS_DIR}/pelec-static-analysis.txt -I ${PELEC_DIR}/include ${PELEC_DIR}/src"
  #  cmd "printf \"%s warnings\n\" \"$(wc -l < ${LOGS_DIR}/pelec-static-analysis.txt | xargs echo -n)\" >> ${LOGS_DIR}/pelec-static-analysis.txt"
  #  CTEST_ARGS="-DHAVE_STATIC_ANALYSIS_OUTPUT:BOOL=TRUE -DSTATIC_ANALYSIS_LOG=${LOGS_DIR}/pelec-static-analysis.txt ${CTEST_ARGS}"
  #fi

  # Unset the TMPDIR variable after building but before testing during ctest nightly script
  if [ "${MACHINE_NAME}" == 'peregrine' ] || [ "${MACHINE_NAME}" == 'eagle' ]; then
    CTEST_ARGS="-DUNSET_TMPDIR_VAR:BOOL=TRUE ${CTEST_ARGS}"
  fi

  # Turn on all warnings unless we're gcc 4.9.4
  #if [ "${COMPILER_ID}" == 'gcc@4.9.4' ]; then
  #  CMAKE_CONFIGURE_ARGS="-DENABLE_ALL_WARNINGS:BOOL=FALSE ${CMAKE_CONFIGURE_ARGS}"
  #else
  #  CMAKE_CONFIGURE_ARGS="-DENABLE_ALL_WARNINGS:BOOL=TRUE ${CMAKE_CONFIGURE_ARGS}"
  #fi

  # Turn on address sanitizer for clang build on rhodes
  if [ "${COMPILER_NAME}" == 'clang' ] && [ "${MACHINE_NAME}" == 'rhodes' ]; then
    printf "\nSetting up address sanitizer in Clang...\n"
    export CXXFLAGS="-fsanitize=address -fno-omit-frame-pointer"
    printf "export CXXFLAGS=${CXX_FLAGS}\n"
    cmd "export ASAN_OPTIONS=detect_container_overflow=0"
    printf "Writing asan.supp file...\n"
    (set -x; printf "leak:libopen-pal\nleak:libmpi" > ${PELEC_DIR}/build/asan.supp)
    cmd "export LSAN_OPTIONS=suppressions=${PELEC_DIR}/build/asan.supp"
    #CMAKE_CONFIGURE_ARGS="-DCMAKE_CXX_FLAGS:STRING=-fsanitize=address\ -fno-omit-frame-pointer ${CMAKE_CONFIGURE_ARGS}"
    #CMAKE_CONFIGURE_ARGS="-DCMAKE_LINKER=clang++ -DCMAKE_CXX_LINK_EXECUTABLE=clang++ -DCMAKE_CXX_FLAGS:STRING=\'-fsanitize=address -fno-omit-frame-pointer\' -DCMAKE_EXE_LINKER_FLAGS:STRING=-fsanitize=address ${CMAKE_CONFIGURE_ARGS}"
    #printf "Disabling OpenMP in PeleC for address sanitizer...\n"
    #CMAKE_CONFIGURE_ARGS="-DENABLE_OPENMP:BOOL=FALSE ${CMAKE_CONFIGURE_ARGS}"
  fi

  # Explicitly set compilers to MPI compilers
  if [ "${COMPILER_NAME}" == 'gcc' ] || [ "${COMPILER_NAME}" == 'clang' ]; then
    MPI_CXX_COMPILER=mpicxx
    MPI_C_COMPILER=mpicc
    MPI_FORTRAN_COMPILER=mpifort
  elif [ "${COMPILER_NAME}" == 'intel' ]; then
    MPI_CXX_COMPILER=mpiicpc
    MPI_C_COMPILER=mpiicc
    MPI_FORTRAN_COMPILER=mpiifort
  fi

  printf "\nListing cmake and compilers that will be used in ctest...\n"
  cmd "which ${MPI_CXX_COMPILER}"
  cmd "which ${MPI_C_COMPILER}"
  cmd "which ${MPI_FORTRAN_COMPILER}"
  cmd "which mpiexec"
  cmd "which cmake"

  CMAKE_CONFIGURE_ARGS="-DCMAKE_CXX_COMPILER:STRING=${MPI_CXX_COMPILER} -DCMAKE_C_COMPILER:STRING=${MPI_C_COMPILER} -DCMAKE_Fortran_COMPILER:STRING=${MPI_FORTRAN_COMPILER} -DMPI_CXX_COMPILER:STRING=${MPI_CXX_COMPILER} -DMPI_C_COMPILER:STRING=${MPI_C_COMPILER} -DMPI_Fortran_COMPILER:STRING=${MPI_FORTRAN_COMPILER} ${CMAKE_CONFIGURE_ARGS}"

  # Set essential arguments for ctest
  CTEST_ARGS="-DTESTING_ROOT_DIR=${PELEC_TESTING_ROOT_DIR} -DPELEC_DIR=${PELEC_TESTING_ROOT_DIR}/pelec -DTEST_LOG=${LOGS_DIR}/pelec-test-log.txt -DHOST_NAME=${HOST_NAME} -DEXTRA_BUILD_NAME=${EXTRA_BUILD_NAME} ${CTEST_ARGS}"

  # Set essential arguments for the ctest cmake configure step
  CMAKE_CONFIGURE_ARGS="-DCMAKE_BUILD_TYPE=RelWithDebInfo ${CMAKE_CONFIGURE_ARGS}"

  # Set looser diff tolerance for GCC 7.3.0 cases that have more optimization flags on
  #if [ "${COMPILER_ID}" == 'gcc@7.3.0' ] && [ "${MACHINE_NAME}" != 'mac' ]; then
  #  CMAKE_CONFIGURE_ARGS="-DTEST_TOLERANCE:STRING=0.00001 ${CMAKE_CONFIGURE_ARGS}"
  #fi

  # Allow OpenMPI to consider hardware threads as cpus and allow for oversubscription
  if [ "${COMPILER_NAME}" != 'intel' ]; then
    CMAKE_CONFIGURE_ARGS="-DMPIEXEC_PREFLAGS:STRING=--oversubscribe ${CMAKE_CONFIGURE_ARGS}"
  fi

  printf "\nRunning CTest at $(date)...\n"
  cmd "cd ${PELEC_DIR}/build"
  if [ "${MACHINE_NAME}" != 'mac' ]; then
    cmd "module list"
  fi
  cmd "ctest ${CTEST_ARGS} -DCMAKE_CONFIGURE_ARGS=\"${CMAKE_CONFIGURE_ARGS}\" -VV -S ${PELEC_DIR}/Testing/CTestNightlyScript.cmake"
  printf "Returned from CTest at $(date)\n"

  printf "\nSaving golds...\n"
  (set -x; find ${PELEC_DIR}/build/Testing/test_files -type d -name *plt00010* | tar -czf ${GOLDS_DIR}/golds${EXTRA_BUILD_NAME}-$(date +%Y-%m-%d-%H-%M).tar.gz -T -)

  printf "\n"
  printf "************************************************************\n"
  printf "Done testing PeleC with:\n"
  printf "${COMPILER_ID}\n"
  printf "MPI_ENABLED: ${MPI_ENABLED}\n"
  printf "OPENMP_ENABLED: ${OPENMP_ENABLED}\n"
  printf "LIST_OF_TPLS: ${LIST_OF_TPLS}\n"
  printf "at $(date)\n"
  printf "************************************************************\n"
}

# Main function for assembling configurations to test
main() {
  printf "============================================================\n"
  printf "$(date)\n"
  printf "============================================================\n"
  printf "Job is running on ${HOSTNAME}\n"
  printf "============================================================\n"

  # Decide what machine we are on
  if [ "${NREL_CLUSTER}" == 'peregrine' ]; then
    MACHINE_NAME=peregrine
  elif [ "${NREL_CLUSTER}" == 'eagle' ]; then
    MACHINE_NAME=eagle
  fi
  if [ $(hostname) == 'rhodes.hpc.nrel.gov' ]; then
    MACHINE_NAME=rhodes
  elif [ $(hostname) == 'jrood-31712s.nrel.gov' ]; then
    MACHINE_NAME=mac
  fi
    
  HOST_NAME="${MACHINE_NAME}.hpc.nrel.gov"
 
  # Set configurations to test for each machine
  declare -a CONFIGURATIONS
  #CONFIGURATION[n]='compiler_name:compiler_version:mpi_enabled:openmp_enabled:list_of_tpls'
  if [ "${MACHINE_NAME}" == 'rhodes' ]; then
    CONFIGURATIONS[0]='gcc:7.3.0:true:false:masa'
    CONFIGURATIONS[1]='gcc:4.9.4:true:false:masa'
    CONFIGURATIONS[2]='intel:18.0.4:true:false:masa'
    CONFIGURATIONS[3]='clang:6.0.1:true:false:masa'
    PELEC_TESTING_ROOT_DIR=/projects/ecp/combustion/pelec-testing2
    INTEL_COMPILER_MODULE=intel-parallel-studio/cluster.2018.4
  elif [ "${MACHINE_NAME}" == 'peregrine' ]; then
    CONFIGURATIONS[0]='gcc:7.3.0:true:false:masa'
    PELEC_TESTING_ROOT_DIR=/projects/ExaCT/pelec-testing2
    INTEL_COMPILER_MODULE=intel-parallel-studio/cluster.2018.4
  elif [ "${MACHINE_NAME}" == 'eagle' ]; then
    CONFIGURATIONS[0]='gcc:7.3.0:true:false:masa'
    PELEC_TESTING_ROOT_DIR=/projects/ExaCT/pelec-testing2
    INTEL_COMPILER_MODULE=intel-parallel-studio/cluster.2018.4
  elif [ "${MACHINE_NAME}" == 'mac' ]; then
    CONFIGURATIONS[0]='gcc:7.3.0:true:false:masa'
    CONFIGURATIONS[1]='intel:18.0.4:true:false:masa'
    PELEC_TESTING_ROOT_DIR=${HOME}/pelec-testing2
  else
    printf "\nMachine name not recognized.\n"
  fi
 
  PELEC_DIR=${PELEC_TESTING_ROOT_DIR}/pelec
  BUILD_TEST_DIR=${PELEC_TESTING_ROOT_DIR}/build-test
  LOGS_DIR=${PELEC_TESTING_ROOT_DIR}/logs
  GOLDS_DIR=${PELEC_TESTING_ROOT_DIR}/golds
  cmd "export SPACK_ROOT=${PELEC_TESTING_ROOT_DIR}/spack"
 
  printf "============================================================\n"
  printf "HOST_NAME: ${HOST_NAME}\n"
  printf "PELEC_TESTING_ROOT_DIR: ${PELEC_TESTING_ROOT_DIR}\n"
  printf "PELEC_DIR: ${PELEC_DIR}\n"
  printf "BUILD_TEST_DIR: ${BUILD_TEST_DIR}\n"
  printf "LOGS_DIR: ${LOGS_DIR}\n"
  printf "GOLDS_DIR: ${GOLDS_DIR}\n"
  printf "SPACK_ROOT: ${SPACK_ROOT}\n"
  printf "Testing configurations:\n"
  printf " compiler_name:compiler_version:mpi_enabled:openmp_enabled:list_of_tpls\n"
  for CONFIGURATION in "${CONFIGURATIONS[@]}"; do
    printf " ${CONFIGURATION}\n"
  done
  printf "============================================================\n"
 
  if [ ! -d "${PELEC_TESTING_ROOT_DIR}" ]; then
    set -e
    printf "============================================================\n"
    printf "Top level testing directory doesn't exist.\n"
    printf "Creating everything from scratch...\n"
    printf "============================================================\n"

    printf "Creating top level testing directory...\n"
    cmd "mkdir -p ${PELEC_TESTING_ROOT_DIR}"
 
    printf "\nCloning Spack repo...\n"
    cmd "git clone https://github.com/spack/spack.git ${SPACK_ROOT}"
 
    printf "\nConfiguring Spack...\n"
    cmd "git clone https://github.com/exawind/build-test.git ${BUILD_TEST_DIR}"
    cmd "cd ${BUILD_TEST_DIR}/configs && ./setup-spack.sh"
 
    # Checkout PeleC and meshes submodule outside of Spack so ctest can build it itself
    printf "\nCloning PeleC repo...\n"
    cmd "git clone --recursive -b ctest https://github.com/AMReX-Combustion/PeleC.git ${PELEC_DIR}"
    #cmd "mkdir -p ${PELEC_DIR}/build"
 
    printf "\nMaking job output directory...\n"
    cmd "mkdir -p ${LOGS_DIR}"

    printf "\nMaking golds archive directory...\n"
    cmd "mkdir -p ${GOLDS_DIR}"
 
    printf "============================================================\n"
    printf "Done setting up testing directory\n"
    printf "============================================================\n"
    set +e
  fi
 
  printf "\nLoading Spack...\n"
  cmd "source ${SPACK_ROOT}/share/spack/setup-env.sh"

  printf "\n"
  printf "============================================================\n"
  printf "Starting testing loops...\n"
  printf "============================================================\n"
 
  # Test PeleC for the list of configurations
  for CONFIGURATION in "${CONFIGURATIONS[@]}"; do
    CONFIG=(${CONFIGURATION//:/ })
    COMPILER_NAME=${CONFIG[0]}
    COMPILER_VERSION=${CONFIG[1]}
    MPI_ENABLED=${CONFIG[2]}
    OPENMP_ENABLED=${CONFIG[3]}
    LIST_OF_TPLS=${CONFIG[4]}
 
    printf "\nRemoving previous test log for uploading to CDash...\n"
    cmd "rm ${LOGS_DIR}/pelec-test-log.txt"
    (test_configuration) 2>&1 | tee -i ${LOGS_DIR}/pelec-test-log.txt
  done

  printf "============================================================\n"
  printf "Done with testing loops\n"
  printf "============================================================\n"
  printf "============================================================\n"
  printf "Final steps\n"
  printf "============================================================\n"
 
  if [ "${MACHINE_NAME}" == 'peregrine' ] || \
     [ "${MACHINE_NAME}" == 'eagle' ] || \
     [ "${MACHINE_NAME}" == 'rhodes' ]; then
    printf "\nSetting permissions...\n"
    cmd "chmod -R a+rX,go-w ${PELEC_TESTING_ROOT_DIR}"
  fi

  if [ "${MACHINE_NAME}" == 'rhodes' ]; then
    printf "\nSetting group...\n"
    cmd "chgrp -R windsim ${PELEC_TESTING_ROOT_DIR}"
  fi

  printf "============================================================\n"
  printf "Done!\n"
  printf "$(date)\n"
  printf "============================================================\n"
}

main "$@"
