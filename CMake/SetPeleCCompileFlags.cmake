# Logic for handling warnings
if(PELEC_ENABLE_ALL_WARNINGS)
  # GCC, Clang, and Intel seem to accept these
  list(APPEND PELEC_CXX_FLAGS "-Wall" "-Wextra" "-pedantic")
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    list(APPEND PELEC_CXX_FLAGS "-diag-disable:11074,11076,10397,15335")
  endif()
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 7.0)
    # Avoid notes about -faligned-new with GCC > 7
    list(APPEND PELEC_CXX_FLAGS "-faligned-new")
  endif()
endif()

# Add our extra flags according to language
separate_arguments(PELEC_CXX_FLAGS)
target_compile_options(${pelec_exe_name} PRIVATE $<$<COMPILE_LANGUAGE:CXX>:${PELEC_CXX_FLAGS}>)

if(PELEC_ENABLE_CUDA)
  list(APPEND PELEC_CUDA_FLAGS "--expt-relaxed-constexpr")
  list(APPEND PELEC_CUDA_FLAGS "--expt-extended-lambda")
  list(APPEND PELEC_CUDA_FLAGS "--Wno-deprecated-gpu-targets")
  list(APPEND PELEC_CUDA_FLAGS "-m64")
  if(ENABLE_CUDA_FASTMATH)
    list(APPEND PELEC_CUDA_FLAGS "--use_fast_math")
  endif()
  separate_arguments(PELEC_CUDA_FLAGS)
  target_compile_options(${pelec_exe_name} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:${PELEC_CUDA_FLAGS}>)
  # Add arch flags to both compile and linker to avoid warnings about missing arch
  set(CMAKE_CUDA_FLAGS ${NVCC_ARCH_FLAGS})
  set_target_properties(
    ${pelec_exe_name} PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_RESOLVE_DEVICE_SYMBOLS OFF)
endif()
