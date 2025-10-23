cmake_minimum_required(VERSION 3.18)
project(OrthrosDLCT LANGUAGES C CXX)

option(ENABLE_CUDA "Build CUDA binaries" ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (ENABLE_CUDA)
    # Check if nvcc is available before trying to enable CUDA
    find_program(NVCC_EXECUTABLE nvcc)
    if (NOT NVCC_EXECUTABLE)
        message(WARNING "nvcc not found in PATH. Disabling CUDA targets.")
        set(ENABLE_CUDA OFF)
    else()
        enable_language(CUDA)
        find_package(CUDAToolkit QUIET)
        if (NOT CUDAToolkit_FOUND)
            message(WARNING "CUDA toolkit not found. Disabling CUDA targets.")
            set(ENABLE_CUDA OFF)
        endif()
    endif()
endif()

if (ENABLE_CUDA)
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)

    # Suppress deprecated GPU architecture warnings
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Wno-deprecated-gpu-targets")

    if (NOT CMAKE_CUDA_ARCHITECTURES)
        set(CMAKE_CUDA_ARCHITECTURES 89)
    endif()
endif()

find_package(OpenMP QUIET)

if (NOT OpenMP_CXX_FOUND)
    if (APPLE AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
        set(_openmp_roots
            /opt/homebrew/opt/libomp
            /usr/local/opt/libomp
            $ENV{OPENMP_ROOT}
        )

        foreach(_root IN LISTS _openmp_roots)
            if (_root)
                find_path(OpenMP_CXX_INCLUDE_DIR
                    NAMES omp.h
                    HINTS "${_root}/include"
                )
                find_library(OpenMP_CXX_LIBRARY
                    NAMES omp
                    HINTS "${_root}/lib"
                )
                if (OpenMP_CXX_INCLUDE_DIR AND OpenMP_CXX_LIBRARY)
                    set(OpenMP_CXX_FLAGS "-Xpreprocessor;-fopenmp")
                    set(OpenMP_CXX_FOUND ON)
                    set(OpenMP_FOUND ON)
                    set(OpenMP_CXX_LIB_NAMES "omp")
                    set(OpenMP_omp_LIBRARY "${OpenMP_CXX_LIBRARY}")
                    set(OpenMP_CXX_INCLUDE_DIRS "${OpenMP_CXX_INCLUDE_DIR}")
                    break()
                endif()
            endif()
        endforeach()
    endif()
endif()

add_library(orthros_analysis STATIC
    src/correlation_analysis.cpp
)

target_include_directories(orthros_analysis PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

if (ENABLE_CUDA)
    add_executable(cuda_dlct
        src/cuda_dlct.cu
    )

    target_include_directories(cuda_dlct PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
    )

    target_link_libraries(cuda_dlct PRIVATE
        CUDA::cudart
        orthros_analysis
    )

    set_target_properties(cuda_dlct PROPERTIES
        CUDA_RUNTIME_LIBRARY Static
        CUDA_RESOLVE_DEVICE_SYMBOLS ON
    )

    # Optimize for register usage and performance
    target_compile_options(cuda_dlct PRIVATE
        $<$<COMPILE_LANGUAGE:CUDA>:-O3>           # Maximum optimization
        $<$<COMPILE_LANGUAGE:CUDA>:--use_fast_math>  # Fast math operations
        $<$<COMPILE_LANGUAGE:CUDA>:--maxrregcount=64> # Limit registers for higher occupancy
        $<$<COMPILE_LANGUAGE:CUDA>:-lineinfo>     # Keep line info for profiling
    )

    install(TARGETS cuda_dlct DESTINATION bin)
endif()

add_executable(openmp_dlct
    src/openmp_dlct.cpp
    src/openmp/random_utils.cpp
    src/openmp/dlct_runner.cpp
)

target_include_directories(openmp_dlct PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_link_libraries(openmp_dlct PRIVATE
    orthros_analysis
)

if (OpenMP_CXX_FOUND)
    if (TARGET OpenMP::OpenMP_CXX)
    target_link_libraries(openmp_dlct PRIVATE OpenMP::OpenMP_CXX)
    else()
        if (OpenMP_CXX_FLAGS)
            target_compile_options(openmp_dlct PRIVATE ${OpenMP_CXX_FLAGS})
        endif()
        if (OpenMP_CXX_INCLUDE_DIRS)
            target_include_directories(openmp_dlct PRIVATE ${OpenMP_CXX_INCLUDE_DIRS})
        endif()
        if (OpenMP_CXX_LIBRARY)
            target_link_libraries(openmp_dlct PRIVATE ${OpenMP_CXX_LIBRARY})
        endif()
    endif()
else()
    message(STATUS "OpenMP not found. The openmp_dlct target will run in serial mode.")
endif()

if (APPLE)
    target_link_libraries(openmp_dlct PRIVATE "-framework Security")
endif()

install(TARGETS openmp_dlct DESTINATION bin)

if (ENABLE_CUDA)
    message(STATUS "CUDA architectures: ${CMAKE_CUDA_ARCHITECTURES}")
else()
    message(STATUS "CUDA architectures: disabled")
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
