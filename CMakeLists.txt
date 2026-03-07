cmake_minimum_required(VERSION 3.16)

project(calc_coupling
    VERSION 2.0
    DESCRIPTION "2-states FMO dimer coupling calculator"
    LANGUAGES Fortran
)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
file(MAKE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY})

# --- 添加选项：是否启用 MKL ---
option(USE_MKL "Enable Intel MKL for linear algebra" ON)

# --- 优先级查找 MKL（仅在 USE_MKL=ON 时尝试） ---
if(USE_MKL)
    find_package(MKL CONFIG QUIET)
    if(NOT MKL_FOUND)
        find_package(MKL QUIET)
    endif()
endif()

if(MKL_FOUND AND USE_MKL)
    message(STATUS "Using Intel MKL")
    set(LINALG_LIBRARIES MKL::MKL)
else()
    message(STATUS "Intel MKL disabled or not found. Using generic LAPACK/BLAS.")
    find_package(LAPACK REQUIRED)
    set(LINALG_LIBRARIES ${LAPACK_LIBRARIES})
endif()

if(MKL_FOUND)
    message(STATUS "Using Intel MKL")
    set(LINALG_LIBRARIES MKL::MKL)
else()
    message(STATUS "Intel MKL not found. Trying generic LAPACK/BLAS.")
    find_package(LAPACK REQUIRED)
    if(LAPACK_FOUND)
        set(LINALG_LIBRARIES ${LAPACK_LIBRARIES})
    else()
        message(FATAL_ERROR "Neither Intel MKL nor generic LAPACK/BLAS found. Please load one of them.")
    endif()
endif()

message(STATUS "Linear algebra library: ${LINALG_LIBRARIES}")

set(SOURCES
    src/coupling_module.f90
    src/calculator_basic.f90
    src/basic_reader.f90
    src/matrix_reader.f90
    src/coupling.f90
    src/main.f90
)

add_executable(${PROJECT_NAME} ${SOURCES})
target_link_libraries(${PROJECT_NAME} PRIVATE ${LINALG_LIBRARIES})
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})

install(TARGETS ${PROJECT_NAME} DESTINATION bin)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

message(STATUS "")
message(STATUS "=== Configuration Summary ===")
message(STATUS "Project: ${PROJECT_NAME}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Fortran compiler: ${CMAKE_Fortran_COMPILER}")
message(STATUS "Linear algebra library: ${LINALG_LIBRARIES}")
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "=============================")
