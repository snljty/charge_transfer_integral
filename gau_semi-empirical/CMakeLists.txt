cmake_minimum_required(VERSION 3.18) # since 3.18 supports module FindLAPACK

enable_language(Fortran)

project(calc_coupling_gau_semi_empirical_use_out)

# list(PREPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/Modules")

option(BUILD_SHARED_LINKED_EXE "Build shared linked executable : ON OFF (BOOL)" OFF)
set(BUILD_SHARED_LINKED_EXE_VALID_OPTIONS ON OFF)
list(FIND BUILD_SHARED_LINKED_EXE_VALID_OPTIONS ${BUILD_SHARED_LINKED_EXE} BUILD_SHARED_LINKED_EXE_VALID_OPTION_INDEX)
if (BUILD_SHARED_LINKED_EXE_VALID_OPTION_INDEX EQUAL -1)
    message(FATAL_ERROR "BUILD_SHARED_LINKED_EXE can only be ON or OFF")
endif ()

if (${BUILD_SHARED_LINKED_EXE} STREQUAL OFF)
    string(APPEND CMAKE_EXE_LINKER_FLAGS " -static")
endif ()

add_executable(${PROJECT_NAME} src/calc_coupling_gau_semi_empirical_use_out.f90)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif ()

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffpe-summary=none")

set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})

if ("${CMAKE_INSTALL_PREFIX}" STREQUAL "${CMAKE_BINARY_DIR}")
    message(FATAL_ERROR "The source directory cannot be the same as install directory.")
endif ()

install(TARGETS ${PROJECT_NAME} DESTINATION bin)

