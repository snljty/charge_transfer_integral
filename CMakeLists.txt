cmake_minimum_required(VERSION 3.18) # since 3.18 supports module FindLAPACK

enable_language(Fortran)

project(calc_coupling)

# list(PREPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/Modules")

option(BUILD_SHARED_LINKED_EXE "Build shared linked executable : ON OFF (BOOL)" OFF)
set(BUILD_SHARED_LINKED_EXE_VALID_OPTIONS ON OFF)
list(FIND BUILD_SHARED_LINKED_EXE_VALID_OPTIONS ${BUILD_SHARED_LINKED_EXE} BUILD_SHARED_LINKED_EXE_VALID_OPTION_INDEX)
if (BUILD_SHARED_LINKED_EXE_VALID_OPTION_INDEX EQUAL -1)
    message(FATAL_ERROR "BUILD_SHARED_LINKED_EXE can only be ON or OFF")
endif ()

find_package(LAPACK REQUIRED)

if (${BUILD_SHARED_LINKED_EXE} STREQUAL OFF)
    string(APPEND CMAKE_EXE_LINKER_FLAGS " -static")
    set(BLA_STATIC ON)
    if (${BLA_STATIC} STREQUAL ON)
        if (${CMAKE_SYSTEM_NAME} MATCHES Windows)
            list(GET LAPACK_LIBRARIES 0 LAPACK_ITSELF_LIBRARY)
            get_filename_component(LAPACK_LIB_DIR ${LAPACK_ITSELF_LIBRARY} DIRECTORY)
            unset(LAPACK_ITSELF_LIBRARY)
            foreach (LAPACK_LIB ${LAPACK_LIBRARIES})
                get_filename_component(LAPACK_LIB_NAME ${LAPACK_LIB} NAME)
                # get_filename_component(LAPACK_LIB_NAME_STATIC ${LAPACK_LIB} NAME_WE)
                # string(APPEND LAPACK_LIB_NAME_STATIC ".a")
                string(REPLACE ".dll.a" ".a" LAPACK_LIB_STATIC_NAME ${LAPACK_LIB_NAME})
                list(APPEND LAPACK_LIBRARIES_SWP ${LAPACK_LIB_DIR}/${LAPACK_LIB_STATIC_NAME})
            endforeach ()
            unset(LAPACK_LIB)
            unset(LAPACK_LIB_NAME)
            unset(LAPACK_LIB_STATIC_NAME)
            set(LAPACK_LIBRARIES ${LAPACK_LIBRARIES_SWP})
            unset(LAPACK_LIBRARIES_SWP)
            message("-- Will use static libraries LAPACK_LIBRARIES ${LAPACK_LIBRARIES} .")
        else ()
            list(GET LAPACK_LIBRARIES 0 LAPACK_ITSELF_LIBRARY)
            get_filename_component(LAPACK_LIB_DIR ${LAPACK_ITSELF_LIBRARY} DIRECTORY)
            unset(LAPACK_ITSELF_LIBRARY)
            foreach (LAPACK_LIB ${LAPACK_LIBRARIES})
                get_filename_component(LAPACK_LIB_NAME ${LAPACK_LIB} NAME)
                # get_filename_component(LAPACK_LIB_NAME_STATIC ${LAPACK_LIB} NAME_WE)
                # string(APPEND LAPACK_LIB_NAME_STATIC ".a")
                string(REPLACE ".so" ".a" LAPACK_LIB_STATIC_NAME ${LAPACK_LIB_NAME})
                list(APPEND LAPACK_LIBRARIES_SWP ${LAPACK_LIB_DIR}/${LAPACK_LIB_STATIC_NAME})
            endforeach ()
            unset(LAPACK_LIB)
            unset(LAPACK_LIB_NAME)
            unset(LAPACK_LIB_STATIC_NAME)
            set(LAPACK_LIBRARIES ${LAPACK_LIBRARIES_SWP})
            unset(LAPACK_LIBRARIES_SWP)
            message("-- Will use static libraries LAPACK_LIBRARIES ${LAPACK_LIBRARIES} .")
        endif ()
    endif ()
endif ()

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffpe-summary=none")

add_executable(${PROJECT_NAME} src/calc_coupling.f90)

target_link_libraries(${PROJECT_NAME} PRIVATE ${LAPACK_LIBRARIES})

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif ()

set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})

if ("${CMAKE_INSTALL_PREFIX}" STREQUAL "${CMAKE_BINARY_DIR}")
    message(FATAL_ERROR "The source directory cannot be the same as install directory.")
endif ()

install(TARGETS ${PROJECT_NAME} DESTINATION bin)

