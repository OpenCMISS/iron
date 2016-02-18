cmake_minimum_required (VERSION 3.3 FATAL_ERROR)
project(Iron VERSION 0.5.0 LANGUAGES C CXX Fortran)# Need all languages as iron's deps cover all (link phase) 

option(WITH_CELLML "Build iron with CellML support" ON)
option(WITH_FIELDML "Build iron with FieldML support" ON)
option(WITH_HYPRE "Build iron with HYPRE integration" ON)
option(WITH_SUNDIALS "Build iron with SUNDIALS integration" ON)
option(WITH_MUMPS "Build iron with MUMPS integration" ON)
option(WITH_SCALAPACK "Build iron with SCALAPACK integration" ON)
option(WITH_PETSC "Build iron with PETSC integration" ON)
option(BUILD_TESTS "${PROJECT_NAME} - Build tests" ON)
option(WITH_PROFILING "${PROJECT_NAME} - Build with profiling flags" OFF)
option(WITH_C_BINDINGS "Build iron C bindings" ON)
option(WITH_Python_BINDINGS "Build iron Python bindings" ON)
set(PACKAGE_CONFIG_DIR "cmake" CACHE STRING "Directory for package config files (relative to CMAKE_INSTALL_PREFIX)")

# Config
if (MSVC)
    set(HAVE_MULTICONFIG_ENV YES)
endif()

# See http://cmake.3232098.n2.nabble.com/How-to-install-Fortran-module-files-td7584399.html
SET(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fortran_modules)
SET(CMAKE_NO_SYSTEM_FROM_IMPORTED TRUE)
# Have cmake find our own scripts
LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# External components
# MPI
find_package(MPI REQUIRED)
 
# Collect libraries to configure for/link against
SET(LINK_LIBS )
if (WITH_CELLML)
    find_package(CELLML ${CELLML_VERSION} CONFIG REQUIRED)
    message(STATUS "Building with CELLML (${LIBCELLML_DIR})")
    LIST(APPEND LINK_LIBS cellml cellml_model_definition)
endif()
if (WITH_FIELDML)
    find_package(FIELDML-API ${FIELDML-API_VERSION} CONFIG REQUIRED)
    message(STATUS "Building with FIELDML (${FIELDML-API_DIR})")
endif()

# This include file is in the main "manage" directory.
# We can include this here as iron is an integrated part of the overall build system.
SET(IRON_UPSTREAM_DEPS PETSC SUNDIALS HYPRE MUMPS SCALAPACK PARMETIS)
foreach(DEP ${IRON_UPSTREAM_DEPS})
    if (WITH_${DEP})
        message(STATUS "Looking for ${DEP}..")
        find_package(${DEP} ${${DEP}_VERSION} REQUIRED)
        if (${DEP}_DIR)
            # This is set if CONFIG mode was used
            set(DEP_LOC ${${DEP}_DIR})
        else()
            # This is set if MODULE mode was used
            if (${DEP}_LIBRARIES)
                set(_ALTERNATIVE ${DEP}_LIBRARIES)
            else()
                set(_ALTERNATIVE ${DEP}_INCLUDE_PATH)
            endif()
            list(GET ${_ALTERNATIVE} 0 DEP_LOC)
            get_filename_component(DEP_LOC ${DEP_LOC} PATH)
        endif()
        #message(STATUS "Looking for ${DEP}.. found at ${DEP_LOC}")
        string(TOLOWER ${DEP} targetname)
        LIST(APPEND LINK_LIBS ${targetname})
    endif()
endforeach()

# Platform-checks
include(CheckSymbolExists)
CHECK_SYMBOL_EXISTS(SA_NODEFER "" HAVE_SA_NODEFER)
CHECK_SYMBOL_EXISTS(sigaction "signal.h" HAVE_SIGACTION_STRUCT)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmiss.in.h ${CMAKE_CURRENT_BINARY_DIR}/cmiss.h)

# Populates the IRON_XX_SRC variables
include(Sources)

if(MSVC)
    # Need separate C/Fortran code due to Project natures in VS
    set(IRON_TARGETS iron iron_static_C)
    
    add_library(iron_static_C STATIC ${IRON_C_SRC})
    
    add_library(iron ${IRON_Fortran_SRC})
    target_link_libraries(iron PRIVATE iron_static_C)
    
    # IDE stuff
    source_group(src FILES ${IRON_C_SRC})
    source_group(src FILES ${IRON_Fortran_SRC})
else()
    set(IRON_TARGETS iron)
    add_library(iron ${IRON_SRC})
    target_compile_options(iron PUBLIC $<BUILD_INTERFACE:-cpp>)
endif()

# Set the same stuff for either the separate targets (windows) or one target (else)
foreach(irontarget ${IRON_TARGETS})
    target_include_directories(${irontarget} PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/src"
        ${CMAKE_CURRENT_BINARY_DIR}
    )
    if (WITH_CELLML)
        target_compile_definitions(${irontarget} PUBLIC USECELLML)
    endif()
    if (WITH_FIELDML)
        target_compile_definitions(${irontarget} PUBLIC USEFIELDML)
    endif()
    # In case the compiler is not already the mpi wrapper,
    # we need all the include paths etc.
    # This is conveniently done by linking against the "mpi" interface target
    # returned by our FindMPI module.
    if (NOT CMAKE_Fortran_COMPILER STREQUAL MPI_Fortran_COMPILER)
        target_link_libraries(${irontarget} PUBLIC $<BUILD_INTERFACE:mpi>)
    endif()
    target_link_libraries(${irontarget} PRIVATE ${LINK_LIBS})
    if(WITH_FIELDML)
        target_link_libraries(${irontarget} PUBLIC $<BUILD_INTERFACE:fieldml-api>)
    endif()
    # Dont use "USE MPI" instructions if the percompiled mpi module files are incompatible.
    # This happens when the fortran compiler used to build mpi differs from the currently used one.
    # This holds on windows (only pre-compiled mpi implementations available; MPICH2 and MSMPI)
    # and on e.g. clusters with inappropriately configured toolchain environments.
    if (NOT MPI_Fortran_MODULE_COMPATIBLE)
        target_compile_definitions(${irontarget} PUBLIC $<BUILD_INTERFACE:NOMPIMOD>)
    endif()
endforeach()

if (BUILD_TESTS)
    enable_testing()
    
    # Composes a native PATH-compatible variable to use for DLL/SO finding.
    # Each extra argument is assumed a path to add. Added in the order specified.
    function(get_library_path OUTPUT_VARIABLE)
        if (WIN32)
            set(PSEP "\\;")
            set(LD_VARNAME "PATH")
        elseif(APPLE)
            set(LD_VARNAME "DYLD_LIBRARY_PATH")
            set(PSEP ":")
        elseif(UNIX)
            set(LD_VARNAME "LD_LIBRARY_PATH")
            set(PSEP ":")
        else()
            message(WARNING "get_library_path not implemented for '${CMAKE_HOST_SYSTEM}'")
        endif()
        # Load system environment - on windows its separated by semicolon, so we need to protect those 
        string(REPLACE ";" "\\;" LD_PATH "$ENV{${LD_VARNAME}}")
        foreach(_PATH ${ARGN})
            if (HAVE_MULTICONFIG_ENV)
                file(TO_NATIVE_PATH "${_PATH}/$<CONFIG>" _PATH)
            else()
                file(TO_NATIVE_PATH "${_PATH}" _PATH)
            endif()
            set(LD_PATH "${_PATH}${PSEP}${LD_PATH}")
        endforeach()
        set(${OUTPUT_VARIABLE} "${LD_VARNAME}=${LD_PATH}" PARENT_SCOPE)
    endfunction()
    
    function(add_test_path TESTNAME)
        get_library_path(LD_PATH ${ARGN})
        set_tests_properties(${TESTNAME} PROPERTIES
            ENVIRONMENT "${LD_PATH}"
        )
    endfunction()
    
    add_subdirectory(tests)
endif()

# Build bindings as selected
add_subdirectory(bindings)

# Install steps
install(TARGETS iron
    EXPORT iron-config
    RUNTIME DESTINATION bin
    COMPONENT Runtime
    ARCHIVE DESTINATION lib
    COMPONENT Runtime
    INCLUDES DESTINATION include/opencmiss
    LIBRARY DESTINATION lib
    COMPONENT Development
    INCLUDES DESTINATION include/opencmiss
)

install(EXPORT iron-config
    DESTINATION ${PACKAGE_CONFIG_DIR}
    COMPONENT Development
)
include(CMakePackageConfigHelpers)
WRITE_BASIC_PACKAGE_VERSION_FILE(${CMAKE_CURRENT_BINARY_DIR}/iron-config-version.cmake COMPATIBILITY AnyNewerVersion)
install(FILES 
    ${CMAKE_CURRENT_BINARY_DIR}/iron-config-version.cmake
    DESTINATION ${PACKAGE_CONFIG_DIR}
    COMPONENT Development
)
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/iron-config-dependencies.cmake
    "include(CMakeFindDependencyMacro)\r\n"
    "set(${PROJECT_NAME}_IMPORT_PREFIX \${_IMPORT_PREFIX})\r\n"
)
foreach(IRON_DEP ${IRON_UPSTREAM_DEPS})
    if (WITH_${IRON_DEP})
        file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/iron-config-dependencies.cmake "find_dependency(${IRON_DEP})\r\n")
    endif()
endforeach()
if (WITH_CELLML)
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/iron-config-dependencies.cmake "find_dependency(CELLML)\r\n")
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/iron-config-dependencies.cmake "find_dependency(LIBCELLML)\r\n")
endif()
if (WITH_FIELDML)
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/iron-config-dependencies.cmake "find_dependency(FIELDML-API)\r\n")
endif()
file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/iron-config-dependencies.cmake "set(_IMPORT_PREFIX \${${PROJECT_NAME}_IMPORT_PREFIX})")
install(FILES 
    ${CMAKE_CURRENT_BINARY_DIR}/iron-config-dependencies.cmake
    DESTINATION ${PACKAGE_CONFIG_DIR}
    COMPONENT DevelopmentSDK)    

# The install directories are "opencmiss" instead of "iron", as those packages are to be seen under one name
# rather than separate entities.
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/
    DESTINATION include/opencmiss
    COMPONENT Development
    FILES_MATCHING PATTERN "opencmiss*.mod"
)
