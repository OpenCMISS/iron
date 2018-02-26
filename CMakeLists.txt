cmake_minimum_required (VERSION 3.3 FATAL_ERROR)
project(Iron VERSION 0.7.0 LANGUAGES C CXX Fortran)# Need all languages as iron's deps cover all (link phase) 
include(CMakeDependentOption)

set(_CORRECT_CMAKE_MODULE_PATH FALSE)
# First check if the CMAKE_MODULE_PATH is already set properly.
foreach(_PATH ${CMAKE_MODULE_PATH})
    if (EXISTS "${_PATH}/OCMiscFunctions.cmake")
        set(_CORRECT_CMAKE_MODULE_PATH TRUE)
        break()
    endif ()
endforeach()

if (NOT _CORRECT_CMAKE_MODULE_PATH)
    # Check bare minimum configuration
    set(IRON_CMAKE_MODULE_PATH "${IRON_CMAKE_MODULE_PATH}" CACHE STRING "Location of OpenCMISS CMake modules.")
    if (EXISTS "${IRON_CMAKE_MODULE_PATH}/OpenCMISS/OCMiscFunctions.cmake")
        list(APPEND CMAKE_MODULE_PATH "${IRON_CMAKE_MODULE_PATH}/FindModuleWrappers")
        list(APPEND CMAKE_MODULE_PATH "${IRON_CMAKE_MODULE_PATH}")
        list(APPEND CMAKE_MODULE_PATH "${IRON_CMAKE_MODULE_PATH}/OpenCMISS")
        set(_CORRECT_CMAKE_MODULE_PATH TRUE)
    endif ()
endif ()

if (NOT _CORRECT_CMAKE_MODULE_PATH)
    message(STATUS "CMAKE_MODULE_PATH: '${CMAKE_MODULE_PATH}'")
    message(STATUS "IRON_CMAKE_MODULE_PATH: '${IRON_CMAKE_MODULE_PATH}'")
    message(FATAL_ERROR "'OCMiscFunctions.cmake' was not found on the module path.  Check the value of CMAKE_MODULE_PATH and/or IRON_CMAKE_MODULE_PATH.")
endif ()

string(TOUPPER "${CMAKE_BUILD_TYPE}" _uppercase_CMAKE_BUILD_TYPE)
if (_uppercase_CMAKE_BUILD_TYPE MATCHES "^(DEBUG|RELWITHDEBINFO)$")
  set(IS_DEBUG_BUILD TRUE)
else()
  set(IS_DEBUG_BUILD FALSE)
endif()

# Options for building iron
option(IRON_BUILD_SHARED_LIBRARY "Build a shared iron library." ON)
option(IRON_BUILD_STATIC_LIBRARY "Build a static iron library." OFF)
set(REAL_PRECISION "double" CACHE STRING "Real precision (single|double)")
option(WITH_CELLML "Build iron with CellML support" ON)
option(WITH_FIELDML "Build iron with FieldML support" ON)
option(WITH_HYPRE "Build iron with HYPRE integration" ON)
option(WITH_SUNDIALS "Build iron with SUNDIALS integration" ON)
option(WITH_MUMPS "Build iron with MUMPS integration" ON)
option(WITH_SCALAPACK "Build iron with SCALAPACK integration" ON)
option(WITH_PETSC "Build iron with PETSC integration" ON)
option(BUILD_TESTS "${PROJECT_NAME} - Build tests" ON)
option(WITH_PROFILING "${PROJECT_NAME} - Build with profiling flags" OFF)
CMAKE_DEPENDENT_OPTION(WITH_DIAGNOSTICS "Build iron with diagnostics support" ON
                       "IS_DEBUG_BUILD" OFF)
option(PRINT_CONFIG_SUMMARY "Show a summary of the configuration." TRUE)
set(PACKAGE_CONFIG_DIR "cmake" CACHE STRING "Directory for package config files (relative to CMAKE_INSTALL_PREFIX)")

# Handle instrumentation choices.
set(IRON_INSTRUMENTATION "none" CACHE STRING "Specify the instrumentation to use.")
if (DEFINED FE_INSTRUMENTATION)
    set(IRON_INSTRUMENTATION "${FE_INSTRUMENTATION}" CACHE STRING "Specify the instrumentation to use.")
endif ()
set(OPENCMISS_INSTRUMENTATION ${IRON_INSTRUMENTATION})

# Bindings options
find_program(VIRTUALENV_EXEC virtualenv)
mark_as_advanced(VIRTUALENV_EXEC)

option(WITH_C_BINDINGS "Build iron C bindings" ON)
option(WITH_Python_BINDINGS "Build iron Python bindings" ON)
set(IRON_VIRTUALENV_INSTALL_PREFIX "${IRON_VIRTUALENV_INSTALL_PREFIX}" CACHE PATH "Install prefix for the virtualenv Iron library.")
if (DEFINED FE_VIRTUALENV_INSTALL_PREFIX)
    set(IRON_VIRTUALENV_INSTALL_PREFIX "${FE_VIRTUALENV_INSTALL_PREFIX}" CACHE PATH "Install prefix for the virtualenv Iron library." FORCE)
endif ()

set(IRON_USE_VIRTUALENV TRUE CACHE BOOL "Install Iron into virtual environment.")
if (DEFINED FE_USE_VIRTUALENV)
    set(IRON_USE_VIRTUALENV ${FE_USE_VIRTUALENV} CACHE BOOL "Install Iron into virtual environment." FORCE)
endif ()
if (NOT EXISTS "${VIRTUALENV_EXEC}")
    set(IRON_USE_VIRTUALENV FALSE CACHE BOOL "Install Iron into virtual environment." FORCE)
    message(STATUS "The program virtualenv was not found disabling it's use.")
endif ()

string(TOUPPER "${REAL_PRECISION}" _uppercase_REAL_PRECISION)
if(_uppercase_REAL_PRECISION STREQUAL "SINGLE")
  set(WORKING_REAL_PRECISION_SINGLE)
elseif(_uppercase_REAL_PRECISION STREQUAL "DOUBLE")
  set(WORKING_REAL_PRECISION_DOUBLE)
else()
  message(STATUS "The specified real precision of '${REAL_PRECISION}' is invalid. Defaulting to double.")
  set(WORKING_REAL_PRECISION_DOUBLE)
endif()

if(NOT IRON_BUILD_SHARED_LIBRARY AND NOT IRON_BUILD_STATIC_LIBRARY )
  message(FATAL_ERROR "You have to build at least one of the types of libraries (shared or static)." )
endif()

include(OCMiscFunctions)
include(OCMultiConfigEnvironment)

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
set(IRON_UPSTREAM_DEPS SCALAPACK PARMETIS MUMPS HYPRE SUNDIALS PETSC)
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
        message(STATUS "Looking for ${DEP}.. found at ${DEP_LOC}")
        string(TOLOWER ${DEP} targetname)
        LIST(APPEND LINK_LIBS ${targetname})
    endif()
endforeach()

# Platform-checks
include(CheckSymbolExists)
CHECK_SYMBOL_EXISTS(SA_NODEFER "signal.h" HAVE_SA_NODEFER)
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
    source_group(headers FILES ${IRON_HEADERS})
    source_group(src FILES ${IRON_C_SRC})
    source_group(src FILES ${IRON_Fortran_SRC})
else()
    set(IRON_TARGETS iron)
    add_library(iron ${IRON_SRC})
    target_compile_options(iron PUBLIC $<BUILD_INTERFACE:-cpp>)
    source_group(headers FILES ${IRON_HEADERS})
    source_group(src FILES ${IRON_Fortran_SRC} ${IRON_C_SRC})
endif()

# Set the same stuff for either the separate targets (windows) or one target (else)
foreach(irontarget ${IRON_TARGETS})
    target_include_directories(${irontarget} PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/src"
        ${CMAKE_CURRENT_BINARY_DIR}
    )
    if(WORKING_REAL_PRECISION_SINGLE)
        target_compile_definitions(${irontarget} PUBLIC SINGLE_REAL_PRECISION)
    else()
        target_compile_definitions(${irontarget} PUBLIC DOUBLE_REAL_PRECISION)
    endif()
    if (IS_DEBUG_BUILD)
        target_compile_definitions(${irontarget} PUBLIC DEBUG)
    endif()
    if (WITH_CELLML)
        target_compile_definitions(${irontarget} PUBLIC WITH_CELLML)
    endif()
    if (WITH_FIELDML)
        target_compile_definitions(${irontarget} PUBLIC WITH_FIELDML)
    endif()
    if (WITH_HYPRE)
        target_compile_definitions(${irontarget} PUBLIC WITH_HYPRE)
    endif()
    if (WITH_SUNDIALS)
        target_compile_definitions(${irontarget} PUBLIC WITH_SUNDIALS)
    endif()
    if (WITH_MUMPS)
        target_compile_definitions(${irontarget} PUBLIC WITH_MUMPS)
    endif()
    if (WITH_SCALAPACK)
        target_compile_definitions(${irontarget} PUBLIC WITH_SCALAPACK)
    endif()
    if (WITH_PETSC)
        target_compile_definitions(${irontarget} PUBLIC WITH_PETSC)
    endif()
    if (WITH_DIAGNOSTICS)
        target_compile_definitions(${irontarget} PUBLIC WITH_DIAGNOSTICS)
    endif()
    
    target_link_libraries(${irontarget} PUBLIC MPI::MPI_Fortran)
    target_link_libraries(${irontarget} PRIVATE ${LINK_LIBS})
    if(WITH_FIELDML)
        target_link_libraries(${irontarget} PUBLIC $<BUILD_INTERFACE:fieldml-api>)
    endif()
    # Dont use "USE MPI" instructions if the precompiled mpi module files are incompatible.
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
            if (OPENCMISS_HAVE_MULTICONFIG_ENV)
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
install(TARGETS ${IRON_TARGETS}
    EXPORT iron-config
    RUNTIME DESTINATION bin
    COMPONENT Runtime
    ARCHIVE DESTINATION lib
    COMPONENT Runtime
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
file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/iron-config-dependencies.cmake "if (${IRON_BUILD_STATIC_LIBRARY})\r\n")
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
file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/iron-config-dependencies.cmake "endif ()\r\n")
file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/iron-config-dependencies.cmake "set(_IMPORT_PREFIX \${${PROJECT_NAME}_IMPORT_PREFIX})\r\n")

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

# Bundle fixup for iron shared library
# We need non-system libraries shipped out with the Runtime libraries.
# This will probably also work on other platforms - need to check
if (WIN32)
    # The wild turkey is out! We need this generator-intermediate step to get a hold of the effective output location
    # of the iron shared library. This in turn in placed inside a directory that is configuration-dependent,
    # and in the cmake_install script the CMAKE_INSTALL_CONFIG_NAME helps us getting the job done.
    file(GENERATE OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/$<CONFIG>/wildturkey.cmake"
        CONTENT "set(IRON_LIB_NAME \"$<TARGET_FILE:iron>\")")
    # Hack: The BundleUtilities dont seem to work unless you give them an executable (just the dll wont work),
    # and as we dont always have the tests around we just use the cmake executable, which will always be around.
    install(FILES ${CMAKE_COMMAND}
        DESTINATION bin)
    get_filename_component(CMAKE_NAME "${CMAKE_COMMAND}" NAME)
    install(CODE "
       include(BundleUtilities)
       include(${CMAKE_CURRENT_BINARY_DIR}/\${CMAKE_INSTALL_CONFIG_NAME}/wildturkey.cmake)
       file(COPY \"${CMAKE_COMMAND}\" DESTINATION \"\${CMAKE_INSTALL_PREFIX}/bin\")
       fixup_bundle(\"\${CMAKE_INSTALL_PREFIX}/bin/${CMAKE_NAME}\" \"\${IRON_LIB_NAME}\" \"\")
       file(REMOVE \"\${CMAKE_INSTALL_PREFIX}/bin/${CMAKE_NAME}\")
       "
       COMPONENT Redist)
endif()

if (PRINT_CONFIG_SUMMARY)
	message(STATUS   "")
	message(STATUS   "=====================================================")
	message(STATUS   "Iron has been configured with the following settings:")
	message(STATUS   "")
	message(STATUS   "          Build type: '${CMAKE_BUILD_TYPE}'")
	if (IS_DEBUG_BUILD)
	  message(STATUS "      Is debug build: TRUE")
	else()
	  message(STATUS "      Is debug build: FALSE")
	endif()
        if(WORKING_REAL_PRECISION_SINGLE)
	  message(STATUS "      Real precision: SINGLE")
	else()
	  message(STATUS "      Real precision: DOUBLE")
	endif()
	if (IRON_BUILD_SHARED_LIBRARY)
	  message(STATUS "        Library type: SHARED")
	elseif(IRON_BUILD_STATIC_LIBRARY)
	  message(STATUS "        Library type: STATIC")
	endif()
	message(STATUS   "         Install dir: '${CMAKE_INSTALL_PREFIX}'")
	if (BUILD_TESTS)
	  message(STATUS "         Build tests: ON")
	else()
	  message(STATUS "         Build tests: OFF")
	endif()
	if (WITH_PROFILING)
	  message(STATUS "           Profiling: ON")
	else()
	  message(STATUS "           Profiling: OFF")
	endif()
	if (WITH_C_BINDINGS)
	  message(STATUS "          C Bindings: ON")
	else()
	  message(STATUS "          C Bindings: OFF")
	endif()
	if (WITH_Python_BINDINGS)
	  message(STATUS "     Python Bindings: ON")
	else()
	  message(STATUS "     Python Bindings: OFF")
	endif()
	if (WITH_DIAGNOSTICS)
	  message(STATUS "         Diagnostics: ON")
	else()
	  message(STATUS "         Diagnostics: OFF")
	endif()
	if (WITH_CELLML)
	  message(STATUS "              CellML: ON")
	else()
	  message(STATUS "              CellML: OFF")
	endif()
	if (WITH_FIELDML)
	  message(STATUS "             FieldML: ON")
	else()
	  message(STATUS "             FieldML: OFF")
	endif()
	if (WITH_HYPRE)
	  message(STATUS "               Hypre: ON")
	else()
	  message(STATUS "               Hypre: OFF")
	endif()
	if (WITH_SUNDIALS)
	  message(STATUS "            SUNDIALS: ON")
	else()
	  message(STATUS "            SUNDIALS: OFF")
	endif()
	if (WITH_MUMPS)
	  message(STATUS "               MUMPS: ON")
	else()
	  message(STATUS "               MUMPS: OFF")
	endif()
	if (WITH_SCALAPACK)
	  message(STATUS "           ScaLAPACK: ON")
	else()
	  message(STATUS "           ScaLAPACK: OFF")
	endif()
	if (WITH_PETSC)
	  message(STATUS "               PETSc: ON")
	else()
	  message(STATUS "               PETSc: OFF")
	endif()
	message(STATUS   "=====================================================")
	message(STATUS   "")
endif()
