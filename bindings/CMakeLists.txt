#########################
# I know there's a UseSWIG.cmake file.
# This did not seem to be very flexible at the current stage.

# Need C bindings for Python bindings in any case (Python wraps around C)
if (WITH_Python_BINDINGS)
    set(WITH_C_BINDINGS YES)
endif()
# Run checks
if (WITH_C_BINDINGS)
    find_package(Python COMPONENTS Interpreter Development)
    if (NOT Python_Interpreter_FOUND)
        message(FATAL_ERROR "No Python found. A working Python interpreter is required to generate any Iron bindings")
    endif()
endif()
# Directory used within tests
set(BINDINGS_BIN_DIR "${CMAKE_CURRENT_BINARY_DIR}")

set(_CLEANUP)
# C-Bindings extra target
if (WITH_C_BINDINGS)
    message(STATUS "Creating C bindings")
    # Generate the c files
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/opencmiss") # need this sub-path as the include directive is 
    set(OPENCMISS_H "${CMAKE_CURRENT_BINARY_DIR}/opencmiss/iron.h")
    set(OPENCMISS_C_F90 "${CMAKE_CURRENT_BINARY_DIR}/iron_c.f90")
    add_custom_command(OUTPUT ${OPENCMISS_C_F90} ${OPENCMISS_H}
        DEPENDS ${Iron_SOURCE_DIR}/src/opencmiss_iron.f90 # Need to re-build if that file changes!
        COMMAND "${Python_EXECUTABLE}" generate_bindings "${Iron_SOURCE_DIR}" C "${OPENCMISS_H}" "${OPENCMISS_C_F90}"
        WORKING_DIRECTORY ${Iron_SOURCE_DIR}/bindings
    )
    list(APPEND _CLEANUP ${OPENCMISS_H} ${OPENCMISS_C_F90})
    
    add_library(iron_c ${OPENCMISS_C_F90})
    target_link_libraries(iron_c PUBLIC iron)
    target_include_directories(iron_c PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
        $<BUILD_INTERFACE:${Iron_SOURCE_DIR}/src>
    )
    install(TARGETS iron_c
        EXPORT iron-config
        RUNTIME DESTINATION bin
        COMPONENT CBindings
        ARCHIVE DESTINATION lib
        COMPONENT CBindings
        LIBRARY DESTINATION lib
        COMPONENT Development
        INCLUDES DESTINATION include
    )
    install(FILES ${OPENCMISS_H}
        DESTINATION include/opencmiss
        COMPONENT CBindings
    )
    
    if (BUILD_TESTS)
        add_subdirectory(c/tests)
    endif()
    
    include(GenerateExportHeader)
    generate_export_header(iron_c)
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/iron_c_export.h"
        DESTINATION include/opencmiss
        COMPONENT CBindings
    )
endif()

# Python-Bindings extra target
if (WITH_Python_BINDINGS)
    message(STATUS "Creating Python bindings")
    find_package(Python COMPONENTS Interpreter Development)
    if(NOT Python_Development_FOUND)
        message(FATAL_ERROR "No Python development libraries found. Unable to generate Python bindings for Iron.")
    endif()
    # I know there's also a UseSWIG.cmake file (http://www.cmake.org/cmake/help/v3.3/module/UseSWIG.html)
    # This did not seem to be very flexible at the time of writing this.
    find_package(SWIG QUIET)
    if (NOT SWIG_FOUND)
        message(FATAL_ERROR "SWIG could not be found. Unable to generate Python bindings for Iron.")
    endif()
    execute_process(COMMAND "${Python_EXECUTABLE}" "${Iron_SOURCE_DIR}/utils/numpy_include.py"
        OUTPUT_VARIABLE NUMPY_INCLUDES
        RESULT_VARIABLE RES_NUMPY_INC
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if (NOT RES_NUMPY_INC EQUAL 0)
        message(FATAL_ERROR "Could not detect NumPy includes. Unable to generate Python bindings for Iron.")
    endif()
        
    set(IRON_PYTHON_MODULE iron_python)

    # Copy interface files to binary dir
    set(SWIG_INTERFACE_FILES ${CMAKE_CURRENT_BINARY_DIR}/numpy.i ${CMAKE_CURRENT_BINARY_DIR}/numpy_extra.i ${CMAKE_CURRENT_BINARY_DIR}/iron.i)
    add_custom_command(OUTPUT ${SWIG_INTERFACE_FILES}
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/python/numpy.i ${CMAKE_CURRENT_BINARY_DIR}
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/python/numpy_extra.i ${CMAKE_CURRENT_BINARY_DIR}
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/python/iron.i ${CMAKE_CURRENT_BINARY_DIR}
    )
    list(APPEND _CLEANUP ${SWIG_INTERFACE_FILES})
 
    # iron_generated.i is included by iron.i later
    set(SWIG_IFACE "${CMAKE_CURRENT_BINARY_DIR}/iron_generated.i")
    add_custom_command(OUTPUT ${SWIG_IFACE}
        DEPENDS ${Iron_SOURCE_DIR}/src/opencmiss_iron.f90 # Need to re-build if that file changes!
        COMMAND "${Python_EXECUTABLE}" generate_bindings "${Iron_SOURCE_DIR}" SWIG "${SWIG_IFACE}"
        COMMENT "Generating swig interface file for Iron"
        WORKING_DIRECTORY ${Iron_SOURCE_DIR}/bindings
    )

    if(Python_VERSION_MAJOR STREQUAL "3")
        set(PYTHONLIBS_SWIG_FLAGS -py3)
    endif()

    set(_PYTHON_DESTINATION_PREFIX lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR})
    set(IRON_PYTHON_PACKAGE_DESTINATION_PREFIX ${_PYTHON_DESTINATION_PREFIX}/$<CONFIG>/opencmiss.iron/opencmiss/iron)
    file(RELATIVE_PATH _REL_INSTALL_PATH ${CMAKE_INSTALL_PREFIX}/${IRON_PYTHON_PACKAGE_DESTINATION_PREFIX} ${CMAKE_INSTALL_PREFIX}/lib)

    # Generate C wrapper
    set(PYTHON_WRAPPER "${CMAKE_CURRENT_BINARY_DIR}/iron_python_wrapper.c")
    add_custom_command(OUTPUT ${PYTHON_WRAPPER}
        DEPENDS ${SWIG_IFACE}
        DEPENDS ${SWIG_INTERFACE_FILES}
        COMMAND ${SWIG_EXECUTABLE} -python ${PYTHONLIBS_SWIG_FLAGS} -o "${PYTHON_WRAPPER}"
           -module ${IRON_PYTHON_MODULE} -outdir ${CMAKE_CURRENT_BINARY_DIR} iron.i
        # we dont need the iron_python.py file generated by swig.
        COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_CURRENT_BINARY_DIR}/${IRON_PYTHON_MODULE}.py
        COMMENT "Generating python wrapper for Iron with swig"
    )
    list(APPEND _CLEANUP ${SWIG_IFACE} ${PYTHON_WRAPPER})

    #Generate wrapper object
    # MODULE: See http://www.cmake.org/cmake/help/v3.3/command/add_library.html
    add_library(${IRON_PYTHON_MODULE} MODULE ${PYTHON_WRAPPER} ${OPENCMISS_H} ${OPENCMISS_C_F90})
    target_link_libraries(${IRON_PYTHON_MODULE} PUBLIC iron_c Python::Python)
    target_include_directories(${IRON_PYTHON_MODULE} PRIVATE
        ${CMAKE_CURRENT_BINARY_DIR}
        ${PYTHON_INCLUDE_DIRS}
        ${NUMPY_INCLUDES}
    )
    set_target_properties(${IRON_PYTHON_MODULE} PROPERTIES PREFIX ""
        DEBUG_POSTFIX ""
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/$<CONFIG>/opencmiss/iron
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/$<CONFIG>/opencmiss/iron
    )
    if (APPLE)
        set_target_properties(${IRON_PYTHON_MODULE} PROPERTIES
            INSTALL_RPATH "@loader_path/${_REL_INSTALL_PATH};@loader_path/")
    elseif (UNIX)
        set_target_properties(${IRON_PYTHON_MODULE} PROPERTIES
            INSTALL_RPATH "\$ORIGIN/${_REL_INSTALL_PATH};\$ORIGIN/")
    endif ()

    # "Bugfix": If we use the same source file name, there will be two processes messing with
    # the same fortran module file which can lead to errors using -j parallel builds.
    # By this we simply force to have iron_c build before the python module.
    if (WITH_C_BINDINGS)
        add_dependencies(${IRON_PYTHON_MODULE} iron_c)
    endif()
    
    # We dont want any libXX or .dll - see UseSWIG.cmake
    set_target_properties(${IRON_PYTHON_MODULE} PROPERTIES PREFIX "_")
    if(WIN32 AND NOT CYGWIN)
        set_target_properties(${IRON_PYTHON_MODULE} PROPERTIES SUFFIX ".pyd")
    endif()
        
    # Generate Python module
    # Step 1: Configuration phase - replaces variables
    set(SETUP_PY ${CMAKE_CURRENT_BINARY_DIR}/setup.py)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/python/setup.in.py ${SETUP_PY} @ONLY)

    if (WIN32)
        file(TO_NATIVE_PATH "${CMAKE_INSTALL_PREFIX}/bin" NATIVE_LIBRARY_PATH)
    else()
        file(TO_NATIVE_PATH "${CMAKE_INSTALL_PREFIX}/lib" NATIVE_LIBRARY_PATH)
    endif()

    set(INIT_PY ${CMAKE_CURRENT_BINARY_DIR}/__init__.py)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/python/opencmiss/iron/__init__.in.py ${INIT_PY} @ONLY)
    
    set(UTILS_PY ${CMAKE_CURRENT_BINARY_DIR}/_utils.py)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/python/opencmiss/iron/_utils.in.py ${UTILS_PY} @ONLY)
            
    # Step 2: Generation phase - replaces remaining generator expressions
    set(IRON_PACKAGE_INIT_PY "${CMAKE_CURRENT_BINARY_DIR}/$<CONFIG>/opencmiss/iron/__init__.py")
    file(GENERATE OUTPUT "${IRON_PACKAGE_INIT_PY}"
        INPUT ${INIT_PY})
    set(IRON_SETUP_PY "${CMAKE_CURRENT_BINARY_DIR}/$<CONFIG>/setup.py")
    file(GENERATE OUTPUT "${IRON_SETUP_PY}"
        INPUT ${SETUP_PY})
    set(IRON_UTILS_PY "${CMAKE_CURRENT_BINARY_DIR}/$<CONFIG>/opencmiss/iron/_utils.py")
    file(GENERATE OUTPUT "${IRON_UTILS_PY}"
        INPUT ${UTILS_PY})
    set(OPENCMISS_PACKAGE_INIT_PY "${CMAKE_CURRENT_BINARY_DIR}/$<CONFIG>/opencmiss/__init__.py")
    file(GENERATE OUTPUT "${OPENCMISS_PACKAGE_INIT_PY}"
        INPUT ${CMAKE_CURRENT_SOURCE_DIR}/python/opencmiss/__init__.py)
    
    set(IRON_PY ${CMAKE_CURRENT_BINARY_DIR}/iron.py) # see python.py script in generate_bindings for iron.py name
    add_custom_command(OUTPUT ${IRON_PY}
        DEPENDS ${Iron_SOURCE_DIR}/src/opencmiss_iron.f90 # Need to re-build if that file changes!
        COMMAND "${Python_EXECUTABLE}" generate_bindings 
            "${Iron_SOURCE_DIR}" Python ${IRON_PYTHON_MODULE} "${CMAKE_CURRENT_BINARY_DIR}"
        COMMENT "Generating Python binding script"
        WORKING_DIRECTORY "${Iron_SOURCE_DIR}/bindings"
    )
    list(APPEND _CLEANUP ${IRON_PY})
    set(IRON_IRON_PY "${CMAKE_CURRENT_BINARY_DIR}/$<CONFIG>/opencmiss/iron/iron.py")
    add_custom_target(collect_python_binding_files ALL
        DEPENDS "${IRON_PY}"
        COMMAND ${CMAKE_COMMAND} -E copy ${IRON_PY} ${IRON_IRON_PY}
    )
    
    install(FILES ${IRON_SETUP_PY}
        DESTINATION ${_PYTHON_DESTINATION_PREFIX}/$<CONFIG>/opencmiss.iron
        COMPONENT PythonBindings)
    install(FILES ${OPENCMISS_PACKAGE_INIT_PY}
        DESTINATION ${_PYTHON_DESTINATION_PREFIX}/$<CONFIG>/opencmiss.iron/opencmiss
        COMPONENT PythonBindings)
    install(FILES ${IRON_PACKAGE_INIT_PY} ${IRON_IRON_PY} ${IRON_UTILS_PY}
        DESTINATION ${IRON_PYTHON_PACKAGE_DESTINATION_PREFIX}
        COMPONENT PythonBindings)
    install(TARGETS ${IRON_PYTHON_MODULE}
        DESTINATION ${IRON_PYTHON_PACKAGE_DESTINATION_PREFIX}
        COMPONENT PythonBindings)
    
    if (IRON_USE_VIRTUALENV)
        # Variables required for OCPythonBindingsVirtualEnv.
        set(VIRTUALENV_INSTALL_PREFIX ${IRON_VIRTUALENV_INSTALL_PREFIX})
        set(PYTHON_PACKAGE_CURRENT_NAME iron)
        # Virtual environment creation target.
        include(OCPythonBindingsVirtualEnv)    
    endif ()

    if (BUILD_TESTS)
        add_subdirectory(python/tests)
    endif()
endif()
if (_CLEANUP)
    set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${_CLEANUP}")
endif()

