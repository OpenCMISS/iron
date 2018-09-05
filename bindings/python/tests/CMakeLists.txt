# Set the correct python path to have it find the opencmiss package
file(TO_NATIVE_PATH "${BINDINGS_BIN_DIR}/$<CONFIG>" NATIVE_PYTHON_BINDINGS_DIR)

# Add the tests
set(PYTHONTESTS )
add_test(NAME Python_Bindings_Import
    COMMAND ${Python_EXECUTABLE} -c "from opencmiss.iron import iron"
)
list(APPEND PYTHONTESTS Python_Bindings_Import)

add_test(NAME Python_Bindings_Cantilever
    COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Cantilever.py
)
list(APPEND PYTHONTESTS Python_Bindings_Cantilever)


add_test(NAME Python_Bindings_CellMLModelIntegration_n98_Euler
    COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/CellMLModelIntegration.py ${Iron_SOURCE_DIR}/tests/CellML/n98.xml euler 0.01
)
list(APPEND PYTHONTESTS Python_Bindings_CellMLModelIntegration_n98_Euler)

add_test(NAME Python_Bindings_CellMLModelIntegration_n98_BDF
    COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/CellMLModelIntegration.py ${Iron_SOURCE_DIR}/tests/CellML/n98.xml bdf 0.01
)
list(APPEND PYTHONTESTS Python_Bindings_CellMLModelIntegration_n98_BDF)

add_test(NAME Python_Bindings_CellMLModelIntegration_n98_Euler_0.005
    COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/CellMLModelIntegration.py ${Iron_SOURCE_DIR}/tests/CellML/n98.xml euler 0.005
)
list(APPEND PYTHONTESTS Python_Bindings_CellMLModelIntegration_n98_Euler_0.005)

add_test(NAME Python_Bindings_Laplace
    COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Laplace.py
)
list(APPEND PYTHONTESTS Python_Bindings_Laplace)

# Does not seem to converge! Needs checking
#add_test(NAME Python_Bindings_Monodomain2DSquare
#    COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Monodomain2DSquare.py ${Iron_SOURCE_DIR}/tests/CellML/n98.xml
#)
#list(APPEND PYTHONTESTS Python_Bindings_Monodomain2DSquare)

get_library_path(PATH_DEF "${Iron_BINARY_DIR}" "${BINDINGS_BIN_DIR}")
foreach(_TEST ${PYTHONTESTS})
    set_tests_properties(${_TEST} PROPERTIES
        TIMEOUT 60 
        ENVIRONMENT "PYTHONPATH=${NATIVE_PYTHON_BINDINGS_DIR};${PATH_DEF}"
    )
endforeach()

if(INSTALL_TO_VIRTUALENV)
    find_program(VIRTUALENV_EXEC virtualenv)
    # Add an extra test ensuring the installed opencmiss iron package works in the virtualenvironment
    add_test(NAME Python_Bindings_Virtualenv
        COMMAND ${CMAKE_COMMAND}
        -DWD=${CMAKE_CURRENT_BINARY_DIR}
        -DVENV_EXEC=${VIRTUALENV_EXEC}
        -DPACKAGE_DIR=${NATIVE_CMAKE_CURRENT_BINARY_DIR}
        -P ${CMAKE_CURRENT_SOURCE_DIR}/TestVirtualenv.cmake
    )
    add_test_path(Python_Bindings_Virtualenv "${Iron_BINARY_DIR}" "${BINDINGS_BIN_DIR}")
endif()

