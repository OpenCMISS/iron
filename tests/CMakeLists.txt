GET_LD_PATH(LD_PATH "${Iron_BINARY_DIR}/$<CONFIG>")
# Defines a simple wrapper for the add_test command that also sets the right environment for test execution
function(oc_add_test TESTNAME)
    add_test(NAME ${TESTNAME} COMMAND ${ARGN})
    set_tests_properties(${TESTNAME} PROPERTIES
        ENVIRONMENT "${LD_VARNAME}=${LD_PATH}"
    )
endfunction()

add_subdirectory(ClassicalField)
add_subdirectory(FiniteElasticity)
add_subdirectory(LinearElasticity)
add_subdirectory(FluidMechanics)
if (WITH_CELLML)
    add_subdirectory(CellML)
endif()
if (WITH_FIELDML)
    add_subdirectory(FieldML)
endif()