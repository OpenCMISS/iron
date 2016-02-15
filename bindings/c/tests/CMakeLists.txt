GET_LD_PATH(LD_PATH "${Iron_BINARY_DIR}/$<CONFIG>" "${BINDINGS_BIN_DIR}/$<CONFIG>")

function(oc_add_test TESTNAME)
    add_test(NAME ${TESTNAME} COMMAND ${ARGN})
    set_tests_properties(${TESTNAME} PROPERTIES
        ENVIRONMENT "${LD_VARNAME}=${LD_PATH}"
    )
endfunction()

add_executable(ComplexMesh complex_mesh.c)
target_link_libraries(ComplexMesh iron_c)
oc_add_test(C_Bindings_ComplexMesh ComplexMesh)

add_executable(Laplace_C laplace.c)
target_link_libraries(Laplace_C iron_c)
oc_add_test(C_Bindings_Laplace Laplace_C)
    
