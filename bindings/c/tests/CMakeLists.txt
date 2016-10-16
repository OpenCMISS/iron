function(oc_add_test TESTNAME)
    add_test(NAME ${TESTNAME} COMMAND ${ARGN})
    add_test_path(${TESTNAME} "${Iron_BINARY_DIR}" "${BINDINGS_BIN_DIR}")
endfunction()

add_executable(ComplexMesh complex_mesh.c)
set_target_properties(ComplexMesh PROPERTIES LINKER_LANGUAGE C)
target_link_libraries(ComplexMesh iron_c)
oc_add_test(C_Bindings_ComplexMesh ComplexMesh)

add_executable(Laplace_C laplace.c)
set_target_properties(Laplace_C PROPERTIES LINKER_LANGUAGE C)
target_link_libraries(Laplace_C iron_c)
oc_add_test(C_Bindings_Laplace Laplace_C)
    
