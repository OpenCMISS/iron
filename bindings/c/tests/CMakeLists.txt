add_executable(ComplexMesh complex_mesh.c)
target_link_libraries(ComplexMesh iron_c)
add_test(C_Bindings_ComplexMesh ComplexMesh)

add_executable(Laplace_C laplace.c)
target_link_libraries(Laplace_C iron_c)
add_test(C_Bindings_Laplace Laplace_C)