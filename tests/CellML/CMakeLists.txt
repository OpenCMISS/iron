add_executable(CellMLModelIntegration CellMLModelIntegration.f90)
target_link_libraries(CellMLModelIntegration iron)
add_test(CellML_ModelIntegration CellMLModelIntegration n98.xml)

add_executable(Monodomain Monodomain.f90)
target_link_libraries(Monodomain iron)
add_test(CellML_Monodomain Monodomain 0.005 0.1001 70 n98.xml)