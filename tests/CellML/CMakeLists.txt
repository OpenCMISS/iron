add_executable(CellMLModelIntegration CellMLModelIntegration.f90)
target_link_libraries(CellMLModelIntegration iron)
oc_add_test(CellML_ModelIntegration CellMLModelIntegration ${CMAKE_CURRENT_SOURCE_DIR}/n98.xml)

add_executable(Monodomain Monodomain.f90)
target_link_libraries(Monodomain iron)
oc_add_test(CellML_Monodomain Monodomain 0.005 0.1001 70 ${CMAKE_CURRENT_SOURCE_DIR}/n98.xml)