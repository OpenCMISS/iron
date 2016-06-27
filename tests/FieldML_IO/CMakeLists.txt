file(COPY input/ DESTINATION input)
add_executable(FieldML_IO fieldml_io.f90 cube.f90 fieldml_arguments.f90 ../Utilities/testframework.f90)
target_link_libraries(FieldML_IO iron)
oc_add_test(FieldML_IO FieldML_IO)
