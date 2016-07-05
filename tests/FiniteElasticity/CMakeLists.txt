add_executable(Cantilever Cantilever.f90)
set_target_properties(Cantilever PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries(Cantilever iron)
# This is reading the default settings that should apply to the file
# - however not giving it any arguments assumes type 9 interpolation opposite to the 
# DisplacementInterpolationType=CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION in line 150. Needs checking.
oc_add_test(FiniteElasticity_Cantilever Cantilever) # 3 2 2 1 3
foreach(ITYPE RANGE 1 9)
    set(TNAME "FiniteElasticity_Cantilever_I${ITYPE}_2_1_1")
    oc_add_test(${TNAME} Cantilever ${ITYPE} 2 1 1)
    set(TNAME "FiniteElasticity_Cantilever_I${ITYPE}_3_2_2")
    oc_add_test(${TNAME} Cantilever ${ITYPE} 3 2 2)
endforeach()

add_executable(SimpleShear SimpleShear.f90)
set_target_properties(SimpleShear PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries(SimpleShear iron)
oc_add_test(FiniteElasticity_SimpleShear SimpleShear)
