set(IRON_C_SRC
    #binary_file_c.c
    cmiss_c.c
    external_dae_solver_routines.c
    FieldExport.c
    timer_c.c
)
set(IRON_HEADERS
    external_dae_solver_routines.h
    FieldExport.h
    FieldExportConstants.h
    macros.h
    dllexport.h
)
set(IRON_Fortran_SRC
    advection_diffusion_equation_routines.F90
    advection_equation_routines.F90
    analytic_analysis_routines.F90
    base_routines.F90
    basis_routines.F90
    basis_access_routines.F90
    #binary_file_f.F90
    biodomain_equation_routines.F90
    bioelectric_finite_elasticity_routines.F90
    bioelectric_routines.F90
    blas.F90
    boundary_condition_routines.F90
    Burgers_equation_routines.F90
    cellml_access_routines.F90
    characteristic_equation_routines.F90
    classical_field_routines.F90
    cmiss_cellml.F90
    cmiss_fortran_c.F90
    cmiss_mpi.F90
    cmiss_parmetis.F90
    cmiss_petsc_types.F90
    cmiss_petsc.F90
    cmiss.F90
    computational_environment.F90
    constants.F90
    control_loop_routines.F90
    control_loop_access_routines.F90
    coordinate_routines.F90
    coordinate_access_routines.F90
    Darcy_equations_routines.F90
    Darcy_pressure_equations_routines.F90
    data_point_routines.F90
    data_point_access_routines.F90
    data_projection_routines.F90
    data_projection_access_routines.F90
    diffusion_advection_diffusion_routines.F90
    diffusion_diffusion_routines.F90
    diffusion_equation_routines.F90
    distributed_matrix_vector_IO.F90
    distributed_matrix_vector.F90
    distributed_matrix_vector_access.F90
    domain_mappings.F90
    elasticity_routines.F90
    electromechanics_routines.F90
    electrophysiology_cell_routines.F90
    equations_routines.F90
    equations_access_routines.F90
    equations_mapping_routines.F90
    equations_mapping_access_routines.F90
    equations_matrices_routines.F90
    equations_matrices_access_routines.F90
    equations_set_constants.F90
    equations_set_routines.F90
    equations_set_access_routines.F90
    field_IO_routines.F90
    field_routines.F90
    field_access_routines.F90
    finite_elasticity_Darcy_routines.F90
    finite_elasticity_fluid_pressure_routines.F90
    finite_elasticity_routines.F90
    fitting_routines.F90
    fluid_mechanics_IO_routines.F90
    fluid_mechanics_routines.F90
    fsi_routines.F90
    generated_mesh_routines.F90
    generated_mesh_access_routines.F90
    Hamilton_Jacobi_equations_routines.F90
    Helmholtz_equations_routines.F90
    #Helmholtz_TEMPLATE_equations_routines.F90
    history_routines.F90
    input_output.F90
    interface_conditions_constants.F90
    interface_conditions_routines.F90
    interface_condition_access_routines.F90
    interface_equations_routines.F90
    interface_equations_access_routines.F90
    interface_mapping_routines.F90
    interface_matrices_constants.F90
    interface_matrices_routines.F90
    interface_operators_routines.F90
    interface_routines.F90
    interface_access_routines.F90
    iso_varying_string.F90
    kinds.F90
    lapack.F90
    Laplace_equations_routines.F90
    linear_elasticity_routines.F90
    linkedlist_routines.F90
    lists.F90
    maths.F90
    matrix_vector.F90
    mesh_routines.F90
    mesh_access_routines.F90
    monodomain_equations_routines.F90
    multi_compartment_transport_routines.F90
    multi_physics_routines.F90
    Navier_Stokes_equations_routines.F90
    node_routines.F90
    opencmiss.F90
    opencmiss_iron.F90
    Poiseuille_equations_routines.F90
    Poisson_equations_routines.F90
    problem_constants.F90
    problem_routines.F90
    problem_access_routines.F90
    profiling_routines.F90
    reaction_diffusion_equation_routines.F90
    reaction_diffusion_IO_routines.F90
    region_routines.F90
    region_access_routines.F90
    solver_mapping_routines.F90
    solver_mapping_access_routines.F90
    solver_matrices_routines.F90
    solver_matrices_access_routines.F90
    solver_routines.F90
    solver_access_routines.F90
    sorting.F90
    Stokes_equations_routines.F90
    stree_equation_routines.F90
    strings.F90
    test_framework_routines.F90
    timer_f.F90
    trees.F90
    types.F90
    util_array.F90
)
# Add platform dependent files
IF(${OPERATING_SYSTEM} MATCHES linux)
    list(APPEND IRON_Fortran_SRC machine_constants_linux.F90)
    #list(INSERT IRON_Fortran_SRC 0 machine_constants_linux.F90)
ELSEIF(${OPERATING_SYSTEM} MATCHES darwin)
    list(APPEND IRON_Fortran_SRC machine_constants_linux.F90)
ELSEIF(${OPERATING_SYSTEM} MATCHES aix)
    list(APPEND IRON_Fortran_SRC machine_constants_aix.F90)
ELSE(${OPERATING_SYSTEM} MATCHES windows)
    list(APPEND IRON_Fortran_SRC machine_constants_win32.F90)
ENDIF()
#machine_constants_irix.F90
#machine_constants_vms.F90
    
# 
set(IRON_FIELDML_SRC)
if (WITH_FIELDML)
    list(APPEND IRON_Fortran_SRC
    #set(IRON_FIELDML_SRC
        fieldml_input_routines.F90
        fieldml_output_routines.F90
        fieldml_types.F90
        fieldml_util_routines.F90
    )
    #list(APPEND IRON_Fortran_SRC ${IRON_FIELDML_SRC})
endif()

# Fix paths to files
set(FIXPATH_VARS IRON_C_SRC IRON_Fortran_SRC)#IRON_FIELDML_SRC
foreach(varname ${FIXPATH_VARS})
    set(_TMP )
    foreach(filename ${${varname}})
        list(APPEND _TMP ${CMAKE_CURRENT_SOURCE_DIR}/src/${filename}) 
    endforeach()
    set(${varname} ${_TMP})
endforeach()

# Set combined sources variable
set(IRON_SRC ${IRON_C_SRC} ${IRON_Fortran_SRC})
