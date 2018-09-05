!> \file
!> \author Chris Bradley
!> \brief This module handles all finite elasticity routines.
!>
!> \section LICENSE
!>
!> Version: MPL 1.1/GPL 2.0/LGPL 2.1
!>
!> The contents of this file are subject to the Mozilla Public License
!> Version 1.1 (the "License"); you may not use this file except in
!> compliance with the License. You may obtain a copy of the License at
!> http://www.mozilla.org/MPL/
!>
!> Software distributed under the License is distributed on an "AS IS"
!> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
!> License for the specific language governing rights and limitations
!> under the License.
!>
!> The Original Code is OpenCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand, the University of Oxford, Oxford, United
!> Kingdom and King's College, London, United Kingdom. Portions created
!> by the University of Auckland, the University of Oxford and King's
!> College, London are Copyright (C) 2007-2010 by the University of
!> Auckland, the University of Oxford and King's College, London.
!> All Rights Reserved.
!>
!> Contributor(s): Kumar Mithraratne, Jack Lee, Alice Hung, Sander Arens
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>

!>This module handles all finite elasticity routines.
MODULE FINITE_ELASTICITY_ROUTINES

  USE BaseRoutines
  USE BasisRoutines
  USE BasisAccessRoutines
  USE BOUNDARY_CONDITIONS_ROUTINES
  USE ComputationEnvironment
  USE Constants
  USE CONTROL_LOOP_ROUTINES
  USE ControlLoopAccessRoutines
  USE COORDINATE_ROUTINES  
  USE CoordinateSystemAccessRoutines
  USE DistributedMatrixVector
  USE DOMAIN_MAPPINGS
  USE EquationsRoutines
  USE EquationsAccessRoutines
  USE EquationsMappingRoutines
  USE EquationsMappingAccessRoutines
  USE EquationsMatricesRoutines
  USE EquationsSetConstants
  USE EquationsSetAccessRoutines
  USE FIELD_ROUTINES
  USE FieldAccessRoutines
  USE FIELD_IO_ROUTINES
  USE FLUID_MECHANICS_IO_ROUTINES
  USE GENERATED_MESH_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE Kinds
  USE Lapack
  USE Maths 
  USE MatrixVector
  USE MESH_ROUTINES
  USE MeshAccessRoutines
#ifndef NOMPIMOD
  USE MPI
#endif
  USE PROBLEM_CONSTANTS
  USE ProfilingRoutines
  USE SOLVER_ROUTINES
  USE SolverAccessRoutines
  USE SolverMappingAccessRoutines
  USE SolverMatricesAccessRoutines
  USE Strings
  USE Timer
  USE Types

#include "macros.h"

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif

  PRIVATE

  !Module parameters

  !> \addtogroup FINITE_ELASTICITY_ROUTINES_AnalyticParamIndices FINITE_ELASTICITY_ROUTINES::AnalyticParamIndices
  !> \brief Indices for EQUATIONS_SET_ANALYTIC_TYPE%ANALYTIC_USER_PARAMS
  !> \see FINITE_ELASTICITY_ROUTINES,OPENCMISS_AnalyticParamIndices
  !>@{
  INTEGER(INTG), PARAMETER :: FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_PIN_IDX=1 !<Inner pressure parameter index \see FINITE_ELASTICITY_ROUTINES_AnalyticParamIndices, FINITE_ELASTICITY_ROUTINES
  INTEGER(INTG), PARAMETER :: FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_POUT_IDX=2 !<Outer pressure parameter index \see FINITE_ELASTICITY_ROUTINES_AnalyticParamIndices, FINITE_ELASTICITY_ROUTINES
  INTEGER(INTG), PARAMETER :: FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_LAMBDA_IDX=3 !<Lambda parameter index \see FINITE_ELASTICITY_ROUTINES_AnalyticParamIndices, FINITE_ELASTICITY_ROUTINES
  INTEGER(INTG), PARAMETER :: FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_TSI_IDX=4 !<Tsi parameter index \see FINITE_ELASTICITY_ROUTINES_AnalyticParamIndices, FINITE_ELASTICITY_ROUTINES
  INTEGER(INTG), PARAMETER :: FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_RIN_IDX=5 !<Inner radius parameter index \see FINITE_ELASTICITY_ROUTINES_AnalyticParamIndices, FINITE_ELASTICITY_ROUTINES
  INTEGER(INTG), PARAMETER :: FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_ROUT_IDX=6 !<Outer radius parameter index \see FINITE_ELASTICITY_ROUTINES_AnalyticParamIndices, FINITE_ELASTICITY_ROUTINES
  INTEGER(INTG), PARAMETER :: FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_C1_IDX=7 !<c1 parameter index \see FINITE_ELASTICITY_ROUTINES_AnalyticParamIndices, FINITE_ELASTICITY_ROUTINES
  INTEGER(INTG), PARAMETER :: FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_C2_IDX=8 !<c2 parameter index \see FINITE_ELASTICITY_ROUTINES_AnalyticParamIndices, FINITE_ELASTICITY_ROUTINES
  !>@}

  !Module types

  !Module variables

  !Interfaces

  PUBLIC FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_PIN_IDX,FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_POUT_IDX, &
    & FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_LAMBDA_IDX, FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_TSI_IDX, &
    & FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_RIN_IDX, FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_ROUT_IDX, &
    & FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_C1_IDX, FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_C2_IDX

  PUBLIC FiniteElasticity_BoundaryConditionsAnalyticCalculate
  
  PUBLIC FiniteElasticity_FiniteElementResidualEvaluate

  PUBLIC FiniteElasticity_FiniteElementPreResidualEvaluate,FiniteElasticity_FiniteElementPostResidualEvaluate

  PUBLIC FiniteElasticity_FiniteElementJacobianEvaluate
  
  PUBLIC FINITE_ELASTICITY_EQUATIONS_SET_SETUP,FiniteElasticity_EquationsSetSolutionMethodSet, &
    & FiniteElasticity_EquationsSetSpecificationSet,FiniteElasticity_ProblemSpecificationSet,FINITE_ELASTICITY_PROBLEM_SETUP, &
    & FiniteElasticity_ContactProblemSpecificationSet,FiniteElasticity_ContactProblemSetup, & 
    & FiniteElasticity_PostSolve,FiniteElasticity_PostSolveOutputData, &
    & FiniteElasticity_PreSolve,FiniteElasticity_ControlTimeLoopPreLoop,FiniteElasticity_ControlLoadIncrementLoopPostLoop, &
    & EVALUATE_CHAPELLE_FUNCTION, GET_DARCY_FINITE_ELASTICITY_PARAMETERS, &
    & FiniteElasticity_GaussDeformationGradientTensor,FINITE_ELASTICITY_LOAD_INCREMENT_APPLY, &
    & FiniteElasticity_StressStrainCalculate
    
  PUBLIC FiniteElasticityEquationsSet_DerivedVariableCalculate
  
  PUBLIC FiniteElasticity_TensorInterpolateGaussPoint
  
  PUBLIC FiniteElasticity_TensorInterpolateXi

CONTAINS

  !
  !================================================================================================================================
  !

  !>Calculates the analytic solution and sets the boundary conditions for an analytic problem
  SUBROUTINE FiniteElasticity_BoundaryConditionsAnalyticCalculate(EQUATIONS_SET,BOUNDARY_CONDITIONS,err,error,*)
    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: node_idx,component_idx,deriv_idx,variable_idx,dim_idx,local_ny,variable_type
    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS,user_node,global_node,local_node
    REAL(DP) :: X(3),DEFORMED_X(3),P
    REAL(DP), POINTER :: GEOMETRIC_PARAMETERS(:)
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN,DOMAIN_PRESSURE
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES,DOMAIN_PRESSURE_NODES
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: NODES_MAPPING
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,GEOMETRIC_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE,GEOMETRIC_VARIABLE
    !BC stuff
    INTEGER(INTG),ALLOCATABLE :: INNER_SURFACE_NODES(:),OUTER_SURFACE_NODES(:),TOP_SURFACE_NODES(:),BOTTOM_SURFACE_NODES(:)
    INTEGER(INTG) :: INNER_NORMAL_XI,OUTER_NORMAL_XI,TOP_NORMAL_XI,BOTTOM_NORMAL_XI,MESH_COMPONENT
    INTEGER(INTG) :: myComputationalNodeNumber, DOMAIN_NUMBER, MPI_IERROR
    REAL(DP) :: PIN,POUT,LAMBDA,DEFORMED_Z
    LOGICAL :: X_FIXED,Y_FIXED,NODE_EXISTS, X_OKAY,Y_OKAY
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    NULLIFY(GEOMETRIC_PARAMETERS)

    ENTERS("FiniteElasticity_BoundaryConditionsAnalyticCalculate",err,error,*999)

    myComputationalNodeNumber=ComputationalEnvironment_NodeNumberGet(err,error)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
        DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
        IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
          GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
          IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN
            CALL FIELD_NUMBER_OF_COMPONENTS_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS,err,error,*999)
            !Get access to geometric coordinates
            NULLIFY(GEOMETRIC_VARIABLE)
            CALL Field_VariableGet(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,GEOMETRIC_VARIABLE,err,error,*999)
            MESH_COMPONENT=GEOMETRIC_VARIABLE%COMPONENTS(1)%MESH_COMPONENT_NUMBER
            CALL Field_ParameterSetDataGet(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS, &
              & err,error,*999)
            !Assign BC here - it's complicated so separate from analytic calculations
            IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
              DECOMPOSITION=>DEPENDENT_FIELD%DECOMPOSITION
              IF(ASSOCIATED(DECOMPOSITION)) THEN
                MESH=>DECOMPOSITION%MESH
                IF(ASSOCIATED(MESH)) THEN
                  GENERATED_MESH=>MESH%GENERATED_MESH
                  IF(ASSOCIATED(GENERATED_MESH)) THEN
                    NODES_MAPPING=>DECOMPOSITION%DOMAIN(1)%ptr%MAPPINGS%NODES   !HACK - ALL CHECKING INTERMEDIATE SKIPPED
                    IF(ASSOCIATED(NODES_MAPPING)) THEN
                      !Get surfaces (hardcoded): fix two nodes on the bottom face, pressure conditions inside & outside
                      CALL GENERATED_MESH_SURFACE_GET(GENERATED_MESH,MESH_COMPONENT,1_INTG, &
                          & INNER_SURFACE_NODES,INNER_NORMAL_XI,err,error,*999) !Inner
                      CALL GENERATED_MESH_SURFACE_GET(GENERATED_MESH,MESH_COMPONENT,2_INTG, &
                          & OUTER_SURFACE_NODES,OUTER_NORMAL_XI,err,error,*999) !Outer
                      CALL GENERATED_MESH_SURFACE_GET(GENERATED_MESH,MESH_COMPONENT,3_INTG, &
                          & TOP_SURFACE_NODES,TOP_NORMAL_XI,err,error,*999) !Top
                      CALL GENERATED_MESH_SURFACE_GET(GENERATED_MESH,MESH_COMPONENT,4_INTG, &
                          & BOTTOM_SURFACE_NODES,BOTTOM_NORMAL_XI,err,error,*999) !Bottom
                      !Set all inner surface nodes to inner pressure (- sign is to make positive P into a compressive force) ?
                      PIN=EQUATIONS_SET%ANALYTIC%ANALYTIC_USER_PARAMS(FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_PIN_IDX)
                      DO node_idx=1,SIZE(INNER_SURFACE_NODES,1)
                        user_node=INNER_SURFACE_NODES(node_idx)
                        !Need to test if this node is in current decomposition
                        CALL DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,user_node,1,DOMAIN_NUMBER,err,error,*999)
                        IF(DOMAIN_NUMBER==myComputationalNodeNumber) THEN
                          !Default to version 1 of each node derivative
                          CALL BOUNDARY_CONDITIONS_SET_NODE(BOUNDARY_CONDITIONS,DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1,1, &
                            & user_node,ABS(INNER_NORMAL_XI),BOUNDARY_CONDITION_PRESSURE_INCREMENTED,PIN,err,error,*999)
                        ENDIF
                      ENDDO
                      !Set all outer surface nodes to outer pressure
                      POUT=EQUATIONS_SET%ANALYTIC%ANALYTIC_USER_PARAMS(FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_POUT_IDX)
                      DO node_idx=1,SIZE(OUTER_SURFACE_NODES,1)
                        user_node=OUTER_SURFACE_NODES(node_idx)
                        !Need to test if this node is in current decomposition
                        CALL DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,user_node,1,DOMAIN_NUMBER,err,error,*999)
                        IF(DOMAIN_NUMBER==myComputationalNodeNumber) THEN
                          !Default to version 1 of each node derivative
                          CALL BOUNDARY_CONDITIONS_SET_NODE(BOUNDARY_CONDITIONS,DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1,1, &
                            & user_node,ABS(OUTER_NORMAL_XI),BOUNDARY_CONDITION_PRESSURE_INCREMENTED,POUT,err,error,*999)
                        ENDIF
                      ENDDO
                      !Set all top nodes fixed in z plane at lambda*height
                      LAMBDA=EQUATIONS_SET%ANALYTIC%ANALYTIC_USER_PARAMS(FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_LAMBDA_IDX)
                      DO node_idx=1,SIZE(TOP_SURFACE_NODES,1)
                        user_node=TOP_SURFACE_NODES(node_idx)
                        !Need to test if this node is in current decomposition
                        CALL DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,user_node,1,DOMAIN_NUMBER,err,error,*999)
                        IF(DOMAIN_NUMBER==myComputationalNodeNumber) THEN
                          CALL MeshTopology_NodeCheckExists(MESH,1,user_node,NODE_EXISTS,global_node,err,error,*999)
                          IF(.NOT.NODE_EXISTS) CYCLE
                          CALL DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET(NODES_MAPPING,global_node,NODE_EXISTS,local_node,err,error,*999)
                          !Default to version 1 of each node derivative
                          local_ny=GEOMETRIC_VARIABLE%COMPONENTS(3)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(local_node)% &
                            & DERIVATIVES(1)%VERSIONS(1)
                          DEFORMED_Z=GEOMETRIC_PARAMETERS(local_ny)*LAMBDA
                          !Default to version 1 of each node derivative
                          CALL BOUNDARY_CONDITIONS_SET_NODE(BOUNDARY_CONDITIONS,DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,1,1, &
                            & user_node,ABS(TOP_NORMAL_XI),BOUNDARY_CONDITION_FIXED,DEFORMED_Z,err,error,*999)
                        ENDIF
                      ENDDO
                      !Set all bottom nodes fixed in z plane
                      DO node_idx=1,SIZE(BOTTOM_SURFACE_NODES,1)
                        user_node=BOTTOM_SURFACE_NODES(node_idx)
                        !Need to check this node exists in the current domain
                        CALL DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,user_node,1,DOMAIN_NUMBER,err,error,*999)
                        IF(DOMAIN_NUMBER==myComputationalNodeNumber) THEN
                          !Default to version 1 of each node derivative
                          CALL BOUNDARY_CONDITIONS_SET_NODE(BOUNDARY_CONDITIONS,DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,1,1, &
                            & user_node,ABS(BOTTOM_NORMAL_XI),BOUNDARY_CONDITION_FIXED,0.0_DP,err,error,*999)
                        ENDIF
                      ENDDO
                      !Set two nodes on the bottom surface to axial displacement only:
                      !Easier for parallel: Fix everything that can be fixed !!!
                      X_FIXED=.FALSE.
                      Y_FIXED=.FALSE.
                      DO node_idx=1,SIZE(BOTTOM_SURFACE_NODES,1)
                        user_node=BOTTOM_SURFACE_NODES(node_idx)
                        CALL DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,user_node,1,DOMAIN_NUMBER,err,error,*999)
                        IF(DOMAIN_NUMBER==myComputationalNodeNumber) THEN
                          CALL MeshTopology_NodeCheckExists(MESH,1,user_node,NODE_EXISTS,global_node,err,error,*999)
                          IF(.NOT.NODE_EXISTS) CYCLE
                          CALL DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET(NODES_MAPPING,global_node,NODE_EXISTS,local_node,err,error,*999)
                          !Default to version 1 of each node derivative
                          local_ny=GEOMETRIC_VARIABLE%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(local_node)% &
                            & DERIVATIVES(1)%VERSIONS(1)
                          X(1)=GEOMETRIC_PARAMETERS(local_ny)
                            CALL MeshTopology_NodeCheckExists(MESH,1,user_node,NODE_EXISTS,global_node,err,error,*999)
                            IF(.NOT.NODE_EXISTS) CYCLE
                            CALL DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET(NODES_MAPPING,global_node,NODE_EXISTS,local_node, &
                              & err,error,*999)
                            !Default to version 1 of each node derivative
                            local_ny=GEOMETRIC_VARIABLE%COMPONENTS(2)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(local_node)% &
                            & DERIVATIVES(1)%VERSIONS(1)
                          X(2)=GEOMETRIC_PARAMETERS(local_ny)
                          IF(ABS(X(1))<1E-7_DP) THEN
                            !Default to version 1 of each node derivative
                            CALL BOUNDARY_CONDITIONS_SET_NODE(BOUNDARY_CONDITIONS,DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,1,1, &
                              & user_node,1,BOUNDARY_CONDITION_FIXED,0.0_DP,err,error,*999)

                            X_FIXED=.TRUE.
                          ENDIF
                          IF(ABS(X(2))<1E-7_DP) THEN
                            !Default to version 1 of each node derivative
                            CALL BOUNDARY_CONDITIONS_SET_NODE(BOUNDARY_CONDITIONS,DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,1,1, &
                              & user_node,2,BOUNDARY_CONDITION_FIXED,0.0_DP,err,error,*999)

                            Y_FIXED=.TRUE.
                          ENDIF
                        ENDIF
                      ENDDO
                      !Check it went well
                      CALL MPI_REDUCE(X_FIXED,X_OKAY,1,MPI_LOGICAL,MPI_LOR,0,computationalEnvironment%mpiCommunicator,MPI_IERROR)
                      CALL MPI_REDUCE(Y_FIXED,Y_OKAY,1,MPI_LOGICAL,MPI_LOR,0,computationalEnvironment%mpiCommunicator,MPI_IERROR)
                      IF(myComputationalNodeNumber==0) THEN
                        IF(.NOT.(X_OKAY.AND.Y_OKAY)) THEN
                          CALL FlagError("Could not fix nodes to prevent rigid body motion",err,error,*999)
                        ENDIF
                      ENDIF
                    ELSE
                      CALL FlagError("Domain nodes mapping is not associated.",err,error,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Generated mesh is not associated. For the Cylinder analytic solution, "// &
                      & "it must be available for automatic boundary condition assignment",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Mesh is not associated",err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("Decomposition is not associated",err,error,*999)
              ENDIF

              !Now calculate analytic solution
              DO variable_idx=1,DEPENDENT_FIELD%NUMBER_OF_VARIABLES
                variable_type=DEPENDENT_FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE
                FIELD_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%ptr
                IF(variable_idx==1) CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Global number of dofs : ", &
                  & FIELD_VARIABLE%NUMBER_OF_GLOBAL_DOFS,err,error,*999)
                IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                  CALL FIELD_PARAMETER_SET_CREATE(DEPENDENT_FIELD,variable_type,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
                  component_idx=1 !Assuming components 1..3 use a common mesh component and 4 uses a different one
                  IF(FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN
                    DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
                    IF(ASSOCIATED(DOMAIN)) THEN
                      IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                        DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                        IF(ASSOCIATED(DOMAIN_NODES)) THEN
                          !Also grab the equivalent pointer for pressure component
                          IF(FIELD_VARIABLE%COMPONENTS(4)%INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN
                            DOMAIN_PRESSURE=>FIELD_VARIABLE%COMPONENTS(4)%DOMAIN
                            IF(ASSOCIATED(DOMAIN_PRESSURE)) THEN
                              IF(ASSOCIATED(DOMAIN_PRESSURE%TOPOLOGY)) THEN
                                DOMAIN_PRESSURE_NODES=>DOMAIN_PRESSURE%TOPOLOGY%NODES
                                  IF(ASSOCIATED(DOMAIN_PRESSURE_NODES)) THEN

                                  !Loop over the local nodes excluding the ghosts.
                                  DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
                                    !!TODO \todo We should interpolate the geometric field here and the node position.
                                    DO dim_idx=1,NUMBER_OF_DIMENSIONS
                                      !Default to version 1 of each node derivative
                                      local_ny=GEOMETRIC_VARIABLE%COMPONENTS(dim_idx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                                        & NODES(node_idx)%DERIVATIVES(1)%VERSIONS(1)
                                      X(dim_idx)=GEOMETRIC_PARAMETERS(local_ny)
                                    ENDDO !dim_idx
                                    !Loop over the derivatives
                                    DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                                      SELECT CASE(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE)
                                      CASE(EQUATIONS_SET_FINITE_ELASTICITY_CYLINDER)
                                        !Cylinder inflation, extension, torsion
                                        SELECT CASE(variable_type)
                                        CASE(FIELD_U_VARIABLE_TYPE)
                                          SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                          CASE(NO_GLOBAL_DERIV)
                                            !Do all components at the same time (r,theta,z)->(x,y,z) & p
                                            CALL FiniteElasticity_CylinderAnalyticCalculate(X, &
                                              & EQUATIONS_SET%ANALYTIC%ANALYTIC_USER_PARAMS,DEFORMED_X,P,err,error,*999)
                                          CASE(GLOBAL_DERIV_S1)
                                            CALL FlagError("Not implemented.",err,error,*999)
                                          CASE(GLOBAL_DERIV_S2)
                                            CALL FlagError("Not implemented.",err,error,*999)
                                          CASE(GLOBAL_DERIV_S1_S2)
                                            CALL FlagError("Not implemented.",err,error,*999)
                                          CASE DEFAULT
                                            LOCAL_ERROR="The global derivative index of "//TRIM(NumberToVString( &
                                              DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                              & err,error))//" is invalid."
                                            CALL FlagError(LOCAL_ERROR,err,error,*999)
                                          END SELECT
                                        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                                          SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                          CASE(NO_GLOBAL_DERIV)
                                            !Not implemented, but don't want to cause an error so do nothing
                                          CASE(GLOBAL_DERIV_S1)
                                            CALL FlagError("Not implemented.",err,error,*999)
                                          CASE(GLOBAL_DERIV_S2)
                                            CALL FlagError("Not implemented.",err,error,*999)
                                          CASE(GLOBAL_DERIV_S1_S2)
                                            CALL FlagError("Not implemented.",err,error,*999)
                                          CASE DEFAULT
                                            LOCAL_ERROR="The global derivative index of "//TRIM(NumberToVString( &
                                              DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                                & err,error))//" is invalid."
                                            CALL FlagError(LOCAL_ERROR,err,error,*999)
                                          END SELECT
                                        CASE DEFAULT
                                          LOCAL_ERROR="The variable type "//TRIM(NumberToVString(variable_type,"*",err,error)) &
                                            & //" is invalid."
                                          CALL FlagError(LOCAL_ERROR,err,error,*999)
                                        END SELECT
                                      CASE DEFAULT
                                        LOCAL_ERROR="The analytic function type of "// &
                                          & TRIM(NumberToVString(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
                                          & " is invalid."
                                        CALL FlagError(LOCAL_ERROR,err,error,*999)
                                      END SELECT
                                      !Set the analytic solution to parameter set
                                      DO component_idx=1,NUMBER_OF_DIMENSIONS
                                        !Default to version 1 of each node derivative
                                        local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                          & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                        CALL Field_ParameterSetUpdateLocalDOF(DEPENDENT_FIELD,variable_type, &
                                          & FIELD_ANALYTIC_VALUES_SET_TYPE,local_ny,DEFORMED_X(component_idx),err,error,*999)
                                      ENDDO
                                      !Don't forget the pressure component
                                      user_node=DOMAIN_NODES%NODES(node_idx)%USER_NUMBER
                                      CALL MeshTopology_NodeCheckExists(MESH,DOMAIN_PRESSURE%MESH_COMPONENT_NUMBER,user_node, &
                                        & NODE_EXISTS,global_node,err,error,*999)
                                      IF(NODE_EXISTS) THEN
                                        CALL DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,user_node, &
                                          & DOMAIN_PRESSURE%MESH_COMPONENT_NUMBER,DOMAIN_NUMBER,err,error,*999)
                                        IF(DOMAIN_NUMBER==myComputationalNodeNumber) THEN
                                          !\todo: test the domain node mappings pointer properly
                                          local_node=DOMAIN_PRESSURE%mappings%nodes%global_to_local_map(global_node)%local_number(1)
                                          !Default to version 1 of each node derivative
                                          local_ny=FIELD_VARIABLE%COMPONENTS(4)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                                            & NODES(local_node)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                          !Because p=2.lambda in this particular constitutive law, we'll assign half the
                                          !hydrostatic pressure to the analytic array
                                          CALL Field_ParameterSetUpdateLocalDOF(DEPENDENT_FIELD,variable_type, &
                                          & FIELD_ANALYTIC_VALUES_SET_TYPE,local_ny,P/2.0_dp,err,error,*999)
                                        ENDIF
                                      ENDIF
                                    ENDDO !deriv_idx
                                  ENDDO !node_idx

                                ELSE
                                  CALL FlagError("Domain for pressure topology node is not associated",err,error,*999)
                                ENDIF
                              ELSE
                                CALL FlagError("Domain for pressure topology is not associated",err,error,*999)
                              ENDIF
                            ELSE
                              CALL FlagError("Domain for pressure component is not associated",err,error,*999)
                            ENDIF
                          ELSE
                            CALL FlagError("Non-nodal based interpolation of pressure cannot be used with analytic solutions", &
                              & err,error,*999)
                          ENDIF
                        ELSE
                          CALL FlagError("Domain topology nodes is not associated.",err,error,*999)
                        ENDIF
                      ELSE
                        CALL FlagError("Domain topology is not associated.",err,error,*999)
                      ENDIF
                    ELSE
                      CALL FlagError("Domain is not associated.",err,error,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Only node based interpolation is implemented.",err,error,*999)
                  ENDIF
                  CALL Field_ParameterSetUpdateStart(DEPENDENT_FIELD,variable_type,FIELD_ANALYTIC_VALUES_SET_TYPE, &
                    & err,error,*999)
                  CALL Field_ParameterSetUpdateFinish(DEPENDENT_FIELD,variable_type,FIELD_ANALYTIC_VALUES_SET_TYPE, &
                    & err,error,*999)
                ELSE
                  CALL FlagError("Field variable is not associated.",err,error,*999)
                ENDIF

              ENDDO !variable_idx
              CALL Field_ParameterSetDataRestore(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & GEOMETRIC_PARAMETERS,err,error,*999)
            ELSE
              CALL FlagError("Boundary conditions is not associated.",err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Equations set geometric field is not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Equations set dependent field is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Equations set analytic is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF


    EXITS("FiniteElasticity_BoundaryConditionsAnalyticCalculate")
    RETURN
999 ERRORS("FiniteElasticity_BoundaryConditionsAnalyticCalculate",err,error)
    EXITS("FiniteElasticity_BoundaryConditionsAnalyticCalculate")
    RETURN 1

  END SUBROUTINE FiniteElasticity_BoundaryConditionsAnalyticCalculate

  !
  !================================================================================================================================
  !

  !>Calcualates the analytic solution (deformed coordinates and hydrostatic pressure) for cylinder inflation+extension+torsion problem
  SUBROUTINE FiniteElasticity_CylinderAnalyticCalculate(X,ANALYTIC_USER_PARAMS,DEFORMED_X,P,err,error,*)
    !Argument variables
    REAL(DP), INTENT(IN) :: X(:)                !<Undeformed coordinates
    REAL(DP), INTENT(IN) :: ANALYTIC_USER_PARAMS(:) !<Array containing the problem parameters
    REAL(DP), INTENT(OUT) :: DEFORMED_X(3)      !<Deformed coordinates
    REAL(DP), INTENT(OUT) :: P                  !<Hydrostatic pressure at the given material coordintae
    INTEGER(INTG), INTENT(OUT) :: ERR           !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR  !<The error string
    !Local variables
    REAL(DP) :: PIN,POUT,LAMBDA,TSI,A1,A2,C1,C2 !A1=external radius, A2=internal radius
    REAL(DP) :: MU1,MU2,MU,K
    REAL(DP) :: F,F2,DF
    REAL(DP) :: R,THETA ! Undeformed coordinates in radial coordinates
    REAL(DP) :: DEFORMED_R,DEFORMED_THETA
    REAL(DP) :: DELTA,RES
    REAL(DP), PARAMETER :: STEP=1E-5_DP, RELTOL=1E-12_DP


    ENTERS("FiniteElasticity_CylinderAnalyticCalculate",err,error,*999)

    !Grab problem parameters
    PIN=ANALYTIC_USER_PARAMS(FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_PIN_IDX)
    POUT=ANALYTIC_USER_PARAMS(FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_POUT_IDX)
    LAMBDA=ANALYTIC_USER_PARAMS(FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_LAMBDA_IDX)
    TSI=ANALYTIC_USER_PARAMS(FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_TSI_IDX)
    A1=ANALYTIC_USER_PARAMS(FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_ROUT_IDX) ! external radius
    A2=ANALYTIC_USER_PARAMS(FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_RIN_IDX) ! internal radius
    C1=ANALYTIC_USER_PARAMS(FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_C1_IDX)
    C2=ANALYTIC_USER_PARAMS(FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_C2_IDX)

    !Solve for MU1 - Newton's method (\todo: Implement here, or separate out for general use?)
    MU1=1.0_DP  !Initial guess - need a better way!
    DO
      !Calculate f(MU1)
      F=FINITE_ELASTICITY_CYLINDER_ANALYTIC_FUNC_EVALUATE(MU1,PIN,POUT,LAMBDA,TSI,A1,A2,C1,C2,err,error)
      IF(ERR/=0) GOTO 999

      !Calculate f'(MU1) by finite differencing
      F2=FINITE_ELASTICITY_CYLINDER_ANALYTIC_FUNC_EVALUATE(MU1+STEP,PIN,POUT,LAMBDA,TSI,A1,A2,C1,C2,err,error)
      IF(ERR/=0) GOTO 999
      DF=(F2-F)/STEP

      !Next increment for MU1
      DELTA=-F/DF

      !Ensure that the step actually reduces residual
      F2=FINITE_ELASTICITY_CYLINDER_ANALYTIC_FUNC_EVALUATE(MU1+DELTA,PIN,POUT,LAMBDA,TSI,A1,A2,C1,C2,err,error)
      IF(ERR/=0) GOTO 999
      DO
        IF (ABS(F2)<ABS(F).OR.ABS(F2)<ZERO_TOLERANCE) THEN    ! PASS
          MU1=MU1+DELTA
          EXIT
        ELSEIF (DELTA<1E-3_DP) THEN ! FAIL: It's likely that the initial guess is too far away
          CALL FlagError("FiniteElasticity_CylinderAnalyticCalculate failed to converge.",err,error,*999)
        ELSE                        ! KEEP GOING
          DELTA=DELTA/2.0_DP
          F2=FINITE_ELASTICITY_CYLINDER_ANALYTIC_FUNC_EVALUATE(MU1+DELTA,PIN,POUT,LAMBDA,TSI,A1,A2,C1,C2,err,error)
          IF(ERR/=0) GOTO 999
        ENDIF
      ENDDO

      !Test for convergence: relative residual
      RES=DELTA/(1.0_DP+MU1)
      IF (RES<RELTOL) EXIT
    ENDDO

    !Calculate MU2
    MU2=SQRT(((A1/A2)**2*(LAMBDA*MU1**2-1.0_DP)+1.0_DP)/LAMBDA)

    !Calculate radius and angle from undeformed coordinates
    R=SQRT(X(1)**2+X(2)**2)
    THETA=ATAN2(X(2),X(1)) ! in radians

    !Calculate deformed coordinates
    K=A1**2*(LAMBDA*MU1**2-1.0_DP)
    MU=SQRT(1.0_DP/LAMBDA*(1.0_DP+K/R**2))
    DEFORMED_R=MU*R
    DEFORMED_THETA=THETA+TSI*LAMBDA*X(3)
    DEFORMED_X(1)=DEFORMED_R*COS(DEFORMED_THETA)
    DEFORMED_X(2)=DEFORMED_R*SIN(DEFORMED_THETA)
    DEFORMED_X(3)=LAMBDA*X(3)

    !Calculate pressure
    P=POUT-(C1/LAMBDA+C2*LAMBDA)*(1.0_DP/LAMBDA/MU1**2-R**2/(R**2+K)+LOG(MU**2/MU1**2))+C1*TSI**2*LAMBDA*(R**2-A1**2) &
      & -2.0_DP*(C1/LAMBDA**2/MU**2+C2*(1.0_DP/LAMBDA**2+1.0_DP/MU**2+TSI**2*R**2))

    EXITS("FiniteElasticity_CylinderAnalyticCalculate")
    RETURN
999 ERRORSEXITS("FiniteElasticity_CylinderAnalyticCalculate",err,error)
    RETURN 1

  END SUBROUTINE FiniteElasticity_CylinderAnalyticCalculate

  !
  !================================================================================================================================
  !

  !>Evaluates the residual function required to solve for MU1, in the cylinder analytic example
  FUNCTION FINITE_ELASTICITY_CYLINDER_ANALYTIC_FUNC_EVALUATE(MU1,PIN,POUT,LAMBDA,TSI,A1,A2,C1,C2,err,error)
    !Argument variables
    REAL(DP) :: FINITE_ELASTICITY_CYLINDER_ANALYTIC_FUNC_EVALUATE
    REAL(DP) :: MU1,PIN,POUT,LAMBDA,TSI,A1,A2,C1,C2
    INTEGER(INTG) :: ERR
    TYPE(VARYING_STRING) :: ERROR
    !Local variables
    REAL(DP) :: MU,K

    ENTERS("FINITE_ELASTICITY_CYLINDER_ANALYTIC_FUNC_EVALUATE",err,error,*999)

    K=A1**2*(LAMBDA*MU1**2-1.0_DP)
    MU=SQRT(1.0_DP/LAMBDA*(1.0_DP+K/A2**2))

    FINITE_ELASTICITY_CYLINDER_ANALYTIC_FUNC_EVALUATE= &
      &  2.0_DP*(C1/LAMBDA**2/MU**2 + C2*(1.0_DP/LAMBDA**2+1.0_DP/MU**2+TSI**2*A2**2))+ &
      & POUT-(C1/LAMBDA+C2*LAMBDA)*(1.0_DP/LAMBDA/MU1**2-A2**2/(A2**2+K)+2*LOG(MU/MU1))+ &
      & C1*TSI**2*LAMBDA*(A2**2-A1**2)-2.0_DP*(C1/LAMBDA**2/MU**2+C2*(1.0_DP/LAMBDA**2+ &
      & 1.0_DP/MU**2+TSI**2*A2**2))+PIN
  
    EXITS("FINITE_ELASTICITY_CYLINDER_ANALYTIC_FUNC_EVALUATE")
    RETURN
999 ERRORS("FINITE_ELASTICITY_CYLINDER_ANALYTIC_FUNC_EVALUATE",err,error)
    EXITS("FINITE_ELASTICITY_CYLINDER_ANALYTIC_FUNC_EVALUATE")
    RETURN
    
  END FUNCTION FINITE_ELASTICITY_CYLINDER_ANALYTIC_FUNC_EVALUATE

  !
  !================================================================================================================================
  !

  !>Evaluates the spatial elasticity and stress tensor in Voigt form at a given Gauss point.
  SUBROUTINE FINITE_ELASTICITY_GAUSS_ELASTICITY_TENSOR(EQUATIONS_SET,DEPENDENT_INTERPOLATED_POINT, &
      & MATERIALS_INTERPOLATED_POINT,ELASTICITY_TENSOR,HYDRO_ELASTICITY_VOIGT,STRESS_TENSOR,DZDNU, &
      & Jznu,ELEMENT_NUMBER,GAUSS_POINT_NUMBER,ERR,ERROR,*)
    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER, INTENT(IN) :: EQUATIONS_SET !<A pointer to the equations set
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: DEPENDENT_INTERPOLATED_POINT,MATERIALS_INTERPOLATED_POINT
    REAL(DP), INTENT(OUT) :: ELASTICITY_TENSOR(:,:) !< Rank 4 elasticity tensor in Voigt notation
    REAL(DP), INTENT(OUT) :: HYDRO_ELASTICITY_VOIGT(:) !<Rank 2 hydrostatic portion of the elasticity tensor in Voigt notation
    REAL(DP), INTENT(OUT) :: STRESS_TENSOR(:) !< Rank 2 stress tensor in Voigt notation
    REAL(DP), INTENT(IN) :: DZDNU(:,:)!< The deformation gradient
    REAL(DP), INTENT(IN) :: Jznu !< The Jacobian
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER,GAUSS_POINT_NUMBER !<Element/Gauss point number
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: PRESSURE_COMPONENT,i,j,dof_idx
    REAL(DP) :: P, I1, I3
    REAL(DP) :: DZDNUT(3,3),AZL(3,3),AZU(3,3),TEMP(3,3)
    REAL(DP) :: AZLv(6), AZUv(6) !<Voigt forms of the C and C^-1 tensors.
    REAL(DP) :: TEMPTERM1,TEMPTERM2,VALUE
    REAL(DP), POINTER :: C(:) !Parameters for constitutive laws
    REAL(DP) :: B(6),E(6),DQ_DE(6),Q
    REAL(DP) :: I3EE(6,6) !<Derivative of I3 wrt E
    REAL(DP) :: ADJCC(6,6) !<Derivative of adj(C) wrt C
    REAL(DP) :: AZUE(6,6) !<Derivative of C^-1 wrt E
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("FINITE_ELASTICITY_GAUSS_ELASTICITY_TENSOR",ERR,ERROR,*999)

    NULLIFY(FIELD_VARIABLE,C)

    !AZL = F'*F (deformed covariant or right cauchy deformation tensor, C)
    !AZU - deformed contravariant tensor; I3 = det(C)
    !E = Green-Lagrange strain tensor = 0.5*(C-I)
    !P is the hydrostatic pressure

    ! Evaluate the Cauchy strain tensor C.
    CALL MatrixTranspose(DZDNU,DZDNUT,ERR,ERROR,*999)
    CALL MatrixProduct(DZDNUT,DZDNU,AZL,ERR,ERROR,*999)
    CALL Invert(AZL,AZU,I3,ERR,ERROR,*999)

    ! Evaluate the derivative of AZU wrt to E (AZUE) for the hydrostatic term. Formulation from Nam-Ho Kim book, pg.198.
    AZLv(1) = AZL(1,1)
    AZLv(2) = AZL(2,2)
    AZLv(3) = AZL(3,3)
    AZLv(4) = AZL(1,2)
    AZLv(5) = AZL(1,3)
    AZLv(6) = AZL(2,3)
    AZUv(1) = AZU(1,1)
    AZUv(2) = AZU(2,2)
    AZUv(3) = AZU(3,3)
    AZUv(4) = AZU(1,2)
    AZUv(5) = AZU(1,3)
    AZUv(6) = AZU(2,3)
    I3EE = RESHAPE([0.0_DP, 4.0_DP*AZLv(3), 4.0_DP*AZLv(2), 0.0_DP,  0.0_DP,-4.0_DP*AZLv(6), &
      & 4.0_DP*AZLv(3), 0.0_DP, 4.0_DP*AZLv(1), 0.0_DP,-4.0_DP*AZLv(5), 0.0_DP,  &
      & 4.0_DP*AZLv(2), 4.0_DP*AZLv(1), 0.0_DP, -2.0_DP*AZLv(4), 0.0_DP, 0.0_DP, &
      & 0.0_DP, 0.0_DP, -4.0_DP*AZLv(4), -2.0_DP*AZLv(3), 2.0_DP*AZLv(6), 2.0_DP*AZLv(5), &
      & 0.0_DP, -4.0_DP*AZLv(5), 0.0_DP, 2.0_DP*AZLv(6), -2.0_DP*AZLv(2), 2.0_DP*AZLv(4), &
      & -4.0_DP*AZLv(6), 0.0_DP, 0.0_DP, 2.0_DP*AZLv(5), 2.0_DP*AZLv(4), -2.0_DP*AZLv(1)], [6,6])
    ADJCC = RESHAPE([0.0_DP, AZLv(3), AZLv(2), 0.0_DP,  0.0_DP,-AZLv(6), &
      & AZLv(3), 0.0_DP, AZLv(1), 0.0_DP,-AZLv(5), 0.0_DP,  &
      & AZLv(2), AZLv(1), 0.0_DP, -AZLv(4), 0.0_DP, 0.0_DP, &
      & 0.0_DP, 0.0_DP, -AZLv(4), -0.5_DP*AZLv(3), 0.5_DP*AZLv(6), 0.5_DP*AZLv(5), &
      & 0.0_DP, -AZLv(5), 0.0_DP,0.5_DP*AZLv(6), -0.5_DP*AZLv(2), 0.5_DP*AZLv(4), &
      & -AZLv(6), 0.0_DP, 0.0_DP, 0.5_DP*AZLv(5), 0.5_DP*AZLv(4), -0.5_DP*AZLv(1)], [6,6])
    !DO i=1,6
    !  DO j=1,6
    !    AZUE(i,j) = -2.0_DP*AZUv(i)*AZUv(j) + 2.0_DP*ADJCC(i,j)/I3
    !  ENDDO
    !ENDDO

    DO i=1,6
      DO j=1,6
        AZUE(i,j) = -2.0_DP*AZUv(i)*AZUv(j) + 0.5_DP*I3EE(i,j)/I3
      ENDDO
    ENDDO

    C=>MATERIALS_INTERPOLATED_POINT%VALUES(:,NO_PART_DERIV)

    ELASTICITY_TENSOR=0.0_DP

    SELECT CASE(EQUATIONS_SET%specification(3))
    CASE(EQUATIONS_SET_MOONEY_RIVLIN_ACTIVECONTRACTION_SUBTYPE, &
      & EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE, &
      & EQUATIONS_SET_MR_AND_GROWTH_LAW_IN_CELLML_SUBTYPE)
      LOCAL_ERROR="Analytic Jacobian has not been validated for the Mooney-Rivlin equations, please use finite differences instead."
      CALL FlagWarning(LOCAL_ERROR,ERR,ERROR,*999)
      PRESSURE_COMPONENT=DEPENDENT_INTERPOLATED_POINT%INTERPOLATION_PARAMETERS%FIELD_VARIABLE%NUMBER_OF_COMPONENTS
      P=DEPENDENT_INTERPOLATED_POINT%VALUES(PRESSURE_COMPONENT,NO_PART_DERIV)
      !Form of constitutive model is:
      ! W=c1*(I1-3)+c2*(I2-3)+p/2*(I3-1)

      ! Calculate isochoric fictitious 2nd Piola tensor (in Voigt form)
      I1=AZL(1,1)+AZL(2,2)+AZL(3,3)
      TEMPTERM1=-2.0_DP*C(2)
      TEMPTERM2=2.0_DP*(C(1)+I1*C(2))
      STRESS_TENSOR(1)=TEMPTERM1*AZL(1,1)+TEMPTERM2
      STRESS_TENSOR(2)=TEMPTERM1*AZL(2,2)+TEMPTERM2
      STRESS_TENSOR(3)=TEMPTERM1*AZL(3,3)+TEMPTERM2
      STRESS_TENSOR(4)=TEMPTERM1*AZL(2,1)
      STRESS_TENSOR(5)=TEMPTERM1*AZL(3,1)
      STRESS_TENSOR(6)=TEMPTERM1*AZL(3,2)
      IF(EQUATIONS_SET%specification(3)==EQUATIONS_SET_MOONEY_RIVLIN_ACTIVECONTRACTION_SUBTYPE) THEN
        
        !add active contraction stress values
        !Be aware for modified DZDNU, should active contraction be added here? Normally should be okay as modified DZDNU and DZDNU
        !converge during the Newton iteration.
        CALL Field_VariableGet(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VARIABLE,ERR,ERROR,*999)
        DO i=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
          dof_idx=FIELD_VARIABLE%COMPONENTS(i)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP% &
            & GAUSS_POINTS(GAUSS_POINT_NUMBER,ELEMENT_NUMBER)
          CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,dof_idx,VALUE,ERR,ERROR,*999)
          STRESS_TENSOR(i)=STRESS_TENSOR(i)+VALUE
        ENDDO
      ENDIF

      ! Calculate material elasticity tensor (in Voigt form) as
      ! this will be compensated for in the push-forward with the modified deformation gradient.
      TEMPTERM1=4.0_DP*C(2)
      TEMPTERM2=-2.0_DP*C(2)
      ELASTICITY_TENSOR(2,1)=TEMPTERM1
      ELASTICITY_TENSOR(3,1)=TEMPTERM1
      ELASTICITY_TENSOR(1,2)=TEMPTERM1
      ELASTICITY_TENSOR(3,2)=TEMPTERM1
      ELASTICITY_TENSOR(1,3)=TEMPTERM1
      ELASTICITY_TENSOR(2,3)=TEMPTERM1
      ELASTICITY_TENSOR(4,4)=TEMPTERM2
      ELASTICITY_TENSOR(5,5)=TEMPTERM2
      ELASTICITY_TENSOR(6,6)=TEMPTERM2
      !Add volumetric part of elasticity tensor - p*d(C^-1)/dE.
      ELASTICITY_TENSOR=ELASTICITY_TENSOR + P*AZUE

      !Hydrostatic portion of the elasticity tensor (dS/dp)
      HYDRO_ELASTICITY_VOIGT = AZUv

      ! Do push-forward of 2nd Piola tensor and the material elasticity tensor.
      CALL FINITE_ELASTICITY_PUSH_STRESS_TENSOR(STRESS_TENSOR,DZDNU,Jznu,ERR,ERROR,*999)
      CALL FINITE_ELASTICITY_PUSH_STRESS_TENSOR(HYDRO_ELASTICITY_VOIGT,DZDNU,Jznu,ERR,ERROR,*999)
      CALL FINITE_ELASTICITY_PUSH_ELASTICITY_TENSOR(ELASTICITY_TENSOR,DZDNU,Jznu,ERR,ERROR,*999)

      ! Add volumetric parts.
      STRESS_TENSOR(1:3)=STRESS_TENSOR(1:3)+P

    CASE(EQUATIONS_SET_TRANSVERSE_ISOTROPIC_GUCCIONE_SUBTYPE,EQUATIONS_SET_GUCCIONE_ACTIVECONTRACTION_SUBTYPE, &
      & EQUATIONS_SET_REFERENCE_STATE_TRANSVERSE_GUCCIONE_SUBTYPE)
      PRESSURE_COMPONENT=DEPENDENT_INTERPOLATED_POINT%INTERPOLATION_PARAMETERS%FIELD_VARIABLE%NUMBER_OF_COMPONENTS
      P=DEPENDENT_INTERPOLATED_POINT%VALUES(PRESSURE_COMPONENT,NO_PART_DERIV)
      B=[2.0_DP*C(2),2.0_DP*C(3),2.0_DP*C(3),C(4),C(4),C(3)] ![2*b_f,2*b_t,2*b_t,b_ft,b_ft,b_t]
      E=[0.5_DP*(AZL(1,1)-1.0_DP),0.5_DP*(AZL(2,2)-1.0_DP),0.5_DP*(AZL(3,3)-1.0_DP),AZL(2,1),AZL(3,1),AZL(3,2)] !(Modified) strain tensor in Voigt form.
      DQ_DE=B*E
      TEMPTERM1=0.5_DP*C(1)*EXP(0.5_DP*DOT_PRODUCT(E,DQ_DE))
      !Calculate 2nd Piola tensor (in Voigt form)
      STRESS_TENSOR=TEMPTERM1*DQ_DE + P*AZUv
      IF(EQUATIONS_SET%specification(3)==EQUATIONS_SET_GUCCIONE_ACTIVECONTRACTION_SUBTYPE) THEN
        !add active contraction stress values
        CALL Field_VariableGet(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VARIABLE,ERR,ERROR,*999)
        DO i=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
          dof_idx=FIELD_VARIABLE%COMPONENTS(i)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP% &
            & GAUSS_POINTS(GAUSS_POINT_NUMBER,ELEMENT_NUMBER)
          CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,dof_idx,VALUE,ERR,ERROR,*999)
          STRESS_TENSOR(i)=STRESS_TENSOR(i)+VALUE
        ENDDO
      ENDIF

      !\todo blas has routines specifically for symmetric matrices, so it would be worth to check if these could give some speedup.

      ! Calculate material elasticity tensor c (in Voigt form).
      ! First calculate lower part of 6X6 matrix
      DO j=1,6
        DO i=j,6
          ELASTICITY_TENSOR(i,j)=TEMPTERM1*DQ_DE(i)*DQ_DE(j)
        ENDDO
      ENDDO
      B=[2.0_DP*C(2),2.0_DP*C(3),2.0_DP*C(3),C(4),C(4),C(3)]
      DO i=1,6
        ELASTICITY_TENSOR(i,i)=ELASTICITY_TENSOR(i,i)+TEMPTERM1*B(i)
      ENDDO
      ! Then calculate upper part.
      DO j=2,6
        DO i=1,j-1
          ELASTICITY_TENSOR(i,j)=ELASTICITY_TENSOR(j,i)
        ENDDO
      ENDDO

      !Add volumetric part of elasticity tensor - p*d(C^-1)/dE.
      ELASTICITY_TENSOR=ELASTICITY_TENSOR + P*AZUE

      !Hydrostatic portion of the elasticity tensor (dS/dp)
      HYDRO_ELASTICITY_VOIGT = AZUv

      !Do push-forward of 2nd Piola tensor and the material elasticity tensor.
      CALL FINITE_ELASTICITY_PUSH_STRESS_TENSOR(STRESS_TENSOR,DZDNU,Jznu,ERR,ERROR,*999)
      CALL FINITE_ELASTICITY_PUSH_STRESS_TENSOR(HYDRO_ELASTICITY_VOIGT,DZDNU,Jznu,ERR,ERROR,*999)
      CALL FINITE_ELASTICITY_PUSH_ELASTICITY_TENSOR(ELASTICITY_TENSOR,DZDNU,Jznu,ERR,ERROR,*999)
    CASE DEFAULT
      LOCAL_ERROR="Analytic Jacobian has not been implemented for the third equations set specification of "// &
        & TRIM(NumberToVString(EQUATIONS_SET%specification(3),"*",ERR,ERROR))
      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
    END SELECT

    EXITS("FINITE_ELASTICITY_GAUSS_ELASTICITY_TENSOR")
    RETURN
999 ERRORSEXITS("FINITE_ELASTICITY_GAUSS_ELASTICITY_TENSOR",ERR,ERROR)
    RETURN 1

  END SUBROUTINE FINITE_ELASTICITY_GAUSS_ELASTICITY_TENSOR

  !
  !================================================================================================================================
  !

  !>Evaluates the element Jacobian matrix for the given element number for a finite elasticity class finite element equation set.
  SUBROUTINE FiniteElasticity_FiniteElementJacobianEvaluate(EQUATIONS_SET,ELEMENT_NUMBER,err,error,*)
    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to evaluate the Jacobian for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: FIELD_VAR_TYPE,ng,nh,ns,nhs,ni,mh,ms,mhs,oh
    INTEGER(INTG) :: PRESSURE_COMPONENT
    INTEGER(INTG) :: SUM_ELEMENT_PARAMETERS,TOTAL_NUMBER_OF_SURFACE_PRESSURE_CONDITIONS
    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS,NUMBER_OF_XI
    INTEGER(INTG) :: ELEMENT_BASE_DOF_INDEX(4),component_idx,component_idx2
    INTEGER(INTG), PARAMETER :: OFF_DIAG_COMP(3)=[0,1,3],OFF_DIAG_DEP_VAR1(3)=[1,1,2],OFF_DIAG_DEP_VAR2(3)=[2,3,3]
    INTEGER(INTG) :: MESH_COMPONENT_NUMBER,NUMBER_OF_ELEMENT_PARAMETERS(4)
    REAL(DP) :: DZDNU(3,3),CAUCHY_TENSOR(3,3),HYDRO_ELASTICITY_TENSOR(3,3)
    REAL(DP) :: JGW_SUB_MAT(3,3)
    REAL(DP) :: TEMPVEC(3)
    REAL(DP) :: STRESS_TENSOR(6),ELASTICITY_TENSOR(6,6),HYDRO_ELASTICITY_VOIGT(6)
    REAL(DP) :: DPHIDZ(3,64,3),DJDZ(64,3)
    REAL(DP) :: JGW_DPHINS_DZ,JGW_DPHIMS_DZ,PHIMS,PHINS,TEMPTERM
    REAL(DP) :: Jznu,JGW,SUM1,SUM2
    TYPE(QUADRATURE_SCHEME_PTR_TYPE) :: QUADRATURE_SCHEMES(4)
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: geometricInterpPoint,fibreInterpPoint, &
      & materialsInterpPoint,dependentInterpPoint
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: geometricInterpPointMetrics, &
      & dependentInterpPointMetrics
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsJacobianType), POINTER :: jacobianMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,GEOMETRIC_FIELD,MATERIALS_FIELD,FIBRE_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: DEPENDENT_QUADRATURE_SCHEME
    INTEGER(INTG) :: EQUATIONS_SET_SUBTYPE

    ENTERS("FiniteElasticity_FiniteElementJacobianEvaluate",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      EQUATIONS_SET_SUBTYPE = EQUATIONS_SET%SPECIFICATION(3)
      IF(ASSOCIATED(EQUATIONS)) THEN
        NULLIFY(vectorEquations)
        CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
        vectorMatrices=>vectorEquations%vectorMatrices
        nonlinearMatrices=>vectorMatrices%nonlinearMatrices
        jacobianMatrix=>nonlinearMatrices%jacobians(1)%ptr
        IF(jacobianMatrix%updateJacobian) THEN
          IF (EQUATIONS_SET_SUBTYPE == EQUATIONS_SET_REFERENCE_STATE_TRANSVERSE_GUCCIONE_SUBTYPE) THEN
            DEPENDENT_FIELD=>equations%interpolation%geometricField
            GEOMETRIC_FIELD=>equations%interpolation%dependentField
          ELSE
            DEPENDENT_FIELD=>equations%interpolation%dependentField
            GEOMETRIC_FIELD=>equations%interpolation%geometricField
          END IF
          MATERIALS_FIELD=>equations%interpolation%materialsField
          FIBRE_FIELD=>equations%interpolation%fibreField

          DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%ptr% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          DEPENDENT_QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
          
          NUMBER_OF_DIMENSIONS=EQUATIONS_SET%REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS
          NUMBER_OF_XI=DEPENDENT_BASIS%NUMBER_OF_XI

          vectorMapping=>vectorEquations%vectorMapping
          nonlinearMapping=>vectorMapping%nonlinearMapping
          
          FIELD_VARIABLE=>nonlinearMapping%residualVariables(1)%ptr
          FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE

          PRESSURE_COMPONENT=FIELD_VARIABLE%NUMBER_OF_COMPONENTS

          BOUNDARY_CONDITIONS=>EQUATIONS_SET%BOUNDARY_CONDITIONS
          CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,EQUATIONS_SET%equations%vectorEquations%vectorMapping% &
            & rhsMapping%rhsVariable,BOUNDARY_CONDITIONS_VARIABLE,err,error,*999)
          TOTAL_NUMBER_OF_SURFACE_PRESSURE_CONDITIONS=BOUNDARY_CONDITIONS_VARIABLE%DOF_COUNTS(BOUNDARY_CONDITION_PRESSURE)+ &
            & BOUNDARY_CONDITIONS_VARIABLE%DOF_COUNTS(BOUNDARY_CONDITION_PRESSURE_INCREMENTED)
        
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
            & dependentInterpParameters(FIELD_VAR_TYPE)%ptr,err,error,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
            & geometricInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
            & materialsInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
          IF(ASSOCIATED(FIBRE_FIELD)) THEN
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
              & fibreInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
          END IF

          !Point interpolation pointer
          geometricInterpPoint=>equations%interpolation%geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr
          geometricInterpPointMetrics=>equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr
          IF (EQUATIONS_SET_SUBTYPE == EQUATIONS_SET_REFERENCE_STATE_TRANSVERSE_GUCCIONE_SUBTYPE) THEN
            dependentInterpPoint=>equations%interpolation%geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr
            dependentInterpPointMetrics=>equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr
            geometricInterpPoint=>equations%interpolation%dependentInterpPoint(FIELD_VAR_TYPE)%ptr
            geometricInterpPointMetrics=>equations%interpolation%dependentInterpPointMetrics(FIELD_VAR_TYPE)%ptr
          ELSE
            geometricInterpPoint=>equations%interpolation%geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr
            geometricInterpPointMetrics=>equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr
            dependentInterpPoint=>equations%interpolation%dependentInterpPoint(FIELD_VAR_TYPE)%ptr
            dependentInterpPointMetrics=>equations%interpolation%dependentInterpPointMetrics(FIELD_VAR_TYPE)%ptr
          END IF
          IF(ASSOCIATED(FIBRE_FIELD)) THEN
            fibreInterpPoint=>equations%interpolation%fibreInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr
          END IF
          materialsInterpPoint=>equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr
          
          SUM_ELEMENT_PARAMETERS=0
          !Loop over geometric dependent basis functions.
          DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
            MESH_COMPONENT_NUMBER=FIELD_VARIABLE%COMPONENTS(nh)%MESH_COMPONENT_NUMBER
            DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%ptr% &
              & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
            QUADRATURE_SCHEMES(nh)%ptr=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
            IF(FIELD_VARIABLE%COMPONENTS(nh)%INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN
              NUMBER_OF_ELEMENT_PARAMETERS(nh)=DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
            ELSEIF(FIELD_VARIABLE%COMPONENTS(nh)%INTERPOLATION_TYPE==FIELD_ELEMENT_BASED_INTERPOLATION) THEN
              NUMBER_OF_ELEMENT_PARAMETERS(nh)=1
            ENDIF
            ELEMENT_BASE_DOF_INDEX(nh)=SUM_ELEMENT_PARAMETERS
            SUM_ELEMENT_PARAMETERS=SUM_ELEMENT_PARAMETERS+NUMBER_OF_ELEMENT_PARAMETERS(nh)
          ENDDO !nh

          !Loop over all Gauss points
          DO ng=1,DEPENDENT_QUADRATURE_SCHEME%NUMBER_OF_GAUSS
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng, &
              & dependentInterpPoint,err,error,*999)
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng, &
              & geometricInterpPoint,err,error,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(COORDINATE_JACOBIAN_VOLUME_TYPE, &
              & geometricInterpPointMetrics,err,error,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(COORDINATE_JACOBIAN_VOLUME_TYPE, &
              & dependentInterpPointMetrics,err,error,*999)
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng, &
              & materialsInterpPoint,err,error,*999)
            IF(ASSOCIATED(FIBRE_FIELD)) THEN
              CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng, &
                & fibreInterpPoint,err,error,*999)
            ENDIF

            Jznu=dependentInterpPointMetrics%JACOBIAN/geometricInterpPointMetrics%JACOBIAN 
            JGW=dependentInterpPointMetrics%JACOBIAN*DEPENDENT_QUADRATURE_SCHEME%GAUSS_WEIGHTS(ng)
            
            !Loop over geometric dependent basis functions.
            DO nh=1,NUMBER_OF_DIMENSIONS
              DO ns=1,NUMBER_OF_ELEMENT_PARAMETERS(nh)
                !Loop over derivative directions.
                SUM2=0.0_DP
                DO mh=1,NUMBER_OF_DIMENSIONS
                  SUM1=0.0_DP
                  DO ni=1,NUMBER_OF_XI
                    SUM1=SUM1+QUADRATURE_SCHEMES(nh)%PTR%GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)* &
                      & dependentInterpPointMetrics%DXI_DX(ni,mh)
                    SUM2=SUM2+QUADRATURE_SCHEMES(mh)%PTR%GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)* &
                      & dependentInterpPointMetrics%DXI_DX(ni,mh)*dependentInterpPointMetrics%GU(ni,mh)
                  ENDDO !mi
                  DPHIDZ(mh,ns,nh)=SUM1
                ENDDO !mh
                DJDZ(ns,nh)=SUM2*dependentInterpPointMetrics%jacobian
              ENDDO !ns
            ENDDO !nh

            CALL FiniteElasticity_GaussDeformationGradientTensor(dependentInterpPointMetrics, &
              & geometricInterpPointMetrics,fibreInterpPoint,dZdNu,err,error,*999)

            CALL FINITE_ELASTICITY_GAUSS_ELASTICITY_TENSOR(EQUATIONS_SET,dependentInterpPoint, &
              & materialsInterpPoint,ELASTICITY_TENSOR,HYDRO_ELASTICITY_VOIGT,STRESS_TENSOR, &
              & DZDNU,Jznu,ELEMENT_NUMBER,ng,ERR,ERROR,*999)

            !Convert from Voigt form to tensor form.
            DO nh=1,NUMBER_OF_DIMENSIONS
              DO mh=1,NUMBER_OF_DIMENSIONS
                CAUCHY_TENSOR(mh,nh)=STRESS_TENSOR(TENSOR_TO_VOIGT3(mh,nh))
                HYDRO_ELASTICITY_TENSOR(mh,nh)=HYDRO_ELASTICITY_VOIGT(TENSOR_TO_VOIGT3(mh,nh))
              ENDDO
            ENDDO

            !1) loop over mh=nh
            !Loop over element columns belonging to geometric dependent variables
            nhs=0
            DO nh=1,NUMBER_OF_DIMENSIONS
              JGW_SUB_MAT=JGW*(ELASTICITY_TENSOR(TENSOR_TO_VOIGT(1:NUMBER_OF_DIMENSIONS,nh,NUMBER_OF_DIMENSIONS), &
                & TENSOR_TO_VOIGT(1:NUMBER_OF_DIMENSIONS,nh,NUMBER_OF_DIMENSIONS))+ &
                & CAUCHY_TENSOR(1:NUMBER_OF_DIMENSIONS,1:NUMBER_OF_DIMENSIONS))              
              DO ns=1,NUMBER_OF_ELEMENT_PARAMETERS(nh)
                TEMPVEC=MATMUL(JGW_SUB_MAT,DPHIDZ(1:NUMBER_OF_DIMENSIONS,ns,nh))
                nhs=nhs+1
                mhs=nhs-1
                !Loop over element rows belonging to geometric dependent variables
                DO ms=ns,NUMBER_OF_ELEMENT_PARAMETERS(nh)
                  mhs=mhs+1
                  jacobianMatrix%elementJacobian%matrix(mhs,nhs)=jacobianMatrix%elementJacobian%matrix(mhs,nhs)+ &
                    & DOT_PRODUCT(dPhiDZ(:,ms,nh),TEMPVEC)
                  DO component_idx=1,NUMBER_OF_DIMENSIONS
                    DO component_idx2=1,NUMBER_OF_DIMENSIONS
                      TEMPTERM=CAUCHY_TENSOR(component_idx,component_idx2)* &
                      & dPhidZ(component_idx2,ms,component_idx)
                    ENDDO
                  ENDDO
                ENDDO !ms
              ENDDO !ns
            ENDDO !nh

            !2) loop over mh>nh
            !Loop over element columns belonging to geometric dependent variables
            DO oh=1,OFF_DIAG_COMP(NUMBER_OF_DIMENSIONS)
              nh=OFF_DIAG_DEP_VAR1(oh)
              mh=OFF_DIAG_DEP_VAR2(oh)
              nhs=ELEMENT_BASE_DOF_INDEX(nh)
              JGW_SUB_MAT=JGW*(ELASTICITY_TENSOR(TENSOR_TO_VOIGT3(1:NUMBER_OF_DIMENSIONS,mh), &
                & TENSOR_TO_VOIGT3(1:NUMBER_OF_DIMENSIONS,nh)))
              
              DO ns=1,NUMBER_OF_ELEMENT_PARAMETERS(nh)
                !Loop over element rows belonging to geometric dependent variables
                TEMPVEC=MATMUL(JGW_SUB_MAT,DPHIDZ(1:NUMBER_OF_DIMENSIONS,ns,nh))
                nhs=nhs+1
                mhs=ELEMENT_BASE_DOF_INDEX(mh)
                DO ms=1,NUMBER_OF_ELEMENT_PARAMETERS(mh)
                  mhs=mhs+1
                  jacobianMatrix%elementJacobian%matrix(mhs,nhs)=jacobianMatrix%elementJacobian%matrix(mhs,nhs)+ &
                    & DOT_PRODUCT(dPhidZ(:,ms,mh),TEMPVEC)
                  DO component_idx=1,NUMBER_OF_DIMENSIONS
                    DO component_idx2=1,NUMBER_OF_DIMENSIONS
                      TEMPTERM=CAUCHY_TENSOR(component_idx,component_idx2)* &
                      & dPhidZ(component_idx2,ms,component_idx)
                    ENDDO
                  ENDDO
                  !JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(mhs,nhs)=JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(mhs,nhs)+ &
                  !  & TEMPTERM*DJDZ(ms,nh)*DEPENDENT_QUADRATURE_SCHEME%GAUSS_WEIGHTS(ng)
                ENDDO !ms
              ENDDO !ns
            ENDDO

            !3) loop over all nh and pressure component
            nhs=0
            IF(FIELD_VARIABLE%COMPONENTS(PRESSURE_COMPONENT)%INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN !node based
              !Loop over element rows belonging to geometric dependent variables
              DO nh=1,NUMBER_OF_DIMENSIONS
                DO ns=1,NUMBER_OF_ELEMENT_PARAMETERS(nh)
                  JGW_DPHINS_DZ=JGW*DPHIDZ(nh,ns,nh)
                  nhs=nhs+1
                 !Loop over element rows belonging to hydrostatic pressure
                  mhs=ELEMENT_BASE_DOF_INDEX(PRESSURE_COMPONENT)
                  DO ms=1,NUMBER_OF_ELEMENT_PARAMETERS(PRESSURE_COMPONENT)
                    mhs=mhs+1
                    PHIMS=QUADRATURE_SCHEMES(PRESSURE_COMPONENT)%ptr%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                    jacobianMatrix%elementJacobian%matrix(mhs,nhs)=jacobianMatrix%elementJacobian%matrix(mhs,nhs)+ &
                      & JGW_DPHINS_DZ*PHIMS
                  ENDDO !ms
                ENDDO !ns
              ENDDO !nh
            ELSEIF(FIELD_VARIABLE%COMPONENTS(PRESSURE_COMPONENT)%INTERPOLATION_TYPE==FIELD_ELEMENT_BASED_INTERPOLATION) THEN !element based
              !Loop over element rows belonging to geometric dependent variables
              DO nh=1,NUMBER_OF_DIMENSIONS
                DO ns=1,NUMBER_OF_ELEMENT_PARAMETERS(nh)
                  JGW_DPHINS_DZ=JGW*DPHIDZ(nh,ns,nh)
                  nhs=nhs+1
                  !Loop over element rows belonging to hydrostatic pressure
                  mhs=ELEMENT_BASE_DOF_INDEX(PRESSURE_COMPONENT)+1
                  jacobianMatrix%elementJacobian%matrix(mhs,nhs)=jacobianMatrix%elementJacobian%matrix(mhs,nhs)+ &
                    & JGW_DPHINS_DZ
                ENDDO !ns
              ENDDO !nh
            ENDIF

            !4) Loop over all mh pressure component
            mhs=0
            IF(FIELD_VARIABLE%COMPONENTS(PRESSURE_COMPONENT)%INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN !node based
              !Loop over element columns belonging to geometric dependent variables.
              DO mh=1,NUMBER_OF_DIMENSIONS
                DO ms=1,NUMBER_OF_ELEMENT_PARAMETERS(mh)
                  TEMPVEC=MATMUL(HYDRO_ELASTICITY_TENSOR,DPHIDZ(:,ms,mh))
                  JGW_DPHIMS_DZ=JGW*TEMPVEC(mh)
                  mhs=mhs+1
                  !Loop over element columns belonging to hydrostatic pressure
                  nhs=ELEMENT_BASE_DOF_INDEX(PRESSURE_COMPONENT)
                  DO ns=1,NUMBER_OF_ELEMENT_PARAMETERS(PRESSURE_COMPONENT)
                    nhs=nhs+1
                    PHINS=QUADRATURE_SCHEMES(PRESSURE_COMPONENT)%PTR%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)
                    jacobianMatrix%elementJacobian%matrix(mhs,nhs)=jacobianMatrix%elementJacobian%matrix(mhs,nhs)+ &
                      & JGW_DPHIMS_DZ*PHINS
                  ENDDO !ns
                ENDDO !ms
              ENDDO !mh
            ELSEIF(FIELD_VARIABLE%COMPONENTS(PRESSURE_COMPONENT)%INTERPOLATION_TYPE==FIELD_ELEMENT_BASED_INTERPOLATION) THEN !element based
              !Loop over element columns belonging to geometric dependent variables.
              DO mh=1,NUMBER_OF_DIMENSIONS
                DO ms=1,NUMBER_OF_ELEMENT_PARAMETERS(mh)
                  TEMPVEC=MATMUL(HYDRO_ELASTICITY_TENSOR,DPHIDZ(:,ms,mh))
                  JGW_DPHIMS_DZ=JGW*TEMPVEC(mh)
                  mhs=mhs+1
                  !Loop over element columns belonging to hydrostatic pressure.
                  nhs=ELEMENT_BASE_DOF_INDEX(PRESSURE_COMPONENT)+1
                  jacobianMatrix%elementJacobian%matrix(mhs,nhs)=jacobianMatrix%elementJacobian%matrix(mhs,nhs) + &
                    & JGW_DPHIMS_DZ
                ENDDO !ms
              ENDDO !mh
            ENDIF
            ! No loop over element columns and rows belonging both to hydrostatic pressure because it is zero.
          ENDDO !ng

          !Scale factor adjustment
          IF(DEPENDENT_FIELD%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
            !Following call is necessary, otherwise wrong face scale factors from function call to surface pressure jacobian are
            !used.
            CALL Field_InterpolationParametersScaleFactorsElementGet(ELEMENT_NUMBER, &
              & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr,err,error,*999) 
            nhs=0          
            ! Loop over element columns
            DO nh=1,NUMBER_OF_DIMENSIONS
              DO ns=1,NUMBER_OF_ELEMENT_PARAMETERS(nh)
                nhs=nhs+1
                mhs=nhs-1
                ! Loop over element rows
                DO ms=ns,NUMBER_OF_ELEMENT_PARAMETERS(nh)
                  mhs=mhs+1
                  jacobianMatrix%elementJacobian%matrix(mhs,nhs)=jacobianMatrix%elementJacobian%matrix(mhs,nhs)* &
                    & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%SCALE_FACTORS(ms,nh)* &
                    & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%SCALE_FACTORS(ns,nh)
                ENDDO !ms
              ENDDO !ns
            ENDDO !nh
            DO oh=1,OFF_DIAG_COMP(NUMBER_OF_DIMENSIONS)
              nh=OFF_DIAG_DEP_VAR1(oh)
              mh=OFF_DIAG_DEP_VAR2(oh)
              nhs=ELEMENT_BASE_DOF_INDEX(nh)
              DO ns=1,NUMBER_OF_ELEMENT_PARAMETERS(nh)
                nhs=nhs+1
                mhs=ELEMENT_BASE_DOF_INDEX(mh)
                !Loop over element rows belonging to geometric dependent variables
                DO ms=1,NUMBER_OF_ELEMENT_PARAMETERS(mh)
                  mhs=mhs+1
                  jacobianMatrix%elementJacobian%matrix(mhs,nhs)=jacobianMatrix%elementJacobian%matrix(mhs,nhs)* &
                    & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%SCALE_FACTORS(ms,mh)* &
                    & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%SCALE_FACTORS(ns,nh)
                ENDDO !ms    
              ENDDO !ns
            ENDDO

            nhs=0
            IF(FIELD_VARIABLE%COMPONENTS(PRESSURE_COMPONENT)%INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN !node based
              !Loop over element rows belonging to geometric dependent variables
              DO nh=1,NUMBER_OF_DIMENSIONS
                DO ns=1,NUMBER_OF_ELEMENT_PARAMETERS(nh)
                  nhs=nhs+1
                  !Loop over element rows belonging to hydrostatic pressure
                  mhs=ELEMENT_BASE_DOF_INDEX(PRESSURE_COMPONENT)
                  DO ms=1,NUMBER_OF_ELEMENT_PARAMETERS(PRESSURE_COMPONENT)
                    mhs=mhs+1
                    jacobianMatrix%elementJacobian%matrix(nhs,mhs)=jacobianMatrix%elementJacobian%matrix(nhs,mhs)* &
                      & EQUATIONS%INTERPOLATION%dependentInterpParameters(FIELD_VAR_TYPE)%PTR% &
                      & SCALE_FACTORS(ms,PRESSURE_COMPONENT)* &
                      & EQUATIONS%INTERPOLATION%dependentInterpParameters(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ns,nh)
                    jacobianMatrix%elementJacobian%matrix(mhs,nhs)=jacobianMatrix%elementJacobian%matrix(mhs,nhs)* &
                      & EQUATIONS%INTERPOLATION%dependentInterpParameters(FIELD_VAR_TYPE)%PTR% &
                      & SCALE_FACTORS(ms,PRESSURE_COMPONENT)* &
                      & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%SCALE_FACTORS(ns,nh)
                  ENDDO !ms    
                ENDDO !ns
              ENDDO !nh
            ELSEIF(FIELD_VARIABLE%COMPONENTS(PRESSURE_COMPONENT)%INTERPOLATION_TYPE==FIELD_ELEMENT_BASED_INTERPOLATION) THEN !element based
              !Loop over element rows belonging to geometric dependent variables
              DO nh=1,NUMBER_OF_DIMENSIONS
                DO ns=1,NUMBER_OF_ELEMENT_PARAMETERS(nh)
                  nhs=nhs+1
                  !Loop over element rows belonging to hydrostatic pressure
                  mhs=ELEMENT_BASE_DOF_INDEX(PRESSURE_COMPONENT)+1
                  jacobianMatrix%elementJacobian%matrix(mhs,nhs)=jacobianMatrix%elementJacobian%matrix(mhs,nhs)* &
                    & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%SCALE_FACTORS(ns,nh)
                ENDDO !ns
              ENDDO !nh
            ENDIF
          ENDIF

          !Mirror the Jacobian matrix except for the hydrostatic rows and columns, which are not necessarily symmetric.
          DO nhs=2,ELEMENT_BASE_DOF_INDEX(PRESSURE_COMPONENT)
            DO mhs=1,nhs-1
              jacobianMatrix%elementJacobian%matrix(mhs,nhs)=jacobianMatrix%elementJacobian%matrix(nhs,mhs)
            ENDDO !mhs
          ENDDO !nhs

          !If unsymmetric pressure Jacobian uncomment this.
          !Call surface pressure term here: should only be executed if THIS element has surface pressure on it (direct or incremented)
          IF(DEPENDENT_FIELD%DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BOUNDARY_ELEMENT.AND. &
            & TOTAL_NUMBER_OF_SURFACE_PRESSURE_CONDITIONS>0) THEN    !
            CALL FiniteElasticity_SurfacePressureJacobianEvaluate(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FlagError("Equations set equations is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

    EXITS("FiniteElasticity_FiniteElementJacobianEvaluate")
    RETURN
999 ERRORS("FiniteElasticity_FiniteElementJacobianEvaluate",err,error)
    EXITS("FiniteElasticity_FiniteElementJacobianEvaluate")
    RETURN 1

  END SUBROUTINE FiniteElasticity_FiniteElementJacobianEvaluate

  !
  !================================================================================================================================
  !

  !>Push-forward the rank 4 elasticity tensor.
  SUBROUTINE FINITE_ELASTICITY_PUSH_ELASTICITY_TENSOR(ELASTICITY_TENSOR,DZDNU,Jznu,err,error,*)
    
    !Argument variables
    REAL(DP), INTENT(INOUT) :: ELASTICITY_TENSOR(6,6)
    REAL(DP), INTENT(IN) :: DZDNU(3,3)
    REAL(DP), INTENT(IN) :: Jznu
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: i,j
    REAL(DP) :: t(6,6),ttrans(6,6) 

    ENTERS("FINITE_ELASTICITY_PUSH_ELASTICITY_TENSOR",err,error,*999)

    DO j=1,3
      DO i=1,6
        t(i,j)=DZDNU(VOIGT_TO_TENSOR3(1,i),VOIGT_TO_TENSOR3(1,j))*DZDNU(VOIGT_TO_TENSOR3(2,i),VOIGT_TO_TENSOR3(2,j))
      ENDDO
    END DO
    DO j=4,6
      DO i=1,6
        t(i,j)=DZDNU(VOIGT_TO_TENSOR3(1,i),VOIGT_TO_TENSOR3(1,j))*DZDNU(VOIGT_TO_TENSOR3(2,i),VOIGT_TO_TENSOR3(2,j))+ &
          & DZDNU(VOIGT_TO_TENSOR3(1,i),VOIGT_TO_TENSOR3(2,j))*DZDNU(VOIGT_TO_TENSOR3(2,i),VOIGT_TO_TENSOR3(1,j))
      ENDDO
    END DO

    CALL MatrixTranspose(t,ttrans,err,error,*999)
    ELASTICITY_TENSOR=MATMUL(MATMUL(t,ELASTICITY_TENSOR),ttrans)/Jznu

    EXITS("FINITE_ELASTICITY_PUSH_ELASTICITY_TENSOR")
    RETURN
999 ERRORSEXITS("FINITE_ELASTICITY_PUSH_ELASTICITY_TENSOR",err,error)
    RETURN 1
  END SUBROUTINE FINITE_ELASTICITY_PUSH_ELASTICITY_TENSOR

  !
  !================================================================================================================================
  !

  !>Push-forward the rank 2 Piola stress tensor.
  SUBROUTINE FINITE_ELASTICITY_PUSH_STRESS_TENSOR(STRESS_TENSOR,DZDNU,Jznu,err,error,*)
    
    !Argument variables
    REAL(DP), INTENT(INOUT) :: STRESS_TENSOR(6)
    REAL(DP), INTENT(IN) :: DZDNU(3,3)
    REAL(DP), INTENT(IN) :: Jznu
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: i,j
    REAL(DP) :: t(6,6)

    ENTERS("FINITE_ELASTICITY_PUSH_STRESS_TENSOR",err,error,*999)

    DO j=1,3
      DO i=1,6
        t(i,j)=DZDNU(VOIGT_TO_TENSOR3(1,i),VOIGT_TO_TENSOR3(1,j))*DZDNU(VOIGT_TO_TENSOR3(2,i),VOIGT_TO_TENSOR3(2,j))
      ENDDO
    END DO
    DO j=4,6
      DO i=1,6
        t(i,j)=DZDNU(VOIGT_TO_TENSOR3(1,i),VOIGT_TO_TENSOR3(1,j))*DZDNU(VOIGT_TO_TENSOR3(2,i),VOIGT_TO_TENSOR3(2,j))+ &
          & DZDNU(VOIGT_TO_TENSOR3(1,i),VOIGT_TO_TENSOR3(2,j))*DZDNU(VOIGT_TO_TENSOR3(2,i),VOIGT_TO_TENSOR3(1,j))
      ENDDO
    END DO

    STRESS_TENSOR=MATMUL(t,STRESS_TENSOR)/Jznu

    EXITS("FINITE_ELASTICITY_PUSH_STRESS_TENSOR")
    RETURN
999 ERRORSEXITS("FINITE_ELASTICITY_PUSH_STRESS_TENSOR",err,error)
    RETURN 1
  END SUBROUTINE FINITE_ELASTICITY_PUSH_STRESS_TENSOR

  !
  !================================================================================================================================
  !

  !>Evaluates the residual and RHS vectors for a finite elasticity finite element equations set.
   SUBROUTINE FiniteElasticity_FiniteElementResidualEvaluate(EQUATIONS_SET,elementNumber,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS,COMPONENT_BASIS,dependentBasis
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,FIBRE_FIELD,GEOMETRIC_FIELD,MATERIALS_FIELD,EQUATIONS_SET_FIELD,SOURCE_FIELD
    TYPE(FIELD_TYPE), POINTER :: INDEPENDENT_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: DEPENDENT_QUADRATURE_SCHEME,COMPONENT_QUADRATURE_SCHEME,quadratureScheme
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: GEOMETRIC_INTERPOLATION_PARAMETERS, &
      & FIBRE_INTERPOLATION_PARAMETERS,MATERIALS_INTERPOLATION_PARAMETERS,DEPENDENT_INTERPOLATION_PARAMETERS, &
      & DARCY_DEPENDENT_INTERPOLATION_PARAMETERS,SOURCE_INTERPOLATION_PARAMETERS,DARCY_MATERIALS_INTERPOLATION_PARAMETERS, &
      & DENSITY_INTERPOLATION_PARAMETERS,INDEPENDENT_INTERPOLATION_PARAMETERS,prevDependentInterpParameters
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: GEOMETRIC_INTERPOLATED_POINT,FIBRE_INTERPOLATED_POINT, &
      & MATERIALS_INTERPOLATED_POINT,DEPENDENT_INTERPOLATED_POINT,DARCY_DEPENDENT_INTERPOLATED_POINT,SOURCE_INTERPOLATED_POINT, &
      & DENSITY_INTERPOLATED_POINT,INDEPENDENT_INTERPOLATED_POINT,DARCY_MATERIALS_INTERPOLATED_POINT, &
      & prevDependentInterpPoint
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: GEOMETRIC_INTERPOLATED_POINT_METRICS, &
      & DEPENDENT_INTERPOLATED_POINT_METRICS,prevDependentInterpPointMetrics
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS_1,GEOMETRIC_BASIS
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_ELEMENT_MAPPING
    TYPE(VARYING_STRING) :: localError,localWarning
    LOGICAL :: DARCY_DENSITY,DARCY_DEPENDENT
    INTEGER(INTG) :: component_idx,component_idx2,parameter_idx,gauss_idx,element_dof_idx,FIELD_VAR_TYPE,DARCY_FIELD_VAR_TYPE
    INTEGER(INTG) :: imatrix,Ncompartments,gaussIdx,rowIdx,columnIdx,componentIdx,elementDofIdx
    INTEGER(INTG) :: numberOfXDimensions,numberOfXiDimensions
    INTEGER(INTG) :: NDOFS,mh,ms,mhs,mi,nh,ns,rowComponentIdx,rowElementParameterIdx,rowElementDofIdx,xiIdx,columnComponentIdx, &
      & columnElementParameterIdx
    INTEGER(INTG) :: DEPENDENT_NUMBER_OF_COMPONENTS
    INTEGER(INTG) :: numberOfDimensions,NUMBER_OF_XI,HYDROSTATIC_PRESSURE_COMPONENT,hydrostaticPressureComponent,numberOfXi
    INTEGER(INTG) :: NUMBER_OF_FIELD_COMPONENT_INTERPOLATION_PARAMETERS
    INTEGER(INTG) :: DEPENDENT_COMPONENT_INTERPOLATION_TYPE,dependentComponentInterpolationType
    INTEGER(INTG) :: DEPENDENT_NUMBER_OF_GAUSS_POINTS       
    INTEGER(INTG) :: MESH_COMPONENT_1,MESH_COMPONENT_NUMBER,meshComponentNumber
    INTEGER(INTG) :: TOTAL_NUMBER_OF_SURFACE_PRESSURE_CONDITIONS
    INTEGER(INTG) :: var1 ! Variable number corresponding to 'U' in single physics case
    INTEGER(INTG) :: var2 ! Variable number corresponding to 'DELUDLEN' in single physics case
    INTEGER(INTG), POINTER :: EQUATIONS_SET_FIELD_DATA(:)
    REAL(DP) :: DZDNU(3,3),DZDNUT(3,3),dzdx(3,3),AZL(3,3),AZU(3,3),Fe(3,3),FeT(3,3),Fg(3,3),C(3,3),f(3,3),E(3,3),I3,P, &
      & piolaTensor(3,3),TEMP(3,3),prevdzdx(3,3),prevdZdNu(3,3),invPrevdZdNu(3,3)
    REAL(DP) :: cauchyTensor(3,3),JGW_CAUCHY_TENSOR(3,3),kirchoffTensor(3,3),STRESS_TENSOR(6),growthValues(3)
    REAL(DP) :: deformationGradientTensor(3,3),growthTensor(3,3),growthTensorInverse(3,3),growthTensorInverseTranspose(3,3), &
      & fibreGrowth,sheetGrowth,normalGrowth,fibreVector(3),sheetVector(3),normalVector(3)
    REAL(DP) :: dNudXi(3,3),dXidNu(3,3)
    REAL(DP) :: DFDZ(64,3,3) !temporary until a proper alternative is found
    REAL(DP) :: dPhidZ(3,64,3) !temporary until a proper alternative is found
    REAL(DP) :: GAUSS_WEIGHT,J,Je,Jg,Jznu,Jxxi,Jzxi,JGw,gaussWeight,prevJZxi
    REAL(DP) :: sum1,tempTerm1
    REAL(DP) :: THICKNESS ! for elastic membrane
    REAL(DP) :: DARCY_MASS_INCREASE,DARCY_VOL_INCREASE,DARCY_RHO_0_F,DENSITY  !coupling with Darcy model
    REAL(DP) :: Mfact, bfact, p0fact
    REAL(DP) :: dt,K,mu,a0,a1,b0,b1,m,Jr,kappa1,kappan,kappas,prevJF
    REAL(DP) :: alpha1,BePrime(3,3),BePrime1(3,3),BePrimeStar(3,3),Br(3,3),c0,c1,c2,Dbar(3,3),deltaEps,detdevBePrimePrime, &
      & devBePrimePrime(3,3),BePrimePrimeStar(3,3),devDbar(3,3),gePrimePrime(3,3),gePrimePrimeStar(3,3),devT(3,3),dtGamma, &
      & dtGamma0,dtGamma1,factor1,factor2,Fr(3,3),FrPrime(3,3),gammaEStar,gamma,hydBePrimeStar(3,3),hydDbar(3,3),ITens(3,3), &
      & kappa,tempTensor(3,3)
    INTEGER(INTG) :: numIterations,maxNumIterations
    REAL(DP) :: mu0,q,kt,kc,k1,k2,k3,k12,ka,Je1,Je2,Jh,GammaM,deltaAlphaNorm,deltaAlphaNormTolerance,detBePrime, &
      & detBePrimePrime,BeDDotH,BeDDotBe,detBe,ktJe,kcJe,macEe11,macEe22,macEe33,QQ,detInvFrT,hstep,hStepAlpha,detJacobian, &
      & H11,H22,H33,H12,H13,H23
    REAL(DP) :: alpha(2),perturbedAlpha(2),deltaAlpha(2),resid(2),perturbedResid(2),H(3,3),HPrimePrime(3,3), &
      & S1(3,3),S2(3,3),S(3,3,3,3),Be(3,3),uniBePrime1(3,3),invBePrime(3,3),Jacobian(2,2),lame(4),Ee(4,4),invBe(3,3), &
      & invBeS33(3,3),S33invBe(3,3),T0(3,3),T1(3,3),T2(3,3),T3(3,3),T4(3,3),T5(3,3),invFrT(3,3),BePrimePrime(3,3), &
      & invJacobian(2,2)
    REAL(DP) :: devH(3,3),HH(3,3),malpha1,B,Bepr(3,3),lame1,lame2,lame3,lamea,Ee11,Ee22,Ee33,Ee12,Eea,TT(3,3)
    REAL(DP) :: statev(13),ddsdde(3,3,3,3),props(8)
    REAL(DP) :: k4,b2,b3,b4,b5,b6,rho0,bigQ,delWdelQ,delQdelJe,delQdelAlpha,delWdelJe,delWdelAlpha
    REAL(DP) :: A1g,B0g,B1g,B2g,alphag,chi1g,chi2g,phig,R1g,R2g,Lg,x,y,z,r,theta,Dg22,time,FgInv(3,3)
    INTEGER(INTG) :: EQUATIONS_SET_SUBTYPE

    ENTERS("FiniteElasticity_FiniteElementResidualEvaluate",err,error,*999)

    NULLIFY(BOUNDARY_CONDITIONS,BOUNDARY_CONDITIONS_VARIABLE)
    NULLIFY(DEPENDENT_BASIS,COMPONENT_BASIS)
    NULLIFY(EQUATIONS,vectorMapping,vectorMatrices,nonlinearMatrices,rhsVector)
    NULLIFY(DEPENDENT_FIELD,FIBRE_FIELD,GEOMETRIC_FIELD,MATERIALS_FIELD,SOURCE_FIELD,INDEPENDENT_FIELD)
    NULLIFY(FIELD_VARIABLE)
    NULLIFY(DEPENDENT_QUADRATURE_SCHEME,COMPONENT_QUADRATURE_SCHEME)
    NULLIFY(GEOMETRIC_INTERPOLATION_PARAMETERS,FIBRE_INTERPOLATION_PARAMETERS,SOURCE_INTERPOLATION_PARAMETERS)
    NULLIFY(MATERIALS_INTERPOLATION_PARAMETERS,DEPENDENT_INTERPOLATION_PARAMETERS,prevDependentInterpParameters)
    NULLIFY(INDEPENDENT_INTERPOLATION_PARAMETERS,DARCY_MATERIALS_INTERPOLATION_PARAMETERS)
    NULLIFY(DARCY_DEPENDENT_INTERPOLATION_PARAMETERS,DENSITY_INTERPOLATION_PARAMETERS)
    NULLIFY(GEOMETRIC_INTERPOLATED_POINT,FIBRE_INTERPOLATED_POINT,SOURCE_INTERPOLATED_POINT)
    NULLIFY(GEOMETRIC_INTERPOLATED_POINT_METRICS,DEPENDENT_INTERPOLATED_POINT_METRICS,prevDependentInterpPointMetrics)
    NULLIFY(MATERIALS_INTERPOLATED_POINT,DEPENDENT_INTERPOLATED_POINT,DARCY_DEPENDENT_INTERPOLATED_POINT)
    NULLIFY(DENSITY_INTERPOLATED_POINT,INDEPENDENT_INTERPOLATED_POINT,prevDependentInterpPoint)
    NULLIFY(DEPENDENT_BASIS_1)
    NULLIFY(DECOMPOSITION)
    NULLIFY(EQUATIONS_SET_FIELD_DATA)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a finite elasticity type equations set.", &
          & err,error,*999)
      END IF
      EQUATIONS_SET_SUBTYPE = EQUATIONS_SET%SPECIFICATION(3)
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        NULLIFY(vectorEquations)
        CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
        !Which variables are we working with - find the variable pair used for this equations set
        !\todo: put in checks for all the objects/mappings below (do we want to do this for every element?)
        var1=vectorEquations%vectorMapping%nonlinearMapping%residualVariables(1)%ptr%VARIABLE_NUMBER ! number for 'U'
        var2=vectorEquations%vectorMapping%rhsMapping%rhsVariable%VARIABLE_NUMBER ! number for 'DELUDELN'

        !Grab pointers: matrices, fields, decomposition, basis
        !\todo: see if we can separate this residual evaluation from the pressure boundary conditions somehow
        !so that the equations set doesn't need to maintain a pointer to the boundary conditions
        BOUNDARY_CONDITIONS=>EQUATIONS_SET%BOUNDARY_CONDITIONS
        CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,EQUATIONS_SET%equations%vectorEquations%vectorMapping% &
          & rhsMapping%rhsVariable,BOUNDARY_CONDITIONS_VARIABLE,err,error,*999)
        TOTAL_NUMBER_OF_SURFACE_PRESSURE_CONDITIONS=BOUNDARY_CONDITIONS_VARIABLE%DOF_COUNTS(BOUNDARY_CONDITION_PRESSURE)+ &
          & BOUNDARY_CONDITIONS_VARIABLE%DOF_COUNTS(BOUNDARY_CONDITION_PRESSURE_INCREMENTED)

        vectorMatrices=>vectorEquations%vectorMatrices
        nonlinearMatrices=>vectorMatrices%nonlinearMatrices
        rhsVector=>vectorMatrices%rhsVector
        vectorMapping =>vectorEquations%vectorMapping

        IF (EQUATIONS_SET_SUBTYPE == EQUATIONS_SET_REFERENCE_STATE_TRANSVERSE_GUCCIONE_SUBTYPE) THEN
          DEPENDENT_FIELD  =>equations%interpolation%geometricField
          GEOMETRIC_FIELD  =>equations%interpolation%dependentField
        ELSE
          GEOMETRIC_FIELD  =>equations%interpolation%geometricField
          DEPENDENT_FIELD  =>equations%interpolation%dependentField
        ENDIF
        FIBRE_FIELD      =>equations%interpolation%fibreField
        MATERIALS_FIELD  =>equations%interpolation%materialsField
        SOURCE_FIELD     =>equations%interpolation%sourceField
        INDEPENDENT_FIELD=>equations%interpolation%independentField

        DECOMPOSITION    =>DEPENDENT_FIELD%DECOMPOSITION
        MESH_COMPONENT_NUMBER = DECOMPOSITION%MESH_COMPONENT_NUMBER

        DOMAIN_ELEMENT_MAPPING=>DECOMPOSITION%DOMAIN(1)%ptr%MAPPINGS%ELEMENTS

        DEPENDENT_BASIS=>DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%ptr%TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS       
        DEPENDENT_QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
        DEPENDENT_NUMBER_OF_GAUSS_POINTS=DEPENDENT_QUADRATURE_SCHEME%NUMBER_OF_GAUSS
        DEPENDENT_NUMBER_OF_COMPONENTS=DEPENDENT_FIELD%VARIABLES(var1)%NUMBER_OF_COMPONENTS
        GEOMETRIC_BASIS=>GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%ptr% &
          & TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS

        numberOfDimensions=EQUATIONS_SET%REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS
        numberOfXi=DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%ptr%TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS%NUMBER_OF_XI

        !Initialise tensors and matrices
        CALL IdentityMatrix(DZDNU,err,error,*999)
        CALL IdentityMatrix(piolaTensor,err,error,*999)
        CALL IdentityMatrix(cauchyTensor,err,error,*999)
        DFDZ=0.0_DP ! (parameter_idx,component_idx)

        !Set flags for coupled finite elasticity and Darcy problems
        !Check if we need Darcy materials field for Density
        IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_STATIC_INRIA_SUBTYPE .OR. &
          & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_HOLMES_MOW_SUBTYPE .OR. &
          & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_ELASTICITY_FLUID_PRES_HOLMES_MOW_ACTIVE_SUBTYPE) THEN
          DARCY_DENSITY=.TRUE.
        ELSE
          DARCY_DENSITY=.FALSE.
        ENDIF
        !Check if we need Darcy dependent field
        IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE .OR. &
          & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE .OR. &
          & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE .OR. &
          & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_STATIC_INRIA_SUBTYPE .OR. &
          & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_HOLMES_MOW_SUBTYPE .OR. &
          & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_ELASTICITY_FLUID_PRES_HOLMES_MOW_ACTIVE_SUBTYPE) THEN
          DARCY_DEPENDENT=.TRUE.
        ELSE
          DARCY_DEPENDENT=.FALSE.
        ENDIF

        !Grab interpolation parameters
        FIELD_VARIABLE=>EQUATIONS_SET%equations%vectorEquations%vectorMapping%nonlinearMapping%residualVariables(1)%ptr
        FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
        IF (EQUATIONS_SET_SUBTYPE == EQUATIONS_SET_REFERENCE_STATE_TRANSVERSE_GUCCIONE_SUBTYPE) THEN
          GEOMETRIC_INTERPOLATION_PARAMETERS=>equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr
          DEPENDENT_INTERPOLATION_PARAMETERS=>equations%interpolation%geometricInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr
        ELSE
          DEPENDENT_INTERPOLATION_PARAMETERS=>equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr
          GEOMETRIC_INTERPOLATION_PARAMETERS=>equations%interpolation%geometricInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr
        ENDIF
        IF(EQUATIONS%timeDependence/=EQUATIONS_STATIC) THEN
          prevDependentInterpParameters=>equations%interpolation%prevDependentInterpParameters(FIELD_VAR_TYPE)%ptr
        ENDIF
        IF(ASSOCIATED(FIBRE_FIELD)) THEN
          FIBRE_INTERPOLATION_PARAMETERS=>equations%interpolation%fibreInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr
        ENDIF
        IF(ASSOCIATED(MATERIALS_FIELD)) THEN
          MATERIALS_INTERPOLATION_PARAMETERS=>equations%interpolation%materialsInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr
!          DENSITY_INTERPOLATION_PARAMETERS=>equations%interpolation%materialsInterpParameters(FIELD_V_VARIABLE_TYPE)%ptr
        ENDIF
        IF(DARCY_DEPENDENT) THEN
          DARCY_DEPENDENT_INTERPOLATION_PARAMETERS=>equations%interpolation%dependentInterpParameters(FIELD_V_VARIABLE_TYPE)%ptr
        ELSE IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_STANDARD_MONODOMAIN_ELASTICITY_SUBTYPE) THEN
          INDEPENDENT_INTERPOLATION_PARAMETERS=>equations%interpolation%independentInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr
        ENDIF
!       IF(ASSOCIATED(SOURCE_FIELD)) THEN
!         SOURCE_INTERPOLATION_PARAMETERS=>equations%interpolation%sourceInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr
!       ENDIF

        CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,elementNumber, &
          & GEOMETRIC_INTERPOLATION_PARAMETERS,err,error,*999)
        IF(ASSOCIATED(FIBRE_FIELD)) THEN
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,elementNumber, &
            & FIBRE_INTERPOLATION_PARAMETERS,err,error,*999)
        END IF
        IF(ASSOCIATED(MATERIALS_FIELD)) THEN
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,elementNumber, &
            & MATERIALS_INTERPOLATION_PARAMETERS,err,error,*999)
!         IF(ASSOCIATED(DENSITY_INTERPOLATION_PARAMETERS)) THEN
!           CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,elementNumber, &
!             & DENSITY_INTERPOLATION_PARAMETERS,err,error,*999)
!         ENDIF
        ENDIF
        CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,elementNumber, &
          & DEPENDENT_INTERPOLATION_PARAMETERS,err,error,*999)
        IF(EQUATIONS%timeDependence/=EQUATIONS_STATIC) THEN
          CALL Field_InterpolationParametersElementGet(FIELD_PREVIOUS_VALUES_SET_TYPE,elementNumber, &
            & prevDependentInterpParameters,err,error,*999)
        ENDIF
        IF(DARCY_DEPENDENT) THEN
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,elementNumber, &
            & DARCY_DEPENDENT_INTERPOLATION_PARAMETERS,err,error,*999)
        ELSE IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_STANDARD_MONODOMAIN_ELASTICITY_SUBTYPE) THEN
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,elementNumber, &
            & INDEPENDENT_INTERPOLATION_PARAMETERS,err,error,*999)
        ENDIF
!       IF(ASSOCIATED(SOURCE_FIELD)) THEN
!         CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,elementNumber, &
!           & SOURCE_INTERPOLATION_PARAMETERS,err,error,*999)
!       END IF

        !Point interpolation pointer
        IF (EQUATIONS_SET_SUBTYPE == EQUATIONS_SET_REFERENCE_STATE_TRANSVERSE_GUCCIONE_SUBTYPE) THEN
          DEPENDENT_INTERPOLATED_POINT=>equations%interpolation%geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr
          DEPENDENT_INTERPOLATED_POINT_METRICS=>equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr
          GEOMETRIC_INTERPOLATED_POINT=>equations%interpolation%dependentInterpPoint(FIELD_VAR_TYPE)%ptr
          GEOMETRIC_INTERPOLATED_POINT_METRICS=>equations%interpolation%dependentInterpPointMetrics(FIELD_VAR_TYPE)%ptr
        ELSE
          GEOMETRIC_INTERPOLATED_POINT=>equations%interpolation%geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr
          GEOMETRIC_INTERPOLATED_POINT_METRICS=>equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr
          DEPENDENT_INTERPOLATED_POINT=>equations%interpolation%dependentInterpPoint(FIELD_VAR_TYPE)%ptr
          DEPENDENT_INTERPOLATED_POINT_METRICS=>equations%interpolation%dependentInterpPointMetrics(FIELD_VAR_TYPE)%ptr
        ENDIF
        IF(EQUATIONS%timeDependence/=EQUATIONS_STATIC) THEN
          prevDependentInterpPoint=>equations%interpolation%prevDependentInterpPoint(FIELD_VAR_TYPE)%ptr
          prevDependentInterpPointMetrics=>equations%interpolation%prevDependentInterpPointMetrics(FIELD_VAR_TYPE)%ptr
        ENDIF
        IF(ASSOCIATED(FIBRE_FIELD)) THEN
          FIBRE_INTERPOLATED_POINT=>equations%interpolation%fibreInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr
        END IF
        IF(ASSOCIATED(MATERIALS_FIELD)) THEN
          MATERIALS_INTERPOLATED_POINT=>equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr
          DENSITY_INTERPOLATED_POINT=>equations%interpolation%materialsInterpPoint(FIELD_V_VARIABLE_TYPE)%ptr
        ENDIF
        IF(DARCY_DEPENDENT) THEN
          DARCY_DEPENDENT_INTERPOLATED_POINT=>equations%interpolation%dependentInterpPoint(FIELD_V_VARIABLE_TYPE)%ptr
        ELSE IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_STANDARD_MONODOMAIN_ELASTICITY_SUBTYPE) THEN
          INDEPENDENT_INTERPOLATED_POINT=>equations%interpolation%independentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr
        ENDIF
        IF(ASSOCIATED(SOURCE_FIELD)) THEN
          SOURCE_INTERPOLATED_POINT=>equations%interpolation%sourceInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr
        ENDIF

        !SELECT: Compressible or incompressible cases, or poro multicompartment
        SELECT CASE(EQUATIONS_SET_SUBTYPE)
        ! ---------------------------------------------------------------
        CASE(EQUATIONS_SET_MOONEY_RIVLIN_ACTIVECONTRACTION_SUBTYPE, &
          & EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE)
          !Loop over gauss points and add residuals
          DO gauss_idx=1,DEPENDENT_NUMBER_OF_GAUSS_POINTS
            !Interpolate dependent, geometric, fibre and materials fields
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & DEPENDENT_INTERPOLATED_POINT,err,error,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(DEPENDENT_BASIS%NUMBER_OF_XI,DEPENDENT_INTERPOLATED_POINT_METRICS, &
              & err,error,*999)
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & GEOMETRIC_INTERPOLATED_POINT,err,error,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,GEOMETRIC_INTERPOLATED_POINT_METRICS, &
              & err,error,*999)
            IF(ASSOCIATED(FIBRE_FIELD)) THEN
              CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
                & FIBRE_INTERPOLATED_POINT,err,error,*999)
            END IF
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & MATERIALS_INTERPOLATED_POINT,err,error,*999)
            
            !Loop over geometric dependent basis functions.
            DO nh=1,numberOfDimensions
              MESH_COMPONENT_NUMBER=FIELD_VARIABLE%COMPONENTS(nh)%MESH_COMPONENT_NUMBER
              DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%ptr% &
                & TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS
              COMPONENT_QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
              DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                !Loop over derivative directions.
                DO mh=1,numberOfDimensions
                  SUM1=0.0_DP
                  DO mi=1,numberOfXi
                    SUM1=SUM1+DEPENDENT_INTERPOLATED_POINT_METRICS%DXI_DX(mi,mh)* &
                      & COMPONENT_QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(mi),gauss_idx)
                  ENDDO !mi
                  DPHIDZ(mh,ns,nh)=SUM1
                ENDDO !mh
              ENDDO !ns
            ENDDO !nh

            CALL FiniteElasticity_GaussDeformationGradientTensor(DEPENDENT_INTERPOLATED_POINT_METRICS, &
              & GEOMETRIC_INTERPOLATED_POINT_METRICS,FIBRE_INTERPOLATED_POINT,DZDNU,ERR,ERROR,*999)

            Jznu=DEPENDENT_INTERPOLATED_POINT_METRICS%JACOBIAN/GEOMETRIC_INTERPOLATED_POINT_METRICS%JACOBIAN
            JGW=DEPENDENT_INTERPOLATED_POINT_METRICS%JACOBIAN*DEPENDENT_QUADRATURE_SCHEME%GAUSS_WEIGHTS(gauss_idx)

           !Calculate the Cauchy stress tensor (in Voigt form) at the gauss point.
           CALL FINITE_ELASTICITY_GAUSS_STRESS_TENSOR(EQUATIONS_SET,DEPENDENT_INTERPOLATED_POINT, &
             & MATERIALS_INTERPOLATED_POINT,GEOMETRIC_INTERPOLATED_POINT,STRESS_TENSOR,DZDNU,Jznu, &
             & elementNumber,gauss_idx,err,error,*999)
            
            ! Convert from Voigt form to tensor form and multiply with Jacobian and Gauss weight.
            DO nh=1,numberOfDimensions
              DO mh=1,numberOfDimensions
                JGW_CAUCHY_TENSOR(mh,nh)=JGW*STRESS_TENSOR(TENSOR_TO_VOIGT3(mh,nh))
              ENDDO
            ENDDO

            !Now add up the residual terms
            mhs=0
            DO mh=1,numberOfDimensions
              MESH_COMPONENT_NUMBER=FIELD_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
              DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%ptr% &
                & TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1
                nonlinearMatrices%elementResidual%vector(mhs)=nonlinearMatrices%elementResidual%vector(mhs)+ &
                  & DOT_PRODUCT(DPHIDZ(1:numberOfDimensions,ms,mh),JGW_CAUCHY_TENSOR(1:numberOfDimensions,mh))
              ENDDO !ms
            ENDDO !mh

            JGW=GEOMETRIC_INTERPOLATED_POINT_METRICS%JACOBIAN*DEPENDENT_QUADRATURE_SCHEME%GAUSS_WEIGHTS(gauss_idx)

            !Hydrostatic pressure component
            MESH_COMPONENT_NUMBER=FIELD_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
            DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%ptr% &
              & TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS
            COMPONENT_QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
            TEMPTERM1=JGW*(Jznu-1.0_DP)
            IF(FIELD_VARIABLE%COMPONENTS(mh)%INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN !node based
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1 
                nonlinearMatrices%elementResidual%vector(mhs)=nonlinearMatrices%elementResidual%vector(mhs)+ &
                  & TEMPTERM1*COMPONENT_QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,gauss_idx)
              ENDDO
            ELSEIF(FIELD_VARIABLE%COMPONENTS(mh)%INTERPOLATION_TYPE==FIELD_ELEMENT_BASED_INTERPOLATION) THEN !element based
              mhs=mhs+1
              nonlinearMatrices%elementResidual%vector(mhs)=nonlinearMatrices%elementResidual%vector(mhs)+TEMPTERM1
            ENDIF

!           !Gravity loading term
!           IF(rhsVector%updateVector) THEN
!             IF(ASSOCIATED(SOURCE_FIELD)) THEN
!               CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
!                 & SOURCE_INTERPOLATED_POINT,err,error,*999)
!               IF(ASSOCIATED(DENSITY_INTERPOLATED_POINT)) THEN
!                 CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
!                   & DENSITY_INTERPOLATED_POINT,err,error,*999)
!                 DENSITY=DENSITY_INTERPOLATED_POINT%VALUES(1,NO_PART_DERIV)
!                 mhs=0
!                 DO mh=1,numberOfDimensions
!                   MESH_COMPONENT_NUMBER=FIELD_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
!                   DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%ptr% &
!                     & TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS
!                   COMPONENT_QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP( &
!                     & BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
!                   G_DENSITY_JGW=SOURCE_INTERPOLATED_POINT%VALUES(mh,NO_PART_DERIV)*DENSITY*JGW
!                   DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
!                     mhs=mhs+1
!                     rhsVector%elementVector%vector(mhs)=rhsVector%elementVector%vector(mhs)+ &
!                       & G_DENSITY_JGW*COMPONENT_QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,gauss_idx)
!                   ENDDO
!                 ENDDO
!               ENDIF
!             ENDIF
!           ENDIF
         ENDDO !gauss_idx


          !Call surface pressure term here: should only be executed if THIS element has surface pressure on it (direct or incremented)
          IF(DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BOUNDARY_ELEMENT.AND. &
            & TOTAL_NUMBER_OF_SURFACE_PRESSURE_CONDITIONS>0) THEN    ! 
            CALL FiniteElasticity_SurfacePressureResidualEvaluate(EQUATIONS_SET,elementNumber,var1,var2,err,error,*999)
          ENDIF

          !Scale factor adjustment
          IF(DEPENDENT_FIELD%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
            ! Following function is necessary, otherwise wrong face scale factors from function call to surface pressure residual are
            ! used.
            CALL Field_InterpolationParametersScaleFactorsElementGet(elementNumber, &
              & DEPENDENT_INTERPOLATION_PARAMETERS,err,error,*999) 
            mhs=0          
            DO mh=1,numberOfDimensions
              MESH_COMPONENT_NUMBER=FIELD_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
              DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%ptr% &
                & TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS
              !Loop over residual vector
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1
                nonlinearMatrices%elementResidual%vector(mhs)=nonlinearMatrices%elementResidual%vector(mhs)* &
                & DEPENDENT_INTERPOLATION_PARAMETERS%SCALE_FACTORS(ms,mh)
              ENDDO !ms
            ENDDO !mh
            IF(FIELD_VARIABLE%COMPONENTS(mh)%INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN !node based
              MESH_COMPONENT_NUMBER=FIELD_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
              DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%ptr% &
                & TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1 
                nonlinearMatrices%elementResidual%vector(mhs)=nonlinearMatrices%elementResidual%vector(mhs)* &
                & DEPENDENT_INTERPOLATION_PARAMETERS%SCALE_FACTORS(ms,mh)
              ENDDO
            ENDIF
          ENDIF

        ! ---------------------------------------------------------------
        CASE(EQUATIONS_SET_NO_SUBTYPE,EQUATIONS_SET_MEMBRANE_SUBTYPE, &
          & EQUATIONS_SET_STVENANT_KIRCHOFF_ACTIVECONTRACTION_SUBTYPE, &
          & EQUATIONS_SET_ISOTROPIC_EXPONENTIAL_SUBTYPE, &
          & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_EXPONENTIAL_SUBTYPE, &
          & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_GUCCIONE_SUBTYPE,EQUATIONS_SET_GUCCIONE_ACTIVECONTRACTION_SUBTYPE, &
          & EQUATIONS_SET_REFERENCE_STATE_TRANSVERSE_GUCCIONE_SUBTYPE, &
          & EQUATIONS_SET_ORTHOTROPIC_MATERIAL_COSTA_SUBTYPE, EQUATIONS_SET_ACTIVECONTRACTION_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE,EQUATIONS_SET_STANDARD_MONODOMAIN_ELASTICITY_SUBTYPE, &
          & EQUATIONS_SET_ORTHOTROPIC_MATERIAL_HOLZAPFEL_OGDEN_SUBTYPE,EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE, &
          & EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE,EQUATIONS_SET_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE, &
          & EQUATIONS_SET_ACTIVE_STRAIN_SUBTYPE,EQUATIONS_SET_MULTISCALE_ACTIVE_STRAIN_SUBTYPE, &
          & EQUATIONS_SET_1D3D_MONODOMAIN_ACTIVE_STRAIN_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE,EQUATIONS_SET_TRANSVERSE_ISOTROPIC_POLYNOMIAL_SUBTYPE, &
          & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_HUMPHREY_YIN_SUBTYPE, &
          & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_ACTIVE_SUBTYPE,EQUATIONS_SET_TRANS_ISOTROPIC_ACTIVE_TRANSITION_SUBTYPE, &
          & EQUATIONS_SET_ANISOTROPIC_POLYNOMIAL_SUBTYPE,EQUATIONS_SET_ANISOTROPIC_POLYNOMIAL_ACTIVE_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE, EQUATIONS_SET_HOLZAPFEL_OGDEN_ACTIVECONTRACTION_SUBTYPE) ! 4 dependent components

          growthValues=[1.0_DP,1.0_DP,1.0_DP]

          !Loop over gauss points and add residuals
          DO gauss_idx=1,DEPENDENT_NUMBER_OF_GAUSS_POINTS
            GAUSS_WEIGHT=DEPENDENT_QUADRATURE_SCHEME%GAUSS_WEIGHTS(gauss_idx)
              !Interpolate dependent, geometric, fibre and materials fields
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & DEPENDENT_INTERPOLATED_POINT,err,error,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(DEPENDENT_BASIS%NUMBER_OF_XI,DEPENDENT_INTERPOLATED_POINT_METRICS, &
              & err,error,*999)
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & GEOMETRIC_INTERPOLATED_POINT,err,error,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,GEOMETRIC_INTERPOLATED_POINT_METRICS, &
              & err,error,*999)
            IF(ASSOCIATED(FIBRE_FIELD)) THEN
              CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
                & FIBRE_INTERPOLATED_POINT,err,error,*999)
            END IF
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & MATERIALS_INTERPOLATED_POINT,err,error,*999)
            IF(DARCY_DEPENDENT) THEN
              CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
                & DARCY_DEPENDENT_INTERPOLATED_POINT,err,error,*999)
            ELSE IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_STANDARD_MONODOMAIN_ELASTICITY_SUBTYPE) THEN
              CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
                & INDEPENDENT_INTERPOLATED_POINT,err,error,*999)
            ENDIF

            !Calculate F=dZ/dNU, the deformation gradient tensor at the gauss point
            CALL FiniteElasticity_GaussDeformationGradientTensor(DEPENDENT_INTERPOLATED_POINT_METRICS, &
              & GEOMETRIC_INTERPOLATED_POINT_METRICS,FIBRE_INTERPOLATED_POINT,DZDNU,ERR,ERROR,*999)
            CALL Determinant(dZdNu,Jznu,err,error,*999)
            IF(Jznu<0.0_DP) THEN
              localWarning="Volume is negative for gauss point "//TRIM(NumberToVString(gauss_idx,"*",err,error))//&
                & " of element "//TRIM(NumberToVString(elementNumber,"*",err,error))//". det(F) = "// &
                & TRIM(NumberToVString(Jznu,"*",err,error))//"."
              CALL FlagWarning(localWarning,err,error,*999) 
            ENDIF

            Jzxi=DEPENDENT_INTERPOLATED_POINT_METRICS%JACOBIAN
            Jxxi=GEOMETRIC_INTERPOLATED_POINT_METRICS%JACOBIAN

            IF(DIAGNOSTICS1) THEN
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  elementNumber = ",elementNumber,err,error,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  gauss_idx = ",gauss_idx,err,error,*999)
            ENDIF

            !Calculate Jacobian of deformation.
            CALL FiniteElasticity_GaussGrowthTensor(EQUATIONS_SET,numberOfDimensions,dZdNu,growthValues,Fg,Fe,Jg,Je, &
              & err,error,*999)

            !Calculate strain tensors
            CALL FiniteElasticity_StrainTensor(Fe,C,f,Jznu,E,err,error,*999)

            !Calculate Sigma=1/Jznu.FTF', the Cauchy stress tensor at the gauss point
            CALL FINITE_ELASTICITY_GAUSS_CAUCHY_TENSOR(EQUATIONS_SET,DEPENDENT_INTERPOLATED_POINT, &
              & MATERIALS_INTERPOLATED_POINT,GEOMETRIC_INTERPOLATED_POINT,DARCY_DEPENDENT_INTERPOLATED_POINT, &
              & INDEPENDENT_INTERPOLATED_POINT,cauchyTensor,Jznu,DZDNU,elementNumber,gauss_idx,err,error,*999)

            IF(DIAGNOSTICS1) THEN
              CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"",err,error,*999)
              CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Stress tensors:",err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Hydrostatic pressure = ",P,err,error,*999)
              CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Second Piola-Kirchoff stress tensor:",err,error,*999)
              CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
                & 3,3,piolaTensor,WRITE_STRING_MATRIX_NAME_AND_INDICES, '("    T','(",I1,",:)','     :",3(X,E13.6))', &
                & '(12X,3(X,E13.6))',err,error,*999)
              CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Cauchy stress tensor:",err,error,*999)
              CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
                & 3,3,cauchyTensor,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("    sigma','(",I1,",:)',' :",3(X,E13.6))', &
                & '(12X,3(X,E13.6))',err,error,*999)
            ENDIF

            IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE) THEN
              !Parameters settings for coupled elasticity Darcy INRIA model:
              CALL GET_DARCY_FINITE_ELASTICITY_PARAMETERS(DARCY_RHO_0_F,Mfact,bfact,p0fact,err,error,*999)
              DARCY_MASS_INCREASE = DARCY_DEPENDENT_INTERPOLATED_POINT%VALUES(4,NO_PART_DERIV) 
              DARCY_VOL_INCREASE = DARCY_MASS_INCREASE / DARCY_RHO_0_F
            ENDIF

            !For membrane theory in 3D space, the final equation is multiplied by thickness. Default to unit thickness if equation set subtype is not membrane
            THICKNESS = 1.0_DP
            IF(EQUATIONS_SET_SUBTYPE == EQUATIONS_SET_MEMBRANE_SUBTYPE) THEN
              IF(numberOfDimensions == 3) THEN
                THICKNESS = MATERIALS_INTERPOLATED_POINT%VALUES(MATERIALS_INTERPOLATED_POINT%INTERPOLATION_PARAMETERS% &
                  & FIELD_VARIABLE%NUMBER_OF_COMPONENTS,1)
              ENDIF
            ENDIF

            !Calculate the combined Jacobian
            JGW=Jzxi*DEPENDENT_QUADRATURE_SCHEME%GAUSS_WEIGHTS(gauss_idx)
            
            !Loop over geometric dependent basis functions and evaluate dPhidZ.
            DO nh=1,numberOfDimensions
              MESH_COMPONENT_NUMBER=FIELD_VARIABLE%COMPONENTS(nh)%MESH_COMPONENT_NUMBER
              DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%ptr% &
                & TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS
              COMPONENT_QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
              DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                !Loop over derivative directions.
                DO mh=1,numberOfDimensions
                  SUM1=0.0_DP
                  DO mi=1,DEPENDENT_BASIS%NUMBER_OF_XI
                    SUM1=SUM1+DEPENDENT_INTERPOLATED_POINT_METRICS%DXI_DX(mi,mh)* &
                      & COMPONENT_QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(mi),gauss_idx)
                  ENDDO !mi
                  DPHIDZ(mh,ns,nh)=SUM1
                ENDDO !mh
              ENDDO !ns
            ENDDO !nh
            
            !Now add up the residual terms
            element_dof_idx=0
            DO component_idx=1,numberOfDimensions
              DEPENDENT_COMPONENT_INTERPOLATION_TYPE=DEPENDENT_FIELD%VARIABLES(var1)%COMPONENTS(component_idx)%INTERPOLATION_TYPE
              IF(DEPENDENT_COMPONENT_INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN !node based
                DEPENDENT_BASIS=>DEPENDENT_FIELD%VARIABLES(var1)%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY% &
                  & ELEMENTS%ELEMENTS(elementNumber)%BASIS
                NUMBER_OF_FIELD_COMPONENT_INTERPOLATION_PARAMETERS=DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                DO parameter_idx=1,NUMBER_OF_FIELD_COMPONENT_INTERPOLATION_PARAMETERS
                  element_dof_idx=element_dof_idx+1
                  nonlinearMatrices%elementResidual%vector(element_dof_idx)= &
                    & nonlinearMatrices%elementResidual%vector(element_dof_idx)+ &
                    & JGW*DOT_PRODUCT(DPhiDZ(1:numberOfDimensions,parameter_idx,component_idx), &
                    & cauchyTensor(1:numberOfDimensions,component_idx))
                ENDDO ! parameter_idx (residual vector loop)
              ELSEIF(DEPENDENT_COMPONENT_INTERPOLATION_TYPE==FIELD_ELEMENT_BASED_INTERPOLATION) THEN
                !Will probably never be used
                CALL FlagError("Finite elasticity with element based interpolation is not implemented.",err,error,*999)
              ENDIF
            ENDDO ! component_idx

            !Hydrostatic pressure component (skip for membrane problems)
            IF (EQUATIONS_SET_SUBTYPE /= EQUATIONS_SET_MEMBRANE_SUBTYPE) THEN
              HYDROSTATIC_PRESSURE_COMPONENT=DEPENDENT_FIELD%VARIABLES(var1)%NUMBER_OF_COMPONENTS
              DEPENDENT_COMPONENT_INTERPOLATION_TYPE=DEPENDENT_FIELD%VARIABLES(var1)%COMPONENTS(HYDROSTATIC_PRESSURE_COMPONENT)% &
                & INTERPOLATION_TYPE
              IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE) THEN
                TEMPTERM1=GAUSS_WEIGHT*(Jzxi-(Jg-DARCY_VOL_INCREASE)*Jxxi)
              ELSE
                TEMPTERM1=GAUSS_WEIGHT*(Jzxi/Jxxi - 1.0_DP)*Jxxi
              ENDIF
              IF(DEPENDENT_COMPONENT_INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN !node based
                COMPONENT_BASIS=>DEPENDENT_FIELD%VARIABLES(var1)%COMPONENTS(HYDROSTATIC_PRESSURE_COMPONENT)%DOMAIN% &
                  & TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS
                COMPONENT_QUADRATURE_SCHEME=>COMPONENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
                NUMBER_OF_FIELD_COMPONENT_INTERPOLATION_PARAMETERS=COMPONENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                DO parameter_idx=1,NUMBER_OF_FIELD_COMPONENT_INTERPOLATION_PARAMETERS
                  element_dof_idx=element_dof_idx+1 
                  nonlinearMatrices%elementResidual%vector(element_dof_idx)= &
                    & nonlinearMatrices%elementResidual%vector(element_dof_idx)+ &
                    & COMPONENT_QUADRATURE_SCHEME%GAUSS_BASIS_FNS(parameter_idx,1,gauss_idx)*TEMPTERM1
                ENDDO
              ELSEIF(DEPENDENT_COMPONENT_INTERPOLATION_TYPE==FIELD_ELEMENT_BASED_INTERPOLATION) THEN !element based
                element_dof_idx=element_dof_idx+1                
                nonlinearMatrices%elementResidual%vector(element_dof_idx)= &
                  & nonlinearMatrices%elementResidual%vector(element_dof_idx)+TEMPTERM1              
              ENDIF
            ENDIF
          ENDDO !gauss_idx

          !Call surface pressure term here: should only be executed if THIS element has surface pressure on it (direct or incremented)
          IF(DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BOUNDARY_ELEMENT.AND. &
            & TOTAL_NUMBER_OF_SURFACE_PRESSURE_CONDITIONS>0) THEN    ! 
            CALL FiniteElasticity_SurfacePressureResidualEvaluate(EQUATIONS_SET,elementNumber,var1,var2,err,error,*999)
          ENDIF

          ! ---------------------------------------------------------------
        CASE(EQUATIONS_SET_MR_AND_GROWTH_LAW_IN_CELLML_SUBTYPE)

          !Loop over gauss points and add residuals
          DO gauss_idx=1,DEPENDENT_NUMBER_OF_GAUSS_POINTS
            
            IF(DIAGNOSTICS1) THEN
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Element number = ",elementNumber,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Gauss index  = ",gauss_idx,err,error,*999)
            ENDIF

            GAUSS_WEIGHT=DEPENDENT_QUADRATURE_SCHEME%GAUSS_WEIGHTS(gauss_idx)
            !Interpolate dependent, geometric, fibre and materials fields
            CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & DEPENDENT_INTERPOLATED_POINT,err,error,*999)
            CALL Field_InterpolatedPointMetricsCalculate(DEPENDENT_BASIS%NUMBER_OF_XI,DEPENDENT_INTERPOLATED_POINT_METRICS, &
              & err,error,*999)
            CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & GEOMETRIC_INTERPOLATED_POINT,err,error,*999)
            CALL Field_InterpolatedPointMetricsCalculate(GEOMETRIC_BASIS%NUMBER_OF_XI,GEOMETRIC_INTERPOLATED_POINT_METRICS, &
              & err,error,*999)
            IF(ASSOCIATED(FIBRE_FIELD)) THEN
              CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
                & FIBRE_INTERPOLATED_POINT,err,error,*999)
            ENDIF
            CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & MATERIALS_INTERPOLATED_POINT,err,error,*999)
            CALL Field_ParameterSetGetLocalGaussPoint(dependent_Field,FIELD_U3_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & gauss_idx,elementNumber,1,growthValues(1),err,error,*999)
            IF(numberofDimensions>1) THEN
              CALL Field_ParameterSetGetLocalGaussPoint(dependent_Field,FIELD_U3_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & gauss_idx,elementNumber,2,growthValues(2),err,error,*999)
              IF(numberOfDimensions>2) THEN
                CALL Field_ParameterSetGetLocalGaussPoint(dependent_Field,FIELD_U3_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & gauss_idx,elementNumber,3,growthValues(3),err,error,*999)
              ENDIF
            ENDIF

            !Calculate F=dZ/dNU, the deformation gradient tensor at the gauss point
            CALL FiniteElasticity_GaussDeformationGradientTensor(DEPENDENT_INTERPOLATED_POINT_METRICS, &
              & GEOMETRIC_INTERPOLATED_POINT_METRICS,FIBRE_INTERPOLATED_POINT,DZDNU,ERR,ERROR,*999)

            Jxxi=GEOMETRIC_INTERPOLATED_POINT_METRICS%JACOBIAN

            Jzxi=DEPENDENT_INTERPOLATED_POINT_METRICS%JACOBIAN
            
            HYDROSTATIC_PRESSURE_COMPONENT=DEPENDENT_INTERPOLATED_POINT%INTERPOLATION_PARAMETERS%FIELD_VARIABLE% &
              & NUMBER_OF_COMPONENTS
            P=DEPENDENT_INTERPOLATED_POINT%VALUES(HYDROSTATIC_PRESSURE_COMPONENT,1)
            
            CALL FiniteElasticity_GaussGrowthTensor(EQUATIONS_SET,numberOfDimensions,dZdNu,growthValues,Fg,Fe,Jg,Je, &
              & err,error,*999)
             
            CALL FiniteElasticity_StrainTensor(Fe,C,f,Jznu,E,err,error,*999)
 
            !Calculate the Cauchy stress tensor (in Voigt form) at the gauss point.
            CALL FINITE_ELASTICITY_GAUSS_STRESS_TENSOR(EQUATIONS_SET,DEPENDENT_INTERPOLATED_POINT, &
             & MATERIALS_INTERPOLATED_POINT,GEOMETRIC_INTERPOLATED_POINT,STRESS_TENSOR,Fe,Jznu, &
             & elementNumber,gauss_idx,err,error,*999)
            ! Convert from Voigt form to tensor form and multiply with Jacobian and Gauss weight.
            
            JGW=Jzxi*DEPENDENT_QUADRATURE_SCHEME%GAUSS_WEIGHTS(gauss_idx)
            DO nh=1,numberOfDimensions
              DO mh=1,numberOfDimensions
                JGW_CAUCHY_TENSOR(mh,nh)=JGW*STRESS_TENSOR(TENSOR_TO_VOIGT3(mh,nh))
              ENDDO
            ENDDO
            
            IF(DIAGNOSTICS1) THEN
              CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"",err,error,*999)
              CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Stress tensors:",err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Hydrostatic pressure = ",P,err,error,*999)
              CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Second Piola-Kirchoff stress tensor:",err,error,*999)
              CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
                & 3,3,piolaTensor,WRITE_STRING_MATRIX_NAME_AND_INDICES, '("    T','(",I1,",:)','     :",3(X,E13.6))', &
                & '(12X,3(X,E13.6))',err,error,*999)
              CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Cauchy stress tensor:",err,error,*999)
              CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
                & 3,3,cauchyTensor,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("    sigma','(",I1,",:)',' :",3(X,E13.6))', &
                & '(12X,3(X,E13.6))',err,error,*999)
            ENDIF

            !For membrane theory in 3D space, the final equation is multiplied by thickness. Default to unit thickness if equation set subtype is not membrane
            !!TODO Maybe have the thickness as a component in the equations set field. Yes, as we don't need a materials field for CellML constituative laws.
            THICKNESS = 1.0_DP
            IF(EQUATIONS_SET_SUBTYPE == EQUATIONS_SET_MEMBRANE_SUBTYPE) THEN
              IF(numberOfDimensions == 3) THEN
                IF(ASSOCIATED(MATERIALS_FIELD)) THEN
                  CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
                    & MATERIALS_INTERPOLATED_POINT,err,error,*999)
                  THICKNESS = MATERIALS_INTERPOLATED_POINT%VALUES(MATERIALS_INTERPOLATED_POINT%INTERPOLATION_PARAMETERS% &
                    & FIELD_VARIABLE%NUMBER_OF_COMPONENTS,1)
                ENDIF
              ENDIF
            ENDIF

            !!Now add up the residual terms

            !Loop over geometric dependent basis functions.
            DO nh=1,numberOfDimensions
              MESH_COMPONENT_NUMBER=FIELD_VARIABLE%COMPONENTS(nh)%MESH_COMPONENT_NUMBER
              DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%ptr% &
                & TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS
              COMPONENT_QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
              DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                !Loop over derivative directions.
                DO mh=1,numberOfDimensions
                  SUM1=0.0_DP
                  DO mi=1,numberOfXi
                    SUM1=SUM1+DEPENDENT_INTERPOLATED_POINT_METRICS%DXI_DX(mi,mh)* &
                      & COMPONENT_QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(mi),gauss_idx)
                  ENDDO !mi
                  DPHIDZ(mh,ns,nh)=SUM1
                ENDDO !mh
              ENDDO !ns
            ENDDO !nh
            !Now add up the residual terms
            mhs=0
            DO mh=1,numberOfDimensions
              MESH_COMPONENT_NUMBER=FIELD_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
              DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%ptr% &
                & TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1
                nonlinearMatrices%elementResidual%vector(mhs)=nonlinearMatrices%elementResidual%vector(mhs)+ &
                  & DOT_PRODUCT(DPHIDZ(1:numberOfDimensions,ms,mh),JGW_CAUCHY_TENSOR(1:numberOfDimensions,mh))
              ENDDO !ms
            ENDDO !mh
            
            !Hydrostatic pressure component (skip for membrane problems)
            IF (EQUATIONS_SET_SUBTYPE /= EQUATIONS_SET_MEMBRANE_SUBTYPE) THEN
              IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_REFERENCE_STATE_TRANSVERSE_GUCCIONE_SUBTYPE) THEN
                HYDROSTATIC_PRESSURE_COMPONENT=GEOMETRIC_FIELD%VARIABLES(var1)%NUMBER_OF_COMPONENTS
                DEPENDENT_COMPONENT_INTERPOLATION_TYPE=GEOMETRIC_FIELD%VARIABLES(var1)%COMPONENTS( &
                  & HYDROSTATIC_PRESSURE_COMPONENT)%INTERPOLATION_TYPE
              ELSE
                HYDROSTATIC_PRESSURE_COMPONENT=DEPENDENT_FIELD%VARIABLES(var1)%NUMBER_OF_COMPONENTS
                DEPENDENT_COMPONENT_INTERPOLATION_TYPE=DEPENDENT_FIELD%VARIABLES(var1)%COMPONENTS( &
                  & HYDROSTATIC_PRESSURE_COMPONENT)%INTERPOLATION_TYPE
              ENDIF
              IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE) THEN
                TEMPTERM1=GAUSS_WEIGHT*Jxxi*(Jznu-(Jg-DARCY_VOL_INCREASE))
              ELSE
                TEMPTERM1=GAUSS_WEIGHT*Jxxi*(Jznu-Jg)
              ENDIF            
              IF(DEPENDENT_COMPONENT_INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN !node based
                IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_REFERENCE_STATE_TRANSVERSE_GUCCIONE_SUBTYPE) THEN
                  COMPONENT_BASIS=>GEOMETRIC_FIELD%VARIABLES(var1)%COMPONENTS(HYDROSTATIC_PRESSURE_COMPONENT)%DOMAIN% &
                    & TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS
                ELSE
                  COMPONENT_BASIS=>DEPENDENT_FIELD%VARIABLES(var1)%COMPONENTS(HYDROSTATIC_PRESSURE_COMPONENT)%DOMAIN% &
                    & TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS
                ENDIF
                COMPONENT_QUADRATURE_SCHEME=>COMPONENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
                NUMBER_OF_FIELD_COMPONENT_INTERPOLATION_PARAMETERS=COMPONENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                DO parameter_idx=1,NUMBER_OF_FIELD_COMPONENT_INTERPOLATION_PARAMETERS
                  mhs=mhs+1 
                  nonlinearMatrices%elementResidual%vector(mhs)= &
                    & nonlinearMatrices%elementResidual%vector(mhs)+ &
                    & COMPONENT_QUADRATURE_SCHEME%GAUSS_BASIS_FNS(parameter_idx,1,gauss_idx)*TEMPTERM1                  
                ENDDO
              ELSEIF(DEPENDENT_COMPONENT_INTERPOLATION_TYPE==FIELD_ELEMENT_BASED_INTERPOLATION) THEN !element based
                mhs=mhs+1
                nonlinearMatrices%elementResidual%vector(mhs)= &
                  & nonlinearMatrices%elementResidual%vector(mhs)+TEMPTERM1
              ENDIF
            ENDIF
          ENDDO !gauss_idx

          !Call surface pressure term here: should only be executed if THIS element has surface pressure on it (direct or incremented)
          IF(DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BOUNDARY_ELEMENT.AND. &
            & TOTAL_NUMBER_OF_SURFACE_PRESSURE_CONDITIONS>0) THEN    ! 
            CALL FiniteElasticity_SurfacePressureResidualEvaluate(EQUATIONS_SET,elementNumber,var1,var2,err,error,*999)
          ENDIF
          
          ! ---------------------------------------------------------------
        CASE(EQUATIONS_SET_GROWTH_LAW_IN_CELLML_SUBTYPE)
          CALL IdentityMatrix(ITens,err,error,*999)
          time=EQUATIONS_SET%currentTime         
          !Loop over gauss points and add residuals
          DO gauss_idx=1,DEPENDENT_NUMBER_OF_GAUSS_POINTS
            
            IF(DIAGNOSTICS1) THEN
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Element number = ",elementNumber,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Gauss index  = ",gauss_idx,err,error,*999)
            ENDIF

            GAUSS_WEIGHT=DEPENDENT_QUADRATURE_SCHEME%GAUSS_WEIGHTS(gauss_idx)
            !Interpolate dependent, geometric, fibre and materials fields
            CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & DEPENDENT_INTERPOLATED_POINT,err,error,*999)
            CALL Field_InterpolatedPointMetricsCalculate(DEPENDENT_BASIS%NUMBER_OF_XI,DEPENDENT_INTERPOLATED_POINT_METRICS, &
              & err,error,*999)
            CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & GEOMETRIC_INTERPOLATED_POINT,err,error,*999)
            CALL Field_InterpolatedPointMetricsCalculate(GEOMETRIC_BASIS%NUMBER_OF_XI,GEOMETRIC_INTERPOLATED_POINT_METRICS, &
              & err,error,*999)
            IF(ASSOCIATED(FIBRE_FIELD)) THEN
              CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
                & FIBRE_INTERPOLATED_POINT,err,error,*999)
            ENDIF
            CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & MATERIALS_INTERPOLATED_POINT,err,error,*999)

            q=MATERIALS_INTERPOLATED_POINT%VALUES(1,NO_PART_DERIV)
            k1=MATERIALS_INTERPOLATED_POINT%VALUES(2,NO_PART_DERIV)
            k2=MATERIALS_INTERPOLATED_POINT%VALUES(3,NO_PART_DERIV)
            k3=MATERIALS_INTERPOLATED_POINT%VALUES(4,NO_PART_DERIV)
            k4=MATERIALS_INTERPOLATED_POINT%VALUES(5,NO_PART_DERIV)
            b1=MATERIALS_INTERPOLATED_POINT%VALUES(6,NO_PART_DERIV)
            b2=MATERIALS_INTERPOLATED_POINT%VALUES(7,NO_PART_DERIV)
            b3=MATERIALS_INTERPOLATED_POINT%VALUES(8,NO_PART_DERIV)
            b4=MATERIALS_INTERPOLATED_POINT%VALUES(9,NO_PART_DERIV)
            b5=MATERIALS_INTERPOLATED_POINT%VALUES(10,NO_PART_DERIV)
            b6=MATERIALS_INTERPOLATED_POINT%VALUES(11,NO_PART_DERIV)
            mu0=MATERIALS_INTERPOLATED_POINT%VALUES(12,NO_PART_DERIV)
            rho0=MATERIALS_INTERPOLATED_POINT%VALUES(13,NO_PART_DERIV)

            A1g=MATERIALS_INTERPOLATED_POINT%VALUES(14,NO_PART_DERIV)
            B0g=MATERIALS_INTERPOLATED_POINT%VALUES(15,NO_PART_DERIV)
            B1g=MATERIALS_INTERPOLATED_POINT%VALUES(16,NO_PART_DERIV)
            B2g=MATERIALS_INTERPOLATED_POINT%VALUES(17,NO_PART_DERIV)
            alphag=MATERIALS_INTERPOLATED_POINT%VALUES(18,NO_PART_DERIV)
            chi1g=MATERIALS_INTERPOLATED_POINT%VALUES(19,NO_PART_DERIV)
            chi2g=MATERIALS_INTERPOLATED_POINT%VALUES(20,NO_PART_DERIV)
            phig=MATERIALS_INTERPOLATED_POINT%VALUES(21,NO_PART_DERIV)
            R1g=MATERIALS_INTERPOLATED_POINT%VALUES(22,NO_PART_DERIV)
            R2g=MATERIALS_INTERPOLATED_POINT%VALUES(23,NO_PART_DERIV)
            Lg=MATERIALS_INTERPOLATED_POINT%VALUES(24,NO_PART_DERIV)

            x=GEOMETRIC_INTERPOLATED_POINT%VALUES(1,NO_PART_DERIV)
            y=GEOMETRIC_INTERPOLATED_POINT%VALUES(2,NO_PART_DERIV)
            z=GEOMETRIC_INTERPOLATED_POINT%VALUES(3,NO_PART_DERIV)

            r=SQRT(x*x+y*y)
            theta=ATAN2(x,y)
              
            !Calculate F=dZ/dNU, the deformation gradient tensor at the gauss point
             CALL MatrixProduct(DEPENDENT_INTERPOLATED_POINT_METRICS%DX_DXI(1:3,1:3), &
              & GEOMETRIC_INTERPOLATED_POINT_METRICS%DXI_DX(1:3,1:3),dZdNu(1:3,1:3),err,error,*999)
            
            Jxxi=GEOMETRIC_INTERPOLATED_POINT_METRICS%JACOBIAN
            
            Jzxi=DEPENDENT_INTERPOLATED_POINT_METRICS%JACOBIAN

            Dg22=A1g*(1.0_DP-EXP(-1.0_DP*alphag*MacaulayBracket(z/lg-chi1g)**2.0_DP)- &
              & EXP(-1.0_DP*alphag*MacaulayBracket(chi2g-z/lg)**2.0_DP))+ &
              & (B0g+(B1g-B0g)*(1.0_DP-EXP(-1.0_DP*alphag*MacaulayBracket(z/lg-chi1g)**2.0_DP))+ &
              & (B2g-B1g)*EXP(-1.0_DP*alphag*MacaulayBracket(chi2g-z/lg)**2.0_DP))*r*SIN(theta-phig)/R2g

            CALL IdentityMatrix(Fg,err,error,*999)
            Fg(3,3)=EXP(Dg22*time)
            
            !Calculate inverse growth deformation tensor, Fg^-1, Jg 
            CALL Invert(Fg,FgInv,Jg,err,error,*999)
            !Calculate elastic deformation tensor, Fe=F.(Fg)^-1.       
            CALL MatrixProduct(dZdNu,FgInv,Fe,err,error,*999)
            CALL Determinant(Fe,Je,err,error,*999)
            
            !Compute the distortional tensor
            CALL MatrixProductTranspose(Fe,Fe,Beprime,err,error,*999)
            Beprime=Beprime*Je**(2.0_DP/3.0_DP)
            CALL DoubleDotProduct(Beprime,Itens,alpha1,err,error,*999)
            bigQ=0.5_DP*(k1*(0.5_DP*(Je*Je-1.0_DP)-LOG(Je))+k2*(alpha1-3.0_DP))
            delWdelQ=mu0*bigQ*EXP(q*bigQ)
            delQdelJe=k1*(Je*Je-1)/(2.0_DP*Je)
            delQdelAlpha=k2
            delWdelJe=delWdelQ*delQdelJe
            delWdelAlpha=delWdelQ*delQdelAlpha

            cauchyTensor=delWdelJe*Itens+(2.0_DP/Je)*delWdelAlpha*(Beprime-(alpha1/3.0_DP)*Itens)

            IF(DIAGNOSTICS1) THEN
              CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"",err,error,*999)
              CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Stress tensors:",err,error,*999)
              CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Cauchy stress tensor:",err,error,*999)
              CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
                & 3,3,cauchyTensor,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("    sigma','(",I1,",:)',' :",3(X,E13.6))', &
                & '(12X,3(X,E13.6))',err,error,*999)
            ENDIF
            
            !Loop over geometric dependent basis functions.
            DO nh=1,numberOfDimensions
              MESH_COMPONENT_NUMBER=FIELD_VARIABLE%COMPONENTS(nh)%MESH_COMPONENT_NUMBER
              DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%ptr% &
                & TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS
              COMPONENT_QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
              DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                !Loop over derivative directions.
                DO mh=1,numberOfDimensions
                  SUM1=0.0_DP
                  DO mi=1,numberOfXi
                    SUM1=SUM1+DEPENDENT_INTERPOLATED_POINT_METRICS%DXI_DX(mi,mh)* &
                      & COMPONENT_QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(mi),gauss_idx)
                  ENDDO !mi
                  DPHIDZ(mh,ns,nh)=SUM1
                ENDDO !mh
              ENDDO !ns
            ENDDO !nh
            JGW=Jzxi*DEPENDENT_QUADRATURE_SCHEME%GAUSS_WEIGHTS(gauss_idx)
            !Now add up the residual terms
            mhs=0
            DO mh=1,numberOfDimensions
              MESH_COMPONENT_NUMBER=FIELD_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
              DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%ptr% &
                & TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1
                nonlinearMatrices%elementResidual%vector(mhs)=nonlinearMatrices%elementResidual%vector(mhs)+ &
                  & JGW*DOT_PRODUCT(DPHIDZ(1:numberOfDimensions,ms,mh),cauchyTensor(1:numberOfDimensions,mh))
              ENDDO !ms
            ENDDO !mh

          ENDDO !gauss_idx
          
          !Call surface pressure term here: should only be executed if THIS element has surface pressure on it (direct or incremented)
          IF(DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BOUNDARY_ELEMENT.AND. &
            & TOTAL_NUMBER_OF_SURFACE_PRESSURE_CONDITIONS>0) THEN    ! 
            CALL FiniteElasticity_SurfacePressureResidualEvaluate(EQUATIONS_SET,elementNumber,var1,var2,err,error,*999)
          ENDIF

        CASE(EQUATIONS_SET_CONSTITUTIVE_LAW_IN_CELLML_EVALUATE_SUBTYPE, &
          & EQUATIONS_SET_CONSTITUTIVE_AND_GROWTH_LAW_IN_CELLML_SUBTYPE)

          !Loop over gauss points and add residuals
          DO gauss_idx=1,DEPENDENT_NUMBER_OF_GAUSS_POINTS
            
            IF(DIAGNOSTICS1) THEN
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Element number = ",elementNumber,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Gauss index  = ",gauss_idx,err,error,*999)
            ENDIF

            GAUSS_WEIGHT=DEPENDENT_QUADRATURE_SCHEME%GAUSS_WEIGHTS(gauss_idx)
            !Interpolate dependent, geometric, fibre and materials fields
            CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & DEPENDENT_INTERPOLATED_POINT,err,error,*999)
            CALL Field_InterpolatedPointMetricsCalculate(DEPENDENT_BASIS%NUMBER_OF_XI,DEPENDENT_INTERPOLATED_POINT_METRICS, &
              & err,error,*999)
            CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & GEOMETRIC_INTERPOLATED_POINT,err,error,*999)
            CALL Field_InterpolatedPointMetricsCalculate(GEOMETRIC_BASIS%NUMBER_OF_XI,GEOMETRIC_INTERPOLATED_POINT_METRICS, &
              & err,error,*999)
            IF(ASSOCIATED(FIBRE_FIELD)) THEN
              CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
                & FIBRE_INTERPOLATED_POINT,err,error,*999)
            ENDIF
            IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTITUTIVE_AND_GROWTH_LAW_IN_CELLML_SUBTYPE) THEN
              CALL Field_ParameterSetGetLocalGaussPoint(dependent_Field,FIELD_U3_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & gauss_idx,elementNumber,1,growthValues(1),err,error,*999)
              IF(numberofDimensions>1) THEN
                CALL Field_ParameterSetGetLocalGaussPoint(dependent_Field,FIELD_U3_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & gauss_idx,elementNumber,2,growthValues(2),err,error,*999)
                IF(numberOfDimensions>2) THEN
                  CALL Field_ParameterSetGetLocalGaussPoint(dependent_Field,FIELD_U3_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                    & gauss_idx,elementNumber,3,growthValues(3),err,error,*999)
                ENDIF
              ENDIF
            ELSE
              growthValues=[1.0_DP,1.0_DP,1.0_DP]
            ENDIF

            !Calculate F=dZ/dNU, the deformation gradient tensor at the gauss point
            CALL FiniteElasticity_GaussDeformationGradientTensor(DEPENDENT_INTERPOLATED_POINT_METRICS, &
              & GEOMETRIC_INTERPOLATED_POINT_METRICS,FIBRE_INTERPOLATED_POINT,dZdNu,err,error,*999)

            Jxxi=GEOMETRIC_INTERPOLATED_POINT_METRICS%JACOBIAN
            
            Jzxi=DEPENDENT_INTERPOLATED_POINT_METRICS%JACOBIAN
            
            HYDROSTATIC_PRESSURE_COMPONENT=DEPENDENT_INTERPOLATED_POINT%INTERPOLATION_PARAMETERS%FIELD_VARIABLE% &
              & NUMBER_OF_COMPONENTS
            P=DEPENDENT_INTERPOLATED_POINT%VALUES(HYDROSTATIC_PRESSURE_COMPONENT,1)
            
            CALL FiniteElasticity_GaussGrowthTensor(EQUATIONS_SET,numberOfDimensions,dZdNu,growthValues,Fg,Fe,Jg,Je, &
              & err,error,*999)
             
            CALL FiniteElasticity_StrainTensor(Fe,C,f,Jznu,E,err,error,*999)
 
            !Get the stress field!!!
            IF(numberOfDimensions==3) THEN
              CALL Field_ParameterSetGetLocalGaussPoint(DEPENDENT_FIELD,FIELD_U2_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & gauss_idx,elementNumber,1,piolaTensor(1,1),err,error,*999)
              CALL Field_ParameterSetGetLocalGaussPoint(DEPENDENT_FIELD,FIELD_U2_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & gauss_idx,elementNumber,2,piolaTensor(1,2),err,error,*999)
              CALL Field_ParameterSetGetLocalGaussPoint(DEPENDENT_FIELD,FIELD_U2_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & gauss_idx,elementNumber,3,piolaTensor(1,3),err,error,*999)
              CALL Field_ParameterSetGetLocalGaussPoint(DEPENDENT_FIELD,FIELD_U2_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & gauss_idx,elementNumber,4,piolaTensor(2,2),err,error,*999)
              CALL Field_ParameterSetGetLocalGaussPoint(DEPENDENT_FIELD,FIELD_U2_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & gauss_idx,elementNumber,5,piolaTensor(2,3),err,error,*999)
              CALL Field_ParameterSetGetLocalGaussPoint(DEPENDENT_FIELD,FIELD_U2_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & gauss_idx,elementNumber,6,piolaTensor(3,3),err,error,*999)
              !CellML computes the deviatoric stress. Add the volumetric component!
              piolaTensor(1,1)=piolaTensor(1,1)+P*f(1,1)
              piolaTensor(2,2)=piolaTensor(2,2)+P*f(2,2)
              piolaTensor(3,3)=piolaTensor(3,3)+P*f(3,3)
              piolaTensor(1,2)=piolaTensor(1,2)+P*f(1,2)
              piolaTensor(1,3)=piolaTensor(1,3)+P*f(1,3)
              piolaTensor(2,3)=piolaTensor(2,3)+P*f(2,3)
              piolaTensor(2,1)=piolaTensor(1,2)
              piolaTensor(3,1)=piolaTensor(1,3)
              piolaTensor(3,2)=piolaTensor(2,3)
            ELSE IF(numberOfDimensions==2) THEN
              CALL Field_ParameterSetGetLocalGaussPoint(DEPENDENT_FIELD,FIELD_U2_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & gauss_idx,elementNumber,1,piolaTensor(1,1),err,error,*999)
              CALL Field_ParameterSetGetLocalGaussPoint(DEPENDENT_FIELD,FIELD_U2_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & gauss_idx,elementNumber,2,piolaTensor(1,2),err,error,*999)
              CALL Field_ParameterSetGetLocalGaussPoint(DEPENDENT_FIELD,FIELD_U2_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & gauss_idx,elementNumber,3,piolaTensor(2,2),err,error,*999)
              !CellML computes the deviatoric stress. Add the volumetric component!
              piolaTensor(1,1)=piolaTensor(1,1)+P*f(1,1)
              piolaTensor(2,2)=piolaTensor(2,2)+P*f(2,2)
              piolaTensor(1,2)=piolaTensor(1,2)+P*f(1,2)
              piolaTensor(2,1)=piolaTensor(1,2)
            ELSE
              CALL Field_ParameterSetGetLocalGaussPoint(DEPENDENT_FIELD,FIELD_U2_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & gauss_idx,elementNumber,1,piolaTensor(1,1),err,error,*999)
              piolaTensor(1,1)=piolaTensor(1,1)+P*f(1,1)
            ENDIF

            !Compute the Kirchoff stress tensor by pushing the 2nd Piola Kirchoff stress tensor forward \tau = F.S.F^T
            CALL MatrixProduct(Fe,piolaTensor,temp,err,error,*999)
            CALL MatrixProductTranspose(temp,Fe,kirchoffTensor,err,error,*999)

            !Calculate the Cauchy stress tensor
            cauchyTensor=kirchoffTensor/Je

            IF(DIAGNOSTICS1) THEN
              CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"",err,error,*999)
              CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Stress tensors:",err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Hydrostatic pressure = ",P,err,error,*999)
              CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Second Piola-Kirchoff stress tensor:",err,error,*999)
              CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
                & 3,3,piolaTensor,WRITE_STRING_MATRIX_NAME_AND_INDICES, '("    T','(",I1,",:)','     :",3(X,E13.6))', &
                & '(12X,3(X,E13.6))',err,error,*999)
              CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Cauchy stress tensor:",err,error,*999)
              CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
                & 3,3,cauchyTensor,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("    sigma','(",I1,",:)',' :",3(X,E13.6))', &
                & '(12X,3(X,E13.6))',err,error,*999)
            ENDIF

            !Calculate dPhi/dZ at the gauss point, Phi is the basis function
            !CALL FINITE_ELASTICITY_GAUSS_DFDZ(DEPENDENT_INTERPOLATED_POINT,elementNumber,gauss_idx,numberOfDimensions, &
            !  & numberOfXi,DFDZ,err,error,*999)

            !For membrane theory in 3D space, the final equation is multiplied by thickness. Default to unit thickness if equation set subtype is not membrane
            !!TODO Maybe have the thickness as a component in the equations set field. Yes, as we don't need a materials field for CellML constituative laws.
            THICKNESS = 1.0_DP
            IF(EQUATIONS_SET_SUBTYPE == EQUATIONS_SET_MEMBRANE_SUBTYPE) THEN
              IF(numberOfDimensions == 3) THEN
                IF(ASSOCIATED(MATERIALS_FIELD)) THEN
                  CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
                    & MATERIALS_INTERPOLATED_POINT,err,error,*999)
                  THICKNESS = MATERIALS_INTERPOLATED_POINT%VALUES(MATERIALS_INTERPOLATED_POINT%INTERPOLATION_PARAMETERS% &
                    & FIELD_VARIABLE%NUMBER_OF_COMPONENTS,1)
                ENDIF
              ENDIF
            ENDIF

            !!Now add up the residual terms
            !element_dof_idx=0
            !DO component_idx=1,numberOfDimensions
            !  DEPENDENT_COMPONENT_INTERPOLATION_TYPE=DEPENDENT_FIELD%VARIABLES(var1)%COMPONENTS(component_idx)%INTERPOLATION_TYPE
            !  IF(DEPENDENT_COMPONENT_INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN !node based
            !    DEPENDENT_BASIS=>DEPENDENT_FIELD%VARIABLES(var1)%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY% &
            !      & ELEMENTS%ELEMENTS(elementNumber)%BASIS
            !    NUMBER_OF_FIELD_COMPONENT_INTERPOLATION_PARAMETERS=DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
            !    DO parameter_idx=1,NUMBER_OF_FIELD_COMPONENT_INTERPOLATION_PARAMETERS
            !      element_dof_idx=element_dof_idx+1
            !      DO component_idx2=1,numberOfDimensions
            !        nonlinearMatrices%elementResidual%vector(element_dof_idx)= &
            !          & nonlinearMatrices%elementResidual%vector(element_dof_idx)+ &
            !          & GAUSS_WEIGHT*Jzxi*THICKNESS*cauchyTensor(component_idx,component_idx2)* &
            !          & DFDZ(parameter_idx,component_idx2,component_idx)
            !      ENDDO ! component_idx2 (inner component index)
            !    ENDDO ! parameter_idx (residual vector loop)
            !  ELSEIF(DEPENDENT_COMPONENT_INTERPOLATION_TYPE==FIELD_ELEMENT_BASED_INTERPOLATION) THEN
            !    !Will probably never be used
            !    CALL FlagError("Finite elasticity with element based interpolation is not implemented.",err,error,*999)
            !  ENDIF
            !ENDDO ! component_idx

            !Loop over geometric dependent basis functions.
            DO nh=1,numberOfDimensions
              MESH_COMPONENT_NUMBER=FIELD_VARIABLE%COMPONENTS(nh)%MESH_COMPONENT_NUMBER
              DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%ptr% &
                & TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS
              COMPONENT_QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
              DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                !Loop over derivative directions.
                DO mh=1,numberOfDimensions
                  SUM1=0.0_DP
                  DO mi=1,numberOfXi
                    SUM1=SUM1+DEPENDENT_INTERPOLATED_POINT_METRICS%DXI_DX(mi,mh)* &
                      & COMPONENT_QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(mi),gauss_idx)
                  ENDDO !mi
                  DPHIDZ(mh,ns,nh)=SUM1
                ENDDO !mh
              ENDDO !ns
            ENDDO !nh
            JGW=Jzxi*DEPENDENT_QUADRATURE_SCHEME%GAUSS_WEIGHTS(gauss_idx)
            !Now add up the residual terms
            elementDofIdx=0
            DO mh=1,numberOfDimensions
              MESH_COMPONENT_NUMBER=FIELD_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
              DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%ptr% &
                & TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                elementDofIdx=elementDofIdx+1
                nonlinearMatrices%elementResidual%vector(elementDofIdx)= &
                  & nonlinearMatrices%elementResidual%vector(elementDofIdx)+ &
                  & JGW*DOT_PRODUCT(DPHIDZ(1:numberOfDimensions,ms,mh),cauchyTensor(1:numberOfDimensions,mh))
              ENDDO !ms
            ENDDO !mh
            
            !Hydrostatic pressure component (skip for membrane problems)
            IF (EQUATIONS_SET_SUBTYPE /= EQUATIONS_SET_MEMBRANE_SUBTYPE) THEN
              IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_REFERENCE_STATE_TRANSVERSE_GUCCIONE_SUBTYPE) THEN
                HYDROSTATIC_PRESSURE_COMPONENT=GEOMETRIC_FIELD%VARIABLES(var1)%NUMBER_OF_COMPONENTS
                DEPENDENT_COMPONENT_INTERPOLATION_TYPE=GEOMETRIC_FIELD%VARIABLES(var1)%COMPONENTS( &
                  & HYDROSTATIC_PRESSURE_COMPONENT)%INTERPOLATION_TYPE
              ELSE
                HYDROSTATIC_PRESSURE_COMPONENT=DEPENDENT_FIELD%VARIABLES(var1)%NUMBER_OF_COMPONENTS
                DEPENDENT_COMPONENT_INTERPOLATION_TYPE=DEPENDENT_FIELD%VARIABLES(var1)%COMPONENTS( &
                  & HYDROSTATIC_PRESSURE_COMPONENT)%INTERPOLATION_TYPE
              ENDIF
              IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE) THEN
                TEMPTERM1=GAUSS_WEIGHT*Jxxi*(Jznu-(Jg-DARCY_VOL_INCREASE))
              ELSE
                TEMPTERM1=GAUSS_WEIGHT*Jxxi*(Jznu-Jg)
              ENDIF            
              IF(DEPENDENT_COMPONENT_INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN !node based
                IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_REFERENCE_STATE_TRANSVERSE_GUCCIONE_SUBTYPE) THEN
                  COMPONENT_BASIS=>GEOMETRIC_FIELD%VARIABLES(var1)%COMPONENTS(HYDROSTATIC_PRESSURE_COMPONENT)%DOMAIN% &
                    & TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS
                ELSE
                  COMPONENT_BASIS=>DEPENDENT_FIELD%VARIABLES(var1)%COMPONENTS(HYDROSTATIC_PRESSURE_COMPONENT)%DOMAIN% &
                    & TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS
                ENDIF
                COMPONENT_QUADRATURE_SCHEME=>COMPONENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
                NUMBER_OF_FIELD_COMPONENT_INTERPOLATION_PARAMETERS=COMPONENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                DO parameter_idx=1,NUMBER_OF_FIELD_COMPONENT_INTERPOLATION_PARAMETERS
                  elementDofIdx=elementDofIdx+1
                  IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE) THEN
                    nonlinearMatrices%elementResidual%vector(elementDofIdx)= &
                      & nonlinearMatrices%elementResidual%vector(elementDofIdx)+ &
                      & GAUSS_WEIGHT*Jzxi*COMPONENT_QUADRATURE_SCHEME%GAUSS_BASIS_FNS(parameter_idx,1,gauss_idx)* &
                      & (Je-1.0_DP-DARCY_VOL_INCREASE)
                  ELSE
                    nonlinearMatrices%elementResidual%vector(elementDofIdx)= &
                      & nonlinearMatrices%elementResidual%vector(elementDofIdx)+ &
                      & GAUSS_WEIGHT*Jzxi*COMPONENT_QUADRATURE_SCHEME%GAUSS_BASIS_FNS(parameter_idx,1,gauss_idx)* &
                      & (Je-1.0_DP)
                  ENDIF
                ENDDO
              ELSE IF(DEPENDENT_COMPONENT_INTERPOLATION_TYPE==FIELD_ELEMENT_BASED_INTERPOLATION) THEN !element based
                elementDofIdx=elementDofIdx+1
                nonlinearMatrices%elementResidual%vector(elementDofIdx)= &
                  & nonlinearMatrices%elementResidual%vector(elementDofIdx)+TEMPTERM1
              ENDIF
            ENDIF
          ENDDO !gauss_idx

          !Call surface pressure term here: should only be executed if THIS element has surface pressure on it (direct or incremented)
          IF(DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BOUNDARY_ELEMENT.AND. &
            & TOTAL_NUMBER_OF_SURFACE_PRESSURE_CONDITIONS>0) THEN    ! 
            CALL FiniteElasticity_SurfacePressureResidualEvaluate(EQUATIONS_SET,elementNumber,var1,var2,err,error,*999)
          ENDIF

        ! ---------------------------------------------------------------
        CASE(EQUATIONS_SET_RATE_BASED_SMOOTH_MODEL_SUBTYPE,EQUATIONS_SET_COMPRESSIBLE_RATE_BASED_SMOOTH_MODEL_SUBTYPE, &
          & EQUATIONS_SET_RATE_BASED_GROWTH_MODEL_SUBTYPE,EQUATIONS_SET_COMPRESSIBLE_RATE_BASED_GROWTH_MODEL_SUBTYPE)

          hydrostaticPressureComponent=DEPENDENT_INTERPOLATED_POINT%INTERPOLATION_PARAMETERS%FIELD_VARIABLE%NUMBER_OF_COMPONENTS
          CALL IdentityMatrix(ITens,err,error,*999)
          !Get time step 
          dt=EQUATIONS_SET%deltaTime         
          !Loop over gauss points and add residuals
          DO gaussIdx=1,DEPENDENT_NUMBER_OF_GAUSS_POINTS
            
            IF(diagnostics1) THEN
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Element number = ",elementNumber,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Gauss index  = ",gaussIdx,err,error,*999)
            ENDIF

            gaussWeight=DEPENDENT_QUADRATURE_SCHEME%GAUSS_WEIGHTS(gaussIdx)
            
            !Interpolate dependent, previous dependent, geometric, fibre and materials fields
            CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussIdx, &
              & DEPENDENT_INTERPOLATED_POINT,err,error,*999)
            CALL Field_InterpolatedPointMetricsCalculate(DEPENDENT_BASIS%NUMBER_OF_XI,DEPENDENT_INTERPOLATED_POINT_METRICS, &
              & err,error,*999)
            CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussIdx, &
              & prevDependentInterpPoint,err,error,*999)
            CALL Field_InterpolatedPointMetricsCalculate(DEPENDENT_BASIS%NUMBER_OF_XI,prevDependentInterpPointMetrics, &
              & err,error,*999)
            CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussIdx, &
              & GEOMETRIC_INTERPOLATED_POINT,err,error,*999)
            CALL Field_InterpolatedPointMetricsCalculate(GEOMETRIC_BASIS%NUMBER_OF_XI,GEOMETRIC_INTERPOLATED_POINT_METRICS, &
              & err,error,*999)
            CALL Field_InterpolateGauss(No_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussIdx, &
              & MATERIALS_INTERPOLATED_POINT,err,error,*999)
            IF(ASSOCIATED(FIBRE_FIELD)) THEN
              CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussIdx, &
                & FIBRE_INTERPOLATED_POINT,err,error,*999)
            ENDIF
            
             !Calculate F=dZ/dNu, the deformation gradient tensor at the gauss point
            CALL FiniteElasticity_GaussDeformationGradientTensor(DEPENDENT_INTERPOLATED_POINT_METRICS, &
              & GEOMETRIC_INTERPOLATED_POINT_METRICS,FIBRE_INTERPOLATED_POINT,dZdNu,err,error,*999)

            !Calculate prevF=prevdZ/dNu, the relative deformation gradient tensor at the gauss point
!!TODO: What happens with the fibres here! They will not be in an orthogonal system
            CALL FiniteElasticity_GaussDeformationGradientTensor(prevDependentInterpPointMetrics, &
              & GEOMETRIC_INTERPOLATED_POINT_METRICS,FIBRE_INTERPOLATED_POINT,prevdZdNu,err,error,*999)

            !Get BePrime from the start of the time step
            DO rowIdx=1,numberOfDimensions
              DO columnIdx=rowIdx,numberOfDimensions
                componentIdx=1+TENSOR_TO_VOIGT(rowIdx,columnIdx,numberOfDimensions)
                CALL Field_ParameterSetGetLocalGaussPoint(DEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & gaussIdx,elementNumber,componentIdx,BePrime1(rowIdx,columnIdx),err,error,*999)
              ENDDO !columnIdx
              DO columnIdx=1,rowIdx-1
                BePrime1(rowIdx,columnIdx)=BePrime1(columnIdx,rowIdx)
              ENDDO !columnIdx
            ENDDO !rowIdx
            

            !Calculate Fr = F.prevF^-1
            CALL Invert(prevdZdNu,invPrevdZdNu,prevJF,err,error,*999)
            CALL MatrixProduct(dZdNu,invPrevdZdNu,Fr,err,error,*999)
            !Calculate Jr
            CALL Determinant(Fr,Jr,err,error,*999)
            
            IF(diagnostics1) THEN
              CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"",err,error,*999)
              CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Deformation gradient tensors:",err,error,*999)
              CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Deformation gradient at time t:",err,error,*999)
              CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
                & 3,3,prevdZdNu,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("            F','(",I1,",:)',' :",3(X,E13.6))', &
                & '(16X,3(X,E13.6))',err,error,*999)
              CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Deformation gradient at time t + dt:",err,error,*999)
              CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
                & 3,3,dZdNu,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("            F','(",I1,",:)',' :",3(X,E13.6))', &
                & '(16X,3(X,E13.6))',err,error,*999)
              CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Relative deformation gradient:",err,error,*999)
              CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
                & 3,3,Fr,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("           Fr','(",I1,",:)',' :",3(X,E13.6))', &
                & '(16X,3(X,E13.6))',err,error,*999)
              CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Jacobians:",err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  J at time t       = ",prevJF,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  J at time t + dt  = ",J,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Jr                = ",Jr,err,error,*999)
            ENDIF
            
            JXXi=GEOMETRIC_INTERPOLATED_POINT_METRICS%JACOBIAN            
            JZxi=DEPENDENT_INTERPOLATED_POINT_METRICS%JACOBIAN
            prevJZxi=prevDependentInterpPointMetrics%JACOBIAN            
                        
            SELECT CASE(EQUATIONS_SET_SUBTYPE)
            CASE(EQUATIONS_SET_RATE_BASED_SMOOTH_MODEL_SUBTYPE,EQUATIONS_SET_COMPRESSIBLE_RATE_BASED_SMOOTH_MODEL_SUBTYPE)
              
              !Equation numbers are from M. Hollenstein, M. Jabareen and M.B. Rubin "Modeling a smoth elastic-inelastic
              !transition with a strongly objective numerical integrator needing no iteration". 2013. Comput Mech. 52:649-667.

              mu=MATERIALS_INTERPOLATED_POINT%VALUES(1,1)
              a0=MATERIALS_INTERPOLATED_POINT%VALUES(2,1)
              a1=MATERIALS_INTERPOLATED_POINT%VALUES(3,1)
              b0=MATERIALS_INTERPOLATED_POINT%VALUES(4,1)
              b1=MATERIALS_INTERPOLATED_POINT%VALUES(5,1)
              m=MATERIALS_INTERPOLATED_POINT%VALUES(6,1)
              kappas=MATERIALS_INTERPOLATED_POINT%VALUES(7,1)
              
              !Get the hardening variable, kappa
              CALL Field_ParameterSetGetLocalGaussPoint(DEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & gaussIdx,elementNumber,1,kappa1,err,error,*999)
              
              !props(1)=K
              !props(2)=mu
              !props(3)=a0
              !props(4)=a1
              !props(5)=b0
              !props(6)=b1
              !props(7)=m
              !props(8)=kappas
              
              !statev(1)=prevJF
              !statev(2)=BePrime1(1,1)
              !statev(3)=BePrime1(2,2)
              !statev(4)=BePrime1(3,3)
              !statev(5)=BePrime1(1,2)
              !statev(6)=BePrime1(2,3)
              !statev(7)=BePrime1(1,3)
              !statev(8)=kappa1
              !statev(9)=0.0_DP
              !statev(10)=0.0_DP
              !statev(11)=0.0_DP
              !statev(12)=statev(2)+statev(3)+statev(4)
              !statev(13)=0.0_DP
              
              !CALL umat(6,13,8,2,1,prevdZdNu,dZdNu,dt,cauchyTensor,statev,ddsdde,props)
                            
              !Eqn (4.6)
              J=Jr*prevJF
              !Eqn (4.4)
              CALL Unimodular(Fr,FrPrime,err,error,*999)
              !Calculate Elastic Trial. Eqn (4.7)
              CALL MatrixProduct(FrPrime,BePrime1,tempTensor,err,error,*999)
              CALL MatrixProductTranspose(tempTensor,FrPrime,BePrimeStar,err,error,*999)
              CALL DecomposeSphericalDeviatoric(BePrimeStar,hydBePrimeStar,BePrimePrimeStar,err,error,*999)
              !Calculate gePrimePrimeStar. Eqn (4.10)
              gePrimePrimeStar=BePrimePrimeStar/2.0_DP
              !Calculate Equivalent Strain, gammaEStar. Use Eqn (4.10) in Eqn (4.17)
              CALL DoubleDotProduct(BePrimePrimeStar,BePrimePrimeStar,gammaEStar,err,error,*999)
              gammaEStar=SQRT(3.0_DP/8.0_DP*gammaEStar)
              !Calculate Total Effective Distortional Strain, Dbar. Eqn (4.13)
              CALL MatrixProductTranspose(Fr,Fr,Br,err,error,*999)
              Dbar=(Br-ITens)/(2.0_DP*dt)
              !Calculate delta epsilon. Eqn (4.14). Note typo in paper. It should be devDbar dot devDbar
              CALL DecomposeSphericalDeviatoric(Dbar,hydDbar,devDbar,err,error,*999)
              CALL DoubleDotProduct(devDbar,devDbar,deltaEps,err,error,*999)
              deltaEps=dt*SQRT(2.0_DP/3.0_DP*deltaEps)
              !Calculate gamma. Eqn (4.18)
              dtGamma0=dt*a0+b0*deltaEps
              dtGamma1=dt*a1+b1*deltaEps
              !Eqn (4.21)
              c0=dtGamma1*(gammaEStar-(1.0_DP+dtGamma0)*kappa1)
              IF(c0<=ZERO_TOLERANCE) THEN
                !Eqn (4.19)
                gamma=0.0_DP
              ELSE
                !Eqn (4.21)
                c1=gammaEStar+dtGamma1*(kappa1-m*(gammaEStar-(1.0_DP+dtGamma0)*kappas))
                IF(ABS(m)<ZERO_TOLERANCE) THEN
                  gamma=c0/c1
                ELSE
                  !Eqn (4.21)
                  c2=m*(gammaEStar+dtGamma1*kappas)
                  !Eqn (4.22)
                  gamma=(-c1+SQRT(c1*c1+4.0_DP*c0*c2))/(2.0_DP*c2)
                ENDIF
              ENDIF
              !Calculate dtGamma. Eqn (4.15)
              dtGamma=dtGamma0+gamma
              !Calculate Hardening. Eqn (4.16)
              kappa=(kappa1+m*kappas*gamma)/(1.0_DP+m*gamma)
              !Calculate devBePrimePrime. Eqn (4.11)
              devBePrimePrime=2.0_DP/(1.0_DP+dtGamma)*gePrimePrimeStar
              !Calculate gePrimePrime (deviatoric part of gePrime). Eqn (3.1) and (4.11)
              gePrimePrime=devBePrimePrime/2.0_DP
              !Calculate alpha
              CALL Determinant(devBePrimePrime,detdevBePrimePrime,err,error,*999)
              CALL DoubleDotProduct(devBePrimePrime,devBePrimePrime,factor1,err,error,*999)
              factor1=2.0_DP/3.0_DP*factor1
              IF(ABS(factor1)<ZERO_TOLERANCE) THEN
                alpha1=3.0_DP
              ELSE
                factor2=(4.0_DP*(1.0_DP-detdevBePrimePrime))/(factor1**1.5_DP)
                IF(factor2>=1.0_DP) THEN
                  alpha1=3.0_DP*SQRT(factor1)*COSH(ACOSH(factor2)/3.0_DP)
                ELSE
                  alpha1=3.0_DP*SQRT(factor1)*COS(ACOS(factor2)/3.0_DP)
                ENDIF
              ENDIF
              !Calculate BePrime. Eqn (4.12)
              BePrime=alpha1/3.0_DP*ITens+devBePrimePrime
              !Calculate deviatoric part of Cauchy stress. Eqn (6.2)
              devT=2.0_DP*mu*gePrimePrime/J
              SELECT CASE(EQUATIONS_SET_SUBTYPE)
              CASE(EQUATIONS_SET_RATE_BASED_SMOOTH_MODEL_SUBTYPE)
                P=DEPENDENT_INTERPOLATED_POINT%VALUES(hydrostaticPressureComponent,1)
              CASE(EQUATIONS_SET_COMPRESSIBLE_RATE_BASED_SMOOTH_MODEL_SUBTYPE)
                P=-K*(J-1.0_DP)
              END SELECT
              !Calculate Cauchy stress. Eqn (2.9)
              cauchyTensor=devT-P*ITens
              
              IF(diagnostics1) THEN
                CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Strain tensors:",err,error,*999)
                CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Be prime at t:",err,error,*999)
                CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
                  & 3,3,BePrime1,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("    Be prime','(",I1,",:)',' :",3(X,E13.6))', &
                  & '(16X,3(X,E13.6))',err,error,*999)
                CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Be prime at t+dt:",err,error,*999)
                CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
                  & 3,3,BePrime,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("    Be prime','(",I1,",:)',' :",3(X,E13.6))', &
                  & '(16X,3(X,E13.6))',err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  gamma = ",gamma,err,error,*999)
                CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Strain hardening:",err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  kappa at time t      = ",kappa1,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  kappa at time t + dt = ",kappa,err,error,*999)
                CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Stress tensors:",err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Hydrostatic pressure = ",P,err,error,*999)
                CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Deviatoric Cauchy stress tensor:",err,error,*999)
                CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
                  & 3,3,devT,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("    dev sigma','(",I1,",:)',' :",3(X,E13.6))', &
                  & '(16X,3(X,E13.6))',err,error,*999)
                CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Cauchy stress tensor:",err,error,*999)
                CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
                  & 3,3,cauchyTensor,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("        sigma','(",I1,",:)',' :",3(X,E13.6))', &
                  & '(16X,3(X,E13.6))',err,error,*999)
              ENDIF

              !Store the current value of the hardening variable, kappa, in the next values parameter set
              CALL Field_ParameterSetUpdateLocalGaussPoint(DEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE,FIELD_NEXT_VALUES_SET_TYPE, &
                & gaussIdx,elementNumber,1,kappa,err,error,*999)            
                           
            CASE(EQUATIONS_SET_COMPRESSIBLE_RATE_BASED_GROWTH_MODEL_SUBTYPE,EQUATIONS_SET_RATE_BASED_GROWTH_MODEL_SUBTYPE)
              
              mu0=MATERIALS_INTERPOLATED_POINT%VALUES(1,1)
              q=MATERIALS_INTERPOLATED_POINT%VALUES(2,1)
              kt=MATERIALS_INTERPOLATED_POINT%VALUES(3,1)
              kc=MATERIALS_INTERPOLATED_POINT%VALUES(4,1)
              k1=MATERIALS_INTERPOLATED_POINT%VALUES(5,1)
              k2=MATERIALS_INTERPOLATED_POINT%VALUES(6,1)
              k3=MATERIALS_INTERPOLATED_POINT%VALUES(7,1)
              k12=MATERIALS_INTERPOLATED_POINT%VALUES(8,1)
              ka=MATERIALS_INTERPOLATED_POINT%VALUES(9,1)
              Gamma=MATERIALS_INTERPOLATED_POINT%VALUES(10,1)
              GammaM=MATERIALS_INTERPOLATED_POINT%VALUES(11,1)
              Jh=MATERIALS_INTERPOLATED_POINT%VALUES(12,1)
              H11=MATERIALS_INTERPOLATED_POINT%VALUES(13,1)
              H22=MATERIALS_INTERPOLATED_POINT%VALUES(14,1)
              H33=MATERIALS_INTERPOLATED_POINT%VALUES(15,1)
              H12=MATERIALS_INTERPOLATED_POINT%VALUES(16,1)
              H13=MATERIALS_INTERPOLATED_POINT%VALUES(17,1)
              H23=MATERIALS_INTERPOLATED_POINT%VALUES(18,1)
              
              !Get the current value of Je at the start of the time step
              CALL Field_ParameterSetGetLocalGaussPoint(DEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & gaussIdx,elementNumber,1,Je1,err,error,*999)
              
              !Get the structural tensor, S1, from the start of the time step
              DO rowIdx=1,numberOfDimensions
                DO columnIdx=1,numberOfDimensions
                  componentIdx=rowIdx+(columnIdx-1)*numberOfDimensions
                  CALL Field_ParameterSetGetLocalGaussPoint(DEPENDENT_FIELD,FIELD_U2_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                    & gaussIdx,elementNumber,componentIdx,S1(rowIdx,columnIdx),err,error,*999)
                ENDDO !columnIdx
              ENDDO !rowIdx

              !CALL clooping(1.0E-10_DP,1.0E-10_DP,100.0_DP,mu0,kc,kt,1.0_DP,k2,k3,k12,ka,q,H11,H22,H12,Jh,Gamma,Gammam,dt,Fr, &
              !  & Jr,Je1,S1(:,1),S1(:,2),S1(:,3),devH,HH,malpha1,B,Bepr,lame1,lame2,lame3,lamea,Ee11,Ee22,Ee33,Ee12,Eea,QQ,T0,T1, &
              !  & T2,T3,T4,T5,TT)

              !Calculate the elastic dilatation. Eqn (6.5) and (6.8)
              Je2 = Jr*Je1
              Je = EXP(LOG(Je2)+dt*GammaM*LOG(Jh))/(1.0_DP+dt*GammaM)
              !Update state variables
              !Update the fibre vectors s1 and s3 to give S at t2 i.e., S2 Eqn (6.4)
              CALL MatrixVectorProduct(Fr,S1(:,1),S2(:,1),err,error,*999)
              CALL Normalise(S2(:,1),S2(:,1),err,error,*999)
              CALL InvertTranspose(Fr,invFrT,detInvFrT,err,error,*999)
              CALL MatrixVectorProduct(invFrT,S1(:,3),S2(:,3),err,error,*999)
              CALL Normalise(S2(:,3),S2(:,3),err,error,*999)
              !Compute s2 = s3 x s1. Eqn (3.1)
              CALL CrossProduct(S2(:,3),S2(:,1),S2(:,2),err,error,*999)
              !Compute the structural tensor(s). Eqn (3.5)
!!TODO: Change this to use a matrix matrix tensor product?/ Indices will change???
              DO rowIdx=1,numberOfDimensions
                DO columnIdx=1,numberOfDimensions
                  CALL TensorProduct(S2(:,rowIdx),S2(:,columnIdx),S(rowIdx,columnIdx,:,:),err,error,*999)
                ENDDO !columnIdx
              ENDDO !rowIdx
              !Compute the homostatis growth tensor. Eqn (7.4) and (???)
              HPrimePrime=H11*S(1,1,:,:)+H22*S(2,2,:,:)-(H11+H22)*S(3,3,:,:)+H12*(S(1,2,:,:)+S(2,1,:,:))
              H=ITens+HPrimePrime
              !Compute alphas via non-linear Newton solve              
              !Calculate BePrimeStar
              CALL MatrixProduct(FrPrime,BePrime1,tempTensor,err,error,*999)
              CALL MatrixProductTranspose(tempTensor,FrPrime,BePrimeStar,err,error,*999)
              !Calculate BePrimePrimeStar (deviatoric part of BePrimeStar)
              CALL DecomposeSphericalDeviatoric(BePrimeStar,hydBePrimeStar,BePrimePrimeStar,err,error,*999)
              !Initialise alpha
              CALL Unimodular(BePrime1,uniBePrime1,err,error,*999)
              CALL Trace(uniBePrime1,alpha(1),err,error,*999)
              alpha(2)=1.0_DP              
              numIterations=0
              deltaAlphaNorm=1.0_DP
              deltaAlphaNormTolerance=10.0E-10
              hStep=10.0E-10
              maxNumIterations=100
              DO WHILE((deltaAlphaNorm>=deltaAlphaNormTolerance).AND.(numIterations<maxNumIterations))
                !Increment iteration count
                numIterations=numIterations+1
                !Calculate deviatoric part of BePrime2, BePrimePrime
                BePrimePrime=(BePrimePrimeStar+alpha(2)*dt*Gamma*HPrimePrime)/(1.0_DP+dt*Gamma)
                !Calculate elastic distortion, Beprime
                BePrime=alpha(1)*ITens/3.0_DP+BePrimePrime
                !Calculate resid(1)
                !Calculate inverse of BePrime
                CALL Invert(BePrime,invBePrime,detBePrime,err,error,*999)
                CALL DoubleDotProduct(invBePrime,H,BeDDotH,err,error,*999)
                resid(1)=alpha(2)-3.0_DP/BeDDotH
                !Calculate resid(2)
                CALL DoubleDotProduct(BePrimePrime,BePrimePrime,BeDDotBe,err,error,*999)
                CALL Determinant(BePrimePrime,detBePrimePrime,err,error,*999)
                resid(2)=(alpha(1)/3.0_DP)**3.0_DP-(BeDDotBe*alpha(1)/6.0_DP)-(1.0_DP-detBePrimePrime)
                perturbedAlpha=alpha
                DO columnIdx=1,2
                  hStepAlpha=MAX(hStep,ABS(hstep*perturbedAlpha(columnIdx)))
                  perturbedAlpha(columnIdx)=perturbedAlpha(columnIdx)+hstepAlpha
                  !Calculate deviatoric part of BePrime2, BePrimePrime
                  BePrimePrime=(BePrimePrimeStar+perturbedAlpha(2)*dt*Gamma*HPrimePrime)/(1.0_DP+dt*Gamma)
                  !Calculate elastic distortion, Beprime
                  BePrime=perturbedAlpha(1)*ITens/3.0_DP+BePrimePrime
                  !Calculate resid(1)
                  !Calculate inverse of BePrime
                  CALL Invert(BePrime,invBePrime,detBePrime,err,error,*999)
                  CALL DoubleDotProduct(invBePrime,H,BeDDotH,err,error,*999)
                  perturbedResid(1)=perturbedAlpha(2)-3.0_DP/BeDDotH
                  !Calculate resid(2)
                  CALL DoubleDotProduct(BePrimePrime,BePrimePrime,BeDDotBe,err,error,*999)
                  CALL Determinant(BePrimePrime,detBePrimePrime,err,error,*999)
                  perturbedResid(2)=(perturbedAlpha(1)/3.0_DP)**3.0_DP-(BeDDotBe*perturbedAlpha(1)/6.0_DP)-(1.0_DP-detBePrimePrime)
                  DO rowIdx=1,2
                    Jacobian(rowIdx,columnIdx)=(perturbedResid(rowIdx)-resid(rowIdx))/hStepAlpha
                  ENDDO !rowIdx
                  perturbedAlpha(columnIdx)=alpha(columnIdx)
                ENDDO !columnIdx
                CALL Invert(Jacobian,invJacobian,detJacobian,err,error,*999)
                CALL MatrixVectorProduct(invJacobian,resid,deltaAlpha,err,error,*999)
                CALL L2Norm(deltaAlpha,deltaAlphaNorm,err,error,*999)
                alpha=alpha-deltaAlpha
              ENDDO
              !Calculate deviatoric part of BePrime2, BePrimePrime
              BePrimePrime=(BePrimePrimeStar+alpha(2)*dt*Gamma*HPrimePrime)/(1.0_DP+dt*Gamma)
              !Calculate elastic distortion, Beprime
              BePrime=alpha(1)*ITens/3.0_DP+BePrimePrime
              !Calculate strain information
              !Calculate Be
              Be=(Je**(2.0_DP/3.0_DP))*BePrime
              !Invert Be
              CALL Invert(Be,invBe,detBe,err,error,*999)
              !Calculate lame and Ee
              Ee=0.0_DP
              DO rowIdx=1,numberOfDimensions
                CALL DoubleDotProduct(invBe,S(rowIdx,rowIdx,:,:),lame(rowIdx),err,error,*999)
                lame(rowIdx)=1.0_DP/SQRT(lame(rowIdx))
                Ee(rowIdx,rowIdx)=(lame(rowIdx)*lame(rowIdx)-1.0_DP)/2.0_DP
              ENDDO !rowIdx
              CALL DoubleDotProduct(invBe,S(1,2,:,:),Ee(1,2),err,error,*999)
              Ee(1,2)=-1.0_DP*lame(1)*lame(2)**Ee(1,2)/2.0_DP
              CALL DoubleDotProduct(Be,S(3,3,:,:),lame(4),err,error,*999)
              lame(4)=Je/SQRT(lame(4))
              Ee(4,4)=lame(4)-1.0_DP
              !Calculate stress information
              ktJe=MacaulayBracket(Je-1.0_DP)
              kcJe=MacaulayBracket(1.0_DP-Je)
              macEe11=MacaulayBracket(Ee(1,1))
              macEe22=MacaulayBracket(Ee(2,2))
              macEe33=MacaulayBracket(Ee(3,3))
              CALL MatrixProduct(invBe,S(3,3,:,:),invBeS33,err,error,*999)
              CALL MatrixProduct(S(3,3,:,:),invBe,S33invBe,err,error,*999)
              
              QQ=0.50_DP*kt*(ktJe**2.0_DP)+0.50_DP*kc*(kcJe**2.0_DP) &
                &  +0.50_DP*(alpha(1)-3.0_DP)+0.50_DP*k1*(macEe11**2.0_DP) &
                &  +0.50_DP*k2*(macEe22**2.0_DP)+0.50_DP*k3*(macEe33**2.0_DP) &
                &  +2.0_DP*k12*(Ee(1,2)**2.0_DP)+0.50_DP*ka*(Ee(4,4)**2.0d0)
              
              mu=mu0*EXP(q*QQ)
              T0=mu*(kt*ktJe-kc*kcJe)*iTens+(mu/Je)*BePrimePrime
              T1=((k1*mu*(lame(1)**2.0_DP)*macEe11)/Je)*S(1,1,:,:)
              T2=((k2*mu*(lame(2)**2.0_DP)*macEe22)/Je)*(S(2,2,:,:)-(2.0_DP*lame(2)/lame(1))*Ee(1,2)*(S(1,2,:,:)+S(2,1,:,:)))
              T3=((k3*mu*(lame(3)**2.0_DP)*macEe33)/Je)*((lame(3)**2.0_DP)*(invBeS33+S33invBe)-S(3,3,:,:))
              T4=2.0_DP*mu*k12*Ee(1,2)*(lame(2)/lame(1)/Je)*(1.0d0-4.0d0*(Ee(1,2)**2.0_DP))*(S(1,2,:,:)+S(2,1,:,:))
              T5=(mu*ka*Ee(4,4)*lame(4)/Je)*(S(1,1,:,:)+S(2,2,:,:))

              cauchyTensor=T0+T1+T2+T3+T4+T5
              
              IF(diagnostics1) THEN
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Je at time t      = ",Je1,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Je at time t + dt = ",Je2,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Je                = ",Je,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Jh                = ",Jh,err,error,*999)
                CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Structural fibre tensors:",err,error,*999)
                CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  S at t:",err,error,*999)
                CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
                  & 3,3,S1,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("           S','(",I1,",:)',' :",3(X,E13.6))', &
                  & '(16X,3(X,E13.6))',err,error,*999)
                CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  S at t + dt:",err,error,*999)
                CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
                  & 3,3,S2,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("           S','(",I1,",:)',' :",3(X,E13.6))', &
                  & '(16X,3(X,E13.6))',err,error,*999)
                CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Growth tensors:",err,error,*999)
                CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  H prime prime:",err,error,*999)
                CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
                  & 3,3,HPrimePrime,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("           H','(",I1,",:)',' :",3(X,E13.6))', &
                  & '(16X,3(X,E13.6))',err,error,*999)
                CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  H:",err,error,*999)
                CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
                  & 3,3,H,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("           H','(",I1,",:)',' :",3(X,E13.6))', &
                  & '(16X,3(X,E13.6))',err,error,*999)
                CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Strain tensors:",err,error,*999)
                CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Be prime at t:",err,error,*999)
                CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
                  & 3,3,BePrime1,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("    Be prime','(",I1,",:)',' :",3(X,E13.6))', &
                  & '(16X,3(X,E13.6))',err,error,*999)
                CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Be prime at t + dt:",err,error,*999)
                CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
                  & 3,3,BePrime,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("    Be prime','(",I1,",:)',' :",3(X,E13.6))', &
                  & '(16X,3(X,E13.6))',err,error,*999)
                CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,2,2,2,alpha,'("",2(X,E13.6))','(2(X,E13.6))',err,error,*999)
                CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Ee at t + dt:",err,error,*999)
                CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,4,1,1,4, &
                  & 4,4,Ee,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("          Ee','(",I1,",:)',' :",4(X,E13.6))', &
                  & '(16X,4(X,E13.6))',err,error,*999)
                CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,4,4,4,lame,'("",4(X,E13.6))','4(X,E13.6))',err,error,*999)
                CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Stress tensors:",err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  QQ = ",QQ,err,error,*999)
                CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  T1-5:",err,error,*999)
                CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3,3,3,T1,WRITE_STRING_MATRIX_NAME_AND_INDICES, &
                  & '("            T1','(",I1,",:)',' :",3(X,E13.6))','(16X,3(X,E13.6))',err,error,*999)
                CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3,3,3,T2,WRITE_STRING_MATRIX_NAME_AND_INDICES, &
                  & '("            T2','(",I1,",:)',' :",3(X,E13.6))','(16X,3(X,E13.6))',err,error,*999)
                CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3,3,3,T3,WRITE_STRING_MATRIX_NAME_AND_INDICES, &
                  & '("            T3','(",I1,",:)',' :",3(X,E13.6))','(16X,3(X,E13.6))',err,error,*999)
                CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3,3,3,T4,WRITE_STRING_MATRIX_NAME_AND_INDICES, &
                  & '("            T4','(",I1,",:)',' :",3(X,E13.6))','(16X,3(X,E13.6))',err,error,*999)
                CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3,3,3,T5,WRITE_STRING_MATRIX_NAME_AND_INDICES, &
                  & '("            T5','(",I1,",:)',' :",3(X,E13.6))','(16X,3(X,E13.6))',err,error,*999)
                CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Cauchy stress tensor:",err,error,*999)
                CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
                  & 3,3,cauchyTensor,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("        sigma','(",I1,",:)',' :",3(X,E13.6))', &
                  & '(16X,3(X,E13.6))',err,error,*999)
              ENDIF
              
              !Store the current value of Je in the next values parameter set
              CALL Field_ParameterSetUpdateLocalGaussPoint(DEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE,FIELD_NEXT_VALUES_SET_TYPE, &
                & gaussIdx,elementNumber,1,Je,err,error,*999)              
            
              !Store the current values of the structural tensor, S2, in the next values parameter set
              DO rowIdx=1,numberOfDimensions
                DO columnIdx=1,numberOfDimensions
                  componentIdx=rowIdx+(columnIdx-1)*numberOfDimensions
                  CALL Field_ParameterSetUpdateLocalGaussPoint(DEPENDENT_FIELD,FIELD_U2_VARIABLE_TYPE,FIELD_NEXT_VALUES_SET_TYPE, &
                    & gaussIdx,elementNumber,componentIdx,S2(rowIdx,columnIdx),err,error,*999)
                ENDDO !columnIdx
              ENDDO !rowIdx
              
           CASE DEFAULT
              localError="The equations set subtype of "//TRIM(NumberToVString(EQUATIONS_SET_SUBTYPE,"*",err,error))// &
                & " is invalid."
              CALL FlagError(localError,err,error,*999)              
            END SELECT
            
            !Store BePrime in the next values parameter set
            DO rowIdx=1,numberOfDimensions
              DO columnIdx=rowIdx,numberOfDimensions
                componentIdx=1+TENSOR_TO_VOIGT(rowIdx,columnIdx,numberOfDimensions)
                CALL Field_ParameterSetUpdateLocalGaussPoint(DEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE,FIELD_NEXT_VALUES_SET_TYPE, &
                  & gaussIdx,elementNumber,componentIdx,BePrime(rowIdx,columnIdx),err,error,*999)            
              ENDDO !columnIdx
            ENDDO !rowIdx

            !Loop over dependent columns directions.
            DO columnComponentIdx=1,numberOfDimensions
              meshComponentNumber=FIELD_VARIABLE%COMPONENTS(columnComponentIdx)%MESH_COMPONENT_NUMBER
              dependentBasis=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(meshComponentNumber)%ptr% &
                & TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS
              quadratureScheme=>dependentBasis%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
              !Loop over dependent columns element parameters.
              DO columnElementParameterIdx=1,dependentBasis%NUMBER_OF_ELEMENT_PARAMETERS
                !Loop over dependent columns directions.
                DO rowComponentIdx=1,numberOfDimensions
                  sum1=0.0_DP
                  DO xiIdx=1,numberOfXi
                    sum1=sum1+DEPENDENT_INTERPOLATED_POINT_METRICS%DXI_DX(xiIdx,rowComponentIdx)* &
                      & quadratureScheme%GAUSS_BASIS_FNS(columnElementParameterIdx, &
                      PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xiIdx),gaussIdx)
                  ENDDO !xiIdx
                  dPhidZ(rowComponentIdx,columnElementParameterIdx,columnComponentIdx)=sum1
                ENDDO !rowComponentIdx
              ENDDO !columnElementParameterIdx
            ENDDO !columnComponentIdx
            JGw=JZxi*gaussWeight
            !Now add up the residual terms
            rowElementDofIdx=0
            DO rowComponentIdx=1,numberOfDimensions
              meshComponentNumber=FIELD_VARIABLE%COMPONENTS(rowComponentIdx)%MESH_COMPONENT_NUMBER
              dependentBasis=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(meshComponentNumber)%ptr% &
                & TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS
              DO rowElementParameterIdx=1,dependentBasis%NUMBER_OF_ELEMENT_PARAMETERS
                rowElementDofIdx=rowElementDofIdx+1
                nonlinearMatrices%elementResidual%vector(rowElementDofIdx)=nonlinearMatrices%elementResidual% &
                  & vector(rowElementDofIdx)+JGw*DOT_PRODUCT(dPhiDz(1:numberOfDimensions,rowElementParameterIdx,rowComponentIdx), &
                  & cauchyTensor(1:numberOfDimensions,rowComponentIdx))
              ENDDO !rowElementParameterIdx
            ENDDO !rowComponentIdx
            
            SELECT CASE(EQUATIONS_SET_SUBTYPE)
            CASE(EQUATIONS_SET_RATE_BASED_SMOOTH_MODEL_SUBTYPE, &
              & EQUATIONS_SET_RATE_BASED_GROWTH_MODEL_SUBTYPE)
              P=DEPENDENT_INTERPOLATED_POINT%VALUES(hydrostaticPressureComponent,1)
              !Hydrostatic pressure component
              hydrostaticPressureComponent=DEPENDENT_FIELD%VARIABLES(var1)%NUMBER_OF_COMPONENTS
              dependentComponentInterpolationType=DEPENDENT_FIELD%VARIABLES(var1)%COMPONENTS(hydrostaticPressureComponent)% &
                & INTERPOLATION_TYPE
              tempTerm1=gaussWeight*(J-1.0_DP)
              SELECT CASE(dependentComponentInterpolationType)
              CASE(FIELD_CONSTANT_INTERPOLATION)
                CALL FlagError("Not implemented.",err,error,*999)
              CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                rowElementDofIdx=rowElementDofIdx+1
                nonlinearMatrices%elementResidual%vector(rowElementDofIdx)= &
                  & nonlinearMatrices%elementResidual%vector(rowElementDofIdx)+tempTerm1
              CASE(FIELD_NODE_BASED_INTERPOLATION)
                dependentBasis=>DEPENDENT_FIELD%VARIABLES(var1)%COMPONENTS(hydrostaticPressureComponent)%DOMAIN% &
                  & TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%basis
                quadratureScheme=>dependentBasis%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
                DO rowElementParameterIdx=1,dependentBasis%NUMBER_OF_ELEMENT_PARAMETERS
                  mhs=mhs+1 
                  nonlinearMatrices%elementResidual%vector(rowElementParameterIdx)= &
                    & nonlinearMatrices%elementResidual%vector(rowElementParameterIdx)+ &
                    & quadratureScheme%GAUSS_BASIS_FNS(rowElementParameterIdx,1,gaussIdx)*tempTerm1
                ENDDO
              CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                CALL FlagError("Not implemented.",err,error,*999)
              CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                CALL FlagError("Not implemented.",err,error,*999)
              CASE(FIELD_DATA_POINT_BASED_INTERPOLATION)
                CALL FlagError("Not implemented.",err,error,*999)
              CASE DEFAULT
                localError="The dependent field interpolation type of "// &
                  & TRIM(NumberToVString(dependentComponentInterpolationType,"*",err,error))//" for field component "// &
                  & TRIM(NumberToVString(hydrostaticPressureComponent,"*",err,error))//" is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE DEFAULT
              !Do nothing
            END SELECT
            
          ENDDO !gaussIdx

          !Call surface pressure term here: should only be executed if THIS element has surface pressure on it (direct or
          !incremented)
          IF(DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BOUNDARY_ELEMENT.AND. &
            & TOTAL_NUMBER_OF_SURFACE_PRESSURE_CONDITIONS>0) THEN    ! 
            CALL FiniteElasticity_SurfacePressureResidualEvaluate(EQUATIONS_SET,elementNumber,var1,var2,err,error,*999)
          ENDIF
          
        CASE(EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
          !keep the multi-compartment case separate for the time being until the formulation has been finalised, then perhaps
          !integrate within the single compartment case
          !Loop over gauss points and add residuals
          EQUATIONS_SET_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
          CALL Field_ParameterSetDataGet(EQUATIONS_SET_FIELD,FIELD_U_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,EQUATIONS_SET_FIELD_DATA,err,error,*999)

          Ncompartments  = EQUATIONS_SET_FIELD_DATA(2)

          DO gauss_idx=1,DEPENDENT_NUMBER_OF_GAUSS_POINTS
            GAUSS_WEIGHT=DEPENDENT_QUADRATURE_SCHEME%GAUSS_WEIGHTS(gauss_idx)
              !Interpolate dependent, geometric, fibre and materials fields
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & DEPENDENT_INTERPOLATED_POINT,err,error,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(DEPENDENT_BASIS%NUMBER_OF_XI,DEPENDENT_INTERPOLATED_POINT_METRICS, &
              & err,error,*999)
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & GEOMETRIC_INTERPOLATED_POINT,err,error,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,GEOMETRIC_INTERPOLATED_POINT_METRICS, &
              & err,error,*999)
             IF(ASSOCIATED(FIBRE_FIELD)) THEN
              CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
                & FIBRE_INTERPOLATED_POINT,err,error,*999)
            END IF
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & MATERIALS_INTERPOLATED_POINT,err,error,*999)

            IF(DIAGNOSTICS1) THEN
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  elementNumber = ",elementNumber,err,error,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  gauss_idx = ",gauss_idx,err,error,*999)
            ENDIF
            IF(DIAGNOSTICS1) THEN
              CALL WRITE_STRING_MATRIX(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
                & 3,3,piolaTensor,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("     Piola Tensor','(",I1,",:)',' :",3(X,E13.6))', &
                & '(17X,3(X,E13.6))',err,error,*999)

              CALL WRITE_STRING_MATRIX(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
                & 3,3,cauchyTensor,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("    Cauchy Tensor','(",I1,",:)',' :",3(X,E13.6))', &
                & '(17X,3(X,E13.6))',err,error,*999)
            ENDIF

              !Parameters settings for coupled elasticity Darcy INRIA model:
              CALL GET_DARCY_FINITE_ELASTICITY_PARAMETERS(DARCY_RHO_0_F,Mfact,bfact,p0fact,err,error,*999)

              DARCY_MASS_INCREASE = 0.0_DP
              DO imatrix=1,Ncompartments
                DARCY_FIELD_VAR_TYPE=FIELD_V_VARIABLE_TYPE+FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(imatrix-1)
                DARCY_DEPENDENT_INTERPOLATION_PARAMETERS=>&
                  & equations%interpolation%dependentInterpParameters(DARCY_FIELD_VAR_TYPE)%ptr

                CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,elementNumber, &
                  & DARCY_DEPENDENT_INTERPOLATION_PARAMETERS,err,error,*999)

                DARCY_DEPENDENT_INTERPOLATED_POINT=>equations%interpolation%dependentInterpPoint(DARCY_FIELD_VAR_TYPE)%ptr
                CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
                  & DARCY_DEPENDENT_INTERPOLATED_POINT,err,error,*999)

                DARCY_MASS_INCREASE = DARCY_MASS_INCREASE + DARCY_DEPENDENT_INTERPOLATED_POINT%VALUES(4,NO_PART_DERIV)
              ENDDO

              DARCY_VOL_INCREASE = DARCY_MASS_INCREASE / DARCY_RHO_0_F

            !Calculate F=dZ/dNU, the deformation gradient tensor at the gauss point
            CALL FiniteElasticity_GaussDeformationGradientTensor(DEPENDENT_INTERPOLATED_POINT_METRICS, &
              & GEOMETRIC_INTERPOLATED_POINT_METRICS,FIBRE_INTERPOLATED_POINT,DZDNU,ERR,ERROR,*999)

            Jxxi=GEOMETRIC_INTERPOLATED_POINT_METRICS%JACOBIAN

            !Calculate Sigma=1/Jznu.FTF', the Cauchy stress tensor at the gauss point
            CALL FINITE_ELASTICITY_GAUSS_CAUCHY_TENSOR(EQUATIONS_SET,DEPENDENT_INTERPOLATED_POINT, &
              & MATERIALS_INTERPOLATED_POINT,GEOMETRIC_INTERPOLATED_POINT,DARCY_DEPENDENT_INTERPOLATED_POINT, &
              & INDEPENDENT_INTERPOLATED_POINT,cauchyTensor,Jznu,DZDNU,elementNumber,gauss_idx,err,error,*999)

            !Calculate dPhi/dZ at the gauss point, Phi is the basis function
            CALL FINITE_ELASTICITY_GAUSS_DFDZ(DEPENDENT_INTERPOLATED_POINT,elementNumber,gauss_idx,numberOfDimensions, &
              & numberOfXi,DFDZ,err,error,*999)

            !For membrane theory in 3D space, the final equation is multiplied by thickness. Default to unit thickness if equation set subtype is not membrane
            THICKNESS = 1.0_DP
            IF(EQUATIONS_SET_SUBTYPE == EQUATIONS_SET_MEMBRANE_SUBTYPE) THEN
              IF(numberOfDimensions == 3) THEN
                THICKNESS = MATERIALS_INTERPOLATED_POINT%VALUES(MATERIALS_INTERPOLATED_POINT%INTERPOLATION_PARAMETERS% &
                  & FIELD_VARIABLE%NUMBER_OF_COMPONENTS,1)
              ENDIF
            ENDIF

            !Now add up the residual terms
            element_dof_idx=0
            DO component_idx=1,numberOfDimensions
              DEPENDENT_COMPONENT_INTERPOLATION_TYPE=DEPENDENT_FIELD%VARIABLES(var1)%COMPONENTS(component_idx)%INTERPOLATION_TYPE
              IF(DEPENDENT_COMPONENT_INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN !node based
                DEPENDENT_BASIS=>DEPENDENT_FIELD%VARIABLES(var1)%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY% &
                  & ELEMENTS%ELEMENTS(elementNumber)%BASIS
                NUMBER_OF_FIELD_COMPONENT_INTERPOLATION_PARAMETERS=DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                DO parameter_idx=1,NUMBER_OF_FIELD_COMPONENT_INTERPOLATION_PARAMETERS
                  element_dof_idx=element_dof_idx+1
                  DO component_idx2=1,numberOfDimensions
                    nonlinearMatrices%elementResidual%vector(element_dof_idx)= &
                      & nonlinearMatrices%elementResidual%vector(element_dof_idx)+ &
                      & GAUSS_WEIGHT*Jxxi*Jznu*THICKNESS*cauchyTensor(component_idx,component_idx2)* &
                      & DFDZ(parameter_idx,component_idx2,component_idx)
                  ENDDO ! component_idx2 (inner component index)
                ENDDO ! parameter_idx (residual vector loop)
              ELSEIF(DEPENDENT_COMPONENT_INTERPOLATION_TYPE==FIELD_ELEMENT_BASED_INTERPOLATION) THEN
                !Will probably never be used
                CALL FlagError("Finite elasticity with element based interpolation is not implemented.",err,error,*999)
              ENDIF
            ENDDO ! component_idx

            !Hydrostatic pressure component (skip for membrane problems)
            IF (EQUATIONS_SET_SUBTYPE /= EQUATIONS_SET_MEMBRANE_SUBTYPE) THEN
              HYDROSTATIC_PRESSURE_COMPONENT=DEPENDENT_FIELD%VARIABLES(var1)%NUMBER_OF_COMPONENTS
              DEPENDENT_COMPONENT_INTERPOLATION_TYPE=DEPENDENT_FIELD%VARIABLES(var1)%COMPONENTS(component_idx)%INTERPOLATION_TYPE
              IF(DEPENDENT_COMPONENT_INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN !node based
                COMPONENT_BASIS=>DEPENDENT_FIELD%VARIABLES(var1)%COMPONENTS(HYDROSTATIC_PRESSURE_COMPONENT)%DOMAIN% &
                  & TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS
                COMPONENT_QUADRATURE_SCHEME=>COMPONENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
                NUMBER_OF_FIELD_COMPONENT_INTERPOLATION_PARAMETERS=COMPONENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                DO parameter_idx=1,NUMBER_OF_FIELD_COMPONENT_INTERPOLATION_PARAMETERS
                  element_dof_idx=element_dof_idx+1 
                    nonlinearMatrices%elementResidual%vector(element_dof_idx)= &
                      & nonlinearMatrices%elementResidual%vector(element_dof_idx)+ &
                      & GAUSS_WEIGHT*Jxxi*COMPONENT_QUADRATURE_SCHEME%GAUSS_BASIS_FNS(parameter_idx,1,gauss_idx)* &
                      & (Jznu-1.0_DP-DARCY_VOL_INCREASE)
                ENDDO
              ELSEIF(DEPENDENT_COMPONENT_INTERPOLATION_TYPE==FIELD_ELEMENT_BASED_INTERPOLATION) THEN !element based
                element_dof_idx=element_dof_idx+1
                  nonlinearMatrices%elementResidual%vector(element_dof_idx)= &
                    & nonlinearMatrices%elementResidual%vector(element_dof_idx)+GAUSS_WEIGHT*Jxxi* &
                    & (Jznu-1.0_DP-DARCY_VOL_INCREASE)
              ENDIF
            ENDIF
          ENDDO !gauss_idx

          !Call surface pressure term here: should only be executed if THIS element has surface pressure on it (direct or incremented)
          IF(DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BOUNDARY_ELEMENT.AND. &
            & TOTAL_NUMBER_OF_SURFACE_PRESSURE_CONDITIONS>0) THEN    ! 
            CALL FiniteElasticity_SurfacePressureResidualEvaluate(EQUATIONS_SET,elementNumber,var1,var2,err,error,*999)
          ENDIF

        CASE (EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE, &
            & EQUATIONS_SET_COMPRESSIBLE_ACTIVECONTRACTION_SUBTYPE, &
            & EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
            & EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_STATIC_INRIA_SUBTYPE, &
            & EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_HOLMES_MOW_SUBTYPE, &
            & EQUATIONS_SET_ELASTICITY_FLUID_PRES_HOLMES_MOW_ACTIVE_SUBTYPE, &
            & EQUATIONS_SET_NEARLY_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE)
          !compressible problem (no pressure component)

          !Loop over gauss points and add up residuals
          DO gauss_idx=1,DEPENDENT_NUMBER_OF_GAUSS_POINTS
            GAUSS_WEIGHT=DEPENDENT_QUADRATURE_SCHEME%GAUSS_WEIGHTS(gauss_idx)

            !Interpolate fields at the gauss points
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & DEPENDENT_INTERPOLATED_POINT,err,error,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(DEPENDENT_BASIS%NUMBER_OF_XI,DEPENDENT_INTERPOLATED_POINT_METRICS, &
              & err,error,*999)
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & GEOMETRIC_INTERPOLATED_POINT,err,error,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,GEOMETRIC_INTERPOLATED_POINT_METRICS, &
              & err,error,*999)
            IF(ASSOCIATED(FIBRE_FIELD)) THEN
               CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
                & FIBRE_INTERPOLATED_POINT,err,error,*999)
            END IF
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & MATERIALS_INTERPOLATED_POINT,err,error,*999)
            IF(DARCY_DEPENDENT) THEN
              CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
                & DARCY_DEPENDENT_INTERPOLATED_POINT,err,error,*999) ! 'FIRST_PART_DERIV' required ???
            ENDIF

            !Calculate F=dZ/dNU at the gauss point
            CALL FiniteElasticity_GaussDeformationGradientTensor(DEPENDENT_INTERPOLATED_POINT_METRICS, &
              & GEOMETRIC_INTERPOLATED_POINT_METRICS,FIBRE_INTERPOLATED_POINT,DZDNU,ERR,ERROR,*999)

            Jxxi=GEOMETRIC_INTERPOLATED_POINT_METRICS%JACOBIAN

            !Calculate Cauchy stress tensor at the gauss point
            CALL FINITE_ELASTICITY_GAUSS_CAUCHY_TENSOR(EQUATIONS_SET,DEPENDENT_INTERPOLATED_POINT, &
              & MATERIALS_INTERPOLATED_POINT,GEOMETRIC_INTERPOLATED_POINT,DARCY_DEPENDENT_INTERPOLATED_POINT, &
              & INDEPENDENT_INTERPOLATED_POINT,cauchyTensor,Jznu,DZDNU,elementNumber,gauss_idx,err,error,*999)

            !Calculate dF/DZ at the gauss point
            CALL FINITE_ELASTICITY_GAUSS_DFDZ(DEPENDENT_INTERPOLATED_POINT,elementNumber,gauss_idx,numberOfDimensions, &
              & numberOfXi,DFDZ,err,error,*999)

            !Add up the residual terms
            element_dof_idx=0
            DO component_idx=1,DEPENDENT_NUMBER_OF_COMPONENTS
              DEPENDENT_COMPONENT_INTERPOLATION_TYPE=DEPENDENT_FIELD%VARIABLES(var1)%COMPONENTS(component_idx)%INTERPOLATION_TYPE
              IF(DEPENDENT_COMPONENT_INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN !node based
                DEPENDENT_BASIS=>DEPENDENT_FIELD%VARIABLES(var1)%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY% &
                  & ELEMENTS%ELEMENTS(elementNumber)%BASIS
                NUMBER_OF_FIELD_COMPONENT_INTERPOLATION_PARAMETERS=DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                DO parameter_idx=1,NUMBER_OF_FIELD_COMPONENT_INTERPOLATION_PARAMETERS  
                  element_dof_idx=element_dof_idx+1    
                  nonlinearMatrices%elementResidual%vector(element_dof_idx)= &
                    & nonlinearMatrices%elementResidual%vector(element_dof_idx)+ &
                    & GAUSS_WEIGHT*Jxxi*Jznu*(cauchyTensor(component_idx,1)*DFDZ(parameter_idx,1,component_idx)+ &
                    & cauchyTensor(component_idx,2)*DFDZ(parameter_idx,2,component_idx)+ &
                    & cauchyTensor(component_idx,3)*DFDZ(parameter_idx,3,component_idx))
                ENDDO
              ELSEIF(DEPENDENT_COMPONENT_INTERPOLATION_TYPE==FIELD_ELEMENT_BASED_INTERPOLATION) THEN
                !Will probably never be used
                CALL FlagError("Finite elasticity with element based interpolation is not implemented.",err,error,*999)
              ENDIF
            ENDDO !component_idx
          ENDDO !gauss_idx

          !Call surface pressure term here: should only be executed if THIS element has surface pressure on it (direct or incremented)
          IF(DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BOUNDARY_ELEMENT.AND. &
            & TOTAL_NUMBER_OF_SURFACE_PRESSURE_CONDITIONS>0) THEN    !
            CALL FiniteElasticity_SurfacePressureResidualEvaluate(EQUATIONS_SET,elementNumber,var1,var2,err,error,*999)
          ENDIF
        END SELECT
        !Gravity loading term
        IF(ASSOCIATED(rhsVector)) THEN
          IF(ASSOCIATED(SOURCE_FIELD)) THEN
            IF(ASSOCIATED(MATERIALS_FIELD%VARIABLE_TYPE_MAP(FIELD_V_VARIABLE_TYPE)%ptr)) THEN
              DENSITY_INTERPOLATION_PARAMETERS=>equations%interpolation%materialsInterpParameters(FIELD_V_VARIABLE_TYPE)%ptr
              CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,elementNumber, &
                & DENSITY_INTERPOLATION_PARAMETERS,err,error,*999)
              DENSITY_INTERPOLATED_POINT=>equations%interpolation%materialsInterpPoint(FIELD_V_VARIABLE_TYPE)%ptr
              IF(DARCY_DENSITY) THEN
                DARCY_MATERIALS_INTERPOLATION_PARAMETERS=>equations%interpolation%materialsInterpParameters( &
                  & FIELD_U1_VARIABLE_TYPE)%ptr
                CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,elementNumber, &
                  & DARCY_MATERIALS_INTERPOLATION_PARAMETERS,err,error,*999)
                DARCY_MATERIALS_INTERPOLATED_POINT=>equations%interpolation%materialsInterpPoint(FIELD_U1_VARIABLE_TYPE)%ptr
              ENDIF
              IF(rhsVector%updateVector) THEN
                SOURCE_INTERPOLATION_PARAMETERS=>equations%interpolation%sourceInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr
                CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,elementNumber, &
                  & SOURCE_INTERPOLATION_PARAMETERS,err,error,*999)
                SOURCE_INTERPOLATED_POINT=>equations%interpolation%sourceInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr

                DO gauss_idx=1,DEPENDENT_NUMBER_OF_GAUSS_POINTS
                  GAUSS_WEIGHT=DEPENDENT_QUADRATURE_SCHEME%GAUSS_WEIGHTS(gauss_idx)
                  CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
                    & SOURCE_INTERPOLATED_POINT,err,error,*999)
                  CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx,equations%interpolation% &
                    & geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                  CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
                    & DENSITY_INTERPOLATED_POINT,err,error,*999)
                  IF(DARCY_DENSITY) THEN
                    CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
                      & DARCY_MATERIALS_INTERPOLATED_POINT,err,error,*999)
                    !Account for separate fluid and solid proportions and densities
                    !Total lagrangian density = m_s + m_f = rho^0_s * (1 - phi^0) + rho_f * phi
                    !By assuming solid incompressibility, phi = (J - 1 + phi^0)
                    !\todo: Think about how this fits in with the constitutive relation, and what happens when the solid
                    !isn't incompressible. Can we assume the solid is incompressible if we aren't enforcing that in the
                    !constitutive relation?
                    DENSITY=DENSITY_INTERPOLATED_POINT%VALUES(1,1)*(1.0_DP-DARCY_MATERIALS_INTERPOLATED_POINT%VALUES(8,1)) + &
                      & DARCY_MATERIALS_INTERPOLATED_POINT%VALUES(7,1)*(Jznu-1.0_DP+DARCY_MATERIALS_INTERPOLATED_POINT%VALUES(8,1))
                  ELSE
                    DENSITY=DENSITY_INTERPOLATED_POINT%VALUES(1,1)
                  ENDIF
                  CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(DEPENDENT_BASIS%NUMBER_OF_XI, &
                    & DEPENDENT_INTERPOLATED_POINT_METRICS,err,error,*999)
                  element_dof_idx=0
                  DO component_idx=1,numberOfDimensions
                    DEPENDENT_BASIS=>DEPENDENT_FIELD%VARIABLES(var1)%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY% &
                      & ELEMENTS%ELEMENTS(elementNumber)%BASIS
                    DO parameter_idx=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                      element_dof_idx=element_dof_idx+1
                      rhsVector%elementVector%vector(element_dof_idx)=rhsVector%elementVector%vector(element_dof_idx) + &
                        & DENSITY*SOURCE_INTERPOLATED_POINT%VALUES(component_idx,1) * &
                        & DEPENDENT_QUADRATURE_SCHEME%GAUSS_BASIS_FNS(parameter_idx,NO_PART_DERIV,gauss_idx)*GAUSS_WEIGHT * &
                        & DEPENDENT_INTERPOLATED_POINT_METRICS%JACOBIAN
                    ENDDO
                  ENDDO
                ENDDO !gauss_idx
              ENDIF
            ENDIF
          ENDIF
        ELSE
          CALL FlagError("RHS vector is not associated.",err,error,*999)
        ENDIF

        !Scale factor adjustment
        IF(DEPENDENT_FIELD%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
          CALL Field_InterpolationParametersScaleFactorsElementGet(elementNumber,equations%interpolation% &
            & dependentInterpParameters(FIELD_VAR_TYPE)%ptr,err,error,*999)
          mhs=0          
          DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
            !Loop over element rows
            DEPENDENT_COMPONENT_INTERPOLATION_TYPE=DEPENDENT_FIELD%VARIABLES(FIELD_VAR_TYPE)%COMPONENTS(mh)%INTERPOLATION_TYPE
            IF(DEPENDENT_COMPONENT_INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN !node based
              DEPENDENT_BASIS=>DEPENDENT_FIELD%VARIABLES(FIELD_VAR_TYPE)%COMPONENTS(mh)%DOMAIN%TOPOLOGY% &
                & ELEMENTS%ELEMENTS(elementNumber)%BASIS
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1
                nonlinearMatrices%elementResidual%vector(mhs)=nonlinearMatrices%elementResidual%vector(mhs)* &
                  & EQUATIONS%INTERPOLATION%dependentInterpParameters(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ms,mh)
                IF(ASSOCIATED(rhsVector)) THEN
                  IF(ASSOCIATED(SOURCE_FIELD)) THEN
                   IF(rhsVector%updateVector) rhsVector%elementVector%vector(mhs)=rhsVector%elementVector%vector(mhs)* &
                      & EQUATIONS%INTERPOLATION%dependentInterpParameters(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ms,mh)
                  ENDIF
                ENDIF
              ENDDO !ms
            ENDIF
          ENDDO !mh
        ENDIF

      ELSE
        CALL FlagError("Equations set equations is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

    IF(DIAGNOSTICS5) THEN
      !Output element residual vector for first element
      IF(elementNumber == 1) THEN
        NDOFS = 0
        FIELD_VARIABLE=>DEPENDENT_FIELD%VARIABLES(var1) ! 'U' variable
        DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
          SELECT CASE(FIELD_VARIABLE%COMPONENTS(mh)%INTERPOLATION_TYPE)
          CASE(FIELD_NODE_BASED_INTERPOLATION)
            MESH_COMPONENT_1 = FIELD_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
            DEPENDENT_BASIS_1 => DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT_1)%ptr% &
              & TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS
            NDOFS = NDOFS + DEPENDENT_BASIS_1%NUMBER_OF_ELEMENT_PARAMETERS
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"EP: ",DEPENDENT_BASIS_1%NUMBER_OF_ELEMENT_PARAMETERS,err,error,*999)
          CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
            NDOFS = NDOFS + 1
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"EP: ",1,err,error,*999)
          CASE DEFAULT
            CALL FlagError("Interpolation type " &
              & //TRIM(NumberToVString(FIELD_VARIABLE%COMPONENTS(mh)%INTERPOLATION_TYPE,"*",err,error))// &
              & " is not valid for a finite elasticity equation.",err,error,*999)
          END SELECT
        END DO
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"NDOFS: ",NDOFS,err,error,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Element Vector for element number * (Fin.Elast.):",err,error,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Element Vector for element number (Fin.Elast.): ", &
          & elementNumber,err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,NDOFS,NDOFS,NDOFS,&
          & nonlinearMatrices%elementResidual%vector(:), &
          & '(4(X,E13.6))','4(4(X,E13.6))',err,error,*999)
      ENDIF
    ENDIF

    EXITS("FiniteElasticity_FiniteElementResidualEvaluate")
    RETURN
999 ERRORS("FiniteElasticity_FiniteElementResidualEvaluate",err,error)
    EXITS("FiniteElasticity_FiniteElementResidualEvaluate")
    RETURN 1

  END SUBROUTINE FiniteElasticity_FiniteElementResidualEvaluate

  !
  !================================================================================================================================
  !

  !>Pre-evaluates the residual for a finite elasticity finite element equations set.
  SUBROUTINE FiniteElasticity_FiniteElementPreResidualEvaluate(equationsSet,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(FIELD_TYPE), POINTER :: dependentField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: dependentVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("FiniteElasticity_FiniteElementPreResidualEvaluate",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(equationsSet%specification)) CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    IF(SIZE(equationsSet%specification,1)<3) &
      & CALL FlagError("Equations set specification must have at least three entries for a finite elasticity type equations set.", &
      & err,error,*999)
    
    SELECT CASE(equationsSet%specification(3))
      
    CASE(EQUATIONS_SET_CONSTITUTIVE_LAW_IN_CELLML_EVALUATE_SUBTYPE, &
      & EQUATIONS_SET_CONSTITUTIVE_AND_GROWTH_LAW_IN_CELLML_SUBTYPE)
      NULLIFY(dependentField)
      CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
      NULLIFY(dependentVariable)
      CALL Field_VariableGet(dependentField,FIELD_U1_VARIABLE_TYPE,dependentVariable,err,error,*999)
      CALL FiniteElasticity_StressStrainCalculate(equationsSet,EQUATIONS_SET_R_CAUCHY_GREEN_DEFORMATION_TENSOR, &
        & dependentVariable,err,error,*999)
    CASE(EQUATIONS_SET_MEMBRANE_SUBTYPE,EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE, &
      & EQUATIONS_SET_MOONEY_RIVLIN_ACTIVECONTRACTION_SUBTYPE, &
      & EQUATIONS_SET_STVENANT_KIRCHOFF_ACTIVECONTRACTION_SUBTYPE, &
      & EQUATIONS_SET_ISOTROPIC_EXPONENTIAL_SUBTYPE, &
      & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_EXPONENTIAL_SUBTYPE,EQUATIONS_SET_TRANSVERSE_ISOTROPIC_POLYNOMIAL_SUBTYPE, &
      & EQUATIONS_SET_ANISOTROPIC_POLYNOMIAL_SUBTYPE,EQUATIONS_SET_ANISOTROPIC_POLYNOMIAL_ACTIVE_SUBTYPE, &
      & EQUATIONS_SET_TRANS_ISOTROPIC_ACTIVE_TRANSITION_SUBTYPE,EQUATIONS_SET_TRANSVERSE_ISOTROPIC_ACTIVE_SUBTYPE, &
      & EQUATIONS_SET_ORTHOTROPIC_MATERIAL_COSTA_SUBTYPE,EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE,&
      & EQUATIONS_SET_COMPRESSIBLE_ACTIVECONTRACTION_SUBTYPE, &
      & EQUATIONS_SET_ACTIVECONTRACTION_SUBTYPE, EQUATIONS_SET_NO_SUBTYPE, &
      & EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
      & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
      & EQUATIONS_SET_ORTHOTROPIC_MATERIAL_HOLZAPFEL_OGDEN_SUBTYPE, &
      & EQUATIONS_SET_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE,EQUATIONS_SET_NEARLY_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE, &
      & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE,EQUATIONS_SET_TRANSVERSE_ISOTROPIC_GUCCIONE_SUBTYPE, &
      & EQUATIONS_SET_GUCCIONE_ACTIVECONTRACTION_SUBTYPE, EQUATIONS_SET_REFERENCE_STATE_TRANSVERSE_GUCCIONE_SUBTYPE, &
      & EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_STATIC_INRIA_SUBTYPE, &
      & EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_HOLMES_MOW_SUBTYPE, &
      & EQUATIONS_SET_ELASTICITY_FLUID_PRES_HOLMES_MOW_ACTIVE_SUBTYPE, &
      & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_HUMPHREY_YIN_SUBTYPE, &
      & EQUATIONS_SET_STANDARD_MONODOMAIN_ELASTICITY_SUBTYPE,EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE, &
      & EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE,EQUATIONS_SET_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE, &
      & EQUATIONS_SET_HOLZAPFEL_OGDEN_ACTIVECONTRACTION_SUBTYPE, &
      & EQUATIONS_SET_GROWTH_LAW_IN_CELLML_SUBTYPE, &
      & EQUATIONS_SET_MR_AND_GROWTH_LAW_IN_CELLML_SUBTYPE, &
      & EQUATIONS_SET_RATE_BASED_SMOOTH_MODEL_SUBTYPE,EQUATIONS_SET_COMPRESSIBLE_RATE_BASED_SMOOTH_MODEL_SUBTYPE, &
      & EQUATIONS_SET_RATE_BASED_GROWTH_MODEL_SUBTYPE,EQUATIONS_SET_COMPRESSIBLE_RATE_BASED_GROWTH_MODEL_SUBTYPE)
      !Do nothing ???
    CASE DEFAULT
      localError="The third equations set specification of "// &
        & TRIM(NumberToVString(equationsSet%specification(3),"*",err,error))// &
        & " is not valid for a finite elasticity type of an elasticity equation set."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("FiniteElasticity_FiniteElementPreResidualEvaluate")
    RETURN
999 ERRORS("FiniteElasticity_FiniteElementPreResidualEvaluate",err,error)
    EXITS("FiniteElasticity_FiniteElementPreResidualEvaluate")
    RETURN 1

  END SUBROUTINE FiniteElasticity_FiniteElementPreResidualEvaluate

  !
  !================================================================================================================================
  !

  !>Post-evaluates the residual for a finite elasticity finite element equations set.
  SUBROUTINE FiniteElasticity_FiniteElementPostResidualEvaluate(EQUATIONS_SET,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("FiniteElasticity_FiniteElementPostResidualEvaluate",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a finite elasticity type equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
      CASE(EQUATIONS_SET_MEMBRANE_SUBTYPE,EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE, &
        & EQUATIONS_SET_MOONEY_RIVLIN_ACTIVECONTRACTION_SUBTYPE, & 
        & EQUATIONS_SET_STVENANT_KIRCHOFF_ACTIVECONTRACTION_SUBTYPE, &
        & EQUATIONS_SET_ISOTROPIC_EXPONENTIAL_SUBTYPE, &
        & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_EXPONENTIAL_SUBTYPE,EQUATIONS_SET_STANDARD_MONODOMAIN_ELASTICITY_SUBTYPE, &
        & EQUATIONS_SET_ORTHOTROPIC_MATERIAL_COSTA_SUBTYPE,EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE,&
        & EQUATIONS_SET_COMPRESSIBLE_ACTIVECONTRACTION_SUBTYPE,EQUATIONS_SET_TRANSVERSE_ISOTROPIC_ACTIVE_SUBTYPE, &
        & EQUATIONS_SET_TRANS_ISOTROPIC_ACTIVE_TRANSITION_SUBTYPE,EQUATIONS_SET_TRANSVERSE_ISOTROPIC_POLYNOMIAL_SUBTYPE, & 
        & EQUATIONS_SET_ANISOTROPIC_POLYNOMIAL_SUBTYPE,EQUATIONS_SET_ANISOTROPIC_POLYNOMIAL_ACTIVE_SUBTYPE, &
        & EQUATIONS_SET_ACTIVECONTRACTION_SUBTYPE,EQUATIONS_SET_NO_SUBTYPE,EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE, &
        & EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE,EQUATIONS_SET_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE, &
        & EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
        & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
        & EQUATIONS_SET_ORTHOTROPIC_MATERIAL_HOLZAPFEL_OGDEN_SUBTYPE, &
        & EQUATIONS_SET_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE,EQUATIONS_SET_NEARLY_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE, & 
        & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE,EQUATIONS_SET_TRANSVERSE_ISOTROPIC_GUCCIONE_SUBTYPE, &
        & EQUATIONS_SET_GUCCIONE_ACTIVECONTRACTION_SUBTYPE,  EQUATIONS_SET_REFERENCE_STATE_TRANSVERSE_GUCCIONE_SUBTYPE,&
        & EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_STATIC_INRIA_SUBTYPE, &
        & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_HUMPHREY_YIN_SUBTYPE,&
        & EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_HOLMES_MOW_SUBTYPE, &
        & EQUATIONS_SET_ELASTICITY_FLUID_PRES_HOLMES_MOW_ACTIVE_SUBTYPE, &
        & EQUATIONS_SET_CONSTITUTIVE_LAW_IN_CELLML_EVALUATE_SUBTYPE, &
        & EQUATIONS_SET_HOLZAPFEL_OGDEN_ACTIVECONTRACTION_SUBTYPE, &
        & EQUATIONS_SET_CONSTITUTIVE_AND_GROWTH_LAW_IN_CELLML_SUBTYPE, &
        & EQUATIONS_SET_GROWTH_LAW_IN_CELLML_SUBTYPE, &
        & EQUATIONS_SET_MR_AND_GROWTH_LAW_IN_CELLML_SUBTYPE, &
        & EQUATIONS_SET_RATE_BASED_SMOOTH_MODEL_SUBTYPE,EQUATIONS_SET_COMPRESSIBLE_RATE_BASED_SMOOTH_MODEL_SUBTYPE, &
        & EQUATIONS_SET_RATE_BASED_GROWTH_MODEL_SUBTYPE,EQUATIONS_SET_COMPRESSIBLE_RATE_BASED_GROWTH_MODEL_SUBTYPE)
         !Do nothing ???
      CASE DEFAULT
        LOCAL_ERROR="The third equations set specification of "// &
          & TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
          & " is not valid for a finite elasticity type of an elasticity equation set."
        CALL FlagError(LOCAL_ERROR,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

    EXITS("FiniteElasticity_FiniteElementPostResidualEvaluate")
    RETURN
999 ERRORS("FiniteElasticity_FiniteElementPostResidualEvaluate",err,error)
    EXITS("FiniteElasticity_FiniteElementPostResidualEvaluate")
    RETURN 1

  END SUBROUTINE FiniteElasticity_FiniteElementPostResidualEvaluate

  !
  !================================================================================================================================
  !

  !>Calculated an output field for a finite elasticity equations set.
  SUBROUTINE FiniteElasticityEquationsSet_DerivedVariableCalculate(equationsSet,derivedType,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER, INTENT(IN) :: equationsSet !<A pointer to the equations set to calculate the output for
    INTEGER(INTG), INTENT(IN) :: derivedType !<The derived field type to calculate. \see EquationsSetConstants_DerivedTypes.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: derivedVariable

    ENTERS("FiniteElasticityEquationsSet_DerivedVariableCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(.NOT.equationsSet%EQUATIONS_SET_FINISHED) CALL FlagError("Equations set has not been finished.",err,error,*999)
    IF(.NOT.ASSOCIATED(equationsSet%equations)) CALL FlagError("Equations set equations are not associated.",err,error,*999)
    
    NULLIFY(derivedVariable)
    CALL Equations_DerivedVariableGet(equationsSet%equations,derivedType,derivedVariable,err,error,*999)
    CALL FiniteElasticity_StressStrainCalculate(equationsSet,derivedType,derivedVariable,err,error,*999)
   
    EXITS("FiniteElasticityEquationsSet_DerivedVariableCalculate")
    RETURN
999 ERRORS("FiniteElasticityEquationsSet_DerivedVariableCalculate",err,error)
    EXITS("FiniteElasticityEquationsSet_DerivedVariableCalculate")
    RETURN 1
    
  END SUBROUTINE FiniteElasticityEquationsSet_DerivedVariableCalculate

  !
  !================================================================================================================================
  !

  !>Calculates the strain field for a finite elasticity finite element equations set.
  SUBROUTINE FiniteElasticity_StressStrainCalculate(equationsSet,derivedType,fieldVariable,err,error,*)

    TYPE(EQUATIONS_SET_TYPE), POINTER, INTENT(IN) :: equationsSet !<A pointer to the equations set to calculate strain for.
    INTEGER(INTG), INTENT(IN) :: derivedType !<The type of derived field to calculate.     
    TYPE(FIELD_VARIABLE_TYPE), POINTER, INTENT(INOUT) :: fieldVariable !<The field variable to store the stress/strain in.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx,dependentNumberOfComponents,elementIdx,elementNumber,fieldVariableType,gaussIdx, &
      & meshComponentNumber,numberOfComponents,numberOfDimensions,numberOfGauss,numberOfTimes,numberOfXi,partIdx, &
      & startIdx,finishIdx,fieldInterpolation,dataPointNumber,numberOfDataPoints,dataPointIdx,residualVariableType, &
      & fieldVarType
    REAL(DP) :: dZdNu(3,3),Fg(3,3),Fe(3,3),J,Jg,Je,C(3,3),f(3,3),E(3,3),growthValues(3),xi(3),values(3,3)
    REAL(SP) :: elementUserElapsed,elementSystemElapsed,systemElapsed,systemTime1(1),systemTime2(1),systemTime3(1), &
      & systemTime4(1),userElapsed,userTime1(1),userTime2(1),userTime3(1),userTime4(1)
    TYPE(BASIS_TYPE), POINTER :: basis
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: coordinateSystem
    TYPE(DataProjectionType), POINTER :: dataProjection
    TYPE(DecompositionDataPointsType), POINTER :: dataPoints
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: decompositionTopology
    TYPE(DOMAIN_TYPE), POINTER :: domain
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: domainElements
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: elementsMappings
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: domainMappings
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: domainTopology
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsInterpolationType), POINTER :: interpolation
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FIELD_TYPE), POINTER :: field,dependentField,geometricField,fibreField,independentField,materialsField
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: quadratureScheme
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: geometricInterpolationParameters,dependentInterpolationParameters, &
      & fibreInterpolationParameters,independentInterpolationParameters,materialsInterpolationParameters
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: geometricInterpolatedPoint,dependentInterpolatedPoint,fibreInterpolatedPoint, &
      & independentInterpolatedPoint,materialsInterpolatedPoint
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER ::geometricInterpolatedPointMetrics,dependentInterpolatedPointMetrics
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: residualVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("FiniteElasticity_StressStrainCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is not associated.",err,error,*999)

    !Get the coordinate system
    NULLIFY(coordinateSystem)
    CALL EquationsSet_CoordinateSystemGet(equationsSet,coordinateSystem,err,error,*999)
    numberOfDimensions=coordinateSystem%NUMBER_OF_DIMENSIONS
    !Check the provided strain field variable has appropriate components and interpolation
    SELECT CASE(numberOfDimensions)
    CASE(3)
      numberOfComponents=6
    CASE(2)
      numberOfComponents=3
    CASE(1)
      numberOfComponents=1
    CASE DEFAULT
      CALL FlagError("The number of dimensions of "//TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
        & " is invalid.",err,error,*999)
    END SELECT
    NULLIFY(field)
    CALL FieldVariable_FieldGet(fieldVariable,field,err,error,*999)
    fieldVarType=fieldVariable%VARIABLE_TYPE
   
    CALL Field_NumberOfComponentsCheck(field,fieldVarType,6,err,error,*999)
    CALL Field_ComponentInterpolationGet(field,fieldVarType,1,fieldInterpolation,err,error,*999)
    !Check the interpolation type
    SELECT CASE(fieldInterpolation)
    CASE(FIELD_CONSTANT_INTERPOLATION)
      CALL FlagError("Can not calculate stress or strain for a field with constant interpolation.",err,error,*999)
    CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
      !OK
    CASE(FIELD_NODE_BASED_INTERPOLATION)
      CALL FlagError("Stress/strain calculation is not implemented for a field with node based interpolation.",err,error,*999)
    CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
      CALL FlagError("Stress/strain calculation is not implemented for a field with grid point based interpolation.", &
        & err,error,*999)
    CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
      !OK
    CASE(FIELD_DATA_POINT_BASED_INTERPOLATION)
      !OK
    CASE DEFAULT
      localError="The field interpolation type for component 1 of "//TRIM(NumberToVString(fieldInterpolation,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    !Check that all the components have the same interpolation type
    DO componentIdx=2,numberOfComponents
      CALL Field_ComponentInterpolationCheck(field,fieldVarType,componentIdx,fieldInterpolation,err,error,*999)
    ENDDO !componentIdx

    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMapping)
    CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
    NULLIFY(nonlinearMapping)
    CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)
    NULLIFY(residualVariable)
    CALL EquationsMappingNonlinear_ResidualVariableGet(nonlinearMapping,1,1,residualVariable,err,error,*999)
  
    NULLIFY(geometricField)
    CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
    NULLIFY(fibreField)
    CALL EquationsSet_FibreFieldExists(equationsSet,fibreField,err,error,*999)
    NULLIFY(materialsField)
    CALL EquationsSet_MaterialsFieldExists(equationsSet,materialsField,err,error,*999)
    NULLIFY(independentField)
    CALL EquationsSet_IndependentFieldExists(equationsSet,independentField,err,error,*999)
    dependentNumberOfComponents=residualVariable%NUMBER_OF_COMPONENTS

    NULLIFY(decomposition)
    CALL Field_DecompositionGet(dependentField,decomposition,err,error,*999)
    NULLIFY(decompositionTopology)
    CALL Decomposition_TopologyGet(decomposition,decompositionTopology,err,error,*999)
    NULLIFY(domain)
    CALL Decomposition_DomainGet(decomposition,0,domain,err,error,*999)
    NULLIFY(domainMappings)
    CALL Domain_MappingsGet(domain,domainMappings,err,error,*999)
    NULLIFY(elementsMappings)
    CALL DomainMappings_ElementsGet(domainMappings,elementsMappings,err,error,*999)
    NULLIFY(domainTopology)
    CALL Domain_TopologyGet(domain,domainTopology,err,error,*999)
    NULLIFY(domainElements)
    CALL DomainTopology_ElementsGet(domainTopology,domainElements,err,error,*999)
    
    !Grab interpolation points
    residualVariableType=residualVariable%VARIABLE_TYPE
    NULLIFY(interpolation)
    CALL Equations_InterpolationGet(equations,interpolation,err,error,*999)    
    NULLIFY(geometricInterpolationParameters)
    CALL EquationsInterpolation_GeometricParametersGet(interpolation,FIELD_U_VARIABLE_TYPE,geometricInterpolationParameters, &
      & err,error,*999)
    NULLIFY(geometricInterpolatedPoint)
    CALL EquationsInterpolation_GeometricPointGet(interpolation,FIELD_U_VARIABLE_TYPE,geometricInterpolatedPoint,err,error,*999)
    NULLIFY(geometricInterpolatedPointMetrics)
    CALL EquationsInterpolation_GeometricPointMetricsGet(interpolation,FIELD_U_VARIABLE_TYPE,geometricInterpolatedPointMetrics, &
      & err,error,*999)
    NULLIFY(dependentInterpolationParameters)
    CALL EquationsInterpolation_DependentParametersGet(interpolation,residualVariableType,dependentInterpolationParameters, &
      & err,error,*999)
    NULLIFY(dependentInterpolatedPoint)
    CALL EquationsInterpolation_DependentPointGet(interpolation,residualVariableType,dependentInterpolatedPoint,err,error,*999)
    NULLIFY(dependentInterpolatedPointMetrics)
    CALL EquationsInterpolation_DependentPointMetricsGet(interpolation,residualVariableType,dependentInterpolatedPointMetrics, &
      & err,error,*999)
    NULLIFY(fibreInterpolationParameters)
    NULLIFY(fibreInterpolatedPoint)
    IF(ASSOCIATED(fibreField)) THEN
      CALL EquationsInterpolation_FibreParametersGet(interpolation,FIELD_U_VARIABLE_TYPE,fibreInterpolationParameters, &
        & err,error,*999)
      CALL EquationsInterpolation_FibrePointGet(interpolation,FIELD_U_VARIABLE_TYPE,fibreInterpolatedPoint, &
        & err,error,*999)
    ENDIF
    NULLIFY(materialsInterpolationParameters)
    NULLIFY(materialsInterpolatedPoint)
    IF(ASSOCIATED(materialsField)) THEN
      CALL EquationsInterpolation_MaterialsParametersGet(interpolation,FIELD_U_VARIABLE_TYPE,materialsInterpolationParameters, &
        & err,error,*999)
      CALL EquationsInterpolation_MaterialsPointGet(interpolation,FIELD_U_VARIABLE_TYPE,materialsInterpolatedPoint, &
        & err,error,*999)
    ENDIF
    NULLIFY(independentInterpolationParameters)
    NULLIFY(independentInterpolatedPoint)
    IF(ASSOCIATED(independentField)) THEN
      CALL EquationsInterpolation_IndependentParametersGet(interpolation,FIELD_U_VARIABLE_TYPE,independentInterpolationParameters, &
        & err,error,*999)
      CALL EquationsInterpolation_IndependentPointGet(interpolation,FIELD_U_VARIABLE_TYPE,independentInterpolatedPoint, &
        & err,error,*999)
    ENDIF
 
    numberOfTimes=0

    !Loop over the two parts: 1 - boundary and ghost elements, 2 - internal
    DO partIdx=1,2          
      
      IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
        CALL Cpu_Timer(USER_CPU,userTime1,err,error,*999)
        CALL Cpu_timer(SYSTEM_CPU,systemTime1,err,error,*999)
      ENDIF
      
      IF(partIdx==1) THEN
        startIdx=elementsMappings%BOUNDARY_START
        finishIdx=elementsMappings%GHOST_FINISH
      ELSE
        startIdx=elementsMappings%INTERNAL_START
        finishIdx=elementsMappings%INTERNAL_FINISH
      ENDIF
      
      !Loop over (1) the boundary and ghost elements, (2) the internal elements
      DO elementIdx=startIdx,finishIdx
        
        numberOfTimes=numberOfTimes+1
        elementNumber=elementsMappings%DOMAIN_LIST(elementIdx)
        
        IF(diagnostics1) THEN
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Element number = ",elementNumber,err,error,*999)
        ENDIF

        NULLIFY(basis)
        CALL DomainElements_BasisGet(domainElements,elementNumber,basis,err,error,*999)
        numberOfXi=basis%NUMBER_OF_XI
                        
        CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,geometricInterpolationParameters, &
          & err,error,*999)
        CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,dependentInterpolationParameters, &
          & err,error,*999)
        IF(ASSOCIATED(fibreField)) THEN
          CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,fibreInterpolationParameters, &
            & err,error,*999)
        ENDIF
        IF(ASSOCIATED(materialsField)) THEN
          CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,materialsInterpolationParameters, &
            & err,error,*999)
        ENDIF
        IF(ASSOCIATED(independentField)) THEN
          CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,independentInterpolationParameters, &
            & err,error,*999)
        ENDIF

        SELECT CASE(fieldInterpolation)
        CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
          !Interpolate dependent, geometric, fibre etc. fields
          xi=[0.5_DP,0.5_DP,0.5_DP]
          CALL Field_InterpolateXi(FIRST_PART_DERIV,xi(1:numberOfXi),geometricInterpolatedPoint,err,error,*999)
          CALL Field_InterpolatedPointMetricsCalculate(numberOfXi,geometricInterpolatedPointMetrics,err,error,*999)
          CALL Field_InterpolateXi(FIRST_PART_DERIV,xi(1:numberOfXi),dependentInterpolatedPoint,err,error,*999)
          CALL Field_InterpolatedPointMetricsCalculate(numberOfXi,dependentInterpolatedPointMetrics,err,error,*999)
          IF(ASSOCIATED(fibreField)) &
            & CALL Field_InterpolateXi(FIRST_PART_DERIV,xi(1:numberOfXi),fibreInterpolatedPoint,err,error,*999)
          IF(ASSOCIATED(materialsField)) &
            & CALL Field_InterpolateXi(FIRST_PART_DERIV,xi(1:numberOfXi),materialsInterpolatedPoint,err,error,*999)
          IF(ASSOCIATED(independentField)) &
            & CALL Field_InterpolateXi(FIRST_PART_DERIV,xi(1:numberOfXi),independentInterpolatedPoint,err,error,*999)
!!\TODO how to get growth values????
          growthValues = [1.0_DP,1.0_DP,1.0_DP]
        
          CALL FiniteElasticity_StressStrainPoint(equationsSet,derivedType,numberOfDimensions,numberOfXi,gaussIdx, &
            & elementNumber,geometricInterpolatedPoint,geometricInterpolatedPointMetrics,dependentInterpolatedPoint, &
            & dependentInterpolatedPointMetrics,fibreInterpolatedPoint,materialsInterpolatedPoint,independentInterpolatedPoint, &
            & growthValues,values,err,error,*999)
          
          !We only want to store the independent components 
          SELECT CASE(numberOfDimensions)
          CASE(3)
            ! 3 dimensional problem
            ! ORDER OF THE COMPONENTS: U_11, U_12, U_13, U_22, U_23, U_33 (upper triangular matrix)
            CALL Field_ParameterSetUpdateLocalElement(field,fieldVarType,FIELD_VALUES_SET_TYPE, &
              & elementNumber,1,values(1,1),err,error,*999)
            CALL Field_ParameterSetUpdateLocalElement(field,fieldVarType,FIELD_VALUES_SET_TYPE, &
              & elementNumber,2,values(1,2),err,error,*999)
            CALL Field_ParameterSetUpdateLocalElement(field,fieldVarType,FIELD_VALUES_SET_TYPE, &
              & elementNumber,3,values(1,3),err,error,*999)
            CALL Field_ParameterSetUpdateLocalElement(field,fieldVarType,FIELD_VALUES_SET_TYPE, &
              & elementNumber,4,values(2,2),err,error,*999)
            CALL Field_ParameterSetUpdateLocalElement(field,fieldVarType,FIELD_VALUES_SET_TYPE, &
              & elementNumber,5,values(2,3),err,error,*999)
            CALL Field_ParameterSetUpdateLocalElement(field,fieldVarType,FIELD_VALUES_SET_TYPE, &
              & elementNumber,6,values(3,3),err,error,*999)
          CASE(2)
            ! 2 dimensional problem
            ! ORDER OF THE COMPONENTS: U_11, U_12, U_22 (upper triangular matrix)
            CALL Field_ParameterSetUpdateLocalElement(field,fieldVarType,FIELD_VALUES_SET_TYPE, &
              & elementNumber,1,values(1,1),err,error,*999)
            CALL Field_ParameterSetUpdateLocalElement(field,fieldVarType,FIELD_VALUES_SET_TYPE, &
              & elementNumber,2,values(1,2),err,error,*999)
            CALL Field_ParameterSetUpdateLocalElement(field,fieldVarType,FIELD_VALUES_SET_TYPE, &
              & elementNumber,3,values(2,2),err,error,*999)
          CASE(1)
            ! 1 dimensional problem
            CALL Field_ParameterSetUpdateLocalElement(field,fieldVarType,FIELD_VALUES_SET_TYPE, &
              & elementNumber,1,values(1,1),err,error,*999)
          CASE DEFAULT
            localError="The number of dimensions of "//TRIM(NumberToVString(numberofDimensions,"*",err,error))// &
                & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
          
        CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)            
            
          NULLIFY(quadratureScheme)               
          CALL Basis_QuadratureSchemeGet(basis,BASIS_DEFAULT_QUADRATURE_SCHEME,quadratureScheme,err,error,*999)
          numberOfGauss=quadratureScheme%NUMBER_OF_GAUSS
          
          !Loop over gauss points        
          DO gaussIdx=1,numberOfGauss
            
            IF(diagnostics1) THEN
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Gauss point number = ",gaussIdx,err,error,*999)
            ENDIF
            
            !Interpolate dependent, geometric, fibre etc. fields
            CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussIdx,dependentInterpolatedPoint, &
              & err,error,*999)
            CALL Field_InterpolatedPointMetricsCalculate(numberOfXi,dependentInterpolatedPointMetrics,err,error,*999)
            CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussIdx,geometricInterpolatedPoint, &
              & err,error,*999)
            CALL Field_InterpolatedPointMetricsCalculate(numberOfXi,geometricInterpolatedPointMetrics,err,error,*999)
            IF(ASSOCIATED(fibreField)) THEN
              CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussIdx,fibreInterpolatedPoint, &
                & err,error,*999)
            ENDIF
            IF(ASSOCIATED(materialsField)) THEN
              CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussIdx,materialsInterpolatedPoint, &
                & err,error,*999)
            ENDIF
            IF(ASSOCIATED(independentField)) THEN
              CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussIdx,independentInterpolatedPoint, &
                & err,error,*999)
            ENDIF
            IF(equationsSet%specification(3)==EQUATIONS_SET_CONSTITUTIVE_AND_GROWTH_LAW_IN_CELLML_SUBTYPE.OR. &
              equationsSet%specification(3)==EQUATIONS_SET_MR_AND_GROWTH_LAW_IN_CELLML_SUBTYPE) THEN
              CALL Field_ParameterSetGetLocalGaussPoint(dependentField,FIELD_U3_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & gaussIdx,elementNumber,1,growthValues(1),err,error,*999)
              IF(numberofDimensions>1) THEN
                CALL Field_ParameterSetGetLocalGaussPoint(dependentField,FIELD_U3_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & gaussIdx,elementNumber,2,growthValues(2),err,error,*999)
                IF(numberOfDimensions>2) THEN
                  CALL Field_ParameterSetGetLocalGaussPoint(dependentField,FIELD_U3_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & gaussIdx,elementNumber,3,growthValues(3),err,error,*999)
                ENDIF
              ENDIF
            ELSE
              growthValues = [1.0_DP,1.0_DP,1.0_DP]
            ENDIF
            
            CALL FiniteElasticity_StressStrainPoint(equationsSet,derivedType,numberOfDimensions,numberOfXi,gaussIdx, &
              & elementNumber,geometricInterpolatedPoint,geometricInterpolatedPointMetrics,dependentInterpolatedPoint, &
              & dependentInterpolatedPointMetrics,fibreInterpolatedPoint,materialsInterpolatedPoint,independentInterpolatedPoint, &
              & growthValues,values,err,error,*999)
            
            !We only want to store the independent components 
            SELECT CASE(numberOfDimensions)
            CASE(3)
              ! 3 dimensional problem
              ! ORDER OF THE COMPONENTS: U_11, U_12, U_13, U_22, U_23, U_33 (upper triangular matrix)
              CALL Field_ParameterSetUpdateLocalGaussPoint(field,fieldVarType,FIELD_VALUES_SET_TYPE, &
                & gaussIdx,elementNumber,1,values(1,1),err,error,*999)
              CALL Field_ParameterSetUpdateLocalGaussPoint(field,fieldVarType,FIELD_VALUES_SET_TYPE, &
                & gaussIdx,elementNumber,2,values(1,2),err,error,*999)
              CALL Field_ParameterSetUpdateLocalGaussPoint(field,fieldVarType,FIELD_VALUES_SET_TYPE, &
                & gaussIdx,elementNumber,3,values(1,3),err,error,*999)
              CALL Field_ParameterSetUpdateLocalGaussPoint(field,fieldVarType,FIELD_VALUES_SET_TYPE, &
                & gaussIdx,elementNumber,4,values(2,2),err,error,*999)
              CALL Field_ParameterSetUpdateLocalGaussPoint(field,fieldVarType,FIELD_VALUES_SET_TYPE, &
                & gaussIdx,elementNumber,5,values(2,3),err,error,*999)
              CALL Field_ParameterSetUpdateLocalGaussPoint(field,fieldVarType,FIELD_VALUES_SET_TYPE, &
                & gaussIdx,elementNumber,6,values(3,3),err,error,*999)
            CASE(2)
              ! 2 dimensional problem
              ! ORDER OF THE COMPONENTS: U_11, U_12, U_22 (upper triangular matrix)
              CALL Field_ParameterSetUpdateLocalGaussPoint(field,fieldVarType,FIELD_VALUES_SET_TYPE, &
                & gaussIdx,elementNumber,1,values(1,1),err,error,*999)
              CALL Field_ParameterSetUpdateLocalGaussPoint(field,fieldVarType,FIELD_VALUES_SET_TYPE, &
                & gaussIdx,elementNumber,2,values(1,2),err,error,*999)
              CALL Field_ParameterSetUpdateLocalGaussPoint(field,fieldVarType,FIELD_VALUES_SET_TYPE, &
                & gaussIdx,elementNumber,3,values(2,2),err,error,*999)
            CASE(1)
              ! 1 dimensional problem
              CALL Field_ParameterSetUpdateLocalGaussPoint(field,fieldVarType,FIELD_VALUES_SET_TYPE, &
                & gaussIdx,elementNumber,1,values(1,1),err,error,*999)
            CASE DEFAULT
              localError="The number of dimensions of "//TRIM(NumberToVString(numberofDimensions,"*",err,error))// &
                & " is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ENDDO !gaussIdx/
        CASE(FIELD_DATA_POINT_BASED_INTERPOLATION)
          
          NULLIFY(dataProjection)
          CALL Field_DataProjectionGet(field,dataProjection,err,error,*999)
          NULLIFY(dataPoints)
          CALL DecompositionTopology_DataPointsGet(decompositionTopology,dataPoints,err,error,*999)

          numberOfDataPoints=dataPoints%elementDataPoint(elementNumber)%numberOfProjectedData

          DO dataPointIdx=1,numberOfDataPoints
            
            dataPointNumber=dataPoints%elementDataPoint(elementNumber)%dataIndices(dataPointIdx)%globalNumber
            xi(1:numberOfXi)=dataProjection%dataProjectionResults(dataPointNumber)%elementXi(1:numberOfXi)
            CALL Field_InterpolateXi(FIRST_PART_DERIV,xi(1:numberOfXi),geometricInterpolatedPoint,err,error,*999)
            CALL Field_InterpolatedPointMetricsCalculate(numberOfXi,geometricInterpolatedPointMetrics,err,error,*999)
            CALL Field_InterpolateXi(FIRST_PART_DERIV,xi(1:numberOfXi),dependentInterpolatedPoint,err,error,*999)
            CALL Field_InterpolatedPointMetricsCalculate(numberOfXi,dependentInterpolatedPointMetrics,err,error,*999)
            IF(ASSOCIATED(fibreField)) &
              & CALL Field_InterpolateXi(FIRST_PART_DERIV,xi(1:numberOfXi),fibreInterpolatedPoint,err,error,*999)
            IF(ASSOCIATED(materialsField)) &
              & CALL Field_InterpolateXi(FIRST_PART_DERIV,xi(1:numberOfXi),materialsInterpolatedPoint,err,error,*999)
            IF(ASSOCIATED(independentField)) &
              & CALL Field_InterpolateXi(FIRST_PART_DERIV,xi(1:numberOfXi),independentInterpolatedPoint,err,error,*999)
!!\TODO how to get growth values????
            growthValues = [1.0_DP,1.0_DP,1.0_DP]
            
            CALL FiniteElasticity_StressStrainPoint(equationsSet,derivedType,numberOfDimensions,numberOfXi,gaussIdx, &
              & elementNumber,geometricInterpolatedPoint,geometricInterpolatedPointMetrics,dependentInterpolatedPoint, &
              & dependentInterpolatedPointMetrics,fibreInterpolatedPoint,materialsInterpolatedPoint,independentInterpolatedPoint, &
              & growthValues,values,err,error,*999)
            
            !We only want to store the independent components 
            SELECT CASE(numberOfDimensions)
            CASE(3)
              ! 3 dimensional problem
              ! ORDER OF THE COMPONENTS: U_11, U_12, U_13, U_22, U_23, U_33 (upper triangular matrix)
              CALL Field_ParameterSetUpdateLocalElement(field,fieldVarType,FIELD_VALUES_SET_TYPE, &
                & elementNumber,1,values(1,1),err,error,*999)
              CALL Field_ParameterSetUpdateLocalElement(field,fieldVarType,FIELD_VALUES_SET_TYPE, &
                & elementNumber,2,values(1,2),err,error,*999)
              CALL Field_ParameterSetUpdateLocalElement(field,fieldVarType,FIELD_VALUES_SET_TYPE, &
                & elementNumber,3,values(1,3),err,error,*999)
              CALL Field_ParameterSetUpdateLocalElement(field,fieldVarType,FIELD_VALUES_SET_TYPE, &
                & elementNumber,4,values(2,2),err,error,*999)
              CALL Field_ParameterSetUpdateLocalElement(field,fieldVarType,FIELD_VALUES_SET_TYPE, &
                & elementNumber,5,values(2,3),err,error,*999)
              CALL Field_ParameterSetUpdateLocalElement(field,fieldVarType,FIELD_VALUES_SET_TYPE, &
              & elementNumber,6,values(3,3),err,error,*999)
            CASE(2)
              ! 2 dimensional problem
              ! ORDER OF THE COMPONENTS: U_11, U_12, U_22 (upper triangular matrix)
              CALL Field_ParameterSetUpdateLocalElement(field,fieldVarType,FIELD_VALUES_SET_TYPE, &
                & elementNumber,1,values(1,1),err,error,*999)
              CALL Field_ParameterSetUpdateLocalElement(field,fieldVarType,FIELD_VALUES_SET_TYPE, &
                & elementNumber,2,values(1,2),err,error,*999)
              CALL Field_ParameterSetUpdateLocalElement(field,fieldVarType,FIELD_VALUES_SET_TYPE, &
                & elementNumber,3,values(2,2),err,error,*999)
            CASE(1)
              ! 1 dimensional problem
              CALL Field_ParameterSetUpdateLocalElement(field,fieldVarType,FIELD_VALUES_SET_TYPE, &
                & elementNumber,1,values(1,1),err,error,*999)
            CASE DEFAULT
              localError="The number of dimensions of "//TRIM(NumberToVString(numberofDimensions,"*",err,error))// &
                & " is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
                      
          ENDDO !dataPointIdx
          
        CASE DEFAULT
          localError="The field interpolation type for component 1 of "// &
            & TRIM(NumberToVString(fieldInterpolation,"*",err,error))//" is invalid or not implemented."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        
      ENDDO !elementIdx
      
      !Output timing information if required
      IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
        CALL Cpu_Timer(USER_CPU,userTime2,err,error,*999)
        CALL Cpu_Timer(SYSTEM_CPU,systemTime2,err,error,*999)
        userElapsed=userTime2(1)-userTime1(1)
        systemElapsed=systemTime2(1)-systemTime1(1)
        elementUserElapsed=elementUserElapsed+userElapsed
        elementSystemElapsed=elementSystemElapsed+systemElapsed
        IF(partIdx==1) THEN
          CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
          CALL Profiling_TimingsOutput(1,"Boundary+ghost elements calculation",userElapsed,systemElapsed,err,error,*999)
       ELSE
          CALL Profiling_TimingsOutput(1,"Internal elements calculation",userElapsed,systemElapsed,err,error,*999)
          IF(numberOfTimes>0) CALL Profiling_TimingsOutput(1,"Average element calculation", &
            & elementUserElapsed/numberOfTimes,elementSystemElapsed/numberOfTimes,err,error,*999)
        ENDIF
      ENDIF !equations%outputType>=EQUATIONS_TIMING_OUTPUT
      
      IF(partIdx==1) THEN
        IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
          CALL Cpu_Timer(USER_CPU,userTime3,err,error,*999)
          CALL Cpu_Timer(SYSTEM_CPU,systemTime3,err,error,*999)
        ENDIF
        !Start to update the field
        CALL Field_ParameterSetUpdateStart(field,fieldVarType,FIELD_VALUES_SET_TYPE,err,error,*999)
      ELSE
        !Finish to update the field
        CALL Field_ParameterSetUpdateFinish(field,fieldVarType,FIELD_VALUES_SET_TYPE,err,error,*999)
        !Output timing information if required
        IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
          CALL Cpu_Timer(USER_CPU,userTime4,err,error,*999)
          CALL Cpu_Timer(SYSTEM_CPU,systemTime4,err,error,*999)
          userElapsed=userTime4(1)-userTime3(1)
          systemElapsed=systemTime4(1)-systemTime3(1)
          CALL Profiling_TimingsOutput(1,"Parameters update transfer",userElapsed,systemElapsed,err,error,*999)
        ENDIF !equations%outputType>=EQUATIONS_TIMING_OUTPUT
      ENDIF
      
    ENDDO !partIdx
        
    EXITS("FiniteElasticity_StressStrainCalculate")
    RETURN
999 ERRORSEXITS("FiniteElasticity_StressStrainCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE FiniteElasticity_StressStrainCalculate

  !
  !================================================================================================================================
  !

  !>Evaluates stress and strain at a point. \TODO merge this with interpolate xi below.
  SUBROUTINE FiniteElasticity_StressStrainPoint(equationsSet,evaluateType,numberOfDimensions,numberOfXi,pointNumber, &
    & elementNumber,geometricInterpolatedPoint,geometricInterpolatedPointMetrics,dependentInterpolatedPoint, &
    & dependentInterpolatedPointMetrics,fibreInterpolatedPoint,materialsInterpolatedPoint,independentInterpolatedPoint, &
    & growthValues,values,err,error,*)
    ! Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER, INTENT(IN) :: equationsSet !<A pointer to the equations set to calculate the tensor for
    INTEGER(INTG), INTENT(IN) :: evaluateType !<The type of tensor to evaluate.
    INTEGER(INTG), INTENT(IN) :: numberOfDimensions !<The number of dimensions
    INTEGER(INTG), INTENT(IN) :: numberOfXi !<The number of xi directions
    INTEGER(INTG), INTENT(IN) :: pointNumber !<The point number to evaluate the tensor for
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The user element number to evaluate the tensor for
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: geometricInterpolatedPoint !<The geometric interpolated point
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: geometricInterpolatedPointMetrics !<The geometric interpolated point metrics
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: dependentInterpolatedPoint !<The dependent interpolated point
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: dependentInterpolatedPointMetrics !<The dependent interpolated point metrics
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: fibreInterpolatedPoint !<The fibre interpolated point
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: materialsInterpolatedPoint !<The materials interpolated point
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: independentInterpolatedPoint !<The independent interpolated point
    REAL(DP), INTENT(IN) :: growthValues(:) !<The growth extension values if any. 
    REAL(DP), INTENT(OUT) :: values(:,:) !<On exit, the interpolated tensor values.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    ! Local variables
    INTEGER(INTG) :: fieldInterpolation,i,nh,mh
    REAL(DP) :: b(3,3),C(3,3),dZdNu(3,3),Fe(3,3),Fg(3,3),Je,Jg,Jznu
    REAL(DP) :: E(3,3),cauchyStressTensor(3,3),cauchyStressVoigt(6)
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: darcyInterpolatedPoint
    TYPE(VARYING_STRING) :: localError

    ENTERS("FiniteElasticity_StressStrainPoint",err,error,*999)

    !Calculate field metrics
    CALL Field_InterpolatedPointMetricsCalculate(numberOfXi,geometricInterpolatedPointMetrics,err,error,*999)
    CALL Field_InterpolatedPointMetricsCalculate(numberOfXi,dependentInterpolatedPointMetrics,err,error,*999)

    !Calculate F=dZ/dNU, the deformation gradient tensor at the xi location
    CALL FiniteElasticity_GaussDeformationGradientTensor(dependentInterpolatedPointMetrics, &
      & geometricInterpolatedPointMetrics,fibreInterpolatedPoint,dZdNu,err,error,*999)

    CALL FiniteElasticity_GaussGrowthTensor(equationsSet,numberOfDimensions,dZdNu,growthValues,Fg,Fe,Jg,Je, &
      & err,error,*999)
    
    IF(evaluateType==EQUATIONS_SET_R_CAUCHY_GREEN_DEFORMATION_TENSOR .OR. &
      & evaluateType==EQUATIONS_SET_GREEN_LAGRANGE_STRAIN_TENSOR) THEN
      CALL MatrixTransposeProduct(Fe(1:numberOfDimensions,1:numberOfXi),Fe(1:numberOfDimensions,1:numberOfXi), &
        & C(1:numberOfDimensions,1:numberOfDimensions),err,error,*999)
    ENDIF
    
    IF(evaluateType==EQUATIONS_SET_L_CAUCHY_GREEN_DEFORMATION_TENSOR) THEN
      CALL MatrixProductTranspose(Fe(1:numberOfDimensions,1:numberOfXi),Fe(1:numberOfDimensions,1:numberOfXi), &
        & b(1:numberOfDimensions,1:numberOfDimensions),err,error,*999)
    ENDIF

    IF(evaluateType==EQUATIONS_SET_GREEN_LAGRANGE_STRAIN_TENSOR) THEN
      !Calculate E
      E(1:numberOfDimensions,1:numberOfDimensions)=0.5_DP*C(1:numberOfDimensions,1:numberOfDimensions)
      DO i=1,numberOfDimensions
        E(i,i)=E(i,i)-0.5_DP
      ENDDO !i
    ENDIF

    IF(evaluateType==EQUATIONS_SET_CAUCHY_STRESS_TENSOR) THEN
!!\TODO the whole stress thing needs to be looked at as the routines below do not take in the deformation gradient that
!! is calculated above but rather they calculate it internally. This will lead to mismatches as things like growth are
!! not taken into account. 
      
      SELECT CASE(equationsSet%specification(3))
      CASE(EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE, &
        & EQUATIONS_SET_MR_AND_GROWTH_LAW_IN_CELLML_SUBTYPE)
        !Calculate the Cauchy stress tensor (in Voigt form) at the gauss point.
        Jznu=dependentInterpolatedPointMetrics%JACOBIAN/geometricInterpolatedPointMetrics%JACOBIAN
        ! Note that some problems, e.g. active contraction, require additonal fields to be evaluated at Gauss points. This is
        ! currently achieved by providing the gausspoint number to the FINITE_ELASTICITY_GAUSS_STRESS_TENSOR routine.
        ! However, the current  routine, FiniteElasticity_TensorInterpolateXi, aims to evaluate tensors as any xi, so the Gauss
        ! point number has been set to 0, which will generate an error for such problems.
        ! To address such issues, the FINITE_ELASTICITY_GAUSS_STRESS_TENSOR routine needs to be generalized to allow calculation
        ! of stress at any xi position and the GaussPoint number argument needs to be replace with a set of xi coordinates.
        CALL FINITE_ELASTICITY_GAUSS_STRESS_TENSOR(equationsSet,dependentInterpolatedPoint, &
          & materialsInterpolatedPoint,geometricInterpolatedPoint,cauchyStressVoigt,dZdNu,Jznu, &
          & elementNumber,0,ERR,ERROR,*999)
        
        !Convert from Voigt form to tensor form. \TODO needs to be generalised for 2D
        DO nh=1,numberOfDimensions
          DO mh=1,numberOfDimensions
            cauchyStressTensor(mh,nh)=cauchyStressVoigt(TENSOR_TO_VOIGT(mh,nh,numberOfDimensions))
          ENDDO
        ENDDO
      CASE(EQUATIONS_SET_ORTHOTROPIC_MATERIAL_COSTA_SUBTYPE, EQUATIONS_SET_TRANSVERSE_ISOTROPIC_GUCCIONE_SUBTYPE)
        CALL FINITE_ELASTICITY_GAUSS_CAUCHY_TENSOR(equationsSet,dependentInterpolatedPoint,materialsInterpolatedPoint, &
          & geometricInterpolatedPoint,darcyInterpolatedPoint,independentInterpolatedPoint, &
          & cauchyStressTensor,Jznu,dZdNu,elementNumber,0,ERR,ERROR,*999)
      CASE DEFAULT
        CALL FlagError("Not implemented ",err,error,*999)
      END SELECT
    END IF

    SELECT CASE(evaluateType)
    CASE(EQUATIONS_SET_DEFORMATION_GRADIENT_TENSOR)
      values(1:numberOfDimensions,1:numberOfDimensions)=Fe(1:numberOfDimensions,1:numberOfDimensions)
    CASE(EQUATIONS_SET_R_CAUCHY_GREEN_DEFORMATION_TENSOR)
      values(1:numberOfDimensions,1:numberOfDimensions)=C(1:numberOfDimensions,1:numberOfDimensions)
    CASE(EQUATIONS_SET_L_CAUCHY_GREEN_DEFORMATION_TENSOR)
      values(1:numberOfDimensions,1:numberOfDimensions)=b(1:numberOfDimensions,1:numberOfDimensions)
    CASE(EQUATIONS_SET_GREEN_LAGRANGE_STRAIN_TENSOR)
      values(1:numberOfDimensions,1:numberOfDimensions)=E(1:numberOfDimensions,1:numberOfDimensions)
    CASE(EQUATIONS_SET_CAUCHY_STRESS_TENSOR)
      values(1:numberOfDimensions,1:numberOfDimensions)=cauchyStressTensor(1:numberOfDimensions,1:numberOfDimensions)
    CASE(EQUATIONS_SET_SECOND_PK_STRESS_TENSOR)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      CALL FlagError("The tensor evalaute type of "//TRIM(NumberToVString(evaluateType,"*",err,error))//" is invalid "// &
        & "for finite elasticity equation sets.",err,error,*999)
    END SELECT
 
    EXITS("FiniteElasticity_StressStrainPoint")
    RETURN
999 ERRORS("FiniteElasticity_StressStrainPoint",err,error)
    EXITS("FiniteElasticity_StressStrainPoint")
    RETURN 1
    
  END SUBROUTINE FiniteElasticity_StressStrainPoint

  !
  !================================================================================================================================
  !

  !>Evaluates a tensor at a given element Gauss point. \TODO merge this with interpolate xi below.
  SUBROUTINE FiniteElasticity_TensorInterpolateGaussPoint(equationsSet,tensorEvaluateType,gaussPointNumber,userElementNumber, &
    & values,err,error,*)
    ! Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER, INTENT(IN) :: equationsSet !<A pointer to the equations set to calculate the tensor for
    INTEGER(INTG), INTENT(IN) :: tensorEvaluateType !<The type of tensor to evaluate.
    INTEGER(INTG), INTENT(IN) :: gaussPointNumber !<The Gauss point number to evaluate the tensor for
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user element number to evaluate the tensor for
    REAL(DP), INTENT(OUT) :: values(:,:) !<On exit, the interpolated tensor values.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    ! Local variables
    INTEGER(INTG) :: dependentVarType,meshComponentNumber
    INTEGER(INTG) :: numberOfDimensions,numberOfXi
    INTEGER(INTG) :: localElementNumber,i,nh,mh
    REAL(DP) :: C(3,3),dZdNu(3,3),E(3,3),cauchyStressTensor(3,3),cauchyStressVoigt(6),Fe(3,3),Fg(3,3),growthValues(3),Je,Jg,Jznu
    LOGICAL :: userElementExists,ghostElement
    TYPE(BASIS_TYPE), POINTER :: elementBasis
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: coordinateSystem
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: decompositionTopology
    TYPE(DOMAIN_TYPE), POINTER :: domain
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: domainTopology
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FIELD_TYPE), POINTER :: dependentField
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: geometricInterpolatedPoint, &
      & fibreInterpolatedPoint,dependentInterpolatedPoint,materialsInterpolatedPoint, &
      & independentInterpolatedPoint,darcyInterpolatedPoint
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: geometricInterpolatedPointMetrics, &
      & dependentInterpolatedPointMetrics
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: residualVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("FiniteElasticity_TensorInterpolateGaussPoint",err,error,*999)

    NULLIFY(geometricInterpolatedPoint)
    NULLIFY(fibreInterpolatedPoint)
    NULLIFY(dependentInterpolatedPoint)
    NULLIFY(materialsInterpolatedPoint)
    NULLIFY(independentInterpolatedPoint)
    NULLIFY(darcyInterpolatedPoint)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
    NULLIFY(coordinateSystem)
    CALL Field_CoordinateSystemGet(dependentField,coordinateSystem,err,error,*999)
    numberOfDimensions=coordinateSystem%number_of_dimensions
    IF(SIZE(values,1)<numberOfDimensions) THEN
      localError="The size of the first dimension of the supplied values array of "// &
        & TRIM(NumberToVString(SIZE(values,1),"*",err,error))//" is too small. The size must be >= "// &
        & TRIM(NumberToVString(numberOfDimensions,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(values,2)<numberOfDimensions) THEN
      localError="The size of the second dimension of the supplied values array of "// &
        & TRIM(NumberToVString(SIZE(values,2),"*",err,error))//" is too small. The size must be >= "// &
        & TRIM(NumberToVString(numberOfDimensions,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF    
    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMapping)
    CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
    NULLIFY(nonlinearMapping)
    CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)
    NULLIFY(residualVariable)
    CALL EquationsMappingNonlinear_ResidualVariableGet(nonlinearMapping,1,1,residualVariable,err,error,*999)
    dependentVarType=residualVariable%VARIABLE_TYPE
    NULLIFY(decomposition)
    CALL Field_DecompositionGet(dependentField,decomposition,err,error,*999)
    NULLIFY(decompositionTopology)
    CALL Decomposition_TopologyGet(decomposition,decompositionTopology,err,error,*999)
    CALL DecompositionTopology_ElementCheckExists(decompositionTopology,userElementNumber,userElementExists, &
      & localElementNumber,ghostElement,err,error,*999)
    IF(.NOT.userElementExists) THEN
      localError="The specified user element number of "//TRIM(NumberToVstring(userElementNumber,"*",err,error))// &
        & " does not exist in the decomposition for the dependent field."
      CALL FlagError(localError,err,error,*999)
    END IF    
    NULLIFY(domain)
    CALL Decomposition_DomainGet(decomposition,0,domain,err,error,*999)
    NULLIFY(domainTopology)
    CALL Domain_TopologyGet(domain,domainTopology,err,error,*999)
    NULLIFY(elementBasis)
    CALL DomainTopology_ElementBasisGet(domainTopology,userElementNumber,elementBasis,err,error,*999)
    numberOfXi=elementBasis%number_of_xi
    
    IF(.NOT.ASSOCIATED(equations%interpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)    

    !Get the interpolation parameters for this element
    CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,localElementNumber, &
      & equations%interpolation%geometricInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
    CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,localElementNumber, &
      & equations%interpolation%dependentInterpParameters(dependentVarType)%ptr,err,error,*999)
    IF(ASSOCIATED(equations%interpolation%fibreInterpParameters)) THEN
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,localElementNumber, &
        & equations%interpolation%fibreInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
    END IF

    !Get interpolated points
    geometricInterpolatedPoint=>equations%interpolation%geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr
    dependentInterpolatedPoint=>equations%interpolation%dependentInterpPoint(dependentVarType)%ptr
    IF(ASSOCIATED(equations%interpolation%fibreInterpPoint)) THEN
      fibreInterpolatedPoint=>equations%interpolation%fibreInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr
    END IF

    !Get interpolated point metrics
    geometricInterpolatedPointMetrics=>equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr
    dependentInterpolatedPointMetrics=>equations%interpolation%dependentInterpPointMetrics(dependentVarType)%ptr

    IF(equationsSet%specification(3)==EQUATIONS_SET_CONSTITUTIVE_AND_GROWTH_LAW_IN_CELLML_SUBTYPE.OR. &
      equationsSet%specification(3)==EQUATIONS_SET_MR_AND_GROWTH_LAW_IN_CELLML_SUBTYPE) THEN
      CALL Field_ParameterSetGetLocalGaussPoint(dependentField,FIELD_U3_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
        & gaussPointNumber,localElementNumber,1,growthValues(1),err,error,*999)
      IF(numberofDimensions>1) THEN
        CALL Field_ParameterSetGetLocalGaussPoint(dependentField,FIELD_U3_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & gaussPointNumber,localElementNumber,2,growthValues(2),err,error,*999)
        IF(numberOfDimensions>2) THEN
          CALL Field_ParameterSetGetLocalGaussPoint(dependentField,FIELD_U3_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & gaussPointNumber,localElementNumber,3,growthValues(3),err,error,*999)
        ENDIF
      ENDIF
    ELSE
      growthValues = [1.0_DP,1.0_DP,1.0_DP]
    ENDIF

    !Interpolate fields at xi position
    CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointNumber,dependentInterpolatedPoint, &
      & err,error,*999)
    CALL Field_interpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointNumber,geometricInterpolatedPoint, &
      & err,error,*999)
    IF(ASSOCIATED(fibreInterpolatedPoint)) &
      & CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointNumber,fibreInterpolatedPoint, &
      & err,error,*999)

    !Calculate field metrics
    CALL Field_InterpolatedPointMetricsCalculate(numberOfXi,geometricInterpolatedPointMetrics,err,error,*999)
    CALL Field_InterpolatedPointMetricsCalculate(numberOfXi,dependentInterpolatedPointMetrics,err,error,*999)

    !Calculate F=dZ/dNU, the deformation gradient tensor at the xi location
    CALL FiniteElasticity_GaussDeformationGradientTensor(dependentInterpolatedPointMetrics, &
      & geometricInterpolatedPointMetrics,fibreInterpolatedPoint,dZdNu,err,error,*999)

    CALL FiniteElasticity_GaussGrowthTensor(equationsSet,numberOfDimensions,dZdNu,growthValues,Fg,Fe,Jg,Je, &
      & err,error,*999)
    
    IF(tensorEvaluateType==EQUATIONS_SET_R_CAUCHY_GREEN_DEFORMATION_TENSOR .OR. &
      & tensorEvaluateType==EQUATIONS_SET_GREEN_LAGRANGE_STRAIN_TENSOR) THEN
      CALL MatrixTransposeProduct(Fe(1:numberOfDimensions,1:numberOfXi),Fe(1:numberOfDimensions,1:numberOfXi), &
        & C(1:numberOfDimensions,1:numberOfDimensions),err,error,*999)
    ENDIF

    IF(tensorEvaluateType==EQUATIONS_SET_GREEN_LAGRANGE_STRAIN_TENSOR) THEN
      !Calculate E
      E(1:numberOfDimensions,1:numberOfDimensions)=0.5_DP*C(1:numberOfDimensions,1:numberOfDimensions)
      DO i=1,numberOfDimensions
        E(i,i)=E(i,i)-0.5_DP
      ENDDO !i
    ENDIF

    IF(tensorEvaluateType==EQUATIONS_SET_CAUCHY_STRESS_TENSOR) THEN
      !Get the interpolation parameters for this element
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,localElementNumber, &
        & equations%interpolation%materialsInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
      IF(ASSOCIATED(equations%interpolation%independentInterpParameters)) THEN
        CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,localElementNumber, &
          & equations%interpolation%independentInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
      END IF

      !Get interpolated points
      materialsInterpolatedPoint=>equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr
      IF(ASSOCIATED(equations%interpolation%independentInterpPoint)) &
        & independentInterpolatedPoint=>equations%interpolation%independentInterpPoint(dependentVarType)%ptr

      !Interpolate fields at Gauss point
      CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointNumber,materialsInterpolatedPoint, &
        & err,error,*999)
      IF(ASSOCIATED(independentInterpolatedPoint)) &
        & CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointNumber, &
        & independentInterpolatedPoint,err,error,*999)

!!\TODO the whole stress thing needs to be looked at as the routines below do not take in the deformation gradient that
!! is calculated above but rather they calculate it internally. This will lead to mismatches as things like growth are
!! not taken into account. 
      
      SELECT CASE(equationsSet%specification(3))
      CASE(EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE, &
        & EQUATIONS_SET_MR_AND_GROWTH_LAW_IN_CELLML_SUBTYPE)
        !Calculate the Cauchy stress tensor (in Voigt form) at the gauss point.
        Jznu=dependentInterpolatedPointMetrics%JACOBIAN/geometricInterpolatedPointMetrics%JACOBIAN
        ! Note that some problems, e.g. active contraction, require additonal fields to be evaluated at Gauss points. This is
        ! currently achieved by providing the gausspoint number to the FINITE_ELASTICITY_GAUSS_STRESS_TENSOR routine.
        ! However, the current  routine, FiniteElasticity_TensorInterpolateXi, aims to evaluate tensors as any xi, so the Gauss
        ! point number has been set to 0, which will generate an error for such problems.
        ! To address such issues, the FINITE_ELASTICITY_GAUSS_STRESS_TENSOR routine needs to be generalized to allow calculation
        ! of stress at any xi position and the GaussPoint number argument needs to be replace with a set of xi coordinates.
        CALL FINITE_ELASTICITY_GAUSS_STRESS_TENSOR(equationsSet,dependentInterpolatedPoint, &
          & materialsInterpolatedPoint,geometricInterpolatedPoint,cauchyStressVoigt,dZdNu,Jznu, &
          & localElementNumber,0,ERR,ERROR,*999)
        
        !Convert from Voigt form to tensor form. \TODO needs to be generalised for 2D
        DO nh=1,numberOfDimensions
          DO mh=1,numberOfDimensions
            cauchyStressTensor(mh,nh)=cauchyStressVoigt(TENSOR_TO_VOIGT3(mh,nh))
          ENDDO
        ENDDO
      CASE(EQUATIONS_SET_ORTHOTROPIC_MATERIAL_COSTA_SUBTYPE, EQUATIONS_SET_TRANSVERSE_ISOTROPIC_GUCCIONE_SUBTYPE)
        CALL FINITE_ELASTICITY_GAUSS_CAUCHY_TENSOR(equationsSet,dependentInterpolatedPoint,materialsInterpolatedPoint, &
          & geometricInterpolatedPoint,darcyInterpolatedPoint,independentInterpolatedPoint, &
          & cauchyStressTensor,Jznu,dZdNu,localElementNumber,0,ERR,ERROR,*999)
      CASE DEFAULT
        CALL FlagError("Not implemented ",err,error,*999)
      END SELECT
    END IF

    SELECT CASE(tensorEvaluateType)
    CASE(EQUATIONS_SET_DEFORMATION_GRADIENT_TENSOR)
      values(1:numberOfDimensions,1:numberOfDimensions)=dZdNu(1:numberOfDimensions,1:numberOfDimensions)
    CASE(EQUATIONS_SET_R_CAUCHY_GREEN_DEFORMATION_TENSOR)
      values(1:numberOfDimensions,1:numberOfDimensions)=C(1:numberOfDimensions,1:numberOfDimensions)
    CASE(EQUATIONS_SET_GREEN_LAGRANGE_STRAIN_TENSOR)
      values(1:numberOfDimensions,1:numberOfDimensions)=E(1:numberOfDimensions,1:numberOfDimensions)
    CASE(EQUATIONS_SET_CAUCHY_STRESS_TENSOR)
      values(1:numberOfDimensions,1:numberOfDimensions)=cauchyStressTensor(1:numberOfDimensions,1:numberOfDimensions)
    CASE(EQUATIONS_SET_SECOND_PK_STRESS_TENSOR)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      CALL FlagError("The tensor evalaute type of "//TRIM(NumberToVString(tensorEvaluateType,"*",err,error))//" is invalid "// &
        & "for finite elasticity equation sets.",err,error,*999)
    END SELECT
 
    EXITS("FiniteElasticity_TensorInterpolateGaussPoint")
    RETURN
999 ERRORS("FiniteElasticity_TensorInterpolateGaussPoint",err,error)
    EXITS("FiniteElasticity_TensorInterpolateGaussPoint")
    RETURN 1
    
  END SUBROUTINE FiniteElasticity_TensorInterpolateGaussPoint

  !
  !================================================================================================================================
  !

  !>Evaluates a tensor at a given element xi location. \TODO merge this with interpolate Gauss above.
  SUBROUTINE FiniteElasticity_TensorInterpolateXi(equationsSet,tensorEvaluateType,userElementNumber,xi,values,err,error,*)
    ! Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER, INTENT(IN) :: equationsSet !<A pointer to the equations set to calculate the tensor for
    INTEGER(INTG), INTENT(IN) :: tensorEvaluateType !<The type of tensor to evaluate.
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user element number to evaluate the tensor for
    REAL(DP), INTENT(IN) :: xi(:) !<The xi location to evaluate the tensor for.
    REAL(DP), INTENT(OUT) :: values(:,:) !<On exit, the interpolated tensor values.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    ! Local variables
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FIELD_TYPE), POINTER :: dependentField
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: geometricInterpolatedPoint, &
      & fibreInterpolatedPoint,dependentInterpolatedPoint,materialsInterpolatedPoint, &
      & independentInterpolatedPoint,darcyInterpolatedPoint
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: geometricInterpolatedPointMetrics, &
      & dependentInterpolatedPointMetrics
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: decompositionTopology
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: domainTopology
    TYPE(BASIS_TYPE), POINTER :: elementBasis
    LOGICAL :: userElementExists,ghostElement
    INTEGER(INTG) :: dependentVarType,meshComponentNumber
    INTEGER(INTG) :: numberOfDimensions,numberOfXi
    INTEGER(INTG) :: localElementNumber,i,nh,mh
    REAL(DP) :: dZdNu(3,3),dZdNuT(3,3),dZdXi(3,3),AZL(3,3),E(3,3),cauchyStressTensor(3,3),cauchyStressVoigt(6),Jznu

    ENTERS("FiniteElasticity_TensorInterpolateXi",err,error,*999)

    NULLIFY(geometricInterpolatedPoint)
    NULLIFY(fibreInterpolatedPoint)
    NULLIFY(dependentInterpolatedPoint)
    NULLIFY(materialsInterpolatedPoint)
    NULLIFY(independentInterpolatedPoint)
    NULLIFY(darcyInterpolatedPoint)
    NULLIFY(decomposition)
    NULLIFY(decompositionTopology)
    NULLIFY(domainTopology)
    NULLIFY(elementBasis)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)

    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMapping)
    CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
    NULLIFY(nonlinearMapping)
    CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)

    dependentVarType=nonlinearMapping%residualVariables(1)%ptr%VARIABLE_TYPE

    IF(.NOT.ASSOCIATED(equations%interpolation)) THEN
      CALL FlagError("Equations interpolation is not associated.",err,error,*999)
    END IF
    dependentField=>equations%interpolation%dependentField
    IF(.NOT.ASSOCIATED(dependentField)) THEN
      CALL FlagError("Equations dependent field is not associated.",err,error,*999)
    END IF
    NULLIFY(decomposition)
    CALL Field_DecompositionGet(dependentField,decomposition,err,error,*999)
    decomposition=>dependentField%decomposition
    CALL DECOMPOSITION_MESH_COMPONENT_NUMBER_GET(decomposition,meshComponentNumber,err,error,*999)
    decompositionTopology=>decomposition%topology
    domainTopology=>decomposition%domain(meshComponentNumber)%ptr%topology
    CALL DECOMPOSITION_TOPOLOGY_ELEMENT_CHECK_EXISTS(decompositionTopology,userElementNumber, &
      & userElementExists,localElementNumber,ghostElement,err,error,*999)
    IF(.NOT.userElementExists) THEN
      CALL FlagError("The specified user element number of "// &
        & TRIM(NumberToVstring(userElementNumber,"*",err,error))// &
        & " does not exist in the decomposition for the dependent field.",err,error,*999)
    END IF
    CALL DomainTopology_ElementBasisGet(domainTopology,userElementNumber,elementBasis,err,error,*999)

    !Get the interpolation parameters for this element
    CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,localElementNumber, &
      & equations%interpolation%geometricInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
    IF(ASSOCIATED(equations%interpolation%fibreInterpParameters)) THEN
      CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,localElementNumber, &
        & equations%interpolation%fibreInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
    END IF
    CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,localElementNumber, &
      & equations%interpolation%dependentInterpParameters(dependentVarType)%ptr,err,error,*999)

    !Get interpolated points
    geometricInterpolatedPoint=>equations%interpolation%geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr
    IF(ASSOCIATED(equations%interpolation%fibreInterpPoint)) THEN
      fibreInterpolatedPoint=>equations%interpolation%fibreInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr
    END IF
    dependentInterpolatedPoint=>equations%interpolation%dependentInterpPoint(dependentVarType)%ptr

    !Get interpolated point metrics
    geometricInterpolatedPointMetrics=>equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr
    dependentInterpolatedPointMetrics=>equations%interpolation%dependentInterpPointMetrics(dependentVarType)%ptr

    !Interpolate fields at xi position
    CALL FIELD_INTERPOLATE_XI(FIRST_PART_DERIV,xi,dependentInterpolatedPoint,err,error,*999)
    CALL FIELD_INTERPOLATE_XI(FIRST_PART_DERIV,xi,geometricInterpolatedPoint,err,error,*999)
    IF(ASSOCIATED(fibreInterpolatedPoint)) THEN
      CALL FIELD_INTERPOLATE_XI(FIRST_PART_DERIV,xi,fibreInterpolatedPoint,err,error,*999)
    END IF

    !Calculate field metrics
    CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE( &
      & elementBasis%number_of_xi,geometricInterpolatedPointMetrics,err,error,*999)
    CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE( &
      & elementBasis%number_of_xi,dependentInterpolatedPointMetrics,err,error,*999)

    !Calculate F=dZ/dNU, the deformation gradient tensor at the xi location
    numberOfDimensions=equationsSet%region%coordinate_system%number_of_dimensions
    numberOfXi=elementBasis%number_of_xi
    CALL FiniteElasticity_GaussDeformationGradientTensor(dependentInterpolatedPointMetrics, &
      & geometricInterpolatedPointMetrics,fibreInterpolatedPoint,dZdNu,err,error,*999)

    IF(tensorEvaluateType==EQUATIONS_SET_R_CAUCHY_GREEN_DEFORMATION_TENSOR .OR. &
      & tensorEvaluateType==EQUATIONS_SET_GREEN_LAGRANGE_STRAIN_TENSOR) THEN
      CALL MatrixTranspose(dZdNu,dZdNuT,err,error,*999)
      CALL MatrixProduct(dZdNuT,dZdNu,AZL,err,error,*999)
     END IF

    IF(tensorEvaluateType==EQUATIONS_SET_GREEN_LAGRANGE_STRAIN_TENSOR) THEN
      !Calculate E
      E=0.5_DP*AZL
      DO i=1,3
        E(i,i)=E(i,i)-0.5_DP
      END DO
    END IF

    IF(tensorEvaluateType==EQUATIONS_SET_CAUCHY_STRESS_TENSOR) THEN
      !Get the interpolation parameters for this element
      CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,localElementNumber, &
        & equations%interpolation%materialsInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
      IF(ASSOCIATED(equations%interpolation%independentInterpParameters)) THEN
        CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,localElementNumber, &
          & equations%interpolation%independentInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
      END IF

      !Get interpolated points
      materialsInterpolatedPoint=>equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr
      IF(ASSOCIATED(equations%interpolation%independentInterpPoint)) THEN
        independentInterpolatedPoint=>equations%interpolation%independentInterpPoint(dependentVarType)%ptr
      END IF

      !Interpolate fields at xi position
      CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,xi,materialsInterpolatedPoint,err,error,*999)
      IF(ASSOCIATED(independentInterpolatedPoint)) THEN
        CALL FIELD_INTERPOLATE_XI(FIRST_PART_DERIV,xi,independentInterpolatedPoint,err,error,*999)
      END IF

      SELECT CASE(equationsSet%specification(3))
      CASE(EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE, &
        & EQUATIONS_SET_MR_AND_GROWTH_LAW_IN_CELLML_SUBTYPE)
        !Calculate the Cauchy stress tensor (in Voigt form) at the gauss point.
        Jznu=dependentInterpolatedPointMetrics%JACOBIAN/geometricInterpolatedPointMetrics%JACOBIAN
        ! Note that some problems, e.g. active contraction, require additonal fields to be evaluated at Gauss points. This is
        ! currently achieved by providing the gausspoint number to the FINITE_ELASTICITY_GAUSS_STRESS_TENSOR routine.
        ! However, the current  routine, FiniteElasticity_TensorInterpolateXi, aims to evaluate tensors as any xi, so the Gauss
        ! point number has been set to 0, which will generate an error for such problems.
        ! To address such issues, the FINITE_ELASTICITY_GAUSS_STRESS_TENSOR routine needs to be generalized to allow calculation
        ! of stress at any xi position and the GaussPoint number argument needs to be replace with a set of xi coordinates.
        CALL FINITE_ELASTICITY_GAUSS_STRESS_TENSOR(equationsSet,dependentInterpolatedPoint, &
          & materialsInterpolatedPoint,geometricInterpolatedPoint,cauchyStressVoigt,dZdNu,Jznu, &
          & localElementNumber,0,ERR,ERROR,*999)
        
        !Convert from Voigt form to tensor form.
        DO nh=1,3
          DO mh=1,3
            cauchyStressTensor(mh,nh)=cauchyStressVoigt(TENSOR_TO_VOIGT3(mh,nh))
          ENDDO
        ENDDO
      CASE(EQUATIONS_SET_ORTHOTROPIC_MATERIAL_COSTA_SUBTYPE, EQUATIONS_SET_TRANSVERSE_ISOTROPIC_GUCCIONE_SUBTYPE)
        CALL FINITE_ELASTICITY_GAUSS_CAUCHY_TENSOR(equationsSet,dependentInterpolatedPoint,materialsInterpolatedPoint, &
          & geometricInterpolatedPoint,darcyInterpolatedPoint,independentInterpolatedPoint, &
          & cauchyStressTensor,Jznu,dZdNu,localElementNumber,0,ERR,ERROR,*999)
      CASE DEFAULT
        CALL FlagError("Not implemented ",err,error,*999)
      END SELECT
    END IF

    SELECT CASE(tensorEvaluateType)
    CASE(EQUATIONS_SET_DEFORMATION_GRADIENT_TENSOR)
      values=dZdNu
    CASE(EQUATIONS_SET_R_CAUCHY_GREEN_DEFORMATION_TENSOR)
      values=AZL
    CASE(EQUATIONS_SET_GREEN_LAGRANGE_STRAIN_TENSOR)
      values=E
    CASE(EQUATIONS_SET_CAUCHY_STRESS_TENSOR)
      values=cauchyStressTensor
    CASE(EQUATIONS_SET_SECOND_PK_STRESS_TENSOR)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      CALL FlagError("The tensor evalaute type of "//TRIM(NumberToVString(tensorEvaluateType,"*",err,error))//" is invalid "// &
        & "for finite elasticity equation sets",err,error,*999)
    END SELECT

    EXITS("FiniteElasticity_TensorInterpolateXi")
    RETURN
999 ERRORSEXITS("FiniteElasticity_TensorInterpolateXi",err,error)
    RETURN 1
  END SUBROUTINE FiniteElasticity_TensorInterpolateXi

  !
  !================================================================================================================================
  !

  !Evaluates the Jacobian surface traction (pressure) term of the equilibrium equation. Here it is assumed that pressure is constant
  !(if not: the jacobian has to be extended to include this) and that along the boundary of the boundary faces (the boundary line)
  !minimal one direction perpendicular to that boundary line is fixed, or that we have no boundary line at all (a closed body). In
  !these cases the jacobian is symmetrical. See Rumpel & Schweizerhof, "Hydrostatic fluid loading in non-linear finite element
  !analysis".
  SUBROUTINE FiniteElasticity_SurfacePressureJacobianEvaluate(EQUATIONS_SET,ELEMENT_NUMBER,err,error,*)
    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS
    TYPE(BASIS_PTR_TYPE) :: BASES(3)
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DECOMPOSITION_ELEMENT_TYPE), POINTER :: ELEMENT
    TYPE(DECOMPOSITION_FACE_TYPE), POINTER :: FACE
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsJacobianType), POINTER :: jacobianMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: DEPENDENT_INTERPOLATION_PARAMETERS
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: PRESSURE_INTERPOLATION_PARAMETERS
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: dependentInterpPoint,PRESSURE_INTERP_POINT
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: dependentInterpPointMetrics
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: DEPENDENT_QUADRATURE_SCHEME
    TYPE(QUADRATURE_SCHEME_PTR_TYPE) :: QUADRATURE_SCHEMES(3)
    INTEGER(INTG) :: FACE_NUMBER,xiDirection(3),orientation
    INTEGER(INTG) :: FIELD_VAR_U_TYPE,FIELD_VAR_DELUDELN_TYPE,MESH_COMPONENT_NUMBER
    INTEGER(INTG) :: oh,mh,ms,mhs,nh,ns,nhs,ng,naf
    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS,NUMBER_OF_LOCAL_FACES
    INTEGER(INTG) :: SUM_ELEMENT_PARAMETERS
    INTEGER(INTG) :: ELEMENT_BASE_DOF_INDEX(3),NUMBER_OF_FACE_PARAMETERS(3)
    INTEGER(INTG), PARAMETER :: OFF_DIAG_COMP(3)=[0,1,3],OFF_DIAG_DEP_VAR1(3)=[1,1,2],OFF_DIAG_DEP_VAR2(3)=[2,3,3]
    REAL(DP) :: PRESSURE_GAUSS,GW_PRESSURE
    REAL(DP) :: NORMAL(3),GW_PRESSURE_W(2),TEMP3, TEMP4
    REAL(DP) :: TEMPVEC1(2),TEMPVEC2(2),TEMPVEC3(3),TEMPVEC4(3),TEMPVEC5(3)
    LOGICAL :: NONZERO_PRESSURE

    ENTERS("FiniteElasticity_SurfacePressureJacobianEvaluate",err,error,*999)

    NULLIFY(DEPENDENT_BASIS)
    NULLIFY(DECOMPOSITION)
    NULLIFY(ELEMENT)
    NULLIFY(EQUATIONS,vectorMapping,vectorMatrices,nonlinearMapping,nonlinearMatrices,jacobianMatrix)
    NULLIFY(DEPENDENT_INTERPOLATION_PARAMETERS,PRESSURE_INTERPOLATION_PARAMETERS)
    NULLIFY(dependentInterpPoint,dependentInterpPointMetrics,PRESSURE_INTERP_POINT)
    NULLIFY(DEPENDENT_FIELD)
    NULLIFY(FIELD_VARIABLE)
    NULLIFY(DEPENDENT_QUADRATURE_SCHEME)

    NUMBER_OF_DIMENSIONS=EQUATIONS_SET%REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS

    EQUATIONS=>EQUATIONS_SET%EQUATIONS
    vectorEquations=>equations%vectorEquations
    vectorMatrices=>vectorEquations%vectorMatrices
    nonlinearMatrices=>vectorMatrices%nonlinearMatrices
    jacobianMatrix=>nonlinearMatrices%jacobians(1)%ptr

    DEPENDENT_FIELD=>equations%interpolation%dependentField
    DECOMPOSITION=>DEPENDENT_FIELD%DECOMPOSITION
    MESH_COMPONENT_NUMBER=DECOMPOSITION%MESH_COMPONENT_NUMBER
    ELEMENT=>DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)
    NUMBER_OF_LOCAL_FACES=DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%ptr% &
      & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS%NUMBER_OF_LOCAL_FACES
    
    FIELD_VARIABLE=>vectorEquations%vectorMapping%nonlinearMapping%residualVariables(1)%ptr
    FIELD_VAR_U_TYPE=vectorEquations%vectorMapping%nonlinearMapping%residualVariables(1)%ptr%VARIABLE_TYPE
    FIELD_VAR_DELUDELN_TYPE=vectorEquations%vectorMapping%rhsMapping%rhsVariableType

    !Surface pressure term calculation: Loop over all faces
    DO naf=1,NUMBER_OF_LOCAL_FACES
      FACE_NUMBER=ELEMENT%ELEMENT_FACES(naf)
      FACE=>DECOMPOSITION%TOPOLOGY%FACES%FACES(FACE_NUMBER)

      !Check if it's a boundary face
      IF(FACE%BOUNDARY_FACE) THEN
        xiDirection(3)=ABS(FACE%XI_DIRECTION)

        PRESSURE_INTERPOLATION_PARAMETERS=>equations%interpolation%dependentInterpParameters(FIELD_VAR_DELUDELN_TYPE)%ptr
        CALL FIELD_INTERPOLATION_PARAMETERS_FACE_GET(FIELD_PRESSURE_VALUES_SET_TYPE,FACE_NUMBER, &
          & PRESSURE_INTERPOLATION_PARAMETERS,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
        PRESSURE_INTERP_POINT=>equations%interpolation%dependentInterpPoint(FIELD_VAR_DELUDELN_TYPE)%ptr

        !Check if nonzero surface pressure is defined on the face
        NONZERO_PRESSURE=ANY(ABS(PRESSURE_INTERPOLATION_PARAMETERS%PARAMETERS(:,xiDirection(3)))>ZERO_TOLERANCE)

        !Nonzero surface pressure found?
        IF(NONZERO_PRESSURE) THEN
          MESH_COMPONENT_NUMBER=DECOMPOSITION%MESH_COMPONENT_NUMBER
          DEPENDENT_BASIS=>DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%ptr%TOPOLOGY%FACES%FACES(FACE_NUMBER)%BASIS
          DEPENDENT_QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr

          DEPENDENT_INTERPOLATION_PARAMETERS=>equations%interpolation%dependentInterpParameters(FIELD_VAR_U_TYPE)%ptr
          CALL FIELD_INTERPOLATION_PARAMETERS_FACE_GET(FIELD_VALUES_SET_TYPE,FACE_NUMBER, &
            & DEPENDENT_INTERPOLATION_PARAMETERS,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
          dependentInterpPoint=>equations%interpolation%dependentInterpPoint(FIELD_VAR_U_TYPE)%ptr
          dependentInterpPointMetrics=>equations%interpolation%dependentInterpPointMetrics(FIELD_VAR_U_TYPE)%ptr

          SUM_ELEMENT_PARAMETERS=0
          !Loop over geometric dependent basis functions.
          DO nh=1,NUMBER_OF_DIMENSIONS
            MESH_COMPONENT_NUMBER=FIELD_VARIABLE%COMPONENTS(nh)%MESH_COMPONENT_NUMBER
            DEPENDENT_BASIS=>DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%ptr%TOPOLOGY%FACES%FACES(FACE_NUMBER)%BASIS
            BASES(nh)%ptr=>DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%ptr% &
              & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
            QUADRATURE_SCHEMES(nh)%ptr=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
            NUMBER_OF_FACE_PARAMETERS(nh)=DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
            ELEMENT_BASE_DOF_INDEX(nh)=SUM_ELEMENT_PARAMETERS
            SUM_ELEMENT_PARAMETERS=SUM_ELEMENT_PARAMETERS+BASES(nh)%ptr%NUMBER_OF_ELEMENT_PARAMETERS
          ENDDO !nh

          xiDirection(1)=OTHER_XI_DIRECTIONS3(xiDirection(3),2,1)
          xiDirection(2)=OTHER_XI_DIRECTIONS3(xiDirection(3),3,1)
          orientation=SIGN(1,OTHER_XI_ORIENTATIONS3(xiDirection(1),xiDirection(2))*FACE%XI_DIRECTION)

          !Loop over all Gauss points
          DO ng=1,DEPENDENT_QUADRATURE_SCHEME%NUMBER_OF_GAUSS
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng, &
              & PRESSURE_INTERP_POINT,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng, &
              & dependentInterpPoint,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(COORDINATE_JACOBIAN_AREA_TYPE, &
              & dependentInterpPointMetrics,err,error,*999)

            CALL CrossProduct(dependentInterpPointMetrics%DX_DXI(:,1), &
              & dependentInterpPointMetrics%DX_DXI(:,2),NORMAL,err,error,*999)
            PRESSURE_GAUSS=PRESSURE_INTERP_POINT%VALUES(xiDirection(3),NO_PART_DERIV)*orientation
            GW_PRESSURE=DEPENDENT_QUADRATURE_SCHEME%GAUSS_WEIGHTS(ng)*PRESSURE_GAUSS

            DO oh=1,OFF_DIAG_COMP(NUMBER_OF_DIMENSIONS)
              nh=OFF_DIAG_DEP_VAR1(oh)
              mh=OFF_DIAG_DEP_VAR2(oh)
              GW_PRESSURE_W(1:2)=(NORMAL(mh)*dependentInterpPointMetrics%DXI_DX(1:2,nh)- &
                & dependentInterpPointMetrics%DXI_DX(1:2,mh)*NORMAL(nh))*GW_PRESSURE
              DO ns=1,NUMBER_OF_FACE_PARAMETERS(nh)
                !Loop over element rows belonging to geometric dependent variables
                nhs=ELEMENT_BASE_DOF_INDEX(nh)+ &
                  & BASES(nh)%ptr%ELEMENT_PARAMETERS_IN_LOCAL_FACE(ns,naf)
                TEMPVEC1(1:2)=GW_PRESSURE_W(1:2)*QUADRATURE_SCHEMES(nh)%ptr% &
                  & GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(1:2),ng)
                DO ms=1,NUMBER_OF_FACE_PARAMETERS(mh)
                  mhs=ELEMENT_BASE_DOF_INDEX(mh)+ &
                    & BASES(mh)%PTR%ELEMENT_PARAMETERS_IN_LOCAL_FACE(ms,naf)
                  TEMPVEC2=QUADRATURE_SCHEMES(mh)%PTR%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                  jacobianMatrix%elementJacobian%matrix(mhs,nhs)=jacobianMatrix%elementJacobian%matrix(mhs,nhs)+ &
                    & DOT_PRODUCT(TEMPVEC1,TEMPVEC2)* &
                    & EQUATIONS%INTERPOLATION%dependentInterpParameters(FIELD_VAR_U_TYPE)%PTR%SCALE_FACTORS(ms,mh)* &
                    & EQUATIONS%INTERPOLATION%dependentInterpParameters(FIELD_VAR_U_TYPE)%PTR%SCALE_FACTORS(ns,nh)
                ENDDO !ms
              ENDDO !ns
            ENDDO !oh

            DO oh=1,OFF_DIAG_COMP(NUMBER_OF_DIMENSIONS)
              nh=OFF_DIAG_DEP_VAR1(oh)
              mh=OFF_DIAG_DEP_VAR2(oh)
              GW_PRESSURE_W(1:2)=(NORMAL(nh)*dependentInterpPointMetrics%DXI_DX(1:2,mh)- &
                & dependentInterpPointMetrics%DXI_DX(1:2,nh)*NORMAL(mh))*GW_PRESSURE
              DO ms=1,NUMBER_OF_FACE_PARAMETERS(mh)
                !Loop over element rows belonging to geometric dependent variables
                mhs=ELEMENT_BASE_DOF_INDEX(mh)+ &
                  & BASES(mh)%PTR%ELEMENT_PARAMETERS_IN_LOCAL_FACE(ms,naf)
                TEMPVEC1(1:2)=GW_PRESSURE_W(1:2)*QUADRATURE_SCHEMES(mh)%PTR% &
                  & GAUSS_BASIS_FNS(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(1:2),ng)
                DO ns=1,NUMBER_OF_FACE_PARAMETERS(nh)
                  nhs=ELEMENT_BASE_DOF_INDEX(nh)+ &
                    & BASES(nh)%PTR%ELEMENT_PARAMETERS_IN_LOCAL_FACE(ns,naf)
                  TEMPVEC2=QUADRATURE_SCHEMES(nh)%PTR%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)
                  jacobianMatrix%elementJacobian%matrix(nhs,mhs)=jacobianMatrix%elementJacobian%matrix(nhs,mhs)+ &
                    & DOT_PRODUCT(TEMPVEC1,TEMPVEC2)* &
                    & EQUATIONS%INTERPOLATION%dependentInterpParameters(FIELD_VAR_U_TYPE)%PTR%SCALE_FACTORS(ms,mh)* &
                    & EQUATIONS%INTERPOLATION%dependentInterpParameters(FIELD_VAR_U_TYPE)%PTR%SCALE_FACTORS(ns,nh)
                ENDDO !ns
              ENDDO !ms
            ENDDO !oh
          ENDDO !ng
        ENDIF !Non-zero pressure on face
      ENDIF !Boundary face
    ENDDO !naf

    EXITS("FiniteElasticity_SurfacePressureJacobianEvaluate")
    RETURN
999 ERRORS("FiniteElasticity_SurfacePressureJacobianEvaluate",err,error)
    EXITS("FiniteElasticity_SurfacePressureJacobianEvaluate")
    RETURN 1

  END SUBROUTINE FiniteElasticity_SurfacePressureJacobianEvaluate

  !
  !================================================================================================================================
  !

  !Evaluates the surface traction (pressure) term of the equilibrium equation
  SUBROUTINE FiniteElasticity_SurfacePressureResidualEvaluate(EQUATIONS_SET,ELEMENT_NUMBER,var1,var2,err,error,*)
    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER
    INTEGER(INTG), INTENT(IN) :: var1 !<'U' variable number in single-physics case
    INTEGER(INTG), INTENT(IN) :: var2 !<'DELUDELN' variable number in single-physics case
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_FACE_BASIS,COMPONENT_FACE_BASIS,COMPONENT_BASIS
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DECOMPOSITION_ELEMENT_TYPE), POINTER :: DECOMP_ELEMENT
    TYPE(DECOMPOSITION_FACE_TYPE), POINTER :: DECOMP_FACE
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: FACE_DEPENDENT_INTERPOLATION_PARAMETERS, &
      & DEPENDENT_INTERPOLATION_PARAMETERS
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: FACE_PRESSURE_INTERPOLATION_PARAMETERS
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: FACE_DEPENDENT_INTERPOLATED_POINT
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: FACE_DEPENDENT_INTERPOLATED_POINT_METRICS
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: FACE_PRESSURE_INTERPOLATED_POINT
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: FACE_QUADRATURE_SCHEME,COMPONENT_FACE_QUADRATURE_SCHEME
    INTEGER(INTG) :: FIELD_VAR_U_TYPE,FIELD_VAR_DUDN_TYPE,MESH_COMPONENT_NUMBER
    INTEGER(INTG) :: element_face_idx,face_number,gauss_idx
    INTEGER(INTG) :: component_idx,element_base_dof_idx,element_dof_idx,parameter_idx,face_parameter_idx
    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS,NUMBER_OF_LOCAL_FACES
    INTEGER(INTG) :: xiDirection(3),orientation
    REAL(DP) :: PRESSURE_GAUSS,GW_PRESSURE,GW_PRESSURE_NORMAL_COMPONENT
    REAL(DP) :: NORMAL(3)
    LOGICAL :: NONZERO_PRESSURE
    INTEGER(INTG) :: EQUATIONS_SET_SUBTYPE

    ENTERS("FiniteElasticity_SurfacePressureResidualEvaluate",err,error,*999)

    NULLIFY(DEPENDENT_FACE_BASIS,COMPONENT_FACE_BASIS,COMPONENT_BASIS)
    NULLIFY(DECOMPOSITION)
    NULLIFY(DECOMP_ELEMENT)
    NULLIFY(DECOMP_FACE)
    NULLIFY(EQUATIONS)
    NULLIFY(EQUATIONS,nonlinearMatrices)
    NULLIFY(DEPENDENT_FIELD,FIELD_VARIABLE)
    NULLIFY(FACE_DEPENDENT_INTERPOLATION_PARAMETERS)
    NULLIFY(FACE_DEPENDENT_INTERPOLATED_POINT,FACE_DEPENDENT_INTERPOLATED_POINT_METRICS)
    NULLIFY(FACE_PRESSURE_INTERPOLATION_PARAMETERS,FACE_PRESSURE_INTERPOLATED_POINT)
    NULLIFY(COMPONENT_FACE_QUADRATURE_SCHEME,FACE_QUADRATURE_SCHEME)

    NUMBER_OF_DIMENSIONS=EQUATIONS_SET%REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS

    !Grab pointers of interest
    EQUATIONS=>EQUATIONS_SET%EQUATIONS
    vectorEquations=>equations%vectorEquations
    nonlinearMatrices=>vectorEquations%vectorMatrices%nonlinearMatrices
    EQUATIONS_SET_SUBTYPE = EQUATIONS_SET%SPECIFICATION(3)

    IF (EQUATIONS_SET_SUBTYPE == EQUATIONS_SET_REFERENCE_STATE_TRANSVERSE_GUCCIONE_SUBTYPE) THEN
      DEPENDENT_FIELD=>equations%interpolation%geometricField
    ELSE
      DEPENDENT_FIELD=>equations%interpolation%dependentField
    END IF
    DECOMPOSITION=>DEPENDENT_FIELD%DECOMPOSITION
    MESH_COMPONENT_NUMBER=DECOMPOSITION%MESH_COMPONENT_NUMBER
    DECOMP_ELEMENT=>DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)

    !Interpolation parameter for metric tensor
    FIELD_VARIABLE=>vectorEquations%vectorMapping%nonlinearMapping%residualVariables(1)%ptr
    FIELD_VAR_U_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
    FIELD_VAR_DUDN_TYPE=vectorEquations%vectorMapping%rhsMapping%rhsVariableType
    NUMBER_OF_LOCAL_FACES=DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%ptr%TOPOLOGY%ELEMENTS% &
      & ELEMENTS(ELEMENT_NUMBER)%BASIS%NUMBER_OF_LOCAL_FACES

    !Surface pressure term calculation: Loop over all faces
    DO element_face_idx=1,NUMBER_OF_LOCAL_FACES
      face_number=DECOMP_ELEMENT%ELEMENT_FACES(element_face_idx)
      DECOMP_FACE=>DECOMPOSITION%TOPOLOGY%FACES%FACES(face_number)

      !Check if it's a boundary face
      IF(DECOMP_FACE%BOUNDARY_FACE) THEN !!temporary until MESH_FACE (or equivalent) is available (decomp face includes ghost faces?)
        xiDirection(3)=ABS(DECOMP_FACE%XI_DIRECTION)  ! if xi=0, this can be a negative number
        !Get pressure interpolation objects (DELUDELN pressure_values_set_type)
        FACE_PRESSURE_INTERPOLATION_PARAMETERS=>equations%interpolation%dependentInterpParameters(FIELD_VAR_DUDN_TYPE)%ptr
        CALL FIELD_INTERPOLATION_PARAMETERS_FACE_GET(FIELD_PRESSURE_VALUES_SET_TYPE,face_number, &
          & FACE_PRESSURE_INTERPOLATION_PARAMETERS,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
        FACE_PRESSURE_INTERPOLATED_POINT=>equations%interpolation%dependentInterpPoint(var2)%ptr

        !Check if nonzero surface pressure is defined on the face
        NONZERO_PRESSURE=ANY(ABS(FACE_PRESSURE_INTERPOLATION_PARAMETERS%PARAMETERS(:,xiDirection(3)))>ZERO_TOLERANCE)

        !Nonzero surface pressure found?
        IF(NONZERO_PRESSURE) THEN
          !Grab some other pointers
          DEPENDENT_FACE_BASIS=>DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%ptr%TOPOLOGY%FACES%FACES(face_number)%BASIS
          FACE_QUADRATURE_SCHEME=>DEPENDENT_FACE_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
          
          IF (EQUATIONS_SET_SUBTYPE == EQUATIONS_SET_REFERENCE_STATE_TRANSVERSE_GUCCIONE_SUBTYPE) THEN
            FACE_DEPENDENT_INTERPOLATION_PARAMETERS=>equations%interpolation%geometricInterpParameters(FIELD_VAR_U_TYPE)%ptr
            FACE_DEPENDENT_INTERPOLATED_POINT=>equations%interpolation%geometricInterpPoint(FIELD_VAR_U_TYPE)%ptr
            FACE_DEPENDENT_INTERPOLATED_POINT_METRICS=>equations%interpolation% &
              & geometricInterpPointMetrics(FIELD_VAR_U_TYPE)%ptr
          ELSE
            FACE_DEPENDENT_INTERPOLATION_PARAMETERS=>equations%interpolation%dependentInterpParameters(FIELD_VAR_U_TYPE)%ptr
            FACE_DEPENDENT_INTERPOLATED_POINT=>equations%interpolation%dependentInterpPoint(FIELD_VAR_U_TYPE)%ptr
            FACE_DEPENDENT_INTERPOLATED_POINT_METRICS=>equations%interpolation% &
              & dependentInterpPointMetrics(FIELD_VAR_U_TYPE)%ptr
          ENDIF
          CALL FIELD_INTERPOLATION_PARAMETERS_FACE_GET(FIELD_VALUES_SET_TYPE,FACE_NUMBER, &
            & FACE_DEPENDENT_INTERPOLATION_PARAMETERS,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)

          xiDirection(1)=OTHER_XI_DIRECTIONS3(xiDirection(3),2,1)
          xiDirection(2)=OTHER_XI_DIRECTIONS3(xiDirection(3),3,1)
          orientation=SIGN(1,OTHER_XI_ORIENTATIONS3(xiDirection(1),xiDirection(2))*DECOMP_FACE%XI_DIRECTION)

          !Start integrating
          ! Note: As the code will look for P(appl) in the *normal* component to the face, the
          !       initial assignment of P(appl) will have to be made appropriately during bc assignment
          DO gauss_idx=1,FACE_QUADRATURE_SCHEME%NUMBER_OF_GAUSS
            !Interpolate p(appl) at gauss point
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & FACE_PRESSURE_INTERPOLATED_POINT,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & FACE_DEPENDENT_INTERPOLATED_POINT,err,error,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(COORDINATE_JACOBIAN_AREA_TYPE, &
              & FACE_DEPENDENT_INTERPOLATED_POINT_METRICS,err,error,*999)
           
            CALL CrossProduct(FACE_DEPENDENT_INTERPOLATED_POINT_METRICS%DX_DXI(:,1), &
              & FACE_DEPENDENT_INTERPOLATED_POINT_METRICS%DX_DXI(:,2),NORMAL,err,error,*999)

            PRESSURE_GAUSS=FACE_PRESSURE_INTERPOLATED_POINT%VALUES(xiDirection(3),NO_PART_DERIV)*orientation
            GW_PRESSURE=FACE_QUADRATURE_SCHEME%GAUSS_WEIGHTS(gauss_idx)*PRESSURE_GAUSS
            element_base_dof_idx=0
            !Loop over 3 components
            DO component_idx=1,NUMBER_OF_DIMENSIONS
              MESH_COMPONENT_NUMBER=FIELD_VARIABLE%COMPONENTS(component_idx)%MESH_COMPONENT_NUMBER
              COMPONENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%ptr% &
                & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
              COMPONENT_FACE_BASIS=>DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%ptr%TOPOLOGY%FACES%FACES(face_number)%BASIS
              COMPONENT_FACE_QUADRATURE_SCHEME=>COMPONENT_FACE_BASIS% &
                & QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
              GW_PRESSURE_NORMAL_COMPONENT=GW_PRESSURE*NORMAL(component_idx)
              DO face_parameter_idx=1,COMPONENT_FACE_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                parameter_idx=COMPONENT_BASIS%ELEMENT_PARAMETERS_IN_LOCAL_FACE(face_parameter_idx,element_face_idx)
                element_dof_idx=element_base_dof_idx+parameter_idx
                nonlinearMatrices%elementResidual%vector(element_dof_idx)= &
                  & nonlinearMatrices%elementResidual%vector(element_dof_idx)+ & ! sign: double -'s. p(appl) always opposite to normal'
                  & GW_PRESSURE_NORMAL_COMPONENT* &
                  & COMPONENT_FACE_QUADRATURE_SCHEME%GAUSS_BASIS_FNS(face_parameter_idx,NO_PART_DERIV,gauss_idx)
              ENDDO !face_parameter_idx
              !Update element_base_dof_idx
              element_base_dof_idx=element_base_dof_idx+COMPONENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
            ENDDO !component_idx
          ENDDO !gauss_idx
        ENDIF !nonzero surface pressure check
      ENDIF !boundary face check
    ENDDO !element_face_idx

    EXITS("FiniteElasticity_SurfacePressureResidualEvaluate")
    RETURN
999 ERRORS("FiniteElasticity_SurfacePressureResidualEvaluate",err,error)
    EXITS("FiniteElasticity_SurfacePressureResidualEvaluate")
    RETURN 1

  END SUBROUTINE FiniteElasticity_SurfacePressureResidualEvaluate

  !
  !================================================================================================================================
  !

  !>Evaluates the deformation gradient tensor at a given Gauss point
  SUBROUTINE FiniteElasticity_GaussDeformationGradientTensor(dependentInterpPointMetrics,geometricInterpPointMetrics,&
    & fibreInterpolatedPoint,dZdNu,err,error,*)

    !Argument variables
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: dependentInterpPointMetrics,geometricInterpPointMetrics
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: fibreInterpolatedPoint
    REAL(DP), INTENT(OUT) :: dZdNu(3,3) !<dZdNu(coordinateIdx,coordianteIdx). On return, the deformation gradient tensor
    INTEGER(INTG), INTENT(OUT) :: err   !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfXDimensions,numberOfXiDimensions,numberOfZDimensions
    REAL(DP) ::detdZdX,detdZdNu,dNudX(3,3),dXdNu(3,3),dNuDXi(3,3),dXidNu(3,3)

    ENTERS("FiniteElasticity_GaussDeformationGradientTensor",err,error,*999)

    IF(ASSOCIATED(dependentInterpPointMetrics)) THEN
      IF(ASSOCIATED(geometricInterpPointMetrics)) THEN
        numberOfXDimensions=geometricInterpPointMetrics%NUMBER_OF_X_DIMENSIONS
        numberOfXiDimensions=geometricInterpPointMetrics%NUMBER_OF_XI_DIMENSIONS
        numberOfZDimensions=dependentInterpPointMetrics%NUMBER_OF_X_DIMENSIONS

        CALL Coordinates_MaterialSystemCalculate(geometricInterpPointMetrics,fibreInterpolatedPoint,dNudX,dXdNu, &
          & dNudXi(1:numberOfXDimensions,1:numberOfXiDimensions), &
          & dXidNu(1:numberOfXiDimensions,1:numberOfXDimensions),err,error,*999)
        !dZ/dNu = dZ/dXi * dXi/dNu  (deformation gradient tensor, F)
        CALL MatrixProduct(dependentInterpPointMetrics%DX_DXI(1:numberOfZDimensions,1:numberOfXiDimensions), &
          & dXiDNu(1:numberOfXiDimensions,1:numberOfXDimensions),dZdNu(1:numberOfZDimensions,1:numberOfXDimensions), &
          & err,error,*999)

        IF(numberOfZDimensions == 2) THEN
          dZdNu(:,3) = [0.0_DP,0.0_DP,1.0_DP]
          dZdNu(3,1:2) = 0.0_DP
        ENDIF

        IF(DIAGNOSTICS1) THEN
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"",err,error,*999)
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Calculated deformation gradient tensor:",err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of Z dimensions  = ",numberOfZDimensions,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of Xi dimensions = ",numberOfXiDimensions,err,error,*999)
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Derivative of X wrt to Nu coordinates:",err,error,*999)
          CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfXiDimensions,1,1,numberOfXDimensions, &
            & numberOfXDimensions,numberOfXDimensions,dXidNu,WRITE_STRING_MATRIX_NAME_AND_INDICES, &
            & '("    dX_dNu','(",I1,",:)','   :",3(X,E13.6))','(19X,3(X,E13.6))',err,error,*999)
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Deformation gradient tensor wrt Nu coordinates:",err,error,*999)
          CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfZDimensions,1,1,numberOfXDimensions, &
            & numberOfXDimensions,numberOfXDimensions,dZdNu,WRITE_STRING_MATRIX_NAME_AND_INDICES, &
            & '("    dZ_dNu','(",I1,",:)','   :",3(X,E13.6))','(19X,3(X,E13.6))',err,error,*999)
          CALL Determinant(dZdNu,detdZdNu,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Determinant dZ_dNu  = ",detdZdNu,err,error,*999)
        ENDIF

      ELSE
        CALL FlagError("Geometric interpolated point metrics is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Dependent interpolated point metrics is not associated.",err,error,*999)
    ENDIF

    EXITS("FiniteElasticity_GaussDeformationGradientTensor")
    RETURN
999 ERRORS("FiniteElasticity_GaussDeformationGradientTensor",err,error)
    EXITS("FiniteElasticity_GaussDeformationGradientTensor")
    RETURN 1

  END SUBROUTINE FiniteElasticity_GaussDeformationGradientTensor

  !
  !================================================================================================================================
  !

  !>Evaluates the Cauchy stress tensor at a given Gauss point
  SUBROUTINE FINITE_ELASTICITY_GAUSS_CAUCHY_TENSOR(EQUATIONS_SET,DEPENDENT_INTERPOLATED_POINT, &
      & MATERIALS_INTERPOLATED_POINT,GEOMETRIC_INTERPOLATED_POINT,DARCY_DEPENDENT_INTERPOLATED_POINT, &
      & INDEPENDENT_INTERPOLATED_POINT,CAUCHY_TENSOR,Jznu,DZDNU,ELEMENT_NUMBER,GAUSS_POINT_NUMBER,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER, INTENT(IN) :: EQUATIONS_SET !<A pointer to the equations set
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: DEPENDENT_INTERPOLATED_POINT,MATERIALS_INTERPOLATED_POINT
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: DARCY_DEPENDENT_INTERPOLATED_POINT,GEOMETRIC_INTERPOLATED_POINT
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: INDEPENDENT_INTERPOLATED_POINT
    REAL(DP), INTENT(OUT) :: CAUCHY_TENSOR(:,:)
    REAL(DP), INTENT(OUT) :: Jznu !Determinant of deformation gradient tensor (AZL)
    REAL(DP), INTENT(IN) :: DZDNU(3,3) !Deformation gradient tensor at the Guass point
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER,GAUSS_POINT_NUMBER !<Element/Gauss point number
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: EQUATIONS_SET_SUBTYPE !<The equation subtype
    INTEGER(INTG) :: i,j,k,PRESSURE_COMPONENT,component_idx,dof_idx
    REAL(DP) :: activation
    REAL(DP) :: AZL(3,3),AZU(3,3),DZDNUT(3,3),PIOLA_TENSOR(3,3),E(3,3),P,IDENTITY(3,3),AZLT(3,3),AZUT(3,3)
    REAL(DP) :: AZL_SQUARED(3,3)
    REAL(DP) :: I1,I2,I3            !Invariants, if needed
    REAL(DP) :: ACTIVE_STRESS_11,ACTIVE_STRESS_22,ACTIVE_STRESS_33 !Active stress to be copied in from independent field.
    REAL(DP) :: TEMP(3,3),TEMPTERM  !Temporary variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    REAL(DP), DIMENSION (:), POINTER :: C !Parameters for constitutive laws
    REAL(DP) :: a, B(3,3), Q !Parameters for orthotropic laws
    REAL(DP) :: ffact,dfdJfact !coupled elasticity Darcy
    INTEGER(INTG) :: DARCY_MASS_INCREASE_ENTRY !position of mass-increase entry in dependent-variable vector
    REAL(DP) :: VALUE,VAL1,VAL2
    REAL(DP) :: WV_PRIME,TOL,TOL1,UP,LOW
    REAL(DP) :: F_e(3,3),F_a(3,3),F_a_inv(3,3),F_a_T(3,3),C_a(3,3),C_a_inv(3,3),lambda_a,C_e(3,3),F_e_T(3,3)
    REAL(DP) :: REFERENCE_VOLUME,XB_STIFFNESS,XB_DISTORTION,V_MAX
    REAL(DP) :: SARCO_LENGTH,FREE_ENERGY,FREE_ENERGY_0,XB_ENERGY_PER_VOLUME,SLOPE,lambda_f,A_1,A_2,x_1,x_2
    REAL(DP) :: MAX_XB_NUMBER_PER_VOLUME,ENERGY_PER_XB,FORCE_LENGTH,I_1e,EVALUES(3),EVECTOR_1(3),EVECTOR_2(3),EVECTOR_3(3)
    REAL(DP) :: EMATRIX_1(3,3),EMATRIX_2(3,3),EMATRIX_3(3,3),TEMP1(3,3),TEMP2(3,3),TEMP3(3,3),N1(3,3),N2(3,3),N3(3,3) 
    REAL(DP), DIMENSION(5) :: PAR
    INTEGER(INTG) :: LWORK,node1,node2
    INTEGER(INTG), PARAMETER :: LWMAX=1000
    REAL(DP) :: WORK(LWMAX),RIGHT_NODE(3),LEFT_NODE(3),delta_t,dist1,dist2,velo
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,INDEPENDENT_FIELD
    REAL(DP) :: ISOMETRIC_FORCE_AT_FULL_ACT,LENGTH_HALF_SARCO
    REAL(DP) :: TITIN_VALUE,TITIN_VALUE_CROSS_FIBRE,TITIN_UNBOUND,TITIN_BOUND
    REAL(DP) :: TITIN_UNBOUND_CROSS_FIBRE,TITIN_BOUND_CROSS_FIBRE

    ENTERS("FINITE_ELASTICITY_GAUSS_CAUCHY_TENSOR",err,error,*999)

    NULLIFY(FIELD_VARIABLE)

    IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
      CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
      CALL FlagError("Equations set specification must have three entries for a finite elasticity type equations set.", &
        & err,error,*999)
    END IF
    EQUATIONS_SET_SUBTYPE=EQUATIONS_SET%SPECIFICATION(3)
    C=>MATERIALS_INTERPOLATED_POINT%VALUES(:,1)

    !AZL = F'*F (deformed covariant or right cauchy deformation tensor, C)
    !AZU - deformed contravariant tensor; I3 = det(C)
    !E = Green-Lagrange strain tensor = 0.5*(C-I)
    !PIOLA_TENSOR is the second Piola-Kirchoff tensor (PK2 or S)
    !P is the actual hydrostatic pressure, not double it

    CALL MatrixTranspose(DZDNU,DZDNUT,err,error,*999)
    CALL MatrixProduct(DZDNUT,DZDNU,AZL,err,error,*999)
    CALL Determinant(DZDNU,Jznu,err,error,*999)

    PRESSURE_COMPONENT=DEPENDENT_INTERPOLATED_POINT%INTERPOLATION_PARAMETERS%FIELD_VARIABLE%NUMBER_OF_COMPONENTS
    P=DEPENDENT_INTERPOLATED_POINT%VALUES(PRESSURE_COMPONENT,1)

    CALL INVERT(AZL,AZU,I3,ERR,ERROR,*999)
    
    E = 0.5_DP*AZL
    DO i=1,3
      E(i,i)=E(i,i)-0.5_DP
    ENDDO
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING_MATRIX(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
        & 3,3,E,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("    E','(",I1,",:)',' :",3(X,E13.6))', &
        & '(17X,3(X,E13.6))',err,error,*999)
    ENDIF
    IDENTITY=0.0_DP
    DO i=1,3
      IDENTITY(i,i)=1.0_DP
    ENDDO

    SELECT CASE(EQUATIONS_SET_SUBTYPE)
    CASE(EQUATIONS_SET_NEARLY_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE)
      !Form of constitutive model is:
      ! W_hat=c1*(I1_hat-3)+c2*(I2_hat-3)+p*J*C^(-1) + W^v(J)
      ! take W^v(J) = 1/2 * kappa * (J-1)^2
      WV_PRIME = C(3)*(Jznu - 1.0_DP)
      !compute the invariants, I3 a few lines up
      I1 = AZL(1,1) + AZL(2,2) + AZL(3,3)
      CALL MatrixProduct(AZL,AZL,AZL_SQUARED,err,error,*999)
      I2 = 0.5_DP * (I1**2 - AZL_SQUARED(1,1) - AZL_SQUARED(2,2) - AZL_SQUARED(3,3))

      PIOLA_TENSOR=2.0_DP*Jznu**(-2.0_DP/3.0_DP)*((C(1)+C(2)*I1)*IDENTITY-C(2)*AZL &
        & -(C(1)*I1+2.0_DP*C(2)*I2-1.5_DP*WV_PRIME*Jznu**(5.0_DP/3.0_DP))/3.0_DP*AZU)

    CASE(EQUATIONS_SET_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE)
      !Form of constitutive model is:
      ! W_hat=c1*(I1_hat-3)+c2*(I2_hat-3)+p*J*C^(-1)

      !compute the invariants, I3 a few lines up
      I1 = AZL(1,1) + AZL(2,2) + AZL(3,3)
      CALL MatrixProduct(AZL,AZL,AZL_SQUARED,err,error,*999)
      I2 = 0.5_DP * (I1**2 - AZL_SQUARED(1,1) - AZL_SQUARED(2,2) - AZL_SQUARED(3,3))

      !compute 2PK
!      PIOLA_TENSOR(1,1) = 2.0_DP * Jznu**(-2.0_DP/3.0_DP) * (C(1) + C(2) * I1 - C(2) * AZL(1,1) &
!                          & - (C(1) * I1 + 2.0_DP * C(2) * I2 - 1.5_DP * P * Jznu**(5.0_DP/3.0_DP)) / 3.0_DP * AZU(1,1))
!      PIOLA_TENSOR(1,2) = 2.0_DP * Jznu**(-2.0_DP/3.0_DP) * (-C(2) * AZL(1,2) &
!                          & - (C(1) * I1 + 2.0_DP * C(2) * I2 - 1.5_DP * P * Jznu**(5.0_DP/3.0_DP)) / 3.0_DP * AZU(1,2))
!      PIOLA_TENSOR(1,3) = 2.0_DP * Jznu**(-2.0_DP/3.0_DP) * (-C(2) * AZL(1,3) &
!                          & - (C(1) * I1 + 2.0_DP * C(2) * I2 - 1.5_DP * P * Jznu**(5.0_DP/3.0_DP)) / 3.0_DP * AZU(1,3))
!      PIOLA_TENSOR(2,1) = PIOLA_TENSOR(1,2)
!      PIOLA_TENSOR(2,2) = 2.0_DP * Jznu**(-2.0_DP/3.0_DP) * (C(1) + C(2) * I1 - C(2) * AZL(2,2) &
!                          & - (C(1) * I1 + 2.0_DP * C(2) * I2 - 1.5_DP * P * Jznu**(5.0_DP/3.0_DP)) / 3.0_DP * AZU(2,2))
!      PIOLA_TENSOR(2,3) = 2.0_DP * Jznu**(-2.0_DP/3.0_DP) * (-C(2) * AZL(2,3) &
!                          & - (C(1) * I1 + 2.0_DP * C(2) * I2 - 1.5_DP * P * Jznu**(5.0_DP/3.0_DP)) / 3.0_DP * AZU(2,3))
!      PIOLA_TENSOR(3,1) = PIOLA_TENSOR(1,3)
!      PIOLA_TENSOR(3,2) = PIOLA_TENSOR(2,3)
!      PIOLA_TENSOR(3,3) = 2.0_DP * Jznu**(-2.0_DP/3.0_DP) * (C(1) + C(2) * I1 - C(2) * AZL(3,3) &
!                          & - (C(1) * I1 + 2.0_DP * C(2) * I2 - 1.5_DP * P * Jznu**(5.0_DP/3.0_DP)) / 3.0_DP * AZU(3,3))
      !????
      PIOLA_TENSOR=2.0_DP*Jznu**(-2.0_DP/3.0_DP)*((C(1)+C(2)*I1)*IDENTITY-C(2)*AZL &
        & -(C(1)*I1+2.0_DP*C(2)*I2-1.5_DP*P*Jznu**(5.0_DP/3.0_DP))/3.0_DP*AZU)

    CASE(EQUATIONS_SET_ACTIVE_STRAIN_SUBTYPE)

      DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
      node1=dependent_field%decomposition%domain(1)%ptr%topology%elements%elements(ELEMENT_NUMBER)%element_nodes(13)
      node2=dependent_field%decomposition%domain(1)%ptr%topology%elements%elements(ELEMENT_NUMBER)%element_nodes(15)

      NULLIFY(FIELD_VARIABLE)
      ! compute the nodal distance of the previous time step
      CALL Field_VariableGet(DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VARIABLE,err,error,*999)
      dof_idx=FIELD_VARIABLE%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node1)%DERIVATIVES(1)%VERSIONS(1)
      CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,dof_idx,LEFT_NODE(1), &
        & err,error,*999)
      dof_idx=FIELD_VARIABLE%COMPONENTS(2)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node1)%DERIVATIVES(1)%VERSIONS(1)
      CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,dof_idx,LEFT_NODE(2), &
        & err,error,*999)
      dof_idx=FIELD_VARIABLE%COMPONENTS(3)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node1)%DERIVATIVES(1)%VERSIONS(1)
      CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,dof_idx,LEFT_NODE(3), &
        & err,error,*999)

      dof_idx=FIELD_VARIABLE%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node2)%DERIVATIVES(1)%VERSIONS(1)
      CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,dof_idx,RIGHT_NODE(1), &
        & err,error,*999)
      dof_idx=FIELD_VARIABLE%COMPONENTS(2)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node2)%DERIVATIVES(1)%VERSIONS(1)
      CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,dof_idx,RIGHT_NODE(2), &
        & err,error,*999)
      dof_idx=FIELD_VARIABLE%COMPONENTS(3)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node2)%DERIVATIVES(1)%VERSIONS(1)
      CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,dof_idx,RIGHT_NODE(3), &
        & err,error,*999)

      dist1=SQRT((RIGHT_NODE(1)-LEFT_NODE(1))*(RIGHT_NODE(1)-LEFT_NODE(1))+ &
               & (RIGHT_NODE(2)-LEFT_NODE(2))*(RIGHT_NODE(2)-LEFT_NODE(2))+ &
               & (RIGHT_NODE(3)-LEFT_NODE(3))*(RIGHT_NODE(3)-LEFT_NODE(3)))

      NULLIFY(FIELD_VARIABLE)
      ! compute the nodal distance of the current time step
      CALL Field_VariableGet(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VARIABLE,err,error,*999)
      dof_idx=FIELD_VARIABLE%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node1)%DERIVATIVES(1)%VERSIONS(1)
      CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,dof_idx,LEFT_NODE(1), &
        & err,error,*999)
      dof_idx=FIELD_VARIABLE%COMPONENTS(2)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node1)%DERIVATIVES(1)%VERSIONS(1)
      CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,dof_idx,LEFT_NODE(2), &
        & err,error,*999)
      dof_idx=FIELD_VARIABLE%COMPONENTS(3)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node1)%DERIVATIVES(1)%VERSIONS(1)
      CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,dof_idx,LEFT_NODE(3), &
        & err,error,*999)

      dof_idx=FIELD_VARIABLE%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node2)%DERIVATIVES(1)%VERSIONS(1)
      CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,dof_idx,RIGHT_NODE(1), &
        & err,error,*999)
      dof_idx=FIELD_VARIABLE%COMPONENTS(2)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node2)%DERIVATIVES(1)%VERSIONS(1)
      CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,dof_idx,RIGHT_NODE(2), &
        & err,error,*999)
      dof_idx=FIELD_VARIABLE%COMPONENTS(3)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node2)%DERIVATIVES(1)%VERSIONS(1)
      CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,dof_idx,RIGHT_NODE(3), &
        & err,error,*999)

      dist2=SQRT((RIGHT_NODE(1)-LEFT_NODE(1))*(RIGHT_NODE(1)-LEFT_NODE(1))+ &
               & (RIGHT_NODE(2)-LEFT_NODE(2))*(RIGHT_NODE(2)-LEFT_NODE(2))+ &
               & (RIGHT_NODE(3)-LEFT_NODE(3))*(RIGHT_NODE(3)-LEFT_NODE(3)))

      delta_t=0.01_DP;
      velo=(dist2-dist1)/delta_t ! velo>0 for lengthening
!      velo=(dist1-dist2)/delta_t ! velo<0 for shortening
      !velo=velo*1.0e-6_DP
      velo=velo*5.0e-8_DP 
      
      !--------------------------------------------------------------------------------------------

      !Force-Velocity-Relation
!      PAR=[1.0_DP,0.5_DP,0.5_DP,0.8_DP,0.2_DP] ! Muscle-Parameters for F-v-Relation
!      IF(velo.GE.0.0_DP) THEN
!        ENERGY_PER_XB=(PAR(1)+PAR(2))*PAR(3)/(velo+PAR(3))-PAR(2)
!      ELSE
!        ENERGY_PER_XB=((2.0_DP*PAR(1)-PAR(4))*velo-PAR(1)*PAR(5))/(velo-PAR(5))
!      ENDIF
      V_MAX=8.9e-8_DP
      XB_DISTORTION=8.0e-9_DP*(1+velo/V_MAX) ! [m]
      
      XB_STIFFNESS=2.2e-3_DP ! [N/m]

      REFERENCE_VOLUME=1.4965e+06_DP ! [nm^3]
      MAX_XB_NUMBER_PER_VOLUME=120.0_DP*2.0_DP/REFERENCE_VOLUME ! [cross-bridges per nm^3]
      ENERGY_PER_XB=0.5_DP*XB_STIFFNESS*XB_DISTORTION**2 ! [J]

      SARCO_LENGTH=DZDNU(1,1)
      
      ! Calculate Filament-Overlap
      IF(SARCO_LENGTH.LE.0.635_DP) THEN
        FORCE_LENGTH=0.0_DP
      ELSE IF(SARCO_LENGTH.LE.0.835_DP) THEN 
        FORCE_LENGTH=4.2_DP*(SARCO_LENGTH-0.635_DP)
      ELSE IF(SARCO_LENGTH.LE.1.0_DP) THEN
        FORCE_LENGTH=0.84_DP+0.9697_DP*(SARCO_LENGTH-0.835_DP)
      ELSE IF(SARCO_LENGTH.LE.1.125_DP) THEN
        FORCE_LENGTH=1.0_DP
      ELSE IF(SARCO_LENGTH.LE.1.825_DP) THEN
        FORCE_LENGTH=1.0_DP-1.4286_DP*(SARCO_LENGTH-1.125_DP)
      ELSE
        FORCE_LENGTH=0.0_DP
      ENDIF
      
      !Mechanical Energy stored in cross-bridges - conversion from J/nm^3 to N/cm^2           
      XB_ENERGY_PER_VOLUME=MAX_XB_NUMBER_PER_VOLUME*FORCE_LENGTH*ENERGY_PER_XB*A_2*10.0_DP**23

      !Initalize lambda_a
      lambda_a=1.0_DP
      
      F_a_inv=0.0_DP
      F_a_inv(1,1)=1.0_DP/lambda_a
      F_a_inv(2,2)=1.0_DP
      F_a_inv(3,3)=1.0_DP

      CALL MatrixProduct(DZDNU,F_a_inv,F_e,err,error,*999)
      CALL MatrixTranspose(F_e,F_e_T,err,error,*999)
      CALL MatrixProduct(F_e_T,F_e,C_e,err,error,*999)

      !Odgen law - 3 terms. Material Parameters C = [mu(1) mu(2) mu(3) alpha(1) alpha(2) alpha(3) mu_0]
      
!      CALL Eigenvalue(C_e,EVALUES,err,error,*999)
      CALL DSYEV('V','U',3,C_e,3,EVALUES,WORK,-1,ERR)
      IF(ERR.NE.0) CALL FlagError("Error in Eigenvalue computation",err,error,*999)
      LWORK=MIN(LWMAX,INT(WORK(1)))
      CALL DSYEV('V','U',3,C_e,3,EVALUES,WORK,LWORK,ERR)
      IF(ERR.NE.0) CALL FlagError("Error in Eigenvalue computation",err,error,*999)
      EVECTOR_1=C_e(:,1)
      EVECTOR_2=C_e(:,2)
      EVECTOR_3=C_e(:,3)

      DO i=1,3
        DO j=1,3
          EMATRIX_1(i,j)=EVECTOR_1(i)*EVECTOR_1(j)
          EMATRIX_2(i,j)=EVECTOR_2(i)*EVECTOR_2(j)
          EMATRIX_3(i,j)=EVECTOR_3(i)*EVECTOR_3(j)
        END DO
      END DO

      CALL MatrixProduct(F_a_inv,EMATRIX_1,N1,err,error,*999)
      CALL MatrixProduct(N1,F_a_inv,N1,err,error,*999) ! F_a_inv=F_a_inv_T
      CALL MatrixProduct(F_a_inv,EMATRIX_2,N2,err,error,*999)
      CALL MatrixProduct(N2,F_a_inv,N2,err,error,*999) ! F_a_inv=F_a_inv_T
      CALL MatrixProduct(F_a_inv,EMATRIX_3,N3,err,error,*999)
      CALL MatrixProduct(N3,F_a_inv,N3,err,error,*999) ! F_a_inv=F_a_inv_T 

      FREE_ENERGY_0=0.0_DP
      DO i=1,3
        FREE_ENERGY_0=FREE_ENERGY_0+C(i)/C(i+3)*( &
          & EVALUES(1)**(C(i+3)/2.0_DP)+ &
          & EVALUES(2)**(C(i+3)/2.0_DP)+ &
          & EVALUES(3)**(C(i+3)/2.0_DP)-3.0_DP)
      END DO
      FREE_ENERGY_0=C(7)*FREE_ENERGY_0

      FREE_ENERGY=FREE_ENERGY_0

      VALUE=XB_ENERGY_PER_VOLUME-(FREE_ENERGY-FREE_ENERGY_0)

      !tolerance for Newton's method
      TOL=0.00001_DP
      !tolerance for the bisection method as preconditioner. Since Newton's method does not converge, we only use the bisection method here
      TOL1=TOL 
      UP=lambda_a
      LOW=0.001_DP
      
!      WRITE(*,*) "VALUE: ", VALUE

      DO WHILE (ABS(VALUE).GE.TOL)

        !bisection method
        IF (ABS(VALUE).GE.TOL1) THEN
          lambda_a=UP-(UP-LOW)/2.0_DP

          F_a_inv=0.0_DP
          IF(lambda_a<TOL) THEN
           CALL FlagWarning("lambda_a is close to zero",err,error,*999)
!            WRITE(*,*) "UP: ", UP
!            WRITE(*,*) "LOW: ", LOW
!            WRITE(*,*) "lambda_a: ", lambda_a
            lambda_a=lambda_a+TOL
          ENDIF
          F_a_inv(1,1)=1.0_DP/lambda_a
          F_a_inv(2,2)=1.0_DP
          F_a_inv(3,3)=1.0_DP

          CALL MatrixProduct(DZDNU,F_a_inv,F_e,err,error,*999)
          CALL MatrixTranspose(F_e,F_e_T,err,error,*999)
          CALL MatrixProduct(F_e_T,F_e,C_e,err,error,*999)

          CALL DSYEV('V','U',3,C_e,3,EVALUES,WORK,-1,ERR)
          IF(ERR.NE.0) CALL FlagError("Error in Eigenvalue computation",err,error,*999)
          LWORK=MIN(LWMAX,INT(WORK(1)))
          CALL DSYEV('V','U',3,C_e,3,EVALUES,WORK,LWORK,ERR)
          IF(ERR.NE.0) CALL FlagError("Error in Eigenvalue computation",err,error,*999)
          EVECTOR_1=C_e(:,1)
          EVECTOR_2=C_e(:,2)
          EVECTOR_3=C_e(:,3)

          DO i=1,3
            DO j=1,3
              EMATRIX_1(i,j)=EVECTOR_1(i)*EVECTOR_1(j)
              EMATRIX_2(i,j)=EVECTOR_2(i)*EVECTOR_2(j)
              EMATRIX_3(i,j)=EVECTOR_3(i)*EVECTOR_3(j)
            END DO
          END DO

          CALL MatrixProduct(F_a_inv,EMATRIX_1,N1,err,error,*999)
          CALL MatrixProduct(N1,F_a_inv,N1,err,error,*999) ! F_a_inv=F_a_inv_T
          CALL MatrixProduct(F_a_inv,EMATRIX_2,N2,err,error,*999)
          CALL MatrixProduct(N2,F_a_inv,N2,err,error,*999) ! F_a_inv=F_a_inv_T
          CALL MatrixProduct(F_a_inv,EMATRIX_3,N3,err,error,*999)
          CALL MatrixProduct(N3,F_a_inv,N3,err,error,*999) ! F_a_inv=F_a_inv_T 

          FREE_ENERGY=0.0_DP
          DO i=1,3
            FREE_ENERGY=FREE_ENERGY+C(i)/C(i+3)*( &
              & EVALUES(1)**(C(i+3)/2.0_DP)+ &
              & EVALUES(2)**(C(i+3)/2.0_DP)+ &
              & EVALUES(3)**(C(i+3)/2.0_DP)-3.0_DP)
          END DO
          FREE_ENERGY=C(7)*FREE_ENERGY

          VALUE=XB_ENERGY_PER_VOLUME-(FREE_ENERGY-FREE_ENERGY_0)

          IF (VALUE.GE.0) THEN
            UP=lambda_a
          ELSE
            LOW=lambda_a
          ENDIF

        ELSE 
          !Newton's method -- needs to be checked TODO

          TEMP=DZDNU+DZDNUT
          CALL MatrixProduct(F_e_T,TEMP,TEMP,err,error,*999)
          CALL MatrixProduct(TEMP,N1,TEMP1,err,error,*999) 
          CALL MatrixProduct(TEMP,N2,TEMP2,err,error,*999) 
          CALL MatrixProduct(TEMP,N3,TEMP3,err,error,*999) 

          TEMP=0.0_DP
          DO i=1,3
            TEMP=TEMP+ &
              & C(i)*EVALUES(1)**(C(i+3)/2.0_DP-1.0_DP)*TEMP1+ &
              & C(i)*EVALUES(2)**(C(i+3)/2.0_DP-1.0_DP)*TEMP2+ &
              & C(i)*EVALUES(3)**(C(i+3)/2.0_DP-1.0_DP)*TEMP3
          END DO
          SLOPE=TEMP(1,1)*C(7)
          lambda_a=lambda_a-VALUE/SLOPE
          !IF (lambda_a.LE.0.0_DP) THEN
          ! lambda_a=0.1_DP
          !END IF
          !lambda_a=lambda_a-0.001

          F_a_inv=0.0_DP
          F_a_inv(1,1)=1.0_DP/lambda_a
          F_a_inv(2,2)=1.0_DP
          F_a_inv(3,3)=1.0_DP

          CALL MatrixProduct(DZDNU,F_a_inv,F_e,err,error,*999)
          CALL MatrixTranspose(F_e,F_e_T,err,error,*999)
          CALL MatrixProduct(F_e_T,F_e,C_e,err,error,*999)

          CALL DSYEV('V','U',3,C_e,3,EVALUES,WORK,-1,ERR)
          IF(ERR.NE.0) CALL FlagError("Error in Eigenvalue computation",err,error,*999)
          LWORK=MIN(LWMAX,INT(WORK(1)))
          CALL DSYEV('V','U',3,C_e,3,EVALUES,WORK,LWORK,ERR)
          IF(ERR.NE.0) CALL FlagError("Error in Eigenvalue computation",err,error,*999)
          EVECTOR_1=C_e(:,1)
          EVECTOR_2=C_e(:,2)
          EVECTOR_3=C_e(:,3)

          DO i=1,3
            DO j=1,3
              EMATRIX_1(i,j)=EVECTOR_1(i)*EVECTOR_1(j)
              EMATRIX_2(i,j)=EVECTOR_2(i)*EVECTOR_2(j)
              EMATRIX_3(i,j)=EVECTOR_3(i)*EVECTOR_3(j)
            END DO
          END DO

          CALL MatrixProduct(F_a_inv,EMATRIX_1,N1,err,error,*999)
          CALL MatrixProduct(N1,F_a_inv,N1,err,error,*999) ! F_a_inv=F_a_inv_T
          CALL MatrixProduct(F_a_inv,EMATRIX_2,N2,err,error,*999)
          CALL MatrixProduct(N2,F_a_inv,N2,err,error,*999) ! F_a_inv=F_a_inv_T
          CALL MatrixProduct(F_a_inv,EMATRIX_3,N3,err,error,*999)
          CALL MatrixProduct(N3,F_a_inv,N3,err,error,*999) ! F_a_inv=F_a_inv_T 

          FREE_ENERGY=0.0_DP
          DO i=1,3
            FREE_ENERGY=FREE_ENERGY+C(i)/C(i+3)*( &
              & EVALUES(1)**(C(i+3)/2.0_DP)+ &
              & EVALUES(2)**(C(i+3)/2.0_DP)+ &
              & EVALUES(3)**(C(i+3)/2.0_DP)-3.0_DP)
          END DO
          FREE_ENERGY=C(7)*FREE_ENERGY

          VALUE=XB_ENERGY_PER_VOLUME-(FREE_ENERGY-FREE_ENERGY_0)
        ENDIF
      ENDDO
    
      PIOLA_TENSOR=0.0_DP
      DO i=1,3
        PIOLA_TENSOR=PIOLA_TENSOR+ &
          & C(i)*EVALUES(1)**(C(i+3)/2.0_DP-1.0_DP)*N1+ &
          & C(i)*EVALUES(2)**(C(i+3)/2.0_DP-1.0_DP)*N2+ &
          & C(i)*EVALUES(3)**(C(i+3)/2.0_DP-1.0_DP)*N3
      END DO
      PIOLA_TENSOR=PIOLA_TENSOR*C(7)+2.0_DP*P*AZU

    CASE(EQUATIONS_SET_MULTISCALE_ACTIVE_STRAIN_SUBTYPE)

      DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
      node1=dependent_field%decomposition%domain(1)%ptr%topology%elements%elements(ELEMENT_NUMBER)%element_nodes(13)
      node2=dependent_field%decomposition%domain(1)%ptr%topology%elements%elements(ELEMENT_NUMBER)%element_nodes(15)

      NULLIFY(FIELD_VARIABLE)
      ! compute the nodal distance of the previous time step
      CALL Field_VariableGet(DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VARIABLE,err,error,*999)
      dof_idx=FIELD_VARIABLE%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node1)%DERIVATIVES(1)%VERSIONS(1)
      CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,dof_idx,LEFT_NODE(1), &
        & err,error,*999)
      dof_idx=FIELD_VARIABLE%COMPONENTS(2)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node1)%DERIVATIVES(1)%VERSIONS(1)
      CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,dof_idx,LEFT_NODE(2), &
        & err,error,*999)
      dof_idx=FIELD_VARIABLE%COMPONENTS(3)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node1)%DERIVATIVES(1)%VERSIONS(1)
      CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,dof_idx,LEFT_NODE(3), &
        & err,error,*999)

      dof_idx=FIELD_VARIABLE%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node2)%DERIVATIVES(1)%VERSIONS(1)
      CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,dof_idx,RIGHT_NODE(1), &
        & err,error,*999)
      dof_idx=FIELD_VARIABLE%COMPONENTS(2)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node2)%DERIVATIVES(1)%VERSIONS(1)
      CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,dof_idx,RIGHT_NODE(2), &
        & err,error,*999)
      dof_idx=FIELD_VARIABLE%COMPONENTS(3)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node2)%DERIVATIVES(1)%VERSIONS(1)
      CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,dof_idx,RIGHT_NODE(3), &
        & err,error,*999)

      dist1=SQRT((RIGHT_NODE(1)-LEFT_NODE(1))*(RIGHT_NODE(1)-LEFT_NODE(1))+ &
               & (RIGHT_NODE(2)-LEFT_NODE(2))*(RIGHT_NODE(2)-LEFT_NODE(2))+ &
               & (RIGHT_NODE(3)-LEFT_NODE(3))*(RIGHT_NODE(3)-LEFT_NODE(3)))

      NULLIFY(FIELD_VARIABLE)
      ! compute the nodal distance of the current time step
      CALL Field_VariableGet(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VARIABLE,err,error,*999)
      dof_idx=FIELD_VARIABLE%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node1)%DERIVATIVES(1)%VERSIONS(1)
      CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,dof_idx,LEFT_NODE(1), &
        & err,error,*999)
      dof_idx=FIELD_VARIABLE%COMPONENTS(2)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node1)%DERIVATIVES(1)%VERSIONS(1)
      CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,dof_idx,LEFT_NODE(2), &
        & err,error,*999)
      dof_idx=FIELD_VARIABLE%COMPONENTS(3)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node1)%DERIVATIVES(1)%VERSIONS(1)
      CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,dof_idx,LEFT_NODE(3), &
        & err,error,*999)

      dof_idx=FIELD_VARIABLE%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node2)%DERIVATIVES(1)%VERSIONS(1)
      CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,dof_idx,RIGHT_NODE(1), &
        & err,error,*999)
      dof_idx=FIELD_VARIABLE%COMPONENTS(2)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node2)%DERIVATIVES(1)%VERSIONS(1)
      CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,dof_idx,RIGHT_NODE(2), &
        & err,error,*999)
      dof_idx=FIELD_VARIABLE%COMPONENTS(3)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node2)%DERIVATIVES(1)%VERSIONS(1)
      CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,dof_idx,RIGHT_NODE(3), &
        & err,error,*999)

      dist2=SQRT((RIGHT_NODE(1)-LEFT_NODE(1))*(RIGHT_NODE(1)-LEFT_NODE(1))+ &
               & (RIGHT_NODE(2)-LEFT_NODE(2))*(RIGHT_NODE(2)-LEFT_NODE(2))+ &
               & (RIGHT_NODE(3)-LEFT_NODE(3))*(RIGHT_NODE(3)-LEFT_NODE(3)))

      delta_t=0.001_DP;
      velo=(dist2-dist1)/delta_t ! velo>0 == lengthening
      !conversion of velocity at the continuum macroscale to the micromechanical cell model half-sarcomere velocity
      velo=velo*5.0e-8_DP 
!      velo=velo*5.0e-2_DP
!      velo=velo*5.0e-7_DP 

      CALL FIELD_PARAMETER_SET_UPDATE_GAUSS_POINT(DEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,GAUSS_POINT_NUMBER, &
        & ELEMENT_NUMBER,2,velo,err,error,*999)

      
      !--------------------------------------------------------------------------------------------
      NULLIFY(INDEPENDENT_FIELD)
      INDEPENDENT_FIELD=>EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD
      NULLIFY(FIELD_VARIABLE)
      CALL Field_VariableGet(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VARIABLE,err,error,*999)

      dof_idx=FIELD_VARIABLE%COMPONENTS(1)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(GAUSS_POINT_NUMBER, &
        & ELEMENT_NUMBER)
      CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,dof_idx,A_1, &
        & err,error,*999)
      dof_idx=FIELD_VARIABLE%COMPONENTS(2)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(GAUSS_POINT_NUMBER, &
        & ELEMENT_NUMBER)
      CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,dof_idx,A_2, &
        & err,error,*999)
      dof_idx=FIELD_VARIABLE%COMPONENTS(3)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(GAUSS_POINT_NUMBER, &
        & ELEMENT_NUMBER)
      CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,dof_idx,x_1, &
        & err,error,*999)
      dof_idx=FIELD_VARIABLE%COMPONENTS(4)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(GAUSS_POINT_NUMBER, &
        & ELEMENT_NUMBER)
      CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,dof_idx,x_2, &
        & err,error,*999)

      !--------------------------------------------------------------------------------------------
      SARCO_LENGTH=DZDNU(1,1)
      ! Calculate Filament-Overlap
      IF(SARCO_LENGTH.LE.0.635_DP) THEN
        FORCE_LENGTH=0.0_DP
      ELSE IF(SARCO_LENGTH.LE.0.835_DP) THEN 
        FORCE_LENGTH=4.2_DP*(SARCO_LENGTH-0.635_DP)
      ELSE IF(SARCO_LENGTH.LE.1.0_DP) THEN
        FORCE_LENGTH=0.84_DP+0.9697_DP*(SARCO_LENGTH-0.835_DP)
      ELSE IF(SARCO_LENGTH.LE.1.125_DP) THEN
        FORCE_LENGTH=1.0_DP
      ELSE IF(SARCO_LENGTH.LE.1.825_DP) THEN
        FORCE_LENGTH=1.0_DP-1.4286_DP*(SARCO_LENGTH-1.125_DP)
      ELSE
        FORCE_LENGTH=0.0_DP
      ENDIF

      REFERENCE_VOLUME=1.4965e+06_DP ! [nm^3]
      MAX_XB_NUMBER_PER_VOLUME=120.0_DP*2.0_DP/REFERENCE_VOLUME ! [cross-bridges per nm^3]
      ENERGY_PER_XB=0.5_DP*x_2**2*C(8) ! joule
      
      !Mechanical Energy stored in cross-bridges - conversion from J/nm^3 to N/cm^2
      XB_ENERGY_PER_VOLUME=MAX_XB_NUMBER_PER_VOLUME*FORCE_LENGTH*ENERGY_PER_XB*A_2*10.0_DP**23

      !Initalize lambda_a
      lambda_a=1.0_DP
      
      F_a_inv=0.0_DP
      F_a_inv(1,1)=1.0_DP/lambda_a
      F_a_inv(2,2)=1.0_DP
      F_a_inv(3,3)=1.0_DP

      CALL MatrixProduct(DZDNU,F_a_inv,F_e,err,error,*999)
      CALL MatrixTranspose(F_e,F_e_T,err,error,*999)
      CALL MatrixProduct(F_e_T,F_e,C_e,err,error,*999)

      !Odgen law - 3 terms. Material Parameters C = [mu(1) mu(2) mu(3) alpha(1) alpha(2) alpha(3) mu_0]
!      CALL Eigenvalue(C_e,EVALUES,err,error,*999)
      CALL DSYEV('V','U',3,C_e,3,EVALUES,WORK,-1,ERR)
      IF(ERR.NE.0) CALL FlagError("Error in Eigenvalue computation",err,error,*999)
      LWORK=MIN(LWMAX,INT(WORK(1)))
      CALL DSYEV('V','U',3,C_e,3,EVALUES,WORK,LWORK,ERR)
      IF(ERR.NE.0) CALL FlagError("Error in Eigenvalue computation",err,error,*999)
      EVECTOR_1=C_e(:,1)
      EVECTOR_2=C_e(:,2)
      EVECTOR_3=C_e(:,3)

      DO i=1,3
        DO j=1,3
          EMATRIX_1(i,j)=EVECTOR_1(i)*EVECTOR_1(j)
          EMATRIX_2(i,j)=EVECTOR_2(i)*EVECTOR_2(j)
          EMATRIX_3(i,j)=EVECTOR_3(i)*EVECTOR_3(j)
        END DO
      END DO

      CALL MatrixProduct(F_a_inv,EMATRIX_1,N1,err,error,*999)
      CALL MatrixProduct(N1,F_a_inv,N1,err,error,*999) ! F_a_inv=F_a_inv_T
      CALL MatrixProduct(F_a_inv,EMATRIX_2,N2,err,error,*999)
      CALL MatrixProduct(N2,F_a_inv,N2,err,error,*999) ! F_a_inv=F_a_inv_T
      CALL MatrixProduct(F_a_inv,EMATRIX_3,N3,err,error,*999)
      CALL MatrixProduct(N3,F_a_inv,N3,err,error,*999) ! F_a_inv=F_a_inv_T 

      FREE_ENERGY_0=0.0_DP
      DO i=1,3
        FREE_ENERGY_0=FREE_ENERGY_0+C(i)/C(i+3)*( &
          & EVALUES(1)**(C(i+3)/2.0_DP)+ &
          & EVALUES(2)**(C(i+3)/2.0_DP)+ &
          & EVALUES(3)**(C(i+3)/2.0_DP)-3.0_DP)
      END DO
      FREE_ENERGY_0=C(7)*FREE_ENERGY_0

      FREE_ENERGY=FREE_ENERGY_0

      VALUE=XB_ENERGY_PER_VOLUME-(FREE_ENERGY-FREE_ENERGY_0)

      !tolerance for Newton's method
      TOL=0.00001_DP
      !tolerance for the bisection method as preconditioner. Since Newton's method does not converge, we only use the bisection method here
      TOL1=TOL 
      UP=lambda_a
      LOW=0.001_DP
      
!      WRITE(*,*) "VALUE: ", VALUE

      DO WHILE (ABS(VALUE).GE.TOL)

        !bisection method
        IF (ABS(VALUE).GE.TOL1) THEN
          lambda_a=UP-(UP-LOW)/2.0_DP

          F_a_inv=0.0_DP
          IF(lambda_a<TOL) THEN
           CALL FlagWarning("lambda_a is close to zero",err,error,*999)
!            WRITE(*,*) "UP: ", UP
!            WRITE(*,*) "LOW: ", LOW
!            WRITE(*,*) "lambda_a: ", lambda_a
            lambda_a=lambda_a+TOL
          ENDIF
          F_a_inv(1,1)=1.0_DP/lambda_a
          F_a_inv(2,2)=1.0_DP
          F_a_inv(3,3)=1.0_DP

          CALL MatrixProduct(DZDNU,F_a_inv,F_e,err,error,*999)
          CALL MatrixTranspose(F_e,F_e_T,err,error,*999)
          CALL MatrixProduct(F_e_T,F_e,C_e,err,error,*999)

          CALL DSYEV('V','U',3,C_e,3,EVALUES,WORK,-1,ERR)
          IF(ERR.NE.0) CALL FlagError("Error in Eigenvalue computation",err,error,*999)
          LWORK=MIN(LWMAX,INT(WORK(1)))
          CALL DSYEV('V','U',3,C_e,3,EVALUES,WORK,LWORK,ERR)
          IF(ERR.NE.0) CALL FlagError("Error in Eigenvalue computation",err,error,*999)
          EVECTOR_1=C_e(:,1)
          EVECTOR_2=C_e(:,2)
          EVECTOR_3=C_e(:,3)

          DO i=1,3
            DO j=1,3
              EMATRIX_1(i,j)=EVECTOR_1(i)*EVECTOR_1(j)
              EMATRIX_2(i,j)=EVECTOR_2(i)*EVECTOR_2(j)
              EMATRIX_3(i,j)=EVECTOR_3(i)*EVECTOR_3(j)
            END DO
          END DO

          CALL MatrixProduct(F_a_inv,EMATRIX_1,N1,err,error,*999)
          CALL MatrixProduct(N1,F_a_inv,N1,err,error,*999) ! F_a_inv=F_a_inv_T
          CALL MatrixProduct(F_a_inv,EMATRIX_2,N2,err,error,*999)
          CALL MatrixProduct(N2,F_a_inv,N2,err,error,*999) ! F_a_inv=F_a_inv_T
          CALL MatrixProduct(F_a_inv,EMATRIX_3,N3,err,error,*999)
          CALL MatrixProduct(N3,F_a_inv,N3,err,error,*999) ! F_a_inv=F_a_inv_T 

          FREE_ENERGY=0.0_DP
          DO i=1,3
            FREE_ENERGY=FREE_ENERGY+C(i)/C(i+3)*( &
              & EVALUES(1)**(C(i+3)/2.0_DP)+ &
              & EVALUES(2)**(C(i+3)/2.0_DP)+ &
              & EVALUES(3)**(C(i+3)/2.0_DP)-3.0_DP)
          END DO
          FREE_ENERGY=C(7)*FREE_ENERGY

          VALUE=XB_ENERGY_PER_VOLUME-(FREE_ENERGY-FREE_ENERGY_0)

          IF (VALUE.GE.0) THEN
            UP=lambda_a
          ELSE
            LOW=lambda_a
          ENDIF

        ELSE 
          !Newton's method -- needs to be checked TODO

          TEMP=DZDNU+DZDNUT
          CALL MatrixProduct(F_e_T,TEMP,TEMP,err,error,*999)
          CALL MatrixProduct(TEMP,N1,TEMP1,err,error,*999) 
          CALL MatrixProduct(TEMP,N2,TEMP2,err,error,*999) 
          CALL MatrixProduct(TEMP,N3,TEMP3,err,error,*999) 

          TEMP=0.0_DP
          DO i=1,3
            TEMP=TEMP+ &
              & C(i)*EVALUES(1)**(C(i+3)/2.0_DP-1.0_DP)*TEMP1+ &
              & C(i)*EVALUES(2)**(C(i+3)/2.0_DP-1.0_DP)*TEMP2+ &
              & C(i)*EVALUES(3)**(C(i+3)/2.0_DP-1.0_DP)*TEMP3
          END DO
          SLOPE=TEMP(1,1)*C(7)
          lambda_a=lambda_a-VALUE/SLOPE
          !IF (lambda_a.LE.0.0_DP) THEN
          ! lambda_a=0.1_DP
          !END IF
          !lambda_a=lambda_a-0.001

          F_a_inv=0.0_DP
          F_a_inv(1,1)=1.0_DP/lambda_a
          F_a_inv(2,2)=1.0_DP
          F_a_inv(3,3)=1.0_DP

          CALL MatrixProduct(DZDNU,F_a_inv,F_e,err,error,*999)
          CALL MatrixTranspose(F_e,F_e_T,err,error,*999)
          CALL MatrixProduct(F_e_T,F_e,C_e,err,error,*999)

          CALL DSYEV('V','U',3,C_e,3,EVALUES,WORK,-1,ERR)
          IF(ERR.NE.0) CALL FlagError("Error in Eigenvalue computation",err,error,*999)
          LWORK=MIN(LWMAX,INT(WORK(1)))
          CALL DSYEV('V','U',3,C_e,3,EVALUES,WORK,LWORK,ERR)
          IF(ERR.NE.0) CALL FlagError("Error in Eigenvalue computation",err,error,*999)
          EVECTOR_1=C_e(:,1)
          EVECTOR_2=C_e(:,2)
          EVECTOR_3=C_e(:,3)

          DO i=1,3
            DO j=1,3
              EMATRIX_1(i,j)=EVECTOR_1(i)*EVECTOR_1(j)
              EMATRIX_2(i,j)=EVECTOR_2(i)*EVECTOR_2(j)
              EMATRIX_3(i,j)=EVECTOR_3(i)*EVECTOR_3(j)
            END DO
          END DO

          CALL MatrixProduct(F_a_inv,EMATRIX_1,N1,err,error,*999)
          CALL MatrixProduct(N1,F_a_inv,N1,err,error,*999) ! F_a_inv=F_a_inv_T
          CALL MatrixProduct(F_a_inv,EMATRIX_2,N2,err,error,*999)
          CALL MatrixProduct(N2,F_a_inv,N2,err,error,*999) ! F_a_inv=F_a_inv_T
          CALL MatrixProduct(F_a_inv,EMATRIX_3,N3,err,error,*999)
          CALL MatrixProduct(N3,F_a_inv,N3,err,error,*999) ! F_a_inv=F_a_inv_T 

          FREE_ENERGY=0.0_DP
          DO i=1,3
            FREE_ENERGY=FREE_ENERGY+C(i)/C(i+3)*( &
              & EVALUES(1)**(C(i+3)/2.0_DP)+ &
              & EVALUES(2)**(C(i+3)/2.0_DP)+ &
              & EVALUES(3)**(C(i+3)/2.0_DP)-3.0_DP)
          END DO
          FREE_ENERGY=C(7)*FREE_ENERGY

          VALUE=XB_ENERGY_PER_VOLUME-(FREE_ENERGY-FREE_ENERGY_0)
        ENDIF
      ENDDO
    
      PIOLA_TENSOR=0.0_DP
      DO i=1,3
        PIOLA_TENSOR=PIOLA_TENSOR+ &
          & C(i)*EVALUES(1)**(C(i+3)/2.0_DP-1.0_DP)*N1+ &
          & C(i)*EVALUES(2)**(C(i+3)/2.0_DP-1.0_DP)*N2+ &
          & C(i)*EVALUES(3)**(C(i+3)/2.0_DP-1.0_DP)*N3
      END DO
      PIOLA_TENSOR=PIOLA_TENSOR*C(7)+2.0_DP*P*AZU
      
      !store lambda_f, so it can be used in the CellML file
      lambda_f=SQRT(AZL(1,1))
      CALL FIELD_PARAMETER_SET_UPDATE_GAUSS_POINT(DEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,GAUSS_POINT_NUMBER,&
        & ELEMENT_NUMBER,1,lambda_f,err,error,*999)

    CASE(EQUATIONS_SET_1D3D_MONODOMAIN_ACTIVE_STRAIN_SUBTYPE)

      NULLIFY(INDEPENDENT_FIELD)
      INDEPENDENT_FIELD=>EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD
      NULLIFY(FIELD_VARIABLE)
      CALL Field_VariableGet(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VARIABLE,err,error,*999)

      dof_idx=FIELD_VARIABLE%COMPONENTS(5)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(GAUSS_POINT_NUMBER, &
        & ELEMENT_NUMBER)
      CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,dof_idx,lambda_a, &
        & err,error,*999)

      F_a_inv=0.0_DP
      F_a_inv(1,1)=1.0_DP/lambda_a
      F_a_inv(2,2)=1.0_DP
      F_a_inv(3,3)=1.0_DP

      CALL MatrixProduct(DZDNU,F_a_inv,F_e,err,error,*999)
      CALL MatrixTranspose(F_e,F_e_T,err,error,*999)
      CALL MatrixProduct(F_e_T,F_e,C_e,err,error,*999)
      
      !Odgen law - 3 terms. Material Parameters C = [mu(1) mu(2) mu(3) alpha(1) alpha(2) alpha(3) mu_0]
!      CALL Eigenvalue(C_e,EVALUES,err,error,*999)
      CALL DSYEV('V','U',3,C_e,3,EVALUES,WORK,-1,ERR)
      IF(ERR.NE.0) CALL FlagError("Error in Eigenvalue computation",err,error,*999)
      LWORK=MIN(LWMAX,INT(WORK(1)))
      CALL DSYEV('V','U',3,C_e,3,EVALUES,WORK,LWORK,ERR)
      IF(ERR.NE.0) CALL FlagError("Error in Eigenvalue computation",err,error,*999)
      EVECTOR_1=C_e(:,1)
      EVECTOR_2=C_e(:,2)
      EVECTOR_3=C_e(:,3)

      DO i=1,3
        DO j=1,3
          EMATRIX_1(i,j)=EVECTOR_1(i)*EVECTOR_1(j)
          EMATRIX_2(i,j)=EVECTOR_2(i)*EVECTOR_2(j)
          EMATRIX_3(i,j)=EVECTOR_3(i)*EVECTOR_3(j)
        END DO
      END DO

      CALL MatrixProduct(F_a_inv,EMATRIX_1,N1,err,error,*999)
      CALL MatrixProduct(N1,F_a_inv,N1,err,error,*999) ! F_a_inv=F_a_inv_T
      CALL MatrixProduct(F_a_inv,EMATRIX_2,N2,err,error,*999)
      CALL MatrixProduct(N2,F_a_inv,N2,err,error,*999) ! F_a_inv=F_a_inv_T
      CALL MatrixProduct(F_a_inv,EMATRIX_3,N3,err,error,*999)
      CALL MatrixProduct(N3,F_a_inv,N3,err,error,*999) ! F_a_inv=F_a_inv_T 

      PIOLA_TENSOR=0.0_DP
      DO i=1,3
        PIOLA_TENSOR=PIOLA_TENSOR+ &
          & C(i)*EVALUES(1)**(C(i+3)/2.0_DP-1.0_DP)*N1+ &
          & C(i)*EVALUES(2)**(C(i+3)/2.0_DP-1.0_DP)*N2+ &
          & C(i)*EVALUES(3)**(C(i+3)/2.0_DP-1.0_DP)*N3
      END DO
      PIOLA_TENSOR=PIOLA_TENSOR*C(7)+2.0_DP*P*AZU
      
    CASE(EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE,EQUATIONS_SET_MOONEY_RIVLIN_ACTIVECONTRACTION_SUBTYPE,EQUATIONS_SET_MEMBRANE_SUBTYPE,&
      & EQUATIONS_SET_NO_SUBTYPE,EQUATIONS_SET_ANISOTROPIC_POLYNOMIAL_SUBTYPE,EQUATIONS_SET_ANISOTROPIC_POLYNOMIAL_ACTIVE_SUBTYPE,&
      & EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE,EQUATIONS_SET_STANDARD_MONODOMAIN_ELASTICITY_SUBTYPE, &
      & EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE,EQUATIONS_SET_TRANSVERSE_ISOTROPIC_POLYNOMIAL_SUBTYPE, &
      & EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE,EQUATIONS_SET_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE, &
      & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_ACTIVE_SUBTYPE,EQUATIONS_SET_MR_AND_GROWTH_LAW_IN_CELLML_SUBTYPE)
      !Form of constitutive model is:
      ! W=c1*(I1-3)+c2*(I2-3)+p*(I3-1)
      !Also assumed I3 = det(AZL) = 1.0
      !  Note that because PIOLA = 2.del{W}/del{C}=[...]+2.lambda.J^2.C^{-1}
      !  lambda here is actually half of hydrostatic pressure -- is this comment still correct?
      !If subtype is membrane, assume Mooney Rivlin constitutive law
      IF (EQUATIONS_SET_SUBTYPE/=EQUATIONS_SET_MEMBRANE_SUBTYPE) THEN
          PIOLA_TENSOR(1,3)=2.0_DP*(C(2)*(-AZL(3,1)))+P*AZU(1,3)
          PIOLA_TENSOR(2,3)=2.0_DP*(C(2)*(-AZL(3,2)))+P*AZU(2,3)
          PIOLA_TENSOR(3,1)=PIOLA_TENSOR(1,3)
          PIOLA_TENSOR(3,2)=PIOLA_TENSOR(2,3)
          PIOLA_TENSOR(3,3)=2.0_DP*(C(1)+C(2)*(AZL(1,1)+AZL(2,2)))+P*AZU(3,3)
      ELSE
        ! Membrane Equations
        ! Assume incompressible => I3 = 1 => C33(C11 x C22 - C12*C21) = 1
        AZL(3,3) = 1.0_DP / ((AZL(1,1) * AZL(2,2)) - (AZL(1,2) * AZL (2,1)))
        ! Assume Mooney-Rivlin constitutive relation
        P = -1.0_DP*((C(1) + C(2) * (AZL(1,1) + AZL(2,2))) * AZL(3,3))
        ! Assume stress normal to the surface is neglible i.e. PIOLA_TENSOR(:,3) = 0,PIOLA_TENSOR(3,:) = 0
        PIOLA_TENSOR(:,3) = 0.0_DP
        PIOLA_TENSOR(3,:) = 0.0_DP
      ENDIF
      PIOLA_TENSOR(1,1)=2.0_DP*(C(1)+C(2)*(AZL(2,2)+AZL(3,3)))+P*AZU(1,1)
      PIOLA_TENSOR(1,2)=2.0_DP*(     C(2)*(-AZL(2,1)))+P*AZU(1,2)
      PIOLA_TENSOR(2,1)=PIOLA_TENSOR(1,2)
      PIOLA_TENSOR(2,2)=2.0_DP*(C(1)+C(2)*(AZL(3,3)+AZL(1,1)))+P*AZU(2,2)


      SELECT CASE(EQUATIONS_SET_SUBTYPE)
      CASE(EQUATIONS_SET_MOONEY_RIVLIN_ACTIVECONTRACTION_SUBTYPE)
        !add active contraction stress value to the trace of the stress tensor - basically adding to hydrostatic pressure.
        !the active stress is stored inside the independent field that has been set up in the user program.
        !for generality we could set up 3 components in independent field for 3 different active stress components
        !1 isotropic value assumed here.
        CALL Field_ParameterSetGetLocalGaussPoint(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, & 
          &  FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,GAUSS_POINT_NUMBER,ELEMENT_NUMBER,1,ACTIVE_STRESS_11, &
          & err,error,*999) ! get the independent field stress value

        CALL Field_ParameterSetGetLocalGaussPoint(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, & 
          &  FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,GAUSS_POINT_NUMBER,ELEMENT_NUMBER,2,ACTIVE_STRESS_22, &
          & err,error,*999) ! get the independent field stress value

        CALL Field_ParameterSetGetLocalGaussPoint(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, & 
          &  FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,GAUSS_POINT_NUMBER,ELEMENT_NUMBER,3,ACTIVE_STRESS_33, &
          & err,error,*999) ! get the independent field stress value

        PIOLA_TENSOR(1,1)=PIOLA_TENSOR(1,1)+ACTIVE_STRESS_11
        PIOLA_TENSOR(2,2)=PIOLA_TENSOR(2,2)+ACTIVE_STRESS_22
        PIOLA_TENSOR(3,3)=PIOLA_TENSOR(3,3)+ACTIVE_STRESS_33

      CASE(EQUATIONS_SET_STANDARD_MONODOMAIN_ELASTICITY_SUBTYPE)
        ! add the active stress component (stored in the independent field) to the 1,1-direction of the 2-PK tensor
        PIOLA_TENSOR(1,1)=PIOLA_TENSOR(1,1)+INDEPENDENT_INTERPOLATED_POINT%VALUES(1,NO_PART_DERIV)

      CASE(EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE)
        !passive anisotropic stiffness -- only in the tension range
        IF(AZL(1,1) > 1.0_DP) THEN
          PIOLA_TENSOR(1,1)=PIOLA_TENSOR(1,1)+C(3)/AZL(1,1)*(AZL(1,1)**(C(4)/2.0_DP)-1.0_DP)
        ENDIF
        !active stress component
        CALL Field_ParameterSetGetLocalGaussPoint(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, & 
          & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,GAUSS_POINT_NUMBER,ELEMENT_NUMBER,1,VALUE, &
          & err,error,*999)
        !divide by lambda and multiply by P_max
        VALUE=VALUE/SQRT(AZL(1,1))*C(5)

        !HINDAWI paper - force-length relation at the continuum level
!        if((SQRT(AZL(1,1))>0.72_DP).AND.(SQRT(AZL(1,1))<1.68_DP)) then
!          VALUE=VALUE*(-25.0_DP/4.0_DP*AZL(1,1)/1.2_DP/1.2_DP + 25.0_DP/2.0_DP*SQRT(AZL(1,1))/1.2_DP - 5.25_DP)
!        else
!          VALUE=0.0_DP
!        endif

        PIOLA_TENSOR(1,1)=PIOLA_TENSOR(1,1)+VALUE

      CASE(EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE)
        !passive anisotropic stiffness -- only in the tension range
        IF(AZL(1,1) > 1.0_DP) THEN
          PIOLA_TENSOR(1,1)=PIOLA_TENSOR(1,1)+C(3)/AZL(1,1)*(AZL(1,1)**(C(4)/2.0_DP)-1.0_DP)
        ENDIF
        !active stress component
        CALL Field_VariableGet(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VARIABLE,err,error,*999)
        dof_idx=FIELD_VARIABLE%COMPONENTS(1)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(GAUSS_POINT_NUMBER, &
          & ELEMENT_NUMBER)
        CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
          & FIELD_VALUES_SET_TYPE,dof_idx,VALUE,err,error,*999)

        IF(VALUE.LT.0.0_DP) VALUE=0.0_DP

        !divide by lambda and multiply by P_max
        VALUE=VALUE/SQRT(AZL(1,1))*C(5)

        PIOLA_TENSOR(1,1)=PIOLA_TENSOR(1,1)+VALUE

        ! unbound Titin-stress
        dof_idx=FIELD_VARIABLE%COMPONENTS(2)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(GAUSS_POINT_NUMBER, &
          & ELEMENT_NUMBER)
        CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,dof_idx,TITIN_UNBOUND,err,error,*999)
        ! bound Titin-stress -> Rode Model
        dof_idx=FIELD_VARIABLE%COMPONENTS(3)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(GAUSS_POINT_NUMBER, &
          & ELEMENT_NUMBER)
        CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
          & FIELD_VALUES_SET_TYPE,dof_idx,TITIN_BOUND,err,error,*999)
        ! activation
        dof_idx=FIELD_VARIABLE%COMPONENTS(6)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(GAUSS_POINT_NUMBER, &
          & ELEMENT_NUMBER)
        CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
          & FIELD_VALUES_SET_TYPE,dof_idx,activation,err,error,*999)

        IF(activation.GT.1.0_DP) activation=1.0_DP
        IF(activation.LT.0.0_DP) activation=0.0_DP

        ! parameter to switch on and off actin-titin interaction
        activation=C(6)*activation
        
        ! normalized Titin-stress -> weighted sum of bound and unbound titin-stress
        TITIN_VALUE=activation*TITIN_BOUND+(1.0_DP-activation)*TITIN_UNBOUND
        !TITIN_VALUE=activation*TITIN_BOUND*0.5_DP+(1.0_DP-activation)*TITIN_UNBOUND !TK Hack
        ! divide by lambda and multiply by P_max
        TITIN_VALUE=TITIN_VALUE/SQRT(AZL(1,1))*C(5)

        PIOLA_TENSOR(1,1)=PIOLA_TENSOR(1,1)+TITIN_VALUE

        ! unbound titin-stress in cross-fibre direction
        dof_idx=FIELD_VARIABLE%COMPONENTS(4)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(GAUSS_POINT_NUMBER, &
          & ELEMENT_NUMBER)
        CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
          & FIELD_VALUES_SET_TYPE,dof_idx,TITIN_UNBOUND_CROSS_FIBRE,err,error,*999)
        ! bound titin-stress in cross-fibre direction
        dof_idx=FIELD_VARIABLE%COMPONENTS(5)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(GAUSS_POINT_NUMBER, &
          & ELEMENT_NUMBER)
        CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
          & FIELD_VALUES_SET_TYPE,dof_idx,TITIN_BOUND_CROSS_FIBRE,err,error,*999)

        ! normalized XF-Titin-stress -> weighted sum of bound and unbound XF-titin-stress
        TITIN_VALUE_CROSS_FIBRE=activation*TITIN_BOUND_CROSS_FIBRE+(1.0_DP-activation)*TITIN_UNBOUND_CROSS_FIBRE
        ! divide by lambda and multiply by P_max
        TITIN_VALUE_CROSS_FIBRE=TITIN_VALUE_CROSS_FIBRE*C(5) !/SQRT(AZL(1,1))
 
        PIOLA_TENSOR(2,2)=PIOLA_TENSOR(2,2)+TITIN_VALUE_CROSS_FIBRE
        PIOLA_TENSOR(3,3)=PIOLA_TENSOR(3,3)+TITIN_VALUE_CROSS_FIBRE

      CASE(EQUATIONS_SET_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE)
        !passive anisotropic stiffness -- only in the tension range
        IF(AZL(1,1) > 1.0_DP) THEN
!tomo
!          PIOLA_TENSOR(1,1)=PIOLA_TENSOR(1,1)+C(3)/AZL(1,1)*(AZL(1,1)**(C(4)/2.0_DP)-1.0_DP)
          PIOLA_TENSOR(1,1)=PIOLA_TENSOR(1,1)+0.355439810963035_DP/AZL(1,1)*(AZL(1,1)**(12.660539325481963_DP/2.0_DP)-1.0_DP)
        ENDIF
!tomo
        IF(AZL(2,2) > 1.0_DP) THEN
          PIOLA_TENSOR(2,2)=PIOLA_TENSOR(2,2)+5316.372204148964_DP/AZL(2,2)*(AZL(2,2)**(0.014991843974911_DP/2.0_DP)-1.0_DP)
        ENDIF
        IF(AZL(3,3) > 1.0_DP) THEN
          PIOLA_TENSOR(3,3)=PIOLA_TENSOR(3,3)+5316.372204148964_DP/AZL(3,3)*(AZL(3,3)**(0.014991843974911_DP/2.0_DP)-1.0_DP)
        ENDIF
!tomo end
        !active stress component
        CALL Field_VariableGet(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VARIABLE,err,error,*999)
        dof_idx=FIELD_VARIABLE%COMPONENTS(1)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(GAUSS_POINT_NUMBER, &
          & ELEMENT_NUMBER)
        CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
          & FIELD_VALUES_SET_TYPE,dof_idx,VALUE,err,error,*999)
        !divide by lambda and multiply by P_max
!tomo RFE
        VAL1=VALUE
!tomo REF end
        VALUE=VALUE/SQRT(AZL(1,1))*C(5)


!tomo RFE
        !alpha*K_rfe*(lambda-lambda_start)/lambda
        !TODO make lambda_start variable --> independent field
!        VAL2=VAL1*100.0_DP*(SQRT(AZL(1,1))-1) !stretch and compression!!!
        VAL2=100.0_DP*(SQRT(AZL(1,1))-1) !stretch and compression!!!
        VALUE=VALUE+VAL2/SQRT(AZL(1,1))
        PIOLA_TENSOR(1,1)=PIOLA_TENSOR(1,1)+VALUE
!tomo REF end

!        PIOLA_TENSOR(1,1)=PIOLA_TENSOR(1,1)+VALUE

      CASE(EQUATIONS_SET_TRANSVERSE_ISOTROPIC_POLYNOMIAL_SUBTYPE)
        !Additional term for transversely isotropic (fibre-reinforced) materials (Markert, B., W. Ehlers, and N. Karajan.
        !A general polyconvex strain-energy function for fiber-reinforced materials.
        !Proceedings in Applied Mathematics and Mechanics 5.1 (2005): 245-246.)

        ! W_aniso=c3*(sqrt(I4)^(c4-2)-1/I4)M
        ! with M being the mapping towards the fibre direction, here: I4=C_11
        !C(3)=c3...polynomial coefficient
        !C(4)=c4...power coefficient
        IF(AZL(1,1) > 1.0_DP) THEN ! only in the tension range
          PIOLA_TENSOR(1,1)=PIOLA_TENSOR(1,1)+C(3)/AZL(1,1)*(AZL(1,1)**(C(4)/2.0_DP)-1.0_DP)
        ENDIF

      CASE(EQUATIONS_SET_TRANSVERSE_ISOTROPIC_ACTIVE_SUBTYPE)
        !Isotropic and anisotropic part from above, additionally an active part in fibre direction
        ! W=W_iso+W_aniso+W_act
        !  with W_act=(1/sqrt(I4)*P_max*f*alpha)M
        !C(5)=alpha...activation parameter [0,1]
        IF(AZL(1,1) > 1.0_DP) THEN ! only in the tension range
          PIOLA_TENSOR(1,1)=PIOLA_TENSOR(1,1)+C(3)/AZL(1,1)*(AZL(1,1)**(C(4)/2.0_DP)-1.0_DP)
        ENDIF
!        IF((SQRT(AZL(1,1))>0.84_DP).AND.(SQRT(AZL(1,1))<1.96_DP)) THEN
        if((SQRT(AZL(1,1))>0.72_DP).AND.(SQRT(AZL(1,1))<1.68_DP)) then
!          VALUE=(-25.0_DP/4.0_DP*AZL(1,1)/1.4_DP/1.4_DP + 25.0_DP/2.0_DP*SQRT(AZL(1,1))/1.4_DP - 5.25_DP) !f
          VALUE=(-25.0_DP/4.0_DP*AZL(1,1)/1.2_DP/1.2_DP + 25.0_DP/2.0_DP*SQRT(AZL(1,1))/1.2_DP - 5.25_DP)
          VALUE=VALUE*(1.0_DP/SQRT(AZL(1,1)))*20.0_DP*C(5)
          PIOLA_TENSOR(1,1)=PIOLA_TENSOR(1,1)+VALUE
        ENDIF

      CASE(EQUATIONS_SET_ANISOTROPIC_POLYNOMIAL_SUBTYPE)
        !Three additional terms for transversely isotropic (fibre-reinforced) materials (Markert, B., W. Ehlers, and N. Karajan.
        !A general polyconvex strain-energy function for fiber-reinforced materials.
        !Proceedings in Applied Mathematics and Mechanics 5.1 (2005): 245-246.)
        ! W_aniso=c3*(sqrt(I4)^(c4-2)-1/I4)M_1 + c5*(sqrt(I4)^(c6-2)-1/I4)M_2 + c7*(sqrt(I4)^(c8-2)-1/I4)M_3
        ! with M_1 being the mapping towards the fibre direction, here: I4=C_11
        !C(3)=c3...polynomial coefficient
        !C(4)=c4...power coefficient
        !C(5)=c5...polynomial coefficient
        !C(6)=c6...power coefficient
        !C(7)=c7...polynomial coefficient
        !C(8)=c8...power coefficient
        IF(AZL(1,1) > 1.0_DP) THEN ! only in the tension range
          PIOLA_TENSOR(1,1)=PIOLA_TENSOR(1,1)+C(3)/AZL(1,1)*(AZL(1,1)**(C(4)/2.0_DP)-1.0_DP)
        ENDIF
        IF(AZL(2,2) > 1.0_DP) THEN
          PIOLA_TENSOR(2,2)=PIOLA_TENSOR(2,2)+C(5)/AZL(2,2)*(AZL(2,2)**(C(6)/2.0_DP)-1.0_DP)
        ENDIF
        IF(AZL(3,3) > 1.0_DP) THEN
          PIOLA_TENSOR(3,3)=PIOLA_TENSOR(3,3)+C(7)/AZL(3,3)*(AZL(3,3)**(C(8)/2.0_DP)-1.0_DP)
        ENDIF

      CASE(EQUATIONS_SET_ANISOTROPIC_POLYNOMIAL_ACTIVE_SUBTYPE)
        !Three additional terms for transversely isotropic (fibre-reinforced) materials (Markert, B., W. Ehlers, and N. Karajan.
        !A general polyconvex strain-energy function for fiber-reinforced materials.
        !Proceedings in Applied Mathematics and Mechanics 5.1 (2005): 245-246.)
        ! W_aniso=c3*(sqrt(I4)^(c4-2)-1/I4)M_1 + c5*(sqrt(I4)^(c6-2)-1/I4)M_2 + c7*(sqrt(I4)^(c8-2)-1/I4)M_3
        ! with M_1 being the mapping towards the fibre direction, here: I4=C_11
        !C(3)=c3...polynomial coefficient
        !C(4)=c4...power coefficient
        !C(5)=c5...polynomial coefficient
        !C(6)=c6...power coefficient
        !C(7)=c7...polynomial coefficient
        !C(8)=c8...power coefficient
        !C(9)=lambda_opt...optimal fibre stretch
        !C(10)=P_max...maximum active tension
        !C(11)=alpha...activation parameter [0 1]
        !C(12)=K_rfe...stiffness of the residual force enhancement
        IF(AZL(1,1) > 1.0_DP) THEN ! only in the tension range
          PIOLA_TENSOR(1,1)=PIOLA_TENSOR(1,1)+C(3)/AZL(1,1)*(AZL(1,1)**(C(4)/2.0_DP)-1.0_DP)
        ENDIF
        IF(AZL(2,2) > 1.0_DP) THEN
          PIOLA_TENSOR(2,2)=PIOLA_TENSOR(2,2)+C(5)/AZL(2,2)*(AZL(2,2)**(C(6)/2.0_DP)-1.0_DP)
        ENDIF
        IF(AZL(3,3) > 1.0_DP) THEN
          PIOLA_TENSOR(3,3)=PIOLA_TENSOR(3,3)+C(7)/AZL(3,3)*(AZL(3,3)**(C(8)/2.0_DP)-1.0_DP)
        ENDIF

        VAL1=SQRT(AZL(1,1))/C(9) !lambda/lambda_opt
        IF((VAL1>0.7_DP).AND.(VAL1<1.3_DP)) THEN
          !active force-length relation
          VALUE=(-11.1111_DP*VAL1*VAL1+22.2222_DP*VAL1-10.1111_DP)
          !multiply by P_max and alpha, divide by lambda
          VALUE=VALUE*C(10)*C(11)/SQRT(AZL(1,1))
        ELSE
          VALUE=0.0_DP
        ENDIF
        !alpha*K_rfe*(lambda-lambda_start)/lambda
        !TODO make lambda_start variable --> independent field
        VAL2=C(11)*C(12)*(SQRT(AZL(1,1))-1) !stretch and compression!!!
        VALUE=VALUE+VAL2/SQRT(AZL(1,1))
        PIOLA_TENSOR(1,1)=PIOLA_TENSOR(1,1)+VALUE

      END SELECT


    CASE(EQUATIONS_SET_TRANS_ISOTROPIC_ACTIVE_TRANSITION_SUBTYPE)
      !Equations set for transversely isotropic (fibre-reinforced), active contractible bodies consitisting of two materials
      ! The local portion between them is defined by the parameter trans
      ! Material 1 is active contractible, material 2 is only passive
      !W=W_iso+W_aniso+W_act
      ! where the three parts are adopted from above (iso Mooney-Rivlin, aniso Markert, active part)
      !Markert, B., W. Ehlers, and N. Karajan.
      !A general polyconvex strain-energy function for fiber-reinforced materials.
      !Proceedings in Applied Mathematics and Mechanics 5.1 (2005): 245-246.)

      !C(1)=c1_m1...Mooney Rivlin parameter material 1
      !C(2)=c2_m1...Mooney Rivlin parameter material 1
      !C(3)=c4_m1...polynomial coefficient (Markert model) material 1
      !C(4)=c5_m1...power coefficient (Markert model) material 1
      !C(5)=c1_m2...Mooney Rivlin parameter material 2
      !C(6)=c2_m2...Mooney Rivlin parameter material 2
      !C(7)=c4_m2...polynomial coefficient (Markert model) material 2
      !C(8)=c5_m2...power coefficient (Markert model) material 2
      !C(9)=alpha...activation parameter [0,1]
      !C(10)=trans...transition parameter [0,1] for the portion between the two materials
      !C(11)=P_max...maximum isometric stress

      !Weighting the Mooney Rivlin parameters and obtaining resulting c1 and c2
      VAL1=C(1)*C(10)+C(5)*(1.0_DP-C(10))
      VAL2=C(2)*C(10)+C(6)*(1.0_DP-C(10))

      !Mooney-Rivlin for the isotropic part
      PIOLA_TENSOR(1,1)=2.0_DP*(VAL1+VAL2*(AZL(2,2)+AZL(3,3))+P*AZU(1,1))
      PIOLA_TENSOR(1,2)=2.0_DP*(     VAL2*(-AZL(2,1))        +P*AZU(1,2))
      PIOLA_TENSOR(1,3)=2.0_DP*(     VAL2*(-AZL(3,1))        +P*AZU(1,3))
      PIOLA_TENSOR(2,1)=PIOLA_TENSOR(1,2)
      PIOLA_TENSOR(2,2)=2.0_DP*(VAL1+VAL2*(AZL(3,3)+AZL(1,1))+P*AZU(2,2))
      PIOLA_TENSOR(2,3)=2.0_DP*(     VAL2*(-AZL(3,2))        +P*AZU(2,3))
      PIOLA_TENSOR(3,1)=PIOLA_TENSOR(1,3)
      PIOLA_TENSOR(3,2)=PIOLA_TENSOR(2,3)
      PIOLA_TENSOR(3,3)=2.0_DP*(VAL1+VAL2*(AZL(1,1)+AZL(2,2))+P*AZU(3,3))

      !passive anisotropic part -- only in the tension range (Markert)
      IF(AZL(1,1) > 1.0_DP) THEN
        VAL1=C(3)/AZL(1,1)*(AZL(1,1)**(C(4)/2.0_DP)-1.0_DP)
        VAL2=C(7)/AZL(1,1)*(AZL(1,1)**(C(8)/2.0_DP)-1.0_DP)
        PIOLA_TENSOR(1,1)=PIOLA_TENSOR(1,1)+(VAL1*C(10)+VAL2*(1.0_DP-C(10)))
      ENDIF

      !active part
      IF((SQRT(AZL(1,1))>0.84_DP).AND.(SQRT(AZL(1,1))<1.96_DP)) THEN
        VALUE=(-25.0_DP/4.0_DP*AZL(1,1)/1.4_DP/1.4_DP + 25.0_DP/2.0_DP*SQRT(AZL(1,1))/1.4_DP - 5.25_DP)
        VALUE=VALUE*(1.0_DP/SQRT(AZL(1,1)))*C(9)*C(10)*C(11)
        PIOLA_TENSOR(1,1)=PIOLA_TENSOR(1,1)+VALUE
      ENDIF

    CASE(EQUATIONS_SET_ISOTROPIC_EXPONENTIAL_SUBTYPE)
      !Form of constitutive model is:
      ! W=c1/2 (e^(c2*(I1-3)) - 1)
      ! S = 2*dW/dC + 2pC^-1
      PIOLA_TENSOR=C(1)*C(2)*EXP(C(2)*(AZL(1,1)+AZL(2,2)+AZL(3,3)-3.0_DP))*IDENTITY+2.0_DP*P*AZU
    CASE(EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_STATIC_INRIA_SUBTYPE)
      !C(1)=Mooney Rivlin parameter
      !C(2)=Mooney Rivlin parameter
      !C(3)=K
      !C(4)=M, Biot modulus
      !C(5)=b, skeleton parameter
      !C(6)=p0, reference pressure

      P=DARCY_DEPENDENT_INTERPOLATED_POINT%VALUES(1,NO_PART_DERIV) !Fluid pressure
      CALL MatrixTranspose(AZL,AZLT,err,error,*999)
      I1=AZL(1,1)+AZL(2,2)+AZL(3,3)
      TEMP=MATMUL(AZL,AZL)
      I2=0.5_DP*(I1**2.0_DP-TEMP(1,1)-TEMP(2,2)-TEMP(3,3))

      CALL EVALUATE_CHAPELLE_FUNCTION(Jznu,ffact,dfdJfact,err,error,*999)

      PIOLA_TENSOR=2.0_DP*C(1)*Jznu**(-2.0_DP/3.0_DP)*(IDENTITY-(1.0_DP/3.0_DP)*I1*AZU)
      PIOLA_TENSOR=PIOLA_TENSOR+2.0_DP*C(2)*Jznu**(-4.0_DP/3.0_DP)*(I1*IDENTITY-AZLT-(2.0_DP/3.0_DP)*I2*AZU)
      PIOLA_TENSOR=PIOLA_TENSOR+(C(3)-C(4)*C(5)**2)*(Jznu-1.0_DP)*AZU
      PIOLA_TENSOR=PIOLA_TENSOR-C(5)*(P-C(6))*Jznu*AZU
      PIOLA_TENSOR=PIOLA_TENSOR+0.5_DP*((P-C(6))**2/C(4))*(dfdJfact/(ffact**2))*Jznu*AZU
    CASE(EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_HOLMES_MOW_SUBTYPE)
      ! See Holmes MH, Mow VC. The nonlinear characteristics of soft gels and hydrated connective tissues in ultrafiltration.
      ! Journal of Biomechanics. 1990;23(11):1145-1156. DOI: 10.1016/0021-9290(90)90007-P
      ! The form of constitutive relation is:
      ! sigma = sigma^s + sigma^f
      ! sigma^f = -phi^f p I
      ! sigma^s = -phi^s p I + rho_0^s sigma^s_E
      ! sigma^s_E is the effective Cauchy stress obtained by differentiating
      ! the free energy function to get the second Piola-Kirchoff stress tensor:
      ! rho_0^s W^s = c0 exp(c1(I1 - 3) + c2(I2 - 3)) / (I_3^(c1 + 2c2))
      ! Rather than add the "phi^s p I" term to the Cauchy stress, we add it here as "phi^s p J C^-1"
      ! We also set rho_0^s = the solid density * initial solidity, and move the solidity
      ! inside the strain energy density function
      !
      ! c0 = C(1)
      ! c1 = C(2)
      ! c2 = C(3)
      ! phi^s_0 = C(4)

      CALL MatrixTranspose(AZL,AZLT,err,error,*999)
      CALL MatrixTranspose(AZU,AZUT,err,error,*999)
      I1=AZL(1,1)+AZL(2,2)+AZL(3,3)
      TEMP=MATMUL(AZL,AZL)
      I2=0.5_DP*(I1**2.0_DP-TEMP(1,1)-TEMP(2,2)-TEMP(3,3))
      !I3 already defined

      TEMPTERM=2.0_DP*C(4)*C(1)*EXP(C(2)*(I1 - 3.0_DP) + C(3)*(I2 - 3.0_DP)) / (I3**(C(2)+2.0_DP*C(3)))
      PIOLA_TENSOR=C(2)*TEMPTERM*IDENTITY + C(3)*TEMPTERM*(I1*IDENTITY-AZLT) - (C(2)+2.0_DP*C(3))*TEMPTERM*AZUT
      PIOLA_TENSOR=PIOLA_TENSOR - DARCY_DEPENDENT_INTERPOLATED_POINT%VALUES(1,NO_PART_DERIV)*Jznu*AZU

    CASE(EQUATIONS_SET_ELASTICITY_FLUID_PRES_HOLMES_MOW_ACTIVE_SUBTYPE)
      ! See Holmes MH, Mow VC. The nonlinear characteristics of soft gels and hydrated connective tissues in ultrafiltration.
      ! Journal of Biomechanics. 1990;23(11):1145-1156. DOI: 10.1016/0021-9290(90)90007-P
      ! The form of constitutive relation is:
      ! sigma = sigma^s + sigma^f
      ! sigma^f = -phi^f p I
      ! sigma^s = -phi^s p I + rho_0^s sigma^s_E
      ! sigma^s_E is the effective Cauchy stress obtained by differentiating
      ! the free energy function to get the second Piola-Kirchoff stress tensor:
      ! rho_0^s W^s = c0 exp(c1(I1 - 3) + c2(I2 - 3)) / (I_3^(c1 + 2c2))
      ! Rather than add the "phi^s p I" term to the Cauchy stress, we add it here as "phi^s p J C^-1"
      ! We also set rho_0^s = the solid density * initial solidity, and move the solidity
      ! inside the strain energy density function
      !
      ! c0 = C(1)
      ! c1 = C(2)
      ! c2 = C(3)
      ! phi^s_0 = C(4)
      ! alpha = C(5) (activation level)
      ! P_max = C(6) (maximum isometric active stress)

      CALL MatrixTranspose(AZL,AZLT,err,error,*999)
      CALL MatrixTranspose(AZU,AZUT,err,error,*999)
      I1=AZL(1,1)+AZL(2,2)+AZL(3,3)
      TEMP=MATMUL(AZL,AZL)
      I2=0.5_DP*(I1**2.0_DP-TEMP(1,1)-TEMP(2,2)-TEMP(3,3))
      !I3 already defined

      TEMPTERM=2.0_DP*C(4)*C(1)*EXP(C(2)*(I1 - 3.0_DP) + C(3)*(I2 - 3.0_DP)) / (I3**(C(2)+2.0_DP*C(3)))
      PIOLA_TENSOR=C(2)*TEMPTERM*IDENTITY + C(3)*TEMPTERM*(I1*IDENTITY-AZLT) - (C(2)+2.0_DP*C(3))*TEMPTERM*AZUT
      PIOLA_TENSOR=PIOLA_TENSOR - DARCY_DEPENDENT_INTERPOLATED_POINT%VALUES(1,NO_PART_DERIV)*Jznu*AZU

      IF((SQRT(AZL(1,1))>0.72_DP).AND.(SQRT(AZL(1,1))<1.68_DP)) THEN
        VALUE=(-25.0_DP/4.0_DP*AZL(1,1)/1.2_DP/1.2_DP + 25.0_DP/2.0_DP*SQRT(AZL(1,1))/1.2_DP - 5.25_DP)
      ELSE
        VALUE=0.0_DP
      END IF

      PIOLA_TENSOR(1,1) = PIOLA_TENSOR(1,1) + 1.0_DP/SQRT(AZL(1,1))*C(5)*C(6)*VALUE

    CASE(EQUATIONS_SET_STVENANT_KIRCHOFF_ACTIVECONTRACTION_SUBTYPE)
    ! For of constitutive model is:
    ! W = 0.5lambda*tr(E)^2 + mu*tr(E^2)
    ! S = dW/dE = lambda*tr(E)Identity + 2muE
      PIOLA_TENSOR(1,3)=(2.0_DP*C(2)*E(1,3))+(2.0_DP*P*AZU(1,3))
      PIOLA_TENSOR(2,3)=(2.0_DP*C(2)*E(2,3))+(2.0_DP*P*AZU(2,3))
      PIOLA_TENSOR(3,1)=PIOLA_TENSOR(1,3)
      PIOLA_TENSOR(3,2)=PIOLA_TENSOR(2,3)
      PIOLA_TENSOR(3,3)=C(1)*(E(1,1)+E(2,2)+E(3,3))+(2.0_DP*E(3,3)*C(2)+(2.0_DP*P*AZU(3,3)))

      PIOLA_TENSOR(1,1)=C(1)*(E(1,1)+E(2,2)+E(3,3))+(2.0_DP*E(1,1)*C(2)+(2.0_DP*P*AZU(1,1)))
      PIOLA_TENSOR(1,2)=(2.0_DP*C(2)*E(1,2))+(2.0_DP*P*AZU(1,2))
      PIOLA_TENSOR(2,1)=PIOLA_TENSOR(1,2)
      PIOLA_TENSOR(2,2)=C(1)*(E(1,1)+E(2,2)+E(3,3))+(2.0_DP*E(2,2)*C(2)+(2.0_DP*P*AZU(2,2)))

      CALL Field_ParameterSetGetLocalGaussPoint(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
        &  FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,GAUSS_POINT_NUMBER,ELEMENT_NUMBER,1,ACTIVE_STRESS_11, &
        & err,error,*999) ! get the independent field stress value

      CALL Field_ParameterSetGetLocalGaussPoint(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
        &  FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,GAUSS_POINT_NUMBER,ELEMENT_NUMBER,2,ACTIVE_STRESS_22, &
        & err,error,*999) ! get the independent field stress value

      CALL Field_ParameterSetGetLocalGaussPoint(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
        &  FIELD_U_VARIABLE_TYPE, FIELD_VALUES_SET_TYPE,GAUSS_POINT_NUMBER,ELEMENT_NUMBER,3,ACTIVE_STRESS_33, &
        & err,error,*999) ! get the independent field stress value

      PIOLA_TENSOR(1,1)=PIOLA_TENSOR(1,1)+ACTIVE_STRESS_11
      PIOLA_TENSOR(2,2)=PIOLA_TENSOR(2,2)+ACTIVE_STRESS_22
      PIOLA_TENSOR(3,3)=PIOLA_TENSOR(3,3)+ACTIVE_STRESS_33

    CASE(EQUATIONS_SET_TRANSVERSE_ISOTROPIC_EXPONENTIAL_SUBTYPE)
      !Form of constitutive model is:
      ! W=c1/2 (e^Q - 1)
      ! where Q=2c2(E11+E22+E33)+c3(E11^2)+c4(E22^2+E33^2+E23^2+E32^2)+c5(E12^2+E21^2+E31^2+E13^2)
      ! with E expressed in fibre coordinates

      TEMPTERM=C(1)*EXP(2.0*C(2)*(E(1,1)+E(2,2)+E(3,3))+C(3)*E(1,1)**2+C(4)*(E(2,2)**2+E(3,3)**2+2.0_DP*E(2,3)**2)+ &
          & C(5)*2.0_DP*(E(1,2)**2+E(1,3)**2))
      PIOLA_TENSOR(1,1)=(C(2)+C(3)*E(1,1))*TEMPTERM+2.0_DP*P*AZU(1,1)
      PIOLA_TENSOR(1,2)=C(5)*E(1,2)*TEMPTERM+2.0_DP*P*AZU(1,2)
      PIOLA_TENSOR(1,3)=C(5)*E(1,3)*TEMPTERM+2.0_DP*P*AZU(1,3)
      PIOLA_TENSOR(2,1)=PIOLA_TENSOR(1,2)
      PIOLA_TENSOR(2,2)=(C(2)+C(4)*E(2,2))*TEMPTERM+2.0_DP*P*AZU(2,2)
      PIOLA_TENSOR(2,3)=C(4)*E(2,3)*TEMPTERM+2.0_DP*P*AZU(2,3)
      PIOLA_TENSOR(3,1)=PIOLA_TENSOR(1,3)
      PIOLA_TENSOR(3,2)=PIOLA_TENSOR(2,3)
      PIOLA_TENSOR(3,3)=(C(2)+C(4)*E(3,3))*TEMPTERM+2.0_DP*P*AZU(3,3)
    CASE(EQUATIONS_SET_TRANSVERSE_ISOTROPIC_GUCCIONE_SUBTYPE,EQUATIONS_SET_GUCCIONE_ACTIVECONTRACTION_SUBTYPE, &
      & EQUATIONS_SET_REFERENCE_STATE_TRANSVERSE_GUCCIONE_SUBTYPE)
      ! W=C1/2*exp*(Q) + p(J-1)
      ! Q=C2*E(1,1)^2 + C3*(E(2,2)^2+E(3,3)^2+2*E(2,3)*E(3,2)) + 2*C4*(E(1,2)*E(2,1)+E(1,3)*E(3,1))
      Q=C(2)*E(1,1)**2 + C(3)*(E(2,2)**2+E(3,3)**2+2.0_DP*E(2,3)**2) + 2.0_DP*C(4)*(E(1,2)**2+E(1,3)**2)
      TEMPTERM=0.5_DP*C(1)*exp(Q) ! iso term
      PIOLA_TENSOR(1,1) = 2.0_DP*C(2) * E(1,1)
      PIOLA_TENSOR(2,2) = 2.0_DP*C(3) * E(2,2)
      PIOLA_TENSOR(3,3) = 2.0_DP*C(3) * E(3,3)
      PIOLA_TENSOR(1,2) = 2.0_DP*C(4) * E(1,2)
      PIOLA_TENSOR(2,1) = PIOLA_TENSOR(1,2)
      PIOLA_TENSOR(1,3) = 2.0_DP*C(4) * E(1,3)
      PIOLA_TENSOR(3,1) = PIOLA_TENSOR(1,3)
      PIOLA_TENSOR(3,2) = 2.0_DP*C(3) * E(2,3)
      PIOLA_TENSOR(2,3) = PIOLA_TENSOR(3,2)
      PIOLA_TENSOR = PIOLA_TENSOR * TEMPTERM
      ! pressure terms
!
! TEMP DURING MERGE
!
!      PIOLA_TENSOR = PIOLA_TENSOR + 2.0_DP*p*Jznu*AZU   ! is Jznu required here, or is it omitted everywhere else?
!
!      IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_GUCCIONE_ACTIVECONTRACTION_SUBTYPE) THEN
!        !the active stress is stored inside the independent field that has been set up in the user program.
!        !for better generality we could set up 3 components in independent field for 3 different active stress components,
!        !but only one component is implemented so far for fibre active tension.
!        CALL Field_VariableGet(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VARIABLE,err,error,*999)
!        DO i=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
!          dof_idx=FIELD_VARIABLE%COMPONENTS(i)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP% &
!            & GAUSS_POINTS(GAUSS_POINT_NUMBER,ELEMENT_NUMBER)
!          CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
!            & FIELD_VALUES_SET_TYPE,dof_idx,VALUE,err,error,*999)
!          PIOLA_TENSOR(1,1)=PIOLA_TENSOR(1,1)+VALUE
!        ENDDO
!      ENDIF
      !PIOLA_TENSOR = PIOLA_TENSOR + 2.0_DP*p*Jznu*AZU   ! is Jznu required here, or is it omitted everywhere else?
      IF(EQUATIONS_SET%specification(3)==EQUATIONS_SET_REFERENCE_STATE_TRANSVERSE_GUCCIONE_SUBTYPE) THEN
        PRESSURE_COMPONENT=GEOMETRIC_INTERPOLATED_POINT%INTERPOLATION_PARAMETERS%FIELD_VARIABLE%NUMBER_OF_COMPONENTS
        P=GEOMETRIC_INTERPOLATED_POINT%VALUES(PRESSURE_COMPONENT,NO_PART_DERIV)
      ELSE
        PRESSURE_COMPONENT=DEPENDENT_INTERPOLATED_POINT%INTERPOLATION_PARAMETERS%FIELD_VARIABLE%NUMBER_OF_COMPONENTS
        P=DEPENDENT_INTERPOLATED_POINT%VALUES(PRESSURE_COMPONENT,NO_PART_DERIV)
      ENDIF 
      PIOLA_TENSOR = PIOLA_TENSOR + P*AZU   ! is Jznu required here, or is it omitted everywhere else?
      IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_GUCCIONE_ACTIVECONTRACTION_SUBTYPE) THEN
      !add active contraction stress value to the trace of the stress tensor - basically adding to hydrostatic pressure.
      !the active stress is stored inside the independent field that has been set up in the user program.
      !for generality we could set up 3 components in independent field for 3 different active stress components
        CALL Field_VariableGet(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VARIABLE,err,error,*999)
        DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
          dof_idx=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP% &
            & GAUSS_POINTS(GAUSS_POINT_NUMBER,ELEMENT_NUMBER)
          CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,dof_idx,VALUE,err,error,*999)
          PIOLA_TENSOR(component_idx,component_idx)=PIOLA_TENSOR(component_idx,component_idx)+VALUE
        ENDDO
      ENDIF
    CASE(EQUATIONS_SET_TRANSVERSE_ISOTROPIC_HUMPHREY_YIN_SUBTYPE)
      ! W=a*(exp(b(I1-3))-1) + c*(exp(d(alpha-1)^2)-1)
      ! a=C(1), b=C(2), c=C(3), d=C(4)
      I1=AZL(1,1)+AZL(2,2)+AZL(3,3)
      PIOLA_TENSOR(1,1)=C(1)*C(2)*EXP(C(2)*(I1-3))+ &
        & C(3)*2.0_DP*(SQRT(AZL(1,1))-1)*C(4)*EXP(C(4)*(SQRT(AZL(1,1))-1)**2)/(2*SQRT(AZL(1,1)))+P*AZU(1,1)
      PIOLA_TENSOR(2,2)=C(1)*C(2)*EXP(C(2)*(I1-3))+P*AZU(2,2)
      PIOLA_TENSOR(3,3)=C(1)*C(2)*EXP(C(2)*(I1-3))+P*AZU(3,3)
      PIOLA_TENSOR(1,2)=P*AZU(1,2)
      PIOLA_TENSOR(1,3)=P*AZU(1,3)
      PIOLA_TENSOR(2,3)=P*AZU(2,3)
      PIOLA_TENSOR(2,1)=PIOLA_TENSOR(1,2)
      PIOLA_TENSOR(3,1)=PIOLA_TENSOR(1,3)
      PIOLA_TENSOR(3,2)=PIOLA_TENSOR(2,3)
      PIOLA_TENSOR=PIOLA_TENSOR*2.0_DP
    CASE(EQUATIONS_SET_ORTHOTROPIC_MATERIAL_COSTA_SUBTYPE,EQUATIONS_SET_ACTIVECONTRACTION_SUBTYPE) !added by Robert 2010-01-23
      !Form of constitutive model is:
      ! W=a/2 (e^Q - 1)
      ! where Q=[b_ff 2b_fs 2b_fn b_ss 2b_sn b_nn]'* [E_ff E_fs E_fn E_ss E_sn E_nn].^2;
      ! f,s,n denotes the fibre sheet and sheet-normal direction
      a = MATERIALS_INTERPOLATED_POINT%VALUES(1,1)
      B(1,1) = MATERIALS_INTERPOLATED_POINT%VALUES(1+1,1)
      B(1,2) = MATERIALS_INTERPOLATED_POINT%VALUES(1+2,1)
      B(1,3) = MATERIALS_INTERPOLATED_POINT%VALUES(1+3,1)
      B(2,1) = B(1,2);
      B(2,2) = MATERIALS_INTERPOLATED_POINT%VALUES(1+4,1)
      B(2,3) = MATERIALS_INTERPOLATED_POINT%VALUES(1+5,1)
      B(3,1) = B(1,3);
      B(3,2) = B(2,3);
      B(3,3) = MATERIALS_INTERPOLATED_POINT%VALUES(1+6,1)
      Q = 0.0_DP;
      DO i=1,3,1
       DO j=1,3,1
         IF (i==j) THEN
              E(i,j) = 0.5_DP * (AZL(i,j)-1);
         ELSE
              E(i,j) = 0.5_DP * AZL(i,j);
         ENDIF
         Q = Q + B(i,j) * E(i,j) * E(i,j)
       ENDDO
      ENDDO
      Q = exp(Q);
      DO i=1,3,1
       DO j=1,3,1
         PIOLA_TENSOR(i,j)=a*B(i,j)*E(i,j)*Q + p*AZU(i,j);
       ENDDO
      ENDDO

      IF(EQUATIONS_SET_SUBTYPE == EQUATIONS_SET_ACTIVECONTRACTION_SUBTYPE) THEN
        CALL FiniteElasticity_PiolaAddActiveContraction(EQUATIONS_SET%equations%interpolation%independentField, &
          & EQUATIONS_SET%equations%interpolation%materialsField,EQUATIONS_SET%currentTime,EQUATIONS_SET%deltaTime, &
          & PIOLA_TENSOR(1,1),E(1,1),ELEMENT_NUMBER,GAUSS_POINT_NUMBER,err,error,*999)
      ENDIF
    CASE (EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
      & EQUATIONS_SET_COMPRESSIBLE_ACTIVECONTRACTION_SUBTYPE, &
      & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
      !Form of constitutive model is:
      ! W=c1*(I1-3)+c2*(I2-3)+c3*(J-1)^2   (this is actually nearly incompressible)
      C(1)=MATERIALS_INTERPOLATED_POINT%VALUES(1,1)
      C(2)=MATERIALS_INTERPOLATED_POINT%VALUES(2,1)

      PIOLA_TENSOR(1,1)=C(1)+C(2)*(AZL(2,2)+AZL(3,3))
      PIOLA_TENSOR(1,2)=C(2)*(-AZL(2,1))
      PIOLA_TENSOR(1,3)=C(2)*(-AZL(3,1))
      PIOLA_TENSOR(2,1)=PIOLA_TENSOR(1,2)
      PIOLA_TENSOR(2,2)=C(1)+C(2)*(AZL(3,3)+AZL(1,1))
      PIOLA_TENSOR(2,3)=C(2)*(-AZL(3,2))
      PIOLA_TENSOR(3,1)=PIOLA_TENSOR(1,3)
      PIOLA_TENSOR(3,2)=PIOLA_TENSOR(2,3)
      PIOLA_TENSOR(3,3)=C(1)+C(2)*(AZL(1,1)+AZL(2,2))
      PIOLA_TENSOR=PIOLA_TENSOR*2.0_DP

      IF(DIAGNOSTICS1) THEN
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  C(1) = ",C(1),err,error,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  C(2) = ",C(2),err,error,*999)
        CALL WRITE_STRING_MATRIX(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
          & 3,3,AZL,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("    AZL','(",I1,",:)',' :",3(X,E13.6))', &
          & '(17X,3(X,E13.6))',err,error,*999)
      ENDIF

      IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_COMPRESSIBLE_ACTIVECONTRACTION_SUBTYPE) THEN

        CALL Field_ParameterSetGetLocalGaussPoint(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
          &  FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,GAUSS_POINT_NUMBER,ELEMENT_NUMBER,1,ACTIVE_STRESS_11, &
          & err,error,*999) ! get the independent field stress value

        CALL Field_ParameterSetGetLocalGaussPoint(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
          &  FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,GAUSS_POINT_NUMBER,ELEMENT_NUMBER,2,ACTIVE_STRESS_22, &
          & err,error,*999) ! get the independent field stress value

        CALL Field_ParameterSetGetLocalGaussPoint(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
          &  FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,GAUSS_POINT_NUMBER,ELEMENT_NUMBER,3,ACTIVE_STRESS_33, &
          & err,error,*999) ! get the independent field stress value

        PIOLA_TENSOR(1,1)=PIOLA_TENSOR(1,1)+ACTIVE_STRESS_11
        PIOLA_TENSOR(2,2)=PIOLA_TENSOR(2,2)+ACTIVE_STRESS_22
        PIOLA_TENSOR(3,3)=PIOLA_TENSOR(3,3)+ACTIVE_STRESS_33
      ENDIF
      IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE .OR. &
        & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_COMPRESSIBLE_ACTIVECONTRACTION_SUBTYPE) THEN
        C(3)=MATERIALS_INTERPOLATED_POINT%VALUES(3,1)
        PIOLA_TENSOR=PIOLA_TENSOR+2.0_DP*C(3)*(I3-SQRT(I3))*AZU
      ELSEIF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE.OR. &
        & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE .OR. &
        & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE) THEN
        SELECT CASE (EQUATIONS_SET_SUBTYPE)
        CASE (EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE) !Nearly incompressible
          C(3)=MATERIALS_INTERPOLATED_POINT%VALUES(3,1)
          !Starting point for this models is above compressible form of 2nd PK tensor
          !Adjust for the modified Ciarlet-Geymonat expression: Eq.(22) of the INRIA paper
          ! Question is: What deviation is to be penalized : (J-1) or (J-1-m/rho) ??? Probably the latter !
          ! However, m/rho is a given 'constant' and, upon differentiation, drops out.
          ! But it is important to retain I3 = J^2, since J ~ 1 + m/rho /= 1
          PIOLA_TENSOR=PIOLA_TENSOR+C(3)*(SQRT(I3)-1.0_DP)*AZU
          DARCY_MASS_INCREASE_ENTRY = 5 !fifth entry
        CASE (EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
           &  EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE) !Incompressible
          !Constitutive model: W=c1*(I1-3)+c2*(I2-3)+p*(I3-1)
          ! The term 'p*(I3-1)' gives rise to: '2p I3 AZU'
          ! Retain I3 = J^2, since J ~ 1 + m/rho /= 1
!         CASE (EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_MR_SUBTYPE)
          !Constitutive model: W=C1*(J1-3)+C2*(J2-3)+C3*(J-1)^2+lambda.(J-1-m/rho)
          !J1 and J2 are the modified invariants, adjusted for volume change (J1=I1*J^(-2/3), J2=I2*J^(-4/3))
          !Strictly speaking this law isn't for an incompressible material, but the fourth equation in the elasticity
          !is used to satisfy a subtly different constraint, which is to require the solid portion of the poroelastic
          !material retains its volume. (This law is applied on the whole pororous body).

          PIOLA_TENSOR=0.0_DP
          TEMP=0.0_DP

          C(1)=MATERIALS_INTERPOLATED_POINT%VALUES(1,1)
          C(2)=MATERIALS_INTERPOLATED_POINT%VALUES(2,1)
          C(3)=MATERIALS_INTERPOLATED_POINT%VALUES(3,1)

          !J1 term: del(J1)/del(C)=J^(-2/3)*I-2/3*I_1*J^(-2/3)*C^-1
          TEMPTERM=Jznu**(-2.0_DP/3.0_DP)
          TEMP(1,1)=TEMPTERM
          TEMP(2,2)=TEMPTERM
          TEMP(3,3)=TEMPTERM
          I1=AZL(1,1)+AZL(2,2)+AZL(3,3)
          PIOLA_TENSOR=C(1)* (TEMP-1.0_DP/3.0_DP*I1*TEMPTERM*AZU)

          !J2 term: del(J2)/del(C)=J^(-4/3)*del(I2)/del(C) -4/3*I_2*J^(-4/3)*C^-1
          TEMP=MATMUL(AZL,AZL)  ! C^2
          I2=0.5_DP*(I1**2.0_DP-(TEMP(1,1)+TEMP(2,2)+TEMP(3,3)))
          TEMPTERM=Jznu**(-4.0_DP/3.0_DP)
          !TEMP is now del(I2)/del(C)
          TEMP(1,1)=AZL(2,2)+AZL(3,3)
!           TEMP(1,2)=-2.0_DP*AZL(1,2)
          TEMP(1,2)=-1.0_DP*AZL(1,2)
!           TEMP(1,3)=-2.0_DP*AZL(1,3)
          TEMP(1,3)=-1.0_DP*AZL(1,3)
          TEMP(2,1)=TEMP(1,2)
          TEMP(2,2)=AZL(1,1)+AZL(3,3)
!           TEMP(2,3)=-2.0_DP*AZL(2,3)
          TEMP(2,3)=-1.0_DP*AZL(2,3)
          TEMP(3,1)=TEMP(1,3)
          TEMP(3,2)=TEMP(2,3)
          TEMP(3,3)=AZL(1,1)+AZL(2,2)
          PIOLA_TENSOR=PIOLA_TENSOR+C(2)* (TEMPTERM*TEMP-2.0_DP/3.0_DP*I2*TEMPTERM*AZU)

          !J (det(F)) term: (2.C3.(J-1)+lambda)*J.C^-1
          PIOLA_TENSOR=PIOLA_TENSOR+(2.0_DP*C(3)*(Jznu-1.0_DP)+P)*Jznu*AZU

          !Don't forget, it's wrt C so there is a factor of 2 - but not for the pressure !!??
          PIOLA_TENSOR=2.0_DP*PIOLA_TENSOR


          DARCY_MASS_INCREASE_ENTRY = 4 !fourth entry

        END SELECT

!         DARCY_MASS_INCREASE = DARCY_DEPENDENT_INTERPOLATED_POINT%VALUES(DARCY_MASS_INCREASE_ENTRY,NO_PART_DERIV)
! 
!         CALL EVALUATE_CHAPELLE_PIOLA_TENSOR_ADDITION(AZL,AZU,DARCY_MASS_INCREASE,PIOLA_TENSOR_ADDITION,err,error,*999)
! 
!         IF(DIAGNOSTICS1) THEN
!           CALL WRITE_STRING_MATRIX(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
!             & 3,3,PIOLA_TENSOR,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("    PIOLA_TENSOR','(",I1,",:)',' :",3(X,E13.6))', &
!             & '(17X,3(X,E13.6))',err,error,*999)
!           CALL WRITE_STRING_MATRIX(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
!             & 3,3,PIOLA_TENSOR_ADDITION, &
!             & WRITE_STRING_MATRIX_NAME_AND_INDICES,'("    PIOLA_TENSOR_ADDITION','(",I1,",:)',' :",3(X,E13.6))', &
!             & '(17X,3(X,E13.6))',err,error,*999)
!         ENDIF
!
!         PIOLA_TENSOR = PIOLA_TENSOR + PIOLA_TENSOR_ADDITION
      ENDIF

    CASE (EQUATIONS_SET_ORTHOTROPIC_MATERIAL_HOLZAPFEL_OGDEN_SUBTYPE, &
        & EQUATIONS_SET_HOLZAPFEL_OGDEN_ACTIVECONTRACTION_SUBTYPE)
      !Form of the constitutive model is:
      ! W = a/(2*b)*exp[b*(I1-3)] + sum_(i=f,s)[H(I4i-1)*a_i/(2*b_i)*(exp[b_i*(I4i-1)^2]-1)] + a_fs/(2*b_fs)*(exp[b_fs*I8fs^2]-1)
      !where H is the Heaviside step function. Fibres only contribute stiffness if in tension.
      !Also assumed I3 = det(AZL) = J^2 = 1.0  -  incompressible material
      !Assume directions: fibre f_0=[1 0 0], sheet s_0=[0 1 0], (sheet) normal n_0=[0 0 1]
      !Based on: Holzapfel, G. A., & Ogden, R. W. (2009). Constitutive modelling of passive myocardium: A structurally based
      !  framework for material characterization. Philosophical Transactions of the Royal Society A: Mathematical, Physical and
      !  Engineering Sciences, 367(1902), 3445-3475. doi:10.1098/rsta.2009.0091
      C(1)=MATERIALS_INTERPOLATED_POINT%VALUES(1,1) !a
      C(2)=MATERIALS_INTERPOLATED_POINT%VALUES(2,1) !b
      C(3)=MATERIALS_INTERPOLATED_POINT%VALUES(3,1) !a_f
      C(4)=MATERIALS_INTERPOLATED_POINT%VALUES(4,1) !a_s
      C(5)=MATERIALS_INTERPOLATED_POINT%VALUES(5,1) !b_f
      C(6)=MATERIALS_INTERPOLATED_POINT%VALUES(6,1) !b_s
      C(7)=MATERIALS_INTERPOLATED_POINT%VALUES(7,1) !a_fs
      C(8)=MATERIALS_INTERPOLATED_POINT%VALUES(8,1) !b_fs
      I1=AZL(1,1)+AZL(2,2)+AZL(3,3)
      TEMPTERM=C(1)*EXP(C(2)*(I1-3.0_DP))
      PIOLA_TENSOR(1,1)=-P*AZU(1,1)+TEMPTERM
      IF(AZL(1,1)>1.0_DP) THEN
        PIOLA_TENSOR(1,1)=PIOLA_TENSOR(1,1)+2.0_DP*C(3)*(AZL(1,1)-1.0_DP)*EXP(C(5)*(AZL(1,1)-1.0_DP)**2.0_DP)
      END IF
      PIOLA_TENSOR(1,2)=-P*AZU(1,2)+C(7)*AZL(1,2)*EXP(C(8)*AZL(1,2)**2.0_DP)
      PIOLA_TENSOR(1,3)=-P*AZU(1,3)
      PIOLA_TENSOR(2,1)=PIOLA_TENSOR(1,2)
      PIOLA_TENSOR(2,2)=-P*AZU(2,2)+TEMPTERM
      IF(AZL(2,2)>1.0_DP) THEN
        PIOLA_TENSOR(2,2)=PIOLA_TENSOR(2,2)+2.0_DP*C(4)*(AZL(2,2)-1.0_DP)*EXP(C(6)*(AZL(2,2)-1.0_DP)**2.0_DP)
      END IF
      PIOLA_TENSOR(2,3)=-P*AZU(2,3)
      PIOLA_TENSOR(3,1)=PIOLA_TENSOR(1,3)
      PIOLA_TENSOR(3,2)=PIOLA_TENSOR(2,3)
      PIOLA_TENSOR(3,3)=-P*AZU(3,3)+TEMPTERM

      IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_HOLZAPFEL_OGDEN_ACTIVECONTRACTION_SUBTYPE) THEN
      !add active contraction stress value to the trace of the stress tensor - basically adding to hydrostatic pressure.
      !the active stress is stored inside the independent field that has been set up in the user program.
      !for generality we could set up 3 components in independent field for 3 different active stress components
        CALL Field_VariableGet(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VARIABLE,err,error,*999)
        DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
          dof_idx=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP% &
            & GAUSS_POINTS(GAUSS_POINT_NUMBER,ELEMENT_NUMBER)
          CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,dof_idx,VALUE,err,error,*999)
          PIOLA_TENSOR(component_idx,component_idx)=PIOLA_TENSOR(component_idx,component_idx)+VALUE
        ENDDO
      ENDIF

    CASE DEFAULT
      LOCAL_ERROR="The third equations set specification of "//TRIM(NumberToVString(EQUATIONS_SET_SUBTYPE,"*",err,error))// &
        & " is not valid for a finite elasticity type of an elasticity equation set."
      CALL FlagError(LOCAL_ERROR,err,error,*999)
    END SELECT

    CALL MatrixProduct(DZDNU,PIOLA_TENSOR,TEMP,err,error,*999)
    CALL MatrixProduct(TEMP,DZDNUT,CAUCHY_TENSOR,err,error,*999)
    
    CAUCHY_TENSOR=CAUCHY_TENSOR/Jznu
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  ELEMENT_NUMBER = ",ELEMENT_NUMBER,err,error,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  gauss_idx = ",GAUSS_POINT_NUMBER,err,error,*999)
      CALL WRITE_STRING_MATRIX(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
        & 3,3,PIOLA_TENSOR,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("    PIOLA_TENSOR','(",I1,",:)',' :",3(X,E13.6))', &
        & '(17X,3(X,E13.6))',err,error,*999)
      CALL WRITE_STRING_MATRIX(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
        & 3,3,CAUCHY_TENSOR,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("    CAUCHY_TENSOR','(",I1,",:)',' :",3(X,E13.6))', &
        & '(17X,3(X,E13.6))',err,error,*999)
    ENDIF
    NULLIFY(C)

    EXITS("FINITE_ELASTICITY_GAUSS_CAUCHY_TENSOR")
    RETURN
999 ERRORSEXITS("FINITE_ELASTICITY_GAUSS_CAUCHY_TENSOR",err,error)
    RETURN 1
  END SUBROUTINE FINITE_ELASTICITY_GAUSS_CAUCHY_TENSOR

  !
  !================================================================================================================================
  !

  !>Evaluates the growth tensor at a given Gauss point and calculates the elastic part of the deformation gradient tensor
  SUBROUTINE FiniteElasticity_GaussGrowthTensor(equationsSet,numberOfDimensions,deformationGradientTensor,growthValues, &
    & growthTensor,elasticDeformationGradientTensor,Jg,Je,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER, INTENT(IN) :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: numberOfDimensions !<The number of dimensions
    REAL(DP), INTENT(IN) :: deformationGradientTensor(3,3) !<The full deformation gradient tensor
    REAL(DP), INTENT(IN) :: growthValues(3) !<The fibre, sheet and normal growth extensions.
    REAL(DP), INTENT(OUT) :: growthTensor(3,3) !<On output, the growth tensor
    REAL(DP), INTENT(OUT) :: elasticDeformationGradientTensor(3,3) !<On output, the elastic part of the deformation gradient tensor
    REAL(DP), INTENT(OUT) :: Jg !<On output, the Jacobian of the growth tensor
    REAL(DP), INTENT(OUT) :: Je !<On output, the Jacobian of the elastic tensor
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    REAL(DP) :: growthTensorInverse(3,3),J
    
    ENTERS("FiniteElasticity_GaussGrowthTensor",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      CALL IdentityMatrix(growthTensor,err,error,*999)
      IF(equationsSet%specification(3)==EQUATIONS_SET_CONSTITUTIVE_AND_GROWTH_LAW_IN_CELLML_SUBTYPE.OR. &
        equationsSet%specification(3)==EQUATIONS_SET_MR_AND_GROWTH_LAW_IN_CELLML_SUBTYPE) THEN
        growthTensor(1,1)=growthValues(1)
        IF(numberofDimensions>1) THEN
          growthTensor(2,2)=growthValues(2)
          IF(numberOfDimensions>2) THEN
            growthTensor(3,3)=growthValues(3)
          ENDIF
        ENDIF
        !Calculate inverse growth deformation tensor, Fg^-1, Jg 
        CALL Invert(growthTensor,growthTensorInverse,Jg,err,error,*999)
        !Calculate elastic deformation tensor, Fe=F.(Fg)^-1.       
        CALL MatrixProduct(deformationGradientTensor,growthTensorInverse,elasticDeformationGradientTensor,err,error,*999)
      ELSE
        Jg=1.0_DP
        elasticDeformationGradientTensor=deformationGradientTensor
      ENDIF
      CALL Determinant(elasticDeformationGradientTensor,Je,err,error,*999)
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

    IF(diagnostics1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"",err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Growth information:",err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Total deformation gradient tensor:",err,error,*999)
      CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3,3,3,deformationGradientTensor, &
        & WRITE_STRING_MATRIX_NAME_AND_INDICES,'("    F','(",I1,",:)','  :",3(X,E13.6))','(13X,3(X,E13.6))',err,error,*999)
      CALL Determinant(deformationGradientTensor,J,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Determinant F, J = ",J,err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Elastic component of the deformation gradient tensor:",err,error,*999)
      CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3,3,3,elasticDeformationGradientTensor, &
        & WRITE_STRING_MATRIX_NAME_AND_INDICES,'("    Fe','(",I1,",:)',' :",3(X,E13.6))','(13X,3(X,E13.6))',err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Determinant Fe, Je = ",Je,err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Growth component of the deformation gradient tensor:",err,error,*999)
      CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3,3,3,growthTensor, &
        & WRITE_STRING_MATRIX_NAME_AND_INDICES,'("    Fg','(",I1,",:)',' :",3(X,E13.6))','(13X,3(X,E13.6))',err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Determinant Fg, Jg = ",Jg,err,error,*999)
    ENDIF
   
    EXITS("FiniteElasticity_GaussGrowthTensor")
    RETURN
    999 ERRORSEXITS("FiniteElasticity_GaussGrowthTensor",err,error)
    RETURN 1
    
  END SUBROUTINE FiniteElasticity_GaussGrowthTensor

 !
  !================================================================================================================================
  !

  !>Evaluates the strain tensor given the deformation gradient tensor
  SUBROUTINE FiniteElasticity_StrainTensor(deformationGradientTensor,rightCauchyDeformationTensor,fingerDeformationTensor, &
    jacobian,greenStrainTensor,err,error,*)

    !Argument variables
    REAL(DP), INTENT(IN) :: deformationGradientTensor(3,3) !<The elastic part of the  deformation gradient tensor
    REAL(DP), INTENT(OUT) :: rightCauchyDeformationTensor(3,3) !<On output, the right Cauchy deformation tensor, C
    REAL(DP), INTENT(OUT) :: fingerDeformationTensor(3,3) !<On output, the finger deformation tensor, f
    REAL(DP), INTENT(OUT) :: jacobian !<On output, the Jacobian of the deformation
    REAL(DP), INTENT(OUT) :: greenStrainTensor(3,3) !<On output, the Green-Lagrange strain tensor
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: i
    REAL(DP) :: I3
    
    ENTERS("FiniteElasticity_StrainTensor",err,error,*999)

    CALL MatrixTransposeProduct(deformationGradientTensor,deformationGradientTensor,rightCauchyDeformationTensor,err,error,*999)
    CALL Invert(rightCauchyDeformationTensor,fingerDeformationTensor,I3,err,error,*999)
    CALL Determinant(deformationGradientTensor,jacobian,err,error,*999)

    greenStrainTensor=0.5_DP*rightCauchyDeformationTensor
    DO i=1,3
      greenStrainTensor(i,i)=greenStrainTensor(i,i)-0.5_DP
    ENDDO !i
    
    IF(diagnostics1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"",err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Strain information:",err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Right Cauchy-Green deformation tensor:",err,error,*999)
      CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
        & 3,3,rightCauchyDeformationTensor,WRITE_STRING_MATRIX_NAME_AND_INDICES, '("    C','(",I1,",:)', &
        & ' :",3(X,E13.6))','(12X,3(X,E13.6))',err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Finger deformation tensor:",err,error,*999)
      CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
        & 3,3,fingerDeformationTensor,WRITE_STRING_MATRIX_NAME_AND_INDICES, '("    f','(",I1,",:)', &
        & ' :",3(X,E13.6))','(12X,3(X,E13.6))',err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Jacobian = ",jacobian,err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Green-Lagrange strain tensor:",err,error,*999)
      CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
        & 3,3,greenStrainTensor,WRITE_STRING_MATRIX_NAME_AND_INDICES, '("    E','(",I1,",:)', &
        & ' :",3(X,E13.6))','(12X,3(X,E13.6))',err,error,*999)
    ENDIF
   
    EXITS("FiniteElasticity_StrainTensor")
    RETURN
    999 ERRORSEXITS("FiniteElasticity_StrainTensor",err,error)
    RETURN 1
    
  END SUBROUTINE FiniteElasticity_StrainTensor

  !
  !================================================================================================================================
  !

  !>Evaluates the Cauchy stress tensor at a given Gauss point
  SUBROUTINE FINITE_ELASTICITY_GAUSS_STRESS_TENSOR(EQUATIONS_SET,DEPENDENT_INTERPOLATED_POINT, &
      & MATERIALS_INTERPOLATED_POINT,GEOMETRIC_INTERPOLATED_POINT,STRESS_TENSOR,DZDNU,Jznu, &
      & ELEMENT_NUMBER,GAUSS_POINT_NUMBER,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER, INTENT(IN) :: EQUATIONS_SET !<A pointer to the equations set
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: DEPENDENT_INTERPOLATED_POINT,MATERIALS_INTERPOLATED_POINT
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: GEOMETRIC_INTERPOLATED_POINT
    REAL(DP), INTENT(OUT) :: STRESS_TENSOR(:)
    REAL(DP), INTENT(IN) :: DZDNU(3,3) !Deformation gradient tensor at the gauss point
    REAL(DP), INTENT(IN) :: Jznu !Determinant of deformation gradient tensor (AZL)
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER,GAUSS_POINT_NUMBER !<Element/Gauss point number
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: PRESSURE_COMPONENT,component_idx,dof_idx
    REAL(DP) :: P
    REAL(DP) :: I1 !Invariants, if needed
    REAL(DP) :: TEMPTERM1,TEMPTERM2,VALUE !Temporary variables
    REAL(DP) :: ONETHIRD_TRACE
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    REAL(DP) :: MOD_DZDNU(3,3),MOD_DZDNUT(3,3),AZL(3,3)
    REAL(DP) :: B(6),E(6),DQ_DE(6)
    REAL(DP), POINTER :: C(:) !Parameters for constitutive laws

    ENTERS("FINITE_ELASTICITY_GAUSS_STRESS_TENSOR",err,error,*999)

    NULLIFY(FIELD_VARIABLE,C)

    !AZL = F'*F (deformed covariant or right cauchy deformation tensor, C)
    !AZU - deformed contravariant tensor; I3 = det(C)

    MOD_DZDNU=DZDNU*Jznu**(-1.0_DP/3.0_DP)
    CALL MatrixTranspose(MOD_DZDNU,MOD_DZDNUT,err,error,*999)
    CALL MatrixProduct(MOD_DZDNUT,MOD_DZDNU,AZL,err,error,*999)
    C=>MATERIALS_INTERPOLATED_POINT%VALUES(:,NO_PART_DERIV)

    SELECT CASE(EQUATIONS_SET%specification(3))
    CASE(EQUATIONS_SET_MOONEY_RIVLIN_ACTIVECONTRACTION_SUBTYPE, &
      & EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE, &
      & EQUATIONS_SET_MR_AND_GROWTH_LAW_IN_CELLML_SUBTYPE)
      PRESSURE_COMPONENT=DEPENDENT_INTERPOLATED_POINT%INTERPOLATION_PARAMETERS%FIELD_VARIABLE%NUMBER_OF_COMPONENTS
      P=DEPENDENT_INTERPOLATED_POINT%VALUES(PRESSURE_COMPONENT,NO_PART_DERIV)
      !Form of constitutive model is:
      !W=c1*(I1-3)+c2*(I2-3)+p/2*(I3-1)

      !Calculate isochoric fictitious 2nd Piola tensor (in Voigt form)
      I1=AZL(1,1)+AZL(2,2)+AZL(3,3)
      TEMPTERM1=-2.0_DP*C(2)
      TEMPTERM2=2.0_DP*(C(1)+I1*C(2))
      STRESS_TENSOR(1)=TEMPTERM1*AZL(1,1)+TEMPTERM2
      STRESS_TENSOR(2)=TEMPTERM1*AZL(2,2)+TEMPTERM2
      STRESS_TENSOR(3)=TEMPTERM1*AZL(3,3)+TEMPTERM2
      STRESS_TENSOR(4)=TEMPTERM1*AZL(2,1)
      STRESS_TENSOR(5)=TEMPTERM1*AZL(3,1)
      STRESS_TENSOR(6)=TEMPTERM1*AZL(3,2)

      IF(EQUATIONS_SET%specification(3)==EQUATIONS_SET_MOONEY_RIVLIN_ACTIVECONTRACTION_SUBTYPE) THEN
        !add active contraction stress values
        !Be aware for modified DZDNU, should active contraction be added here? Normally should be okay as modified DZDNU and DZDNU
        !converge during the Newton iteration.
        CALL Field_VariableGet(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VARIABLE,err,error,*999)
        DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
          dof_idx=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP% &
            & GAUSS_POINTS(GAUSS_POINT_NUMBER,ELEMENT_NUMBER)
          CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,dof_idx,VALUE,err,error,*999)
          STRESS_TENSOR(component_idx)=STRESS_TENSOR(component_idx)+VALUE
        ENDDO
      ENDIF

      !Do push-forward of 2nd Piola tensor. 
      CALL FINITE_ELASTICITY_PUSH_STRESS_TENSOR(STRESS_TENSOR,MOD_DZDNU,Jznu,err,error,*999)
      !Calculate isochoric Cauchy tensor (the deviatoric part) and add the volumetric part (the hydrostatic pressure).
      ONETHIRD_TRACE=SUM(STRESS_TENSOR(1:3))/3.0_DP
      STRESS_TENSOR(1:3)=STRESS_TENSOR(1:3)-ONETHIRD_TRACE+P

    CASE(EQUATIONS_SET_TRANSVERSE_ISOTROPIC_GUCCIONE_SUBTYPE,EQUATIONS_SET_GUCCIONE_ACTIVECONTRACTION_SUBTYPE, &
      & EQUATIONS_SET_REFERENCE_STATE_TRANSVERSE_GUCCIONE_SUBTYPE)
      IF(EQUATIONS_SET%specification(3)==EQUATIONS_SET_REFERENCE_STATE_TRANSVERSE_GUCCIONE_SUBTYPE) THEN
        PRESSURE_COMPONENT=GEOMETRIC_INTERPOLATED_POINT%INTERPOLATION_PARAMETERS%FIELD_VARIABLE%NUMBER_OF_COMPONENTS
        P=GEOMETRIC_INTERPOLATED_POINT%VALUES(PRESSURE_COMPONENT,NO_PART_DERIV)
      ELSE
        PRESSURE_COMPONENT=DEPENDENT_INTERPOLATED_POINT%INTERPOLATION_PARAMETERS%FIELD_VARIABLE%NUMBER_OF_COMPONENTS
        P=DEPENDENT_INTERPOLATED_POINT%VALUES(PRESSURE_COMPONENT,NO_PART_DERIV)
      ENDIF
      B=[2.0_DP*C(2),2.0_DP*C(3),2.0_DP*C(3),C(4),C(4),C(3)] ![2*b_f,2*b_t,2*b_t,b_ft,b_ft,b_t]
      E=[0.5_DP*(AZL(1,1)-1.0_DP),0.5_DP*(AZL(2,2)-1.0_DP),0.5_DP*(AZL(3,3)-1.0_DP),AZL(2,1),AZL(3,1),AZL(3,2)] !(Modified) strain tensor in Voigt form.
      DQ_DE=B*E
      TEMPTERM1=0.5_DP*C(1)*EXP(0.5_DP*DOT_PRODUCT(E,DQ_DE))
      ! Calculate isochoric fictitious 2nd Piola tensor (in Voigt form)
      STRESS_TENSOR=TEMPTERM1*DQ_DE
      IF(EQUATIONS_SET%specification(3)==EQUATIONS_SET_GUCCIONE_ACTIVECONTRACTION_SUBTYPE) THEN
        !add active contraction stress values
        !Be aware for modified DZDNU, should active contraction be added here? Normally should be okay as modified DZDNU and DZDNU
        !converge during the Newton iteration.
        CALL Field_VariableGet(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VARIABLE,err,error,*999)
        DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
          dof_idx=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP% &
            & GAUSS_POINTS(GAUSS_POINT_NUMBER,ELEMENT_NUMBER)
          CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,dof_idx,VALUE,err,error,*999)
          STRESS_TENSOR(component_idx)=STRESS_TENSOR(component_idx)+VALUE
        ENDDO
      ENDIF
      ! Do push-forward of 2nd Piola tensor. 
      CALL FINITE_ELASTICITY_PUSH_STRESS_TENSOR(STRESS_TENSOR,MOD_DZDNU,Jznu,err,error,*999)
      !Calculate isochoric Cauchy tensor (the deviatoric part) and add the volumetric part (the hydrostatic pressure).
      ONETHIRD_TRACE=SUM(STRESS_TENSOR(1:3))/3.0_DP
      STRESS_TENSOR(1:3)=STRESS_TENSOR(1:3)-ONETHIRD_TRACE+P
    CASE DEFAULT
      LOCAL_ERROR="The third equations set specification of "// &
        & TRIM(NumberToVString(EQUATIONS_SET%specification(3),"*",err,error))// &
        & " is not valid for a finite elasticity type of an elasticity equation set."
     CALL FlagError(LOCAL_ERROR,err,error,*999)
    END SELECT

    EXITS("FINITE_ELASTICITY_GAUSS_STRESS_TENSOR")
    RETURN
    999 ERRORSEXITS("FINITE_ELASTICITY_GAUSS_STRESS_TENSOR",err,error)
    RETURN 1
  END SUBROUTINE FINITE_ELASTICITY_GAUSS_STRESS_TENSOR

  !
  !================================================================================================================================
  !

  ! calculates the current active contraction component using the independent field
  ! Uses a hardcoded tension transient based on GPB+NHS with length-dependence for now
  SUBROUTINE FiniteElasticity_PiolaAddActiveContraction(INDEPENDENT_FIELD,MATERIALS_FIELD,currentTime,dt,PIOLA_FF,E_FF,&
             & ELEMENT_NUMBER,GAUSS_POINT_NUMBER,err,error,*)
    !Argument variables
    TYPE(FIELD_TYPE), POINTER, INTENT(IN) :: INDEPENDENT_FIELD
    TYPE(FIELD_TYPE), POINTER, INTENT(IN) :: MATERIALS_FIELD
    REAL(DP), INTENT(IN) :: currentTime !<The time to evaluate at
    REAL(DP), INTENT(IN) :: dt !<The delta time to evaluate at
    REAL(DP), INTENT(INOUT) :: PIOLA_FF  !<The (1,1)=(fiber,fiber) component of the stress tensor
    REAL(DP), INTENT(IN)    :: E_FF !<E(1,1)
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER,GAUSS_POINT_NUMBER !<Element/Gauss point number
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    INTEGER(INTG)  :: I
    REAL(DP) :: S, LAMBDA, ISO_TA, TA, ACTIVTIME, TIME
    REAL(DP), DIMENSION(1:4) :: QL

    REAL(DP), PARAMETER :: PERIOD = 1000 ! 1 Hz
    REAL(DP), PARAMETER, DIMENSION(28) :: TIMES    =    [ 0, 20, 30, 40, 60, 80, 100, 120, 150, 160, 170, 175, 180, 190, 200,&
    & 225, 250, 300, 333, 366, 400, 450, 500, 600, 700, 800, 900,1000 ] ! simple tension curve based on GPB/NHS: times

    REAL(DP), PARAMETER, DIMENSION(28) :: TENSIONFRAC = [ 0.0194, 0.0193, 0.0200, 0.0254, 0.0778, 0.1713, 0.2794, 0.3708,&
    & 0.4472, 0.4578, 0.4624, 0.4627, 0.4618, 0.4567, 0.4478, 0.4121, 0.3614, 0.2326, 0.1471, 0.0920, 0.0681, 0.0526, 0.0438,&
    & 0.0332, 0.0271, 0.0234, 0.0210, 0.0194 ] ! simple isometric tension curve based on GPB/NHS: tension/tref 
    REAL(DP), PARAMETER :: T_REF = 100          ! reference tension
  
    ENTERS("FiniteElasticity_PiolaAddActiveContraction",err,error,*999)

    ! Get Q's
    DO I=1,4
      CALL Field_ParameterSetGetLocalGaussPoint(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,&
        &  FIELD_VALUES_SET_TYPE,GAUSS_POINT_NUMBER,ELEMENT_NUMBER,I,QL(I),err,error,*999)  ! Q(1) Q(2) Q(3) Lambda for prev in 1/2/3/4
    END DO

    ! get activation time from material field
    CALL Field_ParameterSetGetLocalGaussPoint(MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE,&
      &  FIELD_VALUES_SET_TYPE,GAUSS_POINT_NUMBER,ELEMENT_NUMBER,1,ACTIVTIME,err,error,*999)

    LAMBDA = SQRT(2*E_FF + 1)
    TIME =  MAX( MOD(currentTime, PERIOD) - ACTIVTIME, 0.0) ! start activation at this time
   
    I = 1
    DO WHILE (TIMES(I) <= TIME) ! find first I such that times(I) >= time
      I = I+1
    END DO
    S    = (TIME - TIMES(I-1)) /  (TIMES(I) - TIMES(I-1))                     !| linear interpolation of ta/tref
    ISO_TA   = T_REF * (TENSIONFRAC(I-1) * (1-S) + TENSIONFRAC(I) * S)        !/ + multiply by tref
  
    CALL FINITE_ELASTICITY_FMM(TIME,DT,QL(4),LAMBDA,QL,ISO_TA,TA,err,error,*999)

    QL(4) = LAMBDA  ! bounds applied in FMM, Qi integrated
    DO I=1,4
      CALL FIELD_PARAMETER_SET_UPDATE_GAUSS_POINT(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,&
        &  FIELD_VALUES_SET_TYPE,GAUSS_POINT_NUMBER,ELEMENT_NUMBER, 4+I, QL(I),err,error,*999) ! store Q(1) Q(2) Q(3) Lambda for next in 5/6/7/8
    END DO

    PIOLA_FF = PIOLA_FF + TA

    EXITS("FiniteElasticity_PiolaAddActiveContraction")
    RETURN
999 ERRORSEXITS("FiniteElasticity_PiolaAddActiveContraction",err,error)
    RETURN 1

  END SUBROUTINE FiniteElasticity_PiolaAddActiveContraction

  !
  !================================================================================================================================
  !

  ! Implements length and velocity dependence. can be used in both weak and strong coupling
  SUBROUTINE FINITE_ELASTICITY_FMM(TIME,DT,PREV_LAMBDA,CURR_LAMBDA,Q123,ISO_TA,TA,err,error,*)
    ! PARAMETERS FROM Niederer Hunter & Smith 2006
    REAL(DP), PARAMETER, DIMENSION(1:3) :: A     = [-29.0,138.0,129.0]  ! 'A'
    REAL(DP), PARAMETER, DIMENSION(1:3) :: ALPHA = [0.03,0.13,0.625]
    REAL(DP), PARAMETER :: la   = 0.35, BETA_0 = 4.9  ! 'a'

    REAL(DP), INTENT(INOUT), DIMENSION(:) :: Q123
    REAL(DP), INTENT(INOUT) :: CURR_LAMBDA
    REAL(DP), INTENT(IN) :: PREV_LAMBDA, DT, TIME, ISO_TA
    REAL(DP), INTENT(OUT) :: TA

    INTEGER(INTG) :: ERR
    TYPE(VARYING_STRING) :: ERROR

    REAL(DP) :: QFAC, DLAMBDA_DT, Q, OVERLAP
    INTEGER(INTG) :: I

    ENTERS("FINITE_ELASTICITY_FMM",err,error,*999)

    CURR_LAMBDA = MIN(1.15, MAX(0.8, CURR_LAMBDA))  ! inout -> save this

    IF( TIME - 1e-10 <= 0.0) THEN  ! preload / first step -> update method off
      QFAC = 1.0
    ELSE
      DLAMBDA_DT = (CURR_LAMBDA - PREV_LAMBDA) / DT
      DO I=1,3
        Q123(I) = Q123(I) + DT * (A(I) * DLAMBDA_DT - ALPHA(I) * Q123(I))
      END DO
      Q = Q123(1)+Q123(2)+Q123(3)
      IF(Q < 0.0) THEN
        QFAC = (la*Q + 1.0) / (1.0 - Q)
      ELSE
        QFAC = (1.0 + (la+2.0)*Q)/(1.0+Q);
      END IF
    END IF

    OVERLAP= 1.0 + BETA_0 * (CURR_LAMBDA-1.0)
    TA = OVERLAP * QFAC * ISO_TA  ! length dep * vel dep * isometric tension
    
    EXITS("FINITE_ELASTICITY_FMM")
    RETURN
999 ERRORSEXITS("FINITE_ELASTICITY_FMM",err,error)
    RETURN 1

  END SUBROUTINE FINITE_ELASTICITY_FMM


  !
  !================================================================================================================================
  !

  !>Evaluates df/dz (derivative of interpolation function wrt deformed coord) matrix at a given Gauss point
  SUBROUTINE FINITE_ELASTICITY_GAUSS_DFDZ(INTERPOLATED_POINT,ELEMENT_NUMBER,GAUSS_POINT_NUMBER,NUMBER_OF_DIMENSIONS, &
    & NUMBER_OF_XI,DFDZ,err,error,*)

    !Argument variables
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: INTERPOLATED_POINT !<Interpolated point for the dependent field
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number
    INTEGER(INTG), INTENT(IN) :: GAUSS_POINT_NUMBER !<The gauss point number
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_DIMENSIONS !<The number of dimensions
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_XI !<The number of xi directions for the interpolation
    REAL(DP), INTENT(OUT) :: DFDZ(:,:,:) !<On return, a matrix containing the derivatives of the basis functions wrt the deformed coordinates
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(BASIS_TYPE), POINTER :: COMPONENT_BASIS
    TYPE(FIELD_TYPE), POINTER :: FIELD
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME
    INTEGER(INTG) :: derivative_idx,component_idx1,component_idx2,xi_idx,parameter_idx
    REAL(DP) :: DXIDZ(NUMBER_OF_DIMENSIONS,NUMBER_OF_DIMENSIONS),DZDXI(NUMBER_OF_DIMENSIONS,NUMBER_OF_DIMENSIONS)
    REAL(DP) :: Jzxi,DFDXI(NUMBER_OF_DIMENSIONS,64,NUMBER_OF_XI)!temporary until a proper alternative is found
    
    ENTERS("FINITE_ELASTICITY_GAUSS_DFDZ",err,error,*999)

    !Initialise DFDXI array
    DFDXI=0.0_DP  ! DFDXI(component_idx,parameter_idx,xi_idx)
    DFDZ=0.0_DP
    DO component_idx2=1,NUMBER_OF_DIMENSIONS !Always 3 spatial coordinates (3D)
      DO xi_idx=1,NUMBER_OF_XI !Thus always 3 element coordinates
        derivative_idx=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xi_idx)  !2,4,7
        DZDXI(component_idx2,xi_idx)=INTERPOLATED_POINT%VALUES(component_idx2,derivative_idx)  !dz/dxi
      ENDDO
    ENDDO

    ! Populate a 3 x 3 square dzdXi if this is a membrane problem in 3D space
    IF (NUMBER_OF_DIMENSIONS == 3 .AND. NUMBER_OF_XI == 2) THEN
        CALL CrossProduct(DZDXI(:,1),DZDXI(:,2),DZDXI(:,3),err,error,*999)
        CALL Normalise(DZDXI(:,3),DZDXI(:,3),err,error,*999)
    ENDIF

    CALL INVERT(DZDXI,DXIDZ,Jzxi,err,error,*999) !dxi/dz

    FIELD=>INTERPOLATED_POINT%INTERPOLATION_PARAMETERS%FIELD
    DO component_idx1=1,NUMBER_OF_DIMENSIONS
      COMPONENT_BASIS=>FIELD%VARIABLES(1)%COMPONENTS(component_idx1)%DOMAIN%TOPOLOGY%ELEMENTS% &
        & ELEMENTS(ELEMENT_NUMBER)%BASIS
      QUADRATURE_SCHEME=>COMPONENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
      DO parameter_idx=1,COMPONENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
        DO xi_idx=1,NUMBER_OF_XI
          derivative_idx=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xi_idx)
          DFDXI(component_idx1,parameter_idx,xi_idx)=QUADRATURE_SCHEME%GAUSS_BASIS_FNS(parameter_idx,derivative_idx, &
            & GAUSS_POINT_NUMBER)
        ENDDO
      ENDDO
    ENDDO

    DO component_idx1=1,NUMBER_OF_DIMENSIONS
      COMPONENT_BASIS=>FIELD%VARIABLES(1)%COMPONENTS(component_idx1)%DOMAIN%TOPOLOGY%ELEMENTS% &
        & ELEMENTS(ELEMENT_NUMBER)%BASIS
      DO component_idx2=1,NUMBER_OF_DIMENSIONS
        DO parameter_idx=1,COMPONENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
          DO xi_idx=1,NUMBER_OF_XI
            DFDZ(parameter_idx,component_idx2,component_idx1)=DFDZ(parameter_idx,component_idx2,component_idx1) + &
              & DFDXI(component_idx1,parameter_idx,xi_idx) * DXIDZ(xi_idx,component_idx2)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    EXITS("FINITE_ELASTICITY_GAUSS_DFDZ")
    RETURN
999 ERRORSEXITS("FINITE_ELASTICITY_GAUSS_DFDZ",err,error)
    RETURN 1
  END SUBROUTINE FINITE_ELASTICITY_GAUSS_DFDZ

  !
  !================================================================================================================================
  !

  !>Sets up the finite elasticity equation type of an elasticity equations set class.
  SUBROUTINE FINITE_ELASTICITY_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup a Laplace equation on.
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR           !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR  !<The error string
    !Local Variables
    INTEGER(INTG) :: GEOMETRIC_MESH_COMPONENT,GEOMETRIC_SCALING_TYPE,NUMBER_OF_COMPONENTS, &
      & numberOfDimensions,NUMBER_OF_DARCY_COMPONENTS,GEOMETRIC_COMPONENT_NUMBER,NUMBER_OF_COMPONENTS_2,component_idx, &
      & componentIdx,derivedIdx,varIdx,variableType,NUMBER_OF_FLUID_COMPONENTS,numberOfTensorComponents
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: coordinateSystem
    TYPE(DECOMPOSITION_TYPE), POINTER :: GEOMETRIC_DECOMPOSITION
    TYPE(FIELD_TYPE), POINTER :: ANALYTIC_FIELD,DEPENDENT_FIELD,GEOMETRIC_FIELD
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EQUATIONS_SET_EQUATIONS_SET_FIELD_TYPE), POINTER :: EQUATIONS_EQUATIONS_SET_FIELD
    TYPE(FIELD_TYPE), POINTER :: EQUATIONS_SET_FIELD_FIELD
    TYPE(EQUATIONS_SET_MATERIALS_TYPE), POINTER :: EQUATIONS_MATERIALS
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR,localError
    LOGICAL :: IS_HYDROSTATIC_PRESSURE_DEPENDENT_FIELD
    INTEGER(INTG) :: num_var,Ncompartments,DEPENDENT_FIELD_NUMBER_OF_VARIABLES,dimensionIdx
    INTEGER(INTG) :: EQUATIONS_SET_FIELD_NUMBER_OF_VARIABLES,EQUATIONS_SET_FIELD_NUMBER_OF_COMPONENTS    
    INTEGER(INTG), POINTER :: EQUATIONS_SET_FIELD_DATA(:)
    INTEGER(INTG), ALLOCATABLE :: VARIABLE_TYPES(:)
    INTEGER(INTG) :: EQUATIONS_SET_SUBTYPE

    ENTERS("FINITE_ELASTICITY_EQUATIONS_SET_SETUP",err,error,*999)

    NULLIFY(GEOMETRIC_DECOMPOSITION)
    NULLIFY(EQUATIONS)
    NULLIFY(vectorMapping)
    NULLIFY(vectorMatrices)
    NULLIFY(EQUATIONS_MATERIALS)
    NULLIFY(EQUATIONS_EQUATIONS_SET_FIELD)
    NULLIFY(EQUATIONS_SET_FIELD_FIELD)
    NULLIFY(EQUATIONS_SET_FIELD_DATA)

    IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<3) &
      & CALL FlagError("Equations set specification must have at least three entries for a finite elasticity type equations set.", &
        & err,error,*999)
    EQUATIONS_SET_SUBTYPE=EQUATIONS_SET%SPECIFICATION(3)
    IS_HYDROSTATIC_PRESSURE_DEPENDENT_FIELD = EQUATIONS_SET_SUBTYPE/=EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE &
        & .AND. EQUATIONS_SET_SUBTYPE/=EQUATIONS_SET_COMPRESSIBLE_ACTIVECONTRACTION_SUBTYPE &
        & .AND. EQUATIONS_SET_SUBTYPE/=EQUATIONS_SET_MEMBRANE_SUBTYPE &
        & .AND. EQUATIONS_SET_SUBTYPE/=EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE &
        & .AND. EQUATIONS_SET_SUBTYPE/=EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_STATIC_INRIA_SUBTYPE &
        & .AND. EQUATIONS_SET_SUBTYPE/=EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_HOLMES_MOW_SUBTYPE &
        & .AND. EQUATIONS_SET_SUBTYPE/=EQUATIONS_SET_ELASTICITY_FLUID_PRES_HOLMES_MOW_ACTIVE_SUBTYPE &
        & .AND. EQUATIONS_SET_SUBTYPE/=EQUATIONS_SET_NEARLY_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE &
        & .AND. EQUATIONS_SET_SUBTYPE/=EQUATIONS_SET_COMPRESSIBLE_RATE_BASED_SMOOTH_MODEL_SUBTYPE &
        & .AND. EQUATIONS_SET_SUBTYPE/=EQUATIONS_SET_COMPRESSIBLE_RATE_BASED_GROWTH_MODEL_SUBTYPE

    NULLIFY(coordinateSystem)
    CALL EquationsSet_CoordinateSystemGet(EQUATIONS_SET,coordinateSystem,err,error,*999)
    CALL CoordinateSystem_DimensionGet(coordinateSystem,numberOfDimensions,err,error,*999)

    IF(IS_HYDROSTATIC_PRESSURE_DEPENDENT_FIELD) THEN
      NUMBER_OF_COMPONENTS = numberOfDimensions + 1
    ELSE
      NUMBER_OF_COMPONENTS = numberOfDimensions
    ENDIF

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET_SUBTYPE)
      CASE(EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE,EQUATIONS_SET_MEMBRANE_SUBTYPE, &
        & EQUATIONS_SET_MOONEY_RIVLIN_ACTIVECONTRACTION_SUBTYPE, &
        & EQUATIONS_SET_STVENANT_KIRCHOFF_ACTIVECONTRACTION_SUBTYPE, &
        & EQUATIONS_SET_ISOTROPIC_EXPONENTIAL_SUBTYPE,EQUATIONS_SET_TRANSVERSE_ISOTROPIC_ACTIVE_SUBTYPE, &
        & EQUATIONS_SET_TRANS_ISOTROPIC_ACTIVE_TRANSITION_SUBTYPE, &
        & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_EXPONENTIAL_SUBTYPE,EQUATIONS_SET_TRANSVERSE_ISOTROPIC_POLYNOMIAL_SUBTYPE, &
        & EQUATIONS_SET_ANISOTROPIC_POLYNOMIAL_SUBTYPE,EQUATIONS_SET_ANISOTROPIC_POLYNOMIAL_ACTIVE_SUBTYPE, &
        & EQUATIONS_SET_ORTHOTROPIC_MATERIAL_COSTA_SUBTYPE,EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE, &
        & EQUATIONS_SET_COMPRESSIBLE_ACTIVECONTRACTION_SUBTYPE,&
        & EQUATIONS_SET_ACTIVECONTRACTION_SUBTYPE, EQUATIONS_SET_NO_SUBTYPE, &
        & EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
        & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
        & EQUATIONS_SET_ORTHOTROPIC_MATERIAL_HOLZAPFEL_OGDEN_SUBTYPE, &
        & EQUATIONS_SET_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE,EQUATIONS_SET_NEARLY_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE, & 
        & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE,EQUATIONS_SET_TRANSVERSE_ISOTROPIC_GUCCIONE_SUBTYPE, &
        & EQUATIONS_SET_GUCCIONE_ACTIVECONTRACTION_SUBTYPE, EQUATIONS_SET_REFERENCE_STATE_TRANSVERSE_GUCCIONE_SUBTYPE, &
        & EQUATIONS_SET_CONSTITUTIVE_LAW_IN_CELLML_EVALUATE_SUBTYPE, &
        & EQUATIONS_SET_CONSTITUTIVE_AND_GROWTH_LAW_IN_CELLML_SUBTYPE, &
        & EQUATIONS_SET_GROWTH_LAW_IN_CELLML_SUBTYPE, &
        & EQUATIONS_SET_MR_AND_GROWTH_LAW_IN_CELLML_SUBTYPE, &
        & EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_STATIC_INRIA_SUBTYPE, &
        & EQUATIONS_SET_ELASTICITY_FLUID_PRES_HOLMES_MOW_ACTIVE_SUBTYPE, &
        & EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_HOLMES_MOW_SUBTYPE,EQUATIONS_SET_TRANSVERSE_ISOTROPIC_HUMPHREY_YIN_SUBTYPE, &
        & EQUATIONS_SET_STANDARD_MONODOMAIN_ELASTICITY_SUBTYPE,EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE, &
        & EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE,EQUATIONS_SET_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE, &
        & EQUATIONS_SET_HOLZAPFEL_OGDEN_ACTIVECONTRACTION_SUBTYPE, &
        & EQUATIONS_SET_RATE_BASED_SMOOTH_MODEL_SUBTYPE,EQUATIONS_SET_COMPRESSIBLE_RATE_BASED_SMOOTH_MODEL_SUBTYPE, &
        & EQUATIONS_SET_RATE_BASED_GROWTH_MODEL_SUBTYPE,EQUATIONS_SET_COMPRESSIBLE_RATE_BASED_GROWTH_MODEL_SUBTYPE)
        SELECT CASE(EQUATIONS_SET_SETUP%SETUP_TYPE)
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            !Default to FEM solution method
            CALL FiniteElasticity_EquationsSetSolutionMethodSet(EQUATIONS_SET,EQUATIONS_SET_FEM_SOLUTION_METHOD, &
              & err,error,*999)
            CALL EquationsSet_LabelSet(EQUATIONS_SET,"Finite elasticity equations set",err,error,*999)
            IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE) THEN
              !setup equations set field to store number of fluid compartments
              EQUATIONS_SET_FIELD_NUMBER_OF_VARIABLES = 1
              EQUATIONS_SET_FIELD_NUMBER_OF_COMPONENTS = 2
              EQUATIONS_EQUATIONS_SET_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD
              IF(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_AUTO_CREATED) THEN
                !Create the auto created equations set field
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION, &
                  & EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,err,error,*999)
                EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
                CALL FIELD_LABEL_SET(EQUATIONS_SET_FIELD_FIELD,"Equations Set Field",err,error,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,FIELD_GENERAL_TYPE,&
                  & err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,&
                  & FIELD_INDEPENDENT_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_SET(EQUATIONS_SET_FIELD_FIELD, &
                  & EQUATIONS_SET_FIELD_NUMBER_OF_VARIABLES,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,&
                  & [FIELD_U_VARIABLE_TYPE],err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_INTG_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,&
                  & FIELD_U_VARIABLE_TYPE,EQUATIONS_SET_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
              ELSE
                !Check the user specified field
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,EQUATIONS_SET_FIELD_NUMBER_OF_VARIABLES, &
                  & err,error,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_INTG_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, & 
                  & EQUATIONS_SET_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
              ENDIF
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)THEN
              IF(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_AUTO_CREATED) THEN
                CALL FIELD_CREATE_FINISH(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,err,error,*999)
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,&
                  & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, 1, 1_INTG, ERR, ERROR, *999)
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,&
                  & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, 2, 1_INTG, ERR, ERROR, *999)
              ENDIF
            ENDIF
!!TODO: Check valid setup
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a finite elasticity equation."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
          !\todo Check dimension of geometric field
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            !Check whether a fibre field is required, and if so, make sure it has been set
            SELECT CASE(EQUATIONS_SET_SUBTYPE)
            CASE(EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE, &
              & EQUATIONS_SET_MOONEY_RIVLIN_ACTIVECONTRACTION_SUBTYPE, &
              & EQUATIONS_SET_STVENANT_KIRCHOFF_ACTIVECONTRACTION_SUBTYPE, &
              & EQUATIONS_SET_ISOTROPIC_EXPONENTIAL_SUBTYPE, &
              & EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE,&
              & EQUATIONS_SET_COMPRESSIBLE_ACTIVECONTRACTION_SUBTYPE,&
              & EQUATIONS_SET_NO_SUBTYPE, &
              & EQUATIONS_SET_MEMBRANE_SUBTYPE, &
              & EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_STATIC_INRIA_SUBTYPE, &
              & EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_HOLMES_MOW_SUBTYPE, &
              & EQUATIONS_SET_ELASTICITY_FLUID_PRES_HOLMES_MOW_ACTIVE_SUBTYPE, &
              & EQUATIONS_SET_NEARLY_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE, & 
              & EQUATIONS_SET_MR_AND_GROWTH_LAW_IN_CELLML_SUBTYPE)
              ! pass, fibre field isn't required as the constitutive relation is isotropic
            CASE(EQUATIONS_SET_ORTHOTROPIC_MATERIAL_COSTA_SUBTYPE, &
              & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_EXPONENTIAL_SUBTYPE, &
              & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_POLYNOMIAL_SUBTYPE, &
              & EQUATIONS_SET_ANISOTROPIC_POLYNOMIAL_SUBTYPE,EQUATIONS_SET_ANISOTROPIC_POLYNOMIAL_ACTIVE_SUBTYPE, &
              & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_ACTIVE_SUBTYPE, &
              & EQUATIONS_SET_TRANS_ISOTROPIC_ACTIVE_TRANSITION_SUBTYPE, &
              & EQUATIONS_SET_ACTIVECONTRACTION_SUBTYPE, &
              & EQUATIONS_SET_ORTHOTROPIC_MATERIAL_HOLZAPFEL_OGDEN_SUBTYPE, &
              & EQUATIONS_SET_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE, &
              & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_GUCCIONE_SUBTYPE, &
              & EQUATIONS_SET_REFERENCE_STATE_TRANSVERSE_GUCCIONE_SUBTYPE, &
              & EQUATIONS_SET_GUCCIONE_ACTIVECONTRACTION_SUBTYPE, &
              & EQUATIONS_SET_CONSTITUTIVE_LAW_IN_CELLML_EVALUATE_SUBTYPE, &
              & EQUATIONS_SET_CONSTITUTIVE_AND_GROWTH_LAW_IN_CELLML_SUBTYPE, &
              & EQUATIONS_SET_GROWTH_LAW_IN_CELLML_SUBTYPE, &
              & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE, &
              & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
              & EQUATIONS_SET_STANDARD_MONODOMAIN_ELASTICITY_SUBTYPE, &
              & EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE,EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, &
              & EQUATIONS_SET_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE,EQUATIONS_SET_ACTIVE_STRAIN_SUBTYPE, &
                  & EQUATIONS_SET_MULTISCALE_ACTIVE_STRAIN_SUBTYPE,EQUATIONS_SET_1D3D_MONODOMAIN_ACTIVE_STRAIN_SUBTYPE, &
              & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_HUMPHREY_YIN_SUBTYPE, &
              & EQUATIONS_SET_HOLZAPFEL_OGDEN_ACTIVECONTRACTION_SUBTYPE, &
              & EQUATIONS_SET_RATE_BASED_SMOOTH_MODEL_SUBTYPE,EQUATIONS_SET_COMPRESSIBLE_RATE_BASED_SMOOTH_MODEL_SUBTYPE, &
              & EQUATIONS_SET_RATE_BASED_GROWTH_MODEL_SUBTYPE,EQUATIONS_SET_COMPRESSIBLE_RATE_BASED_GROWTH_MODEL_SUBTYPE)
              IF(.NOT.ASSOCIATED(EQUATIONS_SET%GEOMETRY%FIBRE_FIELD)) CALL FlagError( &
                & "Finite elascitiy equations require a fibre field.",err,error,*999)
            CASE DEFAULT
              LOCAL_ERROR="The third equations set specification of "// &
                & TRIM(NumberToVString(EQUATIONS_SET_SUBTYPE,"*",err,error))// &
                & " is invalid for a finite elasticity equation."
              CALL FlagError(LOCAL_ERROR,err,error,*999)
            END SELECT
            IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE .OR. &
                & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE .OR. &
                & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE .OR. &
                & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE) THEN
              ! Set up mesh displacement and equations set field info for elasticity Darcy problems
              FIELD_VARIABLE=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%ptr
              CALL Field_ParameterSetEnsureCreated(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_MESH_DISPLACEMENT_SET_TYPE,err,error,*999)
              IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMPARTMENT_DARCY_SUBTYPE .OR. &
                 EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE) THEN
                !Create the equations set field for multi-compartment Darcy
                EQUATIONS_SET_FIELD_NUMBER_OF_COMPONENTS = 2

                EQUATIONS_EQUATIONS_SET_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD
                EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD

                IF(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_AUTO_CREATED) THEN
                  CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,err,error,*999)
                  CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,&
                    & GEOMETRIC_DECOMPOSITION,err,error,*999)
                  CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,& 
                    & EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 1,GEOMETRIC_COMPONENT_NUMBER,err,error,*999)
                  DO component_idx = 1, EQUATIONS_SET_FIELD_NUMBER_OF_COMPONENTS
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET_AND_LOCK(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,component_idx,GEOMETRIC_COMPONENT_NUMBER,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,component_idx,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                  END DO

                  !Default the field scaling to that of the geometric field
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,GEOMETRIC_SCALING_TYPE, &
                    & err,error,*999)
                ELSE
                  !Do nothing
                ENDIF
              ENDIF
            END IF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            ! do nothing
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a linear diffusion equation."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
          SELECT CASE(EQUATIONS_SET_SUBTYPE)
          !-----------------------------------------------------------------------
          ! Dependent field setup for single-physics
          !-----------------------------------------------------------------------
          CASE(EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE, & 
            & EQUATIONS_SET_MOONEY_RIVLIN_ACTIVECONTRACTION_SUBTYPE, &
            & EQUATIONS_SET_STVENANT_KIRCHOFF_ACTIVECONTRACTION_SUBTYPE, &
            & EQUATIONS_SET_ISOTROPIC_EXPONENTIAL_SUBTYPE,EQUATIONS_SET_TRANSVERSE_ISOTROPIC_ACTIVE_SUBTYPE, &
            & EQUATIONS_SET_TRANS_ISOTROPIC_ACTIVE_TRANSITION_SUBTYPE, &
            & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_EXPONENTIAL_SUBTYPE,EQUATIONS_SET_TRANSVERSE_ISOTROPIC_POLYNOMIAL_SUBTYPE, &
            & EQUATIONS_SET_ANISOTROPIC_POLYNOMIAL_SUBTYPE,EQUATIONS_SET_ANISOTROPIC_POLYNOMIAL_ACTIVE_SUBTYPE, &
            & EQUATIONS_SET_ORTHOTROPIC_MATERIAL_COSTA_SUBTYPE,EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE,&
            & EQUATIONS_SET_COMPRESSIBLE_ACTIVECONTRACTION_SUBTYPE, &
            & EQUATIONS_SET_ACTIVECONTRACTION_SUBTYPE, EQUATIONS_SET_NO_SUBTYPE,EQUATIONS_SET_MEMBRANE_SUBTYPE, &
            & EQUATIONS_SET_ORTHOTROPIC_MATERIAL_HOLZAPFEL_OGDEN_SUBTYPE,EQUATIONS_SET_TRANSVERSE_ISOTROPIC_GUCCIONE_SUBTYPE, &
            & EQUATIONS_SET_GUCCIONE_ACTIVECONTRACTION_SUBTYPE, EQUATIONS_SET_REFERENCE_STATE_TRANSVERSE_GUCCIONE_SUBTYPE, &
            & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_HUMPHREY_YIN_SUBTYPE, &
            & EQUATIONS_SET_STANDARD_MONODOMAIN_ELASTICITY_SUBTYPE,EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE, &
            & EQUATIONS_SET_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE,EQUATIONS_SET_NEARLY_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE, &
            & EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, &
            & EQUATIONS_SET_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE, &
            & EQUATIONS_SET_HOLZAPFEL_OGDEN_ACTIVECONTRACTION_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                !Create the auto created dependent field
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET%DEPENDENT% &
                  & DEPENDENT_FIELD,err,error,*999)
                CALL FIELD_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,"Dependent Field",err,error,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_GEOMETRIC_GENERAL_TYPE,err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,err,error,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & err,error,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY% &
                  & GEOMETRIC_FIELD,err,error,*999)
                IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_ACTIVE_STRAIN_SUBTYPE) THEN
                  CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,3,err,error,*999)
                ELSE IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTISCALE_ACTIVE_STRAIN_SUBTYPE) THEN
                  CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,4,err,error,*999)
                ELSE
                  DEPENDENT_FIELD_NUMBER_OF_VARIABLES=2
                  CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & DEPENDENT_FIELD_NUMBER_OF_VARIABLES,err,error,*999)
                  CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,[FIELD_U_VARIABLE_TYPE, &
                    & FIELD_DELUDELN_VARIABLE_TYPE],err,error,*999)
                  CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & "U",err,error,*999)
                  CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                    & "del U/del n",err,error,*999)
                END IF
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & numberOfDimensions,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_COMPONENTS,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,NUMBER_OF_COMPONENTS,err,error,*999)
                IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_ACTIVE_STRAIN_SUBTYPE) THEN
                  CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                  CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_V_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
                ELSE IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTISCALE_ACTIVE_STRAIN_SUBTYPE) THEN
                  CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                  CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_V_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
                  CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                  CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                    & FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_U1_VARIABLE_TYPE,2,err,error,*999)
                END IF

                !Default to the geometric interpolation setup
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                DO component_idx=1,numberOfDimensions
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_ACTIVE_STRAIN_SUBTYPE) THEN
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_V_VARIABLE_TYPE,component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  ELSE IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTISCALE_ACTIVE_STRAIN_SUBTYPE) THEN
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_V_VARIABLE_TYPE,component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  END IF
                ENDDO !component_idx

                IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTISCALE_ACTIVE_STRAIN_SUBTYPE) THEN
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_U1_VARIABLE_TYPE,1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_U1_VARIABLE_TYPE,2,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                END IF

                IF(IS_HYDROSTATIC_PRESSURE_DEPENDENT_FIELD) THEN
!kmith :09.06.09 - Do we need this ?
                  !Set the hydrostatic component to that of the first geometric component
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_COMPONENTS,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                    & NUMBER_OF_COMPONENTS,GEOMETRIC_MESH_COMPONENT,err,error,*999)
!kmith
                ENDIF
                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  !Set the displacement components to node based interpolation
                  DO component_idx=1,numberOfDimensions
!                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
!                      & component_idx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
!                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
!                      & FIELD_DELUDELN_VARIABLE_TYPE,component_idx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & component_idx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_DELUDELN_VARIABLE_TYPE,component_idx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  ENDDO !component_idx
                  IF(IS_HYDROSTATIC_PRESSURE_DEPENDENT_FIELD) THEN
                    !Set the hydrostatic pressure component to element based interpolation
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & NUMBER_OF_COMPONENTS,FIELD_ELEMENT_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_DELUDELN_VARIABLE_TYPE,NUMBER_OF_COMPONENTS,FIELD_ELEMENT_BASED_INTERPOLATION,err,error,*999)
                  ENDIF
                  !Default the scaling to the geometric field scaling
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(LOCAL_ERROR,err,error,*999)
                END SELECT
              ELSE
                !Check the user specified field
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GEOMETRIC_GENERAL_TYPE,err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
                IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_ACTIVE_STRAIN_SUBTYPE) THEN
                  CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,3,err,error,*999)
                  CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE, &
                    & FIELD_V_VARIABLE_TYPE],err,error,*999)
                  CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                    & err,error,*999)
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                ELSE IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTISCALE_ACTIVE_STRAIN_SUBTYPE) THEN
                  CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,4,err,error,*999)
                  CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE, &
                    & FIELD_V_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE],err,error,*999)
                  CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,&
                    & err,error,*999)
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,&
                    & err,error,*999)
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                ELSE
                  CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,2,err,error,*999)
                  CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE],&
                    & err,error,*999)
                END IF
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & numberOfDimensions,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS, &
                  & err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,NUMBER_OF_COMPONENTS, &
                  & err,error,*999)
                IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_ACTIVE_STRAIN_SUBTYPE) THEN
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE, &
                    & numberOfDimensions,err,error,*999)
                ELSE IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTISCALE_ACTIVE_STRAIN_SUBTYPE) THEN
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE, &
                    & numberOfDimensions,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE, &
                    & 2,err,error,*999)
                END IF
                !Check that the pressure values set type is created here?? (second variable is a DELUDELN type, as checked above)
                !\todo: Decide whether these set_types (previous one as well) is to be created by user or automatically..
                IF(.not.ASSOCIATED(EQUATIONS_SET_SETUP%FIELD%VARIABLES(2)%PARAMETER_SETS% &
                  & SET_TYPE(FIELD_PRESSURE_VALUES_SET_TYPE)%ptr)) THEN
                    LOCAL_ERROR="Variable 2 of type "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%FIELD%VARIABLES(2)% &
                      & VARIABLE_TYPE,"*",err,error))//" does not have a pressure values set type associated."
                ENDIF
                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  DO component_idx=1,numberOfDimensions
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,component_idx, &
                      & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,component_idx, &
                      & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  ENDDO !component_idx
                CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(LOCAL_ERROR,err,error,*999)
                END SELECT
              ENDIF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                CALL FIELD_CREATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,err,error,*999)
              ENDIF
              IF(EQUATIONS_SET%specification(3)==EQUATIONS_SET_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE) THEN
                CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_PREVIOUS_VALUES_SET_TYPE,err,error,*999)
                CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,err,error,*999)
              ENDIF
            CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
                & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a finite elasticity equation"
              CALL FlagError(LOCAL_ERROR,err,error,*999)
            END SELECT

          CASE(EQUATIONS_SET_GROWTH_LAW_IN_CELLML_SUBTYPE)
            
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                !Create the auto created dependent field
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET%DEPENDENT% &
                  & DEPENDENT_FIELD,err,error,*999)
                CALL FIELD_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,"Dependent Field",err,error,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_GEOMETRIC_GENERAL_TYPE,err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,err,error,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & err,error,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY% &
                  & GEOMETRIC_FIELD,err,error,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,2,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,[FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DELUDELN_VARIABLE_TYPE],err,error,*999)
                CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & "Z",err,error,*999)
                CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & "Traction",err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & numberOfDimensions,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & numberOfDimensions,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
                
                !Default to the geometric interpolation setup
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                DO component_idx=1,numberOfDimensions
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                ENDDO !component_idx
                
                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  !Set the displacement components to node based interpolation
                  DO component_idx=1,numberOfDimensions
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & component_idx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_DELUDELN_VARIABLE_TYPE,component_idx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  ENDDO !component_idx

                  IF(IS_HYDROSTATIC_PRESSURE_DEPENDENT_FIELD) THEN
                    !Set the hydrostatic pressure component to element based interpolation
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & component_idx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_DELUDELN_VARIABLE_TYPE,component_idx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  ENDIF
                  !Default the scaling to the geometric field scaling
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%geometric_Field,GEOMETRIC_SCALING_TYPE,err,error,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%DEPENDENT%dependent_Field,GEOMETRIC_SCALING_TYPE,err,error,*999)
                CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(LOCAL_ERROR,err,error,*999)
                END SELECT
              ELSE
                !Check the user specified field
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GEOMETRIC_GENERAL_TYPE,err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,2,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE],&
                  & err,error,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & numberOfDimensions,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,numberOfDimensions, &
                  & err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,numberOfDimensions, &
                  & err,error,*999)
                !Check that the pressure values set type is created here?? (second variable is a DELUDELN type, as checked above)
                !\todo: Decide whether these set_types (previous one as well) is to be created by user or automatically..
                IF(.NOT.ASSOCIATED(EQUATIONS_SET_SETUP%FIELD%VARIABLES(2)%PARAMETER_SETS% &
                  & SET_TYPE(FIELD_PRESSURE_VALUES_SET_TYPE)%ptr)) THEN
                  LOCAL_ERROR="Variable 2 of type "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%FIELD%VARIABLES(2)% &
                    & VARIABLE_TYPE,"*",err,error))//" does not have a pressure values set type associated."
                ENDIF
                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  DO component_idx=1,numberOfDimensions
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,component_idx, &
                      & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,component_idx, &
                      & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  ENDDO !component_idx
                CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(LOCAL_ERROR,err,error,*999)
                END SELECT
              ENDIF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                CALL FIELD_CREATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,err,error,*999)
              ENDIF
            CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
                & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a finite elasticity equation"
              CALL FlagError(LOCAL_ERROR,err,error,*999)
            END SELECT

          CASE(EQUATIONS_SET_MR_AND_GROWTH_LAW_IN_CELLML_SUBTYPE)
            !--------------------------------------------------------------------------------------
            ! Dependent field setup for a code constitutive law with a growth law defined in CellML
            !--------------------------------------------------------------------------------------
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                !Create the auto created dependent field
                CALL Field_CreateStart(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET%DEPENDENT% &
                  & DEPENDENT_FIELD,err,error,*999)
                CALL Field_TypeSetAndLock(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_GEOMETRIC_GENERAL_TYPE,err,error,*999)
                CALL Field_DependentTypeSetAndLock(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
                CALL Field_MeshDecompositionGet(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,err,error,*999)
                CALL Field_MeshDecompositionSetAndLock(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & err,error,*999)
                CALL Field_GeometricFieldSetAndLock(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY% &
                  & GEOMETRIC_FIELD,err,error,*999)
                CALL Field_NumberOfVariablesSetAndLock(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,3,err,error,*999)
                CALL Field_VariableTypesSetAndLock(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,[FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,FIELD_U3_VARIABLE_TYPE],err,error,*999)
                CALL Field_VariableLabelSet(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & "U",err,error,*999)
                CALL Field_VariableLabelSet(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & "del U/del n",err,error,*999)
                CALL Field_VariableLabelSet(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U3_VARIABLE_TYPE, &
                  & "U3",err,error,*999)
                CALL Field_DimensionSetAndLock(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL Field_DimensionSetAndLock(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL Field_DimensionSetAndLock(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U3_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL Field_DataTypeSetAndLock(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL Field_DataTypeSetAndLock(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL Field_DataTypeSetAndLock(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U3_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL Field_NumberOfComponentsGet(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & numberOfDimensions,err,error,*999)
                CALL Field_NumberOfComponentsSetAndLock(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_COMPONENTS,err,error,*999)
                CALL Field_NumberOfComponentsSetAndLock(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,NUMBER_OF_COMPONENTS,err,error,*999)
                CALL Field_NumberOfComponentsSetAndLock(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U3_VARIABLE_TYPE, &
                  & numberOfDimensions,err,error,*999)
               
                !Default to the geometric interpolation setup
                DO component_idx=1,numberOfDimensions
                  CALL Field_ComponentMeshComponentGet(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U3_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                ENDDO !component_idx

                IF(IS_HYDROSTATIC_PRESSURE_DEPENDENT_FIELD) THEN
                  !Set the hydrostatic component to that of the first geometric component
                  CALL Field_ComponentMeshComponentGet(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_COMPONENTS,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                    & NUMBER_OF_COMPONENTS,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                ENDIF

                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  !Set the displacement components to node based interpolation, set the growth to Gauss point
                  DO component_idx=1,numberOfDimensions
                    CALL Field_ComponentInterpolationSetAndLock(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,component_idx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL Field_ComponentInterpolationSetAndLock(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_DELUDELN_VARIABLE_TYPE,component_idx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL Field_ComponentInterpolationSetAndLock(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_U3_VARIABLE_TYPE,component_idx,FIELD_GAUSS_POINT_BASED_INTERPOLATION,err,error,*999)
                  ENDDO !component_idx
                  
                  IF(IS_HYDROSTATIC_PRESSURE_DEPENDENT_FIELD) THEN
                    !Set the hydrostatic pressure component to element based interpolation
                    CALL Field_ComponentInterpolationSet(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & NUMBER_OF_COMPONENTS,FIELD_ELEMENT_BASED_INTERPOLATION,err,error,*999)
                    CALL Field_ComponentInterpolationSet(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_DELUDELN_VARIABLE_TYPE,NUMBER_OF_COMPONENTS,FIELD_ELEMENT_BASED_INTERPOLATION,err,error,*999)
                  ENDIF

                  !Default the scaling to the geometric field scaling
                  CALL Field_ScalingTypeGet(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                  CALL Field_ScalingTypeSet(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(LOCAL_ERROR,err,error,*999)
                END SELECT

              ELSE !EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED

                !Check the user specified field
                CALL Field_TypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_GEOMETRIC_GENERAL_TYPE,err,error,*999)
                CALL Field_DependentTypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
                CALL Field_NumberOfVariablesCheck(EQUATIONS_SET_SETUP%FIELD,4,err,error,*999)
                CALL Field_VariableTypesCheck(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_U3_VARIABLE_TYPE],err,error,*999)
                CALL Field_DimensionCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL Field_DimensionCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL Field_DimensionCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U3_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL Field_DataTypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL Field_DataTypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL Field_DataTypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U3_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL Field_NumberOfComponentsGet(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & numberOfDimensions,err,error,*999)
                CALL Field_NumberOfComponentsCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS, &
                  & err,error,*999)
                CALL Field_NumberOfComponentsCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & NUMBER_OF_COMPONENTS,err,error,*999)
                CALL Field_NumberOfComponentsCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U3_VARIABLE_TYPE,numberOfDimensions, &
                  & err,error,*999)

                !Check that the pressure values set type is created here?? (second variable is a DELUDELN type, as checked above)
                !\todo: Decide whether these set_types (previous one as well) is to be created by user or automatically..
                IF(.not.ASSOCIATED(EQUATIONS_SET_SETUP%FIELD%VARIABLES(2)%PARAMETER_SETS% &
                  & SET_TYPE(FIELD_PRESSURE_VALUES_SET_TYPE)%ptr)) THEN
                    LOCAL_ERROR="Variable 2 of type "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%FIELD%VARIABLES(2)% &
                      & VARIABLE_TYPE,"*",err,error))//" does not have a pressure values set type associated."
                ENDIF
                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  DO component_idx=1,numberOfDimensions
                    CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,component_idx, &
                      & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,component_idx, &
                      & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U3_VARIABLE_TYPE,component_idx, &
                      & FIELD_GAUSS_POINT_BASED_INTERPOLATION,err,error,*999)
                  ENDDO !component_idx
                CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(LOCAL_ERROR,err,error,*999)
                END SELECT
              ENDIF !EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                CALL Field_CreateFinish(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,err,error,*999)
              ENDIF
            CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
                & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a finite elasticity equation"
              CALL FlagError(LOCAL_ERROR,err,error,*999)
            END SELECT
            
          CASE(EQUATIONS_SET_CONSTITUTIVE_LAW_IN_CELLML_EVALUATE_SUBTYPE, &
            & EQUATIONS_SET_CONSTITUTIVE_AND_GROWTH_LAW_IN_CELLML_SUBTYPE)
            !-------------------------------------------------------------------------------
            ! Dependent field setup for elasticity evaluated in CellML
            !-------------------------------------------------------------------------------
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                !Create the auto created dependent field
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET%DEPENDENT% &
                  & DEPENDENT_FIELD,err,error,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_GEOMETRIC_GENERAL_TYPE,err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,err,error,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & err,error,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY% &
                  & GEOMETRIC_FIELD,err,error,*999)
                IF(numberOfDimensions==3) THEN
                  NUMBER_OF_COMPONENTS_2 = 6
                ELSE IF(numberOfDimensions==2) THEN
                  NUMBER_OF_COMPONENTS_2 = 3
                ELSE
                  CALL FlagError("Only 2 and 3 dimensional problems are implemented at the moment",err,error,*999)
                ENDIF !numberOfDimensions
                IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CONSTITUTIVE_AND_GROWTH_LAW_IN_CELLML_SUBTYPE) THEN
                  CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,5,err,error,*999)
                  CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,[FIELD_U_VARIABLE_TYPE, &
                    & FIELD_DELUDELN_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE,FIELD_U2_VARIABLE_TYPE,FIELD_U3_VARIABLE_TYPE], &
                    & err,error,*999)
                ELSE
                  CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,4,err,error,*999)
                  CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,[FIELD_U_VARIABLE_TYPE, &
                    & FIELD_DELUDELN_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE,FIELD_U2_VARIABLE_TYPE],err,error,*999)
                ENDIF
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U2_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U2_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & numberOfDimensions,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_COMPONENTS,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,NUMBER_OF_COMPONENTS,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                  & NUMBER_OF_COMPONENTS_2,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U2_VARIABLE_TYPE, &
                  & NUMBER_OF_COMPONENTS_2,err,error,*999)
                
                IF(EQUATIONS_SET%specification(3)==EQUATIONS_SET_CONSTITUTIVE_AND_GROWTH_LAW_IN_CELLML_SUBTYPE) THEN
                  CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U3_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                  CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U3_VARIABLE_TYPE, &
                    & FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U3_VARIABLE_TYPE, &
                    & numberOfDimensions,err,error,*999)
                ENDIF

                !Default to the geometric interpolation setup
                DO component_idx=1,numberOfDimensions
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                ENDDO !component_idx

                IF(IS_HYDROSTATIC_PRESSURE_DEPENDENT_FIELD) THEN
!kmith :09.06.09 - Do we need this ?
                  !Set the hydrostatic component to that of the first geometric component
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_COMPONENTS,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                    & NUMBER_OF_COMPONENTS,GEOMETRIC_MESH_COMPONENT,err,error,*999)
!kmith
                ENDIF

                !Set the stress and strain components to that of the first geometric component
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                DO component_idx=1,NUMBER_OF_COMPONENTS_2
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U2_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  IF(EQUATIONS_SET%specification(3)==EQUATIONS_SET_CONSTITUTIVE_AND_GROWTH_LAW_IN_CELLML_SUBTYPE) THEN
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U3_VARIABLE_TYPE, &
                      & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  ENDIF
                ENDDO !component_idx

                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  !Set the displacement components to node based interpolation
                  DO component_idx=1,numberOfDimensions
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & component_idx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_DELUDELN_VARIABLE_TYPE,component_idx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  ENDDO !component_idx
                  
                  IF(IS_HYDROSTATIC_PRESSURE_DEPENDENT_FIELD) THEN
                    !Set the hydrostatic pressure component to element based interpolation
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & NUMBER_OF_COMPONENTS,FIELD_ELEMENT_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_DELUDELN_VARIABLE_TYPE,NUMBER_OF_COMPONENTS,FIELD_ELEMENT_BASED_INTERPOLATION,err,error,*999)
                  ENDIF

                  !Set the stress and strain components to gauss point interpolation
                  DO component_idx=1,NUMBER_OF_COMPONENTS_2
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_U1_VARIABLE_TYPE,component_idx,FIELD_GAUSS_POINT_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_U2_VARIABLE_TYPE,component_idx,FIELD_GAUSS_POINT_BASED_INTERPOLATION,err,error,*999)
                    IF(EQUATIONS_SET%specification(3)==EQUATIONS_SET_CONSTITUTIVE_AND_GROWTH_LAW_IN_CELLML_SUBTYPE) THEN
                      CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                        & FIELD_U3_VARIABLE_TYPE,component_idx,FIELD_GAUSS_POINT_BASED_INTERPOLATION,err,error,*999)
                    ENDIF
                  ENDDO !component_idx

                  !Default the scaling to the geometric field scaling
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(LOCAL_ERROR,err,error,*999)
                END SELECT

              ELSE !EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED

                !Check the user specified field
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GEOMETRIC_GENERAL_TYPE,err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
                IF(numberOfDimensions==3) THEN
                  NUMBER_OF_COMPONENTS_2 = 6
                ELSE IF(numberOfDimensions==2) THEN
                  NUMBER_OF_COMPONENTS_2 = 3
                ELSE
                  CALL FlagError("Only 2 and 3 dimensional problems are implemented at the moment",err,error,*999)
                ENDIF !numberOfDimensions
                IF(EQUATIONS_SET%specification(3)==EQUATIONS_SET_CONSTITUTIVE_AND_GROWTH_LAW_IN_CELLML_SUBTYPE) THEN
                  CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,5,err,error,*999)
                  CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE, &
                    & FIELD_U1_VARIABLE_TYPE,FIELD_U2_VARIABLE_TYPE,FIELD_U3_VARIABLE_TYPE],err,error,*999)
                ELSE
                  CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,4,err,error,*999)
                  CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE, &
                    & FIELD_U1_VARIABLE_TYPE,FIELD_U2_VARIABLE_TYPE],err,error,*999)
                ENDIF
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U2_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U2_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & numberOfDimensions,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS, &
                  & err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & NUMBER_OF_COMPONENTS,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE,NUMBER_OF_COMPONENTS_2, &
                  & err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U2_VARIABLE_TYPE,NUMBER_OF_COMPONENTS_2, &
                  & err,error,*999)
                IF(EQUATIONS_SET%specification(3)==EQUATIONS_SET_CONSTITUTIVE_AND_GROWTH_LAW_IN_CELLML_SUBTYPE) THEN
                  CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U3_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                    & err,error,*999)
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U3_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U3_VARIABLE_TYPE,numberOfDimensions, &
                    & err,error,*999)
                ENDIF

                !Check that the pressure values set type is created here?? (second variable is a DELUDELN type, as checked above)
                !\todo: Decide whether these set_types (previous one as well) is to be created by user or automatically..
                IF(.not.ASSOCIATED(EQUATIONS_SET_SETUP%FIELD%VARIABLES(2)%PARAMETER_SETS% &
                  & SET_TYPE(FIELD_PRESSURE_VALUES_SET_TYPE)%ptr)) THEN
                    LOCAL_ERROR="Variable 2 of type "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%FIELD%VARIABLES(2)% &
                      & VARIABLE_TYPE,"*",err,error))//" does not have a pressure values set type associated."
                ENDIF
                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  DO component_idx=1,numberOfDimensions
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,component_idx, &
                      & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,component_idx, &
                      & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    IF(EQUATIONS_SET%specification(3)==EQUATIONS_SET_CONSTITUTIVE_AND_GROWTH_LAW_IN_CELLML_SUBTYPE) THEN
                      CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U3_VARIABLE_TYPE,component_idx, &
                        & FIELD_GAUSS_POINT_BASED_INTERPOLATION,err,error,*999)
                    ENDIF
                  ENDDO !component_idx
                  DO component_idx=1,NUMBER_OF_COMPONENTS_2
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE,component_idx, &
                      & FIELD_GAUSS_POINT_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U2_VARIABLE_TYPE,component_idx, &
                      & FIELD_GAUSS_POINT_BASED_INTERPOLATION,err,error,*999)
                  ENDDO !component_idx

                CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(LOCAL_ERROR,err,error,*999)
                END SELECT
              ENDIF !EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                CALL FIELD_CREATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,err,error,*999)
              ENDIF
            CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
                & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a finite elasticity equation"
              CALL FlagError(LOCAL_ERROR,err,error,*999)
            END SELECT

          CASE(EQUATIONS_SET_RATE_BASED_SMOOTH_MODEL_SUBTYPE,EQUATIONS_SET_COMPRESSIBLE_RATE_BASED_SMOOTH_MODEL_SUBTYPE, &
            & EQUATIONS_SET_RATE_BASED_GROWTH_MODEL_SUBTYPE,EQUATIONS_SET_COMPRESSIBLE_RATE_BASED_GROWTH_MODEL_SUBTYPE)
            !Rubin rate based formulations
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                !Create the auto created dependent field
                CALL Field_CreateStart(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET%DEPENDENT% &
                  & DEPENDENT_FIELD,err,error,*999)
                CALL Field_LabelSet(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,"Dependent Field",err,error,*999)
                CALL Field_TypeSetAndLock(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_GEOMETRIC_GENERAL_TYPE,err,error,*999)
                CALL Field_DependentTypeSetAndLock(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
                CALL Field_MeshDecompositionGet(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,err,error,*999)
                CALL Field_MeshDecompositionSetAndLock(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & err,error,*999)
                CALL Field_GeometricFieldSetAndLock(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY% &
                  & GEOMETRIC_FIELD,err,error,*999)
                SELECT CASE(EQUATIONS_SET_SUBTYPE)
                CASE(EQUATIONS_SET_RATE_BASED_SMOOTH_MODEL_SUBTYPE, &
                  & EQUATIONS_SET_COMPRESSIBLE_RATE_BASED_SMOOTH_MODEL_SUBTYPE)
                  CALL Field_NumberOfVariablesSetAndLock(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,3,err,error,*999)
                  CALL Field_VariableTypesSetAndLock(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,[FIELD_U_VARIABLE_TYPE, &
                    & FIELD_DELUDELN_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE],err,error,*999)
                  CALL Field_VariableLabelSet(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & "U",err,error,*999)
                  CALL Field_VariableLabelSet(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                    & "del U/del n",err,error,*999)
                  CALL Field_VariableLabelSet(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                    & "U1",err,error,*999)
                  CALL Field_DimensionSetAndLock(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                  CALL Field_DimensionSetAndLock(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                  CALL Field_DimensionSetAndLock(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                  CALL Field_DataTypeSetAndLock(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_DP_TYPE,err,error,*999)
                  CALL Field_DataTypeSetAndLock(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                    & FIELD_DP_TYPE,err,error,*999)
                  CALL Field_DataTypeSetAndLock(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                    & FIELD_DP_TYPE,err,error,*999)
                CASE( EQUATIONS_SET_RATE_BASED_GROWTH_MODEL_SUBTYPE, &
                  & EQUATIONS_SET_COMPRESSIBLE_RATE_BASED_GROWTH_MODEL_SUBTYPE)
                  CALL Field_NumberOfVariablesSetAndLock(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,4,err,error,*999)
                  CALL Field_VariableTypesSetAndLock(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,[FIELD_U_VARIABLE_TYPE, &
                    & FIELD_DELUDELN_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE,FIELD_U2_VARIABLE_TYPE],err,error,*999)
                  CALL Field_VariableLabelSet(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & "U",err,error,*999)
                  CALL Field_VariableLabelSet(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                    & "del U/del n",err,error,*999)
                  CALL Field_VariableLabelSet(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                    & "U1",err,error,*999)
                  CALL Field_VariableLabelSet(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                    & "U2",err,error,*999)
                  CALL Field_DimensionSetAndLock(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                  CALL Field_DimensionSetAndLock(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                  CALL Field_DimensionSetAndLock(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                  CALL Field_DimensionSetAndLock(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U2_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                  CALL Field_DataTypeSetAndLock(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_DP_TYPE,err,error,*999)
                  CALL Field_DataTypeSetAndLock(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                    & FIELD_DP_TYPE,err,error,*999)
                  CALL Field_DataTypeSetAndLock(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                    & FIELD_DP_TYPE,err,error,*999)
                  CALL Field_DataTypeSetAndLock(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U2_VARIABLE_TYPE, &
                    & FIELD_DP_TYPE,err,error,*999)
                CASE DEFAULT
                  localError="The equations set subtype of "//TRIM(NumberToVString(EQUATIONS_SET_SUBTYPE,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
                CALL Field_ComponentMeshComponentGet(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                SELECT CASE(EQUATIONS_SET_SUBTYPE)
                CASE(EQUATIONS_SET_RATE_BASED_SMOOTH_MODEL_SUBTYPE, &
                  & EQUATIONS_SET_COMPRESSIBLE_RATE_BASED_SMOOTH_MODEL_SUBTYPE)
                  CALL Field_NumberOfComponentsSetAndLock(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS,err,error,*999)
                  CALL Field_NumberOfComponentsSetAndLock(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_DELUDELN_VARIABLE_TYPE,NUMBER_OF_COMPONENTS,err,error,*999)
                  CALL Field_NumberOfComponentsSetAndLock(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_U1_VARIABLE_TYPE,1+NUMBER_OF_VOIGT(numberOfDimensions),err,error,*999)
                CASE(EQUATIONS_SET_RATE_BASED_GROWTH_MODEL_SUBTYPE, &
                  & EQUATIONS_SET_COMPRESSIBLE_RATE_BASED_GROWTH_MODEL_SUBTYPE)
                  CALL Field_NumberOfComponentsSetAndLock(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS,err,error,*999)
                  CALL Field_NumberOfComponentsSetAndLock(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_DELUDELN_VARIABLE_TYPE,NUMBER_OF_COMPONENTS,err,error,*999)
                  CALL Field_NumberOfComponentsSetAndLock(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_U1_VARIABLE_TYPE,1+NUMBER_OF_VOIGT(numberOfDimensions),err,error,*999)
                  CALL Field_NumberOfComponentsSetAndLock(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_U2_VARIABLE_TYPE,numberOfDimensions*numberOfDimensions,err,error,*999)
                CASE DEFAULT
                  localError="The equations set subtype of "//TRIM(NumberToVString(EQUATIONS_SET_SUBTYPE,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
                !Default to the geometric interpolation setup
                CALL Field_ComponentMeshComponentGet(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                DO component_idx=1,NUMBER_OF_COMPONENTS
                  CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                ENDDO !component_idx
                IF(IS_HYDROSTATIC_PRESSURE_DEPENDENT_FIELD) THEN
                  !Set the hydrostatic component to that of the first geometric component
                  CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & numberOfDimensions+1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                    & numberOfDimensions+1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                END IF
                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  !Set the displacement components to node based interpolation
                  DO component_idx=1,numberOfDimensions
                    CALL Field_ComponentInterpolationSet(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & component_idx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL Field_ComponentInterpolationSet(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_DELUDELN_VARIABLE_TYPE,component_idx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  ENDDO !component_idx
                  IF(IS_HYDROSTATIC_PRESSURE_DEPENDENT_FIELD) THEN
                    !Set the hydrostatic pressure component to element based interpolation
                    CALL Field_ComponentInterpolationSet(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & numberOfDimensions+1,FIELD_ELEMENT_BASED_INTERPOLATION,err,error,*999)
                    CALL Field_ComponentInterpolationSet(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_DELUDELN_VARIABLE_TYPE,numberOfDimensions+1,FIELD_ELEMENT_BASED_INTERPOLATION,err,error,*999)
                  ENDIF
                CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(LOCAL_ERROR,err,error,*999)
                END SELECT
                SELECT CASE(EQUATIONS_SET_SUBTYPE)
                CASE(EQUATIONS_SET_RATE_BASED_SMOOTH_MODEL_SUBTYPE, &
                  & EQUATIONS_SET_COMPRESSIBLE_RATE_BASED_SMOOTH_MODEL_SUBTYPE)
                  !Set the U1 variable components
                  DO componentIdx=1,1+NUMBER_OF_VOIGT(numberOfDimensions)
                    CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                      & componentIdx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                    CALL Field_ComponentInterpolationSet(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                      & componentIdx,FIELD_GAUSS_POINT_BASED_INTERPOLATION,err,error,*999)
                  ENDDO !componentIdx
                CASE(EQUATIONS_SET_RATE_BASED_GROWTH_MODEL_SUBTYPE, &
                  & EQUATIONS_SET_COMPRESSIBLE_RATE_BASED_GROWTH_MODEL_SUBTYPE)
                  !Set the U1 variable components
                  DO componentIdx=1,1+NUMBER_OF_VOIGT(numberOfDimensions)
                    CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                      & componentIdx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                    CALL Field_ComponentInterpolationSet(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                      & componentIdx,FIELD_GAUSS_POINT_BASED_INTERPOLATION,err,error,*999)
                  ENDDO !componentIdx
                  !Set the U2 variable components
                  DO componentIdx=1,numberOfDimensions*numberOfDimensions
                    CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U2_VARIABLE_TYPE, &
                      & componentIdx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                    CALL Field_ComponentInterpolationSet(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U2_VARIABLE_TYPE, &
                      & componentIdx,FIELD_GAUSS_POINT_BASED_INTERPOLATION,err,error,*999)
                  ENDDO !componentIdx
                CASE DEFAULT
                  localError="The equations set subtype of "//TRIM(NumberToVString(EQUATIONS_SET_SUBTYPE,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
                !Default the scaling to the geometric field scaling
                CALL Field_ScalingTypeGet(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                CALL Field_ScalingTypeSet(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
              ELSE
                !Check the user specified field
                CALL Field_TypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_GEOMETRIC_GENERAL_TYPE,err,error,*999)
                CALL Field_DependentTypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
                SELECT CASE(EQUATIONS_SET_SUBTYPE)
                CASE(EQUATIONS_SET_RATE_BASED_SMOOTH_MODEL_SUBTYPE, &
                  & EQUATIONS_SET_COMPRESSIBLE_RATE_BASED_SMOOTH_MODEL_SUBTYPE)
                  CALL Field_NumberOfVariablesCheck(EQUATIONS_SET_SETUP%FIELD,3,err,error,*999)
                  CALL Field_VariableTypesCheck(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE, &
                    & FIELD_U1_VARIABLE_TYPE],err,error,*999)
                  CALL Field_DimensionCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                    & err,error,*999)
                  CALL Field_DimensionCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                    & err,error,*999)
                  CALL Field_DimensionCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                    & err,error,*999)
                  CALL Field_DataTypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                  CALL Field_DataTypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                  CALL Field_DataTypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CASE( EQUATIONS_SET_RATE_BASED_GROWTH_MODEL_SUBTYPE, &
                  & EQUATIONS_SET_COMPRESSIBLE_RATE_BASED_GROWTH_MODEL_SUBTYPE)
                  CALL Field_NumberOfVariablesCheck(EQUATIONS_SET_SETUP%FIELD,4,err,error,*999)
                  CALL Field_VariableTypesCheck(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE, &
                    & FIELD_U1_VARIABLE_TYPE,FIELD_U2_VARIABLE_TYPE],err,error,*999)
                  CALL Field_DimensionCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                    & err,error,*999)
                  CALL Field_DimensionCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                    & err,error,*999)
                  CALL Field_DimensionCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                    & err,error,*999)
                  CALL Field_DimensionCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U2_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                    & err,error,*999)
                  CALL Field_DataTypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                  CALL Field_DataTypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                  CALL Field_DataTypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                  CALL Field_DataTypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U2_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CASE DEFAULT
                  localError="The equations set subtype of "//TRIM(NumberToVString(EQUATIONS_SET_SUBTYPE,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
                CALL Field_NumberOfComponentsCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS, &
                  & err,error,*999)
                CALL Field_NumberOfComponentsCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,NUMBER_OF_COMPONENTS, &
                  & err,error,*999)
                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  DO componentIdx=1,numberOfDimensions
                    CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,componentIdx, &
                      & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,componentIdx, &
                      & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  ENDDO !componentIdx
                CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(LOCAL_ERROR,err,error,*999)
                END SELECT
                SELECT CASE(EQUATIONS_SET_SUBTYPE)
                CASE(EQUATIONS_SET_RATE_BASED_SMOOTH_MODEL_SUBTYPE, &
                  & EQUATIONS_SET_COMPRESSIBLE_RATE_BASED_SMOOTH_MODEL_SUBTYPE)
                  CALL Field_NumberOfComponentsCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE, &
                    & 1+NUMBER_OF_VOIGT(numberOfDimensions),err,error,*999)
                  DO componentIdx=1,1+NUMBER_OF_VOIGT(numberOfDimensions)
                    CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE,componentIdx, &
                      & FIELD_GAUSS_POINT_BASED_INTERPOLATION,err,error,*999)
                  ENDDO !componentIdx
                CASE(EQUATIONS_SET_RATE_BASED_GROWTH_MODEL_SUBTYPE, &
                  & EQUATIONS_SET_COMPRESSIBLE_RATE_BASED_GROWTH_MODEL_SUBTYPE)
                  CALL Field_NumberOfComponentsCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE, &
                    & 1+NUMBER_OF_VOIGT(numberOfDimensions),err,error,*999)
                  CALL Field_NumberOfComponentsCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U2_VARIABLE_TYPE, &
                    & numberOfDimensions*numberOfDimensions,err,error,*999)
                  DO componentIdx=1,1+NUMBER_OF_VOIGT(numberOfDimensions)
                    CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE,componentIdx, &
                      & FIELD_GAUSS_POINT_BASED_INTERPOLATION,err,error,*999)
                  ENDDO !componentIdx
                  DO componentIdx=1,numberOfDimensions*numberOfDimensions
                    CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U2_VARIABLE_TYPE,componentIdx, &
                      & FIELD_GAUSS_POINT_BASED_INTERPOLATION,err,error,*999)
                  ENDDO !componentIdx
                CASE DEFAULT
                  localError="The equations set subtype of "//TRIM(NumberToVString(EQUATIONS_SET_SUBTYPE,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
                !Check that the pressure values set type is created here?? (second variable is a DELUDELN type, as checked above)
                !\todo: Decide whether these set_types (previous one as well) is to be created by user or automatically..
                IF(.NOT.ASSOCIATED(EQUATIONS_SET_SETUP%FIELD%VARIABLES(2)%PARAMETER_SETS% &
                  & SET_TYPE(FIELD_PRESSURE_VALUES_SET_TYPE)%ptr)) THEN
                    LOCAL_ERROR="Variable 2 of type "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%FIELD%VARIABLES(2)% &
                      & VARIABLE_TYPE,"*",err,error))//" does not have a pressure values set type associated."
                ENDIF
              ENDIF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                CALL Field_CreateFinish(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,err,error,*999)
                SELECT CASE(EQUATIONS_SET_SUBTYPE)
                CASE(EQUATIONS_SET_RATE_BASED_SMOOTH_MODEL_SUBTYPE, &
                  & EQUATIONS_SET_COMPRESSIBLE_RATE_BASED_SMOOTH_MODEL_SUBTYPE)
                  !Initialise U1 variables
                  !Initialise Be prime to the identity matrix
                  DO dimensionIdx=1,numberOfDimensions
                    componentIdx=1+TENSOR_TO_VOIGT(dimensionIdx,dimensionIdx,numberOfDimensions)
                    CALL Field_ComponentValuesInitialise(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,componentIdx,1.0_DP,err,error,*999)
                  ENDDO !dimensionIdx
                CASE(EQUATIONS_SET_RATE_BASED_GROWTH_MODEL_SUBTYPE, &
                  & EQUATIONS_SET_COMPRESSIBLE_RATE_BASED_GROWTH_MODEL_SUBTYPE)
                  !Initialise U1 variables
                  !Initialise Je to 1.0
                  CALL Field_ComponentValuesInitialise(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,1,1.0_DP,err,error,*999)
                  !Initialise Be prime to the identity matrix
                  DO dimensionIdx=1,numberOfDimensions
                    componentIdx=1+TENSOR_TO_VOIGT(dimensionIdx,dimensionIdx,numberOfDimensions)
                    CALL Field_ComponentValuesInitialise(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,componentIdx,1.0_DP,err,error,*999)
                  ENDDO !dimensionIdx
                  !Initialise U2 variables
                  !Initialise S to the identity matrix
                  DO dimensionIdx=1,numberOfDimensions
                    componentIdx=dimensionIdx+(dimensionIdx-1)*numberOfDimensions
                    CALL Field_ComponentValuesInitialise(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U2_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,componentIdx,1.0_DP,err,error,*999)
                  ENDDO !dimensionIdx
                CASE DEFAULT
                  localError="The equations set subtype of "//TRIM(NumberToVString(EQUATIONS_SET_SUBTYPE,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              ENDIF
              CALL Field_ParameterSetEnsureCreated(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_PREVIOUS_VALUES_SET_TYPE,err,error,*999)
              CALL Field_ParameterSetEnsureCreated(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,err,error,*999)
              CALL Field_ParameterSetEnsureCreated(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_PREVIOUS_VALUES_SET_TYPE,err,error,*999)
              CALL Field_ParameterSetEnsureCreated(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_PRESSURE_VALUES_SET_TYPE,err,error,*999)
              SELECT CASE(EQUATIONS_SET_SUBTYPE)
              CASE(EQUATIONS_SET_RATE_BASED_SMOOTH_MODEL_SUBTYPE, &
                & EQUATIONS_SET_COMPRESSIBLE_RATE_BASED_SMOOTH_MODEL_SUBTYPE)
                CALL Field_ParameterSetEnsureCreated(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                  & FIELD_NEXT_VALUES_SET_TYPE,err,error,*999)
              CASE(EQUATIONS_SET_RATE_BASED_GROWTH_MODEL_SUBTYPE, &
                & EQUATIONS_SET_COMPRESSIBLE_RATE_BASED_GROWTH_MODEL_SUBTYPE)
                CALL Field_ParameterSetEnsureCreated(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                  & FIELD_NEXT_VALUES_SET_TYPE,err,error,*999)
                CALL Field_ParameterSetEnsureCreated(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U2_VARIABLE_TYPE, &
                  & FIELD_NEXT_VALUES_SET_TYPE,err,error,*999)
              CASE DEFAULT
                localError="The equations set subtype of "//TRIM(NumberToVString(EQUATIONS_SET_SUBTYPE,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
                & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a finite elasticity equation"
              CALL FlagError(LOCAL_ERROR,err,error,*999)
            END SELECT


          CASE(EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
              & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE)
            !-------------------------------------------------------------------------------
            ! Shared Dependent field setup for multi-physics: elasticity coupled with Darcy
            !-------------------------------------------------------------------------------
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                !Create the auto created dependent field
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET%DEPENDENT% &
                  & DEPENDENT_FIELD,err,error,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,err,error,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & err,error,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY% &
                  & GEOMETRIC_FIELD,err,error,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,4,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE,FIELD_DELVDELN_VARIABLE_TYPE],err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELVDELN_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELVDELN_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & numberOfDimensions,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_COMPONENTS,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & NUMBER_OF_COMPONENTS,err,error,*999)

                SELECT CASE(EQUATIONS_SET_SUBTYPE)
                CASE(EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE)
                  NUMBER_OF_DARCY_COMPONENTS=numberOfDimensions+2 !for INRIA model: velocity components, pressure, mass increase
                CASE (EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE)
                  NUMBER_OF_DARCY_COMPONENTS=numberOfDimensions+1 !for standard Darcy: velocity components and pressure
                CASE (EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE)
                  NUMBER_OF_DARCY_COMPONENTS=numberOfDimensions+1 !for Darcy with pressure driven by solid: velocity components and mass increase
                END SELECT

                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & NUMBER_OF_DARCY_COMPONENTS,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELVDELN_VARIABLE_TYPE, &
                  & NUMBER_OF_DARCY_COMPONENTS,err,error,*999)

                !Elasticity: Default to the geometric interpolation setup
                DO component_idx=1,numberOfDimensions
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                ENDDO !component_idx

                IF (IS_HYDROSTATIC_PRESSURE_DEPENDENT_FIELD) THEN
!kmith :09.0.06.09 - Do we need this ?
                  !Set the hydrostatic component to that of the first geometric component
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_COMPONENTS,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                    & NUMBER_OF_COMPONENTS,GEOMETRIC_MESH_COMPONENT,err,error,*999)
!kmith
                ENDIF

                !Darcy: Default to the geometric interpolation setup
                DO component_idx=1,numberOfDimensions
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELVDELN_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                ENDDO !component_idx

                !Darcy: Default pressure and, if present, mass increase to the first geometric component
                DO component_idx=numberOfDimensions+1,NUMBER_OF_DARCY_COMPONENTS
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELVDELN_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                ENDDO !component_idx

              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  !Elasticity: Set the displacement components to node based interpolation
                  DO component_idx=1,numberOfDimensions
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & component_idx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_DELUDELN_VARIABLE_TYPE,component_idx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  ENDDO !component_idx

                  IF (EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE) THEN
                    !Elasticity: Set the hydrostatic pressure component to node based interpolation
                    !as this is used as the pressure field for the Darcy equations
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & NUMBER_OF_COMPONENTS,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_DELUDELN_VARIABLE_TYPE,NUMBER_OF_COMPONENTS,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  ELSE IF (IS_HYDROSTATIC_PRESSURE_DEPENDENT_FIELD) THEN
                    !Elasticity: Set the hydrostatic pressure component to element based interpolation
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & NUMBER_OF_COMPONENTS,FIELD_ELEMENT_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_DELUDELN_VARIABLE_TYPE,NUMBER_OF_COMPONENTS,FIELD_ELEMENT_BASED_INTERPOLATION,err,error,*999)
                  ENDIF

                  !Darcy: Set the velocity, pressure and, if present, mass increase components to node based interpolation
                  DO component_idx=1,NUMBER_OF_DARCY_COMPONENTS
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                      & component_idx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_DELVDELN_VARIABLE_TYPE,component_idx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  ENDDO !component_idx

                  !Default the scaling to the geometric field scaling
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(LOCAL_ERROR,err,error,*999)
                END SELECT
              ELSE
                !Check the user specified field
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,4,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE,&
                  & FIELD_V_VARIABLE_TYPE,FIELD_DELVDELN_VARIABLE_TYPE],err,error,*999)

                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELVDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)

                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELVDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & numberOfDimensions,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS, &
                  & err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,NUMBER_OF_COMPONENTS, &
                  & err,error,*999)

                SELECT CASE(EQUATIONS_SET_SUBTYPE)
                CASE(EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE)
                  NUMBER_OF_DARCY_COMPONENTS=numberOfDimensions+2 !for INRIA model: velocity components, pressure, mass increase
                CASE (EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE)
                  NUMBER_OF_DARCY_COMPONENTS=numberOfDimensions+1 !for standard Darcy: velocity components and pressure
                CASE (EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE)
                  NUMBER_OF_DARCY_COMPONENTS=numberOfDimensions+1 !for Darcy with pressure driven by solid: velocity components and mass increase
                END SELECT

                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,NUMBER_OF_DARCY_COMPONENTS, &
                  & err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELVDELN_VARIABLE_TYPE, &
                  & NUMBER_OF_DARCY_COMPONENTS,err,error,*999)

                !Check that the impermeability flag values set type is created here??
                !\todo: Decide whether these set_types is to be created by user or automatically..
                IF(.not.ASSOCIATED(EQUATIONS_SET_SETUP%FIELD%VARIABLES(4)%PARAMETER_SETS% &
                  & SET_TYPE(FIELD_IMPERMEABLE_FLAG_VALUES_SET_TYPE)%ptr)) THEN
                    LOCAL_ERROR="Variable 4 of type "//TRIM(NumberToVString(EQUATIONS_SET_SETUP% &
                      & FIELD%VARIABLES(4)% &
                      & VARIABLE_TYPE,"*",err,error))//" does not have an impermeable flag values set type associated."
                ENDIF

                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  !Elasticity:
                  DO component_idx=1,numberOfDimensions
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,component_idx, &
                      & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,component_idx, &
                      & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  ENDDO !component_idx
                  IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE) THEN
                    !If solid hydrostatic pressure is driving Darcy flow, check that pressure uses node based interpolation
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,4, &
                      & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,4, &
                      & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  ENDIF
                  !Darcy:
                  DO component_idx=1,NUMBER_OF_DARCY_COMPONENTS
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,component_idx, &
                      & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELVDELN_VARIABLE_TYPE,component_idx, &
                      & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  ENDDO !component_idx

                CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(LOCAL_ERROR,err,error,*999)
                END SELECT
              ENDIF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                CALL FIELD_CREATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,err,error,*999)
              ENDIF
              CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                 & FIELD_INITIAL_VALUES_SET_TYPE,err,error,*999)
              CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                 & FIELD_RELATIVE_VELOCITY_SET_TYPE,err,error,*999)
              CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                 & FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,err,error,*999)

              CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELVDELN_VARIABLE_TYPE, &
                 & FIELD_IMPERMEABLE_FLAG_VALUES_SET_TYPE,err,error,*999)
            CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
                & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a finite elasticity equation"
              CALL FlagError(LOCAL_ERROR,err,error,*999)
            END SELECT
          !---------------------------------------------------------------------------------------------
          ! Shared Dependent field setup for multi-physics: elasticity coupled with Darcy fluid pressure
          !---------------------------------------------------------------------------------------------
          CASE(EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_STATIC_INRIA_SUBTYPE, &
              & EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_HOLMES_MOW_SUBTYPE, &
              & EQUATIONS_SET_ELASTICITY_FLUID_PRES_HOLMES_MOW_ACTIVE_SUBTYPE)
            NUMBER_OF_DARCY_COMPONENTS=1 !Only solving for the fluid pressure at the moment
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                !Create the auto created dependent field
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET%DEPENDENT% &
                  & DEPENDENT_FIELD,err,error,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_GEOMETRIC_GENERAL_TYPE,err,error,*999)
                CALL FIELD_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,"Dependent Field",err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,err,error,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & err,error,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY% &
                  & GEOMETRIC_FIELD,err,error,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,4,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,[FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE,FIELD_DELVDELN_VARIABLE_TYPE],err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELVDELN_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELVDELN_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & numberOfDimensions,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_COMPONENTS,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & NUMBER_OF_COMPONENTS,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & NUMBER_OF_DARCY_COMPONENTS,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELVDELN_VARIABLE_TYPE, &
                  & NUMBER_OF_DARCY_COMPONENTS,err,error,*999)

                !Set labels
                CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,"U",err,error,*999)
                CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,"del U/del n", &
                  & err,error,*999)
                CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE,"V",err,error,*999)
                CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELVDELN_VARIABLE_TYPE,"del V/del n", &
                  & err,error,*999)
                CALL FIELD_COMPONENT_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,1,"x1",err,error,*999)
                CALL FIELD_COMPONENT_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,2,"x2",err,error,*999)
                CALL FIELD_COMPONENT_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,3,"x3",err,error,*999)
                CALL FIELD_COMPONENT_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                  & "del x1/del n",err,error,*999)
                CALL FIELD_COMPONENT_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,2, &
                  & "del x2/del n",err,error,*999)
                CALL FIELD_COMPONENT_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,3, &
                  & "del x3/del n",err,error,*999)
                CALL FIELD_COMPONENT_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,1,"p",err,error,*999)
                CALL FIELD_COMPONENT_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                  & "del p/del n",err,error,*999)

                !Elasticity: Default to the geometric interpolation setup
                DO component_idx=1,numberOfDimensions
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                ENDDO !component_idx
                !Darcy: Default pressure and mass increase to the first geometric component
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                DO component_idx=1,NUMBER_OF_DARCY_COMPONENTS
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELVDELN_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                ENDDO !component_idx

                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  !Elasticity: Set the displacement components to node based interpolation
                  DO component_idx=1,numberOfDimensions
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & component_idx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_DELUDELN_VARIABLE_TYPE,component_idx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  ENDDO !component_idx
                  !Darcy: Set the pressure and mass increase components to node based interpolation
                  DO component_idx=1,NUMBER_OF_DARCY_COMPONENTS
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                      & component_idx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_DELVDELN_VARIABLE_TYPE,component_idx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  ENDDO !component_idx

                  !Default the scaling to the geometric field scaling
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(LOCAL_ERROR,err,error,*999)
                END SELECT
              ELSE
                !Check the user specified field
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GEOMETRIC_GENERAL_TYPE,err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,4,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE, &
                      & FIELD_V_VARIABLE_TYPE,FIELD_DELVDELN_VARIABLE_TYPE] &
                    & ,err,error,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELVDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)

                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELVDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & numberOfDimensions,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_COMPONENTS,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                    & NUMBER_OF_COMPONENTS,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE, &
                    & NUMBER_OF_DARCY_COMPONENTS,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELVDELN_VARIABLE_TYPE, &
                    & NUMBER_OF_DARCY_COMPONENTS,err,error,*999)

                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  !Elasticity:
                  DO component_idx=1,numberOfDimensions
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,component_idx, &
                      & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,component_idx, &
                      & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  ENDDO !component_idx
                  !Darcy:
                  DO component_idx=1,NUMBER_OF_DARCY_COMPONENTS
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,component_idx, &
                      & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELVDELN_VARIABLE_TYPE,component_idx, &
                      & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  ENDDO !component_idx

                CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(LOCAL_ERROR,err,error,*999)
                END SELECT
              ENDIF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                CALL FIELD_CREATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,err,error,*999)
              ENDIF
            CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
                & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a finite elasticity equation"
              CALL FlagError(LOCAL_ERROR,err,error,*999)
            END SELECT
          !-------------------------------------------------------------------------------
          ! Shared Dependent field setup for multi-physics: elasticity coupled with multi-compartment Darcy
          !-------------------------------------------------------------------------------
          CASE(EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                !Create the auto created dependent field
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET%DEPENDENT% &
                  & DEPENDENT_FIELD,err,error,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_GEOMETRIC_GENERAL_TYPE,err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,err,error,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & err,error,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY% &
                  & GEOMETRIC_FIELD,err,error,*999)
                !Get the number of Darcy compartments from the equations set field
                  EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
                  CALL Field_ParameterSetDataGet(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & FIELD_VALUES_SET_TYPE,EQUATIONS_SET_FIELD_DATA,err,error,*999)
                  Ncompartments=EQUATIONS_SET_FIELD_DATA(2)
                !Set number of variables to be 2+2*Ncompartments
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,(2+2*Ncompartments), &
                   & err,error,*999)
                ALLOCATE(VARIABLE_TYPES(2*Ncompartments+2))
                DO num_var=1,Ncompartments+1
                  VARIABLE_TYPES(2*num_var-1)=FIELD_U_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(num_var-1))
                  VARIABLE_TYPES(2*num_var)=FIELD_DELUDELN_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(num_var-1))
                ENDDO
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET_SETUP%FIELD,VARIABLE_TYPES,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & numberOfDimensions,err,error,*999)
                NUMBER_OF_COMPONENTS=numberOfDimensions+1
                NUMBER_OF_DARCY_COMPONENTS=numberOfDimensions+1 !for Darcy with pressure driven by solid: velocity components and mass increase

                DO num_var=1,2*Ncompartments+2
                  CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,VARIABLE_TYPES(num_var), &
                    & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                  CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,VARIABLE_TYPES(num_var), &
                    & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,VARIABLE_TYPES(num_var), &
                  & NUMBER_OF_COMPONENTS,err,error,*999)
                ENDDO

!                 CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
!                   & NUMBER_OF_COMPONENTS,err,error,*999)
!                 CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
!                   & NUMBER_OF_COMPONENTS,err,error,*999)
!                   NUMBER_OF_DARCY_COMPONENTS=numberOfDimensions+1 !for Darcy with pressure driven by solid: velocity components and mass increase
!                 CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
!                   & NUMBER_OF_DARCY_COMPONENTS,err,error,*999)
!                 CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELVDELN_VARIABLE_TYPE, &
!                   & NUMBER_OF_DARCY_COMPONENTS,err,error,*999)

                !Elasticity: Default to the geometric interpolation setup
                DO component_idx=1,numberOfDimensions
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                ENDDO !component_idx

                !Set the hydrostatic component to that of the first geometric component
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_COMPONENTS,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & NUMBER_OF_COMPONENTS,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                DO num_var=3,2*Ncompartments+2
                  !Darcy: Default to the geometric interpolation setup
                  DO component_idx=1,numberOfDimensions
                    CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, & 
                      & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,VARIABLE_TYPES(num_var), &
                      & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  ENDDO !component_idx
                  !Darcy: Default pressure and, if present, mass increase to the first geometric component
                  DO component_idx=numberOfDimensions+1,NUMBER_OF_DARCY_COMPONENTS
                    CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & 1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,VARIABLE_TYPES(num_var), &
                      & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  ENDDO !component_idx
                ENDDO
                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  !Elasticity: Set the displacement components to node based interpolation
                  DO component_idx=1,numberOfDimensions
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & component_idx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_DELUDELN_VARIABLE_TYPE,component_idx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  ENDDO !component_idx

!                   IF (EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE) THEN
                    !Elasticity: Set the hydrostatic pressure component to node based interpolation
                    !as this is used as the pressure field for the Darcy equations
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_COMPONENTS,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_DELUDELN_VARIABLE_TYPE,NUMBER_OF_COMPONENTS,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
!                   ELSE IF (IS_HYDROSTATIC_PRESSURE_DEPENDENT_FIELD) THEN
!                     !Elasticity: Set the hydrostatic pressure component to element based interpolation
!                     CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
!                       & NUMBER_OF_COMPONENTS,FIELD_ELEMENT_BASED_INTERPOLATION,err,error,*999)
!                     CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
!                       & FIELD_DELUDELN_VARIABLE_TYPE,NUMBER_OF_COMPONENTS,FIELD_ELEMENT_BASED_INTERPOLATION,err,error,*999)
!                   ENDIF
                  DO num_var=3,2*Ncompartments+2
                    !Darcy: Set the velocity, pressure and, if present, mass increase components to node based interpolation
                    DO component_idx=1,NUMBER_OF_DARCY_COMPONENTS
                      CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                        & VARIABLE_TYPES(num_var),component_idx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    ENDDO !component_idx
                  ENDDO
                  !Default the scaling to the geometric field scaling
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(LOCAL_ERROR,err,error,*999)
                END SELECT
              ELSE
                !Check the user specified field
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GEOMETRIC_GENERAL_TYPE,err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
                !Get the number of Darcy compartments from the equations set field
                  EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
                  CALL Field_ParameterSetDataGet(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & FIELD_VALUES_SET_TYPE,EQUATIONS_SET_FIELD_DATA,err,error,*999)
                  Ncompartments=EQUATIONS_SET_FIELD_DATA(2)
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,(2+2*Ncompartments),err,error,*999)
                ALLOCATE(VARIABLE_TYPES(2*Ncompartments+2))
                DO num_var=1,Ncompartments+1
                  VARIABLE_TYPES(2*num_var-1)=FIELD_U_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(num_var-1))
                  VARIABLE_TYPES(2*num_var)=FIELD_DELUDELN_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(num_var-1))
                ENDDO
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,VARIABLE_TYPES,err,error,*999)

                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & numberOfDimensions,err,error,*999)
                NUMBER_OF_COMPONENTS=numberOfDimensions+1
                NUMBER_OF_DARCY_COMPONENTS=numberOfDimensions+1

                DO num_var=1,2*Ncompartments+2
                  CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,VARIABLE_TYPES(num_var),FIELD_VECTOR_DIMENSION_TYPE, &
                    & err,error,*999)
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,VARIABLE_TYPES(num_var),FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,VARIABLE_TYPES(num_var),NUMBER_OF_COMPONENTS, &
                    & err,error,*999)

                ENDDO

                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  !Elasticity:
                 DO component_idx=1,numberOfDimensions
                   CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,component_idx, &
                     & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                   CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,component_idx, &
                     & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                 ENDDO !component_idx
                 !If solid hydrostatic pressure is driving Darcy flow, check that pressure uses node based interpolation
                 CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                 CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                    & NUMBER_OF_COMPONENTS,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)

                DO num_var=3,2*Ncompartments+2
                  !Darcy:
                  DO component_idx=1,NUMBER_OF_DARCY_COMPONENTS
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,VARIABLE_TYPES(num_var),component_idx, &
                      & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  ENDDO !component_idx
                ENDDO
                CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(LOCAL_ERROR,err,error,*999)
                END SELECT
                DEALLOCATE(VARIABLE_TYPES)
              ENDIF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                CALL FIELD_CREATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,err,error,*999)
              ENDIF
              EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
              CALL Field_ParameterSetDataGet(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,EQUATIONS_SET_FIELD_DATA,err,error,*999)
              Ncompartments=EQUATIONS_SET_FIELD_DATA(2)
              ALLOCATE(VARIABLE_TYPES(2*Ncompartments+2))
              DO num_var=1,Ncompartments+1
                VARIABLE_TYPES(2*num_var-1)=FIELD_U_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(num_var-1))
                VARIABLE_TYPES(2*num_var)=FIELD_DELUDELN_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(num_var-1))
              ENDDO
              DO num_var=3,2*Ncompartments+2
                CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,VARIABLE_TYPES(num_var), &
                  & FIELD_INITIAL_VALUES_SET_TYPE,err,error,*999)
                CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,VARIABLE_TYPES(num_var), &
                  & FIELD_RELATIVE_VELOCITY_SET_TYPE,err,error,*999)
                CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,VARIABLE_TYPES(num_var), &
                  & FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,err,error,*999)
              ENDDO
              DEALLOCATE(VARIABLE_TYPES)
            CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
                & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a finite elasticity equation"
              CALL FlagError(LOCAL_ERROR,err,error,*999)
            END SELECT
            !end: Dependent field setup for elasticity coupled with Darcy
          CASE DEFAULT
            LOCAL_ERROR="The equation set subtype of "//TRIM(NumberToVString(EQUATIONS_SET_SUBTYPE,"*",err,error))// &
              & " is invalid for a finite elasticity equation"
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT


        !-----------------------------------------------------------------
        ! I n d e p e n d e n t   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_INDEPENDENT_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)

           SELECT CASE(EQUATIONS_SET_SUBTYPE)
           ! ACTIVE CONTRACTION
           CASE(EQUATIONS_SET_ACTIVECONTRACTION_SUBTYPE)
             NUMBER_OF_COMPONENTS = 8  ! Q1 Q2 Q3 lambda    prev Q1 Q2 Q3 lambda
             IF(EQUATIONS_SET%SOLUTION_METHOD /= EQUATIONS_SET_FEM_SOLUTION_METHOD .OR. &
               &.NOT. EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
               CALL FlagError("Not implemented.",err,error,*999)
             END IF
             CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET%INDEPENDENT% &
               & INDEPENDENT_FIELD,err,error,*999)
             CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_GENERAL_TYPE,err,error,*999)

              CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_INDEPENDENT_TYPE, &
                & err,error,*999)
              CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,err,error,*999)
              CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,GEOMETRIC_DECOMPOSITION, &
                & err,error,*999)

             CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY% &
               & GEOMETRIC_FIELD,err,error,*999)
             CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
               & NUMBER_OF_COMPONENTS,err,error,*999)

             DO component_idx=1,NUMBER_OF_COMPONENTS ! other gauss pt based
               CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                 & FIELD_U_VARIABLE_TYPE,component_idx,FIELD_GAUSS_POINT_BASED_INTERPOLATION,err,error,*999)
             END DO

           !Mooney Rivlin, St Venant Kirchoff and Compressible active contraction subtype
           CASE(EQUATIONS_SET_MOONEY_RIVLIN_ACTIVECONTRACTION_SUBTYPE,EQUATIONS_SET_STVENANT_KIRCHOFF_ACTIVECONTRACTION_SUBTYPE, &
             & EQUATIONS_SET_COMPRESSIBLE_ACTIVECONTRACTION_SUBTYPE, EQUATIONS_SET_HOLZAPFEL_OGDEN_ACTIVECONTRACTION_SUBTYPE, &
             & EQUATIONS_SET_GUCCIONE_ACTIVECONTRACTION_SUBTYPE)
             NUMBER_OF_COMPONENTS = 3 !one contractile stress value for each of the three directions
             IF(EQUATIONS_SET%SOLUTION_METHOD /= EQUATIONS_SET_FEM_SOLUTION_METHOD .OR. &
               &.NOT. EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
               CALL FlagError("Not implemented.",err,error,*999)
             END IF
             CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET%INDEPENDENT% &
               & INDEPENDENT_FIELD,err,error,*999)
             CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_GENERAL_TYPE,err,error,*999)

             CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_INDEPENDENT_TYPE, &
               & err,error,*999)
             CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,err,error,*999)
             CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,GEOMETRIC_DECOMPOSITION, &
               & err,error,*999)

             CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY% &
               & GEOMETRIC_FIELD,err,error,*999)
             CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
               & NUMBER_OF_COMPONENTS,err,error,*999)
             
             !Set component to be gauss point based 
             DO component_idx=1,NUMBER_OF_COMPONENTS 
               CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                 & FIELD_U_VARIABLE_TYPE,component_idx,FIELD_GAUSS_POINT_BASED_INTERPOLATION,err,error,*999)
             ENDDO


           ! COUPLED DARCY
           CASE(EQUATIONS_SET_STANDARD_ELASTICITY_DARCY_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE, &
              & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE,&
              & EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
              & EQUATIONS_SET_ELASTICITY_MULTI_COMPARTMENT_DARCY_INRIA_SUBTYPE, &
              & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_MR_SUBTYPE,EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE)
             IF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
               !Create the auto created dependent field
               CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET%INDEPENDENT% &
                 & INDEPENDENT_FIELD,err,error,*999)
               CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_GENERAL_TYPE,err,error,*999)
               CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_INDEPENDENT_TYPE, &
                 & err,error,*999)
               CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,err,error,*999)
               CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,GEOMETRIC_DECOMPOSITION, &
                 & err,error,*999)
               CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY% &
                 & GEOMETRIC_FIELD,err,error,*999)
               CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,2,err,error,*999)
               CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                 & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
               CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                 & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
               CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                 & FIELD_DP_TYPE,err,error,*999)
               CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                 & FIELD_DP_TYPE,err,error,*999)
               CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                 & numberOfDimensions,err,error,*999)
               NUMBER_OF_COMPONENTS=numberOfDimensions !+1 !Include hydrostatic pressure component
               CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                 & NUMBER_OF_COMPONENTS,err,error,*999)
               CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                 & FIELD_DELUDELN_VARIABLE_TYPE,NUMBER_OF_COMPONENTS,err,error,*999)
               !Default to the geometric interpolation setup
               DO component_idx=1,numberOfDimensions
                 CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                 CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                 CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                   & FIELD_DELUDELN_VARIABLE_TYPE,component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
               ENDDO !component_idx

! !kmith :09.06.09 - Do we need this ?
!               !Set the hydrostatic component to that of the first geometric component
!               CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
!                 & 1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
!               CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
!                 & NUMBER_OF_COMPONENTS,GEOMETRIC_MESH_COMPONENT,err,error,*999)
!               CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
!                 & NUMBER_OF_COMPONENTS,GEOMETRIC_MESH_COMPONENT,err,error,*999)
! !kmith

               SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
               CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                 !Set the displacement components to node based interpolation
                 DO component_idx=1,numberOfDimensions
                   CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                     & FIELD_U_VARIABLE_TYPE,component_idx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                   CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                     & FIELD_DELUDELN_VARIABLE_TYPE,component_idx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                 ENDDO !component_idx
!                 !Set the hydrostatic pressure component to element based interpolation
!                 CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
!                   & FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS,FIELD_ELEMENT_BASED_INTERPOLATION,err,error,*999)
!                 CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
!                   & FIELD_DELUDELN_VARIABLE_TYPE,NUMBER_OF_COMPONENTS,FIELD_ELEMENT_BASED_INTERPOLATION,err,error,*999)
                 !Default the scaling to the geometric field scaling
                 CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                 CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
               CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE DEFAULT
                 LOCAL_ERROR="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                   & " is invalid."
                 CALL FlagError(LOCAL_ERROR,err,error,*999)
               END SELECT
             ELSE !INDEPENDENT_FIELD_AUTO_CREATED
               !Check the user specified field
               CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
               CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
               !Question:Better to leave it up for the user to play around?
               CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,2,err,error,*999)
               CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE], &
                 & err,error,*999)
               CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,ERR, &
                 & ERROR,*999)
               CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                 & err,error,*999)
               CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
               CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
               CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                 & numberOfDimensions,err,error,*999)
               CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS, &
                 & err,error,*999)
               CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,NUMBER_OF_COMPONENTS,&
                 & err,error,*999)
               SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
               CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                 !DO component_idx=1,numberOfDimensions
                 !  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,component_idx, &
                 !    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                 !  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,component_idx, &
                 !    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                 !ENDDO !component_idx
               CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE DEFAULT
                 LOCAL_ERROR="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                   & " is invalid."
                 CALL FlagError(LOCAL_ERROR,err,error,*999)
               END SELECT
             ENDIF !INDEPENDENT_FIELD_AUTO_CREATED

           ! BIOELECTRICS COUPLED TO FINITE ELASTICITY
           CASE(EQUATIONS_SET_STANDARD_MONODOMAIN_ELASTICITY_SUBTYPE)
             IF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
               !Create the auto created dependent field
               CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET%INDEPENDENT% &
                 & INDEPENDENT_FIELD,err,error,*999)
               CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_GENERAL_TYPE,err,error,*999)
               CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_INDEPENDENT_TYPE, &
                 & err,error,*999)
               CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,err,error,*999)
               CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,GEOMETRIC_DECOMPOSITION, &
                 & err,error,*999)
               CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY% &
                 & GEOMETRIC_FIELD,err,error,*999)
               CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,1,err,error,*999)
               CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                 & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
               CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                 & FIELD_DP_TYPE,err,error,*999)
               CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                 & 1,err,error,*999)
               !Default to the first component of the geometric interpolation setup
               CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & 1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
               CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                 & 1,GEOMETRIC_MESH_COMPONENT,err,error,*999)

               SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
               CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                 !Set to node based interpolation
                 CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                   & FIELD_U_VARIABLE_TYPE,1,FIELD_GAUSS_POINT_BASED_INTERPOLATION,err,error,*999)
                 !Default the scaling to the geometric field scaling
                 CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                 CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
               CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE DEFAULT
                 LOCAL_ERROR="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                   & " is invalid."
                 CALL FlagError(LOCAL_ERROR,err,error,*999)
               END SELECT
             ELSE !INDEPENDENT_FIELD_AUTO_CREATED
               !Check the user specified field
               CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
               CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
               !Question:Better to leave it up for the user to play around?
               CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,1,err,error,*999)
               CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],err,error,*999)
               CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE,ERR, &
                 & ERROR,*999)
               CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
               CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                 & err,error,*999)
               SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
               CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                 !do/check nothing???
               CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE DEFAULT
                 LOCAL_ERROR="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                   & " is invalid."
                 CALL FlagError(LOCAL_ERROR,err,error,*999)
               END SELECT
             ENDIF !INDEPENDENT_FIELD_AUTO_CREATED

           CASE(EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE,EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, &
             & EQUATIONS_SET_1D3D_MONODOMAIN_ACTIVE_STRAIN_SUBTYPE)
             IF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
               !Create the auto created independent field
               CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET%INDEPENDENT% &
                 & INDEPENDENT_FIELD,err,error,*999)
               CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_GENERAL_TYPE,err,error,*999)
               CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_INDEPENDENT_TYPE, &
                 & err,error,*999)
               CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,err,error,*999)
               CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,GEOMETRIC_DECOMPOSITION, &
                 & err,error,*999)
               CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY% &
                 & GEOMETRIC_FIELD,err,error,*999)
               CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,2,err,error,*999)
               CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,[FIELD_U_VARIABLE_TYPE, &
                 & FIELD_V_VARIABLE_TYPE],err,error,*999)
               IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE) THEN
                 CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
               ENDIF
               CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                 & FIELD_DP_TYPE,err,error,*999)
               CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                 & FIELD_INTG_TYPE,err,error,*999)
               IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE) THEN
                 CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & 1,err,error,*999)
               ELSEIF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE) THEN
                 CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & 6,err,error,*999)
               ELSEIF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_1D3D_MONODOMAIN_ACTIVE_STRAIN_SUBTYPE) THEN
                 CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & 5,err,error,*999)
               ENDIF
               CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                 & numberOfDimensions+1,err,error,*999)
               !Default to the first component of the geometric interpolation setup
               CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                 & 1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
               CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                 & 1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
               IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE) THEN
                 CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & 2,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                 CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & 3,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                 CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & 4,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                 CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & 5,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                 CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & 6,GEOMETRIC_MESH_COMPONENT,err,error,*999)
               ELSEIF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_1D3D_MONODOMAIN_ACTIVE_STRAIN_SUBTYPE) THEN
                 CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & 2,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                 CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & 3,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                 CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & 4,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                 CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & 5,GEOMETRIC_MESH_COMPONENT,err,error,*999)
               ENDIF
               SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
               CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                 !Set to node based interpolation
                 CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                   & FIELD_U_VARIABLE_TYPE,1,FIELD_GAUSS_POINT_BASED_INTERPOLATION,err,error,*999)
                 IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE) THEN
                   CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                     & FIELD_U_VARIABLE_TYPE,2,FIELD_GAUSS_POINT_BASED_INTERPOLATION,err,error,*999)
                   CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                     & FIELD_U_VARIABLE_TYPE,3,FIELD_GAUSS_POINT_BASED_INTERPOLATION,err,error,*999)
                   CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                     & FIELD_U_VARIABLE_TYPE,4,FIELD_GAUSS_POINT_BASED_INTERPOLATION,err,error,*999)
                   CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                     & FIELD_U_VARIABLE_TYPE,5,FIELD_GAUSS_POINT_BASED_INTERPOLATION,err,error,*999)
                   CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                     & FIELD_U_VARIABLE_TYPE,6,FIELD_GAUSS_POINT_BASED_INTERPOLATION,err,error,*999)
                 ENDIF
                 DO component_idx=1,numberOfDimensions
                   CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                     & FIELD_V_VARIABLE_TYPE,component_idx,FIELD_ELEMENT_BASED_INTERPOLATION,err,error,*999)
                 ENDDO
                 !Default the scaling to the geometric field scaling
                 CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                 CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
               CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE DEFAULT
                 LOCAL_ERROR="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                   & " is invalid."
                 CALL FlagError(LOCAL_ERROR,err,error,*999)
               END SELECT
             ELSE !INDEPENDENT_FIELD_AUTO_CREATED
               !Check the user specified field
               CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
               CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
               !Question:Better to leave it up for the user to play around?
               CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,2,err,error,*999)
               CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE],ERR, &
                 & ERROR,*999)
               IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE) THEN
                 CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE,ERR, &
                   & ERROR,*999)
               ENDIF
               CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
               CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_INTG_TYPE,err,error,*999)
               IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE) THEN
                 CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                   & err,error,*999)
               ELSE IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE) THEN
                 CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,6, &
                   & err,error,*999)
               ELSEIF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_1D3D_MONODOMAIN_ACTIVE_STRAIN_SUBTYPE) THEN
                 CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,5, &
                   & err,error,*999)
               ENDIF
               CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,numberOfDimensions+1, &
                 & err,error,*999)
               SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
               CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                 !do/check nothing???
               CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE DEFAULT
                 LOCAL_ERROR="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                   & " is invalid."
                 CALL FlagError(LOCAL_ERROR,err,error,*999)
               END SELECT
             ENDIF !INDEPENDENT_FIELD_AUTO_CREATED

           CASE(EQUATIONS_SET_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE)
             IF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
               !Create the auto created independent field
               CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET%INDEPENDENT% &
                 & INDEPENDENT_FIELD,err,error,*999)
               CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_GENERAL_TYPE,err,error,*999)
               CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_INDEPENDENT_TYPE, &
                 & err,error,*999)
               CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,err,error,*999)
               CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,GEOMETRIC_DECOMPOSITION, &
                 & err,error,*999)
               CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY% &
                 & GEOMETRIC_FIELD,err,error,*999)
               CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,3,err,error,*999)
               CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,[FIELD_U_VARIABLE_TYPE, &
                 & FIELD_V_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE],err,error,*999)
               CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                 & FIELD_DP_TYPE,err,error,*999)
               CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                 & FIELD_INTG_TYPE,err,error,*999)
               CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                 & FIELD_DP_TYPE,err,error,*999)
               CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                 & 1,err,error,*999)
               CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                 & numberOfDimensions+1,err,error,*999)
               CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                 & 3,err,error,*999)
               !Default to the first component of the geometric interpolation setup
               CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                 & 1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
               CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                 & 1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
               CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                 & 2,GEOMETRIC_MESH_COMPONENT,err,error,*999)
               CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                 & 1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
               SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
               CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                 !Set to node based interpolation
                 CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                   & FIELD_U_VARIABLE_TYPE,1,FIELD_GAUSS_POINT_BASED_INTERPOLATION,err,error,*999)
                 CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                   & FIELD_U1_VARIABLE_TYPE,1,FIELD_GAUSS_POINT_BASED_INTERPOLATION,err,error,*999)
                 DO component_idx=1,numberOfDimensions
                   CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                     & FIELD_V_VARIABLE_TYPE,component_idx,FIELD_ELEMENT_BASED_INTERPOLATION,err,error,*999)
                 ENDDO
                 !Default the scaling to the geometric field scaling
                 CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                 CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
               CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE DEFAULT
                 LOCAL_ERROR="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                   & " is invalid."
                 CALL FlagError(LOCAL_ERROR,err,error,*999)
               END SELECT
             ELSE !INDEPENDENT_FIELD_AUTO_CREATED
               !Check the user specified field
               CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
               CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
               !Question:Better to leave it up for the user to play around?
               CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,3,err,error,*999)
               CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE, &
                 & FIELD_U1_VARIABLE_TYPE],err,error,*999)
               CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
               CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_INTG_TYPE,err,error,*999)
               CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
               CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                 & err,error,*999)
               CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,numberOfDimensions+1, &
                 & err,error,*999)
               CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE,3, &
                 & err,error,*999)
               SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
               CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                 !do/check nothing???
               CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE DEFAULT
                 LOCAL_ERROR="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                   & " is invalid."
                 CALL FlagError(LOCAL_ERROR,err,error,*999)
               END SELECT
             ENDIF !INDEPENDENT_FIELD_AUTO_CREATED

           CASE(EQUATIONS_SET_MULTISCALE_ACTIVE_STRAIN_SUBTYPE)
             IF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
               !Create the auto created independent field
               CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET%INDEPENDENT% &
                 & INDEPENDENT_FIELD,err,error,*999)
               CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_GENERAL_TYPE,err,error,*999)
               CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_INDEPENDENT_TYPE, &
                 & err,error,*999)
               CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,err,error,*999)
               CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,GEOMETRIC_DECOMPOSITION, &
                 & err,error,*999)
               CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY% &
                 & GEOMETRIC_FIELD,err,error,*999)
               CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,1,err,error,*999)
               CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,[FIELD_U_VARIABLE_TYPE],ERR, &
                 & ERROR,*999)
               CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                 & FIELD_DP_TYPE,err,error,*999)
               CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                 & 4,err,error,*999)
               !Default to the first component of the geometric interpolation setup
               CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                 & 1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
               CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                 & 1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
               CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                 & 2,GEOMETRIC_MESH_COMPONENT,err,error,*999)
               CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                 & 3,GEOMETRIC_MESH_COMPONENT,err,error,*999)
               CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                 & 4,GEOMETRIC_MESH_COMPONENT,err,error,*999)
               SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
               CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                 !Set to node based interpolation
                 CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                   & FIELD_U_VARIABLE_TYPE,1,FIELD_GAUSS_POINT_BASED_INTERPOLATION,err,error,*999)
                 CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                   & FIELD_U_VARIABLE_TYPE,2,FIELD_GAUSS_POINT_BASED_INTERPOLATION,err,error,*999)
                 !Default the scaling to the geometric field scaling
                 CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                 CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
               CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE DEFAULT
                 LOCAL_ERROR="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                   & " is invalid."
                 CALL FlagError(LOCAL_ERROR,err,error,*999)
               END SELECT
             ELSE !INDEPENDENT_FIELD_AUTO_CREATED
               !Check the user specified field
               CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
               CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
               !Question:Better to leave it up for the user to play around?
               CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,1,err,error,*999)
               CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],ERR, &
                 & ERROR,*999)
               CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)

               CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,4, &
                 & err,error,*999)

               SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
               CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                 !do/check nothing???
               CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                 CALL FlagError("Not implemented.",err,error,*999)
               CASE DEFAULT
                 LOCAL_ERROR="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                   & " is invalid."
                 CALL FlagError(LOCAL_ERROR,err,error,*999)
               END SELECT
             ENDIF !INDEPENDENT_FIELD_AUTO_CREATED

           CASE DEFAULT
             LOCAL_ERROR="The third equations set specification of "// &
               & TRIM(NumberToVString(EQUATIONS_SET_SUBTYPE,"*",err,error))// &
               & " is invalid for an independent field of a finite elasticity equation."
             CALL FlagError(LOCAL_ERROR,err,error,*999)
           END SELECT
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
              CALL FIELD_CREATE_FINISH(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,err,error,*999)
              ! initialize values for active contraction independent field. TODO: actual init for z, trpn, or flag to presolve
              IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_ACTIVECONTRACTION_SUBTYPE) THEN
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_COMPONENTS,err,error,*999)
                DO component_idx=1,NUMBER_OF_COMPONENTS
                  CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & FIELD_VALUES_SET_TYPE,component_idx,0.0_DP,Err,ERROR,*999)
                ENDDO
              ENDIF
            ENDIF
            IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE) THEN
              CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_PREVIOUS_VALUES_SET_TYPE,err,error,*999)
              CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                & FIELD_PREVIOUS_VALUES_SET_TYPE,err,error,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a finite elasticity equation"
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT

        !-----------------------------------------------------------------
        ! M a t e r i a l s   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
            IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
              NUMBER_OF_FLUID_COMPONENTS=0
              SELECT CASE(EQUATIONS_SET_SUBTYPE)
              CASE(EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE,EQUATIONS_SET_NO_SUBTYPE, &
                & EQUATIONS_SET_MOONEY_RIVLIN_ACTIVECONTRACTION_SUBTYPE, &
                & EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
                & EQUATIONS_SET_STANDARD_MONODOMAIN_ELASTICITY_SUBTYPE, &
                & EQUATIONS_SET_MR_AND_GROWTH_LAW_IN_CELLML_SUBTYPE)
                NUMBER_OF_COMPONENTS=2;
              CASE(EQUATIONS_SET_ACTIVE_STRAIN_SUBTYPE)
                NUMBER_OF_COMPONENTS=8;
              CASE(EQUATIONS_SET_MULTISCALE_ACTIVE_STRAIN_SUBTYPE)
                NUMBER_OF_COMPONENTS=8;
              CASE(EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
                & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
                NUMBER_OF_COMPONENTS=3;
              CASE(EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE,EQUATIONS_SET_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE)
                NUMBER_OF_COMPONENTS=5;
              CASE(EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE)
                NUMBER_OF_COMPONENTS=6;
              CASE(EQUATIONS_SET_1D3D_MONODOMAIN_ACTIVE_STRAIN_SUBTYPE)
                NUMBER_OF_COMPONENTS=8;
              CASE(EQUATIONS_SET_MEMBRANE_SUBTYPE)
                !\todo Currently the number of components for a membrane problem's material field has been set to 3 in 3D space or
                ! 2 in 2D space to work with a Mooney Rivlin material (2 material parameters) and a membrane thickness parameter
                ! (only if in 3D space). Extra subtypes will need to be added to use other constitutive relations with
                ! membrane mechanics problems.
                IF (numberOfDimensions==3) THEN
                  NUMBER_OF_COMPONENTS=3
                ELSE
                  NUMBER_OF_COMPONENTS=2
                ENDIF
              CASE(EQUATIONS_SET_STVENANT_KIRCHOFF_ACTIVECONTRACTION_SUBTYPE)
                NUMBER_OF_COMPONENTS=2;
              CASE(EQUATIONS_SET_ISOTROPIC_EXPONENTIAL_SUBTYPE)
                NUMBER_OF_COMPONENTS=2;
              CASE(EQUATIONS_SET_TRANSVERSE_ISOTROPIC_EXPONENTIAL_SUBTYPE)
                NUMBER_OF_COMPONENTS=5;
              CASE(EQUATIONS_SET_TRANSVERSE_ISOTROPIC_POLYNOMIAL_SUBTYPE)
                NUMBER_OF_COMPONENTS=4;
              CASE(EQUATIONS_SET_TRANSVERSE_ISOTROPIC_ACTIVE_SUBTYPE)
                NUMBER_OF_COMPONENTS=5;
              CASE(EQUATIONS_SET_ANISOTROPIC_POLYNOMIAL_SUBTYPE)
                NUMBER_OF_COMPONENTS=8;
              CASE(EQUATIONS_SET_ANISOTROPIC_POLYNOMIAL_ACTIVE_SUBTYPE)
                NUMBER_OF_COMPONENTS=12;
              CASE(EQUATIONS_SET_TRANS_ISOTROPIC_ACTIVE_TRANSITION_SUBTYPE)
                NUMBER_OF_COMPONENTS=11;
              CASE(EQUATIONS_SET_ORTHOTROPIC_MATERIAL_COSTA_SUBTYPE,EQUATIONS_SET_ACTIVECONTRACTION_SUBTYPE)
                NUMBER_OF_COMPONENTS=7;
              CASE(EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE,EQUATIONS_SET_COMPRESSIBLE_ACTIVECONTRACTION_SUBTYPE,&
                & EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE)
                NUMBER_OF_COMPONENTS=3;
              CASE(EQUATIONS_SET_ORTHOTROPIC_MATERIAL_HOLZAPFEL_OGDEN_SUBTYPE, &
                & EQUATIONS_SET_HOLZAPFEL_OGDEN_ACTIVECONTRACTION_SUBTYPE)
                NUMBER_OF_COMPONENTS=8;
              CASE(EQUATIONS_SET_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE)
                NUMBER_OF_COMPONENTS=2;  
              CASE(EQUATIONS_SET_NEARLY_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE)
                NUMBER_OF_COMPONENTS=3;
              CASE(EQUATIONS_SET_TRANSVERSE_ISOTROPIC_GUCCIONE_SUBTYPE, &
                & EQUATIONS_SET_GUCCIONE_ACTIVECONTRACTION_SUBTYPE, &
                & EQUATIONS_SET_REFERENCE_STATE_TRANSVERSE_GUCCIONE_SUBTYPE, &
                & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_HUMPHREY_YIN_SUBTYPE)
                NUMBER_OF_COMPONENTS=4;
              CASE(EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_STATIC_INRIA_SUBTYPE)
                NUMBER_OF_COMPONENTS=6;
                NUMBER_OF_FLUID_COMPONENTS=8
              CASE(EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_HOLMES_MOW_SUBTYPE)
                NUMBER_OF_COMPONENTS=4
                NUMBER_OF_FLUID_COMPONENTS=8
              CASE(EQUATIONS_SET_ELASTICITY_FLUID_PRES_HOLMES_MOW_ACTIVE_SUBTYPE)
                NUMBER_OF_COMPONENTS=6
                NUMBER_OF_FLUID_COMPONENTS=8
              CASE(EQUATIONS_SET_GROWTH_LAW_IN_CELLML_SUBTYPE)
                NUMBER_OF_COMPONENTS=24
              CASE(EQUATIONS_SET_CONSTITUTIVE_LAW_IN_CELLML_EVALUATE_SUBTYPE)
                CALL FlagError("Materials field is not required for CellML based constituative laws.",err,error,*999)
              CASE(EQUATIONS_SET_CONSTITUTIVE_AND_GROWTH_LAW_IN_CELLML_SUBTYPE)
                CALL FlagError("Materials field is not required for CellML based constituative laws.",err,error,*999)
              CASE(EQUATIONS_SET_RATE_BASED_SMOOTH_MODEL_SUBTYPE)
                NUMBER_OF_COMPONENTS=7
              CASE(EQUATIONS_SET_COMPRESSIBLE_RATE_BASED_SMOOTH_MODEL_SUBTYPE)
                NUMBER_OF_COMPONENTS=8
              CASE(EQUATIONS_SET_RATE_BASED_GROWTH_MODEL_SUBTYPE)
                NUMBER_OF_COMPONENTS=9+3+NUMBER_OF_VOIGT(numberOfDimensions)
              CASE(EQUATIONS_SET_COMPRESSIBLE_RATE_BASED_GROWTH_MODEL_SUBTYPE)
                NUMBER_OF_COMPONENTS=9+3+NUMBER_OF_VOIGT(numberOfDimensions)
              CASE DEFAULT
                LOCAL_ERROR="The third equations set specification of "// &
                  & TRIM(NumberToVString(EQUATIONS_SET_SUBTYPE,"*",err,error))// &
                  & " is not valid for a finite elasticity type of an elasticity equation set."
                CALL FlagError(LOCAL_ERROR,err,error,*999)
              END SELECT
              IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                !Create the auto created materials field
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_MATERIALS% &
                  & MATERIALS_FIELD,err,error,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_MATERIAL_TYPE,err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,err,error,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & err,error,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,EQUATIONS_SET%GEOMETRY% &
                  & GEOMETRIC_FIELD,err,error,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_MESH_COMPONENT,err,error,*999)  ! get 1 = x (?) component

                !U variable type is constitutive law parameters
                !V variable type has one component, density
                IF(NUMBER_OF_FLUID_COMPONENTS>0) THEN
                  !If coupled with Darcy pressure equation then a shared material field is used and Darcy material parameters are in U1
                  CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,3,err,error,*999)
                  CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,[FIELD_U_VARIABLE_TYPE, &
                      & FIELD_V_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE],err,error,*999)
                ELSE
                  CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,2,err,error,*999)
                  CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,[FIELD_U_VARIABLE_TYPE, &
                      & FIELD_V_VARIABLE_TYPE],err,error,*999)
                ENDIF
                CALL FIELD_LABEL_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,"Materials",err,error,*999)

                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & NUMBER_OF_COMPONENTS,err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE,"Parameters",err,error,*999)

                IF(EQUATIONS_SET_SUBTYPE == EQUATIONS_SET_ACTIVECONTRACTION_SUBTYPE) THEN
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & 1,err,error,*999) ! just 1 component: activation time
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD, &
                   & FIELD_V_VARIABLE_TYPE,1 ,FIELD_GAUSS_POINT_BASED_INTERPOLATION,err,error,*999) ! gauss pt based interp.
                ELSE
                  !Solid density
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                     & 1,err,error,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & 1,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & 1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE,"Density",err,error,*999)
                ENDIF

                DO component_idx=1,NUMBER_OF_COMPONENTS
                !Default the materials components to the geometric interpolation setup with constant interpolation
                  CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                ENDDO

                IF(NUMBER_OF_FLUID_COMPONENTS>0) THEN
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U1_VARIABLE_TYPE, &
                     & NUMBER_OF_FLUID_COMPONENTS,err,error,*999)
                  CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U1_VARIABLE_TYPE,"Fluid Parameters", &
                    & err,error,*999)
                ENDIF
                DO component_idx=1,NUMBER_OF_FLUID_COMPONENTS
                  CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U1_VARIABLE_TYPE, &
                    & component_idx,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U1_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                ENDDO

                !Default the field scaling to that of the geometric field
                CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                CALL FIELD_SCALING_TYPE_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
              ELSE
                !Check the user specified field
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_MATERIAL_TYPE,err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_GET(EQUATIONS_SET_SETUP%FIELD,EQUATIONS_SET_FIELD_NUMBER_OF_VARIABLES,err,error,*999)
                SELECT CASE(EQUATIONS_SET_FIELD_NUMBER_OF_VARIABLES)
                CASE(1)
                  CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                CASE(2)
                  CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE, &
                      & FIELD_V_VARIABLE_TYPE],err,error,*999)
                CASE(3)
                  CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE, &
                      & FIELD_V_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE],err,error,*999)
                CASE DEFAULT
                  LOCAL_ERROR="Invalid number of variables. The number of variables for field number "// &
                    & TRIM(NumberToVString(EQUATIONS_SET_SETUP%FIELD%USER_NUMBER,"*",err,error))//" is "// &
                    & TRIM(NumberToVString(EQUATIONS_SET_SETUP%FIELD%NUMBER_OF_VARIABLES,"*",err,error))// &
                    & " but should be either 1, 2 or 3"
                  CALL FlagError(LOCAL_ERROR,err,error,*999)
                END SELECT
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_COMPONENTS,err,error,*999)
                IF (EQUATIONS_SET_FIELD_NUMBER_OF_VARIABLES>1) THEN
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE, &
                      & 1,err,error,*999)
                ENDIF
                IF (EQUATIONS_SET_FIELD_NUMBER_OF_VARIABLES>2) THEN
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE, &
                      & NUMBER_OF_FLUID_COMPONENTS,err,error,*999)
                ENDIF
              ENDIF
            ELSE
              CALL FlagError("Equations set materials is not associated.",err,error,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
            IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
              IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                !Finish creating the materials field
                CALL FIELD_CREATE_FINISH(EQUATIONS_MATERIALS%MATERIALS_FIELD,err,error,*999)
                !Set the default values for the materials field
                !Don't bother checking equations types, just default to all componets = 1.0
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_COMPONENTS,err,error,*999)
                DO component_idx=1,NUMBER_OF_COMPONENTS
                  CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,component_idx,1.0_DP,err,error,*999)
                ENDDO
                !Initialise density to 0
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,1,0.0_DP,err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Equations set materials is not associated.",err,error,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a finite elasticity equation."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
          IF(ASSOCIATED(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD)) THEN
            CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
              & numberOfDimensions,err,error,*999)
            NUMBER_OF_COMPONENTS=numberOfDimensions
          ELSE
            CALL FlagError("Equations set geometrc field is not associated",err,error,*999)
          ENDIF
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%SOURCE%SOURCE_FIELD_AUTO_CREATED) THEN
              CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET%SOURCE% &
                & SOURCE_FIELD,err,error,*999)
              CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_GENERAL_TYPE,err,error,*999)
              CALL FIELD_LABEL_SET(EQUATIONS_SET%SOURCE%SOURCE_FIELD,"Source Field",err,error,*999)
              CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_INDEPENDENT_TYPE, &
                & err,error,*999)
              CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,err,error,*999)
              CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD,GEOMETRIC_DECOMPOSITION, &
                & err,error,*999)
              CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD,EQUATIONS_SET%GEOMETRY% &
                & GEOMETRIC_FIELD,err,error,*999)

              CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD,1,err,error,*999)
              CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD,[FIELD_U_VARIABLE_TYPE],err,error,*999)
              CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                & NUMBER_OF_COMPONENTS,err,error,*999)

              CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE,"Gravity",err,error,*999)
              CALL FIELD_COMPONENT_LABEL_SET(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE,1,"g1",err,error,*999)
              CALL FIELD_COMPONENT_LABEL_SET(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE,2,"g2",err,error,*999)
              IF(NUMBER_OF_COMPONENTS==3) THEN
                CALL FIELD_COMPONENT_LABEL_SET(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE,3,"g3",err,error,*999)
              ENDIF

              DO component_idx=1,NUMBER_OF_COMPONENTS
                CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD, &
                  & FIELD_U_VARIABLE_TYPE,component_idx,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
              END DO
            ELSE
              !Check the user specified field
              CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
              CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
              CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,1,err,error,*999)
              CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],err,error,*999)
              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                & err,error,*999)
              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_COMPONENTS,err,error,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(ASSOCIATED(EQUATIONS_SET%SOURCE)) THEN
              IF(EQUATIONS_SET%SOURCE%SOURCE_FIELD_AUTO_CREATED) THEN
                !Finish creating the source field
                CALL FIELD_CREATE_FINISH(EQUATIONS_SET%SOURCE%SOURCE_FIELD,err,error,*999)
                !Set the default values for the field
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_COMPONENTS,err,error,*999)
                DO component_idx=1,NUMBER_OF_COMPONENTS-1
                  CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,component_idx,0.0_DP,err,error,*999)
                ENDDO
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,NUMBER_OF_COMPONENTS,-9.80665_DP,err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Equations set source is not associated.",err,error,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a finite elasticity equation."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
              DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
              IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN
                  SELECT CASE(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE)
                  CASE(EQUATIONS_SET_FINITE_ELASTICITY_CYLINDER)
                    IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE) THEN
                      !Create analytic field if required
                      !Set analtyic function type
                      EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_FINITE_ELASTICITY_CYLINDER
                    ELSE
                      LOCAL_ERROR="The thrid equations set specification of "// &
                        & TRIM(NumberToVString(EQUATIONS_SET_SUBTYPE,"*",err,error))// &
                        & " is invalid. The analytic function type of "// &
                        & TRIM(NumberToVString(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
                        & " requires that the third equations set specification be a Mooney-Rivlin finite elasticity equation."
                      CALL FlagError(LOCAL_ERROR,err,error,*999)
                    ENDIF
                  CASE DEFAULT
                    LOCAL_ERROR="The specified analytic function type of "// &
                      & TRIM(NumberToVString(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
                      & " is invalid for a finite elasticity equation."
                    CALL FlagError(LOCAL_ERROR,err,error,*999)
                  END SELECT
                ELSE
                  CALL FlagError("Equations set geometric field is not associated.",err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("Equations set dependent field is not associated.",err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Equations set dependent field has not been finished.",err,error,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
              ANALYTIC_FIELD=>EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD
              IF(ASSOCIATED(ANALYTIC_FIELD)) THEN
                IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD_AUTO_CREATED) THEN
                  CALL FIELD_CREATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,err,error,*999)
                ENDIF
              ENDIF
            ELSE
              CALL FlagError("Equations set analytic is not associated.",err,error,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a finite elasticity equation."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
              !Start the equations creation
              CALL Equations_CreateStart(EQUATIONS_SET,EQUATIONS,err,error,*999)
              CALL Equations_LinearityTypeSet(EQUATIONS,EQUATIONS_NONLINEAR,err,error,*999)
              SELECT CASE(EQUATIONS_SET_SUBTYPE)
              CASE(EQUATIONS_SET_ACTIVECONTRACTION_SUBTYPE, & 
                & EQUATIONS_SET_RATE_BASED_SMOOTH_MODEL_SUBTYPE,EQUATIONS_SET_COMPRESSIBLE_RATE_BASED_SMOOTH_MODEL_SUBTYPE, & 
                & EQUATIONS_SET_RATE_BASED_GROWTH_MODEL_SUBTYPE,EQUATIONS_SET_COMPRESSIBLE_RATE_BASED_GROWTH_MODEL_SUBTYPE)
                CALL Equations_TimeDependenceTypeSet(EQUATIONS,EQUATIONS_QUASISTATIC,err,error,*999)
              CASE DEFAULT
                CALL Equations_TimeDependenceTypeSet(EQUATIONS,EQUATIONS_STATIC,err,error,*999)
              END SELECT
            ELSE
              CALL FlagError("Equations set dependent field has not been finished.",err,error,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              !Finish the equations creation
              CALL EquationsSet_EquationsGet(EQUATIONS_SET,EQUATIONS,err,error,*999)
              CALL Equations_CreateFinish(equations,err,error,*999)
              NULLIFY(vectorEquations)
              CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
              !Create the equations mapping.              
              CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_DELUDELN_VARIABLE_TYPE,vectorMapping,err,error,*999)
              SELECT CASE(EQUATIONS_SET_SUBTYPE)
              CASE(EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_STATIC_INRIA_SUBTYPE, &
                & EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_HOLMES_MOW_SUBTYPE, &
                & EQUATIONS_SET_ELASTICITY_FLUID_PRES_HOLMES_MOW_ACTIVE_SUBTYPE)
                !Residual vector also depends on the fluid pressure variable
                CALL EquationsMapping_ResidualVariablesNumberSet(vectorMapping,2,err,error,*999)
                CALL EquationsMapping_ResidualVariableTypesSet(vectorMapping, &
                    & [FIELD_U_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE],err,error,*999)
              CASE DEFAULT
                !Single residual variable
                CALL EquationsMapping_ResidualVariableTypesSet(vectorMapping,[FIELD_U_VARIABLE_TYPE],err,error,*999)
              END SELECT
              CALL EquationsMapping_LinearMatricesNumberSet(vectorMapping,0,err,error,*999)
              CALL EquationsMapping_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE,err,error,*999)
              CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
              !Create the equations matrices
              CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
              ! set structure and storage types
              SELECT CASE(EQUATIONS%sparsityType)
                CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                  CALL EquationsMatrices_NonlinearStorageTypeSet(vectorMatrices,MATRIX_BLOCK_STORAGE_TYPE, &
                    & err,error,*999)
                CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                  CALL EquationsMatrices_NonlinearStorageTypeSet(vectorMatrices, & 
                    & MATRIX_COMPRESSED_ROW_STORAGE_TYPE,err,error,*999)
                  CALL EquationsMatrices_NonlinearStructureTypeSet(vectorMatrices, & 
                    & EQUATIONS_MATRIX_FEM_STRUCTURE,err,error,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The equations matrices sparsity type of "// &
                    & TRIM(NumberToVString(EQUATIONS%sparsityType,"*",err,error))//" is invalid."
                  CALL FlagError(LOCAL_ERROR,err,error,*999)
              END SELECT
              !Set Jacobian matrices calculation type to default finite difference. 
              CALL EquationsMatrices_JacobianTypesSet(vectorMatrices,[EQUATIONS_JACOBIAN_FINITE_DIFFERENCE_CALCULATED], &
                & err,error,*999)
              CALL EquationsMatrices_VectorCreateFinish(vectorMatrices,err,error,*999)
            CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE DEFAULT
              LOCAL_ERROR="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                & " is invalid."
              CALL FlagError(LOCAL_ERROR,err,error,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a finite elasticity equation."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_DERIVED_TYPE)
          ! We want to be able to set which derived variables are calculated before finishing the derived
          ! field, so don't create field variables or check the provided field until the finish action.
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%derived%derivedFieldAutoCreated) THEN
              CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET%derived% &
                & derivedField,err,error,*999)
              CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%derived%derivedField,FIELD_GENERAL_TYPE,err,error,*999)
              CALL FIELD_LABEL_SET(EQUATIONS_SET%derived%derivedField,"Derived Field",err,error,*999)
              CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%derived%derivedField,FIELD_DEPENDENT_TYPE, &
                & err,error,*999)
              CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,err,error,*999)
              CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%derived%derivedField,GEOMETRIC_DECOMPOSITION, &
                & err,error,*999)
              CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%derived%derivedField,EQUATIONS_SET%GEOMETRY% &
                & GEOMETRIC_FIELD,err,error,*999)
            END IF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(ASSOCIATED(EQUATIONS_SET%derived)) THEN
              ALLOCATE(VARIABLE_TYPES(EQUATIONS_SET%derived%numberOfVariables),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate derived field variable types.",err,error,*999)
              varIdx=0
              DO derivedIdx=1,EQUATIONS_SET_NUMBER_OF_TENSOR_TYPES
                IF(EQUATIONS_SET%derived%variableTypes(derivedIdx)/=0) THEN
                  varIdx=varIdx+1
                  VARIABLE_TYPES(varIdx)=EQUATIONS_SET%derived%variableTypes(derivedIdx)
                END IF
              END DO
              numberOfTensorComponents=NUMBER_OF_VOIGT(numberOfDimensions)
              IF(EQUATIONS_SET%derived%derivedFieldAutoCreated) THEN
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%derived%derivedField, &
                  & EQUATIONS_SET%derived%numberOfVariables,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%derived%derivedField,VARIABLE_TYPES,err,error,*999)
                DO derivedIdx=1,EQUATIONS_SET_NUMBER_OF_TENSOR_TYPES
                  variableType=EQUATIONS_SET%derived%variableTypes(derivedIdx)
                  IF(variableType/=0) THEN
                    CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%derived%derivedField,variableType, &
                      & FIELD_DP_TYPE,err,error,*999)
                    SELECT CASE(derivedidx)
                    CASE(EQUATIONS_SET_DEFORMATION_GRADIENT_TENSOR)
                      CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%derived%derivedField,variableType, &
                        & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                      CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%derived%derivedField,variableType,"Strain",err,error,*999)
                      CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%derived%derivedField,variableType, &
                        & numberOfTensorComponents,err,error,*999)
                    CASE(EQUATIONS_SET_R_CAUCHY_GREEN_DEFORMATION_TENSOR)
                      CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%derived%derivedField,variableType, &
                        & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                      CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%derived%derivedField,variableType,"Strain",err,error,*999)
                      CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%derived%derivedField,variableType, &
                        & numberOfTensorComponents,err,error,*999)
                    CASE(EQUATIONS_SET_L_CAUCHY_GREEN_DEFORMATION_TENSOR)
                      CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%derived%derivedField,variableType, &
                        & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                      CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%derived%derivedField,variableType,"Strain",err,error,*999)
                      CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%derived%derivedField,variableType, &
                        & numberOfTensorComponents,err,error,*999)
                    CASE(EQUATIONS_SET_GREEN_LAGRANGE_STRAIN_TENSOR)
                      CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%derived%derivedField,variableType, &
                        & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                      CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%derived%derivedField,variableType,"Strain",err,error,*999)
                      CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%derived%derivedField,variableType, &
                        & numberOfTensorComponents,err,error,*999)
                    CASE(EQUATIONS_SET_CAUCHY_STRESS_TENSOR)
                      CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%derived%derivedField,variableType, &
                        & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                      CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%derived%derivedField,variableType,"Stress",err,error,*999)
                      CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%derived%derivedField,variableType, &
                        & numberOfTensorComponents,err,error,*999)
                    CASE(EQUATIONS_SET_FIRST_PK_STRESS_TENSOR)
                      CALL FlagError("Not implemented.",err,error,*999)
                    CASE(EQUATIONS_SET_SECOND_PK_STRESS_TENSOR)
                      CALL FlagError("Not implemented.",err,error,*999)
                    CASE DEFAULT
                      CALL FlagError("The specified derived field type of "//TRIM(NumberToVString(derivedIdx,"*",err,error))// &
                        & " is not supported for a finite elasticity equations set type.",err,error,*999)
                    END SELECT
                  END IF
                END DO
                !Finish creating the derived field
                CALL FIELD_CREATE_FINISH(EQUATIONS_SET%derived%derivedField,err,error,*999)
              ELSE
                !Check the user specified derived field
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
 
                DO derivedIdx=1,EQUATIONS_SET_NUMBER_OF_TENSOR_TYPES
                  variableType=EQUATIONS_SET%derived%variableTypes(derivedIdx)
                  IF(variableType/=0) THEN
                    CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                    SELECT CASE(derivedidx)
                    CASE(EQUATIONS_SET_DEFORMATION_GRADIENT_TENSOR)
                      CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET%derived%derivedField,variableType, &
                        & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                      CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET%derived%derivedField,variableType, &
                        & numberOfTensorComponents,err,error,*999)
                    CASE(EQUATIONS_SET_R_CAUCHY_GREEN_DEFORMATION_TENSOR)
                      CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET%derived%derivedField,variableType, &
                        & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                      CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET%derived%derivedField,variableType, &
                        & numberOfTensorComponents,err,error,*999)
                    CASE(EQUATIONS_SET_L_CAUCHY_GREEN_DEFORMATION_TENSOR)
                      CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET%derived%derivedField,variableType, &
                        & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                      CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET%derived%derivedField,variableType, &
                        & numberOfTensorComponents,err,error,*999)
                    CASE(EQUATIONS_SET_GREEN_LAGRANGE_STRAIN_TENSOR)
                      CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET%derived%derivedField,variableType, &
                        & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                      CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET%derived%derivedField,variableType, &
                        & numberOfTensorComponents,err,error,*999)
                    CASE(EQUATIONS_SET_CAUCHY_STRESS_TENSOR)
                      CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET%derived%derivedField,variableType, &
                        & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                      CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET%derived%derivedField,variableType, &
                        & numberOfTensorComponents,err,error,*999)
                    CASE(EQUATIONS_SET_FIRST_PK_STRESS_TENSOR)
                      CALL FlagError("Not implemented.",err,error,*999)
                    CASE(EQUATIONS_SET_SECOND_PK_STRESS_TENSOR)
                      CALL FlagError("Not implemented.",err,error,*999)
                    CASE DEFAULT
                      CALL FlagError("The specified derived field type of "//TRIM(NumberToVString(derivedIdx,"*",err,error))// &
                        & " is not supported for a finite elasticity equations set type.",err,error,*999)
                    END SELECT
                  END IF
                END DO
              END IF
            ELSE
              CALL FlagError("Equations set derived is not associated.",err,error,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a finite elasticity equation."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
            & " is invalid for a finite elasticity equation."
          CALL FlagError(LOCAL_ERROR,err,error,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="The third equations set specification of "//TRIM(NumberToVString(EQUATIONS_SET_SUBTYPE,"*",err,error))// &
          & " is not valid for a finite elasticity type of an elasticity equation set."
        CALL FlagError(LOCAL_ERROR,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

    EXITS("FINITE_ELASTICITY_EQUATIONS_SET_SETUP")
    RETURN
999 ERRORSEXITS("FINITE_ELASTICITY_EQUATIONS_SET_SETUP",err,error)
    RETURN 1
  END SUBROUTINE FINITE_ELASTICITY_EQUATIONS_SET_SETUP

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a finite elasticity equation type of an elasticity equations set class.
  SUBROUTINE FiniteElasticity_EquationsSetSolutionMethodSet(EQUATIONS_SET,SOLUTION_METHOD,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("FiniteElasticity_EquationsSetSolutionMethodSet",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a finite elasticity type equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
      CASE(EQUATIONS_SET_MEMBRANE_SUBTYPE,EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE, &
          & EQUATIONS_SET_MOONEY_RIVLIN_ACTIVECONTRACTION_SUBTYPE, &
          & EQUATIONS_SET_STVENANT_KIRCHOFF_ACTIVECONTRACTION_SUBTYPE, &
          & EQUATIONS_SET_ISOTROPIC_EXPONENTIAL_SUBTYPE,EQUATIONS_SET_TRANSVERSE_ISOTROPIC_ACTIVE_SUBTYPE, &
          & EQUATIONS_SET_TRANS_ISOTROPIC_ACTIVE_TRANSITION_SUBTYPE, &
          & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_EXPONENTIAL_SUBTYPE,EQUATIONS_SET_TRANSVERSE_ISOTROPIC_POLYNOMIAL_SUBTYPE, &
          & EQUATIONS_SET_ANISOTROPIC_POLYNOMIAL_SUBTYPE,EQUATIONS_SET_ANISOTROPIC_POLYNOMIAL_ACTIVE_SUBTYPE, &
          & EQUATIONS_SET_ORTHOTROPIC_MATERIAL_COSTA_SUBTYPE,EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE, &
          & EQUATIONS_SET_COMPRESSIBLE_ACTIVECONTRACTION_SUBTYPE, &
          & EQUATIONS_SET_ACTIVECONTRACTION_SUBTYPE,EQUATIONS_SET_NO_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
          & EQUATIONS_SET_ORTHOTROPIC_MATERIAL_HOLZAPFEL_OGDEN_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE,EQUATIONS_SET_NEARLY_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE,EQUATIONS_SET_TRANSVERSE_ISOTROPIC_GUCCIONE_SUBTYPE, &
          & EQUATIONS_SET_REFERENCE_STATE_TRANSVERSE_GUCCIONE_SUBTYPE, &
          & EQUATIONS_SET_CONSTITUTIVE_LAW_IN_CELLML_EVALUATE_SUBTYPE, &
          & EQUATIONS_SET_CONSTITUTIVE_AND_GROWTH_LAW_IN_CELLML_SUBTYPE, &
          & EQUATIONS_SET_GROWTH_LAW_IN_CELLML_SUBTYPE, &
          & EQUATIONS_SET_MR_AND_GROWTH_LAW_IN_CELLML_SUBTYPE, &
          & EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_STATIC_INRIA_SUBTYPE, &
          & EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_HOLMES_MOW_SUBTYPE, &
          & EQUATIONS_SET_ELASTICITY_FLUID_PRES_HOLMES_MOW_ACTIVE_SUBTYPE, &
          & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_HUMPHREY_YIN_SUBTYPE, & 
          & EQUATIONS_SET_STANDARD_MONODOMAIN_ELASTICITY_SUBTYPE, &
          & EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE, &
          & EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, &
          & EQUATIONS_SET_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE, &
          & EQUATIONS_SET_HOLZAPFEL_OGDEN_ACTIVECONTRACTION_SUBTYPE, &
          & EQUATIONS_SET_ACTIVE_STRAIN_SUBTYPE, &
          & EQUATIONS_SET_MULTISCALE_ACTIVE_STRAIN_SUBTYPE, &
          & EQUATIONS_SET_1D3D_MONODOMAIN_ACTIVE_STRAIN_SUBTYPE, &
          & EQUATIONS_SET_GUCCIONE_ACTIVECONTRACTION_SUBTYPE, &
          & EQUATIONS_SET_RATE_BASED_SMOOTH_MODEL_SUBTYPE,EQUATIONS_SET_COMPRESSIBLE_RATE_BASED_SMOOTH_MODEL_SUBTYPE, &
          & EQUATIONS_SET_RATE_BASED_GROWTH_MODEL_SUBTYPE,EQUATIONS_SET_COMPRESSIBLE_RATE_BASED_GROWTH_MODEL_SUBTYPE)
        SELECT CASE(SOLUTION_METHOD)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          EQUATIONS_SET%SOLUTION_METHOD=EQUATIONS_SET_FEM_SOLUTION_METHOD
        CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE DEFAULT
          LOCAL_ERROR="The specified solution method of "//TRIM(NumberToVString(SOLUTION_METHOD,"*",err,error))//" is invalid."
          CALL FlagError(LOCAL_ERROR,err,error,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="Equations set subtype of "//TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
          & " is not valid for a finite elasticity equation type of an elasticity equations set class."
        CALL FlagError(LOCAL_ERROR,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

    EXITS("FiniteElasticity_EquationsSetSolutionMethodSet")
    RETURN
999 ERRORS("FiniteElasticity_EquationsSetSolutionMethodSet",err,error)
    EXITS("FiniteElasticity_EquationsSetSolutionMethodSet")
    RETURN 1

  END SUBROUTINE FiniteElasticity_EquationsSetSolutionMethodSet

  !
  !================================================================================================================================
  !

  !>Sets the equation specification for a finite elasticity equation type of an elasticity equations set class.
  SUBROUTINE FiniteElasticity_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to set the equation subtype for
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations set specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: subtype

    ENTERS("FiniteElasticity_EquationsSetSpecificationSet",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(SIZE(specification,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a finite elasticity type equations set.", &
          & err,error,*999)
      END IF
      subtype=specification(3)
      SELECT CASE(subtype)
      CASE(EQUATIONS_SET_MEMBRANE_SUBTYPE,EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE, &
          & EQUATIONS_SET_MOONEY_RIVLIN_ACTIVECONTRACTION_SUBTYPE, &
          & EQUATIONS_SET_STVENANT_KIRCHOFF_ACTIVECONTRACTION_SUBTYPE, &
          & EQUATIONS_SET_ISOTROPIC_EXPONENTIAL_SUBTYPE,EQUATIONS_SET_TRANSVERSE_ISOTROPIC_ACTIVE_SUBTYPE, &
          & EQUATIONS_SET_TRANS_ISOTROPIC_ACTIVE_TRANSITION_SUBTYPE, &
          & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_EXPONENTIAL_SUBTYPE,EQUATIONS_SET_TRANSVERSE_ISOTROPIC_POLYNOMIAL_SUBTYPE, &
          & EQUATIONS_SET_ANISOTROPIC_POLYNOMIAL_SUBTYPE,EQUATIONS_SET_ANISOTROPIC_POLYNOMIAL_ACTIVE_SUBTYPE, &
          & EQUATIONS_SET_ORTHOTROPIC_MATERIAL_COSTA_SUBTYPE,EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE, &
          & EQUATIONS_SET_COMPRESSIBLE_ACTIVECONTRACTION_SUBTYPE, &
          & EQUATIONS_SET_ACTIVECONTRACTION_SUBTYPE, EQUATIONS_SET_NO_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
          & EQUATIONS_SET_ORTHOTROPIC_MATERIAL_HOLZAPFEL_OGDEN_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE,EQUATIONS_SET_NEARLY_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE,EQUATIONS_SET_TRANSVERSE_ISOTROPIC_GUCCIONE_SUBTYPE, &
          & EQUATIONS_SET_GUCCIONE_ACTIVECONTRACTION_SUBTYPE, EQUATIONS_SET_REFERENCE_STATE_TRANSVERSE_GUCCIONE_SUBTYPE, &
          & EQUATIONS_SET_CONSTITUTIVE_LAW_IN_CELLML_EVALUATE_SUBTYPE, &
          & EQUATIONS_SET_CONSTITUTIVE_AND_GROWTH_LAW_IN_CELLML_SUBTYPE, &
          & EQUATIONS_SET_GROWTH_LAW_IN_CELLML_SUBTYPE, &
          & EQUATIONS_SET_MR_AND_GROWTH_LAW_IN_CELLML_SUBTYPE, &
          & EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_STATIC_INRIA_SUBTYPE, &
          & EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_HOLMES_MOW_SUBTYPE, &
          & EQUATIONS_SET_ELASTICITY_FLUID_PRES_HOLMES_MOW_ACTIVE_SUBTYPE, &
          & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_HUMPHREY_YIN_SUBTYPE, &
          & EQUATIONS_SET_STANDARD_MONODOMAIN_ELASTICITY_SUBTYPE,EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE, &
          & EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE,EQUATIONS_SET_ACTIVE_STRAIN_SUBTYPE, &
          & EQUATIONS_SET_MULTISCALE_ACTIVE_STRAIN_SUBTYPE,EQUATIONS_SET_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE, &
          & EQUATIONS_SET_1D3D_MONODOMAIN_ACTIVE_STRAIN_SUBTYPE, &
          & EQUATIONS_SET_HOLZAPFEL_OGDEN_ACTIVECONTRACTION_SUBTYPE, &
          & EQUATIONS_SET_RATE_BASED_SMOOTH_MODEL_SUBTYPE,EQUATIONS_SET_COMPRESSIBLE_RATE_BASED_SMOOTH_MODEL_SUBTYPE, &
          & EQUATIONS_SET_RATE_BASED_GROWTH_MODEL_SUBTYPE,EQUATIONS_SET_COMPRESSIBLE_RATE_BASED_GROWTH_MODEL_SUBTYPE)
        !Set full specification
        IF(ALLOCATED(equationsSet%specification)) THEN
          CALL FlagError("Equations set specification is already allocated.",err,error,*999)
        ELSE
          ALLOCATE(equationsSet%specification(3),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate equations set specification.",err,error,*999)
        END IF
        equationsSet%specification(1:3)=[EQUATIONS_SET_ELASTICITY_CLASS,EQUATIONS_SET_FINITE_ELASTICITY_TYPE,subtype]
      CASE DEFAULT
        localError="Equations set subtype "//TRIM(NumberToVstring(subtype,"*",err,error))// &
          & " is not valid for a finite elasticity equation type of an elasticity equations set class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    END IF

    EXITS("FiniteElasticity_EquationsSetSpecificationSet")
    RETURN
999 ERRORS("FiniteElasticity_EquationsSetSpecificationSet",err,error)
    EXITS("FiniteElasticity_EquationsSetSpecificationSet")
    RETURN 1

  END SUBROUTINE FiniteElasticity_EquationsSetSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the finite elasticity problem.
  SUBROUTINE FINITE_ELASTICITY_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem set to setup a Laplace equation on.
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP,CONTROL_LOOP_ROOT
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_TYPE), POINTER :: CELLML_SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(CELLML_EQUATIONS_TYPE), POINTER :: CELLML_EQUATIONS
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR,localError
    INTEGER(INTG) :: PROBLEM_SUBTYPE

    ENTERS("FINITE_ELASTICITY_PROBLEM_SETUP",err,error,*999)

    NULLIFY(CONTROL_LOOP)
    NULLIFY(SOLVER)
    NULLIFY(CELLML_SOLVER)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVERS)
    NULLIFY(CELLML_EQUATIONS)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(.NOT.ALLOCATED(PROBLEM%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(PROBLEM%SPECIFICATION,1)<3) THEN
        CALL FlagError("Problem specification must have three entries for a finite elasticity problem.",err,error,*999)
      ENDIF
      PROBLEM_SUBTYPE=PROBLEM%SPECIFICATION(3)
      SELECT CASE(PROBLEM_SUBTYPE)
      CASE(PROBLEM_NO_SUBTYPE,PROBLEM_STATIC_FINITE_ELASTICITY_SUBTYPE,PROBLEM_QUASISTATIC_FINITE_ELASTICITY_SUBTYPE, &
        & PROBLEM_QUASISTATIC_FINITE_ELASTICITY_WITH_GROWTH_SUBTYPE,PROBLEM_DYNAMIC_FINITE_ELASTICITY_SUBTYPE, &
        & PROBLEM_MULTISCALE_FINITE_ELASTICITY_SUBTYPE,PROBLEM_FINITE_ELASTICITY_WITH_ACTIVE_SUBTYPE, &
        & PROBLEM_FINITE_ELASTICITY_WITH_CELLML_SUBTYPE,PROBLEM_FINITE_ELASTICITY_WITH_GROWTH_CELLML_SUBTYPE)
        SELECT CASE(PROBLEM_SETUP%SETUP_TYPE)
        CASE(PROBLEM_SETUP_INITIAL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Do nothing????
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Do nothing????
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a finite elasticity problem."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Set up a simple control loop: default is load increment type now
            CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,err,error,*999)
            SELECT CASE(PROBLEM_SUBTYPE)
            CASE(PROBLEM_NO_SUBTYPE,PROBLEM_STATIC_FINITE_ELASTICITY_SUBTYPE,PROBLEM_FINITE_ELASTICITY_WITH_CELLML_SUBTYPE)
              CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE,err,error,*999)
            CASE(PROBLEM_QUASISTATIC_FINITE_ELASTICITY_SUBTYPE,PROBLEM_QUASISTATIC_FINITE_ELASTICITY_WITH_GROWTH_SUBTYPE, &
              & PROBLEM_DYNAMIC_FINITE_ELASTICITY_SUBTYPE,PROBLEM_FINITE_ELASTICITY_WITH_ACTIVE_SUBTYPE, &
              & PROBLEM_FINITE_ELASTICITY_WITH_GROWTH_CELLML_SUBTYPE)
!!TODO: Should we have a load step loop inside the time loop?
              CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,PROBLEM_CONTROL_TIME_LOOP_TYPE,err,error,*999)
              CALL ControlLoop_LabelSet(CONTROL_LOOP,"Time Loop",err,error,*999)
            CASE DEFAULT
              localError="Problem subtype "//TRIM(NumberToVString(PROBLEM_SUBTYPE,"*",err,error))// &
                & " is not valid for a finite elasticity type of an elasticity problem class."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Finish the control loops
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,err,error,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a finite elasticity problem."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVERS_TYPE)
          !Get the control loop
          CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Start the solvers creation
            CALL SOLVERS_CREATE_START(CONTROL_LOOP,SOLVERS,err,error,*999)
            SELECT CASE(PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_NO_SUBTYPE,PROBLEM_STATIC_FINITE_ELASTICITY_SUBTYPE,PROBLEM_DYNAMIC_FINITE_ELASTICITY_SUBTYPE, &
              & PROBLEM_FINITE_ELASTICITY_WITH_ACTIVE_SUBTYPE)
              CALL SOLVERS_NUMBER_SET(SOLVERS,1,err,error,*999)
              !Set the solver to be a nonlinear solver
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              CALL SOLVER_TYPE_SET(SOLVER,SOLVER_NONLINEAR_TYPE,err,error,*999)
              !Set solver defaults
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,err,error,*999)              
            CASE(PROBLEM_QUASISTATIC_FINITE_ELASTICITY_SUBTYPE)
              CALL SOLVERS_NUMBER_SET(SOLVERS,2,err,error,*999)
              !Set the first solver to be an CellML Evaluator for time varying boundary conditions
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              CALL SOLVER_TYPE_SET(SOLVER,SOLVER_CELLML_EVALUATOR_TYPE,err,error,*999)
              CALL Solver_LabelSet(solver,"Boundary Condition Evaluation Solver",err,error,*999)
              !Set solver defaults
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,err,error,*999)
              !Set the second solver to be a nonlinear solver to solve the mechanics
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
              CALL SOLVER_TYPE_SET(SOLVER,SOLVER_NONLINEAR_TYPE,err,error,*999)
              CALL Solver_LabelSet(solver,"Nonlinear Solver",err,error,*999)
              !Set solver defaults
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,err,error,*999)
            CASE(PROBLEM_QUASISTATIC_FINITE_ELASTICITY_WITH_GROWTH_SUBTYPE)
              CALL SOLVERS_NUMBER_SET(SOLVERS,3,err,error,*999)
              !Set the first solver to be an CellML Evaluator for time varying boundary conditions
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              CALL SOLVER_TYPE_SET(SOLVER,SOLVER_CELLML_EVALUATOR_TYPE,err,error,*999)
              CALL Solver_LabelSet(solver,"Boundary Condition Evaluation Solver",err,error,*999)
              !Set solver defaults
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,err,error,*999)
              !Set the second solver to be an ODE integrator for growth
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
              CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DAE_TYPE,err,error,*999)
              !Set solver defaults
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,err,error,*999)
              !Set the third solver to be a nonlinear solver for elasticity
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(SOLVERS,3,SOLVER,err,error,*999)
              CALL SOLVER_TYPE_SET(SOLVER,SOLVER_NONLINEAR_TYPE,err,error,*999)
              !Set solver defaults
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,err,error,*999)
            CASE(PROBLEM_FINITE_ELASTICITY_WITH_GROWTH_CELLML_SUBTYPE) 
              CALL SOLVERS_NUMBER_SET(SOLVERS,2,err,error,*999)
              !Set the first solver to be an ODE integrator for growth
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DAE_TYPE,err,error,*999)
              !Set solver defaults
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,err,error,*999)
              !Set the second solver to be a nonlinear solver for elasticity
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
              CALL SOLVER_TYPE_SET(SOLVER,SOLVER_NONLINEAR_TYPE,err,error,*999)
              !Set solver defaults
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,err,error,*999)
              !Create the CellML evaluator solver for a constituative law via CellML
              NULLIFY(CELLML_SOLVER)
              CALL SOLVER_NEWTON_CELLML_EVALUATOR_CREATE(SOLVER,CELLML_SOLVER,err,error,*999)
              !Link the CellML evaluator solver to the solver
              CALL SOLVER_LINKED_SOLVER_ADD(SOLVER,CELLML_SOLVER,SOLVER_CELLML_EVALUATOR_TYPE,err,error,*999)
            CASE(PROBLEM_FINITE_ELASTICITY_WITH_CELLML_SUBTYPE)
              CALL SOLVERS_NUMBER_SET(SOLVERS,1,err,error,*999)
              !Set the solver to be a nonlinear solver
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              CALL SOLVER_TYPE_SET(SOLVER,SOLVER_NONLINEAR_TYPE,err,error,*999)
              !Set solver defaults
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,err,error,*999)
              !Create the CellML evaluator solver
              NULLIFY(CELLML_SOLVER)
              CALL SOLVER_NEWTON_CELLML_EVALUATOR_CREATE(SOLVER,CELLML_SOLVER,err,error,*999)
              !Link the CellML evaluator solver to the solver
              CALL SOLVER_LINKED_SOLVER_ADD(SOLVER,CELLML_SOLVER,SOLVER_CELLML_EVALUATOR_TYPE,err,error,*999)
            CASE(PROBLEM_MULTISCALE_FINITE_ELASTICITY_SUBTYPE)
              CALL SOLVERS_NUMBER_SET(SOLVERS,2,err,error,*999)
              !Set the first solver to be a DAE solver
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DAE_TYPE,err,error,*999)
              CALL SOLVER_LABEL_SET(SOLVER,"ODE_Solver",err,error,*999)
              !Set solver defaults
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,err,error,*999)
              NULLIFY(SOLVER)
              !Set the second solver to be a nonlinear solver 
              CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
              CALL SOLVER_TYPE_SET(SOLVER,SOLVER_NONLINEAR_TYPE,err,error,*999)
              CALL SOLVER_LABEL_SET(SOLVER,"Nonlinear_Solver",err,error,*999)
              !Set solver defaults
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,err,error,*999)
            CASE DEFAULT
              localError="The third problem specification of "//TRIM(NumberToVString(problem%specification(3),"*",err,error))// &
                & " is not valid for a finite elasticity type of an elasticity problem class."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the solvers
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(SOLVERS,err,error,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a finite elasticity problem."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            !Get the solver
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            SELECT CASE(PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_NO_SUBTYPE,PROBLEM_STATIC_FINITE_ELASTICITY_SUBTYPE,PROBLEM_DYNAMIC_FINITE_ELASTICITY_SUBTYPE, &
              & PROBLEM_FINITE_ELASTICITY_WITH_ACTIVE_SUBTYPE,PROBLEM_FINITE_ELASTICITY_WITH_CELLML_SUBTYPE)
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
            CASE(PROBLEM_QUASISTATIC_FINITE_ELASTICITY_SUBTYPE,PROBLEM_MULTISCALE_FINITE_ELASTICITY_SUBTYPE, &
              & PROBLEM_FINITE_ELASTICITY_WITH_GROWTH_CELLML_SUBTYPE)
              CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
            CASE(PROBLEM_QUASISTATIC_FINITE_ELASTICITY_WITH_GROWTH_SUBTYPE)
              CALL SOLVERS_SOLVER_GET(SOLVERS,3,SOLVER,err,error,*999)
            CASE DEFAULT
              localError="The third problem specification of "//TRIM(NumberToVString(problem%specification(3),"*",err,error))// &
                & " is not valid for a finite elasticity type of an elasticity problem class."
              CALL FlagError(localError,err,error,*999)
            END SELECT
            !Create the solver equations
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,err,error,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
            !Set time dependence
            SELECT CASE(PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_NO_SUBTYPE,PROBLEM_STATIC_FINITE_ELASTICITY_SUBTYPE,PROBLEM_MULTISCALE_FINITE_ELASTICITY_SUBTYPE, &
              & PROBLEM_FINITE_ELASTICITY_WITH_CELLML_SUBTYPE)
              CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,err,error,*999)
            CASE(PROBLEM_QUASISTATIC_FINITE_ELASTICITY_SUBTYPE,PROBLEM_QUASISTATIC_FINITE_ELASTICITY_WITH_GROWTH_SUBTYPE, & 
              & PROBLEM_FINITE_ELASTICITY_WITH_ACTIVE_SUBTYPE,PROBLEM_FINITE_ELASTICITY_WITH_GROWTH_CELLML_SUBTYPE)
              CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_QUASISTATIC,err,error,*999)
            CASE(PROBLEM_DYNAMIC_FINITE_ELASTICITY_SUBTYPE)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE DEFAULT
              localError="The third problem specification of "//TRIM(NumberToVString(problem%specification(3),"*",err,error))// &
                & " is not valid for a finite elasticity type of an elasticity problem class."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            !Get the solver equations
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            SELECT CASE(PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_NO_SUBTYPE,PROBLEM_STATIC_FINITE_ELASTICITY_SUBTYPE,PROBLEM_DYNAMIC_FINITE_ELASTICITY_SUBTYPE, &
              & PROBLEM_FINITE_ELASTICITY_WITH_ACTIVE_SUBTYPE,PROBLEM_FINITE_ELASTICITY_WITH_CELLML_SUBTYPE)
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
            CASE(PROBLEM_QUASISTATIC_FINITE_ELASTICITY_SUBTYPE,PROBLEM_MULTISCALE_FINITE_ELASTICITY_SUBTYPE, &
              & PROBLEM_FINITE_ELASTICITY_WITH_GROWTH_CELLML_SUBTYPE)
              CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
            CASE(PROBLEM_QUASISTATIC_FINITE_ELASTICITY_WITH_GROWTH_SUBTYPE)
              CALL SOLVERS_SOLVER_GET(SOLVERS,3,SOLVER,err,error,*999)
            CASE DEFAULT
              localError="The third problem specification of "//TRIM(NumberToVString(problem%specification(3),"*",err,error))// &
                & " is not valid for a finite elasticity type of an elasticity problem class."
              CALL FlagError(localError,err,error,*999)
            END SELECT
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,err,error,*999)
            !Finish the solver equations creation
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,err,error,*999)
         CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a finite elasticity problem."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CELLML_EQUATIONS_TYPE)
          !Get the control loop
          CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
          CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            SELECT CASE(PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_NO_SUBTYPE,PROBLEM_STATIC_FINITE_ELASTICITY_SUBTYPE,PROBLEM_DYNAMIC_FINITE_ELASTICITY_SUBTYPE, &
              & PROBLEM_FINITE_ELASTICITY_WITH_ACTIVE_SUBTYPE)
              !Do nothing
            CASE(PROBLEM_QUASISTATIC_FINITE_ELASTICITY_SUBTYPE)
              !Get the CellML solver
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,CELLML_SOLVER,err,error,*999)
              !Create the CellML equations
              CALL CELLML_EQUATIONS_CREATE_START(CELLML_SOLVER,CELLML_EQUATIONS,err,error,*999)
              !Set the time dependence
              CALL CellMLEquations_TimeDependenceTypeSet(CELLML_EQUATIONS,CELLML_EQUATIONS_QUASISTATIC,err,error,*999)
            CASE(PROBLEM_QUASISTATIC_FINITE_ELASTICITY_WITH_GROWTH_SUBTYPE)
              !Get the CellML BC solver 
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,CELLML_SOLVER,err,error,*999)
              !Create the CellML equations
              CALL CELLML_EQUATIONS_CREATE_START(CELLML_SOLVER,CELLML_EQUATIONS,err,error,*999)
              !Set the time dependence
              CALL CellMLEquations_TimeDependenceTypeSet(CELLML_EQUATIONS,CELLML_EQUATIONS_QUASISTATIC,err,error,*999)
              !Get the CellML Growth solver
              NULLIFY(CELLML_SOLVER)
              CALL SOLVERS_SOLVER_GET(SOLVERS,2,CELLML_SOLVER,err,error,*999)
              !Create the CellML equations
              NULLIFY(CELLML_EQUATIONS)
              CALL CELLML_EQUATIONS_CREATE_START(CELLML_SOLVER,CELLML_EQUATIONS,err,error,*999)
              !Set the time dependence
              CALL CellMLEquations_TimeDependenceTypeSet(CELLML_EQUATIONS,CELLML_EQUATIONS_QUASISTATIC,err,error,*999)
            CASE(PROBLEM_FINITE_ELASTICITY_WITH_CELLML_SUBTYPE)
              !Get the nonlinear solver
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              !Get the CellML evaluator solver
              CALL SOLVER_NEWTON_CELLML_SOLVER_GET(SOLVER,CELLML_SOLVER,err,error,*999)
              !Create the CellML equations
              CALL CELLML_EQUATIONS_CREATE_START(CELLML_SOLVER,CELLML_EQUATIONS,err,error,*999)
              !Set the time dependence
              CALL CellMLEquations_TimeDependenceTypeSet(CELLML_EQUATIONS,CELLML_EQUATIONS_STATIC,err,error,*999)
            CASE(PROBLEM_FINITE_ELASTICITY_WITH_GROWTH_CELLML_SUBTYPE) 
              !Get the CellML integrator solver
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,CELLML_SOLVER,err,error,*999)
              CALL CELLML_EQUATIONS_CREATE_START(CELLML_SOLVER,CELLML_EQUATIONS,err,error,*999)
              NULLIFY(CELLML_SOLVER)
              NULLIFY(CELLML_EQUATIONS)
              !Get the nonlinear solver
              CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
              !Get the CellML evaluator solver
              CALL SOLVER_NEWTON_CELLML_SOLVER_GET(SOLVER,CELLML_SOLVER,err,error,*999)
              !Create the CellML equations
              CALL CELLML_EQUATIONS_CREATE_START(CELLML_SOLVER,CELLML_EQUATIONS, &
              & err,error,*999)
            CASE(PROBLEM_MULTISCALE_FINITE_ELASTICITY_SUBTYPE)
              !Create the CellML equations for the first DAE solver
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              NULLIFY(CELLML_EQUATIONS)
              CALL CELLML_EQUATIONS_CREATE_START(SOLVER,CELLML_EQUATIONS,err,error,*999)
            CASE DEFAULT
              LOCAL_ERROR="The third problem specification of "//TRIM(NumberToVString(PROBLEM_SUBTYPE,"*",err,error))// &
                & " is not valid for a finite elasticity type of an elasticity problem."
              CALL FlagError(LOCAL_ERROR,err,error,*999)
            END SELECT
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            SELECT CASE(PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_NO_SUBTYPE,PROBLEM_STATIC_FINITE_ELASTICITY_SUBTYPE,PROBLEM_DYNAMIC_FINITE_ELASTICITY_SUBTYPE, &
              & PROBLEM_FINITE_ELASTICITY_WITH_ACTIVE_SUBTYPE)
              !Do nothing
            CASE(PROBLEM_QUASISTATIC_FINITE_ELASTICITY_SUBTYPE)
              !Get the CellML evaluator solver
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,CELLML_SOLVER,err,error,*999)
              !Get the CellML equations for the CellML evaluator solver
              CALL SOLVER_CELLML_EQUATIONS_GET(CELLML_SOLVER,CELLML_EQUATIONS,err,error,*999)
              !Finish the CellML equations creation
              CALL CELLML_EQUATIONS_CREATE_FINISH(CELLML_EQUATIONS,err,error,*999)
            CASE(PROBLEM_QUASISTATIC_FINITE_ELASTICITY_WITH_GROWTH_SUBTYPE)
              !Get the CellML BC evaluator solver
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,CELLML_SOLVER,err,error,*999)
              !Get the CellML equations for the CellML evaluator solver
              CALL SOLVER_CELLML_EQUATIONS_GET(CELLML_SOLVER,CELLML_EQUATIONS,err,error,*999)
              !Finish the CellML equations creation
              CALL CELLML_EQUATIONS_CREATE_FINISH(CELLML_EQUATIONS,err,error,*999)
              !Get the CellML growth integration solver
              NULLIFY(CELLML_SOLVER)
              CALL SOLVERS_SOLVER_GET(SOLVERS,2,CELLML_SOLVER,err,error,*999)
              !Get the CellML equations for the CellML integration solver
              NULLIFY(CELLML_EQUATIONS)
              CALL SOLVER_CELLML_EQUATIONS_GET(CELLML_SOLVER,CELLML_EQUATIONS,err,error,*999)
              !Finish the CellML equations creation
              CALL CELLML_EQUATIONS_CREATE_FINISH(CELLML_EQUATIONS,err,error,*999)
            CASE(PROBLEM_FINITE_ELASTICITY_WITH_CELLML_SUBTYPE)
              !Get the nonlinear solver
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              !Get the CellML evaluator solver
              CALL SOLVER_NEWTON_CELLML_SOLVER_GET(SOLVER,CELLML_SOLVER,err,error,*999)
              !Get the CellML equations for the CellML evaluator solver
              CALL SOLVER_CELLML_EQUATIONS_GET(CELLML_SOLVER,CELLML_EQUATIONS,err,error,*999)
              !Finish the CellML equations creation
              CALL CELLML_EQUATIONS_CREATE_FINISH(CELLML_EQUATIONS,err,error,*999)
            CASE(PROBLEM_FINITE_ELASTICITY_WITH_GROWTH_CELLML_SUBTYPE) 
              !Get the CellML integrator solver
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,CELLML_SOLVER,err,error,*999)
              !Get the CellML equations for the CellML evaluator solver
              CALL SOLVER_CELLML_EQUATIONS_GET(CELLML_SOLVER,CELLML_EQUATIONS,err,error,*999)
              !Finish the CellML equations creation
              CALL CELLML_EQUATIONS_CREATE_FINISH(CELLML_EQUATIONS,err,error,*999)
              NULLIFY(CELLML_SOLVER)
              NULLIFY(CELLML_EQUATIONS)
              !Get the nonlinear solver
              CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
              !Get the CellML evaluator solver
              CALL SOLVER_NEWTON_CELLML_SOLVER_GET(SOLVER,CELLML_SOLVER,err,error,*999)
              !Get the CellML equations for the CellML evaluator solver
              CALL SOLVER_CELLML_EQUATIONS_GET(CELLML_SOLVER,CELLML_EQUATIONS,err,error,*999)
              !Finish the CellML equations creation
              CALL CELLML_EQUATIONS_CREATE_FINISH(CELLML_EQUATIONS,err,error,*999)
            CASE(PROBLEM_MULTISCALE_FINITE_ELASTICITY_SUBTYPE)
              !Get the CellML equations for the first DAE solver
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              CALL SOLVER_CELLML_EQUATIONS_GET(SOLVER,CELLML_EQUATIONS,err,error,*999)
              !Finish the CellML equations creation
              CALL CELLML_EQUATIONS_CREATE_FINISH(CELLML_EQUATIONS,err,error,*999)
            CASE DEFAULT
              localError="The third problem specification of "//TRIM(NumberToVString(problem%specification(3),"*",err,error))// &
                & " is not valid for a finite elasticity type of an elasticity problem class."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a finite elasticity equation."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
            & " is invalid for a finite elasticity problem."
          CALL FlagError(LOCAL_ERROR,err,error,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="Problem subtype "//TRIM(NumberToVString(PROBLEM_SUBTYPE,"*",err,error))// &
          & " is not valid for a finite elasticity type of an elasticity problem class."
        CALL FlagError(LOCAL_ERROR,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    ENDIF

    EXITS("FINITE_ELASTICITY_PROBLEM_SETUP")
    RETURN
999 ERRORSEXITS("FINITE_ELASTICITY_PROBLEM_SETUP",err,error)
    RETURN 1
  END SUBROUTINE FINITE_ELASTICITY_PROBLEM_SETUP

  !
  !================================================================================================================================
  !

  !>Sets up the finite elasticity problem.
  SUBROUTINE FiniteElasticity_ContactProblemSetup(PROBLEM,PROBLEM_SETUP,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem set to setup a Laplace equation on.
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP,CONTROL_LOOP_ROOT
    TYPE(SOLVER_TYPE), POINTER :: nonlinearSolver,transformationSolver
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    INTEGER(INTG) :: PROBLEM_SUBTYPE

    ENTERS("FiniteElasticity_ContactProblemSetup",err,error,*999)

    NULLIFY(CONTROL_LOOP)
    NULLIFY(nonlinearSolver)
    NULLIFY(transformationSolver)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVERS)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(.NOT.ALLOCATED(PROBLEM%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(PROBLEM%SPECIFICATION,1)<3) THEN
        CALL FlagError("Problem specification must have three entries for a finite elasticity problem.",err,error,*999)
      END IF
      PROBLEM_SUBTYPE=PROBLEM%SPECIFICATION(3)
      SELECT CASE(PROBLEM_SUBTYPE)
      CASE(PROBLEM_FE_CONTACT_TRANSFORM_REPROJECT_SUBTYPE,PROBLEM_FE_CONTACT_TRANSFORM_SUBTYPE)
        SELECT CASE(PROBLEM_SETUP%SETUP_TYPE)
        CASE(PROBLEM_SETUP_INITIAL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Do nothing????
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Do nothing????
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a finite elasticity problem."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Set up a simple control loop: default is load increment type now
            CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,err,error,*999)
            CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Finish the control loops
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,err,error,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a finite elasticity problem."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVERS_TYPE)
          !Get the control loop
          CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Start the solvers creation
            CALL SOLVERS_CREATE_START(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL SOLVERS_NUMBER_SET(SOLVERS,2,err,error,*999)
            !Set the first solver to be a geometric transformation solver
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,transformationSolver,err,error,*999)
            CALL SOLVER_TYPE_SET(transformationSolver,SOLVER_GEOMETRIC_TRANSFORMATION_TYPE,err,error,*999)
            !Set the second solver to be a nonlinear solver
            CALL SOLVERS_SOLVER_GET(SOLVERS,2,nonlinearSolver,err,error,*999)
            CALL SOLVER_TYPE_SET(nonlinearSolver,SOLVER_NONLINEAR_TYPE,err,error,*999)
            !Set solver defaults
            CALL SOLVER_LIBRARY_TYPE_SET(nonlinearSolver,SOLVER_PETSC_LIBRARY,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the solvers
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(SOLVERS,err,error,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a finite elasticity problem."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            !Get the solver
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,2,nonlinearSolver,err,error,*999)
            !Create the solver equatgions
            CALL SOLVER_EQUATIONS_CREATE_START(nonlinearSolver,SOLVER_EQUATIONS,err,error,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,err,error,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            !Get the solver equations
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,2,nonlinearSolver,err,error,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(nonlinearSolver,SOLVER_EQUATIONS,err,error,*999)
            !Finish the solver equations creation
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,err,error,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a finite elasticity problem."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
            & " is invalid for a finite elasticity problem."
          CALL FlagError(LOCAL_ERROR,err,error,*999)
        END SELECT
      CASE(PROBLEM_FE_CONTACT_REPROJECT_SUBTYPE)
        SELECT CASE(PROBLEM_SETUP%SETUP_TYPE)
        CASE(PROBLEM_SETUP_INITIAL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Do nothing????
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Do nothing????
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a finite elasticity problem."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Set up a simple control loop: default is load increment type now
            CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,err,error,*999)
            CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Finish the control loops
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,err,error,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a finite elasticity problem."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVERS_TYPE)
          !Get the control loop
          CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Start the solvers creation
            CALL SOLVERS_CREATE_START(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL SOLVERS_NUMBER_SET(SOLVERS,1,err,error,*999)
            !Set the solver to be a nonlinear solver
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,nonlinearSolver,err,error,*999)
            CALL SOLVER_TYPE_SET(nonlinearSolver,SOLVER_NONLINEAR_TYPE,err,error,*999)
            !Set solver defaults
            CALL SOLVER_LIBRARY_TYPE_SET(nonlinearSolver,SOLVER_PETSC_LIBRARY,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the solvers
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(SOLVERS,err,error,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a finite elasticity problem."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            !Get the solver
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,nonlinearSolver,err,error,*999)
            !Create the solver equatgions
            CALL SOLVER_EQUATIONS_CREATE_START(nonlinearSolver,SOLVER_EQUATIONS,err,error,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,err,error,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            !Get the solver equations
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,nonlinearSolver,err,error,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(nonlinearSolver,SOLVER_EQUATIONS,err,error,*999)
            !Finish the solver equations creation
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,err,error,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a finite elasticity problem."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
            & " is invalid for a finite elasticity problem."
          CALL FlagError(LOCAL_ERROR,err,error,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="Problem subtype "//TRIM(NumberToVString(PROBLEM_SUBTYPE,"*",err,error))// &
          & " is not valid for a finite elasticity contact type of an elasticity problem class."
        CALL FlagError(LOCAL_ERROR,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    ENDIF

    EXITS("FiniteElasticity_ContactProblemSetup")
    RETURN
999 ERRORSEXITS("FiniteElasticity_ContactProblemSetup",err,error)
    RETURN 1
  END SUBROUTINE FiniteElasticity_ContactProblemSetup

  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a finite elasticity type problem.
  SUBROUTINE FiniteElasticity_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: problem !<A pointer to the problem to set the problem subtype for
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<the error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<the error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: problemSubtype

    ENTERS("FiniteElasticity_ProblemSpecificationSet",err,error,*999)

    IF(ASSOCIATED(problem)) THEN
      IF(SIZE(problemSpecification,1)<3) THEN
        !Default to no subtype if not set
        problemSubtype=PROBLEM_NO_SUBTYPE
      ELSE IF(SIZE(problemSpecification,1)==3) THEN
        problemSubtype=problemSpecification(3)
        SELECT CASE(problemSubtype)
        CASE(PROBLEM_NO_SUBTYPE)
          !ok
        CASE(PROBLEM_STATIC_FINITE_ELASTICITY_SUBTYPE)
          !ok
        CASE(PROBLEM_QUASISTATIC_FINITE_ELASTICITY_SUBTYPE)
          !ok
        CASE(PROBLEM_QUASISTATIC_FINITE_ELASTICITY_WITH_GROWTH_SUBTYPE)
          !ok
        CASE(PROBLEM_DYNAMIC_FINITE_ELASTICITY_SUBTYPE)
          !ok
        CASE(PROBLEM_FINITE_ELASTICITY_WITH_ACTIVE_SUBTYPE)
          !ok
        CASE(PROBLEM_FINITE_ELASTICITY_WITH_CELLML_SUBTYPE)
          !ok
        CASE(PROBLEM_MULTISCALE_FINITE_ELASTICITY_SUBTYPE)
          !ok
        CASE(PROBLEM_FINITE_ELASTICITY_WITH_GROWTH_CELLML_SUBTYPE)
          !ok
        CASE DEFAULT
          localError="The third problem specification of "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
            & " is not valid for a finite elasticity type of an elasticity problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        CALL FlagError("Finite elasticity problem specification may only have up to 3 entries.",err,error,*999)
      END IF
      IF(ALLOCATED(problem%specification)) THEN
        CALL FlagError("Problem specification is already allocated.",err,error,*999)
      ELSE
        ALLOCATE(problem%specification(3),stat=err)
        IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
      END IF
      problem%specification(1:3)=[PROBLEM_ELASTICITY_CLASS,PROBLEM_FINITE_ELASTICITY_TYPE,problemSubtype]
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    ENDIF

    EXITS("FiniteElasticity_ProblemSpecificationSet")
    RETURN
999 ERRORS("FiniteElasticity_ProblemSpecificationSet",err,error)
    EXITS("FiniteElasticity_ProblemSpecificationSet")
    RETURN 1

  END SUBROUTINE FiniteElasticity_ProblemSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets/changes the problem subtype for a finite elasticity contact type .
  SUBROUTINE FiniteElasticity_ContactProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: problem !<A pointer to the problem to set the problem subtype for
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: problemSubtype

    ENTERS("FiniteElasticity_ContactProblemSpecificationSet",err,error,*999)

    IF(ASSOCIATED(problem)) THEN
      IF(SIZE(problemSpecification,1)<3) THEN
        !Default to no subtype if not set
        problemSubtype=PROBLEM_NO_SUBTYPE
      ELSE IF(SIZE(problemSpecification,1)==3) THEN
        problemSubtype=problemSpecification(3)
        SELECT CASE(problemSubtype)
        CASE(PROBLEM_NO_SUBTYPE)
          !Normal finite elasticity problem subject to contact constraint, no extra solvers required
        CASE(PROBLEM_FE_CONTACT_TRANSFORM_REPROJECT_SUBTYPE)
          !ok
        CASE(PROBLEM_FE_CONTACT_TRANSFORM_SUBTYPE)
          !ok
        CASE(PROBLEM_FE_CONTACT_REPROJECT_SUBTYPE)
          !ok
        CASE DEFAULT
          localError="The third problem specification of "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
            & " is not valid for a finite elasticity contact type of an elasticity problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        CALL FlagError("Finite elasticity contact problem specification may only have up to 3 entries.",err,error,*999)
      END IF
      IF(ALLOCATED(problem%specification)) THEN
        CALL FlagError("Problem specification is already allocated.",err,error,*999)
      ELSE
        ALLOCATE(problem%specification(3),stat=err)
        IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
      END IF
      problem%specification(1:3)=[PROBLEM_ELASTICITY_CLASS,PROBLEM_FINITE_ELASTICITY_CONTACT_TYPE,problemSubtype]
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    ENDIF

    EXITS("FiniteElasticity_ContactProblemSpecificationSet")
    RETURN
999 ERRORS("FiniteElasticity_ContactProblemSpecificationSet",err,error)
    EXITS("FiniteElasticity_ContactProblemSpecificationSet")
    RETURN 1

  END SUBROUTINE FiniteElasticity_ContactProblemSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the finite elasticity problem post solve.
  SUBROUTINE FiniteElasticity_PostSolve(solver,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver!<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx,I
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    TYPE(FIELD_TYPE), POINTER :: dependentField,independentField
    TYPE(PROBLEM_TYPE), POINTER :: problem
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping
 
    ENTERS("FiniteElasticity_PostSolve",err,error,*999)
    
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    IF(.NOT.ALLOCATED(problem%specification)) CALL FlagError("Problem specification is not allocated.",err,error,*999)
    IF(SIZE(problem%specification,1)<3) &
      & CALL FlagError("Problem specification must have three entries for a finite elasticity problem.",err,error,*999)
    
    SELECT CASE(problem%specification(3))
    CASE(PROBLEM_QUASISTATIC_FINITE_ELASTICITY_SUBTYPE)
      IF(solver%GLOBAL_NUMBER==2) THEN
        !Nonlinear solver
        CALL Solver_NonlinearDivergenceExit(solver,err,error,*999)
        !Update the hardening variable if we have converged
        NULLIFY(solverEquations)
        CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
        NULLIFY(solverMapping)
        CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
        DO equationsSetIdx=1,solverMapping%NUMBER_OF_EQUATIONS_SETS
          NULLIFY(equationsSet)
          CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
          IF(equationsSet%specification(1)==EQUATIONS_SET_ELASTICITY_CLASS.AND. &
            & equationsSet%specification(2)==EQUATIONS_SET_FINITE_ELASTICITY_TYPE) THEN
            IF(equationsSet%specification(3)==EQUATIONS_SET_RATE_BASED_SMOOTH_MODEL_SUBTYPE.OR. &
              & equationsSet%specification(3)==EQUATIONS_SET_COMPRESSIBLE_RATE_BASED_SMOOTH_MODEL_SUBTYPE) THEN
              !Rate based smooth model
              NULLIFY(dependentField)
              CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
              !Update U1 variable
              CALL Field_ParameterSetsCopy(dependentField,FIELD_U1_VARIABLE_TYPE,FIELD_NEXT_VALUES_SET_TYPE, &
                & FIELD_VALUES_SET_TYPE,1.0_DP,err,error,*999)
            ELSE IF(equationsSet%specification(3)==EQUATIONS_SET_RATE_BASED_GROWTH_MODEL_SUBTYPE.OR. &
              & equationsSet%specification(3)==EQUATIONS_SET_COMPRESSIBLE_RATE_BASED_GROWTH_MODEL_SUBTYPE) THEN
              !Rate based growh model
              NULLIFY(dependentField)
              CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
              !Update U1 variable
              CALL Field_ParameterSetsCopy(dependentField,FIELD_U1_VARIABLE_TYPE,FIELD_NEXT_VALUES_SET_TYPE, &
                & FIELD_VALUES_SET_TYPE,1.0_DP,err,error,*999)
              !Update U2 variable
              CALL Field_ParameterSetsCopy(dependentField,FIELD_U2_VARIABLE_TYPE,FIELD_NEXT_VALUES_SET_TYPE, &
                & FIELD_VALUES_SET_TYPE,1.0_DP,err,error,*999)
            ENDIF
          ENDIF
        ENDDO !equationsSetIdx
        !Output results
        CALL FiniteElasticity_PostSolveOutputData(solver,err,error,*999)          
      ENDIF
    CASE(PROBLEM_QUASISTATIC_FINITE_ELASTICITY_WITH_GROWTH_SUBTYPE)
      !Add in other quasi-static stuff
      IF(solver%GLOBAL_NUMBER==3) THEN
        !Nonlinear solver
        CALL Solver_NonlinearDivergenceExit(solver,err,error,*999)
        CALL FiniteElasticity_PostSolveOutputData(solver,err,error,*999)          
      ENDIF
    CASE(PROBLEM_FINITE_ELASTICITY_WITH_GROWTH_CELLML_SUBTYPE)
      IF(solver%GLOBAL_NUMBER==2) THEN
        !Nonlinear solver
        CALL Solver_NonlinearDivergenceExit(solver,err,error,*999)
        CALL FiniteElasticity_PostSolveOutputData(solver,err,error,*999)          
      ENDIF
    CASE(PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE,PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE, &
      & PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
      !Call divergence test only if finite element loop: THIS IS NOT A PROPER FIX
      IF(controlLoop%SUB_LOOP_INDEX==1) THEN
        CALL Solver_NonlinearDivergenceExit(solver,err,error,*999)
      ENDIF
      IF(controlLoop%LOOP_TYPE==PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE.AND.solver%GLOBAL_NUMBER==1) THEN
        CALL FiniteElasticity_PostSolveOutputData(solver,err,error,*999)
      END IF
    CASE(PROBLEM_FINITE_ELASTICITY_WITH_ACTIVE_SUBTYPE)
      ! how to check eqn subtype? assume active contraction
      NULLIFY(solverEquations)
      CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
      NULLIFY(independentField)
      CALL EquationsSet_IndependentFieldGet(equationsSet,independentField,err,error,*999)
      ! store lambda Q (5-8) in prev lambda Q (1-4)
      DO I=1,4
        CALL Field_ParametersToFieldParametersCopy(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,I+4, &
          & independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,I,err,error,*999)
      END DO
      ! output data
      CALL FiniteElasticity_PostSolveOutputData(solver,err,error,*999)
    CASE(PROBLEM_MULTISCALE_FINITE_ELASTICITY_SUBTYPE)
      IF(ASSOCIATED(solver%DAE_SOLVER)) THEN
        !do nothing
      ELSE IF(ASSOCIATED(solver%NONLINEAR_SOLVER)) THEN
        CALL Solver_NonlinearDivergenceExit(solver,err,error,*999)
      END IF
    CASE DEFAULT
      !Check that solver converged
      CALL Solver_NonlinearDivergenceExit(solver,err,error,*999)
    END SELECT

    EXITS("FiniteElasticity_PostSolve")
    RETURN
999 ERRORSEXITS("FiniteElasticity_PostSolve",err,error)
    RETURN 1
    
  END SUBROUTINE FiniteElasticity_PostSolve

  !
  !================================================================================================================================
  !

  !>Sets up the finite elasticity problem post solve output data.
  SUBROUTINE FiniteElasticity_PostSolveOutputData(solver,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop !<A pointer to the control loop to solve.
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS,solverEquations  !<A pointer to the solver equations
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET,equationsSet !<A pointer to the equations set
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING,solverMapping !<A pointer to the solver mapping
    TYPE(CONTROL_LOOP_TYPE), POINTER :: TIME_LOOP !<A pointer to the control time loop.
    TYPE(PROBLEM_TYPE), POINTER :: problem
    TYPE(VARYING_STRING) :: LOCAL_ERROR,localError
    TYPE(VARYING_STRING) :: METHOD,outputFile !,FILE
    CHARACTER(14) :: file
    CHARACTER(14) :: OUTPUT_FILE
    CHARACTER(15) :: fileName
    LOGICAL :: EXPORT_FIELD
    REAL(DP) :: startTime,stopTime,currentTime,timeIncrement
    INTEGER(INTG) :: CURRENT_LOOP_ITERATION,currentIteration
    INTEGER(INTG) :: OUTPUT_ITERATION_NUMBER,outputIteration
    INTEGER(INTG) :: equations_set_idx,loop_idx,equationsSetIdx

    ENTERS("FiniteElasticity_PostSolveOutputData",err,error,*999)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
        
    IF(.NOT.ALLOCATED(problem%SPECIFICATION)) &
      & CALL FlagError("Problem specification is not allocated.",err,error,*999)
    IF(SIZE(problem%SPECIFICATION,1)<3) &
      & CALL FlagError("Problem specification must have three entries for a finite elasticity problem.",err,error,*999)
     
    SELECT CASE(problem%SPECIFICATION(3))
    CASE(PROBLEM_NO_SUBTYPE,PROBLEM_STATIC_FINITE_ELASTICITY_SUBTYPE, &
      & PROBLEM_STANDARD_ELASTICITY_FLUID_PRESSURE_SUBTYPE, &
      & PROBLEM_FINITE_ELASTICITY_WITH_CELLML_SUBTYPE, &
      & PROBLEM_FINITE_ELASTICITY_WITH_GROWTH_CELLML_SUBTYPE, &
      & PROBLEM_MULTISCALE_FINITE_ELASTICITY_SUBTYPE)
      SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
      IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
        SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
        IF(ASSOCIATED(SOLVER_MAPPING)) THEN
          !Make sure the equations sets are up to date
          DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
            EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%ptr
            METHOD="FORTRAN"
            !EXPORT_FIELD=.TRUE.
            EXPORT_FIELD=.FALSE.
            IF(EXPORT_FIELD) THEN          
              IF(SOLVER%outputType>=SOLVER_PROGRESS_OUTPUT) THEN
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Finite Elasticity export fields ... ",err,error,*999)
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"STATICSOLUTION",err,error,*999)
              ENDIF
              CALL FLUID_MECHANICS_IO_WRITE_CMGUI(EQUATIONS_SET%REGION,EQUATIONS_SET%GLOBAL_NUMBER, &
                & "STATICSOLIDSOLUTION",err,error,*999)
            ENDIF
          ENDDO
        ENDIF
      ENDIF
    CASE(PROBLEM_QUASISTATIC_FINITE_ELASTICITY_SUBTYPE,PROBLEM_QUASISTATIC_FINITE_ELASTICITY_WITH_GROWTH_SUBTYPE, &
      & PROBLEM_DYNAMIC_FINITE_ELASTICITY_SUBTYPE)
      !Get the current time information
      CALL ControlLoop_CurrentTimeInformationGet(controlLoop,currentTime,timeIncrement,startTime,stopTime,currentIteration, &
        & outputIteration,err,error,*999)
      !See if we should output
      IF(currentTime<=stopTime) THEN
        IF(outputIteration/=0) THEN
          IF(MOD(currentIteration,outputIteration)==0) THEN
            !Output 
            WRITE(fileName,'("Elasticity_",I0)') currentIteration
            outputFile=fileName(1:LEN_TRIM(filename))
            method="FORTRAN"
            !Loop over the equations sets and output
            NULLIFY(solverEquations)
            CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
            NULLIFY(solverMapping)
            CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
            DO equationsSetIdx=1,solverMapping%NUMBER_OF_EQUATIONS_SETS
              NULLIFY(equationsSet)
              CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
              IF(.NOT.ALLOCATED(equationsSet%specification))  &
                & CALL FlagError("Equations set specification is not allocated.",err,error,*999)
              IF(SIZE(equationsSet%specification,1)<2)  &
                & CALL FlagError("Equations set specification does not have a type.",err,error,*999)
              IF(equationsSet%specification(2)==EQUATIONS_SET_FINITE_ELASTICITY_TYPE) THEN
                CALL FIELD_IO_NODES_EXPORT(equationsSet%region%fields,outputFile,method,err,error,*999)
                CALL FIELD_IO_ELEMENTS_EXPORT(equationsSet%region%fields,outputFile,method,err,error,*999)
              ENDIF
            ENDDO !equationsSetIdx
          ENDIF
        ENDIF
      ENDIF
    CASE(PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE,PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE, &
      & PROBLEM_FINITE_ELASTICITY_WITH_ACTIVE_SUBTYPE,PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE, &
      & PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
      SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
      IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
        SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
        IF(ASSOCIATED(SOLVER_MAPPING)) THEN
          !Make sure the equations sets are up to date
          DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
            EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%ptr
            IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
              CALL FlagError("Equations set specification is not allocated.",err,error,*999)
            ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<2) THEN
              CALL FlagError("Equations set specification does not have a type.", &
                & err,error,*999)
            END IF
            IF(EQUATIONS_SET%SPECIFICATION(2)==EQUATIONS_SET_FINITE_ELASTICITY_TYPE) THEN
              TIME_LOOP=>controlLoop !Initialise time loop (load increment loop on first)
              !Move up to find outer time loop
              DO loop_idx=1,controlLoop%CONTROL_LOOP_LEVEL-1
                IF(ASSOCIATED(TIME_LOOP%PARENT_LOOP)) THEN
                  TIME_LOOP=>TIME_LOOP%PARENT_LOOP
                ELSE
                  CALL FlagError("Could not find a time control loop.",err,error,*999)
                ENDIF
              ENDDO
              CURRENT_LOOP_ITERATION=TIME_LOOP%TIME_LOOP%ITERATION_NUMBER
              OUTPUT_ITERATION_NUMBER=TIME_LOOP%TIME_LOOP%OUTPUT_NUMBER
              
              !Write out fields at each timestep
              IF(TIME_LOOP%TIME_LOOP%CURRENT_TIME<=TIME_LOOP%TIME_LOOP%STOP_TIME) THEN
                WRITE(OUTPUT_FILE,'("S_TIMESTP_",I4.4)') CURRENT_LOOP_ITERATION
                FILE=OUTPUT_FILE
                METHOD="FORTRAN"
                EXPORT_FIELD=.TRUE.
                IF(EXPORT_FIELD) THEN
                  IF(OUTPUT_ITERATION_NUMBER/=0.AND.MOD(CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER)==0)  THEN
                    IF(SOLVER%outputType>=SOLVER_PROGRESS_OUTPUT) THEN
                      CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Finite Elasticity export fields ...",err,error,*999)
                      CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,OUTPUT_FILE,err,error,*999)
                    ENDIF
                    CALL FLUID_MECHANICS_IO_WRITE_CMGUI(EQUATIONS_SET%REGION,EQUATIONS_SET%GLOBAL_NUMBER,FILE, &
                      & err,error,*999)
                    IF(SOLVER%outputType>=SOLVER_PROGRESS_OUTPUT) THEN
                      CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Finite Elasticity all fields exported ...",err,error,*999)
                    ENDIF
                    CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,OUTPUT_FILE,err,error,*999)
                  ENDIF
                ENDIF
              ENDIF !stop_time                 
            ENDIF !EQUATIONS_SET_FINITE_ELASTICITY_TYPE
          ENDDO !equations_set_idx
        ENDIF !Solver_mapping
      ENDIF !Solver_equations
    CASE DEFAULT
      LOCAL_ERROR="The third problem specification of "// &
        & TRIM(NumberToVString(problem%SPECIFICATION(3),"*",err,error))// &
        & " is not valid for a finite elasticity problem class."
      CALL FlagError(LOCAL_ERROR,err,error,*999)
    END SELECT
      
    EXITS("FiniteElasticity_PostSolveOutputData")
    RETURN
999 ERRORSEXITS("FiniteElasticity_PostSolveOutputData",err,error)
    RETURN 1
    
  END SUBROUTINE FiniteElasticity_PostSolveOutputData

  !
  !================================================================================================================================
  !

  !>Runs before each time loop for a finite elasticity problem.
  SUBROUTINE FiniteElasticity_ControlTimeLoopPreLoop(CONTROL_LOOP,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the time control loop
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER_SOLID !<A pointer to the solid solver
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP_SOLID
    TYPE(FIELD_TYPE), POINTER :: INDEPENDENT_FIELD

    NULLIFY(SOLVER_SOLID)
    NULLIFY(CONTROL_LOOP_SOLID)

    ENTERS("FiniteElasticity_ControlTimeLoopPreLoop",err,error,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
        IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
          CALL FlagError("Problem specification is not allocated.",err,error,*999)
        ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
          CALL FlagError("Problem specification must have three entries for a finite elasticity problem.",err,error,*999)
        END IF

!!TODO: ALL THE BELOW SHOULD BE RENAMED PRE-LOOP (NOT PRE-SOLVE) AND HAVE THE CONTROL LOOP ONLY PASSED.
        
        SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
        CASE(PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE,PROBLEM_STANDARD_ELASTICITY_FLUID_PRESSURE_SUBTYPE)
          ! could do this in one line with problem_solver_get but the dependence on problem_routines causes a circular dependence
          CALL CONTROL_LOOP_GET(CONTROL_LOOP,[1,CONTROL_LOOP_NODE],CONTROL_LOOP_SOLID,err,error,*999)
          CALL SOLVERS_SOLVER_GET(CONTROL_LOOP_SOLID%SOLVERS,1,SOLVER_SOLID,err,error,*999)
          !--- 3.0 For Standard Elasticity Darcy: Update the boundary conditions of the solid
          CALL FiniteElasticity_PreSolveUpdateBoundaryConditions(SOLVER_SOLID,err,error,*999)
        CASE(PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
          CALL CONTROL_LOOP_GET(CONTROL_LOOP,[1,1,CONTROL_LOOP_NODE],CONTROL_LOOP_SOLID,err,error,*999)
          CALL SOLVERS_SOLVER_GET(CONTROL_LOOP_SOLID%SOLVERS,1,SOLVER_SOLID,err,error,*999)
            !--- 3.0 For Standard Elasticity Darcy: Update the boundary conditions of the solid
          CALL FiniteElasticity_PreSolveUpdateBoundaryConditions(SOLVER_SOLID,err,error,*999)
        CASE(PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE)
          CALL CONTROL_LOOP_GET(CONTROL_LOOP,[1,CONTROL_LOOP_NODE],CONTROL_LOOP_SOLID,err,error,*999)
          CALL SOLVERS_SOLVER_GET(CONTROL_LOOP_SOLID%SOLVERS,1,SOLVER_SOLID,err,error,*999)
          !--- For PGM: Get the displacement field
          CALL FiniteElasticity_PreSolveGetSolidDisplacement(CONTROL_LOOP,SOLVER_SOLID,err,error,*999)
        CASE DEFAULT
          !do nothing
        END SELECT
      ELSE
        CALL FlagError("Problem is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",err,error,*999)
    ENDIF

    EXITS("FiniteElasticity_ControlTimeLoopPreLoop")
    RETURN
999 ERRORSEXITS("FiniteElasticity_ControlTimeLoopPreLoop",err,error)
    RETURN 1

  END SUBROUTINE FiniteElasticity_ControlTimeLoopPreLoop

  !
  !================================================================================================================================
  !

  !>Executes after each loop of a control loop for finite elasticity problems, i.e., after each load increment in a load increment loop
  SUBROUTINE FiniteElasticity_ControlLoadIncrementLoopPostLoop(controlLoop,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop !<A pointer to the control loop
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(SOLVERS_TYPE), POINTER :: solvers
    TYPE(SOLVER_TYPE), POINTER :: solver
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping
    TYPE(REGION_TYPE), POINTER :: region
    TYPE(FIELDS_TYPE), POINTER :: fields
    INTEGER(INTG) :: solverIdx,equationsSetIdx,incrementIdx,outputNumber
    LOGICAL :: dirExist
    TYPE(VARYING_STRING) :: fileName,method,directory

    ENTERS("FiniteElasticity_ControlLoadIncrementLoopPostLoop",err,error,*999)

    IF(ASSOCIATED(controlLoop)) THEN
      IF(controlLoop%LOOP_TYPE==PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE) THEN
        incrementIdx=controlLoop%LOAD_INCREMENT_LOOP%ITERATION_NUMBER
        outputNumber=controlLoop%LOAD_INCREMENT_LOOP%OUTPUT_NUMBER
        IF(outputNumber>0) THEN
          IF(MOD(incrementIdx,outputNumber)==0) THEN
            solvers=>controlLoop%SOLVERS
            IF(ASSOCIATED(solvers)) THEN
              DO solverIdx=1,solvers%NUMBER_OF_SOLVERS
                solver=>solvers%SOLVERS(solverIdx)%ptr
                IF(ASSOCIATED(solver)) THEN
                  solverEquations=>SOLVER%SOLVER_EQUATIONS
                  IF(ASSOCIATED(solverEquations)) THEN
                    solverMapping=>SOLVER%SOLVER_EQUATIONS%SOLVER_MAPPING
                    IF(ASSOCIATED(solverMapping)) THEN
                      DO equationsSetIdx=1,solverMapping%NUMBER_OF_EQUATIONS_SETS
                        region=>solverMapping%EQUATIONS_SETS(equationsSetIdx)%ptr%REGION
                        NULLIFY(fields)
                        fields=>region%FIELDS
                        directory="results_load/"
                        INQUIRE(FILE=CHAR(directory),EXIST=dirExist)
                        IF(.NOT.dirExist) THEN
                          CALL SYSTEM(CHAR("mkdir "//directory))
                        ENDIF
                        fileName=directory//"mesh"//TRIM(NumberToVString(equationsSetIdx,"*",err,error))// &
                          & "_load"//TRIM(NumberToVString(incrementIdx,"*",err,error))
                        method="FORTRAN"
                        CALL FIELD_IO_ELEMENTS_EXPORT(fields,fileName,method,err,error,*999)
                        CALL FIELD_IO_NODES_EXPORT(fields,fileName,method,err,error,*999)
                      ENDDO !equationsSetIdx
                    ENDIF
                  ENDIF
                ENDIF
              ENDDO !solverIdx
            ELSE
              CALL FlagError("Control loop solvers is not associated.",err,error,*999)
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",err,error,*999)
    ENDIF

    EXITS("FiniteElasticity_ControlLoadIncrementLoopPostLoop")
    RETURN
999 ERRORS("FiniteElasticity_ControlLoadIncrementLoopPostLoop",err,error)
    EXITS("FiniteElasticity_ControlLoadIncrementLoopPostLoop")
    RETURN 1

  END SUBROUTINE FiniteElasticity_ControlLoadIncrementLoopPostLoop

  !
  !================================================================================================================================
  !

  !>Sets up the finite elasticity problem pre-solve.
  SUBROUTINE FiniteElasticity_PreSolve(solver,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: solverMatrixIdx,equationsSetIdx
    LOGICAL :: cellMLSolve,nonlinearSolve,validSubType
    REAL(DP) :: currentTime,timeIncrement
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop !<A pointer to the control loop to solve.
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    TYPE(FIELD_TYPE), POINTER :: dependentField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: dependentVariable
    TYPE(PROBLEM_TYPE), POINTER :: problem
    TYPE(SOLVER_TYPE), POINTER :: cellMLSolver
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping !<A pointer to the solver mapping
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: solverMatrices
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: solverMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("FiniteElasticity_PreSolve",err,error,*999)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    
    IF(.NOT.ALLOCATED(problem%specification)) CALL FlagError("Problem specification is not allocated.",err,error,*999)
    IF(SIZE(problem%specification,1)<3) &
      & CALL FlagError("Problem specification must have three entries for a finite elasticity problem.",err,error,*999)
 
    SELECT CASE(problem%specification(3))
    CASE(PROBLEM_NO_SUBTYPE)
      !Do nothing ???
    CASE(PROBLEM_STATIC_FINITE_ELASTICITY_SUBTYPE)
      !Do nothing ???
    CASE(PROBLEM_QUASISTATIC_FINITE_ELASTICITY_SUBTYPE)
      !Do nothing ???
    CASE(PROBLEM_QUASISTATIC_FINITE_ELASTICITY_WITH_GROWTH_SUBTYPE)
      IF(solver%GLOBAL_NUMBER==2) THEN
        CALL ControlLoop_CurrentTimesGet(controlLoop,currentTime,timeIncrement,err,error,*999)
        CALL Solver_DAETimesSet(solver,currentTime,currentTime+timeIncrement,err,error,*999)
      ENDIF
    CASE(PROBLEM_DYNAMIC_FINITE_ELASTICITY_SUBTYPE)
      !Do nothing ???
    CASE(PROBLEM_FINITE_ELASTICITY_WITH_CELLML_SUBTYPE, &
      & PROBLEM_FINITE_ELASTICITY_WITH_GROWTH_CELLML_SUBTYPE)
      IF(problem%specification(3)==PROBLEM_FINITE_ELASTICITY_WITH_GROWTH_CELLML_SUBTYPE) THEN
        IF(solver%GLOBAL_NUMBER==1) THEN
          cellMLSolve=.TRUE.
          nonlinearSolve=.FALSE.
        ELSE
          cellMLSolve=.FALSE.
          nonlinearSolve=.TRUE.
        ENDIF
      ELSE
        cellMLSolve=.FALSE.
        nonlinearSolve=.TRUE.
      ENDIF
      IF(cellMLSolve) THEN
        CALL ControlLoop_CurrentTimesGet(controlLoop,currentTime,timeIncrement,err,error,*999)
        CALL Solver_DAETimesSet(solver,currentTime,currentTime+timeIncrement,err,error,*999)
      ENDIF
      IF(nonlinearSolve) THEN
        NULLIFY(solverEquations)
        CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
        NULLIFY(solverMapping)
        CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
        validSubType=.FALSE.
        DO equationsSetIdx=1,solverMapping%NUMBER_OF_EQUATIONS_SETS
          NULLIFY(equationsSet)
          CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
          IF(equationsSet%specification(3)==EQUATIONS_SET_CONSTITUTIVE_LAW_IN_CELLML_EVALUATE_SUBTYPE.OR. &
            equationsSet%specification(3)==EQUATIONS_SET_CONSTITUTIVE_AND_GROWTH_LAW_IN_CELLML_SUBTYPE) THEN
            validSubtype=.TRUE.
            !compute the strain field
            NULLIFY(dependentField)
            CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
            NULLIFY(dependentVariable)
            CALL Field_VariableGet(dependentField,FIELD_U1_VARIABLE_TYPE,dependentVariable,err,error,*999)
            CALL FiniteElasticity_StressStrainCalculate(equationsSet,EQUATIONS_SET_R_CAUCHY_GREEN_DEFORMATION_TENSOR, &
              & dependentVariable,err,error,*999)
            !check for a linked CellML solver
            cellMLSolver=>solver%NONLINEAR_SOLVER%NEWTON_SOLVER%CELLML_EVALUATOR_SOLVER
            IF(ASSOCIATED(cellMLSolver)) THEN
              !evaluate the constiutive equation in CellML
              CALL Solver_Solve(cellMLSolver,err,error,*999)
            ENDIF
          ENDIF
        ENDDO !equationsSetIdx
      ENDIF
    CASE(PROBLEM_FINITE_ELASTICITY_WITH_ACTIVE_SUBTYPE)
      ! do nothing, time values get updated in CONTROL_TIME_LOOP_PRE_LOOP as there might be 
      ! a load increment loop below the time loop, so we don't want to update times here before
      ! every solve
    CASE(PROBLEM_GUDUNOV_MONODOMAIN_SIMPLE_ELASTICITY_SUBTYPE)
      ! do nothing
    CASE(PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE,PROBLEM_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE)
      ! do nothing
    CASE(EQUATIONS_SET_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE)
      ! do nothing
    CASE(PROBLEM_MULTISCALE_FINITE_ELASTICITY_SUBTYPE)
      ! do nothing
    CASE(PROBLEM_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE)
      !evaluate the evolution law using the cell model variables of the current time step and the deformation gradient tensor of the previous time step
      CALL FINITE_ELASTICITY_EVALUATE_EVOLUTION_LAW(SOLVER,err,error,*999)
    CASE(PROBLEM_STANDARD_ELASTICITY_FLUID_PRESSURE_SUBTYPE)
      ! do nothing
    CASE(PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE,PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE, &
      & PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
      
      !--- Set 'solverMatrix%updateMatrix=.TRUE.'
      NULLIFY(solverEquations)
      CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      NULLIFY(solverMatrices)
      CALL SolverEquations_SolverMatricesGet(solverEquations,solverMatrices,err,error,*999)
      
      DO solverMatrixIdx=1,solverMapping%NUMBER_OF_SOLVER_MATRICES
        NULLIFY(solverMatrix)
        CALL SolverMatrices_SolverMatrixGet(solverMatrices,solverMatrixIdx,solverMatrix,err,error,*999)
        solverMatrix%UPDATE_MATRIX=.TRUE.
      ENDDO !solverMatrixIdx
      
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(problem%specification(3),"*",err,error))// &
        & " is not valid for a finite elasticity problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("FiniteElasticity_PreSolve")
    RETURN
999 ERRORSEXITS("FiniteElasticity_PreSolve",err,error)
    RETURN 1
    
  END SUBROUTINE FiniteElasticity_PreSolve
      
  !   
  !================================================================================================================================
  !

  !>Evaluates the evolution law of a multiscale active strain muscle model
  SUBROUTINE FINITE_ELASTICITY_EVALUATE_EVOLUTION_LAW(SOLVER,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    INTEGER(INTG) :: solver_matrix_idx,equations_set_idx
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,INDEPENDENT_FIELD,FIBRE_FIELD,GEOMETRIC_FIELD,MATERIALS_FIELD
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING
    INTEGER(INTG) :: gauss_idx,element_idx,ne
    INTEGER(INTG) :: DEPENDENT_NUMBER_OF_GAUSS_POINTS
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: DEPENDENT_QUADRATURE_SCHEME
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: GEOMETRIC_INTERPOLATION_PARAMETERS,INDEPENDENT_INTERPOLATION_PARAMETERS, &
      & FIBRE_INTERPOLATION_PARAMETERS,MATERIALS_INTERPOLATION_PARAMETERS,DEPENDENT_INTERPOLATION_PARAMETERS
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: GEOMETRIC_INTERPOLATED_POINT,FIBRE_INTERPOLATED_POINT, &
      & MATERIALS_INTERPOLATED_POINT,DEPENDENT_INTERPOLATED_POINT,INDEPENDENT_INTERPOLATED_POINT
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: GEOMETRIC_INTERPOLATED_POINT_METRICS, &
      & DEPENDENT_INTERPOLATED_POINT_METRICS
    TYPE(BASIS_TYPE), POINTER :: GEOMETRIC_BASIS
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    INTEGER(INTG) :: FIELD_VAR_TYPE
    INTEGER(INTG) :: dof_idx,idx,i,j,LWORK
    INTEGER(INTG) :: MESH_COMPONENT_NUMBER
    REAL(DP) :: DZDXI(3,3),DZDNU(3,3),DZDNUT(3,3),TEMP(3,3)
    REAL(DP), DIMENSION (:), POINTER :: C !Parameters for constitutive laws
    REAL(DP) :: VALUE
    REAL(DP) :: TOL,TOL1,UP,LOW
    REAL(DP) :: F_e(3,3),F_a_inv(3,3),F_a_inv_T(3,3),F_a_T(3,3),C_a(3,3),C_a_inv(3,3),lambda_a,C_e(3,3),F_e_T(3,3)
    REAL(DP) :: REFERENCE_VOLUME,XB_STIFFNESS,XB_DISTORTION
    REAL(DP) :: SARCO_LENGTH,FREE_ENERGY,FREE_ENERGY_0,XB_ENERGY_PER_VOLUME,SLOPE,lambda_f,A_1,A_2,x_1,x_2
    REAL(DP) :: MAX_XB_NUMBER_PER_VOLUME,ENERGY_PER_XB,FORCE_LENGTH,I_1e,EVALUES(3),EVECTOR_1(3),EVECTOR_2(3),EVECTOR_3(3)
    REAL(DP) :: EMATRIX_1(3,3),EMATRIX_2(3,3),EMATRIX_3(3,3),TEMP1(3,3),TEMP2(3,3),TEMP3(3,3),N1(3,3),N2(3,3),N3(3,3) 
    INTEGER(INTG), PARAMETER :: LWMAX=1000
    REAL(DP) :: WORK(LWMAX)
    
    ENTERS("FINITE_ELASTICITY_EVALUATE_EVOLUTION_LAW",err,error,*999)

    NULLIFY(ELEMENTS_MAPPING)
    NULLIFY(DECOMPOSITION)
    NULLIFY(DEPENDENT_BASIS,GEOMETRIC_BASIS)
    NULLIFY(EQUATIONS)
    NULLIFY(DEPENDENT_FIELD,FIBRE_FIELD,GEOMETRIC_FIELD,MATERIALS_FIELD,INDEPENDENT_FIELD)
    NULLIFY(FIELD_VARIABLE)
    NULLIFY(DEPENDENT_QUADRATURE_SCHEME)
    NULLIFY(GEOMETRIC_INTERPOLATION_PARAMETERS,FIBRE_INTERPOLATION_PARAMETERS)
    NULLIFY(MATERIALS_INTERPOLATION_PARAMETERS,DEPENDENT_INTERPOLATION_PARAMETERS)
    NULLIFY(INDEPENDENT_INTERPOLATION_PARAMETERS)
    NULLIFY(GEOMETRIC_INTERPOLATED_POINT,FIBRE_INTERPOLATED_POINT)
    NULLIFY(GEOMETRIC_INTERPOLATED_POINT_METRICS,DEPENDENT_INTERPOLATED_POINT_METRICS)
    NULLIFY(MATERIALS_INTERPOLATED_POINT,DEPENDENT_INTERPOLATED_POINT)
    NULLIFY(INDEPENDENT_INTERPOLATED_POINT)

    !compute the deformation gradient tensor at the Gauss point
    SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
    IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
      SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
      IF(ASSOCIATED(SOLVER_MAPPING)) THEN
        DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
          EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%ptr
          EQUATIONS=>EQUATIONS_SET%EQUATIONS
          NULLIFY(vectorEquations)
          CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
          
          FIBRE_FIELD      =>equations%interpolation%fibreField
          GEOMETRIC_FIELD  =>equations%interpolation%geometricField
          MATERIALS_FIELD  =>equations%interpolation%materialsField
          DEPENDENT_FIELD  =>equations%interpolation%dependentField
          INDEPENDENT_FIELD=>equations%interpolation%independentField
          
          IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_1D3D_MONODOMAIN_ACTIVE_STRAIN_SUBTYPE) THEN

            DECOMPOSITION=>DEPENDENT_FIELD%DECOMPOSITION

            ELEMENTS_MAPPING=>DECOMPOSITION%DOMAIN(DECOMPOSITION%MESH_COMPONENT_NUMBER)% &
              & PTR%MAPPINGS%ELEMENTS

            DO element_idx=ELEMENTS_MAPPING%INTERNAL_START,ELEMENTS_MAPPING%INTERNAL_FINISH
              ne=ELEMENTS_MAPPING%DOMAIN_LIST(element_idx)

              MESH_COMPONENT_NUMBER=DECOMPOSITION%MESH_COMPONENT_NUMBER

              DEPENDENT_BASIS=>DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%ptr%TOPOLOGY%ELEMENTS% &
                & ELEMENTS(ne)%BASIS
              DEPENDENT_QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE% &
                & QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
              DEPENDENT_NUMBER_OF_GAUSS_POINTS=DEPENDENT_QUADRATURE_SCHEME%NUMBER_OF_GAUSS
              GEOMETRIC_BASIS=>GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%&
                & PTR%TOPOLOGY%ELEMENTS%ELEMENTS(ne)%BASIS

              !Initialise tensors and matrices
              DZDNU=0.0_DP
              DO idx=1,3
                DZDNU(idx,idx)=1.0_DP
              ENDDO

              !Grab interpolation parameters
              FIELD_VAR_TYPE=vectorEquations%vectorMapping%nonlinearMapping%residualVariables(1)%ptr%VARIABLE_TYPE
              DEPENDENT_INTERPOLATION_PARAMETERS=>equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr
              GEOMETRIC_INTERPOLATION_PARAMETERS=>equations%interpolation% &
                & geometricInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr
              IF(ASSOCIATED(FIBRE_FIELD)) THEN
                FIBRE_INTERPOLATION_PARAMETERS=>equations%interpolation%fibreInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr
              ENDIF
              IF(ASSOCIATED(MATERIALS_FIELD)) THEN
                MATERIALS_INTERPOLATION_PARAMETERS=>equations%interpolation%materialsInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr
              ENDIF
!              INDEPENDENT_INTERPOLATION_PARAMETERS=>equations%interpolation% &
!                & independentInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr

              CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ne, &
                & GEOMETRIC_INTERPOLATION_PARAMETERS,err,error,*999)
              IF(ASSOCIATED(FIBRE_FIELD)) THEN
                CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ne, &
                  & FIBRE_INTERPOLATION_PARAMETERS,err,error,*999)
              END IF
              IF(ASSOCIATED(MATERIALS_FIELD)) THEN
                CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ne, &
                  & MATERIALS_INTERPOLATION_PARAMETERS,err,error,*999)
              ENDIF
              CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ne, &
                & DEPENDENT_INTERPOLATION_PARAMETERS,err,error,*999)
!              CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ne, &
!                & INDEPENDENT_INTERPOLATION_PARAMETERS,err,error,*999)

              !Point interpolation pointer
              GEOMETRIC_INTERPOLATED_POINT=>equations%interpolation%geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr
              GEOMETRIC_INTERPOLATED_POINT_METRICS=>equations%interpolation% &
                & geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr
              IF(ASSOCIATED(FIBRE_FIELD)) THEN
                FIBRE_INTERPOLATED_POINT=>equations%interpolation%fibreInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr
              END IF
              IF(ASSOCIATED(MATERIALS_FIELD)) THEN
                MATERIALS_INTERPOLATED_POINT=>equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr
              ENDIF
              DEPENDENT_INTERPOLATED_POINT=>equations%interpolation%dependentInterpPoint(FIELD_VAR_TYPE)%ptr
              DEPENDENT_INTERPOLATED_POINT_METRICS=>equations%interpolation%dependentInterpPointMetrics(FIELD_VAR_TYPE)%ptr
!              INDEPENDENT_INTERPOLATED_POINT=>equations%interpolation%independentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr
              
              C=>MATERIALS_INTERPOLATED_POINT%VALUES(:,1)

              !Loop over gauss points
              DO gauss_idx=1,DEPENDENT_NUMBER_OF_GAUSS_POINTS

                !Interpolate dependent, geometric, fibre and material fields
                CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
                  & DEPENDENT_INTERPOLATED_POINT,err,error,*999)
                CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(DEPENDENT_BASIS%NUMBER_OF_XI, &
                  & DEPENDENT_INTERPOLATED_POINT_METRICS,err,error,*999)
                CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
                  & GEOMETRIC_INTERPOLATED_POINT,err,error,*999)
                CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI, &
                  & GEOMETRIC_INTERPOLATED_POINT_METRICS,err,error,*999)
                IF(ASSOCIATED(FIBRE_FIELD)) THEN
                  CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
                    & FIBRE_INTERPOLATED_POINT,err,error,*999)
                END IF
                CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
                  & MATERIALS_INTERPOLATED_POINT,err,error,*999)
!                CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
!                  & INDEPENDENT_INTERPOLATED_POINT,err,error,*999)

                !Calculate F=dZ/dNU, the deformation gradient tensor at the gauss point
                CALL FiniteElasticity_GaussDeformationGradientTensor(DEPENDENT_INTERPOLATED_POINT_METRICS, &
                  & GEOMETRIC_INTERPOLATED_POINT_METRICS,FIBRE_INTERPOLATED_POINT,DZDNU,ERR,ERROR,*999)

                !get A1, A2, x1, x2 at the Gauss point of the 3D finite elasticity element
                NULLIFY(FIELD_VARIABLE)
                CALL Field_VariableGet(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VARIABLE,err,error,*999)

                dof_idx=FIELD_VARIABLE%COMPONENTS(1)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gauss_idx,ne)
                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & dof_idx,A_1,err,error,*999)
                dof_idx=FIELD_VARIABLE%COMPONENTS(2)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gauss_idx,ne)
                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & dof_idx,A_2,err,error,*999)
                dof_idx=FIELD_VARIABLE%COMPONENTS(3)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gauss_idx,ne)
                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & dof_idx,x_1,err,error,*999)
                dof_idx=FIELD_VARIABLE%COMPONENTS(4)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gauss_idx,ne)
                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & dof_idx,x_2,err,error,*999)

                !--------------------------------------------------------------------------------------------
                SARCO_LENGTH=DZDNU(1,1)
                ! Calculate Filament-Overlap
                IF(SARCO_LENGTH.LE.0.635_DP) THEN
                  FORCE_LENGTH=0.0_DP
                ELSE IF(SARCO_LENGTH.LE.0.835_DP) THEN 
                  FORCE_LENGTH=4.2_DP*(SARCO_LENGTH-0.635_DP)
                ELSE IF(SARCO_LENGTH.LE.1.0_DP) THEN
                  FORCE_LENGTH=0.84_DP+0.9697_DP*(SARCO_LENGTH-0.835_DP)
                ELSE IF(SARCO_LENGTH.LE.1.125_DP) THEN
                  FORCE_LENGTH=1.0_DP
                ELSE IF(SARCO_LENGTH.LE.1.825_DP) THEN
                  FORCE_LENGTH=1.0_DP-1.4286_DP*(SARCO_LENGTH-1.125_DP)
                ELSE
                  FORCE_LENGTH=0.0_DP
                ENDIF

                REFERENCE_VOLUME=1.4965e-03_DP ! [micrometer^3]
                MAX_XB_NUMBER_PER_VOLUME=120.0_DP*2.0_DP/REFERENCE_VOLUME ! [cross-bridges per micrometer^3]
                ENERGY_PER_XB=0.5_DP*x_2**2*C(8) ! [Newton times micrometer]
      
                !Mechanical Energy stored in cross-bridges
                XB_ENERGY_PER_VOLUME=MAX_XB_NUMBER_PER_VOLUME*FORCE_LENGTH*ENERGY_PER_XB*A_2 ! [Newton per micrometer^2]
                !Mechanical Energy stored in cross-bridges converted in Newton per cm^2
                XB_ENERGY_PER_VOLUME=XB_ENERGY_PER_VOLUME*1e+08_DP ! [Newton per cm^2]

                !Initalize lambda_a
                lambda_a=1.0_DP
    
                F_a_inv=0.0_DP
                F_a_inv(1,1)=1.0_DP/lambda_a
                F_a_inv(2,2)=1.0_DP
                F_a_inv(3,3)=1.0_DP

                CALL MatrixProduct(DZDNU,F_a_inv,F_e,err,error,*999)
                CALL MatrixTranspose(F_e,F_e_T,err,error,*999)
                CALL MatrixProduct(F_e_T,F_e,C_e,err,error,*999)
    
                !Odgen law - 3 terms. Material Parameters C = [mu(1) mu(2) mu(3) alpha(1) alpha(2) alpha(3) mu_0]
                !CALL Eigenvalue(C_e,EVALUES,err,error,*999)
                CALL DSYEV('V','U',3,C_e,3,EVALUES,WORK,-1,ERR)
                IF(ERR.NE.0) CALL FlagError("Error in Eigenvalue computation",err,error,*999)
                LWORK=MIN(LWMAX,INT(WORK(1)))
                CALL DSYEV('V','U',3,C_e,3,EVALUES,WORK,LWORK,ERR)
                IF(ERR.NE.0) CALL FlagError("Error in Eigenvalue computation",err,error,*999)
                EVECTOR_1=C_e(:,1)
                EVECTOR_2=C_e(:,2)
                EVECTOR_3=C_e(:,3)

                DO i=1,3
                  DO j=1,3
                    EMATRIX_1(i,j)=EVECTOR_1(i)*EVECTOR_1(j)
                    EMATRIX_2(i,j)=EVECTOR_2(i)*EVECTOR_2(j)
                    EMATRIX_3(i,j)=EVECTOR_3(i)*EVECTOR_3(j)
                  END DO
                END DO

                CALL MatrixProduct(F_a_inv,EMATRIX_1,N1,err,error,*999)
                CALL MatrixProduct(N1,F_a_inv,N1,err,error,*999) ! F_a_inv=F_a_inv_T
                CALL MatrixProduct(F_a_inv,EMATRIX_2,N2,err,error,*999)
                CALL MatrixProduct(N2,F_a_inv,N2,err,error,*999) ! F_a_inv=F_a_inv_T
                CALL MatrixProduct(F_a_inv,EMATRIX_3,N3,err,error,*999)
                CALL MatrixProduct(N3,F_a_inv,N3,err,error,*999) ! F_a_inv=F_a_inv_T 

                FREE_ENERGY_0=0.0_DP
                DO i=1,3
                  FREE_ENERGY_0=FREE_ENERGY_0+C(i)/C(i+3)*( &
                    & EVALUES(1)**(C(i+3)/2.0_DP)+ &
                    & EVALUES(2)**(C(i+3)/2.0_DP)+ &
                    & EVALUES(3)**(C(i+3)/2.0_DP)-3.0_DP)
                END DO
                FREE_ENERGY_0=C(7)*FREE_ENERGY_0

                FREE_ENERGY=FREE_ENERGY_0

                VALUE=XB_ENERGY_PER_VOLUME-(FREE_ENERGY-FREE_ENERGY_0)

                !tolerance for Newton's method
                TOL=0.00001_DP
                !tolerance for the bisection method as preconditioner. Since Newton's method does not converge, we only use the bisection method here
                TOL1=TOL 
                UP=lambda_a
                LOW=0.001_DP

                DO WHILE (ABS(VALUE).GE.TOL)

                  !bisection method
                  IF (ABS(VALUE).GE.TOL1) THEN
                    lambda_a=UP-(UP-LOW)/2.0_DP

                    F_a_inv=0.0_DP
                    F_a_inv(1,1)=1.0_DP/lambda_a
                    F_a_inv(2,2)=1.0_DP
                    F_a_inv(3,3)=1.0_DP

                    CALL MatrixProduct(DZDNU,F_a_inv,F_e,err,error,*999)
                    CALL MatrixTranspose(F_e,F_e_T,err,error,*999)
                    CALL MatrixProduct(F_e_T,F_e,C_e,err,error,*999)

                    CALL DSYEV('V','U',3,C_e,3,EVALUES,WORK,-1,ERR)
                    IF(ERR.NE.0) CALL FlagError("Error in Eigenvalue computation",err,error,*999)
                    LWORK=MIN(LWMAX,INT(WORK(1)))
                    CALL DSYEV('V','U',3,C_e,3,EVALUES,WORK,LWORK,ERR)
                    IF(ERR.NE.0) CALL FlagError("Error in Eigenvalue computation",err,error,*999)
                    EVECTOR_1=C_e(:,1)
                    EVECTOR_2=C_e(:,2)
                    EVECTOR_3=C_e(:,3)

                    DO i=1,3
                      DO j=1,3
                        EMATRIX_1(i,j)=EVECTOR_1(i)*EVECTOR_1(j)
                        EMATRIX_2(i,j)=EVECTOR_2(i)*EVECTOR_2(j)
                        EMATRIX_3(i,j)=EVECTOR_3(i)*EVECTOR_3(j)
                      END DO
                    END DO

                    CALL MatrixProduct(F_a_inv,EMATRIX_1,N1,err,error,*999)
                    CALL MatrixProduct(N1,F_a_inv,N1,err,error,*999) ! F_a_inv=F_a_inv_T
                    CALL MatrixProduct(F_a_inv,EMATRIX_2,N2,err,error,*999)
                    CALL MatrixProduct(N2,F_a_inv,N2,err,error,*999) ! F_a_inv=F_a_inv_T
                    CALL MatrixProduct(F_a_inv,EMATRIX_3,N3,err,error,*999)
                    CALL MatrixProduct(N3,F_a_inv,N3,err,error,*999) ! F_a_inv=F_a_inv_T 

                    FREE_ENERGY=0.0_DP
                    DO i=1,3
                      FREE_ENERGY=FREE_ENERGY+C(i)/C(i+3)*( &
                        & EVALUES(1)**(C(i+3)/2.0_DP)+ &
                        & EVALUES(2)**(C(i+3)/2.0_DP)+ &
                        & EVALUES(3)**(C(i+3)/2.0_DP)-3.0_DP)
                    END DO
                    FREE_ENERGY=C(7)*FREE_ENERGY

                    VALUE=XB_ENERGY_PER_VOLUME-(FREE_ENERGY-FREE_ENERGY_0)

                    IF (VALUE.GE.0) THEN
                      UP=lambda_a
                    ELSE
                      LOW=lambda_a
                    ENDIF

                  ELSE 
                    !Newton's method -- needs to be checked TODO

                    TEMP=DZDNU+DZDNUT
                    CALL MatrixProduct(F_e_T,TEMP,TEMP,err,error,*999)
                    CALL MatrixProduct(TEMP,N1,TEMP1,err,error,*999) 
                    CALL MatrixProduct(TEMP,N2,TEMP2,err,error,*999) 
                    CALL MatrixProduct(TEMP,N3,TEMP3,err,error,*999) 

                    TEMP=0.0_DP
                    DO i=1,3
                      TEMP=TEMP+ &
                        & C(i)*EVALUES(1)**(C(i+3)/2.0_DP-1.0_DP)*TEMP1+ &
                        & C(i)*EVALUES(2)**(C(i+3)/2.0_DP-1.0_DP)*TEMP2+ &
                        & C(i)*EVALUES(3)**(C(i+3)/2.0_DP-1.0_DP)*TEMP3
                    END DO
                    
!                    SLOPE=TEMP(1,1)*C(7)
                    lambda_a=lambda_a-VALUE/SLOPE
                    !IF (lambda_a.LE.0.0_DP) THEN
                    ! lambda_a=0.1_DP
                    !END IF
                    !lambda_a=lambda_a-0.001

                    F_a_inv=0.0_DP
                    F_a_inv(1,1)=1.0_DP/lambda_a
                    F_a_inv(2,2)=1.0_DP
                    F_a_inv(3,3)=1.0_DP

                    CALL MatrixProduct(DZDNU,F_a_inv,F_e,err,error,*999)
                    CALL MatrixTranspose(F_e,F_e_T,err,error,*999)
                    CALL MatrixProduct(F_e_T,F_e,C_e,err,error,*999)

                    CALL DSYEV('V','U',3,C_e,3,EVALUES,WORK,-1,ERR)
                    IF(ERR.NE.0) CALL FlagError("Error in Eigenvalue computation",err,error,*999)
                    LWORK=MIN(LWMAX,INT(WORK(1)))
                    CALL DSYEV('V','U',3,C_e,3,EVALUES,WORK,LWORK,ERR)
                    IF(ERR.NE.0) CALL FlagError("Error in Eigenvalue computation",err,error,*999)
                    EVECTOR_1=C_e(:,1)
                    EVECTOR_2=C_e(:,2)
                    EVECTOR_3=C_e(:,3)

                    DO i=1,3
                      DO j=1,3
                        EMATRIX_1(i,j)=EVECTOR_1(i)*EVECTOR_1(j)
                        EMATRIX_2(i,j)=EVECTOR_2(i)*EVECTOR_2(j)
                        EMATRIX_3(i,j)=EVECTOR_3(i)*EVECTOR_3(j)
                      END DO
                    END DO

                    CALL MatrixProduct(F_a_inv,EMATRIX_1,N1,err,error,*999)
                    CALL MatrixProduct(N1,F_a_inv,N1,err,error,*999) ! F_a_inv=F_a_inv_T
                    CALL MatrixProduct(F_a_inv,EMATRIX_2,N2,err,error,*999)
                    CALL MatrixProduct(N2,F_a_inv,N2,err,error,*999) ! F_a_inv=F_a_inv_T
                    CALL MatrixProduct(F_a_inv,EMATRIX_3,N3,err,error,*999)
                    CALL MatrixProduct(N3,F_a_inv,N3,err,error,*999) ! F_a_inv=F_a_inv_T 

                    FREE_ENERGY=0.0_DP
                    DO i=1,3
                      FREE_ENERGY=FREE_ENERGY+C(i)/C(i+3)*( &
                        & EVALUES(1)**(C(i+3)/2.0_DP)+ &
                        & EVALUES(2)**(C(i+3)/2.0_DP)+ &
                        & EVALUES(3)**(C(i+3)/2.0_DP)-3.0_DP)
                    END DO
                    FREE_ENERGY=C(7)*FREE_ENERGY

                    VALUE=XB_ENERGY_PER_VOLUME-(FREE_ENERGY-FREE_ENERGY_0)
                  ENDIF
                ENDDO

                !store lambda_a at the Gauss point
                dof_idx=FIELD_VARIABLE%COMPONENTS(5)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gauss_idx,ne)
                CALL Field_ParameterSetUpdateLocalDOF(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & dof_idx,lambda_a,err,error,*999)

              ENDDO
            ENDDO
          ELSE
            LOCAL_ERROR="This routine is not implemented for the third equations set specification of "// &
              & TRIM(NumberToVString(EQUATIONS_SET%specification(3),"*",err,error))// &
              & " of a finite elasticity type of an elasticity equation set."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          ENDIF
        ENDDO
      ENDIF
    ENDIF

    EXITS("FINITE_ELASTICITY_EVALUATE_EVOLUTION_LAW")
    RETURN
999 ERRORS("FINITE_ELASTICITY_EVALUATE_EVOLUTION_LAW",err,error)
    EXITS("FINITE_ELASTICITY_EVALUATE_EVOLUTION_LAW")
    RETURN 1
  END SUBROUTINE FINITE_ELASTICITY_EVALUATE_EVOLUTION_LAW

  !   
  !================================================================================================================================
  !
  !>Read in the displacement field for a PGM elasticity problem
  SUBROUTINE FiniteElasticity_PreSolveGetSolidDisplacement(CONTROL_LOOP,SOLVER,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER_FINITE_ELASTICITY  !<A pointer to the solvers
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD_FINITE_ELASTICITY
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS_FINITE_ELASTICITY  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING_FINITE_ELASTICITY !<A pointer to the solver mapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET_FINITE_ELASTICITY !<A pointer to the equations set
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_TIME_LOOP

    REAL(DP), POINTER :: MESH_DISPLACEMENT_VALUES(:)
    REAL(DP), POINTER :: DUMMY_VALUES2(:)
    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT

    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS,NDOFS_TO_PRINT
    INTEGER(INTG) :: INPUT_TYPE,INPUT_OPTION
    INTEGER(INTG) :: loop_idx

    ENTERS("FiniteElasticity_PreSolveGetSolidDisplacement",err,error,*999)

!--- \todo : Do we need for each case a Field_ParameterSetUpdateStart / FINISH on FIELD_MESH_DISPLACEMENT_SET_TYPE ?

    NULLIFY(SOLVER_FINITE_ELASTICITY)
    NULLIFY(MESH_DISPLACEMENT_VALUES)
    NULLIFY(DUMMY_VALUES2)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      CONTROL_TIME_LOOP=>CONTROL_LOOP
      DO loop_idx=1,CONTROL_LOOP%CONTROL_LOOP_LEVEL
        IF(CONTROL_TIME_LOOP%LOOP_TYPE==PROBLEM_CONTROL_TIME_LOOP_TYPE) THEN
          CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_TIME_LOOP,CURRENT_TIME,TIME_INCREMENT,err,error,*999)
          EXIT
        ENDIF
        IF (ASSOCIATED(CONTROL_LOOP%PARENT_LOOP)) THEN
          CONTROL_TIME_LOOP=>CONTROL_TIME_LOOP%PARENT_LOOP
        ELSE
          CALL FlagError("Could not find a time control loop.",err,error,*999)
        ENDIF
      ENDDO

      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a finite elasticity problem.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE, &
              & PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE)
              !--- Motion: read in from a file
              IF(SOLVER%GLOBAL_NUMBER==1) THEN
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,1,SOLVER_FINITE_ELASTICITY,err,error,*999)
                SOLVER_EQUATIONS_FINITE_ELASTICITY=>SOLVER_FINITE_ELASTICITY%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS_FINITE_ELASTICITY)) THEN
                  SOLVER_MAPPING_FINITE_ELASTICITY=>SOLVER_EQUATIONS_FINITE_ELASTICITY%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING_FINITE_ELASTICITY)) THEN
                    EQUATIONS_SET_FINITE_ELASTICITY=>SOLVER_MAPPING_FINITE_ELASTICITY%EQUATIONS_SETS(1)%ptr
                    IF(ASSOCIATED(EQUATIONS_SET_FINITE_ELASTICITY)) THEN
                      DEPENDENT_FIELD_FINITE_ELASTICITY=>EQUATIONS_SET_FINITE_ELASTICITY%DEPENDENT%DEPENDENT_FIELD
                    ELSE
                      CALL FlagError("Finite elasticity equations set is not associated.",err,error,*999)
                    END IF
                    CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Finite Elasticity motion read from a file ... ",err,error,*999)

                    CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET_FINITE_ELASTICITY%GEOMETRY%GEOMETRIC_FIELD, & 
                      & FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS,err,error,*999)

                    !Copy input to Finite elasticity's dependent field
                    !\todo: Still need to take into account that we are reading in displacement,
                    !       while dependent field is the absolute position of the structure
                    INPUT_TYPE=42
                    INPUT_OPTION=2
                    CALL Field_ParameterSetDataGet(EQUATIONS_SET_FINITE_ELASTICITY%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,MESH_DISPLACEMENT_VALUES,err,error,*999)
                    CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,MESH_DISPLACEMENT_VALUES, & 
                      & NUMBER_OF_DIMENSIONS,INPUT_TYPE,INPUT_OPTION,CONTROL_LOOP%TIME_LOOP%ITERATION_NUMBER,1.0_DP, &
                      & err,error,*999)
                    CALL Field_ParameterSetUpdateStart(EQUATIONS_SET_FINITE_ELASTICITY%DEPENDENT%DEPENDENT_FIELD, & 
                      & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
                    CALL Field_ParameterSetUpdateFinish(EQUATIONS_SET_FINITE_ELASTICITY%DEPENDENT%DEPENDENT_FIELD, & 
                      & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
                  ELSE
                    CALL FlagError("Finite elasticity solver mapping is not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("Finite elasticity solver equations are not associated.",err,error,*999)
                END IF

               IF(DIAGNOSTICS1) THEN
                 NDOFS_TO_PRINT = SIZE(MESH_DISPLACEMENT_VALUES,1)
                 CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,NDOFS_TO_PRINT,NDOFS_TO_PRINT,NDOFS_TO_PRINT,&
                   & MESH_DISPLACEMENT_VALUES,'(" MESH_DISPLACEMENT_VALUES = ",4(X,E13.6))','4(4(X,E13.6))', &
                   & err,error,*999)
               ENDIF
              ELSE
                ! in case of a solver number different from 3: do nothing ???
              ENDIF
            CASE DEFAULT
              LOCAL_ERROR="Problem subtype "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
                & " is not valid for a Finite elasticity equation fluid type of a fluid mechanics problem class."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",err,error,*999)
    ENDIF

    EXITS("FiniteElasticity_PreSolveGetSolidDisplacement")
    RETURN
999 ERRORS("FiniteElasticity_PreSolveGetSolidDisplacement",err,error)
    EXITS("FiniteElasticity_PreSolveGetSolidDisplacement")
    RETURN 1

  END SUBROUTINE FiniteElasticity_PreSolveGetSolidDisplacement

  !
  !================================================================================================================================
  !

  !>Update boundary conditions for finite elasticity pre solve
  SUBROUTINE FiniteElasticity_PreSolveUpdateBoundaryConditions(solver,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: boundaryConditionsVariable
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: boundaryConditions
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop!<A pointer to the control loop to solve.
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    TYPE(FIELD_TYPE), POINTER :: dependentField, geometricField
    TYPE(PROBLEM_TYPE), POINTER :: problem
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping !<A pointer to the solver mapping
 
    REAL(DP) :: currentTime,timeIncrement,alpha
    REAL(DP), POINTER :: GEOMETRIC_FIELD_VALUES(:) 
    REAL(DP), POINTER :: MESH_POSITION_VALUES(:) 
    REAL(DP), POINTER :: DUMMY_VALUES1(:), CURRENT_PRESSURE_VALUES(:)
    REAL(DP), ALLOCATABLE :: NEW_PRESSURE_VALUES(:)

    INTEGER(INTG) :: BOUNDARY_CONDITION_CHECK_VARIABLE
    INTEGER(INTG) :: dof_number,GEOMETRY_NUMBER_OF_DOFS,DEPENDENT_NUMBER_OF_DOFS
    INTEGER(INTG) :: NDOFS_TO_PRINT
    INTEGER(INTG) :: loop_idx
    INTEGER(INTG) :: SUBITERATION_NUMBER

    ENTERS("FiniteElasticity_PreSolveUpdateBoundaryConditions",err,error,*999)

    NULLIFY( CURRENT_PRESSURE_VALUES, DUMMY_VALUES1 )

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
    
    IF(solver%GLOBAL_NUMBER==1) THEN
      NULLIFY(solverEquations)
      CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      NULLIFY(controlLoop)
      CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
      CALL ControlLoop_CurrentTimesGet(controlLoop,currentTime,timeIncrement,err,error,*999)
      NULLIFY(problem)
      CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
      IF(.NOT.ALLOCATED(problem%specification)) CALL FlagError("Problem specification is not allocated.",err,error,*999)
      IF(SIZE(problem%specification,1)<3) &
        & CALL FlagError("Problem specification must have three entries for a finite elasticity problem.",err,error,*999)
      SELECT CASE(problem%specification(3))
      CASE(PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
        equations=>solverMapping%EQUATIONS_SET_TO_SOLVER_MAP(1)%equations
        IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations are not associated.",err,error,*999)
        NULLIFY(equationsSet)
        CALL Equations_EquationsSetGet(equations,equationsSet,err,error,*999)
        IF(.NOT.ALLOCATED(equationsSet%specification)) &
          & CALL FlagError("Equations set specification is not allocated.",err,error,*999)
        IF(SIZE(equationsSet%specification,1)/=3) &
          & CALL FlagError("Equations set specification must have three entries for a finite elasticity type equations set.", &
          & err,error,*999)
        SELECT CASE(equationsSet%specification(3))
        CASE(EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE)
          NULLIFY(dependentField)
          CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
          NULLIFY(boundaryConditions)
          CALL SolverEquations_BoundaryConditionsGet(solverEquations,boundaryConditions,err,error,*999)
          NULLIFY(vectorEquations)
          CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
          NULLIFY(vectorMapping)
          CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
          CALL Field_VariableGet(dependentField,FIELD_DELUDELN_VARIABLE_TYPE,fieldVariable,err,error,*999)
          CALL BOUNDARY_CONDITIONS_VARIABLE_GET(boundaryConditions,fieldVariable,boundaryConditionsVariable,err,error,*999)
          IF(ASSOCIATED(boundaryConditionsVariable)) THEN
            IF(boundaryConditionsVariable%DOF_COUNTS(BOUNDARY_CONDITION_PRESSURE)>0) THEN
              
              IF(DIAGNOSTICS1) THEN
                CALL Field_ParameterSetDataGet(dependentField,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_PRESSURE_VALUES_SET_TYPE, &
                  & CURRENT_PRESSURE_VALUES,err,error,*999)
                NDOFS_TO_PRINT = SIZE(CURRENT_PRESSURE_VALUES,1)
                CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,NDOFS_TO_PRINT,NDOFS_TO_PRINT, &
                  & NDOFS_TO_PRINT,CURRENT_PRESSURE_VALUES, &
                  & '(" DEP_FIELD,FIELD_U_VAR_TYPE,FIELD_PRESSURE_VAL_SET_TYPE (before) = ",4(X,E13.6))',&
                  & '4(4(X,E13.6))',err,error,*999)
                CALL Field_ParameterSetDataRestore(dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_PRESSURE_VALUES_SET_TYPE,CURRENT_PRESSURE_VALUES,err,error,*999)
              ENDIF
              
              DEPENDENT_NUMBER_OF_DOFS=dependentField%VARIABLE_TYPE_MAP(FIELD_DELUDELN_VARIABLE_TYPE)%PTR%NUMBER_OF_DOFS

              ALLOCATE(NEW_PRESSURE_VALUES(DEPENDENT_NUMBER_OF_DOFS))
              
              !Linear increase of cavity pressure: just a test example prototype
              !\todo: general time-dependent boundary condition input method?
              ALPHA = ( currentTime + timeIncrement ) / currentTime
              NEW_PRESSURE_VALUES(1:DEPENDENT_NUMBER_OF_DOFS) = ALPHA * CURRENT_PRESSURE_VALUES(1:DEPENDENT_NUMBER_OF_DOFS)
              
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Finite Elasticity update pressure BCs",err,error,*999)
              DO dof_number=1,DEPENDENT_NUMBER_OF_DOFS
                CALL Field_ParameterSetUpdateLocalDOF(dependentField,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_PRESSURE_VALUES_SET_TYPE, &
                  & dof_number,NEW_PRESSURE_VALUES(dof_number), err,error,*999)
              ENDDO
              CALL Field_ParameterSetUpdateStart(dependentField,FIELD_DELUDELN_VARIABLE_TYPE, FIELD_PRESSURE_VALUES_SET_TYPE, &
                & err,error,*999)
              CALL Field_ParameterSetUpdateFinish(dependentField,FIELD_DELUDELN_VARIABLE_TYPE, FIELD_PRESSURE_VALUES_SET_TYPE, &
                & err,error,*999)
              
              DEALLOCATE(NEW_PRESSURE_VALUES)
              
              IF(DIAGNOSTICS1) THEN
                NULLIFY( DUMMY_VALUES1 )
                CALL Field_ParameterSetDataGet(dependentField,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_PRESSURE_VALUES_SET_TYPE, &
                  & DUMMY_VALUES1,err,error,*999)
                NDOFS_TO_PRINT = SIZE(DUMMY_VALUES1,1)
                CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,NDOFS_TO_PRINT,NDOFS_TO_PRINT, &
                  & NDOFS_TO_PRINT,DUMMY_VALUES1, &
                  & '(" DEP_FIELD,FIELD_U_VAR_TYPE,FIELD_PRESSURE_VAL_SET_TYPE (after) = ",4(X,E13.6))', &
                  & '4(4(X,E13.6))',err,error,*999)
                CALL Field_ParameterSetDataRestore(dependentField,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_PRESSURE_VALUES_SET_TYPE, &
                  & DUMMY_VALUES1,err,error,*999)
              ENDIF
              CALL Field_ParameterSetDataRestore(dependentField,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_PRESSURE_VALUES_SET_TYPE, &
                & CURRENT_PRESSURE_VALUES,err,error,*999)
            ENDIF !Pressure_condition_used
          ELSE
            CALL FlagError("Boundary condition variable is not associated.",err,error,*999)
          END IF
        CASE DEFAULT
          ! do nothing 
        END SELECT
        
      CASE(PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE)
        equations=>solverMapping%EQUATIONS_SET_TO_SOLVER_MAP(1)%equations
        IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations are not associated.",err,error,*999)
        NULLIFY(equationsSet)
        CALL Equations_EquationsSetGet(equations,equationsSet,err,error,*999)
        IF(.NOT.ALLOCATED(equationsSet%specification)) &
          & CALL FlagError("Equations set specification is not allocated.",err,error,*999)
        IF(SIZE(equationsSet%specification,1)/=3) &
          & CALL FlagError("Equations set specification must have three entries for a finite elasticity type equations set.", &
          & err,error,*999)
        SELECT CASE(equationsSet%specification(3))
        CASE(EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
          & EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE)
          IF(solver%outputType>=SOLVER_PROGRESS_OUTPUT) THEN
            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Finite Elasticity update BCs",err,error,*999)
          ENDIF
          NULLIFY(dependentField)
          CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
          NULLIFY(geometricField)
          CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
          NULLIFY(boundaryConditions)
          CALL SolverEquations_BoundaryConditionsGet(solverEquations,boundaryConditions,err,error,*999)
          NULLIFY(vectorEquations)
          CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
          NULLIFY(vectorMapping)
          CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
          CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,fieldVariable,err,error,*999)
          CALL BOUNDARY_CONDITIONS_VARIABLE_GET(boundaryConditions,fieldVariable,boundaryConditionsVariable,err,error,*999)
          IF(ASSOCIATED(boundaryConditionsVariable)) THEN
            IF(DIAGNOSTICS1) THEN
              NULLIFY( DUMMY_VALUES1 )
              CALL Field_ParameterSetDataGet(dependentField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,DUMMY_VALUES1,err,error,*999)
              NDOFS_TO_PRINT = SIZE(DUMMY_VALUES1,1)
              CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,NDOFS_TO_PRINT,NDOFS_TO_PRINT, &
                & NDOFS_TO_PRINT,DUMMY_VALUES1, &
                & '(" dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE (bef) = ",4(X,E13.6))',&
                & '4(4(X,E13.6))',err,error,*999)
              CALL Field_ParameterSetDataRestore(dependentField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,DUMMY_VALUES1,err,error,*999)
            ENDIF
            
            ! requires solid dependent field and geometry to be interpolated identically !!!
            ! assumes that DOFs for dependent and geometric field are stored in the same order
            ! How does this routine take into account the BC value ???
            ALPHA = 0.10_DP * SIN( 2.0_DP * PI * currentTime / 4.0_DP )
            CALL Field_ParameterSetsCopy(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & FIELD_MESH_DISPLACEMENT_SET_TYPE,ALPHA,err,error,*999)
            
            NULLIFY(GEOMETRIC_FIELD_VALUES)
            CALL Field_ParameterSetDataGet(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,GEOMETRIC_FIELD_VALUES, &
              & err,error,*999)
            
            GEOMETRY_NUMBER_OF_DOFS=geometricField%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR%NUMBER_OF_DOFS
            DO dof_number=1,GEOMETRY_NUMBER_OF_DOFS
              BOUNDARY_CONDITION_CHECK_VARIABLE=boundaryConditionsVariable%CONDITION_TYPES(dof_number)
              IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_MOVED_WALL .OR. &
                & BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_MOVED_WALL_INCREMENTED) THEN
                !--- To obtain absolute positions, add nodal coordinates on top of mesh displacement
                CALL Field_ParameterSetAddLocalDOF(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE, &
                  & dof_number,GEOMETRIC_FIELD_VALUES(dof_number),err,error,*999)
              ELSE
                ! do nothing ???
              END IF
            END DO
            
            NULLIFY(MESH_POSITION_VALUES)
            CALL Field_ParameterSetDataGet(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE, &
              & MESH_POSITION_VALUES,err,error,*999)
            
            DEPENDENT_NUMBER_OF_DOFS=dependentField%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR%NUMBER_OF_DOFS
            DO dof_number=1,DEPENDENT_NUMBER_OF_DOFS
              BOUNDARY_CONDITION_CHECK_VARIABLE=boundaryConditionsVariable%CONDITION_TYPES(dof_number)
              IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_MOVED_WALL .OR. &
                & BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_MOVED_WALL_INCREMENTED) THEN
                
                !Update FIELD_BOUNDARY_CONDITIONS_SET_TYPE or FIELD_VALUES_SET_TYPE
                !(so it is one or the other, but not both) depending on whether or not load increments are used
                IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_MOVED_WALL_INCREMENTED) THEN
                  CALL Field_ParameterSetUpdateLocalDOF(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_BOUNDARY_CONDITIONS_SET_TYPE, &
                    & dof_number,MESH_POSITION_VALUES(dof_number),err,error,*999)
                ELSE
                  !--- Update the dependent field with the new absolute position
                  CALL Field_ParameterSetUpdateLocalDOF(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,dof_number, &
                    & MESH_POSITION_VALUES(dof_number),err,error,*999)
                ENDIF
                
              ELSE
                ! do nothing ???
              END IF
            END DO
            
            IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_MOVED_WALL_INCREMENTED) THEN
              CALL Field_ParameterSetUpdateStart(dependentField,FIELD_U_VARIABLE_TYPE, FIELD_BOUNDARY_CONDITIONS_SET_TYPE, &
                & err,error,*999)
              CALL Field_ParameterSetUpdateFinish(dependentField,FIELD_U_VARIABLE_TYPE, FIELD_BOUNDARY_CONDITIONS_SET_TYPE, &
                & err,error,*999)
            ELSE
              CALL Field_ParameterSetUpdateStart(dependentField,FIELD_U_VARIABLE_TYPE, FIELD_VALUES_SET_TYPE,err,error,*999)
              CALL Field_ParameterSetUpdateFinish(dependentField,FIELD_U_VARIABLE_TYPE, FIELD_VALUES_SET_TYPE,err,error,*999)
            ENDIF
            
            IF(DIAGNOSTICS1) THEN
              NULLIFY( DUMMY_VALUES1 )
              CALL Field_ParameterSetDataGet(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,DUMMY_VALUES1, &
                & err,error,*999)
              NDOFS_TO_PRINT = SIZE(DUMMY_VALUES1,1)
              CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,NDOFS_TO_PRINT,NDOFS_TO_PRINT, &
                & NDOFS_TO_PRINT,DUMMY_VALUES1, &
                & '(" dependentField,FIELD_U_VAR_TYPE,FIELD_VALUES_SET_TYPE (after) = ",4(X,E13.6))', &
                & '4(4(X,E13.6))',err,error,*999)
              CALL Field_ParameterSetDataRestore(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,DUMMY_VALUES1, &
                & err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Boundary condition variable is not associated.",err,error,*999)
          END IF
          CALL Field_ParameterSetUpdateStart(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
          CALL Field_ParameterSetUpdateFinish(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
        CASE DEFAULT
          ! do nothing
        END SELECT
      CASE DEFAULT
        ! do nothing 
      END SELECT
    ENDIF
    
    EXITS("FiniteElasticity_PreSolveUpdateBoundaryConditions")
    RETURN
999 ERRORS("FiniteElasticity_PreSolveUpdateBoundaryConditions",err,error)
    EXITS("FiniteElasticity_PreSolveUpdateBoundaryConditions")
    RETURN 1

  END SUBROUTINE FiniteElasticity_PreSolveUpdateBoundaryConditions

  !
  !================================================================================================================================
  !

  !>Evaluates the functions f(J) and f\'(J);
  !>  Eq.(21) in Chapelle, Gerbeau, Sainte-Marie, Vignon-Clementel, Computational Mechanics (2010)
  SUBROUTINE EVALUATE_CHAPELLE_FUNCTION(Jznu,ffact,dfdJfact,err,error,*)
  
    !Argument variables
    REAL(DP), INTENT(IN) :: Jznu !<Jznu=DETERMINANT(AZL,err,error)**0.5_DP
    REAL(DP), INTENT(OUT) :: ffact !<f(Jznu) of the INRIA model
    REAL(DP), INTENT(OUT) :: dfdJfact !<dfdJfact = f'(Jznu) of the INRIA model
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    ENTERS("EVALUATE_CHAPELLE_FUNCTION",err,error,*999)

!     IF( ABS(Jznu-1.0_DP) > 5.0E-02_DP ) THEN
    IF( ABS(Jznu-1.0_DP) > 1.0E-10_DP ) THEN
      !Eq.(21) of the INRIA paper
      ffact = 2.0_DP * (Jznu - 1.0_DP - log(Jznu)) / (Jznu - 1.0_DP)**2.0_DP
      dfdJfact = ( 2.0_DP * (1.0_DP - 1.0_DP/Jznu) * (Jznu - 1.0_DP)**2.0_DP &
        & - 4.0_DP * (Jznu - 1.0_DP - log(Jznu)) * (Jznu - 1.0_DP) ) / (Jznu - 1.0_DP)**4.0_DP
    ELSE
      ffact = 1.0_DP
      dfdJfact = 0.0_DP
    END IF

    EXITS("EVALUATE_CHAPELLE_FUNCTION")
    RETURN
999 ERRORSEXITS("EVALUATE_CHAPELLE_FUNCTION",err,error)
    RETURN 1
  END SUBROUTINE EVALUATE_CHAPELLE_FUNCTION

  !
  !================================================================================================================================
  !

  !>Evaluates the 2nd Piola-Kirchhoff stress tensor;
  !>  Eq.(13) in Chapelle, Gerbeau, Sainte-Marie, Vignon-Clementel, Computational Mechanics (2010)
  SUBROUTINE EVALUATE_CHAPELLE_PIOLA_TENSOR_ADDITION(AZL,AZU,DARCY_MASS_INCREASE,PIOLA_TENSOR_ADDITION,err,error,*)
  
    !Argument variables
    REAL(DP), INTENT(IN) :: AZL(3,3) !<C=F\'F
    REAL(DP), INTENT(IN) :: AZU(3,3) !<inverse of AZL
    REAL(DP), INTENT(IN) :: DARCY_MASS_INCREASE !<mass increase
    REAL(DP), INTENT(OUT) :: PIOLA_TENSOR_ADDITION(3,3) !<Addition to the 2nd Piola-Kirchhoff tensor
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    REAL(DP) :: Jznu !<Jznu=DETERMINANT(AZL,err,error)**0.5_DP
    REAL(DP) :: ffact !<f(Jznu) of the INRIA model
    REAL(DP) :: dfdJfact !<dfdJfact = f\'(Jznu) of the INRIA model
    REAL(DP) :: Mfact, bfact, p0fact  !<INRIA constitutive law
    REAL(DP) :: DARCY_VOL_INCREASE, DARCY_RHO_0_F
    INTEGER(INTG) :: i,j

    
    ENTERS("EVALUATE_CHAPELLE_PIOLA_TENSOR_ADDITION",err,error,*999)

    !Parameters settings for coupled elasticity Darcy INRIA model:
    CALL GET_DARCY_FINITE_ELASTICITY_PARAMETERS(DARCY_RHO_0_F,Mfact,bfact,p0fact,err,error,*999)

    DARCY_VOL_INCREASE = DARCY_MASS_INCREASE / DARCY_RHO_0_F

    CALL Determinant(AZL,Jznu,err,error,*999)
    Jznu=Jznu**0.5_DP
    IF( ABS(Jznu) < 1.0E-10_DP ) THEN
      CALL FlagError("EVALUATE_CHAPELLE_PIOLA_TENSOR_ADDITION: ABS(Jznu) < 1.0E-10_DP",err,error,*999)
    END IF

    CALL EVALUATE_CHAPELLE_FUNCTION(Jznu,ffact,dfdJfact,err,error,*999)

    DO i=1,3
      DO j=1,3
!         PIOLA_TENSOR_ADDITION(i,j) = - Mfact * bfact * DARCY_VOL_INCREASE * (ffact + (Jznu - 1.0_DP) * dfdJfact) * Jznu * AZU(i,j) &
!           & + 0.5_DP * Mfact * DARCY_VOL_INCREASE**2.0_DP * dfdJfact * Jznu * AZU(i,j)
        PIOLA_TENSOR_ADDITION(i,j) = 0.5_DP * Mfact * DARCY_VOL_INCREASE**2.0_DP * Jznu * AZU(i,j)
!         PIOLA_TENSOR_ADDITION(i,j) = 0.0_DP
      ENDDO
    ENDDO

!     PIOLA_TENSOR_ADDITION = - Mfact * bfact * DARCY_VOL_INCREASE * (ffact + (Jznu - 1.0_DP) * dfdJfact) * Jznu * AZU &
!       & + 0.5_DP * Mfact * DARCY_VOL_INCREASE**2.0_DP * dfdJfact * Jznu * AZU

    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  DARCY_VOL_INCREASE = ",DARCY_VOL_INCREASE,err,error,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Jznu = ",Jznu,err,error,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  ffact = ",ffact,err,error,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  dfdJfact = ",dfdJfact,err,error,*999)
      CALL WRITE_STRING_MATRIX(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
        & 3,3,AZU,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("    AZU','(",I1,",:)',' :",3(X,E13.6))', &
        & '(17X,3(X,E13.6))',err,error,*999)
      CALL WRITE_STRING_MATRIX(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
        & 3,3,PIOLA_TENSOR_ADDITION, &
        & WRITE_STRING_MATRIX_NAME_AND_INDICES,'("    PIOLA_TENSOR_ADDITION','(",I1,",:)',' :",3(X,E13.6))', &
        & '(17X,3(X,E13.6))',err,error,*999)
    ENDIF

    EXITS("EVALUATE_CHAPELLE_PIOLA_TENSOR_ADDITION")
    RETURN
999 ERRORSEXITS("EVALUATE_CHAPELLE_PIOLA_TENSOR_ADDITION",err,error)
    RETURN 1
  END SUBROUTINE EVALUATE_CHAPELLE_PIOLA_TENSOR_ADDITION

  !
  !================================================================================================================================
  !

  !>Sets some data for the coupled Darcy / finite-elasticity model
  SUBROUTINE GET_DARCY_FINITE_ELASTICITY_PARAMETERS(DARCY_RHO_0_F,Mfact,bfact,p0fact,err,error,*)
  
    !Argument variables
    REAL(DP), INTENT(OUT) :: DARCY_RHO_0_F
    REAL(DP), INTENT(OUT) :: Mfact
    REAL(DP), INTENT(OUT) :: bfact
    REAL(DP), INTENT(OUT) :: p0fact
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    ENTERS("GET_DARCY_FINITE_ELASTICITY_PARAMETERS",err,error,*999)

!   DARCY_RHO_0_F = 1.0E-03_DP
    DARCY_RHO_0_F = 1.0_DP
!   Mfact = 2.18E05_DP
    Mfact = 2.18E00_DP
    bfact = 1.0_DP
    p0fact = 0.0_DP

    EXITS("GET_DARCY_FINITE_ELASTICITY_PARAMETERS")
    RETURN
999 ERRORSEXITS("GET_DARCY_FINITE_ELASTICITY_PARAMETERS",err,error)
    RETURN 1
  END SUBROUTINE GET_DARCY_FINITE_ELASTICITY_PARAMETERS

  !
  !================================================================================================================================
  !

  !> Apply load increments to the gravity vector
  SUBROUTINE FINITE_ELASTICITY_LOAD_INCREMENT_APPLY(EQUATIONS_SET,ITERATION_NUMBER,MAXIMUM_NUMBER_OF_ITERATIONS,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    INTEGER(INTG), INTENT(IN) :: ITERATION_NUMBER !<The current load increment iteration index
    INTEGER(INTG), INTENT(IN) :: MAXIMUM_NUMBER_OF_ITERATIONS !<Final index for load increment loop
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local variables
    TYPE(EquationsType), POINTER :: EQUATIONS
    TYPE(FIELD_TYPE), POINTER :: SOURCE_FIELD
    REAL(DP) :: INCREMENT

    ENTERS("FINITE_ELASTICITY_LOAD_INCREMENT_APPLY",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        SOURCE_FIELD=>equations%interpolation%sourceField
        IF(ASSOCIATED(SOURCE_FIELD)) THEN
          IF(MAXIMUM_NUMBER_OF_ITERATIONS>1) THEN
            IF(ITERATION_NUMBER==1) THEN
              !Setup initial values parameter set
              CALL Field_ParameterSetEnsureCreated(SOURCE_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_INITIAL_VALUES_SET_TYPE,err,error,*999)
              CALL Field_ParameterSetsCopy(SOURCE_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & FIELD_INITIAL_VALUES_SET_TYPE,1.0_DP,err,error,*999)
            ENDIF
            INCREMENT=REAL(ITERATION_NUMBER)/REAL(MAXIMUM_NUMBER_OF_ITERATIONS)
            CALL Field_ParameterSetsCopy(SOURCE_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_INITIAL_VALUES_SET_TYPE, &
                & FIELD_VALUES_SET_TYPE,INCREMENT,err,error,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FlagError("Equations set equations is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

    EXITS("FINITE_ELASTICITY_LOAD_INCREMENT_APPLY")
    RETURN
999 ERRORSEXITS("FINITE_ELASTICITY_LOAD_INCREMENT_APPLY",err,error)
    RETURN 1

  END SUBROUTINE FINITE_ELASTICITY_LOAD_INCREMENT_APPLY

  !
  !================================================================================================================================
  !
  
!   Main functions to compute everything
      
  SUBROUTINE clooping(hstep,tol,maxitr,mu0,kc,kt,k1,k2,k3,k12,ka, &
    & q,H1,H2,H12,Jh,Gamma,Gammam,Dt,Fr,Jr,Jen,Beprn,s1n,s2n,s3n, &
    & alpha1n,Bn,Je,s1,s2,s3,devH,HH,alpha1,B,Bepr,lame1,lame2,lame3, &
    & lamea,Ee11,Ee22,Ee33,Ee12,Eea,QQ,T0,T1,T2,T3,T4,T5,TT)
      
    IMPLICIT NONE
      
!   Declarations
    REAL*8 hstep,tol,maxitr,alpha1n,Bn
    REAL*8 mu0,kc,kt,k1,k2,k3,k12,ka,q,H1,H2,H12,Jh,Gamma,Gammam,Dt
    REAL*8 Fr(3,3),Jr,Jen,Beprn(3,3),s1n(3),s2n(3),s3n(3),Je,s1(3)
    REAL*8 s2(3),s3(3),S11(3,3),S22(3,3),S33(3,3),S12(3,3),S21(3,3)
    REAL*8 S13(3,3),S31(3,3),S23(3,3),S32(3,3),devH(3,3),HH(3,3)
    REAL*8 alpha1,B,Bepr(3,3),lame1,lame2,lame3,lamea,Ee11,Ee22,Ee33
    REAL*8 Ee12,Eea,devBepr(3,3),QQ,T0(3,3),T1(3,3),T2(3,3),T3(3,3)
    REAL*8 T4(3,3),T5(3,3),TT(3,3)
    INTEGER i1,i2
    
    CALL ElasticDil(Jr,Jh,Dt,Gammam,Jen,Je)
    
    CALL StateVarsIntegrator(hstep,tol,maxitr,Fr,Dt,Gamma,Beprn, &
      & alpha1n,Bn,s1n,s2n,s3n,H1,H2,H12,s1,s2,s3,S11,S22,S33,S12,S21, &
      & devH,HH,B,alpha1,devBepr,Bepr)
      
    CALL ElasticMeasures(Je,Bepr,S11,S22,S33,S12,lame1,lame2,lame3, &
      & lamea,Ee11,Ee22,Ee33,Ee12,Eea)
     
    CALL Constitutive(mu0,kc,kt,k1,k2,k3,k12,ka,q,Je,devBepr,Bepr, &
      & alpha1,S11,S22,S33,S12,S21,lame1,lame2,lame3,lamea,Ee12,Ee11, &
      & Ee22,Ee33,Eea,QQ,T0,T1,T2,T3,T4,T5,TT)
                      
  END SUBROUTINE clooping
      
!   ====================================================================
      
!   Elastic dilatation at t=t2
      
  SUBROUTINE ElasticDil(Jr,Jh,Dt,Gammam,Jen,Je) 
      
    IMPLICIT NONE
      
!   Declarations
    REAL*8 Jr,Jh,Dt,Gammam,Jen,Jn,Je,Jes
    
    Jes=Jr*Jen
    Je=dexp((dlog(Jes)+Dt*Gammam*dlog(Jh))/(1+Dt*Gammam))
    
  END SUBROUTINE ElasticDil
      
!   ====================================================================
      
  SUBROUTINE StateVarsIntegrator(hstep,tol,maxitr,Fr,Dt,Gamma, &
    & Beprn,alpha1n,Bn,s1n,s2n,s3n,H1,H2,H12,s1,s2,s3,S11,S22,S33,S12, &
    & S21,devH,HH,B,alpha1,devBepr,Bepr)
    
    IMPLICIT NONE
      
!   Declarations
    REAL*8 hstep,tol,maxitr,Fr(3,3),Beprn(3,3),Dt,Gamma,alpha1n,Bn
    REAL*8 s1n(3),s2n(3),s3n(3),H1,H2,H12,alpha1,B,Bepr(3,3),s1(3)
    REAL*8 s2(3),s3(3),S11(3,3),S22(3,3),S33(3,3),S12(3,3),S21(3,3)
    REAL*8 S13(3,3),S31(3,3),S23(3,3),S32(3,3),devH(3,3),devBepr(3,3)
    REAL*8 HH(3,3),BB(3,3),Frt(3,3),aux(3,3),EyeTens(3,3)
    INTEGER i1,i2
      
    CALL IdentityTens(EyeTens)
    CALL UpdateFibersHomeo(s1n,s2n,s3n,Fr,H1,H2,H12,s1,s2,s3,S11, &
      & S22,S33,S12,S21,S13,S31,S23,S32,devH,HH)     
    CALL alpha1B(hstep,tol,maxitr,Fr,Beprn,Dt,Gamma,devH,HH,alpha1n, &
      & Bn,alpha1,B)       
    CALL deviatorBepr(Fr,Beprn,Dt,Gamma,B,devH,devBepr)
    DO i1=1,3
      DO i2=1,3
        Bepr(i1,i2)=(alpha1/3.0d0)*EyeTens(i1,i2)+devBepr(i1,i2)
      END DO
    END DO
    
  END SUBROUTINE StateVarsIntegrator
  
!   ====================================================================
      
!   Compute the deviatoric part Be''

  SUBROUTINE deviatorBepr(Fr,Beprn,Dt,Gamma,B,devH,devBepr)
      
    IMPLICIT NONE
      
    !   Declarations
    REAL*8 Fr(3,3),Beprn(3,3),Dt,Gamma,B,devH(3,3),devBepr(3,3)
    REAL*8 Beprs(3,3),hydBeprs(3,3),devBeprs(3,3),Frpr(3,3)
    REAL*8 FrprT(3,3),aux(3,3)
    INTEGER i1,i2
    
    CALL MUnimodular(Fr,Frpr)
    CALL TransTens(Frpr,FrprT)
    CALL JuxtaTensTens(Frpr,Beprn,aux)
    CALL JuxtaTensTens(aux,FrprT,Beprs)
    CALL DevTens(Beprs,devBeprs,hydBeprs)
    
    DO i1=1,3
      DO i2=1,3
        devBepr(i1,i2)=(1.0d0/(1.0d0+Dt*Gamma)) & 
          & *(devBeprs(i1,i2)+B*Dt*Gamma*devH(i1,i2))
      END DO
    END DO
    
  END SUBROUTINE deviatorBepr
  
!   ====================================================================
      
!   Updated fibers directions and the homeostatic state
      
  SUBROUTINE UpdateFibersHomeo(s1n,s2n,s3n,Fr,H1,H2,H12,s1,s2,s3, &
    & S11,S22,S33,S12,S21,S13,S31,S23,S32,devH,HH)
    
    IMPLICIT NONE
      
!   Declarations
    REAL*8 s1n(3),s2n(3),s3n(3),Fr(3,3),H1,H2,H12,s1(3),s2(3),s3(3)
    REAL*8 S11(3,3),S22(3,3),S33(3,3),S12(3,3),S21(3,3),S13(3,3)
    REAL*8 S31(3,3),S23(3,3),S32(3,3),devH(3,3),HH(3,3),Frt(3,3)
    REAL*8 invFrt(3,3),EyeTens(3,3),ms1,ms3
    INTEGER i,j
    
    CALL IdentityTens(EyeTens)
    CALL TransTens(Fr,Frt)
    CALL InvTens(Frt,invFrt)
    
    CALL JuxtaTensVec(Fr,s1n,s1)
    ms1=dsqrt(s1(1)**2.0d0+s1(2)**2.0d0+s1(3)**2.0d0)
    s1(1)=s1(1)/ms1
    s1(2)=s1(2)/ms1
    s1(3)=s1(3)/ms1
    
    CALL JuxtaTensVec(invFrt,s3n,s3)
    ms3=dsqrt(s3(1)**2.0d0+s3(2)**2.0d0+s3(3)**2.0d0)
    s3(1)=s3(1)/ms3
    s3(2)=s3(2)/ms3
    s3(3)=s3(3)/ms3
    
    CALL CrossVecVec(s3,s1,s2)
    
    DO i=1,3
      DO j=1,3
        S11(i,j)=s1(i)*s1(j)
        S22(i,j)=s2(i)*s2(j)
        S33(i,j)=s3(i)*s3(j)
        S12(i,j)=s1(i)*s2(j)
        S21(i,j)=s2(i)*s1(j)
        S13(i,j)=s1(i)*s3(j)
        S31(i,j)=s3(i)*s1(j)
        S23(i,j)=s2(i)*s3(j)
        S32(i,j)=s3(i)*s2(j)        
        devH(i,j)=H1*S11(i,j)+H2*S22(i,j)-(H1+H2)*S33(i,j) &
          & +H12*(S12(i,j)+S21(i,j))          
        HH(i,j)=EyeTens(i,j)+devH(i,j)
      END DO
    END DO
    
  END SUBROUTINE UpdateFibersHomeo
  
!   ====================================================================
  
!   Compute the nonlinear function to compute {alpha1, B}
      
  SUBROUTINE nonlinearfuncs(Fr,Beprn,Dt,Gamma,devH,HH,alpha1,B, &
    & devBepr,f1,f2)
    
    IMPLICIT NONE
      
!   Declarations
    REAL*8 Fr(3,3),Beprn(3,3),Dt,Gamma,devH(3,3),HH(3,3),alpha1,B
    REAL*8 devBepr(3,3),Bepr(3,3),invBepr(3,3),EyeTens(3,3),f1,f2
    REAL*8 aux1,aux2,detdevBepr
    INTEGER i1,i2
    
    CALL IdentityTens(EyeTens)
    
    DO i1=1,3
      DO i2=1,3
        Bepr(i1,i2)=(alpha1/3.0d0)*EyeTens(i1,i2)+devBepr(i1,i2)
      END DO
    END DO
    CALL InvTens(Bepr,invBepr)
    CALL DotTensTens(invBepr,HH,aux1)
    f1=B-3.0d0/aux1
    
    CALL DotTensTens(devBepr,devBepr,aux2)
    CALL DetTens(devBepr,detdevBepr)
    f2=(alpha1/3.0d0)**3.0d0-(aux2/2.0d0)*(alpha1/3.0d0) &
      & -(1.0d0-detdevBepr)
    
  END SUBROUTINE nonlinearfuncs
      
!   ====================================================================
      
!   Compute the scalars {alpha1, B}
      
  SUBROUTINE alpha1B(hstep,tol,maxitr,Fr,Beprn,Dt,Gamma,devH,HH, &
    & alpha1n,Bn,alpha1,B)   
    
    IMPLICIT NONE
      
!   Declarations
    REAL*8 hstep,tol,maxitr,Fr(3,3),Beprn(3,3),Dt,Gamma,devH(3,3)
    REAL*8 HH(3,3),alpha1n,Bn,alpha1,B,devBepr(3,3),n,mDvars
    REAL*8 varsn(2),xpertb(2),ff(2),ffp(2),Jacob(2,2),varsnew(2)
    REAL*8 invJ(2,2),invJff(2),Dvars(2),hstepx
    INTEGER i1,i2,i3
      
    n=0.0d0
    mDvars=1.0d0
    alpha1=alpha1n
    B=Bn
    
    DO WHILE ((mDvars.GT.tol).AND.(n.LT.maxitr))
      varsn(1)=alpha1
      varsn(2)=B
      CALL deviatorBepr(Fr,Beprn,Dt,Gamma,varsn(2),devH,devBepr)
      CALL nonlinearfuncs(Fr,Beprn,Dt,Gamma,devH,HH,varsn(1), &
        & varsn(2),devBepr,ff(1),ff(2)) 
      xpertb(1)=varsn(1)
      xpertb(2)=varsn(2)
      DO i2=1,2
        hstepx=MAX(hstep,dabs(hstep*xpertb(i2)))
        xpertb(i2)=xpertb(i2)+hstepx
        CALL deviatorBepr(Fr,Beprn,Dt,Gamma,xpertb(2),devH,devBepr)
        CALL nonlinearfuncs(Fr,Beprn,Dt,Gamma,devH,HH,xpertb(1), &
          & xpertb(2),devBepr,ffp(1),ffp(2))         
        DO i3=1,2
          Jacob(i3,i2)=(ffp(i3)-ff(i3))/hstepx
        END DO
        xpertb(i2)=varsn(i2)
      END DO
      CALL InvTens2by2(Jacob,invJ)
      DO i1=1,2
        invJff(i1)=0.0d0
        DO i2=1,2
          invJff(i1)=invJff(i1)+Jacob(i1,i2)*ff(i2)
        END DO
        varsnew(i1)=varsn(i1)-invJff(i1)
        Dvars(i1)=varsnew(i1)-varsn(i1)
      END DO
      mDvars=dsqrt(Dvars(1)**2.0d0+Dvars(2)**2.0d0)
      n=n+1.0d0 
      alpha1=varsnew(1)
      B=varsnew(2)
      
      WRITE(*,*) mDvars
      WRITE(*,*) ''
      WRITE(*,*) n
      WRITE(*,*) '=========================='
      
    END DO
    
  END SUBROUTINE alpha1B
  
!   ====================================================================
      
!   compute the elastic stretches and strains

  SUBROUTINE ElasticMeasures(Je,Bepr,S11,S22,S33,S12,lame1,lame2, &
    & lame3,lamea,Ee11,Ee22,Ee33,Ee12,Eea)
      
    IMPLICIT NONE
      
!   Declarations
    REAL*8 Je,Bepr(3,3),Be(3,3),S11(3,3),S22(3,3),S33(3,3),S12(3,3)
    REAL*8 lame1,lame2,lame3,lamea,Ee11,Ee22,Ee33,Ee12,Eea,invBe(3,3)
    REAL*8 aux3
    INTEGER i1,i2
    
    DO i1=1,3
      DO i2=1,3
        Be(i1,i2)=(Je**(2.0d0/3.0d0))*Bepr(i1,i2)
      END DO
    END DO
    
    CALL InvTens(Be,invBe)
    
    CALL DotTensTens(invBe,S11,lame1)
    lame1=lame1**(-0.50d0)
    Ee11=0.50d0*((lame1**2.0d0)-1.0d0)
    
    CALL DotTensTens(invBe,S22,lame2)
    lame2=lame2**(-0.50d0)
    Ee22=0.50d0*((lame2**2.0d0)-1.0d0)
    
    CALL DotTensTens(invBe,S33,lame3)
    lame3=lame3**(-0.50d0)
    Ee33=0.50d0*((lame3**2.0d0)-1.0d0)
    
    CALL DotTensTens(invBe,S12,Ee12)
    Ee12=-0.50d0*lame1*lame2*Ee12
    
    CALL DotTensTens(Be,S33,lamea)
    lamea=Je*(lamea**(-0.50d0))
    Eea=lamea-1.0d0
    
  END SUBROUTINE ElasticMeasures
  
!   ====================================================================
      
!   Constitutive Equations (Cauchy Stresses)
      
  SUBROUTINE Constitutive(mu0,kc,kt,k1,k2,k3,k12,ka,q,Je,devBepr, &
    & Bepr,alpha1,S11,S22,S33,S12,S21,lame1,lame2,lame3,lamea,Ee12, &
    & Ee11,Ee22,Ee33,Eea,QQ,T0,T1,T2,T3,T4,T5,TT)
    
    IMPLICIT NONE
      
!   Declarations
    REAL*8 mu0,kc,kt,k1,k2,k3,k12,ka,q,mu,Je,Bepr(3,3),Be(3,3),alpha1
    REAL*8 S11(3,3),S22(3,3),S33(3,3),S12(3,3),S21(3,3),lame1,lame2
    REAL*8 lame3,lamea,Ee12,Ee11,Ee22,Ee33,Eea,QQ,T0(3,3),T1(3,3)
    REAL*8 T2(3,3),T3(3,3),T4(3,3),T5(3,3),TT(3,3),ktJe,kcJe
    REAL*8 macEe11,macEe22,macEe33,EyeTens(3,3),devBepr(3,3)
    REAL*8 invBe(3,3),invBeS33(3,3),S33invBe(3,3)
    INTEGER i1,i2
    
    CALL IdentityTens(EyeTens)
    CALL MacBrackets(Je-1,ktJe)
    CALL MacBrackets(1-Je,kcJe)
    CALL MacBrackets(Ee11,macEe11)
    CALL MacBrackets(Ee22,macEe22)
    CALL MacBrackets(Ee33,macEe33)
    DO i1=1,3
      DO i2=1,3
        Be(i1,i2)=(Je**(2.0d0/3.0d0))*Bepr(i1,i2)
      END DO
    END DO
    CALL InvTens(Be,invBe)
    CALL JuxtaTensTens(invBe,S33,invBeS33)
    CALL JuxtaTensTens(S33,invBe,S33invBe)
    
    QQ=0.50d0*kt*(ktJe**2.0d0)+0.50d0*kc*(kcJe**2.0d0) &
      & +0.50d0*(alpha1-3.0d0)+0.50d0*k1*(macEe11**2.0d0) &
      & +0.50d0*k2*(macEe22**2.0d0)+0.50d0*k3*(macEe33**2.0d0) &
      & +2.0d0*k12*(Ee12**2.0d0)+0.50d0*ka*(Eea**2.0d0) 
    
    mu=mu0*dexp(q*QQ)
    DO i1=1,3
      DO i2=1,3
        T0(i1,i2)=mu*(kt*ktJe-kc*kcJe)*EyeTens(i1,i2) &
          & +(mu/Je)*devBepr(i1,i2)
        T1(i1,i2)=((k1*mu*(lame1**2.0d0)*macEe11)/Je)*S11(i1,i2)
        T2(i1,i2)=((k2*mu*(lame2**2.0d0)*macEe22)/Je) & 
          & *(S22(i1,i2)-(2.0d0*lame2/lame1) &
          & *Ee12*(S12(i1,i2)+S21(i1,i2)))
        T3(i1,i2)=((k3*mu*(lame3**2.0d0)*macEe33)/Je) &
          & *((lame3**2.0d0)*(invBeS33(i1,i2)+S33invBe(i1,i2)) &
          & -S33(i1,i2))
        T4(i1,i2)=2.0d0*mu*k12*Ee12*(lame2/lame1/Je) &
          & *(1.0d0-4.0d0*(Ee12**2.0d0)) &
          & *(S12(i1,i2)+S21(i1,i2))
        T5(i1,i2)=(mu*ka*Eea*lamea/Je)*(S11(i1,i2)+S22(i1,i2))
        TT(i1,i2)=T0(i1,i2)+T1(i1,i2)+T2(i1,i2)+T3(i1,i2)+T4(i1,i2) &
          & +T5(i1,i2)
      END DO
    END DO
    
  END SUBROUTINE Constitutive
      
!   ====================================================================
      
  SUBROUTINE umat(ntens,nstatv,nprops,kstep,kinc,dfgrd0,dfgrd1,dtime,stress,statev,ddsdde,props)
    
!  SUBROUTINE umat(stress,statev,ddsdde,sse,spd,scd,rpl,ddsddt, &
!    & drplde,drpldt,stran,dstran,time,dtime,temp,dtemp,predef,dpred, &
!    & materl,ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt, &
!    & celent,dfgrd0,dfgrd1,noel,npt,kslay,kspt,kstep,kinc)
    
    !include 'aba_param.inc'
    
    CHARACTER*8 materl
    INTEGER ntens,nstatv,nprops
    INTEGER ndi,nshr,noel,npt,kslay,kspt,kstep,kinc
    REAL*8 sse,spd,scd,rpl,drpldt,dtime,temp,dtemp,pnewdt,celent
    REAL*8 stress(ntens),statev(nstatv),ddsdde(ntens,ntens), &
      & ddsddt(ntens),drplde(ntens),stran(ntens),dstran(ntens), &
      & dfgrd0(3,3),dfgrd1(3,3),time(2),predef(1),dpred(1), &
      & props(nprops),coords(3),drot(3,3)
    
    REAL*8 EyeTens(3,3),Dt,K,mu,a0,b0,a1,b1,m,kappas,Jn,kappan
    REAL*8 Beprn(3,3),Fr(3,3),J,Beprs(3,3),devBeprs(3,3),games,game
    REAL*8 Br(3,3),Deps,DtGam0,DtGam1,c0,c1,c2,gam,GAMMA
    REAL*8 devBepr(3,3),alpha1,Bepr(3,3),kappa,T(3,3),dJdFrFrT(3,3)
    REAL*8 ddevBeprsdFrFrT(3,3,3,3),dgamesdFrFrT(3,3)
    REAL*8 dDepsdFrFrT(3,3),dDtGam0dFrFrT(3,3),dDtGam1dFrFrT(3,3)
    REAL*8 dgamdFrFrT(3,3),dc0dFrFrT(3,3),dc1dFrFrT(3,3)
    REAL*8 dc2dFrFrT(3,3),dDtGAMdFrFrT(3,3),ddevBeprdFrFrT(3,3,3,3)
    REAL*8 tanmat(3,3,3,3)
    REAL*8 Ft1(3,3),Ft2(3,3),invFt1(3,3)
   
    INTEGER ii1,ii2,ii3,ii4

    INTEGER ind11(6),ind21(6)
    DATA    ind11/1,2,3,1,1,2/
    DATA    ind21/1,2,3,2,3,3/
    INTEGER i1,i2
      
    CALL IdentityTens(EyeTens)

!   -------------------     
!   Material constants
!   -------------------
    K=props(1)
    mu=props(2)
    a0=props(3)
    a1=props(4)
    b0=props(5)
    b1=props(6)
    m=props(7)
    kappas=props(8)
    
!   ----------------------------------------------
!   Set the initial values of the state variables
!   ----------------------------------------------
    IF ((kstep.EQ.1).AND.(kinc.EQ.1)) THEN
      statev(1)=1.0d0
      statev(2)=1.0d0
      statev(3)=1.0d0
      statev(4)=1.0d0
      statev(5)=0.0d0
      statev(6)=0.0d0
      statev(7)=0.0d0
      statev(8)=0.010d0        
      statev(9)=0.0d0
      statev(10)=0.0d0
      statev(11)=0.0d0
      statev(12)=statev(2)+statev(3)+statev(4)
      statev(13)=0.0d0
    END IF
    
!   ------------------------------------------
!   Get the old values of the state variables
!   ------------------------------------------      
    Jn=statev(1)
    Beprn(1,1)=statev(2)
    Beprn(2,2)=statev(3)
    Beprn(3,3)=statev(4)
    Beprn(1,2)=statev(5)
    Beprn(2,3)=statev(6)
    Beprn(1,3)=statev(7)
    kappan=statev(8)
    games=statev(9)
    Deps=statev(10)
    GAMMA=statev(11)
    alpha1=statev(12)
    game=statev(13)
    
    Beprn(2,1)=Beprn(1,2)
    Beprn(3,2)=Beprn(2,3)
    Beprn(3,1)=Beprn(1,3)
    
!   -----------------------------------------
!   Time step; Relative deformation gradient
!   -----------------------------------------
    Dt=dtime
    
    DO i1=1,3
      DO i2=1,3
!           <Ft1=F(t1)>
        Ft1(i1,i2)=dfgrd0(i1,i2)
!           <Ft2=F(t2)>
        Ft2(i1,i2)=dfgrd1(i1,i2)
      END DO
    END DO
      
    CALL InvTens(Ft1,invFt1)
    CALL JuxtaTensTens(Ft2,invFt1,Fr)
    
!   -------------------------------------------------------------------------------------------
!   Calculate the updated values of the state variables, stresses, and Spatial-Tangent Modulus
!   -------------------------------------------------------------------------------------------     
    CALL SmoothMultiPhase(Dt,K,mu,a0,b0,a1,b1,m,kappas,Jn,Beprn, &
      & kappan,Fr,J,Beprs,devBeprs,games,Br,Deps,DtGam0,DtGam1,c0,c1,c2, &
      & gam,GAMMA,devBepr,alpha1,Bepr,kappa,T,dJdFrFrT,ddevBeprsdFrFrT, &
      & dgamesdFrFrT,dDepsdFrFrT,dDtGam0dFrFrT,dDtGam1dFrFrT,dc0dFrFrT, &
      & dc1dFrFrT,dc2dFrFrT,dgamdFrFrT,dDtGAMdFrFrT,ddevBeprdFrFrT, &
      & tanmat)     
      
!   --------------------------- 
!   Update the state variables
!   ---------------------------
    statev(1)=J
    statev(2)=Bepr(1,1)
    statev(3)=Bepr(2,2)
    statev(4)=Bepr(3,3)
    statev(5)=Bepr(1,2)
    statev(6)=Bepr(2,3)
    statev(7)=Bepr(1,3)
    statev(8)=kappa
!      
    statev(9)=games
    statev(10)=Deps
    statev(11)=GAMMA
    statev(12)=alpha1      
    CALL EquivStrain(devBepr,game)
    statev(13)=game   
    
!   ------------------------------------------------
!   Update the Stresses and Spatial-Tangent Modulus
!   ------------------------------------------------
    IF (ntens.EQ.6) THEN
!       !   
!       ! 3D Problem 
!       !
      DO i1=1,6
        ii1=ind11(i1)
        ii2=ind21(i1)
        stress(i1)=T(ii1,ii2)
        DO i2=1,6
          ii3=ind11(i2)
          ii4=ind21(i2)
          ddsdde(i1,i2)=tanmat(ii1,ii2,ii3,ii4)
        END DO ! i2
      END DO ! i1
      
    ELSE IF (ntens.EQ.4) THEN
!       !
!       ! Axisymmetric Problem
!       !
      stress(1)=T(1,1) 
      stress(2)=T(2,2) 
      stress(3)=T(1,2) 
      stress(4)=T(3,3)        
      
      ddsdde(1,1)=tanmat(1,1,1,1)
      ddsdde(1,2)=tanmat(1,1,2,2)
      ddsdde(1,3)=tanmat(1,1,1,2)
      ddsdde(1,4)=tanmat(1,1,3,3)
      
      ddsdde(2,1)=tanmat(2,2,1,1)
      ddsdde(2,2)=tanmat(2,2,2,2)
      ddsdde(2,3)=tanmat(2,2,1,2)
      ddsdde(2,4)=tanmat(2,2,3,3)
        
      ddsdde(3,1)=tanmat(1,2,1,1)
      ddsdde(3,2)=tanmat(1,2,2,2)
      ddsdde(3,3)=tanmat(1,2,1,2)
      ddsdde(3,4)=tanmat(1,2,3,3)
        
      ddsdde(4,1)=tanmat(3,3,1,1)
      ddsdde(4,2)=tanmat(3,3,2,2)
      ddsdde(4,3)=tanmat(3,3,1,2)
      ddsdde(4,4)=tanmat(3,3,3,3)
        
    ELSE IF (ntens.EQ.3) THEN
!       !   
!       ! 2D Problem
!       !
      stress(1)=T(1,1) 
      stress(2)=T(2,2) 
      stress(3)=T(1,2)        
      
      ddsdde(1,1)=tanmat(1,1,1,1)
      ddsdde(1,2)=tanmat(1,1,2,2)
      ddsdde(1,3)=tanmat(1,1,1,2)
      
      ddsdde(2,1)=tanmat(2,2,1,1)
      ddsdde(2,2)=tanmat(2,2,2,2)
      ddsdde(2,3)=tanmat(2,2,1,2)
      
      ddsdde(3,1)=tanmat(1,2,1,1)
      ddsdde(3,2)=tanmat(1,2,2,2)
      ddsdde(3,3)=tanmat(1,2,1,2)
      
    END IF
      
    RETURN
  END SUBROUTINE umat

!   ***************************************************************************
!   ***************************************************************************  
            
!   ================
!   Main Subroutine
!   ================
           
  SUBROUTINE SmoothMultiPhase(Dt,K,mu,a0,b0,a1,b1,m,kappas,Jn, &
    & Beprn,kappan,Fr,J,Beprs,devBeprs,games,Br,Deps,DtGam0, &
    & DtGam1,c0,c1,c2,gam,GAMMA,devBepr,alpha1,Bepr,kappa,T,dJdFrFrT, &
    & ddevBeprsdFrFrT,dgamesdFrFrT,dDepsdFrFrT,dDtGam0dFrFrT, &
    & dDtGam1dFrFrT,dc0dFrFrT,dc1dFrFrT,dc2dFrFrT,dgamdFrFrT, &
    & dDtGAMdFrFrT,ddevBeprdFrFrT,tanmat)
      
    IMPLICIT NONE
      
    REAL*8 EyeTens(3,3),Dt,K,mu,a0,b0,a1,b1,m,kappas,Jn,kappan
    REAL*8 Beprn(3,3),Fr(3,3),FrT(3,3),Frpr(3,3),Jr,J,Beprs(3,3)
    REAL*8 devBeprs(3,3),hydBeprs(3,3),devgeprs(3,3),games
    REAL*8 Dbar(3,3),Br(3,3),Deps,DtGam0,DtGam1,c0,c1,c2,gam
    REAL*8 GAMMA,devBepr(3,3),devgepr(3,3),alpha1,Bepr(3,3)
    REAL*8 kappa,p,devT(3,3),T(3,3),dJdFrFrT(3,3),dgamesdFrFrT(3,3)
    REAL*8 ddevBeprsdFrFrT(3,3,3,3),dDepsdFrFrT(3,3)
    REAL*8 dDtGam0dFrFrT(3,3),dDtGam1dFrFrT(3,3),dgamdFrFrT(3,3)
    REAL*8 dc0dFrFrT(3,3),dc1dFrFrT(3,3),dc2dFrFrT(3,3)
    REAL*8 dDtGAMdFrFrT(3,3),ddevBeprdFrFrT(3,3,3,3),tanmat(3,3,3,3)
    INTEGER i1,i2
      
    CALL IdentityTens(EyeTens)
    CALL DetTens(Fr,Jr)
    J=Jr*Jn
    
    CALL MUnimodular(Fr,Frpr)
    CALL ElasticTrialDist(Frpr,Beprn,Beprs)
    CALL DevTens(Beprs,devBeprs,hydBeprs) 
    DO i1=1,3
      DO i2=1,3
        devgeprs(i1,i2)=0.50d0*devBeprs(i1,i2)
      END DO
    END DO
    
    CALL EquivStrain(devBeprs,games)
    
    CALL TransTens(Fr,FrT)
    CALL JuxtaTensTens(Fr,FrT,Br)
    DO i1=1,3
      DO i2=1,3
        Dbar(i1,i2)=0.50d0/Dt*(Br(i1,i2)-EyeTens(i1,i2))
      END DO
    END DO
    CALL EffDistStrain(Dbar,Dt,Deps) 
    
    CALL DeltatGamma(a0,a1,b0,b1,m,kappas,Dt,kappan,games,Deps, &
      & DtGam0,DtGam1,c0,c1,c2,gam,GAMMA)
    
    CALL Beprpr(devgeprs,Dt,GAMMA,devBepr)
    DO i1=1,3
      DO i2=1,3
        devgepr(i1,i2)=0.50d0*devBepr(i1,i2)
      END DO
    END DO
    CALL OneThirdAlpha(devBepr,alpha1)
    CALL Beprime(devBepr,alpha1,Bepr)
    
    CALL CauchyStress(K,mu,J,devgepr,p,devT,T)      
    CALL hardening(kappan,kappas,m,gam,kappa)
    
    CALL Tangent(K,mu,a0,a1,b0,b1,m,Dt,kappas,kappan,DtGam0, &
      & DtGam1,c0,c1,c2,gam,GAMMA,Fr,J,Beprs,devBeprs,games,devBepr, &
      & Deps,T,dJdFrFrT,ddevBeprsdFrFrT,dgamesdFrFrT,dDepsdFrFrT, &
      & dDtGam0dFrFrT,dDtGam1dFrFrT,dc0dFrFrT,dc1dFrFrT,dc2dFrFrT, &
      & dgamdFrFrT,dDtGAMdFrFrT,ddevBeprdFrFrT,tanmat)
    
    RETURN
  END SUBROUTINE SmoothMultiPhase
            
!   ===========================
!   DISTORTIONAL ELASTIC TRIAL
!   ===========================
      
!   Compute the elastic trial value Be'* of Be'    
     
  SUBROUTINE ElasticTrialDist(Frpr,Beprn,Bepret)
    
    IMPLICIT NONE
    
!   Declarations
    REAL*8 Frpr(3,3),FrprT(3,3),Beprn(3,3),aux(3,3),Bepret(3,3)     
    
    CALL TransTens(Frpr,FrprT)
    CALL JuxtaTensTens(Frpr,Beprn,aux)
    CALL JuxtaTensTens(aux,FrprT,Bepret)
    
    RETURN
  END SUBROUTINE ElasticTrialDist
  
!   ==========================
!   EQUIVALENT STRAIN gamma_e
!   ==========================
      
  SUBROUTINE EquivStrain(devBepr,game)
    
    IMPLICIT NONE
    
!   Declarations
    REAL*8 devBepr(3,3),game
      
    CALL DotTensTens(devBepr,devBepr,game)
    game=dsqrt((3.0d0/8.0d0)*game)
    
    RETURN
  END SUBROUTINE EquivStrain
  
!   ================================================
!   EFFECTIVE TOTAL DISTORTIONAL DEFORMATION D(eps)
!   ================================================
      
  SUBROUTINE EffDistStrain(Dbar,Dt,Deps)
      
    IMPLICIT NONE
    
!   Declarations
    REAL*8 Dbar(3,3),Dt,Deps,devDbar(3,3),hydDbar(3,3)
      
    CALL DevTens(Dbar,devDbar,hydDbar)
    CALL DotTensTens(devDbar,devDbar,Deps)
    Deps=Dt*dsqrt(2.0d0/3.0d0*Deps)
    
    RETURN
  END SUBROUTINE EffDistStrain
            
!   ============================
!   Deviatoric part Be'' of Be'
!   ============================
      
  SUBROUTINE Beprpr(devgeprs,Dt,GAMMA,devBepr)
    
    IMPLICIT NONE
    
!   Declarations
    REAL*8 Dt,GAMMA,devgeprs(3,3),devBepr(3,3),fac1
    INTEGER i1,i2
    
    fac1=2.0d0/(1.0d0+(Dt*GAMMA))
    
    DO i1=1,3
      DO i2=1,3
        devBepr(i1,i2)=fac1*devgeprs(i1,i2)
      END DO
    END DO
      
    RETURN
  END SUBROUTINE Beprpr
            
!   ======================================
!   FIRST NON-TRIVIAL INVARIANT (alpha_1)
!   ======================================
      
!   Calculate alpha=trace(Bepr) by solving a quadratic equation
      
  SUBROUTINE OneThirdAlpha(devBepr,alpha1)
    
    IMPLICIT NONE
    
!   Declarations
    REAL*8 devBepr(3,3),alpha1
    REAL*8 fac1,fac2,detdevBepr,dotprod
    
    CALL DetTens(devBepr,detdevBepr)
    CALL DotTensTens(devBepr,devBepr,dotprod)
      
    fac1=2.0d0*dotprod/3.0d0
    IF (fac1.EQ.0.0d0) THEN
      alpha1=3.0d0
    ELSE
      fac2=(4.0d0*(1.0d0-detdevBepr))/(fac1**1.5d0)
      IF (fac2.GE.1.0d0) THEN
        alpha1=3.0d0*dsqrt(fac1)*dcosh(dacosh(fac2)/3.0d0)
      ELSE
        alpha1=3.0d0*dsqrt(fac1)*dcos(dacos(fac2)/3.0d0)
      END IF
    END IF
    
    RETURN
  END SUBROUTINE OneThirdAlpha
            
!   =========================
!   ELASTIC DISTORTION (Be')
!   =========================
      
  SUBROUTINE Beprime(devBepr,alpha1,Bepr)
      
    IMPLICIT NONE

!   Declarations
    REAL*8 EyeTens(3,3),devBepr(3,3),alpha1,Bepr(3,3)
    INTEGER i1,i2
    
    CALL IdentityTens(EyeTens)
    DO i1=1,3
      DO i2=1,3
        Bepr(i1,i2)=alpha1/3.0d0*EyeTens(i1,i2)+devBepr(i1,i2)
      END DO
    END DO
    
    RETURN
  END SUBROUTINE Beprime
  
!   ======================================
!   AUXILIARY VARIABLES (gamma and GAMMA)
!   ======================================
      
  SUBROUTINE DeltatGamma(a0,a1,b0,b1,m,kappas,Dt,kappan,games, &
    & Deps,DtGam0,DtGam1,c0,c1,c2,gam,GAMMA)
      
    IMPLICIT NONE
    
!   Declarations
    REAL*8 a0,a1,b0,b1,Dt,games,kappan,gam,GAMMA,DtGam0,DtGam1
    REAL*8 c0,c1,c2,Deps,m,kappas
    
    DtGam0=Dt*a0+b0*Deps
    DtGam1=Dt*a1+b1*Deps  
    c0=DtGam1*(games-(1.0d0+DtGam0)*kappan)
    
    IF (c0.LE.0.0d0) THEN
      gam=0.0d0
    ELSE
      c1=games+DtGam1*(kappan-m*(games-(1+DtGam0)*kappas))
      IF (m.EQ.0.0d0) THEN
        gam=c0/c1
      ELSE
        c2=m*(games+DtGam1*kappas)
        gam=(-c1+dsqrt(c1**2.0d0+4.0d0*c0*c2))/2.0d0/c2
      END IF
    END IF
      
    GAMMA=DtGam0/Dt+gam/Dt
      
    RETURN
  END SUBROUTINE DeltatGamma
            
!   =====================
!   CAUCHY STRESS TENSOR
!   =====================
      
  SUBROUTINE CauchyStress(K,mu,J,devgepr,p,devT,T)
    
    IMPLICIT NONE
      
!   Declarations
    REAL*8 K,mu,J,devgepr(3,3),p,devT(3,3),T(3,3)
    REAL*8 EyeTens(3,3)
    INTEGER i1,i2
      
    CALL IdentityTens(EyeTens)
      
    p=-K*(J-1.0d0)
    DO i1=1,3
      DO i2=1,3
        devT(i1,i2)=2.0d0*mu*devgepr(i1,i2)/J
        T(i1,i2)=-p*EyeTens(i1,i2)+devT(i1,i2)
      END DO
    END DO
      
    RETURN
  END SUBROUTINE CauchyStress
            
!   ===================
!   HARDENING VARIABLE
!   ===================
      
  SUBROUTINE hardening(kappan,kappas,m,gam,kappa)
      
    IMPLICIT NONE
      
!   Declarations
    REAL*8 kappan,kappas,m,gam,kappa
    
    kappa=(kappan+m*kappas*gam)/(1.0d0+m*gam)
    
    RETURN
  END SUBROUTINE hardening
  
!   =======================
!   SPATIAL TANGENT MODULI
!   =======================
      
  SUBROUTINE Tangent(K,mu,a0,a1,b0,b1,m,Dt,kappas,kappan,DtGam0, &
    & DtGam1,c0,c1,c2,gam,GAMMA,Fr,J,Beprs,devBeprs,games,devBepr, &
    & Deps,T,dJdFrFrT,ddevBeprsdFrFrT,dgamesdFrFrT,dDepsdFrFrT, &
    & dDtGam0dFrFrT,dDtGam1dFrFrT,dc0dFrFrT,dc1dFrFrT,dc2dFrFrT, &
    & dgamdFrFrT,dDtGAMdFrFrT,ddevBeprdFrFrT,tanmat)
      
    IMPLICIT NONE
      
!   Declarations
    REAL*8 EyeTens(3,3),K,mu,a0,a1,b0,b1,m,Dt,kappas,kappan
    REAL*8 DtGam0,DtGam1,c0,c1,c2,gam,GAMMA,Fr(3,3),J,Beprs(3,3)
    REAL*8 devBeprs(3,3),games,devBepr(3,3),Deps,ex(3,3),hx(3,3)
    REAL*8 Frt(3,3),Br(3,3),T(3,3),devBr(3,3),dJdFrFrT(3,3)
    REAL*8 ddevBeprsdFrFrT(3,3,3,3),dgamesdFrFrT(3,3)
    REAL*8 dDepsdFrFrT(3,3),dDtGam0dFrFrT(3,3),dDtGam1dFrFrT(3,3)
    REAL*8 cf0,dc0dFrFrT(3,3),cf1,dc1dFrFrT(3,3),dc2dFrFrT(3,3)
    REAL*8 dgamdFrFrT(3,3),dDtGAMdFrFrT(3,3),ddevBeprdFrFrT(3,3,3,3)
    REAL*8 cft1,cft2,tanmat(3,3,3,3)
    INTEGER i1,i2,i3,i4
    
    
    CALL IdentityTens(EyeTens)
    CALL JuxtaTensTens(devBeprs,Beprs,ex)
    CALL DevTens(ex,dgamesdFrFrT,hx)
    CALL TransTens(Fr,Frt)
    CALL JuxtaTensTens(Fr,Frt,Br)
    CALL DevTens(Br,devBr,hx)
    CALL JuxtaTensTens(devBr,Br,dDepsdFrFrT)
    
!     calculate diff(J,Fr)*Fr^T     
    DO i1=1,3
      DO i2=1,3
        dJdFrFrT(i1,i2)=J*EyeTens(i1,i2)
      END DO
    END DO
      
!     calculate diff(gammae*,Fr)*Fr^T     
    DO i1=1,3
      DO i2=1,3
        dgamesdFrFrT(i1,i2)=3.0d0/4.0d0/games*dgamesdFrFrT(i1,i2)
      END DO
    END DO
      
!     calculate diff(Depsilon,Fr)*Fr^T 
    DO i1=1,3
      DO i2=1,3
        dDepsdFrFrT(i1,i2)=1.0d0/3.0d0/Deps*dDepsdFrFrT(i1,i2)
      END DO
    END DO
    
!     calculate diff(Dt*Gamma0,Fr)*Fr^T      
    DO i1=1,3
      DO i2=1,3
        dDtGam0dFrFrT(i1,i2)=b0*dDepsdFrFrT(i1,i2)
      END DO
    END DO
      
!     calculate diff(Dt*Gamma1,Fr)*Fr^T   
    DO i1=1,3
      DO i2=1,3
        dDtGam1dFrFrT(i1,i2)=b1*dDepsdFrFrT(i1,i2)
      END DO
    END DO
    
!     calculate diff(gamma,Fr)*Fr^T   
    cf0=games-(1.0d0+DtGam0)*kappan
    cf1=kappan-m*(games-(1.0d0+DtGam0)*kappas)
    DO i1=1,3
      DO i2=1,3
        IF (c0.LE.0.0d0) THEN
          dgamdFrFrT(i1,i2)=0.0d0
        ELSE
          dc0dFrFrT(i1,i2)=cf0*dDtGam1dFrFrT(i1,i2) & 
            & +DtGam1*(dgamesdFrFrT(i1,i2)-kappan*dDtGam0dFrFrT(i1,i2))
          dc1dFrFrT(i1,i2)=(1.0d0-m*DtGam1)*dgamesdFrFrT(i1,i2) &
            & +m*kappas*DtGam1*dDtGam0dFrFrT(i1,i2) &
            & +cf1*dDtGam1dFrFrT(i1,i2)
          dc2dFrFrT(i1,i2)=m*(dgamesdFrFrT(i1,i2) &
            & +kappas*dDtGam1dFrFrT(i1,i2))
          IF (m.EQ.0.0d0) THEN
            dgamdFrFrT(i1,i2)=(1.0d0/c1)*(dc0dFrFrT(i1,i2) &
              & -gam*dc1dFrFrT(i1,i2))
          ELSE
            dgamdFrFrT(i1,i2)=(gam/(2.0d0*c0-c1*gam)) &
              &  *(dc0dFrFrT(i1,i2)-gam*dc1dFrFrT(i1,i2) &
              &  -(gam**2.0d0)*dc2dFrFrT(i1,i2))
          END IF
          dDtGAMdFrFrT(i1,i2)=dDtGam0dFrFrT(i1,i2)+dgamdFrFrT(i1,i2)
        END IF
      END DO
    END DO
      
!     calculate diff(Be''*,Fr)*Fr^T         
    DO i1=1,3
      DO i2=1,3
        DO i3=1,3
          DO i4=1,3
            ddevBeprsdFrFrT(i1,i2,i3,i4)= &
              & EyeTens(i1,i3)*Beprs(i2,i4) &
              & +Beprs(i1,i4)*EyeTens(i2,i3) &
              & -(2.0d0/3.0d0)*Beprs(i1,i2)*EyeTens(i3,i4) &
              & -(2.0d0/3.0d0)*EyeTens(i1,i2)*devBeprs(i3,i4)
          END DO
        END DO
      END DO
    END DO
      
!     calculate diff(Be'',Fr)*Fr^T      
    DO i1=1,3
      DO i2=1,3
        DO i3=1,3
          DO i4=1,3
            ddevBeprdFrFrT(i1,i2,i3,i4)=(1.0d0/(1.0d0+Dt*GAMMA))* &
              & (ddevBeprsdFrFrT(i1,i2,i3,i4)-devBepr(i1,i2)* &
              & dDtGAMdFrFrT(i3,i4))
          END DO
        END DO
      END DO
    END DO

!     calculate the spatial tangent moduli     
    cft1=K/J*(2*J-1.0d0)
    cft2=mu/J
    DO i1=1,3
      DO i2=1,3
        DO i3=1,3
          DO i4=1,3
            tanmat(i1,i2,i3,i4)=cft1*EyeTens(i1,i2)*dJdFrFrT(i3,i4) &
              & +cft2*ddevBeprdFrFrT(i1,i2,i3,i4) &
              & -T(i1,i4)*EyeTens(i2,i3)
          END DO
        END DO
      END DO
    END DO
      
    RETURN
  END SUBROUTINE Tangent

!   ***************************************************************************
!   ***************************************************************************  
      
!   Juxtaposition of a 3x3 tensor and a 3x1 vector
!   Also holds for a Dot Product between a 2nd order tensor and a vector
      
  SUBROUTINE JuxtaTensVec(A,v,Av)
      
    IMPLICIT NONE
      
!   Declarations
    REAL*8 A(3,3),v(3),Av(3)
    INTEGER i,j
    
    DO i=1,3
      Av(i)=0.0d0
      DO j=1,3
        Av(i)=Av(i)+A(i,j)*v(j)
      END DO
    END DO
    
    RETURN
  END SUBROUTINE JuxtaTensVec
      
!   ====================================================================
      
!   Juxtaposition between two 2nd order tensors
      
  SUBROUTINE JuxtaTensTens(A,B,AB)
    
    IMPLICIT NONE
    
!   Declarations
    REAL*8 A(3,3),B(3,3),AB(3,3)
    INTEGER i,j,k

    DO i=1,3
      DO j=1,3        
        AB(i,j)=0.0d0
        DO k=1,3
          AB(i,j)=AB(i,j)+A(i,k)*B(k,j)
        END DO
      END DO
    END DO
    
    RETURN
  END SUBROUTINE JuxtaTensTens
  
!   ====================================================================
      
!   Transpose of a 2nd order tensor
      
  SUBROUTINE TransTens(A,At)
      
    IMPLICIT NONE
      
!   Declarations
    REAL*8 A(3,3),At(3,3)
    INTEGER i,j
    
    DO i=1,3
      DO j=1,3
        At(i,j)=A(j,i)
      END DO
    END DO
    
    RETURN
  END SUBROUTINE TransTens
      
!   ====================================================================
      
!   Determinant of a 2nd order tensor
      
  SUBROUTINE DetTens(A,detA)
      
    IMPLICIT NONE
      
!   Declarations
    REAL*8 A(3,3),A2(3,3),A3(3,3),detA,trA,trA2,trA3
    
    CALL JuxtaTensTens(A,A,A2)
    CALL JuxtaTensTens(A,A2,A3)
    
    CALL TraceTens(A,trA)
    CALL TraceTens(A2,trA2)
    CALL TraceTens(A3,trA3)
    
    detA=1.0d0/3.0d0*(trA3-trA*trA2+0.50d0*(trA**2.0d0-trA2)*trA)
    
    RETURN
  END SUBROUTINE DetTens
  
!   ====================================================================           
      
!   Trace of a 2nd order tensor
      
  SUBROUTINE TraceTens(A,trA)
    
    IMPLICIT NONE

!   Declarations
    REAL*8 A(3,3)
    REAL*8 trA
    
    trA=A(1,1)+A(2,2)+A(3,3)
    
    RETURN
  END SUBROUTINE TraceTens
  
!   ====================================================================    
      
!   Inverse of a 2nd order tensor
      
  SUBROUTINE InvTens(A,invA)
      
    IMPLICIT NONE
      
!   Declarations
    REAL*8 EyeTens(3,3),A(3,3),detA,I11,A2(3,3),trA2,I22,invA(3,3)
    INTEGER i1,i2
    
    CALL IdentityTens(EyeTens)
    CALL DetTens(A,detA)   
    I11=A(1,1)+A(2,2)+A(3,3)      
    CALL JuxtaTensTens(A,A,A2)
    trA2=A2(1,1)+A2(2,2)+A2(3,3)
    I22=0.50d0*(I11**2.0d0-trA2)
    
    DO i1=1,3
      DO i2=1,3
        invA(i1,i2)=1.0d0/detA*(A2(i1,i2)-I11*A(i1,i2) &
          & +I22*EyeTens(i1,i2))
      END DO
    END DO
    
    RETURN
  END SUBROUTINE InvTens
  
!   ====================================================================
      
!   Dot Product between two 2nd order tensors
      
  SUBROUTINE DotTensTens(A,B,AdotB)
      
    IMPLICIT NONE
    
!   Declarations
    REAL*8 A(3,3),B(3,3),Bt(3,3),ABt(3,3)
    REAL*8 trABt,AdotB
    
    CALL TransTens(B,Bt)
    CALL JuxtaTensTens(A,Bt,ABt)
    CALL TraceTens(ABt,trABt)
    
    AdotB=trABt
    
    RETURN
  END SUBROUTINE DotTensTens
      
!   ====================================================================
      
!   2nd order identity tensor
      
  SUBROUTINE IdentityTens(EyeTens)
      
    IMPLICIT NONE
      
!   Declarations
    REAL*8 EyeTens(3,3)
      
    EyeTens(1,1)=1.0d0
    EyeTens(1,2)=0.0d0
    EyeTens(1,3)=0.0d0
    EyeTens(2,1)=0.0d0
    EyeTens(2,2)=1.0d0
    EyeTens(2,3)=0.0d0
    EyeTens(3,1)=0.0d0
    EyeTens(3,2)=0.0d0
    EyeTens(3,3)=1.0d0
    
    RETURN
  END SUBROUTINE IdentityTens
      
!   ====================================================================
      
!   2nd order sparse tensor
      
  SUBROUTINE ZeroTens(zero3)
      
    IMPLICIT NONE
      
!   Declarations
    REAL*8 zero3(3,3)
    
    zero3(1,1)=0.0d0
    zero3(1,2)=0.0d0
    zero3(1,3)=0.0d0
    zero3(2,1)=0.0d0
    zero3(2,2)=0.0d0
    zero3(2,3)=0.0d0
    zero3(3,1)=0.0d0
    zero3(3,2)=0.0d0
    zero3(3,3)=0.0d0
      
    RETURN
  END SUBROUTINE ZeroTens
      
!   ====================================================================
      
!   Deviator of a 2nd order tensor
      
  SUBROUTINE DevTens(A,devA,hydA)
      
    IMPLICIT NONE
      
!   Declarations
    REAL*8 A(3,3),devA(3,3),EyeTens(3,3),hydA(3,3)
    REAL*8 p,trA
    INTEGER i,j
    
    CALL TraceTens(A,trA)
    p=-trA/3.0d0
    CALL IdentityTens(EyeTens)
    
    DO i=1,3
      DO j=1,3
        hydA(i,j)=-p*EyeTens(i,j)
        devA(i,j)=A(i,j)-hydA(i,j)
      END DO
    END DO
      
    RETURN
  END SUBROUTINE DevTens
      
!   ====================================================================
      
!   Unimodular 2nd order tensor
      
  SUBROUTINE MUnimodular(A,Apr)
      
    IMPLICIT NONE
      
!   Declarations
    REAL*8 A(3,3),Apr(3,3)
    REAL*8 detA,power,fac
    INTEGER i,j
    power=-1.0d0/3.0d0
    
    CALL DetTens(A,detA)
    
    DO i=1,3
      DO j=1,3
        Apr(i,j)=SIGN(dabs(detA)**power, detA)*A(i,j) 
      END DO
    END DO
      
!       do i=1,3
!         do j=1,3
!           Apr(i,j)=fac*A(i,j)
!         enddo
!       enddo  
      
    RETURN
  END SUBROUTINE MUnimodular
      
!   ====================================================================
      
!   Cross Product between two vectors
      
  SUBROUTINE CrossVecVec(a,b,axb)
      
    IMPLICIT NONE
      
!   Declarations
    REAL*8 a(3),b(3),axb(3)
      
    axb(1)=a(2)*b(3)-a(3)*b(2)
    axb(2)=a(3)*b(1)-a(1)*b(3)
    axb(3)=a(1)*b(2)-a(2)*b(1)
    
    RETURN
  END SUBROUTINE CrossVecVec
      
!   ====================================================================
      
!   Dot Product between two vectors
      
  SUBROUTINE DotVecVec(a,b,adotb)
      
    IMPLICIT NONE
      
    REAL*8 a(3),b(3)
    REAL*8 adotb
      
    adotb=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
    
    RETURN
  END SUBROUTINE DotVecVec
      
!   ====================================================================
      
!   Tensor Product between two vectors
      
  SUBROUTINE TensProd(a,b,aoxb)
      
    IMPLICIT NONE
      
    REAL*8 a(3),b(3),aoxb(3,3)
    INTEGER i,j
      
    DO i=1,3
      DO j=1,3
        aoxb(i,j)=a(i)*b(j)
      END DO
    END DO
    
    RETURN
  END SUBROUTINE TensProd
      
!   ====================================================================
      
!   Tensor Product berween two 2nd order tensors

  SUBROUTINE TensProd33(A,B,AoxB)
      
    IMPLICIT NONE
      
    REAL*8 A(3,3),B(3,3),AoxB(3,3,3,3)
    INTEGER i,j,m,n
    
    DO i=1,3
      DO j=1,3
        DO m=1,3
          DO n=1,3
            AoxB(i,j,m,n)=A(i,j)*B(m,n)
          END DO
        END DO
      END DO
    END DO
    
    RETURN
  END SUBROUTINE TensProd33
      
!   ====================================================================
      
!   The Operation "oplus" between two 2nd order tensors

  SUBROUTINE oplus(A,B,AopB)
      
    IMPLICIT NONE
      
    REAL*8 A(3,3),B(3,3),AopB(3,3,3,3)
    INTEGER i,j,m,n
    
    DO i=1,3
      DO j=1,3
        DO m=1,3
          DO n=1,3
            AopB(i,j,m,n)=A(i,n)*B(j,m)
          END DO
        END DO
      END DO
    END DO
      
    RETURN
  END SUBROUTINE oplus
      
!   ====================================================================      
!   The Operation "ominus" between two 2nd order tensors

  SUBROUTINE ominus(A,B,AomB)
      
    IMPLICIT NONE
      
    REAL*8 A(3,3),B(3,3),AomB(3,3,3,3)
    INTEGER i,j,m,n
    
    DO i=1,3
      DO j=1,3
        DO m=1,3
          DO n=1,3
            AomB(i,j,m,n)=A(i,m)*B(j,n)
          END DO
        END DO
      END DO
    END DO
    
    RETURN
  END SUBROUTINE ominus
      
!   ====================================================================
      
!   Dot Product between a 2nd order tensor and a 4th order tensor
      
  SUBROUTINE DotTens2Tens4(A2,B4,A2dotB4)
    
    IMPLICIT NONE
    
!   Declarations
    REAL*8 A2(3,3),B4(3,3,3,3),A2dotB4(3,3)
    INTEGER i,j,m,n
    
    DO i=1,3
      DO j=1,3
        A2dotB4(i,j)=0.0d0
        DO m=1,3
          DO n=1,3
            A2dotB4(i,j)=A2dotB4(i,j)+A2(m,n)*B4(m,n,i,j)
          END DO
        END DO
      END DO
    END DO
    
    RETURN
  END SUBROUTINE DotTens2Tens4
      
!   ====================================================================
      
!   Dot Product between a 4th order tensor and a 2nd order tensor
      
  SUBROUTINE DotTens4Tens2(A4,B2,A4dotB2)
    
    IMPLICIT NONE
    
!   Declarations
    REAL*8 A4(3,3,3,3),B2(3,3),A4dotB2(3,3)
    INTEGER i,j,m,n
    
    DO i=1,3
      DO j=1,3
        A4dotB2(i,j)=0.0d0
        DO m=1,3
          DO n=1,3
            A4dotB2(i,j)=A4dotB2(i,j)+A4(i,j,m,n)*B2(m,n)
          END DO
        END DO
      END DO
    END DO
    
    RETURN
  END SUBROUTINE DotTens4Tens2
      
!   ====================================================================    
      
!   Inverse of a 2 by 2 matrix
      
  SUBROUTINE InvTens2by2(A,invA)
      
    IMPLICIT NONE
      
!   Declarations
    REAL*8 A(2,2),detA,invA(2,2)
    
    detA=A(1,1)*A(2,2)-A(1,2)*A(2,1)
    invA(1,1)=1.0d0/detA*A(2,2)
    invA(1,2)=-1.0d0/detA*A(1,2)
    invA(2,1)=-1.0d0/detA*A(2,1)
    invA(2,2)=1.0d0/detA*A(1,1)
    
  END SUBROUTINE InvTens2by2
      
!   ====================================================================  
        
!   Macauley Brackets
      
  SUBROUTINE MacBrackets(x,Macx)
    
    IMPLICIT NONE
      
!   Declarations
    REAL*8 x,Macx
    
    Macx=MAX(x,0.0d0)
      
  END SUBROUTINE MacBrackets
  
END MODULE FINITE_ELASTICITY_ROUTINES
