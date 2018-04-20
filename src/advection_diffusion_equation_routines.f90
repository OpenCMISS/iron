!> \file
!> \author Chris Bradley
!> \brief This module handles all diffusion equation routines.
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
!> Contributor(s):
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

!>This module handles all advection-diffusion equation routines.
MODULE ADVECTION_DIFFUSION_EQUATION_ROUTINES

  USE ANALYTIC_ANALYSIS_ROUTINES
  USE BaseRoutines
  USE BasisRoutines
  USE BasisAccessRoutines
  USE BOUNDARY_CONDITIONS_ROUTINES
  USE CONSTANTS
  USE CONTROL_LOOP_ROUTINES
  USE ControlLoopAccessRoutines
  USE DistributedMatrixVector
  USE DOMAIN_MAPPINGS
  USE EquationsRoutines
  USE EquationsAccessRoutines
  USE EquationsMappingRoutines
  USE EquationsMatricesRoutines
  USE EquationsSetConstants
  USE EquationsSetAccessRoutines
  USE FIELD_ROUTINES
  USE FieldAccessRoutines
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE Kinds
  USE MatrixVector
  USE PROBLEM_CONSTANTS
  USE Strings
  USE SOLVER_ROUTINES
  USE SolverAccessRoutines
  USE Timer
  USE Types
! temporary input for setting velocity field
  USE FLUID_MECHANICS_IO_ROUTINES

#include "macros.h"

  IMPLICIT NONE

  PRIVATE ADVECTION_DIFFUSION_POST_SOLVE_OUTPUT_DATA

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC AdvectionDiffusion_EquationsSetSetup

  PUBLIC AdvectionDiffusion_EquationsSetSolnMethodSet
  
  PUBLIC AdvectionDiffusion_BoundaryConditionsAnalyticCalculate

  PUBLIC AdvectionDiffusion_EquationsSetSpecificationSet
  
  PUBLIC AdvectionDiffusion_FiniteElementCalculate

  PUBLIC AdvectionDiffusion_ProblemSpecificationSet
  
  PUBLIC ADVECTION_DIFFUSION_EQUATION_PROBLEM_SETUP

  PUBLIC AdvectionDiffusion_PreSolveUpdateInputData,AdvectionDiffusion_PreSolveGetSourceValue, &
    & AdvectionDiffusion_PreSolveStoreCurrentSoln

  PUBLIC AdvectionDiffusion_PreSolve,AdvectionDiffusion_PostSolve
  
CONTAINS

  !
  !================================================================================================================================
  !


  !>Calculates the analytic solution and sets the boundary conditions for an analytic problem.
  !>For the advection-diffusion analytic example it is required that the advective velocity
  !>and the source field are set to a particular analytic value, which is performed within this subroutine.
  SUBROUTINE AdvectionDiffusion_BoundaryConditionsAnalyticCalculate(EQUATIONS_SET,BOUNDARY_CONDITIONS,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,deriv_idx,dim_idx,local_ny,node_idx,NUMBER_OF_DIMENSIONS,variable_idx,variable_type
    REAL(DP) :: VALUE,X(3),VALUE_SOURCE,VALUE_INDEPENDENT,VALUE_MATERIAL
    REAL(DP), POINTER :: GEOMETRIC_PARAMETERS(:)
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,GEOMETRIC_FIELD,INDEPENDENT_FIELD,SOURCE_FIELD,MATERIALS_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE,GEOMETRIC_VARIABLE
    TYPE(VARYING_STRING) :: localError
    REAL(DP) :: alpha, phi, Peclet,tanphi

    ENTERS("AdvectionDiffusion_BoundaryConditionsAnalyticCalculate",err,error,*999)

    NULLIFY(GEOMETRIC_VARIABLE)
    NULLIFY(GEOMETRIC_PARAMETERS)

    alpha = 1.0_DP
    phi = 0.2_DP
    tanphi = TAN(phi)
    Peclet= 10.0_DP

    !>Set the analytic boundary conditions 
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
        DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
        IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
          GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
          IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN            
            CALL FIELD_NUMBER_OF_COMPONENTS_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS,err,error,*999)
            CALL Field_VariableGet(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,GEOMETRIC_VARIABLE,err,error,*999)
            CALL FIELD_PARAMETER_SET_DATA_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS, &
              & err,error,*999)
            IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
              DO variable_idx=1,DEPENDENT_FIELD%NUMBER_OF_VARIABLES
                variable_type=DEPENDENT_FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE
                FIELD_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%ptr
                IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                  CALL FIELD_PARAMETER_SET_CREATE(DEPENDENT_FIELD,variable_type,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
                  DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    IF(FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN
                      DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
                      IF(ASSOCIATED(DOMAIN)) THEN
                        IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                          DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                          IF(ASSOCIATED(DOMAIN_NODES)) THEN
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
                                CASE(EQUATIONS_SET_ADVECTION_DIFFUSION_EQUATION_TWO_DIM_1)
                                  !This is a steady-state solution of advection-diffusion equation
                                  !Velocity field takes form v(x,y)=(sin(6y),cos(6x))
                                  !Solution is u(x,y)=tanh(1 - alpha.(x.tan(Phi) - y))
                                 SELECT CASE(variable_type)
                                  CASE(FIELD_U_VARIABLE_TYPE)
                                    SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                    CASE(NO_GLOBAL_DERIV)
                                      VALUE=TANH(1.0-alpha*(X(1)*tanphi-X(2)))
                                    CASE(GLOBAL_DERIV_S1)
                                      CALL FlagError("Not implemented.",err,error,*999)
                                    CASE(GLOBAL_DERIV_S2)
                                      CALL FlagError("Not implemented.",err,error,*999)
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      CALL FlagError("Not implmented.",err,error,*999)
                                    CASE DEFAULT
                                      localError="The global derivative index of "//TRIM(NumberToVString( &
                                        DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                        & err,error))//" is invalid."
                                      CALL FlagError(localError,err,error,*999)
                                    END SELECT
                                  CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                                   SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                    CASE(NO_GLOBAL_DERIV)
                                      VALUE=0.0_DP
                                    CASE(GLOBAL_DERIV_S1)
                                      CALL FlagError("Not implemented.",err,error,*999)
                                    CASE(GLOBAL_DERIV_S2)
                                      CALL FlagError("Not implemented.",err,error,*999)
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      CALL FlagError("Not implemented.",err,error,*999)
                                    CASE DEFAULT
                                      localError="The global derivative index of "//TRIM(NumberToVString( &
                                        & DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                        & err,error))//" is invalid."
                                      CALL FlagError(localError,err,error,*999)
                                    END SELECT
                                  CASE DEFAULT
                                    localError="The variable type of "//TRIM(NumberToVString(variable_type,"*",err,error))// &
                                      & " is invalid."
                                    CALL FlagError(localError,err,error,*999)
                                  END SELECT
                                CASE DEFAULT
                                  localError="The analytic function type of "// &
                                    & TRIM(NumberToVString(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
                                    & " is invalid."
                                  CALL FlagError(localError,err,error,*999)
                                END SELECT
                                !Default to version 1 of each node derivative
                                local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                  & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,variable_type, &
                                  & FIELD_ANALYTIC_VALUES_SET_TYPE,local_ny,VALUE,err,error,*999)
                                IF(variable_type==FIELD_U_VARIABLE_TYPE) THEN
                                  IF(DOMAIN_NODES%NODES(node_idx)%BOUNDARY_NODE) THEN
                                    !If we are a boundary node then set the analytic value on the boundary
                                    CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,DEPENDENT_FIELD,variable_type, &
                                      & local_ny,BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
                                  ENDIF
                                ENDIF
                              ENDDO !deriv_idx
                            ENDDO !node_idx
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
                  ENDDO !component_idx
                  CALL FIELD_PARAMETER_SET_UPDATE_START(DEPENDENT_FIELD,variable_type,FIELD_ANALYTIC_VALUES_SET_TYPE, &
                    & err,error,*999)
                  CALL FIELD_PARAMETER_SET_UPDATE_FINISH(DEPENDENT_FIELD,variable_type,FIELD_ANALYTIC_VALUES_SET_TYPE, &
                    & err,error,*999)
                ELSE
                  CALL FlagError("Field variable is not associated.",err,error,*999)
                ENDIF

              ENDDO !variable_idx
              CALL FIELD_PARAMETER_SET_DATA_RESTORE(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
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
    
    !>Set the independent field (i.e. the advective velocity) to a specified analytical function 
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
        INDEPENDENT_FIELD=>EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD
        IF(ASSOCIATED(INDEPENDENT_FIELD)) THEN
          GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
          IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN            
            CALL FIELD_NUMBER_OF_COMPONENTS_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS,err,error,*999)
            NULLIFY(GEOMETRIC_VARIABLE)
            CALL Field_VariableGet(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,GEOMETRIC_VARIABLE,err,error,*999)
            CALL FIELD_PARAMETER_SET_DATA_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS, &
              & err,error,*999)
            DO variable_idx=1,INDEPENDENT_FIELD%NUMBER_OF_VARIABLES
              variable_type=INDEPENDENT_FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE
              FIELD_VARIABLE=>INDEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%ptr
              IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                  IF(FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN
                    DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
                    IF(ASSOCIATED(DOMAIN)) THEN
                      IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                        DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                        IF(ASSOCIATED(DOMAIN_NODES)) THEN
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
                              CASE(EQUATIONS_SET_ADVECTION_DIFFUSION_EQUATION_TWO_DIM_1)
                                !Velocity field takes form v(x,y)=(sin(6y),cos(6x))
                                IF(component_idx==1) THEN
                                  VALUE_INDEPENDENT=SIN(6*X(1))
                                ELSE
                                  VALUE_INDEPENDENT=COS(6*X(2))           
                                ENDIF
                              CASE DEFAULT
                                localError="The analytic function type of "// &
                                  & TRIM(NumberToVString(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
                                  & " is invalid."
                                CALL FlagError(localError,err,error,*999)
                              END SELECT
                              !Default to version 1 of each node derivative
                              local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                              CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                                & FIELD_VALUES_SET_TYPE,local_ny,VALUE_INDEPENDENT,err,error,*999)
                            ENDDO !deriv_idx
                          ENDDO !node_idx
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
                ENDDO !component_idx
                CALL FIELD_PARAMETER_SET_UPDATE_START(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & err,error,*999)
                CALL FIELD_PARAMETER_SET_UPDATE_FINISH(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & err,error,*999)
              ELSE
                CALL FlagError("Field variable is not associated.",err,error,*999)
              ENDIF

            ENDDO !variable_idx
            CALL FIELD_PARAMETER_SET_DATA_RESTORE(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & GEOMETRIC_PARAMETERS,err,error,*999)
          ELSE
            CALL FlagError("Equations set geometric field is not associated.",err,error,*999)
          ENDIF            
        ELSE
          CALL FlagError("Equations set independent field is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Equations set analytic is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

    !>Set the source field to a specified analytical function
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
        SOURCE_FIELD=>EQUATIONS_SET%SOURCE%SOURCE_FIELD
        IF(ASSOCIATED(SOURCE_FIELD)) THEN
          GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
          IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN            
            CALL FIELD_NUMBER_OF_COMPONENTS_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS,err,error,*999)
            NULLIFY(GEOMETRIC_VARIABLE)
            CALL Field_VariableGet(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,GEOMETRIC_VARIABLE,err,error,*999)
            CALL FIELD_PARAMETER_SET_DATA_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS, &
              & err,error,*999)
            DO variable_idx=1,SOURCE_FIELD%NUMBER_OF_VARIABLES
              variable_type=SOURCE_FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE
              FIELD_VARIABLE=>SOURCE_FIELD%VARIABLE_TYPE_MAP(variable_type)%ptr
              IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                  IF(FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN
                    DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
                    IF(ASSOCIATED(DOMAIN)) THEN
                      IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                        DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                        IF(ASSOCIATED(DOMAIN_NODES)) THEN
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
                              CASE(EQUATIONS_SET_ADVECTION_DIFFUSION_EQUATION_TWO_DIM_1)
                               VALUE_SOURCE= (1.0/Peclet)*(2.0*TANH(-0.1E1+alpha*(tanphi*X(1)-X(2)))*(1.0-(TANH(-0.1E1+alpha*( &
                                  & tanphi*X(1)-X(2)))**2))*alpha*alpha*tanphi*tanphi+2.0*TANH(-0.1E1+alpha*(tanphi*X(1)-X(2)) &
                                  & )*(1.0-(TANH(-0.1E1+alpha*(tanphi*X(1)-X(2)))**2))*alpha*alpha-Peclet*(-SIN(6.0*X(2) &
                                  & )*(1.0-(TANH(-0.1E1+alpha*(tanphi*X(1)-X(2)))**2))*alpha*tanphi+COS(6.0*X(1))*(1.0- &
                                  & (TANH(-0.1E1+alpha*(tanphi*X(1)-X(2)))**2))*alpha))
                                CASE DEFAULT
                                localError="The analytic function type of "// &
                                  & TRIM(NumberToVString(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
                                  & " is invalid."
                                CALL FlagError(localError,err,error,*999)
                              END SELECT
                              !Default to version 1 of each node derivative
                              local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                              CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                                & FIELD_VALUES_SET_TYPE,local_ny,VALUE_SOURCE,err,error,*999)
                            ENDDO !deriv_idx
                          ENDDO !node_idx
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
                ENDDO !component_idx
                CALL FIELD_PARAMETER_SET_UPDATE_START(SOURCE_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & err,error,*999)
                CALL FIELD_PARAMETER_SET_UPDATE_FINISH(SOURCE_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & err,error,*999)
              ELSE
                CALL FlagError("Field variable is not associated.",err,error,*999)
              ENDIF
            ENDDO !variable_idx
            CALL FIELD_PARAMETER_SET_DATA_RESTORE(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & GEOMETRIC_PARAMETERS,err,error,*999)
          ELSE
            CALL FlagError("Equations set geometric field is not associated.",err,error,*999)
          ENDIF            
        ELSE
          CALL FlagError("Equations set source field is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Equations set analytic is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

    !>Set the material field to a specified analytical value
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
        MATERIALS_FIELD=>EQUATIONS_SET%MATERIALS%MATERIALS_FIELD
        IF(ASSOCIATED(MATERIALS_FIELD)) THEN
          GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
          IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN            
            CALL FIELD_NUMBER_OF_COMPONENTS_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS,err,error,*999)
            NULLIFY(GEOMETRIC_VARIABLE)
            CALL Field_VariableGet(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,GEOMETRIC_VARIABLE,err,error,*999)
            CALL FIELD_PARAMETER_SET_DATA_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS, &
              & err,error,*999)
            DO variable_idx=1,MATERIALS_FIELD%NUMBER_OF_VARIABLES
              variable_type=MATERIALS_FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE
              FIELD_VARIABLE=>MATERIALS_FIELD%VARIABLE_TYPE_MAP(variable_type)%ptr
              IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                   SELECT CASE(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE)
                   CASE(EQUATIONS_SET_ADVECTION_DIFFUSION_EQUATION_TWO_DIM_1)
                     VALUE_MATERIAL= (1.0/Peclet)
                   CASE DEFAULT
                     localError="The analytic function type of "// &
                      & TRIM(NumberToVString(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
                      & " is invalid."
                      CALL FlagError(localError,err,error,*999)
                   END SELECT
                      CALL FIELD_PARAMETER_SET_UPDATE_CONSTANT(MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                       & FIELD_VALUES_SET_TYPE,component_idx,VALUE_MATERIAL,err,error,*999)
                ENDDO !component_idx
                CALL FIELD_PARAMETER_SET_UPDATE_START(MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & err,error,*999)
                CALL FIELD_PARAMETER_SET_UPDATE_FINISH(MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & err,error,*999)
              ELSE
                CALL FlagError("Field variable is not associated.",err,error,*999)
              ENDIF
            ENDDO !variable_idx
            CALL FIELD_PARAMETER_SET_DATA_RESTORE(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & GEOMETRIC_PARAMETERS,err,error,*999)
          ELSE
            CALL FlagError("Equations set geometric field is not associated.",err,error,*999)
          ENDIF            
        ELSE
          CALL FlagError("Equations set material field is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Equations set analytic is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF


    EXITS("AdvectionDiffusion_BoundaryConditionsAnalyticCalculate")
    RETURN
999 ERRORS("AdvectionDiffusion_BoundaryConditionsAnalyticCalculate",err,error)
    EXITS("AdvectionDiffusion_BoundaryConditionsAnalyticCalculate")
    RETURN 1
    
  END SUBROUTINE AdvectionDiffusion_BoundaryConditionsAnalyticCalculate


  !
  !================================================================================================================================
  !

  !>Sets up the diffusion equation type of a classical field equations set class.
  SUBROUTINE AdvectionDiffusion_EquationsSetSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup a diffusion equation on.
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("AdvectionDiffusion_EquationsSetSetup",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<3) THEN
        CALL FlagError("Equations set specification does not have a subtype set.",err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
      CASE(EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE)
        CALL AdvectionDiffusion_EquationsSetLinearSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*999)
      CASE(EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE)
        !Need to define the functions diffusion_equation_equations_set_linear_source_setup etc
        CALL AdvectionDiffusion_EquationsSetLinearSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*999)
      CASE(EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE)
        CALL AdvectionDiffusion_EquationsSetLinearSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*999)
      CASE(EQUATIONS_SET_QUADRATIC_SOURCE_ADVECTION_DIFFUSION_SUBTYPE)
         CALL FlagError("Not implemented.",err,error,*999)
!        CALL ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*999)
      CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_ADVECTION_DIFFUSION_SUBTYPE)
         CALL FlagError("Not implemented.",err,error,*999)
!        CALL ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*999)
      CASE(EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
        CALL AdvectionDiffusion_EquationsSetLinearSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*999)
      CASE(EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
        !Need to define the functions diffusion_equation_equations_set_linear_source_setup etc
        CALL AdvectionDiffusion_EquationsSetLinearSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*999)
      CASE(EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
        CALL AdvectionDiffusion_EquationsSetLinearSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*999)
      CASE(EQUATIONS_SET_QUADRATIC_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
         CALL FlagError("Not implemented.",err,error,*999)
!        CALL ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*999)
      CASE(EQUATIONS_SET_EXP_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
         CALL FlagError("Not implemented.",err,error,*999)
!        CALL ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*999)
      CASE(EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE)
        CALL AdvectionDiffusion_EquationsSetLinearSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*999)
      CASE(EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE)
        !Need to define the functions diffusion_equation_equations_set_linear_source_setup etc
        CALL AdvectionDiffusion_EquationsSetLinearSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*999)
      CASE(EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE)
        CALL AdvectionDiffusion_EquationsSetLinearSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*999)
      CASE(EQUATIONS_SET_QUAD_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE)
         CALL FlagError("Not implemented.",err,error,*999)
!        CALL ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*999)
      CASE(EQUATIONS_SET_EXP_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE)
         CALL FlagError("Not implemented.",err,error,*999)
!        CALL ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*999)
      CASE(EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE)
        CALL AdvectionDiffusion_EquationsSetLinearSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*999)
      CASE(EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE)
        !Need to define the functions diffusion_equation_equations_set_linear_source_setup etc
        CALL AdvectionDiffusion_EquationsSetLinearSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*999)
      CASE(EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE)
        CALL AdvectionDiffusion_EquationsSetLinearSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*999)
      CASE(EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE)
        CALL AdvectionDiffusion_EquationsSetLinearSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*999)
      CASE(EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE)
        !Need to define the functions diffusion_equation_equations_set_linear_source_setup etc
        CALL AdvectionDiffusion_EquationsSetLinearSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*999)
      CASE(EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE)
        CALL AdvectionDiffusion_EquationsSetLinearSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*999)
      CASE(EQUATIONS_SET_QUAD_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE)
         CALL FlagError("Not implemented.",err,error,*999)
!        CALL ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*999)
      CASE(EQUATIONS_SET_EXP_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE)
         CALL FlagError("Not implemented.",err,error,*999)
!        CALL ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*999)
      CASE(EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE)
        CALL AdvectionDiffusion_EquationsSetLinearSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*999)
      CASE(EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE)
        !Need to define the functions diffusion_equation_equations_set_linear_source_setup etc
        CALL AdvectionDiffusion_EquationsSetLinearSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*999)
      CASE(EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE)
        CALL AdvectionDiffusion_EquationsSetLinearSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*999)
      CASE(EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
        CALL AdvectionDiffusion_EquationsSetLinearSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*999)
      CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE)
        CALL AdvectionDiffusion_EquationsSetLinearSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*999)
      CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)
        CALL AdvectionDiffusion_EquationsSetLinearSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*999)
      CASE DEFAULT
        localError="Equations set subtype "//TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
          & " is not valid for an advection-diffusion equation type of a classical field equation set class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
       
    EXITS("AdvectionDiffusion_EquationsSetSetup")
    RETURN
999 ERRORS("AdvectionDiffusion_EquationsSetSetup",err,error)
    EXITS("AdvectionDiffusion_EquationsSetSetup")
    RETURN 1
    
  END SUBROUTINE AdvectionDiffusion_EquationsSetSetup

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a diffusion equation type of an classical field equations set class.
  SUBROUTINE AdvectionDiffusion_EquationsSetSolnMethodSet(EQUATIONS_SET,SOLUTION_METHOD,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("AdvectionDiffusion_EquationsSetSolnMethodSet",err,error,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<3) THEN
        CALL FlagError("Equations set specification does not have a subtype set.",err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
      CASE(EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE,EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
         & EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE,EQUATIONS_SET_QUADRATIC_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
         & EQUATIONS_SET_EXPONENTIAL_SOURCE_ADVECTION_DIFFUSION_SUBTYPE,&
         & EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE,&
         & EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
         & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE,&
         & EQUATIONS_SET_QUADRATIC_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
         & EQUATIONS_SET_EXP_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE,&
         & EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE,&
         & EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
         & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE,&
         & EQUATIONS_SET_QUAD_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
         & EQUATIONS_SET_EXP_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE,&
         & EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE, &
         & EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE,EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE, &
         & EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE,EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
         & EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE,EQUATIONS_SET_QUAD_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
         & EQUATIONS_SET_EXP_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE,EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
         & EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, & 
         & EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
         & EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE,EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE, &
         & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)        
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
          localError="The specified solution method of "//TRIM(NumberToVString(SOLUTION_METHOD,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
     CASE DEFAULT
        localError="Equations set subtype of "//TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
          & " is not valid for an advection-diffusion equation type of an classical field equations set class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
       
    EXITS("AdvectionDiffusion_EquationsSetSolnMethodSet")
    RETURN
999 ERRORS("AdvectionDiffusion_EquationsSetSolnMethodSet",err,error)
    EXITS("AdvectionDiffusion_EquationsSetSolnMethodSet")
    RETURN 1
  END SUBROUTINE AdvectionDiffusion_EquationsSetSolnMethodSet

  !
  !================================================================================================================================
  !

  !>Sets the equation specification for a diffusion equation type of a classical field equations set class.
  SUBROUTINE AdvectionDiffusion_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to set the specification for
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations set specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: subtype

    ENTERS("AdvectionDiffusion_EquationsSetSpecificationSet",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(SIZE(specification,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for an advection-diffusion type equations set.", &
          & err,error,*999)
      END IF
      subtype=specification(3)
      SELECT CASE(subtype)
      CASE(EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
          & EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
          & EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
          & EQUATIONS_SET_QUADRATIC_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
          & EQUATIONS_SET_EXPONENTIAL_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
          & EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
          & EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
          & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
          & EQUATIONS_SET_QUADRATIC_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
          & EQUATIONS_SET_EXP_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
          & EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
          & EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
          & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
          & EQUATIONS_SET_QUAD_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
          & EQUATIONS_SET_EXP_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
          & EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE, &
          & EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE, &
          & EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE, &
          & EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
          & EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
          & EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
          & EQUATIONS_SET_QUAD_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
          & EQUATIONS_SET_EXP_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
          & EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
          & EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
          & EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
          & EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE, &
          & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE, &
          & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)
        !ok
      CASE DEFAULT
        localError="The third equations set specification of "//TRIM(NumberToVstring(subtype,"*",err,error))// &
          & " is not valid for an advection-diffusion type of a classical field equations set."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Set full specification
      IF(ALLOCATED(equationsSet%specification)) THEN
        CALL FlagError("Equations set specification is already allocated.",err,error,*999)
      ELSE
        ALLOCATE(equationsSet%specification(3),stat=err)
        IF(err/=0) CALL FlagError("Could not allocate equations set specification.",err,error,*999)
      END IF
      equationsSet%specification(1:3)=[EQUATIONS_SET_CLASSICAL_FIELD_CLASS,EQUATIONS_SET_ADVECTION_DIFFUSION_EQUATION_TYPE,subtype]
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    END IF

    EXITS("AdvectionDiffusion_EquationsSetSpecificationSet")
    RETURN
999 ERRORS("AdvectionDiffusion_EquationsSetSpecificationSet",err,error)
    EXITS("AdvectionDiffusion_EquationsSetSpecificationSet")
    RETURN 1
    
  END SUBROUTINE AdvectionDiffusion_EquationsSetSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the linear advection-diffusion equation.
  SUBROUTINE AdvectionDiffusion_EquationsSetLinearSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,GEOMETRIC_MESH_COMPONENT,GEOMETRIC_SCALING_TYPE,NUMBER_OF_DIMENSIONS, &
      & NUMBER_OF_MATERIALS_COMPONENTS, NUMBER_OF_SOURCE_COMPONENTS, NUMBER_OF_INDEPENDENT_COMPONENTS,imy_matrix,Ncompartments,&
      & GEOMETRIC_COMPONENT_NUMBER,NUMBER_OF_INDEPENDENT_U_VAR_COMPONENTS,NUMBER_OF_INDEPENDENT_V_VAR_COMPONENTS
    TYPE(DECOMPOSITION_TYPE), POINTER :: GEOMETRIC_DECOMPOSITION
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EQUATIONS_SET_EQUATIONS_SET_FIELD_TYPE), POINTER :: EQUATIONS_EQUATIONS_SET_FIELD
    TYPE(EQUATIONS_SET_MATERIALS_TYPE), POINTER :: EQUATIONS_MATERIALS
    TYPE(EQUATIONS_SET_SOURCE_TYPE), POINTER :: EQUATIONS_SOURCE
    TYPE(EQUATIONS_SET_INDEPENDENT_TYPE), POINTER :: EQUATIONS_INDEPENDENT
    TYPE(FIELD_TYPE), POINTER :: ANALYTIC_FIELD,DEPENDENT_FIELD,GEOMETRIC_FIELD,EQUATIONS_SET_FIELD_FIELD
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: num_var,num_var_count,NUMBER_OF_MATERIALS_COUPLING_COMPONENTS    
    INTEGER(INTG) :: EQUATIONS_SET_FIELD_NUMBER_OF_VARIABLES,EQUATIONS_SET_FIELD_NUMBER_OF_COMPONENTS    
    INTEGER(INTG), POINTER :: EQUATIONS_SET_FIELD_DATA(:)
    INTEGER(INTG), ALLOCATABLE :: VARIABLE_TYPES(:),VARIABLE_U_TYPES(:),COUPLING_MATRIX_STORAGE_TYPE(:), &
      & COUPLING_MATRIX_STRUCTURE_TYPE(:)
    INTEGER(INTG) :: EQUATIONS_SET_SUBTYPE
    
    ENTERS("AdvectionDiffusion_EquationsSetLinearSetup",err,error,*999)

    NULLIFY(EQUATIONS)
    NULLIFY(vectorMapping)
    NULLIFY(vectorMatrices)
    NULLIFY(GEOMETRIC_DECOMPOSITION)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<3) THEN
        CALL FlagError("Equations set specification does not have a subtype set.",err,error,*999)
      END IF
      EQUATIONS_SET_SUBTYPE=EQUATIONS_SET%SPECIFICATION(3)
      IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
        & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
        & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
        & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
        & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
        & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
        & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
        & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
        & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
        & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
        & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
        & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
        & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
        & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
        & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
        & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
        & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
        & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
        & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE .OR. &
        & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
        & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
        SELECT CASE(EQUATIONS_SET_SETUP%SETUP_TYPE)
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            CALL AdvectionDiffusion_EquationsSetSolnMethodSet(EQUATIONS_SET, &
              & EQUATIONS_SET_FEM_SOLUTION_METHOD,err,error,*999)
             IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
              & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
              EQUATIONS_SET_FIELD_NUMBER_OF_VARIABLES = 1
              EQUATIONS_SET_FIELD_NUMBER_OF_COMPONENTS = 2
              EQUATIONS_EQUATIONS_SET_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD
              IF(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_AUTO_CREATED) THEN
                !Create the auto created equations set field
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION, &
                  & EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,err,error,*999)
                CALL FIELD_LABEL_SET(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,"Equations Set Field",err,error,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,FIELD_GENERAL_TYPE,&
                  & err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,&
                  & FIELD_INDEPENDENT_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_SET(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                  & EQUATIONS_SET_FIELD_NUMBER_OF_VARIABLES,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,&
                  & [FIELD_U_VARIABLE_TYPE],err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_INTG_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,&
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
             IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
              & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
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
            localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a linear advection-diffusion equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
          !Do nothing???
            SELECT CASE(EQUATIONS_SET_SUBTYPE)
            CASE(EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
               & EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
               & EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE)
            CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE,EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)
              SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
              CASE(EQUATIONS_SET_SETUP_START_ACTION)
                EQUATIONS_SET_FIELD_NUMBER_OF_COMPONENTS = 2
                EQUATIONS_EQUATIONS_SET_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD
                IF(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_AUTO_CREATED) THEN
                  CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,err,error,*999)
                  CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,&
                    & GEOMETRIC_DECOMPOSITION,err,error,*999)
                  CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,& 
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
              CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
                ! do nothing
              CASE DEFAULT
                localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
                  & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                  & " is invalid for a linear advection-diffusion equation."
                CALL FlagError(localError,err,error,*999)
              END SELECT
               !Do nothing???
            CASE(EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
               & EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
               & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
               & EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE)
               CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD, FIELD_U_VARIABLE_TYPE, &
                  & FIELD_MESH_DISPLACEMENT_SET_TYPE, ERR, ERROR, *999)
               CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD, FIELD_U_VARIABLE_TYPE, &
                  & FIELD_MESH_VELOCITY_SET_TYPE, ERR, ERROR, *999)
            END SELECT
        !-----------------------------------------------------------------
        ! D e p e n d e n t   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            SELECT CASE(EQUATIONS_SET_SUBTYPE)
            CASE(EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
               & EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
               & EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
               & EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
               & EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
               & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
               & EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE, &
               & EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE, &
               & EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE, &
               & EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE, &
               & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)
               IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                 !Create the auto created dependent field
                 CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET%DEPENDENT% &
                   & DEPENDENT_FIELD,err,error,*999)
                 CALL FIELD_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,"Dependent Field",err,error,*999)
                 CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                 CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
                 CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,err,error,*999)
                 CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_DECOMPOSITION, &
                   & err,error,*999)
                 CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY% &
                   & GEOMETRIC_FIELD,err,error,*999)
                 CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,2,err,error,*999)
                 CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,[FIELD_U_VARIABLE_TYPE, &
                   & FIELD_DELUDELN_VARIABLE_TYPE],err,error,*999)
                 CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
                 CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                   & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
                 CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & FIELD_DP_TYPE,err,error,*999)
                 CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                   & FIELD_DP_TYPE,err,error,*999)
                 CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                   & err,error,*999)
                 CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                   & FIELD_DELUDELN_VARIABLE_TYPE,1,err,error,*999)
                 !Default to the geometric interpolation setup
                 CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                   & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                 CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                   & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                 CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                   & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                 SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                 CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                   CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                     & FIELD_U_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                   CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                     & FIELD_DELUDELN_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
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
                   localError="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                     & " is invalid."
                   CALL FlagError(localError,err,error,*999)
                 END SELECT
               ELSE
                SELECT CASE(EQUATIONS_SET_SUBTYPE)
                CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE, &
                   & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)
                !uses number of compartments to check that appropriate number and type of variables have been set on the dependent field
                  CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                  CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)                 
                  EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
                  CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & FIELD_VALUES_SET_TYPE,EQUATIONS_SET_FIELD_DATA,err,error,*999)
                  Ncompartments=EQUATIONS_SET_FIELD_DATA(2)
                  CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,2*Ncompartments,err,error,*999)
                  !Create & populate array storing all of the relevant variable types against which to check the field variables
                  ALLOCATE(VARIABLE_TYPES(2*Ncompartments))
                  DO num_var=1,Ncompartments
                    VARIABLE_TYPES(2*num_var-1)=FIELD_U_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(num_var-1))
                    VARIABLE_TYPES(2*num_var)=FIELD_DELUDELN_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(num_var-1))
                  ENDDO
                  CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,VARIABLE_TYPES,err,error,*999)

                  DO num_var=1,2*Ncompartments
                    CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,VARIABLE_TYPES(num_var), & 
                       & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
                    CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,VARIABLE_TYPES(num_var),FIELD_DP_TYPE,err,error,*999)
                    CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,VARIABLE_TYPES(num_var),1, &
                       & err,error,*999)
                  ENDDO
                  CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_DIMENSIONS,err,error,*999)
                  SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                  CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                      component_idx=1
                    DO num_var=1,2*Ncompartments
                      CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,VARIABLE_TYPES(num_var),component_idx, &
                        & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
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
                    localError="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                      & " is invalid."
                    CALL FlagError(localError,err,error,*999)
                  END SELECT
                CASE DEFAULT
                 !Check the user specified field
                 CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                 CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
                 CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,2,err,error,*999)
                 CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE], &
                   & err,error,*999)
                 CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, & 
                   & err,error,*999)
                 CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                   & err,error,*999)
                 CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                 CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                 CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & NUMBER_OF_DIMENSIONS,err,error,*999)
                 CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS, &
                   & err,error,*999)
                 CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                   & NUMBER_OF_DIMENSIONS,err,error,*999)
                 SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                 CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
!                   DO component_idx=1,NUMBER_OF_DIMENSIONS
                     component_idx=1
                     CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,component_idx, &
                       & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                     CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, & 
                       & component_idx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
!                   ENDDO !component_idx
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
                   localError="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                     & " is invalid."
                   CALL FlagError(localError,err,error,*999)
                 END SELECT
                END SELECT 
               ENDIF
            CASE(EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)

               IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                 !Create the auto created dependent field
                 CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET%DEPENDENT% &
                   & DEPENDENT_FIELD,err,error,*999)
                 CALL FIELD_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,"Dependent Field",err,error,*999)
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
                  & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELVDELN_VARIABLE_TYPE, &
                  & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELVDELN_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)

                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                  & err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,1,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & 1,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELVDELN_VARIABLE_TYPE, &
                  & 1,err,error,*999)
                 !Default to the geometric interpolation setup
                 CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                   & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                 CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                   & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                 CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                   & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                 CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                   & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                 CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE,1, &
                   & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                 CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELVDELN_VARIABLE_TYPE,1, &
                   & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                 SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                 CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                   CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                     & FIELD_U_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                   CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                     & FIELD_DELUDELN_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                   CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                     & FIELD_V_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                   CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                     & FIELD_DELVDELN_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
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
                   localError="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                     & " is invalid."
                   CALL FlagError(localError,err,error,*999)
                 END SELECT
               ELSE
                 !Check the user specified field
                 CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                 CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
                 CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,4,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE,&
                  & FIELD_V_VARIABLE_TYPE,FIELD_DELVDELN_VARIABLE_TYPE],err,error,*999)

                 CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, & 
                   & err,error,*999)
                 CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                   & err,error,*999)
                 CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, & 
                   & err,error,*999)
                 CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELVDELN_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                   & err,error,*999)
                 CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                 CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                 CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                 CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELVDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                 CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & NUMBER_OF_DIMENSIONS,err,error,*999)
                 CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                   & err,error,*999)
                 CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, & 
                   & 1,err,error,*999)
                 CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,1, &
                   & err,error,*999)
                 CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELVDELN_VARIABLE_TYPE, & 
                   & 1,err,error,*999)
                 SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                 CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
 !                  DO component_idx=1,NUMBER_OF_DIMENSIONS
                     component_idx=1
                     CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,component_idx, &
                       & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                     CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, & 
                       & component_idx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                     CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,component_idx, &
                       & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                     CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELVDELN_VARIABLE_TYPE, & 
                       & component_idx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
!                   ENDDO !component_idx
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
                   localError="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                     & " is invalid."
                   CALL FlagError(localError,err,error,*999)
                 END SELECT
               ENDIF
              END SELECT
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
              CALL FIELD_CREATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,err,error,*999)
              CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                 & FIELD_BOUNDARY_CONDITIONS_SET_TYPE,err,error,*999)
            ENDIF
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a linear advection-diffusion equation"
            CALL FlagError(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! M a t e r i a l s   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
            IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
              IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
               IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
                & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)THEN
                  CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_MATERIALS% &
                    & MATERIALS_FIELD,err,error,*999)
                  CALL FIELD_LABEL_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,"Materials Field",err,error,*999)
                  CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_MATERIAL_TYPE,err,error,*999)
                  CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                  CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,err,error,*999)
                  CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,GEOMETRIC_DECOMPOSITION, &
                    & err,error,*999)
                  CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,EQUATIONS_SET%GEOMETRY% &
                    & GEOMETRIC_FIELD,err,error,*999)
                  CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,2,err,error,*999)
                  CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,[FIELD_U_VARIABLE_TYPE, &
                    & FIELD_V_VARIABLE_TYPE], &
                    & err,error,*999)
                  CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                  CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                  CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_DIMENSIONS,err,error,*999)
                    NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS
                   !Set the number of materials components
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_MATERIALS_COMPONENTS,err,error,*999)
                  EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
                  CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & FIELD_VALUES_SET_TYPE,EQUATIONS_SET_FIELD_DATA,err,error,*999)
                    Ncompartments=EQUATIONS_SET_FIELD_DATA(2)
                    NUMBER_OF_MATERIALS_COUPLING_COMPONENTS=Ncompartments
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & NUMBER_OF_MATERIALS_COUPLING_COMPONENTS,err,error,*999)
                  !Default the k materials components to the geometric interpolation setup with constant interpolation
                  DO component_idx=1,NUMBER_OF_DIMENSIONS
                    CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & component_idx,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  ENDDO !component_idx
                    CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & 1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  DO component_idx=1,NUMBER_OF_MATERIALS_COUPLING_COMPONENTS
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                      & component_idx,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                      & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  ENDDO 
                    !Default the field scaling to that of the geometric field
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
               ELSE !standard materials field
                !Create the auto created materials field
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_MATERIALS% &
                  & MATERIALS_FIELD,err,error,*999)
                CALL FIELD_LABEL_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,"Materials Field",err,error,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_MATERIAL_TYPE,err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,err,error,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & err,error,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,EQUATIONS_SET%GEOMETRY% &
                  & GEOMETRIC_FIELD,err,error,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,1,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,[FIELD_U_VARIABLE_TYPE], &
                  & err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,err,error,*999)
                IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN
                  NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS
                ELSEIF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN
                  !Linear source. Materials field components are 1 for each dimension and 1 for the linear source
                  !i.e., k and a in div(k.grad(u(x)))=a(x)u(x)+c(x)
                  NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS+1
                ELSE
                  NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS
                ENDIF
                 !Set the number of materials components
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_MATERIALS_COMPONENTS,err,error,*999)
                !Default the k materials components to the geometric interpolation setup with constant interpolation
                DO component_idx=1,NUMBER_OF_DIMENSIONS
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                ENDDO !component_idx
                IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE) THEN
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 1,GEOMETRIC_MESH_COMPONENT,err,error,*999)   
                  DO component_idx=NUMBER_OF_DIMENSIONS+1,NUMBER_OF_MATERIALS_COMPONENTS
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & component_idx,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  ENDDO !component_idx
                ENDIF
                  !Default the field scaling to that of the geometric field
                CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                CALL FIELD_SCALING_TYPE_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
               ENDIF
              ELSE
               IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
                & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)THEN
                  CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_MATERIAL_TYPE,err,error,*999)
                  CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                  CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,2,err,error,*999)
                  CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE, &
                     & FIELD_V_VARIABLE_TYPE],err,error,*999)
                  CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                    & err,error,*999)
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                    & err,error,*999)
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_DIMENSIONS,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS, &
                    & err,error,*999)
                  EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
                  CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & FIELD_VALUES_SET_TYPE,EQUATIONS_SET_FIELD_DATA,err,error,*999)
                    Ncompartments=EQUATIONS_SET_FIELD_DATA(2)
                    NUMBER_OF_MATERIALS_COUPLING_COMPONENTS=Ncompartments
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE, &
                    & NUMBER_OF_MATERIALS_COUPLING_COMPONENTS,err,error,*999)
               ELSE
                !Check the user specified field
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_MATERIAL_TYPE,err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,1,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,err,error,*999)
                IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. & 
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS, &
                    & err,error,*999)
                ELSEIF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS+1, &
                    & err,error,*999)
                ENDIF
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
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,err,error,*999)
                IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
                  NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS             
                ELSEIF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. & 
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. & 
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN
                  !Linear source. Materials field components are 1 for each dimension and 1 for the linear source
                  !i.e., k and a in div(k.grad(u(x)))=a(x)u(x)+c(x)
                  NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS+1
                ELSE
                  NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS
                ENDIF
                !First set the k values to 1.0
                DO component_idx=1,NUMBER_OF_DIMENSIONS
                  !WRITE(*,'("Setting materials components values :")')
                  CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,component_idx,1.0_DP,err,error,*999)
                ENDDO !component_idx
                IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
                & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
                  EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
                  CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & FIELD_VALUES_SET_TYPE,EQUATIONS_SET_FIELD_DATA,err,error,*999)
                    Ncompartments=EQUATIONS_SET_FIELD_DATA(2)
                  DO component_idx=1,Ncompartments
                   CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                     & FIELD_VALUES_SET_TYPE,component_idx,0.0_DP,err,error,*999)
                   ENDDO !component_idx
                ENDIF
                IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. & 
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. & 
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN
                  !Now set the linear source values to 1.0
                  DO component_idx=NUMBER_OF_DIMENSIONS+1,NUMBER_OF_MATERIALS_COMPONENTS
                    CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,component_idx,1.0_DP,err,error,*999)
                  ENDDO !component_idx
                ENDIF
              ENDIF
            ELSE
              CALL FlagError("Equations set materials is not associated.",err,error,*999)
            ENDIF
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a linear advection-diffusion equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! S o u r c e   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
              EQUATIONS_SOURCE=>EQUATIONS_SET%SOURCE
            IF(ASSOCIATED(EQUATIONS_SOURCE)) THEN
              IF(EQUATIONS_SOURCE%SOURCE_FIELD_AUTO_CREATED) THEN
                !Create the auto created source field
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SOURCE% &
                  & SOURCE_FIELD,err,error,*999)
                CALL FIELD_LABEL_SET(EQUATIONS_SOURCE%SOURCE_FIELD,"Source Field",err,error,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,err,error,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & err,error,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,EQUATIONS_SET%GEOMETRY% &
                  & GEOMETRIC_FIELD,err,error,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,1,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,[FIELD_U_VARIABLE_TYPE], &
                  & err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                NUMBER_OF_SOURCE_COMPONENTS=1
                !Set the number of source components
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_SOURCE_COMPONENTS,err,error,*999)
                !Default the source components to the geometric interpolation setup with constant interpolation
                IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &  
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. & 
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
                  DO component_idx=1,NUMBER_OF_SOURCE_COMPONENTS
                    CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                     & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & component_idx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  ENDDO !component_idx
                ENDIF
                  !Default the field scaling to that of the geometric field
                CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SOURCE%SOURCE_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
              ELSE
                !Check the user specified field
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,1,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                  & err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Equations set source is not associated.",err,error,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            EQUATIONS_SOURCE=>EQUATIONS_SET%SOURCE
            IF(ASSOCIATED(EQUATIONS_SOURCE)) THEN
              IF(EQUATIONS_SOURCE%SOURCE_FIELD_AUTO_CREATED) THEN
                !Finish creating the source field
                CALL FIELD_CREATE_FINISH(EQUATIONS_SOURCE%SOURCE_FIELD,err,error,*999)
                !Set the default values for the source field
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,err,error,*999)
                IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. & 
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. & 
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
                    & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
                  NUMBER_OF_SOURCE_COMPONENTS=1
                ELSE
                  NUMBER_OF_SOURCE_COMPONENTS=0
                ENDIF
                !Now set the source values to 1.0
                DO component_idx=1,NUMBER_OF_SOURCE_COMPONENTS
                  CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,component_idx,1.0_DP,err,error,*999)
                ENDDO !component_idx
              ENDIF
            ELSE
              CALL FlagError("Equations set source is not associated.",err,error,*999)
            ENDIF
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a linear advection-diffusion equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! I n d e p e n d e n t   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_INDEPENDENT_TYPE)
          !Setup the equations set for the advective velocity field
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            EQUATIONS_INDEPENDENT=>EQUATIONS_SET%INDEPENDENT
            IF(ASSOCIATED(EQUATIONS_INDEPENDENT)) THEN
             IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
              & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
              IF(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
                !Create the auto created independent field
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_INDEPENDENT% &
                  & INDEPENDENT_FIELD,err,error,*999)
                CALL FIELD_LABEL_SET(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,"Independent Field",err,error,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD, &
                      & FIELD_INDEPENDENT_TYPE,err,error,*999)
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,err,error,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & err,error,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY% &
                  & GEOMETRIC_FIELD,err,error,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,2,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,[FIELD_U_VARIABLE_TYPE, &
                  & FIELD_V_VARIABLE_TYPE],err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,err,error,*999)
                  NUMBER_OF_INDEPENDENT_U_VAR_COMPONENTS=NUMBER_OF_DIMENSIONS
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_INDEPENDENT_U_VAR_COMPONENTS,err,error,*999)
                EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
                CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & FIELD_VALUES_SET_TYPE,EQUATIONS_SET_FIELD_DATA,err,error,*999)
                  Ncompartments=EQUATIONS_SET_FIELD_DATA(2)
                  NUMBER_OF_INDEPENDENT_V_VAR_COMPONENTS=Ncompartments-1
                 !Set the number of independent components
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & NUMBER_OF_INDEPENDENT_V_VAR_COMPONENTS,err,error,*999)
                !Default the k independent components to the geometric interpolation setup with constant interpolation
                DO component_idx=1,NUMBER_OF_INDEPENDENT_U_VAR_COMPONENTS
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,component_idx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                ENDDO !component_idx
                DO component_idx=1,NUMBER_OF_INDEPENDENT_V_VAR_COMPONENTS
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                    & FIELD_V_VARIABLE_TYPE,component_idx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                ENDDO !component_idx
                  !Default the field scaling to that of the geometric field
                CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                CALL FIELD_SCALING_TYPE_SET(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
              ELSE
                !Check the user specified field
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,2,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE],&
                   & err,error,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                   & err,error,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                   & err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS, &
                  & err,error,*999)
                EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
                CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & FIELD_VALUES_SET_TYPE,EQUATIONS_SET_FIELD_DATA,err,error,*999)
                  Ncompartments=EQUATIONS_SET_FIELD_DATA(2)
                  NUMBER_OF_INDEPENDENT_V_VAR_COMPONENTS=Ncompartments-1
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_INDEPENDENT_V_VAR_COMPONENTS,err,error,*999)
                DO component_idx=1,NUMBER_OF_INDEPENDENT_U_VAR_COMPONENTS
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                ENDDO
                DO component_idx=1,NUMBER_OF_INDEPENDENT_V_VAR_COMPONENTS
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,component_idx, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                ENDDO
              ENDIF
             ELSE
              IF(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
                !Create the auto created independent field
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_INDEPENDENT% &
                  & INDEPENDENT_FIELD,err,error,*999)
                CALL FIELD_LABEL_SET(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,"Independent Field",err,error,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD, &
                      & FIELD_INDEPENDENT_TYPE,err,error,*999)
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,err,error,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & err,error,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY% &
                  & GEOMETRIC_FIELD,err,error,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,1,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,[FIELD_U_VARIABLE_TYPE], &
                  & err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,err,error,*999)
                  NUMBER_OF_INDEPENDENT_COMPONENTS=NUMBER_OF_DIMENSIONS
                 !Set the number of independent components
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_INDEPENDENT_COMPONENTS,err,error,*999)
                !Default the k independent components to the geometric interpolation setup with constant interpolation
                DO component_idx=1,NUMBER_OF_DIMENSIONS
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,component_idx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                ENDDO !component_idx
                  !Default the field scaling to that of the geometric field
                CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                CALL FIELD_SCALING_TYPE_SET(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
              ELSE
                !Check the user specified field
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,1,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS, &
                  & err,error,*999)
                CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              ENDIF
             ENDIF
            ELSE
              CALL FlagError("Equations set independent is not associated.",err,error,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            EQUATIONS_INDEPENDENT=>EQUATIONS_SET%INDEPENDENT
            IF(ASSOCIATED(EQUATIONS_INDEPENDENT)) THEN
              IF(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
                !Finish creating the independent field
                CALL FIELD_CREATE_FINISH(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,err,error,*999)
                !Set the default values for the independent field
                CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_INPUT_DATA1_SET_TYPE,err,error,*999)

                 IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
                  & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
!                   CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
!                      & FIELD_INPUT_DATA2_SET_TYPE,err,error,*999)
                 ENDIF
!                 CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
!                   & NUMBER_OF_DIMENSIONS,err,error,*999)
!                   NUMBER_OF_INDEPENDENT_COMPONENTS=NUMBER_OF_DIMENSIONS             
!                 !First set the k values to 1.0
!                 DO component_idx=1,NUMBER_OF_INDEPENDENT_COMPONENTS
!                   CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
!                     & FIELD_VALUES_SET_TYPE,component_idx,1.0_DP,err,error,*999)
!                ENDDO !component_idx
              ENDIF
            ELSE
              CALL FlagError("Equations set independent is not associated.",err,error,*999)
            ENDIF
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a linear advection-diffusion equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! A n a l y t i c   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
              DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
              IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN
                  CALL FIELD_NUMBER_OF_COMPONENTS_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS,err,error,*999)
                  SELECT CASE(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE)
                  CASE(EQUATIONS_SET_ADVECTION_DIFFUSION_EQUATION_TWO_DIM_1)
                      IF(NUMBER_OF_DIMENSIONS/=2) THEN
                        localError="The number of geometric dimensions of "// &
                          & TRIM(NumberToVString(NUMBER_OF_DIMENSIONS,"*",err,error))// &
                          & " is invalid. The analytic function type of "// &
                          & TRIM(NumberToVString(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
                          & " requires that there be 2 geometric dimensions."
                        CALL FlagError(localError,err,error,*999)
                      ENDIF
                      EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE= & 
                       & EQUATIONS_SET_ADVECTION_DIFFUSION_EQUATION_TWO_DIM_1
                  CASE DEFAULT
                    localError="The specified analytic function type of "// &
                      & TRIM(NumberToVString(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
                      & " is invalid for a linear advection-diffusion equation."
                    CALL FlagError(localError,err,error,*999)
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
            localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a linear advection-diffusion equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! E q u a t i o n s    t y p e
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
              CALL Equations_CreateStart(EQUATIONS_SET,EQUATIONS,err,error,*999)
              CALL Equations_LinearityTypeSet(EQUATIONS,EQUATIONS_LINEAR,err,error,*999)
              IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE .OR. & 
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN 
                  CALL Equations_TimeDependenceTypeSet(EQUATIONS,EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
              ELSEIF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN
                  CALL Equations_TimeDependenceTypeSet(EQUATIONS,EQUATIONS_STATIC,err,error,*999)
              ELSE
                CALL FlagError("Equations set subtype not valid.",err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Equations set dependent field has not been finished.",err,error,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
            IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
               & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
               & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
               & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
               & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
               & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. & 
               & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
               & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
               & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
               & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
               & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
               & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
               & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE .OR. & 
               & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
               & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN 
              !Finish the equations
              CALL EquationsSet_EquationsGet(EQUATIONS_SET,EQUATIONS,err,error,*999)
              CALL Equations_CreateFinish(EQUATIONS,err,error,*999)
              NULLIFY(vectorEquations)
              CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
              !Create the equations mapping.
              SELECT CASE(EQUATIONS_SET_SUBTYPE)
              CASE(EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE,EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE,&
               & EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE,EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE,&
               & EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE,&
               & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE,&
               & EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE,&
               & EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE,&
               & EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE,&
               & EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE,&
               & EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE,&
               & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE,&
               & EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
              CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_DELUDELN_VARIABLE_TYPE,vectorMapping,err,error,*999)
              CALL EquationsMapping_DynamicMatricesSet(vectorMapping,.TRUE.,.TRUE.,err,error,*999)
              CALL EquationsMapping_DynamicVariableTypeSet(vectorMapping,FIELD_U_VARIABLE_TYPE,err,error,*999)
              CALL EquationsMapping_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE,err,error,*999)
                IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE) THEN             
                CALL EquationsMapping_SourceVariableTypeSet(vectorMapping,FIELD_U_VARIABLE_TYPE,err,error,*999)
              ENDIF
              CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE, &
                 & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)              
                 EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
                 CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,EQUATIONS_SET_FIELD_DATA,err,error,*999)
                 imy_matrix = EQUATIONS_SET_FIELD_DATA(1)
                 Ncompartments = EQUATIONS_SET_FIELD_DATA(2)    
                 CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_DELUDELN_VARIABLE_TYPE,vectorMapping,err,error,*999)
                 CALL EquationsMapping_DynamicMatricesSet(vectorMapping,.TRUE.,.TRUE.,err,error,*999)
                 CALL EquationsMapping_LinearMatricesNumberSet(vectorMapping,Ncompartments-1,err,error,*999)

                 ALLOCATE(VARIABLE_TYPES(2*Ncompartments))
                 ALLOCATE(VARIABLE_U_TYPES(Ncompartments-1))
                 DO num_var=1,Ncompartments
                   VARIABLE_TYPES(2*num_var-1)=FIELD_U_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(num_var-1))
                   VARIABLE_TYPES(2*num_var)=FIELD_DELUDELN_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(num_var-1))
                 ENDDO
                 num_var_count=0
                 DO num_var=1,Ncompartments
                   IF(num_var/=imy_matrix)THEN
                     num_var_count=num_var_count+1
                     VARIABLE_U_TYPES(num_var_count)=VARIABLE_TYPES(2*num_var-1)
                   ENDIF
                 ENDDO
                 CALL EquationsMapping_DynamicVariableTypeSet(vectorMapping,VARIABLE_TYPES(2*imy_matrix-1),err,error,*999)
                 CALL EquationsMapping_LinearMatricesVariableTypesSet(vectorMapping,VARIABLE_U_TYPES,err,error,*999)
                 CALL EquationsMapping_RHSVariableTypeSet(vectorMapping,VARIABLE_TYPES(2*imy_matrix),err,error,*999)
                 CALL EquationsMapping_SourceVariableTypeSet(vectorMapping,FIELD_U_VARIABLE_TYPE,err,error,*999)
              END SELECT
              CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
              !Create the equations matrices
              CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
              !Set up matrix storage and structure
              IF(EQUATIONS%lumpingType==EQUATIONS_LUMPED_MATRICES) THEN
                !Set up lumping
                CALL EquationsMatrices_DynamicLumpingTypeSet(vectorMatrices, &
                  & [EQUATIONS_MATRIX_UNLUMPED,EQUATIONS_MATRIX_LUMPED],err,error,*999)
                CALL EquationsMatrices_DynamicStorageTypeSet(vectorMatrices, &
                  & [DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE],err,error,*999)
                CALL EquationsMatrices_DynamicStructureTypeSet(vectorMatrices, &
                  [EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_DIAGONAL_STRUCTURE],err,error,*999)
              ELSE
                SELECT CASE(EQUATIONS%sparsityType)
                CASE(EQUATIONS_MATRICES_FULL_MATRICES) 
                  CALL EquationsMatrices_LinearStorageTypeSet(vectorMatrices, &
                    & [DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE,DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE],err,error,*999)
                CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                  CALL EquationsMatrices_DynamicStorageTypeSet(vectorMatrices, &
                    & [DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE], &
                    & err,error,*999)
                  CALL EquationsMatrices_DynamicStructureTypeSet(vectorMatrices, &
                    [EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_FEM_STRUCTURE],err,error,*999)    
                  IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
                   & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)THEN
                    ALLOCATE(COUPLING_MATRIX_STORAGE_TYPE(Ncompartments-1))
                    ALLOCATE(COUPLING_MATRIX_STRUCTURE_TYPE(Ncompartments-1))
                    DO num_var=1,Ncompartments-1
                     COUPLING_MATRIX_STORAGE_TYPE(num_var)=DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE
                     COUPLING_MATRIX_STRUCTURE_TYPE(num_var)=EQUATIONS_MATRIX_FEM_STRUCTURE
                    ENDDO                    
                    CALL EquationsMatrices_LinearStorageTypeSet(vectorMatrices,COUPLING_MATRIX_STORAGE_TYPE,err,error,*999)
                    CALL EquationsMatrices_LinearStructureTypeSet(vectorMatrices,COUPLING_MATRIX_STRUCTURE_TYPE,err,error,*999)
                  ENDIF                          
                CASE DEFAULT
                  localError="The equations matrices sparsity type of "// &
                    & TRIM(NumberToVString(EQUATIONS%sparsityType,"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              ENDIF
              CALL EquationsMatrices_VectorCreateFinish(vectorMatrices,err,error,*999)
            ELSE IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
              & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
              & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
              & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
              & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
              & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN
              !Finish the creation of the equations
              CALL EquationsSet_EquationsGet(EQUATIONS_SET,EQUATIONS,err,error,*999)
              CALL Equations_CreateFinish(EQUATIONS,err,error,*999)
              NULLIFY(vectorEquations)
              CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
              !Create the equations mapping.
              CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_DELUDELN_VARIABLE_TYPE,vectorMapping,err,error,*999)
              CALL EquationsMapping_LinearMatricesNumberSet(vectorMapping,1,err,error,*999)
              CALL EquationsMapping_LinearMatricesVariableTypesSet(vectorMapping,[FIELD_U_VARIABLE_TYPE], &
                & err,error,*999)
              CALL EquationsMapping_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE,err,error,*999)
              IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                   & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                   & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
                   & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN              
                CALL EquationsMapping_SourceVariableTypeSet(vectorMapping,FIELD_U_VARIABLE_TYPE,err,error,*999)
              ENDIF
              CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
              !Create the equations matrices
              CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
              SELECT CASE(EQUATIONS%sparsityType)
              CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                CALL EquationsMatrices_LinearStorageTypeSet(vectorMatrices,[MATRIX_BLOCK_STORAGE_TYPE],err,error,*999)
              CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                CALL EquationsMatrices_LinearStorageTypeSet(vectorMatrices,[MATRIX_COMPRESSED_ROW_STORAGE_TYPE],err,error,*999)
                CALL EquationsMatrices_LinearStructureTypeSet(vectorMatrices,[EQUATIONS_MATRIX_FEM_STRUCTURE],err,error,*999)
              CASE DEFAULT
                localError="The equations matrices sparsity type of "// &
                  & TRIM(NumberToVString(EQUATIONS%sparsityType,"*",err,error))//" is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
              CALL EquationsMatrices_VectorCreateFinish(vectorMatrices,err,error,*999)
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
                localError="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                & " is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a linear advection-diffusion equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
            & " is invalid for a linear advection-diffusion equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        localError="The equations set subtype of "//TRIM(NumberToVString(EQUATIONS_SET_SUBTYPE,"*",err,error))// &
          & " is not a linear advection-diffusion equation subtype."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
       
    EXITS("AdvectionDiffusion_EquationsSetLinearSetup")
    RETURN
999 ERRORS("AdvectionDiffusion_EquationsSetLinearSetup",err,error)
    EXITS("AdvectionDiffusion_EquationsSetLinearSetup")
    RETURN 1
    
  END SUBROUTINE AdvectionDiffusion_EquationsSetLinearSetup

  !
  !===============================================================================================================================
  !
!   SUBROUTINE ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*)
!     !Argument variables
!     TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup
!     TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
!     INTEGER(INTG), INTENT(OUT) :: err !<The error code
!     TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
! 
!     CALL FlagError("Not implemented.",err,error,*999)
!     EXITS("ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP")
!     RETURN
! 999 ERRORSEXITS("ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP",err,error)
!     RETURN 1
! 
!   END SUBROUTINE ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP

  !
  !================================================================================================================================
  !
 
  !>Sets up the diffusion problem.
  SUBROUTINE ADVECTION_DIFFUSION_EQUATION_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem set to setup a diffusion equation on.
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("ADVECTION_DIFFUSION_EQUATION_PROBLEM_SETUP",err,error,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(.NOT.ALLOCATED(PROBLEM%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(PROBLEM%SPECIFICATION,1)<3) THEN
        CALL FlagError("Problem specification does not have a subtype set.",err,error,*999)
      END IF
      SELECT CASE(PROBLEM%SPECIFICATION(3))
      CASE(PROBLEM_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE)
        CALL AdvectionDiffusion_ProblemLinearSetup(PROBLEM,PROBLEM_SETUP,err,error,*999)
      CASE(PROBLEM_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE)
        CALL AdvectionDiffusion_ProblemLinearSetup(PROBLEM,PROBLEM_SETUP,err,error,*999)
      CASE(PROBLEM_NONLINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE)
        CALL FlagError("Not implemented.",err,error,*999)
!        CALL ADVECTION_DIFFUSION_EQUATION_PROBLEM_NONLINEAR_SETUP(PROBLEM,PROBLEM_SETUP,err,error,*999)
      CASE(PROBLEM_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
        CALL AdvectionDiffusion_ProblemLinearSetup(PROBLEM,PROBLEM_SETUP,err,error,*999)
      CASE(PROBLEM_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
        CALL AdvectionDiffusion_ProblemLinearSetup(PROBLEM,PROBLEM_SETUP,err,error,*999)
      CASE(PROBLEM_NONLINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
        CALL FlagError("Not implemented.",err,error,*999)
!        CALL ADVECTION_DIFFUSION_EQUATION_PROBLEM_NONLINEAR_SETUP(PROBLEM,PROBLEM_SETUP,err,error,*999)
      CASE(PROBLEM_NO_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE)
        CALL AdvectionDiffusion_ProblemLinearSetup(PROBLEM,PROBLEM_SETUP,err,error,*999)
      CASE(PROBLEM_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE)
        CALL AdvectionDiffusion_ProblemLinearSetup(PROBLEM,PROBLEM_SETUP,err,error,*999)
      CASE(PROBLEM_NONLINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE)
        CALL FlagError("Not implemented.",err,error,*999)
!        CALL ADVECTION_DIFFUSION_EQUATION_PROBLEM_NONLINEAR_SETUP(PROBLEM,PROBLEM_SETUP,err,error,*999)
      CASE DEFAULT
        localError="Problem subtype "//TRIM(NumberToVString(PROBLEM%SPECIFICATION(3),"*",err,error))// &
          & " is not valid for an advection-diffusion equation type of a classical field problem class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    ENDIF
       
    EXITS("ADVECTION_DIFFUSION_EQUATION_PROBLEM_SETUP")
    RETURN
999 ERRORSEXITS("ADVECTION_DIFFUSION_EQUATION_PROBLEM_SETUP",err,error)
    RETURN 1
  END SUBROUTINE ADVECTION_DIFFUSION_EQUATION_PROBLEM_SETUP
  
  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices and RHS for a diffusion equation finite element equations set.
  SUBROUTINE AdvectionDiffusion_FiniteElementCalculate(EQUATIONS_SET,ELEMENT_NUMBER,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) mh,mhs,ms,ng,nh,nhs,ni,nj,ns,FIELD_VAR_TYPE,my_compartment,Ncompartments,imatrix,num_var_count
    INTEGER(INTG) :: MESH_COMPONENT_1, MESH_COMPONENT_2
    REAL(DP) :: C_PARAM,K_PARAM,RWG,SUM,PGMJ(3),PGNJ(3),ADVEC_VEL,A_PARAM,COUPLING_PARAM,PGM,PGN
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS,GEOMETRIC_BASIS
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS_1, DEPENDENT_BASIS_2
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatrixType), POINTER :: dampingMatrix,stiffnessMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,GEOMETRIC_FIELD,MATERIALS_FIELD,SOURCE_FIELD,INDEPENDENT_FIELD,EQUATIONS_SET_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE,GEOMETRIC_VARIABLE
    TYPE(FIELD_VARIABLE_PTR_TYPE) :: FIELD_VARIABLES(99)
    TYPE(EquationsMatrixPtrType) :: COUPLING_MATRICES(99) 
    INTEGER(INTG) :: FIELD_VAR_TYPES(99)
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME_1, QUADRATURE_SCHEME_2
    TYPE(VARYING_STRING) :: localError
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: DIFFUSION_DEPENDENT_PREVIOUS_INTERPOLATION_PARAMETERS
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: DIFFUSION_DEPENDENT_PREVIOUS_INTERPOLATED_POINT
    INTEGER(INTG), POINTER :: EQUATIONS_SET_FIELD_DATA(:)
    LOGICAL :: updateDampingMatrix,updateStiffnessMatrix,updateRHSVector,updateSourceVector
    INTEGER(INTG) :: EQUATIONS_SET_SUBTYPE

    updateDampingMatrix = .FALSE.
    updateStiffnessMatrix = .FALSE.
    updateRHSVector = .FALSE.
    updateSourceVector = .FALSE.

    ENTERS("AdvectionDiffusion_FiniteElementCalculate",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      
      IF(ASSOCIATED(EQUATIONS)) THEN
        IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
          CALL FlagError("Equations set specification is not allocated.",err,error,*999)
        ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<3) THEN
          CALL FlagError("Equations set specification does not have a subtype set.",err,error,*999)
        END IF
        NULLIFY(vectorEquations)
        CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
        EQUATIONS_SET_SUBTYPE=EQUATIONS_SET%SPECIFICATION(3)
        SELECT CASE(EQUATIONS_SET_SUBTYPE)
        CASE(EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE,EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE,EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE, &
        & EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE, &
        & EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE, &
        & EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
        & EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
        & EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
        & EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
        & EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
        & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
        & EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
        & EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
        & EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
        & EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE, &
        & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)
          !Store all these in equations matrices/somewhere else?????
          DEPENDENT_FIELD=>equations%interpolation%dependentField
          GEOMETRIC_FIELD=>equations%interpolation%geometricField
          MATERIALS_FIELD=>equations%interpolation%materialsField
          INDEPENDENT_FIELD=>equations%interpolation%independentField !Stores the advective velocity field
          IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
             SOURCE_FIELD=>equations%interpolation%sourceField
          ENDIF
          vectorMatrices=>vectorEquations%vectorMatrices
          IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. & 
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
             dynamicMatrices=>vectorMatrices%dynamicMatrices
             stiffnessMatrix=>dynamicMatrices%matrices(1)%ptr
             dampingMatrix=>dynamicMatrices%matrices(2)%ptr
          ELSEIF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN
             linearMatrices=>vectorMatrices%linearMatrices
             stiffnessMatrix=>linearMatrices%matrices(1)%ptr
          ELSE

          ENDIF
          rhsVector=>vectorMatrices%rhsVector
          IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. & 
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN            
             sourceVector=>vectorMatrices%sourceVector
          ENDIF
          IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. & 
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
             IF(ASSOCIATED(dampingMatrix)) updateDampingMatrix=dampingMatrix%updateMatrix
          ENDIF
          IF(ASSOCIATED(stiffnessMatrix)) updateStiffnessMatrix=stiffnessMatrix%updateMatrix
          IF(ASSOCIATED(rhsVector)) updateRHSVector=rhsVector%updateVector
          IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. & 
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN             
             IF(ASSOCIATED(sourceVector)) updateSourceVector=sourceVector%updateVector
          ENDIF
          vectorMapping=>vectorEquations%vectorMapping
          IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
           & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
           EQUATIONS_SET_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
           CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET_FIELD,FIELD_U_VARIABLE_TYPE, &
             & FIELD_VALUES_SET_TYPE,EQUATIONS_SET_FIELD_DATA,err,error,*999)
            my_compartment = EQUATIONS_SET_FIELD_DATA(1)
            Ncompartments  = EQUATIONS_SET_FIELD_DATA(2)
            linearMatrices=>vectorMatrices%linearMatrices
            linearMapping=>vectorMapping%linearMapping
            num_var_count=0
            DO imatrix = 1,Ncompartments
             IF(imatrix/=my_compartment)THEN
              num_var_count=num_var_count+1
              COUPLING_MATRICES(num_var_count)%ptr=>linearMatrices%matrices(num_var_count)%ptr
              FIELD_VARIABLES(num_var_count)%ptr=>linearMapping%equationsMatrixToVarMaps(num_var_count)%VARIABLE
              FIELD_VAR_TYPES(num_var_count)=FIELD_VARIABLES(num_var_count)%ptr%VARIABLE_TYPE
              COUPLING_MATRICES(num_var_count)%ptr%elementMatrix%MATRIX=0.0_DP
             ENDIF
            END DO
          ENDIF
          IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &  
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
             dynamicMapping=>vectorMapping%dynamicMapping
             FIELD_VARIABLE=>dynamicMapping%equationsMatrixToVarMaps(1)%VARIABLE
          ELSEIF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN
             linearMapping=>vectorMapping%linearMapping
             FIELD_VARIABLE=>linearMapping%equationsMatrixToVarMaps(1)%VARIABLE
          ELSE

          ENDIF
          FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
          GEOMETRIC_VARIABLE=>GEOMETRIC_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%ptr
          DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%ptr% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          GEOMETRIC_BASIS=>GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%ptr% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
            & geometricInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
            & materialsInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
          IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE) THEN
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
               & materialsInterpParameters(FIELD_V_VARIABLE_TYPE)%ptr,err,error,*999)
          ENDIF

          IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
              & sourceInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
          ENDIF
          !the following line has been changed to use fieldinputdata1settype
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
            & independentInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
          IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE) THEN  
           DIFFUSION_DEPENDENT_PREVIOUS_INTERPOLATION_PARAMETERS=> &
            & equations%interpolation%dependentInterpParameters(FIELD_V_VARIABLE_TYPE)%ptr
           CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
            & DIFFUSION_DEPENDENT_PREVIOUS_INTERPOLATION_PARAMETERS,err,error,*999)
           DIFFUSION_DEPENDENT_PREVIOUS_INTERPOLATED_POINT=> &
            & equations%interpolation%dependentInterpPoint(FIELD_V_VARIABLE_TYPE)%ptr
          ENDIF


          !Select whether using standard Galerkin scheme, or the stabilised streamwise-upwinding Petrov-Galerkin scheme
          SELECT CASE(EQUATIONS_SET_SUBTYPE)
          CASE(EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE,EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
               & EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
               & EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE,& 
               & EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
               & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
               & EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE, &
               & EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE, &
               & EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE , & 
               & EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE, &
               & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE)
          !Loop over gauss points
          DO ng=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
              & geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,equations%interpolation% &
              & geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
              & materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            !Interpolate to get the advective velocity
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
              & independentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)            
            IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE) THEN
              CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
                & sourceInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)   
            ENDIF
            IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE) THEN
              CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
                 & materialsInterpPoint(FIELD_V_VARIABLE_TYPE)%ptr,err,error,*999)
            ENDIF
            !Calculate RWG.
!!TODO: Think about symmetric problems. 
            RWG=equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%JACOBIAN* &
              & QUADRATURE_SCHEME%GAUSS_WEIGHTS(ng)
            !Loop over field components

            mhs=0          
            DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
              !Loop over element rows
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1
                nhs=0
                IF(updateStiffnessMatrix .OR. updateDampingMatrix) THEN
                  !Loop over element columns
                  DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1
                      IF(updateStiffnessMatrix) THEN
                        SUM=0.0_DP
                        DO nj=1,GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS
                          PGMJ(nj)=0.0_DP
                          PGNJ(nj)=0.0_DP
                          DO ni=1,DEPENDENT_BASIS%NUMBER_OF_XI                          
                            PGMJ(nj)=PGMJ(nj)+ &
                              & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)* &
                              & equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%DXI_DX(ni,nj)
                            PGNJ(nj)=PGNJ(nj)+ &
                              & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)* &
                              & equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%DXI_DX(ni,nj)
                          ENDDO !ni
                          K_PARAM=equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr% &
                            & VALUES(nj,NO_PART_DERIV)
                          SUM=SUM+K_PARAM*PGMJ(nj)*PGNJ(nj)
                          !Advection term is constructed here and then added to SUM for updating the stiffness matrix outside of this loop
                          ADVEC_VEL=equations%interpolation%independentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr% &
                            & VALUES(nj,NO_PART_DERIV)
                          SUM=SUM+ADVEC_VEL*QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)*PGNJ(nj)   
                        ENDDO !nj
                        IF (EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                            & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                            & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                            & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                            & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                            & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE) THEN
                          stiffnessMatrix%elementMatrix%matrix(mhs,nhs)=stiffnessMatrix%elementMatrix%matrix(mhs,nhs)+SUM*RWG
                        ELSEIF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                            & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                            & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. & 
                            & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE) THEN
                          A_PARAM=equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr% &
                            & VALUES(GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS,NO_PART_DERIV)
                          stiffnessMatrix%elementMatrix%matrix(mhs,nhs)=stiffnessMatrix%elementMatrix%matrix(mhs,nhs)+SUM*RWG- &
                            & A_PARAM*QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)* &
                            & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)*RWG
                            ! A_PARAM is the material parameter that multiplies the linear source u
                        ELSEIF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE) THEN
                        !for multi-compartment model must include additional terms into the
                          COUPLING_PARAM=equations%interpolation%materialsInterpPoint(FIELD_V_VARIABLE_TYPE)%ptr% &
                            & VALUES(my_compartment,NO_PART_DERIV)
                          stiffnessMatrix%elementMatrix%matrix(mhs,nhs)=stiffnessMatrix%elementMatrix%matrix(mhs,nhs)+ & 
                            & SUM*RWG + QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)* &
                          & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)*RWG*COUPLING_PARAM
                        ENDIF
                      ENDIF
                    IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE .OR. &
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE) THEN
                      IF(updateDampingMatrix) THEN
                        dampingMatrix%elementMatrix%matrix(mhs,nhs)=dampingMatrix%elementMatrix%matrix(mhs,nhs)+ &
                          & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)* &
                          & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)*RWG
                      ENDIF
                    ENDIF
                    ENDDO !ns
                  ENDDO !nh
                ENDIF
               ENDDO !ms
              ENDDO !mh
              IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE) THEN 
                IF(updateSourceVector) THEN
                    C_PARAM=equations%interpolation%sourceInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(1, NO_PART_DERIV)
                    mhs=0
                    DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    !DO mh=1,DEPENDENT_VARIABLE%NUMBER_OF_COMPONENTS
                     !Loop over element rows
                      DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                       mhs=mhs+1
                       sourceVector%elementVector%vector(mhs)=sourceVector%elementVector%vector(mhs)+ &
                          & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)*C_PARAM*RWG
                      ENDDO !ms
                    ENDDO !mh
                ENDIF
              ELSEIF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE) THEN
                IF(updateSourceVector) THEN
                    CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng, &
                          & DIFFUSION_DEPENDENT_PREVIOUS_INTERPOLATED_POINT,err,error,*999)
                    C_PARAM=DIFFUSION_DEPENDENT_PREVIOUS_INTERPOLATED_POINT%VALUES(1,NO_PART_DERIV)
                    mhs=0
                    DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    !DO mh=1,DEPENDENT_VARIABLE%NUMBER_OF_COMPONENTS
                     !Loop over element rows
                      DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                       mhs=mhs+1
                       sourceVector%elementVector%vector(mhs)=sourceVector%elementVector%vector(mhs)+ &
                          & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)*C_PARAM*RWG
                      ENDDO !ms
                    ENDDO !mh
                ENDIF
              ENDIF
            IF(updateRHSVector) rhsVector%elementVector%vector(mhs)=0.0_DP
  
            IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE) THEN
            !Calculate the coupling matrices

              !Loop over element rows
              mhs=0
              DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS !field_variable is the variable associated with the equations set under consideration

                MESH_COMPONENT_1 = FIELD_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
                DEPENDENT_BASIS_1 => DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT_1)%ptr% &
                  & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                QUADRATURE_SCHEME_1 => DEPENDENT_BASIS_1%QUADRATURE% &
                  & QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
                RWG = equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%JACOBIAN * &
                  & QUADRATURE_SCHEME_1%GAUSS_WEIGHTS(ng)

                DO ms=1,DEPENDENT_BASIS_1%NUMBER_OF_ELEMENT_PARAMETERS
                  mhs=mhs+1

                  num_var_count=0
                  DO imatrix = 1,Ncompartments
                  IF(imatrix/=my_compartment)THEN
                    num_var_count=num_var_count+1

!need to test for the case where imatrix==mycompartment
!the coupling terms then needs to be added into the stiffness matrix
                    IF(COUPLING_MATRICES(num_var_count)%ptr%updateMatrix) THEN

!                       !Loop over element columns
                      nhs=0
! !                       DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                      DO nh=1,FIELD_VARIABLES(num_var_count)%ptr%NUMBER_OF_COMPONENTS

                        MESH_COMPONENT_2 = FIELD_VARIABLE%COMPONENTS(nh)%MESH_COMPONENT_NUMBER
                        DEPENDENT_BASIS_2 => DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT_2)%ptr% &
                          & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                        !--- We cannot use two different quadrature schemes here !!!
                        QUADRATURE_SCHEME_2 => DEPENDENT_BASIS_2%QUADRATURE% &
                         & QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
                        !RWG = equations%interpolation%geometricInterpPointMetrics%JACOBIAN * &
                        !  & QUADRATURE_SCHEME_2%GAUSS_WEIGHTS(ng)

                        DO ns=1,DEPENDENT_BASIS_2%NUMBER_OF_ELEMENT_PARAMETERS
                          nhs=nhs+1

!                           !-------------------------------------------------------------------------------------------------------------
!                           !concentration test function, concentration trial function
!                           !For now, this is only a dummy implementation - this still has to be properly set up.
!                           IF(mh==nh.AND.nh<NUMBER_OF_VEL_PRESS_COMPONENTS) THEN ! don't need this for diffusion equation

!                             SUM = 0.0_DP

                            PGM=QUADRATURE_SCHEME_1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                            PGN=QUADRATURE_SCHEME_2%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)

                            !Get the coupling coefficients 
                              COUPLING_PARAM=equations%interpolation%materialsInterpPoint(FIELD_V_VARIABLE_TYPE)%ptr% &
                                & VALUES(imatrix,NO_PART_DERIV)

!                              SUM = SUM + COUPLING_PARAM * PGM * PGN
 
                             COUPLING_MATRICES(num_var_count)%ptr%elementMatrix%matrix(mhs,nhs) = &
                               & COUPLING_MATRICES(num_var_count)%ptr%elementMatrix%matrix(mhs,nhs) + & 
                               & COUPLING_PARAM * PGM * PGN * RWG
!                           ENDIF
 
                        ENDDO !ns
                      ENDDO !nh
                    ENDIF
                   ENDIF
                  ENDDO !imatrix
                ENDDO !ms
              ENDDO !mh

            ENDIF

          ENDDO !ng
          
          !Scale factor adjustment
          IF(DEPENDENT_FIELD%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
            CALL Field_InterpolationParametersScaleFactorsElementGet(ELEMENT_NUMBER,equations%interpolation% &
              & dependentInterpParameters(FIELD_VAR_TYPE)%ptr,err,error,*999)
            mhs=0          
            DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
              !Loop over element rows
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1                    
                nhs=0
                IF(updateStiffnessMatrix .OR. updateDampingMatrix) THEN
                  !Loop over element columns
                  DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1
                      IF(updateStiffnessMatrix) THEN
                        stiffnessMatrix%elementMatrix%matrix(mhs,nhs)=stiffnessMatrix%elementMatrix%matrix(mhs,nhs)* &
                          & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%SCALE_FACTORS(ms,mh)* &
                          & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%SCALE_FACTORS(ns,nh)
                      ENDIF
                    IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE .OR. & 
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE) THEN
                      IF(updateDampingMatrix) THEN
                        dampingMatrix%elementMatrix%matrix(mhs,nhs)=dampingMatrix%elementMatrix%matrix(mhs,nhs)* &
                          & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%SCALE_FACTORS(ms,mh)* &
                          & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%SCALE_FACTORS(ns,nh)
                      ENDIF
                    ENDIF
                    ENDDO !ns
                  ENDDO !nh
                ENDIF
                IF(updateRHSVector) rhsVector%elementVector%vector(mhs)=rhsVector%elementVector%vector(mhs)* &
                  & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%SCALE_FACTORS(ms,mh)
              IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE) THEN
                IF(updateSourceVector) sourceVector%elementVector%vector(mhs)=sourceVector%elementVector%vector(mhs)* &
                  & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%SCALE_FACTORS(ms,mh)
              ENDIF
              ENDDO !ms
            ENDDO !mh
          ENDIF
          CASE(EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)
            CALL FlagError("Not implemented.",err,error,*999) 
          !Loop over gauss points
          DO ng=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
              & geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,equations%interpolation% &
              & geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
              & materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            !Interpolate to get the advective velocity
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
              & independentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)            
            !Calculate RWG.
!!TODO: Think about symmetric problems. 
            RWG=equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%JACOBIAN* &
              & QUADRATURE_SCHEME%GAUSS_WEIGHTS(ng)
            !Loop over field components
            mhs=0          
            DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
              !Loop over element rows
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1
                nhs=0
                IF(updateStiffnessMatrix .OR. updateDampingMatrix) THEN
                  !Loop over element columns
                  DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1
                      IF(updateStiffnessMatrix) THEN
                        SUM=0.0_DP
                        DO nj=1,GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS
                          PGMJ(nj)=0.0_DP
                          PGNJ(nj)=0.0_DP
                          DO ni=1,DEPENDENT_BASIS%NUMBER_OF_XI                          
                            PGMJ(nj)=PGMJ(nj)+ &
                              & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)* &
                              & equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%DXI_DX(ni,nj)
                            PGNJ(nj)=PGNJ(nj)+ &
                              & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)* &
                              & equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%DXI_DX(ni,nj)
                          ENDDO !ni
                          K_PARAM=equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr% &
                            & VALUES(nj,NO_PART_DERIV)
                          SUM=SUM+K_PARAM*PGMJ(nj)*PGNJ(nj)
                          !Advection term is constructed here and then added to SUM for updating the stiffness matrix outside of this loop
                          ADVEC_VEL=equations%interpolation%independentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr% &
                            & VALUES(nj,NO_PART_DERIV)
                          SUM=SUM+ADVEC_VEL*QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)*PGNJ(nj)   
                        ENDDO !nj
                        IF (EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                            & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                            & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                            & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                            & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
                            & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN
                          stiffnessMatrix%elementMatrix%matrix(mhs,nhs)=stiffnessMatrix%elementMatrix%matrix(mhs,nhs)+SUM*RWG
                        ELSEIF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                            & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                            & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN
                          A_PARAM=equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr% &
                            & VALUES(GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS,NO_PART_DERIV)
                          stiffnessMatrix%elementMatrix%matrix(mhs,nhs)=stiffnessMatrix%elementMatrix%matrix(mhs,nhs)+SUM*RWG- &
                            & A_PARAM*QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)* &
                            & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)*RWG
                            ! A_PARAM is the material parameter that multiplies the linear source u
                        ENDIF
                      ENDIF
                    IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE) THEN
                      IF(updateDampingMatrix) THEN
                        dampingMatrix%elementMatrix%matrix(mhs,nhs)=dampingMatrix%elementMatrix%matrix(mhs,nhs)+ &
                          & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)* &
                          & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)*RWG
                      ENDIF
                    ENDIF
                    ENDDO !ns
                  ENDDO !nh
                ENDIF
               ENDDO !ms
              ENDDO !mh
              IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN 
                IF(updateSourceVector) THEN
                    C_PARAM=equations%interpolation%sourceInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(1, NO_PART_DERIV)
                    mhs=0
                    DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    !DO mh=1,DEPENDENT_VARIABLE%NUMBER_OF_COMPONENTS
                     !Loop over element rows
                      DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                       mhs=mhs+1
                       sourceVector%elementVector%vector(mhs)=sourceVector%elementVector%vector(mhs)+ &
                          & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)*C_PARAM*RWG
                      ENDDO !ms
                    ENDDO !mh
                ENDIF
              ENDIF
            IF(updateRHSVector) rhsVector%elementVector%vector(mhs)=0.0_DP
          ENDDO !ng
          
          !Scale factor adjustment
          IF(DEPENDENT_FIELD%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
            CALL Field_InterpolationParametersScaleFactorsElementGet(ELEMENT_NUMBER,equations%interpolation% &
              & dependentInterpParameters(FIELD_VAR_TYPE)%ptr,err,error,*999)
            mhs=0          
            DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
              !Loop over element rows
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1                    
                nhs=0
                IF(updateStiffnessMatrix .OR. updateDampingMatrix) THEN
                  !Loop over element columns
                  DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1
                      IF(updateStiffnessMatrix) THEN
                        stiffnessMatrix%elementMatrix%matrix(mhs,nhs)=stiffnessMatrix%elementMatrix%matrix(mhs,nhs)* &
                          & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%SCALE_FACTORS(ms,mh)* &
                          & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%SCALE_FACTORS(ns,nh)
                      ENDIF
                    IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                       & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE) THEN
                      IF(updateDampingMatrix) THEN
                        dampingMatrix%elementMatrix%matrix(mhs,nhs)=dampingMatrix%elementMatrix%matrix(mhs,nhs)* &
                          & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%SCALE_FACTORS(ms,mh)* &
                          & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%SCALE_FACTORS(ns,nh)
                      ENDIF
                    ENDIF
                    ENDDO !ns
                  ENDDO !nh
                ENDIF
                IF(updateRHSVector) rhsVector%elementVector%vector(mhs)=rhsVector%elementVector%vector(mhs)* &
                  & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%SCALE_FACTORS(ms,mh)
              IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
                 & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN
                IF(updateSourceVector) sourceVector%elementVector%vector(mhs)=sourceVector%elementVector%vector(mhs)* &
                  & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%SCALE_FACTORS(ms,mh)
              ENDIF
              ENDDO !ms
            ENDDO !mh
          ENDIF
          CASE DEFAULT
            localError="Equations set subtype "//TRIM(NumberToVString(EQUATIONS_SET_SUBTYPE,"*",err,error))// &
              & " is not valid for an advection-diffusion equation type of a classical field equations set class."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_QUADRATIC_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
             & EQUATIONS_SET_EXPONENTIAL_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
             & EQUATIONS_SET_QUADRATIC_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
             & EQUATIONS_SET_EXP_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
             & EQUATIONS_SET_QUAD_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
             & EQUATIONS_SET_EXP_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
             & EQUATIONS_SET_QUAD_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
             & EQUATIONS_SET_EXP_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE)
          CALL FlagError("Can not calculate finite element stiffness matrices for a nonlinear source.",err,error,*999)
        CASE DEFAULT
          localError="Equations set subtype "//TRIM(NumberToVString(EQUATIONS_SET_SUBTYPE,"*",err,error))// &
            & " is not valid for an advection-diffusion equation type of a classical field equations set class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        CALL FlagError("Equations set equations is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
    
    EXITS("AdvectionDiffusion_FiniteElementCalculate")
    RETURN
999 ERRORS("AdvectionDiffusion_FiniteElementCalculate",err,error)
    EXITS("AdvectionDiffusion_FiniteElementCalculate")
    RETURN 1
    
  END SUBROUTINE AdvectionDiffusion_FiniteElementCalculate

  !
  !================================================================================================================================
  !

  !INSERT CODE HERE TO DEAL WITH THE NON-LINEAR SOURCE TERMS - DECIDE HOW TO SOLVE THEM, USE JACOBIAN AS FOR NON-LINEAR POISSON EXAMPLE?

  !
  !================================================================================================================================
  !

  !>Sets the problem specification for an advection diffusion problem type.
  SUBROUTINE AdvectionDiffusion_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: problem !<A pointer to the problem to set the problem specification for
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: problemSubtype

    ENTERS("AdvectionDiffusion_ProblemSpecificationSet",err,error,*999)

    IF(ASSOCIATED(problem)) THEN
      IF(SIZE(problemSpecification,1)==3) THEN
        problemSubtype=problemSpecification(3)
        SELECT CASE(problemSubtype)
        CASE(PROBLEM_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
          & PROBLEM_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
          & PROBLEM_NONLINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
          & PROBLEM_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
          & PROBLEM_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
          & PROBLEM_NONLINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
          & PROBLEM_NO_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE, &
          & PROBLEM_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE, &
          & PROBLEM_NONLINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE)
          !ok
        CASE DEFAULT
          localError="The third problem specification of "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
            & " is not valid for an advection-diffusion type of a classical field problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        IF(ALLOCATED(problem%specification)) THEN
          CALL FlagError("Problem specification is already allocated.",err,error,*999)
        ELSE
          ALLOCATE(problem%specification(3),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
        END IF
        problem%specification(1:3)=[PROBLEM_CLASSICAL_FIELD_CLASS,PROBLEM_ADVECTION_DIFFUSION_EQUATION_TYPE,problemSubtype]
      ELSE
        CALL FlagError("Advection-diffusion problem specification must have three entries.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    END IF

    EXITS("AdvectionDiffusion_ProblemSpecificationSet")
    RETURN
999 ERRORS("AdvectionDiffusion_ProblemSpecificationSet",err,error)
    EXITS("AdvectionDiffusion_ProblemSpecificationSet")
    RETURN 1
    
  END SUBROUTINE AdvectionDiffusion_ProblemSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the diffusion equations.
  SUBROUTINE AdvectionDiffusion_ProblemLinearSetup(PROBLEM,PROBLEM_SETUP,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to setup
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP,CONTROL_LOOP_ROOT
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: PROBLEM_SUBTYPE
    
    ENTERS("AdvectionDiffusion_ProblemLinearSetup",err,error,*999)

    NULLIFY(CONTROL_LOOP)
    NULLIFY(SOLVER)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVERS)
    IF(ASSOCIATED(PROBLEM)) THEN
      IF(.NOT.ALLOCATED(PROBLEM%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(PROBLEM%SPECIFICATION,1)<3) THEN
        CALL FlagError("Problem specification does not have a subtype set.",err,error,*999)
      END IF
      PROBLEM_SUBTYPE=PROBLEM%SPECIFICATION(3)
      IF(PROBLEM_SUBTYPE==PROBLEM_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
         & PROBLEM_SUBTYPE==PROBLEM_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
         & PROBLEM_SUBTYPE==PROBLEM_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
         & PROBLEM_SUBTYPE==PROBLEM_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE) THEN
        SELECT CASE(PROBLEM_SETUP%SETUP_TYPE)
        CASE(PROBLEM_SETUP_INITIAL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Do nothing????
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Do nothing????
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a linear advection-diffusion equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Set up a time control loop
            CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,err,error,*999)
            CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,PROBLEM_CONTROL_TIME_LOOP_TYPE,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Finish the control loops
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,err,error,*999)            
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a linear advection-diffusion equation."
            CALL FlagError(localError,err,error,*999)
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
            !Set the solver to be a first order dynamic solver 
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,err,error,*999)
            CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
            !Set solver defaults
            CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
            CALL SOLVER_DYNAMIC_SCHEME_SET(SOLVER,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,err,error,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the solvers
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(SOLVERS,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a linear advection-diffusion equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            !Get the solver
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
            !Create the solver equations
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,err,error,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,err,error,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            !Get the solver equations
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,err,error,*999)
            !Finish the solver equations creation
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,err,error,*999)             
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a linear advection-diffusion equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
            & " is invalid for a linear advection-diffusion equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSEIF(PROBLEM_SUBTYPE==PROBLEM_NO_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
         & PROBLEM_SUBTYPE==PROBLEM_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE) THEN
        SELECT CASE(PROBLEM_SETUP%SETUP_TYPE)
        CASE(PROBLEM_SETUP_INITIAL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Do nothing????
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Do nothing????
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a linear static advection-diffusion equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Set up a simple control loop
            CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Finish the control loops
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,err,error,*999)            
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a linear static advection-diffusion equation."
            CALL FlagError(localError,err,error,*999)
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
            !Set the solver to be a linear solver
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
             !Start the linear solver creation
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_LINEAR_TYPE,err,error,*999)
            !Set solver defaults
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the solvers
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(SOLVERS,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a linear static advection-diffusion equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            !Get the solver
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
            !Create the solver equations
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,err,error,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,err,error,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,err,error,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            !Get the solver equations
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,err,error,*999)
            !Finish the solver equations creation
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,err,error,*999)             
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a linear static advection-diffusion equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
            & " is invalid for a linear static advection-diffusion equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        localError="The problem subtype of "//TRIM(NumberToVString(PROBLEM_SUBTYPE,"*",err,error))// &
          & " does not equal a linear advection-diffusion equation subtype."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    ENDIF
       
    EXITS("AdvectionDiffusion_ProblemLinearSetup")
    RETURN
999 ERRORS("AdvectionDiffusion_ProblemLinearSetup",err,error)
    EXITS("AdvectionDiffusion_ProblemLinearSetup")
    RETURN 1
    
  END SUBROUTINE AdvectionDiffusion_ProblemLinearSetup
  
  !
  !================================================================================================================================
  !
    !>Sets up the diffusion equations.
!   SUBROUTINE ADVECTION_DIFFUSION_EQUATION_PROBLEM_NONLINEAR_SETUP(PROBLEM,PROBLEM_SETUP,err,error,*)
! 
!     !Argument variables
!     TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to setup
!     TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
!     INTEGER(INTG), INTENT(OUT) :: err !<The error code
!     TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
! 
!     CALL FlagError("Not implemented.",err,error,*999) 
! 
!     EXITS("ADVECTION_DIFFUSION_EQUATION_PROBLEM_NONLINEAR_SETUP")
!     RETURN
! 999 ERRORSEXITS("ADVECTION_DIFFUSION_EQUATION_PROBLEM_NONLINEAR_SETUP",err,error)
!     RETURN 1
! 
!   END SUBROUTINE ADVECTION_DIFFUSION_EQUATION_PROBLEM_NONLINEAR_SETUP

  !
  !================================================================================================================================
  !

  !>Performs pre-solve operations for advection-diffusion problems.
  SUBROUTINE AdvectionDiffusion_PreSolve(solver,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: problemSubType
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop
    TYPE(PROBLEM_TYPE), POINTER :: problem
    TYPE(SOLVER_TYPE), POINTER :: SOLVER2 !<A pointer to the solvers
    TYPE(VARYING_STRING) :: localError

    ENTERS("AdvectionDiffusion_PreSolve",err,error,*999)
    
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    IF(.NOT.ALLOCATED(problem%specification)) CALL FlagError("Problem specification is not allocated.",err,error,*999)
    IF(SIZE(problem%specification,1)<3) &
      & CALL FlagError("Problem specification must have at least three entries for an advection-diffusion problem.",err,error,*999)
    
    problemSubType=problem%specification(3)
    IF(problemSubType==PROBLEM_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
      & problemSubType==PROBLEM_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE) THEN
       CALL WriteString(GENERAL_OUTPUT_TYPE,"Read in vector data... ",err,error,*999)
       !Update independent data fields
       CALL AdvectionDiffusion_PreSolveUpdateInputData(controlLoop,solver,err,error,*999)
       !CALL ADVECTION_DIFFUSION_PRE_SOLVE_UPDATE_BC(controlLoop,solver,err,error,*999)
     ELSE IF(problemSubType==PROBLEM_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
       & problemSubType==PROBLEM_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE) THEN
       CALL WriteString(GENERAL_OUTPUT_TYPE,"Read in vector data... ",err,error,*999)
       !Update independent data fields
       CALL AdvectionDiffusion_PreSolveUpdateInputData(controlLoop,solver,err,error,*999)
       !CALL ADVECTION_DIFFUSION_PRE_SOLVE_UPDATE_BC(controlLoop,solver,err,error,*999)
       CALL AdvectionDiffusion_PreSolveALEUpdateMesh(controlLoop,solver,err,error,*999)
     ELSE IF(problemSubType==PROBLEM_NONLINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
       & problemSubType==PROBLEM_NONLINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE) THEN         
       CALL FlagError("Not implemented.",err,error,*999)
     ELSE IF(problemSubType==PROBLEM_NO_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
       & problemSubType==PROBLEM_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE) THEN
       !do nothing
     ELSE
       localError="The third problem specification of "// &
         & TRIM(NumberToVString(problemSubType,"*",err,error))// &
         & " is not valid for a advection-diffusion type of a classical field problem."
       CALL FlagError(localError,err,error,*999)
     ENDIF

    EXITS("AdvectionDiffusion_PreSolve")
    RETURN
999 ERRORSEXITS("AdvectionDiffusion_PreSolve",err,error)
    RETURN 1
    
  END SUBROUTINE AdvectionDiffusion_PreSolve
  !   
  !================================================================================================================================
  !
  !>Update mesh position and velocity for ALE advection-diffusion problem
  SUBROUTINE AdvectionDiffusion_PreSolveALEUpdateMesh(CONTROL_LOOP,SOLVER,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(FIELD_TYPE), POINTER :: GEOMETRIC_FIELD
    TYPE(SOLVER_TYPE), POINTER :: SOLVER_ALE_DIFFUSION !<A pointer to the solvers
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    TYPE(VARYING_STRING) :: localError

    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT,ALPHA
    REAL(DP), POINTER :: MESH_DISPLACEMENT_VALUES(:)

    INTEGER(INTG) :: dof_number,TOTAL_NUMBER_OF_DOFS,NDOFS_TO_PRINT

    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS
    INTEGER(INTG) :: INPUT_TYPE,INPUT_OPTION
    REAL(DP), POINTER :: INPUT_DATA1(:)

    ENTERS("AdvectionDiffusion_PreSolveALEUpdateMesh",err,error,*999)

    NULLIFY(SOLVER_ALE_DIFFUSION)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVER_MAPPING)
    NULLIFY(EQUATIONS_SET)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%LOOP_TYPE==PROBLEM_CONTROL_TIME_LOOP_TYPE) THEN
        CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CURRENT_TIME,TIME_INCREMENT,err,error,*999)
      ELSE IF(CONTROL_LOOP%CONTROL_LOOP_LEVEL>1) THEN
        CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP%PARENT_LOOP,CURRENT_TIME,TIME_INCREMENT,err,error,*999)
      ENDIF
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification does not have a subtype set.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE,PROBLEM_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
                  & PROBLEM_NONLINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE,PROBLEM_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
                  & PROBLEM_NONLINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
                SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                    EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(1)%ptr
                    IF(ASSOCIATED(EQUATIONS_SET)) THEN
                      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
                        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
                      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<3) THEN
                        CALL FlagError("Equations set specification does not have a subtype set.",err,error,*999)
                      END IF
                      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                        CASE(EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
                          & EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
                          & EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
                          & EQUATIONS_SET_QUADRATIC_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
                          & EQUATIONS_SET_EXPONENTIAL_SOURCE_ADVECTION_DIFFUSION_SUBTYPE)
                          ! do nothing ???
                        CASE(EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
                          & EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
                          & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
                          & EQUATIONS_SET_QUADRATIC_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
                          & EQUATIONS_SET_EXP_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
                          & EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
                          & EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
                          & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
                          & EQUATIONS_SET_QUAD_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
                          & EQUATIONS_SET_EXP_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE)
                          CALL WriteString(GENERAL_OUTPUT_TYPE,"Advection-diffusion update mesh ... ",err,error,*999)
                          GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                          IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN
                            !--- First, read mesh displacement values from file

                           CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                            & NUMBER_OF_DIMENSIONS,err,error,*999)

                           INPUT_TYPE=42
                           INPUT_OPTION=2
                           NULLIFY(INPUT_DATA1)
                           !CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                            !& FIELD_VALUES_SET_TYPE,INPUT_DATA1,err,error,*999)
                           CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,INPUT_DATA1, & 
                            & NUMBER_OF_DIMENSIONS,INPUT_TYPE,INPUT_OPTION,CONTROL_LOOP%TIME_LOOP%ITERATION_NUMBER,1.0_DP, &
                            & err,error,*999)

                            NULLIFY(MESH_DISPLACEMENT_VALUES)
                            CALL FIELD_PARAMETER_SET_DATA_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, & 
                              & FIELD_MESH_DISPLACEMENT_SET_TYPE,MESH_DISPLACEMENT_VALUES,err,error,*999)
                            IF(DIAGNOSTICS1) THEN
                              NDOFS_TO_PRINT = SIZE(MESH_DISPLACEMENT_VALUES,1)
                              CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NDOFS_TO_PRINT,NDOFS_TO_PRINT,NDOFS_TO_PRINT,&
                                & MESH_DISPLACEMENT_VALUES,'(" MESH_DISPLACEMENT_VALUES = ",3(X,E13.6))','3(3(X,E13.6))', &
                                & err,error,*999)
                            ENDIF

                           CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,INPUT_DATA1, & 
                            & NUMBER_OF_DIMENSIONS,INPUT_TYPE,INPUT_OPTION,CONTROL_LOOP%TIME_LOOP%ITERATION_NUMBER,1.0_DP, &
                            & err,error,*999)

                            TOTAL_NUMBER_OF_DOFS = GEOMETRIC_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%ptr% &
                              & TOTAL_NUMBER_OF_DOFS

                            !--- Second, update geometric field
                            DO dof_number=1,TOTAL_NUMBER_OF_DOFS
                              CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(GEOMETRIC_FIELD, & 
                                & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,dof_number, & 
                                & MESH_DISPLACEMENT_VALUES(dof_number), &
                                & err,error,*999)
                            END DO
                            CALL FIELD_PARAMETER_SET_UPDATE_START(GEOMETRIC_FIELD, &
                              & FIELD_U_VARIABLE_TYPE, FIELD_VALUES_SET_TYPE,err,error,*999)
                            CALL FIELD_PARAMETER_SET_UPDATE_FINISH(GEOMETRIC_FIELD, &
                              & FIELD_U_VARIABLE_TYPE, FIELD_VALUES_SET_TYPE,err,error,*999)

                            !--- Third, use displacement values to calculate velocity values
                            ALPHA=1.0_DP/TIME_INCREMENT
                            CALL FIELD_PARAMETER_SETS_COPY(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, & 
                              & FIELD_MESH_DISPLACEMENT_SET_TYPE,FIELD_MESH_VELOCITY_SET_TYPE,ALPHA,err,error,*999)
                            CALL FIELD_PARAMETER_SET_DATA_RESTORE(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, & 
                              & FIELD_MESH_DISPLACEMENT_SET_TYPE,MESH_DISPLACEMENT_VALUES,err,error,*999)
                          ELSE
                            CALL FlagError("Geometric field is not associated.",err,error,*999)
                          ENDIF
                        CASE DEFAULT
                          localError="Equations set subtype " &
                            & //TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
                            & " is not valid for an advection-diffusion equation type of a classical field problem class."
                          CALL FlagError(localError,err,error,*999)
                      END SELECT
                    ELSE
                      CALL FlagError("Equations set is not associated.",err,error,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Solver mapping is not associated.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Solver equations is not associated.",err,error,*999)
                ENDIF
            CASE DEFAULT
              localError="Problem subtype "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
                & " is not valid for an advection-diffusion equation type of a classical field problem class."
              CALL FlagError(localError,err,error,*999)
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

    EXITS("AdvectionDiffusion_PreSolveALEUpdateMesh")
    RETURN
999 ERRORS("AdvectionDiffusion_PreSolveALEUpdateMesh",err,error)
    EXITS("AdvectionDiffusion_PreSolveALEUpdateMesh")
    RETURN 1
    
  END SUBROUTINE AdvectionDiffusion_PreSolveALEUpdateMesh
  !   
  !================================================================================================================================
  !


  SUBROUTINE AdvectionDiffusion_PreSolveStoreCurrentSoln(CONTROL_LOOP,SOLVER,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    !Local Variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER_ADVECTION_DIFFUSION !<A pointer to the solvers
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD_ADVECTION_DIFFUSION
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS_ADVECTION_DIFFUSION !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING_ADVECTION_DIFFUSION !<A pointer to the solver mapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET_ADVECTION_DIFFUSION !<A pointer to the equations set
    TYPE(VARYING_STRING) :: localError

    INTEGER(INTG) :: NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_ADVECTION_DIFFUSION
    INTEGER(INTG) :: I

    ENTERS("AdvectionDiffusion_PreSolveStoreCurrentSoln",err,error,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN

      NULLIFY(SOLVER_ADVECTION_DIFFUSION)

      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification does not have a subtype set.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
               & PROBLEM_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
               & PROBLEM_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_NONLINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
               & PROBLEM_NONLINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
              IF(SOLVER%GLOBAL_NUMBER==1) THEN
                !--- Get the dependent field of the advection-diffusion equations
                CALL WriteString(GENERAL_OUTPUT_TYPE,"Store value of advection-diffusion & 
                  & (dependent field - U variable_type) at time, t ... ",err,error,*999)
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,1,SOLVER_ADVECTION_DIFFUSION,err,error,*999)
                SOLVER_EQUATIONS_ADVECTION_DIFFUSION=>SOLVER_ADVECTION_DIFFUSION%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS_ADVECTION_DIFFUSION)) THEN
                  SOLVER_MAPPING_ADVECTION_DIFFUSION=>SOLVER_EQUATIONS_ADVECTION_DIFFUSION%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING_ADVECTION_DIFFUSION)) THEN
                    EQUATIONS_SET_ADVECTION_DIFFUSION=>SOLVER_MAPPING_ADVECTION_DIFFUSION%EQUATIONS_SETS(1)%ptr
                    IF(ASSOCIATED(EQUATIONS_SET_ADVECTION_DIFFUSION)) THEN
                      DEPENDENT_FIELD_ADVECTION_DIFFUSION=>EQUATIONS_SET_ADVECTION_DIFFUSION%DEPENDENT%DEPENDENT_FIELD
                      IF(ASSOCIATED(DEPENDENT_FIELD_ADVECTION_DIFFUSION)) THEN
                        CALL FIELD_NUMBER_OF_COMPONENTS_GET(DEPENDENT_FIELD_ADVECTION_DIFFUSION, &
                          & FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_ADVECTION_DIFFUSION,err,error,*999)
                      ELSE
                        CALL FlagError("DEPENDENT_FIELD_ADVECTION_DIFFUSIONE is not associated.",err,error,*999)
                      END IF
                    ELSE
                      CALL FlagError("Advection-diffusion equations set is not associated.",err,error,*999)
                    END IF
                  ELSE
                    CALL FlagError("Advection-diffusion solver mapping is not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("Advection-diffusion solver equations are not associated.",err,error,*999)
                END IF

                !--- Copy the current time value parameters set from diffusion-one's dependent field 
                  DO I=1,NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_ADVECTION_DIFFUSION
                    CALL Field_ParametersToFieldParametersCopy(DEPENDENT_FIELD_ADVECTION_DIFFUSION, & 
                      & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,I,DEPENDENT_FIELD_ADVECTION_DIFFUSION, & 
                      & FIELD_U_VARIABLE_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE,I,err,error,*999)
                  END DO

!                 IF(DIAGNOSTICS3) THEN
!                   NULLIFY( DUMMY_VALUES2 )
!                   CALL FIELD_PARAMETER_SET_DATA_GET(DEPENDENT_FIELD_FINITE_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
!                     & FIELD_VALUES_SET_TYPE,DUMMY_VALUES2,err,error,*999)
!                   NDOFS_TO_PRINT = SIZE(DUMMY_VALUES2,1)
!                   CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NDOFS_TO_PRINT,NDOFS_TO_PRINT,NDOFS_TO_PRINT,DUMMY_VALUES2, &
!                     & '(" DEPENDENT_FIELD_FINITE_ELASTICITY,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE = ",4(X,E13.6))',&
!                     & '4(4(X,E13.6))',err,error,*999)
!                   CALL FIELD_PARAMETER_SET_DATA_RESTORE(DEPENDENT_FIELD_FINITE_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
!                     & FIELD_VALUES_SET_TYPE,DUMMY_VALUES2,err,error,*999)
!                 ENDIF

              END IF
            CASE DEFAULT
              localError="Problem subtype "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
                & " is not valid for an advection-diffusion equation type of a classical field problem class."
            CALL FlagError(localError,err,error,*999)
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

    EXITS("AdvectionDiffusion_PreSolveStoreCurrentSoln")
    RETURN
999 ERRORS("AdvectionDiffusion_PreSolveStoreCurrentSoln",err,error)
    EXITS("AdvectionDiffusion_PreSolveStoreCurrentSoln")
    RETURN 1
    
  END SUBROUTINE AdvectionDiffusion_PreSolveStoreCurrentSoln    
  !
  !================================================================================================================================
  !
  SUBROUTINE AdvectionDiffusion_PreSolveGetSourceValue(CONTROL_LOOP,SOLVER,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    !Local Variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER_ADVECTION_DIFFUSION, SOLVER_DIFFUSION  !<A pointer to the solvers
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD_DIFFUSION, SOURCE_FIELD_ADVECTION_DIFFUSION
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS_ADVECTION_DIFFUSION, SOLVER_EQUATIONS_DIFFUSION  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING_ADVECTION_DIFFUSION, SOLVER_MAPPING_DIFFUSION !<A pointer to the solver mapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET_ADVECTION_DIFFUSION, EQUATIONS_SET_DIFFUSION !<A pointer to the equations set
    TYPE(VARYING_STRING) :: localError

    INTEGER(INTG) :: NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_DIFFUSION,NUMBER_OF_COMPONENTS_SOURCE_FIELD_ADVECTION_DIFFUSION
    INTEGER(INTG) :: I


    ENTERS("AdvectionDiffusion_PreSolveGetSourceValue",err,error,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN

      NULLIFY(SOLVER_ADVECTION_DIFFUSION)
      NULLIFY(SOLVER_DIFFUSION)

      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification does not have a subtype set.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
               & PROBLEM_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
               & PROBLEM_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_NONLINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
               & PROBLEM_NONLINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
              IF(SOLVER%GLOBAL_NUMBER==1) THEN
                !--- Get the dependent field of the diffusion equations
                CALL WriteString(GENERAL_OUTPUT_TYPE,"Update advection-diffusion source field ... ",err,error,*999)
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,2,SOLVER_DIFFUSION,err,error,*999)
                SOLVER_EQUATIONS_DIFFUSION=>SOLVER_DIFFUSION%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS_DIFFUSION)) THEN
                  SOLVER_MAPPING_DIFFUSION=>SOLVER_EQUATIONS_DIFFUSION%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING_DIFFUSION)) THEN
                    EQUATIONS_SET_DIFFUSION=>SOLVER_MAPPING_DIFFUSION%EQUATIONS_SETS(1)%ptr
                    IF(ASSOCIATED(EQUATIONS_SET_DIFFUSION)) THEN
                      DEPENDENT_FIELD_DIFFUSION=>EQUATIONS_SET_DIFFUSION%DEPENDENT%DEPENDENT_FIELD
                      IF(ASSOCIATED(DEPENDENT_FIELD_DIFFUSION)) THEN
                        CALL FIELD_NUMBER_OF_COMPONENTS_GET(DEPENDENT_FIELD_DIFFUSION, &
                          & FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_DIFFUSION,err,error,*999)
                      ELSE
                        CALL FlagError("DEPENDENT_FIELD_DIFFUSION is not associated.",err,error,*999)
                      END IF
                    ELSE
                      CALL FlagError("Diffusion equations set is not associated.",err,error,*999)
                    END IF
                  ELSE
                    CALL FlagError("Diffusion solver mapping is not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("Diffusion solver equations are not associated.",err,error,*999)
                END IF


                !--- Get the source field for the advection-diffusion equations
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,1,SOLVER_ADVECTION_DIFFUSION,err,error,*999)
                SOLVER_EQUATIONS_ADVECTION_DIFFUSION=>SOLVER_ADVECTION_DIFFUSION%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS_ADVECTION_DIFFUSION)) THEN
                  SOLVER_MAPPING_ADVECTION_DIFFUSION=>SOLVER_EQUATIONS_ADVECTION_DIFFUSION%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING_ADVECTION_DIFFUSION)) THEN
                    EQUATIONS_SET_ADVECTION_DIFFUSION=>SOLVER_MAPPING_ADVECTION_DIFFUSION%EQUATIONS_SETS(1)%ptr
                    IF(ASSOCIATED(EQUATIONS_SET_ADVECTION_DIFFUSION)) THEN
                      SOURCE_FIELD_ADVECTION_DIFFUSION=>EQUATIONS_SET_ADVECTION_DIFFUSION%SOURCE%SOURCE_FIELD
                      IF(ASSOCIATED(SOURCE_FIELD_ADVECTION_DIFFUSION)) THEN
                        CALL FIELD_NUMBER_OF_COMPONENTS_GET(SOURCE_FIELD_ADVECTION_DIFFUSION, & 
                          & FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS_SOURCE_FIELD_ADVECTION_DIFFUSION,err,error,*999)
                      ELSE
                        CALL FlagError("SOURCE_FIELD_ADVECTION_DIFFUSION is not associated.",err,error,*999)
                      END IF
                    ELSE
                      CALL FlagError("Advection-diffusion equations set is not associated.",err,error,*999)
                    END IF
                  ELSE
                    CALL FlagError("Advection-diffusion solver mapping is not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("Advection-diffusion solver equations are not associated.",err,error,*999)
                END IF

                !--- Copy the result from diffusion's dependent field to advection-diffusion's source field
                IF(NUMBER_OF_COMPONENTS_SOURCE_FIELD_ADVECTION_DIFFUSION==NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_DIFFUSION) THEN
                  DO I=1,NUMBER_OF_COMPONENTS_SOURCE_FIELD_ADVECTION_DIFFUSION
                    CALL Field_ParametersToFieldParametersCopy(DEPENDENT_FIELD_DIFFUSION, & 
                      & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,I,SOURCE_FIELD_ADVECTION_DIFFUSION, & 
                      & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,I,err,error,*999)
                  END DO
                ELSE
                  localError="Number of components of diffusion dependent field "// &
                    & "is not consistent with advection-diffusion equation source field."
                  CALL FlagError(localError,err,error,*999)
                END IF

!                 IF(DIAGNOSTICS3) THEN
!                   NULLIFY( DUMMY_VALUES2 )
!                   CALL FIELD_PARAMETER_SET_DATA_GET(DEPENDENT_FIELD_FINITE_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
!                     & FIELD_VALUES_SET_TYPE,DUMMY_VALUES2,err,error,*999)
!                   NDOFS_TO_PRINT = SIZE(DUMMY_VALUES2,1)
!                   CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NDOFS_TO_PRINT,NDOFS_TO_PRINT,NDOFS_TO_PRINT,DUMMY_VALUES2, &
!                     & '(" DEPENDENT_FIELD_FINITE_ELASTICITY,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE = ",4(X,E13.6))',&
!                     & '4(4(X,E13.6))',err,error,*999)
!                   CALL FIELD_PARAMETER_SET_DATA_RESTORE(DEPENDENT_FIELD_FINITE_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
!                     & FIELD_VALUES_SET_TYPE,DUMMY_VALUES2,err,error,*999)
!                 ENDIF

              END IF
            CASE DEFAULT
              localError="Problem subtype "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
                & " is not valid for an advection-diffusion equation type of a classical field problem class."
            CALL FlagError(localError,err,error,*999)
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

    EXITS("AdvectionDiffusion_PreSolveGetSourceValue")
    RETURN
999 ERRORS("AdvectionDiffusion_PreSolveGetSourceValue",err,error)
    EXITS("AdvectionDiffusion_PreSolveGetSourceValue")
    RETURN 1
    
  END SUBROUTINE AdvectionDiffusion_PreSolveGetSourceValue
  !   
  !================================================================================================================================
  !
  !>Update independent field (velocity) for advection-diffusion pre solve
  SUBROUTINE AdvectionDiffusion_PreSolveUpdateInputData(CONTROL_LOOP,SOLVER,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    TYPE(EquationsType), POINTER :: EQUATIONS

    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS
    INTEGER(INTG) :: INPUT_TYPE,INPUT_OPTION
    REAL(DP), POINTER :: INPUT_DATA1(:)
    INTEGER(INTG) :: PROBLEM_SUBTYPE

    ENTERS("AdvectionDiffusion_PreSolveUpdateInputData",err,error,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification does not have a subtype set.",err,error,*999)
          END IF
          PROBLEM_SUBTYPE=CONTROL_LOOP%PROBLEM%SPECIFICATION(3)
          IF(PROBLEM_SUBTYPE==PROBLEM_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
            & PROBLEM_SUBTYPE==PROBLEM_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
            & PROBLEM_SUBTYPE==PROBLEM_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
            & PROBLEM_SUBTYPE==PROBLEM_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
            & PROBLEM_SUBTYPE==PROBLEM_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE) THEN
            
            CALL WriteString(GENERAL_OUTPUT_TYPE,"Read input data... ",err,error,*999)
            SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
            IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
              SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
              EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
              IF(ASSOCIATED(EQUATIONS)) THEN
                EQUATIONS_SET=>equations%equationsSet
                IF(ASSOCIATED(EQUATIONS_SET)) THEN
                  
                  CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_DIMENSIONS,err,error,*999)
                  
                  INPUT_TYPE=1
                  INPUT_OPTION=1
                  NULLIFY(INPUT_DATA1)
                  CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                    & FIELD_VALUES_SET_TYPE,INPUT_DATA1,err,error,*999)
                  CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,INPUT_DATA1, & 
                    & NUMBER_OF_DIMENSIONS,INPUT_TYPE,INPUT_OPTION,CONTROL_LOOP%TIME_LOOP%ITERATION_NUMBER,1.0_DP, &
                    & err,error,*999)
                  
                ELSE
                  CALL FlagError("Equations set is not associated.",err,error,*999)
                END IF
                
              ELSE
                CALL FlagError("Equations are not associated.",err,error,*999)
              END IF
            ELSE
              CALL FlagError("Solver equations are not associated.",err,error,*999)
            END IF
          ELSE IF(PROBLEM_SUBTYPE==PROBLEM_NONLINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
            & PROBLEM_SUBTYPE==PROBLEM_NONLINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE) THEN
            CALL FlagError("Not implemented.",err,error,*999)
          ENDIF
          
          CALL FIELD_PARAMETER_SET_UPDATE_START(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
            & FIELD_VALUES_SET_TYPE,err,error,*999)
          CALL FIELD_PARAMETER_SET_UPDATE_FINISH(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
            & FIELD_VALUES_SET_TYPE,err,error,*999)
          
        ELSE
          CALL FlagError("Problem is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",err,error,*999)
    ENDIF
    
    EXITS("AdvectionDiffusion_PreSolveUpdateInputData")
    RETURN
999 ERRORS("AdvectionDiffusion_PreSolveUpdateInputData",err,error)
    EXITS("AdvectionDiffusion_PreSolveUpdateInputData")
    RETURN 1
    
  END SUBROUTINE AdvectionDiffusion_PreSolveUpdateInputData

  !
  !================================================================================================================================
  !
  !Update the boundary conditions
  SUBROUTINE ADVECTION_DIFFUSION_PRE_SOLVE_UPDATE_BC(CONTROL_LOOP,SOLVER,err,error,*)
    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    TYPE(EquationsType), POINTER :: EQUATIONS
    TYPE(VARYING_STRING) :: localError
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT

!\todo: Reduce number of variable used
    INTEGER(INTG) :: BOUNDARY_CONDITION_CHECK_VARIABLE,node_idx
    INTEGER(INTG) :: NUMBER_OF_COMPONENTS
    INTEGER(INTG) :: local_ny,global_ny
    REAL(DP), POINTER :: BOUNDARY_VALUES(:)
    INTEGER(INTG), POINTER :: BOUNDARY_NODES(:)

    ENTERS("ADVECTION_DIFFUSION_PRE_SOLVE_UPDATE_BC",err,error,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CURRENT_TIME,TIME_INCREMENT,err,error,*999)
      WRITE (*,*) CONTROL_LOOP%TIME_LOOP%ITERATION_NUMBER
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification does not have a subtype set.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
          CASE (PROBLEM_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE,PROBLEM_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE,&
            & PROBLEM_NONLINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE,PROBLEM_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE,&
            & PROBLEM_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE,PROBLEM_NONLINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
          SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
          IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
            SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
            EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
            IF(ASSOCIATED(EQUATIONS)) THEN
              EQUATIONS_SET=>equations%equationsSet
              IF(ASSOCIATED(EQUATIONS_SET)) THEN
                BOUNDARY_CONDITIONS=>SOLVER_EQUATIONS%BOUNDARY_CONDITIONS
                IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
                  DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                  IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                    CALL Field_VariableGet(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VARIABLE,err,error,*999)
                    IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                      CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,FIELD_VARIABLE,BOUNDARY_CONDITIONS_VARIABLE, &
                        & err,error,*999)
                      IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                        CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                         & NUMBER_OF_COMPONENTS,err,error,*999)
                        NULLIFY(BOUNDARY_VALUES)
                        NULLIFY(BOUNDARY_NODES)
                        CALL FLUID_MECHANICS_IO_READ_BOUNDARY_CONDITIONS_ITERATION(SOLVER_LINEAR_TYPE,BOUNDARY_VALUES, &
                          & BOUNDARY_NODES,NUMBER_OF_COMPONENTS,BOUNDARY_CONDITION_FIXED,CONTROL_LOOP%TIME_LOOP%INPUT_NUMBER, &
                          & CONTROL_LOOP%TIME_LOOP%ITERATION_NUMBER,err,error,*999)
                        WRITE(*,*) SIZE(BOUNDARY_VALUES)
                        DO node_idx=1,SIZE(BOUNDARY_VALUES)
                          !Default to version 1 of each node derivative
                          CALL FIELD_COMPONENT_DOF_GET_USER_NODE(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                            & 1,NO_GLOBAL_DERIV,BOUNDARY_NODES(node_idx), &
                            & 1,local_ny,global_ny,err,error,*999)
                          BOUNDARY_CONDITION_CHECK_VARIABLE=BOUNDARY_CONDITIONS_VARIABLE% &
                            & CONDITION_TYPES(local_ny)
                          IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_FIXED) THEN
                            CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                              & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny, &
                              & BOUNDARY_VALUES(node_idx),err,error,*999)
                          END IF
                        ENDDO
                      ELSE
                        CALL FlagError("Boundary condition variable is not associated.",err,error,*999)
                      END IF
                    ELSE
                      CALL FlagError("Dependent field variable is not associated.",err,error,*999)
                    END IF
                  ELSE
                    CALL FlagError("Equations set dependent variable is not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("Boundary conditions are not associated.",err,error,*999)
                END IF
              ELSE
                CALL FlagError("Equations set is not associated.",err,error,*999)
              END IF
            ELSE
              CALL FlagError("Equations are not associated.",err,error,*999)
            END IF
          ELSE
            CALL FlagError("Solver equations are not associated.",err,error,*999)
          END IF
          CALL FIELD_PARAMETER_SET_UPDATE_START(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
            & FIELD_VALUES_SET_TYPE,err,error,*999)
          CALL FIELD_PARAMETER_SET_UPDATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
            & FIELD_VALUES_SET_TYPE,err,error,*999)
          CASE DEFAULT
            localError="Problem subtype "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
              & " is not valid for an advection-diffusion equation of a classical field problem class."
          CALL FlagError(localError,err,error,*999)
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

    EXITS("ADVECTION_DIFFUSION_PRE_SOLVE_UPDATE_BC")
    RETURN
999 ERRORSEXITS("ADVECTION_DIFFUSION_PRE_SOLVE_UPDATE_BC",err,error)
    RETURN 1

  END SUBROUTINE ADVECTION_DIFFUSION_PRE_SOLVE_UPDATE_BC
  !
  !================================================================================================================================
  !
  SUBROUTINE AdvectionDiffusion_PostSolve(solver,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop
    TYPE(PROBLEM_TYPE), POINTER :: problem
    TYPE(VARYING_STRING) :: localError

    ENTERS("AdvectionDiffusion_PostSolve",err,error,*999)
    
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    IF(.NOT.ALLOCATED(problem%specification)) CALL FlagError("Problem specification is not allocated.",err,error,*999)
    IF(SIZE(problem%specification,1)<3) &
      & CALL FlagError("Problem specification must have at least three entries for an advection-diffusion problem.",err,error,*999)
 
    SELECT CASE(problem%specification(3))
    CASE(PROBLEM_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE,PROBLEM_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
      & PROBLEM_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE,PROBLEM_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
      CALL ADVECTION_DIFFUSION_POST_SOLVE_OUTPUT_DATA(controlLoop,solver,err,error,*999)
    CASE(PROBLEM_NONLINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE,PROBLEM_NONLINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
      & PROBLEM_NO_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE,PROBLEM_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE)
      ! do nothing ???
      !CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The third problem specification of "// &
        & TRIM(NumberToVString(problem%specification(3),"*",err,error))// &
        & " is not valid for an advection-diffusion type of a classical field problem."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("AdvectionDiffusion_PostSolve")
    RETURN
999 ERRORSEXITS("AdvectionDiffusion_PostSolve",err,error)
    RETURN 1
    
  END SUBROUTINE AdvectionDiffusion_PostSolve
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE ADVECTION_DIFFUSION_POST_SOLVE_OUTPUT_DATA(CONTROL_LOOP,SOLVER,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER!<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    TYPE(VARYING_STRING) :: localError

    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT
    INTEGER(INTG) :: EQUATIONS_SET_IDX,CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER

    CHARACTER(14) :: FILE
    CHARACTER(14) :: OUTPUT_FILE

    ENTERS("ADVECTION_DIFFUSION_POST_SOLVE_OUTPUT_DATA",err,error,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification does not have a subtype set.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE,PROBLEM_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE)
              CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CURRENT_TIME,TIME_INCREMENT,err,error,*999)
              SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
              IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                  !Make sure the equations sets are up to date
                  DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                    EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%ptr

                    CURRENT_LOOP_ITERATION=CONTROL_LOOP%TIME_LOOP%ITERATION_NUMBER
                    OUTPUT_ITERATION_NUMBER=CONTROL_LOOP%TIME_LOOP%OUTPUT_NUMBER

                    IF(OUTPUT_ITERATION_NUMBER/=0) THEN
                      IF(CONTROL_LOOP%TIME_LOOP%CURRENT_TIME<=CONTROL_LOOP%TIME_LOOP%STOP_TIME) THEN
                        IF(CURRENT_LOOP_ITERATION<10) THEN
                          WRITE(OUTPUT_FILE,'("TIME_STEP_000",I0)') CURRENT_LOOP_ITERATION
                        ELSE IF(CURRENT_LOOP_ITERATION<100) THEN
                          WRITE(OUTPUT_FILE,'("TIME_STEP_00",I0)') CURRENT_LOOP_ITERATION
                        ELSE IF(CURRENT_LOOP_ITERATION<1000) THEN
                          WRITE(OUTPUT_FILE,'("TIME_STEP_0",I0)') CURRENT_LOOP_ITERATION
                        ELSE IF(CURRENT_LOOP_ITERATION<10000) THEN
                          WRITE(OUTPUT_FILE,'("TIME_STEP_",I0)') CURRENT_LOOP_ITERATION
                        END IF
                        FILE=OUTPUT_FILE
!                        FILE="TRANSIENT_OUTPUT"
!                         METHOD="FORTRAN"
!                         EXPORT_FIELD=.TRUE.
!                         IF(EXPORT_FIELD) THEN          
!                          IF(MOD(CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER)==0)  THEN   
                            CALL WriteString(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
                            CALL WriteString(GENERAL_OUTPUT_TYPE,"Now export fields... ",err,error,*999)
                          CALL FLUID_MECHANICS_IO_WRITE_CMGUI(EQUATIONS_SET%REGION,EQUATIONS_SET%GLOBAL_NUMBER,FILE, &
                              & err,error,*999)
                            CALL WriteString(GENERAL_OUTPUT_TYPE,OUTPUT_FILE,err,error,*999)
                            CALL WriteString(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
!                           ENDIF
!                         ENDIF 
                      ENDIF 
                    ENDIF
                  ENDDO
                ENDIF
              ENDIF
            CASE(PROBLEM_NONLINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE)
              ! do nothing ???
              CALL FlagError("Not implemented.",err,error,*999)
            CASE DEFAULT
              localError="Problem subtype "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
                & " is not valid for an advection-diffusion equation type of a classical field problem class."
            CALL FlagError(localError,err,error,*999)
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

    EXITS("ADVECTION_DIFFUSION_POST_SOLVE_OUTPUT_DATA")
    RETURN
999 ERRORSEXITS("ADVECTION_DIFFUSION_POST_SOLVE_OUTPUT_DATA",err,error)
    RETURN 1
  END SUBROUTINE ADVECTION_DIFFUSION_POST_SOLVE_OUTPUT_DATA
  !
  !================================================================================================================================
  !
END MODULE ADVECTION_DIFFUSION_EQUATION_ROUTINES

