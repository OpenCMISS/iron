!> \file
!> \author Soroush Safaei
!> \brief This module handles pure advection equation routines.
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
!> Contributor(s): Chris Bradley
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

!>This module handles pure advection equation routines.
MODULE ADVECTION_EQUATION_ROUTINES

  USE ANALYTIC_ANALYSIS_ROUTINES
  USE BaseRoutines
  USE BasisRoutines
  USE BasisAccessRoutines
  USE BOUNDARY_CONDITIONS_ROUTINES
  USE Constants
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
  USE Maths
  USE MatrixVector
  USE PROBLEM_CONSTANTS
  USE Strings
  USE SOLVER_ROUTINES
  USE SolverAccessRoutines
  USE Timer
  USE Types

#include "macros.h"

  IMPLICIT NONE

  PUBLIC Advection_EquationsSetSetup
  
  PUBLIC Advection_EquationsSetSolutionMethodSet
  
  PUBLIC Advection_EquationsSetSpecificationSet
  
  PUBLIC Advection_ProblemSpecificationSet
  
  PUBLIC ADVECTION_EQUATION_PROBLEM_SETUP
  
  PUBLIC ADVECTION_EQUATION_FINITE_ELEMENT_CALCULATE
  
  PUBLIC Advection_PreSolve
  
  PUBLIC ADVECTION_PRE_SOLVE_UPDATE_BC
  
  PUBLIC Advection_PostSolve
  
  PUBLIC Advection_Couple1D0D


CONTAINS

  !
  !================================================================================================================================
  !

  !>Sets up the diffusion equation type of a classical field equations set class.
  SUBROUTINE Advection_EquationsSetSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Advection_EquationsSetSetup",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Advection equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%specification(3))
      CASE(EQUATIONS_SET_ADVECTION_SUBTYPE)
        CALL Advection_EquationsSetLinearSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*999)
      CASE DEFAULT
        localError="The third equations set specification of "// &
          & TRIM(NumberToVString(EQUATIONS_SET%specification(3),"*",err,error))// &
          & " is not valid for an advection type of a classical field equation set."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
       
    EXITS("Advection_EquationsSetSetup")
    RETURN
999 ERRORSEXITS("Advection_EquationsSetSetup",err,error)
    RETURN 1
    
  END SUBROUTINE Advection_EquationsSetSetup

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for an advection equation type of an classical field equations set class.
  SUBROUTINE Advection_EquationsSetSolutionMethodSet(EQUATIONS_SET,SOLUTION_METHOD,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Advection_EquationsSetSolutionMethodSet",err,error,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a advection type equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
      CASE(EQUATIONS_SET_ADVECTION_SUBTYPE)
        SELECT CASE(SOLUTION_METHOD)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          EQUATIONS_SET%SOLUTION_METHOD=EQUATIONS_SET_FEM_SOLUTION_METHOD
        CASE DEFAULT
          localError="The specified solution method of "//TRIM(NumberToVString(SOLUTION_METHOD,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The third equations set specification of "// &
          & TRIM(NumberToVString(EQUATIONS_SET%specification(3),"*",err,error))// &
          & " is not valid for an advection type of an classical field equations set."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
       
    EXITS("Advection_EquationsSetSolutionMethodSet")
    RETURN
999 ERRORS("Advection_EquationsSetSolutionMethodSet",err,error)
    EXITS("Advection_EquationsSetSolutionMethodSet")
    RETURN 1
    
  END SUBROUTINE Advection_EquationsSetSolutionMethodSet
  
  !
  !================================================================================================================================
  !

  !>Sets the equation specification for an advection type of a classical field equations set.
  SUBROUTINE Advection_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to set the specification for
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations set specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: subtype
    
    ENTERS("Advection_EquationsSetSpecificationSet",err,error,*999)
    
    IF(ASSOCIATED(equationsSet)) THEN
      IF(SIZE(specification,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Laplace type equations set.", &
          & err,error,*999)
      END IF
      subtype=specification(3)
      SELECT CASE(subtype)
      CASE(EQUATIONS_SET_ADVECTION_SUBTYPE)
        !ok
      CASE DEFAULT
        localError="The third equations set specification of "//TRIM(NumberToVstring(subtype,"*",err,error))// &
          & " is not valid for an advection type of a classical field equations set."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Set full specification
      IF(ALLOCATED(equationsSet%specification)) THEN
        CALL FlagError("Equations set specification is already allocated.",err,error,*999)
      ELSE
        ALLOCATE(equationsSet%specification(3),stat=err)
        IF(err/=0) CALL FlagError("Could not allocate equations set specification.",err,error,*999)
      END IF
      equationsSet%specification(1:3)=[EQUATIONS_SET_CLASSICAL_FIELD_CLASS,EQUATIONS_SET_ADVECTION_EQUATION_TYPE,subtype]
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    END IF
     
    EXITS("Advection_EquationsSetSpecificationSet")
    RETURN
999 ERRORSEXITS("Advection_EquationsSetSpecificationSet",err,error)
    RETURN 1
    
  END SUBROUTINE Advection_EquationsSetSpecificationSet
  
  !
  !================================================================================================================================
  !
  
  !>Sets up the advection equation.
  SUBROUTINE Advection_EquationsSetLinearSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    INTEGER(INTG) :: GEOMETRIC_MESH_COMPONENT,GEOMETRIC_SCALING_TYPE,GEOMETRIC_COMPONENT_NUMBER,NUMBER_OF_DIMENSIONS
    INTEGER(INTG) :: DEPENDENT_FIELD_NUMBER_OF_VARIABLES,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,component_idx
    INTEGER(INTG) :: INDEPENDENT_FIELD_NUMBER_OF_VARIABLES,INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS
    TYPE(DECOMPOSITION_TYPE), POINTER :: GEOMETRIC_DECOMPOSITION
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EQUATIONS_SET_MATERIALS_TYPE), POINTER :: EQUATIONS_MATERIALS
    TYPE(VARYING_STRING) :: localError

    ENTERS("Advection_EquationsSetLinearSetup",err,error,*999)

    NULLIFY(EQUATIONS)
    NULLIFY(vectorMapping)
    NULLIFY(vectorMatrices)
    NULLIFY(GEOMETRIC_DECOMPOSITION)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for an advection equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
      CASE(EQUATIONS_SET_ADVECTION_SUBTYPE)
        SELECT CASE(EQUATIONS_SET_SETUP%SETUP_TYPE)
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            CALL Advection_EquationsSetSolutionMethodSet(EQUATIONS_SET, &
              & EQUATIONS_SET_FEM_SOLUTION_METHOD,err,error,*999)
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for an advection equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT

        CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
          SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
          CASE(EQUATIONS_SET_ADVECTION_SUBTYPE)
            !Do nothing
          END SELECT
        !-----------------------------------------------------------------
        ! D e p e n d e n t   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            DEPENDENT_FIELD_NUMBER_OF_VARIABLES=2
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
              !Create the auto created dependent field
              !start field creation with name 'Concentration'
              CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION, &
                & EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,err,error,*999)
              !start creation of a new field
              CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_GENERAL_TYPE,err,error,*999)
              !label the field
              CALL FIELD_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,"Concentration",err,error,*999)
              !define new created field to be dependent
              CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                & FIELD_DEPENDENT_TYPE,err,error,*999)
              !look for decomposition rule already defined
              CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION, &
                & err,error,*999)
              !apply decomposition rule found on new created field
              CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                & GEOMETRIC_DECOMPOSITION,err,error,*999)
              !point new field to geometric field
              CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY% &
                & GEOMETRIC_FIELD,err,error,*999)
              !set number of variables to 2 (U,delU/delN)
              CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                & DEPENDENT_FIELD_NUMBER_OF_VARIABLES,err,error,*999)
              CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,[FIELD_U_VARIABLE_TYPE, &
                & FIELD_DELUDELN_VARIABLE_TYPE],err,error,*999)
              !set dimension
              CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              !set data type
              CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & NUMBER_OF_DIMENSIONS,err,error,*999)
              !number of components for U,delU/delN (C)
              DEPENDENT_FIELD_NUMBER_OF_COMPONENTS=1
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                & FIELD_U_VARIABLE_TYPE,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                & FIELD_DELUDELN_VARIABLE_TYPE,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
              CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, & 
                & DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,GEOMETRIC_MESH_COMPONENT,err,error,*999)
              !Default to the geometric interpolation setup
              DO component_idx=1,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, & 
                  & FIELD_U_VARIABLE_TYPE,component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
              END DO
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              !Specify fem solution method
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                DO component_idx=1,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,component_idx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_DELUDELN_VARIABLE_TYPE,component_idx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                END DO
                CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE, &
                  & err,error,*999)
                CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE, &
                  & err,error,*999)
              CASE DEFAULT
                localError="The solution method of " &
                  & //TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ELSE 
              !Check the user specified field advection equations
              CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
              CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
              CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,DEPENDENT_FIELD_NUMBER_OF_VARIABLES,err,error,*999)
              CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE], &
                & err,error,*999)
              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, & 
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE, &
                  & err,error,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & NUMBER_OF_DIMENSIONS,err,error,*999)
              DEPENDENT_FIELD_NUMBER_OF_COMPONENTS=1
              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                & DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, & 
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              CASE DEFAULT
                localError="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD, &
                  & "*",err,error))//" is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ENDIF
          !Specify finish action       
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
              CALL FIELD_CREATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,err,error,*999)
              CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                 & FIELD_BOUNDARY_CONDITIONS_SET_TYPE,err,error,*999)
            ENDIF
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for an advection equation"
            CALL FlagError(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! M a t e r i a l s   f i e l d 
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          !Specify start action
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
            IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
              IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                !Create the auto created materials field
                !start field creation with name 'MATERIAL_FIELD'
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION, & 
                  & EQUATIONS_SET%MATERIALS%MATERIALS_FIELD,err,error,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_MATERIAL_TYPE,err,error,*999)
                !label the field
                CALL FIELD_LABEL_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,"Materials Field",err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_INDEPENDENT_TYPE, &
                  & err,error,*999)
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION, & 
                  & err,error,*999)
                !apply decomposition rule found on new created field
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%MATERIALS%MATERIALS_FIELD, & 
                  & GEOMETRIC_DECOMPOSITION,err,error,*999)
                !point new field to geometric field
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,EQUATIONS_SET%GEOMETRY% &
                  & GEOMETRIC_FIELD,err,error,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,1,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE,"Materials",err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,err,error,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD, & 
                  & FIELD_U_VARIABLE_TYPE,1,GEOMETRIC_COMPONENT_NUMBER,err,error,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_COMPONENT_NUMBER,err,error,*999)
                CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                !Default the field scaling to that of the geometric field
                CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE, & 
                  & err,error,*999)
                CALL FIELD_SCALING_TYPE_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
              ELSE
                !Check the user specified field
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_MATERIAL_TYPE,err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,1,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, & 
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE, & 
                  & err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1,err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Equations set materials is not associated.",err,error,*999)
            END IF
          !Specify start action
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
            IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
              IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                !Finish creating the materials field
                CALL FIELD_CREATE_FINISH(EQUATIONS_MATERIALS%MATERIALS_FIELD,err,error,*999)
                !Set the default values for the materials field (Diffusivity)
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,1,1.0_DP,err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Equations set materials is not associated.",err,error,*999)
            ENDIF
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*", & 
              & err,error))//" for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*", & 
              & err,error))//" is invalid for Navier-Stokes equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! I n d e p e n d e n t   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_INDEPENDENT_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          !Set start action
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            INDEPENDENT_FIELD_NUMBER_OF_VARIABLES=1
            INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS=1
            !Create the auto created independent field
            IF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
              !start field creation with name 'INDEPENDENT_FIELD'
              CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION, &
                & EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,err,error,*999)
                !start creation of a new field
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                !label the field
                CALL FIELD_LABEL_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,"Independent Field",err,error,*999)
                !define new created field to be independent
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                  & FIELD_INDEPENDENT_TYPE,err,error,*999)
                !look for decomposition rule already defined
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,err,error,*999)
                !apply decomposition rule found on new created field
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                  & GEOMETRIC_DECOMPOSITION,err,error,*999)
                !point new field to geometric field
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,EQUATIONS_SET% & 
                  & GEOMETRY%GEOMETRIC_FIELD,err,error,*999)
                !set number of variables to 1 (1 for U)
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                  & INDEPENDENT_FIELD_NUMBER_OF_VARIABLES,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,[FIELD_U_VARIABLE_TYPE], &
                  & err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                !calculate number of components with one component for each dimension
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, & 
                  & FIELD_U_VARIABLE_TYPE,INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, & 
                  & 1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                !Default to the geometric interpolation setup
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, & 
                  & FIELD_U_VARIABLE_TYPE,INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                !Specify fem solution method
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE, &
                    & err,error,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE, &
                    & err,error,*999)
                CASE(EQUATIONS_SET_NODAL_SOLUTION_METHOD)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE, &
                    & err,error,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE, &
                    & err,error,*999)
                CASE DEFAULT
                  localError="The solution method of " &
                    & //TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// " is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT 
              ELSE
                !Check the user specified field- Characteristic equation
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              ENDIF    
            !Specify finish action
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
                CALL FIELD_CREATE_FINISH(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,err,error,*999)
              ENDIF
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
                & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a standard Navier-Stokes fluid"
              CALL FlagError(localError,err,error,*999)
            END SELECT
        !-----------------------------------------------------------------
        ! E q u a t i o n s    t y p e
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
              IF(EQUATIONS_SET%specification(3)==EQUATIONS_SET_ADVECTION_SUBTYPE) THEN 
                CALL Equations_CreateStart(EQUATIONS_SET,EQUATIONS,err,error,*999)
                CALL Equations_LinearityTypeSet(EQUATIONS,EQUATIONS_LINEAR,err,error,*999)
                CALL Equations_TimeDependenceTypeSet(EQUATIONS,EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
              ELSE
                CALL FlagError("Equations set subtype not valid.",err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Equations set dependent field has not been finished.",err,error,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              IF(EQUATIONS_SET%specification(3)==EQUATIONS_SET_ADVECTION_SUBTYPE) THEN 
                !Finish the equations
                CALL EquationsSet_EquationsGet(EQUATIONS_SET,EQUATIONS,err,error,*999)
                CALL Equations_CreateFinish(EQUATIONS,err,error,*999)
                NULLIFY(vectorEquations)
                CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
                !Create the equations mapping.
                CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_DELUDELN_VARIABLE_TYPE,vectorMapping,err,error,*999)
                CALL EquationsMapping_DynamicMatricesSet(vectorMapping,.TRUE.,.TRUE.,err,error,*999)
                CALL EquationsMapping_DynamicVariableTypeSet(vectorMapping,FIELD_U_VARIABLE_TYPE,err,error,*999)
                CALL EquationsMapping_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE,err,error,*999)
                CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
                !Create the equations matrices
                CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
                !Set up matrix storage and structure
                SELECT CASE(EQUATIONS%sparsityType)
                CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                  CALL EquationsMatrices_LinearStorageTypeSet(vectorMatrices, &
                    & [DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE,DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE],err,error,*999)
                CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                  CALL EquationsMatrices_DynamicStorageTypeSet(vectorMatrices, &
                    & [DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE, & 
                    & DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE],err,error,*999)
                  CALL EquationsMatrices_DynamicStructureTypeSet(vectorMatrices, &
                    & [EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_FEM_STRUCTURE],err,error,*999)
                CASE DEFAULT
                  localError="The equations matrices sparsity type of "// &
                    & TRIM(NumberToVString(EQUATIONS%sparsityType,"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
                CALL EquationsMatrices_VectorCreateFinish(vectorMatrices,err,error,*999)
              ENDIF
            CASE DEFAULT
                localError="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                & " is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for an advection equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
            & " is invalid for an advection equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The thrid equations set specification of "// &
          & TRIM(NumberToVString(EQUATIONS_SET%specification(3),"*",err,error))// &
          & " does not equal a advection equation set."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
       
    EXITS("Advection_EquationsSetLinearSetup")
    RETURN
999 ERRORSEXITS("Advection_EquationsSetLinearSetup",err,error)
    RETURN 1
    
  END SUBROUTINE Advection_EquationsSetLinearSetup
  
  !
  !================================================================================================================================
  !

  !>Sets the problem specification for an advection problem.
  SUBROUTINE Advection_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: problem !<A pointer to the problem to set the problem specification for
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: problemSubtype
    
    ENTERS("Advection_ProblemSpecificationSet",err,error,*999)
    
    IF(ASSOCIATED(PROBLEM)) THEN
      IF(SIZE(problemSpecification,1)==3) THEN
        problemSubtype=problemSpecification(3)
        SELECT CASE(problemSubtype)
        CASE(PROBLEM_ADVECTION_SUBTYPE)
          !ok
        CASE DEFAULT
          localError="The third problem specification of "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
            & " is not valid for a advection type of a classical field problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        IF(ALLOCATED(problem%specification)) THEN
          CALL FlagError("Problem specification is already allocated.",err,error,*999)
        ELSE
          ALLOCATE(problem%specification(3),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
        END IF
        problem%specification(1:3)=[PROBLEM_CLASSICAL_FIELD_CLASS,PROBLEM_ADVECTION_EQUATION_TYPE,problemSubtype]
      ELSE
        CALL FlagError("Advection problem specification must have three entries.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    END IF
        
    EXITS("Advection_ProblemSpecificationSet")
    RETURN
999 ERRORSEXITS("Advection_ProblemSpecificationSet",err,error)
    RETURN 1
    
  END SUBROUTINE Advection_ProblemSpecificationSet

  !
  !================================================================================================================================
  !
 
  !>Sets up the diffusion problem.
  SUBROUTINE ADVECTION_EQUATION_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("ADVECTION_EQUATION_PROBLEM_SETUP",err,error,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(.NOT.ALLOCATED(PROBLEM%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(PROBLEM%SPECIFICATION,1)<3) THEN
        CALL FlagError("Problem specification must have three entries for a Laplace problem.",err,error,*999)
      END IF
      SELECT CASE(PROBLEM%SPECIFICATION(3))
      CASE(PROBLEM_ADVECTION_SUBTYPE)
        CALL ADVECTION_EQUATION_PROBLEM_LINEAR_SETUP(PROBLEM,PROBLEM_SETUP,err,error,*999)
      CASE DEFAULT
        localError="The third problem specification of "// &
          & TRIM(NumberToVString(PROBLEM%SPECIFICATION(3),"*",err,error))// &
          & " is not valid for an advection type of a classical field problem."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    ENDIF
       
    EXITS("ADVECTION_EQUATION_PROBLEM_SETUP")
    RETURN
999 ERRORSEXITS("ADVECTION_EQUATION_PROBLEM_SETUP",err,error)
    RETURN 1
    
  END SUBROUTINE ADVECTION_EQUATION_PROBLEM_SETUP
  
  !
  !================================================================================================================================
  !

  !>Sets up the advection equations.
  SUBROUTINE ADVECTION_EQUATION_PROBLEM_LINEAR_SETUP(PROBLEM,PROBLEM_SETUP,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP,CONTROL_LOOP_ROOT
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("ADVECTION_EQUATION_PROBLEM_LINEAR_SETUP",err,error,*999)

    NULLIFY(CONTROL_LOOP)
    NULLIFY(SOLVER)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVERS)
    IF(ASSOCIATED(PROBLEM)) THEN
      IF(.NOT.ALLOCATED(PROBLEM%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(PROBLEM%SPECIFICATION,1)<3) THEN
        CALL FlagError("Problem specification must have three entries for an advection problem.",err,error,*999)
      END IF
      IF(PROBLEM%specification(3)==PROBLEM_ADVECTION_SUBTYPE) THEN
        SELECT CASE(PROBLEM_SETUP%SETUP_TYPE)
        CASE(PROBLEM_SETUP_INITIAL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Do nothing
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Do nothing
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for an advection equation."
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
              & " is invalid for an advection equation."
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
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the solvers
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(SOLVERS,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for an advection equation."
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
              & " is invalid for an advection equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
            & " is invalid for an advection equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        localError="The third problem specification of "// &
          & TRIM(NumberToVString(PROBLEM%specification(3),"*",err,error))// &
          & " does not equal an advection equation subtype."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    ENDIF
       
    EXITS("ADVECTION_EQUATION_PROBLEM_LINEAR_SETUP")
    RETURN
999 ERRORSEXITS("ADVECTION_EQUATION_PROBLEM_LINEAR_SETUP",err,error)
    RETURN 1
    
  END SUBROUTINE ADVECTION_EQUATION_PROBLEM_LINEAR_SETUP
  
  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices for an advection equation finite element equations set.
  SUBROUTINE ADVECTION_EQUATION_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS,GEOMETRIC_BASIS
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: ELEMENTS_TOPOLOGY
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatrixType), POINTER :: dampingMatrix,stiffnessMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,INDEPENDENT_FIELD,GEOMETRIC_FIELD,MATERIALS_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: mhs,ms,ng,nhs,ns,xiIdx,coordIdx,MESH_COMPONENT
    REAL(DP) :: JGW,SUM,DXI_DX,DPHIMS_DXI,DPHINS_DXI,PHIMS,PHINS,flow,area,D,Conc,dConc
    LOGICAL :: updateDampingMatrix,updateStiffnessMatrix

    updateDampingMatrix = .FALSE.
    updateStiffnessMatrix = .FALSE.

    ENTERS("ADVECTION_EQUATION_FINITE_ELEMENT_CALCULATE",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
          CALL FlagError("Equations set specification is not allocated.",err,error,*999)
        ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<3) THEN
          CALL FlagError("Equations set specification must have three entries for an advection problem.",err,error,*999)
        END IF
        NULLIFY(vectorEquations)
        CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
        SELECT CASE(EQUATIONS_SET%specification(3))
        CASE(EQUATIONS_SET_ADVECTION_SUBTYPE)
          !Store all these in equations matrices
          DEPENDENT_FIELD=>equations%interpolation%dependentField
          INDEPENDENT_FIELD=>equations%interpolation%independentField
          MATERIALS_FIELD=>equations%interpolation%materialsField
          GEOMETRIC_FIELD=>equations%interpolation%geometricField
          GEOMETRIC_BASIS=>GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%ptr% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%ptr% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
          vectorMatrices=>vectorEquations%vectorMatrices
          dynamicMatrices=>vectorMatrices%dynamicMatrices
          stiffnessMatrix=>dynamicMatrices%matrices(1)%ptr
          dampingMatrix=>dynamicMatrices%matrices(2)%ptr
          vectorMapping=>vectorEquations%vectorMapping
          dynamicMapping=>vectorMapping%dynamicMapping
          FIELD_VARIABLE=>dynamicMapping%equationsMatrixToVarMaps(1)%variable
          stiffnessMatrix%elementMatrix%matrix=0.0_DP
          dampingMatrix%elementMatrix%matrix=0.0_DP
          IF(ASSOCIATED(dampingMatrix)) updateDampingMatrix=dampingMatrix%updateMatrix
          IF(ASSOCIATED(stiffnessMatrix)) updateStiffnessMatrix=stiffnessMatrix%updateMatrix

          SELECT CASE(EQUATIONS_SET%specification(3))
          CASE(EQUATIONS_SET_ADVECTION_SUBTYPE)
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
              & dependentInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
              & independentInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
              & materialsInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
              & geometricInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            !Loop over gauss points
            DO ng=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS
              CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
                & dependentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
                & independentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
                & materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
                & geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,equations%interpolation% &
                & geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)    
   
              Conc=equations%interpolation%dependentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(1,NO_PART_DERIV)
              dConc=equations%interpolation%dependentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(1,FIRST_PART_DERIV)
              flow=equations%interpolation%independentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(1,NO_PART_DERIV)
              area=equations%interpolation%independentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(2,NO_PART_DERIV)
              D=equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(1,NO_PART_DERIV)

              mhs=0         
              MESH_COMPONENT=FIELD_VARIABLE%COMPONENTS(1)%MESH_COMPONENT_NUMBER
              DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT)%ptr% &
                & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
              QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
              JGW=equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%JACOBIAN* &
                & QUADRATURE_SCHEME%GAUSS_WEIGHTS(ng)
              ELEMENTS_TOPOLOGY=>FIELD_VARIABLE%COMPONENTS(1)%DOMAIN%TOPOLOGY%ELEMENTS

              DXI_DX=0.0_DP
              !Calculate dxi_dx in 3D
              DO xiIdx=1,DEPENDENT_BASIS%NUMBER_OF_XI
                DO coordIdx=1,equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE) &
                  & %ptr%NUMBER_OF_X_DIMENSIONS
                  DXI_DX=DXI_DX+(equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)% &
                    & PTR%DXI_DX(xiIdx,coordIdx))**2.0_DP
                END DO !coordIdx
              END DO !xiIdx
              DXI_DX=SQRT(DXI_DX)

              !Loop over element rows
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                PHIMS=QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                DPHIMS_DXI=QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,FIRST_PART_DERIV,ng)
                mhs=mhs+1
                nhs=0
                IF(updateStiffnessMatrix .OR. updateDampingMatrix) THEN
                  !Loop over element columns
                  DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                    PHINS=QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)
                    DPHINS_DXI=QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,FIRST_PART_DERIV,ng)
                    nhs=nhs+1
                        
                    !!!-- D A M P I N G  M A T R I X --!!!
                    IF(updateDampingMatrix) THEN
                      SUM=PHIMS*PHINS
                      dampingMatrix%elementMatrix%matrix(mhs,nhs)=dampingMatrix%elementMatrix%matrix(mhs,nhs)+SUM*JGW
                    ENDIF
                        
                    !!!-- S T I F F N E S S  M A T R I X --!!!
                    IF(updateStiffnessMatrix) THEN
                      SUM=(D*(DPHIMS_DXI*DXI_DX)+(flow/area)*PHIMS)*DPHINS_DXI*DXI_DX
                      stiffnessMatrix%elementMatrix%matrix(mhs,nhs)=stiffnessMatrix%elementMatrix%matrix(mhs,nhs)+SUM*JGW
                    ENDIF

                  ENDDO !ns
                ENDIF

              ENDDO !ms
            ENDDO !ng
            
          CASE DEFAULT
            localError="The third equations set specification of "// &
              & TRIM(NumberToVString(EQUATIONS_SET%specification(3),"*",err,error))// &
              & " is not valid for an advection type of a classical field equations set."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The third equations set specification of "// &
            & TRIM(NumberToVString(EQUATIONS_SET%specification(3),"*",err,error))// &
            & " is not valid for an advection equation type of a classical field equations set class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        CALL FlagError("Equations set equations is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
    
    EXITS("ADVECTION_EQUATION_FINITE_ELEMENT_CALCULATE")
    RETURN
999 ERRORSEXITS("ADVECTION_EQUATION_FINITE_ELEMENT_CALCULATE",err,error)
    RETURN 1
    
  END SUBROUTINE ADVECTION_EQUATION_FINITE_ELEMENT_CALCULATE

  !
  !================================================================================================================================
  !

  !>Sets up the advection problem pre solve.
  SUBROUTINE Advection_PreSolve(solver,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer to the solver which has the pre solver operations
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop
    TYPE(PROBLEM_TYPE), POINTER :: problem
    TYPE(VARYING_STRING) :: localError

    ENTERS("Advection_PreSolve",err,error,*999)
 
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
    
    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    
    IF(.NOT.ALLOCATED(problem%specification)) CALL FlagError("Problem specification is not allocated.",err,error,*999)
    IF(SIZE(problem%specification,1)<3) &
      & CALL FlagError("Problem specification must have three entries for an advection problem.",err,error,*999)
    
    SELECT CASE(problem%specification(3))
    CASE(PROBLEM_ADVECTION_SUBTYPE, &
      & PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
      & PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
      CALL ADVECTION_PRE_SOLVE_UPDATE_BC(solver,err,error,*999)
      !CALL Advection_Couple1D0D(SOLVER,err,error,*999)
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(problem%specification(3),"*",err,error))// &
        & " is not valid for a advection equation type of a classical field problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Advection_PreSolve")
    RETURN
999 ERRORSEXITS("Advection_PreSolve",err,error)
    RETURN 1
    
  END SUBROUTINE Advection_PreSolve
  
  !
  !================================================================================================================================
  !
  
  !>Update the boundary conditions
  SUBROUTINE ADVECTION_PRE_SOLVE_UPDATE_BC(SOLVER,err,error,*)
    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(EquationsType), POINTER :: EQUATIONS
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: localError
    REAL(DP) :: CONC,START_TIME,STOP_TIME,CURRENT_TIME,TIME_INCREMENT,period,delta(300),t(300),c(300),s
    INTEGER(INTG) :: CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER,i,j,n,m

    ENTERS("ADVECTION_PRE_SOLVE_UPDATE_BC",err,error,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      SOLVERS=>SOLVER%SOLVERS
      IF(ASSOCIATED(SOLVERS)) THEN
        CONTROL_LOOP=>SOLVERS%CONTROL_LOOP
        CALL CONTROL_LOOP_TIMES_GET(CONTROL_LOOP,START_TIME,STOP_TIME,CURRENT_TIME,TIME_INCREMENT,CURRENT_LOOP_ITERATION, &
          & OUTPUT_ITERATION_NUMBER,err,error,*999)
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for an advection problem.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%specification(3))
          CASE(PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
             & PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
            SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
            IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
              SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
              IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
                IF(ASSOCIATED(EQUATIONS)) THEN
                  EQUATIONS_SET=>equations%equationsSet
                  DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                  BOUNDARY_CONDITIONS=>SOLVER_EQUATIONS%BOUNDARY_CONDITIONS
                  IF(ASSOCIATED(DEPENDENT_FIELD)) THEN

                    IF(CURRENT_TIME<100)THEN

                      t(1)=0.003228  ; c(1)=0.001513
                      t(2)=0.077482  ; c(2)=0.001513
                      t(3)=0.133979  ; c(3)=0.013616
                      t(4)=0.187248  ; c(4)=0.040847
                      t(5)=0.237288  ; c(5)=0.108926
                      t(6)=0.285714  ; c(6)=0.226929
                      t(7)=0.33414   ; c(7)=0.414523
                      t(8)=0.417272  ; c(8)=0.800303
                      t(9)=0.45117   ; c(9)=0.92587
                      t(10)=0.479419 ; c(10)=0.984871
                      t(11)=0.499596 ; c(11)=0.995461
                      t(12)=0.519774 ; c(12)=0.984871
                      t(13)=0.550444 ; c(13)=0.919818
                      t(14)=0.583535 ; c(14)=0.795764
                      t(15)=0.661017 ; c(15)=0.429652
                      t(16)=0.698951 ; c(16)=0.282905
                      t(17)=0.722357 ; c(17)=0.208775
                      t(18)=0.753834 ; c(18)=0.128593
                      t(19)=0.785311 ; c(19)=0.07413
                      t(20)=0.824052 ; c(20)=0.034796
                      t(21)=0.874899 ; c(21)=0.012103
                      t(22)=0.91364  ; c(22)=0.004539
                      t(23)=0.999193 ; c(23)=0.0

                      !Initialize variables
                      period=100
                      m=1
                      n=23
                      !Compute derivation
                      DO i=1,n-1
                            delta(i)=(c(i+1)-c(i))/(t(i+1)-t(i))
                      END DO
                      delta(n)=delta(n-1)+(delta(n-1)-delta(n-2))/(t(n-1)-t(n-2))*(t(n)-t(n-1))
                      !Find subinterval
                      DO j=1,n-1
                        IF (t(j) <= (CURRENT_TIME/period)) THEN
                         m=j
                        ENDIF
                      END DO
                      !Evaluate interpolant
                      s=(CURRENT_TIME/period)-t(m)
                      CONC=(c(m)+s*delta(m))
                    ELSE 
                      CONC=0.0
                    ENDIF

                    CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                      & 1,CONC,err,error,*999)
                    CALL FIELD_PARAMETER_SET_UPDATE_START(DEPENDENT_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
                    CALL FIELD_PARAMETER_SET_UPDATE_FINISH(DEPENDENT_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
                  ELSE
                    CALL FlagError("Dependent field and/or geometric field is/are not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("Equations are not associated.",err,error,*999)
                END IF
              ELSE
                CALL FlagError("Solver mapping is not associated.",err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Solver equations are not associated.",err,error,*999)
            END IF
          CASE DEFAULT
            localError="Problem subtype "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
              & " is not valid for a advection equation type of a classical field problem class."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Solvers is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Solver is not associated.",err,error,*999)
    ENDIF

    EXITS("ADVECTION_PRE_SOLVE_UPDATE_BC")
    RETURN
999 ERRORSEXITS("ADVECTION_PRE_SOLVE_UPDATE_BC",err,error)
    RETURN 1

  END SUBROUTINE ADVECTION_PRE_SOLVE_UPDATE_BC
  
  !
  !================================================================================================================================
  !
  
  !>Perform post-solve operations for an advection problem
  SUBROUTINE Advection_PostSolve(solver,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer to the solver for post-solve operations
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop
    TYPE(PROBLEM_TYPE), POINTER :: problem
    TYPE(VARYING_STRING) :: localError

    ENTERS("Advection_PostSolve",err,error,*999)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
    
    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    
    IF(.NOT.ALLOCATED(problem%specification)) CALL FlagError("Problem specification is not allocated.",err,error,*999)
    IF(SIZE(problem%specification,1)<3) &
      & CALL FlagError("Problem specification must have three entries for an advection problem.",err,error,*999)
    
    SELECT CASE(problem%specification(3))
    CASE(PROBLEM_ADVECTION_SUBTYPE)
      !Do nothing
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(problem%specification(3),"*",err,error))// &
        & " is not valid for a advection equation type of a classical field problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Advection_PostSolve")
    RETURN
999 ERRORSEXITS("Advection_PostSolve",err,error)
    RETURN 1
    
  END SUBROUTINE Advection_PostSolve
  
  !
  !================================================================================================================================
  !

  !>Update area for boundary nodes.
  SUBROUTINE Advection_Couple1D0D(SOLVER,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    TYPE(FIELD_TYPE), POINTER :: dependentField
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations  
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping 

    ENTERS("Advection_Couple1D0D",err,error,*999)

    IF(ASSOCIATED(solver)) THEN
      solverEquations=>solver%SOLVER_EQUATIONS
      IF(ASSOCIATED(solverEquations)) THEN
        solverMapping=>solverEquations%SOLVER_MAPPING
        IF(ASSOCIATED(solverMapping)) THEN
          equationsSet=>solverMapping%EQUATIONS_SETS(1)%ptr
          IF(ASSOCIATED(equationsSet)) THEN
            dependentField=>equationsSet%DEPENDENT%DEPENDENT_FIELD
          ELSE
            CALL FlagError("Equations set is not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Solver mapping is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver equations is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Solver is not associated.",err,error,*999)
    ENDIF

    EXITS("Advection_Couple1D0D")
    RETURN
999 ERRORSEXITS("Advection_Couple1D0D",err,error)
    RETURN 1

  END SUBROUTINE

  !
  !================================================================================================================================
  !

END MODULE ADVECTION_EQUATION_ROUTINES

