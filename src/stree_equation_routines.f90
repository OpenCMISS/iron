!> \file
!> \author Soroush Safaei
!> \brief This module handles the Stree equation routines. These
!>  equations are often used in concert with 1D fluid modelling to describe
!>  wave propagation phenomena, which is particularly useful for models of
!>  vascular trees. These equations are also often solved using a discontinuous
!>  nodal solution method, rather than FEM.
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

!>This module handles all Stree equation routines.
MODULE Stree_EQUATION_ROUTINES

  USE BaseRoutines
  USE BasisRoutines
  USE BOUNDARY_CONDITIONS_ROUTINES
  USE Constants
  USE CONTROL_LOOP_ROUTINES
  USE CmissMPI
  USE ComputationEnvironment
  USE COORDINATE_ROUTINES
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
  USE FIELD_IO_ROUTINES
  USE FLUID_MECHANICS_IO_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE Kinds
  USE Maths
  USE MatrixVector
  USE MESH_ROUTINES
  USE NODE_ROUTINES
  USE PROBLEM_CONSTANTS
  USE Strings
  USE SOLVER_ROUTINES
  USE Timer
  USE Types

#include "macros.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC Stree_EquationsSetSolutionMethodSet
  PUBLIC Stree_EquationsSetSpecificationSet
  PUBLIC Stree_EquationsSetSetup
  PUBLIC Stree_FINITE_ELEMENT_CALCULATE
  PUBLIC Stree_PRE_SOLVE

CONTAINS

!
!================================================================================================================================
!

  !>Sets/changes the solution method for a Stree equation type of an fluid mechanics equations set class.
  SUBROUTINE Stree_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    INTEGER(INTG), INTENT(IN) :: solutionMethod
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Stree_EquationsSetSolutionMethodSet",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(.NOT.ALLOCATED(equationsSet%specification)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(equationsSet%specification,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a structured tree type equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(equationsSet%specification(3))
      CASE(EQUATIONS_SET_STREE1D0D_SUBTYPE)
        SELECT CASE(solutionMethod)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          equationsSet%SOLUTION_METHOD=EQUATIONS_SET_FEM_SOLUTION_METHOD
        CASE(EQUATIONS_SET_NODAL_SOLUTION_METHOD)
          CALL FLAG_ERROR("Not implemented.",err,error,*999)
        CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
          CALL FLAG_ERROR("Not implemented.",err,error,*999)
        CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
          CALL FLAG_ERROR("Not implemented.",err,error,*999)
        CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
          CALL FLAG_ERROR("Not implemented.",err,error,*999)
        CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
          CALL FLAG_ERROR("Not implemented.",err,error,*999)
        CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
          CALL FLAG_ERROR("Not implemented.",err,error,*999)
        CASE DEFAULT
          localError="The specified solution method of "//TRIM(NUMBER_TO_VSTRING(solutionMethod,"*",err,error))// &
            & " is invalid."
          CALL FLAG_ERROR(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The third equations set specification of "// &
          & TRIM(NUMBER_TO_VSTRING(equationsSet%specification(3),"*",err,error))// &
          & " is not valid for a Stree type of a fluid mechanics equations set."
        CALL FLAG_ERROR(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",err,error,*999)
    ENDIF

    EXITS("Stree_EquationsSetSolutionMethodSet")
    RETURN
999 ERRORSEXITS("Stree_EquationsSetSolutionMethodSet",err,error)
    RETURN 1

  END SUBROUTINE Stree_EquationsSetSolutionMethodSet

!
!================================================================================================================================
!

  !>Sets the equation specification for a Stree type of a fluid mechanics equations set.
  SUBROUTINE Stree_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to set the specification for
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations set specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: subtype
    TYPE(VARYING_STRING) :: localError

    ENTERS("Stree_EquationsSetSpecificationSet",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(SIZE(specification,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Stree type equations set.", &
          & err,error,*999)
      END IF
      subtype=specification(3)
      SELECT CASE(subtype)
      CASE(EQUATIONS_SET_STREE1D0D_SUBTYPE)
        !ok
      CASE DEFAULT
        localError="The third equations set specification of "//TRIM(NumberToVstring(subtype,"*",err,error))// &
          & " is not valid for a Stree type of a fluid mechanics equations set."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Set full specification
      IF(ALLOCATED(equationsSet%specification)) THEN
        CALL FlagError("Equations set specification is already allocated.",err,error,*999)
      ELSE
        ALLOCATE(equationsSet%specification(3),stat=err)
        IF(err/=0) CALL FlagError("Could not allocate equations set specification.",err,error,*999)
      END IF
      equationsSet%specification(1:3)=[EQUATIONS_SET_FLUID_MECHANICS_CLASS,EQUATIONS_SET_STREE_EQUATION_TYPE,subtype]
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    END IF

    EXITS("Stree_EquationsSetSpecificationSet")
    RETURN
999 ERRORSEXITS("Stree_EquationsSetSpecificationSet",err,error)
    RETURN 1

  END SUBROUTINE Stree_EquationsSetSpecificationSet

!
!================================================================================================================================
!

  !>Sets up the Stree equations fluid setup.
  SUBROUTINE Stree_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: equationsSetSetup
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: geometricDecomposition
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EQUATIONS_SET_MATERIALS_TYPE), POINTER :: equationsMaterials
    TYPE(EQUATIONS_SET_EQUATIONS_SET_FIELD_TYPE), POINTER :: equationsEquationsSetField
    TYPE(FIELD_TYPE), POINTER :: equationsSetField
    INTEGER(INTG) :: I,geometricScalingType,geometricMeshComponent,geometricComponentNumber
    INTEGER(INTG) :: dependentFieldNumberOfVariables,dependentFieldNumberOfComponents
    INTEGER(INTG) :: materialsFieldNumberOfVariables,materialsFieldNumberOfComponents
    TYPE(VARYING_STRING) :: localError

    ENTERS("Stree_EquationsSetSetup",err,error,*999)

    NULLIFY(equations)
    NULLIFY(vectorMapping)
    NULLIFY(vectorMatrices)
    NULLIFY(equationsMaterials)
    NULLIFY(geometricDecomposition)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(.NOT.ALLOCATED(equationsSet%specification)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(equationsSet%specification,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Stree type equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(equationsSet%specification(3))
      CASE(EQUATIONS_SET_STREE1D0D_SUBTYPE)
        SELECT CASE(equationsSetSetup%SETUP_TYPE)
        !-----------------------------------------------------------------
        ! I n i t i a l   s e t u p
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(equationsSet%specification(3))
          CASE(EQUATIONS_SET_STREE1D0D_SUBTYPE)
            SELECT CASE(equationsSetSetup%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              CALL Stree_EquationsSetSolutionMethodSet(equationsSet, &
                & EQUATIONS_SET_FEM_SOLUTION_METHOD,err,error,*999)
              equationsSet%SOLUTION_METHOD=EQUATIONS_SET_FEM_SOLUTION_METHOD
              equationsEquationsSetField=>equationsSet%EQUATIONS_SET_FIELD
              IF(equationsEquationsSetField%EQUATIONS_SET_FIELD_AUTO_CREATED) THEN
                !Create the auto created equations set field field for SUPG element metrics
                CALL FIELD_CREATE_START(equationsSetSetup%FIELD_USER_NUMBER,equationsSet%REGION, &
                  & equationsEquationsSetField%EQUATIONS_SET_FIELD_FIELD,err,error,*999)
                equationsSetField=>equationsEquationsSetField%EQUATIONS_SET_FIELD_FIELD
                CALL FIELD_LABEL_SET(equationsSetField,"Equations Set Field",err,error,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(equationsSetField,FIELD_GENERAL_TYPE,&
                  & err,error,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_SET(equationsSetField,1,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(equationsSetField,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                CALL FIELD_VARIABLE_LABEL_SET(equationsSetField,FIELD_U_VARIABLE_TYPE,"Stree",err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(equationsSetField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(equationsSetField,FIELD_U_VARIABLE_TYPE,1,err,error,*999)
              ENDIF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF(equationsSet%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_AUTO_CREATED) THEN
                CALL FIELD_CREATE_FINISH(equationsSet%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,err,error,*999)
                CALL FIELD_COMPONENT_VALUES_INITIALISE(equationsSet%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                 & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1.0_DP,err,error,*999)
              ENDIF
            CASE DEFAULT
              localError="The action type of "//TRIM(NUMBER_TO_VSTRING(equationsSetSetup%ACTION_TYPE, &
                & "*",err,error))// " for a setup type of "//TRIM(NUMBER_TO_VSTRING(equationsSetSetup% &
                & SETUP_TYPE,"*",err,error))// " is not implemented for a Stree equations set."
              CALL FLAG_ERROR(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The third equations set specification of "// &
              & TRIM(NUMBER_TO_VSTRING(equationsSet%specification(3),"*",err,error))// &
              & " is invalid for a Stree equations set."
            CALL FLAG_ERROR(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! G e o m e t r i c   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
          SELECT CASE(equationsSet%specification(3))
          CASE(EQUATIONS_SET_STREE1D0D_SUBTYPE)
            SELECT CASE(equationsSetSetup%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              equationsEquationsSetField=>equationsSet%EQUATIONS_SET_FIELD
              equationsSetField=>equationsEquationsSetField%EQUATIONS_SET_FIELD_FIELD
              IF(equationsEquationsSetField%EQUATIONS_SET_FIELD_AUTO_CREATED) THEN
                CALL FIELD_MESH_DECOMPOSITION_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,geometricDecomposition,err,error,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(equationsSetField,geometricDecomposition,err,error,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(equationsSetField,equationsSet%GEOMETRY%GEOMETRIC_FIELD,err,error,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,geometricComponentNumber,err,error,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET_AND_LOCK(equationsSetField,FIELD_U_VARIABLE_TYPE, &
                  & 1,geometricComponentNumber,err,error,*999)
                CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(equationsSetField,FIELD_U_VARIABLE_TYPE, &
                  & 1,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                !Default the field scaling to that of the geometric field
                CALL FIELD_SCALING_TYPE_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,geometricScalingType,err,error,*999)
                CALL FIELD_SCALING_TYPE_SET(equationsSet%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                  & geometricScalingType,err,error,*999)
              ELSE
                !Do nothing
              ENDIF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              ! do nothing
            CASE DEFAULT
              localError="The action type of "//TRIM(NUMBER_TO_VSTRING(equationsSetSetup%ACTION_TYPE,"*",err,error))// &
                & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a Stree equation."
              CALL FLAG_ERROR(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The third equations set specification of "// &
              & TRIM(NUMBER_TO_VSTRING(equationsSet%specification(3),"*",err,error))// &
              & " is invalid for a Stree equations set."
            CALL FLAG_ERROR(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! D e p e n d e n t   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
          SELECT CASE(equationsSet%specification(3))
          CASE(EQUATIONS_SET_STREE1D0D_SUBTYPE)
            dependentFieldNumberOfVariables=2
            dependentFieldNumberOfComponents=1
            SELECT CASE(equationsSetSetup%ACTION_TYPE)
            !Set start action
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              IF(equationsSet%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                !Create the auto created dependent field
                !start field creation with name 'DEPENDENT_FIELD'
                CALL FIELD_CREATE_START(equationsSetSetup%FIELD_USER_NUMBER,equationsSet%REGION, &
                  & equationsSet%DEPENDENT%DEPENDENT_FIELD,err,error,*999)
                !start creation of a new field
                CALL FIELD_TYPE_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                !label the field
                CALL FIELD_LABEL_SET(equationsSet%DEPENDENT%DEPENDENT_FIELD,"Dependent Field",err,error,*999)
                !define new created field to be dependent
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_DEPENDENT_TYPE,err,error,*999)
                !look for decomposition rule already defined
                CALL FIELD_MESH_DECOMPOSITION_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,geometricDecomposition, &
                  & err,error,*999)
                !apply decomposition rule found on new created field
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                  & geometricDecomposition,err,error,*999)
                !point new field to geometric field
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD,equationsSet%GEOMETRY% &
                  & GEOMETRIC_FIELD,err,error,*999)
                !set number of variables
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                  & dependentFieldNumberOfVariables,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD,[FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DELUDELN_VARIABLE_TYPE],err,error,*999)
                !set dimension
                CALL FIELD_DIMENSION_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                !set data type
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                ! number of components for 1
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_U_VARIABLE_TYPE,dependentFieldNumberOfComponents,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,dependentFieldNumberOfComponents,err,error,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,geometricMeshComponent,err,error,*999)
                !Default to the geometric interpolation setup for U,dUdN
                DO I=1,dependentFieldNumberOfComponents
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,I,geometricMeshComponent,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_DELUDELN_VARIABLE_TYPE,I,geometricMeshComponent,err,error,*999)
                ENDDO
                SELECT CASE(equationsSet%SOLUTION_METHOD)
                !Specify fem solution method
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  DO I=1,dependentFieldNumberOfComponents
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_DELUDELN_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  ENDDO
                  CALL FIELD_SCALING_TYPE_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,geometricScalingType,err,error,*999)
                  CALL FIELD_SCALING_TYPE_SET(equationsSet%DEPENDENT%DEPENDENT_FIELD,geometricScalingType,err,error,*999)
                CASE DEFAULT
                  localError="The solution method of " &
                    & //TRIM(NUMBER_TO_VSTRING(equationsSet%SOLUTION_METHOD,"*",err,error))// " is invalid."
                  CALL FLAG_ERROR(localError,err,error,*999)
                END SELECT
              ELSE
                !Check the user specified field
                CALL FIELD_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(equationsSetSetup%FIELD,dependentFieldNumberOfVariables,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(equationsSetSetup%FIELD,[FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DELUDELN_VARIABLE_TYPE],err,error,*999)
                CALL FIELD_DIMENSION_CHECK(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_CHECK(equationsSetSetup%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE, &
                  & dependentFieldNumberOfComponents,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(equationsSetSetup%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & dependentFieldNumberOfComponents,err,error,*999)
                SELECT CASE(equationsSet%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(equationsSetSetup%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CASE DEFAULT
                  localError="The solution method of "//TRIM(NUMBER_TO_VSTRING(equationsSet%SOLUTION_METHOD, &
                    & "*",err,error))//" is invalid."
                  CALL FLAG_ERROR(localError,err,error,*999)
                END SELECT
              ENDIF
            !Specify finish action
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF(equationsSet%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                CALL FIELD_CREATE_FINISH(equationsSet%DEPENDENT%DEPENDENT_FIELD,err,error,*999)
              ENDIF
            CASE DEFAULT
              localError="The third equations set specification of "// &
                & TRIM(NUMBER_TO_VSTRING(equationsSet%specification(3),"*",err,error))// &
                & " is invalid for a Stree equations set."
              CALL FLAG_ERROR(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The third equations set specification of "// &
              & TRIM(NUMBER_TO_VSTRING(equationsSet%specification(3),"*",err,error))// &
              & " is invalid for a Stree equations set."
            CALL FLAG_ERROR(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! M a t e r i a l s   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
          SELECT CASE(equationsSet%specification(3))
          CASE(EQUATIONS_SET_STREE1D0D_SUBTYPE)
            materialsFieldNumberOfVariables=2
            materialsFieldNumberOfComponents=27
            SELECT CASE(equationsSetSetup%ACTION_TYPE)
            !Specify start action
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              equationsMaterials=>equationsSet%MATERIALS
              IF(ASSOCIATED(equationsMaterials)) THEN
                IF(equationsMaterials%MATERIALS_FIELD_AUTO_CREATED) THEN
                  !Create the auto created materials field
                  !start field creation with name 'MATERIAL_FIELD'
                  CALL FIELD_CREATE_START(equationsSetSetup%FIELD_USER_NUMBER,equationsSet%REGION, &
                    & equationsSet%MATERIALS%MATERIALS_FIELD,err,error,*999)
                  CALL FIELD_TYPE_SET_AND_LOCK(equationsMaterials%MATERIALS_FIELD,FIELD_MATERIAL_TYPE,err,error,*999)
                  !label the field
                  CALL FIELD_LABEL_SET(equationsMaterials%MATERIALS_FIELD,"Materials Field",err,error,*999)
                  CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(equationsMaterials%MATERIALS_FIELD,FIELD_INDEPENDENT_TYPE, &
                    & err,error,*999)
                  CALL FIELD_MESH_DECOMPOSITION_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,geometricDecomposition, &
                    & err,error,*999)
                  !apply decomposition rule found on new created field
                  CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(equationsSet%MATERIALS%MATERIALS_FIELD, &
                    & geometricDecomposition,err,error,*999)
                  !point new field to geometric field
                  CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(equationsMaterials%MATERIALS_FIELD,equationsSet%GEOMETRY% &
                    & GEOMETRIC_FIELD,err,error,*999)
                  CALL FIELD_NUMBER_OF_VARIABLES_SET(equationsMaterials%MATERIALS_FIELD, &
                    & materialsFieldNumberOfVariables,err,error,*999)
                  CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(equationsMaterials%MATERIALS_FIELD, &
                    & [FIELD_U_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE],err,error,*999)
                  CALL FIELD_DIMENSION_SET_AND_LOCK(equationsMaterials%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                  CALL FIELD_DIMENSION_SET_AND_LOCK(equationsMaterials%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                  CALL FIELD_DATA_TYPE_SET_AND_LOCK(equationsMaterials%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_DATA_TYPE_SET_AND_LOCK(equationsMaterials%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(equationsMaterials%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & materialsFieldNumberOfComponents,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(equationsMaterials%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & materialsFieldNumberOfComponents,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 1,geometricComponentNumber,err,error,*999)
                  DO I=1,materialsFieldNumberOfComponents
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(equationsMaterials%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & I,geometricComponentNumber,err,error,*999)
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(equationsMaterials%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                      & I,geometricComponentNumber,err,error,*999)
                  ENDDO
                  DO I=1,materialsFieldNumberOfComponents
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(equationsMaterials%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & I,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(equationsMaterials%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                      & I,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  ENDDO
                  !Default the field scaling to that of the geometric field
                  CALL FIELD_SCALING_TYPE_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,geometricScalingType,err,error,*999)
                  CALL FIELD_SCALING_TYPE_SET(equationsMaterials%MATERIALS_FIELD,geometricScalingType,err,error,*999)
                ELSE
                  !Check the user specified field
                  CALL FIELD_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_MATERIAL_TYPE,err,error,*999)
                  CALL FIELD_DEPENDENT_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                  CALL FIELD_NUMBER_OF_VARIABLES_CHECK(equationsSetSetup%FIELD,materialsFieldNumberOfVariables,err,error,*999)
                  CALL FIELD_VARIABLE_TYPES_CHECK(equationsSetSetup%FIELD,[FIELD_U_VARIABLE_TYPE,FIELD_U_VARIABLE_TYPE], &
                    & err,error,*999)
                  CALL FIELD_DIMENSION_CHECK(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                    & err,error,*999)
                  CALL FIELD_DIMENSION_CHECK(equationsSetSetup%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                    & err,error,*999)
                  CALL FIELD_DATA_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_DATA_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE, &
                    & materialsFieldNumberOfComponents,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(equationsSetSetup%FIELD,FIELD_V_VARIABLE_TYPE, &
                    & materialsFieldNumberOfComponents,err,error,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Equations set materials is not associated.",err,error,*999)
              END IF
              !Specify start action
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              equationsMaterials=>equationsSet%MATERIALS
              IF(ASSOCIATED(equationsMaterials)) THEN
                IF(equationsMaterials%MATERIALS_FIELD_AUTO_CREATED) THEN
                  CALL FIELD_CREATE_FINISH(equationsMaterials%MATERIALS_FIELD,err,error,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Equations set materials is not associated.",err,error,*999)
              ENDIF
            CASE DEFAULT
              localError="The action type of "//TRIM(NUMBER_TO_VSTRING(equationsSetSetup%ACTION_TYPE,"*", &
                & err,error))//" for a setup type of "//TRIM(NUMBER_TO_VSTRING(equationsSetSetup%SETUP_TYPE,"*", &
                & err,error))//" is invalid for Stree equation."
              CALL FLAG_ERROR(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The third equations set specification of "// &
              & TRIM(NUMBER_TO_VSTRING(equationsSet%specification(3),"*",err,error))// &
              & " is invalid for a Stree equation."
            CALL FLAG_ERROR(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! E q u a t i o n s    t y p e
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
          SELECT CASE(equationsSet%specification(3))
          CASE(EQUATIONS_SET_STREE1D0D_SUBTYPE)
            SELECT CASE(equationsSetSetup%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              equationsMaterials=>equationsSet%MATERIALS
              IF(ASSOCIATED(equationsMaterials)) THEN
                IF(equationsMaterials%MATERIALS_FINISHED) THEN
                  CALL Equations_CreateStart(equationsSet,equations,err,error,*999)
                  CALL Equations_LinearityTypeSet(equations,EQUATIONS_LINEAR,err,error,*999)
                  CALL Equations_TimeDependenceTypeSet(equations,EQUATIONS_STATIC,err,error,*999)
                ELSE
                  CALL FLAG_ERROR("Equations set materials has not been finished.",err,error,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Equations materials is not associated.",err,error,*999)
              ENDIF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              SELECT CASE(equationsSet%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                !Finish the creation of the equations
                CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
                CALL Equations_CreateFinish(equations,err,error,*999)
                NULLIFY(vectorEquations)
                CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
                !Create the equations mapping.
                CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_DELUDELN_VARIABLE_TYPE,vectorMapping,err,error,*999)
                CALL EquationsMapping_LinearMatricesNumberSet(vectorMapping,1,err,error,*999)
                CALL EquationsMapping_LinearMatricesVariableTypesSet(vectorMapping,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                CALL EquationsMapping_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE,err,error,*999)
                CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
                !Create the equations matrices
                CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
                SELECT CASE(equations%sparsityType)
                CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                  CALL EquationsMatrices_LinearStorageTypeSet(vectorMatrices,[MATRIX_BLOCK_STORAGE_TYPE],err,error,*999)
                CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                  CALL EquationsMatrices_LinearStorageTypeSet(vectorMatrices,[MATRIX_COMPRESSED_ROW_STORAGE_TYPE],err,error,*999)
                  CALL EquationsMatrices_LinearStructureTypeSet(vectorMatrices,[EQUATIONS_MATRIX_FEM_STRUCTURE],err,error,*999)
                CASE DEFAULT
                  localError="The equations matrices sparsity type of "// &
                    & TRIM(NUMBER_TO_VSTRING(equations%sparsityType,"*",err,error))//" is invalid."
                  CALL FLAG_ERROR(localError,err,error,*999)
                END SELECT
                CALL EquationsMatrices_VectorCreateFinish(vectorMatrices,err,error,*999)
              CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",err,error,*999)
              CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",err,error,*999)
              CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",err,error,*999)
              CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",err,error,*999)
              CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",err,error,*999)
              CASE DEFAULT
                localError="The solution method of "//TRIM(NUMBER_TO_VSTRING(equationsSet%SOLUTION_METHOD, &
                  & "*",err,error))//" is invalid."
                CALL FLAG_ERROR(localError,err,error,*999)
              END SELECT
            CASE DEFAULT
              localError="The action type of "//TRIM(NUMBER_TO_VSTRING(equationsSetSetup%ACTION_TYPE,"*",err,error))// &
                & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a Strees equation."
              CALL FLAG_ERROR(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The third equations set specification of "// &
              & TRIM(NUMBER_TO_VSTRING(equationsSet%specification(3),"*",err,error))// &
              & " is invalid for a Stree equation."
            CALL FLAG_ERROR(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The setup type of "//TRIM(NUMBER_TO_VSTRING(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
            & " is invalid for a Strees equation set."
          CALL FLAG_ERROR(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The third equations set specification of "// &
          & TRIM(NUMBER_TO_VSTRING(equationsSet%specification(3),"*",err,error))// &
          & " does not equal a Strees equation set."
        CALL FLAG_ERROR(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",err,error,*999)
    ENDIF

    EXITS("Stree_EquationsSetSetup")
    RETURN
999 ERRORSEXITS("Stree_EquationsSetSetup",err,error)
    RETURN 1

  END SUBROUTINE Stree_EquationsSetSetup

  !
  !================================================================================================================================
  !

  !>Evaluates the residual nodal stiffness matrices and RHS for a Stree equation nodal equations set.
  SUBROUTINE STREE_FINITE_ELEMENT_CALCULATE(equationsSet,nodeNumber,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    INTEGER(INTG), INTENT(IN) :: nodeNumber
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    TYPE(DOMAIN_NODES_TYPE), POINTER :: domainNodes
    TYPE(DOMAIN_TYPE), POINTER :: domain
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatrixType), POINTER :: stiffnessMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FIELD_TYPE), POINTER :: materialsField,dependentField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    TYPE(VARYING_STRING) :: localError
    REAL(DP), POINTER :: dependentParameters(:),materialsParameters(:)

    ENTERS("STREE_FINITE_ELEMENT_CALCULATE",err,error,*999)

    NULLIFY(equations)
    NULLIFY(vectorMapping)
    NULLIFY(vectorMatrices)
    NULLIFY(linearMapping)
    NULLIFY(linearMatrices)
    NULLIFY(stiffnessMatrix)
    NULLIFY(dependentField)
    NULLIFY(materialsField)
    NULLIFY(domain)
    NULLIFY(domainNodes)
    NULLIFY(dependentParameters)
    NULLIFY(materialsParameters)
    NULLIFY(fieldVariable)

    IF(ASSOCIATED(equationsSet)) THEN
      equations=>equationsSet%EQUATIONS
      IF(ASSOCIATED(equations)) THEN
        NULLIFY(vectorEquations)
        CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
        dependentField=>equations%equationsSet%DEPENDENT%DEPENDENT_FIELD
        IF(ASSOCIATED(dependentField)) THEN
          domain=>dependentField%DECOMPOSITION%DOMAIN(dependentField%DECOMPOSITION%MESH_COMPONENT_NUMBER)%ptr
          IF(ASSOCIATED(domain)) THEN
            domainNodes=>domain%TOPOLOGY%NODES
          ELSE
            CALL FLAG_ERROR("Domain is not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Dependent Field is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations set equations is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",err,error,*999)
    ENDIF

    SELECT CASE(equationsSet%specification(3))
    CASE(EQUATIONS_SET_STREE1D0D_SUBTYPE)
      !Set General and Specific Pointers
      vectorMatrices=>vectorEquations%vectorMatrices
      vectorMapping=>vectorEquations%vectorMapping
      linearMatrices=>vectorMatrices%linearMatrices
      stiffnessMatrix=>linearMatrices%matrices(1)%ptr
      linearMapping=>vectorMapping%linearMapping
      stiffnessMatrix%elementMatrix%matrix=0.0_DP

    CASE DEFAULT
      localError="The third equations set specification of "// &
        & TRIM(NUMBER_TO_VSTRING(equationsSet%specification(3),"*",err,error))// &
        & " is not valid for a Stree type of a fluid mechanics equations set."
      CALL FLAG_ERROR(localError,err,error,*999)
    END SELECT

    EXITS("STREE_FINITE_ELEMENT_CALCULATE")
    RETURN
999 ERRORSEXITS("STREE_FINITE_ELEMENT_CALCULATE",err,error)
    RETURN 1

  END SUBROUTINE STREE_FINITE_ELEMENT_CALCULATE

  !
  !================================================================================================================================
  !

  !>Evaluates the residual nodal stiffness matrices and RHS for a Stree equation nodal equations set.
  SUBROUTINE Stree_PRE_SOLVE(solver,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop,parentLoop,navierstokesLoop
    TYPE(DOMAIN_NODES_TYPE), POINTER :: domainNodes
    TYPE(DOMAIN_TYPE), POINTER :: domain
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet,navierstokesEquationsSet
    TYPE(EquationsType), POINTER :: equations,navierstokesEquations
    TYPE(FIELD_TYPE), POINTER :: materialsField,navierstokesDependentField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: dependentFieldVariable
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations,navierstokesSolverEquations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping,navierstokesSolverMapping
    TYPE(SOLVERS_TYPE), POINTER :: solvers
    TYPE(SOLVER_TYPE), POINTER :: navierstokesSolver
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: variableIdx,BOUNDARY_CONDITION_CHECK_VARIABLE,nodeIdx,componentIdx,derivativeIdx,versionIdx
    INTEGER(INTG) :: dependentDof,dependentVariableType,userNodeNumber,m
    REAL(DP) :: currentTime,timeIncrement,flow

    ENTERS("Stree_PRE_SOLVE",err,error,*999)

    ! Some preliminary sanity checks
    IF(ASSOCIATED(SOLVER)) THEN
      solvers=>SOLVER%SOLVERS
      IF(ASSOCIATED(SOLVERS)) THEN
        controlLoop=>solvers%CONTROL_LOOP
        CALL CONTROL_LOOP_CURRENT_TIMES_GET(controlLoop,currentTime,timeIncrement,err,error,*999)
        parentLoop=>controlLoop%PARENT_LOOP
        navierstokesLoop=>parentLoop%SUB_LOOPS(2)%ptr
        navierstokesSolver=>navierstokesLoop%SOLVERS%SOLVERS(2)%ptr
        IF(ASSOCIATED(controlLoop%PROBLEM)) THEN
          SELECT CASE(controlLoop%PROBLEM%specification(3))
          CASE(PROBLEM_STREE1D0D_NAVIER_STOKES_SUBTYPE, &
             & PROBLEM_STREE1D0D_ADV_NAVIER_STOKES_SUBTYPE)
            solverEquations=>solver%SOLVER_EQUATIONS
            navierstokesSolverEquations=>navierstokesSolver%SOLVER_EQUATIONS
            IF(ASSOCIATED(solverEquations)) THEN
              solverMapping=>solverEquations%SOLVER_MAPPING
              navierstokesSolverMapping=>navierstokesSolverEquations%SOLVER_MAPPING
              IF(ASSOCIATED(solverMapping)) THEN
                equationsSet=>solverMapping%EQUATIONS_SETS(1)%ptr
                navierstokesEquationsSet=>navierstokesSolverMapping%EQUATIONS_SETS(1)%ptr
                IF(ASSOCIATED(equationsSet)) THEN
                  equations=>equationsSet%EQUATIONS
                  navierstokesEquations=>navierstokesEquationsSet%EQUATIONS
                  IF(ASSOCIATED(equations)) THEN
                    materialsField=>equationsSet%MATERIALS%MATERIALS_FIELD
                    navierstokesDependentField=>navierstokesEquationsSet%DEPENDENT%DEPENDENT_FIELD
                    IF(.NOT.ASSOCIATED(materialsField)) THEN
                      CALL FLAG_ERROR("Dependent field is not associated.",err,error,*999)
                    END IF
                  ELSE
                    CALL FLAG_ERROR("Equations set equations is not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FLAG_ERROR("Equations set is not associated.",err,error,*999)
                END IF
              ELSE
                CALL FLAG_ERROR("Solver mapping is not associated.",err,error,*999)
              END IF
            ELSE
              CALL FLAG_ERROR("Solver equations is not associated.",err,error,*999)
            END IF
          CASE DEFAULT
            localError="The third problem specification of "// &
              & TRIM(NUMBER_TO_VSTRING(controlLoop%PROBLEM%specification(3),"*",err,error))// &
              & " is not valid for boundary flux calculation."
            CALL FLAG_ERROR(localError,err,error,*999)
          END SELECT
        ELSE
          CALL FLAG_ERROR("Problem is not associated.",err,error,*999)
        END IF
      ELSE
        CALL FLAG_ERROR("Solvers is not associated.",err,error,*999)
      END IF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",err,error,*999)
    END IF

    BOUNDARY_CONDITIONS=>navierstokesSolverEquations%BOUNDARY_CONDITIONS
    DO variableIdx=1,navierstokesDependentField%NUMBER_OF_VARIABLES
      dependentVariableType=navierstokesDependentField%VARIABLES(variableIdx)%VARIABLE_TYPE
      NULLIFY(dependentFieldVariable)
      CALL Field_VariableGet(navierstokesDependentField,dependentVariableType,dependentFieldVariable,err,error,*999)
      CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS, &
        & dependentFieldVariable,BOUNDARY_CONDITIONS_VARIABLE,err,error,*999)
      IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
        IF(ASSOCIATED(dependentFieldVariable)) THEN
          DO componentIdx=1,dependentFieldVariable%NUMBER_OF_COMPONENTS
            IF(dependentFieldVariable%COMPONENTS(componentIdx)%INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN
              domain=>dependentFieldVariable%COMPONENTS(componentIdx)%DOMAIN
              IF(ASSOCIATED(domain)) THEN
                IF(ASSOCIATED(domain%TOPOLOGY)) THEN
                  domainNodes=>domain%TOPOLOGY%NODES
                  IF(ASSOCIATED(domainNodes)) THEN
                    !Loop over the local nodes excluding the ghosts.
                    DO nodeIdx=1,domainNodes%NUMBER_OF_NODES
                      userNodeNumber=domainNodes%NODES(nodeIdx)%USER_NUMBER
                      DO derivativeIdx=1,domainNodes%NODES(nodeIdx)%NUMBER_OF_DERIVATIVES
                        DO versionIdx=1,domainNodes%NODES(nodeIdx)%DERIVATIVES(derivativeIdx)%numberOfVersions
                          dependentDof = dependentFieldVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP% &
                            & NODE_PARAM2DOF_MAP%NODES(nodeIdx)%DERIVATIVES(derivativeIdx)%VERSIONS(versionIdx)
                          BOUNDARY_CONDITION_CHECK_VARIABLE=BOUNDARY_CONDITIONS_VARIABLE%CONDITION_TYPES(dependentDof)
                          IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_FIXED_STREE) THEN
                            ! Update dependent field value
                            IF(ASSOCIATED(materialsField)) THEN
                              CALL FIELD_PARAMETER_SET_GET_NODE(navierstokesDependentField,FIELD_U_VARIABLE_TYPE, &
                                & FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx,nodeIdx,1,flow,err,error,*999)
                              m=int(currentTime)-800*(int(currentTime)/800)
                              CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(materialsField,FIELD_V_VARIABLE_TYPE, &
                                & FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx,m+1,1,flow,err,error,*999)
                            ENDIF
                          ENDIF
                        ENDDO !versionIdx
                      ENDDO !derivativeIdx
                    ENDDO !nodeIdx
                  ELSE
                    CALL FLAG_ERROR("Domain topology nodes is not associated.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Domain topology is not associated.",err,error,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Domain is not associated.",err,error,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Only node based interpolation is implemented.",err,error,*999)
            ENDIF
          ENDDO !componentIdx
        ELSE
          CALL FLAG_ERROR("Dependent field variable is not associated.",err,error,*999)
        ENDIF
      ENDIF
    ENDDO !variableIdx

    EXITS("Stree_PRE_SOLVE")
    RETURN
999 ERRORSEXITS("Stree_PRE_SOLVE",err,error)
    RETURN 1

  END SUBROUTINE Stree_PRE_SOLVE

  !
  !================================================================================================================================
  !

END MODULE Stree_EQUATION_ROUTINES
