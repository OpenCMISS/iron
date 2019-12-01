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
  USE ControlLoopRoutines
  USE ControlLoopAccessRoutines
  USE CmissMPI
  USE ComputationRoutines
  USE ComputationAccessRoutines
  USE CoordinateSystemRoutines
  USE DistributedMatrixVector
  USE DomainMappings
  USE EquationsRoutines
  USE EquationsAccessRoutines
  USE EquationsMappingRoutines
  USE EquationsMatricesRoutines
  USE EquationsSetConstants
  USE EquationsSetAccessRoutines
  USE FieldRoutines
  USE FieldAccessRoutines
  USE FIELD_IO_ROUTINES
  USE FLUID_MECHANICS_IO_ROUTINES
  USE InputOutput
  USE ISO_VARYING_STRING
  USE Kinds
  USE Maths
  USE MatrixVector
  USE MeshRoutines
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
  
  PUBLIC Stree_FiniteElementCalculate
  
  PUBLIC Stree_PreSolve

CONTAINS

!
!================================================================================================================================
!

  !>Sets/changes the solution method for a Stree equation type of an fluid mechanics equations set class.
  SUBROUTINE Stree_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet
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
          equationsSet%solutionMethod=EQUATIONS_SET_FEM_SOLUTION_METHOD
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
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the specification for
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
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(EquationsSetSetupType), INTENT(INOUT) :: equationsSetSetup
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    TYPE(DecompositionType), POINTER :: geometricDecomposition
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsSetMaterialsType), POINTER :: equationsMaterials
    TYPE(EquationsSetEquationsSetFieldType), POINTER :: equationsEquationsSetField
    TYPE(FieldType), POINTER :: equationsSetField
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
        SELECT CASE(equationsSetSetup%setupType)
        !-----------------------------------------------------------------
        ! I n i t i a l   s e t u p
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(equationsSet%specification(3))
          CASE(EQUATIONS_SET_STREE1D0D_SUBTYPE)
            SELECT CASE(equationsSetSetup%actionType)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              CALL Stree_EquationsSetSolutionMethodSet(equationsSet, &
                & EQUATIONS_SET_FEM_SOLUTION_METHOD,err,error,*999)
              equationsSet%solutionMethod=EQUATIONS_SET_FEM_SOLUTION_METHOD
              equationsEquationsSetField=>equationsSet%equationsSetField
              IF(equationsEquationsSetField%equationsSetFieldAutoCreated) THEN
                !Create the auto created equations set field field for SUPG element metrics
                CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,equationsSet%REGION, &
                  & equationsEquationsSetField%equationsSetFieldField,err,error,*999)
                equationsSetField=>equationsEquationsSetField%equationsSetFieldField
                CALL Field_LabelSet(equationsSetField,"Equations Set Field",err,error,*999)
                CALL Field_TypeSetAndLock(equationsSetField,FIELD_GENERAL_TYPE,&
                  & err,error,*999)
                CALL Field_NumberOfVariablesSet(equationsSetField,1,err,error,*999)
                CALL Field_VariableTypesSetAndLock(equationsSetField,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                CALL Field_VariableLabelSet(equationsSetField,FIELD_U_VARIABLE_TYPE,"Stree",err,error,*999)
                CALL Field_DataTypeSetAndLock(equationsSetField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL Field_NumberOfComponentsSetAndLock(equationsSetField,FIELD_U_VARIABLE_TYPE,1,err,error,*999)
              ENDIF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF(equationsSet%equationsSetField%equationsSetFieldAutoCreated) THEN
                CALL Field_CreateFinish(equationsSet%equationsSetField%equationsSetFieldField,err,error,*999)
                CALL Field_ComponentValuesInitialise(equationsSet%equationsSetField%equationsSetFieldField, &
                 & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1.0_DP,err,error,*999)
              ENDIF
            CASE DEFAULT
              localError="The action type of "//TRIM(NUMBER_TO_VSTRING(equationsSetSetup%actionType, &
                & "*",err,error))// " for a setup type of "//TRIM(NUMBER_TO_VSTRING(equationsSetSetup% &
                & setupType,"*",err,error))// " is not implemented for a Stree equations set."
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
            SELECT CASE(equationsSetSetup%actionType)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              equationsEquationsSetField=>equationsSet%equationsSetField
              equationsSetField=>equationsEquationsSetField%equationsSetFieldField
              IF(equationsEquationsSetField%equationsSetFieldAutoCreated) THEN
                CALL Field_DecompositionGet(equationsSet%GEOMETRY%geometricField,geometricDecomposition,err,error,*999)
                CALL Field_DecompositionSetAndLock(equationsSetField,geometricDecomposition,err,error,*999)
                CALL Field_GeometricFieldSetAndLock(equationsSetField,equationsSet%GEOMETRY%geometricField,err,error,*999)
                CALL Field_ComponentMeshComponentGet(equationsSet%GEOMETRY%geometricField,FIELD_U_VARIABLE_TYPE, &
                  & 1,geometricComponentNumber,err,error,*999)
                CALL Field_ComponentMeshComponentSetAndLock(equationsSetField,FIELD_U_VARIABLE_TYPE, &
                  & 1,geometricComponentNumber,err,error,*999)
                CALL Field_ComponentInterpolationSetAndLock(equationsSetField,FIELD_U_VARIABLE_TYPE, &
                  & 1,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                !Default the field scaling to that of the geometric field
                CALL Field_ScalingTypeGet(equationsSet%GEOMETRY%geometricField,geometricScalingType,err,error,*999)
                CALL Field_ScalingTypeSet(equationsSet%equationsSetField%equationsSetFieldField, &
                  & geometricScalingType,err,error,*999)
              ELSE
                !Do nothing
              ENDIF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              ! do nothing
            CASE DEFAULT
              localError="The action type of "//TRIM(NUMBER_TO_VSTRING(equationsSetSetup%actionType,"*",err,error))// &
                & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(equationsSetSetup%setupType,"*",err,error))// &
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
            SELECT CASE(equationsSetSetup%actionType)
            !Set start action
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              IF(equationsSet%DEPENDENT%dependentFieldAutoCreated) THEN
                !Create the auto created dependent field
                !start field creation with name 'DEPENDENT_FIELD'
                CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,equationsSet%REGION, &
                  & equationsSet%DEPENDENT%dependentField,err,error,*999)
                !start creation of a new field
                CALL Field_TypeSetAndLock(equationsSet%DEPENDENT%dependentField,FIELD_GENERAL_TYPE,err,error,*999)
                !label the field
                CALL Field_LabelSet(equationsSet%DEPENDENT%dependentField,"Dependent Field",err,error,*999)
                !define new created field to be dependent
                CALL Field_DependentTypeSetAndLock(equationsSet%DEPENDENT%dependentField, &
                  & FIELD_DEPENDENT_TYPE,err,error,*999)
                !look for decomposition rule already defined
                CALL Field_DecompositionGet(equationsSet%GEOMETRY%geometricField,geometricDecomposition, &
                  & err,error,*999)
                !apply decomposition rule found on new created field
                CALL Field_DecompositionSetAndLock(equationsSet%DEPENDENT%dependentField, &
                  & geometricDecomposition,err,error,*999)
                !point new field to geometric field
                CALL Field_GeometricFieldSetAndLock(equationsSet%DEPENDENT%dependentField,equationsSet%GEOMETRY% &
                  & geometricField,err,error,*999)
                !set number of variables
                CALL Field_NumberOfVariablesSetAndLock(equationsSet%DEPENDENT%dependentField, &
                  & dependentFieldNumberOfVariables,err,error,*999)
                CALL Field_VariableTypesSetAndLock(equationsSet%DEPENDENT%dependentField,[FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DELUDELN_VARIABLE_TYPE],err,error,*999)
                !set dimension
                CALL Field_DimensionSetAndLock(equationsSet%DEPENDENT%dependentField,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL Field_DimensionSetAndLock(equationsSet%DEPENDENT%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                !set data type
                CALL Field_DataTypeSetAndLock(equationsSet%DEPENDENT%dependentField,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL Field_DataTypeSetAndLock(equationsSet%DEPENDENT%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                ! number of components for 1
                CALL Field_NumberOfComponentsSetAndLock(equationsSet%DEPENDENT%dependentField, &
                  & FIELD_U_VARIABLE_TYPE,dependentFieldNumberOfComponents,err,error,*999)
                CALL Field_NumberOfComponentsSetAndLock(equationsSet%DEPENDENT%dependentField, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,dependentFieldNumberOfComponents,err,error,*999)
                CALL Field_ComponentMeshComponentGet(equationsSet%GEOMETRY%geometricField,FIELD_U_VARIABLE_TYPE, &
                  & 1,geometricMeshComponent,err,error,*999)
                !Default to the geometric interpolation setup for U,dUdN
                DO I=1,dependentFieldNumberOfComponents
                  CALL Field_ComponentMeshComponentSet(equationsSet%DEPENDENT%dependentField, &
                    & FIELD_U_VARIABLE_TYPE,I,geometricMeshComponent,err,error,*999)
                  CALL Field_ComponentMeshComponentSet(equationsSet%DEPENDENT%dependentField, &
                    & FIELD_DELUDELN_VARIABLE_TYPE,I,geometricMeshComponent,err,error,*999)
                ENDDO
                SELECT CASE(equationsSet%solutionMethod)
                !Specify fem solution method
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  DO I=1,dependentFieldNumberOfComponents
                    CALL Field_ComponentInterpolationSetAndLock(equationsSet%DEPENDENT%dependentField, &
                      & FIELD_U_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL Field_ComponentInterpolationSetAndLock(equationsSet%DEPENDENT%dependentField, &
                      & FIELD_DELUDELN_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  ENDDO
                  CALL Field_ScalingTypeGet(equationsSet%GEOMETRY%geometricField,geometricScalingType,err,error,*999)
                  CALL Field_ScalingTypeSet(equationsSet%DEPENDENT%dependentField,geometricScalingType,err,error,*999)
                CASE DEFAULT
                  localError="The solution method of " &
                    & //TRIM(NUMBER_TO_VSTRING(equationsSet%solutionMethod,"*",err,error))// " is invalid."
                  CALL FLAG_ERROR(localError,err,error,*999)
                END SELECT
              ELSE
                !Check the user specified field
                CALL Field_TypeCheck(equationsSetSetup%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                CALL Field_DependentTypeCheck(equationsSetSetup%FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
                CALL Field_NumberOfVariablesCheck(equationsSetSetup%FIELD,dependentFieldNumberOfVariables,err,error,*999)
                CALL Field_VariableTypesCheck(equationsSetSetup%FIELD,[FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DELUDELN_VARIABLE_TYPE],err,error,*999)
                CALL Field_DimensionCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL Field_DimensionCheck(equationsSetSetup%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL Field_DataTypeCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL Field_DataTypeCheck(equationsSetSetup%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL Field_NumberOfComponentsCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE, &
                  & dependentFieldNumberOfComponents,err,error,*999)
                CALL Field_NumberOfComponentsCheck(equationsSetSetup%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & dependentFieldNumberOfComponents,err,error,*999)
                SELECT CASE(equationsSet%solutionMethod)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  CALL Field_ComponentInterpolationCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationCheck(equationsSetSetup%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CASE DEFAULT
                  localError="The solution method of "//TRIM(NUMBER_TO_VSTRING(equationsSet%solutionMethod, &
                    & "*",err,error))//" is invalid."
                  CALL FLAG_ERROR(localError,err,error,*999)
                END SELECT
              ENDIF
            !Specify finish action
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF(equationsSet%DEPENDENT%dependentFieldAutoCreated) THEN
                CALL Field_CreateFinish(equationsSet%DEPENDENT%dependentField,err,error,*999)
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
            SELECT CASE(equationsSetSetup%actionType)
            !Specify start action
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              equationsMaterials=>equationsSet%MATERIALS
              IF(ASSOCIATED(equationsMaterials)) THEN
                IF(equationsMaterials%materialsFieldAutoCreated) THEN
                  !Create the auto created materials field
                  !start field creation with name 'MATERIAL_FIELD'
                  CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,equationsSet%REGION, &
                    & equationsSet%MATERIALS%materialsField,err,error,*999)
                  CALL Field_TypeSetAndLock(equationsMaterials%materialsField,FIELD_MATERIAL_TYPE,err,error,*999)
                  !label the field
                  CALL Field_LabelSet(equationsMaterials%materialsField,"Materials Field",err,error,*999)
                  CALL Field_DependentTypeSetAndLock(equationsMaterials%materialsField,FIELD_INDEPENDENT_TYPE, &
                    & err,error,*999)
                  CALL Field_DecompositionGet(equationsSet%GEOMETRY%geometricField,geometricDecomposition, &
                    & err,error,*999)
                  !apply decomposition rule found on new created field
                  CALL Field_DecompositionSetAndLock(equationsSet%MATERIALS%materialsField, &
                    & geometricDecomposition,err,error,*999)
                  !point new field to geometric field
                  CALL Field_GeometricFieldSetAndLock(equationsMaterials%materialsField,equationsSet%GEOMETRY% &
                    & geometricField,err,error,*999)
                  CALL Field_NumberOfVariablesSet(equationsMaterials%materialsField, &
                    & materialsFieldNumberOfVariables,err,error,*999)
                  CALL Field_VariableTypesSetAndLock(equationsMaterials%materialsField, &
                    & [FIELD_U_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE],err,error,*999)
                  CALL Field_DimensionSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                  CALL Field_DimensionSetAndLock(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                  CALL Field_DataTypeSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_DP_TYPE,err,error,*999)
                  CALL Field_DataTypeSetAndLock(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE, &
                    & FIELD_DP_TYPE,err,error,*999)
                  CALL Field_NumberOfComponentsSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                    & materialsFieldNumberOfComponents,err,error,*999)
                  CALL Field_NumberOfComponentsSetAndLock(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE, &
                    & materialsFieldNumberOfComponents,err,error,*999)
                  CALL Field_ComponentMeshComponentGet(equationsSet%GEOMETRY%geometricField,FIELD_U_VARIABLE_TYPE, &
                    & 1,geometricComponentNumber,err,error,*999)
                  DO I=1,materialsFieldNumberOfComponents
                    CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                      & I,geometricComponentNumber,err,error,*999)
                    CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE, &
                      & I,geometricComponentNumber,err,error,*999)
                  ENDDO
                  DO I=1,materialsFieldNumberOfComponents
                    CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                      & I,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE, &
                      & I,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  ENDDO
                  !Default the field scaling to that of the geometric field
                  CALL Field_ScalingTypeGet(equationsSet%GEOMETRY%geometricField,geometricScalingType,err,error,*999)
                  CALL Field_ScalingTypeSet(equationsMaterials%materialsField,geometricScalingType,err,error,*999)
                ELSE
                  !Check the user specified field
                  CALL Field_TypeCheck(equationsSetSetup%FIELD,FIELD_MATERIAL_TYPE,err,error,*999)
                  CALL Field_DependentTypeCheck(equationsSetSetup%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                  CALL Field_NumberOfVariablesCheck(equationsSetSetup%FIELD,materialsFieldNumberOfVariables,err,error,*999)
                  CALL Field_VariableTypesCheck(equationsSetSetup%FIELD,[FIELD_U_VARIABLE_TYPE,FIELD_U_VARIABLE_TYPE], &
                    & err,error,*999)
                  CALL Field_DimensionCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                    & err,error,*999)
                  CALL Field_DimensionCheck(equationsSetSetup%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                    & err,error,*999)
                  CALL Field_DataTypeCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                  CALL Field_DataTypeCheck(equationsSetSetup%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                  CALL Field_NumberOfComponentsCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE, &
                    & materialsFieldNumberOfComponents,err,error,*999)
                  CALL Field_NumberOfComponentsCheck(equationsSetSetup%FIELD,FIELD_V_VARIABLE_TYPE, &
                    & materialsFieldNumberOfComponents,err,error,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Equations set materials is not associated.",err,error,*999)
              END IF
              !Specify start action
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              equationsMaterials=>equationsSet%MATERIALS
              IF(ASSOCIATED(equationsMaterials)) THEN
                IF(equationsMaterials%materialsFieldAutoCreated) THEN
                  CALL Field_CreateFinish(equationsMaterials%materialsField,err,error,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Equations set materials is not associated.",err,error,*999)
              ENDIF
            CASE DEFAULT
              localError="The action type of "//TRIM(NUMBER_TO_VSTRING(equationsSetSetup%actionType,"*", &
                & err,error))//" for a setup type of "//TRIM(NUMBER_TO_VSTRING(equationsSetSetup%setupType,"*", &
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
            SELECT CASE(equationsSetSetup%actionType)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              CALL EquationsSet_AssertMaterialsIsFinished(equationsSet,err,error,*999)
              CALL Equations_CreateStart(equationsSet,equations,err,error,*999)
              CALL Equations_LinearityTypeSet(equations,EQUATIONS_LINEAR,err,error,*999)
              CALL Equations_TimeDependenceTypeSet(equations,EQUATIONS_STATIC,err,error,*999)
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              SELECT CASE(equationsSet%solutionMethod)
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
                localError="The solution method of "//TRIM(NUMBER_TO_VSTRING(equationsSet%solutionMethod, &
                  & "*",err,error))//" is invalid."
                CALL FLAG_ERROR(localError,err,error,*999)
              END SELECT
            CASE DEFAULT
              localError="The action type of "//TRIM(NUMBER_TO_VSTRING(equationsSetSetup%actionType,"*",err,error))// &
                & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(equationsSetSetup%setupType,"*",err,error))// &
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
          localError="The setup type of "//TRIM(NUMBER_TO_VSTRING(equationsSetSetup%setupType,"*",err,error))// &
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
  SUBROUTINE Stree_FiniteElementCalculate(equationsSet,nodeNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet
    INTEGER(INTG), INTENT(IN) :: nodeNumber
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(DomainType), POINTER :: domain
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatrixType), POINTER :: stiffnessMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: materialsField,dependentField
    TYPE(FieldVariableType), POINTER :: fieldVariable
    TYPE(VARYING_STRING) :: localError
    REAL(DP), POINTER :: dependentParameters(:),materialsParameters(:)

    ENTERS("Stree_FiniteElementCalculate",err,error,*999)

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
        dependentField=>equations%equationsSet%DEPENDENT%dependentField
        IF(ASSOCIATED(dependentField)) THEN
          domain=>dependentField%DECOMPOSITION%DOMAIN(dependentField%decomposition%meshComponentNumber)%ptr
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

    EXITS("Stree_FiniteElementCalculate")
    RETURN
999 ERRORSEXITS("Stree_FiniteElementCalculate",err,error)
    RETURN 1

  END SUBROUTINE Stree_FiniteElementCalculate

  !
  !================================================================================================================================
  !

  !>Evaluates the residual nodal stiffness matrices and RHS for a Stree equation nodal equations set.
  SUBROUTINE Stree_PreSolve(solver,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(ControlLoopType), POINTER :: controlLoop,parentLoop,navierstokesLoop
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(DomainType), POINTER :: domain
    TYPE(EquationsSetType), POINTER :: equationsSet,navierstokesEquationsSet
    TYPE(EquationsType), POINTER :: equations,navierstokesEquations
    TYPE(FieldType), POINTER :: materialsField,navierstokesDependentField
    TYPE(FieldVariableType), POINTER :: dependentFieldVariable
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations,navierstokesSolverEquations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping,navierstokesSolverMapping
    TYPE(SOLVERS_TYPE), POINTER :: solvers
    TYPE(SOLVER_TYPE), POINTER :: navierstokesSolver
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: variableIdx,BOUNDARY_CONDITION_CHECK_VARIABLE,nodeIdx,componentIdx,derivativeIdx,versionIdx
    INTEGER(INTG) :: dependentDof,dependentVariableType,userNodeNumber,m
    REAL(DP) :: currentTime,timeIncrement,flow

    ENTERS("Stree_PreSolve",err,error,*999)

    ! Some preliminary sanity checks
    IF(ASSOCIATED(SOLVER)) THEN
      solvers=>SOLVER%SOLVERS
      IF(ASSOCIATED(SOLVERS)) THEN
        controlLoop=>solvers%controlLoop
        CALL CONTROL_LOOP_CURRENT_TIMES_GET(controlLoop,currentTime,timeIncrement,err,error,*999)
        parentLoop=>controlLoop%parentLoop
        navierstokesLoop=>parentLoop%subLoops(2)%ptr
        navierstokesSolver=>navierstokesLoop%SOLVERS%SOLVERS(2)%ptr
        IF(ASSOCIATED(controlLoop%PROBLEM)) THEN
          SELECT CASE(controlLoop%PROBLEM%specification(3))
          CASE(PROBLEM_STREE1D0D_NAVIER_STOKES_SUBTYPE, &
             & PROBLEM_STREE1D0D_ADV_NAVIER_STOKES_SUBTYPE)
            solverEquations=>solver%SOLVER_EQUATIONS
            navierstokesSolverEquations=>navierstokesSolver%SOLVER_EQUATIONS
            IF(ASSOCIATED(solverEquations)) THEN
              solverMapping=>solverEquations%solverMapping
              navierstokesSolverMapping=>navierstokesSolverEquations%solverMapping
              IF(ASSOCIATED(solverMapping)) THEN
                equationsSet=>solverMapping%equationsSets(1)%ptr
                navierstokesEquationsSet=>navierstokesSolverMapping%equationsSets(1)%ptr
                IF(ASSOCIATED(equationsSet)) THEN
                  equations=>equationsSet%EQUATIONS
                  navierstokesEquations=>navierstokesEquationsSet%EQUATIONS
                  IF(ASSOCIATED(equations)) THEN
                    materialsField=>equationsSet%MATERIALS%materialsField
                    navierstokesDependentField=>navierstokesEquationsSet%DEPENDENT%dependentField
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
    DO variableIdx=1,navierstokesDependentField%numberOfVariables
      dependentVariableType=navierstokesDependentField%VARIABLES(variableIdx)%variableType
      NULLIFY(dependentFieldVariable)
      CALL Field_VariableGet(navierstokesDependentField,dependentVariableType,dependentFieldVariable,err,error,*999)
      CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS, &
        & dependentFieldVariable,BOUNDARY_CONDITIONS_VARIABLE,err,error,*999)
      IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
        IF(ASSOCIATED(dependentFieldVariable)) THEN
          DO componentIdx=1,dependentFieldVariable%numberOfComponents
            IF(dependentFieldVariable%COMPONENTS(componentIdx)%interpolationType==FIELD_NODE_BASED_INTERPOLATION) THEN
              domain=>dependentFieldVariable%COMPONENTS(componentIdx)%DOMAIN
              IF(ASSOCIATED(domain)) THEN
                IF(ASSOCIATED(domain%TOPOLOGY)) THEN
                  domainNodes=>domain%TOPOLOGY%NODES
                  IF(ASSOCIATED(domainNodes)) THEN
                    !Loop over the local nodes excluding the ghosts.
                    DO nodeIdx=1,domainNodes%numberOfNodes
                      userNodeNumber=domainNodes%nodes(nodeIdx)%userNumber
                      DO derivativeIdx=1,domainNodes%nodes(nodeIdx)%numberOfDerivatives
                        DO versionIdx=1,domainNodes%nodes(nodeIdx)%DERIVATIVES(derivativeIdx)%numberOfVersions
                          CALL FieldVariable_LocalNodeDOFGet(dependentFieldVariable,versionIdx,derivativeIdx, &
                            & nodeIdx,componentIdx,dependentDOF,err,error,*999)
                          BOUNDARY_CONDITION_CHECK_VARIABLE=BOUNDARY_CONDITIONS_VARIABLE%CONDITION_TYPES(dependentDof)
                          IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_FIXED_STREE) THEN
                            ! Update dependent field value
                            IF(ASSOCIATED(materialsField)) THEN
                              CALL Field_ParameterSetGetNode(navierstokesDependentField,FIELD_U_VARIABLE_TYPE, &
                                & FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx,nodeIdx,1,flow,err,error,*999)
                              m=int(currentTime)-800*(int(currentTime)/800)
                              CALL Field_ParameterSetUpdateLocalNode(materialsField,FIELD_V_VARIABLE_TYPE, &
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

    EXITS("Stree_PreSolve")
    RETURN
999 ERRORSEXITS("Stree_PreSolve",err,error)
    RETURN 1

  END SUBROUTINE Stree_PreSolve

  !
  !================================================================================================================================
  !

END MODULE Stree_EQUATION_ROUTINES
