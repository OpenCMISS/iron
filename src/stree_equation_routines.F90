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
!> Contributor(s): Soroush Safaei, Chris Bradley
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
MODULE StreeEquationsRoutines

  USE BaseRoutines
  USE BasisRoutines
  USE BoundaryConditionsRoutines
  USE BoundaryConditionAccessRoutines
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
  USE EquationsMatricesAccessRoutines
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
  USE ProblemAccessRoutines
  USE Strings
  USE SolverRoutines
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
    INTEGER(INTG) :: esSpecification(3)
    TYPE(VARYING_STRING) :: localError

    ENTERS("Stree_EquationsSetSolutionMethodSet",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_STREE1D0D_SUBTYPE)
      SELECT CASE(solutionMethod)
      CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
        equationsSet%solutionMethod=EQUATIONS_SET_FEM_SOLUTION_METHOD
      CASE(EQUATIONS_SET_NODAL_SOLUTION_METHOD)
        CALL FlagError("Not implemented.",err,error,*999)
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
        localError="The specified solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The third equations set specification of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is not valid for a Stree type of a fluid mechanics equations set."
      CALL FlagError(localError,err,error,*999)
    END SELECT

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

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(SIZE(specification,1)<3) THEN
      localError="The size of the specified specification array of "// &
        & TRIM(NumberToVString(SIZE(specification,1),"*",err,error))//" is invalid. The size should be >= 3." 
      CALL FlagError(localError,err,error,*999)
    ENDIF
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
    IF(ALLOCATED(equationsSet%specification)) CALL FlagError("Equations set specification is already allocated.",err,error,*999)  
    ALLOCATE(equationsSet%specification(3),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate equations set specification.",err,error,*999)
    equationsSet%specification(1:3)=[EQUATIONS_SET_FLUID_MECHANICS_CLASS,EQUATIONS_SET_STREE_EQUATION_TYPE,subtype]

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
    INTEGER(INTG) :: componentIdx,dependentFieldNumberOfComponents,dependentFieldNumberOfVariables,esSpecification(3), &
      & geometricScalingType,geometricMeshComponent,geometricComponentNumber,materialsFieldNumberOfVariables, &
      & materialsFieldNumberOfComponents,solutionMethod
    TYPE(DecompositionType), POINTER :: geometricDecomposition
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsSetMaterialsType), POINTER :: equationsMaterials
    TYPE(EquationsSetEquationsFieldType), POINTER :: equationsField
    TYPE(FieldType), POINTER :: equationsSetField,geometricField
    TYPE(RegionType), POINTER :: region
    TYPE(VARYING_STRING) :: localError

    ENTERS("Stree_EquationsSetSetup",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)

    SELECT CASE(equationsSet%specification(3))
    CASE(EQUATIONS_SET_STREE1D0D_SUBTYPE)
      !OK
    CASE DEFAULT
      localError="The third equations set specification of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    NULLIFY(region)
    CALL EquationsSet_RegionGet(equationsSet,region,err,error,*999)
    
    SELECT CASE(equationsSetSetup%setupType)
    CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
      !-----------------------------------------------------------------
      ! I n i t i a l   s e t u p
      !-----------------------------------------------------------------
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        CALL Stree_EquationsSetSolutionMethodSet(equationsSet,EQUATIONS_SET_FEM_SOLUTION_METHOD,err,error,*999)
        NULLIFY(equationsField)
        CALL EquationsSet_EquationsFieldGet(equationsSet,equationsField,err,error,*999)
        IF(equationsField%equationsSetFieldAutoCreated) THEN
          !Create the auto created equations set field field for SUPG element metrics
          CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsField%equationsSetField,err,error,*999)
          NULLIFY(equationsSetField)
          CALL EquationsSet_EquationsSetFieldGet(equationsSet,equationsSetField,err,error,*999)
          CALL Field_LabelSet(equationsSetField,"Equations Set Field",err,error,*999)
          CALL Field_TypeSetAndLock(equationsSetField,FIELD_GENERAL_TYPE,err,error,*999)
          CALL Field_NumberOfVariablesSet(equationsSetField,1,err,error,*999)
          CALL Field_VariableTypesSetAndLock(equationsSetField,[FIELD_U_VARIABLE_TYPE],err,error,*999)
          CALL Field_VariableLabelSet(equationsSetField,FIELD_U_VARIABLE_TYPE,"Stree",err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsSetField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_NumberOfComponentsSetAndLock(equationsSetField,FIELD_U_VARIABLE_TYPE,1,err,error,*999)
        ELSE
          !Check the user specified field
          CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
          CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,1,err,error,*999)
          CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE],err,error,*999)
          CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,1,err,error,*999)
        ENDIF
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        IF(equationsField%equationsSetFieldAutoCreated) THEN
          CALL Field_CreateFinish(equationsField%equationsSetField,err,error,*999)
          CALL Field_ComponentValuesInitialise(equationsField%equationsSetField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & 1,1.0_DP,err,error,*999)
        ENDIF
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType, &
          & "*",err,error))// " for a setup type of "//TRIM(NumberToVString(equationsSetSetup% &
          & setupType,"*",err,error))// " is not implemented for a Stree equations set."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
      !-----------------------------------------------------------------
      ! G e o m e t r i c   f i e l d
      !-----------------------------------------------------------------
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        NULLIFY(equationsField)
        CALL EquationsSet_EquationsFieldGet(equationsSet,equationsField,err,error,*999)
        NULLIFY(equationsSetField)
        CALL EquationsSet_EquationsSetFieldGet(equationsSet,equationsSetField,err,error,*999)
        IF(equationsField%equationsSetFieldAutoCreated) THEN
          NULLIFY(geometricDecomposition)
          CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
          CALL Field_DecompositionSetAndLock(equationsSetField,geometricDecomposition,err,error,*999)
          CALL Field_GeometricFieldSetAndLock(equationsSetField,geometricField,err,error,*999)
          CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricComponentNumber,err,error,*999)
          CALL Field_ComponentMeshComponentSetAndLock(equationsSetField,FIELD_U_VARIABLE_TYPE,1,geometricComponentNumber, &
            & err,error,*999)
          CALL Field_ComponentInterpolationSetAndLock(equationsSetField,FIELD_U_VARIABLE_TYPE,1,FIELD_CONSTANT_INTERPOLATION, &
            & err,error,*999)
          !Default the field scaling to that of the geometric field
          CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
          CALL Field_ScalingTypeSet(equationsField%equationsSetField,geometricScalingType,err,error,*999)
        ELSE
          !Do nothing
        ENDIF
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        !Do nothing
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a Stree equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
      !-----------------------------------------------------------------
      ! D e p e n d e n t   f i e l d
      !-----------------------------------------------------------------
      dependentFieldNumberOfVariables=2
      dependentFieldNumberOfComponents=1
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        NULLIFY(geometricField)
        CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
        IF(equationsSet%dependent%dependentFieldAutoCreated) THEN
          !Create the auto created dependent field
          !start field creation with name 'Dependent Field'
          CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsSet%dependent%dependentField,err,error,*999)
          !start creation of a new field
          CALL Field_TypeSetAndLock(equationsSet%dependent%dependentField,FIELD_GENERAL_TYPE,err,error,*999)
          !label the field
          CALL Field_LabelSet(equationsSet%dependent%dependentField,"Dependent Field",err,error,*999)
          !define new created field to be dependent
          CALL Field_DependentTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_DEPENDENT_TYPE,err,error,*999)
          !look for decomposition rule already defined
          CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
          !apply decomposition rule found on new created field
          CALL Field_DecompositionSetAndLock(equationsSet%dependent%dependentField,geometricDecomposition,err,error,*999)
          !point new field to geometric field
          CALL Field_GeometricFieldSetAndLock(equationsSet%dependent%dependentField,geometricField,err,error,*999)
          !set number of variables
          CALL Field_NumberOfVariablesSetAndLock(equationsSet%dependent%dependentField,dependentFieldNumberOfVariables, &
            & err,error,*999)
          CALL Field_VariableTypesSetAndLock(equationsSet%dependent%dependentField,[FIELD_U_VARIABLE_TYPE, &
            & FIELD_DELUDELN_VARIABLE_TYPE],err,error,*999)
          !set dimension
          CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
            & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
            & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
          !set data type
          CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE, &
            & err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE, &
            & err,error,*999)
          ! number of components for 1
          CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
            & dependentFieldNumberOfComponents,err,error,*999)
          CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
            & dependentFieldNumberOfComponents,err,error,*999)
          CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricMeshComponent,err,error,*999)
          !Default to the geometric interpolation setup for U,dUdN
          DO componentIdx=1,dependentFieldNumberOfComponents
            CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,componentIdx, &
              & geometricMeshComponent,err,error,*999)
            CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,componentIdx, &
              & geometricMeshComponent,err,error,*999)
          ENDDO !componentIdx
          CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
          SELECT CASE(solutionMethod)
            !Specify fem solution method
          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
            DO componentIdx=1,dependentFieldNumberOfComponents
              CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
                & componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
                & componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
            ENDDO !componentIdx
            CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
            CALL Field_ScalingTypeSet(equationsSet%dependent%dependentField,geometricScalingType,err,error,*999)
          CASE DEFAULT
            localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))// " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE
          !Check the user specified field
          CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
          CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_DEPENDENT_TYPE,err,error,*999)
          CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,dependentFieldNumberOfVariables,err,error,*999)
          CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE],err,error,*999)
          CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,dependentFieldNumberOfComponents, &
            & err,error,*999)
          CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE, &
            & dependentFieldNumberOfComponents,err,error,*999)
          CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
          SELECT CASE(solutionMethod)
          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
            CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION, &
              & err,error,*999)
            CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,1, &
              & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
          CASE DEFAULT
            localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ENDIF
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        IF(equationsSet%dependent%dependentFieldAutoCreated) THEN
          CALL Field_CreateFinish(equationsSet%dependent%dependentField,err,error,*999)
        ENDIF
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a Stree equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
      !-----------------------------------------------------------------
      ! M a t e r i a l s   f i e l d
      !-----------------------------------------------------------------
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      materialsFieldNumberOfVariables=2
      materialsFieldNumberOfComponents=27
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        NULLIFY(equationsMaterials)
        CALL EquationsSet_EquationsMaterialsGet(equationsSet,equationsMaterials,err,error,*999)
        IF(equationsMaterials%materialsFieldAutoCreated) THEN
          !Create the auto created materials field
          !start field creation with name 'Materials Field'
          CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsSet%materials%materialsField,err,error,*999)
          CALL Field_TypeSetAndLock(equationsMaterials%materialsField,FIELD_MATERIAL_TYPE,err,error,*999)
          !label the field
          CALL Field_LabelSet(equationsMaterials%materialsField,"Materials Field",err,error,*999)
          CALL Field_DependentTypeSetAndLock(equationsMaterials%materialsField,FIELD_INDEPENDENT_TYPE,err,error,*999)
          NULLIFY(geometricDecomposition)
          CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
          !apply decomposition rule found on new created field
          CALL Field_DecompositionSetAndLock(equationsMaterials%materialsField,geometricDecomposition,err,error,*999)
          !point new field to geometric field
          CALL Field_GeometricFieldSetAndLock(equationsMaterials%materialsField,geometricField,err,error,*999)
          CALL Field_NumberOfVariablesSet(equationsMaterials%materialsField,materialsFieldNumberOfVariables,err,error,*999)
          CALL Field_VariableTypesSetAndLock(equationsMaterials%materialsField, &
            & [FIELD_U_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE],err,error,*999)
          CALL Field_DimensionSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
            & err,error,*999)
          CALL Field_DimensionSetAndLock(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
            & err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_NumberOfComponentsSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
            & materialsFieldNumberOfComponents,err,error,*999)
          CALL Field_NumberOfComponentsSetAndLock(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE, &
            & materialsFieldNumberOfComponents,err,error,*999)
          CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricComponentNumber,err,error,*999)
          DO componentIdx=1,materialsFieldNumberOfComponents
            CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,componentIdx, &
              & geometricComponentNumber,err,error,*999)
            CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE,componentIdx, &
              & geometricComponentNumber,err,error,*999)
          ENDDO !componentIdx
          DO componentIdx=1,materialsFieldNumberOfComponents
            CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,componentIdx, &
              & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
            CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE,componentIdx, &
              & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
          ENDDO !componentIdx
          !Default the field scaling to that of the geometric field
          CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
          CALL Field_ScalingTypeSet(equationsMaterials%materialsField,geometricScalingType,err,error,*999)
        ELSE
          !Check the user specified field
          CALL Field_TypeCheck(equationsSetSetup%field,FIELD_MATERIAL_TYPE,err,error,*999)
          CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
          CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,materialsFieldNumberOfVariables,err,error,*999)
          CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE,FIELD_U_VARIABLE_TYPE], &
            & err,error,*999)
          CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
            & err,error,*999)
          CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
            & err,error,*999)
          CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE, &
            & materialsFieldNumberOfComponents,err,error,*999)
          CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE, &
            & materialsFieldNumberOfComponents,err,error,*999)
        ENDIF
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        NULLIFY(equationsMaterials)
        CALL EquationsSet_EquationsMaterialsGet(equationsSet,equationsMaterials,err,error,*999)
        IF(equationsMaterials%materialsFieldAutoCreated) THEN
          CALL Field_CreateFinish(equationsMaterials%materialsField,err,error,*999)
!!TODO: INITIALISE MATERIALS FIELD WITH DEFAULT VALUES
        ENDIF
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*", &
          & err,error))//" for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*", &
          & err,error))//" is invalid for Stree equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
      !-----------------------------------------------------------------
      ! E q u a t i o n s    t y p e
      !-----------------------------------------------------------------
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        CALL EquationsSet_AssertDependentIsFinished(equationsSet,err,error,*999)
        CALL EquationsSet_AssertMaterialsIsFinished(equationsSet,err,error,*999)
        NULLIFY(equations)
        CALL Equations_CreateStart(equationsSet,equations,err,error,*999)
        CALL Equations_LinearityTypeSet(equations,EQUATIONS_LINEAR,err,error,*999)
        CALL Equations_TimeDependenceTypeSet(equations,EQUATIONS_STATIC,err,error,*999)
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
        SELECT CASE(solutionMethod)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          !Finish the creation of the equations
          NULLIFY(equations)
          CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
          CALL Equations_CreateFinish(equations,err,error,*999)
          NULLIFY(vectorEquations)
          CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
          !Create the equations mapping.
          NULLIFY(vectorMapping)
          CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_U_VARIABLE_TYPE,vectorMapping,err,error,*999)
          CALL EquationsMappingVector_NumberOfLinearMatricesSet(vectorMapping,1,err,error,*999)
          CALL EquationsMappingVector_LinearMatricesVariableTypesSet(vectorMapping,[FIELD_U_VARIABLE_TYPE],err,error,*999)
          CALL EquationsMappingVector_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE,err,error,*999)
          CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
          !Create the equations matrices
          NULLIFY(vectorMatrices)
          CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
          SELECT CASE(equations%sparsityType)
          CASE(EQUATIONS_MATRICES_FULL_MATRICES)
            CALL EquationsMatricesVector_LinearStorageTypeSet(vectorMatrices,[MATRIX_BLOCK_STORAGE_TYPE],err,error,*999)
          CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
            CALL EquationsMatricesVector_LinearStorageTypeSet(vectorMatrices,[MATRIX_COMPRESSED_ROW_STORAGE_TYPE],err,error,*999)
            CALL EquationsMatricesVector_LinearStructureTypeSet(vectorMatrices,[EQUATIONS_MATRIX_FEM_STRUCTURE],err,error,*999)
          CASE DEFAULT
            localError="The equations matrices sparsity type of "//TRIM(NumberToVString(equations%sparsityType,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
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
          localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a Strees equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
        & " is invalid for a Strees equation set."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
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
    INTEGER(INTG) :: esSpecification(3)
    LOGICAL :: update,updateMatrix,updateRHS
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatrixType), POINTER :: linearMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldVariableType), POINTER :: rowsVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("Stree_FiniteElementCalculate",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)

    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_STREE1D0D_SUBTYPE)
      !OK
    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is not valid for a Stree equation type of a fluid mechanics equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(equationsInterpolation)
    CALL Equations_InterpolationGet(equations,equationsInterpolation,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMapping)
    CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
    NULLIFY(lhsMapping)
    CALL EquationsMappingVector_LHSMappingGet(vectorMapping,lhsMapping,err,error,*999)
    NULLIFY(rowsVariable)
    CALL EquationsMappingLHS_LHSVariableGet(lhsMapping,rowsVariable,err,error,*999)
    NULLIFY(linearMapping)
    CALL EquationsMappingVector_LinearMappingGet(vectorMapping,linearMapping,err,error,*999)
    NULLIFY(rhsMapping)
    CALL EquationsMappingVector_RHSMappingExists(vectorMapping,rhsMapping,err,error,*999)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
    NULLIFY(linearMatrices)
    CALL EquationsMatricesVector_LinearMatricesGet(vectorMatrices,linearMatrices,err,error,*999)
    NULLIFY(linearMatrix)
    CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,1,linearMatrix,err,error,*999)
    CALL EquationsMatrix_UpdateMatrixGet(linearMatrix,updateMatrix,err,error,*999)
    NULLIFY(rhsVector)
    updateRHS=.FALSE.
    IF(ASSOCIATED(rhsMapping)) THEN
      CALL EquationsMatricesVector_RHSVectorGet(vectorMatrices,rhsVector,err,error,*999)
      CALL EquationsMatricesRHS_UpdateVectorGet(rhsVector,updateRHS,err,error,*999)
    ENDIF

    update=(updateMatrix.OR.updateRHS)

    IF(update) THEN
    ENDIF

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
    TYPE(SolverType), POINTER :: solver
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    INTEGER(INTG) :: boundaryConditionCheckVariable,componentIdx,currentIteration,derivativeIdx,dependentDof, &
      & dependentVariableType,inputIteration,m,nodeIdx,numberOfComponents,numberOfNodes,numberOfNodeDerivatives, &
      & numberOfVariables,numberOfVersions,outputIteration,pSpecification(3),userNodeNumber,variableIdx,versionIdx
    REAL(DP) :: currentTime,flow,startTime,stopTime,timeIncrement
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable
    TYPE(ControlLoopType), POINTER :: controlLoop,parentLoop,navierStokesLoop
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(EquationsSetType), POINTER :: equationsSet,navierStokesEquationsSet
    TYPE(EquationsType), POINTER :: equations,navierStokesEquations
    TYPE(FieldType), POINTER :: materialsField,navierStokesDependentField
    TYPE(FieldVariableType), POINTER :: dependentFieldVariable
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverEquationsType), POINTER :: solverEquations,navierStokesSolverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping,navierStokesSolverMapping
    TYPE(SolversType), POINTER :: navierStokesSolvers,solvers
    TYPE(SolverType), POINTER :: navierStokesSolver
    TYPE(VARYING_STRING) :: localError

    ENTERS("Stree_PreSolve",err,error,*999)

    NULLIFY(solvers)
    CALL Solver_SolversGet(solver,solvers,err,error,*999)
    NULLIFY(controlLoop)
    CALL Solvers_ControlLoopGet(solvers,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    CALL ControlLoop_CurrentTimeInformationGet(controlLoop,currentTime,timeIncrement,startTime,stopTime,currentIteration, &
      & outputIteration,inputIteration,err,error,*999)
    NULLIFY(parentLoop)
    CALL ControlLoop_ParentLoopGet(controlLoop,parentLoop,err,error,*999)
    NULLIFY(navierStokesLoop)
    CALL ControlLoop_SubLoopGet(parentLoop,2,navierStokesLoop,err,error,*999)
    NULLIFY(navierStokesSolvers)
    CALL ControlLoop_SolversGet(navierStokesLoop,navierStokesSolvers,err,error,*999)
    NULLIFY(navierStokesSolver)
    CALL Solvers_SolverGet(navierStokesSolvers,2,navierStokesSolver,err,error,*999)
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STREE1D0D_NAVIER_STOKES_SUBTYPE, &
      & PROBLEM_STREE1D0D_ADV_NAVIER_STOKES_SUBTYPE)
      NULLIFY(solverEquations)
      CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
      NULLIFY(navierStokesSolverEquations)
      CALL Solver_SolverEquationsGet(navierStokesSolver,navierStokesSolverEquations,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      NULLIFY(navierStokesSolverMapping)
      CALL SolverEquations_SolverMappingGet(navierStokesSolverMapping,navierStokesSolverMapping,err,error,*999)
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
      NULLIFY(navierStokesEquationsSet)
      CALL SolverMapping_EquationsSetGet(navierStokesSolverMapping,navierStokesEquationsSet,err,error,*999)
      NULLIFY(materialsField)
      CALL EquationsSet_MaterialsFieldExists(equationsSet,materialsField,err,error,*999)
      NULLIFY(navierStokesDependentField)
      CALL EquationsSet_DependentFieldGet(navierStokesEquationsSet,navierStokesDependentField,err,error,*999)
    CASE DEFAULT
      localError="The third problem specification of "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for boundary flux calculation."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    NULLIFY(boundaryConditions)
    CALL SolverEquations_BoundaryConditionsGet(navierStokesSolverEquations,boundaryConditions,err,error,*999)
    CALL Field_NumberOfVariablesGet(navierStokesDependentField,numberOfVariables,err,error,*999)
    DO variableIdx=1,numberOfVariables
      NULLIFY(dependentFieldVariable)
      CALL Field_VariableIndexGet(navierStokesDependentField,variableIdx,dependentFieldVariable,dependentVariableType, &
        & err,error,*999)
      NULLIFY(boundaryConditionsVariable)
      CALL BoundaryConditions_VariableGet(boundaryConditions,dependentFieldVariable,boundaryConditionsVariable,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(dependentFieldVariable,numberOfComponents,err,error,*999)
      DO componentIdx=1,numberOfComponents
        CALL FieldVariable_ComponentInterpolationCheck(dependentFieldVariable,componentIdx,FIELD_NODE_BASED_INTERPOLATION, &
          & err,error,*999)
        NULLIFY(domain)
        CALL FieldVariable_ComponentDomainGet(dependentFieldVariable,componentIdx,domain,err,error,*999)
        NULLIFY(domainTopology)
        CALL Domain_DomainTopology(domain,domainTopology,err,error,*999)
        NULLIFY(domainNodes)
        CALL DomainTopology_DomainNodes(domainTopology,domainNodes,err,error,*999)
        !Loop over the local nodes excluding the ghosts.
        CALL DomainNodes_NumberOfNodesGet(domainNodes,numberOfNodes,err,error,*999)
        DO nodeIdx=1,numberOfNodes
          CALL DomainNodes_NodeUserNumberGet(domainNodes,nodeIdx,userNodeNumber,err,error,*999)
          CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,nodeIdx,numberOfNodeDerivatives,err,error,*999)
          DO derivativeIdx=1,numberOfNodeDerivatives
            CALL DomainNodes_DerivativeNumberOfVersionsGet(domainNodes,derivativeIdx,nodeIdx,numberOfVersions,err,error,*999)
            DO versionIdx=1,numberOfVersions
              CALL FieldVariable_LocalNodeDOFGet(dependentFieldVariable,versionIdx,derivativeIdx,nodeIdx,componentIdx, &
                & dependentDOF,err,error,*999)
              boundaryConditionCheckVariable=boundaryConditionsVariable%conditionTypes(dependentDof)
              IF(boundaryConditionCheckVariable==BOUNDARY_CONDITION_FIXED_STREE) THEN
                ! Update dependent field value
                IF(ASSOCIATED(materialsField)) THEN
                  CALL Field_ParameterSetGetNode(navierStokesDependentField,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx,nodeIdx,1,flow,err,error,*999)
                  m=INT(currentTime)-800*(INT(currentTime)/800)
                  CALL Field_ParameterSetUpdateLocalNode(materialsField,FIELD_V_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx,m+1,1,flow,err,error,*999)
                ENDIF
              ENDIF
            ENDDO !versionIdx
          ENDDO !derivativeIdx
        ENDDO !nodeIdx
      ENDDO !componentIdx
    ENDDO !variableIdx

    EXITS("Stree_PreSolve")
    RETURN
999 ERRORSEXITS("Stree_PreSolve",err,error)
    RETURN 1

  END SUBROUTINE Stree_PreSolve

  !
  !================================================================================================================================
  !

END MODULE StreeEquationsRoutines
