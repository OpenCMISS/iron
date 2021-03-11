!> \file
!> \author Sebastian Krittian
!> \brief This module handles all Stokes fluid routines.
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
!> Contributor(s): Sebastian Krittian, Chris Bradley
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

!>This module handles all Stokes fluid routines.
MODULE StokesEquationsRoutines

  USE AnalyticAnalysisRoutines
  USE BaseRoutines
  USE BasisRoutines
  USE BasisAccessRoutines
  USE BoundaryConditionsRoutines
  USE BoundaryConditionAccessRoutines
  USE Constants
  USE ControlLoopRoutines
  USE ControlLoopAccessRoutines
  USE DistributedMatrixVector
  USE DistributedMatrixVectorAccessRoutines
  USE DomainMappings
  USE EquationsRoutines
  USE EquationsAccessRoutines
  USE EquationsMappingRoutines
  USE EquationsMatricesAccessRoutines
  USE EquationsMatricesRoutines
  USE EquationsSetAccessRoutines
  USE FieldRoutines
  USE FIELD_IO_ROUTINES
  USE FieldAccessRoutines
  USE FLUID_MECHANICS_IO_ROUTINES
  USE InputOutput
  USE ISO_VARYING_STRING
  USE Kinds
  USE MatrixVector
  USE ProblemAccessRoutines
  USE Strings
  USE SolverRoutines
  USE SolverAccessRoutines
  USE Timer
  USE Types

#include "macros.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC  Stokes_EquationsSetSpecificationSet

  PUBLIC  Stokes_EquationsSetSolutionMethodSet

  PUBLIC  Stokes_EquationsSetSetup

  PUBLIC  Stokes_BoundaryConditionsAnalyticCalculate

  PUBLIC  Stokes_ProblemSpecificationSet

  PUBLIC  Stokes_ProblemSetup

  PUBLIC  Stokes_FiniteElementCalculate

  PUBLIC  Stokes_PostSolve

  PUBLIC  Stokes_PreSolve

  PUBLIC  Stokes_AnalyticFunctions

CONTAINS

  !
  !================================================================================================================================
  !
  
  !>Sets/changes the solution method for a Stokes flow equation type of an fluid mechanics equations set class.
  SUBROUTINE Stokes_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: solutionMethod !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: esSpecification(3)
    TYPE(VARYING_STRING) :: localError

    ENTERS("Stokes_EquationsSetSolutionMethodSet",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
    
    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_STATIC_STOKES_SUBTYPE, &
      & EQUATIONS_SET_LAPLACE_STOKES_SUBTYPE, &
      & EQUATIONS_SET_TRANSIENT_STOKES_SUBTYPE, &
      & EQUATIONS_SET_ALE_STOKES_SUBTYPE)
      SELECT CASE(solutionMethod)
      CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
        equationsSet%solutionMethod=EQUATIONS_SET_FEM_SOLUTION_METHOD
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
      localError="The third equations set specification of "//TRIM(NumberToVstring(esSpecification(3),"*",err,error))// &
        & " is not valid for a Stokes flow equation of a fluid mechanics equations set."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Stokes_EquationsSetSolutionMethodSet")
    RETURN
999 ERRORSEXITS("Stokes_EquationsSetSolutionMethodSet",err,error)
    RETURN 1

  END SUBROUTINE Stokes_EquationsSetSolutionMethodSet

!
!================================================================================================================================
!

  !>Sets the equation specification for a Stokes flow equation of a fluid mechanics equations set.
  SUBROUTINE Stokes_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the specification for
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations set specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: subtype
    TYPE(VARYING_STRING) :: localError

    ENTERS("Stokes_EquationsSetSpecificationSet",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(ALLOCATED(equationsSet%specification)) CALL FlagError("Equations set specification is already allocated.",err,error,*999)    
    IF(SIZE(specification,1)<3) THEN
      localError="The size of the specified specification array of "// &
        & TRIM(NumberToVString(SIZE(specification,1),"*",err,error))//" is invalid. The size should be >= 3."
      CALL FlagError(localError,err,error,*999)
    END IF
    
    subtype=specification(3)
    SELECT CASE(subtype)
    CASE(EQUATIONS_SET_STATIC_STOKES_SUBTYPE, &
      & EQUATIONS_SET_LAPLACE_STOKES_SUBTYPE, &
      & EQUATIONS_SET_TRANSIENT_STOKES_SUBTYPE, &
      & EQUATIONS_SET_ALE_STOKES_SUBTYPE)
      !ok
    CASE(EQUATIONS_SET_OPTIMISED_STOKES_SUBTYPE)
      CALL FlagError("Not implemented yet.",err,error,*999)
    CASE DEFAULT
      localError="The third equations set specification of "//TRIM(NumberToVstring(specification(3),"*",err,error))// &
        & " is not valid for Stokes flow of a fluid mechanics equations set."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    !Set full specification
    ALLOCATE(equationsSet%specification(3),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate equations set specification.",err,error,*999)
    equationsSet%specification(1:3)=[EQUATIONS_SET_FLUID_MECHANICS_CLASS,EQUATIONS_SET_STOKES_EQUATION_TYPE,subtype]

    EXITS("Stokes_EquationsSetSpecificationSet")
    RETURN
999 ERRORS("Stokes_EquationsSetSpecificationSet",err,error)
    EXITS("Stokes_EquationsSetSpecificationSet")
    RETURN 1

  END SUBROUTINE Stokes_EquationsSetSpecificationSet

  !
  !================================================================================================================================
  !
  
  !>Sets up the standard Stokes fluid setup.
  SUBROUTINE Stokes_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to setup
    TYPE(EquationsSetSetupType), INTENT(INOUT) :: equationsSetSetup !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG):: componentIdx,dependentNumberOfVariables,dependentNumberOfComponents,esSpecification(3), &
      & geometricComponentNumber,geometricMeshComponent,geometricScalingType,independentNumberOfVariables, &
      & independentNumberOfComponents,lumpingType,materialsNumberOfComponents,materialsNumberOfVariables,numberOfDimensions, &
      & solutionMethod,sparsityType
    TYPE(DecompositionType), POINTER :: geometricDecomposition
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsSetAnalyticType), POINTER :: equationsAnalytic
    TYPE(EquationsSetIndependentType), POINTER :: equationsIndependent
    TYPE(EquationsSetMaterialsType), POINTER :: equationsMaterials
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField,geometricField
    TYPE(RegionType), POINTER :: region
    TYPE(VARYING_STRING) :: localError

    ENTERS("Stokes_EquationsSetSetup",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
    
    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_STATIC_STOKES_SUBTYPE, &
      & EQUATIONS_SET_TRANSIENT_STOKES_SUBTYPE, &
      & EQUATIONS_SET_LAPLACE_STOKES_SUBTYPE, &
      & EQUATIONS_SET_ALE_STOKES_SUBTYPE)
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
        CALL Stokes_EquationsSetSolutionMethodSet(equationsSet,EQUATIONS_SET_FEM_SOLUTION_METHOD,err,error,*999)
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
!!TODO: Check valid setup
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType, &
          & "*",err,error))// " for a setup type of "//TRIM(NumberToVString(equationsSetSetup% &
          & setupType,"*",err,error))// " is invalid for a standard Stokes fluid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
      !-----------------------------------------------------------------
      ! G e o m e t r i c   f i e l d
      !-----------------------------------------------------------------
      !Do nothing???
    CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
      !-----------------------------------------------------------------
      ! D e p e n d e n t   f i e l d
      !-----------------------------------------------------------------
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
      dependentNumberOfVariables=2 ! U and dUdN
      dependentNumberOfComponents=numberOfDimensions+1 !number of dimensions for u plus one for pressure
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        IF(equationsSet%dependent%dependentFieldAutoCreated) THEN
          !Create the auto created dependent field
          !start field creation with name 'U'
          CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsSet%dependent%dependentField,err,error,*999)
          !start creation of a new field
          CALL Field_TypeSetAndLock(equationsSet%dependent%dependentField,FIELD_GENERAL_TYPE,err,error,*999)
          !label the field
          CALL Field_LabelSetAndLock(equationsSet%dependent%dependentField,"U",err,error,*999)
          !define new created field to be dependent
          CALL Field_DependentTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_DEPENDENT_TYPE,err,error,*999)
          !look for decomposition rule already defined
          CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
          !apply decomposition rule found on new created field
          CALL Field_DecompositionSetAndLock(equationsSet%dependent%dependentField,geometricDecomposition,err,error,*999)
          !point new field to geometric field
          CALL Field_GeometricFieldSetAndLock(equationsSet%dependent%dependentField,geometricField,err,error,*999)
          !set number of variables to 2 (1 for U and one for DELUDELN)
          CALL Field_NumberOfVariablesSetAndLock(equationsSet%dependent%dependentField,dependentNumberOfVariables,err,error,*999)
          CALL Field_VariableTypesSetAndLock(equationsSet%dependent%dependentField,[FIELD_U_VARIABLE_TYPE, &
            & FIELD_DELUDELN_VARIABLE_TYPE],err,error,*999)
          CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
            & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
            & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE, &
            & err,error,*999)
          !calculate number of components with one component for each dimension and one for pressure
          CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
            & dependentNumberOfComponents,err,error,*999)
          CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
            & dependentNumberOfComponents,err,error,*999)
          CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricMeshComponent,err,error,*999)
          !Default to the geometric interpolation setup
          DO componentIdx=1,dependentNumberOfComponents
            CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,componentIdx, &
              & geometricMeshComponent,err,error,*999)
            CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
              & componentIdx,geometricMeshComponent,err,error,*999)
          ENDDO !componentIdx
          CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
          SELECT CASE(solutionMethod)
          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
            DO componentIdx=1,dependentNumberOfComponents
              CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
                & componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
                & componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
            ENDDO !componentIdx
            CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
            CALL Field_ScalingTypeSet(equationsSet%dependent%dependentField,geometricScalingType,err,error,*999)
          CASE DEFAULT
            !Other solutions not defined yet
            localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE
          !Check the user specified field
          CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
          CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_DEPENDENT_TYPE,err,error,*999)
          CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,dependentNumberOfVariables,err,error,*999)
          CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE], &
            & err,error,*999)
          CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
            & err,error,*999)
          CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,dependentNumberOfComponents, &
            & err,error,*999)
          CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,dependentNumberOfComponents, &
            & err,error,*999)
          CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
          SELECT CASE(solutionMethod)
          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
            DO componentIdx=1,numberOfDimensions
              CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,componentIdx, &
                & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,componentIdx, &
                & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
            ENDDO !componentIdx
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
          & " is invalid for a standard Stokes fluid"
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
      !-----------------------------------------------------------------
      ! M a t e r i a l s   f i e l d 
      !-----------------------------------------------------------------
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      NULLIFY(equationsMaterials)
      CALL EquationsSet_MaterialsGet(equationsSet,equationsMaterials,err,error,*999)
      materialsNumberOfVariables=1
      materialsNumberOfComponents=2
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        IF(equationsMaterials%materialsFieldAutoCreated) THEN
          !Create the auto created materials field
          !start field creation with name 'Materials Field'
          CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsMaterials%materialsField,err,error,*999)
          CALL Field_TypeSetAndLock(equationsMaterials%materialsField,FIELD_MATERIAL_TYPE,err,error,*999)
          !label the field
          CALL Field_LabelSetAndLock(equationsMaterials%materialsField,"Materials Field",err,error,*999)
          CALL Field_DependentTypeSetAndLock(equationsMaterials%materialsField,FIELD_INDEPENDENT_TYPE,err,error,*999)
          CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
          !apply decomposition rule found on new created field
          CALL Field_DecompositionSetAndLock(equationsMaterials%materialsField,geometricDecomposition,err,error,*999)
          !point new field to geometric field
          CALL Field_GeometricFieldSetAndLock(equationsMaterials%materialsField,geometricField,err,error,*999)
          CALL Field_NumberOfVariablesSet(equationsMaterials%materialsField,materialsNumberOfVariables,err,error,*999)
          CALL Field_VariableTypesSetAndLock(equationsMaterials%materialsField,[FIELD_U_VARIABLE_TYPE],err,error,*999)
          CALL Field_DimensionSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
            & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_NumberOfComponentsSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
            & materialsNumberOfComponents,err,error,*999)
          CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricComponentNumber,err,error,*999)
          DO componentIdx=1,materialsNumberOfComponents
            CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,componentIdx, &
              & geometricComponentNumber,err,error,*999)
            CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,componentIdx, &
              & FIELD_CONSTANT_INTERPOLATION,err,error,*999)
          ENDDO !componentIdx
          !Default the field scaling to that of the geometric field
          CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
          CALL Field_ScalingTypeSet(equationsMaterials%materialsField,geometricScalingType,err,error,*999)
        ELSE
          !Check the user specified field
          CALL Field_TypeCheck(equationsSetSetup%field,FIELD_MATERIAL_TYPE,err,error,*999)
          CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
          CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,materialsNumberOfVariables,err,error,*999)
          CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE],err,error,*999)
          CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,materialsNumberOfComponents, &
            & err,error,*999)
        ENDIF
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        IF(equationsMaterials%materialsFieldAutoCreated) THEN
          !Finish creating the materials field
          CALL Field_CreateFinish(equationsMaterials%materialsField,err,error,*999)
          !Set the default values for the materials field
          !First set the mu values to 0.001
          !materialsNumberOfComponents
          ! viscosity=1
          CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,1,1.0_DP,err,error,*999)
          ! density=2
          CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,2,100.0_DP,err,error,*999)
        ENDIF
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for Stokes equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_INDEPENDENT_TYPE)
      !-----------------------------------------------------------------
      ! I n d e p e n d e n t   f i e l d
      !-----------------------------------------------------------------
      !define an independent field for ALE information
      SELECT CASE(esSpecification(3))
      CASE(EQUATIONS_SET_ALE_STOKES_SUBTYPE)
        NULLIFY(equationsIndependent)
        CALL EquationsSet_IndependentGet(equationsSet,equationsIndependent,err,error,*999)
        NULLIFY(geometricField)
        CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
        CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
        independentNumberOfVariables=1
        independentNumberOfComponents=numberOfDimensions
        SELECT CASE(equationsSetSetup%actionType)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
          IF(equationsIndependent%independentFieldAutoCreated) THEN
            !Create the auto created independent field
            !start field creation with name 'Independent Field'
            CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsIndependent%independentField,err,error,*999)
            !start creation of a new field
            CALL Field_TypeSetAndLock(equationsIndependent%independentField,FIELD_GENERAL_TYPE,err,error,*999)
            !label the field
            CALL Field_LabelSetAndLock(equationsIndependent%independentField,"Independent Field",err,error,*999)
            !define new created field to be independent
            CALL Field_DependentTypeSetAndLock(equationsIndependent%independentField,FIELD_INDEPENDENT_TYPE,err,error,*999)
            !look for decomposition rule already defined
            CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
            !apply decomposition rule found on new created field
            CALL Field_DecompositionSetAndLock(equationsIndependent%independentField,geometricDecomposition,err,error,*999)
            !point new field to geometric field
            CALL Field_GeometricFieldSetAndLock(equationsIndependent%independentField,geometricField,err,error,*999)
            !set number of variables to 1 (1 for U)
            CALL Field_NumberOfVariablesSetAndLock(equationsIndependent%independentField,independentNumberOfVariables, &
              & err,error,*999)
            CALL Field_VariableTypesSetAndLock(equationsIndependent%independentField,[FIELD_U_VARIABLE_TYPE],err,error,*999)
            CALL Field_DimensionSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DataTypeSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE, &
              & err,error,*999)
            !calculate number of components with one component for each dimension
            CALL Field_NumberOfComponentsSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
              & independentNumberOfComponents,err,error,*999)
            CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricMeshComponent,err,error,*999)
            !Default to the geometric interpolation setup
            DO componentIdx=1,independentNumberOfComponents
              CALL Field_ComponentMeshComponentSet(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
                & componentIdx,geometricMeshComponent,err,error,*999)
            ENDDO !componentIdx
            CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)          
            SELECT CASE(solutionMethod)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              DO componentIdx=1,independentNumberOfComponents
                CALL Field_ComponentInterpolationSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
                  & componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              ENDDO !componentIdx
              CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
              CALL Field_ScalingTypeSet(equationsIndependent%independentField,geometricScalingType,err,error,*999)
              !Other solutions not defined yet
            CASE DEFAULT
              localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ELSE
            !Check the user specified field
            CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
            CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
            CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,independentNumberOfVariables,err,error,*999)
            CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE],err,error,*999)
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,independentNumberOfComponents, &
              & err,error,*999)
            CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)          
            SELECT CASE(solutionMethod)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              DO componentIdx=1,independentNumberOfComponents
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,componentIdx, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              ENDDO !componentIdx
            CASE DEFAULT
              localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ENDIF
        CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
          IF(equationsIndependent%independentFieldAutoCreated) THEN
            CALL Field_CreateFinish(equationsIndependent%independentField,err,error,*999)
            CALL Field_ParameterSetEnsureCreated(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_MESH_DISPLACEMENT_SET_TYPE,err,error,*999)
            CALL Field_ParameterSetEnsureCreated(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_MESH_VELOCITY_SET_TYPE,err,error,*999)
            CALL Field_ParameterSetEnsureCreated(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_BOUNDARY_SET_TYPE,err,error,*999)
          ENDIF
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
            & " is invalid for a standard Stokes fluid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_STATIC_STOKES_SUBTYPE, &
        & EQUATIONS_SET_TRANSIENT_STOKES_SUBTYPE, &
        & EQUATIONS_SET_LAPLACE_STOKES_SUBTYPE)
        !Do nothing
      CASE DEFAULT
        localError="The third equations set specification of "//TRIM(NumberToVstring(esSpecification(3),"*",err,error))// &
          & " is invalid for a Stokes flow equations set."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
      !-----------------------------------------------------------------
      ! S o u r c e   f i e l d 
      !-----------------------------------------------------------------
!!TODO: INCLUDE GRAVITY AS SOURCE TYPE
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        !Do nothing
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        !Do nothing
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a standard Stokes fluid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
      !-----------------------------------------------------------------
      ! A n a l y t i c   f i e l d
      !-----------------------------------------------------------------
      SELECT CASE(esSpecification(3))
      CASE(EQUATIONS_SET_STATIC_STOKES_SUBTYPE,EQUATIONS_SET_LAPLACE_STOKES_SUBTYPE, &
        & EQUATIONS_SET_TRANSIENT_STOKES_SUBTYPE)
        NULLIFY(equationsAnalytic)
        CALL EquationsSet_AnalyticGet(equationsSet,equationsAnalytic,err,error,*999)
        SELECT CASE(equationsSetSetup%actionType)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
          CALL EquationsSet_AssertDependentIsFinished(equationsSet,err,error,*999)
          NULLIFY(geometricField)
          CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
          CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
          NULLIFY(dependentField)
          CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
          SELECT CASE(numberOfDimensions)
          CASE(1)
            CALL FlagError("No analytic functions are defined for a one dimensional Stokes equation.",err,error,*999)
          CASE(2)
            SELECT CASE(equationsSetSetup%analyticFunctionType)
            CASE(EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_1)
              !Set analtyic function type
              equationsAnalytic%analyticFunctionType=EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_1
            CASE(EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_2)
              !Set analtyic function type
              equationsAnalytic%analyticFunctionType=EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_2
            CASE(EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_3)
              !Set analtyic function type
              equationsAnalytic%analyticFunctionType=EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_3
            CASE(EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_4)
              !Set analtyic function type
              equationsAnalytic%analyticFunctionType=EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_4
            CASE(EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_5)
              !Set analtyic function type
              equationsAnalytic%analyticFunctionType=EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_5
            CASE DEFAULT
              localError="The specified analytic function type of "// &
                & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
                & " is invalid for a two-dimensional analytic Stokes equations set."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE(3)
            SELECT CASE(equationsSetSetup%analyticFunctionType)
            CASE(EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_1)
              !Set analtyic function type
              equationsAnalytic%analyticFunctionType=EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_1
            CASE(EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_2)
              !Set analtyic function type
              equationsAnalytic%analyticFunctionType=EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_2
            CASE(EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_3)
              !Set analtyic function type
              equationsAnalytic%analyticFunctionType=EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_3
            CASE(EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_4)
              !Set analtyic function type
              equationsAnalytic%analyticFunctionType=EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_4
            CASE(EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_5)
              !Set analtyic function type
              equationsAnalytic%analyticFunctionType=EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_5
            CASE DEFAULT
              localError="The specified analytic function type of "// &
                & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
                & " is invalid for a three-dimensional analytic Stokes equations set."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The number of dimensions of "//TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
              & " is invalid for an analytic Stokes equations set."
          END SELECT
        CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
          IF(equationsAnalytic%analyticFieldAutoCreated) THEN
            CALL Field_CreateFinish(equationsAnalytic%analyticField,err,error,*999)
          ENDIF
!!TODO: Set default values
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
            & " is invalid for an analytic Stokes problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The third equations set specification of "// &
          & TRIM(NumberToVString(equationsSet%specification(3),"*",err,error))// &
          " is invalid for a Stokes flow equations set."
        CALL FLAG_ERROR(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        CALL EquationsSet_AssertMaterialsIsFinished(equationsSet,err,error,*999)
        NULLIFY(equations)
        CALL Equations_CreateStart(equationsSet,equations,err,error,*999)
        CALL Equations_EquationTypeSet(equations,EQUATIONS_VECTOR_TYPE,err,error,*999)
        CALL Equations_LinearityTypeSet(equations,EQUATIONS_LINEAR,err,error,*999)
        SELECT CASE(esSpecification(3))
        CASE(EQUATIONS_SET_STATIC_STOKES_SUBTYPE,EQUATIONS_SET_LAPLACE_STOKES_SUBTYPE)
          CALL Equations_TimeDependenceTypeSet(equations,EQUATIONS_STATIC,err,error,*999)
        CASE(EQUATIONS_SET_TRANSIENT_STOKES_SUBTYPE,EQUATIONS_SET_ALE_STOKES_SUBTYPE)
          CALL Equations_TimeDependenceTypeSet(equations,EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
        CASE DEFAULT
          localError="The third equations set specification of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
            & " is invalid for a Stokes flow equations set."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)          
        SELECT CASE(solutionMethod)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          !Finish the equations creation
          NULLIFY(equations)
          CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
          CALL Equations_CreateFinish(equations,err,error,*999)
          NULLIFY(vectorEquations)
          CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
          !Create the equations mapping.
          NULLIFY(vectorMapping)
          CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_U_VARIABLE_TYPE,vectorMapping,err,error,*999)
          SELECT CASE(esSpecification(3))
          CASE(EQUATIONS_SET_STATIC_STOKES_SUBTYPE,EQUATIONS_SET_LAPLACE_STOKES_SUBTYPE)
            CALL EquationsMappingVector_NumberOfLinearMatricesSet(vectorMapping,1,err,error,*999)
            CALL EquationsMappingVector_LinearMatricesVariableTypesSet(vectorMapping,[FIELD_U_VARIABLE_TYPE],err,error,*999)
          CASE(EQUATIONS_SET_TRANSIENT_STOKES_SUBTYPE,EQUATIONS_SET_ALE_STOKES_SUBTYPE)
            CALL EquationsMappingVector_DynamicMatricesSet(vectorMapping,.TRUE.,.TRUE.,err,error,*999)
            CALL EquationsMappingVector_DynamicVariableTypeSet(vectorMapping,FIELD_U_VARIABLE_TYPE,err,error,*999)
          CASE DEFAULT
            localError="The third equations set specification of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
              & " is invalid for a Stokes flow equations set."
            CALL FlagError(localError,err,error,*999)
          END SELECT
          CALL EquationsMappingVector_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE,err,error,*999)
          CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
          !Create the equations matrices
          NULLIFY(vectorMatrices)
          CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
          !Set up matrix storage and structure
          CALL Equations_SparsityTypeGet(equations,sparsityType,err,error,*999)
          SELECT CASE(esSpecification(3))
          CASE(EQUATIONS_SET_STATIC_STOKES_SUBTYPE,EQUATIONS_SET_LAPLACE_STOKES_SUBTYPE)
            SELECT CASE(sparsityType)
            CASE(EQUATIONS_MATRICES_FULL_MATRICES)
              CALL EquationsMatricesVector_LinearStorageTypeSet(vectorMatrices,[MATRIX_BLOCK_STORAGE_TYPE], err,error,*999)
            CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
              CALL EquationsMatricesVector_LinearStorageTypeSet(vectorMatrices,[MATRIX_COMPRESSED_ROW_STORAGE_TYPE], &
                & err,error,*999)
              CALL EquationsMatricesVector_LinearStructureTypeSet(vectorMatrices,[EQUATIONS_MATRIX_FEM_STRUCTURE],err,error,*999)
            CASE DEFAULT
              localError="The equations matrices sparsity type of "//TRIM(NumberToVString(sparsityType,"*",err,error))// &
                & " is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE(EQUATIONS_SET_TRANSIENT_STOKES_SUBTYPE,EQUATIONS_SET_ALE_STOKES_SUBTYPE)
            CALL Equations_LumpingTypeGet(equations,lumpingType,err,error,*999)
            IF(lumpingType==EQUATIONS_LUMPED_MATRICES) THEN
              !Set up lumping
              CALL EquationsMatricesVector_DynamicLumpingTypeSet(vectorMatrices,[EQUATIONS_MATRIX_UNLUMPED, &
                & EQUATIONS_MATRIX_LUMPED],err,error,*999)
              CALL EquationsMatricesVector_DynamicStorageTypeSet(vectorMatrices, &
                & [DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE], &
                & err,error,*999)
              CALL EquationsMatricesVector_DynamicStructureTypeSet(vectorMatrices, &
                & [EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_DIAGONAL_STRUCTURE],err,error,*999)
            ELSE
              SELECT CASE(sparsityType)
              CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                CALL EquationsMatricesVector_DynamicStorageTypeSet(vectorMatrices, &
                  & [DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE,DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE],err,error,*999)
              CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                CALL EquationsMatricesVector_DynamicStorageTypeSet(vectorMatrices, &
                  & [DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE, &
                  & DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE],err,error,*999)
                CALL EquationsMatricesVector_DynamicStructureTypeSet(vectorMatrices, &
                  & [EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_FEM_STRUCTURE],err,error,*999)
              CASE DEFAULT
                localError="The equations matrices sparsity type of "//TRIM(NumberToVString(sparsityType,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ENDIF
          CASE DEFAULT
            localError="The third equations set specification of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
              & " is invalid for a Stokes flow equations set."
            CALL FlagError(localError,err,error,*999)
          END SELECT
          CALL EquationsMatrices_VectorCreateFinish(vectorMatrices,err,error,*999)
        CASE DEFAULT
          localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a Stokes flow equation set."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
        & " is invalid for a standard Stokes fluid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("Stokes_EquationsSetSetup")
    RETURN
999 ERRORSEXITS("Stokes_EquationsSetSetup",err,error)
    RETURN 1
    
  END SUBROUTINE Stokes_EquationsSetSetup

  !
  !================================================================================================================================
  !
  
  !>Sets the problem specification for a Stokes fluid problem.
  SUBROUTINE Stokes_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to set the problem specification for
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: problemSubtype
    TYPE(VARYING_STRING) :: localError

    ENTERS("Stokes_ProblemSpecificationSet",err,error,*999)

    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Problem is not associated.",err,error,*999)
    IF(ALLOCATED(problem%specification)) CALL FlagError("Problem specification is already allocated.",err,error,*999)
    IF(SIZE(problemSpecification,1)<3) THEN
      localError="The size of the specified problem specification array of "// &
        & TRIM(NumberToVString(SIZE(problemSpecification,1),"*",err,error))// &
        & " is invalid. The size should be >= 3."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    problemSubtype=problemSpecification(3)
    SELECT CASE(problemSubtype)
    CASE(PROBLEM_STATIC_STOKES_SUBTYPE, &
      & PROBLEM_LAPLACE_STOKES_SUBTYPE, &
      & PROBLEM_TRANSIENT_STOKES_SUBTYPE, &
      & PROBLEM_ALE_STOKES_SUBTYPE)
      !All ok
    CASE(PROBLEM_OPTIMISED_STOKES_SUBTYPE)
      CALL FlagError("Not implemented yet.",err,error,*999)
    CASE DEFAULT
      localError="The third problem specification of "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
        & " is not valid for a Stokes flow fluid mechanics problem."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    ALLOCATE(problem%specification(3),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
    problem%specification(1:3)=[PROBLEM_FLUID_MECHANICS_CLASS,PROBLEM_STOKES_EQUATION_TYPE,problemSubtype]

    EXITS("Stokes_ProblemSpecificationSet")
    RETURN
999 ERRORSEXITS("Stokes_ProblemSpecificationSet",err,error)
    RETURN 1

  END SUBROUTINE Stokes_ProblemSpecificationSet

  !
  !================================================================================================================================
  !
  
  !>Sets up the Stokes problem.
  SUBROUTINE Stokes_ProblemSetup(problem,problemSetup,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem set to setup a Stokes fluid on.
    TYPE(ProblemSetupType), INTENT(INOUT) :: problemSetup !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: pSpecification(3)
    TYPE(ControlLoopType), POINTER :: controlLoop,controlLoopRoot
    TYPE(SolverType), POINTER :: solver, meshSolver
    TYPE(SolverEquationsType), POINTER :: solverEquations,meshSolverEquations
    TYPE(SolversType), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError

    ENTERS("Stokes_ProblemSetup",err,error,*999)

    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)

    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STATIC_STOKES_SUBTYPE, &
      & PROBLEM_LAPLACE_STOKES_SUBTYPE, &
      & PROBLEM_TRANSIENT_STOKES_SUBTYPE, &
      & PROBLEM_ALE_STOKES_SUBTYPE)
      !ok
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a Stokes fluid type of a fluid mechanics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    SELECT CASE(problemSetup%setupType)
    CASE(PROBLEM_SETUP_INITIAL_TYPE)
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)
        !Do nothing????
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Do nothing???
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a standard Stokes fluid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_SETUP_CONTROL_TYPE)
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)
        !Set up a simple control loop
        NULLIFY(controlLoop)
        CALL ControlLoop_CreateStart(problem,controlLoop,err,error,*999)
        IF(pSpecification(3)/=PROBLEM_STATIC_STOKES_SUBTYPE.AND. &
          & pSpecification(3)/=PROBLEM_LAPLACE_STOKES_SUBTYPE) THEN
          CALL ControlLoop_TypeSet(controlLoop,CONTROL_TIME_LOOP_TYPE,err,error,*999)
        ENDIF
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Finish the control loops
        NULLIFY(controlLoopRoot)
        CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
        NULLIFY(controlLoop)
        CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
        CALL ControlLoop_CreateFinish(controlLoop,err,error,*999)
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a standard Stokes fluid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_SETUP_SOLVERS_TYPE)
      !Get the control loop
      NULLIFY(controlLoopRoot)
      CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
      NULLIFY(controlLoop)
      CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)
        !Start the solvers creation
        NULLIFY(solvers)
        CALL Solvers_CreateStart(controlLoop,solvers,err,error,*999)
        SELECT CASE(pSpecification(3))
        CASE(PROBLEM_STATIC_STOKES_SUBTYPE, &
          & PROBLEM_LAPLACE_STOKES_SUBTYPE)
          CALL Solvers_NumberOfSolversSet(solvers,1,err,error,*999)
          !Set the solver to be a linear solver
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          CALL Solver_TypeSet(solver,SOLVER_LINEAR_TYPE,err,error,*999)
          !Set solver defaults
          CALL Solver_LibraryTypeSet(solver,SOLVER_PETSC_LIBRARY,err,error,*999)
        CASE(PROBLEM_TRANSIENT_STOKES_SUBTYPE)
          CALL Solvers_NumberOfSolversSet(solvers,1,err,error,*999)
          !Set the solver to be a first order dynamic solver
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          CALL Solver_TypeSet(solver,SOLVER_DYNAMIC_TYPE,err,error,*999)
          CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
          !Set solver defaults
          CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
          CALL Solver_DynamicSchemeSet(solver,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,err,error,*999)
          CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
        CASE(PROBLEM_ALE_STOKES_SUBTYPE)
          CALL Solvers_NumberOfSolversSet(solvers,2,err,error,*999)
          !Set the first solver to be a linear solver for the Laplace mesh movement problem
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          CALL Solver_TypeSet(solver,SOLVER_LINEAR_TYPE,err,error,*999)
          !Set solver defaults
          CALL Solver_LibraryTypeSet(solver,SOLVER_PETSC_LIBRARY,err,error,*999)
          !Set the second solver to be a first order dynamic solver
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
          CALL Solver_TypeSet(solver,SOLVER_DYNAMIC_TYPE,err,error,*999)
          CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
          !Set solver defaults
          CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
          CALL Solver_DynamicSchemeSet(solver,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,err,error,*999)
          CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
        CASE DEFAULT
          localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
            & " is not valid for a Stokes fluid type of a fluid mechanics problem class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Get the solvers
        NULLIFY(solvers)
        CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
        !Finish the solvers creation
        CALL Solvers_CreateFinish(solvers,err,error,*999)
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a standard Stokes fluid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
      !Get the control loop
      NULLIFY(controlLoopRoot)
      CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
      NULLIFY(controlLoop)
      CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
      NULLIFY(solvers)
      CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)
        SELECT CASE(pSpecification(3))
        CASE(PROBLEM_STATIC_STOKES_SUBTYPE, &
          & PROBLEM_LAPLACE_STOKES_SUBTYPE)
          !Get the solver
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          !Create the solver equations
          NULLIFY(solverEquations)
          CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
          CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_LINEAR,err,error,*999)
          CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_STATIC,err,error,*999)
          CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
        CASE(PROBLEM_TRANSIENT_STOKES_SUBTYPE)
          !Get the solver
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          !Create the solver equations
          CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
          CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_LINEAR,err,error,*999)
          CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,&
            & err,error,*999)
          CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
        CASE(PROBLEM_ALE_STOKES_SUBTYPE)
          !Get the solver
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          !Create the solver equations
          NULLIFY(solverEquations)
          CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
          CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_LINEAR,err,error,*999)
          CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_STATIC,err,error,*999)
          CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
          !Get the solver
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
          !Create the solver equations
          CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
          CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_LINEAR,err,error,*999)
          CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
          CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
            & " is invalid for a standard Stokes fluid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Get the solver equations
        CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
        CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
        !Finish the solver equations creation
        CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)
        IF(pSpecification(3)==PROBLEM_ALE_STOKES_SUBTYPE) THEN
          !Get the solver equations
          CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
          CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
          !Finish the solver equations creation
          CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)
        ENDIF
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a Stokes fluid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
        & " is invalid for a Stokes fluid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Stokes_ProblemSetup")
    RETURN
999 ERRORSEXITS("Stokes_ProblemSetup",err,error)
    RETURN 1
    
  END SUBROUTINE Stokes_ProblemSetup

!
!================================================================================================================================
!

  !>Calculates the element stiffness matrices and RHS for a Stokes fluid finite element equations set.
  SUBROUTINE Stokes_FiniteElementCalculate(equationsSet,elementNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) analyticFunctionType,colsVariableType,columnComponentIdx,columnElementDOFIdx,columnElementParameterIdx, &
      & columnXiIdx,componentIdx,esSpecification(3),gaussPointIdx,maxColumnElementDOFIdx,maxRowElementDOFIdx, &
      & minColumnElementDOFIdx,minRowElementDOFIdx,numberOfColsComponents,numberOfColumnElementParameters,numberOfDimensions, &
      & numberOfGauss,numberOfRowsComponents,numberOfRowElementParameters,numberOfXi,out,rowComponentIdx,rowElementDOFIdx, &
      & rowElementParameterIdx,rowXiIdx,rowsVariableType,scalingType,variableType,xv
    REAL(DP) :: columnPhi,columndPhidXi(3),gaussWeight,jacobian,jacobianGaussWeight,dXidX(3,3),muParam,rhoParam,rowPhi, &
      & rowdPhidXi(3),sum,wValue(3),x(3)
    REAL(DP) :: aMatrix(256,256),dMatrix(256,256),aleMatrix(256,256),bTMatrix(256,256)
    LOGICAL :: update,updateDamping,updateMatrices,updateStiffness,updateRHS,updateSource
    TYPE(BasisType), POINTER :: columnBasis,dependentBasis,geometricBasis,independentBasis,rowBasis
    TYPE(DecompositionType), POINTER :: dependentDecomposition,geometricDecomposition
    TYPE(DomainType), POINTER :: columnDomain,dependentDomain,geometricDomain,rowDomain
    TYPE(DomainElementsType), POINTER :: columnDomainElements,dependentDomainElements,geometricDomainElements,rowDomainElements
    TYPE(DomainTopologyType), POINTER :: columnDomainTopology,dependentDomainTopology,geometricDomainTopology,rowDomainTopology
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping
    TYPE(EquationsMappingSourcesType), POINTER :: sourcesMapping
    TYPE(EquationsMappingSourceType), POINTER :: sourceMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourcesType), POINTER :: sourceVectors
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatrixType), POINTER :: stiffnessMatrix,dampingMatrix
    TYPE(EquationsSetAnalyticType), POINTER :: equationsAnalytic
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField,geometricField,materialsField,independentField
    TYPE(FieldInterpolationParametersType), POINTER :: colsInterpParameters,geometricInterpParameters,independentInterpParameters, &
      & materialsInterpParameters,rowsInterpParameters
    TYPE(FieldInterpolatedPointType), POINTER :: geometricInterpPoint,independentInterpPoint,materialsInterpPoint
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: geometricInterpPointMetrics
    TYPE(FieldVariableType), POINTER :: colsVariable,dependentVariable,geometricVariable,rowsVariable
    TYPE(QuadratureSchemeType), POINTER :: columnQuadratureScheme,dependentQuadratureScheme,geometricQuadratureScheme, &
      & rowQuadratureScheme
    TYPE(VARYING_STRING) :: localError

!\todo: Reduce number of variables and parameters

    ENTERS("Stokes_FiniteElementCalculate",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)

    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_STATIC_STOKES_SUBTYPE, &
      & EQUATIONS_SET_LAPLACE_STOKES_SUBTYPE, &
      & EQUATIONS_SET_TRANSIENT_STOKES_SUBTYPE, &
      & EQUATIONS_SET_ALE_STOKES_SUBTYPE)
      !OK
    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is not valid for a Stokes fluid type of a fluid mechanics equations set class."
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
    NULLIFY(rhsMapping)
    CALL EquationsMappingVector_RHSMappingExists(vectorMapping,rhsMapping,err,error,*999)
    NULLIFY(sourcesMapping)
    CALL EquationsMappingVector_SourcesMappingExists(vectorMapping,sourcesMapping,err,error,*999)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
    NULLIFY(rhsVector)
    updateRHS=.FALSE.
    IF(ASSOCIATED(rhsMapping)) THEN
      CALL EquationsMatricesVector_RHSVectorGet(vectorMatrices,rhsVector,err,error,*999)
      CALL EquationsMatricesRHS_UpdateVectorGet(rhsVector,updateRHS,err,error,*999)
    ENDIF
    NULLIFY(sourceMapping)
    NULLIFY(sourceVectors)
    NULLIFY(sourceVector)
    updateSource=.FALSE.
    IF(ASSOCIATED(sourcesMapping)) THEN
      CALL EquationsMappingSources_SourceMappingGet(sourcesMapping,1,sourceMapping,err,error,*999)
      CALL EquationsMatricesVector_SourceVectorsGet(vectorMatrices,sourceVectors,err,error,*999)
      CALL EquationsMatricesSources_SourceVectorGet(sourceVectors,1,sourceVector,err,error,*999)
      CALL EquationsMatricesSource_UpdateVectorGet(sourceVector,updateSource,err,error,*999)
    ENDIF
    updateStiffness=.FALSE.
    updateDamping=.FALSE.
    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_STATIC_STOKES_SUBTYPE, &
      & EQUATIONS_SET_LAPLACE_STOKES_SUBTYPE)
      NULLIFY(linearMapping)
      CALL EquationsMappingVector_LinearMappingGet(vectorMapping,linearMapping,err,error,*999)
      NULLIFY(colsVariable)
      CALL EquationsMappingLinear_LinearMatrixVariableGet(linearMapping,1,colsVariable,err,error,*999)
      NULLIFY(linearMatrices)
      CALL EquationsMatricesVector_LinearMatricesGet(vectorMatrices,linearMatrices,err,error,*999)
      NULLIFY(stiffnessMatrix)
      CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,1,stiffnessMatrix,err,error,*999)
      NULLIFY(dynamicMapping)
      NULLIFY(dynamicMatrices)
      NULLIFY(dampingMatrix)
      CALL EquationsMatrix_UpdateMatrixGet(stiffnessMatrix,updateStiffness,err,error,*999)
    CASE(EQUATIONS_SET_TRANSIENT_STOKES_SUBTYPE, &
      & EQUATIONS_SET_ALE_STOKES_SUBTYPE)
      NULLIFY(linearMapping)
      NULLIFY(linearMatrices)
      NULLIFY(dynamicMapping)
      CALL EquationsMappingVector_DynamicMappingGet(vectorMapping,dynamicMapping,err,error,*999)
      NULLIFY(colsVariable)
      CALL EquationsMappingDynamic_DynamicVariableGet(dynamicMapping,colsVariable,err,error,*999)
      NULLIFY(dynamicMatrices)
      CALL EquationsMatricesVector_DynamicMatricesGet(vectorMatrices,dynamicMatrices,err,error,*999)
      NULLIFY(stiffnessMatrix)
      CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,1,stiffnessMatrix,err,error,*999)
      NULLIFY(dampingMatrix)
      CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,2,dampingMatrix,err,error,*999)
      CALL EquationsMatrix_UpdateMatrixGet(stiffnessMatrix,updateStiffness,err,error,*999)
      CALL EquationsMatrix_UpdateMatrixGet(dampingMatrix,updateDamping,err,error,*999)
    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is not valid for a Stokes fluid type of a fluid mechanics equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    updateMatrices=(updateStiffness.OR.updateDamping)
    update=(updateMatrices.OR.updateRHS.OR.updateSource)

    IF(update) THEN

      out=0
      wValue=0.0_DP
      dXidX=0.0_DP
      x=0.0_DP
      aMatrix=0.0_DP
      dMatrix=0.0_DP
      aleMatrix=0.0_DP
      bTMatrix=0.0_DP

      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      NULLIFY(dependentField)
      CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
      NULLIFY(materialsField)
      CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)
      NULLIFY(independentField)
      IF(esSpecification(3)==EQUATIONS_SET_ALE_STOKES_SUBTYPE) THEN
        CALL EquationsSet_IndependentFieldGet(equationsSet,independentField,err,error,*999)
      ENDIF

      NULLIFY(geometricVariable)
      CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,geometricVariable,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(geometricVariable,numberOfDimensions,err,error,*999)
      NULLIFY(geometricDecomposition)
      CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
      NULLIFY(geometricDomain)
      CALL Decomposition_DomainGet(geometricDecomposition,0,geometricDomain,err,error,*999)
      NULLIFY(geometricDomainTopology)
      CALL Domain_DomainTopologyGet(geometricDomain,geometricDomainTopology,err,error,*999)
      NULLIFY(geometricDomainElements)
      CALL DomainTopology_DomainElementsGet(geometricDomainTopology,geometricDomainElements,err,error,*999)
      NULLIFY(geometricBasis)
      CALL DomainElements_ElementBasisGet(geometricDomainElements,elementNumber,geometricBasis,err,error,*999)
      CALL Basis_NumberOfXiGet(geometricBasis,numberOfXi,err,error,*999)
      NULLIFY(geometricQuadratureScheme)
      CALL Basis_QuadratureSchemeGet(geometricBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,geometricQuadratureScheme,err,error,*999)
      
      NULLIFY(dependentDecomposition)
      CALL Field_DecompositionGet(dependentField,dependentDecomposition,err,error,*999)
      NULLIFY(dependentDomain)
      CALL Decomposition_DomainGet(dependentDecomposition,0,dependentDomain,err,error,*999)
      NULLIFY(dependentDomainTopology)
      CALL Domain_DomainTopologyGet(dependentDomain,dependentDomainTopology,err,error,*999)
      NULLIFY(dependentDomainElements)
      CALL DomainTopology_DomainElementsGet(dependentDomainTopology,dependentDomainElements,err,error,*999)
      NULLIFY(dependentBasis)
      CALL DomainElements_ElementBasisGet(dependentDomainElements,elementNumber,dependentBasis,err,error,*999)
      NULLIFY(dependentQuadratureScheme)
      CALL Basis_QuadratureSchemeGet(dependentBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,dependentQuadratureScheme,err,error,*999)
      CALL BasisQuadrature_NumberOfGaussGet(dependentQuadratureScheme,numberOfGauss,err,error,*999)
      
      NULLIFY(rowsVariable)
      CALL EquationsMappingLHS_LHSVariableGet(lhsMapping,rowsVariable,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(rowsVariable,numberOfRowsComponents,err,error,*999)
      
      CALL FieldVariable_VariableTypeGet(colsVariable,colsVariableType,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(colsVariable,numberOfColsComponents,err,error,*999)     
      
      NULLIFY(geometricInterpParameters)
      CALL EquationsInterpolation_GeometricParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & geometricInterpParameters,err,error,*999)
      NULLIFY(geometricInterpPoint)
      CALL EquationsInterpolation_GeometricPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,geometricInterpPoint, &
        & err,error,*999)
      NULLIFY(geometricInterpPointMetrics)
      CALL EquationsInterpolation_GeometricPointMetricsGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & geometricInterpPointMetrics,err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,geometricInterpParameters,err,error,*999)
      
      NULLIFY(materialsInterpParameters)
      CALL EquationsInterpolation_MaterialsParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & materialsInterpParameters,err,error,*999)
      NULLIFY(materialsInterpPoint)
      CALL EquationsInterpolation_MaterialsPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,materialsInterpPoint, &
        & err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,materialsInterpParameters,err,error,*999)
      
      NULLIFY(independentInterpParameters)
      NULLIFY(independentInterpPoint)
      IF(esSpecification(3)==EQUATIONS_SET_ALE_STOKES_SUBTYPE) THEN
        CALL EquationsInterpolation_IndependentParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
          & independentInterpParameters,err,error,*999)
        CALL EquationsInterpolation_IndependentPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
          & independentInterpPoint,err,error,*999)
        CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,independentInterpParameters, &
          & err,error,*999)
      ENDIF
      
      NULLIFY(equationsAnalytic)
      CALL EquationsSet_AnalyticExists(equationsSet,equationsAnalytic,err,error,*999)
      IF(ASSOCIATED(equationsAnalytic)) THEN
        CALL EquationsSet_AnalyticFunctionTypeGet(equationsSet,analyticFunctionType,err,error,*999)
      ENDIF
      
      !Start looping over Gauss points
      DO gaussPointIdx=1,numberOfGauss
        
        CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,geometricInterpPoint, &
          & err,error,*999)
        CALL Field_InterpolatedPointMetricsCalculate(numberOfXi,geometricInterpPointMetrics,err,error,*999)
        CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,materialsInterpPoint, &
          & err,error,*999)
        IF(esSpecification(3)==EQUATIONS_SET_ALE_STOKES_SUBTYPE) THEN
          CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,independentInterpPoint, &
            & err,error,*999)
          wValue(1:numberOfDimensions)=independentInterpPoint%values(1:numberOfDimensions,NO_PART_DERIV)
        END IF
        !Define muParam, viscosity=1
        muParam=materialsInterpPoint%values(1,NO_PART_DERIV)
        !Define rhoParam, density=2
        rhoParam=materialsInterpPoint%values(2,NO_PART_DERIV)
        !Calculate partial matrices
        !\todo: Check time spent here

        !Calculate Jacobian and Gauss weight.
!!TODO: Think about symmetric problems. 
        CALL FieldInterpolatedPointsMetrics_JacobianGet(geometricInterpPointMetrics,jacobian,err,error,*999)
        CALL BasisQuadratureScheme_GaussWeightGet(geometricQuadratureScheme,gaussPointIdx,gaussWeight,err,error,*999)
        jacobianGaussWeight=jacobian*gaussWeight

        DO columnComponentIdx=1,numberOfColsComponents
          DO rowXiIdx=1,numberOfXi
            dXidX(rowXiIdx,columnComponentIdx)=geometricInterpPointMetrics%dXidX(rowXiIdx,columnComponentIdx)
          ENDDO !rowXiIdx
        ENDDO !columnXiIdx
 
        !Loop over row components
        rowElementDOFIdx=0          
        DO rowComponentIdx=1,numberOfRowsComponents
          NULLIFY(rowDomain)
          CALL FieldVariable_ComponentDomainGet(rowsVariable,rowComponentIdx,rowDomain,err,error,*999)
          NULLIFY(rowDomainTopology)
          CALL Domain_DomainTopologyGet(rowDomain,rowDomainTopology,err,error,*999)
          NULLIFY(rowDomainElements)
          CALL DomainTopology_DomainElementsGet(rowDomainTopology,rowDomainElements,err,error,*999)
          NULLIFY(rowBasis)
          CALL DomainElements_ElementBasisGet(rowDomainElements,elementNumber,rowBasis,err,error,*999)
          CALL Basis_QuadratureSchemeGet(rowBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,rowQuadratureScheme,err,error,*999)
          CALL BasisQuadratureScheme_GaussWeightGet(rowQuadratureScheme,gaussPointIdx,gaussWeight,err,error,*999)
          jacobianGaussWeight=jacobian*gaussWeight
          CALL Basis_NumberOfElementParametersGet(rowBasis,numberOfRowElementParameters,err,error,*999)
          !Loop over element rows
          DO rowElementParameterIdx=1,numberOfRowElementParameters
            rowElementDOFIdx=rowElementDOFIdx+1
            CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme,rowElementParameterIdx,NO_PART_DERIV, &
              & gaussPointIdx,rowPhi,err,error,*999)
            DO rowXiIdx=1,numberOfXi
              CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme,rowElementParameterIdx, &
                & PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(rowXiIdx),gaussPointIdx,rowdPhidXi(rowXiIdx),err,error,*999)
            ENDDO !rowXiiIdx
            IF(updateMatrices) THEN
              !Loop over element columns
              columnElementDOFIdx=0
              DO columnComponentIdx=1,numberOfColsComponents
                NULLIFY(columnDomain)
                CALL FieldVariable_ComponentDomainGet(colsVariable,columnComponentIdx,columnDomain,err,error,*999)
                NULLIFY(columnDomainTopology)
                CALL Domain_DomainTopologyGet(columnDomain,columnDomainTopology,err,error,*999)
                NULLIFY(columnDomainElements)
                CALL DomainTopology_DomainElementsGet(columnDomainTopology,columnDomainElements,err,error,*999)
                NULLIFY(columnBasis)
                CALL DomainElements_ElementBasisGet(columnDomainElements,elementNumber,columnBasis,err,error,*999)
                CALL Basis_QuadratureSchemeGet(columnBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,columnQuadratureScheme,err,error,*999)
                CALL Basis_NumberOfElementParametersGet(columnBasis,numberOfColumnElementParameters,err,error,*999)
                DO columnElementParameterIdx=1,numberOfColumnElementParameters
                  columnElementDOFIdx=columnElementDOFIdx+1
                  CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme,columnElementParameterIdx,NO_PART_DERIV, &
                    & gaussPointIdx,columnPhi,err,error,*999)
                  DO columnXiIdx=1,numberOfXi
                    CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme,columnElementParameterIdx, &
                      & PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(columnXiIdx),gaussPointIdx,columndPhidXi(columnXiIdx), &
                      & err,error,*999)
                  ENDDO !columnXiIdx
                  
                  IF(updateStiffness) THEN
                    !Laplace type
                    IF(columnComponentIdx==rowComponentIdx) THEN
                      sum=0.0_DP
                      !Calculate sum
                      DO componentIdx=1,numberOfColsComponents
                        DO rowXiIdx=1,numberOfXi
                          DO columnXiIdx=1,numberOfXi
                            sum=sum+muParam*rowdPhidXi(rowXiIdx)*dXidX(rowXiIdx,componentIdx)* &
                              & columndPhidXi(columnXiIdx)*dXidX(columnXiIdx,columnComponentIdx)
                          ENDDO !columnXiIdx
                        ENDDO !rowXiIdx
                      ENDDO !componentIdx
                      !Calculate matrix
                      aMatrix(rowElementDOFIdx,columnElementDOFIdx)=aMatrix(rowElementDOFIdx,columnElementDOFIdx)+ &
                        & sum*jacobianGaussWeight
                    ENDIF !column=row
                    IF(esSpecification(3)/=EQUATIONS_SET_LAPLACE_STOKES_SUBTYPE) THEN
                      IF(columnComponentIdx<numberOfColsComponents) THEN
                        sum=0.0_DP
                        !Calculate sum
                        DO rowXiIdx=1,numberOfXi
                          DO columnXiIdx=1,numberOfXi
                            !note rowComponentIdx/columnComponentIdx derivative in dXidX
                            sum=sum+muParam*rowdPhidXi(columnXiIdx)*dXidX(columnXiIdx,columnComponentIdx)* &
                              & columndPhidXi(rowXiIdx)*dXidX(rowXiIdx,rowComponentIdx)
                          ENDDO !columnXiIdx
                        ENDDO !rowXiIdx
                        !Calculate matrix
                        dMatrix(rowElementDOFIdx,columnElementDOFIdx)=dMatrix(rowElementDOFIdx,columnElementDOFIdx)+ &
                          & sum*jacobianGaussWeight
                      ENDIF
                    ENDIF !Not Laplacce
                    IF(esSpecification(3)==EQUATIONS_SET_ALE_STOKES_SUBTYPE) THEN
                      IF(columnComponentIdx==rowComponentIdx) THEN
                        sum=0.0_DP
                        !Calculate sum
                        DO componentIdx=1,numberOfColsComponents
                          DO columnXiIdx=1,numberOfXi
                            sum=sum-rhoParam*wValue(componentIdx)*rowPhi*columndPhidXi(columnXiIdx)*dXidX(columnXiIdx,componentIdx)
                          ENDDO !columnXiIdx
                        ENDDO !columnComponentIdx
                        !Calculate MATRIX
                        aleMatrix(rowElementDOFIdx,columnElementDOFIdx)=aleMatrix(rowElementDOFIdx,columnElementDOFIdx)+ &
                          & sum*jacobianGaussWeight
                      ENDIF
                    ENDIF !ALE
                    !Calculate pressure contribution (B transpose type)
                    !Laplae type
                    IF(columnComponentIdx==numberOfColsComponents) THEN
                      sum=0.0_DP
                      !Calculate sum
                      DO columnXiIdx=1,numberOfXi
                        sum=sum-rowdPhidXi(columnXiIdx)*dXidX(columnXiIdx,rowComponentIdx)*columnPhi
                      ENDDO !columnXiIdx
                      !Calculate matrix
                      bTMatrix(rowElementDOFIdx,columnElementDOFIdx)=bTMatrix(rowElementDOFIdx,columnElementDOFIdx)+ &
                        & sum*jacobianGaussWeight
                    ENDIF
                  ENDIF !updateStiffness
                  IF(updateDamping) THEN
                    IF(columnComponentIdx==rowComponentIdx) THEN
                      dampingMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                        & dampingMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)+ &
                        & rhoParam*rowPhi*columnPhi*jacobianGaussWeight
                    ENDIF
                  ENDIF !updateDamping
                ENDDO !columnElementParameterIdx
              ENDDO !columnComponentIdx
            ENDIF !update matrices
          ENDDO !rowElementParameterIdx
        ENDDO !rowComponentIdx
        
        !Calculate analytic RHS
        IF(ASSOCIATED(equationsAnalytic)) THEN
          IF(updateRHS) THEN
            IF(analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_1.OR. &
              & analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_2.OR. &
              & analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_3.OR. &
              & analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_4.OR. &
              & analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_5.OR. &
              & analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_1.OR. &
              & analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_2.OR. &
              & analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_3.OR. &
              & analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_4.OR. &
              & analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_5) THEN
              
              !Loop over element rows
              rowElementDOFIdx=0
              DO rowComponentIdx=1,numberOfRowsComponents
                NULLIFY(rowDomain)
                CALL FieldVariable_ComponentDomainGet(rowsVariable,rowComponentIdx,rowDomain,err,error,*999)
                NULLIFY(rowDomainTopology)
                CALL Domain_DomainTopologyGet(rowDomain,rowDomainTopology,err,error,*999)
                NULLIFY(rowDomainElements)
                CALL DomainTopology_DomainElementsGet(rowDomainTopology,rowDomainElements,err,error,*999)
                NULLIFY(rowBasis)
                CALL DomainElements_ElementBasisGet(rowDomainElements,elementNumber,rowBasis,err,error,*999)
                CALL Basis_QuadratureSchemeGet(rowBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,rowQuadratureScheme,err,error,*999)
                CALL BasisQuadratureScheme_GaussWeightGet(rowQuadratureScheme,gaussPointIdx,gaussWeight,err,error,*999)
                jacobianGaussWeight=jacobian*gaussWeight
                CALL Basis_NumberOfElementParametersGet(rowBasis,numberOfRowElementParameters,err,error,*999)
                DO rowElementParameterIdx=1,numberOfRowElementParameters
                  rowElementDOFIdx=rowElementDOFIdx+1
                  CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme,rowElementParameterIdx,NO_PART_DERIV, &
                    & gaussPointIdx,rowPhi,err,error,*999)
                  !note rowComponentIdx value derivative
                  sum=0.0_DP
                  x(1:numberOfDimensions)=geometricInterpPoint%values(1:numberOfDimensions,NO_PART_DERIV)
                  IF(analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_1) THEN
                    IF(rowComponentIdx==1) THEN
                      !Calculate sum
                      sum=0.0_DP
                    ELSE IF(rowComponentIdx==2) THEN
                      !Calculate sum
                      sum=rowPhi*(-2.0_DP*muParam/10.0_DP**2)
                    ENDIF
                  ELSE IF(analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_2) THEN
                    IF(rowComponentIdx==1) THEN
                      !Calculate sum
                      sum=0.0_DP
                    ELSE IF(rowComponentIdx==2) THEN
                      !Calculate sum
                      sum=rowPhi*(-4.0_DP*muParam/100.0_DP*EXP((x(1)-x(2))/10.0_DP))
                    ENDIF
                  ELSE IF(analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_3) THEN
                    IF(rowComponentIdx==1) THEN
                      !Calculate sum
                      sum=0.0_DP
                    ELSE IF(rowComponentIdx==2) THEN
                      !Calculate sum
                      sum=rowPhi*(16.0_DP*muParam*PI*PI/100.0_DP*COS(2.0_DP*PI*x(2)/10.0_DP)*COS(2.0_DP*PI*x(1)/10.0_DP))
                    ENDIF
                  ELSE IF(analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_4) THEN
                    !Do nothing!
                  ELSE IF(analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_5) THEN
                    !Do nothing!
                  ELSE IF(analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_4) THEN
                    !Do nothing!
                  ELSE IF(analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_5) THEN
                    !Do nothing!
                  ELSE IF(analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_1) THEN
                    IF(rowComponentIdx==1) THEN
                      !Calculate sum
                      sum=0.0_DP
                    ELSE IF(rowComponentIdx==2) THEN
                      !Calculate sum
                      sum=rowPhi*(-4.0_DP*muParam/100.0_DP)
                    ELSE IF(rowComponentIdx==3) THEN
                      !Calculate sum
                      sum=rowPhi*(-4.0_DP*muParam/100.0_DP)
                    ENDIF
                  ELSE IF(analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_2) THEN
                    IF(rowComponentIdx==1) THEN
                      !Calculate sum
                      sum=0.0_DP
                    ELSE IF(rowComponentIdx==2) THEN
                      !Calculate sum
                      sum=rowPhi*(-2.0_DP*muParam/100.0_DP*(2.0_DP*EXP((x(1)-x(2))/10.0_DP)+EXP((x(2)-x(3))/10.0_DP)))
                    ELSE IF(rowComponentIdx==3) THEN
                      !Calculate sum
                      sum=rowPhi*(-2.0_DP*muParam/100.0_DP*(2.0_DP*EXP((x(3)-x(1))/10.0_DP)+EXP((x(2)-x(3))/10.0_DP)))
                    ENDIF
                  ELSE IF(analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_3) THEN
                    IF(rowComponentIdx==1) THEN
                      !Calculate sum
                      sum=0.0_DP
                    ELSE IF(rowComponentIdx==2) THEN
                      !Calculate sum
                      sum=rowPhi*(36*muParam*PI**2/100.0_DP*COS(2.0_DP*PI*x(2)/10.0_DP)*SIN(2.0_DP*PI*x(3)/10.0_DP)* &
                        & COS(2.0_DP*PI*x(1)/10.0_DP))
                    ELSE IF(rowComponentIdx==3) THEN
                      !Calculate sum
                      sum=0.0_DP
                    ENDIF
                  ENDIF
                  !Calculate rhs vector
                  rhsVector%elementVector%vector(rowElementDOFIdx)=rhsVector%elementVector%vector(rowElementDOFIdx)+ &
                    & sum*jacobianGaussWeight
                ENDDO !rowElementParameterIdx
              ENDDO !rowComponentIdx
            ELSE
              rhsVector%elementVector%vector(rowElementDOFIdx)=0.0_DP
            ENDIF
          ENDIF
        ENDIF
      ENDDO !gaussPointIdx
    
      !Assemble matrices calculated above
      minRowElementDOFIdx=rowElementDOFIdx
      maxRowElementDOFIdx=columnElementDOFIdx
      minColumnElementDOFIdx=rowElementDOFIdx
      maxColumnElementDOFIdx=columnElementDOFIdx
      IF(updateStiffness) THEN
        stiffnessMatrix%elementMatrix%matrix(1:minRowElementDOFIdx,1:minColumnElementDOFIdx)= &
          & aMatrix(1:minRowElementDOFIdx,1:minColumnElementDOFIdx)+ &
          & dMatrix(1:minRowElementDOFIdx,1:minColumnElementDOFIdx)+ &
          & aleMatrix(1:minRowElementDOFIdx,1:minColumnElementDOFIdx)
        stiffnessMatrix%elementMatrix%matrix(1:minRowElementDOFIdx,minColumnElementDOFIdx+1:maxColumnElementDOFIdx)= &
          & bTmatrix(1:minRowElementDOFIdx,minColumnElementDOFIdx+1:maxColumnElementDOFIdx)
        DO rowElementDOFIdx=minRowElementDOFIdx+1,maxRowElementDOFIdx
          DO columnElementDOFIdx=1,minColumnElementDOFIdx
            !Transpose pressure type entries for mass equation
            stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
              & stiffnessMatrix%elementMatrix%matrix(columnElementDOFIdx,rowElementDOFIdx)
          ENDDO !columnElementDOFIdx
        ENDDO !rowElementDOFIdx
      ENDIF !updateStiffness
      
      !Scale factor adjustment
      CALL Field_ScalingTypeGet(dependentField,scalingType,err,error,*999)
      IF(scalingType/=FIELD_NO_SCALING) THEN
        NULLIFY(rowsInterpParameters)
        CALL EquationsInterpolation_DependentParametersGet(equationsInterpolation,rowsVariableType,rowsInterpParameters, &
          & err,error,*999)
        NULLIFY(colsInterpParameters)
        CALL EquationsInterpolation_DependentParametersGet(equationsInterpolation,colsVariableType,colsInterpParameters, &
          & err,error,*999)
        CALL Field_InterpolationParametersScaleFactorsElementGet(elementNumber,rowsInterpParameters,err,error,*999)
        CALL Field_InterpolationParametersScaleFactorsElementGet(elementNumber,colsInterpParameters,err,error,*999)
        !Loop over element rows
        rowElementDOFIdx=0          
        DO rowComponentIdx=1,numberOfRowsComponents
          NULLIFY(rowDomain)
          CALL FieldVariable_ComponentDomainGet(rowsVariable,rowComponentIdx,rowDomain,err,error,*999)
          NULLIFY(rowDomainTopology)
          CALL Domain_DomainTopologyGet(rowDomain,rowDomainTopology,err,error,*999)
          NULLIFY(rowDomainElements)
          CALL DomainTopology_DomainElementsGet(rowDomainTopology,rowDomainElements,err,error,*999)
          NULLIFY(rowBasis)
          CALL DomainElements_ElementBasisGet(rowDomainElements,elementNumber,rowBasis,err,error,*999)
          CALL Basis_NumberOfElementParametersGet(rowBasis,numberOfRowElementParameters,err,error,*999)
          DO rowElementParameterIdx=1,numberOfRowElementParameters
            rowElementDOFIdx=rowElementDOFIdx+1                    
            IF(updateMatrices) THEN
              !Loop over element columns
              columnElementDOFIdx=0
              DO columnComponentIdx=1,numberOfColsComponents
                NULLIFY(columnDomain)
                CALL FieldVariable_ComponentDomainGet(colsVariable,columnComponentIdx,columnDomain,err,error,*999)
                NULLIFY(columnDomainTopology)
                CALL Domain_DomainTopologyGet(columnDomain,columnDomainTopology,err,error,*999)
                NULLIFY(columnDomainElements)
                CALL DomainTopology_DomainElementsGet(columnDomainTopology,columnDomainElements,err,error,*999)
                NULLIFY(columnBasis)
                CALL DomainElements_ElementBasisGet(columnDomainElements,elementNumber,columnBasis,err,error,*999)
                CALL Basis_NumberOfElementParametersGet(columnBasis,numberOfColumnElementParameters,err,error,*999)
                DO columnElementParameterIdx=1,numberOfColumnElementParameters
                  columnElementDOFIdx=columnElementDOFIdx+1
                  IF(updateStiffness) THEN
                    stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                      & stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)* &
                      & rowsInterpParameters%scaleFactors(rowElementParameterIdx,rowComponentIdx)* &
                      & colsInterpParameters%scaleFactors(columnElementParameterIdx,columnComponentIdx)
                  ENDIF
                  IF(updateDamping) THEN
                    dampingMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                      & dampingMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)* &
                      & rowsInterpParameters%scaleFactors(rowElementParameterIdx,rowComponentIdx)* &
                      & colsInterpParameters%scaleFactors(columnElementParameterIdx,columnComponentIdx)
                  ENDIF
                ENDDO !columnElementParameterIdx
              ENDDO !columnComponentIdx
            ENDIF !update matrices
            IF(updateSource) THEN
              sourceVector%elementVector%vector(rowElementDOFIdx)=sourceVector%elementVector%vector(rowElementDOFIdx)* &
                & rowsInterpParameters%scaleFactors(rowElementParameterIdx,rowComponentIdx)
            ENDIF
            IF(updateRHS) THEN
              rhsVector%elementVector%vector(rowElementDOFIdx)=rhsVector%elementVector%vector(rowElementDOFIdx)* &
                & rowsInterpParameters%scaleFactors(rowElementParameterIdx,rowComponentIdx)
            ENDIF
          ENDDO !rowElementParameterIdx
        ENDDO !rowComponentIdx
      ENDIF !scaling
    ENDIF !update

    EXITS("Stokes_FiniteElementCalculate")
    RETURN
999 ERRORSEXITS("Stokes_FiniteElementCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE Stokes_FiniteElementCalculate

  !
  !================================================================================================================================
  !

  !>Sets up the Stokes problem post solve.
  SUBROUTINE Stokes_PostSolve(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver!<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: outputType,pSpecification(3),solveType
    TYPE(ControlLoopType), POINTER :: controlLoop 
    TYPE(DynamicSolverType), POINTER :: dynamicSolver
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverType), POINTER :: solver2
    TYPE(SolversType), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError

    ENTERS("Stokes_PostSolve",err,error,*999)
    
    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)

    CALL Solver_OutputTypeGet(solver,outputType,err,error,*999)

    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STATIC_STOKES_SUBTYPE,PROBLEM_LAPLACE_STOKES_SUBTYPE)
      CALL Stokes_PostSolveOutputData(solver,err,error,*999)
    CASE(PROBLEM_TRANSIENT_STOKES_SUBTYPE)
      CALL Stokes_PostSolveOutputData(solver,err,error,*999)
    CASE(PROBLEM_ALE_STOKES_SUBTYPE)
      CALL Solver_TypeGet(solver,solveType,err,error,*999)
      IF(solveType==SOLVER_LINEAR_TYPE) THEN
        !Post solve for the linear solver
        IF(outputType>=SOLVER_PROGRESS_OUTPUT) &
          & CALL WriteString(GENERAL_OUTPUT_TYPE,"Mesh movement post solve... ",err,error,*999)
        NULLIFY(solvers)
        CALL Solver_SolversGet(solver,solvers,err,error,*999)
        NULLIFY(solver2)
        CALL Solvers_SolverGet(solvers,2,solver2,err,error,*999)
        NULLIFY(dynamicSolver)
        CALL Solver_DynamicSolverGet(solver,dynamicSolver,err,error,*999)
        dynamicSolver%ale=.TRUE.
      ELSE IF(solveType==SOLVER_DYNAMIC_TYPE) THEN
        !Post solve for the dynamic solver
        IF(outputType>=SOLVER_PROGRESS_OUTPUT) &
          & CALL WriteString(GENERAL_OUTPUT_TYPE,"ALE Stokes post solve... ",err,error,*999)
        CALL Stokes_PostSolveOutputData(solver,err,error,*999)
      END IF
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a Stokes fluid type of a fluid mechanics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Stokes_PostSolve")
    RETURN
999 ERRORSEXITS("Stokes_PostSolve",err,error)
    RETURN 1
    
  END SUBROUTINE Stokes_PostSolve

  !
  !================================================================================================================================
  !

  !>Sets up the Stokes problem pre solve.
  SUBROUTINE Stokes_PreSolve(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: outputType,pSpecification(3),solveType
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(DynamicSolverType), POINTER :: dynamicSolver
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverType), POINTER :: solver2 
    TYPE(SolversType), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError

    ENTERS("Stokes_PreSolve",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)

    CALL Solver_OutputTypeGet(solver,outputType,err,error,*999)
    
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STATIC_STOKES_SUBTYPE,PROBLEM_LAPLACE_STOKES_SUBTYPE)
      ! do nothing ???
    CASE(PROBLEM_TRANSIENT_STOKES_SUBTYPE)
      CALL Stokes_PreSolveUpdateBoundaryConditions(solver,err,error,*999)
    CASE(PROBLEM_ALE_STOKES_SUBTYPE)
      CALL Solver_TypeGet(solver,solveType,err,error,*999)
      !Pre solve for the linear solver
      IF(solveType==SOLVER_LINEAR_TYPE) THEN
        IF(outputType>=SOLVER_PROGRESS_OUTPUT) &
          & CALL WriteString(GENERAL_OUTPUT_TYPE,"Mesh movement pre solve... ",err,error,*999)
        !Update boundary conditions for mesh-movement
        CALL Stokes_PreSolveUpdateBoundaryConditions(solver,err,error,*999)
        NULLIFY(solvers)
        CALL Solver_SolversGet(solver,solvers,err,error,*999)
        NULLIFY(solver2)
        CALL Solvers_SolverGet(solvers,2,solver2,err,error,*999)
        NULLIFY(dynamicSolver)
        CALL Solver_DynamicSolverGet(solver2,dynamicSolver,err,error,*999)
        !\todo: Avoid ALE flag in future
        dynamicSolver%ale=.FALSE.
        !Update material properties for Laplace mesh movement
        CALL Stokes_PreSolveALEUpdateParameters(solver,err,error,*999)
      ELSE IF(solveType==SOLVER_DYNAMIC_TYPE) THEN
        !Pre solve for the dynamic solver
        IF(outputType>=SOLVER_PROGRESS_OUTPUT) &
          & CALL WriteString(GENERAL_OUTPUT_TYPE,"ALE Stokes pre solve... ",err,error,*999)
        NULLIFY(dynamicSolver)
        CALL Solver_DynamicSolverGet(solver2,dynamicSolver,err,error,*999)
        IF(dynamicSolver%ale) THEN
          !First update mesh and calculates boundary velocity values
          CALL Stokes_PreSolveALEUpdateMesh(solver,err,error,*999)
          !Then apply both normal and moving mesh boundary conditions
          CALL Stokes_PreSolveUpdateBoundaryConditions(solver,err,error,*999)
        ELSE
          CALL FlagError("Mesh motion calculation not successful for ALE problem.",err,error,*999)
        END IF
      ELSE
        CALL FlagError("Solver type is not associated for ALE problem.",err,error,*999)
      END IF
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a Stokes fluid type of a fluid mechanics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
 
    EXITS("Stokes_PreSolve")
    RETURN
999 ERRORSEXITS("Stokes_PreSolve",err,error)
    RETURN 1
    
  END SUBROUTINE Stokes_PreSolve
  !
  !================================================================================================================================
  !

  !>Update boundary conditions for Stokes flow pre solve
  SUBROUTINE Stokes_PreSolveUpdateBoundaryConditions(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: analyticFunctionType,boundaryConditionCheckVariable,componentIdx,currentIteration,derivativeIdx,dimensionIdx, &
      & elementIdx,globalDerivativeIndex,i,inputIteration,j,k,localDOFIdx,localNodeIdx,localNodeNumber, &
      & maximumNumberOfElementParameters,nodeIdx,numberOfComponents,numberOfDimensions,numberOfNodes,numberOfNodeDerivatives, &
      & numberOfNodesXic(3),numberOfVariables,outputIteration,outputType,pSpecification(3),solveType,variableIdx,variableType
    REAL(DP) :: analyticValue,currentTime,displacementValue,muParam,rhoParam,startTime,stopTime,tCoordinates(20,3), &
      & timeIncrement,x(3),xiCoordinates(4)
    REAL(DP), POINTER :: boundaryValues(:),geometricParameters(:),materialsParameters(:),meshVelocityValues(:)
    TYPE(BasisType), POINTER :: basis
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainElementsType), POINTER :: domainElements
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(EquationsSetAnalyticType), POINTER :: equationsAnalytic
    TYPE(FieldType), POINTER :: dependentField,geometricField,independentField,materialsField
    TYPE(FieldInterpolatedPointType), POINTER :: interpolatedPoint
    TYPE(FieldInterpolationParametersType), POINTER :: interpolationParameters
    TYPE(FieldVariableType), POINTER :: dependentVariable,uDependentVariable,geometricVariable,independentVariable, &
      & materialsVariable
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(VARYING_STRING) :: localError

    ENTERS("Stokes_PreSolveUpdateBoundaryConditions",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)

    CALL Solver_TypeGet(solver,solveType,err,error,*999)
    CALL Solver_OutputTypeGet(solver,outputType,err,error,*999)
    CALL ControlLoop_CurrentTimeInformationGet(controlLoop,currentTime,timeIncrement,startTime,stopTime,currentIteration, &
      & outputIteration,inputIteration,err,error,*999)

    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STATIC_STOKES_SUBTYPE,PROBLEM_LAPLACE_STOKES_SUBTYPE)
      ! do nothing ???
    CASE(PROBLEM_TRANSIENT_STOKES_SUBTYPE)
      NULLIFY(solverEquations)
      CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      NULLIFY(boundaryConditions)
      CALL SolverEquations_BoundaryConditionsGet(solverEquations,boundaryConditions,err,error,*999)
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
      NULLIFY(equationsAnalytic)
      CALL EquationsSet_AnalyticExists(equationsSet,equationsAnalytic,err,error,*999)
      IF(ASSOCIATED(equationsAnalytic)) THEN
        CALL EquationsSet_AnalyticFunctionTypeGet(equationsSet,analyticFunctionType,err,error,*999)
        IF(analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_4.OR. &
          & analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_5.OR. &
          & analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_1.OR. &
          & analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_4.OR. &
          & analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_5) THEN
          NULLIFY(geometricField)
          CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
          NULLIFY(geometricVariable)
          CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,geometricVariable,err,error,*999)
          CALL FieldVariable_NumberOfComponentsGet(geometricVariable,numberOfDimensions,err,error,*999)
          NULLIFY(geometricParameters)
          CALL FieldVariable_ParameterSetDataGet(geometricVariable,FIELD_VALUES_SET_TYPE,geometricParameters,err,error,*999)
          NULLIFY(interpolationParameters)
          CALL FieldVariable_InterpolationParametersInitialise(geometricVariable,interpolationParameters,err,error,*999)
          NULLIFY(interpolatedPoint)
          CALL Field_InterpolatedPointInitialise(interpolationParameters,interpolatedPoint,err,error,*999)
          NULLIFY(materialsField)
          CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)
          NULLIFY(materialsVariable)
          CALL Field_VariableGet(materialsField,FIELD_U_VARIABLE_TYPE,materialsVariable,err,error,*999)
          NULLIFY(materialsParameters)
          CALL FieldVariable_ParameterSetDataGet(materialsVariable,FIELD_VALUES_SET_TYPE,materialsParameters,err,error,*999)
          !Define muParam, density=1
          muParam=materialsParameters(1)
          !Define rhoParam, density=2
          IF(analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_4.OR. &
            & analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_5.OR. &
            & analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_4.OR. &
            & analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_5) THEN
            rhoParam=materialsParameters(2)
          ELSE
            rhoParam=0.0_DP
          ENDIF
          NULLIFY(dependentField)
          CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
          NULLIFY(uDependentVariable)
          CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,uDependentVariable,err,error,*999)
          CALL Field_NumberOfVariablesGet(dependentField,numberOfVariables,err,error,*999)
          DO variableIdx=1,numberOfVariables
            NULLIFY(dependentVariable)
            CALL Field_VariableIndexGet(dependentField,variableIdx,dependentVariable,variableType,err,error,*999)
            CALL FieldVariable_NumberOfComponentsGet(dependentVariable,numberOfComponents,err,error,*999)
            DO componentIdx=1,numberOfComponents
              CALL FieldVariable_ComponentInterpolationCheck(dependentVariable,componentIdx,FIELD_NODE_BASED_INTERPOLATION, &
                & err,error,*999)
              NULLIFY(domain)
              CALL FieldVariable_ComponentDomainGet(dependentVariable,componentIdx,domain,err,error,*999)
              NULLIFY(domainTopology)
              CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
              NULLIFY(domainNodes)
              CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
              NULLIFY(domainElements)
              CALL DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*999)
              !Should be replaced by boundary node flag
              CALL DomainNodes_NumberOfNodesGet(domainNodes,numberOfNodes,err,error,*999)
              DO nodeIdx=1,numberOfNodes
                CALL DomainNodes_NodeSurroundingElementGet(domainNodes,1,nodeIdx,elementIdx,err,error,*999)
                NULLIFY(basis)
                CALL DomainElements_ElementBasisGet(domainElements,elementIdx,basis,err,error,*999)
                CALL DomainElements_MaxElementParametersGet(domainElements,maximumNumberOfElementParameters,err,error,*999)
                CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementIdx,interpolationParameters, &
                  & err,error,*999)
                localNodeIdx=0
                xiCoordinates=0.0_DP
                numberOfNodesXic=1
                CALL Basis_NumberOfNodesXiCGet(basis,numberOfNodesXiC,err,error,*999)
                !\todo: change definitions as soon as adjacent elements / boundary elements calculation works for simplex
                IF(maximumNumberOfElementParameters==4.OR. &
                  & maximumNumberOfElementParameters==9.OR. &
                  & maximumNumberOfElementParameters==16.OR. &
                  & maximumNumberOfElementParameters==8.OR. &
                  & maximumNumberOfElementParameters==27.OR. &
                  & maximumNumberOfElementParameters==64) THEN
                  DO k=1,numberOfNodesXic(3)
                    DO j=1,numberOfNodesXic(2)
                      DO i=1,numberOfNodesXic(1)
                        localNodeIdx=localNodeIdx+1
                        CALL DomainElements_ElementNodeGet(domainElements,localNodeIdx,elementIdx,localNodeNumber,err,error,*999)
                        IF(localNodeNumber==nodeIdx) EXIT
                        xiCoordinates(1)=xiCoordinates(1)+(1.0_DP/(numberOfNodesXic(1)-1))
                      ENDDO !i
                      IF(localNodeNumber==nodeIdx) EXIT
                      xiCoordinates(1)=0.0_DP
                      xiCoordinates(2)=xiCoordinates(2)+(1.0_DP/(numberOfNodesXic(2)-1))
                    ENDDO !j
                    IF(localNodeNumber==nodeIdx) EXIT
                    xiCoordinates(1)=0.0_DP
                    xiCoordinates(2)=0.0_DP
                    IF(numberOfNodesXic(3)/=1) xiCoordinates(3)=xiCoordinates(3)+(1.0_DP/(numberOfNodesXic(3)-1))
                  ENDDO !k
                  CALL Field_InterpolateXi(NO_PART_DERIV,xiCoordinates,interpolatedPoint,err,error,*999)
                ELSE
                  !\todo: Use boundary flag
                  IF(maximumNumberOfElementParameters==3) THEN
                    tCoordinates(1,1:2)=[0.0_DP,1.0_DP]
                    tCoordinates(2,1:2)=[1.0_DP,0.0_DP]
                    tCoordinates(3,1:2)=[1.0_DP,1.0_DP]
                  ELSE IF(maximumNumberOfElementParameters==6) THEN
                    tCoordinates(1,1:2)=[0.0_DP,1.0_DP]
                    tCoordinates(2,1:2)=[1.0_DP,0.0_DP]
                    tCoordinates(3,1:2)=[1.0_DP,1.0_DP]
                    tCoordinates(4,1:2)=[0.5_DP,0.5_DP]
                    tCoordinates(5,1:2)=[1.0_DP,0.5_DP]
                    tCoordinates(6,1:2)=[0.5_DP,1.0_DP]
                  ELSE IF(maximumNumberOfElementParameters==10.AND.numberOfDimensions==2) THEN
                    tCoordinates(1,1:2)=[0.0_DP,1.0_DP]
                    tCoordinates(2,1:2)=[1.0_DP,0.0_DP]
                    tCoordinates(3,1:2)=[1.0_DP,1.0_DP]
                    tCoordinates(4,1:2)=[1.0_DP/3.0_DP,2.0_DP/3.0_DP]
                    tCoordinates(5,1:2)=[2.0_DP/3.0_DP,1.0_DP/3.0_DP]
                    tCoordinates(6,1:2)=[1.0_DP,1.0_DP/3.0_DP]
                    tCoordinates(7,1:2)=[1.0_DP,2.0_DP/3.0_DP]
                    tCoordinates(8,1:2)=[2.0_DP/3.0_DP,1.0_DP]
                    tCoordinates(9,1:2)=[1.0_DP/3.0_DP,1.0_DP]
                    tCoordinates(10,1:2)=[2.0_DP/3.0_DP,2.0_DP/3.0_DP]
                  ELSE IF(maximumNumberOfElementParameters==4) THEN
                    tCoordinates(1,1:3)=[0.0_DP,1.0_DP,1.0_DP]
                    tCoordinates(2,1:3)=[1.0_DP,0.0_DP,1.0_DP]
                    tCoordinates(3,1:3)=[1.0_DP,1.0_DP,0.0_DP]
                    tCoordinates(4,1:3)=[1.0_DP,1.0_DP,1.0_DP]
                  ELSE IF(maximumNumberOfElementParameters==10.AND.numberOfDimensions==3) THEN
                    tCoordinates(1,1:3)=[0.0_DP,1.0_DP,1.0_DP]
                    tCoordinates(2,1:3)=[1.0_DP,0.0_DP,1.0_DP]
                    tCoordinates(3,1:3)=[1.0_DP,1.0_DP,0.0_DP]
                    tCoordinates(4,1:3)=[1.0_DP,1.0_DP,1.0_DP]
                    tCoordinates(5,1:3)=[0.5_DP,0.5_DP,1.0_DP]
                    tCoordinates(6,1:3)=[0.5_DP,1.0_DP,0.5_DP]
                    tCoordinates(7,1:3)=[0.5_DP,1.0_DP,1.0_DP]
                    tCoordinates(8,1:3)=[1.0_DP,0.5_DP,0.5_DP]
                    tCoordinates(9,1:3)=[1.0_DP,1.0_DP,0.5_DP]
                    tCoordinates(10,1:3)=[1.0_DP,0.5_DP,1.0_DP]
                  ELSE IF(maximumNumberOfElementParameters==20) THEN
                    tCoordinates(1,1:3)=[0.0_DP,1.0_DP,1.0_DP]
                    tCoordinates(2,1:3)=[1.0_DP,0.0_DP,1.0_DP]
                    tCoordinates(3,1:3)=[1.0_DP,1.0_DP,0.0_DP]
                    tCoordinates(4,1:3)=[1.0_DP,1.0_DP,1.0_DP]
                    tCoordinates(5,1:3)=[1.0_DP/3.0_DP,2.0_DP/3.0_DP,1.0_DP]
                    tCoordinates(6,1:3)=[2.0_DP/3.0_DP,1.0_DP/3.0_DP,1.0_DP]
                    tCoordinates(7,1:3)=[1.0_DP/3.0_DP,1.0_DP,2.0_DP/3.0_DP]
                    tCoordinates(8,1:3)=[2.0_DP/3.0_DP,1.0_DP,1.0_DP/3.0_DP]
                    tCoordinates(9,1:3)=[1.0_DP/3.0_DP,1.0_DP,1.0_DP]
                    tCoordinates(10,1:3)=[2.0_DP/3.0_DP,1.0_DP,1.0_DP]
                    tCoordinates(11,1:3)=[1.0_DP,1.0_DP/3.0_DP,2.0_DP/3.0_DP]
                    tCoordinates(12,1:3)=[1.0_DP,2.0_DP/3.0_DP,1.0_DP/3.0_DP]
                    tCoordinates(13,1:3)=[1.0_DP,1.0_DP,1.0_DP/3.0_DP]
                    tCoordinates(14,1:3)=[1.0_DP,1.0_DP,2.0_DP/3.0_DP]
                    tCoordinates(15,1:3)=[1.0_DP,1.0_DP/3.0_DP,1.0_DP]
                    tCoordinates(16,1:3)=[1.0_DP,2.0_DP/3.0_DP,1.0_DP]
                    tCoordinates(17,1:3)=[2.0_DP/3.0_DP,2.0_DP/3.0_DP,2.0_DP/3.0_DP]
                    tCoordinates(18,1:3)=[2.0_DP/3.0_DP,2.0_DP/3.0_DP,1.0_DP]
                    tCoordinates(19,1:3)=[2.0_DP/3.0_DP,1.0_DP,2.0_DP/3.0_DP]
                    tCoordinates(20,1:3)=[1.0_DP,2.0_DP/3.0_DP,2.0_DP/3.0_DP]
                  ENDIF
                  DO k=1,maximumNumberOfElementParameters
                    CALL DomainElements_ElementNodeGet(domainElements,k,elementIdx,localNodeNumber,err,error,*999)
                    IF(localNodeNumber==nodeIdx) EXIT
                  ENDDO !k
                  CALL Field_InterpolateXi(NO_PART_DERIV,tCoordinates(K,1:numberOfDimensions),interpolatedPoint,err,error,*999)
                ENDIF
                x=0.0_DP
                DO dimensionIdx=1,numberOfDimensions
                  x(dimensionIdx)=interpolatedPoint%values(dimensionIdx,NO_PART_DERIV)
                ENDDO !dimensionIdx
                NULLIFY(boundaryConditionsVariable)
                CALL BoundaryConditions_VariableGet(boundaryConditions,dependentVariable,boundaryConditionsVariable,err,error,*999)
                !Loop over the derivatives
                CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,nodeIdx,numberOfNodeDerivatives,err,error,*999)
                DO derivativeIdx=1,numberOfNodeDerivatives
                  CALL DomainNodes_DerivativeGlobalIndexGet(domainNodes,derivativeIdx,nodeIdx,globalDerivativeIndex,err,error,*999)
                  CALL Stokes_AnalyticFunctions(analyticValue,x,muParam,rhoParam,currentTime,variableType, &
                    & globalDerivativeIndex,analyticFunctionType,numberOfDimensions,numberOfComponents,componentIdx, &
                    & err,error,*999)
                  !Default to version 1 of each node derivative
                  CALL FieldVariable_LocalNodeDOFGet(dependentVariable,1,derivativeIdx,nodeIdx,componentIdx,localDOFIdx, &
                    & err,error,*999)
                  CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,localDOFIdx, &
                    & analyticValue,err,error,*999)
                  CALL BoundaryConditionsVariable_ConditionTypeGet(boundaryConditionsVariable,localDOFIdx, &
                    & boundaryConditionCheckVariable,err,error,*999)
                  IF(boundaryConditionCheckVariable==BOUNDARY_CONDITION_FIXED) THEN
                    CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_VALUES_SET_TYPE,localDOFIdx, &
                      & analyticValue,err,error,*999)
                  ENDIF
                ENDDO !derivativeIdx
              ENDDO !nodeIdx
            ENDDO !componentIdx
            CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
            CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
            CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
            CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
          ENDDO !variableIdx
          CALL FieldVariable_ParameterSetDataRestore(geometricVariable,FIELD_VALUES_SET_TYPE,geometricParameters,err,error,*999)
          CALL FieldVariable_ParameterSetDataRestore(materialsVariable,FIELD_VALUES_SET_TYPE,materialsParameters,err,error,*999)
        ENDIF
      ENDIF
      CALL FieldVariable_ParameterSetUpdateStart(uDependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
      CALL FieldVariable_ParameterSetUpdateFinish(uDependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
    CASE(PROBLEM_ALE_STOKES_SUBTYPE)
      !Pre solve for the linear solver
      IF(solveType==SOLVER_LINEAR_TYPE) THEN
        IF(outputType>=SOLVER_PROGRESS_OUTPUT) &
          & CALL WriteString(GENERAL_OUTPUT_TYPE,"Mesh movement change boundary conditions... ",err,error,*999)
        NULLIFY(solverEquations)
        CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
        NULLIFY(solverMapping)
        CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
        NULLIFY(boundaryConditions)
        CALL SolverEquations_BoundaryConditionsGet(solverEquations,boundaryConditions,err,error,*999)
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverEquations,1,equationsSet,err,error,*999)
        NULLIFY(geometricField)
        CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
        NULLIFY(geometricVariable)
        CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,geometricVariable,err,error,*999)
        CALL FieldVariable_NumberOfComponentsGet(geometricVariable,numberOfDimensions,err,error,*999)
        NULLIFY(dependentField)
        CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
        NULLIFY(dependentVariable)
        CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,dependentVariable,err,error,*999)
        NULLIFY(boundaryConditionsVariable)
        CALL BoundaryConditions_VariableGet(boundaryConditions,dependentVariable,boundaryConditionsVariable,err,error,*999)
        NULLIFY(independentField)
        CALL EquationsSet_IndependentFieldGet(equationsSet,independentField,err,error,*999)
        NULLIFY(independentVariable)
        CALL Field_VariableGet(independentField,FIELD_U_VARIABLE_TYPE,independentVariable,err,error,*999)
        NULLIFY(boundaryValues)
        CALL FieldVariable_ParameterSetDataGet(independentVariable,FIELD_BOUNDARY_SET_TYPE,boundaryValues,err,error,*999)
        CALL FLUID_MECHANICS_IO_READ_boundaryConditions(SOLVER_LINEAR_TYPE,boundaryValues,numberOfDimensions, &
          & BOUNDARY_CONDITION_MOVED_WALL,inputIteration,currentIteration,,currentTime,1.0_DP,err,error,*999)
        CALL Field_NumberOfVariablesGet(dependentField,numberOfVariables,err,error,*999)
        DO variableIdx=1,numberOfVariables
          NULLIFY(dependentVariable)
          CALL Field_VariableIndexGet(dependentField,variableIdx,dependentVariable,variableType,err,error,*999)
          CALL FieldVariable_NumberOfComponentsGet(dependentVariable,numberOfComponents,err,error,*999)
          DO componentIdx=1,numberOfComponents
            NULLIFY(domain)
            CALL FieldVariable_ComponentDomainGet(dependentVariable,componentIdx,domain,err,error,*999)
            NULLIFY(domainTopology)
            CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
            NULLIFY(domainNodes)
            CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
            !Loop over the local nodes excluding the ghosts.
            CALL DomainNodes_NumberOfNodesGet(domainNodes,numberOfNodes,err,error,*999)
            DO nodeIdx=1,numberOfNodes
              CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,nodeIdx,numberOfNodeDerivatives,err,error,*999)
              DO derivativeIdx=1,numberOfNodeDerivatives
                !Default to version 1 of each node derivative
                CALL FieldVariable_LocalNodeDOFGet(dependentVariable,1,derivativeIdx,nodeIdx,componentIdx,localDOFIdx, &
                  & err,error,*999)
                CALL BoundaryConditionsVariable_ConditionTypeGet(boundaryConditionsVariable,localDOFIdx, &
                  & boundaryConditionCheckVariable,err,error,*999)
                IF(boundaryConditionCheckVariable==BOUNDARY_CONDITION_MOVED_WALL) THEN
                  CALL Field_ParameterSetUpdateLocalDOF(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                    & localDOFIdx,boundaryValues(localDOFIdx),err,error,*999)
                ENDIF
              ENDDO !derivativeIdx
            ENDDO !nodeIdx
          ENDDO !componentIdx
        ENDDO !variableIdx
        CALL FieldVariable_ParameterSetDataRestore(independentVariable,FIELD_BOUNDARY_SET_TYPE,boundaryValues,err,error,*999)
        !\todo: This part should be read in out of a file eventually
        CALL Field_ParameterSetUpdateStart(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
        CALL Field_ParameterSetUpdateFinish(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
        !Pre solve for the dynamic solver
      ELSE IF(solveType==SOLVER_DYNAMIC_TYPE) THEN
        IF(outputType>=SOLVER_PROGRESS_OUTPUT) &
          & CALL WriteString(GENERAL_OUTPUT_TYPE,"Mesh movement change boundary conditions... ",err,error,*999)
        NULLIFY(solverEquations)
        CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
        NULLIFY(solverMapping)
        CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
        NULLIFY(boundaryConditions)
        CALL SolverEquations_BoundaryConditionsGet(solverEquations,boundaryConditions,err,error,*999)
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverEquations,1,equationsSet,err,error,*999)
        NULLIFY(geometricField)
        CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
        NULLIFY(geometricVariable)
        CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,geometricVariable,err,error,*999)
        CALL FieldVariable_NumberOfComponentsGet(geometricVariable,numberOfDimensions,err,error,*999)
        NULLIFY(dependentField)
        CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
        NULLIFY(dependentVariable)
        CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,dependentVariable,err,error,*999)
        NULLIFY(boundaryConditionsVariable)
        CALL BoundaryConditions_VariableGet(boundaryConditions,dependentVariable,boundaryConditionsVariable,err,error,*999)
        NULLIFY(independentField)
        CALL EquationsSet_IndependentFieldGet(equationsSet,independentField,err,error,*999)
        NULLIFY(independentVariable)
        CALL Field_VariableGet(independentField,FIELD_U_VARIABLE_TYPE,independentVariable,err,error,*999)
        NULLIFY(meshVelocityValues)
        CALL FieldVariable_ParameterSetDataGet(independentVariable,FIELD_MESH_VELOCITY_SET_TYPE,meshVelocityValues,err,error,*999)
        NULLIFY(boundaryValues)
        CALL FieldVariable_ParameterSetDataGet(independentVariable,FIELD_BOUNDARY_SET_TYPE,boundaryValues,err,error,*999)
        CALL FLUID_MECHANICS_IO_READ_boundaryConditions(SOLVER_LINEAR_TYPE,boundaryValues,numberOfDimensions, &
          & BOUNDARY_CONDITION_FIXED_INLET,inputIteration,currentIteration,currentTime,1.0_DP,err,error,*999)
        CALL Field_NumberOfVariablesGet(dependentField,numberOfVariables,err,error,*999)
        DO variableIdx=1,numberOfVariables
          NULLIFY(dependentVariable)
          CALL FieldVariable_NumberOfComponentsGet(dependentVariable,numberOfComponents,err,error,*999)
          DO componentIdx=1,numberOfComponents
            NULLIFY(domain)
            CALL FieldVariable_ComponentDomainGet(dependentVariable,componentIdx,domain,err,error,*999)
            NULLIFY(domainTopology)
            CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
            NULLIFY(domainNodes)
            CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
            !Loop over the local nodes excluding the ghosts.
            CALL DomainNodes_NumberOfNodesGet(domainNodes,numberOfNodes,err,error,*999)
            DO nodeIdx=1,numberOfNodes
              CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,nodeIdx,numberOfNodeDerivatives,err,error,*999)
              DO derivativeIdx=1,numberOfNodeDerivatives
                !Default to version 1 of each node derivative
                CALL FieldVariable_LocalNodeDOFGet(dependentVariable,1,derivativeIdx,nodeIdx,componentIdx,localDOFIdx, &
                  & err,error,*999)
                displacementValue=0.0_DP
                CALL BoundaryConditionsVariable_ConditionTypeGet(boundaryConditionsVariable,localDOFIdx, &
                  & boundaryConditionCheckVariable,err,error,*999)
                IF(boundaryConditionCheckVariable==BOUNDARY_CONDITION_MOVED_WALL) THEN
                  CALL Field_ParameterSetUpdateLocalDOF(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,localDOFIdx, &
                    & meshVelocityValues(localDOFIdx),err,error,*999)
                ELSE IF(boundaryConditionCheckVariable==BOUNDARY_CONDITION_FIXED_INLET) THEN
                  CALL Field_ParameterSetUpdateLocalDOF(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,localDOFIdx, &
                    & boundaryValues(localDOFIdx),err,error,*999)
                ENDIF
              ENDDO !derivativeIdx
            ENDDO !nodeIdx
          ENDDO !componentIdx
        ENDDO !variableIdx
        CALL FieldVariable_ParameterSetDataRestore(independentVariable,FIELD_MESH_VELOCITY_SET_TYPE,meshVelocityValues, &
          & err,error,*999)
        CALL FieldVariable_ParameterSetDataRestore(independentVariable,FIELD_BOUNDARY_SET_TYPE,boundaryValues,err,error,*999)
        CALL Field_ParameterSetUpdateStart(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
        CALL Field_ParameterSetUpdateFinish(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
      END IF
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a Stokes equation fluid type of a fluid mechanics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
     
    EXITS("Stokes_PreSolveUpdateBoundaryConditions")
    RETURN
999 ERRORSEXITS("Stokes_PreSolveUpdateBoundaryConditions",err,error)
    RETURN 1
    
  END SUBROUTINE Stokes_PreSolveUpdateBoundaryConditions

  !
  !================================================================================================================================
  !
  
  !>Update mesh velocity and move mesh for ALE Stokes problem
  SUBROUTINE Stokes_PreSolveALEUpdateMesh(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx,currentIteration,derivativeIdx,geometricMeshComponent,inputIteration,inputOption,inputType, &
      & localDOFIdx,numberOfComponents,numberOfDimensionsLaplace,numberOfDimensionsALEStokes,numberOfNodes, &
      & numberOfNodeDerivatives,numberOfVariables,nodeIdx,outputIteration,pSpecification(3),solveType,variableIdx,variableType
    REAL(DP) :: alpha,currentTime,startTime,stopTime,timeIncrement
    REAL(DP), POINTER :: meshDisplacementValues(:)
    TYPE(ControlLoopType), POINTER :: controlLoop 
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(DynamicSolverType), POINTER :: dynamicSolver
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsSetType), POINTER :: equationsSetLaplace,equationsSetALEStokes 
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentFieldALEStokes,dependentFieldLaplace,geometricFieldALEStokes,geometricFieldLaplace, &
      & independentFieldALEStokes
    TYPE(FieldVariableType), POINTER :: dependentVariableALEStokes,dependentVariableLaplace,geometricVariableALEStokes, &
      & geometricVariableLaplace,independentVariableALEStokes,independentVariableLaplace
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverType), POINTER :: solverALEStokes,solverLaplace
    TYPE(SolverEquationsType), POINTER :: solverEquationsLaplace,solverEquationsALEStokes 
    TYPE(SolverMappingType), POINTER :: solverMappingLaplace,solverMappingALEStokes
    TYPE(SolversType), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError

    ENTERS("Stokes_PreSolveALEUpdateMesh",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)

    CALL Solver_TypeGet(solver,solveType,err,error,*999)
    CALL ControlLoop_CurrentTimeInformationGet(controlLoop,currentTime,timeIncrement,startTime,stopTime,currentIteration, &
      & outputIteration,inputIteration,err,error,*999)

    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STATIC_STOKES_SUBTYPE,PROBLEM_LAPLACE_STOKES_SUBTYPE)
      ! do nothing ???
    CASE(PROBLEM_TRANSIENT_STOKES_SUBTYPE)
      ! do nothing ???
    CASE(PROBLEM_ALE_STOKES_SUBTYPE)
      !Update mesh within the dynamic solver
      IF(solveType/=SOLVER_DYNAMIC_TYPE) CALL FlagError("Mesh update is not defined for non-dynamic problems.",err,error,*999)
      NULLIFY(dynamicSolver)
      CALL Solver_DynamicSolverGet(solver,dynamicSolver,err,error,*999)
      IF(.NOT.dynamicSolver%ale) CALL FlagError("Mesh motion calculation not successful for ALE problem.",err,error,*999)
      !Get the dependent field for the three component Laplace problem
      NULLIFY(solvers)
      CALL Solver_SolversGet(solver,solvers,err,error,*999)
      NULLIFY(solverLaplace)
      CALL Solvers_SolverGet(solvers,1,solverLaplace,err,error,*999)
      NULLIFY(solverEquationsLaplace)
      CALL Solver_SolverEquations(solverLaplace,solverEquationsLaplace,err,error,*999)
      NULLIFY(solverMappingLaplace)
      CALL SolverEquations_SolverMappingGet(solverEquationsLaplace,solverMappingLaplace,err,error,*999)
      NULLIFY(equationsSetLaplace)
      CALL SolverMapping_EquationsSetGet(solverMappingLaplace,equationsSetLaplace,err,error,*999)
      NULLIFY(geometricFieldLaplace)
      CALL EquationsSet_GeometricFieldGet(equationsSetLaplace,geometricFieldLaplace,err,error,*999)
      NULLIFY(geometricVariableLaplace)
      CALL Field_VariableGet(geometricFieldLaplace,FIELD_U_VARIABLE_TYPE,geometricVariableLaplace,err,error,*999)
      CALL FieldVarible_NumberOfComponentsGet(geometricVariableLaplace,numberOfDimensionsLaplace,err,error,*999)
      NULLIFY(dependentFieldLaplace)
      CALL EquationsSet_DependentFieldGet(equationsSetLaplace,dependentFieldLaplace,err,error,*999)
      !Get the independent field for the ALE Stokes problem
      NULLIFY(solvers)
      CALL Solver_SolversGet(solver,solvers,err,error,*999)
      NULLIFY(solverALEStokes)
      CALL Solvers_SolverGet(solver%solvers,2,solverALEStokes,err,error,*999)
      NULLIFY(solverEquationsALEStokes)
      CALL Solver_SolverEquationsGet(solverALEStokes,solverEquationsALEStokes,err,error,*999)
      NULLIFY(solverMappingALEStokes)
      CALL SolverEquations_SolverMappingGet(solverEquationsALEStokes,solverMappingALEStokes,err,error,*999)
      NULLIFY(equationsSetALEStokes)
      CALL SolverMapping_EquationsSetGet(solverMappingALEStokes,1,equationsSetALEStokes,err,error,*999)
      NULLIFY(geometricFieldALEStokes)
      CALL EquationsSet_GeometricFieldGet(equationsSetALEStokes,geometricFieldALEStokes,err,error,*999)
      NULLIFY(geometricVariableALEStokes)
      CALL Field_VariableGet(geometricFieldALEStokes,FIELD_U_VARIABLE_TYPE,geometricVariableALEStokes,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(geometricVariableALEStokes,numberOfDimensionsALEStokes,err,error,*999)
      NULLIFY(independentFieldALEStokes)
      CALL EquationsSet_IndependentFieldGet(equationsSetALEStokes,independentFieldALEStokes,err,error,*999)
      NULLIFY(independentVariableALEStokes)
      CALL Field_VariableGet(independentFieldALEStokes,FIELD_U_VARIABLE_TYPE,independentVariableALEStokes,err,error,*999)
      !Copy result from Laplace mesh movement to Stokes' independent field
      IF(numberOfDimensionsALEStokes/=numberOfDimensionsLaplace) THEN
        localError="The number of dimensions in the Laplace problem of "// &
          & TRIM(NumberToVString(numberOfDimensionsLaplace,"*",err,error))// &
          & " does not match the number of dimensions in the ALE Stokes problem of "// &
          & TRIM(NumberToVString(numberOfDimensionsALEStokes,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      DO componentIdx=1,numberOfDimensionsALEStokes
        CALL Field_ParametersToFieldParametersCopy(dependentFieldLaplace,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & componentIdx,independentFieldALEStokes,FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE,componentIdx, &
          & err,error,*999)
      ENDDO !componentIdx
      !Use calculated values to update mesh
      NULLIFY(meshDisplacementValues)
      CALL FieldVariable_ParameterSetDataGet(independentVariableALEStokes,FIELD_MESH_DISPLACEMENT_SET_TYPE, &
        & meshDisplacementValues,err,error,*999)
      NULLIFY(equations)
      CALL EquationsSet_EquationsGet(equationsSetLaplace,equations,err,error,*999)
      NULLIFY(vectorEquations)
      CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
      NULLIFY(vectorMapping)
      CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
      CALL Field_NumberOfVariablesGet(dependentFieldALEStokes,numberOfVariables,err,error,*999)
      DO variableIdx=1,numberOfVariables
        NULLIFY(dependentVariableALEStokes)
        CALL Field_VariableIndexGet(dependentFieldALEStokes,variableIdx,dependentVariableALEStokes,variableType,err,error,*999)
        CALL FieldVariable_NumberOfComponents(dependentVariableALEStokes,numberOfComponents,err,error,*999)
        DO componentIdx=1,numberOfComponents
          NULLIFY(domain)
          CALL FieldVariable_ComponentDomainGet(dependentVariableALEStokes,componentIdx,domain,err,error,*999)
          NULLIFY(domainTopology)
          CALL Domain_DomainTopology(domain,domainTopology,err,error,*999)
          NULLIFY(domainNodes)
          CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
          !Loop over the local nodes excluding the ghosts.
          CALL DomainNodes_NumberOfNodesGet(domainNodes,numberOfNodes,err,error,*999)
          DO nodeIdx=1,numberOfNodes
            CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,nodeIdx,numberOfNodeDerivatives,err,error,*999)
            DO derivativeIdx=1,numberOfNodeDerivatives
              !Default to version 1 of each node derivative
              CALL FieldVariable_LocalNodeDOFGet(dependentVariableALEStokes,1,derivativeIdx,nodeIdx,componentIdx, &
                & localDOFIdx,err,error,*999)
              CALL FieldVariable_ParameterSetAddLocalDOF(geometricVariableALEStokes,FIELD_VALUES_SET_TYPE,localDOFIdx, &
                & meshDisplacementValues(localDOFIdx),err,error,*999)
            ENDDO !derivativeIdx
          ENDDO !nodeIdx
        ENDDO !componentIdx
      ENDDO !variableIdx
      CALL FieldVariable_ParameterSetDataRestore(independentVariableALEStokes,FIELD_MESH_DISPLACEMENT_SET_TYPE, &
        & meshDisplacementValues,err,error,*999)
      CALL FieldVariable_ParameterSetUpdateStart(geometricVariableALEStokes,FIELD_VALUES_SET_TYPE,err,error,*999)
      CALL FieldVariable_ParameterSetUpdateFinish(geometricVariableALEStokes,FIELD_VALUES_SET_TYPE,err,error,*999)
      !Now use displacement values to calculate velocity values
      alpha=1.0_DP/timeIncrement
      CALL FieldVariable_ParameterSetsCopy(independentVariableALEStokes,FIELD_MESH_DISPLACEMENT_SET_TYPE, &
        & FIELD_MESH_VELOCITY_SET_TYPE,alpha,err,error,*999)
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a Stokes equation fluid type of a fluid mechanics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("Stokes_PreSolveALEUpdateMesh")    
    RETURN
999 ERRORSEXITS("Stokes_PreSolveALEUpdateMesh",err,error)
    RETURN 1
    
  END SUBROUTINE Stokes_PreSolveALEUpdateMesh

  !
  !================================================================================================================================
  !
  
  !>Update mesh parameters for three component Laplace problem
  SUBROUTINE Stokes_PreSolveALEUpdateParameters(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx,nodeIdx,derivativeIdx,localDOFIdx,numberOfComponents,numberOfNodes,numberOfNodeDerivatives, &
      & numberOfVariables,pSpecification(3),solveType,variableIdx,variableType
    REAL(DP) :: currentTime,timeIncrement
    REAL(DP), POINTER :: meshStiffValues(:)
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(FieldType), POINTER :: dependentField,independentField
    TYPE(FieldVariableType), POINTER :: dependentVariable
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping 
    TYPE(VARYING_STRING) :: localError

    ENTERS("Stokes_PreSolveALEUpdateParameters",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_Problem(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    
    CALL ControlLoop_CurrentTimesGet(controlLoop,currentTime,timeIncrement,err,error,*999)
    
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STATIC_STOKES_SUBTYPE,PROBLEM_LAPLACE_STOKES_SUBTYPE)
      ! do nothing ???
    CASE(PROBLEM_TRANSIENT_STOKES_SUBTYPE)
      ! do nothing ???
    CASE(PROBLEM_ALE_STOKES_SUBTYPE)
      CALL Solver_TypeGet(solver,solveType,err,error,*999)
      IF(solveType==SOLVER_LINEAR_TYPE) THEN
        !Get the independent field for the ALE Stokes problem
        NULLIFY(solverEquations)
        CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
        NULLIFY(solverMapping)
        CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
        NULLIFY(dependentField)
        CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
        NULLIFY(independentField)
        CALL EquationsSet_IndependentFieldGet(equationsSet,independentField,err,error,*999)
        NULLIFY(meshStiffValues)
        CALL Field_ParameterSetDataGet(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,meshStiffValues,err,error,*999)
        CALL Field_NumberOfVariablesGet(dependentField,numberOfVariables,err,error,*999)
        DO variableIdx=1,numberOfVariables
          NULLIFY(dependentVariable)
          CALL Field_VariableIndexGet(dependentField,variableIdx,dependentVariable,variableType,err,error,*999)
          CALL FieldVariable_NumberOfComponentsGet(dependentVariable,numberOfComponents,err,error,*999)
          DO componentIdx=1,numberOfComponents
            NULLIFY(domain)
            CALL FieldVariable_ComponentDomainGet(dependentVariable,componentIdx,domain,err,error,*999)
            NULLIFY(domainTopology)
            CALL Doman_DomainTopologyGet(domain,domainTopology,err,error,*999)
            NULLIFY(domainNodes)
            CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
            !Loop over the local nodes excluding the ghosts.
            CALL DomainNodes_NumberOfNodesGet(domainNodes,numberOfNodes,err,error,*999)
            DO nodeIdx=1,numberOfNodes
              CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,nodeIdx,numberOfNodeDerivatives,err,error,*999)
              DO derivativeIdx=1,numberOfNodeDerivatives
                !Default to version 1 of each node derivative
                CALL FieldVariable_LocalNodeDOFGet(dependentVariable,1,derivativeIdx,nodeIdx,componentIdx,localDOFIdx, &
                  & err,error,*999)
                !Calculation of K values dependent on current mesh topology
                meshStiffValues(localDOFIdx)=1.0_DP
                CALL Field_ParameterSetUpdateLocalDOF(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,localDOFIdx, &
                  & meshStiffValues(localDOFIdx),err,error,*999)
              ENDDO !derivativeIdx
            ENDDO !nodeIdx
          ENDDO !componentIdx
        ENDDO !variableIdx
        CALL Field_ParameterSetDataRestore(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,meshStiffValues, &
          & err,error,*999)
      ELSE IF(solveType==SOLVER_DYNAMIC_TYPE) THEN
        CALL FlagError("Mesh motion calculation not successful for ALE problem.",err,error,*999)
      ENDIF
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a Stokes equation fluid type of a fluid mechanics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("Stokes_PreSolveALEUpdateParameters")
    RETURN
999 ERRORSEXITS("Stokes_PreSolveALEUpdateParameters",err,error)
    RETURN 1
    
  END SUBROUTINE Stokes_PreSolveALEUpdateParameters

  !
  !================================================================================================================================
  !

  !>Output data post solve
  SUBROUTINE Stokes_PostSolveOutputData(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: analyticFunctionType,currentIteration,equationsSetIdx,inputIteration,numberOfDimensions, &
      & numberOfEquationsSets,outputIteration,pSpecification(3),outputType
    REAL(DP) :: currentTime,startTime,stopTime,timeIncrement
    LOGICAL :: exportField
    CHARACTER(14) :: outputFile
    TYPE(ControlLoopType), POINTER :: controlLoop 
    TYPE(EquationsSetType), POINTER :: equationsSet 
    TYPE(EquationsSetAnalyticType), POINTER :: equationsAnalytic
    TYPE(FieldType), POINTER :: dependentField
    TYPE(FieldsType), POINTER :: fields
    TYPE(ProblemType), POINTER :: problem
    TYPE(RegionType), POINTER :: region
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(VARYING_STRING) :: localError,method,filename

    ENTERS("Stokes_PostSolveOutputData",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)

    CALL Solver_OutputTypeGet(solver,outputType,err,error,*999)
    CALL System('mkdir -p ./output')
    
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STATIC_STOKES_SUBTYPE,PROBLEM_LAPLACE_STOKES_SUBTYPE)
      NULLIFY(solverEquations)
      CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      !Make sure the equations sets are up to date
      CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
      DO equationsSetIdx=1,numberOfEquationsSets
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
        NULLIFY(region)
        CALL EquationsSet_RegionGet(equationsSet,region,err,error,*999)
        NULLIFY(fields)
        CALL Region_FieldsGet(region,fields,err,error,*999)
        filename="./output/"//"StaticSolution"
        method="FORTRAN"
        IF(outputType>=SOLVER_PROGRESS_OUTPUT) THEN
          CALL WriteString(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
          CALL WriteString(GENERAL_OUTPUT_TYPE,"Now export fields... ",err,error,*999)
        ENDIF
        CALL FIELD_IO_NODES_EXPORT(fields,filename,method,err,error,*999)
        CALL FIELD_IO_ELEMENTS_EXPORT(fields,filename,method,err,error,*999)
      ENDDO !equationsSetIdx
    CASE(PROBLEM_TRANSIENT_STOKES_SUBTYPE,PROBLEM_ALE_STOKES_SUBTYPE)
      CALL ControlLoop_CurrentTimeInformationGet(controlLoop,currentTime,timeIncrement,startTime,stopTime,currentIteration, &
        & outputIteration,inputIteration,err,error,*999)
      NULLIFY(solverEquations)
      CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      !Make sure the equations sets are up to date
      CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
      DO equationsSetIdx=1,numberOfEquationsSets
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
        IF(outputIteration/=0) THEN
          IF(currentTime<=stopTime) THEN
            IF(currentIteration<10) THEN
              WRITE(outputFile,'("TimeStep_000",I0)') currentIteration
            ELSE IF(currentIteration<100) THEN
              WRITE(outputFile,'("TimeStep_00",I0)') currentIteration
            ELSE IF(currentIteration<1000) THEN
              WRITE(outputFile,'("TimeStep_0",I0)') currentIteration
            ELSE IF(currentIteration<10000) THEN
              WRITE(outputFile,'("TimeStep_",I0)') currentIteration
            END IF
            filename="./output/"//"MainTime_"//TRIM(NumberToVString(currentIteration,"*",err,error))
            method="FORTRAN"
            IF(MOD(currentIteration,outputIteration)==0)  THEN
              IF(outputtype >= SOLVER_PROGRESS_OUTPUT) THEN
                CALL WriteString(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
                CALL WriteString(GENERAL_OUTPUT_TYPE,"Now export fields... ",err,error,*999)
              ENDIF
              NULLIFY(region)
              CALL EquationsSet_RegionGet(equationsSet,region,err,error,*999)
              NULLIFY(fields)
              CALL Region_FieldsGet(region,fields,err,error,*999)
              CALL FIELD_IO_NODES_EXPORT(fields,filename,method,err,error,*999)
              CALL FIELD_IO_ELEMENTS_EXPORT(fields,filename,method,err,error,*999)
              IF(outputType >= SOLVER_PROGRESS_OUTPUT) THEN
                CALL WriteString(GENERAL_OUTPUT_TYPE,filename,err,error,*999)
                CALL WriteString(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
              ENDIF
            ENDIF
            NULLIFY(equationsAnalytic)
            CALL EquationsSet_AnalyticExists(equationsSet,equationsAnalytic,err,error,*999)            
            IF(ASSOCIATED(equationsAnalytic)) THEN
              CALL EquationsSet_AnalyticFunctionTypeGet(equationsSet,analyticFunctionType,err,error,*999)
              IF(analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_4.OR. &
                & analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_5.OR. &
                & analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_4.OR. &
                & analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_5.OR. &
                & analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_1) THEN
                NULLIFY(dependentField)
                CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
                CALL AnalyticAnalysis_Output(dependentField,outputFile,err,error,*999)
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDDO !equationsSetIdx
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a Stokes equation fluid type of a fluid mechanics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("Stokes_PostSolveOutputData")
    RETURN
999 ERRORSEXITS("Stokes_PostSolveOutputData",err,error)
    RETURN 1
    
  END SUBROUTINE Stokes_PostSolveOutputData

  !
  !================================================================================================================================
  !

  !>Calculates the analytic solution and sets the boundary conditions for an analytic problem.
  SUBROUTINE Stokes_BoundaryConditionsAnalyticCalculate(equationsSet,boundaryConditions,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
!\todo: Reduce number of variables used
    INTEGER(INTG) :: analyticFunctionType,boundCount,componentIdx,derivativeIdx,dimensionIdx,elementIdx,globalDerivativeIndex, &
      & i,j,k,localDOFIdx,localNodeIdx,localNodeNumber,maximumNumberOfElementParameters,nodeIdx,numberOfComponents, &
      & numberOfDimensions,numberOfNodes,numberOfNodeDerivatives,numberOfNodesXiC(4),numberOfVariables,variableIdx,variableType
    REAL(DP) :: analyticValue,currentTime,muParam,rhoParam,tCoordinates(20,3),x(3),xiCoordinates(3)
    !REAL(DP) :: boundaryTolerance, boundaryX(3,2),muParam,L
    REAL(DP), POINTER :: geometricParameters(:),materialsParameters(:)
    LOGICAL :: boundaryNode
    TYPE(BasisType), POINTER :: basis
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainElementsType), POINTER :: domainElements
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(EquationsSetAnalyticType), POINTER :: equationsAnalytic
    TYPE(FieldType), POINTER :: dependentField,geometricField,materialsField
    TYPE(FieldVariableType), POINTER :: dependentVariable,geometricVariable,materialsVariable
    TYPE(FieldInterpolatedPointType), POINTER :: interpolatedPoint
    TYPE(FieldInterpolationParametersType), POINTER :: interpolationParameters
    !TYPE(VARYING_STRING) :: localError

    ENTERS("Stokes_BoundaryConditionsAnalyticCalculate",err,error,*999)
    
!\todo: Introduce user call to set parameters
    boundCount=0
    !L=10.0_DP
    xiCoordinates(3)=0.0_DP    
    !boundaryTolerance=0.000000001_DP

    IF(.NOT.ASSOCIATED(boundaryConditions)) CALL FlagError("Boundary conditions is not associated.",err,error,*999)
    NULLIFY(equationsAnalytic)
    CALL EquationsSet_AnalyticGet(equationsSet,equationsAnalytic,err,error,*999)
    CALL EquationsSet_AnalyticFunctionTypeGet(equationsSet,analyticFunctionType,err,error,*999)
    NULLIFY(geometricField)
    CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
    NULLIFY(geometricVariable)
    CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,geometricVariable,err,error,*999)
    CALL FieldVariable_NumberOfComponentsGet(geometricVariable,numberOfDimensions,err,error,*999)
    NULLIFY(geometricParameters)
    CALL FieldVariable_ParameterSetDataGet(geometricVariable,FIELD_VALUES_SET_TYPE,geometricParameters,err,error,*999)
    NULLIFY(interpolationParameters)
    CALL FieldVariable_InterpolationParametersInitialise(geometricVariable,interpolationParameters,err,error,*999)
    NULLIFY(interpolatedPoint)
    CALL Field_InterpolatedPointInitialise(interpolationParameters,interpolatedPoint,err,error,*999)
    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
    NULLIFY(materialsField)
    CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)
    NULLIFY(materialsVariable)
    CALL Field_VariableGet(materialsField,FIELD_U_VARIABLE_TYPE,materialsVariable,err,error,*999)
    NULLIFY(materialsParameters)
    CALL FieldVariable_ParameterSetDataGet(materialsVariable,FIELD_VALUES_SET_TYPE,materialsParameters,err,error,*999)
    !Define muParam, density=1
    muParam=materialsParameters(1)
    !Define rhoParam, density=2
    IF(analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_4.OR. &
      & analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_5.OR. &
      & analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_4.OR. &
      & analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_5) THEN
      rhoParam=materialsParameters(2)
    ELSE
      rhoParam=0.0_DP
    ENDIF

    !\todo: Check adjacent element calculation / use boundary node flag instead / didn't work for simplex
    !IF(numberOfDimensions==2) THEN
    !  boundaryX(1,1)=0.0_DP
    !  boundaryX(1,2)=10.0_DP
    !  boundaryX(2,1)=0.0_DP
    !  boundaryX(2,2)=10.0_DP
    !ELSE IF(numberOfDimensions==3) THEN
    !  boundaryX(1,1)=-5.0_DP
    !  boundaryX(1,2)=5.0_DP
    !  boundaryX(2,1)=-5.0_DP
    !  boundaryX(2,2)=5.0_DP
    !  boundaryX(3,1)=-5.0_DP
    !  boundaryX(3,2)=5.0_DP
    !ENDIF

    CALL Field_NumberOfVariablesGet(dependentField,numberOfVariables,err,error,*999)
    DO variableIdx=1,numberOfVariables
      NULLIFY(dependentVariable)
      CALL Field_VariableIndexGet(dependentField,variableIdx,dependentVariable,variableType,err,error,*999)
      CALL FieldVariable_ParameterSetEnsureCreated(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(dependentVariable,numberOfComponents,err,error,*999)
      DO componentIdx=1,numberOfComponents
        boundCount=0
        CALL FieldVariable_ComponentInterpolationCheck(dependentVariable,componentIdx,FIELD_NODE_BASED_INTERPOLATION, &
          & err,error,*999)
        NULLIFY(domain)
        CALL FieldVariable_ComponentDomainGet(dependentVariable,componentIdx,domain,err,error,*999)
        NULLIFY(domainTopology)
        CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
        NULLIFY(domainNodes)
        CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
        NULLIFY(domainElements)
        CALL DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*999)
        !Loop over the local nodes excluding the ghosts.
        CALL DomainNodes_NumberOfNodesGet(domainNodes,numberOfNodes,err,error,*999)
        DO nodeIdx=1,numberOfNodes
          CALL DomainNodes_NodeBoundaryNodeGet(domainNodes,nodeIdx,boundaryNode,err,error,*999)
          CALL DomainNodes_NodeSurroundingElementGet(domainNodes,1,nodeIdx,elementIdx,err,error,*999)
          NULLIFY(basis)
          CALL DomainElements_ElementBasisGet(domainElements,elementIdx,basis,err,error,*999)
          CALL DomainElements_MaxElementParametersGet(domainElements,maximumNumberOfElementParameters,err,error,*999)
          CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementIdx,interpolationParameters,err,error,*999)
          localNodeIdx=0
          xiCoordinates=0.0_DP
          numberOfNodesXiC=1
          CALL Basis_NumberOfNodesXiCGet(basis,numberOfNodesXiC,err,error,*999)
          !\todo: Use boundary flag
          IF(maximumNumberOfElementParameters==4.AND.numberOfDimensions==2 .OR. &
            & maximumNumberOfElementParameters==9.OR. &
            & maximumNumberOfElementParameters==16.OR. &
            & maximumNumberOfElementParameters==8.OR. &
            & maximumNumberOfElementParameters==27.OR. &
            & maximumNumberOfElementParameters==64) THEN
            DO k=1,numberOfNodesXic(3)
              DO j=1,numberOfNodesXic(2)
                DO i=1,numberOfNodesXic(1)
                  localNodeIdx=localNodeIdx+1
                  CALL DomainElements_ElementNodeGet(domainElements,localNodeIdx,elementIdx,localNodeNumber,err,error,*999)
                  IF(localNodeNumber==nodeIdx) EXIT
                  xiCoordinates(1)=xiCoordinates(1)+(1.0_DP/(numberOfNodesXic(1)-1))
                ENDDO !i
                IF(localNodeNumber==nodeIdx) EXIT
                xiCoordinates(1)=0.0_DP
                xiCoordinates(2)=xiCoordinates(2)+(1.0_DP/(numberOfNodesXic(2)-1))
              ENDDO !j
              IF(localNodeNumber==nodeIdx) EXIT
              xiCoordinates(1)=0.0_DP
              xiCoordinates(2)=0.0_DP
              IF(numberOfNodesXic(3)/=1) THEN
                xiCoordinates(3)=xiCoordinates(3)+(1.0_DP/(numberOfNodesXic(3)-1))
              ENDIF
            ENDDO !k
            CALL Field_InterpolateXi(NO_PART_DERIV,xiCoordinates,interpolatedPoint,err,error,*999)
          ELSE
            !\todo: Use boundary flag
            IF(maximumNumberOfElementParameters==3) THEN
              tCoordinates(1,1:2)=[0.0_DP,1.0_DP]
              tCoordinates(2,1:2)=[1.0_DP,0.0_DP]
              tCoordinates(3,1:2)=[1.0_DP,1.0_DP]
            ELSE IF(maximumNumberOfElementParameters==6) THEN
              tCoordinates(1,1:2)=[0.0_DP,1.0_DP]
              tCoordinates(2,1:2)=[1.0_DP,0.0_DP]
              tCoordinates(3,1:2)=[1.0_DP,1.0_DP]
              tCoordinates(4,1:2)=[0.5_DP,0.5_DP]
              tCoordinates(5,1:2)=[1.0_DP,0.5_DP]
              tCoordinates(6,1:2)=[0.5_DP,1.0_DP]
            ELSE IF(maximumNumberOfElementParameters==10.AND.numberOfDimensions==2) THEN
              tCoordinates(1,1:2)=[0.0_DP,1.0_DP]
              tCoordinates(2,1:2)=[1.0_DP,0.0_DP]
              tCoordinates(3,1:2)=[1.0_DP,1.0_DP]
              tCoordinates(4,1:2)=[1.0_DP/3.0_DP,2.0_DP/3.0_DP]
              tCoordinates(5,1:2)=[2.0_DP/3.0_DP,1.0_DP/3.0_DP]
              tCoordinates(6,1:2)=[1.0_DP,1.0_DP/3.0_DP]
              tCoordinates(7,1:2)=[1.0_DP,2.0_DP/3.0_DP]
              tCoordinates(8,1:2)=[2.0_DP/3.0_DP,1.0_DP]
              tCoordinates(9,1:2)=[1.0_DP/3.0_DP,1.0_DP]
              tCoordinates(10,1:2)=[2.0_DP/3.0_DP,2.0_DP/3.0_DP]
            ELSE IF(maximumNumberOfElementParameters==4) THEN
              tCoordinates(1,1:3)=[0.0_DP,1.0_DP,1.0_DP]
              tCoordinates(2,1:3)=[1.0_DP,0.0_DP,1.0_DP]
              tCoordinates(3,1:3)=[1.0_DP,1.0_DP,0.0_DP]
              tCoordinates(4,1:3)=[1.0_DP,1.0_DP,1.0_DP]
            ELSE IF(maximumNumberOfElementParameters==10.AND.numberOfDimensions==3) THEN
              tCoordinates(1,1:3)=[0.0_DP,1.0_DP,1.0_DP]
              tCoordinates(2,1:3)=[1.0_DP,0.0_DP,1.0_DP]
              tCoordinates(3,1:3)=[1.0_DP,1.0_DP,0.0_DP]
              tCoordinates(4,1:3)=[1.0_DP,1.0_DP,1.0_DP]
              tCoordinates(5,1:3)=[0.5_DP,0.5_DP,1.0_DP]
              tCoordinates(6,1:3)=[0.5_DP,1.0_DP,0.5_DP]
              tCoordinates(7,1:3)=[0.5_DP,1.0_DP,1.0_DP]
              tCoordinates(8,1:3)=[1.0_DP,0.5_DP,0.5_DP]
              tCoordinates(9,1:3)=[1.0_DP,1.0_DP,0.5_DP]
              tCoordinates(10,1:3)=[1.0_DP,0.5_DP,1.0_DP]
            ELSE IF(maximumNumberOfElementParameters==20) THEN
              tCoordinates(1,1:3)=[0.0_DP,1.0_DP,1.0_DP]
              tCoordinates(2,1:3)=[1.0_DP,0.0_DP,1.0_DP]
              tCoordinates(3,1:3)=[1.0_DP,1.0_DP,0.0_DP]
              tCoordinates(4,1:3)=[1.0_DP,1.0_DP,1.0_DP]
              tCoordinates(5,1:3)=[1.0_DP/3.0_DP,2.0_DP/3.0_DP,1.0_DP]
              tCoordinates(6,1:3)=[2.0_DP/3.0_DP,1.0_DP/3.0_DP,1.0_DP]
              tCoordinates(7,1:3)=[1.0_DP/3.0_DP,1.0_DP,2.0_DP/3.0_DP]
              tCoordinates(8,1:3)=[2.0_DP/3.0_DP,1.0_DP,1.0_DP/3.0_DP]
              tCoordinates(9,1:3)=[1.0_DP/3.0_DP,1.0_DP,1.0_DP]
              tCoordinates(10,1:3)=[2.0_DP/3.0_DP,1.0_DP,1.0_DP]
              tCoordinates(11,1:3)=[1.0_DP,1.0_DP/3.0_DP,2.0_DP/3.0_DP]
              tCoordinates(12,1:3)=[1.0_DP,2.0_DP/3.0_DP,1.0_DP/3.0_DP]
              tCoordinates(13,1:3)=[1.0_DP,1.0_DP,1.0_DP/3.0_DP]
              tCoordinates(14,1:3)=[1.0_DP,1.0_DP,2.0_DP/3.0_DP]
              tCoordinates(15,1:3)=[1.0_DP,1.0_DP/3.0_DP,1.0_DP]
              tCoordinates(16,1:3)=[1.0_DP,2.0_DP/3.0_DP,1.0_DP]
              tCoordinates(17,1:3)=[2.0_DP/3.0_DP,2.0_DP/3.0_DP,2.0_DP/3.0_DP]
              tCoordinates(18,1:3)=[2.0_DP/3.0_DP,2.0_DP/3.0_DP,1.0_DP]
              tCoordinates(19,1:3)=[2.0_DP/3.0_DP,1.0_DP,2.0_DP/3.0_DP]
              tCoordinates(20,1:3)=[1.0_DP,2.0_DP/3.0_DP,2.0_DP/3.0_DP]
            ENDIF
            DO k=1,maximumNumberOfElementParameters
              CALL DomainElements_ElementNodeGet(domainElements,k,elementIdx,localNodeNumber,err,error,*999)
              IF(localNodeNumber==nodeIdx) EXIT
            ENDDO
            CALL Field_InterpolateXi(NO_PART_DERIV,tCoordinates(k,1:numberOfDimensions),interpolatedPoint,err,error,*999)
          ENDIF
          x=0.0_DP
          DO dimensionIdx=1,numberOfDimensions
            x(dimensionIdx)=interpolatedPoint%values(dimensionIdx,NO_PART_DERIV)
          ENDDO !dimensionIdx          
          !Loop over the derivatives
          CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,nodeIdx,numberOfNodeDerivatives,err,error,*999)
          DO derivativeIdx=1,numberOfNodeDerivatives
            CALL DomainNodes_DerivativeGlobalIndexGet(domainNodes,derivativeIdx,nodeIdx,globalDerivativeIndex,err,error,*999)
            currentTime=0.0_DP
            CALL Stokes_AnalyticFunctions(analyticValue,x,muParam,rhoParam,currentTime,variableType, &
              & globalDerivativeIndex,analyticFunctionType,numberOfDimensions,numberOfComponents,componentIdx,err,error,*999)
            !Default to version 1 of each node derivative
            CALL FieldVariable_LocalNodeDOFGet(dependentVariable,1,derivativeIdx,nodeIdx,componentIdx,localDOFIdx,err,error,*999)
            CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,localDOFIdx, &
              & analyticValue,err,error,*999)
            IF(variableType==FIELD_U_VARIABLE_TYPE) THEN
              ! \todo: This part should work even for simplex elements as soon as adjacent element calculation has been fixed
              IF(boundaryNode) THEN
                !If we are a boundary node then set the analytic value on the boundary
                IF(componentIdx<=numberOfDimensions) THEN
                  CALL BoundaryConditions_SetLocalDOF(boundaryConditions,dependentVariable,localDOFIdx,BOUNDARY_CONDITION_FIXED, &
                    & analyticValue,err,error,*999)
                  boundCount=boundCount+1
                ELSE
                  ! \todo: This is just a workaround for linear pressure fields in simplex element components
                  IF(maximumNumberOfElementParameters==3) THEN
                    IF(analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_1.OR. &
                      & analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_2.OR. &
                      & analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_3.OR. &
                      & analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_4.OR. &
                      & analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_5) THEN
                      IF(-0.001_DP<X(1).AND.X(1)<0.001_DP.AND.-0.001_DP<X(2).AND.X(2)<0.001_DP.OR. &
                        & 10.0_DP-0.001_DP<X(1).AND.X(1)<10.0_DP+0.001_DP.AND.-0.001_DP<X(2).AND. &
                        & X(2)<0.001_DP.OR. &
                        & 10.0_DP-0.001_DP<X(1).AND.X(1)<10.0_DP+0.001_DP.AND.10.0_DP-0.001_DP<X(2).AND. &
                        & X(2)<10.0_DP+0.001_DP.OR. &
                        & -0.001_DP<X(1).AND.X(1)<0.001_DP.AND.10.0_DP-0.001_DP<X(2).AND. &
                        & X(2)<10.0_DP+0.001_DP) THEN
                        CALL BoundaryConditions_SetLocalDOF(boundaryConditions,dependentVariable,localDOFIdx, &
                          & BOUNDARY_CONDITION_FIXED,analyticValue,err,error,*999)
                        boundCount=boundCount+1
                      ENDIF
                    ENDIF
                  ELSE IF(maximumNumberOfElementParameters==4.AND.numberOfDimensions==3) THEN
                    IF(analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_1.OR. &
                      & analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_2.OR. &
                      & analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_3.OR. &
                      & analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_4.OR. &
                      & analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_5) THEN
                      IF(-5.0_DP-0.001_DP<X(1).AND.X(1)<-5.0_DP+0.001_DP.AND.-5.0_DP-0.001_DP<X(2).AND. &
                        & X(2)<-5.0_DP+0.001_DP.AND.-5.0_DP-0.001_DP<X(3).AND.X(3)<-5.0_DP+0.001_DP.OR. &
                        & -5.0_DP-0.001_DP<X(1).AND.X(1)<-5.0_DP+0.001_DP.AND.5.0_DP-0.001_DP<X(2).AND. &
                        & X(2)<5.0_DP+0.001_DP.AND.-5.0_DP-0.001_DP<X(3).AND.X(3)<-5.0_DP+0.001_DP.OR. &
                        & 5.0_DP-0.001_DP<X(1).AND.X(1)<5.0_DP+0.001_DP.AND.5.0_DP-0.001_DP<X(2).AND. &
                        & X(2)<5.0_DP+0.001_DP.AND.-5.0_DP-0.001_DP<X(3).AND.X(3)<-5.0_DP+0.001_DP.OR. &
                        & 5.0_DP-0.001_DP<X(1).AND.X(1)<5.0_DP+0.001_DP.AND.-5.0_DP-0.001_DP<X(2).AND. &
                        & X(2)<-5.0_DP+0.001_DP.AND.-5.0_DP-0.001_DP<X(3).AND.X(3)<-5.0_DP+0.001_DP.OR. &
                        & -5.0_DP-0.001_DP<X(1).AND.X(1)<-5.0_DP+0.001_DP.AND.-5.0_DP-0.001_DP<X(2).AND. &
                        & X(2)<-5.0_DP+0.001_DP.AND.5.0_DP-0.001_DP<X(3).AND.X(3)<5.0_DP+0.001_DP.OR. &
                        & -5.0_DP-0.001_DP<X(1).AND.X(1)<-5.0_DP+0.001_DP.AND.5.0_DP-0.001_DP<X(2).AND. &
                        & X(2)<5.0_DP+0.001_DP.AND.5.0_DP-0.001_DP<X(3).AND.X(3)<5.0_DP+0.001_DP.OR. &
                        & 5.0_DP-0.001_DP<X(1).AND.X(1)<5.0_DP+0.001_DP.AND.5.0_DP-0.001_DP<X(2).AND. &
                        & X(2)<5.0_DP+0.001_DP.AND.5.0_DP-0.001_DP<X(3).AND.X(3)<5.0_DP+0.001_DP.OR. &
                        & 5.0_DP-0.001_DP<X(1).AND.X(1)<5.0_DP+0.001_DP.AND.-5.0_DP-0.001_DP<X(2).AND. &
                        & X(2)<-5.0_DP+ 0.001_DP.AND.5.0_DP-0.001_DP<X(3).AND.X(3)<5.0_DP+0.001_DP) THEN
                        CALL BoundaryConditions_SetLocalDOF(boundaryConditions,dependentVariable,localDOFIdx, &
                          & BOUNDARY_CONDITION_FIXED,analyticValue,err,error,*999)
                        boundCount=boundCount+1
                      ENDIF
                    ENDIF
                    ! \todo: This is how it should be if adjacent elements would be working
                  ELSE IF(boundCount==0) THEN
                    CALL BoundaryConditions_SetLocalDOF(boundaryConditions,dependentVariable,localDOFIdx, &
                      & BOUNDARY_CONDITION_FIXED,analyticValue,err,error,*999)
                    boundCount=boundCount+1
                  ENDIF
                ENDIF
              ELSE
                IF(componentIdx<=numberOfDimensions) THEN
                  CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_VALUES_SET_TYPE,localDOFIdx, &
                    & analyticValue,err,error,*999)
                ENDIF
              ENDIF
            ENDIF
          ENDDO !derivativeIdx
        ENDDO !nodeIdx
      ENDDO !componentIdx
      CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
      CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
      CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
      CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
    ENDDO !variableIdx
    CALL FieldVariable_ParameterSetDataRestore(geometricVariable,FIELD_VALUES_SET_TYPE,geometricParameters,err,error,*999)
    CALL FieldVariable_ParameterSetDataRestore(materialsVariable,FIELD_VALUES_SET_TYPE,materialsParameters,err,error,*999)
    CALL Field_InterpolatedPointFinalise(interpolatedPoint,err,error,*999)
    CALL FieldVariable_InterpolationParametersFinalise(interpolationParameters,err,error,*999)

    EXITS("Stokes_BoundaryConditionsAnalyticCalculate")
    RETURN
999 ERRORSEXITS("Stokes_BoundaryConditionsAnalyticCalculate",err,error)
    RETURN 1

  END SUBROUTINE Stokes_BoundaryConditionsAnalyticCalculate

  !
  !================================================================================================================================
  !
  !>Calculates the various analytic solutions given X and time, can be called from within analytic calculate or elsewhere if needed
  SUBROUTINE Stokes_AnalyticFunctions(analyticValue,X,muParam,rhoParam,currentTime,variableType,globalDerivativeIndex, &
    & analyticFunctionType,numberOfDimensions,numberOfComponents,componentIdx,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    REAL(DP), INTENT(OUT) :: analyticValue
    REAL(DP) :: muParam,rhoParam
    REAL(DP), INTENT(IN) :: currentTime
    REAL(DP), INTENT(IN), DIMENSION(3) :: X
    INTEGER(INTG), INTENT(IN) :: numberOfDimensions,numberOfComponents,componentIdx
    !Local variables
    INTEGER(INTG) :: variableType,globalDerivativeIndex,analyticFunctionType
    REAL(DP) :: internalTime
    !TYPE(DomainType), POINTER :: domain
    !TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(VARYING_STRING) :: localError

    ENTERS("Stokes_AnalyticFunctions",err,error,*999)

    !\todo: Introduce user-defined or default values instead for density and viscosity
    internalTime=currentTime
    SELECT CASE(analyticFunctionType)
    CASE(EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_1)
      IF(numberOfDimensions==2.AND.numberOfComponents==3) THEN
        !Polynomial function
        SELECT CASE(variableType)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(globalDerivativeIndex)
          CASE(NO_GLOBAL_DERIV)
            IF(componentIdx==1) THEN
              !calculate u
              analyticValue=x(2)**2/10.0_DP**2
            ELSE IF(componentIdx==2) THEN
              !calculate v
              analyticValue=x(1)**2/10.0_DP**2
            ELSE IF(componentIdx==3) THEN
              !calculate p
              analyticValue=2.0_DP*muParam/10.0_DP**2*x(1)
            ELSE
              CALL FlagError("Not implemented.",err,error,*999)
            ENDIF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(globalDerivativeIndex)
          CASE(NO_GLOBAL_DERIV)
            analyticValue= 0.0_DP
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The variable type of "//TRIM(NumberToVString(variableType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        localError="The number of components does not correspond to the number of dimensions."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE(EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_2)
      IF(numberOfDimensions==2.AND.numberOfComponents==3) THEN
        !Exponential function
        SELECT CASE(variableType)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(globalDerivativeIndex)
          CASE(NO_GLOBAL_DERIV)
            IF(componentIdx==1) THEN
              !calculate u
              analyticValue= EXP((x(1)-x(2))/10.0_DP)
            ELSE IF(componentIdx==2) THEN
              !calculate v
              analyticValue= EXP((x(1)-x(2))/10.0_DP)
            ELSE IF(componentIdx==3) THEN
              !calculate p
              analyticValue= 2.0_DP*muParam/10.0_DP*EXP((x(1)-x(2))/10.0_DP)
            ELSE
              CALL FlagError("Not implemented.",err,error,*999)
            ENDIF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(globalDerivativeIndex)
          CASE(NO_GLOBAL_DERIV)
            IF(componentIdx==1) THEN
              !calculate u
              analyticValue= 0.0_DP
            ELSE IF(componentIdx==2) THEN
              !calculate v
              analyticValue= 0.0_DP
            ELSE IF(componentIdx==3) THEN
              !calculate p
              analyticValue= 0.0_DP
            ELSE
              CALL FlagError("Not implemented.",err,error,*999)
            ENDIF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The variable type of "//TRIM(NumberToVString(variableType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        localError="The number of components does not correspond to the number of dimensions."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE(EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_3)
      IF(numberOfDimensions==2.AND.numberOfComponents==3) THEN
        !Sine and cosine functions
        SELECT CASE(variableType)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(globalDerivativeIndex)
          CASE(NO_GLOBAL_DERIV)
            IF(componentIdx==1) THEN
              !calculate u
              analyticValue=SIN(2.0_DP*PI*x(1)/10.0_DP)*SIN(2.0_DP*PI*x(2)/10.0_DP)
            ELSE IF(componentIdx==2) THEN
              !calculate v
              analyticValue=COS(2.0_DP*PI*x(1)/10.0_DP)*COS(2.0_DP*PI*x(2)/10.0_DP)
            ELSE IF(componentIdx==3) THEN
              !calculate p
              analyticValue=4.0_DP*muParam*PI/10.0_DP*SIN(2.0_DP*PI*x(2)/10.0_DP)*COS(2.0_DP*PI*x(1)/10.0_DP)
            ELSE
              CALL FlagError("Not implemented.",err,error,*999)
            ENDIF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(globalDerivativeIndex)
          CASE(NO_GLOBAL_DERIV)
            IF(componentIdx==1) THEN
              !calculate u
              analyticValue=0.0_DP
            ELSE IF(componentIdx==2) THEN
              !calculate v
              analyticValue=16.0_DP*muParam*PI**2/10.0_DP**2*cos(2.0_DP*PI*x(2)/10.0_DP)*cos(2.0_DP*PI*x(1)/10.0_DP)
            ELSE IF(componentIdx==3) THEN
              !calculate p
              analyticValue=0.0_DP
            ELSE
              CALL FlagError("Not implemented.",err,error,*999)
            ENDIF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString( &
              & globalDerivativeIndex,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        localError="The number of components does not correspond to the number of dimensions."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE(EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_4)
      IF(numberOfDimensions==2.AND.numberOfComponents==3) THEN
        !Reduced Taylor-Green solution for Stokes
        CALL FlagError("Not implemented.",err,error,*999)
      ENDIF
    CASE(EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_5)
      IF(numberOfDimensions==2.AND.numberOfComponents==3) THEN
        !Stokes-Taylor-Green dynamic
        SELECT CASE(variableType)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(globalDerivativeIndex)
          CASE(NO_GLOBAL_DERIV)
            IF(componentIdx==1) THEN
              !calculate u
              analyticValue=x(2)*exp(-(2.0_DP*muParam/rhoParam*currentTime))
            ELSE IF(componentIdx==2) THEN
              !calculate v
              analyticValue=x(1)*exp(-(2.0_DP*muParam/rhoParam*currentTime))
            ELSE IF(componentIdx==3) THEN
              !calculate p
              analyticValue=2.0_DP*x(2)*muParam*exp(-(2.0_DP*muParam/rhoParam*currentTime))*x(1)
            ELSE
              CALL FlagError("Not implemented.",err,error,*999)
            ENDIF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString( &
              & globalDerivativeIndex,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(globalDerivativeIndex)
          CASE(NO_GLOBAL_DERIV)
            IF(componentIdx==1) THEN
              !calculate u
              analyticValue=0.0_DP
            ELSE IF(componentIdx==2) THEN
              !calculate v
              analyticValue=0.0_DP
            ELSE IF(componentIdx==3) THEN
              !calculate p
              analyticValue=0.0_DP
            ELSE
              CALL FlagError("Not implemented.",err,error,*999)
            ENDIF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString( &
              & globalDerivativeIndex,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The variable type of "//TRIM(NumberToVString(variableType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        localError="The number of components does not correspond to the number of dimensions."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE(EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_1)
      IF(numberOfDimensions==3.AND.numberOfComponents==4) THEN
        !POLYNOM
        SELECT CASE(variableType)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(globalDerivativeIndex)
          CASE(NO_GLOBAL_DERIV)
            IF(componentIdx==1) THEN
              !calculate u
              analyticValue=x(2)**2/10.0_DP**2+x(3)**2/10.0_DP**2
            ELSE IF(componentIdx==2) THEN
              !calculate v
              analyticValue=x(1)**2/10.0_DP**2+x(3)**2/10.0_DP**2
            ELSE IF(componentIdx==3) THEN
              !calculate w
              analyticValue=x(1)**2/10.0_DP**2+x(2)**2/10.0_DP**2
            ELSE IF(componentIdx==4) THEN
              !calculate p
              analyticValue=4.0_DP*muParam/10.0_DP**2*x(1)
            ELSE
              CALL FlagError("Not implemented.",err,error,*999)
            ENDIF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(globalDerivativeIndex)
          CASE(NO_GLOBAL_DERIV)
            analyticValue=0.0_DP
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString( &
              & globalDerivativeIndex,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The variable type of "//TRIM(NumberToVString(variableType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        localError="The number of components does not correspond to the number of dimensions."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE(EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_2)
      IF(numberOfDimensions==3.AND.numberOfComponents==4) THEN
        !Exponential function
        SELECT CASE(variableType)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(globalDerivativeIndex)
          CASE(NO_GLOBAL_DERIV)
            IF(componentIdx==1) THEN
              !calculate u
              analyticValue=EXP((x(1)-x(2))/10.0_DP)+EXP((x(3)-x(1))/10.0_DP)
            ELSE IF(componentIdx==2) THEN
              !calculate v
              analyticValue=EXP((x(1)-x(2))/10.0_DP)+EXP((x(2)-x(3))/10.0_DP)
            ELSE IF(componentIdx==3) THEN
              !calculate w
              analyticValue=EXP((x(3)-x(1))/10.0_DP)+EXP((x(2)-x(3))/10.0_DP)
            ELSE IF(componentIdx==4) THEN
              !calculate p
              analyticValue=2.0_DP*muParam/10.0_DP*(EXP((x(1)-x(2))/10.0_DP)-EXP((x(3)-x(1))/10.0_DP))
            ELSE
              CALL FlagError("Not implemented.",err,error,*999)
            ENDIF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(globalDerivativeIndex)
          CASE(NO_GLOBAL_DERIV)
            IF(componentIdx==1) THEN
              !calculate u
              analyticValue=0.0_DP
            ELSE IF(componentIdx==2) THEN
              !calculate v
              analyticValue=-2.0_DP*muParam*(2.0_DP*EXP(x(1)-x(2))+EXP(x(2)-x(3)))
            ELSE IF(componentIdx==3) THEN
              !calculate w
              analyticValue=-2.0_DP*muParam*(2.0_DP*EXP(x(3)-x(1))+EXP(x(2)-x(3)))
            ELSE IF(componentIdx==4) THEN
              !calculate p
              analyticValue=0.0_DP
            ELSE
              CALL FlagError("Not implemented.",err,error,*999)
            ENDIF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The variable type of "//TRIM(NumberToVString(variableType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        localError="The number of components does not correspond to the number of dimensions."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE(EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_3)
      IF(numberOfDimensions==3.AND.numberOfComponents==4) THEN
        !Sine and cosine functions
        SELECT CASE(variableType)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(globalDerivativeIndex)
          CASE(NO_GLOBAL_DERIV)
            IF(componentIdx==1) THEN
              !calculate u
              analyticValue=sin(2.0_DP*PI*x(1)/10.0_DP)*sin(2.0_DP*PI*x(2)/10.0_DP)*sin(2.0_DP*PI*x(3)/10.0_DP)
            ELSE IF(componentIdx==2) THEN
              !calculate v
              analyticValue=2.0_DP*cos(2.0_DP*PI*x(1)/10.0_DP)*sin(2.0_DP*PI*x(3)/10.0_DP)*cos(2.0_DP*PI*x(2)/10.0_DP)
            ELSE IF(componentIdx==3) THEN
              !calculate w
              analyticValue=-cos(2.0_DP*PI*x(1)/10.0_DP)*sin(2.0_DP*PI*x(2)/10.0_DP)*cos(2.0_DP*PI*x(3)/10.0_DP)
            ELSE IF(componentIdx==4) THEN
              !calculate p
              analyticValue=6.0_DP*muParam*PI/10.0_DP*sin(2.0_DP*PI*x(2)/10.0_DP)*sin(2.0_DP*PI*x(3)/10.0_DP)* &
                & cos(2.0_DP*PI*x(1)/10.0_DP)
            ELSE
              CALL FlagError("Not implemented.",err,error,*999)
            ENDIF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(globalDerivativeIndex)
          CASE(NO_GLOBAL_DERIV)
            IF(componentIdx==1) THEN
              !calculate u
              analyticValue=0.0_DP
            ELSE IF(componentIdx==2) THEN
              !calculate v
              analyticValue=36*muParam*PI**2/10.0_DP**2*cos(2.0_DP*PI*x(2)/10.0_DP)*sin(2.0_DP*PI*x(3)/10.0_DP)* &
                & cos(2.0_DP*PI*x(1)/10.0_DP)
            ELSE IF(componentIdx==3) THEN
              !calculate w
              analyticValue=0.0_DP
            ELSE IF(componentIdx==4) THEN
              !calculate p
              analyticValue=0.0_DP
            ELSE
              CALL FlagError("Not implemented.",err,error,*999)
            ENDIF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The variable type of "//TRIM(NumberToVString(variableType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        localError="The number of components does not correspond to the number of dimensions."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE(EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_4)
      IF(numberOfDimensions==3.AND.numberOfComponents==4) THEN
        !Reduced Taylor-Green solution for Stokes
        CALL FlagError("Not implemented.",err,error,*999)
      ENDIF
    CASE(EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_5)
      IF(numberOfDimensions==3.AND.numberOfComponents==4) THEN
        !Stokes-Taylor-Green dynamic
        SELECT CASE(variableType)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(globalDerivativeIndex)
          CASE(NO_GLOBAL_DERIV)
            IF(componentIdx==1) THEN
              !calculate u
              analyticValue=x(2)*exp(-(2.0_DP*muParam/rhoParam*currentTime))
            ELSE IF(componentIdx==2) THEN
              !calculate v
              analyticValue=x(1)*exp(-(2.0_DP*muParam/rhoParam*currentTime))
            ELSE IF(componentIdx==3) THEN
              !calculate v
              analyticValue=0.0_DP
            ELSE IF(componentIdx==4) THEN
              !calculate p
              analyticValue=2.0_DP*x(2)*muParam*exp(-(2.0_DP*muParam/rhoParam*currentTime))*x(1)
            ELSE
              CALL FlagError("Not implemented.",err,error,*999)
            ENDIF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(globalDerivativeIndex)
          CASE(NO_GLOBAL_DERIV)
            IF(componentIdx==1) THEN
              !calculate u
              analyticValue=0.0_DP
            ELSE IF(componentIdx==2) THEN
              !calculate v
              analyticValue=0.0_DP
            ELSE IF(componentIdx==3) THEN
              !calculate p
              analyticValue=0.0_DP
            ELSE IF(componentIdx==4) THEN
              !calculate p
              analyticValue=0.0_DP
            ELSE
              CALL FlagError("Not implemented.",err,error,*999)
            ENDIF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The variable type of "//TRIM(NumberToVString(variableType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        localError="The number of components does not correspond to the number of dimensions."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE DEFAULT
      localError="The analytic function type of "//TRIM(NumberToVString(analyticFunctionType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Stokes_AnalyticFunctions")
    RETURN
999 ERRORSEXITS("Stokes_AnalyticFunctions",err,error)
    RETURN 1

  END SUBROUTINE Stokes_AnalyticFunctions

  !
  !================================================================================================================================
  !

END MODULE StokesEquationsRoutines
