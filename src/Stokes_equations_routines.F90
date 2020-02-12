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
  USE Constants
  USE ControlLoopRoutines
  USE ControlLoopAccessRoutines
  USE DistributedMatrixVector
  USE DistributedMatrixVectorAccessRoutines
  USE DomainMappings
  USE EquationsRoutines
  USE EquationsAccessRoutines
  USE EquationsMappingRoutines
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

  PUBLIC  STOKES_EQUATION_ANALYTIC_FUNCTIONS

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
      & EQUATIONS_SET_ALE_STOKES_SUBTYPE, &
      & EQUATIONS_SET_PGM_STOKES_SUBTYPE)
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
      & EQUATIONS_SET_ALE_STOKES_SUBTYPE, &
      & EQUATIONS_SET_PGM_STOKES_SUBTYPE)
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
    INTEGER(INTG):: dependentNumberOfVariables,dependentNumberOfComponents
    INTEGER(INTG):: independentFieldNumberOfVariables,independentNumberOfComponents
    INTEGER(INTG):: numberOfDimensions,geometricComponentNumber
    INTEGER(INTG):: materialsNumberOfVariables,materialsNumberOfComponents,componentIdx
    INTEGER(INTG) :: geometricScalingType,geometricMeshComponent
    TYPE(DecompositionType), POINTER :: geometricDecomposition
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsSetMaterialsType), POINTER :: equationsMaterials
    TYPE(VARYING_STRING) :: localError

    ENTERS("Stokes_EquationsSetSetup",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
    
    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_STATIC_STOKES_SUBTYPE, &
      & EQUATIONS_SET_TRANSIENT_STOKES_SUBTYPE, &
      & EQUATIONS_SET_LAPLACE_STOKES_SUBTYPE, &
      & EQUATIONS_SET_ALE_STOKES_SUBTYPE, &
      & EQUATIONS_SET_PGM_STOKES_SUBTYPE)
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
      CALL Field_NumberOfComponentsGetgeometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
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
          CALL Field_ComponentMeshComponentGetgeometricField,FIELD_U_VARIABLE_TYPE,1,geometricMeshComponent,err,error,*999)
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
      CASE(EQUATIONS_SET_ALE_STOKES_SUBTYPE,EQUATIONS_SET_PGM_STOKES_SUBTYPE)
        NULLIFY(geometricField)
        CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
        CALL Field_NumberOfComponentsGetgeometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
        NULLIFY(equationsIndependent)
        CALL EquationsSet_IndependentGet(equationsSet,equationsIndependent,err,error,*999)
        independentFieldNumberOfVariables=1
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
            CALL Field_NumberOfVariablesSetAndLock(equationsIndependent%independentField,independentFieldNumberOfVariables, &
              & err,error,*999)
            CALL Field_VariableTypesSetAndLock(equationsIndependent%independentField,[FIELD_U_VARIABLE_TYPE],err,error,*999)
            CALL Field_DimensionSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DataTypeSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE, &
              & err,error,*999)
            CALL Field_NumberOfComponentsGet(equationsSet%GEOMETRY%geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions, &
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
          CALL Field_NumberOfComponentsGetgeometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
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
      CASE(EQUATIONS_SET_ALE_STOKES_SUBTYPE,EQUATIONS_SET_PGM_STOKES_SUBTYPE)
        !Do nothing
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
        CASE(EQUATIONS_SET_TRANSIENT_STOKES_SUBTYPE,EQUATIONS_SET_ALE_STOKES_SUBTYPE,EQUATIONS_SET_PGM_STOKES_SUBTYPE)
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
          CASE(EQUATIONS_SET_TRANSIENT_STOKES_SUBTYPE,EQUATIONS_SET_ALE_STOKES_SUBTYPE,EQUATIONS_SET_PGM_STOKES_SUBTYPE)
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
          CASE(EQUATIONS_SET_TRANSIENT_STOKES_SUBTYPE,EQUATIONS_SET_ALE_STOKES_SUBTYPE,EQUATIONS_SET_PGM_STOKES_SUBTYPE)
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
      & PROBLEM_ALE_STOKES_SUBTYPE, &
      & PROBLEM_PGM_STOKES_SUBTYPE)
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
    TYPE(ControlLoopType), POINTER :: controlLoop,controLoopRoot
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
      & PROBLEM_PGM_STOKES_SUBTYPE, &
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
        IF(esSpecification(3)/=PROBLEM_STATIC_STOKES_SUBTYPE.AND. &
          & esSpecification(3)/=PROBLEM_LAPLACE_STOKES_SUBTYPE) THEN
          CALL ControlLoop_TypeSet(controlLoop,CONTROL_TIME_LOOP_TYPE,err,error,*999)
        ENDIF
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Finish the control loops
        NULLIFY(controlLoopRoot)
        CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
        NULLIFY(controlLoop)
        CALL ControlLoop_Get(controLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
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
      CALL ControlLoop_Get(controLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)
        !Start the solvers creation
        NULLIFY(solvers)
        CALL Solvers_CreateStart(controlLoop,solvers,err,error,*999)
        SELECT CASE(esSpecification(3))
        CASE(PROBLEM_STATIC_STOKES_SUBTYPE, &
          & PROBLEM_LAPLACE_STOKES_SUBTYPE)
          CALL Solvers_NumberOfSolversSet(solvers,1,err,error,*999)
          !Set the solver to be a linear solver
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          CALL Solver_TypeSet(solver,SOLVER_LINEAR_TYPE,err,error,*999)
          !Set solver defaults
          CALL Solver_LibraryTypeSet(solver,SOLVER_PETSC_LIBRARY,err,error,*999)
        CASE(PROBLEM_TRANSIENT_STOKES_SUBTYPE, &
          & PROBLEM_PGM_STOKES_SUBTYPE)
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
      CALL ControlLoop_Get(controLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
      NULLIFY(solvers)
      CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)
        SELECT CASE(esSpecification(3))
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
        CASE(PROBLEM_TRANSIENT_STOKES_SUBTYPE, &
          & PROBLEM_PGM_STOKES_SUBTYPE)
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
          CALL SolverEquations_LinearityTypeSet(solverEquations,solverEquations_LINEAR,err,error,*999)
          CALL SolverEquations_TimeDependenceTypeSet(solverEquations,solverEquations_STATIC,err,error,*999)
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
        IF(esSpecification(3)==PROBLEM_ALE_STOKES_SUBTYPE) THEN
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
    INTEGER(INTG) variableType,gaussPointIdx,rowComponentIdx,rowElelemntDOFIdx,rowXiIdx,rowElementParameterIdx,columnComponentIdx,columnElementDOFIdx,columnXiIdx,columnElementParameterIdx,meshComponent1,meshComponent2, columnElementDOFIdx_max, rowElelemntDOFIdx_max, columnElementDOFIdx_min, rowElelemntDOFIdx_min
    REAL(DP) :: JGW,SUM,dXidX(3,3),rowsPhi,colsPhi,muParam,rhoParam,rowsdPhidXi(3),colsdPhidXi(3)
    REAL(DP) :: AG_MATRIX(256,256) ! "A" Matrix ("G"radient part) - maximum size allocated
    REAL(DP) :: AL_MATRIX(256,256) ! "A" Matrix ("L"aplace part) - maximum size allocated
    REAL(DP) :: BT_MATRIX(256,256) ! "B" "T"ranspose Matrix - maximum size allocated
    REAL(DP) :: MT_MATRIX(256,256) ! "M"ass "T"ime Matrix - maximum size allocated
    REAL(DP) :: CT_MATRIX(256,256) ! "C"onvective "T"erm Matrix - maximum size allocated
    REAL(DP) :: ALE_MATRIX(256,256) ! "A"rbitrary "L"agrangian "E"ulerian Matrix - maximum size allocated
    REAL(DP) :: RH_VECTOR(256) ! "R"ight "H"and vector - maximum size allocated
    REAL(DP) :: wValue(3)
    REAL(DP)::  X(3)
    LOGICAL :: updateStiffnessMatrix, updateDampingMatrix,updateRHSVector
    TYPE(BasisType), POINTER :: dependentBasis,dependentBasis1,dependentBasis2,geometricBasis,independentBasis
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatrixType), POINTER :: stiffnessMatrix, dampingMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField,geometricField,materialsField,independentField
    TYPE(FieldVariableType), POINTER :: fieldVariable
    TYPE(QuadratureSchemeType), POINTER :: quadratureScheme,quadratureScheme1,quadratureScheme2
    TYPE(VARYING_STRING) :: localError
    INTEGER:: xv,out

!\todo: Reduce number of variables and parameters

    ENTERS("Stokes_FiniteElementCalculate",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)

    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_STATIC_STOKES_SUBTYPE, &
      & EQUATIONS_SET_LAPLACE_STOKES_SUBTYPE, &
      & EQUATIONS_SET_TRANSIENT_STOKES_SUBTYPE, &
      & EQUATIONS_SET_ALE_STOKES_SUBTYPE, &
      & EQUATIONS_SET_PGM_STOKES_SUBTYPE)
    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is not valid for a Stokes fluid type of a fluid mechanics equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
      
    out=0
    AG_MATRIX=0.0_DP
    AL_MATRIX=0.0_DP
    BT_MATRIX=0.0_DP
    MT_MATRIX=0.0_DP
    CT_MATRIX=0.0_DP
    ALE_MATRIX=0.0_DP
    RH_VECTOR=0.0_DP
    X=0.0_DP
!     L=10.0_DP

    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*9999)
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
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
    NULLIFY(rhsVector)
    CALL EquationsMatricesVector_RHSVectorGet(vectorMatrices,rhsVector,err,error,*999)
    updateRHS=rhsVector%updateVector
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
      CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,1,stiffessMatrix,err,error,*999)
      NULLIFY(dynamicMapping)
      NULLIFY(dynamicMatrices)
      NULLIFY(dampingMatrix)
      updateStiffnes=stiffnessMatrix%updateMatrix
      updateDamping=.FALSE.
    CASE(EQUATIONS_SET_TRANSIENT_STOKES_SUBTYPE, &
      & EQUATIONS_SET_ALE_STOKES_SUBTYPE, &
      & EQUATIONS_SET_PGM_STOKES_SUBTYPE)
      NULLIFY(linearMapping)
      NULLIFY(linearMatrices)
      NULLIFY(dynamicMapping)
      CALL EquationsMappingVector_DynamicMappingGet(vectorMapping,dynamicMapping,err,error,*999)
      NULLIFY(colsVariable)
      CALL EquationsMappingDynamic_DynamicVariableGet(dynamicMapping,colsVariable,err,error,*999)
      NULLIFY(dynamicMatrices)
      CALL EquationsMatricesVector_DynamicMatricesGet(vectorMatrices,dynamicMatrices,err,error,*999)
      NULLIFY(stiffnessMatrix)
      CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,1,stiffessMatrix,err,error,*999)
      NULLIFY(dampingMatrix)
      CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,2,dampingMatrix,err,error,*999)
      updateStiffnes=stiffnessMatrix%updateMatrix
      updateDamping=dampingMatrix%updateMatrix
    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is not valid for a Stokes fluid type of a fluid mechanics equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    updateMatrices=(updateStiffness.OR.updateDamping)
    update=(updateMatrices.OR.updateRHS)

    IF(update) THEN
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      NULLIFY(dependentField)
      CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
      NULLIFY(materialsField)
      CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)
      NULLIFY(independentField)
      IF(esSpecification(3)==EQUATIONS_SET_ALE_STOKES_SUBTYPE.OR. &
        & esSpecification(3)==EQUATIONS_SET_PGM_STOKES_SUBTYPE) THEN
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
      CALL DomainElements_BasisGet(geometricDomainElements,elementNumber,geometricBasis,err,error,*999)
      CALL Basis_NumberOfXiGet(geometricBasis,numberOfXi,err,error,*999)
      
      NULLIFY(dependentDecomposition)
      CALL Field_DecompositionGet(dependentField,dependentDecomposition,err,error,*999)
      NULLIFY(colsDomain)
      CALL Decomposition_DomainGet(dependentDecomposition,0,colsDomain,err,error,*999)
      NULLIFY(colsDomainTopology)
      CALL Domain_DomainTopologyGet(colsDomain,colsDomainTopology,err,error,*999)
      NULLIFY(colsDomainElements)
      CALL DomainTopology_DomainElementsGet(colsDomainTopology,colsDomainElements,err,error,*999)
      NULLIFY(colsBasis)
      CALL DomainElements_BasisGet(colsDomainElements,elementNumber,colsBasis,err,error,*999)
      
      NULLIFY(rowsVariable)
      CALL EquationsMappingLHS_LHSVariableGet(lhsMapping,rowsVariable,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(rowsVariable,numberOfRowsComponents,err,error,*999)
      
      CALL FieldVariable_VariableTypeGet(colsVariable,colsVariableType,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(colsVariable,numberOfColsComponents,err,error,*999)
      
      NULLIFY(geometricQuadratureScheme)
      CALL Basis_QuadratureSchemeGet(geometricBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,geometricQuadratureScheme,err,error,*999)
      NULLIFY(colsQuadratureScheme)
      CALL Basis_QuadratureSchemeGet(colsBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,colsQuadratureScheme,err,error,*999)
      CALL BasisQuadrature_NumberOfGaussGet(colsQuadratureScheme,numberOfGauss,err,error,*999)
      
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
      
      NULLIFY(materialsInterpParamters)
      CALL EquationsInterpolation_MaterialsParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & materialsInterpParameters,err,error,*999)
      NULLIFY(materialsInterpPoint)
      CALL EquationsInterpolation_MaterialsPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,materialsInterpPoint, &
        & err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,materialsInterpParameters,err,error,*999)
      
      NULLIFY(independentInterpParameters)
      NULLIFY(independentInterpPoint)
      IF(esSpecification(3)==EQUATIONS_SET_ALE_STOKES_SUBTYPE.OR. &
        & esSpecification(3)==EQUATIONS_SET_PGM_STOKES_SUBTYPE) THEN
        CALL EquationsInterpolation_IndependentParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
          & independentInterpParameters,err,error,*999)
        CALL EquationsInterpolation_IndependentPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
          & independentInterpPoint,err,error,*999)
        CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,independentInterpParameters, &
          & err,error,*999)
      ENDIF
      

      !Start looping over Gauss points
      DO gaussPointIdx=1,numberOfGauss
        CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,geometricInterpPoint, &
          & err,error,*999)
        CALL Field_InterpolatedPointMetricsCalculate(numberOfXi,geometricInterpPointMetrics,err,error,*999)
        CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,materialsInterpPoint, &
          & err,error,*999)
        IF(equationsSet%specification(3)==EQUATIONS_SET_ALE_STOKES_SUBTYPE.OR. &
          & equationsSet%specification(3)==EQUATIONS_SET_PGM_STOKES_SUBTYPE) THEN
          CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,independentInterpPoint, &
            & err,error,*999)
          wValue(1)=independentInterpPoint%values(1,NO_PART_DERIV)
          wValue(2)=independentInterpPoint%values(2,NO_PART_DERIV)
          IF(numberOfDimensions==3) THEN
            wValue(3)=independentInterpPoint%values(3,NO_PART_DERIV)
          END IF
        ELSE
          wValue=0.0_DP
        END IF
        !Define muParam, viscosity=1
        muParam=materialsInterpPoint%values(1,NO_PART_DERIV)
        !Define rhoParam, density=2
        rhoParam=materialsInterpPointr%values(2,NO_PART_DERIV)
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
 
        !Loop over field components
        rowElementDOFIdx=0          
        DO rowComponentIdx=1,numberOfRowsComponents
          NULLIFY(rowsDomain)
          CALL FieldVariable_ComponentDomainGet(rowsVariable,rowComponentIdx,rowsDomain,err,error,*999)
          NULLIFY(rowsDomainTopology)
          CALL Domain_DomainTopologyGet(rowsDomain,rowsDomainTopology,err,error,*999)
          NULLIFY(rowsDomainElements)
          CALL DomainTopology_DomainElementsGet(rowsDomainTopology,rowsDomainElements,err,error,*999)
          NULLIFY(rowsBasis)
          CALL DomainElements_ElementBasisGet(rowsDomainElements,elementNumber,rowsBasis,err,error,*999)
          CALL Basis_QuadratureSchemeGet(rowsBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,rowsQuadratureScheme,err,error,*999)
          CALL Basis_NumberOfElementParametersGet(rowsBasis,numberOfRowsElementParameters,err,erorr,*999)
          !Loop over element rows
          DO rowElementParameterIdx=1,numberOfRowsElementParameters
            rowElementDOFIdx=rowElementDOFIdx+1
            CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowsQuadratureScheme,rowElementParameterIdx,NO_PART_DERIV, &
              & rowsPhi,err,error,*999)
            DO rowXiIdx=1,numberOfXi
              CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowsQuadratureScheme,rowElementParameterIdx, &
                & PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(rowXiIdx1),gaussPointIdx,rowsdPhidXi(rowsXiIdx),err,error,*999)
            ENDDO !rowXiiIdx
            IF(updateMatrices) THEN
              !Loop over element columns
              columnElementDOFIdx=0
              DO columnComponentIdx=1,numberOfColsComponents
                NULLIFY(colsDomain)
                CALL FieldVariable_ComponentDomainGet(colsVariable,columnComponentIdx,colsDomain,err,error,*999)
                NULLIFY(colsDomainTopology)
                CALL Domain_DomainTopologyGet(colsDomain,colsDomainTopology,err,error,*999)
                NULLIFY(colsDomainElements)
                CALL DomainTopology_DomainElementsGet(colsDomainTopology,colsDomainElements,err,error,*999)
                NULLIFY(colsBasis)
                CALL DomainElements_ElementBasisGet(colsDomainElements,elementNumber,colsBasis,err,error,*999)
                CALL Basis_QuadratureSchemeGet(colsBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,colsQuadratureScheme,err,error,*999)
                CALL Basis_NumberOfElementParametersGet(colsBasis,numberOfColsElementParameters,err,error,*999)
                DO columnElementParameterIdx=1,numberOfColsElementParameters
                  columnElementDOFIdx=columnElementDOFIdx+1
                  CALL BasisQuadratureScheme_GaussBasisFunctionGet(colsQuadratureScheme,columnElementParameterIdx,NO_PART_DERIV, &
                    & colsPhi,err,error,*999)
                  DO columnXiIdx1=1,numberOfXi
                    CALL BasisQuadratureScheme_GaussBasisFunctionGet(colsQuadratureScheme,columnElementParameterIdx, &
                      & PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(columnXiIdx),gaussPointIdx,colsdPhidXi(columnXiIdx), &
                      & err,error,*999)
                  ENDDO !columnXiIdx
                  IF(updateStiffness) THEN
                    !LAPLACE TYPE
                    IF(columnComponentIdx==rowComponentIdx) THEN
                      sum=0.0_DP
                      !Calculate SUM
                      DO componentIdx=1,numberOfColsComponents
                        DO rowXiIdx=1,numberOfXi
                          DO columnXiIdx=1,numberOfXi
                            sum=sum+muParam*rowsdPhidXi(rowXiIdx)*dXidX(rowXiIdx,componentIdx)* &
                              & colsdPhidXi(columnXiIdx)*dXidX(columnXiIdx,componentIdx)
                          ENDDO !columnXiIdx
                        ENDDO !rowXiIdx
                      ENDDO !x
                      !Calculate MATRIX
                      AL_MATRIX(rowElelemntDOFIdx,columnElementDOFIdx)=AL_MATRIX(rowElelemntDOFIdx,columnElementDOFIdx)+ &
                        & sum*jacobianGaussWeight
                    ENDIF
                    IF(equationsSet%specification(3)/=EQUATIONS_SET_LAPLACE_STOKES_SUBTYPE) THEN
                      IF(columnComponentIdx<numberOfColsComponents) THEN
                        sum=0.0_DP
                        !Calculate sum
                        DO rowXiIdx=1,numberOfXi
                          DO columnXiIdx=1,dependentBasis2%numberOfXi
                            !note rowComponentIdx/columnComponentIdx derivative in dXidX
                            sum=sum+muParam*rowsdPhidXi(columnXiIdx)*dXidX(columnXiIdx,columnComponentIdx)* &
                              & colsdPhidXi(rowXiIdx)*dXidX(rowXiIdx,rowComponentIdx)
                          ENDDO !columnXiIdx
                        ENDDO !rowXiIdx
                        !Calculate MATRIX
                        AG_MATRIX(rowElelemntDOFIdx,columnElementDOFIdx)= &
                          & AG_MATRIX(rowElelemntDOFIdx,columnElementDOFIdx)+ &
                          & sum*jacobianGaussWeight
                      ENDIF
                    ENDIF
                    IF(equationsSet%specification(3)==EQUATIONS_SET_ALE_STOKES_SUBTYPE.OR. &
                      & equationsSet%specification(3)==EQUATIONS_SET_PGM_STOKES_SUBTYPE) THEN
                      IF(columnComponentIdx==rowComponentIdx) THEN
                        sum=0.0_DP
                        !Calculate SUM
                        DO componentIdx=1,numberOfColsComponents
                          DO columnXiIdx=1,numberOfXi
                            sum=sum-rhoParam*wValue(componentIdx)*rowsPhi*colsdPhidXi(columnXiIdx)*dXidX(columnXiIdx,componentIdx)
                          ENDDO !columnXiIdx
                        ENDDO !rowXiIdx
                        !Calculate MATRIX
                        ALE_MATRIX(rowElelemntDOFIdx,columnElementDOFIdx)= &
                          & ALE_MATRIX(rowElelemntDOFIdx,columnElementDOFIdx)+ &
                          & sum*jacobianGaussWeight
                      ENDIF
                    ENDIF
                    !Calculate pressure contribution (B transpose type)
                    !LAPLACE TYPE
                    IF(columnComponentIdx==numberOfColsComponents) THEN
                      sum=0.0_DP
                      !Calculate SUM
                      DO columnXiIdx=1,dependentBasis1%numberOfXi
                        sum=sum-rowsdPhidXi(columnXiIdx)*dXidX(columnXiIdx,rowComponentIdx)*colsPhi
                      ENDDO !columnXiIdx
                      !Calculate MATRIX
                      BT_MATRIX(rowElelemntDOFIdx,columnElementDOFIdx)= &
                        & BT_MATRIX(rowElelemntDOFIdx,columnElementDOFIdx)+ &
                        & sum*jacobianGaussWeight
                    ENDIF
                  ENDIF
                  IF(updateDamping) THEN
                    IF(columnComponentIdx==rowComponentIdx) THEN
                      dampingMatrix%elementMatrix%matrix(rowElelemntDOFIdx,columnElementDOFIdx)= &
                        & dampingMatrix%elementMatrix%matrix(rowElelemntDOFIdx,columnElementDOFIdx)+ &
                        & rowParam*rowsPhi*colsPhi*jacobianGaussWeight
                    ENDIF
                  ENDIF
                ENDDO !columnElementParameterIdx
              ENDDO !columnComponentIdx
            ENDDO !rowElementParameterIdx
          ENDDO !rowComponentIdx
          
                    IF(updateStiffnessMatrix.OR.updateDampigaussPointIdxMatrix) THEN
                      !Loop over element columns
                      DO columnComponentIdx=1,(fieldVariable%numberOfComponents)

                        meshComponent2=fieldVariable%COMPONENTS(columnComponentIdx)%meshComponentNumber
                        dependentBasis2=>dependentField%DECOMPOSITION%DOMAIN(meshComponent2)%ptr% &
                          & TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS
                        quadratureScheme2=>dependentBasis2%QUADRATURE%quadratureSchemeMap &
                          & (BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
                        ! JGW=equations%interpolation%geometricInterpPointMetrics%jacobian*quadratureScheme2%&
                        ! &gaussWeights(gaussPointIdx)
                        DO columnElementParameterIdx=1,dependentBasis2%numberOfElementParameters
                          columnElementDOFIdx=columnElementDOFIdx+1
                        !Calculate some variables used later on
                          DO columnXiIdx=1,dependentBasis2%numberOfXi
                            DO rowXiIdx=1,dependentBasis1%numberOfXi
                              dXidX(rowXiIdx,columnXiIdx)=equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr% &
                                & dXidX(rowXiIdx,columnXiIdx)
                            END DO
                            rowsdPhidXi(columnXiIdx)=quadratureScheme1%gaussBasisFunctiocolumnElementParameterIdx(rowElementParameterIdx,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(columnXiIdx),gaussPointIdx)
                            colsdPhidXi(columnXiIdx)=quadratureScheme2%gaussBasisFunctions(columnElementParameterIdx,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(columnXiIdx),gaussPointIdx)
                          END DO !columnXiIdx
                          rowsPhi=quadratureScheme1%gaussBasisFunctions(rowElementParameterIdx,NO_PART_DERIV,gaussPointIdx)
                          colsPhi=quadratureScheme2%gaussBasisFunctions(columnElementParameterIdx,NO_PART_DERIV,gaussPointIdx)
                        !                         DO rowXiIdx=1,dependentBasis1%numberOfXi
                        !                           DO columnXiIdx=1,dependentBasis2%numberOfXi
                        !                             SUM=SUM-muParam*DrowsPhiS_DXI(rowXiIdx)*DPHINSS_DXI(columnXiIdx)*equations%interpolation%geometricInterpPointMetrics%GU(rowXiIdx,columnXiIdx)
                        !                           ENDDO !columnXiIdx
                        !                         ENDDO !rowXiIdx
                          IF(updateStiffnessMatrix) THEN

                            !LAPLACE TYPE
                            IF(columnComponentIdx==rowComponentIdx) THEN
                              SUM=0.0_DP
                              !Calculate SUM
                              DO xv=1,dependentBasis1%numberOfXi
                                DO rowXiIdx=1,dependentBasis1%numberOfXi
                                  DO columnXiIdx=1,dependentBasis2%numberOfXi
                                    SUM=SUM+muParam*colsdPhidXi(columnXiIdx)*dXidX(columnXiIdx,xv)*rowsdPhidXi(rowXiIdx)*dXidX(rowXiIdx,xv)
                                  ENDDO !columnXiIdx
                                ENDDO !rowXiIdx
                              ENDDO !x
                              !Calculate MATRIX
                              AL_MATRIX(rowElelemntDOFIdx,columnElementDOFIdx)=AL_MATRIX(rowElelemntDOFIdx,columnElementDOFIdx)+SUM*JGW
                            END IF

                          END IF
                          !Calculate standard matrix (gradient transpose type)
                          IF(updateStiffnessMatrix) THEN

                            IF(equationsSet%specification(3)/=EQUATIONS_SET_LAPLACE_STOKES_SUBTYPE) THEN
                              IF(columnComponentIdx<fieldVariable%numberOfComponents) THEN
                                SUM=0.0_DP
                                !Calculate SUM
                                DO rowXiIdx=1,dependentBasis1%numberOfXi
                                  DO columnXiIdx=1,dependentBasis2%numberOfXi
                                    !note rowComponentIdx/columnComponentIdx derivative in dXidX
                                    SUM=SUM+muParam*colsdPhidXi(rowXiIdx)*dXidX(rowXiIdx,rowComponentIdx)*rowsdPhidXi(columnXiIdx)*dXidX(columnXiIdx,columnComponentIdx)
                                  ENDDO !columnXiIdx
                                ENDDO !rowXiIdx
                                !Calculate MATRIX
                                AG_MATRIX(rowElelemntDOFIdx,columnElementDOFIdx)=AG_MATRIX(rowElelemntDOFIdx,columnElementDOFIdx)+SUM*JGW
                              END IF
                            END IF

                          END IF
                          !Calculate ALE matric contribution
                          IF(updateStiffnessMatrix) THEN

                            IF(equationsSet%specification(3)==EQUATIONS_SET_ALE_STOKES_SUBTYPE.OR. &
                              & equationsSet%specification(3)==EQUATIONS_SET_PGM_STOKES_SUBTYPE) THEN
                              IF(columnComponentIdx==rowComponentIdx) THEN
                                SUM=0.0_DP
                                !Calculate SUM
                                DO rowXiIdx=1,dependentBasis1%numberOfXi
                                  DO columnXiIdx=1,dependentBasis1%numberOfXi
                                    SUM=SUM-rhoParam*wValue(rowXiIdx)*colsdPhidXi(columnXiIdx)*dXidX(columnXiIdx,rowXiIdx)*rowsPhi
                                  ENDDO !columnXiIdx
                                ENDDO !rowXiIdx
                                !Calculate MATRIX
                                ALE_MATRIX(rowElelemntDOFIdx,columnElementDOFIdx)=ALE_MATRIX(rowElelemntDOFIdx,columnElementDOFIdx)+SUM*JGW
                              END IF
                            END IF

                          END IF
                          !Calculate pressure contribution (B transpose type)
                          IF(updateStiffnessMatrix) THEN

                            !LAPLACE TYPE
                            IF(columnComponentIdx==fieldVariable%numberOfComponents) THEN
                              SUM=0.0_DP
                              !Calculate SUM
                              DO columnXiIdx=1,dependentBasis1%numberOfXi
                                SUM=SUM-colsPhi*rowsdPhidXi(columnXiIdx)*dXidX(columnXiIdx,rowComponentIdx)
                              ENDDO !columnXiIdx
                              !Calculate MATRIX
                              BT_MATRIX(rowElelemntDOFIdx,columnElementDOFIdx)=BT_MATRIX(rowElelemntDOFIdx,columnElementDOFIdx)+SUM*JGW
                            END IF

                          END IF
                          !Calculate mass matrix if needed
                          IF(equationsSet%specification(3)==EQUATIONS_SET_TRANSIENT_STOKES_SUBTYPE.OR. &
                            & equationsSet%specification(3)==EQUATIONS_SET_ALE_STOKES_SUBTYPE.OR. &
                            & equationsSet%specification(3)==EQUATIONS_SET_PGM_STOKES_SUBTYPE) THEN
                            IF(updateDampingMatrix) THEN
                              IF(columnComponentIdx==rowComponentIdx) THEN
                                SUM=0.0_DP
                                !Calculate SUM
                                SUM=rhoParam*rowsPhi*colsPhi
                                !Calculate MATRIX
                                MT_MATRIX(rowElelemntDOFIdx,columnElementDOFIdx)=MT_MATRIX(rowElelemntDOFIdx,columnElementDOFIdx)+SUM*JGW
                              END IF
                            END IF
                          END IF
                        ENDDO !columnElementParameterIdx
                      ENDDO !columnComponentIdx
                    ENDIF
                  ENDDO !rowElementParameterIdx
                ENDDO !rowComponentIdx

                !Calculate analytic RHS
                IF(ASSOCIATED(equationsSet%ANALYTIC)) THEN
                  IF(equationsAnalytic%analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_1.OR. &
                    & equationsAnalytic%analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_2.OR. &
                    & equationsAnalytic%analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_3.OR. &
                    & equationsAnalytic%analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_4.OR. &
                    & equationsAnalytic%analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_5.OR. &
                    & equationsAnalytic%analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_1.OR. &
                    & equationsAnalytic%analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_2.OR. &
                    & equationsAnalytic%analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_3.OR. &
                    & equationsAnalytic%analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_4.OR. &
                    & equationsAnalytic%analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_5) THEN

                    rowElelemntDOFIdx=0
                    DO rowComponentIdx=1,(fieldVariable%numberOfComponents-1)
                      meshComponent1=fieldVariable%COMPONENTS(rowComponentIdx)%meshComponentNumber
                      dependentBasis1=>dependentField%DECOMPOSITION%DOMAIN(meshComponent1)%ptr% &
                        & TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS
                      quadratureScheme1=>dependentBasis1%QUADRATURE%quadratureSchemeMap(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
                      JGW=equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%jacobian* &
                        & quadratureScheme1%gaussWeights(gaussPointIdx)
                      DO rowElementParameterIdx=1,dependentBasis1%numberOfElementParameters
                        rowElelemntDOFIdx=rowElelemntDOFIdx+1
                        rowsPhi=quadratureScheme1%gaussBasisFunctions(rowElementParameterIdx,NO_PART_DERIV,gaussPointIdx)
                        !note rowComponentIdx value derivative
                        SUM=0.0_DP
                        X(1) = equations%interpolation%geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(1,1)
                        X(2) = equations%interpolation%geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(2,1)
                        IF(dependentBasis1%numberOfXi==3) THEN
                          X(3) = equations%interpolation%geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(3,1)
                        END IF
                        IF(equationsAnalytic%analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_1) THEN
                          IF(rowComponentIdx==1) THEN
                            !Calculate SUM
                            SUM=0.0_DP
                          ELSE IF(rowComponentIdx==2) THEN
                            !Calculate SUM
                            SUM=rowsPhi*(-2.0_DP*muParam/10.0_DP**2)
                          ENDIF
                        ELSE IF(equationsAnalytic%analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_2) THEN
                          IF(rowComponentIdx==1) THEN
                            !Calculate SUM
                            SUM=0.0_DP
                          ELSE IF(rowComponentIdx==2) THEN
                            !Calculate SUM
                            SUM=rowsPhi*(-4.0_DP*muParam/100.0_DP*EXP((X(1)-X(2))/10.0_DP))
                          ENDIF
                        ELSE IF(equationsAnalytic%analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_3) THEN
                          IF(rowComponentIdx==1) THEN
                            !Calculate SUM
                            SUM=0.0_DP
                          ELSE IF(rowComponentIdx==2) THEN
                            !Calculate SUM
                            SUM=rowsPhi*(16.0_DP*muParam*PI*PI/100.0_DP*COS(2.0_DP*PI*X(2)/10.0_DP)*COS(2.0_DP*PI*X(1)/10.0_DP))
                          ENDIF
                        ELSE IF(equationsAnalytic%analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_4) THEN
!                           do nothing!
                        ELSE IF(equationsAnalytic%analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_5) THEN
!                           do nothing!
                        ELSE IF(equationsAnalytic%analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_4) THEN
!                           do nothing!
                        ELSE IF(equationsAnalytic%analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_5) THEN
!                           do nothing!
                        ELSE IF(equationsAnalytic%analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_1) THEN
                          IF(rowComponentIdx==1) THEN
                            !Calculate SUM
                            SUM=0.0_DP
                          ELSE IF(rowComponentIdx==2) THEN
                            !Calculate SUM
                            SUM=rowsPhi*(-4.0_DP*muParam/100.0_DP)
                          ELSE IF(rowComponentIdx==3) THEN
                            !Calculate SUM
                            SUM=rowsPhi*(-4.0_DP*muParam/100.0_DP)
                          ENDIF
                        ELSE IF(equationsAnalytic%analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_2) THEN
                          IF(rowComponentIdx==1) THEN
                            !Calculate SUM
                            SUM=0.0_DP
                          ELSE IF(rowComponentIdx==2) THEN
                            !Calculate SUM
                            SUM=rowsPhi*(-2.0_DP*muParam/100.0_DP*(2.0_DP*EXP((X(1)-X(2))/10.0_DP)+EXP((X(2)-X(3))/10.0_DP)))
                          ELSE IF(rowComponentIdx==3) THEN
                            !Calculate SUM
                            SUM=rowsPhi*(-2.0_DP*muParam/100.0_DP*(2.0_DP*EXP((X(3)-X(1))/10.0_DP)+EXP((X(2)-X(3))/10.0_DP)))
                          ENDIF
                        ELSE IF(equationsAnalytic%analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_3) THEN
                          IF(rowComponentIdx==1) THEN
                            !Calculate SUM
                            SUM=0.0_DP
                          ELSE IF(rowComponentIdx==2) THEN
                            !Calculate SUM
                            SUM=rowsPhi*(36*muParam*PI**2/100.0_DP*COS(2.0_DP*PI*X(2)/10.0_DP)*SIN(2.0_DP*PI*X(3)/10.0_DP)* &
                              & COS(2.0_DP*PI*X(1)/10.0_DP))
                          ELSE IF(rowComponentIdx==3) THEN
                            !Calculate SUM
                            SUM=0.0_DP
                          ENDIF
                        ENDIF
                        !Calculate RH VECTOR
                        RH_VECTOR(rowElelemntDOFIdx)=RH_VECTOR(rowElelemntDOFIdx)+SUM*JGW
                      ENDDO !rowElementParameterIdx
                    ENDDO !rowComponentIdx
                  ELSE
                    RH_VECTOR(rowElelemntDOFIdx)=0.0_DP
                  ENDIF
                ENDIF
              END IF
            ENDDO !gaussPointIdx
            !Assemble matrices calculated above
            rowElelemntDOFIdx_min=rowElelemntDOFIdx
            rowElelemntDOFIdx_max=columnElementDOFIdx
            columnElementDOFIdx_min=rowElelemntDOFIdx
            columnElementDOFIdx_max=columnElementDOFIdx
            IF(equationsSet%specification(3)==EQUATIONS_SET_STATIC_STOKES_SUBTYPE.OR.  &
              & equationsSet%specification(3)==EQUATIONS_SET_LAPLACE_STOKES_SUBTYPE.OR. &
              & equationsSet%specification(3)==EQUATIONS_SET_ALE_STOKES_SUBTYPE.OR. &
              & equationsSet%specification(3)==EQUATIONS_SET_PGM_STOKES_SUBTYPE.OR. &
              & equationsSet%specification(3)==EQUATIONS_SET_TRANSIENT_STOKES_SUBTYPE) THEN
              IF(updateStiffnessMatrix) THEN
                stiffnessMatrix%elementMatrix%matrix(1:rowElelemntDOFIdx_min,1:columnElementDOFIdx_min)=AL_MATRIX(1:rowElelemntDOFIdx_min,1:columnElementDOFIdx_min)+AG_MATRIX(1:rowElelemntDOFIdx_min, &
                  & 1:columnElementDOFIdx_min)+ALE_MATRIX(1:rowElelemntDOFIdx_min,1:columnElementDOFIdx_min)
                stiffnessMatrix%elementMatrix%matrix(1:rowElelemntDOFIdx_min,columnElementDOFIdx_min+1:columnElementDOFIdx_max)=BT_MATRIX(1:rowElelemntDOFIdx_min,columnElementDOFIdx_min+1:columnElementDOFIdx_max)
                DO rowElelemntDOFIdx=rowElelemntDOFIdx_min+1,rowElelemntDOFIdx_max
                  DO columnElementDOFIdx=1,columnElementDOFIdx_min
                    !Transpose pressure type entries for mass equation
                    stiffnessMatrix%elementMatrix%matrix(rowElelemntDOFIdx,columnElementDOFIdx)=stiffnessMatrix%elementMatrix%matrix(columnElementDOFIdx,rowElelemntDOFIdx)
                  END DO
                END DO
              END IF
            END IF
            IF(equationsSet%specification(3)==EQUATIONS_SET_TRANSIENT_STOKES_SUBTYPE.OR. &
              & equationsSet%specification(3)==EQUATIONS_SET_ALE_STOKES_SUBTYPE.OR. &
              & equationsSet%specification(3)==EQUATIONS_SET_PGM_STOKES_SUBTYPE) THEN
              IF(updateDampingMatrix) THEN
                dampingMatrix%elementMatrix%matrix(1:rowElelemntDOFIdx_min,1:columnElementDOFIdx_min)=MT_MATRIX(1:rowElelemntDOFIdx_min,1:columnElementDOFIdx_min)
              END IF
            END IF
          !Assemble RHS vector
          IF(rhsVector%firstAssembly) THEN
            IF(updateRHSVector) THEN
              rhsVector%elementVector%vector(1:rowElelemntDOFIdx_max)=RH_VECTOR(1:rowElelemntDOFIdx_max)
            ENDIF
          ENDIF
          !Scale factor adjustment
            IF(dependentField%SCALINGS%scalingType/=FIELD_NO_SCALING) THEN
              CALL Field_InterpolationParametersScaleFactorsElementGet(elementNumber,equations%interpolation% &
                & dependentInterpParameters(variableType)%ptr,err,error,*999)
              rowElelemntDOFIdx=0
              DO rowComponentIdx=1,fieldVariable%numberOfComponents
                !Loop over element rows
                meshComponent1=fieldVariable%COMPONENTS(rowComponentIdx)%meshComponentNumber
                dependentBasis1=>dependentField%DECOMPOSITION%DOMAIN(meshComponent1)%ptr% &
                  & TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS
                DO rowElementParameterIdx=1,dependentBasis1%numberOfElementParameters
                  rowElelemntDOFIdx=rowElelemntDOFIdx+1
                  columnElementDOFIdx=0
                   IF(updateStiffnessMatrix.OR.updateDampingMatrix) THEN
                    !Loop over element columns
                    DO columnComponentIdx=1,fieldVariable%numberOfComponents
                      meshComponent2=fieldVariable%COMPONENTS(columnComponentIdx)%meshComponentNumber
                      dependentBasis2=>dependentField%DECOMPOSITION%DOMAIN(meshComponent2)%ptr% &
                        & TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS
                      DO columnElementParameterIdx=1,dependentBasis2%numberOfElementParameters
                        columnElementDOFIdx=columnElementDOFIdx+1
                        IF(updateStiffnessMatrix)THEN
                          stiffnessMatrix%elementMatrix%matrix(rowElelemntDOFIdx,columnElementDOFIdx)=stiffnessMatrix%elementMatrix%matrix(rowElelemntDOFIdx,columnElementDOFIdx)* &
                            & equations%interpolation%dependentInterpParameters(variableType)%ptr%scaleFactors(rowElementParameterIdx,rowComponentIdx)* &
                            & equations%interpolation%dependentInterpParameters(variableType)%ptr%scaleFactors(columnElementParameterIdx,columnComponentIdx)
                        END IF
                        IF(updateDampingMatrix)THEN
                          dampingMatrix%elementMatrix%matrix(rowElelemntDOFIdx,columnElementDOFIdx)=dampingMatrix%elementMatrix%matrix(rowElelemntDOFIdx,columnElementDOFIdx)* &
                            & equations%interpolation%dependentInterpParameters(variableType)%ptr%scaleFactors(rowElementParameterIdx,rowComponentIdx)* &
                            & equations%interpolation%dependentInterpParameters(variableType)%ptr%scaleFactors(columnElementParameterIdx,columnComponentIdx)
                        END IF
                      ENDDO !columnElementParameterIdx
                    ENDDO !columnComponentIdx
                  ENDIF
                  IF(updateRHSVector) rhsVector%elementVector%vector(rowElelemntDOFIdx)=rhsVector%elementVector%vector(rowElelemntDOFIdx)* &
                    & equations%interpolation%dependentInterpParameters(variableType)%ptr%scaleFactors(rowElementParameterIdx,rowComponentIdx)
                ENDDO !rowElementParameterIdx
              ENDDO !rowComponentIdx
            ENDIF
          CASE DEFAULT
            localError="Equations set subtype "//TRIM(NumberToVString(equationsSet%specification(3),"*",err,error))// &
              & " is not valid for a Stokes fluid type of a fluid mechanics equations set class."
            CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        CALL FlagError("Equations set equations is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

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
    TYPE(ControlLoopType), POINTER :: controlLoop 
    TYPE(SolverType), POINTER :: solver2 !<A pointer to the solver
    TYPE(VARYING_STRING) :: localError

    ENTERS("Stokes_PostSolve",err,error,*999)
    NULLIFY(solver2)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)

    IF(ASSOCIATED(controlLoop)) THEN
      IF(ASSOCIATED(solver)) THEN
        IF(ASSOCIATED(controlLoop%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(control_loop%problem%specification)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(control_loop%problem%specification,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a Stokes problem.",err,error,*999)
          END IF
          SELECT CASE(controlLoop%PROBLEM%specification(3))
            CASE(PROBLEM_STATIC_STOKES_SUBTYPE,PROBLEM_LAPLACE_STOKES_SUBTYPE)
              CALL STOKES_POST_SOLVE_OUTPUT_DATA(controlLoop,solver,err,error,*999)
            CASE(PROBLEM_PGM_STOKES_SUBTYPE)
              CALL STOKES_POST_SOLVE_OUTPUT_DATA(controlLoop,solver,err,error,*999)
            CASE(PROBLEM_TRANSIENT_STOKES_SUBTYPE)
              CALL STOKES_POST_SOLVE_OUTPUT_DATA(controlLoop,solver,err,error,*999)
            CASE(PROBLEM_ALE_STOKES_SUBTYPE)
              !Post solve for the linear solver
              IF(solver%solveType==SOLVER_LINEAR_TYPE) THEN
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Mesh movement post solve... ",err,error,*999)
                CALL Solvers_SolverGet(solver%solvers,2,solver2,err,error,*999)
                IF(ASSOCIATED(solver2%dynamicSolver)) THEN
                  solver2%dynamicSolver%ale=.TRUE.
                ELSE
                  CALL FlagError("Dynamic solver is not associated for ALE problem.",err,error,*999)
                END IF
              !Post solve for the linear solver
              ELSE IF(solver%solveType==SOLVER_DYNAMIC_TYPE) THEN
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"ALE Stokes post solve... ",err,error,*999)
                CALL STOKES_POST_SOLVE_OUTPUT_DATA(controlLoop,solver,err,error,*999)
              END IF
            CASE DEFAULT
              localError="Problem subtype "//TRIM(NumberToVString(controlLoop%PROBLEM%specification(3),"*",err,error))// &
                & " is not valid for a Stokes fluid type of a fluid mechanics problem class."
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
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverType), POINTER :: solver2 
    TYPE(VARYING_STRING) :: localError

    ENTERS("Stokes_PreSolve",err,error,*999)
    NULLIFY(solver2)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    IF(.NOT.ALLOCATED(problem%specification)) &
      & CALL FlagError("Problem specification is not allocated.",err,error,*999)
    IF(SIZE(problem%specification,1)<3) &
      & CALL FlagError("Problem specification must have three entries for a Stokes problem.",err,error,*999)
    
    SELECT CASE(problem%specification(3))
    CASE(PROBLEM_STATIC_STOKES_SUBTYPE,PROBLEM_LAPLACE_STOKES_SUBTYPE)
      ! do nothing ???
    CASE(PROBLEM_TRANSIENT_STOKES_SUBTYPE)
      ! do nothing ???
      CALL Stokes_PreSolveUpdateBoundaryConditions(solver,err,error,*999)
    CASE(PROBLEM_PGM_STOKES_SUBTYPE)
      ! do nothing ???
      !First update mesh and calculates boundary velocity values
      CALL Stokes_PreSolveALEUpdateMesh(solver,err,error,*999)
      !Then apply both normal and moving mesh boundary conditions
      CALL Stokes_PreSolveUpdateBoundaryConditions(solver,err,error,*999)
    CASE(PROBLEM_ALE_STOKES_SUBTYPE)
      !Pre solve for the linear solver
      IF(solver%solveType==SOLVER_LINEAR_TYPE) THEN
        CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Mesh movement pre solve... ",err,error,*999)
        !Update boundary conditions for mesh-movement
        CALL Stokes_PreSolveUpdateBoundaryConditions(solver,err,error,*999)
        CALL Solvers_SolverGet(solver%solvers,2,solver2,err,error,*999)
        IF(ASSOCIATED(solver2%dynamicSolver)) THEN
          !\todo: Avoid ALE flag in future
          solver2%dynamicSolver%ale=.FALSE.
        ELSE
          CALL FlagError("Dynamic solver is not associated for ALE problem.",err,error,*999)
        END IF
        !Update material properties for Laplace mesh movement
        CALL Stokes_PreSolveALEUpdateParameters(solver,err,error,*999)
        !Pre solve for the linear solver
      ELSE IF(solver%solveType==SOLVER_DYNAMIC_TYPE) THEN
        CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"ALE Stokes pre solve... ",err,error,*999)
        IF(solver%dynamicSolver%ale) THEN
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
      localError="Problem subtype "//TRIM(NumberToVString(controlLoop%PROBLEM%specification(3),"*",err,error))// &
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
    TYPE(SolverEquationsType), POINTER :: solverEquations  !<A pointer to the solver equations
    TYPE(SolverMappingType), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set
    TYPE(EquationsType), POINTER :: EQUATIONS
    TYPE(VARYING_STRING) :: localError
    TYPE(BoundaryConditionVariableType), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(BoundaryConditionsType), POINTER :: BOUNDARY_CONDITIONS
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the control loop to solve.
    TYPE(FieldType), POINTER :: dependentField,geometricField,materialsField
    TYPE(FieldVariableType), POINTER :: fieldVariable,GEOMETRIC_VARIABLE
    TYPE(DomainType), POINTER :: DOMAIN
    TYPE(DomainNodesType), POINTER :: DOMAIN_NODES
    TYPE(FieldInterpolatedPointPtrType), POINTER :: INTERPOLATED_POINT(:)
    TYPE(FieldInterpolationParametersPtrType), POINTER :: INTERPOLATION_PARAMETERS(:)
    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT,DISPLACEMENT_VALUE,VALUE,XI_COORDINATES(3)
    REAL(DP) :: T_COORDINATES(20,3)
    INTEGER(INTG) :: numberOfDimensions,BOUNDARY_CONDITION_CHECK_VARIABLE,GLOBAL_DERIV_INDEX,node_idx,variable_type
    INTEGER(INTG) :: variable_idx,local_ny,ANALYTIC_FUNCTION_TYPE,component_idx,deriv_idx,dim_idx
    INTEGER(INTG) :: element_idx,en_idx,I,J,K,number_of_nodes_xic(3)
    REAL(DP) :: X(3),muParam,rhoParam
    REAL(DP), POINTER :: MESH_VELOCITY_VALUES(:), GEOMETRIC_PARAMETERS(:)
    REAL(DP), POINTER :: BOUNDARY_VALUES(:)


    ENTERS("Stokes_PreSolveUpdateBoundaryConditions",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    CALL CONTROL_LOOP_CURRENT_TIMES_GET(controlLoop,CURRENT_TIME,TIME_INCREMENT,err,error,*999)
    IF(ASSOCIATED(controlLoop%PROBLEM)) THEN
      IF(.NOT.ALLOCATED(control_loop%problem%specification)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(control_loop%problem%specification,1)<3) THEN
        CALL FlagError("Problem specification must have three entries for a Stokes problem.",err,error,*999)
      END IF
      SELECT CASE(controlLoop%PROBLEM%specification(3))
      CASE(PROBLEM_STATIC_STOKES_SUBTYPE,PROBLEM_LAPLACE_STOKES_SUBTYPE)
        ! do nothing ???
      CASE(PROBLEM_TRANSIENT_STOKES_SUBTYPE)
        solverEquations=>SOLVER%solverEquations
        IF(ASSOCIATED(solverEquations)) THEN
          SOLVER_MAPPING=>SOLVER_equations%solverMapping
          EQUATIONS=>SOLVER_MAPPING%equationsSetToSolverMatricesMap(1)%EQUATIONS
          IF(ASSOCIATED(EQUATIONS)) THEN
            equationsSet=>equations%equationsSet
            IF(ASSOCIATED(equationsSet%ANALYTIC)) THEN
              IF(equationsAnalytic%analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_4 .OR. &
                & equationsAnalytic%analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_1 .OR. &
                & equationsAnalytic%analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_4 .OR. &
                & equationsAnalytic%analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_5 .OR. &
                & equationsAnalytic%analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_5) THEN
                IF(ASSOCIATED(equationsSet)) THEN
                  IF(ASSOCIATED(equationsSet%ANALYTIC)) THEN
                    dependentField=>equationsSet%dependent%dependentField
                    IF(ASSOCIATED(dependentField)) THEN
                      geometricField=>equationsSet%GEOMETRY%geometricField
                      IF(ASSOCIATED(geometricField)) THEN
                        CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,&
                          & numberOfDimensions,err,error,*999)
                        NULLIFY(INTERPOLATION_PARAMETERS)
                        NULLIFY(INTERPOLATED_POINT)
                        CALL Field_InterpolationParametersInitialise(geometricField,INTERPOLATION_PARAMETERS,err,error, &
                          & *999)
                        CALL Field_InterpolatedPointsInitialise(INTERPOLATION_PARAMETERS,INTERPOLATED_POINT,err,error,*999)
                        NULLIFY(GEOMETRIC_VARIABLE)
                        CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,GEOMETRIC_VARIABLE,err,error,*999)
                        NULLIFY(GEOMETRIC_PARAMETERS)
                        CALL Field_ParameterSetDataGet(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,&
                          & GEOMETRIC_PARAMETERS,err,error,*999)
                        DO variable_idx=1,dependentField%numberOfVariables
                          variable_type=dependentField%VARIABLES(variable_idx)%variableType
                          fieldVariable=>dependentField%variableTypeMap(variable_type)%ptr
                          IF(ASSOCIATED(fieldVariable)) THEN
                            DO component_idx=1,fieldVariable%numberOfComponents
                              IF(fieldVariable%COMPONENTS(component_idx)%interpolationType== &
                                & FIELD_NODE_BASED_INTERPOLATION) THEN
                                DOMAIN=>fieldVariable%COMPONENTS(component_idx)%DOMAIN
                                IF(ASSOCIATED(DOMAIN)) THEN
                                  IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                                    DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                                    IF(ASSOCIATED(DOMAIN_NODES)) THEN
                                      !Should be replaced by boundary node flag
                                      DO node_idx=1,DOMAIN_NODES%numberOfNodes
                                        element_idx=DOMAIN%topology%nodes%nodes(node_idx)%surroundingElements(1)
                                        CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,element_idx, &
                                          & INTERPOLATION_PARAMETERS(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                                        en_idx=0
                                        XI_COORDINATES=0.0_DP
                                        number_of_nodes_xic(1)=DOMAIN%topology%elements%elements(element_idx)% &
                                          & basis%numberOfNodesXiC(1)
                                        number_of_nodes_xic(2)=DOMAIN%topology%elements%elements(element_idx)% &
                                          & basis%numberOfNodesXiC(2)
                                        IF(numberOfDimensions==3) THEN
                                          number_of_nodes_xic(3)=DOMAIN%topology%elements%elements(element_idx)%basis% &
                                            & numberOfNodesXiC(3)
                                        ELSE
                                          number_of_nodes_xic(3)=1
                                        ENDIF
                                        !\todo: change definitions as soon as adjacent elements / boundary elements calculation works for simplex
                                        IF(DOMAIN%topology%elements%maximumNumberOfElementParameters==4.OR. &
                                          & DOMAIN%topology%elements%maximumNumberOfElementParameters==9.OR. &
                                          & DOMAIN%topology%elements%maximumNumberOfElementParameters==16.OR. &
                                          & DOMAIN%topology%elements%maximumNumberOfElementParameters==8.OR. &
                                          & DOMAIN%topology%elements%maximumNumberOfElementParameters==27.OR. &
                                          & DOMAIN%topology%elements%maximumNumberOfElementParameters==64) THEN
                                          DO K=1,number_of_nodes_xic(3)
                                            DO J=1,number_of_nodes_xic(2)
                                              DO I=1,number_of_nodes_xic(1)
                                                en_idx=en_idx+1
                                                IF(DOMAIN%topology%elements%elements(element_idx)% &
                                                  & elementNodes(en_idx)==node_idx) EXIT
                                                XI_COORDINATES(1)=XI_COORDINATES(1)+(1.0_DP/(number_of_nodes_xic(1)-1))
                                              ENDDO
                                              IF(DOMAIN%topology%elements%elements(element_idx)% &
                                                & elementNodes(en_idx)==node_idx) EXIT
                                              XI_COORDINATES(1)=0.0_DP
                                              XI_COORDINATES(2)=XI_COORDINATES(2)+(1.0_DP/(number_of_nodes_xic(2)-1))
                                            ENDDO
                                            IF(DOMAIN%topology%elements%elements(element_idx)% &
                                              & elementNodes(en_idx)==node_idx) EXIT
                                            XI_COORDINATES(1)=0.0_DP
                                            XI_COORDINATES(2)=0.0_DP
                                            IF(number_of_nodes_xic(3)/=1) THEN
                                              XI_COORDINATES(3)=XI_COORDINATES(3)+(1.0_DP/(number_of_nodes_xic(3)-1))
                                            ENDIF
                                          ENDDO
                                          CALL Field_InterpolateXi(NO_PART_DERIV,XI_COORDINATES, &
                                            & INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                                        ELSE
                                          !\todo: Use boundary flag
                                          IF(DOMAIN%topology%elements%maximumNumberOfElementParameters==3) THEN
                                            T_COORDINATES(1,1:2)=[0.0_DP,1.0_DP]
                                            T_COORDINATES(2,1:2)=[1.0_DP,0.0_DP]
                                            T_COORDINATES(3,1:2)=[1.0_DP,1.0_DP]
                                          ELSE IF(DOMAIN%topology%elements%maximumNumberOfElementParameters==6) THEN
                                            T_COORDINATES(1,1:2)=[0.0_DP,1.0_DP]
                                            T_COORDINATES(2,1:2)=[1.0_DP,0.0_DP]
                                            T_COORDINATES(3,1:2)=[1.0_DP,1.0_DP]
                                            T_COORDINATES(4,1:2)=[0.5_DP,0.5_DP]
                                            T_COORDINATES(5,1:2)=[1.0_DP,0.5_DP]
                                            T_COORDINATES(6,1:2)=[0.5_DP,1.0_DP]
                                          ELSE IF(DOMAIN%topology%elements%maximumNumberOfElementParameters==10.AND. &
                                            & numberOfDimensions==2) THEN
                                            T_COORDINATES(1,1:2)=[0.0_DP,1.0_DP]
                                            T_COORDINATES(2,1:2)=[1.0_DP,0.0_DP]
                                            T_COORDINATES(3,1:2)=[1.0_DP,1.0_DP]
                                            T_COORDINATES(4,1:2)=[1.0_DP/3.0_DP,2.0_DP/3.0_DP]
                                            T_COORDINATES(5,1:2)=[2.0_DP/3.0_DP,1.0_DP/3.0_DP]
                                            T_COORDINATES(6,1:2)=[1.0_DP,1.0_DP/3.0_DP]
                                            T_COORDINATES(7,1:2)=[1.0_DP,2.0_DP/3.0_DP]
                                            T_COORDINATES(8,1:2)=[2.0_DP/3.0_DP,1.0_DP]
                                            T_COORDINATES(9,1:2)=[1.0_DP/3.0_DP,1.0_DP]
                                            T_COORDINATES(10,1:2)=[2.0_DP/3.0_DP,2.0_DP/3.0_DP]
                                          ELSE IF(DOMAIN%topology%elements%maximumNumberOfElementParameters==4) THEN
                                            T_COORDINATES(1,1:3)=[0.0_DP,1.0_DP,1.0_DP]
                                            T_COORDINATES(2,1:3)=[1.0_DP,0.0_DP,1.0_DP]
                                            T_COORDINATES(3,1:3)=[1.0_DP,1.0_DP,0.0_DP]
                                            T_COORDINATES(4,1:3)=[1.0_DP,1.0_DP,1.0_DP]
                                          ELSE IF(DOMAIN%topology%elements%maximumNumberOfElementParameters==10.AND. &
                                            & numberOfDimensions==3) THEN
                                            T_COORDINATES(1,1:3)=[0.0_DP,1.0_DP,1.0_DP]
                                            T_COORDINATES(2,1:3)=[1.0_DP,0.0_DP,1.0_DP]
                                            T_COORDINATES(3,1:3)=[1.0_DP,1.0_DP,0.0_DP]
                                            T_COORDINATES(4,1:3)=[1.0_DP,1.0_DP,1.0_DP]
                                            T_COORDINATES(5,1:3)=[0.5_DP,0.5_DP,1.0_DP]
                                            T_COORDINATES(6,1:3)=[0.5_DP,1.0_DP,0.5_DP]
                                            T_COORDINATES(7,1:3)=[0.5_DP,1.0_DP,1.0_DP]
                                            T_COORDINATES(8,1:3)=[1.0_DP,0.5_DP,0.5_DP]
                                            T_COORDINATES(9,1:3)=[1.0_DP,1.0_DP,0.5_DP]
                                            T_COORDINATES(10,1:3)=[1.0_DP,0.5_DP,1.0_DP]
                                          ELSE IF(DOMAIN%topology%elements%maximumNumberOfElementParameters==20) THEN
                                            T_COORDINATES(1,1:3)=[0.0_DP,1.0_DP,1.0_DP]
                                            T_COORDINATES(2,1:3)=[1.0_DP,0.0_DP,1.0_DP]
                                            T_COORDINATES(3,1:3)=[1.0_DP,1.0_DP,0.0_DP]
                                            T_COORDINATES(4,1:3)=[1.0_DP,1.0_DP,1.0_DP]
                                            T_COORDINATES(5,1:3)=[1.0_DP/3.0_DP,2.0_DP/3.0_DP,1.0_DP]
                                            T_COORDINATES(6,1:3)=[2.0_DP/3.0_DP,1.0_DP/3.0_DP,1.0_DP]
                                            T_COORDINATES(7,1:3)=[1.0_DP/3.0_DP,1.0_DP,2.0_DP/3.0_DP]
                                            T_COORDINATES(8,1:3)=[2.0_DP/3.0_DP,1.0_DP,1.0_DP/3.0_DP]
                                            T_COORDINATES(9,1:3)=[1.0_DP/3.0_DP,1.0_DP,1.0_DP]
                                            T_COORDINATES(10,1:3)=[2.0_DP/3.0_DP,1.0_DP,1.0_DP]
                                            T_COORDINATES(11,1:3)=[1.0_DP,1.0_DP/3.0_DP,2.0_DP/3.0_DP]
                                            T_COORDINATES(12,1:3)=[1.0_DP,2.0_DP/3.0_DP,1.0_DP/3.0_DP]
                                            T_COORDINATES(13,1:3)=[1.0_DP,1.0_DP,1.0_DP/3.0_DP]
                                            T_COORDINATES(14,1:3)=[1.0_DP,1.0_DP,2.0_DP/3.0_DP]
                                            T_COORDINATES(15,1:3)=[1.0_DP,1.0_DP/3.0_DP,1.0_DP]
                                            T_COORDINATES(16,1:3)=[1.0_DP,2.0_DP/3.0_DP,1.0_DP]
                                            T_COORDINATES(17,1:3)=[2.0_DP/3.0_DP,2.0_DP/3.0_DP,2.0_DP/3.0_DP]
                                            T_COORDINATES(18,1:3)=[2.0_DP/3.0_DP,2.0_DP/3.0_DP,1.0_DP]
                                            T_COORDINATES(19,1:3)=[2.0_DP/3.0_DP,1.0_DP,2.0_DP/3.0_DP]
                                            T_COORDINATES(20,1:3)=[1.0_DP,2.0_DP/3.0_DP,2.0_DP/3.0_DP]
                                          ENDIF
                                          DO K=1,DOMAIN%topology%elements%maximumNumberOfElementParameters
                                            IF(DOMAIN%topology%elements%elements(element_idx)%elementNodes(K)==node_idx) EXIT
                                          ENDDO
                                          IF(numberOfDimensions==2) THEN
                                            CALL Field_InterpolateXi(NO_PART_DERIV,T_COORDINATES(K,1:2), &
                                              & INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                                          ELSE IF(numberOfDimensions==3) THEN
                                            CALL Field_InterpolateXi(NO_PART_DERIV,T_COORDINATES(K,1:3), &
                                              & INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                                          ENDIF
                                        ENDIF
                                        X=0.0_DP
                                        DO dim_idx=1,numberOfDimensions
                                          X(dim_idx)=INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(dim_idx,1)
                                        ENDDO !dim_idx
                                        !Loop over the derivatives
                                        NULLIFY(BOUNDARY_CONDITIONS_VARIABLE)
                                        CALL BoundaryConditions_VariableGet(SOLVER_equations%boundaryConditions, &
                                          & dependentField%variableTypeMap(FIELD_U_VARIABLE_TYPE)%ptr, &
                                          & BOUNDARY_CONDITIONS_VARIABLE,err,error,*999)
                                        DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%numberOfDerivatives
                                          ANALYTIC_FUNCTION_TYPE=equationsAnalytic%analyticFunctionType
                                          GLOBAL_DERIV_INDEX=DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)% &
                                            & globalDerivativeIndex
                                          materialsField=>equationsSet%MATERIALS%materialsField
                                          !Define muParam, density=1
                                          muParam=materialsField%variables(1)%parameterSets%parameterSets(1)%ptr% &
                                            & parameters%cmiss%dataDP(1)
                                          !Define rhoParam, density=2
                                          IF(ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_4.OR. &
                                            & ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_5.OR. &
                                            & ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_4.OR. &
                                            & ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_5) THEN
                                            rhoParam=materialsField%variables(1)%parameterSets%parameterSets(1)%ptr% &
                                              & parameters%cmiss%dataDP(2)
                                          ELSE
                                            rhoParam=0.0_DP
                                          ENDIF
                                          CALL STOKES_EQUATION_ANALYTIC_FUNCTIONS(VALUE,X,muParam,rhoParam,CURRENT_TIME, &
                                            & variable_type, &
                                            & GLOBAL_DERIV_INDEX,ANALYTIC_FUNCTION_TYPE,numberOfDimensions, &
                                            & fieldVariable%numberOfComponents,component_idx,err,error,*999)
                                          !Default to version 1 of each node derivative
                                          CALL FieldVariable_LocalNodeDOFGet(fieldVariable,1,deriv_idx,node_idx, &
                                            & component_idx,local_ny,err,error,*999)
                                          CALL Field_ParameterSetUpdateLocalDOF(dependentField,variable_type, &
                                            & FIELD_ANALYTIC_VALUES_SET_TYPE,local_ny,VALUE,err,error,*999)
                                          BOUNDARY_CONDITION_CHECK_VARIABLE=BOUNDARY_CONDITIONS_VARIABLE% &
                                            & conditionTypes(local_ny)
                                          IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_FIXED) THEN
                                            CALL Field_ParameterSetUpdateLocalDOF(dependentField, &
                                              & variable_type,FIELD_VALUES_SET_TYPE,local_ny, &
                                              & VALUE,err,error,*999)
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
                            CALL Field_ParameterSetUpdateStart(dependentField,variable_type, &
                              & FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
                            CALL Field_ParameterSetUpdateFinish(dependentField,variable_type, &
                              & FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
                            CALL Field_ParameterSetUpdateStart(dependentField,variable_type, &
                              & FIELD_VALUES_SET_TYPE,err,error,*999)
                            CALL Field_ParameterSetUpdateFinish(dependentField,variable_type, &
                              & FIELD_VALUES_SET_TYPE,err,error,*999)
                          ELSE
                            CALL FlagError("Field variable is not associated.",err,error,*999)
                          ENDIF
                        ENDDO !variable_idx
                        CALL Field_ParameterSetDataRestore(geometricField,FIELD_U_VARIABLE_TYPE,&
                          & FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS,err,error,*999)
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
              ENDIF
            ENDIF
          ELSE
            CALL FlagError("Equations are not associated.",err,error,*999)
          END IF
        ELSE
          CALL FlagError("Solver equations are not associated.",err,error,*999)
        END IF
        CALL Field_ParameterSetUpdateStart(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
          & FIELD_VALUES_SET_TYPE,err,error,*999)
        CALL Field_ParameterSetUpdateFinish(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
          & FIELD_VALUES_SET_TYPE,err,error,*999)
      CASE(PROBLEM_PGM_STOKES_SUBTYPE)
        !Pre solve for the dynamic solver
        IF(SOLVER%solveType==SOLVER_DYNAMIC_TYPE) THEN
          CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Mesh movement change boundary conditions... ",err,error,*999)
          solverEquations=>SOLVER%solverEquations
          IF(ASSOCIATED(solverEquations)) THEN
            SOLVER_MAPPING=>SOLVER_equations%solverMapping
            EQUATIONS=>SOLVER_MAPPING%equationsSetToSolverMatricesMap(1)%EQUATIONS
            IF(ASSOCIATED(EQUATIONS)) THEN
              equationsSet=>equations%equationsSet
              IF(ASSOCIATED(equationsSet)) THEN
                BOUNDARY_CONDITIONS=>SOLVER_equations%boundaryConditions
                IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
                  NULLIFY(BOUNDARY_CONDITIONS_VARIABLE)
                  CALL BoundaryConditions_VariableGet(BOUNDARY_CONDITIONS,equationsSet%dependent%dependentField% &
                    & variableTypeMap(FIELD_U_VARIABLE_TYPE)%ptr,BOUNDARY_CONDITIONS_VARIABLE,err,error,*999)
                  CALL Field_NumberOfComponentsGet(equationsSet%GEOMETRY%geometricField,FIELD_U_VARIABLE_TYPE, &
                    & numberOfDimensions,err,error,*999)
                  NULLIFY(MESH_VELOCITY_VALUES)
                  CALL Field_ParameterSetDataGet(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_MESH_VELOCITY_SET_TYPE,MESH_VELOCITY_VALUES,err,error,*999)
                  NULLIFY(BOUNDARY_VALUES)
                  CALL Field_ParameterSetDataGet(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_BOUNDARY_SET_TYPE,BOUNDARY_VALUES,err,error,*999)
                  CALL FLUID_MECHANICS_IO_READ_BOUNDARY_CONDITIONS(SOLVER_LINEAR_TYPE,BOUNDARY_VALUES, &
                    & numberOfDimensions,BOUNDARY_CONDITION_FIXED_INLET,controlLoop%timeLoop%inputNumber, &
                    & controlLoop%timeLoop%iterationNumber,CURRENT_TIME,1.0_DP,err,error,*999)
                  !                           DO equations_row_number=1,vectorEquations%vectorMapping%totalNumberOfRows
                  ! xxxxxxxxxxxxxxxxxxxxxx
                  DO variable_idx=1,equationsSet%dependent%dependentField%numberOfVariables
                    variable_type=equationsSet%dependent%dependentField%VARIABLES(variable_idx)%variableType
                    fieldVariable=>equationsSet%dependent%dependentField%variableTypeMap(variable_type)%ptr
                    IF(ASSOCIATED(fieldVariable)) THEN
                      DO component_idx=1,fieldVariable%numberOfComponents
                        DOMAIN=>fieldVariable%COMPONENTS(component_idx)%DOMAIN
                        IF(ASSOCIATED(DOMAIN)) THEN
                          IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                            DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                            IF(ASSOCIATED(DOMAIN_NODES)) THEN
                              !Loop over the local nodes excluding the ghosts.
                              DO node_idx=1,DOMAIN_NODES%numberOfNodes
                                DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%numberOfDerivatives
                                  !Default to version 1 of each node derivative
                                  CALL FieldVariable_LocalNodeDOFGet(fieldVariable,1,deriv_idx,node_idx, &
                                    & component_idx,local_ny,err,error,*999)
                                  ! xxxxxxxxxxxxxxxxxxxxxxxxx
                                  DISPLACEMENT_VALUE=0.0_DP
                                  BOUNDARY_CONDITION_CHECK_VARIABLE=BOUNDARY_CONDITIONS_VARIABLE% &
                                    & conditionTypes(local_ny)
                                  IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_MOVED_WALL) THEN
                                    CALL Field_ParameterSetUpdateLocalDOF(equationsSet%dependent%dependentField, &
                                      & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny, &
                                      & MESH_VELOCITY_VALUES(local_ny),err,error,*999)
                                  ELSE IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_FIXED_INLET) THEN
                                    CALL Field_ParameterSetUpdateLocalDOF(equationsSet%dependent%dependentField, &
                                      & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny, &
                                      & BOUNDARY_VALUES(local_ny),err,error,*999)
                                  END IF
                                ENDDO !deriv_idx
                              ENDDO !node_idx
                            ENDIF
                          ENDIF
                        ENDIF
                      ENDDO !component_idx
                    ENDIF
                  ENDDO !variable_idx
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
          CALL Field_ParameterSetUpdateStart(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,err,error,*999)
          CALL Field_ParameterSetUpdateFinish(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,err,error,*999)
        END IF
      CASE(PROBLEM_ALE_STOKES_SUBTYPE)
        !Pre solve for the linear solver
        IF(SOLVER%solveType==SOLVER_LINEAR_TYPE) THEN
          CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Mesh movement change boundary conditions... ",err,error,*999)
          solverEquations=>SOLVER%solverEquations
          IF(ASSOCIATED(solverEquations)) THEN
            SOLVER_MAPPING=>SOLVER_equations%solverMapping
            EQUATIONS=>SOLVER_MAPPING%equationsSetToSolverMatricesMap(1)%EQUATIONS
            IF(ASSOCIATED(EQUATIONS)) THEN
              equationsSet=>equations%equationsSet
              IF(ASSOCIATED(equationsSet)) THEN
                BOUNDARY_CONDITIONS=>SOLVER_equations%boundaryConditions
                IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
                  NULLIFY(BOUNDARY_CONDITIONS_VARIABLE)
                  CALL BoundaryConditions_VariableGet(BOUNDARY_CONDITIONS,equationsSet%dependent%dependentField% &
                    & variableTypeMap(FIELD_U_VARIABLE_TYPE)%ptr,BOUNDARY_CONDITIONS_VARIABLE,err,error,*999)
                  CALL Field_NumberOfComponentsGet(equationsSet%GEOMETRY%geometricField,FIELD_U_VARIABLE_TYPE, &
                    & numberOfDimensions,err,error,*999)
                  NULLIFY(BOUNDARY_VALUES)
                  CALL Field_ParameterSetDataGet(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_BOUNDARY_SET_TYPE,BOUNDARY_VALUES,err,error,*999)
                  CALL FLUID_MECHANICS_IO_READ_BOUNDARY_CONDITIONS(SOLVER_LINEAR_TYPE,BOUNDARY_VALUES, &
                    & numberOfDimensions,BOUNDARY_CONDITION_MOVED_WALL,controlLoop%timeLoop%inputNumber, &
                    & controlLoop%timeLoop%iterationNumber,CURRENT_TIME,1.0_DP,err,error,*999)
                  DO variable_idx=1,equationsSet%dependent%dependentField%numberOfVariables
                    variable_type=equationsSet%dependent%dependentField%VARIABLES(variable_idx)%variableType
                    fieldVariable=>equationsSet%dependent%dependentField%variableTypeMap(variable_type)%ptr
                    IF(ASSOCIATED(fieldVariable)) THEN
                      DO component_idx=1,fieldVariable%numberOfComponents
                        DOMAIN=>fieldVariable%COMPONENTS(component_idx)%DOMAIN
                        IF(ASSOCIATED(DOMAIN)) THEN
                          IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                            DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                            IF(ASSOCIATED(DOMAIN_NODES)) THEN
                              !Loop over the local nodes excluding the ghosts.
                              DO node_idx=1,DOMAIN_NODES%numberOfNodes
                                DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%numberOfDerivatives
                                  !Default to version 1 of each node derivative
                                  CALL FieldVariable_LocalNodeDOFGet(fieldVariable,1,deriv_idx,node_idx, &
                                    & component_idx,local_ny,err,error,*999)
                                  BOUNDARY_CONDITION_CHECK_VARIABLE=BOUNDARY_CONDITIONS_VARIABLE% &
                                    & conditionTypes(local_ny)
                                  IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_MOVED_WALL) THEN
                                    CALL Field_ParameterSetUpdateLocalDOF(equationsSet%dependent%dependentField, &
                                      & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny, &
                                      & BOUNDARY_VALUES(local_ny),err,error,*999)
                                  END IF
                                END DO !deriv_idx
                              ENDDO !node_idx
                            ENDIF
                          ENDIF
                        ENDIF
                      ENDDO !component_idx
                    ENDIF
                  ENDDO !variable_idx
                  CALL Field_ParameterSetDataRestore(equationsIndependent%independentField, &
                    & FIELD_U_VARIABLE_TYPE,FIELD_BOUNDARY_SET_TYPE,BOUNDARY_VALUES,err,error,*999)
                  !\todo: This part should be read in out of a file eventually
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
          CALL Field_ParameterSetUpdateStart(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,err,error,*999)
          CALL Field_ParameterSetUpdateFinish(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,err,error,*999)
          !Pre solve for the dynamic solver
        ELSE IF(SOLVER%solveType==SOLVER_DYNAMIC_TYPE) THEN
          CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Mesh movement change boundary conditions... ",err,error,*999)
          solverEquations=>SOLVER%solverEquations
          IF(ASSOCIATED(solverEquations)) THEN
            SOLVER_MAPPING=>SOLVER_equations%solverMapping
            EQUATIONS=>SOLVER_MAPPING%equationsSetToSolverMatricesMap(1)%EQUATIONS
            IF(ASSOCIATED(EQUATIONS)) THEN
              equationsSet=>equations%equationsSet
              IF(ASSOCIATED(equationsSet)) THEN
                BOUNDARY_CONDITIONS=>SOLVER_equations%boundaryConditions
                IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
                  NULLIFY(BOUNDARY_CONDITIONS_VARIABLE)
                  CALL BoundaryConditions_VariableGet(BOUNDARY_CONDITIONS,equationsSet%dependent%dependentField% &
                    & variableTypeMap(FIELD_U_VARIABLE_TYPE)%ptr,BOUNDARY_CONDITIONS_VARIABLE,err,error,*999)
                  CALL Field_NumberOfComponentsGet(equationsSet%GEOMETRY%geometricField,FIELD_U_VARIABLE_TYPE, &
                    & numberOfDimensions,err,error,*999)
                  NULLIFY(MESH_VELOCITY_VALUES)
                  CALL Field_ParameterSetDataGet(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_MESH_VELOCITY_SET_TYPE,MESH_VELOCITY_VALUES,err,error,*999)
                  NULLIFY(BOUNDARY_VALUES)
                  CALL Field_ParameterSetDataGet(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_BOUNDARY_SET_TYPE,BOUNDARY_VALUES,err,error,*999)
                  CALL FLUID_MECHANICS_IO_READ_BOUNDARY_CONDITIONS(SOLVER_LINEAR_TYPE,BOUNDARY_VALUES, &
                    & numberOfDimensions,BOUNDARY_CONDITION_FIXED_INLET,controlLoop%timeLoop%inputNumber, &
                    & controlLoop%timeLoop%iterationNumber,CURRENT_TIME,1.0_DP,err,error,*999)
                  DO variable_idx=1,equationsSet%dependent%dependentField%numberOfVariables
                    variable_type=equationsSet%dependent%dependentField%VARIABLES(variable_idx)%variableType
                    fieldVariable=>equationsSet%dependent%dependentField%variableTypeMap(variable_type)%ptr
                    IF(ASSOCIATED(fieldVariable)) THEN
                      DO component_idx=1,fieldVariable%numberOfComponents
                        DOMAIN=>fieldVariable%COMPONENTS(component_idx)%DOMAIN
                        IF(ASSOCIATED(DOMAIN)) THEN
                          IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                            DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                            IF(ASSOCIATED(DOMAIN_NODES)) THEN
                              !Loop over the local nodes excluding the ghosts.
                              DO node_idx=1,DOMAIN_NODES%numberOfNodes
                                DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%numberOfDerivatives
                                  !Default to version 1 of each node derivative
                                  CALL FieldVariable_LocalNodeDOFGet(fieldVariable,1,deriv_idx,node_idx, &
                                    & component_idx,local_ny,err,error,*999)
                                  DISPLACEMENT_VALUE=0.0_DP
                                  BOUNDARY_CONDITION_CHECK_VARIABLE=BOUNDARY_CONDITIONS_VARIABLE% &
                                    & conditionTypes(local_ny)
                                  IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_MOVED_WALL) THEN
                                    CALL Field_ParameterSetUpdateLocalDOF(equationsSet%dependent%dependentField, &
                                      & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny, &
                                      & MESH_VELOCITY_VALUES(local_ny),err,error,*999)
                                  ELSE IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_FIXED_INLET) THEN
                                    CALL Field_ParameterSetUpdateLocalDOF(equationsSet%dependent%dependentField, &
                                      & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny, &
                                      & BOUNDARY_VALUES(local_ny),err,error,*999)
                                  END IF
                                END DO !deriv_idx
                              ENDDO !node_idx
                            ENDIF
                          ENDIF
                        ENDIF
                      ENDDO !component_idx
                    ENDIF
                  ENDDO !variable_idx
                  CALL Field_ParameterSetDataRestore(equationsIndependent%independentField, &
                    & FIELD_U_VARIABLE_TYPE,FIELD_MESH_VELOCITY_SET_TYPE,MESH_VELOCITY_VALUES,err,error,*999)
                  CALL Field_ParameterSetDataRestore(equationsIndependent%independentField, &
                    & FIELD_U_VARIABLE_TYPE,FIELD_BOUNDARY_SET_TYPE,BOUNDARY_VALUES,err,error,*999)
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
          CALL Field_ParameterSetUpdateStart(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,err,error,*999)
          CALL Field_ParameterSetUpdateFinish(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,err,error,*999)
        END IF
        ! do nothing ???
      CASE DEFAULT
        localError="Problem subtype "//TRIM(NumberToVString(controlLoop%PROBLEM%specification(3),"*",err,error))// &
          & " is not valid for a Stokes equation fluid type of a fluid mechanics problem class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    ENDIF
      
    EXITS("Stokes_PreSolveUpdateBoundaryConditions")
    RETURN
999 ERRORSEXITS("Stokes_PreSolveUpdateBoundaryConditions",err,error)
    RETURN 1
    
  END SUBROUTINE Stokes_PreSolveUpdateBoundaryConditions

  !
  !================================================================================================================================
  !
  !>Update mesh velocity and move mesh for ALE Stokes problem
  SUBROUTINE Stokes_PreSolveALEUpdateMesh(SOLVER,err,error,*)

    !Argument variables
   TYPE(SolverType), POINTER :: SOLVER !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ControlLoopType), POINTER :: controlLoop 
    TYPE(SolverType), POINTER :: SOLVER_ALE_STOKES, SOLVER_LAPLACE !<A pointer to the solvers
    TYPE(FieldType), POINTER :: DEPENDENT_FIELD_LAPLACE, INDEPENDENT_FIELD_ALE_STOKES
    TYPE(SolverEquationsType), POINTER :: SOLVER_EQUATIONS_LAPLACE, SOLVER_EQUATIONS_ALE_STOKES  !<A pointer to the solver equations
    TYPE(SolverMappingType), POINTER :: SOLVER_MAPPING_LAPLACE, SOLVER_MAPPING_ALE_STOKES !<A pointer to the solver mapping
    TYPE(EquationsSetType), POINTER :: EQUATIONS_SET_LAPLACE, EQUATIONS_SET_ALE_STOKES !<A pointer to the equations set
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(VARYING_STRING) :: localError
    TYPE(FieldVariableType), POINTER :: fieldVariable
    TYPE(DomainType), POINTER :: DOMAIN
    TYPE(DomainNodesType), POINTER :: DOMAIN_NODES

    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT,ALPHA
    REAL(DP), POINTER :: MESH_DISPLACEMENT_VALUES(:)
    INTEGER(INTG) :: I,numberOfDimensions_LAPLACE,numberOfDimensions_ALE_STOKES,geometricMeshComponent
    INTEGER(INTG) :: INPUT_TYPE,INPUT_OPTION,component_idx,deriv_idx,local_ny,node_idx,variable_idx,variable_type

    ENTERS("Stokes_PreSolveALEUpdateMesh",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(SOLVER,controlLoop,err,error,*999)
    CALL CONTROL_LOOP_CURRENT_TIMES_GET(controlLoop,CURRENT_TIME,TIME_INCREMENT,err,error,*999)
    NULLIFY(SOLVER_LAPLACE)
    NULLIFY(SOLVER_ALE_STOKES)
    IF(ASSOCIATED(SOLVER)) THEN
      IF(ASSOCIATED(controlLoop%PROBLEM)) THEN
        IF(.NOT.ALLOCATED(control_loop%problem%specification)) THEN
          CALL FlagError("Problem specification is not allocated.",err,error,*999)
        ELSE IF(SIZE(control_loop%problem%specification,1)<3) THEN
          CALL FlagError("Problem specification must have three entries for a Stokes problem.",err,error,*999)
        END IF
        SELECT CASE(controlLoop%PROBLEM%specification(3))
        CASE(PROBLEM_STATIC_STOKES_SUBTYPE,PROBLEM_LAPLACE_STOKES_SUBTYPE)
          ! do nothing ???
        CASE(PROBLEM_TRANSIENT_STOKES_SUBTYPE)
          ! do nothing ???
        CASE(PROBLEM_PGM_STOKES_SUBTYPE)
          !Update mesh within the dynamic solver
          IF(SOLVER%solveType==SOLVER_DYNAMIC_TYPE) THEN
            !Get the independent field for the ALE Stokes problem
            CALL Solvers_SolverGet(SOLVER%solvers,1,SOLVER_ALE_STOKES,err,error,*999)
            SOLVER_EQUATIONS_ALE_STOKES=>SOLVER_ALE_STOKES%solverEquations
            IF(ASSOCIATED(SOLVER_EQUATIONS_ALE_STOKES)) THEN
              SOLVER_MAPPING_ALE_STOKES=>SOLVER_EQUATIONS_ALE_STOKES%solverMapping
              IF(ASSOCIATED(SOLVER_MAPPING_ALE_STOKES)) THEN
                EQUATIONS_SET_ALE_STOKES=>SOLVER_MAPPING_ALE_STOKES%equationsSets(1)%ptr
                IF(ASSOCIATED(EQUATIONS_SET_ALE_STOKES)) THEN
                  INDEPENDENT_FIELD_ALE_STOKES=>EQUATIONS_SET_ALE_STOKES%INDEPENDENT%independentField
                ELSE
                  CALL FlagError("ALE Stokes equations set is not associated.",err,error,*999)
                END IF
                !Get the data
                CALL Field_NumberOfComponentsGet(EQUATIONS_SET_ALE_STOKES%GEOMETRY%geometricField, &
                  & FIELD_U_VARIABLE_TYPE,numberOfDimensions_ALE_STOKES,err,error,*999)
                !\todo: Introduce flags set by the user (42/1 only for testings purpose)
                !Copy input to Stokes' independent field
                INPUT_TYPE=42
                INPUT_OPTION=1
                NULLIFY(MESH_DISPLACEMENT_VALUES)
                CALL Field_ParameterSetDataGet(EQUATIONS_SET_ALE_STOKES%INDEPENDENT%independentField, &
                  & FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE,MESH_DISPLACEMENT_VALUES,err,error,*999)
                CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,MESH_DISPLACEMENT_VALUES, &
                  & numberOfDimensions_ALE_STOKES,INPUT_TYPE,INPUT_OPTION,controlLoop%timeLoop%iterationNumber,1.0_DP, &
                  & err,error,*999)
                CALL Field_ParameterSetUpdateStart(EQUATIONS_SET_ALE_STOKES%INDEPENDENT%independentField, &
                  & FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE,err,error,*999)
                CALL Field_ParameterSetUpdateFinish(EQUATIONS_SET_ALE_STOKES%INDEPENDENT%independentField, &
                  & FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE,err,error,*999)
              ELSE
                CALL FlagError("ALE Stokes solver mapping is not associated.",err,error,*999)
              END IF
            ELSE
              CALL FlagError("ALE Stokes solver equations are not associated.",err,error,*999)
            END IF
            !Use calculated values to update mesh
            CALL Field_ComponentMeshComponentGet(EQUATIONS_SET_ALE_STOKES%GEOMETRY%geometricField, &
              & FIELD_U_VARIABLE_TYPE,1,geometricMeshComponent,err,error,*999)
            !                 CALL Field_ParameterSetDataGet(INDEPENDENT_FIELD_ALE_STOKES,FIELD_U_VARIABLE_TYPE, &
            !                   & FIELD_MESH_DISPLACEMENT_SET_TYPE,MESH_DISPLACEMENT_VALUES,err,error,*999)
            EQUATIONS=>SOLVER_MAPPING_ALE_STOKES%equationsSetToSolverMatricesMap(1)%EQUATIONS
            IF(ASSOCIATED(EQUATIONS)) THEN
              NULLIFY(vectorEquations)
              CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
              vectorMapping=>vectorEquations%vectorMapping
              IF(ASSOCIATED(vectorMapping)) THEN
                DO variable_idx=1,EQUATIONS_SET_ALE_STOKES%dependent%dependentField%numberOfVariables
                  variable_type=EQUATIONS_SET_ALE_STOKES%dependent%dependentField%VARIABLES(variable_idx)%variableType
                  fieldVariable=>EQUATIONS_SET_ALE_STOKES%GEOMETRY%geometricField%variableTypeMap(variable_type)%ptr
                  IF(ASSOCIATED(fieldVariable)) THEN
                    DO component_idx=1,fieldVariable%numberOfComponents
                      DOMAIN=>fieldVariable%COMPONENTS(component_idx)%DOMAIN
                      IF(ASSOCIATED(DOMAIN)) THEN
                        IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                          DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                          IF(ASSOCIATED(DOMAIN_NODES)) THEN
                            !Loop over the local nodes excluding the ghosts.
                            DO node_idx=1,DOMAIN_NODES%numberOfNodes
                              DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%numberOfDerivatives
                                !Default to version 1 of each node derivative
                                CALL FieldVariable_LocalNodeDOFGet(fieldVariable,1,deriv_idx,node_idx, &
                                  & component_idx,local_ny,err,error,*999)
                                CALL Field_ParameterSetAddLocalDOF(EQUATIONS_SET_ALE_STOKES%GEOMETRY%geometricField, &
                                  & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny, &
                                  & MESH_DISPLACEMENT_VALUES(local_ny),err,error,*999)
                              ENDDO !deriv_idx
                            ENDDO !node_idx
                          ENDIF
                        ENDIF
                      ENDIF
                    ENDDO !component_idx
                  ENDIF
                ENDDO !variable_idx
              ELSE
                CALL FlagError("Equations mapping is not associated.",err,error,*999)
              END IF
            ELSE
              CALL FlagError("Equations are not associated.",err,error,*999)
            END IF
            CALL Field_ParameterSetUpdateStart(EQUATIONS_SET_ALE_STOKES%GEOMETRY%geometricField, &
              & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
            CALL Field_ParameterSetUpdateFinish(EQUATIONS_SET_ALE_STOKES%GEOMETRY%geometricField, &
              & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
            !Now use displacement values to calculate velocity values
            TIME_INCREMENT=controlLoop%timeLoop%timeIncrement
            ALPHA=1.0_DP/TIME_INCREMENT
            CALL Field_ParameterSetsCopy(INDEPENDENT_FIELD_ALE_STOKES,FIELD_U_VARIABLE_TYPE, &
              & FIELD_MESH_DISPLACEMENT_SET_TYPE,FIELD_MESH_VELOCITY_SET_TYPE,ALPHA,err,error,*999)
          ELSE
            CALL FlagError("Mesh motion calculation not successful for ALE problem.",err,error,*999)
          END IF
        CASE(PROBLEM_ALE_STOKES_SUBTYPE)
          !Update mesh within the dynamic solver
          IF(SOLVER%solveType==SOLVER_DYNAMIC_TYPE) THEN
            IF(SOLVER%dynamicSolver%ale) THEN
              !Get the dependent field for the three component Laplace problem
              CALL Solvers_SolverGet(SOLVER%solvers,1,SOLVER_LAPLACE,err,error,*999)
              SOLVER_EQUATIONS_LAPLACE=>SOLVER_LAPLACE%solverEquations
              IF(ASSOCIATED(SOLVER_EQUATIONS_LAPLACE)) THEN
                SOLVER_MAPPING_LAPLACE=>SOLVER_EQUATIONS_LAPLACE%solverMapping
                IF(ASSOCIATED(SOLVER_MAPPING_LAPLACE)) THEN
                  EQUATIONS_SET_LAPLACE=>SOLVER_MAPPING_LAPLACE%equationsSets(1)%ptr
                  IF(ASSOCIATED(EQUATIONS_SET_LAPLACE)) THEN
                    DEPENDENT_FIELD_LAPLACE=>EQUATIONS_SET_LAPLACE%dependent%dependentField
                  ELSE
                    CALL FlagError("Laplace equations set is not associated.",err,error,*999)
                  END IF
                  CALL Field_NumberOfComponentsGet(EQUATIONS_SET_LAPLACE%GEOMETRY%geometricField,FIELD_U_VARIABLE_TYPE, &
                    & numberOfDimensions_LAPLACE,err,error,*999)
                ELSE
                  CALL FlagError("Laplace solver mapping is not associated.",err,error,*999)
                END IF
              ELSE
                CALL FlagError("Laplace solver equations are not associated.",err,error,*999)
              END IF
              !Get the independent field for the ALE Stokes problem
              CALL Solvers_SolverGet(SOLVER%solvers,2,SOLVER_ALE_STOKES,err,error,*999)
              SOLVER_EQUATIONS_ALE_STOKES=>SOLVER_ALE_STOKES%solverEquations
              IF(ASSOCIATED(SOLVER_EQUATIONS_ALE_STOKES)) THEN
                SOLVER_MAPPING_ALE_STOKES=>SOLVER_EQUATIONS_ALE_STOKES%solverMapping
                IF(ASSOCIATED(SOLVER_MAPPING_ALE_STOKES)) THEN
                  EQUATIONS_SET_ALE_STOKES=>SOLVER_MAPPING_ALE_STOKES%equationsSets(1)%ptr
                  IF(ASSOCIATED(EQUATIONS_SET_ALE_STOKES)) THEN
                    INDEPENDENT_FIELD_ALE_STOKES=>EQUATIONS_SET_ALE_STOKES%INDEPENDENT%independentField
                  ELSE
                    CALL FlagError("ALE Stokes equations set is not associated.",err,error,*999)
                  END IF
                  CALL Field_NumberOfComponentsGet(EQUATIONS_SET_ALE_STOKES%GEOMETRY%geometricField, &
                    & FIELD_U_VARIABLE_TYPE,numberOfDimensions_ALE_STOKES,err,error,*999)
                ELSE
                  CALL FlagError("ALE Stokes solver mapping is not associated.",err,error,*999)
                END IF
              ELSE
                CALL FlagError("ALE Stokes solver equations are not associated.",err,error,*999)
              END IF
              !Copy result from Laplace mesh movement to Stokes' independent field
              IF(numberOfDimensions_ALE_STOKES==numberOfDimensions_LAPLACE) THEN
                DO I=1,numberOfDimensions_ALE_STOKES
                  CALL Field_ParametersToFieldParametersCopy(DEPENDENT_FIELD_LAPLACE, &
                    & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,I,INDEPENDENT_FIELD_ALE_STOKES, &
                    & FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE,I,err,error,*999)
                END DO
              ELSE
                CALL FlagError("Dimension of Laplace and ALE Stokes equations set is not consistent.",err,error,*999)
              END IF
              !Use calculated values to update mesh
              CALL Field_ComponentMeshComponentGet(EQUATIONS_SET_ALE_STOKES%GEOMETRY%geometricField, &
                & FIELD_U_VARIABLE_TYPE,1,geometricMeshComponent,err,error,*999)
              NULLIFY(MESH_DISPLACEMENT_VALUES)
              CALL Field_ParameterSetDataGet(INDEPENDENT_FIELD_ALE_STOKES,FIELD_U_VARIABLE_TYPE, &
                & FIELD_MESH_DISPLACEMENT_SET_TYPE,MESH_DISPLACEMENT_VALUES,err,error,*999)
              EQUATIONS=>SOLVER_MAPPING_LAPLACE%equationsSetToSolverMatricesMap(1)%EQUATIONS
              IF(ASSOCIATED(EQUATIONS)) THEN
                NULLIFY(vectorEquations)
                CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
                vectorMapping=>vectorEquations%vectorMapping
                IF(ASSOCIATED(vectorMapping)) THEN
                  DO variable_idx=1,EQUATIONS_SET_ALE_STOKES%dependent%dependentField%numberOfVariables
                    variable_type=EQUATIONS_SET_ALE_STOKES%dependent%dependentField%VARIABLES(variable_idx)%variableType
                    fieldVariable=>EQUATIONS_SET_ALE_STOKES%GEOMETRY%geometricField%variableTypeMap(variable_type)%ptr
                    IF(ASSOCIATED(fieldVariable)) THEN
                      DO component_idx=1,fieldVariable%numberOfComponents
                        DOMAIN=>fieldVariable%COMPONENTS(component_idx)%DOMAIN
                        IF(ASSOCIATED(DOMAIN)) THEN
                          IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                            DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                            IF(ASSOCIATED(DOMAIN_NODES)) THEN
                              !Loop over the local nodes excluding the ghosts.
                              DO node_idx=1,DOMAIN_NODES%numberOfNodes
                                DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%numberOfDerivatives
                                  !Default to version 1 of each node derivative
                                  CALL FieldVariable_LocalNodeDOFGet(fieldVariable,1,deriv_idx,node_idx, &
                                    & component_idx,local_ny,err,error,*999)
                                  CALL Field_ParameterSetAddLocalDOF(EQUATIONS_SET_ALE_STOKES%GEOMETRY%geometricField, &
                                    & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny, &
                                    & MESH_DISPLACEMENT_VALUES(local_ny),err,error,*999)
                                ENDDO !deriv_idx
                              ENDDO !node_idx
                            ENDIF
                          ENDIF
                        ENDIF
                      ENDDO !component_idx
                    ENDIF
                  ENDDO !variable_idx
                ELSE
                  CALL FlagError("Equations mapping is not associated.",err,error,*999)
                ENDIF
                CALL Field_ParameterSetDataRestore(INDEPENDENT_FIELD_ALE_STOKES,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_MESH_DISPLACEMENT_SET_TYPE,MESH_DISPLACEMENT_VALUES,err,error,*999)
              ELSE
                CALL FlagError("Equations are not associated.",err,error,*999)
              END IF
              CALL Field_ParameterSetUpdateStart(EQUATIONS_SET_ALE_STOKES%GEOMETRY%geometricField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,err,error,*999)
              CALL Field_ParameterSetUpdateFinish(EQUATIONS_SET_ALE_STOKES%GEOMETRY%geometricField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,err,error,*999)
              !Now use displacement values to calculate velocity values
              TIME_INCREMENT=controlLoop%timeLoop%timeIncrement
              ALPHA=1.0_DP/TIME_INCREMENT
              CALL Field_ParameterSetsCopy(INDEPENDENT_FIELD_ALE_STOKES,FIELD_U_VARIABLE_TYPE, &
                & FIELD_MESH_DISPLACEMENT_SET_TYPE,FIELD_MESH_VELOCITY_SET_TYPE,ALPHA,err,error,*999)
            ELSE
              CALL FlagError("Mesh motion calculation not successful for ALE problem.",err,error,*999)
            END IF
          ELSE
            CALL FlagError("Mesh update is not defined for non-dynamic problems.",err,error,*999)
          END IF
        CASE DEFAULT
          localError="Problem subtype "//TRIM(NumberToVString(controlLoop%PROBLEM%specification(3),"*",err,error))// &
            & " is not valid for a Stokes equation fluid type of a fluid mechanics problem class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        CALL FlagError("Problem is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Solver is not associated.",err,error,*999)
    ENDIF
    
    EXITS("Stokes_PreSolveALEUpdateMesh")    
    RETURN
999 ERRORSEXITS("Stokes_PreSolveALEUpdateMesh",err,error)
    RETURN 1
    
  END SUBROUTINE Stokes_PreSolveALEUpdateMesh

  !
  !================================================================================================================================
  !
  !>Update mesh parameters for three component Laplace problem
  SUBROUTINE Stokes_PreSolveALEUpdateParameters(SOLVER,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(FieldType), POINTER :: independentField
    TYPE(SolverEquationsType), POINTER :: SOLVER_EQUATIONS  !<A pointer to the solver equations
    TYPE(SolverMappingType), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    TYPE(EquationsSetType), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    TYPE(EquationsType), POINTER :: EQUATIONS
    TYPE(VARYING_STRING) :: localError
    TYPE(DomainType), POINTER :: DOMAIN
    TYPE(FieldVariableType), POINTER :: fieldVariable
    TYPE(DomainNodesType), POINTER :: DOMAIN_NODES

    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT
    INTEGER(INTG) :: component_idx,node_idx,deriv_idx,local_ny,variable_idx,variable_type
    REAL(DP), POINTER :: MESH_STIFF_VALUES(:)


    ENTERS("Stokes_PreSolveALEUpdateParameters",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(SOLVER,controlLoop,ERR,ERROR,*999)
    IF(ASSOCIATED(controlLoop)) THEN
      CALL CONTROL_LOOP_CURRENT_TIMES_GET(controlLoop,CURRENT_TIME,TIME_INCREMENT,err,error,*999)
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(control_loop%problem%specification)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(control_loop%problem%specification,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a Stokes problem.",err,error,*999)
          END IF
          SELECT CASE(controlLoop%PROBLEM%specification(3))
            CASE(PROBLEM_STATIC_STOKES_SUBTYPE,PROBLEM_LAPLACE_STOKES_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_TRANSIENT_STOKES_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_ALE_STOKES_SUBTYPE)
              IF(SOLVER%solveType==SOLVER_LINEAR_TYPE) THEN
                !Get the independent field for the ALE Stokes problem
                SOLVER_EQUATIONS=>SOLVER%solverEquations
                IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  SOLVER_MAPPING=>SOLVER_equations%solverMapping
                  IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                    EQUATIONS_SET=>SOLVER_MAPPING%equationsSets(1)%ptr
                    NULLIFY(MESH_STIFF_VALUES)
                    CALL Field_ParameterSetDataGet(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,MESH_STIFF_VALUES,err,error,*999)
                    IF(ASSOCIATED(EQUATIONS_SET)) THEN
                      EQUATIONS=>SOLVER_MAPPING%equationsSetToSolverMatricesMap(1)%EQUATIONS
                      IF(ASSOCIATED(EQUATIONS)) THEN
                        independentField=>EQUATIONS_SET%INDEPENDENT%independentField
                        IF(ASSOCIATED(independentField)) THEN
                          DO variable_idx=1,EQUATIONS_SET%dependent%dependentField%numberOfVariables
                            variable_type=EQUATIONS_SET%dependent%dependentField%VARIABLES(variable_idx)%variableType
                            fieldVariable=>EQUATIONS_SET%dependent%dependentField%variableTypeMap(variable_type)%ptr
                            IF(ASSOCIATED(fieldVariable)) THEN
                              DO component_idx=1,fieldVariable%numberOfComponents
                                DOMAIN=>fieldVariable%COMPONENTS(component_idx)%DOMAIN
                                IF(ASSOCIATED(DOMAIN)) THEN
                                  IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                                    DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                                    IF(ASSOCIATED(DOMAIN_NODES)) THEN
                                      !Loop over the local nodes excluding the ghosts.
                                      DO node_idx=1,DOMAIN_NODES%numberOfNodes
                                        DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%numberOfDerivatives
                                          !Default to version 1 of each node derivative
                                          CALL FieldVariable_LocalNodeDOFGet(fieldVariable,1,deriv_idx,node_idx, &
                                            & component_idx,local_ny,err,error,*999)
           !                             !Calculation of K values dependent on current mesh topology
                                          MESH_STIFF_VALUES(local_ny)=1.0_DP
                                          CALL Field_ParameterSetUpdateLocalDOF(EQUATIONS_SET%INDEPENDENT% &
                                            & independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny, &
                                            & MESH_STIFF_VALUES(local_ny),err,error,*999)
                                        ENDDO !deriv_idx
                                      ENDDO !node_idx
                                    ENDIF
                                  ENDIF
                                ENDIF
                              ENDDO !component_idx
                            ENDIF
                          ENDDO !variable_idx
                        ELSE
                          CALL FlagError("Independent field is not associated.",err,error,*999)
                        END IF
                      ELSE
                        CALL FlagError("Equations are not associated.",err,error,*999)
                      END IF
                    ELSE
                      CALL FlagError("Equations set is not associated.",err,error,*999)
                    ENDIF
                    CALL Field_ParameterSetDataRestore(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,MESH_STIFF_VALUES,err,error,*999)
                  ELSE
                    CALL FlagError("Solver mapping is not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("Solver equations are not associated.",err,error,*999)
                END IF
              ELSE IF(SOLVER%solveType==SOLVER_DYNAMIC_TYPE) THEN
                CALL FlagError("Mesh motion calculation not successful for ALE problem.",err,error,*999)
              END IF
            CASE DEFAULT
              localError="Problem subtype "//TRIM(NumberToVString(controlLoop%PROBLEM%specification(3),"*",err,error))// &
                & " is not valid for a Stokes equation fluid type of a fluid mechanics problem class."
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
    
    EXITS("Stokes_PreSolveALEUpdateParameters")
    RETURN
999 ERRORSEXITS("Stokes_PreSolveALEUpdateParameters",err,error)
    RETURN 1
    
  END SUBROUTINE Stokes_PreSolveALEUpdateParameters

  !
  !================================================================================================================================
  !

  !>Output data post solve
  SUBROUTINE STOKES_POST_SOLVE_OUTPUT_DATA(controlLoop,SOLVER,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the control loop to solve.
    TYPE(SolverType), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(SolverEquationsType), POINTER :: SOLVER_EQUATIONS  !<A pointer to the solver equations
    TYPE(SolverMappingType), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    TYPE(EquationsSetType), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    TYPE(FieldsType), POINTER :: Fields
    TYPE(VARYING_STRING) :: localError,METHOD,FILENAME

    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT
    INTEGER(INTG) :: EQUATIONS_SET_IDX,CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER,numberOfDimensions
    LOGICAL :: EXPORT_FIELD
    CHARACTER(14) :: OUTPUT_FILE

    ENTERS("STOKES_POST_SOLVE_OUTPUT_DATA",err,error,*999)

    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVER_MAPPING)
    NULLIFY(EQUATIONS_SET)

    IF(ASSOCIATED(controlLoop)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(controlLoop%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(control_loop%problem%specification)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(control_loop%problem%specification,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a Stokes problem.",err,error,*999)
          END IF
          CALL SYSTEM('mkdir -p ./output')
          SELECT CASE(controlLoop%PROBLEM%specification(3))
          CASE(PROBLEM_STATIC_STOKES_SUBTYPE,PROBLEM_LAPLACE_STOKES_SUBTYPE)
            SOLVER_EQUATIONS=>SOLVER%solverEquations
            IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
              SOLVER_MAPPING=>SOLVER_equations%solverMapping
              IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                !Make sure the equations sets are up to date
                DO equations_set_idx=1,SOLVER_MAPPING%numberOfEquationsSets
                  EQUATIONS_SET=>SOLVER_MAPPING%equationsSets(equations_set_idx)%ptr
                  FILENAME="./output/"//"STATIC_SOLUTION"
                  METHOD="FORTRAN"
                  IF(SOLVER%outputType>=SOLVER_PROGRESS_OUTPUT) THEN
                    CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
                    CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Now export fields... ",err,error,*999)
                  ENDIF
                  Fields=>EQUATIONS_SET%REGION%fieldS
                  CALL FIELD_IO_NODES_EXPORT(Fields,FILENAME,METHOD,err,error,*999)
                  CALL FIELD_IO_ELEMENTS_EXPORT(Fields,FILENAME,METHOD,err,error,*999)
                  NULLIFY(Fields)
                ENDDO
              ENDIF
            ENDIF
          CASE(PROBLEM_TRANSIENT_STOKES_SUBTYPE,PROBLEM_ALE_STOKES_SUBTYPE,PROBLEM_PGM_STOKES_SUBTYPE)
            CALL CONTROL_LOOP_CURRENT_TIMES_GET(controlLoop,CURRENT_TIME,TIME_INCREMENT,err,error,*999)
            SOLVER_EQUATIONS=>SOLVER%solverEquations
            IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
              SOLVER_MAPPING=>SOLVER_equations%solverMapping
              IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                !Make sure the equations sets are up to date
                DO equations_set_idx=1,SOLVER_MAPPING%numberOfEquationsSets
                  EQUATIONS_SET=>SOLVER_MAPPING%equationsSets(equations_set_idx)%ptr
                  CURRENT_LOOP_ITERATION=controlLoop%timeLoop%iterationNumber
                  OUTPUT_ITERATION_NUMBER=controlLoop%timeLoop%outputNumber
                  IF(OUTPUT_ITERATION_NUMBER/=0) THEN
                    IF(controlLoop%timeLoop%currentTime<=controlLoop%timeLoop%stopTime) THEN
                      IF(CURRENT_LOOP_ITERATION<10) THEN
                        WRITE(OUTPUT_FILE,'("TIME_STEP_000",I0)') CURRENT_LOOP_ITERATION
                      ELSE IF(CURRENT_LOOP_ITERATION<100) THEN
                        WRITE(OUTPUT_FILE,'("TIME_STEP_00",I0)') CURRENT_LOOP_ITERATION
                      ELSE IF(CURRENT_LOOP_ITERATION<1000) THEN
                        WRITE(OUTPUT_FILE,'("TIME_STEP_0",I0)') CURRENT_LOOP_ITERATION
                      ELSE IF(CURRENT_LOOP_ITERATION<10000) THEN
                        WRITE(OUTPUT_FILE,'("TIME_STEP_",I0)') CURRENT_LOOP_ITERATION
                      END IF
                      FILENAME="./output/"//"MainTime_"//TRIM(NumberToVString(CURRENT_LOOP_ITERATION,"*",err,error))
                      METHOD="FORTRAN"
                      IF(MOD(CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER)==0)  THEN
                        IF(controlLoop%outputtype >= CONTROL_LOOP_PROGRESS_OUTPUT) THEN
                          CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
                          CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Now export fields... ",err,error,*999)
                        ENDIF
                        Fields=>EQUATIONS_SET%REGION%fieldS
                        CALL FIELD_IO_NODES_EXPORT(Fields,FILENAME,METHOD,err,error,*999)
                        CALL FIELD_IO_ELEMENTS_EXPORT(Fields,FILENAME,METHOD,err,error,*999)
                        NULLIFY(Fields)
                        IF(controlLoop%outputtype >= CONTROL_LOOP_PROGRESS_OUTPUT) THEN
                          CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,FILENAME,err,error,*999)
                          CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
                        ENDIF
                      END IF
                      IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
                        IF(EQUATIONS_SET%ANALYTIC%analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_4.OR. &
                          & EQUATIONS_SET%ANALYTIC%analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_5.OR. &
                          & EQUATIONS_SET%ANALYTIC%analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_4.OR. &
                          & EQUATIONS_SET%ANALYTIC%analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_5.OR. &
                          & EQUATIONS_SET%ANALYTIC%analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_1) THEN
                          CALL AnalyticAnalysis_Output(EQUATIONS_SET%dependent%dependentField,OUTPUT_FILE,err,error,*999)
                        ENDIF
                      ENDIF
                    ENDIF
                  ENDIF
                ENDDO
              ENDIF
            ENDIF
          CASE DEFAULT
            localError="Problem subtype "//TRIM(NumberToVString(controlLoop%PROBLEM%specification(3),"*",err,error))// &
              & " is not valid for a Stokes equation fluid type of a fluid mechanics problem class."
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
    EXITS("STOKES_POST_SOLVE_OUTPUT_DATA")
    RETURN
999 ERRORSEXITS("STOKES_POST_SOLVE_OUTPUT_DATA",err,error)
    RETURN 1
  END SUBROUTINE STOKES_POST_SOLVE_OUTPUT_DATA

  !
  !================================================================================================================================
  !

  !>Calculates the analytic solution and sets the boundary conditions for an analytic problem.
  SUBROUTINE Stokes_BoundaryConditionsAnalyticCalculate(EQUATIONS_SET,BOUNDARY_CONDITIONS,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: EQUATIONS_SET
    TYPE(BoundaryConditionsType), POINTER :: BOUNDARY_CONDITIONS
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
!\todo: Reduce number of variables used
    INTEGER(INTG) :: component_idx,deriv_idx,dim_idx,local_ny,node_idx,numberOfDimensions,variable_idx,variable_type,I,J,K
    INTEGER(INTG) :: number_of_nodes_xic(3),element_idx,en_idx,BOUND_COUNT,ANALYTIC_FUNCTION_TYPE,GLOBAL_DERIV_INDEX
    REAL(DP) :: VALUE,X(3),XI_COORDINATES(3)
!     REAL(DP) :: BOUNDARY_TOLERANCE, BOUNDARY_X(3,2),muParam,L
    REAL(DP) :: T_COORDINATES(20,3),CURRENT_TIME,muParam,rhoParam
    REAL(DP), POINTER :: GEOMETRIC_PARAMETERS(:)
    TYPE(DomainType), POINTER :: DOMAIN
    TYPE(DomainNodesType), POINTER :: DOMAIN_NODES
    TYPE(FieldType), POINTER :: dependentField,geometricField,materialsField
    TYPE(FieldVariableType), POINTER :: fieldVariable,GEOMETRIC_VARIABLE
    TYPE(FieldInterpolatedPointPtrType), POINTER :: INTERPOLATED_POINT(:)
    TYPE(FieldInterpolationParametersPtrType), POINTER :: INTERPOLATION_PARAMETERS(:)
!     TYPE(VARYING_STRING) :: localError

! ! !     !Temp variables
! ! !     INTEGER(INTG) :: number_of_element_nodes,temp_local_ny,temp_node_number,velocity_DOF_check,temp_local_node_number

    ENTERS("Stokes_BoundaryConditionsAnalyticCalculate",err,error,*999)
!\todo: Introduce user call to set parameters
    BOUND_COUNT=0
! ! ! !     L=10.0_DP
    XI_COORDINATES(3)=0.0_DP
!     BOUNDARY_TOLERANCE=0.000000001_DP
! ! !     BOUNDARY_X=0.0_DP
! ! !     T_COORDINATES=0.0_DP
! ! !     number_of_element_nodes=0
! ! !     temp_local_node_number=0
! ! !     temp_local_ny=0
! ! !     temp_node_number=0
! ! !     velocity_DOF_check=0
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
        dependentField=>EQUATIONS_SET%dependent%dependentField
        IF(ASSOCIATED(dependentField)) THEN
          geometricField=>EQUATIONS_SET%GEOMETRY%geometricField
          IF(ASSOCIATED(geometricField)) THEN
            NULLIFY(INTERPOLATION_PARAMETERS)
            NULLIFY(INTERPOLATED_POINT)
            CALL Field_InterpolationParametersInitialise(geometricField,INTERPOLATION_PARAMETERS,err,error,*999)
            CALL Field_InterpolatedPointsInitialise(INTERPOLATION_PARAMETERS,INTERPOLATED_POINT,err,error,*999)
            CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
! ! ! !\todo: Check adjacent element calculation / use boundary node flag instead / didn't work for simplex
! ! !             IF(numberOfDimensions==2) THEN
! ! !               BOUNDARY_X(1,1)=0.0_DP
! ! !               BOUNDARY_X(1,2)=10.0_DP
! ! !               BOUNDARY_X(2,1)=0.0_DP
! ! !               BOUNDARY_X(2,2)=10.0_DP
! ! !             ELSE IF(numberOfDimensions==3) THEN
! ! !               BOUNDARY_X(1,1)=-5.0_DP
! ! !               BOUNDARY_X(1,2)=5.0_DP
! ! !               BOUNDARY_X(2,1)=-5.0_DP
! ! !               BOUNDARY_X(2,2)=5.0_DP
! ! !               BOUNDARY_X(3,1)=-5.0_DP
! ! !               BOUNDARY_X(3,2)=5.0_DP
! ! !             ENDIF
            NULLIFY(GEOMETRIC_VARIABLE)
            CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,GEOMETRIC_VARIABLE,err,error,*999)
            NULLIFY(GEOMETRIC_PARAMETERS)
            CALL Field_ParameterSetDataGet(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS, &
              & err,error,*999)
            IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
              DO variable_idx=1,dependentField%numberOfVariables
                variable_type=dependentField%VARIABLES(variable_idx)%variableType
                fieldVariable=>dependentField%variableTypeMap(variable_type)%ptr
                IF(ASSOCIATED(fieldVariable)) THEN
                  CALL Field_ParameterSetCreate(dependentField,variable_type,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
                  DO component_idx=1,fieldVariable%numberOfComponents
                    BOUND_COUNT=0
                    IF(fieldVariable%COMPONENTS(component_idx)%interpolationType==FIELD_NODE_BASED_INTERPOLATION) THEN
                      DOMAIN=>fieldVariable%COMPONENTS(component_idx)%DOMAIN
                      IF(ASSOCIATED(DOMAIN)) THEN
                        IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                          DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                          IF(ASSOCIATED(DOMAIN_NODES)) THEN
                            !Loop over the local nodes excluding the ghosts.
                            DO node_idx=1,DOMAIN_NODES%numberOfNodes
                              element_idx=DOMAIN%topology%nodes%nodes(node_idx)%surroundingElements(1)
                              CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,element_idx, &
                                & INTERPOLATION_PARAMETERS(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                              en_idx=0
                              XI_COORDINATES=0.0_DP
                              number_of_nodes_xic(1)=DOMAIN%topology%elements%elements(element_idx)%basis%numberOfNodesXiC(1)
                              number_of_nodes_xic(2)=DOMAIN%topology%elements%elements(element_idx)%basis%numberOfNodesXiC(2)
                              IF(numberOfDimensions==3) THEN
                                number_of_nodes_xic(3)=DOMAIN%topology%elements%elements(element_idx)%basis%numberOfNodesXiC(3)
                              ELSE
                                number_of_nodes_xic(3)=1
                              ENDIF
  !\todo: Use boundary flag
                              IF(DOMAIN%topology%elements%maximumNumberOfElementParameters==4.AND.numberOfDimensions==2 .OR. &
                                & DOMAIN%topology%elements%maximumNumberOfElementParameters==9.OR. &
                                & DOMAIN%topology%elements%maximumNumberOfElementParameters==16.OR. &
                                & DOMAIN%topology%elements%maximumNumberOfElementParameters==8.OR. &
                                & DOMAIN%topology%elements%maximumNumberOfElementParameters==27.OR. &
                                & DOMAIN%topology%elements%maximumNumberOfElementParameters==64) THEN
                                DO K=1,number_of_nodes_xic(3)
                                  DO J=1,number_of_nodes_xic(2)
                                    DO I=1,number_of_nodes_xic(1)
                                      en_idx=en_idx+1
                                      IF(DOMAIN%topology%elements%elements(element_idx)%elementNodes(en_idx)==node_idx) EXIT
                                      XI_COORDINATES(1)=XI_COORDINATES(1)+(1.0_DP/(number_of_nodes_xic(1)-1))
                                    ENDDO
                                      IF(DOMAIN%topology%elements%elements(element_idx)%elementNodes(en_idx)==node_idx) EXIT
                                      XI_COORDINATES(1)=0.0_DP
                                      XI_COORDINATES(2)=XI_COORDINATES(2)+(1.0_DP/(number_of_nodes_xic(2)-1))
                                  ENDDO
                                  IF(DOMAIN%topology%elements%elements(element_idx)%elementNodes(en_idx)==node_idx) EXIT
                                  XI_COORDINATES(1)=0.0_DP
                                  XI_COORDINATES(2)=0.0_DP
                                  IF(number_of_nodes_xic(3)/=1) THEN
                                    XI_COORDINATES(3)=XI_COORDINATES(3)+(1.0_DP/(number_of_nodes_xic(3)-1))
                                  ENDIF
                                ENDDO
                                CALL Field_InterpolateXi(NO_PART_DERIV,XI_COORDINATES, &
                                  & INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                              ELSE
  !\todo: Use boundary flag
                                IF(DOMAIN%topology%elements%maximumNumberOfElementParameters==3) THEN
                                  T_COORDINATES(1,1:2)=[0.0_DP,1.0_DP]
                                  T_COORDINATES(2,1:2)=[1.0_DP,0.0_DP]
                                  T_COORDINATES(3,1:2)=[1.0_DP,1.0_DP]
                                ELSE IF(DOMAIN%topology%elements%maximumNumberOfElementParameters==6) THEN
                                  T_COORDINATES(1,1:2)=[0.0_DP,1.0_DP]
                                  T_COORDINATES(2,1:2)=[1.0_DP,0.0_DP]
                                  T_COORDINATES(3,1:2)=[1.0_DP,1.0_DP]
                                  T_COORDINATES(4,1:2)=[0.5_DP,0.5_DP]
                                  T_COORDINATES(5,1:2)=[1.0_DP,0.5_DP]
                                  T_COORDINATES(6,1:2)=[0.5_DP,1.0_DP]
                                ELSE IF(DOMAIN%topology%elements%maximumNumberOfElementParameters==10.AND. &
                                  & numberOfDimensions==2) THEN
                                  T_COORDINATES(1,1:2)=[0.0_DP,1.0_DP]
                                  T_COORDINATES(2,1:2)=[1.0_DP,0.0_DP]
                                  T_COORDINATES(3,1:2)=[1.0_DP,1.0_DP]
                                  T_COORDINATES(4,1:2)=[1.0_DP/3.0_DP,2.0_DP/3.0_DP]
                                  T_COORDINATES(5,1:2)=[2.0_DP/3.0_DP,1.0_DP/3.0_DP]
                                  T_COORDINATES(6,1:2)=[1.0_DP,1.0_DP/3.0_DP]
                                  T_COORDINATES(7,1:2)=[1.0_DP,2.0_DP/3.0_DP]
                                  T_COORDINATES(8,1:2)=[2.0_DP/3.0_DP,1.0_DP]
                                  T_COORDINATES(9,1:2)=[1.0_DP/3.0_DP,1.0_DP]
                                  T_COORDINATES(10,1:2)=[2.0_DP/3.0_DP,2.0_DP/3.0_DP]
                                ELSE IF(DOMAIN%topology%elements%maximumNumberOfElementParameters==4) THEN
                                  T_COORDINATES(1,1:3)=[0.0_DP,1.0_DP,1.0_DP]
                                  T_COORDINATES(2,1:3)=[1.0_DP,0.0_DP,1.0_DP]
                                  T_COORDINATES(3,1:3)=[1.0_DP,1.0_DP,0.0_DP]
                                  T_COORDINATES(4,1:3)=[1.0_DP,1.0_DP,1.0_DP]
                                ELSE IF(DOMAIN%topology%elements%maximumNumberOfElementParameters==10.AND. &
                                  & numberOfDimensions==3) THEN
                                  T_COORDINATES(1,1:3)=[0.0_DP,1.0_DP,1.0_DP]
                                  T_COORDINATES(2,1:3)=[1.0_DP,0.0_DP,1.0_DP]
                                  T_COORDINATES(3,1:3)=[1.0_DP,1.0_DP,0.0_DP]
                                  T_COORDINATES(4,1:3)=[1.0_DP,1.0_DP,1.0_DP]
                                  T_COORDINATES(5,1:3)=[0.5_DP,0.5_DP,1.0_DP]
                                  T_COORDINATES(6,1:3)=[0.5_DP,1.0_DP,0.5_DP]
                                  T_COORDINATES(7,1:3)=[0.5_DP,1.0_DP,1.0_DP]
                                  T_COORDINATES(8,1:3)=[1.0_DP,0.5_DP,0.5_DP]
                                  T_COORDINATES(9,1:3)=[1.0_DP,1.0_DP,0.5_DP]
                                  T_COORDINATES(10,1:3)=[1.0_DP,0.5_DP,1.0_DP]
                                ELSE IF(DOMAIN%topology%elements%maximumNumberOfElementParameters==20) THEN
                                  T_COORDINATES(1,1:3)=[0.0_DP,1.0_DP,1.0_DP]
                                  T_COORDINATES(2,1:3)=[1.0_DP,0.0_DP,1.0_DP]
                                  T_COORDINATES(3,1:3)=[1.0_DP,1.0_DP,0.0_DP]
                                  T_COORDINATES(4,1:3)=[1.0_DP,1.0_DP,1.0_DP]
                                  T_COORDINATES(5,1:3)=[1.0_DP/3.0_DP,2.0_DP/3.0_DP,1.0_DP]
                                  T_COORDINATES(6,1:3)=[2.0_DP/3.0_DP,1.0_DP/3.0_DP,1.0_DP]
                                  T_COORDINATES(7,1:3)=[1.0_DP/3.0_DP,1.0_DP,2.0_DP/3.0_DP]
                                  T_COORDINATES(8,1:3)=[2.0_DP/3.0_DP,1.0_DP,1.0_DP/3.0_DP]
                                  T_COORDINATES(9,1:3)=[1.0_DP/3.0_DP,1.0_DP,1.0_DP]
                                  T_COORDINATES(10,1:3)=[2.0_DP/3.0_DP,1.0_DP,1.0_DP]
                                  T_COORDINATES(11,1:3)=[1.0_DP,1.0_DP/3.0_DP,2.0_DP/3.0_DP]
                                  T_COORDINATES(12,1:3)=[1.0_DP,2.0_DP/3.0_DP,1.0_DP/3.0_DP]
                                  T_COORDINATES(13,1:3)=[1.0_DP,1.0_DP,1.0_DP/3.0_DP]
                                  T_COORDINATES(14,1:3)=[1.0_DP,1.0_DP,2.0_DP/3.0_DP]
                                  T_COORDINATES(15,1:3)=[1.0_DP,1.0_DP/3.0_DP,1.0_DP]
                                  T_COORDINATES(16,1:3)=[1.0_DP,2.0_DP/3.0_DP,1.0_DP]
                                  T_COORDINATES(17,1:3)=[2.0_DP/3.0_DP,2.0_DP/3.0_DP,2.0_DP/3.0_DP]
                                  T_COORDINATES(18,1:3)=[2.0_DP/3.0_DP,2.0_DP/3.0_DP,1.0_DP]
                                  T_COORDINATES(19,1:3)=[2.0_DP/3.0_DP,1.0_DP,2.0_DP/3.0_DP]
                                  T_COORDINATES(20,1:3)=[1.0_DP,2.0_DP/3.0_DP,2.0_DP/3.0_DP]
                                ENDIF
                                DO K=1,DOMAIN%topology%elements%maximumNumberOfElementParameters
                                  IF(DOMAIN%topology%elements%elements(element_idx)%elementNodes(K)==node_idx) EXIT
                                ENDDO
                                IF(numberOfDimensions==2) THEN
                                  CALL Field_InterpolateXi(NO_PART_DERIV,T_COORDINATES(K,1:2), &
                                    & INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                                ELSE IF(numberOfDimensions==3) THEN
                                  CALL Field_InterpolateXi(NO_PART_DERIV,T_COORDINATES(K,1:3), &
                                    & INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                                ENDIF
                              ENDIF
                              X=0.0_DP
                              DO dim_idx=1,numberOfDimensions
                                X(dim_idx)=INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(dim_idx,1)
                              ENDDO !dim_idx

                              !Loop over the derivatives
                              DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%numberOfDerivatives
                                ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET%ANALYTIC%analyticFunctionType
                                GLOBAL_DERIV_INDEX=DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%globalDerivativeIndex
                                CURRENT_TIME=0.0_DP
                                materialsField=>EQUATIONS_SET%MATERIALS%materialsField
                                !Define muParam, density=1
                                muParam=materialsField%variables(1)%parameterSets%parameterSets(1)%ptr% &
                                  & parameters%cmiss%dataDP(1)
                                !Define rhoParam, density=2
                                IF(ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_4.OR. &
                                  & ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_5.OR. &
                                  & ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_4.OR. &
                                  & ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_5) THEN
                                  rhoParam=materialsField%variables(1)%parameterSets%parameterSets(1)%ptr% &
                                    & parameters%cmiss%dataDP(2)
                                ELSE
                                  rhoParam=0.0_DP
                                ENDIF
                                CALL STOKES_EQUATION_ANALYTIC_FUNCTIONS(VALUE,X,muParam,rhoParam,CURRENT_TIME,variable_type, &
                                  & GLOBAL_DERIV_INDEX,ANALYTIC_FUNCTION_TYPE,numberOfDimensions, &
                                  & fieldVariable%numberOfComponents,component_idx,err,error,*999)
                                !Default to version 1 of each node derivative
                                CALL FieldVariable_LocalNodeDOFGet(fieldVariable,1,deriv_idx,node_idx, &
                                  & component_idx,local_ny,err,error,*999)
                                CALL Field_ParameterSetUpdateLocalDOF(dependentField,variable_type, &
                                  & FIELD_ANALYTIC_VALUES_SET_TYPE,local_ny,VALUE,err,error,*999)
                                IF(variable_type==FIELD_U_VARIABLE_TYPE) THEN
  ! \todo: This part should work even for simplex elements as soon as adjacent element calculation has been fixed
                                  IF(DOMAIN_NODES%NODES(node_idx)%boundaryNode) THEN
                                    !If we are a boundary node then set the analytic value on the boundary
                                    IF(component_idx<=numberOfDimensions) THEN
                                      CALL BoundaryConditions_SetLocalDOF(BOUNDARY_CONDITIONS,dependentField,variable_type, &
                                        & local_ny,BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
                                      BOUND_COUNT=BOUND_COUNT+1
                                    ELSE
  ! \todo: This is just a workaround for linear pressure fields in simplex element components
                                      IF(DOMAIN%topology%elements%maximumNumberOfElementParameters==3) THEN
                                        IF(ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_1.OR. &
                                          & ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_2.OR. &
                                          & ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_3.OR. &
                                          & ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_4.OR. &
                                          & ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_5) THEN
                                          IF(-0.001_DP<X(1).AND.X(1)<0.001_DP.AND.-0.001_DP<X(2).AND.X(2)<0.001_DP.OR. &
                                            &  10.0_DP-0.001_DP<X(1).AND.X(1)<10.0_DP+0.001_DP.AND.-0.001_DP<X(2).AND. &
                                            & X(2)<0.001_DP.OR. &
                                            &  10.0_DP-0.001_DP<X(1).AND.X(1)<10.0_DP+0.001_DP.AND.10.0_DP-0.001_DP<X(2).AND. &
                                            & X(2)<10.0_DP+0.001_DP.OR. &
                                            &  -0.001_DP<X(1).AND.X(1)<0.001_DP.AND.10.0_DP-0.001_DP<X(2).AND. &
                                            & X(2)<10.0_DP+0.001_DP) THEN
                                              CALL BoundaryConditions_SetLocalDOF(BOUNDARY_CONDITIONS,dependentField, &
                                                & variable_type,local_ny,BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
                                              BOUND_COUNT=BOUND_COUNT+1
                                          ENDIF
                                        ENDIF
                                      ELSE IF(DOMAIN%topology%elements%maximumNumberOfElementParameters==4.AND. &
                                        & numberOfDimensions==3) THEN
                                        IF(ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_1.OR. &
                                          & ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_2.OR. &
                                          & ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_3.OR. &
                                          & ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_4.OR. &
                                          & ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_5) THEN
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
                                            CALL BoundaryConditions_SetLocalDOF(BOUNDARY_CONDITIONS,dependentField, &
                                              & variable_type,local_ny,BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
                                            BOUND_COUNT=BOUND_COUNT+1
                                          ENDIF
                                        ENDIF
  ! \todo: This is how it should be if adjacent elements would be working
                                      ELSE IF(BOUND_COUNT==0) THEN
                                        CALL BoundaryConditions_SetLocalDOF(BOUNDARY_CONDITIONS,dependentField,variable_type, &
                                          & local_ny,BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
                                        BOUND_COUNT=BOUND_COUNT+1
                                      ENDIF


                                    ENDIF
                                  ELSE
                                    IF(component_idx<=numberOfDimensions) THEN
                                      CALL Field_ParameterSetUpdateLocalDOF(dependentField,variable_type, &
                                        & FIELD_VALUES_SET_TYPE,local_ny,VALUE,err,error,*999)
                                    ENDIF
                                  ENDIF
  ! \todo: Use boundary node flag
  ! ! !                                 !If we are a boundary node then set the analytic value on the boundary
  ! ! !                                 IF(numberOfDimensions==2) THEN
  ! ! !                                   IF(X(1)<BOUNDARY_X(1,1)+BOUNDARY_TOLERANCE.AND.X(1)>BOUNDARY_X(1,1)-BOUNDARY_TOLERANCE.OR. &
  ! ! !                                     & X(1)<BOUNDARY_X(1,2)+BOUNDARY_TOLERANCE.AND.X(1)>BOUNDARY_X(1,2)-BOUNDARY_TOLERANCE.OR. &
  ! ! !                                     & X(2)<BOUNDARY_X(2,1)+BOUNDARY_TOLERANCE.AND.X(2)>BOUNDARY_X(2,1)-BOUNDARY_TOLERANCE.OR. &
  ! ! !                                     & X(2)<BOUNDARY_X(2,2)+BOUNDARY_TOLERANCE.AND.X(2)>BOUNDARY_X(2,2)-BOUNDARY_TOLERANCE) THEN
  ! ! !                                     IF(component_idx<=numberOfDimensions) THEN
  ! ! !                                       CALL BoundaryConditions_SetLocalDOF(BOUNDARY_CONDITIONS,variable_type,local_ny, &
  ! ! !                                         & BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
  ! ! !                                     BOUND_COUNT=BOUND_COUNT+1
  ! ! !                                     !Apply boundary conditions check for pressure nodes
  ! ! !                                     ELSE IF(component_idx>numberOfDimensions) THEN
  ! ! !                                       IF(DOMAIN%topology%elements%maximumNumberOfElementParameters==4) THEN
  ! ! !                                       IF(X(1)<BOUNDARY_X(1,1)+BOUNDARY_TOLERANCE.AND.X(1)>BOUNDARY_X(1,1)-BOUNDARY_TOLERANCE.AND. &
  ! ! !                                         & X(2)<BOUNDARY_X(2,1)+BOUNDARY_TOLERANCE.AND.X(2)>BOUNDARY_X(2,1)-BOUNDARY_TOLERANCE) &
  ! ! !                                         & THEN
  ! ! !                                            ! Commented out for testing purposes
  ! ! !                                           CALL BoundaryConditions_SetLocalDOF(BOUNDARY_CONDITIONS,variable_type,local_ny, &
  ! ! !                                             & BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
  ! ! !                                           BOUND_COUNT=BOUND_COUNT+1
  ! ! !                                       ENDIF
  ! ! !                                       ENDIF
  ! ! ! !\todo: Again, ...
  ! ! !                                       IF(DOMAIN%topology%elements%maximumNumberOfElementParameters==3.OR. &
  ! ! !                                         & DOMAIN%topology%elements%maximumNumberOfElementParameters==6.OR. &
  ! ! !                                         & DOMAIN%topology%elements%maximumNumberOfElementParameters==10) THEN
  ! ! !                                       IF(X(1)<BOUNDARY_X(1,1)+BOUNDARY_TOLERANCE.AND.X(1)>BOUNDARY_X(1,1)-BOUNDARY_TOLERANCE.AND. &
  ! ! !                                         & X(2)<BOUNDARY_X(2,1)+BOUNDARY_TOLERANCE.AND.X(2)>BOUNDARY_X(2,1)-BOUNDARY_TOLERANCE.OR. &
  ! ! !                                         & X(1)<BOUNDARY_X(1,1)+BOUNDARY_TOLERANCE.AND.X(1)>BOUNDARY_X(1,1)-BOUNDARY_TOLERANCE.AND.&
  ! ! !                                         & X(2)<BOUNDARY_X(2,2)+BOUNDARY_TOLERANCE.AND.X(2)>BOUNDARY_X(2,2)-BOUNDARY_TOLERANCE.OR. &
  ! ! !                                         & X(1)<BOUNDARY_X(1,2)+BOUNDARY_TOLERANCE.AND.X(1)>BOUNDARY_X(1,2)-BOUNDARY_TOLERANCE.AND.&
  ! ! !                                         & X(2)<BOUNDARY_X(2,1)+BOUNDARY_TOLERANCE.AND.X(2)>BOUNDARY_X(2,1)-BOUNDARY_TOLERANCE.OR. &
  ! ! !                                         & X(1)<BOUNDARY_X(1,2)+BOUNDARY_TOLERANCE.AND.X(1)>BOUNDARY_X(1,2)-BOUNDARY_TOLERANCE.AND.&
  ! ! !                                         & X(2)<BOUNDARY_X(2,2)+BOUNDARY_TOLERANCE.AND.X(2)>BOUNDARY_X(2,2)-BOUNDARY_TOLERANCE) &
  ! ! !                                         & THEN
  ! ! !                                           CALL BoundaryConditions_SetLocalDOF(BOUNDARY_CONDITIONS,variable_type,local_ny, &
  ! ! !                                             & BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
  ! ! !                                           BOUND_COUNT=BOUND_COUNT+1
  ! ! !                                       ENDIF
  ! ! !                                       ENDIF
  ! ! !                                     ENDIF
  ! ! !                                   ENDIF
  ! ! !                                     IF(component_idx<=numberOfDimensions+1) THEN
  ! ! !                                       CALL Field_ParameterSetUpdateLocalDOF(dependentField,variable_type, &
  ! ! !                                         & FIELD_VALUES_SET_TYPE,local_ny,VALUE,err,error,*999)
  ! ! !                                     ENDIF
  ! ! !                                 ELSE IF(numberOfDimensions==3) THEN
  ! ! !                                   IF(X(1)<BOUNDARY_X(1,1)+BOUNDARY_TOLERANCE.AND.X(1)>BOUNDARY_X(1,1)-BOUNDARY_TOLERANCE.OR. &
  ! ! !                                     & X(1)<BOUNDARY_X(1,2)+BOUNDARY_TOLERANCE.AND.X(1)>BOUNDARY_X(1,2)-BOUNDARY_TOLERANCE.OR. &
  ! ! !                                     & X(2)<BOUNDARY_X(2,1)+BOUNDARY_TOLERANCE.AND.X(2)>BOUNDARY_X(2,1)-BOUNDARY_TOLERANCE.OR. &
  ! ! !                                     & X(2)<BOUNDARY_X(2,2)+BOUNDARY_TOLERANCE.AND.X(2)>BOUNDARY_X(2,2)-BOUNDARY_TOLERANCE.OR. &
  ! ! !                                     & X(3)<BOUNDARY_X(3,1)+BOUNDARY_TOLERANCE.AND.X(3)>BOUNDARY_X(3,1)-BOUNDARY_TOLERANCE.OR. &
  ! ! !                                     & X(3)<BOUNDARY_X(3,2)+BOUNDARY_TOLERANCE.AND.X(3)>BOUNDARY_X(3,2)-BOUNDARY_TOLERANCE) THEN
  ! ! !                                     IF(component_idx<=numberOfDimensions) THEN
  ! ! !                                       CALL BoundaryConditions_SetLocalDOF(BOUNDARY_CONDITIONS,variable_type,local_ny, &
  ! ! !                                         & BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
  ! ! !                                     BOUND_COUNT=BOUND_COUNT+1
  ! ! !                                     !Apply boundary conditions check for pressure nodes
  ! ! !                                     ELSE IF(component_idx>numberOfDimensions) THEN
  ! ! !                                       IF(DOMAIN%topology%elements%maximumNumberOfElementParameters==4.OR. &
  ! ! !                                         & DOMAIN%topology%elements%maximumNumberOfElementParameters==10.OR. &
  ! ! !                                         & DOMAIN%topology%elements%maximumNumberOfElementParameters==20) THEN
  ! ! !                                       IF(X(1)<BOUNDARY_X(1,1)+BOUNDARY_TOLERANCE.AND.X(1)>BOUNDARY_X(1,1)-BOUNDARY_TOLERANCE.AND. &
  ! ! !                                        & X(2)<BOUNDARY_X(2,1)+BOUNDARY_TOLERANCE.AND.X(2)>BOUNDARY_X(2,1)-BOUNDARY_TOLERANCE.AND. &
  ! ! !                                        & X(3)<BOUNDARY_X(3,1)+BOUNDARY_TOLERANCE.AND.X(3)>BOUNDARY_X(3,1)-BOUNDARY_TOLERANCE.OR. &
  ! ! !                                        & X(1)<BOUNDARY_X(1,1)+BOUNDARY_TOLERANCE.AND.X(1)>BOUNDARY_X(1,1)-BOUNDARY_TOLERANCE.AND. &
  ! ! !                                        & X(2)<BOUNDARY_X(2,1)+BOUNDARY_TOLERANCE.AND.X(2)>BOUNDARY_X(2,1)-BOUNDARY_TOLERANCE.AND. &
  ! ! !                                        & X(3)<BOUNDARY_X(3,2)+BOUNDARY_TOLERANCE.AND.X(3)>BOUNDARY_X(3,2)-BOUNDARY_TOLERANCE.OR. &
  ! ! !                                        & X(1)<BOUNDARY_X(1,1)+BOUNDARY_TOLERANCE.AND.X(1)>BOUNDARY_X(1,1)-BOUNDARY_TOLERANCE.AND. &
  ! ! !                                        & X(2)<BOUNDARY_X(2,2)+BOUNDARY_TOLERANCE.AND.X(2)>BOUNDARY_X(2,2)-BOUNDARY_TOLERANCE.AND. &
  ! ! !                                        & X(3)<BOUNDARY_X(3,1)+BOUNDARY_TOLERANCE.AND.X(3)>BOUNDARY_X(3,1)-BOUNDARY_TOLERANCE.OR. &
  ! ! !                                        & X(1)<BOUNDARY_X(1,1)+BOUNDARY_TOLERANCE.AND.X(1)>BOUNDARY_X(1,1)-BOUNDARY_TOLERANCE.AND. &
  ! ! !                                        & X(2)<BOUNDARY_X(2,2)+BOUNDARY_TOLERANCE.AND.X(2)>BOUNDARY_X(2,2)-BOUNDARY_TOLERANCE.AND. &
  ! ! !                                        & X(3)<BOUNDARY_X(3,2)+BOUNDARY_TOLERANCE.AND.X(3)>BOUNDARY_X(3,2)-BOUNDARY_TOLERANCE.OR. &
  ! ! !                                        & X(1)<BOUNDARY_X(1,2)+BOUNDARY_TOLERANCE.AND.X(1)>BOUNDARY_X(1,2)-BOUNDARY_TOLERANCE.AND. &
  ! ! !                                        & X(2)<BOUNDARY_X(2,1)+BOUNDARY_TOLERANCE.AND.X(2)>BOUNDARY_X(2,1)-BOUNDARY_TOLERANCE.AND. &
  ! ! !                                        & X(3)<BOUNDARY_X(3,1)+BOUNDARY_TOLERANCE.AND.X(3)>BOUNDARY_X(3,1)-BOUNDARY_TOLERANCE.OR. &
  ! ! !                                        & X(1)<BOUNDARY_X(1,2)+BOUNDARY_TOLERANCE.AND.X(1)>BOUNDARY_X(1,2)-BOUNDARY_TOLERANCE.AND. &
  ! ! !                                        & X(2)<BOUNDARY_X(2,1)+BOUNDARY_TOLERANCE.AND.X(2)>BOUNDARY_X(2,1)-BOUNDARY_TOLERANCE.AND. &
  ! ! !                                        & X(3)<BOUNDARY_X(3,2)+BOUNDARY_TOLERANCE.AND.X(3)>BOUNDARY_X(3,2)-BOUNDARY_TOLERANCE.OR. &
  ! ! !                                        & X(1)<BOUNDARY_X(1,2)+BOUNDARY_TOLERANCE.AND.X(1)>BOUNDARY_X(1,2)-BOUNDARY_TOLERANCE.AND. &
  ! ! !                                        & X(2)<BOUNDARY_X(2,2)+BOUNDARY_TOLERANCE.AND.X(2)>BOUNDARY_X(2,2)-BOUNDARY_TOLERANCE.AND. &
  ! ! !                                        & X(3)<BOUNDARY_X(3,1)+BOUNDARY_TOLERANCE.AND.X(3)>BOUNDARY_X(3,1)-BOUNDARY_TOLERANCE.OR. &
  ! ! !                                        & X(1)<BOUNDARY_X(1,2)+BOUNDARY_TOLERANCE.AND.X(1)>BOUNDARY_X(1,2)-BOUNDARY_TOLERANCE.AND. &
  ! ! !                                        & X(2)<BOUNDARY_X(2,2)+BOUNDARY_TOLERANCE.AND.X(2)>BOUNDARY_X(2,2)-BOUNDARY_TOLERANCE.AND. &
  ! ! !                                        & X(3)<BOUNDARY_X(3,2)+BOUNDARY_TOLERANCE.AND.X(3)>BOUNDARY_X(3,2)-BOUNDARY_TOLERANCE) THEN
  ! ! !                                          CALL BoundaryConditions_SetLocalDOF(BOUNDARY_CONDITIONS,variable_type,local_ny, &
  ! ! !                                            & BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
  ! ! !                                          BOUND_COUNT=BOUND_COUNT+1
  ! ! !                                       ENDIF
  ! ! !                                       ENDIF
  ! ! !                                     ENDIF
  ! ! !                                   ELSE
  ! ! !                                     IF(component_idx<=numberOfDimensions+1) THEN
  ! ! !                                       CALL Field_ParameterSetUpdateLocalDOF(dependentField,variable_type, &
  ! ! !                                         & FIELD_VALUES_SET_TYPE,local_ny,VALUE,err,error,*999)
  ! ! !                                     ENDIF
  ! ! !                                   ENDIF
  ! ! !                                 ENDIF
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
                  CALL Field_ParameterSetUpdateStart(dependentField,variable_type,FIELD_ANALYTIC_VALUES_SET_TYPE, &
                    & err,error,*999)
                  CALL Field_ParameterSetUpdateFinish(dependentField,variable_type,FIELD_ANALYTIC_VALUES_SET_TYPE, &
                    & err,error,*999)
                  CALL Field_ParameterSetUpdateStart(dependentField,variable_type,FIELD_VALUES_SET_TYPE, &
                    & err,error,*999)
                  CALL Field_ParameterSetUpdateFinish(dependentField,variable_type,FIELD_VALUES_SET_TYPE, &
                    & err,error,*999)
                ELSE
                  CALL FlagError("Field variable is not associated.",err,error,*999)
                ENDIF
              ENDDO !variable_idx
              CALL Field_ParameterSetDataRestore(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & GEOMETRIC_PARAMETERS,err,error,*999)
              CALL Field_InterpolatedPointsFinalise(INTERPOLATED_POINT,err,error,*999)
              CALL Field_InterpolationParametersFinalise(INTERPOLATION_PARAMETERS,err,error,*999)
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

    EXITS("Stokes_BoundaryConditionsAnalyticCalculate")
    RETURN
999 ERRORSEXITS("Stokes_BoundaryConditionsAnalyticCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE Stokes_BoundaryConditionsAnalyticCalculate

  !
  !================================================================================================================================
  !
  !>Calculates the various analytic solutions given X and time, can be called from within analytic calculate or elsewhere if needed
  SUBROUTINE STOKES_EQUATION_ANALYTIC_FUNCTIONS(VALUE,X,muParam,rhoParam,CURRENT_TIME,VARIABLE_TYPE, &
    & GLOBAL_DERIV_INDEX,ANALYTIC_FUNCTION_TYPE,numberOfDimensions,NUMBER_OF_COMPONENTS,COMPONENT_IDX,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    REAL(DP), INTENT(OUT) :: VALUE
    REAL(DP) :: muParam,rhoParam
    REAL(DP), INTENT(IN) :: CURRENT_TIME
    REAL(DP), INTENT(IN), DIMENSION(3) :: X
    INTEGER(INTG), INTENT(IN) :: numberOfDimensions,NUMBER_OF_COMPONENTS,COMPONENT_IDX
    !Local variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: variable_type,GLOBAL_DERIV_INDEX,ANALYTIC_FUNCTION_TYPE
    !TYPE(DomainType), POINTER :: DOMAIN
    !TYPE(DomainNodesType), POINTER :: DOMAIN_NODES
    REAL(DP) :: INTERNAL_TIME

    ENTERS("STOKES_EQUATION_ANALYTIC_FUNCTIONS",err,error,*999)

!\todo: Introduce user-defined or default values instead for density and viscosity
    INTERNAL_TIME=CURRENT_TIME
     SELECT CASE(ANALYTIC_FUNCTION_TYPE)
       CASE(EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_1)
         IF(numberOfDimensions==2.AND.NUMBER_OF_COMPONENTS==3) THEN
           !Polynomial function
           SELECT CASE(variable_type)
             CASE(FIELD_U_VARIABLE_TYPE)
               SELECT CASE(GLOBAL_DERIV_INDEX)
                 CASE(NO_GLOBAL_DERIV)
                   IF(component_idx==1) THEN
                     !calculate u
                     VALUE=X(2)**2/10.0_DP**2
                   ELSE IF(component_idx==2) THEN
                     !calculate v
                     VALUE=X(1)**2/10.0_DP**2
                   ELSE IF(component_idx==3) THEN
                     !calculate p
                     VALUE=2.0_DP*muParam/10.0_DP**2*X(1)
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
                     & GLOBAL_DERIV_INDEX,"*",err,error))// &
                     & " is invalid."
                   CALL FlagError(localError,err,error,*999)
               END SELECT
             CASE(FIELD_DELUDELN_VARIABLE_TYPE)
               SELECT CASE(GLOBAL_DERIV_INDEX)
                 CASE(NO_GLOBAL_DERIV)
                   VALUE= 0.0_DP
                 CASE(GLOBAL_DERIV_S1)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE(GLOBAL_DERIV_S2)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE(GLOBAL_DERIV_S1_S2)
                   CALL FlagError("Not implemented.",err,error,*999)
                 CASE DEFAULT
                   localError="The global derivative index of "//TRIM(NumberToVString( &
                     & GLOBAL_DERIV_INDEX,"*",err,error))// &
                     & " is invalid."
                   CALL FlagError(localError,err,error,*999)
               END SELECT
             CASE DEFAULT
               localError="The variable type of "//TRIM(NumberToVString(variable_type,"*",err,error))// &
                 & " is invalid."
               CALL FlagError(localError,err,error,*999)
           END SELECT
         ELSE
           localError="The number of components does not correspond to the number of dimensions."
           CALL FlagError(localError,err,error,*999)
         ENDIF
       CASE(EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_2)
         IF(numberOfDimensions==2.AND.NUMBER_OF_COMPONENTS==3) THEN
           !Exponential function
           SELECT CASE(variable_type)
             CASE(FIELD_U_VARIABLE_TYPE)
               SELECT CASE(GLOBAL_DERIV_INDEX)
                 CASE(NO_GLOBAL_DERIV)
                   IF(component_idx==1) THEN
                     !calculate u
                     VALUE= EXP((X(1)-X(2))/10.0_DP)
                   ELSE IF(component_idx==2) THEN
                     !calculate v
                     VALUE= EXP((X(1)-X(2))/10.0_DP)
                   ELSE IF(component_idx==3) THEN
                     !calculate p
                     VALUE= 2.0_DP*muParam/10.0_DP*EXP((X(1)-X(2))/10.0_DP)
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
                     & GLOBAL_DERIV_INDEX,"*",err,error))// &
                     & " is invalid."
                   CALL FlagError(localError,err,error,*999)
               END SELECT
             CASE(FIELD_DELUDELN_VARIABLE_TYPE)
               SELECT CASE(GLOBAL_DERIV_INDEX)
                 CASE(NO_GLOBAL_DERIV)
                   IF(component_idx==1) THEN
                     !calculate u
                     VALUE= 0.0_DP
                   ELSE IF(component_idx==2) THEN
                     !calculate v
                     VALUE= 0.0_DP
                   ELSE IF(component_idx==3) THEN
                     !calculate p
                     VALUE= 0.0_DP
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
                     & GLOBAL_DERIV_INDEX,"*",err,error))// &
                     & " is invalid."
                   CALL FlagError(localError,err,error,*999)
               END SELECT
             CASE DEFAULT
               localError="The variable type of "//TRIM(NumberToVString(variable_type,"*",err,error))// &
                 & " is invalid."
               CALL FlagError(localError,err,error,*999)
           END SELECT
         ELSE
           localError="The number of components does not correspond to the number of dimensions."
           CALL FlagError(localError,err,error,*999)
         ENDIF
       CASE(EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_3)
         IF(numberOfDimensions==2.AND.NUMBER_OF_COMPONENTS==3) THEN
           !Sine and cosine functions
           SELECT CASE(variable_type)
             CASE(FIELD_U_VARIABLE_TYPE)
               SELECT CASE(GLOBAL_DERIV_INDEX)
                 CASE(NO_GLOBAL_DERIV)
                   IF(component_idx==1) THEN
                     !calculate u
                     VALUE=SIN(2.0_DP*PI*X(1)/10.0_DP)*SIN(2.0_DP*PI*X(2)/10.0_DP)
                   ELSE IF(component_idx==2) THEN
                     !calculate v
                     VALUE=COS(2.0_DP*PI*X(1)/10.0_DP)*COS(2.0_DP*PI*X(2)/10.0_DP)
                   ELSE IF(component_idx==3) THEN
                     !calculate p
                     VALUE=4.0_DP*muParam*PI/10.0_DP*SIN(2.0_DP*PI*X(2)/10.0_DP)*COS(2.0_DP*PI*X(1)/10.0_DP)
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
                     & GLOBAL_DERIV_INDEX,"*",err,error))// &
                     & " is invalid."
                   CALL FlagError(localError,err,error,*999)
               END SELECT
             CASE(FIELD_DELUDELN_VARIABLE_TYPE)
               SELECT CASE(GLOBAL_DERIV_INDEX)
                 CASE(NO_GLOBAL_DERIV)
                   IF(component_idx==1) THEN
                     !calculate u
                     VALUE=0.0_DP
                   ELSE IF(component_idx==2) THEN
                     !calculate v
                     VALUE=16.0_DP*muParam*PI**2/10.0_DP**2*cos(2.0_DP*PI*X(2)/10.0_DP)*cos(2.0_DP*PI*X(1)/10.0_DP)
                   ELSE IF(component_idx==3) THEN
                     !calculate p
                     VALUE=0.0_DP
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
                     & GLOBAL_DERIV_INDEX,"*",err,error))// &
                     & " is invalid."
                   CALL FlagError(localError,err,error,*999)
               END SELECT
             CASE DEFAULT
               localError="The variable type of "//TRIM(NumberToVString(variable_type,"*",err,error))// &
                 & " is invalid."
               CALL FlagError(localError,err,error,*999)
           END SELECT
         ELSE
           localError="The number of components does not correspond to the number of dimensions."
           CALL FlagError(localError,err,error,*999)
         ENDIF
       CASE(EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_4)
         IF(numberOfDimensions==2.AND.NUMBER_OF_COMPONENTS==3) THEN
           !Reduced Taylor-Green solution for Stokes
           CALL FlagError("Not implemented.",err,error,*999)
         ENDIF
       CASE(EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_5)
         IF(numberOfDimensions==2.AND.NUMBER_OF_COMPONENTS==3) THEN
           !Stokes-Taylor-Green dynamic
           SELECT CASE(variable_type)
             CASE(FIELD_U_VARIABLE_TYPE)
               SELECT CASE(GLOBAL_DERIV_INDEX)
                 CASE(NO_GLOBAL_DERIV)
                   IF(component_idx==1) THEN
                     !calculate u
                     VALUE=X(2)*exp(-(2.0_DP*muParam/rhoParam*CURRENT_TIME))
                   ELSE IF(component_idx==2) THEN
                     !calculate v
                     VALUE=X(1)*exp(-(2.0_DP*muParam/rhoParam*CURRENT_TIME))
                   ELSE IF(component_idx==3) THEN
                     !calculate p
                     VALUE=2.0_DP*X(2)*muParam*exp(-(2.0_DP*muParam/rhoParam*CURRENT_TIME))*X(1)
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
                     & GLOBAL_DERIV_INDEX,"*",err,error))// &
                     & " is invalid."
                   CALL FlagError(localError,err,error,*999)
               END SELECT
             CASE(FIELD_DELUDELN_VARIABLE_TYPE)
               SELECT CASE(GLOBAL_DERIV_INDEX)
                 CASE(NO_GLOBAL_DERIV)
                   IF(component_idx==1) THEN
                     !calculate u
                     VALUE=0.0_DP
                   ELSE IF(component_idx==2) THEN
                     !calculate v
                     VALUE=0.0_DP
                   ELSE IF(component_idx==3) THEN
                     !calculate p
                     VALUE=0.0_DP
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
                     & GLOBAL_DERIV_INDEX,"*",err,error))// &
                     & " is invalid."
                   CALL FlagError(localError,err,error,*999)
               END SELECT
             CASE DEFAULT
               localError="The variable type of "//TRIM(NumberToVString(variable_type,"*",err,error))// &
                 & " is invalid."
               CALL FlagError(localError,err,error,*999)
           END SELECT
         ELSE
           localError="The number of components does not correspond to the number of dimensions."
           CALL FlagError(localError,err,error,*999)
         ENDIF
       CASE(EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_1)
         IF(numberOfDimensions==3.AND.NUMBER_OF_COMPONENTS==4) THEN
!POLYNOM
           SELECT CASE(variable_type)
             CASE(FIELD_U_VARIABLE_TYPE)
               SELECT CASE(GLOBAL_DERIV_INDEX)
                 CASE(NO_GLOBAL_DERIV)
                   IF(component_idx==1) THEN
                     !calculate u
                     VALUE=X(2)**2/10.0_DP**2+X(3)**2/10.0_DP**2
                   ELSE IF(component_idx==2) THEN
                     !calculate v
                     VALUE=X(1)**2/10.0_DP**2+X(3)**2/10.0_DP**2
                   ELSE IF(component_idx==3) THEN
                     !calculate w
                     VALUE=X(1)**2/10.0_DP**2+X(2)**2/10.0_DP**2
                   ELSE IF(component_idx==4) THEN
                     !calculate p
                     VALUE=4.0_DP*muParam/10.0_DP**2*X(1)
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
                     & GLOBAL_DERIV_INDEX,"*",err,error))// &
                     & " is invalid."
                   CALL FlagError(localError,err,error,*999)
               END SELECT
             CASE(FIELD_DELUDELN_VARIABLE_TYPE)
               SELECT CASE(GLOBAL_DERIV_INDEX)
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
                     & GLOBAL_DERIV_INDEX,"*",err,error))// &
                     & " is invalid."
                   CALL FlagError(localError,err,error,*999)
                END SELECT
              CASE DEFAULT
                localError="The variable type of "//TRIM(NumberToVString(variable_type,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
           END SELECT
         ELSE
           localError="The number of components does not correspond to the number of dimensions."
           CALL FlagError(localError,err,error,*999)
         ENDIF
       CASE(EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_2)
         IF(numberOfDimensions==3.AND.NUMBER_OF_COMPONENTS==4) THEN
           !Exponential function
           SELECT CASE(variable_type)
             CASE(FIELD_U_VARIABLE_TYPE)
               SELECT CASE(GLOBAL_DERIV_INDEX)
                 CASE(NO_GLOBAL_DERIV)
                   IF(component_idx==1) THEN
                     !calculate u
                     VALUE=EXP((X(1)-X(2))/10.0_DP)+EXP((X(3)-X(1))/10.0_DP)
                   ELSE IF(component_idx==2) THEN
                     !calculate v
                     VALUE=EXP((X(1)-X(2))/10.0_DP)+EXP((X(2)-X(3))/10.0_DP)
                   ELSE IF(component_idx==3) THEN
                     !calculate w
                     VALUE=EXP((X(3)-X(1))/10.0_DP)+EXP((X(2)-X(3))/10.0_DP)
                   ELSE IF(component_idx==4) THEN
                     !calculate p
                     VALUE=2.0_DP*muParam/10.0_DP*(EXP((X(1)-X(2))/10.0_DP)-EXP((X(3)-X(1))/10.0_DP))
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
                     & GLOBAL_DERIV_INDEX,"*",err,error))// &
                     & " is invalid."
                   CALL FlagError(localError,err,error,*999)
               END SELECT
             CASE(FIELD_DELUDELN_VARIABLE_TYPE)
               SELECT CASE(GLOBAL_DERIV_INDEX)
                 CASE(NO_GLOBAL_DERIV)
                   IF(component_idx==1) THEN
                     !calculate u
                     VALUE=0.0_DP
                   ELSE IF(component_idx==2) THEN
                     !calculate v
                     VALUE=-2.0_DP*muParam*(2.0_DP*EXP(X(1)-X(2))+EXP(X(2)-X(3)))
                   ELSE IF(component_idx==3) THEN
                     !calculate w
                     VALUE=-2.0_DP*muParam*(2.0_DP*EXP(X(3)-X(1))+EXP(X(2)-X(3)))
                   ELSE IF(component_idx==4) THEN
                     !calculate p
                     VALUE=0.0_DP
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
                     & GLOBAL_DERIV_INDEX,"*",err,error))// &
                     & " is invalid."
                   CALL FlagError(localError,err,error,*999)
               END SELECT
             CASE DEFAULT
               localError="The variable type of "//TRIM(NumberToVString(variable_type,"*",err,error))// &
                 & " is invalid."
               CALL FlagError(localError,err,error,*999)
           END SELECT
         ELSE
           localError="The number of components does not correspond to the number of dimensions."
           CALL FlagError(localError,err,error,*999)
         ENDIF
       CASE(EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_3)
         IF(numberOfDimensions==3.AND.NUMBER_OF_COMPONENTS==4) THEN
           !Sine and cosine functions
           SELECT CASE(variable_type)
             CASE(FIELD_U_VARIABLE_TYPE)
               SELECT CASE(GLOBAL_DERIV_INDEX)
                 CASE(NO_GLOBAL_DERIV)
                   IF(component_idx==1) THEN
                     !calculate u
                     VALUE=sin(2.0_DP*PI*X(1)/10.0_DP)*sin(2.0_DP*PI*X(2)/10.0_DP)*sin(2.0_DP*PI*X(3)/10.0_DP)
                   ELSE IF(component_idx==2) THEN
                     !calculate v
                     VALUE=2.0_DP*cos(2.0_DP*PI*x(1)/10.0_DP)*sin(2.0_DP*PI*X(3)/10.0_DP)*cos(2.0_DP*PI*X(2)/10.0_DP)
                   ELSE IF(component_idx==3) THEN
                     !calculate w
                     VALUE=-cos(2.0_DP*PI*X(1)/10.0_DP)*sin(2.0_DP*PI*X(2)/10.0_DP)*cos(2.0_DP*PI*X(3)/10.0_DP)
                   ELSE IF(component_idx==4) THEN
                     !calculate p
                     VALUE=6.0_DP*muParam*PI/10.0_DP*sin(2.0_DP*PI*X(2)/10.0_DP)*sin(2.0_DP*PI*X(3)/10.0_DP)* &
                       & cos(2.0_DP*PI*X(1)/10.0_DP)
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
                     & GLOBAL_DERIV_INDEX,"*",err,error))// &
                     & " is invalid."
                   CALL FlagError(localError,err,error,*999)
                END SELECT
              CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                SELECT CASE(GLOBAL_DERIV_INDEX)
                  CASE(NO_GLOBAL_DERIV)
                    IF(component_idx==1) THEN
                      !calculate u
                      VALUE=0.0_DP
                    ELSE IF(component_idx==2) THEN
                      !calculate v
                      VALUE=36*muParam*PI**2/10.0_DP**2*cos(2.0_DP*PI*X(2)/10.0_DP)*sin(2.0_DP*PI*X(3)/10.0_DP)* &
                        & cos(2.0_DP*PI*X(1)/10.0_DP)
                    ELSE IF(component_idx==3) THEN
                      !calculate w
                      VALUE=0.0_DP
                    ELSE IF(component_idx==4) THEN
                      !calculate p
                      VALUE=0.0_DP
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
                      & GLOBAL_DERIV_INDEX,"*",err,error))// &
                      & " is invalid."
                    CALL FlagError(localError,err,error,*999)
                END SELECT
              CASE DEFAULT
                localError="The variable type of "//TRIM(NumberToVString(variable_type,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
            END SELECT
          ELSE
            localError="The number of components does not correspond to the number of dimensions."
            CALL FlagError(localError,err,error,*999)
          ENDIF
       CASE(EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_4)
         IF(numberOfDimensions==3.AND.NUMBER_OF_COMPONENTS==4) THEN
           !Reduced Taylor-Green solution for Stokes
           CALL FlagError("Not implemented.",err,error,*999)
         ENDIF
       CASE(EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_5)
         IF(numberOfDimensions==3.AND.NUMBER_OF_COMPONENTS==4) THEN
           !Stokes-Taylor-Green dynamic
           SELECT CASE(variable_type)
             CASE(FIELD_U_VARIABLE_TYPE)
               SELECT CASE(GLOBAL_DERIV_INDEX)
                 CASE(NO_GLOBAL_DERIV)
                   IF(component_idx==1) THEN
                     !calculate u
                     VALUE=X(2)*exp(-(2.0_DP*muParam/rhoParam*CURRENT_TIME))
                   ELSE IF(component_idx==2) THEN
                     !calculate v
                     VALUE=X(1)*exp(-(2.0_DP*muParam/rhoParam*CURRENT_TIME))
                   ELSE IF(component_idx==3) THEN
                     !calculate v
                     VALUE=0.0_DP
                   ELSE IF(component_idx==4) THEN
                     !calculate p
                     VALUE=2.0_DP*X(2)*muParam*exp(-(2.0_DP*muParam/rhoParam*CURRENT_TIME))*X(1)
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
                     & GLOBAL_DERIV_INDEX,"*",err,error))// &
                     & " is invalid."
                   CALL FlagError(localError,err,error,*999)
               END SELECT
             CASE(FIELD_DELUDELN_VARIABLE_TYPE)
               SELECT CASE(GLOBAL_DERIV_INDEX)
                 CASE(NO_GLOBAL_DERIV)
                   IF(component_idx==1) THEN
                     !calculate u
                     VALUE=0.0_DP
                   ELSE IF(component_idx==2) THEN
                     !calculate v
                     VALUE=0.0_DP
                   ELSE IF(component_idx==3) THEN
                     !calculate p
                     VALUE=0.0_DP
                   ELSE IF(component_idx==4) THEN
                     !calculate p
                     VALUE=0.0_DP
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
                     & GLOBAL_DERIV_INDEX,"*",err,error))// &
                     & " is invalid."
                   CALL FlagError(localError,err,error,*999)
               END SELECT
             CASE DEFAULT
               localError="The variable type of "//TRIM(NumberToVString(variable_type,"*",err,error))// &
                 & " is invalid."
               CALL FlagError(localError,err,error,*999)
           END SELECT
         ELSE
           localError="The number of components does not correspond to the number of dimensions."
           CALL FlagError(localError,err,error,*999)
         ENDIF
        CASE DEFAULT
          localError="The analytic function type of "// &
            & TRIM(NumberToVString(ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT

    EXITS("STOKES_EQUATION_ANALYTIC_FUNCTIONS")
    RETURN
999 ERRORSEXITS("STOKES_EQUATION_ANALYTIC_FUNCTIONS",err,error)
    RETURN 1

  END SUBROUTINE STOKES_EQUATION_ANALYTIC_FUNCTIONS

  !
  !================================================================================================================================
  !

END MODULE StokesEquationsRoutines
