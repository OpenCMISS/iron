!> \file
!> \brief This module handles all Darcy equations routines.
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
!> Contributor(s): Adam Reeve, Chris Bradley
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

!>This module handles all Darcy equations routines.
MODULE DarcyEquationsRoutines  

  USE BaseRoutines
  USE BasisRoutines
  USE BasisAccessRoutines
  USE BoundaryConditionsRoutines
  USE BoundaryConditionAccessRoutines
  USE Constants
  USE ControlLoopRoutines
  USE ControlLoopAccessRoutines
  USE ComputationRoutines
  USE ComputationAccessRoutines
  USE CoordinateSystemAccessRoutines
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
  USE FiniteElasticityRoutines
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
  USE SolverAccessRoutines
  USE Timer
  USE Types

#include "macros.h"

  IMPLICIT NONE

  PUBLIC Darcy_EquationsSetSetup
  
  PUBLIC Darcy_EquationsSetSpecificationSet
  
  PUBLIC Darcy_EquationsSetSolutionMethodSet
  
  PUBLIC Darcy_BoundaryConditionsAnalyticCalculate

  PUBLIC Darcy_ProblemSetup
  
  PUBLIC Darcy_ProblemSpecificationSet

  PUBLIC Darcy_FiniteElementCalculate

  PUBLIC Darcy_PreSolve
  
  PUBLIC Darcy_PostSolve
  
  PUBLIC Darcy_PostSolveOutputData

  PUBLIC Darcy_PreLoop

  PUBLIC Darcy_PreSolveStorePreviousIterate

  PUBLIC Darcy_MonitorConvergence

  LOGICAL :: idebug1, idebug2, idebug3

CONTAINS

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a Darcy equation type of a fluid mechanics equations set class.
  SUBROUTINE Darcy_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: solutionMethod !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: esSpecification(3)
    TYPE(VARYING_STRING) :: localError

    ENTERS("Darcy_EquationsSetSolutionMethodSet",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
    
    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE, &
      & EQUATIONS_SET_QUASISTATIC_DARCY_SUBTYPE, &
      & EQUATIONS_SET_ALE_DARCY_SUBTYPE, &
      & EQUATIONS_SET_TRANSIENT_DARCY_SUBTYPE, &
      & EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
      & EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE, &
      & EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
      & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
      & EQUATIONS_SET_MULTI_COMPARTMENT_DARCY_SUBTYPE, &
      & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
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
      localError="Equations set subtype of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is not valid for a Darcy equation type of a fluid mechanics equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Darcy_EquationsSetSolutionMethodSet")
    RETURN
999 ERRORSEXITS("Darcy_EquationsSetSolutionMethodSet",err,error)
    RETURN 1
    
  END SUBROUTINE Darcy_EquationsSetSolutionMethodSet

  !
  !================================================================================================================================
  !

  !>Sets up the Darcy equation.
  SUBROUTINE Darcy_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to setup
    TYPE(EquationsSetSetupType), INTENT(INOUT) :: equationsSetSetup !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx,darcyDependentNumberOfComponents,darcyNumberOfComponents,dependentNumberOfComponents, &
      & dependentNumberOfVariables,elasticityDependentNumberOfComponents,equationsSetNumberOfComponents, &
      & equationsSetNumberOfVariables,esSpecification(3),geometricComponentNumber,geometricMeshComponent,geometricScalingType, &
      & independentNumberOfComponents,independentNumberOfVariables,lumpingType,materialFieldNumberOfComponents, &
      & materialsNumberOfUComponents,materialsNumberOfU1Components,materialsNumberOfVComponents,materialsNumberOfVariables, &
      & meshComponent,myMatrixIdx,numberOfCompartments,numberOfComponents,numberOfDimensions,solutionMethod, &
      & sourceNumberOfComponents,sparsityType,variableCount,variableIdx
    INTEGER(INTG), POINTER :: equationsSetFieldData(:)
    INTEGER(INTG), ALLOCATABLE :: couplingMatrixStorageType(:),couplingMatrixStructureType(:),variableTypes(:),variableUTypes(:)
    TYPE(DecompositionType), POINTER :: geometricDecomposition
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsSetAnalyticType), POINTER :: equationsAnalytic
    TYPE(EquationsSetEquationsFieldType), POINTER :: equationsField
    TYPE(EquationsSetIndependentType), POINTER :: equationsIndependent
    TYPE(EquationsSetMaterialsType), POINTER :: equationsMaterials
    TYPE(EquationsSetSourceType), POINTER :: equationsSource
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsSetEquationsFieldType), POINTER :: eqsEquationsSetField
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: equationsSetField,geometricField
    TYPE(FieldVariableType), POINTER :: fieldVariable,geometricVariable
    TYPE(RegionType), POINTER :: region
    TYPE(VARYING_STRING) :: localError

    ENTERS("Darcy_EquationsSetSetup",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)

    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE, EQUATIONS_SET_QUASISTATIC_DARCY_SUBTYPE, EQUATIONS_SET_ALE_DARCY_SUBTYPE, &
      & EQUATIONS_SET_TRANSIENT_DARCY_SUBTYPE, EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
      & EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE, EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
      & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE,EQUATIONS_SET_MULTI_COMPARTMENT_DARCY_SUBTYPE, &
      & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
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
      SELECT CASE(esSpecification(3))
      CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE,EQUATIONS_SET_QUASISTATIC_DARCY_SUBTYPE,EQUATIONS_SET_ALE_DARCY_SUBTYPE, &
        & EQUATIONS_SET_TRANSIENT_DARCY_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
        & EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
        & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE)
        SELECT CASE(equationsSetSetup%actionType)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
          CALL Darcy_EquationsSetSolutionMethodSet(equationsSet,EQUATIONS_SET_FEM_SOLUTION_METHOD,err,error,*999)
        CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
          !do nothing
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
            & " is invalid for a standard or quasistatic Darcy equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_MULTI_COMPARTMENT_DARCY_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
        equationsSetNumberOfVariables = 1
        equationsSetNumberOfComponents = 2
        NULLIFY(equationsField)
        CALL EquationsSet_EquationsFieldGet(equationsSet,equationsField,err,error,*999)
        SELECT CASE(equationsSetSetup%actionType)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
          CALL Darcy_EquationsSetSolutionMethodSet(equationsSet,EQUATIONS_SET_FEM_SOLUTION_METHOD,err,error,*999)
          IF(equationsField%equationsSetFieldAutoCreated) THEN
            !Create the auto created equations set field
            CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsField%equationsSetField,err,error,*999)
            NULLIFY(equationsSetField)
            CALL EquationsSet_EquationsSetFieldGet(equationsSet,equationsSetField,err,error,*999)
            CALL Field_LabelSet(equationsSetField,"Equations Set Field",err,error,*999)
            CALL Field_TypeSetAndLock(equationsSetField,FIELD_GENERAL_TYPE,err,error,*999)
            CALL Field_DependentTypeSetAndLock(equationsSetField,FIELD_INDEPENDENT_TYPE,err,error,*999)
            CALL Field_NumberOfVariablesSet(equationsSetField,equationsSetNumberOfVariables,err,error,*999)
            CALL Field_VariableTypesSetAndLock(equationsSetField,[FIELD_U_VARIABLE_TYPE],err,error,*999)
            CALL Field_DimensionSetAndLock(equationsSetField,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DataTypeSetAndLock(equationsSetField,FIELD_U_VARIABLE_TYPE,FIELD_INTG_TYPE,err,error,*999)
            CALL Field_NumberOfComponentsSetAndLock(equationsSetField,FIELD_U_VARIABLE_TYPE,equationsSetNumberOfComponents, &
              & err,error,*999)
          ELSE
            !Check the user specified field
            CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
            CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
            CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,equationsSetNumberOfVariables,err,error,*999)
            CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE],err,error,*999)
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_INTG_TYPE,err,error,*999)
            CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE, &
              & equationsSetNumberOfComponents,err,error,*999)
          ENDIF
        CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
          IF(equationsField%equationsSetFieldAutoCreated) THEN
            CALL Field_CreateFinish(equationsField%equationsSetField,err,error,*999)
            CALL Field_ComponentValuesInitialise(equationsField%equationsSetField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & 1,1_INTG,err,error,*999)
            CALL Field_ComponentValuesInitialise(equationsField%equationsSetField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & 2,1_INTG,err,error,*999)
          ENDIF
!!TODO: Check valid setup
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
            & " is invalid for a standard or quasistatic Darcy equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      END SELECT

    CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
      !-----------------------------------------------------------------
      ! g e o m e t r y   f i e l d
      !-----------------------------------------------------------------
      SELECT CASE(esSpecification(3))
      CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE,EQUATIONS_SET_QUASISTATIC_DARCY_SUBTYPE, &
        & EQUATIONS_SET_TRANSIENT_DARCY_SUBTYPE)
        !Do nothing
      CASE(EQUATIONS_SET_ALE_DARCY_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
        & EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
        & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE,EQUATIONS_SET_MULTI_COMPARTMENT_DARCY_SUBTYPE, &
        & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)        
        SELECT CASE(equationsSetSetup%actionType)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
          NULLIFY(geometricField)
          CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
          NULLIFY(geometricVariable)
          CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,geometricVariable,err,error,*999)
          CALL FieldVariable_ParameterSetEnsureCreated(geometricVariable,FIELD_INITIAL_VALUES_SET_TYPE,err,error,*999)
          CALL FieldVariable_ParameterSetEnsureCreated(geometricVariable,FIELD_PREVIOUS_VALUES_SET_TYPE,err,error,*999)
          CALL FieldVariable_ParameterSetEnsureCreated(geometricVariable,FIELD_MESH_DISPLACEMENT_SET_TYPE,err,error,*999)
          CALL FieldVariable_ParameterSetEnsureCreated(geometricVariable,FIELD_MESH_VELOCITY_SET_TYPE,err,error,*999)
          CALL FieldVariable_ParameterSetEnsureCreated(geometricVariable,FIELD_NEGATIVE_MESH_VELOCITY_SET_TYPE,err,error,*999)
          IF(esSpecification(3)==EQUATIONS_SET_MULTI_COMPARTMENT_DARCY_SUBTYPE .OR. &
            esSpecification(3)==EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE) THEN
            !Create the equations set field for multi-compartment Darcy
            equationsSetNumberOfComponents = 2
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
              DO componentIdx=1,equationsSetNumberOfComponents
                CALL Field_ComponentMeshComponentSetAndLock(equationsSetField,FIELD_U_VARIABLE_TYPE,componentIdx, &
                  & geometricComponentNumber,err,error,*999)
                CALL Field_ComponentInterpolationSetAndLock(equationsSetField,FIELD_U_VARIABLE_TYPE,componentIdx, &
                  & FIELD_CONSTANT_INTERPOLATION,err,error,*999)
              ENDDO !componentIdx
              !Default the field scaling to that of the geometric field
              CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
              CALL Field_ScalingTypeSet(equationsSetField,geometricScalingType,err,error,*999)
            ELSE
              !Do nothing
            ENDIF
          ENDIF
        CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
          ! do nothing
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
            & " is invalid for a linear diffusion equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      END SELECT
    CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
      !-----------------------------------------------------------------
      ! d e p e n d e n t   f i e l d
      !-----------------------------------------------------------------
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
      SELECT CASE(esSpecification(3))
      CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE, EQUATIONS_SET_QUASISTATIC_DARCY_SUBTYPE, EQUATIONS_SET_ALE_DARCY_SUBTYPE, &
        & EQUATIONS_SET_TRANSIENT_DARCY_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
        & EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
        & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
        & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
        dependentNumberOfVariables = 2  ! U and the normal component of its flux
        dependentNumberOfComponents = numberOfDimensions + 1
        SELECT CASE(equationsSetSetup%actionType)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
          IF(equationsSet%dependent%dependentFieldAutoCreated) THEN
            !Create the auto created dependent field
            CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsSet%dependent%dependentField,err,error,*999)
            CALL Field_TypeSetAndLock(equationsSet%dependent%dependentField,FIELD_GEOMETRIC_GENERAL_TYPE,err,error,*999)
            CALL Field_DependentTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_DEPENDENT_TYPE,err,error,*999)
            NULLIFY(geometricDecomposition)
            CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
            CALL Field_DecompositionSetAndLock(equationsSet%dependent%dependentField,geometricDecomposition, err,error,*999)
            CALL Field_GeometricFieldSetAndLock(equationsSet%dependent%dependentField,geometricField,err,error,*999)
            CALL Field_NumberOfVariablesSetAndLock(equationsSet%dependent%dependentField,dependentNumberOfVariables,err,error,*999)
            CALL Field_VariableTypesSetAndLock(equationsSet%dependent%dependentField,[FIELD_U_VARIABLE_TYPE, &
              & FIELD_DELUDELN_VARIABLE_TYPE],err,error,*999)
            CALL Field_VariableLabelSet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,"U",err,error,*999)
            CALL Field_VariableLabelSet(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,"del U/del n", &
              & err,error,*999)
            CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
              & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_DP_TYPE,err,error,*999)
            CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
              & FIELD_DP_TYPE,err,error,*999)            
            CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField, FIELD_U_VARIABLE_TYPE, &
              & dependentNumberOfComponents,err,error,*999)
            CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField, &
              & FIELD_DELUDELN_VARIABLE_TYPE, dependentNumberOfComponents,err,error,*999)
            !Default to the geometric interpolation setup
            CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricMeshComponent,err,error,*999)
            DO componentIdx=1,dependentNumberOfComponents
              IF(componentIdx<dependentNumberOfComponents ) THEN
                !Set velocity mesh component (default to the geometric one)
                meshComponent = geometricMeshComponent
                CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
                  & componentIdx,meshComponent,err,error,*999)
                CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & componentIdx,meshComponent,err,error,*999)
              ELSE
                !Set pressure mesh component (default to the geometric one)
                meshComponent = geometricMeshComponent
                CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
                  & componentIdx,meshComponent,err,error,*999)
                CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & componentIdx,meshComponent,err,error,*999)
              ENDIF
            ENDDO !componentIdx
            CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
            SELECT CASE(solutionMethod)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              DO componentIdx=1,dependentNumberOfComponents
                CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField, &
                  & FIELD_U_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              ENDDO !componentIdx
              !Default the scaling to the geometric field scaling
              CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
              CALL Field_ScalingTypeSet(equationsSet%dependent%dependentField,geometricScalingType,err,error,*999)
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
              localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))// &
                & " is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ELSE
            SELECT CASE(esSpecification(3))
            CASE(EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
              & EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
              & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE)
              !-----------------------------------------------------------------------
              ! Check the shared dependent field set up in finite elasticity routines
              !-----------------------------------------------------------------------
              CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GEOMETRIC_GENERAL_TYPE,err,error,*999)
              CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_DEPENDENT_TYPE,err,error,*999)
              CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,4,err,error,*999)
              CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_V_VARIABLE_TYPE,FIELD_DELVDELN_VARIABLE_TYPE],err,error,*999)
              CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                & err,error,*999)
              CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_DELVDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                & err,error,*999)
              CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_DELVDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              SELECT CASE(esSpecification(3))
              CASE(EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE)  !compressible elasticity
                elasticityDependentNumberOfComponents = numberOfDimensions
                darcyDependentNumberOfComponents = numberOfDimensions + 2  !(u,v,w,p,m)
              CASE(EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE)
                elasticityDependentNumberOfComponents = numberOfDimensions + 1
                darcyDependentNumberOfComponents = numberOfDimensions + 1
              CASE(EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE)
                elasticityDependentNumberOfComponents = numberOfDimensions + 1 !(u1,u2,u3,p)
                darcyDependentNumberOfComponents = numberOfDimensions + 1 !(u,v,w,m)
              END SELECT
              CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE, &
                & elasticityDependentNumberOfComponents,err,error,*999)
              CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE, &
                & elasticityDependentNumberOfComponents,err,error,*999)
              CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE, &
                & darcyDependentNumberOfComponents,err,error,*999)
              CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_DELVDELN_VARIABLE_TYPE, &
                & darcyDependentNumberOfComponents,err,error,*999)
              CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)              
              SELECT CASE(solutionMethod)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                !Mind that elastic hydrostatic pressure might be interpolated element-wise
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_DELVDELN_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
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
            CASE(EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
              !-----------------------------------------------------------------------
              ! Check the shared dependent field set up in finite elasticity routines
              ! Must have 2+2*numberOfCompartments number of variable types
              !-----------------------------------------------------------------------
              CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GEOMETRIC_GENERAL_TYPE,err,error,*999)
              CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_DEPENDENT_TYPE,err,error,*999)
              !Get the number of Darcy compartments from the equations set field
              NULLIFY(equationsSetField)
              CALL EquationsSet_EquationsSetFieldFieldGet(equationsSet,equationsSetField,err,error,*999)
              CALL Field_ParameterSetDataGet(equationsSetField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,equationsSetFieldData,err,error,*999)
              numberOfCompartments=equationsSetFieldData(2)
              CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,(2+2*numberOfCompartments),err,error,*999)
              ALLOCATE(variableTypes(2*numberOfCompartments+2))
              DO variableIdx=1,numberOfCompartments+1
                variableTypes(2*variableIdx-1)=FIELD_U_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(variableIdx-1))
                variableTypes(2*variableIdx)=FIELD_DELUDELN_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(variableIdx-1))
              ENDDO !variableIdx
              CALL Field_VariableTypesCheck(equationsSetSetup%field,variableTypes,err,error,*999)
              numberOfComponents=numberOfDimensions+1
              darcyNumberOfComponents=numberOfDimensions+1              
              DO variableIdx=1,2*numberOfCompartments+2
                CALL Field_DimensionCheck(equationsSetSetup%field,variableTypes(variableIdx),FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL Field_DataTypeCheck(equationsSetSetup%field,variableTypes(variableIdx),FIELD_DP_TYPE,err,error,*999)
                CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,variableTypes(variableIdx),numberOfComponents, &
                  & err,error,*999)
              ENDDO !variableIdx
              CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
              SELECT CASE(solutionMethod)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                !Elasticity:
                DO componentIdx=1,numberOfDimensions
                  CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,componentIdx, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,componentIdx,&
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                ENDDO !componentIdx
                !If solid hydrostatic pressure is driving Darcy flow, check that pressure uses node based interpolation
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfComponents,&
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & numberOfComponents,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                DO variableIdx=3,2*numberOfCompartments+2
                  !Darcy:
                  DO componentIdx=1,darcyNumberOfComponents
                    CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,variableTypes(variableIdx),componentIdx, &
                      & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  ENDDO !componentIdx
                ENDDO !variableIdx
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
                localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE(EQUATIONS_SET_MULTI_COMPARTMENT_DARCY_SUBTYPE)
              !Check the field created by Darcy routines for the multi-compartment model
              CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
              CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_DEPENDENT_TYPE,err,error,*999)
              NULLIFY(equationsSetField)
              CALL EquationsSet_EquationSetFieldGet(equationsSet,equationsSetField,err,error,*999)
              CALL Field_ParameterSetDataGet(equationsSetField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & equationsSetFieldData,err,error,*999)
              numberOfCompartments=equationsSetFieldData(2)
              CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,2*numberOfCompartments,err,error,*999)
              !Create & populate array storing all of the relevant variable types against which to check the field variables
              ALLOCATE(variableTypes(2*numberOfCompartments))
              DO variableIdx=1,numberOfCompartments
                variableTypes(2*variableIdx-1)=FIELD_U_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(variableIdx-1))
                variableTypes(2*variableIdx)=FIELD_DELUDELN_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(variableIdx-1))
              ENDDO !variableIdx
              CALL Field_VariableTypesCheck(equationsSetSetup%field,variableTypes,err,error,*999)
              DO variableIdx=1,2*numberOfCompartments
                CALL Field_DimensionCheck(equationsSetSetup%field,variableTypes(variableIdx),FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL Field_DataTypeCheck(equationsSetSetup%field,variableTypes(variableIdx),FIELD_DP_TYPE,err,error,*999)
                CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,variableTypes(variableIdx),numberOfDimensions+1, &
                  & err,error,*999)
              ENDDO !variableIdx
              CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
              SELECT CASE(solutionMethod)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                componentIdx=1
                DO variableIdx=1,2*numberOfCompartments
                  DO componentIdx=1,numberOfDimensions+1
                    CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,variableTypes(variableIdx),componentIdx, &
                      & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    !NOTE-pressure might use element based interpolation - need to account for this
                  ENDDO !componentIdx
                ENDDO !variableIdx
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
              !--------------------------------
              ! Check the user specified field
              !--------------------------------
              CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
              CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_DEPENDENT_TYPE,err,error,*999)
              CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,2,err,error,*999)
              CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE],&
                & err,error,*999)
              CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                & err,error,*999)
              CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                & err,error,*999)
              CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_NumberOfComponentsGet(equationsSet%GEOMETRY%geometricField, FIELD_U_VARIABLE_TYPE, &
                & numberOfDimensions,err,error,*999)
              dependentNumberOfComponents = numberOfDimensions + 1
              CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE, &
                & dependentNumberOfComponents,err,error,*999)
              CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE, &
                & dependentNumberOfComponents,err,error,*999)

              SELECT CASE(equationsSet%solutionMethod)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
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
                localError="The solution method of "//TRIM(NumberToVString(equationsSet%solutionMethod,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            END SELECT ! on (equationsSet%SPECIFICATION(3))
          ENDIF ! on (equationsSet%dependent%dependentFieldAutoCreated)
        CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
          IF(equationsSet%dependent%dependentFieldAutoCreated) THEN
            CALL Field_CreateFinish(equationsSet%dependent%dependentField,err,error,*999)
            CALL Field_ParameterSetEnsureCreated(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_INITIAL_VALUES_SET_TYPE,err,error,*999)
          ENDIF
          IF(esSpecification(3)/=EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)THEN
            !Actually, only needed for PGM (for elasticity_Darcy defined in elasticity V var):
            CALL Field_ParameterSetEnsureCreated(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_RELATIVE_VELOCITY_SET_TYPE,err,error,*999)
          ENDIF
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
            & " is invalid for a standard, quasistatic or ALE Darcy equation"
          CALL FlagError(localError,err,error,*999)
        END SELECT
      END SELECT
    CASE(EQUATIONS_SET_SETUP_INDEPENDENT_TYPE)
      !-----------------------------------------------------------------
      ! I N d e p e n d e n t   f i e l d
      !-----------------------------------------------------------------
      SELECT CASE(esSpecification(3))
      CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE, EQUATIONS_SET_QUASISTATIC_DARCY_SUBTYPE, EQUATIONS_SET_ALE_DARCY_SUBTYPE, &
        & EQUATIONS_SET_TRANSIENT_DARCY_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
        & EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
        & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
        & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
        !\todo: revise: do they all need an independent field ?
        NULLIFY(equationsIndependent)
        CALL EquationsSet_IndependentGet(equationsSet,equationsIndependent,err,error,*999)
        NULLIFY(geometricField)
        CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
        CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
        independentNumberOfVariables = 2  ! U and the normal component of its flux
        independentNumberOfComponents = numberOfDimensions !+ 1
        SELECT CASE(equationsSetSetup%actionType)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
          IF(equationsIndependent%independentFieldAutoCreated) THEN
            !Create the auto created INdependent field
            CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsIndependent%independentField,err,error,*999)
            CALL Field_TypeSetAndLock(equationsIndependent%independentField,FIELD_GENERAL_TYPE,err,error,*999)
            CALL Field_DependentTypeSetAndLock(equationsIndependent%independentField,FIELD_INDEPENDENT_TYPE,err,error,*999)
            NULLIFY(geometricDecomposition)
            CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
            CALL Field_DecompositionSetAndLock(equationsIndependent%independentField,geometricDecomposition,err,error,*999)
            CALL Field_GeometricFieldSetAndLock(equationsIndependent%independentField,geometricField,err,error,*999)
            CALL Field_NumberOfVariablesSetAndLock(equationsIndependent%independentField, &
              & independentNumberOfVariables,err,error,*999)
            CALL Field_VariableTypesSetAndLock(equationsIndependent%independentField,[FIELD_U_VARIABLE_TYPE, &
              & FIELD_DELUDELN_VARIABLE_TYPE],err,error,*999)
            CALL Field_VariableLabelSet(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE,"Independent U", &
              & err,error,*999)
            CALL Field_VariableLabelSet(equationsIndependent%independentField,FIELD_DELUDELN_VARIABLE_TYPE, &
              & "Independent del U/del n",err,error,*999)
            CALL Field_DimensionSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DimensionSetAndLock(equationsIndependent%independentField,FIELD_DELUDELN_VARIABLE_TYPE, &
              & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DataTypeSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_DP_TYPE,err,error,*999)
            CALL Field_DataTypeSetAndLock(equationsIndependent%independentField,FIELD_DELUDELN_VARIABLE_TYPE, &
              & FIELD_DP_TYPE,err,error,*999)
            CALL Field_NumberOfComponentsSetAndLock(equationsIndependent%independentField, FIELD_U_VARIABLE_TYPE, &
              & independentNumberOfComponents,err,error,*999)
            CALL Field_NumberOfComponentsSetAndLock(equationsIndependent%independentField, &
              & FIELD_DELUDELN_VARIABLE_TYPE,independentNumberOfComponents,err,error,*999)
            !Default to the geometric interpolation setup
            CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricMeshComponent,err,error,*999)
            DO componentIdx=1,independentNumberOfComponents
              IF(componentIdx<independentNumberOfComponents) THEN
                !Set velocity mesh component (default to the geometric one)
                meshComponent = geometricMeshComponent
                CALL Field_ComponentMeshComponentSet(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
                  & componentIdx,meshComponent,err,error,*999)
                CALL Field_ComponentMeshComponentSet(equationsIndependent%independentField,&
                  & FIELD_DELUDELN_VARIABLE_TYPE,componentIdx,meshComponent,err,error,*999)
              ELSE
                !Set pressure mesh component (default to the geometric one)
                meshComponent = geometricMeshComponent
                CALL Field_ComponentMeshComponentSet(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
                  & componentIdx,meshComponent,err,error,*999)
                CALL Field_ComponentMeshComponentSet(equationsIndependent%independentField,&
                  & FIELD_DELUDELN_VARIABLE_TYPE,componentIdx,meshComponent,err,error,*999)
              ENDIF
            ENDDO !componentIdx
            CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)            
            SELECT CASE(solutionMethod)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              DO componentIdx=1,independentNumberOfComponents
                CALL Field_ComponentInterpolationSetAndLock(equationsIndependent%independentField, &
                  & FIELD_U_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationSetAndLock(equationsIndependent%independentField, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              ENDDO !componentIdx
              !Default the scaling to the geometric field scaling
              CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
              CALL Field_ScalingTypeSet(equationsIndependent%independentField,geometricScalingType,err,error,*999)
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
              localError="The solution method of "//TRIM(NumberToVString(equationsSet%solutionMethod,"*",err,error))// &
                & " is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ELSE
            !Check the user specified field
            CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)            
            CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)            
            CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,independentNumberOfVariables,err,error,*999)
            CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE], &
              & err,error,*999)
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
              & err,error,*999)
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
              & err,error,*999)
            CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE, &
              & independentNumberOfComponents,err,error,*999)
            CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE, &
              & independentNumberOfComponents,err,error,*999)
            CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)            
            SELECT CASE(solutionMethod)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,1, &
                & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
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
              localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))// &
                & " is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ENDIF
        CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
          IF(equationsIndependent%independentFieldAutoCreated) THEN
            CALL Field_CreateFinish(equationsIndependent%independentField,err,error,*999)
          ENDIF
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
            & " is invalid for a standard, quasistatic or ALE Darcy equation"
          CALL FlagError(localError,err,error,*999)
        END SELECT
      END SELECT
    CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
      !-----------------------------------------------------------------
      !   m a t e r i a l   f i e l d
      !-----------------------------------------------------------------
      NULLIFY(equationsMaterials)
      CALL EquationsSet_MaterialsGet(equationsSet,equationsMaterials,err,error,*999)
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      SELECT CASE(esSpecification(3))
      CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE, EQUATIONS_SET_QUASISTATIC_DARCY_SUBTYPE, EQUATIONS_SET_ALE_DARCY_SUBTYPE, &
        & EQUATIONS_SET_TRANSIENT_DARCY_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
        & EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
        & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, EQUATIONS_SET_MULTI_COMPARTMENT_DARCY_SUBTYPE)
        materialsNumberOfVariables = 1
        SELECT CASE(esSpecification(3))
        CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE,EQUATIONS_SET_QUASISTATIC_DARCY_SUBTYPE,EQUATIONS_SET_ALE_DARCY_SUBTYPE)
          !Porosity + scalar permeability/viscosity
          materialFieldNumberOfComponents = 2
        CASE DEFAULT
          !Porosity + symmetric permeability/viscosity tensor
          materialFieldNumberOfComponents = 7
        END SELECT
        SELECT CASE(equationsSetSetup%actionType)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
          IF(equationsMaterials%materialsFieldAutoCreated) THEN
            !Create the auto created materials field
            CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsMaterials%materialsField,err,error,*999)
            CALL Field_TypeSetAndLock(equationsMaterials%materialsField,FIELD_MATERIAL_TYPE,err,error,*999)
            CALL Field_DependentTypeSetAndLock(equationsMaterials%materialsField,FIELD_INDEPENDENT_TYPE,err,error,*999)
            NULLIFY(geometricDecomposition)
            CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
            CALL Field_DecompositionSetAndLock(equationsMaterials%materialsField,geometricDecomposition,err,error,*999)
            CALL Field_GeometricFieldSetAndLock(equationsMaterials%materialsField,geometricField,err,error,*999)
            CALL Field_NumberOfVariablesSet(equationsMaterials%materialsField,materialsNumberOfVariables,err,error,*999)
            CALL Field_VariableTypesSetAndLock(equationsMaterials%materialsField,[FIELD_U_VARIABLE_TYPE],err,error,*999)
            CALL Field_VariableLabelSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,"Material",err,error,*999)
            CALL Field_DimensionSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DataTypeSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_DP_TYPE,err,error,*999)
            CALL Field_NumberOfComponentsSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
              & materialFieldNumberOfComponents,err,error,*999)
            CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricComponentNumber,err,error,*999)
            !Auto-created / default is node_based_interpolation: that's an expensive default ...
            !Maybe default should be constant; node_based should be requested by the user \todo
            DO componentIdx=1,materialFieldNumberOfComponents
              CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                & componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                & componentIdx,geometricComponentNumber,err,error,*999)
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
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
              & err,error,*999)
            CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE, &
              & materialFieldNumberOfComponents,err,error,*999)
          ENDIF
        CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
          IF( equationsMaterials%materialsFieldAutoCreated ) THEN
            CALL Field_CreateFinish(equationsMaterials%materialsField,err,error,*999)
            !Set the default values for the materials field
            DO componentIdx=1,materialFieldNumberOfComponents
              CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,componentIdx,1.0_DP,err,error,*999)
            ENDDO !componentIdx
          ENDIF
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
            & " is invalid for a standard, quasistatic or ALE Darcy equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
        !Materials field needs two extra variable types
        !The V variable type stores the Darcy coupling coefficients that govern flux between compartments
        !The U1 variable type stores the parameters for the constitutive laws that determine the partial pressure in
        !each compartment
        !For a first attempt at this, it will be assumed that the functional form of this law is the same for each
        !compartment, with only the paramenters varying (default will be three components)
        NULLIFY(equationsSetField)
        CALL EquationsSet_EquationsSetFieldFieldGet(equationsSet,equationsSetField,err,error,*999)
        CALL Field_ParameterSetDataGet(equationsSetField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,equationsSetFieldData, &
          & err,error,*999)
        numberOfCompartments=equationsSetFieldData(2)
        materialsNumberOfVariables = 3
        materialsNumberOfUComponents = 2
        materialsNumberOfVComponents = numberOfCompartments
        materialsNumberOfU1Components = 3
        SELECT CASE(equationsSetSetup%actionType)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
          IF(equationsMaterials%materialsFieldAutoCreated) THEN
            !Create the auto created materials field
            CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsMaterials%materialsField,err,error,*999)
            CALL Field_TypeSetAndLock(equationsMaterials%materialsField,FIELD_MATERIAL_TYPE,err,error,*999)
            CALL Field_DependentTypeSetAndLock(equationsMaterials%materialsField,FIELD_INDEPENDENT_TYPE,err,error,*999)
            NULLIFY(geometricDecomposition)
            CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
            CALL Field_DecompositionSetAndLock(equationsMaterials%materialsField,geometricDecomposition,err,error,*999)
            CALL Field_GeometricFieldSetAndLock(equationsMaterials%materialsField,geometricField,err,error,*999)
            CALL Field_NumberOfVariablesSet(equationsMaterials%materialsField,materialsNumberOfVariables,err,error,*999)
            CALL Field_VariableTypesSetAndLock(equationsMaterials%materialsField,[FIELD_U_VARIABLE_TYPE, &
              & FIELD_V_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE],err,error,*999)
            CALL Field_DimensionSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DataTypeSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            CALL Field_DimensionSetAndLock(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE, &
              & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DataTypeSetAndLock(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            CALL Field_DimensionSetAndLock(equationsMaterials%materialsField,FIELD_U1_VARIABLE_TYPE, &
              & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DataTypeSetAndLock(equationsMaterials%materialsField,FIELD_U1_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            CALL Field_NumberOfComponentsSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
              & materialsNumberOfUComponents,err,error,*999)
            CALL Field_NumberOfComponentsSetAndLock(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE, &
              & materialsNumberOfVComponents,err,error,*999)
            CALL Field_NumberOfComponentsSetAndLock(equationsMaterials%materialsField,FIELD_U1_VARIABLE_TYPE, &
              & materialsNumberOfU1Components,err,error,*999)
            CALL Field_ComponentMeshComponentGet(equationsSet%GEOMETRY%geometricField,FIELD_U_VARIABLE_TYPE, &
              & 1,geometricComponentNumber,err,error,*999)            
            !Auto-created / default is node_based_interpolation: that's an expensive default ...
            !Maybe default should be constant; node_based should be requested by the user \todo
            DO componentIdx = 1, materialsNumberOfUComponents
              CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                & componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                & componentIdx,geometricComponentNumber,err,error,*999)
            ENDDO !componentIdx
            DO componentIdx = 1, materialsNumberOfVComponents
              CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE, &
                & componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE, &
                & componentIdx,geometricComponentNumber,err,error,*999)
            ENDDO !componentIdx
            DO componentIdx = 1, materialsNumberOfU1Components
              CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_U1_VARIABLE_TYPE, &
                & componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_U1_VARIABLE_TYPE, &
                & componentIdx,geometricComponentNumber,err,error,*999)
            ENDDO !componentIdx
            !Default the field scaling to that of the geometric field
            CALL Field_ScalingTypeGet(equationsSet%GEOMETRY%geometricField,geometricScalingType,err,error,*999)
            CALL Field_ScalingTypeSet(equationsMaterials%materialsField,geometricScalingType,err,error,*999)
          ELSE
            !Check the user specified field
            CALL Field_TypeCheck(equationsSetSetup%field,FIELD_MATERIAL_TYPE,err,error,*999)
            CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
            CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,materialsNumberOfVariables,err,error,*999)
            CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE, &
              & FIELD_V_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE],err,error,*999)
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U1_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U1_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,materialsNumberOfUComponents, &
              & err,error,*999)
            CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,materialsNumberOfVComponents, &
              & err,error,*999)
            CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U1_VARIABLE_TYPE,materialsNumberOfU1Components, &
              & err,error,*999)
          ENDIF
        CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
          IF( equationsMaterials%materialsFieldAutoCreated ) THEN
            CALL Field_CreateFinish(equationsMaterials%materialsField,err,error,*999)
            !Set the default values for the materials field
            DO componentIdx=1,materialsNumberOfUComponents
              CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,componentIdx,1.0_DP,err,error,*999)
            ENDDO !componentIdx
            DO componentIdx=1,materialsNumberOfVComponents
              CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,componentIdx,0.0_DP,err,error,*999)
            ENDDO !componentIdx
            DO componentIdx=1,materialsNumberOfU1Components
              CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U1_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,componentIdx,0.0_DP,err,error,*999)
            ENDDO !componentIdx
          ENDIF
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
            & " is invalid for a standard, quasistatic or ALE Darcy equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      END SELECT
    CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
      !-----------------------------------------------------------------
      !   a n a l y t i c   f i e l d
      !-----------------------------------------------------------------
      NULLIFY(equationsAnalytic)
      CALL EquationsSet_AnalyticGet(equationsSet,equationsAnalytic,err,error,*999)
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
      SELECT CASE(esSpecification(3))
      CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE)
        SELECT CASE(equationsSetSetup%actionType)
          !Set start action
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
          CALL EquationsSet_AssertDependentIsFinished(equationsSet,err,error,*999)
          SELECT CASE(equationsSetSetup%analyticFunctionType)
          CASE(EQUATIONS_SET_DARCY_EQUATION_TWO_DIM_1)
            !Set analytic function type
            equationsSet%analytic%analyticFunctionType=EQUATIONS_SET_DARCY_EQUATION_TWO_DIM_1
          CASE(EQUATIONS_SET_DARCY_EQUATION_TWO_DIM_2)
            !Set analytic function type
            equationsSet%analytic%analyticFunctionType=EQUATIONS_SET_DARCY_EQUATION_TWO_DIM_2
          CASE(EQUATIONS_SET_DARCY_EQUATION_TWO_DIM_3)
            !Set analytic function type
            equationsSet%analytic%analyticFunctionType=EQUATIONS_SET_DARCY_EQUATION_TWO_DIM_3
          CASE(EQUATIONS_SET_DARCY_EQUATION_THREE_DIM_1)
            !Set analytic function type
            equationsSet%analytic%analyticFunctionType=EQUATIONS_SET_DARCY_EQUATION_THREE_DIM_1
          CASE(EQUATIONS_SET_DARCY_EQUATION_THREE_DIM_2)
            !Set analytic function type
            equationsSet%analytic%analyticFunctionType=EQUATIONS_SET_DARCY_EQUATION_THREE_DIM_2
          CASE(EQUATIONS_SET_DARCY_EQUATION_THREE_DIM_3)
            !Set analytic function type
            equationsSet%analytic%analyticFunctionType=EQUATIONS_SET_DARCY_EQUATION_THREE_DIM_3
          CASE(EQUATIONS_SET_INCOMP_ELAST_DARCY_ANALYTIC_DARCY)
            !Set analytic function type
            equationsSet%analytic%analyticFunctionType=EQUATIONS_SET_INCOMP_ELAST_DARCY_ANALYTIC_DARCY
          CASE DEFAULT
            localError="The specified analytic function type of "// &
              & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
              & " is invalid for an analytic Darcy problem."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
          IF(equationsSet%analytic%analyticFieldAutoCreated) THEN
            !--- Why finish the dependent field and not the analytic one ???
            !CALL Field_CreateFinish(equationsSet%dependent%dependentField,err,error,*999)
          ENDIF
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
            & " is invalid for an analytic Darcy problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE)
        SELECT CASE(equationsSetSetup%actionType)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
          CALL EquationsSet_AssertDependentIsFinished(equationsSet,err,error,*999)
          !Initialise analytic parameter which stores value of time to zero - need to update this somewhere in a pre_solve routine
          SELECT CASE(equationsSetSetup%analyticFunctionType)
          CASE(EQUATIONS_SET_INCOMP_ELAST_DARCY_ANALYTIC_DARCY)
            !Set analytic function type
            equationsSet%analytic%analyticFunctionType=EQUATIONS_SET_INCOMP_ELAST_DARCY_ANALYTIC_DARCY
          CASE DEFAULT
            localError="The specified analytic function type of "// &
              & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
              & " is invalid for an analytic Darcy problem."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
          IF(equationsSet%analytic%analyticFieldAutoCreated) THEN
            !--- Why finish the dependent field and not the analytic one ???
            !CALL Field_CreateFinish(equationsSet%dependent%dependentField,err,error,*999)
          ENDIF
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
            & " is invalid for an analytic Darcy problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The equation set subtype of "//TRIM(NumberToVString(equationsSet%SPECIFICATION(3),"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a Darcy equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
      !-----------------------------------------------------------------
      !   s o u r c e   t y p e   -   include gravity at some point
      !-----------------------------------------------------------------
      NULLIFY(equationsSource)
      CALL EquationsSet_SourceGet(equationsSet,equationsSource,err,error,*999)
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
      SELECT CASE(esSpecification(3))
      CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE, EQUATIONS_SET_QUASISTATIC_DARCY_SUBTYPE, EQUATIONS_SET_ALE_DARCY_SUBTYPE, &
        & EQUATIONS_SET_TRANSIENT_DARCY_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
        & EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
        & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
        & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
        sourceNumberOfComponents = numberOfDimensions + 1
        SELECT CASE(equationsSetSetup%actionType)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
          IF(equationsSource%sourceFieldAutoCreated) THEN
            CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsSource%sourceField,err,error,*999)
            CALL Field_LabelSet(equationsSource%sourceField,"Source Field",err,error,*999)
            CALL Field_TypeSetAndLock(equationsSource%sourceField,FIELD_GENERAL_TYPE,err,error,*999)
            CALL Field_DependentTypeSetAndLock(equationsSource%sourceField,FIELD_INDEPENDENT_TYPE,err,error,*999)
            NULLIFY(geometricDecomposition)
            CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
            CALL Field_DecompositionSetAndLock(equationsSource%sourceField,geometricDecomposition,err,error,*999)
            CALL Field_GeometricFieldSetAndLock(equationsSource%sourceField,geometricField,err,error,*999)
            CALL Field_NumberOfVariablesSetAndLock(equationsSource%sourceField,1,err,error,*999)
            CALL Field_VariableTypesSetAndLock(equationsSource%sourceField,[FIELD_U_VARIABLE_TYPE],err,error,*999)
            CALL Field_VariableLabelSet(equationsSource%sourceField,FIELD_U_VARIABLE_TYPE,"Source",err,error,*999)
            CALL Field_DimensionSetAndLock(equationsSource%sourceField,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
              & err,error,*999)
            CALL Field_DataTypeSetAndLock(equationsSource%sourceField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            CALL Field_NumberOfComponentsSetAndLock(equationsSource%sourceField,FIELD_U_VARIABLE_TYPE,sourceNumberOfComponents, &
              & err,error,*999)
            !Default the source components to the geometric interpolation setup with nodal interpolation
            IF(esSpecification(3)==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE .OR. &
              & esSpecification(3)==EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE) THEN
              !nodal / mesh based
              DO componentIdx=1,numberOfDimensions !sourceNumberOfComponents
                CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,componentIdx,geometricMeshComponent, &
                  & err,error,*999)
                CALL Field_ComponentInterpolationSet(equationsSource%sourceField,FIELD_U_VARIABLE_TYPE,componentIdx, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentMeshComponentSet(equationsSource%sourceField,FIELD_U_VARIABLE_TYPE,componentIdx, &
                  & geometricMeshComponent,err,error,*999)
              ENDDO !componentIdx
              !Set source component 'numberOfDimensions + 1' according to geometricMeshComponent 'numberOfDimensions'
              CALL Field_ComponentInterpolationSet(equationsSource%sourceField,FIELD_U_VARIABLE_TYPE, &
                & numberOfDimensions+1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              CALL Field_ComponentMeshComponentSet(equationsSource%sourceField,FIELD_U_VARIABLE_TYPE, &
                & numberOfDimensions+1,geometricMeshComponent,err,error,*999)
            ENDIF
            !Default the field scaling to that of the geometric field
            CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
            CALL Field_ScalingTypeSet(equationsSource%sourceField,geometricScalingType,err,error,*999)
          ELSE
            !Check the user specified field
            CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
            CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
            CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,1,err,error,*999)
            CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE],err,error,*999)
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            CALL Field_NumberOfComponentsGet(equationsSet%GEOMETRY%geometricField, FIELD_U_VARIABLE_TYPE,numberOfDimensions, &
              & err,error,*999)
            CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,sourceNumberOfComponents, &
              & err,error,*999)
          ENDIF
        CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
          IF(equationsSource%sourceFieldAutoCreated) THEN
            !Finish creating the source field
            CALL Field_CreateFinish(equationsSource%sourceField,err,error,*999)
            !Set the default values for the source field
            IF(esSpecification(3)==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE .OR. &
              & esSpecification(3)==EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE) THEN
              sourceNumberOfComponents=numberOfDimensions + 1
            ELSE
              sourceNumberOfComponents=0
            ENDIF
            !Now set the source values to 0.0
            DO componentIdx=1,sourceNumberOfComponents
              CALL Field_ComponentValuesInitialise(equationsSource%sourceField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,componentIdx,0.0_DP,err,error,*999)
            ENDDO !componentIdx
          ENDIF
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
            & " is invalid for a standard, quasistatic or ALE Darcy equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      END SELECT
    CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
      !-----------------------------------------------------------------
      !   e q u a t i o n s   t y p e
      !-----------------------------------------------------------------
      SELECT CASE(equationsSet%SPECIFICATION(3))
      CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE,EQUATIONS_SET_MULTI_COMPARTMENT_DARCY_SUBTYPE)
        !-----------------------------------------------------------------
        !   s t a t i c
        !-----------------------------------------------------------------
        SELECT CASE(equationsSetSetup%actionType)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
          CALL EquationsSet_AssertMaterialsIsFinished(equationsSet,err,error,*999)
          NULLIFY(equations)
          CALL Equations_CreateStart(equationsSet,equations,err,error,*999)
          CALL Equations_LinearityTypeSet(equations,EQUATIONS_LINEAR,err,error,*999)
          CALL Equations_TimeDependenceTypeSet(equations,EQUATIONS_STATIC,err,error,*999)
        CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
          CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
          SELECT CASE(solutionMethod)
          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
            SELECT CASE(esSpecification(3))
            CASE(EQUATIONS_SET_MULTI_COMPARTMENT_DARCY_SUBTYPE)
!!!!!THE FOLLOWING IF STATEMENT IS ILLUSTRATIVE ONLY - need to implement the equation set field thing, and make a generalised case statement
              NULLIFY(equationsSetField)
              CALL EquationsSet_EquationsSetFieldGet(equationsSet,equationsSetField,err,error,*999)
              NULLIFY(equations)
              CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
              CALL Equations_CreateFinish(equations,err,error,*999)
              NULLIFY(vectorEquations)
              CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
              !Create the equations mapping.
              NULLIFY(vectorMapping)
              CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_U_VARIABLE_TYPE,vectorMapping,err,error,*999)
              CALL EquationsMappingVector_NumberOfLinearMatricesSet(vectorMapping,1,err,error,*999)
              CALL Field_ParameterSetDataGet(equationsSetField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & equationsSetFieldData,err,error,*999)
              myMatrixIdx = equationsSetFieldData(1)
              numberOfCompartments = equationsSetFieldData(2)
              ALLOCATE(variableTypes(2*numberOfCompartments))
              DO variableIdx=1,numberOfCompartments
                variableTypes(2*variableIdx-1)=FIELD_U_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(variableIdx-1))
                variableTypes(2*variableIdx)=FIELD_DELUDELN_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(variableIdx-1))
              ENDDO !variableIdx
              CALL EquationsMappingVector_LinearMatricesVariableTypesSet(vectorMapping,[variableTypes(2*myMatrixIdx-1)], &
                & err,error,*999)
              CALL EquationsMappingVector_RHSVariableTypeSet(vectorMapping,variableTypes(2*myMatrixIdx),err,error,*999)
              CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
              !Create the equations matrices
              NULLIFY(vectorMatrices)
              CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
              CALL Equations_SparsityTypeGet(equations,sparsityType,err,error,*999)
              SELECT CASE(sparsityType)
              CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                CALL EquationsMatricesVector_LinearStorageTypeSet(vectorMatrices,[MATRIX_BLOCK_STORAGE_TYPE],err,error,*999)
              CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                CALL EquationsMatricesVector_LinearStorageTypeSet(vectorMatrices,[MATRIX_COMPRESSED_ROW_STORAGE_TYPE], &
                  & err,error,*999)
                CALL EquationsMatricesVector_LinearStructureTypeSet(vectorMatrices,[EQUATIONS_MATRIX_FEM_STRUCTURE],err,error,*999)
              CASE DEFAULT
                localError="The equations matrices sparsity type of "// &
                  & TRIM(NumberToVString(sparsityType,"*",err,error))//" is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
              CALL EquationsMatrices_VectorCreateFinish(vectorMatrices,err,error,*999)
            CASE DEFAULT
              !Finish the equations creation
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
              CALL Equations_SparsityTypeGet(equations,sparsityType,err,error,*999)
              SELECT CASE(sparsityType)
              CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                CALL EquationsMatricesVector_LinearStorageTypeSet(vectorMatrices,[MATRIX_BLOCK_STORAGE_TYPE],err,error,*999)
              CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                CALL EquationsMatricesVector_LinearStorageTypeSet(vectorMatrices,[MATRIX_COMPRESSED_ROW_STORAGE_TYPE], &
                  & err,error,*999)
                CALL EquationsMatricesVector_LinearStructureTypeSet(vectorMatrices,[EQUATIONS_MATRIX_FEM_STRUCTURE],err,error,*999)
              CASE DEFAULT
                localError="The equations matrices sparsity type of "//TRIM(NumberToVString(sparsityType,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
              CALL EquationsMatrices_VectorCreateFinish(vectorMatrices,err,error,*999)
            END SELECT
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
            localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
            & " is invalid for a standard Darcy equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_QUASISTATIC_DARCY_SUBTYPE, EQUATIONS_SET_ALE_DARCY_SUBTYPE, &
        & EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE)
        !-----------------------------------------------------------------
        !   q u a s i s t a t i c   and    A L E
        !-----------------------------------------------------------------
        SELECT CASE(equationsSetSetup%actionType)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
          CALL EquationsSet_AssertMaterialsIsFinished(equationsSet,err,error,*999)
          NULLIFY(equations)
          CALL Equations_CreateStart(equationsSet,equations,err,error,*999)
          CALL Equations_LinearityTypeSet(equations,EQUATIONS_LINEAR,err,error,*999)
          CALL Equations_TimeDependenceTypeSet(equations,EQUATIONS_QUASISTATIC,err,error,*999)
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
            CALL EquationsMappingVector_NumberOfLinearMatricesSet(vectorMapping,1,err,error,*999)
            SELECT CASE(esSpecification(3))
            CASE(EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE)
              CALL EquationsMappingVector_LinearMatricesVariableTypesSet(vectorMapping,[FIELD_V_VARIABLE_TYPE],err,error,*999)
              CALL EquationsMappingVector_RHSVariableTypeSet(vectorMapping,FIELD_DELVDELN_VARIABLE_TYPE,err,error,*999)
            CASE DEFAULT
              CALL EquationsMappingVector_LinearMatricesVariableTypesSet(vectorMapping,[FIELD_U_VARIABLE_TYPE],err,error,*999)
              CALL EquationsMappingVector_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE,err,error,*999)
            END SELECT
            CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
            !Create the equations matrices
            NULLIFY(vectorMatrices)
            CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
            CALL Equations_SparsityTypeGet(equations,sparsityType,err,error,*999)
            SELECT CASE(sparsityType)
            CASE(EQUATIONS_MATRICES_FULL_MATRICES)
              CALL EquationsMatricesVector_LinearStorageTypeSet(vectorMatrices,[MATRIX_BLOCK_STORAGE_TYPE],err,error,*999)
            CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
              CALL EquationsMatricesVector_LinearStorageTypeSet(vectorMatrices,[MATRIX_COMPRESSED_ROW_STORAGE_TYPE], &
                & err,error,*999)
              CALL EquationsMatricesVector_LinearStructureTypeSet(vectorMatrices,[EQUATIONS_MATRIX_FEM_STRUCTURE],err,error,*999)
            CASE DEFAULT
              localError="The equations matrices sparsity type of "//TRIM(NumberToVString(sparsityType,"*",err,error))// &
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
            localError="The solution method of "//TRIM(NumberToVString(equationsSet%solutionMethod,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
            & " is invalid for a quasistatic Darcy equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_TRANSIENT_DARCY_SUBTYPE,EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE, &
        & EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
        & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
        !-----------------------------------------------------------------
        !   d y n a m i c
        !-----------------------------------------------------------------
        SELECT CASE(equationsSetSetup%actionType)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
          CALL EquationsSet_AssertMaterialsIsFinished(equationsSet,err,error,*999)
          NULLIFY(equations)
          CALL Equations_CreateStart(equationsSet,equations,err,error,*999)
          CALL Equations_LinearityTypeSet(equations,EQUATIONS_LINEAR,err,error,*999)
          CALL Equations_TimeDependenceTypeSet(equations,EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
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
            IF(esSpecification(3)==EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE .OR. &
              & esSpecification(3)==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE) THEN
              CALL EquationsMappingVector_NumberOfLinearMatricesSet(vectorMapping,0,err,error,*999)
            ENDIF
            CALL EquationsMappingVector_DynamicMatricesSet(vectorMapping,.TRUE.,.TRUE.,err,error,*999)
            SELECT CASE(esSpecification(3))
            CASE(EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
              & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE)
              CALL EquationsMappingVector_DynamicVariableTypeSet(vectorMapping,FIELD_V_VARIABLE_TYPE,err,error,*999)
              CALL EquationsMappingVector_RHSVariableTypeSet(vectorMapping,FIELD_DELVDELN_VARIABLE_TYPE,err,error,*999)
              IF(esSpecification(3)==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE) THEN
                CALL EquationsMappingVector_NumberOfSourcesSet(vectorMapping,1,err,error,*999)
                CALL EquationsMappingVector_SourceVariableTypeSet(vectorMapping,1,FIELD_U_VARIABLE_TYPE,err,error,*999)
              ENDIF
            CASE(EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
              NULLIFY(equationsSetField)
              CALL EquationsSet_EquationsSetFieldFieldGet(equationsSet,equationsSetField,err,error,*999)
              NULLIFY(equationsSetFieldData)
              CALL Field_ParameterSetDataGet(equationsSetField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & equationsSetFieldData,err,error,*999)
              myMatrixIdx = equationsSetFieldData(1)
              numberOfCompartments = equationsSetFieldData(2)
              CALL EquationsMappingVector_NumberOfLinearMatricesSet(vectorMapping,numberOfCompartments-1,err,error,*999)
              ALLOCATE(variableTypes(2*numberOfCompartments+2))
              ALLOCATE(variableUTypes(numberOfCompartments-1))
              DO variableIdx=1,numberOfCompartments+1
                variableTypes(2*variableIdx-1)=FIELD_U_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(variableIdx-1))
                variableTypes(2*variableIdx)=FIELD_DELUDELN_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(variableIdx-1))
              ENDDO !variableIdx
              variableCount=0
              DO variableIdx=2,numberOfCompartments+1
                IF((variableIdx-1)/=myMatrixIdx)THEN
                  variableCount=variableCount+1
                  variableUTypes(variableCount)=variableTypes(2*variableIdx-1)
                ENDIF
              ENDDO
              CALL EquationsMappingVector_DynamicVariableTypeSet(vectorMapping,variableTypes(2*myMatrixIdx+1),err,error,*999)
              CALL EquationsMappingVector_LinearMatricesVariableTypesSet(vectorMapping,variableUTypes,err,error,*999)
              CALL EquationsMappingVector_RHSVariableTypeSet(vectorMapping,variableTypes(2*myMatrixIdx+2),err,error,*999)
              CALL EquationsMappingVector_NumberOfSourcesSet(vectorMapping,1,err,error,*999)
              CALL EquationsMappingVector_SourceVariableTypeSet(vectorMapping,1,FIELD_U_VARIABLE_TYPE,err,error,*999)
            CASE DEFAULT
              CALL EquationsMappingVector_DynamicVariableTypeSet(vectorMapping,FIELD_U_VARIABLE_TYPE,err,error,*999)
              CALL EquationsMappingVector_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE, &
                & err,error,*999)
            END SELECT
            CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
            !Create the equations matrices
            NULLIFY(vectorMatrices)
            CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
            !Set up matrix storage and structure
            CALL Equations_LumpingTypeGet(equations,lumpingType,err,error,*999)
            IF(lumpingType==EQUATIONS_LUMPED_MATRICES) THEN
              !Set up lumping
              CALL EquationsMatricesVector_DynamicLumpingTypeSet(vectorMatrices,[EQUATIONS_MATRIX_UNLUMPED, &
                & EQUATIONS_MATRIX_LUMPED],err,error,*999)
              CALL EquationsMatricesVector_DynamicStorageTypeSet(vectorMatrices,[DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE, &
                & DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE],err,error,*999)
              CALL EquationsMatricesVector_DynamicStructureTypeSet(vectorMatrices,[EQUATIONS_MATRIX_FEM_STRUCTURE, &
                & EQUATIONS_MATRIX_DIAGONAL_STRUCTURE],err,error,*999)
            ELSE
              CALL Equations_SparsityTypeGet(equations,sparsityType,err,error,*999)
              SELECT CASE(sparsityType)
              CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                CALL EquationsMatricesVector_DynamicStorageTypeSet(vectorMatrices,[DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE, &
                  & DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE],err,error,*999)
              CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                CALL EquationsMatricesVector_DynamicStorageTypeSet(vectorMatrices, &
                  & [DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE, &
                  & DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE],err,error,*999)
                CALL EquationsMatricesVector_DynamicStructureTypeSet(vectorMatrices, &
                  & [EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_FEM_STRUCTURE],err,error,*999)
                IF(equationsSet%SPECIFICATION(3)==EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE) THEN
                  ALLOCATE(couplingMatrixStorageType(numberOfCompartments-1))
                  ALLOCATE(couplingMatrixStructureType(numberOfCompartments-1))
                  DO variableIdx=1,numberOfCompartments-1
                    couplingMatrixStorageType(variableIdx)=DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE
                    couplingMatrixStructureType(variableIdx)=EQUATIONS_MATRIX_FEM_STRUCTURE
                  ENDDO
                  CALL EquationsMatricesVector_LinearStorageTypeSet(vectorMatrices,couplingMatrixStorageType, &
                    & err,error,*999)
                  CALL EquationsMatricesVector_LinearStructureTypeSet(vectorMatrices,couplingMatrixStructureType, &
                    & err,error,*999)
                ENDIF
              CASE DEFAULT
                localError="The equations matrices sparsity type of "//TRIM(NumberToVString(sparsityType,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ENDIF
            CALL EquationsMatrices_VectorCreateFinish(vectorMatrices,err,error,*999)
          CASE DEFAULT
            localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
            & " is invalid for a Darcy equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The equation set subtype of "//TRIM(NumberToVString(equationsSet%SPECIFICATION(3),"*",err,error))// &
          & " for a setup of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a Darcy equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
        & " is invalid for a standard, quasistatic, ALE or dynamic Darcy equation."
      CALL FlagError(localError,err,error,*999)      
    END SELECT

    EXITS("Darcy_EquationsSetSetup")
    RETURN
999 ERRORSEXITS("Darcy_EquationsSetSetup",err,error)
    RETURN 1
    
  END SUBROUTINE Darcy_EquationsSetSetup

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices and RHS for a Darcy equation finite element equations set.
  SUBROUTINE Darcy_FiniteElementCalculate(equationsSet,elementNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: analyticFunctionType,componentIdx,columnComponentIdx,columnElementDOFIdx,columnElementParameterIdx, &
      & columnXiIdx,colsVariableType,derivativeIdx,dimensionIdx,equationsSetSubtype,esSpecification(3), &
      & fieldVarTypes(FIELD_NUMBER_OF_VARIABLE_TYPES),gaussPointIdx,globalElementIdx,matrixIdx,meshComponentNumber, &
      & meshComponent1,meshComponent2,myCompartment,numberOfColsComponents,numberOfColsElementParameters,numberOfCompartments, &
      & numberOfDimensions,numberOfDOFs,numberOfGauss,numberOfRowsComponents,numberOfRowsElementParameters, &
      & numberOfVelPressComponents,numberOfXi,rowComponentIdx,rowElementDOFIdx,rowXiIdx,rowElementParameterIdx, &
      & rowsVariableType,scalingType,variableCount,xiIdx
    INTEGER(INTG), POINTER :: equationsSetFieldData(:)
    REAL(DP) :: arg(3),betaParameter,bfact,colsdPhidXi(3),colsPhi,couplingParameter,cParameter,darcyRho0F,dfdJfact,dXdXi(3,3), &
      & dXdY(3,3),dXidX(3,3),dXidY(3,3),dYdXi(3,3),fact,ffact,gaussWeight,gradientLMPressure(3),interComparmentSource, &
      & interCompartmentPermeability1,interCompartmentPermeability2,jacobian,jacobianGaussWeight,Jmat,Jxxi,Jxy,Jyxi,L,lmPressure, &
      & Mfact,p0fact,permeabilityOverViscosity(3,3),permeabilityOverViscosityParameter,pSinkParameter,rowsdPhidXi(3),rowsPhi, &
      & source,sum,viscosityOverPermeability(3,3),x(3)
    REAL(DP), ALLOCATABLE :: pressureCoefficient(:),pressure(:),pressureGradient(:,:)
    LOGICAL :: boundaryElement,stabilized,update,updateCoupling,updateDamping,updateMatrices,updateMatrix,updateRHS, &
      & updateSource,updateStiffness
    TYPE(BasisType), POINTER :: colsBasis,dependentBasis,dependentBasis1,dependentBasis2,geometricBasis,rowsBasis
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable
    TYPE(DecompositionType), POINTER :: dependentDecomposition,geometricDecomposition
    TYPE(DomainType), POINTER :: colsDomain,geometricDomain,rowsDomain
    TYPE(DomainElementsType), POINTER :: colsDomainElements,geometricDomainElements,rowsDomainElements
    TYPE(DomainTopologyType), POINTER :: colsDomainTopology,geometricDomainTopology,rowsDomainTopology
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping
    TYPE(EquationsMappingSourceType), POINTER :: sourceMapping
    TYPE(EquationsMappingSourcesType), POINTER :: sourcesMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatricesSourcesType), POINTER :: sourceVectors
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatrixType), POINTER :: stiffnessMatrix, dampingMatrix
    TYPE(EquationsMatrixPtrType) :: couplingMatrices(99)
    TYPE(EquationsSetAnalyticType), POINTER :: equationsAnalytic
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField,equationsSetField,geometricField,materialsField,sourceField
    TYPE(FieldInterpolatedPointType), POINTER :: elasticityDependentInterpPoint,geometricInterpPoint,materialsInterpPoint, &
      & materialsUInterpPoint,materialsU1InterpPoint,materialsVInterpPoint,referenceGeometricInterpPoint,sourceInterpPoint
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: geometricInterpPointMetrics
    TYPE(FieldInterpolationParametersType), POINTER :: colsInterpParameters,elasticityDependentInterpParameters, &
      & geometricInterpParameters,materialsUInterpParameters,materialsU1InterpParameters,materialsVInterpParameters, &
      & referenceGeometricInterpParameters,rowsInterpParameters,sourceInterpParameters
    TYPE(FieldVariableType), POINTER :: colsVariable,fieldVariable,geometricVariable,rowsVariable
    TYPE(FieldVariablePtrType) :: fieldVariables(99)
    TYPE(MeshElementType), POINTER :: meshElement
    TYPE(QuadratureSchemeType), POINTER :: colsQuadratureScheme,geometricQuadratureScheme,quadratureScheme,quadratureScheme1, &
      & quadratureScheme2,rowsQuadratureScheme
    TYPE(VARYING_STRING) :: localError

    !--- Parameter settings concerning the Finite Element implementation
    stabilized = .TRUE.
    darcy%length = 10.0_DP
    L = darcy%length

    !--- testcase: default
    darcy%testcase = 0
    darcy%analytic = .FALSE.

    ENTERS("Darcy_FiniteElementCalculate",err,error,*999)

    !Parameters settings for coupled elasticity Darcy INRIA model:
    CALL FiniteElasticity_GetDarcyParameters(darcyRho0F,Mfact,bfact,p0fact,err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)

    equationsSetSubtype=esSpecification(3)
    SELECT CASE(equationsSetSubtype)
    CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE, EQUATIONS_SET_QUASISTATIC_DARCY_SUBTYPE, EQUATIONS_SET_ALE_DARCY_SUBTYPE, &
      & EQUATIONS_SET_TRANSIENT_DARCY_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
      & EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
      & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE,EQUATIONS_SET_MULTI_COMPARTMENT_DARCY_SUBTYPE, &
      & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
      !OK
    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NumberToVString(equationsSetSubtype,"*",err,error))// &
        & " is not valid for a Darcy equation type of a fluid mechanics equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
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
    updateRHS=.FALSE.
    IF(ASSOCIATED(rhsMapping)) THEN
      CALL EquationsMatricesVector_RHSVectorGet(vectorMatrices,rhsVector,err,error,*999)
      CALL EquationsMatricesRHS_UpdateVectorGet(rhsVector,updateRHS,err,error,*999)
    ENDIF

    NULLIFY(dynamicMapping)
    NULLIFY(dynamicMatrices)
    NULLIFY(linearMapping)
    NULLIFY(linearMatrices)
    updateStiffness=.FALSE.
    updateDamping=.FALSE.
    updateCoupling=.FALSE.
    NULLIFY(sourcesMapping)
    NULLIFY(sourceVectors)
    NULLIFY(sourceVector)
    updateSource=.FALSE.
    SELECT CASE(equationsSetSubtype)
    CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE, &
      & EQUATIONS_SET_QUASISTATIC_DARCY_SUBTYPE, &
      & EQUATIONS_SET_ALE_DARCY_SUBTYPE, &
      & EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE)
      CALL EquationsMappingVector_LinearMappingGet(vectorMapping,linearMapping,err,error,*999)
      CALL EquationsMatricesVector_LinearMatricesGet(vectorMatrices,linearMatrices,err,error,*999)
      CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,1,stiffnessMatrix,err,error,*999)
      CALL EquationsMatrix_UpdateMatrixGet(stiffnessMatrix,updateStiffness,err,error,*999)
      NULLIFY(colsVariable)
      CALL EquationsMappingLinear_MatrixVariableGet(linearMapping,1,colsVariable,err,error,*999)
    CASE(EQUATIONS_SET_TRANSIENT_DARCY_SUBTYPE, &
      & EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE, &
      & EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
      & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE)
      CALL EquationsMappingVector_DynamicMappingGet(vectorMapping,dynamicMapping,err,error,*999)
      CALL EquationsMatricesVector_DynamicMatricesGet(vectorMatrices,dynamicMatrices,err,error,*999)
      CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,1,stiffnessMatrix,err,error,*999)
      CALL EquationsMatrix_UpdateMatrixGet(stiffnessMatrix,updateStiffness,err,error,*999)
      CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,2,dampingMatrix,err,error,*999)
      CALL EquationsMatrix_UpdateMatrixGet(dampingMatrix,updateDamping,err,error,*999)
      NULLIFY(colsVariable)
      CALL EquationsMappingDynamic_DynamicVariableGet(dynamicMapping,colsVariable,err,error,*999)
      updateSource=.FALSE.
      IF(equationsSetSubtype==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE) THEN
        CALL EquationsMappingVector_SourcesMappingExists(vectorMapping,sourcesMapping,err,error,*999)
        IF(ASSOCIATED(sourcesMapping)) THEN
          CALL EquationsMappingSources_SourceMappingGet(sourcesMapping,1,sourceMapping,err,error,*999)
          CALL EquationsMatricesVector_SourceVectorsGet(vectorMatrices,sourceVectors,err,error,*999)
          CALL EquationsMatricesSources_SourceVectorGet(sourceVectors,1,sourceVector,err,error,*999)
          CALL EquationsMatricesSource_UpdateVectorGet(sourceVector,updateSource,err,error,*999)
        ENDIF
      ENDIF
    CASE(EQUATIONS_SET_MULTI_COMPARTMENT_DARCY_SUBTYPE)
      NULLIFY(linearMapping)
      CALL EquationsMappingVector_LinearMappingGet(vectorMapping,linearMapping,err,error,*999)
      NULLIFY(linearMatrices)
      CALL EquationsMatricesVector_LinearMatricesGet(vectorMatrices,linearMatrices,err,error,*999)
      NULLIFY(stiffnessMatrix)
      CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,1,stiffnessMatrix,err,error,*999)
      CALL EquationsMatrix_UpdateMatrixGet(stiffnessMatrix,updateStiffness,err,error,*999)
      updateDamping=.FALSE.
      updateSource=.FALSE.
      NULLIFY(colsVariable)
      CALL EquationsMappingLinear_LinearMatrixVariableGet(linearMapping,1,colsVariable,err,error,*999)
      NULLIFY(equationsSetField)
      CALL EquationsSet_EquationsSetFieldGet(equationsSet,equationsSetField,err,error,*999)
      NULLIFY(equationsSetFieldData)
      CALL Field_ParameterSetDataGet(equationsSetField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,equationsSetFieldData, &
        & err,error,*999)
      myCompartment = equationsSetFieldData(1)
      numberOfCompartments  = equationsSetFieldData(2) 
    CASE(EQUATIONS_SET_ELASTICITY_MULTI_COMPARTMENT_DARCY_INRIA_SUBTYPE)
      NULLIFY(dynamicMapping)
      CALL EquationsMappingVector_DynamicMappingGet(vectorMapping,dynamicMapping,err,error,*999)
      NULLIFY(dynamicMatrices)
      CALL EquationsMatricesVector_DynamicMatricesGet(vectorMatrices,dynamicMatrices,err,error,*999)
      NULLIFY(stiffnessMatrix)
      CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,1,stiffnessMatrix,err,error,*999)
      CALL EquationsMatrix_UpdateMatrixGet(stiffnessMatrix,updateStiffness,err,error,*999)
      NULLIFY(dampingMatrix)
      CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,2,dampingMatrix,err,error,*999)
      CALL EquationsMatrix_UpdateMatrixGet(dampingMatrix,updateDamping,err,error,*999)
      updateSource=.FALSE.
      NULLIFY(colsVariable)
      CALL EquationsMappingDynamic_DynamicVariableGet(dynamicMapping,colsVariable,err,error,*999)
      NULLIFY(linearMapping)
      CALL EquationsMappingVector_LinearMappingGet(vectorMapping,linearMapping,err,error,*999)
      NULLIFY(linearMatrices)
      CALL EquationsMatricesVector_LinearMatricesGet(vectorMatrices,linearMatrices,err,error,*999)
      variableCount=0
      DO matrixIdx=1,numberOfCompartments
        IF(matrixIdx/=myCompartment)THEN
          variableCount=variableCount+1
          NULLIFY(couplingMatrices(variableCount)%ptr)
          CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,variableCount,couplingMatrices(variableCount)%ptr, &
            & err,error,*999)
          CALL EquationsMatrix_UpdateMatrixGet(couplingMatrices(variableCount)%ptr,updateMatrix,err,error,*999)
          updateCoupling=updateCoupling.OR.updateMatrix
          NULLIFY(fieldVariables(variableCount)%ptr)
          CALL EquationsMappingLinear_LinearMatrixVariableGet(linearMapping,variableCount,fieldVariables(variableCount)%ptr, &
            & err,error,*999)
          CALL EquationsMappingLinear_LinearMatrixVariableTypeGet(linearMapping,variableCount,fieldVarTypes(variableCount), &
            & err,error,*999)          
        ENDIF
      ENDDO !matrixIdx
    CASE(EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
      CALL EquationsMappingVector_DynamicMappingGet(vectorMapping,dynamicMapping,err,error,*999)
      CALL EquationsMatricesVector_DynamicMatricesGet(vectorMatrices,dynamicMatrices,err,error,*999)
      CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,1,stiffnessMatrix,err,error,*999)
      CALL EquationsMatrix_UpdateMatrixGet(stiffnessMatrix,updateStiffness,err,error,*999)
      CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,2,dampingMatrix,err,error,*999)
      CALL EquationsMatrix_UpdateMatrixGet(dampingMatrix,updateDamping,err,error,*999)
      NULLIFY(colsVariable)
      CALL EquationsMappingDynamic_DynamicVariableGet(dynamicMapping,colsVariable,err,error,*999)
      !These linear matrices are actually only required if we are coupling the momentum terms too
      !If it is just a mass coupling, then all of the additional terms are placed in the RHS of the mass-increase equation      
      CALL EquationsMappingVector_LinearMappingGet(vectorMapping,linearMapping,err,error,*999)
      CALL EquationsMatricesVector_LinearMatricesGet(vectorMatrices,linearMatrices,err,error,*999)
      variableCount=0
      DO matrixIdx=1,numberOfCompartments
        IF(matrixIdx/=myCompartment)THEN
          variableCount=variableCount+1
          NULLIFY(couplingMatrices(variableCount)%ptr)
          CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,variableCount,couplingMatrices(variableCount)%ptr, &
            & err,error,*999)
          CALL EquationsMatrix_UpdateMatrixGet(couplingMatrices(variableCount)%ptr,updateMatrix,err,error,*999)
          updateCoupling=updateCoupling.OR.updateMatrix
          NULLIFY(fieldVariables(variableCount)%ptr)
          CALL EquationsMappingLinear_LinearMatrixVariableGet(linearMapping,variableCount,fieldVariables(variableCount)%ptr, &
            & err,error,*999)
          CALL EquationsMappingLinear_LinearMatrixVariableTypeGet(linearMapping,variableCount,fieldVarTypes(variableCount), &
            & err,error,*999)          
        ENDIF
      ENDDO !matrixIdx
      CALL EquationsMappingVector_SourcesMappingExists(vectorMapping,sourcesMapping,err,error,*999)
      IF(ASSOCIATED(sourcesMapping)) THEN
        CALL EquationsMappingSources_SourceMappingGet(sourcesMapping,1,sourceMapping,err,error,*999)
        CALL EquationsMatricesVector_SourceVectorsGet(vectorMatrices,sourceVectors,err,error,*999)
        CALL EquationsMatricesSources_SourceVectorGet(sourceVectors,1,sourceVector,err,error,*999)
        CALL EquationsMatricesSource_UpdateVectorGet(sourceVector,updateSource,err,error,*999)
      ENDIF
      ALLOCATE(pressureCoefficient(numberOfCompartments),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate pressure coefficients.",err,error,*999)
      ALLOCATE(pressure(numberOfCompartments),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate pressure.",err,error,*999)
      ALLOCATE(pressureGradient(3,numberOfCompartments),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate pressure gradient.",err,error,*999)
      pressure = 0.0_DP
      pressureGradient = 0.0_DP
      pressureCoefficient(1)=0.25_DP
      pressureCoefficient(2)=0.25_DP
      pressureCoefficient(3)=0.25_DP
      pressureCoefficient(4)=0.25_DP
    END SELECT

    updateMatrices=(updateStiffness.OR.updateDamping)
    update=(updateMatrices.OR.updateCoupling.OR.updateSource.OR.updateRHS)

    IF(update) THEN

      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      NULLIFY(dependentField)
      CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
      NULLIFY(materialsField)
      CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)
      NULLIFY(sourceField)
      IF(equationsSetSubtype==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE .OR. &
        & equationsSetSubtype==EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE) THEN      
        CALL EquationsSet_SourceFieldGet(equationsSet,sourceField,err,error,*999)
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
      CALL DomainElements_ElementBoundaryElementGet(geometricDomainElements,elementNumber,boundaryElement,err,error,*999)
      NULLIFY(geometricBasis)
      CALL DomainElements_ElementBasisGet(geometricDomainElements,elementNumber,geometricBasis,err,error,*999)
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
      CALL DomainElements_ElementBasisGet(colsDomainElements,elementNumber,colsBasis,err,error,*999)

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

      NULLIFY(equationsInterpolation)
      CALL Equations_InterpolationGet(equations,equationsInterpolation,err,error,*999)

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
      NULLIFY(referenceGeometricInterpParameters)
      NULLIFY(referenceGeometricInterpPoint)
      IF(equationsSetSubtype==EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE.OR. &
        & equationsSetSubtype==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE .OR. &
        & equationsSetSubtype==EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE) THEN
        CALL FieldVariable_InterpolationParameterInitialise(geometricVariable,referenceGeometricInterpParameters,err,error,*999)
        CALL Field_InterpolatedPointInitialise(referenceGeometricInterpParameters,referenceGeometricInterpPoint,err,error,*999)
        CALL Field_InterpolationParametersElementGet(FIELD_INITIAL_VALUES_SET_TYPE,elementNumber, &
          & referenceGeometricInterpParameters,err,error,*999)
      ENDIF

      NULLIFY(elasticityDependentInterpParameters)
      NULLIFY(elasticityDependentInterpPoint)
      IF(equationsSetSubtype==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE .OR. &
        & equationsSetSubtype==EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE) THEN
        CALL EquationsInterpolation_DependentParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
          & elasticityDependentInterpParameters,err,error,*999)
        CALL EquationsInterpolation_DependentPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,elasticityDependentInterpPoint, &
          & err,error,*999)
        CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,elasticityDependentInterpParameters, &
          & err,error,*999)
      ENDIF

      NULLIFY(materialsUInterpParameters)
      CALL EquationsInterpolation_MaterialsParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & materialsUInterpParameters,err,error,*999)
      NULLIFY(materialsUInterpPoint)
      CALL EquationsInterpolation_MaterialsPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,materialsUInterpPoint, &
        & err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,materialsUInterpParameters,err,error,*999)
      IF(equationsSetSubtype==EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE) THEN
        NULLIFY(materialsVInterpParameters)
        CALL EquationsInterpolation_MaterialsParametersGet(equationsInterpolation,FIELD_V_VARIABLE_TYPE, &
          & materialsVInterpParameters,err,error,*999)
        NULLIFY(materialsVInterpPoint)
        CALL EquationsInterpolation_MaterialsPointGet(equationsInterpolation,FIELD_V_VARIABLE_TYPE,materialsVInterpPoint, &
          & err,error,*999)
        CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,materialsVInterpParameters,err,error,*999)
        NULLIFY(materialsU1InterpParameters)
        CALL EquationsInterpolation_MaterialsParametersGet(equationsInterpolation,FIELD_U1_VARIABLE_TYPE, &
          & materialsU1InterpParameters,err,error,*999)
        NULLIFY(materialsU1InterpPoint)
        CALL EquationsInterpolation_MaterialsPointGet(equationsInterpolation,FIELD_U1_VARIABLE_TYPE,materialsU1InterpPoint, &
          & err,error,*999)
        CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,materialsU1InterpParameters,err,error,*999)
      ENDIF

      NULLIFY(sourceInterpParameters)
      NULLIFY(sourceInterpPoint)
      IF(equationsSetSubtype==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE .OR. &
        & equationsSetSubtype==EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE) THEN
        CALL EquationsInterpolation_MaterialsParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
          & sourceInterpParameters,err,error,*999)
        CALL EquationsInterpolation_MaterialsPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,sourceInterpPoint, &
          & err,error,*999)
        CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,sourceInterpParameters,err,error,*999)
      ENDIF

      SELECT CASE(equationsSetSubtype)
      CASE(EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
        & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)        
        numberOfVelPressComponents = numberOfColsComponents - 1  !last component: mass increase
      CASE DEFAULT
        numberOfVelPressComponents = numberOfColsComponents
      END SELECT

      !---------------------------------------------------------------------------------------------------------
      !Invoke penalty term to enforce impermeable BC. 
      !This sshould only be executed if THIS element lies on the surface
      !(within the routine we check whether the element nodes have actually been set impermeable)
      SELECT CASE(equationsSetSubtype)
      CASE(EQUATIONS_SET_TRANSIENT_DARCY_SUBTYPE,EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE, &
        & EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
        & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE)
        IF(boundaryElement) CALL Darcy_ImpermeableBCViaPenalty(equationsSet,elementNumber,err,error,*999)
      END SELECT
      !---------------------------------------------------------------------------------------------------------

      !--- Loop over gauss points
      !    Given that also materials field is interpolated, ensure sufficient number of Gauss points !!!
      DO gaussPointIdx=1,numberOfGauss

        CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,geometricInterpPoint, &
          & err,error,*999)
        CALL Field_InterpolatedPointMetricsCalculate(numberOfXi,geometricInterpPointMetrics,err,error,*999)
        !--- Retrieve local map dXdXi
        DO componentIdx=1,numberOfDimensions
          DO xiIdx=1,numberOfXi
            derivativeIdx=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xiIdx) !2,4,7
            dXdXi(componentIdx,xiIdx)=geometricInterpPoint%values(componentIdx,derivativeIdx) !dx/dxi
          ENDDO !xiIdx
        ENDDO !componentIdx
        CALL Invert(dXdXi,dXidX,Jxxi,err,error,*999) !dy/dxi -> dxi/dy

        IF(equationsSetSubtype==EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE.OR. &
          & equationsSetSubtype==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE .OR. &
          & equationsSetSubtype==EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE) THEN
          !------------------------------------------------------------------------------
          !--- begin: Compute the Jacobian of the mapping

          !--- Interpolation of Reference Geometry
          CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx, &
            & referenceGeometricInterpPoint,err,error,*999)
          !--- Retrieve local map dYdXi
          DO componentIdx=1,numberOfDimensions
            DO xiIdx=1,numberOfXi
              derivativeIdx=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xiIdx) !2,4,7
              dYdXi(componentIdx,xiIdx)=referenceGeometricInterpPoint%values(componentIdx,derivativeIdx) !dy/dxi (y = referential)
            ENDDO !xiIdx
          ENDDO !componentIdx

          !--- Compute deformation gradient tensor dXdY and its Jacobian Jxy
          CALL Invert(dYdXi,dXidY,Jyxi,err,error,*999) !dy/dxi -> dxi/dy
          CALL MatrixProduct(dXdXi,dXidY,dXdY,err,error,*999) !dx/dxi * dxi/dy = dx/dy (deformation gradient tensor, F)
          CALL Determinant(dXdY,Jxy,err,error,*999)

          IF( ABS(Jxy) < 1.0E-10_DP ) THEN
            localError="Darcy_FiniteElementCalculate: Jacobian Jxy is smaller than 1.0E-10_DP."
            CALL FlagError(localError,err,error,*999)
          END IF

          !ffact = f(Jxy) of the INRIA model, dfdJfact is not relevant here
          CALL FiniteElasticity_EvaluateChapelleFunction(Jxy,ffact,dfdJfact,err,error,*999)

          !--- end: Compute the Jacobian of the mapping
          !------------------------------------------------------------------------------
        END IF

        !--- Material Settings ---!
        !*** If material is variable, need to account for this in deriving the variational statement ***!        

        !--- Interpolate materials field
        !Get the Darcy permeability
        CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx, &
          & materialsUInterpPoint,err,error,*999)
        IF(equationsSetSubtype==EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE) THEN
          !Get the intercompartmental permeabilities
          CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,materialsVInterpPoint, &
            & err,error,*999)
          !Get the material parameters for the constitutive law for each Darcy compartment (for determining the partial pressures)
          CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,materialsU1InterpPoint, &
            & err,error,*999)
        ENDIF

        IF(equationsSetSubtype==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE .OR. &
          & equationsSetSubtype==EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE) THEN
          CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,sourceInterpPoint, &
            & err,error,*999)
        END IF

        SELECT CASE(equationsSetSubtype)
        CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE,EQUATIONS_SET_QUASISTATIC_DARCY_SUBTYPE,EQUATIONS_SET_ALE_DARCY_SUBTYPE)
          !scalar permeability/viscosity
          permeabilityOverViscosity=0.0_DP
          permeabilityOverViscosity(1,1) = materialsUInterpPoint%values(2,NO_PART_DERIV)
          permeabilityOverViscosity(2,2) = materialsUInterpPoint%values(2,NO_PART_DERIV)
          permeabilityOverViscosity(3,3) = materialsUInterpPoint%values(2,NO_PART_DERIV)
          !Multiply by porosity
          permeabilityOverViscosity=permeabilityOverViscosity*materialsUInterpPoint%values(1,NO_PART_DERIV)
        CASE DEFAULT
          !symmetric permeability/viscosity tensor
          permeabilityOverViscosity(1,1) = materialsUInterpPoint%values(2,NO_PART_DERIV)
          permeabilityOverViscosity(1,2) = materialsUInterpPoint%values(3,NO_PART_DERIV)
          permeabilityOverViscosity(1,3) = materialsUInterpPoint%values(4,NO_PART_DERIV)
          permeabilityOverViscosity(2,2) = materialsUInterpPoint%values(5,NO_PART_DERIV)
          permeabilityOverViscosity(2,3) = materialsUInterpPoint%values(6,NO_PART_DERIV)
          permeabilityOverViscosity(3,3) = materialsUInterpPoint%values(7,NO_PART_DERIV)

          permeabilityOverViscosity(2,1) = permeabilityOverViscosity(1,2)
          permeabilityOverViscosity(3,1) = permeabilityOverViscosity(1,3)
          permeabilityOverViscosity(3,2) = permeabilityOverViscosity(2,3)
        END SELECT

        IF(diagnostics3) THEN
          IF(idebug2) THEN
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"materialsUInterpPoint%values(1,NO_PART_DERIV) = ", &
              & materialsUInterpPoint%values(1,NO_PART_DERIV),err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"materialsUInterpPoint%values(2,NO_PART_DERIV) = ", &
              & materialsUInterpPoint%values(2,NO_PART_DERIV),err,error,*999)
            CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"",err,error,*999)
            idebug2 = .FALSE.
          ENDIF
        ENDIF

        CALL Determinant(permeabilityOverViscosity,Jmat,err,error,*999)
        IF(Jmat>ZERO_TOLERANCE) THEN
          CALL Invert(permeabilityOverViscosity,viscosityOverPermeability,Jmat,err,error,*999)
        ELSE
          viscosityOverPermeability = 0.0_DP
          DO componentIdx=1,3
            viscosityOverPermeability(componentIdx,componentIdx) = 1.0e10_DP
          ENDDO !componentIdx
        ENDIF

        !Two parameters that are used only for TESTCASE==3: VenousCompartment problem: Exclude this, too specific ???
        betaParameter   = - darcy%PERM_OVER_VIS * (2.0_DP * PI / darcy%length) * (2.0_DP * PI / darcy%length)
        pSinkParameter = darcy%P_SINK

        IF(equationsSetSubtype==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE) THEN
          CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx, &
            & elasticityDependentInterpPoint,err,error,*999)
          !Mind the sign !!!
          !The minus sign derives from the convention of using "+ P * Jznu * AZU(i,j)"
          ! in the constitutive law in FiniteElasticity_GaussCauchyTensor
          lmPressure = -elasticityDependentInterpPoint%values(4,NO_PART_DERIV)
          DO xiIdx=1,numberOfXi
            derivativeIdx=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xiIdx) !2,4,7
            !gradient wrt. element coordinates xi
            gradientLMPressure(xiIdx) = -elasticityDependentInterpPoint%values(4,derivativeIdx)
          ENDDO !xiIdx
        ENDIF

        !For multi-compartment model - determine pressure from partial derivative of constitutive law
        IF(equationsSetSubtype==EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE) THEN
          CALL FlagWarning("NEED CONSTITUTIVE LAWS HERE!!!! THE FOLLOWING IS PLACEHOLDER ONLY!",err,error,*999)
          !BEGIN PLACEHOLDER
          CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx, &
            & elasticityDependentInterpPoint,err,error,*999)
          !Mind the sign !!!
          lmPressure = -elasticityDependentInterpPoint%values(4,NO_PART_DERIV)
          DO xiIdx=1,dependentBasis%numberOfXi
            derivativeIdx=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xiIdx) !2,4,7
            !gradient wrt. element coordinates xi
            gradientLMPressure(xiIdx) = -elasticityDependentInterpPoint%values(4,derivativeIdx)
          ENDDO !xiIdx
          !loop over compartments to determine the pressure in each one - this could be quite inefficient, as it will be
          !calculated several times over

          !unless calculate the pressures in a pre-solve and store them in extra components/variables of the dependent field
          !these pressures should really be known immediately after the finite elasticity solve and not determined here
          !END PLACEHOLDER
          !The following pressure_coeff matrix is just for testing purposes and ultimately will be replaced with functions
          !and materials field parameters (for present, sum of coefficients should be 1).

          DO matrixIdx=1,numberOfCompartments
            pressure(matrixIdx) =  pressureCoefficient(matrixIdx)*lmPressure
            DO xiIdx=1,numberOfXi
              pressureGradient(xiIdx,matrixIdx) = pressureCoefficient(matrixIdx)*gradientLMPressure(xiIdx)
            ENDDO !xiIdx
          ENDDO !matrixIdx
        ENDIF

        CALL FieldInterpolatedPointsMetrics_JacobianGet(geometricInterpPointMetrics,jacobian,err,error,*999)
        CALL BasisQuadratureScheme_GaussWeightGet(geometricQuadratureScheme,gaussPointIdx,gaussWeight,err,error,*999)
        jacobianGaussWeight=jacobian*gaussWeight

        !Loop over element rows
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
          CALL Basis_NumberOfElementParametersGet(rowsBasis,numberOfRowsElementParameters,err,error,*999)
          !Loop over element rows
          DO rowElementParameterIdx=1,numberOfRowsElementParameters
            rowElementDOFIdx=rowElementDOFIdx+1
            CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowsQuadratureScheme,rowElementParameterIdx,NO_PART_DERIV, &
              & gaussPointIdx,rowsPhi,err,error,*999)
            DO xiIdx=1,numberOfXi
              CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowsQuadratureScheme,rowElementParameterIdx, &
                & PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xiIdx),gaussPointIdx,rowsdPhidXi(xiIdx),err,error,*999)
            ENDDO !xiIdx
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
                    & gaussPointIdx,colsPhi,err,error,*999)
                  DO xiIdx=1,numberOfXi
                    CALL BasisQuadratureScheme_GaussBasisFunctionGet(colsQuadratureScheme,columnElementParameterIdx, &
                      & PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xiIdx),gaussPointIdx,colsdPhidXi(xiIdx),err,error,*999)
                  ENDDO !xiIdx
                  SELECT CASE(equationsSetSubtype)
                  CASE(EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE)
                    !----------------------------------------------------------------------------------------------------
                    !  i n c o m p r e s s i b l e   e l a s t i c i t y   d r i v e n   D a r c y   :   M A T R I C E S
                    !----------------------------------------------------------------------------------------------------
                    IF(updateStiffness) THEN
                      !velocity test function, velocity trial function
                      IF(rowComponentIdx==columnComponentIdx.AND.columnComponentIdx<numberOfColsComponents) THEN

                        sum = viscosityOverPermeability( rowComponentIdx, columnComponentIdx ) * rowsPhi * colsPhi
                        !MIND: double check the matrix index order: (rowComponentIdx, columnComponentIdx)
                        !or (columnComponentIdx, rowComponentIdx) !within this conditional: rowComponentIdx==columnComponentIdx
                        !anyway

                        stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) = &
                          & stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) + &
                          & sum * jacobianGaussWeight
                        !--------------------------------------------------------------------------------------------------------
                        !mass-increase test function, velocity trial function
                      ELSE IF(rowComponentIdx==numberOfRowsComponents.AND. &
                        & columnComponentIdx<numberOfColsComponents) THEN

                        sum = 0.0_DP
                        DO columnXiIdx=1,numberOfXi
                          sum = sum + rowsPhi * colsdPhidXi(columnXiIdx) * dXidX(columnXiIdx,columnComponentIdx)
                        ENDDO !columnXiIdx

                        stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) = &
                          & stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) + &
                          & sum * jacobianGaussWeight
                      ENDIF
                    ENDIF
                    IF(updateDamping) THEN
                      !MASS-INCREASE test function, mass-increase trial function
                      IF(rowComponentIdx==columnComponentIdx.AND.columnComponentIdx==numberOfColsComponents) THEN

                        !To integrate the mass-increase term in the reference configuration, we divide by Jxy.
                        sum = rowsPhi * colsPhi / (Jxy * darcyRho0F)

                        dampingMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) = &
                          & dampingMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) + &
                          & sum * jacobianGaussWeight
                      END IF

                      ! !Try out adding the inertia term ...
                      ! IF(rowComponentIdx==columnComponentIdx.AND.rowComponentIdx<fieldVariable%numberOfComponents) THEN
                      !   rowsPhi=quadratureScheme1%gaussBasisFunctions(rowElementParameterIdx,NO_PART_DERIV,gaussPointIdx)
                      !   colsPhi=quadratureScheme2%gaussBasisFunctions(columnElementParameterIdx,NO_PART_DERIV,gaussPointIdx)

                      !   sum = 0.0_DP

                      !   sum = rowsPhi*colsPhi*darcyRho0F

                      !   dampingMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) = &
                      !     & dampingMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) + &
                      !     & sum * jacobianGaussWeight
                      ! END IF

                    END IF

                  CASE(EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
                    ! matrices for multi-compartment poroelastic equations
                    IF(updateStiffness) THEN

                      !velocity test function, velocity trial function
                      IF(rowComponentIdx==columnComponentIdx.AND.columnComponentIdx<numberOfColsComponents) THEN

                        sum = viscosityOverPermeability( rowComponentIdx, columnComponentIdx ) * rowsPhi * colsPhi
                        !MIND: double check the matrix index order: (rowComponentIdx, columnComponentIdx) or
                        !(columnComponentIdx, rowComponentIdx) within this conditional: rowComponentIdx==columnComponentIdx anyway

                        stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) = &
                          & stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) + &
                          & sum * jacobianGaussWeight

                        !--------------------------------------------------------------------------------------------------------
                        !mass-increase test function, velocity trial function
                      ELSE IF(rowComponentIdx==numberOfRowsComponents.AND.columnComponentIdx<numberOfColsComponents) THEN

                        sum = 0.0_DP
                        DO columnXiIdx=1,dependentBasis2%numberOfXi
                          sum = sum + rowsPhi * colsdPhidXi(columnXiIdx) * dXidX(columnXiIdx,columnComponentIdx)
                        ENDDO !columnXiIdx

                        stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) = &
                          & stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) + &
                          & sum * jacobianGaussWeight
                      ENDIF
                    ENDIF !update stiffness

                    IF(updateDamping) THEN
                      !MASS-INCREASE test function, mass-increase trial function
                      IF(rowComponentIdx==columnComponentIdx.AND.columnComponentIdx==numberOfColsComponents) THEN

                        !To integrate the mass-increase term in the reference configuration, we divide by Jxy.
                        sum = rowsPhi * colsPhi / (Jxy * darcyRho0F)
                        
                        dampingMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) = &
                          & dampingMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) + &
                          & sum * jacobianGaussWeight
                      ENDIF

!!Try out adding the inertia term ...
                      !IF(rowComponentIdx==columnComponentIdx.AND.rowComponentIdx<fieldVariable%numberOfComponents) THEN
                      !  rowsPhi=quadratureScheme1%gaussBasisFunctions(rowElementParameterIdx,NO_PART_DERIV,gaussPointIdx)
                      !  colsPhi=quadratureScheme2%gaussBasisFunctions(columnElementParameterIdx,NO_PART_DERIV,gaussPointIdx)
                      !
                      !  sum = 0.0_DP
                      !
                      !  sum = rowsPhi*colsPhi*darcyRho0F
                      !
                      !  dampingMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) = &
                      !    dampingMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) + &
                      !    & sum * jacobianGaussWeight
                      !ENDIF
                      
                    ENDIF !update damping


                    !=================================================================================
                    !    d e f a u l t   :   M A T R I C E S
                  CASE DEFAULT

                    IF(updateStiffness) THEN
                      !---------------------------------------------------------------------------------------------------
                      !velocity test function, velocity trial function
                      IF(rowComponentIdx==columnComponentIdx.AND.columnComponentIdx<numberOfVelPressComponents) THEN
                        
                        sum=viscosityOverPermeability( rowComponentIdx, columnComponentIdx ) * rowsPhi * colsPhi
                        !MIND: double check the matrix index order: (rowComponentIdx, columnComponentIdx) or
                        !(columnComponentIdx, rowComponentIdx) within this conditional: rowComponentIdx==columnComponentIdx
                        !anyway
                        
                        IF( stabilized ) THEN
                          sum=sum-0.5_DP*viscosityOverPermeability(rowComponentIdx,columnComponentIdx)*rowsPhi*colsPhi
                        END IF
                        
                        stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) = &
                          & stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) + &
                          & sum * jacobianGaussWeight
                        
                        !--------------------------------------------------------------------------------------------------------
                        !velocity test function, pressure trial function
                      ELSE IF(rowComponentIdx<numberOfVelPressComponents.AND.columnComponentIdx==numberOfVelPressComponents) THEN
                        
                        sum = 0.0_DP
                        DO rowXiIdx=1,numberOfXi
                          sum = sum - rowsdPhidXi(rowXiIdx)*colsPhi*dXidX(rowXiIdx,rowComponentIdx)
                        ENDDO !rowXiIdx
                        
                        IF( stabilized ) THEN
                          DO columnXiIdx=1,numberOfXi
                            sum = sum - 0.5_DP * rowsPhi * colsdPhidXi(columnXiIdx) * dXidX(columnXiIdx,rowComponentIdx)
                          ENDDO !columnXiIdx
                        END IF

                        stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) = &
                          & stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) + &
                          & sum * jacobianGaussWeight

                        !---------------------------------------------------------------------------------------------------------
                        !pressure test function, velocity trial function
                      ELSE IF(rowComponentIdx==numberOfVelPressComponents.AND.columnComponentIdx<numberOfVelPressComponents) THEN
                        
                        sum = 0.0_DP
                        DO columnXiIdx=1,dependentBasis2%numberOfXi
                          sum = sum + rowsPhi * colsdPhidXi(columnXiIdx) * dXidX(columnXiIdx,columnComponentIdx)
                        ENDDO !columnXiIdx
                        
                        IF( stabilized ) THEN
                          DO rowXiIdx=1,dependentBasis1%numberOfXi
                            sum = sum + 0.5_DP * rowsdPhidXi(rowXiIdx) * colsPhi * dXidX(rowXiIdx,columnComponentIdx)
                          ENDDO !rowXiIdx
                        END IF
                        
                        stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) = &
                          & stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) + &
                          & sum * jacobianGaussWeight
                        
                        !---------------------------------------------------------------------------------------------------------
                        !pressure test function, pressure trial function
                      ELSE IF(rowComponentIdx==columnComponentIdx.AND.columnComponentIdx==numberOfVelPressComponents) THEN
                        
                        sum = 0.0_DP

                        IF( stabilized ) THEN
                          DO dimensionIdx =1,numberOfDimensions !number space dimension equiv. number of xi
                            DO rowXiIdx=1,numberOfXi
                              DO columnXiIdx=1,numberOfXi
                                sum = sum + 0.5_DP * permeabilityOverViscosity( dimensionIdx, dimensionIdx ) * &
                                  & rowsdPhidXi(rowXiIdx) * colsdPhidXi(columnXiIdx) * dXidX(rowXiIdx,dimensionIdx) *  &
                                  & dXidX(columnXiIdx,dimensionIdx)
                              ENDDO !columnXiIdx
                            ENDDO !rowXiIdx
                          ENDDO !dimensionIdx
                        END IF

                        IF( darcy%testcase == 3 ) THEN
                          !This forms part of the pressure-dependent source term,
                          !thus it enters the LHS
                          
                          sum = sum + betaParameter * rowsPhi * colsPhi
                        END IF
                        
                        stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) = &
                          & stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) + &
                          & sum * jacobianGaussWeight
                        
                        !--------------------------------------------------------------------------------------------------------
                        !For the INRIA model, and: mass-increase test function, pressure trial function
                      ELSE IF(equationsSetSubtype==EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE.AND. &
                        & rowComponentIdx==numberOfRowsComponents.AND. &
                        & columnComponentIdx==numberOfVelPressComponents) THEN


                        sum = -rowsPhi * colsPhi / (Mfact * ffact)

                        stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) = &
                          & stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) + &
                          & sum * jacobianGaussWeight
                        
                        !--------------------------------------------------------------------------------------------------------
                        !For the INRIA model, and: mass-increase test function, mass-increase trial function
                      ELSE IF(equationsSetSubtype==EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE.AND. &
                        & rowComponentIdx==columnComponentIdx.AND.columnComponentIdx==numberOfColsComponents) THEN


                        sum = rowsPhi * colsPhi / darcyRho0F
                        
                        stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) = &
                          & stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) + &
                          & sum * jacobianGaussWeight
                        
                        !-------------------------------------------------------------------------------------------------------
                      ELSE

                        stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) = 0.0_DP
                        
                      ENDIF
                      
                    ENDIF !update stiffness
                    
                    !=======================================================================================================
                    ! dampingMatrix
                    
                    IF(equationsSetSubtype==EQUATIONS_SET_TRANSIENT_DARCY_SUBTYPE) THEN
                      IF(updateDamping) THEN
                        IF(rowComponentIdx==columnComponentIdx.AND.rowComponentIdx<numberOfVelPressComponents) THEN
                          
                          sum = rowsPhi*colsPhi*darcyRho0F
                          
                          dampingMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) = &
                            & dampingMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) + &
                            & sum * jacobianGaussWeight
                        ENDIF
                      ENDIF !update damping
                    ELSE IF(equationsSetSubtype==EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE) THEN
                      IF(updateDamping) THEN
                        !pressure test function, mass-increase trial function
                        IF(rowComponentIdx==numberOfVelPressComponents.AND.columnComponentIdx==numberOfColsComponents) THEN
                          
                          sum = rowsPhi * colsPhi / (Jxy * darcyRho0F)
                          
                          dampingMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) = &
                            & dampingMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) + &
                            & sum * jacobianGaussWeight
                        ENDIF
                      ENDIF !update damping
                    ENDIF

                  END SELECT
                  !   e n d   s e l e c t   equationsSetSubtype
                  !=================================================================================
                  
                ENDDO !columnElementParameterIdx
              ENDDO !columnComponentIdx
            ENDIF !udpate matrices
            !===================================================================================================================
            !rhsVector
            IF(updateRHS) THEN
              
              SELECT CASE(equationsSetSubtype)
                !==========================================================================================
                !  i n c o m p r e s s i b l e   e l a s t i c i t y   d r i v e n   D a r c y   :   R H S
              CASE(EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE)
                
                !---------------------------------------------------------------------------------------------------------
                !velocity test function
                IF( rowComponentIdx<numberOfRowsComponents ) THEN
                  
                  sum = 0.0_DP
                  
                  !Term arising from the pressure / Lagrange Multiplier of elasticity (given):
                  DO rowXiIdx=1,dependentBasis1%numberOfXi
                    sum = sum - rowsPhi * gradientLMPressure(rowXiIdx) * dXidX(rowXiIdx,rowComponentIdx)
                  ENDDO !rowXiIdx
                  
                  rhsVector%elementVector%vector(rowElementDOFIdx) = &
                    & rhsVector%elementVector%vector(rowElementDOFIdx) + &
                    & sum * jacobianGaussWeight
                  
                  !-------------------------------------------------------------------------------------------------------
                  !mass-increase test function
                ELSE IF( rowComponentIdx==numberOfRowsComponents ) THEN
                  
                  ! + possible source AND SINK TERMS
                  source = 0.0_DP
                  
                  sum = rowsPhi * source

                  rhsVector%elementVector%vector(rowElementDOFIdx) = &
                    & rhsVector%elementVector%vector(rowElementDOFIdx) + &
                    & sum * jacobianGaussWeight

                ELSE
                  
                  rhsVector%elementVector%vector(rowElementDOFIdx) = 0.0_DP
                  
                ENDIF
                !------------------------------------------------------------------------------------------------------
              CASE(EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
                !------------------------------------------------------------------------------------------------------
                !velocity test function
                IF( rowComponentIdx<fieldVariable%numberOfComponents ) THEN

                  sum = 0.0_DP
                  !Term arising from the pressure / Lagrange Multiplier of elasticity (given):
                  !TO DO- need to read different grad p depending on the compartment of interest
                  DO rowXiIdx=1,numberOfXi
                    !sum = sum - rowsPhi * gradientLMPressure(rowXiIdx) * dXidX(rowXiIdx,rowComponentIdx)
                    !this is the pressure gradient for the appropriate compartment
                    sum = sum - rowsPhi * pressureGradient(rowXiIdx,myCompartment) * dXidX(rowXiIdx,rowComponentIdx)
                  ENDDO !rowXiIdx

                  rhsVector%elementVector%vector(rowElementDOFIdx) = &
                    & rhsVector%elementVector%vector(rowElementDOFIdx) + &
                    & sum * jacobianGaussWeight

                  !-----------------------------------------------------------------------------------------------------
                  !mass-increase test function
                ELSE IF( rowComponentIdx==fieldVariable%numberOfComponents ) THEN

                  ! n o   s o u r c e
                  !source terms need to be converted to use source field & vector
                  source = 0.0_DP                      

                  !Add in the source/sink terms due to the pressure difference between compartments
                  DO matrixIdx=1,numberOfCompartments
                    IF(matrixIdx/=myCompartment) THEN
                      !Interpolate the coupling material parameter from the V variable type of the materials field
                      interCompartmentPermeability1=materialsVInterpPoint%values(myCompartment,NO_PART_DERIV)
                      interCompartmentPermeability2=materialsVInterpPoint%values(matrixIdx,NO_PART_DERIV)
                      !Source term is coefficient*(p(myCompartment) - p(matrixIdx))
                      interComparmentSource=-interCompartmentPermeability1*pressure(myCompartment) + &
                        & interCompartmentPermeability2*pressure(matrixIdx)
                    ENDIF
                  ENDDO !matrixIdx                      

                  sum = rowsPhi * source + rowsPhi * interComparmentSource

                  rhsVector%elementVector%vector(rowElementDOFIdx) = &
                    & rhsVector%elementVector%vector(rowElementDOFIdx) + sum * jacobianGaussWeight

                ELSE

                  rhsVector%elementVector%vector(rowElementDOFIdx) = 0.0_DP

                ENDIF
                !=================================================================================
                !    d e f a u l t   :   R H S
              CASE DEFAULT
                !-------------------------------------------------------------------------------------------------
                !velocity test function
                IF( rowComponentIdx<numberOfVelPressComponents ) THEN

                  sum = 0.0_DP

                  rhsVector%elementVector%vector(rowElementDOFIdx) = &
                    & rhsVector%elementVector%vector(rowElementDOFIdx) + sum * jacobianGaussWeight

                  !-------------------------------------------------------------------------------------------------
                  !pressure test function
                ELSE IF( rowComponentIdx==numberOfVelPressComponents ) THEN

                  ! n o   s o u r c e
                  source = 0.0_DP

                  sum = rowsPhi * source

                  rhsVector%elementVector%vector(rowElementDOFIdx) = &
                    & rhsVector%elementVector%vector(rowElementDOFIdx) + sum * jacobianGaussWeight

                  !-------------------------------------------------------------------------------------------------------------
                  !For the INRIA model, and: mass-increase test function
                ELSE IF(equationsSetSubtype==EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE.AND. &
                  & rowComponentIdx==numberOfRowsComponents) THEN

                  sum = -rowsPhi * bfact * (1.0_DP - Jxy)

                  sum = sum - rowsPhi * p0fact / (Mfact * ffact)

                  rhsVector%elementVector%vector(rowElementDOFIdx) = &
                    & rhsVector%elementVector%vector(rowElementDOFIdx) + sum * jacobianGaussWeight

                ELSE

                  rhsVector%elementVector%vector(rowElementDOFIdx) = 0.0_DP

                END IF
                !-------------------------------------------------------------------------------------------------------------
              END SELECT
              !   e n d   s e l e c t   equationsSetSubtype
              !=================================================================================

            ENDIF !update RHS

            IF(updateSource) THEN
              IF(equationsSetSubtype==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE .OR. &
                & equationsSetSubtype==EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE) THEN

                cParameter=sourceInterpPoint%values(rowComponentIdx, NO_PART_DERIV)

                !IF(ABS(cParameter)>1.0E-08) WRITE(*,*)'cParameter = ',cParameter

                sum = rowsPhi * cParameter
                sourceVector%elementVector%vector(rowElementDOFIdx) = &
                  & sourceVector%elementVector%vector(rowElementDOFIdx) + sum * jacobianGaussWeight
              ENDIF
            ENDIF !update source
          ENDDO !rowElementParameterIdx
        ENDDO !rowComponentIdx

        IF(equationsSetSubtype==EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE) THEN
          !Calculate the momentum coupling matrices

          !Loop over element rows
          rowElementDOFIdx=0
          DO rowComponentIdx=1,numberOfRowsComponents !field_variable is the variable associated with the equations set under consideration
            NULLIFY(rowsDomain)
            CALL FieldVariable_ComponentDomainGet(rowsVariable,rowComponentIdx,rowsDomain,err,error,*999)
            NULLIFY(rowsDomainTopology)
            CALL Domain_DomainTopologyGet(rowsDomain,rowsDomainTopology,err,error,*999)
            NULLIFY(rowsDomainElements)
            CALL DomainTopology_DomainElementsGet(rowsDomainTopology,rowsDomainElements,err,error,*999)
            NULLIFY(rowsBasis)
            CALL DomainElements_ElementBasisGet(rowsDomainElements,elementNumber,rowsBasis,err,error,*999)
            CALL Basis_QuadratureSchemeGet(rowsBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,rowsQuadratureScheme,err,error,*999)
            CALL Basis_NumberOfElementParametersGet(rowsBasis,numberOfRowsElementParameters,err,error,*999)

            DO rowElementParameterIdx=1,numberOfRowsElementParameters
              rowElementDOFIdx=rowElementDOFIdx+1
              CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowsQuadratureScheme,rowElementParameterIdx,NO_PART_DERIV, &
                & gaussPointIdx,rowsPhi,err,error,*999)

              variableCount=0
              DO matrixIdx = 1,numberOfCompartments
                IF(matrixIdx/=myCompartment)THEN
                  variableCount=variableCount+1

                  !need to test for the case where matrixIdx==mycompartment
                  !the coupling terms then needs to be added into the stiffness matrix
                  CALL EquationsMatrix_UpdateMatrixGet(couplingMatrices(variableCount)%ptr,updateMatrix,err,error,*999)
                  IF(updateMatrix) THEN

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
                      CALL Basis_QuadratureSchemeGet(colsBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,colsQuadratureScheme, &
                        & err,error,*999)
                      CALL Basis_NumberOfElementParametersGet(colsBasis,numberOfColsElementParameters,err,error,*999)
                      DO columnElementParameterIdx=1,numberOfColsElementParameters
                        columnElementDOFIdx=columnElementDOFIdx+1
                        CALL BasisQuadratureScheme_GaussBasisFunctionGet(colsQuadratureScheme,columnElementParameterIdx, &
                          & NO_PART_DERIV,gaussPointIdx,colsPhi,err,error,*999)

                        !---------------------------------------------------------------------------------------------------
                        !concentration test function, concentration trial function
                        !For now, this is only a dummy implementation - this still has to be properly set up.
                        !IF(rowComponentIdx==columnComponentIdx.AND.columnComponentIdx<numberOfVelPressComponents) THEN ! don't need this for diffusion equation

                        !                             sum = 0.0_DP

                        !Get the coupling coefficients
                        couplingParameter=materialsVInterpPoint%values(matrixIdx,NO_PART_DERIV)

                        !                              sum = sum + couplingParameter * rowsPhi * PGN

                        couplingMatrices(variableCount)%ptr%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) = &
                          & couplingMatrices(variableCount)%ptr%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) + &
                          & couplingParameter * rowsPhi * colsPhi * jacobianGaussWeight
                        !                           ENDIF

                      ENDDO !columnElementParameterIdx
                    ENDDO !columnComponentIdx
                  ENDIF
                ENDIF
              ENDDO !matrixIdx
            ENDDO !rowElementParameterIdx
          ENDDO !rowComponentIdx
        ENDIF ! multi component

        !---------------------------------------------------------------------------------------------------------------
        ! RIGHT HAND SIDE FOR ANALYTIC SOLUTION
        !---------------------------------------------------------------------------------------------------------------

        NULLIFY(equationsAnalytic)
        CALL EquationsSet_AnalyticExists(equationsSet,equationsAnalytic,err,error,*999)
        IF(ASSOCIATED(equationsAnalytic)) THEN
          CALL EquationsSet_AnalyticFunctionTypeGet(equationsSet,analyticFunctionType,err,error,*999)
          IF(analyticFunctionType==EQUATIONS_SET_DARCY_EQUATION_TWO_DIM_1.OR. &
            & analyticFunctionType==EQUATIONS_SET_DARCY_EQUATION_TWO_DIM_2.OR. &
            & analyticFunctionType==EQUATIONS_SET_DARCY_EQUATION_TWO_DIM_3.OR. &
            & analyticFunctionType==EQUATIONS_SET_DARCY_EQUATION_THREE_DIM_1.OR. &
            & analyticFunctionType==EQUATIONS_SET_DARCY_EQUATION_THREE_DIM_2.OR. &
            & analyticFunctionType==EQUATIONS_SET_DARCY_EQUATION_THREE_DIM_3) THEN

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
              CALL Basis_NumberOfElementParametersGet(rowsBasis,numberOfRowsElementParameters,err,error,*999)
              !Loop over element rows
              DO rowElementParameterIdx=1,numberOfRowsElementParameters
                rowElementDOFIdx=rowElementDOFIdx+1
                CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowsQuadratureScheme,rowElementParameterIdx,NO_PART_DERIV, &
                  & gaussPointIdx,rowsPhi,err,error,*999)
                !note rowComponentIdx value derivative
                sum=0.0_DP

                x(1) = geometricInterpPoint%values(1,1)
                x(2) = geometricInterpPoint%values(2,1)
                IF(numberOfDimensions==3) x(3) = geometricInterpPoint%values(3,1)
                IF(analyticFunctionType==EQUATIONS_SET_DARCY_EQUATION_TWO_DIM_1) THEN
                  sum=0.0_DP
                ELSE IF(analyticFunctionType==EQUATIONS_SET_DARCY_EQUATION_TWO_DIM_2) THEN
                  IF(rowComponentIdx==3) THEN
                    fact   = permeabilityOverViscosityParameter / L
                    arg(1) = x(1) / L
                    arg(2) = x(2) / L
                    source = -2.0_DP / L * fact * EXP( arg(1) ) * EXP( arg(2) )
                    sum = rowsPhi * source
                  ELSE
                    sum = 0.0_DP
                  ENDIF
                ELSE IF(analyticFunctionType==EQUATIONS_SET_DARCY_EQUATION_TWO_DIM_3) THEN
                  IF(rowComponentIdx==3) THEN
                    fact   = 2.0_DP * PI * permeabilityOverViscosityParameter / L
                    arg(1) = 2.0_DP * PI * x(1) / L
                    arg(2) = 2.0_DP * PI * x(2) / L
                    source = +2.0_DP * (2.0_DP * PI / L) * fact * SIN( arg(1) ) * SIN( arg(2) )
                    sum = rowsPhi * source
                  ELSE
                    sum = 0.0_DP
                  ENDIF
                ELSE IF(analyticFunctionType==EQUATIONS_SET_DARCY_EQUATION_THREE_DIM_1) THEN
                  sum=0.0_DP
                ELSE IF(analyticFunctionType==EQUATIONS_SET_DARCY_EQUATION_THREE_DIM_2) THEN
                  IF(rowComponentIdx==4) THEN
                    fact   = permeabilityOverViscosityParameter / L
                    arg(1) = x(1) / L
                    arg(2) = x(2) / L
                    arg(3) = x(3) / L
                    source = -3.0_DP / L * fact * EXP( arg(1) ) * EXP( arg(2) ) * EXP( arg(3) )
                    sum = rowsPhi * source
                  ELSE
                    sum = 0.0_DP
                  ENDIF
                ELSE IF(analyticFunctionType==EQUATIONS_SET_DARCY_EQUATION_THREE_DIM_3) THEN
                  IF(rowComponentIdx==4) THEN
                    fact   = 2.0_DP * PI * permeabilityOverViscosityParameter / L
                    arg(1) = 2.0_DP * PI * x(1) / L
                    arg(2) = 2.0_DP * PI * x(2) / L
                    arg(3) = 2.0_DP * PI * x(3) / L
                    source = +3.0_DP * ( 2.0_DP * PI / L ) * fact * SIN( arg(1) ) * SIN( arg(2) ) * SIN( arg(3) )
                    sum = rowsPhi * source
                  ELSE
                    sum = 0.0_DP
                  END IF
                ENDIF

                !Calculate RHS VECTOR
                rhsVector%elementVector%vector(rowElementDOFIdx)= &
                  & rhsVector%elementVector%vector(rowElementDOFIdx)+sum*jacobianGaussWeight
              ENDDO !rowElementParameterIdx
            ENDDO !rowComponentIdx
          ELSE
            rhsVector%elementVector%vector(rowElementDOFIdx)=0.0_DP
          ENDIF
        ENDIF !analytic

        ! end: RIGHT HAND SIDE FOR ANALYTIC SOLUTION
        !-------------------------------------------------------------------------------------------------------------

        ! !===================================================================================================================
        ! !couplingMatrices
        ! SELECT CASE(equationsSetSubtype)
        ! CASE(EQUATIONS_SET_MULTI_COMPARTMENT_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_MULTI_COMPARTMENT_DARCY_INRIA_SUBTYPE)

        !   !Create fieldVariables type, couplingMatrices type

        !   !Loop over element rows
        !   rowElementDOFIdx=0
        !   DO rowComponentIdx=1,fieldVariable%numberOfComponents

        !     meshComponent1 = fieldVariable%COMPONENTS(rowComponentIdx)%meshComponentNumber
        !     dependentBasis1 => dependentField%decomposition%DOMAIN(meshComponent1)%ptr% &
        !       & TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS
        !     quadratureScheme1 => dependentBasis1%QUADRATURE% &
        !       & quadratureSchemeMap(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
        !     jacobianGaussWeight = equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%jacobian * &
        !       & quadratureScheme1%gaussWeights(gaussPointIdx)

        !     DO rowElementParameterIdx=1,dependentBasis1%numberOfElementParameters
        !       rowElementDOFIdx=rowElementDOFIdx+1

        !       DO matrixIdx=1,numberOfCompartments

        !         IF(couplingMatrices(matrixIdx)%ptr%updateMatrix) THEN

        !           !Loop over element columns
        !           columnElementDOFIdx=0
        !           !                       DO columnComponentIdx=1,fieldVariable%numberOfComponents
        !           DO columnComponentIdx=1,fieldVariables(matrixIdx)%ptr%numberOfComponents

        !             meshComponent2 = fieldVariable%COMPONENTS(columnComponentIdx)%meshComponentNumber
        !             dependentBasis2 => dependentField%decomposition%DOMAIN(meshComponent2)%ptr% &
        !               & TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS
        !             !--- We cannot use two different quadrature schemes here !!!
        !             quadratureScheme2 => dependentBasis2%QUADRATURE% &
        !               & quadratureSchemeMap(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
        !             !jacobianGaussWeight = equations%interpolation%geometricInterpPointMetrics%jacobian * &
        !             !  & quadratureScheme2%gaussWeights(gaussPointIdx)

        !             DO columnElementParameterIdx=1,dependentBasis2%numberOfElementParameters
        !               columnElementDOFIdx=columnElementDOFIdx+1

        !               !-------------------------------------------------------------------------------------------------------------
        !               !velocity test function, velocity trial function
        !               !For now, this is only a dummy implementation - this still has to be properly set up.
        !               IF(rowComponentIdx==columnComponentIdx.AND.columnComponentIdx<numberOfVelPressComponents) THEN

        !                 sum = 0.0_DP

        !                 rowsPhi=quadratureScheme1%gaussBasisFunctions(rowElementParameterIdx,NO_PART_DERIV,gaussPointIdx)
        !                 colsPhi=quadratureScheme2%gaussBasisFunctions(columnElementParameterIdx,NO_PART_DERIV,gaussPointIdx)

        !                 sum = sum + viscosityOverPermeability( rowComponentIdx, columnComponentIdx ) * rowsPhi * PGN

        !                 couplingMatrices(matrixIdx)%ptr%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) = &
        !                   & couplingMatrices(matrixIdx)%ptr%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) + sum * jacobianGaussWeight
        !               ENDIF

        !             ENDDO !columnElementParameterIdx
        !           ENDDO !columnComponentIdx
        !         ENDIF
        !       ENDDO !matrixIdx
        !     ENDDO !rowElementParameterIdx
        !   ENDDO !rowComponentIdx
        ! CASE DEFAULT
        !   !Do nothing
        ! END SELECT

      ENDDO !gaussPointIdx

      IF(updateRHS) THEN
        ! Integrate pressure over faces, and add to RHS vector
        CALL Darcy_FiniteElementFaceIntegrate(equationsSet,elementNumber,fieldVariable,err,error,*999)
      ENDIF !update RHS

      !Scale factor adjustment
      CALL Field_ScalingTypeGet(dependentField,scalingType,err,error,*999)
      IF(scalingType/=FIELD_NO_SCALING) THEN
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
          NULLIFY(rowsDomain)
          CALL FieldVariable_ComponentDomainGet(rowsVariable,rowComponentIdx,rowsDomain,err,error,*999)
          NULLIFY(rowsDomainTopology)
          CALL Domain_DomainTopologyGet(rowsDomain,rowsDomainTopology,err,error,*999)
          NULLIFY(rowsDomainElements)
          CALL DomainTopology_DomainElementsGet(rowsDomainTopology,rowsDomainElements,err,error,*999)
          NULLIFY(rowsBasis)
          CALL DomainElements_ElementBasisGet(rowsDomainElements,elementNumber,rowsBasis,err,error,*999)
          CALL Basis_NumberOfElementParametersGet(rowsBasis,numberOfRowsElementParameters,err,error,*999)
          DO rowElementParameterIdx=1,numberOfRowsElementParameters
            rowElementDOFIdx=rowElementDOFIdx+1                    
            IF(updateMatrices) THEN
              columnElementDOFIdx=0
              !Loop over element columns
              DO columnComponentIdx=1,numberOfColsComponents
                NULLIFY(colsDomain)
                CALL FieldVariable_ComponentDomainGet(colsVariable,columnComponentIdx,colsDomain,err,error,*999)
                NULLIFY(colsDomainTopology)
                CALL Domain_DomainTopologyGet(colsDomain,colsDomainTopology,err,error,*999)
                NULLIFY(colsDomainElements)
                CALL DomainTopology_DomainElementsGet(colsDomainTopology,colsDomainElements,err,error,*999)
                NULLIFY(colsBasis)
                CALL DomainElements_ElementBasisGet(colsDomainElements,elementNumber,colsBasis,err,error,*999)
                CALL Basis_NumberOfElementParametersGet(colsBasis,numberOfColsElementParameters,err,error,*999)
                DO columnElementParameterIdx=1,numberOfColsElementParameters
                  columnElementDOFIdx=columnElementDOFIdx+1
                  IF(updateStiffness)THEN
                    stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                      & stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)* &
                      & rowsInterpParameters%scaleFactors(rowElementParameterIdx,rowComponentIdx)* &
                      & colsInterpParameters%scaleFactors(columnElementParameterIdx,columnComponentIdx)
                  ENDIF !update stiffness
                  IF(updateDamping)THEN
                    dampingMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                      & dampingMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)* &
                      & rowsInterpParameters%scaleFactors(rowElementParameterIdx,rowComponentIdx)* &
                      & colsInterpParameters%scaleFactors(columnElementParameterIdx,columnComponentIdx)
                  ENDIF !update damping
                ENDDO !columnElementParameterIdx
              ENDDO !columnComponentIdx
            ENDIF !update matrices
            IF(updateRHS) THEN
              rhsVector%elementVector%vector(rowElementDOFIdx)= &
                & rhsVector%elementVector%vector(rowElementDOFIdx)* &
                & rowsInterpParameters%scaleFactors(rowElementParameterIdx,rowComponentIdx)
            ENDIF !update RHS
            IF(updateSource) THEN
              sourceVector%elementVector%vector(rowElementDOFIdx)= &
                & sourceVector%elementVector%vector(rowElementDOFIdx)* &
                & rowsInterpParameters%scaleFactors(rowElementParameterIdx,rowComponentIdx)
            ENDIF !update source
          ENDDO !rowElementParameterIdx
        ENDDO !rowComponentIdx
      ENDIF

      ! RESTORE ALL POINTERS CALL PARAMATER_SET_FIELD_DATA_RESTORE

    ENDIF !update

    EXITS("Darcy_FiniteElementCalculate")
    RETURN
999 ERRORSEXITS("Darcy_FiniteElementCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE Darcy_FiniteElementCalculate

  !
  !================================================================================================================================
  !

  !>Calculates the face integration term of the finite element formulation for Darcy's equation,
  !>required for pressure boundary conditions.
  SUBROUTINE Darcy_FiniteElementFaceIntegrate(equationsSet,elementNumber,dependentVariable,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<The equations set to calculate the RHS term for
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to calculat the RHS term for
    TYPE(FieldVariableType), POINTER :: dependentVariable
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: componentIdx,elementBaseDOFIdx,elementDOFIdx,elementNodeIdx,esSpecification(3),faceIdx, &
      & faceNodeDerivativeIdx,faceNodeIdx,faceNumber,faceParameterIdx,gaussPointIdx,meshComponentNumber,nodeDerivativeIdx, &
      & normalComponentIdx,numberOfDependentComponents,numberOfElementParameters,numberOfFaceNodes,numberOfGauss, &
      & numberOfLocalFaces,numberOfNodeDerivatives,parameterIdx,variableType
    REAL(DP) :: facePhi,gaussWeight,jacobian,normalProjection,pressureGauss
    LOGICAL :: boundaryFace,calculateFaces,updateRHS
    TYPE(BasisType), POINTER :: dependentBasis,faceBasis
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionElementsType), POINTER :: decompositionElements
    TYPE(DecompositionFacesType), POINTER :: decompositionFaces
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainElementsType), POINTER :: domainElements
    TYPE(DomainFacesType), POINTER :: domainFaces
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField
    TYPE(FieldInterpolationParametersType), POINTER :: dependentInterpParameters,geometricInterpParameters
    TYPE(FieldInterpolatedPointType), POINTER :: dependentInterpPoint,geometricInterpPoint
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: geometricInterpPointMetrics
    TYPE(QuadratureSchemeType), POINTER :: faceQuadratureScheme

    ENTERS("Darcy_FiniteElementFaceIntegrate",err,error,*999)

    IF(.NOT.ASSOCIATED(dependentVariable)) CALL FlagError("Dependent variable is not associated.",err,error,*999)
    
    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMapping)
    CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
    NULLIFY(rhsMapping)
    CALL EquationsMappingVector_RHSMappingExists(vectorMapping,rhsMapping,err,error,*999)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
    NULLIFY(rhsVector)
    updateRHS=.FALSE.
    IF(ASSOCIATED(rhsMapping)) THEN
      CALL EquationsMatricesVector_RHSVectorGet(vectorMatrices,rhsVector,err,error,*999)
      CALL EquationsMatricesRHS_UpdateVectorGet(rhsVector,updateRHS,err,error,*999)
    ENDIF

    IF(updateRHS) THEN

      CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)

      SELECT CASE(esSpecification(3))
      CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE, &
        & EQUATIONS_SET_QUASISTATIC_DARCY_SUBTYPE, &
        & EQUATIONS_SET_ALE_DARCY_SUBTYPE, &
        & EQUATIONS_SET_TRANSIENT_DARCY_SUBTYPE)

        !Get the mesh decomposition and basis for this element
        NULLIFY(dependentField)
        CALL FieldVariable_FieldGet(dependentVariable,dependentField,err,error,*999)
        NULLIFY(decomposition)
        CALL Field_DecompositionGet(dependentField,decomposition,err,error,*999)
        !Only add RHS terms if the face geometric parameters are calculated
        CALL Decomposition_CalculateFacesGet(decomposition,calculateFaces,err,error,*999)
        IF(calculateFaces) THEN
          !These RHS terms are associated with the equations for the three velocity components,
          !rather than the pressure term
          CALL FieldVariable_VariableTypeGet(dependentVariable,variableType,err,error,*999)
          CALL FieldVariable_NumberOfComponentsGet(dependentVariable,numberOfDependentComponents,err,error,*999)
          NULLIFY(domain)
          CALL FieldVariable_ComponentDomainGet(dependentVariable,1,domain,err,error,*999)
          NULLIFY(domainTopology)
          CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
          NULLIFY(domainElements)
          CALL DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*999)
          NULLIFY(domainFaces)
          CALL DomainTopology_DomainFacesGet(domainTopology,domainFaces,err,error,*999)
          NULLIFY(dependentBasis)
          CALL DomainElements_ElementBasisGet(domainElements,elementNumber,dependentBasis,err,error,*999)
          CALL Basis_NumberOfElementParametersGet(dependentBasis,numberOfElementParameters,err,error,*999)
          NULLIFY(decompositionTopology)
          CALL Decomposition_DecompositionTopologyGet(decomposition,decompositionTopology,err,error,*999)
          NULLIFY(decompositionElements)
          CALL DecompositionTopology_DecompositionElementsGet(decompositionTopology,decompositionElements,err,error,*999)
          NULLIFY(decompositionFaces)
          CALL DecompositionTopology_DecompositionFacesGet(decompositionTopology,decompositionFaces,err,error,*999)
          !Get interpolation parameters and point for Darcy pressure
          NULLIFY(equationsInterpolation)
          CALL Equations_InterpolationGet(equations,equationsInterpolation,err,error,*999)
          NULLIFY(geometricInterpParameters)
          CALL EquationsInterpolation_GeometricParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
            & geometricInterpParameters,err,error,*999)
          NULLIFY(geometricInterpPoint)
          CALL EquationsInterpolation_GeometricPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
            & geometricInterpPoint,err,error,*999)
          NULLIFY(geometricInterpPointMetrics)
          CALL EquationsInterpolation_GeometricPointMetricsGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
            & geometricInterpPointMetrics,err,error,*999)
          CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,geometricInterpParameters, &
            & err,error,*999)
          NULLIFY(dependentInterpParameters)
          CALL EquationsInterpolation_DependentParametersGet(equationsInterpolation,variableType,dependentInterpParameters, &
            & err,error,*999)
          NULLIFY(dependentInterpPoint)
          CALL EquationsInterpolation_DependentPointGet(equationsInterpolation,variableType,dependentInterpPoint,err,error,*999)
          CALL Basis_NumberOfLocalFacesGet(dependentBasis,numberOfLocalFaces,err,error,*999)
          DO faceIdx=1,numberOfLocalFaces
            !Get the face normal and quadrature information
            CALL DecompositionElements_ElementFaceNumberGet(decompositionElements,faceIdx,faceNumber,err,error,*999)
            !This speeds things up but is also important, as non-boundary faces have an XI_DIRECTION that might
            !correspond to the other element.
            CALL DecompositionFaces_FaceBoundaryFaceGet(decompositionFaces,faceIdx,boundaryFace,err,error,*999)
            IF(.NOT.boundaryFace) CYCLE
            CALL Field_InterpolationParametersFaceGet(FIELD_VALUES_SET_TYPE,faceNumber,dependentInterpParameters, &
              & err,error,*999)
            CALL DecompositionFaces_FaceXiNormalDirectionGet(decompositionFaces,faceIdx,normalComponentIdx,err,error,*999)
            normalComponentIdx=ABS(normalComponentIdx)
            NULLIFY(faceBasis)
            CALL DomainFaces_FaceBasisGet(domainFaces,faceNumber,faceBasis,err,error,*999)            
            NULLIFY(faceQuadratureScheme)
            CALL Basis_QuadratureSchemeGet(faceBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,faceQuadratureScheme,err,error,*999)
            CALL BasisQuadrature_NumberOfGaussGet(faceQuadratureScheme,numberOfGauss,err,error,*999)
            CALL Basis_NumberOfNodesGet(faceBasis,numberOfFaceNodes,err,error,*999)
            DO gaussPointIdx=1,numberOfGauss
              !Get interpolated Darcy pressure
              CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,dependentInterpPoint, &
                & err,error,*999)
              pressureGauss=dependentInterpPoint%values(4,NO_PART_DERIV) !(component,derivative)
              !Use the geometric field to find the face normal and the Jacobian for the face integral
              CALL Field_InterpolateLocalFaceGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,faceIdx,gaussPointIdx, &
                & geometricInterpPoint,err,error,*999)
              !Calculate the metric tensors and Jacobian
              CALL Field_InterpolatedPointMetricsCalculate(COORDINATE_JACOBIAN_VOLUME_TYPE,geometricInterpPointMetrics, &
                & err,error,*999)
              CALL FieldInterpolatedPointsMetrics_JacobianGet(geometricInterpPointMetrics,jacobian,err,error,*999)
              CALL BasisQuadratureScheme_GaussWeightGet(faceQuadratureScheme,gaussPointIdx,gaussWeight,err,error,*999)
              DO componentIdx=1,numberOfDependentComponents-1
                normalProjection=DOT_PRODUCT(geometricInterpPointMetrics%gu(normalComponentIdx,:), &
                  & geometricInterpPointMetrics%dXdXi(componentIdx,:))
                IF(normalComponentIdx<0) normalProjection=-normalProjection
                IF(ABS(normalProjection)<ZERO_TOLERANCE) CYCLE
                !Work out the first index of the rhs vector for this element - 1
                elementBaseDofIdx=numberOfElementParameters*(componentIdx-1)
                DO faceNodeIdx=1,numberOfFaceNodes
                  CALL Basis_FaceNodeNumberGet(dependentBasis,faceNodeIdx,faceIdx,elementNodeIdx,err,error,*999)
                  CALL Basis_FaceNodeNumberOfDerivativesGet(dependentBasis,faceNodeIdx,faceIdx,numberOfNodeDerivatives, &
                    & err,error,*999)
                  DO faceNodeDerivativeIdx=1,numberOfNodeDerivatives
                    CALL Basis_FaceNodeDerivativeNumberGet(dependentBasis,faceNodeDerivativeIdx,faceNodeIdx,faceIdx, &
                      & nodeDerivativeIdx,err,error,*999)
                    CALL Basis_ElementParameterGet(dependentBasis,nodeDerivativeIdx,elementNodeIdx,parameterIdx,err,error,*999)
                    CALL Basis_ElementParameterGet(faceBasis,faceNodeDerivativeIdx,faceNodeIdx,faceParameterIdx,err,error,*999)
                    CALL BasisQuadratureScheme_GaussBasisFunctionGet(faceQuadratureScheme,faceParameterIdx,NO_PART_DERIV, &
                      & gaussPointIdx,facePhi,err,error,*999)
                    elementDofIdx=elementBaseDofIdx+parameterIdx
                    rhsVector%elementVector%vector(elementDofIdx) = rhsVector%elementVector%vector(elementDofIdx) - &
                      & gaussWeight*pressureGauss*normalProjection*facePhi*jacobian
                  ENDDO !nodeDerivativeIdx
                ENDDO !faceNodeIdx
              ENDDO !componentIdx
            ENDDO !gaussPointIdx
          END DO !faceIdx
        END IF !decomposition%calculateFaces

      CASE DEFAULT
        ! Do nothing for other equation set subtypes
      END SELECT

    ENDIF
      
    EXITS("Darcy_FiniteElementFaceIntegrate")
    RETURN
999 ERRORSEXITS("Darcy_FiniteElementFaceIntegrate",err,error)
    RETURN 1
    
  END SUBROUTINE Darcy_FiniteElementFaceIntegrate

  !
  !================================================================================================================================
  !

  !>Sets the equation specification for a Darcy equation type of a fluid mechanics equations set class.
  SUBROUTINE Darcy_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the specification for
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations set specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: subtype
    TYPE(VARYING_STRING) :: localError

    ENTERS("Darcy_EquationsSetSpecificationSet",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(ALLOCATED(equationsSet%specification)) CALL FlagError("Equations set specification is already allocated.",err,error,*999)    
    IF(SIZE(specification,1)<3) THEN
      localError="The size of the specified specification array of "// &
        & TRIM(NumberToVString(SIZE(specification,1),"*",err,error))//" is invalid. The size should be >= 3."
      CALL FlagError(localError,err,error,*999)
    END IF
    
    subtype=specification(3)
    SELECT CASE(subtype)
    CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE, &
      & EQUATIONS_SET_QUASISTATIC_DARCY_SUBTYPE, &
      & EQUATIONS_SET_ALE_DARCY_SUBTYPE, &
      & EQUATIONS_SET_TRANSIENT_DARCY_SUBTYPE, &
      & EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
      & EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
      & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
      & EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE, &
      & EQUATIONS_SET_MULTI_COMPARTMENT_DARCY_SUBTYPE, &
      & EQUATIONS_SET_ELASTICITY_MULTI_COMPARTMENT_DARCY_INRIA_SUBTYPE, &
      & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
      !ok
    CASE DEFAULT
      localError="The third equations set specification of "//TRIM(NumberToVstring(specification(3),"*",err,error))// &
        & " is not valid for a Darcy type of a fluid mechanics equations set."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    !Set full specification
    ALLOCATE(equationsSet%specification(3),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate equations set specification.",err,error,*999)
    equationsSet%specification(1:3)=[EQUATIONS_SET_FLUID_MECHANICS_CLASS,EQUATIONS_SET_DARCY_EQUATION_TYPE,subtype]

    EXITS("Darcy_EquationsSetSpecificationSet")
    RETURN
999 ERRORS("Darcy_EquationsSetSpecificationSet",err,error)
    EXITS("Darcy_EquationsSetSpecificationSet")
    RETURN 1

  END SUBROUTINE Darcy_EquationsSetSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a Darcy problem.
  SUBROUTINE Darcy_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to set the problem specification for
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: problemSubtype
    TYPE(VARYING_STRING) :: localError

    ENTERS("Darcy_ProblemSpecificationSet",err,error,*998)

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
    CASE(PROBLEM_STANDARD_DARCY_SUBTYPE, &
      & PROBLEM_QUASISTATIC_DARCY_SUBTYPE, &
      & PROBLEM_ALE_DARCY_SUBTYPE, &
      & PROBLEM_TRANSIENT_DARCY_SUBTYPE, &
      & PROBLEM_PGM_DARCY_SUBTYPE, &
      & PROBLEM_PGM_TRANSIENT_DARCY_SUBTYPE)
      !All ok
    CASE DEFAULT
      localError="The third problem subtype of "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
        & " is not valid for a Darcy type of a fluid mechanics problem."
      CALL FlagError(localError,err,error,*998)
    END SELECT
    ALLOCATE(problem%specification(3),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
    problem%specification(1:3)=[PROBLEM_FLUID_MECHANICS_CLASS,PROBLEM_DARCY_EQUATION_TYPE,problemSubtype]

    EXITS("Darcy_ProblemSpecificationSet")
    RETURN
999 IF(ALLOCATED(problem%specification)) DEALLOCATE(problem%specification)
998 ERRORS("Darcy_ProblemSpecificationSet",err,error)
    EXITS("Darcy_ProblemSpecificationSet")
    RETURN 1

  END SUBROUTINE Darcy_ProblemSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the Darcy equations problem.
  SUBROUTINE Darcy_ProblemSetup(problem,problemSetup,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to setup
    TYPE(ProblemSetupType), INTENT(INOUT) :: problemSetup !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfSolvers,pSpecification(3),problemSubType,solverIdx
    TYPE(ControlLoopType), POINTER :: controlLoop,controlLoopRoot
    TYPE(SolverType), POINTER :: solver, solverMatProperties
    TYPE(SolverEquationsType), POINTER :: solverEquations, solverEquationsMatProperties
    TYPE(SolversType), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError

    ENTERS("Darcy_ProblemSetup",err,error,*999)

    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)

    problemSubType=pSpecification(3)
    SELECT CASE(problemSubType)
    CASE(PROBLEM_STANDARD_DARCY_SUBTYPE, &
      & PROBLEM_QUASISTATIC_DARCY_SUBTYPE, &
      & PROBLEM_ALE_DARCY_SUBTYPE, &
      & PROBLEM_TRANSIENT_DARCY_SUBTYPE, &
      & PROBLEM_PGM_DARCY_SUBTYPE, &
      & PROBLEM_PGM_TRANSIENT_DARCY_SUBTYPE)
      !All ok
    CASE DEFAULT
      localError="The third problem subtype of "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
        & " is not valid for a Darcy type of a fluid mechanics problem."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
      !-----------------------------------------------------------------
      !   s t a n d a r d   D a r c y
      !-----------------------------------------------------------------
      !CASE(PROBLEM_STANDARD_DARCY_SUBTYPE)
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
          & " is invalid for a standard Darcy equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_SETUP_CONTROL_TYPE)
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)
        !Set up a simple control loop
        NULLIFY(controlLoop)
        CALL ControlLoop_CreateStart(problem,controlLoop,err,error,*999)
        IF(problemSubType/=PROBLEM_STANDARD_DARCY_SUBTYPE) &
          & CALL ControlLoop_TypeSet(controlLoop,CONTROL_TIME_LOOP_TYPE,err,error,*999)
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
          & " is invalid for a standard Darcy equation."
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
        SELECT CASE(problemSubType)
        CASE(PROBLEM_STANDARD_DARCY_SUBTYPE, &
          & PROBLEM_QUASISTATIC_DARCY_SUBTYPE)
          CALL Solvers_NumberOfSolversSet(solvers,1,err,error,*999)
          !Set the solver to be a linear solver
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          CALL Solver_TypeSet(solver,SOLVER_LINEAR_TYPE,err,error,*999)
          !Set solver defaults
          CALL Solver_LibraryTypeSet(solver,SOLVER_PETSC_LIBRARY,err,error,*999)
        CASE(PROBLEM_ALE_DARCY_SUBTYPE, &
          & PROBLEM_PGM_DARCY_SUBTYPE)
          CALL Solvers_NumberOfSolversSet(solvers,2,err,error,*999)
          !Set the first solver to be a linear solver for the material update
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          CALL Solver_TypeSet(solver,SOLVER_LINEAR_TYPE,err,error,*999)
          CALL Solver_LibraryTypeSet(solver,SOLVER_PETSC_LIBRARY,err,error,*999)
          !Set the second solver to be a linear solver for the ALE Darcy
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
          CALL Solver_TypeSet(solver,SOLVER_LINEAR_TYPE,err,error,*999)
          CALL Solver_LibraryTypeSet(solver,SOLVER_PETSC_LIBRARY,err,error,*999)
        CASE(PROBLEM_PGM_TRANSIENT_DARCY_SUBTYPE)
          CALL Solvers_NumberOfSolversSet(solvers,2,err,error,*999)
          !Set the first solver to be a linear solver for the material update
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          CALL Solver_TypeSet(solver,SOLVER_LINEAR_TYPE,err,error,*999)
          CALL Solver_LibraryTypeSet(solver,SOLVER_PETSC_LIBRARY,err,error,*999)
          !Set the second solver to be a first order dynamic solver for the ALE Darcy
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
          CALL Solver_TypeSet(solver,SOLVER_DYNAMIC_TYPE,err,error,*999)
          CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
          !Set solver defaults
          CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
          CALL Solver_DynamicSchemeSet(solver,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,err,error,*999)
          CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
        CASE(PROBLEM_TRANSIENT_DARCY_SUBTYPE)
          CALL Solvers_NumberOfSolversSet(solvers,1,err,error,*999)
          !Set the solver to be a first order dynamic solver
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          CALL Solver_TypeSet(solver,SOLVER_DYNAMIC_TYPE,err,error,*999)
          CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
          !Set solver defaults
          CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
          CALL Solver_DynamicSchemeSet(solver,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,err,error,*999)
          CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
        CASE DEFAULT
          localError="The problem subtype of "//TRIM(NumberToVString(problemSubType,"*",err,error))// &
            & " is invalid for a Darcy equation subtype."
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
          & " is invalid for a standard Darcy equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT

    CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)
        !Get the control loop
        NULLIFY(controlLoopRoot)
        CALL Problem_ControlLoopRoot(problem,controlLoopRoot,err,error,*999)
        NULLIFY(controlLoop)
        CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
        !Get the solver
        NULLIFY(solvers)
        CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
        SELECT CASE(problemSubType)
        CASE(PROBLEM_STANDARD_DARCY_SUBTYPE, &
          & PROBLEM_QUASISTATIC_DARCY_SUBTYPE)
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          !Create the solver equations
          NULLIFY(solverEquations)
          CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
          CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_LINEAR,err,error,*999)
          IF(problemSubType==PROBLEM_QUASISTATIC_DARCY_SUBTYPE) THEN
            CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_QUASISTATIC,err,error,*999)
          ELSE
            CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_STATIC,err,error,*999)
          ENDIF
          CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
        CASE(PROBLEM_ALE_DARCY_SUBTYPE, &
          & PROBLEM_PGM_DARCY_SUBTYPE)
          DO solverIdx=1,2
            !solverIdx=1 is material-properties solver and solverIdx=2 is the Darcy-ALE solver equations
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,solverIdx,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
            CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_LINEAR,err,error,*999)
            CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_QUASISTATIC,err,error,*999)
            CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
          ENDDO !solverIdx
        CASE(PROBLEM_PGM_TRANSIENT_DARCY_SUBTYPE)
          !Get the material-properties solver and create the material-properties solver equations
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          NULLIFY(solverEquations)
          CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
          CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_LINEAR,err,error,*999)
          CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_QUASISTATIC,err,error,*999)
          CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
          !Get the Darcy-ALE solver and create the Darcy-ALE solver equations
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
          NULLIFY(solverEquations)
          CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
          CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_LINEAR,err,error,*999)
          CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
          CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
        CASE(PROBLEM_TRANSIENT_DARCY_SUBTYPE)
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          !Create the solver equations
          NULLIFY(solverEquations)
          CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
          CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_LINEAR,err,error,*999)
          CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
          CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
        CASE DEFAULT
          localError="The problem subtype of "//TRIM(NumberToVString(problemSubType,"*",err,error))// &
            & " is invalid for a Darcy equation subtype."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Get the control loop
        NULLIFY(controlLoopRoot)
        CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
        NULLIFY(controlLoop)
        CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
        !Get the solver equations
        NULLIFY(solvers)
        CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
        CALL Solvers_NumberOfSolversGet(solvers,numberOfSolvers,err,error,*999)
        DO solverIdx=1,numberOfSolvers
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,solverIdx,solver,err,error,*999)
          NULLIFY(solverEquations)
          CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
          !Finish the solver equations creation
          CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)
        ENDDO !solverIdx
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a standard Darcy equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
        & " is invalid for a standard Darcy equation."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Darcy_ProblemSetup")
    RETURN
999 ERRORSEXITS("Darcy_ProblemSetup",err,error)
    RETURN 1
    
  END SUBROUTINE Darcy_ProblemSetup

  !
  !=============================================================================================================================
  !

  !>Sets up the Darcy problem pre-solve.
  SUBROUTINE Darcy_PreSolve(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: analyticFunctionType,loopType,matrixNumber,numberOfSolverMatrices,pSpecification(3),problemSubType, &
      & solveNumber,solverMatrixIdx,solverNumber,solverNumberDarcy,solverNumberMatProperties,solverNumberSolid
    TYPE(ControlLoopType), POINTER :: controlLoop 
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(EquationsSetAnalyticType), POINTER :: equationsAnalytic
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverType), POINTER :: solverALEDarcy
    TYPE(SolversType), POINTER :: solvers
    TYPE(SolverEquationsType), POINTER :: solverEquations,solverEquationsALEDarcy
    TYPE(SolverMappingType), POINTER :: solverMapping,solverMappingALEDarcy
    TYPE(SolverMatricesType), POINTER :: solverMatrices
    TYPE(SolverMatrixType), POINTER :: solverMatrix
    TYPE(VARYING_STRING) :: localError

    ENTERS("Darcy_PreSolve",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    problemSubType=pSpecification(3)
    
    SELECT CASE(problemSubType)
    CASE(PROBLEM_STANDARD_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_DARCY_SUBTYPE,PROBLEM_TRANSIENT_DARCY_SUBTYPE)
      solverNumberDarcy=1
    CASE(PROBLEM_ALE_DARCY_SUBTYPE,PROBLEM_PGM_DARCY_SUBTYPE,PROBLEM_PGM_TRANSIENT_DARCY_SUBTYPE)
      solverNumberMatProperties=1
      solverNumberDarcy=2
    CASE(PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE)
      solverNumberSolid=1
      solverNumberDarcy=1
    CASE(PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE,PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE, &
      & PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
      solverNumberSolid=1
      solverNumberMatProperties=1
      solverNumberDarcy=2
    END SELECT

    !--- Set explicitly 'solverMatrix%updateMatrix=.TRUE.'
    NULLIFY(solverEquations)
    CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
    NULLIFY(solverMapping)
    CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
    NULLIFY(solverMatrices)
    CALL SolverEquations_SolverMatricesGet(solverEquations,solverMatrices,err,error,*999)
    CALL SolverMatrices_NumberOfMatricesGet(solverMatrices,numberOfSolverMatrices,err,error,*999)
    DO solverMatrixIdx=1,numberOfSolverMatrices
      NULLIFY(solverMatrix)
      CALL SolverMatrices_SolverMatrixGet(solverMatrices,solverMatrixIdx,solverMatrix,err,error,*999)
      solverMatrix%updateMatrix=.TRUE.
    ENDDO !solverMatrixIdx
 
    !--- pre_solve calls for various actions
    SELECT CASE(problemSubType)
    CASE(PROBLEM_STANDARD_DARCY_SUBTYPE)
      ! do nothing
    CASE(PROBLEM_QUASISTATIC_DARCY_SUBTYPE)
      ! do nothing
    CASE(PROBLEM_TRANSIENT_DARCY_SUBTYPE)
      CALL Solver_NumberGet(solver,solverNumber,err,error,*999)
      IF(solveNumber==solverNumberDarcy) CALL Darcy_PreSolveUpdateBoundaryConditions(solver,err,error,*999)
    CASE(PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE)
      CALL ControlLoop_TypeGet(controlLoop,loopType,err,error,*999)
      CALL Solver_NumberGet(solver,solverNumber,err,error,*999)
      IF((loopType==CONTROL_SIMPLE_TYPE.OR.loopType==CONTROL_TIME_LOOP_TYPE).AND.solverNumber==solverNumberDarcy) THEN
!!TODO remove these debug flags
        !--- flags to ensure once-per-time-step output in conjunction with diagnostics
        idebug1 = .TRUE.
        idebug2 = .TRUE.
        idebug3 = .TRUE.

        NULLIFY(solvers)
        CALL Solver_SolversGet(solver,solvers,err,error,*999)
        NULLIFY(solverALEDarcy)
        CALL Solvers_SolverGet(solvers,solverNumberDarcy,solverALEDarcy,err,error,*999)
        NULLIFY(solverEquationsALEDarcy)
        CALL Solver_SolverEquationsGet(solverALEDarcy,solverEquationsALEDarcy,err,error,*999)
        NULLIFY(solverMappingALEDarcy)
        CALL SolverEquations_SolverMappingGet(solverEquationsALEDarcy,solverMappingALEDarcy,err,error,*999)
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMappingALEDarcy,1,equationsSet,err,error,*999)
        NULLIFY(equationsAnalytic)
        CALL EquationsSet_AnalyticExists(equationsSet,equationsAnalytic,err,error,*999)
        IF(ASSOCIATED(equationsAnalytic)) THEN
          CALL EquationsSet_AnalyticFunctionTypeGet(equationsSet,analyticFunctionType,err,error,*999)
          IF(analyticFunctionType==EQUATIONS_SET_INCOMP_ELAST_DARCY_ANALYTIC_DARCY)THEN
            !call only analytic update and DO NOT call the other standard pre-solve routines as the mesh does not require deform
            CALL Darcy_PreSolveUpdateAnalyticValues(solverALEDarcy,err,error,*999)
          ENDIF
        ELSE
          !default
          !--- 1.1 Transfer solid displacement to Darcy geometric field
          CALL Darcy_PreSolveGetSolidDisplacement(solverALEDarcy,err,error,*999)
          
          !--- 1.2 Update the mesh (and calculate boundary velocities) PRIOR to solving for new material properties
          IF(solverNumber==solverNumberDarcy) CALL Darcy_PreSolveALEUpdateMesh(solverALEDarcy,err,error,*999)
          
          ! ! ! ! i n   p r i n c i p l e   c u r r e n t l y   d o   n o t   n e e d   t o   u p d a t e   B C s
          ! ! ! !unless:
          ! ! ! !--- 1.3 Apply both normal and moving mesh boundary conditions, OR:
          ! ! ! !--- 1.3 (Iteratively) Render the boundary impermeable (ellipsoid, general curvilinear mesh)
          ! ! ! CALL Darcy_PreSolveUpdateBoundaryConditions(solverALEDarcy,err,error,*999)
        ENDIF
      ENDIF
    CASE(PROBLEM_ALE_DARCY_SUBTYPE,PROBLEM_PGM_DARCY_SUBTYPE,PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE, &
      & PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE,PROBLEM_PGM_TRANSIENT_DARCY_SUBTYPE, &
      & PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
      CALL ControlLoop_TypeGet(controlLoop,loopType,err,error,*999)
      CALL SolverMatrix_NumberGet(solverMatrix,matrixNumber,err,error,*999)
      IF((loopType==CONTROL_SIMPLE_TYPE.OR.loopType==CONTROL_TIME_LOOP_TYPE).AND.solverNumber==solverNumberMatProperties) THEN
!!TODO remove these debug flags
        !--- flags to ensure once-per-time-step output in conjunction with diagnostics
        idebug1 = .TRUE.
        idebug2 = .TRUE.
        idebug3 = .TRUE.
        
        NULLIFY(solvers)
        CALL Solver_SolversGet(solver,solvers,err,error,*999)
        NULLIFY(solverALEDarcy)
        CALL Solvers_SolverGet(solvers,solverNumberDarcy,solverALEDarcy,err,error,*999)
        NULLIFY(solverEquationsALEDarcy)
        CALL Solver_SolverEquationsGet(solverALEDarcy,solverEquationsALEDarcy,err,error,*999)
        NULLIFY(solverMappingALEDarcy)
        CALL SolverEquations_SolverMappingGet(solverEquationsALEDarcy,solverMappingALEDarcy,err,error,*999)
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMappingALEDarcy,1,equationsSet,err,error,*999)
        NULLIFY(equationsAnalytic)
        CALL EquationsSet_AnalyticExists(equationsSet,equationsAnalytic,err,error,*999)
        IF(ASSOCIATED(equationsAnalytic)) THEN
          CALL EquationsSet_AnalyticFunctionTypeGet(equationsSet,analyticFunctionType,err,error,*999)
          IF(analyticFunctionType==EQUATIONS_SET_INCOMP_ELAST_DARCY_ANALYTIC_DARCY) THEN
            !call only analytic update and DO NOT call the other standard pre-solve routines as the mesh does not require deform
            CALL Darcy_PreSolveUpdateAnalyticValues(solverALEDarcy,err,error,*999)
          ENDIF
        ELSE
          !default
          !--- 1.1 Transfer solid displacement to Darcy geometric field
          CALL Darcy_PreSolveGetSolidDisplacement(solverALEDarcy,err,error,*999)

          !--- 1.2 Update the mesh (and calculate boundary velocities) PRIOR to solving for new material properties
          CALL Darcy_PreSolveALEUpdateMesh(solverALEDarcy,err,error,*999)

          ! ! ! ! i n   p r i n c i p l e   c u r r e n t l y   d o   n o t   n e e d   t o   u p d a t e   B C s
          ! ! ! !--- 1.3 Apply both normal and moving mesh boundary conditions
          ! ! ! CALL Darcy_PreSolveUpdateBoundaryConditions(solverALEDarcy,err,error,*999)
        ENDIF
      ELSE IF((loopType==CONTROL_SIMPLE_TYPE.OR.loopType==CONTROL_TIME_LOOP_TYPE).AND.solverNumber==solverNumberDarcy) THEN
        ! ! ! !n o t   f o r   n o w   ! ! !
        ! ! ! !--- 2.1 Update the material field
        ! ! ! CALL Darcy_PreSolveUpdateMatrixProperties(solver,solverNumberDarcy,solverNumberMatProperties,err,error,*999)
      END IF
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(problemSubType,"*",err,error))// &
        & " is not valid for a Darcy fluid type of a fluid mechanics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Darcy_PreSolve")
    RETURN
999 ERRORSEXITS("Darcy_PreSolve",err,error)
    RETURN 1
    
  END SUBROUTINE Darcy_PreSolve

  !
  !================================================================================================================================
  !

  SUBROUTINE Darcy_PreLoop(controlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the time control loop for the Darcy  problem
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: iterationNumber,outputType,pSpecification(3),solverNumber,solverNumberDarcy,solverNumberMatProperties, &
      & solverNumberSolid
    TYPE(ControlLoopType), POINTER :: controlLoopDarcy
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverType), POINTER :: solverDarcy
    TYPE(SolverSType), POINTER :: solvers

    ENTERS("Darcy_PreLoop",err,error,*999)

    !Get the solver for the Darcy problem
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STANDARD_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_DARCY_SUBTYPE,PROBLEM_TRANSIENT_DARCY_SUBTYPE)
      !solverNumberDarcy has to be set here so that store_reference_data and store_previous_data have access to it
      solverNumberDarcy=1
      NULLIFY(controlLoopDarcy)
      CALL ControlLoop_Get(controlLoop,[CONTROL_LOOP_NODE],controlLoopDarcy,err,error,*999)
      NULLIFY(solvers)
      CALL ControlLoop_SolversGet(controlLoopDarcy,solvers,err,error,*999)
      NULLIFY(solverDarcy)
      CALL Solvers_SolverGet(solvers,solverNumberDarcy,solverDarcy,err,error,*999)
    CASE(PROBLEM_ALE_DARCY_SUBTYPE,PROBLEM_PGM_DARCY_SUBTYPE,PROBLEM_PGM_TRANSIENT_DARCY_SUBTYPE)
      solverNumberMatProperties=1
      solverNumberDarcy=2
      NULLIFY(controlLoopDarcy)
      CALL ControlLoop_Get(controlLoop,[CONTROL_LOOP_NODE],controlLoopDarcy,err,error,*999)
      NULLIFY(solvers)
      CALL ControlLoop_SolversGet(controlLoopDarcy,solvers,err,error,*999)
      NULLIFY(solverDarcy)
      CALL Solvers_SolverGet(solvers,solverNumberDarcy,solverDarcy,err,error,*999)
    CASE(PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE,PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE)
      solverNumberSolid=1
      solverNumberMatProperties=1
      solverNumberDarcy=2
      NULLIFY(controlLoopDarcy)
      CALL ControlLoop_Get(controlLoop,[2,CONTROL_LOOP_NODE],controlLoopDarcy,err,error,*999)
      NULLIFY(solvers)
      CALL ControlLoop_SolversGet(controlLoopDarcy,solvers,err,error,*999)
      NULLIFY(solverDarcy)
      CALL Solvers_SolverGet(solvers,solverNumberDarcy,solverDarcy,err,error,*999)
    CASE(PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE)
      solverNumberSolid=1
      solverNumberDarcy=1
      NULLIFY(controlLoopDarcy)
      CALL ControlLoop_Get(controlLoop,[1,2,CONTROL_LOOP_NODE],controlLoopDarcy,err,error,*999)
      NULLIFY(solvers)
      CALL ControlLoop_SolversGet(controlLoopDarcy,solvers,err,error,*999)
      NULLIFY(solverDarcy)
      CALL Solvers_SolverGet(solvers,solverNumberDarcy,solverDarcy,err,error,*999)
    CASE(PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
      solverNumberSolid=1
      solverNumberMatProperties=1
      solverNumberDarcy=2
      NULLIFY(controlLoopDarcy)
      CALL ControlLoop_Get(controlLoop,[1,2,CONTROL_LOOP_NODE],controlLoopDarcy,err,error,*999)
      NULLIFY(solvers)
      CALL ControlLoop_SolversGet(controlLoopDarcy,solvers,err,error,*999)
      NULLIFY(solverDarcy)
      CALL Solvers_SolverGet(solvers,solverNumberDarcy,solverDarcy,err,error,*999)
    END SELECT

    CALL Solver_NumberGet(solverDarcy,solverNumber,err,error,*999)
    CALL ControlLoop_IterationNumberGet(controlLoop,iterationNumber,err,error,*999)
    CALL ControlLoop_OutputTypeGet(controlLoop,outputType,err,error,*999)
    !If this is the first time step then store reference data
    IF(iterationNumber==0) THEN
      IF(outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
        CALL WriteString(GENERAL_OUTPUT_TYPE,'== Storing reference data',err,error,*999)
      ENDIF
      IF(solverNumber==solverNumberDarcy) CALL Darcy_PreSolveStoreReferenceData(solverDarcy,err,error,*999)
    ENDIF

    !Store data of previous time step (mesh position); executed once per time step before subiteration
    IF(outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
      CALL WriteString(GENERAL_OUTPUT_TYPE,'== Storing previous data',err,error,*999)
    ENDIF
    IF(solverNumber==solverNumberDarcy) CALL Darcy_PreSolveStorePreviousData(solverDarcy,err,error,*999)

    EXITS("Darcy_PreLoop")
    RETURN
999 ERRORSEXITS("Darcy_PreLoop",err,error)
    RETURN 1
    
  END SUBROUTINE Darcy_PreLoop

  !
  !================================================================================================================================
  !

  !>Store some reference data for ALE Darcy problem
  SUBROUTINE Darcy_PreSolveStoreReferenceData(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx,esSpecification(3),numberDOFsToPrint,numberOfEquationsSets,pSpecification(3),variableType
    REAL(DP) :: alpha
    REAL(DP), POINTER :: initialValues(:)
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the control loop to solve.
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField, geometricField
    TYPE(FieldVariableType), POINTER :: fieldVariable
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverEquationsType), POINTER :: solverEquations  !<A pointer to the solver equations
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping
    TYPE(VARYING_STRING) :: localError

    ENTERS("Darcy_PreSolveStoreReferenceData",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)

    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STANDARD_DARCY_SUBTYPE)
      ! do nothing
    CASE(PROBLEM_QUASISTATIC_DARCY_SUBTYPE)
      ! do nothing
    CASE(PROBLEM_TRANSIENT_DARCY_SUBTYPE)
      ! do nothing
    CASE(PROBLEM_ALE_DARCY_SUBTYPE,PROBLEM_PGM_DARCY_SUBTYPE,PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE, &
      & PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE,PROBLEM_PGM_TRANSIENT_DARCY_SUBTYPE, &
      & PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
      NULLIFY(solverEquations)
      CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
      DO equationsSetIdx=1,numberOfEquationsSets
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
        CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
        SELECT CASE(esSpecification(3))
        CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE)
          ! do nothing
        CASE(EQUATIONS_SET_QUASISTATIC_DARCY_SUBTYPE)
          ! do nothing
        CASE(EQUATIONS_SET_TRANSIENT_DARCY_SUBTYPE)
          ! do nothing
        CASE(EQUATIONS_SET_ALE_DARCY_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
          NULLIFY(geometricField)
          CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
          NULLIFY(dependentField)
          CALL EquationstSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
          !--- Store the initial (= reference) geometry field values
          alpha = 1.0_DP
          CALL Field_ParameterSetsCopy(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & FIELD_INITIAL_VALUES_SET_TYPE,alpha,err,error,*999)
          NULLIFY(equations)
          CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
          NULLIFY(vectorEquations)
          CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
          NULLIFY(vectorMapping)
          CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
          SELECT CASE(esSpecification(3))
          CASE(EQUATIONS_SET_ALE_DARCY_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE)
            NULLIFY(linearMapping)
            CALL EquationsMappingVector_LinearMappingGet(vectorMapping,linearMapping,err,error,*999)
            NULLIFY(fieldVariable)
            CALL EquationsMappingLinear_LinearMatrixVariableGet(linearMapping,1,fieldVariable,err,error,*999)
            ! '1' associated with linear matrix
          CASE(EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
            & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
            & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
            NULLIFY(dynamicMapping)
            CALL EquationsMappingVector_DynamicMappingGet(vectorMapping,dynamicMapping,err,error,*999)
            NULLIFY(fieldVariable)
            CALL EquationsMappingDynamic_DynamicVariableGet(dynamicMapping,fieldVariable,err,error,*999)
          END SELECT

          CALL FieldVariable_VariableTypeGet(fieldVariable,variableType,err,error,*999)
          !--- Store the initial DEPENDENT field values
          alpha = 1.0_DP
          CALL FieldVariable_ParameterSetsCopy(fieldVariable,FIELD_VALUES_SET_TYPE,FIELD_INITIAL_VALUES_SET_TYPE, &
            & alpha,err,error,*999)

          IF(diagnostics1) THEN
            NULLIFY(initialValues)
            CALL FieldVariable_ParameterSetDataGet(fieldVariable,FIELD_INITIAL_VALUES_SET_TYPE,initialValues,err,error,*999)
            numberDOFsToPrint = SIZE(initialValues,1)
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberDOFsToPrint,numberDOFsToPrint,numberDOFsToPrint, &
              & initialValues,'(" dependentField,FIELD_U_VARIABLE_TYPE,FIELD_INITIAL_VALUES_SET_TYPE = ",4(X,E13.6))', &
              & '4(4(X,E13.6))',err,error,*999)
            CALL FieldVariable_ParameterSetDataRestore(fieldVariable,FIELD_INITIAL_VALUES_SET_TYPE,initialValues,err,error,*999)
          ENDIF
        CASE DEFAULT
          localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
            & " is not valid for a Darcy equation fluid type of a fluid mechanics problem class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDDO !equationsSetIdx
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a Darcy equation fluid type of a fluid mechanics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Darcy_PreSolveStoreReferenceData")
    RETURN
999 ERRORSEXITS("Darcy_PreSolveStoreReferenceData",err,error)
    RETURN 1

  END SUBROUTINE Darcy_PreSolveStoreReferenceData

  !
  !================================================================================================================================
  !

  !>Store data of previous time step (mesh position) for ALE Darcy problem
  SUBROUTINE Darcy_PreSolveStorePreviousData(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx,esSpecification(3),numberOfEquationsSets,pSpecification(3)
    REAL(DP) :: alpha
    TYPE(ControlLoopType), POINTER :: controlLoop 
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set
    TYPE(FieldType), POINTER :: geometricField
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverEquationsType), POINTER :: solverEquations  !<A pointer to the solver equations
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping
    TYPE(VARYING_STRING) :: localError

    ENTERS("Darcy_PreSolveStorePreviousData",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)

    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STANDARD_DARCY_SUBTYPE)
      ! do nothing
    CASE(PROBLEM_QUASISTATIC_DARCY_SUBTYPE)
      ! do nothing
    CASE(PROBLEM_TRANSIENT_DARCY_SUBTYPE)
      ! do nothing
    CASE(PROBLEM_ALE_DARCY_SUBTYPE,PROBLEM_PGM_DARCY_SUBTYPE,PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE, &
      & PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE,PROBLEM_PGM_TRANSIENT_DARCY_SUBTYPE, &
      & PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
      NULLIFY(solverEquations)
      CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
      DO equationsSetIdx=1,numberOfEquationsSets
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
        CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
        SELECT CASE(esSpecification(3))
        CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE)
          ! do nothing
        CASE(EQUATIONS_SET_QUASISTATIC_DARCY_SUBTYPE)
          ! do nothing
        CASE(EQUATIONS_SET_TRANSIENT_DARCY_SUBTYPE)
          ! do nothing
        CASE(EQUATIONS_SET_ALE_DARCY_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
          NULLIFY(geometricField)
          CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
          !--- Store the geometry field values of the previous time step
          alpha = 1.0_DP
          CALL Field_ParameterSetsCopy(geometricField,FIELD_U_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE,alpha,err,error,*999)
        CASE DEFAULT
          localError="Equations set subtype " &
            & //TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
            & " is not valid for a Darcy equation fluid type of a fluid mechanics problem class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDDO !equationsSetIdx
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a Darcy equation fluid type of a fluid mechanics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Darcy_PreSolveStorePreviousData")
    RETURN
999 ERRORSEXITS("Darcy_PreSolveStorePreviousData",err,error)
    RETURN 1

  END SUBROUTINE Darcy_PreSolveStorePreviousData

  !
  !================================================================================================================================
  !

  !>Update mesh position and velocity for ALE Darcy problem
  SUBROUTINE Darcy_PreSolveALEUpdateMesh(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dofNumber,esSpecification(3),loopIdx,numberOfDofs,numberDOFsToPrint,problemSubtype,pSpecification(3), &
      & outputType
    REAL(DP) :: currentTime,timeIncrement,alpha
    REAL(DP), POINTER :: meshDisplacementValues(:)
    TYPE(ControlLoopType), POINTER :: controlLoop,controlTimeLoop
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set
    TYPE(FieldType), POINTER :: geometricField
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverType), POINTER :: solverALEDarcy !<A pointer to the solvers
    TYPE(SolverEquationsType), POINTER :: solverEquations  !<A pointer to the solver equations
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping
    TYPE(VARYING_STRING) :: localError

    ENTERS("Darcy_PreSolveALEUpdateMesh",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    
    CALL ControlLoop_CurrentTimesGet(controlLoop,currentTime,timeIncrement,err,error,*999)
    problemSubtype=pSpecification(3)
    
    SELECT CASE(problemSubtype)
    CASE(PROBLEM_STANDARD_DARCY_SUBTYPE)
      ! do nothing
    CASE(PROBLEM_QUASISTATIC_DARCY_SUBTYPE)
      ! do nothing
    CASE(PROBLEM_TRANSIENT_DARCY_SUBTYPE)
      ! do nothing
    CASE(PROBLEM_ALE_DARCY_SUBTYPE,PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE,PROBLEM_PGM_DARCY_SUBTYPE, &
      & PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE,PROBLEM_PGM_TRANSIENT_DARCY_SUBTYPE, &
      & PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
      NULLIFY(solverEquations)
      CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
      CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
      CALL Solver_OutputTypeGet(solver,outputType,err,error,*999)
      SELECT CASE(esSpecification(3))
      CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE)
        ! do nothing
      CASE(EQUATIONS_SET_QUASISTATIC_DARCY_SUBTYPE)
        ! do nothing
      CASE(EQUATIONS_SET_ALE_DARCY_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
        & EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
        & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
        & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
        IF(outputType>=SOLVER_PROGRESS_OUTPUT) THEN
          CALL WriteString(GENERAL_OUTPUT_TYPE,"Darcy update mesh ... ",err,error,*999)
        ENDIF
        NULLIFY(geometricField)
        CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
        !--- First, get pointer to mesh displacement values
        NULLIFY(meshDisplacementValues)
        CALL Field_ParameterSetDataGet(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE, &
          & meshDisplacementValues,err,error,*999)
        IF(diagnostics1) THEN
          numberDOFsToPrint = SIZE(meshDisplacementValues,1)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberDOFsToPrint,numberDOFsToPrint,numberDOFsToPrint,&
            & meshDisplacementValues,'(" meshDisplacementValues = ",3(X,E13.6))','3(3(X,E13.6))', &
            & err,error,*999)
        ENDIF

        CALL Field_NumberOfDOFsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDOFs,err,error,*999)
        
        IF(problemSubtype==PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE &
          & .OR. problemSubtype==PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE &
          & .OR. problemSubtype==PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE) THEN
          !--- Don't update geometric field here, this is done in
          !    darcy_equation_pre_solve_get_solid_displacement for these problems, but
          !    needs to be made consistent between the different problem types
        ELSE
          !--- Second, update geometric field
          DO dofNumber=1,numberOfDofs
            CALL Field_ParameterSetAddLocalDOF(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,dofNumber, &
              & meshDisplacementValues(dofNumber),err,error,*999)
          END DO
          CALL Field_ParameterSetUpdateStart(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
          CALL Field_ParameterSetUpdateFinish(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
        ENDIF
        
        !--- Third, use displacement values to calculate velocity values
        alpha=1.0_DP/timeIncrement
        CALL Field_ParameterSetsCopy(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE, &
          & FIELD_MESH_VELOCITY_SET_TYPE,alpha,err,error,*999)
        CALL Field_ParameterSetDataRestore(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE, &
          & meshDisplacementValues,err,error,*999)
      CASE DEFAULT
        localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
          & " is not valid for a Darcy equation fluid type of a fluid mechanics equations set."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(problemSubtype,"*",err,error))// &
        & " is not valid for a Darcy equation fluid type of a fluid mechanics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Darcy_PreSolveALEUpdateMesh")
    RETURN
999 ERRORSEXITS("Darcy_PreSolveALEUpdateMesh",err,error)
    RETURN 1
    
  END SUBROUTINE Darcy_PreSolveALEUpdateMesh

  !
  !================================================================================================================================
  !

  !>Update boundary conditions for Darcy equation pre solve
  SUBROUTINE Darcy_PreSolveUpdateBoundaryConditions(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: boundaryConditionCheckVariable,dofNumber,esSpecification(3),loopIdx,numberOfDofs,numberDOFsToPrint, &
      & outputType,problemSubType,pSpecification(3),variableType
    REAL(DP) :: currentTime, pressure,timeIncrement
    REAL(DP), POINTER :: dummyValues1(:),initialValues(:),meshVelocityValues(:)
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable
    TYPE(ControlLoopType), POINTER :: controlLoop,controlTimeLoop
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set
    TYPE(FieldType), POINTER :: dependentField, geometricField
    TYPE(FieldVariableType), POINTER :: fieldVariable
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverEquationsType), POINTER :: solverEquations  !<A pointer to the solver equations
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping
    TYPE(VARYING_STRING) :: localError

    ENTERS("Darcy_PreSolveUpdateBoundaryConditions",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    
    CALL ControlLoop_CurrentTimesGet(controlLoop,currentTime,timeIncrement,err,error,*999)
    problemSubtype=pSpecification(3)
    SELECT CASE(problemSubType)
    CASE(PROBLEM_STANDARD_DARCY_SUBTYPE)
      ! do nothing
    CASE(PROBLEM_QUASISTATIC_DARCY_SUBTYPE)
      ! do nothing
    CASE(PROBLEM_TRANSIENT_DARCY_SUBTYPE)
      ! do nothing
    CASE(PROBLEM_ALE_DARCY_SUBTYPE,PROBLEM_PGM_DARCY_SUBTYPE,PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE, &
      & PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE,PROBLEM_PGM_TRANSIENT_DARCY_SUBTYPE, &
      & PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
      NULLIFY(solverEquations)
      CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
      CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
      CALL Solver_OutputTypeGet(solver,outputType,err,error,*999)
      SELECT CASE(esSpecification(3))
      CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE)
        ! do nothing
      CASE(EQUATIONS_SET_QUASISTATIC_DARCY_SUBTYPE)
        ! do nothing
      CASE(EQUATIONS_SET_ALE_DARCY_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
        & EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
        & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
        & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
        IF(outputType>=SOLVER_PROGRESS_OUTPUT) THEN
          CALL WriteString(GENERAL_OUTPUT_TYPE,"Darcy update boundary conditions ... ",err,error,*999)
        ENDIF
        NULLIFY(geometricField)
        CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
        NULLIFY(dependentField)
        CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
        NULLIFY(boundaryConditions)
        CALL SolverEquations_BoundaryConditionsGet(solverEquations,boundaryConditions,err,error,*999)
        NULLIFY(equations)
        CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
        NULLIFY(vectorEquations)
        CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
        NULLIFY(vectorMapping)
        CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
        SELECT CASE(esSpecification(3))
        CASE(EQUATIONS_SET_ALE_DARCY_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE)
          NULLIFY(linearMapping)
          CALL EquationsMappingVector_LinearMappingGet(vectorMapping,linearMapping,err,error,*999)
          NULLIFY(fieldVariable)
          CALL EquationsMappingLinear_LinearMatrixVariableGet(linearMapping,1,fieldVariable,err,error,*999)
        CASE(EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE, &
          & EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
          NULLIFY(dynamicMapping)
          CALL EquationsMappingVector_DynamicMappingGet(vectorMapping,dynamicMapping,err,error,*999)
          NULLIFY(fieldVariable)
          CALL EquationsMappingDynamic_DynamicVariableGet(dynamicMapping,fieldVariable,err,error,*999)
        END SELECT

        CALL FieldVariable_VariableTypeGet(fieldVariable,variableType,err,error,*999)
        
        NULLIFY(boundaryConditionsVariable)
        CALL BoundaryConditions_VariableGet(boundaryConditions,fieldVariable,boundaryConditionsVariable,err,error,*999)
        NULLIFY(meshVelocityValues)
        CALL Field_ParameterSetDataGet(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_MESH_VELOCITY_SET_TYPE,meshVelocityValues, &
          & err,error,*999)
        NULLIFY(initialValues)
        CALL Field_ParameterSetDataGet(dependentField,variableType,FIELD_INITIAL_VALUES_SET_TYPE,initialValues,err,error,*999)
        IF(diagnostics1) THEN
          NULLIFY(dummyValues1)
          CALL Field_ParameterSetDataGet(dependentField,variableType,FIELD_VALUES_SET_TYPE,dummyValues1,err,error,*999)
          numberDOFsToPrint = SIZE(dummyValues1,1)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberDOFsToPrint,numberDOFsToPrint, &
            & numberDOFsToPrint,dummyValues1, &
            & '(" dependentField,variableType,FIELD_VALUES_SET_TYPE (before) = ",4(X,E13.6))', &
            & '4(4(X,E13.6))',err,error,*999)
        ENDIF
        CALL Field_NumberOfDOFsGet(dependentField,variableType,numberOfDOFs,err,error,*999)
        DO dofNumber=1,numberOfDofs
          boundaryConditionCheckVariable=boundaryConditionsVariable%conditionTypes(dofNumber)
          IF(boundaryConditionCheckVariable==BOUNDARY_CONDITION_MOVED_WALL) THEN
            !--- Reset boundary condition to the initial normal-velocity boundary condition
            CALL Field_ParameterSetUpdateLocalDOF(dependentField,variableType,FIELD_VALUES_SET_TYPE,dofNumber, &
              & initialValues(dofNumber),err,error,*999)
            !--- Add the velocity of the moving boundary on top of the initial boundary condition
!! === If we solve in terms of Darcy flow vector, then do not add mesh velocity === !!
!! === The BC is kept to the initial BC, for instance: null-flux                === !!
            !                                          CALL Field_ParameterSetAddLocalDOF(dependentField, &
            !                                            & variableType,FIELD_VALUES_SET_TYPE,dofNumber, &
            !                                            & meshVelocityValues(dofNumber),err,error,*999)
            !                                            ! dependent field      ( V_u, V_v, V_w, P_p )
            !                                            ! meshVelocityValues ( V_u, V_v, V_w )
            
          ELSE IF( boundaryConditionCheckVariable==BOUNDARY_CONDITION_FIXED .AND. &
            & esSpecification(3)==EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE) THEN
            !\ToDo: Check component number; this way we can also apply it to velocity
            !--- Set the time-dependent pressure BC
            pressure = initialValues(dofNumber) * (1.0_DP - EXP(- currentTime**2.0_DP / 0.25_DP))
            
            CALL Field_ParameterSetUpdateLocalDOF(dependentField,variableType,FIELD_VALUES_SET_TYPE,dofNumber, &
              & pressure,err,error,*999)
          ELSE
            ! do nothing
          END IF
        ENDDO !dofNumber
        CALL Field_ParameterSetUpdateStart(dependentField,variableType,FIELD_VALUES_SET_TYPE,err,error,*999)
        CALL Field_ParameterSetUpdateFinish(dependentField,variableType,FIELD_VALUES_SET_TYPE,err,error,*999)
        IF(diagnostics1) THEN
          numberDOFsToPrint = SIZE(meshVelocityValues,1)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberDOFsToPrint,numberDOFsToPrint, &
            & numberDOFsToPrint,meshVelocityValues,'(" meshVelocityValues = ",4(X,E13.6))','4(4(X,E13.6))',err,error,*999)
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE," ",err,error,*999)
          !
          numberDOFsToPrint = SIZE(initialValues,1)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberDOFsToPrint,numberDOFsToPrint, &
            & numberDOFsToPrint,initialValues,'(" initialValues = ",4(X,E13.6))','4(4(X,E13.6))',err,error,*999)
          !
          NULLIFY( dummyValues1 )
          CALL Field_ParameterSetDataGet(dependentField,variableType,FIELD_VALUES_SET_TYPE,dummyValues1,err,error,*999)
          numberDOFsToPrint = SIZE(dummyValues1,1)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberDOFsToPrint,numberDOFsToPrint, &
            & numberDOFsToPrint,dummyValues1,'(" dependentField,variableType,FIELD_VALUES_SET_TYPE (after) = ",4(X,E13.6))', &
            & '4(4(X,E13.6))',err,error,*999)
          CALL Field_ParameterSetDataRestore(dependentField,variableType,FIELD_VALUES_SET_TYPE,dummyValues1,err,error,*999)
        ENDIF
        CALL Field_ParameterSetDataRestore(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_MESH_VELOCITY_SET_TYPE, &
          & meshVelocityValues,err,error,*999)
        CALL Field_ParameterSetDataRestore(dependentField,variableType,FIELD_INITIAL_VALUES_SET_TYPE,initialValues, &
          & err,error,*999)        
        CALL Field_ParameterSetUpdateStart(dependentField,variableType,FIELD_VALUES_SET_TYPE,err,error,*999)
        CALL Field_ParameterSetUpdateFinish(dependentField,variableType,FIELD_VALUES_SET_TYPE,err,error,*999)
      CASE DEFAULT
        localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
          & " is not valid for a Darcy equation fluid type of a fluid mechanics problem class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(problemSubType,"*",err,error))// &
        & " is not valid for a Darcy equation fluid type of a fluid mechanics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Darcy_PreSolveUpdateBoundaryConditions")
    RETURN
999 ERRORS("Darcy_PreSolveUpdateBoundaryConditions",err,error)
    EXITS("Darcy_PreSolveUpdateBoundaryConditions")
    RETURN 1

  END SUBROUTINE Darcy_PreSolveUpdateBoundaryConditions

  !
  !================================================================================================================================
  !

  !>Update materials field for ALE Darcy problem
  SUBROUTINE Darcy_PreSolveUpdateMatrixProperties(solver,solverNumberDarcy,solverNumberMatProperties,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solvers
    INTEGER(INTG), INTENT(IN) :: solverNumberDarcy !<The solver number corresponding to the Darcy solver
    INTEGER(INTG), INTENT(IN) :: solverNumberMatProperties !<The solver number corresponding to the material properties solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx,numberOfComponentsMatProperties,numberOfComponentsALEDarcy,numberDOFsToPrint, &
      & pSpecification(3),outputType
    REAL(DP), POINTER :: dummyValues2(:)
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(EquationsSetType), POINTER :: equationsSetMatProperties,equationsSetALEDarcy
    TYPE(FieldType), POINTER :: dependentFieldMatProperties,materialsFieldALEDarcy
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverType), POINTER :: solverMatProperties,solverALEDarcy
    TYPE(SolversType), POINTER :: solvers
    TYPE(SolverEquationsType), POINTER :: solverEquationsMatProperties,solverEquationsALEDarcy
    TYPE(SolverMappingType), POINTER :: solverMappingMatProperties,solverMappingALEDarcy
    TYPE(VARYING_STRING) :: localError

    ENTERS("Darcy_PreSolveUpdateMatrixProperties",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STANDARD_DARCY_SUBTYPE)
      ! do nothing
    CASE(PROBLEM_QUASISTATIC_DARCY_SUBTYPE)
      ! do nothing
    CASE(PROBLEM_TRANSIENT_DARCY_SUBTYPE)
      ! do nothing
    CASE(PROBLEM_ALE_DARCY_SUBTYPE,PROBLEM_PGM_DARCY_SUBTYPE,PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE, &
      & PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE,PROBLEM_PGM_TRANSIENT_DARCY_SUBTYPE, &
      & PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
      IF(controlLoop%loopType==CONTROL_SIMPLE_TYPE.OR.controlLoop%loopType==CONTROL_TIME_LOOP_TYPE) THEN
        !--- Get the dependent field of the Material-Properties Galerkin-Projection equations
        CALL Solver_OutputTypeGet(solver,outputType,err,error,*999)
        IF(outputType>=SOLVER_PROGRESS_OUTPUT) THEN
          CALL WriteString(GENERAL_OUTPUT_TYPE,"Darcy update materials ... ",err,error,*999)
        ENDIF
        NULLIFY(solvers)
        CALL Solver_SolversGet(solver,solvers,err,error,*999)
        NULLIFY(solverMatProperties)
        CALL Solvers_SolverGet(solvers,solverNumberMatProperties,solverMatProperties,err,error,*999)
        NULLIFY(solverEquationsMatProperties)
        CALL Solver_SolverEquationsGet(solverMatProperties,solverEquationsMatProperties,err,error,*999)
        NULLIFY(solverMappingMatProperties)
        CALL SolverEquations_SolverMappingGet(solverEquationsMatProperties,solverMappingMatProperties,err,error,*999)
        NULLIFY(equationsSetMatProperties)
        CALL SolverMapping_EquationsSetGet(solverMappingMatProperties,1,equationsSetMatProperties,err,error,*999)
        NULLIFY(dependentFieldMatProperties)
        CALL EquationsSet_DependentFieldGet(equationsSetMatProperties,dependentFieldMatProperties,err,error,*999)
        CALL Field_NumberOfComponentsGet(dependentFieldMatProperties,FIELD_U_VARIABLE_TYPE,numberOfComponentsMatProperties, &
          & err,error,*999)
 
        !--- Get the materials field for the ALE Darcy equations
        NULLIFY(solverALEDarcy)
        CALL Solvers_SolverGet(solvers,solverNumberDarcy,solverALEDarcy,err,error,*999)
        NULLIFY(solverEquationsALEDarcy)
        CALL Solver_SolverEquationGet(solverALEDarcy,solverEquationsALEDarcy,err,error,*999)
        NULLIFY(solverMappingALEDarcy)
        CALL SolverEquations_SolverMappingGet(solverEquationsALEDarcy,solverMappingALEDarcy,err,error,*999)
        NULLIFY(equationsSetALEDarcy)
        CALL SolverMapping_EquationsSetGet(solverMappingALEDarcy,1,equationsSetALEDarcy,err,error,*999)
        NULLIFY(materialsFieldALEDarcy)
        CALL EquationsSet_MaterialsFieldGet(equationsSetALEDarcy,materialsFieldALEDarcy,err,error,*999)
        CALL Field_NumberOfComponentsGet(materialsFieldALEDarcy,FIELD_U_VARIABLE_TYPE,numberOfComponentsALEDarcy,err,error,*999)
        
        !--- Copy the result from Galerkin-Projection's dependent field to ALE Darcy's material field
        IF(numberOfComponentsALEDarcy/=numberOfComponentsMatProperties) THEN
         localError="Number of components of Galerkin-Projection dependent field "// &
            & "is not consistent with ALE-Darcy-equation material field."
          CALL FlagError(localError,err,error,*999)
        END IF
        DO componentIdx=1,numberOfComponentsALEDarcy
          CALL Field_ParametersToFieldParametersCopy(dependentFieldMatProperties,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & componentIdx,materialsFieldALEDarcy,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,componentIdx,err,error,*999)
        ENDDO !componentIdx
 
        IF(diagnostics3) THEN
          NULLIFY( dummyValues2 )
          CALL Field_ParameterSetDataGet(dependentFieldMatProperties,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,dummyValues2, &
            & err,error,*999)
          numberDOFsToPrint = SIZE(dummyValues2,1)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberDOFsToPrint,numberDOFsToPrint,numberDOFsToPrint,dummyValues2, &
            & '(" DEPENDENT_FIELD_MAT_PROPERTIES,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE = ",4(X,E13.6))',&
            & '4(4(X,E13.6))',err,error,*999)
          CALL Field_ParameterSetDataRestore(dependentFieldMatProperties,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & dummyValues2,err,error,*999)
        ENDIF
      ENDIF
    CASE DEFAULT
      localError="The third problem specification of "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for Darcy_PreSolveUpdateMatrixProperties."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Darcy_PreSolveUpdateMatrixProperties")
    RETURN
999 ERRORS("Darcy_PreSolveUpdateMatrixProperties",err,error)
    EXITS("Darcy_PreSolveUpdateMatrixProperties")
    RETURN 1

  END SUBROUTINE Darcy_PreSolveUpdateMatrixProperties

  !
  !================================================================================================================================
  !

  !>Sets up the Darcy problem post solve.
  SUBROUTINE Darcy_PostSolve(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver!<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: pSpecification(3),solverNumber
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the control loop to solve.
    TYPE(ProblemType), POINTER :: problem
    TYPE(VARYING_STRING) :: localError

    ENTERS("Darcy_PostSolve",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    CALL Solver_GlobalNumberGet(solver,solverNumber,err,error,*999)

    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STANDARD_DARCY_SUBTYPE)
      CALL Darcy_PostSolveOutputData(solver,err,error,*999)
    CASE(PROBLEM_QUASISTATIC_DARCY_SUBTYPE)
      CALL Darcy_PostSolveOutputData(solver,err,error,*999)
    CASE(PROBLEM_ALE_DARCY_SUBTYPE,PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE)
      IF(solverNumber==2) CALL Darcy_PostSolveOutputData(solver,err,error,*999)
    CASE(PROBLEM_PGM_DARCY_SUBTYPE,PROBLEM_PGM_TRANSIENT_DARCY_SUBTYPE,PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE, &
      & PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
      IF(solverNumber==2) THEN
        CALL Darcy_PostSolveOutputData(solver,err,error,*999)
        
        ! The following command only when setting the Darcy mass increase explicitly to test finite elasticity !!!
        ! ! ! CALL Darcy_PostSolveSetMassIncrease(solver,err,error,*999)
        
      ENDIF
    CASE(PROBLEM_TRANSIENT_DARCY_SUBTYPE)
      CALL Darcy_PostSolveOutputData(solver,err,error,*999)
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a Darcy fluid type of a fluid mechanics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Darcy_PostSolve")
    RETURN
999 ERRORSEXITS("Darcy_PostSolve",err,error)
    RETURN 1
    
  END SUBROUTINE Darcy_PostSolve

  !
  !================================================================================================================================
  !

  !>Sets up the Darcy problem post solve output data.
  SUBROUTINE Darcy_PostSolveOutputData(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: controlLoopOutputType,currentIteration,equationsSetIdx,esSpecification(3),inputIteration,loopIdx, &
      & maximumNumberOfIterations,numberOfEquationsSets,parentLoopType,pSpecification(3),outputIteration,solverOutputType, &
      & subIterationNumber
    REAL(DP) :: absoluteTolerance,currentTime,relativeTolerance,startTime,stopTime,timeIncrement
    CHARACTER(14) :: outputFile
    LOGICAL :: continueLoop,exportField
    TYPE(ControlLoopType), POINTER :: controlLoop,controlTimeLoop,parentLoop
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(FieldsType), POINTER :: fields
    TYPE(ProblemType), POINTER :: problem
    TYPE(RegionType), POINTER :: region
    TYPE(SolverType), POINTER :: solverALEDarcy,solverMatProperties
    TYPE(SolversType), POINTER :: solvers
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(VARYING_STRING) :: localError,method,filename
  
    ENTERS("Darcy_PostSolveOutputData",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    CALL Solver_OutputTypeGet(solver,solverOutputType,err,error,*999)
    
    CALL SYSTEM('mkdir -p ./output')
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STANDARD_DARCY_SUBTYPE)
      NULLIFY(solverEquations)
      CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
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
        IF(solverOutputType>=SOLVER_PROGRESS_OUTPUT) THEN
          CALL WriteString(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
          CALL WriteString(GENERAL_OUTPUT_TYPE,"Now export fields... ",err,error,*999)
        ENDIF
        CALL FIELD_IO_NODES_EXPORT(Fields,filename,method,err,error,*999)
        CALL FIELD_IO_ELEMENTS_EXPORT(Fields,filename,method,err,error,*999)
      ENDDO !equationsSetIdx
    CASE(PROBLEM_QUASISTATIC_DARCY_SUBTYPE, PROBLEM_ALE_DARCY_SUBTYPE, PROBLEM_PGM_DARCY_SUBTYPE, &
      & PROBLEM_TRANSIENT_DARCY_SUBTYPE, PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE, &
      & PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE,PROBLEM_PGM_TRANSIENT_DARCY_SUBTYPE, &
      & PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
      CALL ControlLoop_CurrentTimeInformationGet(controlLoop,currentTime,timeIncrement,startTime,stopTime,currentIteration, &
        & outputIteration,inputIteration,err,error,*999)
      NULLIFY(solverEquations)
      CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
      !Make sure the equations sets are up to date
      DO equationsSetIdx=1,numberOfEquationsSets
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
        CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
        IF(esSpecification(2)==EQUATIONS_SET_DARCY_EQUATION_TYPE) THEN
          IF(esSpecification(3)==EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE) THEN
            NULLIFY(parentLoop)
            CALL ControlLoop_ParentLoopGet(controlLoop,parentLoop,err,error,*999)
            CALL ControlLoop_TypeGet(parentLoop,parentLoopType,err,error,*999)
            IF(parentLoopType==CONTROL_WHILE_LOOP_TYPE) &
              & CALL ControlLoop_WhileInformationGet(parentLoop,subIterationNumber,,maximumNumberOfIterations,absoluteTolerance, &
              & relativeTolerance,continueLoop,err,error,*999)            
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
                NULLIFY(region)
                CALL EquationsSet_RegionGet(equationsSet,region,err,error,*999)
                NULLIFY(fields)
                CALL Region_FieldsGet(region,fields,err,error,*999)
                filename="./output/"//"MainTime_"//TRIM(NumberToVString(currentIteration,"*",err,error))
                method="FORTRAN"
                CALL ControlLoop_OutputTypeGet(parentLoop,controlLoopOutputType,err,error,*999)
                IF(MOD(currentIteration,outputIteration)==0)  THEN
                  IF(controlLoopOutputType >= CONTROL_LOOP_PROGRESS_OUTPUT) THEN
                    CALL WriteString(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
                    CALL WriteString(GENERAL_OUTPUT_TYPE,"Now export fields... ",err,error,*999)
                  ENDIF
                  CALL FIELD_IO_NODES_EXPORT(fields,filename,method,err,error,*999)
                  CALL FIELD_IO_ELEMENTS_EXPORT(fields,filename,method,err,error,*999)
                  IF(controlLoopOutputType >= CONTROL_LOOP_PROGRESS_OUTPUT) THEN
                    CALL WriteString(GENERAL_OUTPUT_TYPE,filename,err,error,*999)
                    CALL WriteString(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
                  ENDIF
                ENDIF
              ENDIF
            ENDIF

!!Subiteration intermediate solutions / iterates output:
            !IF(controlLoop%parentLoop%loopType==CONTROL_WHILE_LOOP_TYPE) THEN  !subiteration exists
            !  IF(currentIteration<10) THEN
            !    IF(subIterationNumber<10) THEN
            !      WRITE(outputFile,'("T_00",I0,"_SB_0",I0,"_C",I0)') currentIteration,subIterationNumber,equationsSetIdx
            !    ELSE IF(subIterationNumber<100) THEN
            !      WRITE(outputFile,'("T_00",I0,"_SB_",I0,"_C",I0)') currentIteration,subIterationNumber,equationsSetIdx
            !    END IF
            !    FILE=outputFile
            !    method="FORTRAN"
            !    exportField=.TRUE.
            !    IF(exportField) THEN
            !      CALL WriteString(GENERAL_OUTPUT_TYPE,"Darcy export subiterates ...",err,error,*999)
            !      CALL FLUID_MECHANICS_IO_WRITE_CMGUI(equationsSet%region,equationsSet%globalNumber,file,err,error,*999)
            !      CALL WriteString(GENERAL_OUTPUT_TYPE,outputFile,err,error,*999)
            !    ENDIF
            !  ENDIF
            !ENDIF

          ELSE !for single compartment (i.e. standary Darcy flow) equations sets
            !Find the time loop
            !If coupled with finite elasticity and using subiterations, get the while loop iteration number
            NULLIFY(parentLoop)
            CALL ControlLoop_ParentLoopGet(controlLoop,parentLoop,err,error,*999)
            CALL ControlLoop_TypeGet(parentLoop,parentLoopType,err,error,*999)
            IF(parentLoopType==CONTROL_WHILE_LOOP_TYPE) &
              & CALL ControlLoop_WhileInformationGet(parentLoop,subIterationNumber,maximumNumberOfIterations,absoluteTolerance, &
              & relativeTolerance,continueLoop,err,error,*999)            

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
                NULLIFY(region)
                CALL EquationsSet_RegionGet(equationsSet,region,err,error,*999)
                NULLIFY(fields)
                CALL Region_FieldsGet(region,fields,err,error,*999)              
                filename="./output/"//"MainTime_"//TRIM(NumberToVString(currentIteration,"*",err,error))
                method="FORTRAN"
                IF(MOD(currentIteration,outputIteration)==0)  THEN
                  CALL ControlLoop_OutputTypeGet(controlLoop,controlLoopOutputType,err,error,*999)
                  IF(controlLoopOutputType >= CONTROL_LOOP_PROGRESS_OUTPUT) THEN
                    CALL WriteString(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
                    CALL WriteString(GENERAL_OUTPUT_TYPE,"Now export fields... ",err,error,*999)
                  ENDIF
                  CALL FIELD_IO_NODES_EXPORT(Fields,filename,method,err,error,*999)
                  CALL FIELD_IO_ELEMENTS_EXPORT(Fields,filename,method,err,error,*999)
                  IF(controlLoopOutputType >= CONTROL_LOOP_PROGRESS_OUTPUT) THEN
                    CALL WriteString(GENERAL_OUTPUT_TYPE,filename,err,error,*999)
                    CALL WriteString(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
                  ENDIF
                ENDIF
              ENDIF
            ENDIF

!!Subiteration intermediate solutions / iterates output:
            !IF(controlLoop%parentLoop%loopType==CONTROL_WHILE_LOOP_TYPE) THEN  !subiteration exists
            !  IF(currentIteration<10) THEN
            !    IF(subIterationNumber<10) THEN
            !      WRITE(outputFile,'("T_00",I0,"_SUB_000",I0)') currentIteration,subIterationNumber
            !    ELSE IF(subIterationNumber<100) THEN
            !      WRITE(outputFile,'("T_00",I0,"_SUB_00",I0)') currentIteration,subIterationNumber
            !    END IF
            !    FILE=outputFile
            !    method="FORTRAN"
            !    exportField=.TRUE.
            !    IF(exportField) THEN
            !      CALL WriteString(GENERAL_OUTPUT_TYPE,"Darcy export subiterates ...",err,error,*999)
            !      CALL FLUID_MECHANICS_IO_WRITE_CMGUI(equationsSet%region,equationsSet%lobalNumber,file,err,error,*999)
            !      CALL WriteString(GENERAL_OUTPUT_TYPE,outputFile,err,error,*999)
            !    ENDIF
            !  ENDIF
            !ENDIF

          ENDIF
        ENDIF
      ENDDO !equationsSetIdx
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a Darcy fluid type of a fluid mechanics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Darcy_PostSolveOutputData")
    RETURN
999 ERRORSEXITS("Darcy_PostSolveOutputData",err,error)
    RETURN 1
    
  END SUBROUTINE Darcy_PostSolveOutputData

  !
  !================================================================================================================================
  !

  !>Calculates the analytic solution and sets the boundary conditions for an analytic problem.
  SUBROUTINE Darcy_BoundaryConditionsAnalyticCalculate(equationsSet,boundaryConditions,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: analyticFunctionType,boundCount,componentIdx,derivativeIdx,dimensionIdx,elementIdx,elementNode, &
      & elementNodeIdx,globalDerivativeIndex,i,j,k,localDOFIdx,maxElementParameters,nodeIdx,numberOfComponents,numberOfDimensions, &
      & numberOfElementNodes,numberOfNodes,numberOfNodeDerivatives,numberOfNodesXic(3),numberOfVariables,tempLocalDOFIdx, &
      & tempLocalNodeNumber,tempNodeNumber,variableIdx,variableType,velocityDOFCheck
    REAL(DP) ::  arg(3),boundaryTolerance,boundaryX(3,2),currentTime,fact,L,permeabilityOverViscosityParameter, &
      & tCoordinates(20,3),VALUE,x(3),xiCoordinates(3)
    REAL(DP), POINTER :: geometricParameters(:)
    LOGICAL :: boundaryNode
    TYPE(BasisType), POINTER :: basis
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainElementsType), POINTER :: domainElements
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(EquationsSetAnalyticType), POINTER :: equationsAnalytic
    TYPE(FieldType), POINTER :: dependentField,geometricField
    TYPE(FieldVariableType), POINTER :: dependentVariable,fieldVariable,geometricVariable
    TYPE(FieldInterpolatedPointPtrType), POINTER :: interpolatedPoint(:)
    TYPE(FieldInterpolationParametersPtrType), POINTER :: interpolationParameters(:)
    TYPE(VARYING_STRING) :: localError
    !Temp variables
 
    ENTERS("Darcy_BoundaryConditionsAnalyticCalculate",err,error,*999)

    boundCount=0

    permeabilityOverViscosityParameter = 1.0_DP  !temporarily hard-coded: Should rather be determined by interpolating materials field
    
    L=10.0_DP
    xiCoordinates(3)=0.0_DP
    boundaryTolerance=0.000000001_DP
    boundaryX=0.0_DP
    tCoordinates=0.0_DP

    numberOfElementNodes=0
    tempLocalNodeNumber=0
    tempLocalDOFIdx=0
    tempNodeNumber=0
    velocityDOFCheck=0

    IF(.NOT.ASSOCIATED(boundaryConditions)) CALL FlagError("Boundary conditions is not associated.",err,error,*999)
    NULLIFY(equationsAnalytic)
    CALL EquationsSet_AnalyticGet(equationsSet,equationsAnalytic,err,error,*999)
    CALL EquationsSet_AnalyticFunctionTypeGet(equationsSet,analyticFunctionType,err,error,*999)
    SELECT CASE(analyticFunctionType)
    CASE(EQUATIONS_SET_INCOMP_ELAST_DARCY_ANALYTIC_DARCY)
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      NULLIFY(dependentField)
      CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
      NULLIFY(geometricVariable)
      CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,geometricVariable,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(geometricVariable,numberOfDimensions,err,error,*999)
      NULLIFY(geometricParameters)
      CALL FieldVariable_ParameterSetDataGet(geometricVariable,FIELD_VALUES_SET_TYPE,geometricParameters,err,error,*999)
      CALL Field_NumberOfVariablesGet(dependentField,numberOfVariables,err,error,*999)
      CALL EquationsSet_AnalyticTimeGet(equationsSet,currentTime,err,error,*999)
      DO variableIdx=3,numberOfVariables
        NULLIFY(fieldVariable)
        CALL Field_VariableIndexGet(dependentField,variableIdx,fieldVariable,variableType,err,error,*999)
        CALL FieldVariable_ParameterSetEnsureCreated(fieldVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
        CALL FieldVariable_NumberOfComponentsGet(fieldVariable,numberOfComponents,err,error,*999)
        DO componentIdx=1,numberOfComponents
          CALL FieldVariable_ComponentInterpolationCheck(fieldVariable,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
          NULLIFY(domain)
          CALL FieldVariable_ComponentDomainGet(fieldVariable,componentIdx,domain,err,error,*999)
          NULLIFY(domainTopology)
          CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
          NULLIFY(domainNodes)
          CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
          !Loop over the local nodes excluding the ghosts.
          CALL DomainNodes_NumberOfNodesGet(domainNodes,numberOfNodes,err,error,*999)
          DO nodeIdx=1,numberOfNodes
!!TODO \todo We should interpolate the geometric field here and the node position.
            DO dimensionIdx=1,numberOfDimensions
              !Default to version 1 of each node derivative
              CALL FieldVariable_LocalNodeDOFGet(geometricVariable,1,1,nodeIdx,dimensionIdx,localDOFIdx,err,error,*999)
              x(dimensionIdx)=geometricParameters(localDOFIdx)
            ENDDO !dimensionIdx
            CALL DomainNodes_NodeBoundaryNodeGet(domainNodes,nodeIdx,boundaryNode,err,error,*999)
            !Loop over the derivatives
            CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,nodeIdx,numberOfNodeDerivatives,err,error,*999)
            DO derivativeIdx=1,numberOfNodeDerivatives
              CALL DomainNodes_DerivativeGlobalIndexGet(domainNodes,derivativeIdx,nodeIdx,globalDerivativeIndex,err,error,*999)
              !CALL DIFFUSION_EQUATION_ANALYTIC_FUNCTIONS(VALUE,X,currentTime,variableType, &
              !  & globalDerivativeIndex,analyticFunctionType,err,error,*999)
!!!!!!!!!!!!NEED TO SET APPROPRIATE VALUE DEPENDING ON WHETHER IT IS A VELOCITY COMPONENT OR THE MASS INCREASE COMPONENT
              VALUE=0.0_DP
              !Default to version 1 of each node derivative
              CALL FieldVariable_LocalNodeDOFGet(fieldVariable,1,derivativeIdx,nodeIdx,componentIdx,localDOFIdx,err,error,*999)
              CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,localDOFIdx, &
                & VALUE,err,error,*999)
              IF(variableType==FIELD_V_VARIABLE_TYPE) THEN
                IF(boundaryNode) THEN
                  !If we are a boundary node then set the analytic value on the boundary
                  CALL BoundaryConditions_SetLocalDOF(boundaryConditions,dependentVariable,localDOFIdx, &
                    & BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
                ELSE
                  CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_VALUES_SET_TYPE,localDOFIdx,VALUE, &
                    & err,error,*999)
                ENDIF
              ENDIF
            ENDDO !derivativeIdx
          ENDDO !nodeIdx
        ENDDO !componentIdx
        CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
        CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
      ENDDO !variableIdx
      CALL FieldVariable_ParameterSetDataRestore(geometricVariable,FIELD_VALUES_SET_TYPE,geometricParameters,err,error,*999)
    CASE DEFAULT
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      NULLIFY(geometricVariable)
      CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,geometricVariable,err,error,*999)
      NULLIFY(dependentField)
      CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
      NULLIFY(interpolationParameters)
      CALL Field_InterpolationParametersInitialise(geometricField,interpolationParameters,err,error,*999)
      NULLIFY(interpolatedPoint)
      CALL Field_InterpolatedPointsInitialise(interpolationParameters,interpolatedPoint,err,error,*999)

      CALL FieldVariable_NumberOfComponentsGet(geometricVariable,numberOfDimensions,err,error,*999)

      IF(numberOfDimensions==2) THEN
        boundaryX(1,1)=0.0_DP
        boundaryX(1,2)=10.0_DP
        boundaryX(2,1)=0.0_DP
        boundaryX(2,2)=10.0_DP
      ELSE IF(numberOfDimensions==3) THEN
        boundaryX(1,1)=-5.0_DP
        boundaryX(1,2)=5.0_DP
        boundaryX(2,1)=-5.0_DP
        boundaryX(2,2)=5.0_DP
        boundaryX(3,1)=-5.0_DP
        boundaryX(3,2)=5.0_DP
      ENDIF

      NULLIFY(geometricParameters)
      CALL FieldVariable_ParameterSetDataGet(geometricVariable,FIELD_VALUES_SET_TYPE,geometricParameters,err,error,*999)
      CALL Field_NumberOfVariablesGet(dependentField,numberOfVariables,err,error,*999)
      DO variableIdx=1,numberOfVariables
        NULLIFY(fieldVariable)
        CALL Field_VariableIndexGet(dependentField,variableIdx,fieldVariable,variableType,err,error,*999)
        CALL FieldVariable_ParameterSetEnsureCreated(fieldVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
        CALL FieldVariable_NumberOfCompnentsGet(fieldVariable,numberOfComponents,err,error,*999)
        DO componentIdx=1,numberOfComponents
          boundCount=0
          CALL FieldVariable_ComponentInterpolationCheck(fieldVariable,componentIdx,FIELD_NODE_BASED_INTERPOLATION, &
            & err,error,*999)
          NULLIFY(domain)
          CALL FieldVariable_ComponentDomainGet(fieldVariable,componentIdx,domain,err,error,*999)
          NULLIFY(domainTopology)
          CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
          NULLIFY(domainNodes)
          CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
          NULLIFY(domainElements)
          CALL DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*999)
          CALL DomainElements_MaxElementParametersGet(domainElements,maxElementParameters,err,error,*999)
          !Loop over the local nodes excluding the ghosts.
          CALL DomainNodes_NumberOfNodesGet(domainNodes,numberOfNodes,err,error,*999)
          DO nodeIdx=1,numberOfNodes
            CALL DomainNodes_NodeSurroundingElementGet(domainNodes,1,nodeIdx,elementIdx,err,error,*999)
            NULLIFY(basis)
            CALL DomainElements_ElementBasisGet(domainElements,elementIdx,basis,err,error,*999)
            CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementIdx, &
              & interpolationParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)

            !DO I=1,DOMAIN%topology%elements%maximumNumberOfElementParameters
            !  IF(DOMAIN%topology%elements%elements(elementIdx)%elementNodes(I)=nodeIdx THEN

            elementNodeIdx=0
            xiCoordinates=0.0_DP
            numberOfNodesXiC=1
            CALL Basis_NumberOfNodesXiCGet(basis,numberOfNodesXiC,err,error,*999)

            IF(maxElementParameters==4.AND.numberOfDimensions==2 .OR. &
              & maxElementParameters==9.OR. &
              & maxElementParameters==16.OR. &
              & maxElementParameters==8.OR. &
              & maxElementParameters==27.OR. &
              & maxElementParameters==64) THEN

              DO K=1,numberOfNodesXic(3)
                DO J=1,numberOfNodesXic(2)
                  DO I=1,numberOfNodesXic(1)
                    elementNodeIdx=elementNodeIdx+1
                    CALL DomainElement_ElementNodeGet(domainElements,elementNodeIdx,elementIdx,elementNode,err,error,*999)
                    IF(elementNode==nodeIdx) EXIT
                    xiCoordinates(1)=xiCoordinates(1)+(1.0_DP/(numberOfNodesXic(1)-1))
                  ENDDO
                  IF(elementNode==nodeIdx) EXIT
                  xiCoordinates(1)=0.0_DP
                  xiCoordinates(2)=xiCoordinates(2)+(1.0_DP/(numberOfNodesXic(2)-1))
                ENDDO
                IF(elementNode==nodeIdx) EXIT
                xiCoordinates(1)=0.0_DP
                xiCoordinates(2)=0.0_DP
                IF(numberOfNodesXic(3)/=1) THEN
                  xiCoordinates(3)=xiCoordinates(3)+(1.0_DP/(numberOfNodesXic(3)-1))
                ENDIF
              ENDDO
              CALL Field_InterpolateXi(NO_PART_DERIV,xiCoordinates,interpolatedPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            ELSE
              IF(maxElementParameters==3) THEN
                tCoordinates(1,1:2)=[0.0_DP,1.0_DP]
                tCoordinates(2,1:2)=[1.0_DP,0.0_DP]
                tCoordinates(3,1:2)=[1.0_DP,1.0_DP]
              ELSE IF(maxElementParameters==6) THEN
                tCoordinates(1,1:2)=[0.0_DP,1.0_DP]
                tCoordinates(2,1:2)=[1.0_DP,0.0_DP]
                tCoordinates(3,1:2)=[1.0_DP,1.0_DP]
                tCoordinates(4,1:2)=[0.5_DP,0.5_DP]
                tCoordinates(5,1:2)=[1.0_DP,0.5_DP]
                tCoordinates(6,1:2)=[0.5_DP,1.0_DP]
              ELSE IF(maxElementParameters==10.AND.numberOfDimensions==2) THEN
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
              ELSE IF(maxElementParameters==4) THEN
                tCoordinates(1,1:3)=[0.0_DP,1.0_DP,1.0_DP]
                tCoordinates(2,1:3)=[1.0_DP,0.0_DP,1.0_DP]
                tCoordinates(3,1:3)=[1.0_DP,1.0_DP,0.0_DP]
                tCoordinates(4,1:3)=[1.0_DP,1.0_DP,1.0_DP]
              ELSE IF(maxElementParameters==10.AND.numberOfDimensions==3) THEN
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
              ELSE IF(maxElementParameters==20) THEN
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

              DO K=1,maxElementParameters
                CALL DomainElement_ElementNodeGet(domainElements,K,elementIdx,elementNode,err,error,*999)
                IF(elementNode==nodeIdx) EXIT
              ENDDO !K

              IF(numberOfDimensions==2) THEN
                CALL Field_InterpolateXi(NO_PART_DERIV,tCoordinates(K,1:2),interpolatedPoint(FIELD_U_VARIABLE_TYPE)%ptr, &
                  & err,error,*999)
              ELSE IF(numberOfDimensions==3) THEN
                CALL Field_InterpolateXi(NO_PART_DERIV,tCoordinates(K,1:3),interpolatedPoint(FIELD_U_VARIABLE_TYPE)%ptr, &
                  & err,error,*999)
              ENDIF
            ENDIF

            X=0.0_DP
            DO dimensionIdx=1,numberOfDimensions
              x(dimensionIdx)=interpolatedPoint(FIELD_U_VARIABLE_TYPE)%ptr%values(dimensionIdx,1)
            ENDDO !dimensionIdx

            !Loop over the derivatives
            CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,nodeIdx,numberOfNodeDerivatives,err,error,*999)
            DO derivativeIdx=1,numberOfNodeDerivatives
              CALL DomainNodes_DerivativeGlobalIndexGet(domainNodes,derivativeIdx,nodeIdx,globalDerivativeIndex,err,error,*999)
              SELECT CASE(analyticFunctionType)
              CASE(EQUATIONS_SET_DARCY_EQUATION_TWO_DIM_1)
                IF(numberOfDimensions==2.AND.numberOfComponents==3) THEN
                  !POLYNOM
                  SELECT CASE(variableType)
                  CASE(FIELD_U_VARIABLE_TYPE)
                    SELECT CASE(globalDerivativeIndex)
                    CASE(NO_GLOBAL_DERIV)
                      fact = permeabilityOverViscosityParameter
                      IF(componentIdx==1) THEN
                        !calculate u
                        VALUE = - fact * ( 2.0_DP*x(1) + 2.0_DP*x(2) )
                      ELSE IF(componentIdx==2) THEN
                        !calculate v
                        VALUE = - fact * ( 2.0_DP*x(1) - 2.0_DP*x(2) )
                      ELSE IF(componentIdx==3) THEN
                        !calculate p
                        VALUE = x(1)**2.0_DP + 2.0_DP*x(1)*x(2) - x(2)**2.0_DP
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
                      localError="The global derivative index of "// &
                        & TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))//" is invalid."
                      CALL FlagError(localError,err,error,*999)
                    END SELECT
                  CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                    SELECT CASE(globalDerivativeIndex)
                    CASE(NO_GLOBAL_DERIV)
                      VALUE= 0.0_DP
                    CASE(GLOBAL_DERIV_S1)
                      CALL FlagError("Not implemented.",err,error,*999)
                    CASE(GLOBAL_DERIV_S2)
                      CALL FlagError("Not implemented.",err,error,*999)
                    CASE(GLOBAL_DERIV_S1_S2)
                      CALL FlagError("Not implemented.",err,error,*999)
                    CASE DEFAULT
                      localError="The global derivative index of "// &
                        & TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))//" is invalid."
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


              CASE(EQUATIONS_SET_DARCY_EQUATION_TWO_DIM_2)
                IF(numberOfDimensions==2.AND.numberOfComponents==3) THEN
                  !EXPONENTIAL
                  SELECT CASE(variableType)
                  CASE(FIELD_U_VARIABLE_TYPE)
                    SELECT CASE(globalDerivativeIndex)
                    CASE(NO_GLOBAL_DERIV)
                      fact   = permeabilityOverViscosityParameter / L
                      arg(1) = x(1) / L
                      arg(2) = x(2) / L
                      IF(componentIdx==1) THEN
                        !calculate u
                        VALUE = - fact * EXP( arg(1) ) * EXP( arg(2) )
                      ELSE IF(componentIdx==2) THEN
                        !calculate v
                        VALUE = - fact * EXP( arg(1) ) * EXP( arg(2) )
                      ELSE IF(componentIdx==3) THEN
                        !calculate p
                        VALUE =          EXP( arg(1) ) * EXP( arg(2) )
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
                      localError="The global derivative index of "// &
                        & TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))//" is invalid."
                      CALL FlagError(localError,err,error,*999)
                    END SELECT
                  CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                    SELECT CASE(globalDerivativeIndex)
                    CASE(NO_GLOBAL_DERIV)
                      IF(componentIdx==1) THEN
                        !calculate u
                        VALUE= 0.0_DP
                      ELSE IF(componentIdx==2) THEN
                        !calculate v
                        VALUE= 0.0_DP
                      ELSE IF(componentIdx==3) THEN
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
                      localError="The global derivative index of "// &
                        & TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))//" is invalid."
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
              CASE(EQUATIONS_SET_DARCY_EQUATION_TWO_DIM_3)
                IF(numberOfDimensions==2.AND.numberOfComponents==3) THEN
                  !SINUS/COSINUS
                  SELECT CASE(variableType)
                  CASE(FIELD_U_VARIABLE_TYPE)
                    SELECT CASE(globalDerivativeIndex)
                    CASE(NO_GLOBAL_DERIV)
                      fact = 2.0_DP * PI * permeabilityOverViscosityParameter / L
                      arg(1) = 2.0_DP * PI * x(1) / L
                      arg(2) = 2.0_DP * PI * x(2) / L
                      IF(componentIdx==1) THEN
                        !calculate u
                        VALUE = - fact * COS( arg(1) ) * SIN( arg(2) )
                      ELSE IF(componentIdx==2) THEN
                        !calculate v
                        VALUE = - fact * SIN( arg(1) ) * COS( arg(2) )
                      ELSE IF(componentIdx==3) THEN
                        !calculate p
                        VALUE =          SIN( arg(1) ) * SIN( arg(2) )
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
                      localError="The global derivative index of "// &
                        & TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))//" is invalid."
                      CALL FlagError(localError,err,error,*999)
                    END SELECT
                  CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                    SELECT CASE(globalDerivativeIndex)
                    CASE(NO_GLOBAL_DERIV)
                      IF(componentIdx==1) THEN
                        !calculate u
                        VALUE=0.0_DP
                      ELSE IF(componentIdx==2) THEN
                        !calculate v
                        VALUE=0.0_DP
                      ELSE IF(componentIdx==3) THEN
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
                      localError="The global derivative index of "// &
                        & TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))//" is invalid."
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

              CASE(EQUATIONS_SET_DARCY_EQUATION_THREE_DIM_1)
                IF(numberOfDimensions==3.AND.fieldVariable%numberOfComponents==4) THEN
                  !POLYNOM
                  SELECT CASE(variableType)
                  CASE(FIELD_U_VARIABLE_TYPE)
                    SELECT CASE(globalDerivativeIndex)
                    CASE(NO_GLOBAL_DERIV)
                      fact = permeabilityOverViscosityParameter
                      IF(componentIdx==1) THEN
                        !calculate u
                        VALUE = - fact * ( 2.0_DP*x(1) + 2.0_DP*x(2) + x(3) )
                      ELSE IF(componentIdx==2) THEN
                        !calculate v
                        VALUE = - fact * ( 2.0_DP*x(1) - 2.0_DP*x(2) + x(3) )
                      ELSE IF(componentIdx==3) THEN
                        !calculate w
                        VALUE = - fact * ( 3.0_DP + x(1) + x(2) )
                      ELSE IF(componentIdx==4) THEN
                        !calculate p
                        VALUE = x(1)**2.0_DP + 2.0_DP*x(1)*x(2) - x(2)**2.0_DP + &
                          & 3.0_DP*x(3) + x(3)*x(1) + x(3)*x(2)
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
                        & domainNodes%NODES(nodeIdx)%DERIVATIVES(derivativeIdx)%globalDerivativeIndex,"*", &
                        & err,error))//" is invalid."
                      CALL FlagError(localError,err,error,*999)
                    END SELECT
                  CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                    SELECT CASE(globalDerivativeIndex)
                    CASE(NO_GLOBAL_DERIV)
                      VALUE=0.0_DP
                    CASE(GLOBAL_DERIV_S1)
                      CALL FlagError("Not implemented.",err,error,*999)
                    CASE(GLOBAL_DERIV_S2)
                      CALL FlagError("Not implemented.",err,error,*999)
                    CASE(GLOBAL_DERIV_S1_S2)
                      CALL FlagError("Not implemented.",err,error,*999)
                    CASE DEFAULT
                      localError="The global derivative index of "// &
                        & TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))//" is invalid."
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
              CASE(EQUATIONS_SET_DARCY_EQUATION_THREE_DIM_2)
                IF(numberOfDimensions==3.AND.fieldVariable%numberOfComponents==4) THEN
                  !EXPONENTIAL
                  SELECT CASE(variableType)
                  CASE(FIELD_U_VARIABLE_TYPE)
                    SELECT CASE(globalDerivativeIndex)
                    CASE(NO_GLOBAL_DERIV)
                      fact   = permeabilityOverViscosityParameter / L
                      arg(1) = x(1) / L
                      arg(2) = x(2) / L
                      arg(3) = x(3) / L
                      IF(componentIdx==1) THEN
                        !calculate u
                        VALUE = - fact * EXP( arg(1) ) * EXP( arg(2) ) * EXP( arg(3) )
                      ELSE IF(componentIdx==2) THEN
                        !calculate v
                        VALUE = - fact * EXP( arg(1) ) * EXP( arg(2) ) * EXP( arg(3) )
                      ELSE IF(componentIdx==3) THEN
                        !calculate w
                        VALUE = - fact * EXP( arg(1) ) * EXP( arg(2) ) * EXP( arg(3) )
                      ELSE IF(componentIdx==4) THEN
                        !calculate p
                        VALUE =          EXP( arg(1) ) * EXP( arg(2) ) * EXP( arg(3) )
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
                      localError="The global derivative index of "// &
                        & TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))//" is invalid."
                      CALL FlagError(localError,err,error,*999)
                    END SELECT
                  CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                    SELECT CASE(globalDerivativeIndex)
                    CASE(NO_GLOBAL_DERIV)
                      IF(componentIdx==1) THEN
                        !calculate u
                        VALUE=0.0_DP
                      ELSE IF(componentIdx==2) THEN
                        !calculate v
                        VALUE=0.0_DP
                      ELSE IF(componentIdx==3) THEN
                        !calculate w
                        VALUE=0.0_DP
                      ELSE IF(componentIdx==4) THEN
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
                      localError="The global derivative index of "// &
                        & TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))//" is invalid."
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

              CASE(EQUATIONS_SET_DARCY_EQUATION_THREE_DIM_3)
                IF(numberOfDimensions==3.AND.fieldVariable%numberOfComponents==4) THEN
                  !SINE/COSINE
                  SELECT CASE(variableType)
                  CASE(FIELD_U_VARIABLE_TYPE)
                    SELECT CASE(globalDerivativeIndex)
                    CASE(NO_GLOBAL_DERIV)
                      fact = 2.0_DP * PI * permeabilityOverViscosityParameter / L
                      arg(1) = 2.0_DP * PI * x(1) / L
                      arg(2) = 2.0_DP * PI * x(2) / L
                      arg(3) = 2.0_DP * PI * x(3) / L
                      IF(componentIdx==1) THEN
                        !calculate u
                        VALUE = - fact * COS( arg(1) ) * SIN( arg(2) )  * SIN( arg(3) )
                      ELSE IF(componentIdx==2) THEN
                        !calculate v
                        VALUE = - fact * SIN( arg(1) ) * COS( arg(2) )  * SIN( arg(3) )
                      ELSE IF(componentIdx==3) THEN
                        !calculate w
                        VALUE = - fact * SIN( arg(1) ) * SIN( arg(2) )  * COS( arg(3) )
                      ELSE IF(componentIdx==4) THEN
                        !calculate p
                        VALUE =          SIN( arg(1) ) * SIN( arg(2) )  * SIN( arg(3) )
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
                      localError="The global derivative index of "// &
                        & TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))//" is invalid."
                      CALL FlagError(localError,err,error,*999)
                    END SELECT
                  CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                    SELECT CASE(globalDerivativeIndex)
                    CASE(NO_GLOBAL_DERIV)
                      IF(componentIdx==1) THEN
                        !calculate u
                        VALUE=0.0_DP
                      ELSE IF(componentIdx==2) THEN
                        !calculate v
                        VALUE=0.0_DP
                      ELSE IF(componentIdx==3) THEN
                        !calculate w
                        VALUE=0.0_DP
                      ELSE IF(componentIdx==4) THEN
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
                      localError="The global derivative index of "// &
                        & TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))//" is invalid."
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
              !Default to version 1 of each node derivative
              CALL FieldVariable_LocalNodeDOFGet(fieldVariable,1,derivativeIdx,nodeIdx,componentIdx,localDOFIdx,err,error,*999)
              CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,localDOFIdx, &
                & VALUE,err,error,*999)
              IF(variableType==FIELD_U_VARIABLE_TYPE) THEN                  
                ! ! !                                 IF(domainNodes%NODES(nodeIdx)%boundaryNode) THEN
                ! ! !                                   !If we are a boundary node then set the analytic value on the boundary
                ! ! !                                   IF(componentIdx<=numberOfDimensions) THEN
                ! ! !                                     CALL BoundaryConditions_SetLocalDOF(boundaryConditions,variableType,localDOFIdx, &
                ! ! !                                       & BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
                ! ! !                                   boundCount=boundCount+1
                ! ! !                                   ENDIF
                ! ! !                                 ELSE
                ! ! !                                   IF(componentIdx<=numberOfDimensions) THEN
                ! ! !                                     CALL Field_ParameterSetUpdateLocalDOF(dependentField,variableType, &
                ! ! !                                       & FIELD_VALUES_SET_TYPE,localDOFIdx,VALUE,err,error,*999)
                ! ! !                                   ENDIF
                ! ! !                                 ENDIF                                                      

                !If we are a boundary node then set the analytic value on the boundary
                IF(numberOfDimensions==2) THEN
                  IF(ABS(x(1)-boundaryX(1,1))<boundaryTolerance.OR. &
                    & ABS(x(1)-boundaryX(1,2))<boundaryTolerance.OR. &
                    & ABS(x(2)-boundaryX(2,1))<boundaryTolerance.OR. &
                    & ABS(x(2)-boundaryX(2,2))<boundaryTolerance) THEN
                    IF(componentIdx<=numberOfDimensions) THEN
                      CALL BoundaryConditions_SetLocalDOF(boundaryConditions,dependentVariable,localDOFIdx, &
                        & BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
                      boundCount=boundCount+1
                      !Apply boundary conditions check for pressure nodes
                    ELSE IF(componentIdx>numberOfDimensions) THEN
                      IF(maxElementParameters==3.OR. &
                        & maxElementParameters==6.OR. &
                        & maxElementParameters==10) THEN                          
                        IF(ABS(x(1)-boundaryX(1,1))<boundaryTolerance.AND. &
                          & ABS(x(2)-boundaryX(2,1))<boundaryTolerance.OR. &
                          & ABS(x(1)-boundaryX(1,1))<boundaryTolerance.AND.&
                          & ABS(x(2)-boundaryX(2,2))<boundaryTolerance.OR. &
                          & ABS(x(1)-boundaryX(1,2))<boundaryTolerance.AND.&
                          & ABS(x(2)-boundaryX(2,1))<boundaryTolerance.OR. &
                          & ABS(x(1)-boundaryX(1,2))<boundaryTolerance.AND.&
                          & ABS(x(2)-boundaryX(2,2))<boundaryTolerance) &
                          & THEN
                          CALL BoundaryConditions_SetLocalDOF(boundaryConditions,dependentVariable,localDOFIdx, &
                            & BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
                          boundCount=boundCount+1
                        ENDIF
                      ENDIF
                    ENDIF
                  ELSE
                    IF(componentIdx<=numberOfDimensions) THEN
                      CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_VALUES_SET_TYPE,localDOFIdx, &
                        & VALUE,err,error,*999)
                    ENDIF
                  ENDIF
                ELSE IF(numberOfDimensions==3) THEN
                  IF(ABS(x(1)-boundaryX(1,1))<boundaryTolerance.OR. &
                    & ABS(x(1)-boundaryX(1,2))<boundaryTolerance.OR. &
                    & ABS(x(2)-boundaryX(2,1))<boundaryTolerance.OR. &
                    & ABS(x(2)-boundaryX(2,2))<boundaryTolerance.OR. &
                    & ABS(x(3)-boundaryX(3,1))<boundaryTolerance.OR. &
                    & ABS(x(3)-boundaryX(3,2))<boundaryTolerance) THEN
                    IF(componentIdx<=numberOfDimensions) THEN
                      CALL BoundaryConditions_SetLocalDOF(boundaryConditions,dependentVariable,localDOFIdx, &
                        & BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
                      boundCount=boundCount+1
                      !Apply boundary conditions check for pressure nodes
                    ELSE IF(componentIdx>numberOfDimensions) THEN
                      IF(maxElementParameters==4.OR. &
                        & maxElementParameters==10.OR. &
                        & maxElementParameters==20) THEN
                        IF(ABS(x(1)-boundaryX(1,1))<boundaryTolerance.AND. &
                          & ABS(x(2)-boundaryX(2,1))<boundaryTolerance.AND. &
                          & ABS(x(3)-boundaryX(3,1))<boundaryTolerance.OR. &
                          & ABS(x(1)-boundaryX(1,1))<boundaryTolerance.AND. &
                          & ABS(x(2)-boundaryX(2,1))<boundaryTolerance.AND. &
                          & ABS(x(3)-boundaryX(3,2))<boundaryTolerance.OR. &
                          & ABS(x(1)-boundaryX(1,1))<boundaryTolerance.AND. &
                          & ABS(x(2)-boundaryX(2,2))<boundaryTolerance.AND. &
                          & ABS(x(3)-boundaryX(3,1))<boundaryTolerance.OR. &
                          & ABS(x(1)-boundaryX(1,1))<boundaryTolerance.AND. &
                          & ABS(x(2)-boundaryX(2,2))<boundaryTolerance.AND. &
                          & ABS(x(3)-boundaryX(3,2))<boundaryTolerance.OR. &
                          & ABS(x(1)-boundaryX(1,2))<boundaryTolerance.AND. &
                          & ABS(x(2)-boundaryX(2,1))<boundaryTolerance.AND. &
                          & ABS(x(3)-boundaryX(3,1))<boundaryTolerance.OR. &
                          & ABS(x(1)-boundaryX(1,2))<boundaryTolerance.AND. &
                          & ABS(x(2)-boundaryX(2,1))<boundaryTolerance.AND. &
                          & ABS(x(3)-boundaryX(3,2))<boundaryTolerance.OR. &
                          & ABS(x(1)-boundaryX(1,2))<boundaryTolerance.AND. &
                          & ABS(x(2)-boundaryX(2,2))<boundaryTolerance.AND. &
                          & ABS(x(3)-boundaryX(3,1))<boundaryTolerance.OR. &
                          & ABS(x(1)-boundaryX(1,2))<boundaryTolerance.AND. &
                          & ABS(x(2)-boundaryX(2,2))<boundaryTolerance.AND. &
                          & ABS(x(3)-boundaryX(3,2))<boundaryTolerance) THEN
                          CALL BoundaryConditions_SetLocalDOF(boundaryConditions,dependentVariable,localDOFIdx, &
                            & BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
                          boundCount=boundCount+1
                        ENDIF
                      ENDIF
                    ENDIF
                  ELSE
                    IF(componentIdx<=numberOfDimensions) THEN
                      CALL FieldVariable_ParameterSetUpdateLocalDOF(fieldVariable,FIELD_VALUES_SET_TYPE,localDOFIdx, &
                        & VALUE,err,error,*999)
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
            ENDDO !derivativeIdx
          ENDDO !nodeIdx
        ENDDO !componentIdx
        CALL FieldVariable_ParameterSetUpdateStart(fieldVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
        CALL FieldVariable_ParameterSetUpdateFinish(fieldVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
        CALL FieldVariable_ParameterSetUpdateStart(fieldVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
        CALL FieldVariable_ParameterSetUpdateFinish(fieldVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
      ENDDO !variableIdx
      CALL FieldVariable_ParameterSetDataRestore(geometricVariable,FIELD_VALUES_SET_TYPE,geometricParameters,err,error,*999)
    END SELECT

    EXITS("Darcy_BoundaryConditionsAnalyticCalculate")
    RETURN
999 ERRORSEXITS("Darcy_BoundaryConditionsAnalyticCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE Darcy_BoundaryConditionsAnalyticCalculate

  !
  !================================================================================================================================
  !

  !>Update geometric field for ALE Darcy problem
  SUBROUTINE Darcy_PreSolveGetSolidDisplacement(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: currentIteration,dofNumber,esSpecification(3),inputIteration,inputOption,inputType,loopIdx, &
      & numberOfDimensions,numberOfDofs,numberDOFsToPrint,outputIteration,pSpecification(3),solverOutputType
    REAL(DP) :: alpha,currentTime,startTime,stopTime,timeIncrement
    REAL(DP), POINTER :: dummyValues2(:),meshDisplacementValues(:),solutionValuesSolid(:)
    TYPE(ControlLoopType), POINTER :: controlLoop,controlTimeLoop,rootControlLoop,solidControlLoop
    TYPE(EquationsSetType), POINTER :: equationsSetFiniteElasticity, equationsSetDarcy
    TYPE(FieldType), POINTER :: dependentFieldFiniteElasticity, geometricFieldDarcy
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverType), POINTER :: solverFiniteElasticity, solverDarcy
    TYPE(SolversType), POINTER :: solvers
    TYPE(SolverEquationsType), POINTER :: solverEquationsFiniteElasticity, solverEquationsDarcy 
    TYPE(SolverMappingType), POINTER :: solverMappingFiniteElasticity, solverMappingDarcy
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Darcy_PreSolveGetSolidDisplacement",err,error,*999)

!--- \todo : Do we need for each case a Field_ParameterSetUpdateStart / FINISH on FIELD_MESH_DISPLACEMENT_SET_TYPE ?

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)

    CALL ControlLoop_CurrentTimeInformationGet(controlLoop,currentTime,timeIncrement,startTime,stopTime,currentIteration, &
      & outputIteration,inputIteration,err,error,*999)
    
    IF(diagnostics1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE, &
        & "*******************************************************************************************************", &
        & err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Current Time   = ",currentTime,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Time Increment = ",timeIncrement,err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE, &
        & "*******************************************************************************************************", &
        & err,error,*999)
    ENDIF

    NULLIFY(rootControlLoop)
    CALL Problem_RootControlLoopGet(problem,rootControlLoop,err,error,*999)
    CALL Solver_OutputTypeGet(solver,solverOutputType,err,error,*999)
    NULLIFY(solvers)
    CALL Solver_SolversGet(solver,solvers,err,error,*999)
    
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STANDARD_DARCY_SUBTYPE)
      ! do nothing
    CASE(PROBLEM_QUASISTATIC_DARCY_SUBTYPE)
      ! do nothing
    CASE(PROBLEM_TRANSIENT_DARCY_SUBTYPE)
      ! do nothing
    CASE(PROBLEM_ALE_DARCY_SUBTYPE)
      !--- Motion: specified
      NULLIFY(solverEquationsDarcy)
      CALL Solver_SolverEquationsGet(solver,solverEquationsDarcy,err,error,*999)
      NULLIFY(solverMappingDarcy)
      CALL SolverEquations_SolverMappingGet(solverEquationsDarcy,solverMappingDarcy,err,error,*999)
      NULLIFY(equationsSetDarcy)
      CALL SolverMapping_EquationsSetGet(solverMappingDarcy,1,equationsSetDarcy,err,error,*999)
      CALL EquationsSet_SpecificationGet(equationsSetDarcy,3,esSpecification,err,error,*999)
      SELECT CASE(esSpecification(3))
      CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE)
        ! do nothing
      CASE(EQUATIONS_SET_QUASISTATIC_DARCY_SUBTYPE)
        ! do nothing
      CASE(EQUATIONS_SET_ALE_DARCY_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
        & EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
        & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE)
        IF(solverOutputType>=SOLVER_PROGRESS_OUTPUT) &
          & CALL WriteString(GENERAL_OUTPUT_TYPE,"Darcy motion specified ... ",err,error,*999)
        NULLIFY(geometricFieldDarcy)
        CALL EquationsSet_GeometricFieldGet(equationsSetDarcy,geometricFieldDarcy,err,error,*999)
        alpha = 0.085_DP * SIN( 2.0_DP * PI * currentTime / 4.0_DP )        
        CALL Field_ParameterSetsCopy(geometricFieldDarcy,FIELD_U_VARIABLE_TYPE,FIELD_INITIAL_VALUES_SET_TYPE, &
          & FIELD_MESH_DISPLACEMENT_SET_TYPE,alpha,err,error,*999)
      CASE DEFAULT
        localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
          & " is not valid for a Darcy equation fluid type of a fluid mechanics problem class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_PGM_DARCY_SUBTYPE,PROBLEM_PGM_TRANSIENT_DARCY_SUBTYPE,PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE)
      !--- Motion: read in from a file
      NULLIFY(solverEquationsDarcy)
      CALL Solver_SolverEquationsGet(solver,solverEquationsDarcy,err,error,*999)
      NULLIFY(solverMappingDarcy)
      CALL SolverEquations_SolverMappingGet(solverEquationsDarcy,solverMappingDarcy,err,error,*999)
      NULLIFY(equationsSetDarcy)
      CALL SolverMapping_EquationsSetGet(solverMappingDarcy,1,equationsSetDarcy,err,error,*999)
      CALL EquationsSet_SpecificationGet(equationsSetDarcy,3,esSpecification,err,error,*999)
      NULLIFY(geometricFieldDarcy)
      CALL EquationsSet_GeometricFieldGet(equationsSetDarcy,geometricFieldDarcy,err,error,*999)
      IF(solverOutputType>=SOLVER_PROGRESS_OUTPUT) &
        & CALL WriteString(GENERAL_OUTPUT_TYPE,"Darcy motion read from a file ... ",err,error,*999)
      CALL Field_NumberOfComponentsGet(geometricFieldDarcy,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
      !Copy input to Darcy' geometric field
      inputType=42
      inputOption=2
      CALL Field_ParameterSetDataGet(geometricFieldDarcy,FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE, &
        & meshDisplacementValues,err,error,*999)
      CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,meshDisplacementValues, &
        & numberOfDimensions,inputType,inputOption,controlLoop%timeLoop%iterationNumber,1.0_DP, &
        & err,error,*999)
      CALL Field_ParameterSetUpdateStart(geometricFieldDarcy,FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE,err,error,*999)
      CALL Field_ParameterSetUpdateFinish(geometricFieldDarcy,FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE,err,error,*999)
      
      IF(diagnostics1) THEN
        numberDOFsToPrint = SIZE(meshDisplacementValues,1)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberDOFsToPrint,numberDOFsToPrint,numberDOFsToPrint,&
          & meshDisplacementValues,'(" meshDisplacementValues = ",4(X,E13.6))','4(4(X,E13.6))', &
          & err,error,*999)
      ENDIF
    CASE(PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE, &
      & PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
      !--- Motion: defined by fluid-solid interaction (thus read from solid's dependent field)
      !--- Get the dependent field of the finite elasticity equations
      IF(solverOutputType>=SOLVER_PROGRESS_OUTPUT) &
        & CALL WriteString(GENERAL_OUTPUT_TYPE,"Darcy motion read from solid's dependent field ... ",err,error,*999)
      SELECT CASE(pSpecification(3))
      CASE(PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE)
        CALL ControlLoop_Get(rootControlLoop,[1,CONTROL_LOOP_NODE],solidControlLoop,err,error,*999)
      CASE(PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
        CALL ControlLoop_Get(rootControlLoop,[1,1,CONTROL_LOOP_NODE],solidControlLoop,err,error,*999)
      END SELECT
      CALL Solvers_SolverGet(solidControlLoop%solvers,1,solverFiniteElasticity,err,error,*999)
      NULLIFY(solverEquationsFiniteElasticity)
      CALL Solver_SolverEquationsGet(solverFiniteElasticity,solverEquationsFiniteElasticity,err,error,*999)
      NULLIFY(solverMappingFiniteElasticity)
      CALL SolverEquations_SolverMappingGet(solverEquationsFiniteElasticity,solverMappingFiniteElasticity,err,error,*999)
      NULLIFY(equationsSetFiniteElasticity)
      CALL SolverMapping_EquationsSetGet(solverMappingFiniteElasticity,equationsSetFiniteElasticity,err,error,*999)
      NULLIFY(dependentFieldFiniteElasticity)
      CALL EquationsSet_DependentFieldGet(equationsSetFiniteElasticity,dependentFieldFiniteElasticity,err,error,*999)
      !No longer needed, since no more 'Field_ParametersToFieldParametersCopy'
      !                         CALL Field_NumberOfComponentsGet(dependentFieldFiniteElasticity, &
      !                           & FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS_dependentFieldFiniteElasticity,err,error,*999)
      
      !--- Get the geometric field for the ALE Darcy equations
      IF(pSpecification(3)==PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE) THEN
        CALL Solvers_SolverGet(solvers,1,solverDarcy,err,error,*999)
      ELSE
        CALL Solvers_SolverGet(solvers,2,solverDarcy,err,error,*999)
      ENDIF
      NULLIFY(solverEquationsDarcy)
      CALL Solver_SolverEquations(solverDarcy,solverEquationsDarcy,err,error,*999)
      NULLIFY(solverMappingDarcy)
      CALL SolverEquations_SolverMappingGet(solverEquationsDarcy,solverMappingDarcy,err,error,*999)
      NULLIFY(equationsSetDarcy)
      CALL SolverMapping_EquationsSetGet(solverMappingDarcy,1,equationsSetDarcy,err,error,*999)
      NULLIFY(geometricFieldDarcy)
      CALL EquationsSet_GeometricFieldGet(equationsSetDarcy,geometricFieldDarcy,err,error,*999)
      !No longer needed, since no more 'Field_ParametersToFieldParametersCopy'
      !                         CALL Field_NumberOfComponentsGet(geometricFieldDarcy, &
      !                           & FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS_geometricFieldDarcy,err,error,*999)
      
      !--- Copy the result from Finite-elasticity's dependent field to ALE Darcy's geometric field
      !--- First: FIELD_MESH_DISPLACEMENT_SET_TYPE = - FIELD_PREVIOUS_VALUES_SET_TYPE
      alpha=-1.0_DP
      CALL Field_ParameterSetsCopy(geometricFieldDarcy,FIELD_U_VARIABLE_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE, &
        & FIELD_MESH_DISPLACEMENT_SET_TYPE,alpha,err,error,*999)

      ! Write 'FIELD_PREVIOUS_VALUES_SET_TYPE'
      IF(diagnostics3) THEN
        NULLIFY( dummyValues2 )
        CALL Field_ParameterSetDataGet(geometricFieldDarcy,FIELD_U_VARIABLE_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE, &
          & dummyValues2,err,error,*999)
        numberDOFsToPrint = SIZE(dummyValues2,1)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberDOFsToPrint,numberDOFsToPrint,numberDOFsToPrint,dummyValues2, &
          & '(" geometricFieldDarcy,FIELD_U_VARIABLE_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE = ",4(X,E13.6))',&
          & '4(4(X,E13.6))',err,error,*999)
        CALL Field_ParameterSetDataRestore(geometricFieldDarcy,FIELD_U_VARIABLE_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE, &
          & dummyValues2,err,error,*999)
      ENDIF
      
      !--- Second: Get a pointer to the solution values of the solid
      !    (deformed absolute positions in x, y, z; possibly solid pressure)
      CALL Field_ParameterSetDataGet(dependentFieldFiniteElasticity,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
        & solutionValuesSolid,err,error,*999)
      !                 CALL Field_ParameterSetDataRestore(dependentFieldFiniteElasticity,FIELD_U_VARIABLE_TYPE, &
      !                   & FIELD_VALUES_SET_TYPE,solutionValuesSolid,err,error,*999) ! necessary ???
      
      ! Write 'dependentFieldFiniteElasticity'
      IF(diagnostics3) THEN
        NULLIFY( dummyValues2 )
        CALL Field_ParameterSetDataGet(dependentFieldFiniteElasticity,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & dummyValues2,err,error,*999)
        numberDOFsToPrint = SIZE(dummyValues2,1)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberDOFsToPrint,numberDOFsToPrint,numberDOFsToPrint,dummyValues2, &
          & '(" dependentFieldFiniteElasticity,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE = ",4(X,E13.6))',&
          & '4(4(X,E13.6))',err,error,*999)
        CALL Field_ParameterSetDataRestore(dependentFieldFiniteElasticity,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & dummyValues2,err,error,*999)
      ENDIF
      
      !--- Third: FIELD_MESH_DISPLACEMENT_SET_TYPE += Deformed absolute position of solid
      CALL Field_NumberOfDOFsGet(geometricFieldDarcy,FIELD_U_VARIABLE_TYPE,numberOfDOFs,err,error,*999)
      DO dofNumber=1,numberOfDofs
        ! assumes fluid-geometry and solid-dependent mesh are identical \todo: introduce check
        CALL Field_ParameterSetAddLocalDOF(geometricFieldDarcy,FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE, &
          & dofNumber,solutionValuesSolid(dofNumber),err,error,*999)
        
        !---              !!! Why not directly do the mesh update here ??? !!!
        CALL Field_ParameterSetUpdateLocalDOF(geometricFieldDarcy,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,dofNumber, &
          & solutionValuesSolid(dofNumber),err,error,*999)
        !---
        
      ENDDO !dofNumber
      CALL Field_ParameterSetUpdateStart(geometricFieldDarcy,FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE,err,error,*999)
      CALL Field_ParameterSetUpdateStart(geometricFieldDarcy,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
      
      CALL Field_ParameterSetUpdateFinish(geometricFieldDarcy,FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE,err,error,*999)
      CALL Field_ParameterSetUpdateFinish(geometricFieldDarcy,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
      
      ! Write 'FIELD_MESH_DISPLACEMENT_SET_TYPE'
      IF(diagnostics3) THEN
        NULLIFY( dummyValues2 )
        CALL Field_ParameterSetDataGet(geometricFieldDarcy,FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE, &
          & dummyValues2,err,error,*999)
        numberDOFsToPrint = SIZE(dummyValues2,1)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberDOFsToPrint,numberDOFsToPrint,numberDOFsToPrint,dummyValues2, &
          & '(" geometricFieldDarcy,FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE = ",4(X,E13.6))',&
          & '4(4(X,E13.6))',err,error,*999)
        CALL Field_ParameterSetDataRestore(geometricFieldDarcy,FIELD_U_VARIABLE_TYPE, &
          & FIELD_MESH_DISPLACEMENT_SET_TYPE,dummyValues2,err,error,*999)
      ENDIF
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a Darcy equation fluid type of a fluid mechanics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Darcy_PreSolveGetSolidDisplacement")
    RETURN
999 ERRORSEXITS("Darcy_PreSolveGetSolidDisplacement",err,error)
    RETURN 1

  END SUBROUTINE Darcy_PreSolveGetSolidDisplacement

  !
  !================================================================================================================================

  !>Store solution of previous subiteration iterate
  SUBROUTINE Darcy_PreSolveStorePreviousIterate(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx,esSpecification(3),numberOfEquationsSets,pSpecification(3),solverOutputType,variableType
    REAL(DP) :: alpha
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField
    TYPE(FieldVariableType), POINTER :: fieldVariable
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(VARYING_STRING) :: localError

    ENTERS("Darcy_PreSolveStorePreviousIterate",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    CALL Solver_OutputTypeGet(solver,solverOutputType,err,error,*999)

    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STANDARD_DARCY_SUBTYPE)
      ! do nothing
    CASE(PROBLEM_QUASISTATIC_DARCY_SUBTYPE)
      ! do nothing
    CASE(PROBLEM_ALE_DARCY_SUBTYPE,PROBLEM_PGM_DARCY_SUBTYPE,PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE, &
      & PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE,PROBLEM_PGM_TRANSIENT_DARCY_SUBTYPE, &
      & PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
      NULLIFY(solverEquations)
      CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      !loop over the equations sets
      CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
      DO equationsSetIdx=1,numberOfEquationsSets
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
        CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
        SELECT CASE(esSpecification(3))
        CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE)
          ! do nothing
        CASE(EQUATIONS_SET_QUASISTATIC_DARCY_SUBTYPE)
          ! do nothing
        CASE(EQUATIONS_SET_ALE_DARCY_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
          NULLIFY(dependentField)
          CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
          IF(solverOutputType>=SOLVER_PROGRESS_OUTPUT) THEN
            CALL WriteStirng(GENERAL_OUTPUT_TYPE,'-------------------------------------------------------',err,error,*999)
            CALL WriteStirng(GENERAL_OUTPUT_TYPE,'+++     Storing previous subiteration iterate       +++',err,error,*999)
            CALL WriteStirng(GENERAL_OUTPUT_TYPE,'-------------------------------------------------------',err,error,*999)
          ENDIF
          !--- Store the DEPENDENT field values of the previous subiteration iterate
          NULLIFY(equations)
          CALL EquationSet_EquationsGet(equationsSet,equations,err,error,*999)
          NULLIFY(vectorEquations)
          CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
          NULLIFY(vectorMapping)
          CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
          SELECT CASE(esSpecification(3))
          CASE(EQUATIONS_SET_ALE_DARCY_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE)
            NULLIFY(linearMapping)
            CALL EquationsMapping_LinearMappingGet(vectorMapping,linearMapping,err,error,*999)
            NULLIFY(fieldVariable)
            CALL EquationsMappingLinear_LinearMatrixVariableGet(linearMapping,1,fieldVariable,err,error,*999)
          CASE(EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
            & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
            & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
            NULLIFY(dynamicMapping)
            CALL EquationsMappingVector_DynamicMappingGet(vectorMapping,dynamicMapping,err,error,*999)
            NULLIFY(fieldVariable)
            CALL EquationsMappingDynamic_DynamicVariableGet(dynamicMapping,fieldVariable,err,error,*999)
          END SELECT
          alpha = 1.0_DP
          CALL FieldVariable_ParameterSetsCopy(fieldVariable,FIELD_VALUES_SET_TYPE,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE, &
            & alpha,err,error,*999)
        CASE DEFAULT
          localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
            & " is not valid for a Darcy equation fluid type of a fluid mechanics problem class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDDO !equationsSetIdx
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a Darcy equation fluid type of a fluid mechanics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Darcy_PreSolveStorePreviousIterate")
    RETURN
999 ERRORSEXITS("Darcy_PreSolveStorePreviousIterate",err,error)
    RETURN 1

  END SUBROUTINE Darcy_PreSolveStorePreviousIterate

  !
  !================================================================================================================================
  !
  !updates the boundary conditions etc to the required analytic values
  !for the case EquationsSetIncompElastDarcyAnalyticDarcy the pressure field obtained from the finite elasticity solve is overwritten
  !by the appropriate mass increase for that time step
  SUBROUTINE Darcy_PreSolveUpdateAnalyticValues(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: analyticFunctionType,equationsSetIdx,loopIdx,numberOfDimensions,numberOfEquationsSets,numberOfVariables, &
      & pSpecification(3),variableIdx,variableType
    REAL(DP) :: A1,currentTime,D1,timeIncrement
    REAL(DP), POINTER :: geometricParameters(:)
    TYPE(ControlLoopType), POINTER :: controlLoop,controlTimeLoop
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(EquationsSetAnalyticType), POINTER :: equationsAnalytic
    TYPE(FieldType), POINTER :: dependentField,geometricField
    TYPE(FieldVariableType), POINTER :: fieldVariable,geometricVariable
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(VARYING_STRING) :: localError

    ENTERS("Darcy_PreSolveUpdateAnalyticValues",err,error,*999)

    A1=0.4_DP
    D1=1.0_DP
    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    
    CALL ControlLoop_CurrentTimesGet(controlLoop,currentTime,timeIncrement,err,error,*999)
    
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
      NULLIFY(solverEquations)
      CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      !loop over all the equation sets and set the appropriate field variable type BCs and
      !the source field associated with each equation set
      CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
      DO equationsSetIdx=1,numberOfEquationsSets
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
        NULLIFY(equationsAnalytic)
        CALL EquationsSet_AnalyticExists(equationsSet,equationsAnalytic,err,error,*999)
        IF(ASSOCIATED(equationsAnalytic)) THEN
          CALL EquationsSet_AnalyticFunctionTypeGet(equationsSet,analyticFunctionType,err,error,*999)
          IF(analyticFunctionType==EQUATIONS_SET_INCOMP_ELAST_DARCY_ANALYTIC_DARCY)THEN
            !for this analytic case we copy the mass variable to the pressure variable
            NULLIFY(geometricField)
            CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
            NULLIFY(geometricVariable)
            CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,geometricVariable,err,error,*999)
            NULLIFY(geometricParameters)
            CALL FieldVariable_ParameterSetDataGet(geometricVariable,FIELD_VALUES_SET_TYPE,geometricParameters,err,error,*999)
            CALL FieldVariable_NumberOfComponentsGet(geometricVariable,numberOfDimensions,err,error,*999)
            NULLIFY(dependentField)
            CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
            CALL Field_NumberOfVariablesGet(dependentField,numberOfVariables,err,error,*999)
            !DO variableIdx=1,numberOfVariables
              
            variableType=FIELD_V_VARIABLE_TYPE
            NULLIFY(fieldVariable)
            CALL Field_VariableGet(dependentField,variableType,fieldVariable,err,error,*999)
            !DO componentIdx=4,fieldVariable%numberOfComponents

            
            CALL Field_ParametersToFieldParametersCopy(dependentField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,4, &
              & dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,4,err,error,*999)

            !CALL FieldVariable_ComponentInterpolationCheck(fieldVariable,componentIdx,FIELD_NODE_BASED_INTERPOLATION, &
            !  & err,error,*999)
            !NULLIFY(domain)
            !CALL FieldVariable_ComponentDomainGet(fieldVariable,componentIdx,domain,err,error,*999)
            !NULLIFY(domainTopology)
            !CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
            !NULLIFY(domainNodes)
            !CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
            !!Loop over the local nodes excluding the ghosts.
            !CALL DomainNodes_NumberOfNodesGet(domainNodes,numberOfNodes,err,error,*999)
            !DO nodeIdx=1,numberOfNodes
            !  CALL Field_ParameterSetGetNode(dependentField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,nodeIdx,4, &
            !    & massIncrease,err,error,*999)
            !  CALL Field_ParameterSetUpdateNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,nodeIdx,4, &
            !    & 0.1*massIncrease,err,error,*999)                      
            !  !TODO \todo We should interpolate the geometric field here and the node position.
            !  DO dimensionIdx=1,numberOfDimensions
            !    CALL FieldVariable_LocalDOFIdxGet(geometricVariable,1,1,nodeIdx,componentIdx,localDOFIdx,err,error,*999)
            !    x(dimensionIdx)=geometricParameters(localDOFIdx)
            !  ENDDO !dimensionIdx
            !  !Loop over the derivatives
            !  CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,nodeIdx,numberOfNodeDerivatives,err,error,*999)
            !  DO derivativeIdx=1,numberOfNodeDerivatives
            !    CALL DomainNodes_DerivativeGlobalIndexGet(domainNodes,derivativeIdx,nodeIdx,globalDerivativeIndex,err,error,*999)
            !    CALL Diffusion_AnalyticFunctions(VALUE,X,currentTime,variableType,globalDerivativeIndex,analyticFunctionType, &
            !      & err,error,*999)
            !    CALL FieldVariable_LocalDOFIdxGet(fieldVariable,1,derivativeIdx,nodeIdx,componentIdx,localDOFIdx,err,error,*999)
            !    CALL Field_ParameterSetUpdateLocalDOF(dependentField,variableType,FIELD_ANALYTIC_VALUES_SET_TYPE,localDOFIdx, &
            !      & VALUE,err,error,*999)
            !    boundaryConditionCheckVariable=solverEquations%%boundaryConditions% &
            !      & boundaryConditionsVariableTypeMap(variableType)%ptr%conditionTypes(localDOFIdx)
            !    IF(boundaryConditionCheckVariable==BOUNDARY_CONDITION_FIXED) THEN
            !      CALL Field_ParameterSetUpdateLocalDOF(dependentFieldvariableType,FIELD_VALUES_SET_TYPE,localDOFIdx,VALUE, &
            !        & err,error,*999)
            !    ENDIF                        
            !    IF(variableType==FIELD_U_VARIABLE_TYPE) THEN
            !      IF(domainNodes%NODES(nodeIdx)%boundaryNode) THEN
            !        !If we are a boundary node then set the analytic value on the boundary
            !        CALL BoundaryConditions_SetLocalDOF(boundaryConditions,variableType,localDOFIdx, &
            !          & BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
            !      ENDIF
            !    ENDIF
            !  ENDDO !derivativeIdx
            !ENDDO !nodeIdx
          !ENDDO !componentIdx
                      
            CALL Field_ParameterSetUpdateStart(dependentField,variableType,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
            CALL Field_ParameterSetUpdateStart(dependentField,variableType,FIELD_VALUES_SET_TYPE,err,error,*999)
            CALL Field_ParameterSetUpdateFinish(dependentField,variableType,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
            CALL Field_ParameterSetUpdateFinish(dependentField,variableType,FIELD_VALUES_SET_TYPE,err,error,*999)
          ENDIF
        ENDIF
        !ENDDO !variableIdx
        CALL Field_ParameterSetDataRestore(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & geometricParameters,err,error,*999)
        CALL Field_ParameterSetUpdateStart(dependentField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
        CALL Field_ParameterSetUpdateFinish(dependentField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)

      ENDDO !equationsSetIdx
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a Darcy equation type of a fluid mechanics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("Darcy_PreSolveUpdateAnalyticValues")
    RETURN
999 ERRORSEXITS("Darcy_PreSolveUpdateAnalyticValues",err,error)
    RETURN 1

  END SUBROUTINE Darcy_PreSolveUpdateAnalyticValues

  !
  !================================================================================================================================
  !
  
  !> Monitor convergence of the Darcy solution
  SUBROUTINE Darcy_MonitorConvergence(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: computationNodeNumber,dofNumber,equationsSetIdx,esSpecification(3),fileunitN,fileunitN1,iterationNumber, &
      & loopType,numberOfDofs,numberOfEquationsSets,pSpecification(3),solverOutputType,variableType
    REAL(DP) :: residualNorm,residualNorm_0
    REAL(DP), POINTER :: iterationValuesN(:),iterationValuesN1(:)
    REAL(DP), PARAMETER :: RESIDUAL_TOLERANCE_RELATIVE=1.0E-05_DP
    REAL(DP), PARAMETER :: RESIDUAL_TOLERANCE_ABSOLUTE=1.0E-10_DP
    CHARACTER(25) :: filename
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(FieldType), POINTER :: dependentField
    TYPE(FieldVariableType), POINTER :: fieldVariable
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(VARYING_STRING) :: filepath,localError
    TYPE(WorkGroupType), POINTER :: workGroup

    ENTERS("Darcy_MonitorConvergence",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    
    NULLIFY(workGroup)
    CALL Solver_WorkGroupGet(solver,workGroup,err,error,*999)
    CALL WorkGroup_GroupNodeNumberGet(workGroup,computationNodeNumber,err,error,*999)
    WRITE(filename,'("Darcy_",I3.3,".conv")') computationNodeNumber
    filepath = "./output/"//filename
    OPEN(UNIT=23, FILE=CHAR(filepath),STATUS='unknown',ACCESS='append')

    CALL Solver_OutputTypeGet(solver,solverOutputType,err,error,*999)

    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STANDARD_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_DARCY_SUBTYPE,PROBLEM_TRANSIENT_DARCY_SUBTYPE, &
      & PROBLEM_ALE_DARCY_SUBTYPE,PROBLEM_PGM_DARCY_SUBTYPE,PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE, &
      & PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE,PROBLEM_PGM_TRANSIENT_DARCY_SUBTYPE)
      ! do nothing
    CASE(PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
      NULLIFY(solverEquations)
      CALL Solver_SolverEquationGet(solver,solverEquations,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
      DO equationsSetIdx=1,numberOfEquationsSets
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
        CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
        SELECT CASE(esSpecification(3))
        CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE,EQUATIONS_SET_QUASISTATIC_DARCY_SUBTYPE, &
          & EQUATIONS_SET_ALE_DARCY_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE)
          ! do nothing
        CASE(EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
          IF(solverOutputType>=SOLVER_PROGRESS_OUTPUT) THEN
            CALL WriteString(GENERAL_OUTPUT_TYPE,"Darcy monitor convergence ... ",err,error,*999)
          ENDIF
          NULLIFY(dependentField)
          CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
          NULLIFY(equations)
          CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
          NULLIFY(vectorEquations)
          CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
          NULLIFY(vectorMapping)
          CALL EquationsVector_VectorMapping(vectorEquations,vectorMapping,err,error,*999)
          SELECT CASE(esSpecification(3))
          CASE(EQUATIONS_SET_ALE_DARCY_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE)
            NULLIFY(linearMapping)
            CALL EquationsMappingVector_LinearMapping(vectorMapping,linearMapping,err,error,*999)
            NULLIFY(fieldVariable)
            CALL EquationsMappingLinear_LinearMatrixVariableGet(linearMapping,1,fieldVariable,err,error,*999)
          CASE(EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
            & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
            & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
            NULLIFY(dynamicMapping)
            CALL EquationsMappingVector_DynamicMappingGet(vectorMapping,dynamicMapping,err,error,*999)
            NULLIFY(fieldVariable)
            CALL EquationsMappingDynamic_DynamicVariableGet(dynamicMapping,fieldVariable,err,error,*999)
          END SELECT
          CALL FieldVariable_VariableTypeGet(fieldVariable,variableType,err,error,*999)
          !iter 1
          NULLIFY(iterationValuesN)
          CALL FieldVariable_ParameterSetDataGet(fieldVariable,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,iterationValuesN, &
            & err,error,*999)          
          !iter 2
          NULLIFY(iterationValuesN1)
          CALL FieldVariable_ParameterSetDataGet(fieldVariable,FIELD_VALUES_SET_TYPE,iterationValuesN1,err,error,*999)
          CALL FieldVariable_NumberOfDOFsGet(fieldVariable,numberOfDOFs,err,error,*999)
          residualNorm = 0.0_DP
          DO dofNumber=1,numberOfDofs
            residualNorm = residualNorm + ( iterationValuesN1(dofNumber) - iterationValuesN(dofNumber) )**2.0_DP
          ENDDO !dofNumber
          residualNorm = SQRT(residualNorm / numberOfDofs)
          CALL ControlLoop_TypeGet(controlLoop,loopType,err,error,*999)
          IF(loopType==CONTROL_WHILE_LOOP_TYPE) THEN
            CALL ControlLoop_IterationNumberGet(controlLoop,iterationNumber,err,error,*999)
            IF(iterationNumber>=2) THEN !Omit initialised solution
              IF(iterationNumber==2) THEN
                residualNorm_0 = residualNorm
                WRITE(23,*) 'residualNorm_0 = ',residualNorm_0
                WRITE(23,*) 'R / R0 :'
                WRITE(23,*) residualNorm / residualNorm_0
              ENDIF
              IF(solverOutputType>=SOLVER_PROGRESS_OUTPUT) THEN
                CALL WriteString(GENERAL_OUTPUT_TYPE,'-------------------------------------------------------',err,error,*999)
                CALL WriteStringValue(GENERAL_OUTPUT_TYPE,'+++     residualNorm   =        +++',residualNorm,err,error,*999)
                CALL WriteStringValue(GENERAL_OUTPUT_TYPE,'+++     residualNorm_0 =        +++',residualNorm_0,err,error,*999)
                CALL WriteStringValue(GENERAL_OUTPUT_TYPE,'+++     R / R_0         =        +++',residualNorm/residualNorm_0, &
                  & err,error,*999)
                CALL WriteString(GENERAL_OUTPUT_TYPE,'-------------------------------------------------------',err,error,*999)
              ENDIF
              !End subiteration loop if residual is small relative to residual in first step
              IF((residualNorm/residualNorm_0)<=RESIDUAL_TOLERANCE_RELATIVE .OR. &
                & residualNorm<=RESIDUAL_TOLERANCE_ABSOLUTE ) THEN
                IF(solverOutputType>=SOLVER_PROGRESS_OUTPUT) THEN
                  CALL WriteString(GENERAL_OUTPUT_TYPE,'++++++++++++++++++++++++++++++++++++',err,error,*999)
                  CALL WriteString(GENERAL_OUTPUT_TYPE,'+++    SUBITERATION CONVERGED    +++',err,error,*999)
                  CALL WriteString(GENERAL_OUTPUT_TYPE,'++++++++++++++++++++++++++++++++++++',err,error,*999)
                ENDIF
                controlLoop%whileLoop%continueLoop=.FALSE.
              ELSE IF(iterationNumber==controlLoop%whileLoop%maximumNumberOfIterations) THEN
                CALL FlagWarning("Subiterations between solid and fluid equations did not converge.",err,error,*999)
              ENDIF
            ENDIF
          ELSE
            CALL FlagError("Darcy_MonitorConvergence must be called with a while control loop.",err,error,*999)
          ENDIF


          !                                   subIterationNumber = controlLoop%whileLoop%iterationNumber
          !
          !                                   WRITE(filename,'("Darcy_DOFs_N_",I2.2,".dat")') subIterationNumber
          !                                   filepath = "./output/"//filename
          !                                   fileunitN = 7777 + 2*subIterationNumber
          !                                   OPEN(UNIT=fileunitN,FILE=CHAR(filepath),STATUS='unknown',ACCESS='append')
          !                                   DO dofNumber=1,numberOfDofs
          !                                     WRITE(fileunitN,*) iterationValuesN(dofNumber)
          !                                   END DO
          !
          !
          !                                   WRITE(filename,'("Darcy_DOFs_N1_",I2.2,".dat")') subIterationNumber
          !                                   filepath = "./output/"//filename
          !                                   fileunitN1 = 7777 + 2*subIterationNumber+1
          !                                   OPEN(UNIT=fileunitN1,FILE=CHAR(filepath),STATUS='unknown',ACCESS='append')
          !                                   DO dofNumber=1,numberOfDofs
          !                                     WRITE(fileunitN1,*) iterationValuesN1(dofNumber)
          !                                   END DO
          
          
          CALL FieldVariable_ParameterSetDataRestore(fieldVariable,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,iterationValuesN, &
            & err,error,*999)
          CALL FieldVariable_ParameterSetDataRestore(fieldVariable,FIELD_VALUES_SET_TYPE,iterationValuesN1,err,error,*999)
        CASE DEFAULT
          localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
            & " is not valid for a Darcy equation fluid type of a fluid mechanics problem class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDDO !equationsSetIdx
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a Darcy equation fluid type of a fluid mechanics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    CLOSE(23)
    CLOSE(fileunitN)
    CLOSE(fileunitN1)

    EXITS("Darcy_MonitorConvergence")
    RETURN
999 ERRORSEXITS("Darcy_MonitorConvergence",err,error)
    RETURN 1
    
  END SUBROUTINE Darcy_MonitorConvergence

  !
  !================================================================================================================================
  !

  !> Accelerate convergence of the Darcy solution
  SUBROUTINE Darcy_AccelerateConvergence(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dofNumber,equationsSetIdx,esSpecification(3),numberOfDofs,numberOfEquationsSets,pSpecification(3), &
      & solverOutputType,variableType
    REAL(DP), POINTER :: iterationValuesN(:),iterationValuesN1(:)
    REAL(DP) :: relaxationParam,acceleratedValue
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField
    TYPE(FieldVariableType), POINTER :: fieldVariable
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(VARYING_STRING) :: localError

    ENTERS("Darcy_AccelerateConvergence",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)

    CALL Solver_OutputTypeGet(solver,solverOutputType,err,error,*999)
    
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STANDARD_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_DARCY_SUBTYPE,PROBLEM_TRANSIENT_DARCY_SUBTYPE, &
      & PROBLEM_ALE_DARCY_SUBTYPE,PROBLEM_PGM_DARCY_SUBTYPE,PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE, &
      & PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE,PROBLEM_PGM_TRANSIENT_DARCY_SUBTYPE)
      ! do nothing
    CASE(PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
      NULLIFY(solverEquations)
      CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
      DO equationsSetIdx=1,numberOfEquationsSets
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
        CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
        SELECT CASE(esSpecification(3))
        CASE(EQUATIONS_SET_STANDARD_DARCY_SUBTYPE,EQUATIONS_SET_QUASISTATIC_DARCY_SUBTYPE, &
          & EQUATIONS_SET_ALE_DARCY_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE)
          ! do nothing
        CASE(EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
          NULLIFY(dependentField)
          CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
          NULLIFY(equations)
          CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
          NULLIFY(vectorEquations)
          CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
          NULLIFY(vectorMapping)
          CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
          SELECT CASE(esSpecification(3))
          CASE(EQUATIONS_SET_ALE_DARCY_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE)
            NULLIFY(linearMapping)
            CALL EquationsMappingVector_LinearMappingGet(vectorMapping,linearMapping,err,error,*999)
            NULLIFY(fieldVariable)
            CALL EquationsMappingLinear_LinearMatrixVariableGet(linearMapping,1,fieldVariable,err,error,*999)
          CASE(EQUATIONS_SET_TRANSIENT_ALE_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
            & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
            & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
            NULLIFY(dynamicMapping)
            CALL EquationsMappingVector_DynamicMappingGet(vectorMapping,dynamicMapping,err,error,*999)
            NULLIFY(fieldVariable)
            CALL EquationsMappingDynamic_DynamicVariableGet(dynamicMapping,fieldVariable,err,error,*999)
          END SELECT
          CALL FieldVariable_VariableTypeGet(fieldVariable,variableType,err,error,*999)
          !iter 1
          NULLIFY(iterationValuesN)
          CALL FieldVariable_ParameterSetDataGet(fieldVariable,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,iterationValuesN, &
            & err,error,*999)
          !iter 2
          NULLIFY(iterationValuesN1)
          CALL FieldVariable_ParameterSetDataGet(fieldVariable,FIELD_VALUES_SET_TYPE,iterationValuesN1,err,error,*999)

          !residualNorm = 0.0_DP
          CALL FieldVariable_NumberOfDOFsGet(fieldVariable,numberOfDOFs,err,error,*999)
          
          !DO dofNumber=1,numberOfDofs
          !  residualNorm = residualNorm + &
          !    & ( iterationValuesN1(dofNumber) - iterationValuesN(dofNumber) )**2.0_DP
          !END DO
          !residualNorm = SQRT(residualNorm / numberOfDofs)
          
          relaxationParam = 2.0_DP  !\ToDo Devise better way of determining optimal Aitken parameter

          IF( controlLoop%whileLoop%iterationNumber>2 )THEN
            IF(solverOutputType>=SOLVER_PROGRESS_OUTPUT) THEN
              CALL WriteString(GENERAL_OUTPUT_TYPE,"Darcy accelerate convergence ... ",err,error,*999)
            ENDIF
            DO dofNumber=1,numberOfDofs
              acceleratedValue = iterationValuesN(dofNumber) &
                & + relaxationParam * ( iterationValuesN1(dofNumber) - iterationValuesN(dofNumber) )
              CALL FieldVariable_ParameterSetUpdateLocalDOF(fieldVariable,FIELD_VALUES_SET_TYPE,dofNumber,acceleratedValue, &
                & err,error,*999)
            ENDDO !dofNumber
            CALL FieldVariable_ParameterSetUpdateStart(fieldVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
            CALL FieldCariable_ParameterSetUpdateFinish(fieldVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
          END IF
          CALL FieldVariable_ParameterSetDataRestore(fieldVariable,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,iterationValuesN, &
            & err,error,*999)
          CALL FieldVariable_ParameterSetDataRestore(fieldVariable,FIELD_VALUES_SET_TYPE,iterationValuesN1,err,error,*999)
        CASE DEFAULT
          localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
            & " is not valid for a Darcy equation fluid type of a fluid mechanics problem class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDDO !equationsSetIdx
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a Darcy equation fluid type of a fluid mechanics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Darcy_AccelerateConvergence")
    RETURN
999 ERRORSEXITS("Darcy_AccelerateConvergence",err,error)
    RETURN 1
    
  END SUBROUTINE Darcy_AccelerateConvergence

  !
  !================================================================================================================================
  !

  !> Allows to set an explicit Darcy mass increase to test finite elasticity
  !> (and only then this function is called, but not for the coupled problem)
  SUBROUTINE Darcy_PostSolveSetMassIncrease(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dofNumber,loopIdx,numberOfDofs,pSpecification(3)
    REAL(DP) :: currentTime,timeIncrement,alpha
    REAL(DP), POINTER :: meshDisplacementValues(:),solutionValuesSolid(:)
    TYPE(ControlLoopType), POINTER :: controlLoop,controlTimeLoop,rootControlLoop,solidControlLoop
    TYPE(EquationsSetType), POINTER :: equationsSetDarcy
    TYPE(FieldType), POINTER :: dependentFieldDarcy
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverType), POINTER :: solverFiniteElasticity,solverDarcy
    TYPE(SolversType), POINTER :: solvers
    TYPE(SolverEquationsType), POINTER :: solverEquationsDarcy
    TYPE(SolverMappingType), POINTER :: solverMappingDarcy
    TYPE(VARYING_STRING) :: localError

    ENTERS("Darcy_PostSolveSetMassIncrease",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    NULLIFY(rootControlLoop)
    CALL Problem_RootControlLoopGet(problem,rootControlLoop,err,error,*999)
    
    CALL ControlLoop_CurrentTimesGet(controlLoop,currentTime,timeIncrement,err,error,*999)
      
    IF(diagnostics1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE, &
        & "*******************************************************************************************************", &
        & err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Current Time   = ",currentTime,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Time Increment = ",timeIncrement,err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE, &
        & "*******************************************************************************************************", &
        & err,error,*999)
    ENDIF
    
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STANDARD_DARCY_SUBTYPE)
      ! do nothing
    CASE(PROBLEM_QUASISTATIC_DARCY_SUBTYPE)
      ! do nothing
    CASE(PROBLEM_TRANSIENT_DARCY_SUBTYPE)
      ! do nothing
    CASE(PROBLEM_ALE_DARCY_SUBTYPE)
      ! do nothing
    CASE(PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE, &
      & PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
      !--- Mass increase specified
      !--- Get the dependent field of the Darcy equations
      NULLIFY(solvers)
      CALL Solver_SolversGet(solver,solvers,err,error,*999)
      NULLIFY(solverDarcy)
      IF(pSpecification(3)==PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE) THEN
        CALL Solvers_SolverGet(solvers,1,solverDarcy,err,error,*999)
      ELSE
        CALL Solvers_SolverGet(solvers,2,solverDarcy,err,error,*999)
      ENDIF
      NULLIFY(solverEquationsDarcy)
      CALL Solver_SolverEquationsGet(solverDarcy,solverEquationsDarcy,err,error,*999)
      NULLIFY(solverMappingDarcy)
      CALL SolverEquations_SolverMappingGet(solverEquationsDarcy,solverMappingDarcy,err,error,*999)
      NULLIFY(equationsSetDarcy)
      CALL SolverMapping_EquationsSetGet(solverMappingDarcy,1,equationsSetDarcy,err,error,*999)
      NULLIFY(dependentFieldDarcy)
      CALL EquationsSet_DependentFieldGet(equationsSetDarcy,dependentFieldDarcy,err,error,*999)
      ! do nothing
      ! Set the mass increase for Darcy dependent field (u, v, w; m)
      
      !                 alpha = 2.0E-03_DP
      
      !                 alpha = 5.0E-04_DP * currentTime / timeIncrement
      
      alpha = 5.0E-04_DP * SIN(2.0_DP * PI * currentTime / timeIncrement / 20.0_DP)

      CALL Field_NumberOfDOFsGet(dependentFieldDarcy,FIELD_V_VARIABLE_TYPE,numberOfDOFs,err,error,*999)
 
      DO dofNumber = NINT(3.0/4.0*numberOfDofs) + 1, numberOfDofs
        !'3/4' only works for equal order interpolation in (u,v,w) and p
        CALL Field_ParameterSetUpdateLocalDOF(dependentFieldDarcy,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,dofNumber, &
          & alpha,err,error,*999)
      ENDDO !dofNumber
      CALL Field_ParameterSetUpdateStart(dependentFieldDarcy,FIELD_V_VARIABLE_TYPE, FIELD_VALUES_SET_TYPE,err,error,*999)
      CALL Field_ParameterSetUpdateFinish(dependentFieldDarcy,FIELD_V_VARIABLE_TYPE, FIELD_VALUES_SET_TYPE,err,error,*999)
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a Darcy equation fluid type of a fluid mechanics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Darcy_PostSolveSetMassIncrease")
    RETURN
999 ERRORSEXITS("Darcy_PostSolveSetMassIncrease",err,error)
    RETURN 1
    
  END SUBROUTINE Darcy_PostSolveSetMassIncrease

  !
  !================================================================================================================================
  !

  !\ToDo: enable this penalty formulation also for (quasi-)static; as made available in solver_routines

  !Adds a penalty term to the equilibrium equations to enforce impermeability at certain boundaries
  ! derived from: "FINITE_ELASTICITY_SURFACE_PRESSURE_RESIDUAL_EVALUATE"; same restrictions apply
  SUBROUTINE Darcy_ImpermeableBCViaPenalty(equationsSet,elementNumber,err,error,*)
    
    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: elementNumber
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: componentIdx1,componentIdx2,elementBaseDOFIdx1,elementBaseDOFIdx2,elementDOFIdx1,elementDOFIdx2, &
      & elementFaceIdx1,elementFaceIdx2,elementNodeDerivativeIdx1,elementNodeDerivativeIdx2,elementNodeIdx1,elementNodeIdx2, &
      & faceNodeIdx1,faceNodeIdx2,faceNodeDerivativeIdx1,faceNodeDerivativeIdx2,faceNumber,faceNumberOfGauss, &
      & faceParameterIdx1,faceParameterIdx2,gaussPointIdx,meshComponentNumber,normalComponentIdx,numberOfElementParameters, &
      & numberOfLocalFaces,numberOfNodes,numberOfNodeDerivatives1,numberOfNodeDerivatives2,parameterIdx1,parameterIdx2, &
      & xiNormalDirection
    REAL(DP) :: columnBasis,dZdXi(3,3),dZdXiT(3,3),gaussWeight,G,gFlat(3,3),gSharp(3,3),normalProjection1,normalProjection2, &
      & penaltyParameter,rowBasis,sqrtG,sum
    LOGICAL :: boundaryFace,impermeableBC
    TYPE(BasisType), POINTER :: faceBasis,dependentBasis
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionElementsType), POINTER :: decompositionElements
    TYPE(DecompositionFacesType), POINTER :: decompositionFaces
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainElementsType), POINTER :: domainElements
    TYPE(DomainFacesType), POINTER :: domainFaces
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatrixType), POINTER :: stiffnessMatrix
    TYPE(FieldType), POINTER :: dependentField,independentField
    TYPE(FieldVariableType), POINTER :: dynamicVariable
    TYPE(EquationsSetType), POINTER :: equationsSetDarcy
    TYPE(FieldInterpolationParametersType), POINTER :: faceVelocityInterpParameters,geometricInterpParameters, &
      & independentInterpParameters
    TYPE(FieldInterpolatedPointType), POINTER :: faceInterpPoint,geometricInterpPoint,independentInterpPoint
    TYPE(QuadratureSchemeType), POINTER :: faceQuadratureScheme

    ENTERS("Darcy_ImpermeableBCViaPenalty",err,error,*999)

    !Make this routine conditional on (stiffnessMatrix%updateMatrix)

    penaltyParameter = 1.0e04_DP

    !Grab pointers of interest
    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
    NULLIFY(vectorMapping)
    CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
    NULLIFY(dynamicMapping)
    CALL EquationsMappingVector_DynamicMappingGet(vectorMapping,dynamicMapping,err,error,*999)
    NULLIFY(dynamicVariable)
    CALL EquationMappingDynamic_DynamicVariableGet(dynamicMapping,dynamicVariable,err,error,*999)
    NULLIFY(dynamicMatrices)
    CALL EquationsMatricesVector_DynamicMatricesGet(vectorMatrices,dynamicMatrices,err,error,*999)
    NULLIFY(stiffnessMatrix)
    CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,1,stiffnessMatrix,err,error,*999)
    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
    NULLIFY(decomposition)
    CALL Field_DecompositionGet(dependentField,decomposition,err,error,*999)
    NULLIFY(decompositionTopology)
    CALL Decomposition_DecompositionTopologyGet(decomposition,decompositionTopology,err,error,*999)
    NULLIFY(decompositionElements)
    CALL DecompositionTopology_DecompositionElementsGet(decompositionTopology,decompositionElements,err,error,*999)
    NULLIFY(decompositionFaces)
    CALL DecompositionTopology_DecompositionFacesGet(decompositionTopology,decompositionFaces,err,error,*999)
    NULLIFY(domain)
    CALL FieldVariable_ComponentDomainGet(dynamicVariable,1,domain,err,error,*999)
    NULLIFY(domainTopology)
    CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
    NULLIFY(domainElements)
    CALL DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*999)
    NULLIFY(domainFaces)
    CALL DomainTopology_DomainFacesGet(domainTopology,domainFaces,err,error,*999)
    
    NULLIFY(independentField)
    CALL EquationsSet_IndependentFieldGet(equationsSet,independentField,err,error,*999)
    
    NULLIFY(equationsInterpolation)
    CALL Equations_InterpolationGet(equations,equationsInterpolation,err,error,*999)
    NULLIFY(geometricInterpParameters)
    CALL EquationsInterpolation_GeometricInterpParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
      & geometricInterpParameters,err,error,*999)
    NULLIFY(geometricInterpPoint)
    CALL EquationsInterpolation_GeometricInterpPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
      & geometricInterpPoint,err,error,*999)
    NULLIFY(independentInterpParameters)
    CALL EquationsInterpolation_IndependentInterpParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
      & independentInterpParameters,err,error,*999)
    NULLIFY(independentInterpPoint)
    CALL EquationsInterpolation_IndependentInterpPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
      & independentInterpPoint,err,error,*999)
    NULLIFY(dependentBasis)
    CALL DomainElements_ElementBasisGet(domainElements,elementNumber,dependentBasis,err,error,*999)

    CALL Basis_NumberOfElementParametersGet(dependentBasis,numberOfElementParameters,err,error,*999)
    CALL Basis_NumberOfLocalFacesGet(dependentBasis,numberOfLocalFaces,err,error,*999)
 
    CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,geometricInterpParameters,err,error,*999)
    
     !Calculate penalty term to render surfaces impermeable: Loop over all faces
    DO elementFaceIdx1=1,numberOfLocalFaces
      CALL DecompositionElements_ElementFaceNumberGet(decompositionElements,elementFaceIdx1,elementNumber,faceNumber,err,error,*999)
      CALL DecompositionFaces_FaceBoundaryFaceGet(decompositionFaces,faceNumber,boundaryFace,err,error,*999)

      !Check if it's a boundary face
      IF(boundaryFace) THEN !!temporary until MESH_FACE (or equivalent) is available (decomp face includes ghost faces?)

        !Grab normal xi direction of the face and the other two xi directions
        CALL DecompositionFaces_FaceXiNormalDirectionGet(decompositionFaces,faceNumber,xiNormalDirection,err,error,*999)
        normalComponentIdx=ABS(xiNormalDirection)  ! if xi=0, this can be a negative number
        !         FACE_COMPONENTS=OTHER_XI_DIRECTIONS3(normalComponentIdx,2:3,1)  !Two xi directions for the current face
        !\todo: will FACE_COMPONENTS be a problem with sector elements? Check this.

        ! To find out which faces are set impermeable:
        NULLIFY(faceVelocityInterpParameters)
        CALL EquationsInterpolation_IndependentInterpParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,  &
          & faceVelocityInterpParameters,err,error,*999)
        NULLIFY(faceInterpPoint)
        CALL EquationsInterpolation_IndependentInterpPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,faceInterpPoint, &
          & err,error,*999)
        CALL Field_InterpolationParametersFaceGet(FIELD_VALUES_SET_TYPE,faceNumber,independentInterpParameters,err,error,*999)

        !Check if impermeable boundary condition is defined on the face
        impermeableBC=.FALSE.
        IF(ANY(ABS(independentInterpParameters%parameters(:,normalComponentIdx))>ZERO_TOLERANCE)) impermeableBC=.TRUE.

        IF(impermeableBC) THEN

          !Grab some other pointers
          NULLIFY(faceBasis)
          CALL DomainFaces_FaceBasisGet(domainFaces,faceNumber,faceBasis,err,error,*999)
          NULLIFY(faceQuadratureScheme)
          CALL Basis_QuadratureSchemeGet(faceBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,faceQuadratureScheme,err,error,*999)
          CALL BasisQuadratureScheme_NumberOfGaussGet(faceQuadratureScheme,faceNumberOfGauss,err,error,*999)

          !A single faceBasis and dependentBasis should suffice, since we only deal with terms
          !  deriving from velocity test AND trial functions, and moreover use Galerkin,
          !  i.e. same basis functions for test and trial functions

          !Start integrating
!\todo: hopefully all quadrature stuff will always match up between face basis and local face stuff.
! Annoying issue here that p(appl) is interpolated using the face_basis, while dZdXI has to be evaluated
! using the 3D face interpolation... many variables are shared, probably supposed to be the same but I
! can't guarantee it and checking every single thing will be a fair bit of overhead
          DO gaussPointIdx=1,faceNumberOfGauss

            CALL BasisQuadratureScheme_GaussWeightGet(faceQuadratureScheme,gaussPointIdx,gaussWeight,err,error,*999)
            !What happens with surface Jacobian ? sqrtG ? - Apparently contained in normal calculation

            !Use (deformed) Geometric field to obtain delx_j/delxi_M = dZdxi at the face gauss point
            CALL Field_InterpolateLocalFaceGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,elementFaceIdx1,gaussPointIdx, &
              & geometricInterpPoint,err,error,*999)
            
            dZdXi=geometricInterpPoint%values(1:3,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(1:3)) !(component,derivative)

            !Calculate covariant metric tensor
            CALL MatrixTranspose(dZdXi,dZdXiT,err,error,*999)
            CALL MatrixProduct(dZdXiT,dZdXi,gFlat,err,error,*999) !g_ij = dZdXI' * dZdXI
            CALL Invert(gFlat,gSharp,G,err,error,*999) !g^ij = inv(g_ij), G=DET(gFlat)
            sqrtG=SQRT(G)

            !--- L o o p   1 : over element rows (3 velocity components) -----------------------------------
            DO componentIdx1=1,3
              !Calculate g^{normalComponentIdx}M*dZ_j/dxi_M; this apparently includes the face Jacobian
              CALL DotProduct(gSharp(normalComponentIdx,:),dZdXi(componentIdx1,:),normalProjection1,err,error,*999)

              IF(xiNormalDirection<0) normalProjection1=-normalProjection1  !always outward normal

              IF(ABS(normalProjection1)<ZERO_TOLERANCE) CYCLE !Makes it a bit quicker

              elementBaseDOFIdx1 = (componentIdx1-1) * numberOfElementParameters

              CALL Basis_NumberOfNodesGet(faceBasis,numberOfNodes,err,error,*999)
              DO faceNodeIdx1=1,numberOfNodes
                CALL Basis_FaceNodeNumberGet(dependentBasis,faceNodeIdx1,elementFaceIdx1,elementNodeIdx1,err,error,*999)
                CALL Basis_FaceNodeNumberOfDerivativesGet(faceBasis,faceNodeIdx1,elementFaceIdx1,numberOfNodeDerivatives1, &
                  & err,error,*999)
                DO faceNodeDerivativeIdx1=1,numberOfNodeDerivatives1
                  CALL Basis_FaceNodeDerivativeNumberGet(dependentBasis,faceNodeDerivativeIdx1,faceNodeIdx1,elementFaceIdx1, &
                    & elementNodeDerivativeIdx1,err,error,*999)
                  CALL Basis_ElementParameterGet(dependentBasis,elementNodeDerivativeIdx1,elementNodeIdx1,parameterIdx1, &
                    & err,error,*999)
                  CALL Basis_ElementParameterGet(faceBasis,faceNodeDerivativeIdx1,faceNodeIdx1,faceParameterIdx1, &
                    & err,error,*999)
                  elementDOFIdx1=elementBaseDOFIdx1+parameterIdx1

                  CALL BasisQuadratureScheme_GaussBasisFunctionGet(faceQuadratureScheme,faceParameterIdx1,NO_PART_DERIV, &
                    & gaussPointIdx,rowBasis,err,error,*999)

                  !--- L o o p   2 : over element columns (3 velocity components) -----------------------------------
                  DO componentIdx2=1,3
                    !Calculate g^3M*dZ_j/dxi_M
                    CALL DotProduct(gSharp(normalComponentIdx,:),dZdXi(componentIdx2,:),normalProjection2,err,error,*999)

                    IF(xiNormalDirection<0) normalProjection2=-normalProjection2  !always outward normal

                    IF(ABS(normalProjection2)<ZERO_TOLERANCE) CYCLE !Makes it a bit quicker

                    elementBaseDOFIdx2 = (componentIdx2-1) * numberOfElementParameters

                    DO faceNodeIdx2=1,numberOfNodes !nnf
                      CALL Basis_FaceNodeNumberGet(dependentBasis,faceNodeIdx2,elementFaceIdx1,elementNodeIdx2,err,error,*999)
                      CALL Basis_FaceNodeNumberOfDerivativesGet(faceBasis,faceNodeIdx2,elementFaceIdx1,numberOfNodeDerivatives2, &
                        & err,error,*999)
                      DO faceNodeDerivativeIdx2=1,numberOfNodeDerivatives2
                        CALL Basis_FaceNodeDerivativeNumberGet(dependentBasis,faceNodeDerivativeIdx2,faceNodeIdx2, &
                          & elementFaceIdx2,elementNodeDerivativeIdx2,err,error,*999)
                        CALL Basis_ElementParameterGet(dependentBasis,elementNodeDerivativeIdx2,elementNodeIdx2, &
                          & parameterIdx2,err,error,*999)
                        CALL Basis_ElementParameterGet(faceBasis,faceNodeDerivativeIdx2,faceNodeIdx2, &
                          & faceParameterIdx2,err,error,*999)
                        elementDOFIdx2=elementBaseDOFIdx2+parameterIdx2

                        sum = 0.0_DP

                        CALL BasisQuadratureScheme_GaussBasisFunctionGet(faceQuadratureScheme,faceParameterIdx2,NO_PART_DERIV, &
                          & gaussPointIdx,columnBasis,err,error,*999)

                        sum = sum + penaltyParameter * rowBasis * normalProjection1 * sqrtG * &
                                                  & columnBasis * normalProjection2 * sqrtG

                        stiffnessMatrix%elementMatrix%matrix(elementDOFIdx1,elementDOFIdx2) = &
                          & stiffnessMatrix%elementMatrix%matrix(elementDOFIdx1,elementDOFIdx2) + &
                          & sum * gaussWeight

                      ENDDO !elementNodeDerivativeIdx2
                    ENDDO !faceNodeIdx2
                  ENDDO !componentIdx2
                ENDDO !elementNodeDerivativeIdx1
              ENDDO !faceNodeIdx1
            ENDDO !componentIdx1
          ENDDO !gaussPointIdx
        ENDIF !impermeableBC
      ENDIF !boundary face check
    ENDDO !elementFaceIdx1

    EXITS("Darcy_ImpermeableBCViaPenalty")
    RETURN

999 ERRORSEXITS("Darcy_ImpermeableBCViaPenalty",err,error)
    RETURN 1
    
  END SUBROUTINE Darcy_ImpermeableBCViaPenalty

  !
  !================================================================================================================================
  !


END MODULE DarcyEquationsRoutines
