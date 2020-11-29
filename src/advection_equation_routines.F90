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
MODULE AdvectionEquationsRoutines

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
  USE EquationsMatricesAccessRoutines
  USE EquationsSetAccessRoutines
  USE FieldRoutines
  USE FieldAccessRoutines
  USE InputOutput
  USE ISO_VARYING_STRING
  USE Kinds
  USE Maths
  USE MatrixVector
  USE ProblemAccessRoutines
  USE Strings
  USE SolverRoutines
  USE SolverAccessRoutines
  USE Timer
  USE Types

#include "macros.h"

  IMPLICIT NONE

  PUBLIC Advection_EquationsSetSetup
  
  PUBLIC Advection_EquationsSetSolutionMethodSet
  
  PUBLIC Advection_EquationsSetSpecificationSet
  
  PUBLIC Advection_ProblemSpecificationSet
  
  PUBLIC Advection_ProblemSetup
  
  PUBLIC Advection_FiniteElementCalculate
  
  PUBLIC Advection_PreSolve
  
  PUBLIC Advection_PostSolve
  
  PUBLIC Advection_Couple1D0D


CONTAINS

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for an advection equation type of an classical field equations set class.
  SUBROUTINE Advection_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet
    INTEGER(INTG), INTENT(IN) :: solutionMethod
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    INTEGER(INTG) :: esSpecification(3)
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Advection_EquationsSetSolutionMethodSet",err,error,*999)

    CALL EquationSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_ADVECTION_SUBTYPE)
      SELECT CASE(solutionMethod)
      CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
        equationsSet%solutionMethod=EQUATIONS_SET_FEM_SOLUTION_METHOD
      CASE DEFAULT
        localError="The specified solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The third equations set specification of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is not valid for an advection type of an classical field equations set."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
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
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the specification for
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations set specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: subtype
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Advection_EquationsSetSpecificationSet",err,error,*999)
    
    IF(.NOT.ASSOCIATED(equationsSet))  CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(SIZE(specification,1)<3) THEN
      localError="The size of the specified specification array of "// &
        & TRIM(NumberToVString(SIZE(specification,1),"*",err,error))//" is invalid. The size must be >= 3."
      CALL FlagError(localError,err,error,*999)
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
    IF(ALLOCATED(equationsSet%specification)) CALL FlagError("Equations set specification is already allocated.",err,error,*999)
    ALLOCATE(equationsSet%specification(3),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate equations set specification.",err,error,*999)
    equationsSet%specification(1:3)=[EQUATIONS_SET_CLASSICAL_FIELD_CLASS,EQUATIONS_SET_ADVECTION_EQUATION_TYPE,subtype]
     
    EXITS("Advection_EquationsSetSpecificationSet")
    RETURN
999 ERRORSEXITS("Advection_EquationsSetSpecificationSet",err,error)
    RETURN 1
    
  END SUBROUTINE Advection_EquationsSetSpecificationSet
  
  !
  !================================================================================================================================
  !
  
  !>Sets up the advection equation.
  SUBROUTINE Advection_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(EquationsSetSetupType), INTENT(INOUT) :: equationsSetSetup
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    INTEGER(INTG) :: componentIdx,esSpecification(3),geometricMeshComponent,geometricScalingType,geometricComponentNumber, &
      & numberOfDependentComponents,numberOfDependentVariables,numberOfDimensions,numberOfIndependentComponents,solutionMethod
    TYPE(DecompositionType), POINTER :: geometricDecomposition
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsSetIndependentType), POINTER :: equationsIndependent
    TYPE(EquationsSetMaterialsType), POINTER :: equationsMaterials
    TYPE(FieldType), POINTER :: geometricField
    TYPE(RegionType), POINTER :: region
    TYPE(VARYING_STRING) :: localError

    ENTERS("Advection_EquationsSetSetup",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)

    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_ADVECTION_SUBTYPE)
      !OK
    CASE DEFAULT
      localError="The thrid equations set specification of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is invalid for an advection equation set."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    SELECT CASE(equationsSetSetup%setupType)
    CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        CALL Advection_EquationsSetSolutionMethodSet(equationsSet,EQUATIONS_SET_FEM_SOLUTION_METHOD,err,error,*999)
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        !Do nothing
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for an advection equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        !Do nothing
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        !Do nothing
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for an advection equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
      !-----------------------------------------------------------------
      ! D e p e n d e n t   f i e l d
      !-----------------------------------------------------------------
      numberOfDependentVariables=2
      numberOfDependentComponents=1
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        NULLIFY(geometricField)
        CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
        IF(equationsSet%dependent%dependentFieldAutoCreated) THEN
          !Create the auto created dependent field
          !start field creation with name 'Concentration'
          NULLIFY(region)
          CALL EquationsSet_RegionGet(equationsSet,region,err,error,*999)
          CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsSet%dependent%dependentField,err,error,*999)
          !start creation of a new field
          CALL Field_TypeSetAndLock(equationsSet%dependent%dependentField,FIELD_GENERAL_TYPE,err,error,*999)
          !label the field
          CALL Field_LabelSet(equationsSet%dependent%dependentField,"Concentration",err,error,*999)
          !define new created field to be dependent
          CALL Field_DependentTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_DEPENDENT_TYPE,err,error,*999)
          !look for decomposition rule already defined
          CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
          !apply decomposition rule found on new created field
          CALL Field_DecompositionSetAndLock(equationsSet%dependent%dependentField,geometricDecomposition,err,error,*999)
          !point new field to geometric field
          CALL Field_GeometricFieldSetAndLock(equationsSet%dependent%dependentField,geometricField,err,error,*999)
          !set number of variables to 2 (U,delU/delN)
          CALL Field_NumberOfVariablesSetAndLock(equationsSet%dependent%dependentField,numberOfDependentVariables,err,error,*999)
          CALL Field_VariableTypesSetAndLock(equationsSet%dependent%dependentField,[FIELD_U_VARIABLE_TYPE, &
            & FIELD_DELUDELN_VARIABLE_TYPE],err,error,*999)
          !set dimension
          CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
            & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
            & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
          !set data type
          CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
            & FIELD_DP_TYPE,err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
            & FIELD_DP_TYPE,err,error,*999)
          CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
          !number of components for U,delU/delN (C)
          CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
            & numberOfDependentComponents,err,error,*999)
          CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
            & numberOfDependentComponents,err,error,*999)
          CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricMeshComponent,err,error,*999)
          !Default to the geometric interpolation setup
          DO componentIdx=1,numberOfDependentComponents
            CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,componentIdx, &
              & geometricMeshComponent,err,error,*999)
            CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,componentIdx, &
              & geometricMeshComponent,err,error,*999)
          ENDDO !componentIdx
          CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
          SELECT CASE(solutionMethod)
          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
            DO componentIdx=1,numberOfDependentComponents
              CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
                & componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
                & componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
            ENDDO !componentIdx
            CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
            CALL Field_ScalingTypeSet(equationsSet%dependent%dependentField,geometricScalingType,err,error,*999)
          CASE(EQUATIONS_SET_NODAL_SOLUTION_METHOD)
            DO componentIdx=1,numberOfDependentComponents
              CALL Field_ComponentInterpolationSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
                & componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
                & componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
            ENDDO !componentIdx
            CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
            CALL Field_ScalingTypeSet(equationsIndependent%independentField,geometricScalingType,err,error,*999)
          CASE DEFAULT
            localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE 
          !Check the user specified field advection equations
          CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
          CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_DEPENDENT_TYPE,err,error,*999)
          CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,numberOfDependentVariables,err,error,*999)
          CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE],err,error,*999)
          CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
          CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfDependentComponents, &
            & err,error,*999)
          CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,numberOfDependentComponents, &
            & err,error,*999)
          CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
          SELECT CASE(solutionMethod)
          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
            DO componentIdx=1,numberOfDependentComponents
              CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,componentIdx, &
                & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,componentIdx, & 
                & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
            ENDDO !componentIdx
          CASE(EQUATIONS_SET_NODAL_SOLUTION_METHOD)
            DO componentIdx=1,numberOfDependentComponents
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
          CALL Field_ParameterSetCreate(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
            & FIELD_BOUNDARY_CONDITIONS_SET_TYPE,err,error,*999)
        ENDIF
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for an advection equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
      !-----------------------------------------------------------------
      ! M a t e r i a l s   f i e l d 
      !-----------------------------------------------------------------
      NULLIFY(equationsMaterials)
      CALL EquationsSet_MaterialsGet(equationsSet,equationsMaterials,err,error,*999)
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        NULLIFY(geometricField)
        CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
        IF(equationsMaterials%materialsFieldAutoCreated) THEN
          !Create the auto created materials field
          !start field creation with name 'Materials Field'
          NULLIFY(region)
          CALL EquationsSet_RegionGet(equationsSet,region,err,error,*999)
          CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsMaterials%materialsField,err,error,*999)
          CALL Field_TypeSetAndLock(equationsMaterials%materialsField,FIELD_MATERIAL_TYPE,err,error,*999)
          !label the field
          CALL Field_LabelSet(equationsMaterials%materialsField,"Materials Field",err,error,*999)
          CALL Field_DependentTypeSetAndLock(equationsMaterials%materialsField,FIELD_INDEPENDENT_TYPE,err,error,*999)
          CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
          !apply decomposition rule found on new created field
          CALL Field_DecompositionSetAndLock(equationsMaterials%materialsField,geometricDecomposition,err,error,*999)
          !point new field to geometric field
          CALL Field_GeometricFieldSetAndLock(equationsMaterials%materialsField,geometricField,err,error,*999)
          CALL Field_NumberOfVariablesSet(equationsMaterials%materialsField,1,err,error,*999)
          CALL Field_VariableTypesSetAndLock(equationsMaterials%materialsField,[FIELD_U_VARIABLE_TYPE],err,error,*999)
          CALL Field_VariableLabelSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,"Materials",err,error,*999)
          CALL Field_DimensionSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
            & err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_NumberOfComponentsSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,1,err,error,*999)
          CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricComponentNumber,err,error,*999)
          CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,1, &
            & geometricComponentNumber,err,error,*999)
          CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,1, &
            & FIELD_CONSTANT_INTERPOLATION,err,error,*999)
          !Default the field scaling to that of the geometric field
          CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
          CALL Field_ScalingTypeSet(equationsMaterials%materialsField,geometricScalingType,err,error,*999)
        ELSE
          !Check the user specified field
          CALL Field_TypeCheck(equationsSetSetup%field,FIELD_MATERIAL_TYPE,err,error,*999)
          CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
          CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,1,err,error,*999)
          CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE],err,error,*999)
          CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
          CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,1,err,error,*999)
        ENDIF
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        IF(equationsMaterials%materialsFieldAutoCreated) THEN
          !Finish creating the materials field
          CALL Field_CreateFinish(equationsMaterials%materialsField,err,error,*999)
          !Set the default values for the materials field (Diffusivity)
          CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & 1,1.0_DP,err,error,*999)
        ENDIF
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*", & 
          & err,error))//" for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*", & 
          & err,error))//" is invalid for Navier-Stokes equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_INDEPENDENT_TYPE)
      !-----------------------------------------------------------------
      ! I n d e p e n d e n t   f i e l d
      !-----------------------------------------------------------------
      NULLIFY(equationsIndependent)
      CALL EquationsSet_IndependentGet(equationsSet,equationsIndependent,err,error,*999)
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        NULLIFY(geometricField)
        CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
        !Create the auto created independent field
        IF(equationsIndependent%independentFieldAutoCreated) THEN
          NULLIFY(region)
          CALL EquationsSet_RegionGet(equationsSet,region,err,error,*999)
          !start field creation with name 'Independent Field'
          CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsIndependent%independentField,err,error,*999)
          !start creation of a new field
          CALL Field_TypeSetAndLock(equationsIndependent%independentField,FIELD_GENERAL_TYPE,err,error,*999)
          !label the field
          CALL Field_LabelSet(equationsIndependent%independentField,"Independent Field",err,error,*999)
          !define new created field to be independent
          CALL Field_DependentTypeSetAndLock(equationsIndependent%independentField,FIELD_INDEPENDENT_TYPE,err,error,*999)
          !look for decomposition rule already defined
          CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
          !apply decomposition rule found on new created field
          CALL Field_DecompositionSetAndLock(equationsIndependent%independentField,geometricDecomposition,err,error,*999)
          !point new field to geometric field
          CALL Field_GeometricFieldSetAndLock(equationsIndependent%independentField,geometricField,err,error,*999)
          !set number of variables to 1 (1 for U)
          CALL Field_NumberOfVariablesSetAndLock(equationsIndependent%independentField,1,err,error,*999)
          CALL Field_VariableTypesSetAndLock(equationsIndependent%independentField,[FIELD_U_VARIABLE_TYPE],err,error,*999)
          CALL Field_DimensionSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
            & err,error,*999)
          !Set the number of components
          numberOfIndependentComponents=2
          CALL Field_NumberOfComponentsSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
            & numberOfIndependentComponents,err,error,*999)
          CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricMeshComponent,err,error,*999)
          !Default to the geometric interpolation setup
          DO componentIdx=1,numberOfIndependentComponents
            CALL Field_ComponentMeshComponentSet(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE,componentIdx, &
              & geometricMeshComponent,err,error,*999)
          ENDDO !componentIdx
          CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
          SELECT CASE(equationsSet%solutionMethod)
          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
            DO componentIdx=1,numberOfIndependentComponents
              CALL Field_ComponentInterpolationSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
                & componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
            ENDDO !componentIdx
            CALL Field_ScalingTypeGet(equationsSet%geometry%geometricField,geometricScalingType,err,error,*999)
            CALL Field_ScalingTypeSet(equationsIndependent%independentField,geometricScalingType,err,error,*999)
          CASE(EQUATIONS_SET_NODAL_SOLUTION_METHOD)
            DO componentIdx=1,numberOfIndependentComponents
              CALL Field_ComponentInterpolationSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
                & componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
            ENDDO !componentIdx
            CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
            CALL Field_ScalingTypeSet(equationsIndependent%independentField,geometricScalingType,err,error,*999)
          CASE DEFAULT
            localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE
          !Check the user specified field- 
          CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
          CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
          CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,1,err,error,*999)
          CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE],err,error,*999)
          CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
          numberOfIndependentComponents=2
          CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfIndependentComponents, &
            & err,error,*999)
          CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
          SELECT CASE(solutionMethod)
          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
            DO componentIdx=1,numberOfIndependentComponents
              CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,componentIdx, &
                & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
            ENDDO !componentIdx
          CASE(EQUATIONS_SET_NODAL_SOLUTION_METHOD)
            DO componentIdx=1,numberOfIndependentComponents
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
        ENDIF
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a standard Navier-Stokes fluid"
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
      !-----------------------------------------------------------------
      ! E q u a t i o n s    t y p e
      !-----------------------------------------------------------------
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        CALL EquationsSet_AssertDependentIsFinished(equationsSet,err,error,*999)
        NULLIFY(equations)
        CALL Equations_CreateStart(equationsSet,equations,err,error,*999)
        CALL Equations_LinearityTypeSet(equations,EQUATIONS_LINEAR,err,error,*999)
        CALL Equations_TimeDependenceTypeSet(equations,EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
        SELECT CASE(solutionMethod)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          !Finish the equations
          CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
          CALL Equations_CreateFinish(equations,err,error,*999)
          NULLIFY(vectorEquations)
          CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
          !Create the equations mapping.
          CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_U_VARIABLE_TYPE,vectorMapping,err,error,*999)
          CALL EquationsMappingVector_DynamicMatricesSet(vectorMapping,.TRUE.,.TRUE.,err,error,*999)
          CALL EquationsMappingVector_DynamicVariableTypeSet(vectorMapping,FIELD_U_VARIABLE_TYPE,err,error,*999)
          CALL EquationsMappingVector_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE,err,error,*999)
          CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
          !Create the equations matrices
          CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
          !Set up matrix storage and structure
          SELECT CASE(equations%sparsityType)
          CASE(EQUATIONS_MATRICES_FULL_MATRICES)
            CALL EquationsMatricesVector_LinearStorageTypeSet(vectorMatrices, &
              & [DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE,DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE],err,error,*999)
          CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
            CALL EquationsMatricesVector_DynamicStorageTypeSet(vectorMatrices, &
              & [DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE, & 
              & DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE],err,error,*999)
            CALL EquationsMatricesVector_DynamicStructureTypeSet(vectorMatrices, &
              & [EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_FEM_STRUCTURE],err,error,*999)
          CASE DEFAULT
            localError="The equations matrices sparsity type of "//TRIM(NumberToVString(equations%sparsityType,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
          CALL EquationsMatrices_VectorCreateFinish(vectorMatrices,err,error,*999)
        CASE DEFAULT
          localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for an advection equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
        & " is invalid for an advection equation."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("Advection_EquationsSetSetup")
    RETURN
999 ERRORSEXITS("Advection_EquationsSetSetup",err,error)
    RETURN 1
    
  END SUBROUTINE Advection_EquationsSetSetup
  
  !
  !================================================================================================================================
  !

  !>Sets the problem specification for an advection problem.
  SUBROUTINE Advection_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to set the problem specification for
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: problemSubtype
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Advection_ProblemSpecificationSet",err,error,*999)
    
    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Problem is not associated.",err,error,*999)
    IF(SIZE(problemSpecification,1)<3) THEN
      localError="The size of the specified problem specification array of "// &
        & TRIM(NumberToVString(SIZE(problemSpecification,1),"*",err,error))//" is invalid. The size should be >= 3."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    problemSubtype=problemSpecification(3)
    SELECT CASE(problemSubtype)
    CASE(PROBLEM_ADVECTION_SUBTYPE)
      !ok
    CASE DEFAULT
      localError="The third problem specification of "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
        & " is not valid for a advection type of a classical field problem."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    IF(ALLOCATED(problem%specification)) CALL FlagError("Problem specification is already allocated.",err,error,*999)      
    ALLOCATE(problem%specification(3),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
    problem%specification(1:3)=[PROBLEM_CLASSICAL_FIELD_CLASS,PROBLEM_ADVECTION_EQUATION_TYPE,problemSubtype]
        
    EXITS("Advection_ProblemSpecificationSet")
    RETURN
999 ERRORSEXITS("Advection_ProblemSpecificationSet",err,error)
    RETURN 1
    
  END SUBROUTINE Advection_ProblemSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the advection equations.
  SUBROUTINE Advection_ProblemSetup(problem,problemSetup,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem
    TYPE(ProblemSetupType), INTENT(INOUT) :: problemSetup
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    INTEGER(INTG) :: pSpecification(3)
    TYPE(ControlLoopType), POINTER :: controlLoop,controlLoopRoot
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverType), POINTER :: solver
    TYPE(SolversType), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Advection_ProblemSetup",err,error,*999)

    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_ADVECTION_SUBTYPE)
      !ok
    CASE DEFAULT
      localError="The third problem specification of "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    SELECT CASE(problemSetup%setupType)
    CASE(PROBLEM_SETUP_INITIAL_TYPE)
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)
        !Do nothing
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Do nothing
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for an advection equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_SETUP_CONTROL_TYPE)
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)
        !Set up a time control loop
        NULLIFY(controlLoop)
        CALL ControlLoop_CreateStart(problem,controlLoop,err,error,*999)
        CALL ControlLoop_TypeSet(controlLoop,CONTROL_TIME_LOOP_TYPE,err,error,*999)
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Finish the control loops
        NULLIFY(controlLoopRoot)
        CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
        CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
        CALL ControlLoop_CreateFinish(controlLoop,err,error,*999)            
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for an advection equation."
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
        CALL Solvers_NumberOfSolversSet(solvers,1,err,error,*999)
        !Set the solver to be a first order dynamic solver
        NULLIFY(solver)
        CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
        CALL Solver_TypeSet(solver,SOLVER_DYNAMIC_TYPE,err,error,*999)
        CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
        !Set solver defaults
        CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
        CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Get the solvers
        NULLIFY(solvers)
        CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
        !Finish the solvers creation
        CALL Solvers_CreateFinish(solvers,err,error,*999)
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for an advection equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)
        !Get the control loop
        NULLIFY(controlLoopRoot)
        CALL Problem_ControlLoopRootGet(problem,controlLoop,err,error,*999)
        NULLIFY(controlLoop)
        CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
        !Get the solver
        NULLIFY(solvers)
        CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
        NULLIFY(solver)
        CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
        !Create the solver equations
        NULLIFY(solverEquations)
        CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
        CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_LINEAR,err,error,*999)
        CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
        CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Get the control loop
        NULLIFY(controlLoopRoot)
        CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
        NULLIFY(controlLoop)
        CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
        !Get the solver equations
        NULLIFY(solvers)
        CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
        NULLIFY(solver)
        CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
        NULLIFY(solverEquations)
        CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
        !Finish the solver equations creation
        CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)             
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for an advection equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
        & " is invalid for an advection equation."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("Advection_ProblemSetup")
    RETURN
999 ERRORSEXITS("Advection_ProblemSetup",err,error)
    RETURN 1
    
  END SUBROUTINE Advection_ProblemSetup
  
  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices for an advection equation finite element equations set.
  SUBROUTINE Advection_FiniteElementCalculate(equationsSet,elementNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet
    INTEGER(INTG), INTENT(IN) :: elementNumber
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    INTEGER(INTG) :: columnComponentIdx,columnElementDOFIdx,columnElementParameterIdx,colsVariableType,componentIdx, &
      & dynamicVariableType,esSpecification(3),gaussPointIdx,numberOfColumnComponents,numberOfColumnElementParameters, &
      & numberOfDependentXi,numberOfDimensions,numberOfGauss,numberOfGeometricXi,numberOfRowComponents, &
      & numberOfRowElementParameters,rowComponentIdx,rowElementDOFIdx,rowElementParameterIdx,rowsVariableType,scalingType,xiIdx     
    REAL(DP) :: area,columndPhidXi(3),columnPhi,conc,D,dConc,dXidX,flow,gaussWeight,jacobian,jacobianGaussWeight,rowdPhidXi(3), &
      & rowPhi,sum
    LOGICAL :: update,updateDamping,updateMatrices,updateStiffness,updateRHS
    TYPE(BasisType), POINTER :: columnBasis,dependentBasis,geometricBasis,rowBasis
    TYPE(DecompositionType), POINTER :: dependentDecomposition,geometricDecomposition
    TYPE(DomainType), POINTER :: columnDomain,dependentDomain,geometricDomain,rowDomain
    TYPE(DomainElementsType), POINTER :: columnDomainElements,dependentDomainElements,geometricDomainElements,rowDomainElements
    TYPE(DomainTopologyType), POINTER :: columnDomainTopology,dependentDomainTopology,geometricDomainTopology,rowDomainTopology
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatrixType), POINTER :: dampingMatrix,stiffnessMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField,independentField,geometricField,materialsField
    TYPE(FieldInterpolationParametersType), POINTER :: colsInterpParameters,dependentInterpParameters,geometricInterpParameters, &
      & independentInterpParameters,materialsInterpParameters,rowsInterpParameters
    TYPE(FieldInterpolatedPointType), POINTER :: dependentInterpPoint,geometricInterpPoint,independentInterpPoint, &
      & materialsInterpPoint
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: geometricInterpPointMetrics
    TYPE(FieldVariableType), POINTER :: colsVariable,geometricVariable,rowsVariable
    TYPE(QuadratureSchemeType), POINTER :: columnQuadratureScheme,dependentQuadratureScheme,geometricQuadratureScheme, &
      & rowQuadratureScheme
    TYPE(RegionType), POINTER :: region
    TYPE(VARYING_STRING) :: localError

    ENTERS("Advection_FiniteElementCalculate",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)

    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_ADVECTION_SUBTYPE)
      !ok
    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is not valid for and advection equation type of a classical field equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(equationsInterpolation)
    CALL Equations_InterpolationGet(equations,equationsInterpolation,err,error,*999)
    NULLIFY(vectorMapping)
    CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
    NULLIFY(lhsMapping)
    CALL EquationsMappingVector_LHSMappingGet(vectorMapping,lhsMapping,err,error,*999)
    NULLIFY(rowsVariable)
    CALL EquationsMappingLHS_LHSVariableGet(lhsMapping,rowsVariable,err,error,*999)
    NULLIFY(dynamicMapping)
    CALL EquationsMappingVector_DynamicMappingGet(vectorMapping,dynamicMapping,err,error,*999)
    NULLIFY(rhsMapping)
    CALL EquationsMappingVector_RHSMappingExists(vectorMapping,rhsMapping,err,error,*999)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
    NULLIFY(dynamicMatrices)
    CALL EquationsMatricesVector_DynamicMatricesGet(vectorMatrices,dynamicMatrices,err,error,*999)
    NULLIFY(stiffnessMatrix)
    CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,1,stiffnessMatrix,err,error,*999)
    CALL EquationsMatrix_UpdateMatrixGet(stiffnessMatrix,updateStiffness,err,error,*999)
    CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,2,dampingMatrix,err,error,*999)
    CALL EquationsMatrix_UpdateMatrixGet(dampingMatrix,updateDamping,err,error,*999)
    NULLIFY(rhsVector)
    updateRHS=.FALSE.
    IF(ASSOCIATED(rhsMapping)) THEN
      CALL EquationsMatricesVector_RHSVectorGet(vectorMatrices,rhsVector,err,error,*999)
      CALL EquationsMatricesRHS_UpdateVectorGet(rhsVector,updateRHS,err,error,*999)
    ENDIF

    updateMatrices=(updateStiffness.OR.updateDamping)
    update=(updateMatrices.OR.updateRHS)

    IF(update) THEN
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
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
      CALL Basis_NumberOfXiGet(geometricBasis,numberOfGeometricXi,err,error,*999)
      
      NULLIFY(dependentField)
      CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
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
      CALL Basis_NumberOfXiGet(geometricBasis,numberOfDependentXi,err,error,*999)
      NULLIFY(dependentQuadratureScheme)
      CALL Basis_QuadratureSchemeGet(dependentBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,dependentQuadratureScheme,err,error,*999)
      CALL BasisQuadrature_NumberOfGaussGet(dependentQuadratureScheme,numberOfGauss,err,error,*999)
       
      NULLIFY(materialsField)
      CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)
      
      NULLIFY(independentField)
      CALL EquationsSet_IndependentFieldGet(equationsSet,independentField,err,error,*999)
    
      NULLIFY(rowsVariable)
      CALL EquationsMappingLHS_LHSVariableGet(lhsMapping,rowsVariable,err,error,*999)
      CALL FieldVariable_VariableTypeGet(rowsVariable,rowsVariableType,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(rowsVariable,numberOfRowComponents,err,error,*999)
      
      NULLIFY(colsVariable)
      CALL EquationsMappingDynamic_DynamicVariableGet(dynamicMapping,colsVariable,err,error,*999)
      CALL FieldVariable_VariableTypeGet(colsVariable,dynamicVariableType,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(colsVariable,numberOfColumnComponents,err,error,*999)
            
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
      
      NULLIFY(dependentInterpParameters)
      CALL EquationsInterpolation_DependentParametersGet(equationsInterpolation,dynamicVariableType,dependentInterpParameters, &
        & err,error,*999)
      NULLIFY(dependentInterpPoint)
      CALL EquationsInterpolation_DependentPointGet(equationsInterpolation,dynamicVariableType,dependentInterpPoint, &
        & err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,dependentInterpParameters,err,error,*999)
      
      NULLIFY(materialsInterpParameters)
      CALL EquationsInterpolation_MaterialsParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & materialsInterpParameters,err,error,*999)
      NULLIFY(materialsInterpPoint)
      CALL EquationsInterpolation_MaterialsPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,materialsInterpPoint, &
        & err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,materialsInterpParameters,err,error,*999)
      
      NULLIFY(independentInterpParameters)
      CALL EquationsInterpolation_IndependentParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & independentInterpParameters,err,error,*999)
      NULLIFY(independentInterpPoint)
      CALL EquationsInterpolation_IndependentPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & independentInterpPoint,err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,independentInterpParameters, &
        & err,error,*999)
      
      !Loop over gauss points
      DO gaussPointIdx=1,numberOfGauss
        
        CALL BasisQuadratureScheme_GaussWeightGet(dependentQuadratureScheme,gaussPointIdx,gaussWeight,err,error,*999)
        
        CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,geometricInterpPoint, &
          & err,error,*999)
        CALL Field_InterpolatedPointMetricsCalculate(numberOfGeometricXi,geometricInterpPointMetrics,err,error,*999)    
        CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,dependentInterpPoint, &
          & err,error,*999)
        CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,independentInterpPoint, &
          & err,error,*999)
        CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,materialsInterpPoint, &
          & err,error,*999)
   
        conc=dependentInterpPoint%values(1,NO_PART_DERIV)
        dConc=dependentInterpPoint%values(1,FIRST_PART_DERIV)
        flow=independentInterpPoint%values(1,NO_PART_DERIV)
        area=independentInterpPoint%values(2,NO_PART_DERIV)
        D=materialsInterpPoint%values(1,NO_PART_DERIV)
        
        !Calculate jacobianGaussWeight.
        CALL FieldInterpolatedPointsMetrics_JacobianGet(geometricInterpPointMetrics,jacobian,err,error,*999)
        jacobianGaussWeight=jacobian*gaussWeight
        
        !Loop over field components
        rowElementDOFIdx=0
        DO rowComponentIdx=1,numberOfRowComponents
          NULLIFY(rowDomain)
          CALL FieldVariable_ComponentDomainGet(rowsVariable,rowComponentIdx,rowDomain,err,error,*999)
          NULLIFY(rowDomainTopology)
          CALL Domain_DomainTopologyGet(rowDomain,rowDomainTopology,err,error,*999)
          NULLIFY(rowDomainElements)
          CALL DomainTopology_DomainElementsGet(rowDomainTopology,rowDomainElements,err,error,*999)
          NULLIFY(rowBasis)
          CALL DomainElements_ElementBasisGet(rowDomainElements,elementNumber,rowBasis,err,error,*999)
          CALL Basis_NumberOfElementParametersGet(rowBasis,numberOfRowElementParameters,err,error,*999)
          NULLIFY(rowQuadratureScheme)
          CALL Basis_QuadratureSchemeGet(rowBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,rowQuadratureScheme,err,error,*999)
          DO rowElementParameterIdx=1,numberOfRowElementParameters
            rowElementDOFIdx=rowElementDOFIdx+1                    
            CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme,rowElementParameterIdx,NO_PART_DERIV, &
              & gaussPointIdx,rowPhi,err,error,*999)
            DO xiIdx=1,numberOfDependentXi
              CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme,rowElementParameterIdx, &
                & PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xiIdx),gaussPointIdx,rowdPhidXi(xiIdx),err,error,*999)
            ENDDO !xiIdx
            IF(updateMatrices) THEN
              !Loop over element columns
              columnElementDOFIdx=0
              DO columnComponentIdx=1,numberOfColumnComponents
                NULLIFY(columnDomain)
                CALL FieldVariable_ComponentDomainGet(colsVariable,columnComponentIdx,columnDomain,err,error,*999)
                NULLIFY(columnDomainTopology)
                CALL Domain_DomainTopologyGet(columnDomain,columnDomainTopology,err,error,*999)
                NULLIFY(columnDomainElements)
                CALL DomainTopology_DomainElementsGet(columnDomainTopology,columnDomainElements,err,error,*999)
                NULLIFY(columnBasis)
                CALL DomainElements_ElementBasisGet(columnDomainElements,elementNumber,columnBasis,err,error,*999)
                CALL Basis_NumberOfElementParametersGet(columnBasis,numberOfColumnElementParameters,err,error,*999)
                NULLIFY(columnQuadratureScheme)
                CALL Basis_QuadratureSchemeGet(columnBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,columnQuadratureScheme,err,error,*999)
                DO columnElementParameterIdx=1,numberOfColumnElementParameters
                  columnElementDOFIdx=columnElementDOFIdx+1
                  CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme,columnElementParameterIdx,NO_PART_DERIV, &
                    & gaussPointIdx,columnPhi,err,error,*999)
                  !Calculate dXidX in 3D
                  dXidX=0.0_DP
                  DO componentIdx=1,numberOfDimensions
                    DO xiIdx=1,numberOfGeometricXi
                      dXidX=dXidX+geometricInterpPointMetrics%dXidX(xiIdx,componentIdx)**2.0_DP
                    ENDDO !xiIdx
                  ENDDO !componentIdx
                  dXidX=SQRT(dXidX)
                  
                  !-- S T I F F N E S S  M A T R I X --!
                  IF(updateStiffness) THEN
                    CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme,columnElementParameterIdx, &
                      & PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(1),gaussPointIdx,columndPhidXi(1),err,error,*999)
                    sum=(D*(rowdPhidXi(1)*dXidX)+(flow/area)*rowPhi)*columndPhidXi(1)*dXidX
                    stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                      & stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)+ &
                      & sum*jacobianGaussWeight
                  ENDIF
                  
                  !-- D A M P I N G  M A T R I X --!
                  IF(updateDamping) THEN
                    sum=rowPhi*columnPhi
                    dampingMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                      & dampingMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)+ &
                      & sum*jacobianGaussWeight
                  ENDIF                  
                  
                ENDDO !columnElementParameterIdx
              ENDDO !columnComponentIdx
            ENDIF !update matrices
            IF(updateRHS) rhsVector%elementVector%vector(rowElementDOFIdx)=0.0_DP
          ENDDO !rowElementParameterIdx
        ENDDO !rowComponentIdx
        
      ENDDO !gaussPointIdx

      !Scale factor adjustment
      CALL Field_ScalingTypeGet(dependentField,scalingType,err,error,*999)
      IF(scalingType/=FIELD_NO_SCALING) THEN
        NULLIFY(rowsInterpParameters)
        CALL EquationsInterpolation_DependentParametersGet(equationsInterpolation,rowsVariableType,rowsInterpParameters, &
          & err,error,*999)
        NULLIFY(colsInterpParameters)
        CALL EquationsInterpolation_DependentParametersGet(equationsInterpolation,dynamicVariableType,colsInterpParameters, &
          & err,error,*999)
        CALL Field_InterpolationParametersScaleFactorsElementGet(elementNumber,rowsInterpParameters,err,error,*999)
        CALL Field_InterpolationParametersScaleFactorsElementGet(elementNumber,colsInterpParameters,err,error,*999)
        !Loop over element rows
        rowElementDOFIdx=0          
        DO rowComponentIdx=1,numberOfRowComponents
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
            columnElementDOFIdx=0
            IF(updateMatrices) THEN
              !Loop over element columns
              DO columnComponentIdx=1,numberOfColumnComponents
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
            IF(updateRHS) THEN
              rhsVector%elementVector%vector(rowElementDOFIdx)= &
                & rhsVector%elementVector%vector(rowElementDOFIdx)* &
                & rowsInterpParameters%scaleFactors(rowElementParameterIdx,rowComponentIdx)
            ENDIF
          ENDDO !rowElementParameterIdx
        ENDDO !rowComponentIdx
      ENDIF !scaling

    ENDIF !update
    
    EXITS("Advection_FiniteElementCalculate")
    RETURN
999 ERRORSEXITS("Advection_FiniteElementCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE Advection_FiniteElementCalculate

  !
  !================================================================================================================================
  !

  !>Sets up the advection problem pre solve.
  SUBROUTINE Advection_PreSolve(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver which has the pre solver operations
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: pSpecification(3)
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemType), POINTER :: problem
    TYPE(VARYING_STRING) :: localError

    ENTERS("Advection_PreSolve",err,error,*999)
 
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
    
    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_ADVECTION_SUBTYPE, &
      & PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
      & PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
      CALL Advection_PreSolveUpdateBC(solver,err,error,*999)
      !CALL Advection_Couple1D0D(SOLVER,err,error,*999)
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
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
  SUBROUTINE Advection_PreSolveUpdateBC(SOLVER,err,error,*)
    !Argument variables
    TYPE(SolverType), POINTER :: SOLVER
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    INTEGER(INTG) :: currentIteration,i,inputIteration,j,n,m,outputIteration,pSpecification(3)
    REAL(DP) :: c(300),conc,currentTime,delta(300),period,s,startTime,stopTime,t(300),timeIncrement
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(FieldType), POINTER :: dependentField
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolversType), POINTER :: solvers
    TYPE(ProblemType), POINTER :: problem
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Advection_PreSolveUpdateBC",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
      & PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
      CALL ControlLoop_CurrentTimeInformationGet(controlLoop,currentTime,timeIncrement,startTime,stopTime,currentIteration, &
        & outputIteration,inputIteration,err,error,*999)
      NULLIFY(solverEquations)
      CALL Solver_SolverEquations(solver,solverEquations,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
      NULLIFY(dependentField)
      CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
      NULLIFY(boundaryConditions)
      CALL SolverEquations_BoundaryConditionsGet(solverEquations,boundaryConditions,err,error,*999)

      IF(currentTime<100.0_DP)THEN

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
        ENDDO !i
        delta(n)=delta(n-1)+(delta(n-1)-delta(n-2))/(t(n-1)-t(n-2))*(t(n)-t(n-1))
        !Find subinterval
        DO j=1,n-1
          IF (t(j) <= (currentTime/period)) THEN
            m=j
          ENDIF
        ENDDO !j
        !Evaluate interpolant
        s=(currentTime/period)-t(m)
        conc=(c(m)+s*delta(m))
      ELSE 
        conc=0.0
      ENDIF
      
      CALL Field_ParameterSetUpdateLocalDOF(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,conc,err,error,*999)
      CALL Field_ParameterSetUpdateStart(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
      CALL Field_ParameterSetUpdateFinish(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
      
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a advection equation type of a classical field problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Advection_PreSolveUpdateBC")
    RETURN
999 ERRORSEXITS("Advection_PreSolveUpdateBC",err,error)
    RETURN 1

  END SUBROUTINE Advection_PreSolveUpdateBC
  
  !
  !================================================================================================================================
  !
  
  !>Perform post-solve operations for an advection problem
  SUBROUTINE Advection_PostSolve(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver for post-solve operations
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    INTEGER(INTG) :: pSpecification(3)
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemType), POINTER :: problem
    TYPE(VARYING_STRING) :: localError

    ENTERS("Advection_PostSolve",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
     
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_ADVECTION_SUBTYPE)
      !Do nothing
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
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
  SUBROUTINE Advection_Couple1D0D(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(FieldType), POINTER :: dependentField
    TYPE(SolverEquationsType), POINTER :: solverEquations  
    TYPE(SolverMappingType), POINTER :: solverMapping 

    ENTERS("Advection_Couple1D0D",err,error,*999)

    NULLIFY(solverEquations)
    CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
    NULLIFY(solverMapping)
    CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
    NULLIFY(equationsSet)
    CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)

    EXITS("Advection_Couple1D0D")
    RETURN
999 ERRORSEXITS("Advection_Couple1D0D",err,error)
    RETURN 1

  END SUBROUTINE

  !
  !================================================================================================================================
  !

END MODULE AdvectionEquationsRoutines

