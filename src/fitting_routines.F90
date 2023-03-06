!> \file
!> \author Sebastian Krittian
!> \brief This module handles all fitting routines.
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
!> Contributor(s): Chris Bradley, Sebastian Krittian, Prasad Babarenda Gamage
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

!>This module handles all fitting routines.
MODULE FittingRoutines

  USE BaseRoutines
  USE BasisRoutines
  USE BasisAccessRoutines
  USE BoundaryConditionsRoutines
  USE Constants
  USE ControlLoopRoutines
  USE ControlLoopAccessRoutines
  USE DataProjectionAccessRoutines
  USE DecompositionAccessRoutines
  USE DistributedMatrixVector
  USE DomainMappings
  USE EquationsRoutines
  USE EquationsAccessRoutines
  USE EquationsMappingRoutines
  USE EquationsMappingAccessRoutines
  USE EquationsMatricesRoutines
  USE EquationsMatricesAccessRoutines
  USE EquationsSetAccessRoutines
  USE FieldRoutines
  USE FieldAccessRoutines
  USE FiniteElasticityUtilityRoutines
  USE FLUID_MECHANICS_IO_ROUTINES
  USE InputOutput
  USE ISO_VARYING_STRING
  USE Kinds
  USE MatrixVector
  USE Maths
  USE ProblemAccessRoutines
  USE Strings
  USE SolverRoutines
  USE SolverAccessRoutines
  USE SolverMappingAccessRoutines
  USE Timer
  USE Types

#include "macros.h"

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC Fitting_EquationsSetSetup
  
  PUBLIC Fitting_EquationsSetSpecificationSet
  
  PUBLIC Fitting_EquationsSetSolutionMethodSet

  PUBLIC Fitting_FiniteElementCalculate

  PUBLIC Fitting_FiniteElementResidualEvaluate

  PUBLIC Fitting_ProblemSetup
  
  PUBLIC Fitting_ProblemSpecificationSet

  PUBLIC Fitting_PreSolve
  
  PUBLIC Fitting_PostSolve

CONTAINS

  !
  !================================================================================================================================
  !

  !>Sets up a fitting equations set.
  SUBROUTINE Fitting_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to setup
    TYPE(EquationsSetSetupType), INTENT(INOUT) :: equationsSetSetup !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx,esSpecification(4),geometricMeshComponent,geometricScalingType,numberOfComponents, &
      & numberOfComponents2,numberOfDecompositionDimensions,numberOfDependentComponents,numberOfDimensions, &
      & numberOfIndependentComponents,numberOfMaterialsComponents,solutionMethod,sparsityType
    TYPE(DecompositionType), POINTER :: geometricDecomposition
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsSetIndependentType), POINTER :: equationsIndependent
    TYPE(EquationsSetMaterialsType), POINTER :: equationsMaterials
    TYPE(EquationsSetSourceType), POINTER :: equationsSource
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField,geometricField,independentField
    TYPE(RegionType), POINTER :: region
    TYPE(VARYING_STRING) :: localError

    ENTERS("Fitting_EquationsSetSetup",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,4,esSpecification,err,error,*999)
    
    SELECT CASE(esSpecification(2))
    CASE(EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE)
      SELECT CASE(esSpecification(3))
      CASE(EQUATIONS_SET_GENERALISED_DATA_FITTING_SUBTYPE, &
        & EQUATIONS_SET_DIFFUSION_TENSOR_FIBRE_DATA_FITTING_SUBTYPE)
        !OK        
      CASE DEFAULT
        localError="The equations set subtype of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
          & " is not valid for a data fitting equations type."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_GAUSS_FITTING_EQUATION_TYPE)
      SELECT CASE(esSpecification(3))
      CASE(EQUATIONS_SET_GENERALISED_GAUSS_FITTING_SUBTYPE, &
        & EQUATIONS_SET_MAT_PROPERTIES_GAUSS_FITTING_SUBTYPE, &
        & EQUATIONS_SET_MAT_PROPERTIES_INRIA_GAUSS_FITTING_SUBTYPE, &
        & EQUATIONS_SET_DIV_FREE_GAUSS_FITTING_SUBTYPE)
        !OK
      CASE DEFAULT
        localError="The equations set subtype of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
          & " is not valid for a Gauss fitting equations type."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The equations set type of "//TRIM(NumberToVString(esSpecification(2),"*",err,error))// &
        & " is not valid for a fitting equations type."
      CALL FlagError(localError,err,error,*999)
    END SELECT
      
    NULLIFY(region)
    CALL EquationsSet_RegionGet(equationsSet,region,err,error,*999)
    
    SELECT CASE(equationsSetSetup%setupType)
    CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
      !-----------------------------------------------------------------
      ! I n i t i a l   S e t u p
      !-----------------------------------------------------------------
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        CALL Fitting_EquationsSetSolutionMethodSet(equationsSet,EQUATIONS_SET_FEM_SOLUTION_METHOD,err,error,*999)
        SELECT CASE(esSpecification(2))
        CASE(EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE)
          SELECT CASE(esSpecification(3))
          CASE(EQUATIONS_SET_GENERALISED_DATA_FITTING_SUBTYPE)
            CALL EquationsSet_LabelSet(equationsSet,"Generalised data fitting equations set",err,error,*999)
          CASE(EQUATIONS_SET_DIFFUSION_TENSOR_FIBRE_DATA_FITTING_SUBTYPE)
            CALL EquationsSet_LabelSet(equationsSet,"Diffusion tensor fibre data fitting equations set",err,error,*999)
          CASE DEFAULT
            localError="The equations set subtype of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
              & " is not valid for a data fitting equations type."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_GAUSS_FITTING_EQUATION_TYPE)
          SELECT CASE(esSpecification(3))
          CASE(EQUATIONS_SET_GENERALISED_GAUSS_FITTING_SUBTYPE)
            CALL EquationsSet_LabelSet(equationsSet,"Generalised Gauss fitting equations set",err,error,*999)            
          CASE(EQUATIONS_SET_MAT_PROPERTIES_GAUSS_FITTING_SUBTYPE)
            CALL EquationsSet_LabelSet(equationsSet,"Material properties Gauss fitting equations set",err,error,*999)            
          CASE(EQUATIONS_SET_MAT_PROPERTIES_INRIA_GAUSS_FITTING_SUBTYPE)
            CALL EquationsSet_LabelSet(equationsSet,"INRIA material properties Gauss fitting equations set",err,error,*999)
          CASE(EQUATIONS_SET_DIV_FREE_GAUSS_FITTING_SUBTYPE)
            CALL EquationsSet_LabelSet(equationsSet,"Divergence free Gauss fitting equations set",err,error,*999)            
          CASE DEFAULT
            localError="The equations set subtype of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
              & " is not valid for a Gauss fitting equations type."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The equations set type of "//TRIM(NumberToVString(esSpecification(2),"*",err,error))// &
            & " is not valid for a fitting equations type."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        !Do nothing
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a Gauss point fitting equations set."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
      !-----------------------------------------------------------------
      ! G e o m e t r i c   S e t u p
      !-----------------------------------------------------------------
      !Do nothing
    CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
      !-----------------------------------------------------------------
      ! D e p e n d e n t   S e t u p
      !-----------------------------------------------------------------
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
      NULLIFY(independentField)
      CALL EquationsSet_IndependentFieldExists(equationsSet,independentField,err,error,*999)
      IF(ASSOCIATED(independentField)) THEN
        CALL Field_NumberOfComponentsGet(independentField,FIELD_U_VARIABLE_TYPE,numberOfIndependentComponents,err,error,*999)
      ELSE
        IF(esSpecification(2)==EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE.AND. &
          & esSpecification(3)==EQUATIONS_SET_DIFFUSION_TENSOR_FIBRE_DATA_FITTING_SUBTYPE) THEN
          !Set the number of components for a diffusion tensor
          numberOfIndependentComponents=NUMBER_OF_VOIGT(numberOfDimensions)
        ELSE
          !Default to 1.
          numberOfIndependentComponents=1
        ENDIF
      ENDIF
      !Set the number of components.
      SELECT CASE(esSpecification(2))
      CASE(EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE)
        SELECT CASE(esSpecification(3))
        CASE(EQUATIONS_SET_GENERALISED_DATA_FITTING_SUBTYPE)
          !If the independent field has been defined use that number of components
          numberOfDependentComponents=numberOfIndependentComponents
        CASE(EQUATIONS_SET_DIFFUSION_TENSOR_FIBRE_DATA_FITTING_SUBTYPE)
          !If the independent field has been defined use that number of components
          numberOfDependentComponents=numberOfIndependentComponents
        CASE DEFAULT
          localError="The equations set subtype of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
            & " is not valid for a data fitting equations type."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_GAUSS_FITTING_EQUATION_TYPE)
        SELECT CASE(esSpecification(3))
        CASE(EQUATIONS_SET_GENERALISED_GAUSS_FITTING_SUBTYPE)
          !If the independent field has been defined use that number of components
          numberOfDependentComponents=numberOfIndependentComponents
        CASE(EQUATIONS_SET_MAT_PROPERTIES_GAUSS_FITTING_SUBTYPE, &
          & EQUATIONS_SET_MAT_PROPERTIES_INRIA_GAUSS_FITTING_SUBTYPE)
          !Component 1: dependent porosity variable, Component 2: dependent permeability variable
          numberOfDependentComponents=2
        CASE(EQUATIONS_SET_DIV_FREE_GAUSS_FITTING_SUBTYPE)
          numberOfDependentComponents=numberOfDimensions+1
        CASE DEFAULT
          localError="The equations set subtype of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
            & " is not valid for a Gauss fitting equations type."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The equations set type of "//TRIM(NumberToVString(esSpecification(2),"*",err,error))// &
          & " is not valid for a fitting equations type."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        !Set start action
        CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
        IF(equationsSet%dependent%dependentFieldAutoCreated) THEN
          !Create the auto created dependent field
          CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsSet%dependent%dependentField,err,error,*999)
          CALL Field_LabelSet(equationsSet%dependent%dependentField,"Dependent Field",err,error,*999)
          !start creation of a new field
          CALL Field_TypeSetAndLock(equationsSet%dependent%dependentField,FIELD_GENERAL_TYPE,err,error,*999)
          !define new created field to be dependent
          CALL Field_DependentTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_DEPENDENT_TYPE,err,error,*999)
          !look for decomposition rule already defined
          NULLIFY(geometricDecomposition)
          CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
          !apply decomposition rule found on new created field
          CALL Field_DecompositionSetAndLock(equationsSet%dependent%dependentField,geometricDecomposition,err,error,*999)
          !point new field to geometric field
          CALL Field_GeometricFieldSetAndLock(equationsSet%dependent%dependentField,geometricField,err,error,*999)
          !set number of variables to 1 (U)
          CALL Field_NumberOfVariablesSetAndLock(equationsSet%dependent%dependentField,1,err,error,*999)
          CALL Field_VariableTypesSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,err,error,*999)
          CALL Field_VariableLabelSet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,"U",err,error,*999)
          CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
            & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
            & FIELD_DP_TYPE,err,error,*999)
          !Set number of components
          CALL Field_NumberOfComponentsSet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
            & numberOfDependentComponents,err,error,*999)
          !Default to the first geometric interpolation setup
          CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricMeshComponent,err,error,*999)
          DO componentIdx=1,numberOfDependentComponents
            CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,componentIdx, &
              & geometricMeshComponent,err,error,*999)
          ENDDO !componentIdx
          SELECT CASE(solutionMethod)
          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
            DO componentIdx=1,numberOfDependentComponents
              CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
                & componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
            ENDDO !componentIdx
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
            localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE
          !Check the user specified field
          CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
          CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_DEPENDENT_TYPE,err,error,*999)
          CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,1,err,error,*999)
          CALL Field_VariableTypesCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,err,error,*999)
          CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_NumberOfComponentsGet(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfComponents,err,error,*999)
          !If the independent field has been defined check the number of components is the same
          IF(ASSOCIATED(independentField)) THEN
            IF((esSpecification(2)==EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE.AND. &
              & esSpecification(3)==EQUATIONS_SET_GENERALISED_DATA_FITTING_SUBTYPE).OR. &
              & (esSpecification(2)==EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE.AND. &
              & esSpecification(3)==EQUATIONS_SET_DIFFUSION_TENSOR_FIBRE_DATA_FITTING_SUBTYPE).OR. &
              & (esSpecification(2)==EQUATIONS_SET_GAUSS_FITTING_EQUATION_TYPE.AND. &
              & esSpecification(3)==EQUATIONS_SET_GENERALISED_GAUSS_FITTING_SUBTYPE)) THEN
              IF(numberOfComponents /= numberOfIndependentComponents) THEN
                localError="The number of components for the specified dependent field of "// &
                  & TRIM(NumberToVString(numberOfComponents,"*",err,error))// &
                  & " does not match the number of components for the independent field of "// &
                  & TRIM(NumberToVString(numberOfIndependentComponents,"*",err,error))//"."
                CALL FlagError(localError,err,error,*999)
              ENDIF
            ELSE
              CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfDependentComponents, &
                & err,error,*999)
            ENDIF
          ELSE
            CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfDependentComponents, &
              & err,error,*999)
          ENDIF
          SELECT CASE(solutionMethod)
          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
            DO componentIdx=1,numberOfComponents
              CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,componentIdx, &
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
            localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ENDIF
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        IF(equationsSet%dependent%dependentFieldAutoCreated) THEN
          !Check that we have the same number of components as the independent field
          CALL Field_NumberOfComponentsGet(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfComponents, &
            & err,error,*999)
          IF(ASSOCIATED(independentField)) THEN
            IF((esSpecification(2)==EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE.AND. &
              & esSpecification(3)==EQUATIONS_SET_GENERALISED_DATA_FITTING_SUBTYPE).OR. &
              & (esSpecification(2)==EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE.AND. &
              & esSpecification(3)==EQUATIONS_SET_DIFFUSION_TENSOR_FIBRE_DATA_FITTING_SUBTYPE).OR. &
              & (esSpecification(2)==EQUATIONS_SET_GAUSS_FITTING_EQUATION_TYPE.AND. &
              & esSpecification(3)==EQUATIONS_SET_GENERALISED_GAUSS_FITTING_SUBTYPE)) THEN
              IF(numberOfComponents /= numberOfIndependentComponents) THEN
                localError="The number of components for the specified dependent field of "// &
                  & TRIM(NumberToVString(numberOfComponents,"*",err,error))// &
                  & " does not match the number of components for the independentt field of "// &
                  & TRIM(NumberToVString(numberOfIndependentComponents,"*",err,error))//"."
                CALL FlagError(localError,err,error,*999)
              ENDIF
            ENDIF
          ENDIF
          !Finish creating the field
          CALL Field_CreateFinish(equationsSet%dependent%dependentField,err,error,*999)
          !Default the values (to zero)
        ENDIF
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for an update-materials Galerkin projection"
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
      IF(((esSpecification(2)==EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE.AND. &
        & esSpecification(3)==EQUATIONS_SET_GENERALISED_DATA_FITTING_SUBTYPE).OR. &
        & (esSpecification(2)==EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE.AND. &
        & esSpecification(3)==EQUATIONS_SET_DIFFUSION_TENSOR_FIBRE_DATA_FITTING_SUBTYPE).OR. &
        & (esSpecification(2)==EQUATIONS_SET_GAUSS_FITTING_EQUATION_TYPE.AND. &
        & esSpecification(3)==EQUATIONS_SET_GENERALISED_GAUSS_FITTING_SUBTYPE).AND. &
        & (esSpecification(4)==EQUATIONS_SET_FITTING_SOBOLEV_VALUE_SMOOTHING.OR. &
        & esSpecification(4)==EQUATIONS_SET_FITTING_SOBOLEV_DIFFERENCE_SMOOTHING)).OR. &
        & (esSpecification(2)==EQUATIONS_SET_GAUSS_FITTING_EQUATION_TYPE.AND. &
        & (esSpecification(3)==EQUATIONS_SET_MAT_PROPERTIES_GAUSS_FITTING_SUBTYPE.OR. &
        & esSpecification(3)==EQUATIONS_SET_MAT_PROPERTIES_INRIA_GAUSS_FITTING_SUBTYPE))) THEN
        !-----------------------------------------------------------------
        ! M a t e r i a l s   S e t u p
        !-----------------------------------------------------------------
        NULLIFY(equationsMaterials)
        CALL EquationsSet_MaterialsGet(equationsSet,equationsMaterials,err,error,*999)
        NULLIFY(geometricField)
        CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
        CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
        NULLIFY(geometricDecomposition)
        CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
        CALL Decomposition_NumberOfDimensionsGet(geometricDecomposition,numberOfDecompositionDimensions,err,error,*999)
        !Sobolev smoothing material parameters
        SELECT CASE(numberOfDecompositionDimensions)
        CASE(1)
          numberOfMaterialsComponents=2
        CASE(2)
          numberOfMaterialsComponents=5
        CASE(3)
          numberOfMaterialsComponents=9
        CASE DEFAULT
          localError="The number of geometric decomposition dimensions of "// &
            & TRIM(NumberToVString(numberOfDecompositionDimensions,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        SELECT CASE(equationsSetSetup%actionType)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
          IF(equationsMaterials%materialsFieldAutoCreated) THEN
            !Create the auto created materials field
            CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsMaterials%materialsField,err,error,*999)
            CALL Field_TypeSetAndLock(equationsMaterials%materialsField,FIELD_MATERIAL_TYPE,err,error,*999)
            CALL Field_DependentTypeSetAndLock(equationsMaterials%materialsField,FIELD_INDEPENDENT_TYPE,err,error,*999)
            !apply decomposition rule found on new created field
            CALL Field_DecompositionSetAndLock(equationsMaterials%materialsField,geometricDecomposition,err,error,*999)
            !point new field to geometric field
            CALL Field_GeometricFieldSetAndLock(equationsMaterials%materialsField,geometricField,err,error,*999)
            CALL Field_NumberOfVariablesSet(equationsMaterials%materialsField,1,err,error,*999)
            CALL Field_VariableTypesSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,err,error,*999)
            CALL Field_DimensionSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
              & err,error,*999)
            CALL Field_DataTypeSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            !Sobolev smoothing material parameters- tau and kappa            
            CALL Field_NumberOfComponentsSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
              & numberOfMaterialsComponents,err,error,*999)
            CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricMeshComponent,err,error,*999)
            DO componentIdx=1,numberOfMaterialsComponents
              CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,componentIdx, &
                & geometricMeshComponent,err,error,*999)
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
            CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,1,err,error,*999)
            CALL Field_VariableTypesCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,err,error,*999)
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
              & err,error,*999)
            CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfMaterialsComponents, &
              & err,error,*999)
          ENDIF
        CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
          IF(equationsMaterials%materialsFieldAutoCreated) THEN
            !Finish creating the materials field
            CALL Field_CreateFinish(equationsMaterials%materialsField,err,error,*999)
            !Set the default values for the materials field
            DO componentIdx=1,numberOfMaterialsComponents
              CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,componentIdx,1.0_DP,err,error,*999)
            ENDDO !componentIdx
          ENDIF
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
            & " is invalid for an update-materials Galerkin projection."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDIF
    CASE(EQUATIONS_SET_SETUP_INDEPENDENT_TYPE)
      IF((esSpecification(2)==EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE.AND. &
        & esSpecification(3)==EQUATIONS_SET_GENERALISED_DATA_FITTING_SUBTYPE).OR. &
        & (esSpecification(2)==EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE.AND. &
        & esSpecification(3)==EQUATIONS_SET_DIFFUSION_TENSOR_FIBRE_DATA_FITTING_SUBTYPE).OR. &
        & (esSpecification(2)==EQUATIONS_SET_GAUSS_FITTING_EQUATION_TYPE.AND. &
        & esSpecification(3)==EQUATIONS_SET_GENERALISED_GAUSS_FITTING_SUBTYPE)) THEN
        !-----------------------------------------------------------------
        ! I n d e p e n d e n t   S e t u p
        ! (this field holds the data point based field of vectors/tensors to map to the dependent field)
        !-----------------------------------------------------------------
        NULLIFY(equationsIndependent)
        CALL EquationsSet_IndependentGet(equationsSet,equationsIndependent,err,error,*999)
        NULLIFY(geometricField)
        CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
        CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
        NULLIFY(dependentField)
        CALL EquationsSet_DependentFieldExists(equationsSet,dependentField,err,error,*999)
        IF(ASSOCIATED(dependentField)) THEN
          CALL Field_NumberOfComponentsGet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
            & numberOfDependentComponents,err,error,*999)
        ELSE
          IF(esSpecification(2)==EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE.AND. &
            & esSpecification(3)==EQUATIONS_SET_DIFFUSION_TENSOR_FIBRE_DATA_FITTING_SUBTYPE) THEN
            numberOfDependentComponents=NUMBER_OF_VOIGT(numberOfDimensions)
          ELSE
            numberOfDependentComponents=1
          ENDIF
        ENDIF
        !Default the number of independent to the number of dependent
        numberOfIndependentComponents=numberOfDependentComponents
        SELECT CASE(equationsSetSetup%actionType)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
          !Set start action
          IF(equationsIndependent%independentFieldAutoCreated) THEN
            !Create the auto created independent field
            !start field creation with name 'INDEPENDENT_FIELD'
            CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsIndependent%independentField,err,error,*999)
            !start creation of a new field
            CALL Field_TypeSetAndLock(equationsIndependent%independentField,FIELD_GENERAL_TYPE,err,error,*999)
            !label the field
            CALL Field_LabelSetAndLock(equationsIndependent%independentField,"Independent Field",err,error,*999)
            !define new created field to be independent
            CALL Field_DependentTypeSetAndLock(equationsIndependent%independentField,FIELD_INDEPENDENT_TYPE,err,error,*999)
            !look for decomposition rule already defined
            NULLIFY(geometricDecomposition)
            CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
            !apply decomposition rule found on new created field
            CALL Field_DecompositionSetAndLock(equationsIndependent%independentField,geometricDecomposition,err,error,*999)
            !point new field to geometric field
            CALL Field_GeometricFieldSetAndLock(equationsIndependent%independentField,geometricField,err,error,*999)
            !Create two variables: U for data and V for weights
            CALL Field_NumberOfVariablesSetAndLock(equationsIndependent%independentField, 2,err,error,*999)
            CALL Field_VariableTypesSetAndLock(equationsIndependent%independentField,[FIELD_U_VARIABLE_TYPE, &
              & FIELD_V_VARIABLE_TYPE],err,error,*999)
            CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricMeshComponent,err,error,*999)
            ! U Variable: data points
            CALL Field_DimensionSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DataTypeSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE, &
              & err,error,*999)
            CALL Field_NumberOfComponentsSet(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
              & numberOfIndependentComponents,err,error,*999)
            !Default to the geometric interpolation setup
            ! V Variable: data point weights
            CALL Field_DimensionSetAndLock(equationsIndependent%independentField,FIELD_V_VARIABLE_TYPE, &
              & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DataTypeSetAndLock(equationsIndependent%independentField,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE, &
              & err,error,*999)
            CALL Field_NumberOfComponentsSet(equationsIndependent%independentField,FIELD_V_VARIABLE_TYPE, &
              & numberOfIndependentComponents,err,error,*999)
            !Default to the geometric interpolation setup
            DO componentIdx=1,numberOfIndependentComponents
              CALL Field_ComponentMeshComponentSet(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE,componentIdx, &
                & geometricMeshComponent,err,error,*999)
              CALL Field_ComponentMeshComponentSet(equationsIndependent%independentField,FIELD_V_VARIABLE_TYPE,componentIdx, &
                & geometricMeshComponent,err,error,*999)
            ENDDO !componentIdx
            CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
            SELECT CASE(solutionMethod)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              !Specify fem solution method
              IF((esSpecification(2)==EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE.AND. &
                & esSpecification(3)==EQUATIONS_SET_GENERALISED_DATA_FITTING_SUBTYPE).OR. &
                & (esSpecification(2)==EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE.AND. &
                & esSpecification(3)==EQUATIONS_SET_DIFFUSION_TENSOR_FIBRE_DATA_FITTING_SUBTYPE)) THEN
                !Data point based
                DO componentIdx = 1,numberOfIndependentComponents
                  CALL Field_ComponentInterpolationSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
                    & componentIdx,FIELD_DATA_POINT_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationSetAndLock(equationsIndependent%independentField,FIELD_V_VARIABLE_TYPE, &
                    & componentIdx,FIELD_DATA_POINT_BASED_INTERPOLATION,err,error,*999)
                ENDDO !componentIdx
              ELSE
                !Gauss point based
                DO componentIdx=1,numberOfIndependentComponents
                  CALL Field_ComponentInterpolationSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
                    & componentIdx,FIELD_GAUSS_POINT_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationSetAndLock(equationsIndependent%independentField,FIELD_V_VARIABLE_TYPE, &
                    & componentIdx,FIELD_GAUSS_POINT_BASED_INTERPOLATION,err,error,*999)
                ENDDO !componentIdx
              ENDIF
              CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
              CALL Field_ScalingTypeSet(equationsIndependent%independentField,geometricScalingType,err,error,*999)
            CASE DEFAULT
              localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ELSE
            !Check the user specified field
            CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
            ! U (vector) variable
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            CALL Field_NumberOfComponentsGet(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfComponents, &
              & err,error,*999)
            !If the dependent field has been defined use that number of components            
            IF(ASSOCIATED(dependentField)) THEN
              IF(numberOfComponents /= numberOfDependentComponents) THEN
                localError="The number of components for the specified independent field of "// &
                  & TRIM(NumberToVString(numberOfComponents,"*",err,error))// &
                  & " does not match the number of components for the dependent field of "// &
                  & TRIM(NumberToVString(numberOfDependentComponents,"*",err,error))//"."
                CALL FlagError(localError,err,error,*999)
              ENDIF
            ENDIF
            ! V (weight) variable
            CALL Field_DimensionCheck(equationsSetSetup%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DataTypeCheck(equationsSetSetup%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,numberOfComponents,err,error,*999)
            CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
            SELECT CASE(solutionMethod)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              IF((esSpecification(2)==EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE.AND. &
                & esSpecification(3)==EQUATIONS_SET_GENERALISED_DATA_FITTING_SUBTYPE).OR. &
                & (esSpecification(2)==EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE.AND. &
                & esSpecification(3)==EQUATIONS_SET_DIFFUSION_TENSOR_FIBRE_DATA_FITTING_SUBTYPE)) THEN
                !Data point based
                DO componentIdx=1,numberOfComponents
                  CALL Field_ComponentInterpolationCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,componentIdx, &
                    & FIELD_DATA_POINT_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationCheck(equationsSetSetup%FIELD,FIELD_V_VARIABLE_TYPE,componentIdx, &
                    & FIELD_DATA_POINT_BASED_INTERPOLATION,err,error,*999)
                ENDDO !componentIdx
              ELSE
                !Gauss point based
                DO componentIdx=1,numberOfComponents
                  CALL Field_ComponentInterpolationCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,componentIdx, &
                    & FIELD_GAUSS_POINT_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationCheck(equationsSetSetup%FIELD,FIELD_V_VARIABLE_TYPE,componentIdx, &
                    & FIELD_GAUSS_POINT_BASED_INTERPOLATION,err,error,*999)
                ENDDO !componentIdx
              ENDIF
            CASE DEFAULT
              localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ENDIF
        CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
          IF(equationsIndependent%independentFieldAutoCreated) THEN
            !Check that we have the same number of components as the dependent field
            IF(ASSOCIATED(dependentField)) THEN
              CALL Field_NumberOfComponentsGet(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfComponents,err,error,*999)
              IF(numberOfComponents /= numberOfDependentComponents) THEN
                localError="The number of components for the specified independent field of "// &
                  & TRIM(NumberToVString(numberOfComponents,"*",err,error))// &
                  & " does not match the number of components for the dependent field of "// &
                  & TRIM(NumberToVString(numberOfDependentComponents,"*",err,error))//"."
                CALL FlagError(localError,err,error,*999)
              ENDIF
            ENDIF
            !Specify finish action
            CALL Field_CreateFinish(equationsIndependent%independentField,err,error,*999)
            !Dafult values
          ENDIF
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
            & " is invalid for a fitting equations set."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDIF
    CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
      IF((esSpecification(2)==EQUATIONS_SET_GAUSS_FITTING_EQUATION_TYPE.AND. &
        & esSpecification(3)==EQUATIONS_SET_DIV_FREE_GAUSS_FITTING_SUBTYPE)) THEN
        !-----------------------------------------------------------------
        ! S o u r c e   S e t u p
        !-----------------------------------------------------------------
        NULLIFY(equationsSource)
        CALL EquationsSet_SourceGet(equationsSet,equationsSource,err,error,*999)
        SELECT CASE(equationsSetSetup%actionType)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
          NULLIFY(geometricField)
          CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
          CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
          IF(equationsSource%sourceFieldAutoCreated) THEN
            !Create the auto created source field
            !start field creation with name 'sourceField'
            CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsSource%sourceField,err,error,*999)
            !start creation of a new field
            CALL Field_TypeSetAndLock(equationsSource%sourceField,FIELD_GENERAL_TYPE,err,error,*999)
            !label the field
            CALL Field_LabelSetAndLock(equationsSource%sourceField,"Source Field",err,error,*999)
            !define new created field to be source
            CALL Field_DependentTypeSetAndLock(equationsSource%sourceField,FIELD_INDEPENDENT_TYPE,err,error,*999)
            !look for decomposition rule already defined
            NULLIFY(geometricDecomposition)
            CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
            !apply decomposition rule found on new created field
            CALL Field_DecompositionSetAndLock(equationsSource%sourceField,geometricDecomposition,err,error,*999)
            !point new field to geometric field
            CALL Field_GeometricFieldSetAndLock(equationsSource%sourceField,geometricField,err,error,*999)
            CALL Field_NumberOfVariablesSetAndLock(equationsSource%sourceField,1,err,error,*999)
            CALL Field_VariableTypesSetAndLock(equationsSource%sourceField,FIELD_U_VARIABLE_TYPE,err,error,*999)
            CALL Field_DimensionSetAndLock(equationsSource%sourceField,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
              & err,error,*999)
            CALL Field_DataTypeSetAndLock(equationsSource%sourceField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            CALL Field_NumberOfComponentsSetAndLock(equationsSource%sourceField,FIELD_U_VARIABLE_TYPE,numberOfDimensions+1, &
              & err,error,*999)
            CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
            SELECT CASE(solutionMethod)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              !Specify fem solution method
              DO componentIdx=1,numberOfDimensions
                CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,componentIdx,geometricMeshComponent, &
                  & err,error,*999)
                !Default to the geometric interpolation setup
                CALL Field_ComponentMeshComponentSet(equationsSource%sourceField,FIELD_U_VARIABLE_TYPE,componentIdx, &
                  & geometricMeshComponent,err,error,*999)
                !Default to node based
                CALL Field_ComponentInterpolationSetAndLock(equationsSource%sourceField,FIELD_U_VARIABLE_TYPE,componentIdx, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              ENDDO !componentIdx
              CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricMeshComponent, &
                & err,error,*999)
              !Default to the geometric interpolation setup
              CALL Field_ComponentMeshComponentSet(equationsSource%sourceField,FIELD_U_VARIABLE_TYPE,numberOfDimensions+1, &
                & geometricMeshComponent,err,error,*999)
              !Default to node based
              CALL Field_ComponentInterpolationSet(equationsSource%sourceField,FIELD_U_VARIABLE_TYPE,numberOfDimensions+1, &
                & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
              CALL Field_ScalingTypeSet(equationsSource%sourceField,geometricScalingType,err,error,*999)
            CASE DEFAULT
              !Other solutions not defined yet
              localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ELSE
            !Check the user specified field
            CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
            CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
            CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,1,err,error,*999)
            CALL Field_VariableTypesCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,err,error,*999)
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            !calculate number of components with one component for each dimension and one for pressure
            CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfDimensions+1,err,error,*999)
            CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
            SELECT CASE(solutionMethod)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              DO componentIdx=1,numberOfDimensions+1
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,componentIdx, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              ENDDO !componentIdx
            CASE DEFAULT
              localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ENDIF
        CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
          IF(equationsSource%sourceFieldAutoCreated) THEN
            CALL Field_CreateFinish(equationsSource%sourceField,err,error,*999)
          ENDIF
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
            & " is invalid for a standard Galerkin projection."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDIF
    CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
      !-----------------------------------------------------------------
      ! E q u a t i o n s   S e t u p
      !-----------------------------------------------------------------
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        CALL EquationsSet_AssertDependentIsFinished(equationsSet,err,error,*999)
        IF(((esSpecification(2)==EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE.AND. &
          & esSpecification(3)==EQUATIONS_SET_GENERALISED_DATA_FITTING_SUBTYPE).OR. &
          & (esSpecification(2)==EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE.AND. &
          & esSpecification(3)==EQUATIONS_SET_DIFFUSION_TENSOR_FIBRE_DATA_FITTING_SUBTYPE).OR. &
          & (esSpecification(2)==EQUATIONS_SET_GAUSS_FITTING_EQUATION_TYPE.AND. &
          & esSpecification(3)==EQUATIONS_SET_GENERALISED_GAUSS_FITTING_SUBTYPE).AND. &
          & (esSpecification(4)==EQUATIONS_SET_FITTING_SOBOLEV_VALUE_SMOOTHING.OR. &
          & esSpecification(4)==EQUATIONS_SET_FITTING_SOBOLEV_DIFFERENCE_SMOOTHING)).OR. &
          & (esSpecification(2)==EQUATIONS_SET_GAUSS_FITTING_EQUATION_TYPE.AND. &
          & (esSpecification(3)==EQUATIONS_SET_MAT_PROPERTIES_GAUSS_FITTING_SUBTYPE.OR. &
          & esSpecification(3)==EQUATIONS_SET_MAT_PROPERTIES_INRIA_GAUSS_FITTING_SUBTYPE))) THEN
          CALL EquationsSet_AssertMaterialsIsFinished(equationsSet,err,error,*999)
        ENDIF
        IF((esSpecification(2)==EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE.AND. &
          & esSpecification(3)==EQUATIONS_SET_GENERALISED_DATA_FITTING_SUBTYPE).OR. &
          & (esSpecification(2)==EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE.AND. &
          & esSpecification(3)==EQUATIONS_SET_DIFFUSION_TENSOR_FIBRE_DATA_FITTING_SUBTYPE).OR. &
          & (esSpecification(2)==EQUATIONS_SET_GAUSS_FITTING_EQUATION_TYPE.AND. &
          & esSpecification(3)==EQUATIONS_SET_GENERALISED_GAUSS_FITTING_SUBTYPE)) THEN
          CALL EquationsSet_AssertIndependentIsFinished(equationsSet,err,error,*999)
        ENDIF
        IF((esSpecification(2)==EQUATIONS_SET_GAUSS_FITTING_EQUATION_TYPE.AND. &
          & esSpecification(3)==EQUATIONS_SET_DIV_FREE_GAUSS_FITTING_SUBTYPE)) THEN
          CALL EquationsSet_AssertSourceIsFinished(equationsSet,err,error,*999)
        ENDIF
        NULLIFY(equations)
        CALL Equations_CreateStart(equationsSet,equations,err,error,*999)
        IF(esSpecification(2)==EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE.AND. &
          & esSpecification(3)==EQUATIONS_SET_DIFFUSION_TENSOR_FIBRE_DATA_FITTING_SUBTYPE) THEN
          CALL Equations_LinearityTypeSet(equations,EQUATIONS_NONLINEAR,err,error,*999)
        ELSE
          CALL Equations_LinearityTypeSet(equations,EQUATIONS_LINEAR,err,error,*999)
        ENDIF
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
          IF(esSpecification(2)==EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE.AND. &
            & esSpecification(3)==EQUATIONS_SET_DIFFUSION_TENSOR_FIBRE_DATA_FITTING_SUBTYPE) THEN
            !Create the equations mapping.
            NULLIFY(vectorMapping)
            CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_U_VARIABLE_TYPE,vectorMapping,err,error,*999)
            CALL EquationsMappingVector_NumberOfLinearMatricesSet(vectorMapping,0,err,error,*999)
            CALL EquationsMappingVector_NumberOfResidualsSet(vectorMapping,1,err,error,*999)
            CALL EquationsMappingVector_ResidualNumberOfVariablesSet(vectorMapping,1,1,err,error,*999)
            CALL EquationsMappingVector_ResidualVariableTypesSet(vectorMapping,1,FIELD_U_VARIABLE_TYPE,err,error,*999)            
            CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
            !Create the equations matrices
            NULLIFY(vectorMatrices)
            CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
            CALL Equations_SparsityTypeGet(equations,sparsityType,err,error,*999)
            SELECT CASE(sparsityType)
            CASE(EQUATIONS_MATRICES_FULL_MATRICES)
              CALL EquationsMatricesVector_NonlinearStorageTypeSet(vectorMatrices,1,MATRIX_BLOCK_STORAGE_TYPE,err,error,*999)
            CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
              CALL EquationsMatricesVector_NonlinearStorageTypeSet(vectorMatrices,1,MATRIX_COMPRESSED_ROW_STORAGE_TYPE, &
                & err,error,*999)
              CALL EquationsMatricesVector_NonlinearStructureTypeSet(vectorMatrices,1,EQUATIONS_MATRIX_FEM_STRUCTURE,err,error,*999)
            CASE DEFAULT
              localError="The equations matrices sparsity type of "//TRIM(NumberToVString(sparsityType,"*",err,error))// &
                & " is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
            CALL EquationsMatrices_VectorCreateFinish(vectorMatrices,err,error,*999)
          ELSE
            !Create the equations mapping.
            NULLIFY(vectorMapping)
            CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_U_VARIABLE_TYPE,vectorMapping,err,error,*999)
            CALL EquationsMappingVector_NumberOfLinearMatricesSet(vectorMapping,1,err,error,*999)
            CALL EquationsMappingVector_LinearMatricesVariableTypesSet(vectorMapping,FIELD_U_VARIABLE_TYPE,err,error,*999)
            CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
            !Create the equations matrices
            NULLIFY(vectorMatrices)
            CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
            CALL Equations_SparsityTypeGet(equations,sparsityType,err,error,*999)
            SELECT CASE(sparsityType)
            CASE(EQUATIONS_MATRICES_FULL_MATRICES)
              CALL EquationsMatricesVector_LinearStorageTypeSet(vectorMatrices,MATRIX_BLOCK_STORAGE_TYPE,err,error,*999)
            CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
              CALL EquationsMatricesVector_LinearStorageTypeSet(vectorMatrices,MATRIX_COMPRESSED_ROW_STORAGE_TYPE, &
                & err,error,*999)
              CALL EquationsMatricesVector_LinearStructureTypeSet(vectorMatrices,EQUATIONS_MATRIX_FEM_STRUCTURE,err,error,*999)
            CASE DEFAULT
              localError="The equations matrices sparsity type of "//TRIM(NumberToVString(sparsityType,"*",err,error))// &
                & " is invalid."
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
          localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a fitting equations set."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
        & " is invalid for a fitting equations set."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Fitting_EquationsSetSetup")
    RETURN
999 ERRORSEXITS("Fitting_EquationsSetSetup",err,error)
    RETURN 1

  END SUBROUTINE Fitting_EquationsSetSetup

  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a data fitting equation set class.
  SUBROUTINE Fitting_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetType,equationsSetSubtype,equationsSetSmoothing
    TYPE(VARYING_STRING) :: localError

    ENTERS("Fitting_EquationsSetSpecificationSet",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(SIZE(specification,1)/=4) THEN
        CALL FlagError("Equations set specification must have four entries for a fitting class equations set.", &
          & err,error,*999)
      ENDIF
      equationsSetType=specification(2)
      SELECT CASE(equationsSetType)
      CASE(EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE)
        equationsSetSubtype=specification(3)
        SELECT CASE(equationsSetSubtype)
        CASE(EQUATIONS_SET_GENERALISED_DATA_FITTING_SUBTYPE, &
          & EQUATIONS_SET_DIFFUSION_TENSOR_FIBRE_DATA_FITTING_SUBTYPE)
          !OK
        CASE DEFAULT
          localError="The third equations set specifiction of "//TRIM(NumberToVstring(equationsSetSubtype,"*",err,error))// &
            & " is not valid for a data fitting equations set."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_GAUSS_FITTING_EQUATION_TYPE)
        equationsSetSubtype=specification(3)
        SELECT CASE(equationsSetSubtype)
        CASE(EQUATIONS_SET_GENERALISED_GAUSS_FITTING_SUBTYPE, &
          & EQUATIONS_SET_MAT_PROPERTIES_GAUSS_FITTING_SUBTYPE, &
          & EQUATIONS_SET_MAT_PROPERTIES_INRIA_GAUSS_FITTING_SUBTYPE, &
          & EQUATIONS_SET_DIV_FREE_GAUSS_FITTING_SUBTYPE)
          !OK
        CASE DEFAULT
          localError="The third equations set specifiction of "//TRIM(NumberToVstring(equationsSetSubtype,"*",err,error))// &
            & " is not valid for a Gauss fitting equations set."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The second equations set specification of "//TRIM(NumberToVstring(equationsSetType,"*",err,error))// &
          & " is not valid for a data fitting equations set."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      equationsSetSmoothing=specification(4)
      SELECT CASE(equationsSetSmoothing)
      CASE(EQUATIONS_SET_FITTING_NO_SMOOTHING, &
        & EQUATIONS_SET_FITTING_SOBOLEV_VALUE_SMOOTHING, &
        & EQUATIONS_SET_FITTING_SOBOLEV_DIFFERENCE_SMOOTHING, &
        & EQUATIONS_SET_FITTING_STRAIN_ENERGY_SMOOTHING)
        !ok
      CASE DEFAULT
        localError="The fourth equations set specification of "//TRIM(NumberToVstring(equationsSetSmoothing,"*",err,error))// &
          & " is not valid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Set full specification
      IF(ALLOCATED(equationsSet%specification)) THEN
        CALL FlagError("Equations set specification is already allocated.",err,error,*999)
      ELSE
        ALLOCATE(equationsSet%specification(4),stat=err)
        IF(err/=0) CALL FlagError("Could not allocate equations set specification.",err,error,*999)
      END IF
      equationsSet%specification(1:4)=[EQUATIONS_SET_FITTING_CLASS,equationsSetType,equationsSetSubtype,equationsSetSmoothing]
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    END IF

    EXITS("Fitting_EquationsSetSpecificationSet")
    RETURN
999 ERRORS("Fitting_EquationsSetSpecificationSet",err,error)
    EXITS("Fitting_EquationsSetSpecificationSet")
    RETURN 1

  END SUBROUTINE Fitting_EquationsSetSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a Galerkin projection type of an data fitting equations set class.
  SUBROUTINE Fitting_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: solutionMethod !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Fitting_EquationsSetSolutionMethodSet",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(.NOT.ALLOCATED(equationsSet%specification)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE
        IF(SIZE(equationsSet%specification,1)/=4) &
          & CALL FlagError("Equations set specification must have four entries for a fitting type equations set.", &
          & err,error,*999)
      ENDIF
      SELECT CASE(equationsSet%specification(2))
      CASE(EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE)
        SELECT CASE(equationsSet%specification(3))
        CASE(EQUATIONS_SET_GENERALISED_DATA_FITTING_SUBTYPE, &
          & EQUATIONS_SET_DIFFUSION_TENSOR_FIBRE_DATA_FITTING_SUBTYPE)
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
          localError="Equations set subtype of "//TRIM(NumberToVString(equationsSet%specification(3),"*",err,error))// &
            & " is not valid for a data fitting equations set class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_GAUSS_FITTING_EQUATION_TYPE)
        SELECT CASE(equationsSet%specification(3))
        CASE(EQUATIONS_SET_GENERALISED_GAUSS_FITTING_SUBTYPE, &
          & EQUATIONS_SET_MAT_PROPERTIES_GAUSS_FITTING_SUBTYPE, &
          & EQUATIONS_SET_MAT_PROPERTIES_INRIA_GAUSS_FITTING_SUBTYPE, &
          & EQUATIONS_SET_DIV_FREE_GAUSS_FITTING_SUBTYPE)
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
          localError="Equations set subtype of "//TRIM(NumberToVString(equationsSet%specification(3),"*",err,error))// &
            & " is not valid for a Gauss fitting equations set class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="Equations set type of "//TRIM(NumberToVString(equationsSet%specification(2),"*",err,error))// &
          & " is not valid for a fitting equations set class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

    EXITS("Fitting_EquationsSetSolutionMethodSet")
    RETURN
999 ERRORSEXITS("Fitting_EquationsSetSolutionMethodSet",err,error)
    RETURN 1

  END SUBROUTINE Fitting_EquationsSetSolutionMethodSet

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices and RHS for fitting with a finite element equations set.
   SUBROUTINE Fitting_FiniteElementCalculate(equationsSet,elementNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG), PARAMETER :: MAXIMUM_DATA_COMPONENTS=99
    INTEGER(INTG) :: colsVariableType,columnComponentIdx,columnElementDOFIdx,columnElementParameterIdx,columnXiIdx,componentIdx, &
      & dataPointGlobalNumber,dataPointIdx,dataPointLocalNumber,dataPointUserNumber,derivativeIdx,esSpecification(4), &
      & gaussPointIdx,geometricDerivative,localDOFIdx,numberOfColsComponents, &
      & numberOfColumnElementParameters(MAXIMUM_DATA_COMPONENTS),numberOfDataComponents,numberOfDimensions,numberOfDOFs, &
      & numberOfElementDataPoints,numberOfElementXi,numberOfGauss,numberOfRowsComponents, &
      & numberOfRowElementParameters(MAXIMUM_DATA_COMPONENTS),numberOfXi,rowComponentIdx,rowElementDOFIdx,rowElementParameterIdx, &
      & rowXiIdx,rowsVariableType,scalingType,smoothingType,xiIdx
    REAL(DP) :: columnPhi,columndPhidXi(3),columndPhidXi1,columndPhidXi2,columndPhidXi3,columnd2PhidXi1dXi1,columnd2PhidXi2dXi2, &
      & columnd2PhidXi3dXi3,columnd2PhidXi1dXi2,columnd2PhidXi1dXi3,columnd2PhidXi2dXi3,curvature, &
      & dataPointVector(MAXIMUM_DATA_COMPONENTS),dataPointWeight(MAXIMUM_DATA_COMPONENTS),dXdY(3,3),dXdXi(3,3),dYdXi(3,3), &
      & dXidY(3,3),dXidX(3,3),elementXi(3),gaussWeight,jacobian,jacobianGaussWeight,Jxy,Jyxi,kappa11,kappa22,kappa33,kappa12, &
      & kappa13,kappa23,kappaParam,materialFact,permOverVisParam0,permOverVisParam,porosity0,porosity,projectionXi(3), &
      & rhsCurvature,rhsTension,rowPhi, &
      & rowdPhidXi(3),rowdPhidXi1,rowdPhidXi2,rowdPhidXi3,rowd2PhidXi1dXi1,rowd2PhidXi2dXi2,rowd2PhidXi3dXi3,rowd2PhidXi1dXi2, &
      & rowd2PhidXi1dXi3,rowd2PhidXi2dXi3,sum,tau1,tau2,tau3,tauParam,tension,uValue(3)
    REAL(DP), POINTER :: independentVectorParameters(:),independentWeightParameters(:)
    LOGICAL :: update,updateMatrix,updateRHS,updateSource
    TYPE(BasisType), POINTER :: dependentBasis,geometricBasis
    TYPE(BasisPtrType) :: columnBasis(MAXIMUM_DATA_COMPONENTS),rowBasis(MAXIMUM_DATA_COMPONENTS)
    TYPE(DataProjectionType), POINTER :: dataProjection
    TYPE(DecompositionType), POINTER :: dependentDecomposition,independentDecomposition,geometricDecomposition
    TYPE(DecompositionDataPointsType), POINTER :: decompositionDataPoints
    TYPE(DecompositionTopologyType), POINTER :: independentDecompositionTopology
    TYPE(DomainType), POINTER :: columnDomain,dependentDomain,geometricDomain,rowDomain
    TYPE(DomainElementsType), POINTER :: columnDomainElements,dependentDomainElements,geometricDomainElements,rowDomainElements
    TYPE(DomainTopologyType), POINTER :: columnDomainTopology,dependentDomainTopology,geometricDomainTopology,rowDomainTopology
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping
    TYPE(EquationsMappingSourceType), POINTER :: sourceMapping
    TYPE(EquationsMappingSourcesType), POINTER :: sourcesMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatricesSourcesType), POINTER :: sourceVectors
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatrixType), POINTER :: linearMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField,geometricField,independentField,materialsField,sourceField
    TYPE(FieldInterpolationParametersType), POINTER :: colsInterpParameters,dataInterpParameters,dataWeightInterpParameters, &
      & dependentInterpParameters,independentInterpParameters,initialGeometricInterpParameters,geometricInterpParameters, &
      & materialsInterpParameters,rowsInterpParameters,sourceInterpParameters
    TYPE(FieldInterpolatedPointType), POINTER :: dataInterpPoint,dataWeightInterpPoint,dependentInterpPoint, &
      & independentInterpPoint,initialGeometricInterpPoint,geometricInterpPoint,materialsInterpPoint,sourceInterpPoint
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: initialGeometricInterpPointMetrics,geometricInterpPointMetrics
    TYPE(FieldVariableType), POINTER :: colsVariable,dataVariable,dataWeightVariable,dependentVariable,independentVariable, &
      & geometricVariable,materialsVariable,rowsVariable,sourceVariable
    TYPE(QuadratureSchemeType), POINTER :: dependentQuadratureScheme,geometricQuadratureScheme
    TYPE(QuadratureSchemePtrType) :: columnQuadratureScheme(MAXIMUM_DATA_COMPONENTS),rowQuadratureScheme(MAXIMUM_DATA_COMPONENTS)
    TYPE(VARYING_STRING) :: localError

    ENTERS("Fitting_FiniteElementCalculate",err,error,*999)

    dataPointVector = 0.0_DP
    dataPointWeight = 0.0_DP

    CALL EquationsSet_SpecificationGet(equationsSet,4,esSpecification,err,error,*999)

    SELECT CASE(esSpecification(2))
    CASE(EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE)
      SELECT CASE(esSpecification(3))
      CASE(EQUATIONS_SET_GENERALISED_DATA_FITTING_SUBTYPE)
        !OK
      CASE(EQUATIONS_SET_DIFFUSION_TENSOR_FIBRE_DATA_FITTING_SUBTYPE)
        CALL FlagError("Diffusion tensor fibre data fitting is nonlinear.",err,error,*999)
      CASE DEFAULT
        localError="Equations set subtype "//TRIM(NumberToVString(equationsSet%specification(3),"*",err,error))// &
          & " is not valid for a data fitting equations type."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_GAUSS_FITTING_EQUATION_TYPE)
      SELECT CASE(esSpecification(3))
      CASE(EQUATIONS_SET_GENERALISED_GAUSS_FITTING_SUBTYPE)
        !OK
      CASE(EQUATIONS_SET_MAT_PROPERTIES_GAUSS_FITTING_SUBTYPE, &
        & EQUATIONS_SET_MAT_PROPERTIES_INRIA_GAUSS_FITTING_SUBTYPE)
        !OK
      CASE(EQUATIONS_SET_DIV_FREE_GAUSS_FITTING_SUBTYPE)
        !OK
      CASE DEFAULT
        localError="Equations set subtype "//TRIM(NumberToVString(equationsSet%specification(3),"*",err,error))// &
          & " is not valid for a Gauss fitting equations type."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="Equations set type "//TRIM(NumberToVString(esSpecification(2),"*",err,error))// &
          & " is not valid for a fitting equations class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
        
    smoothingType=esSpecification(4)
    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMapping)
    CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
    NULLIFY(linearMatrices)
    CALL EquationsMatricesVector_LinearMatricesGet(vectorMatrices,linearMatrices,err,error,*999)
    NULLIFY(linearMatrix)
    CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,1,linearMatrix,err,error,*999)
    CALL EquationsMatrix_UpdateMatrixGet(linearMatrix,updateMatrix,err,error,*999)
    NULLIFY(rhsMapping)
    CALL EquationsMappingVector_RHSMappingExists(vectorMapping,rhsMapping,err,error,*999)
    NULLIFY(rhsVector)
    updateRHS=.FALSE.
    IF(ASSOCIATED(rhsMapping)) THEN
      CALL EquationsMatricesVector_RHSVectorGet(vectorMatrices,rhsVector,err,error,*999)
      CALL EquationsMatricesRHS_UpdateVectorGet(rhsVector,updateRHS,err,error,*999)
    ENDIF
    NULLIFY(sourcesMapping)
    CALL EquationsMappingVector_SourcesMappingExists(vectorMapping,sourcesMapping,err,error,*999)
    NULLIFY(sourceVectors)
    NULLIFY(sourceVector)
    updateSource=.FALSE.
    IF(ASSOCIATED(sourcesMapping)) THEN
      CALL EquationsMatricesVector_SourceVectorsGet(vectorMatrices,sourceVectors,err,error,*999)
      CALL EquationsMatricesSources_SourceVectorGet(sourceVectors,1,sourceVector,err,error,*999)
      CALL EquationsMatricesSource_UpdateVectorGet(sourceVector,updateSource,err,error,*999)
    ENDIF

    update=(updateMatrix.OR.updateRHS.OR.updateSource)

    IF(update) THEN
    
      NULLIFY(lhsMapping)
      CALL EquationsMappingVector_LHSMappingGet(vectorMapping,lhsMapping,err,error,*999)
      
      NULLIFY(linearMapping)
      CALL EquationsMappingVector_LinearMappingGet(vectorMapping,linearMapping,err,error,*999)
      
      NULLIFY(equationsInterpolation)
      CALL Equations_InterpolationGet(equations,equationsInterpolation,err,error,*999)

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
      CALL Basis_NumberOfXiGet(geometricBasis,numberOfXi,err,error,*999)
      NULLIFY(geometricInterpParameters)
      CALL EquationsInterpolation_GeometricParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & geometricInterpParameters,err,error,*999)
      NULLIFY(geometricInterpPoint)
      CALL EquationsInterpolation_GeometricPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & geometricInterpPoint,err,error,*999)
      NULLIFY(geometricInterpPointMetrics)
      CALL EquationsInterpolation_GeometricPointMetricsGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & geometricInterpPointMetrics,err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,geometricInterpParameters,err,error,*999)
      NULLIFY(initialGeometricInterpParameters)
      NULLIFY(initialGeometricInterpPoint)
      NULLIFY(initialGeometricInterpPointMetrics)

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
      NULLIFY(dependentQuadratureScheme)
      CALL Basis_QuadratureSchemeGet(dependentBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,dependentQuadratureScheme,err,error,*999)
      CALL BasisQuadratureScheme_NumberOfGaussGet(dependentQuadratureScheme,numberOfGauss,err,error,*999)
      NULLIFY(dependentInterpParameters)
      CALL EquationsInterpolation_DependentParametersGet(equationsInterpolation,colsVariableType, &
        & dependentInterpParameters,err,error,*999)
      NULLIFY(dependentInterpPoint)
      CALL EquationsInterpolation_DependentPointGet(equationsInterpolation,colsVariableType,dependentInterpPoint,err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,dependentInterpParameters,err,error,*999)
       
      NULLIFY(rowsVariable)
      CALL EquationsMappingLHS_LHSVariableGet(lhsMapping,rowsVariable,err,error,*999)
      CALL FieldVariable_VariableTypeGet(rowsVariable,rowsVariableType,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(rowsVariable,numberOfRowsComponents,err,error,*999)
      
      NULLIFY(colsVariable)
      CALL EquationsMappingLinear_LinearMatrixVariableGet(linearMapping,1,colsVariable,err,error,*999)
      CALL FieldVariable_VariableTypeGet(colsVariable,colsVariableType,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(colsVariable,numberOfColsComponents,err,error,*999)
      
      !Cache row and column bases and quadrature schemes to avoid repeated calculations
      IF(numberOfRowsComponents>MAXIMUM_DATA_COMPONENTS) THEN
        localError="The number of rows components of "//TRIM(NumberToVString(numberOfRowsComponents,"*",err,error))// &
          & " is greater than the maximum number of components of "//TRIM(NumberToVString(MAXIMUM_DATA_COMPONENTS,"*", &
          & err,error))//". Increase MAXIMUM_DATA_COMPONENTS."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      DO rowComponentIdx=1,numberOfRowsComponents
        NULLIFY(rowDomain)
        CALL FieldVariable_ComponentDomainGet(rowsVariable,rowComponentIdx,rowDomain,err,error,*999)
        NULLIFY(rowDomainTopology)
        CALL Domain_DomainTopologyGet(rowDomain,rowDomainTopology,err,error,*999)
        NULLIFY(rowDomainElements)
        CALL DomainTopology_DomainElementsGet(rowDomainTopology,rowDomainElements,err,error,*999)
        NULLIFY(rowBasis(rowComponentIdx)%ptr)
        CALL DomainElements_ElementBasisGet(rowDomainElements,elementNumber,rowBasis(rowComponentIdx)%ptr,err,error,*999)
        CALL Basis_NumberOfElementParametersGet(rowBasis(rowComponentIdx)%ptr,numberOfRowElementParameters(rowComponentIdx), &
          & err,error,*999)
        NULLIFY(rowQuadratureScheme(rowComponentIdx)%ptr)
        CALL Basis_QuadratureSchemeGet(rowBasis(rowComponentIdx)%ptr,BASIS_DEFAULT_QUADRATURE_SCHEME, &
          & rowQuadratureScheme(rowComponentIdx)%ptr,err,error,*999)
      ENDDO !rowComponentIdx
      IF(updateMatrix) THEN
        IF(numberOfColsComponents>MAXIMUM_DATA_COMPONENTS) THEN
          localError="The number of columns components of "//TRIM(NumberToVString(numberOfColsComponents,"*",err,error))// &
            & " is greater than the maximum number of components of "//TRIM(NumberToVString(MAXIMUM_DATA_COMPONENTS,"*", &
            & err,error))//". Increase MAXIMUM_DATA_COMPONENTS."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        DO columnComponentIdx=1,numberOfColsComponents
          NULLIFY(columnDomain)
          CALL FieldVariable_ComponentDomainGet(colsVariable,columnComponentIdx,columnDomain,err,error,*999)
          NULLIFY(columnDomainTopology)
          CALL Domain_DomainTopologyGet(columnDomain,columnDomainTopology,err,error,*999)
          NULLIFY(columnDomainElements)
          CALL DomainTopology_DomainElementsGet(columnDomainTopology,columnDomainElements,err,error,*999)
          NULLIFY(columnBasis(columnComponentIdx)%ptr)
          CALL DomainElements_ElementBasisGet(columnDomainElements,elementNumber,columnBasis(columnComponentIdx)%ptr,err,error,*999)
          CALL Basis_NumberOfElementParametersGet(columnBasis(columnComponentIdx)%ptr, &
            & numberOfColumnElementParameters(columnComponentIdx),err,error,*999)
          NULLIFY(columnQuadratureScheme(columnComponentIdx)%ptr)
          CALL Basis_QuadratureSchemeGet(columnBasis(columnComponentIdx)%ptr,BASIS_DEFAULT_QUADRATURE_SCHEME, &
            & columnQuadratureScheme(columnComponentIdx)%ptr,err,error,*999)
        ENDDO !columnComponentIdx
      ENDIF
      
      NULLIFY(materialsField)
      NULLIFY(materialsVariable)
      NULLIFY(materialsInterpParameters)
      NULLIFY(materialsInterpPoint)
      IF(((esSpecification(2)==EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE.AND. &
        & esSpecification(3)==EQUATIONS_SET_GENERALISED_DATA_FITTING_SUBTYPE).OR. &
        & (esSpecification(2)==EQUATIONS_SET_GAUSS_FITTING_EQUATION_TYPE.AND. &
        & esSpecification(3)==EQUATIONS_SET_GENERALISED_GAUSS_FITTING_SUBTYPE).AND. &
        & (esSpecification(4)==EQUATIONS_SET_FITTING_SOBOLEV_VALUE_SMOOTHING.OR. &
        & esSpecification(4)==EQUATIONS_SET_FITTING_SOBOLEV_DIFFERENCE_SMOOTHING)).OR. &
        & (esSpecification(2)==EQUATIONS_SET_GAUSS_FITTING_EQUATION_TYPE.AND. &
        & (esSpecification(3)==EQUATIONS_SET_MAT_PROPERTIES_GAUSS_FITTING_SUBTYPE.OR. &
        & esSpecification(3)==EQUATIONS_SET_MAT_PROPERTIES_INRIA_GAUSS_FITTING_SUBTYPE))) THEN
        CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)
        CALL Field_VariableGet(materialsField,FIELD_U_VARIABLE_TYPE,materialsVariable,err,error,*999)
        CALL EquationsInterpolation_MaterialsParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
          & materialsInterpParameters,err,error,*999)
        CALL EquationsInterpolation_MaterialsPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
          & materialsInterpPoint,err,error,*999)
        CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,materialsInterpParameters, &
          & err,error,*999)
      ENDIF

      NULLIFY(independentField)
      NULLIFY(dataVariable)
      NULLIFY(dataWeightVariable)
      NULLIFY(independentVectorParameters)
      NULLIFY(independentWeightParameters)
      NULLIFY(dataInterpParameters)
      NULLIFY(dataInterpPoint)
      NULLIFY(dataWeightInterpParameters)
      NULLIFY(dataWeightInterpPoint)

      NULLIFY(sourceField)
      NULLIFY(sourceVector)
      NULLIFY(sourceInterpParameters)
      NULLIFY(sourceInterpPoint)

      SELECT CASE(esSpecification(2))
      CASE(EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE)
        
        !============================
        ! D a t a   P o i n t   F i t
        !============================
               
        CALL EquationsSet_IndependentFieldGet(equationsSet,independentField,err,error,*999)
        CALL Field_VariableGet(independentField,FIELD_U_VARIABLE_TYPE,dataVariable,err,error,*999)
        CALL FieldVariable_NumberOfComponentsGet(dataVariable,numberOfDataComponents,err,error,*999)
        CALL Field_VariableGet(independentField,FIELD_V_VARIABLE_TYPE,dataWeightVariable,err,error,*999)
        IF(numberOfDataComponents>MAXIMUM_DATA_COMPONENTS) THEN
          localError="Increase the size of the data point vectors. The data variable has "// &
            & TRIM(NumberToVString(numberOfDataComponents,"*",err,error))//" and the data point vectors only have "// &
            & TRIM(NumberToVString(MAXIMUM_DATA_COMPONENTS,"*",err,error))//" components."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        !Get data point vector parameters
        CALL FieldVariable_ParameterSetDataGet(dataVariable,FIELD_VALUES_SET_TYPE,independentVectorParameters, &
          & err,error,*999)
        !Get data point weight parameters
        CALL FieldVariable_ParameterSetDataGet(dataWeightVariable,FIELD_VALUES_SET_TYPE,independentWeightParameters, &
          & err,error,*999)
        
        SELECT CASE(esSpecification(3))
        CASE(EQUATIONS_SET_GENERALISED_DATA_FITTING_SUBTYPE)
          
          !====================================================
          ! G e n e r a l i s e d   D a t a   P o i n t   F i t
          !====================================================
          
          NULLIFY(dataProjection)
          CALL Field_DataProjectionGet(independentField,dataProjection,err,error,*999)
          NULLIFY(independentDecomposition)
          CALL Field_DecompositionGet(independentField,independentDecomposition,err,error,*999)
          NULLIFY(independentDecompositionTopology)
          CALL Decomposition_DecompositionTopologyGet(independentDecomposition,independentDecompositionTopology,err,error,*999)
          NULLIFY(decompositionDataPoints)
          CALL DecompositionTopology_DecompositionDataPointsGet(independentDecompositionTopology,decompositionDataPoints, &
            & err,error,*999)
          !Loop over data points
          CALL DecompositionDataPoints_ElementNumberOfDataPointsGet(decompositionDataPoints,elementNumber, &
            & numberOfElementDataPoints,err,error,*999)
          DO dataPointIdx=1,numberOfElementDataPoints
            CALL DecompositionDataPoints_ElementDataNumbersGet(decompositionDataPoints,dataPointIdx,elementNumber, &
              & dataPointLocalNumber,dataPointGlobalNumber,dataPointUserNumber,err,error,*999)
            !Need to use global number to get the correct projection results
            CALL DataProjection_ResultElementXiGet(dataProjection,dataPointGlobalNumber,numberOfElementXi,elementXi,err,error,*999)
            !Get data point vector value and weight
            DO componentIdx=1,numberOfDataComponents
              CALL FieldVariable_LocalDataPointDOFGet(dataVariable,dataPointLocalNumber,componentIdx,localDOFIdx, &
                & err,error,*999)
              dataPointVector(componentIdx)=independentVectorParameters(localDOFIdx)
              CALL FieldVariable_LocalDataPointDOFGet(dataWeightVariable,dataPointLocalNumber,componentIdx,localDOFIdx, &
                & err,error,*999)
              dataPointWeight(componentIdx)=independentWeightParameters(localDOFIdx)
            ENDDO !componentIdx

            rowElementDOFIdx=0
            !Loop over element rows
            DO rowComponentIdx=1,numberOfRowsComponents
              DO rowElementParameterIdx=1,numberOfRowElementParameters(rowComponentIdx)
                rowElementDOFIdx=rowElementDOFIdx+1
                CALL Basis_EvaluateXi(rowBasis(rowComponentIdx)%ptr,rowElementParameterIdx,NO_PART_DERIV, &
                  & elementXi(1:numberOfXi),rowPhi,err,error,*999)
                IF(updateMatrix) THEN
                  columnElementDOFIdx=0
                  !Loop over element columns
                  DO columnComponentIdx=1,numberOfColsComponents
                    !Treat each component as separate and independent so only calculate the diagonal blocks
                    IF(columnComponentIdx==rowComponentIdx) THEN
                      DO columnElementParameterIdx=1,numberOfColumnElementParameters(columnComponentIdx)
                        columnElementDOFIdx=columnElementDOFIdx+1
                        CALL Basis_EvaluateXi(columnBasis(columnComponentIdx)%ptr,columnElementParameterIdx,NO_PART_DERIV, &
                          & elementXi(1:numberOfXi),columnPhi,err,error,*999)
                        sum = rowPhi*columnPhi*dataPointWeight(rowComponentIdx)
                        linearMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                          & linearMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)+sum
                      ENDDO !columnElementParameterIdx
                    ENDIF !column = row
                  ENDDO !columnComponentIdx
                ENDIF !updateMatrix
                IF(updateRHS) THEN
                  sum = rowPhi*dataPointVector(rowComponentIdx)*dataPointWeight(rowComponentIdx)
                  rhsVector%elementVector%vector(rowElementDOFIdx)= &
                    & rhsVector%elementVector%vector(rowElementDOFIdx)+sum
                ENDIF !updateRHS
              ENDDO !rowElementParameterIdx
            ENDDO !rowComponentIdx
          ENDDO !dataPointIdx

          SELECT CASE(smoothingType)
          CASE(EQUATIONS_SET_FITTING_NO_SMOOTHING)
            !Do nothing
          CASE(EQUATIONS_SET_FITTING_SOBOLEV_VALUE_SMOOTHING, &
            & EQUATIONS_SET_FITTING_SOBOLEV_DIFFERENCE_SMOOTHING)
            IF(updateMatrix) THEN
              IF(smoothingType==EQUATIONS_SET_FITTING_SOBOLEV_DIFFERENCE_SMOOTHING) THEN
                geometricDerivative=SECOND_PART_DERIV
              ELSE
                geometricDerivative=FIRST_PART_DERIV
              ENDIF
              !Loop over Gauss points
              DO gaussPointIdx=1,numberOfGauss
                !Interpolate fields
                CALL Field_InterpolateGauss(geometricDerivative,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx, &
                  & geometricInterpPoint,err,error,*999)
                CALL Field_InterpolateGauss(SECOND_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx, &
                  & dependentInterpPoint,err,error,*999)
                CALL Field_InterpolatedPointMetricsCalculate(numberOfXi, geometricInterpPointMetrics,err,error,*999)
                !Get Sobolev smoothing parameters from interpolated material field
                CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx, &
                  & materialsInterpPoint,err,error,*999)
                tau1=materialsInterpPoint%values(1,NO_PART_DERIV)
                kappa11=materialsInterpPoint%values(2,NO_PART_DERIV)
                IF(numberOfXi>1) THEN
                  tau2=materialsInterpPoint%values(3,NO_PART_DERIV)
                  kappa22=materialsInterpPoint%values(4,NO_PART_DERIV)
                  kappa12=materialsInterpPoint%values(5,NO_PART_DERIV)
                  IF(numberOfXi>2) THEN
                    tau3=materialsInterpPoint%values(6,NO_PART_DERIV)
                    kappa33=materialsInterpPoint%values(7,NO_PART_DERIV)
                    kappa13=materialsInterpPoint%values(8,NO_PART_DERIV)
                    kappa23=materialsInterpPoint%values(9,NO_PART_DERIV)
                  ENDIF
                ENDIF
               
                !Calculate Jacobian and Gauss weight.
                CALL BasisQuadratureScheme_GaussWeightGet(dependentQuadratureScheme,gaussPointIdx,gaussWeight,err,error,*999)
                CALL FieldInterpolatedPointMetrics_JacobianGet(geometricInterpPointMetrics,jacobian,err,error,*999)
                jacobianGaussWeight=jacobian*gaussWeight
                
                !Loop over field components
                rowElementDOFIdx=0
                DO rowComponentIdx=1,numberOfRowsComponents
                  !Loop over element rows
                  DO rowElementParameterIdx=1,numberOfRowElementParameters(rowComponentIdx)
                    rowElementDOFIdx=rowElementDOFIdx+1
                    CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme(rowComponentIdx)%ptr, &
                      & rowElementParameterIdx,PART_DERIV_S1,gaussPointIdx,rowdPhidXi1,err,error,*999)
                    CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme(rowComponentIdx)%ptr, &
                      & rowElementParameterIdx,PART_DERIV_S1_S1,gaussPointIdx,rowd2PhidXi1dXi1,err,error,*999)
                    IF(numberOfXi>1) THEN
                      CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme(rowComponentIdx)%ptr, &
                        & rowElementParameterIdx,PART_DERIV_S2,gaussPointIdx,rowdPhidXi2,err,error,*999)
                      CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme(rowComponentIdx)%ptr, &
                        & rowElementParameterIdx,PART_DERIV_S2_S2,gaussPointIdx,rowd2PhidXi2dXi2,err,error,*999)
                      CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme(rowComponentIdx)%ptr, &
                        & rowElementParameterIdx,PART_DERIV_S1_S2,gaussPointIdx,rowd2PhidXi1dXi2,err,error,*999)
                      IF(numberOfXi>2) THEN
                        CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme(rowComponentIdx)%ptr, &
                          & rowElementParameterIdx,PART_DERIV_S3,gaussPointIdx,rowdPhidXi3,err,error,*999)
                        CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme(rowComponentIdx)%ptr, &
                          & rowElementParameterIdx,PART_DERIV_S3_S3,gaussPointIdx,rowd2PhidXi3dXi3,err,error,*999)
                        CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme(rowComponentIdx)%ptr, &
                          & rowElementParameterIdx,PART_DERIV_S1_S3,gaussPointIdx,rowd2PhidXi1dXi3,err,error,*999)
                        CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme(rowComponentIdx)%ptr, &
                          & rowElementParameterIdx,PART_DERIV_S2_S3,gaussPointIdx,rowd2PhidXi2dXi3,err,error,*999)
                      ENDIF
                    ENDIF
                    columnElementDOFIdx=0
                    !Loop over element columns
                    DO columnComponentIdx=1,numberOfColsComponents
                      DO columnElementParameterIdx=1,numberOfColumnElementParameters(columnComponentIdx)
                        columnElementDOFIdx=columnElementDOFIdx+1
                        !Calculate Sobolev surface tension and curvature smoothing terms
                        CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme(columnComponentIdx)%ptr, &
                          & columnElementParameterIdx,PART_DERIV_S1,gaussPointIdx,columndPhidXi1,err,error,*999)
                        CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme(columnComponentIdx)%ptr, &
                          & columnElementParameterIdx,PART_DERIV_S1_S1,gaussPointIdx,columnd2PhidXi1dXi1,err,error,*999)
                        tension = tau1*rowdPhidXi1*columndPhidXi1
                        curvature = kappa11*rowd2PhidXi1dXi1*columnd2PhidXi1dXi1
                        IF(numberOfXi > 1) THEN
                          CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme(columnComponentIdx)%ptr, &
                            & columnElementParameterIdx,PART_DERIV_S2,gaussPointIdx,columndPhidXi2,err,error,*999)
                          CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme(columnComponentIdx)%ptr, &
                            & columnElementParameterIdx,PART_DERIV_S2_S2,gaussPointIdx,columnd2PhidXi2dXi2,err,error,*999)
                          CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme(columnComponentIdx)%ptr, &
                            & columnElementParameterIdx,PART_DERIV_S1_S2,gaussPointIdx,columnd2PhidXi1dXi2,err,error,*999)
                          tension = tension + tau2*rowdPhidXi2*columndPhidXi2
                          curvature = curvature + kappa22*rowd2PhidXi2dXi2*columnd2PhidXi2dXi2+ &
                            & kappa12*rowd2PhidXi1dXi2*columnd2PhidXi1dXi2
                          IF(numberOfXi > 2) THEN
                            CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme(columnComponentIdx)%ptr, &
                              & columnElementParameterIdx,PART_DERIV_S3,gaussPointIdx,columndPhidXi3,err,error,*999)
                            CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme(columnComponentIdx)%ptr, &
                              & columnElementParameterIdx,PART_DERIV_S3_S3,gaussPointIdx,columnd2PhidXi3dXi3,err,error,*999)
                            CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme(columnComponentIdx)%ptr, &
                              & columnElementParameterIdx,PART_DERIV_S1_S3,gaussPointIdx,columnd2PhidXi1dXi3,err,error,*999)
                            CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme(columnComponentIdx)%ptr, &
                              & columnElementParameterIdx,PART_DERIV_S2_S3,gaussPointIdx,columnd2PhidXi2dXi3,err,error,*999)
                            tension = tension + tau3*rowdPhidXi3*columndPhidXi3
                            curvature = curvature + kappa33*rowd2PhidXi3dXi3*columnd2PhidXi3dXi3+ &
                              & kappa13*rowd2PhidXi1dXi3*columnd2PhidXi1dXi3+kappa23*rowd2PhidXi2dXi3*columnd2PhidXi2dXi3
                          ENDIF ! 3D
                        ENDIF ! 2 or 3D
                        sum = 2.0_DP*(tension + curvature)*jacobianGaussWeight
                                                
                        linearMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                          linearMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)+sum
                        
                      ENDDO !columnElementParameterIdx
                    ENDDO !columnComponentIdx
                    IF(smoothingType==EQUATIONS_SET_FITTING_SOBOLEV_DIFFERENCE_SMOOTHING) THEN
                      IF(updateRHS) THEN
                        !Calculate Sobolev surface tension and curvature smoothing terms
                        rhsTension=0.0_DP
                        rhsCurvature=0.0_DP
                        DO columnComponentIdx=1,numberOfColsComponents
                          rhsTension = rhsTension+ &
                            & tau1*rowdPhidXi1*geometricInterpPoint%values(columnComponentIdx,PART_DERIV_S1)
                          rhsCurvature = rhsCurvature + &
                            & kappa11*rowd2PhidXi1dXi1*geometricInterpPoint%values(columnComponentIdx,PART_DERIV_S1_S1)
                          IF(numberOfXi > 1) THEN
                            rhsTension = rhsTension + &
                              & tau2*rowdPhidXi2*geometricInterpPoint%values(columnComponentIdx,PART_DERIV_S2)
                            rhsCurvature = rhsCurvature + &
                              & kappa22*rowd2PhidXi2dXi2*geometricInterpPoint%values(columnComponentIdx,PART_DERIV_S2_S2) + &
                              & kappa12*rowd2PhidXi1dXi2*geometricInterpPoint%values(columnComponentIdx,PART_DERIV_S1_S2) 
                            IF(numberOfXi > 2) THEN
                              rhsTension = rhsTension + &
                                & tau3*rowdPhidXi3*geometricInterpPoint%values(columnComponentIdx,PART_DERIV_S3) 
                              rhsCurvature = rhsCurvature + &
                                & kappa33*rowd2PhidXi3dXi3*geometricInterpPoint%values(columnComponentIdx,PART_DERIV_S3_S3) + &
                                & kappa13*rowd2PhidXi1dXi3*geometricInterpPoint%values(columnComponentIdx,PART_DERIV_S1_S3) + &
                                & kappa23*rowd2PhidXi2dXi3*geometricInterpPoint%values(columnComponentIdx,PART_DERIV_S2_S3)
                            ENDIF ! 3D
                          ENDIF ! 2 or 3D
                        ENDDO !columnComponentIdx
                        sum = 2.0_DP*(rhsTension + rhsCurvature)*jacobianGaussWeight
                        rhsVector%elementVector%vector(rowElementDOFIdx)= &
                          & rhsVector%elementVector%vector(rowElementDOFIdx)+sum
                      ENDIF !update RHS
                    ENDIF !Sobolev difference smoothing                                     
                  ENDDO !rowElementParameterIdx
                ENDDO !rowComponentIdx
              ENDDO !gaussPointIdx
            ENDIF !updateMatrix
            
          CASE(EQUATIONS_SET_FITTING_STRAIN_ENERGY_SMOOTHING)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The fitting smoothing type of "//TRIM(NumberToVString(smoothingType,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
          
        CASE DEFAULT
          localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
            & " is not valid for a data fitting equations set class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        
        !Restore data point vector parameters
        CALL FieldVariable_ParameterSetDataRestore(dataVariable,FIELD_VALUES_SET_TYPE,independentVectorParameters,err,error,*999)
        !Restore data point weight parameters
        CALL FieldVariable_ParameterSetDataRestore(dataWeightVariable,FIELD_VALUES_SET_TYPE,independentWeightParameters, &
          & err,error,*999)
        
      CASE(EQUATIONS_SET_GAUSS_FITTING_EQUATION_TYPE)
        
        !==============================
        ! G a u s s   P o i n t   F i t
        !==============================

        NULLIFY(independentField)
        CALL EquationsSet_IndependentFieldGet(equationsSet,independentField,err,error,*999)
        NULLIFY(dataVariable)
        CALL Field_VariableGet(independentField,FIELD_U_VARIABLE_TYPE,dataVariable,err,error,*999)
        CALL FieldVariable_NumberOfComponentsGet(independentVariable,numberOfDataComponents,err,error,*999)
        NULLIFY(dataWeightVariable)
        CALL Field_VariableGet(independentField,FIELD_V_VARIABLE_TYPE,dataWeightVariable,err,error,*999)
        NULLIFY(dataInterpParameters)
        CALL EquationsInterpolation_IndependentParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
          & dataInterpParameters,err,error,*999)
        NULLIFY(dataInterpPoint)
        CALL EquationsInterpolation_IndependentPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
          & dataInterpPoint,err,error,*999)
        CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,dataInterpParameters,err,error,*999)
        NULLIFY(dataWeightInterpParameters)
        CALL EquationsInterpolation_IndependentParametersGet(equationsInterpolation,FIELD_V_VARIABLE_TYPE, &
          & dataWeightInterpParameters,err,error,*999)
        NULLIFY(dataWeightInterpParameters)
        CALL EquationsInterpolation_IndependentPointGet(equationsInterpolation,FIELD_V_VARIABLE_TYPE, &
          & dataWeightInterpPoint,err,error,*999)
        CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,dataWeightInterpParameters, &
          & err,error,*999)
        
        SELECT CASE(esSpecification(3))
        CASE(EQUATIONS_SET_GENERALISED_GAUSS_FITTING_SUBTYPE)
          
          !======================================================
          ! G e n e r a l i s e d   G a u s s   P o i n t   F i t
          !======================================================
            
          !Loop over Gauss points
          DO gaussPointIdx=1,numberOfGauss
            !Interpolate fields
            CALL Field_InterpolateGauss(SECOND_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,geometricInterpPoint, &
              & err,error,*999)
            CALL Field_InterpolatedPointMetricsCalculate(numberOfXi,geometricInterpPointMetrics,err,error,*999)
            CALL Field_InterpolateGauss(SECOND_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,dependentInterpPoint, &
              & err,error,*999)
            CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,dataInterpPoint, &
              & err,error,*999)
            CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,dataWeightInterpPoint, &
              & err,error,*999)
            
            !Get fitting data from interpolated fields
            DO componentIdx=1,numberOfDataComponents
              dataPointVector(componentIdx)=dataInterpPoint%values(componentIdx,NO_PART_DERIV)
              dataPointWeight(componentIdx)=dataWeightInterpPoint%values(componentIdx,NO_PART_DERIV)
            ENDDO !componentIdx
            
            SELECT CASE(smoothingType)
            CASE(EQUATIONS_SET_FITTING_NO_SMOOTHING)
              !Do nothing
            CASE(EQUATIONS_SET_FITTING_SOBOLEV_VALUE_SMOOTHING, &
              & EQUATIONS_SET_FITTING_SOBOLEV_DIFFERENCE_SMOOTHING)
              !Get Sobolev smoothing data from interpolated fields
              CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,materialsInterpPoint, &
                & err,error,*999)
              tau1=materialsInterpPoint%values(1,NO_PART_DERIV)
              kappa11=materialsInterpPoint%values(2,NO_PART_DERIV)
              IF(numberOfXi>1) THEN
                tau2=materialsInterpPoint%values(3,NO_PART_DERIV)
                kappa22=materialsInterpPoint%values(4,NO_PART_DERIV)
                kappa12=materialsInterpPoint%values(5,NO_PART_DERIV)
                IF(numberOfXi>2) THEN
                  tau3=materialsInterpPoint%values(6,NO_PART_DERIV)
                  kappa33=materialsInterpPoint%values(7,NO_PART_DERIV)
                  kappa13=materialsInterpPoint%values(8,NO_PART_DERIV)
                  kappa23=materialsInterpPoint%values(9,NO_PART_DERIV)
                ENDIF
              ENDIF
            CASE(EQUATIONS_SET_FITTING_STRAIN_ENERGY_SMOOTHING)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE DEFAULT
              localError="The fitting smoothing type of "//TRIM(NumberToVString(smoothingType,"*",err,error))// &
                & " is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
            
            !Calculate Jacobian and Gauss weight.
            CALL BasisQuadratureScheme_GaussWeightGet(dependentQuadratureScheme,gaussPointIdx,gaussWeight,err,error,*999)
            CALL FieldInterpolatedPointMetrics_JacobianGet(geometricInterpPointMetrics,jacobian,err,error,*999)
            jacobianGaussWeight=jacobian*gaussWeight
            
            rowElementDOFIdx=0
            !Loop over element rows
            DO rowComponentIdx=1,numberOfRowsComponents
               DO rowElementParameterIdx=1,numberOfRowElementParameters(rowComponentIdx)
                rowElementDOFIdx=rowElementDOFIdx+1
                CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme(rowComponentIdx)%ptr, &
                  & rowElementParameterIdx,NO_PART_DERIV,gaussPointIdx,rowPhi,err,error,*999)
                IF(updateMatrix) THEN
                  IF(smoothingType==EQUATIONS_SET_FITTING_SOBOLEV_VALUE_SMOOTHING.OR. &
                    & smoothingType==EQUATIONS_SET_FITTING_SOBOLEV_DIFFERENCE_SMOOTHING) THEN
                    CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme(rowComponentIdx)%ptr, &
                      & rowElementParameterIdx,PART_DERIV_S1,gaussPointIdx,rowdPhidXi1,err,error,*999)
                    CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme(rowComponentIdx)%ptr, &
                      & rowElementParameterIdx,PART_DERIV_S1_S1,gaussPointIdx,rowd2PhidXi1dXi1,err,error,*999)
                    IF(numberOfXi>1) THEN
                      CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme(rowComponentIdx)%ptr, &
                        & rowElementParameterIdx,PART_DERIV_S2,gaussPointIdx,rowdPhidXi2,err,error,*999)
                      CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme(rowComponentIdx)%ptr, &
                        & rowElementParameterIdx,PART_DERIV_S2_S2,gaussPointIdx,rowd2PhidXi2dXi2,err,error,*999)
                      CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme(rowComponentIdx)%ptr, &
                        & rowElementParameterIdx,PART_DERIV_S1_S2,gaussPointIdx,rowd2PhidXi1dXi2,err,error,*999)
                      IF(numberOfXi>2) THEN
                        CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme(rowComponentIdx)%ptr, &
                          & rowElementParameterIdx,PART_DERIV_S3,gaussPointIdx,rowdPhidXi3,err,error,*999)
                        CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme(rowComponentIdx)%ptr, &
                          & rowElementParameterIdx,PART_DERIV_S3_S3,gaussPointIdx,rowd2PhidXi3dXi3,err,error,*999)
                        CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme(rowComponentIdx)%ptr, &
                          & rowElementParameterIdx,PART_DERIV_S1_S3,gaussPointIdx,rowd2PhidXi1dXi3,err,error,*999)
                        CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme(rowComponentIdx)%ptr, &
                          & rowElementParameterIdx,PART_DERIV_S2_S3,gaussPointIdx,rowd2PhidXi2dXi3,err,error,*999)
                      ENDIF
                    ENDIF
                  ENDIF !Sobolev smoothing
                  columnElementDOFIdx=0
                  !Loop over element columns
                  DO columnComponentIdx=1,numberOfColsComponents
                    !Treat each component as separate and independent so only calculate the diagonal blocks
                    IF(columnComponentIdx==rowComponentIdx) THEN
                      DO columnElementParameterIdx=1,numberOfColumnElementParameters(columnComponentIdx)
                        columnElementDOFIdx=columnElementDOFIdx+1
                        CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme(columnComponentIdx)%ptr, &
                          & columnElementParameterIdx,NO_PART_DERIV,gaussPointIdx,columnPhi,err,error,*999)
                        sum=rowPhi*columnPhi*dataPointWeight(rowComponentIdx)
                        SELECT CASE(smoothingType)
                        CASE(EQUATIONS_SET_FITTING_NO_SMOOTHING)
                          !Do nothing
                        CASE(EQUATIONS_SET_FITTING_SOBOLEV_VALUE_SMOOTHING)
                          !Calculate Sobolev surface tension and curvature smoothing terms
                          CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme(columnComponentIdx)%ptr, &
                            & columnElementParameterIdx,PART_DERIV_S1,gaussPointIdx,columndPhidXi1,err,error,*999)
                          CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme(columnComponentIdx)%ptr, &
                            & columnElementParameterIdx,PART_DERIV_S1_S1,gaussPointIdx,columnd2PhidXi1dXi1,err,error,*999)
                          tension = tau1*rowdPhidXi1*columndPhidXi1
                          curvature =  kappa11*rowd2PhidXi1dXi1*columnd2PhidXi1dXi1
                          IF(numberOfXi > 1) THEN
                            CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme(columnComponentIdx)%ptr, &
                              & columnElementParameterIdx,PART_DERIV_S2,gaussPointIdx,columndPhidXi2,err,error,*999)
                            CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme(columnComponentIdx)%ptr, &
                              & columnElementParameterIdx,PART_DERIV_S2_S2,gaussPointIdx,columnd2PhidXi2dXi2,err,error,*999)
                            CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme(columnComponentIdx)%ptr, &
                              & columnElementParameterIdx,PART_DERIV_S1_S2,gaussPointIdx,columnd2PhidXi1dXi2,err,error,*999)
                            tension = tension + tau2*rowdPhidXi2*columndPhidXi2
                            curvature = curvature + kappa22*rowd2PhidXi2dXi2*columnd2PhidXi2dXi2+ &
                              & kappa12*rowd2PhidXi1dXi2*columnd2PhidXi1dXi2
                            IF(numberOfXi > 2) THEN
                              CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme(columnComponentIdx)%ptr, &
                                & columnElementParameterIdx,PART_DERIV_S3,gaussPointIdx,columndPhidXi3,err,error,*999)
                              CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme(columnComponentIdx)%ptr, &
                                & columnElementParameterIdx,PART_DERIV_S3_S3,gaussPointIdx,columnd2PhidXi3dXi3,err,error,*999)
                              CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme(columnComponentIdx)%ptr, &
                                & columnElementParameterIdx,PART_DERIV_S1_S3,gaussPointIdx,columnd2PhidXi1dXi3,err,error,*999)
                              CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme(columnComponentIdx)%ptr, &
                                & columnElementParameterIdx,PART_DERIV_S2_S3,gaussPointIdx,columnd2PhidXi2dXi3,err,error,*999)
                              tension = tension + tau3*rowdPhidXi3*columndPhidXi3
                              curvature = curvature + kappa33*rowd2PhidXi3dXi3*columnd2PhidXi3dXi3+ &
                                & kappa13*rowd2PhidXi1dXi3*columnd2PhidXi1dXi3+kappa23*rowd2PhidXi2dXi3*columnd2PhidXi2dXi3
                            ENDIF ! 3D
                          ENDIF ! 2 or 3D
                          sum = sum + 2.0_DP*(tension + curvature)*jacobianGaussWeight
                        CASE(EQUATIONS_SET_FITTING_SOBOLEV_DIFFERENCE_SMOOTHING)
                          CALL FlagError("Not implemented.",err,error,*999)
                        CASE(EQUATIONS_SET_FITTING_STRAIN_ENERGY_SMOOTHING)
                          CALL FlagError("Not implemented.",err,error,*999)
                        CASE DEFAULT
                          localError="The fitting smoothing type of "//TRIM(NumberToVString(smoothingType,"*",err,error))// &
                            & " is invalid."
                          CALL FlagError(localError,err,error,*999)
                        END SELECT
                        
                        linearMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                          linearMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)+sum
                      ENDDO !columnElementParameterIdx
                    ENDIF
                  ENDDO !columnComponentIdx
                ENDIF
                IF(updateRHS) THEN
                  rhsVector%elementVector%vector(rowElementDOFIdx)= &
                    & rhsVector%elementVector%vector(rowElementDOFIdx) + &
                    & rowPhi*dataPointVector(rowComponentIdx)*dataPointWeight(rowComponentIdx)
                ENDIF
              ENDDO !rowElementParameterIdx
            ENDDO !rowComponentIdx
          ENDDO !gaussPointIdx
          
        CASE(EQUATIONS_SET_MAT_PROPERTIES_GAUSS_FITTING_SUBTYPE, &
          & EQUATIONS_SET_MAT_PROPERTIES_INRIA_GAUSS_FITTING_SUBTYPE)
          
          !========================================================================
          ! M a t e r i a l s   P r o p e r t i e s   G a u s s   P o i n t   F i t
          !========================================================================
          
          CALL FieldVariable_InterpolationParameterInitialise(geometricVariable,initialGeometricInterpParameters, &
            & err,error,*999)
          CALL Field_InterpolatedPointInitialise(initialGeometricInterpParameters,initialGeometricInterpPoint, &
            & err,error,*999)
          CALL Field_InterpolatedPointMetricsInitialise(initialGeometricInterpPoint,initialGeometricInterpPointMetrics, &
            & err,error,*999)
          CALL Field_InterpolationParametersElementGet(FIELD_INITIAL_VALUES_SET_TYPE,elementNumber, &
            & initialGeometricInterpParameters,err,error,*999)
          
          !--- Loop over gauss points
          DO gaussPointIdx=1,numberOfGauss
            
            !--- Interpolation of (actual) Geometry and Metrics
            CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx, &
              & geometricInterpPoint,err,error,*999)
            CALL Field_InterpolatedPointMetricsCalculate(numberOfXi,geometricInterpPointMetrics,err,error,*999)
            !--- Interpolation of Reference Geometry
            CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx, &
              & initialGeometricInterpPoint,err,error,*999)
            !--- Retrieve local dXdXi and dYdXi
            DO xiIdx=1,numberOfXi
              derivativeIdx=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xiIdx)
              DO componentIdx=1,numberOfDimensions
                !dy/dxi (y = referential)
                dYdXi(componentIdx,xiIdx)=initialGeometricInterpPoint%values(componentIdx,derivativeIdx)
                !dx/dxi
                dXdXi(componentIdx,xiIdx)=geometricInterpPoint%values(componentIdx,derivativeIdx)
              ENDDO !componentIdx
            ENDDO !xiIdx
            
            !--- Compute deformation gradient tensor dXdY and its Jacobian Jxy
            CALL Invert(dYdXi,dXidY,Jyxi,err,error,*999) !dy/dxi -> dxi/dy
            CALL MatrixProduct(dXdXi,dXidY,dXdY,err,error,*999) !dx/dxi * dxi/dy = dx/dy (deformation gradient tensor, F)
            CALL Determinant(dXdY,Jxy,err,error,*999)
            IF(ABS(Jxy)<=ZERO_TOLERANCE) THEN
              localError="Jacobian Jxy="//TRIM(NumberToVString(Jxy,"*",err,error))//" is <= 0.0 for Gauss point "// &
                & TRIM(NumberToVString(gaussPointIdx,"*",err,error))//" of element number "// &
                & TRIM(NumberToVString(elementNumber,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
            
            !--- Interpolation of Materials Field
            CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx, &
              & materialsInterpPoint,err,error,*999)
            
            !--- Retrieve reference material parameters:
            porosity0=materialsInterpPoint%values(1,NO_PART_DERIV)
            permOverVisParam0=materialsInterpPoint%values(2,NO_PART_DERIV)
            
            !--- Material dependence on structural deformation
            porosity = 1.0_DP - ( 1.0_DP - porosity0 ) / Jxy
            
            IF(esSpecification(3)==EQUATIONS_SET_MAT_PROPERTIES_INRIA_GAUSS_FITTING_SUBTYPE) THEN
              permOverVisParam = permOverVisParam0
            ELSE
              !material modeling could use gradient information, or solve some PDE
              materialFact = ( Jxy * porosity / porosity0 )**2.0_DP
              permOverVisParam = materialFact * permOverVisParam0
            END IF
            
            !Calculate Jacobian and Gauss weight.
            CALL BasisQuadratureScheme_GaussWeightGet(dependentQuadratureScheme,gaussPointIdx,gaussWeight,err,error,*999)
            CALL FieldInterpolatedPointMetrics_JacobianGet(geometricInterpPointMetrics,jacobian,err,error,*999)
            jacobianGaussWeight=jacobian*gaussWeight
            
            IF(diagnostics2) THEN
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Gauss point : ",gaussPointIdx,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Jacobian = ",jacobian,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Jxy = ",Jxy,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  porosity = ",porosity,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  permOverVisParam = ",permOverVisParam,err,error,*999)
              CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE," ",err,error,*999)
            ENDIF
            
            !Loop over field components
            rowElementDOFIdx=0
            DO rowComponentIdx=1,numberOfRowsComponents
              !Loop over element rows
              DO rowElementParameterIdx=1,numberOfRowElementParameters(rowComponentIdx)
                rowElementDOFIdx=rowElementDOFIdx+1
                CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme(rowComponentIdx)%ptr,rowElementParameterIdx, &
                  & NO_PART_DERIV,gaussPointIdx,rowPhi,err,error,*999)
                IF(updateMatrix) THEN
                  columnElementDOFIdx=0
                  !Loop over element columns
                  DO columnComponentIdx=1,numberOfColsComponents
                    IF(rowComponentIdx==columnComponentIdx) THEN
                      DO columnElementParameterIdx=1,numberOfColumnElementParameters(columnComponentIdx)
                        columnElementDOFIdx=columnElementDOFIdx+1
                        CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme(columnComponentIdx)%ptr, &
                          & columnElementParameterIdx,NO_PART_DERIV,gaussPointIdx,columnPhi,err,error,*999)
                        sum=rowPhi*columnPhi
                        linearMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) = &
                          & linearMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) + &
                          & sum*jacobianGaussWeight                         
                      ENDDO !columnElementParameterIdx
                    ENDIF !row = column
                  ENDDO !columnComponentIdx
                ENDIF !updateMatrix
                IF(updateRHS) THEN
                  sum=0.0_DP
                  IF(rowComponentIdx==1) THEN
                    sum=sum+rowPhi*porosity
                  ELSE IF(rowComponentIdx==2) THEN
                    sum=sum+rowPhi*permOverVisParam
                  ENDIF
                  rhsVector%elementVector%vector(rowElementDOFIdx)= &
                    & rhsVector%elementVector%vector(rowElementDOFIdx)+sum*jacobianGaussWeight
                ENDIF !update RHS
              ENDDO !rowElementParameterIdx
            ENDDO !rowComponentIdx
          ENDDO !gaussPointIdx
          
          CALL Field_InterpolatedPointMetricsFinalise(initialGeometricInterpPointMetrics,err,error,*999)
          CALL Field_InterpolatedPointFinalise(initialGeometricInterpPoint,err,error,*999)
          CALL FieldVariable_InterpolationParameterFinalise(initialGeometricInterpParameters,err,error,*999)
          
        CASE(EQUATIONS_SET_DIV_FREE_GAUSS_FITTING_SUBTYPE)
          
          !==============================================================
          ! D i v e r g e n c e   F r e e   G a u s s   P o i n t   F i t
          !==============================================================
          
          CALL EquationsSet_MaterialsFieldGet(equationsSet,sourceField,err,error,*999)
          CALL Field_VariableGet(sourceField,FIELD_U_VARIABLE_TYPE,sourceVariable,err,error,*999)
          CALL EquationsInterpolation_SourceParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
            & sourceInterpParameters,err,error,*999)
          CALL EquationsInterpolation_SourcePointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
            & sourceInterpPoint,err,error,*999)
          CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,sourceInterpParameters, &
            & err,error,*999)

          !Loop over gauss points
          DO gaussPointIdx=1,numberOfGauss
            
            CALL Field_InterpolateGauss(SECOND_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx, &
              & geometricInterpPoint,err,error,*999)
            CALL Field_InterpolatedPointMetricsCalculate(numberOfXi,geometricInterpPointMetrics,err,error,*999)
            CALL Field_InterpolateGauss(SECOND_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx, &
              & dependentInterpPoint,err,error,*999)
            
            DO componentIdx=1,numberOfDimensions
              DO xiIdx=1,numberOfXi
                dXidX(xiIdx,componentIdx)=geometricInterpPointMetrics%dXidX(xiIdx,componentIdx)
              ENDDO !xiIdx
            ENDDO !componentIdx
            
            IF(smoothingType==EQUATIONS_SET_FITTING_SOBOLEV_VALUE_SMOOTHING.OR. &
              & smoothingType==EQUATIONS_SET_FITTING_SOBOLEV_DIFFERENCE_SMOOTHING) THEN
              CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx, &
                & materialsInterpPoint,err,error,*999)
              tau1=materialsInterpPoint%values(1,NO_PART_DERIV)
              kappa11=materialsInterpPoint%values(2,NO_PART_DERIV)
              IF(numberOfXi>1) THEN
                tau2=materialsInterpPoint%values(3,NO_PART_DERIV)
                kappa22=materialsInterpPoint%values(4,NO_PART_DERIV)
                kappa12=materialsInterpPoint%values(5,NO_PART_DERIV)
                IF(numberOfXi>2) THEN
                  tau3=materialsInterpPoint%values(6,NO_PART_DERIV)
                  kappa33=materialsInterpPoint%values(7,NO_PART_DERIV)
                  kappa13=materialsInterpPoint%values(8,NO_PART_DERIV)
                  kappa23=materialsInterpPoint%values(9,NO_PART_DERIV)
                ENDIF
              ENDIF
            ENDIF

            uValue=0.0_DP
            IF(updateSource) THEN
              CALL Field_InterpolateGauss(SECOND_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx, &
                & sourceInterpPoint,err,error,*999)
              uValue(1:numberOfDimensions)=sourceInterpPoint%values(1:numberOfDimensions,NO_PART_DERIV)
            ENDIF
            
            !Calculate Jacobian and Gauss weight.
            CALL BasisQuadratureScheme_GaussWeightGet(dependentQuadratureScheme,gaussPointIdx,gaussWeight,err,error,*999)
            CALL FieldInterpolatedPointMetrics_JacobianGet(geometricInterpPointMetrics,jacobian,err,error,*999)
            jacobianGaussWeight=jacobian*gaussWeight
            
            !Loop over field components
            rowElementDOFIdx=0
            DO rowComponentIdx=1,numberOfRowsComponents
              DO rowElementParameterIdx=1,numberOfRowElementParameters(rowComponentIdx)
                rowElementDOFIdx=rowElementDOFIdx+1
                CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme(rowComponentIdx)%ptr, &
                  & rowElementParameterIdx,NO_PART_DERIV,gaussPointIdx,rowPhi,err,error,*999)
                DO rowXiIdx=1,numberOfXi
                  CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme(rowComponentIdx)%ptr, &
                    & rowElementParameterIdx,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(rowXiIdx),gaussPointIdx, &
                    & rowdPhidXi(rowXiIdx),err,error,*999)
                ENDDO !rowXiIdx
                IF(updateMatrix) THEN
                  !Loop over element columns
                  columnElementDOFIdx=0
                  DO columnComponentIdx=1,numberOfColsComponents
                    DO columnElementParameterIdx=1,numberOfColumnElementParameters(columnComponentIdx)
                      columnElementDOFIdx=columnElementDOFIdx+1
                      CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme(columnComponentIdx)%ptr, &
                        & columnElementParameterIdx,NO_PART_DERIV,gaussPointIdx,columnPhi,err,error,*999)
                      DO columnXiIdx=1,numberOfXi
                        CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme(rowComponentIdx)%ptr, &
                          & columnElementParameterIdx,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(columnXiIdx),gaussPointIdx, &
                          & columndPhidXi(columnXiIdx),err,error,*999)
                      ENDDO !columnXiIdx
                      sum = 0.0_DP
                      !Calculate sum
                      IF(rowComponentIdx==columnComponentIdx.AND.rowComponentIdx<=numberOfDimensions) &
                        & sum=sum+rowPhi*columnPhi
                      
                      SELECT CASE(smoothingType)
                      CASE(EQUATIONS_SET_FITTING_NO_SMOOTHING)
                        !Do nothing
                      CASE(EQUATIONS_SET_FITTING_SOBOLEV_VALUE_SMOOTHING, &
                        & EQUATIONS_SET_FITTING_SOBOLEV_DIFFERENCE_SMOOTHING)
                        CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme(columnComponentIdx)%ptr, &
                          & columnElementParameterIdx,PART_DERIV_S1,gaussPointIdx,columndPhidXi1,err,error,*999)
                        CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme(columnComponentIdx)%ptr, &
                          & columnElementParameterIdx,PART_DERIV_S1_S1,gaussPointIdx,columnd2PhidXi1dXi1,err,error,*999)
                        tension = tau1*rowdPhidXi1*columndPhidXi1
                        curvature =  kappa11*rowd2PhidXi1dXi1*columnd2PhidXi1dXi1
                        IF(numberOfXi > 1) THEN
                          CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme(columnComponentIdx)%ptr, &
                            & columnElementParameterIdx,PART_DERIV_S2,gaussPointIdx,columndPhidXi2,err,error,*999)
                          CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme(columnComponentIdx)%ptr, &
                            & columnElementParameterIdx,PART_DERIV_S2_S2,gaussPointIdx,columnd2PhidXi2dXi2,err,error,*999)
                          CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme(columnComponentIdx)%ptr, &
                            & columnElementParameterIdx,PART_DERIV_S1_S2,gaussPointIdx,columnd2PhidXi1dXi2,err,error,*999)
                          tension = tension + tau2*rowdPhidXi2*columndPhidXi2
                          curvature = curvature + kappa22*rowd2PhidXi2dXi2*columnd2PhidXi2dXi2+ &
                            & kappa12*rowd2PhidXi1dXi2*columnd2PhidXi1dXi2
                          IF(numberOfXi > 2) THEN
                            CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme(columnComponentIdx)%ptr, &
                              & columnElementParameterIdx,PART_DERIV_S3,gaussPointIdx,columndPhidXi3,err,error,*999)
                            CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme(columnComponentIdx)%ptr, &
                              & columnElementParameterIdx,PART_DERIV_S3_S3,gaussPointIdx,columnd2PhidXi3dXi3,err,error,*999)
                            CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme(columnComponentIdx)%ptr, &
                              & columnElementParameterIdx,PART_DERIV_S1_S3,gaussPointIdx,columnd2PhidXi1dXi3,err,error,*999)
                            CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme(columnComponentIdx)%ptr, &
                              & columnElementParameterIdx,PART_DERIV_S2_S3,gaussPointIdx,columnd2PhidXi2dXi3,err,error,*999)
                            tension = tension + tau3*rowdPhidXi3*columndPhidXi3
                            curvature = curvature + kappa33*rowd2PhidXi3dXi3*columnd2PhidXi3dXi3+ &
                              & kappa13*rowd2PhidXi1dXi3*columnd2PhidXi1dXi3+kappa23*rowd2PhidXi2dXi3*columnd2PhidXi2dXi3
                          ENDIF ! 3D
                        ENDIF ! 2 or 3D
                        sum = sum + 2.0_DP*(tension + curvature)*jacobianGaussWeight
                      CASE(EQUATIONS_SET_FITTING_STRAIN_ENERGY_SMOOTHING)
                        CALL FlagError("Not implemented.",err,error,*999)
                      CASE DEFAULT
                        localError="The fitting smoothing type of "//TRIM(NumberToVString(smoothingType,"*",err,error))// &
                          & " is invalid."
                        CALL FlagError(localError,err,error,*999)
                      END SELECT
                      
                      linearMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) = &
                        & linearMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) + sum
                      
                      IF(columnComponentIdx==numberOfColsComponents.AND.rowComponentIdx<=numberOfDimensions) THEN
                        sum=0.0_DP
                        !Calculate sum
                        DO rowXiIdx=1,numberOfXi
                          sum=sum+columnPhi*rowdPhidXi(rowXiIdx)*dXidX(rowXiIdx,rowComponentIdx)
                        ENDDO !columnXiIdx
                        linearMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) = &
                          & linearMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) + &
                          & sum*jacobianGaussWeight
                        linearMatrix%elementMatrix%matrix(columnElementDOFIdx,rowElementDOFIdx) = &
                          & linearMatrix%elementMatrix%matrix(columnElementDOFIdx,rowElementDOFIdx) + &
                          & sum*jacobianGaussWeight
                      ENDIF
                    ENDDO !columnElementParameterIdx
                  ENDDO !columnComponentIdx
                ENDIF !update matrix
                IF(updateRHS) THEN
                  rhsVector%elementVector%vector(rowElementDOFIdx)=0.0_DP
                  IF(smoothingType==EQUATIONS_SET_FITTING_SOBOLEV_DIFFERENCE_SMOOTHING) THEN
                    !Calculate Sobolev surface tension and curvature smoothing terms
                    rhsTension=0.0_DP
                    rhsCurvature=0.0_DP
                    DO columnComponentIdx=1,numberOfColsComponents
                      rhsTension = rhsTension+ &
                        & tau1*rowdPhidXi1*geometricInterpPoint%values(columnComponentIdx,PART_DERIV_S1)
                      rhsCurvature = rhsCurvature + &
                        & kappa11*rowd2PhidXi1dXi1*geometricInterpPoint%values(columnComponentIdx,PART_DERIV_S1_S1)
                      IF(numberOfXi > 1) THEN
                        rhsTension = rhsTension + &
                          & tau2*rowdPhidXi2*geometricInterpPoint%values(columnComponentIdx,PART_DERIV_S2)
                        rhsCurvature = rhsCurvature + &
                          & kappa22*rowd2PhidXi2dXi2*geometricInterpPoint%values(columnComponentIdx,PART_DERIV_S2_S2) + &
                          & kappa12*rowd2PhidXi1dXi2*geometricInterpPoint%values(columnComponentIdx,PART_DERIV_S1_S2) 
                        IF(numberOfXi > 2) THEN
                          rhsTension = rhsTension + &
                            & tau3*rowdPhidXi3*geometricInterpPoint%values(columnComponentIdx,PART_DERIV_S3) 
                          rhsCurvature = rhsCurvature + &
                            & kappa33*rowd2PhidXi3dXi3*geometricInterpPoint%values(columnComponentIdx,PART_DERIV_S3_S3) + &
                            & kappa13*rowd2PhidXi1dXi3*geometricInterpPoint%values(columnComponentIdx,PART_DERIV_S1_S3) + &
                            & kappa23*rowd2PhidXi2dXi3*geometricInterpPoint%values(columnComponentIdx,PART_DERIV_S2_S3)
                        ENDIF ! 3D
                      ENDIF ! 2 or 3D
                    ENDDO !columnComponentIdx
                    sum = 2.0_DP*(rhsTension + rhsCurvature)*jacobianGaussWeight
                    rhsVector%elementVector%vector(rowElementDOFIdx)= &
                      & rhsVector%elementVector%vector(rowElementDOFIdx)+sum
                  ENDIF !Sobolev difference smoothing                                     
                ENDIF !update RHS
                IF(updateSource) THEN
                  IF(rowComponentIdx<=numberOfDimensions) THEN
                    sum=uValue(rowComponentIdx)*rowPhi
                    sourceVector%elementVector%vector(rowElementDOFIdx)=sourceVector%elementVector%vector(rowElementDOFIdx)+sum
                  ENDIF
                ENDIF !update source
              ENDDO !rowElementParameterIdx
            ENDDO !rowComponentIdx
          ENDDO !gaussPointIdx
          
        CASE DEFAULT
          localError="Equations set subtype "//TRIM(NumberToVString(equationsSet%specification(3),"*",err,error))// &
            & " is not valid for a Gauss fitting equations set class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        
      CASE DEFAULT
        localError="Equations set type "//TRIM(NumberToVString(equationsSet%specification(2),"*",err,error))// &
          & " is not valid for a fitting equations set class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      
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
          DO rowElementParameterIdx=1,numberOfRowElementParameters(rowComponentIdx)
            rowElementDOFIdx=rowElementDOFIdx+1                    
            columnElementDOFIdx=0
            IF(updateMatrix) THEN
              !Loop over element columns
              DO columnComponentIdx=1,numberOfColsComponents
                DO columnElementParameterIdx=1,numberOfColumnElementParameters(columnComponentIdx)
                  columnElementDOFIdx=columnElementDOFIdx+1
                  linearMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                    & linearMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)* &
                    & rowsInterpParameters%scaleFactors(rowElementParameterIdx,rowComponentIdx)* &
                    & colsInterpParameters%scaleFactors(columnElementParameterIdx,columnComponentIdx)
                ENDDO !columnElementParameterIdx
              ENDDO !columnComponentIdx
            ENDIF !update matrix
            IF(updateRHS) THEN
              rhsVector%elementVector%vector(rowElementDOFIdx)= &
                & rhsVector%elementVector%vector(rowElementDOFIdx)* &
                & rowsInterpParameters%scaleFactors(rowElementParameterIdx,rowComponentIdx)
            ENDIF !update RHS
            IF(updateSource) THEN
              sourceVector%elementVector%vector(rowElementDOFIdx)= &
                & sourceVector%elementVector%vector(rowElementDOFIdx)* &
                & rowsInterpParameters%scaleFactors(rowElementParameterIdx,rowComponentIdx)
            ENDIF !update RHS
          ENDDO !rowElementParameterIdx
        ENDDO !rowComponentIdx
      ENDIF !scale factors

    ENDIF !update
      
    EXITS("Fitting_FiniteElementCalculate")
    RETURN
999 ERRORSEXITS("Fitting_FiniteElementCalculate",err,error)
    RETURN 1

  END SUBROUTINE Fitting_FiniteElementCalculate

  !
  !================================================================================================================================
  !

  !>Evaluates the element residual and RHS vectors for the given element number for a fitting class finite element equation set.
  SUBROUTINE Fitting_FiniteElementResidualEvaluate(equationsSet,elementNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG), PARAMETER :: MAXIMUM_DATA_COMPONENTS=99
    INTEGER(INTG) ::  colsVariableType,columnComponentIdx,columnElementDOFIdx,columnElementParameterIdx,componentIdx, &
      & dataPointGlobalNumber,dataPointIdx,dataPointLocalNumber,dataPointUserNumber,dependentComponentColumnIdx, &
      & dependentComponentRowIdx,dependentElementParameterColumnIdx,dependentElementParameterRowIdx,dependentParameterColumnIdx, &
      & dependentParameterRowIdx,dependentVariableType,esSpecification(4),gaussPointIdx,geometricDerivative,localDofIdx, &
      & meshComponentRow,meshComponentColumn,numberOfColsComponents,numberOfColumnElementParameters(MAXIMUM_DATA_COMPONENTS), &
      & numberOfDataComponents,numberOfDimensions,numberOfElementDataPoints,numberOfElementXi,numberOfGauss, &
      & numberOfRowsComponents,numberOfRowElementParameters(MAXIMUM_DATA_COMPONENTS),numberOfXi,rowComponentIdx, &
      & rowElementDOFIdx,rowElementParameterIdx,rowsVariableType,scalingType,smoothingType
    REAL(DP) :: columnPhi,curvature,dataPointWeight(MAXIMUM_DATA_COMPONENTS),dataPointVector(MAXIMUM_DATA_COMPONENTS), &
      & elementXi(3),jacobianGaussWeight,kappaParam,rowPhi,sum,tauParam,tension
    REAL(DP), POINTER :: independentVectorParameters(:),independentWeightParameters(:)
    LOGICAL :: update,updateMatrix,updateResidual,updateRHS
    TYPE(BasisType), POINTER :: dependentBasis,geometricBasis,dependentBasisRow,dependentBasisColumn
    TYPE(BasisPtrType) :: columnBasis(MAXIMUM_DATA_COMPONENTS),rowBasis(MAXIMUM_DATA_COMPONENTS)
    TYPE(DataProjectionType), POINTER :: dataProjection
    TYPE(DecompositionType), POINTER :: dependentDecomposition,geometricDecomposition,independentDecomposition
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology,independentDecompositionTopology
    TYPE(DecompositionDataPointsType), POINTER :: decompositionDataPoints
    TYPE(DomainType), POINTER :: columnDomain,dependentDomain,geometricDomain,rowDomain
    TYPE(DomainElementsType), POINTER :: columnDomainElements,dependentDomainElements,geometricDomainElements,rowDomainElements
    TYPE(DomainTopologyType), POINTER :: columnDomainTopology,dependentDomainTopology,geometricDomainTopology,rowDomainTopology
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField,geometricField,independentField,materialsField
    TYPE(FieldInterpolationParametersType), POINTER :: dataInterpParameters,dataWeightInterpParameters,dependentInterpParameters, &
      & geometricInterpParameters,materialsInterpParameters,rowsInterpParameters
    TYPE(FieldInterpolatedPointType), POINTER :: dataInterpPoint,dataWeightInterpPoint,dependentInterpPoint,geometricInterpPoint, &
      & materialsInterpPoint
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: geometricInterpPointMetrics
    TYPE(FieldVariableType), POINTER :: colsVariable,dataVariable,dataWeightVariable,dependentVariable,geometricVariable, &
      & materialsVariable,rowsVariable
    TYPE(QuadratureSchemeType), POINTER :: dependentQuadratureScheme
    TYPE(QuadratureSchemePtrType) :: columnQuadratureScheme(MAXIMUM_DATA_COMPONENTS), &
      & rowQuadratureScheme(MAXIMUM_DATA_COMPONENTS)
    TYPE(VARYING_STRING) :: localError

    ENTERS("Fitting_FiniteElementResidualEvaluate",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,4,esSpecification,err,error,*999)
    
    SELECT CASE(esSpecification(2))
    CASE(EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE)
      SELECT CASE(esSpecification(3))
      CASE(EQUATIONS_SET_GENERALISED_DATA_FITTING_SUBTYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(EQUATIONS_SET_DIFFUSION_TENSOR_FIBRE_DATA_FITTING_SUBTYPE)
        !OK
      CASE DEFAULT
        localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
          & " is not valid for a data fitting equation type of a fitting equations set class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_GAUSS_FITTING_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="Equations set type "//TRIM(NumberToVString(esSpecification(2),"*",err,error))// &
        & " is not valid for a fitting equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    smoothingType=esSpecification(4)
    
    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMapping)
    CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
    NULLIFY(lhsMapping)
    CALL EquationsMappingVector_LHSMappingGet(vectorMapping,lhsMapping,err,error,*999)
    NULLIFY(nonlinearMapping)
    CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)
    NULLIFY(residualMapping)
    CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,1,residualMapping,err,error,*999)
    NULLIFY(rhsMapping)
    CALL EquationsMappingVector_RHSMappingExists(vectorMapping,rhsMapping,err,error,*999)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
    NULLIFY(nonlinearMatrices)
    CALL EquationsMatricesVector_NonlinearMatricesGet(vectorMatrices,nonlinearMatrices,err,error,*999)
    NULLIFY(residualVector)
    CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,1,residualVector,err,error,*999)
    CALL EquationsMatricesResidual_UpdateVectorGet(residualVector,updateResidual,err,error,*999)
    NULLIFY(rhsVector)
    updateRHS=.FALSE.
    IF(ASSOCIATED(rhsMapping)) THEN
      CALL EquationsMatricesVector_RHSVectorGet(vectorMatrices,rhsVector,err,error,*999)
      CALL EquationsMatricesRHS_UpdateVectorGet(rhsVector,updateRHS,err,error,*999)
    ENDIF

    update=(updateResidual.OR.updateRHS)

    IF(update) THEN
      
      NULLIFY(equationsInterpolation)
      CALL Equations_InterpolationGet(equations,equationsInterpolation,err,error,*999)

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
      CALL Basis_NumberOfXiGet(geometricBasis,numberOfXi,err,error,*999)
      NULLIFY(geometricInterpParameters)
      CALL EquationsInterpolation_GeometricParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & geometricInterpParameters,err,error,*999)
      NULLIFY(geometricInterpPoint)
      CALL EquationsInterpolation_GeometricPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & geometricInterpPoint,err,error,*999)
      NULLIFY(geometricInterpPointMetrics)
      CALL EquationsInterpolation_GeometricPointMetricsGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & geometricInterpPointMetrics,err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,geometricInterpParameters,err,error,*999)
 
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
      NULLIFY(dependentQuadratureScheme)
      CALL Basis_QuadratureSchemeGet(dependentBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,dependentQuadratureScheme,err,error,*999)
      CALL BasisQuadratureScheme_NumberOfGaussGet(dependentQuadratureScheme,numberOfGauss,err,error,*999)
      NULLIFY(dependentInterpParameters)
      CALL EquationsInterpolation_DependentParametersGet(equationsInterpolation,colsVariableType, &
        & dependentInterpParameters,err,error,*999)
      NULLIFY(dependentInterpPoint)
      CALL EquationsInterpolation_DependentPointGet(equationsInterpolation,colsVariableType,dependentInterpPoint,err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,dependentInterpParameters,err,error,*999)
    
      NULLIFY(rowsVariable)
      CALL EquationsMappingLHS_LHSVariableGet(lhsMapping,rowsVariable,err,error,*999)
      CALL FieldVariable_VariableTypeGet(rowsVariable,rowsVariableType,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(rowsVariable,numberOfRowsComponents,err,error,*999)
      
      NULLIFY(colsVariable)
      CALL EquationsMappingResidual_VariableGet(residualMapping,1,colsVariable,err,error,*999)
      CALL FieldVariable_VariableTypeGet(colsVariable,colsVariableType,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(colsVariable,numberOfColsComponents,err,error,*999)
      
      !Cache row and column bases and quadrature schemes to avoid repeated calculations
      IF(numberOfRowsComponents>MAXIMUM_DATA_COMPONENTS) THEN
        localError="The number of rows components of "//TRIM(NumberToVString(numberOfRowsComponents,"*",err,error))// &
          & " is greater than the maximum number of components of "//TRIM(NumberToVString(MAXIMUM_DATA_COMPONENTS,"*", &
          & err,error))//". Increase MAXIMUM_DATA_COMPONENTS."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      DO rowComponentIdx=1,numberOfRowsComponents
        NULLIFY(rowDomain)
        CALL FieldVariable_ComponentDomainGet(rowsVariable,rowComponentIdx,rowDomain,err,error,*999)
        NULLIFY(rowDomainTopology)
        CALL Domain_DomainTopologyGet(rowDomain,rowDomainTopology,err,error,*999)
        NULLIFY(rowDomainElements)
        CALL DomainTopology_DomainElementsGet(rowDomainTopology,rowDomainElements,err,error,*999)
        NULLIFY(rowBasis(rowComponentIdx)%ptr)
        CALL DomainElements_ElementBasisGet(rowDomainElements,elementNumber,rowBasis(rowComponentIdx)%ptr,err,error,*999)
        CALL Basis_NumberOfElementParametersGet(rowBasis(rowComponentIdx)%ptr,numberOfRowElementParameters(rowComponentIdx), &
          & err,error,*999)
        NULLIFY(rowQuadratureScheme(rowComponentIdx)%ptr)
        CALL Basis_QuadratureSchemeGet(rowBasis(rowComponentIdx)%ptr,BASIS_DEFAULT_QUADRATURE_SCHEME, &
          & rowQuadratureScheme(rowComponentIdx)%ptr,err,error,*999)
      ENDDO !rowComponentIdx
      IF(updateMatrix) THEN
        IF(numberOfColsComponents>MAXIMUM_DATA_COMPONENTS) THEN
          localError="The number of columns components of "//TRIM(NumberToVString(numberOfColsComponents,"*",err,error))// &
            & " is greater than the maximum number of components of "//TRIM(NumberToVString(MAXIMUM_DATA_COMPONENTS,"*", &
            & err,error))//". Increase MAXIMUM_DATA_COMPONENTS."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        DO columnComponentIdx=1,numberOfColsComponents
          NULLIFY(columnDomain)
          CALL FieldVariable_ComponentDomainGet(colsVariable,columnComponentIdx,columnDomain,err,error,*999)
          NULLIFY(columnDomainTopology)
          CALL Domain_DomainTopologyGet(columnDomain,columnDomainTopology,err,error,*999)
          NULLIFY(columnDomainElements)
          CALL DomainTopology_DomainElementsGet(columnDomainTopology,columnDomainElements,err,error,*999)
          NULLIFY(columnBasis(columnComponentIdx)%ptr)
          CALL DomainElements_ElementBasisGet(columnDomainElements,elementNumber,columnBasis(columnComponentIdx)%ptr,err,error,*999)
          CALL Basis_NumberOfElementParametersGet(columnBasis(columnComponentIdx)%ptr, &
            & numberOfColumnElementParameters(columnComponentIdx),err,error,*999)
          NULLIFY(columnQuadratureScheme(columnComponentIdx)%ptr)
          CALL Basis_QuadratureSchemeGet(columnBasis(columnComponentIdx)%ptr,BASIS_DEFAULT_QUADRATURE_SCHEME, &
            & columnQuadratureScheme(columnComponentIdx)%ptr,err,error,*999)
        ENDDO !columnComponentIdx
      ENDIF
      
      NULLIFY(materialsField)
      NULLIFY(materialsVariable)
      NULLIFY(materialsInterpParameters)
      NULLIFY(materialsInterpPoint)
      IF(smoothingType==EQUATIONS_SET_FITTING_SOBOLEV_VALUE_SMOOTHING.OR. &
        & smoothingType==EQUATIONS_SET_FITTING_SOBOLEV_DIFFERENCE_SMOOTHING) THEN
        CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)
        CALL Field_VariableGet(materialsField,FIELD_U_VARIABLE_TYPE,materialsVariable,err,error,*999)
        CALL EquationsInterpolation_MaterialsParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
          & materialsInterpParameters,err,error,*999)
        CALL EquationsInterpolation_MaterialsPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
          & materialsInterpPoint,err,error,*999)
        CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,materialsInterpParameters, &
          & err,error,*999)
      ENDIF

      NULLIFY(independentField)
      NULLIFY(dataVariable)
      NULLIFY(dataWeightVariable)
      NULLIFY(independentVectorParameters)
      NULLIFY(independentWeightParameters)
      NULLIFY(dataInterpParameters)
      NULLIFY(dataInterpPoint)
      NULLIFY(dataWeightInterpParameters)
      NULLIFY(dataWeightInterpPoint)
      
      SELECT CASE(esSpecification(2))
      CASE(EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE)
        
        !============================
        ! D a t a   P o i n t   F i t
        !============================
        
        CALL EquationsSet_IndependentFieldGet(equationsSet,independentField,err,error,*999)
        CALL Field_VariableGet(independentField,FIELD_U_VARIABLE_TYPE,dataVariable,err,error,*999)
        CALL FieldVariable_NumberOfComponentsGet(dataVariable,numberOfDataComponents,err,error,*999)
        CALL Field_VariableGet(independentField,FIELD_V_VARIABLE_TYPE,dataWeightVariable,err,error,*999)
        IF(numberOfDataComponents>MAXIMUM_DATA_COMPONENTS) THEN
          localError="Increase the size of the data point vectors. The data variable has "// &
            & TRIM(NumberToVString(numberOfDataComponents,"*",err,error))//" and the data point vectors only have "// &
            & TRIM(NumberToVString(MAXIMUM_DATA_COMPONENTS,"*",err,error))//" components."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        !Get data point vector parameters
        CALL FieldVariable_ParameterSetDataGet(dataVariable,FIELD_VALUES_SET_TYPE,independentVectorParameters, &
          & err,error,*999)
        !Get data point weight parameters
        CALL FieldVariable_ParameterSetDataGet(dataWeightVariable,FIELD_VALUES_SET_TYPE,independentWeightParameters, &
          & err,error,*999)

        NULLIFY(dataProjection)
        CALL Field_DataProjectionGet(independentField,dataProjection,err,error,*999)
        NULLIFY(independentDecomposition)
        CALL Field_DecompositionGet(independentField,independentDecomposition,err,error,*999)
        NULLIFY(independentDecompositionTopology)
        CALL Decomposition_DecompositionTopologyGet(independentDecomposition,independentDecompositionTopology,err,error,*999)
        NULLIFY(decompositionDataPoints)
        CALL DecompositionTopology_DecompositionDataPointsGet(independentDecompositionTopology,decompositionDataPoints, &
          & err,error,*999)
        CALL DecompositionDataPoints_ElementNumberOfDataPointsGet(decompositionDataPoints,elementNumber, &
          & numberOfElementDataPoints,err,error,*999)
        
        SELECT CASE(esSpecification(3))
        CASE(EQUATIONS_SET_DIFFUSION_TENSOR_FIBRE_DATA_FITTING_SUBTYPE)
          
          !=======================================================================
          ! D i f f u s i o n  T e n s o r  F i b r e  D a t a   P o i n t   F i t
          !=======================================================================
          
          !Loop over data points
          DO dataPointIdx=1,numberOfElementDataPoints
            CALL DecompositionDataPoints_ElementDataNumbersGet(decompositionDataPoints,dataPointIdx,elementNumber, &
              & dataPointLocalNumber,dataPointGlobalNumber,dataPointUserNumber,err,error,*999)
            !Need to use global number to get the correct projection results
            CALL DataProjection_ResultElementXiGet(dataProjection,dataPointGlobalNumber,numberOfElementXi,elementXi,err,error,*999)
            !Get data point vector value and weight
            DO componentIdx=1,numberOfDataComponents
              CALL FieldVariable_LocalDataPointDOFGet(dataVariable,dataPointLocalNumber,componentIdx,localDOFIdx, &
                & err,error,*999)
              dataPointVector(componentIdx)=independentVectorParameters(localDOFIdx)
              CALL FieldVariable_LocalDataPointDOFGet(dataWeightVariable,dataPointLocalNumber,componentIdx,localDOFIdx, &
                & err,error,*999)
              dataPointWeight(componentIdx)=independentWeightParameters(localDOFIdx)
            ENDDO !componentIdx
            CALL Field_InterpolateXi(FIRST_PART_DERIV,elementXi,geometricInterpPoint,err,error,*999)
            CALL Field_InterpolateXi(FIRST_PART_DERIV,elementXi,dependentInterpPoint,err,error,*999)
            CALL Field_InterpolatedPointMetricsCalculate(numberOfXi,geometricInterpPointMetrics,err,error,*999)
            
            rowElementDOFIdx=0
            !Loop over element rows
            DO rowComponentIdx=1,numberOfRowsComponents
              DO rowElementParameterIdx=1,numberOfRowElementParameters(rowComponentIdx)
                rowElementDOFIdx=rowElementDOFIdx+1
                CALL Basis_EvaluateXi(rowBasis(rowComponentIdx)%ptr,rowElementParameterIdx,NO_PART_DERIV, &
                  & elementXi(1:numberOfXi),rowPhi,err,error,*999)
                IF(updateResidual) THEN
                  columnElementDOFIdx=0
                  !Loop over element columns
                  DO columnComponentIdx=1,numberOfColsComponents
                    sum = 0.0_DP
                    !Treat each component as separate and independent so only calculate the diagonal blocks
                    IF(columnComponentIdx==rowComponentIdx) THEN
                      DO columnElementParameterIdx=1,numberOfColumnElementParameters(columnComponentIdx)
                        columnElementDOFIdx=columnElementDOFIdx+1
                        CALL Basis_EvaluateXi(columnBasis(columnComponentIdx)%ptr,columnElementParameterIdx,NO_PART_DERIV, &
                          & elementXi(1:numberOfXi),columnPhi,err,error,*999)
                        !sum = rowPhi*columnPhi*dataPointWeight(columnComponentIdx)* &
                        !  & dependentInterpPoint%values(rowComponentIdx,NO_PART_DERIV)
                      ENDDO !columnElementParameterIdx
                    ENDIF
                  ENDDO !columnComponentIdx
                  residualVector%elementResidual%vector(rowElementParameterIdx) = &
                    & residualVector%elementResidual%vector(rowElementParameterIdx) + sum
                ENDIF !update residual
                IF(updateRHS) THEN
                  sum = rowPhi*dataPointVector(rowComponentIdx)*dataPointWeight(rowComponentIdx)
                  rhsVector%elementVector%vector(rowElementParameterIdx)= &
                    & rhsVector%elementVector%vector(rowElementParameterIdx)+sum
                ENDIF !update RHS
                ENDDO !dependentElementParameterRowIdx
              ENDDO !dependentComponentRowIdx
            ENDDO !dataPointIdx

            !Restore data point vector parameters
            CALL Field_ParameterSetDataRestore(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & independentVectorParameters,err,error,*999)
            !Restore data point weight parameters
            CALL Field_ParameterSetDataRestore(independentField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & independentWeightParameters,err,error,*999)

            SELECT CASE(smoothingType)
            CASE(EQUATIONS_SET_FITTING_NO_SMOOTHING)
              !Do nothing
            CASE(EQUATIONS_SET_FITTING_SOBOLEV_VALUE_SMOOTHING)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(EQUATIONS_SET_FITTING_SOBOLEV_DIFFERENCE_SMOOTHING)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(EQUATIONS_SET_FITTING_STRAIN_ENERGY_SMOOTHING)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE DEFAULT
              localError="The fitting smoothing type of "//TRIM(NumberToVString(smoothingType,"*",err,error))// &
                & " is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT

          CASE DEFAULT
            localError="Equations set subtype "//TRIM(NumberToVString(equationsSet%specification(3),"*",err,error))// &
              & " is not valid for a data fitting equations set class."
            CALL FlagError(localError,err,error,*999)
          END SELECT

        CASE DEFAULT
          localError="Equations set type "//TRIM(NumberToVString(equationsSet%specification(2),"*",err,error))// &
            & " is not valid for a fitting equations set class."
          CALL FlagError(localError,err,error,*999)
        END SELECT

      !Scale factor adjustment
      CALL Field_ScalingTypeGet(dependentField,scalingType,err,error,*999)
      IF(scalingType/=FIELD_NO_SCALING) THEN
        NULLIFY(rowsInterpParameters)
        CALL EquationsInterpolation_DependentParametersGet(equationsInterpolation,rowsVariableType,rowsInterpParameters, &
          & err,error,*999)
        CALL Field_InterpolationParametersScaleFactorsElementGet(elementNumber,rowsInterpParameters,err,error,*999)
        !Loop over element rows
        rowElementDOFIdx=0          
        DO rowComponentIdx=1,numberOfRowsComponents
          DO rowElementParameterIdx=1,numberOfRowElementParameters(rowComponentIdx)
            rowElementDOFIdx=rowElementDOFIdx+1                    
            IF(updateResidual) THEN
              residualVector%elementResidual%vector(rowElementDOFIdx)= &
                & residualVector%elementResidual%vector(rowElementDOFIdx)* &
                & rowsInterpParameters%scaleFactors(rowElementParameterIdx,rowComponentIdx)
            ENDIF !update RHS
            IF(updateRHS) THEN
              rhsVector%elementVector%vector(rowElementDOFIdx)= &
                & rhsVector%elementVector%vector(rowElementDOFIdx)* &
                & rowsInterpParameters%scaleFactors(rowElementParameterIdx,rowComponentIdx)
            ENDIF !update RHS
         ENDDO !rowElementParameterIdx
        ENDDO !rowComponentIdx
      ENDIF !scale factors
      
    ENDIF !update

    EXITS("Fitting_FiniteElementResidualEvaluate")
    RETURN
999 ERRORS("Fitting_FiniteElementResidualEvaluate",err,error)
    EXITS("Fitting_FiniteElementResidualEvaluate")
    RETURN 1

  END SUBROUTINE Fitting_FiniteElementResidualEvaluate

  !
  !================================================================================================================================
  !

  SUBROUTINE Fitting_ProblemSetup(problem,problemSetup,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to setup
    TYPE(ProblemSetupType), INTENT(INOUT) :: problemSetup !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: pSpecification(3)
    TYPE(ControlLoopType), POINTER :: controlLoop, controlLoopRoot
    TYPE(SolverType), POINTER :: solver
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolversType), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError

    ENTERS("Fitting_ProblemSetup",err,error,*999)

    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)

    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STATIC_LINEAR_FITTING_SUBTYPE, &
      & PROBLEM_STATIC_NONLINEAR_FITTING_SUBTYPE, &
      & PROBLEM_QUASISTATIC_LINEAR_FITTING_SUBTYPE, &
      & PROBLEM_QUASISTATIC_NONLINEAR_FITTING_SUBTYPE, &
      & PROBLEM_DIV_FREE_VELOCITY_FITTING_SUBTYPE)
      !OK
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a fitting problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
      
    SELECT CASE(problemSetup%setupType)
    CASE(PROBLEM_SETUP_INITIAL_TYPE)
      !
      ! Initial Setup
      !
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)
        !Do nothing????
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Do nothing???
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a static fitting problem."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_SETUP_CONTROL_TYPE)
      !
      ! Control loop setup
      !
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)
        !Set up a simple control loop
        NULLIFY(controlLoop)
        CALL ControlLoop_CreateStart(problem,controlLoop,err,error,*999)
        IF(pSpecification(3)==PROBLEM_STATIC_LINEAR_FITTING_SUBTYPE.OR. &
          & pSpecification(3)==PROBLEM_STATIC_NONLINEAR_FITTING_SUBTYPE) THEN
          CALL ControlLoop_TypeSet(controlLoop,CONTROL_SIMPLE_TYPE,err,error,*999)
        ELSE
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
          & " is invalid for a static fitting problem."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_SETUP_SOLVERS_TYPE)
      !
      ! Solvers setup
      !
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
        IF(pSpecification(3)==PROBLEM_STATIC_NONLINEAR_FITTING_SUBTYPE.OR. &
          & pSpecification(3)==PROBLEM_QUASISTATIC_NONLINEAR_FITTING_SUBTYPE) THEN
          !Set the solver to be a nonlinear solver
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          CALL Solver_TypeSet(solver,SOLVER_NONLINEAR_TYPE,err,error,*999)
          !Set solver defaults
          CALL Solver_LibraryTypeSet(solver,SOLVER_PETSC_LIBRARY,err,error,*999)
        ELSE
          !Set the solver to be a linear solver
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          CALL Solver_TypeSet(solver,SOLVER_LINEAR_TYPE,err,error,*999)
          !Set solver defaults
          CALL Solver_LibraryTypeSet(solver,SOLVER_PETSC_LIBRARY,err,error,*999)
        ENDIF
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Get the solvers
        NULLIFY(solvers)
        CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
        !Finish the solvers creation
        CALL Solvers_CreateFinish(solvers,err,error,*999)
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a static fitting problem."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
      !
      ! Solver equations setup
      !
      !Get the control loop
      NULLIFY(controlLoopRoot)
      CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
      NULLIFY(controlLoop)
      CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
      !Get the solver
      NULLIFY(solvers)
      CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
      NULLIFY(solver)
      CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)
        !Create the solver equations
        NULLIFY(solverEquations)
        CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
        IF(pSpecification(3)==PROBLEM_STATIC_NONLINEAR_FITTING_SUBTYPE.OR. &
          & pSpecification(3)==PROBLEM_QUASISTATIC_NONLINEAR_FITTING_SUBTYPE ) THEN
          CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
        ELSE
          CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_LINEAR,err,error,*999)
        ENDIF
        IF(pSpecification(3)==PROBLEM_STATIC_LINEAR_FITTING_SUBTYPE.OR. &
          & pSpecification(3)==PROBLEM_STATIC_NONLINEAR_FITTING_SUBTYPE ) THEN
          CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_STATIC,err,error,*999)
        ELSE
          CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_QUASISTATIC,err,error,*999)
        ENDIF
        CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        NULLIFY(solverEquations)
        CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
        !Finish the solver equations creation
        CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a static fitting problem."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
        & " is invalid for a static fitting problem."
      CALL FlagError(localError,err,error,*999)
    END SELECT
 
    EXITS("Fitting_ProblemSetup")
    RETURN
999 ERRORSEXITS("Fitting_ProblemSetup",err,error)
    RETURN 1

  END SUBROUTINE Fitting_ProblemSetup

  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a data fitting problem class.
  SUBROUTINE Fitting_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to set the specification for
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The proboem specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: problemType,problemSubtype
    TYPE(VARYING_STRING) :: localError

    ENTERS("Fitting_ProblemSpecificationSet",err,error,*999)

    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Problem is not associated.",err,error,*999)
    IF(ALLOCATED(problem%specification)) CALL FlagError("Problem specification is already allocated.",err,error,*999)
    IF(SIZE(problemSpecification,1)<3) THEN
      localError="The size of the specified problem specification array of "// &
        & TRIM(NumberToVString(SIZE(problemSpecification,1),"*",err,error))// &
        & " is invalid. The size should be >= 3."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    problemType=problemSpecification(2)
    problemSubtype=problemSpecification(3)
    SELECT CASE(problemType)
    CASE(PROBLEM_FITTING_TYPE)
      SELECT CASE(problemSubtype)
      CASE(PROBLEM_STATIC_LINEAR_FITTING_SUBTYPE, &
        & PROBLEM_STATIC_NONLINEAR_FITTING_SUBTYPE, &
        & PROBLEM_QUASISTATIC_LINEAR_FITTING_SUBTYPE, &
        & PROBLEM_QUASISTATIC_NONLINEAR_FITTING_SUBTYPE, &
        & PROBLEM_DIV_FREE_VELOCITY_FITTING_SUBTYPE)
        !ok
      CASE DEFAULT
        localError="The third problem specification of "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
          & " is not valid for a Galerkin projection type of a data fitting problem."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The second problem specification of "//TRIM(NumberToVstring(problemType,"*",err,error))// &
        & " is not valid for a data fitting problem."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    ALLOCATE(problem%specification(3),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
    problem%specification(1:3)=[PROBLEM_FITTING_CLASS,problemType,problemSubtype]

    EXITS("Fitting_ProblemSpecificationSet")
    RETURN
999 ERRORS("Fitting_ProblemSpecificationSet",err,error)
    EXITS("Fitting_ProblemSpecificationSet")
    RETURN 1

  END SUBROUTINE Fitting_ProblemSpecificationSet
  
  !
  !================================================================================================================================
  !

 !>Sets up the output type for a data fitting problem class.
  SUBROUTINE Fitting_PreSolve(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: pSpecification(3)
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolversType), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError

    ENTERS("Fitting_PreSolve",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)

    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STATIC_LINEAR_FITTING_SUBTYPE)
      !Do nothing
    CASE(PROBLEM_STATIC_NONLINEAR_FITTING_SUBTYPE)
      !Do nothing
    CASE(PROBLEM_QUASISTATIC_LINEAR_FITTING_SUBTYPE)
      !Do nothing
    CASE(PROBLEM_QUASISTATIC_NONLINEAR_FITTING_SUBTYPE)
      !Do nothing
    CASE(PROBLEM_DIV_FREE_VELOCITY_FITTING_SUBTYPE)
      !Update indpendent data fields
      CALL Fitting_PreSolveUpdateInputData(solver,err,error,*999)
    CASE DEFAULT
      localError="The third problem specification of "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a data fitting problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Fitting_PreSolve")
    RETURN
999 ERRORSEXITS("Fitting_PreSolve",err,error)
    RETURN 1

  END SUBROUTINE Fitting_PreSolve

  !
  !================================================================================================================================
  !

  !>Sets up the output type for a data fitting problem class.
  SUBROUTINE Fitting_PostSolve(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: pSpecification(3)
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemType), POINTER :: problem
    TYPE(VARYING_STRING) :: localError

    ENTERS("Fitting_PostSolve",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)

    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STATIC_LINEAR_FITTING_SUBTYPE)
      !Do nothing
    CASE(PROBLEM_STATIC_NONLINEAR_FITTING_SUBTYPE)
      !Do nothing
    CASE(PROBLEM_QUASISTATIC_LINEAR_FITTING_SUBTYPE)
      !Do nothing
    CASE(PROBLEM_QUASISTATIC_NONLINEAR_FITTING_SUBTYPE)
      !Do nothing
    CASE(PROBLEM_DIV_FREE_VELOCITY_FITTING_SUBTYPE)
      CALL Fitting_PostSolveOutputData(solver,err,error,*999)
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a fitting problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Fitting_PostSolve")
    RETURN
999 ERRORSEXITS("Fitting_PostSolve",err,error)
    RETURN 1

  END SUBROUTINE Fitting_PostSolve

  !
  !================================================================================================================================
  !

  !>Output data post solve
  SUBROUTINE Fitting_PostSolveOutputData(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: currentIteration,equationsSetIdx,globalNumber,inputIteration,numberOfEquationsSets,outputIteration, &
      & outputType,pSpecification(3)
    REAL(DP) :: currentTime,startTime,stopTime,timeIncrement
    LOGICAL :: exportField
    CHARACTER(14) :: outputFile
    TYPE(ControlLoopType), POINTER :: controlLoop 
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(ProblemType), POINTER :: problem 
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(RegionType), POINTER :: region
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Fitting_PostSolveOutputData",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)

    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STATIC_LINEAR_FITTING_SUBTYPE)
      !Do nothing
    CASE(PROBLEM_STATIC_NONLINEAR_FITTING_SUBTYPE)
      !Do nothing
    CASE(PROBLEM_QUASISTATIC_LINEAR_FITTING_SUBTYPE)
      !Do nothing
    CASE(PROBLEM_QUASISTATIC_NONLINEAR_FITTING_SUBTYPE)
      !Do nothing
    CASE(PROBLEM_DIV_FREE_VELOCITY_FITTING_SUBTYPE)
      CALL Solver_OutputTypeGet(solver,outputType,err,error,*999)
      CALL ControlLoop_CurrentTimeInformationGet(controlLoop,currentTime,timeIncrement,startTime,stopTime,currentIteration, &
        & outputIteration,inputIteration,err,error,*999)
      NULLIFY(solverEquations)
      CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
      DO equationsSetIdx=1,numberOfEquationsSets
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
        IF(outputIteration/=0) THEN
          IF(currentTime<stopTime) THEN
            IF(currentIteration<10) THEN
              WRITE(outputFile,'("FittedData_0",I0)') currentIteration
            ELSE IF(currentIteration<100) THEN
              WRITE(outputFile,'("FittedData_",I0)') currentIteration
            ENDIF
            exportField=.TRUE.
            IF(exportField) THEN
              IF(MOD(currentIteration,outputIteration)==0)  THEN
                IF(outputType>=SOLVER_PROGRESS_OUTPUT) THEN
                  CALL WriteString(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
                  CALL WriteString(GENERAL_OUTPUT_TYPE,"Now export fields... ",err,error,*999)
                ENDIF
                NULLIFY(region)
                CALL EquationsSet_RegionGet(equationsSet,region,err,error,*999)
                CALL EquationsSet_GlobalNumberGet(equationsSet,globalNumber,err,error,*999)
                CALL FLUID_MECHANICS_IO_WRITE_FITTED_FIELD(region,globalNumber,outputFile,err,error,*999)
                IF(outputType>=SOLVER_PROGRESS_OUTPUT) THEN
                  CALL WriteString(GENERAL_OUTPUT_TYPE,outputFile,err,error,*999)
                  CALL WriteString(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDDO !equationsSetIdx
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(problem%specification(3),"*",err,error))// &
        & " is not valid for a fitting equation of a classical field problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Fitting_PostSolveOutputData")
    RETURN
999 ERRORSEXITS("Fitting_PostSolveOutputData",err,error)
    RETURN 1

  END SUBROUTINE Fitting_PostSolveOutputData

  !
  !================================================================================================================================
  !

  !>Update input data conditions for field fitting
  SUBROUTINE Fitting_PreSolveUpdateInputData(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: currentIteration,inputIteration,inputOption,inputType,numberOfDimensions,outputIteration,outputType, &
      & pSpecification(3)
    REAL(DP) :: currentTime,startTime,stopTime,timeIncrement
    REAL(DP), POINTER :: inputVelNewData(:)
    TYPE(ControlLoopType), POINTER :: controlLoop,controlTimeLoop
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(FieldType), POINTER :: geometricField,sourceField
    TYPE(FieldVariableType), POINTER :: sourceVariable
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolversType), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError

    ENTERS("Fitting_PreSolveUpdateInputData",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)

    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STATIC_LINEAR_FITTING_SUBTYPE)
      !Do nothing
    CASE(PROBLEM_STATIC_NONLINEAR_FITTING_SUBTYPE)
      !Do nothing
    CASE(PROBLEM_QUASISTATIC_LINEAR_FITTING_SUBTYPE)
      !Do nothing
    CASE(PROBLEM_QUASISTATIC_NONLINEAR_FITTING_SUBTYPE)
      !Do nothing
    CASE(PROBLEM_DIV_FREE_VELOCITY_FITTING_SUBTYPE)
      CALL Solver_OutputTypeGet(solver,outputType,err,error,*999)
      CALL ControlLoop_CurrentTimeInformationGet(controlLoop,currentTime,timeIncrement,startTime,stopTime,currentIteration, &
        & outputIteration,inputIteration,err,error,*999)
      IF(outputType>=SOLVER_PROGRESS_OUTPUT) CALL WriteString(GENERAL_OUTPUT_TYPE,"Read input data... ",err,error,*999)
      NULLIFY(solverEquations)
      CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
      NULLIFY(sourceField)
      CALL EquationsSet_SourceFieldGet(equationsSet,sourceField,err,error,*999)
      NULLIFY(sourceVariable)
      CALL Field_VariableGet(sourceField,FIELD_U_VARIABLE_TYPE,sourceVariable,err,error,*999)
      !this is the current time step
      !\todo: Provide possibility for user to define input type and option (that's more or less an IO question)
      inputType=1
      inputOption=1
      NULLIFY(inputVelNewData)
      CALL FieldVariable_ParameterSetDataGet(sourceVariable,FIELD_VALUES_SET_TYPE,inputVelNewData,err,error,*999)
      CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,inputVelNewData,numberOfDimensions,inputType,inputOption, &
        & currentIteration,1.0_DP,err,error,*999)
      CALL FieldVariable_ParameterSetUpdateStart(sourceVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
      CALL FieldVariable_ParameterSetUpdateFinish(sourceVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a fitting problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Fitting_PreSolveUpdateInputData")
    RETURN
999 ERRORSEXITS("Fitting_PreSolveUpdateInputData",err,error)
    RETURN 1

  END SUBROUTINE Fitting_PreSolveUpdateInputData

  !
  !================================================================================================================================
  !

END MODULE FittingRoutines
