!> \file
!> \author Vijay Rajagopal
!> \brief This module handles all reaction diffusion equation routines.
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
!> Contributor(s): Vijay Rajagopal, Chris Bradley
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

!>This module handles all reaction diffusion equation routines.
MODULE ReactionDiffusionEquationsRoutines

  USE BaseRoutines
  USE BasisRoutines
  USE BasisAccessRoutines
  USE BoundaryConditionsRoutines
  USE ComputationRoutines
  USE ComputationAccessRoutines
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
  USE FieldAccessRoutines
  USE InputOutput
  USE ISO_VARYING_STRING
  USE Kinds
  USE MatrixVector
#ifndef NOMPIMOD
  USE MPI
#endif
  USE ProblemAccessRoutines
  USE Strings
  USE SolverRoutines
  USE SolverAccessRoutines
  USE Timer
  USE Types
  
  USE REACTION_DIFFUSION_IO_ROUTINES

#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

#ifdef NOMPIMOD
#include "mpif.h"
#endif

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC ReactionDiffusion_EquationsSetSetup

  PUBLIC ReactionDiffusion_EquationsSetSolutionMethodSet
  
  PUBLIC ReactionDiffusion_EquationsSetSpecificationSet

  PUBLIC ReactionDiffusion_FiniteElementCalculate
  
  PUBLIC ReactionDiffusion_PreSolve
  
  PUBLIC ReactionDiffusion_ProblemSetup

  PUBLIC ReactionDiffusion_ProblemSpecificationSet

  PUBLIC ReactionDiffusion_PostSolve

  PUBLIC ReactionDiffusion_PostLoop


CONTAINS

  !
  !================================================================================================================================
  !

  !>Sets up the reaction diffusion equation type of a classical equations set class.
  SUBROUTINE ReactionDiffusion_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to setup a bioelectric domain equation on.
    TYPE(EquationsSetSetupType), INTENT(INOUT) :: equationsSetSetup !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx,dimensionMultiplier,geometricComponentNumber,geometricScalingType, &
      & numberOfDimensions,numberOfMaterialsComponents,geometricMeshComponent
    TYPE(DecompositionType), POINTER :: geometricDecomposition
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsSetMaterialsType), POINTER :: equationsMaterials
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(FieldType), POINTER :: geometricField
    TYPE(RegionType), POINTER :: regionType
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("ReactionDiffusion_EquationsSetSetup",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
    
    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_CONSTANT_REAC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_CELLML_REAC_NO_SPLIT_REAC_DIFF_SUBTYPE)
      !OK
    CASE DEFAULT
      localError="The third equations set specification of "//TRIM(NumberToVstring(esSpecification(3),"*",err,error))// &
        & " is not valid for a reaction-diffusion type of a classical field equations set."
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
      CASE(equationsSet_SETUP_START_ACTION)
        CALL ReactionDiffusion_EquationsSetSolutionMethodSet(equationsSet,EQUATIONS_SET_FEM_SOLUTION_METHOD,err,error,*999)
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
!!Todo: CHECK VALID SETUP
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a reaction diffusion domain equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
      !-----------------------------------------------------------------
      ! G e o m e t r i c   f i e l d
      !-----------------------------------------------------------------
      !\todo Check geometric dimension
    CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
      !-----------------------------------------------------------------
      ! D e p e n d e n t   f i e l d
      !-----------------------------------------------------------------
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        SELECT CASE(esSpecification(3))
        CASE(EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE, &
          & EQUATIONS_SET_CONSTANT_REAC_DIFF_SUBTYPE)
          IF(equationsSet%dependent%dependentFieldAutoCreated) THEN
              !Create the auto created dependent field
            CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsSet%dependent%dependentField,err,error,*999)
            CALL Field_LabelSet(equationsSet%dependent%dependentField,"Dependent Field",err,error,*999)
            CALL Field_TypeSetAndLock(equationsSet%dependent%dependentField,FIELD_GENERAL_TYPE,err,error,*999)
            CALL Field_DependentTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_DEPENDENT_TYPE,err,error,*999)
            CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
            CALL Field_DecompositionSetAndLock(equationsSet%dependent%dependentField,geometricDecomposition,err,error,*999)
            CALL Field_GeometricFieldSetAndLock(equationsSet%dependent%dependentField,geometricField,err,error,*999)
            CALL Field_NumberOfVariablesSetAndLock(equationsSet%dependent%dependentField,2,err,error,*999)
            CALL Field_VariableTypesSetAndLock(equationsSet%dependent%dependentField,[FIELD_U_VARIABLE_TYPE, &
              & FIELD_DELUDELN_VARIABLE_TYPE],err,error,*999)
            !CALL Field_VariableLabelSet
            CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
              & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE, &
              & err,error,*999)
            CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,1, &
              & err,error,*999)
            CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,1, &
              & err,error,*999)
            !Default to the geometric interpolation setup
            CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricMeshComponent,err,error,*999)
            CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,1, &
              & geometricMeshComponent,err,error,*999)
            CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,1, &
              & geometricMeshComponent,err,error,*999)
            CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
            SELECT CASE(solutionMethod)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,1, &
                & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
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
              localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ELSE
            !Check the user specified field
            CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
            CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_DEPENDENT_TYPE,err,error,*999)
            CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,2,err,error,*999)
            CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE], &
              & err,error,*999)
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
              & err,error,*999)
            CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,1,err,error,*999)
            CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,1,err,error,*999)
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
              localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ENDIF
        CASE(EQUATIONS_SET_CELLML_REAC_NO_SPLIT_REAC_DIFF_SUBTYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE DEFAULT
          localError="The equation set subtype of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
            & " is invalid for a reaction diffusion equation set class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        IF(equationsSet%dependent%dependentFieldAutoCreated) THEN
          CALL Field_CreateFinish(equationsSet%dependent%dependentField,err,error,*999)
        ENDIF
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a reaction diffusion equation"
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
      !-----------------------------------------------------------------
      ! M a t e r i a l s   f i e l d
      !-----------------------------------------------------------------
      NULLIFY(equationsMaterials)
      CALL EquationsSet_MaterialsGet(equationsSet,equationsMaterials,err,error,*999)
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
      IF(esSpecification(3)==EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE) THEN
        !a.\delby{u}{t}+div(\sigma.\grad u)+b.cellmlRC=0
        !Reaction Diffusion. Materials field components are a + diffusion coeff tensor components for each dimension
        numberOfMaterialsComponents=1+NUMBER_OF_VOIGT(numberOfDimensions)
        dimensionMultiplier=1
      ELSE IF(esSpecification(3)==EQUATIONS_SET_CONSTANT_REAC_DIFF_SUBTYPE) THEN
        !Constant reaction + diffusion. Materials field has 1 diffuse coeff for each dimension
        !plus one for the storage coefficient om  alpha(delC/delt) = Div(-kgradC)+const(x)_source
        numberOfMaterialsComponents=1+NUMBER_OF_VOIGT(numberOfDimensions)
        dimensionMultiplier=1
      ENDIF
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        CALL EquationsSet_AssertDependentIsFinished(equationsSet,err,error,*999)
        IF(equationsMaterials%materialsFieldAutoCreated) THEN
          !Create the auto created materials field
          CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsMaterials%materialsField,err,error,*999)
          CALL Field_TypeSetAndLock(equationsMaterials%materialsField,FIELD_MATERIAL_TYPE,err,error,*999)
          CALL Field_DependentTypeSetAndLock(equationsMaterials%materialsField,FIELD_INDEPENDENT_TYPE,err,error,*999)
          CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
          CALL Field_DecompositionSetAndLock(equationsMaterials%materialsField,geometricDecomposition,err,error,*999)
          CALL Field_GeometricFieldSetAndLock(equationsMaterials%materialsField,geometricField,err,error,*999)
          CALL Field_NumberOfVariablesSetAndLock(equationsMaterials%materialsField,1,err,error,*999)
          CALL Field_VariableTypesSetAndLock(equationsMaterials%materialsField,[FIELD_U_VARIABLE_TYPE], err,error,*999)
          CALL Field_DimensionSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
            & err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          !Set the number of materials components
          CALL Field_NumberOfComponentsSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
            & numberOfMaterialsComponents,err,error,*999)
          !Default the materials components for diffusivity param to the first component geometric interpolation
          !with const interpolation
          CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricComponentNumber,err,error,*999)
          DO componentIdx=1,numberOfDimensions
            CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
              & componentIdx,geometricComponentNumber,err,error,*999)
            CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
              & componentIdx,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
          ENDDO !components_idx
          !Default the storage co-efficient to the first geometric interpolation setup with constant interpolation
          IF(esSpecification(3)==EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE) THEN
            componentIdx=numberOfMaterialsComponents
            CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricComponentNumber,err,error,*999)
            CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,componentIdx, &
              & geometricComponentNumber,err,error,*999)
            CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,componentIdx, &
              & FIELD_CONSTANT_INTERPOLATION,err,error,*999)
          ENDIF
          !Default the field scaling to that of the geometric field
          CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
          CALL Field_ScalingTypeSet(equationsMaterials%materialsField,geometricScalingType,err,error,*999)
        ELSE
          !Check the user specified field
          CALL Field_TypeCheck(equationsSetSetup%field,FIELD_MATERIAL_TYPE,err,error,*999)
          CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
          CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,1,err,error,*999)
          CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE],err,error,*999)
          CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
            & err,error,*999)
          CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          IF(esSpecification(3)==EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE .OR. &
            & esSpecification(3)==EQUATIONS_SET_CONSTANT_REAC_DIFF_SUBTYPE) THEN
            !Reaction Diffusion with cellml. Materials field components are 1 for storage coeff plus one for each dimension i.e., k
            CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfDimensions+1, &
              & err,error,*999)
          ENDIF
        ENDIF
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        IF(equationsMaterials%materialsFieldAutoCreated) THEN
          !Finish creating the materials field
          CALL Field_CreateFinish(equationsMaterials%materialsField,err,error,*999)
          !Set the default values for the materials field
          IF(esSpecification(3)==EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE .OR. &
            & esSpecification(3)==EQUATIONS_SET_CONSTANT_REAC_DIFF_SUBTYPE) THEN
            !Reaction Diffusion with cellml. Materials field components are 1 plus one for each dimension i.e.,storage coeff, and k.
            numberOfMaterialsComponents=numberOfDimensions+1
            dimensionMultiplier=1
          ENDIF
          !set the diffusion coefficients to be 1.0
          DO componentIdx=1,numberOfDimensions
            CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_VALUES_SET_TYPE,componentIdx,1.0_DP,err,error,*999)
          ENDDO !componentIdx
          !Now set storage-coefficient
          componentIdx=numberOfDimensions+1
          CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,componentIdx,1.0_DP,err,error,*999)
        ENDIF
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a reaction diffusion equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
      NULLIFY(equationsSource)
      CALL EquationsSet_SourceGet(equationsSet,equationsSource,err,error,*999)
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        CALL EquationsSet_AssertMaterialsIsFinished(equationsSet,err,error,*999)
        IF(equationsSource%sourceFieldAutoCreated) THEN
          !Create the auto created source field
          !Start field creation with name 'Source Field'
          CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsSource%sourceField,err,error,*999)
          !Create a general field
          CALL Field_TypeSetAndLock(equationsSource%sourceField,FIELD_GENERAL_TYPE,err,error,*999)
          !Label the field
          CALL Field_LabelSetAndLock(equationsSource%sourceField,"Source Field",err,error,*999)
          !Set the dependent type
          CALL Field_DependentTypeSetAndLock(equationsSource%sourceField,FIELD_INDEPENDENT_TYPE,err,error,*999)
          !Set the field decomposition to be that of the geometric decomposition
          CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
          CALL Field_DecompositionSetAndLock(equationsSource%sourceField,geometricDecomposition,err,error,*999)
          !Set the geometric field
          CALL Field_GeometricFieldSetAndLock(equationsSource%sourceField,geometricField,err,error,*999)
          !Set the field variables.
          CALL Field_NumberOfVariablesSetAndLock(equationsSource%sourceField,1,err,error,*999)
          CALL Field_VariableTypesSetAndLock(equationsSource%sourceField,[FIELD_U_VARIABLE_TYPE],err,error,*999)
          !Set the dimension
          CALL Field_DimensionSetAndLock(equationsSource%sourceField,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
            & err,error,*999)
          !Set the data type
          CALL Field_DataTypeSetAndLock(equationsSource%sourceField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          !Set the number of components to one
          CALL Field_NumberOfComponentsSetAndLock(equationsSource%sourceField,FIELD_U_VARIABLE_TYPE,1,err,error,*999)
          !Get the geometric mesh component
          CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricMeshComponent,err,error,*999)
          !Default to the geometric interpolation setup
          CALL Field_ComponentMeshComponentSet(equationsSource%sourceField,FIELD_U_VARIABLE_TYPE,1,geometricMeshComponent, &
            & err,error,*999)
          !Specify the interpolation to be same as geometric interpolation
          CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
          SELECT CASE(solutionMethod)
          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
            CALL Field_ComponentInterpolationSetAndLock(equationsSource%sourceField,FIELD_U_VARIABLE_TYPE,1, &
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
            localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))// " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
          !Set the scaling to be the same as the geometric field
          CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
          CALL Field_ScalingTypeSet(equationsSource%sourceField,geometricScalingType,err,error,*999)
        ELSE
          !Check the user specified field
          CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
          CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
          CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,1,err,error,*999)
          CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE],err,error,*999)
          CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,1,err,error,*999)
          CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
          SELECT CASE(solutionMethod)
          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
            CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,1, &
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
        ENDIF
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        IF(equationsSource%sourceFieldAutoCreated) THEN
          CALL Field_CreateFinish(equationsSource%sourceField,err,error,*999)
        ENDIF
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a reaction diffusion equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a reaction diffusion equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        CALL EquationsSet_AssertSourceIsFinished(equationsSet,err,error,*999)
        !Create the equations
        NULLIFY(equations)
        CALL Equations_CreateStart(equationsSet,equations,err,error,*999)
        SELECT CASE(esSpecification(3))
        CASE(EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE)
          CALL Equations_LinearityTypeSet(equations,EQUATIONS_LINEAR,err,error,*999)
        CASE(EQUATIONS_SET_CONSTANT_REAC_DIFF_SUBTYPE)
          CALL Equations_LinearityTypeSet(equations,EQUATIONS_LINEAR,err,error,*999)
        CASE(EQUATIONS_SET_CELLML_REAC_NO_SPLIT_REAC_DIFF_SUBTYPE)
          CALL Equations_LinearityTypeSet(equations,EQUATIONS_NONLINEAR,err,error,*999)
        CASE DEFAULT
          localError="The third equations set specification of "//TRIM(NumberToVstring(esSpecification(3),"*",err,error))// &
            & " is not valid for a reaction-diffusion type of a classical field equations set."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        CALL Equations_TimeDependenceTypeSet(equations,EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
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
          CALL EquationsMappingVector_DynamicMatricesSet(vectorMapping,.TRUE.,.TRUE.,err,error,*999)
          CALL EquationsMappingVector_DynamicVariableTypeSet(vectorMapping,FIELD_U_VARIABLE_TYPE,err,error,*999)
          CALL EquationsMappingVector_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE,err,error,*999)
          CALL EquationsMappingVector_NumberOfSourcesSet(vectorMapping,1,err,error,*999)
          CALL EquationsMappingVector_SourceVariableTypeSet(vectorMapping,1,FIELD_U_VARIABLE_TYPE,err,error,*999)
          CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
          !Create the equations matrices
          NULLIFY(vectorMatrices)
          CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
          !Set up matrix storage and structure
          CALL EquationsSet_LumpingTypeGet(equations,lumpingType,err,error,*999)
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
              CALL EquationsMatricesVector_DynamicStorageTypeSet(vectorMatrices,[DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE, &
                & DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE],err,error,*999)
              CALL EquationsMatricesVector_DynamicStructureTypeSet(vectorMatrices,[EQUATIONS_MATRIX_FEM_STRUCTURE, &
                & EQUATIONS_MATRIX_FEM_STRUCTURE],err,error,*999)
            CASE DEFAULT
              localError="The equations matrices sparsity type of "//TRIM(NumberToVString(sparsityType,"*",err,error))// &
                & " is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ENDIF
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
          & " is invalid for a bioelectric domain equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
        & " is invalid for reaction diffusion equation."
      CALL FlagError(localError,err,error,*999)
    END SELECT
  
    EXITS("ReactionDiffusion_EquationsSetSetup")
    RETURN
999 ERRORSEXITS("ReactionDiffusion_EquationsSetSetup",err,error)
    RETURN 1
    
  END SUBROUTINE ReactionDiffusion_EquationsSetSetup

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a reaction diffusion equation type of a classical equations set class.
  SUBROUTINE ReactionDiffusion_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: solutionMethod !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("ReactionDiffusion_EquationsSetSolutionMethodSet",err,error,*999)
    
    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
    
    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_CELLML_REAC_NO_SPLIT_REAC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_CONSTANT_REAC_DIFF_SUBTYPE)        
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
        & " is not valid for a reaction diffusion equation type of classical equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("ReactionDiffusion_EquationsSetSolutionMethodSet")
    RETURN
999 ERRORS("ReactionDiffusion_EquationsSetSolutionMethodSet",err,error)
    EXITS("ReactionDiffusion_EquationsSetSolutionMethodSet")
    RETURN 1
    
  END SUBROUTINE ReactionDiffusion_EquationsSetSolutionMethodSet

  !
  !================================================================================================================================ 
  !

  !>Sets the equation specification for a reaction diffusion equation type of a classical equations set class.
  SUBROUTINE ReactionDiffusion_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the specification for
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations set specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: subtype
    TYPE(VARYING_STRING) :: localError

    ENTERS("ReactionDiffusion_EquationsSetSpecificationSet",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(ALLOCATED(equationsSet%specification)) CALL FlagError("Equations set specification is already allocated.",err,error,*999)    
    IF(SIZE(specification,1)<3) THEN
      localError="The size of the specified specification array of "// &
        & TRIM(NumberToVString(SIZE(specification,1),"*",err,error))//" is invalid. The size should be >= 3."
      CALL FlagError(localError,err,error,*999)
    END IF
    
    subtype=specification(3)
    SELECT CASE(subtype)
    CASE(EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE)
      !ok
    CASE(EQUATIONS_SET_CELLML_REAC_NO_SPLIT_REAC_DIFF_SUBTYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_CONSTANT_REAC_DIFF_SUBTYPE)
      !ok
    CASE DEFAULT
      localError="The specified equations set subtype of "//TRIM(NumberToVstring(subtype,"*",err,error))// &
        & " is not valid for reaction diffusion equation type of a classical equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    !Set full specification
    ALLOCATE(equationsSet%specification(3),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate equations set specification.",err,error,*999)
    equationsSet%specification(1:3)=[EQUATIONS_SET_CLASSICAL_FIELD_CLASS,EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE,subtype]

    EXITS("ReactionDiffusion_EquationsSetSpecificationSet")
    RETURN
999 ERRORS("ReactionDiffusion_EquationsSetSpecificationSet",err,error)
    EXITS("ReactionDiffusion_EquationsSetSpecificationSet")
    RETURN 1
    
  END SUBROUTINE ReactionDiffusion_EquationsSetSpecificationSet

  !
  !================================================================================================================================
  !
  !>Calculates the element stiffness matrices and RHS for a reaction diffusion equation finite element equations set.
  SUBROUTINE ReactionDiffusion_FiniteElementCalculate(equationsSet,ELEMENT_NUMBER,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) FIELD_VAR_TYPE,mh,mhs,ms,ng,nh,nhs,ni,nj,ns,componentIdx
    LOGICAL :: USE_FIBRES
    REAL(DP) :: DIFFUSIVITY(3,3),DPHIDX(3,64),RWG,SUM,STORAGE_COEFFICIENT,C_PARAM
    TYPE(BasisType), POINTER :: DEPENDENT_BASIS,GEOMETRIC_BASIS,FIBRE_BASIS
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatrixType), POINTER :: dampingMatrix,stiffnessMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField,geometricField,fibreField,materialsField,sourceField
    TYPE(FieldVariableType), POINTER :: FIELD_VARIABLE,GEOMETRIC_VARIABLE
    TYPE(QuadratureSchemeType), POINTER :: QUADRATURE_SCHEME
    
    ENTERS("ReactionDiffusion_FiniteElementCalculate",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(.NOT.ALLOCATED(equationsSet%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(equationsSet%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a reaction-diffusion type equations set.", &
          & err,error,*999)
      END IF
      EQUATIONS=>equationsSet%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        NULLIFY(vectorEquations)
        CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
        dependentField=>equations%interpolation%dependentField
        geometricField=>equations%interpolation%geometricField
        materialsField=>equations%interpolation%materialsField
        IF(equationsSet%SPECIFICATION(3)==EQUATIONS_SET_CELLML_REAC_NO_SPLIT_REAC_DIFF_SUBTYPE .OR. &
            & equationsSet%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_REAC_DIFF_SUBTYPE) THEN
           sourceField=>equations%interpolation%sourceField
        ENDIF
        fibreField=>equations%interpolation%fibreField
        USE_FIBRES=ASSOCIATED(fibreField)
        vectorMapping=>vectorEquations%vectorMapping
        vectorMatrices=>vectorEquations%vectorMatrices
        DEPENDENT_BASIS=>dependentField%DECOMPOSITION%DOMAIN(dependentField%decomposition%meshComponentNumber)%ptr% &
          & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
        GEOMETRIC_BASIS=>geometricField%DECOMPOSITION%DOMAIN(geometricField%decomposition%meshComponentNumber)%ptr% &
          & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
        GEOMETRIC_VARIABLE=>geometricField%variableTypeMap(FIELD_U_VARIABLE_TYPE)%ptr
        IF(USE_FIBRES) FIBRE_BASIS=>fibreField%DECOMPOSITION%DOMAIN(geometricField%decomposition%meshComponentNumber)%ptr% &
          & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
        QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%quadratureSchemeMap(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
        CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
          & geometricInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
        CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
          & materialsInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
        IF(USE_FIBRES) CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations% &
          & interpolation%fibreInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
        IF(equationsSet%SPECIFICATION(3)==EQUATIONS_SET_CELLML_REAC_NO_SPLIT_REAC_DIFF_SUBTYPE .OR. &
            & equationsSet%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_REAC_DIFF_SUBTYPE) THEN
           CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
             & sourceInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
        ENDIF                  
        dynamicMatrices=>vectorMatrices%dynamicMatrices
        stiffnessMatrix=>dynamicMatrices%matrices(1)%ptr
        dampingMatrix=>dynamicMatrices%matrices(2)%ptr
        rhsVector=>vectorMatrices%rhsVector
        IF(equationsSet%SPECIFICATION(3)==EQUATIONS_SET_CELLML_REAC_NO_SPLIT_REAC_DIFF_SUBTYPE .OR. &
            & equationsSet%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_REAC_DIFF_SUBTYPE) THEN
           sourceVector=>vectorMatrices%sourceVector
        ENDIF
        dynamicMapping=>vectorMapping%dynamicMapping
        FIELD_VARIABLE=>dynamicMapping%equationsMatrixToVarMaps(1)%variable
        FIELD_VAR_TYPE=FIELD_VARIABLE%variableType
        IF(stiffnessMatrix%updateMatrix.OR.dampingMatrix%updateMatrix.OR.rhsVector%updateVector) THEN
          !Loop over gauss points
          DO ng=1,QUADRATURE_SCHEME%numberOfGauss
            CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
              & geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            CALL Field_InterpolatedPointMetricsCalculate(GEOMETRIC_BASIS%numberOfXi,equations%interpolation% &
              & geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
              & materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            IF(USE_FIBRES) THEN
              CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
                & fibreInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              CALL Field_InterpolatedPointMetricsCalculate(FIBRE_BASIS%numberOfXi,equations%interpolation% &
                & fibreInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            ENDIF
            IF(equationsSet%SPECIFICATION(3)==EQUATIONS_SET_CELLML_REAC_NO_SPLIT_REAC_DIFF_SUBTYPE .OR. &
            & equationsSet%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_REAC_DIFF_SUBTYPE) THEN
              CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
                & sourceInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            ENDIF
            !Calculate RWG.
            RWG=equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%jacobian* &
              & QUADRATURE_SCHEME%gaussWeights(ng)
            !Calculate the diffusivity tensor
            DIFFUSIVITY=0.0_DP
            IF(USE_FIBRES) THEN
              !Calculate the diffusivity tensor in fibre coordinates
              CALL FlagError("Not implemented.",err,error,*999)
            ELSE
              !Use the diffusivity tensor in geometric coordinates
              DO nj=1,GEOMETRIC_VARIABLE%numberOfComponents !first three components of material field are the diffusivities
                DIFFUSIVITY(nj,nj)=equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(nj,1)
              ENDDO !nj
            ENDIF
            !Get the storage Coefficient, stored in the component after the diffusivities for each dimension
            componentIdx=GEOMETRIC_VARIABLE%numberOfComponents+1
            STORAGE_COEFFICIENT=equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(componentIdx,1)
            !Compute basis dPhi/dx terms
            DO nj=1,GEOMETRIC_VARIABLE%numberOfComponents
              DO ms=1,DEPENDENT_BASIS%numberOfElementParameters
                DPHIDX(nj,ms)=0.0_DP
                DO ni=1,DEPENDENT_BASIS%numberOfXi
                  DPHIDX(nj,ms)=DPHIDX(nj,ms)+ &
                    & QUADRATURE_SCHEME%gaussBasisFunctions(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)* &
                    & equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%dXidX(ni,nj)
                ENDDO !ni
              ENDDO !ms
            ENDDO !nj            
            !Loop over field components
            mhs=0          
            DO mh=1,FIELD_VARIABLE%numberOfComponents
              !Loop over element rows
              DO ms=1,DEPENDENT_BASIS%numberOfElementParameters
                mhs=mhs+1
                nhs=0
                !Loop over element columns
                DO nh=1,FIELD_VARIABLE%numberOfComponents
                  DO ns=1,DEPENDENT_BASIS%numberOfElementParameters
                    nhs=nhs+1
                    SUM=0.0_DP
                    IF(stiffnessMatrix%updateMatrix) THEN
                      DO ni=1,GEOMETRIC_VARIABLE%numberOfComponents
                        DO nj=1,GEOMETRIC_VARIABLE%numberOfComponents
                          SUM=SUM+DIFFUSIVITY(ni,nj)*DPHIDX(ni,mhs)*DPHIDX(nj,nhs)
                       ENDDO !nj
                      ENDDO !ni
                      stiffnessMatrix%elementMatrix%matrix(mhs,nhs)=stiffnessMatrix%elementMatrix%matrix(mhs,nhs)+(SUM*RWG)
                    ENDIF
                    IF(dampingMatrix%updateMatrix) THEN
                      dampingMatrix%elementMatrix%matrix(mhs,nhs)=dampingMatrix%elementMatrix%matrix(mhs,nhs)+ &
                        & QUADRATURE_SCHEME%gaussBasisFunctions(ms,NO_PART_DERIV,ng)* &
                        & QUADRATURE_SCHEME%gaussBasisFunctions(ns,NO_PART_DERIV,ng)*STORAGE_COEFFICIENT*RWG
                    ENDIF
                  ENDDO !ns
                ENDDO !nh

                IF(rhsVector%updateVector) rhsVector%elementVector%vector(mhs)=0.0_DP
              ENDDO !ms
            ENDDO !mh
            IF(equationsSet%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_REAC_DIFF_SUBTYPE) THEN
              IF(sourceVector%updateVector) THEN
                C_PARAM=equations%interpolation%sourceInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(1,NO_PART_DERIV)
                mhs=0
                DO mh=1,FIELD_VARIABLE%numberOfComponents
                !Loop over element rows
                  DO ms=1,DEPENDENT_BASIS%numberOfElementParameters
                    mhs=mhs+1
                    sourceVector%elementVector%vector(mhs)=sourceVector%elementVector%vector(mhs)+ &
                      & QUADRATURE_SCHEME%gaussBasisFunctions(ms,NO_PART_DERIV,ng)*C_PARAM*RWG
                  ENDDO
                ENDDO
              ENDIF
            ENDIF
            IF(rhsVector%updateVector) rhsVector%elementVector%vector(mhs)=0.0_DP
          ENDDO !ng
        ENDIF
        !Scale factor adjustment
        IF(dependentField%SCALINGS%scalingType/=FIELD_NO_SCALING) THEN
          CALL Field_InterpolationParametersScaleFactorsElementGet(ELEMENT_NUMBER,equations%interpolation% &
            & dependentInterpParameters(FIELD_VAR_TYPE)%ptr,err,error,*999)
          mhs=0          
          DO mh=1,FIELD_VARIABLE%numberOfComponents
            !Loop over element rows
            DO ms=1,DEPENDENT_BASIS%numberOfElementParameters
              mhs=mhs+1                    
              nhs=0
              IF(stiffnessMatrix%updateMatrix.OR.dampingMatrix%updateMatrix) THEN
                !Loop over element columns
                DO nh=1,FIELD_VARIABLE%numberOfComponents
                  DO ns=1,DEPENDENT_BASIS%numberOfElementParameters
                    nhs=nhs+1
                    IF(stiffnessMatrix%updateMatrix) THEN
                      stiffnessMatrix%elementMatrix%matrix(mhs,nhs)=stiffnessMatrix%elementMatrix%matrix(mhs,nhs)* &
                        & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%scaleFactors(ms,mh)* &
                        & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%scaleFactors(ns,nh)
                    ENDIF
                    IF(dampingMatrix%updateMatrix) THEN
                      dampingMatrix%elementMatrix%matrix(mhs,nhs)=dampingMatrix%elementMatrix%matrix(mhs,nhs)* &
                        & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%scaleFactors(ms,mh)* &
                        & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%scaleFactors(ns,nh)
                    ENDIF
                  ENDDO !ns
                ENDDO !nh
              ENDIF
              IF(rhsVector%updateVector) rhsVector%elementVector%vector(mhs)=rhsVector%elementVector%vector(mhs)* &
                & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%scaleFactors(ms,mh)
              IF(equationsSet%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_REAC_DIFF_SUBTYPE) THEN
                 IF(sourceVector%updateVector) sourceVector%elementVector%vector(mhs)= & 
                    & sourceVector%elementVector%vector(mhs)* &
                    & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%scaleFactors(ms,mh)
              ENDIF
            ENDDO !ms
          ENDDO !mh
        ENDIF
      ELSE
        CALL FlagError("Equations set equations is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
       
    EXITS("ReactionDiffusion_FiniteElementCalculate")
    RETURN
999 ERRORS("ReactionDiffusion_FiniteElementCalculate",err,error)
    EXITS("ReactionDiffusion_FiniteElementCalculate")
    RETURN 1
    
  END SUBROUTINE ReactionDiffusion_FiniteElementCalculate

  !
  !================================================================================================================================
  !
  !>Sets the problem specification for a reaction-diffusion problem.
  SUBROUTINE ReactionDiffusion_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER, INTENT(IN) :: problem !<A pointer to the problem to set the problem specification for.
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: problemSubtype
    TYPE(VARYING_STRING) :: localError

    ENTERS("ReactionDiffusion_ProblemSpecificationSet",err,error,*999)

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
    CASE(PROBLEM_CELLML_REAC_INTEG_REAC_DIFF_STRANG_SPLIT_SUBTYPE, &
      & PROBLEM_CELLML_REAC_EVAL_REAC_DIFF_NO_SPLIT_SUBTYPE, &
      & PROBLEM_CONSTANT_REAC_DIFF_NO_SPLIT_SUBTYPE)
      !ok
    CASE DEFAULT
      localError="The specified problem subtype of "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
        & " is not valid for a reaction-diffusion problem type of a classical problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    ALLOCATE(problem%specification(3),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
    problem%specification(1:3)=[PROBLEM_CLASSICAL_FIELD_CLASS,PROBLEM_REACTION_DIFFUSION_EQUATION_TYPE,problemSubtype]
 
    EXITS("ReactionDiffusion_ProblemSpecificationSet")
    RETURN
999 ERRORS("ReactionDiffusion_ProblemSpecificationSet",err,error)
    EXITS("ReactionDiffusion_ProblemSpecificationSet")
    RETURN 1
    
  END SUBROUTINE ReactionDiffusion_ProblemSpecificationSet

  !
  !================================================================================================================================
  !
  
  !>Sets up the reaction-diffusion problem.
  SUBROUTINE ReactionDiffusion_ProblemSetup(problem,problemSetup,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem set to setup a bioelectric domain equation on.
    TYPE(ProblemSetupType), INTENT(INOUT) :: problemSetup !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: pSpecification(3)
    TYPE(CellMLEquationsType), POINTER :: cellMLEquations
    TYPE(ControlLoopType), POINTER :: controlLoop,controlLoopRoot
    TYPE(SolverType), POINTER :: solver
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolversType), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError
   
    ENTERS("ReactionDiffusion_ProblemSetup",err,error,*999)

    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_CELLML_REAC_INTEG_REAC_DIFF_STRANG_SPLIT_SUBTYPE)
      !OK
    CASE(PROBLEM_CELLML_REAC_EVAL_REAC_DIFF_NO_SPLIT_SUBTYPE)
      !OK
    CASE(PROBLEM_CONSTANT_REAC_DIFF_NO_SPLIT_SUBTYPE)
      !OK
    CASE DEFAULT
      localError="The problem subtype of "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is invalid for a reaction-diffusion problem type of a classical problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
  
    SELECT CASE(problemSetup%setupType)
    CASE(PROBLEM_SETUP_INITIAL_TYPE)
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)
        !Do nothing????
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Do nothing????
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a reaction diffusion equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_SETUP_CONTROL_TYPE)
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)
        !Set up a time control loop
        NULLIFY(controlLoop)
        CALL ControlLoop_CreateStart(problem,controlLoop,err,error,*999)
        CALL ControlLoop_TypeSet(controlLoop,CONTROL_TIME_LOOP_TYPE,err,error,*999)
        CALL ControlLoop_LabelSet(controlLoop,"Time Loop",err,error,*999)
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
          & " is invalid for a reaction-diffusion equation."
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
        CASE(PROBLEM_CELLML_REAC_INTEG_REAC_DIFF_STRANG_SPLIT_SUBTYPE)
          CALL Solvers_NumberOfSolversSet(solvers,3,err,error,*999)
          !Set the first solver to be a differential-algebraic equations solver
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          CALL Solver_TypeSet(solver,SOLVER_DAE_TYPE,err,error,*999)
          CALL Solver_LabelSet(solver,"First ODE solver",err,error,*999)
          !Set solver defaults
          CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
          !Set the second solver to be a dynamic solver 
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
          CALL Solver_TypeSet(solver,SOLVER_DYNAMIC_TYPE,err,error,*999)
          CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
          CALL Solver_LabelSet(solver,"Parabolic solver",err,error,*999)
          !Set solver defaults
          CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
          CALL Solver_DynamicSchemeSet(solver,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,err,error,*999)
          CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
          !Set the third solver to be a differential-algebraic equations solver
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,3,solver,err,error,*999)
          CALL Solver_TypeSet(solver,SOLVER_DAE_TYPE,err,error,*999)
          CALL Solver_LabelSet(solver,"Second ODE solver",err,error,*999)
          !Set solver defaults
          CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
        CASE(PROBLEM_CELLML_REAC_EVAL_REAC_DIFF_NO_SPLIT_SUBTYPE)
          CALL Solvers_NumberOfSolversSet(solvers,2,err,error,*999)
          !Set the first solver to be a CELLML evaluator equations solver
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          CALL Solver_TypeSet(solver,SOLVER_CELLML_EVALUATOR_TYPE,err,error,*999)
          CALL Solver_LabelSet(solver,"Evaluator solver",err,error,*999)
          !Set solver defaults
          CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
          !Set the second solver to be a dynamic solver 
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
          CALL Solver_TypeSet(solver,SOLVER_DYNAMIC_TYPE,err,error,*999)
          CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
          CALL Solver_LabelSet(solver,"Parabolic solver",err,error,*999)
          !Set solver defaults
          CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
          CALL Solver_DynamicSchemeSet(solver,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,err,error,*999)
          CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
        CASE(PROBLEM_CONSTANT_REAC_DIFF_NO_SPLIT_SUBTYPE)
          CALL Solvers_NumberOfSolversSet(solvers,1,err,error,*999)
          !Set the solver to be a dynamic solver 
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          CALL Solver_TypeSet(solver,SOLVER_DYNAMIC_TYPE,err,error,*999)
          CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
          CALL Solver_LabelSet(solver,"Parabolic solver",err,error,*999)
          !Set solver defaults
          CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
          CALL Solver_DynamicSchemeSet(solver,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,err,error,*999)
          CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
        CASE DEFAULT
          localError="The problem subtype of "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
            & " is invalid for a reaction-diffusion problem type of a classical problem class."
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
          & " is invalid for a classical equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)
        !Get the control loop
        NULLIFY(controlLoopRoot)
        CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
        NULLIFY(controlLoop)
        CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
        SELECT CASE(pSpecification(3))
        CASE(PROBLEM_CELLML_REAC_INTEG_REAC_DIFF_STRANG_SPLIT_SUBTYPE)
          !Get the solver
          NULLIFY(solvers)
          CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
          !Create the solver equations for the second (parabolic) solver
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
          NULLIFY(solverEquations)
          CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
          CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_LINEAR,err,error,*999)
          CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
          CALL SolverEquations_TypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
        CASE(PROBLEM_CELLML_REAC_EVAL_REAC_DIFF_NO_SPLIT_SUBTYPE)
          !Get the solver
          NULLIFY(solvers)
          CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
          !Create the solver equations for the parabolic solver
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
          NULLIFY(solverEquations)
          CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
          CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
          CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
          CALL SolverEquations_TypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
        CASE(PROBLEM_CONSTANT_REAC_DIFF_NO_SPLIT_SUBTYPE)
          !Get the solver
          NULLIFY(solvers)
          CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
          !Create the solver equations for the parabolic solver
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          NULLIFY(solverEquations)
          CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
          CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_LINEAR,err,error,*999)
          CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
          CALL SolverEquations_TypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
        CASE DEFAULT
          localError="The problem subtype of "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
            & " is invalid for a reaction-diffusion problem type of a classical problem class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
     CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Get the control loop
        NULLIFY(controlLoopRoot)
        CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
        NULLIFY(controlLoop)
        CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
        NULLIFY(solvers)
        CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
        SELECT CASE(pSpecification(3))
        CASE(PROBLEM_CELLML_REAC_INTEG_REAC_DIFF_STRANG_SPLIT_SUBTYPE)
          !Get the solver equations for the second (parabolic) solver
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
          NULLIFY(solverEquations)
          CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
          !Finish the solver equations creation
          CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)             
        CASE(PROBLEM_CELLML_REAC_EVAL_REAC_DIFF_NO_SPLIT_SUBTYPE)
          !Get the solver equations for the parabolic solver
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
          NULLIFY(solverEquations)
          CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
          !Finish the solver equations creation
          CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)             
        CASE(PROBLEM_CONSTANT_REAC_DIFF_NO_SPLIT_SUBTYPE)
          !Get the solver equations for thE solver
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          NULLIFY(solverEquations)
          CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
          !Finish the solver equations creation
          CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)             
        CASE DEFAULT
          localError="The problem subtype of "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
            & " is invalid for a reaction-diffusion problem type of a classical problem class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a reaction-diffusion equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_SETUP_CELLML_EQUATIONS_TYPE)
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)
        !Get the control loop
        NULLIFY(controlLoopRoot)
        CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
        NULLIFY(controlLoop)
        CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
        !Get the solver
        NULLIFY(solvers)
        CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)        
        SELECT CASE(pSpecification(3))
        CASE(PROBLEM_CELLML_REAC_INTEG_REAC_DIFF_STRANG_SPLIT_SUBTYPE)
          !Create the CellML equations for the first DAE solver
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          NULLIFY(cellMLEquations)
          CALL CellMLEquations_CreateStart(solver,cellMLEquations,err,error,*999)
          !Set the time dependence
          CALL CellMLEquations_TimeDependenceTypeSet(cellMLEquations,CELLML_EQUATIONS_DYNAMIC,err,error,*999)
          !Set the linearity
          CALL CellMLEquations_LinearityTypeSet(cellMLEquations,CELLML_EQUATIONS_NONLINEAR,err,error,*999)
          !Create the CellML equations for the second DAE solver
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,3,solver,err,error,*999)
          NULLIFY(cellMLEquations)
          CALL CellMLEquations_CreateStart(solver,cellMLEquations,err,error,*999)
          !Set the time dependence
          CALL CellMLEquations_TimeDependenceTypeSet(cellMLEquations,CELLML_EQUATIONS_DYNAMIC,err,error,*999)
          !Set the linearity
          CALL CellMLEquations_LinearityTypeSet(cellMLEquations,CELLML_EQUATIONS_NONLINEAR,err,error,*999)
        CASE(PROBLEM_CELLML_REAC_EVAL_REAC_DIFF_NO_SPLIT_SUBTYPE)
          !Create the CellML equations for the first cellml evaluator solver
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          NULLIFY(cellMLEquations)
          CALL CellMLEquations_CreateStart(solver,cellMLEquations,err,error,*999)
          !Set the time dependence
          CALL CellMLEquations_TimeDependenceTypeSet(cellMLEquations,CELLML_EQUATIONS_DYNAMIC,err,error,*999)
          !Set the linearity
          CALL CellMLEquations_LinearityTypeSet(cellMLEquations,CELLML_EQUATIONS_NONLINEAR,err,error,*999)
        ENDIF
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Get the control loop
        NULLIFY(controlLoopRoot)
        CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
        NULLIFY(controlLoop)
        CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
        NULLIFY(solvers)
        CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
        SELECT CASE(pSpecification(3))
        CASE(PROBLEM_CELLML_REAC_INTEG_REAC_DIFF_STRANG_SPLIT_SUBTYPE)
          !Get the CellML equations for the first DAE solver
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          NULLIFY(cellMLEquations)
          CALL Solver_CellMLEquationsGet(solver,cellMLEquations,err,error,*999)
          !Finish the CellML equations creation
          CALL CellMLEquations_CreateFinish(cellMLEquations,err,error,*999)
          !Get the CellML equations for the second DAE solver
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,3,solver,err,error,*999)
          NULLIFY(cellMLEquations)
          CALL Solver_CellMLEquationsGet(solver,cellMLEquations,err,error,*999)
          !Finish the CellML equations creation
          CALL CellMLEquations_CreateFinish(cellMLEquations,err,error,*999)
        CASE(PROBLEM_CELLML_REAC_EVAL_REAC_DIFF_NO_SPLIT_SUBTYPE)
          !Get the CellML equations for the first evaluator solver
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          NULLIFY(cellMLEquations)
          CALL Solver_CellMLEquationsGet(solver,cellMLEquations,err,error,*999)
          !Finish the CellML equations creation
          CALL CellMLEquations_CreateFinish(cellMLEquations,err,error,*999)
        CASE DEFAULT
          localError="The problem subtype of "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
            & " is invalid for a reaction-diffusion problem type of a classical problem class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for reaction-diffusion equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
        & " is invalid for areaction-diffusion equation."
      CALL FlagError(localError,err,error,*999)
    END SELECT
     
    EXITS("ReactionDiffusion_ProblemSetup")
    RETURN
999 ERRORSEXITS("ReactionDiffusion_ProblemSetup",err,error)
    RETURN 1
    
  END SUBROUTINE ReactionDiffusion_ProblemSetup
  
  !
  !================================================================================================================================
  !
  
  !>Performs pre-solve actions for reaction-diffusion problems.
  SUBROUTINE ReactionDiffusion_PreSolve(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver to perform the pre-solve actions for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: pSpecification(3),solverNumber
    REAL(DP) :: currentTime,timeIncrement
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemType), POINTER :: problem
     TYPE(VARYING_STRING) :: localError

    ENTERS("ReactionDiffusion_PreSolve",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
  
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_CELLML_REAC_INTEG_REAC_DIFF_STRANG_SPLIT_SUBTYPE)
      CALL ControlLoop_CurrentTimesGet(controlLoop,currentTime,timeIncrement,err,error,*999)
      CALL Solver_GlobalNumberGet(solver,solverNumber,err,error,*999)
      SELECT CASE(solverNumber)
      CASE(1)
        CALL Solver_DAETimesSet(solver,currentTime,currentTime+timeIncrement/2.0_DP,err,error,*999)
      CASE(2)
        !Do nothing
      CASE(3)
        CALL Solver_DAETimesSet(solver,currentTime,currentTime+timeIncrement/2.0_DP,err,error,*999)
      CASE DEFAULT
        localError="The solver global number of "//TRIM(NumberToVString(solverNumber,"*",err,error))// &
          & " is invalid for a Strang split reaction-diffusion problem."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_CONSTANT_REAC_DIFF_NO_SPLIT_SUBTYPE)
      !No splitting, therefore entire problem is solved as a dynamic one, with 1 solver nothing to do.
    CASE(PROBLEM_CELLML_REAC_EVAL_REAC_DIFF_NO_SPLIT_SUBTYPE)
      !No splitting, therefore entire problem is solved as a dynamic one, with 1 solver nothing to do.
    CASE DEFAULT
      localError="The problem subtype of "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is invalid for a reaction-diffusion problem type."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("ReactionDiffusion_PreSolve")
    RETURN
999 ERRORSEXITS("ReactionDiffusion_PreSolve",err,error)
    RETURN 1
    
  END SUBROUTINE ReactionDiffusion_PreSolve

  !
  !================================================================================================================================
  !

  !>Performs post-solve operations for reaction-diffusion problems. 
  SUBROUTINE ReactionDiffusion_PostSolve(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: pSpecification(3),solverNumber
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolversType), POINTER :: solvers
    TYPE(SolverType), POINTER :: pdeSolver
    TYPE(VARYING_STRING) :: localError

    ENTERS("ReactionDiffusion_PostSolve",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
     
    SELECT CASE(pSpecification(3))      
    CASE(PROBLEM_CELLML_REAC_INTEG_REAC_DIFF_STRANG_SPLIT_SUBTYPE)
      CALL Solver_GlobalNumberGet(solver,solverNumber,err,error,*999)
      SELECT CASE(solverNumber)
      CASE(1)        
        !do nothing
      CASE(2)        
        !do nothing
      CASE(3)
        !OUTPUT SOLUTIONS AT EACH TIME STEP - should probably change this bit below to output 
        !mesh solutions directly from the 3rd solver itself rather than by getting the 2nd solver.
        !I just don't know how to work with cellml_equations to do this.
        NULLIFY(solvers)
        CALL Solver_SolversGet(solver,solvers,err,error,*999)
        NULLIFY(pdeSolver)
        CALL Solvers_SolverGet(solvers,2,pdeSolver,err,error,*999)
        CALL ReactionDiffusion_PostSolveOutputData(controlLoop,pdeSolver,err,error,*999)
      CASE DEFAULT
        localError="The solver global number of "//TRIM(NumberToVString(solverNumber,"*",err,error))// &
          & " is invalid for a Strang split reaction-diffusion problem."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE (PROBLEM_CELLML_REAC_EVAL_REAC_DIFF_NO_SPLIT_SUBTYPE)
      !do nothing - time output not implemented
    CASE (PROBLEM_CONSTANT_REAC_DIFF_NO_SPLIT_SUBTYPE)
      !OUTPUT SOLUTIONS AT TIME STEP
      CALL ReactionDiffusion_PostSolveOutputData(controlLoop,solver,err,error,*999)
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a reaction diffusion type of a classical field problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("ReactionDiffusion_PostSolve")
    RETURN
999 ERRORSEXITS("ReactionDiffusion_PostSolve",err,error)
    RETURN 1
    
  END SUBROUTINE ReactionDiffusion_PostSolve
  
  !   
  !================================================================================================================================
  !

  SUBROUTINE ReactionDiffusion_PostSolveOutputData(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver!<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables 
    INTEGER(INTG) :: equationsSetIdx,currentIteration,esNumber,inputFrequency,maxDigits,myGroupComputationNodeNumber, &
      & outputFrequency,pSpecification(3)
    REAL(DP) :: currentTime,startTime,stopTime,timeIncrement
    LOGICAL :: exportExelem
    CHARACTER(30) :: file
    CHARACTER(30) :: outputFile
    CHARACTER(100) :: fmt, tempFmt
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(ProblemType), POINTER :: problem
    TYPE(RegionType), POINTER :: region
    TYPE(VARYING_STRING) :: localError
    TYPE(WorkGroupType), POINTER :: workGroup

    ENTERS("ReactionDiffusion_PostSolveOutputData",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_CELLML_REAC_INTEG_REAC_DIFF_STRANG_SPLIT_SUBTYPE, &
      & PROBLEM_CONSTANT_REAC_DIFF_NO_SPLIT_SUBTYPE)
      CALL ControlLoop_CurrentTimeInformationGet(controlLoop,currentTime,timeIncrement,startTime,stopTime,currentIteration, &
        & outputFrequency,inputFrequency,err,error,*999)
      CALL Solver_OutputTypeGet(solver,outputType,err,error,*999)
      NULLIFY(solverEquations)
      CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      !Make sure the equations sets are up to date
      CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
      DO equationsSetIdx=1,numberOfEquationsSets
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,solverMapping,err,error,*999)
        NULLIFY(region)
        CALL EquationsSet_RegionGet(equationsSet,region,err,error,*999)
        CALL EquationsSet_GlobalNumberGet(equationsSet,esNumber,err,error,*999)
        NULLIFY(workGroup)
        CALL Solver_WorkGroupGet(solver,workGroup,err,error,*999)
        CALL WorkGroup_GroupNodeNumberGet(workGroup,myGroupComputationNodeNumber,err,error,*999)
        maxDigits=FLOOR(LOG10((stopTime-startTime)/timeIncrement))+1
        IF(outputFrequency>0) THEN
          IF(MOD(currentIteration,outputFrequency)==0) THEN
            IF(currentTime<=stopTime) THEN
              IF(numberOfEquationsSets==1) THEN
                WRITE(tempFmt,'("I",I0,".",I0)') maxDigits,maxDigits
                !100 FORMAT 
                fmt = TRIM(tempFmt)
                WRITE(tempFmt,'(A2,A38,A20,A2)') "(", '"TimeStepSpec_1.part",I2.2,".",',fmt,")"
                fmt = TRIM(tempFmt)
                WRITE(outputFile,fmt) myGroupComputationNodeNumber,currentIteration
              ELSE
                WRITE(tempFmt,'("I",I0,".",I0)') maxDigits,maxDigits
                !200 FORMAT 
                fmt = TRIM(tempFmt)
                WRITE(tempFmt,'(A2,A38,A20,A2)') "(", '"TimeStepSpec_",I0,".part",I2.2,".",',fmt,")"
                fmt = TRIM(tempFmt)
                WRITE(outputFile,fmt) equationsSetIdx, myGroupComputationNodeNumber,currentIteration
              ENDIF
              WRITE(*,*) outputFile
              file=TRIM(outputFile)
              IF(outputType>=SOLVER_PROGRESS_OUTPUT) THEN
                CALL WriteString(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
                CALL WriteString(GENERAL_OUTPUT_TYPE,"Now export fields... ",err,error,*999)
              ENDIF
              exportExelem = .FALSE.
              IF(currentIteration==0) THEN
                IF(equationsSetIdx==1) exportExelem = .TRUE.
              ENDIF
              CALL REACTION_DIFFUSION_IO_WRITE_CMGUI(region,esNumber,file,exportExelem,err,error,*999)
            ENDIF
          ENDIF
        ENDIF
      ENDDO !equationsSetIdx
    CASE(PROBLEM_CELLML_REAC_EVAL_REAC_DIFF_NO_SPLIT_SUBTYPE)
      ! do nothing ???
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for an advection-diffusion equation type of a classical field problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("ReactionDiffusion_PostSolveOutputData")
    RETURN
999 ERRORSEXITS("ReactionDiffusion_PostSolveOutputData",err,error)
    RETURN 1
    
  END SUBROUTINE ReactionDiffusion_PostSolveOutputData
  
  !
  !================================================================================================================================
  !

  SUBROUTINE ReactionDiffusion_PostLoop(controlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatrixType), POINTER :: dampingMatrix,stiffnessMatrix
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(EquationsVectorType), POINTER :: vectorEquations 
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverType), POINTER :: solver 
    TYPE(SolversType), POINTER :: solvers
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
   
    ENTERS("ReactionDiffusion_PostLoop",err,error,*999)

    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,pSpecification,err,error,*999)
    
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_CELLML_REAC_INTEG_REAC_DIFF_STRANG_SPLIT_SUBTYPE)
      NULLIFY(solvers)
      CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
      NULLIFY(solver)
      CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
      NULLIFY(solverEquations)
      CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
      NULLIFY(equations)
      CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
      NULLIFY(vectorEquations)
      CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
      NULLIFY(vectorMatrices)
      CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
      NULLIFY(dynamicMatrices)
      CALL EquationsMatricesVector_DynamicMatrices(vectorMatrices,dynamicMatrices,err,error,*999)
      NULLIFY(stiffnessMatrix)
      CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,1,stiffnessMatrix,err,error,*999)
      NULLIFY(dampingMatrix)
      CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,2,dampingMatrix,err,error,*999)
      stiffnessMatrix%updateMatrix = .FALSE.
      dampingMatrix%updateMatrix = .FALSE.
    CASE DEFAULT
      !do nothing
    END SELECT
    
    EXITS("ReactionDiffusion_PostLoop")
    RETURN
999 ERRORSEXITS("ReactionDiffusion_PostLoop",err,error)
    RETURN 1
  END SUBROUTINE ReactionDiffusion_PostLoop
  
  !
  !================================================================================================================================
  !
  
END MODULE ReactionDiffusionEquationsRoutines
