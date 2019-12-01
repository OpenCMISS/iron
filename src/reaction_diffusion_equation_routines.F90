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
!> Contributor(s): Vijay Rajagopal,Chris Bradley
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
MODULE REACTION_DIFFUSION_EQUATION_ROUTINES

  USE BaseRoutines
  USE BasisRoutines
  USE BasisAccessRoutines
  USE BOUNDARY_CONDITIONS_ROUTINES
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
  USE EquationsSetConstants
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
  USE PROBLEM_CONSTANTS
  USE Strings
  USE SOLVER_ROUTINES
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
  
  PUBLIC REACTION_DIFFUSION_EQUATION_PROBLEM_SETUP

  PUBLIC ReactionDiffusion_ProblemSpecificationSet

  PUBLIC ReactionDiffusion_PostSolve

  PUBLIC REACTION_DIFFUSION_CONTROL_LOOP_POST_LOOP


CONTAINS

  !
  !================================================================================================================================
  !

  !>Sets up the reaction diffusion equation type of a classical equations set class.
  SUBROUTINE ReactionDiffusion_EquationsSetSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup a bioelectric domain equation on.
    TYPE(EquationsSetSetupType), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,DIMENSION_MULTIPLIER,GEOMETRIC_COMPONENT_NUMBER,GEOMETRIC_SCALING_TYPE, &
      & numberOfDimensions,NUMBER_OF_MATERIALS_COMPONENTS,GEOMETRIC_MESH_COMPONENT
    TYPE(DecompositionType), POINTER :: GEOMETRIC_DECOMPOSITION
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsSetMaterialsType), POINTER :: EQUATIONS_MATERIALS
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("ReactionDiffusion_EquationsSetSetup",err,error,*999)

    NULLIFY(EQUATIONS)
    NULLIFY(vectorMapping)
    NULLIFY(vectorMatrices)
    NULLIFY(GEOMETRIC_DECOMPOSITION)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a reaction-diffusion type equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET_SETUP%setupType)
      CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
        SELECT CASE(EQUATIONS_SET_SETUP%actionType)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
          CALL ReactionDiffusion_EquationsSetSolutionMethodSet(EQUATIONS_SET, &
           & EQUATIONS_SET_FEM_SOLUTION_METHOD,err,error,*999)
        CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
!!Todo: CHECK VALID SETUP
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%setupType,"*",err,error))// &
            & " is invalid for a reaction diffusion domain equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
        !\todo Check geometric dimension
      CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
        SELECT CASE(EQUATIONS_SET_SETUP%actionType)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
          SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
          CASE(EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE, &
           & EQUATIONS_SET_CONSTANT_REAC_DIFF_SUBTYPE)
            IF(EQUATIONS_SET%DEPENDENT%dependentFieldAutoCreated) THEN
              !Create the auto created dependent field
              CALL Field_CreateStart(EQUATIONS_SET_SETUP%fieldUserNumber,EQUATIONS_SET%REGION,EQUATIONS_SET%DEPENDENT% &
                & dependentField,err,error,*999)
              CALL Field_LabelSet(EQUATIONS_SET%DEPENDENT%dependentField,"Dependent Field",err,error,*999)
              CALL Field_TypeSetAndLock(EQUATIONS_SET%DEPENDENT%dependentField,FIELD_GENERAL_TYPE,err,error,*999)
              CALL Field_DependentTypeSetAndLock(EQUATIONS_SET%DEPENDENT%dependentField,FIELD_DEPENDENT_TYPE,err,error,*999)
              CALL Field_DecompositionGet(EQUATIONS_SET%GEOMETRY%geometricField,GEOMETRIC_DECOMPOSITION,err,error,*999)
              CALL Field_DecompositionSetAndLock(EQUATIONS_SET%DEPENDENT%dependentField,GEOMETRIC_DECOMPOSITION, &
                & err,error,*999)
              CALL Field_GeometricFieldSetAndLock(EQUATIONS_SET%DEPENDENT%dependentField,EQUATIONS_SET%GEOMETRY% &
                & geometricField,err,error,*999)
              CALL Field_NumberOfVariablesSetAndLock(EQUATIONS_SET%DEPENDENT%dependentField,2,err,error,*999)
              CALL Field_VariableTypesSetAndLock(EQUATIONS_SET%DEPENDENT%dependentField,[FIELD_U_VARIABLE_TYPE, &
                & FIELD_DELUDELN_VARIABLE_TYPE],err,error,*999)
              !CALL Field_VariableLabelSet
              CALL Field_DimensionSetAndLock(EQUATIONS_SET%DEPENDENT%dependentField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DimensionSetAndLock(EQUATIONS_SET%DEPENDENT%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DataTypeSetAndLock(EQUATIONS_SET%DEPENDENT%dependentField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              CALL Field_DataTypeSetAndLock(EQUATIONS_SET%DEPENDENT%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              CALL Field_NumberOfComponentsSetAndLock(EQUATIONS_SET%DEPENDENT%dependentField,FIELD_U_VARIABLE_TYPE,1, &
                & err,error,*999)
              CALL Field_NumberOfComponentsSetAndLock(EQUATIONS_SET%DEPENDENT%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                & err,error,*999)
              !Default to the geometric interpolation setup
              CALL Field_ComponentMeshComponentGet(EQUATIONS_SET%GEOMETRY%geometricField,FIELD_U_VARIABLE_TYPE,1, &
                & GEOMETRIC_MESH_COMPONENT,err,error,*999)
              CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%DEPENDENT%dependentField,FIELD_U_VARIABLE_TYPE,1, &
                & GEOMETRIC_MESH_COMPONENT,err,error,*999)
              CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%DEPENDENT%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                & GEOMETRIC_MESH_COMPONENT,err,error,*999)
              SELECT CASE(EQUATIONS_SET%solutionMethod)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL Field_ComponentInterpolationSetAndLock(EQUATIONS_SET%DEPENDENT%dependentField, &
                  & FIELD_U_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationSetAndLock(EQUATIONS_SET%DEPENDENT%dependentField, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                !Default the scaling to the geometric field scaling
                CALL Field_ScalingTypeGet(EQUATIONS_SET%GEOMETRY%geometricField,GEOMETRIC_SCALING_TYPE,err,error,*999)
                CALL Field_ScalingTypeSet(EQUATIONS_SET%DEPENDENT%dependentField,GEOMETRIC_SCALING_TYPE,err,error,*999)
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
                localError="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%solutionMethod,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ELSE
              !Check the user specified field
              CALL Field_TypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
              CALL Field_DependentTypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
              CALL Field_NumberOfVariablesCheck(EQUATIONS_SET_SETUP%FIELD,2,err,error,*999)
              CALL Field_VariableTypesCheck(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE], &
                & err,error,*999)
              CALL Field_DimensionCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DimensionCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                & err,error,*999)
              CALL Field_DataTypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_DataTypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_NumberOfComponentsCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1,err,error,*999)
              CALL Field_NumberOfComponentsCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1,err,error,*999)
              SELECT CASE(EQUATIONS_SET%solutionMethod)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
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
                localError="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%solutionMethod,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ENDIF
          CASE(EQUATIONS_SET_CELLML_REAC_NO_SPLIT_REAC_DIFF_SUBTYPE)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The equation set subtype of "//TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
              & " is invalid for a reaction diffusion equation set class."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
          IF(EQUATIONS_SET%DEPENDENT%dependentFieldAutoCreated) THEN
            CALL Field_CreateFinish(EQUATIONS_SET%DEPENDENT%dependentField,err,error,*999)
          ENDIF
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%setupType,"*",err,error))// &
            & " is invalid for a reaction diffusion equation"
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
        SELECT CASE(EQUATIONS_SET_SETUP%actionType)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
          CALL EquationsSet_AssertDependentIsFinished(EQUATIONS_SET,err,error,*999)
          EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
          IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
            IF(EQUATIONS_MATERIALS%materialsFieldAutoCreated) THEN
              !Create the auto created materials field
              CALL Field_CreateStart(EQUATIONS_SET_SETUP%fieldUserNumber,EQUATIONS_SET%REGION,EQUATIONS_MATERIALS% &
                & materialsField,err,error,*999)
              CALL Field_TypeSetAndLock(EQUATIONS_MATERIALS%materialsField,FIELD_MATERIAL_TYPE,err,error,*999)
              CALL Field_DependentTypeSetAndLock(EQUATIONS_MATERIALS%materialsField,FIELD_INDEPENDENT_TYPE,err,error,*999)
              CALL Field_DecompositionGet(EQUATIONS_SET%GEOMETRY%geometricField,GEOMETRIC_DECOMPOSITION,err,error,*999)
              CALL Field_DecompositionSetAndLock(EQUATIONS_MATERIALS%materialsField,GEOMETRIC_DECOMPOSITION, &
                & err,error,*999)
              CALL Field_GeometricFieldSetAndLock(EQUATIONS_MATERIALS%materialsField,EQUATIONS_SET%GEOMETRY% &
                & geometricField,err,error,*999)
              CALL Field_NumberOfVariablesSetAndLock(EQUATIONS_MATERIALS%materialsField,1,err,error,*999)
              CALL Field_VariableTypesSetAndLock(EQUATIONS_MATERIALS%materialsField,[FIELD_U_VARIABLE_TYPE], &
                & err,error,*999)
              CALL Field_DimensionSetAndLock(EQUATIONS_MATERIALS%materialsField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DataTypeSetAndLock(EQUATIONS_MATERIALS%materialsField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              CALL Field_NumberOfComponentsGet(EQUATIONS_SET%GEOMETRY%geometricField,FIELD_U_VARIABLE_TYPE, &
                & numberOfDimensions,err,error,*999)
              IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE) THEN
                !Reaction Diffusion. Materials field components are 1 diffusion coeff for each dimension
                !plus one for the storage coefficient in alpha(delC/delt) = Div(-kgradC)+cellmlRC
                NUMBER_OF_MATERIALS_COMPONENTS=numberOfDimensions+1
                DIMENSION_MULTIPLIER=1
              ELSEIF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_REAC_DIFF_SUBTYPE) THEN
                !Constant reaction + diffusion. Materials field has 1 diffuse coeff for each dimension
                !plus one for the storage coefficient om  alpha(delC/delt) = Div(-kgradC)+const(x)_source
                NUMBER_OF_MATERIALS_COMPONENTS=numberOfDimensions+1
                DIMENSION_MULTIPLIER=1
              ENDIF
              !Set the number of materials components
              CALL Field_NumberOfComponentsSetAndLock(EQUATIONS_MATERIALS%materialsField,FIELD_U_VARIABLE_TYPE, &
                & NUMBER_OF_MATERIALS_COMPONENTS,err,error,*999)
              !Default the first three materials components for diffusivity param to the first component geometric interpolation with const interpolation
              DO component_idx=1,numberOfDimensions
                CALL Field_ComponentMeshComponentGet(EQUATIONS_SET%GEOMETRY%geometricField,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_COMPONENT_NUMBER,err,error,*999)
                CALL Field_ComponentMeshComponentSet(EQUATIONS_MATERIALS%materialsField,FIELD_U_VARIABLE_TYPE, &
                  & component_idx,GEOMETRIC_COMPONENT_NUMBER,err,error,*999)
                CALL Field_ComponentInterpolationSet(EQUATIONS_MATERIALS%materialsField,FIELD_U_VARIABLE_TYPE, &
                  & component_idx,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
              ENDDO !components_idx
              !Default the storage co-efficient to the first geometric interpolation setup with constant interpolation
              IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE) THEN
                component_idx=NUMBER_OF_MATERIALS_COMPONENTS
                CALL Field_ComponentMeshComponentGet(EQUATIONS_SET%GEOMETRY%geometricField,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_COMPONENT_NUMBER,err,error,*999)
                CALL Field_ComponentMeshComponentSet(EQUATIONS_MATERIALS%materialsField,FIELD_U_VARIABLE_TYPE, &
                  & component_idx,GEOMETRIC_COMPONENT_NUMBER,err,error,*999)
                CALL Field_ComponentInterpolationSet(EQUATIONS_MATERIALS%materialsField,FIELD_U_VARIABLE_TYPE, &
                  & component_idx,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
              ENDIF
              !Default the field scaling to that of the geometric field
              CALL Field_ScalingTypeGet(EQUATIONS_SET%GEOMETRY%geometricField,GEOMETRIC_SCALING_TYPE,err,error,*999)
              CALL Field_ScalingTypeSet(EQUATIONS_MATERIALS%materialsField,GEOMETRIC_SCALING_TYPE,err,error,*999)
            ELSE
              !Check the user specified field
              CALL Field_TypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_MATERIAL_TYPE,err,error,*999)
              CALL Field_DependentTypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
              CALL Field_NumberOfVariablesCheck(EQUATIONS_SET_SETUP%FIELD,1,err,error,*999)
              CALL Field_VariableTypesCheck(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],err,error,*999)
              CALL Field_DimensionCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                & err,error,*999)
              CALL Field_DataTypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_NumberOfComponentsGet(EQUATIONS_SET%GEOMETRY%geometricField,FIELD_U_VARIABLE_TYPE, &
                & numberOfDimensions,err,error,*999)
              IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE .OR. &
                & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_REAC_DIFF_SUBTYPE) THEN
                !Reaction Diffusion with cellml. Materials field components are 1 for storage coeff plus one for each dimension i.e., k
                CALL Field_NumberOfComponentsCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,numberOfDimensions+1, &
                  & err,error,*999)
              ENDIF
            ENDIF
          ELSE
            CALL FlagError("Equations set materials is not associated.",err,error,*999)
          ENDIF
        CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
          EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
          IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
            IF(EQUATIONS_MATERIALS%materialsFieldAutoCreated) THEN
              !Finish creating the materials field
              CALL Field_CreateFinish(EQUATIONS_MATERIALS%materialsField,err,error,*999)
              !Set the default values for the materials field
              CALL Field_NumberOfComponentsGet(EQUATIONS_SET%GEOMETRY%geometricField,FIELD_U_VARIABLE_TYPE, &
                & numberOfDimensions,err,error,*999)
              IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE .OR. &
               & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_REAC_DIFF_SUBTYPE) THEN
                !Reaction Diffusion with cellml. Materials field components are 1 plus one for each dimension i.e.,storage coeff, and k.
                NUMBER_OF_MATERIALS_COMPONENTS=numberOfDimensions+1
                DIMENSION_MULTIPLIER=1
              ENDIF
              !set the diffusion coefficients to be 1.0
              DO component_idx=1,numberOfDimensions
                CALL Field_ComponentValuesInitialise(EQUATIONS_MATERIALS%materialsField,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,component_idx,1.0_DP,err,error,*999)
              ENDDO !component_idx
               !Now set storage-coefficient
              component_idx=numberOfDimensions+1
              CALL Field_ComponentValuesInitialise(EQUATIONS_MATERIALS%materialsField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,component_idx,1.0_DP,err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Equations set materials is not associated.",err,error,*999)
          ENDIF
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%setupType,"*",err,error))// &
            & " is invalid for a reaction diffusion equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
        SELECT CASE(EQUATIONS_SET_SETUP%actionType)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
          IF(ASSOCIATED(EQUATIONS_SET%MATERIALS)) THEN
            CALL EquationsSet_AssertMaterialsIsFinished(EQUATIONS_SET,err,error,*999)
            IF(ASSOCIATED(EQUATIONS_SET%SOURCE)) THEN
              IF(EQUATIONS_SET%SOURCE%sourceFieldAutoCreated) THEN
                !Create the auto created source field
                !Start field creation with name 'SOURCE_FIELD'
                CALL Field_CreateStart(EQUATIONS_SET_SETUP%fieldUserNumber,EQUATIONS_SET%REGION, &
                  & EQUATIONS_SET%SOURCE%sourceField,err,error,*999)
                !Create a general field
                CALL Field_TypeSetAndLock(EQUATIONS_SET%SOURCE%sourceField,FIELD_GENERAL_TYPE,err,error,*999)
                !Label the field
                CALL Field_LabelSetAndLock(EQUATIONS_SET%SOURCE%sourceField,"Source Field",err,error,*999)
                !Set the dependent type
                CALL Field_DependentTypeSetAndLock(EQUATIONS_SET%SOURCE%sourceField,FIELD_INDEPENDENT_TYPE,err,error,*999)
                !Set the field decomposition to be that of the geometric decomposition
                CALL Field_DecompositionGet(EQUATIONS_SET%GEOMETRY%geometricField,GEOMETRIC_DECOMPOSITION,err,error,*999)
                CALL Field_DecompositionSetAndLock(EQUATIONS_SET%SOURCE%sourceField,GEOMETRIC_DECOMPOSITION, &
                  & err,error,*999)
                !Set the geometric field
                CALL Field_GeometricFieldSetAndLock(EQUATIONS_SET%SOURCE%sourceField,EQUATIONS_SET% & 
                  & GEOMETRY%geometricField,err,error,*999)
                !Set the field variables.
                CALL Field_NumberOfVariablesSetAndLock(EQUATIONS_SET%SOURCE%sourceField,1,err,error,*999)
                CALL Field_VariableTypesSetAndLock(EQUATIONS_SET%SOURCE%sourceField,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                !Set the dimension
                CALL Field_DimensionSetAndLock(EQUATIONS_SET%SOURCE%sourceField,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
                !Set the data type
                CALL Field_DataTypeSetAndLock(EQUATIONS_SET%SOURCE%sourceField,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                !Set the number of components to one
                CALL Field_NumberOfComponentsSetAndLock(EQUATIONS_SET%SOURCE%sourceField,FIELD_U_VARIABLE_TYPE,1, &
                  & err,error,*999)
                !Get the geometric mesh component
                CALL Field_ComponentMeshComponentGet(EQUATIONS_SET%GEOMETRY%geometricField,FIELD_U_VARIABLE_TYPE, & 
                  & 1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                !Default to the geometric interpolation setup
                CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%SOURCE%sourceField,FIELD_U_VARIABLE_TYPE,1, &
                  & GEOMETRIC_MESH_COMPONENT,err,error,*999)            
                !Specify the interpolation to be same as geometric interpolation
                SELECT CASE(EQUATIONS_SET%solutionMethod)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  CALL Field_ComponentInterpolationSetAndLock(EQUATIONS_SET%SOURCE%sourceField,FIELD_U_VARIABLE_TYPE,1, &
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
                  localError="The solution method of " &
                    & //TRIM(NumberToVString(EQUATIONS_SET%solutionMethod,"*",err,error))// " is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
                !Set the scaling to be the same as the geometric field
                CALL Field_ScalingTypeGet(EQUATIONS_SET%GEOMETRY%geometricField,GEOMETRIC_SCALING_TYPE, &
                  & err,error,*999)
                CALL Field_ScalingTypeSet(EQUATIONS_SET%SOURCE%sourceField,GEOMETRIC_SCALING_TYPE, &
                  & err,error,*999)
              ELSE
                !Check the user specified field
                CALL Field_TypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                CALL Field_DependentTypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                CALL Field_NumberOfVariablesCheck(EQUATIONS_SET_SETUP%FIELD,1,err,error,*999)
                CALL Field_VariableTypesCheck(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                CALL Field_DimensionCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL Field_DataTypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL Field_NumberOfComponentsCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1,err,error,*999)
                
                SELECT CASE(EQUATIONS_SET%solutionMethod)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
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
                  localError="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%solutionMethod, &
                    &"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
                
              ENDIF
            ELSE
              CALL FlagError("Equations set source is not associated.",err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Equations set materials is not associated.",err,error,*999)
          ENDIF
        CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
          IF(EQUATIONS_SET%SOURCE%sourceFieldAutoCreated) THEN
            CALL Field_CreateFinish(EQUATIONS_SET%SOURCE%sourceField,err,error,*999)
          ENDIF
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%setupType,"*",err,error))// &
            & " is invalid for a reaction diffusion equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
        SELECT CASE(EQUATIONS_SET_SETUP%actionType)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%setupType,"*",err,error))// &
            & " is invalid for a reaction diffusion equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
        SELECT CASE(EQUATIONS_SET_SETUP%actionType)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
          CALL EquationsSet_AssertSourceIsFinished(EQUATIONS_SET,err,error,*999)
          !Create the equations
          CALL Equations_CreateStart(EQUATIONS_SET,EQUATIONS,err,error,*999)
          SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
          CASE(EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE)
            CALL Equations_LinearityTypeSet(EQUATIONS,EQUATIONS_LINEAR,err,error,*999)
          CASE(EQUATIONS_SET_CONSTANT_REAC_DIFF_SUBTYPE)
            CALL Equations_LinearityTypeSet(EQUATIONS,EQUATIONS_LINEAR,err,error,*999)
          CASE(EQUATIONS_SET_CELLML_REAC_NO_SPLIT_REAC_DIFF_SUBTYPE)
            CALL Equations_LinearityTypeSet(EQUATIONS,EQUATIONS_NONLINEAR,err,error,*999)
          CASE DEFAULT
            localError="The equations matrices linearity set up of "// &
              & TRIM(NumberToVString(EQUATIONS%sparsityType,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
          CALL Equations_TimeDependenceTypeSet(EQUATIONS,EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
        CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
          SELECT CASE(EQUATIONS_SET%solutionMethod)
          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
            !Finish the creation of the equations
            CALL EquationsSet_EquationsGet(EQUATIONS_SET,equations,err,error,*999)
            CALL Equations_CreateFinish(equations,err,error,*999)
            NULLIFY(vectorEquations)
            CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
            !Create the equations mapping.
            CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_DELUDELN_VARIABLE_TYPE,vectorMapping,err,error,*999)
            CALL EquationsMapping_DynamicMatricesSet(vectorMapping,.TRUE.,.TRUE.,err,error,*999)
            CALL EquationsMapping_DynamicVariableTypeSet(vectorMapping,FIELD_U_VARIABLE_TYPE,err,error,*999)
            CALL EquationsMapping_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE,err,error,*999)

            CALL EquationsMapping_SourceVariableTypeSet(vectorMapping,FIELD_U_VARIABLE_TYPE,err,error,*999)
            CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
            !Create the equations matrices
            CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
            !Set up matrix storage and structure
            IF(EQUATIONS%lumpingType==EQUATIONS_LUMPED_MATRICES) THEN
              !Set up lumping
              CALL EquationsMatrices_DynamicLumpingTypeSet(vectorMatrices, &
                & [EQUATIONS_MATRIX_UNLUMPED,EQUATIONS_MATRIX_LUMPED],err,error,*999)
              CALL EquationsMatrices_DynamicStorageTypeSet(vectorMatrices, &
                & [DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE], &
                & err,error,*999)
              CALL EquationsMatrices_DynamicStructureTypeSet(vectorMatrices, &
                [EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_DIAGONAL_STRUCTURE],err,error,*999)
            ELSE
              SELECT CASE(EQUATIONS%sparsityType)
              CASE(EQUATIONS_MATRICES_FULL_MATRICES) 
                CALL EquationsMatrices_DynamicStorageTypeSet(vectorMatrices, &
                  & [DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE,DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE],err,error,*999)
              CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                CALL EquationsMatrices_DynamicStorageTypeSet(vectorMatrices, &
                  & [DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE], &
                  & err,error,*999)
                CALL EquationsMatrices_DynamicStructureTypeSet(vectorMatrices, &
                  [EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_FEM_STRUCTURE],err,error,*999)                  
              CASE DEFAULT
                localError="The equations matrices sparsity type of "// &
                  & TRIM(NumberToVString(EQUATIONS%sparsityType,"*",err,error))//" is invalid."
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
            localError="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%solutionMethod,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%setupType,"*",err,error))// &
            & " is invalid for a bioelectric domain equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%setupType,"*",err,error))// &
          & " is invalid for reaction diffusion equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
    
    EXITS("ReactionDiffusion_EquationsSetSetup")
    RETURN
999 ERRORSEXITS("ReactionDiffusion_EquationsSetSetup",err,error)
    RETURN 1
    
  END SUBROUTINE ReactionDiffusion_EquationsSetSetup

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a reaction diffusion equation type of a classical equations set class.
  SUBROUTINE ReactionDiffusion_EquationsSetSolutionMethodSet(EQUATIONS_SET,SOLUTION_METHOD,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("ReactionDiffusion_EquationsSetSolutionMethodSet",err,error,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a reaction-diffusion type equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
      CASE(EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE, &
         & EQUATIONS_SET_CELLML_REAC_NO_SPLIT_REAC_DIFF_SUBTYPE, &
         & EQUATIONS_SET_CONSTANT_REAC_DIFF_SUBTYPE)        
        SELECT CASE(SOLUTION_METHOD)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          EQUATIONS_SET%solutionMethod=EQUATIONS_SET_FEM_SOLUTION_METHOD
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
          localError="The specified solution method of "//TRIM(NumberToVString(SOLUTION_METHOD,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="Equations set subtype of "//TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
          & " is not valid for a reaction diffusion equation type of classical equations set class."
        CALL FlagError(localError,err,error,*999)
      END SELECT

    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
       
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
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: subtype

    ENTERS("ReactionDiffusion_EquationsSetSpecificationSet",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(SIZE(specification,1)>3) THEN
        CALL FlagError("Equations set specification must have 3 entries for a reaction-diffusion type equations set.", &
          & err,error,*999)
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
      IF(ALLOCATED(equationsSet%specification)) THEN
        CALL FlagError("Equations set specification is already allocated.",err,error,*999)
      ELSE
        ALLOCATE(equationsSet%specification(3),stat=err)
        IF(err/=0) CALL FlagError("Could not allocate equations set specification.",err,error,*999)
      END IF
      equationsSet%specification(1:3)=[EQUATIONS_SET_CLASSICAL_FIELD_CLASS,EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE,subtype]
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    END IF

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
  SUBROUTINE ReactionDiffusion_FiniteElementCalculate(EQUATIONS_SET,ELEMENT_NUMBER,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) FIELD_VAR_TYPE,mh,mhs,ms,ng,nh,nhs,ni,nj,ns,component_idx
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

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a reaction-diffusion type equations set.", &
          & err,error,*999)
      END IF
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        NULLIFY(vectorEquations)
        CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
        dependentField=>equations%interpolation%dependentField
        geometricField=>equations%interpolation%geometricField
        materialsField=>equations%interpolation%materialsField
        IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CELLML_REAC_NO_SPLIT_REAC_DIFF_SUBTYPE .OR. &
            & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_REAC_DIFF_SUBTYPE) THEN
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
        IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CELLML_REAC_NO_SPLIT_REAC_DIFF_SUBTYPE .OR. &
            & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_REAC_DIFF_SUBTYPE) THEN
           CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
             & sourceInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
        ENDIF                  
        dynamicMatrices=>vectorMatrices%dynamicMatrices
        stiffnessMatrix=>dynamicMatrices%matrices(1)%ptr
        dampingMatrix=>dynamicMatrices%matrices(2)%ptr
        rhsVector=>vectorMatrices%rhsVector
        IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CELLML_REAC_NO_SPLIT_REAC_DIFF_SUBTYPE .OR. &
            & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_REAC_DIFF_SUBTYPE) THEN
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
            IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CELLML_REAC_NO_SPLIT_REAC_DIFF_SUBTYPE .OR. &
            & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_REAC_DIFF_SUBTYPE) THEN
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
            component_idx=GEOMETRIC_VARIABLE%numberOfComponents+1
            STORAGE_COEFFICIENT=equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(component_idx,1)
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
            IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_REAC_DIFF_SUBTYPE) THEN
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
              IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_REAC_DIFF_SUBTYPE) THEN
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
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: problemSubtype

    ENTERS("ReactionDiffusion_ProblemSpecificationSet",err,error,*999)

    IF(ASSOCIATED(problem)) THEN
      IF(SIZE(problemSpecification,1)>=3) THEN
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
        IF(ALLOCATED(problem%specification)) THEN
          CALL FlagError("Problem specification is already allocated.",err,error,*999)
        ELSE
          ALLOCATE(problem%specification(3),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
        END IF
        problem%specification(1:3)=[PROBLEM_CLASSICAL_FIELD_CLASS,PROBLEM_REACTION_DIFFUSION_EQUATION_TYPE,problemSubtype]
      ELSE
        CALL FlagError("Reaction-diffusion problem specification must have >=3 entries.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    END IF

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
  SUBROUTINE REACTION_DIFFUSION_EQUATION_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: PROBLEM !<A pointer to the problem set to setup a bioelectric domain equation on.
    TYPE(ProblemSetupType), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(CELLML_EQUATIONS_TYPE), POINTER :: CELLML_EQUATIONS
    TYPE(ControlLoopType), POINTER :: CONTROL_LOOP,CONTROL_LOOP_ROOT
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: localError

    NULLIFY(CELLML_EQUATIONS)
    NULLIFY(CONTROL_LOOP)
    NULLIFY(SOLVER)
    NULLIFY(SOLVERS)
    NULLIFY(SOLVER_EQUATIONS)
    
    ENTERS("REACTION_DIFFUSION_EQUATION_PROBLEM_SETUP",err,error,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(.NOT.ALLOCATED(PROBLEM%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(PROBLEM%SPECIFICATION,1)<3) THEN
        CALL FlagError("Problem specification must have three entries for a reaction diffusion problem.",err,error,*999)
      END IF
      SELECT CASE(PROBLEM_SETUP%setupType)
      CASE(PROBLEM_SETUP_INITIAL_TYPE)
        SELECT CASE(PROBLEM_SETUP%actionType)
        CASE(PROBLEM_SETUP_START_ACTION)
          !Do nothing????
        CASE(PROBLEM_SETUP_FINISH_ACTION)
          !Do nothing????
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%setupType,"*",err,error))// &
            & " is invalid for a reaction diffusion equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PROBLEM_SETUP_CONTROL_TYPE)
        SELECT CASE(PROBLEM_SETUP%actionType)
        CASE(PROBLEM_SETUP_START_ACTION)
          !Set up a time control loop
          CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,err,error,*999)
          CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,CONTROL_TIME_LOOP_TYPE,err,error,*999)
          CALL CONTROL_LOOP_LABEL_SET(CONTROL_LOOP,"Time Loop",err,error,*999)
        CASE(PROBLEM_SETUP_FINISH_ACTION)
          !Finish the control loops
          CONTROL_LOOP_ROOT=>PROBLEM%controlLoop
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
          CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,err,error,*999)            
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%setupType,"*",err,error))// &
            & " is invalid for a reaction-diffusion equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PROBLEM_SETUP_SOLVERS_TYPE)
        !Get the control loop
        CONTROL_LOOP_ROOT=>PROBLEM%controlLoop
        CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
        SELECT CASE(PROBLEM_SETUP%actionType)
        CASE(PROBLEM_SETUP_START_ACTION)
          !Start the solvers creation
          CALL SOLVERS_CREATE_START(CONTROL_LOOP,SOLVERS,err,error,*999)
          SELECT CASE(PROBLEM%SPECIFICATION(3))
          CASE(PROBLEM_CELLML_REAC_INTEG_REAC_DIFF_STRANG_SPLIT_SUBTYPE)
            CALL SOLVERS_NUMBER_SET(SOLVERS,3,err,error,*999)
            !Set the first solver to be a differential-algebraic equations solver
            NULLIFY(SOLVER)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DAE_TYPE,err,error,*999)
            CALL SOLVER_LABEL_SET(SOLVER,"First ODE solver",err,error,*999)
            !Set solver defaults
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,err,error,*999)
            !Set the second solver to be a dynamic solver 
            NULLIFY(SOLVER)
            CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,err,error,*999)
            CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
            CALL SOLVER_LABEL_SET(SOLVER,"Parabolic solver",err,error,*999)
            !Set solver defaults
            CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
            CALL SOLVER_DYNAMIC_SCHEME_SET(SOLVER,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,err,error,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,err,error,*999)
            !Set the third solver to be a differential-algebraic equations solver
            NULLIFY(SOLVER)
            CALL SOLVERS_SOLVER_GET(SOLVERS,3,SOLVER,err,error,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DAE_TYPE,err,error,*999)
            CALL SOLVER_LABEL_SET(SOLVER,"Second ODE solver",err,error,*999)
            !Set solver defaults
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,err,error,*999)
          CASE(PROBLEM_CELLML_REAC_EVAL_REAC_DIFF_NO_SPLIT_SUBTYPE)
            CALL SOLVERS_NUMBER_SET(SOLVERS,2,err,error,*999)
            !Set the first solver to be a CELLML evaluator equations solver
            NULLIFY(SOLVER)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_CELLML_EVALUATOR_TYPE,err,error,*999)
            CALL SOLVER_LABEL_SET(SOLVER,"Evaluator solver",err,error,*999)
            !Set solver defaults
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,err,error,*999)
            !Set the second solver to be a dynamic solver 
            NULLIFY(SOLVER)
            CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,err,error,*999)
            CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
            CALL SOLVER_LABEL_SET(SOLVER,"Parabolic solver",err,error,*999)
            !Set solver defaults
            CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
            CALL SOLVER_DYNAMIC_SCHEME_SET(SOLVER,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,err,error,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,err,error,*999)
          CASE(PROBLEM_CONSTANT_REAC_DIFF_NO_SPLIT_SUBTYPE)
            CALL SOLVERS_NUMBER_SET(SOLVERS,1,err,error,*999)
            !Set the solver to be a dynamic solver 
            NULLIFY(SOLVER)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,err,error,*999)
            CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
            CALL SOLVER_LABEL_SET(SOLVER,"Parabolic solver",err,error,*999)
            !Set solver defaults
            CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
            CALL SOLVER_DYNAMIC_SCHEME_SET(SOLVER,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,err,error,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,err,error,*999)
          CASE DEFAULT
            localError="The problem subtype of "//TRIM(NumberToVString(PROBLEM%SPECIFICATION(3),"*",err,error))// &
              & " is invalid for a reaction-diffusion problem type of a classical problem class."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_FINISH_ACTION)
          !Get the solvers
          CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
          !Finish the solvers creation
          CALL SOLVERS_CREATE_FINISH(SOLVERS,err,error,*999)
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%setupType,"*",err,error))// &
            & " is invalid for a classical equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
        SELECT CASE(PROBLEM_SETUP%actionType)
        CASE(PROBLEM_SETUP_START_ACTION)
          !Get the control loop
          CONTROL_LOOP_ROOT=>PROBLEM%controlLoop
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
          SELECT CASE(PROBLEM%SPECIFICATION(3))
          CASE(PROBLEM_CELLML_REAC_INTEG_REAC_DIFF_STRANG_SPLIT_SUBTYPE)
            !Get the solver
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            !Create the solver equations for the second (parabolic) solver
            CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,err,error,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,err,error,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
          CASE(PROBLEM_CELLML_REAC_EVAL_REAC_DIFF_NO_SPLIT_SUBTYPE)
            !Get the solver
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            !Create the solver equations for the parabolic solver
            CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,err,error,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
          CASE(PROBLEM_CONSTANT_REAC_DIFF_NO_SPLIT_SUBTYPE)
            !Get the solver
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            !Create the solver equations for the parabolic solver
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,err,error,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,err,error,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%actionType,"*",err,error))// &
              & " for a setup subtype of "//TRIM(NumberToVString(PROBLEM%SPECIFICATION(3),"*",err,error))// &
              & " is invalid for a reaction-diffusion equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT

        CASE(PROBLEM_SETUP_FINISH_ACTION)
          !Get the control loop
          CONTROL_LOOP_ROOT=>PROBLEM%controlLoop
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
          CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
          SELECT CASE(PROBLEM%SPECIFICATION(3))
          CASE(PROBLEM_CELLML_REAC_INTEG_REAC_DIFF_STRANG_SPLIT_SUBTYPE)
            !Get the solver equations for the second (parabolic) solver
            CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,err,error,*999)
            !Finish the solver equations creation
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,err,error,*999)             
          CASE(PROBLEM_CELLML_REAC_EVAL_REAC_DIFF_NO_SPLIT_SUBTYPE)
            !Get the solver equations for the parabolic solver
            CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,err,error,*999)
            !Finish the solver equations creation
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,err,error,*999)             
          CASE(PROBLEM_CONSTANT_REAC_DIFF_NO_SPLIT_SUBTYPE)
            !Get the solver equations for thE solver
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,err,error,*999)
            !Finish the solver equations creation
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,err,error,*999)             
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%actionType,"*",err,error))// &
              & " for a setup subtype of "//TRIM(NumberToVString(PROBLEM%SPECIFICATION(3),"*",err,error))// &
              & " is invalid for a reaction-diffusion equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%setupType,"*",err,error))// &
            & " is invalid for a reaction-diffusion equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PROBLEM_SETUP_CELLML_EQUATIONS_TYPE)
        SELECT CASE(PROBLEM_SETUP%actionType)
        CASE(PROBLEM_SETUP_START_ACTION)
          !Get the control loop
          CONTROL_LOOP_ROOT=>PROBLEM%controlLoop
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
          !Get the solver
          CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)

          IF(PROBLEM%SPECIFICATION(3)==PROBLEM_CELLML_REAC_INTEG_REAC_DIFF_STRANG_SPLIT_SUBTYPE) THEN
            NULLIFY(SOLVER)
            NULLIFY(CELLML_EQUATIONS)
            !Create the CellML equations for the first DAE solver
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
            CALL CELLML_EQUATIONS_CREATE_START(SOLVER,CELLML_EQUATIONS,err,error,*999)
            !Set the time dependence
            CALL CellMLEquations_TimeDependenceTypeSet(CELLML_EQUATIONS,CELLML_EQUATIONS_DYNAMIC,err,error,*999)
            !Set the linearity
            CALL CellMLEquations_LinearityTypeSet(CELLML_EQUATIONS,CELLML_EQUATIONS_NONLINEAR,err,error,*999)
            !Create the CellML equations for the second DAE solver
            NULLIFY(SOLVER)
            NULLIFY(CELLML_EQUATIONS)
            CALL SOLVERS_SOLVER_GET(SOLVERS,3,SOLVER,err,error,*999)
            CALL CELLML_EQUATIONS_CREATE_START(SOLVER,CELLML_EQUATIONS,err,error,*999)
            !Set the time dependence
            CALL CellMLEquations_TimeDependenceTypeSet(CELLML_EQUATIONS,CELLML_EQUATIONS_DYNAMIC,err,error,*999)
            !Set the linearity
            CALL CellMLEquations_LinearityTypeSet(CELLML_EQUATIONS,CELLML_EQUATIONS_NONLINEAR,err,error,*999)
          ELSEIF(PROBLEM%SPECIFICATION(3)== &
           & PROBLEM_CELLML_REAC_EVAL_REAC_DIFF_NO_SPLIT_SUBTYPE) THEN
            !CREATE the CellML equations for the first evaluator solver
            NULLIFY(SOLVER)
            NULLIFY(CELLML_EQUATIONS)
            !Create the CellML equations for the first cellml evaluator solver
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
            CALL CELLML_EQUATIONS_CREATE_START(SOLVER,CELLML_EQUATIONS,err,error,*999)
            !Set the time dependence
            CALL CellMLEquations_TimeDependenceTypeSet(CELLML_EQUATIONS,CELLML_EQUATIONS_DYNAMIC,err,error,*999)
            !Set the linearity
            CALL CellMLEquations_LinearityTypeSet(CELLML_EQUATIONS,CELLML_EQUATIONS_NONLINEAR,err,error,*999)
          ENDIF
        CASE(PROBLEM_SETUP_FINISH_ACTION)
          !Get the control loop
          CONTROL_LOOP_ROOT=>PROBLEM%controlLoop
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
          CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
          SELECT CASE(PROBLEM%SPECIFICATION(3))
          CASE(PROBLEM_CELLML_REAC_INTEG_REAC_DIFF_STRANG_SPLIT_SUBTYPE)
            !Get the CellML equations for the first DAE solver
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
            CALL SOLVER_CELLML_EQUATIONS_GET(SOLVER,CELLML_EQUATIONS,err,error,*999)
            !Finish the CellML equations creation
            CALL CELLML_EQUATIONS_CREATE_FINISH(CELLML_EQUATIONS,err,error,*999)
            !Get the CellML equations for the second DAE solver
            NULLIFY(SOLVER)
            NULLIFY(CELLML_EQUATIONS)
            CALL SOLVERS_SOLVER_GET(SOLVERS,3,SOLVER,err,error,*999)
            CALL SOLVER_CELLML_EQUATIONS_GET(SOLVER,CELLML_EQUATIONS,err,error,*999)
            !Finish the CellML equations creation
            CALL CELLML_EQUATIONS_CREATE_FINISH(CELLML_EQUATIONS,err,error,*999)
          CASE(PROBLEM_CELLML_REAC_EVAL_REAC_DIFF_NO_SPLIT_SUBTYPE)
            !Get the CellML equations for the first evaluator solver
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
            CALL SOLVER_CELLML_EQUATIONS_GET(SOLVER,CELLML_EQUATIONS,err,error,*999)
            !Finish the CellML equations creation
            CALL CELLML_EQUATIONS_CREATE_FINISH(CELLML_EQUATIONS,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%actionType,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM%SPECIFICATION(3),"*",err,error))// &
              & " is invalid for reaction-diffusion equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%setupType,"*",err,error))// &
            & " is invalid for reaction-diffusion equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%setupType,"*",err,error))// &
          & " is invalid for areaction-diffusion equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    ENDIF
      
    EXITS("REACTION_DIFFUSION_EQUATION_PROBLEM_SETUP")
    RETURN
999 ERRORSEXITS("REACTION_DIFFUSION_EQUATION_PROBLEM_SETUP",err,error)
    RETURN 1
  END SUBROUTINE REACTION_DIFFUSION_EQUATION_PROBLEM_SETUP
  
  !
  !================================================================================================================================
  !
  
  !>Performs pre-solve actions for reaction-diffusion problems.
  SUBROUTINE ReactionDiffusion_PreSolve(solver,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer to the solver to perform the pre-solve actions for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    REAL(DP) :: currentTime,timeIncrement
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemType), POINTER :: problem
     TYPE(VARYING_STRING) :: localError

    ENTERS("ReactionDiffusion_PreSolve",err,error,*999)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    IF(.NOT.ALLOCATED(problem%specification)) CALL FlagError("Problem specification is not allocated.",err,error,*999)
    IF(SIZE(problem%specification,1)<3) &
      & CALL FlagError("Problem specification must have at least three entries for an advection-diffusion problem.",err,error,*999)
    
    SELECT CASE(problem%specification(3))
    CASE(PROBLEM_CELLML_REAC_INTEG_REAC_DIFF_STRANG_SPLIT_SUBTYPE)
      CALL ControlLoop_CurrentTimesGet(controlLoop,currentTime,timeIncrement,err,error,*999)    
      SELECT CASE(solver%globalNumber)
      CASE(1)
        CALL Solver_DAETimesSet(solver,currentTime,currentTime+timeIncrement/2.0_DP,err,error,*999)
      CASE(2)
        !Do nothing
      CASE(3)
        CALL Solver_DAETimesSet(solver,currentTime,currentTime+timeIncrement/2.0_DP,err,error,*999)
      CASE DEFAULT
        localError="The solver global number of "//TRIM(NumberToVString(SOLVER%globalNumber,"*",err,error))// &
          & " is invalid for a Strang split reaction-diffusion problem."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_CONSTANT_REAC_DIFF_NO_SPLIT_SUBTYPE)
      !No splitting, therefore entire problem is solved as a dynamic one, with 1 solver nothing to do.
    CASE(PROBLEM_CELLML_REAC_EVAL_REAC_DIFF_NO_SPLIT_SUBTYPE)
      !No splitting, therefore entire problem is solved as a dynamic one, with 1 solver nothing to do.
    CASE DEFAULT
      localError="The problem subtype of "//TRIM(NumberToVString(problem%specification(3),"*",err,error))// &
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
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemType), POINTER :: problem
    TYPE(SOLVERS_TYPE), POINTER :: solvers
    TYPE(SOLVER_TYPE), POINTER :: pdeSolver
    TYPE(VARYING_STRING) :: localError

    ENTERS("ReactionDiffusion_PostSolve",err,error,*999)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    IF(.NOT.ALLOCATED(problem%specification)) CALL FlagError("Problem specification is not allocated.",err,error,*999)
    IF(SIZE(problem%specification,1)<3) &
      & CALL FlagError("Problem specification must have at least three entries for an reaction-diffusion problem.",err,error,*999)
    
    SELECT CASE(problem%specification(3))
      
    CASE(PROBLEM_CELLML_REAC_INTEG_REAC_DIFF_STRANG_SPLIT_SUBTYPE)
      SELECT CASE(solver%globalNumber)
      CASE(1)        
        !do nothing
      CASE(2)        
        !do nothing
      CASE(3)
        !OUTPUT SOLUTIONS AT EACH TIME STEP - should probably change this bit below to output 
        !mesh solutions directly from the 3rd solver itself rather than by getting the 2nd solver.
        !I just don't know how to work with cellml_equations to do this.
        solvers=>solver%solvers
        NULLIFY(pdeSolver)
        CALL SOLVERS_SOLVER_GET(solvers,2,pdeSolver,err,error,*999)
        CALL REACTION_DIFFUSION_POST_SOLVE_OUTPUT_DATA(controlLoop,pdeSolver,err,error,*999)
      CASE DEFAULT
        localError="The solver global number of "//TRIM(NumberToVString(solver%globalNumber,"*",err,error))// &
          & " is invalid for a Strang split reaction-diffusion problem."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE (PROBLEM_CELLML_REAC_EVAL_REAC_DIFF_NO_SPLIT_SUBTYPE)
      !do nothing - time output not implemented
    CASE (PROBLEM_CONSTANT_REAC_DIFF_NO_SPLIT_SUBTYPE)
      !OUTPUT SOLUTIONS AT TIME STEP
      CALL REACTION_DIFFUSION_POST_SOLVE_OUTPUT_DATA(controlLoop,solver,err,error,*999)
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(problem%specification(3),"*",err,error))// &
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

  SUBROUTINE REACTION_DIFFUSION_POST_SOLVE_OUTPUT_DATA(CONTROL_LOOP,SOLVER,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER!<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    TYPE(EquationsSetType), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    TYPE(VARYING_STRING) :: localError
    TYPE(WorkGroupType), POINTER :: workGroup

    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT
    INTEGER(INTG) :: EQUATIONS_SET_IDX,CURRENT_LOOP_ITERATION,OUTPUT_FREQUENCY,MAX_DIGITS
    INTEGER(INTG) :: myGroupComputationNodeNumber
    LOGICAL :: exportExelem

    CHARACTER(30) :: FILE
    CHARACTER(30) :: OUTPUT_FILE
    CHARACTER(100) :: FMT, TEMP_FMT

    ENTERS("REACTION_DIFFUSION_POST_SOLVE_OUTPUT_DATA",err,error,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a reaction diffusion problem.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_CELLML_REAC_INTEG_REAC_DIFF_STRANG_SPLIT_SUBTYPE, &
              & PROBLEM_CONSTANT_REAC_DIFF_NO_SPLIT_SUBTYPE)
              CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CURRENT_TIME,TIME_INCREMENT,err,error,*999)
              SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
              IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                SOLVER_MAPPING=>SOLVER_EQUATIONS%solverMapping
                IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                  !Make sure the equations sets are up to date
                  DO equations_set_idx=1,SOLVER_MAPPING%numberOfEquationsSets
                    EQUATIONS_SET=>SOLVER_MAPPING%equationsSets(equations_set_idx)%PTR
                    CURRENT_LOOP_ITERATION=CONTROL_LOOP%timeLoop%iterationNumber
                    OUTPUT_FREQUENCY=CONTROL_LOOP%timeLoop%outputNumber
                    NULLIFY(workGroup)
                    CALL Solver_WorkGroupGet(solver,workGroup,err,error,*999)
                    CALL WorkGroup_GroupNodeNumberGet(workGroup,myGroupComputationNodeNumber,err,error,*999)
                    MAX_DIGITS=FLOOR(LOG10((CONTROL_LOOP%timeLoop%stopTime-CONTROL_LOOP%timeLoop%startTime)/ &
                      & CONTROL_LOOP%timeLoop%timeIncrement))+1
                    IF(OUTPUT_FREQUENCY>0) THEN
                      IF(MOD(CURRENT_LOOP_ITERATION,OUTPUT_FREQUENCY)==0) THEN
                        IF(CONTROL_LOOP%timeLoop%currentTime<=CONTROL_LOOP%timeLoop%stopTime) THEN
                          IF(SOLVER_MAPPING%numberOfEquationsSets.EQ.1) THEN
                            WRITE(TEMP_FMT,'("I",I0,".",I0)') MAX_DIGITS,MAX_DIGITS
                            !100 FORMAT 
                            FMT = TRIM(TEMP_FMT)
                            WRITE(TEMP_FMT,'(A2,A38,A20,A2)') "(", '"TIME_STEP_SPEC_1.part",I2.2,".",',FMT,")"
                            FMT = TRIM(TEMP_FMT)
                            WRITE(OUTPUT_FILE,FMT) &
                              & myGroupComputationNodeNumber,CURRENT_LOOP_ITERATION
                          ELSE
                            WRITE(TEMP_FMT,'("I",I0,".",I0)') MAX_DIGITS,MAX_DIGITS
                            !200 FORMAT 
                            FMT = TRIM(TEMP_FMT)
                            WRITE(TEMP_FMT,'(A2,A38,A20,A2)') "(", '"TIME_STEP_SPEC_",I0,".part",I2.2,".",',FMT,")"
                            FMT = TRIM(TEMP_FMT)
                            WRITE(OUTPUT_FILE,FMT) &
                              & equations_set_idx, myGroupComputationNodeNumber,CURRENT_LOOP_ITERATION
                          ENDIF
                          WRITE(*,*) OUTPUT_FILE
                          FILE=TRIM(OUTPUT_FILE)
                          CALL WriteString(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
                          CALL WriteString(GENERAL_OUTPUT_TYPE,"Now export fields... ",err,error,*999)
                          exportExelem = .False.
                          IF (CURRENT_LOOP_ITERATION==0) THEN
                             IF (equations_set_idx==1) exportExelem = .True.
                          ENDIF
                          CALL REACTION_DIFFUSION_IO_WRITE_CMGUI(EQUATIONS_SET%REGION,EQUATIONS_SET%globalNumber,FILE, &
                            & exportExelem,err,error,*999)
                        ENDIF
                      ENDIF 
                    ENDIF
                  ENDDO
                ENDIF
              ENDIF
            CASE(PROBLEM_CELLML_REAC_EVAL_REAC_DIFF_NO_SPLIT_SUBTYPE)
              ! do nothing ???
              CALL FlagError("Not implemented.",err,error,*999)
            CASE DEFAULT
              localError="Problem subtype "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
                & " is not valid for an advection-diffusion equation type of a classical field problem class."
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

    EXITS("REACTION_DIFFUSION_POST_SOLVE_OUTPUT_DATA")
    RETURN
999 ERRORSEXITS("REACTION_DIFFUSION_POST_SOLVE_OUTPUT_DATA",err,error)
    RETURN 1
    
  END SUBROUTINE REACTION_DIFFUSION_POST_SOLVE_OUTPUT_DATA
  !
  !================================================================================================================================
  !

  SUBROUTINE REACTION_DIFFUSION_CONTROL_LOOP_POST_LOOP(CONTROL_LOOP,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(ProblemType), POINTER :: PROBLEM
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    TYPE(EquationsType), POINTER :: equations !<A pointer to the equations
    TYPE(EquationsVectorType), POINTER :: vectorEquations 
    TYPE(EquationsSetType), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS !<A pointer to the solvers
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatrixType), POINTER :: dampingMatrix,stiffnessMatrix
    
    ENTERS("REACTION_DIFFUSION_CONTROL_LOOP_POST_LOOP",err,error,*999)
    
    NULLIFY(SOLVER)
    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      PROBLEM=>CONTROL_LOOP%PROBLEM
      IF(ASSOCIATED(PROBLEM)) THEN
        SELECT CASE(PROBLEM%SPECIFICATION(3))
        CASE(PROBLEM_CELLML_REAC_INTEG_REAC_DIFF_STRANG_SPLIT_SUBTYPE)
          SOLVERS=>CONTROL_LOOP%SOLVERS
          IF(ASSOCIATED(SOLVERS)) THEN
            CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
            SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
            IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
              SOLVER_MAPPING=>SOLVER_EQUATIONS%solverMapping
              IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                EQUATIONS_SET=>SOLVER_MAPPING%equationsSets(1)%ptr
                IF(ASSOCIATED(EQUATIONS_SET)) THEN
                  EQUATIONS=>EQUATIONS_SET%EQUATIONS
                  IF(ASSOCIATED(EQUATIONS)) THEN
                    NULLIFY(vectorEquations)
                    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
                    vectorMatrices=>vectorEquations%vectorMatrices
                    dynamicMatrices=>vectorMatrices%dynamicMatrices
                    stiffnessMatrix=>dynamicMatrices%matrices(1)%ptr
                    dampingMatrix=>dynamicMatrices%matrices(2)%ptr
                    stiffnessMatrix%updateMatrix = .FALSE.
                    dampingMatrix%updateMatrix = .FALSE.
                  ELSE
                    CALL FlagError("Equations not associated.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Equations Set not associated.",err,error,*999)
                ENDIF
      
              ELSE
                CALL FlagError("Solver Mapping not associated.",err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Solver Equations not associated.", err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Solvers is not associated.", err,error,*999)
          ENDIF


        CASE DEFAULT
          !do nothing
        END SELECT
      ELSE
        CALL FlagError("Problem is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Control Loop is not associated.",err,error,*999)
    ENDIF
    EXITS("REACTION_DIFFUSION_CONTROL_LOOP_POST_LOOP")
    RETURN
999 ERRORSEXITS("REACTION_DIFFUSION_CONTROL_LOOP_POST_LOOP",err,error)
    RETURN 1
  END SUBROUTINE REACTION_DIFFUSION_CONTROL_LOOP_POST_LOOP
END MODULE REACTION_DIFFUSION_EQUATION_ROUTINES
