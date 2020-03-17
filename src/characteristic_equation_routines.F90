!> \file  
!> \author David Ladd
!> \brief This module handles the characteristic equation routines. These 
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
!> Contributor(s): David Ladd, Soroush Safaei, Chris Bradley
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

!>This module handles all characteristic equation routines.
MODULE CharacteristicEquationsRoutines

  USE BaseRoutines
  USE BasisRoutines
  USE BoundaryConditionsRoutines
  USE Constants
  USE ControlLoopRoutines
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

  PUBLIC Characteristic_EquationsSetSolutionMethodSet
  
  PUBLIC Characteristic_EquationsSetSpecificationSet
  
  PUBLIC Characteristic_EquationsSetSetup
  
  PUBLIC Characteristic_NodalResidualEvaluate
  
  PUBLIC Characteristic_NodalJacobianEvaluate
  
  PUBLIC Characteristic_Extrapolate
  
  PUBLIC Characteristic_PrimitiveToCharacteristic

CONTAINS 

!
!================================================================================================================================
!

  !>Sets/changes the solution method for a Characteristic equation type of an fluid mechanics equations set class.
  SUBROUTINE Characteristic_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet
    INTEGER(INTG), INTENT(IN) :: solutionMethod
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    INTEGER(INTG) :: esSpecification(3)
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Characteristic_EquationsSetSolutionMethodSet",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
    
    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_CHARACTERISTIC_SUBTYPE)                                
      SELECT CASE(solutionMethod)
      CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(EQUATIONS_SET_NODAL_SOLUTION_METHOD)
        equationsSet%solutionMethod=EQUATIONS_SET_NODAL_SOLUTION_METHOD
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
        & " is not valid for a characteristic type of a fluid mechanics equations set."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("Characteristic_EquationsSetSolutionMethodSet")
    RETURN
999 ERRORS("Characteristic_EquationsSetSolutionMethodSet",err,error)
    EXITS("Characteristic_EquationsSetSolutionMethodSet")
    RETURN 1
    
  END SUBROUTINE Characteristic_EquationsSetSolutionMethodSet

!
!================================================================================================================================
!

  !>Sets the equation specification for a Characteristic type of a fluid mechanics equations set class.
  SUBROUTINE Characteristic_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the specification for
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations set specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: subtype

    ENTERS("Characteristic_EquationsSetSpecificationSet",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(ALLOCATED(equationsSet%specification)) CALL FlagError("Equations set specification is already allocated.",err,error,*999)
    IF(SIZE(specification,1)<3) THEN
      localError="The size of the specified specification array of "// &
        & TRIM(NumberToVString(SIZE(specification,1),"*",err,error))//" is invalid. The size should be >= 3."
      CALL FlagError(localError,err,error,*999)
    END IF
    
    subtype=specification(3)
    SELECT CASE(subtype)
    CASE(EQUATIONS_SET_CHARACTERISTIC_SUBTYPE)
      !ok
    CASE DEFAULT
      localError="The third equations set specification of "//TRIM(NumberToVstring(subtype,"*",err,error))// &
        & " is not valid for a characteristic type of a fluid mechanics equations set."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    !Set full specification
    ALLOCATE(equationsSet%specification(3),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate equations set specification.",err,error,*999)
    equationsSet%specification(1:3)=[EQUATIONS_SET_FLUID_MECHANICS_CLASS,EQUATIONS_SET_CHARACTERISTIC_EQUATION_TYPE,subtype]

    EXITS("Characteristic_EquationsSetSpecificationSet")
    RETURN
999 ERRORS("Characteristic_EquationsSetSpecificationSet",err,error)
    EXITS("Characteristic_EquationsSetSpecificationSet")
    RETURN 1
    
  END SUBROUTINE Characteristic_EquationsSetSpecificationSet

!
!================================================================================================================================
!

  !>Sets up the Characteristic equations fluid setup.
  SUBROUTINE Characteristic_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(EquationsSetSetupType), INTENT(INOUT) :: equationsSetSetup
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    INTEGER(INTG) :: componentIdx,geometricScalingType,geometricMeshComponent,geometricComponentNumber
    INTEGER(INTG) :: dependentFieldNumberOfVariables,dependentFieldNumberOfComponents
    INTEGER(INTG) :: independentFieldNumberOfVariables,independentFieldNumberOfComponents
    INTEGER(INTG) :: materialsFieldNumberOfVariables,materialsFieldNumberOfComponents1,materialsFieldNumberOfComponents2
    TYPE(DecompositionType), POINTER :: geometricDecomposition
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsSetMaterialsType), POINTER :: equationsMaterials
    TYPE(EquationsSetEquationsSetFieldType), POINTER :: equationsEquationsSetField
    TYPE(FieldType), POINTER :: equationsSetField
    TYPE(VARYING_STRING) :: localError

    ENTERS("Characteristic_EquationsSetSetup",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)

    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_CHARACTERISTIC_SUBTYPE)
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
      SELECT CASE(equationsSet%specification(3))
      CASE(EQUATIONS_SET_CHARACTERISTIC_SUBTYPE)
        SELECT CASE(equationsSetSetup%actionType)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
          CALL Characteristic_EquationsSetSolutionMethodSet(equationsSet,EQUATIONS_SET_NODAL_SOLUTION_METHOD,err,error,*999)
          NULLIFY(equationsEquationsSetField)
          CALL EquationsSet_EquationsSetFieldGet(equationsSet,equationsEquationsSetField,err,error,*999)
          IF(equationsEquationsSetField%equationsSetFieldAutoCreated) THEN
            !Create the auto created equations set field field for SUPG element metrics
            CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsEquationsSetField%equationsSetFieldField, &
              & err,error,*999)
            NULLIFY(equationsSetField)
            CALL EquationsSet_EquationsSetFieldFieldGet(equationsSet,equationsSetField,err,error,*999)
            CALL Field_LabelSet(equationsSetField,"Equations Set Field",err,error,*999)
            CALL Field_TypeSetAndLock(equationsSetField,FIELD_GENERAL_TYPE,err,error,*999)
            CALL Field_NumberOfVariablesSet(equationsSetField,1,err,error,*999)
            CALL Field_VariableTypesSetAndLock(equationsSetField,[FIELD_U_VARIABLE_TYPE],err,error,*999)
            CALL Field_VariableLabelSet(equationsSetField,FIELD_U_VARIABLE_TYPE,"W2Initialise",err,error,*999)
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
          localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType, &
            & "*",err,error))// " for a setup type of "//TRIM(NumberToVString(equationsSetSetup% &
            & setupType,"*",err,error))// " is not implemented for a characteristic equations set."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
        !-----------------------------------------------------------------
        ! G e o m e t r i c   f i e l d
        !-----------------------------------------------------------------
        SELECT CASE(equationsSetSetup%actionType)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
          NULLIFY(equationsEquationsSetField)
          CALL EquationsSet_EquationsSetFieldGet(equationsSet,equationsEquationsSetField,err,error,*999)
          NULLIFY(equationsSetField)
          CALL EquationsSet_EquationsSetFieldFieldGet(equationsSet,equationsSetField,err,error,*999)
          NULLIFY(geometricField)
          CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
          IF(equationsEquationsSetField%equationsSetFieldAutoCreated) THEN
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
            CALL Field_ScalingTypeSet(equationsSet%equationsSetField%equationsSetFieldField,geometricScalingType,err,error,*999)
          ELSE
            !Do nothing
          ENDIF
        CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
          ! do nothing
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
            & " is invalid for a characteristic equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
        !-----------------------------------------------------------------
        ! D e p e n d e n t   f i e l d
        !-----------------------------------------------------------------
        NULLIFY(geometricField)
        CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
        dependentFieldNumberOfVariables=5
        dependentFieldNumberOfComponents=2
        SELECT CASE(equationsSetSetup%actionType)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
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
            !set number of variables to 5 (U,DELUDELN,V,U1,U2)=>(Q,A;dQ,dA;W;pCellML;Pressure)
            CALL Field_NumberOfVariablesSetAndLock(equationsSet%dependent%dependentField,dependentFieldNumberOfVariables, &
              & err,error,*999)
            CALL Field_VariableTypesSetAndLock(equationsSet%dependent%dependentField,[FIELD_U_VARIABLE_TYPE, &
              & FIELD_DELUDELN_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE,FIELD_U2_VARIABLE_TYPE], &
              & err,error,*999)
            ! set dimension
            CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
              & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_V_VARIABLE_TYPE, &
              & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_U1_VARIABLE_TYPE, &
              & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_U2_VARIABLE_TYPE, &
              & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            ! set data type
            CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_DP_TYPE,err,error,*999)
            CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
              & FIELD_DP_TYPE,err,error,*999)
            CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_V_VARIABLE_TYPE, &
              & FIELD_DP_TYPE,err,error,*999)
            CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_U1_VARIABLE_TYPE, &
              & FIELD_DP_TYPE,err,error,*999)
            CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_U2_VARIABLE_TYPE, &
              & FIELD_DP_TYPE,err,error,*999)
            ! number of components for U,DELUDELN=2 (Q,A)
            CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
              & dependentFieldNumberOfComponents,err,error,*999)
            CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
              & dependentFieldNumberOfComponents,err,error,*999)
            CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_V_VARIABLE_TYPE, &
              & dependentFieldNumberOfComponents,err,error,*999)
            CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_U1_VARIABLE_TYPE, &
              & dependentFieldNumberOfComponents,err,error,*999)
            CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_U2_VARIABLE_TYPE, &
              & dependentFieldNumberOfComponents,err,error,*999)
            CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricMeshComponent,err,error,*999)
            !Default to the geometric interpolation setup for U,dUdN
            DO componentIdx=1,dependentFieldNumberOfComponents
              CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,componentIdx, &
                & geometricMeshComponent,err,error,*999)
              CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
                & componentIdx,geometricMeshComponent,err,error,*999)
              CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_V_VARIABLE_TYPE,componentIdx, &
                & geometricMeshComponent,err,error,*999)
              CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_U1_VARIABLE_TYPE,componentIdx, &
                & geometricMeshComponent,err,error,*999)
              CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_U2_VARIABLE_TYPE,componentIdx, &
                & geometricMeshComponent,err,error,*999)
            ENDDO !componentIdx
            CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
            SELECT CASE(solutionMethod)
            CASE(EQUATIONS_SET_NODAL_SOLUTION_METHOD)
              !Specify nodal solution method
              ! (U, dUdN); 2 components (Q,A)
              DO componentIdx=1,dependentFieldNumberOfComponents
                CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField, &
                  & FIELD_U_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField, &
                  & FIELD_V_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField, &
                  & FIELD_U1_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField, &
                  & FIELD_U2_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              ENDDO
              CALL Field_ScalingTypeGet(equationsSet%GEOMETRY%geometricField,geometricScalingType, &
                & err,error,*999)
              CALL Field_ScalingTypeSet(equationsSet%dependent%dependentField,geometricScalingType, &
                & err,error,*999)
            CASE DEFAULT
              localError="The solution method of " &
                & //TRIM(NumberToVString(equationsSet%solutionMethod,"*",err,error))// " is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ELSE 
            !Check the user specified field
            CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
            CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_DEPENDENT_TYPE,err,error,*999)
            CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,dependentFieldNumberOfVariables,err,error,*999)
            CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE, &
              & FIELD_DELUDELN_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE,FIELD_U2_VARIABLE_TYPE] &
              & ,err,error,*999)
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE, & 
              & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE, &
              & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE, & 
              & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U1_VARIABLE_TYPE, & 
              & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U2_VARIABLE_TYPE, & 
              & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U1_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U2_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            !calculate number of components (Q,A) for U and dUdN
            CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE, &
              & dependentFieldNumberOfComponents,err,error,*999)
            CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE, &
              & dependentFieldNumberOfComponents,err,error,*999)
            CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE, &
              & dependentFieldNumberOfComponents,err,error,*999)
            CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U1_VARIABLE_TYPE, &
              & dependentFieldNumberOfComponents,err,error,*999)
            CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U2_VARIABLE_TYPE, &
              & dependentFieldNumberOfComponents,err,error,*999)
            SELECT CASE(equationsSet%solutionMethod)
            CASE(EQUATIONS_SET_NODAL_SOLUTION_METHOD)
              CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,1, &
                & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,1, &
                & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U1_VARIABLE_TYPE,1, &
                & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U2_VARIABLE_TYPE,1, &
                & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
            CASE DEFAULT
              localError="The solution method of "//TRIM(NumberToVString(equationsSet%solutionMethod, &
                & "*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ENDIF
          !Specify finish action
        CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
          IF(equationsSet%dependent%dependentFieldAutoCreated) THEN
            CALL Field_CreateFinish(equationsSet%dependent%dependentField,err,error,*999)
          ENDIF
        CASE DEFAULT
          localError="The third equations set specification of "// &
            & TRIM(NumberToVString(equationsSet%specification(3),"*",err,error))// &
            & " is invalid for a characteristic equations set."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_SETUP_INDEPENDENT_TYPE)
        !-----------------------------------------------------------------
        ! I n d e p e n d e n t   f i e l d
        !-----------------------------------------------------------------
        NULLIFY(geometricField)
        CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
        independentFieldNumberOfVariables=1 !set number of variables to 1 (W)
        independentFieldNumberOfComponents=2 !normalDirection for wave relative to node for W1,W2
        SELECT CASE(equationsSetSetup%actionType)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
          IF(equationsSet%INDEPENDENT%independentFieldAutoCreated) THEN
            !Create the auto created independent field
            !start field creation with name 'Independent Field'
            CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsSet%INDEPENDENT%independentField, &
              & err,error,*999)
            !start creation of a new field
            CALL Field_TypeSetAndLock(equationsSet%INDEPENDENT%independentField,FIELD_GENERAL_TYPE,err,error,*999)
            !label the field
            CALL Field_LabelSet(equationsSet%INDEPENDENT%independentField,"Independent Field",err,error,*999)
            !define new created field to be independent
            CALL Field_DependentTypeSetAndLock(equationsSet%INDEPENDENT%independentField,FIELD_INDEPENDENT_TYPE,err,error,*999)
            !look for decomposition rule already defined
            CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
            !apply decomposition rule found on new created field
            CALL Field_DecompositionSetAndLock(equationsSet%INDEPENDENT%independentField,geometricDecomposition,err,error,*999)
            !point new field to geometric field
            CALL Field_GeometricFieldSetAndLock(equationsSet%INDEPENDENT%independentField,geometricField,err,error,*999)
            CALL Field_NumberOfVariablesSetAndLock(equationsSet%INDEPENDENT%independentField,independentFieldNumberOfVariables, &
              & err,error,*999)
            CALL Field_VariableTypesSetAndLock(equationsSet%INDEPENDENT%independentField,[FIELD_U_VARIABLE_TYPE],err,error,*999)
            CALL Field_DimensionSetAndLock(equationsSet%INDEPENDENT%independentField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            !characteristic normal direction (normalWave) is +/- 1
            CALL Field_DataTypeSetAndLock(equationsSet%INDEPENDENT%independentField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE, &
              & err,error,*999)
            !calculate number of components with one component for each dimension
            CALL Field_NumberOfComponentsSetAndLock(equationsSet%INDEPENDENT%independentField,FIELD_U_VARIABLE_TYPE, &
              & independentFieldNumberOfComponents,err,error,*999)
            CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricMeshComponent,err,error,*999)
            !Default to the geometric interpolation setup
            DO componentIdx=1,independentFieldNumberOfComponents
              CALL Field_ComponentMeshComponentSet(equationsSet%INDEPENDENT%independentField,FIELD_U_VARIABLE_TYPE, &
                & componentIdx,geometricMeshComponent,err,error,*999)
            ENDDO !componentIdx
            CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)          
            SELECT CASE(solutionMethod)
            CASE(EQUATIONS_SET_NODAL_SOLUTION_METHOD)
              DO componentIdx=1,independentFieldNumberOfComponents
                CALL Field_ComponentInterpolationSetAndLock(equationsSet%INDEPENDENT%independentField, &
                  & FIELD_U_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              ENDDO !componentIdx
              CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
            CASE DEFAULT
              localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))// " is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ELSE
            !Check the user specified field
            CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
            CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
            CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,independentFieldNumberOfVariables,err,error,*999)
            CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE],err,error,*999)
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE, &
              & independentFieldNumberOfComponents,err,error,*999)
          ENDIF
        CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
          IF(equationsSet%INDEPENDENT%independentFieldAutoCreated) THEN
            CALL Field_CreateFinish(equationsSet%INDEPENDENT%independentField,err,error,*999)
          ENDIF
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
            & " is invalid for a standard characteristic equations set"
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The third equations set specification of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
          & " is invalid for a standard characteristic equations set."
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
      materialsFieldNumberOfVariables=2 ! U type-7 constant / V type-3 variable
      materialsFieldNumberOfComponents1=8
      materialsFieldNumberOfComponents2=3
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        IF(equationsMaterials%materialsFieldAutoCreated) THEN
          !Create the auto created materials field
          !start field creation with name 'Materials Field'
          CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsSet%materials%materialsField,err,error,*999)
          CALL Field_TypeSetAndLock(equationsMaterials%materialsField,FIELD_MATERIAL_TYPE,err,error,*999)
          !label the field
          CALL Field_LabelSet(equationsMaterials%materialsField,"Materials Field",err,error,*999)
          CALL Field_DependentTypeSetAndLock(equationsMaterials%materialsField,FIELD_INDEPENDENT_TYPE,err,error,*999)
          CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
          !apply decomposition rule found on new created field
          CALL Field_DecompositionSetAndLock(equationsMaterials%materialsField,geometricDecomposition,err,error,*999)
          !point new field to geometric field
          CALL Field_GeometricFieldSetAndLock(equationsMaterials%materialsField,geometricField,err,error,*999)
          CALL Field_NumberOfVariablesSet(equationsMaterials%materialsField,materialsFieldNumberOfVariables,err,error,*999)
          CALL Field_VariableTypesSetAndLock(equationsMaterials%materialsField,[FIELD_U_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE], &
            & err,error,*999)
          CALL Field_DimensionSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
            & err,error,*999)
          CALL Field_DimensionSetAndLock(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
            & err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_NumberOfComponentsSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
            & materialsFieldNumberOfComponents1,err,error,*999)
          CALL Field_NumberOfComponentsSetAndLock(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE, &
            & materialsFieldNumberOfComponents2,err,error,*999)
          CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricComponentNumber,err,error,*999)
          CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,1, &
            & geometricComponentNumber,err,error,*999)
          CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE,1, &
            & geometricComponentNumber,err,error,*999)
          DO componentIdx=1,materialsFieldNumberOfComponents1 !(mu,rho,alpha,pressureExternal,LengthScale,TimeScale,MassScale)
            CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
              & componentIdx,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
          ENDDO !componentIdx
          DO componentIdx=1,materialsFieldNumberOfComponents2 !(A0,E,H0)
            CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE, &
              & componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
          ENDDO !componentIdx
          !Default the field scaling to that of the geometric field
          CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
          CALL Field_ScalingTypeSet(equationsMaterials%materialsField,geometricScalingType,err,error,*999)
        ELSE
          !Check the user specified field
          CALL Field_TypeCheck(equationsSetSetup%field,FIELD_MATERIAL_TYPE,err,error,*999)
          CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
          CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,materialsFieldNumberOfVariables,err,error,*999)
          CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE],err,error,*999)
          ! U-variable
          CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,materialsFieldNumberOfComponents1, &
            & err,error,*999)
          CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,materialsFieldNumberOfComponents2, &
            & err,error,*999)
        ENDIF
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        IF(equationsMaterials%materialsFieldAutoCreated) THEN
          !Finish creating the materials field
          CALL Field_CreateFinish(equationsMaterials%materialsField,err,error,*999)
          !Should be initialized from example file \TODO: No, should be initialised here.
        ENDIF
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*", & 
          & err,error))//" for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*", & 
          & err,error))//" is invalid for characteristic equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
      !-----------------------------------------------------------------
      ! E q u a t i o n s 
      !-----------------------------------------------------------------
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        CALL EquationsSet_AssertMaterialsIsFinished(equationsSet,err,error,*999)
        NULLIFY(equations)
        CALL Equations_CreateStart(equationsSet,equations,err,error,*999)
        CALL Equations_LinearityTypeSet(equations,EQUATIONS_NONLINEAR,err,error,*999)
        CALL Equations_TimeDependenceTypeSet(equations,EQUATIONS_STATIC,err,error,*999)
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)          
        SELECT CASE(solutionMethod)
        CASE(EQUATIONS_SET_NODAL_SOLUTION_METHOD)
          !Finish the creation of the equations
          NULLIFY(equations)
          CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
          CALL Equations_CreateFinish(equations,err,error,*999)
          NULLIFY(vectorEquations)
          CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
          !Create the equations mapping.
          NULLIFY(vectorMapping)
          CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_U_VARIABLE_TYPE,vectorMapping,err,error,*999)
          CALL EquationsMappingVector_NumberOfResidualsSet(vectorMapping,1,err,error,*999)
          CALL EquationsMappingVector_ResidualNumberOfVariablesSet(vectorMapping,1,1,err,error,*999)
          CALL EquationsMappingVector_ResidualVariableTypesSet(vectorMapping,1,[FIELD_U_VARIABLE_TYPE],err,error,*999)
          CALL EquationsMappingVector_NumberOfLinearMatricesSet(vectorMapping,1,err,error,*999)
          CALL EquationsMappingVector_LinearMatricesVariableTypesSet(vectorMapping,[FIELD_U_VARIABLE_TYPE],err,error,*999)
          CALL EquationsMappingVector_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE,err,error,*999)
          CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
          !Create the equations matrices
          CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
          CALL Equations_SparsityTypeGet(equations,sparsityType,err,error,*999)
          SELECT CASE(sparsityType)
          CASE(EQUATIONS_MATRICES_FULL_MATRICES)
            CALL EquationsMatricesVector_LinearStorageTypeSet(vectorMatrices,[MATRIX_BLOCK_STORAGE_TYPE],err,error,*999)
            CALL EquationsMatricesVector_NonlinearStorageTypeSet(vectorMatrices,MATRIX_BLOCK_STORAGE_TYPE,err,error,*999)
          CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
            CALL EquationsMatricesVector_LinearStorageTypeSet(vectorMatrices,[MATRIX_COMPRESSED_ROW_STORAGE_TYPE],err,error,*999)
            CALL EquationsMatricesVector_LinearStructureTypeSet(vectorMatrices,[EQUATIONS_MATRIX_NODAL_STRUCTURE],err,error,*999)
            CALL EquationsMatricesVector_NonlinearStorageTypeSet(vectorMatrices,MATRIX_COMPRESSED_ROW_STORAGE_TYPE,err,error,*999)
            CALL EquationsMatricesVector_NonlinearStructureTypeSet(vectorMatrices,EQUATIONS_MATRIX_NODAL_STRUCTURE,err,error,*999)
          CASE DEFAULT
            localError="The equations matrices sparsity type of "//TRIM(NumberToVString(sparsityType,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
          CALL EquationsMatrices_VectorCreateFinish(vectorMatrices,err,error,*999)
          !Use the analytic Jacobian calculation
          CALL EquationsMatricesVector_JacobianCalculationTypeSet(vectorMatrices,FIELD_U_VARIABLE_TYPE,1, &
            & EQUATIONS_JACOBIAN_ANALYTIC_CALCULATED,err,error,*999)
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
          & " is invalid for a characteristics equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
        & " is invalid for a characteristics equation set."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("Characteristic_EquationsSetSetup")
    RETURN
999 ERRORSEXITS("Characteristic_EquationsSetSetup",err,error)
    RETURN 1
    
  END SUBROUTINE Characteristic_EquationsSetSetup

  !
  !================================================================================================================================
  !

  !>Evaluates the residual nodal stiffness matrices and RHS for a characteristic equation nodal equations set.
  SUBROUTINE Characteristic_NodalResidualEvaluate(equationsSet,nodeNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet
    INTEGER(INTG), INTENT(IN) :: nodeNumber
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    INTEGER(INTG) :: derivativeIdx,versionIdx,versionIdx2,componentIdx,rowIdx,columnIdx,componentIdx2,numberOfVersions
    REAL(DP) :: qBif(4),aBif(4),a0Param(4),eParam(4),h0Param(4),beta(4),w(2,4),normalWave(2,4),sum,rho
    REAL(DP), POINTER :: dependentParameters(:),independentParameters(:),materialsParameters(:)
    LOGICAL :: boundaryNode,update,updateMatrix,updateResidual
    TYPE(DomainType), POINTER :: dependentDomain
    TYPE(DomainNodesType), POINTER :: dependentDomainNodes
    TYPE(DomainTopologyType), POINTER :: dependentDomainTopology
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector
    TYPE(EquationsMatrixType), POINTER :: stiffnessMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: materialsField,dependentField,independentField
    TYPE(FieldVariableType), POINTER :: dependentUVariable,dependentVVariable,independentVariable,materialsUVariable, &
      & materialsVVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("Characteristic_NodalResidualEvaluate",err,error,*999)
  
    CALL EquationsSet_SpecificationGet(esSpecification,3,esSpecification,err,error,*999)

    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_CHARACTERISTIC_SUBTYPE)
      !OK
    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is not valid for a Poisson equation type of a classical field equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    NULLIFY(equations)
    CALL EquationSet_EquationsGet(equationsSet,equations,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMapping)
    CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
    NULLIFY(nonlinearMapping)
    CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)
    NULLIFY(linearMapping)
    CALL EquationsMappingVector_LinearMappingGet(vectorMapping,linearMapping,err,error,*999)
    NULLIFY(residualMapping)
    CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,1,residualMapping,err,error,*999)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
    NULLIFY(nonlinearMatrices)
    CALL EquationsMatricesVector_NonlinearMatricesGet(vectorMatrices,nonlinearMatrices,err,error,*999)
    NULLIFY(residualVector)
    CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,1,residualVector,err,error,*999)
    updateResidual=residualVector%updateVector
    NULLIFY(linearMatrices)
    CALL EquationsMatricesVector_LinearMatricesGet(vectorMatrices,linearMatrices,err,error,*999)
    NULLIFY(stiffnessMatrix)
    CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,1,stiffnessMatrix,err,error,*999)
    updateMatrix=stiffnessMatrix%updateMatrix

    update=(updateMatrix.OR.updateResidual)

    IF(update) THEN

      NULLIFY(dependentField)
      CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
      NULLIFY(materialsField)
      CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)
      NULLIFY(independentField)
      CALL EquationsSet_IndependentFieldGet(equationsSet,independentField,err,error,*999)

      NULLIFY(dependentUVariable)
      CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,dependentUVariable,err,error,*999)
      NULLIFY(dependentVVariable)
      CALL Field_VariableGet(dependentField,FIELD_V_VARIABLE_TYPE,dependentVVariable,err,error,*999)
      NULLIFY(dependentDecomposition)
      CALL Field_DecompositionGet(dependentField,dependentDecomposition,err,error,*999)
      NULLIFY(dependentDecompositionTopology)
      CALL Decomposition_DecompositionTopologyGet(dependentDecomposition,dependentDecompositionTopology,err,error,*999)
      NULLIFY(dependentDecompositionElements)
      CALL DecompositionTopology_DecompositionElementsGet(dependentDecompositionTopology,dependentDecompositionElements, &
        & err,error,*999)
      NULLIFY(dependentDomain)
      CALL Decomposition_DomainGet(dependentDecomposition,0,dependentDomain,err,error,*999)
      NULLIFY(dependentDomainTopology)
      CALL Domain_DomainTopologyGet(dependentDomain,dependentDomainTopology,err,error,*999)
      NULLIFY(dependentDomainNodes)
      CALL DomainTopology_DomainNodesGet(dependentDomainTopology,dependentDomainNodes,err,error,*999)
       
      NULLIFY(materialsUVariable)
      CALL Field_VariableGet(materialsField,FIELD_U_VARIABLE_TYPE,materialsUVariable,err,error,*999)
      NULLIFY(materialsVVariable)
      CALL Field_VariableGet(materialsField,FIELD_V_VARIABLE_TYPE,materialsVVariable,err,error,*999)

      NULLIFY(independentVariable)
      CALL Field_VariableGet(independentField,FIELD_U_VARIABLE_TYPE,independentVariable,err,error,*999)

      derivativeIdx=1
      normalWave=0.0_DP
      CALL DomainNodes_DerivativeNumberOfVersionsGet(dependentDomainNodes,derivativeIdx,nodeNumber,numberOfVersions,err,error,*999)
      CALL DomainNodes_NodeBoundaryNodeGet(dependentDomainNodes,nodeNumber,boundaryNode,err,error,*999)
 
      !Get normal wave direction for nodes
      DO componentIdx=1,2
        DO versionIdx=1,numberOfVersions
          CALL FieldVariable_ParameterSetGetLocalNode(independentVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
            & nodeNumber,componentIdx,normalWave(componentIdx,versionIdx),err,error,*999)  
        ENDDO !versionIdx
      ENDDO !componentIdx

      !!!-- F i n d   B r a n c h   N o d e s --!!!
      IF(ABS(normalWave(1,1))>0 .OR. ABS(normalWave(2,1))>0) THEN
        IF(.NOT. boundaryNode) THEN

          !Get material constants
          CALL FieldVariable_ParameterSetGetConstant(materialsUVariable,FIELD_VALUES_SET_TYPE,2,rho,err,error,*999)
          !Get node-based material parameters
          DO versionIdx=1,numberOfVersions
            CALL FieldVariable_ParameterSetGetLocalNode(materialsVVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
              & nodeNumber,1,a0Param(versionIdx),err,error,*999)  
            CALL FieldVariable_ParameterSetGetLocalNode(materialsVVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
              & nodeNumber,2,eParam(versionIdx),err,error,*999)                
            CALL FieldVariable_ParameterSetGetLocalNode(materialsVVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
              & nodeNumber,3,h0Param(versionIdx),err,error,*999)                
            beta(versionIdx)=(4.0_DP*SQRT(PI)*eParam(versionIdx)*h0Param(versionIdx))/(3.0_DP*a0Param(versionIdx))     
          ENDDO !versionIdx

          DO versionIdx=1,numberOfVersions
            !Get current Q & A Values at the node
            CALL FieldVariable_ParameterSetGetLocalNode(dependentUVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
              & nodeNumber,1,qBif(versionIdx),err,error,*999)                
            CALL FieldVariable_ParameterSetGetLocalNode(dependentUVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
              & nodeNumber,2,aBif(versionIdx),err,error,*999)                
            !Set as upwind field values
            CALL FieldVariable_ParameterSetUpdateLocalNode(dependentUVariable,FIELD_UPWIND_VALUES_SET_TYPE,versionIdx,1, &
              & nodeNumber,1,qBif(versionIdx),err,error,*999)            
            CALL FieldVariable_ParameterSetUpdateLocalNode(dependentUVariable,FIELD_UPWIND_VALUES_SET_TYPE,versionIdx,1, &
              & nodeNumber,2,aBif(versionIdx),err,error,*999)
            ! If A goes negative during nonlinear iteration, set to A0
            IF (aBif(versionIdx) < a0Param(versionIdx)*0.001_DP) aBif(versionIdx) = a0Param(versionIdx)*0.001_DP
          ENDDO !versionIdx

          !Get extrapolated W for the node
          DO componentIdx=1,2
            DO versionIdx=1,numberOfVersions
              CALL FieldVariable_ParameterSetGetLocalNode(dependentVVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
                & nodeNumber,componentIdx,w(componentIdx,versionIdx),err,error,*999)                
            ENDDO !versionIdx
          ENDDO !componentIdx

          !!!-- S T I F F N E S S  M A T R I X  --!!!
          IF(updateMatrix) THEN
            !Conservation of Mass
            rowIdx=numberOfVersions+1
            columnIdx=0
            DO componentIdx=1,2
              DO versionIdx=1,numberOfVersions
                IF(ABS(normalWave(componentIdx,versionIdx))>ZERO_TOLERANCE) THEN
                  columnIdx=columnIdx+1
                  stiffnessMatrix%nodalMatrix%matrix(rowIdx,columnIdx)=normalWave(componentIdx,versionIdx)
                ENDIF
              ENDDO !versionIdx
            ENDDO !componentIdx
          ENDIF

          !!!-- N O N L I N E A R   V E C T O R --!!!
          IF(updateResidual) THEN
            rowIdx=0
            !Characteristics Equations
            DO componentIdx=1,2
              DO versionIdx=1,numberOfVersions
                IF(ABS(normalWave(componentIdx,versionIdx))>ZERO_TOLERANCE) THEN
                  rowIdx=rowIdx+1
                  residualVector%nodalResidual%vector(rowIdx)=(qBif(versionIdx)/aBif(versionIdx)) &
                    & +normalWave(componentIdx,versionIdx)*4.0_DP*SQRT(beta(versionIdx)/(2.0_DP*rho))* &
                    & (aBif(versionIdx)**0.25_DP - a0Param(versionIdx)**0.25_DP)-w(componentIdx,versionIdx)
                ENDIF
              ENDDO !versionIdx
            ENDDO !componentIdx
            !Continuity of Total Pressure
            DO componentIdx=1,2
              DO versionIdx=1,numberOfVersions
                IF(ABS(normalWave(componentIdx,versionIdx))>ZERO_TOLERANCE) THEN
                  rowIdx=rowIdx+1
                  IF(versionIdx==1) THEN
                    sum=0.0_DP
                    DO componentIdx2=1,2
                      DO versionIdx2=1,numberOfVersions
                        sum=sum+normalWave(componentIdx2,versionIdx2)*qBif(versionIdx2)
                      ENDDO !versionIdx
                    ENDDO !componentIdx
                    residualVector%nodalResidual%vector(rowIdx)=sum
                  ELSE
                    residualVector%nodalResidual%vector(rowIdx)= &
                      & (rho/2.0_DP*((qBif(1)/aBif(1))**2.0_DP) + beta(1)*(SQRT(aBif(1)) - SQRT(a0Param(1)))) - &
                      & (rho/2.0_DP*((qBif(versionIdx)/aBif(versionIdx))**2.0_DP) + &
                      & beta(versionIdx)*(SQRT(aBif(versionIdx)) - SQRT(a0Param(versionIdx))))
                  ENDIF
                ENDIF
              ENDDO !versionIdx
            ENDDO !componentIdx
          ENDIF
        ENDIF
      ENDIF !Find branch nodes

    ENDIF !update
    
    EXITS("Characteristic_NodalResidualEvaluate")
    RETURN
999 ERRORSEXITS("Characteristic_NodalResidualEvaluate",err,error)
    RETURN 1
    
  END SUBROUTINE Characteristic_NodalResidualEvaluate

  !
  !================================================================================================================================
  !

  !>Evaluates the Jacobian nodal matrix for a characteristic equation nodal equations set.
  SUBROUTINE Characteristic_NodalJacobianEvaluate(equationsSet,nodeNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet
    INTEGER(INTG), INTENT(IN) :: nodeNumber
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    INTEGER(INTG) :: columnIdx,columnIdx2,componentIdx,derivativeIdx,endRow,esSpecification(3),numberOfVersions,localDOFIdx, &
      & rowIdx,startColumn2,startRow,versionIdx
    REAL(DP) :: qBif(4),aBif(4),a0Param(4),eParam(4),h0Param(4),beta(4),w(2,4),normalWave(2,4),rho
    REAL(DP), POINTER :: independentParameters(:)
    LOGICAL :: update,updateMatrix,boundaryNode
    TYPE(DomainType), POINTER :: dependentDomain
    TYPE(DomainNodesType), POINTER :: dependentDomainNodes
    TYPE(DomainTopologyType), POINTER :: dependentDomainNodes
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix
    TYPE(FieldType), POINTER :: materialsField,dependentField,independentField
    TYPE(FieldVariableType), POINTER :: dependentUVariable,dependentVVariable,independentVariable,materialsUVariable, &
      & materialsVVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("Characteristic_NodalJacobianEvaluate",err,error,*999)
  
    CALL EquationsSet_SpecificationGet(esSpecification,3,esSpecification,err,error,*999)

    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_CHARACTERISTIC_SUBTYPE)
      !OK
    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is not valid for a Poisson equation type of a classical field equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    NULLIFY(equations)
    CALL EquationSet_EquationsGet(equationsSet,equations,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMapping)
    CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
    NULLIFY(nonlinearMapping)
    CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)
    NULLIFY(residualMapping)
    CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,1,residualMapping,err,error,*999)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
    NULLIFY(nonlinearMatrices)
    CALL EquationsMatricesVector_NonlinearMatricesGet(vectorMatrices,nonlinearMatrices,err,error,*999)
    NULLIFY(residualVector)
    CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,1,residualVector,err,error,*999)
    NULLIFY(jacobianMatrix)
    CALL EquationsMatricesResidual_JacobianMatrixGet(residualVector,1,jacobianMatrix,err,error,*999)
    updateMatrix=jacobianMatrix%updateJacobian
 
    update=(updateMatrix)

    IF(update) THEN

      NULLIFY(dependentField)
      CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
      NULLIFY(materialsField)
      CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)
      NULLIFY(independentField)
      CALL EquationsSet_IndependentFieldGet(equationsSet,independentField,err,error,*999)

      NULLIFY(dependentUVariable)
      CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,dependentUVariable,err,error,*999)
      NULLIFY(dependentVVariable)
      CALL Field_VariableGet(dependentField,FIELD_V_VARIABLE_TYPE,dependentVVariable,err,error,*999)
      NULLIFY(dependentDecomposition)
      CALL Field_DecompositionGet(dependentField,dependentDecomposition,err,error,*999)
      NULLIFY(dependentDomain)
      CALL Decomposition_DomainGet(dependentDecomposition,0,dependentDomain,err,error,*999)
      NULLIFY(dependentDomainTopology)
      CALL Domain_DomainTopologyGet(dependentDomain,dependentDomainTopology,err,error,*999)
      NULLIFY(dependentDomainNodes)
      CALL DomainTopology_DomainNodesGet(dependentDomainTopology,dependentDomainNodes,err,error,*999)
      
      NULLIFY(materialsUVariable)
      CALL Field_VariableGet(materialsField,FIELD_U_VARIABLE_TYPE,materialsUVariable,err,error,*999)
      NULLIFY(materialsVVariable)
      CALL Field_VariableGet(materialsField,FIELD_V_VARIABLE_TYPE,materialsVVariable,err,error,*999)

      NULLIFY(independentVariable)
      CALL Field_VariableGet(independentField,FIELD_U_VARIABLE_TYPE,independentVariable,err,error,*999)

      derivativeIdx=1
      normalWave=0.0_DP
      CALL DomainNodes_DerivativeNumberOfVersionsGet(dependentDomainNodes,derivativeIdx,nodeNumber,numberOfVersions,err,error,*999)
      CALL DomainNodes_NodeBoundaryNodeGet(dependentDomainNodes,nodeNumber,boundaryNode,err,error,*999)

      !Get normal wave direction for nodes
      CALL FieldVariable_ParameterSetDataGet(independentVariable,FIELD_VALUES_SET_TYPE,independentParameters,err,error,*999)
      DO componentIdx=1,2
        DO versionIdx=1,numberOfVersions
          CALL FieldVariable_LocalNodeDOFGet(dependentUVariable,versionIdx,derivativeIdx,nodeNumber,componentIdx, &
            & localDOFIdx,err,error,*999)
          normalWave(componentIdx,versionIdx)=independentParameters(localDOFIdx)
        ENDDO !versionIdx
      ENDDO !componentIdx
      CALL FieldVariable_ParameterSetDataRestore(independentVariable,FIELD_VALUES_SET_TYPE,independentParameters,err,error,*999)

      !!!-- F i n d   B r a n c h   N o d e s --!!!
      IF(ABS(normalWave(1,1))>0 .OR. ABS(normalWave(2,1))>0) THEN
        IF(.NOT. boundaryNode) THEN

          !Get material constants
          CALL FieldVariable_ParameterSetGetConstant(materialsUVariable,FIELD_VALUES_SET_TYPE,2,rho,err,error,*999)
          !Get node-based material parameters
          DO versionIdx=1,numberOfVersions
            CALL FieldVariable_ParameterSetGetLocalNode(materialsVVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
              & nodeNumber,1,a0Param(versionIdx),err,error,*999)  
            CALL FieldVariable_ParameterSetGetLocalNode(materialsVVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
              & nodeNumber,2,eParam(versionIdx),err,error,*999)
            CALL FieldVariable_ParameterSetGetLocalNode(materialsVVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
              & nodeNumber,3,h0Param(versionIdx),err,error,*999)
            beta(versionIdx)=(4.0_DP*SQRT(PI)*eParam(versionIdx)*h0Param(versionIdx))/(3.0_DP*a0Param(versionIdx))     
          ENDDO !versionIdx

          !Get current Q & A Values at the node
          DO versionIdx=1,numberOfVersions
            CALL FieldVariable_ParameterSetGetLocalNode(dependentUVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
              & nodeNumber,1,qBif(versionIdx),err,error,*999)
            CALL FieldVariable_ParameterSetGetLocalNode(dependentUVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
              & nodeNumber,2,aBif(versionIdx),err,error,*999)                
            ! If A goes negative during nonlinear iteration, set to A0
            IF (aBif(versionIdx) < a0Param(versionIdx)*0.001_DP) aBif(versionIdx) = a0Param(versionIdx)*0.001_DP
          ENDDO !versionIdx

          !Get extrapolated W for the node
          DO componentIdx=1,2
            DO versionIdx=1,numberOfVersions
              CALL FieldVariable_ParameterSetGetLocalNode(dependentVVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
                & nodeNumber,componentIdx,w(componentIdx,versionIdx),err,error,*999)                
            ENDDO !versionIdx
          ENDDO !componentIdx

          !!!--  J A C O B I A N   M A T R I X  --!!!
          IF(updateMatrix) THEN
            ! Characteristic equations (dW/dU)
            columnIdx=0
            rowIdx=0
            DO componentIdx=1,2
              DO versionIdx=1,numberOfVersions
                IF(ABS(normalWave(componentIdx,versionIdx))>ZERO_TOLERANCE) THEN
                  columnIdx=columnIdx+1
                  rowIdx=rowIdx+1
                  jacobianMatrix%nodalJacobian%matrix(rowIdx,columnIdx)=1.0_DP/aBif(versionIdx)
                ENDIF                
              ENDDO !versionIdx
            ENDDO !componentIdx
            rowIdx=0
            DO componentIdx=1,2
              DO versionIdx=1,numberOfVersions
                IF(ABS(normalWave(componentIdx,versionIdx))>ZERO_TOLERANCE) THEN
                  columnIdx=columnIdx+1
                  rowIdx=rowIdx+1
                  jacobianMatrix%nodalJacobian%matrix(rowIdx,columnIdx)=(-qBif(versionIdx)/(aBif(versionIdx)**2)) &
                    & +normalWave(componentIdx,versionIdx)*SQRT(beta(versionIdx)/(2.0_DP*rho))*(aBif(versionIdx)**(-0.75_DP))
                ENDIF
              ENDDO !versionIdx
            ENDDO !componentIdx

            !Conservation of Mass
            rowIdx=numberOfVersions+1
            columnIdx = 0
            DO componentIdx=1,2
              DO versionIdx=1,numberOfVersions
                IF(ABS(normalWave(componentIdx,versionIdx))>ZERO_TOLERANCE) THEN
                  columnIdx=columnIdx+1
                  jacobianMatrix%nodalJacobian%matrix(rowIdx,columnIdx)=normalWave(componentIdx,versionIdx)
                ENDIF                
              ENDDO !versionIdx
            ENDDO !componentIdx

            !Continuity of Total Pressure (dP/dU)
            startRow=numberOfVersions+2
            endRow=numberOfVersions*2
            startColumn2=numberOfVersions+1
            DO rowIdx=startRow,endRow
              columnIdx=1
              columnIdx2=startColumn2
              DO componentIdx=1,2
                DO versionIdx=1,numberOfVersions
                  IF(ABS(normalWave(componentIdx,versionIdx))>ZERO_TOLERANCE) THEN
                    IF(columnIdx==1) THEN
                      ! dP/dQ
                      jacobianMatrix%nodalJacobian%matrix(rowIdx,columnIdx)=rho*(qBif(1)/(aBif(1)**2.0_DP))
                      ! dP/dA
                      jacobianMatrix%nodalJacobian%matrix(rowIdx,columnIdx2)=beta(1)/(2.0_DP*SQRT(aBif(1))) - &
                        & (rho)*((qBif(1)**2.0_DP)/(aBif(1)**3.0_DP))
                    ELSE IF(columnIdx2==rowIdx) THEN
                      ! dP/dQ
                      jacobianMatrix%nodalJacobian%matrix(rowIdx,columnIdx)=-rho*(qBif(versionIdx)/(aBif(versionIdx)**2.0_DP))
                      ! dP/dA
                      jacobianMatrix%nodalJacobian%matrix(rowIdx,columnIdx2)=-beta(versionIdx)/(2.0_DP*SQRT(aBif(versionIdx))) + &
                        & (rho)*((qBif(versionIdx)**2.0_DP)/(aBif(versionIdx)**3.0_DP))
                    ELSE
                      jacobianMatrix%nodalJacobian%matrix(rowIdx,versionIdx)=0.0_DP
                      jacobianMatrix%nodalJacobian%matrix(rowIdx,columnIdx)=0.0_DP
                    ENDIF
                    columnIdx=columnIdx+1
                    columnIdx2=columnIdx2+1
                  ENDIF
                ENDDO !versionIdx
              ENDDO !componentIdx
            ENDDO !rowIdx

          ENDIF
        ENDIF
      ENDIF !Find branch nodes

    ENDIF !update
       
    EXITS("Characteristic_NodalJacobianEvaluate")
    RETURN
999 ERRORSEXITS("Characteristic_NodalJacobianEvaluate",err,error)
    RETURN 1
    
  END SUBROUTINE Characteristic_NodalJacobianEvaluate

  !
  !================================================================================================================================
  !

  !>Extrapolate W for branch nodes and boundaries .
  SUBROUTINE Characteristic_Extrapolate(solver,currentTime,timeIncrement,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver 
    REAL(DP), INTENT(IN) :: currentTime
    REAL(DP), INTENT(IN) :: timeIncrement
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    TYPE(BasisType), POINTER :: dependentBasis,materialsBasis
    TYPE(DomainType), POINTER :: dependentDomain,materialsDomain
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(EquationsType), POINTER :: equations
    TYPE(FieldType), POINTER ::  dependentField,materialsField,independentField,geometricField
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    REAL(DP) :: w(2,4),qEx(4),aEx(4),xi(1),a0Param(4),h0Param(4),eParam(4),beta(4),normalWave(2,4),elementLengths(4)
    REAL(DP) :: a0Ex(4),h0Ex(4),eEx(4),betaEx(4),f(4),l,friction
    REAL(DP) :: qPrevious,aPrevious,rho,lambda(4)
    REAL(DP) :: elementLength,extrapolationDistance
    INTEGER(INTG) :: nodeIdx,versionIdx,derivativeIdx,elementIdx,elementNumber,versionElementNumber(4),lineNumber
    INTEGER(INTG) :: elementNodeIdx,elementNodeNumber,elementNodeVersion,numberOfVersions,componentIdx,numberOfLocalNodes
    LOGICAL :: overExtrapolated

    ENTERS("Characteristic_Extrapolate",err,error,*999)

    NULLIFY(solverEquations)
    CALL Solver_SolverEquations(solver,solverEquations,err,error,*999)
    NULLIFY(solverMapping)
    CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
    NULLIFY(equationsSet)
    CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    NULLIFY(geometricField)
    CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
    NULLIFY(independentField)
    CALL EquationsSet_IndependentFieldGet(equationsSet,independentField,err,error,*999)
    NULLIFY(materialsField)
    CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)

    NULLIFY(geometricParameters)
    CALL Field_GeometricParametersGet(geometricField,geometricParameters,err,error,*999)
    NULLIFY(geometricVariable)
    CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,geometricVariable,err,error,*999)
    NULLIFY(geometricDecomposition)
    CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
    NULLIFY(geometricDecompositionTopology)
    CALL Decomposition_DecompositionTopologyGet(geometricDecomposition,geometricDecompositionTopology,err,error,*999)
    NULLIFY(geometricDecompositionElements)
    CALL DecompositionTopology_DecompositionElementsGet(geometricDecompositionTopology,geometricDecompositionElements, &
        & err,error,*999)
    
    NULLIFY(dependentUVariable)
    CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,dependentUVariable,err,error,*999)
    NULLIFY(dependentVVariable)
    CALL Field_VariableGet(dependentField,FIELD_V_VARIABLE_TYPE,dependentVVariable,err,error,*999)
    NULLIFY(dependentDecomposition)
    CALL Field_DecompositionGet(dependentField,dependentDecomposition,err,error,*999)
    NULLIFY(dependentDomain)
    CALL Decomposition_DomainGet(dependentDecomposition,0,dependentDomain,err,error,*999)
    NULLIFY(dependentDomainTopology)
    CALL Domain_DomainTopologyGet(dependentDomain,dependentDomainTopology,err,error,*999)
    NULLIFY(dependentDomainNodes)
    CALL DomainTopology_DomainNodesGet(dependentDomainTopology,dependentDomainNodes,err,error,*999)
    CALL DomainNodes_NumberOfNodesGet(dependentDomainNodes,numberOfLocalNodes,err,error,*999)
    NULLIFY(dependentDomainElements)
    CALL DomainTopology_DomainElementsGet(dependentDomainTopology,dependentDomainElements,err,error,*999)

    NULLIFY(materialsUVariable)
    CALL Field_VariableGet(materialsField,FIELD_U_VARIABLE_TYPE,materialsUVariable,err,error,*999)
    NULLIFY(materialsVVariable)
    CALL Field_VariableGet(materialsField,FIELD_V_VARIABLE_TYPE,materialsVVariable,err,error,*999)
    NULLIFY(materialsDecomposition)
    CALL Field_DecompositionGet(materialsField,materialsDecomposition,err,error,*999)
    NULLIFY(materialsDomain)
    CALL Decomposition_DomainGet(materialsDecomposition,0,materialsDomain,err,error,*999)
    NULLIFY(materialsDomainTopology)
    CALL Domain_DomainTopologyGet(materialsDomain,materialsDomainTopology,err,error,*999)
    NULLIFY(materialsDomainNodes)
    CALL DomainTopology_DomainNodesGet(materialsDomainTopology,materialsDomainNodes,err,error,*999)
    NULLIFY(materialsDomainElements)
    CALL DomainTopology_DomainElementsGet(materialsDomainTopology,materialsDomainElements,err,error,*999)
    
    NULLIFY(independentVariable)
    CALL Field_VariableGet(independentField,FIELD_U_VARIABLE_TYPE,independentVariable,err,error,*999)

    NULLIFY(equationsInterpolation)
    CALL Equations_EquationsInterpolationGet(equations,equationsInterpolation,err,error,*999)
    NULLIFY(dependentUInterpParams)
    CALL EquationsInterpolation_DependentParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
      & dependentUInterpParameters,err,error,*999)
    NULLIFY(dependentUInterpPoint)
    CALL EquationsInterpolation_DependentPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,dependentUInterpPoint, &
      & err,error,*999)
    NULLIFY(materialsVInterpParamters)
    CALL EquationsInterpolation_MaterialsParametersGet(equationsInterpolation,FIELD_V_VARIABLE_TYPE, &
      & materialsVInterpParameters,err,error,*999)
    NULLIFY(materialsVInterpPoint)
    CALL EquationsInterpolation_MaterialsPointGet(equationsInterpolation,FIELD_V_VARIABLE_TYPE,materialsVInterpPoint, &
      & err,error,*999)
    
    derivativeIdx=1
    
    !--  L o o p   O v e r   L o c a l  N o d e s  --!
    DO nodeIdx=1,numberOfLocalNodes

      CALL DomainNodes_DerivativeNumberOfVersionsGet(dependentDomainNodes,derivativeIdx,nodeIdx,numberOfVersions,err,error,*999)
      CALL DomainNodes_NodeNumberOfSurroundingElementsGet(dependentDomainNodes,nodeIdx,numberOfSurroundingElements,err,error,*999)

      !Get normal wave direction
      normalWave=0.0_DP
      DO componentIdx=1,2
        DO versionIdx=1,numberOfVersions
          CALL FieldVariable_ParameterSetGetLocalNode(independentVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
            & nodeIdx,componentIdx,normalWave(componentIdx,versionIdx),err,error,*999)
        ENDDO !versionIdx
      ENDDO !componentIdx

      !-- F i n d   B r a n c h   a n d   B o u n d a r y    N o d e s --!
      IF(ABS(normalWave(1,1)) > ZERO_TOLERANCE .OR. ABS(normalWave(2,1))> ZERO_TOLERANCE) THEN
        !Get constant material parameters
        CALL FieldVariable_ParameterSetGetConstant(materialsUVariable,FIELD_VALUES_SET_TYPE,2,rho,err,error,*999)
        
        overExtrapolated = .FALSE.
        !-- G e t   E l e m e n t   L e n g t h s --!
        elementLengths = 0.0_DP
        DO elementIdx=1,numberOfSurroundingElements
          CALL DomainNodes_NodeSurroundingElementGet(dependentDomainNodes,elementIdx,nodesIdx,elementNumber,err,error,*999)
          CALL DecompositionElements_ElementLineNumberGet(geometricDecompositionElements,1,elementNumber,lineNumber,err,error,*999)
          ! Get the line lengths to extrapolate at equidistant points from the branch node
          CALL FieldGeometricParameters_LineLengthGet(geometricParameters,lineNumber,elementLength,err,error,*999)
          !Loop over the nodes on this (surrounding) element
          NULLIFY(dependentBasis)
          CALL DomainElements_ElementBasisGet(dependentDomainElements,elementNumber,dependentBasis,err,error,*999)
          NULLIFY(materialBasis)
          CALL DomainElements_ElementBasisGet(materialsDomainElements,elementNumber,materialsBasis,err,error,*999)
          CALL Basis_NumberOfNodesGet(dependentBasis,numberOfElementNodes,err,error,*999)
          DO elementNodeIdx=1,numberOfElementNodes
            CALL DomainElements_ElementNodeGet(dependentDomainElements,elementNodeIdx,elementNumber,elementNodeNumber, &
              & err,error,*999)
            !Check that this node is the same as the current iterative node
            IF(elementNodeNumber==nodeIdx) THEN
              !Loop over the versions to find the element index that matches the version
              DO versionIdx=1,numberOfVersions
                !Version number for the local element node
                CALL DomainElements_ElementVersionGet(dependentDomainElements,1,elementNodeIdx,elementNumber,elementNodeVersion, &
                  & err,error,*999)
                IF(elementNodeVersion==versionIdx) THEN
                  versionElementNumber(versionIdx)=elementNumber
                  elementLengths(versionIdx) = elementLength
                ENDIF
              ENDDO !versionIDx
            ENDIF
          ENDDO !elementNodeIdx
        ENDDO !elementIdx

        !-- E x t r a p o l a t e   Q   a n d   A    V a l u e s --!
        ! --------------------------------------------------------------
        ! Extrapolate along the characteristic curve a distance x - lambda*dt from node location (x) to get 
        ! values for w(t) from Q,A(t-delta(t)). Note that since the characteristic solver runs before the 
        ! Navier-Stokes solver, 'previous' values are still in the 'current' field at this time-step as the 
        ! time integration occurs as part of the Navier-Stokes solution.
        DO componentIdx=1,2
          DO versionIdx=1,numberOfVersions                         
            IF(ABS(normalWave(componentIdx,versionIdx))> ZERO_TOLERANCE) THEN
              
              ! Get materials values at node
              CALL FieldVariable_ParameterSetGetLocalNode(materialsVVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
                & nodeIdx,1,a0Param(versionIdx),err,error,*999)
              CALL FieldVariable_ParameterSetGetLocalNode(materialsVVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
                & nodeIdx,2,eParam(versionIdx),err,error,*999)            
              CALL FieldVariable_ParameterSetGetLocalNode(materialsVVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
                & nodeIdx,3,h0Param(versionIdx),err,error,*999)
              beta(versionIdx) = (4.0_DP*SQRT(PI)*eParam(versionIdx)*h0Param(versionIdx))/(3.0_DP*a0Param(versionIdx))

              ! Get previous Q,A values at node
              CALL FieldVariable_ParameterSetGetLocalNode(dependentUVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
                & nodeIdx,1,qPrevious,err,error,*999)
              CALL FieldVariable_ParameterSetGetLocalNode(dependentUVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
                & nodeIdx,2,aPrevious,err,error,*999)
              
              ! Calculate wave speed
              lambda(versionIdx) = qPrevious/aPrevious + normalWave(componentIdx,versionIdx)*(aPrevious**0.25)* &
                & SQRT(beta(versionIdx)/(2.0_DP*rho))
              ! Check that lambda(1) > 0, lambda(2) < 0
              IF (lambda(versionIdx)*normalWave(componentIdx,versionIdx) < 0.0_DP) THEN
                CALL FlagError("Subcritical 1D system violated.",err,error,*999)
              ENDIF
              
              ! Calculate extrapolation distance and xi location
              extrapolationDistance = (timeIncrement)*lambda(versionIdx)
              !  Convert to xi-space within the element
              IF((normalWave(componentIdx,versionIdx)>ZERO_TOLERANCE)) THEN
                ! Parent branch / outlet boundary
                xi(1)=1.0_DP - extrapolationDistance/(elementLengths(versionIdx))
              ELSE
                ! Daughter branch / inlet boundary
                xi(1)=0.0_DP - extrapolationDistance/(elementLengths(versionIdx))
              ENDIF
              IF (xi(1) > 1.0_DP .OR. xi(1) < 0.0_DP) THEN
                CALL FlagWarning("1D extrapolation location outside of element xi space. Reduce time increment.", &
                  & err,error,*999)
                overExtrapolated = .TRUE.
              ENDIF
              
              ! Get Q,A values at extrapolated xi locations
              CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,versionElementNumber(versionIdx), &
                & dependentUInterpParameters,err,error,*999)
              CALL Field_InterpolateXi(NO_PART_DERIV,xi,dependentUInterpPoint,err,error,*999)
              qEx(versionIdx)=dependentUInterpParameters%values(1,NO_PART_DERIV)
              aEx(versionIdx)=dependentUInterpParameters%values(2,NO_PART_DERIV)
              ! Get spatially varying material values at extrapolated xi locations
              CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,versionElementNumber(versionIdx), &
                & materialsVInterpParameters,err,error,*999)
              CALL Field_InterpolateXi(NO_PART_DERIV,xi,matrieralsVInterpPoint,err,error,*999)
              a0Ex(versionIdx)=materialsVInterpParameters%values(1,NO_PART_DERIV)
              eEx(versionIdx)=materialsVInterpParameters%values(2,NO_PART_DERIV)
              h0Ex(versionIdx)=materialsVInterpParameters%values(3,NO_PART_DERIV)
              betaEx(versionIdx) = (4.0_DP*SQRT(PI)*eEx(versionIdx)*h0Ex(versionIdx))/(3.0_DP*a0Ex(versionIdx))
              ! Calculate friction term if necessary
              f(versionIdx) = -qEx(versionIdx)/(aEx(versionIdx)**2.0_DP)
            ENDIF
          ENDDO !versionIdx
        ENDDO !componentIdx

        !Calculate W
        w(:,:)=0.0_DP
        DO componentIdx=1,2
          DO versionIdx=1,numberOfVersions
            IF(ABS(normalWave(componentIdx,versionIdx))>ZERO_TOLERANCE) THEN
              ! w(t+delta(t)) = W_extrap(t)
              w(componentIdx,versionIdx)= ((qEx(versionIdx)/aEx(versionIdx))+ &
                & normalWave(componentIdx,versionIdx)*4.0_DP*SQRT(betaEx(versionIdx)/(2.0_DP*rho))* &
                & (aEx(versionIdx)**(0.25_DP) - (a0Ex(versionIdx))**(0.25_DP)))
              
              ! Add friction term if not neglected
              l = (1.0_DP/(qEx(versionIdx)/aEx(versionIdx) +  &
                & normalWave(componentIdx,versionIdx)*aEx(versionIdx)**0.25_DP*SQRT(betaEx(versionIdx)/(2.0_DP*rho))))
              friction = timeIncrement*l*f(versionIdx)
              !w(componentIdx,versionIdx)= w(componentIdx,versionIdx) + friction
              
              ! Check extrapolated wave speed is coherent
              lambda(versionIdx) = qEx(versionIdx)/aEx(versionIdx) + normalWave(componentIdx,versionIdx)* &
                & (aEx(versionIdx)**0.25)*SQRT(beta(versionIdx)/(2.0_DP*rho))
              IF (lambda(versionIdx)*normalWave(componentIdx,versionIdx) < -ZERO_TOLERANCE ) THEN
                CALL FlagError("Subcritical 1D system violated.",err,error,*999)
              ENDIF
              
              IF (.NOT. overExtrapolated) THEN
                CALL FieldVariable_ParameterSetUpdateLocalNode(dependentVVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
                  & nodeIdx,componentIdx,w(componentIdx,versionIdx),err,error,*999)
              ENDIF
            ENDIF
          ENDDO !versionIdx
        ENDDO !componentIdx
      ENDIF ! branch or boundary node
    ENDDO !Loop over nodes

    EXITS("Characteristic_Extrapolate")
    RETURN
999 ERRORSEXITS("Characteristic_Extrapolate",err,error)
    RETURN 1
    
  END SUBROUTINE Characteristic_Extrapolate

  !
  !================================================================================================================================
  !

  !>Calculate Characteristic (W) values based on dependent field values
  SUBROUTINE Characteristic_PrimitiveToCharacteristic(equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer the equations set
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: componentIdx,derivativeIdx,esSpecification(3),nodeIdx,numberOfNodes,numberOfVersions,versionIdx
    REAL(DP) :: a0Param,aCurrent(4),beta,eParam,h0Param,normalWave,qCurrent(4),w(2,4)
    LOGICAL :: boundaryNode
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(DomainTopologyType), POINTER :: domainToplogy
    TYPE(FieldType), POINTER ::  dependentField,materialsField,independentField
    TYPE(FieldVariableType), POINTER :: dependentUVariable,dependentVVariable,independentVariable,materialsVVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("Characteristic_PrimitiveToCharacteristic",err,error,*999)

    CALL EquationsSet_SpecificationGet(esSpecification,3,esSpecification,err,error,*999)

    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE)
      !OK
    CASE DEFAULT
      localError="The third equations set specification of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is not valid for a call to Characteristic_PrimitiveToCharacteristic."
    END SELECT
    
    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
    NULLIFY(materialsField)
    CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)
    NULLIFY(independentField)
    CALL EquationsSet_IndependentFieldGet(equationsSet,independentField,err,error,*999)

    NULLIFY(dependentUVariable)
    CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,dependentUVariable,err,error,*999)
    NULLIFY(dependentVVariable)
    CALL Field_VariableGet(dependentField,FIELD_V_VARIABLE_TYPE,dependentVVariable,err,error,*999)
    NULLIFY(decomposition)
    CALL Field_DecompositionGet(dependentField,decomposition,err,error,*999)
    NULLIFY(domain)
    CALL Decomposition_DomainGet(decomposition,0,domain,err,error,*999)
    NULLIFY(domainTopology)
    CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
    NULLIFY(domainNodes)
    CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)

    NULLIFY(materialsVVariable)
    CALL Field_VariableGet(materialsField,FIELD_V_VARIABLE_TYPE,materialsVVariable,err,error,*999)
    
    NULLIFY(independentVariable)
    CALL Field_VariableGet(independentField,FIELD_U_VARIABLE_TYPE,independentVariable,err,error,*999)
    
    !--  L o o p   O v e r   L o c a l  N o d e s  --!!!
    CALL DomainNodes_NumberOfNodesGet(domainNodes,numberOfNodes,err,error,*999)
    DO nodeIdx=1,numberOfNodes
      derivativeIdx = 1
      CALL DomainNodes_DerivativeNumberOfVersionsGet(domainNodes,derivativeIdx,nodeIdx,numberOfVersions,err,error,*999)
      CALL DomainNodes_NodeBoundaryNodeGet(domainNodes,nodeIdx,boundaryNode,err,error,*999)
      !-- F i n d    B r a n c h    N o d e s --!!!
      IF(numberOfVersions > 1 .AND. .NOT. boundaryNode) THEN
        DO componentIdx=1,2
          DO versionIdx=1,numberOfVersions
            CALL FieldVariable_ParameterSetGetLocalNode(independentVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
              & nodeIdx,componentIdx,normalWave,err,error,*999)
            IF(ABS(normalWave)>ZERO_TOLERANCE) THEN
              !Get material parameters
              CALL FieldVariable_ParameterSetGetLocalNode(materialsVVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
                & nodeIdx,1,a0Param,err,error,*999)  
              CALL FieldVariable_ParameterSetGetLocalNode(materialsVVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
                & nodeIdx,2,eParam,err,error,*999)                
              CALL FieldVariable_ParameterSetGetLocalNode(materialsVVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
                & nodeIdx,3,h0Param,err,error,*999)
              beta=(4.0_DP*SQRT(PI)*eParam*h0Param)/(3.0_DP*a0Param)     

              ! Get current Q,A values
              CALL FieldVariable_ParameterSetGetLocalNode(dependentUVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
                & nodeIdx,1,qCurrent(versionIdx),err,error,*999)         
              CALL FieldVariable_ParameterSetGetLocalNode(dependentUVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
                & nodeIdx,2,aCurrent(versionIdx),err,error,*999)

              ! Calculate the characteristic based on current Q,A values
              w(componentIdx,versionIdx)= ((qCurrent(versionIdx)/aCurrent(versionIdx))+ &
                & normalWave*4.0_DP*SQRT(((beta)))*(aCurrent(versionIdx)**(0.25_DP) - (a0Param)**(0.25_DP)))
              
              !Update W values
              CALL FieldVariable_ParameterSetUpdateLocalNode(dependentVVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
                & nodeIdx,componentIdx,w(componentIdx,versionIdx),err,error,*999)
            ENDIF
          ENDDO !versionIdx
        ENDDO !componentIdx
      ENDIF ! branch check
    ENDDO !nodeIdx

    EXITS("Characteristic_PrimitiveToCharacteristic")
    RETURN
999 ERRORSEXITS("Characteristic_PrimitiveToCharacteristic",err,error)
    RETURN 1
    
  END SUBROUTINE Characteristic_PrimitiveToCharacteristic

  !
  !================================================================================================================================
  !      

END MODULE CharacteristicEquationsRoutines
