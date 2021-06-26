!> \file
!> \author Chris Bradley
!> \brief This module handles all equations mapping routines.
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

!>This module handles all equations mapping routines.
MODULE EquationsMappingRoutines

  USE BaseRoutines
  USE DomainMappings
  USE EquationsAccessRoutines
  USE EquationsMappingAccessRoutines
  USE EquationsSetAccessRoutines
  USE FieldRoutines
  USE FieldAccessRoutines
  USE InputOutput
  USE ISO_VARYING_STRING
  USE Kinds
  USE Lists
  USE Strings
  USE Types

#include "macros.h"  
 
  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  INTERFACE EquationsMappingVector_DynamicMatricesSet
    MODULE PROCEDURE EquationsMappingVector_DynamicMatricesSet1
    MODULE PROCEDURE EquationsMappingVector_DynamicMatricesSet2
  END INTERFACE EquationsMappingVector_DynamicMatricesSet
  
  INTERFACE EquationsMappingVector_DynamicMatricesCoefficientsSet
    MODULE PROCEDURE EquationsMappingVector_DynamicMatricesCoefficientsSet1
    MODULE PROCEDURE EquationsMappingVector_DynamicMatricesCoefficientsSet2
  END INTERFACE EquationsMappingVector_DynamicMatricesCoefficientsSet
  
  INTERFACE EquationsMappingVector_LinearMatricesCoefficientsSet
    MODULE PROCEDURE EquationsMappingVector_LinearMatricesCoefficientsSet0
    MODULE PROCEDURE EquationsMappingVector_LinearMatricesCoefficientsSet1
  END INTERFACE EquationsMappingVector_LinearMatricesCoefficientsSet
  
  INTERFACE EquationsMappingVector_LinearMatricesVariableTypesSet
    MODULE PROCEDURE EquationsMappingVector_LinearMatricesVariableTypesSet0
    MODULE PROCEDURE EquationsMappingVector_LinearMatricesVariableTypesSet1
  END INTERFACE EquationsMappingVector_LinearMatricesVariableTypesSet
  
  INTERFACE EquationsMappingVector_ResidualCoefficientsSet
    MODULE PROCEDURE EquationsMappingVector_ResidualCoefficientsSet0
    MODULE PROCEDURE EquationsMappingVector_ResidualCoefficientsSet1
  END INTERFACE EquationsMappingVector_ResidualCoefficientsSet
  
  INTERFACE EquationsMappingVector_ResidualVariableTypesSet
    MODULE PROCEDURE EquationsMappingVector_ResidualVariableTypesSet0
    MODULE PROCEDURE EquationsMappingVector_ResidualVariableTypesSet1
  END INTERFACE EquationsMappingVector_ResidualVariableTypesSet
  
  INTERFACE EquationsMappingVector_SourcesCoefficientsSet
    MODULE PROCEDURE EquationsMappingVector_SourcesCoefficientsSet0
    MODULE PROCEDURE EquationsMappingVector_SourcesCoefficientsSet1
  END INTERFACE EquationsMappingVector_SourcesCoefficientsSet
  
  INTERFACE EquationsMappingVector_SourcesVariableTypesSet
    MODULE PROCEDURE EquationsMappingVector_SourcesVariableTypesSet0
    MODULE PROCEDURE EquationsMappingVector_SourcesVariableTypesSet1
  END INTERFACE EquationsMappingVector_SourcesVariableTypesSet
  
  PUBLIC EquationsMapping_ScalarDestroy

  PUBLIC EquationsMapping_VectorCreateFinish,EquationsMapping_VectorCreateStart

  PUBLIC EquationsMapping_VectorDestroy

  PUBLIC EquationsMappingVector_DynamicMatricesSet

  PUBLIC EquationsMappingVector_DynamicMatricesCoefficientsSet

  PUBLIC EquationsMappingVector_DynamicVariableTypeSet

  PUBLIC EquationsMappingVector_LinearMatricesCoefficientsSet

  PUBLIC EquationsMappingVector_LinearMatricesVariableTypesSet
  
  PUBLIC EquationsMappingVector_NumberOfLinearMatricesSet

  PUBLIC EquationsMappingVector_NumberOfResidualsSet

  PUBLIC EquationsMappingVector_NumberOfSourcesSet

  PUBLIC EquationsMappingVector_ResidualCoefficientsSet

  PUBLIC EquationsMappingVector_ResidualNumberOfVariablesSet
  
  PUBLIC EquationsMappingVector_ResidualVariableTypesSet

  PUBLIC EquationsMappingVector_RHSCoefficientSet

  PUBLIC EquationsMappingVector_RHSVariableTypeSet

  PUBLIC EquationsMappingVector_SourcesCoefficientsSet

  PUBLIC EquationsMappingVector_SourcesVariableTypesSet

CONTAINS

  !
  !================================================================================================================================
  !

  !>Destroy an scalar equations mapping.
  SUBROUTINE EquationsMapping_ScalarDestroy(scalarMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingScalarType), POINTER :: scalarMapping !<A pointer the scalar equations mapping to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsMapping_ScalarDestroy",err,error,*999)

    IF(.NOT.ASSOCIATED(scalarMapping)) CALL FlagError("Scalar equations mapping is not associated.",err,error,*999)
    
    CALL EquationsMapping_ScalarFinalise(scalarMapping,err,error,*999)
        
    EXITS("EquationsMapping_ScalarDestroy")
    RETURN
999 ERRORSEXITS("EquationsMapping_ScalarDestroy",err,error)    
    RETURN 1
   
  END SUBROUTINE EquationsMapping_ScalarDestroy

  !
  !================================================================================================================================
  !

  !>Finalises the scalar equations mapping and deallocates all memory.
  SUBROUTINE EquationsMapping_ScalarFinalise(scalarMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingScalarType), POINTER :: scalarMapping !<A pointer to the scalar equations mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsMapping_ScalarFinalise",err,error,*999)

    IF(ASSOCIATED(scalarMapping)) THEN
      CALL EquationsMappingScalar_FunctionsFinalise(scalarMapping%functionMappings,err,error,*999)
      CALL EquationsMappingScalar_NormsFinalise(scalarMapping%normMappings,err,error,*999)
      CALL EquationsMappingScalar_DotProductsFinalise(scalarMapping%dotProductMappings,err,error,*999)
      CALL EquationsMappingScalar_QuadraticsFinalise(scalarMapping%quadraticMappings,err,error,*999)
      DEALLOCATE(scalarMapping)
    ENDIF
    
    EXITS("EquationsMapping_ScalarFinalise")
    RETURN
999 ERRORSEXITS("EquationsMapping_ScalarFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMapping_ScalarFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the scalar equations mapping.
  SUBROUTINE EquationsMapping_ScalarInitialise(scalarEquations,err,error,*)

    !Argument variables
    TYPE(EquationsScalarType), POINTER :: scalarEquations !<A pointer to the scalar equations to initialise the scalar equations mapping for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("EquationsMapping_ScalarInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(scalarEquations)) CALL FlagError("Scalar equations is not associated.",err,error,*998)
    IF(ASSOCIATED(scalarEquations%scalarMapping)) CALL FlagError("Scalar equations mapping is already associated.",err,error,*998)
     
    ALLOCATE(scalarEquations%scalarMapping,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate scalar equations scalar mapping.",err,error,*999)
    scalarEquations%scalarMapping%scalarEquations=>scalarEquations
    scalarEquations%scalarMapping%scalarMappingFinished=.FALSE.
    NULLIFY(scalarEquations%scalarMapping%scalarMatrices)
    NULLIFY(scalarEquations%scalarMapping%functionMappings)
    NULLIFY(scalarEquations%scalarMapping%normMappings)
    NULLIFY(scalarEquations%scalarMapping%dotProductMappings)
    NULLIFY(scalarEquations%scalarMapping%quadraticMappings)
    NULLIFY(scalarEquations%scalarMapping%createValuesCache)
    !CALL EquationsMapping_ScalarCreateValuesCacheInitialise(scalarEquations%scalarMapping,err,error,*999)        
       
    EXITS("EquationsMapping_ScalarInitialise")
    RETURN
999 CALL EquationsMapping_ScalarFinalise(scalarEquations%scalarMapping,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMapping_ScalarInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMapping_ScalarInitialise

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating an equations vector mapping
  SUBROUTINE EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations vector mapping
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dynamicVariableType,linearMatrixVariableType,linearVariableType,matrixIdx,numberOfResidualVariables, &
      & residualIdx,residualVariableType,rhsVariableType,sourceIdx,sourceVariableType,variableIdx
    LOGICAL :: isResidualType
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(FieldType), POINTER :: dependentField,sourceField
    TYPE(FieldVariableType), POINTER :: dependentVariable,dynamicVariable,linearVariable,residualVariable,rhsVariable, &
      & sourceVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsMapping_VectorCreateFinish",err,error,*999)

    CALL EquationsMappingVector_AssertNotFinished(vectorMapping,err,error,*999)
 
    NULLIFY(createValuesCache)
    CALL EquationsMappingVector_CreateValuesCacheGet(vectorMapping,createValuesCache,err,error,*999)
    NULLIFY(vectorEquations)
    CALL EquationsMappingVector_VectorEquationsGet(vectorMapping,vectorEquations,err,error,*999)
    NULLIFY(equations)
    CALL EquationsVector_EquationsGet(vectorEquations,equations,err,error,*999)
    NULLIFY(equationsSet)
    CALL Equations_EquationsSetGet(equations,equationsSet,err,error,*999)
    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)

    !Check that all the variables have been mapped properly
    SELECT CASE(equations%timeDependence)
    CASE(EQUATIONS_STATIC,EQUATIONS_QUASISTATIC)
      SELECT CASE(equations%linearity)
      CASE(EQUATIONS_LINEAR,EQUATIONS_NONLINEAR_BCS)
        !Static, linear equations
        !Check the number of linear matrices
        IF(createValuesCache%numberOfLinearMatrices<1) THEN
          localError="The number of linear matrices of "// &
            & TRIM(NumberToVString(createValuesCache%numberOfLinearMatrices,"*",err,error))// &
            & " is invalid. The number of linear equations matrices must be >= 1 for static linear equations."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        !Check we have the linear variable types in the dependent field
        DO matrixIdx=1,createValuesCache%numberOfLinearMatrices
          CALL EquationsMappingVectorCVC_LinearMatrixVariableTypeGet(createValuesCache,matrixIdx,linearVariableType,err,error,*999)
          NULLIFY(linearVariable)
          CALL Field_VariableGet(dependentField,linearVariableType,linearVariable,err,error,*999)
        ENDDO !matrixIdx
      CASE(EQUATIONS_NONLINEAR)
        !Static, nonlinear equations
        !Check the number of residuals
        IF(createValuesCache%numberOfResiduals<1) THEN
          localError="The number of residuals of "//TRIM(NumberToVString(createValuesCache%numberOfResiduals,"*",err,error))// &
            & " is invalid. The number of residuals must be >= 1 for static nonlinear equations."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        !Check the residuals
        DO residualIdx=1,createValuesCache%numberOfResiduals
          CALL EquationsMappingVectorCVC_NumberOfResidualVariablesGet(createValuesCache,residualIdx,numberOfResidualVariables, &
            & err,error,*999)
          !Check we have some variables
          IF(numberOfResidualVariables<1) THEN
            localError="The number of residual variables of "// &
              & TRIM(NumberToVString(numberOfResidualVariables,"*",err,error))//" for residual number "// &
              & TRIM(NumberToVString(residualIdx,"*",err,error))// &
              & " is invalid. The number of residual variables must be >= 1 for residuals in static nonlinear equations."
            CALL FlagError(localError,err,error,*999)
          ENDIF
          DO variableIdx=1,numberOfResidualVariables
            !Check that we have the variables in the dependent field
            CALL EquationsMappingVectorCVC_ResidualVariableTypeGet(createValuesCache,variableIdx,residualIdx, &
              & residualVariableType,err,error,*999)
            NULLIFY(dependentVariable)
            CALL Field_VariableGet(dependentField,residualVariableType,dependentVariable,err,error,*999)
          ENDDO !variableIdx
        ENDDO !residualIdx
        !Check any linear variable types are in the dependent field
        DO matrixIdx=1,createValuesCache%numberOfLinearMatrices
          CALL EquationsMappingVectorCVC_LinearMatrixVariableTypeGet(createValuesCache,matrixIdx,linearVariableType,err,error,*999)
          NULLIFY(linearVariable)
          CALL Field_VariableGet(dependentField,linearVariableType,linearVariable,err,error,*999)
        ENDDO !matrixIdx
      CASE DEFAULT
        localError="The equations linearity type of "//TRIM(NumberToVString(equations%linearity,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_FIRST_ORDER_DYNAMIC,EQUATIONS_SECOND_ORDER_DYNAMIC)
      !Check that we have a dynamic variable
      CALL EquationsMappingVectorCVC_DynamicVariableTypeGet(createValuesCache,dynamicVariableType,err,error,*999)
      NULLIFY(dynamicVariable)
      CALL Field_VariableGet(dependentField,dynamicVariableType,dynamicVariable,err,error,*999)
      SELECT CASE(equations%linearity)
      CASE(EQUATIONS_LINEAR,EQUATIONS_NONLINEAR_BCS)
        !Dynamic, linear equations
        !Check any linear variable types are in the dependent field
        DO matrixIdx=1,createValuesCache%numberOfLinearMatrices
          CALL EquationsMappingVectorCVC_LinearMatrixVariableTypeGet(createValuesCache,matrixIdx,linearVariableType,err,error,*999)
          NULLIFY(linearVariable)
          CALL Field_VariableGet(dependentField,linearVariableType,linearVariable,err,error,*999)
        ENDDO !matrixIdx
      CASE(EQUATIONS_NONLINEAR)
        !Check the residuals
        IF(createValuesCache%numberOfResiduals<1) THEN
          localError="The number of residuals of "//TRIM(NumberToVString(createValuesCache%numberOfResiduals,"*",err,error))// &
            & " is invalid. The number of residuals must be >= 1 for dynamic nonlinear equations."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        isResidualType=.FALSE.
        DO residualIdx=1,createValuesCache%numberOfResiduals
          CALL EquationsMappingVectorCVC_NumberOfResidualVariablesGet(createValuesCache,residualIdx,numberOfResidualVariables, &
            & err,error,*999)
          !Check we have some variables
          IF(numberOfResidualVariables<1) THEN
            localError="The number of residual variables of "// &
              & TRIM(NumberToVString(numberOfResidualVariables,"*",err,error))//" for residual number "// &
              & TRIM(NumberToVString(residualIdx,"*",err,error))// &
              & " is invalid. The number of residual variables must be >= 1 for residuals in static nonlinear equations."
            CALL FlagError(localError,err,error,*999)
          ENDIF
          DO variableIdx=1,numberOfResidualVariables
            !Check that we have the variables in the dependent field
            CALL EquationsMappingVectorCVC_ResidualVariableTypeGet(createValuesCache,variableIdx,residualIdx, &
              & residualVariableType,err,error,*999)
            IF(residualVariableType==dynamicVariableType) isResidualType=.TRUE. 
            NULLIFY(residualVariable)
            CALL Field_VariableGet(dependentField,residualVariableType,residualVariable,err,error,*999)
          ENDDO !variableIdx          
        ENDDO !residualIdx
        !Check that at least one of the residual variables correspond to the dynamic variable.
        IF(.NOT.isResidualType) THEN
          localError="The specified dynamic variable type of "//TRIM(NumberToVString(dynamicVariableType,"*",err,error))// &
            & " is not the same as any residual variable type in dynamic nonlinear equations."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        !Check any linear variable types are in the dependent field
        DO matrixIdx=1,createValuesCache%numberOfLinearMatrices
          CALL EquationsMappingVectorCVC_LinearMatrixVariableTypeGet(createValuesCache,matrixIdx,linearVariableType,err,error,*999)
          NULLIFY(linearVariable)
          CALL Field_VariableGet(dependentField,linearVariableType,linearVariable,err,error,*999)
        ENDDO !matrixIdx
     CASE DEFAULT
        localError="The equations linearity type of "//TRIM(NumberToVString(equations%linearity,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The equations time dependence type of "//TRIM(NumberToVString(equations%timeDependence,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    !Check the linear matrices variable types
    DO matrixIdx=1,createValuesCache%numberOfLinearMatrices
      IF(createValuesCache%linearMatrixVariableTypes(matrixIdx)==0) THEN
        localError="Invalid equations mapping. The linear matrix variable type is not set for linear matrix number "//&
          & TRIM(NumberToVString(matrixIdx,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDDO !matrixIdx
    !Check any RHS variables
    IF(createValuesCache%rhsVariableType/=0) THEN
      !Check that the RHS variable exists in the dependent field
      rhsVariableType=createValuesCache%rhsVariableType
      NULLIFY(rhsVariable)
      CALL Field_VariableGet(dependentField,rhsVariableType,rhsVariable,err,error,*999)      
      !Check the RHS variable is not mapped to any LHS matrices or vectors
      IF(rhsVariableType==createValuesCache%dynamicVariableType) THEN
        localError="Invalid equations mapping. The dependent field variable type of "// &
          & TRIM(NumberToVString(rhsVariableType,"*",err,error))// &
          & " is mapped to both the RHS and the dynamic matrices."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      DO matrixIdx=1,createValuesCache%numberOfLinearMatrices
        CALL EquationsMappingVectorCVC_LinearMatrixVariableTypeGet(createValuesCache,matrixIdx,linearMatrixVariableType, &
          & err,error,*999)
        IF(rhsVariableType==linearMatrixVariableType) THEN
          localError="Invalid equations mapping. The dependent field variable type of "// &
            & TRIM(NumberToVString(rhsVariableType,"*",err,error))// &
            & " is mapped to both the RHS and linear matrix number "//TRIM(NumberToVString(matrixIdx,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDDO !matrixIdx
      DO residualIdx=1,createValuesCache%numberOfResiduals
        CALL EquationsMappingVectorCVC_NumberOfResidualVariablesGet(createValuesCache,residualIdx,numberOfResidualVariables, &
          & err,error,*999)
        DO variableIdx=1,numberOfResidualVariables
          CALL EquationsMappingVectorCVC_ResidualVariableTypeGet(createValuesCache,variableIdx,residualIdx, &
            & residualVariableType,err,error,*999)
          IF(rhsVariableType==residualVariableType) THEN
            localError="Invalid equations mapping. The dependent field variable type of "// &
              & TRIM(NumberToVString(rhsVariableType,"*",err,error))// &
              & " is mapped to both the RHS and variable number "//TRIM(NumberToVString(variableIdx,"*",err,error))// &
              & " of residual number "//TRIM(NumberToVString(residualIdx,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ENDDO !variableIdx          
      ENDDO !residualIdx
    ENDIF
    IF(createValuesCache%numberOfSources/=0) THEN
      !Check the source variables exists     
      NULLIFY(sourceField)
      CALL EquationsSet_SourceFieldGet(equationsSet,sourceField,err,error,*999)
      DO sourceIdx=1,createValuesCache%numberOfSources
        CALL EquationsMappingVectorCVC_SourceVariableTypeGet(createValuesCache,sourceIdx,sourceVariableType,err,error,*999)
        NULLIFY(sourceVariable)
        CALL Field_VariableGet(sourceField,sourceVariableType,sourceVariable,err,error,*999)
      ENDDO !sourceIdx
    ENDIF
    !Now calculate the equations mapping and clean up
    CALL EquationsMappingVector_Calculate(vectorMapping,err,error,*999)
    CALL EquationsMappingVector_CreateValuesCacheFinalise(vectorMapping%createValuesCache,err,error,*999)
    vectorMapping%vectorMappingFinished=.TRUE.
       
    EXITS("EquationsMapping_VectorCreateFinish")
    RETURN
999 ERRORSEXITS("EquationsMapping_VectorCreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMapping_VectorCreateFinish

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating an equations vector mapping for an equations
  SUBROUTINE EquationsMapping_VectorCreateStart(vectorEquations,lhsVariableType,vectorMapping,err,error,*)

    !Argument variables
    TYPE(EquationsVectorType), POINTER :: vectorEquations !<A pointer to the equations to create the equations vector mapping from.
    INTEGER(INTG), INTENT(IN) :: lhsVariableType !<The variable type associated with the equations left hand side (LHS). The LHS variable controls the row mappings for vector equations sets. 
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<On return, a pointer to the equations vector mapping. This must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dofIdx,rowIdx
    TYPE(DomainMappingType), POINTER :: lhsDomainMapping
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(FieldType), POINTER :: dependentField
    TYPE(FieldVariableType), POINTER :: lhsVariable
    
    ENTERS("EquationsMapping_VectorCreateStart",err,error,*998)

    IF(ASSOCIATED(vectorMapping)) CALL FlagError("Vector equations mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorEquations)) CALL FlagError("Vector equations is not associated.",err,error,*999)
    NULLIFY(equations)
    CALL EquationsVector_EquationsGet(vectorEquations,equations,err,error,*999)
    CALL Equations_AssertIsFinished(equations,err,error,*999)
    NULLIFY(equationsSet)
    CALL Equations_EquationsSetGet(equations,equationsSet,err,error,*999)
    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)

    NULLIFY(lhsVariable)
    CALL Field_VariableGet(dependentField,lhsVariableType,lhsVariable,err,error,*999)
    NULLIFY(lhsDomainMapping)
    CALL FieldVariable_DomainMappingGet(lhsVariable,lhsDomainMapping,err,error,*999)
     
    CALL EquationsMapping_VectorInitialise(vectorEquations,err,error,*999)
    CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
    NULLIFY(lhsMapping)
    CALL EquationsMappingVector_LHSMappingGet(vectorMapping,lhsMapping,err,error,*999)
    lhsMapping%lhsVariableType=lhsVariableType
    lhsMapping%lhsVariable=>lhsVariable
    lhsMapping%rowDofsMapping=>lhsDomainMapping
    CALL FieldVariable_NumberOfDOFsGet(lhsVariable,lhsMapping%numberOfRows,err,error,*999)
    CALL FieldVariable_TotalNumberOfDOFsGet(lhsVariable,lhsMapping%totalNumberOfRows,err,error,*999)
    CALL FieldVariable_NumberOfGlobalDOFsGet(lhsVariable,lhsMapping%numberOfGlobalRows,err,error,*999)
   !Allocate and set up the row mappings
    ALLOCATE(lhsMapping%lhsDOFToEquationsRowMap(lhsVariable%totalNumberOfDOFs),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate LHS DOF to equations row map.",err,error,*999)
    ALLOCATE(lhsMapping%equationsRowToLHSDOFMap(lhsMapping%totalNumberOfRows),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate equations row to LHS DOF map.",err,error,*999)
    DO dofIdx=1,lhsMapping%totalNumberOfRows
      !1-1 mapping for now
      rowIdx=dofIdx
      lhsMapping%lhsDOFToEquationsRowMap(dofIdx)=rowIdx
    ENDDO !dofIdx
    DO rowIdx=1,lhsMapping%totalNumberOfRows
      !1-1 mapping for now
      dofIdx=rowIdx
      lhsMapping%equationsRowToLHSDOFMap(rowIdx)=dofIdx
    ENDDO !rowIdx
   
    EXITS("EquationsMapping_VectorCreateStart")
    RETURN
999 CALL EquationsMapping_VectorFinalise(vectorMapping,err,error,*998)
998 ERRORSEXITS("EquationsMapping_VectorCreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMapping_VectorCreateStart

  !
  !================================================================================================================================
  !

  !>Destroy an vector equations mapping.
  SUBROUTINE EquationsMapping_VectorDestroy(vectorMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer the vector equations mapping to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsMapping_VectorDestroy",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector equations mapping is not associated.",err,error,*999)
    
    CALL EquationsMapping_VectorFinalise(vectorMapping,err,error,*999)
        
    EXITS("EquationsMapping_VectorDestroy")
    RETURN
999 ERRORSEXITS("EquationsMapping_VectorDestroy",err,error)    
    RETURN 1
   
  END SUBROUTINE EquationsMapping_VectorDestroy

  !
  !================================================================================================================================
  !

  !>Finalises the equations mapping and deallocates all memory.
  SUBROUTINE EquationsMapping_VectorFinalise(vectorMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsMapping_VectorFinalise",err,error,*999)

    IF(ASSOCIATED(vectorMapping)) THEN
       !Row dofs mappings are linked to the field mapping therefore do not deallocate here
       !NULLIFY(vectorMapping%rowDofsMapping)
       CALL EquationsMappingVector_LHSMappingFinalise(vectorMapping%lhsMapping,err,error,*999)
       CALL EquationsMappingVector_DynamicMappingFinalise(vectorMapping%dynamicMapping,err,error,*999)
       CALL EquationsMappingVector_LinearMappingFinalise(vectorMapping%linearMapping,err,error,*999)
       CALL EquationsMappingVector_NonlinearMappingFinalise(vectorMapping%nonlinearMapping,err,error,*999)
       CALL EquationsMappingVector_SourcesMappingFinalise(vectorMapping%sourcesMapping,err,error,*999)      
       CALL EquationsMappingVector_RHSMappingFinalise(vectorMapping%rhsMapping,err,error,*999)      
       CALL EquationsMappingVector_CreateValuesCacheFinalise(vectorMapping%createValuesCache,err,error,*999)
       DEALLOCATE(vectorMapping)
    ENDIF
       
    EXITS("EquationsMapping_VectorFinalise")
    RETURN
999 ERRORSEXITS("EquationsMapping_VectorFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMapping_VectorFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the vector equations mapping.
  SUBROUTINE EquationsMapping_VectorInitialise(vectorEquations,err,error,*)

    !Argument variables
    TYPE(EquationsVectorType), POINTER :: vectorEquations !<A pointer to the vector equations to initialise the vector equations mapping for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("EquationsMapping_VectorInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(vectorEquations)) CALL FlagError("Vector equations is not associated.",err,error,*998)
    IF(ASSOCIATED(vectorEquations%vectorMapping)) CALL FlagError("Vector equations mapping is already associated.",err,error,*998)
     
    ALLOCATE(vectorEquations%vectorMapping,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate vector equations vector mapping.",err,error,*999)
    vectorEquations%vectorMapping%vectorEquations=>vectorEquations
    vectorEquations%vectorMapping%vectorMappingFinished=.FALSE.
    NULLIFY(vectorEquations%vectorMapping%vectorMatrices)
    NULLIFY(vectorEquations%vectorMapping%lhsMapping)
    NULLIFY(vectorEquations%vectorMapping%dynamicMapping)
    NULLIFY(vectorEquations%vectorMapping%linearMapping)
    NULLIFY(vectorEquations%vectorMapping%nonlinearMapping)
    NULLIFY(vectorEquations%vectorMapping%sourcesMapping)
    NULLIFY(vectorEquations%vectorMapping%rhsMapping)
    NULLIFY(vectorEquations%vectorMapping%createValuesCache)
    
    CALL EquationsMappingVector_LHSMappingInitialise(vectorEquations%vectorMapping,err,error,*999)    
    CALL EquationsMappingVector_CreateValuesCacheInitialise(vectorEquations%vectorMapping,err,error,*999)        
       
    EXITS("EquationsMapping_VectorInitialise")
    RETURN
999 CALL EquationsMapping_VectorFinalise(vectorEquations%vectorMapping,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMapping_VectorInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMapping_VectorInitialise

  !
  !================================================================================================================================
  !

  !>Finalises the scalar equations mapping dot product mapping and deallocates all memory
  SUBROUTINE EquationsMappingDotProduct_Finalise(dotProductMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingDotProductType) :: dotProductMapping !<A pointer to the dot product mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingDotProduct_Finalise",err,error,*999)

    dotProductMapping%dotProductNumber=0
    NULLIFY(dotProductMapping%dotProductVariables(1)%ptr)
    NULLIFY(dotProductMapping%dotProductVariables(2)%ptr)    
       
    EXITS("EquationsMappingDotProduct_Finalise")
    RETURN
999 ERRORSEXITS("EquationsMappingDotProduct_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingDotProduct_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises the scalar equations mapping dot product mapping
  SUBROUTINE EquationsMappingDotProduct_Initialise(dotProductMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingDotProductType) :: dotProductMapping !<The dot product mapping to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsMappingDotProduct_Initialise",err,error,*999)

    dotProductMapping%dotProductNumber=0
    NULLIFY(dotProductMapping%dotProductVariables(1)%ptr)
    NULLIFY(dotProductMapping%dotProductVariables(2)%ptr)
    dotProductMapping%dotProductCoefficient=1.0_DP
    
    EXITS("EquationsMappingDotProduct_Initialise")
    RETURN
999 ERRORS("EquationsMappingDotProduct_Initialise",err,error)
    EXITS("EquationsMappingDotProduct_Initialise")
    RETURN 1
    
  END SUBROUTINE EquationsMappingDotProduct_Initialise

  !
  !================================================================================================================================
  !

  !>Finalise an equations matrix to variable maps and deallocate all memory.
  SUBROUTINE EquationsMappingEMToVMap_Finalise(equationsMatrixToVarMap,err,error,*)

    !Argument variables
    TYPE(EquationsMatrixToVarMapType), POINTER :: equationsMatrixToVarMap !<A pointer to the equations matrix to variable map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsMappingEMToVMap_Finalise",err,error,*999)

    IF(ASSOCIATED(equationsMatrixToVarMap)) THEN
      IF(ASSOCIATED(equationsMatrixToVarMap%columnToDOFMap)) DEALLOCATE(equationsMatrixToVarMap%columnToDOFMap)
      DEALLOCATE(equationsMatrixToVarMap)
    ENDIF
    
    EXITS("EquationsMappingEMToVMap_Finalise")
    RETURN
999 ERRORS("EquationsMappingEMToVMap_Finalise",err,error)    
    EXITS("EquationsMappingEMToVMap_Finalise")    
    RETURN 1
   
  END SUBROUTINE EquationsMappingEMToVMap_Finalise

  !
  !================================================================================================================================
  !

  !>Initialise an equations matrix to variable maps.
  SUBROUTINE EquationsMappingEMToVMap_Initialise(equationsMatrixToVarMap,err,error,*)

    !Argument variables
    TYPE(EquationsMatrixToVarMapType), POINTER :: equationsMatrixToVarMap !<A pointer to the equations matrix to variable map to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("EquationsMappingEMToVMap_Initialise",err,error,*998)

    IF(ASSOCIATED(equationsMatrixToVarMap)) &
      & CALL FlagError("Equations matrix to variable map is already associated.",err,error,*998)

    ALLOCATE(equationsMatrixToVarMap,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate equations matrix to variable map.",err,error,*999)
    equationsMatrixToVarMap%matrixNumber=0
    NULLIFY(equationsMatrixToVarMap%equationsMatrix)
    equationsMatrixToVarMap%variableIndex=0
    equationsMatrixToVarMap%variableType=0
    NULLIFY(equationsMatrixToVarMap%variable)
    equationsMatrixToVarMap%numberOfColumns=0
    equationsMatrixToVarMap%matrixCoefficient=1.0_DP !Matrices in an equation set are added by default
    NULLIFY(equationsMatrixToVarMap%columnToDOFMap)
    NULLIFY(equationsMatrixToVarMap%columnDOFSMapping)
    
    EXITS("EquationsMappingEMToVMap_Initialise")
    RETURN
999 CALL EquationsMappingEMToVMap_Finalise(equationsMatrixToVarMap,dummyErr,dummyError,*998)
998 ERRORS("EquationsMappingEMToVMap_Initialise",err,error)    
    EXITS("EquationsMappingEMToVMap_Initialise")    
    RETURN 1
   
  END SUBROUTINE EquationsMappingEMToVMap_Initialise

  !
  !================================================================================================================================
  !

  !>Finalises the scalar equations mapping function mapping and deallocates all memory
  SUBROUTINE EquationsMappingFunction_Finalise(functionMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingFunctionType) :: functionMapping !<A pointer to the function mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingFunction_Finalise",err,error,*999)

    functionMapping%functionNumber=0
       
    EXITS("EquationsMappingFunction_Finalise")
    RETURN
999 ERRORSEXITS("EquationsMappingFunction_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingFunction_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises the scalar equations mapping function mapping
  SUBROUTINE EquationsMappingFunction_Initialise(functionMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingFunctionType) :: functionMapping !<The function mapping to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsMappingFunction_Initialise",err,error,*999)

    functionMapping%functionNumber=0
    
    EXITS("EquationsMappingFunction_Initialise")
    RETURN
999 ERRORS("EquationsMappingFunction_Initialise",err,error)
    EXITS("EquationsMappingFunction_Initialise")
    RETURN 1
    
  END SUBROUTINE EquationsMappingFunction_Initialise

  !
  !================================================================================================================================
  !

  !>Finalises a variable to equations Jacobian map and deallocates all memory.
  SUBROUTINE EquationsMappingJMToVMap_Finalise(jacobianMatrixToVarMap,err,error,*)

    !Argument variables
    TYPE(JacobianMatrixToVarMapType), POINTER :: jacobianMatrixToVarMap !<A pointer to the equations Jacobian to variable map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsMappingJMToVMap_Finalise",err,error,*999)

    IF(ASSOCIATED(jacobianMatrixToVarMap)) THEN
      IF(ASSOCIATED(jacobianMatrixToVarMap%equationsColumnToDOFVariableMap)) &
        & DEALLOCATE(jacobianMatrixToVarMap%equationsColumnToDOFVariableMap)
      DEALLOCATE(jacobianMatrixToVarMap)
    ENDIF
    
    EXITS("EquationsMappingJMToVMap_Finalise")
    RETURN
999 ERRORS("EquationsMappingJMToVMap_Finalise",err,error)    
    EXITS("EquationsMappingJMToVMap_Finalise")    
    RETURN 1
   
  END SUBROUTINE EquationsMappingJMToVMap_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises a variable to equations Jacobian map.
  SUBROUTINE EquationsMappingJMToVMap_Initialise(jacobianMatrixToVarMap,err,error,*)

    !Argument variables
    TYPE(JacobianMatrixToVarMapType), POINTER :: jacobianMatrixToVarMap !<A pointer to the equations Jacobian to variable map to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("EquationsMappingJMToVMap_Initialise",err,error,*998)

    IF(ASSOCIATED(jacobianMatrixToVarMap)) &
      & CALL FlagError("The Jacobian matrix to variable map is already associated.",err,error,*998)

    ALLOCATE(jacobianMatrixToVarMap,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate the Jacobian matrix to variable map.",err,error,*999)
    jacobianMatrixToVarMap%jacobianNumber=0
    jacobianMatrixToVarMap%variableType=0
    NULLIFY(jacobianMatrixToVarMap%variable)
    NULLIFY(jacobianMatrixToVarMap%jacobian)
    jacobianMatrixToVarMap%numberOfColumns=0
    jacobianMatrixToVarMap%jacobianCoefficient=0.0_DP
    NULLIFY(jacobianMatrixToVarMap%equationsColumnToDOFVariableMap)
    NULLIFY(jacobianMatrixToVarMap%columnDOFSMapping)    
    
    EXITS("EquationsMappingJMToVMap_Initialise")
    RETURN
999 CALL EquationsMappingJMToVMap_Finalise(jacobianMatrixToVarMap,dummyErr,dummyError,*998)
998 ERRORS("EquationsMappingJMToVMap_Initialise",err,error)    
    EXITS("EquationsMappingJMToVMap_Initialise")    
    RETURN 1
   
  END SUBROUTINE EquationsMappingJMToVMap_Initialise

  !
  !================================================================================================================================
  !

  !>Finalises a residual mapping and deallocates all memory
  SUBROUTINE EquationsMappingNonlinear_ResidualMappingFinalise(residualMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping !<A pointer to the residual mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: jacobianIdx,variableIdx
 
    ENTERS("EquationsMappingNonlinear_ResidualMappingFinalise",err,error,*999)

    IF(ASSOCIATED(residualMapping)) THEN
      IF(ALLOCATED(residualMapping%variableTypes)) DEALLOCATE(residualMapping%variableTypes)
      IF(ALLOCATED(residualMapping%variableTypesMap)) DEALLOCATE(residualMapping%variableTypesMap)
      IF(ALLOCATED(residualMapping%variables)) DEALLOCATE(residualMapping%variables)
      IF(ALLOCATED(residualMapping%varToJacobianMatrixMaps)) THEN
        DO variableIdx=1,SIZE(residualMapping%varToJacobianMatrixMaps,1)
          CALL EquationsMappingVToJMMap_Finalise(residualMapping%varToJacobianMatrixMaps(variableIdx)%ptr,err,error,*999)
        ENDDO !variableIdx
        DEALLOCATE(residualMapping%varToJacobianMatrixMaps)
      ENDIF
      IF(ALLOCATED(residualMapping%jacobianMatrixToVarMaps)) THEN
        DO jacobianIdx=1,SIZE(residualMapping%jacobianMatrixToVarMaps,1)
          CALL EquationsMappingJMToVMap_Finalise(residualMapping%jacobianMatrixToVarMaps(jacobianIdx)%ptr,err,error,*999)
        ENDDO !jacobianIdx
        DEALLOCATE(residualMapping%jacobianMatrixToVarMaps)
      ENDIF
      IF(ASSOCIATED(residualMapping%equationsRowToResidualDOFMap)) DEALLOCATE(residualMapping%equationsRowToResidualDOFMap)
      DEALLOCATE(residualMapping)
    ENDIF
       
    EXITS("EquationsMappingNonlinear_ResidualMappingFinalise")
    RETURN
999 ERRORSEXITS("EquationsMappingNonlinear_ResidualMappingFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingNonlinear_ResidualMappingFinalise

  !
  !================================================================================================================================
  !

  !>Allocates and initialises a residual mapping for a nonlinear mapping.
  SUBROUTINE EquationsMappingNonlinear_ResidualMappingInitialise(nonlinearMapping,residualIdx,err,error,*)

    !Argument variables
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping !<A pointer to the nonlinear mapping to initialise the residual mapping for.
    INTEGER(INTG), INTENT(IN) :: residualIdx !<The residual index to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError,localError
 
    ENTERS("EquationsMappingNonlinear_ResidualMappingInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(nonlinearMapping)) CALL FlagError("Nonlinear mapping is not associated.",err,error,*998)
    IF(residualIdx<1.OR.residualIdx>nonlinearMapping%numberOfResiduals) THEN
      localError="The specified residual index of "//TRIM(NumberToVString(residualIdx,"*",err,error))// &
        & " is invalid. The residual index should be >= 1 and <= "// &
        & TRIM(NumberToVString(nonlinearMapping%numberOfResiduals,"*",err,error))//"."
      CALL FlagError(localError,err,error,*998)      
    ENDIF
    IF(.NOT.ALLOCATED(nonlinearMapping%residuals)) &
      & CALL FlagError("The residuals array is not allocated for the nonlinear mapping.",err,error,*998)
    IF(ASSOCIATED(nonlinearMapping%residuals(residualIdx)%ptr)) THEN
      localError="The residual mapping is already associated for residual index "// &
        & TRIM(NumberToVString(residualIdx,"*",err,error))//"."
      CALL FlagError(localError,err,error,*998)
    ENDIF

    ALLOCATE(nonlinearMapping%residuals(residualIdx)%ptr,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate residual mapping.",err,error,*999)
    nonlinearMapping%residuals(residualIdx)%ptr%nonlinearMapping=>nonlinearMapping
    nonlinearMapping%residuals(residualIdx)%ptr%residualNumber=residualIdx
    nonlinearMapping%residuals(residualIdx)%ptr%numberOfVariables=0
    nonlinearMapping%residuals(residualIdx)%ptr%numberOfJacobianMatrices=0
    nonlinearMapping%residuals(residualIdx)%ptr%residualCoefficient=1.0_DP
    NULLIFY(nonlinearMapping%residuals(residualIdx)%ptr%equationsRowToResidualDOFMap)
        
    EXITS("EquationsMappingNonlinear_ResidualMappingInitialise")
    RETURN
999 CALL EquationsMappingNonlinear_ResidualMappingFinalise(nonlinearMapping%residuals(residualIdx)%ptr,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMappingNonlinear_ResidualMappingInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingNonlinear_ResidualMappingInitialise

  !
  !================================================================================================================================
  !

  !>Finalises the scalar equations mapping norm mapping and deallocates all memory
  SUBROUTINE EquationsMappingNormMapping_Finalise(normMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingNormType) :: normMapping !<A pointer to the norm mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingNormMapping_Finalise",err,error,*999)

    normMapping%normNumber=0
    NULLIFY(normMapping%normVariable)
       
    EXITS("EquationsMappingNormMapping_Finalise")
    RETURN
999 ERRORSEXITS("EquationsMappingNormMapping_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingNormMapping_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises the scalar equations mapping norm mapping
  SUBROUTINE EquationsMappingNormMapping_Initialise(normMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingNormType) :: normMapping !<The norm mapping to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsMappingNormMapping_Initialise",err,error,*999)

    normMapping%normNumber=0
    NULLIFY(normMapping%normVariable)
    normMapping%normCoefficient=1.0_DP
    
    EXITS("EquationsMappingNormMapping_Initialise")
    RETURN
999 ERRORSEXITS("EquationsMappingNormMapping_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingNormMapping_Initialise

  !
  !================================================================================================================================
  !

  !>Finalises the scalar equations mapping quadratic form mapping and deallocates all memory
  SUBROUTINE EquationsMappingQuadratic_Finalise(quadraticMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingQuadraticType) :: quadraticMapping !<A pointer to the quadratic mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingQuadratic_Finalise",err,error,*999)

    quadraticMapping%quadraticNumber=0
    NULLIFY(quadraticMapping%quadraticVariables(1)%ptr)
    NULLIFY(quadraticMapping%quadraticVariables(2)%ptr)
    quadraticMapping%quadraticCoefficient=1.0_DP
       
    EXITS("EquationsMappingQuadratic_Finalise")
    RETURN
999 ERRORSEXITS("EquationsMappingQuadratic_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingQuadratic_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises the scalar equations mapping quadratic form mapping
  SUBROUTINE EquationsMappingQuadratic_Initialise(quadraticMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingQuadraticType) :: quadraticMapping !<The quadratic form mapping to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsMappingQuadratic_Initialise",err,error,*999)

    quadraticMapping%quadraticNumber=0
    NULLIFY(quadraticMapping%quadraticVariables(1)%ptr)
    NULLIFY(quadraticMapping%quadraticVariables(2)%ptr)
    quadraticMapping%quadraticCoefficient=1.0_DP
    
    EXITS("EquationsMappingQuadratic_Initialise")
    RETURN
999 ERRORSEXITS("EquationsMappingQuadratic_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingQuadratic_Initialise

  !
  !================================================================================================================================
  !

  !>Finalises the scalar equations mapping dot product mappings and deallocates all memory
  SUBROUTINE EquationsMappingScalar_DotProductsFinalise(dotProductMappings,err,error,*)

    !Argument variables
    TYPE(EquationsMappingDotProductsType), POINTER :: dotProductMappings !<A pointer to the dot product mappings to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dotProductIdx

    ENTERS("EquationsMappingScalar_DotProductsFinalise",err,error,*999)

    IF(ASSOCIATED(dotProductMappings)) THEN
      IF(ALLOCATED(dotProductMappings%dotProducts)) THEN
        DO dotProductIdx=1,SIZE(dotProductMappings%dotProducts,1)
          CALL EquationsMappingDotProduct_Finalise(dotProductMappings%dotProducts(dotProductIdx),err,error,*999)
        ENDDO !dotProductIdx
        DEALLOCATE(dotProductMappings%dotProducts)
      ENDIF
      DEALLOCATE(dotProductMappings)
    ENDIF
       
    EXITS("EquationsMappingScalar_DotProductsFinalise")
    RETURN
999 ERRORS("EquationsMappingScalar_DotProductsFinalise",err,error)
    EXITS("EquationsMappingScalar_DotProductsFinalise")
    RETURN 1
    
  END SUBROUTINE EquationsMappingScalar_DotProductsFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the scalar equations mapping dot product mappings
  SUBROUTINE EquationsMappingScalar_DotProductsInitialise(scalarMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingScalarType), POINTER :: scalarMapping !<A pointer to the scalar equations mapping to initialise the dot product mappings for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("EquationsMappingScalar_DotProductsInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(scalarMapping)) CALL FlagError("Scalar equations mapping is not associated.",err,error,*998)
    IF(ASSOCIATED(scalarMapping%dotProductMappings)) &
      & CALL FlagError("Scalar equations mapping dot product mappings is already associated.",err,error,*998)
    
    ALLOCATE(scalarMapping%dotProductMappings,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate scalar equations mapping dot product mappings.",err,error,*999)
    scalarMapping%dotProductMappings%scalarMapping=>scalarMapping
    scalarMapping%dotProductMappings%numberOfDotProducts=0
    
    EXITS("EquationsMappingScalar_DotProductsInitialise")
    RETURN
999 CALL EquationsMappingScalar_DotProductsFinalise(scalarMapping%dotProductMappings,dummyErr,dummyError,*998)
998 ERRORS("EquationsMappingScalar_DotProductsInitialise",err,error)
    EXITS("EquationsMappingScalar_DotProductsInitialise")
    RETURN 1
    
  END SUBROUTINE EquationsMappingScalar_DotProductsInitialise

  !
  !================================================================================================================================
  !

  !>Finalises the scalar equations mapping function mappings and deallocates all memory
  SUBROUTINE EquationsMappingScalar_FunctionsFinalise(functionMappings,err,error,*)

    !Argument variables
    TYPE(EquationsMappingFunctionsType), POINTER :: functionMappings !<A pointer to the function mappings to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: functionIdx

    ENTERS("EquationsMappingScalar_FunctionsFinalise",err,error,*999)

    IF(ASSOCIATED(functionMappings)) THEN
      IF(ALLOCATED(functionMappings%functions)) THEN
        DO functionIdx=1,SIZE(functionMappings%functions,1)
          CALL EquationsMappingFunction_Finalise(functionMappings%functions(functionIdx),err,error,*999)
        ENDDO !functionIdx
        DEALLOCATE(functionMappings%functions)
      ENDIF
      DEALLOCATE(functionMappings)
    ENDIF
       
    EXITS("EquationsMappingScalar_FunctionsFinalise")
    RETURN
999 ERRORS("EquationsMappingScalar_FunctionsFinalise",err,error)
    EXITS("EquationsMappingScalar_FunctionsFinalise")
    RETURN 1
    
  END SUBROUTINE EquationsMappingScalar_FunctionsFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the scalar equations mapping function mappings
  SUBROUTINE EquationsMappingScalar_FunctionsInitialise(scalarMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingScalarType), POINTER :: scalarMapping !<A pointer to the scalar equations mapping to initialise the function mappings for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("EquationsMappingScalar_FunctionsInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(scalarMapping)) CALL FlagError("Scalar equations mapping is not associated.",err,error,*998)
    IF(ASSOCIATED(scalarMapping%functionMappings)) &
      & CALL FlagError("Scalar equations mapping function mappings is already associated.",err,error,*998)
    
    ALLOCATE(scalarMapping%functionMappings,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate scalar equations mapping function mappings.",err,error,*999)
    scalarMapping%functionMappings%scalarMapping=>scalarMapping
    scalarMapping%functionMappings%numberOfFunctions=0
    
    EXITS("EquationsMappingScalar_FunctionsInitialise")
    RETURN
999 CALL EquationsMappingScalar_FunctionsFinalise(scalarMapping%functionMappings,dummyErr,dummyError,*998)
998 ERRORS("EquationsMappingScalar_FunctionsInitialise",err,error)
    EXITS("EquationsMappingScalar_FunctionsInitialise")
    RETURN 1
    
  END SUBROUTINE EquationsMappingScalar_FunctionsInitialise

  !
  !================================================================================================================================
  !

  !>Finalises the scalar equations mapping norm mappings and deallocates all memory
  SUBROUTINE EquationsMappingScalar_NormsFinalise(normMappings,err,error,*)

    !Argument variables
    TYPE(EquationsMappingNormsType), POINTER :: normMappings !<A pointer to the norm mappings to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: normIdx

    ENTERS("EquationsMappingScalar_NormsFinalise",err,error,*999)

    IF(ASSOCIATED(normMappings)) THEN
      IF(ALLOCATED(normMappings%norms)) THEN
        DO normIdx=1,SIZE(normMappings%norms,1)
          CALL EquationsMappingNormMapping_Finalise(normMappings%norms(normIdx),err,error,*999)
        ENDDO !normIdx
        DEALLOCATE(normMappings%norms)
      ENDIF
      DEALLOCATE(normMappings)
    ENDIF
       
    EXITS("EquationsMappingScalar_NormsFinalise")
    RETURN
999 ERRORSEXITS("EquationsMappingScalar_NormsFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingScalar_NormsFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the scalar equations mapping norm mappings
  SUBROUTINE EquationsMappingScalar_NormsInitialise(scalarMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingScalarType), POINTER :: scalarMapping !<A pointer to the scalar equations mapping to initialise the norm mappings for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("EquationsMappingScalar_NormsInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(scalarMapping)) CALL FlagError("Scalar equations mapping is not associated.",err,error,*998)
    IF(ASSOCIATED(scalarMapping%normMappings)) &
      & CALL FlagError("Scalar equations mapping norm mappings is already associated.",err,error,*998)
    
    ALLOCATE(scalarMapping%normMappings,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate scalar equations mapping norm mappings.",err,error,*999)
    scalarMapping%normMappings%scalarMapping=>scalarMapping
    scalarMapping%normMappings%numberOfNorms=0
    
    EXITS("EquationsMappingScalar_NormsInitialise")
    RETURN
999 CALL EquationsMappingScalar_NormsFinalise(scalarMapping%normMappings,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMappingScalar_NormsInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingScalar_NormsInitialise

  !
  !================================================================================================================================
  !

  !>Finalises the scalar equations mapping quadratic form mappings and deallocates all memory
  SUBROUTINE EquationsMappingScalar_QuadraticsFinalise(quadraticMappings,err,error,*)

    !Argument variables
    TYPE(EquationsMappingQuadraticsType), POINTER :: quadraticMappings !<A pointer to the quadratic forms mappings to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: quadraticIdx

    ENTERS("EquationsMappingScalar_QuadraticsFinalise",err,error,*999)

    IF(ASSOCIATED(quadraticMappings)) THEN
      IF(ALLOCATED(quadraticMappings%quadratics)) THEN
        DO quadraticIdx=1,SIZE(quadraticMappings%quadratics,1)
          CALL EquationsMappingQuadratic_Finalise(quadraticMappings%quadratics(quadraticIdx),err,error,*999)
        ENDDO !quadraticIdx
        DEALLOCATE(quadraticMappings%quadratics)
      ENDIF
      DEALLOCATE(quadraticMappings)
    ENDIF
       
    EXITS("EquationsMappingScalar_QuadraticsFinalise")
    RETURN
999 ERRORSEXITS("EquationsMappingScalar_QuadraticsFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingScalar_QuadraticsFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the scalar equations mapping quadratic form mappings
  SUBROUTINE EquationsMappingScalar_QuadraticsInitialise(scalarMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingScalarType), POINTER :: scalarMapping !<A pointer to the scalar equations mapping to initialise the quadratic form mappings for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("EquationsMappingScalar_QuadraticsInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(scalarMapping)) CALL FlagError("Scalar equations mapping is not associated.",err,error,*998)
    IF(ASSOCIATED(scalarMapping%quadraticMappings)) &
      & CALL FlagError("Scalar equations mapping quadratic mappings is already associated.",err,error,*998)
    
    ALLOCATE(scalarMapping%quadraticMappings,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate scalar equations mapping quadratic mappings.",err,error,*999)
    scalarMapping%quadraticMappings%scalarMapping=>scalarMapping
    scalarMapping%quadraticMappings%numberOfQuadratics=0
    
    EXITS("EquationsMappingScalar_QuadraticsInitialise")
    RETURN
999 CALL EquationsMappingScalar_QuadraticsFinalise(scalarMapping%quadraticMappings,dummyErr,dummyError,*998)
998 ERRORS("EquationsMappingScalar_QuadraticsInitialise",err,error)
    EXITS("EquationsMappingScalar_QuadraticsInitialise")
    RETURN 1
    
  END SUBROUTINE EquationsMappingScalar_QuadraticsInitialise

  !
  !================================================================================================================================
  !

  !>Finalises the equations mapping source mapping and deallocates all memory
  SUBROUTINE EquationsMappingSources_SourceMappingFinalise(sourceMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingSourceType), POINTER :: sourceMapping !<A pointer to the SOURCE mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingSources_SourceMappingFinalise",err,error,*999)

    IF(ASSOCIATED(sourceMapping)) THEN
      IF(ASSOCIATED(sourceMapping%sourceDOFToEquationsRowMap)) DEALLOCATE(sourceMapping%sourceDOFToEquationsRowMap)
      IF(ASSOCIATED(sourceMapping%equationsRowToSourceDOFMap)) DEALLOCATE(sourceMapping%equationsRowToSourceDOFMap)
      DEALLOCATE(sourceMapping)
    ENDIF
       
    EXITS("EquationsMappingSources_SourceMappingFinalise")
    RETURN
999 ERRORSEXITS("EquationsMappingSources_SourceMappingFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingSources_SourceMappingFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the equations mapping source mapping
  SUBROUTINE EquationsMappingSources_SourceMappingInitialise(sourcesMapping,sourceIdx,err,error,*)

    !Argument variables
    TYPE(EquationsMappingSourcesType), POINTER :: sourcesMapping !<A pointer to the sources mapping to initialise the source mapping for
    INTEGER(INTG), INTENT(IN) :: sourceIdx !<The source index to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError,localError

    ENTERS("EquationsMappingSources_SourceMappingInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(sourcesMapping)) CALL FlagError("Sources mapping is not associated.",err,error,*998)
    IF(sourceIdx<1.OR.sourceIdx>sourcesMapping%numberOfSources) THEN
      localError="The specified source index of "//TRIM(NumberToVString(sourceIdx,"*",err,error))// &
        & " is invalid. The source index should be >= 1 and <= "// &
        & TRIM(NumberToVString(sourcesMapping%numberOfSources,"*",err,error))//"."
      CALL FlagError(localError,err,error,*998)      
    ENDIF
    IF(.NOT.ALLOCATED(sourcesMapping%sources)) &
      & CALL FlagError("The sources array is not allocated for the sources mapping.",err,error,*998)
    IF(ASSOCIATED(sourcesMapping%sources(sourceIdx)%ptr)) THEN
      localError="The source mapping is already associated for source index "// &
        & TRIM(NumberToVString(sourceIdx,"*",err,error))//"."
      CALL FlagError(localError,err,error,*998)
    ENDIF

    ALLOCATE(sourcesMapping%sources(sourceIdx)%ptr,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate source mapping.",err,error,*999)
    sourcesMapping%sources(sourceIdx)%ptr%sourcesMapping=>sourcesMapping
    sourcesMapping%sources(sourceIdx)%ptr%sourceNumber=sourceIdx
    sourcesMapping%sources(sourceIdx)%ptr%sourceVariableType=0
    NULLIFY(sourcesMapping%sources(sourceIdx)%ptr%sourceVariable)
    NULLIFY(sourcesMapping%sources(sourceIdx)%ptr%sourceVariableMapping)
    sourcesMapping%sources(sourceIdx)%ptr%sourceCoefficient=1.0_DP
    NULLIFY(sourcesMapping%sources(sourceIdx)%ptr%sourceDOFToEquationsRowMap)
    NULLIFY(sourcesMapping%sources(sourceIdx)%ptr%equationsRowToSourceDOFMap)
       
    EXITS("EquationsMappingSources_SourceMappingInitialise")
    RETURN
999 CALL EquationsMappingSources_SourceMappingFinalise(sourcesMapping%sources(sourceIdx)%ptr,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMappingSources_SourceMappingInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingSources_SourceMappingInitialise

  !
  !================================================================================================================================
  !

  !>Calculates the equations vector mapping.
  SUBROUTINE EquationsMappingVector_Calculate(vectorMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations vector mapping to calculate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: columnIdx,dofIdx,lhsVariableType,linearMatrixStart,matrixIdx,numberOfDependentDOFs, &
      & numberOfEquationsMatrices,numberOfGlobal,numberOfLHSDOFs,numberOfRows,numberOfGlobalRows,numberOfVariables, &
      & rowIdx,residualIdx,sourceIdx,totalNumberOfDependentDOFs,totalNumberOfLHSDOFs,totalNumberOfRows,variableIdx, &
      & variableType,variableTypeIdx
    INTEGER(INTG), ALLOCATABLE :: variableTypes(:)
    INTEGER(INTG), POINTER :: equationsRowToLHSDOFMap(:),equationsRowToResidualDOFMap(:),equationsRowToVariableDOFMap(:), &
      & lhsDOFToEquationsRowMap(:),linEquationsRowToVariableDOFMap(:,:)
    REAL(DP) :: dynamicMatrixCoefficient,residualCoefficient
    TYPE(DomainMappingType), POINTER :: columnDomainMapping,rowDomainMapping
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsMappingVectorCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping
    TYPE(EquationsMappingSourceType), POINTER :: sourceMapping
    TYPE(EquationsMappingSourcesType), POINTER :: sourcesMapping
    TYPE(EquationsMatrixToVarMapType), POINTER :: equationsMatrixToVarMap
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(FieldType), POINTER :: dependentField,sourceField
    TYPE(FieldVariableType), POINTER :: dependentVariable,dynamicVariable,jacobianVariable,lhsVariable,linearVariable, &
      & rhsVariable,rowVariable,sourceVariable
    TYPE(JacobianMatrixToVarMapType), POINTER :: jacobianMatrixToVarMap
    TYPE(ListType), POINTER :: variableList
    TYPE(VarToEquationsMatricesMapType), POINTER :: varToEquationsMatricesMap
    TYPE(VarToJacobianMatrixMapType), POINTER :: varToJacobianMatrixMap
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsMappingVector_Calculate",err,error,*999)

    NULLIFY(vectorEquations)
    CALL EquationsMappingVector_VectorEquationsGet(vectorMapping,vectorEquations,err,error,*999)
    NULLIFY(equations)
    CALL EquationsVector_EquationsGet(vectorEquations,equations,err,error,*999)
    NULLIFY(equationsSet)
    CALL Equations_EquationsSetGet(equations,equationsSet,err,error,*999)
    NULLIFY(createValuesCache)
    CALL EquationsMappingVector_CreateValuesCacheGet(vectorMapping,createValuesCache,err,error,*999)

    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)

    !Get the number of rows in the equations set from the LHS mapping
    NULLIFY(lhsMapping)
    CALL EquationsMappingVector_LHSMappingGet(vectorMapping,lhsMapping,err,error,*999)
    NULLIFY(lhsVariable)
    CALL EquationsMappingLHS_LHSVariableGet(lhsMapping,lhsVariable,err,error,*999)
    CALL FieldVariable_NumberOfDOFsGet(lhsVariable,numberOfLHSDOFs,err,error,*999)
    numberOfRows=numberOfLHSDOFs
    CALL FieldVariable_TotalNumberOfDOFsGet(lhsVariable,totalNumberOfLHSDOFs,err,error,*999)
    totalNumberOfRows=totalNumberOfLHSDOFs
    NULLIFY(rowDomainMapping)
    CALL EquationsMappingLHS_RowDOFsMappingGet(lhsMapping,rowDomainMapping,err,error,*999)

    !Calculate dynamic mappings
    IF(createValuesCache%dynamicVariableType/=0) THEN
      CALL EquationsMappingVector_DynamicMappingInitialise(vectorMapping,err,error,*999)
      NULLIFY(dynamicMapping)
      CALL EquationsMappingVector_DynamicMappingGet(vectorMapping,dynamicMapping,err,error,*999)
      dynamicMapping%numberOfDynamicMatrices=createValuesCache%numberOfDynamicMatrices
      dynamicMapping%stiffnessMatrixNumber=createValuesCache%dynamicStiffnessMatrixNumber
      dynamicMapping%dampingMatrixNumber=createValuesCache%dynamicDampingMatrixNumber
      dynamicMapping%massMatrixNumber=createValuesCache%dynamicMassMatrixNumber
      !Initialise the variable map
      NULLIFY(dependentVariable)
      CALL Field_VariableGet(dependentField,createValuesCache%dynamicVariableType,dependentVariable,err,error,*999)
      NULLIFY(columnDomainMapping)
      CALL FieldVariable_DomainMappingGet(dependentVariable,columnDomainMapping,err,error,*999)
      CALL FieldVariable_NumberOfDOFsGet(lhsVariable,numberOfDependentDOFs,err,error,*999)
      CALL FieldVariable_TotalNumberOfDOFsGet(dependentVariable,totalNumberOfDependentDOFs,err,error,*999)
      CALL DomainMapping_NumberOfGlobalGet(columnDomainMapping,numberOfGlobal,err,error,*999)
      dynamicMapping%dynamicVariableType=createValuesCache%dynamicVariableType
      dynamicMapping%dynamicVariable=>dependentVariable
      NULLIFY(dynamicMapping%varToEquationsMatricesMap)
      CALL EquationsMappingVToEMSMap_Initialise(dynamicMapping%varToEquationsMatricesMap,err,error,*999)
      variableIdx=1
      dynamicMapping%varToEquationsMatricesMap%variableIndex=variableIdx
      dynamicMapping%varToEquationsMatricesMap%variableType=createValuesCache%dynamicVariableType
      dynamicMapping%varToEquationsMatricesMap%variable=>dependentVariable
      dynamicMapping%varToEquationsMatricesMap%numberOfEquationsMatrices=createValuesCache%numberOfDynamicMatrices
      !Allocate and initialise the variable to equations matrices maps
      ALLOCATE(dynamicMapping%varToEquationsMatricesMap%equationsMatrixNumbers( &
        & dynamicMapping%varToEquationsMatricesMap%numberOfEquationsMatrices),STAT=err)      
      IF(err/=0) CALL FlagError("Could not allocate variable to equations matrices maps equations matrix numbers.",err,error,*999)
      ALLOCATE(dynamicMapping%varToEquationsMatricesMap%dofToColumnsMaps(dynamicMapping%varToEquationsMatricesMap% &
        & numberOfEquationsMatrices),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate variable to equations matrices maps DOF to columns map.",err,error,*999)
      dynamicMapping%varToEquationsMatricesMap%equationsMatrixNumbers=0
      IF(createValuesCache%dynamicStiffnessMatrixNumber/=0) &
        & dynamicMapping%varToEquationsMatricesMap%equationsMatrixNumbers(createValuesCache%dynamicStiffnessMatrixNumber)= &
        & createValuesCache%dynamicStiffnessMatrixNumber
      IF(createValuesCache%dynamicDampingMatrixNumber/=0) &
        & dynamicMapping%varToEquationsMatricesMap%equationsMatrixNumbers(createValuesCache%dynamicDampingMatrixNumber)= &
        & createValuesCache%dynamicDampingMatrixNumber
      IF(createValuesCache%dynamicMassMatrixNumber/=0) &
        & dynamicMapping%varToEquationsMatricesMap%equationsMatrixNumbers(createValuesCache%dynamicMassMatrixNumber)= &
        & createValuesCache%dynamicMassMatrixNumber
      DO matrixIdx=1,dynamicMapping%varToEquationsMatricesMap%numberOfEquationsMatrices
        ALLOCATE(dynamicMapping%varToEquationsMatricesMap%dofToColumnsMaps(matrixIdx)%columnDOF(totalNumberOfDependentDOFs), &
          & STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate variable DOF to columns map column DOF.",err,error,*999)
        DO dofIdx=1,totalNumberOfDependentDofs
          !1-1 mapping for now
          CALL DomainMapping_LocalToGlobalGet(columnDomainMapping,dofIdx,columnIdx,err,error,*999)
          dynamicMapping%varToEquationsMatricesMap%dofToColumnsMaps(matrixIdx)%columnDOF(dofIdx)=columnIdx
        ENDDO !dofIdx                          
      ENDDO !matrixIdx
      !Allocate and initialise the equations matrix to variable maps types
      ALLOCATE(dynamicMapping%equationsMatrixToVarMaps(dynamicMapping%numberOfDynamicMatrices),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate equations mapping equations matrix to variable maps.",err,error,*999)
      !Create the individual matrix maps and column maps
      DO matrixIdx=1,dynamicMapping%numberOfDynamicMatrices
        NULLIFY(dynamicMapping%equationsMatrixToVarMaps(matrixIdx)%ptr)
        CALL EquationsMappingEMToVMap_Initialise(dynamicMapping%equationsMatrixToVarMaps(matrixIdx)%ptr,err,error,*999)
        dynamicMapping%equationsMatrixToVarMaps(matrixIdx)%ptr%matrixNumber=matrixIdx
        dynamicMapping%equationsMatrixToVarMaps(matrixIdx)%ptr%variableType=createValuesCache%dynamicVariableType
        dynamicMapping%equationsMatrixToVarMaps(matrixIdx)%ptr%variable=>dependentVariable
        dynamicMapping%equationsMatrixToVarMaps(matrixIdx)%ptr%numberOfColumns=numberOfGlobal
        CALL EquationsMappingVectorCVC_DynamicMatrixCoefficientGet(createValuesCache,matrixIdx,dynamicMatrixCoefficient, &
          & err,error,*999)
        dynamicMapping%equationsMatrixToVarMaps(matrixIdx)%ptr%matrixCoefficient=dynamicMatrixCoefficient
        ALLOCATE(dynamicMapping%equationsMatrixToVarMaps(matrixIdx)%ptr%columnToDOFMap(numberOfGlobal),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate equation matrix to variable map column to dof map.",err,error,*999)
        dynamicMapping%equationsMatrixToVarMaps(matrixIdx)%ptr%columnToDOFMap=0
        DO dofIdx=1,totalNumberOfDependentDofs
          !1-1 mapping for now
          CALL DomainMapping_LocalToGlobalGet(columnDomainMapping,dofIdx,columnIdx,err,error,*999)
          dynamicMapping%equationsMatrixToVarMaps(matrixIdx)%ptr%columnToDOFMap(columnIdx)=dofIdx
        ENDDO !dofIdx
        dynamicMapping%equationsMatrixToVarMaps(matrixIdx)%ptr%columnDOFSMapping=>columnDomainMapping
      ENDDO !matrixIdx
      !Set up row <-> DOF mappings. For now just have 1-1 so check the variable against the LHS variable
      IF(numberOfDependentDOFs/=numberOfLHSDOFs) THEN
        localError="The dynamic dependent variable is not compatible with the LHS variable for 1-1 mapping. "// &
          & "The number of DOFs for the dynamic depdendent variable is "// &
          & TRIM(NumberToVString(numberOfDependentDOFs,"*",err,error))// &
          & " and the number of DOFs for the LHS variable is "// &
          & TRIM(NumberToVString(numberOfLHSDOFs,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      IF(totalNumberOfDependentDOFs/=totalNumberOfLHSDOFs) THEN
        localError="The dynamic dependent variable is not compatible with the LHS variable for 1-1 mapping. "// &
          & "The total number of DOFs for the dynamic depdendent variable is "// &
          & TRIM(NumberToVString(totalNumberOfDependentDOFs,"*",err,error))// &
          & " and the total number of DOFs for the LHS variable is "// &
          & TRIM(NumberToVString(totalNumberOfLHSDOFs,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      !Allocate the row mappings
      ALLOCATE(dynamicMapping%equationsRowToVariableDOFMaps(totalNumberOfRows),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate equations row to variable DOF maps.",err,error,*999)
      !Set up the row mappings
      DO rowIdx=1,totalNumberOfRows
        !1-1 mapping for now
        dynamicMapping%equationsRowToVariableDOFMaps(rowIdx)=rowIdx
      ENDDO !rowIdx
    ENDIF
    !Calculate linear mappings
    IF(createValuesCache%numberOfLinearMatrices>0) THEN                  
      CALL EquationsMappingVector_LinearMappingInitialise(vectorMapping,err,error,*999)
      NULLIFY(linearMapping)
      CALL EquationsMappingVector_LinearMappingGet(vectorMapping,linearMapping,err,error,*999)
      linearMapping%numberOfLinearMatrices=createValuesCache%numberOfLinearMatrices
      !Find the list of individual variables
      NULLIFY(variableList)
      CALL List_CreateStart(variableList,err,error,*999)
      CALL List_DataTypeSet(variableList,LIST_INTG_TYPE,err,error,*999)
      CALL List_CreateFinish(variableList,err,error,*999)
      DO matrixIdx=1,linearMapping%numberOfLinearMatrices
        CALL EquationsMappingVectorCVC_LinearMatrixVariableTypeGet(createValuesCache,matrixIdx,variableType,err,error,*999)
        CALL List_ItemAdd(variableList,variableType,err,error,*999)
      ENDDO !matrixIdx
      CALL List_RemoveDuplicates(variableList,err,error,*999)
      CALL List_DetachAndDestroy(variableList,numberOfVariables,variableTypes,err,error,*999)
      !Allocate and initialise the independent variables maps
      linearMapping%numberOfLinearVariables=numberOfVariables
      ALLOCATE(linearMapping%linearVariableTypes(numberOfVariables),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate linear mapping linear variable types.",err,error,*999)
      ALLOCATE(linearMapping%linearVariableTypesMap(FIELD_NUMBER_OF_VARIABLE_TYPES),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate linear mapping linear variable types map.",err,error,*999)
      ALLOCATE(linearMapping%linearVariables(numberOfVariables),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate linear mapping linear variables.",err,error,*999)
      linearMapping%linearVariableTypesMap=0
      DO variableIdx=1,numberOfVariables
        linearMapping%linearVariableTypes(variableIdx)=variableTypes(variableIdx)
        linearMapping%linearVariableTypesMap(variableTypes(variableIdx))=variableIdx
        NULLIFY(linearMapping%linearVariables(variableIdx)%ptr)
        CALL Field_VariableGet(dependentField,variableTypes(variableIdx),linearMapping%linearVariables(variableIdx)%ptr, &
          & err,error,*999)
      ENDDO !variableIdx
      IF(ALLOCATED(variableTypes)) DEALLOCATE(variableTypes)
      !Allocate and initialise the variable type maps
      ALLOCATE(linearMapping%varToEquationsMatricesMaps(numberOfVariables),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate equations mapping variable to equations map.",err,error,*999)
      DO variableIdx=1,numberOfVariables
        NULLIFY(linearMapping%varToEquationsMatricesMaps(variableIdx)%ptr)
        CALL EquationsMappingVToEMSMap_Initialise(linearMapping%varToEquationsMatricesMaps(variableIdx)%ptr,err,error,*999)
        linearMapping%varToEquationsMatricesMaps(variableIdx)%ptr%variableIndex=variableIdx
        linearMapping%varToEquationsMatricesMaps(variableIdx)%ptr%variableType=linearMapping%linearVariableTypes(variableIdx)
        linearMapping%varToEquationsMatricesMaps(variableIdx)%ptr%variable=>linearMapping%linearVariables(variableIdx)%ptr
      ENDDO !variableIdx
      !Calculate the number of variable type maps and initialise
      DO matrixIdx=1,linearMapping%numberOfLinearMatrices
        CALL EquationsMappingVectorCVC_LinearMatrixVariableTypeGet(createValuesCache,matrixIdx,variableType,err,error,*999)
        variableIdx=linearMapping%linearVariableTypesMap(variableType)
        linearMapping%varToEquationsMatricesMaps(variableIdx)%ptr%numberOfEquationsMatrices= &
          & linearMapping%varToEquationsMatricesMaps(variableIdx)%ptr%numberOfEquationsMatrices+1
      ENDDO !matrixIdx
      !Allocate and initialise the variable to equations matrices maps
      DO variableIdx=1,numberOfVariables      
        NULLIFY(dependentVariable)
        CALL EquationsMappingLinear_LinearVariableGet(linearMapping,variableIdx,dependentVariable,err,error,*999)
        NULLIFY(columnDomainMapping)
        CALL FieldVariable_DomainMappingGet(dependentVariable,columnDomainMapping,err,error,*999)
        CALL FieldVariable_TotalNumberOfDOFsGet(dependentVariable,totalNumberOfDependentDOFs,err,error,*999)
        IF(linearMapping%varToEquationsMatricesMaps(variableIdx)%ptr%numberOfEquationsMatrices>0) THEN
          ALLOCATE(linearMapping%varToEquationsMatricesMaps(variableIdx)%ptr%equationsMatrixNumbers( &
            & linearMapping%varToEquationsMatricesMaps(variableIdx)%ptr%numberOfEquationsMatrices),STAT=err)
          IF(err/=0) &
            & CALL FlagError("Could not allocate variable to equations matrices maps equations matrix numbers.",err,error,*999)
          ALLOCATE(linearMapping%varToEquationsMatricesMaps(variableIdx)%ptr%dofToColumnsMaps( &
            & linearMapping%varToEquationsMatricesMaps(variableIdx)%ptr%numberOfEquationsMatrices),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate variable to equations matrices maps DOF to columns map.",err,error,*999)
          linearMapping%varToEquationsMatricesMaps(variableIdx)%ptr%equationsMatrixNumbers=0
          linearMapping%varToEquationsMatricesMaps(variableIdx)%ptr%numberOfEquationsMatrices=0
          DO matrixIdx=1,linearMapping%numberOfLinearMatrices
            CALL EquationsMappingVectorCVC_LinearMatrixVariableTypeGet(createValuesCache,matrixIdx,variableType,err,error,*999)
            IF(linearMapping%linearVariableTypesMap(variableType)==variableIdx) THEN
              linearMapping%varToEquationsMatricesMaps(variableIdx)%ptr%numberOfEquationsMatrices= &
                & linearMapping%varToEquationsMatricesMaps(variableIdx)%ptr%numberOfEquationsMatrices+1
              linearMapping%varToEquationsMatricesMaps(variableIdx)%ptr%equationsMatrixNumbers( &
                & linearMapping%varToEquationsMatricesMaps(variableIdx)%ptr%numberOfEquationsMatrices)=matrixIdx
              ALLOCATE(linearMapping%varToEquationsMatricesMaps(variableIdx)%ptr%dofToColumnsMaps( &
                & linearMapping%varToEquationsMatricesMaps(variableIdx)%ptr% &
                & numberOfEquationsMatrices)%columnDOF(totalNumberOfDependentDOFs),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate variable DOF to columns map column DOF.",err,error,*999)
              DO dofIdx=1,totalNumberOfDependentDofs
                !1-1 mapping for now
                CALL DomainMapping_LocalToGlobalGet(columnDomainMapping,dofIdx,columnIdx,err,error,*999)
                linearMapping%varToEquationsMatricesMaps(variableIdx)%ptr%dofToColumnsMaps( &
                  & linearMapping%varToEquationsMatricesMaps(variableIdx)%ptr% &
                  & numberOfEquationsMatrices)%columnDOF(dofIdx)=columnIdx
              ENDDO !dofIdx
            ENDIF
          ENDDO !matrixIdx
        ENDIF
      ENDDO !variableIdx
      !Allocate and initialise the equations matrix to variable maps types
      ALLOCATE(linearMapping%equationsMatrixToVarMaps(linearMapping%numberOfLinearMatrices),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate equations mapping equations matrix to variable maps.",err,error,*999)
      !Create the individual matrix maps and column maps
      DO matrixIdx=1,linearMapping%numberOfLinearMatrices
        CALL EquationsMappingVectorCVC_LinearMatrixVariableTypeGet(createValuesCache,matrixIdx,variableType,err,error,*999)
        variableIdx=linearMapping%linearVariableTypesMap(variableType)
        NULLIFY(dependentVariable)
        CALL EquationsMappingLinear_LinearVariableGet(linearMapping,variableIdx,dependentVariable,err,error,*999)
        NULLIFY(columnDomainMapping)
        CALL FieldVariable_DomainMappingGet(dependentVariable,columnDomainMapping,err,error,*999)
        CALL DomainMapping_NumberOfGlobalGet(columnDomainMapping,numberOfGlobal,err,error,*999)
        NULLIFY(linearMapping%equationsMatrixToVarMaps(matrixIdx)%ptr)
        CALL EquationsMappingEMToVMap_Initialise(linearMapping%equationsMatrixToVarMaps(matrixIdx)%ptr,err,error,*999)
        linearMapping%equationsMatrixToVarMaps(matrixIdx)%ptr%matrixNumber=matrixIdx
        linearMapping%equationsMatrixToVarMaps(matrixIdx)%ptr%variableIndex=variableIdx
        linearMapping%equationsMatrixToVarMaps(matrixIdx)%ptr%variableType=variableType
        linearMapping%equationsMatrixToVarMaps(matrixIdx)%ptr%variable=>dependentVariable
        linearMapping%equationsMatrixToVarMaps(matrixIdx)%ptr%numberOfColumns=numberOfGlobal
        CALL EquationsMappingVectorCVC_LinearMatrixCoefficientGet(createValuesCache,matrixIdx,linearMapping% &
          & equationsMatrixToVarMaps(matrixIdx)%ptr%matrixCoefficient,err,error,*999)
        ALLOCATE(linearMapping%equationsMatrixToVarMaps(matrixIdx)%ptr%columnToDOFMap(numberOfGlobal),STAT=err)                  
        IF(err/=0) CALL FlagError("Could not allocate equation matrix to variable map column to dof map.",err,error,*999)
        linearMapping%equationsMatrixToVarMaps(matrixIdx)%ptr%columnToDOFMap=0
        DO dofIdx=1,totalNumberOfDependentDOFs
          !1-1 mapping for now
          CALL DomainMapping_LocalToGlobalGet(columnDomainMapping,dofIdx,columnIdx,err,error,*999)
          linearMapping%equationsMatrixToVarMaps(matrixIdx)%ptr%columnToDOFMap(columnIdx)=dofIdx
        ENDDO !dofIdx
        linearMapping%equationsMatrixToVarMaps(matrixIdx)%ptr%columnDOFSMapping=>columnDomainMapping
      ENDDO !matrixIdx
      !Set up row <-> DOF mappings. For now just have 1-1 so check the variable against the LHS variable
      DO variableIdx=1,numberOfVariables      
        NULLIFY(dependentVariable)
        CALL EquationsMappingLinear_LinearVariableGet(linearMapping,variableIdx,dependentVariable,err,error,*999)
        NULLIFY(columnDomainMapping)
        CALL FieldVariable_DomainMappingGet(dependentVariable,columnDomainMapping,err,error,*999)
        CALL FieldVariable_NumberOfDOFsGet(dependentVariable,numberOfDependentDOFs,err,error,*999)
        CALL FieldVariable_TotalNumberOfDOFsGet(dependentVariable,totalNumberOfDependentDOFs,err,error,*999)
        IF(numberOfDependentDOFs/=numberOfLHSDOFs) THEN
          localError="The linear dependent variable number "//TRIM(NumberToVString(variableIdx,"*",err,error))// &
            & "is not compatible with the LHS variable for 1-1 mapping. The number of DOFs for the linear depdendent "// &
            & "variable is "//TRIM(NumberToVString(numberOfDependentDOFs,"*",err,error))// &
            & " and the number of DOFs for the LHS variable is "// &
            & TRIM(NumberToVString(numberOfLHSDOFs,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        IF(totalNumberOfDependentDOFs/=totalNumberOfLHSDOFs) THEN
          localError="The linear dependent variable number "//TRIM(NumberToVString(variableIdx,"*",err,error))// &
            & "is not compatible with the LHS variable for 1-1 mapping. The total number of DOFs for the linear depdendent "// &
            & "variable is "//TRIM(NumberToVString(totalNumberOfDependentDOFs,"*",err,error))// &
            & " and the total number of DOFs for the LHS variable is "// &
            & TRIM(NumberToVString(totalNumberOfLHSDOFs,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDDO !variableIdx
      !Allocate the row mappings
      ALLOCATE(linearMapping%equationsRowToVariableDOFMaps(totalNumberOfRows,linearMapping%numberOfLinearVariables), &
        & STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate equations row to variable DOF maps.",err,error,*999)
      !Set up the row mappings
      DO variableIdx=1,linearMapping%numberOfLinearVariables
        DO rowIdx=1,totalNumberOfRows
          !1-1 mapping for now
          linearMapping%equationsRowToVariableDOFMaps(rowIdx,variableIdx)=rowIdx
        ENDDO !rowIdx
      ENDDO !variableIdx
    ENDIF
    !Calculate non-linear mappings
    IF(createValuesCache%numberOfResiduals>0) THEN
      CALL EquationsMappingVector_NonlinearMappingInitialise(vectorMapping,err,error,*999)
      NULLIFY(nonlinearMapping)
      CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)
      ALLOCATE(nonlinearMapping%residuals(createValuesCache%numberOfResiduals),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate residuals.",err,error,*999)
      nonlinearMapping%numberOfresiduals=createValuesCache%numberOfResiduals
      DO residualIdx=1,createValuesCache%numberOfResiduals
        CALL EquationsMappingVectorCVC_NumberOfResidualVariablesGet(createValuesCache,residualIdx,numberOfVariables,err,error,*999)
        IF(numberOfVariables<1) THEN
          localError="The number of variables of "//TRIM(NumberToVString(numberOfVariables,"*",err,error))// &
            & " mapped to residual number "//TRIM(NumberToVString(residualIdx,"*",err,error))// &
            & " is invalid. The number of variables should be >= 1."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        NULLIFY(nonlinearMapping%residuals(residualIdx)%ptr)
        CALL EquationsMappingNonlinear_ResidualMappingInitialise(nonlinearMapping,residualIdx,err,error,*999)
        NULLIFY(residualMapping)
        CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,residualIdx,residualMapping,err,error,*999)
        !Find the list of individual residual variables
        NULLIFY(variableList)
        CALL List_CreateStart(variableList,err,error,*999)
        CALL List_DataTypeSet(variableList,LIST_INTG_TYPE,err,error,*999)
        CALL List_CreateFinish(variableList,err,error,*999)
        DO variableIdx=1,numberOfVariables
          CALL EquationsMappingVectorCVC_ResidualVariableTypeGet(createValuesCache,variableIdx,residualIdx,variableType, &
            & err,error,*999)
          CALL List_ItemAdd(variableList,variableType,err,error,*999)
        ENDDO !matrixIdx
        CALL List_RemoveDuplicates(variableList,err,error,*999)
        CALL List_DetachAndDestroy(variableList,numberOfVariables,variableTypes,err,error,*999)        
        residualMapping%numberOfVariables=numberOfVariables
        ALLOCATE(residualMapping%variableTypes(numberOfVariables),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate residual mapping variable types.",err,error,*999)
        ALLOCATE(residualMapping%variableTypesMap(FIELD_NUMBER_OF_VARIABLE_TYPES),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate residual mapping variable types map.",err,error,*999)
        ALLOCATE(residualMapping%variables(numberOfVariables),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate residual mapping variables.",err,error,*999)
        residualMapping%variableTypesMap=0
        DO variableIdx=1,residualMapping%numberOfVariables
          residualMapping%variableTypes(variableIdx)=variableTypes(variableIdx)
          residualMapping%variableTypesMap(variableTypes(variableIdx))=variableIdx
          NULLIFY(residualMapping%variables(variableIdx)%ptr)
          CALL Field_VariableGet(dependentField,variableTypes(variableIdx),residualMapping%variables(variableIdx)%ptr, &
            & err,error,*999)
        ENDDO !variableIdx
        IF(ALLOCATED(variableTypes)) DEALLOCATE(variableTypes)
        residualMapping%numberOfJacobianMatrices=residualMapping%numberOfVariables
        ALLOCATE(residualMapping%varToJacobianMatrixMaps(residualMapping%numberOfVariables),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate variable to Jacobian matrix maps.",err,error,*999)
        ALLOCATE(residualMapping%jacobianMatrixToVarMaps(residualMapping%numberOfVariables),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate Jacobian matrix to variable maps.",err,error,*999)
        CALL EquationsMappingVectorCVC_ResidualCoefficientGet(createValuesCache,residualIdx,residualCoefficient,err,error,*999)
        DO variableIdx=1,residualMapping%numberOfVariables
          NULLIFY(residualMapping%varToJacobianMatrixMaps(variableIdx)%ptr)
          CALL EquationsMappingVToJMMap_Initialise(residualMapping%varToJacobianMatrixMaps(variableIdx)%ptr,err,error,*999)
          dependentVariable=>residualMapping%variables(variableIdx)%ptr
          CALL FieldVariable_NumberOfDOFsGet(dependentVariable,numberOfDependentDOFs,err,error,*999)
          CALL FieldVariable_TotalNumberOfDOFsGet(dependentVariable,totalNumberOfDependentDOFs,err,error,*999)
          NULLIFY(columnDomainMapping)
          CALL FieldVariable_DomainMappingGet(dependentVariable,columnDomainMapping,err,error,*999)
          CALL DomainMapping_NumberOfGlobalGet(columnDomainMapping,numberOfGlobal,err,error,*999)
          residualMapping%varToJacobianMatrixMaps(variableIdx)%ptr%jacobianNumber=variableIdx
          residualMapping%varToJacobianMatrixMaps(variableIdx)%ptr%variableType=residualMapping%variableTypes(variableIdx)
          residualMapping%varToJacobianMatrixMaps(variableIdx)%ptr%variable=>dependentVariable
          !Allocate and set dof to Jacobian columns map
          ALLOCATE(residualMapping%varToJacobianMatrixMaps(variableIdx)%ptr%dofToColumnsMap(totalNumberOfDependentDOFs),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate variable to Jacobian matrix map DOF to columns map.",err,error,*999)
          DO dofIdx=1,totalNumberOfDependentDOFs
            !1-1 mapping for now
            CALL DomainMapping_LocalToGlobalGet(columnDomainMapping,dofIdx,columnIdx,err,error,*999)
            residualMapping%varToJacobianMatrixMaps(variableIdx)%ptr%dofToColumnsMap(dofIdx)=columnIdx
          ENDDO !dofIdx
          !Set up row <-> DOF mappings. For now just have 1-1 so check the variable against the LHS variable
          IF(numberOfDependentDOFs/=numberOfLHSDOFs) THEN
            localError="The residual dependent variable number "//TRIM(NumberToVString(variableIdx,"*",err,error))// &
              & " for residual number "//TRIM(NumberToVString(residualIdx,"*",err,error))//" is not compatible with the "// &
              & "LHS variable for 1-1 mapping. The number of DOFs for the residual depdendent variable is "// &
              & TRIM(NumberToVString(numberOfDependentDOFs,"*",err,error))// &
              & " and the number of DOFs for the LHS variable is "// &
              & TRIM(NumberToVString(numberOfLHSDOFs,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
          IF(totalNumberOfDependentDOFs/=totalNumberOfLHSDOFs) THEN
            localError="The residual dependent variable number "//TRIM(NumberToVString(variableIdx,"*",err,error))// &
              & " for residual number "//TRIM(NumberToVString(residualIdx,"*",err,error))//" is not compatible with the "// &
              & "is not compatible with the LHS variable for 1-1 mapping. The total number of DOFs for the residual depdendent"// &
              & " variable is "//TRIM(NumberToVString(totalNumberOfDependentDOFs,"*",err,error))// &
              & " and the total number of DOFs for the LHS variable is "// &
              & TRIM(NumberToVString(totalNumberOfLHSDOFs,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF                    
          NULLIFY(residualMapping%jacobianMatrixToVarMaps(variableIdx)%ptr)
          CALL EquationsMappingJMToVMap_Initialise(residualMapping%jacobianMatrixToVarMaps(variableIdx)%ptr,err,error,*999)
          residualMapping%jacobianMatrixToVarMaps(variableIdx)%ptr%jacobianNumber=variableIdx          
          residualMapping%jacobianMatrixToVarMaps(variableIdx)%ptr%variableType=residualMapping%variableTypes(variableIdx)
          residualMapping%jacobianMatrixToVarMaps(variableIdx)%ptr%variable=>dependentVariable
          residualMapping%jacobianMatrixToVarMaps(variableIdx)%ptr%numberOfColumns=numberOfGlobal
          residualMapping%jacobianMatrixToVarMaps(variableIdx)%ptr%jacobianCoefficient=residualCoefficient
          ALLOCATE(residualMapping%jacobianMatrixToVarMaps(variableIdx)%ptr%equationsColumnToDOFVariableMap(numberOfGlobal), &
            & STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate equations column to DOF variable map.",err,error,*999)
          residualMapping%jacobianMatrixToVarMaps(variableIdx)%ptr%equationsColumnToDOFVariableMap=0
          DO dofIdx=1,totalNumberOfDependentDOFs
            !1-1 mapping for now
            CALL DomainMapping_LocalToGlobalGet(columnDomainMapping,dofIdx,columnIdx,err,error,*999)
            residualMapping%jacobianMatrixToVarMaps(variableIdx)%ptr%equationsColumnToDOFVariableMap(columnIdx)=dofIdx
          ENDDO !dofIdx
          residualMapping%jacobianMatrixToVarMaps(variableIdx)%ptr%columnDOFSMapping=>columnDomainMapping
        ENDDO !variableIdx
        residualMapping%residualCoefficient=residualCoefficient
      ENDDO !residualIdx
    ENDIF    
    !Calculate RHS mappings
    IF(createValuesCache%rhsVariableType/=0) THEN                  
      CALL EquationsMappingVector_RHSMappingInitialise(vectorMapping,err,error,*999)
      NULLIFY(rhsMapping)
      CALL EquationsMappingVector_RHSMappingGet(vectorMapping,rhsMapping,err,error,*999)
      rhsMapping%rhsVariableType=createValuesCache%rhsVariableType
      NULLIFY(dependentVariable)
      CALL Field_VariableGet(dependentField,createValuesCache%rhsVariableType,dependentVariable,err,error,*999)
      CALL FieldVariable_DomainMappingGet(dependentVariable,rhsMapping%rhsVariableMapping,err,error,*999)
      CALL FieldVariable_TotalNumberOfDOFsGet(dependentVariable,totalNumberOfDependentDOFs,err,error,*999)
      rhsMapping%rhsVariable=>dependentVariable
      rhsMapping%rhsCoefficient=createValuesCache%rhsCoefficient
      !Allocate and set up the row mappings
      ALLOCATE(rhsMapping%rhsDOFToEquationsRowMap(totalNumberOfDependentDOFs),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate RHS to equations row map.",err,error,*999)
      ALLOCATE(rhsMapping%equationsRowToRHSDOFMap(totalNumberOfRows),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate equations row to DOF map.",err,error,*999)
      DO dofIdx=1,totalNumberOfDependentDOFs
        !1-1 mapping for now
        rowIdx=dofIdx
        rhsMapping%rhsDOFToEquationsRowMap(dofIdx)=rowIdx
      ENDDO !dofIdx
      DO rowIdx=1,totalNumberOfRows
        !1-1 mapping for now
        dofIdx=rowIdx
        rhsMapping%equationsRowToRHSDOFMap(rowIdx)=dofIdx
      ENDDO !rowIdx
    ENDIF
    !Calcuate the source mappings
    IF(createValuesCache%numberOfSources>0) THEN
      NULLIFY(sourceField)
      CALL EquationsSet_SourceFieldGet(equationsSet,sourceField,err,error,*999)
      CALL EquationsMappingVector_SourcesMappingInitialise(vectorMapping,err,error,*999)
      NULLIFY(sourcesMapping)
      CALL EquationsMappingVector_SourcesMappingGet(vectorMapping,sourcesMapping,err,error,*999)
      ALLOCATE(sourcesMapping%sources(createValuesCache%numberOfSources),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate sources.",err,error,*999)
      sourcesMapping%numberOfSources=createValuesCache%numberOfSources
      DO sourceIdx=1,createValuesCache%numberOfSources
        NULLIFY(sourcesMapping%sources(sourceIdx)%ptr)
        CALL EquationsMappingSources_SourceMappingInitialise(sourcesMapping,sourceIdx,err,error,*999)
        NULLIFY(sourceMapping)
        CALL EquationsMappingSources_SourceMappingGet(sourcesMapping,sourceIdx,sourceMapping,err,error,*999)
        CALL EquationsMappingVectorCVC_SourceVariableTypeGet(createValuesCache,sourceIdx,sourceMapping%sourceVariableType, &
          & err,error,*999)
        NULLIFY(sourceVariable)
        CALL Field_VariableGet(sourceField,sourceMapping%sourceVariableType,sourceVariable,err,error,*999)
        sourceMapping%sourceVariable=>sourceVariable
        CALL EquationsMappingVectorCVC_SourceCoefficientGet(createValuesCache,sourceIdx,sourceMapping%sourceCoefficient, &
          & err,error,*999)
        !sourceMapping%sourceVariableMapping=>sourceVariable%domainMapping
        !!Allocate and set up the row mappings
        !ALLOCATE(sourceMapping%sourceDOFToEquationsRowMap(sourceVariable%totalNumberOfDofs),STAT=err)
        !IF(err/=0) CALL FlagError("Could not allocate source dof to equations row map.",err,error,*999)
        !ALLOCATE(sourceMapping%equationsRowToSourceDOFMap(lhsMapping%totalNumberOfRows),STAT=err)
        !IF(err/=0) CALL FlagError("Could not allocate equations row to source map.",err,error,*999)
        !DO dofIdx=1,sourceVariable%totalNumberOfDofs
        !  !1-1 mapping for now
        !  rowIdx=dofIdx
        !  sourceMapping%sourceDOFToEquationsRowMap(dofIdx)=rowIdx
        !ENDDO !dofIdx
        !DO rowIdx=1,lhsMapping%totalNumberOfRows
        !  !1-1 mapping for now
        !  dofIdx=rowIdx
        !  sourceMapping%equationsRowToSourceDOFMap(rowIdx)=dofIdx
        !ENDDO !rowIdx
      ENDDO !sourceIdx
    ENDIF
    
    IF(diagnostics1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Equations vector mappings:",err,error,*999)
      NULLIFY(lhsMapping)
      CALL EquationsMappingVector_LHSMappingGet(vectorMapping,lhsMapping,err,error,*999)
      CALL EquationsMappingLHS_NumberOfRowsGet(lhsMapping,numberOfRows,err,error,*999)
      CALL EquationsMappingLHS_TotalNumberOfRowsGet(lhsMapping,totalNumberOfRows,err,error,*999)
      CALL EquationsMappingLHS_NumberOfGlobalRowsGet(lhsMapping,numberOfGlobalRows,err,error,*999)
      NULLIFY(lhsVariable)
      CALL EquationsMappingLHS_LHSVariableGet(lhsMapping,lhsVariable,err,error,*999)
      CALL FieldVariable_VariableTypeGet(lhsVariable,lhsVariableType,err,error,*999)
      CALL FieldVariable_TotalNumberOfDOFsGet(lhsVariable,totalNumberOfLHSDOFs,err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  LHS mappings:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of rows = ",numberOfRows,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Total number of rows = ",totalNumberOfRows,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of global rows = ",numberOfGlobalRows,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    LHS variable type = ",lhsVariableType,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Total number of LHS DOFs = ",totalNumberOfLHSDOFs,err,error,*999)
      NULLIFY(lhsDOFToEquationsRowMap)
      CALL EquationsMappingLHS_LHSDOFToEquationsRowMapGet(lhsMapping,lhsDOFToEquationsRowMap,err,error,*999)
      NULLIFY(equationsRowToLHSDOFMap)
      CALL EquationsMappingLHS_EquationsRowToLHSDOFMapGet(lhsMapping,equationsRowToLHSDOFMap,err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Row mappings:",err,error,*999)
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,totalNumberOfLHSDOFs,5,5, &
        & lhsDOFToEquationsRowMap,'("    DOF to row mappings :",5(X,I13))','(25X,5(X,I13))',err,error,*999) 
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,lhsMapping%totalNumberOfRows,5,5, &
        & equationsRowToLHSDOFMap,'("    Row to DOF mappings :",5(X,I13))','(25X,5(X,I13))',err,error,*999) 
      NULLIFY(dynamicMapping)
      CALL EquationsMappingVector_DynamicMappingExists(vectorMapping,dynamicMapping,err,error,*999)
      IF(ASSOCIATED(dynamicMapping)) THEN
        NULLIFY(dynamicVariable)
        CALL EquationsMappingDynamic_DynamicVariableGet(dynamicMapping,dynamicVariable,err,error,*999)
        CALL FieldVariable_VariableTypeGet(dynamicVariable,totalNumberOfDependentDOFs,err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Dynamic mappings:",err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of dynamic equations matrices = ",dynamicMapping% &
          & numberOfDynamicMatrices,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Dynamic stiffness matrix number = ",dynamicMapping% &
          & stiffnessMatrixNumber,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Dynamic damping matrix number = ",dynamicMapping% &
          & dampingMatrixNumber,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Dynamic mass matrix number = ",dynamicMapping% &
          & massMatrixNumber,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Dynamic variable type = ",dynamicMapping%dynamicVariableType, &
          & err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Total number of DOFs = ",totalNumberOfDependentDOFs,err,error,*999)
        NULLIFY(varToEquationsMatricesMap)
        CALL EquationsMappingDynamic_VariableToEquationsMatricesMapGet(dynamicMapping,varToEquationsMatricesMap,err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Variable to matrices mappings:",err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of equations matrices = ",varToEquationsMatricesMap% &
          & numberOfEquationsMatrices,err,error,*999)
        IF(varToEquationsMatricesMap%numberOfEquationsMatrices>0) THEN
          NULLIFY(dependentVariable)
          CALL EquationsMappingVectorVToEMSMap_VariableGet(varToEquationsMatricesMap,dependentVariable,err,error,*999)
          CALL FieldVariable_TotalNumberOfDOFsGet(dependentVariable,totalNumberOfDependentDOFs,err,error,*999)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,varToEquationsMatricesMap%numberOfEquationsMatrices,4,4, &
            & varToEquationsMatricesMap%equationsMatrixNumbers,'("      Matrix numbers :",4(X,I12))','(22X,4(X,I12))', &
            & err,error,*999) 
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"      DOF to column maps :",err,error,*999)
          DO matrixIdx=1,varToEquationsMatricesMap%numberOfEquationsMatrices
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Matrix number : ",matrixIdx,err,error,*999)
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,totalNumberOfDependentDOFs,5,5, &
              & varToEquationsMatricesMap%dofToColumnsMaps(matrixIdx)%columnDOF, &
              & '("        Column numbers :",5(X,I13))','(24X,5(X,I13))',err,error,*999) 
          ENDDO !matrixIdx
        ENDIF
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Matrix to variable mappings:",err,error,*999)
        DO matrixIdx=1,dynamicMapping%numberOfDynamicMatrices
          NULLIFY(equationsMatrixToVarMap)
          CALL EquationsMappingDynamic_EquationsMatrixToVarMapGet(dynamicMapping,matrixIdx,equationsMatrixToVarMap,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Matrix number : ",matrixIdx,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Variable type = ",equationsMatrixToVarMap%variableType, &
            & err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Number of columns = ",equationsMatrixToVarMap%numberOfColumns, &
            & err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Matrix coefficient = ",equationsMatrixToVarMap%matrixCoefficient, &
            & err,error,*999)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,equationsMatrixToVarMap%numberOfColumns,5,5, &
            & equationsMatrixToVarMap%columnToDOFMap,'("        Column to DOF maps :",5(X,I13))','(28X,5(X,I13))',err,error,*999) 
        ENDDO !matrixIdx
      ENDIF
      NULLIFY(linearMapping)
      CALL EquationsMappingVector_LinearMappingExists(vectorMapping,linearMapping,err,error,*999)
      IF(ASSOCIATED(linearMapping)) THEN
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Linear mappings:",err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of linear equations matrices = ",linearMapping% &
          & numberOfLinearMatrices,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of linear variables = ",linearMapping% &
          & numberOfLinearVariables,err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,linearMapping%numberOfLinearVariables,4,4, &
          & linearMapping%linearVariableTypes,'("    Linear variable types :",4(X,I12))','(20X,4(X,I12))', &
          & err,error,*999) 
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Variable to matrices mappings:",err,error,*999)
        DO variableIdx=1,linearMapping%numberOfLinearVariables
          NULLIFY(varToEquationsMatricesMap)
          CALL EquationsMappingLinear_VariableToEquationsMatricesMapGet(linearMapping,variableIdx,varToEquationsMatricesMap, &
            & err,error,*999)
          NULLIFY(linearVariable)
          CALL EquationsMappingVectorVToEMSMap_VariableGet(varToEquationsMatricesMap,linearVariable,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Variable index : ",variableIdx,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Total number of DOFs = ",linearVariable%totalNumberOfDOFs, &
            & err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of equations matrices = ",varToEquationsMatricesMap% &
            & numberOfEquationsMatrices,err,error,*999)
          IF(varToEquationsMatricesMap%numberOfEquationsMatrices>0) THEN
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,varToEquationsMatricesMap%numberOfEquationsMatrices,4,4, &
              & varToEquationsMatricesMap%equationsMatrixNumbers,'("      Matrix numbers :",4(X,I12))','(22X,4(X,I12))', &
              & err,error,*999)
            CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"      DOF to column maps :",err,error,*999)
            DO matrixIdx=1,varToEquationsMatricesMap%numberOfEquationsMatrices
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Matrix number : ",matrixIdx,err,error,*999)
              CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,linearVariable%totalNumberOfDOFs,5,5, &
                & varToEquationsMatricesMap%dofToColumnsMaps(matrixIdx)%columnDOF, &
                & '("        Column numbers :",5(X,I13))','(24X,5(X,I13))',err,error,*999) 
            ENDDO !matrixIdx
          ENDIF
        ENDDO !variableIdx
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Matrix to variable mappings:",err,error,*999)
        DO matrixIdx=1,linearMapping%numberOfLinearMatrices
          NULLIFY(equationsMatrixToVarMap)
          CALL EquationsMappingLinear_EquationsMatrixToVarMapGet(linearMapping,matrixIdx,equationsMatrixToVarMap,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Matrix number : ",matrixIdx,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Variable index = ",equationsMatrixToVarMap%variableIndex, &
            & err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Variable type = ",equationsMatrixToVarMap%variableType, &
            & err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Number of columns = ",equationsMatrixToVarMap%numberOfColumns, &
            & err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Matrix coefficient = ",equationsMatrixToVarMap%matrixCoefficient, &
            & err,error,*999)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,equationsMatrixToVarMap%numberOfColumns,5,5, &
            & equationsMatrixToVarMap%columnToDOFMap,'("        Column to DOF maps :",5(X,I13))','(28X,5(X,I13))', &
            & err,error,*999)          
        ENDDO !matrixIdx
      ENDIF
      NULLIFY(nonlinearMapping)
      CALL EquationsMappingVector_NonlinearMappingExists(vectorMapping,nonlinearMapping,err,error,*999)
      IF(ASSOCIATED(nonlinearMapping)) THEN
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Nonlinear mappings:",err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of residuals = ",nonlinearMapping%numberOfResiduals,err,error,*999)
        DO residualIdx=1,nonlinearMapping%numberOfResiduals
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Residual number : ",residualIdx,err,error,*999)
          NULLIFY(residualMapping)
          CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,residualIdx,residualMapping,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Residual coefficient = ",residualMapping%residualCoefficient, &
            & err,error,*999)
          NULLIFY(equationsRowToResidualDOFMap)
          CALL EquationsMappingResidual_EquationsRowToResidualDOFMapGet(residualMapping,equationsRowToResidualDOFMap, &
            & err,error,*999)
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"      Residual row mappings:",err,error,*999)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,lhsMapping%totalNumberOfRows,5,5, &
            & equationsRowToResidualDOFMap,'("      Row to DOF mappings :",5(X,I13))','(27X,5(X,I13))',err,error,*999) 
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of variables = ",residualMapping%numberOfVariables, &
            & err,error,*999)
          DO variableIdx=1,residualMapping%numberOfVariables
            NULLIFY(jacobianMatrixToVarMap)
            CALL EquationsMappingResidual_JacobianMatrixToVarMapGet(residualMapping,variableIdx,jacobianMatrixToVarMap, &
              & err,error,*999)
            NULLIFY(varToJacobianMatrixMap)
            CALL EquationsMappingResidual_VarToJacobianMatrixMapGet(residualMapping,variableIdx,varToJacobianMatrixMap, &
              & err,error,*999)
            NULLIFY(jacobianVariable)
            CALL EquationsMappingVectorJMToVMap_VariableGet(jacobianMatrixToVarMap,jacobianVariable,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Variable number : ",variableIdx,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Residual variable type = ",jacobianMatrixToVarMap%variableType, &
              & err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Total number of variable DOFs = ",jacobianVariable% &
              & totalNumberOfDofs,err,error,*999)
            CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"        Jacobian mappings:",err,error,*999)
            CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"          Variable to Jacobian mappings:",err,error,*999)
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,jacobianVariable%totalNumberOfDofs,5,5, &
              & varToJacobianMatrixMap%dofToColumnsMap,'("            DOF to column map :",5(X,I13))','(32X,5(X,I13))', &
              & err,error,*999) 
            CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"          Jacobian to variable mappings:",err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Number of columns = ",jacobianMatrixToVarMap% &
              & numberOfColumns,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Jacobian coefficient = ",jacobianMatrixToVarMap% &
              & jacobianCoefficient,err,error,*999)
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,jacobianMatrixToVarMap%numberOfColumns,5,5, &
              & jacobianMatrixToVarMap%equationsColumnToDOFVariableMap, &
              & '("            Column to DOF map :",5(X,I13))','(32X,5(X,I13))',err,error,*999) 
          ENDDO !variableIdx
        ENDDO !residualIdx
      ENDIF
      NULLIFY(rhsMapping)
      CALL EquationsMappingVector_RHSMappingExists(vectorMapping,rhsMapping,err,error,*999)
      IF(ASSOCIATED(rhsMapping)) THEN
        NULLIFY(rhsVariable)
        CALL EquationsMappingRHS_RHSVariableGet(rhsMapping,rhsVariable,err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  RHS mappings:",err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    RHS variable type = ",rhsMapping%rhsVariableType,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Total number of RHS DOFs = ",rhsVariable%totalNumberOfDofs, &
          & err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    RHS coefficient = ",rhsMapping%rhsCoefficient,err,error,*999)
      ENDIF
      NULLIFY(sourcesMapping)
      CALL EquationsMappingVector_SourcesMappingExists(vectorMapping,sourcesMapping,err,error,*999)
      IF(ASSOCIATED(sourcesMapping)) THEN
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Source mappings:",err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of sources = ",sourcesMapping%numberOfSources,err,error,*999)
        DO sourceIdx=1,sourcesMapping%numberOfSources
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Source number : ",sourceIdx,err,error,*999)
          NULLIFY(sourceMapping)
          CALL EquationsMappingSources_SourceMappingGet(sourcesMapping,sourceIdx,sourceMapping,err,error,*999)
          NULLIFY(sourceVariable)
          CALL EquationsMappingSource_SourceVariableGet(sourceMapping,sourceVariable,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Source variable type = ",sourceMapping%sourceVariableType, &
            & err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Total number of source DOFs = ",sourceVariable% &
            & totalNumberOfDofs,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Source coefficient = ",sourceMapping%sourceCoefficient, &
            & err,error,*999)
        ENDDO !sourceIdx
      ENDIF
    ENDIF
       
    EXITS("EquationsMappingVector_Calculate")
    RETURN
999 IF(ALLOCATED(variableTypes)) DEALLOCATE(variableTypes)
    ERRORSEXITS("EquationsMappingVector_Calculate",err,error)    
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_Calculate

  !
  !================================================================================================================================
  !

  !>Finalises an equations mapping create values cache and deallocates all memory
  SUBROUTINE EquationsMappingVector_CreateValuesCacheFinalise(createValuesCache,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorCreateValuesCacheType), POINTER :: createValuesCache !<A pointer to the create values cache to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsMappingVector_CreateValuesCacheFinalise",err,error,*999)

    IF(ASSOCIATED(createValuesCache)) THEN
      IF(ALLOCATED(createValuesCache%dynamicMatrixCoefficients)) DEALLOCATE(createValuesCache%dynamicMatrixCoefficients)
      IF(ALLOCATED(createValuesCache%linearMatrixVariableTypes)) DEALLOCATE(createValuesCache%linearMatrixVariableTypes)
      IF(ALLOCATED(createValuesCache%linearMatrixCoefficients)) DEALLOCATE(createValuesCache%linearMatrixCoefficients)
      IF(ALLOCATED(createValuesCache%numberOfResidualVariables)) DEALLOCATE(createValuesCache%numberOfResidualVariables)
      IF(ALLOCATED(createValuesCache%residualVariableTypes)) DEALLOCATE(createValuesCache%residualVariableTypes)
      IF(ALLOCATED(createValuesCache%residualCoefficients)) DEALLOCATE(createValuesCache%residualCoefficients)
      IF(ALLOCATED(createValuesCache%sourceVariableTypes)) DEALLOCATE(createValuesCache%sourceVariableTypes)
      IF(ALLOCATED(createValuesCache%sourceCoefficients)) DEALLOCATE(createValuesCache%sourceCoefficients)
      DEALLOCATE(createValuesCache)
    ENDIF
       
    EXITS("EquationsMappingVector_CreateValuesCacheFinalise")
    RETURN
999 ERRORS("EquationsMappingVector_CreateValuesCacheFinalise",err,error)
    EXITS("EquationsMappingVector_CreateValuesCacheFinalise")
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_CreateValuesCacheFinalise

  !
  !================================================================================================================================
  !

  !>Initialises a vector equations mapping create values cache 
  SUBROUTINE EquationsMappingVector_CreateValuesCacheInitialise(vectorMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the vector equations mapping to initialise the create values cache for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,dynamicVariableType,fluxVariableType,linearVariableType,matrixIdx,matrixIdx2,residualIdx, &
      & residualVariableType,variableIdx,variableNumber
    LOGICAL :: isResidualType
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(FieldType), POINTER :: dependentField
    TYPE(FieldVariableType), POINTER :: dynamicVariable,fluxVariable,linearVariable,residualVariable
    TYPE(VARYING_STRING) :: dummyError,localError

    ENTERS("EquationsMappingVector_CreateValuesCacheInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector equations mapping is not associated.",err,error,*998)
    IF(ASSOCIATED(vectorMapping%createValuesCache)) &
      & CALL FlagError("Equations mapping create values cache is already associated.",err,error,*998)
    
    NULLIFY(vectorEquations)
    CALL EquationsMappingVector_VectorEquationsGet(vectorMapping,vectorEquations,err,error,*999)
    NULLIFY(equations)
    CALL EquationsVector_EquationsGet(vectorEquations,equations,err,error,*999)
    NULLIFY(equationsSet)
    CALL Equations_EquationsSetGet(equations,equationsSet,err,error,*999)
    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
    
    !Allocate and initialise the create values cache
    ALLOCATE(vectorMapping%createValuesCache,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate equations mapping create values cache.",err,error,*999)
    vectorMapping%createValuesCache%numberOfDynamicMatrices=0
    vectorMapping%createValuesCache%dynamicVariableType=0
    vectorMapping%createValuesCache%dynamicStiffnessMatrixNumber=0
    vectorMapping%createValuesCache%dynamicDampingMatrixNumber=0
    vectorMapping%createValuesCache%dynamicMassMatrixNumber=0
    vectorMapping%createValuesCache%numberOfLinearMatrices=0    
    vectorMapping%createValuesCache%numberOfResiduals=0
    vectorMapping%createValuesCache%rhsVariableType=0
    vectorMapping%createValuesCache%rhsCoefficient=1.0_DP
    vectorMapping%createValuesCache%numberOfSources=0
!!TODO: SHOULD WE HAVE ANY DEFAULTS HERE? JUST LET THE USER SET EXACTLY WHAT THEY USE. THEY DON'T HAVE TO USE ALL VARIABLES
!!      IN THE DEPENDENT FIELD    
    !Set the default equations mapping in the create values cache.
    !First calculate how many linear and dynamic matrices we have and set the variable types for the dynamic, residual
    !and RHS variables
    IF(dependentField%numberOfVariables<1) THEN
      localError="The number of variables of "//TRIM(NumberToVString(dependentField%numberOfVariables,"*",err,error))// &
        & " for dependent field number "//TRIM(NumberToVString(dependentField%userNumber,"*",err,error))// &
        & " is invalid. The number of variables should be >= 1."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    SELECT CASE(equations%timeDependence)
      !Static equations
    CASE(EQUATIONS_STATIC,EQUATIONS_QUASISTATIC)
      SELECT CASE(equations%linearity)
      CASE(EQUATIONS_LINEAR,EQUATIONS_NONLINEAR_BCS)
        !Static, linear equations
        IF(dependentField%numberOfVariables==1) THEN
          !Only one variable so map it to a linear matrix and have no RHS.
          vectorMapping%createValuesCache%numberOfLinearMatrices=1
        ELSE
          !Map the second variable to the RHS and all the other variables to linear matrices
          vectorMapping%createValuesCache%numberOfLinearMatrices=dependentField%numberOfVariables-1
          NULLIFY(fluxVariable)
          CALL Field_VariableIndexGet(dependentField,2,fluxVariable,fluxVariableType,err,error,*999)
          vectorMapping%createValuesCache%rhsVariableType=fluxVariableType
        ENDIF
      CASE(EQUATIONS_NONLINEAR)
        !Static, nonlinear equations
        IF(dependentField%numberOfVariables==1) THEN
          !Only one variable so map it to a residual and have no RHS.
          vectorMapping%createValuesCache%numberOfResidualVariables=1
        ELSE         
          !Map first variable to a residual and the second variable to the RHS.
          vectorMapping%createValuesCache%numberOfResidualVariables=1
          NULLIFY(fluxVariable)
          CALL Field_VariableIndexGet(dependentField,2,fluxVariable,fluxVariableType,err,error,*999)
          vectorMapping%createValuesCache%rhsVariableType=fluxVariableType
        ENDIF
      CASE DEFAULT
        localError="The equations linearity type of "//TRIM(NumberToVString(equations%linearity,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_FIRST_ORDER_DYNAMIC,EQUATIONS_SECOND_ORDER_DYNAMIC)
      !Dynamic equations
      IF(equations%timeDependence==EQUATIONS_FIRST_ORDER_DYNAMIC) THEN
        vectorMapping%createValuesCache%numberOfDynamicMatrices=2
        vectorMapping%createValuesCache%dynamicStiffnessMatrixNumber=1
        vectorMapping%createValuesCache%dynamicDampingMatrixNumber=2
      ELSE
        vectorMapping%createValuesCache%numberOfDynamicMatrices=3
        vectorMapping%createValuesCache%dynamicStiffnessMatrixNumber=1
        vectorMapping%createValuesCache%dynamicDampingMatrixNumber=2
        vectorMapping%createValuesCache%dynamicMassMatrixNumber=3
      ENDIF
      !Map the first variable to dynamic matrices
      NULLIFY(dynamicVariable)
      CALL Field_VariableIndexGet(dependentField,1,dynamicVariable,dynamicVariableType,err,error,*999)
      vectorMapping%createValuesCache%dynamicVariableType=dynamicVariableType
      SELECT CASE(equations%linearity)
      CASE(EQUATIONS_LINEAR,EQUATIONS_NONLINEAR_BCS)
        !Dynamic linear equations
        IF(dependentField%numberOfVariables>1) THEN
          !Map the second variable to the RHS
          NULLIFY(fluxVariable)
          CALL Field_VariableIndexGet(dependentField,2,fluxVariable,fluxVariableType,err,error,*999)
          vectorMapping%createValuesCache%rhsVariableType=fluxVariableType
        ENDIF
      CASE(EQUATIONS_NONLINEAR)
        !Dynamic, nonlinear equations.
        !Have one residual (same as the dynamic variable)
        vectorMapping%createValuesCache%numberOfResidualVariables=1
        IF(dependentField%numberOfVariables>1) THEN
          !Map the second variable to the RHS
          NULLIFY(fluxVariable)
          CALL Field_VariableIndexGet(dependentField,2,fluxVariable,fluxVariableType,err,error,*999)
          vectorMapping%createValuesCache%rhsVariableType=fluxVariableType
        ENDIF
      CASE DEFAULT
        localError="The equations linearity type of "//TRIM(NumberToVString(equations%linearity,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The equations time dependence type of "//TRIM(NumberToVString(equations%timeDependence,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    !Allocate the dynamic matrix coefficients and set their values
    IF(vectorMapping%createValuesCache%numberOfDynamicMatrices>0) THEN
      ALLOCATE(vectorMapping%createValuesCache%dynamicMatrixCoefficients(vectorMapping% &
        & createValuesCache%numberOfDynamicMatrices),STAT=err)
      IF(err/=0) &
        & CALL FlagError("Could not allocate equations mapping create values cache dynamic matrix coefficients.",err,error,*999)
      vectorMapping%createValuesCache%dynamicMatrixCoefficients=1.0_DP !Equations matrices are added by default
    ENDIF
    !Allocate the residual variable types
    IF(vectorMapping%createValuesCache%numberOfResiduals>0) THEN
      ALLOCATE(vectorMapping%createValuesCache%numberOfResidualVariables(vectorMapping% &
        & createValuesCache%numberOfResiduals),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate create values cache number of residual variables.",err,error,*999)
      !Just have one variable per residual
      vectorMapping%createValuesCache%numberOfResidualVariables=1
      ALLOCATE(vectorMapping%createValuesCache%residualVariableTypes(1,vectorMapping%createValuesCache%numberOfResiduals), &
        & STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate the create values cache residual variable types.",err,error,*999)
      !Whether the equations are static or dynamic the residual variable will be the first dependent variable
      NULLIFY(residualVariable)
      CALL Field_VariableIndexGet(dependentField,1,residualVariable,residualVariableType,err,error,*999)
      vectorMapping%createValuesCache%residualVariableTypes=residualVariableType
      !Set up the coefficients
      ALLOCATE(vectorMapping%createValuesCache%residualCoefficients(vectorMapping%createValuesCache%numberOfResiduals), &
        & STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate create values cache residual coefficients.",err,error,*999)
      vectorMapping%createValuesCache%residualCoefficients=1.0_DP
    ENDIF
    !Allocate the linear matrix variable types and linear matrix coefficients and set their values
    IF(vectorMapping%createValuesCache%numberOfLinearMatrices>0) THEN
      ALLOCATE(vectorMapping%createValuesCache%linearMatrixVariableTypes(vectorMapping%createValuesCache%numberOfLinearMatrices), &
        & STAT=err)
      IF(err/=0) & 
        & CALL FlagError("Could not allocate equations mapping create values cache linear matrix variable types.",err,error,*999)
      ALLOCATE(vectorMapping%createValuesCache%linearMatrixCoefficients(vectorMapping% &
        & createValuesCache%numberOfLinearMatrices),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate equations mapping create values cache linear matrix coefficients.", &
        & err,error,*999)
      !Set up the matrices variable types
      vectorMapping%createValuesCache%linearMatrixVariableTypes=0
      !The first linear matrix will be the first variable
      NULLIFY(linearVariable)
      CALL Field_VariableIndexGet(dependentField,1,linearVariable,linearVariableType,err,error,*999)
      vectorMapping%createValuesCache%linearMatrixVariableTypes(1)=linearVariableType
      !Subsequent linear matrices will be the third (and onwards) variable
      DO matrixIdx=2,vectorMapping%createValuesCache%numberOfLinearMatrices
        NULLIFY(linearVariable)
        CALL Field_VariableIndexGet(dependentField,matrixIdx+1,linearVariable,linearVariableType,err,error,*999)
        vectorMapping%createValuesCache%linearMatrixVariableTypes(matrixIdx)=linearVariableType
      ENDDO !matrixIdx
      vectorMapping%createValuesCache%linearMatrixCoefficients=1.0_DP !Equations matrices are added by default
    ENDIF
      
    EXITS("EquationsMappingVector_CreateValuesCacheInitialise")
    RETURN
999 CALL EquationsMappingVector_CreateValuesCacheFinalise(vectorMapping%createValuesCache,dummyErr,dummyError,*998)
998 ERRORS("EquationsMappingVector_CreateValuesCacheInitialise",err,error)
    EXITS("EquationsMappingVector_CreateValuesCacheInitialise")
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_CreateValuesCacheInitialise

  !
  !================================================================================================================================
  !

  !>Finalises the equations mapping dynamic mapping and deallocates all memory
  SUBROUTINE EquationsMappingVector_DynamicMappingFinalise(dynamicMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping !<A pointer to the dynamic mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx,variableType
 
    ENTERS("EquationsMappingVector_DynamicMappingFinalise",err,error,*999)

    IF(ASSOCIATED(dynamicMapping)) THEN
      CALL EquationsMappingVToEMSMap_Finalise(dynamicMapping%varToEquationsMatricesMap,err,error,*999)
      IF(ALLOCATED(dynamicMapping%equationsMatrixToVarMaps)) THEN
        DO matrixIdx=1,SIZE(dynamicMapping%equationsMatrixToVarMaps,1)
          CALL EquationsMappingEMToVMap_Finalise(dynamicMapping%equationsMatrixToVarMaps(matrixIdx)%ptr, &
            & err,error,*999)
        ENDDO !matrixIdx
        DEALLOCATE(dynamicMapping%equationsMatrixToVarMaps)
      ENDIF
      IF(ASSOCIATED(dynamicMapping%equationsRowToVariableDOFMaps)) DEALLOCATE(dynamicMapping%equationsRowToVariableDOFMaps)
      DEALLOCATE(dynamicMapping)
    ENDIF
       
    EXITS("EquationsMappingVector_DynamicMappingFinalise")
    RETURN
999 ERRORSEXITS("EquationsMappingVector_DynamicMappingFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_DynamicMappingFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the equations mapping dynamic mapping
  SUBROUTINE EquationsMappingVector_DynamicMappingInitialise(vectorMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations mapping to initialise the dynamic mapping for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("EquationsMappingVector_DynamicMappingInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector equations mapping is not associated.",err,error,*998)
    IF(ASSOCIATED(vectorMapping%dynamicMapping)) &
      & CALL FlagError("Equations mapping dynamic mapping is already associated.",err,error,*998)
    
    ALLOCATE(vectorMapping%dynamicMapping,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate equations mapping dynamic mapping.",err,error,*999)
    vectorMapping%dynamicMapping%vectorMapping=>vectorMapping
    vectorMapping%dynamicMapping%numberOfDynamicMatrices=0
    vectorMapping%dynamicMapping%stiffnessMatrixNumber=0
    vectorMapping%dynamicMapping%dampingMatrixNumber=0
    vectorMapping%dynamicMapping%massMatrixNumber=0
    vectorMapping%dynamicMapping%dynamicVariableType=0
    NULLIFY(vectorMapping%dynamicMapping%dynamicVariable)
    NULLIFY(vectorMapping%dynamicMapping%varToEquationsMatricesMap)
    NULLIFY(vectorMapping%dynamicMapping%equationsRowToVariableDOFMaps)
    
    EXITS("EquationsMappingVector_DynamicMappingInitialise")
    RETURN
999 CALL EquationsMappingVector_DynamicMappingFinalise(vectorMapping%dynamicMapping,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMappingVector_DynamicMappingInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_DynamicMappingInitialise

  !
  !================================================================================================================================
  !

  !>Sets/changes the matrices involved in dynamic equations mapping
  SUBROUTINE EquationsMappingVector_DynamicMatricesSetAll(vectorMapping,massMatrix,dampingMatrix,stiffnessMatrix,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the vector equations mapping to set the atrices for
    LOGICAL, INTENT(IN) :: massMatrix !<Is .TRUE. if the mass matrix is in the vector equations mapping, .FALSE. if not
    LOGICAL, INTENT(IN) :: dampingMatrix !<Is .TRUE. if the damping matrix is in the vector equations mapping, .FALSE. if not
    LOGICAL, INTENT(IN) :: stiffnessMatrix !<Is .TRUE. if the stiffness matrix is in the vector equations mapping, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: newDynamicDampingMatrixNumber,newDynamicMassMatrixNumber,newDynamicStiffnessMatrixNumber, &
      & numberOfDynamicMatrices
    REAL(DP), ALLOCATABLE :: newDynamicMatrixCoefficients(:)
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsMappingVector_DynamicMatricesSetAll",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector equations mapping is not associated.",err,error,*999)
    
    NULLIFY(createValuesCache)
    CALL EquationsMappingVector_CreateValuesCacheGet(vectorMapping,createValuesCache,err,error,*999)
    NULLIFY(vectorEquations)
    CALL EquationsMappingVector_VectorEquationsGet(vectorMapping,vectorEquations,err,error,*999)
    NULLIFY(equations)
    CALL EquationsVector_EquationsGet(vectorEquations,equations,err,error,*999)
    CALL Equations_AssertIsDynamic(equations,err,error,*999)
     
    numberOfDynamicMatrices=0
    newDynamicStiffnessMatrixNumber=0
    newDynamicDampingMatrixNumber=0
    newDynamicMassMatrixNumber=0
    IF(stiffnessMatrix) THEN
      numberOfDynamicMatrices=numberOfDynamicMatrices+1
      newDynamicStiffnessMatrixNumber=numberOfDynamicMatrices
    ENDIF
    IF(dampingMatrix) THEN
      numberOfDynamicMatrices=numberOfDynamicMatrices+1
      newDynamicDampingMatrixNumber=numberOfDynamicMatrices
    ENDIF
    IF(massMatrix) THEN
      numberOfDynamicMatrices=numberOfDynamicMatrices+1
      newDynamicMassMatrixNumber=numberOfDynamicMatrices
    ENDIF
    IF(numberOfDynamicMatrices<=0) &
      & CALL FlagError("Invalid dynamic matrices set up. There are no dynamic equations matrices.",err,error,*999)
    ALLOCATE(newDynamicMatrixCoefficients(numberOfDynamicMatrices),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new dynamic matrix coefficients.",err,error,*999)
    IF(newDynamicStiffnessMatrixNumber/=0) THEN
      IF(createValuesCache%dynamicStiffnessMatrixNumber==0) THEN
        newDynamicMatrixCoefficients(newDynamicStiffnessMatrixNumber)=1.0_DP
      ELSE
        newDynamicMatrixCoefficients(newDynamicStiffnessMatrixNumber)= &
          & createValuesCache%dynamicMatrixCoefficients(createValuesCache%dynamicStiffnessMatrixNumber)
      ENDIF
    ENDIF
    IF(newDynamicDampingMatrixNumber/=0) THEN
      IF(createValuesCache%dynamicDampingMatrixNumber==0) THEN
        newDynamicMatrixCoefficients(newDynamicDampingMatrixNumber)=1.0_DP
      ELSE
        newDynamicMatrixCoefficients(newDynamicDampingMatrixNumber)= &
          & createValuesCache%dynamicMatrixCoefficients(createValuesCache%dynamicDampingMatrixNumber)
      ENDIF
    ENDIF
    IF(newDynamicMassMatrixNumber/=0) THEN
      IF(createValuesCache%dynamicMassMatrixNumber==0) THEN
        newDynamicMatrixCoefficients(newDynamicMassMatrixNumber)=1.0_DP
      ELSE
        newDynamicMatrixCoefficients(newDynamicMassMatrixNumber)= &
          & createValuesCache%dynamicMatrixCoefficients(createValuesCache%dynamicMassMatrixNumber)
      ENDIF
    ENDIF
    createValuesCache%numberOfDynamicMatrices=numberOfDynamicMatrices
    createValuesCache%dynamicStiffnessMatrixNumber=newDynamicStiffnessMatrixNumber
    createValuesCache%dynamicDampingMatrixNumber=newDynamicDampingMatrixNumber
    createValuesCache%dynamicMassMatrixNumber=newDynamicMassMatrixNumber
    CALL MOVE_ALLOC(newDynamicMatrixCoefficients,createValuesCache%dynamicMatrixCoefficients)
    
    EXITS("EquationsMappingVector_DynamicMatricesSetAll")
    RETURN
999 IF(ALLOCATED(newDynamicMatrixCoefficients)) DEALLOCATE(newDynamicMatrixCoefficients)    
    ERRORSEXITS("EquationsMappingVector_DynamicMatricesSetAll",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_DynamicMatricesSetAll

  !
  !================================================================================================================================
  !

  !>Sets/changes the matrices involved in a first order dynamic vector equations mapping
  SUBROUTINE EquationsMappingVector_DynamicMatricesSet1(vectorMapping,dampingMatrix,stiffnessMatrix,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the vector equations mapping to set the first order matrices for
    LOGICAL, INTENT(IN) :: dampingMatrix !<Is .TRUE. if the damping matrix is in the vector equations mapping, .FALSE. if not
    LOGICAL, INTENT(IN) :: stiffnessMatrix !<Is .TRUE. if the stiffness matrix is in the vector equations mapping, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsMappingVector_DynamicMatricesSet1",err,error,*999)

    CALL EquationsMappingVector_AssertNotFinished(vectorMapping,err,error,*999)    
    NULLIFY(vectorEquations)
    CALL EquationsMappingVector_VectorEquationsGet(vectorMapping,vectorEquations,err,error,*999)
    NULLIFY(equations)
    CALL EquationsVector_EquationsGet(vectorEquations,equations,err,error,*999)
    CALL Equations_AssertIsDynamic(equations,err,error,*999)
    
    IF(equations%timeDependence==EQUATIONS_FIRST_ORDER_DYNAMIC) THEN
      IF(.NOT.dampingMatrix) CALL FlagWarning("No damping matrix for first order dynamic equations.",err,error,*999)
      CALL EquationsMappingVector_DynamicMatricesSetAll(vectorMapping,.FALSE.,dampingMatrix,stiffnessMatrix,err,error,*999)
    ELSE
      CALL FlagError("Need to specify three matrices to set for second order dynamic equations.",err,error,*999)
    ENDIF
    
    EXITS("EquationsMappingVector_DynamicMatricesSet1")
    RETURN
999 ERRORSEXITS("EquationsMappingVector_DynamicMatricesSet1",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_DynamicMatricesSet1

  !
  !================================================================================================================================
  !

  !>Sets/changes the matrices involved in a second order dynamic equations mapping
  SUBROUTINE EquationsMappingVector_DynamicMatricesSet2(vectorMapping,massMatrix,dampingMatrix,stiffnessMatrix,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the vector equations mapping to set the first order matrices for
    LOGICAL, INTENT(IN) :: massMatrix !<Is .TRUE. if the mass matrix is in the vector equations mapping, .FALSE. if not
    LOGICAL, INTENT(IN) :: dampingMatrix !<Is .TRUE. if the damping matrix is in the vector equations mapping, .FALSE. if not
    LOGICAL, INTENT(IN) :: stiffnessMatrix !<Is .TRUE. if the stiffness matrix is in the vector equations mapping, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsMappingVector_DynamicMatricesSet2",err,error,*999)

    CALL EquationsMappingVector_AssertNotFinished(vectorMapping,err,error,*999)    
    NULLIFY(vectorEquations)
    CALL EquationsMappingVector_VectorEquationsGet(vectorMapping,vectorEquations,err,error,*999)
    NULLIFY(equations)
    CALL EquationsVector_EquationsGet(vectorEquations,equations,err,error,*999)
    CALL Equations_AssertIsDynamic(equations,err,error,*999)
 
    IF(equations%timeDependence==EQUATIONS_FIRST_ORDER_DYNAMIC) THEN
      IF(massMatrix) THEN
        CALL FlagError("The mass matrix cannot be present for first order dynamic equations.",err,error,*999)
      ELSE
        IF(.NOT.dampingMatrix) CALL FlagWarning("No damping matrix for a first order dynamic system.",err,error,*999)
        CALL EquationsMappingVector_DynamicMatricesSetAll(vectorMapping,.FALSE.,dampingMatrix,stiffnessMatrix,err,error,*999)
      ENDIF
    ELSE
      IF(.NOT.massMatrix) CALL FlagWarning("No mass matrix for a second order dynamic system.",err,error,*999)
      CALL EquationsMappingVector_DynamicMatricesSetAll(vectorMapping,massMatrix,dampingMatrix,stiffnessMatrix,err,error,*999)
    ENDIF
    
    EXITS("EquationsMappingVector_DynamicMatricesSet2")
    RETURN
999 ERRORSEXITS("EquationsMappingVector_DynamicMatricesSet2",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_DynamicMatricesSet2

  !
  !================================================================================================================================
  !

  !>Sets/changes the matrix coefficients in a first order dynamic equations mapping
  SUBROUTINE EquationsMappingVector_DynamicMatricesCoefficientsSet1(vectorMapping,dampingMatrixCoefficient, &
    & stiffnessMatrixCoefficient,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations mapping to set 
    REAL(DP), INTENT(IN) :: dampingMatrixCoefficient !<The damping matrix coefficient
    REAL(DP), INTENT(IN) :: stiffnessMatrixCoefficient !<The stiffness matrix coefficient
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsMappingVector_DynamicMatricesCoefficientsSet1",err,error,*999)
    
    CALL EquationsMappingVector_AssertNotFinished(vectorMapping,err,error,*999)     
    NULLIFY(createValuesCache)
    CALL EquationsMappingVector_CreateValuesCacheGet(vectorMapping,createValuesCache,err,error,*999)
    NULLIFY(vectorEquations)
    CALL EquationsMappingVector_VectorEquationsGet(vectorMapping,vectorEquations,err,error,*999)
    NULLIFY(equations)
    CALL EquationsVector_EquationsGet(vectorEquations,equations,err,error,*999)
    CALL Equations_AssertIsDynamic(equations,err,error,*999)
     
    IF(equations%timeDependence==EQUATIONS_FIRST_ORDER_DYNAMIC) THEN
      IF(createValuesCache%dynamicStiffnessMatrixNumber/=0) THEN
        createValuesCache%dynamicMatrixCoefficients(createValuesCache%dynamicStiffnessMatrixNumber)=stiffnessMatrixCoefficient
      ENDIF
      IF(createValuesCache%dynamicDampingMatrixNumber/=0) THEN
        createValuesCache%dynamicMatrixCoefficients(createValuesCache%dynamicDampingMatrixNumber)=dampingMatrixCoefficient
      ENDIF
    ELSE
      CALL FlagError("Need to specify three matrix coefficients for second order dynamic equations.", &
        & err,error,*999)
    ENDIF
    
    EXITS("EquationsMappingVector_DynamicMatricesCoefficientsSet1")
    RETURN
999 ERRORSEXITS("EquationsMappingVector_DynamicMatricesCoefficientsSet1",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_DynamicMatricesCoefficientsSet1

  !
  !================================================================================================================================
  !

  !>Sets/changes the matrix coefficients in a second order dynamic equations mapping
  SUBROUTINE EquationsMappingVector_DynamicMatricesCoefficientsSet2(vectorMapping,massMatrixCoefficient, &
    & dampingMatrixCoefficient,stiffnessMatrixCoefficient,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations mapping to set 
    REAL(DP), INTENT(IN) :: massMatrixCoefficient !<The mass matrix coefficient
    REAL(DP), INTENT(IN) :: dampingMatrixCoefficient !<The damping matrix coefficient
    REAL(DP), INTENT(IN) :: stiffnessMatrixCoefficient !<The stiffness matrix coefficient
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsMappingVector_DynamicMatricesCoefficientsSet2",err,error,*999)

    CALL EquationsMappingVector_AssertNotFinished(vectorMapping,err,error,*999)     
    NULLIFY(createValuesCache)
    CALL EquationsMappingVector_CreateValuesCacheGet(vectorMapping,createValuesCache,err,error,*999)
    NULLIFY(vectorEquations)
    CALL EquationsMappingVector_VectorEquationsGet(vectorMapping,vectorEquations,err,error,*999)
    NULLIFY(equations)
    CALL EquationsVector_EquationsGet(vectorEquations,equations,err,error,*999)
    CALL Equations_AssertIsDynamic(equations,err,error,*999)
    
    IF(equations%timeDependence==EQUATIONS_FIRST_ORDER_DYNAMIC) THEN
      CALL FlagError("Need to specify two matrix coefficients for second order dynamic equations.",err,error,*999)
    ELSE
      IF(createValuesCache%dynamicStiffnessMatrixNumber/=0) &
        & createValuesCache%dynamicMatrixCoefficients(createValuesCache%dynamicStiffnessMatrixNumber)=stiffnessMatrixCoefficient
      IF(createValuesCache%dynamicDampingMatrixNumber/=0) &
        & createValuesCache%dynamicMatrixCoefficients(createValuesCache%dynamicDampingMatrixNumber)=dampingMatrixCoefficient
      IF(createValuesCache%dynamicMassMatrixNumber/=0) &
        & createValuesCache%dynamicMatrixCoefficients(createValuesCache%dynamicMassMatrixNumber)=massMatrixCoefficient
    ENDIF
    
    EXITS("EquationsMappingVector_DynamicMatricesCoefficientsSet2")
    RETURN
999 ERRORSEXITS("EquationsMappingVector_DynamicMatricesCoefficientsSet2",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_DynamicMatricesCoefficientsSet2

  !
  !================================================================================================================================
  !

  !>Sets the mapping between a dependent field variable and the equations set dynamic matrices
  SUBROUTINE EquationsMappingVector_DynamicVariableTypeSet(vectorMapping,dynamicVariableType,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the vector equations mapping to set
    INTEGER(INTG), INTENT(IN) :: dynamicVariableType !<The variable type associated with the vector equations dynamic matrices.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx,numberOfLinearMatrixVariables,numberOfResidualVariables,residualIdx,residualVariableType, &
      & variableIdx
    LOGICAL :: isResidualType
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(FieldType), POINTER :: dependentField
    TYPE(FieldVariableType), POINTER :: dynamicVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsMappingVector_DynamicVariableTypeSet",err,error,*999)

    CALL EquationsMappingVector_AssertNotFinished(vectorMapping,err,error,*999)     
    NULLIFY(createValuesCache)
    CALL EquationsMappingVector_CreateValuesCacheGet(vectorMapping,createValuesCache,err,error,*999)
    NULLIFY(vectorEquations)
    CALL EquationsMappingVector_VectorEquationsGet(vectorMapping,vectorEquations,err,error,*999)
    NULLIFY(equations)
    CALL EquationsVector_EquationsGet(vectorEquations,equations,err,error,*999)
    CALL Equations_AssertIsDynamic(equations,err,error,*999)
   
    IF(dynamicVariableType==0) THEN
      createValuesCache%dynamicVariableType=0
    ELSE
      NULLIFY(equationsSet)
      CALL Equations_EquationsSetGet(equations,equationsSet,err,error,*999)      
      NULLIFY(dependentField)
      CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
      !Check the dynamic variable type is not being by other equations matrices or vectors
      IF(equations%linearity==EQUATIONS_NONLINEAR) THEN
        IF(createValuesCache%numberOfResiduals>0) THEN
          isResidualType=.FALSE.
          DO residualIdx=1,createValuesCache%numberOfResiduals
            CALL EquationsMappingVectorCVC_NumberOfResidualVariablesGet(createValuesCache,residualIdx,numberOfResidualVariables, &
              & err,error,*999)
            DO variableIdx=1,numberOfResidualVariables
              CALL EquationsMappingVectorCVC_ResidualVariableTypeGet(createValuesCache,variableIdx,residualIdx, &
                & residualVariableType,err,error,*999)
              IF(residualVariableType==dynamicVariableType) isResidualType=.TRUE.
            ENDDO !variableIdx
          ENDDO !residualIdx
          IF(.NOT.isResidualType) THEN
            localError="The specified dynamic variable type of "//TRIM(NumberToVString(dynamicVariableType,"*",err,error))// &
              & " is not the same as any residual variable type in dynamic nonlinear equations."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ENDIF
      ENDIF
      IF(createValuesCache%rhsVariableType==dynamicVariableType) THEN
        localError="The specified dynamic variable type of "//TRIM(NumberToVString(dynamicVariableType,"*",err,error))// &
          & " is the same as the variable type for the RHS vector."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      NULLIFY(dynamicVariable)
      CALL Field_VariableGet(dependentField,dynamicVariableType,dynamicVariable,err,error,*999)
      vectorMapping%createValuesCache%dynamicVariableType=dynamicVariableType
    ENDIF
    
    EXITS("EquationsMappingVector_DynamicVariableTypeSet")
    RETURN
999 ERRORSEXITS("EquationsMappingVector_DynamicVariableTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_DynamicVariableTypeSet
  
  !
  !================================================================================================================================
  !

  !>Finalises the vector equations mapping LHS mapping and deallocates all memory
  SUBROUTINE EquationsMappingVector_LHSMappingFinalise(lhsMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping !<A pointer to the LHS mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingVector_LHSMappingFinalise",err,error,*999)

    IF(ASSOCIATED(lhsMapping)) THEN
      IF(ASSOCIATED(lhsMapping%lhsDOFToEquationsRowMap)) DEALLOCATE(lhsMapping%lhsDOFToEquationsRowMap)
      IF(ASSOCIATED(lhsMapping%equationsRowToLHSDOFMap)) DEALLOCATE(lhsMapping%equationsRowToLHSDOFMap)
      DEALLOCATE(lhsMapping)
    ENDIF
       
    EXITS("EquationsMappingVector_LHSMappingFinalise")
    RETURN
999 ERRORSEXITS("EquationsMappingVector_LHSMappingFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_LHSMappingFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the equations mapping LHS mapping
  SUBROUTINE EquationsMappingVector_LHSMappingInitialise(vectorMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations mapping to initialise the LHS mapping for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("EquationsMappingVector_LHSMappingInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector equations mapping is not associated.",err,error,*998)
    IF(ASSOCIATED(vectorMapping%lhsMapping)) &
      & CALL FlagError("Vector equations mapping LHS mapping is already associated.",err,error,*998)
    
    ALLOCATE(vectorMapping%lhsMapping,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate equations mapping LHS mapping.",err,error,*999)
    vectorMapping%lhsMapping%vectorMapping=>vectorMapping        
    vectorMapping%lhsMapping%lhsVariableType=0
    NULLIFY(vectorMapping%lhsMapping%lhsVariable)
    vectorMapping%lhsMapping%numberOfRows=0
    vectorMapping%lhsMapping%totalNumberOfRows=0
    vectorMapping%lhsMapping%numberOfGlobalRows=0
    NULLIFY(vectorMapping%lhsMapping%rowDofsMapping)
    NULLIFY(vectorMapping%lhsMapping%lhsDOFToEquationsRowMap)
    NULLIFY(vectorMapping%lhsMapping%equationsRowToLHSDOFMap)
       
    EXITS("EquationsMappingVector_LHSMappingInitialise")
    RETURN
999 CALL EquationsMappingVector_LHSMappingFinalise(vectorMapping%lhsMapping,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMappingVector_LHSMappingInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_LHSMappingInitialise

  !
  !================================================================================================================================
  !

  !>Finalises the equations mapping linear mapping and deallocates all memory
  SUBROUTINE EquationsMappingVector_LinearMappingFinalise(linearMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping !<A pointer to the linear mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx,variableIdx
 
    ENTERS("EquationsMappingVector_LinearMappingFinalise",err,error,*999)

    IF(ASSOCIATED(linearMapping)) THEN
      IF(ALLOCATED(linearMapping%linearVariableTypes)) DEALLOCATE(linearMapping%linearVariableTypes)
      IF(ALLOCATED(linearMapping%linearVariableTypesMap)) DEALLOCATE(linearMapping%linearVariableTypesMap)
      IF(ALLOCATED(linearMapping%linearVariables)) DEALLOCATE(linearMapping%linearVariables)
      IF(ALLOCATED(linearMapping%varToEquationsMatricesMaps)) THEN
        DO variableIdx=1,SIZE(linearMapping%varToEquationsMatricesMaps,1)
          CALL EquationsMappingVToEMSMap_Finalise(linearMapping%varToEquationsMatricesMaps(variableIdx)%ptr,err,error,*999)
        ENDDO !variableIdx
        DEALLOCATE(linearMapping%varToEquationsMatricesMaps)        
      ENDIF
      IF(ALLOCATED(linearMapping%equationsMatrixToVarMaps)) THEN
        DO matrixIdx=1,SIZE(linearMapping%equationsMatrixToVarMaps,1)
          CALL EquationsMappingEMToVMap_Finalise(linearMapping%equationsMatrixToVarMaps(matrixIdx)%ptr,err,error,*999)
        ENDDO !matrixIdx
        DEALLOCATE(linearMapping%equationsMatrixToVarMaps)
      ENDIF
      IF(ALLOCATED(linearMapping%equationsRowToVariableDOFMaps)) DEALLOCATE(linearMapping%equationsRowToVariableDOFMaps)
      DEALLOCATE(linearMapping)
    ENDIF
    
    EXITS("EquationsMappingVector_LinearMappingFinalise")
    RETURN
999 ERRORSEXITS("EquationsMappingVector_LinearMappingFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_LinearMappingFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the equations mapping linear mapping
  SUBROUTINE EquationsMappingVector_LinearMappingInitialise(vectorMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations mapping to initialise the linear mapping for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("EquationsMappingVector_LinearMappingInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector equations mapping is not associated.",err,error,*998)
    IF(ASSOCIATED(vectorMapping%linearMapping)) &
      & CALL FlagError("Vector equations mapping linear mapping is already associated.",err,error,*998)
    
    ALLOCATE(vectorMapping%linearMapping,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate equations mapping linear mapping.",err,error,*999)
    vectorMapping%linearMapping%vectorMapping=>vectorMapping       
    vectorMapping%linearMapping%numberOfLinearMatrices=0
    vectorMapping%linearMapping%numberOfLinearVariables=0
       
    EXITS("EquationsMappingVector_LinearMappingInitialise")
    RETURN
999 CALL EquationsMappingVector_LinearMappingFinalise(vectorMapping%linearMapping,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMappingVector_LinearMappingInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_LinearMappingInitialise

  !
  !================================================================================================================================
  !

  !>Sets the coefficient for a linear equations matrix in an equation set. 
  SUBROUTINE EquationsMappingVector_LinearMatricesCoefficientsSet0(vectorMapping,linearMatrixCoefficient,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations mapping.
    REAL(DP), INTENT(IN) :: linearMatrixCoefficient !<The linear matrix coefficient
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(EquationsMappingVectorCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsMappingVector_LinearMatricesCoefficientsSet0",err,error,*999)

    CALL EquationsMappingVector_LinearMatricesCoefficientsSet1(vectorMapping,[linearMatrixCoefficient],err,error,*999)
    
    EXITS("EquationsMappingVector_LinearMatricesCoefficientsSet0")
    RETURN
999 ERRORSEXITS("EquationsMappingVector_LinearMatricesCoefficientsSet0",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_LinearMatricesCoefficientsSet0

  !
  !================================================================================================================================
  !

  !>Sets the coefficients for the linear equations matrices in an equation set. 
  SUBROUTINE EquationsMappingVector_LinearMatricesCoefficientsSet1(vectorMapping,linearMatricesCoefficients,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations mapping.
    REAL(DP), INTENT(IN) :: linearMatricesCoefficients(:) !<linearMatricesCoefficients(matrixIdx). The linear matrices coefficients
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(EquationsMappingVectorCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsMappingVector_LinearMatricesCoefficientsSet1",err,error,*999)

    CALL EquationsMappingVector_AssertNotFinished(vectorMapping,err,error,*999)
    NULLIFY(createValuesCache)
    CALL EquationsMappingVector_CreateValuesCacheGet(vectorMapping,createValuesCache,err,error,*999)
    IF(SIZE(linearMatricesCoefficients,1)<createValuesCache%numberOfLinearMatrices) THEN
      localError="The size of the specified linear matrices coefficients array of "// &
        & TRIM(NumberToVString(SIZE(linearMatricesCoefficients,1),"*",err,error))// &
        & " is invalid. The size of linear matrix coefficients array should be >= "// &
        & TRIM(NumberToVString(vectorMapping%createValuesCache%numberOfLinearMatrices,"*",err,error))//")."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    createValuesCache%linearMatrixCoefficients(1:createValuesCache%numberOfLinearMatrices)= &
      & linearMatricesCoefficients(1:createValuesCache%numberOfLinearMatrices)
    
    EXITS("EquationsMappingVector_LinearMatricesCoefficientsSet1")
    RETURN
999 ERRORSEXITS("EquationsMappingVector_LinearMatricesCoefficientsSet1",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_LinearMatricesCoefficientsSet1

  !
  !================================================================================================================================
  !

  !>Sets the mapping between a dependent field variable type and a linear equations matrix
  SUBROUTINE EquationsMappingVector_LinearMatricesVariableTypesSet0(vectorMapping,linearMatrixVariableType,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations mapping
    INTEGER(INTG), INTENT(IN) :: linearMatrixVariableType !<The matrix variable typs to map to the linear equations matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("EquationsMappingVector_LinearMatricesVariableTypesSet0",err,error,*999)

    CALL EquationsMappingVector_LinearMatricesVariableTypesSet1(vectorMapping,[linearMatrixVariableType],err,error,*999)
    
    EXITS("EquationsMappingVector_LinearMatricesVariableTypesSet0")
    RETURN
999 ERRORS("EquationsMappingVector_LinearMatricesVariableTypesSet0",err,error)
    EXITS("EquationsMappingVector_LinearMatricesVariableTypesSet0")
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_LinearMatricesVariableTypesSet0

  !
  !================================================================================================================================
  !

  !>Sets the mapping between the dependent field variable types and the linear equations matrices
  SUBROUTINE EquationsMappingVector_LinearMatricesVariableTypesSet1(vectorMapping,linearMatricesVariableTypes,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations mapping
    INTEGER(INTG), INTENT(IN) :: linearMatricesVariableTypes(:) !<linearMatricesVariableTypes(matrixIdx). The matrix variable types to map to each linear equations matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(FieldType), POINTER :: dependentField
    TYPE(FieldVariableType), POINTER :: linearVariable
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("EquationsMappingVector_LinearMatricesVariableTypesSet1",err,error,*999)

    CALL EquationsMappingVector_AssertNotFinished(vectorMapping,err,error,*999)
    NULLIFY(createValuesCache)
    CALL EquationsMappingVector_CreateValuesCacheGet(vectorMapping,createValuesCache,err,error,*999)
    NULLIFY(vectorEquations)
    CALL EquationsMappingVector_VectorEquationsGet(vectorMapping,vectorEquations,err,error,*999)
    NULLIFY(equations)
    CALL EquationsVector_EquationsGet(vectorEquations,equations,err,error,*999)
    NULLIFY(equationsSet)
    CALL Equations_EquationsSetGet(equations,equationsSet,err,error,*999)
    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)    
    !Check input values
    IF(SIZE(linearMatricesVariableTypes,1)<createValuesCache%numberOfLinearMatrices) THEN
      localError="The size of supplied linear matrix variable types of "// &
        & TRIM(NumberToVString(SIZE(linearMatricesVariableTypes,1),"*",err,error))// &
        & " is invalid. The size must be >= "// &
        & TRIM(NumberToVString(vectorMapping%createValuesCache%numberOfLinearMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    DO matrixIdx=1,createValuesCache%numberOfLinearMatrices
      IF(linearMatricesVariableTypes(matrixIdx)==0) THEN
        localError="The specified linear matrix variable type of "// &
          & TRIM(NumberToVString(linearMatricesVariableTypes(matrixIdx),"*",err,error))// &
          & " for linear matrix number "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
          & " is invalid. The variable type should be >= 1 and <= "// &
          & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      !Check the linear variable type is not the same as the RHS variable type
      IF(createValuesCache%rhsVariableType==linearMatricesVariableTypes(matrixIdx)) THEN
        localError="The specified linear matrix variable type of "// &
          & TRIM(NumberToVString(linearMatricesVariableTypes(matrixIdx),"*",err,error))// &
          & " for linear matrix number "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
          & " is the same as the variable type for the RHS vector."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      !Check to see if the linear matrix variable numbers are defined on the dependent field
      NULLIFY(linearVariable)
      CALL Field_VariableExists(dependentField,linearMatricesVariableTypes(matrixIdx),linearVariable,err,error,*999)
      IF(.NOT.ASSOCIATED(linearVariable)) THEN
        localError="The linear matrix variable type of "// &
          & TRIM(NumberToVString(linearMatricesVariableTypes(matrixIdx),"*",err,error))// &
            & " for linear matrix number "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
            & " is not defined in the dependent field."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDDO !matrixIdx
    
    createValuesCache%linearMatrixVariableTypes(1:createValuesCache%numberOfLinearMatrices)= &
      & linearMatricesVariableTypes(1:createValuesCache%numberOfLinearMatrices)
    
    EXITS("EquationsMappingVector_LinearMatricesVariableTypesSet1")
    RETURN
999 ERRORS("EquationsMappingVector_LinearMatricesVariableTypesSet1",err,error)
    EXITS("EquationsMappingVector_LinearMatricesVariableTypesSet1")
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_LinearMatricesVariableTypesSet1

  !
  !================================================================================================================================
  !

  !>Finalises the equations mapping nonlinear mapping and deallocates all memory
  SUBROUTINE EquationsMappingVector_NonlinearMappingFinalise(nonlinearMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping !<A pointer to the nonlinear mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: residualIdx
 
    ENTERS("EquationsMappingVector_NonlinearMappingFinalise",err,error,*999)

    IF(ASSOCIATED(nonlinearMapping)) THEN
      IF(ALLOCATED(nonlinearMapping%nonlinearVariableTypes)) DEALLOCATE(nonlinearMapping%nonlinearVariableTypes)
      IF(ALLOCATED(nonlinearMapping%nonlinearVariableTypesMap)) DEALLOCATE(nonlinearMapping%nonlinearVariableTypesMap)
      IF(ALLOCATED(nonlinearMapping%nonlinearVariables)) DEALLOCATE(nonlinearMapping%nonlinearVariables)
      IF(ALLOCATED(nonlinearMapping%residuals)) THEN
        DO residualIdx=1,SIZE(nonlinearMapping%residuals)
          CALL EquationsMappingNonlinear_ResidualMappingFinalise(nonlinearMapping%residuals(residualIdx)%ptr,err,error,*999)
        ENDDO !residualIdx
      ENDIF
    ENDIF
    
    EXITS("EquationsMappingVector_NonlinearMappingFinalise")
    RETURN
999 ERRORSEXITS("EquationsMappingVector_NonlinearMappingFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_NonlinearMappingFinalise
  
  !
  !================================================================================================================================
  !

  !>Initialises the equations mapping nonlinear mapping
  SUBROUTINE EquationsMappingVector_NonlinearMappingInitialise(vectorMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations mapping to initialise the nonlinear mapping for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("EquationsMappingVector_NonlinearMappingInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector equations mapping is not associated.",err,error,*998)
    IF(ASSOCIATED(vectorMapping%nonlinearMapping)) &
      & CALL FlagError("Equations mapping nonlinear mapping is already associated.",err,error,*998)
     
    ALLOCATE(vectorMapping%nonlinearMapping,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate equations mapping nonlinear mapping.",err,error,*999)
    vectorMapping%nonlinearMapping%vectorMapping=>vectorMapping
    vectorMapping%nonlinearMapping%numberOfResiduals=1
    vectorMapping%nonlinearMapping%numberOfNonlinearVariables=1

    EXITS("EquationsMappingVector_NonlinearMappingInitialise")
    RETURN
999 CALL EquationsMappingVector_NonlinearMappingFinalise(vectorMapping%nonlinearMapping,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMappingVector_NonlinearMappingInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_NonlinearMappingInitialise

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of linear equations matrices
  SUBROUTINE EquationsMappingVector_NumberOfLinearMatricesSet(vectorMapping,numberOfLinearMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations mapping to set the number of matrices for.
    INTEGER(INTG), INTENT(IN) :: numberOfLinearMatrices !<The number of linear equations matrices for the mapping.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: linearVariableType,matrixIdx,minNumberOfLinearMatrices,numberOfVariables,variableIdx
    INTEGER(INTG), ALLOCATABLE :: newLinearMatrixVariableTypes(:)
    REAL(DP), ALLOCATABLE :: newLinearMatrixCoefficients(:)
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(FieldType), POINTER :: dependentField
    TYPE(FieldVariableType), POINTER :: linearVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsMappingVector_NumberOfLinearMatricesSet",err,error,*999)

    CALL EquationsMappingVector_AssertNotFinished(vectorMapping,err,error,*999)
    NULLIFY(createValuesCache)
    CALL EquationsMappingVector_CreateValuesCacheGet(vectorMapping,createValuesCache,err,error,*999)
    NULLIFY(vectorEquations)
    CALL EquationsMappingVector_VectorEquationsGet(vectorMapping,vectorEquations,err,error,*999)
    NULLIFY(equations)
    CALL EquationsVector_EquationsGet(vectorEquations,equations,err,error,*999)

    !Check number of matrices to create is valid
    SELECT CASE(equations%timeDependence)
    CASE(EQUATIONS_STATIC,EQUATIONS_QUASISTATIC)
      SELECT CASE(equations%linearity)
      CASE(EQUATIONS_LINEAR,EQUATIONS_NONLINEAR_BCS)
        !Static, linear equations
        IF(numberOfLinearMatrices<1) THEN
          localError="The specified number of linear matrices of "//TRIM(NumberToVString(numberOfLinearMatrices,"*", &
            & err,error))//" is invalid for static linear problems. The number of linear matrices should be >= 1."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      CASE(EQUATIONS_NONLINEAR)
        IF(numberOfLinearMatrices<0) THEN
          localError="The specified number of linear matrices of "//TRIM(NumberToVString(numberOfLinearMatrices,"*", &
            & err,error))//" is invalid for static nonlinear problems. The number of linear matrices should be >= 0."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      CASE DEFAULT
        localError="The equations linearity type of "//TRIM(NumberToVString(equations%linearity,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_FIRST_ORDER_DYNAMIC,EQUATIONS_SECOND_ORDER_DYNAMIC)
      SELECT CASE(equations%linearity)
      CASE(EQUATIONS_LINEAR,EQUATIONS_NONLINEAR_BCS)
        IF(numberOfLinearMatrices<0) THEN
          localError="The specified number of linear matrices of "//TRIM(NumberToVString(numberOfLinearMatrices,"*", &
            & err,error))//" is invalid for dynamic linear problems. The number of linear matrices should be >= 0."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      CASE(EQUATIONS_NONLINEAR)
        IF(numberOfLinearMatrices<0) THEN
          localError="The specified number of linear matrices of "//TRIM(NumberToVString(numberOfLinearMatrices,"*", &
            & err,error))//" is invalid for dynamic nonlinear problems. The number of linear matrices should be >= 0."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      CASE DEFAULT
        localError="The equations linearity type of "//TRIM(NumberToVString(equations%linearity,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The equations time dependence type of "//TRIM(NumberToVString(equations%timeDependence,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    !See if we have any linear matrices left.
    IF(numberOfLinearMatrices==0) THEN
      IF(createValuesCache%numberOfLinearMatrices>0) THEN
        IF(ALLOCATED(createValuesCache%linearMatrixVariableTypes)) DEALLOCATE(createValuesCache%linearMatrixVariableTypes)
        IF(ALLOCATED(createValuesCache%linearMatrixCoefficients)) DEALLOCATE(createValuesCache%linearMatrixCoefficients)        
      ENDIF
      createValuesCache%numberOfLinearMatrices=0
    ELSE
      !If we need to reallocate and reset all the create values cache arrays and change the number of matrices
      IF(numberOfLinearMatrices/=createValuesCache%numberOfLinearMatrices) THEN
        ALLOCATE(newLinearMatrixVariableTypes(numberOfLinearMatrices),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate new linear matrix variable types.",err,error,*999)
        ALLOCATE(newLinearMatrixCoefficients(numberOfLinearMatrices),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate old linear matrix coefficients.",err,error,*999)
        newLinearMatrixVariableTypes=0
        newLinearMatrixCoefficients=1.0_DP
        minNumberOfLinearMatrices=MIN(numberOfLinearMatrices,createValuesCache%numberOfLinearMatrices)
        newLinearMatrixVariableTypes(1:minNumberOfLinearMatrices)= &
          & createValuesCache%linearMatrixVariableTypes(1:minNumberOfLinearMatrices)
        newLinearMatrixCoefficients(1:minNumberOfLinearMatrices)= &
          & createValuesCache%linearMatrixCoefficients(1:numberOfLinearMatrices)
        IF(numberOfLinearMatrices>createValuesCache%numberOfLinearMatrices) THEN
          linearVariableType=0
          IF(createValuesCache%numberOfLinearMatrices==0) THEN
            !Set the linear variable type to the first variable that is not the rhs variable
            NULLIFY(equationsSet)
            CALL Equations_EquationsSetGet(equations,equationsSet,err,error,*999)
            NULLIFY(dependentField)
            CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
            CALL Field_NumberOfVariablesGet(dependentField,numberOfVariables,err,error,*999)
            DO variableIdx=1,numberOfVariables
              NULLIFY(linearVariable)
              CALL Field_VariableIndexGet(dependentField,1,linearVariable,linearVariableType,err,error,*999)
              IF(linearVariableType/=createValuesCache%rhsVariableType) EXIT
            ENDDO !variableIdx
          ELSE
            linearVariableType=createValuesCache%linearMatrixVariableTypes(1)
          ENDIF
          IF(linearVariableType==0) CALL FlagError("Could not find a linear variable type.",err,error,*999)
          DO matrixIdx=minNumberOfLinearMatrices+1,numberOfLinearMatrices
            newLinearMatrixVariableTypes(matrixIdx)=linearVariableType
          ENDDO !matrixIdx
        ENDIF
        CALL MOVE_ALLOC(newLinearMatrixVariableTypes,createValuesCache%linearMatrixVariableTypes)
        CALL MOVE_ALLOC(newLinearMatrixCoefficients,createValuesCache%linearMatrixCoefficients)
        createValuesCache%numberOfLinearMatrices=numberOfLinearMatrices
      ENDIF
    ENDIF
    
    EXITS("EquationsMappingVector_NumberOfLinearMatricesSet")
    RETURN
999 IF(ALLOCATED(newLinearMatrixVariableTypes)) DEALLOCATE(newLinearMatrixVariableTypes)    
    IF(ALLOCATED(newLinearMatrixCoefficients)) DEALLOCATE(newLinearMatrixCoefficients)    
    ERRORSEXITS("EquationsMappingVector_NumberOfLinearMatricesSet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_NumberOfLinearMatricesSet

  !
  !================================================================================================================================
  !

  !>Sets the number of residuals in the vector equations mapping.
  SUBROUTINE EquationsMappingVector_NumberOfResidualsSet(vectorMapping,numberOfResiduals,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the vector equations mapping to set
    INTEGER(INTG), INTENT(IN) :: numberOfResiduals !<The number of residual variables for this vector equations mapping.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: previousNumberOfResiduals,residualIdx,variableIdx
    INTEGER(INTG), ALLOCATABLE :: newNumberOfResidualVariables(:),newResidualVariableTypes(:,:)
    REAL(DP), ALLOCATABLE :: newResidualCoefficients(:)
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsMappingVectorCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsMappingVector_NumberOfResidualsSet",err,error,*999)

    CALL EquationsMappingVector_AssertNotFinished(vectorMapping,err,error,*999)
    NULLIFY(createValuesCache)
    CALL EquationsMappingVector_CreateValuesCacheGet(vectorMapping,createValuesCache,err,error,*999)
    NULLIFY(vectorEquations)
    CALL EquationsMappingVector_VectorEquationsGet(vectorMapping,vectorEquations,err,error,*999)
    NULLIFY(equations)
    CALL EquationsVector_EquationsGet(vectorEquations,equations,err,error,*999)

    !Check that the number of residuals is valid
    SELECT CASE(equations%timeDependence)
    CASE(EQUATIONS_STATIC,EQUATIONS_QUASISTATIC)
      SELECT CASE(equations%linearity)
      CASE(EQUATIONS_LINEAR,EQUATIONS_NONLINEAR_BCS)
        !Linear equations can't have any residuals
        IF(numberOfResiduals/=0) THEN
          localError="The specified number of residuals of "//TRIM(NumberToVString(numberOfResiduals,"*",err,error))// &
            & " is invalid for static linear equations. The number of residuals should be 0."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      CASE(EQUATIONS_NONLINEAR)
        !Must have at least one residual for a nonlinear equations set
        IF(numberOfResiduals<1) THEN
          localError="The specified number of residuals of "//TRIM(NumberToVString(numberOfResiduals,"*",err,error))// &
            & " is invalid for static nonlinear equations. The number of residuals should be >= 1."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      CASE DEFAULT
        localError="The equations linearity type of "//TRIM(NumberToVString(equations%linearity,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_FIRST_ORDER_DYNAMIC,EQUATIONS_SECOND_ORDER_DYNAMIC)
      SELECT CASE(equations%linearity)
      CASE(EQUATIONS_LINEAR,EQUATIONS_NONLINEAR_BCS)
        !Linear equations can't have any residuals
        IF(numberOfResiduals/=0) THEN
          localError="The specified number of residuals of "//TRIM(NumberToVString(numberOfResiduals,"*",err,error))// &
            & " is invalid for dynamic linear equations. The number of residuals should be 0."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      CASE(EQUATIONS_NONLINEAR)
        !Must have at least one residual for a nonlinear equations set
        IF(numberOfResiduals<1) THEN
          localError="The specified number of residuals of "//TRIM(NumberToVString(numberOfResiduals,"*",err,error))// &
            & " is invalid for dynamic nonlinear equations. The number of residuals should be >= 1."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      CASE DEFAULT
        localError="The equations linearity type of "//TRIM(NumberToVString(equations%linearity,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The equations time dependence type of "//TRIM(NumberToVString(equations%timeDependence,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
     
    IF(numberOfResiduals/=createValuesCache%numberOfResiduals) THEN
      !Create new number of residual variables, residual variable types array and coeficients arrays
      ALLOCATE(newNumberOfResidualVariables(numberOfResiduals),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate new number of residuals variables array.",err,error,*999)
      ALLOCATE(newResidualVariableTypes(SIZE(createValuesCache%residualVariableTypes,1),numberOfResiduals),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate new residual variables types array.",err,error,*999)
      ALLOCATE(newResidualCoefficients(numberOfResiduals),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate new residual coefficients array.",err,error,*999)
      newNumberOfResidualVariables=0
      newResidualVariableTypes=0
      newResidualCoefficients=1.0_DP
      IF(previousNumberOfResiduals>0) THEN
        DO residualIdx=1,MIN(numberOfResiduals,previousNumberOfResiduals)
          newNumberOfResidualVariables(residualIdx)=createValuesCache%numberOfResidualVariables(residualIdx)
          DO variableIdx=1,createValuesCache%numberOfResidualVariables(residualIdx)
            newResidualVariableTypes(variableIdx,residualIdx)=createValuesCache%residualVariableTypes(variableIdx,residualIdx)
          ENDDO !variableIdx
          newResidualCoefficients(residualIdx)=createValuesCache%residualCoefficients(residualIdx)
        ENDDO !residualIdx
        !Default any new residuals to be the same as the first residual
        DO residualIdx=createValuesCache%numberOfResiduals+1,numberOfResiduals
          newNumberOfResidualVariables(residualIdx)=createValuesCache%numberOfResidualVariables(1)
          DO variableIdx=1,createValuesCache%numberOfResidualVariables(1)
            newResidualVariableTypes(variableIdx,residualIdx)=createValuesCache%residualVariableTypes(variableIdx,1)
          ENDDO !variableIdx
        ENDDO !residualIdx
      ENDIF
      CALL MOVE_ALLOC(newNumberOfResidualVariables,createValuesCache%numberOfResidualVariables)
      CALL MOVE_ALLOC(newResidualVariableTypes,createValuesCache%residualVariableTypes)
      CALL MOVE_ALLOC(newResidualCoefficients,createValuesCache%residualCoefficients)
      !Set number of residual variables
      createValuesCache%numberOfResiduals=numberOfResiduals
    ENDIF

    EXITS("EquationsMappingVector_NumberOfResidualsSet")
    RETURN
999 IF(ALLOCATED(newNumberOfResidualVariables)) DEALLOCATE(newNumberOfResidualVariables)
    IF(ALLOCATED(newResidualVariableTypes)) DEALLOCATE(newResidualVariableTypes)
    IF(ALLOCATED(newResidualCoefficients)) DEALLOCATE(newResidualCoefficients)
    ERRORSEXITS("EquationsMappingVector_NumberOfResidualsSet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_NumberOfResidualsSet

  !
  !================================================================================================================================
  !

  !>Sets the number of sources in the vector equations mapping.
  SUBROUTINE EquationsMappingVector_NumberOfSourcesSet(vectorMapping,numberOfSources,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the vector equations mapping to set
    INTEGER(INTG), INTENT(IN) :: numberOfSources !<The number of sources for this vector equations mapping.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: sourceIdx,sourceVariableType
    INTEGER(INTG), ALLOCATABLE :: newSourceVariableTypes(:)
    REAL(DP), ALLOCATABLE :: newSourceCoefficients(:)
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsMappingVectorCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(FieldType), POINTER :: sourceField
    TYPE(FieldVariableType), POINTER :: sourceVariable    
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsMappingVector_NumberOfSourcesSet",err,error,*999)

    CALL EquationsMappingVector_AssertNotFinished(vectorMapping,err,error,*999)
    NULLIFY(createValuesCache)
    CALL EquationsMappingVector_CreateValuesCacheGet(vectorMapping,createValuesCache,err,error,*999)
    NULLIFY(vectorEquations)
    CALL EquationsMappingVector_VectorEquationsGet(vectorMapping,vectorEquations,err,error,*999)
    NULLIFY(equations)
    CALL EquationsVector_EquationsGet(vectorEquations,equations,err,error,*999)
    NULLIFY(equationsSet)
    CALL Equations_EquationsSetGet(equations,equationsSet,err,error,*999)
    CALL EquationsSet_AssertSourceIsFinished(equationsSet,err,error,*999)
    NULLIFY(sourceField)
    CALL EquationsSet_SourceFieldExists(equationsSet,sourceField,err,error,*999)
    IF(.NOT.ASSOCIATED(sourceField)) THEN
      localError="The source field does not exist for equations set "// &
        & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))//" does not exist."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    !Check that the number of sources is valid
    IF(numberOfSources<0) THEN
      localError="The specified number of sources of "//TRIM(NumberToVString(numberOfSources,"*",err,error))// &
        & " is invalid. The number of sources should be >= 0."
      CALL FlagError(localError,err,error,*999)
    ENDIF
      
    IF(numberOfSources/=createValuesCache%numberOfSources) THEN
      IF(numberOfSources==0) THEN
        IF(ALLOCATED(createValuesCache%sourceVariableTypes)) DEALLOCATE(createValuesCache%sourceVariableTypes)
        IF(ALLOCATED(createValuesCache%sourceCoefficients)) DEALLOCATE(createValuesCache%sourceCoefficients)
      ELSE
        !Create new source variable types array and coeficients arrays
        ALLOCATE(newSourceVariableTypes(numberOfSources),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate new source variable types array.",err,error,*999)
        ALLOCATE(newSourceCoefficients(numberOfSources),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate new source coefficients array.",err,error,*999)
        newSourceVariableTypes=0
        newSourceCoefficients=1.0_DP
        DO sourceIdx=1,MIN(numberOfSources,createValuesCache%numberOfSources)
          newSourceVariableTypes(sourceIdx)=createValuesCache%sourceVariableTypes(sourceIdx)
          newSourceCoefficients(sourceIdx)=createValuesCache%sourceCoefficients(sourceIdx)
        ENDDO !sourceIdx
        IF(ALLOCATED(createValuesCache%sourceVariableTypes)) THEN
          !Default any new sources to be the previous first source vriable type
          sourceVariableType=createValuesCache%sourceVariableTypes(1)                    
        ELSE 
          !Default any new sources to be the same as the first source field variable
          NULLIFY(sourceVariable)
          CALL Field_VariableIndexGet(sourceField,1,sourceVariable,sourceVariableType,err,error,*999)
        ENDIF
        DO sourceIdx=createValuesCache%numberOfSources+1,numberOfSources
          newSourceVariableTypes(sourceIdx)=sourceVariableType
        ENDDO !sourceIdx
        CALL MOVE_ALLOC(newSourceVariableTypes,createValuesCache%sourceVariableTypes)
        CALL MOVE_ALLOC(newSourceCoefficients,createValuesCache%sourceCoefficients)
      ENDIF
      !Set number of source variables
      createValuesCache%numberOfSources=numberOfSources
    ENDIF

    EXITS("EquationsMappingVector_NumberOfSourcesSet")
    RETURN
999 IF(ALLOCATED(newSourceVariableTypes)) DEALLOCATE(newSourceVariableTypes)
    IF(ALLOCATED(newSourceCoefficients)) DEALLOCATE(newSourceCoefficients)
    ERRORSEXITS("EquationsMappingVector_NumberOfSourcesSet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_NumberOfSourcesSet

  !
  !================================================================================================================================
  !

  !>Sets the coefficient applied to the equations set residual vector.
  SUBROUTINE EquationsMappingVector_ResidualCoefficientsSet0(vectorMapping,residualCoefficient,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations vector mapping to set   
    REAL(DP), INTENT(IN) :: residualCoefficient !<The coefficient applied to the residual vector.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsMappingVector_ResidualCoefficientsSet0",err,error,*999)
    
    CALL EquationsMappingVector_ResidualCoefficientsSet1(vectorMapping,[residualCoefficient],err,error,*999)
      
    EXITS("EquationsMappingVector_ResidualCoefficientsSet0")
    RETURN
999 ERRORSEXITS("EquationsMappingVector_ResidualCoefficientsSet0",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_ResidualCoefficientsSet0

  !
  !================================================================================================================================
  !

  !>Sets the coefficient applied to the equations set residual vector.
  SUBROUTINE EquationsMappingVector_ResidualCoefficientsSet1(vectorMapping,residualCoefficients,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations mapping to set   
    REAL(DP), INTENT(IN) :: residualCoefficients(:) !<The coefficients applied to the residual vectors.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(EquationsMappingVectorCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsMappingVector_ResidualCoefficientsSet1",err,error,*999)
    
    CALL EquationsMappingVector_AssertNotFinished(vectorMapping,err,error,*999)
    NULLIFY(createValuesCache)
    CALL EquationsMappingVector_CreateValuesCacheGet(vectorMapping,createValuesCache,err,error,*999)
    IF(SIZE(residualCoefficients,1)<createValuesCache%numberOfResiduals) THEN
      localError="The size of the specified residual coefficients array of "// &
        & TRIM(NumberToVString(SIZE(residualCoefficients,1),"*",err,error))// &
        & " is invalid. The size of residual coefficients array should be >= "// &
        & TRIM(NumberToVString(createValuesCache%numberOfResiduals,"*",err,error))//")."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    createValuesCache%residualCoefficients(1:createValuesCache%numberOfResiduals)= &
      & residualCoefficients(1:createValuesCache%numberOfResiduals)
      
    EXITS("EquationsMappingVector_ResidualCoefficientsSet1")
    RETURN
999 ERRORSEXITS("EquationsMappingVector_ResidualCoefficientsSet1",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_ResidualCoefficientsSet1

  !
  !================================================================================================================================
  !

  !>Sets the number of dependent field variable that are mapped to a specified vector equations residual vector.
  SUBROUTINE EquationsMappingVector_ResidualNumberOfVariablesSet(vectorMapping,residualIdx,numberOfVariables,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the vector equations mapping to set the number of residual variables for
    INTEGER(INTG), INTENT(IN) :: residualIdx !<The index of the residual to set the number of variables for
    INTEGER(INTG), INTENT(IN) :: numberOfVariables !<The number of residual variables for the specified residual.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: maxNumberOfVariables,minNumberOfVariables,previousNumberOfVariables
    INTEGER(INTG), ALLOCATABLE :: newResidualVariableTypes(:,:)
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsMappingVectorCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsMappingVector_ResidualNumberOfVariablesSet",err,error,*999)

    CALL EquationsMappingVector_AssertNotFinished(vectorMapping,err,error,*999)
    NULLIFY(vectorEquations)
    CALL EquationsMappingVector_VectorEquationsGet(vectorMapping,vectorEquations,err,error,*999)
    NULLIFY(equations)
    CALL EquationsVector_EquationsGet(vectorEquations,equations,err,error,*999)
    CALL Equations_AssertIsNonlinear(equations,err,error,*999)
    IF(numberOfVariables<1) THEN
      localError="The specified number of variables of "//TRIM(NumberToVString(numberOfVariables,"*",err,error))// &
        & " is invalid. The number of variables for a residual should be >= 1."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    NULLIFY(createValuesCache)
    CALL EquationsMappingVector_CreateValuesCacheGet(vectorMapping,createValuesCache,err,error,*999)
    CALL EquationsMappingVectorCVC_NumberOfResidualVariablesGet(createValuesCache,residualIdx,previousNumberOfVariables, &
      & err,error,*999)
    IF(numberOfVariables/=previousNumberOfVariables) THEN
      !Set number of residual variables
      createValuesCache%numberOfResidualVariables(residualIdx)=numberOfVariables
      !Create new residual variable types array if required and copy over previous values
      maxNumberOfVariables=MAXVAL(createValuesCache%numberOfResidualVariables)
      IF(maxNumberOfVariables>SIZE(createValuesCache%residualVariableTypes,1)) THEN
        minNumberOfVariables=MIN(maxNumberOfVariables,SIZE(createValuesCache%residualVariableTypes,1))
        ALLOCATE(newResidualVariableTypes(maxNumberOfVariables,createValuesCache%numberOfResiduals),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate new residual variable types.",err,error,*999)
        newResidualVariableTypes=0
        newResidualVariableTypes(1:minNumberOfVariables,1:createValuesCache%numberOfResiduals)= &
          & createValuesCache%residualVariableTypes(1:minNumberOfVariables,1:createValuesCache%numberOfResiduals)
        !Set extra variables to be the same as the first variable
        newResidualVariableTypes(minNumberOfVariables+1:numberOfVariables,residualIdx)= &
          & createValuesCache%residualVariableTypes(1,residualIdx)
        CALL MOVE_ALLOC(newResidualVariableTypes,createValuesCache%residualVariableTypes)
      ENDIF
    ENDIF

    EXITS("EquationsMappingVector_ResidualNumberOfVariablesSet")
    RETURN
999 IF(ALLOCATED(newResidualVariableTypes)) DEALLOCATE(newResidualVariableTypes)
    ERRORSEXITS("EquationsMappingVector_ResidualNumberOfVariablesSet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_ResidualNumberOfVariablesSet

  !
  !================================================================================================================================
  !

  !>Sets the mapping between a dependent field variable and the equations residual vector.
  SUBROUTINE EquationsMappingVector_ResidualVariableTypesSet0(vectorMapping,residualIdx,residualVariableType,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations mapping to set
    INTEGER(INTG), INTENT(IN) :: residualIdx !<The residual index to set the variable types for.
    INTEGER(INTG), INTENT(IN) :: residualVariableType !<The variable type associated with the residualIdx'th residual vector. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingVector_ResidualVariableTypesSet0",err,error,*999)

    CALL EquationsMappingVector_ResidualVariableTypesSet1(vectorMapping,residualIdx,[residualVariableType],err,error,*999) 
      
    EXITS("EquationsMappingVector_ResidualVariableTypesSet0")
    RETURN
999 ERRORSEXITS("EquationsMappingVector_ResidualVariableTypesSet0",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_ResidualVariableTypesSet0
  
  !
  !================================================================================================================================
  !

  !>Sets the mapping between a dependent field variable and the equations residual vector.
  SUBROUTINE EquationsMappingVector_ResidualVariableTypesSet1(vectorMapping,residualIdx,residualVariableTypes,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations mapping to set
    INTEGER(INTG), INTENT(IN) :: residualIdx !<The residual index to set the variable types for.
    INTEGER(INTG), INTENT(IN) :: residualVariableTypes(:) !<residualVariableTypes(variableIdx). The variableIdx'th variable type associated with the residualIdx'th residual vector. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx,variableIdx,numberOfResidualVariables,residualVariableType
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(FieldType), POINTER :: dependentField
    TYPE(FieldVariableType), POINTER :: residualVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsMappingVector_ResidualVariableTypesSet1",err,error,*999)

    CALL EquationsMappingVector_AssertNotFinished(vectorMapping,err,error,*999)
    NULLIFY(vectorEquations)
    CALL EquationsMappingVector_VectorEquationsGet(vectorMapping,vectorEquations,err,error,*999)
    NULLIFY(equations)
    CALL EquationsVector_EquationsGet(vectorEquations,equations,err,error,*999)
    CALL Equations_AssertIsNonlinear(equations,err,error,*999)   
    NULLIFY(createValuesCache)
    CALL EquationsMappingVector_CreateValuesCacheGet(vectorMapping,createValuesCache,err,error,*999)
    NULLIFY(equationsSet)
    CALL Equations_EquationsSetGet(equations,equationsSet,err,error,*999)
    IF(residualIdx<1.OR.residualIdx>createValuesCache%numberOfResiduals) THEN
      localError="The specified residual index of "//TRIM(NumberToVString(residualIdx,"*",err,error))// &
        & " is invalid. The residual index should be >= 1 and <= "// &
        & TRIM(NumberToVString(createValuesCache%numberOfResiduals,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF    
    numberOfResidualVariables=SIZE(residualVariableTypes,1)
    IF(numberOfResidualVariables/=createValuesCache%numberOfResidualVariables(residualIdx)) THEN
      localError="Invalid number of variables. The number of residual variables " &
        & //TRIM(NumberToVString(numberOfResidualVariables,"*",err,error))// &
        & " for residual index "//TRIM(NumberToVString(residualIdx,"*",err,error))//" should be "// &
        & TRIM(NumberToVString(createValuesCache%numberOfResidualVariables(residualIdx),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)      
    !Check the residual variable types
    DO variableIdx=1,numberOfResidualVariables
      residualVariableType=residualVariableTypes(variableIdx)
      IF(createValuesCache%rhsVariableType==residualVariableType) THEN
        localError="The specified residual variable type of "//TRIM(NumberToVString(residualVariableType,"*",err,error))// &
          & " is the same as the variable type for the RHS vector."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      !Check the residual variable number is defined on the dependent field
      IF(residualVariableType<1.OR.residualVariableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
        localError="The specified residual variable type of "//TRIM(NumberToVString(residualVariableType,"*",err,error))// &
          & " is invalid. The variable type must either be zero or >= 1 and <= "// &
          & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      NULLIFY(residualVariable)
      CALL Field_VariableExists(dependentField,residualVariableType,residualVariable,err,error,*999)
      IF(.NOT.ASSOCIATED(residualVariable)) THEN
        localError="The specified residual variable type of "//TRIM(NumberToVString(residualVariableType,"*",err,error))// &
          & " is not defined on the dependent field."
        CALL FlagError(localError,err,error,*999)
      ENDIF      
    ENDDO !variableIdx
    
    createValuesCache%residualVariableTypes(1:numberOfResidualVariables,residualIdx)= &
      & residualVariableTypes(1:numberOfResidualVariables)
      
    EXITS("EquationsMappingVector_ResidualVariableTypesSet1")
    RETURN
999 ERRORSEXITS("EquationsMappingVector_ResidualVariableTypesSet1",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_ResidualVariableTypesSet1
  
  !
  !================================================================================================================================
  !

  !>Sets the coefficient applied to the equations set RHS vector.
  SUBROUTINE EquationsMappingVector_RHSCoefficientSet(vectorMapping,rhsCoefficient,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations mapping to set
    REAL(DP), INTENT(IN) :: rhsCoefficient!<The coefficient applied to the equations set RHS vector.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(EquationsMappingVectorCreateValuesCacheType), POINTER :: createValuesCache

    ENTERS("EquationsMappingVector_RHSCoefficientSet",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector equations mapping is not associated",err,error,*999)
    IF(vectorMapping%vectorMappingFinished) CALL FlagError("Vector equations mapping has been finished.",err,error,*999)

    NULLIFY(createValuesCache)
    CALL EquationsMappingVector_CreateValuesCacheGet(vectorMapping,createValuesCache,err,error,*999)
    
    IF(createValuesCache%rhsVariableType==0) &
      & CALL FlagError("The equations mapping RHS variable type has not been set.",err,error,*999)

    vectorMapping%createValuesCache%rhsCoefficient=rhsCoefficient
       
    EXITS("EquationsMappingVector_RHSCoefficientSet")
    RETURN
999 ERRORSEXITS("EquationsMappingVector_RHSCoefficientSet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_RHSCoefficientSet
  
  !
  !================================================================================================================================
  !

  !>Finalises the equations mapping RHS mapping and deallocates all memory
  SUBROUTINE EquationsMappingVector_RHSMappingFinalise(rhsMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping !<A pointer to the RHS mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingVector_RHSMappingFinalise",err,error,*999)

    IF(ASSOCIATED(rhsMapping)) THEN
      IF(ASSOCIATED(rhsMapping%rhsDOFToEquationsRowMap)) DEALLOCATE(rhsMapping%rhsDOFToEquationsRowMap)
      IF(ASSOCIATED(rhsMapping%equationsRowToRHSDOFMap)) DEALLOCATE(rhsMapping%equationsRowToRHSDOFMap)
      DEALLOCATE(rhsMapping)
    ENDIF
       
    EXITS("EquationsMappingVector_RHSMappingFinalise")
    RETURN
999 ERRORSEXITS("EquationsMappingVector_RHSMappingFinalise",err,error)
    RETURN 1
  END SUBROUTINE EquationsMappingVector_RHSMappingFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the equations mapping RHS mapping
  SUBROUTINE EquationsMappingVector_RHSMappingInitialise(vectorMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations mapping to initialise the RHS mapping for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("EquationsMappingVector_RHSMappingInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector equations mapping is not associated.",err,error,*998)
    IF(ASSOCIATED(vectorMapping%rhsMapping)) &
      & CALL FlagError("Vector equations mapping RHS mapping is already associated.",err,error,*998)
    
    ALLOCATE(vectorMapping%rhsMapping,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate equations mapping RHS mapping.",err,error,*999)
    vectorMapping%rhsMapping%vectorMapping=>vectorMapping        
    vectorMapping%rhsMapping%rhsVariableType=0
    NULLIFY(vectorMapping%rhsMapping%rhsVariable)
    NULLIFY(vectorMapping%rhsMapping%rhsVariableMapping)
    vectorMapping%rhsMapping%rhsCoefficient=1.0_DP
       
    EXITS("EquationsMappingVector_RHSMappingInitialise")
    RETURN
999 CALL EquationsMappingVector_RHSMappingFinalise(vectorMapping%rhsMapping,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMappingVector_RHSMappingInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_RHSMappingInitialise

  !
  !================================================================================================================================
  !

  !>Sets the mapping between a dependent field variable and the equations set rhs vector.
  SUBROUTINE EquationsMappingVector_RHSVariableTypeSet(vectorMapping,rhsVariableType,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations mapping to set
    INTEGER(INTG), INTENT(IN) :: rhsVariableType !<The variable type associated with the equations set rhs vector. If the problem does not have a rhs vector then the variable type on input should be zero.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx,residualIdx
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(FieldType), POINTER :: dependentField
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsMappingVector_RHSVariableTypeSet",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector equations mapping is not associated",err,error,*999)
    IF(vectorMapping%vectorMappingFinished) CALL FlagError("Vector equations mapping has been finished.",err,error,*999)

    NULLIFY(createValuesCache)
    CALL EquationsMappingVector_CreateValuesCacheGet(vectorMapping,createValuesCache,err,error,*999)
    NULLIFY(vectorEquations)
    CALL EquationsMappingVector_VectorEquationsGet(vectorMapping,vectorEquations,err,error,*999)
    NULLIFY(equations)
    CALL EquationsVector_EquationsGet(vectorEquations,equations,err,error,*999)
    NULLIFY(equationsSet)
    CALL Equations_EquationsSetGet(equations,equationsSet,err,error,*999)

    IF(rhsVariableType==0) THEN
      createValuesCache%rhsVariableType=0
    ELSE
      NULLIFY(dependentField)
      CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
      
      !Check the RHS variable type is not being by other equations matrices or vectors
      IF(createValuesCache%dynamicVariableType==rhsVariableType) THEN
        localError="The specified RHS variable type of "//TRIM(NumberToVString(rhsVariableType,"*",err,error))// &
          & " is the same as the variable type for the dynamic matrices."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      DO residualIdx=1,createValuesCache%numberOfResiduals
        DO matrixIdx=1,createValuesCache%numberOfResidualVariables(residualIdx)
          IF(createValuesCache%residualVariableTypes(matrixIdx,residualIdx)==rhsVariableType) THEN
            localError="The specified RHS variable type of "//TRIM(NumberToVString(rhsVariableType,"*",err,error))// &
              & " is the same as the variable type for variable index "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
              & " of residual index "//TRIM(NumberToVString(residualIdx,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ENDDO !matrixIdx
      ENDDO !residualIdx
      DO matrixIdx=1,createValuesCache%numberOfLinearMatrices
        IF(createValuesCache%linearMatrixVariableTypes(matrixIdx)==rhsVariableType) THEN
          localError="The specified RHS variable type of "//TRIM(NumberToVString(rhsVariableType,"*",err,error))// &
            & " is the same as the variable type for linear matrix number "//TRIM(NumberToVString(matrixIdx,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDDO !matrixIdx
      !Check the RHS variable number is defined on the dependent field
      IF(rhsVariableType<1.OR.rhsVariableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
        localError="The specified RHS variable type of "//TRIM(NumberToVString(rhsVariableType,"*",err,error))// &
          & " is invalid. The number must either be zero or >= 1 and <= "// &
          & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      IF(.NOT.ASSOCIATED(dependentField%variableTypeMap(rhsVariableType)%ptr)) THEN
        localError="The specified RHS variable type of "//TRIM(NumberToVString(rhsVariableType,"*",err,error))// &
          & " is not defined on the dependent field."
        CALL FlagError(localError,err,error,*999)
      ENDIF

      createValuesCache%rhsVariableType=rhsVariableType
      
    ENDIF
       
    EXITS("EquationsMappingVector_RHSVariableTypeSet")
    RETURN
999 ERRORSEXITS("EquationsMappingVector_RHSVariableTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_RHSVariableTypeSet
  
  !
  !================================================================================================================================
  !

  !>Sets the coefficient applied to the equations set source vector.
  SUBROUTINE EquationsMappingVector_SourcesCoefficientsSet0(vectorMapping,sourceCoefficient,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the vector equations mapping to set the source coefficient for
    REAL(DP), INTENT(IN) :: sourceCoefficient !<The coefficient applied to the equations set source vector.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsMappingVector_SourcesCoefficientsSet0",err,error,*999)

    CALL EquationsMappingVector_SourcesCoefficientsSet1(vectorMapping,[sourceCoefficient],err,error,*999)
      
    EXITS("EquationsMappingVector_SourcesCoefficientsSet0")
    RETURN
999 ERRORSEXITS("EquationsMappingVector_SourcesCoefficientsSet0",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_SourcesCoefficientsSet0
  
  !
  !================================================================================================================================
  !

  !>Sets the coefficients applied to the equations set source vectors.
  SUBROUTINE EquationsMappingVector_SourcesCoefficientsSet1(vectorMapping,sourceCoefficients,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the vector equations mapping to set the source coefficients for
    REAL(DP), INTENT(IN) :: sourceCoefficients(:) !<sourceCoefficients(sourceIdx). The coefficients applied to the equations set sourceIdx'th source vector.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(EquationsMappingVectorCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsMappingVector_SourcesCoefficientSet1",err,error,*999)

    CALL EquationsMappingVector_AssertNotFinished(vectorMapping,err,error,*999)
    NULLIFY(createValuesCache)
    CALL EquationsMappingVector_CreateValuesCacheGet(vectorMapping,createValuesCache,err,error,*999)
    IF(SIZE(sourceCoefficients,1)<createValuesCache%numberOfSources) THEN
      localError="The size of the specified source coefficients array of "// &
        & TRIM(NumberToVString(SIZE(sourceCoefficients,1),"*",err,error))// &
        & " is invalid. The size of source coefficients array should be >= "// &
        & TRIM(NumberToVString(createValuesCache%numberOfSources,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    createValuesCache%sourceCoefficients(1:createValuesCache%numberOfSources)= &
      & sourceCoefficients(1:createValuesCache%numberOfSources)
      
    EXITS("EquationsMappingVector_SourcesCoefficientsSet1")
    RETURN
999 ERRORSEXITS("EquationsMappingVector_SourcesCoefficientsSet1",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_SourcesCoefficientsSet1
  
  !
  !================================================================================================================================
  !

  !>Finalises the equations mapping sources mapping and deallocates all memory
  SUBROUTINE EquationsMappingVector_SourcesMappingFinalise(sourcesMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingSourcesType), POINTER :: sourcesMapping !<A pointer to the sources mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: sourceIdx
 
    ENTERS("EquationsMappingVector_SourcesMappingFinalise",err,error,*999)

    IF(ASSOCIATED(sourcesMapping)) THEN
      IF(ALLOCATED(sourcesMapping%sources)) THEN
        DO sourceIdx=1,SIZE(sourcesMapping%sources,1)
          CALL EquationsMappingSources_SourceMappingFinalise(sourcesMapping%sources(sourceIdx)%ptr,err,error,*999)
        ENDDO !sourceIdx
        DEALLOCATE(sourcesMapping%sources)
      ENDIF
      DEALLOCATE(sourcesMapping)
    ENDIF
       
    EXITS("EquationsMappingVector_SourcesMappingFinalise")
    RETURN
999 ERRORSEXITS("EquationsMappingVector_SourcesMappingFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_SourcesMappingFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the equations mapping sources mapping
  SUBROUTINE EquationsMappingVector_SourcesMappingInitialise(vectorMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations mapping to initialise the sources mapping for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("EquationsMappingVector_SourcesMappingInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector equations mapping is not associated.",err,error,*998)
    IF(ASSOCIATED(vectorMapping%sourcesMapping)) &
      & CALL FlagError("Equations mapping sources mapping is already associated.",err,error,*998)

    ALLOCATE(vectorMapping%sourcesMapping,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate equations mapping sources mapping.",err,error,*999)
    vectorMapping%sourcesMapping%vectorMapping=>vectorMapping        
    vectorMapping%sourcesMapping%numberOfSources=0
       
    EXITS("EquationsMappingVector_SourcesMappingInitialise")
    RETURN
999 CALL EquationsMappingVector_SourcesMappingFinalise(vectorMapping%sourcesMapping,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMappingVector_SourcesMappingInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_SourcesMappingInitialise

  !
  !================================================================================================================================
  !

  !>Sets the mapping between a source field variable and the equations set source vector.
  SUBROUTINE EquationsMappingVector_SourcesVariableTypesSet0(vectorMapping,sourceVariableType,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the vector equations mapping to set
    INTEGER(INTG), INTENT(IN) :: sourceVariableType !<The variable type associated with the equations set source vector.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsMappingVector_SourcesVariableTypesSet0",err,error,*999)

    CALL EquationsMappingVector_SourcesVariableTypesSet1(vectorMapping,[sourceVariableType],err,error,*999)
    
    EXITS("EquationsMappingVector_SourcesVariableTypesSet0")
    RETURN
999 ERRORSEXITS("EquationsMappingVector_SourcesVariableTypesSet0",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_SourcesVariableTypesSet0
  
  !
  !================================================================================================================================
  !

  !>Sets the mapping between a source field variables and the equations set source vectors.
  SUBROUTINE EquationsMappingVector_SourcesVariableTypesSet1(vectorMapping,sourcesVariableTypes,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the vector equations mapping to set
    INTEGER(INTG), INTENT(IN) :: sourcesVariableTypes(:) !<sourcesVariableTypes(sourceIdx). The variable type associated with the equations set sourceIdx'th source vector. If the problem does not have a source vector then the variable type on input should be zero.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: sourceIdx,sourceVariableType
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(FieldType), POINTER :: sourceField
    TYPE(FieldVariableType), POINTER :: sourceVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsMappingVector_SourcesVariableTypesSet1",err,error,*999)

    CALL EquationsMappingVector_AssertNotFinished(vectorMapping,err,error,*999)
    NULLIFY(createValuesCache)
    CALL EquationsMappingVector_CreateValuesCacheGet(vectorMapping,createValuesCache,err,error,*999)
    IF(SIZE(sourcesVariableTypes,1)<createValuesCache%numberOfSources) THEN
      localError="The size of the specified sources coefficients array of "// &
        & TRIM(NumberToVString(SIZE(sourcesVariableTypes,1),"*",err,error))// &
        & " is invalid. The size of sources coefficients array should be >= "// &
        & TRIM(NumberToVString(createValuesCache%numberOfSources,"*",err,error))//")."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    NULLIFY(vectorEquations)
    CALL EquationsMappingVector_VectorEquationsGet(vectorMapping,vectorEquations,err,error,*999)
    NULLIFY(equations)
    CALL EquationsVector_EquationsGet(vectorEquations,equations,err,error,*999)
    NULLIFY(equationsSet)
    CALL Equations_EquationsSetGet(equations,equationsSet,err,error,*999)
    NULLIFY(sourceField)
    CALL EquationsSet_SourceFieldGet(equationsSet,sourceField,err,error,*999)
    !Check the source variables
    DO sourceIdx=1,createValuesCache%numberOfSources
      sourceVariableType=sourcesVariableTypes(sourceIdx)
      IF(sourceVariableType<1.OR.sourceVariableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
        localError="The specified source variable type of "//TRIM(NumberToVString(sourceVariableType,"*",err,error))// &
          & " for source index "//TRIM(NumberToVString(sourceIdx,"*",err,error))// &
          & " is invalid. The source variable type must be >= 1 and <= "// &
          & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      NULLIFY(sourceVariable)
      CALL Field_VariableExists(sourceField,sourceVariableType,sourceVariable,err,error,*999)
      IF(.NOT.ASSOCIATED(sourceVariable)) THEN
        localError="The specified source variable type of "//TRIM(NumberToVString(sourceVariableType,"*",err,error))// &
          & " for source index "//TRIM(NumberToVString(sourceIdx,"*",err,error))//" is not defined on the source field."
        CALL FlagError(localError,err,error,*999)
      ENDIF      
    ENDDO !sourceIdx

    createValuesCache%sourceVariableTypes(1:createValuesCache%numberOfSources)= &
      & sourcesVariableTypes(1:createValuesCache%numberOfSources)
    
    EXITS("EquationsMappingVector_SourcesVariableTypesSet1")
    RETURN
999 ERRORSEXITS("EquationsMappingVector_SourcesVariableTypesSet1",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_SourcesVariableTypesSet1
  
  !
  !================================================================================================================================
  !

  !>Finalise an equations mapping equations matrix map.
  SUBROUTINE EquationsMappingVToECM_Finalise(varToEquationsColumnMap,err,error,*)

    !Argument variables
    TYPE(VarToEquationsColumnMapType) :: varToEquationsColumnMap !<The variable dof to equations column map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsMappingVToECM_Finalise",err,error,*999)
    
    IF(ALLOCATED(varToEquationsColumnMap%columnDOF)) DEALLOCATE(varToEquationsColumnMap%columnDOF)
    
    EXITS("EquationsMappingVToECM_Finalise")
    RETURN
999 ERRORS("EquationsMappingVToECM_Finalise",err,error)    
    EXITS("EquationsMappingVToECM_Finalise")    
    RETURN 1
   
  END SUBROUTINE EquationsMappingVToECM_Finalise

  !
  !================================================================================================================================
  !

  !>Finalises a variable to equations matrices map and deallocates all memory.
  SUBROUTINE EquationsMappingVToEMSMap_Finalise(varToEquationsMatricesMap,err,error,*)

    !Argument variables
    TYPE(VarToEquationsMatricesMapType), POINTER :: varToEquationsMatricesMap !<A pointer to the variable to equations matrices map to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx

    ENTERS("EquationsMappingVToEMSMap_Finalise",err,error,*999)

    IF(ASSOCIATED(varToEquationsMatricesMap)) THEN
      IF(ALLOCATED(varToEquationsMatricesMap%equationsMatrixNumbers)) DEALLOCATE(varToEquationsMatricesMap%equationsMatrixNumbers)
      IF(ALLOCATED(varToEquationsMatricesMap%dofToColumnsMaps)) THEN
        DO matrixIdx=1,SIZE(varToEquationsMatricesMap%dofToColumnsMaps,1)
          CALL EquationsMappingVToECM_Finalise(varToEquationsMatricesMap%dofToColumnsMaps(matrixIdx),err,error,*999)
        ENDDO !matrixIdx
        DEALLOCATE(varToEquationsMatricesMap%dofToColumnsMaps)
      ENDIF
    ENDIF
    
    EXITS("EquationsMappingVToEMSMap_Finalise")
    RETURN
999 ERRORS("EquationsMappingVToEMSMap_Finalise",err,error)    
    EXITS("EquationsMappingVToEMSMap_Finalise")    
    RETURN 1
   
  END SUBROUTINE EquationsMappingVToEMSMap_Finalise

  !
  !================================================================================================================================
  !

  !>Allocates and initialises an equations mapping equations matrix map.
  SUBROUTINE EquationsMappingVToEMSMap_Initialise(varToEquationsMatricesMap,err,error,*)

    !Argument variables
    TYPE(VarToEquationsMatricesMapType), POINTER :: varToEquationsMatricesMap !<A pointer to the variable to equations matrices map to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("EquationsMappingVToEMSMap_Initialise",err,error,*998)

    IF(ASSOCIATED(varToEquationsMatricesMap)) &
      & CALL FlagError("Variable to equations matrices map is already associated.",err,error,*998)

    ALLOCATE(varToEquationsMatricesMap,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate variable to equations matrices map.",err,error,*999)
    varToEquationsMatricesMap%variableIndex=0
    varToEquationsMatricesMap%variableType=0
    NULLIFY(varToEquationsMatricesMap%variable)
    varToEquationsMatricesMap%numberOfEquationsMatrices=0
    
    EXITS("EquationsMappingVToEMSMap_Initialise")
    RETURN
999 CALL EquationsMappingVToEMSMap_Finalise(varToEquationsMatricesMap,dummyErr,dummyError,*998)
998 ERRORS("EquationsMappingVToEMSMap_Initialise",err,error)    
    EXITS("EquationsMappingVToEMSMap_Initialise")    
    RETURN 1
   
  END SUBROUTINE EquationsMappingVToEMSMap_Initialise

  !
  !================================================================================================================================
  !

  !>Finalises a variable to equations Jacobian map and deallocates all memory.
  SUBROUTINE EquationsMappingVToJMMap_Finalise(varToEquationsJacobianMap,err,error,*)

    !Argument variables
    TYPE(VarToJacobianMatrixMapType), POINTER :: varToEquationsJacobianMap !<A pointer to the variable to equations Jacobian map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsMappingVToJMMap_Finalise",err,error,*999)

    IF(ASSOCIATED(varToEquationsJacobianMap)) THEN
      IF(ALLOCATED(varToEquationsJacobianMap%dofToColumnsMap)) DEALLOCATE(varToEquationsJacobianMap%dofToColumnsMap)
      DEALLOCATE(varToEquationsJacobianMap)
    ENDIF
    
    EXITS("EquationsMappingVToJMMap_Finalise")
    RETURN
999 ERRORS("EquationsMappingVToJMMap_Finalise",err,error)    
    EXITS("EquationsMappingVToJMMap_Finalise")    
    RETURN 1
   
  END SUBROUTINE EquationsMappingVToJMMap_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises a variable to equations Jacobian map
  SUBROUTINE EquationsMappingVToJMMap_Initialise(varToEquationsJacobianMap,err,error,*)

    !Argument variables
    TYPE(VarToJacobianMatrixMapType), POINTER :: varToEquationsJacobianMap !<A pointer to the variable to equations Jacobian map to initialise. Must not be allocated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("EquationsMappingVToJMMap_Initialise",err,error,*998)

    IF(ASSOCIATED(varToEquationsJacobianMap)) &
      &  CALL FlagError("Variable to equations Jacobian map is already associated.",err,error,*998)

    ALLOCATE(varToEquationsJacobianMap,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate the variable to equations Jacobian map.",err,error,*999)
    varToEquationsJacobianMap%variableType=0
    NULLIFY(varToEquationsJacobianMap%variable)
    
    EXITS("EquationsMappingVToJMMap_Initialise")
    RETURN
999 CALL EquationsMappingVToJMMap_Finalise(varToEquationsJacobianMap,dummyErr,dummyError,*998)
998 ERRORS("EquationsMappingVToJMMap_Initialise",err,error)    
    EXITS("EquationsMappingVToJMMap_Initialise")    
    RETURN 1
   
  END SUBROUTINE EquationsMappingVToJMMap_Initialise

  !
  !================================================================================================================================
  !
  
END MODULE EquationsMappingRoutines
