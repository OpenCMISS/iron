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
!> Contributor(s):
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
  USE DOMAIN_MAPPINGS
  USE EquationsAccessRoutines
  USE EquationsMappingAccessRoutines
  USE EquationsSetConstants
  USE EquationsSetAccessRoutines
  USE FIELD_ROUTINES
  USE FieldAccessRoutines
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE Kinds
  USE Strings
  USE Types

#include "macros.h"  
 
  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  INTERFACE EquationsMapping_DynamicMatricesSet
    MODULE PROCEDURE EquationsMapping_DynamicMatricesSet1
    MODULE PROCEDURE EquationsMapping_DynamicMatricesSet2
  END INTERFACE EquationsMapping_DynamicMatricesSet
  
  INTERFACE EquationsMapping_DynamicMatricesCoeffsSet
    MODULE PROCEDURE EquationsMapping_DynamicMatricesCoeffsSet1
    MODULE PROCEDURE EquationsMapping_DynamicMatricesCoeffsSet2
  END INTERFACE EquationsMapping_DynamicMatricesCoeffsSet
  
  PUBLIC EquationsMapping_DynamicMatricesSet

  PUBLIC EquationsMapping_DynamicMatricesCoeffsSet

  PUBLIC EquationsMapping_DynamicVariableTypeSet

  PUBLIC EquationsMapping_LinearMatricesCoeffsSet

  PUBLIC EquationsMapping_LinearMatricesNumberSet

  PUBLIC EquationsMapping_LinearMatricesVariableTypesSet
  
  PUBLIC EquationsMapping_ResidualCoeffSet

  PUBLIC EquationsMapping_ResidualVariableTypesSet

  PUBLIC EquationsMapping_ResidualVariablesNumberSet
  
  PUBLIC EquationsMapping_RHSCoeffSet

  PUBLIC EquationsMapping_RHSVariableTypeSet

  PUBLIC EquationsMapping_ScalarDestroy

  PUBLIC EquationsMapping_SourceCoeffSet

  PUBLIC EquationsMapping_SourceVariableTypeSet

  PUBLIC EquationsMapping_VectorCreateFinish,EquationsMapping_VectorCreateStart

  PUBLIC EquationsMapping_VectorDestroy

CONTAINS

  !
  !================================================================================================================================
  !

  !>Calculates the equations vector mapping.
  SUBROUTINE EquationsMapping_VectorCalculate(vectorMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations vector mapping to calculate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: columnIdx,dofIdx,linearMatrixStart,matrixIdx,numberOfRows,numberOfGlobalRows,rowIdx, &
      & totalNumberOfRows,variableIdx,variableType,variableTypeIdx
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsMappingVectorCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping
    TYPE(EquationsMappingSourceType), POINTER :: sourceMapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    TYPE(FIELD_TYPE), POINTER :: dependentField,sourceField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: dependentVariable,sourceVariable,rowVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsMapping_VectorCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Equations mapping is not associated.",err,error,*999)
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
    NULLIFY(sourceField)
    IF(createValuesCache%sourceVariableType/=0) CALL EquationsSet_SourceFieldGet(equationsSet,sourceField,err,error,*999)

    !Calculate the number of rows in the equations set
    linearMatrixStart=1
    SELECT CASE(equations%timeDependence)
    CASE(EQUATIONS_STATIC,EQUATIONS_QUASISTATIC)
      SELECT CASE(equations%linearity)
      CASE(EQUATIONS_LINEAR,EQUATIONS_NONLINEAR_BCS)
        !Static linear equations set
        IF(createValuesCache%numberOfLinearMatrices>=1) THEN
          linearMatrixStart=2
          dependentVariable=>dependentField%VARIABLE_TYPE_MAP(createValuesCache%linearMatrixVariableTypes(1))%ptr
        ELSE
          localError="The create values cache number of linear matrices of "// &
            & TRIM(NumberToVString(createValuesCache%numberOfLinearMatrices,"*",err,error))// &
            & " is invalid. The number of linear equations matrices must be >= 1 for a linear equations set."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      CASE(EQUATIONS_NONLINEAR)
        !Static nonlinear equations set
        !Use first listed nonlinear variable
        IF(createValuesCache%numberOfResidualVariables>=1) THEN
          dependentVariable=>dependentField%VARIABLE_TYPE_MAP(createValuesCache%residualVariableTypes(1))%ptr
        ELSE
          localError="The create values cache number of residual vectors of "// &
            & TRIM(NumberToVString(createValuesCache%numberOfResidualVariables,"*",err,error))// &
            & " is invalid. The number of Jacobian matrices must be >= 1 for a nonlinear equations set."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      CASE DEFAULT
        localError="The equations linearity type of "//TRIM(NumberToVString(equations%linearity,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_FIRST_ORDER_DYNAMIC,EQUATIONS_SECOND_ORDER_DYNAMIC)
      SELECT CASE(equations%linearity)
      CASE(EQUATIONS_LINEAR,EQUATIONS_NONLINEAR_BCS)
        !Dynamic linear equations set
        dependentVariable=>dependentField%VARIABLE_TYPE_MAP(createValuesCache%dynamicVariableType)%ptr
      CASE(EQUATIONS_NONLINEAR)
        !Dynamic nonlinear equations set
        dependentVariable=>dependentField%VARIABLE_TYPE_MAP(createValuesCache%residualVariableTypes(1))%ptr
      CASE DEFAULT
        localError="The equations linearity type of "//TRIM(NumberToVString(equations%linearity,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_TIME_STEPPING)
      !Time stepping DAE equations set
!!NOTE: The time stepping variable type doesn't have to come from the dependent field, it could come from, say, the source field.
      !dependentVariable=>dependentField%VARIABLE_TYPE_MAP(createValuesCache%TIME_STEPPING_VARIABLE_TYPE)%ptr
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The equations time dependence type of "//TRIM(NumberToVString(equations%timeDependence,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    IF(.NOT.ASSOCIATED(dependentVariable)) CALL FlagError("The dependent variable is not associated.",err,error,*999)
    numberOfRows=dependentVariable%NUMBER_OF_DOFS
    totalNumberOfRows=dependentVariable%TOTAL_NUMBER_OF_DOFS
    vectorMapping%rowDofsMapping=>dependentVariable%DOMAIN_MAPPING
    IF(.NOT.ASSOCIATED(vectorMapping%rowDofsMapping)) &
      & CALL FlagError("Dependent variable domain mapping is not associated.",err,error,*999)
    numberOfGlobalRows=vectorMapping%rowDofsMapping%NUMBER_OF_GLOBAL
    
    !Check that the number of rows is consistent across the remaining linear matrices
    DO matrixIdx=linearMatrixStart,createValuesCache%numberOfLinearMatrices
      dependentVariable=>dependentField%VARIABLE_TYPE_MAP(createValuesCache%linearMatrixVariableTypes(matrixIdx))%ptr
      IF(ASSOCIATED(dependentVariable)) THEN
        IF(dependentVariable%NUMBER_OF_DOFS/=numberOfRows) THEN
          localError="Invalid equations set up. The number of rows in the equations set ("// &
            & TRIM(NumberToVString(numberOfRows,"*",err,error))// &
            & ") does not match the number of rows in equations linear matrix number "// &
            & TRIM(NumberToVString(matrixIdx,"*",err,error))//" ("// &
                      & TRIM(NumberToVString(dependentVariable%NUMBER_OF_DOFS,"*",err,error))//")."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        IF(dependentVariable%TOTAL_NUMBER_OF_DOFS/=totalNumberOfRows) THEN
          localError="Invalid equations set up. The total number of rows in the equations set ("// &
            & TRIM(NumberToVString(totalNumberOfRows,"*",err,error))// &
            & ") does not match the total number of rows in equations matrix number "// &
            & TRIM(NumberToVString(matrixIdx,"*",err,error))//" ("// &
            & TRIM(NumberToVString(dependentVariable%TOTAL_NUMBER_OF_DOFS,"*",err,error))//")."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        localError="The dependent variable mapped to linear matrix number "// &
          & TRIM(NumberToVString(matrixIdx,"*",err,error))//" is not associated."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDDO !matrixIdx
    !Check the Jacobian matrices
    !Can't check the number of rows now as Jacobian's might not be square so just check variables are associated
    DO matrixIdx=1,createValuesCache%numberOfResidualVariables
      dependentVariable=>dependentField%VARIABLE_TYPE_MAP(createValuesCache%residualVariableTypes(matrixIdx))%ptr
      IF(.NOT.ASSOCIATED(dependentVariable)) THEN
        localError="The dependent variable mapped to Jacobian matrix number "// &
          & TRIM(NumberToVString(matrixIdx,"*",err,error))//" is not associated."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDDO !matrixIdx
    !Check that the number of rows are consistent with the RHS vector if it exists
    IF(createValuesCache%rhsVariableType/=0) THEN
      dependentVariable=>dependentField%VARIABLE_TYPE_MAP(createValuesCache%rhsVariableType)%ptr
      IF(ASSOCIATED(dependentVariable)) THEN
        IF(dependentVariable%NUMBER_OF_DOFS/=numberOfRows) THEN
          localError="Invalid equations set up. The number of rows in the equations set ("// &
            & TRIM(NumberToVString(numberOfRows,"*",err,error))// &
            & ") does not match the number of rows in the RHS vector ("// &
            & TRIM(NumberToVString(dependentVariable%NUMBER_OF_DOFS,"*",err,error))//")."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        IF(dependentVariable%TOTAL_NUMBER_OF_DOFS/=totalNumberOfRows) THEN
          localError="Invalid equations set up. The total number of rows in the equations set ("// &
            & TRIM(NumberToVString(totalNumberOfRows,"*",err,error))// &
            & ") does not match the total number of rows in the RHS vector ("// &
            & TRIM(NumberToVString(dependentVariable%TOTAL_NUMBER_OF_DOFS,"*",err,error))//")."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("The dependent variable mapped to the RHS vector is not associated.",err,error,*999)
      ENDIF
    ENDIF
    !Check that the number of rows are consistent with the source vector if it exists
    !IF(createValuesCache%sourceVariableType/=0) THEN
    !  sourceVariable=>sourceField%VARIABLE_TYPE_MAP(createValuesCache%sourceVariableType)%ptr
    !  IF(ASSOCIATED(sourceVariable)) THEN
    !    IF(sourceVariable%NUMBER_OF_DOFS/=numberOfRows) THEN
    !      localError="Invalid equations set up. The number of rows in the equations set ("// &
    !        & TRIM(NumberToVString(numberOfRows,"*",err,error))// &
    !        & ") does not match the number of rows in the source vector ("// &
    !        & TRIM(NumberToVString(sourceVariable%NUMBER_OF_DOFS,"*",err,error))//")."
    !      CALL FlagError(localError,err,error,*999)
    !    ENDIF
    !    IF(sourceVariable%TOTAL_NUMBER_OF_DOFS/=totalNumberOfRows) THEN
    !      localError="Invalid equations set up. The total number of rows in the equations set ("// &
    !        & TRIM(NumberToVString(totalNumberOfRows,"*",err,error))// &
    !        & ") does not match the total number of rows in the source vector ("// &
    !        & TRIM(NumberToVString(sourceVariable%TOTAL_NUMBER_OF_DOFS,"*",err,error))//")."
    !      CALL FlagError(localError,err,error,*999)
    !    ENDIF
    !  ELSE
    !    CALL FlagError("The source variable mapped to the source vector is not associated.",err,error,*999)
    !  ENDIF
    !ENDIF
    vectorMapping%numberOfRows=numberOfRows
    vectorMapping%totalNumberOfRows=totalNumberOfRows
    vectorMapping%numberOfGlobalRows=numberOfGlobalRows
    !Calculate dynamic mappings
    IF(createValuesCache%dynamicVariableType/=0) THEN
      CALL EquationsMapping_DynamicMappingInitialise(vectorMapping,err,error,*999)
      NULLIFY(dynamicMapping)
      CALL EquationsMappingVector_DynamicMappingGet(vectorMapping,dynamicMapping,err,error,*999)
      dynamicMapping%numberOfDynamicMatrices=createValuesCache%numberOfDynamicMatrices
      dynamicMapping%stiffnessMatrixNumber=createValuesCache%dynamicStiffnessMatrixNumber
      dynamicMapping%dampingMatrixNumber=createValuesCache%dynamicDampingMatrixNumber
      dynamicMapping%massMatrixNumber=createValuesCache%dynamicMassMatrixNumber
      !Initialise the variable type maps
      ALLOCATE(dynamicMapping%varToEquationsMatricesMaps(FIELD_NUMBER_OF_VARIABLE_TYPES),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate equations mapping variable to equations map.",err,error,*999)
      DO variableTypeIdx=1,FIELD_NUMBER_OF_VARIABLE_TYPES
        CALL EquationsMapping_VarToEquatsMatricesMapInitialise(dynamicMapping%varToEquationsMatricesMaps(variableTypeIdx), &
          & err,error,*999)
        dynamicMapping%varToEquationsMatricesMaps(variableTypeIdx)%variableIndex=variableTypeIdx
        dynamicMapping%varToEquationsMatricesMaps(variableTypeIdx)%variableType=variableTypeIdx
        dynamicMapping%varToEquationsMatricesMaps(variableTypeIdx)%variable=>dependentField%VARIABLE_TYPE_MAP(variableTypeIdx)%ptr
      ENDDO !variableTypeIdx
      dynamicMapping%varToEquationsMatricesMaps(createValuesCache%dynamicVariableType)% &
        & numberOfEquationsMatrices=createValuesCache%numberOfDynamicMatrices
      IF(createValuesCache%rhsVariableType/=0) &
        & dynamicMapping%varToEquationsMatricesMaps(createValuesCache%rhsVariableType)%numberOfEquationsMatrices=-1
      !Allocate and initialise the variable to equations matrices maps
      DO variableTypeIdx=1,FIELD_NUMBER_OF_VARIABLE_TYPES
        dependentVariable=>dependentField%VARIABLE_TYPE_MAP(variableTypeIdx)%ptr
        IF(ASSOCIATED(dependentVariable)) THEN
          IF(dynamicMapping%varToEquationsMatricesMaps(variableTypeIdx)%numberOfEquationsMatrices==-1) THEN
!!TODO: check if this can be removed and just allocate those variables that are actually used
            ALLOCATE(dynamicMapping%varToEquationsMatricesMaps(variableTypeIdx)%dofToRowsMap( &
              & dependentVariable%TOTAL_NUMBER_OF_DOFS),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate variable to equations matrices maps dof to rows map.",err,error,*999)
            dynamicMapping%varToEquationsMatricesMaps(variableTypeIdx)%dofToRowsMap=0
          ELSE IF(dynamicMapping%varToEquationsMatricesMaps(variableTypeIdx)%numberOfEquationsMatrices>0) THEN
            ALLOCATE(dynamicMapping%varToEquationsMatricesMaps(variableTypeIdx)%equationsMatrixNumbers( &
              & dynamicMapping%varToEquationsMatricesMaps(variableTypeIdx)%numberOfEquationsMatrices),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate variable to equations matrices maps equations matrix numbers.", &
              & err,error,*999)
            ALLOCATE(dynamicMapping%varToEquationsMatricesMaps(variableTypeIdx)%dofToColumnsMaps( &
              & dynamicMapping%varToEquationsMatricesMaps(variableTypeIdx)%numberOfEquationsMatrices),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate variable to equations matrices maps dof to columns map.",err,error,*999)
            dynamicMapping%varToEquationsMatricesMaps(variableTypeIdx)%equationsMatrixNumbers=0
            IF(createValuesCache%dynamicVariableType==variableTypeIdx) THEN
              IF(createValuesCache%dynamicStiffnessMatrixNumber/=0) &
                & dynamicMapping%varToEquationsMatricesMaps(variableTypeIdx)%equationsMatrixNumbers( &
                & createValuesCache%dynamicStiffnessMatrixNumber)=createValuesCache%dynamicStiffnessMatrixNumber
              IF(createValuesCache%dynamicDampingMatrixNumber/=0) &
                & dynamicMapping%varToEquationsMatricesMaps(variableTypeIdx)%equationsMatrixNumbers( &
                & createValuesCache%dynamicDampingMatrixNumber)=createValuesCache%dynamicDampingMatrixNumber
              IF(createValuesCache%dynamicMassMatrixNumber/=0) &
                & dynamicMapping%varToEquationsMatricesMaps(variableTypeIdx)%equationsMatrixNumbers( &
                & createValuesCache%dynamicMassMatrixNumber)=createValuesCache%dynamicMassMatrixNumber
              DO matrixIdx=1,dynamicMapping%varToEquationsMatricesMaps(variableTypeIdx)%numberOfEquationsMatrices
                ALLOCATE(dynamicMapping%varToEquationsMatricesMaps(variableTypeIdx)%dofToColumnsMaps( &
                  & matrixIdx)%columnDOF(dependentVariable%TOTAL_NUMBER_OF_DOFS),STAT=err)
                IF(err/=0) CALL FlagError("Could not allocate variable dof to columns map column dof.",err,error,*999)
                DO dofIdx=1,dependentVariable%TOTAL_NUMBER_OF_DOFS
                  !1-1 mapping for now
                  columnIdx=dependentVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(dofIdx)
                  dynamicMapping%varToEquationsMatricesMaps(variableTypeIdx)%dofToColumnsMaps(matrixIdx)%columnDOF(dofIdx)=columnIdx
                ENDDO !dofIdx                          
              ENDDO !matrixIdx
              ALLOCATE(dynamicMapping%varToEquationsMatricesMaps(variableTypeIdx)%dofToRowsMap( &
                & dependentVariable%TOTAL_NUMBER_OF_DOFS),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate variable to equations matrices maps dof to rows map.",err,error,*999)
              DO dofIdx=1,dependentVariable%TOTAL_NUMBER_OF_DOFS
                !1-1 mappings for now.
                rowIdx=dofIdx
                dynamicMapping%varToEquationsMatricesMaps(variableTypeIdx)%dofToRowsMap(dofIdx)=rowIdx
              ENDDO !dofIdx
            ENDIF
          ENDIF
        ENDIF
      ENDDO !variableTypeIdx
      !Allocate and initialise the equations matrix to variable maps types
      ALLOCATE(dynamicMapping%equationsMatrixToVarMaps(dynamicMapping%numberOfDynamicMatrices),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate equations mapping equations matrix to variable maps.",err,error,*999)
      !Create the individual matrix maps and column maps
      variableType=createValuesCache%dynamicVariableType
      dependentVariable=>dependentField%VARIABLE_TYPE_MAP(variableType)%ptr
      dynamicMapping%dynamicVariableType=variableType
      dynamicMapping%dynamicVariable=>dependentVariable
      DO matrixIdx=1,dynamicMapping%numberOfDynamicMatrices
        CALL EquationsMapping_EquatsMatrixToVarMapInitialise(dynamicMapping%equationsMatrixToVarMaps(matrixIdx),err,error,*999)
        dynamicMapping%equationsMatrixToVarMaps(matrixIdx)%matrixNumber=matrixIdx
        dynamicMapping%equationsMatrixToVarMaps(matrixIdx)%variableType=variableType
        dynamicMapping%equationsMatrixToVarMaps(matrixIdx)%variable=>dependentVariable
        dynamicMapping%equationsMatrixToVarMaps(matrixIdx)%numberOfColumns=dependentVariable%DOMAIN_MAPPING%NUMBER_OF_GLOBAL
        dynamicMapping%equationsMatrixToVarMaps(matrixIdx)%matrixCoefficient=createValuesCache%dynamicMatrixCoefficients(matrixIdx)
        ALLOCATE(dynamicMapping%equationsMatrixToVarMaps(matrixIdx)%columnToDOFMap(dependentVariable%DOMAIN_MAPPING% &
          & NUMBER_OF_GLOBAL),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate equation matrix to variable map column to dof map.",err,error,*999)
        dynamicMapping%equationsMatrixToVarMaps(matrixIdx)%columnToDOFMap=0
        DO dofIdx=1,dependentVariable%TOTAL_NUMBER_OF_DOFS
          !1-1 mapping for now
          columnIdx=dependentVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(dofIdx)
          dynamicMapping%equationsMatrixToVarMaps(matrixIdx)%columnToDOFMap(columnIdx)=dofIdx
        ENDDO !dofIdx
        dynamicMapping%equationsMatrixToVarMaps(matrixIdx)%columnDOFSMapping=>dependentVariable%DOMAIN_MAPPING
      ENDDO !matrixIdx
      !Allocate the row mappings
      ALLOCATE(dynamicMapping%equationsRowToVariableDOFMaps(vectorMapping%totalNumberOfRows),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate equations row to variable dof maps.",err,error,*999)
      !Set up the row mappings
      DO rowIdx=1,vectorMapping%totalNumberOfRows
        !1-1 mapping for now
        dynamicMapping%equationsRowToVariableDOFMaps(rowIdx)=rowIdx
      ENDDO !rowIdx
    ENDIF
    !Calculate linear mappings
    IF(createValuesCache%numberOfLinearMatrices>0) THEN                  
      CALL EquationsMapping_LinearMappingInitialise(vectorMapping,err,error,*999)
      NULLIFY(linearMapping)
      CALL EquationsMappingVector_LinearMappingGet(vectorMapping,linearMapping,err,error,*999)
      linearMapping%numberOfLinearMatrices=createValuesCache%numberOfLinearMatrices
      !Allocate and initialise the variable type maps
      ALLOCATE(linearMapping%varToEquationsMatricesMaps(FIELD_NUMBER_OF_VARIABLE_TYPES),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate equations mapping variable to equations map.",err,error,*999)
      DO variableTypeIdx=1,FIELD_NUMBER_OF_VARIABLE_TYPES
        CALL EquationsMapping_VarToEquatsMatricesMapInitialise(linearMapping%varToEquationsMatricesMaps(variableTypeIdx), &
          & err,error,*999)
        linearMapping%varToEquationsMatricesMaps(variableTypeIdx)%variableIndex=variableTypeIdx
        linearMapping%varToEquationsMatricesMaps(variableTypeIdx)%variableType=variableTypeIdx
        linearMapping%varToEquationsMatricesMaps(variableTypeIdx)%variable=>dependentField%VARIABLE_TYPE_MAP(variableTypeIdx)%ptr
      ENDDO !variableTypeIdx
      !Calculate the number of variable type maps and initialise
      DO matrixIdx=1,linearMapping%numberOfLinearMatrices
        variableType=createValuesCache%linearMatrixVariableTypes(matrixIdx)
        linearMapping%varToEquationsMatricesMaps(variableType)%numberOfEquationsMatrices= &
          & linearMapping%varToEquationsMatricesMaps(variableType)%numberOfEquationsMatrices+1
      ENDDO !matrixIdx
      IF(createValuesCache%rhsVariableType/=0) &
        & linearMapping%varToEquationsMatricesMaps(createValuesCache%rhsVariableType)%numberOfEquationsMatrices=-1
      linearMapping%numberOfLinearMatrixVariables=0
      !Allocate and initialise the variable to equations matrices maps
      DO variableTypeIdx=1,FIELD_NUMBER_OF_VARIABLE_TYPES        
        dependentVariable=>dependentField%VARIABLE_TYPE_MAP(variableTypeIdx)%ptr
        IF(ASSOCIATED(dependentVariable)) THEN
          IF(linearMapping%varToEquationsMatricesMaps(variableTypeIdx)%numberOfEquationsMatrices==-1) THEN
!!TODO: check if this can be removed and just allocate those variables that are actually used
            ALLOCATE(linearMapping%varToEquationsMatricesMaps(variableTypeIdx)%dofToRowsMap( &
              & dependentVariable%TOTAL_NUMBER_OF_DOFS),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate variable to equations matrices maps dof to rows map.",err,error,*999)
            linearMapping%varToEquationsMatricesMaps(variableTypeIdx)%dofToRowsMap=0
            linearMapping%numberOfLinearMatrixVariables=linearMapping%numberOfLinearMatrixVariables+1
          ELSE IF(linearMapping%varToEquationsMatricesMaps(variableTypeIdx)%numberOfEquationsMatrices>0) THEN
            ALLOCATE(linearMapping%varToEquationsMatricesMaps(variableTypeIdx)%equationsMatrixNumbers( &
              & linearMapping%varToEquationsMatricesMaps(variableTypeIdx)%numberOfEquationsMatrices),STAT=err)
            IF(err/=0) &
              & CALL FlagError("Could not allocate variable to equations matrices maps equations matrix numbers.",err,error,*999)
            ALLOCATE(linearMapping%varToEquationsMatricesMaps(variableTypeIdx)%dofToColumnsMaps( &
              & linearMapping%varToEquationsMatricesMaps(variableTypeIdx)%numberOfEquationsMatrices),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate variable to equations matrices maps dof to columns map.",err,error,*999)
            linearMapping%varToEquationsMatricesMaps(variableTypeIdx)%equationsMatrixNumbers=0
            linearMapping%varToEquationsMatricesMaps(variableTypeIdx)%numberOfEquationsMatrices=0
            DO matrixIdx=1,linearMapping%numberOfLinearMatrices
              IF(createValuesCache%linearMatrixVariableTypes(matrixIdx)==variableTypeIdx) THEN
                linearMapping%varToEquationsMatricesMaps(variableTypeIdx)%numberOfEquationsMatrices= &
                  & linearMapping%varToEquationsMatricesMaps(variableTypeIdx)%numberOfEquationsMatrices+1
                linearMapping%varToEquationsMatricesMaps(variableTypeIdx)%equationsMatrixNumbers( &
                  & linearMapping%varToEquationsMatricesMaps(variableTypeIdx)%numberOfEquationsMatrices)=matrixIdx
                ALLOCATE(linearMapping%varToEquationsMatricesMaps(variableTypeIdx)%dofToColumnsMaps( &
                  & linearMapping%varToEquationsMatricesMaps(variableTypeIdx)% &
                  & numberOfEquationsMatrices)%columnDOF(dependentVariable%TOTAL_NUMBER_OF_DOFS),STAT=err)
                IF(err/=0) CALL FlagError("Could not allocate variable dof to columns map column dof.",err,error,*999)
                DO dofIdx=1,dependentVariable%TOTAL_NUMBER_OF_DOFS
                  !1-1 mapping for now
                  columnIdx=dependentVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(dofIdx)
                  linearMapping%varToEquationsMatricesMaps(variableTypeIdx)%dofToColumnsMaps( &
                    & linearMapping%varToEquationsMatricesMaps(variableTypeIdx)% &
                    & numberOfEquationsMatrices)%columnDOF(dofIdx)=columnIdx
                ENDDO !dofIdx
              ENDIF
            ENDDO !matrixIdx
            ALLOCATE(linearMapping%varToEquationsMatricesMaps(variableTypeIdx)%dofToRowsMap( &
              & dependentVariable%TOTAL_NUMBER_OF_DOFS),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate variable to equations matrices maps dof to rows map.",err,error,*999)
            DO dofIdx=1,dependentVariable%TOTAL_NUMBER_OF_DOFS
              !1-1 mappings for now.
              rowIdx=dofIdx
              linearMapping%varToEquationsMatricesMaps(variableTypeIdx)%dofToRowsMap(dofIdx)=rowIdx
            ENDDO !dofIdx
            linearMapping%numberOfLinearMatrixVariables=linearMapping%numberOfLinearMatrixVariables+1
          ENDIF
        ENDIF
      ENDDO !variableTypeIdx
      !Allocate and initialise the variable types
      ALLOCATE(linearMapping%linearMatrixVariableTypes(linearMapping%numberOfLinearMatrixVariables),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate equations mapping matrix variable types.",err,error,*999)
      linearMapping%numberOfLinearMatrixVariables=0
      DO variableTypeIdx=1,FIELD_NUMBER_OF_VARIABLE_TYPES
        IF(linearMapping%varToEquationsMatricesMaps(variableTypeIdx)%numberOfEquationsMatrices>0) THEN
          linearMapping%numberOfLinearMatrixVariables=linearMapping%numberOfLinearMatrixVariables+1
          linearMapping%linearMatrixVariableTypes(linearMapping%numberOfLinearMatrixVariables)=variableTypeIdx
        ENDIF
      ENDDO !variableTypeIdx
      !Allocate and initialise the equations matrix to variable maps types
      ALLOCATE(linearMapping%equationsMatrixToVarMaps(linearMapping%numberOfLinearMatrices),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate equations mapping equations matrix to variable maps.",err,error,*999)
      !Create the individual matrix maps and column maps
      DO matrixIdx=1,linearMapping%numberOfLinearMatrices
        variableType=createValuesCache%linearMatrixVariableTypes(matrixIdx)
        dependentVariable=>dependentField%VARIABLE_TYPE_MAP(variableType)%ptr
        CALL EquationsMapping_EquatsMatrixToVarMapInitialise(linearMapping%equationsMatrixToVarMaps(matrixIdx),err,error,*999)
        linearMapping%equationsMatrixToVarMaps(matrixIdx)%matrixNumber=matrixIdx
        linearMapping%equationsMatrixToVarMaps(matrixIdx)%variableType=variableType
        linearMapping%equationsMatrixToVarMaps(matrixIdx)%variable=>dependentVariable
        linearMapping%equationsMatrixToVarMaps(matrixIdx)%numberOfColumns=dependentVariable%DOMAIN_MAPPING%NUMBER_OF_GLOBAL
        linearMapping%equationsMatrixToVarMaps(matrixIdx)%matrixCoefficient=vectorMapping% &
          createValuesCache%linearMatrixCoefficients(matrixIdx)
        ALLOCATE(linearMapping%equationsMatrixToVarMaps(matrixIdx)%columnToDOFMap( &
          & dependentVariable%DOMAIN_MAPPING%NUMBER_OF_GLOBAL),STAT=err)                  
        IF(err/=0) CALL FlagError("Could not allocate equation matrix to variable map column to dof map.",err,error,*999)
        linearMapping%equationsMatrixToVarMaps(matrixIdx)%columnToDOFMap=0
        DO dofIdx=1,dependentVariable%TOTAL_NUMBER_OF_DOFS
          !1-1 mapping for now
          columnIdx=dependentVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(dofIdx)
          linearMapping%equationsMatrixToVarMaps(matrixIdx)%columnToDOFMap(columnIdx)=dofIdx
        ENDDO !dofIdx
        linearMapping%equationsMatrixToVarMaps(matrixIdx)%columnDOFSMapping=>dependentVariable%DOMAIN_MAPPING
      ENDDO !matrixIdx
      !Allocate the row mappings
      ALLOCATE(linearMapping%equationsRowToVariableDOFMaps(vectorMapping%totalNumberOfRows, &
        & linearMapping%numberOfLinearMatrixVariables),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate equations row to variable dof maps.",err,error,*999)
      !Set up the row mappings
      DO variableIdx=1,linearMapping%numberOfLinearMatrixVariables
        DO rowIdx=1,vectorMapping%totalNumberOfRows
          !1-1 mapping for now
          linearMapping%equationsRowToVariableDOFMaps(rowIdx,variableIdx)=rowIdx
        ENDDO !rowIdx
      ENDDO !variableIdx
    ENDIF
    !Calculate non-linear mappings
    IF(createValuesCache%numberOfResidualVariables/=0) THEN
      CALL EquationsMapping_NonlinearMappingInitialise(vectorMapping,err,error,*999)
      NULLIFY(nonlinearMapping)
      CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)
      nonlinearMapping%numberOfResidualVariables=createValuesCache%numberOfResidualVariables
      ALLOCATE(nonlinearMapping%varToJacobianMap(nonlinearMapping%numberOfResidualVariables),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate variable to Jacobian maps.",err,error,*999)
      ALLOCATE(nonlinearMapping%jacobianToVarMap(nonlinearMapping%numberOfResidualVariables),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate Jacobian to variable maps.",err,error,*999)
      ALLOCATE(nonlinearMapping%residualVariables(nonlinearMapping%numberOfResidualVariables),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate nonlinear mapping residual variables.",err,error,*999)
      DO matrixIdx=1,nonlinearMapping%numberOfResidualVariables
        CALL EquationsMapping_VarToEquatsJacobianMapInitialise(nonlinearMapping%varToJacobianMap(matrixIdx),err,error,*999)
        nonlinearMapping%varToJacobianMap(matrixIdx)%jacobianNumber=matrixIdx
        nonlinearMapping%varToJacobianMap(matrixIdx)%variableType=createValuesCache%residualVariableTypes(matrixIdx)
        dependentVariable=>dependentField%VARIABLE_TYPE_MAP(createValuesCache%residualVariableTypes(matrixIdx))%ptr
        nonlinearMapping%varToJacobianMap(matrixIdx)%variable=>dependentVariable
        nonlinearMapping%residualVariables(matrixIdx)%ptr=>dependentVariable
        !Row variable is RHS if set, otherwise first nonlinear variable
        IF(createValuesCache%rhsVariableType/=0) THEN
          rowVariable=>dependentField%VARIABLE_TYPE_MAP(createValuesCache%rhsVariableType)%ptr
        ELSE
          rowVariable=>dependentField%VARIABLE_TYPE_MAP(createValuesCache%residualVariableTypes(1))%ptr
        ENDIF
        !Allocate and set dof to Jacobian columns map
        ALLOCATE(nonlinearMapping%varToJacobianMap(matrixIdx)%dofToColumnsMap(dependentVariable%TOTAL_NUMBER_OF_DOFS),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate variable to Jacobian map dof to columns map.",err,error,*999)
        DO dofIdx=1,dependentVariable%TOTAL_NUMBER_OF_DOFS
          !1-1 mapping for now
          columnIdx=dependentVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(dofIdx)
          nonlinearMapping%varToJacobianMap(matrixIdx)%dofToColumnsMap(dofIdx)=columnIdx
        ENDDO !dofIdx
        !Allocate and set dof to Jacobian rows map
        ALLOCATE(nonlinearMapping%varToJacobianMap(matrixIdx)%dofToRowsMap(rowVariable%TOTAL_NUMBER_OF_DOFS),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate variable to Jacobian map dof to columns map.",err,error,*999)
        DO dofIdx=1,rowVariable%TOTAL_NUMBER_OF_DOFS
          !1-1 mapping for now
          rowIdx=dofIdx
          nonlinearMapping%varToJacobianMap(matrixIdx)%dofToRowsMap(dofIdx)=rowIdx
        ENDDO !dofIdx
        CALL EquationsMapping_EquatsJacobianToVarMapInitialise(nonlinearMapping%jacobianToVarMap(matrixIdx),err,error,*999)
        nonlinearMapping%jacobianToVarMap(matrixIdx)%jacobianNumber=matrixIdx
        nonlinearMapping%jacobianToVarMap(matrixIdx)%variableType=createValuesCache%residualVariableTypes(matrixIdx)
        nonlinearMapping%jacobianToVarMap(matrixIdx)%variable=>dependentVariable
        nonlinearMapping%jacobianToVarMap(matrixIdx)%numberOfColumns=dependentVariable%DOMAIN_MAPPING%NUMBER_OF_GLOBAL
        nonlinearMapping%jacobianToVarMap(matrixIdx)%jacobianCoefficient=createValuesCache%residualCoefficient
        ALLOCATE(nonlinearMapping%jacobianToVarMap(matrixIdx)%equationsColumnToDOFVariableMap( &
          & dependentVariable%DOMAIN_MAPPING%NUMBER_OF_GLOBAL),STAT=err)
        nonlinearMapping%jacobianToVarMap(matrixIdx)%equationsColumnToDOFVariableMap=0
        DO dofIdx=1,dependentVariable%TOTAL_NUMBER_OF_DOFS
          !1-1 mapping for now
          columnIdx=dependentVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(dofIdx)
          nonlinearMapping%jacobianToVarMap(matrixIdx)%equationsColumnToDOFVariableMap(columnIdx)=dofIdx
        ENDDO !dofIdx
        nonlinearMapping%jacobianToVarMap(matrixIdx)%columnDOFSMapping=>dependentVariable%DOMAIN_MAPPING
      ENDDO !matrixIdx
      !Set up the row mappings
      ALLOCATE(nonlinearMapping%equationsRowToResidualDOFMap(totalNumberOfRows),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate equations row to residual dof map.",err,error,*999)
      DO rowIdx=1,totalNumberOfRows
        !1-1 mapping for now
        dofIdx=rowIdx
        nonlinearMapping%equationsRowToResidualDOFMap(rowIdx)=dofIdx
      ENDDO !rowIdx
    ENDIF
    !Calculate RHS mappings
    IF(createValuesCache%rhsVariableType/=0) THEN                  
      CALL EquationsMapping_RHSMappingInitialise(vectorMapping,err,error,*999)
      NULLIFY(rhsMapping)
      CALL EquationsMappingVector_RHSMappingGet(vectorMapping,rhsMapping,err,error,*999)
      rhsMapping%rhsVariableType=createValuesCache%rhsVariableType
      dependentVariable=>dependentField%VARIABLE_TYPE_MAP(createValuesCache%rhsVariableType)%ptr
      rhsMapping%rhsVariable=>dependentVariable
      rhsMapping%rhsVariableMapping=>dependentVariable%DOMAIN_MAPPING
      rhsMapping%rhsCoefficient=createValuesCache%rhsCoefficient
      !Allocate and set up the row mappings
      ALLOCATE(rhsMapping%rhsDOFToEquationsRowMap(dependentVariable%TOTAL_NUMBER_OF_DOFS),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate rhs dof to equations row map.",err,error,*999)
      ALLOCATE(rhsMapping%equationsRowToRHSDOFMap(totalNumberOfRows),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate equations row to dof map.",err,error,*999)
      DO dofIdx=1,dependentVariable%TOTAL_NUMBER_OF_DOFS
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
    IF(createValuesCache%sourceVariableType/=0) THEN                  
      CALL EquationsMapping_SourceMappingInitialise(vectorMapping,err,error,*999)
      NULLIFY(sourceMapping)
      CALL EquationsMappingVector_SourceMappingGet(vectorMapping,sourceMapping,err,error,*999)
      sourceMapping%sourceVariableType=createValuesCache%sourceVariableType
      sourceVariable=>sourceField%VARIABLE_TYPE_MAP(createValuesCache%sourceVariableType)%ptr
      sourceMapping%sourceVariable=>sourceVariable
      !    sourceMapping%sourceVariableMapping=>sourceVariable%DOMAIN_MAPPING
      !    sourceMapping%sourceCoefficient=createValuesCache%sourceCoefficient
      !    !Allocate and set up the row mappings
      !    ALLOCATE(sourceMapping%sourceDOFToEquationsRowMap(sourceVariable%TOTAL_NUMBER_OF_DOFS),STAT=err)
      !    IF(err/=0) CALL FlagError("Could not allocate source dof to equations row map.",err,error,*999)
      !    ALLOCATE(sourceMapping%equationsRowToSourceDOFMap(totalNumberOfRows),STAT=err)
      !    IF(err/=0) CALL FlagError("Could not allocate equations row to source map.",err,error,*999)
      !    DO dofIdx=1,sourceVariable%TOTAL_NUMBER_OF_DOFS
      !      !1-1 mapping for now
      !      rowIdx=dofIdx
      !      sourceMapping%sourceDOFToEquationsRowMap(dofIdx)=rowIdx
      !    ENDDO !dofIdx
      !    DO rowIdx=1,totalNumberOfRows
      !      !1-1 mapping for now
      !      dofIdx=rowIdx
      !      sourceMapping%equationsRowToSourceDOFMap(rowIdx)=dofIdx
      !    ENDDO !rowIdx
    ENDIF
    
    IF(diagnostics1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Equations vector mappings:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of rows = ",vectorMapping%numberOfRows,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Total umber of rows = ",vectorMapping%totalNumberOfRows,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of global rows = ",vectorMapping%numberOfGlobalRows,err,error,*999)
      dynamicMapping=>vectorMapping%dynamicMapping
      IF(ASSOCIATED(dynamicMapping)) THEN
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
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Variable to matrices mappings:",err,error,*999)
        DO variableTypeIdx=1,FIELD_NUMBER_OF_VARIABLE_TYPES
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Variable type : ",variableTypeIdx,err,error,*999)
          IF(ASSOCIATED(dynamicMapping%varToEquationsMatricesMaps(variableTypeIdx)%variable)) THEN
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Total number of DOFs = ",dynamicMapping% &
              & varToEquationsMatricesMaps(variableTypeIdx)%variable%TOTAL_NUMBER_OF_DOFS,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of equations matrices = ",dynamicMapping% &
              & varToEquationsMatricesMaps(variableTypeIdx)%numberOfEquationsMatrices,err,error,*999)
            IF(dynamicMapping%varToEquationsMatricesMaps(variableTypeIdx)%numberOfEquationsMatrices>0) THEN
              CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,dynamicMapping%varToEquationsMatricesMaps(variableTypeIdx)% &
                & numberOfEquationsMatrices,4,4,dynamicMapping%varToEquationsMatricesMaps(variableTypeIdx)% &
                & equationsMatrixNumbers,'("      Matrix numbers :",4(X,I12))','(22X,4(X,I12))',err,error,*999) 
              CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"      DOF to column maps :",err,error,*999)
              DO matrixIdx=1,dynamicMapping%varToEquationsMatricesMaps(variableTypeIdx)%numberOfEquationsMatrices
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Matrix number : ",matrixIdx,err,error,*999)
                CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,dynamicMapping%varToEquationsMatricesMaps( &
                  & variableTypeIdx)%variable%TOTAL_NUMBER_OF_DOFS,5,5,dynamicMapping%varToEquationsMatricesMaps( &
                  & variableTypeIdx)%dofToColumnsMaps(matrixIdx)%columnDOF, &
                  & '("        Column numbers :",5(X,I13))','(24X,5(X,I13))',err,error,*999) 
              ENDDO !matrixIdx
              CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,dynamicMapping%varToEquationsMatricesMaps(variableTypeIdx)% &
                & VARIABLE%TOTAL_NUMBER_OF_DOFS,5,5,dynamicMapping%varToEquationsMatricesMaps(variableTypeIdx)% &
                & dofToRowsMap,'("      DOF to row maps  :",5(X,I13))','(24X,5(X,I13))',err,error,*999)
            ENDIF
          ENDIF
        ENDDO !variableTypeIdx
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Matrix to variable mappings:",err,error,*999)
        DO matrixIdx=1,dynamicMapping%numberOfDynamicMatrices
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Matrix number : ",matrixIdx,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Variable type = ",dynamicMapping% &
            & equationsMatrixToVarMaps(matrixIdx)%variableType,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Number of columns = ",dynamicMapping% &
            & equationsMatrixToVarMaps(matrixIdx)%numberOfColumns,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Matrix coefficient = ",dynamicMapping% &
            & equationsMatrixToVarMaps(matrixIdx)%matrixCoefficient,err,error,*999)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,dynamicMapping%equationsMatrixToVarMaps(matrixIdx)% &
            & numberOfColumns,5,5,dynamicMapping%equationsMatrixToVarMaps(matrixIdx)%columnToDOFMap, &
            & '("        Column to DOF maps :",5(X,I13))','(28X,5(X,I13))',err,error,*999) 
        ENDDO !matrixIdx
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Row mappings:",err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,vectorMapping%totalNumberOfRows,5,5, &
          & dynamicMapping%equationsRowToVariableDOFMaps,'("    Row to DOF maps :",5(X,I13))','(21X,5(X,I13))', &
          & err,error,*999) 
      ENDIF
      linearMapping=>vectorMapping%linearMapping
      IF(ASSOCIATED(linearMapping)) THEN
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Linear mappings:",err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of linear equations matrices = ",linearMapping% &
          & numberOfLinearMatrices,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of linear matrix variables = ",linearMapping% &
          & numberOfLinearMatrixVariables,err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,linearMapping%numberOfLinearMatrixVariables,4,4, &
          & linearMapping%linearMatrixVariableTypes,'("    Linear matrix variable types :",4(X,I12))','(27X,4(X,I12))', &
          & err,error,*999) 
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Variable to matrices mappings:",err,error,*999)
        DO variableTypeIdx=1,FIELD_NUMBER_OF_VARIABLE_TYPES
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Variable type : ",variableTypeIdx,err,error,*999)
          IF(ASSOCIATED(linearMapping%varToEquationsMatricesMaps(variableTypeIdx)%variable)) THEN
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Total number of DOFs = ",linearMapping% &
              & varToEquationsMatricesMaps(variableTypeIdx)%variable%TOTAL_NUMBER_OF_DOFS,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of equations matrices = ",linearMapping% &
              & varToEquationsMatricesMaps(variableTypeIdx)%numberOfEquationsMatrices,err,error,*999)
            IF(linearMapping%varToEquationsMatricesMaps(variableTypeIdx)%numberOfEquationsMatrices>0) THEN
              CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,linearMapping%varToEquationsMatricesMaps(variableTypeIdx)% &
                & numberOfEquationsMatrices,4,4,linearMapping%varToEquationsMatricesMaps(variableTypeIdx)% &
                & equationsMatrixNumbers,'("      Matrix numbers :",4(X,I12))','(22X,4(X,I12))',err,error,*999) 
              CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"      DOF to column maps :",err,error,*999)
              DO matrixIdx=1,linearMapping%varToEquationsMatricesMaps(variableTypeIdx)%numberOfEquationsMatrices
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Matrix number : ",matrixIdx,err,error,*999)
                CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,linearMapping%varToEquationsMatricesMaps( &
                  & variableTypeIdx)%variable%TOTAL_NUMBER_OF_DOFS,5,5,linearMapping%varToEquationsMatricesMaps( &
                  & variableTypeIdx)%dofToColumnsMaps(matrixIdx)%columnDOF, &
                  & '("        Column numbers :",5(X,I13))','(24X,5(X,I13))',err,error,*999) 
              ENDDO !matrixIdx
              CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,linearMapping%varToEquationsMatricesMaps(variableTypeIdx)% &
                & VARIABLE%TOTAL_NUMBER_OF_DOFS,5,5,linearMapping%varToEquationsMatricesMaps(variableTypeIdx)% &
                & dofToRowsMap,'("      DOF to row maps  :",5(X,I13))','(24X,5(X,I13))',err,error,*999)
            ENDIF
          ENDIF
        ENDDO !variableTypeIdx
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Matrix to variable mappings:",err,error,*999)
        DO matrixIdx=1,linearMapping%numberOfLinearMatrices
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Matrix number : ",matrixIdx,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Variable type = ",linearMapping% &
            & equationsMatrixToVarMaps(matrixIdx)%variableType,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Number of columns = ",linearMapping% &
            & equationsMatrixToVarMaps(matrixIdx)%numberOfColumns,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Matrix coefficient = ",linearMapping% &
            & equationsMatrixToVarMaps(matrixIdx)%matrixCoefficient,err,error,*999)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,linearMapping%equationsMatrixToVarMaps(matrixIdx)% &
            & numberOfColumns,5,5,linearMapping%equationsMatrixToVarMaps(matrixIdx)%columnToDOFMap, &
            & '("        Column to DOF maps :",5(X,I13))','(28X,5(X,I13))',err,error,*999) 
        ENDDO !matrixIdx
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Row mappings:",err,error,*999)
        DO variableIdx=1,linearMapping%numberOfLinearMatrixVariables
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Variable number : ",variableIdx,err,error,*999)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,vectorMapping%totalNumberOfRows,5,5, &
            & linearMapping%equationsRowToVariableDOFMaps(:,variableIdx), &
            & '("    Row to DOF maps :",5(X,I13))','(21X,5(X,I13))',err,error,*999) 
        ENDDO !variableIdx
      ENDIF
      nonlinearMapping=>vectorMapping%nonlinearMapping
      IF(ASSOCIATED(nonlinearMapping)) THEN
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Nonlinear mappings:",err,error,*999)
        DO matrixIdx=1,nonlinearMapping%numberOfResidualVariables
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Matrix number : ",matrixIdx,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Residual variable type = ",nonlinearMapping% &
            & jacobianToVarMap(matrixIdx)%variableType,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Total number of residual DOFs = ",nonlinearMapping% &
            & jacobianToVarMap(matrixIdx)%variable%TOTAL_NUMBER_OF_DOFS,err,error,*999)
        ENDDO
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Residual coefficient = ",nonlinearMapping%residualCoefficient, &
          & err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Residual row mappings:",err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,vectorMapping%totalNumberOfRows,5,5, &
          & nonlinearMapping%equationsRowToResidualDOFMap,'("    Row to DOF mappings :",5(X,I13))','(25X,5(X,I13))', &
          & err,error,*999) 
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Jacobian mappings:",err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Variable to Jacobian mappings:",err,error,*999)
        DO matrixIdx=1,nonlinearMapping%numberOfResidualVariables
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Matrix number : ",matrixIdx,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Jacobian variable type = ",nonlinearMapping% &
            & varToJacobianMap(matrixIdx)%variableType,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Total number of Jacobain DOFs = ",nonlinearMapping% &
            & varToJacobianMap(matrixIdx)%variable%TOTAL_NUMBER_OF_DOFS,err,error,*999)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,nonlinearMapping%varToJacobianMap(matrixIdx)%variable% &
            & TOTAL_NUMBER_OF_DOFS,5,5,nonlinearMapping%varToJacobianMap(matrixIdx)%dofToColumnsMap, &
            & '("      DOF to column map :",5(X,I13))','(26X,5(X,I13))',err,error,*999) 
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,nonlinearMapping%varToJacobianMap(matrixIdx)%variable% &
            & TOTAL_NUMBER_OF_DOFS,5,5,nonlinearMapping%varToJacobianMap(matrixIdx)%dofToRowsMap, &
            & '("      DOF to row map    :",5(X,I13))','(26X,5(X,I13))',err,error,*999) 
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Jacobian to variable mappings:",err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Jacobian variable type = ",nonlinearMapping% &
            & jacobianToVarMap(matrixIdx)%variableType,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of columns = ",nonlinearMapping% &
            & jacobianToVarMap(matrixIdx)%numberOfColumns,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Jacobian coefficient = ",nonlinearMapping% &
            & jacobianToVarMap(matrixIdx)%jacobianCoefficient,err,error,*999)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,nonlinearMapping%jacobianToVarMap(matrixIdx)%numberOfColumns, &
            & 5,5,nonlinearMapping%jacobianToVarMap(matrixIdx)%equationsColumnToDOFVariableMap, &
            & '("      Column to DOF map :",5(X,I13))','(26X,5(X,I13))',err,error,*999) 
        ENDDO
      ENDIF
      rhsMapping=>vectorMapping%rhsMapping
      IF(ASSOCIATED(rhsMapping)) THEN
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  RHS mappings:",err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    RHS variable type = ",rhsMapping%rhsVariableType,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Total number of RHS DOFs = ",rhsMapping%rhsVariable% &
          & TOTAL_NUMBER_OF_DOFS,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    RHS coefficient = ",rhsMapping%rhsCoefficient,err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Row mappings:",err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,rhsMapping%rhsVariable%TOTAL_NUMBER_OF_DOFS,5,5, &
          & rhsMapping%rhsDOFToEquationsRowMap,'("    DOF to row mappings :",5(X,I13))','(25X,5(X,I13))',err,error,*999) 
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,vectorMapping%totalNumberOfRows,5,5, &
          & rhsMapping%equationsRowToRHSDOFMap,'("    Row to DOF mappings :",5(X,I13))','(25X,5(X,I13))',err,error,*999) 
       ENDIF
      sourceMapping=>vectorMapping%sourceMapping
      IF(ASSOCIATED(sourceMapping)) THEN
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Source mappings:",err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Source variable type = ",sourceMapping%sourceVariableType, &
          & err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Total number of source DOFs = ",sourceMapping%sourceVariable% &
          & TOTAL_NUMBER_OF_DOFS,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Source coefficient = ",sourceMapping%sourceCoefficient,err,error,*999)
        !CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Row mappings:",err,error,*999)
        !CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,sourceMapping%sourceVariable%TOTAL_NUMBER_OF_DOFS,5,5, &
        !  & sourceMapping%sourceDOFToEquationsRowMap,'("    DOF to row mappings :",5(X,I13))','(25X,5(X,I13))', &
        !  & err,error,*999) 
        !CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,vectorMapping%totalNumberOfRows,5,5, &
        !  & sourceMapping%equationsRowToSourceDOFMap,'("    Row to DOF mappings :",5(X,I13))','(25X,5(X,I13))', &
        !  & err,error,*999) 
      ENDIF
    ENDIF
       
    EXITS("EquationsMapping_VectorCalculate")
    RETURN
999 ERRORSEXITS("EquationsMapping_VectorCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMapping_VectorCalculate

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
    INTEGER(INTG) :: matrixIdx
    LOGICAL :: isResidualType
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsMapping_VectorCreateFinish",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector equations mapping is not associated.",err,error,*999)
    IF(vectorMapping%vectorMappingFinished) CALL FlagError("Equations mapping has already been finished.",err,error,*999)

    NULLIFY(createValuesCache)
    CALL EquationsMappingVector_CreateValuesCacheGet(vectorMapping,createValuesCache,err,error,*999)
    NULLIFY(vectorEquations)
    CALL EquationsMappingVector_VectorEquationsGet(vectorMapping,vectorEquations,err,error,*999)
    NULLIFY(equations)
    CALL EquationsVector_EquationsGet(vectorEquations,equations,err,error,*999)
    NULLIFY(equationsSet)
    CALL Equations_EquationsSetGet(equations,equationsSet,err,error,*999)

    !Check that all the variables have been mapped properly
    SELECT CASE(equations%timeDependence)
    CASE(EQUATIONS_STATIC,EQUATIONS_QUASISTATIC)
      SELECT CASE(equations%linearity)
      CASE(EQUATIONS_LINEAR,EQUATIONS_NONLINEAR_BCS)
        IF(createValuesCache%rhsVariableType==0.AND.createValuesCache%numberOfLinearMatrices==0) &
          & CALL FlagError("Invalid equations mapping. The RHS variable type must be set if there are no "// &
          & "linear matrices.",err,error,*999)
      CASE(EQUATIONS_NONLINEAR)
        DO matrixIdx=1,createValuesCache%numberOfResidualVariables
          IF(createValuesCache%residualVariableTypes(matrixIdx)==0) THEN
            localError="Invalid equations mapping. The residual variable type is not set for Jacobian number "// &
              & TRIM(NumberToVString(matrixIdx,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ENDDO !matrixIdx
        IF(createValuesCache%rhsVariableType==0.AND.createValuesCache%numberOfLinearMatrices==0) &
          & CALL FlagError("Invalid equations mapping. The RHS variable type must be set if there are no "// &
          & "linear matrices.",err,error,*999)
      CASE DEFAULT
        localError="The equations linearity type of "//TRIM(NumberToVString(equations%linearity,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_FIRST_ORDER_DYNAMIC,EQUATIONS_SECOND_ORDER_DYNAMIC)
      SELECT CASE(equations%linearity)
      CASE(EQUATIONS_LINEAR,EQUATIONS_NONLINEAR_BCS)
        IF(createValuesCache%dynamicVariableType==0) CALL FlagError("Invalid equations mapping. "// &
          & "The dynamic variable type must be set for dynamic equations.", err,error,*999)
        IF(createValuesCache%rhsVariableType==0.AND.createValuesCache%numberOfLinearMatrices==0) &
          & CALL FlagError("Invalid equations mapping. The RHS variable type must be set if there are no "// &
          & "linear matrices.",err,error,*999)
      CASE(EQUATIONS_NONLINEAR)
        ! SEBK 19/08/2009 not sure about mapping here
        !|
        IF(createValuesCache%dynamicVariableType==0) CALL FlagError("Invalid equations mapping. "// &
          & "The dynamic variable type must be set for dynamic equations.", err,error,*999)
        IF(createValuesCache%rhsVariableType==0.AND.createValuesCache%numberOfLinearMatrices==0) &
          & CALL FlagError("Invalid equations mapping. The RHS variable type must be set if there are no "// &
          & "linear matrices.",err,error,*999)
        isResidualType=.FALSE.
        DO matrixIdx=1,createValuesCache%numberOfResidualVariables
          IF(createValuesCache%residualVariableTypes(matrixIdx)==0) THEN
            localError="Invalid equations mapping. The residual variable type is not set for Jacobian number "// &
              & TRIM(NumberToVString(matrixIdx,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
          IF(createValuesCache%residualVariableTypes(matrixIdx)==createValuesCache%dynamicVariableType) THEN
            isResidualType=.TRUE.
          ENDIF
        ENDDO !matrixIdx
        IF(.NOT.isResidualType) THEN
          CALL FlagError("Invalid equations mapping. The residual variable type must correspond to the "// &
            & "dynamic variable type for nonlinear dynamic equations.", err,error,*999)
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
    !Check the linear matrices variable types
    DO matrixIdx=1,createValuesCache%numberOfLinearMatrices
      IF(createValuesCache%linearMatrixVariableTypes(matrixIdx)==0) THEN
        localError="Invalid equations mapping. The linear matrix variable type is not set for linear matrix number "//&
          & TRIM(NumberToVString(matrixIdx,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDDO !matrixIdx
    !Now calculate the equations mapping and clean up
    CALL EquationsMapping_VectorCalculate(vectorMapping,err,error,*999)
    CALL EquationsMapping_VectorCreateValuesCacheFinalise(vectorMapping%createValuesCache,err,error,*999)
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
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    TYPE(FIELD_TYPE), POINTER :: dependentField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: lhsVariable
    
    ENTERS("EquationsMapping_VectorCreateStart",err,error,*998)

    IF(ASSOCIATED(vectorMapping)) CALL FlagError("Vector equations mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorEquations)) CALL FlagError("Vector equations is not associated.",err,error,*999)
    NULLIFY(equations)
    CALL EquationsVector_EquationsGet(vectorEquations,equations,err,error,*999)
    NULLIFY(equationsSet)
    CALL Equations_EquationsSetGet(equations,equationsSet,err,error,*999)
    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
    IF(.NOT.equations%equationsFinished) CALL FlagError("Vector equations equations has not been finished.",err,error,*999)

    NULLIFY(lhsVariable)
    CALL Field_VariableGet(dependentField,lhsVariableType,lhsVariable,err,error,*999)
    IF(.NOT.ASSOCIATED(lhsVariable%DOMAIN_MAPPING))  &
      & CALL FlagError("LHS variable domain mapping is not associated.",err,error,*999)
     
    CALL EquationsMapping_VectorInitialise(vectorEquations,err,error,*999)
    CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
    NULLIFY(lhsMapping)
    CALL EquationsMappingVector_LHSMappingGet(vectorMapping,lhsMapping,err,error,*999)
    lhsMapping%lhsVariableType=lhsVariableType
    lhsMapping%lhsVariable=>lhsVariable
    lhsMapping%rowDofsMapping=>lhsVariable%DOMAIN_MAPPING
    lhsMapping%numberOfRows=lhsVariable%NUMBER_OF_DOFS
    lhsMapping%totalNumberOfRows=lhsVariable%TOTAL_NUMBER_OF_DOFS
    lhsMapping%numberOfGlobalRows=lhsVariable%NUMBER_OF_GLOBAL_DOFS
   
    EXITS("EquationsMapping_VectorCreateStart")
    RETURN
999 NULLIFY(vectorMapping)
998 ERRORSEXITS("EquationsMapping_VectorCreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMapping_VectorCreateStart

  !
  !================================================================================================================================
  !

  !>Finalises an equations mapping create values cache and deallocates all memory
  SUBROUTINE EquationsMapping_VectorCreateValuesCacheFinalise(createValuesCache,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorCreateValuesCacheType), POINTER :: createValuesCache !<A pointer to the create values cache to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsMapping_VectorCreateValuesCacheFinalise",err,error,*999)

    IF(ASSOCIATED(createValuesCache)) THEN
      IF(ALLOCATED(createValuesCache%dynamicMatrixCoefficients)) DEALLOCATE(createValuesCache%dynamicMatrixCoefficients)
      IF(ALLOCATED(createValuesCache%linearMatrixVariableTypes)) DEALLOCATE(createValuesCache%linearMatrixVariableTypes)
      IF(ALLOCATED(createValuesCache%linearMatrixCoefficients)) DEALLOCATE(createValuesCache%linearMatrixCoefficients)
      IF(ALLOCATED(createValuesCache%residualVariableTypes)) DEALLOCATE(createValuesCache%residualVariableTypes)
      DEALLOCATE(createValuesCache)
    ENDIF
       
    EXITS("EquationsMapping_VectorCreateValuesCacheFinalise")
    RETURN
999 ERRORS("EquationsMapping_VectorCreateValuesCacheFinalise",err,error)
    EXITS("EquationsMapping_VectorCreateValuesCacheFinalise")
    RETURN 1
    
  END SUBROUTINE EquationsMapping_VectorCreateValuesCacheFinalise

  !
  !================================================================================================================================
  !

  !>Initialises a vector equations mapping create values cache 
  SUBROUTINE EquationsMapping_VectorCreateValuesCacheInitialise(vectorMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the vector equations mapping to initialise the create values cache for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,matrixIdx,matrixIdx2,variableNumber
    LOGICAL :: isResidualType
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    TYPE(FIELD_TYPE), POINTER :: dependentField
    TYPE(VARYING_STRING) :: dummyError,localError

    ENTERS("EquationsMapping_VectorCreateValuesCacheInitialise",err,error,*998)

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
    vectorMapping%createValuesCache%numberOfResidualVariables=0
    vectorMapping%createValuesCache%residualCoefficient=1.0_DP
    vectorMapping%createValuesCache%rhsVariableType=0
    vectorMapping%createValuesCache%rhsCoefficient=1.0_DP
    vectorMapping%createValuesCache%sourceVariableType=0
    vectorMapping%createValuesCache%sourceCoefficient=1.0_DP
    !Set the default equations mapping in the create values cache
    !First calculate how many linear and dynamic matrices we have and set the variable types for the dynamic, residual
    !and RHS variables
    IF(dependentField%NUMBER_OF_VARIABLES==1) THEN
      SELECT CASE(equations%linearity)
      CASE(EQUATIONS_LINEAR,EQUATIONS_NONLINEAR_BCS)
        CALL FlagError("Dependent field only has one variable which cannot be mapped to both an equations matrix and RHS vector.", &
          & err,error,*999)
      CASE(EQUATIONS_NONLINEAR)
        CALL FlagError("Dependent field only has one variable which cannot be mapped to both the residual and RHS vector.", &
          & err,error,*999)
      CASE DEFAULT
        localError="The equations linearity type of "//TRIM(NumberToVString(equations%linearity,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE IF(dependentField%NUMBER_OF_VARIABLES>1) THEN
      SELECT CASE(equations%timeDependence)
      CASE(EQUATIONS_STATIC,EQUATIONS_QUASISTATIC)
        SELECT CASE(equations%linearity)
        CASE(EQUATIONS_LINEAR,EQUATIONS_NONLINEAR_BCS)
          vectorMapping%createValuesCache%numberOfLinearMatrices=dependentField%NUMBER_OF_VARIABLES-1
          IF(ASSOCIATED(dependentField%VARIABLE_TYPE_MAP(FIELD_DELUDELN_VARIABLE_TYPE)%ptr)) THEN
            vectorMapping%createValuesCache%rhsVariableType=dependentField%VARIABLE_TYPE_MAP(FIELD_DELUDELN_VARIABLE_TYPE)% &
              & PTR%VARIABLE_TYPE
          ELSE
            CALL FlagError("Not implemented.",err,error,*999)
          ENDIF
        CASE(EQUATIONS_NONLINEAR)
          vectorMapping%createValuesCache%numberOfLinearMatrices=0
          vectorMapping%createValuesCache%numberOfResidualVariables=1
          IF(ASSOCIATED(dependentField%VARIABLE_TYPE_MAP(FIELD_DELUDELN_VARIABLE_TYPE)%ptr)) THEN
            vectorMapping%createValuesCache%rhsVariableType=dependentField%VARIABLE_TYPE_MAP(FIELD_DELUDELN_VARIABLE_TYPE)% &
              & ptr%VARIABLE_TYPE
          ELSE
            CALL FlagError("Not implemented.",err,error,*999)
          ENDIF
        CASE DEFAULT
          localError="The equations linearity type of "//TRIM(NumberToVString(equations%linearity,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_FIRST_ORDER_DYNAMIC,EQUATIONS_SECOND_ORDER_DYNAMIC)
        SELECT CASE(equations%linearity)
        CASE(EQUATIONS_LINEAR,EQUATIONS_NONLINEAR_BCS)
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
          !vectorMapping%createValuesCache%numberOfLinearMatrices=dependentField%NUMBER_OF_VARIABLES-2
          vectorMapping%createValuesCache%numberOfLinearMatrices=0
          IF(ASSOCIATED(dependentField%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%ptr)) THEN
            vectorMapping%createValuesCache%dynamicVariableType=dependentField% &
              & VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%ptr%VARIABLE_TYPE
          ELSE
            CALL FlagError("Not implemented.",err,error,*999)
          ENDIF
          IF(ASSOCIATED(dependentField%VARIABLE_TYPE_MAP(FIELD_DELUDELN_VARIABLE_TYPE)%ptr)) THEN
            vectorMapping%createValuesCache%rhsVariableType=dependentField% &
              & VARIABLE_TYPE_MAP(FIELD_DELUDELN_VARIABLE_TYPE)%ptr%VARIABLE_TYPE
          ELSE
            CALL FlagError("Not implemented.",err,error,*999)
          ENDIF
        CASE(EQUATIONS_NONLINEAR)
          ! SEBK 19/08/2009 not sure about mapping here
          !|
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
          vectorMapping%createValuesCache%numberOfLinearMatrices=0
          vectorMapping%createValuesCache%numberOfResidualVariables=1
          IF(ASSOCIATED(dependentField%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%ptr)) THEN
            vectorMapping%createValuesCache%dynamicVariableType=dependentField% &
              & VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%ptr%VARIABLE_TYPE
          ELSE
            CALL FlagError("Not implemented.",err,error,*999)
          ENDIF
          IF(ASSOCIATED(dependentField%VARIABLE_TYPE_MAP(FIELD_DELUDELN_VARIABLE_TYPE)%ptr)) THEN
            vectorMapping%createValuesCache%rhsVariableType=dependentField% &
              & VARIABLE_TYPE_MAP(FIELD_DELUDELN_VARIABLE_TYPE)%ptr%VARIABLE_TYPE
          ELSE
            CALL FlagError("Not implemented.",err,error,*999)
          ENDIF
        CASE DEFAULT
          localError="The equations linearity type of "//TRIM(NumberToVString(equations%linearity,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        !|
        ! SEBK 19/08/2009 not sure about mapping here
      CASE DEFAULT
        localError="The equations time dependence type of "//TRIM(NumberToVString(equations%timeDependence,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      localError="The number of dependent field variables of "// &
        & TRIM(NumberToVString(dependentField%NUMBER_OF_VARIABLES,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    !Allocate the dynamic matrix coefficients and set their values
    IF(vectorMapping%createValuesCache%numberOfDynamicMatrices>0) THEN
      ALLOCATE(vectorMapping%createValuesCache%dynamicMatrixCoefficients(vectorMapping% &
        & createValuesCache%numberOfDynamicMatrices),STAT=err)
      IF(err/=0) &
        & CALL FlagError("Could not allocate equations mapping create values cache dynamic matrix coefficients.",err,error,*999)
      vectorMapping%createValuesCache%dynamicMatrixCoefficients=1.0_DP !Equations matrices are added by default
    ENDIF
    !Allocate the residual variable types
    IF(vectorMapping%createValuesCache%numberOfResidualVariables>0) THEN
      ALLOCATE(vectorMapping%createValuesCache%residualVariableTypes(vectorMapping% &
        & createValuesCache%numberOfResidualVariables),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate equations mapping create values cache residual variable types.",err,error,*999)
      vectorMapping%createValuesCache%residualVariableTypes=0
      DO matrixIdx=1,vectorMapping%createValuesCache%numberOfResidualVariables
        variableNumber=1
        DO WHILE(vectorMapping%createValuesCache%residualVariableTypes(matrixIdx)==0.AND. &
          & variableNumber<=FIELD_NUMBER_OF_VARIABLE_TYPES)
          IF(ASSOCIATED(dependentField%VARIABLE_TYPE_MAP(variableNumber)%ptr)) THEN
            IF(dependentField%VARIABLE_TYPE_MAP(variableNumber)%ptr%VARIABLE_TYPE/= &
              & vectorMapping%createValuesCache%dynamicVariableType) THEN
              vectorMapping%createValuesCache%residualVariableTypes(matrixIdx)= &
                & dependentField%VARIABLE_TYPE_MAP(variableNumber)%ptr%VARIABLE_TYPE
            ENDIF
          ENDIF
          variableNumber=variableNumber+1
        ENDDO
      ENDDO !matrixIdx
      IF(vectorMapping%createValuesCache%residualVariableTypes(vectorMapping%createValuesCache%numberOfResidualVariables)==0) &
        & CALL FlagError("Invalid setup. All Jacobian matrices do not have a mapped dependent field variable.", &
        & err,error,*999)
    ENDIF
    !Allocate the linear matrix variable types and linear matrix coefficients and set their values
    IF(vectorMapping%createValuesCache%numberOfLinearMatrices>0) THEN
      ALLOCATE(vectorMapping%createValuesCache%linearMatrixVariableTypes(vectorMapping% &
        & createValuesCache%numberOfLinearMatrices),STAT=err)
      IF(err/=0) CALL  & 
        & FLAG_ERROR("Could not allocate equations mapping create values cache linear matrix variable types.",err,error,*999)
      ALLOCATE(vectorMapping%createValuesCache%linearMatrixCoefficients(vectorMapping% &
        & createValuesCache%numberOfLinearMatrices),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate equations mapping create values cache linear matrix coefficients.", &
        & err,error,*999)
      !Set up the matrices variable types
      vectorMapping%createValuesCache%linearMatrixVariableTypes=0
      variableNumber=1
      DO matrixIdx=1,vectorMapping%createValuesCache%numberOfLinearMatrices
        DO WHILE(vectorMapping%createValuesCache%linearMatrixVariableTypes(matrixIdx)==0.AND. &
          & variableNumber<=FIELD_NUMBER_OF_VARIABLE_TYPES)
          IF(ASSOCIATED(dependentField%VARIABLE_TYPE_MAP(variableNumber)%ptr)) THEN
            IF(dependentField%VARIABLE_TYPE_MAP(variableNumber)%ptr%VARIABLE_TYPE/= &
              & vectorMapping%createValuesCache%dynamicVariableType) THEN
              isResidualType=.FALSE.
              DO matrixIdx2=1,vectorMapping%createValuesCache%numberOfResidualVariables
                IF(dependentField%VARIABLE_TYPE_MAP(variableNumber)%ptr%VARIABLE_TYPE== &
                  & vectorMapping%createValuesCache%residualVariableTypes(matrixIdx2)) THEN
                  isResidualType=.TRUE.
                ENDIF
              ENDDO
              IF(.NOT.isResidualType) THEN
                vectorMapping%createValuesCache%linearMatrixVariableTypes(matrixIdx)= &
                  & dependentField%VARIABLE_TYPE_MAP(variableNumber)%ptr%VARIABLE_TYPE
              ENDIF
            ENDIF
          ENDIF
          variableNumber=variableNumber+1
        ENDDO
      ENDDO !matrixIdx
      IF(vectorMapping%createValuesCache%linearMatrixVariableTypes(vectorMapping% &
        & createValuesCache%numberOfLinearMatrices)==0) &
        & CALL FlagError("Invalid setup. All linear matrices do not have a mapped dependent field variable.", &
        & err,error,*999)
      vectorMapping%createValuesCache%linearMatrixCoefficients=1.0_DP !Equations matrices are added by default
    ENDIF
      
    EXITS("EquationsMapping_VectorCreateValuesCacheInitialise")
    RETURN
999 CALL EquationsMapping_VectorCreateValuesCacheFinalise(vectorMapping%createValuesCache,dummyErr,dummyError,*998)
998 ERRORS("EquationsMapping_VectorCreateValuesCacheInitialise",err,error)
    EXITS("EquationsMapping_VectorCreateValuesCacheInitialise")
    RETURN 1
    
  END SUBROUTINE EquationsMapping_VectorCreateValuesCacheInitialise

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

  !>Finalises the scalar equations mapping dot product mapping and deallocates all memory
  SUBROUTINE EquationsMapping_DotProductMappingFinalise(dotProductMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingDotProductType) :: dotProductMapping !<A pointer to the dot product mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMapping_DotProductMappingFinalise",err,error,*999)

    dotProductMapping%dotProductNumber=0
    NULLIFY(dotProductMapping%dotProductVariables(1)%ptr)
    NULLIFY(dotProductMapping%dotProductVariables(2)%ptr)    
       
    EXITS("EquationsMapping_DotProductMappingFinalise")
    RETURN
999 ERRORSEXITS("EquationsMapping_DotProductMappingFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMapping_DotProductMappingFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the scalar equations mapping dot product mapping
  SUBROUTINE EquationsMapping_DotProductMappingInitialise(dotProductMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingDotProductType) :: dotProductMapping !<The dot product mapping to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsMapping_DotProductMappingInitialise",err,error,*999)

    dotProductMapping%dotProductNumber=0
    NULLIFY(dotProductMapping%dotProductVariables(1)%ptr)
    NULLIFY(dotProductMapping%dotProductVariables(2)%ptr)
    dotProductMapping%dotProductCoefficient=1.0_DP
    
    EXITS("EquationsMapping_DotProductMappingInitialise")
    RETURN
999 ERRORS("EquationsMapping_DotProductMappingInitialise",err,error)
    EXITS("EquationsMapping_DotProductMappingInitialise")
    RETURN 1
    
  END SUBROUTINE EquationsMapping_DotProductMappingInitialise

  !
  !================================================================================================================================
  !

  !>Finalises the scalar equations mapping dot product mappings and deallocates all memory
  SUBROUTINE EquationsMapping_DotProductMappingsFinalise(dotProductMappings,err,error,*)

    !Argument variables
    TYPE(EquationsMappingDotProductsType), POINTER :: dotProductMappings !<A pointer to the dot product mappings to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dotProductIdx

    ENTERS("EquationsMapping_DotProductMappingsFinalise",err,error,*999)

    IF(ASSOCIATED(dotProductMappings)) THEN
      IF(ALLOCATED(dotProductMappings%dotProducts)) THEN
        DO dotProductIdx=1,SIZE(dotProductMappings%dotProducts,1)
          CALL EquationsMapping_DotProductMappingFinalise(dotProductMappings%dotProducts(dotProductIdx),err,error,*999)
        ENDDO !dotProductIdx
        DEALLOCATE(dotProductMappings%dotProducts)
      ENDIF
      DEALLOCATE(dotProductMappings)
    ENDIF
       
    EXITS("EquationsMapping_DotProductMappingsFinalise")
    RETURN
999 ERRORS("EquationsMapping_DotProductMappingsFinalise",err,error)
    EXITS("EquationsMapping_DotProductMappingsFinalise")
    RETURN 1
    
  END SUBROUTINE EquationsMapping_DotProductMappingsFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the scalar equations mapping dot product mappings
  SUBROUTINE EquationsMapping_DotProductMappingsInitialise(scalarMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingScalarType), POINTER :: scalarMapping !<A pointer to the scalar equations mapping to initialise the dot product mappings for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("EquationsMapping_DotProductMappingsInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(scalarMapping)) CALL FlagError("Scalar equations mapping is not associated.",err,error,*998)
    IF(ASSOCIATED(scalarMapping%dotProductMappings)) &
      & CALL FlagError("Scalar equations mapping dot product mappings is already associated.",err,error,*998)
    
    ALLOCATE(scalarMapping%dotProductMappings,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate scalar equations mapping dot product mappings.",err,error,*999)
    scalarMapping%dotProductMappings%scalarMapping=>scalarMapping
    scalarMapping%dotProductMappings%numberOfDotProducts=0
    
    EXITS("EquationsMapping_DotProductMappingsInitialise")
    RETURN
999 CALL EquationsMapping_DotProductMappingsFinalise(scalarMapping%dotProductMappings,dummyErr,dummyError,*998)
998 ERRORS("EquationsMapping_DotProductMappingsInitialise",err,error)
    EXITS("EquationsMapping_DotProductMappingsInitialise")
    RETURN 1
    
  END SUBROUTINE EquationsMapping_DotProductMappingsInitialise

  !
  !================================================================================================================================
  !

  !>Finalises the equations mapping dynamic mapping and deallocates all memory
  SUBROUTINE EquationsMapping_DynamicMappingFinalise(dynamicMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping !<A pointer to the dynamic mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx,variableType
 
    ENTERS("EquationsMapping_DynamicMappingFinalise",err,error,*999)

    IF(ASSOCIATED(dynamicMapping)) THEN
      IF(ALLOCATED(dynamicMapping%varToEquationsMatricesMaps)) THEN
        DO variableType=1,SIZE(dynamicMapping%varToEquationsMatricesMaps,1)
          CALL EquationsMapping_VarToEquatsMatricesMapFinalise(dynamicMapping% &
            & varToEquationsMatricesMaps(variableType),err,error,*999)
        ENDDO !variableType
        DEALLOCATE(dynamicMapping%varToEquationsMatricesMaps)        
      ENDIF
      IF(ALLOCATED(dynamicMapping%equationsMatrixToVarMaps)) THEN
        DO matrixIdx=1,SIZE(dynamicMapping%equationsMatrixToVarMaps,1)
          CALL EquationsMapping_EquationsMatrixToVarMapFinalise(dynamicMapping%equationsMatrixToVarMaps(matrixIdx), &
            & err,error,*999)
        ENDDO !matrixIdx
        DEALLOCATE(dynamicMapping%equationsMatrixToVarMaps)
      ENDIF
      IF(ALLOCATED(dynamicMapping%equationsRowToVariableDOFMaps)) &
        & DEALLOCATE(dynamicMapping%equationsRowToVariableDOFMaps)
      DEALLOCATE(dynamicMapping)
    ENDIF
       
    EXITS("EquationsMapping_DynamicMappingFinalise")
    RETURN
999 ERRORSEXITS("EquationsMapping_DynamicMappingFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMapping_DynamicMappingFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the equations mapping dynamic mapping
  SUBROUTINE EquationsMapping_DynamicMappingInitialise(vectorMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations mapping to initialise the dynamic mapping for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("EquationsMapping_DynamicMappingInitialise",err,error,*998)

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
    
    EXITS("EquationsMapping_DynamicMappingInitialise")
    RETURN
999 CALL EquationsMapping_DynamicMappingFinalise(vectorMapping%dynamicMapping,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMapping_DynamicMappingInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMapping_DynamicMappingInitialise

  !
  !================================================================================================================================
  !

  !>Sets/changes the matrices involved in dynamic equations mapping
  SUBROUTINE EquationsMapping_DynamicMatricesSetAll(vectorMapping,massMatrix,dampingMatrix,stiffnessMatrix,err,error,*)

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
    REAL(DP), ALLOCATABLE :: oldDynamicMatrixCoefficients(:)
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsMapping_DynamicMatricesSetAll",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector equations mapping is not associated.",err,error,*999)
    NULLIFY(vectorEquations)
    CALL EquationsMappingVector_VectorEquationsGet(vectorMapping,vectorEquations,err,error,*999)
    NULLIFY(equations)
    CALL EquationsVector_EquationsGet(vectorEquations,equations,err,error,*999)
    NULLIFY(equationsSet)
    CALL Equations_EquationsSetGet(equations,equationsSet,err,error,*999)
    NULLIFY(createValuesCache)
    CALL EquationsMappingVector_CreateValuesCacheGet(vectorMapping,createValuesCache,err,error,*999)
    
    SELECT CASE(equations%linearity)
    CASE(EQUATIONS_LINEAR,EQUATIONS_NONLINEAR_BCS)
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
      IF(numberOfDynamicMatrices>0) THEN
        ALLOCATE(oldDynamicMatrixCoefficients(createValuesCache%numberOfDynamicMatrices),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate old dynamic matrix coefficients.",err,error,*999)
        oldDynamicMatrixCoefficients(1:createValuesCache%numberOfDynamicMatrices)= &
          & createValuesCache%dynamicMatrixCoefficients(1:createValuesCache%numberOfDynamicMatrices)
        DEALLOCATE(createValuesCache%dynamicMatrixCoefficients)
        ALLOCATE(createValuesCache%dynamicMatrixCoefficients(numberOfDynamicMatrices),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate dynamic matrix coefficients.",err,error,*999)
        IF(newDynamicStiffnessMatrixNumber/=0) THEN
          IF(createValuesCache%dynamicStiffnessMatrixNumber==0) THEN
            createValuesCache%dynamicMatrixCoefficients(newDynamicStiffnessMatrixNumber)=1.0_DP
          ELSE
            createValuesCache%dynamicMatrixCoefficients(newDynamicStiffnessMatrixNumber)= &
              & oldDynamicMatrixCoefficients(createValuesCache%dynamicStiffnessMatrixNumber)
          ENDIF
        ENDIF
        IF(newDynamicDampingMatrixNumber/=0) THEN
          IF(createValuesCache%dynamicDampingMatrixNumber==0) THEN
            createValuesCache%dynamicMatrixCoefficients(newDynamicDampingMatrixNumber)=1.0_DP
          ELSE
            createValuesCache%dynamicMatrixCoefficients(newDynamicDampingMatrixNumber)= &
              & oldDynamicMatrixCoefficients(createValuesCache%dynamicDampingMatrixNumber)
          ENDIF
        ENDIF
        IF(newDynamicMassMatrixNumber/=0) THEN
          IF(createValuesCache%dynamicMassMatrixNumber==0) THEN
            createValuesCache%dynamicMatrixCoefficients(newDynamicMassMatrixNumber)=1.0_DP
          ELSE
            createValuesCache%dynamicMatrixCoefficients(newDynamicMassMatrixNumber)= &
              & oldDynamicMatrixCoefficients(createValuesCache%dynamicMassMatrixNumber)
          ENDIF
        ENDIF
        createValuesCache%numberOfDynamicMatrices=numberOfDynamicMatrices
        createValuesCache%dynamicStiffnessMatrixNumber=newDynamicStiffnessMatrixNumber
        createValuesCache%dynamicDampingMatrixNumber=newDynamicDampingMatrixNumber
        createValuesCache%dynamicMassMatrixNumber=newDynamicMassMatrixNumber
        IF(ALLOCATED(oldDynamicMatrixCoefficients)) DEALLOCATE(oldDynamicMatrixCoefficients)
      ELSE
        CALL FlagError("Invalid dynamic matrices set up. There are no dynamic equations matrices.",err,error,*999)
      ENDIF
    CASE(EQUATIONS_NONLINEAR)
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
      IF(numberOfDynamicMatrices>0) THEN
        ALLOCATE(oldDynamicMatrixCoefficients(createValuesCache%numberOfDynamicMatrices),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate old dynamic matrix coefficients.",err,error,*999)
        oldDynamicMatrixCoefficients(1:createValuesCache%numberOfDynamicMatrices)= &
          & createValuesCache%dynamicMatrixCoefficients(1:createValuesCache%numberOfDynamicMatrices)
        DEALLOCATE(createValuesCache%dynamicMatrixCoefficients)
        ALLOCATE(createValuesCache%dynamicMatrixCoefficients(numberOfDynamicMatrices),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate dynamic matrix coefficients.",err,error,*999)
        IF(newDynamicStiffnessMatrixNumber/=0) THEN
          IF(createValuesCache%dynamicStiffnessMatrixNumber==0) THEN
            createValuesCache%dynamicMatrixCoefficients(newDynamicStiffnessMatrixNumber)=1.0_DP
          ELSE
            createValuesCache%dynamicMatrixCoefficients(newDynamicStiffnessMatrixNumber)= &
              & oldDynamicMatrixCoefficients(createValuesCache%dynamicStiffnessMatrixNumber)
          ENDIF
        ENDIF
        IF(newDynamicDampingMatrixNumber/=0) THEN
          IF(createValuesCache%dynamicDampingMatrixNumber==0) THEN
            createValuesCache%dynamicMatrixCoefficients(newDynamicDampingMatrixNumber)=1.0_DP
          ELSE
            createValuesCache%dynamicMatrixCoefficients(newDynamicDampingMatrixNumber)= &
              & oldDynamicMatrixCoefficients(createValuesCache%dynamicDampingMatrixNumber)
          ENDIF
        ENDIF
        IF(newDynamicMassMatrixNumber/=0) THEN
          IF(createValuesCache%dynamicMassMatrixNumber==0) THEN
            createValuesCache%dynamicMatrixCoefficients(newDynamicMassMatrixNumber)=1.0_DP
          ELSE
            createValuesCache%dynamicMatrixCoefficients(newDynamicMassMatrixNumber)= &
              & oldDynamicMatrixCoefficients(createValuesCache%dynamicMassMatrixNumber)
          ENDIF
        ENDIF
        createValuesCache%numberOfDynamicMatrices=numberOfDynamicMatrices
        createValuesCache%dynamicStiffnessMatrixNumber=newDynamicStiffnessMatrixNumber
        createValuesCache%dynamicDampingMatrixNumber=newDynamicDampingMatrixNumber
        createValuesCache%dynamicMassMatrixNumber=newDynamicMassMatrixNumber
        IF(ALLOCATED(oldDynamicMatrixCoefficients)) DEALLOCATE(oldDynamicMatrixCoefficients)
      ELSE
        CALL FlagError("Invalid dynamic matrices set up. There are no dynamic equations matrices.",err,error,*999)
      ENDIF
    CASE DEFAULT
      localError="The equations linearity type of "//TRIM(NumberToVString(equations%linearity,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
   
    EXITS("EquationsMapping_DynamicMatricesSetAll")
    RETURN
999 IF(ALLOCATED(oldDynamicMatrixCoefficients)) DEALLOCATE(oldDynamicMatrixCoefficients)    
    ERRORSEXITS("EquationsMapping_DynamicMatricesSetAll",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMapping_DynamicMatricesSetAll

  !
  !================================================================================================================================
  !

  !>Sets/changes the matrices involved in a first order dynamic vector equations mapping
  SUBROUTINE EquationsMapping_DynamicMatricesSet1(vectorMapping,dampingMatrix,stiffnessMatrix,err,error,*)

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

    ENTERS("EquationsMapping_DynamicMatricesSet1",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector equations mapping is not associated.",err,error,*999)
    IF(vectorMapping%vectorMappingFinished) CALL FlagError("Vector equations mapping has already been finished.",err,error,*999)
    
    NULLIFY(vectorEquations)
    CALL EquationsMappingVector_VectorEquationsGet(vectorMapping,vectorEquations,err,error,*999)
    NULLIFY(equations)
    CALL EquationsVector_EquationsGet(vectorEquations,equations,err,error,*999)

    SELECT CASE(equations%timeDependence)
    CASE(EQUATIONS_STATIC)
      CALL FlagError("Can not set dynamic matrices for static equations.",err,error,*999)
    CASE(EQUATIONS_QUASISTATIC)
      CALL FlagError("Can not set dynamic matrices for quasi-static equations.",err,error,*999)
    CASE(EQUATIONS_FIRST_ORDER_DYNAMIC)
      IF(.NOT.dampingMatrix) CALL FlagWarning("No damping matrix for first order dynamic equations.",err,error,*999)
      CALL EquationsMapping_DynamicMatricesSetAll(vectorMapping,.FALSE.,dampingMatrix,stiffnessMatrix,err,error,*999)
    CASE(EQUATIONS_SECOND_ORDER_DYNAMIC)
      CALL FlagError("Need to specify three matrices to set for second order dynamic equations.",err,error,*999)
    CASE DEFAULT
      localError="The equations time dependence type of "//TRIM(NumberToVString(equations%timeDependence,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("EquationsMapping_DynamicMatricesSet1")
    RETURN
999 ERRORSEXITS("EquationsMapping_DynamicMatricesSet1",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMapping_DynamicMatricesSet1

  !
  !================================================================================================================================
  !

  !>Sets/changes the matrices involved in a second order dynamic equations mapping
  SUBROUTINE EquationsMapping_DynamicMatricesSet2(vectorMapping,massMatrix,dampingMatrix,stiffnessMatrix,err,error,*)

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

    ENTERS("EquationsMapping_DynamicMatricesSet2",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector equations mapping is not associated.",err,error,*999)
    IF(vectorMapping%vectorMappingFinished) CALL FlagError("Vector equations mapping has already been finished.",err,error,*999)

    NULLIFY(vectorEquations)
    CALL EquationsMappingVector_VectorEquationsGet(vectorMapping,vectorEquations,err,error,*999)
    NULLIFY(equations)
    CALL EquationsVector_EquationsGet(vectorEquations,equations,err,error,*999)

    SELECT CASE(equations%timeDependence)
    CASE(EQUATIONS_STATIC)
      CALL FlagError("Can not set dynamic matrices for static equations.",err,error,*999)
    CASE(EQUATIONS_QUASISTATIC)
      CALL FlagError("Can not set dynamic matrices for quasi-static equations.",err,error,*999)
    CASE(EQUATIONS_FIRST_ORDER_DYNAMIC)
      IF(massMatrix) THEN
        CALL FlagError("The mass matrix cannot be present for first order dynamic equations.",err,error,*999)
      ELSE
        IF(.NOT.dampingMatrix) CALL FlagWarning("No damping matrix for a first order dynamic system.",err,error,*999)
        CALL EquationsMapping_DynamicMatricesSetAll(vectorMapping,.FALSE.,dampingMatrix,stiffnessMatrix,err,error,*999)
      ENDIF
    CASE(EQUATIONS_SECOND_ORDER_DYNAMIC)
      IF(.NOT.massMatrix) CALL FlagWarning("No mass matrix for a second order dynamic system.",err,error,*999)
      CALL EquationsMapping_DynamicMatricesSetAll(vectorMapping,massMatrix,dampingMatrix,stiffnessMatrix,err,error,*999)
    CASE DEFAULT
      localError="The equations time dependence type of "//TRIM(NumberToVString(equations%timeDependence,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("EquationsMapping_DynamicMatricesSet2")
    RETURN
999 ERRORSEXITS("EquationsMapping_DynamicMatricesSet2",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMapping_DynamicMatricesSet2

  !
  !================================================================================================================================
  !

  !>Sets/changes the matrix coefficients in a first order dynamic equations mapping
  SUBROUTINE EquationsMapping_DynamicMatricesCoeffsSet1(vectorMapping,dampingMatrixCoefficient, &
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

    ENTERS("EquationsMapping_DynamicMatricesCoeffsSet1",err,error,*999)
    
    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector equations mapping is not associated.",err,error,*999)
    IF(vectorMapping%vectorMappingFinished) CALL FlagError("Vector equations mapping has already been finished.",err,error,*999)
 
    NULLIFY(createValuesCache)
    CALL EquationsMappingVector_CreateValuesCacheGet(vectorMapping,createValuesCache,err,error,*999)
    NULLIFY(vectorEquations)
    CALL EquationsMappingVector_VectorEquationsGet(vectorMapping,vectorEquations,err,error,*999)
    NULLIFY(equations)
    CALL EquationsVector_EquationsGet(vectorEquations,equations,err,error,*999)
    
    SELECT CASE(equations%timeDependence)
    CASE(EQUATIONS_STATIC)
      CALL FlagError("Can not set dynamic matrix coefficients for static equations.",err,error,*999)
    CASE(EQUATIONS_QUASISTATIC)
      CALL FlagError("Can not set dynamic matrix coefficients for quasi-static equations.",err,error,*999)
    CASE(EQUATIONS_FIRST_ORDER_DYNAMIC)
      IF(createValuesCache%dynamicStiffnessMatrixNumber/=0) THEN
        createValuesCache%dynamicMatrixCoefficients(createValuesCache%dynamicStiffnessMatrixNumber)=stiffnessMatrixCoefficient
      ENDIF
      IF(createValuesCache%dynamicDampingMatrixNumber/=0) THEN
        createValuesCache%dynamicMatrixCoefficients(createValuesCache%dynamicDampingMatrixNumber)=dampingMatrixCoefficient
      ENDIF
    CASE(EQUATIONS_SECOND_ORDER_DYNAMIC)
      CALL FlagError("Need to specify three matrix coefficients for second order dynamic equations.", &
        & err,error,*999)
    CASE DEFAULT
      localError="The equations time dependence type of "//TRIM(NumberToVString(equations%timeDependence,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("EquationsMapping_DynamicMatricesCoeffsSet1")
    RETURN
999 ERRORSEXITS("EquationsMapping_DynamicMatricesCoeffsSet1",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMapping_DynamicMatricesCoeffsSet1

  !
  !================================================================================================================================
  !

  !>Sets/changes the matrix coefficients in a second order dynamic equations mapping
  SUBROUTINE EquationsMapping_DynamicMatricesCoeffsSet2(vectorMapping,massMatrixCoefficient, &
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

    ENTERS("EquationsMapping_DynamicMatricesCoeffsSet2",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector equations mapping is not associated.",err,error,*999)
    IF(vectorMapping%vectorMappingFinished) CALL FlagError("Vector equations mapping has already been finished.",err,error,*999)

    NULLIFY(createValuesCache)
    CALL EquationsMappingVector_CreateValuesCacheGet(vectorMapping,createValuesCache,err,error,*999)
    NULLIFY(vectorEquations)
    CALL EquationsMappingVector_VectorEquationsGet(vectorMapping,vectorEquations,err,error,*999)
    NULLIFY(equations)
    CALL EquationsVector_EquationsGet(vectorEquations,equations,err,error,*999)
    
    SELECT CASE(equations%timeDependence)
    CASE(EQUATIONS_STATIC)
      CALL FlagError("Can not set dynamic matrices for static equations.",err,error,*999)
    CASE(EQUATIONS_QUASISTATIC)
      CALL FlagError("Can not set dynamic matrices for quasi-static equations.",err,error,*999)
    CASE(EQUATIONS_FIRST_ORDER_DYNAMIC)
      CALL FlagError("Need to specify two matrix coefficients for second order dynamic equations.",err,error,*999)
    CASE(EQUATIONS_SECOND_ORDER_DYNAMIC)
      IF(createValuesCache%dynamicStiffnessMatrixNumber/=0) &
        & createValuesCache%dynamicMatrixCoefficients(createValuesCache%dynamicStiffnessMatrixNumber)=stiffnessMatrixCoefficient
      IF(createValuesCache%dynamicDampingMatrixNumber/=0) &
        & createValuesCache%dynamicMatrixCoefficients(createValuesCache%dynamicDampingMatrixNumber)=dampingMatrixCoefficient
      IF(createValuesCache%dynamicMassMatrixNumber/=0) &
        & createValuesCache%dynamicMatrixCoefficients(createValuesCache%dynamicMassMatrixNumber)=massMatrixCoefficient
    CASE DEFAULT
      localError="The equations time dependence type of "//TRIM(NumberToVString(equations%timeDependence,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("EquationsMapping_DynamicMatricesCoeffsSet2")
    RETURN
999 ERRORSEXITS("EquationsMapping_DynamicMatricesCoeffsSet2",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMapping_DynamicMatricesCoeffsSet2

  !
  !================================================================================================================================
  !

  !>Sets the mapping between a dependent field variable and the equations set dynamic matrices
  SUBROUTINE EquationsMapping_DynamicVariableTypeSet(vectorMapping,dynamicVariableType,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the vector equations mapping to set
    INTEGER(INTG), INTENT(IN) :: dynamicVariableType !<The variable type associated with the vector equations dynamic matrices.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    TYPE(FIELD_TYPE), POINTER :: dependentField
    TYPE(VARYING_STRING) :: localError
    LOGICAL :: isResidualType

    ENTERS("EquationsMapping_DynamicVariableTypeSet",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector equations mapping is not associated",err,error,*999)
    IF(vectorMapping%vectorMappingFinished) CALL FlagError("Vector equations mapping has been finished.",err,error,*999)

    NULLIFY(createValuesCache)
    CALL EquationsMappingVector_CreateValuesCacheGet(vectorMapping,createValuesCache,err,error,*999)
    NULLIFY(vectorEquations)
    CALL EquationsMappingVector_VectorEquationsGet(vectorMapping,vectorEquations,err,error,*999)
    NULLIFY(equations)
    CALL EquationsVector_EquationsGet(vectorEquations,equations,err,error,*999)
    
    IF(dynamicVariableType==0) THEN
      createValuesCache%dynamicVariableType=0
    ELSE
      IF(equations%timeDependence==EQUATIONS_FIRST_ORDER_DYNAMIC.OR. &
        equations%timeDependence==EQUATIONS_SECOND_ORDER_DYNAMIC) THEN
        NULLIFY(equationsSet)
        CALL Equations_EquationsSetGet(equations,equationsSet,err,error,*999)      
        NULLIFY(dependentField)
        CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
        !Check the dynamic variable type is not being by other equations matrices or vectors
        IF(equations%linearity==EQUATIONS_NONLINEAR) THEN
          isResidualType=.FALSE.
          DO matrixIdx=1,createValuesCache%numberOfResidualVariables
            IF(createValuesCache%residualVariableTypes(matrixIdx)==dynamicVariableType) isResidualType=.TRUE.
          ENDDO !matrixIdx
          IF(.NOT.isResidualType) THEN
            localError="The specified dynamic variable type of "//TRIM(NumberToVString(dynamicVariableType,"*",err,error))// &
              & " is not the same as any residual variable type."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        END IF
        IF(createValuesCache%rhsVariableType==dynamicVariableType) THEN
          localError="The specified dynamic variable type of "//TRIM(NumberToVString(dynamicVariableType,"*",err,error))// &
            & " is the same as the variable type for the RHS vector."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        DO matrixIdx=1,createValuesCache%numberOfLinearMatrices
          IF(createValuesCache%linearMatrixVariableTypes(matrixIdx)==dynamicVariableType) THEN
            localError="The specified dynamic variable type of "//TRIM(NumberToVString(dynamicVariableType,"*",err,error))// &
              & " is the same as the variable type for linear matrix number "//TRIM(NumberToVString(matrixIdx,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ENDDO !matrixIdx
        !Check the dynamic variable type is defined on the dependent field
        IF(dynamicVariableType>=1.AND.dynamicVariableType<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          IF(ASSOCIATED(dependentField%VARIABLE_TYPE_MAP(dynamicVariableType)%ptr)) THEN
            vectorMapping%createValuesCache%dynamicVariableType=dynamicVariableType
          ELSE
            localError="The specified dynamic variable type of "//TRIM(NumberToVString(dynamicVariableType,"*",err,error))// &
              & " is not defined on the dependent field."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ELSE
          localError="The specified dynamic variable type of "//TRIM(NumberToVString(dynamicVariableType,"*",err,error))// &
            & " is invalid. The number must either be zero or >= 1 and <= "// &
            & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("The equations are not dynamic equations.",err,error,*999)
      ENDIF
    ENDIF
    
    EXITS("EquationsMapping_DynamicVariableTypeSet")
    RETURN
999 ERRORSEXITS("EquationsMapping_DynamicVariableTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMapping_DynamicVariableTypeSet
  
  !
  !================================================================================================================================
  !

  !>Finalises a variable to equations Jacobian map and deallocates all memory.
  SUBROUTINE EquationsMapping_EquatsJacobianToVarMapFinalise(equationsJacobianToVarMap,err,error,*)

    !Argument variables
    TYPE(EquationsJacobianToVarMapType) :: equationsJacobianToVarMap !<The equations Jacobian to variable map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsMapping_EquatsJacobianToVarMapFinalise",err,error,*999)
    
    IF(ALLOCATED(equationsJacobianToVarMap%equationsColumnToDOFVariableMap)) &
      & DEALLOCATE(equationsJacobianToVarMap%equationsColumnToDOFVariableMap)
    
    EXITS("EquationsMapping_EquatsJacobianToVarMapFinalise")
    RETURN
999 ERRORS("EquationsMapping_EquatsJacobianToVarMapFinalise",err,error)    
    EXITS("EquationsMapping_EquatsJacobianToVarMapFinalise")    
    RETURN 1
   
  END SUBROUTINE EquationsMapping_EquatsJacobianToVarMapFinalise

  !
  !================================================================================================================================
  !

  !>Initialises a variable to equations Jacobian map.
  SUBROUTINE EquationsMapping_EquatsJacobianToVarMapInitialise(equationsJacobianToVarMap,err,error,*)

    !Argument variables
    TYPE(EquationsJacobianToVarMapType) :: equationsJacobianToVarMap !<The equations Jacobian to variable map to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsMapping_EquatsJacobianToVarMapInitialise",err,error,*999)
    
    equationsJacobianToVarMap%variableType=0
    NULLIFY(equationsJacobianToVarMap%variable)
    NULLIFY(equationsJacobianToVarMap%jacobian)
    equationsJacobianToVarMap%numberOfColumns=0
    equationsJacobianToVarMap%jacobianCoefficient=0
    NULLIFY(equationsJacobianToVarMap%columnDOFSMapping)    
    
    EXITS("EquationsMapping_EquatsJacobianToVarMapInitialise")
    RETURN
999 ERRORS("EquationsMapping_EquatsJacobianToVarMapInitialise",err,error)    
    EXITS("EquationsMapping_EquatsJacobianToVarMapInitialise")    
    RETURN 1
   
  END SUBROUTINE EquationsMapping_EquatsJacobianToVarMapInitialise

  !
  !================================================================================================================================
  !

  !>Finalise an equations matrix to variable maps and deallocate all memory.
  SUBROUTINE EquationsMapping_EquationsMatrixToVarMapFinalise(equationsMatrixToVarMap,err,error,*)

    !Argument variables
    TYPE(EquationsMatrixToVarMapType) :: equationsMatrixToVarMap !<The equations matrix to variable map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsMapping_EquationsMatrixToVarMapFinalise",err,error,*999)

    IF(ALLOCATED(equationsMatrixToVarMap%columnToDOFMap)) &
      & DEALLOCATE(equationsMatrixToVarMap%columnToDOFMap)
    
    EXITS("EquationsMapping_EquationsMatrixToVarMapFinalise")
    RETURN
999 ERRORS("EquationsMapping_EquationsMatrixToVarMapFinalise",err,error)    
    EXITS("EquationsMapping_EquationsMatrixToVarMapFinalise")    
    RETURN 1
   
  END SUBROUTINE EquationsMapping_EquationsMatrixToVarMapFinalise

  !
  !================================================================================================================================
  !

  !>Initialise an equations matrix to variable maps.
  SUBROUTINE EquationsMapping_EquatsMatrixToVarMapInitialise(equationsMatrixToVarMap,err,error,*)

    !Argument variables
    TYPE(EquationsMatrixToVarMapType) :: equationsMatrixToVarMap !<The equations matrix to variable map to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsMapping_EquatsMatrixToVarMapInitialise",err,error,*999)

    equationsMatrixToVarMap%matrixNumber=0
    equationsMatrixToVarMap%variableType=0
    NULLIFY(equationsMatrixToVarMap%variable)
    equationsMatrixToVarMap%numberOfColumns=0
    equationsMatrixToVarMap%matrixCoefficient=1.0_DP !Matrices in an equation set are added by default
    NULLIFY(equationsMatrixToVarMap%columnDOFSMapping)
    
    EXITS("EquationsMapping_EquatsMatrixToVarMapInitialise")
    RETURN
999 ERRORS("EquationsMapping_EquatsMatrixToVarMapInitialise",err,error)    
    EXITS("EquationsMapping_EquatsMatrixToVarMapInitialise")    
    RETURN 1
   
  END SUBROUTINE EquationsMapping_EquatsMatrixToVarMapInitialise

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
       NULLIFY(vectorMapping%rowDofsMapping)
       CALL EquationsMapping_DynamicMappingFinalise(vectorMapping%dynamicMapping,err,error,*999)
       CALL EquationsMapping_LinearMappingFinalise(vectorMapping%linearMapping,err,error,*999)
       CALL EquationsMapping_NonlinearMappingFinalise(vectorMapping%nonlinearMapping,err,error,*999)
       CALL EquationsMapping_RHSMappingFinalise(vectorMapping%rhsMapping,err,error,*999)      
       CALL EquationsMapping_SourceMappingFinalise(vectorMapping%sourceMapping,err,error,*999)      
       CALL EquationsMapping_VectorCreateValuesCacheFinalise(vectorMapping%createValuesCache,err,error,*999)
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
    NULLIFY(vectorEquations%vectorMapping%rowDofsMapping)
    NULLIFY(vectorEquations%vectorMapping%dynamicMapping)
    NULLIFY(vectorEquations%vectorMapping%linearMapping)
    NULLIFY(vectorEquations%vectorMapping%nonlinearMapping)
    NULLIFY(vectorEquations%vectorMapping%rhsMapping)
    NULLIFY(vectorEquations%vectorMapping%sourceMapping)
    NULLIFY(vectorEquations%vectorMapping%createValuesCache)
    CALL EquationsMapping_LHSMappingInitialise(vectorEquations%vectorMapping,err,error,*999)
    
    CALL EquationsMapping_VectorCreateValuesCacheInitialise(vectorEquations%vectorMapping,err,error,*999)        
       
    EXITS("EquationsMapping_VectorInitialise")
    RETURN
999 CALL EquationsMapping_VectorFinalise(vectorEquations%vectorMapping,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMapping_VectorInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMapping_VectorInitialise

  !
  !================================================================================================================================
  !

  !>Finalises the scalar equations mapping function mapping and deallocates all memory
  SUBROUTINE EquationsMapping_FunctionMappingFinalise(functionMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingFunctionType) :: functionMapping !<A pointer to the function mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMapping_FunctionMappingFinalise",err,error,*999)

    functionMapping%functionNumber=0
       
    EXITS("EquationsMapping_FunctionMappingFinalise")
    RETURN
999 ERRORSEXITS("EquationsMapping_FunctionMappingFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMapping_FunctionMappingFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the scalar equations mapping function mapping
  SUBROUTINE EquationsMapping_FunctionMappingInitialise(functionMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingFunctionType) :: functionMapping !<The function mapping to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsMapping_FunctionMappingInitialise",err,error,*999)

    functionMapping%functionNumber=0
    
    EXITS("EquationsMapping_FunctionMappingInitialise")
    RETURN
999 ERRORS("EquationsMapping_FunctionMappingInitialise",err,error)
    EXITS("EquationsMapping_FunctionMappingInitialise")
    RETURN 1
    
  END SUBROUTINE EquationsMapping_FunctionMappingInitialise

  !
  !================================================================================================================================
  !

  !>Finalises the scalar equations mapping function mappings and deallocates all memory
  SUBROUTINE EquationsMapping_FunctionMappingsFinalise(functionMappings,err,error,*)

    !Argument variables
    TYPE(EquationsMappingFunctionsType), POINTER :: functionMappings !<A pointer to the function mappings to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: functionIdx

    ENTERS("EquationsMapping_FunctionMappingsFinalise",err,error,*999)

    IF(ASSOCIATED(functionMappings)) THEN
      IF(ALLOCATED(functionMappings%functions)) THEN
        DO functionIdx=1,SIZE(functionMappings%functions,1)
          CALL EquationsMapping_FunctionMappingFinalise(functionMappings%functions(functionIdx),err,error,*999)
        ENDDO !functionIdx
        DEALLOCATE(functionMappings%functions)
      ENDIF
      DEALLOCATE(functionMappings)
    ENDIF
       
    EXITS("EquationsMapping_FunctionMappingsFinalise")
    RETURN
999 ERRORS("EquationsMapping_FunctionMappingsFinalise",err,error)
    EXITS("EquationsMapping_FunctionMappingsFinalise")
    RETURN 1
    
  END SUBROUTINE EquationsMapping_FunctionMappingsFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the scalar equations mapping function mappings
  SUBROUTINE EquationsMapping_FunctionsMappingsInitialise(scalarMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingScalarType), POINTER :: scalarMapping !<A pointer to the scalar equations mapping to initialise the function mappings for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("EquationsMapping_FunctionMappingsInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(scalarMapping)) CALL FlagError("Scalar equations mapping is not associated.",err,error,*998)
    IF(ASSOCIATED(scalarMapping%functionMappings)) &
      & CALL FlagError("Scalar equations mapping function mappings is already associated.",err,error,*998)
    
    ALLOCATE(scalarMapping%functionMappings,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate scalar equations mapping function mappings.",err,error,*999)
    scalarMapping%functionMappings%scalarMapping=>scalarMapping
    scalarMapping%functionMappings%numberOfFunctions=0
    
    EXITS("EquationsMapping_FunctionMappingsInitialise")
    RETURN
999 CALL EquationsMapping_FunctionMappingsFinalise(scalarMapping%functionMappings,dummyErr,dummyError,*998)
998 ERRORS("EquationsMapping_FunctionMappingsInitialise",err,error)
    EXITS("EquationsMapping_FunctionMappingsInitialise")
    RETURN 1
    
  END SUBROUTINE EquationsMapping_FunctionsMappingsInitialise

  !
  !================================================================================================================================
  !

  !>Sets the mapping between a dependent field variable and the equations set residual vector.
  SUBROUTINE EquationsMapping_ResidualVariablesNumberSet(vectorMapping,numberOfVariables,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations mapping to set
    INTEGER(INTG), INTENT(IN) :: numberOfVariables !<The number of residual variables for this equations set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: previousNumber,minNumber
    INTEGER(INTG), ALLOCATABLE :: newResidualVariableTypes(:)
    TYPE(EquationsMappingVectorCreateValuesCacheType), POINTER :: createValuesCache

    ENTERS("EquationsMapping_ResidualVariablesNumberSet",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector equations mapping is not associated",err,error,*999)
    IF(vectorMapping%vectorMappingFinished) CALL FlagError("Vector equations mapping has been finished.",err,error,*999)

    NULLIFY(createValuesCache)
    CALL EquationsMappingVector_CreateValuesCacheGet(vectorMapping,createValuesCache,err,error,*999)
    previousNumber=createValuesCache%numberOfResidualVariables
    IF(numberOfVariables/=previousNumber) THEN
      !Create new residual_variable_types array and copy over previous values
      minNumber=MIN(numberOfVariables,previousNumber)
      ALLOCATE(newResidualVariableTypes(numberOfVariables),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate new residual variable types.",err,error,*999)
      newResidualVariableTypes=0
      newResidualVariableTypes(1:minNumber)=createValuesCache%residualVariableTypes(1:minNumber)
      CALL MOVE_ALLOC(newResidualVariableTypes,createValuesCache%residualVariableTypes)
      !Set number of residual variables
      createValuesCache%numberOfResidualVariables=numberOfVariables
    ENDIF

    EXITS("EquationsMapping_ResidualVariablesNumberSet")
    RETURN
999 ERRORSEXITS("EquationsMapping_ResidualVariablesNumberSet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMapping_ResidualVariablesNumberSet

  !
  !================================================================================================================================
  !

  !>Finalises the vector equations mapping LHS mapping and deallocates all memory
  SUBROUTINE EquationsMapping_LHSMappingFinalise(lhsMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping !<A pointer to the LHS mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMapping_LHSMappingFinalise",err,error,*999)

    IF(ASSOCIATED(lhsMapping)) THEN
      IF(ALLOCATED(lhsMapping%lhsDOFToEquationsRowMap)) DEALLOCATE(lhsMapping%lhsDOFToEquationsRowMap)
      IF(ALLOCATED(lhsMapping%equationsRowToLHSDOFMap)) DEALLOCATE(lhsMapping%equationsRowToLHSDOFMap)
      DEALLOCATE(lhsMapping)
    ENDIF
       
    EXITS("EquationsMapping_LHSMappingFinalise")
    RETURN
999 ERRORSEXITS("EquationsMapping_LHSMappingFinalise",err,error)
    RETURN 1
  END SUBROUTINE EquationsMapping_LHSMappingFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the equations mapping LHS mapping
  SUBROUTINE EquationsMapping_LHSMappingInitialise(vectorMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations mapping to initialise the LHS mapping for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("EquationsMapping_LHSMappingInitialise",err,error,*998)

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
       
    EXITS("EquationsMapping_LHSMappingInitialise")
    RETURN
999 CALL EquationsMapping_LHSMappingFinalise(vectorMapping%lhsMapping,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMapping_LHSMappingInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMapping_LHSMappingInitialise

  !
  !================================================================================================================================
  !

  !>Finalises the equations mapping linear mapping and deallocates all memory
  SUBROUTINE EquationsMapping_LinearMappingFinalise(linearMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping !<A pointer to the linear mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx,variableTypeIdx
 
    ENTERS("EquationsMapping_LinearMappingFinalise",err,error,*999)

    IF(ASSOCIATED(linearMapping)) THEN
      IF(ALLOCATED(linearMapping%linearMatrixVariableTypes)) DEALLOCATE(linearMapping%linearMatrixVariableTypes)
      IF(ALLOCATED(linearMapping%varToEquationsMatricesMaps)) THEN
        DO variableTypeIdx=1,SIZE(linearMapping%varToEquationsMatricesMaps,1)
          CALL EquationsMapping_VarToEquatsMatricesMapFinalise(linearMapping%varToEquationsMatricesMaps(variableTypeIdx), &
            & err,error,*999)
        ENDDO !variableTypeIdx
        DEALLOCATE(linearMapping%varToEquationsMatricesMaps)        
      ENDIF
      IF(ALLOCATED(linearMapping%equationsMatrixToVarMaps)) THEN
        DO matrixIdx=1,SIZE(linearMapping%equationsMatrixToVarMaps,1)
          CALL EquationsMapping_EquationsMatrixToVarMapFinalise(linearMapping% equationsMatrixToVarMaps(matrixIdx), &
            & err,error,*999)
        ENDDO !matrixIdx
        DEALLOCATE(linearMapping%equationsMatrixToVarMaps)
      ENDIF
      IF(ALLOCATED(linearMapping%equationsRowToVariableDOFMaps)) DEALLOCATE(linearMapping%equationsRowToVariableDOFMaps)
      DEALLOCATE(linearMapping)
    ENDIF
       
    EXITS("EquationsMapping_LinearMappingFinalise")
    RETURN
999 ERRORSEXITS("EquationsMapping_LinearMappingFinalise",err,error)
    RETURN 1
  END SUBROUTINE EquationsMapping_LinearMappingFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the equations mapping linear mapping
  SUBROUTINE EquationsMapping_LinearMappingInitialise(vectorMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations mapping to initialise the linear mapping for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("EquationsMapping_LinearMappingInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector equations mapping is not associated.",err,error,*998)
    IF(ASSOCIATED(vectorMapping%linearMapping)) &
      & CALL FlagError("Vector equations mapping linear mapping is already associated.",err,error,*998)
    
    ALLOCATE(vectorMapping%linearMapping,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate equations mapping linear mapping.",err,error,*999)
    vectorMapping%linearMapping%vectorMapping=>vectorMapping       
    vectorMapping%linearMapping%numberOfLinearMatrices=0
    vectorMapping%linearMapping%numberOfLinearMatrixVariables=0
       
    EXITS("EquationsMapping_LinearMappingInitialise")
    RETURN
999 CALL EquationsMapping_LinearMappingFinalise(vectorMapping%linearMapping,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMapping_LinearMappingInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMapping_LinearMappingInitialise

  !
  !================================================================================================================================
  !

  !>Sets the coefficients for the linear equations matrices in an equation set. 
  SUBROUTINE EquationsMapping_LinearMatricesCoeffsSet(vectorMapping,linearMatrixCoefficients,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations mapping.
    REAL(DP), INTENT(IN) :: linearMatrixCoefficients(:) !<The linear matrix coefficients
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(EquationsMappingVectorCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsMapping_LinearMatricesCoeffsSet",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector equations mapping is not associated.",err,error,*999)
    IF(vectorMapping%vectorMappingFinished) CALL FlagError("Vector equations mapping has been finished.",err,error,*999)
    
    NULLIFY(createValuesCache)
    CALL EquationsMappingVector_CreateValuesCacheGet(vectorMapping,createValuesCache,err,error,*999)

    IF(SIZE(linearMatrixCoefficients,1)==createValuesCache%numberOfLinearMatrices) THEN
      createValuesCache%linearMatrixCoefficients(1:createValuesCache%numberOfLinearMatrices)= &
        & linearMatrixCoefficients(1:createValuesCache%numberOfLinearMatrices)
    ELSE
      localError="Invalid size of linear matrix coefficeints. The size of the supplied array ("// &
        & TRIM(NumberToVString(SIZE(linearMatrixCoefficients,1),"*",err,error))// &
        & ") must match the number of linear equations matrices ("// &
        & TRIM(NumberToVString(vectorMapping%createValuesCache%numberOfLinearMatrices,"*",err,error))//")."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("EquationsMapping_LinearMatricesCoeffsSet")
    RETURN
999 ERRORSEXITS("EquationsMapping_LinearMatricesCoeffsSet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMapping_LinearMatricesCoeffsSet

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of linear equations matrices
  SUBROUTINE EquationsMapping_LinearMatricesNumberSet(vectorMapping,numberOfLinearMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations mapping to set the number of matrices for.
    INTEGER(INTG), INTENT(IN) :: numberOfLinearMatrices !<The number of linear equations matrices for the mapping.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    INTEGER(INTG), ALLOCATABLE :: oldLinearMatrixVariableTypes(:)
    REAL(DP), ALLOCATABLE :: oldLinearMatrixCoefficients(:)
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    TYPE(FIELD_TYPE), POINTER :: dependentField
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsMapping_LinearMatricesNumberSet",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector equations mapping is not associated",err,error,*999)
    IF(vectorMapping%vectorMappingFinished) CALL FlagError("Vector equations mapping has been finished",err,error,*999)
      

    NULLIFY(createValuesCache)
    CALL EquationsMappingVector_CreateValuesCacheGet(vectorMapping,createValuesCache,err,error,*999)
    NULLIFY(vectorEquations)
    CALL EquationsMappingVector_VectorEquationsGet(vectorMapping,vectorEquations,err,error,*999)
    NULLIFY(equations)
    CALL EquationsVector_EquationsGet(vectorEquations,equations,err,error,*999)
    NULLIFY(equationsSet)
    CALL Equations_EquationsSetGet(equations,equationsSet,err,error,*999)

    !Check number of matrices to create is valid
    SELECT CASE(equations%timeDependence)
    CASE(EQUATIONS_STATIC,EQUATIONS_QUASISTATIC)
      SELECT CASE(equations%linearity)
      CASE(EQUATIONS_LINEAR,EQUATIONS_NONLINEAR_BCS)
        IF(createValuesCache%rhsVariableType==0) THEN                  
          IF(numberOfLinearMatrices<1) THEN
            localError="The specified number of linear matrices of "//TRIM(NumberToVString(numberOfLinearMatrices,"*", &
              & err,error))//" is invalid. For non-dynamic linear problems without a equations set RHS the number must be >= 1."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ELSE
          IF(numberOfLinearMatrices<1) THEN
            localError="The specified number of linear matrices of "//TRIM(NumberToVString(numberOfLinearMatrices,"*", &
              & err,error))//" is invalid. For non-dynamic linear problems with a equations set RHS the number must be >= 1."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ENDIF
      CASE(EQUATIONS_NONLINEAR)
        IF(numberOfLinearMatrices<0.OR.numberOfLinearMatrices>FIELD_NUMBER_OF_VARIABLE_TYPES-2) THEN
          localError="The specified number of linear matrices of "//TRIM(NumberToVString(numberOfLinearMatrices,"*", &
            & err,error))//" is invalid. For non-dynamic non-linear problems the number must be between >= 0 and <= "// &
            & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES-2,"*",err,error))
          CALL FlagError(localError,err,error,*999)
        ENDIF
      CASE DEFAULT
        localError="The equations linearity type of "//TRIM(NumberToVString(equations%linearity,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_FIRST_ORDER_DYNAMIC,EQUATIONS_SECOND_ORDER_DYNAMIC)
      SELECT CASE(equations%linearity)
      CASE(EQUATIONS_LINEAR,EQUATIONS_NONLINEAR_BCS)
        IF(createValuesCache%rhsVariableType==0) THEN                  
          IF(numberOfLinearMatrices<1.OR.numberOfLinearMatrices>FIELD_NUMBER_OF_VARIABLE_TYPES-1) THEN
            localError="The specified number of linear matrices of "//TRIM(NumberToVString(numberOfLinearMatrices,"*", &
              & err,error))//" is invalid. For dynamic linear problems without a equations set RHS the number must be "// &
              & "between >= 1 and <= "//TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES-1,"*",err,error))
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ELSE
          IF(numberOfLinearMatrices<0) THEN
            localError="The specified number of linear matrices of "//TRIM(NumberToVString(numberOfLinearMatrices,"*", &
              & err,error))//" is invalid. For dynamic linear problems with a equations set RHS the number must be >= 0."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ENDIF
      CASE(EQUATIONS_NONLINEAR)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The equations linearity type of "//TRIM(NumberToVString(equations%linearity,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The equations time dependence type of "//TRIM(NumberToVString(equations%timeDependence,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    !If we need to reallocate and reset all the create_values cache arrays and change the number of matrices
    IF(numberOfLinearMatrices/=createValuesCache%numberOfLinearMatrices) THEN
      IF(createValuesCache%numberOfLinearMatrices>0) THEN                  
        ALLOCATE(oldLinearMatrixVariableTypes(createValuesCache%numberOfLinearMatrices),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate old linear matrix variable types.",err,error,*999)
        ALLOCATE(oldLinearMatrixCoefficients(createValuesCache%numberOfLinearMatrices),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate old linear matrix coefficients.",err,error,*999)
        oldLinearMatrixVariableTypes(1:createValuesCache%numberOfLinearMatrices)= &
          & createValuesCache%linearMatrixVariableTypes(1:createValuesCache%numberOfLinearMatrices)
        oldLinearMatrixCoefficients(1:createValuesCache%numberOfLinearMatrices)= &
          & createValuesCache%linearMatrixCoefficients(1:createValuesCache%numberOfLinearMatrices)
      ENDIF
      IF(ALLOCATED(createValuesCache%linearMatrixVariableTypes)) DEALLOCATE(createValuesCache%linearMatrixVariableTypes)
      IF(ALLOCATED(createValuesCache%linearMatrixCoefficients)) DEALLOCATE(createValuesCache%linearMatrixCoefficients)
      ALLOCATE(createValuesCache%linearMatrixVariableTypes(numberOfLinearMatrices),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate linear matrix variable types.",err,error,*999)
      ALLOCATE(createValuesCache%linearMatrixCoefficients(numberOfLinearMatrices),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate linear matrix coefficients.",err,error,*999)
      IF(createValuesCache%numberOfLinearMatrices>0) THEN                  
        IF(numberOfLinearMatrices>createValuesCache%numberOfLinearMatrices) THEN
          createValuesCache%linearMatrixVariableTypes(1:createValuesCache%numberOfLinearMatrices)= &
            & oldLinearMatrixVariableTypes
          createValuesCache%linearMatrixVariableTypes(createValuesCache%numberOfLinearMatrices+1: &
            & numberOfLinearMatrices)=oldLinearMatrixVariableTypes(1)
          createValuesCache%linearMatrixCoefficients(1:createValuesCache%numberOfLinearMatrices)= &
            & oldLinearMatrixCoefficients
          createValuesCache%linearMatrixCoefficients(createValuesCache%numberOfLinearMatrices+1: &
            & numberOfLinearMatrices)=oldLinearMatrixCoefficients(1)
        ELSE
          createValuesCache%linearMatrixVariableTypes(1:numberOfLinearMatrices)= &
            & oldLinearMatrixVariableTypes(1:numberOfLinearMatrices)
          createValuesCache%linearMatrixCoefficients(1:numberOfLinearMatrices)= &
            & oldLinearMatrixCoefficients(1:numberOfLinearMatrices)
        ENDIF
      ELSE
        NULLIFY(dependentField)
        CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
        createValuesCache%linearMatrixVariableTypes=0
        SELECT CASE(equations%timeDependence)
        CASE(EQUATIONS_STATIC,EQUATIONS_QUASISTATIC)
          SELECT CASE(equations%linearity)
          CASE(EQUATIONS_LINEAR,EQUATIONS_NONLINEAR_BCS)
            IF(ASSOCIATED(dependentField%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%ptr)) THEN
              createValuesCache%linearMatrixVariableTypes(1)=dependentField%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%ptr% &
                & VARIABLE_TYPE
            ELSE
              CALL FlagError("Not implemented.",err,error,*999)
            ENDIF
            DO matrixIdx=2,createValuesCache%numberOfLinearMatrices
              IF(ASSOCIATED(dependentField%VARIABLE_TYPE_MAP(matrixIdx+1)%ptr)) THEN
                createValuesCache%linearMatrixVariableTypes(matrixIdx)=dependentField%VARIABLE_TYPE_MAP(matrixIdx+1)%ptr% &
                  & VARIABLE_TYPE
              ELSE
                CALL FlagError("Not implemented.",err,error,*999)
              ENDIF
            ENDDO !matrixIdx
          CASE(EQUATIONS_NONLINEAR)
            DO matrixIdx=1,createValuesCache%numberOfLinearMatrices
              IF(ASSOCIATED(dependentField%VARIABLE_TYPE_MAP(matrixIdx+2)%ptr)) THEN
                createValuesCache%linearMatrixVariableTypes(matrixIdx)=dependentField%VARIABLE_TYPE_MAP(matrixIdx+2)%ptr% &
                  & VARIABLE_TYPE
              ELSE
                CALL FlagError("Not implemented.",err,error,*999)
              ENDIF
            ENDDO !matrixIdx
          CASE DEFAULT
            localError="The equations linearity type of "//TRIM(NumberToVString(equations%linearity,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(EQUATIONS_FIRST_ORDER_DYNAMIC,EQUATIONS_SECOND_ORDER_DYNAMIC)
          SELECT CASE(equations%linearity)
          CASE(EQUATIONS_LINEAR,EQUATIONS_NONLINEAR_BCS)
            DO matrixIdx=1,createValuesCache%numberOfLinearMatrices
              IF(ASSOCIATED(dependentField%VARIABLE_TYPE_MAP(matrixIdx+2)%ptr)) THEN
                createValuesCache%linearMatrixVariableTypes(matrixIdx)=dependentField%VARIABLE_TYPE_MAP(matrixIdx+2)%ptr% &
                  & VARIABLE_TYPE
              ELSE
                CALL FlagError("Not implemented.",err,error,*999)
              ENDIF
            ENDDO !matrixIdx
          CASE(EQUATIONS_NONLINEAR)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The equations linearity type of "//TRIM(NumberToVString(equations%linearity,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The equations time dependence type of "//TRIM(NumberToVString(equations%timeDependence,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        createValuesCache%linearMatrixCoefficients=1.0_DP !Equations matrices are added by default
      ENDIF
      createValuesCache%numberOfLinearMatrices=numberOfLinearMatrices
      IF(ALLOCATED(oldLinearMatrixVariableTypes)) DEALLOCATE(oldLinearMatrixVariableTypes)
      IF(ALLOCATED(oldLinearMatrixCoefficients)) DEALLOCATE(oldLinearMatrixCoefficients)
    ENDIF
    
    EXITS("EquationsMapping_LinearMatricesNumberSet")
    RETURN
999 IF(ALLOCATED(oldLinearMatrixVariableTypes)) DEALLOCATE(oldLinearMatrixVariableTypes)    
    IF(ALLOCATED(oldLinearMatrixCoefficients)) DEALLOCATE(oldLinearMatrixCoefficients)    
    ERRORSEXITS("EquationsMapping_LinearMatricesNumberSet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMapping_LinearMatricesNumberSet

  !
  !================================================================================================================================
  !

  !>Sets the mapping between the dependent field variable types and the linear equations matrices
  SUBROUTINE EquationsMapping_LinearMatricesVariableTypesSet(vectorMapping,linearMatrixVariableTypes,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations mapping
    INTEGER(INTG), INTENT(IN) :: linearMatrixVariableTypes(:) !<The matrix variable types to map to each linear equations matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    TYPE(FIELD_TYPE), POINTER :: dependentField
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("EquationsMapping_LinearMatricesVariableTypesSet",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector equations mapping is not associated.",err,error,*999)
    IF(vectorMapping%vectorMappingFinished) CALL FlagError("Vector equations mapping has been finished.",err,error,*999)

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
    IF(SIZE(linearMatrixVariableTypes,1)/=createValuesCache%numberOfLinearMatrices) THEN
      localError="Invalid size of linear matrix variable types. The size of the supplied array ("// &
        & TRIM(NumberToVString(SIZE(linearMatrixVariableTypes,1),"*",err,error))// &
        & ") must match the number of linear equations matrices ("// &
        & TRIM(NumberToVString(vectorMapping%createValuesCache%numberOfLinearMatrices,"*",err,error))//")."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    DO matrixIdx=1,createValuesCache%numberOfLinearMatrices
      IF(linearMatrixVariableTypes(matrixIdx)/=0) THEN
        !Check the residual variable type is not being by other equations matrices or vectors
        !Don't check against the residual variable as we can have linear parts of nonlinear equations
        IF(createValuesCache%dynamicVariableType==linearMatrixVariableTypes(matrixIdx)) THEN
          localError="The specified linear matrix variable type of "// &
            & TRIM(NumberToVString(linearMatrixVariableTypes(matrixIdx),"*",err,error))// &
            & " for linear matrix number "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
            & " is the same as the variable type for the dynamic matrices."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        IF(createValuesCache%rhsVariableType==linearMatrixVariableTypes(matrixIdx)) THEN
          localError="The specified linear matrix variable type of "// &
            & TRIM(NumberToVString(linearMatrixVariableTypes(matrixIdx),"*",err,error))// &
            & " for linear matrix number "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
            & " is the same as the variable type for the RHS vector."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        !Check to see if the linear matrix variable numbers are defined on the dependent field
        IF(linearMatrixVariableTypes(matrixIdx)>=1.OR. &
          & linearMatrixVariableTypes(matrixIdx)<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          IF(.NOT.ASSOCIATED(dependentField%VARIABLE_TYPE_MAP(linearMatrixVariableTypes(matrixIdx))%ptr)) THEN
            localError="The linear matrix variable type of "// &
              & TRIM(NumberToVString(linearMatrixVariableTypes(matrixIdx),"*",err,error))// &
              & " for linear matrix NUMBER "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
              & " is not defined on the dependent field."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ELSE
          localError="The linear matrix variable type of "// &
            & TRIM(NumberToVString(linearMatrixVariableTypes(matrixIdx),"*",err,error))// &
            & " for linear matrix number "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
            & " is invalid. The variable types must be either zero or >= 1 and <= "// &
            & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDIF
    ENDDO !matrixIdx
    
    createValuesCache%linearMatrixVariableTypes(1:SIZE(linearMatrixVariableTypes))= &
      & linearMatrixVariableTypes(1:SIZE(linearMatrixVariableTypes))
       
    EXITS("EquationsMapping_LinearMatricesVariableTypesSet")
    RETURN
999 ERRORS("EquationsMapping_LinearMatricesVariableTypesSet",err,error)
    EXITS("EquationsMapping_LinearMatricesVariableTypesSet")
    RETURN 1
    
  END SUBROUTINE EquationsMapping_LinearMatricesVariableTypesSet

  !
  !================================================================================================================================
  !

  !>Finalises the equations mapping nonlinear mapping and deallocates all memory
  SUBROUTINE EquationsMapping_NonlinearMappingFinalise(nonlinearMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping !<A pointer to the nonlinear mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) matrixIdx
 
    ENTERS("EquationsMapping_NonlinearMappingFinalise",err,error,*999)

    IF(ASSOCIATED(nonlinearMapping)) THEN
      DO matrixIdx=1,nonlinearMapping%numberOfResidualVariables
        CALL EquationsMapping_VarToEquatsJacobianMapFinalise(nonlinearMapping%varToJacobianMap(matrixIdx),err,error,*999)
        CALL EquationsMapping_EquatsJacobianToVarMapFinalise(nonlinearMapping%jacobianToVarMap(matrixIdx),err,error,*999)
      ENDDO !matrixIdx
      IF(ALLOCATED(nonlinearMapping%equationsRowToResidualDOFMap)) &
        & DEALLOCATE(nonlinearMapping%equationsRowToResidualDOFMap)
      IF(ALLOCATED(nonlinearMapping%residualVariables)) DEALLOCATE(nonlinearMapping%residualVariables)
      IF(ALLOCATED(nonlinearMapping%varToJacobianMap)) DEALLOCATE(nonlinearMapping%varToJacobianMap)
      IF(ALLOCATED(nonlinearMapping%jacobianToVarMap)) DEALLOCATE(nonlinearMapping%jacobianToVarMap)
      DEALLOCATE(nonlinearMapping)
    ENDIF
    
    EXITS("EquationsMapping_NonlinearMappingFinalise")
    RETURN
999 ERRORSEXITS("EquationsMapping_NonlinearMappingFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMapping_NonlinearMappingFinalise
  
  !
  !================================================================================================================================
  !

  !>Initialises the equations mapping nonlinear mapping
  SUBROUTINE EquationsMapping_NonlinearMappingInitialise(vectorMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations mapping to initialise the nonlinear mapping for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("EquationsMapping_NonlinearMappingInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector equations mapping is not associated.",err,error,*998)
    IF(ASSOCIATED(vectorMapping%nonlinearMapping)) &
      & CALL FlagError("Equations mapping nonlinear mapping is already associated.",err,error,*998)
     
    ALLOCATE(vectorMapping%nonlinearMapping,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate equations mapping nonlinear mapping.",err,error,*999)
    vectorMapping%nonlinearMapping%vectorMapping=>vectorMapping
    vectorMapping%nonlinearMapping%numberOfResiduals=1
    vectorMapping%nonlinearMapping%numberOfResidualVariables=0
    vectorMapping%nonlinearMapping%residualCoefficient=1.0_DP

    EXITS("EquationsMapping_NonlinearMappingInitialise")
    RETURN
999 CALL EquationsMapping_NonlinearMappingFinalise(vectorMapping%nonlinearMapping,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMapping_NonlinearMappingInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMapping_NonlinearMappingInitialise

  !
  !================================================================================================================================
  !

  !>Finalises the scalar equations mapping norm mapping and deallocates all memory
  SUBROUTINE EquationsMapping_NormMappingFinalise(normMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingNormType) :: normMapping !<A pointer to the norm mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMapping_NormMappingFinalise",err,error,*999)

    normMapping%normNumber=0
    NULLIFY(normMapping%normVariable)
       
    EXITS("EquationsMapping_NormMappingFinalise")
    RETURN
999 ERRORSEXITS("EquationsMapping_NormMappingFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMapping_NormMappingFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the scalar equations mapping norm mapping
  SUBROUTINE EquationsMapping_NormMappingInitialise(normMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingNormType) :: normMapping !<The norm mapping to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsMapping_NormMappingInitialise",err,error,*999)

    normMapping%normNumber=0
    NULLIFY(normMapping%normVariable)
    normMapping%normCoefficient=1.0_DP
    
    EXITS("EquationsMapping_NormMappingInitialise")
    RETURN
999 ERRORSEXITS("EquationsMapping_NormMappingInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMapping_NormMappingInitialise

  !
  !================================================================================================================================
  !

  !>Finalises the scalar equations mapping norm mappings and deallocates all memory
  SUBROUTINE EquationsMapping_NormMappingsFinalise(normMappings,err,error,*)

    !Argument variables
    TYPE(EquationsMappingNormsType), POINTER :: normMappings !<A pointer to the norm mappings to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: normIdx

    ENTERS("EquationsMapping_NormMappingsFinalise",err,error,*999)

    IF(ASSOCIATED(normMappings)) THEN
      IF(ALLOCATED(normMappings%norms)) THEN
        DO normIdx=1,SIZE(normMappings%norms,1)
          CALL EquationsMapping_NormMappingFinalise(normMappings%norms(normIdx),err,error,*999)
        ENDDO !normIdx
        DEALLOCATE(normMappings%norms)
      ENDIF
      DEALLOCATE(normMappings)
    ENDIF
       
    EXITS("EquationsMapping_NormMappingsFinalise")
    RETURN
999 ERRORSEXITS("EquationsMapping_NormMappingsFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMapping_NormMappingsFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the scalar equations mapping norm mappings
  SUBROUTINE EquationsMapping_NormMappingsInitialise(scalarMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingScalarType), POINTER :: scalarMapping !<A pointer to the scalar equations mapping to initialise the norm mappings for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("EquationsMapping_NormMappingsInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(scalarMapping)) CALL FlagError("Scalar equations mapping is not associated.",err,error,*998)
    IF(ASSOCIATED(scalarMapping%normMappings)) &
      & CALL FlagError("Scalar equations mapping norm mappings is already associated.",err,error,*998)
    
    ALLOCATE(scalarMapping%normMappings,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate scalar equations mapping norm mappings.",err,error,*999)
    scalarMapping%normMappings%scalarMapping=>scalarMapping
    scalarMapping%normMappings%numberOfNorms=0
    
    EXITS("EquationsMapping_NormMappingsInitialise")
    RETURN
999 CALL EquationsMapping_NormMappingsFinalise(scalarMapping%normMappings,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMapping_NormMappingsInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMapping_NormMappingsInitialise

  !
  !================================================================================================================================
  !

  !>Finalises the scalar equations mapping quadratic form mapping and deallocates all memory
  SUBROUTINE EquationsMapping_QuadraticMappingFinalise(quadraticMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingQuadraticType) :: quadraticMapping !<A pointer to the quadratic mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMapping_QuadraticMappingFinalise",err,error,*999)

    quadraticMapping%quadraticNumber=0
    NULLIFY(quadraticMapping%quadraticVariables(1)%ptr)
    NULLIFY(quadraticMapping%quadraticVariables(2)%ptr)
    quadraticMapping%quadraticCoefficient=1.0_DP
       
    EXITS("EquationsMapping_QuadraticMappingFinalise")
    RETURN
999 ERRORSEXITS("EquationsMapping_QuadraticMappingFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMapping_QuadraticMappingFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the scalar equations mapping quadratic form mapping
  SUBROUTINE EquationsMapping_QuadraticMappingInitialise(quadraticMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingQuadraticType) :: quadraticMapping !<The quadratic form mapping to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsMapping_QuadraticMappingInitialise",err,error,*999)

    quadraticMapping%quadraticNumber=0
    NULLIFY(quadraticMapping%quadraticVariables(1)%ptr)
    NULLIFY(quadraticMapping%quadraticVariables(2)%ptr)
    quadraticMapping%quadraticCoefficient=1.0_DP
    
    EXITS("EquationsMapping_QauadraticMappingInitialise")
    RETURN
999 ERRORSEXITS("EquationsMapping_QuadraticMappingInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMapping_QuadraticMappingInitialise

  !
  !================================================================================================================================
  !

  !>Finalises the scalar equations mapping quadratic form mappings and deallocates all memory
  SUBROUTINE EquationsMapping_QuadraticMappingsFinalise(quadraticMappings,err,error,*)

    !Argument variables
    TYPE(EquationsMappingQuadraticsType), POINTER :: quadraticMappings !<A pointer to the quadratic forms mappings to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: quadraticIdx

    ENTERS("EquationsMapping_QuadraticMappingsFinalise",err,error,*999)

    IF(ASSOCIATED(quadraticMappings)) THEN
      IF(ALLOCATED(quadraticMappings%quadratics)) THEN
        DO quadraticIdx=1,SIZE(quadraticMappings%quadratics,1)
          CALL EquationsMapping_QuadraticMappingFinalise(quadraticMappings%quadratics(quadraticIdx),err,error,*999)
        ENDDO !quadraticIdx
        DEALLOCATE(quadraticMappings%quadratics)
      ENDIF
      DEALLOCATE(quadraticMappings)
    ENDIF
       
    EXITS("EquationsMapping_QuadraticMappingsFinalise")
    RETURN
999 ERRORSEXITS("EquationsMapping_QuadraticMappingsFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMapping_QuadraticMappingsFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the scalar equations mapping quadratic form mappings
  SUBROUTINE EquationsMapping_QuadraticMappingsInitialise(scalarMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingScalarType), POINTER :: scalarMapping !<A pointer to the scalar equations mapping to initialise the quadratic form mappings for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("EquationsMapping_QuadraticMappingsInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(scalarMapping)) CALL FlagError("Scalar equations mapping is not associated.",err,error,*998)
    IF(ASSOCIATED(scalarMapping%quadraticMappings)) &
      & CALL FlagError("Scalar equations mapping quadratic mappings is already associated.",err,error,*998)
    
    ALLOCATE(scalarMapping%quadraticMappings,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate scalar equations mapping quadratic mappings.",err,error,*999)
    scalarMapping%quadraticMappings%scalarMapping=>scalarMapping
    scalarMapping%quadraticMappings%numberOfQuadratics=0
    
    EXITS("EquationsMapping_QuadraticMappingsInitialise")
    RETURN
999 CALL EquationsMapping_QuadraticMappingsFinalise(scalarMapping%quadraticMappings,dummyErr,dummyError,*998)
998 ERRORS("EquationsMapping_QuadraticMappingsInitialise",err,error)
    EXITS("EquationsMapping_QuadraticMappingsInitialise")
    RETURN 1
    
  END SUBROUTINE EquationsMapping_QuadraticMappingsInitialise

  !
  !================================================================================================================================
  !

  !>Sets the coefficient applied to the equations set residual vector.
  SUBROUTINE EquationsMapping_ResidualCoeffSet(vectorMapping,residualCoefficient,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations mapping to set
    REAL(DP), INTENT(IN) :: residualCoefficient!<The coefficient applied to the equations set residual vector.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(EquationsMappingVectorCreateValuesCacheType), POINTER :: createValuesCache

    ENTERS("EquationsMapping_ResidualCoeffSet",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector equations mapping is not associated",err,error,*999)
    IF(vectorMapping%vectorMappingFinished) CALL FlagError("Vector equations mapping has been finished.",err,error,*999)

    NULLIFY(createValuesCache)
    CALL EquationsMappingVector_CreateValuesCacheGet(vectorMapping,createValuesCache,err,error,*999)

    createValuesCache%residualCoefficient=residualCoefficient
      
    EXITS("EquationsMapping_ResidualCoeffSet")
    RETURN
999 ERRORSEXITS("EquationsMapping_ResidualCoeffSet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMapping_ResidualCoeffSet

  !
  !================================================================================================================================
  !

  !>Sets the mapping between a dependent field variable and the equations residual vector.
  SUBROUTINE EquationsMapping_ResidualVariableTypesSet(vectorMapping,residualVariableTypes,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations mapping to set
    INTEGER(INTG), INTENT(IN) :: residualVariableTypes(:) !<residualVariableTypes(variableIdx). The variableIdx'th variable type associated with the equations set residual vector. The first variable type must correspond to the diagonal terms in the full solver Jacobian so that the solver mapping can use boundary conditions on this first variable to decide whether to keep rows.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx,variableIdx,numberOfResidualVariables,residualVariableType
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    TYPE(FIELD_TYPE), POINTER :: dependentField
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsMapping_ResidualVariableTypesSet",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector equations mapping is not associated",err,error,*999)
    IF(vectorMapping%vectorMappingFinished) CALL FlagError("Vector equations mapping have been finished.",err,error,*999)

    NULLIFY(createValuesCache)
    CALL EquationsMappingVector_CreateValuesCacheGet(vectorMapping,createValuesCache,err,error,*999)
    NULLIFY(vectorEquations)
    CALL EquationsMappingVector_VectorEquationsGet(vectorMapping,vectorEquations,err,error,*999)
    NULLIFY(equations)
    CALL EquationsVector_EquationsGet(vectorEquations,equations,err,error,*999)
    NULLIFY(equationsSet)
    CALL Equations_EquationsSetGet(equations,equationsSet,err,error,*999)
    
    numberOfResidualVariables=SIZE(residualVariableTypes,1)
    IF(numberOfResidualVariables/=createValuesCache%numberOfResidualVariables) THEN
      localError="Invalid number of variables. The number of residual variables " &
        & //TRIM(NumberToVString(numberOfResidualVariables,"*",err,error))//" should be "// &
        & TRIM(NumberToVString(createValuesCache%numberOfResidualVariables,"*",err,error))
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(equations%linearity/=EQUATIONS_NONLINEAR) CALL FlagError("The equations are not nonlinear equations.",err,error,*999)
    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
      
    !Check the residual variable types are not being used by other equations matrices or vectors
    DO variableIdx=1,numberOfResidualVariables
      residualVariableType=residualVariableTypes(variableIdx)
      SELECT CASE(equations%timeDependence)
      CASE(EQUATIONS_STATIC,EQUATIONS_QUASISTATIC)
        IF(createValuesCache%dynamicVariableType==residualVariableType) THEN
          localError="The specified residual variable type of "//TRIM(NumberToVString(residualVariableType,"*",err,error))// &
            & " is the same as the variable type for the dynamic matrices."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      CASE(EQUATIONS_FIRST_ORDER_DYNAMIC,EQUATIONS_SECOND_ORDER_DYNAMIC)
        IF(createValuesCache%dynamicVariableType/=residualVariableType) THEN
          localError="The specified residual variable type of "//TRIM(NumberToVString(residualVariableType,"*",err,error))// &
            & " is not the same as the variable type for the dynamic matrices."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      CASE DEFAULT
        localError="The equations time dependence type of "//TRIM(NumberToVString(equations%timeDependence,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      IF(createValuesCache%rhsVariableType==residualVariableType) THEN
        localError="The specified residual variable type of "//TRIM(NumberToVString(residualVariableType,"*",err,error))// &
          & " is the same as the variable type for the RHS vector."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      DO matrixIdx=1,createValuesCache%numberOfLinearMatrices
        IF(createValuesCache%linearMatrixVariableTypes(matrixIdx)==residualVariableType) THEN
          localError="The specified residual variable type of "//TRIM(NumberToVString(residualVariableType,"*",err,error))// &
            & " is the same as the variable type for linear matrix number "//TRIM(NumberToVString(matrixIdx,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDDO !matrixIdx
      !Check the residual variable number is defined on the dependent field
      IF(residualVariableType<1.OR.residualVariableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
        localError="The specified residual variable type of "//TRIM(NumberToVString(residualVariableType,"*",err,error))// &
          & " is invalid. The variable type must either be zero or >= 1 and <= "// &
          & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      IF(.NOT.ASSOCIATED(dependentField%VARIABLE_TYPE_MAP(residualVariableType)%ptr)) THEN
        localError="The specified residual variable type of "//TRIM(NumberToVString(residualVariableType,"*",err,error))// &
          & " is not defined on the dependent field."
        CALL FlagError(localError,err,error,*999)
      END IF      
    ENDDO !variableIdx
    
    createValuesCache%residualVariableTypes(1:numberOfResidualVariables)=residualVariableTypes(1:numberOfResidualVariables)
      
    EXITS("EquationsMapping_ResidualVariableTypesSet")
    RETURN
999 ERRORSEXITS("EquationsMapping_ResidualVariableTypesSet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMapping_ResidualVariableTypesSet
  
  !
  !================================================================================================================================
  !

  !>Sets the coefficient applied to the equations set RHS vector.
  SUBROUTINE EquationsMapping_RHSCoeffSet(vectorMapping,rhsCoefficient,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations mapping to set
    REAL(DP), INTENT(IN) :: rhsCoefficient!<The coefficient applied to the equations set RHS vector.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(EquationsMappingVectorCreateValuesCacheType), POINTER :: createValuesCache

    ENTERS("EquationsMapping_RHSCoeffSet",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector equations mapping is not associated",err,error,*999)
    IF(vectorMapping%vectorMappingFinished) CALL FlagError("Vector equations mapping has been finished.",err,error,*999)

    NULLIFY(createValuesCache)
    CALL EquationsMappingVector_CreateValuesCacheGet(vectorMapping,createValuesCache,err,error,*999)
    
    IF(createValuesCache%rhsVariableType==0) &
      & CALL FlagError("The equations mapping RHS variable type has not been set.",err,error,*999)

    vectorMapping%createValuesCache%rhsCoefficient=rhsCoefficient
       
    EXITS("EquationsMapping_RHSCoeffSet")
    RETURN
999 ERRORSEXITS("EquationsMapping_RHSCoeffSet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMapping_RHSCoeffSet
  
  !
  !================================================================================================================================
  !

  !>Finalises the equations mapping RHS mapping and deallocates all memory
  SUBROUTINE EquationsMapping_RHSMappingFinalise(rhsMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping !<A pointer to the RHS mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMapping_RHSMappingFinalise",err,error,*999)

    IF(ASSOCIATED(rhsMapping)) THEN
      IF(ALLOCATED(rhsMapping%rhsDOFToEquationsRowMap)) DEALLOCATE(rhsMapping%rhsDOFToEquationsRowMap)
      IF(ALLOCATED(rhsMapping%equationsRowToRHSDOFMap)) DEALLOCATE(rhsMapping%equationsRowToRHSDOFMap)
      DEALLOCATE(rhsMapping)
    ENDIF
       
    EXITS("EquationsMapping_RHSMappingFinalise")
    RETURN
999 ERRORSEXITS("EquationsMapping_RHSMappingFinalise",err,error)
    RETURN 1
  END SUBROUTINE EquationsMapping_RHSMappingFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the equations mapping RHS mapping
  SUBROUTINE EquationsMapping_RHSMappingInitialise(vectorMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations mapping to initialise the RHS mapping for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("EquationsMapping_RHSMappingInitialise",err,error,*998)

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
       
    EXITS("EquationsMapping_RHSMappingInitialise")
    RETURN
999 CALL EquationsMapping_RHSMappingFinalise(vectorMapping%rhsMapping,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMapping_RHSMappingInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMapping_RHSMappingInitialise

  !
  !================================================================================================================================
  !

  !>Sets the mapping between a dependent field variable and the equations set rhs vector.
  SUBROUTINE EquationsMapping_RHSVariableTypeSet(vectorMapping,rhsVariableType,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations mapping to set
    INTEGER(INTG), INTENT(IN) :: rhsVariableType !<The variable type associated with the equations set rhs vector. If the problem does not have a rhs vector then the variable type on input should be zero.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    TYPE(FIELD_TYPE), POINTER :: dependentField
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsMapping_RHSVariableTypeSet",err,error,*999)

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
      DO matrixIdx=1,createValuesCache%numberOfResidualVariables
        IF(createValuesCache%residualVariableTypes(matrixIdx)==rhsVariableType) THEN
          localError="The specified RHS variable type of "//TRIM(NumberToVString(rhsVariableType,"*",err,error))// &
            & " is the same as the variable type for the residual vector."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDDO !matrixIdx
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
      IF(.NOT.ASSOCIATED(dependentField%VARIABLE_TYPE_MAP(rhsVariableType)%ptr)) THEN
        localError="The specified RHS variable type of "//TRIM(NumberToVString(rhsVariableType,"*",err,error))// &
          & " is not defined on the dependent field."
        CALL FlagError(localError,err,error,*999)
      ENDIF

      createValuesCache%rhsVariableType=rhsVariableType
      
    ENDIF
       
    EXITS("EquationsMapping_RHSVariableTypeSet")
    RETURN
999 ERRORSEXITS("EquationsMapping_RHSVariableTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMapping_RHSVariableTypeSet
  
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
      CALL EquationsMapping_FunctionMappingsFinalise(scalarMapping%functionMappings,err,error,*999)
      CALL EquationsMapping_NormMappingsFinalise(scalarMapping%normMappings,err,error,*999)
      CALL EquationsMapping_DotProductMappingsFinalise(scalarMapping%dotProductMappings,err,error,*999)
      CALL EquationsMapping_QuadraticMappingsFinalise(scalarMapping%quadraticMappings,err,error,*999)
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

  !>Sets the coefficient applied to the equations set source vectorQ.
  SUBROUTINE EquationsMapping_SourceCoeffSet(vectorMapping,sourceCoefficient,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations mapping to set
    REAL(DP), INTENT(IN) :: sourceCoefficient!<The coefficient applied to the equations set source vector.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(EquationsMappingVectorCreateValuesCacheType), POINTER :: createValuesCache

    ENTERS("EquationsMapping_SourceCoeffSet",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector equations mapping is not associated",err,error,*999)
    IF(vectorMapping%vectorMappingFinished) CALL FlagError("Vector equations mapping has been finished.",err,error,*999)

    NULLIFY(createValuesCache)
    CALL EquationsMappingVector_CreateValuesCacheGet(vectorMapping,createValuesCache,err,error,*999)
    
    IF(createValuesCache%sourceVariableType==0) &
      & CALL FlagError("The equations mapping source variable type has not been set.",err,error,*999)
    
    createValuesCache%sourceCoefficient=sourceCoefficient
       
    EXITS("EquationsMapping_SourceCoeffSet")
    RETURN
999 ERRORSEXITS("EquationsMapping_SourceCoeffSet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMapping_SourceCoeffSet
  
  !
  !================================================================================================================================
  !

  !>Finalises the equations mapping source mapping and deallocates all memory
  SUBROUTINE EquationsMapping_SourceMappingFinalise(sourceMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingSourceType), POINTER :: sourceMapping !<A pointer to the SOURCE mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMapping_SourceMappingFinalise",err,error,*999)

    IF(ASSOCIATED(sourceMapping)) THEN
      IF(ALLOCATED(sourceMapping%sourceDOFToEquationsRowMap)) DEALLOCATE(sourceMapping%sourceDOFToEquationsRowMap)
      IF(ALLOCATED(sourceMapping%equationsRowToSourceDOFMap)) DEALLOCATE(sourceMapping%equationsRowToSourceDOFMap)
      DEALLOCATE(sourceMapping)
    ENDIF
       
    EXITS("EquationsMapping_SourceMappingFinalise")
    RETURN
999 ERRORSEXITS("EquationsMapping_SourceMappingFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMapping_SourceMappingFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the equations mapping source mapping
  SUBROUTINE EquationsMapping_SourceMappingInitialise(vectorMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations mapping to initialise the source mapping for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("EquationsMapping_SourceMappingInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector equations mapping is not associated.",err,error,*998)
    IF(ASSOCIATED(vectorMapping%sourceMapping)) &
      & CALL FlagError("Equations mapping source mapping is already associated.",err,error,*998)

    ALLOCATE(vectorMapping%sourceMapping,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate equations mapping source mapping.",err,error,*999)
    vectorMapping%sourceMapping%vectorMapping=>vectorMapping        
    vectorMapping%sourceMapping%sourceVariableType=0
    NULLIFY(vectorMapping%sourceMapping%sourceVariable)
    NULLIFY(vectorMapping%sourceMapping%sourceVariableMapping)
    vectorMapping%sourceMapping%sourceCoefficient=1.0_DP
       
    EXITS("EquationsMapping_SourceMappingInitialise")
    RETURN
999 CALL EquationsMapping_SourceMappingFinalise(vectorMapping%sourceMapping,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMapping_SourceMappingInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMapping_SourceMappingInitialise

  !
  !================================================================================================================================
  !

  !>Sets the mapping between a source field variable and the equations set source vector.
  SUBROUTINE EquationsMapping_SourceVariableTypeSet(vectorMapping,sourceVariableType,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the vector equations mapping to set
    INTEGER(INTG), INTENT(IN) :: sourceVariableType !<The variable type associated with the equations set source vector. If the problem does not have a source vector then the variable type on input should be zero.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    TYPE(FIELD_TYPE), POINTER :: sourceField
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsMapping_SourceVariableTypeSet",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector equations mapping is not associated.",err,error,*999)
    IF(vectorMapping%vectorMappingFinished) CALL FlagError("Vector equations mapping has been finished.",err,error,*999)

    NULLIFY(createValuesCache)
    CALL EquationsMappingVector_CreateValuesCacheGet(vectorMapping,createValuesCache,err,error,*999)
    NULLIFY(vectorEquations)
    CALL EquationsMappingVector_VectorEquationsGet(vectorMapping,vectorEquations,err,error,*999)
    NULLIFY(equations)
    CALL EquationsVector_EquationsGet(vectorEquations,equations,err,error,*999)
    NULLIFY(equationsSet)
    CALL Equations_EquationsSetGet(equations,equationsSet,err,error,*999)

    IF(sourceVariableType==0) THEN
      createValuesCache%sourceVariableType=0
    ELSE
      NULLIFY(sourceField)
      CALL EquationsSet_SourceFieldGet(equationsSet,sourceField,err,error,*999)
      
      !Check the source variable type is defined on the source field
      IF(sourceVariableType<1.OR.sourceVariableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
        localError="The specified source variable type of "//TRIM(NumberToVString(sourceVariableType,"*",err,error))// &
          & " is invalid. The number must either be zero or >= 1 and <= "// &
          & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      IF(.NOT.ASSOCIATED(sourceField%VARIABLE_TYPE_MAP(sourceVariableType)%ptr)) THEN
        localError="The specified source variable type of "//TRIM(NumberToVString(sourceVariableType,"*",err,error))// &
          & " is not defined on the source field."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      
      createValuesCache%sourceVariableType=sourceVariableType
      
    ENDIF
    
    EXITS("EquationsMapping_SourceVariableTypeSet")
    RETURN
999 ERRORSEXITS("EquationsMapping_SourceVariableTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMapping_SourceVariableTypeSet
  
  !
  !================================================================================================================================
  !
!
  !>Finalise an equations mapping equations matrix map.
  SUBROUTINE EquationsMapping_VarToEquatsColumnMapFinalise(varToEquationsColumnMap,err,error,*)

    !Argument variables
    TYPE(VarToEquationsColumnMapType) :: varToEquationsColumnMap !<The variable dof to equations column map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsMapping_VarToEquatsColumnMapFinalise",err,error,*999)
    
    IF(ALLOCATED(varToEquationsColumnMap%columnDOF)) DEALLOCATE(varToEquationsColumnMap%columnDOF)
    
    EXITS("EquationsMapping_VarToEquatsColumnMapFinalise")
    RETURN
999 ERRORS("EquationsMapping_VarToEquatsColumnMapFinalise",err,error)    
    EXITS("EquationsMapping_VarToEquatsColumnMapFinalise")    
    RETURN 1
   
  END SUBROUTINE EquationsMapping_VarToEquatsColumnMapFinalise

  !
  !================================================================================================================================
  !

  !>Finalises a variable to equations Jacobian map and deallocates all memory.
  SUBROUTINE EquationsMapping_VarToEquatsJacobianMapFinalise(varToEquationsJacobianMap,err,error,*)

    !Argument variables
    TYPE(VarToEquationsJacobianMapType) :: varToEquationsJacobianMap !<The variable to equations Jacobian map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsMapping_VarToEquatsJacobianMapFinalise",err,error,*999)
    
    IF(ALLOCATED(varToEquationsJacobianMap%dofToColumnsMap)) DEALLOCATE(varToEquationsJacobianMap%dofToColumnsMap)
    IF(ALLOCATED(varToEquationsJacobianMap%dofToRowsMap)) DEALLOCATE(varToEquationsJacobianMap%dofToRowsMap)
    
    EXITS("EquationsMapping_VarToEquatsJacobianMapFinalise")
    RETURN
999 ERRORS("EquationsMapping_VarToEquatsJacobianMapFinalise",err,error)    
    EXITS("EquationsMapping_VarToEquatsJacobianMapFinalise")    
    RETURN 1
   
  END SUBROUTINE EquationsMapping_VarToEquatsJacobianMapFinalise

  !
  !================================================================================================================================
  !

  !>Initialises a variable to equations Jacobian map
  SUBROUTINE EquationsMapping_VarToEquatsJacobianMapInitialise(varToEquationsJacobianMap,err,error,*)

    !Argument variables
    TYPE(VarToEquationsJacobianMapType) :: varToEquationsJacobianMap !<The variable to equations Jacobian map to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsMapping_VarToEquatsJacobianMapInitialise",err,error,*999)
    
    varToEquationsJacobianMap%variableType=0
    NULLIFY(varToEquationsJacobianMap%variable)
    
    EXITS("EquationsMapping_VarToEquatsJacobianMapInitialise")
    RETURN
999 ERRORS("EquationsMapping_VarToEquatsJacobianMapInitialise",err,error)    
    EXITS("EquationsMapping_VarToEquatsJacobianMapInitialise")    
    RETURN 1
   
  END SUBROUTINE EquationsMapping_VarToEquatsJacobianMapInitialise

  !
  !================================================================================================================================
  !

  !>Finalises a variable to equations matrices map and deallocates all memory.
  SUBROUTINE EquationsMapping_VarToEquatsMatricesMapFinalise(varToEquationsMatricesMap,err,error,*)

    !Argument variables
    TYPE(VarToEquationsMatricesMapType) :: varToEquationsMatricesMap !<The variable to equations matrices map to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx

    ENTERS("EquationsMapping_VarToEquatsMatricesMapFinalise",err,error,*999)
    
    IF(ALLOCATED(varToEquationsMatricesMap%equationsMatrixNumbers)) DEALLOCATE(varToEquationsMatricesMap%equationsMatrixNumbers)
    IF(ALLOCATED(varToEquationsMatricesMap%dofToColumnsMaps)) THEN
      DO matrixIdx=1,SIZE(varToEquationsMatricesMap%dofToColumnsMaps,1)
        CALL EquationsMapping_VarToEquatsColumnMapFinalise(varToEquationsMatricesMap%dofToColumnsMaps(matrixIdx),err,error,*999)
      ENDDO !matrixIdx
      DEALLOCATE(varToEquationsMatricesMap%dofToColumnsMaps)
    ENDIF
    IF(ALLOCATED(varToEquationsMatricesMap%dofToRowsMap)) DEALLOCATE(varToEquationsMatricesMap%dofToRowsMap)
    
    EXITS("EquationsMapping_VarToEquatsMatricesMapFinalise")
    RETURN
999 ERRORS("EquationsMapping_VarToEquatsMatricesMapFinalise",err,error)    
    EXITS("EquationsMapping_VarToEquatsMatricesMapFinalise")    
    RETURN 1
   
  END SUBROUTINE EquationsMapping_VarToEquatsMatricesMapFinalise

  !
  !================================================================================================================================
  !

  !>Initialise an equations mapping equations matrix map.
  SUBROUTINE EquationsMapping_VarToEquatsMatricesMapInitialise(varToEquationsMatricesMap,err,error,*)

    !Argument variables
    TYPE(VarToEquationsMatricesMapType) :: varToEquationsMatricesMap !<The variable to equations matrices map to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsMapping_VarToEquatsMatricesMapInitialise",err,error,*999)

    varToEquationsMatricesMap%variableIndex=0
    varToEquationsMatricesMap%variableType=0
    NULLIFY(varToEquationsMatricesMap%variable)
    varToEquationsMatricesMap%numberOfEquationsMatrices=0
    
    EXITS("EquationsMapping_VarToEquatsMatricesMapInitialise")
    RETURN
999 ERRORS("EquationsMapping_VarToEquatsMatricesMapInitialise",err,error)    
    EXITS("EquationsMapping_VarToEquatsMatricesMapInitialise")    
    RETURN 1
   
  END SUBROUTINE EquationsMapping_VarToEquatsMatricesMapInitialise

  !
  !================================================================================================================================
  !
  
END MODULE EquationsMappingRoutines
