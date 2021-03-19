!> \file
!> \author Chris Bradley
!> \brief This module contains all interface mapping routines.
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

!>This module contains all interface mapping routines.
MODULE InterfaceMappingRoutines

  USE BaseRoutines
  USE FieldRoutines
  USE FieldAccessRoutines
  USE InputOutput
  USE InterfaceConditionAccessRoutines
  USE InterfaceEquationsAccessRoutines
  USE InterfaceMappingAccessRoutines
  USE InterfaceMatricesAccessRoutines
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

  PUBLIC InterfaceMapping_CreateFinish,InterfaceMapping_CreateStart

  PUBLIC InterfaceMapping_Destroy

  PUBLIC InterfaceMapping_LagrangeVariableTypeSet

  PUBLIC InterfaceMapping_MatrixCoefficientsSet

  PUBLIC InterfaceMapping_TransposeMatrixCoefficientsSet

  PUBLIC InterfaceMapping_MatricesColumnMeshIndicesSet,InterfaceMapping_MatricesRowMeshIndicesSet

  PUBLIC InterfaceMapping_NumberOfMatricesSet

  PUBLIC InterfaceMapping_MatricesTransposeSet

  PUBLIC InterfaceMapping_RHSCoefficientSet

  PUBLIC InterfaceMapping_RHSVariableTypeSet

CONTAINS

  !
  !================================================================================================================================
  !

  !>Calculates interface mapping
  SUBROUTINE InterfaceMapping_Calculate(interfaceMapping,err,error,*)
    
    !Argument variables
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping !<A pointer to the interface mapping to calculate the mapping for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: columnIdx,dofIdx,interfaceConditionMethod,matrixIdx,meshIdx,variableIdx,numberOfDOFs,numberOfGlobalDOFs, &
      & numberOfInterfaceMatrices,totalNumberOfDOFs,variableType
    REAL(DP) :: matrixCoefficient,transposeMatrixCoefficient
    LOGICAL :: hasTranspose
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainMappingType), POINTER :: domainMapping
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(FieldType), POINTER :: lagrangeField
    TYPE(FieldVariableType), POINTER :: fieldVariable,lagrangeVariable
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition
    TYPE(InterfaceDependentType), POINTER :: interfaceDependent
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations
    TYPE(InterfaceLagrangeType), POINTER :: interfaceLagrange
    TYPE(InterfaceMappingCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(InterfaceMappingRHSType), POINTER :: rhsMapping
    TYPE(InterfaceMatrixToVarMapType), POINTER :: interfaceMatrixToVarMap
    TYPE(VARYING_STRING) :: localError

    ENTERS("InterfaceMapping_Calculate",err,error,*999)

    IF(.NOT.ASSOCIATED(interfaceMapping)) CALL FlagError("Interface mapping is not associated.",err,error,*999)

    NULLIFY(createValuesCache)
    CALL InterfaceMapping_CreateValuesCacheGet(interfaceMapping,createValuesCache,err,error,*999)
    NULLIFY(interfaceEquations)
    CALL InterfaceMapping_InterfaceEquationsGet(interfaceMapping,interfaceEquations,err,error,*999)
    NULLIFY(interfaceCondition)
    CALL InterfaceEquations_InterfaceConditionGet(interfaceEquations,interfaceCondition,err,error,*999)
    CALL InterfaceCondition_MethodGet(interfaceCondition,interfaceConditionMethod,err,error,*999)
    SELECT CASE(interfaceConditionMethod)
    CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
      NULLIFY(interfaceLagrange)
      CALL InterfaceCondition_InterfaceLagrangeGet(interfaceCondition,interfaceLagrange,err,error,*999)
      NULLIFY(interfaceDependent)
      CALL InterfaceCondition_InterfaceDependentGet(interfaceCondition,interfaceDependent,err,error,*999)
      NULLIFY(lagrangeField)
      CALL InterfaceLagrange_LagrangeFieldGet(interfaceLagrange,lagrangeField,err,error,*999)
      !Set the Lagrange variable information
      NULLIFY(lagrangeVariable)
      CALL Field_VariableGet(lagrangeField,createValuesCache%lagrangeVariableType,lagrangeVariable,err,error,*999)
      interfaceMapping%lagrangeVariableType=createValuesCache%lagrangeVariableType
      interfaceMapping%lagrangeVariable=>lagrangeVariable
      !Set the number of columns in the interface matrices
      CALL FieldVariable_NumberOfDOFsGet(lagrangeVariable,interfaceMapping%numberOfColumns,err,error,*999)
      CALL FieldVariable_TotalNumberOfDOFsGet(lagrangeVariable,interfaceMapping%totalNumberOfColumns,err,error,*999)
      CALL FieldVariable_NumberOfGlobalDOFsGet(lagrangeVariable,interfaceMapping%numberOfGlobalColumns,err,error,*999)
      !Set the column dofs mapping
      CALL FieldVariable_DomainMappingGet(lagrangeVariable,interfaceMapping%columnDOFsMapping,err,error,*999)
      ALLOCATE(interfaceMapping%lagrangeDOFToColumnMap(interfaceMapping%totalNumberOfColumns),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate Lagrange DOF to column map.",err,error,*999)
      !1-1 mapping for now
      DO dofIdx=1,interfaceMapping%totalNumberOfColumns
        CALL DomainMapping_LocalToGlobalGet(interfaceMapping%columnDOFsMapping,dofIdx,columnIdx,err,error,*999)
        interfaceMapping%lagrangeDOFToColumnMap(dofIdx)=columnIdx
      ENDDO !dofIdx
      !Set the number of interface matrices
      interfaceMapping%numberOfInterfaceMatrices=createValuesCache%numberOfInterfaceMatrices
      ALLOCATE(interfaceMapping%interfaceMatrixToVarMaps(interfaceMapping%numberOfInterfaceMatrices),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate interface matrix rows to variable maps.",err,error,*999)
      !Loop over the interface matrices and calculate the row mappings
      SELECT CASE(interfaceConditionMethod)
      CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
        numberOfInterfaceMatrices=interfaceMapping%numberOfInterfaceMatrices
      CASE(INTERFACE_CONDITION_PENALTY_METHOD)
        !Number of interface matrices whose rows/columns are related to Dependent/Lagrange variables and not
        !Lagrange/Lagrange variables (last interface matrix is Lagrange/Lagrange (Penalty matrix)
        numberOfInterfaceMatrices=interfaceMapping%numberOfInterfaceMatrices-1 
      END SELECT
      DO matrixIdx=1,numberOfInterfaceMatrices
        !Initialise and setup the interface matrix
        NULLIFY(interfaceMapping%interfaceMatrixToVarMaps(matrixIdx)%ptr)
        CALL InterfaceMapping_MatrixToVarMapInitialise(interfaceMapping,matrixIdx,err,error,*999)
        CALL InterfaceMappingCVC_RowVariableIndexGet(createValuesCache,matrixIdx,meshIdx,err,error,*999)
        NULLIFY(equationsSet)
        NULLIFY(fieldVariable)
        DO variableIdx=1,interfaceDependent%numberOfDependentVariables
          IF(interfaceDependent%variableMeshIndices(variableIdx)==meshIdx) THEN
            CALL InterfaceDependent_EquationsSetGet(interfaceDependent,variableIdx,equationsSet,err,error,*999)
            CALL InterfaceDependent_DependentVariableGet(interfaceDependent,variableIdx,fieldVariable,err,error,*999)
            EXIT
          ENDIF
        ENDDO !variableIdx
        IF(.NOT.ASSOCIATED(equationsSet)) THEN
          localError="Equations set for mesh index "//TRIM(NumberToVString(meshIdx,"*",err,error))//" could not be found."
          CALL FlagError(localError,err,error,*999)
        ENDIF          
        IF(.NOT.ASSOCIATED(fieldVariable)) THEN
          localError="Dependent variable for mesh index "//TRIM(NumberToVString(meshIdx,"*",err,error))//" could not be found."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        CALL FieldVariable_VariableTypeGet(fieldVariable,variableType,err,error,*999)
        CALL FieldVariable_NumberOfDOFsGet(fieldVariable,numberOfDOFs,err,error,*999)
        CALL FieldVariable_TotalNumberOfDOFsGet(fieldVariable,totalNumberOfDOFs,err,error,*999)
        CALL FieldVariable_NumberOfGlobalDOFsGet(fieldVariable,numberOfGlobalDOFs,err,error,*999)
        NULLIFY(domainMapping)
        CALL FieldVariable_DomainMappingGet(fieldVariable,domainMapping,err,error,*999)
        CALL InterfaceMappingCVC_MatrixCoefficientGet(createValuesCache,matrixIdx,matrixCoefficient,err,error,*999)
        CALL InterfaceMappingCVC_TransposeMatrixCoefficientGet(createValuesCache,matrixIdx,transposeMatrixCoefficient, &
          & err,error,*999)
        CALL InterfaceMappingCVS_HasTransposeGet(createValuesCache,matrixIdx,hasTranspose,err,error,*999)
        NULLIFY(interfaceMatrixToVarMap)
        CALL InterfaceMapping_InterfaceMatrixToVarMapGet(interfaceMapping,matrixIdx,interfaceMatrixToVarMap,err,error,*999)
        interfaceMatrixToVarMap%equationsSet=>equationsSet
        interfaceMatrixToVarMap%variableType=variableType
        interfaceMatrixToVarMap%variable=>fieldVariable
        interfaceMatrixToVarMap%meshIndex=meshIdx
        interfaceMatrixToVarMap%matrixCoefficient=matrixCoefficient
        interfaceMatrixToVarMap%hasTranspose=hasTranspose
        IF(hasTranspose) interfaceMatrixToVarMap%transposeMatrixCoefficient=transposeMatrixCoefficient
        !Set the number of rows
        interfaceMatrixToVarMap%numberOfRows=numberOfDOFs
        interfaceMatrixToVarMap%totalNumberOfRows=totalNumberOfDOFs
        interfaceMatrixToVarMap%numberOfGlobalRows=numberOfGlobalDOFs
        !Set the row mapping
        interfaceMatrixToVarMap%rowDOFsMapping=>domainMapping
        ALLOCATE(interfaceMatrixToVarMap%variableDOFToRowMap(totalNumberOfDOFs),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate variable DOF to row map.",err,error,*999)
        !1-1 mapping for now
        DO dofIdx=1,totalNumberOfDOFs
          interfaceMatrixToVarMap%variableDOFToRowMap(dofIdx)=dofIdx
        ENDDO !dofIdx
      ENDDO !matrixIdx
      IF(interfaceConditionMethod==INTERFACE_CONDITION_PENALTY_METHOD) THEN
        !Sets up the Lagrange-(Penalty) interface matrix mapping and calculate the row mappings
        matrixIdx = interfaceMapping%numberOfInterfaceMatrices !last of the interface matrices
        !Initialise and setup the interface matrix
        NULLIFY(interfaceMapping%interfaceMatrixToVarMaps(matrixIdx)%ptr)
        CALL InterfaceMapping_MatrixToVarMapInitialise(interfaceMapping,matrixIdx,err,error,*999)
        CALL InterfaceMappingCVC_RowVariableIndexGet(createValuesCache,matrixIdx,meshIdx,err,error,*999)
        NULLIFY(lagrangeVariable)
        CALL Field_VariableGet(lagrangeField,createValuesCache%lagrangeVariableType,lagrangeVariable,err,error,*999)
        CALL FieldVariable_VariableTypeGet(lagrangeVariable,variableType,err,error,*999)
        CALL FieldVariable_NumberOfDOFsGet(lagrangeVariable,numberOfDOFs,err,error,*999)
        CALL FieldVariable_TotalNumberOfDOFsGet(lagrangeVariable,totalNumberOfDOFs,err,error,*999)
        CALL FieldVariable_NumberOfGlobalDOFsGet(fieldVariable,numberOfGlobalDOFs,err,error,*999)
        NULLIFY(domainMapping)
        CALL FieldVariable_DomainMappingGet(lagrangeVariable,domainMapping,err,error,*999)
        CALL InterfaceMappingCVC_MatrixCoefficientGet(createValuesCache,matrixIdx,matrixCoefficient,err,error,*999)
        CALL InterfaceMappingCVC_TransposeMatrixCoefficientGet(createValuesCache,matrixIdx,transposeMatrixCoefficient, &
          & err,error,*999)
        CALL InterfaceMappingCVS_HasTransposeGet(createValuesCache,matrixIdx,hasTranspose,err,error,*999)
        interfaceMatrixToVarMap%interfaceEquations=>interfaceEquations
        interfaceMatrixToVarMap%variableType=variableType
        interfaceMatrixToVarMap%variable=>lagrangeVariable
        interfaceMatrixToVarMap%meshIndex=meshIdx
        interfaceMatrixToVarMap%matrixCoefficient=matrixCoefficient
        interfaceMatrixToVarMap%hasTranspose=hasTranspose
        !Set the number of rows
        interfaceMatrixToVarMap%numberOfRows=numberOfDOFs
        interfaceMatrixToVarMap%totalNumberOfRows=totalNumberOfDOFs
        interfaceMatrixToVarMap%numberOfGlobalRows=numberOfGlobalDOFs
        !Set the row mapping
        interfaceMatrixToVarMap%rowDOFsMapping=>domainMapping
        ALLOCATE(interfaceMatrixToVarMap%variableDOFToRowMap(totalNumberOfDOFs),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate variable DOF to row map.",err,error,*999)
        !1-1 mapping for now
        DO dofIdx=1,totalNumberOfDOFs
          interfaceMatrixToVarMap%variableDOFToRowMap(dofIdx)=dofIdx
        ENDDO !dofIdx
      ENDIF
      !Calculate RHS mappings
      IF(createValuesCache%rhsLagrangeVariableType/=0) THEN
        CALL InterfaceMapping_RHSMappingInitialise(interfaceMapping,err,error,*999)
        NULLIFY(rhsMapping)
        CALL InterfaceMapping_RHSMappingGet(interfaceMapping,rhsMapping,err,error,*999)
        rhsMapping%rhsVariableType=createValuesCache%rhsLagrangeVariableType
        NULLIFY(lagrangeVariable)
        CALL Field_VariableGet(lagrangeField,createValuesCache%rhsLagrangeVariableType,lagrangeVariable,err,error,*999)
        CALL FieldVariable_NumberOfDOFsGet(lagrangeVariable,numberOfDOFs,err,error,*999)
        CALL FieldVariable_TotalNumberOfDOFsGet(lagrangeVariable,totalNumberOfDOFs,err,error,*999)
        CALL FieldVariable_NumberOfGlobalDOFsGet(fieldVariable,numberOfGlobalDOFs,err,error,*999)
        NULLIFY(domainMapping)
        CALL FieldVariable_DomainMappingGet(lagrangeVariable,domainMapping,err,error,*999)         
        rhsMapping%rhsVariable=>lagrangeVariable
        rhsMapping%rhsVariableMapping=>domainMapping
        rhsMapping%rhsCoefficient=createValuesCache%rhsCoefficient
        !Allocate and set up the row mappings
        ALLOCATE(rhsMapping%rhsDOFToInterfaceRowMap(totalNumberOfDOFs),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate RHS DOF to interface row map.",err,error,*999)
        ALLOCATE(rhsMapping%interfaceRowToRHSDOFMap(interfaceMapping%totalNumberOfColumns),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate interface row to RHS DOF map.",err,error,*999)
        DO dofIdx=1,totalNumberOfDOFs
          !1-1 mapping for now
          columnIdx=dofIdx
          rhsMapping%rhsDOFToInterfaceRowMap(dofIdx)=columnIdx
        ENDDO !dofIdx
        DO columnIdx=1,interfaceMapping%totalNumberOfColumns
          !1-1 mapping for now
          dofIdx=columnIdx
          rhsMapping%interfaceRowToRHSDOFMap(columnIdx)=dofIdx
        ENDDO !columnIdx
      ENDIF
    CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The interface condition method of "//TRIM(NumberToVString(interfaceConditionMethod,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("InterfaceMapping_Calculate")
    RETURN
999 ERRORSEXITS("InterfaceMapping_Calculate",err,error)
    RETURN 1
   
  END SUBROUTINE InterfaceMapping_Calculate

  !
  !================================================================================================================================
  !

  !>Finishes the creation of interface mapping
  SUBROUTINE InterfaceMapping_CreateFinish(interfaceMapping,err,error,*)
    
    !Argument variables
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping !<A pointer to the interface mapping to finish the creation of.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("InterfaceMapping_CreateFinish",err,error,*999)

    CALL InterfaceMapping_AssertNotFinished(interfaceMapping,err,error,*999)

    !Calculate the equations mapping and clean up
    CALL InterfaceMapping_Calculate(interfaceMapping,err,error,*999)
    CALL InterfaceMapping_CreateValuesCacheFinalise(interfaceMapping%createValuesCache,err,error,*999)
    interfaceMapping%interfaceMappingFinished=.TRUE.
    
    EXITS("InterfaceMapping_CreateFinish")
    RETURN
999 ERRORSEXITS("InterfaceMapping_CreateFinish",err,error)    
    RETURN 1
   
  END SUBROUTINE InterfaceMapping_CreateFinish

  !
  !================================================================================================================================
  !

  !>Starts the process of creating the interface mapping for interface equations.
  SUBROUTINE InterfaceMapping_CreateStart(interfaceEquations,interfaceMapping,err,error,*)
    
    !Argument variables
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations !<A pointer to the interface equations to create the mapping for.
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping !<On exit, a pointer to the created interface mapping. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("InterfaceMapping_CreateStart",err,error,*998)

    IF(ASSOCIATED(interfaceMapping)) CALL FlagError("Interface mapping is already associated.",err,error,*998)
    CALL InterfaceEquations_AssertIsFinished(interfaceEquations,err,error,*999)

    CALL InterfaceMapping_Initialise(interfaceEquations,err,error,*999)
    interfaceMapping=>interfaceEquations%interfaceMapping
   
    EXITS("InterfaceMapping_CreateStart")
    RETURN
999 NULLIFY(interfaceMapping)
998 ERRORSEXITS("InterfaceMapping_CreateStart",err,error)    
    RETURN 1
   
  END SUBROUTINE InterfaceMapping_CreateStart

  !
  !================================================================================================================================
  !

  !>Finalises an interface mapping create values cache and deallocates all memory
  SUBROUTINE InterfaceMapping_CreateValuesCacheFinalise(createValuesCache,err,error,*)

    !Argument variables
    TYPE(InterfaceMappingCreateValuesCacheType), POINTER :: createValuesCache !<A pointer to the create values cache to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("InterfaceMapping_CreateValuesCacheFinalise",err,error,*999)

    IF(ASSOCIATED(createValuesCache)) THEN
      IF(ALLOCATED(createValuesCache%matrixCoefficients)) DEALLOCATE(createValuesCache%matrixCoefficients)
      IF(ALLOCATED(createValuesCache%transposeMatrixCoefficients)) DEALLOCATE(createValuesCache%transposeMatrixCoefficients)
      IF(ALLOCATED(createValuesCache%hasTranspose)) DEALLOCATE(createValuesCache%hasTranspose)
      IF(ALLOCATED(createValuesCache%matrixRowFieldVariableIndices)) DEALLOCATE(createValuesCache%matrixRowFieldVariableIndices)
      IF(ALLOCATED(createValuesCache%matrixColFieldVariableIndices)) DEALLOCATE(createValuesCache%matrixColFieldVariableIndices)
      DEALLOCATE(createValuesCache)
    ENDIF
       
    EXITS("InterfaceMapping_CreateValuesCacheFinalise")
    RETURN
999 ERRORSEXITS("InterfaceMapping_CreateValuesCacheFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMapping_CreateValuesCacheFinalise

  !
  !================================================================================================================================
  !

  !>Initialises an interface mapping create values cache 
  SUBROUTINE InterfaceMapping_CreateValuesCacheInitialise(interfaceMapping,err,error,*)

    !Argument variables
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping !<A pointer to the interface mapping to create the create values cache for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,interfaceConditionMethod,variableIdx,variableTypeIdx,variableTypeIdx2
    TYPE(FieldType), POINTER :: lagrangeField
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition
    TYPE(InterfaceDependentType), POINTER :: interfaceDependent
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations
    TYPE(VARYING_STRING) :: dummyError,localError

    ENTERS("InterfaceMapping_CreateValuesCacheInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(interfaceMapping)) CALL FlagError("Interface mapping is not associated.",err,error,*998)
    IF(ASSOCIATED(interfaceMapping%createValuesCache)) &
      & CALL FlagError("Interface mapping create values cache is already associated.",err,error,*998)

    NULLIFY(interfaceEquations)
    CALL InterfaceMapping_InterfaceEquationsGet(interfaceMapping,interfaceEquations,err,error,*999)
    NULLIFY(interfaceCondition)
    CALL InterfaceEquations_InterfaceConditionGet(interfaceEquations,interfaceCondition,err,error,*999)
    CALL InterfaceCondition_MethodGet(interfaceCondition,interfaceConditionMethod,err,error,*999)
    !Allocate and initialise the create values cache
    ALLOCATE(interfaceMapping%createValuesCache,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate interface mapping create values cache.",err,error,*999)
    interfaceMapping%createValuesCache%numberOfInterfaceMatrices=0
    interfaceMapping%createValuesCache%lagrangeVariableType=0
    interfaceMapping%createValuesCache%rhsLagrangeVariableType=0
    interfaceMapping%createValuesCache%rhsCoefficient=0.0_DP
    !Set the default interface mapping in the create values cache
    !First calculate how many interface matrices we have and set the variable types
    SELECT CASE(interfaceConditionMethod)
    CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
      NULLIFY(lagrangeField)
      CALL InterfaceCondition_LagrangeFieldGet(interfaceCondition,lagrangeField,err,error,*999)
      NULLIFY(interfaceDependent)
      CALL InterfaceCondition_InterfaceDependentGet(interfaceCondition,interfaceDependent,err,error,*999)
      SELECT CASE(interfaceConditionMethod)
      CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
        !Default the number of interface matrices to the number of added dependent variables
        interfaceMapping%createValuesCache%numberOfInterfaceMatrices=interfaceDependent%numberOfDependentVariables
      CASE(INTERFACE_CONDITION_PENALTY_METHOD)
        !Default the number of interface matrices to the number of added dependent variables plus a single Lagrange variable
        interfaceMapping%createValuesCache%numberOfInterfaceMatrices=interfaceDependent%numberOfDependentVariables+1
      END SELECT
      !Default the Lagrange variable to the first Lagrange variable
      interfaceMapping%createValuesCache%lagrangeVariableType=0
      DO variableTypeIdx=1,FIELD_NUMBER_OF_VARIABLE_TYPES
        IF(ASSOCIATED(lagrangeField%variableTypeMap(variableTypeIdx)%ptr)) THEN
          interfaceMapping%createValuesCache%lagrangeVariableType=variableTypeIdx
          EXIT
        ENDIF
      ENDDO !variableTypeIdx
      IF(interfaceMapping%createValuesCache%lagrangeVariableType==0) &
        & CALL FlagError("Could not find a Lagrange variable type in the Lagrange field.",err,error,*999)
      !Default the RHS Lagrange variable to zero
      interfaceMapping%createValuesCache%rhsLagrangeVariableType=0
      IF(interfaceMapping%createValuesCache%rhsLagrangeVariableType==0) &
        & CALL FlagError("Could not find a RHS Lagrange variable type in the Lagrange field.",err,error,*999)
      ALLOCATE(interfaceMapping%createValuesCache%matrixCoefficients(interfaceMapping%createValuesCache% &
        & numberOfInterfaceMatrices),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate create values cache matrix coefficients.",err,error,*999)
      ALLOCATE(interfaceMapping%createValuesCache%transposeMatrixCoefficients(interfaceMapping%createValuesCache% &
        & numberOfInterfaceMatrices),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate create values cache transpose matrix coefficients.",err,error,*999)
       !Default the interface matrices coefficients to add.
      interfaceMapping%createValuesCache%matrixCoefficients=1.0_DP
      interfaceMapping%createValuesCache%transposeMatrixCoefficients=1.0_DP
      interfaceMapping%createValuesCache%rhsCoefficient=1.0_DP
      ALLOCATE(interfaceMapping%createValuesCache%hasTranspose(interfaceMapping%createValuesCache% &
        & numberOfInterfaceMatrices),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate create values cache has transpose.",err,error,*999)
      !Default the interface matrices to all have a transpose
      interfaceMapping%createValuesCache%hasTranspose=.TRUE.
      ALLOCATE(interfaceMapping%createValuesCache%matrixRowFieldVariableIndices(interfaceMapping%createValuesCache% &
        & numberOfInterfaceMatrices),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate create values cache matrix row field variable indexes.",err,error,*999)
      !Default the interface matrices to be in mesh index order.
      DO variableIdx=1,interfaceDependent%numberOfDependentVariables
        interfaceMapping%createValuesCache%matrixRowFieldVariableIndices(variableIdx)=variableIdx
      ENDDO !variableIdx
      !The pointers below have been checked for association above.
      IF(interfaceConditionMethod==INTERFACE_CONDITION_PENALTY_METHOD) THEN
        !Default the interface matrix (Penalty) to have no transpose
        interfaceMapping%createValuesCache%hasTranspose(interfaceMapping% &
          & createValuesCache%numberOfInterfaceMatrices)=.FALSE.
        !Default the interface matrices to be in mesh index order (and set Penalty matrix (last interface matrix) to be the 
        !first Lagrange variable).
        interfaceMapping%createValuesCache%matrixRowFieldVariableIndices(interfaceDependent%numberOfDependentVariables+1)=1
      ENDIF
    CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The interface equations method of "//TRIM(NumberToVString(interfaceCondition%method,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("InterfaceMapping_CreateValuesCacheInitialise")
    RETURN
999 CALL InterfaceMapping_CreateValuesCacheFinalise(interfaceMapping%createValuesCache,dummyErr,dummyError,*998)
998 ERRORS("InterfaceMapping_CreateValuesCacheInitialise",err,error)
    EXITS("InterfaceMapping_CreateValuesCacheInitialise")
    RETURN 1
    
  END SUBROUTINE InterfaceMapping_CreateValuesCacheInitialise

  !
  !================================================================================================================================
  !

  !>Destroys an interface mapping.
  SUBROUTINE InterfaceMapping_Destroy(interfaceMapping,err,error,*)

    !Argument variables
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping !<A pointer the interface mapping to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("InterfaceMapping_Destroy",err,error,*999)

    IF(.NOT.ASSOCIATED(interfaceMapping)) CALL FlagError("Interface mapping is not associated.",err,error,*999)
    
    CALL InterfaceMapping_Finalise(interfaceMapping,err,error,*999)
        
    EXITS("InterfaceMapping_Destroy")
    RETURN
999 ERRORSEXITS("InterfaceMapping_Destroy",err,error)
    RETURN 1
   
  END SUBROUTINE InterfaceMapping_Destroy

  !
  !================================================================================================================================
  !

  !>Finalises the interface mapping and deallocates all memory.
  SUBROUTINE InterfaceMapping_Finalise(interfaceMapping,err,error,*)

    !Argument variables
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping !<A pointer to the interface mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx

    ENTERS("InterfaceMapping_Finalise",err,error,*999)

    IF(ASSOCIATED(interfaceMapping)) THEN
      IF(ASSOCIATED(interfaceMapping%lagrangeDOFToColumnMap)) DEALLOCATE(interfaceMapping%lagrangeDOFToColumnMap)
      IF(ALLOCATED(interfaceMapping%interfaceMatrixToVarMaps)) THEN
        DO matrixIdx=1,SIZE(interfaceMapping%interfaceMatrixToVarMaps,1)
          CALL InterfaceMapping_MatrixToVarMapFinalise(interfaceMapping%interfaceMatrixToVarMaps(matrixIdx)%ptr, &
            & err,error,*999)
        ENDDO !matrixIdx
        DEALLOCATE(interfaceMapping%interfaceMatrixToVarMaps)
      ENDIF
      CALL InterfaceMapping_RHSMappingFinalise(interfaceMapping%rhsMapping,err,error,*999)
      CALL InterfaceMapping_CreateValuesCacheFinalise(interfaceMapping%createValuesCache,err,error,*999)
      DEALLOCATE(interfaceMapping)
    ENDIF
       
    EXITS("InterfaceMapping_Finalise")
    RETURN
999 ERRORSEXITS("InterfaceMapping_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMapping_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises the interface mapping and deallocates all memory.
  SUBROUTINE InterfaceMapping_Initialise(interfaceEquations,err,error,*)

    !Argument variables
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations !<A pointer to the interface equations to initialise the interface mapping for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("InterfaceMapping_Initialise",err,error,*998)

    IF(.NOT.ASSOCIATED(interfaceEquations)) CALL FlagError("Interface equations is not associated.",err,error,*998)
    IF(ASSOCIATED(interfaceEquations%interfaceMapping)) CALL FlagError("Interface mapping is already associated.",err,error,*998)
    
    ALLOCATE(interfaceEquations%interfaceMapping,STAT=err)    
    IF(err/=0) CALL FlagError("Could not allocate interface equations interface mapping.",err,error,*999)
    interfaceEquations%interfaceMapping%interfaceEquations=>interfaceEquations
    interfaceEquations%interfaceMapping%interfaceMappingFinished=.FALSE.
    interfaceEquations%interfaceMapping%lagrangeVariableType=0
    NULLIFY(interfaceEquations%interfaceMapping%lagrangeVariable)
    interfaceEquations%interfaceMapping%numberOfColumns=0
    interfaceEquations%interfaceMapping%totalNumberOfColumns=0
    interfaceEquations%interfaceMapping%numberOfGlobalColumns=0
    NULLIFY(interfaceEquations%interfaceMapping%columnDOFSMapping)
    NULLIFY(interfaceEquations%interfaceMapping%lagrangeDOFToColumnMap)
    interfaceEquations%interfaceMapping%numberOfInterfaceMatrices=0
    NULLIFY(interfaceEquations%interfaceMapping%rhsMapping)
    NULLIFY(interfaceEquations%interfaceMapping%createValuesCache)
    CALL InterfaceMapping_CreateValuesCacheInitialise(interfaceEquations%interfaceMapping,err,error,*999)
        
    EXITS("InterfaceMapping_Initialise")
    RETURN
999 CALL InterfaceMapping_Finalise(interfaceEquations%interfaceMapping,dummyErr,dummyError,*998)
998 ERRORSEXITS("InterfaceMapping_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMapping_Initialise

  !
  !================================================================================================================================
  !

  !>Sets the Lagrange variable type for the interface mapping. 
  SUBROUTINE InterfaceMapping_LagrangeVariableTypeSet(interfaceMapping,lagrangeVariableType,err,error,*)

    !Argument variables
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping !<A pointer to the interface mapping.
    INTEGER(INTG), INTENT(IN) :: lagrangeVariableType !<The Lagrange variable type to set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: interfaceConditionMethod
    TYPE(FieldType), POINTER :: lagrangeField
    TYPE(FieldVariableType), POINTER :: lagrangeVariable
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations
    TYPE(InterfaceLagrangeType), POINTER :: interfaceLagrange
    TYPE(InterfaceMappingCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(VARYING_STRING) :: localError

    ENTERS("InterfaceMapping_LagrangeVariableTypeSet",err,error,*999)

    CALL InterfaceMapping_AssertNotFinished(interfaceMapping,err,error,*999)
    
    NULLIFY(createValuesCache)
    CALL InterfaceMapping_CreateValuesCacheGet(interfaceMapping,createValuesCache,err,error,*999)
    NULLIFY(interfaceEquations)
    CALL InterfaceMapping_InterfaceEquationsGet(interfaceMapping,interfaceEquations,err,error,*999)
    NULLIFY(interfaceCondition)
    CALL InterfaceEquations_InterfaceConditionGet(interfaceEquations,interfaceCondition,err,error,*999)
    CALL InterfaceCondition_MethodGet(interfaceCondition,interfaceConditionMethod,err,error,*999)
    SELECT CASE(interfaceConditionMethod)
    CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
      NULLIFY(interfaceLagrange)
      CALL InterfaceCondition_InterfaceLagrangeGet(interfaceCondition,interfaceLagrange,err,error,*999)
      CALL InterfaceLagrange_AssertIsFinished(interfaceLagrange,err,error,*999)
      NULLIFY(lagrangeField)
      CALL InterfaceCondition_LagrangeFieldGet(interfaceCondition,lagrangeField,err,error,*999)
      NULLIFY(lagrangeVariable)
      CALL Field_VariableGet(lagrangeField,lagrangeVariableType,lagrangeVariable,err,error,*999)
      createValuesCache%lagrangeVariableType=lagrangeVariableType
    CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The interface condition method of "//TRIM(NumberToVString(interfaceConditionMethod,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("InterfaceMapping_LagrangeVariableTypeSet")
    RETURN
999 ERRORSEXITS("InterfaceMapping_LagrangeVariableTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMapping_LagrangeVariableTypeSet

  !
  !================================================================================================================================
  !

  !>Finalises an interface matrix to variable map and deallocates all memory.
  SUBROUTINE InterfaceMapping_MatrixToVarMapFinalise(interfaceMatrixToVarMap,err,error,*)

    !Argument variables
    TYPE(InterfaceMatrixToVarMapType), POINTER :: interfaceMatrixToVarMap !<The interface matrix to var map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
  
    ENTERS("InterfaceMapping_MatrixToVarMapFinalise",err,error,*999)

    IF(ASSOCIATED(interfaceMatrixToVarMap)) THEN
      IF(ASSOCIATED(interfaceMatrixToVarMap%variableDOFToRowMap)) DEALLOCATE(interfaceMatrixToVarMap%variableDOFToRowMap)
      DEALLOCATE(interfaceMatrixToVarMap)
    ENDIF
      
    EXITS("InterfaceMapping_MatrixToVarMapFinalise")
    RETURN
999 ERRORSEXITS("InterfaceMapping_MatrixToVarMapFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMapping_MatrixToVarMapFinalise

  !
  !================================================================================================================================
  !

  !>Initialises an interface matrix to variable map.
  SUBROUTINE InterfaceMapping_MatrixToVarMapInitialise(interfaceMapping,matrixIdx,err,error,*)

    !Argument variables
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping !<A pointer to the interface mapping to initialise the matrix to variable map for a given matrix index.
    INTEGER(INTG), INTENT(IN) :: matrixIdx !<The matrix index to intialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
  
    ENTERS("InterfaceMapping_MatrixToVarMapInitialise",err,error,*999)

    IF(.NOT.ASSOCIATED(interfaceMapping)) CALL FlagError("Interface mapping is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(interfaceMapping%interfaceMatrixToVarMaps)) &
      & CALL FlagError("Interface mapping matrix rows to var maps is not allocated.",err,error,*999)
    IF(matrixIdx<=0.OR.matrixIdx>interfaceMapping%numberOfInterfaceMatrices) THEN
      localError="The specified matrix index of "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
        & " is invalid. The index must be >= 1 and <= "// &
        & TRIM(NumberToVString(interfaceMapping%numberOfInterfaceMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(ASSOCIATED(interfaceMapping%interfaceMatrixToVarMaps(matrixIdx)%ptr)) THEN
      localError="The interface matrix rows to variable map is already associated for matrix index "// &
        & TRIM(NumberToVString(matrixIdx,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    ALLOCATE(interfaceMapping%interfaceMatrixToVarMaps(matrixIdx)%ptr,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate interface matrix rows to variable map.",err,error,*999)
    interfaceMapping%interfaceMatrixToVarMaps(matrixIdx)%ptr%matrixNumber=matrixIdx
    NULLIFY(interfaceMapping%interfaceMatrixToVarMaps(matrixIdx)%ptr%interfaceMatrix)
    NULLIFY(interfaceMapping%interfaceMatrixToVarMaps(matrixIdx)%ptr%equationsSet)
    interfaceMapping%interfaceMatrixToVarMaps(matrixIdx)%ptr%variableType=0
    NULLIFY(interfaceMapping%interfaceMatrixToVarMaps(matrixIdx)%ptr%variable)
    interfaceMapping%interfaceMatrixToVarMaps(matrixIdx)%ptr%meshIndex=0
    interfaceMapping%interfaceMatrixToVarMaps(matrixIdx)%ptr%matrixCoefficient=0.0_DP
    interfaceMapping%interfaceMatrixToVarMaps(matrixIdx)%ptr%hasTranspose=.FALSE.
    interfaceMapping%interfaceMatrixToVarMaps(matrixIdx)%ptr%numberOfRows=0
    interfaceMapping%interfaceMatrixToVarMaps(matrixIdx)%ptr%totalNumberOfRows=0
    interfaceMapping%interfaceMatrixToVarMaps(matrixIdx)%ptr%numberOfGlobalRows=0
    NULLIFY(interfaceMapping%interfaceMatrixToVarMaps(matrixIdx)%ptr%rowDOFsMapping)          
    NULLIFY(interfaceMapping%interfaceMatrixToVarMaps(matrixIdx)%ptr%variableDOFToRowMap)          
    
    EXITS("InterfaceMapping_MatrixToVarMapInitialise")
    RETURN
999 ERRORSEXITS("InterfaceMapping_MatrixToVarMapInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMapping_MatrixToVarMapInitialise

  !
  !================================================================================================================================
  !

  !>Sets the coefficients for the interface matrices. 
  SUBROUTINE InterfaceMapping_MatrixCoefficientsSet(interfaceMapping,matrixCoefficients,err,error,*)

    !Argument variables
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping !<A pointer to the interface mapping.
    REAL(DP), INTENT(IN) :: matrixCoefficients(:) !<The interface matrix coefficients
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: interfaceConditionMethod
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations
    TYPE(InterfaceMappingCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(VARYING_STRING) :: localError

    ENTERS("InterfaceMapping_MatrixCoefficientsSet",err,error,*999)

    CALL InterfaceMapping_AssertNotFinished(interfaceMapping,err,error,*999)

    NULLIFY(createValuesCache)
    CALL InterfaceMapping_CreateValuesCacheGet(interfaceMapping,createValuesCache,err,error,*999)
    NULLIFY(interfaceEquations)
    CALL InterfaceMapping_InterfaceEquationsGet(interfaceMapping,interfaceEquations,err,error,*999)
    NULLIFY(interfaceCondition)
    CALL InterfaceEquations_InterfaceConditionGet(interfaceEquations,interfaceCondition,err,error,*999)
    CALL InterfaceCondition_MethodGet(interfaceCondition,interfaceConditionMethod,err,error,*999)
    SELECT CASE(interfaceConditionMethod)
    CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
      !Check that the number of supplied coefficients matches the number of interface matrices
      IF(SIZE(matrixCoefficients,1)<createValuesCache%numberOfInterfaceMatrices) THEN
        localError="The size of the supplied array of matrix coefficients of "// &
          & TRIM(NumberToVString(SIZE(matrixCoefficients,1),"*",err,error))// &
          & "is invalid. The size must be >="// &
          & TRIM(NumberToVString(createValuesCache%numberOfInterfaceMatrices,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      createValuesCache%matrixCoefficients(1:createValuesCache%numberOfInterfaceMatrices)= &
        & matrixCoefficients(1:createValuesCache%numberOfInterfaceMatrices)
    CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The interface condition method of "//TRIM(NumberToVString(interfaceConditionMethod,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
   
    EXITS("InterfaceMapping_MatrixCoefficientsSet")
    RETURN
999 ERRORSEXITS("InterfaceMapping_MatrixCoefficientsSet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMapping_MatrixCoefficientsSet

  !
  !================================================================================================================================
  !

  !>Sets the coefficients for the transpose interface matrices. 
  SUBROUTINE InterfaceMapping_TransposeMatrixCoefficientsSet(interfaceMapping,transposeMatrixCoefficients,err,error,*)

    !Argument variables
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping !<A pointer to the interface mapping.
    REAL(DP), INTENT(IN) :: transposeMatrixCoefficients(:) !<The transpose interface matrix coefficients
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: interfaceConditionMethod
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations
    TYPE(InterfaceMappingCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(VARYING_STRING) :: localError

    ENTERS("InterfaceMapping_TransposeMatrixCoefficientsSet",err,error,*999)

    CALL InterfaceMapping_AssertNotFinished(interfaceMapping,err,error,*999)

    NULLIFY(createValuesCache)
    CALL InterfaceMapping_CreateValuesCacheGet(interfaceMapping,createValuesCache,err,error,*999)
    NULLIFY(interfaceEquations)
    CALL InterfaceMapping_InterfaceEquationsGet(interfaceMapping,interfaceEquations,err,error,*999)
    NULLIFY(interfaceCondition)
    CALL InterfaceEquations_InterfaceConditionGet(interfaceEquations,interfaceCondition,err,error,*999)
    CALL InterfaceCondition_MethodGet(interfaceCondition,interfaceConditionMethod,err,error,*999)
    SELECT CASE(interfaceConditionMethod)
    CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
      !Check that the number of supplied coefficients matches the number of interface matrices
      IF(SIZE(transposeMatrixCoefficients,1)<createValuesCache%numberOfInterfaceMatrices) THEN
        localError="The size of the supplied transpose matrix coefficients array of "// &
          & TRIM(NumberToVString(SIZE(transposeMatrixCoefficients,1),"*",err,error))// &
          & "is invalid. The size should be >= "// &
          & TRIM(NumberToVString(createValuesCache%numberOfInterfaceMatrices,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      createValuesCache%transposeMatrixCoefficients(1:createValuesCache%numberOfInterfaceMatrices)= &
        & transposeMatrixCoefficients(1:createValuesCache%numberOfInterfaceMatrices)
    CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The interface condition method of "//TRIM(NumberToVString(interfaceConditionMethod,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
   
    EXITS("InterfaceMapping_TransposeMatrixCoefficientsSet")
    RETURN
999 ERRORSEXITS("InterfaceMapping_TransposeMatrixCoefficientsSet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMapping_TransposeMatrixCoefficientsSet

  !
  !================================================================================================================================
  !

  !>Sets the column mesh indices for the interface matrices. 
  SUBROUTINE InterfaceMapping_MatricesColumnMeshIndicesSet(interfaceMapping,columnMeshIndices,err,error,*)

    !Argument variables
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping !<A pointer to the interface mapping.
    INTEGER(INTG), INTENT(IN) :: columnMeshIndices(:) !<The interface matrix column mesh indices
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: interfaceConditionMethod
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations
    TYPE(InterfaceMappingCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(VARYING_STRING) :: localError

    ENTERS("InterfaceMapping_MatricesColumnMeshIndicesSet",err,error,*999)

    CALL InterfaceMapping_AssertNotFinished(interfaceMapping,err,error,*999)
    
    NULLIFY(createValuesCache)
    CALL InterfaceMapping_CreateValuesCacheGet(interfaceMapping,createValuesCache,err,error,*999)
    NULLIFY(interfaceEquations)
    CALL InterfaceMapping_InterfaceEquationsGet(interfaceMapping,interfaceEquations,err,error,*999)
    NULLIFY(interfaceCondition)
    CALL InterfaceEquations_InterfaceConditionGet(interfaceEquations,interfaceCondition,err,error,*999)
    CALL InterfaceCondition_MethodGet(interfaceCondition,interfaceConditionMethod,err,error,*999)
    SELECT CASE(interfaceConditionMethod)
    CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
      CALL FlagError("Can not set the column mesh indices when using the Lagrange multipliers interface condition method.", &
        & err,error,*999)
    CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(INTERFACE_CONDITION_PENALTY_METHOD)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The interface condition method of "//TRIM(NumberToVString(interfaceConditionMethod,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("InterfaceMapping_MatricesColumnMeshIndicesSet")
    RETURN
999 ERRORS("InterfaceMapping_MatricesColumnMeshIndicesSet",err,error)
    EXITS("InterfaceMapping_MatricesColumnMeshIndicesSet")
    RETURN 1
    
  END SUBROUTINE InterfaceMapping_MatricesColumnMeshIndicesSet

  !
  !================================================================================================================================
  !

  !>Sets the row mesh indices for the interface matrices. 
  SUBROUTINE InterfaceMapping_MatricesRowMeshIndicesSet(interfaceMapping,rowMeshIndices,err,error,*)

    !Argument variables
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping !<A pointer to the interface mapping.
    INTEGER(INTG), INTENT(IN) :: rowMeshIndices(:) !<The interface matrix mesh indices
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: interfaceConditionMethod,meshIdx,meshIdx2,meshIdx3
    LOGICAL :: found
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition
    TYPE(InterfaceDependentType), POINTER :: interfaceDependent
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations
    TYPE(InterfaceMappingCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(VARYING_STRING) :: localError

    ENTERS("InterfaceMapping_MatricesRowMeshIndicesSet",err,error,*999)

    CALL InterfaceMapping_AssertNotFinished(interfaceMapping,err,error,*999)
    
    NULLIFY(createValuesCache)
    CALL InterfaceMapping_CreateValuesCacheGet(interfaceMapping,createValuesCache,err,error,*999)
    NULLIFY(interfaceEquations)
    CALL InterfaceMapping_InterfaceEquationsGet(interfaceMapping,interfaceEquations,err,error,*999)
    NULLIFY(interfaceCondition)
    CALL InterfaceEquations_InterfaceConditionGet(interfaceEquations,interfaceCondition,err,error,*999)
    CALL InterfaceCondition_MethodGet(interfaceCondition,interfaceConditionMethod,err,error,*999)
    SELECT CASE(interfaceConditionMethod)
    CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
      !Check the size of the mesh indicies array
      IF(SIZE(rowMeshIndices,1)/=createValuesCache%numberOfInterfaceMatrices) THEN
        localError="Invalid size of row mesh indices. The size of the supplied array of "// &
          & TRIM(NumberToVString(SIZE(rowMeshIndices,1),"*",err,error))// &
          & " must match the number of interface matrices of "// &
          & TRIM(NumberToVString(createValuesCache%numberOfInterfaceMatrices,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      !Check that mesh indices are valid.
      NULLIFY(interfaceDependent)
      CALL InterfaceCondition_InterfaceDependentGet(interfaceCondition,interfaceDependent,err,error,*999)
      DO meshIdx=1,createValuesCache%numberOfInterfaceMatrices
        found=.FALSE.
        DO meshIdx2=1,interfaceDependent%numberOfDependentVariables
          IF(rowMeshIndices(meshIdx)==interfaceDependent%variableMeshIndices(meshIdx2)) THEN
            found=.TRUE.
            EXIT
          ENDIF
        ENDDO !meshIdx2
        IF(.NOT.found) THEN
          localError="The supplied mesh index of "//TRIM(NumberToVString(rowMeshIndices(meshIdx),"*",err,error))// &
            & " at position "//TRIM(NumberToVString(meshIdx,"*",err,error))// &
            & " has not been added as a dependent variable to the interface condition."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        !Check that the mesh index has not been repeated.
        DO meshIdx3=meshIdx+1,createValuesCache%numberOfInterfaceMatrices
          IF(rowMeshIndices(meshIdx)==rowMeshIndices(meshIdx3)) THEN
            localError="The supplied mesh index of "//TRIM(NumberToVString(rowMeshIndices(meshIdx),"*",err,error))// &
              & " at position "//TRIM(NumberToVString(meshIdx,"*",err,error))// &
              & " has been repeated at position "//TRIM(NumberToVString(meshIdx3,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ENDDO !meshIdx3
        !Set the mesh indices
        createValuesCache%matrixRowFieldVariableIndices(1:createValuesCache%numberOfInterfaceMatrices)= &
          & rowMeshIndices(1:createValuesCache%numberOfInterfaceMatrices)
      ENDDO !meshIdx
    CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The interface condition method of "//TRIM(NumberToVString(interfaceConditionMethod,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
   
    EXITS("InterfaceMapping_MatricesRowMeshIndicesSet")
    RETURN
999 ERRORSEXITS("InterfaceMapping_MatricesRowMeshIndicesSet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMapping_MatricesRowMeshIndicesSet

  !
  !================================================================================================================================
  !

  !>Sets the number of interface matrices for an interface mapping.
  SUBROUTINE InterfaceMapping_NumberOfMatricesSet(interfaceMapping,numberOfInterfaceMatrices,err,error,*)

    !Argument variables
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping !<A pointer to the interface mapping to set the number of linear matrices for.
    INTEGER(INTG), INTENT(IN) :: numberOfInterfaceMatrices !<The number of interface matrices to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: interfaceConditionMethod,matrixIdx,matrixIdx2,meshIdx,numberOfDependentVariables,variableIdx
    INTEGER(INTG), ALLOCATABLE :: newMatrixRowFieldVariableIndices(:)
    REAL(DP), ALLOCATABLE :: newMatrixCoefficients(:),newTransposeMatrixCoefficients(:)
    LOGICAL :: found
    LOGICAL, ALLOCATABLE :: newHasTranspose(:)
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition
    TYPE(InterfaceDependentType), POINTER :: interfaceDependent
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations
    TYPE(InterfaceMappingCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(VARYING_STRING) :: localError

    ENTERS("InterfaceMapping_NumberOfMatricesSet",err,error,*999)

    CALL InterfaceMapping_AssertNotFinished(interfaceMapping,err,error,*999)

    NULLIFY(createValuesCache)
    CALL InterfaceMapping_CreateValuesCacheGet(interfaceMapping,createValuesCache,err,error,*999)
    NULLIFY(interfaceEquations)
    CALL InterfaceMapping_InterfaceEquationsGet(interfaceMapping,interfaceEquations,err,error,*999)
    NULLIFY(interfaceCondition)
    CALL InterfaceEquations_InterfaceConditionGet(interfaceEquations,interfaceCondition,err,error,*999)
    CALL InterfaceCondition_MethodGet(interfaceCondition,interfaceConditionMethod,err,error,*999)
    SELECT CASE(interfaceConditionMethod)
    CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
      !Check the number of interface matrices
      IF(numberOfInterfaceMatrices<=0) THEN
        localError="The specified number of interface matrices of "// &
          & TRIM(NumberToVString(numberOfInterfaceMatrices,"*",err,error))//" is invalid. The number must be >= 1."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      NULLIFY(interfaceDependent)
      CALL InterfaceCondition_InterfaceDependentGet(interfaceCondition,interfaceDependent,err,error,*999)
      SELECT CASE(interfaceConditionMethod)
      CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
        numberOfDependentVariables=interfaceDependent%numberOfDependentVariables
      CASE(INTERFACE_CONDITION_PENALTY_METHOD)
        numberOfDependentVariables=interfaceDependent%numberOfDependentVariables+1
      END SELECT
      IF(numberOfInterfaceMatrices>numberOfDependentVariables) THEN
        localError="The specified number of interface matrices of "// &
          & TRIM(NumberToVString(numberOfInterfaceMatrices,"*",err,error))// &
          & " is invalid. The number must be <= the number of added dependent variables of "// &
          & TRIM(NumberToVString(numberOfDependentVariables,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      !If we need to reallocate and reset all the create values cache arrays and change the number of matrices
      IF(numberOfInterfaceMatrices/=createValuesCache%numberOfInterfaceMatrices) THEN
        ALLOCATE(newMatrixCoefficients(createValuesCache%numberOfInterfaceMatrices),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate new matrix coefficients.",err,error,*999)
        ALLOCATE(newTransposeMatrixCoefficients(createValuesCache%numberOfInterfaceMatrices),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate new transpose matrix coefficients.",err,error,*999)
        ALLOCATE(newHasTranspose(createValuesCache%numberOfInterfaceMatrices),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate new matrix transpose.",err,error,*999)
        ALLOCATE(newMatrixRowFieldVariableIndices(createValuesCache%numberOfInterfaceMatrices),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate new matrix row field indexes.",err,error,*999)
        newMatrixCoefficients(1:createValuesCache%numberOfInterfaceMatrices)= &
          createValuesCache%matrixCoefficients(1:createValuesCache%numberOfInterfaceMatrices)
        newTransposeMatrixCoefficients(1:createValuesCache%numberOfInterfaceMatrices)= &
          createValuesCache%transposeMatrixCoefficients(1:createValuesCache%numberOfInterfaceMatrices)
        newHasTranspose(1:createValuesCache%numberOfInterfaceMatrices)= &
          & createValuesCache%hasTranspose(1:createValuesCache%numberOfInterfaceMatrices)
        newMatrixRowFieldVariableIndices(1:createValuesCache%numberOfInterfaceMatrices)= &
          & createValuesCache%matrixRowFieldVariableIndices(1:createValuesCache%numberOfInterfaceMatrices)
        IF(numberOfInterfaceMatrices>createValuesCache%numberOfInterfaceMatrices) THEN
          !Loop through in mesh index order and set the default matrix to variable map to be in mesh index order
          DO matrixIdx=createValuesCache%numberOfInterfaceMatrices+1,numberOfInterfaceMatrices
            createValuesCache%matrixRowFieldVariableIndices(matrixIdx)=0
            DO variableIdx=1,interfaceDependent%numberOfDependentVariables
              found=.FALSE.
              DO matrixIdx2=1,createValuesCache%numberOfInterfaceMatrices
                CALL InterfaceMappingCVC_MatrixRowVariableIndexGet(createValuesCache,matrixIdx2,meshIdx,err,error,*999)
                IF(interfaceDependent%variableMeshIndices(variableIdx)==meshIdx) THEN
                  found=.TRUE.
                  EXIT
                ENDIF
              ENDDO !matrixIdx2
              IF(.NOT.found) newMatrixRowFieldVariableIndices(matrixIdx)=interfaceDependent%variableMeshIndices(variableIdx)
            ENDDO !variableIdx2
            IF(newMatrixRowFieldVariableIndices(matrixIdx)==0) THEN
              localError="Could not map an interface mesh index for interface matrix "// &
                & TRIM(NumberToVString(matrixIdx,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ENDDO !matrixIdx
        ENDIF
        CALL MOVE_ALLOC(newMatrixCoefficients,createValuesCache%matrixCoefficients)
        CALL MOVE_ALLOC(newTransposeMatrixCoefficients,createValuesCache%transposeMatrixCoefficients)
        CALL MOVE_ALLOC(newHasTranspose,createValuesCache%hasTranspose)
        CALL MOVE_ALLOC(newMatrixRowFieldVariableIndices,createValuesCache%matrixRowFieldVariableIndices)
        createValuesCache%numberOfInterfaceMatrices=numberOfInterfaceMatrices
      ENDIF
    CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The interface condition method of "// &
        & TRIM(NumberToVString(interfaceConditionMethod,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("InterfaceMapping_NumberOfMatricesSet")
    RETURN
999 IF(ALLOCATED(newMatrixCoefficients)) DEALLOCATE(newMatrixCoefficients)
    IF(ALLOCATED(newTransposeMatrixCoefficients)) DEALLOCATE(newTransposeMatrixCoefficients)
    IF(ALLOCATED(newHastranspose)) DEALLOCATE(newHastranspose)
    IF(ALLOCATED(newMatrixRowFieldVariableIndices)) DEALLOCATE(newMatrixRowFieldVariableIndices)
    ERRORSEXITS("InterfaceMapping_NumberOfMatricesSet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMapping_NumberOfMatricesSet

  !
  !================================================================================================================================
  !

  !>Sets the transpose flag for the interface matrices. 
  SUBROUTINE InterfaceMapping_MatricesTransposeSet(interfaceMapping,matrixTranspose,err,error,*)

    !Argument variables
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping !<A pointer to the interface mapping.
    LOGICAL, INTENT(IN) :: matrixTranspose(:) !<matrixTranspose(matrixIdx). The interface matrix transpose flag for the matrixIdx'th interface matrix.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: interfaceConditionMethod
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations
    TYPE(InterfaceMappingCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(VARYING_STRING) :: localError

    ENTERS("InterfaceMapping_MatricesTransposeSet",err,error,*999)

    CALL InterfaceMapping_AssertNotFinished(interfaceMapping,err,error,*999)
    
    NULLIFY(createValuesCache)
    CALL InterfaceMapping_CreateValuesCacheGet(interfaceMapping,createValuesCache,err,error,*999)
    NULLIFY(interfaceEquations)
    CALL InterfaceMapping_InterfaceEquationsGet(interfaceMapping,interfaceEquations,err,error,*999)
    NULLIFY(interfaceCondition)
    CALL InterfaceEquations_InterfaceConditionGet(interfaceEquations,interfaceCondition,err,error,*999)
    CALL InterfaceCondition_MethodGet(interfaceCondition,interfaceConditionMethod,err,error,*999)
    SELECT CASE(interfaceConditionMethod)
    CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
      !Check that the number of supplied coefficients matches the number of interface matrices
      IF(SIZE(matrixTranspose,1)<createValuesCache%numberOfInterfaceMatrices) THEN
        localError="The size of the supplied matrix transpose array of "// &
          & TRIM(NumberToVString(SIZE(matrixTranspose,1),"*",err,error))// &
          & "is invalid. The size should be >= "// &
          & TRIM(NumberToVString(createValuesCache%numberOfInterfaceMatrices,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      createValuesCache%hasTranspose(1:createValuesCache%numberOfInterfaceMatrices)= &
        matrixTranspose(1:createValuesCache%numberOfInterfaceMatrices)
    CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The interface condition method of "//TRIM(NumberToVString(interfaceConditionMethod,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("InterfaceMapping_MatricesTransposeSet")
    RETURN
999 ERRORSEXITS("InterfaceMapping_MatricesTransposeSet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMapping_MatricesTransposeSet

  !
  !================================================================================================================================
  !

  !>Sets the coefficient applied to the interface RHS vector.
  SUBROUTINE InterfaceMapping_RHSCoefficientSet(interfaceMapping,rhsCoefficient,err,error,*)

    !Argument variables
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping !<A pointer to the interface mapping to set the RHS coefficent for
    REAL(DP), INTENT(IN) :: rhsCoefficient !<The coefficient applied to the interface RHS vector.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(InterfaceMappingCreateValuesCacheType), POINTER :: createValuesCache

    ENTERS("InterfaceMapping_RHSCoefficientSet",err,error,*999)

    CALL InterfaceMapping_AssertNotFinished(interfaceMapping,err,error,*999)
   
    NULLIFY(createValuesCache)
    CALL InterfaceMapping_CreateValuesCacheGet(interfaceMapping,createValuesCache,err,error,*999)
    IF(createValuesCache%rhsLagrangeVariableType==0) &
      & CALL FlagError("The interface mapping RHS Lagrange variable type has not been set.",err,error,*999)
    
    createValuesCache%rhsCoefficient=rhsCoefficient
       
    EXITS("InterfaceMapping_RHSCoefficientSet")
    RETURN
999 ERRORSEXITS("InterfaceMapping_RHSCoefficientSet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMapping_RHSCoefficientSet
  
  !
  !================================================================================================================================
  !

  !>Finalises the interface mapping RHS mapping and deallocates all memory
  SUBROUTINE InterfaceMapping_RHSMappingFinalise(rhsMapping,err,error,*)

    !Argument variables
    TYPE(InterfaceMappingRHSType), POINTER :: rhsMapping !<A pointer to the RHS mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMapping_RHSMappingFinalise",err,error,*999)

    IF(ASSOCIATED(rhsMapping)) THEN
      IF(ASSOCIATED(rhsMapping%rhsDOFToInterfaceRowMap)) DEALLOCATE(rhsMapping%rhsDOFToInterfaceRowMap)
      IF(ASSOCIATED(rhsMapping%interfaceRowToRHSDOFMap)) DEALLOCATE(rhsMapping%interfaceRowToRHSDOFMap)
      DEALLOCATE(rhsMapping)
    ENDIF
       
    EXITS("InterfaceMapping_RHSMappingFinalise")
    RETURN
999 ERRORSEXITS("InterfaceMapping_RHSMappingFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMapping_RHSMappingFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the interface mapping RHS mapping
  SUBROUTINE InterfaceMapping_RHSMappingInitialise(interfaceMapping,err,error,*)

    !Argument variables
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping !<A pointer to the interface mapping to initialise the RHS mapping for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("InterfaceMapping_RHSMappingInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(interfaceMapping)) CALL FlagError("Interface mapping is not associated.",err,error,*998)
    IF(ASSOCIATED(interfaceMapping%rhsMapping)) &
      & CALL FlagError("Interface mapping RHS mapping is already associated.",err,error,*998)
    
    ALLOCATE(interfaceMapping%rhsMapping,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate interface mapping RHS mapping.",err,error,*999)
    interfaceMapping%rhsMapping%interfaceMapping=>interfaceMapping      
    interfaceMapping%rhsMapping%rhsVariableType=0
    NULLIFY(interfaceMapping%rhsMapping%rhsVariable)
    NULLIFY(interfaceMapping%rhsMapping%rhsVariableMapping)
    interfaceMapping%rhsMapping%rhsCoefficient=1.0_DP
    NULLIFY(interfaceMapping%rhsMapping%rhsDOFToInterfaceRowMap)
    NULLIFY(interfaceMapping%rhsMapping%interfaceRowToRHSDOFMap)
       
    EXITS("InterfaceMapping_RHSMappingInitialise")
    RETURN
999 CALL InterfaceMapping_RHSMappingFinalise(interfaceMapping%rhsMapping,dummyErr,dummyError,*998)
998 ERRORSEXITS("InterfaceMapping_RHSMappingInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMapping_RHSMappingInitialise

  !
  !================================================================================================================================
  !

  !>Sets the mapping between a Lagrange field variable and the interface rhs vector.
  SUBROUTINE InterfaceMapping_RHSVariableTypeSet(interfaceMapping,rhsVariableType,err,error,*)

    !Argument variables
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping !<A pointer to the interface mapping to set the RHS variable type for.
    INTEGER(INTG), INTENT(IN) :: rhsVariableType !<The variable type associated with the interface rhs vector. If the interface condition does not have a rhs vector then the variable type on input should be zero.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: interfaceConditionMethod
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations
    TYPE(InterfaceMappingCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(FieldType), POINTER :: lagrangeField
    TYPE(FieldVariableType), POINTER :: lagrangeVariable,rhsVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("InterfaceMapping_RHSVariableTypeSet",err,error,*999)

    CALL InterfaceMapping_AssertNotFinished(interfaceMapping,err,error,*999)

    NULLIFY(createValuesCache)
    CALL InterfaceMapping_CreateValuesCacheGet(interfaceMapping,createValuesCache,err,error,*999)
    IF(rhsVariableType==0) THEN
      createValuesCache%rhsLagrangeVariableType=0
    ELSE
      NULLIFY(interfaceEquations)
      CALL InterfaceMapping_InterfaceEquationsGet(interfaceMapping,interfaceEquations,err,error,*999)
      NULLIFY(interfaceCondition)
      CALL InterfaceEquations_InterfaceConditionGet(interfaceEquations,interfaceCondition,err,error,*999)
      CALL InterfaceCondition_MethodGet(interfaceCondition,interfaceConditionMethod,err,error,*999)
      SELECT CASE(interfaceConditionMethod)
      CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
        NULLIFY(lagrangeField)
        CALL InterfaceCondition_LagrangeFieldGet(interfaceCondition,lagrangeField,err,error,*999)
        !Check the RHS variable type is not being by the interface matrices
        IF(createValuesCache%lagrangeVariableType==rhsVariableType) THEN
          localError="The specified RHS variable type of "//TRIM(NumberToVString(rhsVariableType,"*",err,error))// &
            & " is the same as the Lagrange variable type for the interface matrices."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        !Check the RHS variable number is defined on the Lagrange field
        NULLIFY(rhsVariable)
        CALL Field_VariableGet(lagrangeField,rhsVariableType,rhsVariable,err,error,*999)
        createValuesCache%rhsLagrangeVariableType=rhsVariableType
      CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The interface condition method of "//TRIM(NumberToVString(interfaceConditionMethod,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDIF
       
    EXITS("InterfaceMapping_RHSVariableTypeSet")
    RETURN
999 ERRORSEXITS("InterfaceMapping_RHSVariableTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMapping_RHSVariableTypeSet
  
  !
  !================================================================================================================================
  !

END MODULE InterfaceMappingRoutines
