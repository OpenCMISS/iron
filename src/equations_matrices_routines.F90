!> \file
!> \author Chris Bradley
!> \brief This module handles all equations matrix and rhs routines.
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
!> Contributor(s): David Ladd
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

!> This module handles all equations matrix and rhs routines.
MODULE EquationsMatricesRoutines

  USE BaseRoutines
  USE DistributedMatrixVector
  USE DistributedMatrixVectorAccessRoutines
  USE EquationsAccessRoutines
  USE EquationsMappingAccessRoutines
  USE EquationsMatricesAccessRoutines
  USE EquationsSetAccessRoutines
  USE FieldRoutines
  USE FieldAccessRoutines
  USE ISO_VARYING_STRING
  USE Kinds
  USE Lists
  USE MatrixVector
  USE Strings
  USE Types
  USE LINKEDLIST_ROUTINES

#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters


  !Module types

  !Module variables

  !Interfaces

  !>Sets the storage type (sparsity) of the Linear equations matrices
  INTERFACE EquationsMatricesVector_LinearStorageTypeSet
    MODULE PROCEDURE EquationsMatricesVector_LinearStorageTypeSet0
    MODULE PROCEDURE EquationsMatricesVector_LinearStorageTypeSet1
  END INTERFACE EquationsMatricesVector_LinearStorageTypeSet

  !>Sets the structure (sparsity) of the Linear equations matrices
  INTERFACE EquationsMatricesVector_LinearStructureTypeSet
    MODULE PROCEDURE EquationsMatricesVector_LinearStructureTypeSet0
    MODULE PROCEDURE EquationsMatricesVector_LinearStructureTypeSet1
  END INTERFACE EquationsMatricesVector_LinearStructureTypeSet

  !>Sets the storage type (sparsity) of the nonlinear (Jacobian) equations matrices
  INTERFACE EquationsMatricesVector_NonlinearStorageTypeSet
    MODULE PROCEDURE EquationsMatricesVector_NonlinearStorageTypeSet0
    MODULE PROCEDURE EquationsMatricesVector_NonlinearStorageTypeSet1
  END INTERFACE EquationsMatricesVector_NonlinearStorageTypeSet

  !>Sets the structure (sparsity) of the nonlinear (Jacobian) equations matrices
  INTERFACE EquationsMatricesVector_NonlinearStructureTypeSet
    MODULE PROCEDURE EquationsMatricesVector_NonlinearStructureTypeSet0
    MODULE PROCEDURE EquationsMatricesVector_NonlinearStructureTypeSet1
  END INTERFACE EquationsMatricesVector_NonlinearStructureTypeSet

  PUBLIC EquationsMatricesVector_DynamicLumpingTypeSet

  PUBLIC EquationsMatricesVector_DynamicStorageTypeSet

  PUBLIC EquationsMatricesVector_DynamicStructureTypeSet

  !!TODO check if the elements should be create/destroy rather than initialise/finalise
  PUBLIC EquationsMatricesVector_ElementInitialise,EquationsMatricesVector_ElementFinalise

  PUBLIC EquationsMatricesVector_ElementAdd

  PUBLIC EquationsMatricesVector_JacobianElementAdd

  PUBLIC EquationsMatricesVector_ElementCalculate

  PUBLIC EquationsMatrices_ElementMatrixFinalise,EquationsMatrices_ElementMatrixInitialise
  
  PUBLIC EquationsMatrices_ElementMatrixCalculate

  PUBLIC EquationsMatrices_ElementMatrixSetup

  PUBLIC EquationsMatrices_ElementVectorFinalise,EquationsMatrices_ElementVectorInitialise

  PUBLIC EquationsMatrices_ElementVectorCalculate

  PUBLIC EquationsMatrices_ElementVectorSetup

  PUBLIC EquationsMatricesResidual_DistributedVectorEnsureCreated

  PUBLIC EquationsMatricesRHS_DistributedVectorEnsureCreated

  PUBLIC EquationsMatricesSource_DistributedVectorEnsureCreated

  PUBLIC EquationsMatricesVector_HessianOutput

  PUBLIC EquationsMatricesVector_JacobianNodeAdd

  PUBLIC EquationsMatricesVector_JacobianOutput

  PUBLIC EquationsMatricesVector_JacobianCalculationTypeSet

  PUBLIC EquationsMatricesVector_JacobianFiniteDifferenceStepSizeSet

  PUBLIC EquationsMatricesVector_LinearStorageTypeSet

  PUBLIC EquationsMatricesVector_LinearStructureTypeSet

  PUBLIC EquationsMatricesVector_NodalInitialise,EquationsMatricesVector_NodalFinalise

  PUBLIC EquationsMatricesVector_NodeAdd

  PUBLIC EquationsMatricesVector_NodalCalculate

  PUBLIC EquationsMatrices_NodalMatrixFinalise,EquationsMatrices_NodalMatrixInitialise

  PUBLIC EquationsMatrices_NodalMatrixCalculate

  PUBLIC EquationsMatrices_NodalMatrixSetup

  PUBLIC EquationsMatrices_NodalVectorFinalise,EquationsMatrices_NodalVectorInitialise

  PUBLIC EquationsMatrices_NodalVectorCalculate

  PUBLIC EquationsMatrices_NodalVectorSetup

  PUBLIC EquationsMatricesVector_NonlinearStorageTypeSet

  PUBLIC EquationsMatricesVector_NonlinearStructureTypeSet

  PUBLIC EquationsMatrices_ScalarDestroy

  PUBLIC EquationsMatrices_VectorCreateFinish,EquationsMatrices_VectorCreateStart

  PUBLIC EquationsMatrices_VectorDestroy

  PUBLIC EquationsMatricesVector_Output

  PUBLIC EquationsMatricesVector_VectorValuesInitialise

 CONTAINS

  !
  !================================================================================================================================
  !

  !>Finalise the equations Hessian and deallocate all memory
  SUBROUTINE EquationsMatrices_HessianFinalise(hessianMatrix,err,error,*)

    !Argument variables
    TYPE(HessianMatrixType), POINTER :: hessianMatrix !<A pointer to the equations Hessian to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("EquationsMatrices_HessianFinalise",err,error,*999)

    IF(.NOT.ASSOCIATED(hessianMatrix)) THEN
      IF(ASSOCIATED(hessianMatrix%hessian)) CALL DistributedMatrix_Destroy(hessianMatrix%hessian,err,error,*999)
      CALL EquationsMatrices_ElementMatrixFinalise(hessianMatrix%elementHessian,err,error,*999)
      DEALLOCATE(hessianMatrix)
    ENDIF
    
    EXITS("EquationsMatrices_HessianFinalise")
    RETURN
999 ERRORSEXITS("EquationsMatrices_HessianFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_HessianFinalise

  !
  !================================================================================================================================
  !

  !>Finishes the creation of the vector equations matrices and RHS for the vector equations
  SUBROUTINE EquationsMatrices_VectorCreateFinish(vectorMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<The pointer to the vector equations matrices
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string  
    !Local Variables
    INTEGER(INTG) :: dummyErr,matrixIdx,numberOfNonZeros,residualIdx,sourceIdx
    INTEGER(INTG), POINTER :: rowIndices(:),columnIndices(:)
    TYPE(DomainMappingType), POINTER :: rowDomainMap,columnDomainMap
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatricesSourcesType), POINTER :: sourceVectors
    TYPE(EquationsMatrixType), POINTER :: dynamicMatrix,linearMatrix
    TYPE(EquationsMatrixToVarMapType), POINTER :: equationsMatrixToVarMap
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix
    TYPE(JacobianMatrixToVarMapType), POINTER :: jacobianMatrixToVarMap
    TYPE(VARYING_STRING) :: dummyError,localError
    TYPE(LinkedList), POINTER :: list(:)
    
    NULLIFY(rowIndices)
    NULLIFY(columnIndices)

    ENTERS("EquationsMatrices_VectorCreateFinish",err,error,*998)

    CALL EquationsMatricesVector_AssertNotFinished(vectorMatrices,err,error,*999)
    NULLIFY(vectorMapping)
    CALL EquationsMatricesVector_VectorMappingGet(vectorMatrices,vectorMapping,err,error,*999)
    NULLIFY(lhsMapping)
    CALL EquationsMappingVector_LHSMappingGet(vectorMapping,lhsMapping,err,error,*999)
    NULLIFY(rowDomainMap)
    CALL EquationsMappingLHS_RowDOFsMappingGet(lhsMapping,rowDomainMap,err,error,*999)

    NULLIFY(dynamicMatrices)
    CALL EquationsMatricesVector_DynamicMatricesExists(vectorMatrices,dynamicMatrices,err,error,*999)
    IF(ASSOCIATED(dynamicMatrices)) THEN
      !Dynamic matrices
      NULLIFY(dynamicMapping)
      CALL EquationsMappingVector_DynamicMappingGet(vectorMapping,dynamicMapping,err,error,*999)
      !Now create the individual dynamic equations matrices
      DO matrixIdx=1,dynamicMatrices%numberOfDynamicMatrices
        NULLIFY(dynamicMatrix)
        CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,matrixIdx,dynamicMatrix,err,error,*999)
        NULLIFY(equationsMatrixToVarMap)
        CALL EquationsMappingDynamic_EquationsMatrixToVarMapGet(dynamicMapping,matrixIdx,equationsMatrixToVarMap,err,error,*999)
        NULLIFY(columnDomainMap)
        CALL EquationsMappingEMToVMap_ColumnDOFsMappingGet(equationsMatrixToVarMap,columnDomainMap,err,error,*999)
        !Create the distributed equations matrix
        NULLIFY(dynamicMatrix%matrix)
        CALL DistributedMatrix_CreateStart(rowDomainMap,columnDomainMap,dynamicMatrix%matrix,err,error,*999)
        CALL DistributedMatrix_DataTypeSet(dynamicMatrix%matrix,MATRIX_VECTOR_DP_TYPE,err,error,*999)
        CALL DistributedMatrix_StorageTypeSet(dynamicMatrix%matrix,dynamicMatrix%storageType,err,error,*999)
        CALL DistributedMatrix_TransposeTypeSet(dynamicMatrix%matrix,DISTRIBUTED_MATRIX_PARTIAL_TRANSPOSE_REQUIRED,err,error,*999)
        !Calculate and set the matrix structure/sparsity pattern
        IF(dynamicMatrix%storageType/=DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE.AND. &
          & dynamicMatrix%storageType/=DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE) THEN
          NULLIFY(list)
          CALL EquationsMatrix_StructureCalculate(dynamicMatrix,numberOfNonZeros,rowIndices,columnIndices,list,err,error,*999)
          CALL DistributedMatrix_LinkListSet(dynamicMatrix%matrix,list,err,error,*999)
          CALL DistributedMatrix_NumberOfNonZerosSet(dynamicMatrix%matrix,numberOfNonZeros,err,error,*999)
          CALL DistributedMatrix_StorageLocationsSet(dynamicMatrix%matrix,rowIndices,columnIndices,err,error,*999)
          IF(ASSOCIATED(rowIndices)) DEALLOCATE(rowIndices)
          IF(ASSOCIATED(columnIndices)) DEALLOCATE(columnIndices)
        ENDIF
        CALL DistributedMatrix_CreateFinish(dynamicMatrix%matrix,err,error,*999)
      ENDDO !matrixIdx                
    ENDIF
    NULLIFY(linearMatrices)
    CALL EquationsMatricesVector_LinearMatricesExists(vectorMatrices,linearMatrices,err,error,*999)
    IF(ASSOCIATED(linearMatrices)) THEN
      !Linear matrices
      NULLIFY(linearMapping)
      CALL EquationsMappingVector_LinearMappingGet(vectorMapping,linearMapping,err,error,*999)
      !Now create the individual linear equations matrices
      DO matrixIdx=1,linearMatrices%numberOfLinearMatrices
        NULLIFY(linearMatrix)
        CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,matrixIdx,linearMatrix,err,error,*999)
        NULLIFY(equationsMatrixToVarMap)
        CALL EquationsMappingLinear_EquationsMatrixToVarMapGet(linearMapping,matrixIdx,equationsMatrixToVarMap,err,error,*999)
        NULLIFY(columnDomainMap)
        CALL EquationsMappingEMToVMap_ColumnDOFsMappingGet(equationsMatrixToVarMap,columnDomainMap,err,error,*999)
        !Create the distributed equations matrix
        NULLIFY(linearMatrix%matrix)
        CALL DistributedMatrix_CreateStart(rowDomainMap,columnDomainMap,linearMatrix%matrix,err,error,*999)
        CALL DistributedMatrix_DataTypeSet(linearMatrix%matrix,MATRIX_VECTOR_DP_TYPE,err,error,*999)
        CALL DistributedMatrix_StorageTypeSet(linearMatrix%matrix,linearMatrix%storageType,err,error,*999)
        CALL DistributedMatrix_TransposeTypeSet(linearMatrix%matrix,DISTRIBUTED_MATRIX_PARTIAL_TRANSPOSE_REQUIRED,err,error,*999)
        !Calculate and set the matrix structure/sparsity pattern
        IF(linearMatrix%storageType/=DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE.AND. &
          & linearMatrix%storageType/=DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE) THEN
          NULLIFY(list)
          CALL EquationsMatrix_StructureCalculate(linearMatrix,numberOfNonZeros,rowIndices,columnIndices,list,err,error,*999)
          CALL DistributedMatrix_LinkListSet(linearMatrix%matrix,list,err,error,*999)
          CALL DistributedMatrix_NumberOfNonZerosSet(linearMatrix%matrix,numberOfNonZeros,err,error,*999)
          CALL DistributedMatrix_StorageLocationsSet(linearMatrix%matrix,rowIndices,columnIndices,err,error,*999)
          IF(ASSOCIATED(rowIndices)) DEALLOCATE(rowIndices)
          IF(ASSOCIATED(columnIndices)) DEALLOCATE(columnIndices)
        ENDIF
        CALL DistributedMatrix_CreateFinish(linearMatrix%matrix,err,error,*999)
      ENDDO !matrixIdx
    ENDIF
    NULLIFY(nonlinearMatrices)
    CALL EquationsMatricesVector_NonlinearMatricesExists(vectorMatrices,nonlinearMatrices,err,error,*999)
    IF(ASSOCIATED(nonlinearMatrices)) THEN
      !Nonlinear matrices
      NULLIFY(nonlinearMapping)
      CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)
      DO residualIdx=1,nonlinearMapping%numberOfResiduals
        NULLIFY(residualMapping)
        CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,residualIdx,residualMapping,err,error,*999)
        NULLIFY(residualVector)
        CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,residualIdx,residualVector,err,error,*999)
        !Set up the residual vector
        NULLIFY(residualVector%residual)
        CALL DistributedVector_CreateStart(rowDomainMap,residualVector%residual,err,error,*999)
        CALL DistributedVector_DataTypeSet(residualVector%residual,MATRIX_VECTOR_DP_TYPE,err,error,*999)
        CALL DistributedVector_CreateFinish(residualVector%residual,err,error,*999)
        !Initialise the residual vector to zero for time dependent problems so that the previous residual is set to zero
        CALL DistributedVector_AllValuesSet(residualVector%residual,0.0_DP,err,error,*999)
        !Set up the Jacobian matrices
        DO matrixIdx=1,residualVector%numberOfJacobians
          NULLIFY(jacobianMatrix)
          CALL EquationsMatricesResidual_JacobianMatrixGet(residualVector,matrixIdx,jacobianMatrix,err,error,*999)
          NULLIFY(jacobianMatrixToVarMap)
          CALL EquationsMappingResidual_JacobianMatrixToVarMapGet(residualMapping,matrixIdx,jacobianMatrixToVarMap,err,error,*999)
          NULLIFY(columnDomainMap)
          CALL EquationsMappingJMToVMap_ColumnDOFsMapGet(jacobianMatrixToVarMap,columnDomainMap,err,error,*999)
!!TODO: Set the distributed matrix not to allocate the data if the Jacobian is not calculated.
          !Create the distributed Jacobian matrix
          NULLIFY(jacobianMatrix%jacobian)
          CALL DistributedMatrix_CreateStart(rowDomainMap,columnDomainMap,jacobianMatrix%jacobian,err,error,*999)
          CALL DistributedMatrix_DataTypeSet(jacobianMatrix%jacobian,MATRIX_VECTOR_DP_TYPE,err,error,*999)
          CALL DistributedMatrix_StorageTypeSet(jacobianMatrix%jacobian,jacobianMatrix%storageType,err,error,*999)
          CALL DistributedMatrix_TransposeTypeSet(jacobianMatrix%jacobian,DISTRIBUTED_MATRIX_PARTIAL_TRANSPOSE_REQUIRED, &
            & err,error,*999)
          !Calculate and set the matrix structure/sparsity pattern
          IF(jacobianMatrix%storageType/=DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE.AND. &
            & jacobianMatrix%storageType/=DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE) THEN
            CALL JacobianMatrix_StructureCalculate(jacobianMatrix,numberOfNonZeros,rowIndices,columnIndices,err,error,*999)
            CALL DistributedMatrix_NumberOfNonZerosSet(jacobianMatrix%jacobian,numberOfNonZeros,err,error,*999)
            CALL DistributedMatrix_StorageLocationsSet(jacobianMatrix%jacobian,rowIndices,columnIndices,err,error,*999)
            IF(ASSOCIATED(rowIndices)) DEALLOCATE(rowIndices)
            IF(ASSOCIATED(columnIndices)) DEALLOCATE(columnIndices)
          ENDIF
          CALL DistributedMatrix_CreateFinish(jacobianMatrix%jacobian,err,error,*999)
        ENDDO !matrixIdx
      ENDDO !residualIdx
    ENDIF
    NULLIFY(rhsVector)
    CALL EquationsMatricesVector_RHSVectorExists(vectorMatrices,rhsVector,err,error,*999)
    IF(ASSOCIATED(rhsVector)) THEN      
      !Set up the equations RHS vector
      NULLIFY(rhsVector%vector)
      CALL DistributedVector_CreateStart(rowDomainMap,rhsVector%vector,err,error,*999)
      CALL DistributedVector_DataTypeSet(rhsVector%vector,MATRIX_VECTOR_DP_TYPE,err,error,*999)
      CALL DistributedVector_CreateFinish(rhsVector%vector,err,error,*999)
    ENDIF
    NULLIFY(sourceVectors)
    CALL EquationsMatricesVector_SourceVectorsExists(vectorMatrices,sourceVectors,err,error,*999)
    IF(ASSOCIATED(sourceVectors)) THEN
      DO sourceIdx=1,sourceVectors%numberOfSources
        NULLIFY(sourceVector)
        CALL EquationsMatricesSources_SourceVectorGet(sourceVectors,sourceIdx,sourceVector,err,error,*999)
        !Set up the equations source vector
        NULLIFY(sourceVector%vector)
        CALL DistributedVector_CreateStart(rowDomainMap,sourceVector%vector,err,error,*999)
        CALL DistributedVector_DataTypeSet(sourceVector%vector,MATRIX_VECTOR_DP_TYPE,err,error,*999)
        CALL DistributedVector_CreateFinish(sourceVector%vector,err,error,*999)
      ENDDO !sourceIdx
    ENDIF
    !Finish up
    vectorMatrices%vectorMatricesFinished=.TRUE.
        
    EXITS("EquationsMatrices_VectorCreateFinish")
    RETURN
999 IF(ASSOCIATED(rowIndices)) DEALLOCATE(rowIndices)
    IF(ASSOCIATED(columnIndices)) DEALLOCATE(columnIndices)
    CALL EquationsMatrices_VectorFinalise(vectorMatrices,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMatrices_VectorCreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_VectorCreateFinish

  !
  !================================================================================================================================
  !

  !>Starts the creation of the vector equations matrices and rhs for the vector equations
  SUBROUTINE EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsVectorType), POINTER :: vectorEquations !<The pointer to the vector equations to create the vector equations matrices for
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<On return, a pointer to the vector equations matrices being created.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string  
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(EquationsType), POINTER :: equations
    TYPE(VARYING_STRING) :: dummyError    

    ENTERS("EquationsMatrices_VectorCreateStart",err,error,*998)

    IF(ASSOCIATED(vectorMatrices)) CALL FlagError("Vector equations matrices is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorEquations)) CALL FlagError("Vector equations is not associated.",err,error,*998)   
    NULLIFY(equations)
    CALL EquationsVector_EquationsGet(vectorEquations,equations,err,error,*999)
    CALL Equations_AssertIsFinished(equations,err,error,*999)

    !Initialise the equations matrices
    CALL EquationsMatrices_VectorInitialise(vectorEquations,err,error,*999)
    vectorMatrices=>vectorEquations%vectorMatrices
    
    EXITS("EquationsMatrices_VectorCreateStart")
    RETURN
999 CALL EquationsMatrices_VectorFinalise(vectorEquations%vectorMatrices,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMatrices_VectorCreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_VectorCreateStart

  !
  !================================================================================================================================
  !

  !>Destroy the vector equations matrices
  SUBROUTINE EquationsMatrices_VectorDestroy(vectorMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer the vector equations matrices to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsMatrices_VectorDestroy",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector equations matrices is not associated",err,error,*999)
    
    CALL EquationsMatrices_VectorFinalise(vectorMatrices,err,error,*999)
        
    EXITS("EquationsMatrices_VectorDestroy")
    RETURN
999 ERRORSEXITS("EquationsMatrices_VectorDestroy",err,error)    
    RETURN 1
   
  END SUBROUTINE EquationsMatrices_VectorDestroy

  !
  !================================================================================================================================
  !

  !>Calculate the positions in the equations matrices of the element matrix. Old CMISS name MELGE.
  SUBROUTINE EquationsMatrices_ElementMatrixCalculate(elementMatrix,updateMatrix,rowElementNumbers,columnElementNumbers, &
    & rowsFieldVariable,colsFieldVariable,err,error,*)

    !Argument variables
    TYPE(ElementMatrixType) :: elementMatrix !<The element matrix to calculate
    LOGICAL :: updateMatrix !<Is .TRUE. if the element matrix is to be updated, .FALSE. if not.
    INTEGER(INTG), INTENT(IN) :: rowElementNumbers(:) !<The row element number to calculate
    INTEGER(INTG), INTENT(IN) :: columnElementNumbers(:) !<The column element number to calculate
    TYPE(FieldVariableType), POINTER :: rowsFieldVariable !<A pointer to the field variable associated with the rows
    TYPE(FieldVariableType), POINTER :: colsFieldVariable !<A pointer to the field variable associated with the columns
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: colElementNumber,colsVariableType,componentIdx,dataPointIdx,derivativeNumber,derivativeIdx, &
      & elementIdx,globalDOFIdx,interpolationType,localDataPointNumber,localDOFIdx,localNodeIdx,nodeNumber, &
      & numberOfComponents,numberOfElementDataPoints,numberOfLocalNodes,numberOfNodeDerivatives,rowElementNumber, &
      & rowsVariableType,totalNumberOfElements,versionNumber    
    TYPE(BasisType), POINTER :: basis
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionDataPointsType), POINTER :: dataPoints
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainElementsType), POINTER :: domainElements
    TYPE(DomainMappingType), POINTER :: domainMapping
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsMatrices_ElementMatrixCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(rowsFieldVariable)) CALL FlagError("Rows field variable is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(colsFieldVariable)) CALL FlagError("Columns field variable is not associated.",err,error,*999)
    CALL FieldVariable_VariableTypeGet(rowsFieldVariable,rowsVariableType,err,error,*999)
    CALL FieldVariable_VariableTypeGet(colsFieldVariable,colsVariableType,err,error,*999)
    
    elementMatrix%numberOfRows=0
    elementMatrix%numberOfColumns=0
    IF(updateMatrix) THEN
      IF(ASSOCIATED(rowsFieldVariable,colsFieldVariable)) THEN
        !Row and columns variable is the same.
        CALL FieldVariable_NumberOfComponentsGet(rowsFieldVariable,numberOfComponents,err,error,*999)
        DO componentIdx=1,numberOfComponents
          NULLIFY(domain)
          CALL FieldVariable_ComponentDomainGet(rowsFieldVariable,componentIdx,domain,err,error,*999)
          NULLIFY(domainTopology)
          CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
          NULLIFY(domainElements)
          CALL DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*999)
          CALL FieldVariable_ComponentInterpolationGet(rowsFieldVariable,componentIdx,interpolationType,err,error,*999)
          NULLIFY(domainMapping)
          CALL FieldVariable_DomainMappingGet(rowsFieldVariable,domainMapping,err,error,*999)
          CALL DomainElements_TotalNumberOfElementsGet(domainElements,totalNumberOfElements,err,error,*999)
          DO elementIdx=1,SIZE(rowElementNumbers)
            rowElementNumber=rowElementNumbers(elementIdx)
            IF(rowElementNumber<1.OR.rowElementNumber>totalNumberOfElements) THEN
              localError="Element number "//TRIM(NumberToVString(rowElementNumber,"*",err,error))// &
                & " is invalid for component number "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
                & " of rows field variable type "//TRIM(NumberToVString(rowsVariableType,"*",err,error))// &
                & ". The element number must be >= 1 and <= "//TRIM(NumberToVString(totalNumberOfElements,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
            SELECT CASE(interpolationType)
            CASE(FIELD_CONSTANT_INTERPOLATION)
              CALL FieldVariable_ConstantDOFGet(rowsFieldVariable,componentIdx,localDOFIdx,err,error,*999)
              CALL DomainMapping_LocalToGlobalGet(domainMapping,localDOFIdx,globalDOFIdx,err,error,*999)
              elementMatrix%numberOfRows=elementMatrix%numberOfRows+1
              elementMatrix%numberOfColumns=elementMatrix%numberOfColumns+1
              elementMatrix%rowDOFS(elementMatrix%numberOfRows)=localDOFIdx
              elementMatrix%columnDOFS(elementMatrix%numberOfColumns)=globalDOFIdx
            CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
              CALL FieldVariable_LocalElementDOFGet(rowsFieldVariable,rowElementNumber,componentIdx,localDOFIdx,err,error,*999)
              CALL DomainMapping_LocalToGlobalGet(domainMapping,localDOFIdx,globalDOFIdx,err,error,*999)
              elementMatrix%numberOfRows=elementMatrix%numberOfRows+1
              elementMatrix%numberOfColumns=elementMatrix%numberOfColumns+1
              elementMatrix%rowDOFS(elementMatrix%numberOfRows)=localDOFIdx
              elementMatrix%columnDOFS(elementMatrix%numberOfColumns)=globalDOFIdx
            CASE(FIELD_NODE_BASED_INTERPOLATION)
              NULLIFY(basis)
              CALL DomainElements_ElementBasisGet(domainElements,rowElementNumber,basis,err,error,*999)
              CALL Basis_NumberOfLocalNodesGet(basis,numberOfLocalNodes,err,error,*999)
              DO localNodeIdx=1,numberOfLocalNodes
                CALL DomainElements_ElementNodeGet(domainElements,localNodeIdx,rowElementNumber,nodeNumber,err,error,*999)
                CALL Basis_NodeNumberOfDerivativesGet(basis,localNodeIdx,numberOfNodeDerivatives,err,error,*999)
                DO derivativeIdx=1,numberOfNodeDerivatives
                  CALL DomainElements_ElementDerivativeGet(domainElements,derivativeIdx,localNodeIdx,rowElementNumber, &
                    & derivativeNumber,err,error,*999)
                  CALL DomainElements_ElementVersionGet(domainElements,derivativeIdx,localNodeIdx,rowElementNumber, &
                    & versionNumber,err,error,*999)
                  CALL FieldVariable_LocalNodeDOFGet(rowsFieldVariable,versionNumber,derivativeNumber,nodeNumber, &
                    & componentIdx,localDOFIdx,err,error,*999)
                  CALL DomainMapping_LocalToGlobalGet(domainMapping,localDOFIdx,globalDOFIdx,err,error,*999)
                  elementMatrix%numberOfRows=elementMatrix%numberOfRows+1
                  elementMatrix%numberOfColumns=elementMatrix%numberOfColumns+1
                  elementMatrix%rowDOFS(elementMatrix%numberOfRows)=localDOFIdx
                  elementMatrix%columnDOFS(elementMatrix%numberOfColumns)=globalDOFIdx
                ENDDO !derivativeIdx
              ENDDO !localNodeIdx
            CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(FIELD_DATA_POINT_BASED_INTERPOLATION)
              NULLIFY(decomposition)
              CALL Domain_DecompositionGet(domain,decomposition,err,error,*999)
              NULLIFY(decompositionTopology)
              CALL Decomposition_DecompositionTopologyGet(decomposition,decompositionTopology,err,error,*999)
              NULLIFY(dataPoints)
              CALL DecompositionTopology_DataPointsGet(decompositionTopology,dataPoints,err,error,*999)
              CALL DecompositionDataPoints_ElementNumberOfDataPointsGet(dataPoints,rowElementNumber,numberOfElementDataPoints, &
                & err,error,*999)
              DO dataPointIdx=1,numberOfElementDataPoints
                CALL DecompositionDataPoints_ElementDataLocalNumberGet(dataPoints,dataPointIdx,rowElementNumber, &
                  & localDataPointNumber,err,error,*999)
                CALL FieldVariable_LocalDataPointDOFGet(rowsFieldVariable,localDataPointNumber,componentIdx,localDOFIdx, &
                  & err,error,*999)
                CALL DomainMapping_LocalToGlobalGet(domainMapping,localDOFIdx,globalDOFIdx,err,error,*999)
                elementMatrix%numberOfRows=elementMatrix%numberOfRows+1
                elementMatrix%numberOfColumns=elementMatrix%numberOfColumns+1
                elementMatrix%rowDOFS(elementMatrix%numberOfRows)=localDOFIdx
                elementMatrix%columnDOFS(elementMatrix%numberOfColumns)=globalDOFIdx
              ENDDO !dataPointIdx
            CASE DEFAULT
              localError="The interpolation type of "// &
                & TRIM(NumberToVString(interpolationType,"*",err,error))// &
                & " is invalid for component number "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
                & " of rows field variable type "//TRIM(NumberToVString(rowsFieldVariable%variableType,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)          
            END SELECT
          ENDDO !elementIdx
        ENDDO !componentIdx
      ELSE
        !Row and column variables are different
        !Row mapping
        CALL FieldVariable_NumberOfComponentsGet(rowsFieldVariable,numberOfComponents,err,error,*999)
        DO componentIdx=1,numberOfComponents
          NULLIFY(domain)
          CALL FieldVariable_ComponentDomainGet(rowsFieldVariable,componentIdx,domain,err,error,*999)
          NULLIFY(domainTopology)
          CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
          NULLIFY(domainElements)
          CALL DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*999)
          CALL FieldVariable_ComponentInterpolationGet(rowsFieldVariable,componentIdx,interpolationType,err,error,*999)
          NULLIFY(domainMapping)
          CALL FieldVariable_DomainMappingGet(rowsFieldVariable,domainMapping,err,error,*999)
          CALL DomainElements_TotalNumberOfElementsGet(domainElements,totalNumberOfElements,err,error,*999)
          DO elementIdx=1,SIZE(rowElementNumbers)
            rowElementNumber=rowElementNumbers(elementIdx)
            IF(rowElementNumber<1.OR.rowElementNumber>totalNumberOfElements) THEN
              localError="Element number "//TRIM(NumberToVString(rowElementNumber,"*",err,error))// &
                & " is invalid for component number "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
                & " of rows field variable type "//TRIM(NumberToVString(rowsVariableType,"*",err,error))// &
                & ". The element number must be >= 1 and <= "//TRIM(NumberToVString(totalNumberOfElements,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
           SELECT CASE(interpolationType)
            CASE(FIELD_CONSTANT_INTERPOLATION)
              CALL FieldVariable_ConstantDOFGet(rowsFieldVariable,componentIdx,localDOFIdx,err,error,*999)
              elementMatrix%numberOfRows=elementMatrix%numberOfRows+1
              elementMatrix%rowDOFS(elementMatrix%numberOfRows)=localDOFIdx
            CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
              CALL FieldVariable_LocalElementDOFGet(rowsFieldVariable,rowElementNumber,componentIdx,localDOFIdx,err,error,*999)
              elementMatrix%numberOfRows=elementMatrix%numberOfRows+1
              elementMatrix%rowDOFS(elementMatrix%numberOfRows)=localDOFIdx
            CASE(FIELD_NODE_BASED_INTERPOLATION)
              NULLIFY(basis)
              CALL DomainElements_ElementBasisGet(domainElements,rowElementNumber,basis,err,error,*999)
              CALL Basis_NumberOfLocalNodesGet(basis,numberOfLocalNodes,err,error,*999)
              DO localNodeIdx=1,numberOfLocalNodes
                CALL DomainElements_ElementNodeGet(domainElements,localNodeIdx,rowElementNumber,nodeNumber,err,error,*999)
                CALL Basis_NodeNumberOfDerivativesGet(basis,localNodeIdx,numberOfNodeDerivatives,err,error,*999)
                DO derivativeIdx=1,numberOfNodeDerivatives
                  CALL DomainElements_ElementDerivativeGet(domainElements,derivativeIdx,localNodeIdx,rowElementNumber, &
                    & derivativeNumber,err,error,*999)
                  CALL DomainElements_ElementVersionGet(domainElements,derivativeIdx,localNodeIdx,rowElementNumber, &
                    & versionNumber,err,error,*999)
                  CALL FieldVariable_LocalNodeDOFGet(rowsFieldVariable,versionNumber,derivativeNumber,nodeNumber, &
                    & componentIdx,localDOFIdx,err,error,*999)
                  elementMatrix%numberOfRows=elementMatrix%numberOfRows+1
                  elementMatrix%rowDOFS(elementMatrix%numberOfRows)=localDOFIdx
                ENDDO !derivativeIdx
              ENDDO !localNodeIdx
            CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(FIELD_DATA_POINT_BASED_INTERPOLATION)
              NULLIFY(decomposition)
              CALL Domain_DecompositionGet(domain,decomposition,err,error,*999)
              NULLIFY(decompositionTopology)
              CALL Decomposition_DecompositionTopologyGet(decomposition,decompositionTopology,err,error,*999)
              NULLIFY(dataPoints)
              CALL DecompositionTopology_DataPointsGet(decompositionTopology,dataPoints,err,error,*999)
              CALL DecompositionDataPoints_ElementNumberOfDataPointsGet(dataPoints,rowElementNumber,numberOfElementDataPoints, &
                & err,error,*999)
               DO dataPointIdx=1,numberOfElementDataPoints
                CALL DecompositionDataPoints_ElementDataLocalNumberGet(dataPoints,dataPointIdx,rowElementNumber, &
                  & localDataPointNumber,err,error,*999)
                CALL FieldVariable_LocalDataPointDOFGet(rowsFieldVariable,localDataPointNumber,componentIdx,localDOFIdx, &
                  & err,error,*999)
                elementMatrix%numberOfRows=elementMatrix%numberOfRows+1
                elementMatrix%rowDOFS(elementMatrix%numberOfRows)=localDOFIdx
              ENDDO !dataPointIdx
            CASE DEFAULT
              localError="The interpolation type of "// &
                & TRIM(NumberToVString(interpolationType,"*",err,error))// &
                & " is invalid for component number "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
                & " of rows field variable type "//TRIM(NumberToVString(rowsVariableType,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)          
            END SELECT
          ENDDO !elementIdx
        ENDDO !componentIdx
        !Column mapping
        CALL FieldVariable_NumberOfComponentsGet(colsFieldVariable,numberOfComponents,err,error,*999)
        DO componentIdx=1,numberOfComponents
          NULLIFY(domain)
          CALL FieldVariable_ComponentDomainGet(colsFieldVariable,componentIdx,domain,err,error,*999)
          NULLIFY(domainTopology)
          CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
          NULLIFY(domainElements)
          CALL DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*999)
          CALL FieldVariable_ComponentInterpolationGet(colsFieldVariable,componentIdx,interpolationType,err,error,*999)
          NULLIFY(domainMapping)
          CALL FieldVariable_DomainMappingGet(colsFieldVariable,domainMapping,err,error,*999)
          CALL DomainElements_TotalNumberOfElementsGet(domainElements,totalNumberOfElements,err,error,*999)
          DO elementIdx=1,SIZE(columnElementNumbers)
            colElementNumber=columnElementNumbers(elementIdx)
            IF(colElementNumber<1.AND.colElementNumber>totalNumberOfElements) THEN
              localError="Column element number "//TRIM(NumberToVString(colElementNumber,"*",err,error))// &
                & " is invalid for component number "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
                & " of column field variable type "//TRIM(NumberToVString(colsVariableType,"*",err,error))// &
                & ". The element number must be between 1 and "//TRIM(NumberToVString(totalNumberOfElements,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
            SELECT CASE(interpolationType)
            CASE(FIELD_CONSTANT_INTERPOLATION)
              CALL FieldVariable_ConstantDOFGet(colsFieldVariable,componentIdx,localDOFIdx,err,error,*999)
              CALL DomainMapping_LocalToGlobalGet(domainMapping,localDOFIdx,globalDOFIdx,err,error,*999)
              elementMatrix%numberOfColumns=elementMatrix%numberOfColumns+1
              elementMatrix%columnDOFS(elementMatrix%numberOfColumns)=globalDOFIdx
            CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
              CALL FieldVariable_LocalElementDOFGet(colsFieldVariable,colElementNumber,componentIdx,localDOFIdx,err,error,*999)
              CALL DomainMapping_LocalToGlobalGet(domainMapping,localDOFIdx,globalDOFIdx,err,error,*999)
              elementMatrix%numberOfColumns=elementMatrix%numberOfColumns+1
              elementMatrix%columnDOFS(elementMatrix%numberOfColumns)=globalDOFIdx
            CASE(FIELD_NODE_BASED_INTERPOLATION)
              NULLIFY(basis)
              CALL DomainElements_ElementBasisGet(domainElements,colElementNumber,basis,err,error,*999)
              CALL Basis_NumberOfLocalNodesGet(basis,numberOfLocalNodes,err,error,*999)
              DO localNodeIdx=1,numberOfLocalNodes
                CALL DomainElements_ElementNodeGet(domainElements,localNodeIdx,rowElementNumber,nodeNumber,err,error,*999)
                CALL Basis_NodeNumberOfDerivativesGet(basis,localNodeIdx,numberOfNodeDerivatives,err,error,*999)
                DO derivativeIdx=1,numberOfNodeDerivatives
                  CALL DomainElements_ElementDerivativeGet(domainElements,derivativeIdx,localNodeIdx,rowElementNumber, &
                    & derivativeNumber,err,error,*999)
                  CALL DomainElements_ElementVersionGet(domainElements,derivativeIdx,localNodeIdx,rowElementNumber, &
                    & versionNumber,err,error,*999)
                  CALL FieldVariable_LocalNodeDOFGet(colsFieldVariable,versionNumber,derivativeNumber,nodeNumber, &
                    & componentIdx,localDOFIdx,err,error,*999)
                  CALL DomainMapping_LocalToGlobalGet(domainMapping,localDOFIdx,globalDOFIdx,err,error,*999)
                  elementMatrix%numberOfColumns=elementMatrix%numberOfColumns+1
                  elementMatrix%columnDOFS(elementMatrix%numberOfColumns)=globalDOFIdx
                ENDDO !derivativeIdx
              ENDDO !localNodeIdx
            CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(FIELD_DATA_POINT_BASED_INTERPOLATION)
              NULLIFY(decomposition)
              CALL Domain_DecompositionGet(domain,decomposition,err,error,*999)
              NULLIFY(decompositionTopology)
              CALL Decomposition_DecompositionTopologyGet(decomposition,decompositionTopology,err,error,*999)
              NULLIFY(dataPoints)
              CALL DecompositionTopology_DataPointsGet(decompositionTopology,dataPoints,err,error,*999)
              CALL DecompositionDataPoints_ElementNumberOfDataPointsGet(dataPoints,rowElementNumber,numberOfElementDataPoints, &
                & err,error,*999)
              DO dataPointIdx=1,numberOfElementDataPoints
                CALL DecompositionDataPoints_ElementDataLocalNumberGet(dataPoints,dataPointIdx,rowElementNumber, &
                  & localDataPointNumber,err,error,*999)
                CALL FieldVariable_LocalDataPointDOFGet(colsFieldVariable,localDataPointNumber,componentIdx,localDOFIdx, &
                  & err,error,*999)
                CALL DomainMapping_LocalToGlobalGet(domainMapping,localDOFIdx,globalDOFIdx,err,error,*999)
                elementMatrix%numberOfColumns=elementMatrix%numberOfColumns+1
                elementMatrix%columnDOFS(elementMatrix%numberOfColumns)=globalDOFIdx
              ENDDO !dataPointIdx
            CASE DEFAULT
              localError="The interpolation type of "// &
                & TRIM(NumberToVString(interpolationType,"*",err,error))// &
                & " is invalid for component number "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
                & " of column field variable type "//TRIM(NumberToVString(colsVariableType,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)          
            END SELECT
          ENDDO !elementIdx
        ENDDO !componentIdx
      ENDIF
      elementMatrix%matrix=0.0_DP
    ENDIF
    
    EXITS("EquationsMatrices_ElementMatrixCalculate")
    RETURN
999 ERRORSEXITS("EquationsMatrices_ElementMatrixCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_ElementMatrixCalculate

  !
  !================================================================================================================================
  !

  !>Finalise an element matrix and deallocate all memory
  SUBROUTINE EquationsMatrices_ElementMatrixFinalise(elementMatrix,err,error,*)

    !Argument variables
    TYPE(ElementMatrixType):: elementMatrix!<The element matrix to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("EquationsMatrices_ElementMatrixFinalise",err,error,*999)

    elementMatrix%equationsMatrixNumber=0
    elementMatrix%numberOfRows=0
    elementMatrix%numberOfColumns=0
    elementMatrix%maxNumberOfRows=0
    elementMatrix%maxNumberOfColumns=0
    IF(ALLOCATED(elementMatrix%rowDOFS)) DEALLOCATE(elementMatrix%rowDOFS)
    IF(ALLOCATED(elementMatrix%columnDOFS)) DEALLOCATE(elementMatrix%columnDOFS)
    IF(ALLOCATED(elementMatrix%matrix)) DEALLOCATE(elementMatrix%matrix)
    
    EXITS("EquationsMatrices_ElementMatrixFinalise")
    RETURN
999 ERRORSEXITS("EquationsMatrices_ElementMatrixFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_ElementMatrixFinalise

  !
  !================================================================================================================================
  !

  !>Initialise the element matrix.
  SUBROUTINE EquationsMatrices_ElementMatrixInitialise(elementMatrix,err,error,*)

    !Argument variables
    TYPE(ElementMatrixType) :: elementMatrix !The element matrix to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsMatrices_ElementMatrixInitialise",err,error,*999)

    elementMatrix%equationsMatrixNumber=0
    elementMatrix%numberOfRows=0
    elementMatrix%numberOfColumns=0
    elementMatrix%maxNumberOfRows=0
    elementMatrix%maxNumberOfColumns=0
       
    EXITS("EquationsMatrices_ElementMatrixInitialise")
    RETURN
999 ERRORSEXITS("EquationsMatrices_ElementMatrixInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_ElementMatrixInitialise

  !
  !================================================================================================================================
  !

  !>Sets up the element matrix for the row and column field variables.
  SUBROUTINE EquationsMatrices_ElementMatrixSetup(elementMatrix,rowsFieldVariable,colsFieldVariable,rowsNumberOfElements, &
    & colsNumberOfElements,err,error,*)

    !Argument variables
    TYPE(ElementMatrixType) :: elementMatrix !<The element matrix to setup
    TYPE(FieldVariableType), POINTER :: rowsFieldVariable !<A pointer to the field variable associated with the rows
    TYPE(FieldVariableType), POINTER :: colsFieldVariable !<A pointer to the field variable associated with the columns
    INTEGER(INTG), INTENT(IN)  :: rowsNumberOfElements !<Number of elements in the row variables whose dofs are present in this element matrix
    INTEGER(INTG), INTENT(IN)  :: colsNumberOfElements !<Number of elements in the col variables whose dofs are present in this element matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,componentIdx,maxElementInterpParameters,numberOfComponents
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("EquationsMatrices_ElementMatrixSetup",err,error,*998)

    IF(.NOT.ASSOCIATED(rowsFieldVariable)) CALL FlagError("Rows field variable is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(colsFieldVariable)) CALL FlagError("Columns field variable is not associated.",err,error,*999)
    IF(ALLOCATED(elementMatrix%rowDOFS)) CALL FlagError("Element matrix row dofs already allocated.",err,error,*999)
    IF(ALLOCATED(elementMatrix%columnDOFS)) CALL FlagError("Element matrix column dofs already allocated.",err,error,*999)
    IF(ALLOCATED(elementMatrix%matrix)) CALL FlagError("Element matrix already allocated.",err,error,*999)
     
    elementMatrix%maxNumberOfRows=0
    CALL FieldVariable_NumberOfComponentsGet(rowsFieldVariable,numberOfComponents,err,error,*999)
    DO componentIdx=1,numberOfComponents
      CALL FieldVariable_ComponentMaxElementInterpParametersGet(rowsFieldVariable,componentIdx,maxElementInterpParameters, &
        & err,error,*999)
      elementMatrix%maxNumberOfRows=elementMatrix%maxNumberOfRows+maxElementInterpParameters
    ENDDO !componentIdx
    elementMatrix%maxNumberOfRows=elementMatrix%maxNumberOfRows*rowsNumberOfElements
    elementMatrix%maxNumberOfColumns=0
    CALL FieldVariable_NumberOfComponentsGet(colsFieldVariable,numberOfComponents,err,error,*999)
    DO componentIdx=1,numberOfComponents
      CALL FieldVariable_ComponentMaxElementInterpParametersGet(colsFieldVariable,componentIdx,maxElementInterpParameters, &
        & err,error,*999)
      elementMatrix%maxNumberOfColumns=elementMatrix%maxNumberOfColumns+maxElementInterpParameters
    ENDDO !componentIdx
    elementMatrix%maxNumberOfColumns=elementMatrix%maxNumberOfColumns*colsNumberOfElements
    ALLOCATE(elementMatrix%rowDOFS(elementMatrix%maxNumberOfRows),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate element matrix row dofs.",err,error,*999)
    ALLOCATE(elementMatrix%columnDOFS(elementMatrix%maxNumberOfColumns),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate element matrix column dofs.",err,error,*999)
    ALLOCATE(elementMatrix%matrix(elementMatrix%maxNumberOfRows,elementMatrix%maxNumberOfColumns),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate element matrix.",err,error,*999)
    
    EXITS("EquationsMatrices_ElementMatrixSetup")
    RETURN
999 CALL EquationsMatrices_ElementMatrixFinalise(elementMatrix,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMatrices_ElementMatrixSetup",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_ElementMatrixSetup

  !
  !================================================================================================================================
  !

  !>Calculate the positions in the equations rhs of the element rhs vector. Old CMISS name MELGE.
  SUBROUTINE EquationsMatrices_ElementVectorCalculate(elementVector,updateVector,elementNumber,rowsFieldVariable,err,error,*)

    !Argument variables
    TYPE(ElementVectorType) :: elementVector !<The element vector to calculate.
    LOGICAL :: updateVector !<Is .TRUE. if the element vector is to be updated, .FALSE. if not.
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to calculate
    TYPE(FieldVariableType), POINTER :: rowsFieldVariable !<A pointer to the field variable associated with the rows
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx,dataPointIdx,derivativeNumber,derivativeIdx,interpolationType,localDataPointNumber, &
      & localDOFIdx,localNodeIdx,nodeNumber,numberOfComponents,numberOfNodeDerivatives,numberOfElementDataPoints, &
      & numberOfLocalNodes,rowsVariableType,totalNumberOfElements,versionNumber
    TYPE(BasisType), POINTER :: basis
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionDataPointsType), POINTER :: dataPoints
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainElementsType), POINTER :: domainElements
    TYPE(DomainMappingType), POINTER :: domainMapping
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(DomainElementsType), POINTER :: elementsTopology
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("EquationsMatrices_ElementVectorCalculate",err,error,*999)

    CALL FieldVariable_VariableTypeGet(rowsFieldVariable,rowsVariableType,err,error,*999)
    
    !Calculate the rows for the element vector
    elementVector%numberOfRows=0
    IF(updateVector) THEN
      CALL FieldVariable_NumberOfComponentsGet(rowsFieldVariable,numberOfComponents,err,error,*999)
      DO componentIdx=1,numberOfComponents
        NULLIFY(domain)
        CALL FieldVariable_ComponentDomainGet(rowsFieldVariable,componentIdx,domain,err,error,*999)
        NULLIFY(domainTopology)
        CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
        NULLIFY(domainElements)
        CALL DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*999)
        CALL FieldVariable_ComponentInterpolationGet(rowsFieldVariable,componentIdx,interpolationType,err,error,*999)
        NULLIFY(domainMapping)
        CALL FieldVariable_DomainMappingGet(rowsFieldVariable,domainMapping,err,error,*999)
        CALL DomianMapping_TotalNumberOfElementsGet(domain,totalNumberOfElements,err,error,*999)
        IF(elementNumber<1.OR.elementNumber>totalNumberOfElements) THEN
          localError="Element number "//TRIM(NumberToVString(elementNumber,"*",err,error))// &
            & " is invalid for component number "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
            & " of rows field variable type "//TRIM(NumberToVString(rowsVariableType,"*",err,error))// &
            & ". The element number must be >= 1 and <= "//TRIM(NumberToVString(totalNumberOfElements,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        SELECT CASE(interpolationType)
        CASE(FIELD_CONSTANT_INTERPOLATION)
          CALL FieldVariable_ConstantDOFGet(rowsFieldVariable,componentIdx,localDOFIdx,err,error,*999)
          elementVector%numberOfRows=elementVector%numberOfRows+1
          elementVector%rowDOFS(elementVector%numberOfRows)=localDOFIdx
        CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
          CALL FieldVariable_LocalElementDOFGet(rowsFieldVariable,elementNumber,componentIdx,localDOFIdx,err,error,*999)
          elementVector%numberOfRows=elementVector%numberOfRows+1
          elementVector%rowDOFS(elementVector%numberOfRows)=localDOFIdx
        CASE(FIELD_NODE_BASED_INTERPOLATION)
          NULLIFY(basis)
          CALL DomainElements_ElementBasisGet(domainElements,elementNumber,basis,err,error,*999)
          CALL Basis_NumberOfLocalNodesGet(basis,numberOfLocalNodes,err,error,*999)
          DO localNodeIdx=1,numberOfLocalNodes
            CALL DomainElements_ElementNodeGet(domainElements,localNodeIdx,elementNumber,nodeNumber,err,error,*999)
            CALL Basis_NodeNumberOfDerivativesGet(basis,localNodeIdx,numberOfNodeDerivatives,err,error,*999)
            DO derivativeIdx=1,numberOfNodeDerivatives
              CALL DomainElements_ElementDerivativeGet(domainElements,derivativeIdx,localNodeIdx,elementNumber, &
                & derivativeNumber,err,error,*999)
              CALL DomainElements_ElementVersionGet(domainElements,derivativeIdx,localNodeIdx,elementNumber, &
                & versionNumber,err,error,*999)
              CALL FieldVariable_LocalNodeDOFGet(rowsFieldVariable,versionNumber,derivativeNumber,nodeNumber, &
                & componentIdx,localDOFIdx,err,error,*999)
              elementVector%numberOfRows=elementVector%numberOfRows+1
              elementVector%rowDOFS(elementVector%numberOfRows)=localDOFIdx
            ENDDO !derivativeIdx
          ENDDO !localNodeIdx
        CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(FIELD_DATA_POINT_BASED_INTERPOLATION)
          NULLIFY(decomposition)
          CALL Domain_DecompositionGet(domain,decomposition,err,error,*999)
          NULLIFY(decompositionTopology)
          CALL Decomposition_DecompositionTopologyGet(decomposition,decompositionTopology,err,error,*999)
          NULLIFY(dataPoints)
          CALL DecompositionTopology_DataPointsGet(decompositionTopology,dataPoints,err,error,*999)
          CALL DecompositionDataPoints_ElementNumberOfDataPointsGet(dataPoints,elementNumber,numberOfElementDataPoints, &
            & err,error,*999)
           DO dataPointIdx=1,numberOfElementDataPoints
             CALL DecompositionDataPoints_ElementDataLocalNumberGet(dataPoints,dataPointIdx,elementNumber, &
               & localDataPointNumber,err,error,*999)
            CALL FieldVariable_LocalDataPointDOFGet(rowsFieldVariable,localDataPointNumber,componentIdx,localDOFIdx,err,error,*999)
            elementVector%numberOfRows=elementVector%numberOfRows+1
            elementVector%rowDOFS(elementVector%numberOfRows)=localDOFIdx
          ENDDO !dataPointIdx
        CASE DEFAULT
          localError="The interpolation type of "// &
            & TRIM(NumberToVString(interpolationType,"*",err,error))// &
            & " is invalid for component number "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
            & " of rows field variable type "//TRIM(NumberToVString(rowsVariableType,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)          
        END SELECT
      ENDDO !componentIdx
      elementVector%vector=0.0_DP
    ENDIF
    
    EXITS("EquationsMatrices_ElementVectorCalculate")
    RETURN
999 ERRORSEXITS("EquationsMatrices_ElementVectorCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_ElementVectorCalculate

  !
  !================================================================================================================================
  !

  !>Finalise an element vector and deallocate all memory
  SUBROUTINE EquationsMatrices_ElementVectorFinalise(elementVector,err,error,*)

    !Argument variables
    TYPE(ElementVectorType) :: elementVector !<The element vector to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("EquationsMatrices_ElementVectorFinalise",err,error,*999)

    elementVector%numberOfRows=0
    elementVector%maxNumberOfRows=0    
    IF(ALLOCATED(elementVector%rowDOFS)) DEALLOCATE(elementVector%rowDOFS)    
    IF(ALLOCATED(elementVector%vector)) DEALLOCATE(elementVector%vector)
    
    EXITS("EquationsMatrices_ElementVectorFinalise")
    RETURN
999 ERRORSEXITS("EquationsMatrices_ElementVectorFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_ElementVectorFinalise

  !
  !================================================================================================================================
  !

  !>Initialise the element vector
  SUBROUTINE EquationsMatrices_ElementVectorInitialise(elementVector,err,error,*)

    !Argument variables
    TYPE(ElementVectorType) :: elementVector !The element vector to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsMatrices_ElementVectorInitialise",err,error,*999)

    elementVector%numberOfRows=0
    elementVector%maxNumberOfRows=0
       
    EXITS("EquationsMatrices_ElementVectorInitialise")
    RETURN
999 ERRORSEXITS("EquationsMatrices_ElementVectorInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_ElementVectorInitialise

  !
  !================================================================================================================================
  !

  !>Sets up the element vector for the row field variables.
  SUBROUTINE EquationsMatrices_ElementVectorSetup(elementVector,rowsFieldVariable,err,error,*)

    !Argument variables
    TYPE(ElementVectorType) :: elementVector !<The element vector to setup
    TYPE(FieldVariableType), POINTER :: rowsFieldVariable !<A pointer to the field variable associated with the rows
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx,dummyErr,maxElementInterpParameters,numberOfComponents
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("EquationsMatrices_ElementVectorSetup",err,error,*998)

    IF(.NOT.ASSOCIATED(rowsFieldVariable)) CALL FlagError("Rows field variable is not associated.",err,error,*999)
    IF(ALLOCATED(elementVector%rowDOFS)) CALL FlagError("Element vector row dofs is already allocated.",err,error,*999)
    IF(ALLOCATED(elementVector%vector)) CALL FlagError("Element vector vector already allocated.",err,error,*999)
   
    elementVector%maxNumberOfRows = 0
    CALL FieldVariable_NumberOfComponentsGet(rowsFieldVariable,numberOfComponents,err,error,*999)
    DO componentIdx=1,numberOfComponents
      CALL FieldVariable_ComponentMaxElementInterpParametersGet(rowsFieldVariable,componentIdx,maxElementInterpParameters, &
        & err,error,*999)
      elementVector%maxNumberOfRows=elementVector%maxNumberOfRows+maxElementInterpParameters
    ENDDO !componentIdx
    ALLOCATE(elementVector%rowDOFS(elementVector%maxNumberOfRows),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate element vector row dofs.",err,error,*999)
    ALLOCATE(elementVector%vector(elementVector%maxNumberOfRows),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate element vector vector.",err,error,*999)
    
    EXITS("EquationsMatrices_ElementVectorSetup")
    RETURN
999 CALL EquationsMatrices_ElementVectorFinalise(elementVector,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMatrices_ElementVectorSetup",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_ElementVectorSetup

  !
  !================================================================================================================================
  !

  !>Ensure that a particular temporal type distributed vector has been created for a equations matrices residual vector.
  SUBROUTINE EquationsMatricesResidual_DistributedVectorEnsureCreated(residualVector,temporalType,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector !<The residual vector to ensure that the distributed vector is created for.
    INTEGER(INTG), INTENT(IN) :: temporalType !<The equations matrices vector temporal type to ensure that it has been created. \see EquationsMatricesRoutines_VectorTemporalTypes
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DomainMappingType), POINTER :: rowsDomainMap
    TYPE(DistributedVectorType), POINTER :: residualDistributedVector
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(FieldVariableType), POINTER :: lhsVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsMatricesResidual_DistributedVectorEnsureCreated",err,error,*999)

    NULLIFY(residualDistributedVector)
    CALL EquationsMatricesResidual_DistributedVectorExists(residualVector,temporalType,residualDistributedVector,err,error,*999)
    IF(.NOT.ASSOCIATED(residualDistributedVector)) THEN
      NULLIFY(nonlinearMatrices)
      CALL EquationsMatricesResidual_NonlinearMatricesGet(residualVector,nonlinearMatrices,err,error,*999)
      NULLIFY(vectorMatrices)
      CALL EquationsMatricesNonlinear_VectorMatricesGet(nonlinearMatrices,vectorMatrices,err,error,*999)
      NULLIFY(vectorMapping)
      CALL EquationsMatricesVector_VectorMappingGet(vectorMatrices,vectorMapping,err,error,*999)
      NULLIFY(lhsMapping)
      CALL EquationsMappingVector_LHSMappingGet(vectorMapping,lhsMapping,err,error,*999)
      NULLIFY(lhsVariable)
      CALL EquationsMappingLHS_LHSVariableGet(lhsMapping,lhsVariable,err,error,*999)
      NULLIFY(rowsDomainMap)
      CALL FieldVariable_DomainMappingGet(lhsVariable,rowsDomainMap,err,error,*999)
      SELECT CASE(temporalType)
      CASE(EQUATIONS_MATRICES_CURRENT_VECTOR)
        NULLIFY(residualVector%residual)
        CALL DistributedVector_CreateStart(rowsDomainMap,residualVector%residual,err,error,*999)
        CALL DistributedVector_DataTypeSet(residualVector%residual,MATRIX_VECTOR_DP_TYPE,err,error,*999)
        CALL DistributedVector_CreateFinish(residualVector%residual,err,error,*999)
        CALL DistributedVector_AllValuesSet(residualVector%residual,0.0_DP,err,error,*999)
      CASE(EQUATIONS_MATRICES_PREVIOUS_VECTOR)
        NULLIFY(residualVector%previousResidual)
        CALL DistributedVector_CreateStart(rowsDomainMap,residualVector%previousResidual,err,error,*999)
        CALL DistributedVector_DataTypeSet(residualVector%previousResidual,MATRIX_VECTOR_DP_TYPE,err,error,*999)
        CALL DistributedVector_CreateFinish(residualVector%previousResidual,err,error,*999)
        CALL DistributedVector_AllValuesSet(residualVector%previousResidual,0.0_DP,err,error,*999)
      CASE(EQUATIONS_MATRICES_PREVIOUS2_VECTOR)
        NULLIFY(residualVector%previous2Residual)
        CALL DistributedVector_CreateStart(rowsDomainMap,residualVector%previous2Residual,err,error,*999)
        CALL DistributedVector_DataTypeSet(residualVector%previous2Residual,MATRIX_VECTOR_DP_TYPE,err,error,*999)
        CALL DistributedVector_CreateFinish(residualVector%previous2Residual,err,error,*999)
        CALL DistributedVector_AllValuesSet(residualVector%previous2Residual,0.0_DP,err,error,*999)
      CASE(EQUATIONS_MATRICES_PREVIOUS3_VECTOR)
        NULLIFY(residualVector%previous3Residual)
        CALL DistributedVector_CreateStart(rowsDomainMap,residualVector%previous3Residual,err,error,*999)
        CALL DistributedVector_DataTypeSet(residualVector%previous3Residual,MATRIX_VECTOR_DP_TYPE,err,error,*999)
        CALL DistributedVector_CreateFinish(residualVector%previous3Residual,err,error,*999)
        CALL DistributedVector_AllValuesSet(residualVector%previous3Residual,0.0_DP,err,error,*999)
      CASE DEFAULT
        localError="The specified RHS vector temporal type of "//TRIM(NumberToVString(temporalType,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDIF
    
    EXITS("EquationsMatricesResidual_DistributedVectorEnsureCreated")
    RETURN
999 ERRORS("EquationsMatricesResidual_DistributedVectorEnsureCreated",err,error)
    EXITS("EquationsMatricesResidual_DistributedVectorEnsureCreated")
    RETURN 1
    
  END SUBROUTINE EquationsMatricesResidual_DistributedVectorEnsureCreated

  !
  !================================================================================================================================
  !

  !>Ensure that a particular temporal type distributed vector has been created for a equations matrices RHS vector.
  SUBROUTINE EquationsMatricesRHS_DistributedVectorEnsureCreated(rhsVector,temporalType,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector !<The RHS vector to ensure that the distributed vector is created for.
    INTEGER(INTG), INTENT(IN) :: temporalType !<The equations matrices vector temporal type to ensure that it has been created. \see EquationsMatricesRoutines_VectorTemporalTypes
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DomainMappingType), POINTER :: rowsDomainMap
    TYPE(DistributedVectorType), POINTER :: rhsDistributedVector
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(FieldVariableType), POINTER :: lhsVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsMatricesRHS_DistributedVectorEnsureCreated",err,error,*999)

    NULLIFY(rhsDistributedVector)
    CALL EquationsMatricesRHS_DistributedVectorExists(rhsVector,temporalType,rhsDistributedVector,err,error,*999)
    IF(.NOT.ASSOCIATED(rhsDistributedVector)) THEN
      NULLIFY(vectorMatrices)
      CALL EquationsMatricesRHS_VectorMatricesGet(rhsVector,vectorMatrices,err,error,*999)
      NULLIFY(vectorMapping)
      CALL EquationsMatricesVector_VectorMappingGet(vectorMatrices,vectorMapping,err,error,*999)
      NULLIFY(lhsMapping)
      CALL EquationsMappingVector_LHSMappingGet(vectorMapping,lhsMapping,err,error,*999)
      NULLIFY(lhsVariable)
      CALL EquationsMappingLHS_LHSVariableGet(lhsMapping,lhsVariable,err,error,*999)
      NULLIFY(rowsDomainMap)
      CALL FieldVariable_DomainMappingGet(lhsVariable,rowsDomainMap,err,error,*999)
      SELECT CASE(temporalType)
      CASE(EQUATIONS_MATRICES_CURRENT_VECTOR)
        NULLIFY(rhsVector%vector)
        CALL DistributedVector_CreateStart(rowsDomainMap,rhsVector%vector,err,error,*999)
        CALL DistributedVector_DataTypeSet(rhsVector%vector,MATRIX_VECTOR_DP_TYPE,err,error,*999)
        CALL DistributedVector_CreateFinish(rhsVector%vector,err,error,*999)
        CALL DistributedVector_AllValuesSet(rhsVector%vector,0.0_DP,err,error,*999)
      CASE(EQUATIONS_MATRICES_PREVIOUS_VECTOR)
        NULLIFY(rhsVector%previousRHSVector)
        CALL DistributedVector_CreateStart(rowsDomainMap,rhsVector%previousRHSVector,err,error,*999)
        CALL DistributedVector_DataTypeSet(rhsVector%previousRHSVector,MATRIX_VECTOR_DP_TYPE,err,error,*999)
        CALL DistributedVector_CreateFinish(rhsVector%previousRHSVector,err,error,*999)
        CALL DistributedVector_AllValuesSet(rhsVector%previousRHSVector,0.0_DP,err,error,*999)
      CASE(EQUATIONS_MATRICES_PREVIOUS2_VECTOR)
        NULLIFY(rhsVector%previous2RHSVector)
        CALL DistributedVector_CreateStart(rowsDomainMap,rhsVector%previous2RHSVector,err,error,*999)
        CALL DistributedVector_DataTypeSet(rhsVector%previous2RHSVector,MATRIX_VECTOR_DP_TYPE,err,error,*999)
        CALL DistributedVector_CreateFinish(rhsVector%previous2RHSVector,err,error,*999)
        CALL DistributedVector_AllValuesSet(rhsVector%previous2RHSVector,0.0_DP,err,error,*999)
      CASE(EQUATIONS_MATRICES_PREVIOUS3_VECTOR)
        NULLIFY(rhsVector%previous3RHSVector)
        CALL DistributedVector_CreateStart(rowsDomainMap,rhsVector%previous3RHSVector,err,error,*999)
        CALL DistributedVector_DataTypeSet(rhsVector%previous3RHSVector,MATRIX_VECTOR_DP_TYPE,err,error,*999)
        CALL DistributedVector_CreateFinish(rhsVector%previous3RHSVector,err,error,*999)
        CALL DistributedVector_AllValuesSet(rhsVector%previous3RHSVector,0.0_DP,err,error,*999)
      CASE DEFAULT
        localError="The specified RHS vector temporal type of "//TRIM(NumberToVString(temporalType,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDIF
    
    EXITS("EquationsMatricesRHS_DistributedVectorEnsureCreated")
    RETURN
999 ERRORS("EquationsMatricesRHS_DistributedVectorEnsureCreated",err,error)
    EXITS("EquationsMatricesRHS_DistributedVectorEnsureCreated")
    RETURN 1
    
  END SUBROUTINE EquationsMatricesRHS_DistributedVectorEnsureCreated

  !
  !================================================================================================================================
  !

  !>Ensure that a particular temporal type distributed vector has been created for a equations matrices source vector.
  SUBROUTINE EquationsMatricesSource_DistributedVectorEnsureCreated(sourceVector,temporalType,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector !<The source vector to ensure that the distributed vector is created for.
    INTEGER(INTG), INTENT(IN) :: temporalType !<The equations matrices vector temporal type to ensure that it has been created. \see EquationsMatricesRoutines_VectorTemporalTypes
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DomainMappingType), POINTER :: rowsDomainMap
    TYPE(DistributedVectorType), POINTER :: sourceDistributedVector
    TYPE(EquationsMatricesSourcesType), POINTER :: sourceVectors
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(FieldVariableType), POINTER :: lhsVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsMatricesSource_DistributedVectorEnsureCreated",err,error,*999)

    NULLIFY(sourceDistributedVector)
    CALL EquationsMatricesSource_DistributedVectorExists(sourceVector,temporalType,sourceDistributedVector,err,error,*999)
    IF(.NOT.ASSOCIATED(sourceDistributedVector)) THEN
      NULLIFY(sourceVectors)
      CALL EquationsMatricesSource_SourceVectorsGet(sourceVector,sourceVectors,err,error,*999)
      NULLIFY(vectorMatrices)
      CALL EquationsMatricesSources_VectorMatricesGet(sourceVectors,vectorMatrices,err,error,*999)
      NULLIFY(vectorMapping)
      CALL EquationsMatricesVector_VectorMappingGet(vectorMatrices,vectorMapping,err,error,*999)
      NULLIFY(lhsMapping)
      CALL EquationsMappingVector_LHSMappingGet(vectorMapping,lhsMapping,err,error,*999)
      NULLIFY(lhsVariable)
      CALL EquationsMappingLHS_LHSVariableGet(lhsMapping,lhsVariable,err,error,*999)
      NULLIFY(rowsDomainMap)
      CALL FieldVariable_DomainMappingGet(lhsVariable,rowsDomainMap,err,error,*999)
      SELECT CASE(temporalType)
      CASE(EQUATIONS_MATRICES_CURRENT_VECTOR)
        NULLIFY(sourceVector%vector)
        CALL DistributedVector_CreateStart(rowsDomainMap,sourceVector%vector,err,error,*999)
        CALL DistributedVector_DataTypeSet(sourceVector%vector,MATRIX_VECTOR_DP_TYPE,err,error,*999)
        CALL DistributedVector_CreateFinish(sourceVector%vector,err,error,*999)
        CALL DistributedVector_AllValuesSet(sourceVector%vector,0.0_DP,err,error,*999)
      CASE(EQUATIONS_MATRICES_PREVIOUS_VECTOR)
        NULLIFY(sourceVector%previousSourceVector)
        CALL DistributedVector_CreateStart(rowsDomainMap,sourceVector%previousSourceVector,err,error,*999)
        CALL DistributedVector_DataTypeSet(sourceVector%previousSourceVector,MATRIX_VECTOR_DP_TYPE,err,error,*999)
        CALL DistributedVector_CreateFinish(sourceVector%previousSourceVector,err,error,*999)
        CALL DistributedVector_AllValuesSet(sourceVector%previousSourceVector,0.0_DP,err,error,*999)
      CASE(EQUATIONS_MATRICES_PREVIOUS2_VECTOR)
        NULLIFY(sourceVector%previous2SourceVector)
        CALL DistributedVector_CreateStart(rowsDomainMap,sourceVector%previous2SourceVector,err,error,*999)
        CALL DistributedVector_DataTypeSet(sourceVector%previous2SourceVector,MATRIX_VECTOR_DP_TYPE,err,error,*999)
        CALL DistributedVector_CreateFinish(sourceVector%previous2SourceVector,err,error,*999)
        CALL DistributedVector_AllValuesSet(sourceVector%previous2SourceVector,0.0_DP,err,error,*999)
      CASE(EQUATIONS_MATRICES_PREVIOUS3_VECTOR)
        NULLIFY(sourceVector%previous3SourceVector)
        CALL DistributedVector_CreateStart(rowsDomainMap,sourceVector%previous3SourceVector,err,error,*999)
        CALL DistributedVector_DataTypeSet(sourceVector%previous3SourceVector,MATRIX_VECTOR_DP_TYPE,err,error,*999)
        CALL DistributedVector_CreateFinish(sourceVector%previous3SourceVector,err,error,*999)
        CALL DistributedVector_AllValuesSet(sourceVector%previous3SourceVector,0.0_DP,err,error,*999)
      CASE DEFAULT
        localError="The specified RHS vector temporal type of "//TRIM(NumberToVString(temporalType,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDIF
    
    EXITS("EquationsMatricesSource_DistributedVectorEnsureCreated")
    RETURN
999 ERRORS("EquationsMatricesSource_DistributedVectorEnsureCreated",err,error)
    EXITS("EquationsMatricesSource_DistributedVectorEnsureCreated")
    RETURN 1
    
  END SUBROUTINE EquationsMatricesSource_DistributedVectorEnsureCreated

  !
  !================================================================================================================================
  !

  !>Adds the element matrices and rhs vector into the vector equations matrices and rhs vector.
  SUBROUTINE EquationsMatricesVector_ElementAdd(vectorMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the vector equations matrices
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: columnIdx,matrixIdx,residualIdx,rowIdx,sourceIdx
    REAL(DP) :: sum
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatricesSourcesType), POINTER :: sourceVectors
    TYPE(EquationsMatrixType), POINTER :: dynamicMatrix,linearMatrix

    ENTERS("EquationsMatricesVector_ElementAdd",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector equations matrices is not allocated.",err,error,*999)

    NULLIFY(dynamicMatrices)
    CALL EquationsMatricesVector_DynamicMatricesExists(vectorMatrices,dynamicMatrices,err,error,*999)
    IF(ASSOCIATED(dynamicMatrices)) THEN
      !Add the element matrices
      DO matrixIdx=1,dynamicMatrices%numberOfDynamicMatrices
        NULLIFY(dynamicMatrix)
        CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,matrixIdx,dynamicMatrix,err,error,*999)
        IF(dynamicMatrix%updateMatrix) THEN
          !Handle lumped matrices
          IF(dynamicMatrix%lumped) THEN
            DO rowIdx=1,dynamicMatrix%elementMatrix%numberOfRows
              sum=0.0_DP
              DO columnIdx=1,dynamicMatrix%elementMatrix%numberOfColumns
                sum=sum+dynamicMatrix%elementMatrix%matrix(rowIdx,columnIdx)
                dynamicMatrix%elementMatrix%matrix(rowIdx,columnIdx)=0.0_DP
              ENDDO !columnIdx
              dynamicMatrix%elementMatrix%matrix(rowIdx,rowIdx)=sum
              !Add the element matrices into the distributed equations matrix
              CALL DistributedMatrix_ValuesAdd(dynamicMatrix%matrix,dynamicMatrix%elementMatrix%rowDOFS(rowIdx), &
                & dynamicMatrix%elementMatrix%columnDOFS(rowIdx),dynamicMatrix%elementMatrix%matrix(rowIdx,rowIdx), &
                & err,error,*999)
            ENDDO !rowIdx
          ELSE
            !Add the element matrice into the distributed equations matrix
            CALL DistributedMatrix_ValuesAdd(dynamicMatrix%matrix,dynamicMatrix%elementMatrix%rowDOFS(1: &
              & dynamicMatrix%elementMatrix%numberOfRows),dynamicMatrix%elementMatrix%columnDOFS(1: &
              & dynamicMatrix%elementMatrix%numberOfColumns),dynamicMatrix%elementMatrix%matrix(1: &
              & dynamicMatrix%elementMatrix%numberOfRows,1:dynamicMatrix%elementMatrix%numberOfColumns), &
              & err,error,*999)
          ENDIF
        ENDIF
      ENDDO !matrixIdx
    ENDIF
    NULLIFY(linearMatrices)
    CALL EquationsMatricesVector_LinearMatricesExists(vectorMatrices,linearMatrices,err,error,*999)
    IF(ASSOCIATED(linearMatrices)) THEN
      !Add the element matrices
      DO matrixIdx=1,linearMatrices%numberOfLinearMatrices
        NULLIFY(linearMatrix)
        CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,matrixIdx,linearMatrix,err,error,*999)
        IF(linearMatrix%updateMatrix) THEN
          !Handle lumped matrices
          IF(linearMatrix%lumped) THEN
            DO rowIdx=1,linearMatrix%elementMatrix%numberOfRows
              sum=0.0_DP
              DO columnIdx=1,linearMatrix%elementMatrix%numberOfColumns
                sum=sum+linearMatrix%elementMatrix%matrix(rowIdx,columnIdx)
                linearMatrix%elementMatrix%matrix(rowIdx,columnIdx)=0.0_DP
              ENDDO !columnIdx
              linearMatrix%elementMatrix%matrix(rowIdx,rowIdx)=sum
              !Add the element matrice into the distributed equations matrix
              CALL DistributedMatrix_ValuesAdd(linearMatrix%matrix,linearMatrix%elementMatrix%rowDOFS(rowIdx), &
                & linearMatrix%elementMatrix%columnDOFS(rowIdx),linearMatrix%elementMatrix%matrix(rowIdx, &
                & rowIdx),err,error,*999)
            ENDDO !rowIdx
          ELSE
            !Add the element matrice into the distributed equations matrix
            CALL DistributedMatrix_ValuesAdd(linearMatrix%matrix,linearMatrix%elementMatrix%rowDOFS(1: &
              & linearMatrix%elementMatrix%numberOfRows),linearMatrix%elementMatrix%columnDOFS(1: &
              & linearMatrix%elementMatrix%numberOfColumns),linearMatrix%elementMatrix%matrix(1: &
              & linearMatrix%elementMatrix%numberOfRows,1:linearMatrix%elementMatrix%numberOfColumns), &
              & err,error,*999)
          ENDIF
        ENDIF
      ENDDO !matrixIdx
    ENDIF
    NULLIFY(nonlinearMatrices)
    CALL EquationsMatricesVector_NonlinearMatricesExists(vectorMatrices,nonlinearMatrices,err,error,*999)
    IF(ASSOCIATED(nonlinearMatrices)) THEN
      DO residualIdx=1,nonlinearMatrices%numberOfResiduals
        NULLIFY(residualVector)
        CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,residualIdx,residualVector,err,error,*999)
        IF(residualVector%updateResidual) THEN
          !Add the residual element vector
          CALL DistributedVector_ValuesAdd(residualVector%residual,residualVector%elementResidual%rowDOFS(1: &
            & residualVector%elementResidual%numberOfRows),residualVector%elementResidual%vector(1:residualVector% &
            & elementResidual%numberOfRows),err,error,*999)
        ENDIF
      ENDDO !residualIdx
    ENDIF
    NULLIFY(rhsVector)
    CALL EquationsMatricesVector_RHSVectorExists(vectorMatrices,rhsVector,err,error,*999)
    IF(ASSOCIATED(rhsVector)) THEN
      IF(rhsVector%updateVector) THEN
        !Add the rhs element vector
        CALL DistributedVector_ValuesAdd(rhsVector%vector,rhsVector%elementVector%rowDOFS(1: &
          & rhsVector%elementVector%numberOfRows),rhsVector%elementVector%vector(1:rhsVector% &
          & elementVector%numberOfRows),err,error,*999)
      ENDIF
    ENDIF
    NULLIFY(sourceVectors)
    CALL EquationsMatricesVector_SourceVectorsExists(vectorMatrices,sourceVectors,err,error,*999)
    IF(ASSOCIATED(sourceVectors)) THEN
      DO sourceIdx=1,sourceVectors%numberOfSources
        NULLIFY(sourceVector)
        CALL EquationsMatricesSources_SourceVectorGet(sourceVectors,sourceIdx,sourceVector,err,error,*999)
        IF(sourceVector%updateVector) THEN
          !Add the source element vector
          CALL DistributedVector_ValuesAdd(sourceVector%vector,sourceVector%elementVector%rowDOFS(1: &
            & sourceVector%elementVector%numberOfRows),sourceVector%elementVector%vector(1:sourceVector% &
            & elementVector%numberOfRows),err,error,*999)
        ENDIF
      ENDDO !sourceIdx
    ENDIF
   
    EXITS("EquationsMatricesVector_ElementAdd")
    RETURN
999 ERRORSEXITS("EquationsMatricesVector_ElementAdd",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_ElementAdd

  !
  !================================================================================================================================
  !

  !>Calculate the positions in the vector equations matrices and rhs of the element matrices and rhs vector. Old CMISS name MELGE.
  SUBROUTINE EquationsMatricesVector_ElementCalculate(vectorMatrices,elementNumber,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the vector equations matrices
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to calculate the mappings for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx,residualIdx,sourceIdx
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatricesSourcesType), POINTER :: sourceVectors
    TYPE(EquationsMatrixType), POINTER :: dynamicMatrix,linearMatrix
    TYPE(FieldVariableType), POINTER :: colVariable,rowVariable
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix

    ENTERS("EquationsMatricesVector_ElementCalculate",err,error,*999)

    NULLIFY(vectorMapping)
    CALL EquationsMatricesVector_VectorMappingGet(vectorMatrices,vectorMapping,err,error,*999)
    NULLIFY(lhsMapping)
    CALL EquationsMappingVector_LHSMappingGet(vectorMapping,lhsMapping,err,error,*999)
    NULLIFY(rowVariable)
    CALL EquationsMappingLHS_LHSVariableGet(lhsMapping,rowVariable,err,error,*999)
    NULLIFY(dynamicMatrices)
    CALL EquationsMatricesVector_DynamicMatricesExists(vectorMatrices,dynamicMatrices,err,error,*999)
    IF(ASSOCIATED(dynamicMatrices)) THEN
      !Calculate the row and columns for the dynamic equations matrices
      NULLIFY(dynamicMapping)
      CALL EquationsMappingVector_DynamicMappingGet(vectorMapping,dynamicMapping,err,error,*999)
      NULLIFY(colVariable)
      CALL EquationsMappingDynamic_DynamicVariableGet(dynamicMapping,colVariable,err,error,*999)
      DO matrixIdx=1,dynamicMatrices%numberOfDynamicMatrices
        NULLIFY(dynamicMatrix)
        CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,matrixIdx,dynamicMatrix,err,error,*999)
        CALL EquationsMatrices_ElementMatrixCalculate(dynamicMatrix%elementMatrix,dynamicMatrix%updateMatrix, &
          & [elementNumber],[elementNumber],rowVariable,colVariable,err,error,*999)
      ENDDO !matrixIdx
    ENDIF
    NULLIFY(linearMatrices)
    CALL EquationsMatricesVector_LinearMatricesExists(vectorMatrices,linearMatrices,err,error,*999)
    IF(ASSOCIATED(linearMatrices)) THEN
      !Calculate the row and columns for the linear equations matrices
      NULLIFY(linearMapping)
      CALL EquationsMappingVector_LinearMappingGet(vectorMapping,linearMapping,err,error,*999)
      DO matrixIdx=1,linearMatrices%numberOfLinearMatrices
        NULLIFY(linearMatrix)
        CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,matrixIdx,linearMatrix,err,error,*999)
        NULLIFY(colVariable)
        CALL EquationsMappingLinear_LinearMatrixVariableGet(linearMapping,matrixIdx,colVariable,err,error,*999)
        CALL EquationsMatrices_ElementMatrixCalculate(linearMatrix%elementMatrix,linearMatrix%updateMatrix, &
          & [elementNumber],[elementNumber],rowVariable,colVariable,err,error,*999)
      ENDDO !matrixIdx
    ENDIF
    NULLIFY(nonlinearMatrices)
    CALL EquationsMatricesVector_NonlinearMatricesExists(vectorMatrices,nonlinearMatrices,err,error,*999)
    IF(ASSOCIATED(nonlinearMatrices)) THEN
      NULLIFY(nonlinearMapping)
      CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)
      DO residualIdx=1,nonlinearMapping%numberOfResiduals
        NULLIFY(residualVector)
        CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,residualIdx,residualVector,err,error,*999)
        NULLIFY(residualMapping)
        CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,residualIdx,residualMapping,err,error,*999)
        !Calculate the rows of the equations residual
        CALL EquationsMatrices_ElementVectorCalculate(residualVector%elementResidual,residualVector%updateResidual, &
          & elementNumber,rowVariable,err,error,*999)
        residualVector%elementResidualCalculated=0
        DO matrixIdx=1,residualVector%numberOfJacobians
          !Calculate the rows and columns of the Jacobian
          NULLIFY(jacobianMatrix)
          CALL EquationsMatricesResidual_JacobianMatrixGet(residualVector,matrixIdx,jacobianMatrix,err,error,*999)
          NULLIFY(colVariable)
          CALL EquationsMappingResidual_JacobianMatrixVariableGet(residualMapping,matrixIdx,colVariable,err,error,*999)
          CALL EquationsMatrices_ElementMatrixCalculate(jacobianMatrix%elementJacobian,jacobianMatrix%updateJacobian, &
            & [elementNumber],[elementNumber],rowVariable,colVariable,err,error,*999)
        ENDDO !matrixIdx
      ENDDO !residualIdx
    ENDIF
    NULLIFY(rhsVector)
    CALL EquationsMatricesVector_RHSVectorExists(vectorMatrices,rhsVector,err,error,*999)
    IF(ASSOCIATED(rhsVector)) THEN
      !Calculate the rows for the vector equations RHS
      CALL EquationsMatrices_ElementVectorCalculate(rhsVector%elementVector,rhsVector%updateVector,elementNumber,rowVariable, &
        & err,error,*999)
    ENDIF
    NULLIFY(sourceVectors)
    CALL EquationsMatricesVector_SourceVectorsExists(vectorMatrices,sourceVectors,err,error,*999)
    IF(ASSOCIATED(sourceVectors)) THEN
      DO sourceIdx=1,sourceVectors%numberOfSources
        NULLIFY(sourceVector)
        CALL EquationsMatricesSources_SourceVectorGet(sourceVectors,sourceIdx,sourceVector,err,error,*999)
        !Calculate the rows the equations source.
        CALL EquationsMatrices_ElementVectorCalculate(sourceVector%elementVector,sourceVector%updateVector,elementNumber, &
          & rowVariable,err,error,*999)
      ENDDO !sourceIdx
    ENDIF
    
    EXITS("EquationsMatricesVector_ElementCalculate")
    RETURN
999 ERRORSEXITS("EquationsMatricesVector_ElementCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_ElementCalculate

  !
  !================================================================================================================================
  !

  !>Calculate the positions in the vector equations matrices and rhs of the nodal matrices and rhs vector. Old CMISS name MELGE.
  SUBROUTINE EquationsMatricesVector_NodalCalculate(vectorMatrices,nodeNumber,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the vector equations matrices
    INTEGER(INTG), INTENT(IN) :: nodeNumber !<The nodal number to calculate the mappings for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx,residualIdx,sourceIdx
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatricesSourcesType), POINTER :: sourceVectors
    TYPE(EquationsMatrixType), POINTER :: dynamicMatrix,linearMatrix
    TYPE(FieldVariableType), POINTER :: colVariable,rowVariable
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix

    ENTERS("EquationsMatricesVector_NodalCalculate",err,error,*999)

    NULLIFY(vectorMapping)
    CALL EquationsMatricesVector_VectorMappingGet(vectorMatrices,vectorMapping,err,error,*999)
    NULLIFY(lhsMapping)
    CALL EquationsMappingVector_LHSMappingGet(vectorMapping,lhsMapping,err,error,*999)
    NULLIFY(rowVariable)
    CALL EquationsMappingLHS_LHSVariableGet(lhsMapping,rowVariable,err,error,*999)
    NULLIFY(dynamicMatrices)
    CALL EquationsMatricesVector_DynamicMatricesExists(vectorMatrices,dynamicMatrices,err,error,*999)
    IF(ASSOCIATED(dynamicMatrices)) THEN
      !Calculate the row and columns for the dynamic equations matrices
      NULLIFY(dynamicMapping)
      CALL EquationsMappingVector_DynamicMappingGet(vectorMapping,dynamicMapping,err,error,*999)
      NULLIFY(colVariable)
      CALL EquationsMappingDynamic_DynamicVariableGet(dynamicMapping,colVariable,err,error,*999)
      DO matrixIdx=1,dynamicMatrices%numberOfDynamicMatrices
        NULLIFY(dynamicMatrix)
        CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,matrixIdx,dynamicMatrix,err,error,*999)
        CALL EquationsMatrices_NodalMatrixCalculate(dynamicMatrix%nodalMatrix,dynamicMatrix%updateMatrix, &
          & nodeNumber,nodeNumber,rowVariable,colVariable,err,error,*999)
      ENDDO !matrixIdx
    ENDIF
    NULLIFY(linearMatrices)
    CALL EquationsMatricesVector_LinearMatricesExists(vectorMatrices,linearMatrices,err,error,*999)
    IF(ASSOCIATED(linearMatrices)) THEN
      !Calculate the row and columns for the linear equations matrices
      NULLIFY(linearMapping)
      CALL EquationsMappingVector_LinearMappingGet(vectorMapping,linearMapping,err,error,*999)
      DO matrixIdx=1,linearMatrices%numberOfLinearMatrices
        NULLIFY(linearMatrix)
        CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,matrixIdx,linearMatrix,err,error,*999)
        NULLIFY(colVariable)
        CALL EquationsMappingLinear_LinearMatrixVariableGet(linearMapping,matrixIdx,colVariable,err,error,*999)
        CALL EquationsMatrices_NodalMatrixCalculate(linearMatrix%nodalMatrix,linearMatrix%updateMatrix, &
          & nodeNumber,nodeNumber,rowVariable,colVariable,err,error,*999)
      ENDDO !matrixIdx
    ENDIF
    NULLIFY(nonlinearMatrices)
    CALL EquationsMatricesVector_NonlinearMatricesExists(vectorMatrices,nonlinearMatrices,err,error,*999)
    IF(ASSOCIATED(nonlinearMatrices)) THEN
      !Calculate the rows and columns of the Jacobian
      NULLIFY(nonlinearMapping)
      CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)
      DO residualIdx=1,nonlinearMapping%numberOfResiduals
        NULLIFY(residualVector)
        CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,residualIdx,residualVector,err,error,*999)
        NULLIFY(residualMapping)
        CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,residualIdx,residualMapping,err,error,*999)
        !Calculate the rows of the equations residual
        CALL EquationsMatrices_NodalVectorCalculate(residualVector%nodalResidual,residualVector%updateResidual, &
          & nodeNumber,rowVariable,err,error,*999)
        residualVector%nodalResidualCalculated=0
        DO matrixIdx=1,residualVector%numberOfJacobians
          NULLIFY(jacobianMatrix)
          CALL EquationsMatricesResidual_JacobianMatrixGet(residualVector,matrixIdx,jacobianMatrix,err,error,*999)
          NULLIFY(colVariable)
          CALL EquationsMappingResidual_JacobianMatrixVariableGet(residualMapping,matrixIdx,colVariable,err,error,*999)
          CALL EquationsMatrices_NodalMatrixCalculate(jacobianMatrix%nodalJacobian,jacobianMatrix%updateJacobian, &
            & nodeNumber,nodeNumber,rowVariable,colVariable,err,error,*999)
        ENDDO !matrixIdx
      ENDDO !residualIdx
    ENDIF
    NULLIFY(rhsVector)
    CALL EquationsMatricesVector_RHSVectorExists(vectorMatrices,rhsVector,err,error,*999)
    IF(ASSOCIATED(rhsVector)) THEN
      !Calculate the rows for the equations RHS
      CALL EquationsMatrices_NodalVectorCalculate(rhsVector%nodalVector,rhsVector%updateVector,nodeNumber, &
        & rowVariable,err,error,*999)
    ENDIF
    NULLIFY(sourceVectors)
    CALL EquationsMatricesVector_SourceVectorExists(vectorMatrices,sourceVectors,err,error,*999)
    IF(ASSOCIATED(sourceVectors)) THEN
      DO sourceIdx=1,sourceVectors%numberOfSources
        NULLIFY(sourceVector)
        CALL EquationsMatricesSources_SourceVectorGet(sourceVectors,sourceIdx,sourceVector,err,error,*999)
        !Calculate the rows the equations source.
        CALL EquationsMatrices_NodalVectorCalculate(sourceVector%nodalVector,sourceVector%updateVector,nodeNumber, &
          & rowVariable,err,error,*999)
      ENDDO !sourceIdx
    ENDIF
    
    EXITS("EquationsMatricesVector_NodalCalculate")
    RETURN
999 ERRORSEXITS("EquationsMatricesVector_NodalCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_NodalCalculate

  !
  !================================================================================================================================
  !

  !>Calculate the positions in the equations matrices of the nodal matrix.
  SUBROUTINE EquationsMatrices_NodalMatrixCalculate(nodalMatrix,updateMatrix,rowNodeNumber,columnNodeNumber, &
    & rowsFieldVariable,colsFieldVariable,err,error,*)

    !Argument variables
    TYPE(NodalMatrixType) :: nodalMatrix !<The nodal matrix to calculate
    LOGICAL :: updateMatrix !<Is .TRUE. if the nodal matrix is to be updated, .FALSE. if not.
    INTEGER(INTG), INTENT(IN) :: rowNodeNumber !<The row nodal number to calculate
    INTEGER(INTG), INTENT(IN) :: columnNodeNumber !<The column nodal number to calculate
    TYPE(FieldVariableType), POINTER :: rowsFieldVariable !<A pointer to the field variable associated with the rows
    TYPE(FieldVariableType), POINTER :: colsFieldVariable !<A pointer to the field variable associated with the columns
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: colsVariableType,componentIdx,derivativeIdx,globalColumn,globalRow,interpolationType,localColumn,localRow, &
      & numberOfComponents,numberOfDerivatives,numberOfVersions,rowsVariableType,totalNumberOfNodes,versionIdx
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionDataPointsType), POINTER :: dataPoints
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainMappingType), POINTER :: domainMapping
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsMatrices_NodalMatrixCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(rowsFieldVariable)) CALL FlagError("Rows field variable is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(colsFieldVariable)) CALL FlagError("Columns field variable is not associated.",err,error,*999)
    CALL FieldVariable_VariableTypeGet(rowsFieldVariable,rowsVariableType,err,error,*999)
    CALL FieldVariable_VariableTypeGet(colsFieldVariable,colsVariableType,err,error,*999)
    
    nodalMatrix%numberOfRows=0
    nodalMatrix%numberOfColumns=0
    IF(updateMatrix) THEN
      IF(ASSOCIATED(rowsFieldVariable,colsFieldVariable)) THEN
        !Row and columns variable is the same.
        CALL FieldVariable_NumberOfComponentsGet(rowsFieldVariable,numberOfComponents,err,error,*999)
        DO componentIdx=1,numberOfComponents
          NULLIFY(domain)
          CALL FieldVariable_ComponentDomainGet(rowsFieldVariable,componentIdx,domain,err,error,*999)
          NULLIFY(domainTopology)
          CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
          NULLIFY(domainNodes)
          CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
          CALL FieldVariable_ComponentInterpolationGet(rowsFieldVariable,componentIdx,interpolationType,err,error,*999)
          NULLIFY(domainMapping)
          CALL FieldVariable_DomainMappingGet(rowsFieldVariable,domainMapping,err,error,*999)
          CALL DomainNodes_TotalNumberOfNodesGet(domainNodes,totalNumberOfNodes,err,error,*999)
          IF(rowNodeNumber<1.OR.rowNodeNumber>totalNumberOfNodes) THEN
            localError="Nodal number "//TRIM(NumberToVString(rowNodeNumber,"*",err,error))// &
              & " is invalid for component number "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
              & " of rows field variable type "//TRIM(NumberToVString(rowsVariableType,"*",err,error))// &
              & ". The nodal number must be >= 1 and <= "// &
              & TRIM(NumberToVString(totalNumberOfNodes,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
          SELECT CASE(interpolationType)
          CASE(FIELD_CONSTANT_INTERPOLATION)
            CALL FieldVariable_ConstantDOFGet(rowsFieldVariable,componentIdx,localRow,err,error,*999)
            CALL DomainMapping_LocalToGlobalGet(domainMapping,localRow,globalColumn,err,error,*999)
            nodalMatrix%numberOfRows=nodalMatrix%numberOfRows+1
            nodalMatrix%numberOfColumns=nodalMatrix%numberOfColumns+1
            nodalMatrix%rowDofs(nodalMatrix%numberOfRows)=localRow
            nodalMatrix%columnDofs(nodalMatrix%numberOfColumns)=globalColumn
          CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
            CALL FieldVariable_LocalElementDOFGet(rowsFieldVariable,rowNodeNumber,componentIdx,localRow,err,error,*999)
            CALL DomainMapping_LocalToGlobalGet(domainMapping,localRow,globalColumn,err,error,*999)
            nodalMatrix%numberOfRows=nodalMatrix%numberOfRows+1
            nodalMatrix%numberOfColumns=nodalMatrix%numberOfColumns+1
            nodalMatrix%rowDofs(nodalMatrix%numberOfRows)=localRow
            nodalMatrix%columnDofs(nodalMatrix%numberOfColumns)=globalColumn
          CASE(FIELD_NODE_BASED_INTERPOLATION)
            CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,rowNodeNumber,numberOfDerivatives,err,error,*999)
            DO derivativeIdx=1,numberOfDerivatives
              CALL DomainNodes_DerivativeNumberOfVersionsGet(domainNodes,derivativeIdx,rowNodeNumber,numberOfVersions, &
                & err,error,*999)
              DO versionIdx=1,numberOfVersions
                CALL FieldVariable_LocalNodeDOFGet(rowsFieldVariable,versionIdx,derivativeIdx,rowNodeNumber,componentIdx, &
                  & localRow,err,error,*999)
                CALL DomainMapping_LocalToGlobalGet(domainMapping,localRow,globalColumn,err,error,*999)
                nodalMatrix%numberOfRows=nodalMatrix%numberOfRows+1
                nodalMatrix%numberOfColumns=nodalMatrix%numberOfColumns+1
                nodalMatrix%rowDofs(nodalMatrix%numberOfRows)=localRow
                nodalMatrix%columnDofs(nodalMatrix%numberOfColumns)=globalColumn
              ENDDO !versionIdx
            ENDDO !derivativeIdx
          CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The interpolation type of "// &
              & TRIM(NumberToVString(interpolationType,"*",err,error))// &
              & " is invalid for component number "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
              & " of rows field variable type "//TRIM(NumberToVString(rowsVariableType,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)          
          END SELECT
        ENDDO !componentIdx
      ELSE
        !Row and column variables are different
        !Row mapping
        CALL FieldVariable_NumberOfComponentsGet(rowsFieldVariable,numberOfComponents,err,error,*999)
        DO componentIdx=1,numberOfComponents
          NULLIFY(domain)
          CALL FieldVariable_ComponentDomainGet(rowsFieldVariable,componentIdx,domain,err,error,*999)
          NULLIFY(domainTopology)
          CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
          NULLIFY(domainNodes)
          CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
          CALL FieldVariable_ComponentInterpolationGet(rowsFieldVariable,componentIdx,interpolationType,err,error,*999)
          CALL DomainNodes_TotalNumberOfNodesGet(domainNodes,totalNumberOfNodes,err,error,*999)
          IF(rowNodeNumber<1.OR.rowNodeNumber>totalNumberOfNodes) THEN
            localError="Row nodal number "//TRIM(NumberToVString(rowNodeNumber,"*",err,error))// &
              & " is invalid for component number "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
              & " of rows field variable type "//TRIM(NumberToVString(rowsVariableType,"*",err,error))// &
              & ". The nodal number must be between 1 and "//TRIM(NumberToVString(totalNumberOfNodes,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
          SELECT CASE(interpolationType)
          CASE(FIELD_CONSTANT_INTERPOLATION)
            CALL FieldVariable_ConstantDOFGet(rowsFieldVariable,componentIdx,localRow,err,error,*999)
            nodalMatrix%numberOfRows=nodalMatrix%numberOfRows+1
            nodalMatrix%rowDofs(nodalMatrix%numberOfRows)=localRow
          CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
!!TODO: why is the node number being used as an element number here?            
            CALL FieldVariable_LocalElementDOFGet(rowsFieldVariable,rowNodeNumber,componentIdx,localRow,err,error,*999)
            nodalMatrix%numberOfRows=nodalMatrix%numberOfRows+1
            nodalMatrix%rowDofs(nodalMatrix%numberOfRows)=localRow
          CASE(FIELD_NODE_BASED_INTERPOLATION)
            CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,rowNodeNumber,numberOfDerivatives,err,error,*999)
            DO derivativeIdx=1,numberOfDerivatives
              CALL DomainNodes_DerivativeNumberOfVersionsGet(domainNodes,derivativeIdx,rowNodeNumber,numberOfVersions, &
                & err,error,*999)
              DO versionIdx=1,numberOfVersions
                CALL FieldVariable_LocalNodeDOFGet(rowsFieldVariable,versionIdx,derivativeIdx,rowNodeNumber,componentIdx, &
                  & localRow,err,error,*999)
                nodalMatrix%numberOfRows=nodalMatrix%numberOfRows+1
                nodalMatrix%rowDofs(nodalMatrix%numberOfRows)=localRow
              ENDDO !versionIdx
            ENDDO !derivativeIdx
          CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The interpolation type of "// &
              & TRIM(NumberToVString(interpolationType,"*",err,error))// &
              & " is invalid for component number "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
              & " of rows field variable type "//TRIM(NumberToVString(rowsVariableType,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)          
          END SELECT
        ENDDO !componentIdx
        !Column mapping
        CALL FieldVariable_NumberOfComponentsGet(colsFieldVariable,numberOfComponents,err,error,*999)
        DO componentIdx=1,numberOfComponents
          NULLIFY(domain)
          CALL FieldVariable_ComponentDomainGet(colsFieldVariable,componentIdx,domain,err,error,*999)
          NULLIFY(domainTopology)
          CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
          NULLIFY(domainNodes)
          CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
          CALL FieldVariable_ComponentInterpolationGet(colsFieldVariable,componentIdx,interpolationType,err,error,*999)
          NULLIFY(domainMapping)
          CALL FieldVariable_DomainMappingGet(colsFieldVariable,domainMapping,err,error,*999)
          CALL DomainNodes_TotalNumberOfNodesGet(domainNodes,totalNumberOfNodes,err,error,*999)
          IF(columnNodeNumber<1.OR.columnNodeNumber>totalNumberOfNodes) THEN
            localError="Column nodal number "//TRIM(NumberToVString(columnNodeNumber,"*",err,error))// &
              & " is invalid for component number "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
              & " of column field variable type "//TRIM(NumberToVString(colsVariableType,"*",err,error))// &
              & ". The nodal number must be between 1 and "//TRIM(NumberToVString(totalNumberOfNodes,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
          SELECT CASE(interpolationType)
          CASE(FIELD_CONSTANT_INTERPOLATION)
            CALL FieldVariable_ConstantDOFGet(colsFieldVariable,componentIdx,localColumn,err,error,*999)
            CALL DomainMapping_LocalToGlobalGet(domainMapping,localColumn,globalColumn,err,error,*999)
            nodalMatrix%numberOfColumns=nodalMatrix%numberOfColumns+1
            nodalMatrix%columnDofs(nodalMatrix%numberOfColumns)=globalColumn
          CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
!!TODO: why is a node number being used as an element number here?            
            CALL FieldVariable_LocalElementDOFGet(colsFieldVariable,columnNodeNumber,componentIdx,localColumn,err,error,*999)
            CALL DomainMapping_LocalToGlobalGet(domainMapping,localColumn,globalColumn,err,error,*999)
            nodalMatrix%numberOfColumns=nodalMatrix%numberOfColumns+1
            nodalMatrix%columnDofs(nodalMatrix%numberOfColumns)=globalColumn
          CASE(FIELD_NODE_BASED_INTERPOLATION)
            CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,columnNodeNumber,numberOfDerivatives,err,error,*999)
            DO derivativeIdx=1,numberOfDerivatives
              CALL DomainNodes_DerivativeNumberOfVersionsGet(domainNodes,derivativeIdx,columnNodeNumber,numberOfVersions, &
                & err,error,*999)
              DO versionIdx=1,numberOfVersions
                CALL FieldVariable_LocalNodeDOFGet(colsFieldVariable,versionIdx,derivativeIdx,columnNodeNumber,componentIdx, &
                  & localColumn,err,error,*999)
                CALL DomainMapping_LocalToGlobalGet(domainMapping,localColumn,globalColumn,err,error,*999)
                nodalMatrix%numberOfColumns=nodalMatrix%numberOfColumns+1
                nodalMatrix%columnDofs(nodalMatrix%numberOfColumns)=globalColumn
              ENDDO !versionIdx
            ENDDO !derivativeIdx
          CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The interpolation type of "// &
              & TRIM(NumberToVString(interpolationType,"*",err,error))// &
              & " is invalid for component number "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
              & " of column field variable type "//TRIM(NumberToVString(colsVariableType,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)          
          END SELECT
        ENDDO !componentIdx
      ENDIF
      nodalMatrix%matrix=0.0_DP
    ENDIF
    
    EXITS("EquationsMatrices_NodalMatrixCalculate")
    RETURN
999 ERRORSEXITS("EquationsMatrices_NodalMatrixCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_NodalMatrixCalculate
  
  !
  !================================================================================================================================
  !

  !>Calculate the positions in the vector equations rhs of the nodal rhs vector.
  SUBROUTINE EquationsMatrices_NodalVectorCalculate(nodalVector,updateVector,rowNodeNumber,rowsFieldVariable,err,error,*)

    !Argument variables
    TYPE(NodalVectorType) :: nodalVector !<The nodal vector to calculate.
    LOGICAL :: updateVector !<Is .TRUE. if the nodal vector is to be updated, .FALSE. if not.
    INTEGER(INTG), INTENT(IN) :: rowNodeNumber !<The nodal number to calculate
    TYPE(FieldVariableType), POINTER :: rowsFieldVariable !<A pointer to the field variable associated with the rows
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx,derivativeIdx,interpolationType,localRow,numberOfComponents,numberOfDerivatives, &
      & numberOfVersions,rowsVariableType,totalNumberOfNodes,versionIdx
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("EquationsMatrices_NodalVectorCalculate",err,error,*999)
    
    IF(.NOT.ASSOCIATED(rowsFieldVariable)) CALL FlagError("Rows field variable is not associated.",err,error,*999)
    CALL FieldVariable_VariableTypeGet(rowsFieldVariable,rowsVariableType,err,error,*999)
    
    !Calculate the rows for the nodal vector
    nodalVector%numberOfRows=0
    IF(updateVector) THEN
      CALL FieldVariable_NumberOfComponentsGet(rowsFieldVariable,numberOfComponents,err,error,*999)
      DO componentIdx=1,numberOfComponents
        NULLIFY(domain)
        CALL FieldVariable_ComponentDomainGet(rowsFieldVariable,componentIdx,domain,err,error,*999)
        NULLIFY(domainTopology)
        CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
        NULLIFY(domainNodes)
        CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
        CALL FieldVariable_ComponentInterpolationGet(rowsFieldVariable,componentIdx,interpolationType,err,error,*999)
        CALL DomainNodes_TotalNumberOfNodesGet(domainNodes,totalNumberOfNodes,err,error,*999)
        IF(rowNodeNumber<1.OR.rowNodeNumber>totalNumberOfNodes) THEN
          localError="Node number "//TRIM(NumberToVString(rowNodeNumber,"*",err,error))// &
            & " is invalid for component number "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
            & " of rows field variable type "//TRIM(NumberToVString(rowsVariableType,"*",err,error))// &
            & ". The nodal number must be between 1 and "//TRIM(NumberToVString(totalNumberOfNodes,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        SELECT CASE(interpolationType)
        CASE(FIELD_CONSTANT_INTERPOLATION)
          CALL FieldVariable_ConstantDOFGet(rowsFieldVariable,componentIdx,localRow,err,error,*999)
          nodalVector%numberOfRows=nodalVector%numberOfRows+1
          nodalVector%rowDofs(nodalVector%numberOfRows)=localRow
        CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
!!TODO: Why is the node number being used as an element number here?          
          CALL FieldVariable_LocalElementDOFGet(rowsFieldVariable,rowNodeNumber,componentIdx,localRow,err,error,*999)
          nodalVector%numberOfRows=nodalVector%numberOfRows+1
          nodalVector%rowDofs(nodalVector%numberOfRows)=localRow
        CASE(FIELD_NODE_BASED_INTERPOLATION)
          CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,rowNodeNumber,numberOfDerivatives,err,error,*999)
          DO derivativeIdx=1,numberOfDerivatives
            CALL DomainNodes_DerivativeNumberOfVersionsGet(domainNodes,derivativeIdx,rowNodeNumber,numberOfVersions, &
              & err,error,*999)
            DO versionIdx=1,numberOfVersions
              CALL FieldVariable_LocalNodeDOFGet(rowsFieldVariable,versionIdx,derivativeIdx,rowNodeNumber,componentIdx, &
                & localRow,err,error,*999)
              nodalVector%numberOfRows=nodalVector%numberOfRows+1
              nodalVector%rowDofs(nodalVector%numberOfRows)=localRow
            ENDDO !versionIdx
          ENDDO !derivativeIdx
        CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE DEFAULT
          localError="The interpolation type of "//TRIM(NumberToVString(interpolationType,"*",err,error))// &
            & " is invalid for component number "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
            & " of rows field variable type "//TRIM(NumberToVString(rowsVariableType,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)          
        END SELECT
      ENDDO !componentIdx
      nodalVector%vector=0.0_DP
    ENDIF
    
    EXITS("EquationsMatrices_NodalVectorCalculate")
    RETURN
999 ERRORSEXITS("EquationsMatrices_NodalVectorCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_NodalVectorCalculate

  !
  !================================================================================================================================
  !

  !>Adds the nodal matrices and rhs vector into the equations matrices and rhs vector.
  SUBROUTINE EquationsMatricesVector_NodeAdd(vectorMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the equations matrices
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: columnIdx,matrixIdx,residualIdx,rowIdx,sourceIdx
    REAL(DP) :: sum
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatricesSourcesType), POINTER :: sourceVectors
    TYPE(EquationsMatrixType), POINTER :: dynamicMatrix,linearMatrix

    ENTERS("EquationsMatricesVector_NodeAdd",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector equations matrices is not associated.",err,error,*999)

    NULLIFY(dynamicMatrices)
    CALL EquationsMatricesVector_DynamicMatricesExists(vectorMatrices,dynamicMatrices,err,error,*999)
    IF(ASSOCIATED(dynamicMatrices)) THEN
      !Add the nodal matrices
      DO matrixIdx=1,dynamicMatrices%numberOfDynamicMatrices
        NULLIFY(dynamicMatrix)
        CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,matrixIdx,dynamicMatrix,err,error,*999)
        IF(dynamicMatrix%updateMatrix) THEN
          !Handle lumped matrices
          IF(dynamicMatrix%lumped) THEN
            DO rowIdx=1,dynamicMatrix%nodalMatrix%numberOfRows
              sum=0.0_DP
              DO columnIdx=1,dynamicMatrix%nodalMatrix%numberOfColumns
                sum=sum+dynamicMatrix%nodalMatrix%matrix(rowIdx,columnIdx)
                dynamicMatrix%nodalMatrix%matrix(rowIdx,columnIdx)=0.0_DP
              ENDDO !columnIdx
              dynamicMatrix%nodalMatrix%matrix(rowIdx,rowIdx)=sum
              !Add the nodal matrice into the distributed equations matrix
              CALL DistributedMatrix_ValuesAdd(dynamicMatrix%matrix,dynamicMatrix%nodalMatrix%rowDofs(rowIdx), &
                & dynamicMatrix%nodalMatrix%columnDofs(rowIdx),dynamicMatrix%nodalMatrix%matrix(rowIdx,rowIdx), &
                & err,error,*999)
            ENDDO !rowIdx
          ELSE
            !Add the nodal matrice into the distributed equations matrix
            CALL DistributedMatrix_ValuesAdd(dynamicMatrix%matrix,dynamicMatrix%nodalMatrix%rowDofs(1: &
              & dynamicMatrix%nodalMatrix%numberOfRows),dynamicMatrix%nodalMatrix%columnDofs(1: &
              & dynamicMatrix%nodalMatrix%numberOfColumns),dynamicMatrix%nodalMatrix%matrix(1: &
              & dynamicMatrix%nodalMatrix%numberOfRows,1:dynamicMatrix%nodalMatrix%numberOfColumns), &
              & err,error,*999)
          ENDIF
        ENDIF
      ENDDO !matrixIdx
    ENDIF
    NULLIFY(linearMatrices)
    CALL EquationsMatricesVector_LinearMatricesExists(vectorMatrices,linearMatrices,err,error,*999)
    IF(ASSOCIATED(linearMatrices)) THEN
      !Add the nodal matrices
      DO matrixIdx=1,linearMatrices%numberOfLinearMatrices
        NULLIFY(linearMatrix)
        CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,matrixIdx,linearMatrix,err,error,*999)
        IF(linearMatrix%updateMatrix) THEN
          !Handle lumped matrices
          IF(linearMatrix%lumped) THEN
            DO rowIdx=1,linearMatrix%nodalMatrix%numberOfRows
              sum=0.0_DP
              DO columnIdx=1,linearMatrix%nodalMatrix%numberOfColumns
                sum=sum+linearMatrix%nodalMatrix%matrix(rowIdx,columnIdx)
                linearMatrix%nodalMatrix%matrix(rowIdx,columnIdx)=0.0_DP
              ENDDO !columnIdx
              linearMatrix%nodalMatrix%matrix(rowIdx,rowIdx)=sum
              !Add the nodal matrice into the distributed equations matrix
              CALL DistributedMatrix_ValuesAdd(linearMatrix%matrix,linearMatrix%nodalMatrix%rowDofs(rowIdx), &
                & linearMatrix%nodalMatrix%columnDofs(rowIdx),linearMatrix%nodalMatrix%matrix(rowIdx,rowIdx), &
                & err,error,*999)
            ENDDO !rowIdx
          ELSE
            !Add the nodal matrice into the distributed equations matrix
            CALL DistributedMatrix_ValuesAdd(linearMatrix%matrix,linearMatrix%nodalMatrix%rowDofs(1: &
              & linearMatrix%nodalMatrix%numberOfRows),linearMatrix%nodalMatrix%columnDofs(1: &
              & linearMatrix%nodalMatrix%numberOfColumns),linearMatrix%nodalMatrix%matrix(1: &
              & linearMatrix%nodalMatrix%numberOfRows,1:linearMatrix%nodalMatrix%numberOfColumns), &
              & err,error,*999)
          ENDIF
        ENDIF
      ENDDO !matrixIdx
    ENDIF
    NULLIFY(nonlinearMatrices)
    CALL EquationsMatricesVector_NonlinearMatricesExists(vectorMatrices,nonlinearMatrices,err,error,*999)
    IF(ASSOCIATED(nonlinearMatrices)) THEN
      DO residualIdx=1,nonlinearMatrices%numberOfResiduals
        NULLIFY(residualVector)
        CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,residualIdx,residualVector,err,error,*999)
        IF(residualVector%updateResidual) THEN
          !Add the residual nodal vector
          CALL DistributedVector_ValuesAdd(residualVector%residual,residualVector%nodalResidual%rowDofs(1: &
            & residualVector%nodalResidual%numberOfRows),residualVector%nodalResidual%vector(1:residualVector% &
            & nodalResidual%numberOfRows),err,error,*999)
        ENDIF
      ENDDO !residualIdx
    ENDIF
    NULLIFY(rhsVector)
    CALL EquationsMatricesVector_RHSVectorExists(vectorMatrices,rhsVector,err,error,*999)
    IF(ASSOCIATED(rhsVector)) THEN
      IF(rhsVector%updateVector) THEN
        !Add the rhs nodal vector
        CALL DistributedVector_ValuesAdd(rhsVector%vector,rhsVector%nodalVector%rowDofs(1: &
          & rhsVector%nodalVector%numberOfRows),rhsVector%nodalVector%vector(1:rhsVector% &
          & NodalVector%numberOfRows),err,error,*999)
      ENDIF
    ENDIF
    NULLIFY(sourceVectors)
    CALL EquationsMatricesVector_SourceVectorsExists(vectorMatrices,sourceVectors,err,error,*999)
    IF(ASSOCIATED(sourceVectors)) THEN
      DO sourceIdx=1,sourceVectors%numberOfSources
        NULLIFY(sourceVector)
        CALL EquationsMatricesSources_SourceVectorGet(sourceVectors,sourceIdx,sourceVector,err,error,*999)
        IF(sourceVector%updateVector) THEN
          !Add the rhs nodal vector
          CALL DistributedVector_ValuesAdd(sourceVector%vector,sourceVector%nodalVector%rowDofs(1: &
            & sourceVector%nodalVector%numberOfRows),sourceVector%nodalVector%vector(1:sourceVector% &
            & NodalVector%numberOfRows),err,error,*999)
        ENDIF
      ENDDO !sourceIdx
    ENDIF
    
    EXITS("EquationsMatricesVector_NodeAdd")
    RETURN
999 ERRORSEXITS("EquationsMatricesVector_NodeAdd",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_NodeAdd

  !
  !================================================================================================================================
  !

  !>Initialise the nodal calculation information for the equations matrices
  SUBROUTINE EquationsMatricesVector_NodalInitialise(vectorMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !The equations matrices to initialise the nodal information for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx,residualIdx,sourceIdx
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatricesSourcesType), POINTER :: sourceVectors
    TYPE(EquationsMatrixType), POINTER :: dynamicMatrix,linearMatrix
    TYPE(FieldVariableType), POINTER :: colsVariable,rowsVariable
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix
    
    ENTERS("EquationsMatricesVector_NodalInitialise",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector equations matrices is not associated.",err,error,*999)

    NULLIFY(vectorMapping)
    CALL EquationsMatricesVector_VectorMappingGet(vectorMatrices,vectorMapping,err,error,*999)
    NULLIFY(lhsMapping)
    CALL EquationsMappingVector_LHSMappingGet(vectorMapping,lhsMapping,err,error,*999)
    NULLIFY(rowsVariable)
    CALL EquationsMappingLHS_LHSVariableGet(lhsMapping,rowsVariable,err,error,*999)
    NULLIFY(dynamicMatrices)
    CALL EquationsMatricesVector_DynamicMatricesExists(vectorMatrices,dynamicMatrices,err,error,*999)
    IF(ASSOCIATED(dynamicMatrices)) THEN
      !Initialise the dynamic nodal matrices
      NULLIFY(dynamicMapping)
      CALL EquationsMappingVector_DynamicMappingGet(vectorMapping,dynamicMapping,err,error,*999)
      NULLIFY(colsVariable)
      CALL EquationsMappingDynamic_DynamicVariableGet(dynamicMapping,colsVariable,err,error,*999)
      DO matrixIdx=1,dynamicMatrices%numberOfDynamicMatrices
        NULLIFY(dynamicMatrix)
        CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,matrixIdx,dynamicMatrix,err,error,*999)
        CALL EquationsMatrices_NodalMatrixSetup(dynamicMatrix%nodalMatrix,rowsVariable,colsVariable,err,error,*999)
      ENDDO !matrixIdx
    ENDIF
    NULLIFY(linearMatrices)
    CALL EquationsMatricesVector_LinearMatricesExists(vectorMatrices,linearMatrices,err,error,*999)
    IF(ASSOCIATED(linearMatrices)) THEN
      !Initialise the linear nodal matrices
      NULLIFY(linearMapping)
      CALL EquationsMappingVector_LinearMappingGet(vectorMapping,linearMapping,err,error,*999)
      DO matrixIdx=1,linearMatrices%numberOfLinearMatrices
        NULLIFY(linearMatrix)
        CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,matrixIdx,linearMatrix,err,error,*999)
        NULLIFY(colsVariable)
        CALL EquationsMatricesLinear_LinearMatrixVariableGet(linearMatrices,matrixIdx,colsVariable,err,error,*999)
        CALL EquationsMatrices_NodalMatrixSetup(linearMatrix%nodalMatrix,rowsVariable,colsVariable,err,error,*999)
      ENDDO !matrixIdx
    ENDIF
    NULLIFY(nonlinearMatrices)
    CALL EquationsMatricesVector_NonlinearMatricesExists(vectorMatrices,nonlinearMatrices,err,error,*999)
    IF(ASSOCIATED(nonlinearMatrices)) THEN
      !Initialise the Jacobian nodal matrices
      NULLIFY(nonlinearMapping)
      CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)
      DO residualIdx=1,nonlinearMapping%numberOfResiduals
        !Initialise the residual nodal arrays
        NULLIFY(residualMapping)
        CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,residualIdx,residualMapping,err,error,*999)
        NULLIFY(residualVector)
        CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,residualIdx,residualVector,err,error,*999)
        CALL EquationsMatrices_NodalVectorSetup(residualVector%nodalResidual,rowsVariable,err,error,*999)
        residualVector%nodalResidualCalculated=0
        DO matrixIdx=1,residualVector%numberOfJacobians
          NULLIFY(jacobianMatrix)
          CALL EquationsMatricesResidual_JacobianMatrixGet(residualVector,matrixIdx,jacobianMatrix,err,error,*999)
          NULLIFY(colsVariable)
          CALL EquationsMappingResidual_JacobianMatrixVariableGet(residualMapping,matrixIdx,colsVariable,err,error,*999)
          CALL EquationsMatrices_NodalMatrixSetup(jacobianMatrix%nodalJacobian,rowsVariable,colsVariable,err,error,*999)
        ENDDO !matrixIdx
      ENDDO !residualIdx
    ENDIF
    NULLIFY(rhsVector)
    CALL EquationsMatricesVector_RHSVectorExists(vectorMatrices,rhsVector,err,error,*999)
    IF(ASSOCIATED(rhsVector)) THEN
      !Initialise the RHS nodal vector
      CALL EquationsMatrices_NodalVectorSetup(rhsVector%nodalVector,rowsVariable,err,error,*999)
    ENDIF
    NULLIFY(sourceVectors)
    CALL EquationsMatricesVector_SourceVectorsExists(vectorMatrices,sourceVectors,err,error,*999)
    IF(ASSOCIATED(sourceVectors)) THEN
      DO sourceIdx=1,sourceVectors%numberOfSources
        NULLIFY(sourceVector)
        CALL EquationsMatricesSources_SourceVectorGet(sourceVectors,sourceIdx,sourceVector,err,error,*999)
        !Initialise the source nodal vector.
        CALL EquationsMatrices_NodalVectorSetup(sourceVector%nodalVector,rowsVariable,err,error,*999)
      ENDDO !sourceIdx
    ENDIF
    
    EXITS("EquationsMatricesVector_NodalInitialise")
    RETURN
999 ERRORSEXITS("EquationsMatricesVector_NodalInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_NodalInitialise

  !
  !================================================================================================================================
  !

  !>Sets up the nodal matrix for the row and column field variables.
  SUBROUTINE EquationsMatrices_NodalMatrixSetup(nodalMatrix,rowsFieldVariable,colsFieldVariable,err,error,*)

    !Argument variables
    TYPE(NodalMatrixType) :: nodalMatrix !<The nodal matrix to setup
    TYPE(FieldVariableType), POINTER :: rowsFieldVariable !<A pointer to the field variable associated with the rows
    TYPE(FieldVariableType), POINTER :: colsFieldVariable !<A pointer to the field variable associated with the columns
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx,dummyErr,maxNodeInterpParameters,numberOfComponents
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("EquationsMatrices_NodalMatrixSetup",err,error,*998)

    IF(.NOT.ASSOCIATED(rowsFieldVariable)) CALL FlagError("Rows field variable is not associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(colsFieldVariable)) CALL FlagError("Columns field variable is not associated.",err,error,*998)
    IF(ALLOCATED(nodalMatrix%rowDofs)) CALL FlagError("Nodal matrix row dofs already allocated.",err,error,*998)
    IF(ALLOCATED(nodalMatrix%columnDofs)) CALL FlagError("Nodal matrix column dofs already allocated.",err,error,*998)
    IF(ALLOCATED(nodalMatrix%matrix)) CALL FlagError("Nodal matrix already allocated.",err,error,*998)
    
    nodalMatrix%maxNumberOfRows=0
    CALL FieldVariable_NumberOfComponentsGet(rowsFieldVariable,numberOfComponents,err,error,*999)
    DO componentIdx=1,numberOfComponents
      CALL FieldVariable_ComponentMaxNodeInterpParametersGet(rowsFieldVariable,componentIdx,maxNodeInterpParameters, &
        & err,error,*999)
      nodalMatrix%maxNumberOfRows=nodalMatrix%maxNumberOfRows+maxNodeInterpParameters
    ENDDO !componentIdx
    nodalMatrix%maxNumberOfColumns=0
    CALL FieldVariable_NumberOfComponentsGet(colsFieldVariable,numberOfComponents,err,error,*999)
    DO componentIdx=1,numberOfComponents
      CALL FieldVariable_ComponentMaxNodalInterpParametersGet(colsFieldVariable,componentIdx,maxNodeInterpParameters, &
        & err,error,*999)
      nodalMatrix%maxNumberOfColumns=nodalMatrix%maxNumberOfColumns+maxNodeInterpParameters
    ENDDO !componentIdx
    ALLOCATE(nodalMatrix%rowDofs(nodalMatrix%maxNumberOfRows),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate nodal matrix row dofs.",err,error,*999)
    ALLOCATE(nodalMatrix%columnDofs(nodalMatrix%maxNumberOfColumns),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate nodal matrix column dofs.",err,error,*999)
    ALLOCATE(nodalMatrix%matrix(nodalMatrix%maxNumberOfRows,nodalMatrix%maxNumberOfColumns),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate nodal matrix.",err,error,*999)
    
    EXITS("EquationsMatrices_NodalMatrixSetup")
    RETURN
999 CALL EquationsMatrices_NodalMatrixFinalise(nodalMatrix,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMatrices_NodalMatrixSetup",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_NodalMatrixSetup

  !
  !================================================================================================================================
  !

  !>Sets up the nodal vector for the row field variables.
  SUBROUTINE EquationsMatrices_NodalVectorSetup(nodalVector,rowsFieldVariable,err,error,*)

    !Argument variables
    TYPE(NodalVectorType) :: nodalVector !<The nodal vector to setup
    TYPE(FieldVariableType), POINTER :: rowsFieldVariable !<A pointer to the field variable associated with the rows
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx,dummyErr,maxNodeInterpParameters,numberOfComponents
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("EquationsMatrices_NodalVectorSetup",err,error,*998)

    IF(.NOT.ASSOCIATED(rowsFieldVariable)) CALL FlagError("Rows field variable is not associated.",err,error,*998)
    IF(ALLOCATED(nodalVector%rowDofs)) CALL FlagError("Nodal vector row dofs is already allocated.",err,error,*998)
    IF(ALLOCATED(nodalVector%vector)) CALL FlagError("Nodal vector vector already allocated.",err,error,*998)
    
    nodalVector%maxNumberOfRows = 0
    CALL FieldVariable_NumberOfComponentsGet(rowsFieldVariable,numberOfComponents,err,error,*999)
    DO componentIdx=1,numberOfComponents
      CALL FieldVariable_ComponentMaxNodeInterpParametersGet(rowsFieldVariable,componentIdx,maxNodeInterpParameters, &
        & err,error,*999)
      nodalVector%maxNumberOfRows=nodalVector%maxNumberOfRows+maxNodeInterpParameters
    ENDDO !componentIdx
        
    ALLOCATE(nodalVector%rowDofs(nodalVector%maxNumberOfRows),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate nodal vector row dofs.",err,error,*999)
    ALLOCATE(nodalVector%vector(nodalVector%maxNumberOfRows),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate nodal vector vector.",err,error,*999)
    
    EXITS("EquationsMatrices_NodalVectorSetup")
    RETURN
999 CALL EquationsMatrices_NodalVectorFinalise(nodalVector,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMatrices_NodalVectorSetup",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_NodalVectorSetup

  !
  !================================================================================================================================
  !

  !>Finalise the nodal calculation information and deallocate all memory
  SUBROUTINE EquationsMatricesVector_NodalFinalise(vectorMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<The equations matrices for which to finalise the nodals
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx,residualIdx,sourceIdx
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatricesSourcesType), POINTER :: sourceVectors
    TYPE(EquationsMatrixType), POINTER :: dynamicMatrix,linearMatrix
    
    ENTERS("EquationsMatricesVector_NodalFinalise",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector equations matrices is not associated.",err,error,*999)

    NULLIFY(dynamicMatrices)
    CALL EquationsMatricesVector_DynamicMatricesExists(vectorMatrices,dynamicMatrices,err,error,*999)
    IF(ASSOCIATED(dynamicMatrices)) THEN
      !Finalise the dynamic nodal matrices
      DO matrixIdx=1,dynamicMatrices%numberOfDynamicMatrices
        NULLIFY(dynamicMatrix)
        CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,matrixIdx,dynamicMatrix,err,error,*999)
        CALL EquationsMatrices_NodalMatrixFinalise(dynamicMatrix%nodalMatrix,err,error,*999)
       ENDDO !matrixIdx
    ENDIF
    NULLIFY(linearMatrices)
    CALL EquationsMatricesVector_LinearMatricesExists(vectorMatrices,linearMatrices,err,error,*999)
    IF(ASSOCIATED(linearMatrices)) THEN
      !Finalise the linear nodal matrices
      DO matrixIdx=1,linearMatrices%numberOfLinearMatrices
        NULLIFY(linearMatrix)
        CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,matrixIdx,linearMatrix,err,error,*999)
        CALL EquationsMatrices_NodalMatrixFinalise(linearMatrix%nodalMatrix,err,error,*999)
      ENDDO !matrixIdx
    ENDIF
    NULLIFY(nonlinearMatrices)
    CALL EquationsMatricesVector_NonlinearMatricesExists(vectorMatrices,nonlinearMatrices,err,error,*999)
    IF(ASSOCIATED(nonlinearMatrices)) THEN
      DO residualIdx=1,nonlinearMatrices%numberOfResiduals
        NULLIFY(residualVector)
        CALL EquationsMatricesVector_ResidualVectorGet(nonlinearMatrices,residualIdx,residualVector,err,error,*999)
        CALL EquationsMatrices_NodalVectorFinalise(residualVector%nodalResidual,err,error,*999)
        DO matrixIdx=1,residualVector%numberOfJacobians
          NULLIFY(jacobianMatrix)
          CALL EquationsMatricesResidual_JacobianMatrixGet(residualVector,matrixIdx,jacobianMatrix,err,error,*999)
          CALL EquationsMatrices_NodalMatrixFinalise(jacobianMatrix%nodalJacobian,err,error,*999)
        ENDDO !matrixIdx
      ENDDO !residualIdx
    ENDIF
    NULLIFY(rhsVector)
    CALL EquationsMatricesVector_RHSVectorExists(vectorMatrices,rhsVector,err,error,*999)
    IF(ASSOCIATED(rhsVector)) THEN
      !Finalise the nodal vector
      CALL EquationsMatrices_NodalVectorFinalise(rhsVector%nodalVector,err,error,*999)
    ENDIF
    NULLIFY(sourceVectors)
    CALL EquationsMatricesVector_SourceVectorsExists(vectorMatrices,sourceVectors,err,error,*999)
    IF(ASSOCIATED(sourceVectors)) THEN
      DO sourceIdx=1,sourceVectors%numberOfSources
        NULLIFY(sourceVector)
        CALL EquationsMatricesSources_SourceVectorGet(sourceVectors,sourceIdx,sourceVector,err,error,*999)
        !Finalise the nodal source vector
        CALL EquationsMatrices_NodalVectorFinalise(sourceVector%nodalVector,err,error,*999)
      ENDDO !sourceIdx
    ENDIF
   
    EXITS("EquationsMatricesVector_NodalFinalise")
    RETURN
999 ERRORSEXITS("EquationsMatricesVector_NodalFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_NodalFinalise

  !
  !================================================================================================================================
  !

  !>Initialise the nodal matrix.
  SUBROUTINE EquationsMatrices_NodalMatrixInitialise(nodalMatrix,err,error,*)

    !Argument variables
    TYPE(NodalMatrixType) :: nodalMatrix !The nodal matrix to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsMatrices_NodalMatrixInitialise",err,error,*999)

    nodalMatrix%equationsMatrixNumber=0
    nodalMatrix%numberOfRows=0
    nodalMatrix%numberOfColumns=0
    nodalMatrix%maxNumberOfRows=0
    nodalMatrix%maxNumberOfColumns=0
       
    EXITS("EquationsMatrices_NodalMatrixInitialise")
    RETURN
999 ERRORSEXITS("EquationsMatrices_NodalMatrixInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_NodalMatrixInitialise

  !
  !================================================================================================================================
  !

  !>Finalise an nodal matrix and deallocate all memory
  SUBROUTINE EquationsMatrices_NodalMatrixFinalise(nodalMatrix,err,error,*)

    !Argument variables
    TYPE(NodalMatrixType):: nodalMatrix !<The nodal matrix to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("EquationsMatrices_NodalMatrixFinalise",err,error,*999)

    nodalMatrix%equationsMatrixNumber=0
    nodalMatrix%numberOfRows=0
    nodalMatrix%numberOfColumns=0
    nodalMatrix%maxNumberOfRows=0
    nodalMatrix%maxNumberOfColumns=0
    IF(ALLOCATED(nodalMatrix%rowDofs)) DEALLOCATE(nodalMatrix%rowDofs)    
    IF(ALLOCATED(nodalMatrix%columnDofs)) DEALLOCATE(nodalMatrix%columnDofs)
    IF(ALLOCATED(nodalMatrix%matrix)) DEALLOCATE(nodalMatrix%matrix)
    
    EXITS("EquationsMatrices_NodalMatrixFinalise")
    RETURN
999 ERRORSEXITS("EquationsMatrices_NodalMatrixFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_NodalMatrixFinalise

  !
  !================================================================================================================================
  !

  !>Initialise the nodal vector
  SUBROUTINE EquationsMatrices_NodalVectorInitialise(nodalVector,err,error,*)

    !Argument variables
    TYPE(NodalVectorType) :: nodalVector !The nodal vector to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsMatrices_NodalVectorInitialise",err,error,*999)

    nodalVector%numberOfRows=0
    nodalVector%maxNumberOfRows=0
       
    EXITS("EquationsMatrices_NodalVectorInitialise")
    RETURN
999 ERRORSEXITS("EquationsMatrices_NodalVectorInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_NodalVectorInitialise
  
  !
  !================================================================================================================================
  !

  !>Finalise an nodal vector and deallocate all memory
  SUBROUTINE EquationsMatrices_NodalVectorFinalise(nodalVector,err,error,*)

    !Argument variables
    TYPE(NodalVectorType):: nodalVector !<The nodal vector to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("EquationsMatrices_NodalVectorFinalise",err,error,*999)

    nodalVector%numberOfRows=0
    nodalVector%maxNumberOfRows=0
    IF(ALLOCATED(nodalVector%rowDofs)) DEALLOCATE(nodalVector%rowDofs)
    IF(ALLOCATED(nodalVector%vector)) DEALLOCATE(nodalVector%vector)
    
    EXITS("EquationsMatrices_NodalVectorFinalise")
    RETURN
999 ERRORSEXITS("EquationsMatrices_NodalVectorFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_NodalVectorFinalise

  !
  !================================================================================================================================
  !

  !>Adds the Jacobian matrices into the equations Jacobian.
  SUBROUTINE EquationsMatricesVector_JacobianNodeAdd(vectorMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the equations matrices
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: jacobianMatrixIdx,residualIdx
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix

    ENTERS("EquationsMatricesVector_JacobianNodeAdd",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector equations matrices is not allocated.",err,error,*999)

    NULLIFY(nonlinearMatrices)
    CALL EquationsMatricesVector_NonlinearMatricesExists(vectorMatrices,nonlinearMatrices,err,error,*999)
    IF(ASSOCIATED(nonlinearMatrices)) THEN
      DO residualIdx=1,nonlinearMatrices%numberOfResiduals
        NULLIFY(residualVector)
        CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,residualIdx,residualVector,err,error,*999)
        DO jacobianMatrixIdx=1,residualVector%numberOfJacobians
          NULLIFY(jacobianMatrix)
          CALL EquationsMatricesResidual_JacobianMatrixGet(residualVector,jacobianMatrixIdx,jacobianMatrix,err,error,*999)
          IF(jacobianMatrix%updateJacobian) THEN
            !Add in Jacobian element matrices
            CALL DistributedMatrix_ValuesAdd(jacobianMatrix%jacobian,jacobianMatrix%nodalJacobian%rowDofs(1: &
              & jacobianMatrix%nodalJacobian%numberOfRows),jacobianMatrix%nodalJacobian%columnDofs(1: &
              & jacobianMatrix%nodalJacobian%numberOfColumns),jacobianMatrix%nodalJacobian%matrix(1: &
              & jacobianMatrix%nodalJacobian%numberOfRows,1:jacobianMatrix%nodalJacobian%numberOfColumns), &
              & err,error,*999)
          ENDIF
        ENDDO !jacobianMatrixIdx
      ENDDO !residualIdx
    ENDIF
   
    EXITS("EquationsMatricesVector_JacobianNodeAdd")
    RETURN
999 ERRORSEXITS("EquationsMatricesVector_JacobianNodeAdd",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_JacobianNodeAdd

  !
  !================================================================================================================================
  !

  !>Finalise the element calculation information and deallocate all memory
  SUBROUTINE EquationsMatricesVector_ElementFinalise(vectorMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<The equations matrices for which to finalise the elements
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx,residualIdx,sourceIdx
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatricesSourcesType), POINTER :: sourceVectors
    TYPE(EquationsMatrixType), POINTER :: dynamicMatrix,linearMatrix
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix
    
    ENTERS("EquationsMatricesVector_ElementFinalise",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector equations matrices is not associated.",err,error,*999)

    NULLIFY(dynamicMatrices)
    CALL EquationsMatricesVector_DynamicMatricesExists(vectorMatrices,dynamicMatrices,err,error,*999)
    IF(ASSOCIATED(dynamicMatrices)) THEN
      !Finalise the dynamic element matrices
      DO matrixIdx=1,dynamicMatrices%numberOfDynamicMatrices
        NULLIFY(dynamicMatrix)
        CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,matrixIdx,dynamicMatrix,err,error,*999)
        CALL EquationsMatrices_ElementMatrixFinalise(dynamicMatrix%elementMatrix,err,error,*999)
      ENDDO !matrixIdx
    ENDIF
    NULLIFY(linearMatrices)
    CALL EquationsMatricesVector_LinearMatricesExists(vectorMatrices,linearMatrices,err,error,*999)
    IF(ASSOCIATED(linearMatrices)) THEN
      !Finalise the linear element matrices
      DO matrixIdx=1,linearMatrices%numberOfLinearMatrices
        NULLIFY(linearMatrix)
        CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,matrixIdx,linearMatrix,err,error,*999)
        CALL EquationsMatrices_ElementMatrixFinalise(linearMatrix%elementMatrix,err,error,*999)
      ENDDO !matrixIdx
    ENDIF
    NULLIFY(nonlinearMatrices)
    CALL EquationsMatricesVector_NonlinearMatricesExists(vectorMatrices,nonlinearMatrices,err,error,*999)
    IF(ASSOCIATED(nonlinearMatrices)) THEN
      DO residualIdx=1,nonlinearMatrices%numberOfResiduals
        NULLIFY(residualVector)
        CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,residualIdx,residualVector,err,error,*999)
        CALL EquationsMatrices_ElementVectorFinalise(residualVector%elementResidual,err,error,*999)
        DO matrixIdx=1,residualVector%numberOfJacobians
          NULLIFY(jacobianMatrix)
          CALL EquationsMatricesResidual_JacobianMatrixGet(residualVector,matrixIdx,jacobianMatrix,err,error,*999)
          CALL EquationsMatrices_ElementMatrixFinalise(jacobianMatrix%elementJacobian,err,error,*999)
        ENDDO !matrixIdx
      ENDDO !residualIdx
    ENDIF
    NULLIFY(rhsVector)
    CALL EquationsMatricesVector_RHSVectorExists(vectorMatrices,rhsVector,err,error,*999)
    IF(ASSOCIATED(rhsVector)) THEN
      !Finalise the element vector
      CALL EquationsMatrices_ElementVectorFinalise(rhsVector%elementVector,err,error,*999)
    ENDIF
    NULLIFY(sourceVectors)
    CALL EquationsMatricesVector_SourceVectorsExists(vectorMatrices,sourceVectors,err,error,*999)
    IF(ASSOCIATED(sourceVectors)) THEN
      DO sourceIdx=1,sourceVectors%numberOfSources
        NULLIFY(sourceVector)
        CALL EquationsMatricesSources_SourceVectorGet(sourceVectors,sourceIdx,sourceVector,err,error,*999)
        !Finalise the element source vector
        CALL EquationsMatrices_ElementVectorFinalise(sourceVector%elementVector,err,error,*999)
      ENDDO !sourceIdx
    ENDIF
    
    EXITS("EquationsMatricesVector_ElementFinalise")
    RETURN
999 ERRORSEXITS("EquationsMatricesVector_ElementFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_ElementFinalise

  !
  !================================================================================================================================
  !

  !>Initialise the element calculation information for the equations matrices
  SUBROUTINE EquationsMatricesVector_ElementInitialise(vectorMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !The equations matrices to initialise the element information for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx,residualIdx,sourceIdx
    INTEGER(INTG) :: rowsNumberOfElements,colsNumberOfElements !Number of elements in the row and col variables whose dofs are present in the element matrix
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatricesSourcesType), POINTER :: sourceVectors
    TYPE(EquationsMatrixType), POINTER :: dynamicMatrix,linearMatrix
    TYPE(FieldVariableType), POINTER :: colsVariable,rowsVariable
    
    ENTERS("EquationsMatricesVector_ElementInitialise",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector equations matrices is not associated.",err,error,*999)
    
    rowsNumberOfElements=1
    colsNumberOfElements=1
    NULLIFY(vectorMapping)
    CALL EquationsMatricesVector_VectorMappingGet(vectorMatrices,vectorMapping,err,error,*999)
    NULLIFY(lhsMapping)
    CALL EquationsMappingVector_LHSMappingGet(vectorMapping,lhsMapping,err,error,*999)
    NULLIFY(rowsVariable)
    CALL EquationsMappingLHS_LHSVariableGet(lhsMapping,rowsVariable,err,error,*999)
    NULLIFY(dynamicMatrices)
    CALL EquationsMappingVector_DynamicMatricesExists(vectorMatrices,dynamicMatrices,err,err,*999)
    IF(ASSOCIATED(dynamicMatrices)) THEN
      !Initialise the dynamic element matrices
      NULLIFY(dynamicMapping)
      CALL EquationsMappingVector_DynamicMappingGet(vectorMapping,dynamicMapping,err,error,*999)
      NULLIFY(colsVariable)
      CALL EquationsMappingDynamic_DynamicVariableGet(dynamicMapping,colsVariable,err,error,*999)
      DO matrixIdx=1,dynamicMatrices%numberOfDynamicMatrices
        NULLIFY(dynamicMatrix)
        CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,matrixIdx,dynamicMatrix,err,error,*999)        
        CALL EquationsMatrices_ElementMatrixSetup(dynamicMatrix%elementMatrix,rowsVariable,colsVariable, &
          & rowsNumberOfElements,colsNumberOfElements,err,error,*999)
      ENDDO !matrixIdx
    ENDIF
    NULLIFY(LinearMatrices)
    CALL EquationsMappingVector_LinearMatricesExists(vectorMatrices,linearMatrices,err,err,*999)
    IF(ASSOCIATED(linearMatrices)) THEN
      !Initialise the linear element matrices
      NULLIFY(linearMapping)
      CALL EquationsMappingVector_LinearMappingGet(vectorMapping,linearMapping,err,error,*999)
      DO matrixIdx=1,linearMatrices%numberOfLinearMatrices
        NULLIFY(linearMatrix)
        CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,matrixIdx,linearMatrix,err,error,*999)
        NULLIFY(colsVariable)
        CALL EquationsMappingLinear_LinearMatrixVariableGet(linearMapping,matrixIdx,colsVariable,err,error,*999)
        CALL EquationsMatrices_ElementMatrixSetup(linearMatrix%elementMatrix,rowsVariable,colsVariable, &
          & rowsNumberOfElements,colsNumberOfElements,err,error,*999)
      ENDDO !matrixIdx
    ENDIF
    NULLIFY(nonlinearMatrices)
    CALL EquationsMappingVector_NonlinearMatricesExists(vectorMatrices,nonlinearMatrices,err,err,*999)
    IF(ASSOCIATED(nonlinearMatrices)) THEN
      NULLIFY(nonlinearMapping)
      CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)
      DO residualIdx=1,nonlinearMatrices%numberOfResiduals
        !Initialise the residual element vectors
        NULLIFY(residualMapping)
        CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,residualIdx,residualMapping,err,error,*999)
        NULLIFY(residualVector)
        CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,residualIdx,residualVector,err,error,*999)
        CALL EquationsMatrices_ElementVectorSetup(residualVector%elementResidual,rowsVariable,err,error,*999)
        residualVector%elementResidualCalculated=0
        DO matrixIdx=1,residualVector%numberOfJacobians
          !Initialise the Jacobian element matrices
          NULLIFY(jacobianMatrix)
          CALL EquationsMatricesResidual_JacobianMatrixGet(residualVector,matrixIdx,jacobianMatrix,err,error,*999)
          NULLIFY(colsVariable)
          CALL EquationsMappingResidual_JacobianMatrixVariableGet(residualMapping,matrixIdx,colsVariable,err,error,*999)
          CALL EquationsMatrices_ElementMatrixSetup(jacobianMatrix%elementJacobian,rowsVariable,colsVariable, &
            & rowsNumberOfElements,colsNumberOfElements,err,error,*999)
        ENDDO !matrixIdx
      ENDDO !residualIdx
    ENDIF
    NULLIFY(rhsVector)
    CALL EquationsMatricesVector_RHSVectorExists(vectorMatrices,rhsVector,err,error,*999)
    IF(ASSOCIATED(rhsVector)) THEN
      !Initialise the RHS element vector
      CALL EquationsMatrices_ElementVectorSetup(rhsVector%elementVector,rowsVariable,err,error,*999)
    ENDIF
    NULLIFY(sourceVectors)
    CALL EquationsMatricesVector_SourceVectorsExists(vectorMatrices,sourceVectors,err,error,*999)
    IF(ASSOCIATED(sourceVectors)) THEN
      DO sourceIdx=1,sourceVectors%numberOfSources
        NULLIFY(sourceVector)
        CALL EquationsMatricesSources_SourceVectorGet(sourceVectors,sourceIdx,sourceVector,err,error,*999)
        !Initialise the source element vector. 
        CALL EquationsMatrices_ElementVectorSetup(sourceVector%elementVector,rowsVariable,err,error,*999)
      ENDDO !sourceIdx
    ENDIF
   
    EXITS("EquationsMatricesVector_ElementInitialise")
    RETURN
999 ERRORSEXITS("EquationsMatricesVector_ElementInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_ElementInitialise

  !
  !================================================================================================================================
  !

  !>Finalise a equations matrix and deallocate all memory
  SUBROUTINE EquationsMatrices_EquationsMatrixFinalise(equationsMatrix,err,error,*)

    !Argument variables
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix !<A pointer to the equations matrix to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("EquationsMatrices_EquationsMatrixFinalise",err,error,*999)

    IF(ASSOCIATED(equationsMatrix)) THEN
      IF(ASSOCIATED(equationsMatrix%matrix)) CALL DistributedMatrix_Destroy(equationsMatrix%matrix,err,error,*999)
      CALL EquationsMatrices_ElementMatrixFinalise(equationsMatrix%elementMatrix,err,error,*999)
      CALL EquationsMatrices_NodalMatrixFinalise(equationsMatrix%nodalMatrix,err,error,*999)
      IF(ASSOCIATED(equationsMatrix%tempVector)) CALL DistributedVector_Destroy(equationsMatrix%tempVector,err,error,*999)
      DEALLOCATE(equationsMatrix)
    ENDIF
    
    EXITS("EquationsMatrices_EquationsMatrixFinalise")
    RETURN
999 ERRORSEXITS("EquationsMatrices_EquationsMatrixFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_EquationsMatrixFinalise

  !
  !================================================================================================================================
  !

  !>Initialise the dynamic equations matrix.
  SUBROUTINE EquationsMatrices_EquationsMatrixDynamicInitialise(dynamicMatrices,matrixNumber,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices !<A pointer to the dynamic matrices to initialise the dynamic equations matrix for
    INTEGER(INTG) :: matrixNumber !<The dynamic matrix number in the dynamic equations matrices to initialise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix
    TYPE(EquationsMatrixToVarMapType), POINTER :: equationsMatrixToVarMap
    TYPE(VARYING_STRING) :: dummyError,localError

    ENTERS("EquationsMatrices_EquationsMatrixDynamicInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(dynamicMatrices)) CALL FlagError("Dynamic matrices is not associated.",err,error,*998)
    IF(matrixNumber<1.OR.matrixNumber>dynamicMatrices%numberOfDynamicMatrices) THEN
      localError="The specified dynamic matrix number of "//TRIM(NumberToVString(matrixNumber,"*",err,error))// &
        & " is invalid. The matrix number must be >= 1 and <= "// &
        & TRIM(NumberToVString(dynamicMatrices%numberOfDynamicMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*998)
    ENDIF
    IF(.NOT.ALLOCATED(dynamicMatrices%matrices)) CALL FlagError("Dynamic matrices matrices is not allocated.",err,error,*998)
    IF(ASSOCIATED(dynamicMatrices%matrices(matrixNumber)%ptr)) THEN
      localError="Equations matrix for dynamic matrix number "//TRIM(NumberToVString(matrixNumber,"*",err,error))// &
        & " is already associated."
      CALL FlagError(localError,err,error,*998)
    ENDIF
    NULLIFY(dynamicMapping)
    CALL EquationsMatricesDynamic_DynamicMappingGet(dynamicMatrices,dynamicMapping,err,error,*999)
    NULLIFY(equationsMatrixToVarMap)
    CALL EquationsMappingDynamic_EquationsMatrixToVarMapGet(dynamicMapping,matrixNumber,equationsMatrixToVarMap,err,error,*999)
    
    ALLOCATE(dynamicMatrices%matrices(matrixNumber)%ptr,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate equations matrix.",err,error,*999)
    equationsMatrix=>dynamicMatrices%matrices(matrixNumber)%ptr
    equationsMatrix%matrixNumber=matrixNumber
    equationsMatrix%matrixMatricesType=EQUATIONS_MATRIX_DYNAMIC
    equationsMatrix%dynamicMatrices=>dynamicMatrices
    NULLIFY(equationsMatrix%linearMatrices)
    equationsMatrix%matrixCoefficient=equationsMatrixToVarMap%matrixCoefficient
    equationsMatrix%storageType=MATRIX_BLOCK_STORAGE_TYPE
    equationsMatrix%structureType=EQUATIONS_MATRIX_NO_STRUCTURE
    equationsMatrix%lumped=.FALSE.
    equationsMatrix%symmetric=.FALSE.
    equationsMatrix%updateMatrix=.TRUE.
    equationsMatrix%firstAssembly=.TRUE.
    equationsMatrix%numberOfColumns=equationsMatrixToVarMap%numberOfColumns
    equationsMatrixToVarMap%equationsMatrix=>equationsMatrix
    NULLIFY(equationsMatrix%matrix)
    CALL EquationsMatrices_ElementMatrixInitialise(equationsMatrix%elementMatrix,err,error,*999)
    CALL EquationsMatrices_NodalMatrixInitialise(equationsMatrix%nodalMatrix,err,error,*999)
    NULLIFY(equationsMatrix%tempVector)
    
    EXITS("EquationsMatrices_EquationsMatrixDynamicInitialise")
    RETURN
999 CALL EquationsMatrices_EquationsMatrixFinalise(dynamicMatrices%matrices(matrixNumber)%ptr,dummyErr,dummyError,*998)
998 ERRORS("EquationsMatrices_EquationsMatrixDynamicInitialise",err,error)
    EXITS("EquationsMatrices_EquationsMatrixDynamicInitialise")
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_EquationsMatrixDynamicInitialise

  !
  !================================================================================================================================
  !

  !>Initialise the linear equations matrix.
  SUBROUTINE EquationsMatrices_EquationsMatrixLinearInitialise(linearMatrices,matrixNumber,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices !<A pointer to the linear matrices to initialise the linear equations matrix for
    INTEGER(INTG) :: matrixNumber !<The linear matrix number in the linear equations matrices to initialise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix
    TYPE(EquationsMatrixToVarMapType), POINTER :: equationsMatrixToVarMap
    TYPE(VARYING_STRING) :: dummyError,localError

    ENTERS("EquationsMatrices_EquationsMatrixLinearInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(linearMatrices)) CALL FlagError("Linear matrices is not associated.",err,error,*998)    
    IF(matrixNumber<1.OR.matrixNumber>linearMatrices%numberOfLinearMatrices) THEN
      localError="The specified linear matrix number of "//TRIM(NumberToVString(matrixNumber,"*",err,error))// &
        & " is invalid. The matrix number must be >= 1 and <= "// &
        & TRIM(NumberToVString(linearMatrices%numberOfLinearMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*998)
    ENDIF
    IF(.NOT.ALLOCATED(linearMatrices%matrices)) CALL FlagError("Linear matrices matrices is not allocated.",err,error,*998)
    IF(ASSOCIATED(linearMatrices%matrices(matrixNumber)%ptr)) THEN
      localError="Equations matrix for linear matrix number "//TRIM(NumberToVString(matrixNumber,"*",err,error))// &
        & " is already associated."
      CALL FlagError(localError,err,error,*998)
    ENDIF
    NULLIFY(linearMapping)
    CALL EquationsMatricesLinear_LinearMappingGet(linearMatrices,linearMapping,err,error,*999)
    NULLIFY(equationsMatrixToVarMap)
    CALL EquationsMappingLinear_EquationsMatrixToVarGet(linearMapping,matrixNumber,equationsMatrixToVarMap,err,error,*999)
    
    ALLOCATE(linearMatrices%matrices(matrixNumber)%ptr,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate equations matrix.",err,error,*999)
    equationsMatrix=>linearMatrices%matrices(matrixNumber)%ptr
    equationsMatrix%matrixNumber=matrixNumber
    equationsMatrix%matrixMatricesType=EQUATIONS_MATRIX_LINEAR
    NULLIFY(equationsMatrix%dynamicMatrices)
    equationsMatrix%linearMatrices=>linearMatrices
    equationsMatrix%matrixCoefficient=1.0_DP
    equationsMatrix%storageType=MATRIX_BLOCK_STORAGE_TYPE
    equationsMatrix%structureType=EQUATIONS_MATRIX_NO_STRUCTURE
    equationsMatrix%lumped=.FALSE.
    equationsMatrix%symmetric=.FALSE.
    equationsMatrix%updateMatrix=.TRUE.
    equationsMatrix%firstAssembly=.TRUE.
    equationsMatrix%numberOfColumns=equationsMatrixToVarMap%numberOfColumns
    equationsMatrixToVarMap%equationsMatrix=>equationsMatrix
    NULLIFY(equationsMatrix%matrix)
    CALL EquationsMatrices_ElementMatrixInitialise(equationsMatrix%elementMatrix,err,error,*999)
    CALL EquationsMatrices_NodalMatrixInitialise(equationsMatrix%nodalMatrix,err,error,*999)
    NULLIFY(equationsMatrix%tempVector)
   
    EXITS("EquationsMatrices_EquationsMatrixLinearInitialise")
    RETURN
999 CALL EquationsMatrices_EquationsMatrixFinalise(linearMatrices%matrices(matrixNumber)%ptr,dummyErr,dummyError,*998)
998 ERRORS("EquationsMatrices_EquationsMatrixLinearInitialise",err,error)
    EXITS("EquationsMatrices_EquationsMatrixLinearInitialise")
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_EquationsMatrixLinearInitialise

  !
  !================================================================================================================================
  !

  !>Finalises the vector equations matrices dynamic matrices and deallocates all memory
  SUBROUTINE EquationsMatricesVector_DynamicFinalise(dynamicMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices !<A pointer to the vector equation matrices dynamic matrices to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
     
    ENTERS("EquationsMatricesVector_DynamicFinalise",err,error,*999)

    IF(ASSOCIATED(dynamicMatrices)) THEN
      IF(ALLOCATED(dynamicMatrices%matrices)) THEN
        DO matrixIdx=1,SIZE(dynamicMatrices%matrices,1)
          CALL EquationsMatrices_EquationsMatrixFinalise(dynamicMatrices%matrices(matrixIdx)%ptr,err,error,*999)
        ENDDO !matrixIdx
        DEALLOCATE(dynamicMatrices%matrices)
      ENDIF
      IF(ASSOCIATED(dynamicMatrices%tempVector)) CALL DistributedVector_Destroy(dynamicMatrices%tempVector,err,error,*999)
      DEALLOCATE(dynamicMatrices)
    ENDIF
    
    EXITS("EquationsMatricesVector_DynamicFinalise")
    RETURN
999 ERRORSEXITS("EquationsMatricesVector_DynamicFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_DynamicFinalise
  
  !
  !================================================================================================================================
  !

  !>Initialises the vector equations matrices dynamic matrices
  SUBROUTINE EquationsMatricesVector_DynamicInitialise(vectorMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the vector equation matrices to initialise the dynamic matrices for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,matrixIdx
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMatrixToVarMapType), POINTER :: equationsMatrixToVarMap
    TYPE(VARYING_STRING) :: dummyError
     
    ENTERS("EquationsMatricesVector_DynamicInitialise",err,error,*998)
    
    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector equations matrices is not associated.",err,error,*998)
    IF(ASSOCIATED(vectorMatrices%dynamicMatrices)) &
      & CALL FlagError("Equations matrices dynamic matrices is already associated.",err,error,*998)
    NULLIFY(vectorMapping)
    CALL EquationsMatricesVector_VectorMappingGet(vectorMatrices,vectorMapping,err,error,*998)
    NULLIFY(dynamicMapping)
    CALL EquationsMapping_DynamicMappingExists(vectorMapping,dynamicMapping,err,error,*998)
 
    IF(ASSOCIATED(dynamicMapping)) THEN
      ALLOCATE(vectorMatrices%dynamicMatrices,STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate equations matrices dynamic matrices.",err,error,*999)
      vectorMatrices%dynamicMatrices%vectorMatrices=>vectorMatrices
      vectorMatrices%dynamicMatrices%numberOfDynamicMatrices=dynamicMapping%numberOfDynamicMatrices
      ALLOCATE(vectorMatrices%dynamicMatrices%matrices(dynamicMapping%numberOfDynamicMatrices),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate equations matrices dynamic matrices matrices.",err,error,*999)
      DO matrixIdx=1,dynamicMapping%numberOfDynamicMatrices
        NULLIFY(vectorMatrices%dynamicMatrices%matrices(matrixIdx)%ptr)
        CALL EquationsMatrices_EquationsMatrixDynamicInitialise(vectorMatrices%dynamicMatrices,matrixIdx,err,error,*999)
        NULLIFY(equationsMatrixToVarMap)
        CALL EquationsMappingDynamic_EquationsMatrixToVarMapGet(dynamicMapping,matrixIdx,equationsMatrixToVarMap,err,error,*999)
        vectorMatrices%dynamicMatrices%matrices(matrixIdx)%ptr%matrixCoefficient=equationsMatrixToVarMap%matrixCoefficient
      ENDDO !matrixIdx
      NULLIFY(vectorMatrices%dynamicMatrices%tempVector)
    ENDIF
    
    EXITS("EquationsMatricesVector_DynamicInitialise")
    RETURN
999 CALL EquationsMatricesVector_DynamicFinalise(vectorMatrices%dynamicMatrices,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMatricesVector_DynamicInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_DynamicInitialise
  
  !
  !================================================================================================================================
  !

  !>Adds the Hessian elmental matrices into the equations Hessian.
  SUBROUTINE EquationsMatricesVector_HessianElementAdd(vectorMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the equations matrices
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: hessianMatrixIdx
    TYPE(HessianMatrixType), POINTER :: hessianMatrix
    TYPE(EquationsMatricesOptimisationType), POINTER :: optimisationMatrices
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsMatricesVector_HessianElementAdd",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Equations matrices is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(vectorMatrices%optimisationMatrices)) &
      & CALL FlagError("Equations matrices optimisation matrices is not associated.",err,error,*999)
    NULLIFY(optimisationMatrices)
    CALL EquationsMatricesVector_OptimisationMatricesGet(vectorMatrices,optimisationMatrices,err,error,*999)
    DO hessianMatrixIdx=1,optimisationMatrices%numberOfHessians
      NULLIFY(hessianMatrix)
      CALL EquationsMatricesOptimisation_HessianMatrixGet(optimisationMatrices,hessianMatrixIdx,hessianMatrix,err,error,*999)
      IF(hessianMatrix%updateHessian) THEN
        !Add in Hessian element matrices
        CALL DistributedMatrix_ValuesAdd(hessianMatrix%hessian,hessianMatrix%elementHessian%rowDOFS( &
          & 1:hessianMatrix%elementHessian%numberOfRows),hessianMatrix%elementHessian%columnDOFS( &
          & 1:hessianMatrix%elementHessian%numberOfColumns),hessianMatrix%elementHessian%matrix(1: &
          & hessianMatrix%elementHessian%numberOfRows,1:hessianMatrix%elementHessian%numberOfColumns), &
          & err,error,*999)
      ENDIF
    ENDDO !hessianMatrixIdx
    
    EXITS("EquationsMatricesVector_HessianElementAdd")
    RETURN
999 ERRORSEXITS("EquationsMatricesVector_HessianElementAdd",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_HessianElementAdd

  !
  !================================================================================================================================
  !

  !>Outputs the equations Hessian matrices
  SUBROUTINE EquationsMatricesVector_HessianOutput(id,vectorMatrices,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the ouptut stream
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the equations Hessian matrices
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: hessianMatrixIdx
    TYPE(HessianMatrixType), POINTER :: hessianMatrix
    TYPE(EquationsMatricesOptimisationType), POINTER :: optimisationMatrices
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("EquationsMatricesVector_HessianOutput",err,error,*999)

    CALL EquationsMatricesVector_AssertIsFinished(vectorMatrices,err,error,*999)
    NULLIFY(optimisationMatrices)
    CALL EquationsMatricesVector_OptimisationMatricesGet(vectorMatrices,optimisationMatrices,err,error,*999)
    
    CALL WriteString(id,"",err,error,*999)
    CALL WriteString(id,"Hessian matrices:",err,error,*999)
    DO hessianMatrixIdx=1,optimisationMatrices%numberOfHessians
      NULLIFY(hessianMatrix)
      CALL EquationsMatricesOptimisation_HessianMatrixGet(optimisationMatrices,hessianMatrixIdx,hessianMatrix,err,error,*999)
      CALL WriteStringValue(id,"Hessian matrix: ",hessianMatrixIdx,err,error,*999)
      CALL DistributedMatrix_Output(id,hessianMatrix%hessian,err,error,*999)
    ENDDO !hessianMatrixIdx
    
    EXITS("EquationsMatricesVector_HessianOutput")
    RETURN
999 ERRORSEXITS("EquationsMatricesVector_HessianOutput",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_HessianOutput
  
  !
  !================================================================================================================================
  !

  !>Adds the Jacobain matrices into the equations Jacobian.
  SUBROUTINE EquationsMatricesVector_JacobianElementAdd(vectorMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the equations matrices
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: jacobianMatrixIdx,residualIdx
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix

    ENTERS("EquationsMatricesVector_JacobianElementAdd",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector equations matrices is not associated.",err,error,*999)
    NULLIFY(nonlinearMatrices)
    CALL EquationsMatricesVector_NonlinearMatricesGet(vectorMatrices,nonlinearMatrices,err,error,*999)
    DO residualIdx=1,nonlinearMatrices%numberOfResiduals
      NULLIFY(residualVector)
      CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,residualIdx,residualVector,err,error,*999)
      DO jacobianMatrixIdx=1,residualVector%numberOfJacobians
        NULLIFY(jacobianMatrix)
        CALL EquationsMatricesResidual_JacobianMatrixGet(residualVector,jacobianMatrixIdx,jacobianMatrix,err,error,*999)
        IF(jacobianMatrix%updateJacobian) THEN
          !Add in Jacobian element matrices
          CALL DistributedMatrix_ValuesAdd(jacobianMatrix%jacobian,jacobianMatrix%elementJacobian%rowDOFS(1: &
            & jacobianMatrix%elementJacobian%numberOfRows),jacobianMatrix%elementJacobian%columnDOFS(1: &
            & jacobianMatrix%elementJacobian%numberOfColumns),jacobianMatrix%elementJacobian%matrix(1: &
            & jacobianMatrix%elementJacobian%numberOfRows,1:jacobianMatrix%elementJacobian%numberOfColumns), &
            & err,error,*999)
        ENDIF
      ENDDO !jacobianMatrixIdx
    ENDDO !residualIdx
      
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_STOP("EquationsMatricesVector_JacobianElementAdd()")
#endif
    
    EXITS("EquationsMatricesVector_JacobianElementAdd")
    RETURN
999 ERRORSEXITS("EquationsMatricesVector_JacobianElementAdd",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_JacobianElementAdd

  !
  !================================================================================================================================
  !

  !>Outputs the equations Jacobian matrices
  SUBROUTINE EquationsMatricesVector_JacobianOutput(id,vectorMatrices,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the ouptut stream
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the equations Jacobian matrices
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: jacobianMatrixIdx,residualIdx
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix
    
    ENTERS("EquationsMatricesVector_JacobianOutput",err,error,*999)    

    CALL EquationsMatricesVector_AssertIsFinished(vectorMatrices,err,error,*999)
    NULLIFY(nonlinearMatrices)
    CALL EquationsMatricesVector_NonlinearMatricesExists(vectorMatrices,nonlinearMatrices,err,error,*999)
    IF(ASSOCIATED(nonlinearMatrices)) THEN
      CALL WriteString(id,"",err,error,*999)
      CALL WriteString(id,"Jacobian matrices:",err,error,*999)
      CALL WriteStringValue(id,"Number of residual vectors = ",nonlinearMatrices%numberOfResiduals,err,error,*999)
      DO residualIdx=1,nonlinearMatrices%numberOfResiduals
        CALL WriteStringValue(id,"Residual : ",residualIdx,err,error,*999)
        NULLIFY(residualVector)
        CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,residualIdx,residualVector,err,error,*999)
        CALL WriteStringValue(id,"Number of Jacobian matrices = ",residualVector%numberOfJacobians,err,error,*999)
        DO jacobianMatrixIdx=1,residualVector%numberOfJacobians
          NULLIFY(jacobianMatrix)
          CALL EquationsMatricesResidual_JacobianMatrixGet(residualVector,jacobianMatrixIdx,jacobianMatrix,err,error,*999)
          CALL WriteStringValue(id,"Jacobian matrix : ",jacobianMatrixIdx,err,error,*999)
          CALL DistributedMatrix_Output(id,jacobianMatrix%jacobian,err,error,*999)
        ENDDO !jacobianMatrixIdx
      ENDDO !residualIdx
    ENDIF
    
    EXITS("EquationsMatricesVector_JacobianOutput")
    RETURN
999 ERRORSEXITS("EquationsMatricesVector_JacobianOutput",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_JacobianOutput
  
  !
  !================================================================================================================================
  !

  !>Sets the Jacobian calculation types of the residual variables
  SUBROUTINE EquationsMatricesVector_JacobianCalculationTypeSet(vectorMatrices,variableType,residualIndex,jacobianCalculationType, &
    & err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the vector equations matrices to set the Jacobian type for.
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type that the residual is differentiated with respect to for this Jacobian
    INTEGER(INTG), INTENT(IN) :: residualIndex !<The index of the residual vector to get the Jacobian matrix for
    INTEGER(INTG), INTENT(IN) :: jacobianCalculationType !<The Jacobian calculation type for the matrixIdx'th Jacobian matrix.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx,variableIdx
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsMatricesVector_JacobianCalculationTypeSet",err,error,*999)

    CALL EquationsMatricesVector_AssertIsFinished(vectorMatrices,err,error,*999)
    NULLIFY(nonlinearMatrices)
    CALL EquationsMatricesVector_NonlinearMatricesGet(vectorMatrices,nonlinearMatrices,err,error,*999)
    NULLIFY(residualVector)
    CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,residualIndex,err,error,*999)
    NULLIFY(vectorMapping)
    CALL EquationsMatricesVector_VectorMappingGet(vectorMatrices,vectorMapping,err,error,*999)
    NULLIFY(nonlinearMapping)
    CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)
    NULLIFY(residualMapping)
    CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,residualIndex,residualMapping,err,error,*999)
    !Find Jacobian matrix index for the variable type
    CALL EquationsMappingResidual_VariableIndexGet(residualMapping,variableType,matrixIdx,err,error,*999)
    IF(matrixIdx==0) THEN
      localError="Equations do not have a Jacobian matrix for residual index "// &
        & TRIM(NumberToVstring(residualIndex,"*",err,error))//" and variable type "// &
        & TRIM(NumberToVstring(variableType,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    !Now get Jacobian matrix using the matrix index
    NULLIFY(jacobianMatrix)
    CALL EquationsMatricesResidual_JacobianMatrixGet(residualVector,matrixIdx,jacobianMatrix,err,error,*999)
    
    SELECT CASE(jacobianCalculationType)
    CASE(EQUATIONS_JACOBIAN_FINITE_DIFFERENCE_CALCULATED)
      jacobianMatrix%jacobianCalculationType=EQUATIONS_JACOBIAN_FINITE_DIFFERENCE_CALCULATED
    CASE(EQUATIONS_JACOBIAN_ANALYTIC_CALCULATED)
      jacobianMatrix%jacobianCalculationType=EQUATIONS_JACOBIAN_ANALYTIC_CALCULATED
    CASE DEFAULT
      localError="The specified Jacobian calculation type of "//TRIM(NumberToVString(jacobianCalculationType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("EquationsMatricesVector_JacobianCalculationTypeSet")
    RETURN
999 ERRORS("EquationsMatricesVector_JacobianCalculationTypeSet",err,error)
    EXITS("EquationsMatricesVector_JacobianCalculationTypeSet")
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_JacobianCalculationTypeSet

  !
  !================================================================================================================================
  !

  !>Sets/changes the finite difference step size used for calculating the Jacobian
  SUBROUTINE EquationsMatricesVector_JacobianFiniteDifferenceStepSizeSet(vectorMatrices,variableType,residualIndex, &
    & jacobianStepSize,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the vector equations matrices to set the Jacobian step size for.
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type that the residual is differentiated with respect to for this Jacobian
    INTEGER(INTG), INTENT(IN) :: residualIndex !<The index of the residual vector to set the Jacobian step size for
    REAL(DP), INTENT(IN) :: jacobianStepSize !<The finite difference step size for the Jacobian matrix. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx,variableIdx
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsMatrices_JacobianFiniteDifferenceStepSizesSet",err,error,*999)

    CALL EquationsMatricesVector_AssertIsFinished(vectorMatrices,err,error,*999)
    NULLIFY(nonlinearMatrices)
    CALL EquationsMatricesVector_NonlinearMatricesGet(vectorMatrices,nonlinearMatrices,err,error,*999)
    NULLIFY(residualVector)
    CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,residualIndex,err,error,*999)
    NULLIFY(vectorMapping)
    CALL EquationsMatricesVector_VectorMappingGet(vectorMatrices,vectorMapping,err,error,*999)
    NULLIFY(nonlinearMapping)
    CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)
    NULLIFY(residualMapping)
    CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,residualIndex,residualMapping,err,error,*999)
    !Find Jacobian matrix index for the variable type
    CALL EquationsMappingResidual_VariableIndexGet(residualMapping,variableType,matrixIdx,err,error,*999)
    IF(matrixIdx==0) THEN
      localError="Equations do not have a Jacobian matrix for residual index "// &
        & TRIM(NumberToVstring(residualIndex,"*",err,error))//" and variable type "// &
        & TRIM(NumberToVstring(variableType,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    !Now get Jacobian matrix using the matrix index
    NULLIFY(jacobianMatrix)
    CALL EquationsMatricesResidual_JacobianMatrixGet(residualVector,matrixIdx,jacobianMatrix,err,error,*999)
    
    IF(jacobianStepSize<=ZERO_TOLERANCE) THEN
      localError="The specified Jacobian step size of "//TRIM(NumberToVString(jacobianStepSize,"*",err,error))// &
        & " is invalid. The step size must be >= "//TRIM(NumberToVString(ZERO_TOLERANCE,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    jacobianMatrix%jacobianFiniteDifferenceStepSize=jacobianStepSize

    EXITS("EquationsMatricesVector_JacobianFiniteDifferenceStepSizeSet")
    RETURN
999 ERRORS("EquationsMatricesVector_JacobianFiniteDifferenceStepSizeSet",err,error)
    EXITS("EquationsMatricesVector_JacobianFiniteDifferenceStepSizeSet")
    RETURN 1

  END SUBROUTINE EquationsMatricesVector_JacobianFiniteDifferenceStepSizeSet

  !
  !================================================================================================================================
  !

  !>Finalises the vector equations matrices linear matrices and deallocates all memory
  SUBROUTINE EquationsMatricesVector_LinearFinalise(linearMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices !<A pointer to the vector equation matrices linear matrices to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
     
    ENTERS("EquationsMatricesVector_LinearFinalise",err,error,*999)

    IF(ASSOCIATED(linearMatrices)) THEN
      IF(ALLOCATED(linearMatrices%matrices)) THEN
        DO matrixIdx=1,SIZE(linearMatrices%matrices,1)
          CALL EquationsMatrices_EquationsMatrixFinalise(linearMatrices%matrices(matrixIdx)%ptr,err,error,*999)
        ENDDO !matrixIdx
        DEALLOCATE(linearMatrices%matrices)
      ENDIF
      IF(ASSOCIATED(linearMatrices%tempVector)) CALL DistributedVector_Destroy(linearMatrices%tempVector,err,error,*999)
      DEALLOCATE(linearMatrices)
    ENDIF
    
    EXITS("EquationsMatricesVector_LinearFinalise")
    RETURN
999 ERRORSEXITS("EquationsMatricesVector_LinearFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_LinearFinalise
  
  !
  !================================================================================================================================
  !

  !>Initialises the equations matrices linear matrices
  SUBROUTINE EquationsMatricesVector_LinearInitialise(vectorMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the equation matrices to initialise the linear matrices for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,matrixIdx
    REAL(DP) :: matrixCoefficient
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMatrixToVarMapType), POINTER :: equationsMatrixToVarMap
    TYPE(VARYING_STRING) :: dummyError
     
    ENTERS("EquationsMatricesVector_LinearInitialise",err,error,*998)
    
    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector equations matrices is not associated.",err,error,*998)
    IF(ASSOCIATED(vectorMatrices%linearMatrices)) &
      & CALL FlagError("Vector equations matrices linear matrices is already associated.",err,error,*998)
    NULLIFY(vectorMapping)
    CALL EquationsMatricesVector_VectorMappingGet(vectorMatrices,vectorMapping,err,error,*999)
    NULLIFY(linearMapping)
    CALL EquationsMappingVector_LinearMappingExists(vectorMapping,linearMapping,err,error,*999)

    IF(ASSOCIATED(linearMapping)) THEN
      ALLOCATE(vectorMatrices%linearMatrices,STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate equations matrices linear matrices.",err,error,*999)
      vectorMatrices%linearMatrices%vectorMatrices=>vectorMatrices
      vectorMatrices%linearMatrices%numberOfLinearMatrices=linearMapping%numberOfLinearMatrices
      ALLOCATE(vectorMatrices%linearMatrices%matrices(linearMapping%numberOfLinearMatrices),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate equations matrices linear matrices matrices.",err,error,*999)
      DO matrixIdx=1,linearMapping%numberOfLinearMatrices
        NULLIFY(vectorMatrices%linearMatrices%matrices(matrixIdx)%ptr)
        CALL EquationsMatrices_EquationsMatrixLinearInitialise(vectorMatrices%linearMatrices,matrixIdx,err,error,*999)
        NULLIFY(equationsMatrixToVarMap)
        CALL EquationsMappingLinear_EquationsMatrixToVarMapGet(linearMapping,matrixIdx,equationsMatrixToVarMap,err,error,*999)
        CALL EquationsMappingVectorEMToVMap_MatrixCoefficientGet(equationsMatrixToVarMap,matrixCoefficient,err,error,*999)
        vectorMatrices%dynamicMatrices%matrices(matrixIdx)%ptr%matrixCoefficient=matrixCoefficient
      ENDDO !matrixIdx
      NULLIFY(vectorMatrices%linearMatrices%tempVector)
    ENDIF
    
    EXITS("EquationsMatricesVector_LinearInitialise")
    RETURN
999 CALL EquationsMatricesVector_LinearFinalise(vectorMatrices%linearMatrices,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMatricesVector_LinearInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_LinearInitialise
  
  !
  !================================================================================================================================
  !

  !>Finalises the equations matrices nonlinear matrices and deallocates all memory
  SUBROUTINE EquationsMatricesVector_NonlinearFinalise(nonlinearMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices !<A pointer to the equation matrices nonlinear matrices to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: residualIdx
     
    ENTERS("EquationsMatricesVector_NonlinearFinalise",err,error,*999)

    IF(ASSOCIATED(nonlinearMatrices)) THEN
      IF(ALLOCATED(nonlinearMatrices%residuals)) THEN
        DO residualIdx=1,SIZE(nonlinearMatrices%residuals,1)
          CALL EquationsMatricesNonlinear_ResidualFinalise(nonlinearMatrices%residuals(residualIdx)%ptr,err,error,*999)
        ENDDO !residualIdx
        DEALLOCATE(nonlinearMatrices%residuals)
      ENDIF
      IF(ASSOCIATED(nonlinearMatrices%tempVector)) CALL DistributedVector_Destroy(nonlinearMatrices%tempVector,err,error,*999)
      DEALLOCATE(nonlinearMatrices)
    ENDIF
    
    EXITS("EquationsMatricesVector_NonlinearFinalise")
    RETURN
999 ERRORSEXITS("EquationsMatricesVector_NonlinearFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_NonlinearFinalise
  
  !
  !================================================================================================================================
  !

  !>Initialises the equations matrices nonlinear matrices
  SUBROUTINE EquationsMatricesVector_NonlinearInitialise(vectorMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the equation matrices to initialise the nonlinear matrices for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,matrixIdx,numberOfResiduals,residualIdx
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("EquationsMatricesVector_NonlinearInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector equations matrices is not associated.",err,error,*998)
    IF(ASSOCIATED(vectorMatrices%nonlinearMatrices)) &
      & CALL FlagError("Vector equations matrices nonlinear matrices is already associated.",err,error,*998)

    NULLIFY(vectorMapping)
    CALL EquationsMatricesVector_VectorMappingGet(vectorMatrices,vectorMapping,err,error,*999)
    NULLIFY(nonlinearMapping)
    CALL EquationsMappingVector_NonlinearMappingExists(vectorMapping,nonlinearMapping,err,error,*999)
    IF(ASSOCIATED(nonlinearMapping)) THEN
      CALL EquationsMappingNonlinear_NumberOfResidualsGet(nonlinearMapping,numberOfResiduals,err,error,*999)
      ALLOCATE(vectorMatrices%nonlinearMatrices,STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate equations matrices nonlinear matrices.",err,error,*999)
      vectorMatrices%nonlinearMatrices%vectorMatrices=>vectorMatrices
      ALLOCATE(vectorMatrices%nonlinearMatrices%residuals(numberOfResiduals),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate nonlinear matrices residuals.",err,error,*999)
      vectorMatrices%nonlinearMatrices%numberOfResiduals=numberOfResiduals
      DO residualIdx=1,numberOfResiduals
        NULLIFY(vectorMatrices%nonlinearMatrices%residuals(residualIdx)%ptr)
        CALL EquationsMatricesNonlinear_ResidualInitialise(vectorMatrices%nonlinearMatrices,residualIdx,err,error,*999)
      ENDDO !residualIdx
      NULLIFY(vectorMatrices%nonlinearMatrices%tempVector)
   ENDIF
    
    EXITS("EquationsMatricesVector_NonlinearInitialise")
    RETURN
999 CALL EquationsMatricesVector_NonlinearFinalise(vectorMatrices%nonlinearMatrices,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMatricesVector_NonlinearInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_NonlinearInitialise
  
  !
  !================================================================================================================================
  !

  !>Finalises the equations matrices optimisation matrices and deallocates all memory
  SUBROUTINE EquationsMatrices_OptimisationFinalise(optimisationMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesOptimisationType), POINTER :: optimisationMatrices !<A pointer to the equation matrices optimisation matrices to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
     
    ENTERS("EquationsMatrices_OptimisationFinalise",err,error,*999)

    IF(ASSOCIATED(optimisationMatrices)) THEN
      IF(ALLOCATED(optimisationMatrices%hessians)) THEN
        DO matrixIdx=1,SIZE(optimisationMatrices%hessians,1)
          CALL EquationsMatrices_HessianFinalise(optimisationMatrices%hessians(matrixIdx)%ptr,err,error,*999)
        ENDDO !matrixIdx
        DEALLOCATE(optimisationMatrices%hessians)
      ENDIF
      IF(ASSOCIATED(optimisationMatrices%gradient)) CALL DistributedVector_Destroy(optimisationMatrices%gradient,err,error,*999)
      CALL EquationsMatrices_ElementVectorFinalise(optimisationMatrices%elementGradient,err,error,*999)
      IF(ASSOCIATED(optimisationMatrices%constraints)) CALL DistributedVector_Destroy(optimisationMatrices%constraints, &
        & err,error,*999)
      CALL EquationsMatrices_ElementVectorFinalise(optimisationMatrices%elementConstraints,err,error,*999)
      IF(ASSOCIATED(optimisationMatrices%lowerBounds)) CALL DistributedVector_Destroy(optimisationMatrices%lowerBounds, &
        & err,error,*999)
      CALL EquationsMatrices_ElementVectorFinalise(optimisationMatrices%elementLowerBounds,err,error,*999)
      IF(ASSOCIATED(optimisationMatrices%upperBounds)) CALL DistributedVector_Destroy(optimisationMatrices%upperBounds, &
        & err,error,*999)
      CALL EquationsMatrices_ElementVectorFinalise(optimisationMatrices%elementUpperBounds,err,error,*999)
      IF(ASSOCIATED(optimisationMatrices%residual)) CALL DistributedVector_Destroy(optimisationMatrices%residual,err,error,*999)
      CALL EquationsMatrices_ElementVectorFinalise(optimisationMatrices%elementResidual,err,error,*999)
      DEALLOCATE(optimisationMatrices)
    ENDIF
    
    EXITS("EquationsMatrices_OptimisationFinalise")
    RETURN
999 ERRORSEXITS("EquationsMatrices_OptimisationFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_OptimisationFinalise
  
  !
  !================================================================================================================================
  !

  !>Finalises the vector equations matrices RHS vector and deallocates all memory
  SUBROUTINE EquationsMatricesVector_RHSFinalise(rhsVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector !<A pointer to the vector equation matrices RHS vector to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
     
    ENTERS("EquationsMatricesVector_RHSFinalise",err,error,*999)

    IF(ASSOCIATED(rhsVector)) THEN
      IF(ASSOCIATED(rhsVector%vector)) CALL DistributedVector_Destroy(rhsVector%vector,err,error,*999)
      IF(ASSOCIATED(rhsVector%previousRHSVector)) CALL DistributedVector_Destroy(rhsVector%previousRHSVector,err,error,*999)
      IF(ASSOCIATED(rhsVector%previous2RHSVector)) CALL DistributedVector_Destroy(rhsVector%previous2RHSVector,err,error,*999)
      IF(ASSOCIATED(rhsVector%previous3RHSVector)) CALL DistributedVector_Destroy(rhsVector%previous3RHSVector,err,error,*999)
      CALL EquationsMatrices_ElementVectorFinalise(rhsVector%elementVector,err,error,*999)
      CALL EquationsMatrices_NodalVectorFinalise(rhsVector%nodalVector,err,error,*999)
      DEALLOCATE(rhsVector)
    ENDIF      
     
    EXITS("EquationsMatricesVector_RHSFinalise")
    RETURN
999 ERRORSEXITS("EquationsMatricesVector_RHSFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_RHSFinalise
  
  !
  !================================================================================================================================
  !

  !>Initialises the vector equations matrices RHS vector
  SUBROUTINE EquationsMatricesVector_RHSInitialise(vectorMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the vector equation matrices to initialise the rhs vector for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    REAL(DP) :: rhsCoefficient
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("EquationsMatricesVector_RHSInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector equations matrices is not associated.",err,error,*998)
    IF(ASSOCIATED(vectorMatrices%rhsVector)) &
      & CALL FlagError("Vector equations matrices RHS vector is already associated.",err,error,*998)
    NULLIFY(vectorMapping)
    CALL EquationsMatricesVector_VectorMappingGet(vectorMatrices,vectorMapping,err,error,*998)
    NULLIFY(rhsMapping)
    CALL EquationsMappingVector_RHSMappingExists(vectorMapping,rhsMapping,err,error,*998)
    IF(ASSOCIATED(rhsMapping)) THEN
      CALL EquationsMappingRHS_VectorCoefficientGet(rhsMapping,rhsCoefficient,err,error,*999)
      ALLOCATE(vectorMatrices%rhsVector,STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate equations matrices RHS vector.",err,error,*999)
      vectorMatrices%rhsVector%rhsCoefficient=rhsCoefficient
      vectorMatrices%rhsVector%updateVector=.TRUE.
      vectorMatrices%rhsVector%firstAssembly=.TRUE.
      NULLIFY(vectorMatrices%rhsVector%vector)
      NULLIFY(vectorMatrices%rhsVector%previousRHSVector)
      NULLIFY(vectorMatrices%rhsVector%previous2RHSVector)
      NULLIFY(vectorMatrices%rhsVector%previous3RHSVector)
      CALL EquationsMatrices_ElementVectorInitialise(vectorMatrices%rhsVector%elementVector,err,error,*999)
      CALL EquationsMatrices_NodalVectorInitialise(vectorMatrices%rhsVector%nodalVector,err,error,*999)
    ENDIF
    
    EXITS("EquationsMatricesVector_RHSInitialise")
    RETURN
999 CALL EquationsMatricesVector_RHSFinalise(vectorMatrices%rhsVector,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMatricesVector_RHSInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_RHSInitialise
  
   !
  !================================================================================================================================
  !

  !>Destroy the scalar equations matrices
  SUBROUTINE EquationsMatrices_ScalarDestroy(scalarMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesScalarType), POINTER :: scalarMatrices !<A pointer the scalar equations matrices to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsMatrices_ScalarDestroy",err,error,*999)

    IF(.NOT.ASSOCIATED(scalarMatrices)) CALL FlagError("Scalar equations matrices is not associated",err,error,*999)
    
    CALL EquationsMatrices_ScalarFinalise(scalarMatrices,err,error,*999)
        
    EXITS("EquationsMatrices_ScalarDestroy")
    RETURN
999 ERRORSEXITS("EquationsMatrices_ScalarDestroy",err,error)    
    RETURN 1
   
  END SUBROUTINE EquationsMatrices_ScalarDestroy

  !
  !================================================================================================================================
  !

  !>Finalises the vector equations matrices source vectors and deallocates all memory
  SUBROUTINE EquationsMatricesVector_SourcesFinalise(sourceVectors,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesSourcesType), POINTER :: sourceVectors !<A pointer to the vector equation matrices source vectors to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: sourceIdx
     
    ENTERS("EquationsMatricesVector_SourcesFinalise",err,error,*999)

    IF(ASSOCIATED(sourceVectors)) THEN
      IF(ALLOCATED(sourceVectors%sources)) THEN
        DO sourceIdx=1,SIZE(sourceVectors%sources,1)
          CALL EquationsMatricesSources_SourceVectorFinalise(sourceVectors%sources(sourceIdx)%ptr,err,error,*999)
        ENDDO !sourceIdx
        DEALLOCATE(sourceVectors%sources)
      ENDIF
      IF(ASSOCIATED(sourceVectors%tempVector)) CALL DistributedVector_Destroy(sourceVectors%tempVector,err,error,*999)
      DEALLOCATE(sourceVectors)
    ENDIF      
     
    EXITS("EquationsMatricesVector_SourcesFinalise")
    RETURN
999 ERRORSEXITS("EquationsMatricesVector_SourcesFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_SourcesFinalise
  
  !
  !================================================================================================================================
  !

  !>Initialises the equations matrices source vector
  SUBROUTINE EquationsMatricesVector_SourcesInitialise(vectorMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the equation matrices to initialise the source vector for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,numberOfSources,sourceIdx
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingSourcesType), POINTER :: sourcesMapping
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("EquationsMatricesVector_SourcesInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector equations matrices is not associated.",err,error,*998)
    IF(ASSOCIATED(vectorMatrices%sourceVectors)) &
      & CALL FlagError("Vector equations matrices source vectors is already associated.",err,error,*998)
    NULLIFY(vectorMapping)
    CALL EquationsMatricesVector_VectorMappingGet(vectorMatrices,vectorMapping,err,error,*998)
    NULLIFY(sourcesMapping)
    CALL EquationsMappingVector_SourcesMappingExists(vectorMapping,sourcesMapping,err,error,*998)
    IF(ASSOCIATED(sourcesMapping)) THEN
      CALL EquationsMappingSources_NumberOfSourcesGet(sourcesMapping,numberOfSources,err,error,*999)
      ALLOCATE(vectorMatrices%sourceVectors,STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate equations matrices source vectors.",err,error,*999)
      vectorMatrices%sourceVectors%vectorMatrices=>vectorMatrices
      ALLOCATE(vectorMatrices%sourceVectors%sources(numberOfSources),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate the source vectors sources.",err,error,*999)
      vectorMatrices%sourceVectors%numberOfSources=numberOfSources
      DO sourceIdx=1,numberOfSources
        NULLIFY(vectorMatrices%sourceVectors%sources(sourceIdx)%ptr)
        CALL EquationsMatricesSources_SourceVectorInitialise(vectorMatrices%sourceVectors,sourceIdx,err,error,*999)
      ENDDO !sourceIdx
      NULLIFY(vectorMatrices%sourceVectors%tempVector)
    ENDIF
    
    EXITS("EquationsMatricesVector_SourcesInitialise")
    RETURN
999 CALL EquationsMatricesVector_SourcesFinalise(vectorMatrices%sourceVectors,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMatricesVector_SourcesInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_SourcesInitialise
  
  !
  !================================================================================================================================
  !

  !>Sets the lumping of the dynamic equations matrices
  SUBROUTINE EquationsMatricesVector_DynamicLumpingTypeSet(vectorMatrices,lumpingTypes,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the vector eqautions matrices
    INTEGER(INTG), INTENT(IN) :: lumpingTypes(:) !<lumpingTypes(matrixIdx). The lumping type for the matrixIdx'th dynamic equations matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatrixType), POINTER :: dynamicMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("EquationsMatricesVector_DynamicLumpingTypeSet",err,error,*999)

    CALL EquationsMatricesVector_AssertNotFinished(vectorMatrices,err,error,*999)
    NULLIFY(dynamicMatrices)
    CALL EquationsMatricesVector_DynamicMatricesGet(vectorMatrices,dynamicMatrices,err,error,*999)
    IF(SIZE(lumpingTypes,1)<dynamicMatrices%numberOfDynamicMatrices) THEN
      localError="The size of the lumping type array of "//TRIM(NumberToVString(SIZE(lumpingTypes,1),"*",err,error))// &
        & " is invalid. The size should be >= "//TRIM(NumberToVString(dynamicMatrices%numberOfDynamicMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
      
    DO matrixIdx=1,dynamicMatrices%numberOfDynamicMatrices
      NULLIFY(dynamicMatrix)
      CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,matrixIdx,dynamicMatrix,err,error,*999)
      SELECT CASE(lumpingTypes(matrixIdx))
      CASE(EQUATIONS_MATRIX_UNLUMPED)
        dynamicMatrix%lumped=.FALSE.
      CASE(EQUATIONS_MATRIX_LUMPED)
        dynamicMatrix%lumped=.TRUE.        
      CASE DEFAULT
        localError="The specified lumping type of "//TRIM(NumberToVString(lumpingTypes(matrixIdx),"*",err,error))// &
          & " for dynamic matrix number "//TRIM(NumberToVString(matrixIdx,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDDO !matrixIdx
    
    EXITS("EquationsMatricesVector_DynamicLumpingTypeSet")
    RETURN
999 ERRORSEXITS("EquationsMatricesVector_DynamicLumpingTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_DynamicLumpingTypeSet

  !
  !================================================================================================================================
  !

  !>Sets the storage type (sparsity) of the dynamic equations matrices
  SUBROUTINE EquationsMatricesVector_DynamicStorageTypeSet(vectorMatrices,storageTypes,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the vector equations matrices
    INTEGER(INTG), INTENT(IN) :: storageTypes(:) !<storageTypes(matrixIdx). The storage type for the matrixIdx'th dynamic equations matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatrixType), POINTER :: dynamicMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("EquationsMatricesVector_DynamicStorageTypeSet",err,error,*999)

    CALL EquationsMatricesVector_AssertNotFinished(vectorMatrices,err,error,*999)
    NULLIFY(dynamicMatrices)
    CALL EquationsMatricesVector_DynamicMatricesGet(vectorMatrices,dynamicMatrices,err,error,*999)
    IF(SIZE(storageTypes,1)<dynamicMatrices%numberOfDynamicMatrices) THEN
      localError="The size of the storage type array of "//TRIM(NumberToVString(SIZE(storageTypes,1),"*",err,error))// &
        & " is invalid. The size should be >= "//TRIM(NumberToVString(dynamicMatrices%numberOfDynamicMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    DO matrixIdx=1,dynamicMatrices%numberOfDynamicMatrices
      NULLIFY(dynamicMatrix)
      CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,matrixIdx,dynamicMatrix,err,error,*999)
      SELECT CASE(storageTypes(matrixIdx))
      CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
        dynamicMatrix%storageType=DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE
      CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
        dynamicMatrix%storageType=DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE        
      CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
        dynamicMatrix%storageType=DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE
      CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
        dynamicMatrix%storageType=DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE
      CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
        dynamicMatrix%storageType=DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE
      CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
        dynamicMatrix%storageType=DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE
      CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
        dynamicMatrix%storageType=DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE
      CASE DEFAULT
        localError="The specified storage type of "//TRIM(NumberToVString(storageTypes(matrixIdx),"*",err,error))// &
          & " for dynamic matrix number "//TRIM(NumberToVString(matrixIdx,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDDO !matrixIdx
    
    EXITS("EquationsMatricesVector_DynamicStorageTypeSet")
    RETURN
999 ERRORSEXITS("EquationsMatricesVector_DynamicStorageTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_DynamicStorageTypeSet

  !
  !================================================================================================================================
  !

  !>Sets the storage type (sparsity) of a linear equations matrix
  SUBROUTINE EquationsMatricesVector_LinearStorageTypeSet0(vectorMatrices,storageType,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the vector equations matrices
    INTEGER(INTG), INTENT(IN) :: storageType !<The storage type for the linear equations matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
     
    ENTERS("EquationsMatricesVector_LinearStorageTypeSet0",err,error,*999)

    CALL EquationsMatricesVector_LinearStorageTypeSet1(vectorMatrices,[storageType],err,error,*999)
    
    EXITS("EquationsMatricesVector_LinearStorageTypeSet0")
    RETURN
999 ERRORSEXITS("EquationsMatricesVector_LinearStorageTypeSet0",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_LinearStorageTypeSet0

  !
  !================================================================================================================================
  !

  !>Sets the storage type (sparsity) of the linear equations matrices
  SUBROUTINE EquationsMatricesVector_LinearStorageTypeSet1(vectorMatrices,storageTypes,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the vector equations matrices
    INTEGER(INTG), INTENT(IN) :: storageTypes(:) !<storageTypes(matrixIdx). The storage type for the matrixIdx'th linear equations matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatrixType), POINTER :: linearMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("EquationsMatricesVector_LinearStorageTypeSet1",err,error,*999)

    CALL EquationsMatricesVector_AssertNotFinished(vectorMatrices,err,error,*999)
    NULLIFY(linearMatrices)
    CALL EquationsMatricesVector_LinearMatricesGet(vectorMatrices,linearMatrices,err,error,*999)
    IF(SIZE(storageTypes,1)<linearMatrices%numberOfLinearMatrices) THEN
      localError="The size of the storage type array of "//TRIM(NumberToVString(SIZE(storageTypes,1),"*",err,error))// &
        & " is not equal to the number of linear matrices of "// &
        & TRIM(NumberToVString(linearMatrices%numberOfLinearMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    DO matrixIdx=1,linearMatrices%numberOfLinearMatrices
      NULLIFY(linearMatrix)
      CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,matrixIdx,linearMatrix,err,error,*999)
      SELECT CASE(storageTypes(matrixIdx))
      CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
        linearMatrix%storageType=DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE
      CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
        linearMatrix%storageType=DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE        
      CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
        linearMatrix%storageType=DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE
      CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
        linearMatrix%storageType=DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE
      CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
        linearMatrix%storageType=DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE
      CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
        linearMatrix%storageType=DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE
      CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
        linearMatrix%storageType=DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE
      CASE DEFAULT
        localError="The specified storage type of "//TRIM(NumberToVString(storageTypes(matrixIdx),"*",err,error))// &
          & " for the linear matrix number "//TRIM(NumberToVString(matrixIdx,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDDO !matrixIdx
    
    EXITS("EquationsMatricesVector_LinearStorageTypeSet1")
    RETURN
999 ERRORSEXITS("EquationsMatricesVector_LinearStorageTypeSet1",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_LinearStorageTypeSet1

  !
  !================================================================================================================================
  !

  !>Sets the storage type (sparsity) of all nonlinear (Jacobian) equations matrices
  SUBROUTINE EquationsMatricesVector_NonlinearStorageTypeSet0(vectorMatrices,residualIdx,storageType,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the vector equations matrices
    INTEGER(INTG), INTENT(IN) :: residualIdx !<The residual index to set the Jacobian matrices storage types for
    INTEGER(INTG), INTENT(IN) :: storageType !<storageType. The storage type for all Jacobian equations matrices
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG), ALLOCATABLE :: storageTypes(:)
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector

    ENTERS("EquationsMatricesVector_NonlinearStorageTypeSet0",err,error,*998)

    CALL EquationsMatricesVector_AssertNotFinished(vectorMatrices,err,error,*999)
    NULLIFY(nonlinearMatrices)
    CALL EquationsMatricesVector_NonlinearMatricesGet(vectorMatrices,nonlinearMatrices,err,error,*998)
    NULLIFY(residualVector)
    CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,residualIdx,residualVector,err,error,*998)

    ALLOCATE(storageTypes(residualVector%numberOfJacobians),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate storage types.",err,error,*999)
    storageTypes=storageType
    CALL EquationsMatricesVector_NonlinearStorageTypeSet1(vectorMatrices,residualIdx,storageTypes,err,error,*999)
    DEALLOCATE(storageTypes)

    EXITS("EquationsMatricesVector_NonlinearStorageTypeSet0")
    RETURN
999 IF(ALLOCATED(storageTypes)) DEALLOCATE(storageTypes)
998 ERRORSEXITS("EquationsMatricesVector_NonlinearStorageTypeSet0",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_NonlinearStorageTypeSet0

  !
  !================================================================================================================================
  !

  !>Sets the storage type (sparsity) of the nonlinear (Jacobian) equations matrices
  SUBROUTINE EquationsMatricesVector_NonlinearStorageTypeSet1(vectorMatrices,residualIdx,storageTypes,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the vector eqautions matrices
    INTEGER(INTG), INTENT(IN) :: residualIdx !<The residual index to set the nonlinear storage types for
    INTEGER(INTG), INTENT(IN) :: storageTypes(:) !<storageTypes(matrixIdx). The storage type for the matrixIdx'th Jacobian equations matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("EquationsMatricesVector_NonlinearStorageTypeSet1",err,error,*999)

    CALL EquationsMatricesVector_AssertNotFinished(vectorMatrices,err,error,*999)
    NULLIFY(nonlinearMatrices)
    CALL EquationsMatricesVector_NonlinearMatricesGet(vectorMatrices,nonlinearMatrices,err,error,*999)
    NULLIFY(residualVector)
    CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,residualIdx,residualVector,err,error,*999)
    IF(SIZE(storageTypes,1)<residualVector%numberOfJacobians) THEN
      localError="The size of the storage type array of "//TRIM(NumberToVString(SIZE(storageTypes,1),"*",err,error))// &
        & " is invalid. The size should be >= "// &
        & TRIM(NumberToVString(residualVector%numberOfJacobians,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    DO matrixIdx=1,residualVector%numberOfJacobians
      NULLIFY(jacobianMatrix)
      CALL EquationsMatricesResidual_JacobianMatrixGet(residualVector,matrixIdx,jacobianMatrix,err,error,*999)
      SELECT CASE(storageTypes(matrixIdx))
      CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
        jacobianMatrix%storageType=DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE
      CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
        jacobianMatrix%storageType=DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE        
      CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
        jacobianMatrix%storageType=DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE
      CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
        jacobianMatrix%storageType=DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE
      CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
        jacobianMatrix%storageType=DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE
      CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
        jacobianMatrix%storageType=DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE
      CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
        jacobianMatrix%storageType=DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE
      CASE DEFAULT
        localError="The specified storage type of "//TRIM(NumberToVString(storageTypes(matrixIdx),"*",err,error))// &
          & " for Jacobian matrix number "//TRIM(NumberToVString(matrixIdx,"*",err,error))//" of residual number "// &
          & TRIM(NumberToVString(residualIdx,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDDO !matrixIdx
    
    EXITS("EquationsMatricesVector_NonlinearStorageTypeSet1")
    RETURN
999 ERRORSEXITS("EquationsMatricesVector_NonlinearStorageTypeSet1",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_NonlinearStorageTypeSet1

  !
  !================================================================================================================================
  !

  !>Sets the structure (sparsity) of the dynamic equations matrices
  SUBROUTINE EquationsMatricesVector_DynamicStructureTypeSet(vectorMatrices,structureTypes,err,error,*)
    
    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the vector equations matrices
    INTEGER(INTG), INTENT(IN) :: structureTypes(:) !<structureTypes(matrixIdx). The storage type for the matrixIdx'th dynamic equations matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatrixType), POINTER :: dynamicMatrix
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsMatricesVector_DynamicStructureTypeSet",err,error,*999)

    CALL EqutionsMatricesVector_AssertNotFinished(vectorMatrices,err,error,*999)
    NULLIFY(dynamicMatrices)
    CALL EquationsMatricesVector_DynamicMatricesGet(vectorMatrices,dynamicMatrices,err,error,*999)
    IF(SIZE(structureTypes,1)<dynamicMatrices%numberOfDynamicMatrices) THEN
      localError="The size of the structure types array of "//TRIM(NumberToVString(SIZE(structureTypes,1),"*",err,error))// &
        & " is invalid. The size should be >= "// &
        & TRIM(NumberToVString(dynamicMatrices%numberOfDynamicMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    DO matrixIdx=1,dynamicMatrices%numberOfDynamicMatrices
      NULLIFY(dynamicMatrix)
      CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,matrixIdx,dynamicMatrix,err,error,*999)
      SELECT CASE(structureTypes(matrixIdx))
      CASE(EQUATIONS_MATRIX_NO_STRUCTURE)
        dynamicMatrix%structureType=EQUATIONS_MATRIX_NO_STRUCTURE
      CASE(EQUATIONS_MATRIX_FEM_STRUCTURE)
        dynamicMatrix%structureType=EQUATIONS_MATRIX_FEM_STRUCTURE
      CASE(EQUATIONS_MATRIX_DIAGONAL_STRUCTURE)
        dynamicMatrix%structureType=EQUATIONS_MATRIX_DIAGONAL_STRUCTURE
      CASE(EQUATIONS_MATRIX_NODAL_STRUCTURE)
        dynamicMatrix%structureType=EQUATIONS_MATRIX_NODAL_STRUCTURE
      CASE DEFAULT
        localError="The specified strucutre type of "//TRIM(NumberToVString(structureTypes(matrixIdx),"*",err,error))// &
          & " for dynamic matrix number "//TRIM(NumberToVString(matrixIdx,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDDO !matrixIdx
     
    EXITS("EquationsMatricesVector_DynamicStructureTypeSet")
    RETURN
999 ERRORSEXITS("EquationsMatricesVector_DynamicStructureTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_DynamicStructureTypeSet
  
  !
  !================================================================================================================================
  !

  !>Sets the structure (sparsity) of a linear equations matrix
  SUBROUTINE EquationsMatricesVector_LinearStructureTypeSet0(vectorMatrices,structureType,err,error,*)
    
    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the vector equations matrices
    INTEGER(INTG), INTENT(IN) :: structureType !<The storage type for the linear equations matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsMatricesVector_LinearStructureTypeSet0",err,error,*999)

    CALL EquationsMatricesVector_LinearStructureTypeSet1(vectorMatrices,[structureType],err,error,*999)
    
    EXITS("EquationsMatricesVector_LinearStructureTypeSet0")
    RETURN
999 ERRORSEXITS("EquationsMatricesVector_LinearStructureTypeSet0",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_LinearStructureTypeSet0
  
  !
  !================================================================================================================================
  !

  !>Sets the structure (sparsity) of the linear equations matrices
  SUBROUTINE EquationsMatricesVector_LinearStructureTypeSet1(vectorMatrices,structureTypes,err,error,*)
    
    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the vector equations matrices
    INTEGER(INTG), INTENT(IN) :: structureTypes(:) !<structureTypes(matrixIdx). The storage type for the matrixIdx'th linear equations matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatrixType), POINTER :: linearMatrix
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsMatricesVector_LinearStructureTypeSet1",err,error,*999)

    CALL EquationsMatricesVector_AssertNotFinished(vectorMatrices,err,error,*999)
    NULLIFY(linearMatrices)
    CALL EquationsMatricesVector_LinearMatricesGet(vectorMatrices,linearMatrices,err,error,*999)
    IF(SIZE(structureTypes,1)<linearMatrices%numberOfLinearMatrices) THEN
      localError="The size of the structure types array of "//TRIM(NumberToVString(SIZE(structureTypes,1),"*",err,error))// &
        & " is invalid. The size should be >= "// &
        & TRIM(NumberToVString(linearMatrices%numberOfLinearMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    DO matrixIdx=1,linearMatrices%numberOfLinearMatrices
      NULLIFY(linearMatrix)
      CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,matrixIdx,linearMatrix,err,error,*999)
      SELECT CASE(structureTypes(matrixIdx))
      CASE(EQUATIONS_MATRIX_NO_STRUCTURE)
        linearMatrix%structureType=EQUATIONS_MATRIX_NO_STRUCTURE
      CASE(EQUATIONS_MATRIX_FEM_STRUCTURE)
        linearMatrix%structureType=EQUATIONS_MATRIX_FEM_STRUCTURE
      CASE(EQUATIONS_MATRIX_DIAGONAL_STRUCTURE)
        linearMatrix%structureType=EQUATIONS_MATRIX_DIAGONAL_STRUCTURE
      CASE(EQUATIONS_MATRIX_NODAL_STRUCTURE)
        linearMatrix%structureType=EQUATIONS_MATRIX_NODAL_STRUCTURE
      CASE DEFAULT
        localError="The specified strucutre type of "//TRIM(NumberToVString(structureTypes(matrixIdx),"*",err,error))// &
          & " for linear matrix number "//TRIM(NumberToVString(matrixIdx,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDDO !matrixIdx
    
    EXITS("EquationsMatricesVector_LinearStructureTypeSet1")
    RETURN
999 ERRORSEXITS("EquationsMatricesVector_LinearStructureTypeSet1",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_LinearStructureTypeSet1
  
  !
  !================================================================================================================================
  !

  !>Sets the structure (sparsity) of all nonlinear (Jacobian) equations matrices
  SUBROUTINE EquationsMatricesVector_NonlinearStructureTypeSet0(vectorMatrices,residualIdx,structureType,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the vector equations matrices
    INTEGER(INTG), INTENT(IN) :: residualIdx !<The residual index to set the Jacobian structure types for
    INTEGER(INTG), INTENT(IN) :: structureType !<The structure type for all Jacobian equations matrices
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG), ALLOCATABLE :: structureTypes(:)
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector

    ENTERS("EquationsMatricesVector_NonlinearStructureTypeSet0",err,error,*998)

    CALL EquationsMatricesVector_AssertNotFinished(vectorMatrices,err,error,*998)
    NULLIFY(nonlinearMatrices)
    CALL EquationsMatricesVector_NonlinearMatricesGet(vectorMatrices,nonlinearMatrices,err,error,*998)
    NULLIFY(residualVector)
    CALL EquationsMatricesNonlinear_ResisdualVectorGet(nonlinearMatrices,residualVector,err,error,*998)

    ALLOCATE(structureTypes(residualVector%numberOfJacobians),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate storage types.",err,error,*999)
    structureTypes=structureType
    CALL EquationsMatricesVector_NonlinearStructureTypeSet1(vectorMatrices,residualIdx,structureTypes,err,error,*999)
    DEALLOCATE(structureTypes)

    EXITS("EquationsMatricesVector_NonlinearStructureTypeSet0")
    RETURN
999 IF(ALLOCATED(structureTypes)) DEALLOCATE(structureTypes)
998 ERRORS("EquationsMatricesVector_NonlinearStructureTypeSet0",err,error)
    EXITS("EquationsMatricesVector_NonlinearStructureTypeSet0")
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_NonlinearStructureTypeSet0

  !
  !================================================================================================================================
  !

  !>Sets the structure (sparsity) of the nonlinear (Jacobian) equations matrices
  SUBROUTINE EquationsMatricesVector_NonlinearStructureTypeSet1(vectorMatrices,residualIdx,structureTypes,err,error,*)
    
    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the vector equations matrices
    INTEGER(INTG), INTENT(IN) :: residualIdx !<The residual index to set the Jacobian structure types for
    INTEGER(INTG), INTENT(IN) :: structureTypes(:) !<structureTypes(matrixIdx). The structure type for the matrixIdx'th Jacobian equations matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsMatricesVector_NonlinearStructureTypeSet1",err,error,*999)

    CALL EquationsMatricesVector_AssertNotFinished(vectorMatrices,err,error,*999)
    NULLIFY(nonlinearMatrices)
    CALL EquationsMatricesVector_NonlinearMatricesGet(vectorMatrices,nonlinearMatrices,err,error,*999)
    NULLIFY(residualVector)
    CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,residualIdx,residualVector,err,error,*999)
    IF(SIZE(structureTypes,1)<residualVector%numberOfJacobians) THEN
      localError="The size of the structure types array of "//TRIM(NumberToVString(SIZE(structureTypes,1),"*",err,error))// &
        & " is invalid. The size should be >= "// TRIM(NumberToVString(residualVector%numberOfJacobians,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    DO matrixIdx=1,residualVector%numberOfJacobians
      NULLIFY(jacobianMatrix)
      CALL EquationsMatricesResidual_JacobianMatrixGet(residualVector,matrixIdx,jacobianMatrix,err,error,*999)
      SELECT CASE(structureTypes(matrixIdx))
      CASE(EQUATIONS_MATRIX_NO_STRUCTURE)
        jacobianMatrix%structureType=EQUATIONS_MATRIX_NO_STRUCTURE
      CASE(EQUATIONS_MATRIX_FEM_STRUCTURE)
        jacobianMatrix%structureType=EQUATIONS_MATRIX_FEM_STRUCTURE
      CASE(EQUATIONS_MATRIX_DIAGONAL_STRUCTURE)
        jacobianMatrix%structureType=EQUATIONS_MATRIX_DIAGONAL_STRUCTURE
      CASE(EQUATIONS_MATRIX_NODAL_STRUCTURE)
        jacobianMatrix%structureType=EQUATIONS_MATRIX_NODAL_STRUCTURE
      CASE DEFAULT
        localError="The specified strucutre type of "//TRIM(NumberToVString(structureTypes(matrixIdx),"*",err,error))// &
          & " for Jacobian matrix number "//TRIM(NumberToVString(matrixIdx,"*",err,error))//" of residual number "// &
          & TRIM(NumberToVString(residualIdx,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDDO !matrixIdx
    
    EXITS("EquationsMatricesVector_NonlinearStructureTypeSet1")
    RETURN
999 ERRORS("EquationsMatricesVector_NonlinearStructureTypeSet1",err,error)
    EXITS("EquationsMatricesVector_NonlinearStructureTypeSet0")
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_NonlinearStructureTypeSet1

  !
  !================================================================================================================================
  !

  !>Finalise the scalar equations matrices and deallocate all memory.
  SUBROUTINE EquationsMatrices_ScalarFinalise(scalarMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesScalarType), POINTER :: scalarMatrices !<A pointer to the scalar equations matrices
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
   
    ENTERS("EquationsMatrices_ScalarFinalise",err,error,*999)

    IF(ASSOCIATED(scalarMatrices)) THEN
      DEALLOCATE(scalarMatrices)
    ENDIF
       
    EXITS("EquationsMatrices_ScalarFinalise")
    RETURN
999 ERRORSEXITS("EquationsMatrices_ScalarFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_ScalarFinalise

  !
  !================================================================================================================================
  !

  !>Initialise the scalar equations matrices for the scalar equations.
  SUBROUTINE EquationsMatrices_ScalarInitialise(scalarEquations,err,error,*)
    
    !Argument variables
    TYPE(EquationsScalarType), POINTER :: scalarEquations !<A pointer to the scalar equations to initialise the scalar equations matrices for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(EquationsMappingScalarType), POINTER :: scalarMapping
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("EquationsMatrices_ScalarInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(scalarEquations)) CALL FlagError("Scalar equations is not associated.",err,error,*998)
    IF(ASSOCIATED(scalarEquations%scalarMatrices)) &
      & CALL FlagError("Scalar equations matrices is already associated for this equations.",err,error,*998)
    NULLIFY(scalarMapping)
    CALL EquationsScalar_ScalarMappingGet(scalarEquations,scalarMapping,err,error,*998)
    CALL EquationsMappingScalar_AssertIsFinished(scalarMapping,err,error,*998)
    
    ALLOCATE(scalarEquations%scalarMatrices,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate scalar equations scalar equations matrices.",err,error,*999)
    scalarEquations%scalarMatrices%scalarEquations=>scalarEquations
    scalarEquations%scalarMatrices%scalarMatricesFinished=.FALSE.
    scalarEquations%scalarMatrices%scalarMapping=>scalarMapping
    NULLIFY(scalarEquations%scalarMatrices%solverMapping)
    NULLIFY(scalarEquations%scalarMatrices%functions)
    NULLIFY(scalarEquations%scalarMatrices%normMatrices)
    NULLIFY(scalarEquations%scalarMatrices%dotProductMatrices)
    NULLIFY(scalarEquations%scalarMatrices%quadraticMatrices)
       
    EXITS("EquationsMatrices_ScalarInitialise")
    RETURN
999 CALL EquationsMatrices_ScalarFinalise(scalarEquations%scalarMatrices,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMatrices_ScalarInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_ScalarInitialise

  !
  !================================================================================================================================
  !

  !>Finalise the vector equations matrices and deallocate all memory.
  SUBROUTINE EquationsMatrices_VectorFinalise(vectorMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the vector equations matrices
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
   
    ENTERS("EquationsMatrices_VectorFinalise",err,error,*999)

    IF(ASSOCIATED(vectorMatrices)) THEN
      CALL EquationsMatricesVector_DynamicFinalise(vectorMatrices%dynamicMatrices,err,error,*999)
      CALL EquationsMatricesVector_LinearFinalise(vectorMatrices%linearMatrices,err,error,*999)
      CALL EquationsMatricesVector_NonlinearFinalise(vectorMatrices%nonlinearMatrices,err,error,*999)
      CALL EquationsMatricesVector_RHSFinalise(vectorMatrices%rhsVector,err,error,*999)      
      CALL EquationsMatricesVector_SourcesFinalise(vectorMatrices%sourceVectors,err,error,*999)      
      DEALLOCATE(vectorMatrices)
    ENDIF
       
    EXITS("EquationsMatrices_VectorFinalise")
    RETURN
999 ERRORSEXITS("EquationsMatrices_VectorFinalise",err,error)
    RETURN 1
  END SUBROUTINE EquationsMatrices_VectorFinalise

  !
  !================================================================================================================================
  !

  !>Initialise the vector equations matrices for the vector equations.
  SUBROUTINE EquationsMatrices_VectorInitialise(vectorEquations,err,error,*)
    
     !Argument variables
    TYPE(EquationsVectorType), POINTER :: vectorEquations !<A pointer to the vector equations to initialise the vector equations matrices for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,numberOfGlobalRows,numberOfRows,totalNumberOfRows
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("EquationsMatrices_VectorInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(vectorEquations)) CALL FlagError("Vector equations is not associated.",err,error,*998)
    IF(ASSOCIATED(vectorEquations%vectorMatrices)) &
      & CALL FlagError("Vector equations matrices is already associated for this equations.",err,error,*998)
    NULLIFY(vectorMapping)
    CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*998)
    CALL EquationsMappingVector_AssertIsFinished(vectorMapping,err,error,*998)
    NULLIFY(lhsMapping)
    CALL EquationsMappingVector_LHSMappingGet(vectorMapping,lhsMapping,err,error,*998)
    CALL EquationsMappingLHS_NumberOfRowsGet(lhsMapping,numberOfRows,err,error,*998)
    CALL EquationsMappingLHS_TotalNumberOfRowsGet(lhsMapping,totalNumberOfRows,err,error,*998)
    CALL EquationsMappingLHS_NumberOfGlobalRowsGet(lhsMapping,numberOfGlobalRows,err,error,*998)
    
    ALLOCATE(vectorEquations%vectorMatrices,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate vector equations vector equations matrices.",err,error,*999)
    vectorEquations%vectorMatrices%vectorEquations=>vectorEquations
    vectorEquations%vectorMatrices%vectorMatricesFinished=.FALSE.
    vectorEquations%vectorMatrices%vectorMapping=>vectorMapping
    NULLIFY(vectorEquations%vectorMatrices%solverMapping)
    vectorEquations%vectorMatrices%numberOfRows=numberOfRows
    vectorEquations%vectorMatrices%totalNumberOfRows=totalNumberOfRows
    vectorEquations%vectorMatrices%numberOfGlobalRows=numberOfGlobalRows
    NULLIFY(vectorEquations%vectorMatrices%dynamicMatrices)
    NULLIFY(vectorEquations%vectorMatrices%linearMatrices)
    NULLIFY(vectorEquations%vectorMatrices%nonlinearMatrices)
    NULLIFY(vectorEquations%vectorMatrices%rhsVector)
    NULLIFY(vectorEquations%vectorMatrices%sourceVectors)            
    CALL EquationsMatricesVector_DynamicInitialise(vectorEquations%vectorMatrices,err,error,*999)            
    CALL EquationsMatricesVector_LinearInitialise(vectorEquations%vectorMatrices,err,error,*999)            
    CALL EquationsMatricesVector_NonlinearInitialise(vectorEquations%vectorMatrices,err,error,*999)            
    CALL EquationsMatricesVector_RHSInitialise(vectorEquations%vectorMatrices,err,error,*999)            
    CALL EquationsMatricesVector_SourcesInitialise(vectorEquations%vectorMatrices,err,error,*999)            
       
    EXITS("EquationsMatrices_VectorInitialise")
    RETURN
999 CALL EquationsMatrices_VectorFinalise(vectorEquations%vectorMatrices,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMatrices_VectorInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_VectorInitialise

  !
  !================================================================================================================================
  !

  !>Initialise the values of the equations matrices and vectors to the given value e.g., 0.0_DP
  SUBROUTINE EquationsMatricesVector_VectorValuesInitialise(vectorMatrices,selectionType,value,err,error,*)
    
    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the equations matrices to initialise the values for
    INTEGER(INTG), INTENT(IN) :: selectionType !<The selection of equations matrices to be initialised \see EquationsMatricesRoutines::SelectMatricesTypes,EquationsMatricesRoutines
    REAL(DP), INTENT(IN) :: value !<The value to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx,residualIdx,sourceIdx
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatricesSourcesType), POINTER :: sourceVectors
    TYPE(EquationsMatrixType), POINTER :: dynamicMatrix,linearMatrix
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix
    
    ENTERS("EquationsMatricesVector_VectorValuesInitialise",err,error,*999)
    
    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector equations matrices is not associated.",err,error,*999)
    
    IF(selectionType==EQUATIONS_MATRICES_ALL.OR. &
      & selectionType==EQUATIONS_MATRICES_DYNAMIC_ONLY.OR. &
      & selectionType==EQUATIONS_MATRICES_LINEAR_ONLY.OR. &
      & selectionType==EQUATIONS_MATRICES_NONLINEAR_ONLY) THEN
      NULLIFY(dynamicMatrices)
      CALL EquationsMatricesVector_DynamicMatricesExists(vectorMatrices,dynamicMatrices,err,error,*999)
      IF(ASSOCIATED(dynamicMatrices)) THEN
        DO matrixIdx=1,dynamicMatrices%numberOfDynamicMatrices
          NULLIFY(dynamicMatrix)
          CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,matrixIdx,dynamicMatrix,err,error,*999)
          IF(dynamicMatrix%updateMatrix) CALL DistributedMatrix_AllValuesSet(dynamicMatrix%matrix,value,err,error,*999)
        ENDDO !matrixIdx
      ENDIF
    ENDIF
    IF(selectionType==EQUATIONS_MATRICES_ALL.OR. &
      & selectionType==EQUATIONS_MATRICES_DYNAMIC_ONLY.OR. &
      & selectionType==EQUATIONS_MATRICES_LINEAR_ONLY.OR. &
      & selectionType==EQUATIONS_MATRICES_NONLINEAR_ONLY) THEN
      NULLIFY(linearMatrices)
      CALL EquationsMatricesVector_LinearMatricesExists(vectorMatrices,linearMatrices,err,error,*999)
      IF(ASSOCIATED(linearMatrices)) THEN
        DO matrixIdx=1,linearMatrices%numberOfLinearMatrices
          NULLIFY(linearMatrix)
          CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,matrixIdx,linearMatrix,err,error,*999)
          IF(linearMatrix%updateMatrix) CALL DistributedMatrix_AllValuesSet(linearMatrix%matrix,value,err,error,*999)
       ENDDO !matrixIdx
      ENDIF
    ENDIF
    IF(selectionType==EQUATIONS_MATRICES_ALL.OR. &
      & selectionType==EQUATIONS_MATRICES_NONLINEAR_ONLY.OR. &
      & selectionType==EQUATIONS_MATRICES_JACOBIAN_ONLY) THEN
      NULLIFY(nonlinearMatrices)
      CALL EquationsMatricesVector_NonlinearMatricesExists(vectorMatrices,nonlinearMatrices,err,error,*999)
      IF(ASSOCIATED(nonlinearMatrices)) THEN
        DO residualIdx=1,nonlinearMatrices%numberOfResiduals
          NULLIFY(residualVector)
          CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,residualIdx,residualVector,err,error,*999)
          DO matrixIdx=1,residualVector%numberOfJacobians
            NULLIFY(jacobianMatrix)
            CALL EquationsMatricesResidual_JacobianMatrixGet(residualVector,matrixIdx,jacobianMatrix,err,error,*999)
            IF(jacobianMatrix%updateJacobian) CALL DistributedMatrix_AllValuesSet(jacobianMatrix%jacobian,VALUE,err,error,*999)
          ENDDO !matrixIdx
        ENDDO !residualIdx
      ENDIF
    ENDIF
    IF(selectionType==EQUATIONS_MATRICES_ALL.OR. &
      & selectionType==EQUATIONS_MATRICES_NONLINEAR_ONLY.OR. &
      & selectionType==EQUATIONS_MATRICES_RESIDUAL_ONLY.OR. &
      & selectionType==EQUATIONS_MATRICES_RHS_RESIDUAL_ONLY.OR. &
      & selectionType==EQUATIONS_MATRICES_RESIDUAL_SOURCE_ONLY.OR. &
      & selectionType==EQUATIONS_MATRICES_VECTORS_ONLY) THEN
      NULLIFY(nonlinearMatrices)
      CALL EquationsMatricesVector_NonlinearMatricesExists(vectorMatrices,nonlinearMatrices,err,error,*999)
      IF(ASSOCIATED(nonlinearMatrices)) THEN
        DO residualIdx=1,nonlinearMatrices%numberOfResiduals
          NULLIFY(residualVector)
          CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,residualIdx,residualVector,err,error,*999)
          IF(residualVector%updateResidual) CALL DistributedVector_AllValuesSet(residualVector%residual,VALUE,err,error,*999)
        ENDDO !residualIdx
      ENDIF
    ENDIF
    IF(selectionType==EQUATIONS_MATRICES_ALL.OR. &
      & selectionType==EQUATIONS_MATRICES_DYNAMIC_ONLY.OR. &
      & selectionType==EQUATIONS_MATRICES_LINEAR_ONLY.OR. &
      & selectionType==EQUATIONS_MATRICES_NONLINEAR_ONLY.OR. &
      & selectionType==EQUATIONS_MATRICES_RHS_ONLY.OR. &
      & selectionType==EQUATIONS_MATRICES_RHS_RESIDUAL_ONLY.OR. &
      & selectionType==EQUATIONS_MATRICES_RHS_SOURCE_ONLY.OR. &
      & selectionType==EQUATIONS_MATRICES_VECTORS_ONLY) THEN
      NULLIFY(rhsVector)
      CALL EquationsMatricesVector_RHSVectorExists(vectorMatrices,rhsVector,err,error,*999)
      IF(ASSOCIATED(rhsVector)) THEN
        IF(rhsVector%updateVector) CALL DistributedVector_AllValuesSet(rhsVector%vector,value,err,error,*999)
      ENDIF
    ENDIF
    IF(selectionType==EQUATIONS_MATRICES_ALL.OR. &
      & selectionType==EQUATIONS_MATRICES_DYNAMIC_ONLY.OR. &
      & selectionType==EQUATIONS_MATRICES_LINEAR_ONLY.OR. &
      & selectionType==EQUATIONS_MATRICES_NONLINEAR_ONLY.OR. &
      & selectionType==EQUATIONS_MATRICES_SOURCE_ONLY.OR. &
      & selectionType==EQUATIONS_MATRICES_RHS_SOURCE_ONLY.OR. &
      & selectionType==EQUATIONS_MATRICES_RESIDUAL_SOURCE_ONLY.OR. &
      & selectionType==EQUATIONS_MATRICES_VECTORS_ONLY) THEN    
      NULLIFY(sourceVector)
      CALL EquationsMatricesVector_SourceVectorsExists(vectorMatrices,sourceVectors,err,error,*999)
      IF(ASSOCIATED(sourceVectors)) THEN
        DO sourceIdx=1,sourceVectors%numberOfSources
          NULLIFY(sourceVector)
          CALL EquationsMatricesSources_SourceVectorGet(sourceVectors,sourceIdx,sourceVector,err,error,*999)
          IF(sourceVector%updateVector) CALL DistributedVector_AllValuesSet(sourceVector%vector,VALUE,err,error,*999)
        ENDDO !sourceIdx
      ENDIF
    ENDIF
 
    EXITS("EquationsMatricesVector_VectorValuesInitialise")
    RETURN
999 ERRORSEXITS("EquationsMatricesVector_VectorValuesInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_VectorValuesInitialise

  !
  !================================================================================================================================
  !

  !>Finalise a residual vector and deallocate all memory.
  SUBROUTINE EquationsMatricesNonlinear_ResidualFinalise(residualVector,err,error,*)
    
     !Argument variables
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector !<A pointer to the residual vector to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    
    ENTERS("EquationsMatricesNonlinear_ResidualFinalise",err,error,*999)

    IF(ASSOCIATED(residualVector)) THEN
      IF(ASSOCIATED(residualVector%residual)) CALL DistributedVector_Destroy(residualVector%residual,err,error,*999)
      IF(ASSOCIATED(residualVector%previousResidual)) CALL DistributedVector_Destroy(residualVector%previousResidual, &
        & err,error,*999)
      IF(ASSOCIATED(residualVector%previous2Residual)) CALL DistributedVector_Destroy(residualVector%previous2Residual, &
        & err,error,*999)
      IF(ASSOCIATED(residualVector%previous3Residual)) CALL DistributedVector_Destroy(residualVector%previous3Residual, &
        & err,error,*999)
      CALL EquationsMatrices_ElementVectorFinalise(residualVector%elementResidual,err,error,*999)
      CALL EquationsMatrices_NodalVectorFinalise(residualVector%nodalResidual,err,error,*999)
      IF(ALLOCATED(residualVector%jacobians)) THEN
        DO matrixIdx=1,SIZE(residualVector%jacobians,1)
          CALL EquationsMatricesResidual_JacobianFinalise(residualVector%jacobians(matrixIdx)%ptr,err,error,*999)
        ENDDO !matrixIdx
        DEALLOCATE(residualVector%jacobians)
      ENDIF
      DEALLOCATE(residualVector)
    ENDIF
       
    EXITS("EquationsMatricesNonlinear_ResidualFinalise")
    RETURN
999 ERRORSEXITS("EquationsMatricesNonlinear_ResidualFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesNonlinear_ResidualFinalise

  !
  !================================================================================================================================
  !

  !>Allocate and initialise a residual vector for the nonlinear matrices
  SUBROUTINE EquationsMatricesNonlinear_ResidualInitialise(nonlinearMatrices,residualIdx,err,error,*)
    
     !Argument variables
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices !<A pointer to the nonlinear matrices to initialise the residual vector for
    INTEGER(INTG), INTENT(IN) :: residualIdx !<The residual index to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,matrixIdx,numberOfJacobians
    REAL(DP) :: residualCoefficient
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector
    TYPE(VARYING_STRING) :: dummyError,localError
    
    ENTERS("EquationsMatricesNonlinear_ResidualInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(nonlinearMatrices)) CALL FlagError("Nonlinear matrices is not associated.",err,error,*998)
    IF(residualIdx<1.OR.residualIdx>nonlinearMatrices%numberOfResiduals) THEN
      localError="The specified residual index of "//TRIM(NumberToVString(residualIdx,"*",err,error))// &
        & " is invalid. The residual index should be >= 1 and <= "// &
        & TRIM(NumberToVString(nonlinearMatrices%numberOfResiduals,"*",err,error))//"."
      CALL FlagError(localError,err,error,*998)
    ENDIF
    IF(.NOT.ALLOCATED(nonlinearMatrices%residuals)) &
      & CALL FlagError("The residuals are not allocated for the nonlinear matrices.",err,error,*998)
    IF(ASSOCIATED(nonlinearMatrices%residuals(residualIdx)%ptr)) THEN
      localError="Residual vector is already associated for residual index "// &
        & TRIM(NumberToVString(residualIdx,"*",err,error))//" of the nonlinear matrices."
      CALL FlagError(localError,err,error,*998)
    ENDIF
    NULLIFY(nonlinearMapping)
    CALL EquationsMatricesNonlinear_NonlinearMappingGet(nonlinearMatrices,nonlinearMapping,err,error,*998)
    NULLIFY(residualMapping)
    CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,residualIdx,residualMapping,err,error,*998)
    CALL EquationsMappingResidual_VectorCoefficientGet(residualMapping,residualCoefficient,err,error,*998)
    !Allocate and initialise the residual vector
    ALLOCATE(nonlinearMatrices%residuals(residualIdx)%ptr,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate nonlinear matrices residual vector.",err,error,*999)
    NULLIFY(residualVector)
    CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,residualIdx,residualVector,err,error,*999)
    residualVector%nonlinearMatrices=>nonlinearMatrices
    residualVector%residualNumber=residualIdx
    residualVector%residualCoefficient=residualCoefficient
    residualVector%updateResidual=.TRUE.
    residualVector%firstAssembly=.TRUE.
    NULLIFY(residualVector%residual)
    NULLIFY(residualVector%previousResidual)
    NULLIFY(residualVector%previous2Residual)
    NULLIFY(residualVector%previous3Residual)
    CALL EquationsMatrices_ElementVectorInitialise(residualVector%elementResidual,err,error,*999)
    CALL EquationsMatrices_NodalVectorInitialise(residualVector%nodalResidual,err,error,*999)
    !Allocate and initialise the Jacobians
    CALL EquationsMappingResidual_NumberOfJacobianMatricesGet(residualMapping,numberOfJacobians,err,error,*999)
    ALLOCATE(residualVector%jacobians(numberOfJacobians),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate Jacobians for the residual vector.",err,error,*999)
    residualVector%numberOfJacobians=numberOfJacobians
    DO matrixIdx=1,numberOfJacobians
      NULLIFY(residualVector%jacobians(matrixIdx)%ptr)
      CALL EquationsMatricesResidual_JacobianInitialise(residualVector,matrixIdx,err,error,*999)
    ENDDO !matrixIdx
       
    EXITS("EquationsMatricesNonlinear_ResidualInitialise")
    RETURN
999 CALL EquationsMatricesNonlinear_ResidualFinalise(nonlinearMatrices%residuals(residualIdx)%ptr,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMatricesNonlinear_ResidualInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesNonlinear_ResidualInitialise

  !
  !================================================================================================================================
  !

  !>Finalise the equations Jacobian and deallocate all memory
  SUBROUTINE EquationsMatricesResidual_JacobianFinalise(jacobianMatrix,err,error,*)

    !Argument variables
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix !<A pointer to the equations Jacobian to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("EquationsMatricesResidual_JacobianFinalise",err,error,*999)

    IF(ASSOCIATED(jacobianMatrix)) THEN
      IF(ASSOCIATED(jacobianMatrix%jacobian)) CALL DistributedMatrix_Destroy(jacobianMatrix%jacobian,err,error,*999)
      CALL EquationsMatrices_ElementMatrixFinalise(jacobianMatrix%elementJacobian,err,error,*999)
      CALL EquationsMatrices_NodalMatrixFinalise(jacobianMatrix%nodalJacobian,err,error,*999)
      DEALLOCATE(jacobianMatrix)
    ENDIF
    
    EXITS("EquationsMatricesResidual_JacobianFinalise")
    RETURN
999 ERRORSEXITS("EquationsMatricesResidual_JacobianFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesResidual_JacobianFinalise

  !
  !================================================================================================================================
  !

  !>Initialise the equations Jacobian matrix.
  SUBROUTINE EquationsMatricesResidual_JacobianInitialise(residualVector,matrixIdx,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector !<A pointer to the residual vector to initialise the Jacobian matri for
    INTEGER(INTG), INTENT(IN) :: matrixIdx !<The matrix index to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,numberOfColumns
    REAL(DP) :: jacobianCoefficient
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix
    TYPE(JacobianMatrixToVarMapType), POINTER :: jacobianMatrixToVarMap
    TYPE(VARYING_STRING) :: dummyError,localError

    ENTERS("EquationsMatricesResidual_JacobianInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(residualVector)) CALL FlagError("ResidualVector is not associated.",err,error,*998)
    IF(matrixIdx<1.OR.matrixIdx>residualVector%numberOfJacobians) THEN
      localError="The specified Jacobian matrix index of "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
        & " is invalid. The Jacobian matrix index should be >= 1 and <= "// &
        & TRIM(NumberToVString(residualVector%numberOfJacobians,"*",err,error))//"."
      CALL FlagError(localError,err,error,*998)
    ENDIF
    IF(.NOT.ALLOCATED(residualVector%jacobians)) &
      & CALL FlagError("The Jacobians are not allocated for the residual vector.",err,error,*998)
    IF(ASSOCIATED(residualVector%jacobians(matrixIdx)%ptr)) THEN
      localError="Jacobian matri is already associated for matrix index "// &
        & TRIM(NumberToVString(matrixIdx,"*",err,error))//" of the residual vector."
      CALL FlagError(localError,err,error,*998)
    ENDIF
    NULLIFY(residualMapping)
    CALL EquationsMatricesResidual_ResidualMappingGet(residualVector,residualMapping,err,error,*998)
    NULLIFY(jacobianMatrixToVarMap)
    CALL EquationsMappingResidual_JacobianMatrixToVarMapGet(residualMapping,matrixIdx,jacobianMatrixToVarMap,err,error,*999)
    CALL EquationsMappingVectorJMToVMap_MatrixCoefficientGet(jacobianMatrixToVarMap,jacobianCoefficient,err,error,*999)
    CALL EquationsMappingVectorJMToVMap_NumberOfColumnsGet(jacobianMatrixToVarMap,numberOfColumns,err,error,*999)

    ALLOCATE(residualVector%jacobians(matrixIdx)%ptr,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate Jacobian matrix.",err,error,*999)
    NULLIFY(jacobianMatrix)
    CALL EquationsMatricesResidual_JacobianMatrixGet(residualVector,matrixIdx,jacobianMatrix,err,error,*999)
    jacobianMatrix%residualVector=>residualVector
    jacobianMatrix%jacobianNumber=0
    jacobianMatrix%jacobianCoefficient=jacobianCoefficient
    jacobianMatrix%storageType=MATRIX_BLOCK_STORAGE_TYPE
    jacobianMatrix%structureType=EQUATIONS_MATRIX_NO_STRUCTURE
    jacobianMatrix%symmetric=.FALSE.
    jacobianMatrix%numberOfColumns=numberOfColumns
    jacobianMatrix%updateJacobian=.TRUE.
    jacobianMatrix%firstAssembly=.TRUE.
    NULLIFY(jacobianMatrix%jacobian)
    CALL EquationsMatrices_ElementMatrixInitialise(jacobianMatrix%elementJacobian,err,error,*999)
    CALL EquationsMatrices_NodalMatrixInitialise(jacobianMatrix%nodalJacobian,err,error,*999)
    jacobianMatrix%jacobianCalculationType=EQUATIONS_JACOBIAN_FINITE_DIFFERENCE_CALCULATED
    jacobianMatrix%jacobianFiniteDifferenceStepSize=1.0E-6_DP
    jacobianMatrixToVarMap%jacobian=>jacobianMatrix
    
    EXITS("EquationsMatricesResidual_JacobianInitialise")
    RETURN
999 CALL EquationsMatricesResidual_JacobianFinalise(residualVector%jacobians(matrixIdx)%ptr,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMatricesResidual_JacobianInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesResidual_JacobianInitialise

  !
  !================================================================================================================================
  !

  !>Finalise a source vector and deallocate all memory.
  SUBROUTINE EquationsMatricesSources_SourceVectorFinalise(sourceVector,err,error,*)
    
     !Argument variables
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector !<A pointer to the source vector to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    
    ENTERS("EquationsMatricesSources_SourceVectorFinalise",err,error,*999)

    IF(ASSOCIATED(sourceVector)) THEN
      IF(ASSOCIATED(sourceVector%vector)) CALL DistributedVector_Destroy(sourceVector%vector,err,error,*999)
      IF(ASSOCIATED(sourceVector%previousSourceVector)) CALL DistributedVector_Destroy(sourceVector%previousSourceVector, &
        & err,error,*999)
      IF(ASSOCIATED(sourceVector%previous2SourceVector)) CALL DistributedVector_Destroy(sourceVector%previous2SourceVector, &
        & err,error,*999)
      IF(ASSOCIATED(sourceVector%previous3SourceVector)) CALL DistributedVector_Destroy(sourceVector%previous3SourceVector, &
        & err,error,*999)
      CALL EquationsMatrices_ElementVectorFinalise(sourceVector%elementVector,err,error,*999)
      CALL EquationsMatrices_NodalVectorFinalise(sourceVector%nodalVector,err,error,*999)
      DEALLOCATE(sourceVector)
    ENDIF
       
    EXITS("EquationsMatricesSources_SourceVectorFinalise")
    RETURN
999 ERRORSEXITS("EquationsMatricesSources_SourceVectorFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesSources_SourceVectorFinalise

  !
  !================================================================================================================================
  !

  !>Allocate and initialise a source vector for the source vectors
  SUBROUTINE EquationsMatricesSources_SourceVectorInitialise(sourceVectors,sourceIdx,err,error,*)
    
     !Argument variables
    TYPE(EquationsMatricesSourcesType), POINTER :: sourceVectors !<A pointer to the source vectors to initialise the source vector for
    INTEGER(INTG), INTENT(IN) :: sourceIdx !<The source index to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    REAL(DP) :: sourceCoefficient
    TYPE(EquationsMappingSourceType), POINTER :: sourceMapping
    TYPE(EquationsMappingSourcesType), POINTER :: sourcesMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(VARYING_STRING) :: dummyError,localError
    
    ENTERS("EquationsMatricesSources_SourceVectorInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(sourceVectors)) CALL FlagError("Source vectors is not associated.",err,error,*998)
    IF(sourceIdx<1.OR.sourceIdx>sourceVectors%numberOfSources) THEN
      localError="The specified source index of "//TRIM(NumberToVString(sourceIdx,"*",err,error))// &
        & " is invalid. The source index should be >= 1 and <= "// &
        & TRIM(NumberToVString(sourceVectors%numberOfSources,"*",err,error))//"."
      CALL FlagError(localError,err,error,*998)
    ENDIF
    IF(.NOT.ALLOCATED(sourceVectors%sources)) &
      & CALL FlagError("The sources are not allocated for the source vectors.",err,error,*998)
    IF(ASSOCIATED(sourceVectors%sources(sourceIdx)%ptr)) THEN
      localError="Source vector is already associated for source index "// &
        & TRIM(NumberToVString(sourceIdx,"*",err,error))//" of the source vectors."
      CALL FlagError(localError,err,error,*998)
    ENDIF
    NULLIFY(sourcesMapping)
    CALL EquationsMatricesSources_SourcesMappingGet(sourceVectors,sourcesMapping,err,error,*998)
    NULLIFY(sourceMapping)
    CALL EquationsMappingSources_SourceMappingGet(sourcesMapping,sourceIdx,sourceMapping,err,error,*998)
    CALL EquationsMappingSource_VectorCoefficientGet(sourceMapping,sourceCoefficient,err,error,*999)
    !Allocate and initialise the residual vector
    ALLOCATE(sourceVectors%sources(sourceIdx)%ptr,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate source vectors source vector.",err,error,*999)
    NULLIFY(sourceVector)
    CALL EquationsMatricesSources_SourceVectorGet(sourceVectors,sourceIdx,sourceVector,err,error,*999)
    sourceVectors%sources(sourceIdx)%ptr%sources=>sourceVectors
    sourceVectors%sources(sourceIdx)%ptr%sourceNumber=sourceIdx
    sourceVectors%sources(sourceIdx)%ptr%sourceCoefficient=sourceCoefficient
    sourceVectors%sources(sourceIdx)%ptr%updateVector=.TRUE.
    sourceVectors%sources(sourceIdx)%ptr%firstAssembly=.TRUE.
    NULLIFY(sourceVectors%sources(sourceIdx)%ptr%vector)
    CALL EquationsMatrices_ElementVectorInitialise(sourceVectors%sources(sourceIdx)%ptr%elementVector,err,error,*999)
    CALL EquationsMatrices_NodalVectorInitialise(sourceVectors%sources(sourceIdx)%ptr%nodalVector,err,error,*999)
       
    EXITS("EquationsMatricesSources_SourceVectorInitialise")
    RETURN
999 CALL EquationsMatricesSources_SourceVectorFinalise(sourceVectors%sources(sourceIdx)%ptr,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMatricesSources_SourceVectorInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesSources_SourceVectorInitialise

  !
  !================================================================================================================================
  !

  !>Outputs the equations matrices
  SUBROUTINE EquationsMatricesVector_Output(id,vectorMatrices,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the ouptut stream
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the equations matrices
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx,residualIdx,sourceIdx
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatricesSourcesType), POINTER :: sourceVectors
    TYPE(EquationsMatrixType), POINTER :: dynamicMatrix,linearMatrix
    
    ENTERS("EquationsMatricesVector_Output",err,error,*999)

    CALL EquationsMatricesVector_AssertIsFinalised(vectorMatrices,err,error,*999)
     
    CALL WriteString(id,"",err,error,*999)
    CALL WriteString(id,"Equations matrices:",err,error,*999)
    NULLIFY(dynamicMatrices)
    CALL EquationsMatricesVector_DynamicMatricesExists(vectorMatrices,dynamicMatrices,err,error,*999)
    IF(ASSOCIATED(dynamicMatrices)) THEN
      CALL WriteString(id,"Dynamic matrices:",err,error,*999)
      CALL WriteStringValue(id,"Number of dynamic matrices = ",dynamicMatrices%numberOfDynamicMatrices,err,error,*999)
      DO matrixIdx=1,dynamicMatrices%numberOfDynamicMatrices
        NULLIFY(dynamicMatrix)
        CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,matrixIdx,dynamicMatrix,err,error,*999)
        CALL WriteStringValue(id,"Dynamic matrix : ",matrixIdx,err,error,*999)
        CALL DistributedMatrix_Output(id,dynamicMatrix%matrix,err,error,*999)
      ENDDO !matrixIdx
    ENDIF
    NULLIFY(linearMatrices)
    CALL EquationsMatricesVector_LinearMatricesExists(vectorMatrices,linearMatrices,err,error,*999)
    IF(ASSOCIATED(linearMatrices)) THEN
      CALL WriteString(id,"Linear matrices:",err,error,*999)
      CALL WriteStringValue(id,"Number of linear matrices = ",linearMatrices%numberOfLinearMatrices,err,error,*999)
      DO matrixIdx=1,linearMatrices%numberOfLinearMatrices
        NULLIFY(linearMatrix)
        CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,matrixIdx,linearMatrix,err,error,*999)
        CALL WriteStringValue(id,"Linear matrix : ",matrixIdx,err,error,*999)
        CALL DistributedMatrix_Output(id,linearMatrix%matrix,err,error,*999)
      ENDDO !matrixIdx
    ENDIF
    NULLIFY(nonlinearMatrices)
    CALL EquationsMatricesVector_NonlinearMatricesExists(vectorMatrices,nonlinearMatrices,err,error,*999)
    IF(ASSOCIATED(nonlinearMatrices)) THEN
      CALL WriteString(id,"Nonlinear vectors:",err,error,*999)
      DO residualIdx=1,nonlinearMatrices%numberOfResiduals
        NULLIFY(residualVector)
        CALL EqutionsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,residualIdx,residualVector,err,error,*999)
        CALL WriteStringValue(id,"Residual vector : ",residualIdx,err,error,*999)
        CALL DistributedVector_Output(id,residualVector%residual,err,error,*999)
      ENDDO !residualIdx
    ENDIF
    NULLIFY(rhsVector)
    CALL EquationsMatricesVector_RHSVectorExists(vectorMatrices,rhsVector,err,error,*999)
    IF(ASSOCIATED(rhsVector)) THEN
      CALL WriteString(id,"RHS vector:",err,error,*999)
      CALL DistributedVector_Output(id,rhsVector%vector,err,error,*999)
    ENDIF
    NULLIFY(sourceVectors)
    CALL EquationsMatricesVector_SourceVectorsExists(vectorMatrices,sourceVectors,err,error,*999)
    IF(ASSOCIATED(sourceVectors)) THEN
      CALL WriteString(id,"Source vectors:",err,error,*999)
      DO sourceIdx=1,sourceVectors%numberOfSources
        NULLIFY(sourceVector)
        CALL EquationsMatricesSources_SourceVectorGet(sourceVectors,sourceIdx,sourceVector,err,error,*999)
        CALL WriteStringValue(id,"Source vector : ",sourceIdx,err,error,*999)
        CALL DistributedVector_Output(id,sourceVector%vector,err,error,*999)
      ENDDO !sourceIdx
    ENDIF
    
    EXITS("EquationsMatricesVector_Output")
    RETURN
999 ERRORSEXITS("EquationsMatricesVector_Output",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_Output
  
  !
  !================================================================================================================================
  !

  !>Caclulates the matrix structure (sparsity) for an equations matrix.
  SUBROUTINE EquationsMatrix_StructureCalculate(equationsMatrix,numberOfNonZeros,rowIndices,columnIndices,list,err,error,*)

    !Argument variables
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix !<A pointer to the equations matrix to calculate the structure for
    INTEGER(INTG), INTENT(OUT) :: numberOfNonZeros !<On return the number of non-zeros in the matrix
    INTEGER(INTG), POINTER :: rowIndices(:) !<On return a pointer to row location indices in compressed row format. The pointer must be NULL on entry and the calling routine is responsible for deallocation.
    INTEGER(INTG), POINTER :: columnIndices(:) !<On return a pointer to the column location indices in compressed row format. The pointer must be NULL on entry and the calling routine is responsible for deallocation.
    TYPE(LinkedList), POINTER :: list(:) 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: columnIdx,componentIdx,componentNumber,derivativeIdx,derivativeNumber,derivativeNumber2,dofIdx,dofType, &
      & dofTypeIdx,dummyErr,elementIdx,elementNumber,globalColumn,interpolationType,localColumn,localDOF,localDOFIdx, &
      & localNodeIdx,matrixNumber,maxElementInterpParameters,nodeNumber,nodeNumber2,numberOfColumns,numberOfColumnComponents, &
      & numberOfDerivatives,numberOfGlobalColumns,numberOfLocalNodes,numberOfNodeDerivatives,numberOfRowComponents, &
      & numberOfSurroundingElements,numberOfVersions,residualNumber,totalNumberOfColumns,totalNumberOfRows,versionIdx, &
      & versionNumber,versionNumber2
    INTEGER(INTG), ALLOCATABLE :: columns(:)
    REAL(DP) :: sparsity
    LOGICAL :: boundaryElement,boundaryNode
    TYPE(BasisType), POINTER :: colsBasis,rowsBasis
    TYPE(DomainType), POINTER :: colsDomain,rowsDomain
    TYPE(DomainElementsType), POINTER :: colsDomainElements,rowsDomainElements
    TYPE(DomainMappingType), POINTER :: colsDOFsDomainMapping,rowsDOFsDomainMapping
    TYPE(DomainNodesType), POINTER :: colsDomainNodes,rowsDomainNodes
    TYPE(DomainTopologyType), POINTER :: colsDomainTopology,rowsDomainTopology
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector
    TYPE(EquationsMatrixToVarMapType), POINTER :: equationsMatrixToVarMap
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(FieldType), POINTER :: dependentField
    TYPE(FieldDOFToParamMapType), POINTER :: dependentDofsParamMapping
    TYPE(FieldVariableType), POINTER :: colsVariable,rowsVariable
    TYPE(ListPtrType), ALLOCATABLE :: columnIndicesLists(:)
    TYPE(VARYING_STRING) :: dummyError,localError
    
    ENTERS("EquationsMatrix_StructureCalculate",err,error,*998)

    numberOfNonZeros=0
    
    IF(.NOT.ASSOCIATED(equationsMatrix)) CALL FlagError("Equations matrix is not associated.",err,error,*998)
    IF(ASSOCIATED(rowIndices)) CALL FlagError("Row indices is already associated.",err,error,*998)
    IF(ASSOCIATED(columnIndices)) CALL FlagError("Column indices is already associated.",err,error,*998)
    NULLIFY(linearMatrices)
    CALL EquationsMatrix_LinearMatricesExists(equationsMatrix,linearMatrices,err,error,*999)
    NULLIFY(dynamicMatrices)
    CALL EquationsMatrix_DynamicMatricesExists(equationsMatrix,dynamicMatrices,err,error,*999)
    IF(.NOT.ASSOCIATED(dynamicMatrices).AND..NOT.ASSOCIATED(linearMatrices)) &
      & CALL FlagError("Neither of equations matrix dynamic or linear matrices is associated.",err,error,*998)
    NULLIFY(vectorMatrices)
    IF(ASSOCIATED(dynamicMatrices)) THEN
      CALL EquationsMatricesDynamic_VectorMatricesGet(dynamicMatrices,vectorMatrices,err,error,*999)
    ELSE
      CALL EquationsMatricesLinear_VectorMatricesGet(linearMatrices,vectorMatrices,err,error,*999)
    ENDIF
    NULLIFY(vectorMapping)
    CALL EquationsMatricesVector_VectorMappingGet(vectorMatrices,vectorMapping,err,error,*998)
    NULLIFY(vectorEquations)
    CALL EquationsMatricesVector_VectorEquationsGet(vectorMatrices,vectorEquations,err,error,*998)
    NULLIFY(equations)
    CALL EquationsVector_EquationsGet(vectorEquations,equations,err,error,*998)
    NULLIFY(equationsSet)
    CALL Equations_EquationsSetGet(equations,equationsSet,err,error,*998)
    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*998)
    CALL EquationsMatrix_MatrixNumberGet(equationsMatrix,matrixNumber,err,error,*998)
    NULLIFY(lhsMapping)
    CALL EquationsMappingVector_LHSMappingGet(vectorMapping,lhsMapping,err,error,*998)
    NULLIFY(rowsVariable)
    CALL EquationsMappingLHS_LHSVariableGet(lhsMapping,rowsVariable,err,error,*998)
    NULLIFY(dynamicMapping)
    NULLIFY(linearMapping)
    NULLIFY(colsVariable)
    NULLIFY(equationsMatrixToVarMap)
    IF(ASSOCIATED(dynamicMatrices)) THEN
      CALL EquationsMappingVector_DynamicMappingGet(vectorMapping,dynamicMapping,err,error,*999)
      CALL EquationsMappingDynamic_DynamicVariableGet(dynamicMapping,colsVariable,err,error,*999)
      CALL EquationsMappingDynamic_EquationsMatrixToVarMapGet(dynamicMapping,matrixNumber,equationsMatrixToVarMap, err,error,*999)
    ELSE
      CALL EquationsMappingVector_LinearMappingGet(vectorMapping,linearMapping,err,error,*999)
      CALL EquationsMappingLinear_LinearMatrixVariableGet(linearMapping,matrixNumber,colsVariable,err,error,*999)
      CALL EquationsMappingLinear_EquationsMatrixToVarMapGet(linearMapping,matrixNumber,equationsMatrixToVarMap,err,error,*999)
    ENDIF
    CALL FieldVariable_NumberOfComponentsGet(rowsVariable,numberOfRowComponents,err,error,*999)
    CALL FieldVariable_NumberOfComponentsGet(colsVariable,numberOfColumnComponents,err,error,*999)
    NULLIFY(rowsDOFsDomainMapping)
    CALL FieldVariable_DomainMappingGet(rowsVariable,rowsDOFsDomainMapping,err,error,*999)
    CALL DomainMapping_TotalNumberOfLocalGet(rowsDOFsDomainMapping,totalNumberOfRows,err,error,*999)
    NULLIFY(colsDOFsDomainMapping)
    CALL FieldVariable_DomainMappingGet(colsVariable,colsDOFsDomainMapping,err,error,*999)
    CALL DomainMapping_TotalNumberOfLocalGet(colsDOFsDomainMapping,totalNumberOfColumns,err,error,*999)
    CALL DomainMapping_NumberOfGlobalGet(colsDOFsDomainMapping,numberOfGlobalColumns,err,error,*999)
    
    SELECT CASE(equationsMatrix%structureType)
    CASE(EQUATIONS_MATRIX_NO_STRUCTURE)
      CALL FlagError("There is no structure to calculate for a matrix with no structure.",err,error,*998)
    CASE(EQUATIONS_MATRIX_FEM_STRUCTURE)
      SELECT CASE(equationsMatrix%storageType)
      CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
        
        !Allocate lists
        ALLOCATE(columnIndicesLists(totalNumberOfRows),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate column indices lists.",err,error,*999)
        !Allocate row indices
        ALLOCATE(rowIndices(totalNumberOfRows+1),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate row indices.",err,error,*999)
        rowIndices(1)=1
        
        !First, loop over the rows and calculate the number of non-zeros
        numberOfNonZeros=0
        DO localDOFIdx=1,totalNumberOfRows
          CALL FieldVariable_DOFTypeGet(rowsVariable,localDOFIdx,dofType,dofTypeIdx,err,error,*999)          
          SELECT CASE(dofType)
          CASE(FIELD_CONSTANT_DOF_TYPE)
            CALL FlagError("Constant based DOFs is not implemented.",err,error,*999)
          CASE(FIELD_NODE_DOF_TYPE)
            CALL FieldVariable_DOFParameterGetNode(rowsVariable,dofTypeIdx,versionNumber,derivativeNumber,nodeNumber, &
              & componentNumber,err,error,*999)
            NULLIFY(rowsDomain)
            CALL FieldVariable_ComponentDomainGet(rowsVariable,componentNumber,rowsDomain,err,error,*999)
            CALL FieldVariable_ComponentMaxElementInterpParametersGet(rowsVariable,componentNumber,maxElementInterpParameters, &
              & err,error,*999)
            NULLIFY(rowsDomainTopology)
            CALL Domain_DomainTopologyGet(rowsDomain,rowsDomainTopology,err,error,*999)
            NULLIFY(rowsDomainNodes)
            CALL DomainTopology_DomainNodesGet(rowsDomainTopology,rowsDomainNodes,err,error,*999)
            CALL DomainNodes_NodeNumberOfSurroundingElementsGet(rowsDomainNodes,nodeNumber,numberOfSurroundingElements, &
              & err,error,*999)
            !Set up list
            NULLIFY(columnIndicesLists(localDOFIdx)%ptr)
            CALL List_CreateStart(columnIndicesLists(localDOFIdx)%ptr,err,error,*999)
            CALL List_DataTypeSet(columnIndicesLists(localDOFIdx)%ptr,LIST_INTG_TYPE,err,error,*999)
            CALL List_InitialSizeSet(columnIndicesLists(localDOFIdx)%ptr,numberOfSurroundingElements*maxElementInterpParameters, &
              & err,error,*999)
            CALL List_CreateFinish(columnIndicesLists(localDOFIdx)%ptr,err,error,*999)
            !Loop over all elements containing the dof
            DO elementIdx=1,numberOfSurroundingElements
              CALL DomainNodes_NodeSurroundingElementGet(rowsDomainNodes,elementIdx,nodeNumber,elementNumber,err,error,*999)
              DO componentIdx=1,numberOfColumnComponents
                CALL FieldVariable_ComponentInterpolationGet(colsVariable,componentIdx,interpolationType,err,error,*999)
                SELECT CASE(interpolationType)
                CASE(FIELD_CONSTANT_INTERPOLATION)
                  ! do nothing? this will probably never be encountered...?
                  CALL FlagError("Not implemented?",err,error,*999)
                CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                  CALL FieldVariable_LocalElementDOFGet(colsVariable,elementNumber,componentIdx,localColumn,err,error,*999)
                  CALL DomainMapping_LocalToGlobalGet(colsDOFsDomainMapping,localColumn,globalColumn,err,error,*999)
                  CALL List_ItemAdd(columnIndicesLists(localDOFIdx)%ptr,globalColumn,err,error,*999)
                CASE(FIELD_NODE_BASED_INTERPOLATION)
                  NULLIFY(colsDomain)
                  CALL FieldVariable_ComponentDomainGet(colsVariable,componentIdx,colsDomain,err,error,*999)
                  NULLIFY(colsDomainTopology)
                  CALL Domain_DomainTopologyGet(colsDomain,colsDomainTopology,err,error,*999)
                  NULLIFY(colsDomainElements)
                  CALL DomainTopology_DomainElementsGet(colsDomainTopology,colsDomainElements,err,error,*999)
                  NULLIFY(colsBasis)
                  CALL DomainElements_ElementBasisGet(colsDomainElements,elementNumber,colsBasis,err,error,*999)
                  CALL Basis_NumberOfLocalNodesGet(colsBasis,numberOfLocalNodes,err,error,*999)
                  DO localNodeIdx=1,numberOfLocalNodes
                    CALL DomainElements_ElementNodeGet(colsDomainElements,localNodeIdx,elementNumber,nodeNumber2,err,error,*999)
                    CALL Basis_NodeNumberOfDerivativesGet(colsBasis,localNodeIdx,numberOfNodeDerivatives,err,error,*999)
                    DO derivativeIdx=1,numberOfNodeDerivatives
                      CALL DomainElemennts_ElementDerivativeGet(colsDomainElements,derivativeIdx,localNodeIdx,elementNumber,&
                        & derivativeNumber2,err,error,*999)
                      CALL DomainElemennts_ElementVersionGet(colsDomainElements,derivativeIdx,localNodeIdx,elementNumber, &
                        & versionNumber2,err,error,*999)
                      !Find the local and global column and add the global column to the indices list
                      CALL FieldVariable_LocalNodeDOFGet(colsVariable,versionNumber2,derivativeNumber2,nodeNumber2,componentIdx, &
                        & localColumn,err,error,*999)
                      CALL DomainMapping_LocalToGlobalGet(colsDOFsDomainMapping,localColumn,globalColumn,err,error,*999)
                      CALL List_ItemAdd(columnIndicesLists(localDOFIdx)%ptr,globalColumn,err,error,*999)
                    ENDDO !derivativeIdx
                  ENDDO !localNodeIdx
                CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                  CALL FlagError("Grid point based interpolation is not implemented.",err,error,*999)
                CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                  CALL FlagError("Gauss point based interpolation is not implemented.",err,error,*999)
                CASE(FIELD_DATA_POINT_BASED_INTERPOLATION)
                  CALL FlagError("Data point based interpolation is not implemented.",err,error,*999)
                CASE DEFAULT
                  localError="The columns interpolation type of "//TRIM(NumberToVString(interpolationType,"*",err,error))// &
                    & " for component number "//TRIM(NumberToVString(componentIdx,"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              ENDDO !componentIdx
            ENDDO !elementIdx
            CALL List_RemoveDuplicates(columnIndicesLists(localDOFIdx)%ptr,err,error,*999)
            CALL List_NumberOfItemsGet(columnIndicesLists(localDOFIdx)%ptr,numberOfColumns,err,error,*999)
            numberOfNonZeros=numberOfNonZeros+numberOfColumns
            rowIndices(localDOFIdx+1)=numberOfNonZeros+1
          CASE(FIELD_ELEMENT_DOF_TYPE)
            CALL FieldVariable_DOFParameterGetElement(rowsVariable,dofTypeIdx,elementNumber,componentNumber,err,error,*999)
            CALL FieldVariable_ComponentMaxElementInterpParametersGet(rowsVariable,componentNumber,maxElementInterpParameters, &
              & err,error,*999)
            NULLIFY(rowsDomain)
            CALL FieldVariable_ComponentDomainGet(rowsVariable,componentNumber,rowsDomain,err,error,*999)
            NULLIFY(rowsDomainTopology)
            CALL Domain_DomainTopologyGet(rowsDomain,rowsDomainTopology,err,error,*999)
            NULLIFY(rowsDomainElements)
            CALL DomainTopology_DomainElementsGet(rowsDomainTopology,rowsDomainElements,err,error,*999)
            NULLIFY(rowsBasis)
            CALL DomainElements_ElementBasisGet(rowsDomainElements,elementNumber,rowsBasis,err,error,*999)
            !Set up list
            NULLIFY(columnIndicesLists(localDOFIdx)%ptr)
            CALL List_CreateStart(columnIndicesLists(localDOFIdx)%ptr,err,error,*999)
            CALL List_DataTypeSet(columnIndicesLists(localDOFIdx)%ptr,LIST_INTG_TYPE,err,error,*999)
            CALL List_InitialSizeSet(columnIndicesLists(localDOFIdx)%ptr,maxElementInterpParameters+1,err,error,*999) ! size = all nodal dofs + itself
            CALL List_CreateFinish(columnIndicesLists(localDOFIdx)%ptr,err,error,*999)
            DO componentIdx=1,numberOfColumnComponents
              NULLIFY(colsDomain)
              CALL FieldVariable_ComponentDomainGet(colsVariable,componentIdx,colsDomain,err,error,*999)
              NULLIFY(colsDomainTopology)
              CALL Domain_DomainTopologyGet(colsDomain,colsDomainTopology,err,error,*999)
              NULLIFY(colsDomainElements)
              CALL DomainTopology_DomainElementsGet(colsDomainTopology,colsDomainElements,err,error,*999)
              NULLIFY(colsBasis)
              CALL DomainElements_ElementBasisGet(colsDomainElements,elementNumber,colsBasis,err,error,*999)
              CALL FieldVariable_ComponentInterpolationGet(colsVariable,componentIdx,interpolationType,err,error,*999)
              SELECT CASE(interpolationType)
              CASE(FIELD_CONSTANT_INTERPOLATION)
                CALL FlagError("Constant interpolation is not implemented yet.",err,error,*999)
              CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                !It's assumed that element-based variables arne't directly coupled put a diagonal entry
                CALL FieldVariable_LocalElementDOFGet(colsVariable,elementNumber,componentIdx,localColumn,err,error,*999)
                CALL DomainMapping_LocalToGlobalGet(colsDOFsDomainMapping,localColumn,globalColumn,err,error,*999)
                CALL List_ItemAdd(columnIndicesLists(localDOFIdx)%ptr,globalColumn,err,error,*999)
              CASE(FIELD_NODE_BASED_INTERPOLATION)
                !Loop over all nodes in the element (and dofs belonging to them)
                CALL Basis_NumberOfLocalNodesGet(colsBasis,numberOfLocalNodes,err,error,*999)               
                DO localNodeIdx=1,numberOfLocalNodes
                  CALL DomainElements_ElementNodeGet(colsDomainElements,localNodeIdx,elementNumber,nodeNumber2,err,error,*999)
                  CALL Basis_NodeNumberOfDerivativesGet(colsBasis,localNodeIdx,numberOfNodeDerivatives,err,error,*999)
                  DO derivativeIdx=1,numberOfNodeDerivatives
                    CALL DomainElemennts_ElementDerivativeGet(colsDomainElements,derivativeIdx,localNodeIdx,elementNumber, &
                      & derivativeNumber2,err,error,*999)
                    CALL DomainElemennts_ElementVersionGet(colsDomainElements,derivativeIdx,localNodeIdx,elementNumber, &
                      & versionNumber2,err,error,*999)
                    !Find the local and global column and add the global column to the indices list
                    CALL FieldVariable_LocalNodeDOFGet(colsVariable,versionNumber2,derivativeNumber2,nodeNumber2,componentIdx, &
                      & localColumn,err,error,*999)
                    CALL DomainMapping_LocalToGlobalGet(colsDOFsDomainMapping,localColumn,globalColumn,err,error,*999)
                    CALL List_ItemAdd(columnIndicesLists(localDOFIdx)%ptr,globalColumn,err,error,*999)
                  ENDDO !derivativeIdx
                ENDDO !localNodeIdx
              CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                CALL FlagError("Grid point based interpolation is not implemented.",err,error,*999)
              CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                CALL FlagError("Gauss point based interpolation is not implemented.",err,error,*999)
              CASE(FIELD_DATA_POINT_BASED_INTERPOLATION)
                CALL FlagError("Data point based interpolation is not implemented.",err,error,*999)
              CASE DEFAULT
                localError="The columns interpolation type of "//TRIM(NumberToVString(interpolationType,"*",err,error))// &
                  & " for component number "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ENDDO !componentIdx
            !Clean up the list
            CALL List_RemoveDuplicates(columnIndicesLists(localDOFIdx)%ptr,err,error,*999)
            CALL List_NumberOfItemsGet(columnIndicesLists(localDOFIdx)%ptr,numberOfColumns,err,error,*999)
            numberOfNonZeros=numberOfNonZeros+numberOfColumns
            rowIndices(localDOFIdx+1)=numberOfNonZeros+1
          CASE(FIELD_GRID_POINT_DOF_TYPE)
            CALL FlagError("Grid point based DOFs is not implemented yet.",err,error,*999)
          CASE(FIELD_GAUSS_POINT_DOF_TYPE)
            CALL FlagError("Gauss point based DOFs is not implemented yet.",err,error,*999)
          CASE(FIELD_DATA_POINT_DOF_TYPE)
            CALL FlagError("Data point based DOFs is not implemented yet.",err,error,*999)
          CASE DEFAULT
            localError="Local DOF number "//TRIM(NumberToVString(localDOFIdx,"*",err,error))//" has an invalid type."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ENDDO !localDOFIdx
        
      CASE DEFAULT
        localError="The matrix storage type of "//TRIM(NumberToVString(equationsMatrix%storageType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      
    CASE(EQUATIONS_MATRIX_NODAL_STRUCTURE)
      SELECT CASE(equationsMatrix%storageType)
      CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
        !Allocate lists
        ALLOCATE(columnIndicesLists(totalNumberOfRows),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate column indices lists.",err,error,*999)
        !Allocate row indices
        ALLOCATE(rowIndices(totalNumberOfRows+1),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate row indices.",err,error,*999)
        rowIndices(1)=1
        
        !First, loop over the rows and calculate the number of non-zeros
        numberOfNonZeros=0
        DO localDofIdx=1,totalNumberOfRows
          CALL FieldVariable_DOFTypeGet(rowsVariable,localDOFIdx,dofType,dofTypeIdx,err,error,*999)          
          IF(dofType/=FIELD_NODE_DOF_TYPE) THEN
            localError="Local DOF number "//TRIM(NumberToVString(localDOFIdx,"*",err,error))//" is not a node based DOF."
            CALL FlagError(localError,err,error,*999)
          ENDIF
          CALL FieldVariable_DOFParameterGetNode(rowsVariable,dofTypeIdx,versionNumber,derivativeNumber,nodeNumber, &
            & componentNumber,err,error,*999)
          NULLIFY(rowsDomain)
          CALL FieldVariable_ComponentDomainGet(rowsVariable,componentNumber,rowsDomain,err,error,*999)
          CALL FieldVariable_ComponentMaxElementInterpParametersGet(rowsVariable,componentNumber,maxElementInterpParameters, &
            & err,error,*999)
          NULLIFY(rowsDomainTopology)
          CALL Domain_DomainTopologyGet(rowsDomain,rowsDomainTopology,err,error,*999)
          NULLIFY(rowsDomainNodes)
          CALL DomainTopology_DomainNodesGet(rowsDomainTopology,rowsDomainNodes,err,error,*999)
          !Set up list
          NULLIFY(columnIndicesLists(localDofIdx)%ptr)
          CALL List_CreateStart(columnIndicesLists(localDofIdx)%ptr,err,error,*999)
          CALL List_DataTypeSet(columnIndicesLists(localDofIdx)%ptr,LIST_INTG_TYPE,err,error,*999)          
          CALL List_InitialSizeSet(columnIndicesLists(localDofIdx)%ptr,numberOfRowComponents*maxElementInterpParameters, &
            & err,error,*999)
          CALL List_CreateFinish(columnIndicesLists(localDofIdx)%ptr,err,error,*999)
          !Loop over all components, nodes, derivatives and versions
          DO componentIdx=1,numberOfColumnComponents
            NULLIFY(colsDomain)
            CALL FieldVariable_ComponentDomainGet(colsVariable,componentIdx,colsDomain,err,error,*999)
            NULLIFY(colsDomainTopology)
            CALL Domain_DomainTopologyGet(colsDomain,colsDomainTopology,err,error,*999)
            NULLIFY(colsDomainNodes)
            CALL DomainTopology_DomainNodesGet(colsDomainTopology,colsDomainNodes,err,error,*999)
            CALL DomainNodes_NodeNumberOfDerivativesGet(colsDomainNodes,nodeNumber,numberOfDerivatives,err,error,*999)
            DO derivativeIdx=1,numberOfDerivatives
              CALL DomainNodes_DerivativeNumberOfVersionsGet(colsDomainNodes,derivativeIdx,nodeNumber,numberOfVersions, &
                & err,error,*999)              
              DO versionIdx=1,numberOfVersions
                CALL FieldVariable_LocalNodeDOFGet(colsVariable,versionIdx,derivativeIdx,nodeNumber,componentIdx,localColumn, &
                  & err,error,*999)
                CALL DomainMapping_LocalToGlobalGet(colsDOFsDomainMapping,localColumn,globalColumn,err,error,*999)
                CALL List_ItemAdd(columnIndicesLists(localDofIdx)%ptr,globalColumn,err,error,*999)                
              ENDDO !versionIdx
            ENDDO !derivativeIdx
          ENDDO !componentIdx            
          CALL List_RemoveDuplicates(columnIndicesLists(localDofIdx)%ptr,err,error,*999)
          CALL List_NumberOfItemsGet(columnIndicesLists(localDofIdx)%ptr,numberOfColumns,err,error,*999)
          numberOfNonZeros=numberOfNonZeros+numberOfColumns
          rowIndices(localDofIdx+1)=numberOfNonZeros+1
        ENDDO !localDofIdx
        
      CASE DEFAULT
        localError="The matrix storage type of "//TRIM(NumberToVString(equationsMatrix%storageType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      
    CASE(EQUATIONS_MATRIX_DIAGONAL_STRUCTURE)
      CALL FlagError("There is not structure to calculate for a diagonal matrix.",err,error,*998)
    CASE DEFAULT
      localError="The matrix structure type of "//TRIM(NumberToVString(equationsMatrix%structureType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*998)
    END SELECT

    IF(numberOfNonZeros>0) THEN
      !Allocate and setup the column locations
      ALLOCATE(columnIndices(numberOfNonZeros),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate column indices.",err,error,*999)        
      ALLOCATE(list(colsDOFsDomainMapping%numberOfGlobal),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate list.",err,error,*999)
      
      DO localDOFIdx=1,totalNumberOfRows 
        CALL List_DetachAndDestroy(columnIndicesLists(localDOFIdx)%ptr,numberOfColumns,columns,err,error,*999)        
        DO columnIdx=1,numberOfColumns
          !columns stores the list of nonzero column indices for each local row (localDOFIdx)
          columnIndices(rowIndices(localDOFIdx)+columnIdx-1)=columns(columnIdx)             
          !global to local columns
          CALL DomainMapping_LocalNumberFromGlobalGet(colsDOFsDomainMapping,columns(columnIdx),1,localColumn,err,error,*999)
          localDOF=localColumn
          CALL FieldVariable_DOFTypeGet(colsVariable,localDOF,dofType,dofTypeIdx,err,error,*999)
          SELECT CASE(dofType)
          CASE(FIELD_CONSTANT_DOF_TYPE)
            CALL FlagError("Constant based DOFs is not implemented.",err,error,*999)
          CASE(FIELD_ELEMENT_DOF_TYPE)
            CALL FieldVariable_DOFParameterGetElement(rowsVariable,dofTypeIdx,elementNumber,componentNumber,err,error,*999)
            NULLIFY(colsDomain)
            CALL FieldVariable_ComponentDomainGet(colsVariable,componentNumber,colsDomain,err,error,*999)
            NULLIFY(colsDomainTopology)
            CALL Domain_DomainTopologyGet(colsDomain,colsDomainTopology,err,error,*999)
            NULLIFY(colsDomainElements)
            CALL DomainTopology_DomainElementsGet(colsDomainTopology,colsDomainElements,err,error,*999)
            CALL DomainElements_ElementBoundaryElementGet(colsDomainElements,elementNumber,boundaryElement,err,error,*999)
            !Check whether boundary element    
            IF(boundaryElement) CALL LinkedList_Add(list(columns(columnIdx)),localDOFIdx,err,error,*999)
          CASE(FIELD_NODE_DOF_TYPE)
            CALL FieldVariable_DOFParameterGetNode(colsVariable,dofTypeIdx,versionNumber,derivativeNumber,nodeNumber, &
              & componentNumber,err,error,*999)
            NULLIFY(colsDomain)
            CALL FieldVariable_ComponentDomainGet(colsVariable,componentNumber,colsDomain,err,error,*999)
            NULLIFY(colsDomainTopology)
            CALL Domain_DomainTopologyGet(colsDomain,colsDomainTopology,err,error,*999)
            NULLIFY(colsDomainNodes)
            CALL DomainTopology_DomainNodesGt(colsDomainTopology,colsDomainNodes,err,error,*999)
            CALL DomainNodes_NodeBoundaryNodeGet(colsDomainNodes,nodeNumber,boundaryNode,err,error,*999)            
            !Check whether boundary node    
            IF(boundaryNode) CALL LinkedList_Add(list(columns(columnIdx)),localDOFIdx,err,error,*999)
          CASE(FIELD_GRID_POINT_DOF_TYPE)
            CALL FlagError("Grid point based DOFs is not implemented yet.",err,error,*999)
          CASE(FIELD_GAUSS_POINT_DOF_TYPE)
            CALL FlagError("Gauss point based DOFs is not implemented yet.",err,error,*999)
          CASE(FIELD_DATA_POINT_DOF_TYPE)
            CALL FlagError("Data point based DOFs is not implemented yet.",err,error,*999)
          CASE DEFAULT
            localError="Local DOF number "//TRIM(NumberToVString(localDOF,"*",err,error))//" has an invalid type."
            CALL FlagError(localError,err,error,*999)
          END SELECT          
        ENDDO !columnIdx
        DEALLOCATE(columns)                                    
      ENDDO !localDOFIdx
    ENDIF
      
    IF(diagnostics1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Equations matrix structure:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Equations matrix number : ",matrixNumber,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of rows = ",totalNumberOfRows,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of columns = ",numberOfGlobalColumns, err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of non zeros = ",numberOfNonZeros,err,error,*999)
      IF(totalNumberOfRows*numberOfGlobalColumns/=0) THEN
        sparsity=(1.0_DP-REAL(numberOfNonZeros,DP)/REAL(totalNumberOfRows*numberOfGlobalColumns,DP))*100.0_DP
        CALL WriteStringFmtValue(DIAGNOSTIC_OUTPUT_TYPE,"  Sparsity (% of zeros) = ",sparsity,"F6.2",err,error,*999)
      ENDIF
      IF(numberOfNonZeros>0) THEN
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,totalNumberOfRows+1,8,8, &
          & rowIndices,'("  Row indices    :",8(X,I13))','(18X,8(X,I13))',err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfNonZeros,8,8,columnIndices, &
          & '("  Column indices :",8(X,I13))','(18X,8(X,I13))', err,error,*999)
      ENDIF
    ENDIF
      
    EXITS("EquationsMatrix_StructureCalculate")
    RETURN
999 IF(ASSOCIATED(rowIndices)) DEALLOCATE(rowIndices)
    IF(ASSOCIATED(columnIndices)) DEALLOCATE(columnIndices)
    IF(ALLOCATED(columns)) DEALLOCATE(columns)
    IF(ALLOCATED(columnIndicesLists)) THEN
      DO localDofIdx=1,SIZE(columnIndicesLists,1)
        IF(ASSOCIATED(columnIndicesLists(localDofIdx)%ptr)) &
          & CALL List_Destroy(columnIndicesLists(localDofIdx)%ptr,dummyErr,dummyError,*998)
      ENDDO !localDofIdx
      DEALLOCATE(columnIndicesLists)
    ENDIF
998 ERRORSEXITS("EquationsMatrix_StructureCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrix_StructureCalculate

  !
  !================================================================================================================================
  !

  !>Caclulates the matrix structure (sparsity) for a Jacobian matrix.
  SUBROUTINE JacobianMatrix_StructureCalculate(jacobianMatrix,numberOfNonZeros,rowIndices,columnIndices,err,error,*)

    !Argument variables
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix !<A pointer to the Jacobian matrix to calculate the strucute for
    INTEGER(INTG), INTENT(OUT) :: numberOfNonZeros !<On return the number of non-zeros in the matrix
    INTEGER(INTG), POINTER :: rowIndices(:) !<On return a pointer to row location indices in compressed row format. The pointer must be NULL on entry and the calling routine is responsible for deallocation.
    INTEGER(INTG), POINTER :: columnIndices(:) !<On return a pointer to the column location indices in compressed row format. The pointer must be NULL on entry and the calling routine is responsible for deallocation.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: columnIdx,componentIdx,componentNumber,derivativeIdx,derivativeNumber,derivativeNumber2,dofIdx,dofType, &
      & dofTypeIdx,dummyErr,elementIdx,elementNumber,globalColumn,interpolationType,localColumn,localDOFIdx,localNodeIdx, &
      & matrixNumber,maxElementInterpParameters,nodeNumber,nodeNumber2,numberOfColumns,numberOfColumnComponents, &
      & numberOfDerivatives,numberOfGlobalColumns,numberOfLocalNodes,numberOfNodeDerivatives,numberOfRowComponents, &
      & numberOfSurroundingElements,numberOfVersions,residualNumber,totalNumberOfColumns,totalNumberOfRows,versionIdx, &
      & versionNumber,versionNumber2
    INTEGER(INTG), ALLOCATABLE :: columns(:)
    REAL(DP) :: sparsity
    TYPE(BasisType), POINTER :: colsBasis,rowsBasis
    TYPE(DomainType), POINTER :: colsDomain,rowsDomain
    TYPE(DomainElementsType), POINTER :: colsDomainElements,rowsDomainElements
    TYPE(DomainMappingType), POINTER :: colsDOFsDomainMapping,rowsDOFsDomainMapping
    TYPE(DomainNodesType), POINTER :: colsDomainNodes,rowsDomainNodes
    TYPE(DomainTopologyType), POINTER :: colsDomainTopology,rowsDomainTopology
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField
    TYPE(FieldDOFToParamMapType), POINTER :: colsDOFsParamMapping,rowDOFsParamMapping
    TYPE(FieldVariableType), POINTER :: colsVariable,rowsVariable
    TYPE(ListPtrType), ALLOCATABLE :: columnIndicesLists(:)
    TYPE(VARYING_STRING) :: dummyError,localError

    ENTERS("JacobianMatrix_StructureCalculate",err,error,*998)

    numberOfNonZeros=0
    IF(.NOT.ASSOCIATED(jacobianMatrix)) CALL FlagError("Jacobian matrix is not associated.",err,error,*998)
    IF(ASSOCIATED(rowIndices)) CALL FlagError("Row indices is already associated.",err,error,*998)
    IF(ASSOCIATED(columnIndices)) CALL FlagError("Column indices is already associated.",err,error,*998)
    matrixNumber=jacobianMatrix%jacobianNumber
    NULLIFY(residualVector)
    CALL JacobianMatrix_ResidualVectorGet(jacobianMatrix,residualVector,err,error,*999)
    residualNumber=residualVector%residualNumber
    NULLIFY(nonlinearMatrices)
    CALL EquationsMatricesResidual_NonlinearMatricesGet(residualVector,nonlinearMatrices,err,error,*999)
    NULLIFY(vectorMatrices)
    CALL EquationsMatricesNonlinear_VectorMatricesGet(nonlinearMatrices,vectorMatrices,err,error,*998)
    NULLIFY(vectorMapping)
    CALL EquationsMatricesVector_VectorMappingGet(vectorMatrices,vectorMapping,err,error,*998)
    NULLIFY(lhsMapping)
    CALL EquationsMappingVector_LHSMappingGet(vectorMapping,lhsMapping,err,error,*999)
    NULLIFY(rowsVariable)
    CALL EquationsMappingLHS_LHSVariableGet(lhsMapping,rowsVariable,err,error,*999)
    NULLIFY(vectorEquations)
    CALL EquationsMatricesVector_VectorEquationsGet(vectorMatrices,vectorEquations,err,error,*998)
    NULLIFY(equations)
    CALL EquationsVector_EquationsGet(vectorEquations,equations,err,error,*998)
    NULLIFY(equationsSet)
    CALL Equations_EquationsSetGet(equations,equationsSet,err,error,*998)
    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*998)
    NULLIFY(nonlinearMapping)
    CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*998)
    NULLIFY(residualMapping)
    CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,residualNumber,residualMapping,err,error,*999)
    NULLIFY(colsVariable)
    CALL EquationsMappingResidual_JacobianMatrixVariableGet(residualMapping,matrixNumber,colsVariable,err,error,*999)
    CALL FieldVariable_NumberOfComponentsGet(rowsVariable,numberOfRowComponents,err,error,*999)
    CALL FieldVariable_NumberOfComponentsGet(colsVariable,numberOfColumnComponents,err,error,*999)
    NULLIFY(rowsDOFsDomainMapping)
    CALL FieldVariable_DomainMappingGet(rowsVariable,rowsDOFsDomainMapping,err,error,*999)
    CALL DomainMapping_TotalNumberOfLocalGet(rowsDOFsDomainMapping,totalNumberOfRows,err,error,*999)
    NULLIFY(colsDOFsDomainMapping)
    CALL FieldVariable_DomainMappingGet(colsVariable,colsDOFsDomainMapping,err,error,*999)
    CALL DomainMapping_TotalNumberOfLocalGet(colsDOFsDomainMapping,totalNumberOfColumns,err,error,*999)
    CALL DomainMapping_NumberOfGlobalGet(colsDOFsDomainMapping,numberOfGlobalColumns,err,error,*999)
    
    SELECT CASE(jacobianMatrix%structureType)
    CASE(EQUATIONS_MATRIX_NO_STRUCTURE)
      CALL FlagError("Not implemented.",err,error,*998)
    CASE(EQUATIONS_MATRIX_FEM_STRUCTURE)
      SELECT CASE(jacobianMatrix%storageType)
      CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
        !Allocate lists
        ALLOCATE(columnIndicesLists(totalNumberOfRows),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate column indices lists.",err,error,*999)
        !Allocate row indices
        ALLOCATE(rowIndices(totalNumberOfRows+1),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate row indices.",err,error,*999)
        rowIndices(1)=1
        !First, loop over the rows and calculate the number of non-zeros
        numberOfNonZeros=0
        DO localDOFIdx=1,totalNumberOfRows
          CALL FieldVariable_DOFTypeGet(rowsVariable,localDOFIdx,dofType,dofTypeIdx,err,error,*999)
          SELECT CASE(dofType)
          CASE(FIELD_CONSTANT_DOF_TYPE)
            CALL FlagError("Constant based DOFs is not implemented.",err,error,*999)
          CASE(FIELD_NODE_DOF_TYPE)
            CALL FieldVariable_DOFParameterGetNode(rowsVariable,dofTypeIdx,versionNumber,derivativeNumber,nodeNumber, &
              & componentNumber,err,error,*999)
            NULLIFY(rowsDomain)
            CALL FieldVariable_ComponentDomainGet(rowsVariable,componentNumber,rowsDomain,err,error,*999)
            CALL FieldVariable_ComponentMaxElementInterpParametersGet(rowsVariable,componentNumber,maxElementInterpParameters, &
              & err,error,*999)
            NULLIFY(rowsDomainTopology)
            CALL Domain_DomainTopologyGet(rowsDomain,rowsDomainTopology,err,error,*999)
            NULLIFY(rowsDomainNodes)
            CALL DomainTopology_DomainNodesGet(rowsDomainTopology,rowsDomainNodes,err,error,*999)
            CALL DomainNodes_NodeNumberOfSurroundingElementsGet(rowsDomainNodes,nodeNumber,numberOfSurroundingElements, &
              & err,error,*999)
            !Set up list
            NULLIFY(columnIndicesLists(localDOFIdx)%ptr)
            CALL List_CreateStart(columnIndicesLists(localDOFIdx)%ptr,err,error,*999)
            CALL List_DataTypeSet(columnIndicesLists(localDOFIdx)%ptr,LIST_INTG_TYPE,err,error,*999)
            CALL List_InitialSizeSet(columnIndicesLists(localDOFIdx)%ptr,numberOfSurroundingElements*maxElementInterpParameters, &
              & err,error,*999)
            CALL List_CreateFinish(columnIndicesLists(localDOFIdx)%ptr,err,error,*999)
            !Loop over all elements containing the dof
            DO elementIdx=1,numberOfSurroundingElements
              CALL DomainNodes_NodeSurroundingElementGet(rowsDomainNodes,elementIdx,nodeNumber,elementNumber,err,error,*999)
              DO componentIdx=1,numberOfColumnComponents
                CALL FieldVariable_ComponentInterpolationGet(colsVariable,componentIdx,interpolationType,err,error,*999)
                SELECT CASE(interpolationType)
                CASE(FIELD_CONSTANT_INTERPOLATION)
                  ! do nothing? this will probably never be encountered...?
                  CALL FlagError("Not implemented?",err,error,*999)
                CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                  CALL FieldVariable_LocalElementDOFGet(colsVariable,elementNumber,componentIdx,localColumn,err,error,*999)
                  CALL DomainMapping_LocalToGlobalGet(colsDOFsDomainMapping,localColumn,globalColumn,err,error,*999)
                  CALL List_ItemAdd(columnIndicesLists(localDOFIdx)%ptr,globalColumn,err,error,*999)
                CASE(FIELD_NODE_BASED_INTERPOLATION)
                  NULLIFY(colsDomain)
                  CALL FieldVariable_ComponentDomainGet(colsVariable,componentIdx,colsDomain,err,error,*999)
                  NULLIFY(colsDomainTopology)
                  CALL Domain_DomainTopologyGet(colsDomain,colsDomainTopology,err,error,*999)
                  NULLIFY(colsDomainElements)
                  CALL DomainTopology_DomainElementsGet(colsDomainTopology,colsDomainElements,err,error,*999)
                  NULLIFY(colsBasis)
                  CALL DomainElements_ElementBasisGet(colsDomainElements,elementNumber,colsBasis,err,error,*999)
                  CALL Basis_NumberOfLocalNodesGet(colsBasis,numberOfLocalNodes,err,error,*999)
                  DO localNodeIdx=1,numberOfLocalNodes
                    CALL DomainElements_ElementNodeGet(colsDomainElements,localNodeIdx,elementNumber,nodeNumber2, &
                      & err,error,*999)
                    CALL Basis_NodeNumberOfDerivativesGet(colsBasis,localNodeIdx,numberOfNodeDerivatives,err,error,*999)
                    DO derivativeIdx=1,numberOfNodeDerivatives
                      CALL DomainElemennts_ElementDerivativeGet(colsDomainElements,derivativeIdx,localNodeIdx,elementNumber,&
                        & derivativeNumber2,err,error,*999)
                      CALL DomainElemennts_ElementVersionGet(colsDomainElements,derivativeIdx,localNodeIdx,elementNumber, &
                        & versionNumber2,err,error,*999)
                      !Find the local and global column and add the global column to the indices list
                      CALL FieldVariable_LocalNodeDOFGet(colsVariable,versionNumber2,derivativeNumber2,nodeNumber2,componentIdx, &
                        & localColumn,err,error,*999)
                      CALL DomainMapping_LocalToGlobalGet(colsDOFsDomainMapping,localColumn,globalColumn,err,error,*999)
                      CALL List_ItemAdd(columnIndicesLists(localDOFIdx)%ptr,globalColumn,err,error,*999)
                    ENDDO !derivative
                  ENDDO !localNodeIdx
                CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                  CALL FlagError("Grid point based interpolation is not implemented.",err,error,*999)
                CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                  CALL FlagError("Gauss point based interpolation is not implemented.",err,error,*999)
                CASE(FIELD_DATA_POINT_BASED_INTERPOLATION)
                  CALL FlagError("Data point based interpolation is not implemented.",err,error,*999)
                CASE DEFAULT
                  localError="The columns interpolation type of "//TRIM(NumberToVString(interpolationType,"*",err,error))// &
                    & " for component number "//TRIM(NumberToVString(componentIdx,"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              ENDDO !componentIdx
            ENDDO !elementIdx
            CALL List_RemoveDuplicates(columnIndicesLists(localDOFIdx)%ptr,err,error,*999)
            CALL List_NumberOfItemsGet(columnIndicesLists(localDOFIdx)%ptr,numberOfColumns,err,error,*999)
            numberOfNonZeros=numberOfNonZeros+numberOfColumns
            rowIndices(localDOFIdx+1)=numberOfNonZeros+1
          CASE(FIELD_ELEMENT_DOF_TYPE)
            CALL FieldVariable_DOFParameterGetElement(rowsVariable,dofTypeIdx,elementNumber,componentNumber,err,error,*999)
            CALL FieldVariable_ComponentMaxElementInterpParametersGet(rowsVariable,componentNumber,maxElementInterpParameters, &
              & err,error,*999)
            NULLIFY(rowsDomain)
            CALL FieldVariable_ComponentDomainGet(rowsVariable,componentNumber,rowsDomain,err,error,*999)
            NULLIFY(rowsDomainTopology)
            CALL Domain_DomainTopologyGet(rowsDomain,rowsDomainTopology,err,error,*999)
            NULLIFY(rowsDomainElements)
            CALL DomainTopology_DomainElementsGet(rowsDomainTopology,rowsDomainElements,err,error,*999)
            NULLIFY(rowsBasis)
            CALL DomainElements_ElementBasisGet(rowsDomainElements,elementNumber,rowsBasis,err,error,*999)
            !Set up list
            NULLIFY(columnIndicesLists(localDOFIdx)%ptr)
            CALL List_CreateStart(columnIndicesLists(localDOFIdx)%ptr,err,error,*999)
            CALL List_DataTypeSet(columnIndicesLists(localDOFIdx)%ptr,LIST_INTG_TYPE,err,error,*999)
            CALL List_InitialSizeSet(columnIndicesLists(localDOFIdx)%ptr,maxElementInterpParameters+1,err,error,*999) ! size = all nodal dofs + itself
            CALL List_CreateFinish(columnIndicesLists(localDOFIdx)%ptr,err,error,*999)
            DO componentIdx=1,numberOfColumnComponents
              NULLIFY(colsDomain)
              CALL FieldVariable_ComponentDomainGet(colsVariable,componentIdx,colsDomain,err,error,*999)
              NULLIFY(colsDomainTopology)
              CALL Domain_DomainTopologyGet(colsDomain,colsDomainTopology,err,error,*999)
              NULLIFY(colsDomainElements)
              CALL DomainTopology_DomainElementsGet(colsDomainTopology,colsDomainElements,err,error,*999)
              NULLIFY(colsBasis)
              CALL DomainElements_ElementBasisGet(colsDomainElements,elementNumber,colsBasis,err,error,*999)
              CALL FieldVariable_ComponentInterpolationGet(colsVariable,componentIdx,interpolationType,err,error,*999)
              SELECT CASE(interpolationType)
              CASE(FIELD_CONSTANT_INTERPOLATION)
                CALL FlagError("Constant interpolation is not implemented.",err,error,*999)
              CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                !It's assumed that element-based variables arne't directly coupled put a diagonal entry
                CALL FieldVariable_LocalElementDOFGet(colsVariable,elementNumber,componentIdx,localColumn,err,error,*999)
                CALL DomainMapping_LocalToGlobalGet(colsDOFsDomainMapping,localColumn,globalColumn,err,error,*999)
                CALL List_ItemAdd(columnIndicesLists(localDOFIdx)%ptr,globalColumn,err,error,*999)
              CASE(FIELD_NODE_BASED_INTERPOLATION)
                !Loop over all nodes in the element (and dofs belonging to them)
                CALL Basis_NumberOfLocalNodesGet(colsBasis,numberOfLocalNodes,err,error,*999)               
                DO localNodeIdx=1,numberOfLocalNodes
                  CALL DomainElements_ElementNodeGet(colsDomainElements,localNodeIdx,elementNumber,nodeNumber2,err,error,*999)
                  CALL Basis_NodeNumberOfDerivativesGet(colsBasis,localNodeIdx,numberOfNodeDerivatives,err,error,*999)
                  DO derivativeIdx=1,numberOfNodeDerivatives
                    CALL DomainElemennts_ElementDerivativeGet(colsDomainElements,derivativeIdx,localNodeIdx,elementNumber, &
                      & derivativeNumber2,err,error,*999)
                    CALL DomainElemennts_ElementVersionGet(colsDomainElements,derivativeIdx,localNodeIdx,elementNumber, &
                      & versionNumber2,err,error,*999)
                    !Find the local and global column and add the global column to the indices list
                    CALL FieldVariable_LocalNodeDOFGet(colsVariable,versionNumber2,derivativeNumber2,nodeNumber2,componentIdx, &
                      & localColumn,err,error,*999)
                    CALL DomainMapping_LocalToGlobalGet(colsDOFsDomainMapping,localColumn,globalColumn,err,error,*999)
                    CALL List_ItemAdd(columnIndicesLists(localDOFIdx)%ptr,globalColumn,err,error,*999)
                  ENDDO !derivativeIdx
                ENDDO !localNodeIdx
              CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                CALL FlagError("Grid point based interpolation is not implemented.",err,error,*999)
              CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                CALL FlagError("Gauss point based interpolation is not implemented.",err,error,*999)
              CASE(FIELD_DATA_POINT_BASED_INTERPOLATION)
                CALL FlagError("Data point based interpolation is not implemented.",err,error,*999)
              CASE DEFAULT
                localError="The columns interpolation type of "//TRIM(NumberToVString(interpolationType,"*",err,error))// &
                  & " for component number "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ENDDO !componentIdx
            !Clean up the list
            CALL List_RemoveDuplicates(columnIndicesLists(localDOFIdx)%ptr,err,error,*999)
            CALL List_NumberOfItemsGet(columnIndicesLists(localDOFIdx)%ptr,numberOfColumns,err,error,*999)
            numberOfNonZeros=numberOfNonZeros+numberOfColumns
            rowIndices(localDOFIdx+1)=numberOfNonZeros+1
          CASE(FIELD_GRID_POINT_DOF_TYPE)
            CALL FlagError("Grid point based DOFs is not implemented yet.",err,error,*999)
          CASE(FIELD_GAUSS_POINT_DOF_TYPE)
            CALL FlagError("Gauss point based DOFs is not implemented yet.",err,error,*999)
          CASE(FIELD_DATA_POINT_DOF_TYPE)
            CALL FlagError("Data point based DOFs is not implemented yet.",err,error,*999)
         CASE DEFAULT
            localError="Local dof number "//TRIM(NumberToVString(localDOFIdx,"*",err,error))//" has an invalid type."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ENDDO !localDOFIdx
      CASE DEFAULT
        localError="The matrix storage type of "//TRIM(NumberToVString(jacobianMatrix%storageType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_MATRIX_NODAL_STRUCTURE)
      SELECT CASE(jacobianMatrix%storageType)
      CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
        !Allocate lists
        ALLOCATE(columnIndicesLists(totalNumberOfRows),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate column indices lists.",err,error,*999)
        !Allocate row indices
        ALLOCATE(rowIndices(totalNumberOfRows+1),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate row indices.",err,error,*999)
        rowIndices(1)=1
        !First, loop over the rows and calculate the number of non-zeros
        numberOfNonZeros=0
        DO localDofIdx=1,totalNumberOfRows
          CALL FieldVariable_DOFTypeGet(rowsVariable,localDOFIdx,dofType,dofTypeIdx,err,error,*999)          
          SELECT CASE(dofType)
          CASE(FIELD_CONSTANT_DOF_TYPE)
            CALL FlagError("Constant based DOFs is not implemented yet.",err,error,*999)
          CASE(FIELD_NODE_DOF_TYPE)
            CALL FieldVariable_DOFParameterGetNode(rowsVariable,dofTypeIdx,versionNumber,derivativeNumber,nodeNumber, &
              & componentNumber,err,error,*999)
            NULLIFY(rowsDomain)
            CALL FieldVariable_ComponentDomainGet(rowsVariable,componentNumber,rowsDomain,err,error,*999)
            CALL FieldVariable_ComponentMaxElementInterpParametersGet(rowsVariable,componentNumber,maxElementInterpParameters, &
              & err,error,*999)
            NULLIFY(rowsDomainTopology)
            CALL Domain_DomainTopologyGet(rowsDomain,rowsDomainTopology,err,error,*999)
            NULLIFY(rowsDomainNodes)
            CALL DomainTopology_DomainNodesGet(rowsDomainTopology,rowsDomainNodes,err,error,*999)
            !Set up list
            NULLIFY(columnIndicesLists(localDofIdx)%ptr)
            CALL List_CreateStart(columnIndicesLists(localDofIdx)%ptr,err,error,*999)
            CALL List_DataTypeSet(columnIndicesLists(localDofIdx)%ptr,LIST_INTG_TYPE,err,error,*999)            
            CALL List_InitialSizeSet(columnIndicesLists(localDofIdx)%ptr,numberOfRowComponents*maxElementInterpParameters, &
              & err,error,*999)
            CALL List_CreateFinish(columnIndicesLists(localDofIdx)%ptr,err,error,*999)
            !Loop over all components,nodes,derivatives, and versions
            DO componentIdx=1,numberOfColumnComponents
              CALL FieldVariable_ComponentInterpolationGet(colsVariable,componentIdx,interpolationType,err,error,*999)
              SELECT CASE(interpolationType)
              CASE(FIELD_NODE_BASED_INTERPOLATION)
                NULLIFY(colsDomain)
                CALL FieldVariable_ComponentDomainGet(colsVariable,componentIdx,colsDomain,err,error,*999)
                NULLIFY(colsDomainTopology)
                CALL Domain_DomainTopologyGet(colsDomain,colsDomainTopology,err,error,*999)
                NULLIFY(colsDomainNodes)
                CALL DomainTopology_DomainNodesGet(colsDomainTopology,colsDomainNodes,err,error,*999)
                CALL DomainNodes_NodeNumberOfDerivativesGet(colsDomainNodes,nodeNumber,numberOfDerivatives,err,error,*999)
                DO derivativeIdx=1,numberOfDerivatives
                  CALL DomainNodes_DerivativeNumberOfVersionsGet(colsDomainNodes,derivativeIdx,nodeNumber,numberOfVersions, &
                    & err,error,*999)              
                  DO versionIdx=1,numberOfVersions
                    CALL FieldVariable_LocalNodeDOFGet(colsVariable,versionIdx,derivativeIdx,nodeNumber,componentIdx,localColumn, &
                      & err,error,*999)
                    CALL DomainMapping_LocalToGlobalGet(colsDOFsDomainMapping,localColumn,globalColumn,err,error,*999)
                    CALL List_ItemAdd(columnIndicesLists(localDofIdx)%ptr,globalColumn,err,error,*999)                    
                  ENDDO !versionIdx
                ENDDO !derivativeIdx
              CASE DEFAULT
                localError="Local dof number "//TRIM(NumberToVString(localDofIdx,"*",err,error))//" has invalid interpolation type."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ENDDO !componentIdx            
            CALL List_RemoveDuplicates(columnIndicesLists(localDofIdx)%ptr,err,error,*999)
            CALL List_NumberOfItemsGet(columnIndicesLists(localDofIdx)%ptr,numberOfColumns,err,error,*999)
            numberOfNonZeros=numberOfNonZeros+numberOfColumns
            rowIndices(localDofIdx+1)=numberOfNonZeros+1
          CASE(FIELD_ELEMENT_DOF_TYPE)
            CALL FlagError("Element based DOFs is not implemented.",err,error,*999)
          CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
            CALL FlagError("Grid point based DOFs is not implemented.",err,error,*999)
          CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
            CALL FlagError("Gauss point based DOFs is not implemented.",err,error,*999)
          CASE DEFAULT
            localError="Local dof number "//TRIM(NumberToVString(localDofIdx,"*",err,error))//" has an invalid type."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ENDDO !localDofIdx
      CASE DEFAULT
        localError="The matrix storage type of "//TRIM(NumberToVString(jacobianMatrix%storageType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The matrix structure type of "//TRIM(NumberToVString(jacobianMatrix%structureType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*998)
    END SELECT

    IF(numberOfNonZeros>0) THEN
      !Allocate and setup the column locations
      ALLOCATE(columnIndices(numberOfNonZeros),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate column indices.",err,error,*999)
      DO localDOFIdx=1,rowsDOFsDomainMapping%totalNumberOfLocal
        CALL List_DetachAndDestroy(columnIndicesLists(localDOFIdx)%ptr,numberOfColumns,columns,err,error,*999)
        DO columnIdx=1,numberOfColumns
          columnIndices(rowIndices(localDOFIdx)+columnIdx-1)=columns(columnIdx)
        ENDDO !columnIdx
        DEALLOCATE(columns)
      ENDDO !localDOFIdx
    ENDIF
    
    IF(diagnostics1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Jacobian matrix structure:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of rows = ",totalNumberOfRows,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of columns = ",numberOfGlobalColumns,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of non zeros = ",numberOfNonZeros,err,error,*999)
      IF(totalNumberOfRows*numberOfGlobalColumns/=0) THEN
        sparsity=(1.0_DP-REAL(numberOfNonZeros,DP)/REAL(totalNumberOfRows*numberOfGlobalColumns,DP))*100.0_DP
        CALL WriteStringFmtValue(DIAGNOSTIC_OUTPUT_TYPE,"  Sparsity (% of zeros) = ",sparsity,"F6.2",err,error,*999)
      ENDIF
      IF(numberOfNonZeros>0) THEN
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,totalNumberOfRows+1,8,8, &
          & rowIndices,'("  Row indices    :",8(X,I13))','(18X,8(X,I13))',err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfNonZeros,8,8,columnIndices,&
          & '("  Column indices :",8(X,I13))','(18X,8(X,I13))', err,error,*999)
      ENDIF
    ENDIF
    
    EXITS("JacobianMatrix_StructureCalculate")
    RETURN
999 IF(ASSOCIATED(rowIndices)) DEALLOCATE(rowIndices)
    IF(ASSOCIATED(columnIndices)) DEALLOCATE(columnIndices)
    IF(ALLOCATED(columns)) DEALLOCATE(columns)
    IF(ALLOCATED(columnIndicesLists)) THEN
      DO localDofIdx=1,colsDOFsDomainMapping%totalNumberOfLocal
        IF(ASSOCIATED(columnIndicesLists(localDofIdx)%ptr)) &
          & CALL List_Destroy(columnIndicesLists(localDofIdx)%ptr,dummyErr,dummyError,*998)
      ENDDO !localDofIdx
      DEALLOCATE(columnIndicesLists)
    ENDIF
998 ERRORSEXITS("JacobianMatrix_StructureCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE JacobianMatrix_StructureCalculate
  
  !
  !================================================================================================================================
  !
 
END MODULE EquationsMatricesRoutines
