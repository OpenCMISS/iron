!> \file
!> \author Chris Bradley
!> \brief This module contains all interface matrices routines.
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

!>This module contains all interface matrices routines.
MODULE InterfaceMatricesRoutines

  USE BaseRoutines
  USE BasisAccessRoutines
  USE DecompositionAccessRoutines
  USE DistributedMatrixVector
  USE DistributedMatrixVectorAccessRoutines
  USE DomainMappings
  USE EquationsMatricesRoutines
  USE FieldRoutines
  USE FieldAccessRoutines
  USE InputOutput
  USE InterfaceAccessRoutines
  USE InterfaceConditionAccessRoutines
  USE InterfaceEquationsAccessRoutines
  USE InterfaceMappingAccessRoutines
  USE InterfaceMatricesAccessRoutines
  USE ISO_VARYING_STRING
  USE Kinds
  USE MatrixVector
  USE Strings
  USE Types

#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  INTERFACE InterfaceMatrices_MatrixCoefficientSet
    MODULE PROCEDURE InterfaceMatrices_MatrixCoefficientSet0
    MODULE PROCEDURE InterfaceMatrices_MatrixCoefficientSet1
  END INTERFACE InterfaceMatrices_MatrixCoefficientSet

  INTERFACE InterfaceMatrices_StorageTypeSet
    MODULE PROCEDURE InterfaceMatrices_StorageTypeSet0
    MODULE PROCEDURE InterfaceMatrices_StorageTypeSet1    
  END INTERFACE InterfaceMatrices_StorageTypeSet

  INTERFACE InterfaceMatrices_StructureTypeSet
    MODULE PROCEDURE InterfaceMatrices_StructureTypeSet0
    MODULE PROCEDURE InterfaceMatrices_StructureTypeSet1    
  END INTERFACE InterfaceMatrices_StructureTypeSet

  INTERFACE InterfaceMatrices_TransposeMatrixCoefficientSet
    MODULE PROCEDURE InterfaceMatrices_TransposeMatrixCoefficientSet0
    MODULE PROCEDURE InterfaceMatrices_TransposeMatrixCoefficientSet1
  END INTERFACE InterfaceMatrices_TransposeMatrixCoefficientSet

  PUBLIC InterfaceMatrices_CreateFinish,InterfaceMatrices_CreateStart

  PUBLIC InterfaceMatrices_Destroy

  PUBLIC InterfaceMatrices_ElementAdd

  PUBLIC InterfaceMatrices_ElementCalculate

  PUBLIC InterfaceMatrices_ElementFinalise,InterfaceMatrices_ElementInitialise

  PUBLIC InterfaceMatrices_Output

  PUBLIC InterfaceMatrices_MatrixCoefficientSet

  PUBLIC InterfaceMatrices_StorageTypeSet

  PUBLIC InterfaceMatrices_StructureTypeSet

  PUBLIC InterfaceMatrices_TransposeMatrixCoefficientSet

  PUBLIC InterfaceMatrices_ValueInitialise

  PUBLIC InterfaceMatrix_MatrixCoefficientSet

  PUBLIC InterfaceMatrix_StorageTypeSet

  PUBLIC InterfaceMatrix_StructureTypeSet

  PUBLIC InterfaceMatrix_TransposeMatrixCoefficientSet

CONTAINS

  !
  !================================================================================================================================
  !

  !>Adds the element matrices into the interface matrices.
  SUBROUTINE InterfaceMatrices_ElementAdd(interfaceMatrices,err,error,*)

    !Argument variables
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices !<A pointer to the interface matrices
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(InterfaceMatrixType), POINTER :: interfaceMatrix
    TYPE(InterfaceRHSType), POINTER :: rhsVector

    ENTERS("InterfaceMatrices_ElementAdd",err,error,*999)

    IF(.NOT.ASSOCIATED(interfaceMatrices)) CALL FlagError("Interface matrices is not allocated.",err,error,*999)
    
    !Add the element matrices
    DO matrixIdx=1,interfaceMatrices%numberOfInterfaceMatrices
      NULLIFY(interfaceMatrix)
      CALL InterfaceMatrices_InterfaceMatrixGet(interfaceMatrices,matrixIdx,interfaceMatrix,err,error,*999)
      IF(interfaceMatrix%updateMatrix) THEN
        !Add the element matrix into the distributed interface equations matrix
        CALL DistributedMatrix_ValuesAdd(interfaceMatrix%matrix,interfaceMatrix%elementMatrix%rowDOFS(1: &
          & interfaceMatrix%elementMatrix%numberOfRows),interfaceMatrix%elementMatrix%columnDOFS(1: &
          & interfaceMatrix%elementMatrix%numberOfColumns),interfaceMatrix%elementMatrix%matrix(1: &
          & interfaceMatrix%elementMatrix%numberOfRows,1:interfaceMatrix%elementMatrix%numberOfColumns), &
          & err,error,*999)
        !If the interface matrix has a transpose add it
        IF(interfaceMatrix%hasTranspose) THEN
          CALL DistributedMatrix_ValuesAdd(interfaceMatrix%matrixTranspose,interfaceMatrix%elementMatrix%columnDOFS(1: &
            & interfaceMatrix%elementMatrix%numberOfColumns),interfaceMatrix%elementMatrix%rowDOFS(1: &
            & interfaceMatrix%elementMatrix%numberOfRows),TRANSPOSE(interfaceMatrix%elementMatrix%matrix(1: &
            & interfaceMatrix%elementMatrix%numberOfRows,1:interfaceMatrix%elementMatrix%numberOfColumns)), &
            & err,error,*999)
        ENDIF
      ENDIF
    ENDDO !matrixIdx
    NULLIFY(rhsVector)
    CALL InterfaceMatrices_RHSVectorExists(interfaceMatrices,rhsVector,err,error,*999)
    IF(ASSOCIATED(rhsVector)) THEN
      IF(rhsVector%updateVector) THEN
        !Add the rhs element vector
        CALL DistributedVector_ValuesAdd(rhsVector%rhsVector,rhsVector%elementVector%rowDOFS(1: &
          & rhsVector%elementVector%numberOfRows),rhsVector%elementVector%vector(1:rhsVector% &
          & elementVector%numberOfRows),err,error,*999)
      ENDIF
    ENDIF
    
    EXITS("InterfaceMatrices_ElementAdd")
    RETURN
999 ERRORSEXITS("InterfaceMatrices_ElementAdd",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMatrices_ElementAdd

  !
  !================================================================================================================================
  !

  !>Calculate the positions of the element matrices in the interface matrices. 
  SUBROUTINE InterfaceMatrices_ElementCalculate(interfaceMatrices,interfaceElementNumber,err,error,*)

    !Argument variables
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices !<A pointer to the interface matrices
    INTEGER(INTG), INTENT(IN) :: interfaceElementNumber !<The element number to calculate the mappings for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: interfaceConditionMethod,integrationType,matrixIdx,rowsElementNumber,rowsMeshIdx
    TYPE(FieldVariableType), POINTER :: colsFieldVariable,rowsFieldVariable
    TYPE(InterfaceType), POINTER :: interface
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping
    TYPE(InterfaceMappingRHSType), POINTER :: rhsMapping
    TYPE(InterfaceMatrixType), POINTER :: interfaceMatrix
    TYPE(InterfaceMatrixToVarMapType), POINTER :: interfaceMatrixToVarMap
    TYPE(InterfaceMeshConnectivityType), POINTER :: meshConnectivity
    TYPE(InterfacePointsConnectivityType), POINTER :: pointsConnectivity
    TYPE(InterfaceRHSType), POINTER :: rhsVector
    TYPE(VARYING_STRING) :: localError

    ENTERS("InterfaceMatrices_ElementCalculate",err,error,*999)

    NULLIFY(interfaceMapping)
    CALL InterfaceMatrices_InterfaceMappingGet(interfaceMatrices,interfaceMapping,err,error,*999)
    NULLIFY(interfaceEquations)
    CALL InterfaceMatrices_InterfaceEquationsGet(interfaceMatrices,interfaceEquations,err,error,*999)
    NULLIFY(interfaceCondition)
    CALL InterfaceEquations_InterfaceConditionGet(interfaceEquations,interfaceCondition,err,error,*999)
    NULLIFY(INTERFACE)
    CALL InterfaceCondition_InterfaceGet(interfaceCondition,INTERFACE,err,error,*999)
    CALL InterfaceCondition_MethodGet(interfaceCondition,interfaceConditionMethod,err,error,*999)
    CALL InterfaceCondition_IntegrationTypeGet(interfaceCondition,integrationType,err,error,*999)
    
    SELECT CASE(integrationType)
    CASE(INTERFACE_CONDITION_GAUSS_INTEGRATION)
      NULLIFY(meshConnectivity)
      CALL Interface_MeshConnectivityGet(INTERFACE,meshConnectivity,err,error,*999)
      NULLIFY(colsFieldVariable)
!!\todo: TEMPORARY: Needs generalising
      CALL InterfaceMapping_LagrangeVariableGet(interfaceMapping,colsFieldVariable,err,error,*999)       
      !Calculate the row and columns for the interface equations matrices
      DO matrixIdx=1,interfaceMatrices%numberOfInterfaceMatrices
        NULLIFY(interfaceMatrix)
        CALL InterfaceMatrices_InterfaceMatrixGet(interfaceMatrices,matrixIdx,interfaceMatrix,err,error,*999)
        NULLIFY(interfaceMatrixToVarMap)
        CALL InterfaceMapping_InterfaceMatrixToVarMapGet(interfaceMapping,matrixIdx,interfaceMatrixToVarMap,err,error,*999)
        NULLIFY(rowsFieldVariable)
        CALL InterfaceMappingIMToVMap_VariableGet(interfaceMatrixToVarMap,rowsFieldVariable,err,error,*999)
        CALL InterfaceMappingIMToVMap_MeshIndexGet(interfaceMatrixToVarMap,rowsMeshIdx,err,error,*999)
        IF(ASSOCIATED(rowsFieldVariable,colsFieldVariable)) THEN
          !If the rows and column variables are both the Lagrange variable (this is the diagonal matrix)
          rowsElementNumber=interfaceElementNumber
        ELSE
          CALL InterfaceMeshConnectivity_CoupledElementNumberGet(meshConnectivity,interfaceElementNumber,rowsMeshIdx, &
            & rowsElementNumber,err,error,*999)
         ENDIF
        CALL EquationsMatrices_ElementMatrixCalculate(interfaceMatrix%elementMatrix, &
          & interfaceMatrix%updateMatrix,[rowsElementNumber],[interfaceElementNumber],rowsFieldVariable, &
          & colsFieldVariable,err,error,*999)
      ENDDO !matrixIdx
    CASE(INTERFACE_CONDITION_DATA_POINTS_INTEGRATION)
      NULLIFY(pointsConnectivity)
      CALL Interface_PointsConnectivityGet(INTERFACE,pointsConnectivity,err,error,*999)
      DO matrixIdx=1,interfaceMatrices%numberOfInterfaceMatrices
        NULLIFY(interfaceMatrix)
        CALL InterfaceMatrices_InterfaceMatrixGet(interfaceMatrices,matrixIdx,interfaceMatrix,err,error,*999)
        IF(interfaceConditionmethod==INTERFACE_CONDITION_PENALTY_METHOD.AND. &
          matrixIdx==interfaceMatrices%numberOfInterfaceMatrices) THEN
          NULLIFY(rowsFieldVariable)
          CALL InterfaceMapping_LagrangeVariableGet(interfaceMapping,rowsFieldVariable,err,error,*999)
          NULLIFY(colsFieldVariable)
          CALL InterfaceMapping_LagrangeVariableGet(interfaceMapping,colsFieldVariable,err,error,*999)
          CALL EquationsMatrices_ElementMatrixCalculate(interfaceMatrix%elementMatrix, &
            & interfaceMatrix%updateMatrix,[interfaceElementNumber],[interfaceElementNumber], &
            & rowsFieldVariable,colsFieldVariable,err,error,*999)
        ELSE
          NULLIFY(interfaceMatrixToVarMap)
          CALL InterfaceMapping_InterfaceMatrixToVarMapGet(interfaceMapping,matrixIdx,interfaceMatrixToVarMap,err,error,*999)
          NULLIFY(rowsFieldVariable)
          CALL InterfaceMappingIMToVMap_VariableGet(interfaceMatrixToVarMap,rowsFieldVariable,err,error,*999)
          CALL InterfaceMappingIMToVMap_MeshIndexGet(interfaceMatrixToVarMap,rowsMeshIdx,err,error,*999)
          NULLIFY(colsFieldVariable)
          CALL InterfaceMapping_LagrangeVariableGet(interfaceMapping,colsFieldVariable,err,error,*999)
          CALL EquationsMatrices_ElementMatrixCalculate(interfaceMatrix%elementMatrix, &
            & interfaceMatrix%updateMatrix,pointsConnectivity%coupledElements(InterfaceElementNumber, &
            & rowsMeshIdx)%elementNumbers,[interfaceElementNumber],rowsFieldVariable,colsFieldVariable, &
            & err,error,*999)
        ENDIF
      ENDDO !matrixIdx
    CASE DEFAULT
      localError="The interface condition integration type of "// &
        & TRIM(NumberToVString(integrationType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    !RHS element matrix dofs are the same for both mesh and points connectivity, right now
    NULLIFY(rhsVector)
    CALL InterfaceMatrices_RHSVectorGet(interfaceMatrices,rhsVector,err,error,*999)
    NULLIFY(rhsMapping)
    CALL InterfaceMapping_RHSMappingGet(interfaceMapping,rhsMapping,err,error,*999)
    !Calculate the rows for the equations RHS
    NULLIFY(rowsFieldVariable)
    CALL InterfaceMappingRHS_RHSVariableGet(rhsMapping,rowsFieldVariable,err,error,*999)
    CALL EquationsMatrices_ElementVectorCalculate(rhsVector%elementVector,rhsVector%updateVector, &
      & interfaceElementNumber,rowsFieldVariable,err,error,*999)
    
    EXITS("InterfaceMatrices_ElementCalculate")
    RETURN
999 ERRORSEXITS("InterfaceMatrices_ElementCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMatrices_ElementCalculate

  !
  !================================================================================================================================
  !

  !>Finalise the element calculation information for interface matrices and deallocate all memory
  SUBROUTINE InterfaceMatrices_ElementFinalise(interfaceMatrices,err,error,*)

    !Argument variables
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices !<The interface matrices for which to finalise the elements
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(InterfaceMatrixType), POINTER :: interfaceMatrix
    TYPE(InterfaceRHSType), POINTER :: rhsVector
    
    ENTERS("InterfaceMatrices_ElementFinalise",err,error,*999)

    IF(.NOT.ASSOCIATED(interfaceMatrices)) CALL FlagError("Interface matrices is not associated.",err,error,*999)
    DO matrixIdx=1,interfaceMatrices%numberOfInterfaceMatrices
      NULLIFY(interfaceMatrix)
      CALL InterfaceMatrices_InterfaceMatrixGet(interfaceMatrices,matrixIdx,interfaceMatrix,err,error,*999)
      CALL EquationsMatrices_ElementMatrixFinalise(interfaceMatrix%elementMatrix,err,error,*999)
    ENDDO !matrixIdx
    rhsVector=>interfaceMatrices%rhsVector
    IF(ASSOCIATED(rhsVector)) THEN
      !Finalise the interface element vector
      rhsVector%elementVector%maxNumberOfRows=0
      IF(ALLOCATED(rhsVector%elementVector%rowDOFS)) DEALLOCATE(rhsVector%elementVector%rowDOFS)
      IF(ALLOCATED(rhsVector%elementVector%vector)) DEALLOCATE(rhsVector%elementVector%vector)
    ENDIF
    
    EXITS("InterfaceMatrices_ElementFinalise")
    RETURN
999 ERRORSEXITS("InterfaceMatrices_ElementFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMatrices_ElementFinalise

  !
  !================================================================================================================================
  !

  !>Initialise the element calculation information for the interface matrices
  SUBROUTINE InterfaceMatrices_ElementInitialise(interfaceMatrices,err,error,*)

    !Argument variables
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices !The interface matrices to initialise the element information for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: colsNumberOfElements,integrationType,interfaceConditionMethod,matrixIdx,rowsMeshIdx,rowsNumberOfElements
    TYPE(InterfaceType), POINTER :: INTERFACE
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping
    TYPE(InterfaceMappingRHSType), POINTER :: rhsMapping
    TYPE(InterfaceMatrixType), POINTER :: interfaceMatrix
    TYPE(InterfaceMatrixToVarMapType), POINTER :: interfaceMatrixToVarMap
    TYPE(InterfacePointsConnectivityType), POINTER :: pointsConnectivity
    TYPE(InterfaceRHSType), POINTER :: rhsVector
    TYPE(FieldVariableType), POINTER :: colsFieldVariable,rowsFieldVariable
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("InterfaceMatrices_ElementInitialise",err,error,*999)

    IF(.NOT.ASSOCIATED(interfaceMatrices)) CALL FlagError("Interface matrices is not associated.",err,error,*999)

    NULLIFY(interfaceMapping)
    CALL InterfaceMatrices_InterfaceMappingGet(interfaceMatrices,interfaceMapping,err,error,*999)
    NULLIFY(interfaceEquations)
    CALL InterfaceMatrices_InterfaceEquationsGet(interfaceMatrices,interfaceEquations,err,error,*999)
    NULLIFY(interfaceCondition)
    CALL InterfaceEquations_InterfaceConditionGet(interfaceEquations,interfaceCondition,err,error,*999)
    CALL InterfaceCondition_MethodGet(interfaceCondition,interfaceConditionMethod,err,error,*999)
    CALL InterfaceCondition_IntegrationTypeGet(interfaceCondition,integrationType,err,error,*999)
    
    SELECT CASE(integrationType)
    CASE(INTERFACE_CONDITION_GAUSS_INTEGRATION)
      DO matrixIdx=1,interfaceMatrices%numberOfInterfaceMatrices
        NULLIFY(interfaceMatrix)
        CALL InterfaceMatrices_InterfaceMatrixGet(interfaceMatrices,matrixIdx,interfaceMatrix,err,error,*999)
        NULLIFY(interfaceMatrixToVarMap)
        CALL InterfaceMapping_InterfaceMatrixToVarMapGet(interfaceMapping,matrixIdx,interfaceMatrixToVarMap,err,error,*999)
        NULLIFY(rowsFieldVariable)
        CALL InterfaceMappingIMToVMap_VariableGet(interfaceMatrixToVarMap,rowsFieldVariable,err,error,*999)
        NULLIFY(colsFieldVariable)
        CALL InterfaceMapping_LagrangeVariableGet(interfaceMapping,colsFieldVariable,err,error,*999)
        rowsNumberOfElements=1
        colsNumberOfElements=1
        CALL EquationsMatrices_ElementMatrixSetup(interfaceMatrix%elementMatrix,rowsFieldVariable, &
          & colsFieldVariable,rowsNumberOfElements,colsNumberOfElements,err,error,*999)
      ENDDO !matrixIdx
    CASE(INTERFACE_CONDITION_DATA_POINTS_INTEGRATION)
      NULLIFY(INTERFACE)
      CALL InterfaceCondition_InterfaceGet(interfaceCondition,INTERFACE,err,error,*999)
      NULLIFY(pointsConnectivity)
      CALL Interface_PointsConnectivityGet(INTERFACE,pointsConnectivity,err,error,*999)
      DO matrixIdx=1,interfaceMatrices%numberOfInterfaceMatrices !\todo: Need to separate the case for penalty matrix        
        NULLIFY(interfaceMatrix)
        CALL InterfaceMatrices_InterfaceMatrixGet(interfaceMatrices,matrixIdx,interfaceMatrix,err,error,*999)
        NULLIFY(interfaceMatrixToVarMap)
        CALL InterfaceMapping_InterfaceMatrixToVarMapGet(interfaceMapping,matrixIdx,interfaceMatrixToVarMap,err,error,*999)
        NULLIFY(rowsFieldVariable)
        CALL InterfaceMappingIMToVMap_VariableGet(interfaceMatrixToVarMap,rowsFieldVariable,err,error,*999)
        CALL InterfaceMappingIMToVMap_MeshIndexGet(interfaceMatrixToVarMap,rowsMeshIdx,err,error,*999)
        NULLIFY(colsFieldVariable)
        CALL InterfaceMapping_LagrangeVariableGet(interfaceMapping,colsFieldVariable,err,error,*999)
        CALL InterfacePointsConnectivity_MaximumCoupledElementsGet(pointsConnectivity,rowsMeshIdx,rowsNumberOfElements, &
          & err,error,*999)
        colsNumberOfElements=1        
        CALL EquationsMatrices_ElementMatrixSetup(interfaceMatrix%elementMatrix,rowsFieldVariable, &
          & colsFieldVariable,rowsNumberOfElements,colsNumberOfElements,err,error,*999)
      ENDDO !matrixIdx
    CASE DEFAULT
      localError="The interface condition integration type of "// &
        & TRIM(NumberToVString(integrationType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    rhsVector=>interfaceMatrices%rhsVector
    IF(ASSOCIATED(rhsVector)) THEN
      !Initialise the RHS element vector
      NULLIFY(rhsMapping)
      CALL InterfaceMapping_RHSMappingGet(interfaceMapping,rhsMapping,err,error,*999)
      NULLIFY(rowsFieldVariable)
      CALL InterfaceMappingRHS_RHSVariableGet(rhsMapping,rowsFieldVariable,err,error,*999)
      CALL EquationsMatrices_ElementVectorSetup(rhsVector%elementVector,rowsFieldVariable,err,error,*999)      
    ENDIF
    
    EXITS("InterfaceMatrices_ElementInitialise")
    RETURN
999 ERRORSEXITS("InterfaceMatrices_ElementInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMatrices_ElementInitialise

  !
  !================================================================================================================================
  !

  !>Caclulates the matrix structure (sparsity) for an interface matrix.
  SUBROUTINE InterfaceMatrix_StructureCalculate(interfaceMatrix,numberOfNonZeros,rowIndices,columnIndices, &
    & transposeRowIndices,transposeColumnIndices,err,error,*)

    !Argument variables
    TYPE(InterfaceMatrixType), POINTER :: interfaceMatrix !<A pointer to the interface matrix to calculate the strucute for
    INTEGER(INTG), INTENT(OUT) :: numberOfNonZeros !<On return, the number of non-zeros in the matrix
    INTEGER(INTG), POINTER :: rowIndices(:) !<On return, a pointer to row location indices in compressed row format. The pointer must be NULL on entry and the calling routine is responsible for deallocation.
    INTEGER(INTG), POINTER :: columnIndices(:) !<On return, a pointer to the column location indices in compressed row format. The pointer must be NULL on entry and the calling routine is responsible for deallocation.
    INTEGER(INTG), POINTER :: transposeRowIndices(:) !<On return, if the interface matrix has a transpose a pointer to transpose row location indices in compressed row format. The pointer must be NULL on entry and the calling routine is responsible for deallocation.
    INTEGER(INTG), POINTER :: transposeColumnIndices(:) !<On return, if the interface matrix has a transpose a pointer to the transpose column location indices in compressed row format. The pointer must be NULL on entry and the calling routine is responsible for deallocation.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: columnVersion,columnDerivative,columnIdx,columnComponentIdx,columnLocalDerivativeIdx, &
     & columnInterpolationType,columnLocalNodeIdx,columnNode,dummyErr,domainElement,globalColumn,globalRow,interfaceElementIdx, &
     & interfaceMeshIdx,localColumn,localRow,matrixNumber,numberOfColumns,numberOfColsComponents,numberOfGlobalCols, &
     & numberOfLocalColNodes,numberOfColNodeDerivatives,numberOfLocalRowNodes,numberOfRowNodeDerivatives,numberOfRows, &
     & numberOfRowsComponents,rowComponentIdx,rowInterpolationType,rowVersion,rowDerivative,rowLocalDerivativeIdx,rowIdx, &
     & rowLocalNodeIdx,rowNode,totalNumberOfElements,totalNumberOfLocalCols,totalNumberOfLocalRows,transposeNumberOfNonZeros
    INTEGER(INTG), ALLOCATABLE :: columns(:),transposeColumns(:)
    REAL(DP) :: sparsity
    TYPE(BasisType), POINTER :: columnBasis,rowBasis
    TYPE(DomainType), POINTER :: columnDomain,rowDomain
    TYPE(DomainElementsType), POINTER :: columnDomainElements,rowDomainElements
    TYPE(DomainMappingType), POINTER :: columnDOFSDomainMapping,rowDOFSDomainMapping
    TYPE(DomainTopologyType), POINTER :: columnDomainTopology,rowDomainTopology
    TYPE(InterfaceType), POINTER :: INTERFACE
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices
    TYPE(InterfaceMatrixToVarMapType), POINTER :: interfaceMatrixToVarMap
    TYPE(InterfaceMeshConnectivityType), POINTER :: meshConnectivity
    TYPE(FieldVariableType), POINTER :: columnVariable,rowVariable
    TYPE(ListPtrType), ALLOCATABLE :: columnIndicesLists(:)
    TYPE(ListPtrType), ALLOCATABLE :: transposeColumnIndicesLists(:)
    TYPE(VARYING_STRING) :: dummyError,localError

    ENTERS("InterfaceMatrix_StructureCalculate",err,error,*999)

    numberOfNonZeros=0
    IF(ASSOCIATED(rowIndices)) CALL FlagError("Row indices is already associated.",err,error,*998)
    IF(ASSOCIATED(columnIndices)) CALL FlagError("Column indices is already associated.",err,error,*998)
    IF(ASSOCIATED(transposeRowIndices)) CALL FlagError("Transpose row indices is already associated.",err,error,*998)
    IF(ASSOCIATED(transposeColumnIndices)) CALL FlagError("Transpose column indices is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMatrix)) CALL FlagError("Interface matrix is not associated.",err,error,*998)
    
    matrixNumber=interfaceMatrix%matrixNumber
    SELECT CASE(interfaceMatrix%structureType)
    CASE(INTERFACE_MATRIX_NO_STRUCTURE)
      CALL FlagError("There is no structure to calculate for a matrix with no structure.",err,error,*998)
    CASE(INTERFACE_MATRIX_FEM_STRUCTURE)
      SELECT CASE(interfaceMatrix%storageType)
      CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
        NULLIFY(interfaceMatrices)
        CALL InterfaceMatrix_InterfaceMatricesGet(interfaceMatrix,interfaceMatrices,err,error,*999)
        NULLIFY(interfaceMapping)
        CALL InterfaceMatrices_InterfaceMappingGet(interfaceMatrices,interfaceMapping,err,error,*999)
        NULLIFY(interfaceMatrixToVarMap)
        CALL InterfaceMapping_InterfaceMatrixToVarMapGet(interfaceMapping,matrixNumber,interfaceMatrixToVarMap,err,error,*999)
        NULLIFY(interfaceEquations)
        CALL InterfaceMatrices_InterfaceEquationsGet(interfaceMatrices,interfaceEquations,err,error,*999)
        NULLIFY(interfaceCondition)
        CALL InterfaceEquations_InterfaceConditionGet(interfaceEquations,interfaceCondition,err,error,*999)
        NULLIFY(INTERFACE)
        CALL InterfaceCondition_InterfaceGet(interfaceCondition,INTERFACE,err,error,*999)
        NULLIFY(meshConnectivity)
        CALL Interface_MeshConnectivityGet(INTERFACE,meshConnectivity,err,error,*999)
        NULLIFY(rowVariable)
        CALL InterfaceMapping_MatrixVariableGet(interfaceMapping,matrixNumber,rowVariable,err,error,*999)
        NULLIFY(columnVariable)
        CALL InterfaceMapping_LagrangeVariableGet(interfaceMapping,columnVariable,err,error,*999)
        CALL InterfaceMappingIMToVMap_MeshIndexGet(interfaceMatrixToVarMap,interfaceMeshIdx,err,error,*999)
        NULLIFY(rowDOFSDomainMapping)
        CALL FieldVariable_DomainMappingGet(rowVariable,rowDOFSDomainMapping,err,error,*999)
        NULLIFY(columnDOFSDomainMapping)
        CALL FieldVariable_DomainMappingGet(columnVariable,columnDOFSDomainMapping,err,error,*999)
        
        !Allocate lists
        CALL DomainMapping_TotalNumberOfLocalGet(rowDOFSDomainMapping,totalNumberOfLocalRows,err,error,*999)
        ALLOCATE(columnIndicesLists(totalNumberOfLocalRows),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate column indices lists.",err,error,*999)
        DO localRow=1,totalNumberOfLocalRows
          !Set up list
          NULLIFY(columnIndicesLists(localRow)%ptr)
          CALL List_CreateStart(columnIndicesLists(localRow)%ptr,err,error,*999)
          CALL List_DataTypeSet(columnIndicesLists(localRow)%ptr,LIST_INTG_TYPE,err,error,*999)
          CALL List_InitialSizeSet(columnIndicesLists(localRow)%ptr,50,err,error,*999)
          CALL List_CreateFinish(columnIndicesLists(localRow)%ptr,err,error,*999)
        ENDDO !localRow
        !Allocate row indices
        ALLOCATE(rowIndices(totalNumberOfLocalRows+1),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate row indices.",err,error,*999)
        IF(interfaceMatrix%hasTranspose) THEN
          !Allocate transpose lists
          CALL DomainMapping_TotalNumberOfLocalGet(columnDOFsDomainMapping,totalNumberOfLocalCols,err,error,*999)
          ALLOCATE(transposeColumnIndicesLists(totalNumberOfLocalCols),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate transpose column indices lists.",err,error,*999)
          DO localColumn=1,totalNumberOfLocalCols
            !Set up list
            NULLIFY(transposeColumnIndicesLists(localColumn)%ptr)
            CALL List_CreateStart(transposeColumnIndicesLists(localColumn)%ptr,err,error,*999)
            CALL List_DataTypeSet(transposeColumnIndicesLists(localColumn)%ptr,LIST_INTG_TYPE,err,error,*999)
            CALL List_InitialSizeSet(transposeColumnIndicesLists(localColumn)%ptr,50,err,error,*999)
            CALL List_CreateFinish(transposeColumnIndicesLists(localColumn)%ptr,err,error,*999)
          ENDDO !localColumn
          !Allocate transpose row indices
          ALLOCATE(transposeRowIndices(totalNumberOfLocalCols+1),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate transpose row indices.",err,error,*999)
        ENDIF
        !Loop over the number of components in the Lagrange multipler variable
        CALL FieldVariable_NumberOfComponentsGet(rowVariable,numberOfRowsComponents,err,error,*999)
        CALL FieldVariable_NumberOfComponentsGet(columnVariable,numberOfColsComponents,err,error,*999)
        DO columnComponentIdx=1,numberOfColsComponents          
          CALL FieldVariable_ComponentInterpolationGet(columnVariable,columnComponentIdx,columnInterpolationType,err,error,*999)
          IF(columnInterpolationType/=FIELD_NODE_BASED_INTERPOLATION) &
            & CALL FlagError("Only node based fields implemented.",err,error,*999)
          !Loop over the elements in the interface mesh
          NULLIFY(columnDomain)
          CALL FieldVariable_DomainGet(columnVariable,columnComponentIdx,columnDomain,err,error,*999)
          NULLIFY(columnDomainTopology)
          CALL Domain_DomainTopologyGet(columnDomain,columnDomainTopology,err,error,*999)
          NULLIFY(columnDomainElements)
          CALL DomainTopology_DomainElementsGet(columnDomainTopology,columnDomainElements,err,error,*999)
          CALL DomainElements_TotalNumberOfElementsGet(columnDomainElements,totalNumberOfElements,err,error,*999)
          DO interfaceElementIdx=1,totalNumberOfElements
            NULLIFY(columnBasis)
            CALL DomainElements_ElementBasisGet(columnDomainElements,interfaceElementIdx,columnBasis,err,error,*999)
            !Loop over the column DOFs in the element
            CALL Basis_NumberOfLocalNodesGet(columnBasis,numberOfLocalColNodes,err,error,*999)
            DO columnLocalNodeIdx=1,numberOfLocalColNodes
              CALL DomainElements_ElementNodeGet(columnDomainElements,columnLocalNodeIdx,interfaceElementIdx,columnNode, &
                & err,error,*999)
              CALL Basis_NodeNumberOfDerivativesGet(columnBasis,columnLocalNodeIdx,numberOfColNodeDerivatives,err,error,*999)
              DO columnLocalDerivativeIdx=1,numberOfColNodeDerivatives
                CALL DomainElements_ElementDerivativeGet(columnDomainElements,columnLocalDerivativeIdx,columnLocalNodeIdx, &
                  & interfaceElementIdx,columnDerivative,err,error,*999)
                CALL DomainElements_ElementVersionGet(columnDomainElements,columnLocalDerivativeIdx,columnLocalNodeIdx, &
                  & interfaceElementIdx,columnVersion,err,error,*999)
                CALL FieldVariable_LocalNodeDOFGet(columnVariable,columnVersion,columnDerivative,columnNode, &
                  & columnComponentIdx,localColumn,err,error,*999)
                CALL DomainMapping_LocalToGlobalGet(columnDOFsDomainMapping,localColumn,globalColumn,err,error,*999)
               !Loop over the components in the dependent variable
                DO rowComponentIdx=1,numberOfRowsComponents
                  CALL FieldVariable_ComponentInterpolationGet(rowVariable,rowComponentIdx,rowInterpolationType,err,error,*999)
                  SELECT CASE(rowInterpolationType)
                  CASE(FIELD_CONSTANT_INTERPOLATION)
                    CALL FieldVariable_ConstantDOFGet(rowVariable,rowComponentIdx,localRow,err,error,*999)
                    CALL List_ItemAdd(columnIndicesLists(localRow)%ptr,globalColumn,err,error,*999)
                    IF(interfaceMatrix%hasTranspose) THEN
                      CALL DomainMapping_LocalToGlobalGet(rowDOFsDomainMapping,localRow,globalRow,err,error,*999)
                      CALL List_ItemAdd(transposeColumnIndicesLists(localColumn)%ptr,globalRow,err,error,*999)
                    ENDIF
                  CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                    CALL InterfaceMeshConnectivity_CoupledElementNumberGet(meshConnectivity,interfaceElementIdx,interfaceMeshIdx, &
                      & domainElement,err,error,*999)
                    CALL FieldVariable_LocalElementDOFGet(rowVariable,domainElement,rowComponentIdx,localRow,err,error,*999)
                    CALL List_ItemAdd(columnIndicesLists(localRow)%ptr,globalColumn,err,error,*999)
                    IF(interfaceMatrix%hasTranspose) THEN
                      CALL DomainMapping_LocalToGlobalGet(rowDOFsDomainMapping,localRow,globalRow,err,error,*999)
                      CALL List_ItemAdd(transposeColumnIndicesLists(localColumn)%ptr,globalRow,err,error,*999)
                    ENDIF
                  CASE(FIELD_NODE_BASED_INTERPOLATION)
                    NULLIFY(rowDomain)
                    CALL FieldVariable_DomainGet(rowVariable,rowComponentIdx,rowDomain,err,error,*999)
                    NULLIFY(rowDomainTopology)
                    CALL Domain_DomainTopologyGet(rowDomain,rowDomainTopology,err,error,*999)
                    NULLIFY(rowDomainElements)
                    CALL DomainTopology_DomainElementsGet(rowDomainTopology,rowDomainElements,err,error,*999)
                    CALL InterfaceMeshConnectivity_CoupledElementNumberGet(meshConnectivity,interfaceElementIdx,interfaceMeshIdx, &
                      & domainElement,err,error,*999)
                    NULLIFY(rowBasis)
                    CALL DomainElements_ElementBasisGet(rowDomainElements,domainElement,rowBasis,err,error,*999)
                    !Loop over the row DOFs in the domain mesh element
                    CALL Basis_NumberOfLocalNodesGet(rowBasis,numberOfLocalRowNodes,err,error,*999)
                    DO rowLocalNodeIdx=1,numberOfLocalRowNodes
                      CALL DomainElements_ElementNodeGet(rowDomainElements,rowLocalNodeIdx,domainElement,rowNode, &
                        & err,error,*999)
                      CALL Basis_NodeNumberOfDerivativesGet(rowBasis,rowLocalNodeIdx,numberOfRowNodeDerivatives,err,error,*999)
                      DO rowLocalDerivativeIdx=1,numberOfRowNodeDerivatives
                        CALL DomainElements_ElementDerivativeGet(rowDomainElements,rowLocalDerivativeIdx,rowLocalNodeIdx, &
                          & domainElement,rowDerivative,err,error,*999)
                        CALL DomainElements_ElementVersionGet(rowDomainElements,rowLocalDerivativeIdx,rowLocalNodeIdx, &
                          & domainElement,rowVersion,err,error,*999)
                        CALL FieldVariable_LocalNodeDOFGet(rowVariable,rowVersion,rowDerivative,rowNode,rowComponentIdx, &
                          & localRow,err,error,*999)
                        CALL List_ItemAdd(columnIndicesLists(localRow)%ptr,globalColumn,err,error,*999)
                        IF(interfaceMatrix%hasTranspose) THEN
                          CALL DomainMapping_LocalToGlobalGet(rowDOFsDomainMapping,localRow,globalRow,err,error,*999)
                          CALL List_ItemAdd(transposeColumnIndicesLists(localColumn)%ptr,globalRow,err,error,*999)
                        ENDIF
                      ENDDO !rowLocalDerivativeIdx
                    ENDDO !rowLocalNodeIdx
                  CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                    CALL FlagError("Not implemented.",err,error,*999)
                  CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                    CALL FlagError("Not implemented.",err,error,*999)
                  CASE DEFAULT
                    localError="The row variable interpolation type of "// &
                      & TRIM(NumberToVString(rowInterpolationType,"*",err,error))//" is invalid."
                    CALL FlagError(localError,err,error,*999)
                  END SELECT
                ENDDO !rowComponentIdx
              ENDDO !columnLocalDerivativeIdx
            ENDDO !columnLocalNodeIdx
          ENDDO !interfaceElementIdx
        ENDDO !columnComponentIdx
        rowIndices(1)=1
        DO localRow=1,totalNumberOfLocalRows
          CALL List_RemoveDuplicates(columnIndicesLists(localRow)%ptr,err,error,*999)
          CALL List_NumberOfItemsGet(columnIndicesLists(localRow)%ptr,numberOfColumns,err,error,*999)
          numberOfNonZeros=numberOfNonZeros+numberOfColumns
          rowIndices(localRow+1)=numberOfNonZeros+1
        ENDDO !localRow
        IF(interfaceMatrix%hasTranspose) THEN
          transposeNumberOfNonZeros=0
          transposeRowIndices(1)=1
          DO localColumn=1,totalNumberOfLocalCols
            CALL List_RemoveDuplicates(transposeColumnIndicesLists(localColumn)%ptr,err,error,*999)
            CALL List_NumberOfItemsGet(transposeColumnIndicesLists(localColumn)%ptr,numberOfColumns,err,error,*999)
            transposeNumberOfNonZeros=transposeNumberOfNonZeros+numberOfColumns
            transposeRowIndices(localColumn+1)=transposeNumberOfNonZeros+1
          ENDDO !localColumn
          !Sanity check - the number of non-zeros should be the same
          IF(transposeNumberOfNonZeros/=numberOfNonZeros) THEN
            localError="Invalid number of non-zeros. The number of non-zeros in the "// &
              & "transposed matrix ("//TRIM(NumberToVString(transposeNumberOfNonZeros, &
              & "*",err,error))//") does not match the number of non-zeros in the interface "// &
              & "matrix ("//TRIM(NumberToVString(numberOfNonZeros,"*",err,error))//")."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ENDIF
        !Allocate and setup the column locations
        ALLOCATE(columnIndices(numberOfNonZeros),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate column indices.",err,error,*999)
        DO localRow=1,totalNumberOfLocalRows
          CALL List_DetachAndDestroy(columnIndicesLists(localRow)%ptr,numberOfColumns,columns,err,error,*999)
          DO columnIdx=1,numberOfColumns
            columnIndices(rowIndices(localRow)+columnIdx-1)=columns(columnIdx)
          ENDDO !columnIdx
          DEALLOCATE(columns)
        ENDDO !localRow
        IF(interfaceMatrix%hasTranspose) THEN
          !Allocate and setup the column locations
          ALLOCATE(transposeColumnIndices(numberOfNonZeros),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate transpose column indices.",err,error,*999)
          DO localColumn=1,totalNumberOfLocalCols
            CALL List_DetachAndDestroy(transposeColumnIndicesLists(localColumn)%ptr,numberOfRows,transposeColumns,err,error,*999)
            DO rowIdx=1,numberOfRows
              transposeColumnIndices(transposeRowIndices(localColumn)+rowIdx-1)=transposeColumns(rowIdx)
            ENDDO !rowIdx
            DEALLOCATE(transposeColumns)
          ENDDO !localColumn
        ENDIF
        
        IF(diagnostics1) THEN
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Interface matrix structure:",err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Interface matrix number : ",matrixNumber,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of rows = ",totalNumberOfLocalRows,err,error,*999)
          CALL DomainMapping_NumberOfGlobalGet(columnDOFsDomainMapping,numberOfGlobalCols,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of columns = ",numberOfGlobalCols,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of non zeros = ",numberOfNonZeros,err,error,*999)
          IF(totalNumberOfLocalRows*numberOfGlobalCols/=0) THEN
            sparsity=REAL(numberOfNonZeros,DP)/REAL(totalNumberOfLocalRows*numberOfGlobalCols,DP)*100.0_DP
            CALL WriteStringFmtValue(DIAGNOSTIC_OUTPUT_TYPE,"  Sparsity (%) = ",sparsity,"F6.2",err,error,*999)
          ENDIF
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,totalNumberOfLocalRows+1,5,5,rowIndices, &
            & '("  Row indices              :",5(X,I13))','(28X,5(X,I13))',err,error,*999)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfNonZeros,8,8,columnIndices, &
            & '("  Column indices           :",5(X,I13))','(28X,5(X,I13))',err,error,*999)
          IF(interfaceMatrix%hasTranspose) THEN 
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,totalNumberOfLocalCols+1,5,5, &
              & transposeRowIndices,'("  Transpose row indices    :",5(X,I13))','(28X,5(X,I13))',err,error,*999)
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfNonZeros,8,8,transposeColumnIndices, &
              & '("  Transpose column indices :",5(X,I13))','(28X,5(X,I13))',err,error,*999)
          ENDIF
        ENDIF
        
      CASE DEFAULT
        localError="The matrix storage type of "//TRIM(NumberToVString(interfaceMatrix%storageType,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The matrix structure type of "//TRIM(NumberToVString(interfaceMatrix%structureType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*998)
    END SELECT
    
    EXITS("InterfaceMatrix_StructureCalculate")
    RETURN
999 IF(ASSOCIATED(rowIndices)) DEALLOCATE(rowIndices)
    IF(ASSOCIATED(columnIndices)) DEALLOCATE(columnIndices)
    IF(ASSOCIATED(transposeRowIndices)) DEALLOCATE(transposeRowIndices)
    IF(ASSOCIATED(transposeColumnIndices)) DEALLOCATE(transposeColumnIndices)
    IF(ALLOCATED(columns)) DEALLOCATE(columns)
    IF(ALLOCATED(transposeColumns)) DEALLOCATE(transposeColumns)
    IF(ALLOCATED(columnIndicesLists)) THEN
      DO localRow=1,SIZE(columnIndicesLists,1)
        IF(ASSOCIATED(columnIndicesLists(localRow)%ptr)) &
          & CALL List_Destroy(columnIndicesLists(localRow)%ptr,dummyErr,dummyError,*998)
      ENDDO !localRow
      DEALLOCATE(columnIndicesLists)
    ENDIF
    IF(ALLOCATED(transposeColumnIndicesLists)) THEN
      DO localColumn=1,SIZE(transposeColumnIndicesLists,1)
        IF(ASSOCIATED(transposeColumnIndicesLists(localColumn)%ptr)) &
          & CALL List_Destroy(transposeColumnIndicesLists(localColumn)%ptr,dummyErr,dummyError,*998)
      ENDDO !localRow
      DEALLOCATE(transposeColumnIndicesLists)
    ENDIF
998 ERRORSEXITS("InterfaceMatrix_StructureCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMatrix_StructureCalculate

  !
  !================================================================================================================================
  !

  !>Finishes the creation of the interface matrices for the interface equations
  SUBROUTINE InterfaceMatrices_CreateFinish(interfaceMatrices,err,error,*)

    !Argument variables
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices !<The pointer to the interface matrices
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string  
    !Local Variables
    INTEGER(INTG) :: dummyErr,matrixIdx,numberOfNonZeros
    INTEGER(INTG), POINTER :: rowIndices(:),columnIndices(:),transposeRowIndices(:),transposeColumnIndices(:)
    TYPE(DomainMappingType), POINTER :: rowDomainMapping,columnDomainMapping
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping
    TYPE(InterfaceMatrixType), POINTER :: interfaceMatrix
    TYPE(InterfaceMatrixToVarMapType), POINTER :: interfaceMatrixToVarMap
    TYPE(InterfaceRHSType), POINTER :: rhsVector
    TYPE(VARYING_STRING) :: dummyError,localError

    ENTERS("InterfaceMatrices_CreateFinish",err,error,*998)

    CALL InterfaceMatrices_AssertNotFinished(interfaceMatrices,err,error,*999)

    NULLIFY(interfaceMapping)
    CALL InterfaceMatrices_InterfaceMappingGet(interfaceMatrices,interfaceMapping,err,error,*999)
    NULLIFY(columnDomainMapping)
    CALL InterfaceMapping_ColumnDOFsMappingGet(interfaceMapping,columnDomainMapping,err,error,*999)
    !Now create the individual interface matrices
    DO matrixIdx=1,interfaceMatrices%numberOfInterfaceMatrices
      NULLIFY(interfaceMatrix)
      CALL InterfaceMatrices_InterfaceMatrixGet(interfaceMatrices,matrixIdx,interfaceMatrix,err,error,*999)
      NULLIFY(interfaceMatrixToVarMap)
      CALL InterfaceMapping_InterfaceMatrixToVarMapGet(interfaceMapping,matrixIdx,interfaceMatrixToVarMap,err,error,*999)
      NULLIFY(rowDomainMapping)
      CALL InterfaceMappingIMToVMap_RowDOFsMappingGet(interfaceMatrixToVarMap,rowDomainMapping,err,error,*999)

      !Create the distributed equations matrix
      CALL DistributedMatrix_CreateStart(rowDomainMapping,columnDomainMapping,interfaceMatrices%matrices(matrixIdx)%ptr%matrix, &
        & err,error,*999)
      CALL DistributedMatrix_DataTypeSet(interfaceMatrix%matrix,DISTRIBUTED_MATRIX_VECTOR_DP_TYPE,err,error,*999)
      CALL DistributedMatrix_StorageTypeSet(interfaceMatrix%matrix,interfaceMatrix%storageType,err,error,*999)
      IF(interfaceMatrix%hasTranspose) THEN
        CALL DistributedMatrix_TransposeTypeSet(interfaceMatrix%matrix,DISTRIBUTED_MATRIX_FULL_TRANSPOSE_REQUIRED,err,error,*999)
        CALL DistributedMatrix_CreateStart(columnDomainMapping,rowDomainMapping,interfaceMatrices%matrices(matrixIdx)%ptr% &
          & matrixTranspose,err,error,*999)
        CALL DistributedMatrix_DataTypeSet(interfaceMatrix%matrixTranspose,DISTRIBUTED_MATRIX_VECTOR_DP_TYPE,err,error,*999)
        CALL DistributedMatrix_StorageTypeSet(interfaceMatrix%matrixTranspose,interfaceMatrix%storageType,err,error,*999)
      ENDIF
      !Calculate and set the matrix structure/sparsity pattern
      IF(interfaceMatrix%storageType/=DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE.AND. &
        & interfaceMatrix%storageType/=DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE) THEN
        CALL InterfaceMatrix_StructureCalculate(interfaceMatrix,numberOfNonZeros,rowIndices,columnIndices, &
          & transposeRowIndices,transposeColumnIndices,err,error,*999)
        CALL DistributedMatrix_NumberOfNonZerosSet(interfaceMatrix%matrix,numberOfNonZeros,err,error,*999)
        CALL DistributedMatrix_StorageLocationsSet(interfaceMatrix%matrix,rowIndices,columnIndices,err,error,*999)
        IF(interfaceMatrix%hasTranspose) THEN
          CALL DistributedMatrix_NumberOfNonZerosSet(interfaceMatrix%matrixTranspose,numberOfNonZeros,err,error,*999)
          CALL DistributedMatrix_StorageLocationsSet(interfaceMatrix%matrixTranspose,transposeRowIndices, &
            & transposeColumnIndices,err,error,*999)
        ENDIF
        IF(ASSOCIATED(rowIndices)) DEALLOCATE(rowIndices)
        IF(ASSOCIATED(columnIndices)) DEALLOCATE(columnIndices)
        IF(ASSOCIATED(transposeRowIndices)) DEALLOCATE(transposeRowIndices)
        IF(ASSOCIATED(transposeColumnIndices)) DEALLOCATE(transposeColumnIndices)
      ENDIF
      CALL DistributedMatrix_CreateFinish(interfaceMatrix%matrix,err,error,*999)
      !Create the temp vector
      CALL DistributedVector_CreateStart(rowDomainMapping,interfaceMatrix%tempVector,err,error,*999)
      CALL DistributedVector_DataTypeSet(interfaceMatrix%tempVector,MATRIX_VECTOR_DP_TYPE,err,error,*999)
      CALL DistributedVector_CreateFinish(interfaceMatrix%tempVector,err,error,*999)
      IF(interfaceMatrix%hasTranspose) THEN
        CALL DistributedMatrix_CreateFinish(interfaceMatrix%matrixTranspose,err,error,*999)
        !Create the temp vector
        CALL DistributedVector_CreateStart(columnDomainMapping,interfaceMatrix%tempTransposeVector,err,error,*999)
        CALL DistributedVector_DataTypeSet(interfaceMatrix%tempTransposeVector,MATRIX_VECTOR_DP_TYPE,err,error,*999)
        CALL DistributedVector_CreateFinish(interfaceMatrix%tempTransposeVector,err,error,*999)
      ENDIF
    ENDDO !matrixIdx
    rhsVector=>interfaceMatrices%rhsVector
    IF(ASSOCIATED(rhsVector)) THEN
      !Set up the interface RHS vector
      CALL DistributedVector_CreateStart(columnDomainMapping,interfaceMatrices%rhsVector%rhsVector,err,error,*999)
      CALL DistributedVector_DataTypeSet(rhsVector%rhsVector,MATRIX_VECTOR_DP_TYPE,err,error,*999)
      CALL DistributedVector_CreateFinish(rhsVector%rhsVector,err,error,*999)
    ENDIF
    !Finish up
    interfaceMatrices%interfaceMatricesFinished=.TRUE.
       
    EXITS("InterfaceMatrices_CreateFinish")
    RETURN
999 IF(ASSOCIATED(rowIndices)) DEALLOCATE(rowIndices)
    IF(ASSOCIATED(columnIndices)) DEALLOCATE(columnIndices)
    IF(ASSOCIATED(transposeRowIndices)) DEALLOCATE(transposeRowIndices)
    IF(ASSOCIATED(transposeColumnIndices)) DEALLOCATE(transposeColumnIndices)
    CALL InterfaceMatrices_Finalise(interfaceMatrices,dummyErr,dummyError,*998)
998 ERRORSEXITS("InterfaceMatrices_CreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMatrices_CreateFinish

  !
  !================================================================================================================================
  !

  !>Starts the creation of the interface matrices and rhs for the interface equations
  SUBROUTINE InterfaceMatrices_CreateStart(interfaceEquations,interfaceMatrices,err,error,*)

    !Argument variables
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations !<The pointer to the interface equations to create the interface equations matrices for
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices !<On return, a pointer to the interface matrices being created. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string  
    !Local Variables

    ENTERS("InterfaceMatrices_CreateStart",err,error,*998)

    IF(ASSOCIATED(interfaceMatrices)) CALL FlagError("Interface matrices is already associated.",err,error,*998)
    CALL InterfaceEquations_AssertIsFinished(interfaceEquations,err,error,*999)
    
    !Initialise the interface matrices
    CALL InterfaceMatrices_Initialise(interfaceEquations,err,error,*999)
    !Return the pointer
    interfaceMatrices=>interfaceEquations%interfaceMatrices
   
    EXITS("InterfaceMatrices_CreateStart")
    RETURN
999 NULLIFY(interfaceMatrices)
998 ERRORSEXITS("InterfaceMatrices_CreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMatrices_CreateStart

  !
  !================================================================================================================================
  !

  !>Destroy the interface matrices
  SUBROUTINE InterfaceMatrices_Destroy(interfaceMatrices,err,error,*)

    !Argument variables
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices !<A pointer the interface matrices to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("InterfaceMatrices_Destroy",err,error,*999)

    IF(.NOT.ASSOCIATED(interfaceMatrices)) CALL FlagError("Interface matrices is not associated.",err,error,*999)
    
    CALL InterfaceMatrices_Finalise(interfaceMatrices,err,error,*999)
        
    EXITS("InterfaceMatrices_Destroy")
    RETURN
999 ERRORSEXITS("InterfaceMatrices_Destroy",err,error)    
    RETURN 1
   
  END SUBROUTINE InterfaceMatrices_Destroy

  !
  !================================================================================================================================
  !

  !>Finalise the interface matrices and deallocate all memory.
  SUBROUTINE InterfaceMatrices_Finalise(interfaceMatrices,err,error,*)

    !Argument variables
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices !<A pointer to the interface matrices
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
   
    ENTERS("InterfaceMatrices_Finalise",err,error,*999)

    IF(ASSOCIATED(interfaceMatrices)) THEN
      IF(ALLOCATED(interfaceMatrices%matrices)) THEN
        DO matrixIdx=1,SIZE(interfaceMatrices%matrices,1)
          CALL InterfaceMatrix_Finalise(interfaceMatrices%matrices(matrixIdx)%ptr,err,error,*999)
        ENDDO !matrixIdx
        DEALLOCATE(interfaceMatrices%matrices)
      ENDIF
      CALL InterfaceMatrices_RHSFinalise(interfaceMatrices%rhsVector,err,error,*999)
      DEALLOCATE(interfaceMatrices)
    ENDIF
       
    EXITS("InterfaceMatrices_Finalise")
    RETURN
999 ERRORSEXITS("InterfaceMatrices_Finalise",err,error)
    RETURN 1
  END SUBROUTINE InterfaceMatrices_Finalise

  !
  !================================================================================================================================
  !

  !>Initialise the interface matrices for the interface equations.
  SUBROUTINE InterfaceMatrices_Initialise(interfaceEquations,err,error,*)
    
     !Argument variables
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations !<A pointer to the interface equations to initialise the interface matrices for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,matrixIdx
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping
    TYPE(InterfaceMappingRHSType), POINTER :: rhsInterfaceMapping
    TYPE(InterfaceMatrixToVarMapType), POINTER :: interfaceMatrixToVarMap
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("InterfaceMatrices_Initialise",err,error,*998)

    IF(.NOT.ASSOCIATED(interfaceEquations)) CALL FlagError("Interface equations is not associated.",err,error,*998)
    IF(ASSOCIATED(interfaceEquations%interfaceMatrices)) &
      & CALL FlagError("Interface matrices is already associated for this interface equations.",err,error,*998)

    NULLIFY(interfaceMapping)
    CALL InterfaceEquations_InterfaceMappingGet(interfaceEquations,interfaceMapping,err,error,*999)
    CALL InterfaceMapping_AssertIsFinished(interfaceMapping,err,error,*999)
    
    ALLOCATE(interfaceEquations%interfaceMatrices,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate interface equations interface matrices.",err,error,*999)
    interfaceEquations%interfaceMatrices%interfaceEquations=>interfaceEquations
    interfaceEquations%interfaceMatrices%interfaceMatricesFinished=.FALSE.
    interfaceEquations%interfaceMatrices%interfaceMapping=>interfaceMapping
    NULLIFY(interfaceEquations%interfaceMatrices%solverMapping)
    CALL InterfaceMapping_NumberOfColumnsGet(interfaceMapping,interfaceEquations%interfaceMatrices%numberOfColumns,err,error,*999)
    CALL InterfaceMapping_TotalNumberOfColumnsGet(interfaceMapping,interfaceEquations%interfaceMatrices%totalNumberOfColumns, &
      & err,error,*999)
    CALL InterfaceMapping_NumberOfGlobalColumnsGet(interfaceMapping,interfaceEquations%interfaceMatrices%numberOfGlobalColumns, &
      & err,error,*999)
    NULLIFY(interfaceEquations%interfaceMatrices%rhsVector)
    !Allocate and initialise the matrices
    CALL InterfaceMapping_NumberOfInterfaceMatricesGet(interfaceMapping,interfaceEquations%interfaceMatrices% &
      & numberOfInterfaceMatrices,err,error,*999)
    ALLOCATE(interfaceEquations%interfaceMatrices%matrices(interfaceEquations%interfaceMatrices% &
      & numberOfInterfaceMatrices),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate interface matrices matrices.",err,error,*999)
    DO matrixIdx=1,interfaceEquations%interfaceMatrices%numberOfInterfaceMatrices
      NULLIFY(interfaceEquations%interfaceMatrices%matrices(matrixIdx)%ptr)
      CALL InterfaceMatrix_Initialise(interfaceEquations%interfaceMatrices,matrixIdx,err,error,*999)
      NULLIFY(interfaceMatrixToVarMap)
      CALL InterfaceMapping_InterfaceMatrixToVarMapGet(interfaceMapping,matrixIdx,interfaceMatrixToVarMap,err,error,*999)
      CALL InterfaceMappingIMToVMap_MatrixCoefficientGet(interfaceMatrixToVarMap,interfaceEquations%interfaceMatrices% &
        & matrices(matrixIdx)%ptr%matrixCoefficient,err,error,*999)
      CALL InterfaceMappingIMToVMap_TransposeMatrixCoefficientGet(interfaceMatrixToVarMap,interfaceEquations%interfaceMatrices% &
        & matrices(matrixIdx)%ptr%transposeMatrixCoefficient,err,error,*999)
    ENDDO !matrixIdx
    CALL InterfaceMatrices_RHSInitialise(interfaceEquations%interfaceMatrices,err,error,*999)
    CALL InterfaceMapping_RHSMappingGet(interfaceMapping,rhsInterfaceMapping,err,error,*999)
    CALL InterfaceMappingRHS_VectorCoefficientGet(rhsInterfaceMapping,interfaceEquations%interfaceMatrices%rhsVector% &
      & rhsCoefficient,err,error,*999)
       
    EXITS("InterfaceMatrices_Initialise")
    RETURN
999 CALL InterfaceMatrices_Finalise(interfaceEquations%interfaceMatrices,dummyErr,dummyError,*998)
998 ERRORSEXITS("InterfaceMatrices_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMatrices_Initialise

  !
  !================================================================================================================================
  !

  !>Outputs the interface matrices
  SUBROUTINE InterfaceMatrices_Output(id,interfaceMatrices,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the ouptut stream
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices !<A pointer to the interface matrices to output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(InterfaceMatrixType), POINTER :: interfaceMatrix
    TYPE(InterfaceRHSType), POINTER :: rhsVector
    
    ENTERS("InterfaceMatrices_Output",err,error,*999)

    CALL InterfaceMatrices_AssertIsFinished(interfaceMatrices,err,error,*999)

    CALL WriteString(id,"",err,error,*999)
    CALL WriteString(id,"Interface matrices:",err,error,*999)
    CALL WriteStringValue(id,"Number of interface matrices = ",interfaceMatrices%numberOfInterfaceMatrices,err,error,*999)
    DO matrixIdx=1,interfaceMatrices%numberOfInterfaceMatrices
      NULLIFY(interfaceMatrix)
      CALL InterfaceMatrices_InterfaceMatrixGet(interfaceMatrices,matrixIdx,interfaceMatrix,err,error,*999)
      CALL WriteStringValue(id,"Interface matrix : ",matrixIdx,err,error,*999)
      CALL WriteString(id,"Standard matrix:",err,error,*999)
      CALL DistributedMatrix_Output(id,interfaceMatrix%matrix,err,error,*999)
      IF(interfaceMatrix%hasTranspose) THEN
        CALL WriteString(id,"Transposed matrix:",err,error,*999)
        CALL DistributedMatrix_Output(id,interfaceMatrix%matrixTranspose,err,error,*999)
      ENDIF
    ENDDO !matrixIdx
    rhsVector=>interfaceMatrices%rhsVector
    IF(ASSOCIATED(rhsVector)) THEN
      CALL WriteString(id,"Interface RHS vector:",err,error,*999)
      CALL DistributedVector_Output(id,rhsVector%rhsVector,err,error,*999)
    ENDIF
   
    EXITS("InterfaceMatrices_Output")
    RETURN
999 ERRORSEXITS("InterfaceMatrices_Output",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMatrices_Output
  
  !
  !================================================================================================================================
  !

  !>Finalises the interface matrices RHS vector and deallocates all memory
  SUBROUTINE InterfaceMatrices_RHSFinalise(rhsVector,err,error,*)

    !Argument variables
    TYPE(InterfaceRHSType), POINTER :: rhsVector !<A pointer to the equation matrices RHS vector to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
     
    ENTERS("InterfaceMatrices_RHSFinalise",err,error,*999)

    IF(ASSOCIATED(rhsVector)) THEN
      IF(ASSOCIATED(rhsVector%rhsVector)) CALL DistributedVector_Destroy(rhsVector%rhsVector,err,error,*999)
      CALL EquationsMatrices_ElementVectorFinalise(rhsVector%elementVector,err,error,*999)
      DEALLOCATE(rhsVector)
    ENDIF      
     
    EXITS("InterfaceMatrices_RHSFinalise")
    RETURN
999 ERRORSEXITS("InterfaceMatrices_RHSFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMatrices_RHSFinalise
  
  !
  !================================================================================================================================
  !

  !>Initialises the interface matrices RHS vector
  SUBROUTINE InterfaceMatrices_RHSInitialise(interfaceMatrices,err,error,*)

    !Argument variables
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices !<A pointer to the equation matrices to initialise the rhs vector for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping
    TYPE(InterfaceMappingRHSType), POINTER :: rhsMapping
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("InterfaceMatrices_RHSInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(interfaceMatrices)) CALL FlagError("Interface matrices is not associated.",err,error,*998)

    NULLIFY(interfaceMapping)
    CALL InterfaceMatrices_InterfaceMappingGet(interfaceMatrices,interfaceMapping,err,error,*999)
    NULLIFY(rhsMapping)
    CALL InterfaceMapping_RHSMappingGet(interfaceMapping,rhsMapping,err,error,*999)
    IF(ASSOCIATED(interfaceMatrices%rhsVector)) &
      & CALL FlagError("Interface matrices RHS vector is already associated.",err,error,*998)
    ALLOCATE(interfaceMatrices%rhsVector,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate interface matrices RHS vector.",err,error,*999)
    interfaceMatrices%rhsVector%rhsCoefficient=1.0_DP
    interfaceMatrices%rhsVector%updateVector=.TRUE.
    interfaceMatrices%rhsVector%firstAssembly=.TRUE.
    NULLIFY(interfaceMatrices%rhsVector%rhsVector)
    CALL EquationsMatrices_ElementVectorInitialise(interfaceMatrices%rhsVector%elementVector,err,error,*999)
    
    EXITS("InterfaceMatrices_RHSInitialise")
    RETURN
999 CALL InterfaceMatrices_RHSFinalise(interfaceMatrices%rhsVector,dummyErr,dummyError,*998)
998 ERRORSEXITS("InterfaceMatrices_RHSInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMatrices_RHSInitialise
  
  !
  !================================================================================================================================
  !

  !>Sets the matrix coefficients of the interface matrices
  SUBROUTINE InterfaceMatrices_MatrixCoefficientSet0(interfaceMatrices,matrixCoefficient,err,error,*)

    !Argument variables
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices !<A pointer to the interface matrices
    REAL(DP), INTENT(IN) :: matrixCoefficient !<The multiplicative matrix coefficient for the interface matrices. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("InterfaceMatrices_MatrixCoefficientSet0",err,error,*999)

    CALL InterfaceMatrices_MatrixCoefficientSet1(interfaceMatrices,[matrixCoefficient],err,error,*999)
     
    EXITS("InterfaceMatrices_MatrixCoefficientSet0")
    RETURN
999 ERRORSEXITS("InterfaceMatrices_MatrixCoefficientSet0",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMatrices_MatrixCoefficientSet0

  !
  !================================================================================================================================
  !

  !>Sets the matrix coefficients of the interface matrices
  SUBROUTINE InterfaceMatrices_MatrixCoefficientSet1(interfaceMatrices,matrixCoefficients,err,error,*)

    !Argument variables
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices !<A pointer to the interface matrices
    REAL(DP), INTENT(IN) :: matrixCoefficients(:) !<matrixCoefficients(matrixIdx). The multiplicative matrix coefficient for the matrixIdx'th inteface matrices.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(InterfaceMatrixType), POINTER :: interfaceMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("InterfaceMatrices_MatrixCoefficientSet1",err,error,*999)

    IF(.NOT.ASSOCIATED(interfaceMatrices)) CALL FlagError("Interface matrices is not associated.",err,error,*999)
    IF(SIZE(matrixCoefficients,1)/=interfaceMatrices%numberOfInterfaceMatrices) THEN
      localError="The size of the matrix coefficients array of "// &
        & TRIM(NumberToVString(SIZE(matrixCoefficients,1),"*",err,error))// &
        & " is not equal to the number of interface matrices of"// &
        & TRIM(NumberToVString(interfaceMatrices%numberOfInterfaceMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    DO matrixIdx=1,interfaceMatrices%numberOfInterfaceMatrices
      NULLIFY(interfaceMatrix)
      CALL InterfaceMatrices_InterfaceMatrixGet(interfaceMatrices,matrixIdx,interfaceMatrix,err,error,*999)
      CALL InterfaceMatrix_MatrixCoefficientSet(interfaceMatrix,matrixCoefficients(matrixIdx),err,error,*999)
    ENDDO !matrixIdx
     
    EXITS("InterfaceMatrices_MatrixCoefficientSet1")
    RETURN
999 ERRORSEXITS("InterfaceMatrices_MatrixCoefficientSet1",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMatrices_MatrixCoefficientSet1

  !
  !================================================================================================================================
  !

  !>Sets the storage type (sparsity) of the interface matrices
  SUBROUTINE InterfaceMatrices_StorageTypeSet0(interfaceMatrices,storageType,err,error,*)

    !Argument variables
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices !<A pointer to the interface matrices
    INTEGER(INTG), INTENT(IN) :: storageType !<storageType. The storage type for the inteface matrices. \see InterfaceMatricesRoutines_InterfaceMatricesSparsityTypes,InterfaceMatricesRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("InterfaceMatrices_StorageTypeSet0",err,error,*999)

    CALL InterfaceMatrices_StorageTypeSet1(interfaceMatrices,[storageType],err,error,*999)
     
    EXITS("InterfaceMatrices_StorageTypeSet0")
    RETURN
999 ERRORSEXITS("InterfaceMatrices_StorageTypeSet0",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMatrices_StorageTypeSet0

  !
  !================================================================================================================================
  !

  !>Sets the storage type (sparsity) of the interface matrices
  SUBROUTINE InterfaceMatrices_StorageTypeSet1(interfaceMatrices,storageTypes,err,error,*)

    !Argument variables
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices !<A pointer to the interface matrices
    INTEGER(INTG), INTENT(IN) :: storageTypes(:) !<storageTypes(matrixIdx). The storage type for the matrixIdx'th inteface matrices. \see InterfaceMatricesRoutines_InterfaceMatricesSparsityTypes,InterfaceMatricesRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(InterfaceMatrixType), POINTER :: interfaceMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("InterfaceMatrices_StorageTypeSet1",err,error,*999)

    CALL InterfaceMatrices_AssertNotFinished(interfaceMatrices,err,error,*999)

    IF(SIZE(storageTypes,1)/=interfaceMatrices%numberOfInterfaceMatrices) THEN
      localError="The size of the storage type array of "//TRIM(NumberToVString(SIZE(storageTypes,1),"*",err,error))// &
        & " is not equal to the number of interface matrices of"// &
        & TRIM(NumberToVString(interfaceMatrices%numberOfInterfaceMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    DO matrixIdx=1,interfaceMatrices%numberOfInterfaceMatrices
      NULLIFY(interfaceMatrix)
      CALL InterfaceMatrices_InterfaceMatrixGet(interfaceMatrices,matrixIdx,interfaceMatrix,err,error,*999)
      CALL InterfaceMatrix_StorageTypeSet(interfaceMatrix,storageTypes(matrixIdx),err,error,*999)
    ENDDO !matrixIdx
     
    EXITS("InterfaceMatrices_StorageTypeSet1")
    RETURN
999 ERRORSEXITS("InterfaceMatrices_StorageTypeSet1",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMatrices_StorageTypeSet1

  !
  !================================================================================================================================
  !

  !>Sets the structure (sparsity) of the interface matrices.
  SUBROUTINE InterfaceMatrices_StructureTypeSet0(interfaceMatrices,structureType,err,error,*)
    
    !Argument variables
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices !<A pointer to the interface matrices
    INTEGER(INTG), INTENT(IN) :: structureType !<The structure type for the interface matrix \see InterfaceMatricesRoutines_InterfaceMatrixStructureTypes,InterfaceMatricesRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("InterfaceMatrices_StructureTypeSet0",err,error,*999)

    CALL InterfaceMatrices_StructureTypeSet1(interfaceMatrices,[structureType],err,error,*999)
    
    EXITS("InterfaceMatrices_StructureTypeSet0")
    RETURN
999 ERRORSEXITS("InterfaceMatrices_StructureTypeSet0",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMatrices_StructureTypeSet0
  
  !
  !================================================================================================================================
  !

  !>Sets the structure (sparsity) of the interface matrices.
  SUBROUTINE InterfaceMatrices_StructureTypeSet1(interfaceMatrices,structureTypes,err,error,*)
    
    !Argument variables
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices !<A pointer to the interface matrices
    INTEGER(INTG), INTENT(IN) :: structureTypes(:) !<structureTypes(matrixIdx). The structure type for the  matrixIdx'th interface matrix \see InterfaceMatricesRoutines_InterfaceMatrixStructureTypes,InterfaceMatricesRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(InterfaceMatrixType), POINTER :: interfaceMatrix
    TYPE(VARYING_STRING) :: localError

    ENTERS("InterfaceMatrices_StructureTypeSet1",err,error,*999)

    CALL InterfaceMatrices_AssertNotFinished(interfaceMatrices,err,error,*999)
    IF(SIZE(structureTypes,1)/=interfaceMatrices%numberOfInterfaceMatrices) THEN
      localError="The size of the structure type array of "//TRIM(NumberToVString(SIZE(structureTypes,1),"*",err,error))// &
        & " is not equal to the number of interface matrices of "// &
        & TRIM(NumberToVString(interfaceMatrices%numberOfInterfaceMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    DO matrixIdx=1,interfaceMatrices%numberOfInterfaceMatrices
      NULLIFY(interfaceMatrix)
      CALL InterfaceMatrices_InterfaceMatrixGet(interfaceMatrices,matrixIdx,interfaceMatrix,err,error,*999)
      CALL InterfaceMatrix_StructureTypeSet(interfaceMatrix,structureTypes(matrixIdx),err,error,*999)
    ENDDO !matrixIdx
    
    EXITS("InterfaceMatrices_StructureTypeSet1")
    RETURN
999 ERRORSEXITS("InterfaceMatrices_StructureTypeSet1",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMatrices_StructureTypeSet1
  
  !
  !================================================================================================================================
  !

  !>Sets the transpose matrix coefficients of the interface matrices
  SUBROUTINE InterfaceMatrices_TransposeMatrixCoefficientSet0(interfaceMatrices,transposeMatrixCoefficient,err,error,*)

    !Argument variables
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices !<A pointer to the interface matrices
    REAL(DP), INTENT(IN) :: transposeMatrixCoefficient !<The multiplicative matrix coefficient for the transposed interface matrices. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("InterfaceMatrices_TransposeMatrixCoefficientSet0",err,error,*999)

    CALL InterfaceMatrices_TransposeMatrixCoefficientSet1(interfaceMatrices,[transposeMatrixCoefficient],err,error,*999)
     
    EXITS("InterfaceMatrices_TranposeMatrixCoefficientSet0")
    RETURN
999 ERRORSEXITS("InterfaceMatrices_TransposeMatrixCoefficientSet0",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMatrices_TransposeMatrixCoefficientSet0

  !
  !================================================================================================================================
  !

  !>Sets the matrix coefficients of the transposed interface matrices
  SUBROUTINE InterfaceMatrices_TransposeMatrixCoefficientSet1(interfaceMatrices,transposeMatrixCoefficients,err,error,*)

    !Argument variables
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices !<A pointer to the interface matrices
    REAL(DP), INTENT(IN) :: transposeMatrixCoefficients(:) !<transposeMatrixCoefficients(matrixIdx). The multiplicative matrix coefficient for the matrixIdx'th transposed inteface matrices.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(InterfaceMatrixType), POINTER :: interfaceMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("InterfaceMatrices_TransposeMatrixCoefficientSet1",err,error,*999)

    IF(.NOT.ASSOCIATED(interfaceMatrices)) CALL FlagError("Interface matrices is not associated.",err,error,*999)
    IF(SIZE(transposeMatrixCoefficients,1)/=interfaceMatrices%numberOfInterfaceMatrices) THEN
      localError="The size of the transpose matrix coefficients array of "// &
        & TRIM(NumberToVString(SIZE(transposeMatrixCoefficients,1),"*",err,error))// &
        & " is not equal to the number of interface matrices of"// &
        & TRIM(NumberToVString(interfaceMatrices%numberOfInterfaceMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    DO matrixIdx=1,interfaceMatrices%numberOfInterfaceMatrices
      NULLIFY(interfaceMatrix)
      CALL InterfaceMatrices_InterfaceMatrixGet(interfaceMatrices,matrixIdx,interfaceMatrix,err,error,*999)
      CALL InterfaceMatrix_TransposeMatrixCoefficientSet(interfaceMatrix,transposeMatrixCoefficients(matrixIdx),err,error,*999)
    ENDDO !matrixIdx
     
    EXITS("InterfaceMatrices_TransposeMatrixCoefficientSet1")
    RETURN
999 ERRORSEXITS("InterfaceMatrices_TransposeMatrixCoefficientSet1",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMatrices_TransposeMatrixCoefficientSet1

  !
  !================================================================================================================================
  !

  !>Initialise the values of the interface matrices to the given value e.g., 0.0_DP
  SUBROUTINE InterfaceMatrices_ValueInitialise(interfaceMatrices,value,err,error,*)
    
    !Argument variables
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices !<A pointer to the interface matrices to initialise the values for
    REAL(DP), INTENT(IN) :: value !<The value to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(InterfaceMatrixType), POINTER :: interfaceMatrix
    TYPE(InterfaceRHSType), POINTER :: rhsVector
    
    ENTERS("InterfaceMatrices_ValueInitialise",err,error,*999)
    
    IF(.NOT.ASSOCIATED(interfaceMatrices)) CALL FlagError("Interface matrices is not associated.",err,error,*999)
    
    DO matrixIdx=1,interfaceMatrices%numberOfInterfaceMatrices
      NULLIFY(interfaceMatrix)
      CALL InterfaceMatrices_InterfaceMatrixGet(interfaceMatrices,matrixIdx,interfaceMatrix,err,error,*999)
      IF(interfaceMatrix%updateMatrix) THEN
        CALL DistributedMatrix_AllValuesSet(interfaceMatrix%matrix,VALUE,err,error,*999)
        IF(interfaceMatrix%hasTranspose) CALL DistributedMatrix_AllValuesSet(interfaceMatrix%matrixTranspose,VALUE,err,error,*999)
      ENDIF
    ENDDO !matrixIdx
    rhsVector=>interfaceMatrices%rhsVector
    IF(ASSOCIATED(rhsVector)) THEN
      IF(rhsVector%updateVector) CALL DistributedVector_AllValuesSet(rhsVector%rhsVector,VALUE,err,error,*999)
    ENDIF

    EXITS("InterfaceMatrices_ValueInitialise")
    RETURN
999 ERRORSEXITS("InterfaceMatrices_ValueInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMatrices_ValueInitialise

  !
  !================================================================================================================================
  !

  !>Finalise a interface matrix and deallocate all memory
  SUBROUTINE InterfaceMatrix_Finalise(interfaceMatrix,err,error,*)

    !Argument variables
    TYPE(InterfaceMatrixType), POINTER :: interfaceMatrix !<A pointer to the interface matrix to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("InterfaceMatrix_Finalise",err,error,*999)

    IF(ASSOCIATED(interfaceMatrix)) THEN
      IF(ASSOCIATED(interfaceMatrix%matrix)) CALL DistributedMatrix_Destroy(interfaceMatrix%matrix,err,error,*999)
      IF(ASSOCIATED(interfaceMatrix%matrixTranspose)) CALL DistributedMatrix_Destroy(interfaceMatrix%matrixTranspose, &
        & err,error,*999)
      CALL EquationsMatrices_ElementMatrixFinalise(interfaceMatrix%elementMatrix,err,error,*999)
      DEALLOCATE(interfaceMatrix)
    ENDIF
    
    EXITS("InterfaceMatrix_Finalise")
    RETURN
999 ERRORSEXITS("InterfaceMatrix_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMatrix_Finalise

  !
  !================================================================================================================================
  !

  !>Initialise an interface matrix.
  SUBROUTINE InterfaceMatrix_Initialise(interfaceMatrices,matrixNumber,err,error,*)

    !Argument variables
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices !<A pointer to the interface matrices to initialise the interface matrix for
    INTEGER(INTG) :: matrixNumber !<The matrix number in the interface matrices to initialise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping
    TYPE(InterfaceMatrixType), POINTER :: interfaceMatrix
    TYPE(InterfaceMatrixToVarMapType), POINTER :: interfaceMatrixToVarMap
    TYPE(VARYING_STRING) :: dummyError,localError

    ENTERS("InterfaceMatrix_Initialise",err,error,*998)

    IF(.NOT.ASSOCIATED(interfaceMatrices)) CALL FlagError("Interface matrices is not associated.",err,error,*998)
    IF(matrixNumber<=0.OR.matrixNumber>interfaceMatrices%numberOfInterfaceMatrices) THEN
      localError="The specified interface matrix number of "//TRIM(NumberToVString(matrixNumber,"*",err,error))// &
        & " is invalid. The matrix number must be >= 1 and <= "// &
        & TRIM(NumberToVString(interfaceMatrices%numberOfInterfaceMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*998)
    ENDIF
    IF(ASSOCIATED(interfaceMatrices%matrices(matrixNumber)%ptr)) THEN
      localError="Interface matrix for matrix number "//TRIM(NumberToVString(matrixNumber,"*",err,error))// &
        & " is already associated."
      CALL FlagError(localError,err,error,*998)
    ENDIF

    NULLIFY(interfaceMapping)
    CALL InterfaceMatrices_InterfaceMappingGet(interfaceMatrices,interfaceMapping,err,error,*999)
    NULLIFY(interfaceMatrixToVarMap)
    CALL InterfaceMapping_InterfaceMatrixToVarMapGet(interfaceMapping,matrixNumber,interfaceMatrixToVarMap,err,error,*999)
    ALLOCATE(interfaceMatrices%matrices(matrixNumber)%ptr,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate interface matrix.",err,error,*999)
    interfaceMatrix=>interfaceMatrices%matrices(matrixNumber)%ptr
    interfaceMatrix%matrixNumber=matrixNumber
    interfaceMatrix%interfaceMatrices=>interfaceMatrices
    interfaceMatrix%matrixCoefficient=1.0_DP
    interfaceMatrix%transposeMatrixCoefficient=1.0_DP
    interfaceMatrix%storageType=MATRIX_BLOCK_STORAGE_TYPE
    interfaceMatrix%structureType=INTERFACE_MATRIX_NO_STRUCTURE
    interfaceMatrix%updateMatrix=.TRUE.
    interfaceMatrix%firstAssembly=.TRUE.
    CALL InterfaceMappingIMToVMap_HasTransposeGet(interfaceMatrixToVarMap,interfaceMatrix%hasTranspose,err,error,*999)
    CALL InterfaceMappingIMToVMap_NumberOfRowsGet(interfaceMatrixToVarMap,interfaceMatrix%numberOfRows,err,error,*999)
    CALL InterfaceMappingIMToVMap_TotalNumberOfRowsGet(interfaceMatrixToVarMap,interfaceMatrix%totalNumberOfRows,err,error,*999)
    interfaceMatrixToVarMap%interfaceMatrix=>interfaceMatrix
    NULLIFY(interfaceMatrix%matrix)
    NULLIFY(interfaceMatrix%matrixTranspose)
    NULLIFY(interfaceMatrix%tempVector)
    NULLIFY(interfaceMatrix%tempTransposeVector)
    CALL EquationsMatrices_ElementMatrixInitialise(interfaceMatrix%elementMatrix,err,error,*999)
   
    EXITS("InterfaceMatrix_Initialise")
    RETURN
999 CALL InterfaceMatrix_Finalise(interfaceMatrices%matrices(matrixNumber)%ptr,dummyErr,dummyError,*998)
998 ERRORSEXITS("InterfaceMatrix_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMatrix_Initialise

  !
  !================================================================================================================================
  !

  !>Sets/changes the matrix coefficient for an interface matrix.
  SUBROUTINE InterfaceMatrix_MatrixCoefficientSet(interfaceMatrix,matrixCoefficient,err,error,*)

    !Argument variables
    TYPE(InterfaceMatrixType), POINTER :: interfaceMatrix !<A pointer to the interface matrix to set the matrix coefficient for
    REAL(DP), INTENT(IN) :: matrixCoefficient !<The matrix coefficient to set for the interface matrix.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("InterfaceMatrix_MatrixCoefficientSet",err,error,*999)
    
    IF(.NOT.ASSOCIATED(interfaceMatrix)) CALL FlagError("Interface matrix is not associated.",err,error,*999)
    IF(ABS(matrixCoefficient)<=ZERO_TOLERANCE) THEN
      localError="The specified matrix coefficient of "//TRIM(NumberToVString(matrixCoefficient,"*",err,error))// &
        & " for interface matrix "//TRIM(NumberToVString(interfaceMatrix%matrixNumber,"*",err,error))// &
        & " is invalid. The matrix coefficient should be /= 0.0"
      CALL FlagError(localError,err,error,*999)
    ENDIF

    interfaceMatrix%matrixCoefficient=matrixCoefficient
       
    EXITS("InterfaceMatrix_MatrixCoefficientSet")
    RETURN
999 ERRORSEXITS("InterfaceMatrix_MatrixCoefficientSet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMatrix_MatrixCoefficientSet

  !
  !================================================================================================================================
  !

  !>Sets/changes the storage (sparsity) type for an interface matrix.
  SUBROUTINE InterfaceMatrix_StorageTypeSet(interfaceMatrix,storageType,err,error,*)

    !Argument variables
    TYPE(InterfaceMatrixType), POINTER :: interfaceMatrix !<A pointer to the interface matrix to set the storage type for
    INTEGER(INTG), INTENT(IN) :: storageType !<The storage type to set for the interface matrix. \see InterfaceMatricesRoutines_InterfaceMatricesSparsityTypes,InterfaceMatricesRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("InterfaceMatrix_StorageTypeSet",err,error,*999)
    
    IF(.NOT.ASSOCIATED(interfaceMatrix)) CALL FlagError("Interface matrix is not associated.",err,error,*999)

    SELECT CASE(storageType)
    CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
      interfaceMatrix%storageType=DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE
    CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
      interfaceMatrix%storageType=DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE        
    CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
      interfaceMatrix%storageType=DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE
    CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
      interfaceMatrix%storageType=DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE
    CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
      interfaceMatrix%storageType=DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE
    CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
      interfaceMatrix%storageType=DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE
    CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
      interfaceMatrix%storageType=DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE
    CASE DEFAULT
      localError="The specified storage type of "//TRIM(NumberToVString(storageType,"*",err,error))// &
        & " for interface matrix number "//TRIM(NumberToVString(interfaceMatrix%matrixNumber,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("InterfaceMatrix_StorageTypeSet")
    RETURN
999 ERRORSEXITS("InterfaceMatrix_StorageTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMatrix_StorageTypeSet

  !
  !================================================================================================================================
  !

  !>Sets/changes the structure (sparsity) type for an interface matrix.
  SUBROUTINE InterfaceMatrix_StructureTypeSet(interfaceMatrix,structureType,err,error,*)

    !Argument variables
    TYPE(InterfaceMatrixType), POINTER :: interfaceMatrix !<A pointer to the interface matrix to set the structure type for
    INTEGER(INTG), INTENT(IN) :: structureType !<The structure type to set for the interface matrix. \see InterfaceMatricesRoutines_InterfaceMatrixStructureTypes,InterfaceMatricesRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("InterfaceMatrix_StructureTypeSet",err,error,*999)
    
    IF(.NOT.ASSOCIATED(interfaceMatrix)) CALL FlagError("Interface matrix is not associated.",err,error,*999)
    
    SELECT CASE(structureType)
    CASE(INTERFACE_MATRIX_NO_STRUCTURE)
      interfaceMatrix%structureType=INTERFACE_MATRIX_NO_STRUCTURE
    CASE(INTERFACE_MATRIX_FEM_STRUCTURE)
      interfaceMatrix%structureType=INTERFACE_MATRIX_FEM_STRUCTURE
    CASE DEFAULT
      localError="The specified strucutre type of "// &
        & TRIM(NumberToVString(structureType,"*",err,error))//" for interface matrix number "// &
        & TRIM(NumberToVString(interfaceMatrix%matrixNumber,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("InterfaceMatrix_StructureTypeSet")
    RETURN
999 ERRORSEXITS("InterfaceMatrix_StructureTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMatrix_StructureTypeSet

  !
  !================================================================================================================================
  !

  !>Sets/changes the transpose matrix coefficient for an interface matrix.
  SUBROUTINE InterfaceMatrix_TransposeMatrixCoefficientSet(interfaceMatrix,transposeMatrixCoefficient,err,error,*)

    !Argument variables
    TYPE(InterfaceMatrixType), POINTER :: interfaceMatrix !<A pointer to the interface matrix to set the transpose matrix coefficient for
    REAL(DP), INTENT(IN) :: transposeMatrixCoefficient !<The transpose matrix coefficient to set for the interface matrix.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("InterfaceMatrix_TransposeMatrixCoefficientSet",err,error,*999)
    
    IF(.NOT.ASSOCIATED(interfaceMatrix)) CALL FlagError("Interface matrix is not associated.",err,error,*999)
    IF(ABS(transposeMatrixCoefficient)<=ZERO_TOLERANCE) THEN
      localError="The specified transpose matrix coefficient of "// &
        & TRIM(NumberToVString(transposeMatrixCoefficient,"*",err,error))// &
        & " for interface matrix "//TRIM(NumberToVString(interfaceMatrix%matrixNumber,"*",err,error))// &
        & " is invalid. The transpose matrix coefficient should be /= 0.0"
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.interfaceMatrix%hasTranspose) THEN
      localError="Can not set the matrix transpose coefficient for interface matrix "// &
        & TRIM(NumberToVString(interfaceMatrix%matrixNumber,"*",err,error))//" as it does not have a transpose."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    interfaceMatrix%transposeMatrixCoefficient=transposeMatrixCoefficient
       
    EXITS("InterfaceMatrix_TransposeMatrixCoefficientSet")
    RETURN
999 ERRORSEXITS("InterfaceMatrix_TranposeMatrixCoefficientSet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMatrix_TransposeMatrixCoefficientSet

  !
  !================================================================================================================================
  !

END MODULE InterfaceMatricesRoutines
