!> \file
!> \author Chris Bradley
!> \brief This module handles all solver matrix and rhs routines.
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

!> This module handles all solver matrix and rhs routines.
MODULE SolverMatricesRoutines

  USE BaseRoutines
  USE DistributedMatrixVector
  USE DistributedMatrixVectorAccessRoutines
  USE EquationsMatricesAccessRoutines
  USE InterfaceConditionAccessRoutines
  USE InterfaceMatricesAccessRoutines
  USE ISO_VARYING_STRING
  USE Kinds
  USE MatrixVector
  USE ProblemAccessRoutines
  USE SolverAccessRoutines
  USE SolverMappingAccessRoutines
  USE SolverMatricesAccessRoutines
  USE Strings
  USE Types

#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  INTERFACE SolverMatrices_StorageTypesSet
    MODULE PROCEDURE SolverMatrices_StorageTypesSet0
    MODULE PROCEDURE SolverMatrices_StorageTypesSet1
  END INTERFACE SolverMatrices_StorageTypesSet
  
  INTERFACE SolverMatrices_SymmetryTypesSet
    MODULE PROCEDURE SolverMatrices_SymmetryTypesSet0
    MODULE PROCEDURE SolverMatrices_SymmetryTypesSet1
  END INTERFACE SolverMatrices_SymmetryTypesSet
  
  !PUBLIC SolverMatrix_EquationsMatrixAdd

  !PUBLIC SolverMatrix_InterfaceMatrixAdd

  !PUBLIC SolverMatrix_JacobianMatrixAdd

  PUBLIC SolverMatrices_CreateFinish,SolverMatrices_CreateStart

  PUBLIC SolverMatrices_Destroy

  PUBLIC SolverMatrices_LibraryTypeGet,SolverMatrices_LibraryTypeSet

  PUBLIC SolverMatrices_Output

  PUBLIC SolverMatrices_StorageTypesSet

  PUBLIC SolverMatrices_SymmetryTypesSet

CONTAINS

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating the solver matrices
  SUBROUTINE SolverMatrices_CreateFinish(solverMatrices,err,error,*)

    !Argument variables
    TYPE(SolverMatricesType), POINTER :: solverMatrices !<A pointer to the solver matrices
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,matrixIdx,numberOfNonZeros
    INTEGER(INTG), POINTER :: columnIndices(:),rowIndices(:)
    TYPE(DomainMappingType), POINTER :: rowDomainMap,columnDomainMap
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolverMatrixType), POINTER :: solverMatrix
    TYPE(VARYING_STRING) :: dummyError

    NULLIFY(columnIndices)
    NULLIFY(rowIndices)
    
    ENTERS("SolverMatrices_CreateFinish",err,error,*998)

    CALL SolverMatrices_AssertNotFinished(solverMatrices,err,error,*999)
    NULLIFY(solverEquations)
    CALL SolverMatrices_SolverEquationsGet(solverMatrices,solverEquations,err,error,*999)
    NULLIFY(solverMapping)
    CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
    !Now create the individual solver matrices
    NULLIFY(rowDomainMap)
    CALL SolverMapping_RowDOFSMappingGet(solverMapping,rowDomainMap,err,error,*999)
    DO matrixIdx=1,solverMatrices%numberOfMatrices
      NULLIFY(solverMatrix)
      CALL SolverMatrices_SolverMatrixGet(solverMatrices,matrixIdx,solverMatrix,err,error,*999)
      NULLIFY(columnDomainMap)
      CALL SolverMapping_ColumnDOFSMappingGet(solverMapping,matrixIdx,columnDomainMap,err,error,*999)
      !Create the distributed solver matrix
      CALL DistributedMatrix_CreateStart(rowDomainMap,columnDomainMap,solverMatrices%matrices(matrixIdx)%ptr%matrix,err,error,*999)
      CALL DistributedMatrix_LibraryTypeSet(solverMatrix%matrix,solverMatrices%matrixLibraryType,err,error,*999)
      CALL DistributedMatrix_DataTypeSet(solverMatrix%matrix,MATRIX_VECTOR_DP_TYPE,err,error,*999)
      CALL DistributedMatrix_StorageTypeSet(solverMatrix%matrix,solverMatrix%storageType,err,error,*999)
      CALL DistributedMatrix_TransposeTypeSet(solverMatrix%matrix,DISTRIBUTED_MATRIX_NO_TRANSPOSE_REQUIRED,err,error,*999)
      !Calculate and set the matrix structure/sparsity pattern
      IF(solverMatrix%storageType/=DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE.AND. &
        & solverMatrix%storageType/=DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE) THEN
        CALL SolverMatrix_StructureCalculate(solverMatrix,numberOfNonZeros,rowIndices,columnIndices,err,error,*999) 
        CALL DistributedMatrix_NumberOfNonZerosSet(solverMatrix%matrix,numberOfNonZeros,err,error,*999)
        CALL DistributedMatrix_StorageLocationsSet(solverMatrix%matrix,rowIndices,columnIndices,err,error,*999)
        IF(ASSOCIATED(rowIndices)) DEALLOCATE(rowIndices)
        IF(ASSOCIATED(columnIndices)) DEALLOCATE(columnIndices)
      ENDIF
      CALL DistributedMatrix_SymmetryTypeSet(solverMatrix%matrix,solverMatrix%symmetryType,err,error,*999)
      CALL DistributedMatrix_CreateFinish(solverMatrix%matrix,err,error,*999)
      !Allocate the distributed solver vector
      CALL DistributedVector_CreateStart(columnDomainMap,solverMatrices%matrices(matrixIdx)%ptr%solverVector,err,error,*999)
      CALL DistributedVector_LibraryTypeSet(solverMatrix%solverVector,solverMatrices%matrixLibraryType,err,error,*999)
      CALL DistributedVector_DataTypeSet(solverMatrix%solverVector,MATRIX_VECTOR_DP_TYPE,err,error,*999)
      CALL DistributedVector_CreateFinish(solverMatrix%solverVector,err,error,*999)
    ENDDO !matrixIdx
    IF(solverEquations%linearity==PROBLEM_SOLVER_NONLINEAR) THEN
      !Allocate the nonlinear matrices and vectors                  
      !Allocate the distributed residual vector
      CALL DistributedVector_CreateStart(rowDomainMap,solverMatrices%residual,err,error,*999)
      CALL DistributedVector_LibraryTypeSet(solverMatrices%residual,solverMatrices%matrixLibraryType,err,error,*999)
      CALL DistributedVector_DataTypeSet(solverMatrices%residual,MATRIX_VECTOR_DP_TYPE,err,error,*999)
      CALL DistributedVector_CreateFinish(solverMatrices%residual,err,error,*999)                  
    ENDIF
!!TODO: what to do if there is no RHS
    !Allocate the distributed rhs vector
    CALL DistributedVector_CreateStart(rowDomainMap,solverMatrices%rhsVector,err,error,*999)
    CALL DistributedVector_LibraryTypeSet(solverMatrices%rhsVector,solverMatrices%matrixLibraryType,err,error,*999)
    CALL DistributedVector_DataTypeSet(solverMatrices%rhsVector,MATRIX_VECTOR_DP_TYPE,err,error,*999)
    CALL DistributedVector_CreateFinish(solverMatrices%rhsVector,err,error,*999)
    !Finish up
    solverMatrices%solverMatricesFinished=.TRUE.
    
    EXITS("SolverMatrices_CreateFinish")
    RETURN
999 IF(ASSOCIATED(rowIndices)) DEALLOCATE(rowIndices)
    IF(ASSOCIATED(columnIndices)) DEALLOCATE(columnIndices)
    CALL SolverMatrices_Finalise(solverMatrices,dummyErr,dummyError,*998)
998 ERRORSEXITS("SolverMatrices_CreateFinish",err,error)    
    RETURN 1
    
  END SUBROUTINE SolverMatrices_CreateFinish
        
  !
  !================================================================================================================================
  !

  !>Starts the process of creating the solver matrices
  SUBROUTINE SolverMatrices_CreateStart(solverEquations,solverMatrices,err,error,*)

    !Argument variables
    TYPE(SolverEquationsType), POINTER :: solverEquations !<A pointer to the solver equations to create the solver matrices for
    TYPE(SolverMatricesType), POINTER :: solverMatrices !<On return, a pointer to the solver matrices. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("SolverMatrices_CreateStart",err,error,*998)

    IF(ASSOCIATED(solverMatrices)) CALL FlagError("Solver matrices is already associated.",err,error,*998)
    CALL SolverEquations_AssertIsFinished(solverEquations,err,error,*999)
    IF(ASSOCIATED(solverEquations%solverMatrices)) &
      & CALL FlagError("The solver equations already have solver matrices associated.",err,error,*999)
    
    CALL SolverMatrices_Initialise(solverEquations,err,error,*999)
    solverMatrices=>solverEquations%solverMatrices
        
    EXITS("SolverMatrices_CreateStart")
    RETURN
999 CALL SolverMatrices_Finalise(solverEquations%solverMatrices,dummyErr,dummyError,*998)
998 ERRORSEXITS("SolverMatrices_CreateStart",err,error)    
    RETURN 1
   
  END SUBROUTINE SolverMatrices_CreateStart
        
  !
  !================================================================================================================================
  !

  !>Destroy the solver matrices
  SUBROUTINE SolverMatrices_Destroy(solverMatrices,err,error,*)

    !Argument variables
    TYPE(SolverMatricesType), POINTER :: solverMatrices !<A pointer the solver matrices to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverMatrices_Destroy",err,error,*999)

    IF(.NOT.ASSOCIATED(solverMatrices)) CALL FlagError("Solver matrices is not associated",err,error,*999)
    
    CALL SolverMatrices_Finalise(solverMatrices,err,error,*999)
        
    EXITS("SolverMatrices_Destroy")
    RETURN
999 ERRORSEXITS("SolverMatrices_Destroy",err,error)    
    RETURN 1
   
  END SUBROUTINE SolverMatrices_Destroy

  !
  !================================================================================================================================
  !

  !>Finalises the solver matrices and deallocates all memory
  SUBROUTINE SolverMatrices_Finalise(solverMatrices,err,error,*)

    !Argument variables
    TYPE(SolverMatricesType), POINTER :: solverMatrices !<A pointer to the solver matrices to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx

    ENTERS("SolverMatrices_Finalise",err,error,*999)

    IF(ASSOCIATED(solverMatrices)) THEN
      IF(ALLOCATED(solverMatrices%matrices)) THEN
        DO matrixIdx=1,SIZE(solverMatrices%matrices,1)
          CALL SolverMatrix_Finalise(solverMatrices%matrices(matrixIdx)%ptr,err,error,*999)
        ENDDO !matrixIdx
        DEALLOCATE(solverMatrices%matrices)
      ENDIF
      IF(ASSOCIATED(solverMatrices%residual)) CALL DistributedVector_Destroy(solverMatrices%residual,err,error,*999)
      IF(ASSOCIATED(solverMatrices%rhsVector)) CALL DistributedVector_Destroy(solverMatrices%rhsVector,err,error,*999)
      DEALLOCATE(solverMatrices)
    ENDIF
        
    EXITS("SolverMatrices_Finalise")
    RETURN
999 ERRORSEXITS("SolverMatrices_Finalise",err,error)    
    RETURN 1
   
  END SUBROUTINE SolverMatrices_Finalise
        
  !
  !================================================================================================================================
  !

  !>Initialises the solver matrices for solver equations
  SUBROUTINE SolverMatrices_Initialise(solverEquations,err,error,*)

    !Argument variables
    TYPE(SolverEquationsType), POINTER :: solverEquations !<A pointer to the solver equations to initialise the solver matrices for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,equationsMatrixIdx,equationsSetIdx,numberOfDynamicMatrices,numberOfEquationsSets, &
      & numberOfGlobalRows,numberOfJacobianMatrices,numberOfLinearMatrices,numberOfRows,numberOfSolverMatrices, &
      & solverEquationsLinearity,solverMatrixIdx
    TYPE(EquationsMatricesToSolverMatrixMapType), POINTER :: equationsMatricesToSolverMatrixMap
    TYPE(EquationsMatrixToSolverMatrixMapType), POINTER :: dynamicEquationsMatrixToSolverMatrixMap, &
      & linearEquationsMatrixToSolverMatrixMap
    TYPE(EquationsSetToSolverMatricesMapType), POINTER :: equationsSetToSolverMatricesMap
    TYPE(JacobianMatrixToSolverMatrixMapType), POINTER :: jacobianMatrixToSolverMatrixMap
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("SolverMatrices_Initialise",err,error,*998)

    IF(.NOT.ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is not associated",err,error,*998)
    IF(ASSOCIATED(solverEquations%solverMatrices)) &
      & CALL FlagError("Solver matrices is already associated for this solver equations.",err,error,*998)

    NULLIFY(solverMapping)
    CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
    CALL SolverMapping_NumberOfSolverMatricesGet(solverMapping,numberOfSolverMatrices,err,error,*999)
    CALL SolverMapping_NumberOfRowsGet(solverMapping,numberOfRows,err,error,*999)
    CALL SolverMapping_NumberOfGlobalRowsGet(solverMapping,numberOfGlobalRows,err,error,*999)
    CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
    
    ALLOCATE(solverEquations%solverMatrices,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate solver matrices.",err,error,*999)
    solverEquations%solverMatrices%solverEquations=>solverEquations
    solverEquations%solverMatrices%solverMatricesFinished=.FALSE.
    solverEquations%solverMatrices%solverMapping=>solverMapping
    solverEquations%solverMatrices%numberOfRows=numberOfRows
    solverEquations%solverMatrices%numberOfGlobalRows=numberOfGlobalRows
    solverEquations%solverMatrices%solverLibraryType=0
    solverEquations%solverMatrices%matrixLibraryType=0
    solverEquations%solverMatrices%numberOfMatrices=numberOfSolverMatrices
    ALLOCATE(solverEquations%solverMatrices%matrices(numberOfSolverMatrices),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate solver matrices matrices.",err,error,*999)
    DO solverMatrixIdx=1,numberOfSolverMatrices
      NULLIFY(solverEquations%solverMatrices%matrices(solverMatrixIdx)%ptr)
      CALL SolverMatrix_Initialise(solverEquations%solverMatrices,solverMatrixIdx,err,error,*999)
      DO equationsSetIdx=1,numberOfEquationsSets
        NULLIFY(equationsSetToSolverMatricesMap)
        CALL SolverMapping_EquationsSetToSolverMatricesMapGet(solverMapping,equationsSetIdx,equationsSetToSolverMatricesMap, &
          & err,error,*999)
        NULLIFY(equationsMatricesToSolverMatrixMap)
        CALL SolverMappingESToSMSMap_EquationsMatricesToSolverMatrixMapGet(equationsSetToSolverMatricesMap,solverMatrixIdx, &
          & equationsMatricesToSolverMatrixMap,err,error,*999)
        CALL SolverMappingEMSToSMMap_NumberOfDynamicMatricesGet(equationsMatricesToSolverMatrixMap,numberOfDynamicMatrices, &
          & err,error,*999)
        IF(numberOfDynamicMatrices>0) THEN
          DO equationsMatrixIdx=1,numberOfDynamicMatrices
            !Add the solver matrix to the solvers mapping
            NULLIFY(dynamicEquationsMatrixToSolverMatrixMap)
            CALL SolverMappingEMSToSMMap_DynamicMatrixToSolverMatrixMapGet(equationsMatricesToSolverMatrixMap, &
              & equationsMatrixIdx,dynamicEquationsMatrixToSolverMatrixMap,err,error,*999)
            dynamicEquationsMatrixToSolverMatrixMap%solverMatrix=>solverEquations%solverMatrices%matrices(solverMatrixIdx)%ptr
          ENDDO !equationsMatrixIdx
        ELSE
          CALL SolverMappingEMSToSMMap_NumberOfJacobianMatricesGet(equationsMatricesToSolverMatrixMap,numberOfJacobianMatrices, &
            & err,error,*999)
          IF(numberOfJacobianMatrices>0) THEN
            DO equationsMatrixIdx=1,numberOfJacobianMatrices
              NULLIFY(jacobianMatrixToSolverMatrixMap)
              CALL SolverMappingEMSToSMMap_JacobianMatrixToSolverMatrixMapGet(equationsMatricesToSolverMatrixMap, &
                & equationsMatrixIdx,jacobianMatrixToSolverMatrixMap,err,error,*999)
              jacobianMatrixToSolverMatrixMap%solverMatrix=>solverEquations%solverMatrices%matrices(solverMatrixIdx)%ptr
            ENDDO !equationsMatrixIdx
          ELSE
            CALL SolverMappingEMSToSMMap_NumberOfLinearMatricesGet(equationsMatricesToSolverMatrixMap,numberOfLinearMatrices, &
            & err,error,*999)
            DO equationsMatrixIdx=1,numberOfLinearMatrices
              !Add the solver matrix to the solvers mapping
              NULLIFY(linearEquationsMatrixToSolverMatrixMap)
              CALL SolverMappingEMSToSMMap_LinearMatrixToSolverMatrixMapGet(equationsMatricesToSolverMatrixMap, &
              & equationsMatrixIdx,linearEquationsMatrixToSolverMatrixMap,err,error,*999)
              linearEquationsMatrixToSolverMatrixMap%solverMatrix=>solverEquations%solverMatrices%matrices(solverMatrixIdx)%ptr
            ENDDO !equationsMatrixIdx
          ENDIF
        ENDIF
      ENDDO !equationsSetIdx
    ENDDO !solverMatrixIdx
    CALL SolverEquations_LinearityTypeGet(solverEquations,solverEquationsLinearity,err,error,*999)
    IF(solverEquationsLinearity==PROBLEM_SOLVER_NONLINEAR) THEN
      solverEquations%solverMatrices%updateResidual=.TRUE.
    ELSE
      solverEquations%solverMatrices%updateResidual=.FALSE.
    ENDIF
    NULLIFY(solverEquations%solverMatrices%residual)
    solverEquations%solverMatrices%updateRHSVector=.TRUE.
    NULLIFY(solverEquations%solverMatrices%rhsVector)
       
    EXITS("SolverMatrices_Initialise")
    RETURN
999 CALL SolverMatrices_Finalise(solverEquations%solverMatrices,dummyErr,dummyError,*998)
998 ERRORSEXITS("SolverMatrices_Initialise",err,error)    
    RETURN 1
   
  END SUBROUTINE SolverMatrices_Initialise
  
  !
  !================================================================================================================================
  !

  !>Sets the library type for the solver matrices (and vectors)
  SUBROUTINE SolverMatrices_LibraryTypeSet(solverMatrices,libraryType,err,error,*)

    !Argument variables
    TYPE(SolverMatricesType), POINTER :: solverMatrices !<A pointer to the solver matrices.
    INTEGER(INTG), INTENT(IN) :: libraryType !<The library type to set \see SolverRoutines_SolverLibraries
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverMatrices_LibraryTypeSet",err,error,*999)

    CALL SolverMatrices_AssertNotFinished(solverMatrices,err,error,*999)
    
    SELECT CASE(libraryType)
    CASE(LIBRARY_CMISS_TYPE)
      solverMatrices%matrixLibraryType=DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE
    CASE(LIBRARY_PETSC_TYPE)
      solverMatrices%matrixLibraryType=DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE
    CASE(LIBRARY_MUMPS_TYPE)
      solverMatrices%matrixLibraryType=DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE
    CASE(LIBRARY_SUPERLU_TYPE)
      solverMatrices%matrixLibraryType=DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE
    CASE(LIBRARY_SPOOLES_TYPE)
    CASE(LIBRARY_UMFPACK_TYPE)
      solverMatrices%matrixLibraryType=DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE
    CASE(LIBRARY_LUSOL_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(LIBRARY_ESSL_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(LIBRARY_LAPACK_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(LIBRARY_HYPRE_TYPE)
      solverMatrices%matrixLibraryType=DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE
    CASE(LIBRARY_PASTIX_TYPE)
      solverMatrices%matrixLibraryType=DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE
    CASE DEFAULT
      localError="The solver library type of "// TRIM(NumberToVString(libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    solverMatrices%solverLibraryType=libraryType
    
    EXITS("SolverMatrices_LibraryTypeSet")
    RETURN
999 ERRORSEXITS("SolverMatrices_LibraryTypeSet",err,error)
    RETURN 1
     
  END SUBROUTINE SolverMatrices_LibraryTypeSet

  !
  !================================================================================================================================
  !

  !>Outputs the solver matrices
  SUBROUTINE SolverMatrices_Output(id,selectionType,solverMatrices,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The id of the ouptut stream
    INTEGER(INTG), INTENT(IN) :: selectionType !<The type of matrix selection \see SolverMatricesRoutines_SelectMatricesTypes,SolverMatricesRoutines
    TYPE(SolverMatricesType), POINTER :: solverMatrices !<A pointer to the solver matrices
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(SolverMatrixType), POINTER :: solverMatrix
    
    ENTERS("SolverMatrices_Output",err,error,*999)

    CALL SolverMatrices_AssertIsFinished(solverMatrices,err,error,*999)
    CALL WriteString(id,"",err,error,*999)
    IF(selectionType==SOLVER_MATRICES_ALL.OR. &
      !& selectionType==SOLVER_MATRICES_DYNAMIC_ONLY.OR. &
      & selectionType==SOLVER_MATRICES_LINEAR_ONLY.OR. &
      & selectionType==SOLVER_MATRICES_NONLINEAR_ONLY.OR. &
      & selectionType==SOLVER_MATRICES_JACOBIAN_ONLY) THEN
      CALL WriteString(id,"Solver matrices:",err,error,*999)
      CALL WriteStringValue(id,"Number of matrices = ",solverMatrices%numberOfMatrices,err,error,*999)
      DO matrixIdx=1,solverMatrices%numberOfMatrices
        NULLIFY(solverMatrix)
        CALL SolverMatrices_SolverMatrixGet(solverMatrices,matrixIdx,solverMatrix,err,error,*999)
        CALL WriteStringValue(id,"Solver matrix : ",matrixIdx,err,error,*999)
        CALL DistributedMatrix_Output(id,solverMatrix%matrix,err,error,*999)
      ENDDO !matrixIdx
    ENDIF
    IF(selectionType==SOLVER_MATRICES_ALL.OR. &
      & selectionType==SOLVER_MATRICES_NONLINEAR_ONLY.OR. &
      & selectionType==SOLVER_MATRICES_RESIDUAL_ONLY.OR. &
      & selectionType==SOLVER_MATRICES_RHS_RESIDUAL_ONLY) THEN
      IF(ASSOCIATED(solverMatrices%residual)) THEN
        CALL WriteString(id,"Solver residual vector:",err,error,*999)     
        CALL DistributedVector_Output(id,solverMatrices%residual,err,error,*999)  
      ENDIF
    ENDIF
    IF(selectionType==SOLVER_MATRICES_ALL.OR. &
      !& selectionType==SOLVER_MATRICES_DYNAMIC_ONLY.OR. &
      & selectionType==SOLVER_MATRICES_LINEAR_ONLY.OR. &
      & selectionType==SOLVER_MATRICES_NONLINEAR_ONLY.OR. &
      & selectionType==SOLVER_MATRICES_RHS_ONLY.OR. &
      & selectionType==SOLVER_MATRICES_RHS_RESIDUAL_ONLY) THEN
      IF(ASSOCIATED(solverMatrices%rhsVector)) THEN
        CALL WriteString(id,"Solver RHS vector:",err,error,*999)     
        CALL DistributedVector_Output(id,solverMatrices%rhsVector,err,error,*999)
      ENDIF
    ENDIF
    
    EXITS("SolverMatrices_Output")
    RETURN
999 ERRORSEXITS("SolverMatrices_Output",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMatrices_Output
  
  !
  !================================================================================================================================
  !

  !>Sets the storage type (sparsity) of the solver matrices
  SUBROUTINE SolverMatrices_StorageTypesSet0(solverMatrices,storageType,err,error,*)

    !Argument variables
    TYPE(SolverMatricesType), POINTER :: solverMatrices !<A pointer to the solver matrices
    INTEGER(INTG), INTENT(IN) :: storageType !<The storage type for the solver matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(SolverMatrixType), POINTER :: solverMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("SolverMatrices_StorageTypesSet0",err,error,*999)

    CALL SolverMatrices_StorageTypesSet1(solverMatrices,[storageType],err,error,*999)
   
    EXITS("SolverMatrices_StorageTypesSet0")
    RETURN
999 ERRORSEXITS("SolverMatrices_StorageTypesSet0",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMatrices_StorageTypesSet0

  !
  !================================================================================================================================
  !

  !>Sets the storage type (sparsity) of the solver matrices
  SUBROUTINE SolverMatrices_StorageTypesSet1(solverMatrices,storageTypes,err,error,*)

    !Argument variables
    TYPE(SolverMatricesType), POINTER :: solverMatrices !<A pointer to the solver matrices
    INTEGER(INTG), INTENT(IN) :: storageTypes(:) !<storageTypes(matrixIdx). The storage type for the matrixIdx'th solver matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(SolverMatrixType), POINTER :: solverMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("SolverMatrices_StorageTypesSet1",err,error,*999)

    CALL SolverMatrices_AssertNotFinished(solverMatrices,err,error,*999)
    IF(SIZE(storageTypes,1)<solverMatrices%numberOfMatrices) THEN
      localError="The size of the solver storage type array of "//TRIM(NumberToVString(SIZE(storageTypes,1),"*",err,error))// &
        & " is less than the number of solver matrices of "// &
        & TRIM(NumberToVString(solverMatrices%numberOfMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    DO matrixIdx=1,solverMatrices%numberOfMatrices
      NULLIFY(solverMatrix)
      CALL SolverMatrices_SolverMatrixGet(solverMatrices,matrixIdx,solverMatrix,err,error,*999)
      SELECT CASE(storageTypes(matrixIdx))
      CASE(MATRIX_BLOCK_STORAGE_TYPE)
        solverMatrix%storageType=MATRIX_BLOCK_STORAGE_TYPE
      CASE(MATRIX_DIAGONAL_STORAGE_TYPE)
        solverMatrix%storageType=MATRIX_DIAGONAL_STORAGE_TYPE        
      CASE(MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
        solverMatrix%storageType=MATRIX_COLUMN_MAJOR_STORAGE_TYPE
      CASE(MATRIX_ROW_MAJOR_STORAGE_TYPE)
        solverMatrix%storageType=MATRIX_ROW_MAJOR_STORAGE_TYPE
      CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
        solverMatrix%storageType=MATRIX_COMPRESSED_ROW_STORAGE_TYPE
      CASE(MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
        solverMatrix%storageType=MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE
      CASE(MATRIX_ROW_COLUMN_STORAGE_TYPE)
        solverMatrix%storageType=MATRIX_ROW_COLUMN_STORAGE_TYPE
      CASE DEFAULT
        localError="The specified storage type of "//TRIM(NumberToVString(storageTypes(matrixIdx),"*",err,error))// &
          & " for the matrix number "//TRIM(NumberToVString(matrixIdx,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDDO !matrixIdx
   
    EXITS("SolverMatrices_StorageTypesSet1")
    RETURN
999 ERRORSEXITS("SolverMatrices_StorageTypesSet1",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMatrices_StorageTypesSet1

  !
  !================================================================================================================================
  !

  !>Sets the symmetry types of the solver matrices
  SUBROUTINE SolverMatrices_SymmetryTypesSet0(solverMatrices,symmetryType,err,error,*)

    !Argument variables
    TYPE(SolverMatricesType), POINTER :: solverMatrices !<A pointer to the solver matrices to set the symmetry types for
    INTEGER(INTG), INTENT(IN) :: symmetryType !<The symmetry type for the solver matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("SolverMatrices_SymmetryTypesSet0",err,error,*999)

    CALL SolverMatrices_SymmetryTypesSet1(solverMatrices,[symmetryType],err,error,*999)
    
    EXITS("SolverMatrices_SymmetryTypesSet0")
    RETURN
999 ERRORSEXITS("SolverMatrices_SymmetryTypesSet0",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMatrices_SymmetryTypesSet0

  !
  !================================================================================================================================
  !

  !>Sets the symmetry type of the solver matrices
  SUBROUTINE SolverMatrices_SymmetryTypesSet1(solverMatrices,symmetryTypes,err,error,*)

    !Argument variables
    TYPE(SolverMatricesType), POINTER :: solverMatrices !<A pointer to the solver matrices to set the symmetry types for
    INTEGER(INTG), INTENT(IN) :: symmetryTypes(:) !<symmetryTypes(matrixIdx). The symmetry type for the matrixIdx'th solver matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(SolverMatrixType), POINTER :: solverMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("SolverMatrices_SymmetryTypesSet1",err,error,*999)

    CALL SolverMatrices_AssertNotFinished(solverMatrices,err,error,*999)
    IF(SIZE(symmetryTypes,1)/=solverMatrices%numberOfMatrices) THEN
      localError="The size of the symmetry types array of "//TRIM(NumberToVString(SIZE(symmetryTypes,1),"*",err,error))// &
        & " is not equal to the number of matrices of "//TRIM(NumberToVString(solverMatrices%numberOfMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
      
    DO matrixIdx=1,solverMatrices%numberOfMatrices
      NULLIFY(solverMatrix)
      CALL SolverMatrices_SolverMatrixGet(solverMatrices,matrixIdx,solverMatrix,err,error,*999)
      solverMatrix%symmetryType=symmetryTypes(matrixIdx)
    ENDDO !matrixIdx
    
    EXITS("SolverMatrices_SymmetryTypesSet1")
    RETURN
999 ERRORSEXITS("SolverMatrices_SymmetryTypesSet1",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMatrices_SymmetryTypesSet1

  !
  !================================================================================================================================
  !

  !>Adds alpha times the equations matrix into the solver matrix
  SUBROUTINE SolverMatrix_EquationsMatrixAdd(solverMatrix,equationsSetIdx,alpha,equationsMatrix,err,error,*)

    !Argument variables
    TYPE(SolverMatrixType), POINTER :: solverMatrix !<A pointer to the solver matrix
    INTEGER(INTG), INTENT(IN) :: equationsSetIdx !<The equations set index in the solver mapping that contains the equations matrix to add
    REAL(DP), INTENT(IN) :: alpha !<The multiplicative factor for the equations matrix
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix !<A pointer to the equations matrix to add    
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsMatrixNumber,solverMatrixNumber
    TYPE(DistributedMatrixType), POINTER :: equationsDistributedMatrix,solverDistributedMatrix
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatrixToSolverMatricesMapType), POINTER :: equationsMatrixToSolverMatricesMap
    TYPE(EquationsMatrixToSolverMatrixMapType), POINTER :: equationsMatrixToSolverMatrixMap
    TYPE(EquationsSetToSolverMatricesMapType), POINTER :: equationsSetToSolverMatricesMap
    TYPE(MatrixRowColCouplingType), POINTER :: equationsColToSolverColsMap(:),equationsRowToSolverRowsMap(:)
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolverMatricesType), POINTER :: solverMatrices
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("SolverMatrix_EquationsMatrixAdd",err,error,*999)

    IF(ABS(alpha)>ZERO_TOLERANCE) THEN
      NULLIFY(solverMatrices)
      CALL SolverMatrix_SolverMatricesGet(solverMatrix,solverMatrices,err,error,*999)
      CALL SolverMatrices_AssertIsFinished(solverMatrices,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverMatrices_SolverMappingGet(solverMatrices,solverMapping,err,error,*999)
      NULLIFY(linearMatrices)
      CALL EquationsMatrix_LinearMatricesExists(equationsMatrix,linearMatrices,err,error,*999)
      NULLIFY(dynamicMatrices)
      CALL EquationsMatrix_DynamicMatricesExists(equationsMatrix,dynamicMatrices,err,error,*999)
      IF(.NOT.ASSOCIATED(dynamicMatrices).AND..NOT.ASSOCIATED(linearMatrices)) &
        & CALL FlagError("Equations matrix dynamic or linear matrices is not associated.",err,error,*999)
      NULLIFY(vectorMatrices)
      IF(ASSOCIATED(dynamicMatrices)) THEN
        CALL EquationsMatricesDynamic_VectorMatricesGet(dynamicMatrices,vectorMatrices,err,error,*999)
      ELSE
        CALL EquationsMatricesLinear_VectorMatricesGet(linearMatrices,vectorMatrices,err,error,*999)
      ENDIF
      CALL EquationsMatricesVector_AssertIsFinished(vectorMatrices,err,error,*999)
      NULLIFY(solverDistributedMatrix)
      CALL SolverMatrix_SolverDistributedMatrixGet(solverMatrix,solverDistributedMatrix,err,error,*999)
      NULLIFY(equationsDistributedMatrix)
      CALL EquationsMatrix_DistributedMatrixGet(equationsMatrix,equationsDistributedMatrix,err,error,*999)
      NULLIFY(equationsSetToSolverMatricesMap)
      CALL SolverMapping_EquationsSetToSolverMatricesMapGet(solverMapping,equationsSetIdx,equationsSetToSolverMatricesMap, &
        & err,error,*999)
      NULLIFY(equationsRowToSolverRowsMap)
      CALL SolverMappingESToSMSMap_EquationsRowToSolverRowsMapGet(equationsSetToSolverMatricesMap,equationsRowToSolverRowsMap, &
        & err,error,*999)
      CALL EquationsMatrix_MatrixNumberGet(equationsMatrix,equationsMatrixNumber,err,error,*999)
      NULLIFY(equationsMatrixToSolverMatricesMap)
      CALL SolverMappingESToSMSMap_EquationsMatrixToSolverMatricesMapGet(equationsSetToSolverMatricesMap,equationsMatrixNumber, &
        & equationsMatrixToSolverMatricesMap,err,error,*999)
      CALL SolverMatrix_MatrixNumberGet(solverMatrix,solverMatrixNumber,err,error,*999)
      NULLIFY(equationsMatrixToSolverMatrixMap)
      CALL SolverMappingEMToSMSMap_EquationsMatrixToSolverMatrixMapGet(equationsMatrixToSolverMatricesMap,solverMatrixNumber, &
        & equationsMatrixToSolverMatrixMap,err,error,*999)
      NULLIFY(equationsColToSolverColsMap)
      CALL SolverMappingEMToSMMap_EquationsColToSolverColsMapGet(equationsMatrixToSolverMatrixMap,equationsColToSolverColsMap, &
        & err,error,*999)
      CALL DistributedMatrix_MatrixCoupleAdd(solverDistributedMatrix,DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE, &
        & equationsRowToSolverRowsMap,equationsColToSolverColsMap,alpha,equationsDistributedMatrix,.FALSE.,err,error,*999)
    ENDIF
   
    EXITS("SolverMatrix_EquationsMatrixAdd")
    RETURN
999 ERRORSEXITS("SolverMatrix_EquationsMatrixAdd",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMatrix_EquationsMatrixAdd

  !
  !================================================================================================================================
  !

  !>Adds alpha times the interface matrix into the solver matrix
  SUBROUTINE SolverMatrix_InterfaceMatrixAdd(solverMatrix,interfaceConditionIdx,alpha,interfaceMatrix,err,error,*)

    !Argument variables
    TYPE(SolverMatrixType), POINTER :: solverMatrix !<A pointer to the solver matrix
    INTEGER(INTG), INTENT(IN) :: interfaceConditionIdx !<The interfaceConditionIdx index in the solver mapping that contains the interface matrix to add
    REAL(DP), INTENT(IN) :: alpha(2) !<The multiplicative factor for the interface matrix
    TYPE(InterfaceMatrixType), POINTER :: interfaceMatrix !<A pointer to the interface matrix to add    
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: interfaceMatrixNumber,solverMatrixNumber
    INTEGER(INTG), POINTER :: columnIndices(:),rowIndices(:)
    REAL(DP), POINTER :: interfaceMatrixData(:)
    LOGICAL :: hasTranspose
    TYPE(DistributedMatrixType), POINTER :: interfaceDistributedMatrix,solverDistributedMatrix,transposeDistributedMatrix
    TYPE(InterfaceConditionToSolverMatricesMapType), POINTER :: interfaceConditionToSolverMatricesMap
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices
    TYPE(InterfaceMatricesToSolverMatrixMapType), POINTER :: interfaceMatricesToSolverMatrixMap
    TYPE(InterfaceMatrixToSolverMatricesMapType), POINTER :: interfaceMatrixToSolverMatricesMap
    TYPE(InterfaceMatrixToSolverMatrixMapType), POINTER :: interfaceMatrixToSolverMatrixMap
    TYPE(MatrixRowColCouplingType), POINTER :: interfaceColToSolverColsMap(:),interfaceColToSolverRowsMap(:), &
      & interfaceRowToSolverColsMap(:),interfaceRowToSolverRowsMap(:)
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolverMatricesType), POINTER :: solverMatrices
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("SolverMatrix_InterfaceMatrixAdd",err,error,*999)

    IF(ABS(alpha(1))>ZERO_TOLERANCE) THEN
      NULLIFY(solverMatrices)
      CALL SolverMatrix_SolverMatricesGet(solverMatrix,solverMatrices,err,error,*999)
      CALL SolverMatrices_AssertIsFinished(solverMatrices,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverMatrices_SolverMappingGet(solverMatrices,solverMapping,err,error,*999)      
      NULLIFY(interfaceMatrices)
      CALL InterfaceMatrix_InterfaceMatricesGet(interfaceMatrix,interfaceMatrices,err,error,*999)
      CALL InterfaceMatrices_AssertIsFinished(interfaceMatrices,err,error,*999)
      NULLIFY(solverDistributedMatrix)
      CALL SolverMatrix_SolverDistributedMatrixGet(solverMatrix,solverDistributedMatrix,err,error,*999)
      NULLIFY(interfaceDistributedMatrix)
      CALL InterfaceMatrix_DistributedMatrixGet(interfaceMatrix,interfaceDistributedMatrix,err,error,*999)
      NULLIFY(interfaceConditionToSolverMatricesMap)
      CALL SolverMapping_InterfaceConditionToSolverMatricesMapGet(solverMapping,interfaceConditionIdx, &
        & interfaceConditionToSolverMatricesMap,err,error,*999)
      CALL InterfaceMatrix_MatrixNumberGet(interfaceMatrix,interfaceMatrixNumber,err,error,*999)
      NULLIFY(interfaceMatrixToSolverMatricesMap)
      CALL SolverMappingICToSMSMap_InterfaceMatrixToSolverMatricesMapGet(interfaceConditionToSolverMatricesMap, &
        & interfaceMatrixNumber,interfaceMatrixToSolverMatricesMap,err,error,*999)
      CALL SolverMatrix_MatrixNumberGet(solverMatrix,solverMatrixNumber,err,error,*999)
      NULLIFY(interfaceMatricesToSolverMatrixMap)
      CALL SolverMappingICToSMSMap_InterfaceMatricesToSolverMatrixMapGet(interfaceConditionToSolverMatricesMap, &
        & solverMatrixNumber,interfaceMatricesToSolverMatrixMap,err,error,*999)
      NULLIFY(interfaceRowToSolverRowsMap)
      CALL SolverMappingIMToSMSMap_InterfaceRowToSolverRowsMapGet(interfaceMatrixToSolverMatricesMap,interfaceRowToSolverRowsMap, &
        & err,error,*999)
      NULLIFY(interfaceColToSolverColsMap)
      CALL SolverMappingIMSToSMMap_InterfaceColToSolverColsMapGet(interfaceMatricesToSolverMatrixMap,interfaceColToSolverColsMap, &
        & err,error,*999)
      CALL DistributedMatrix_MatrixCoupleAdd(solverDistributedMatrix,DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE, &
        & interfaceRowToSolverRowsMap,interfaceColToSolverColsMap,alpha(1),interfaceDistributedMatrix,.FALSE.,err,error,*999)
      CALL InterfaceMatrix_HasTransposeGet(interfaceMatrix,hasTranspose,err,error,*999)
      IF(hasTranspose) THEN
        IF(ABS(alpha(2))>ZERO_TOLERANCE) THEN
          NULLIFY(transposeDistributedMatrix)
          CALL InterfaceMatrix_TransposeDistributedMatrixGet(interfaceMatrix,transposeDistributedMatrix,err,error,*999)
          NULLIFY(interfaceMatrixToSolverMatrixMap)
          CALL SolverMappingIMToSMSMap_InterfaceMatrixToSolverMatrixMapGet(interfaceMatrixToSolverMatricesMap,solverMatrixNumber, &
            & interfaceMatrixToSolverMatrixMap,err,error,*999)
          NULLIFY(interfaceRowToSolverColsMap)
          CALL SolverMappingIMToSMMap_InterfaceRowToSolverColsMapGet(interfaceMatrixToSolverMatrixMap,interfaceRowToSolverColsMap, &
            & err,error,*999)
          NULLIFY(interfaceColToSolverRowsMap)
          CALL SolverMappingICToSMSMap_InterfaceColToSolverRowsMapGet(interfaceConditionToSolverMatricesMap, &
            & interfaceColToSolverRowsMap,err,error,*999)
          CALL DistributedMatrix_MatrixCoupleAdd(solverDistributedMatrix,DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE, &
            & interfaceColToSolverRowsMap,interfaceRowToSolverColsMap,alpha(2),transposeDistributedMatrix,.FALSE., &
            & err,error,*999)
        ENDIF
      ENDIF !Interface matrix transpose
    ENDIF
    
    EXITS("SolverMatrix_InterfaceMatrixAdd")
    RETURN
999 ERRORSEXITS("SolverMatrix_InterfaceMatrixAdd",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMatrix_InterfaceMatrixAdd

  !
  !================================================================================================================================
  !

  !>Adds alpha times the Jacobian matrix into the solver matrix
  SUBROUTINE SolverMatrix_JacobianMatrixAdd(solverMatrix,equationsSetIdx,alpha,jacobianMatrix,err,error,*)

    !Argument variables
    TYPE(SolverMatrixType), POINTER :: solverMatrix !<A pointer to the solver matrix
    INTEGER(INTG), INTENT(IN) :: equationsSetIdx !<The equations set index in the solver mapping that contains the Jacobian matrix to add
    REAL(DP), INTENT(IN) :: alpha !<The multiplicative factor for the Jacobian matrix
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix !<A pointer to the Jacobian matrix to add    
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: jacobianMatrixNumber
    REAL(DP), POINTER :: jacobianMatrixData(:)
    TYPE(DistributedMatrixType), POINTER :: jacobianDistributedMatrix,solverDistributedMatrix
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector
    TYPE(EquationsSetToSolverMatricesMapType), POINTER :: equationsSetToSolverMatricesMap
    TYPE(JacobianMatrixToSolverMatrixMapType), POINTER :: jacobianMatrixToSolverMatrixMap
    TYPE(MatrixRowColCouplingType), POINTER :: equationsRowToSolverRowsMap(:),jacobianColToSolverColsMap(:)
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolverMatricesType), POINTER :: solverMatrices
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("SolverMatrix_JacobianMatrixAdd",err,error,*999)

    IF(.NOT.ASSOCIATED(solverMatrix)) CALL FlagError("Solver matrix is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(jacobianMatrix)) CALL FlagError("Jacobian matrix is not associated.",err,error,*999)
    
    IF(ABS(alpha)>ZERO_TOLERANCE) THEN
      NULLIFY(solverMatrices)
      CALL SolverMatrix_SolverMatricesGet(solverMatrix,solverMatrices,err,error,*999)
      CALL SolverMatrices_AssertIsFinished(solverMatrices,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverMatrices_SolverMappingGet(solverMatrices,solverMapping,err,error,*999)
      NULLIFY(equationsSetToSolverMatricesMap)
      CALL SolverMapping_EquationsSetToSolverMatricesMapGet(solverMapping,equationsSetIdx,equationsSetToSolverMatricesMap, &
        & err,error,*999)
      NULLIFY(equationsRowToSolverRowsMap)
      CALL SolverMappingESToSMSMap_EquationsRowToSolverRowsMapGet(equationsSetToSolverMatricesMap,equationsRowToSolverRowsMap, &
        & err,error,*999)
      CALL JacobianMatrix_MatrixNumberGet(jacobianMatrix,jacobianMatrixNumber,err,error,*999)
      NULLIFY(jacobianMatrixToSolverMatrixMap)
      CALL SolverMappingESToSMSMap_JacobianMatrixToSolverMatrixMapGet(equationsSetToSolverMatricesMap, &
        & jacobianMatrixNumber,jacobianMatrixToSolverMatrixMap,err,error,*999)
      NULLIFY(jacobianColToSolverColsMap)
      CALL SolverMappingJMToSMMap_JacobianColToSolverColsMapGet(jacobianMatrixToSolverMatrixMap,jacobianColToSolverColsMap, &
        & err,error,*999)
      NULLIFY(residualVector)
      CALL JacobianMatrix_ResidualVectorGet(jacobianMatrix,residualVector,err,error,*999)
      NULLIFY(nonlinearMatrices)
      CALL EquationsMatricesResidual_NonlinearMatricesGet(residualVector,nonlinearMatrices,err,error,*999)
      NULLIFY(vectorMatrices)
      CALL EquationsMatricesNonlinear_VectorMatricesGet(nonlinearMatrices,vectorMatrices,err,error,*999)
      CALL EquationsMatricesVector_AssertIsFinished(vectorMatrices,err,error,*999)
      NULLIFY(solverDistributedMatrix)
      CALL SolverMatrix_SolverDistributedMatrixGet(solverMatrix,solverDistributedMatrix,err,error,*999)
      NULLIFY(jacobianDistributedMatrix)
      CALL JacobianMatrix_DistributedMatrixGet(jacobianMatrix,jacobianDistributedMatrix,err,error,*999)
      CALL DistributedMatrix_MatrixCoupleAdd(solverDistributedMatrix,DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE, &
        & equationsRowToSolverRowsMap,jacobianColToSolverColsMap,alpha,jacobianDistributedMatrix,.FALSE.,err,error,*999)
    ENDIF
    
    EXITS("SolverMatrix_JacobianMatrixAdd")
    RETURN
999 ERRORSEXITS("SolverMatrix_JacobianMatrixAdd",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMatrix_JacobianMatrixAdd

  !
  !================================================================================================================================
  !

  !>Calculates the structure (sparsity) of the solver matrix from the soluton mapping.
  SUBROUTINE SolverMatrix_StructureCalculate(solverMatrix,numberOfNonZeros,rowIndices,columnIndices,err,error,*)

    !Argument variables
    TYPE(SolverMatrixType), POINTER :: solverMatrix !<A pointer to the solver matrix to calculate the structure for
    INTEGER(INTG), INTENT(OUT) :: numberOfNonZeros !<On return the number of non-zeros in the solver matrix
    INTEGER(INTG), POINTER :: rowIndices(:) !<On return a pointer to row location indices in compressed row format. The pointers must be NULL on entry and the calling routine is responsible for deallocation.
    INTEGER(INTG), POINTER :: columnIndices(:) !<On return a pointer to the column location indices in compressed row format. The pointers must be NULL on entry and the calling routine is responsible for deallocation.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG)  :: dummyErr,equationsMatrixIdx,equationsSetIdx,equationsSolverMatrixNumber,interfaceConditionIdx, &
      & interfaceConditionMethod,interfaceMatrixIdx,maximumColumnIndices,maximumColumnsPerRow,maximumTransposeColumnsPerRow, &
      & numberOfColumns,numberOfDynamicMatrices,numberOfEquationsSets,numberOfInterfaceConditions,numberOfInterfaceMatrices, &
      & numberOfJacobianMatrices,numberOfLinearMatrices,numberOfMatrices,numberOfRows,numberOfSolverMatrices,solverColumnIdx, &
      & solverMatrixIdx,solverMatrixNumber,solverRowNumber,storageType
    INTEGER(INTG), ALLOCATABLE :: columns(:)
    REAL(DP) :: sparsity
    LOGICAL :: hasTranspose
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix,solverDistributedMatrix,transposeDistributedMatrix
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesToSolverMatrixMapType), POINTER :: equationsMatricesToSolverMatrixMap
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatrixToSolverMatrixMapType), POINTER :: equationsMatrixToSolverMatrixMap
    TYPE(EquationsSetToSolverMatricesMapType), POINTER :: equationsSetToSolverMatricesMap
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix
    TYPE(JacobianMatrixToSolverMatrixMapType), POINTER :: jacobianMatrixToSolverMatrixMap
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition
    TYPE(InterfaceConditionToSolverMatricesMapType), POINTER :: interfaceConditionToSolverMatricesMap
    TYPE(InterfaceMatrixType), POINTER :: interfaceMatrix
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices
    TYPE(InterfaceMatricesToSolverMatrixMapType), POINTER :: interfaceMatricesToSolverMatrixMap
    TYPE(InterfaceMatrixToSolverMatricesMapType), POINTER :: interfaceMatrixToSolverMatricesMap
    TYPE(InterfaceMatrixToSolverMatrixMapType), POINTER :: interfaceMatrixToSolverMatrixMap
    TYPE(ListPtrType), ALLOCATABLE :: columnIndicesLists(:)
    TYPE(MatrixRowColCouplingType), POINTER :: equationsColToSolverColsMap(:),equationsRowToSolverRowsMap(:), &
      & jacobianColToSolverColsMap(:),interfaceColToSolverColsMap(:),interfaceColToSolverRowsMap(:), &
      & interfaceRowToSolverColsMap(:),interfaceRowToSolverRowsMap(:)
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolverMatricesType), POINTER :: solverMatrices
    TYPE(SolverMatrixToEquationsMapType), POINTER :: solverMatrixToEquationsMap
    TYPE(VARYING_STRING) :: dummyError,localError
    
    ENTERS("SolverMatrix_StructureCalculate",err,error,*999)

    numberOfNonZeros=0
    IF(ASSOCIATED(rowIndices)) CALL FlagError("Row indices is already associated.",err,error,*998)
    IF(ASSOCIATED(columnIndices)) CALL FlagError("Row indices is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverMatrix)) CALL FlagError("Solver matrix is not associated.",err,error,*999)
    
    CALL SolverMatrix_MatrixNumberGet(solverMatrix,solverMatrixNumber,err,error,*999)
    NULLIFY(solverDistributedMatrix)
    CALL SolverMatrix_SolverDistributedMatrixGet(solverMatrix,solverDistributedMatrix,err,error,*999)
    CALL DistributedMatrix_AssertNotFinished(solverDistributedMatrix,err,error,*999)
    NULLIFY(solverMatrices)
    CALL SolverMatrix_SolverMatricesGet(solverMatrix,solverMatrices,err,error,*999)
    NULLIFY(solverMapping)
    CALL SolverMatrices_SolverMappingGet(solverMatrices,solverMapping,err,error,*999)
    NULLIFY(solverMatrixToEquationsMap)
    CALL SolverMapping_SolverMatrixToEquationsMapGet(solverMapping,solverMatrixNumber,solverMatrixToEquationsMap,err,error,*999)
    CALL SolverMappingSMToEQSMap_NumberOfEquationsSetsGet(solverMatrixToEquationsMap,numberOfEquationsSets,err,error,*999)
    CALL SolverMappingSMToEQSMap_NumberOfInterfaceConditionsGet(solverMatrixToEquationsMap,numberOfInterfaceConditions, &
      & err,error,*999)    
    CALL SolverMatrix_StorageTypeGet(solverMatrix,storageType,err,error,*999)
    SELECT CASE(storageType)
    CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
      CALL FlagError("Can not calculate the structure for a block storage matrix.",err,error,*999)
    CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
      CALL FlagError("Can not calcualte the structure for a diagonal matrix.",err,error,*999)
    CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
      !Find the maximum number of column indices
      maximumColumnIndices=0
      DO equationsSetIdx=1,numberOfEquationsSets
         NULLIFY(equationsSetToSolverMatricesMap)
        CALL SolverMapping_EquationsSetToSolverMatricesMapGet(solverMapping,equationsSetIdx,equationsSetToSolverMatricesMap, &
          & err,error,*999)
        CALL SolverMappingESToSMSMap_NumberOfSolverMatricesGet(equationsSetToSolverMatricesMap,numberOfSolverMatrices, &
          & err,error,*999)
        DO solverMatrixIdx=1,numberOfSolverMatrices
          NULLIFY(equationsMatricesToSolverMatrixMap)
          CALL SolverMappingESToSMSMap_EquationsMatricesToSolverMatrixMapGet(equationsSetToSolverMatricesMap,solverMatrixIdx, &
            & equationsMatricesToSolverMatrixMap,err,error,*999)
          CALL SolverMappingEMSToSMMap_SolverMatrixNumberGet(equationsMatricesToSolverMatrixMap,equationsSolverMatrixNumber, &
            & err,error,*999)
          IF(equationsSolverMatrixNumber==solverMatrixNumber) THEN
            !Loop over dynamic matrices mapped to the solver matrix
            CALL SolverMappingEMSToSMMap_NumberOfDynamicMatricesGet(equationsMatricesToSolverMatrixMap,numberOfDynamicMatrices, &
              & err,error,*999)
            DO equationsMatrixIdx=1,numberOfDynamicMatrices
              NULLIFY(equationsMatrixToSolverMatrixMap)
              CALL SolverMappingEMSToSMMap_DynamicMatrixToSolverMatrixMapGet(equationsMatricesToSolverMatrixMap, &
                & equationsMatrixIdx,equationsMatrixToSolverMatrixMap,err,error,*999)
              CALL SolverMappingEMToSMMap_EquationsMatrixGet(equationsMatrixToSolverMatrixMap,equationsMatrix,err,error,*999)
              NULLIFY(distributedMatrix)
              CALL EquationsMatrix_DistributedMatrixGet(equationsMatrix,distributedMatrix,err,error,*999)
              CALL DistributedMatrix_MaxColumnsPerRowGet(distributedMatrix,maximumColumnsPerRow,err,error,*999)
              maximumColumnIndices=maximumColumnIndices+maximumColumnsPerRow
            ENDDO !equationsMatrixIdx
            !Loop over linear matrices mapped to the solver matrix
            CALL SolverMappingEMSToSMMap_NumberOfLinearMatricesGet(equationsMatricesToSolverMatrixMap,numberOfLinearMatrices, &
              & err,error,*999)
            DO equationsMatrixIdx=1,numberOfLinearMatrices
              NULLIFY(equationsMatrixToSolverMatrixMap)
              CALL SolverMappingEMSToSMMap_LinearMatrixToSolverMatrixMapGet(equationsMatricesToSolverMatrixMap, &
                & equationsMatrixIdx,equationsMatrixToSolverMatrixMap,err,error,*999)
              CALL SolverMappingEMToSMMap_EquationsMatrixGet(equationsMatrixToSolverMatrixMap,equationsMatrix,err,error,*999)
              NULLIFY(distributedMatrix)
              CALL EquationsMatrix_DistributedMatrixGet(equationsMatrix,distributedMatrix,err,error,*999)
              CALL DistributedMatrix_MaxColumnsPerRowGet(distributedMatrix,maximumColumnsPerRow,err,error,*999)
              maximumColumnIndices=maximumColumnIndices+maximumColumnsPerRow
            ENDDO !equationsMatrixIdx
            !Loop over Jacobian matrices mapped to the solver matrix
            CALL SolverMappingEMSToSMMap_NumberOfJacobianMatricesGet(equationsMatricesToSolverMatrixMap,numberOfJacobianMatrices, &
              & err,error,*999)
            DO equationsMatrixIdx=1,numberOfJacobianMatrices
              NULLIFY(jacobianMatrixToSolverMatrixMap)
              CALL SolverMappingEMSToSMMap_JacobianMatrixToSolverMatrixMapGet(equationsMatricesToSolverMatrixMap, &
                & equationsMatrixIdx,jacobianMatrixToSolverMatrixMap,err,error,*999)
              CALL SolverMappingJMToSMMap_JacobianMatrixGet(jacobianMatrixToSolverMatrixMap,jacobianMatrix,err,error,*999)
              NULLIFY(distributedMatrix)
              CALL JacobianMatrix_DistributedMatrixGet(jacobianMatrix,distributedMatrix,err,error,*999)
              CALL DistributedMatrix_MaxColumnsPerRowGet(distributedMatrix,maximumColumnsPerRow,err,error,*999)
              maximumColumnIndices=maximumColumnIndices+maximumColumnsPerRow
            ENDDO !equationsMatrixIdx
          ENDIF
        ENDDO !solverMatrixIdx
      ENDDO !equationsSetIdx
      !Loop over any interface conditions
      CALL SolverMapping_NumberOfInterfaceConditionsGet(solverMapping,numberOfInterfaceConditions,err,error,*999)
      DO interfaceConditionIdx=1,numberOfInterfaceConditions
        NULLIFY(interfaceCondition)
        CALL SolverMapping_InterfaceConditionGet(solverMapping,interfaceConditionIdx,interfaceCondition,err,error,*999)
        CALL InterfaceCondition_MethodGet(interfaceCondition,interfaceConditionMethod,err,error,*999)
        SELECT CASE(interfaceConditionMethod)
        CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
          NULLIFY(interfaceConditionToSolverMatricesMap)
          CALL SolverMapping_InterfaceConditionToSolverMatricesMapGet(solverMapping,interfaceConditionIdx, &
            & interfaceConditionToSolverMatricesMap,err,error,*999)
          CALL SolverMappingICToSMSMap_NumberOfSolverMatricesGet(interfaceConditionToSolverMatricesMap,numberOfSolverMatrices, &
            & err,error,*999)
          DO solverMatrixIdx=1,numberOfSolverMatrices
            NULLIFY(interfaceMatricesToSolverMatrixMap)
            CALL SolverMappingICToSMSMap_InterfaceMatricesToSolverMatrixMapGet(interfaceConditionToSolverMatricesMap, &
              & solverMatrixIdx,interfaceMatricesToSolverMatrixMap,err,error,*999)
            CALL SolverMappingIMSToSMMap_SolverMatrixNumberGet(interfaceMatricesToSolverMatrixMap,equationsSolverMatrixNumber, &
              & err,error,*999)
            IF(equationsSolverMatrixNumber==solverMatrixNumber) THEN
              !Loop over interface matrices mapped to the solver matrix
              CALL SolverMappingIMSToSMMap_NumberOfInterfaceMatricesGet(interfaceMatricesToSolverMatrixMap, &
                & numberOfInterfaceMatrices,err,error,*999)
              DO interfaceMatrixIdx=1,numberOfInterfaceMatrices
                NULLIFY(interfaceMatrixToSolverMatrixMap)
                CALL SolverMappingIMSToSMMap_InterfaceMatrixToSolverMatrixMapGet(interfaceMatricesToSolverMatrixMap, &
                  & interfaceMatrixIdx,interfaceMatrixToSolverMatrixMap,err,error,*999)
                NULLIFY(interfaceMatrix)
                CALL SolverMappingIMToSMMap_InterfaceMatrixGet(interfaceMatrixToSolverMatrixMap,interfaceMatrix,err,error,*999)
                NULLIFY(distributedMatrix)
                CALL InterfaceMatrix_DistributedMatrixGet(interfaceMatrix,distributedMatrix,err,error,*999)
                CALL DistributedMatrix_MaxColumnsPerRowGet(distributedMatrix,maximumColumnsPerRow,err,error,*999)
                maximumTransposeColumnsPerRow=0
                CALL InterfaceMatrix_HasTransposeGet(interfaceMatrix,hasTranspose,err,error,*999)
                IF(hasTranspose) THEN
                  NULLIFY(transposeDistributedMatrix)
                  CALL InterfaceMatrix_TransposeDistributedMatrixGet(interfaceMatrix,transposeDistributedMatrix,err,error,*999)
                  CALL DistributedMatrix_MaxColumnsPerRowGet(transposeDistributedMatrix,maximumTransposeColumnsPerRow, &
                    & err,error,*999)
                ENDIF
                maximumColumnIndices=maximumColumnIndices+MAX(maximumColumnsPerRow,maximumTransposeColumnsPerRow)
              ENDDO !interfaceMatrixIdx
            ENDIF
          ENDDO !solverMatrixIdx          
        CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE DEFAULT
          localError="The interface condition method of "//TRIM(NumberToVString(interfaceConditionMethod,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDDO !interfaceConditionIdx
      !Allocate lists
      CALL SolverMapping_NumberOfRowsGet(solverMapping,numberOfRows,err,error,*999)
      ALLOCATE(columnIndicesLists(numberOfRows),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate column indices lists.",err,error,*999)
      !Allocate row indices
      ALLOCATE(rowIndices(numberOfRows+1),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate row indices.",err,error,*999)
      rowIndices(1)=1
      !Set up the column indicies lists
      DO solverRowNumber=1,numberOfRows
        NULLIFY(columnIndicesLists(solverRowNumber)%ptr)
        CALL List_CreateStart(columnIndicesLists(solverRowNumber)%ptr,err,error,*999)
        CALL List_DataTypeSet(columnIndicesLists(solverRowNumber)%ptr,LIST_INTG_TYPE,err,error,*999)
        CALL List_InitialSizeSet(columnIndicesLists(solverRowNumber)%ptr,maximumColumnIndices,err,error,*999)
        CALL List_CreateFinish(columnIndicesLists(solverRowNumber)%ptr,err,error,*999)
      ENDDO !solverRowNumber
      !Loop over the equations sets
      DO equationsSetIdx=1,numberOfEquationsSets
         NULLIFY(equationsSetToSolverMatricesMap)
        CALL SolverMapping_EquationsSetToSolverMatricesMapGet(solverMapping,equationsSetIdx,equationsSetToSolverMatricesMap, &
          & err,error,*999)
        CALL SolverMappingESToSMSMap_NumberOfSolverMatricesGet(equationsSetToSolverMatricesMap,numberOfSolverMatrices, &
          & err,error,*999)
        DO solverMatrixIdx=1,numberOfSolverMatrices
          NULLIFY(equationsMatricesToSolverMatrixMap)
          CALL SolverMappingESToSMSMap_EquationsMatricesToSolverMatrixMapGet(equationsSetToSolverMatricesMap,solverMatrixIdx, &
            & equationsMatricesToSolverMatrixMap,err,error,*999)
          CALL SolverMappingEMSToSMMap_SolverMatrixNumberGet(equationsMatricesToSolverMatrixMap,equationsSolverMatrixNumber, &
            & err,error,*999)
          IF(equationsSolverMatrixNumber==solverMatrixNumber) THEN
            NULLIFY(equationsRowToSolverRowsMap)
            CALL SolverMappingESToSMSMap_EquationsRowToSolverRowsMapGet(equationsSetToSolverMatricesMap, &
              & equationsRowToSolverRowsMap,err,error,*999)
            !Loop over dynamic matrices mapped to the solver matrix
            CALL SolverMappingEMSToSMMap_NumberOfDynamicMatricesGet(equationsMatricesToSolverMatrixMap,numberOfDynamicMatrices, &
              & err,error,*999)
            DO equationsMatrixIdx=1,numberOfDynamicMatrices
              NULLIFY(equationsMatrixToSolverMatrixMap)
              CALL SolverMappingEMSToSMMap_DynamicMatrixToSolverMatrixMapGet(equationsMatricesToSolverMatrixMap, &
                & equationsMatrixIdx,equationsMatrixToSolverMatrixMap,err,error,*999)
              NULLIFY(equationsMatrix)
              CALL SolverMappingEMToSMMap_EquationsMatrixGet(equationsMatrixToSolverMatrixMap,equationsMatrix,err,error,*999)
              NULLIFY(distributedMatrix)
              CALL EquationsMatrix_DistributedMatrixGet(equationsMatrix,distributedMatrix,err,error,*999)
              NULLIFY(equationsColToSolverColsMap)
              CALL SolverMappingEMToSMMap_EquationsColToSolverColsMapGet(equationsMatrixToSolverMatrixMap, &
                & equationsColToSolverColsMap,err,error,*999)
              CALL DistributedMatrix_MatrixStructureCoupleCalculate(distributedMatrix,.FALSE.,equationsRowToSolverRowsMap, &
                & equationsColToSolverColsMap,columnIndicesLists,err,error,*999)          
            ENDDO !equationsMatrixIdx
            !Loop over linear matrices mapped to the solver matrix
            CALL SolverMappingEMSToSMMap_NumberOfLinearMatricesGet(equationsMatricesToSolverMatrixMap,numberOfLinearMatrices, &
              & err,error,*999)
            DO equationsMatrixIdx=1,numberOfLinearMatrices
              NULLIFY(equationsMatrixToSolverMatrixMap)
              CALL SolverMappingEMSToSMMap_LinearMatrixToSolverMatrixMapGet(equationsMatricesToSolverMatrixMap, &
                & equationsMatrixIdx,equationsMatrixToSolverMatrixMap,err,error,*999)
              NULLIFY(equationsMatrix)
              CALL SolverMappingEMToSMMap_EquationsMatrixGet(equationsMatrixToSolverMatrixMap,equationsMatrix,err,error,*999)
              NULLIFY(distributedMatrix)
              CALL EquationsMatrix_DistributedMatrixGet(equationsMatrix,distributedMatrix,err,error,*999)
              NULLIFY(equationsColToSolverColsMap)
              CALL SolverMappingEMToSMMap_EquationsColToSolverColsMapGet(equationsMatrixToSolverMatrixMap, &
                & equationsColToSolverColsMap,err,error,*999)
              CALL DistributedMatrix_MatrixStructureCoupleCalculate(distributedMatrix,.FALSE.,equationsRowToSolverRowsMap, &
                & equationsColToSolverColsMap,columnIndicesLists,err,error,*999)          
            ENDDO !equationsMatrixIdx
            !Loop over Jacobian matrices mapped to the solver matrix
            CALL SolverMappingEMSToSMMap_NumberOfJacobianMatricesGet(equationsMatricesToSolverMatrixMap,numberOfJacobianMatrices, &
              & err,error,*999)
            DO equationsMatrixIdx=1,numberOfJacobianMatrices
              NULLIFY(jacobianMatrixToSolverMatrixMap)
              CALL SolverMappingEMSToSMMap_JacobianMatrixToSolverMatrixMapGet(equationsMatricesToSolverMatrixMap, &
                & equationsMatrixIdx,jacobianMatrixToSolverMatrixMap,err,error,*999)
              NULLIFY(jacobianMatrix)
              CALL SolverMappingJMToSMMap_JacobianMatrixGet(jacobianMatrixToSolverMatrixMap,jacobianMatrix,err,error,*999)
              NULLIFY(distributedMatrix)
              CALL JacobianMatrix_DistributedMatrixGet(jacobianMatrix,distributedMatrix,err,error,*999)
              NULLIFY(jacobianColToSolverColsMap)
              CALL SolverMappingJMToSMMap_JacobianColToSolverColsMapGet(jacobianMatrixToSolverMatrixMap, &
                & jacobianColToSolverColsMap,err,error,*999)            
              CALL DistributedMatrix_MatrixStructureCoupleCalculate(distributedMatrix,.FALSE.,equationsRowToSolverRowsMap, &
                & jacobianColToSolverColsMap,columnIndicesLists,err,error,*999)
            ENDDO !equationsMatrixIdx
          ENDIF
        ENDDO !solverMatrixIdx
      ENDDO !equationsSetIdx
      !Loop over any interface conditions
      CALL SolverMapping_NumberOfInterfaceConditionsGet(solverMapping,numberOfInterfaceConditions,err,error,*999)
      DO interfaceConditionIdx=1,numberOfInterfaceConditions
        NULLIFY(interfaceCondition)
        CALL SolverMapping_InterfaceConditionGet(solverMapping,interfaceConditionIdx,interfaceCondition,err,error,*999)
        CALL InterfaceCondition_MethodGet(interfaceCondition,interfaceConditionMethod,err,error,*999)
        SELECT CASE(interfaceConditionMethod)
        CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
          NULLIFY(interfaceConditionToSolverMatricesMap)
          CALL SolverMapping_InterfaceConditionToSolverMatricesMapGet(solverMapping,interfaceConditionIdx, &
            & interfaceConditionToSolverMatricesMap,err,error,*999)
          NULLIFY(interfaceColToSolverRowsMap)
          CALL SolverMappingICToSMSMap_InterfaceColToSolverRowsMapGet(interfaceConditionToSolverMatricesMap, &
            & interfaceColToSolverRowsMap,err,error,*999)
          CALL SolverMappingICToSMSMap_NumberOfSolverMatricesGet(interfaceConditionToSolverMatricesMap,numberOfSolverMatrices, &
            & err,error,*999)
          DO solverMatrixIdx=1,numberOfSolverMatrices
            NULLIFY(interfaceMatricesToSolverMatrixMap)
            CALL SolverMappingICToSMSMap_InterfaceMatricesToSolverMatrixMapGet(interfaceConditionToSolverMatricesMap, &
              & solverMatrixIdx,interfaceMatricesToSolverMatrixMap,err,error,*999)
            CALL SolverMappingIMSToSMMap_SolverMatrixNumberGet(interfaceMatricesToSolverMatrixMap,equationsSolverMatrixNumber, &
              & err,error,*999)
            IF(equationsSolverMatrixNumber==solverMatrixNumber) THEN
              NULLIFY(interfaceColToSolverColsMap)
              CALL SolverMappingIMSToSMMap_InterfaceColToSolverColsMapGet(interfaceMatricesToSolverMatrixMap, &
                & interfaceColToSolverColsMap,err,error,*999)
              !Loop over interface matrices mapped to the solver matrix
              CALL SolverMappingIMSToSMMap_NumberOfInterfaceMatricesGet(interfaceMatricesToSolverMatrixMap, &
                & numberOfInterfaceMatrices,err,error,*999)
              DO interfaceMatrixIdx=1,numberOfInterfaceMatrices
                NULLIFY(interfaceMatrixToSolverMatricesMap)
                CALL SolverMappingICToSMSMap_InterfaceMatrixToSolverMatricesMapGet(interfaceConditionToSolverMatricesMap, &
                  & interfaceMatrixIdx,interfaceMatrixToSolverMatricesMap,err,error,*999)
                NULLIFY(interfaceRowToSolverRowsMap)
                CALL SolverMappingIMToSMSMap_InterfaceRowToSolverRowsMapGet(interfaceMatrixToSolverMatricesMap, &
                  & interfaceRowToSolverRowsMap,err,error,*999)            
                NULLIFY(interfaceMatrixToSolverMatrixMap)
                CALL SolverMappingIMSToSMMap_InterfaceMatrixToSolverMatrixMapGet(interfaceMatricesToSolverMatrixMap, &
                  & interfaceMatrixIdx,interfaceMatrixToSolverMatrixMap,err,error,*999)
                NULLIFY(interfaceMatrix)
                CALL SolverMappingIMToSMMap_InterfaceMatrixGet(interfaceMatrixToSolverMatrixMap,interfaceMatrix,err,error,*999)
                NULLIFY(distributedMatrix)
                CALL InterfaceMatrix_DistributedMatrixGet(interfaceMatrix,distributedMatrix,err,error,*999)
                CALL DistributedMatrix_MatrixStructureCoupleCalculate(distributedMatrix,.FALSE.,interfaceRowToSolverRowsMap, &
                  &  interfaceColToSolverColsMap,columnIndicesLists,err,error,*999)
                CALL InterfaceMatrix_HasTransposeGet(interfaceMatrix,hasTranspose,err,error,*999)
                IF(hasTranspose) THEN
                  NULLIFY(transposeDistributedMatrix)
                  CALL InterfaceMatrix_TransposeDistributedMatrixGet(interfaceMatrix,transposeDistributedMatrix,err,error,*999)
                  NULLIFY(interfaceRowToSolverColsMap)
                  CALL SolverMappingIMToSMMap_InterfaceRowToSolverColsMapGet(interfaceMatrixToSolverMatrixMap, &
                    & interfaceRowToSolverColsMap,err,error,*999)
                  CALL DistributedMatrix_MatrixStructureCoupleCalculate(transposeDistributedMatrix,.FALSE., &
                    & interfaceColToSolverRowsMap,interfaceRowToSolverColsMap,columnIndicesLists,err,error,*999)              
                ENDIF
              ENDDO !interfaceMatrixIdx
            ENDIF
          ENDDO !solverMatrixIdx          
        CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE DEFAULT
          localError="The interface condition method of "//TRIM(NumberToVString(interfaceConditionMethod,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDDO !interfaceConditionIdx
      !Loop over the rows to calculate the number of non-zeros and setup the row indicces
      DO solverRowNumber=1,numberOfRows
        CALL List_RemoveDuplicates(columnIndicesLists(solverRowNumber)%ptr,err,error,*999)
        CALL List_NumberOfItemsGet(columnIndicesLists(solverRowNumber)%ptr,numberOfColumns,err,error,*999)
        numberOfNonZeros=numberOfNonZeros+numberOfColumns
        rowIndices(solverRowNumber+1)=numberOfNonZeros+1
      ENDDO !solverRowNumber
      !Allocate and setup the column locations
      ALLOCATE(columnIndices(numberOfNonZeros),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate column indices.",err,error,*999)
      DO solverRowNumber=1,numberOfRows
        CALL List_DetachAndDestroy(columnIndicesLists(solverRowNumber)%ptr,numberOfColumns,columns,err,error,*999)
        DO solverColumnIdx=1,numberOfColumns
          columnIndices(rowIndices(solverRowNumber)+solverColumnIdx-1)=columns(solverColumnIdx)
        ENDDO !solverColumnIdx
        DEALLOCATE(columns)
      ENDDO !solver_row_idx
    CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)                        
    CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)                      
    CASE DEFAULT
      localError="The matrix storage type of "// &
        & TRIM(NumberToVString(solverMatrix%storageType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    IF(diagnostics1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Solver matrix structure:",err,error,*999)
      CALL SolverMatrix_MatrixNumberGet(solverMatrix,solverMatrixIdx,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Solver matrix number : ",solverMatrixIdx,err,error,*999)
      CALL SolverMatrices_NumberOfRowsGet(solverMatrices,numberOfRows,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of rows = ",numberOfRows,err,error,*999)
      CALL SolverMatrix_NumberOfColumnsGet(solverMatrix,numberOfColumns,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of columns = ",numberOfColumns,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of non zeros = ",numberOfNonZeros,err,error,*999)
      IF(numberOfRows*numberOfColumns/=0) THEN
        sparsity=REAL(numberOfNonZeros,DP)/REAL(numberOfRows*numberOfColumns,DP)*100.0_DP
        CALL WriteStringFmtValue(DIAGNOSTIC_OUTPUT_TYPE,"  Sparsity (%) = ",sparsity,"F6.2", err,error,*999)
      ENDIF
      IF(diagnostics2) THEN
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfRows+1,8,8,rowIndices, &
          & '("  Row indices    :",8(X,I13))','(18X,8(X,I13))',err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfNonZeros,8,8,columnIndices, &
          & '("  Column indices :",8(X,I13))','(18X,8(X,I13))', err,error,*999)
      ENDIF
    ENDIF

    EXITS("SolverMatrix_StructureCalculate")
    RETURN
999 IF(ASSOCIATED(rowIndices)) DEALLOCATE(rowIndices)
    IF(ASSOCIATED(columnIndices)) DEALLOCATE(columnIndices)
    IF(ALLOCATED(columns)) DEALLOCATE(columns)
    IF(ALLOCATED(columnIndicesLists)) THEN
      DO solverRowNumber=1,SIZE(columnIndicesLists,1)
        IF(ASSOCIATED(columnIndicesLists(solverRowNumber)%ptr)) &
          & CALL List_Destroy(columnIndicesLists(solverRowNumber)%ptr,dummyErr,dummyError,*998)
      ENDDO !solverRowNumber
      DEALLOCATE(columnIndicesLists)
    ENDIF
998 ERRORSEXITS("SolverMatrix_StructureCalculate",err,error)    
    RETURN 1
   
  END SUBROUTINE SolverMatrix_StructureCalculate
        
  !
  !================================================================================================================================
  !

  !>Finalises a solver matrix and deallocates all memory
  SUBROUTINE SolverMatrix_Finalise(solverMatrix,err,error,*)

    !Argument variables
    TYPE(SolverMatrixType), POINTER :: solverMatrix !<A pointer to the solver matrix to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("SolverMatrix_Finalise",err,error,*999)

    IF(ASSOCIATED(solverMatrix)) THEN
      IF(ASSOCIATED(solverMatrix%matrix)) CALL DistributedMatrix_Destroy(solverMatrix%matrix,err,error,*999)
      IF(ASSOCIATED(solverMatrix%solverVector)) CALL DistributedVector_Destroy(solverMatrix%solverVector,err,error,*999)
      DEALLOCATE(solverMatrix)
    ENDIF
    
    EXITS("SolverMatrix_Finalise")
    RETURN
999 ERRORSEXITS("SolverMatrix_Finalise",err,error)    
    RETURN 1
   
  END SUBROUTINE SolverMatrix_Finalise
        
  !
  !================================================================================================================================
  !

  !>Forms a solver matrix by initialising the structure of the matrix to zero.
  SUBROUTINE SolverMatrix_Form(solverMatrix,err,error,*)

    !Argument variables
    TYPE(SolverMatrixType), POINTER :: solverMatrix !<A pointer to the solver matrix to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix
    
    ENTERS("SolverMatrix_Form",err,error,*999)

    NULLIFY(distributedMatrix)
    CALL SolverMatrix_SolverDistributedMatrixGet(solverMatrix,distributedMatrix,err,error,*999)
    CALL DistributedMatrix_Form(distributedMatrix,err,error,*999)
   
    EXITS("SolverMatrix_Form")
    RETURN
999 ERRORSEXITS("SolverMatrix_Form",err,error)    
    RETURN 1
   
  END SUBROUTINE SolverMatrix_Form
        
    !
  !================================================================================================================================
  !

  !>Initialises a solver matrix
  SUBROUTINE SolverMatrix_Initialise(solverMatrices,matrixNumber,err,error,*)

    !Argument variables
    TYPE(SolverMatricesType), POINTER :: solverMatrices !<A pointer to the solver matrices to initialise
    INTEGER(INTG), INTENT(IN) :: matrixNumber !<The matrix number in the solver matrices to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,numberOfColumns,numberOfMatrices
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolverMatrixType), POINTER :: solverMatrix
    TYPE(SolverMatrixToEquationsMapType), POINTER :: solverMatrixToEquationsMap
    TYPE(VARYING_STRING) :: dummyError,localError
    
    ENTERS("SolverMatrix_Initialise",err,error,*998)

    CALL SolverMatrices_NumberOfSolverMatricesGet(solverMatrices,numberOfMatrices,err,error,*999)
    IF(matrixNumber<1.OR.matrixNumber>numberOfMatrices) THEN
      localError="The specified matrix number of "//TRIM(NumberToVString(matrixNumber,"*",err,error))// &
        & " is invalid. The matrix number must be >= 1 and <= "// &
        & TRIM(NumberToVString(numberOfMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*998)
    ENDIF
    IF(ASSOCIATED(solverMatrices%matrices(matrixNumber)%ptr)) THEN
      localError="The solver matrix for matrix index "//TRIM(NumberToVString(matrixNumber,"*",err,error))// &
        & " is already associated."
      CALL FlagError(localError,err,error,*998)
    ENDIF
    
    NULLIFY(solverMapping)
    CALL SolverMatrices_SolverMappingGet(solverMatrices,solverMapping,err,error,*999)
    NULLIFY(solverMatrixToEquationsMap)
    CALL SolverMapping_SolverMatrixToEquationsMapGet(solverMapping,matrixNumber,solverMatrixToEquationsMap,err,error,*999)
    CALL SolverMappingSMToEQSMap_NumberOfColumnsGet(solverMatrixToEquationsMap,numberOfColumns,err,error,*999)
    
    ALLOCATE(solverMatrices%matrices(matrixNumber)%ptr,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate solver matrix.",err,error,*999)
    solverMatrix=>solverMatrices%matrices(matrixNumber)%ptr
    solverMatrix%matrixNumber=matrixNumber
    solverMatrix%solverMatrices=>solverMatrices
    solverMatrix%storageType=MATRIX_BLOCK_STORAGE_TYPE
    solverMatrix%updateMatrix=.TRUE.
    solverMatrix%numberOfColumns=numberOfColumns
    solverMatrixToEquationsMap%solverMatrix=>solverMatrix
    NULLIFY(solverMatrix%solverVector)
    NULLIFY(solverMatrix%matrix)
    
    EXITS("SolverMatrix_Initialise")
    RETURN
999 CALL SolverMatrix_Finalise(solverMatrices%matrices(matrixNumber)%ptr,dummyErr,dummyError,*998)
998 ERRORSEXITS("SolverMatrix_Initialise",err,error)    
    RETURN 1
   
  END SUBROUTINE SolverMatrix_Initialise
        
  !
  !================================================================================================================================
  !

END MODULE SolverMatricesRoutines
