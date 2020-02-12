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
  USE InterfaceConditionAccessRoutines
  USE InterfaceMatricesAccessRoutines
  USE ISO_VARYING_STRING
  USE Kinds
  USE MatrixVector
  USE ProblemAccessRoutines
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
 
  PUBLIC SolverMatrix_EquationsMatrixAdd

  PUBLIC SolverMatrix_InterfaceMatrixAdd

  PUBLIC SolverMatrix_JacobianMatrixAdd

  PUBLIC SolverMatrices_CreateFinish,SolverMatrices_CreateStart

  PUBLIC SolverMatrices_Destroy

  PUBLIC SolverMatrices_LibraryTypeGet,SolverMatrices_LibraryTypeSet

  PUBLIC SolverMatrices_Output

  PUBLIC SolverMatrices_StorageTypeSet

  PUBLIC SolverMatrices_SymmetryTypeGet,SolverMatrices_SymmetryTypeSet

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
      CALL SolverMatrices_SolverMatrixGet(solvrMatrices,matrixIdx,solverMatrix,err,error,*999)
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
          CALL SOLVER_MATRIX_FINALISE(solverMatrices%matrices(matrixIdx)%ptr,err,error,*999)
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
    INTEGER(INTG) :: dummyErr,equationsMatrixIdx,equationsSetIdx,matrixIdx
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("SolverMatrices_Initialise",err,error,*998)

    IF(.NOT.ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is not associated",err,error,*998)
    IF(ASSOCIATED(solverEquations%solverMatrices)) &
      & CALL FlagError("Solver matrices is already associated for this solver equations.",err,error,*998)

    NULLIFY(solverMapping)
    CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
    
    ALLOCATE(solverEquations%solverMatrices,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate solver matrices.",err,error,*999)
    solverEquations%solverMatrices%solverEquations=>solverEquations
    solverEquations%solverMatrices%solverMatricesFinished=.FALSE.
    solverEquations%solverMatrices%solverMapping=>solverMapping
    solverEquations%solverMatrices%numberOfRows=solverMapping%numberOfRows
    solverEquations%solverMatrices%numberOfGlobalRows=solverMapping%numberOfGlobalRows
    solverEquations%solverMatrices%solverLibraryType=0
    solverEquations%solverMatrices%matrixLibraryType=0
    solverEquations%solverMatrices%numberOfMatrices=solverMapping%numberOfSolverMatrices
    ALLOCATE(solverEquations%solverMatrices%matrices(solverMapping%numberOfSolverMatrices),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate solver matrices matrices.",err,error,*999)
    DO matrixIdx=1,solverMapping%numberOfSolverMatrices
      NULLIFY(solverEquations%solverMatrices%matrices(matrixIdx)%ptr)
      CALL SOLVER_MATRIX_INITIALISE(solverEquations%solverMatrices,matrixIdx,err,error,*999)
      DO equationsSetIdx=1,SOLVER_MAPPING%numberOfEquationsSets
        IF(ALLOCATED(solverMapping%equationsSetToSolverMatricesMap(equationsSetIdx)%equationsMatricesToSolverMatrixMaps( &
          & matrixIdx)%dynamicEquationsMatrixToSolverMatrixMaps)) THEN
          DO equationsMatrixIdx=1,SOLVER_MAPPING%equationsSetToSolverMatricesMap(equationsSetIdx)% &
            & equationsMatricesToSolverMatrixMaps(matrixIdx)%numberOfDynamicEquationsMatrices
            !Add the solver matrix to the solvers mapping
            SOLVER_MAPPING%equationsSetToSolverMatricesMap(equationsSetIdx)%equationsMatricesToSolverMatrixMaps( &
              & matrixIdx)%dynamicEquationsMatrixToSolverMatrixMaps(equationsMatrixIdx)%ptr%solverMatrix=> &
              & solverEquations%solverMatrices%matrices(matrixIdx)%ptr
          ENDDO !equationsMatrixIdx
        ELSE
          IF(ALLOCATED(SOLVER_MAPPING%equationsSetToSolverMatricesMap(equationsSetIdx)% &
            & equationsMatricesToSolverMatrixMaps(matrixIdx)%JACOBIAN_TO_SOLVER_MATRIX_MAPS)) THEN
            DO equationsMatrixIdx=1,SOLVER_MAPPING%equationsSetToSolverMatricesMap(equationsSetIdx)% &
              & equationsMatricesToSolverMatrixMaps(matrixIdx)%numberOfEquationsJacobians
              SOLVER_MAPPING%equationsSetToSolverMatricesMap(equationsSetIdx)% &
                & equationsMatricesToSolverMatrixMaps(matrixIdx)%JACOBIAN_TO_SOLVER_MATRIX_MAPS(equationsMatrixIdx)%ptr% &
                & solverMatrix=>solverEquations%solverMatrices%matrices(matrixIdx)%ptr
            ENDDO !equationsMatrixIdx
          ELSE
            DO equationsMatrixIdx=1,SOLVER_MAPPING%equationsSetToSolverMatricesMap(equationsSetIdx)% &
              & equationsMatricesToSolverMatrixMaps(matrixIdx)%numberOfLinearEquationsMatrices
              !Add the solver matrix to the solvers mapping
              SOLVER_MAPPING%equationsSetToSolverMatricesMap(equationsSetIdx)%equationsMatricesToSolverMatrixMaps( &
                & matrixIdx)%linearEquationsMatrixToSolverMatrixMaps(equationsMatrixIdx)%ptr%solverMatrix=> &
                & solverEquations%solverMatrices%matrices(matrixIdx)%ptr
            ENDDO !equationsMatrixIdx
          ENDIF
        ENDIF
      ENDDO !equationsSetIdx
    ENDDO !matrixIdx
    IF(solverEquations%linearity==PROBLEM_SOLVER_NONLINEAR) THEN
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
  SUBROUTINE SolverMatrices_StorageTypeSet(solverMatrices,storageType,err,error,*)

    !Argument variables
    TYPE(SolverMatricesType), POINTER :: solverMatrices !<A pointer to the solver matrices
    INTEGER(INTG), INTENT(IN) :: storageType(:) !<storageType(matrixIdx). The storage type for the matrixIdx'th solver matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(SolverMatrixType), POINTER :: solverMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("SolverMatrices_StorageTypeSet",err,error,*999)

    CALL SolverMatrices_AssertNotFinished(solverMatrices,err,error,*999)
    IF(SIZE(storageType,1)<solverMatrices%numberOfMatrices) THEN
      localError="The size of the solver storage type array of "//TRIM(NumberToVString(SIZE(storageType,1),"*",err,error))// &
        & " is less than the number of solver matrices of "// &
        & TRIM(NumberToVString(solverMatrices%numberOfMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    DO matrixIdx=1,solverMatrices%numberOfMatrices
      NULLIFY(solverMatrix)
      CALL SolverMatrices_SolverMatrixGet(solverMatrices,matrixIdx,solverMatrix,err,error,*999)
      SELECT CASE(storageType(matrixIdx))
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
        localError="The specified storage type of "//TRIM(NumberToVString(storageType(matrixIdx),"*",err,error))// &
          & " for the matrix number "//TRIM(NumberToVString(matrixIdx,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDDO !matrixIdx
   
    EXITS("SolverMatrices_StorageTypeSet")
    RETURN
999 ERRORSEXITS("SolverMatrices_StorageTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMatrices_StorageTypeSet

  !
  !================================================================================================================================
  !

  !>Sets the symmetry type of the solver matrices
  SUBROUTINE SolverMatrices_SymmetryTypeSet(solverMatrices,symmetryTypes,err,error,*)

    !Argument variables
    TYPE(SolverMatricesType), POINTER :: solverMatrices !<A pointer to the solver matrices to set the symmetry types for
    INTEGER(INTG), INTENT(IN) :: symmetryTypes(:) !<symmetryTypes(matrixIdx). The symmetry type for the matrixIdx'th solver matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(SolverMatrixType), POINTER :: solverMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("SolverMatrices_SymmetryTypeSet",err,error,*999)

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
    
    EXITS("SolverMatrices_SymmetryTypeSet")
    RETURN
999 ERRORSEXITS("SolverMatrices_SymmetryTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMatrices_SymmetryTypeSet

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
    TYPE(DistributedMatrixType), POINTER :: equationsDistributedMatrix,solverDistributedMatrix
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatrixToSolverMatrixMapType), POINTER :: equationsToSolverMap
    TYPE(MatrixRowColCouplingType), POINTER :: equationsColToSolverColsMap(:)
    TYPE(MatrixRowColCouplingType), POINTER :: equationsRowToSolverRowsMap(:)
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolverMatricesType), POINTER :: solverMatrices
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("SolverMatrix_EquationsMatrixAdd",err,error,*999)

    IF(.NOT.ASSOCIATED(solverMatrix)) CALL FlagError("Solver matrix is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(equationsMatrix)) CALL FlagError("Equations matrix is not associated.",err,error,*999)
    IF(ABS(alpha)>ZERO_TOLERANCE) THEN
      NULLIFY(solverMatrices)
      CALL SolverMatrix_SolverMatricesGet(solverMatrix,solverMatrices,err,error,*999)
      CALL SolverMatrices_AssertIsFinished(solverMatrices,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverMatrices_SolverMappingGet(solverMatrices,solverMapping,err,error,*999)      
      linearMatrices=>equationsMatrix%linearMatrices
      dynamicMatrices=>equationsMatrix%dynamicMatrices
      IF(.NOT.ASSOCIATED(dynamicMatrices).AND..NOT.ASSOCIATED(linearMatrices)) &
        & CALL FlagError("Equations matrix dynamic or linear matrices is not associated.",err,error,*999)
      NULLIFY(vectorMatrices)
      IF(ASSOCIATED(dynamicMatrices)) THEN
        CALL EquationsMatricesDynamic_VectorMatricesGet(dynamicMatrices,vectorMatrices,err,error,*999)
      ELSE
        CALL EquationsMatricesLinear_VectorMatricesGet(linearMatrices,vectorMatrices,err,error,*999)
      ENDIF
      CALL EquationsMatricesVector_AssertIsFInished(vectorMatrices,err,error,*999)
      NULLIFY(solverDistributedMatrix)
      CALL SolverMatrix_SolverDistributedMatrixGet(solverMatrix,solverDistributedMatrix,err,error,*999)
      NULLIFY(equationsDistributedMatrix)
      CALL EquationsMatrix_DistributedMatrixGet(equationsMatrix,equationsDistributedMatrix,err,error,*999)
      NULLIFY(equationsToSolverMap)
      CALL SolverMapping_EquationsToSolverMapGet(solverMapping,equationsSetIdx,equationsMatrix%matrixNumber, &
        & solverMatrix%matrixNumber,equationsToSolverMap,err,error,*999)
      NULLIFY(equationsRowToSolverRowsMap)
      CALL SolverMapping_EquationsRowToSolverRowsMapGet(solverMapping,equationsSetIdx,equationsRowToSolverRowsMap,err,error,*999)
      NULLIFY(equationsColToSolverColsMap)
      CALL SolverMappingEquationsToSolverMap_EquationsColToSolverColsMapGet(equationsToSolverMap,equationsColToSolverColsMap, &
        & err,error,*999)
      CALL DistributedMatrix_MatrixCoupleAdd(solverDistributedMatrix,equationsRowToSolverRowsMap, &
        & equationsColToSolverColsMap,alpha,equationsDistributedMatrix,.FALSE.,err,error,*999)                              
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
    INTEGER(INTG), POINTER :: columnIndices(:),rowIndices(:)
    REAL(DP), POINTER :: interfaceMatrixData(:)
    TYPE(DistributedMatrixType), POINTER :: interfaceDistributedMatrix,solverDistributedMatrix
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices
    TYPE(InterfaceMatrixToSolverMatrixMapType), POINTER :: interfaceToSolverMap
    TYPE(MatrixRowColCouplingType), POINTER :: interfaceColToSolverColsMap(:)
    TYPE(MatrixRowColCouplingType), POINTER :: interfaceRowToSolverColsMap(:)
    TYPE(MatrixRowColCouplingType), POINTER :: interfaceRowToSolverColsMap(:)
    TYPE(MatrixRowColCouplingType), POINTER :: interfaceColToSolverRowsMap(:)
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolverMatricesType), POINTER :: solverMatrices
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("SolverMatrix_InterfaceMatrixAdd",err,error,*999)

    IF(ASSOCIATED(solverMatrix)) CALL FlagError("Solver matrix is not associated.",err,error,*999)
    IF(ASSOCIATED(interfaceMatrix)) CALL FlagError("Interface matrix is not associated.",err,error,*999)
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
      NULLIFY(interfaceToSolverMap)
      CALL SolverMapping_InterfaceToSolverMap(solverMapping,interfaceConditionIdx,interfaceMatrix%matrixNumber, &
        & solverMatrix%matrixNumber,interfaceToSolverMap,err,error,*999)
      NULLIFY(interfaceRowToSolverRowsMap)
      CALL SolverMapping_InterfaceRowToSolverRowsMapGet(solverMapping,interfaceConditionIdx,interfaceMatrix%matrixNumber, &
        & interfaceRowToSolverRowsMap,err,error,*999)
      NULLIFY(interfaceColToSolverColsMap)
      CALL SolverMapping_InterfaceColToSolverColsMapGet(solverMapping,interfaceConditionIdx,solverMatrix%matrixNumber, &
        & interfaceColToSolverColsMap,err,error,*999)
      CALL DistributedMatrix_MatrixCoupleAdd(solverDistributedMatrix,interfaceRowToSolverRowsMap, &
        & interfaceColToSolverColsMap,alpha(1),interfaceDistributedMatrix,.FALSE.,err,error,*999)                           
      IF(interfaceMatrix%hasTranspose) THEN
        IF(ABS(alpha(2))>ZERO_TOLERANCE) THEN
          NULLIFY(interfaceDistributedMatrix)
          CALL InterfaceMatrix_TransposeDistributeMatrixGet(interfaceMatrix,interfaceDistributedMatrix,err,error,*999)
          NULLIFY(interfaceColToSolverRowsMap)
          CALL SolverMapping_InterfaceColToSolverRowsMap(solverMapping,interfaceConditionIdx,interfaceColToSolverRowsMap, &
            & err,error,*999)
          NULLIFY(interfaceRowToSolverColsMap)
          CALL SolverMappingInterfaceToSolverMap_InterfaceRowToSolverColsMap(interfaceToSolverMap,interfaceRowToSolverColsMap, &
            & err,error,*999)
          CALL DistributedMatrix_MatrixCoupleAdd(solverDistributedMatrix,interfaceColToSolverRowsMap,interfaceRowToSolverColsMap, &
            & alpha(2),interfaceDistributedMatrix,.FALSE.,err,error,*999)
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
    REAL(DP), POINTER :: jacobianMatrixData(:)
    TYPE(DistributedMatrixType), POINTER :: jacobianDistributedMatrix,solverDistributedMatrix
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(JacobianToSolverMapType), POINTER :: jacobianToSolverMap
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolverMatricesType), POINTER :: solverMatrices
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("SolverMatrix_JacobianMatrixAdd",err,error,*999)

    IF(ASSOCIATED(solverMatrix)) CALL FlagError("Solver matrix is not associated.",err,error,*999)
    IF(ASSOCIATED(jacobianMatrix)) CALL FlagError("Jacobian matrix is not associated.",err,error,*999)
    IF(ABS(alpha)>ZERO_TOLERANCE) THEN
      NULLIFY(solverMatrices)
      CALL SolverMatrix_SolverMatricesGet(solverMatrix,solverMatrices,err,error,*999)
      CALL SolverMatrices_AssertIsFinished(solverMatrices,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverMatrices_SolverMappingGet(solverMatrices,solverMapping,err,error,*999)
      NULLIFY(residualVector)
      CALL JacobianMatrix_ResidualVectorGet(jacobianMatrix,residualVector,err,error,*999)
      NULLIFY(nonlinearMatrices)
      CALL EquationsMatricesResidual_NonlinearMatricesGet(residualVector,nonlinearMatrices,err,error,*999)
      NULLIFY(vectorMatrices)
      CALL EquationsMatricesNonlinear_VectorMatricesGet(nonlinearMatrices,vectorMatrices,err,error,*999)
      CALL EquationsMatricesVector_AssertIsFinished(vectorMatrices,err,error,*999)
      NULLIFY(jacobianToSolverMap)
      CALL SolverMapping_JacobianToSolverMapGet(solverMapping,equationsSetIdx,jacobianMatrix%jacobianNumber, &
        & solverMatrix%matrixNumber,jacobianToSolverMap,err,error,*999)
      NULLIFY(solverDistributedMatrix)
      CALL SolverMatrix_SolverDistributedMatrixGet(solverMatrix,solverDistributedMatrix,err,error,*999)
      NULLIFY(jacobianDistributedMatrix)
      CALL JacobianMatrix_DistributedMatrixGet(jacobianMatrix,jacobianDistributedMatrix,err,error,*999)
      CALL DistributedMatrix_MatrixCoupleAdd(solverDistributedMatrix,solverMapping% &
        & equationsSetToSolverMatricesMap(equationsSetIdx)%equationsRowToSolverRowsMap, &
        & jacobianToSolverMap%jacobianColToSolverColsMap,alpha,jacobianDistributedMatrix, &
        & .FALSE.,err,error,*999)
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
    INTEGER(INTG)  :: dummyErr,equationsMatrixIdx,equationsSetIdx,interfaceConditionIdx,interfaceMatrixIdx, &
      & maximumColumnIndices,maximumColumnsPerRow,maximumTransposeColumnsPerRow, &
      & numberOfColumns,solverColumnIdx,solverMatrixIdx,solverRowNumber
    INTEGER(INTG), ALLOCATABLE :: columns(:)
    REAL(DP) :: sparsity
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix,solverDistributedMatrix
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatrixToSolverMatrixMapType), POINTER :: equationsToSolverMap
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition
    TYPE(InterfaceMatrixType), POINTER :: interfaceMatrix
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices    
    TYPE(InterfaceMatrixToSolverMatrixMapType), POINTER :: interfaceToSolverMap
    TYPE(JacobianToSolverMapType), POINTER :: jacobianToSolverMap
    TYPE(ListPtrType), ALLOCATABLE :: columnIndicesLists(:)
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolverMatricesType), POINTER :: solverMatrices
    TYPE(VARYING_STRING) :: dummyError,localError
    
    ENTERS("SolverMatrix_StructureCalculate",err,error,*999)

    numberOfNonZeros=0
    IF(ASSOCIATED(rowIndices)) CALL FlagError("Row indices is already associated.",err,error,*998)
    IF(ASSOCIATED(columnIndices)) CALL FlagError("Row indices is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverMatrix)) CALL FlagError("Solver matrix is not associated.",err,error,*999)
    NULLIFY(solverDistributedMatrix)
    CALL SolverMatrix_SolverDistributedMatrixGet(solverMatrix,solverDistributedMatrix,err,error,*999)
    CALL DistributedMatrix_AssertNotFinished(solverDistributedMatrix,err,error,*999)
    NULLIFY(solverMatrices)
    CALL SolverMatix_SolverMatricesGet(solverMatrix,solverMatrices,err,error,*999)
    NULLIFY(solverMapping)
    CALL SolverMatrices_SolverMappingGet(solverMatrices,solverMapping,err,error,*999)
    SELECT CASE(solverMatrix%storageType)
    CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
      CALL FlagError("Can not calculate the structure for a block storage matrix.",err,error,*999)
    CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
      CALL FlagError("Can not calcualte the structure for a diagonal matrix.",err,error,*999)
    CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
      solverMatrixIdx=solverMatrix%matrixNumber
      !Find the maximum number of column indices
      maximumColumnIndices=0
      DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
        !Loop over dynamic matrices mapped to the solver matrix
        CALL SolverMapping_SolverNumberOfDynamicMatricesGet(solverMapping,solverMatrixIdx,equationsSetIdx, &
          & numberOfDynamicMatrices,err,error,*999)
        DO equationsMatrixIdx=1,numberOfDynamicMatrices
          NULLIFY(equationsMatrix)
          CALL SolverMapping_SolverDynamicMatrixGet(solverMapping,solverMatrixIdx,equationsSetIdx,equationsMatrixIdx, &
            & equationsMatrix,err,error,*999)
          NULLIFY(distributedMatrix)
          CALL EquationsMatrix_DistributedMatrixGet(equationsMatrix,distributedMatrix,err,error,*999)
          CALL DistributedMatrix_MaxColumnsPerRowGet(distributedMatrix,maximumColumnsPerRow,err,error,*999)
          maximumColumnIndices=maximumColumnIndices+maximumColumnsPerRow
        ENDDO !equationsMatrixIdx
        !Loop over linear matrices mapped to the solver matrix
        CALL SolverMapping_SolverNumberOfLinearMatricesGet(solverMapping,solverMatrixIdx,equationsSetIdx, &
          & numberOfLinearMatrices,err,error,*999)
        DO equationsMatrixIdx=1,numberOfLinearMatrices
          NULLIFY(equationsMatrix)
          CALL SolverMapping_SolverLinearMatrixGet(solverMapping,solverMatrixIdx,equationsSetIdx,equationsMatrixIdx, &
            & equationsMatrix,err,error,*999)
          NULLIFY(distributedMatrix)
          CALL EquationsMatrix_DistributedMatrixGet(equationsMatrix,distributedMatrix,err,error,*999)
          CALL DistributedMatrix_MaxColumnsPerRowGet(distributedMatrix,maximumColumnsPerRow,err,error,*999)
          maximumColumnIndices=maximumColumnIndices+maximumColumnsPerRow
        ENDDO !equationsMatrixIdx
        !Loop over Jacobian matrices mapped to the solver matrix
        CALL SolverMappping_SolverNumberOfJacobianMatricesGet(solverMapping,solverMatrixIdx,equationsSetIdx, &
          & numberOfJacobianMatrices,err,error,*999)
        DO equationsMatrixIdx=1,numberOfJacobianMatrices
          NULLIFY(jacobianMatrix)
          CALL SolverMapping_SolverJacobianMatrixGet(solverMapping,solverMatrixIdx,equationsSetIdx,equationsMatrixIdx, &
            & jacobianMatrix,err,error,*999)
          NULLIFY(distributedMatrix)
          CALL JacobianMatrix_DistributedMatrixGet(jacobianMatrix,distributedMatrixGet,err,error,*999)
          CALL DistributedMatrix_MaxColumnsPerRowGet(distributedMatrix,maximumColumnsPerRow,err,error,*999)
          maximumColumnIndices=maximumColumnIndices+maximumColumnsPerRow
        ENDDO !equationsMatrixIdx
      ENDDO !equationsSetIdx
      !Loop over any interface conditions
      DO interfaceConditionIdx=1,solverMapping%numberOfInterfaceConditions
        NULLIFY(interfaceCondition)
        CALL SolverMapping_InterfaceConditionGet(solverMapping,interfaceConditionIdx,interfaceCondition,err,error,*999)
        SELECT CASE(interfaceCondition%method)
        CASE(interfaceCondition_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
          CALL SolverMapping_SolverNumberOfInterfaceMatricesGet(solverMapping,solverMatrixIdx,interfaceConditionIdx, &
            & numberOfInterfaceMatrices,err,error,*999)
          DO interfaceMatrixIdx=1,numberOfInterfaceMatrices
            NULLIFY(interfaceMatrix)
            CALL SolverMapping_SolverInterfaceMatrixGet(solverMapping,solverMatrixIdx,interfaceConditionIdx,interfaceMatrixIdx, &
              & interfaceMatrix,err,error,*999)
            NULLIFY(distributedMatrix)
            CALL InterfaceMatrix_DistributedMatrixGet(interfaceMatrix,distributedMatrix,err,error,*999)
            CALL DistributedMatrix_MaxColumnsPerRowGet(distributedMatrix,maximumColumnsPerRow,err,error,*999)
            maximumTransposeColumnsPerRow=0
            IF(interfaceMatrix%hasTranspose) THEN
              NULLIFY(distributedMatrix)
              CALL InterfaceMatrix_TransposeDistributedMatrixGet(interfaceMatrix,distributedMatrix,err,error,*999)
              CALL DistributedMatrix_MaxColumnsPerRowGet(distributedMatrix,maximumTransposeColumnsPerRow,err,error,*999)
            ENDIF
            maximumColumnIndices=maximumColumnIndices+MAX(maximumColumnsPerRow,maximumTransposeColumnsPerRow)
          ENDDO !interfaceMatrixIdx
        CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE DEFAULT
          localError="The interface condition method of "//TRIM(NumberToVString(interfaceCondition%method,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDDO !interfaceConditionIdx
      !Allocate lists
      ALLOCATE(columnIndicesLists(solverMapping%numberOfRows),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate column indices lists.",err,error,*999)
      !Allocate row indices
      ALLOCATE(rowIndices(solverMapping%numberOfRows+1),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate row indices.",err,error,*999)
      rowIndices(1)=1
      !Set up the column indicies lists
      DO solverRowNumber=1,solverMapping%numberOfRows
        NULLIFY(columnIndicesLists(solverRowNumber)%ptr)
        CALL List_CreateStart(columnIndicesLists(solverRowNumber)%ptr,err,error,*999)
        CALL List_DataTypeSet(columnIndicesLists(solverRowNumber)%ptr,LIST_INTG_TYPE,err,error,*999)
        CALL List_InitialSizeSet(columnIndicesLists(solverRowNumber)%ptr,maximumColumnIndices,err,error,*999)
        CALL List_CreateFinish(columnIndicesLists(solverRowNumber)%ptr,err,error,*999)
      ENDDO !solverRowNumber
      !Loop over the equations sets
      DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
        !Loop over the dynamic equations matrices mapped to the solver matrix and calculate the col indices by row.
        CALL SolverMapping_SolverNumberOfDynamicMatricesGet(solverMapping,solverMatrixIdx,equationsSetIdx, &
          & numberOfDynamicMatrices,err,error,*999)
        DO equationsMatrixIdx=1,numberOfDynamicMatrices
          NULLIFY(equationsMatrix)
          CALL SolverMapping_SolverDynamicMatrixGet(solverMapping,solverMatrixIdx,equationSetIdx,equationsMatrixIdx, &
            & equationsMatrix,err,error,*999)
          NULLIFY(distributedMatrix)
          CALL EquationsMatrix_DistributedMatrixGet(equationsMatrix,distributedMatrix,err,error,*999)
          NULLIFY(equationsToSolverMap)
          !Note: pointers have been checked above
          equationsToSolverMap=>solverMapping%equationsSetToSolverMatricesMap(equationsSetIdx)% &
            & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%dynamicEquationsMatrixToSolverMatrixMaps( &
            & equationsMatrixIdx)%ptr
          equationsMatrix=>equationsToSolverMap%equationsMatrix
          dynamicMatrices=>equationsMatrix%dynamicMatrices
          vectorMatrices=>dynamicMatrices%vectorMatrices
          distributedMatrix=>equationsMatrix%matrix
          
          CALL DistributedMatrix_MatrixStructureCoupleCalculate(distributedMatrix,.FALSE.,solverMapping% &
            & equationsSetToSolverMatricesMap(equationsSetIdx)%equationsRowToSolverRowsMap, &
            & equationsToSolverMap%equationsColToSolverColsMap,columnIndicesLists,err,error,*999)
          
        ENDDO !equationsMatrixIdx
        !Loop over the linear equations matrices mapped to the solver matrix and calculate the col indices by row.
        DO equationsMatrixIdx=1,solverMapping%equationsSetToSolverMatricesMap(equationsSetIdx)% &
          & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%numberOfLinearEquationsMatrices
          !Note: pointers have been checked above
          equationsToSolverMap=>solverMapping%equationsSetToSolverMatricesMap(equationsSetIdx)% &
            & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%linearEquationsMatrixToSolverMatrixMaps( &
            & equationsMatrixIdx)%ptr
          equationsMatrix=>equationsToSolverMap%equationsMatrix
          linearMatrices=>equationsMatrix%linearMatrices
          vectorMatrices=>linearMatrices%vectorMatrices
          distributedMatrix=>equationsMatrix%matrix
          
          CALL DistributedMatrix_MatrixStructureCoupleCalculate(distributedMatrix,.FALSE.,solverMapping% &
            & equationsSetToSolverMatricesMap(equationsSetIdx)%equationsRowToSolverRowsMap, &
            & equationsToSolverMap%equationsColToSolverColsMap,columnIndicesLists,err,error,*999)
          
        ENDDO !equationsMatrixIdx
        !Now add any columns from the Jacobians
        DO equationsMatrixIdx=1,solverMapping%equationsSetToSolverMatricesMap(equationsSetIdx)% &
          & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%numberOfEquationsJacobians
          jacobianToSolverMap=>solverMapping%equationsSetToSolverMatricesMap(equationsSetIdx)% &
            & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%JACOBIAN_TO_SOLVER_MATRIX_MAPS( &
            & equationsMatrixIdx)%ptr
          IF(ASSOCIATED(jacobianToSolverMap)) THEN
            !Note: pointers have been checked above
            jacobianMatrix=>jacobianToSolverMap%jacobianMatrix
            nonlinearMatrices=>jacobianMatrix%nonlinearMatrices
            vectorMatrices=>nonlinearMatrices%vectorMatrices
            distributedMatrix=>jacobianMatrix%JACOBIAN
            
            CALL DistributedMatrix_MatrixStructureCoupleCalculate(distributedMatrix,.FALSE.,solverMapping% &
              & equationsSetToSolverMatricesMap(equationsSetIdx)%equationsRowToSolverRowsMap, &
              & jacobianToSolverMap%jacobianColToSolverColsMap,columnIndicesLists,err,error,*999)
            
          ENDIF
        ENDDO !equationsMatrixIdx
        !Now add in any interface matrices columns
        DO interfaceConditionIdx=1,solverMapping%equationsSetToSolverMatricesMap(equationsSetIdx)% &
          & numberOfInterfaceConditions
        ENDDO !interfaceConditionIdx
      ENDDO !equationsSetIdx
      !Loop over any equations sets
      DO interfaceConditionIdx=1,solverMapping%numberOfInterfaceConditions
        INTERFACE_CONDITION=>solverMapping%interfaceConditions(interfaceConditionIdx)%ptr
        SELECT CASE(interfaceCondition%METHOD)
        CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
          DO interfaceMatrixIdx=1,solverMapping%interfaceConditionToSolverMatricesMap(interfaceConditionIdx)% &
            & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%numberOfInterfaceMatrices
            interfaceToSolverMap=>solverMapping%interfaceConditionToSolverMatricesMap(interfaceConditionIdx)% &
              & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%interfaceMatrixToSolverMatrixMaps( &
              & interfaceMatrixIdx)%ptr
            interfaceMatrix=>interfaceToSolverMap%interfaceMatrix
            interfaceMatrices=>interfaceMatrix%interfaceMatrices
            distributedMatrix=>interfaceMatrix%matrix
            
            CALL DistributedMatrix_MatrixStructureCoupleCalculate(distributedMatrix,.FALSE.,solverMapping% &
              & interfaceConditionToSolverMatricesMap(interfaceConditionIdx)%interfaceMatrixToSolverMatricesMaps( &
              & interfaceMatrixIdx)%interfaceRowToSolverRowsMap,solverMapping%interfaceConditionToSolverMatricesMap( &
              & interfaceConditionIdx)%interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)% &
              & interfaceColToSolverColsMap,columnIndicesLists,err,error,*999)
            
            IF(interfaceMatrix%hasTranspose) THEN
              distributedMatrix=>interfaceMatrix%matrixTranspose
              
              CALL DistributedMatrix_MatrixStructureCoupleCalculate(distributedMatrix,.FALSE.,solverMapping% &
                & interfaceConditionToSolverMatricesMap(interfaceConditionIdx)%interfaceColToSolverRowsMap, &
                & interfaceToSolverMap%interfaceRowToSolverColsMap,columnIndicesLists,err,error,*999)
              
            ENDIF
          ENDDO !interfaceMatrixIdx
        CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE DEFAULT
          localError="The interface condition method of "// &
            & TRIM(NumberToVString(interfaceCondition%method,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDDO !interfaceConditionIdx
      !Loop over the rows to calculate the number of non-zeros and setup the row indicces
      DO solverRowNumber=1,solverMapping%numberOfRows
        CALL List_RemoveDuplicates(columnIndicesLists(solverRowNumber)%ptr,err,error,*999)
        CALL List_NumberOfItemsGet(columnIndicesLists(solverRowNumber)%ptr,numberOfColumns,err,error,*999)
        numberOfNonZeros=numberOfNonZeros+numberOfColumns
        rowIndices(solverRowNumber+1)=numberOfNonZeros+1
      ENDDO !solverRowNumber
      !Allocate and setup the column locations
      ALLOCATE(columnIndices(numberOfNonZeros),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate column indices.",err,error,*999)
      DO solverRowNumber=1,solverMapping%numberOfRows
        CALL List_DetachAndDestroy(columnIndicesLists(solverRowNumber)%ptr,numberOfColumns,columns, &
          & err,error,*999)
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
    
    IF(DIAGNOSTICS1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Solver matrix structure:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Solver matrix number : ",solverMatrix%matrixNumber, &
        & err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of rows = ",solverMatrices%numberOfRows, &
        & err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of columns = ",solverMatrix%numberOfColumns, &
        & err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of non zeros = ",numberOfNonZeros,err,error,*999)
      IF(solverMatrices%numberOfRows*solverMatrix%numberOfColumns/=0) THEN
        sparsity=REAL(numberOfNonZeros,DP)/REAL(solverMatrices%numberOfRows* &
          & solverMatrix%numberOfColumns,DP)*100.0_DP
        CALL WriteStringFmtValue(DIAGNOSTIC_OUTPUT_TYPE,"  Sparsity (%) = ",sparsity,"F6.2", err,error,*999)
      ENDIF
      IF(DIAGNOSTICS2) THEN
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMatrices%numberOfRows+1,8,8,rowIndices, &
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
      DO solverRowNumber=1,solverMapping%numberOfRows
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
  SUBROUTINE SOLVER_MATRIX_FINALISE(solverMatrix,err,error,*)

    !Argument variables
    TYPE(SolverMatrixType), POINTER :: solverMatrix !<A pointer to the solver matrix to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("SOLVER_MATRIX_FINALISE",err,error,*999)

    IF(ASSOCIATED(solverMatrix)) THEN
      IF(ASSOCIATED(solverMatrix%matrix)) CALL DistributedMatrix_Destroy(solverMatrix%matrix,err,error,*999)
      IF(ASSOCIATED(solverMatrix%solverVector)) CALL DistributedVector_Destroy(solverMatrix%solverVector,err,error,*999)
      DEALLOCATE(solverMatrix)
    ENDIF
    
    EXITS("SOLVER_MATRIX_FINALISE")
    RETURN
999 ERRORSEXITS("SOLVER_MATRIX_FINALISE",err,error)    
    RETURN 1
   
  END SUBROUTINE SOLVER_MATRIX_FINALISE
        
  !
  !================================================================================================================================
  !

  !>Forms a solver matrix by initialising the structure of the matrix to zero.
  SUBROUTINE SOLVER_MATRIX_FORM(solverMatrix,err,error,*)

    !Argument variables
    TYPE(SolverMatrixType), POINTER :: solverMatrix !<A pointer to the solver matrix to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("SOLVER_MATRIX_FORM",err,error,*999)

    IF(ASSOCIATED(solverMatrix)) THEN
      CALL DistributedMatrix_Form(solverMatrix%matrix,err,error,*999)
    ELSE
      CALL FlagError("Solver matrix is not associated.",err,error,*999)
    ENDIF
    
    EXITS("SOLVER_MATRIX_FORM")
    RETURN
999 ERRORSEXITS("SOLVER_MATRIX_FORM",err,error)    
    RETURN 1
   
  END SUBROUTINE SOLVER_MATRIX_FORM
        
    !
  !================================================================================================================================
  !

  !>Initialises a solver matrix
  SUBROUTINE SOLVER_MATRIX_INITIALISE(solverMatrices,MATRIX_NUMBER,err,error,*)

    !Argument variables
    TYPE(SolverMatricesType), POINTER :: solverMatrices !<A pointer to the solver matrices to initialise
    INTEGER(INTG), INTENT(IN) :: MATRIX_NUMBER !<The matrix number in the solver matrices to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolverMatrixType), POINTER :: solverMatrix
    TYPE(VARYING_STRING) :: dummyError,localError
    
    ENTERS("SOLVER_MATRIX_INITIALISE",err,error,*998)

    IF(ASSOCIATED(solverMatrices)) THEN
      IF(MATRIX_NUMBER>0.AND.MATRIX_NUMBER<=solverMatrices%numberOfMatrices) THEN
        solverMapping=>solverMatrices%solverMapping
        IF(ASSOCIATED(SOLVER_MAPPING)) THEN
          IF(ASSOCIATED(solverMatrices%matrices(MATRIX_NUMBER)%ptr)) THEN
            CALL FlagError("Solver matrix is already associated.",err,error,*998)
          ELSE
            ALLOCATE(solverMatrices%matrices(MATRIX_NUMBER)%ptr,STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate solver matrix.",err,error,*999)
            solverMatrix=>solverMatrices%matrices(MATRIX_NUMBER)%ptr
            solverMatrix%matrixNumber=MATRIX_NUMBER
            solverMatrix%solverMatrices=>solverMatrices
            solverMatrix%storageType=MATRIX_BLOCK_STORAGE_TYPE
            solverMatrix%updateMatrix=.TRUE.
            solverMatrix%numberOfColumns=SOLVER_MAPPING%solverColToEquationsColsMap(MATRIX_NUMBER)%numberOfColumns
            SOLVER_MAPPING%solverColToEquationsColsMap(MATRIX_NUMBER)%solverMatrix=>solverMatrix
            NULLIFY(solverMatrix%solverVector)
            NULLIFY(solverMatrix%matrix)
          ENDIF
        ELSE
          CALL FlagError("Solver mapping is not associated.",err,error,*998)
        ENDIF
      ELSE
        localError="The specified matrix number of "//TRIM(NumberToVString(MATRIX_NUMBER,"*",err,error))// &
          & " is invalid. The number must be > 0 and <= "// &
          & TRIM(NumberToVString(solverMatrices%numberOfMatrices,"*",err,error))//"."
        CALL FlagError(localError,err,error,*998)
      ENDIF
    ELSE
      CALL FlagError("Solver matrices is not associated.",err,error,*998)
    ENDIF
    
    EXITS("SOLVER_MATRIX_INITIALISE")
    RETURN
999 CALL SOLVER_MATRIX_FINALISE(solverMatrices%matrices(MATRIX_NUMBER)%ptr,dummyErr,dummyError,*998)
998 ERRORSEXITS("SOLVER_MATRIX_INITIALISE",err,error)    
    RETURN 1
   
  END SUBROUTINE SOLVER_MATRIX_INITIALISE
        
  !
  !================================================================================================================================
  !

END MODULE SolverMatricesRoutines
