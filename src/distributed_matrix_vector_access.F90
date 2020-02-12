!> \file
!> \author Chris Bradley
!> \brief This module contains all distributed matrix vector access method routines.
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

!> This module contains all distributed matrix vector access method routines.
MODULE DistributedMatrixVectorAccessRoutines
  
  USE BaseRoutines
  USE Constants
  USE Kinds
  USE ISO_VARYING_STRING
  USE MatrixVectorAccessRoutines
  USE Strings
  USE Types
  
#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup DistributedMatrixVector_Constants OpenCMISS::Iron::DistributedMatrixVector::Constants
  !> \brief Distributed matrix vector constants.
  !>@{
  !> \addtogroup DistributedMatrixVector_LibraryTypes DistributedMatrixVector::LibraryTypes
  !> \brief Distributed matrix-vector library types
  !> \see DistributedMatrixVector
  !>@{
  INTEGER(INTG), PARAMETER :: DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE=LIBRARY_CMISS_TYPE !<CMISS distributed matrix-vector library type \see DistributedMatrixVector_LibraryTypes,DistributedMatrixVector
  INTEGER(INTG), PARAMETER :: DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE=LIBRARY_PETSC_TYPE !<PETSc distributed matrix-vector library type \see DistributedMatrixVector_LibraryTypes,DistributedMatrixVector
  !>@}
  
  !> \addtogroup DistributedMatrixVector_DataTypes DistributedMatrixVector::DataTypes
  !> \brief Distributed matrix-vector data types
  !> \see DistributedMatrixVector
  !>@{
  INTEGER(INTG), PARAMETER :: DISTRIBUTED_MATRIX_VECTOR_INTG_TYPE=MATRIX_VECTOR_INTG_TYPE !<Integer distributed matrix-vector data type \see DistributedMatrixVector_DataTypes,DistributedMatrixVector
  INTEGER(INTG), PARAMETER :: DISTRIBUTED_MATRIX_VECTOR_SP_TYPE=MATRIX_VECTOR_SP_TYPE !<Single precision real distributed matrix-vector data type \see DistributedMatrixVector_DataTypes,DistributedMatrixVector
  INTEGER(INTG), PARAMETER :: DISTRIBUTED_MATRIX_VECTOR_DP_TYPE=MATRIX_VECTOR_DP_TYPE !<Double precision real distributed matrix-vector data type \see DistributedMatrixVector_DataTypes,DistributedMatrixVector
  INTEGER(INTG), PARAMETER :: DISTRIBUTED_MATRIX_VECTOR_L_TYPE=MATRIX_VECTOR_L_TYPE !<Logical distributed matrix-vector data type \see DistributedMatrixVector_DataTypes,DistributedMatrixVector
  INTEGER(INTG), PARAMETER :: DISTRIBUTED_MATRIX_VECTOR_SPC_TYPE=MATRIX_VECTOR_SPC_TYPE !<Single precision complex distributed matrix-vector data type \see DistributedMatrixVector_DataTypes,DistributedMatrixVector
  INTEGER(INTG), PARAMETER :: DISTRIBUTED_MATRIX_VECTOR_DPC_TYPE=MATRIX_VECTOR_DPC_TYPE !<Double precision complex distributed matrix-vector data type \see DistributedMatrixVector_DataTypes,DistributedMatrixVector
  !>@}

  !> \addtogroup DistributedMatrixVector_StorageTypes DistributedMatrixVector::StorageTypes
  !> \brief Distributed matrix-vector storage type parameters
  !>@{
  INTEGER(INTG), PARAMETER :: DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE=MATRIX_BLOCK_STORAGE_TYPE !<Distributed matrix block storage type \see DistributedMatrixVector_StorageTypes,DistributedMatrixVector
  INTEGER(INTG), PARAMETER :: DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE=MATRIX_DIAGONAL_STORAGE_TYPE !<Distributed matrix diagonal storage type \see DistributedMatrixVector_StorageTypes,DistributedMatrixVector
  INTEGER(INTG), PARAMETER :: DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE=MATRIX_COLUMN_MAJOR_STORAGE_TYPE !<Distributed matrix column major storage type \see DistributedMatrixVector_StorageTypes,DistributedMatrixVector
  INTEGER(INTG), PARAMETER :: DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE=MATRIX_ROW_MAJOR_STORAGE_TYPE !<Distributed matrix row major storage type \see DistributedMatrixVector_StorageTypes,DistributedMatrixVector
  INTEGER(INTG), PARAMETER :: DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE=MATRIX_COMPRESSED_ROW_STORAGE_TYPE !<Distributed matrix compressed row storage type \see DistributedMatrixVector_StorageTypes,DistributedMatrixVector
  INTEGER(INTG), PARAMETER :: DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE=MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE !<Distributed matrix compressed column storage type \see DistributedMatrixVector_StorageTypes,DistributedMatrixVector
  INTEGER(INTG), PARAMETER :: DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE=MATRIX_ROW_COLUMN_STORAGE_TYPE !<Distributed matrix row-column storage type \see DistributedMatrixVector_StorageTypes,MATRIX_VECTOR
  INTEGER(INTG), PARAMETER :: DISTRIBUTED_MATRIX_BLOCK_COMPRESSED_ROW_STORAGE_TYPE=MATRIX_BLOCK_COMPRESSED_ROW_STORAGE_TYPE !<Distributed matrix block compressed row storage type \see DistributedMatrixVector_StorageTypes,DistributedMatrixVector
  INTEGER(INTG), PARAMETER :: DISTRIBUTED_MATRIX_MODIFIED_KRM_STORAGE_TYPE=MATRIX_MODIFIED_KRM_STORAGE_TYPE !<Distributed matrix modified Knuth-Rheinbolt-Mesztenyi storage type \see DistributedMatrixVector_StorageTypes,DistributedMatrixVector
  !>@}
  
  !> \addtogroup DistributedMatrixVector_SymmetryTypes DistributedMatrixVector::SymmetryTypes
  !> \brief Distributed matrix symmetry type parameters
  !> \see DistributedMatrixVector_Symmetry,DistributedMatrixVector
  !>@{
  INTEGER(INTG), PARAMETER :: DISTRIBUTED_MATRIX_SYMMETRIC_TYPE=MATRIX_SYMMETRIC_TYPE !<Distributed matrix is symmetric \see DistributedMatrixVector_SymmetryTypes,DistributedMatrixVector
  INTEGER(INTG), PARAMETER :: DISTRIBUTED_MATRIX_HERMITIAN_TYPE=MATRIX_HERMITIAN_TYPE !<Distributed matrix is Hermitian (complex symmetric) \see DistributedMatrixVector_SymmetryTypes,DistributedMatrixVector
  INTEGER(INTG), PARAMETER :: DISTRIBUTED_MATRIX_SKEW_SYMMETRIC_TYPE=MATRIX_SKEW_SYMMETRIC_TYPE !<Distributed matrix is skew-symmetric \see DistributedMatrixVector_SymmetryTypes,DistributedMatrixVector
  INTEGER(INTG), PARAMETER :: DISTRIBUTED_MATRIX_UNSYMMETRIC_TYPE=MATRIX_UNSYMMETRIC_TYPE !<Distributed matrix is unsymmetric \see DistributedMatrixVector_SymmetryTypes,DistributedMatrixVector
  !>@}

  !> \addtogroup DistributedMatrixVector_GhostingTypes DistributedMatrixVector::GhostingTypes
  !> \brief Distributed matrix-vector ghosting types
  !> \see DistributedMatrixVector
  !>@{
  INTEGER(INTG), PARAMETER :: DISTRIBUTED_MATRIX_VECTOR_INCLUDE_GHOSTS_TYPE=1 !<Include ghost values in the distributed matrix/vector \see DistributedMatrixVector_GhostingTypes,DistributedMatrixVector
  INTEGER(INTG), PARAMETER :: DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE=2 !<Do not include ghost values/rows in the distributed matrix/vector \see DistributedMatrixVector_GhostingTypes,DistributedMatrixVector
  !>@}
  
  !> \addtogroup DistributedMatrixVector_TransposeTypes DistributedMatrixVector::TransposeTypes
  !> \brief Distributed matrix-vector transpose type parameters
  !> \see DistributedMatrixVector_MatrixTransposeTypes,DistributedatrixVector
  !>@{
  INTEGER(INTG), PARAMETER :: DISTRIBUTED_MATRIX_NO_TRANSPOSE_REQUIRED=MATRIX_NO_TRANSPOSE_REQUIRED !<Distributed matrix will not have a transpose required \see DistributedMatrixVector_TransposeTypes,DistributedMatrixVector
  INTEGER(INTG), PARAMETER :: DISTRIBUTED_MATRIX_PARTIAL_TRANSPOSE_REQUIRED=MATRIX_PARTIAL_TRANSPOSE_REQUIRED !<Distributed matrix will require a partial transpose for selected rows/columns \see DistributedMatrixVector_TransposeTypes,DistributedMatrixVector
  INTEGER(INTG), PARAMETER :: DISTRIBUTED_MATRIX_FULL_TRANSPOSE_REQUIRED=MATRIX_FULL_TRANSPOSE_REQUIRED !<Distributed matrix will require a full transpose  \see DistributedMatrixVector_TransposeTypes,DistributedMatrixVector
  !>@}
  !>@}
  
  !Module types

  !Module variables

  !Interfaces

  INTERFACE DISTRIBUTED_MATRIX_MAX_COLUMNS_PER_ROW_GET
    MODULE PROCEDURE DistributedMatrix_MaxColumnsPerRowGet
  END INTERFACE DISTRIBUTED_MATRIX_MAX_COLUMNS_PER_ROW_GET

  INTERFACE DISTRIBUTED_MATRIX_NUMBER_NON_ZEROS_GET
    MODULE PROCEDURE DistributedMatrix_NumberOfNonZerosGet
  END INTERFACE DISTRIBUTED_MATRIX_NUMBER_NON_ZEROS_GET

  INTERFACE DISTRIBUTED_MATRIX_STORAGE_TYPE_GET
    MODULE PROCEDURE DistributedMatrix_StorageTypeGet
  END INTERFACE DISTRIBUTED_MATRIX_STORAGE_TYPE_GET

  PUBLIC DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE,DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE

  PUBLIC DISTRIBUTED_MATRIX_VECTOR_INTG_TYPE,DISTRIBUTED_MATRIX_VECTOR_SP_TYPE,DISTRIBUTED_MATRIX_VECTOR_DP_TYPE, &
    & DISTRIBUTED_MATRIX_VECTOR_L_TYPE,DISTRIBUTED_MATRIX_VECTOR_SPC_TYPE,DISTRIBUTED_MATRIX_VECTOR_DPC_TYPE

  PUBLIC DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE,DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE, &
    & DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE,DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE, &
    & DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE, &
    & DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE,DISTRIBUTED_MATRIX_BLOCK_COMPRESSED_ROW_STORAGE_TYPE, &
    & DISTRIBUTED_MATRIX_MODIFIED_KRM_STORAGE_TYPE

  PUBLIC DISTRIBUTED_MATRIX_SYMMETRIC_TYPE,DISTRIBUTED_MATRIX_HERMITIAN_TYPE,DISTRIBUTED_MATRIX_SKEW_SYMMETRIC_TYPE, &
    & DISTRIBUTED_MATRIX_UNSYMMETRIC_TYPE

  PUBLIC DISTRIBUTED_MATRIX_VECTOR_INCLUDE_GHOSTS_TYPE,DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE

  PUBLIC DISTRIBUTED_MATRIX_NO_TRANSPOSE_REQUIRED,DISTRIBUTED_MATRIX_PARTIAL_TRANSPOSE_REQUIRED, &
    & DISTRIBUTED_MATRIX_FULL_TRANSPOSE_REQUIRED

  PUBLIC DistributedMatrix_AssertIsFinished,DistributedMatrix_AssertNotFinished

  PUBLIC DistributedMatrix_AssertIsIntgData,DistributedMatrix_AssertIsSPData,DistributedMatrix_AssertIsDPData, &
    & DistributedMatrix_AssertIsLData
  
  PUBLIC DistributedMatrix_BlockSizeGet
  
  PUBLIC DistributedMatrix_CMISSMatrixGet
  
  PUBLIC DistributedMatrix_ColumnMappingGet

  PUBLIC DistributedMatrix_DataTypeGet

  PUBLIC DistributedMatrix_DimensionsGet

  PUBLIC DISTRIBUTED_MATRIX_MAX_COLUMNS_PER_ROW_GET

  PUBLIC DistributedMatrix_MaxColumnsPerRowGet

  PUBLIC DistributedMatrix_MaxNumberOfRowsGet

  PUBLIC DistributedMatrix_NumberOfBlocksGet
  
  PUBLIC DistributedMatrix_NumberOfColumnsGet
  
  PUBLIC DistributedMatrix_NumberOfGlobalRowsGet
  
  PUBLIC DistributedMatrix_NumberOfLocalRowsGet

  PUBLIC DISTRIBUTED_MATRIX_NUMBER_NON_ZEROS_GET
  
  PUBLIC DistributedMatrix_NumberOfNonZerosGet
  
  PUBLIC DistributedMatrix_PETScMatrixGet
    
  PUBLIC DistributedMatrix_RowMappingGet

  PUBLIC DISTRIBUTED_MATRIX_STORAGE_TYPE_GET

  PUBLIC DistributedMatrix_StorageTypeGet

  PUBLIC DistributedMatrix_SymmetryTypeGet

  PUBLIC DistributedMatrix_TransposeTypeGet

  PUBLIC DistributedMatrixCMISS_DistributedMatrixGet

  PUBLIC DistributedMatrixPetsc_DistributedMatrixGet

  PUBLIC DistributedVector_AssertIsFinished,DistributedVector_AssertNotFinished

  PUBLIC DistributedVector_AssertIsIntgData,DistributedVector_AssertIsSPData,DistributedVector_AssertIsDPData, &
    & DistributedVector_AssertIsLData
  
  PUBLIC DistributedVector_CMISSVectorGet

  PUBLIC DistributedVector_DataTypeGet

  PUBLIC DistributedVector_NumberOfGlobalRowsGet
  
  PUBLIC DistributedVector_NumberOfLocalRowsGet
  
  PUBLIC DistributedVector_PETScVectorGet
  
  PUBLIC DistributedVector_RowMappingGet

  PUBLIC DistributedVector_TotalNumberOfLocalRowsGet
  
  PUBLIC DistributedVectorCMISS_DistributedVectorGet
  
  PUBLIC DistributedVectorPETSc_DistributedVectorGet  

CONTAINS

  !
  !================================================================================================================================
  !

  !>Assert that a distributed matrix has been finished
  SUBROUTINE DistributedMatrix_AssertIsFinished(distributedMatrix,err,error,*)

    !Argument Variables
    TYPE(DistributedMatrixType), POINTER, INTENT(IN) :: distributedMatrix !<The distributed matrix to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("DistributedMatrix_AssertIsFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
#endif    

    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("Distributed matrix has not been finished.",err,error,*999)
    
    EXITS("DistributedMatrix_AssertIsFinished")
    RETURN
999 ERRORSEXITS("DistributedMatrix_AssertIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_AssertIsFinished

  !
  !================================================================================================================================
  !

  !>Assert that a distribute matrix has not been finished
  SUBROUTINE DistributedMatrix_AssertNotFinished(distributedMatrix,err,error,*)

    !Argument Variables
    TYPE(DistributedMatrixType), POINTER, INTENT(IN) :: distributedMatrix !<The distributed matrix to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("DistributedMatrix_AssertNotFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
#endif    

    IF(distributedMatrix%matrixFinished) CALL FlagError("Distributed matrix has already been finished.",err,error,*999)
    
    EXITS("DistributedMatrix_AssertNotFinished")
    RETURN
999 ERRORSEXITS("DistributedMatrix_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_AssertNotFinished

  !
  !================================================================================================================================
  !

  !>Assert that a distributed matrix has an integer data type
  SUBROUTINE DistributedMatrix_AssertIsIntgData(distributedMatrix,err,error,*)

    !Argument Variables
    TYPE(DistributedMatrixType), POINTER, INTENT(IN) :: distributedMatrix !<The distributed matrix to assert the integer data type for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("DistributedMatrix_AssertIsIntgData",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
#endif    

    IF(distributedMatrix%dataType/=DISTRIBUTED_MATRIX_VECTOR_INTG_TYPE) THEN
      localError="The distributed matrix data type of "//TRIM(NumberToVString(distributedMatrix%dataType,"*",err,error))// &
        & " does not correspond to the required integer data type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("DistributedMatrix_AssertIsIntgData")
    RETURN
999 ERRORSEXITS("DistributedMatrix_AssertIsIntgData",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_AssertIsIntgData

  !
  !================================================================================================================================
  !

  !>Assert that a distributed matrix has a single precision real data type
  SUBROUTINE DistributedMatrix_AssertIsSPData(distributedMatrix,err,error,*)

    !Argument Variables
    TYPE(DistributedMatrixType), POINTER, INTENT(IN) :: distributedMatrix !<The distributed matrix to assert the single precision real data type for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("DistributedMatrix_AssertIsSPData",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
#endif    

    IF(distributedMatrix%dataType/=DISTRIBUTED_MATRIX_VECTOR_SP_TYPE) THEN
      localError="The distributed matrix data type of "//TRIM(NumberToVString(distributedMatrix%dataType,"*",err,error))// &
        & " does not correspond to the required single precision real data type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("DistributedMatrix_AssertIsSPData")
    RETURN
999 ERRORSEXITS("DistributedMatrix_AssertIsSPData",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_AssertIsSPData

  !
  !================================================================================================================================
  !

  !>Assert that a distributed matrix has a double precision real data type
  SUBROUTINE DistributedMatrix_AssertIsDPData(distributedMatrix,err,error,*)

    !Argument Variables
    TYPE(DistributedMatrixType), POINTER, INTENT(IN) :: distributedMatrix !<The distributed matrix to assert the double precision real data type for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("DistributedMatrix_AssertIsDPData",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
#endif    

    IF(distributedMatrix%dataType/=DISTRIBUTED_MATRIX_VECTOR_DP_TYPE) THEN
      localError="The distributed matrix data type of "//TRIM(NumberToVString(distributedMatrix%dataType,"*",err,error))// &
        & " does not correspond to the required double precision real data type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("DistributedMatrix_AssertIsDPData")
    RETURN
999 ERRORSEXITS("DistributedMatrix_AssertIsDPData",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_AssertIsDPData

  !
  !================================================================================================================================
  !

  !>Assert that a distributed matrix has a logical data type
  SUBROUTINE DistributedMatrix_AssertIsLData(distributedMatrix,err,error,*)

    !Argument Variables
    TYPE(DistributedMatrixType), POINTER, INTENT(IN) :: distributedMatrix !<The distributed matrix to assert the logical data type for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("DistributedMatrix_AssertIsLData",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
#endif    

    IF(distributedMatrix%dataType/=DISTRIBUTED_MATRIX_VECTOR_L_TYPE) THEN
      localError="The distributed matrix data type of "//TRIM(NumberToVString(distributedMatrix%dataType,"*",err,error))// &
        & " does not correspond to the required logical data type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("DistributedMatrix_AssertIsLData")
    RETURN
999 ERRORSEXITS("DistributedMatrix_AssertIsLData",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_AssertIsLData

  !  
  !================================================================================================================================
  !
  
  !>Gets the block size for the distributed matrix
  SUBROUTINE DistributedMatrix_BlockSizeGet(distributedMatrix,blockSize,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix to get the block size for
    INTEGER(INTG), INTENT(OUT) :: blockSize !<On return, the block size in the distributed matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedMatrix_BlockSizeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
#endif    

    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_BlockSizeGet(cmissMatrix%matrix,blockSize,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot get the block size for a PETSc matrix.",err,error,*999)
    CASE DEFAULT
      localError="The distributed matrix library type of "//TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
         
    EXITS("DistributedMatrix_BlockSizeGet")
    RETURN
999 ERRORSEXITS("DistributedMatrix_BlockSizeGet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_BlockSizeGet
  
  !  
  !================================================================================================================================
  !
  
  !>Gets the CMISS matrix for the distributed matrix
  SUBROUTINE DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix to get the CMISS matrix for
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix !<On return, the CMISS matrix for the distributed matrix. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DistributedMatrix_CMISSMatrixGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(cmissMatrix)) CALL FlagError("CMISS matrix is associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
#endif    
    
    cmissMatrix=>distributedMatrix%cmiss

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(cmissMatrix))  CALL FlagError("CMISS matrix is not associated for the distributed matrix.",err,error,*999)
#endif    
     
    EXITS("DistributedMatrix_CMISSMatrixGet")
    RETURN
999 NULLIFY(cmissMatrix)
998 ERRORSEXITS("DistributedMatrix_CMISSMatrixGet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_CMISSMatrixGet
  
  !  
  !================================================================================================================================
  !
  
  !>Gets the column mapping for the distributed matrix
  SUBROUTINE DistributedMatrix_ColumnMappingGet(distributedMatrix,columnMapping,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix to get the column mapping for
    TYPE(DomainMappingType), POINTER :: columnMapping !<On return, the column mapping for the distributed matrix. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DistributedMatrix_ColumnMappingGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(columnMapping)) CALL FlagError("Column mapping is associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
#endif    
    
    columnMapping=>distributedMatrix%columnDomainMapping

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(columnMapping)) CALL FlagError("Column mapping is not associated for the distributed matrix.",err,error,*999)
#endif    
     
    EXITS("DistributedMatrix_ColumnMappingGet")
    RETURN
999 NULLIFY(columnMapping)
998 ERRORSEXITS("DistributedMatrix_ColumnMappingGet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_ColumnMappingGet
  
  !
  !================================================================================================================================
  !

  !>Gets the data type of a distributed matrix.
  SUBROUTINE DistributedMatrix_DataTypeGet(distributedMatrix,dataType,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(OUT) :: dataType !<On return, the data type of the matrix. \see DistributedMatrixVector_DataTypes,DistributedMatrixVector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("DistributedMatrix_DataTypeGet",err,error,*999)

    CALL DistributedMatrix_AssertIsFinished(distributedMatrix,err,error,*999)
     
    dataType=distributedMatrix%dataType
 
    EXITS("DistributedMatrix_DataTypeGet")
    RETURN
999 ERRORSEXITS("DistributedMatrix_DataTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_DataTypeGet

  !
  !================================================================================================================================
  !

  !>Gets the dimensions of a matrix on this computation node.
  SUBROUTINE DistributedMatrix_DimensionsGet(distributedMatrix,m,n,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix to get dimensions for
    INTEGER(INTG), INTENT(OUT) :: m !<On return, the number of rows in the matrix for this domain
    INTEGER(INTG), INTENT(OUT) :: n !<On return, the number of columns in the matrix for this domain
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(DistributedMatrixPETScType), POINTER :: petscMatrix
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedMatrix_DimensionsGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
#endif    
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_NumberOfRowsGet(cmissMatrix%matrix,m,err,error,*999)
      CALL Matrix_NumberOfColumnsGet(cmissMatrix%matrix,n,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      NULLIFY(petscMatrix)
      CALL DistributedMatrix_PETScMatrixGet(distributedMatrix,petscMatrix,err,error,*999)
      m=petscMatrix%m
      n=petscMatrix%n
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("DistributedMatrix_DimensionsGet")
    RETURN
999 ERRORSEXITS("DistributedMatrix_DimensionsGet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_DimensionsGet

  !
  !================================================================================================================================
  !

  !>Gets the maximum number of columns in each row of a distributed matrix.
  SUBROUTINE DistributedMatrix_MaxColumnsPerRowGet(distributedMatrix,maxColumnsPerRow,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(OUT) :: maxColumnsPerRow !<On return, the maximum number of columns in each row
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(DistributedMatrixPETScType), POINTER :: petscMatrix
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedMatrix_MaxColumnsPerRowGet",err,error,*999)

    CALL DistributedMatrix_AssertIsFinished(distributedMatrix,err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_MaxColumnsPerRowGet(cmissMatrix%matrix,maxColumnsPerRow,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      NULLIFY(petscMatrix)
      CALL DistributedMatrix_PETScMatrixGet(distributedMatrix,petscMatrix,err,error,*999)
      maxColumnsPerRow=petscMatrix%maximumColumnIndicesPerRow
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
 
    EXITS("DistributedMatrix_MaxColumnsPerRowGet")
    RETURN
999 ERRORSEXITS("DistributedMatrix_MaxColumnsPerRowGet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_MaxColumnsPerRowGet

  !  
  !================================================================================================================================
  !
  
  !>Gets the number of columns for the distributed matrix
  SUBROUTINE DistributedMatrix_NumberOfColumnsGet(distributedMatrix,numberOfColumns,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix to get the number of columns for
    INTEGER(INTG), INTENT(OUT) :: numberOfColumns !<On return, the number of columns in the distributed matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(DistributedMatrixPETScType), POINTER :: petscMatrix
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedMatrix_NumberOfColumnsGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
#endif    

    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_NumberOfColumnsGet(cmissMatrix%matrix,numberOfColumns,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      NULLIFY(petscMatrix)
      CALL DistributedMatrix_PETScMatrixGet(distributedMatrix,petscMatrix,err,error,*999)
      numberOfColumns=petscMatrix%globalN
    CASE DEFAULT
      localError="The distributed matrix library type of "//TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
         
    EXITS("DistributedMatrix_NumberOfColumnsGet")
    RETURN
999 ERRORSEXITS("DistributedMatrix_NumberOfColumnsGet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_NumberOfColumnsGet
  
  !
  !================================================================================================================================
  !

  !>Gets the number of non zeros for a distributed matrix.
  SUBROUTINE DistributedMatrix_NumberOfNonZerosGet(distributedMatrix,numberOfNonZeros,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(OUT) :: numberOfNonZeros !<On return, the number of non zeros in the matrix to get
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(DistributedMatrixPETScType), POINTER :: petscMatrix
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedMatrix_NumberOfNonZerosGet",err,error,*999)

    CALL DistributedMatrix_AssertIsFinished(distributedMatrix,err,error,*999)
     
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_NumberOfNonZerosGet(cmissMatrix%matrix,numberOfNonZeros,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      NULLIFY(petscMatrix)
      CALL DistributedMatrix_PETScMatrixGet(distributedMatrix,petscMatrix,err,error,*999)
      numberOfNonZeros=petscMatrix%numberOfNonZeros
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("DistributedMatrix_NumberOfNonZerosGet")
    RETURN
999 ERRORSEXITS("DistributedMatrix_NumberOfNonZerosGet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_NumberOfNonZerosGet

  !  
  !================================================================================================================================
  !
  
  !>Gets the maximum number of rows for the distributed matrix
  SUBROUTINE DistributedMatrix_MaxNumberOfRowsGet(distributedMatrix,maxNumberOfRows,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix to get the maximum number of rows for
    INTEGER(INTG), INTENT(OUT) :: maxNumberOfRows !<On return, the maximum number of rows in the distributed matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedMatrix_MaxNumberOfRowsGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
#endif    

    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_MaxNumberOfRowsGet(cmissMatrix%matrix,maxNumberOfRows,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot get the maximum number of rows for a PETSc distributed matrix.",err,error,*999)
    CASE DEFAULT
      localError="The distributed matrix library type of "//TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
         
    EXITS("DistributedMatrix_MaxNumberOfRowsGet")
    RETURN
999 ERRORSEXITS("DistributedMatrix_MaxNumberOfRowsGet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_MaxNumberOfRowsGet
  
  !  
  !================================================================================================================================
  !
  
  !>Gets the number of blocks for the distributed matrix
  SUBROUTINE DistributedMatrix_NumberOfBlocksGet(distributedMatrix,numberOfRowBlocks,numberOfColumnBlocks,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix to get the number of blocks for
    INTEGER(INTG), INTENT(OUT) :: numberOfRowBlocks !<On return, the number of row blocks in the distributed matrix
    INTEGER(INTG), INTENT(OUT) :: numberOfColumnBlocks !<On return, the number of column blocks in the distributed matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedMatrix_NumberOfBlocksGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
#endif    

    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_NumberOfBlocksGet(cmissMatrix%matrix,numberOfRowBlocks,numberOfColumnBlocks,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot get the number of blocks for a PETSc matrix.",err,error,*999)
    CASE DEFAULT
      localError="The distributed matrix library type of "//TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
         
    EXITS("DistributedMatrix_NumberOfBlocksGet")
    RETURN
999 ERRORSEXITS("DistributedMatrix_NumberOfBlocksGet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_NumberOfBlocksGet
  
  !  
  !================================================================================================================================
  !
  
  !>Gets the number of global rows for the distributed matrix
  SUBROUTINE DistributedMatrix_NumberOfGlobalRowsGet(distributedMatrix,numberOfGlobalRows,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix to get the number of global rows for
    INTEGER(INTG), INTENT(OUT) :: numberOfGlobalRows !<On return, the number of global rows in the distributed matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(DistributedMatrixPETScType), POINTER :: petscMatrix
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedMatrix_NumberOfGlobalRowsGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
#endif    

    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      numberOfGlobalRows=cmissMatrix%globalM
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      NULLIFY(petscMatrix)
      CALL DistributedMatrix_PETScMatrixGet(distributedMatrix,petscMatrix,err,error,*999)
      numberOfGlobalRows=petscMatrix%globalM
    CASE DEFAULT
      localError="The distributed matrix library type of "//TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
         
    EXITS("DistributedMatrix_NumberOfGlobalRowsGet")
    RETURN
999 ERRORSEXITS("DistributedMatrix_NumberOfGlobalRowsGet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_NumberOfGlobalRowsGet
  
  !  
  !================================================================================================================================
  !
  
  !>Gets the number of local rows for the distributed matrix
  SUBROUTINE DistributedMatrix_NumberOfLocalRowsGet(distributedMatrix,numberOfLocalRows,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix to get the number of local rows for
    INTEGER(INTG), INTENT(OUT) :: numberOfLocalRows !<On return, the number of local rows in the distributed matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(DistributedMatrixPETScType), POINTER :: petscMatrix
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedMatrix_NumberOfLocalRowsGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
#endif    

    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_NumberOfRowsGet(cmissMatrix%matrix,numberOfLocalRows,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      NULLIFY(petscMatrix)
      CALL DistributedMatrix_PETScMatrixGet(distributedMatrix,petscMatrix,err,error,*999)
      numberOfLocalRows=petscMatrix%m
    CASE DEFAULT
      localError="The distributed matrix library type of "//TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
         
    EXITS("DistributedMatrix_NumberOfLocalRowsGet")
    RETURN
999 ERRORSEXITS("DistributedMatrix_NumberOfLocalRowsGet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_NumberOfLocalRowsGet
  
  !  
  !================================================================================================================================
  !
  
  !>Gets the PETSc matrix for the distributed matrix
  SUBROUTINE DistributedMatrix_PETScMatrixGet(distributedMatrix,petscMatrix,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix to get the PETSc matrix for
    TYPE(DistributedMatrixPETScType), POINTER :: petscMatrix !<On return, the PETSc matrix for the distributed matrix. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DistributedMatrix_PETScMatrixGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(petscMatrix)) CALL FlagError("PETSc matrix is associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
#endif    
    
    petscMatrix=>distributedMatrix%petsc

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(petscMatrix))  CALL FlagError("PETSc matrix is not associated for the distributed matrix.",err,error,*999)
#endif    
     
    EXITS("DistributedMatrix_PETScMatrixGet")
    RETURN
999 NULLIFY(petscMatrix)
998 ERRORSEXITS("DistributedMatrix_PETScMatrixGet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_PETScMatrixGet
  
  !
  !================================================================================================================================
  !
  
  !>Gets the row mapping for the distributed matrix
  SUBROUTINE DistributedMatrix_RowMappingGet(distributedMatrix,rowMapping,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix to get the row mapping for
    TYPE(DomainMappingType), POINTER :: rowMapping !<On return, the row mapping for the distributed matrix. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DistributedMatrix_RowMappingGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(rowMapping)) CALL FlagError("Row mapping is associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
#endif    
    
    rowMapping=>distributedMatrix%rowDomainMapping

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(rowMapping))  CALL FlagError("Row mapping is not associated for the distributed matrix.",err,error,*999)
#endif    
     
    EXITS("DistributedMatrix_RowMappingGet")
    RETURN
999 NULLIFY(rowMapping)
998 ERRORSEXITS("DistributedMatrix_RowMappingGet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_RowMappingGet
  
  !
  !================================================================================================================================
  !

  !>Gets the storage type of a distributed matrix.
  SUBROUTINE DistributedMatrix_StorageTypeGet(distributedMatrix,storageType,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(OUT) :: storageType !<On return, the storage (sparsity) type of the distributed matrix. \see MATRIX_VECTOR_StorageTypes,MATRIX_VECTOR
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(DistributedMatrixPETScType), POINTER :: petscMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedMatrix_StorageTypeGet",err,error,*999)

    CALL DistributedMatrix_AssertIsFinished(distributedMatrix,err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_StorageTypeGet(cmissMatrix%matrix,storageType,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      NULLIFY(petscMatrix)
      CALL DistributedMatrix_PETScMatrixGet(distributedMatrix,petscMatrix,err,error,*999)
      storageType=petscMatrix%storageType
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_StorageTypeGet")
    RETURN
999 ERRORSEXITS("DistributedMatrix_StorageTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_StorageTypeGet

  !
  !================================================================================================================================
  !

  !>Gets the symetry type of a distributed matrix.
  SUBROUTINE DistributedMatrix_SymmetryTypeGet(distributedMatrix,symmetryType,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix to get the symmetry type for
    INTEGER(INTG), INTENT(OUT) :: symmetryType !<On return, the symmetry type of the distributed matrix. \see DistributedMatrixVector_SymmetryTypes,DistributedMatrixVector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(DistributedMatrixPETScType), POINTER :: petscMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedMatrix_SymmetryTypeGet",err,error,*999)

    CALL DistributedMatrix_AssertIsFinished(distributedMatrix,err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_SymmetryTypeGet(cmissMatrix%matrix,symmetryType,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      NULLIFY(petscMatrix)
      CALL DistributedMatrix_PETScMatrixGet(distributedMatrix,petscMatrix,err,error,*999)
      symmetryType=petscMatrix%symmetryType
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_SymmetryTypeGet")
    RETURN
999 ERRORSEXITS("DistributedMatrix_SymmetryTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_SymmetryTypeGet

  !
  !================================================================================================================================
  !

  !>Gets the transpose type of a distributed matrix.
  SUBROUTINE DistributedMatrix_TransposeTypeGet(distributedMatrix,transposeType,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix to get the transpose type for
    INTEGER(INTG), INTENT(OUT) :: transposeType !<On return, the transpose type of the distributed matrix. \see DistributedMatrixVector_TransposeTypes,DistributedMatrixVector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(DistributedMatrixPETScType), POINTER :: petscMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedMatrix_TransposeTypeGet",err,error,*999)

    CALL DistributedMatrix_AssertIsFinished(distributedMatrix,err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_TransposeTypeGet(cmissMatrix%matrix,transposeType,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      NULLIFY(petscMatrix)
      CALL DistributedMatrix_PETScMatrixGet(distributedMatrix,petscMatrix,err,error,*999)
      transposeType=petscMatrix%transposeType
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_TransposeTypeGet")
    RETURN
999 ERRORSEXITS("DistributedMatrix_TransposeTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_TransposeTypeGet

  !  
  !================================================================================================================================
  !
  
  !>Gets the distributed matrix for the CMISS matrix
  SUBROUTINE DistributedMatrixCMISS_DistributedMatrixGet(cmissMatrix,distributedMatrix,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix !<A pointer to the CMISS matrix to get the distributed matrix for. 
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<On return, a pointer to the distributed matrix for the CMISS matrix for. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DistributedMatrixCMISS_DistributedMatrixGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is associated.",err,error,*998)    
    IF(.NOT.ASSOCIATED(cmissMatrix)) CALL FlagError("CMISS matrix is not associated.",err,error,*999)
#endif    
    
    distributedMatrix=>cmissMatrix%distributedMatrix

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(distributedMatrix)) &
      & CALL FlagError("Distributed matrix is not associated for the CMISS matrix.",err,error,*999)
#endif    
     
    EXITS("DistributedMatrixCMISS_DistributedMatrixGet")
    RETURN
999 NULLIFY(distributedMatrix)
998 ERRORSEXITS("DistributedMatrixCMISS_DistributedMatrixGet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrixCMISS_DistributedMatrixGet
  
  !  
  !================================================================================================================================
  !
  
  !>Gets the distributed matrix for the PETSc matrix
  SUBROUTINE DistributedMatrixPETSc_DistributedMatrixGet(petscMatrix,distributedMatrix,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixPETScType), POINTER :: petscMatrix !<A pointer to the PETSc matrix to get the distributed matrix for. 
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<On return, a pointer to the distributed matrix to get the PETSc matrix for. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DistributedMatrixPETSc_DistributedMatrixGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is associated.",err,error,*998)    
    IF(.NOT.ASSOCIATED(petscMatrix)) CALL FlagError("PETSc matrix is not associated.",err,error,*999)
#endif    
    
    distributedMatrix=>petscMatrix%distributedMatrix

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(distributedMatrix)) &
      & CALL FlagError("Distributed matrix is not associated for the PETSc matrix.",err,error,*999)
#endif
     
    EXITS("DistributedMatrixPETSc_DistributedMatrixGet")
    RETURN
999 NULLIFY(distributedMatrix)
998 ERRORSEXITS("DistributedMatrixPETSc_DistributedMatrixGet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrixPETSc_DistributedMatrixGet
  
  !
  !================================================================================================================================
  !

  !>Assert that a distributed vector has been finished
  SUBROUTINE DistributedVector_AssertIsFinished(distributedVector,err,error,*)

    !Argument Variables
    TYPE(DistributedVectorType), POINTER, INTENT(IN) :: distributedVector !<The distributed vector to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("DistributedVector_AssertIsFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
#endif    

    IF(.NOT.distributedVector%vectorFinished) CALL FlagError("Distributed vector has not been finished.",err,error,*999)
    
    EXITS("DistributedVector_AssertIsFinished")
    RETURN
999 ERRORSEXITS("DistributedVector_AssertIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_AssertIsFinished

  !
  !================================================================================================================================
  !

  !>Assert that a distributed vector has not been finished
  SUBROUTINE DistributedVector_AssertNotFinished(distributedVector,err,error,*)

    !Argument Variables
    TYPE(DistributedVectorType), POINTER, INTENT(IN) :: distributedVector !<The distributed vector to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("DistributedVector_AssertNotFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
#endif    

    IF(distributedVector%vectorFinished) CALL FlagError("Distributed vector has already been finished.",err,error,*999)
    
    EXITS("DistributedVector_AssertNotFinished")
    RETURN
999 ERRORSEXITS("DistributedVector_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_AssertNotFinished

  !
  !================================================================================================================================
  !

  !>Assert that a distributed vector has an integer data type
  SUBROUTINE DistributedVector_AssertIsIntgData(distributedVector,err,error,*)

    !Argument Variables
    TYPE(DistributedVectorType), POINTER, INTENT(IN) :: distributedVector !<The distributed vector to assert the integer data type for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("DistributedVector_AssertIsIntgData",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
#endif    

    IF(distributedVector%dataType/=DISTRIBUTED_MATRIX_VECTOR_INTG_TYPE) THEN
      localError="The distributed vector data type of "//TRIM(NumberToVString(distributedVector%dataType,"*",err,error))// &
        & " does not correspond to the required integer data type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("DistributedVector_AssertIsIntgData")
    RETURN
999 ERRORSEXITS("DistributedVector_AssertIsIntgData",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_AssertIsIntgData

  !
  !================================================================================================================================
  !

  !>Assert that a distributed vector has a single precision real data type
  SUBROUTINE DistributedVector_AssertIsSPData(distributedVector,err,error,*)

    !Argument Variables
    TYPE(DistributedVectorType), POINTER, INTENT(IN) :: distributedVector !<The distributed vector to assert the single precision real data type for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("DistributedVector_AssertIsSPData",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
#endif    

    IF(distributedVector%dataType/=DISTRIBUTED_MATRIX_VECTOR_SP_TYPE) THEN
      localError="The distributed vector data type of "//TRIM(NumberToVString(distributedVector%dataType,"*",err,error))// &
        & " does not correspond to the required single precision real data type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("DistributedVector_AssertIsSPData")
    RETURN
999 ERRORSEXITS("DistributedVector_AssertIsSPData",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_AssertIsSPData

  !
  !================================================================================================================================
  !

  !>Assert that a distributed vector has a double precision real data type
  SUBROUTINE DistributedVector_AssertIsDPData(distributedVector,err,error,*)

    !Argument Variables
    TYPE(DistributedVectorType), POINTER, INTENT(IN) :: distributedVector !<The distributed vector to assert the double precision real data type for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("DistributedVector_AssertIsDPData",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
#endif    

    IF(distributedVector%dataType/=DISTRIBUTED_MATRIX_VECTOR_DP_TYPE) THEN
      localError="The distributed vector data type of "//TRIM(NumberToVString(distributedVector%dataType,"*",err,error))// &
        & " does not correspond to the required double precision real data type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("DistributedVector_AssertIsDPData")
    RETURN
999 ERRORSEXITS("DistributedVector_AssertIsDPData",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_AssertIsDPData

  !
  !================================================================================================================================
  !

  !>Assert that a distributed vector has a logical data type
  SUBROUTINE DistributedVector_AssertIsLData(distributedVector,err,error,*)

    !Argument Variables
    TYPE(DistributedVectorType), POINTER, INTENT(IN) :: distributedVector !<The distributed vector to assert the logical data type for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("DistributedVector_AssertIsLData",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
#endif    

    IF(distributedVector%dataType/=DISTRIBUTED_MATRIX_VECTOR_L_TYPE) THEN
      localError="The distributed vector data type of "//TRIM(NumberToVString(distributedVector%dataType,"*",err,error))// &
        & " does not correspond to the required logical data type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("DistributedVector_AssertIsLData")
    RETURN
999 ERRORSEXITS("DistributedVector_AssertIsLData",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_AssertIsLData

  !  
  !================================================================================================================================
  !
  
  !>Gets the CMISS vector for the distributed vector
  SUBROUTINE DistributedVector_CMISSVectorGet(distributedVector,cmissVector,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector to get the CMISS vector for
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector !<On return, the CMISS vector for the distributed vector. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DistributedVector_CMISSVectorGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(cmissVector)) CALL FlagError("CMISS vector is associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
#endif    
    
    cmissVector=>distributedVector%cmiss

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(cmissVector))  CALL FlagError("CMISS vector is not associated for the distributed vector.",err,error,*999)
#endif    
     
    EXITS("DistributedVector_CMISSVectorGet")
    RETURN
999 NULLIFY(cmissVector)
998 ERRORSEXITS("DistributedVector_CMISSVectorGet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_CMISSVectorGet
  
  !
  !================================================================================================================================
  !

  !>Gets the data type of a distributed vector.
  SUBROUTINE DistributedVector_DataTypeGet(distributedVector,dataType,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    INTEGER(INTG), INTENT(OUT) :: dataType !<On return, the data type of the vector. \see DistributedMatrixVector_DataTypes,DistributedMatrixVector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("DistributedVector_DataTypeGet",err,error,*999)

    CALL DistributedVector_AssertIsFinished(distributedVector,err,error,*999)
    
    dataType=distributedVector%dataType

    EXITS("DistributedVector_DataTypeGet")
    RETURN
999 ERRORSEXITS("DistributedVector_DataTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_DataTypeGet

  !
  !================================================================================================================================
  !

  !>Gets the number of global rows in a distributed vector.
  SUBROUTINE DistributedVector_NumberOfGlobalRowsGet(distributedVector,numberOfGlobalRows,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    INTEGER(INTG), INTENT(OUT) :: numberOfGlobalRows !<On return, the number of global rows in the vector.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DomainMappingType), POINTER :: rowMapping

    ENTERS("DistributedVector_NumberOfGlobalRowsGet",err,error,*999)

    CALL DistributedVector_AssertIsFinished(distributedVector,err,error,*999)
#ifdef WITH_PRECHECKS    
    rowMapping=>distributedVector%domainMapping    
    IF(.NOT.ASSOCIATED(rowMapping)) CALL FlagError("Row mapping is not associated for the distributed vector.",err,error,*999)
#endif    
    
    numberOfGlobalRows=rowMapping%numberOfGlobal

    EXITS("DistributedVector_NumberOfGlobalRowsGet")
    RETURN
999 ERRORSEXITS("DistributedVector_NumberOfGlobalRowsGet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_NumberOfGlobalRowsGet

  !
  !================================================================================================================================
  !

  !>Gets the number of local rows in a distributed vector.
  SUBROUTINE DistributedVector_NumberOfLocalRowsGet(distributedVector,numberOfLocalRows,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    INTEGER(INTG), INTENT(OUT) :: numberOfLocalRows !<On return, the number of local rows in the vector.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DomainMappingType), POINTER :: rowMapping

    ENTERS("DistributedVector_NumberOfLocalRowsGet",err,error,*999)

    CALL DistributedVector_AssertIsFinished(distributedVector,err,error,*999)
#ifdef WITH_PRECHECKS    
    rowMapping=>distributedVector%domainMapping
    IF(.NOT.ASSOCIATED(rowMapping)) CALL FlagError("Row mapping is not associated for the distributed vector.",err,error,*999)
#endif    
    
    numberOfLocalRows=rowMapping%numberOfLocal

    EXITS("DistributedVector_NumberOfLocalRowsGet")
    RETURN
999 ERRORSEXITS("DistributedVector_NumberOfLocalRowsGet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_NumberOfLocalRowsGet

  !  
  !================================================================================================================================
  !
  
  !>Gets the PETSc vector for the distributed vector
  SUBROUTINE DistributedVector_PETScVectorGet(distributedVector,petscVector,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector to get the PETSc vector for
    TYPE(DistributedVectorPETScType), POINTER :: petscVector !<On return, the PETSc vector for the distributed vector. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DistributedVector_PETScVectorGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(petscVector)) CALL FlagError("PETSc vector is associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
#endif    
    
    petscVector=>distributedVector%petsc

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(petscVector))  CALL FlagError("PETSc vector is not associated for the distributed vector.",err,error,*999)
#endif    
     
    EXITS("DistributedVector_PETScVectorGet")
    RETURN
999 NULLIFY(petscVector)
998 ERRORSEXITS("DistributedVector_PETScVectorGet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_PETScVectorGet
  
  !
  !================================================================================================================================
  !
  
  !>Gets the row mapping for the distributed vector
  SUBROUTINE DistributedVector_RowMappingGet(distributedVector,rowMapping,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector to get the row mapping for
    TYPE(DomainMappingType), POINTER :: rowMapping !<On return, the row mapping for the distributed vector. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DistributedVector_RowMappingGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(rowMapping)) CALL FlagError("Row mapping is associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
#endif    
    
    rowMapping=>distributedVector%domainMapping

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(rowMapping)) CALL FlagError("Row mapping is not associated for the distributed vector.",err,error,*999)
#endif    
     
    EXITS("DistributedVector_RowMappingGet")
    RETURN
999 NULLIFY(rowMapping)
998 ERRORSEXITS("DistributedVector_RowMappingGet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_RowMappingGet
  
  !
  !================================================================================================================================
  !

  !>Gets the total number of local rows in a distributed vector.
  SUBROUTINE DistributedVector_TotalNumberOfLocalRowsGet(distributedVector,totalNumberOfLocalRows,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    INTEGER(INTG), INTENT(OUT) :: totalNumberOfLocalRows !<On return, the total number of local rows in the vector.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DomainMappingType), POINTER :: rowMapping

    ENTERS("DistributedVector_TotalNumberOfLocalRowsGet",err,error,*999)

    CALL DistributedVector_AssertIsFinished(distributedVector,err,error,*999)
#ifdef WITH_PRECHECKS    
    rowMapping=>distributedVector%domainMapping
    IF(.NOT.ASSOCIATED(rowMapping)) CALL FlagError("Row mapping is not associated for the distributed vector.",err,error,*999)
#endif    
    
    totalNumberOfLocalRows=rowMapping%totalNumberOfLocal

    EXITS("DistributedVector_TotalNumberOfLocalRowsGet")
    RETURN
999 ERRORSEXITS("DistributedVector_TotalNumberOfLocalRowsGet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_TotalNumberOfLocalRowsGet

  !  
  !================================================================================================================================
  !
  
  !>Gets the distributed vector for the CMISS vector
  SUBROUTINE DistributedVectorCMISS_DistributedVectorGet(cmissVector,distributedVector,err,error,*)

    !Argument variables
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector !<A pointer to the CMISS vector to get the distributed vector for.
    TYPE(DistributedVectorType), POINTER :: distributedVector !<On return, the distributed vector for the CMISS vector. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DistributedVectorCMISS_DistributedVectorGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(cmissVector)) CALL FlagError("CMISS vector is not associated.",err,error,*999)
#endif    
    
    distributedVector=>cmissVector%distributedVector

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(distributedVector)) &
      & CALL FlagError("Distributed vector is not associated for the CMISS vector.",err,error,*999)
#endif    
     
    EXITS("DistributedVectorCMISS_DistributedVectorGet")
    RETURN
999 NULLIFY(distributedVector)
998 ERRORSEXITS("DistributedVectorCMISS_DistributedVectorGet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVectorCMISS_DistributedVectorGet
  
  !  
  !================================================================================================================================
  !
  
  !>Gets the distributed vector for the PETSc vector
  SUBROUTINE DistributedVectorPETSc_DistributedVectorGet(petscVector,distributedVector,err,error,*)

    !Argument variables
    TYPE(DistributedVectorPETScType), POINTER :: petscVector !<A pointer to the PETSc vector to get the distributed vector for.
    TYPE(DistributedVectorType), POINTER :: distributedVector !<On return, the distributed vector for the PETSc vector. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DistributedVectorPETSc_DistributedVectorGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(petscVector)) CALL FlagError("PETSc vector is not associated.",err,error,*999)
#endif    
    
    distributedVector=>petscVector%distributedVector

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(distributedVector)) &
      & CALL FlagError("Distributed vector is not associated for the PETSc vector.",err,error,*999)
#endif    
     
    EXITS("DistributedVectorPETSc_DistributedVectorGet")
    RETURN
999 NULLIFY(distributedVector)
998 ERRORSEXITS("DistributedVectorPETSc_DistributedVectorGet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVectorPETSc_DistributedVectorGet
  
  !
  !================================================================================================================================
  !
           
END MODULE DistributedMatrixVectorAccessRoutines
