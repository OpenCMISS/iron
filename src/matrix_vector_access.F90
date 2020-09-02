!> \file
!> \author Chris Bradley
!> \brief This module contains all matrix vector access method routines.
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

!> This module contains all matrix vector access method routines.
MODULE MatrixVectorAccessRoutines
  
  USE BaseRoutines
  USE Constants
  USE Kinds
  USE ISO_VARYING_STRING
  USE Strings
  USE Types
  
#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters
  

  !> \addtogroup OpenCMISS_MatrixVectorConstants OpenCMISS::Iron::MatrixVector::Constants
  !> \brief Matrix vector constants.
  !>@{
  !> \addtogroup MatrixVector_DataTypes MatrixVector::Constants::DataTypes
  !> \brief Matrix vector data types
  !> \see MatrixVector
  !>@{
  INTEGER(INTG), PARAMETER :: MATRIX_VECTOR_INTG_TYPE=INTEGER_TYPE !<Integer matrix-vector data type \see MatrixVector_DataTypes,MatrixVector
  INTEGER(INTG), PARAMETER :: MATRIX_VECTOR_SP_TYPE=SINGLE_REAL_TYPE !<Single precision real matrix-vector data type \see MatrixVector_DataTypes,MatrixVector
  INTEGER(INTG), PARAMETER :: MATRIX_VECTOR_DP_TYPE=DOUBLE_REAL_TYPE !<Double precision real matrix-vector data type \see MatrixVector_DataTypes,MatrixVector
  INTEGER(INTG), PARAMETER :: MATRIX_VECTOR_L_TYPE=LOGICAL_TYPE !<Logical matrix-vector data type \see MatrixVector_DataTypes,MatrixVector
  INTEGER(INTG), PARAMETER :: MATRIX_VECTOR_SPC_TYPE=SINGLE_COMPLEX_TYPE !<Single precision complex matrix-vector data type \see MatrixVector_DataTypes,MatrixVector
  INTEGER(INTG), PARAMETER :: MATRIX_VECTOR_DPC_TYPE=DOUBLE_COMPLEX_TYPE !<Double precision complex matrix-vector data type \see MatrixVector_DataTypes,MatrixVector
   !>@}
  
  !> \addtogroup MatrixVector_StorageTypes MatrixVector::Constants::StorageTypes
  !> \brief Matrix-vector storage type parameters
  !> \see MatrixVector_MatrixStorageStructures,MatrixVector
  !>@{
  INTEGER(INTG), PARAMETER :: MATRIX_BLOCK_STORAGE_TYPE=0 !<Matrix block storage type \see MatrixVector_StorageTypes,MatrixVector
  INTEGER(INTG), PARAMETER :: MATRIX_DIAGONAL_STORAGE_TYPE=1 !<Matrix diagonal storage type \see MatrixVector_StorageTypes,MatrixVector
  INTEGER(INTG), PARAMETER :: MATRIX_COLUMN_MAJOR_STORAGE_TYPE=2 !<Matrix column major storage type \see MatrixVector_StorageTypes,MatrixVector
  INTEGER(INTG), PARAMETER :: MATRIX_ROW_MAJOR_STORAGE_TYPE=3 !<Matrix row major storage type \see MatrixVector_StorageTypes,MatrixVector
  INTEGER(INTG), PARAMETER :: MATRIX_COMPRESSED_ROW_STORAGE_TYPE=4 !<Matrix compressed row storage type \see MatrixVector_StorageTypes,MatrixVector
  INTEGER(INTG), PARAMETER :: MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE=5 !<Matrix compressed column storage type \see MatrixVector_StorageTypes,MatrixVector
  INTEGER(INTG), PARAMETER :: MATRIX_ROW_COLUMN_STORAGE_TYPE=6 !<Matrix row-column storage type \see MatrixVector_StorageTypes,MatrixVector
  INTEGER(INTG), PARAMETER :: MATRIX_BLOCK_COMPRESSED_ROW_STORAGE_TYPE=7 !<Matrix block compressed row storage type \see MatrixVector_StorageTypes,MatrixVector
  INTEGER(INTG), PARAMETER :: MATRIX_MODIFIED_KRM_STORAGE_TYPE=8 !<Matrix modified Knuth-Rheinboldt-Mesztenyi storage type \see MatrixVector_StorageTypes,MatrixVector
  !>@}
  
  !> \addtogroup MatrixVector_SymmetryTypes MatrixVector::Constants::SymmetryTypes
  !> \brief Matrix-vector storage type parameters
  !> \see MatrixVector_MatrixStorageStructures,MatrixVector
  !>@{
  INTEGER(INTG), PARAMETER :: MATRIX_SYMMETRIC_TYPE=0 !<Matrix is symmetric \see MatrixVector_SymmetryTypes,MatrixVector
  INTEGER(INTG), PARAMETER :: MATRIX_HERMITIAN_TYPE=1 !<Matrix is Hermitian (complex symmetric) \see MatrixVector_SymmetryTypes,MatrixVector
  INTEGER(INTG), PARAMETER :: MATRIX_SKEW_SYMMETRIC_TYPE=2 !<Matrix is skew-symmetric \see MatrixVector_SymmetryTypes,MatrixVector
  INTEGER(INTG), PARAMETER :: MATRIX_UNSYMMETRIC_TYPE=3 !<Matrix is unsymmetric \see MatrixVector_SymmetryTypes,MatrixVector
  !>@}

  
  !> \addtogroup MatrixVector_TransposeTypes MatrixVector::Constants::TransposeTypes
  !> \brief Matrix-vector transpose type parameters
  !> \see MatrixVector_MatrixStorageStructures,MatrixVector
  !>@{
  INTEGER(INTG), PARAMETER :: MATRIX_NO_TRANSPOSE_REQUIRED=0 !<Matrix will not have a transpose required \see MatrixVector_TransposeTypes,MatrixVector
  INTEGER(INTG), PARAMETER :: MATRIX_PARTIAL_TRANSPOSE_REQUIRED=1 !<Matrix will require a partial transpose for selected rows/columns \see MatrixVector_TransposeTypes,MatrixVector
  INTEGER(INTG), PARAMETER :: MATRIX_FULL_TRANSPOSE_REQUIRED=2 !<Matrix will require a full transpose  \see MatrixVector_TransposeTypes,MatrixVector
  !>@}
  !>@}
  
  !Module types

  !Module variables

  !Interfaces
  
  INTERFACE Matrix_DataGet
    MODULE PROCEDURE Matrix_DataGetIntg
    MODULE PROCEDURE Matrix_DataGetSP
    MODULE PROCEDURE Matrix_DataGetDP
    MODULE PROCEDURE Matrix_DataGetL
  END INTERFACE Matrix_DataGet

  INTERFACE Vector_DataGet
    MODULE PROCEDURE Vector_DataGetIntg
    MODULE PROCEDURE Vector_DataGetSP
    MODULE PROCEDURE Vector_DataGetDP
    MODULE PROCEDURE Vector_DataGetL
  END INTERFACE Vector_DataGet
  
  PUBLIC MATRIX_VECTOR_INTG_TYPE,MATRIX_VECTOR_SP_TYPE,MATRIX_VECTOR_DP_TYPE,MATRIX_VECTOR_L_TYPE, &
    MATRIX_VECTOR_SPC_TYPE,MATRIX_VECTOR_DPC_TYPE

  PUBLIC MATRIX_BLOCK_STORAGE_TYPE,MATRIX_DIAGONAL_STORAGE_TYPE,MATRIX_COLUMN_MAJOR_STORAGE_TYPE,MATRIX_ROW_MAJOR_STORAGE_TYPE, &
    & MATRIX_COMPRESSED_ROW_STORAGE_TYPE,MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE,MATRIX_ROW_COLUMN_STORAGE_TYPE, &
    & MATRIX_BLOCK_COMPRESSED_ROW_STORAGE_TYPE,MATRIX_MODIFIED_KRM_STORAGE_TYPE

  PUBLIC MATRIX_SYMMETRIC_TYPE,MATRIX_HERMITIAN_TYPE,MATRIX_SKEW_SYMMETRIC_TYPE,MATRIX_UNSYMMETRIC_TYPE

  PUBLIC MATRIX_NO_TRANSPOSE_REQUIRED,MATRIX_PARTIAL_TRANSPOSE_REQUIRED,MATRIX_FULL_TRANSPOSE_REQUIRED

  PUBLIC Matrix_AssertIsFinished,Matrix_AssertNotFinished

  PUBLIC Matrix_AssertIsIntgData,Matrix_AssertIsSPData,Matrix_AssertIsDPData,Matrix_AssertIsLData

  PUBLIC Matrix_BlockSizeGet

  PUBLIC Matrix_DataGet

  PUBLIC Matrix_DataTypeGet

  PUBLIC Matrix_MaxColumnsPerRowGet

  PUBLIC Matrix_MaxNumberOfRowsGet

  PUBLIC Matrix_NumberOfBlocksGet

  PUBLIC Matrix_NumberOfColumnsGet

  PUBLIC Matrix_NumberOfNonZerosGet

  PUBLIC Matrix_NumberOfRowsGet

  PUBLIC Matrix_StorageLocationsGet

  PUBLIC Matrix_StorageTransposeLocationsGet

  PUBLIC Matrix_StorageTypeGet

  PUBLIC Matrix_SymmetryTypeGet

  PUBLIC Matrix_TransposeTypeGet

  PUBLIC MatrixRowColCoupling_RowColCouplingInfoGet

  PUBLIC MatrixRowColCoupling_NumberOfRowColsGet
  
  PUBLIC Vector_AssertIsFinished,Vector_AssertNotFinished
 
  PUBLIC Vector_AssertIsIntgData,Vector_AssertIsSPData,Vector_AssertIsDPData,Vector_AssertIsLData
  
  PUBLIC Vector_DataGet

  PUBLIC Vector_DataTypeGet
 
CONTAINS

  !
  !================================================================================================================================
  !

  !>Assert that a matrix has been finished
  SUBROUTINE Matrix_AssertIsFinished(matrix,err,error,*)

    !Argument Variables
    TYPE(MatrixType), POINTER, INTENT(IN) :: matrix !<The matrix to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("Matrix_AssertIsFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(matrix)) CALL FlagError("Matrix is not associated.",err,error,*999)
#endif 

    IF(.NOT.matrix%matrixFinished) CALL FlagError("Matrix has not been finished.",err,error,*999)
    
    EXITS("Matrix_AssertIsFinished")
    RETURN
999 ERRORSEXITS("Matrix_AssertIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_AssertIsFinished

  !
  !=================================================================================================================================
  !

  !>Assert that a matrix has not been finished
  SUBROUTINE Matrix_AssertNotFinished(matrix,err,error,*)

    !Argument Variables
    TYPE(MatrixType), POINTER, INTENT(IN) :: matrix !<The matrix to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Matrix_AssertNotFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(matrix)) CALL FlagError("Matrix is not associated.",err,error,*999)
#endif    

    IF(matrix%matrixFinished) CALL FlagError("Matrix has already been finished.",err,error,*999)
    
    EXITS("Matrix_AssertNotFinished")
    RETURN
999 ERRORSEXITS("Matrix_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_AssertNotFinished

  !
  !================================================================================================================================
  !

  !>Assert that a matrix has an integer data type
  SUBROUTINE Matrix_AssertIsINTGData(matrix,err,error,*)

    !Argument Variables
    TYPE(MatrixType), POINTER, INTENT(INOUT) :: matrix !<The matrix to assert the integer data type for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Matrix_AssertIsINTGData",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(matrix)) CALL FlagError("Matrix is not associated.",err,error,*999)
#endif    

    IF(matrix%dataType/=MATRIX_VECTOR_INTG_TYPE) THEN
      localError="The matrix data type of "//TRIM(NumberToVString(matrix%dataType,"*",err,error))// &
        & " does not correspond to the required integer data type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Matrix_AssertIsINTGData")
    RETURN
999 ERRORSEXITS("Matrix_AssertIsINTGData",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_AssertIsINTGData

  !
  !=================================================================================================================================
  !

  !>Assert that a matrix has a single precision real data type
  SUBROUTINE Matrix_AssertIsSPData(matrix,err,error,*)

    !Argument Variables
    TYPE(MatrixType), POINTER, INTENT(INOUT) :: matrix !<The matrix to assert the single precision real data type for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Matrix_AssertIsSPData",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(matrix)) CALL FlagError("Matrix is not associated.",err,error,*999)
#endif    

    IF(matrix%dataType/=MATRIX_VECTOR_SP_TYPE) THEN
      localError="The matrix data type of "//TRIM(NumberToVString(matrix%dataType,"*",err,error))// &
        & " does not correspond to the required single precision real data type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Matrix_AssertIsSPData")
    RETURN
999 ERRORSEXITS("Matrix_AssertIsSPData",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_AssertIsSPData

  !
  !=================================================================================================================================
  !

  !>Assert that a matrix has a double precision real data type
  SUBROUTINE Matrix_AssertIsDPData(matrix,err,error,*)

    !Argument Variables
    TYPE(MatrixType), POINTER, INTENT(INOUT) :: matrix !<The matrix to assert the double precision real data type for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Matrix_AssertIsDPData",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(matrix)) CALL FlagError("Matrix is not associated.",err,error,*999)
#endif    

    IF(matrix%dataType/=MATRIX_VECTOR_DP_TYPE) THEN
      localError="The matrix data type of "//TRIM(NumberToVString(matrix%dataType,"*",err,error))// &
        & " does not correspond to the required double precision real data type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Matrix_AssertIsDPData")
    RETURN
999 ERRORSEXITS("Matrix_AssertIsDPData",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_AssertIsDPData

  !
  !=================================================================================================================================
  !

  !>Assert that a matrix has a logical data type
  SUBROUTINE Matrix_AssertIsLData(matrix,err,error,*)

    !Argument Variables
    TYPE(MatrixType), POINTER, INTENT(INOUT) :: matrix !<The matrix to assert the logical data type for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Matrix_AssertIsLData",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(matrix)) CALL FlagError("Matrix is not associated.",err,error,*999)
#endif    

    IF(matrix%dataType/=MATRIX_VECTOR_L_TYPE) THEN
      localError="The matrix data type of "//TRIM(NumberToVString(matrix%dataType,"*",err,error))// &
        & " does not correspond to the required logical data type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Matrix_AssertIsLData")
    RETURN
999 ERRORSEXITS("Matrix_AssertIsLData",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_AssertIsLData

  !
  !================================================================================================================================
  !

  !>Gets the block size of a matrix.
  SUBROUTINE Matrix_BlockSizeGet(matrix,blockSize,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix to get the block size for
    INTEGER(INTG), INTENT(OUT) :: blockSize !<On return, the block size of the matrix.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("Matrix_BlockSizeGet",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    
    blockSize=matrix%blockSize

    EXITS("Matrix_BlockSizeGet")
    RETURN
999 ERRORSEXITS("Matrix_BlockSizeGet",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_BlockSizeGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the data of an integer matrix. Note: the values can be used for read operations but a Matrix_ValuesSet call must be used to change any values. The pointer should not be deallocated and should be returned with a Matrix_DataRestore call.
  SUBROUTINE Matrix_DataGetIntg(matrix,matrixData,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), POINTER :: matrixData(:) !<On return a pointer to the matrix data. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Matrix_DataGetIntg",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(matrixData)) CALL FlagError("Matrix data is already associated.",err,error,*998)
#endif    
    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    CALL Matrix_AssertIsINTGData(matrix,err,error,*999)
    
    matrixData=>matrix%dataIntg
    
    EXITS("Matrix_DataGetIntg")
    RETURN
999 NULLIFY(matrixData)
998 ERRORSEXITS("Matrix_DataGetIntg",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_DataGetIntg

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the data of a single precision real matrix. Note: the values can be used for read operations but a Matrix_ValuesSet call must be used to change any values. The pointer should not be deallocated and should be returned with a Matrix_DataRestore call.
  SUBROUTINE Matrix_DataGetSP(matrix,matrixData,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    REAL(SP), POINTER :: matrixData(:) !<On return a pointer to the matrix data. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Matrix_DataGetSP",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(matrixData)) CALL FlagError("Matrix data is already associated.",err,error,*998)
#endif    
    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    CALL Matrix_AssertIsSPData(matrix,err,error,*999)
    
    matrixData=>matrix%dataSP
    
    EXITS("Matrix_DataGetSP")
    RETURN
999 NULLIFY(matrixData)
998 ERRORSEXITS("Matrix_DataGetSP",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_DataGetSP

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the data of a double precision real matrix. Note: the values can be used for read operations but a Matrix_ValuesSet call must be used to change any values. The pointer should not be deallocated and should be returned with a Matrix_DataRestore call.
  SUBROUTINE Matrix_DataGetDP(matrix,matrixData,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    REAL(DP), POINTER :: matrixData(:) !<On return a pointer to the matrix data. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Matrix_DataGetDP",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(matrixData)) CALL FlagError("Matrix data is already associated.",err,error,*998)
#endif    
    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    CALL Matrix_AssertIsDPData(matrix,err,error,*999)
    
    matrixData=>matrix%dataDP
    
    EXITS("Matrix_DataGetDP")
    RETURN
999 NULLIFY(matrixData)
998 ERRORSEXITS("Matrix_DataGetDP",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_DataGetDP

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the data of a logical matrix. Note: the values can be used for read operations but a Matrix_ValuesSet call must be used to change any values. The pointer should not be deallocated and should be returned with a Matrix_DataRestore call.
  SUBROUTINE Matrix_DataGetL(matrix,matrixData,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    LOGICAL, POINTER :: matrixData(:) !<On return a pointer to the matrix data. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Matrix_DataGetL",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(matrixData)) CALL FlagError("Matrix data is already associated.",err,error,*998)
#endif    
    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    CALL Matrix_AssertIsLData(matrix,err,error,*999)
    
    matrixData=>matrix%dataL
    
    EXITS("Matrix_DataGetL")
    RETURN
999 NULLIFY(matrixData)
998 ERRORSEXITS("Matrix_DataGetL",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_DataGetL

  !
  !================================================================================================================================
  !

  !>Gets the data type of a matrix.
  SUBROUTINE Matrix_DataTypeGet(matrix,dataType,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), INTENT(OUT) :: dataType !<On return, the data type of the matrix. \see MatrixVector_DataTypes,MatrixVector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("Matrix_DataTypeGet",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    
    dataType=matrix%dataType

    EXITS("Matrix_DataTypeGet")
    RETURN
999 ERRORSEXITS("Matrix_DataTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_DataTypeGet

  !
  !================================================================================================================================
  !

  !>Gets the maximum number of columns in each row of a matrix.
  SUBROUTINE Matrix_MaxColumnsPerRowGet(matrix,maxColumnsPerRow,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix to get the maximum number of columns per row for
    INTEGER(INTG), INTENT(OUT) :: maxColumnsPerRow !<On return, the maximum number of columns in each row
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("Matrix_MaxColumnsPerRowGet",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)

    maxColumnsPerRow=matrix%maximumColumnIndicesPerRow
    
    EXITS("Matrix_MaxColumnsPerRowGet")
    RETURN
999 ERRORSEXITS("Matrix_MaxColumnsPerRowGet",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_MaxColumnsPerRowGet

  !
  !================================================================================================================================
  !

  !>Gets the maximum number of rows for a matrix.
  SUBROUTINE Matrix_MaxNumberOfRowsGet(matrix,maxNumberOfRows,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix to get the maximum number of rows for
    INTEGER(INTG), INTENT(OUT) :: maxNumberOfRows !<On return, the maximum number of rows in the matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Matrix_MaxNumberOfRowsGet",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)

    maxNumberOfRows=matrix%maxM

    EXITS("Matrix_MaxNumberOfRowsGet")
    RETURN
999 ERRORSEXITS("Matrix_MaxNumberOfRowsGet",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_MaxNumberOfRowsGet

  !
  !================================================================================================================================
  !

  !>Gets the number of columns for a matrix.
  SUBROUTINE Matrix_NumberOfColumnsGet(matrix,numberOfColumns,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix to get the number of columns for
    INTEGER(INTG), INTENT(OUT) :: numberOfColumns !<On return, The number of columns in the matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Matrix_NumberOfColumnsGet",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)

    numberOfColumns=matrix%n

    EXITS("Matrix_NumberOfColumnsGet")
    RETURN
999 ERRORSEXITS("Matrix_NumberOfColumnsGet",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_NumberOfColumnsGet

  !
  !================================================================================================================================
  !

  !>Gets the number of blocks for a matrix.
  SUBROUTINE Matrix_NumberOfBlocksGet(matrix,numberOfRowBlocks,numberOfColumnBlocks,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix to get the number of blocks for
    INTEGER(INTG), INTENT(OUT) :: numberOfRowBlocks !<On return, The number of row blocks in the matrix
    INTEGER(INTG), INTENT(OUT) :: numberOfColumnBlocks !<On return, The number of column blocks in the matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Matrix_NumberOfBlocksGet",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)

    numberOfRowBlocks=matrix%numberOfRowBlocks
    numberOfColumnBlocks=matrix%numberOfColumnBlocks
    
    EXITS("Matrix_NumberOfBlocksGet")
    RETURN
999 ERRORSEXITS("Matrix_NumberOfBlocksGet",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_NumberOfBlocksGet

  !
  !================================================================================================================================
  !

  !>Gets the number of non zeros for a matrix.
  SUBROUTINE Matrix_NumberOfNonZerosGet(matrix,numberOfNonZeros,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix to get the number of non zeros for
    INTEGER(INTG), INTENT(OUT) :: numberOfNonZeros !<On return, The number of non zeros in the matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Matrix_NumberOfNonZerosGet",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)

    numberOfNonZeros=matrix%numberOfNonZeros

    EXITS("Matrix_NumberOfNonZerosGet")
    RETURN
999 ERRORSEXITS("Matrix_NumberOfNonZerosGet",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_NumberOfNonZerosGet

  !
  !================================================================================================================================
  !

  !>Gets the number of rows for a matrix.
  SUBROUTINE Matrix_NumberOfRowsGet(matrix,numberOfRows,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix to get the number of rows for
    INTEGER(INTG), INTENT(OUT) :: numberOfRows !<On return, The number of rows in the matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Matrix_NumberOfRowsGet",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)

    numberOfRows=matrix%m

    EXITS("Matrix_NumberOfRowsGet")
    RETURN
999 ERRORSEXITS("Matrix_NumberOfRowsGet",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_NumberOfRowsGet

  !
  !================================================================================================================================
  !

  !>Gets the storage locations (sparsity pattern) of a matrix.
  SUBROUTINE Matrix_StorageLocationsGet(matrix,rowIndices,columnIndices,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), POINTER :: rowIndices(:) !<rowIndices(i). On return, the row index values for the matrix. Must not be associated on entry.
    INTEGER(INTG), POINTER :: columnIndices(:) !<columnIndices(i). On return, the column index values for the matrix. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Matrix_StorageLocationsGet",err,error,*997)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(rowIndices)) CALL FlagError("Row indicies is already associated.",err,error,*997)
    IF(ASSOCIATED(columnIndices)) CALL FlagError("Column indicies is already associated.",err,error,*998)
#endif    
    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    
    SELECT CASE(matrix%storageType)
    CASE(MATRIX_BLOCK_STORAGE_TYPE)
      CALL FlagError("Can not get matrix locations for a block storage matrix.",err,error,*999)
    CASE(MATRIX_DIAGONAL_STORAGE_TYPE)
      CALL FlagError("Can not get matrix locations for a diagonal storage matrix.",err,error,*999)
    CASE(MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
      CALL FlagError("Can not get matrix locations for a column major storage matrix.",err,error,*999)
    CASE(MATRIX_ROW_MAJOR_STORAGE_TYPE)
      CALL FlagError("Can not get matrix locations for a row major storage matrix.",err,error,*999)
    CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)          
      rowIndices=>matrix%rowIndices
      columnIndices=>matrix%columnIndices
    CASE(MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
      rowIndices=>matrix%rowIndices
      columnIndices=>matrix%columnIndices          
    CASE(MATRIX_ROW_COLUMN_STORAGE_TYPE)
      rowIndices=>matrix%rowIndices
      columnIndices=>matrix%columnIndices
    CASE(MATRIX_BLOCK_COMPRESSED_ROW_STORAGE_TYPE)          
      rowIndices=>matrix%rowIndices
      columnIndices=>matrix%columnIndices
    CASE(MATRIX_MODIFIED_KRM_STORAGE_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The matrix storage type of "//TRIM(NumberToVString(matrix%storageType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Matrix_StorageLocationsGet")
    RETURN
999 NULLIFY(columnIndices)
998 NULLIFY(rowIndices)
997 ERRORSEXITS("Matrix_StorageLocationsGet",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_StorageLocationsGet

  !
  !================================================================================================================================
  !

  !>Gets the transpose storage locations (sparsity pattern) of a matrix.
  SUBROUTINE Matrix_StorageTransposeLocationsGet(matrix,rowIndicesT,columnIndicesT,dataSwivelT,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), POINTER :: rowIndicesT(:) !<rowIndicesT(i). On return, the row index values for the transpose matrix. Must not be associated on entry.
    INTEGER(INTG), POINTER :: columnIndicesT(:) !<columnIndicesT(i). On return, the column index values for the transpose matrix. Must not be associated on entry.
    INTEGER(INTG), POINTER :: dataSwivelT(:) !<dataSwivelT(i). On return, the data swivel accessing for the transpose matrix value from the original matrix data. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Matrix_StorageTransposeLocationsGet",err,error,*997)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(rowIndicesT)) CALL FlagError("Row indicies is already associated.",err,error,*997)
    IF(ASSOCIATED(columnIndicesT)) CALL FlagError("Column indicies is already associated.",err,error,*998)
    IF(ASSOCIATED(dataSwivelT)) CALL FlagError("Data swivel is already associated.",err,error,*998)
#endif    
    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    
    SELECT CASE(matrix%storageType)
    CASE(MATRIX_BLOCK_STORAGE_TYPE)
      CALL FlagError("Can not get matrix locations for a block storage matrix.",err,error,*999)
    CASE(MATRIX_DIAGONAL_STORAGE_TYPE)
      CALL FlagError("Can not get matrix locations for a diagonal storage matrix.",err,error,*999)
    CASE(MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
      CALL FlagError("Can not get matrix locations for a column major storage matrix.",err,error,*999)
    CASE(MATRIX_ROW_MAJOR_STORAGE_TYPE)
      CALL FlagError("Can not get matrix locations for a row major storage matrix.",err,error,*999)
    CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)          
      rowIndicesT=>matrix%rowIndicesT
      columnIndicesT=>matrix%columnIndicesT
      dataSwivelT=>matrix%transposeDataSwivel
    CASE(MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
      rowIndicesT=>matrix%rowIndicesT
      columnIndicesT=>matrix%columnIndicesT
      dataSwivelT=>matrix%transposeDataSwivel
    CASE(MATRIX_ROW_COLUMN_STORAGE_TYPE)
      rowIndicesT=>matrix%rowIndicesT
      columnIndicesT=>matrix%columnIndicesT
    CASE(MATRIX_BLOCK_COMPRESSED_ROW_STORAGE_TYPE)          
      rowIndicesT=>matrix%rowIndicesT
      columnIndicesT=>matrix%columnIndicesT
    CASE(MATRIX_MODIFIED_KRM_STORAGE_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The matrix storage type of "//TRIM(NumberToVString(matrix%storageType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Matrix_StorageTransposeLocationsGet")
    RETURN
999 NULLIFY(columnIndicesT)
998 NULLIFY(rowIndicesT)
997 ERRORSEXITS("Matrix_StorageTransposeLocationsGet",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_StorageTransposeLocationsGet

  !
  !================================================================================================================================
  !

  !>Gets the storage type for a matrix.
  SUBROUTINE Matrix_StorageTypeGet(matrix,storageType,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix to get the storage type for
    INTEGER(INTG), INTENT(OUT) :: storageType !<On return, the storage type of the matrix. \see MatrixVector_StorageTypes,MatrixVector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Matrix_StorageTypeGet",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)

    storageType=matrix%storageType

    EXITS("Matrix_StorageTypeGet")
    RETURN
999 ERRORSEXITS("Matrix_StorageTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_StorageTypeGet

  !
  !================================================================================================================================
  !

  !>Gets the symmetry type for a matrix.
  SUBROUTINE Matrix_SymmetryTypeGet(matrix,symmetryType,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix to get the symmetry type for
    INTEGER(INTG), INTENT(OUT) :: symmetryType !<On return, the symmetry type of the matrix. \see MatrixVector_SymmetryTypes,MatrixVector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Matrix_SymmetryTypeGet",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)

    symmetryType=matrix%symmetryType

    EXITS("Matrix_SymmetryTypeGet")
    RETURN
999 ERRORSEXITS("Matrix_SymmetryTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_SymmetryTypeGet

  !
  !================================================================================================================================
  !

  !>Gets the transpose type for a matrix.
  SUBROUTINE Matrix_TransposeTypeGet(matrix,transposeType,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix to get the transpose type for
    INTEGER(INTG), INTENT(OUT) :: transposeType !<On return, the transpose type of the matrix. \see MatrixVector_TransposeTypes,MatrixVector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Matrix_TransposeTypeGet",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)

    transposeType=matrix%transposeType

    EXITS("Matrix_TransposeTypeGet")
    RETURN
999 ERRORSEXITS("Matrix_TransposeTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_TransposeTypeGet

  !
  !================================================================================================================================
  !

  !>Gets the row/col coupling information in the coupling for a row/column in a matrix row/column coupling.
  SUBROUTINE MatrixRowColCoupling_RowColCouplingInfoGet(matrixRowColCoupling,rowColumnNumber,coupledRowColIdx, &
    & coupledRowColNumber,couplingCoefficient,err,error,*)

    !Argument variables
    TYPE(MatrixRowColCouplingType), POINTER :: matrixRowColCoupling(:) !<A pointer to the matrix row/column coupling
    INTEGER(INTG), INTENT(IN) :: rowColumnNumber !<The row/column number to get the coupling for
    INTEGER(INTG), INTENT(IN) :: coupledRowColIdx !<The index of the coupled row/col to get the coupling information for
    INTEGER(INTG), INTENT(OUT) :: coupledRowColNumber !<On return, the row/col number the specified row/column is coupled to
    REAL(DP), INTENT(OUT) :: couplingCoefficient !<On return, the row/col coupling coefficient
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("MatrixRowColCoupling_RowColCouplingInfoGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(matrixRowColCoupling)) CALL FlagError("The matrix row column coupling is not associated.",err,error,*999)
    IF(rowColumnNumber<1.OR.rowColumnNumber>SIZE(matrixRowColCoupling,1)) THEN
      localError="The specified row/column number of "//TRIM(NumberToVString(rowColumnNumber,"*",err,error))// &
        & " is invalid. The row/column number should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(matrixRowColCoupling,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(coupledRowColIdx<1.OR.coupledRowColIdx>matrixRowColCoupling(rowColumnNumber)%numberOfRowCols) THEN
      localError="The specified row/column index of "//TRIM(NumberToVString(coupledRowColIdx,"*",err,error))// &
        & " is invalid. The row/column index should be >= 1 and <= "// &
        & TRIM(NumberToVString(matrixRowColCoupling(rowColumnNumber)%numberOfRowCols,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(matrixRowColCoupling(rowColumnNumber)%rowCols)) THEN
      localError="The row columns array is not allocated for row/column number "// &
        & TRIM(NumberToVString(rowColumnNumber,"*",err,error))//" of the matrix row column coupling."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(matrixRowColCoupling(rowColumnNumber)%couplingCoefficients)) THEN
      localError="The coupling coefficients array is not allocated for row/column number "// &
        & TRIM(NumberToVString(rowColumnNumber,"*",err,error))//" of the matrix row column coupling."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    coupledRowColNumber=matrixRowColCoupling(rowColumnNumber)%rowCols(coupledRowColIdx)
    couplingCoefficient=matrixRowColCoupling(rowColumnNumber)%couplingCoefficients(coupledRowColIdx)

    EXITS("MatrixRowColCoupling_RowColCouplingInfoGet")
    RETURN
999 ERRORSEXITS("MatrixRowColCoupling_RowColCouplingInfoGet",err,error)
    RETURN 1
    
  END SUBROUTINE MatrixRowColCoupling_RowColCouplingInfoGet

  !
  !================================================================================================================================
  !

  !>Gets the number of row/cols in the coupling for a row/column in a matrix row/column coupling.
  SUBROUTINE MatrixRowColCoupling_NumberOfRowColsGet(matrixRowColCoupling,rowColumnNumber,numberOfRowCols,err,error,*)

    !Argument variables
    TYPE(MatrixRowColCouplingType), POINTER :: matrixRowColCoupling(:) !<A pointer to the matrix row/column coupling
    INTEGER(INTG), INTENT(IN) :: rowColumnNumber !<The row/column number to get the coupling for
    INTEGER(INTG), INTENT(OUT) :: numberOfRowCols !<On return, the number of row/cols the specified row/column is coupled to
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("MatrixRowColCoupling_NumberOfRowColsGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(matrixRowColCoupling)) CALL FlagError("The matrix row column coupling is not associated.",err,error,*999)
    IF(rowColumnNumber<1.OR.rowColumnNumber>SIZE(matrixRowColCoupling,1)) THEN
      localError="The specified row/column number of "//TRIM(NumberToVString(rowColumnNumber,"*",err,error))// &
        & " is invalid. The row/column number should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(matrixRowColCoupling,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    numberOfRowCols=matrixRowColCoupling(rowColumnNumber)%numberOfRowCols

    EXITS("MatrixRowColCoupling_NumberOfRowColsGet")
    RETURN
999 ERRORSEXITS("MatrixRowColCoupling_NumberOfRowColsGet",err,error)
    RETURN 1
    
  END SUBROUTINE MatrixRowColCoupling_NumberOfRowColsGet

  !
  !=================================================================================================================================
  !

  !>Assert that a vector has been finished
  SUBROUTINE Vector_AssertIsFinished(vector,err,error,*)

    !Argument Variables
    TYPE(VectorType), POINTER, INTENT(IN) :: vector !<The vector to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Vector_AssertIsFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(vector)) CALL FlagError("Vector is not associated.",err,error,*999)
#endif    

    IF(.NOT.vector%vectorFinished) CALL FlagError("Vector has not been finished.",err,error,*999)
    
    EXITS("Vector_AssertIsFinished")
    RETURN
999 ERRORSEXITS("Vector_AssertIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE Vector_AssertIsFinished

  !
  !=================================================================================================================================
  !
  !>Assert that a vector has not been finished
  SUBROUTINE Vector_AssertNotFinished(vector,err,error,*)

    !Argument Variables
    TYPE(VectorType), POINTER, INTENT(IN) :: vector !<The vector to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Vector_AssertNotFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(vector)) CALL FlagError("Vector is not associated.",err,error,*999)
#endif    

    IF(vector%vectorFinished) CALL FlagError("Vector has already been finished.",err,error,*999)
    
    EXITS("Vector_AssertNotFinished")
    RETURN
999 ERRORSEXITS("Vector_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE Vector_AssertNotFinished

  !
  !=================================================================================================================================
  !

  !>Assert that a vector has an integer data type
  SUBROUTINE Vector_AssertIsINTGData(vector,err,error,*)

    !Argument Variables
    TYPE(VectorType), POINTER, INTENT(INOUT) :: vector !<The vector to assert the integer data type for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Vector_AssertIsINTGData",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(vector)) CALL FlagError("Vector is not associated.",err,error,*999)
#endif    

    IF(vector%dataType/=MATRIX_VECTOR_INTG_TYPE) THEN
      localError="The vector data type of "//TRIM(NumberToVString(vector%dataType,"*",err,error))// &
        & " does not correspond to the required integer data type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Vector_AssertIsINTGData")
    RETURN
999 ERRORSEXITS("Vector_AssertIsINTGData",err,error)
    RETURN 1
    
  END SUBROUTINE Vector_AssertIsINTGData

  !
  !=================================================================================================================================
  !

  !>Assert that a vector has a single precision real data type
  SUBROUTINE Vector_AssertIsSPData(vector,err,error,*)

    !Argument Variables
    TYPE(VectorType), POINTER, INTENT(INOUT) :: vector !<The vector to assert the single precision real data type for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Vector_AssertIsSPData",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(vector)) CALL FlagError("Vector is not associated.",err,error,*999)
#endif    

    IF(vector%dataType/=MATRIX_VECTOR_SP_TYPE) THEN
      localError="The vector data type of "//TRIM(NumberToVString(vector%dataType,"*",err,error))// &
        & " does not correspond to the required single precision real data type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Vector_AssertIsSPData")
    RETURN
999 ERRORSEXITS("Vector_AssertIsSPData",err,error)
    RETURN 1
    
  END SUBROUTINE Vector_AssertIsSPData

  !
  !=================================================================================================================================
  !

  !>Assert that a vector has a double precision real data type
  SUBROUTINE Vector_AssertIsDPData(vector,err,error,*)

    !Argument Variables
    TYPE(VectorType), POINTER, INTENT(INOUT) :: vector !<The vector to assert the double precision real data type for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Vector_AssertIsDPData",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(vector)) CALL FlagError("Vector is not associated.",err,error,*999)
#endif    

    IF(vector%dataType/=MATRIX_VECTOR_DP_TYPE) THEN
      localError="The vector data type of "//TRIM(NumberToVString(vector%dataType,"*",err,error))// &
        & " does not correspond to the required double precision real data type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Vector_AssertIsDPData")
    RETURN
999 ERRORSEXITS("Vector_AssertIsDPData",err,error)
    RETURN 1
    
  END SUBROUTINE Vector_AssertIsDPData

  !
  !=================================================================================================================================
  !

  !>Assert that a vector has a logical data type
  SUBROUTINE Vector_AssertIsLData(vector,err,error,*)

    !Argument Variables
    TYPE(VectorType), POINTER, INTENT(INOUT) :: vector !<The vector to assert the logical data type for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Vector_AssertIsLData",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(vector)) CALL FlagError("Vector is not associated.",err,error,*999)
#endif    

    IF(vector%dataType/=MATRIX_VECTOR_L_TYPE) THEN
      localError="The vector data type of "//TRIM(NumberToVString(vector%dataType,"*",err,error))// &
        & " does not correspond to the required logical data type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Vector_AssertIsLData")
    RETURN
999 ERRORSEXITS("Vector_AssertIsLData",err,error)
    RETURN 1
    
  END SUBROUTINE Vector_AssertIsLData

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the data of an integer vector. Note: the values can be used for read operations but a Vector_ValuesSet call must be used to change any values. The pointer should not be deallocated and a Vector_DataRestore call should be used to return the data.
  SUBROUTINE Vector_DataGetIntg(vector,vectorData,err,error,*)

    !Argument variables
    TYPE(VectorType), POINTER :: vector !<A pointer to the vector to get the vector data for
    INTEGER(INTG), POINTER :: vectorData(:) !<On return, a pointer to the vector data. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Vector_DataGetIntg",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(vectorData)) CALL FlagError("Vector data is already associated.",err,error,*998)
#endif    
    CALL Vector_AssertIsFinished(vector,err,error,*999)
    CALL Vector_AssertIsINTGData(vector,err,error,*999)
    
    vectorData=>vector%dataIntg
       
    EXITS("Vector_DataGetIntg")
    RETURN
999 NULLIFY(vectorData)
998 ERRORSEXITS("Vector_DataGetIntg",err,error)
    RETURN 1
    
  END SUBROUTINE Vector_DataGetIntg

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the data of a single precision real vector. Note: the values can be used for read operations but a Vector_ValuesSet call must be used to change any values. The pointer should not be deallocated and a Vector_DataRestore call should be used to return the data.
  SUBROUTINE Vector_DataGetSP(vector,vectorData,err,error,*)

    !Argument variables
    TYPE(VectorType), POINTER :: vector !<A pointer to the vector to get the vector data for
    REAL(SP), POINTER :: vectorData(:) !<On return, a pointer to the vector data. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Vector_DataGetSP",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(vectorData)) CALL FlagError("Vector data is already associated.",err,error,*998)
#endif    
    CALL Vector_AssertIsFinished(vector,err,error,*999)
    CALL Vector_AssertIsSPData(vector,err,error,*999)
    
    vectorData=>vector%dataSP
       
    EXITS("Vector_DataGetSP")
    RETURN
999 NULLIFY(vectorData)
998 ERRORSEXITS("Vector_DataGetSP",err,error)
    RETURN 1
    
  END SUBROUTINE Vector_DataGetSP

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the data of a double precision real vector. Note: the values can be used for read operations but a Vector_ValuesSet call must be used to change any values. The pointer should not be deallocated and a Vector_DataRestore call should be used to return the data.
  SUBROUTINE Vector_DataGetDP(vector,vectorData,err,error,*)

    !Argument variables
    TYPE(VectorType), POINTER :: vector !<A pointer to the vector to get the vector data for
    REAL(DP), POINTER :: vectorData(:) !<On return, a pointer to the vector data. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Vector_DataGetDP",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(vectorData)) CALL FlagError("Vector data is already associated.",err,error,*998)
#endif    
    CALL Vector_AssertIsFinished(vector,err,error,*999)
    CALL Vector_AssertIsDPData(vector,err,error,*999)
    
    vectorData=>vector%dataDP
       
    EXITS("Vector_DataGetDP")
    RETURN
999 NULLIFY(vectorData)
998 ERRORSEXITS("Vector_DataGetDP",err,error)
    RETURN 1
    
  END SUBROUTINE Vector_DataGetDP

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the data of a logical vector. Note: the values can be used for read operations but a Vector_ValuesSet call must be used to change any values. The pointer should not be deallocated and a Vector_DataRestore call should be used to return the data.
  SUBROUTINE Vector_DataGetL(vector,vectorData,err,error,*)

    !Argument variables
    TYPE(VectorType), POINTER :: vector !<A pointer to the vector to get the vector data for
    LOGICAL, POINTER :: vectorData(:) !<On return, a pointer to the vector data. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Vector_DataGetL",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(vectorData)) CALL FlagError("Vector data is already associated.",err,error,*998)
#endif    
    CALL Vector_AssertIsFinished(vector,err,error,*999)
    CALL Vector_AssertIsLData(vector,err,error,*999)
    
    vectorData=>vector%dataL
       
    EXITS("Vector_DataGetL")
    RETURN
999 NULLIFY(vectorData)
998 ERRORSEXITS("Vector_DataGetL",err,error)
    RETURN 1
    
  END SUBROUTINE Vector_DataGetL

  !
  !================================================================================================================================
  !

  !>Gets the data type of a vector.
  SUBROUTINE Vector_DataTypeGet(vector,dataType,err,error,*)

    !Argument variables
    TYPE(VectorType), POINTER :: vector !<A pointer to the vector
    INTEGER(INTG), INTENT(OUT) :: dataType !<On return, the data type of the vector. \see MatrixVector_DataTypes,MatrixVector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("Vector_DataTypeGet",err,error,*999)

    CALL Vector_AssertIsFinished(vector,err,error,*999)
    
    dataType=vector%dataType

    EXITS("Vector_DataTypeGet")
    RETURN
999 ERRORSEXITS("Vector_DataTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE Vector_DataTypeGet

  !  
  !================================================================================================================================
  !
  
END MODULE MatrixVectorAccessRoutines
