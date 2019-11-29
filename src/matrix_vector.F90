!> \file
!> \author Chris Bradley
!> \brief This module contains all routines dealing with (non-distributed) matrix and vectors types.
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
!> \section MatrixVector_MatrixStorageStructures MATRIX STORAGE STRUCTURES
!>    The matrix storage structures used are governed by the STORAGE parameter associated with the array. If STORAGE is
!>    MATRIX_BLOCK_STORAGE_TYPE the matrix is not sparse and the the non-sparse matrix dimension M is used to calculate
!>    the matrix storage locations. If storage is MATRIX_DIAGONAL_STORAGE_TYPE then only the matrix diagonal is stored.
!>    If storage is MATRIX_COLUMN_MAJOR_STORAGE_TYPE the matrix is not sparse and the non-sparse
!>    matrix dimension MAX_M (>=M) is used to calcualte the matrix storage locations. If storage is
!>    MATRIX_ROW_MAJOR_STORAGE_TYPE the matrix is not sparse and the non-sparse matrix dimension MAX_N (>=N) is used to
!>    calcualte the matrix storage locations. If STORAGE is MATRIX_COMPRESSED_ROW_STORAGE_TYPE the matrix has compressed row
!>    storage/sparsity (see below) and the sparsity structure arrays rowIndices and columnIndices are used for the
!>    storage location calculation. If STORAGE is MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE the matrix has compressed column
!>    storage/sparsity (see below) and the sparsity structure arrays rowIndices and columnIndices are used for the
!>    storage location calculation. If STORAGE is MATRIX_ROW_COLUMN_STORAGE_TYPE the matrix has row column
!>    storage/sparsity (see below) and the sparsity structure arrays rowIndices and columnIndices are used for the
!>    storage location calculation. 
!>    
!>    \subsection MatrixVector_CompressedRowStorage COMPRESSED-ROW STORAGE:
!>    
!>      The storage structure scheme is based on storing a MxN matrix as a one dimensional array of length SIZE
!>      (=numberOfNonZeros) (where numberOfNonZeros=sxMxN, s is the sparsity of the array) that stores only the non-zero
!>      elements of the matrix. Two additional arrays rowIndices and columnIndices store the positions of the non-zero
!>      elements. rowIndices is of length M+1 and COLUMN is of length numberOfNonZeros. rowIndices(i) stores the position
!>      in columnIndices of  the start of row i. The M+1 position of rowIndices stores the size of columnIndices+1 i.e.,
!>      numberOfNonZeros+1. The number of non-zero elements in row i can be found from rowIndices(i+1)-rowIndices(i).
!>      columnIndices(nz) gives the column number for non-zero element nz. See also COMPRESSED-COLUMN storage.
!>   
!>      Example of the compressed-row storage scheme on a NxN matrix (N=6). Here the sparsity is 8/36 or 22%
!>      \verbatim
!>      
!>      GX  1 2 3 4 5 6        
!>         ____________        DATA(nz)                        
!>       1| 0 A 0 B 0 0          A B C D E F G H  
!>       2| 0 0 C 0 0 0          
!>       3| 0 0 0 0 D E        rowIndices(i)
!>       4| F 0 0 0 0 0          1 3 4 6 7 8 9
!>       5| 0 0 G 0 0 0        columnIndices(i)  
!>       6| 0 0 0 0 0 H          2 4 3 5 6 1 3 6
!>      
!>      \endverbatim
!>
!>    \subsection MatrixVector_CompressedColumnStorage COMPRESSED-COLUMN STORAGE:
!>    
!>      The storage structure scheme is based on storing a MxN matrix as a one dimensional array of length SIZE
!>      (=numberOfNonZeros) (where numberOfNonZeros=sxMxN, s is the sparsity of the array) that stores only the non-zero
!>      elements of the matrix. Two additional arrays rowIndices and columnIndices store the positions of the non-zero
!>      elements. rowIndices is of length numberOfNonZeros and COLUMN is of length N+1. columnIndices(j) stores the position
!>      in rowIndices of  the start of column j. The N+1 position of columnIndices stores the size of rowIndices+1 i.e.,
!>      numberOfNonZeros+1. The number of non-zero elements in column j can be found from columnIndices(j+1)-columnIndices(j).
!>      rowIndices(nz) gives the row number for non-zero element nz. See also COMPRESSED-ROW storage.
!>   
!>      Example of compressed-column storage scheme on a NxN matrix (N=6). Here the sparsity is 8/36 or 22%
!>      \verbatim
!>      
!>      GX  1 2 3 4 5 6        
!>         ____________        DATA(nz)                        
!>       1| 0 A 0 B 0 0          F A C G B D E H  
!>       2| 0 0 C 0 0 0          
!>       3| 0 0 0 0 D E        rowIndices(i)
!>       4| F 0 0 0 0 0          4 1 2 5 1 3 3 6
!>       5| 0 0 G 0 0 0        columnIndices(i)  
!>       6| 0 0 0 0 0 H          1 2 3 5 6 7 9
!>      
!>      \endverbatim
!>
!>    \subsection MatrixVector_RowColumnStorage ROW-COLUMN STORAGE:
!>    
!>      The storage structure scheme is based on storing a MxN matrix as a one dimensional array of length SIZE
!>      (=numberOfNonZeros) (where numberOfNonZeros=sxMxN, s is the sparsity of the array) that stores only the non-zero
!>      elements of the matrix. Two additional arrays rowIndices and columnIndices store the positions of the non-zero
!>      elements. Both rowIndices and columnIndices are of length numberOfNonZeros. rowIndices(nz) gives the row number for
!>      non-zero element nz and columnIndices(nz) gives the column number for non-zero element nz.
!>  
!>      Example of row-column storage scheme on a NxN matrix (N=6). Here the sparsity is 8/36 or 22%
!>      \verbatim
!>   
!>      GX  1 2 3 4 5 6        
!>         ____________        DATA(nz)                        
!>       1| 0 A 0 B 0 0          A B C D E F G H  
!>       2| 0 0 C 0 0 0          
!>       3| 0 0 0 0 D E        rowIndices(i)
!>       4| F 0 0 0 0 0          1 1 2 3 3 4 5 6 
!>       5| 0 0 G 0 0 0        columnIndices(i)  
!>       6| 0 0 0 0 0 H          2 4 3 5 6 1 3 6
!>
!>      \endverbatim
!>
!>    \subsection MatrixVector_BlockCompressedRowStorage BLOCK COMPRESSED-ROW STORAGE:
!>    
!>      The storage structure scheme is based on storing a MxN matrix as a number of blocks of size blockSize. Blocks are
!>      only stored if they have one non-zero element. All entries in the block (zero's included) are stored in the data
!>      array in row-major order. The row and column indices are then used to locate the non-zero blocks inside the matrix.
!>      rowIndicies is of size numberOfRowBlocks + 1 and columnIndicies is of size numberOfBlocks. The row and column indices
!>      use the compressed row storage scheme to locate the non-zero blocks. 
!>   
!>      Example of the compressed-row storage scheme on a NxN matrix (N=6). Here the sparsity is 8/36 or 22%. The matrix is
!>      split into 3x3 blocks of 2x2 entries.
!>      \verbatim
!>                             BLOCKS
!>      GX  1 2 3 4 5 6                                       data(nz) 0 0 A 0 0 C B 0 0 F 0 0 D 0 E 0 G 0 0 0 0 0 0 H
!>         ____________        0 A   0 B   0 0                       
!>       1| 0 A 0 B 0 0        0 0   C 0   0 0                rowIndices(i)
!>       2| 0 0 C 0 0 0                                        1 3 5 7     
!>       3| 0 0 0 0 D E        0 0   0 0   D E                columnIndices(i)
!>       4| F 0 0 0 0 0        F 0   0 0   0 0                 1 2 1 3 2 3
!>       5| 0 0 G 0 0 0                                       
!>       6| 0 0 0 0 0 H        0 0   G 0   0 0                 
!>                             0 0   0 0   0 H                               
!>                                                            
!>      \endverbatim
!>
!>    \subsection MatrixVector_ModifiedKRMStorage MODIFIED KRM STORAGE:
!>    
!>      The storage structure scheme is based on storing a MxN matrix in a modified Knuth-Rheinboldt-Mesztenyi format. This
!>      is also known as Duff-Reid or Harwell. The format allows for the matrix to be indexed by row or by column. The non-zero
!>      entries are stored in a DATA array of length number of zeros. The row and column indicies arrays are also of number of
!>      non-zero size. Two additional arrays give the start for each row (size M) and the start for each column (size N).
!>      The row and column indicies index using a circular linked-list. The end of each row or column is indicated by a negative
!>      row or column number. For example consider the matrix below. To loop across row 3 you go to rowStart(3). This gives the
!>      non-zero index of 4. The first element in the row is thus data(4) or D. The next element in the row is given by
!>      rowIndices(4) which is 5. As the entry is non-negative the next entry in the row is given by data(5) which is E. The
!>      next entry is given by rowIndices(5) which is -3. The value is the negative of the row number and so this indicates the
!>      last entry in the row. To loop down column 3 we go to columnStart(3) which is 3. The first entry in the column is given
!>      by data(3) which is C. The next entry in the column is given by columnIndices(3) which is 7. The next value in the column
!>      is given by data(7) which is G. The next entry in the column is given by columnIndices(7) which is -3. The value is the
!>      negative of the column number and so this indicates the last entry in the column. To find out the row and column number
!>      for a particular non-zero number you need to scan te row and column indicies array for a negative number. For example
!>      consider the non-zero number 4 which has a value of data(4) or D. To find the row number you look at rowIndices(4) which
!>      has value 5. As this is not negative you move along to rowIndicies(5) which is a value of -3. This is negative and so
!>      D has a row number of 3. To find the column number you look at columnIndices(4) which has a value of -5. This is negative
!>      and so D has a column number of 5.
!>  
!>      \verbatim
!>                             data(nz)         
!>      GX  1 2 3 4 5 6          A  B  C  D  E  F  G  H          
!>         ____________     
!>       1| 0 A 0 B 0 0        rowIndices(nz)        
!>       2| 0 0 C 0 0 0          2 -1 -2  5 -3 -4 -5 -6         
!>       3| 0 0 0 0 D E        columnIndices(nz)        
!>       4| F 0 0 0 0 0         -2 -4  7 -5  8 -1 -3 -6      
!>       5| 0 0 G 0 0 0        rowStart(i)                                        
!>       6| 0 0 0 0 0 H          1 3 4 6 7 8         
!>                             columnStart(j)                                       
!>                               6 1 3 2 4 5                                   
!>      \endverbatim
!>

!>This module contains all routines dealing with (non-distributed) matrix and vectors types.
MODULE MatrixVector

  USE BaseRoutines
  USE Constants
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE Kinds
  USE Lists
  USE MatrixVectorAccessRoutines
  USE Strings
  USE Types
  USE LINKEDLIST_ROUTINES

#include "macros.h"  

  IMPLICIT NONE

  !PRIVATE

  !Module parameters

  INTEGER(INTG), PARAMETER :: bisectionToLinearSearchThreshold=10 !<Threshold for transition from bisection to linear search.

  !Module types

  !Matrix types
  
  !Module variables

  INTEGER(INTG), SAVE :: MATRIX_VECTOR_ID=1

  !Interfaces

  INTERFACE MATRIX_ALL_VALUES_SET
    MODULE PROCEDURE Matrix_AllValuesSetIntg
    MODULE PROCEDURE Matrix_AllValuesSetSP
    MODULE PROCEDURE Matrix_AllValuesSetDP
    MODULE PROCEDURE Matrix_AllValuesSetL
  END INTERFACE MATRIX_ALL_VALUES_SET

  INTERFACE Matrix_AllValuesSet
    MODULE PROCEDURE Matrix_AllValuesSetIntg
    MODULE PROCEDURE Matrix_AllValuesSetSP
    MODULE PROCEDURE Matrix_AllValuesSetDP
    MODULE PROCEDURE Matrix_AllValuesSetL
  END INTERFACE Matrix_AllValuesSet

  INTERFACE MATRIX_CREATE_FINISH
    MODULE PROCEDURE Matrix_CreateFinish
  END INTERFACE MATRIX_CREATE_FINISH

  INTERFACE MATRIX_CREATE_START
    MODULE PROCEDURE Matrix_CreateStart
  END INTERFACE MATRIX_CREATE_START
  
  INTERFACE MATRIX_DATA_TYPE_SET
    MODULE PROCEDURE Matrix_DataTypeSet
  END INTERFACE MATRIX_DATA_TYPE_SET

  INTERFACE MATRIX_NUMBER_NON_ZEROS_SET
    MODULE PROCEDURE Matrix_NumberOfNonZerosSet
  END INTERFACE MATRIX_NUMBER_NON_ZEROS_SET

  INTERFACE MATRIX_MAX_SIZE_SET
    MODULE PROCEDURE Matrix_MaxSizeSet
  END INTERFACE MATRIX_MAX_SIZE_SET

  INTERFACE MATRIX_SIZE_SET
    MODULE PROCEDURE Matrix_SizeSet
  END INTERFACE MATRIX_SIZE_SET

  INTERFACE MATRIX_STORAGE_LOCATION_FIND
    MODULE PROCEDURE Matrix_StorageLocationFind
  END INTERFACE MATRIX_STORAGE_LOCATION_FIND

  INTERFACE MATRIX_STORAGE_LOCATIONS_SET
    MODULE PROCEDURE Matrix_StorageLocationsSet
  END INTERFACE MATRIX_STORAGE_LOCATIONS_SET
  
  INTERFACE MATRIX_STORAGE_TYPE_SET
    MODULE PROCEDURE Matrix_StorageTypeSet
  END INTERFACE MATRIX_STORAGE_TYPE_SET

  INTERFACE Matrix_TransposeRowsColumnsSet
    MODULE PROCEDURE Matrix_TransposeRowsColumnsSet0
    MODULE PROCEDURE Matrix_TransposeRowsColumnsSet1
  END INTERFACE Matrix_TransposeRowsColumnsSet

  INTERFACE MATRIX_VALUES_ADD
    MODULE PROCEDURE Matrix_ValuesAddIntg
    MODULE PROCEDURE Matrix_ValuesAddIntg1
    MODULE PROCEDURE Matrix_ValuesAddIntg2
    MODULE PROCEDURE Matrix_ValuesAddSP
    MODULE PROCEDURE Matrix_ValuesAddSP1
    MODULE PROCEDURE Matrix_ValuesAddSP2
    MODULE PROCEDURE Matrix_ValuesAddDP
    MODULE PROCEDURE Matrix_ValuesAddDP1
    MODULE PROCEDURE Matrix_ValuesAddDP2
    MODULE PROCEDURE Matrix_ValuesAddL
    MODULE PROCEDURE Matrix_ValuesAddL1
    MODULE PROCEDURE Matrix_ValuesAddL2
  END INTERFACE MATRIX_VALUES_ADD

  INTERFACE Matrix_ValuesAdd
    MODULE PROCEDURE Matrix_ValuesAddIntg
    MODULE PROCEDURE Matrix_ValuesAddIntg1
    MODULE PROCEDURE Matrix_ValuesAddIntg2
    MODULE PROCEDURE Matrix_ValuesAddSP
    MODULE PROCEDURE Matrix_ValuesAddSP1
    MODULE PROCEDURE Matrix_ValuesAddSP2
    MODULE PROCEDURE Matrix_ValuesAddDP
    MODULE PROCEDURE Matrix_ValuesAddDP1
    MODULE PROCEDURE Matrix_ValuesAddDP2
    MODULE PROCEDURE Matrix_ValuesAddL
    MODULE PROCEDURE Matrix_ValuesAddL1
    MODULE PROCEDURE Matrix_ValuesAddL2
  END INTERFACE Matrix_ValuesAdd

  INTERFACE MATRIX_VALUES_GET
    MODULE PROCEDURE Matrix_ValuesGetIntg
    MODULE PROCEDURE Matrix_ValuesGetIntg1
    MODULE PROCEDURE Matrix_ValuesGetIntg2
    MODULE PROCEDURE Matrix_ValuesGetSP
    MODULE PROCEDURE Matrix_ValuesGetSP1
    MODULE PROCEDURE Matrix_ValuesGetSP2
    MODULE PROCEDURE Matrix_ValuesGetDP
    MODULE PROCEDURE Matrix_ValuesGetDP1
    MODULE PROCEDURE Matrix_ValuesGetDP2
    MODULE PROCEDURE Matrix_ValuesGetL
    MODULE PROCEDURE Matrix_ValuesGetL1
    MODULE PROCEDURE Matrix_ValuesGetL2
  END INTERFACE MATRIX_VALUES_GET
  
  INTERFACE Matrix_ValuesGet
    MODULE PROCEDURE Matrix_ValuesGetIntg
    MODULE PROCEDURE Matrix_ValuesGetIntg1
    MODULE PROCEDURE Matrix_ValuesGetIntg2
    MODULE PROCEDURE Matrix_ValuesGetSP
    MODULE PROCEDURE Matrix_ValuesGetSP1
    MODULE PROCEDURE Matrix_ValuesGetSP2
    MODULE PROCEDURE Matrix_ValuesGetDP
    MODULE PROCEDURE Matrix_ValuesGetDP1
    MODULE PROCEDURE Matrix_ValuesGetDP2
    MODULE PROCEDURE Matrix_ValuesGetL
    MODULE PROCEDURE Matrix_ValuesGetL1
    MODULE PROCEDURE Matrix_ValuesGetL2
  END INTERFACE Matrix_ValuesGet
  
  INTERFACE MATRIX_VALUES_SET
    MODULE PROCEDURE Matrix_ValuesSetIntg
    MODULE PROCEDURE Matrix_ValuesSetIntg1
    MODULE PROCEDURE Matrix_ValuesSetIntg2
    MODULE PROCEDURE Matrix_ValuesSetSP
    MODULE PROCEDURE Matrix_ValuesSetSP1
    MODULE PROCEDURE Matrix_ValuesSetSP2
    MODULE PROCEDURE Matrix_ValuesSetDP
    MODULE PROCEDURE Matrix_ValuesSetDP1
    MODULE PROCEDURE Matrix_ValuesSetDP2
    MODULE PROCEDURE Matrix_ValuesSetL
    MODULE PROCEDURE Matrix_ValuesSetL1
    MODULE PROCEDURE Matrix_ValuesSetL2
  END INTERFACE MATRIX_VALUES_SET

  INTERFACE Matrix_ValuesSet
    MODULE PROCEDURE Matrix_ValuesSetIntg
    MODULE PROCEDURE Matrix_ValuesSetIntg1
    MODULE PROCEDURE Matrix_ValuesSetIntg2
    MODULE PROCEDURE Matrix_ValuesSetSP
    MODULE PROCEDURE Matrix_ValuesSetSP1
    MODULE PROCEDURE Matrix_ValuesSetSP2
    MODULE PROCEDURE Matrix_ValuesSetDP
    MODULE PROCEDURE Matrix_ValuesSetDP1
    MODULE PROCEDURE Matrix_ValuesSetDP2
    MODULE PROCEDURE Matrix_ValuesSetL
    MODULE PROCEDURE Matrix_ValuesSetL1
    MODULE PROCEDURE Matrix_ValuesSetL2
  END INTERFACE Matrix_ValuesSet

  INTERFACE VECTOR_ALL_VALUES_SET
    MODULE PROCEDURE Vector_AllValuesSetIntg
    MODULE PROCEDURE Vector_AllValuesSetSP
    MODULE PROCEDURE Vector_AllValuesSetDP
    MODULE PROCEDURE Vector_AllValuesSetL
  END INTERFACE VECTOR_ALL_VALUES_SET

  INTERFACE Vector_AllValuesSet
    MODULE PROCEDURE Vector_AllValuesSetIntg
    MODULE PROCEDURE Vector_AllValuesSetSP
    MODULE PROCEDURE Vector_AllValuesSetDP
    MODULE PROCEDURE Vector_AllValuesSetL
  END INTERFACE Vector_AllValuesSet

  INTERFACE VECTOR_CREATE_FINISH
    MODULE PROCEDURE Vector_CreateFinish
  END INTERFACE VECTOR_CREATE_FINISH

  INTERFACE VECTOR_CREATE_START
    MODULE PROCEDURE Vector_CreateStart
  END INTERFACE VECTOR_CREATE_START

  INTERFACE VECTOR_DATA_TYPE_SET
    MODULE PROCEDURE Vector_DataTypeSet
  END INTERFACE VECTOR_DATA_TYPE_SET

  INTERFACE VECTOR_SIZE_SET
    MODULE PROCEDURE Vector_SizeSet
  END INTERFACE VECTOR_SIZE_SET

  INTERFACE VECTOR_VALUES_GET
    MODULE PROCEDURE Vector_ValuesGetIntg
    MODULE PROCEDURE Vector_ValuesGetIntg1
    MODULE PROCEDURE Vector_ValuesGetSP
    MODULE PROCEDURE Vector_ValuesGetSP1
    MODULE PROCEDURE Vector_ValuesGetDP
    MODULE PROCEDURE Vector_ValuesGetDP1
    MODULE PROCEDURE Vector_ValuesGetL
    MODULE PROCEDURE Vector_ValuesGetL1
  END INTERFACE VECTOR_VALUES_GET
  
  INTERFACE Vector_ValuesGet
    MODULE PROCEDURE Vector_ValuesGetIntg
    MODULE PROCEDURE Vector_ValuesGetIntg1
    MODULE PROCEDURE Vector_ValuesGetSP
    MODULE PROCEDURE Vector_ValuesGetSP1
    MODULE PROCEDURE Vector_ValuesGetDP
    MODULE PROCEDURE Vector_ValuesGetDP1
    MODULE PROCEDURE Vector_ValuesGetL
    MODULE PROCEDURE Vector_ValuesGetL1
  END INTERFACE Vector_ValuesGet
  
  INTERFACE VECTOR_VALUES_SET
    MODULE PROCEDURE Vector_ValuesSetIntg
    MODULE PROCEDURE Vector_ValuesSetIntg1
    MODULE PROCEDURE Vector_ValuesSetSP
    MODULE PROCEDURE Vector_ValuesSetSP1
    MODULE PROCEDURE Vector_ValuesSetDP
    MODULE PROCEDURE Vector_ValuesSetDP1
    MODULE PROCEDURE Vector_ValuesSetL
    MODULE PROCEDURE Vector_ValuesSetL1
  END INTERFACE VECTOR_VALUES_SET

  INTERFACE Vector_ValuesSet
    MODULE PROCEDURE Vector_ValuesSetIntg
    MODULE PROCEDURE Vector_ValuesSetIntg1
    MODULE PROCEDURE Vector_ValuesSetSP
    MODULE PROCEDURE Vector_ValuesSetSP1
    MODULE PROCEDURE Vector_ValuesSetDP
    MODULE PROCEDURE Vector_ValuesSetDP1
    MODULE PROCEDURE Vector_ValuesSetL
    MODULE PROCEDURE Vector_ValuesSetL1
  END INTERFACE Vector_ValuesSet

  INTERFACE MATRIX_LINKLIST_GET
    MODULE PROCEDURE Matrix_LinklistGet
  END INTERFACE MATRIX_LINKLIST_GET

  INTERFACE MATRIX_LINKLIST_SET
    MODULE PROCEDURE Matrix_LinklistSet
  END INTERFACE MATRIX_LINKLIST_SET

  PUBLIC MATRIX_VECTOR_INTG_TYPE,MATRIX_VECTOR_SP_TYPE,MATRIX_VECTOR_DP_TYPE,MATRIX_VECTOR_L_TYPE, &
    MATRIX_VECTOR_SPC_TYPE,MATRIX_VECTOR_DPC_TYPE

  PUBLIC MATRIX_BLOCK_STORAGE_TYPE,MATRIX_DIAGONAL_STORAGE_TYPE,MATRIX_COLUMN_MAJOR_STORAGE_TYPE,MATRIX_ROW_MAJOR_STORAGE_TYPE, &
    & MATRIX_COMPRESSED_ROW_STORAGE_TYPE,MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE,MATRIX_ROW_COLUMN_STORAGE_TYPE

  PUBLIC MATRIX_SYMMETRIC_TYPE,MATRIX_HERMITIAN_TYPE,MATRIX_SKEW_SYMMETRIC_TYPE,MATRIX_UNSYMMETRIC_TYPE

  PUBLIC MATRIX_ALL_VALUES_SET

  PUBLIC Matrix_AllValuesSet

  PUBLIC Matrix_BlockSizeSet

  PUBLIC MATRIX_CREATE_FINISH,MATRIX_CREATE_START

  PUBLIC Matrix_CreateStart,Matrix_CreateFinish

  PUBLIC Matrix_DataTypeSet

  PUBLIC Matrix_Destroy
   
  PUBLIC Matrix_Duplicate

  PUBLIC MATRIX_NUMBER_NON_ZEROS_SET

  PUBLIC Matrix_NumberOfNonZerosSet

  PUBLIC MATRIX_MAX_SIZE_SET

  PUBLIC Matrix_MaxSizeSet

  PUBLIC Matrix_Output

  PUBLIC MatrixRowColCoupling_Finalise,MatrixRowColCoupling_Initialise

  PUBLIC MATRIX_SIZE_SET

  PUBLIC Matrix_SizeSet

  PUBLIC MATRIX_STORAGE_LOCATION_FIND

  PUBLIC Matrix_StorageLocationFind

  PUBLIC MATRIX_STORAGE_LOCATIONS_SET

  PUBLIC Matrix_StorageLocationsSet

  PUBLIC MATRIX_STORAGE_TYPE_SET

  PUBLIC Matrix_StorageTypeSet

  PUBLIC Matrix_SymmetryTypeSet

  PUBLIC Matrix_TransposeRowColumnPositionGet

  PUBLIC Matrix_TransposeRowsColumnsSet

  PUBLIC Matrix_TransposeTypeSet
  
  PUBLIC MATRIX_VALUES_ADD

  PUBLIC Matrix_ValuesAdd

  PUBLIC MATRIX_VALUES_GET,MATRIX_VALUES_SET

  PUBLIC Matrix_ValuesGet,Matrix_ValuesSet

  PUBLIC VECTOR_ALL_VALUES_SET

  PUBLIC Vector_AllValuesSet

  PUBLIC VECTOR_CREATE_FINISH,VECTOR_CREATE_START

  PUBLIC Vector_CreateFinish,Vector_CreateStart

  PUBLIC VECTOR_DATA_TYPE_SET

  PUBLIC Vector_DataTypeSet

  PUBLIC Vector_Destroy

  PUBLIC Vector_Duplicate

  PUBLIC VECTOR_SIZE_SET

  PUBLIC Vector_SizeSet

  PUBLIC VECTOR_VALUES_GET,VECTOR_VALUES_SET

  PUBLIC Vector_ValuesGet,Vector_ValuesSet

  PUBLIC MATRIX_LINKLIST_GET,MATRIX_LINKLIST_SET

  PUBLIC Matrix_LinklistGet,Matrix_LinklistSet
  
CONTAINS
  
  !
  !================================================================================================================================
  !

  !>Sets all values in an integer matrix to the specified value.
  SUBROUTINE Matrix_AllValuesSetIntg(matrix,matrixValue,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix to set all values for
    INTEGER(INTG), INTENT(IN) :: matrixValue !<The value to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables

    ENTERS("Matrix_AllValuesSetIntg",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    CALL Matrix_AssertIsIntgData(matrix,err,error,*999)

    matrix%dataIntg=matrixValue

    EXITS("Matrix_AllValuesSetIntg")
    RETURN
999 ERRORSEXITS("Matrix_AllValuesSetIntg",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_AllValuesSetIntg

  !
  !================================================================================================================================
  !

  !>Sets all values in a single precision real matrix to the specified value.
  SUBROUTINE Matrix_AllValuesSetSP(matrix,matrixValue,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix to set all values for
    REAL(SP), INTENT(IN) :: matrixValue !<The value to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables

    ENTERS("Matrix_AllValuesSetSP",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    CALL Matrix_AssertIsSPData(matrix,err,error,*999)

    matrix%dataSP=matrixValue

    EXITS("Matrix_AllValuesSetSP")
    RETURN
999 ERRORSEXITS("Matrix_AllValuesSetSP",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_AllValuesSetSP

  !
  !================================================================================================================================
  !

  !>Sets all values in a double precision real matrix to the specified value.
  SUBROUTINE Matrix_AllValuesSetDP(matrix,matrixValue,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix to set all values for
    REAL(DP), INTENT(IN) :: matrixValue !<The value to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables

    ENTERS("Matrix_AllValuesSetDP",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    CALL Matrix_AssertIsDPData(matrix,err,error,*999)

    matrix%dataDP=matrixValue

    EXITS("Matrix_AllValuesSetDP")
    RETURN
999 ERRORSEXITS("Matrix_AllValuesSetDP",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_AllValuesSetDP

  !
  !================================================================================================================================
  !

  !>Sets all values in a logical matrix to the specified value.
  SUBROUTINE Matrix_AllValuesSetL(matrix,matrixValue,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix to set all values for
    LOGICAL, INTENT(IN) :: matrixValue !<The value to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables

    ENTERS("Matrix_AllValuesSetL",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    CALL Matrix_AssertIsLData(matrix,err,error,*999)

    matrix%dataL=matrixValue

    EXITS("Matrix_AllValuesSetL")
    RETURN
999 ERRORSEXITS("Matrix_AllValuesSetL",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_AllValuesSetL

  !
  !================================================================================================================================
  !

  !>Sets/changes the block size for a matrix.
  SUBROUTINE Matrix_BlockSizeSet(matrix,blockSize,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: blockSize !<The block size of the matrix to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Matrix_BlockSizeSet",err,error,*999)

    CALL Matrix_AssertNotFinished(matrix,err,error,*999)
    
    SELECT CASE(matrix%storageType)
    CASE(MATRIX_BLOCK_STORAGE_TYPE)
      CALL FlagError("Can not set the block size for a matrix with block storage.",err,error,*999)
    CASE(MATRIX_DIAGONAL_STORAGE_TYPE)
      CALL FlagError("Can not set the block size for a matrix with diagonal storage.",err,error,*999)
    CASE(MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
      CALL FlagError("Can not set the block size for a matrix with column major storage.",err,error,*999)          
    CASE(MATRIX_ROW_MAJOR_STORAGE_TYPE)
      CALL FlagError("Can not set the block size for a matrix with row major storage.",err,error,*999)          
    CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
      CALL FlagError("Can not set the block size for a matrix with compressed row storage.",err,error,*999)
    CASE(MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
      CALL FlagError("Can not set the block size for a matrix with compressed column storage.",err,error,*999)
    CASE(MATRIX_ROW_COLUMN_STORAGE_TYPE)
      CALL FlagError("Can not set the block size for a matrix with row column storage.",err,error,*999)      
    CASE(MATRIX_BLOCK_COMPRESSED_ROW_STORAGE_TYPE)
      IF(blockSize<1.OR.blockSize>MIN(matrix%m,matrix%n)) THEN
        localError="The specified block size of "//TRIM(NumberToVString(blockSize,"*",err,error))// &
          & " is invalid. The number must be >= 1 and <= "//TRIM(NumberToVString(MIN(matrix%M,matrix%n),"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      matrix%blockSize=blockSize
      matrix%numberOfRowBlocks=CEILING(REAL(matrix%m,DP)/REAL(matrix%blockSize,DP),INTG)
      matrix%numberOfColumnBlocks=CEILING(REAL(matrix%n,DP)/REAL(matrix%blockSize,DP),INTG)
    CASE(MATRIX_MODIFIED_KRM_STORAGE_TYPE)
      CALL FlagError("Can not set the block size for a matrix with modified KRM storage.",err,error,*999)      
    CASE DEFAULT      
      localError="The matrix storage type of "//TRIM(NumberToVString(matrix%storageType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Matrix_BlockSizeSet")
    RETURN
999 ERRORSEXITS("Matrix_BlockSizeSet",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_BlockSizeSet

  !
  !================================================================================================================================
  !

  !>Finishes the creation a matrix.
  SUBROUTINE Matrix_CreateFinish(matrix,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix to finish
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: columnIdx,count,rowIdx,rowIdx2
    TYPE(VARYING_STRING) :: localError

    ENTERS("Matrix_CreateFinish",err,error,*999)

    CALL Matrix_AssertNotFinished(matrix,err,error,*999)
    
    SELECT CASE(matrix%storageType)
    CASE(MATRIX_BLOCK_STORAGE_TYPE)
      IF(matrix%maxM==-1) matrix%maxM=matrix%m
      IF(matrix%maxN==-1) matrix%maxN=matrix%n
      matrix%size=matrix%m*matrix%n
      matrix%numberOfNonZeros=matrix%m*matrix%n
      matrix%maximumColumnIndicesPerRow=matrix%n
      matrix%blockSize=MAX(matrix%m,matrix%n)
      matrix%numberOfBlocks=1
      matrix%numberOfRowBlocks=1
      matrix%numberOfColumnBlocks=1
    CASE(MATRIX_DIAGONAL_STORAGE_TYPE)
      IF(matrix%maxM==-1) matrix%maxM=matrix%m
      IF(matrix%maxN==-1) matrix%maxN=matrix%n
      matrix%size=matrix%m
      matrix%numberOfNonZeros=matrix%m
      matrix%maximumColumnIndicesPerRow=1
      matrix%blockSize=MAX(matrix%m,matrix%n)
      matrix%numberOfBlocks=1
      matrix%numberOfRowBlocks=1
      matrix%numberOfColumnBlocks=1
    CASE(MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
      IF(matrix%maxM==-1) CALL FlagError("Maximum number of rows has not been set for this matrix.",err,error,*999)
      IF(matrix%maxN==-1) CALL FlagError("Maximum number of columns has not been set for this matrix.",err,error,*999)
      matrix%size=matrix%maxM*matrix%n
      matrix%numberOfNonZeros=matrix%m*matrix%n
      matrix%maximumColumnIndicesPerRow=matrix%n
      matrix%blockSize=MAX(matrix%m,matrix%n)
      matrix%numberOfBlocks=1
      matrix%numberOfRowBlocks=1
      matrix%numberOfColumnBlocks=1
    CASE(MATRIX_ROW_MAJOR_STORAGE_TYPE)
      IF(matrix%maxM==-1) CALL FlagError("Maximum number of rows has not been set for this matrix.",err,error,*999)
      IF(matrix%maxN==-1) CALL FlagError("Maximum number of columns has not been set for this matrix.",err,error,*999)
      matrix%size=matrix%m*matrix%maxN
      matrix%numberOfNonZeros=matrix%m*matrix%n
      matrix%maximumColumnIndicesPerRow=matrix%n
      matrix%blockSize=MAX(matrix%m,matrix%n)
      matrix%numberOfBlocks=1
      matrix%numberOfRowBlocks=1
      matrix%numberOfColumnBlocks=1
    CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
      IF(matrix%numberOfNonZeros==-1) CALL FlagError("Number of non-zeros has not been set for this matrix.",err,error,*999)
      IF(matrix%maxM==-1) matrix%maxM=matrix%m
      IF(matrix%maxN==-1) matrix%maxN=matrix%n
      matrix%size=matrix%numberOfNonZeros
      IF(.NOT.ALLOCATED(matrix%columnIndices))  &
        & CALL FlagError("Matrix storage locations column indices have not been set.",err,error,*999)
      IF(.NOT.ALLOCATED(matrix%rowIndices))  &
        & CALL FlagError("Matrix storage locations row indices have not been set.",err,error,*999)
      matrix%maximumColumnIndicesPerRow=0
      DO rowIdx=1,matrix%m
        IF((matrix%rowIndices(rowIdx+1)-matrix%rowIndices(rowIdx))>matrix%maximumColumnIndicesPerRow) &
          & matrix%maximumColumnIndicesPerRow=matrix%rowIndices(rowIdx+1)-matrix%rowIndices(rowIdx)
      ENDDO !rowIdx
      matrix%blockSize=MAX(matrix%m,matrix%n)
      matrix%numberOfBlocks=1
      matrix%numberOfRowBlocks=1
      matrix%numberOfColumnBlocks=1
    CASE(MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
      IF(matrix%numberOfNonZeros==-1) CALL FlagError("Number of non-zeros has not been set for this matrix.",err,error,*999)
      IF(matrix%maxM==-1) matrix%maxM=matrix%m
      IF(matrix%maxN==-1) matrix%maxN=matrix%n
      matrix%size=matrix%numberOfNonZeros
      IF(.NOT.ALLOCATED(matrix%columnIndices))  &
        & CALL FlagError("Matrix storage locations column indices have not been set.",err,error,*999)
      IF(.NOT.ALLOCATED(matrix%rowIndices))  &
        & CALL FlagError("Matrix storage locations row indices have not been set.",err,error,*999)
      matrix%maximumColumnIndicesPerRow=0
      DO rowIdx=1,matrix%m
        count=0
        DO columnIdx=1,matrix%n
          DO rowIdx2=matrix%columnIndices(columnIdx),matrix%columnIndices(columnIdx+1)-1
            IF(matrix%rowIndices(rowIdx2)==rowIdx) count=count+1
          ENDDO !rowIdx2
        ENDDO !columnIdx
        IF(count>matrix%maximumColumnIndicesPerRow) matrix%maximumColumnIndicesPerRow=count
      ENDDO !rowIdx
      matrix%blockSize=MAX(matrix%m,matrix%n)
      matrix%numberOfBlocks=1
      matrix%numberOfRowBlocks=1
      matrix%numberOfColumnBlocks=1
    CASE(MATRIX_ROW_COLUMN_STORAGE_TYPE)
      IF(matrix%numberOfNonZeros==-1) CALL FlagError("Number of non-zeros has not been set for this matrix.",err,error,*999)
      IF(matrix%maxM==-1) matrix%maxM=matrix%m
      IF(matrix%maxN==-1) matrix%maxN=matrix%n
      matrix%size=matrix%numberOfNonZeros  
      IF(.NOT.ALLOCATED(matrix%columnIndices))  &
        & CALL FlagError("Matrix storage locations column indices have not been set.",err,error,*999)
      IF(.NOT.ALLOCATED(matrix%rowIndices))  &
        & CALL FlagError("Matrix storage locations row indices have not been set.",err,error,*999)
      matrix%maximumColumnIndicesPerRow=0
      DO rowIdx=1,matrix%m
        count=0
        DO rowIdx2=1,matrix%numberOfNonZeros
          IF(matrix%rowIndices(rowIdx2)==rowIdx) count=count+1
        ENDDO !rowIdx2            
        IF(count>matrix%maximumColumnIndicesPerRow) matrix%maximumColumnIndicesPerRow=count
      ENDDO !rowIdx          
      matrix%blockSize=MAX(matrix%m,matrix%n)
      matrix%numberOfBlocks=1
      matrix%numberOfRowBlocks=1
      matrix%numberOfColumnBlocks=1
    CASE(MATRIX_BLOCK_COMPRESSED_ROW_STORAGE_TYPE)
      IF(.NOT.ALLOCATED(matrix%columnIndices))  &
        & CALL FlagError("Matrix storage locations column indices have not been set.",err,error,*999)
      IF(.NOT.ALLOCATED(matrix%rowIndices))  &
        & CALL FlagError("Matrix storage locations row indices have not been set.",err,error,*999)
      matrix%numberOfRowBlocks=CEILING(REAL(matrix%m,DP)/REAL(matrix%blockSize,DP),INTG)
      matrix%numberOfColumnBlocks=CEILING(REAL(matrix%n,DP)/REAL(matrix%blockSize,DP),INTG)
      matrix%maxM=matrix%numberOfRowBlocks*matrix%blockSize
      matrix%maxN=matrix%numberOfColumnBlocks*matrix%blockSize
      matrix%numberOfNonZeros=matrix%numberOfBlocks*matrix%blockSize*matrix%blockSize
      matrix%size=matrix%numberOfNonZeros
      matrix%maximumColumnIndicesPerRow=0
      DO rowIdx=1,matrix%numberOfRowBlocks
        IF((matrix%rowIndices(rowIdx+1)-matrix%rowIndices(rowIdx))>matrix%maximumColumnIndicesPerRow) &
          & matrix%maximumColumnIndicesPerRow=matrix%rowIndices(rowIdx+1)-matrix%rowIndices(rowIdx)
      ENDDO !rowIdx
      matrix%maximumColumnIndicesPerRow=matrix%maximumColumnIndicesPerRow*matrix%blockSize
    CASE(MATRIX_MODIFIED_KRM_STORAGE_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The matrix storage type of "//TRIM(NumberToVString(matrix%storageType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    IF(matrix%size>0) THEN
      SELECT CASE(matrix%dataType)
      CASE(MATRIX_VECTOR_INTG_TYPE)
        ALLOCATE(matrix%dataIntg(matrix%size),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate matrix integer data.",err,error,*999)
      CASE(MATRIX_VECTOR_SP_TYPE)
        ALLOCATE(matrix%dataSP(matrix%size),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate matrix single precision data.",err,error,*999)
      CASE(MATRIX_VECTOR_DP_TYPE)
        ALLOCATE(matrix%dataDP(matrix%size),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate matrix double precision data.",err,error,*999)
      CASE(MATRIX_VECTOR_L_TYPE)
        ALLOCATE(matrix%dataL(matrix%size),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate matrix logical data.",err,error,*999)
      CASE(MATRIX_VECTOR_SPC_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(MATRIX_VECTOR_DPC_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The matrix data type of "//TRIM(NumberToVString(matrix%dataType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDIF
    matrix%ID=MATRIX_VECTOR_ID
    MATRIX_VECTOR_ID=MATRIX_VECTOR_ID+1
    matrix%matrixFinished=.TRUE.
    IF(matrix%transposeType==MATRIX_FULL_TRANSPOSE_REQUIRED) CALL Matrix_TransposeLocationsCalculate(matrix,err,error,*999)
    
    EXITS("Matrix_CreateFinish")
    RETURN
!!TODO: deallocate on error
999 ERRORSEXITS("Matrix_CreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_CreateFinish

  !
  !================================================================================================================================
  !

  !>Starts the creation a matrix.
  SUBROUTINE Matrix_CreateStart(matrix,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Matrix_CreateStart",err,error,*998)

    IF(ASSOCIATED(matrix)) CALL FlagError("Matrix is already associated.",err,error,*998)
    
    ALLOCATE(matrix,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate the matrix.",err,error,*999)
    CALL Matrix_Initialise(matrix,err,error,*999)
    !Set the defaults
    matrix%dataType=MATRIX_VECTOR_DP_TYPE
    matrix%storageType=MATRIX_BLOCK_STORAGE_TYPE
    
    EXITS("Matrix_CreateStart")
    RETURN
999 IF(ASSOCIATED(matrix)) CALL Matrix_Finalise(matrix,err,error,*998)
998 ERRORSEXITS("Matrix_CreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_CreateStart

  !
  !================================================================================================================================
  !

  !>Sets/changes the data type of a matrix.
  SUBROUTINE Matrix_DataTypeSet(matrix,dataType,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix to set the data type for
    INTEGER(INTG), INTENT(IN) :: dataType !<The data type to set for the matrix. \see MatrixVector_DataTypes,MatrixVector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Matrix_DataTypeSet",err,error,*999)

    CALL Matrix_AssertNotFinished(matrix,err,error,*999)
    
    SELECT CASE(dataType)
    CASE(MATRIX_VECTOR_INTG_TYPE)
      matrix%dataType=MATRIX_VECTOR_INTG_TYPE
    CASE(MATRIX_VECTOR_SP_TYPE)
      matrix%dataType=MATRIX_VECTOR_SP_TYPE
    CASE(MATRIX_VECTOR_DP_TYPE)
      matrix%dataType=MATRIX_VECTOR_DP_TYPE
    CASE(MATRIX_VECTOR_L_TYPE)
      matrix%dataType=MATRIX_VECTOR_L_TYPE
    CASE(MATRIX_VECTOR_SPC_TYPE)
      matrix%dataType=MATRIX_VECTOR_SPC_TYPE
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(MATRIX_VECTOR_DPC_TYPE)
      matrix%dataType=MATRIX_VECTOR_DPC_TYPE
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The matrix vector data type of "//TRIM(NumberToVString(dataType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Matrix_DataTypeSet")
    RETURN
999 ERRORSEXITS("Matrix_DataTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_DataTypeSet

  !
  !================================================================================================================================
  !

  !>Destroys a matrix 
  SUBROUTINE Matrix_Destroy(matrix,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("Matrix_Destroy",err,error,*999)

    IF(.NOT.ASSOCIATED(matrix)) CALL FlagError("Matrix is not associated.",err,error,*999)
    
    CALL Matrix_Finalise(matrix,err,error,*999)

    EXITS("Matrix_Destroy")
    RETURN
999 ERRORSEXITS("Matrix_Destroy",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_Destroy

  !
  !================================================================================================================================
  !

  !>Duplicates the matrix and returns a pointer to the duplicated matrix in newMatrix.
  SUBROUTINE Matrix_Duplicate(matrix,newMatrix,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix to duplicate
    TYPE(MatrixType), POINTER :: newMatrix !<On return a pointer to a new duplicated matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError,localError
    
    ENTERS("Matrix_Duplicate",err,error,*998)

    IF(ASSOCIATED(newMatrix)) CALL FlagError("New matrix is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(matrix)) CALL FlagError("Matrix is not associated.",err,error,*998)
    
    CALL Matrix_CreateStart(newMatrix,err,error,*999)
    CALL Matrix_DataTypeSet(newMatrix,matrix%dataType,err,error,*999)
    CALL Matrix_SizeSet(newMatrix,matrix%m,matrix%n,err,error,*999)
    CALL Matrix_StorageTypeSet(newMatrix,matrix%storageType,err,error,*999)
    CALL Matrix_SymmetryTypeSet(newMatrix,matrix%symmetryType,err,error,*999)
    CALL Matrix_TransposeTypeSet(newMatrix,matrix%transposeType,err,error,*999)
    SELECT CASE(matrix%storageType)
    CASE(MATRIX_BLOCK_STORAGE_TYPE,MATRIX_DIAGONAL_STORAGE_TYPE)
      !Do nothing
    CASE(MATRIX_COLUMN_MAJOR_STORAGE_TYPE,MATRIX_ROW_MAJOR_STORAGE_TYPE)
      CALL Matrix_MaxSizeSet(newMatrix,matrix%maxM,matrix%maxN,err,error,*999)          
    CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE,MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE,MATRIX_ROW_COLUMN_STORAGE_TYPE)
      CALL Matrix_NumberOfNonZerosSet(newMatrix,matrix%numberOfNonZeros,err,error,*999)
      CALL Matrix_StorageLocationsSet(newMatrix,matrix%rowIndices,matrix%columnIndices,err,error,*999)
    CASE(MATRIX_BLOCK_COMPRESSED_ROW_STORAGE_TYPE)
      CALL Matrix_BlockSizeSet(newMatrix,matrix%blockSize,err,error,*999)
      CALL Matrix_StorageLocationsSet(newMatrix,matrix%rowIndices,matrix%columnIndices,err,error,*999)
    CASE(MATRIX_MODIFIED_KRM_STORAGE_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The matrix storage type of "//TRIM(NumberToVString(matrix%storageType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    CALL Matrix_CreateFinish(newMatrix,err,error,*999)

    EXITS("Matrix_Duplicate")
    RETURN
999 CALL Matrix_Finalise(newMatrix,dummyErr,dummyError,*998)
998 ERRORSEXITS("Matrix_Duplicate",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_Duplicate

  !
  !================================================================================================================================
  !

  !>Finalises a matrix and deallocates all memory.
  SUBROUTINE Matrix_Finalise(matrix,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Matrix_Finalise",err,error,*999)

    IF(ASSOCIATED(matrix)) THEN
      IF(ALLOCATED(matrix%rowIndices)) DEALLOCATE(matrix%rowIndices)
      IF(ALLOCATED(matrix%columnIndices)) DEALLOCATE(matrix%columnIndices)
      IF(ALLOCATED(matrix%rowIndicesT)) DEALLOCATE(matrix%rowIndicesT)
      IF(ALLOCATED(matrix%columnIndicesT)) DEALLOCATE(matrix%columnIndicesT)
      IF(ALLOCATED(matrix%dataIntg)) DEALLOCATE(matrix%dataIntg)
      IF(ALLOCATED(matrix%dataSP)) DEALLOCATE(matrix%dataSP)
      IF(ALLOCATED(matrix%dataDP)) DEALLOCATE(matrix%dataDP)
      IF(ALLOCATED(matrix%dataL)) DEALLOCATE(matrix%dataL)
      DEALLOCATE(matrix)
    ENDIF

    EXITS("Matrix_Finalise")
    RETURN
999 ERRORSEXITS("Matrix_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises a matrix.
  SUBROUTINE Matrix_Initialise(matrix,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Matrix_Initialise",err,error,*999)

    IF(.NOT.ASSOCIATED(matrix)) CALL FlagError("Matrix is not associated.",err,error,*999)
    
!!TODO: have a matrix user number etc.
    matrix%ID=0
    matrix%matrixFinished=.FALSE.
    matrix%m=0
    matrix%n=0
    matrix%maxM=-1
    matrix%maxN=-1
    matrix%dataType=0
    matrix%storageType=0
    matrix%symmetryType=MATRIX_UNSYMMETRIC_TYPE !Should this default to symmetric???
    matrix%transposeType=MATRIX_NO_TRANSPOSE_REQUIRED
    matrix%numberOfNonZeros=0
    matrix%size=0      
    matrix%maximumColumnIndicesPerRow=0
    matrix%blockSize=2
    matrix%numberOfBlocks=0
    matrix%numberOfRowBlocks=0
    matrix%numberOfColumnBlocks=0

    EXITS("Matrix_Initialise")
    RETURN
999 ERRORSEXITS("Matrix_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_Initialise

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of non zeros for a matrix.
  SUBROUTINE Matrix_NumberOfNonZerosSet(matrix,numberOfNonZeros,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: numberOfNonZeros !<The number of non zeros in the matrix to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Matrix_NumberOfNonZerosSet",err,error,*999)

    CALL Matrix_AssertNotFinished(matrix,err,error,*999)
    
    SELECT CASE(matrix%storageType)
    CASE(MATRIX_BLOCK_STORAGE_TYPE)
      CALL FlagError("Can not set the number of non-zeros for a matrix with block storage.",err,error,*999)
    CASE(MATRIX_DIAGONAL_STORAGE_TYPE)
      CALL FlagError("Can not set the number of non-zeros for a matrix with diagonal storage.",err,error,*999)
    CASE(MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
      CALL FlagError("Can not set the number of non-zeros for a matrix with column major storage.",err,error,*999)          
    CASE(MATRIX_ROW_MAJOR_STORAGE_TYPE)
      CALL FlagError("Can not set the number of non-zeros for a matrix with row major storage.",err,error,*999)          
    CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE,MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE,MATRIX_ROW_COLUMN_STORAGE_TYPE)
      IF(numberOfNonZeros>=0) THEN
        matrix%numberOfNonZeros=numberOfNonZeros
      ELSE
        localError="The number of non-zeros ("//TRIM(NumberToVString(numberOfNonZeros,"*",err,error))// &
          & ") is invalid. The number must be greater than or equal to zero."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE(MATRIX_BLOCK_COMPRESSED_ROW_STORAGE_TYPE)
      CALL FlagError("Can not set the number of non-zeros for a matrix with block compressed row storage.",err,error,*999)
    CASE(MATRIX_MODIFIED_KRM_STORAGE_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The matrix storage type of "//TRIM(NumberToVString(matrix%storageType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Matrix_NumberOfNonZerosSet")
    RETURN
999 ERRORSEXITS("Matrix_NumberOfNonZerosSet",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_NumberOfNonZerosSet

  !
  !================================================================================================================================
  !
  
  !>Gets the maximum number of columns in each row of a distributed matrix.
  SUBROUTINE Matrix_LinkListSet(matrix,list,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    TYPE(LinkedList),pointer :: list(:) 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("Matrix_LinkListSet",err,error,*999)

    CALL Matrix_AssertNotFinished(matrix,err,error,*999)
    
    matrix%list=>list      
    
    EXITS("Matrix_LinkListSet")
    RETURN
999 ERRORSEXITS("Matrix_LinkListSet",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_LinkListSet

   !
  !================================================================================================================================
  !!>Gets the maximum number of columns in each row of a distributed matrix.
  SUBROUTINE Matrix_LinkListGet(matrix,list,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    type(LinkedList),pointer :: list(:) 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("Matrix_LinkListGet",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    
    list=>matrix%list
     
    EXITS("Matrix_LinkListGet")
    RETURN
999 ERRORSEXITS("Matrix_LinkListGet",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_LinkListGet

  !
  !================================================================================================================================
  !
  
  !>Sets/changes the maximum size of a matrix.
  SUBROUTINE Matrix_MaxSizeSet(matrix,maxM,maxN,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: maxM !<The maximum number of rows to set
    INTEGER(INTG), INTENT(IN) :: maxN !<The maximum number of columns to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Matrix_MaxSizeSet",err,error,*999)

    CALL Matrix_AssertNotFinished(matrix,err,error,*999)
    IF(maxM<=0) THEN
      localError="The maximum number of matrix rows of "//TRIM(NumberToVString(maxM,"*",err,error))// &
        & " is invalid. The number must be > 0."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(maxN<=0) THEN
      localError="The maximum number of matrix columns of "//TRIM(NumberToVString(maxN,"*",err,error))// &
        & " is invalid. The number must be > 0."
      CALL FlagError(localError,err,error,*999)
    ENDIF  
    IF(maxM<matrix%m) THEN
      localError="The maximum number of matrix columns of "//TRIM(NumberToVString(maxM,"*",err,error))// &
        & " must be >= the number of matrix columns of "//TRIM(NumberToVString(matrix%m,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(maxN<matrix%n) THEN
      localError="The maximum number of matrix rows of "//TRIM(NumberToVString(maxN,"*",err,error))// &
        & " must be >= the number of matrix rows of "//TRIM(NumberToVString(matrix%n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    matrix%maxM=maxM
    matrix%maxN=maxN

    EXITS("Matrix_MaxSizeSet")
    RETURN
999 ERRORSEXITS("Matrix_MaxSizeSet",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_MaxSizeSet
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the size of a matrix.
  SUBROUTINE Matrix_Output(ID,matrix,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID to output to
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: columnIdx,columnBlockIdx,columnBlockNumber,dataIdx,entryIdx,rowIdx,rowBlockIdx,rowNumber
    INTEGER(INTG), ALLOCATABLE :: rowEntriesIntg(:)
    REAL(SP), ALLOCATABLE :: rowEntriesSP(:)
    REAL(DP), ALLOCATABLE :: rowEntriesDP(:)
    LOGICAL, ALLOCATABLE :: rowEntriesL(:)
    CHARACTER(LEN=9) :: rowString,columnString
    CHARACTER(LEN=39) :: initialString
    TYPE(VARYING_STRING) :: localError

    ENTERS("Matrix_Output",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    
    SELECT CASE(matrix%storageType)
    CASE(MATRIX_BLOCK_STORAGE_TYPE)
      SELECT CASE(matrix%dataType)
      CASE(MATRIX_VECTOR_INTG_TYPE)
        CALL WriteStringMatrix(ID,1,1,matrix%m,1,1,matrix%n,8,8,RESHAPE(matrix%dataIntg,[matrix%maxM,matrix%maxN]), &
          & WRITE_STRING_MATRIX_NAME_AND_INDICES,'("Matrix','(",I9,",:)',':",8(X,I13))','(20X,8(X,I13))', &
          & err,error,*999)
      CASE(MATRIX_VECTOR_SP_TYPE)
        CALL WriteStringMatrix(ID,1,1,matrix%m,1,1,matrix%n,8,8,RESHAPE(matrix%dataSP,[matrix%maxM,matrix%maxN]), &
          & WRITE_STRING_MATRIX_NAME_AND_INDICES,'("Matrix','(",I9,",:)',':",8(X,E13.6))','(20X,8(X,E13.6))', &
          & err,error,*999)
      CASE(MATRIX_VECTOR_DP_TYPE)
        CALL WriteStringMatrix(ID,1,1,matrix%m,1,1,matrix%n,8,8,RESHAPE(matrix%dataDP,[matrix%maxM,matrix%maxN]), &
          & WRITE_STRING_MATRIX_NAME_AND_INDICES,'("Matrix','(",I9,",:)',':",8(X,E13.6))','(20X,8(X,E13.6))', &
          & err,error,*999)
      CASE(MATRIX_VECTOR_L_TYPE)            
        CALL WriteStringMatrix(ID,1,1,matrix%m,1,1,matrix%n,8,8,RESHAPE(matrix%dataL,[matrix%maxM,matrix%maxN]), &
          & WRITE_STRING_MATRIX_NAME_AND_INDICES,'("Matrix','(",I9,",:)',':",8(X,L13))','(20X,8(X,L13))', &
          & err,error,*999)
      CASE(MATRIX_VECTOR_SPC_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(MATRIX_VECTOR_DPC_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The matrix data type of "//TRIM(NumberToVString(matrix%dataType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(MATRIX_DIAGONAL_STORAGE_TYPE)
      DO rowIdx=1,matrix%m
        rowString=NumberToCharacter(rowIdx,"I9",err,error)
        initialString='Matrix('//rowString//','//rowString//'):'
        SELECT CASE(matrix%dataType)              
        CASE(MATRIX_VECTOR_INTG_TYPE)
          CALL WriteStringFmtValue(ID,initialString(1:LEN_TRIM(initialString)),matrix%dataIntg(rowIdx),'(7X,I13)', &
            & err,error,*999)
        CASE(MATRIX_VECTOR_SP_TYPE)
          CALL WriteStringFmtValue(ID,initialString(1:LEN_TRIM(initialString)),matrix%dataSP(rowIdx),'(7X,E13.6)', &
            & err,error,*999)
        CASE(MATRIX_VECTOR_DP_TYPE)
          CALL WriteStringFmtValue(ID,initialString(1:LEN_TRIM(initialString)),matrix%dataDP(rowIdx),'(7X,E13.6)', &
            & err,error,*999)
        CASE(MATRIX_VECTOR_L_TYPE)            
          CALL WriteStringFmtValue(ID,initialString(1:LEN_TRIM(initialString)),matrix%dataL(rowIdx),'(7X,L13)', &
            & err,error,*999)
        CASE(MATRIX_VECTOR_SPC_TYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(MATRIX_VECTOR_DPC_TYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE DEFAULT
          localError="The matrix data type of "//TRIM(NumberToVString(matrix%dataType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDDO !rowIdx
    CASE(MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(MATRIX_ROW_MAJOR_STORAGE_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
      DO rowIdx=1,matrix%m
        rowString=NumberToCharacter(rowIdx,"I9",err,error)
        SELECT CASE(matrix%dataType)              
        CASE(MATRIX_VECTOR_INTG_TYPE)
          initialString='("Matrix('//rowString//',:):",8(X,I13))'
          CALL WriteStringVector(ID,matrix%rowIndices(rowIdx),1,matrix%rowIndices(rowIdx+1)-1,8,8,matrix%dataIntg,initialString, &
            & '(20X,8(X,I13))',err,error,*999)
        CASE(MATRIX_VECTOR_SP_TYPE)
          initialString='("Matrix('//rowString//',:):",8(X,E13.6))'
          CALL WriteStringVector(ID,matrix%rowIndices(rowIdx),1,matrix%rowIndices(rowIdx+1)-1,8,8,matrix%dataSP,initialString, &
            & '(20X,8(X,E13.6))',err,error,*999)
        CASE(MATRIX_VECTOR_DP_TYPE)
          initialString='("Matrix('//rowString//',:):",8(X,E13.6))'
          CALL WriteStringVector(ID,matrix%rowIndices(rowIdx),1,matrix%rowIndices(rowIdx+1)-1,8,8,matrix%dataDP,initialString, &
            & '(20X,8(X,E13.6))',err,error,*999)
        CASE(MATRIX_VECTOR_L_TYPE)            
          initialString='("Matrix('//rowString//',:):",8(X,L13))'
          CALL WriteStringVector(ID,matrix%rowIndices(rowIdx),1,matrix%rowIndices(rowIdx+1)-1,8,8,matrix%dataL,initialString, &
            & '(20X,8(X,L13))',err,error,*999)
        CASE(MATRIX_VECTOR_SPC_TYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(MATRIX_VECTOR_DPC_TYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE DEFAULT
          localError="The matrix data type of "//TRIM(NumberToVString(matrix%dataType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDDO !rowIdx
    CASE(MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
      DO columnIdx=1,matrix%n
        columnString=NumberToCharacter(columnIdx,"I9",err,error)
        SELECT CASE(matrix%dataType)              
        CASE(MATRIX_VECTOR_INTG_TYPE)
          initialString='("Matrix(:,'//columnString//'):",8(X,I13))'
          CALL WriteStringVector(ID,matrix%columnIndices(columnIdx),1,matrix%columnIndices(columnIdx+1)-1,8,8,matrix%dataIntg, &
            & initialString,'(20X,8(X,I13))',err,error,*999)
        CASE(MATRIX_VECTOR_SP_TYPE)
          initialString='("Matrix(:,'//columnString//'):",8(X,E13.6))'
          CALL WriteStringVector(ID,matrix%columnIndices(columnIdx),1,matrix%columnIndices(columnIdx+1)-1,8,8,matrix%dataSP, &
            & initialString,'(20X,8(X,E13.6))',err,error,*999)
        CASE(MATRIX_VECTOR_DP_TYPE)
          initialString='("Matrix(:,'//columnString//'):",8(X,E13.6))'
          CALL WriteStringVector(ID,matrix%columnIndices(columnIdx),1,matrix%columnIndices(columnIdx+1)-1,8,8,matrix%dataDP, &
            & initialString,'(20X,8(X,E13.6))',err,error,*999)
        CASE(MATRIX_VECTOR_L_TYPE)            
          initialString='("Matrix(:,'//columnString//'):",8(X,L13))'
          CALL WriteStringVector(ID,matrix%columnIndices(columnIdx),1,matrix%columnIndices(columnIdx+1)-1,8,8,matrix%dataL, &
            & initialString,'(20X,8(X,L13))',err,error,*999)
        CASE(MATRIX_VECTOR_SPC_TYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(MATRIX_VECTOR_DPC_TYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE DEFAULT
          localError="The matrix data type of "//TRIM(NumberToVString(matrix%dataType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDDO !columnIdx
    CASE(MATRIX_ROW_COLUMN_STORAGE_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(MATRIX_BLOCK_COMPRESSED_ROW_STORAGE_TYPE)
      !Create temporary array to store the row
      SELECT CASE(matrix%dataType)
      CASE(MATRIX_VECTOR_INTG_TYPE)
        initialString='("Matrix(:,'//columnString//'):",8(X,I13))'
        ALLOCATE(rowEntriesIntg(matrix%maximumColumnIndicesPerRow),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate integer row entries.",err,error,*999)
      CASE(MATRIX_VECTOR_SP_TYPE)
        initialString='("Matrix(:,'//columnString//'):",8(X,E13.6))'
        ALLOCATE(rowEntriesSP(matrix%maximumColumnIndicesPerRow),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate single precision real row entries.",err,error,*999)
      CASE(MATRIX_VECTOR_DP_TYPE)
        initialString='("Matrix(:,'//columnString//'):",8(X,E13.6))'
        ALLOCATE(rowEntriesDP(matrix%maximumColumnIndicesPerRow),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate double precision real row entries.",err,error,*999)
      CASE(MATRIX_VECTOR_L_TYPE)            
        initialString='("Matrix(:,'//columnString//'):",8(X,L13))'
        ALLOCATE(rowEntriesL(matrix%maximumColumnIndicesPerRow),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate logical row entries.",err,error,*999)
      CASE(MATRIX_VECTOR_SPC_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(MATRIX_VECTOR_DPC_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The matrix data type of "//TRIM(NumberToVString(matrix%dataType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT  
      DO rowBlockIdx=1,matrix%numberOfRowBlocks
        DO rowIdx=1,matrix%blockSize
          rowNumber=rowIdx+(rowBlockIdx-1)*matrix%blockSize
          rowString=NumberToCharacter(rowNumber,"I9",err,error)
          entryIdx=0
          DO columnBlockIdx=matrix%rowIndices(rowBlockIdx),matrix%rowIndices(rowBlockIdx+1)-1
            columnBlockNumber=matrix%columnIndices(columnBlockIdx)
            DO columnIdx=1,matrix%blockSize
              dataIdx=columnIdx+(columnBlockNumber-1)*matrix%blockSize
              entryIdx=entryIdx+1
              SELECT CASE(matrix%dataType)              
              CASE(MATRIX_VECTOR_INTG_TYPE)
                rowEntriesIntg(entryIdx)=matrix%dataIntg(dataIdx)
              CASE(MATRIX_VECTOR_SP_TYPE)
                rowEntriesSP(entryIdx)=matrix%dataSP(dataIdx)
              CASE(MATRIX_VECTOR_DP_TYPE)
                rowEntriesDP(entryIdx)=matrix%dataDP(dataIdx)
              CASE(MATRIX_VECTOR_L_TYPE)            
                rowEntriesL(entryIdx)=matrix%dataL(dataIdx)
              CASE(MATRIX_VECTOR_SPC_TYPE)
                CALL FlagError("Not implemented.",err,error,*999)
              CASE(MATRIX_VECTOR_DPC_TYPE)
                CALL FlagError("Not implemented.",err,error,*999)
              CASE DEFAULT
                localError="The matrix data type of "//TRIM(NumberToVString(matrix%dataType,"*",err,error))//" is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ENDDO !columnIdx
          ENDDO !columnBlockIdx
          SELECT CASE(matrix%dataType)              
          CASE(MATRIX_VECTOR_INTG_TYPE)
            CALL WriteStringVector(ID,1,1,entryIdx,8,8,rowEntriesIntg,initialString,'(20X,8(X,I13))',err,error,*999)
          CASE(MATRIX_VECTOR_SP_TYPE)
            CALL WriteStringVector(ID,1,1,entryIdx,8,8,rowEntriesSP,initialString,'(20X,8(X,E13.6))',err,error,*999)
          CASE(MATRIX_VECTOR_DP_TYPE)
            CALL WriteStringVector(ID,1,1,entryIdx,8,8,rowEntriesDP,initialString,'(20X,8(X,E13.6))',err,error,*999)
          CASE(MATRIX_VECTOR_L_TYPE)            
            CALL WriteStringVector(ID,1,1,entryIdx,8,8,rowEntriesL,initialString,'(20X,8(X,L13))',err,error,*999)
          CASE(MATRIX_VECTOR_SPC_TYPE)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(MATRIX_VECTOR_DPC_TYPE)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The matrix data type of "//TRIM(NumberToVString(matrix%dataType,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT          
        ENDDO ! rowIdx
      ENDDO !rowBlockIdx
      IF(ALLOCATED(rowEntriesIntg)) DEALLOCATE(rowEntriesIntg)
      IF(ALLOCATED(rowEntriesSP)) DEALLOCATE(rowEntriesSP)
      IF(ALLOCATED(rowEntriesDP)) DEALLOCATE(rowEntriesDP)
      IF(ALLOCATED(rowEntriesL)) DEALLOCATE(rowEntriesL)
    CASE(MATRIX_MODIFIED_KRM_STORAGE_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The matrix storage type of "//TRIM(NumberToVString(matrix%storageType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Matrix_Output")
    RETURN
999 IF(ALLOCATED(rowEntriesIntg)) DEALLOCATE(rowEntriesIntg)
    IF(ALLOCATED(rowEntriesSP)) DEALLOCATE(rowEntriesSP)
    IF(ALLOCATED(rowEntriesDP)) DEALLOCATE(rowEntriesDP)
    IF(ALLOCATED(rowEntriesL)) DEALLOCATE(rowEntriesL)
    ERRORSEXITS("Matrix_Output",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_Output

  !
  !================================================================================================================================
  !

  !>Finalise the matrix row/column coupling 
  SUBROUTINE MatrixRowColCoupling_Finalise(rowColCoupling,err,error,*)

    !Argument variables
    TYPE(MatrixRowColCouplingType), INTENT(INOUT) :: rowColCoupling !<The row or column couplings to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("MatrixRowColCoupling_Finalise",err,error,*999)

    IF(ALLOCATED(rowColCoupling%rowCols)) DEALLOCATE(rowColCoupling%rowCols)
    IF(ALLOCATED(rowColCoupling%couplingCoefficients)) DEALLOCATE(rowColCoupling%couplingCoefficients)
    rowColCoupling%numberOfRowCols=0

    EXITS("MatrixRowColCoupling_Finalise")
    RETURN
999 ERRORSEXITS("MatrixRowColCoupling_Finalise",err,error)
    RETURN 1

  END SUBROUTINE MatrixRowColCoupling_Finalise

  !
  !================================================================================================================================
  !

  !>Initialise the matrix row/column coupling
  SUBROUTINE MatrixRowColCoupling_Initialise(rowColCoupling,err,error,*)

    !Argument variables
    TYPE(MatrixRowColCouplingType), INTENT(INOUT) :: rowColCoupling !<The row/column coupling to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("MatrixRowColCoupling_Initialise",err,error,*999)

    rowColCoupling%numberOfRowCols=0

    EXITS("MatrixRowColCoupling_Initialise")
    RETURN
999 ERRORSEXITS("MatrixRowColCoupling_Initialise",err,error)
    RETURN 1

  END SUBROUTINE MatrixRowColCoupling_Initialise
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the size of a matrix.
  SUBROUTINE Matrix_SizeSet(matrix,m,n,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: m !<The number of rows to set
    INTEGER(INTG), INTENT(IN) :: n !<The number of columns to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Matrix_SizeSet",err,error,*999)

    CALL Matrix_AssertNotFinished(matrix,err,error,*999)    
    IF(m<=0) THEN
      localError="The number of matrix rows of "//TRIM(NumberToVString(m,"*",err,error))// &
        & " is invalid. The number must be >0."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(n<=0) THEN
      localError="The number of matrix columns of "//TRIM(NumberToVString(n,"*",err,error))// &
        & " is invalid. The number must be >0."
      CALL FlagError(localError,err,error,*999)
    ENDIF
      
    matrix%m=m
    matrix%n=n

    IF(matrix%storageType==MATRIX_BLOCK_COMPRESSED_ROW_STORAGE_TYPE) THEN
      matrix%numberOfRowBlocks=CEILING(REAL(matrix%m,DP)/REAL(matrix%blockSize,DP),INTG)
      matrix%numberOfColumnBlocks=CEILING(REAL(matrix%n,DP)/REAL(matrix%blockSize,DP),INTG)
    ENDIF

    EXITS("Matrix_SizeSet")
    RETURN
999 ERRORSEXITS("Matrix_SizeSet",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_SizeSet

  !
  !================================================================================================================================
  !

  !>Returns the storage location in the data array of a matrix that correponds to location i,j. If the location does not exist the routine returns zero.
  SUBROUTINE Matrix_StorageLocationFind(matrix,i,j,location,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: i !<The row number of the location to find
    INTEGER(INTG), INTENT(IN) :: j !<The column number of the location to find
    INTEGER(INTG), INTENT(OUT) :: location !<On return the location of the specified row and column in the matrix data. If the row and column does not exist in the matrix then zero is returned.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: columnBlock,columnBlockNumber,columnOffset,k,lowLimit,midPoint,rowBlockNumber,rowOffset,upLimit
    LOGICAL :: foundColumn, foundRow
    TYPE(VARYING_STRING) :: localError

    ENTERS("Matrix_StorageLocationFind",err,error,*999)

    location=0
    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    IF(i<1.OR.i>matrix%m) THEN
      localError="Row number "//TRIM(NumberToVString(i,"*",err,error))//" is outside the matrix range of 1 to "// &
        & TRIM(NumberToVString(matrix%m,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(j<1.OR.j>matrix%n) THEN
      localError="Column number "//TRIM(NumberToVString(j,"*",err,error))//" is outside the matrix range of 1 to "// &
        & TRIM(NumberToVString(matrix%n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
      
    SELECT CASE(matrix%storageType)
    CASE(MATRIX_BLOCK_STORAGE_TYPE)
      location=i+(j-1)*matrix%m
    CASE(MATRIX_DIAGONAL_STORAGE_TYPE)
      IF(i==j) location=i
    CASE(MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
      location=i+(j-1)*matrix%maxM
    CASE(MATRIX_ROW_MAJOR_STORAGE_TYPE)
      location=(i-1)*matrix%maxN+j
    CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
      !Search for the column number in the sparsity list using the bisection (binary search) algorithm
      lowLimit=matrix%rowIndices(i)
      IF(lowLimit<=matrix%numberOfNonZeros) THEN
        IF(j>=matrix%columnIndices(lowLimit)) THEN
          upLimit=matrix%rowIndices(i+1)
          IF(upLimit>matrix%numberOfNonZeros) upLimit=upLimit-1
          IF(upLimit>lowLimit) THEN
            IF(j<=matrix%columnIndices(upLimit-1)) THEN
              DO WHILE((upLimit-lowLimit)>bisectionToLinearSearchThreshold)
                midPoint=(upLimit+lowLimit)/2
                IF(matrix%columnIndices(midPoint)>j) THEN
                  upLimit=midPoint
                ELSE
                  lowLimit=midPoint
                ENDIF
              ENDDO
              DO k=lowLimit,upLimit
                IF(matrix%columnIndices(k)==j) THEN
                  location=k
                  EXIT
                ENDIF
              ENDDO !k
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    CASE(MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
      !Search for the row number in the sparsity list using the bisection (binary search) algorithm
      lowLimit=matrix%columnIndices(j)
      IF(lowLimit<=matrix%numberOfNonZeros) THEN
        IF(i>=matrix%rowIndices(lowLimit)) THEN
          upLimit=matrix%columnIndices(j+1)
          IF(upLimit>matrix%numberOfNonZeros) upLimit=upLimit-1
          IF(upLimit>lowLimit) THEN
            IF(i<=matrix%rowIndices(upLimit-1)) THEN
              DO WHILE((upLimit-lowLimit)>bisectionToLinearSearchThreshold)
                midPoint=(upLimit+lowLimit)/2
                IF(matrix%rowIndices(midPoint)>i) THEN
                  upLimit=midPoint
                ELSE
                  lowLimit=midPoint
                ENDIF
              ENDDO
              DO k=lowLimit,upLimit
                IF(matrix%rowIndices(k)==i) THEN
                  location=k
                  EXIT
                ENDIF
              ENDDO !k
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    CASE(MATRIX_ROW_COLUMN_STORAGE_TYPE)
      foundRow=.FALSE.
      location=1
      DO WHILE(.NOT.foundColumn.AND.location<=matrix%size)
        IF(matrix%rowIndices(location)==i) THEN
          DO WHILE(.NOT.foundColumn.AND.location<=matrix%size)
            IF(matrix%columnIndices(location)==j.AND.matrix%rowIndices(location)==i) THEN
              foundColumn=.TRUE.
            ELSE IF(matrix%rowIndices(location)/=i) THEN
              location=matrix%size+1
            ELSE
              location=location+1
            ENDIF
          ENDDO
        ELSE
          location=location+1
        ENDIF
      ENDDO
      IF(.NOT.(foundRow.AND.foundColumn)) location=0
    CASE(MATRIX_BLOCK_COMPRESSED_ROW_STORAGE_TYPE)
      rowBlockNumber=FLOOR(REAL(i,DP)/REAL(matrix%blockSize,DP),INTG)
      rowOffset=MOD(i,matrix%blockSize)
      columnBlockNumber=FLOOR(REAL(j,DP)/REAL(matrix%blockSize,DP),INTG)
      columnOffset=MOD(j,matrix%blockSize)
      !Search for the column block number in the sparsity list using the bisection (binary search) algorithm
      columnBlock=0
      lowLimit=matrix%rowIndices(rowBlockNumber)
      IF(columnBlockNumber>=matrix%columnIndices(lowLimit)) THEN
        upLimit=matrix%rowIndices(rowBlockNumber+1)
        IF(upLimit>lowLimit) THEN
          IF(columnBlockNumber<=matrix%columnIndices(upLimit-1)) THEN
            DO WHILE((upLimit-lowLimit)>bisectionToLinearSearchThreshold)
              midPoint=(upLimit+lowLimit)/2
              IF(matrix%columnIndices(midPoint)>columnBlockNumber) THEN
                upLimit=midPoint
              ELSE
                lowLimit=midPoint
              ENDIF
            ENDDO
            DO k=lowLimit,upLimit
              IF(matrix%columnIndices(k)==columnBlockNumber) THEN
                columnBlock=k
                EXIT
              ENDIF
            ENDDO !k
          ENDIF
        ENDIF
      ENDIF
      IF(columnBlock/=0) location=columnBlock*matrix%blockSize*matrix%blockSize+(columnOffset-1)*matrix%blockSize+rowOffset
    CASE(MATRIX_MODIFIED_KRM_STORAGE_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The matrix storage type of "//TRIM(NumberToVString(matrix%storageType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Matrix_StorageLocationFind")
    RETURN
999 ERRORSEXITS("Matrix_StorageLocationFind",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_StorageLocationFind

  !
  !================================================================================================================================
  !

  !>Returns the storage location in a matrix corresponding to the specified row and column. If the storage location corresponding to the row and column does not exist then an error is flagged.
  SUBROUTINE Matrix_StorageLocationGet(matrix,row,column,location,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix to get the storage location for
    INTEGER(INTG), INTENT(IN) :: row !<The row to get
    INTEGER(INTG), INTENT(IN) :: column !<The column to get
    INTEGER(INTG), INTENT(OUT) :: location !<On return, the storage location corresponding to the row and column.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Matrix_StorageLocationGet",err,error,*999)


    CALL Matrix_StorageLocationFind(matrix,row,column,location,err,error,*999)
    IF(location==0) THEN
      localError="Row "//TRIM(NumberToVString(row,"*",err,error))//" and column "// &
        & TRIM(NumberToVString(column,"*",err,error))//" does not exist in the matrix."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    EXITS("Matrix_StorageLocationGet")
    RETURN
999 ERRORSEXITS("Matrix_StorageLocationGet",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_StorageLocationGet

  !
  !================================================================================================================================
  !

  !>Sets the storage locations (sparsity pattern) in a matrix to that specified by the row and column indices.
  SUBROUTINE Matrix_StorageLocationsSet(matrix,rowIndices,columnIndices,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The row index values for the sparisty pattern.
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(i). The column index values for the sparsity pattern.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,j,k
    INTEGER(INTG), ALLOCATABLE :: rowCounts(:)
    TYPE(VARYING_STRING) :: localError

    ENTERS("Matrix_StorageLocationsSet",err,error,*999)

    CALL Matrix_AssertNotFinished(matrix,err,error,*999)
    
    SELECT CASE(matrix%storageType)
    CASE(MATRIX_BLOCK_STORAGE_TYPE)
      CALL FlagError("Can not set matrix locations for a block storage matrix.",err,error,*999)
    CASE(MATRIX_DIAGONAL_STORAGE_TYPE)
      CALL FlagError("Can not set matrix locations for a diagonal storage matrix.",err,error,*999)
    CASE(MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
      CALL FlagError("Can not set matrix locations for a column major storage matrix.",err,error,*999)
    CASE(MATRIX_ROW_MAJOR_STORAGE_TYPE)
      CALL FlagError("Can not set matrix locations for a row major storage matrix.",err,error,*999)
    CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
      IF(SIZE(rowIndices,1)/=matrix%m+1) THEN
        localError="The supplied number of row indices of "//TRIM(NumberToVString(SIZE(rowIndices,1),"*",err,error))// &
          & " does not match the number of rows in the matrix + 1 of "// &
          & TRIM(NumberToVString(matrix%m+1,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF        
      IF(SIZE(columnIndices,1)/=matrix%numberOfNonZeros) THEN
        localError="The supplied number of column indices of "//TRIM(NumberToVString(SIZE(columnIndices,1),"*",err,error))// &
          & " does not match the number of non-zeros in the matrix of "// &
          & TRIM(NumberToVString(matrix%numberOfNonZeros,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF        
      IF(rowIndices(1)/=1) THEN
        localError="Invalid row indices. The first row index of "//TRIM(NumberToVString(rowIndices(1),"*",err,error))// &
          & " does not equal 1."
        CALL FlagError(localError,err,error,*999)
      ENDIF        
      IF(rowIndices(matrix%m+1)/=matrix%numberOfNonZeros+1) THEN
        localError="Invalid row indices. The last row index of "//TRIM(NumberToVString(rowIndices(matrix%m+1),"*",err,error))// &
          & " does not equal the number of non-zeros + 1 of "// &
          & TRIM(NumberToVString(matrix%numberOfNonZeros+1,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF        
      DO i=2,matrix%m+1
        IF(rowIndices(i)<rowIndices(i-1)) THEN
          localError="Invalid row indices. Row "//TRIM(NumberToVString(i,"*",err,error))//" index number "// &
            & TRIM(NumberToVString(rowIndices(i),"*",err,error))//" is less than row "// &
            & TRIM(NumberToVString(i-1,"*",err,error))//" index number "// &
            & TRIM(NumberToVString(rowIndices(i-1),"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDDO !i
      DO i=1,matrix%m
        DO j=rowIndices(i),rowIndices(i+1)-1
          k=columnIndices(j)
          IF(k>0) THEN
            IF(k>matrix%n) THEN
              localError="Invalid column indices. Column index "//TRIM(NumberToVString(j,"*",err,error))//" with value "// &
                & TRIM(NumberToVString(k,"*",err,error))//" is greater than the number of columns of "// &
                & TRIM(NumberToVString(matrix%n,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ELSE
            localError="Invalid column indices. Column index "//TRIM(NumberToVString(j,"*",err,error))//" with value "// &
              & TRIM(NumberToVString(k,"*",err,error))//" is less than zero."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ENDDO !j
      ENDDO !i
      
      IF(ALLOCATED(matrix%rowIndices)) DEALLOCATE(matrix%rowIndices)
      IF(ALLOCATED(matrix%columnIndices)) DEALLOCATE(matrix%columnIndices)
      ALLOCATE(matrix%rowIndices(matrix%m+1),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate matrix row indices.",err,error,*999)
      ALLOCATE(matrix%columnIndices(matrix%numberOfNonZeros),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate matrix column indices.",err,error,*999)                  
      matrix%rowIndices(1:matrix%m+1)=rowIndices(1:matrix%m+1)
      matrix%columnIndices(1:matrix%numberOfNonZeros)=columnIndices(1:matrix%numberOfNonZeros)
      !Don't really need this???
      !DO i=1,matrix%m
      !  CALL List_Sort(matrix%columnIndices(matrix%rowIndices(i):matrix%rowIndices(i+1)-1),err,error,*999)
      !ENDDO !i
    CASE(MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
      IF(SIZE(columnIndices,1)/=matrix%n+1) THEN
        localError="The supplied number of column indices of "// &
          & TRIM(NumberToVString(SIZE(columnIndices,1),"*",err,error))// &
          & " does not match the number of columns in the matrix + 1 of "// &
          & TRIM(NumberToVString(matrix%n+1,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF        
      IF(SIZE(rowIndices,1)==matrix%numberOfNonZeros) THEN
        localError="The supplied number of row indices of "// &
          & TRIM(NumberToVString(SIZE(rowIndices,1),"*",err,error))// &
          & " does not match the number of non-zeros in the matrix of "// &
          & TRIM(NumberToVString(matrix%numberOfNonZeros,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF        
      IF(columnIndices(1)/=1) THEN
        localError="Invalid column indices. The first column index of "// &
          & TRIM(NumberToVString(columnIndices(1),"*",err,error))//" does not equal 1."
        CALL FlagError(localError,err,error,*999)
      ENDIF        
      IF(columnIndices(matrix%n+1)==matrix%numberOfNonZeros+1) THEN
        localError="Invalid column indices. The last column index of "// &
          & TRIM(NumberToVString(columnIndices(matrix%n+1),"*",err,error))// &
          & " does not equal the number of non-zeros + 1 of "// &
          & TRIM(NumberToVString(matrix%numberOfNonZeros+1,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      DO j=2,matrix%n+1
        IF(columnIndices(j)<columnIndices(j-1)) THEN
          localError="Invalid column indices. Column "//TRIM(NumberToVString(j,"*",err,error))// &
            & " index number "//TRIM(NumberToVString(columnIndices(j),"*",err,error))//" is less than column "// &
            & TRIM(NumberToVString(j-1,"*",err,error))//" index number "// &
            & TRIM(NumberToVString(columnIndices(j-1),"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        END IF
        IF(columnIndices(j)<0.OR.columnIndices(j)>matrix%numberOfNonZeros+1) THEN
          localError="Invalid column indices. Column index "//TRIM(NumberToVString(j,"*",err,error))//" with value "// &
            & TRIM(NumberToVString(columnIndices(j),"*",err,error))//" "// &
            & " should be in the range of one to the number of non-zeros + 1 of "// &
            & TRIM(NumberToVString(matrix%numberOfNonZeros+1,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        END IF
      ENDDO !i
      DO j=1,matrix%n
        DO i=columnIndices(j),columnIndices(j+1)-1
          k=rowIndices(i)
          IF(k>0) THEN
            IF(k>matrix%m) THEN
              localError="Invalid row indices. Row index "//TRIM(NumberToVString(i,"*",err,error))//" with value "// &
                & TRIM(NumberToVString(k,"*",err,error))//" is greater than the number of rows of "// &
                & TRIM(NumberToVString(matrix%m,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ELSE
            localError="Invalid row indices. Row index "//TRIM(NumberToVString(i,"*",err,error))//" with value "// &
              & TRIM(NumberToVString(k,"*",err,error))//" is less than zero."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ENDDO !i
      ENDDO !j
      IF(ALLOCATED(matrix%rowIndices)) DEALLOCATE(matrix%rowIndices)
      IF(ALLOCATED(matrix%columnIndices)) DEALLOCATE(matrix%columnIndices)
      ALLOCATE(matrix%rowIndices(matrix%numberOfNonZeros),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate matrix row indices.",err,error,*999)
      ALLOCATE(matrix%columnIndices(matrix%n+1),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate matrix column indices.",err,error,*999)                  
      matrix%rowIndices(1:matrix%numberOfNonZeros)=rowIndices(1:matrix%numberOfNonZeros)
      matrix%columnIndices(1:matrix%n+1)=columnIndices(1:matrix%n+1)
      !Don't really need this???
      !DO j=1,matrix%n                    
      !  CALL List_Sort(matrix%rowIndices(matrix%columnIndices(j):matrix%columnIndices(j+1)-1),err,error,*999)
      !ENDDO !j
    CASE(MATRIX_ROW_COLUMN_STORAGE_TYPE)
      IF(SIZE(rowIndices,1)/=matrix%numberOfNonZeros) THEN
        localError="The supplied number of row indices of "// &
          & TRIM(NumberToVString(SIZE(rowIndices,1),"*",err,error))// &
          & " does not match the number of non-zeros in the matrix of "// &
          & TRIM(NumberToVString(matrix%numberOfNonZeros,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF        
      IF(SIZE(columnIndices,1)/=matrix%numberOfNonZeros) THEN
        localError="The supplied number of column indices of "// &
          & TRIM(NumberToVString(SIZE(columnIndices,1),"*",err,error))// &
          & " does not match the number of non-zeros in the matrix of "// &
          & TRIM(NumberToVString(matrix%numberOfNonZeros,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
        
      DO k=1,matrix%numberOfNonZeros
        IF(rowIndices(k)<1.OR.rowIndices(k)>matrix%m) THEN
          localError="Invalid row indices. Row index number "//TRIM(NumberToVString(k,"*",err,error))//" with value "// &
            & TRIM(NumberToVString(rowIndices(k),"*",err,error))// &
            & " is out of range. The row index must be between 1 and "// &
            & TRIM(NumberToVString(matrix%m,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ELSE IF(columnIndices(k)<1.OR.columnIndices(k)>matrix%n) THEN
          localError="Invalid column indices. Column index number "//TRIM(NumberToVString(k,"*",err,error))//" with value "// &
            & TRIM(NumberToVString(columnIndices(k),"*",err,error))// &
            & " is out of range. The column index must be between 1 and "// &
            & TRIM(NumberToVString(matrix%n,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDDO !k
      matrix%rowIndices(1:matrix%numberOfNonZeros)=rowIndices(1:matrix%numberOfNonZeros)
      matrix%columnIndices(1:matrix%numberOfNonZeros)=columnIndices(1:matrix%numberOfNonZeros)
!!TODO: sort the row and colum indices!!!!!
    CASE(MATRIX_BLOCK_COMPRESSED_ROW_STORAGE_TYPE)
      IF(SIZE(rowIndices,1)/=matrix%numberOfRowBlocks+1) THEN
        localError="The supplied number of row indices of "//TRIM(NumberToVString(SIZE(rowIndices,1),"*",err,error))// &
          & " does not match the number of row blocks in the matrix + 1 of "// &
          & TRIM(NumberToVString(matrix%numberOfRowBlocks+1,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF        
      IF(SIZE(columnIndices,1)>matrix%numberOfRowBlocks*matrix%numberOfColumnBlocks) THEN
        localError="The supplied number of column indices of "//TRIM(NumberToVString(SIZE(columnIndices,1),"*",err,error))// &
          & " is greater than the number of blocks in the matrix of "// &
          & TRIM(NumberToVString(matrix%numberOfRowBlocks*matrix%numberOfColumnBlocks,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      IF(rowIndices(1)/=1) THEN
        localError="Invalid row indices. The first row index of "//TRIM(NumberToVString(rowIndices(1),"*",err,error))// &
          & " does not equal 1."
        CALL FlagError(localError,err,error,*999)
      ENDIF        
      IF(rowIndices(matrix%numberOfRowBlocks+1)/=SIZE(columnIndices,1)) THEN
        localError="Invalid row indices. The last row index of "// &
          & TRIM(NumberToVString(rowIndices(matrix%numberOfRowBlocks+1),"*",err,error))// &
          & " does not equal the number of column blocks of "// &
          & TRIM(NumberToVString(SIZE(columnIndices,1),"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF        
      DO i=2,matrix%numberOfRowBlocks+1
        IF(rowIndices(i)<rowIndices(i-1)) THEN
          localError="Invalid row indices. Row "//TRIM(NumberToVString(i,"*",err,error))//" index number "// &
            & TRIM(NumberToVString(rowIndices(i),"*",err,error))//" is less than row "// &
            & TRIM(NumberToVString(i-1,"*",err,error))//" index number "// &
            & TRIM(NumberToVString(rowIndices(i-1),"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDDO !i
      DO i=1,matrix%m
        DO j=rowIndices(i),rowIndices(i+1)-1
          k=columnIndices(j)
          IF(k>0) THEN
            IF(k>matrix%numberOfColumnBlocks) THEN
              localError="Invalid column indices. Column index "//TRIM(NumberToVString(j,"*",err,error))//" with value "// &
                & TRIM(NumberToVString(k,"*",err,error))//" is greater than the number of columns blocks of "// &
                & TRIM(NumberToVString(matrix%numberOfColumnBlocks,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ELSE
            localError="Invalid column indices. Column index "//TRIM(NumberToVString(j,"*",err,error))//" with value "// &
              & TRIM(NumberToVString(k,"*",err,error))//" is less than zero."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ENDDO !j
      ENDDO !i
      matrix%numberOfBlocks=SIZE(columnIndices,1)
      IF(ALLOCATED(matrix%rowIndices)) DEALLOCATE(matrix%rowIndices)
      IF(ALLOCATED(matrix%columnIndices)) DEALLOCATE(matrix%columnIndices)
      ALLOCATE(matrix%rowIndices(matrix%numberOfRowBlocks+1),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate matrix row indices.",err,error,*999)
      ALLOCATE(matrix%columnIndices(matrix%numberOfBlocks),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate matrix column indices.",err,error,*999)                  
      matrix%rowIndices(1:matrix%numberOfRowBlocks+1)=rowIndices(1:matrix%numberOfRowBlocks+1)
      matrix%columnIndices(1:matrix%numberOfBlocks)=columnIndices(1:matrix%numberOfBlocks)
    CASE(MATRIX_MODIFIED_KRM_STORAGE_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The matrix storage type of "//TRIM(NumberToVString(matrix%storageType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("Matrix_StorageLocationsSet")
    RETURN
999 IF(ALLOCATED(rowCounts)) DEALLOCATE(rowCounts)
    ERRORSEXITS("Matrix_StorageLocationsSet",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_StorageLocationsSet

  !
  !================================================================================================================================
  !

  !>Sets/changes the storage type for a matrix.
  SUBROUTINE Matrix_StorageTypeSet(matrix,storageType,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: storageType !<The storage type to set. \see MatrixVector_StorageTypes,MatrixVector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Matrix_StorageTypeSet",err,error,*999)

    CALL Matrix_AssertNotFinished(matrix,err,error,*999)
    
    SELECT CASE(storageType)
    CASE(MATRIX_BLOCK_STORAGE_TYPE)
      matrix%storageType=MATRIX_BLOCK_STORAGE_TYPE
    CASE(MATRIX_DIAGONAL_STORAGE_TYPE)
      matrix%storageType=MATRIX_DIAGONAL_STORAGE_TYPE
    CASE(MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
      matrix%storageType=MATRIX_COLUMN_MAJOR_STORAGE_TYPE
    CASE(MATRIX_ROW_MAJOR_STORAGE_TYPE)
      matrix%storageType=MATRIX_ROW_MAJOR_STORAGE_TYPE
    CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
      matrix%storageType=MATRIX_COMPRESSED_ROW_STORAGE_TYPE
    CASE(MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
      matrix%storageType=MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE
    CASE(MATRIX_ROW_COLUMN_STORAGE_TYPE)
      matrix%storageType=MATRIX_ROW_COLUMN_STORAGE_TYPE
    CASE(MATRIX_BLOCK_COMPRESSED_ROW_STORAGE_TYPE)
      matrix%storageType=MATRIX_BLOCK_COMPRESSED_ROW_STORAGE_TYPE
    CASE(MATRIX_MODIFIED_KRM_STORAGE_TYPE)
      matrix%storageType=MATRIX_MODIFIED_KRM_STORAGE_TYPE
    CASE DEFAULT
      localError="The matrix storage type of "//TRIM(NumberToVString(storageType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Matrix_StorageTypeSet")
    RETURN
999 ERRORSEXITS("Matrix_StorageTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_StorageTypeSet

  !
  !================================================================================================================================
  !

  !>Sets/changes the symmetry type for a matrix.
  SUBROUTINE Matrix_SymmetryTypeSet(matrix,symmetryType,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix to set the symmetry type for.
    INTEGER(INTG), INTENT(IN) :: symmetryType !<The symmetry type to set. \see MatrixVector_SymmetryTypes,MatrixVector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Matrix_SymmetryTypeSet",err,error,*999)

    CALL Matrix_AssertNotFinished(matrix,err,error,*999)
   
    SELECT CASE(symmetryType)
    CASE(MATRIX_SYMMETRIC_TYPE)
      matrix%symmetryType=MATRIX_SYMMETRIC_TYPE
    CASE(MATRIX_HERMITIAN_TYPE)
      IF(matrix%dataType==MATRIX_VECTOR_SPC_TYPE.OR.matrix%dataType==MATRIX_VECTOR_DPC_TYPE) THEN
        matrix%symmetryType=MATRIX_HERMITIAN_TYPE
      ELSE
        CALL FlagError("Cannot set the matrix symmetry type to Hermitian as the matrix does not have a complex data type.", &
          & err,error,*999)
      ENDIF
    CASE(MATRIX_SKEW_SYMMETRIC_TYPE)
      matrix%symmetryType=MATRIX_SKEW_SYMMETRIC_TYPE
    CASE(MATRIX_UNSYMMETRIC_TYPE)
      matrix%symmetryType=MATRIX_UNSYMMETRIC_TYPE
    CASE DEFAULT
      localError="The matrix symmetry type of "//TRIM(NumberToVString(symmetryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("Matrix_SymmetryTypeSet")
    RETURN
999 ERRORSEXITS("Matrix_SymmetryTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_SymmetryTypeSet

  !
  !================================================================================================================================
  !

  !>Calculates the transpose rows/columns for a matrix.
  SUBROUTINE Matrix_TransposeLocationsCalculate(matrix,err,error,*,transposeRowsColumns)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix to set the transpose rows/columns for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    INTEGER(INTG), INTENT(IN), OPTIONAL :: transposeRowsColumns(:) !<transposeRowsColumns(i). The list of transpose rows/columns to calculate. If this parameter is not present then all rows/columns are included in the transpose
    !Local Variables
    INTEGER(INTG) :: columnIdx,dummyErr,listItem(2),location,nonZeroIdx,numberOfNonZeros,numberOfRowsColumns, &
      & numberOfRows,rowColumnIdx,rowIdx
    INTEGER(INTG), ALLOCATABLE :: newColumnIndicesT(:),newRowIndicesT(:),newTransposeDataSwivel(:),newTransposeRowsColumns(:)
    INTEGER(INTG), ALLOCATABLE :: rowColumnData(:,:)
    LOGICAL :: fullTranspose
    TYPE(ListPtrType), ALLOCATABLE :: rowColumnLists(:)
    TYPE(VARYING_STRING) :: dummyError,localError

    ENTERS("Matrix_TransposeLocationsCalculate",err,error,*999)
 
    SELECT CASE(matrix%storageType)
    CASE(MATRIX_BLOCK_STORAGE_TYPE)
      !Do nothing, transpose is trivial to obtain
    CASE(MATRIX_DIAGONAL_STORAGE_TYPE)
      !Do nothing, transpose is trivial to obtain
    CASE(MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
      !Do nothing, transpose is trivial to obtain
    CASE(MATRIX_ROW_MAJOR_STORAGE_TYPE)
      !Do nothing, transpose is trivial to obtain
    CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
      !Store the transpose as a compressed column scheme      
      IF(PRESENT(transposeRowsColumns)) THEN
        numberOfRowsColumns=SIZE(transposeRowsColumns,1)
        ALLOCATE(newTransposeRowsColumns(numberOfRowsColumns),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate new transpose rows columns array.",err,error,*999)
        newTransposeRowsColumns(1:numberOfRowsColumns)=transposeRowsColumns(1:numberOfRowsColumns)
        CALL List_SortHeap(newTransposeRowsColumns,err,error,*999)
        fullTranspose=.FALSE.
      ELSE
        numberOfRowsColumns=matrix%n
        fullTranspose=.TRUE.
      ENDIF
      !Create lists of the column entries for the specified columns
      ALLOCATE(rowColumnLists(numberOfRowsColumns),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate row column lists.",err,error,*999)
      DO rowColumnIdx=1,numberOfRowsColumns
        NULLIFY(rowColumnLists(rowColumnIdx)%ptr)
      ENDDO !rowColumnIdx
      numberOfNonZeros=0
      DO rowColumnIdx=1,numberOfRowsColumns
        IF(fullTranspose) THEN
          columnIdx=rowColumnIdx
        ELSE
          columnIdx=newTransposeRowsColumns(rowColumnIdx)
        ENDIF
        CALL List_CreateStart(rowColumnLists(rowColumnIdx)%ptr,err,error,*999)
        CALL List_DataDimensionSet(rowColumnLists(rowColumnIdx)%ptr,2,err,error,*999)
        CALL List_CreateFinish(rowColumnLists(rowColumnIdx)%ptr,err,error,*999)
        DO rowIdx=1,matrix%m
          CALL Matrix_StorageLocationFind(matrix,rowIdx,columnIdx,location,err,error,*999)
          IF(location/=0) THEN
            numberOfNonZeros=numberOfNonZeros+1
            listItem(1)=rowIdx
            listItem(2)=location
            CALL List_ItemAdd(rowColumnLists(rowColumnIdx)%ptr,listItem,err,error,*999)
          ENDIF
        ENDDO !rowIdx
      ENDDO !rowColumnIdx
      ALLOCATE(newColumnIndicesT(numberOfRowsColumns+1),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate new column indices transpose.",err,error,*999)
      ALLOCATE(newRowIndicesT(numberOfNonZeros),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate new row indices transpose.",err,error,*999)
      ALLOCATE(newTransposeDataSwivel(numberOfNonZeros),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate new transpose data swivel.",err,error,*999)
      newColumnIndicesT=0
      newColumnIndicesT(1)=1
      numberOfNonZeros=0
      DO rowColumnIdx=1,numberOfRowsColumns
        IF(fullTranspose) THEN
          columnIdx=rowColumnIdx
        ELSE
          columnIdx=newTransposeRowsColumns(rowColumnIdx)
        ENDIF
        CALL List_DetachAndDestroy(rowColumnLists(rowColumnIdx)%ptr,numberOfRows,rowColumnData,err,error,*999)
        numberOfNonZeros=numberOfNonZeros+numberOfRows
        IF(numberOfRows>matrix%maximumRowIndicesPerColumn) matrix%maximumRowIndicesPerColumn=numberOfRows
        DO nonZeroIdx=1,numberOfRows
          rowIdx=rowColumnData(1,nonZeroIdx)
          location=rowColumnData(2,nonZeroIdx)
          newRowIndicesT(newColumnIndicesT(rowColumnIdx)+nonZeroIdx-1)=rowIdx
          newTransposeDataSwivel(newColumnIndicesT(rowColumnIdx)+nonZeroIdx-1)=location
        ENDDO !nonZeroIdx
        newColumnIndicesT(rowColumnIdx+1)=numberOfNonZeros+1
        IF(ALLOCATED(rowColumnData)) DEALLOCATE(rowColumnData)
      ENDDO !rowColumnIdx
      IF(.NOT.fullTranspose) CALL MOVE_ALLOC(newTransposeRowsColumns,matrix%transposeRowsColumns)
      CALL MOVE_ALLOC(newColumnIndicesT,matrix%columnIndicesT)
      CALL MOVE_ALLOC(newRowIndicesT,matrix%rowIndicesT)
      CALL MOVE_ALLOC(newTransposeDataSwivel,matrix%transposeDataSwivel)
    CASE(MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
      !Store the transpose as a compressed row scheme
      !Create lists of the row entries for the specified rows
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(MATRIX_ROW_COLUMN_STORAGE_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(MATRIX_MODIFIED_KRM_STORAGE_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The matrix storage type of "//TRIM(NumberToVString(matrix%storageType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
   
    EXITS("Matrix_TransposeLocationsCalculate")
    RETURN
999 IF(ALLOCATED(newColumnIndicesT)) DEALLOCATE(newColumnIndicesT)
    IF(ALLOCATED(newRowIndicesT)) DEALLOCATE(newRowIndicesT)
    IF(ALLOCATED(newTransposeDataSwivel)) DEALLOCATE(newTransposeDataSwivel)
    IF(ALLOCATED(rowColumnData)) DEALLOCATE(rowColumnData)
    IF(ALLOCATED(rowColumnLists)) THEN
      DO rowColumnIdx=1,SIZE(rowColumnLists,1)
        IF(ASSOCIATED(rowColumnLists(rowColumnIdx)%ptr)) &
          & CALL List_Destroy(rowColumnLists(rowColumnIdx)%ptr,dummyErr,dummyError,*998)
      ENDDO !rowColumnIdx
    ENDIF
998 ERRORSEXITS("Matrix_TranposeLocationsCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_TransposeLocationsCalculate

  !
  !================================================================================================================================
  !

  !>Returns the position in the transpose indices list for the specified row/column number.
  SUBROUTINE Matrix_TransposeRowColumnPositionGet(matrix,rowColumnNumber,transposePosition,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix to find the transpose position for.
    INTEGER(INTG), INTENT(IN) :: rowColumnNumber !<The row/column number to find the position for
    INTEGER(INTG), INTENT(OUT) :: transposePosition !<On return, the position in the transfer indices for the row/column. If the row/column is not in the transfer indices then the transpose position will be zero. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Matrix_TransposeRowColumnPositionGet",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
   
    SELECT CASE(matrix%transposeType)
    CASE(MATRIX_NO_TRANSPOSE_REQUIRED)
      transposePosition=0
    CASE(MATRIX_PARTIAL_TRANSPOSE_REQUIRED)
      SELECT CASE(matrix%storageType)
      CASE(MATRIX_BLOCK_STORAGE_TYPE)
        transposePosition=0
      CASE(MATRIX_DIAGONAL_STORAGE_TYPE)
        transposePosition=0
      CASE(MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
        transposePosition=0
      CASE(MATRIX_ROW_MAJOR_STORAGE_TYPE)
        transposePosition=0
      CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
        IF(rowColumnNumber<=0.OR.rowColumnNumber>matrix%n) THEN
          localError="The specified column number of "//TRIM(NumberToVString(rowColumnNumber,"*",err,error))// &
            & " is invalid. The column number should be >=1 and <= "//TRIM(NumberToVString(matrix%n,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        CALL List_SearchBinary(matrix%transposeRowsColumns,rowColumnNumber,transposePosition,err,error,*999)
      CASE(MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
        IF(rowColumnNumber<=0.OR.rowColumnNumber>matrix%M) THEN
          localError="The specified row number of "//TRIM(NumberToVString(rowColumnNumber,"*",err,error))// &
            & " is invalid. The row number should be >=1 and <= "//TRIM(NumberToVString(matrix%m,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        CALL List_SearchBinary(matrix%transposeRowsColumns,rowColumnNumber,transposePosition,err,error,*999)
      CASE(MATRIX_ROW_COLUMN_STORAGE_TYPE)
        transposePosition=0
      CASE(MATRIX_BLOCK_COMPRESSED_ROW_STORAGE_TYPE)
        transposePosition=0
      CASE(MATRIX_MODIFIED_KRM_STORAGE_TYPE)
        transposePosition=0
      CASE DEFAULT
        localError="The matrix storage type of "//TRIM(NumberToVString(matrix%storageType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(MATRIX_FULL_TRANSPOSE_REQUIRED)
      SELECT CASE(matrix%storageType)
      CASE(MATRIX_BLOCK_STORAGE_TYPE)
        transposePosition=0
      CASE(MATRIX_DIAGONAL_STORAGE_TYPE)
        transposePosition=0
      CASE(MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
        transposePosition=0
      CASE(MATRIX_ROW_MAJOR_STORAGE_TYPE)
        transposePosition=0
      CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
        IF(rowColumnNumber<=0.OR.rowColumnNumber>matrix%n) THEN
          localError="The specified column number of "//TRIM(NumberToVString(rowColumnNumber,"*",err,error))// &
            & " is invalid. The column number should be >=1 and <= "//TRIM(NumberToVString(matrix%n,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        transposePosition=rowColumnNumber
      CASE(MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
        IF(rowColumnNumber<=0.OR.rowColumnNumber>matrix%M) THEN
          localError="The specified row number of "//TRIM(NumberToVString(rowColumnNumber,"*",err,error))// &
            & " is invalid. The row number should be >=1 and <= "//TRIM(NumberToVString(matrix%m,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        transposePosition=rowColumnNumber
      CASE(MATRIX_ROW_COLUMN_STORAGE_TYPE)
        transposePosition=0
      CASE(MATRIX_BLOCK_COMPRESSED_ROW_STORAGE_TYPE)
        transposePosition=0
      CASE(MATRIX_MODIFIED_KRM_STORAGE_TYPE)
        transposePosition=0
      CASE DEFAULT
        localError="The matrix storage type of "//TRIM(NumberToVString(matrix%storageType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The matrix transpose type of "//TRIM(NumberToVString(matrix%transposeType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("Matrix_TransposeRowColumnPositionGet")
    RETURN
999 ERRORSEXITS("Matrix_TransposeRowColumnPositionGet",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_TransposeRowColumnPositionGet

  !
  !================================================================================================================================
  !

  !>Sets/changes the tranpose type for a matrix.
  SUBROUTINE Matrix_TransposeTypeSet(matrix,transposeType,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix to set the transpose type for.
    INTEGER(INTG), INTENT(IN) :: transposeType !<The transpose type to set. \see MatrixVector_TransposeTypes,MatrixVector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Matrix_TransposeTypeSet",err,error,*999)

    CALL Matrix_AssertNotFinished(matrix,err,error,*999)
   
    SELECT CASE(transposeType)
    CASE(MATRIX_NO_TRANSPOSE_REQUIRED)
      matrix%transposeType=MATRIX_NO_TRANSPOSE_REQUIRED
    CASE(MATRIX_PARTIAL_TRANSPOSE_REQUIRED)
      matrix%transposeType=MATRIX_PARTIAL_TRANSPOSE_REQUIRED
    CASE(MATRIX_FULL_TRANSPOSE_REQUIRED)
      matrix%transposeType=MATRIX_FULL_TRANSPOSE_REQUIRED
    CASE DEFAULT
      localError="The matrix transpose type of "//TRIM(NumberToVString(transposeType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("Matrix_TransposeTypeSet")
    RETURN
999 ERRORSEXITS("Matrix_TranposeTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_TransposeTypeSet

  !
  !================================================================================================================================
  !

  !>Sets/changes the transpose row/column for a matrix.
  SUBROUTINE Matrix_TransposeRowsColumnsSet0(matrix,transposeRowColumn,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix to set the transpose row/column for.
    INTEGER(INTG), INTENT(IN) :: transposeRowColumn !<The transpose row/column to set. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Matrix_TransposeRowsColumnsSet0",err,error,*999)

    CALL Matrix_TransposeRowsColumnsSet(matrix,[transposeRowColumn],err,error,*999)
    
    EXITS("Matrix_TransposeRowsColumnsSet0")
    RETURN
999 ERRORSEXITS("Matrix_TranposeRowsColumnsSet0",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_TransposeRowsColumnsSet0

  !
  !================================================================================================================================
  !

  !>Sets/changes the transpose rows/columns for a matrix.
  SUBROUTINE Matrix_TransposeRowsColumnsSet1(matrix,transposeRowsColumns,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix to set the transpose rows/columns for.
    INTEGER(INTG), INTENT(IN) :: transposeRowsColumns(:) !<transposeRowsColumns(i). The list of transpose rows/columns to set. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: columnIdx,numberOfRowsColumns
    INTEGER(INTG), ALLOCATABLE :: newTransposeRowsColumns(:)
    TYPE(VARYING_STRING) :: localError

    ENTERS("Matrix_TransposeRowsColumnsSet1",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    IF(SIZE(transposeRowsColumns,1)<1.OR.SIZE(transposeRowsColumns,1)>matrix%n) THEN
      localError="The size of the transpose rows columns array of "// &
        & TRIM(NumberToVString(SIZE(transposeRowsColumns,1),"*",err,error))//" is invalid. The size must be >= 1 and <= "// &
        & TRIM(NumberToVString(matrix%n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    DO columnIdx=1,SIZE(transposeRowsColumns,1)
      IF(transposeRowsColumns(columnIdx)<1.OR.transposeRowsColumns(columnIdx)>matrix%n) THEN
        localError="The column number of "//TRIM(NumberToVString(transposeRowsColumns(columnIdx),"*",err,error))// &
          & " at column index position "//TRIM(NumberToVString(columnIdx,"*",err,error))// &
          & " is invalid. The column number must be >= 1 and <= "//TRIM(NumberToVString(matrix%n,"*",err,error))// &
          & "."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDDO !columnIdx

    SELECT CASE(matrix%storageType)
    CASE(MATRIX_BLOCK_STORAGE_TYPE)
      !Do nothing, transpose is trivial to obtain
    CASE(MATRIX_DIAGONAL_STORAGE_TYPE)
      !Do nothing, transpose is trivial to obtain
    CASE(MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
      !Do nothing, transpose is trivial to obtain
    CASE(MATRIX_ROW_MAJOR_STORAGE_TYPE)
      !Do nothing, transpose is trivial to obtain
    CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
      !Store the transpose as a compressed column scheme
      numberOfRowsColumns=SIZE(transposeRowsColumns,1)
      ALLOCATE(newTransposeRowsColumns(numberOfRowsColumns),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate new transpose rows columns array.",err,error,*999)
      newTransposeRowsColumns(1:numberOfRowsColumns)=transposeRowsColumns(1:numberOfRowsColumns)
      CALL List_SortHeap(newTransposeRowsColumns,err,error,*999)
      IF(err/=0) CALL FlagError("Could not allocate new transpose data swivel.",err,error,*999)
      CALL MOVE_ALLOC(newTransposeRowsColumns,matrix%transposeRowsColumns)
      CALL Matrix_TransposeLocationsCalculate(matrix,err,error,*999,matrix%transposeRowsColumns)
    CASE(MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
      !Store the transpose as a compressed row scheme
      !Create lists of the row entries for the specified rows
    CASE(MATRIX_ROW_COLUMN_STORAGE_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(MATRIX_MODIFIED_KRM_STORAGE_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The matrix storage type of "//TRIM(NumberToVString(matrix%storageType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("Matrix_TransposeRowsColumnsSet1")
    RETURN
999 ERRORSEXITS("Matrix_TranposeRowsColumnsSet1",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_TransposeRowsColumnsSet1

  !
  !================================================================================================================================
  !

  !>Adds a value to an integer matrix at the location specified by the row and column index i.e., matrix(i,j)=matrix(i,j)+value
  SUBROUTINE Matrix_ValuesAddIntg(matrix,rowIndex,columnIndex,value,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: rowIndex !<The row index for the value to add
    INTEGER(INTG), INTENT(IN) :: columnIndex !<The column index for the value to add
    INTEGER(INTG), INTENT(IN) :: value !<The value to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: location

    ENTERS("Matrix_ValuesAddIntg",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    CALL Matrix_AssertIsIntgData(matrix,err,error,*999)

    CALL Matrix_StorageLocationGet(matrix,rowIndex,columnIndex,location,err,error,*999)
    matrix%dataIntg(location)=matrix%dataIntg(location)+value

    EXITS("Matrix_ValuesAddIntg")
    RETURN
999 ERRORSEXITS("Matrix_ValuesAddIntg",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_ValuesAddIntg

  !
  !================================================================================================================================
  !

  !>Adds values to an integer matrix at the location specified by the row and column indices i.e., matrix(i,j)=matrix(i,j)+value
  SUBROUTINE Matrix_ValuesAddIntg1(matrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The row index for the i'th value to add
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(i). The column index for the i'th value to add
    INTEGER(INTG), INTENT(IN) :: values(:) !<values(i). The value of the i'th value to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: k,location
    TYPE(VARYING_STRING) :: localError

    ENTERS("Matrix_ValuesAddIntg1",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    CALL Matrix_AssertIsIntgData(matrix,err,error,*999)
    IF(SIZE(rowIndices,1)/=SIZE(values,1)) THEN
      localError="The size of the row indices array of "//TRIM(NumberToVString(SIZE(rowIndices,1),"*",err,error))// &
        & " does not conform to the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF      
    IF(SIZE(columnIndices,1)/=SIZE(values,1)) THEN
      localError="The size of the column indices array of "//TRIM(NumberToVString(SIZE(columnIndices,1),"*",err,error))// &
        & " does not conform to the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    DO k=1,SIZE(rowIndices,1)
      CALL Matrix_StorageLocationGet(matrix,rowIndices(k),columnIndices(k),location,err,error,*999)
      matrix%dataIntg(location)=matrix%dataIntg(location)+values(k)
    ENDDO !k
    
    EXITS("Matrix_ValuesAddIntg1")
    RETURN
999 ERRORSEXITS("Matrix_ValuesAddIntg1",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_ValuesAddIntg1

  !
  !================================================================================================================================
  !

  !>Adds a matrix of values to an integer matrix at the location specified by the row and column indices i.e., matrix(i,j)=matrix(i,j)+value
  SUBROUTINE Matrix_ValuesAddIntg2(matrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The row index for the ij'th value to add
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(j). The column index for the ij'th value to add
    INTEGER(INTG), INTENT(IN) :: values(:,:) !<values(i,j). The value of the ij'th value to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: blockLocation,i,j,k,rowIndex,previousRowIndex,columnIndex,columnNumber,columnBlockNumber,columnOffset, &
      & previousColumnIndex,location,lowLimit,midPoint,previousColumnBlockNumber,rowBlockNumber,rowNumber,rowOffset,upLimit
    TYPE(VARYING_STRING) :: localError

    ENTERS("Matrix_ValuesAddIntg2",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    CALL Matrix_AssertIsIntgData(matrix,err,error,*999)    
    IF(SIZE(rowIndices,1)/=SIZE(values,1)) THEN
      localError="The size of the row indices array of "//TRIM(NumberToVString(SIZE(rowIndices,1),"*",err,error))// &
        & " does not conform to the number of rows the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(columnIndices,1)/=SIZE(values,2)) THEN
      localError="The size of the column indices array of "//TRIM(NumberToVString(SIZE(columnIndices,1),"*",err,error))// &
        & " does not conform to the number of columns in the values array of "// &
        & TRIM(NumberToVString(SIZE(values,2),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(matrix%storageType)
    CASE(MATRIX_BLOCK_STORAGE_TYPE,MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
      DO j=1,SIZE(columnIndices,1)
        DO i=1,SIZE(rowIndices,1)
          CALL Matrix_StorageLocationGet(matrix,rowIndices(i),columnIndices(j),location,err,error,*999)
          matrix%dataIntg(location)=matrix%dataIntg(location)+values(i,j)
        ENDDO !i
      ENDDO !j
    CASE(MATRIX_DIAGONAL_STORAGE_TYPE,MATRIX_ROW_MAJOR_STORAGE_TYPE,MATRIX_ROW_COLUMN_STORAGE_TYPE)
      DO i=1,SIZE(rowIndices,1)
        DO j=1,SIZE(columnIndices,1)
          CALL Matrix_StorageLocationGet(matrix,rowIndices(i),columnIndices(j),location,err,error,*999)
          matrix%dataIntg(location)=matrix%dataIntg(location)+values(i,j)
        ENDDO !j
      ENDDO !i
    CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
      !Search for the column number in the sparsity list using the bisection (binary search) algorithm in order
      !to take advantage of previous searches.
      previousColumnIndex=-1
      DO i=1,SIZE(rowIndices,1)
        rowIndex=rowIndices(i)
        lowLimit=matrix%rowIndices(rowIndex)
        upLimit=matrix%rowIndices(rowIndex+1)
        DO j=1,SIZE(columnIndices,1)
          location=0
          columnIndex=columnIndices(j)
          IF(columnIndex<=previousColumnIndex) THEN
            lowLimit=matrix%rowIndices(rowIndex)
          ELSE
            upLimit=matrix%rowIndices(rowIndex+1)
          ENDIF
          previousColumnIndex=columnIndex
          DO WHILE((upLimit-lowLimit)>bisectionToLinearSearchThreshold)
            midPoint=(upLimit+lowLimit)/2
            IF(matrix%columnIndices(midPoint)>columnIndex) THEN
              upLimit=midPoint
            ELSE
              lowLimit=midPoint
            ENDIF
          ENDDO
          DO k=lowLimit,upLimit
            IF(matrix%columnIndices(k)==columnIndex) THEN
              location=k
              lowLimit=k+1
              EXIT
            ENDIF
          ENDDO !k
          IF(location==0) THEN
            localError="Row "//TRIM(NumberToVString(rowIndices(i),"*",err,error))//" and column "// &
              & TRIM(NumberToVString(columnIndices(j),"*",err,error))//" does not exist in the matrix."
            CALL FlagError(localError,err,error,*999)
          ELSE
            matrix%dataIntg(location)=matrix%dataIntg(location)+values(i,j)
          ENDIF
        ENDDO !j
      ENDDO !i
    CASE(MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
      !Search for the row number in the sparsity list using the bisection (binary search) algorithm in order
      !to take advantage of previous searches.
      previousRowIndex=-1
      DO j=1,SIZE(columnIndices,1)
        columnIndex=columnIndices(j)
        lowLimit=matrix%columnIndices(columnIndex)
        upLimit=matrix%columnIndices(columnIndex+1)
        DO i=1,SIZE(rowIndices,1)
          location=0
          rowIndex=rowIndices(i)
          IF(rowIndex<=previousRowIndex) THEN
            lowLimit=matrix%columnIndices(columnIndex)
          ELSE
            upLimit=matrix%columnIndices(columnIndex+1)
          ENDIF
          previousRowIndex=rowIndex
          DO WHILE((upLimit-lowLimit)>bisectionToLinearSearchThreshold)
            midPoint=(upLimit+lowLimit)/2
            IF(matrix%rowIndices(midPoint)>rowIndex) THEN
              upLimit=midPoint
            ELSE
              lowLimit=midPoint
            ENDIF
          ENDDO
          DO k=lowLimit,upLimit
            IF(matrix%rowIndices(k)==rowIndex) THEN
              location=k
              lowLimit=k+1
              EXIT
            ENDIF
          ENDDO !k
          IF(location==0) THEN
            localError="Row "//TRIM(NumberToVString(rowIndices(i),"*",err,error))//" and column "// &
              & TRIM(NumberToVString(columnIndices(j),"*",err,error))//" does not exist in the matrix."
            CALL FlagError(localError,err,error,*999)
          ELSE
            matrix%dataIntg(location)=matrix%dataIntg(location)+values(i,j)
          ENDIF
        ENDDO !i
      ENDDO !j
    CASE(MATRIX_BLOCK_COMPRESSED_ROW_STORAGE_TYPE)
      !Search for the column number in the sparsity list using the bisection (binary search) algorithm in order
      !to take advantage of previous searches.
      previousColumnBlockNumber=-1
      DO i=1,SIZE(rowIndices,1)
        rowNumber=rowIndices(i)
        rowBlockNumber=FLOOR(REAL(rowNumber,DP)/REAL(matrix%blockSize,DP),INTG)
        rowOffset=MOD(rowNumber,matrix%blockSize)
        lowLimit=matrix%rowIndices(rowBlockNumber)
        upLimit=matrix%rowIndices(rowBlockNumber+1)
        DO j=1,SIZE(columnIndices,1)
          columnNumber=columnIndices(j)
          columnBlockNumber=FLOOR(REAL(columnNumber,DP)/REAL(matrix%blockSize,DP),INTG)
          columnOffset=MOD(columnNumber,matrix%blockSize)
          blockLocation=0
          IF(columnBlockNumber<=previousColumnBlockNumber) THEN
            lowLimit=matrix%rowIndices(rowBlockNumber)
          ELSE
            upLimit=matrix%rowIndices(rowBlockNumber+1)
          ENDIF
          previousColumnBlockNumber=columnBlockNumber
          DO WHILE((upLimit-lowLimit)>bisectionToLinearSearchThreshold)
            midPoint=(upLimit+lowLimit)/2
            IF(matrix%columnIndices(midPoint)>columnBlockNumber) THEN
              upLimit=midPoint
            ELSE
              lowLimit=midPoint
            ENDIF
          ENDDO
          DO k=lowLimit,upLimit
            IF(matrix%columnIndices(k)==columnBlockNumber) THEN
              blockLocation=k
              lowLimit=k+1
              EXIT
            ENDIF
          ENDDO !k
          IF(blockLocation==0) THEN
            localError="Row "//TRIM(NumberToVString(rowIndices(i),"*",err,error))//" and column "// &
              & TRIM(NumberToVString(columnIndices(j),"*",err,error))//" does not exist in the matrix."
            CALL FlagError(localError,err,error,*999)
          ELSE
            location=columnBlockNumber*matrix%blockSize*matrix%blockSize+(columnOffset-1)*matrix%blockSize+rowOffset
            matrix%dataIntg(location)=matrix%dataIntg(location)+values(i,j)
          ENDIF
        ENDDO !j
      ENDDO !i
    CASE(MATRIX_MODIFIED_KRM_STORAGE_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The matrix storage type of "//TRIM(NumberToVString(matrix%storageType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Matrix_ValuesAddIntg2")
    RETURN
999 ERRORSEXITS("Matrix_ValuesAddIntg2",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_ValuesAddIntg2

  !
  !================================================================================================================================
  !

  !>Adds a value to a single precision real matrix at the location specified by the row and column index i.e., matrix(i,j)=matrix(i,j)+value
  SUBROUTINE Matrix_ValuesAddSP(matrix,rowIndex,columnIndex,VALUE,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: rowIndex !<The row index for the value to add
    INTEGER(INTG), INTENT(IN) :: columnIndex !<The column index for the value to add
    REAL(SP), INTENT(IN) :: value !<The value to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: location

    ENTERS("Matrix_ValuesAddSP",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    CALL Matrix_AssertIsSPData(matrix,err,error,*999)

    CALL Matrix_StorageLocationGet(matrix,rowIndex,columnIndex,location,err,error,*999)
    matrix%dataSP(location)=matrix%dataSP(location)+value

    EXITS("Matrix_ValuesAddSP")
    RETURN
999 ERRORSEXITS("Matrix_ValuesAddSP",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_ValuesAddSP

  !
  !================================================================================================================================
  !

  !>Adds values to a single precision real matrix at the location specified by the row and column indices i.e., matrix(i,j)=matrix(i,j)+value
  SUBROUTINE Matrix_ValuesAddSP1(matrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The row index for the i'th value to add
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(i). The column index for the i'th value to add
    REAL(SP), INTENT(IN) :: values(:) !<values(i). The value of the i'th value to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: k,location
    TYPE(VARYING_STRING) :: localError

    ENTERS("Matrix_ValuesAddSP1",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    CALL Matrix_AssertIsSPData(matrix,err,error,*999)
    IF(SIZE(rowIndices,1)/=SIZE(values,1)) THEN
      localError="The size of the row indices array of "//TRIM(NumberToVString(SIZE(rowIndices,1),"*",err,error))// &
        & " does not conform to the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF      
    IF(SIZE(columnIndices,1)/=SIZE(values,1)) THEN
      localError="The size of the column indices array of "//TRIM(NumberToVString(SIZE(columnIndices,1),"*",err,error))// &
        & " does not conform to the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    DO k=1,SIZE(rowIndices,1)
      CALL Matrix_StorageLocationGet(matrix,rowIndices(k),columnIndices(k),location,err,error,*999)
      matrix%dataSP(location)=matrix%dataSP(location)+values(k)
    ENDDO !k
    
    EXITS("Matrix_ValuesAddSP1")
    RETURN
999 ERRORSEXITS("Matrix_ValuesAddSP1",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_ValuesAddSP1

  !
  !================================================================================================================================
  !

  !>Adds a matrix of values to a single precision real matrix at the location specified by the row and column indices i.e., matrix(i,j)=matrix(i,j)+value
  SUBROUTINE Matrix_ValuesAddSP2(matrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The row index for the ij'th value to add
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(j). The column index for the ij'th value to add
    REAL(SP), INTENT(IN) :: values(:,:) !<values(i,j). The value of the ij'th value to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: blockLocation,i,j,k,rowIndex,previousRowIndex,columnIndex,columnNumber,columnBlockNumber,columnOffset, &
      & previousColumnIndex,location,lowLimit,midPoint,previousColumnBlockNumber,rowBlockNumber,rowNumber,rowOffset,upLimit
    TYPE(VARYING_STRING) :: localError

    ENTERS("Matrix_ValuesAddSP2",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    CALL Matrix_AssertIsSPData(matrix,err,error,*999)    
    IF(SIZE(rowIndices,1)/=SIZE(values,1)) THEN
      localError="The size of the row indices array of "//TRIM(NumberToVString(SIZE(rowIndices,1),"*",err,error))// &
        & " does not conform to the number of rows the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(columnIndices,1)/=SIZE(values,2)) THEN
      localError="The size of the column indices array of "//TRIM(NumberToVString(SIZE(columnIndices,1),"*",err,error))// &
        & " does not conform to the number of columns in the values array of "// &
        & TRIM(NumberToVString(SIZE(values,2),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(matrix%storageType)
    CASE(MATRIX_BLOCK_STORAGE_TYPE,MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
      DO j=1,SIZE(columnIndices,1)
        DO i=1,SIZE(rowIndices,1)
          CALL Matrix_StorageLocationGet(matrix,rowIndices(i),columnIndices(j),location,err,error,*999)
          matrix%dataSP(location)=matrix%dataSP(location)+values(i,j)
        ENDDO !i
      ENDDO !j
    CASE(MATRIX_DIAGONAL_STORAGE_TYPE,MATRIX_ROW_MAJOR_STORAGE_TYPE,MATRIX_ROW_COLUMN_STORAGE_TYPE)
      DO i=1,SIZE(rowIndices,1)
        DO j=1,SIZE(columnIndices,1)
          CALL Matrix_StorageLocationGet(matrix,rowIndices(i),columnIndices(j),location,err,error,*999)
          matrix%dataSP(location)=matrix%dataSP(location)+values(i,j)
        ENDDO !j
      ENDDO !i
    CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
      !Search for the column number in the sparsity list using the bisection (binary search) algorithm in order
      !to take advantage of previous searches.
      previousColumnIndex=-1
      DO i=1,SIZE(rowIndices,1)
        rowIndex=rowIndices(i)
        lowLimit=matrix%rowIndices(rowIndex)
        upLimit=matrix%rowIndices(rowIndex+1)
        DO j=1,SIZE(columnIndices,1)
          location=0
          columnIndex=columnIndices(j)
          IF(columnIndex<=previousColumnIndex) THEN
            lowLimit=matrix%rowIndices(rowIndex)
          ELSE
            upLimit=matrix%rowIndices(rowIndex+1)
          ENDIF
          previousColumnIndex=columnIndex
          DO WHILE((upLimit-lowLimit)>bisectionToLinearSearchThreshold)
            midPoint=(upLimit+lowLimit)/2
            IF(matrix%columnIndices(midPoint)>columnIndex) THEN
              upLimit=midPoint
            ELSE
              lowLimit=midPoint
            ENDIF
          ENDDO
          DO k=lowLimit,upLimit
            IF(matrix%columnIndices(k)==columnIndex) THEN
              location=k
              lowLimit=k+1
              EXIT
            ENDIF
          ENDDO !k
          IF(location==0) THEN
            localError="Row "//TRIM(NumberToVString(rowIndices(i),"*",err,error))//" and column "// &
              & TRIM(NumberToVString(columnIndices(j),"*",err,error))//" does not exist in the matrix."
            CALL FlagError(localError,err,error,*999)
          ELSE
            matrix%dataSP(location)=matrix%dataSP(location)+values(i,j)
          ENDIF
        ENDDO !j
      ENDDO !i
    CASE(MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
      !Search for the row number in the sparsity list using the bisection (binary search) algorithm in order
      !to take advantage of previous searches.
      previousRowIndex=-1
      DO j=1,SIZE(columnIndices,1)
        columnIndex=columnIndices(j)
        lowLimit=matrix%columnIndices(columnIndex)
        upLimit=matrix%columnIndices(columnIndex+1)
        DO i=1,SIZE(rowIndices,1)
          location=0
          rowIndex=rowIndices(i)
          IF(rowIndex<=previousRowIndex) THEN
            lowLimit=matrix%columnIndices(columnIndex)
          ELSE
            upLimit=matrix%columnIndices(columnIndex+1)
          ENDIF
          previousRowIndex=rowIndex
          DO WHILE((upLimit-lowLimit)>bisectionToLinearSearchThreshold)
            midPoint=(upLimit+lowLimit)/2
            IF(matrix%rowIndices(midPoint)>rowIndex) THEN
              upLimit=midPoint
            ELSE
              lowLimit=midPoint
            ENDIF
          ENDDO
          DO k=lowLimit,upLimit
            IF(matrix%rowIndices(k)==rowIndex) THEN
              location=k
              lowLimit=k+1
              EXIT
            ENDIF
          ENDDO !k
          IF(location==0) THEN
            localError="Row "//TRIM(NumberToVString(rowIndices(i),"*",err,error))//" and column "// &
              & TRIM(NumberToVString(columnIndices(j),"*",err,error))//" does not exist in the matrix."
            CALL FlagError(localError,err,error,*999)
          ELSE
            matrix%dataSP(location)=matrix%dataSP(location)+values(i,j)
          ENDIF
        ENDDO !i
      ENDDO !j
    CASE(MATRIX_BLOCK_COMPRESSED_ROW_STORAGE_TYPE)
      !Search for the column number in the sparsity list using the bisection (binary search) algorithm in order
      !to take advantage of previous searches.
      previousColumnBlockNumber=-1
      DO i=1,SIZE(rowIndices,1)
        rowNumber=rowIndices(i)
        rowBlockNumber=FLOOR(REAL(rowNumber,DP)/REAL(matrix%blockSize,DP),INTG)
        rowOffset=MOD(rowNumber,matrix%blockSize)
        lowLimit=matrix%rowIndices(rowBlockNumber)
        upLimit=matrix%rowIndices(rowBlockNumber+1)
        DO j=1,SIZE(columnIndices,1)
          columnNumber=columnIndices(j)
          columnBlockNumber=FLOOR(REAL(columnNumber,DP)/REAL(matrix%blockSize,DP),INTG)
          columnOffset=MOD(columnNumber,matrix%blockSize)
          blockLocation=0
          IF(columnBlockNumber<=previousColumnBlockNumber) THEN
            lowLimit=matrix%rowIndices(rowBlockNumber)
          ELSE
            upLimit=matrix%rowIndices(rowBlockNumber+1)
          ENDIF
          previousColumnBlockNumber=columnBlockNumber
          DO WHILE((upLimit-lowLimit)>bisectionToLinearSearchThreshold)
            midPoint=(upLimit+lowLimit)/2
            IF(matrix%columnIndices(midPoint)>columnBlockNumber) THEN
              upLimit=midPoint
            ELSE
              lowLimit=midPoint
            ENDIF
          ENDDO
          DO k=lowLimit,upLimit
            IF(matrix%columnIndices(k)==columnBlockNumber) THEN
              blockLocation=k
              lowLimit=k+1
              EXIT
            ENDIF
          ENDDO !k
          IF(blockLocation==0) THEN
            localError="Row "//TRIM(NumberToVString(rowIndices(i),"*",err,error))//" and column "// &
              & TRIM(NumberToVString(columnIndices(j),"*",err,error))//" does not exist in the matrix."
            CALL FlagError(localError,err,error,*999)
          ELSE
            location=columnBlockNumber*matrix%blockSize*matrix%blockSize+(columnOffset-1)*matrix%blockSize+rowOffset
            matrix%dataSP(location)=matrix%dataSP(location)+values(i,j)
          ENDIF
        ENDDO !j
      ENDDO !i
    CASE(MATRIX_MODIFIED_KRM_STORAGE_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The matrix storage type of "//TRIM(NumberToVString(matrix%storageType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Matrix_ValuesAddSP2")
    RETURN
999 ERRORSEXITS("Matrix_ValuesAddSP2",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_ValuesAddSP2

  !
  !================================================================================================================================
  !

  !>Adds a value to a double precision real matrix at the location specified by the row and column index i.e., matrix(i,j)=matrix(i,j)+value
  SUBROUTINE Matrix_ValuesAddDP(matrix,rowIndex,columnIndex,VALUE,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: rowIndex !<The row index for the value to add
    INTEGER(INTG), INTENT(IN) :: columnIndex !<The column index for the value to add
    REAL(DP), INTENT(IN) :: value !<The value to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: location

    ENTERS("Matrix_ValuesAddDP",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    CALL Matrix_AssertIsDPData(matrix,err,error,*999)

    CALL Matrix_StorageLocationGet(matrix,rowIndex,columnIndex,location,err,error,*999)
    matrix%dataDP(location)=matrix%dataDP(location)+value

    EXITS("Matrix_ValuesAddDP")
    RETURN
999 ERRORSEXITS("Matrix_ValuesAddDP",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_ValuesAddDP

  !
  !================================================================================================================================
  !

  !>Adds values to a double precision real matrix at the location specified by the row and column indices i.e., matrix(i,j)=matrix(i,j)+value
  SUBROUTINE Matrix_ValuesAddDP1(matrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The row index for the i'th value to add
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(i). The column index for the i'th value to add
    REAL(DP), INTENT(IN) :: values(:) !<values(i). The value of the i'th value to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: k,location
    TYPE(VARYING_STRING) :: localError

    ENTERS("Matrix_ValuesAddDP1",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    CALL Matrix_AssertIsDPData(matrix,err,error,*999)
    IF(SIZE(rowIndices,1)/=SIZE(values,1)) THEN
      localError="The size of the row indices array of "//TRIM(NumberToVString(SIZE(rowIndices,1),"*",err,error))// &
        & " does not conform to the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF      
    IF(SIZE(columnIndices,1)/=SIZE(values,1)) THEN
      localError="The size of the column indices array of "//TRIM(NumberToVString(SIZE(columnIndices,1),"*",err,error))// &
        & " does not conform to the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    DO k=1,SIZE(rowIndices,1)
      CALL Matrix_StorageLocationGet(matrix,rowIndices(k),columnIndices(k),location,err,error,*999)
      matrix%dataDP(location)=matrix%dataDP(location)+values(k)
    ENDDO !k
    
    EXITS("Matrix_ValuesAddDP1")
    RETURN
999 ERRORSEXITS("Matrix_ValuesAddDP1",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_ValuesAddDP1

  !
  !================================================================================================================================
  !

  !>Adds a matrix of values to a double precision real matrix at the location specified by the row and column indices i.e., matrix(i,j)=matrix(i,j)+value
  SUBROUTINE Matrix_ValuesAddDP2(matrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The row index for the ij'th value to add
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(j). The column index for the ij'th value to add
    REAL(DP), INTENT(IN) :: values(:,:) !<values(i,j). The value of the ij'th value to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: blockLocation,i,j,k,rowIndex,previousRowIndex,columnIndex,columnNumber,columnBlockNumber,columnOffset, &
      & previousColumnIndex,location,lowLimit,midPoint,previousColumnBlockNumber,rowBlockNumber,rowNumber,rowOffset,upLimit
    TYPE(VARYING_STRING) :: localError

    ENTERS("Matrix_ValuesAddDP2",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    CALL Matrix_AssertIsDPData(matrix,err,error,*999)    
    IF(SIZE(rowIndices,1)/=SIZE(values,1)) THEN
      localError="The size of the row indices array of "//TRIM(NumberToVString(SIZE(rowIndices,1),"*",err,error))// &
        & " does not conform to the number of rows the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(columnIndices,1)/=SIZE(values,2)) THEN
      localError="The size of the column indices array of "//TRIM(NumberToVString(SIZE(columnIndices,1),"*",err,error))// &
        & " does not conform to the number of columns in the values array of "// &
        & TRIM(NumberToVString(SIZE(values,2),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(matrix%storageType)
    CASE(MATRIX_BLOCK_STORAGE_TYPE,MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
      DO j=1,SIZE(columnIndices,1)
        DO i=1,SIZE(rowIndices,1)
          CALL Matrix_StorageLocationGet(matrix,rowIndices(i),columnIndices(j),location,err,error,*999)
          matrix%dataDP(location)=matrix%dataDP(location)+values(i,j)
        ENDDO !i
      ENDDO !j
    CASE(MATRIX_DIAGONAL_STORAGE_TYPE,MATRIX_ROW_MAJOR_STORAGE_TYPE,MATRIX_ROW_COLUMN_STORAGE_TYPE)
      DO i=1,SIZE(rowIndices,1)
        DO j=1,SIZE(columnIndices,1)
          CALL Matrix_StorageLocationGet(matrix,rowIndices(i),columnIndices(j),location,err,error,*999)
          matrix%dataDP(location)=matrix%dataDP(location)+values(i,j)
        ENDDO !j
      ENDDO !i
    CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
      !Search for the column number in the sparsity list using the bisection (binary search) algorithm in order
      !to take advantage of previous searches.
      previousColumnIndex=-1
      DO i=1,SIZE(rowIndices,1)
        rowIndex=rowIndices(i)
        lowLimit=matrix%rowIndices(rowIndex)
        upLimit=matrix%rowIndices(rowIndex+1)
        DO j=1,SIZE(columnIndices,1)
          location=0
          columnIndex=columnIndices(j)
          IF(columnIndex<=previousColumnIndex) THEN
            lowLimit=matrix%rowIndices(rowIndex)
          ELSE
            upLimit=matrix%rowIndices(rowIndex+1)
          ENDIF
          previousColumnIndex=columnIndex
          DO WHILE((upLimit-lowLimit)>bisectionToLinearSearchThreshold)
            midPoint=(upLimit+lowLimit)/2
            IF(matrix%columnIndices(midPoint)>columnIndex) THEN
              upLimit=midPoint
            ELSE
              lowLimit=midPoint
            ENDIF
          ENDDO
          DO k=lowLimit,upLimit
            IF(matrix%columnIndices(k)==columnIndex) THEN
              location=k
              lowLimit=k+1
              EXIT
            ENDIF
          ENDDO !k
          IF(location==0) THEN
            localError="Row "//TRIM(NumberToVString(rowIndices(i),"*",err,error))//" and column "// &
              & TRIM(NumberToVString(columnIndices(j),"*",err,error))//" does not exist in the matrix."
            CALL FlagError(localError,err,error,*999)
          ELSE
            matrix%dataDP(location)=matrix%dataDP(location)+values(i,j)
          ENDIF
        ENDDO !j
      ENDDO !i
    CASE(MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
      !Search for the row number in the sparsity list using the bisection (binary search) algorithm in order
      !to take advantage of previous searches.
      previousRowIndex=-1
      DO j=1,SIZE(columnIndices,1)
        columnIndex=columnIndices(j)
        lowLimit=matrix%columnIndices(columnIndex)
        upLimit=matrix%columnIndices(columnIndex+1)
        DO i=1,SIZE(rowIndices,1)
          location=0
          rowIndex=rowIndices(i)
          IF(rowIndex<=previousRowIndex) THEN
            lowLimit=matrix%columnIndices(columnIndex)
          ELSE
            upLimit=matrix%columnIndices(columnIndex+1)
          ENDIF
          previousRowIndex=rowIndex
          DO WHILE((upLimit-lowLimit)>bisectionToLinearSearchThreshold)
            midPoint=(upLimit+lowLimit)/2
            IF(matrix%rowIndices(midPoint)>rowIndex) THEN
              upLimit=midPoint
            ELSE
              lowLimit=midPoint
            ENDIF
          ENDDO
          DO k=lowLimit,upLimit
            IF(matrix%rowIndices(k)==rowIndex) THEN
              location=k
              lowLimit=k+1
              EXIT
            ENDIF
          ENDDO !k
          IF(location==0) THEN
            localError="Row "//TRIM(NumberToVString(rowIndices(i),"*",err,error))//" and column "// &
              & TRIM(NumberToVString(columnIndices(j),"*",err,error))//" does not exist in the matrix."
            CALL FlagError(localError,err,error,*999)
          ELSE
            matrix%dataDP(location)=matrix%dataDP(location)+values(i,j)
          ENDIF
        ENDDO !i
      ENDDO !j
    CASE(MATRIX_BLOCK_COMPRESSED_ROW_STORAGE_TYPE)
      !Search for the column number in the sparsity list using the bisection (binary search) algorithm in order
      !to take advantage of previous searches.
      previousColumnBlockNumber=-1
      DO i=1,SIZE(rowIndices,1)
        rowNumber=rowIndices(i)
        rowBlockNumber=FLOOR(REAL(rowNumber,DP)/REAL(matrix%blockSize,DP),INTG)
        rowOffset=MOD(rowNumber,matrix%blockSize)
        lowLimit=matrix%rowIndices(rowBlockNumber)
        upLimit=matrix%rowIndices(rowBlockNumber+1)
        DO j=1,SIZE(columnIndices,1)
          columnNumber=columnIndices(j)
          columnBlockNumber=FLOOR(REAL(columnNumber,DP)/REAL(matrix%blockSize,DP),INTG)
          columnOffset=MOD(columnNumber,matrix%blockSize)
          blockLocation=0
          IF(columnBlockNumber<=previousColumnBlockNumber) THEN
            lowLimit=matrix%rowIndices(rowBlockNumber)
          ELSE
            upLimit=matrix%rowIndices(rowBlockNumber+1)
          ENDIF
          previousColumnBlockNumber=columnBlockNumber
          DO WHILE((upLimit-lowLimit)>bisectionToLinearSearchThreshold)
            midPoint=(upLimit+lowLimit)/2
            IF(matrix%columnIndices(midPoint)>columnBlockNumber) THEN
              upLimit=midPoint
            ELSE
              lowLimit=midPoint
            ENDIF
          ENDDO
          DO k=lowLimit,upLimit
            IF(matrix%columnIndices(k)==columnBlockNumber) THEN
              blockLocation=k
              lowLimit=k+1
              EXIT
            ENDIF
          ENDDO !k
          IF(blockLocation==0) THEN
            localError="Row "//TRIM(NumberToVString(rowIndices(i),"*",err,error))//" and column "// &
              & TRIM(NumberToVString(columnIndices(j),"*",err,error))//" does not exist in the matrix."
            CALL FlagError(localError,err,error,*999)
          ELSE
            location=columnBlockNumber*matrix%blockSize*matrix%blockSize+(columnOffset-1)*matrix%blockSize+rowOffset
            matrix%dataDP(location)=matrix%dataDP(location)+values(i,j)
          ENDIF
        ENDDO !j
      ENDDO !i
    CASE(MATRIX_MODIFIED_KRM_STORAGE_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The matrix storage type of "//TRIM(NumberToVString(matrix%storageType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Matrix_ValuesAddDP2")
    RETURN
999 ERRORSEXITS("Matrix_ValuesAddDP2",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_ValuesAddDP2

  !
  !================================================================================================================================
  !

  !>Adds a value to a logical matrix at the location specified by the row and column index i.e., matrix(i,j)=matrix(i,j)+value
  SUBROUTINE Matrix_ValuesAddL(matrix,rowIndex,columnIndex,VALUE,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: rowIndex !<The row index for the value to add
    INTEGER(INTG), INTENT(IN) :: columnIndex !<The column index for the value to add
    LOGICAL, INTENT(IN) :: value !<The value to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: location

    ENTERS("Matrix_ValuesAddL",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    CALL Matrix_AssertIsLData(matrix,err,error,*999)

    CALL Matrix_StorageLocationGet(matrix,rowIndex,columnIndex,location,err,error,*999)
    matrix%dataL(location)=matrix%dataL(location).OR.value

    EXITS("Matrix_ValuesAddL")
    RETURN
999 ERRORSEXITS("Matrix_ValuesAddL",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_ValuesAddL

  !
  !================================================================================================================================
  !

  !>Adds values to a logical matrix at the location specified by the row and column indices i.e., matrix(i,j)=matrix(i,j)+value
  SUBROUTINE Matrix_ValuesAddL1(matrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The row index for the i'th value to add
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(i). The column index for the i'th value to add
    LOGICAL, INTENT(IN) :: values(:) !<values(i). The value of the i'th value to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: k,location
    TYPE(VARYING_STRING) :: localError

    ENTERS("Matrix_ValuesAddL1",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    CALL Matrix_AssertIsLData(matrix,err,error,*999)
    IF(SIZE(rowIndices,1)/=SIZE(values,1)) THEN
      localError="The size of the row indices array of "//TRIM(NumberToVString(SIZE(rowIndices,1),"*",err,error))// &
        & " does not conform to the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF      
    IF(SIZE(columnIndices,1)/=SIZE(values,1)) THEN
      localError="The size of the column indices array of "//TRIM(NumberToVString(SIZE(columnIndices,1),"*",err,error))// &
        & " does not conform to the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    DO k=1,SIZE(rowIndices,1)
      CALL Matrix_StorageLocationGet(matrix,rowIndices(k),columnIndices(k),location,err,error,*999)
      matrix%dataL(location)=matrix%dataL(location).OR.values(k)
    ENDDO !k
    
    EXITS("Matrix_ValuesAddL1")
    RETURN
999 ERRORSEXITS("Matrix_ValuesAddL1",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_ValuesAddL1

  !
  !================================================================================================================================
  !

  !>Adds a matrix of values to a logical matrix at the location specified by the row and column indices i.e., matrix(i,j)=matrix(i,j)+value
  SUBROUTINE Matrix_ValuesAddL2(matrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The row index for the ij'th value to add
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(j). The column index for the ij'th value to add
    LOGICAL, INTENT(IN) :: values(:,:) !<values(i,j). The value of the ij'th value to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: blockLocation,i,j,k,rowIndex,previousRowIndex,columnIndex,columnNumber,columnBlockNumber,columnOffset, &
      & previousColumnIndex,location,lowLimit,midPoint,previousColumnBlockNumber,rowBlockNumber,rowNumber,rowOffset,upLimit
    TYPE(VARYING_STRING) :: localError

    ENTERS("Matrix_ValuesAddL2",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    CALL Matrix_AssertIsLData(matrix,err,error,*999)    
    IF(SIZE(rowIndices,1)/=SIZE(values,1)) THEN
      localError="The size of the row indices array of "//TRIM(NumberToVString(SIZE(rowIndices,1),"*",err,error))// &
        & " does not conform to the number of rows the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(columnIndices,1)/=SIZE(values,2)) THEN
      localError="The size of the column indices array of "//TRIM(NumberToVString(SIZE(columnIndices,1),"*",err,error))// &
        & " does not conform to the number of columns in the values array of "// &
        & TRIM(NumberToVString(SIZE(values,2),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(matrix%storageType)
    CASE(MATRIX_BLOCK_STORAGE_TYPE,MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
      DO j=1,SIZE(columnIndices,1)
        DO i=1,SIZE(rowIndices,1)
          CALL Matrix_StorageLocationGet(matrix,rowIndices(i),columnIndices(j),location,err,error,*999)
          matrix%dataL(location)=matrix%dataL(location).OR.values(i,j)
        ENDDO !i
      ENDDO !j
    CASE(MATRIX_DIAGONAL_STORAGE_TYPE,MATRIX_ROW_MAJOR_STORAGE_TYPE,MATRIX_ROW_COLUMN_STORAGE_TYPE)
      DO i=1,SIZE(rowIndices,1)
        DO j=1,SIZE(columnIndices,1)
          CALL Matrix_StorageLocationGet(matrix,rowIndices(i),columnIndices(j),location,err,error,*999)
          matrix%dataL(location)=matrix%dataL(location).OR.values(i,j)
        ENDDO !j
      ENDDO !i
    CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
      !Search for the column number in the sparsity list using the bisection (binary search) algorithm in order
      !to take advantage of previous searches.
      previousColumnIndex=-1
      DO i=1,SIZE(rowIndices,1)
        rowIndex=rowIndices(i)
        lowLimit=matrix%rowIndices(rowIndex)
        upLimit=matrix%rowIndices(rowIndex+1)
        DO j=1,SIZE(columnIndices,1)
          location=0
          columnIndex=columnIndices(j)
          IF(columnIndex<=previousColumnIndex) THEN
            lowLimit=matrix%rowIndices(rowIndex)
          ELSE
            upLimit=matrix%rowIndices(rowIndex+1)
          ENDIF
          previousColumnIndex=columnIndex
          DO WHILE((upLimit-lowLimit)>bisectionToLinearSearchThreshold)
            midPoint=(upLimit+lowLimit)/2
            IF(matrix%columnIndices(midPoint)>columnIndex) THEN
              upLimit=midPoint
            ELSE
              lowLimit=midPoint
            ENDIF
          ENDDO
          DO k=lowLimit,upLimit
            IF(matrix%columnIndices(k)==columnIndex) THEN
              location=k
              lowLimit=k+1
              EXIT
            ENDIF
          ENDDO !k
          IF(location==0) THEN
            localError="Row "//TRIM(NumberToVString(rowIndices(i),"*",err,error))//" and column "// &
              & TRIM(NumberToVString(columnIndices(j),"*",err,error))//" does not exist in the matrix."
            CALL FlagError(localError,err,error,*999)
          ELSE
            matrix%dataL(location)=matrix%dataL(location).OR.values(i,j)
          ENDIF
        ENDDO !j
      ENDDO !i
    CASE(MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
      !Search for the row number in the sparsity list using the bisection (binary search) algorithm in order
      !to take advantage of previous searches.
      previousRowIndex=-1
      DO j=1,SIZE(columnIndices,1)
        columnIndex=columnIndices(j)
        lowLimit=matrix%columnIndices(columnIndex)
        upLimit=matrix%columnIndices(columnIndex+1)
        DO i=1,SIZE(rowIndices,1)
          location=0
          rowIndex=rowIndices(i)
          IF(rowIndex<=previousRowIndex) THEN
            lowLimit=matrix%columnIndices(columnIndex)
          ELSE
            upLimit=matrix%columnIndices(columnIndex+1)
          ENDIF
          previousRowIndex=rowIndex
          DO WHILE((upLimit-lowLimit)>bisectionToLinearSearchThreshold)
            midPoint=(upLimit+lowLimit)/2
            IF(matrix%rowIndices(midPoint)>rowIndex) THEN
              upLimit=midPoint
            ELSE
              lowLimit=midPoint
            ENDIF
          ENDDO
          DO k=lowLimit,upLimit
            IF(matrix%rowIndices(k)==rowIndex) THEN
              location=k
              lowLimit=k+1
              EXIT
            ENDIF
          ENDDO !k
          IF(location==0) THEN
            localError="Row "//TRIM(NumberToVString(rowIndices(i),"*",err,error))//" and column "// &
              & TRIM(NumberToVString(columnIndices(j),"*",err,error))//" does not exist in the matrix."
            CALL FlagError(localError,err,error,*999)
          ELSE
            matrix%dataL(location)=matrix%dataL(location).OR.values(i,j)
          ENDIF
        ENDDO !i
      ENDDO !j
    CASE(MATRIX_BLOCK_COMPRESSED_ROW_STORAGE_TYPE)
      !Search for the column number in the sparsity list using the bisection (binary search) algorithm in order
      !to take advantage of previous searches.
      previousColumnBlockNumber=-1
      DO i=1,SIZE(rowIndices,1)
        rowNumber=rowIndices(i)
        rowBlockNumber=FLOOR(REAL(rowNumber,DP)/REAL(matrix%blockSize,DP),INTG)
        rowOffset=MOD(rowNumber,matrix%blockSize)
        lowLimit=matrix%rowIndices(rowBlockNumber)
        upLimit=matrix%rowIndices(rowBlockNumber+1)
        DO j=1,SIZE(columnIndices,1)
          columnNumber=columnIndices(j)
          columnBlockNumber=FLOOR(REAL(columnNumber,DP)/REAL(matrix%blockSize,DP),INTG)
          columnOffset=MOD(columnNumber,matrix%blockSize)
          blockLocation=0
          IF(columnBlockNumber<=previousColumnBlockNumber) THEN
            lowLimit=matrix%rowIndices(rowBlockNumber)
          ELSE
            upLimit=matrix%rowIndices(rowBlockNumber+1)
          ENDIF
          previousColumnBlockNumber=columnBlockNumber
          DO WHILE((upLimit-lowLimit)>bisectionToLinearSearchThreshold)
            midPoint=(upLimit+lowLimit)/2
            IF(matrix%columnIndices(midPoint)>columnBlockNumber) THEN
              upLimit=midPoint
            ELSE
              lowLimit=midPoint
            ENDIF
          ENDDO
          DO k=lowLimit,upLimit
            IF(matrix%columnIndices(k)==columnBlockNumber) THEN
              blockLocation=k
              lowLimit=k+1
              EXIT
            ENDIF
          ENDDO !k
          IF(blockLocation==0) THEN
            localError="Row "//TRIM(NumberToVString(rowIndices(i),"*",err,error))//" and column "// &
              & TRIM(NumberToVString(columnIndices(j),"*",err,error))//" does not exist in the matrix."
            CALL FlagError(localError,err,error,*999)
          ELSE
            location=columnBlockNumber*matrix%blockSize*matrix%blockSize+(columnOffset-1)*matrix%blockSize+rowOffset
            matrix%dataL(location)=matrix%dataL(location).OR.values(i,j)
          ENDIF
        ENDDO !j
      ENDDO !i
    CASE(MATRIX_MODIFIED_KRM_STORAGE_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The matrix storage type of "//TRIM(NumberToVString(matrix%storageType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Matrix_ValuesAddL2")
    RETURN
999 ERRORSEXITS("Matrix_ValuesAddL2",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_ValuesAddL2

  !
  !================================================================================================================================
  !

  !>Gets a value in an integer matrix at the location specified by the row and column index i.e., value=matrix(i,j)
  SUBROUTINE Matrix_ValuesGetIntg(matrix,rowIndex,columnIndex,value,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: rowIndex !<The row index of the value to get
    INTEGER(INTG), INTENT(IN) :: columnIndex !<The column index of the value to get
    INTEGER(INTG), INTENT(OUT) :: value !<On return the value in the matrix at the specified row and column
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: location

    ENTERS("Matrix_ValuesGetIntg",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    CALL Matrix_AssertIsIntgData(matrix,err,error,*999)

    CALL Matrix_StorageLocationFind(matrix,rowIndex,columnIndex,location,err,error,*999)
    IF(location==0) THEN
      value=0
    ELSE
      value=matrix%dataIntg(location)
    ENDIF

    EXITS("Matrix_ValuesGetIntg")
    RETURN
999 ERRORSEXITS("Matrix_ValuesGetIntg",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_ValuesGetIntg

  !
  !================================================================================================================================
  !

  !>Gets the values in an integer matrix at the location specified by the row and column indices i.e., value=matrix(i,j)
  SUBROUTINE Matrix_ValuesGetIntg1(matrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The row index for the i'th value to get
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(i). The column index for the i'th value to get
    INTEGER(INTG), INTENT(OUT) :: values(:) !<values(i). On return the value of the i'th value to get
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: k,location
    TYPE(VARYING_STRING) :: localError

    ENTERS("Matrix_ValuesGetIntg1",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    CALL Matrix_AssertIsIntgData(matrix,err,error,*999)
    IF(SIZE(rowIndices,1)/=SIZE(values,1)) THEN
      localError="The size of the row indices array of "//TRIM(NumberToVString(SIZE(rowIndices,1),"*",err,error))// &
        & " does not conform to the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(columnIndices,1)/=SIZE(values,1)) THEN
      localError="The size of the column indices array of "//TRIM(NumberToVString(SIZE(columnIndices,1),"*",err,error))// &
        & " does not conform to the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
      
    DO k=1,SIZE(rowIndices,1)
      CALL Matrix_StorageLocationFind(matrix,rowIndices(k),columnIndices(k),location,err,error,*999)
      IF(location==0) THEN
        values(k)=0
      ELSE
        values(k)=matrix%dataIntg(location)
      ENDIF
    ENDDO !k

    EXITS("Matrix_ValuesGetIntg1")
    RETURN
999 ERRORSEXITS("Matrix_ValuesGetIntg1",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_ValuesGetIntg1

  !
  !================================================================================================================================
  !

  !>Gets the matrix of values in an integer matrix at the location specified by the row and column indices i.e., value=matrix(i,j)
  SUBROUTINE Matrix_ValuesGetIntg2(matrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The row index for the ij'th value to get
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(j). The column index for the ij'th value to get
    INTEGER(INTG), INTENT(OUT) :: values(:,:) !<values(i,j). On return the value of the ij'th value to get
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,j,location
    TYPE(VARYING_STRING) :: localError

    ENTERS("Matrix_ValuesGetIntg2",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    CALL Matrix_AssertIsIntgData(matrix,err,error,*999)
    IF(SIZE(rowIndices,1)/=SIZE(values,1)) THEN
      localError="The size of the row indices array of "//TRIM(NumberToVString(SIZE(rowIndices,1),"*",err,error))// &
        & " does not conform to the number of rows in the values array of"// &
        & TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF      
    IF(SIZE(columnIndices,1)/=SIZE(values,2)) THEN
      localError="The size of the row indices array of "//TRIM(NumberToVString(SIZE(rowIndices,1),"*",err,error))// &
        & " does not conform to the number of rows in the values array of "// &
        & TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
      
    DO i=1,SIZE(rowIndices,1)
      DO j=1,SIZE(columnIndices,1)
        CALL Matrix_StorageLocationFind(matrix,rowIndices(i),columnIndices(j),location,err,error,*999)
        IF(location==0) THEN
          values(i,j)=0
        ELSE
          values(i,j)=matrix%dataIntg(location)
        ENDIF
      ENDDO !j
    ENDDO !i
    
    EXITS("Matrix_ValuesGetIntg2")
    RETURN
999 ERRORSEXITS("Matrix_ValuesGetIntg2",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_ValuesGetIntg2

  !
  !================================================================================================================================
  !

  !>Gets a value in a single precision real matrix at the location specified by the row and column index i.e., value=matrix(i,j)
  SUBROUTINE Matrix_ValuesGetSP(matrix,rowIndex,columnIndex,VALUE,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: rowIndex !<The row index of the value to get
    INTEGER(INTG), INTENT(IN) :: columnIndex !<The column index of the value to get
    REAL(SP), INTENT(OUT) :: value !<On return the value in the matrix at the specified row and column
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: location

    ENTERS("Matrix_ValuesGetSP",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    CALL Matrix_AssertIsSPData(matrix,err,error,*999)

    CALL Matrix_StorageLocationFind(matrix,rowIndex,columnIndex,location,err,error,*999)
    IF(location==0) THEN
      value=0.0_SP
    ELSE
      value=matrix%dataSP(location)
    ENDIF

    EXITS("Matrix_ValuesGetSP")
    RETURN
999 ERRORSEXITS("Matrix_ValuesGetSP",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_ValuesGetSP

  !
  !================================================================================================================================
  !

  !>Gets the values in a single precision real matrix at the location specified by the row and column indices i.e., value=matrix(i,j)
  SUBROUTINE Matrix_ValuesGetSP1(matrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The row index for the i'th value to get
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(i). The column index for the i'th value to get
    REAL(SP), INTENT(OUT) :: values(:) !<values(i). On return the value of the i'th value to get
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: k,location
    TYPE(VARYING_STRING) :: localError

    ENTERS("Matrix_ValuesGetSP1",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    CALL Matrix_AssertIsSPData(matrix,err,error,*999)
    IF(SIZE(rowIndices,1)/=SIZE(values,1)) THEN
      localError="The size of the row indices array of "//TRIM(NumberToVString(SIZE(rowIndices,1),"*",err,error))// &
        & " does not conform to the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(columnIndices,1)/=SIZE(values,1)) THEN
      localError="The size of the column indices array of "//TRIM(NumberToVString(SIZE(columnIndices,1),"*",err,error))// &
        & " does not conform to the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
      
    DO k=1,SIZE(rowIndices,1)
      CALL Matrix_StorageLocationFind(matrix,rowIndices(k),columnIndices(k),location,err,error,*999)
      IF(location==0) THEN
        values(k)=0.0_SP
      ELSE
        values(k)=matrix%dataSP(location)
      ENDIF
    ENDDO !k

    EXITS("Matrix_ValuesGetSP1")
    RETURN
999 ERRORSEXITS("Matrix_ValuesGetSP1",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_ValuesGetSP1
  
  !
  !================================================================================================================================
  !

  !>Gets the matrix of values in a single precision real matrix at the location specified by the row and column indices i.e., value=matrix(i,j)
  SUBROUTINE Matrix_ValuesGetSP2(matrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The row index for the ij'th value to get
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(j). The column index for the ij'th value to get
    REAL(SP), INTENT(OUT) :: values(:,:) !<values(i,j). On return the value of the ij'th value to get
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,j,location
    TYPE(VARYING_STRING) :: localError

    ENTERS("Matrix_ValuesGetSP2",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    CALL Matrix_AssertIsSPData(matrix,err,error,*999)
    IF(SIZE(rowIndices,1)/=SIZE(values,1)) THEN
      localError="The size of the row indices array of "//TRIM(NumberToVString(SIZE(rowIndices,1),"*",err,error))// &
        & " does not conform to the number of rows in the values array of"// &
        & TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF      
    IF(SIZE(columnIndices,1)/=SIZE(values,2)) THEN
      localError="The size of the row indices array of "//TRIM(NumberToVString(SIZE(rowIndices,1),"*",err,error))// &
        & " does not conform to the number of rows in the values array of "// &
        & TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
      
    DO i=1,SIZE(rowIndices,1)
      DO j=1,SIZE(columnIndices,1)
        CALL Matrix_StorageLocationFind(matrix,rowIndices(i),columnIndices(j),location,err,error,*999)
        IF(location==0) THEN
          values(i,j)=0.0_SP
        ELSE
          values(i,j)=matrix%dataSP(location)
        ENDIF
      ENDDO !j
    ENDDO !i
    
    EXITS("Matrix_ValuesGetSP2")
    RETURN
999 ERRORSEXITS("Matrix_ValuesGetSP2",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_ValuesGetSP2

  !
  !================================================================================================================================
  !

  !>Gets a value in a double precision real matrix at the location specified by the row and column index i.e., value=matrix(i,j)
  SUBROUTINE Matrix_ValuesGetDP(matrix,rowIndex,columnIndex,VALUE,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: rowIndex !<The row index of the value to get
    INTEGER(INTG), INTENT(IN) :: columnIndex !<The column index of the value to get
    REAL(DP), INTENT(OUT) :: value !<On return the value in the matrix at the specified row and column
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: location

    ENTERS("Matrix_ValuesGetDP",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    CALL Matrix_AssertIsDPData(matrix,err,error,*999)

    CALL Matrix_StorageLocationFind(matrix,rowIndex,columnIndex,location,err,error,*999)
    IF(location==0) THEN
      value=0.0_DP
    ELSE
      value=matrix%dataDP(location)
    ENDIF

    EXITS("Matrix_ValuesGetDP")
    RETURN
999 ERRORSEXITS("Matrix_ValuesGetDP",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_ValuesGetDP

  !
  !================================================================================================================================
  !

  !>Gets the values in a double precision real matrix at the location specified by the row and column indices i.e., value=matrix(i,j)
  SUBROUTINE Matrix_ValuesGetDP1(matrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The row index for the i'th value to get
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(i). The column index for the i'th value to get
    REAL(DP), INTENT(OUT) :: values(:) !<values(i). On return the value of the i'th value to get
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: k,location
    TYPE(VARYING_STRING) :: localError

    ENTERS("Matrix_ValuesGetDP1",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    CALL Matrix_AssertIsDPData(matrix,err,error,*999)
    IF(SIZE(rowIndices,1)/=SIZE(values,1)) THEN
      localError="The size of the row indices array of "//TRIM(NumberToVString(SIZE(rowIndices,1),"*",err,error))// &
        & " does not conform to the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(columnIndices,1)/=SIZE(values,1)) THEN
      localError="The size of the column indices array of "//TRIM(NumberToVString(SIZE(columnIndices,1),"*",err,error))// &
        & " does not conform to the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
      
    DO k=1,SIZE(rowIndices,1)
      CALL Matrix_StorageLocationFind(matrix,rowIndices(k),columnIndices(k),location,err,error,*999)
      IF(location==0) THEN
        values(k)=0.0_DP
      ELSE
        values(k)=matrix%dataDP(location)
      ENDIF
    ENDDO !k

    EXITS("Matrix_ValuesGetDP1")
    RETURN
999 ERRORSEXITS("Matrix_ValuesGetDP1",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_ValuesGetDP1
  
  !
  !================================================================================================================================
  !

  !>Gets the matrix of values in a double precision real matrix at the location specified by the row and column indices i.e., value=matrix(i,j)
  SUBROUTINE Matrix_ValuesGetDP2(matrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The row index for the ij'th value to get
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(j). The column index for the ij'th value to get
    REAL(DP), INTENT(OUT) :: values(:,:) !<values(i,j). On return the value of the ij'th value to get
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,j,location
    TYPE(VARYING_STRING) :: localError

    ENTERS("Matrix_ValuesGetDP2",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    CALL Matrix_AssertIsDPData(matrix,err,error,*999)
    IF(SIZE(rowIndices,1)/=SIZE(values,1)) THEN
      localError="The size of the row indices array of "//TRIM(NumberToVString(SIZE(rowIndices,1),"*",err,error))// &
        & " does not conform to the number of rows in the values array of"// &
        & TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF      
    IF(SIZE(columnIndices,1)/=SIZE(values,2)) THEN
      localError="The size of the row indices array of "//TRIM(NumberToVString(SIZE(rowIndices,1),"*",err,error))// &
        & " does not conform to the number of rows in the values array of "// &
        & TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
      
    DO i=1,SIZE(rowIndices,1)
      DO j=1,SIZE(columnIndices,1)
        CALL Matrix_StorageLocationFind(matrix,rowIndices(i),columnIndices(j),location,err,error,*999)
        IF(location==0) THEN
          values(i,j)=0.0_DP
        ELSE
          values(i,j)=matrix%dataDP(location)
        ENDIF
      ENDDO !j
    ENDDO !i
    
    EXITS("Matrix_ValuesGetDP2")
    RETURN
999 ERRORSEXITS("Matrix_ValuesGetDP2",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_ValuesGetDP2

  !
  !================================================================================================================================
  !

  !>Gets a value in a logical matrix at the location specified by the row and column index i.e., value=matrix(i,j)
  SUBROUTINE Matrix_ValuesGetL(matrix,rowIndex,columnIndex,VALUE,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: rowIndex !<The row index of the value to get
    INTEGER(INTG), INTENT(IN) :: columnIndex !<The column index of the value to get
    LOGICAL, INTENT(OUT) :: value !<On return the value in the matrix at the specified row and column
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: location

    ENTERS("Matrix_ValuesGetL",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    CALL Matrix_AssertIsLData(matrix,err,error,*999)

    CALL Matrix_StorageLocationFind(matrix,rowIndex,columnIndex,location,err,error,*999)
    IF(location==0) THEN
      value=.FALSE.
    ELSE
      value=matrix%dataL(location)
    ENDIF

    EXITS("Matrix_ValuesGetL")
    RETURN
999 ERRORSEXITS("Matrix_ValuesGetL",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_ValuesGetL

  !
  !================================================================================================================================
  !

  !>Gets the values in a logical matrix at the location specified by the row and column indices i.e., value=matrix(i,j)
  SUBROUTINE Matrix_ValuesGetL1(matrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The row index for the i'th value to get
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(i). The column index for the i'th value to get
    LOGICAL, INTENT(OUT) :: values(:) !<values(i). On return the value of the i'th value to get
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: k,location
    TYPE(VARYING_STRING) :: localError

    ENTERS("Matrix_ValuesGetL1",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    CALL Matrix_AssertIsLData(matrix,err,error,*999)
    IF(SIZE(rowIndices,1)/=SIZE(values,1)) THEN
      localError="The size of the row indices array of "//TRIM(NumberToVString(SIZE(rowIndices,1),"*",err,error))// &
        & " does not conform to the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(columnIndices,1)/=SIZE(values,1)) THEN
      localError="The size of the column indices array of "//TRIM(NumberToVString(SIZE(columnIndices,1),"*",err,error))// &
        & " does not conform to the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
      
    DO k=1,SIZE(rowIndices,1)
      CALL Matrix_StorageLocationFind(matrix,rowIndices(k),columnIndices(k),location,err,error,*999)
      IF(location==0) THEN
        values(k)=.FALSE.
      ELSE
        values(k)=matrix%dataL(location)
      ENDIF
    ENDDO !k

    EXITS("Matrix_ValuesGetL1")
    RETURN
999 ERRORSEXITS("Matrix_ValuesGetL1",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_ValuesGetL1
  
  !
  !================================================================================================================================
  !

  !>Gets the matrix of values in a logical matrix at the location specified by the row and column indices i.e., value=matrix(i,j)
  SUBROUTINE Matrix_ValuesGetL2(matrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The row index for the ij'th value to get
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(j). The column index for the ij'th value to get
    LOGICAL, INTENT(OUT) :: values(:,:) !<values(i,j). On return the value of the ij'th value to get
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,j,location
    TYPE(VARYING_STRING) :: localError

    ENTERS("Matrix_ValuesGetL2",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    CALL Matrix_AssertIsLData(matrix,err,error,*999)
    IF(SIZE(rowIndices,1)/=SIZE(values,1)) THEN
      localError="The size of the row indices array of "//TRIM(NumberToVString(SIZE(rowIndices,1),"*",err,error))// &
        & " does not conform to the number of rows in the values array of"// &
        & TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF      
    IF(SIZE(columnIndices,1)/=SIZE(values,2)) THEN
      localError="The size of the row indices array of "//TRIM(NumberToVString(SIZE(rowIndices,1),"*",err,error))// &
        & " does not conform to the number of rows in the values array of "// &
        & TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
      
    DO i=1,SIZE(rowIndices,1)
      DO j=1,SIZE(columnIndices,1)
        CALL Matrix_StorageLocationFind(matrix,rowIndices(i),columnIndices(j),location,err,error,*999)
        IF(location==0) THEN
          values(i,j)=.FALSE.
        ELSE
          values(i,j)=matrix%dataL(location)
        ENDIF
      ENDDO !j
    ENDDO !i
    
    EXITS("Matrix_ValuesGetL2")
    RETURN
999 ERRORSEXITS("Matrix_ValuesGetL2",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_ValuesGetL2

  !
  !================================================================================================================================
  !

  !>Sets a value in an integer matrix at the location specified by the row and column index i.e., matrix(i,j)=value
  SUBROUTINE Matrix_ValuesSetIntg(matrix,rowIndex,columnIndex,value,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: rowIndex !<The row index of the value to set
    INTEGER(INTG), INTENT(IN) :: columnIndex !<The column index of the value to set
    INTEGER(INTG), INTENT(IN) :: value !<The value to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: location

    ENTERS("Matrix_ValuesSetIntg",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    CALL Matrix_AssertIsIntgData(matrix,err,error,*999)
    
    CALL Matrix_StorageLocationGet(matrix,rowIndex,columnIndex,location,err,error,*999)
    matrix%dataIntg(location)=value

    EXITS("Matrix_ValuesSetIntg")
    RETURN
999 ERRORSEXITS("Matrix_ValuesSetIntg",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_ValuesSetIntg

  !
  !================================================================================================================================
  !

  !>Sets the values in an integer matrix at the location specified by the row and column indices i.e., matrix(i,j)=value
  SUBROUTINE Matrix_ValuesSetIntg1(matrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The row index for the i'th value to set
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(i). The column index for the i'th value to set
    INTEGER(INTG), INTENT(IN) :: values(:) !<values(i). The value of the i'th value to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: k,location
    TYPE(VARYING_STRING) :: localError

    ENTERS("Matrix_ValuesSetIntg1",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    CALL Matrix_AssertIsIntgData(matrix,err,error,*999)
    IF(SIZE(rowIndices,1)/=SIZE(values,1)) THEN
      localError="The size of the row indices array of "//TRIM(NumberToVString(SIZE(rowIndices,1),"*",err,error))// &
        & " does not conform to the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(columnIndices,1)/=SIZE(values,1)) THEN
      localError="The size of the column indices array of "//TRIM(NumberToVString(SIZE(columnIndices,1),"*",err,error))// &
        & " does not conform to the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    DO k=1,SIZE(rowIndices,1)
      CALL Matrix_StorageLocationGet(matrix,rowIndices(k),columnIndices(k),location,err,error,*999)
      matrix%dataIntg(location)=values(k)
    ENDDO !k

    EXITS("Matrix_ValuesSetIntg1")
    RETURN
999 ERRORSEXITS("Matrix_ValuesSetIntg1",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_ValuesSetIntg1

  !
  !================================================================================================================================
  !

  !>Sets the matrix of values in an integer matrix at the location specified by the row and column indices i.e., matrix(i,j)=value
  SUBROUTINE Matrix_ValuesSetIntg2(matrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The row index for the ij'th value to set
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(j). The column index for the ij'th value to set
    INTEGER(INTG), INTENT(IN) :: values(:,:) !<values(i,j). The value of the i,j'th value to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,j,location
    TYPE(VARYING_STRING) :: localError

    ENTERS("Matrix_ValuesSetIntg2",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    CALL Matrix_AssertIsIntgData(matrix,err,error,*999)
    IF(SIZE(rowIndices,1)/=SIZE(values,1)) THEN
      localError="The size of the row indices array of "//TRIM(NumberToVString(SIZE(rowIndices,1),"*",err,error))// &
        & " does not conform to the number of rows in the values array of "// &
        & TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF      
    IF(SIZE(columnIndices,1)/=SIZE(values,2)) THEN
      localError="The size of the column indices array of "//TRIM(NumberToVString(SIZE(columnIndices,1),"*",err,error))// &
        & " does not conform to the number of columns in the values array of "// &
        & TRIM(NumberToVString(SIZE(values,2),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
      
    DO i=1,SIZE(rowIndices,1)
      DO j=1,SIZE(columnIndices,1)
        CALL Matrix_StorageLocationGet(matrix,rowIndices(i),columnIndices(j),location,err,error,*999)
        matrix%dataIntg(location)=values(i,j)
      ENDDO !j
    ENDDO !i
   
    EXITS("Matrix_ValuesSetIntg2")
    RETURN
999 ERRORSEXITS("Matrix_ValuesSetIntg2",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_ValuesSetIntg2
  
  !
  !================================================================================================================================
  !

  !>Sets a value in a single precision real matrix at the location specified by the row and column index i.e., matrix(i,j)=value
  SUBROUTINE Matrix_ValuesSetSP(matrix,rowIndex,columnIndex,VALUE,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: rowIndex !<The row index of the value to set
    INTEGER(INTG), INTENT(IN) :: columnIndex !<The column index of the value to set
    REAL(SP), INTENT(IN) :: value !<The value to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: location

    ENTERS("Matrix_ValuesSetSP",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    CALL Matrix_AssertIsSPData(matrix,err,error,*999)
    
    CALL Matrix_StorageLocationGet(matrix,rowIndex,columnIndex,location,err,error,*999)
    matrix%dataSP(location)=value

    EXITS("Matrix_ValuesSetSP")
    RETURN
999 ERRORSEXITS("Matrix_ValuesSetSP",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_ValuesSetSP

  !
  !================================================================================================================================
  !

  !>Sets the values in a single precision real matrix at the location specified by the row and column indices i.e., matrix(i,j)=value
  SUBROUTINE Matrix_ValuesSetSP1(matrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The row index for the i'th value to set
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(i). The column index for the i'th value to set
    REAL(SP), INTENT(IN) :: values(:) !<values(i). The value of the i'th value to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: k,location
    TYPE(VARYING_STRING) :: localError

    ENTERS("Matrix_ValuesSetSP1",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    CALL Matrix_AssertIsSPData(matrix,err,error,*999)
    IF(SIZE(rowIndices,1)/=SIZE(values,1)) THEN
      localError="The size of the row indices array of "//TRIM(NumberToVString(SIZE(rowIndices,1),"*",err,error))// &
        & " does not conform to the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(columnIndices,1)/=SIZE(values,1)) THEN
      localError="The size of the column indices array of "//TRIM(NumberToVString(SIZE(columnIndices,1),"*",err,error))// &
        & " does not conform to the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    DO k=1,SIZE(rowIndices,1)
      CALL Matrix_StorageLocationGet(matrix,rowIndices(k),columnIndices(k),location,err,error,*999)
      matrix%dataSP(location)=values(k)
    ENDDO !k

    EXITS("Matrix_ValuesSetSP1")
    RETURN
999 ERRORSEXITS("Matrix_ValuesSetSP1",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_ValuesSetSP1

  !
  !================================================================================================================================
  !

  !>Sets the matrix of values in a single precision real matrix at the location specified by the row and column indices i.e., matrix(i,j)=value
  SUBROUTINE Matrix_ValuesSetSP2(matrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The row index for the ij'th value to set
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(j). The column index for the ij'th value to set
    REAL(SP), INTENT(IN) :: values(:,:) !<values(i,j). The value of the i,j'th value to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,j,location
    TYPE(VARYING_STRING) :: localError

    ENTERS("Matrix_ValuesSetSP2",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    CALL Matrix_AssertIsSPData(matrix,err,error,*999)
    IF(SIZE(rowIndices,1)/=SIZE(values,1)) THEN
      localError="The size of the row indices array of "//TRIM(NumberToVString(SIZE(rowIndices,1),"*",err,error))// &
        & " does not conform to the number of rows in the values array of "// &
        & TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF      
    IF(SIZE(columnIndices,1)/=SIZE(values,2)) THEN
      localError="The size of the column indices array of "//TRIM(NumberToVString(SIZE(columnIndices,1),"*",err,error))// &
        & " does not conform to the number of columns in the values array of "// &
        & TRIM(NumberToVString(SIZE(values,2),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
      
    DO i=1,SIZE(rowIndices,1)
      DO j=1,SIZE(columnIndices,1)
        CALL Matrix_StorageLocationGet(matrix,rowIndices(i),columnIndices(j),location,err,error,*999)
        matrix%dataSP(location)=values(i,j)
      ENDDO !j
    ENDDO !i
   
    EXITS("Matrix_ValuesSetSP2")
    RETURN
999 ERRORSEXITS("Matrix_ValuesSetSP2",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_ValuesSetSP2
  
  !
  !================================================================================================================================
  !

  !>Sets a value in a double precision real matrix at the location specified by the row and column index i.e., matrix(i,j)=value
  SUBROUTINE Matrix_ValuesSetDP(matrix,rowIndex,columnIndex,VALUE,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: rowIndex !<The row index of the value to set
    INTEGER(INTG), INTENT(IN) :: columnIndex !<The column index of the value to set
    REAL(DP), INTENT(IN) :: value !<The value to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: location

    ENTERS("Matrix_ValuesSetDP",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    CALL Matrix_AssertIsDPData(matrix,err,error,*999)
    
    CALL Matrix_StorageLocationGet(matrix,rowIndex,columnIndex,location,err,error,*999)
    matrix%dataDP(location)=value

    EXITS("Matrix_ValuesSetDP")
    RETURN
999 ERRORSEXITS("Matrix_ValuesSetDP",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_ValuesSetDP

  !
  !================================================================================================================================
  !

  !>Sets the values in a double precision real matrix at the location specified by the row and column indices i.e., matrix(i,j)=value
  SUBROUTINE Matrix_ValuesSetDP1(matrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The row index for the i'th value to set
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(i). The column index for the i'th value to set
    REAL(DP), INTENT(IN) :: values(:) !<values(i). The value of the i'th value to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: k,location
    TYPE(VARYING_STRING) :: localError

    ENTERS("Matrix_ValuesSetDP1",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    CALL Matrix_AssertIsDPData(matrix,err,error,*999)
    IF(SIZE(rowIndices,1)/=SIZE(values,1)) THEN
      localError="The size of the row indices array of "//TRIM(NumberToVString(SIZE(rowIndices,1),"*",err,error))// &
        & " does not conform to the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(columnIndices,1)/=SIZE(values,1)) THEN
      localError="The size of the column indices array of "//TRIM(NumberToVString(SIZE(columnIndices,1),"*",err,error))// &
        & " does not conform to the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    DO k=1,SIZE(rowIndices,1)
      CALL Matrix_StorageLocationGet(matrix,rowIndices(k),columnIndices(k),location,err,error,*999)
      matrix%dataDP(location)=values(k)
    ENDDO !k

    EXITS("Matrix_ValuesSetDP1")
    RETURN
999 ERRORSEXITS("Matrix_ValuesSetDP1",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_ValuesSetDP1

  !
  !================================================================================================================================
  !

  !>Sets the matrix of values in a double precision real matrix at the location specified by the row and column indices i.e., matrix(i,j)=value
  SUBROUTINE Matrix_ValuesSetDP2(matrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The row index for the ij'th value to set
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(j). The column index for the ij'th value to set
    REAL(DP), INTENT(IN) :: values(:,:) !<values(i,j). The value of the i,j'th value to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,j,location
    TYPE(VARYING_STRING) :: localError

    ENTERS("Matrix_ValuesSetDP2",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    CALL Matrix_AssertIsDPData(matrix,err,error,*999)
    IF(SIZE(rowIndices,1)/=SIZE(values,1)) THEN
      localError="The size of the row indices array of "//TRIM(NumberToVString(SIZE(rowIndices,1),"*",err,error))// &
        & " does not conform to the number of rows in the values array of "// &
        & TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF      
    IF(SIZE(columnIndices,1)/=SIZE(values,2)) THEN
      localError="The size of the column indices array of "//TRIM(NumberToVString(SIZE(columnIndices,1),"*",err,error))// &
        & " does not conform to the number of columns in the values array of "// &
        & TRIM(NumberToVString(SIZE(values,2),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
      
    DO i=1,SIZE(rowIndices,1)
      DO j=1,SIZE(columnIndices,1)
        CALL Matrix_StorageLocationGet(matrix,rowIndices(i),columnIndices(j),location,err,error,*999)
        matrix%dataDP(location)=values(i,j)
      ENDDO !j
    ENDDO !i
   
    EXITS("Matrix_ValuesSetDP2")
    RETURN
999 ERRORSEXITS("Matrix_ValuesSetDP2",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_ValuesSetDP2
  
  !
  !================================================================================================================================
  !

  !>Sets a value in a logical matrix at the location specified by the row and column index i.e., matrix(i,j)=value
  SUBROUTINE Matrix_ValuesSetL(matrix,rowIndex,columnIndex,VALUE,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: rowIndex !<The row index of the value to set
    INTEGER(INTG), INTENT(IN) :: columnIndex !<The column index of the value to set
    LOGICAL, INTENT(IN) :: value !<The value to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: location

    ENTERS("Matrix_ValuesSetL",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    CALL Matrix_AssertIsLData(matrix,err,error,*999)
    
    CALL Matrix_StorageLocationGet(matrix,rowIndex,columnIndex,location,err,error,*999)
    matrix%dataL(location)=value

    EXITS("Matrix_ValuesSetL")
    RETURN
999 ERRORSEXITS("Matrix_ValuesSetL",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_ValuesSetL

  !
  !================================================================================================================================
  !

  !>Sets the values in a logical matrix at the location specified by the row and column indices i.e., matrix(i,j)=value
  SUBROUTINE Matrix_ValuesSetL1(matrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The row index for the i'th value to set
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(i). The column index for the i'th value to set
    LOGICAL, INTENT(IN) :: values(:) !<values(i). The value of the i'th value to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: k,location
    TYPE(VARYING_STRING) :: localError

    ENTERS("Matrix_ValuesSetL1",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    CALL Matrix_AssertIsLData(matrix,err,error,*999)
    IF(SIZE(rowIndices,1)/=SIZE(values,1)) THEN
      localError="The size of the row indices array of "//TRIM(NumberToVString(SIZE(rowIndices,1),"*",err,error))// &
        & " does not conform to the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(columnIndices,1)/=SIZE(values,1)) THEN
      localError="The size of the column indices array of "//TRIM(NumberToVString(SIZE(columnIndices,1),"*",err,error))// &
        & " does not conform to the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    DO k=1,SIZE(rowIndices,1)
      CALL Matrix_StorageLocationGet(matrix,rowIndices(k),columnIndices(k),location,err,error,*999)
      matrix%dataL(location)=values(k)
    ENDDO !k

    EXITS("Matrix_ValuesSetL1")
    RETURN
999 ERRORSEXITS("Matrix_ValuesSetL1",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_ValuesSetL1

  !
  !================================================================================================================================
  !

  !>Sets the matrix of values in a logical matrix at the location specified by the row and column indices i.e., matrix(i,j)=value
  SUBROUTINE Matrix_ValuesSetL2(matrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The row index for the ij'th value to set
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(j). The column index for the ij'th value to set
    LOGICAL, INTENT(IN) :: values(:,:) !<values(i,j). The value of the i,j'th value to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,j,location
    TYPE(VARYING_STRING) :: localError

    ENTERS("Matrix_ValuesSetL2",err,error,*999)

    CALL Matrix_AssertIsFinished(matrix,err,error,*999)
    CALL Matrix_AssertIsLData(matrix,err,error,*999)
    IF(SIZE(rowIndices,1)/=SIZE(values,1)) THEN
      localError="The size of the row indices array of "//TRIM(NumberToVString(SIZE(rowIndices,1),"*",err,error))// &
        & " does not conform to the number of rows in the values array of "// &
        & TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF      
    IF(SIZE(columnIndices,1)/=SIZE(values,2)) THEN
      localError="The size of the column indices array of "//TRIM(NumberToVString(SIZE(columnIndices,1),"*",err,error))// &
        & " does not conform to the number of columns in the values array of "// &
        & TRIM(NumberToVString(SIZE(values,2),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
      
    DO i=1,SIZE(rowIndices,1)
      DO j=1,SIZE(columnIndices,1)
        CALL Matrix_StorageLocationGet(matrix,rowIndices(i),columnIndices(j),location,err,error,*999)
        matrix%dataL(location)=values(i,j)
      ENDDO !j
    ENDDO !i
   
    EXITS("Matrix_ValuesSetL2")
    RETURN
999 ERRORSEXITS("Matrix_ValuesSetL2",err,error)
    RETURN 1
    
  END SUBROUTINE Matrix_ValuesSetL2
  
  !
  !================================================================================================================================
  !

  !>Sets all values in an integer vector to the specified value.
  SUBROUTINE Vector_AllValuesSetIntg(vector,value,err,error,*)

    !Argument variables
    TYPE(VectorType), POINTER :: vector !<A pointer to the vector to set
    INTEGER(INTG), INTENT(IN) :: value !<The value to set the vector to
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables

    ENTERS("Vector_AllValuesSetIntg",err,error,*999)

    CALL Vector_AssertIsFinished(vector,err,error,*999)
    CALL Vector_AssertIsIntgData(vector,err,error,*999)

    vector%dataIntg=value

    EXITS("Vector_AllValuesSetIntg")
    RETURN
999 ERRORSEXITS("Vector_AllValuesSetIntg",err,error)
    RETURN 1
    
  END SUBROUTINE Vector_AllValuesSetIntg

  !
  !================================================================================================================================
  !

  !>Sets all values in a single precision real vector to the specified value.
  SUBROUTINE Vector_AllValuesSetSP(vector,VALUE,err,error,*)

    !Argument variables
    TYPE(VectorType), POINTER :: vector !<A pointer to the vector to set
    REAL(SP), INTENT(IN) :: value !<The value to set the vector to
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables

    ENTERS("Vector_AllValuesSetSP",err,error,*999)

    CALL Vector_AssertIsFinished(vector,err,error,*999)
    CALL Vector_AssertIsSPData(vector,err,error,*999)

    vector%dataSP=value

    EXITS("Vector_AllValuesSetSP")
    RETURN
999 ERRORSEXITS("Vector_AllValuesSetSP",err,error)
    RETURN 1
    
  END SUBROUTINE Vector_AllValuesSetSP

  !
  !================================================================================================================================
  !

  !>Sets all values in a double precision real vector to the specified value.
  SUBROUTINE Vector_AllValuesSetDP(vector,VALUE,err,error,*)

    !Argument variables
    TYPE(VectorType), POINTER :: vector !<A pointer to the vector to set
    REAL(DP), INTENT(IN) :: value !<The value to set the vector to
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables

    ENTERS("Vector_AllValuesSetDP",err,error,*999)

    CALL Vector_AssertIsFinished(vector,err,error,*999)
    CALL Vector_AssertIsDPData(vector,err,error,*999)

    vector%dataDP=value

    EXITS("Vector_AllValuesSetDP")
    RETURN
999 ERRORSEXITS("Vector_AllValuesSetDP",err,error)
    RETURN 1
    
  END SUBROUTINE Vector_AllValuesSetDP

  !
  !================================================================================================================================
  !

  !>Sets all values in a logical vector to the specified value.
  SUBROUTINE Vector_AllValuesSetL(vector,VALUE,err,error,*)

    !Argument variables
    TYPE(VectorType), POINTER :: vector !<A pointer to the vector to set
    LOGICAL, INTENT(IN) :: value !<The value to set the vector to
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables

    ENTERS("Vector_AllValuesSetL",err,error,*999)

    CALL Vector_AssertIsFinished(vector,err,error,*999)
    CALL Vector_AssertIsLData(vector,err,error,*999)

    vector%dataL=value

    EXITS("Vector_AllValuesSetL")
    RETURN
999 ERRORSEXITS("Vector_AllValuesSetL",err,error)
    RETURN 1
    
  END SUBROUTINE Vector_AllValuesSetL

  !
  !================================================================================================================================
  !

  !>Finihses the creation of a vector. 
  SUBROUTINE Vector_CreateFinish(vector,err,error,*)

    !Argument variables
    TYPE(VectorType), POINTER :: vector !<A pointer to the vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Vector_CreateFinish",err,error,*999)

    CALL Vector_AssertNotFinished(vector,err,error,*999)
    IF(vector%size>0) THEN
      SELECT CASE(vector%dataType)
      CASE(MATRIX_VECTOR_INTG_TYPE)
        ALLOCATE(vector%dataIntg(vector%size),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate vector integer data.",err,error,*999)
      CASE(MATRIX_VECTOR_SP_TYPE)
        ALLOCATE(vector%dataSP(vector%size),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate vector single precision data.",err,error,*999)
      CASE(MATRIX_VECTOR_DP_TYPE)
        ALLOCATE(vector%dataDP(vector%size),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate vector double precision data.",err,error,*999)
      CASE(MATRIX_VECTOR_L_TYPE)
        ALLOCATE(vector%dataL(vector%size),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate vector logical data.",err,error,*999)
      CASE(MATRIX_VECTOR_SPC_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(MATRIX_VECTOR_DPC_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The vector data type of "//TRIM(NumberToVString(vector%dataType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDIF
    vector%ID=MATRIX_VECTOR_ID
    MATRIX_VECTOR_ID=MATRIX_VECTOR_ID+1
    vector%vectorFinished=.TRUE.
   
    EXITS("Vector_CreateFinish")
    RETURN
999 ERRORSEXITS("Vector_CreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE Vector_CreateFinish

  !
  !================================================================================================================================
  !

  !>Starts the creation a vector. 
  SUBROUTINE Vector_CreateStart(vector,err,error,*)

    !Argument variables
    TYPE(VectorType), POINTER :: vector !<A pointer to the vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Vector_CreateStart",err,error,*998)

    IF(ASSOCIATED(vector)) CALL FlagError("Vector is already associated.",err,error,*998)
    
    ALLOCATE(vector,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate the vector.",err,error,*999)
    CALL Vector_Initialise(vector,err,error,*999)
    !Set the defaults
    vector%dataType=MATRIX_VECTOR_DP_TYPE
 
    EXITS("Vector_CreateStart")
    RETURN
999 CALL Vector_Finalise(vector,err,error,*998)
998 ERRORSEXITS("Vector_CreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE Vector_CreateStart

  !
  !================================================================================================================================
  !

  !>Sets/changes the data type of a vector.
  SUBROUTINE Vector_DataTypeSet(vector,dataType,err,error,*)

    !Argument variables
    TYPE(VectorType), POINTER :: vector !<A pointer to the vector.
    INTEGER(INTG), INTENT(IN) :: dataType !<The data type to set. \see MatrixVector_DataTypes,MatrixVector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Vector_DataTypeSet",err,error,*999)

    CALL Vector_AssertNotFinished(vector,err,error,*999)

    SELECT CASE(dataType)
    CASE(MATRIX_VECTOR_INTG_TYPE)
      vector%dataType=MATRIX_VECTOR_INTG_TYPE
    CASE(MATRIX_VECTOR_SP_TYPE)
      vector%dataType=MATRIX_VECTOR_SP_TYPE
    CASE(MATRIX_VECTOR_DP_TYPE)
      vector%dataType=MATRIX_VECTOR_DP_TYPE
    CASE(MATRIX_VECTOR_L_TYPE)
      vector%dataType=MATRIX_VECTOR_L_TYPE
    CASE(MATRIX_VECTOR_SPC_TYPE)
      vector%dataType=MATRIX_VECTOR_SPC_TYPE
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(MATRIX_VECTOR_DPC_TYPE)
      vector%dataType=MATRIX_VECTOR_DPC_TYPE
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The vector data type of "//TRIM(NumberToVString(dataType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Vector_DataTypeSet")
    RETURN
999 ERRORSEXITS("Vector_DataTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE Vector_DataTypeSet

  !
  !================================================================================================================================
  !

  !>Destroys a vector
  SUBROUTINE Vector_Destroy(vector,err,error,*)

    !Argument variables
    TYPE(VectorType), POINTER :: vector !<A pointer to the vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Vector_Destroy",err,error,*999)

    IF(.NOT.ASSOCIATED(vector)) CALL FlagError("Vector is not associated.",err,error,*999)
    
    CALL Vector_Finalise(vector,err,error,*999)

    EXITS("Vector_Destroy")
    RETURN
999 ERRORSEXITS("Vector_Destroy",err,error)
    RETURN 1
    
  END SUBROUTINE Vector_Destroy

  !
  !================================================================================================================================
  !

  !>Duplicates a vector structure and returns a pointer to the new vector in newVector.
  SUBROUTINE Vector_Duplicate(vector,newVector,err,error,*)

    !Argument variables
    TYPE(VectorType), POINTER :: vector !<A pointer to the vector to duplicate
    TYPE(VectorType), POINTER :: newVector !<On return a pointer to the new duplicated vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Vector_Duplicate",err,error,*998)

    IF(ASSOCIATED(newVector)) CALL FlagError("New vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vector)) CALL FlagError("Vector is not associated.",err,error,*999)
    
    CALL Vector_CreateStart(newVector,err,error,*999)
    CALL Vector_DataTypeSet(newVector,vector%dataType,err,error,*999)
    CALL Vector_SizeSet(newVector,vector%n,err,error,*999)
    CALL Vector_CreateFinish(newVector,err,error,*999)
 
    EXITS("Vector_Duplicate")
    RETURN
999 CALL Vector_Finalise(newVector,err,error,*998)
998 ERRORSEXITS("Vector_Duplicate",err,error)
    RETURN 1
    
  END SUBROUTINE Vector_Duplicate

  !
  !================================================================================================================================
  !

  !>Finalises a vector and deallocates all memory
  SUBROUTINE Vector_Finalise(vector,err,error,*)

    !Argument variables
    TYPE(VectorType), POINTER :: vector !<A pointer to the vector to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Vector_Finalise",err,error,*999)

    IF(ASSOCIATED(vector)) THEN
      IF(ALLOCATED(vector%dataIntg)) DEALLOCATE(vector%dataIntg)
      IF(ALLOCATED(vector%dataSP)) DEALLOCATE(vector%dataSP)
      IF(ALLOCATED(vector%dataDP)) DEALLOCATE(vector%dataDP)
      IF(ALLOCATED(vector%dataL)) DEALLOCATE(vector%dataL)
      DEALLOCATE(vector)
    ENDIF

    EXITS("Vector_Finalise")
    RETURN
999 ERRORSEXITS("Vector_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE Vector_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises a vector
  SUBROUTINE Vector_Initialise(vector,err,error,*)

    !Argument variables
    TYPE(VectorType), POINTER :: vector !<A pointer to the vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Vector_Initialise",err,error,*999)

    IF(.NOT.ASSOCIATED(vector)) CALL FlagError("Vector is not associated.",err,error,*999)
    
!!TODO: have a vector user number etc.
    vector%ID=0
    vector%vectorFinished=.FALSE.
    vector%n=0
    vector%dataType=0
    vector%size=0      

    EXITS("Vector_Initialise")
    RETURN
999 ERRORSEXITS("Vector_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE Vector_Initialise

  !
  !================================================================================================================================
  !

  !>Sets/changes the size of a vector
  SUBROUTINE Vector_SizeSet(vector,n,err,error,*)

    !Argument variables
    TYPE(VectorType), POINTER :: vector !<A pointer to the vector
    INTEGER(INTG), INTENT(IN) :: n !<The size of the vector to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Vector_SizeSet",err,error,*999)

    CALL Vector_AssertNotFinished(vector,err,error,*999)

    IF(n<=0) THEN
      localError="The size of the vector of "//TRIM(NumberToVString(n,"*",err,error))//" is invalid. The number must be >0."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    vector%n=n
    
    EXITS("Vector_SizeSet")
    RETURN
999 ERRORSEXITS("Vector_SizeSet",err,error)
    RETURN 1
    
  END SUBROUTINE Vector_SizeSet

  !
  !================================================================================================================================
  !

  !>Gets a value in an integer vector at the location specified by the index
  SUBROUTINE Vector_ValuesGetIntg(vector,index,value,err,error,*)

    !Argument variables
    TYPE(VectorType), POINTER :: vector !<A pointer to the vector
    INTEGER(INTG), INTENT(IN) :: index !<The index of the vector to get
    INTEGER(INTG), INTENT(OUT) :: value !<On return the value of the vector at the specified index
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Vector_ValuesGetIntg",err,error,*999)

    CALL Vector_AssertIsFinished(vector,err,error,*999)
    CALL Vector_AssertIsIntgData(vector,err,error,*999)
    IF(index<1.OR.index>vector%n) THEN
      localError="The index value of "//TRIM(NumberToVString(index,"*",err,error))// &
        & " is invalid. The index must be between 1 and "//TRIM(NumberToVString(vector%n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    value=vector%dataIntg(index)

    EXITS("Vector_ValuesGetIntg")
    RETURN
999 ERRORSEXITS("Vector_ValuesGetIntg",err,error)
    RETURN 1
    
  END SUBROUTINE Vector_ValuesGetIntg

  !
  !================================================================================================================================
  !

  !>Gets the values in an integer vector at the indices specified.
  SUBROUTINE Vector_ValuesGetIntg1(vector,indices,values,err,error,*)

    !Argument variables
    TYPE(VectorType), POINTER :: vector !<A pointer to the vector
    INTEGER(INTG), INTENT(IN) :: indices(:) !<indices(i). The i'th index to get
    INTEGER(INTG), INTENT(OUT) :: values(:) !<values(i). On return the i'th value to get
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,k
    TYPE(VARYING_STRING) :: localError

    ENTERS("Vector_ValuesGetIntg1",err,error,*999)

    CALL Vector_AssertIsFinished(vector,err,error,*999)
    CALL Vector_AssertIsIntgData(vector,err,error,*999)
    IF(SIZE(indices,1)/=SIZE(values,1)) THEN
      localError="The size of the indices array of "//TRIM(NumberToVString(SIZE(INDICES,1),"*",err,error))// &
        & " does not conform to the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    DO i=1,SIZE(indices,1)
      k=indices(i)
      IF(k<1.OR.k>vector%n) THEN
        localError="Index number "//TRIM(NumberToVString(i,"*",err,error))//" is invalid. The index is "// &
          & TRIM(NumberToVString(k,"*",err,error))//" and it must be between 1 and "// &
          & TRIM(NumberToVString(vector%n,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      values(i)=vector%dataIntg(k)
    ENDDO !i

    EXITS("Vector_ValuesGetIntg1")
    RETURN
999 ERRORSEXITS("Vector_ValuesGetIntg1",err,error)
    RETURN 1
    
  END SUBROUTINE Vector_ValuesGetIntg1

  !
  !================================================================================================================================
  !

  !>Gets a value in a single precision real vector at the location specified by the index
  SUBROUTINE Vector_ValuesGetSP(vector,index,VALUE,err,error,*)

    !Argument variables
    TYPE(VectorType), POINTER :: vector !<A pointer to the vector
    INTEGER(INTG), INTENT(IN) :: index !<The index of the vector to get
    REAL(SP), INTENT(OUT) :: value !<On return the value of the vector at the specified index
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Vector_ValuesGetSP",err,error,*999)

    CALL Vector_AssertIsFinished(vector,err,error,*999)
    CALL Vector_AssertIsSPData(vector,err,error,*999)
    IF(index<1.OR.index>vector%n) THEN
      localError="The index value of "//TRIM(NumberToVString(index,"*",err,error))// &
        & " is invalid. The index must be between 1 and "//TRIM(NumberToVString(vector%n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    value=vector%dataSP(index)

    EXITS("Vector_ValuesGetSP")
    RETURN
999 ERRORSEXITS("Vector_ValuesGetSP",err,error)
    RETURN 1
    
  END SUBROUTINE Vector_ValuesGetSP

  !
  !================================================================================================================================
  !

  !>Gets the values in a single precision real vector at the indices specified.
  SUBROUTINE Vector_ValuesGetSP1(vector,indices,values,err,error,*)

    !Argument variables
    TYPE(VectorType), POINTER :: vector !<A pointer to the vector
    INTEGER(INTG), INTENT(IN) :: indices(:) !<indices(i). The i'th index to get
    REAL(SP), INTENT(OUT) :: values(:) !<values(i). On return the i'th value to get
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,k
    TYPE(VARYING_STRING) :: localError

    ENTERS("Vector_ValuesGetSP1",err,error,*999)

    CALL Vector_AssertIsFinished(vector,err,error,*999)
    CALL Vector_AssertIsSPData(vector,err,error,*999)
    IF(SIZE(indices,1)/=SIZE(values,1)) THEN
      localError="The size of the indices array of "//TRIM(NumberToVString(SIZE(INDICES,1),"*",err,error))// &
        & " does not conform to the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    DO i=1,SIZE(indices,1)
      k=indices(i)
      IF(k<1.OR.k>vector%n) THEN
        localError="Index number "//TRIM(NumberToVString(i,"*",err,error))//" is invalid. The index is "// &
          & TRIM(NumberToVString(k,"*",err,error))//" and it must be between 1 and "// &
          & TRIM(NumberToVString(vector%n,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      values(i)=vector%dataSP(k)
    ENDDO !i

    EXITS("Vector_ValuesGetSP1")
    RETURN
999 ERRORSEXITS("Vector_ValuesGetSP1",err,error)
    RETURN 1
    
  END SUBROUTINE Vector_ValuesGetSP1

  !
  !================================================================================================================================
  !

  !>Gets a value in a double precision real vector at the location specified by the index
  SUBROUTINE Vector_ValuesGetDP(vector,index,VALUE,err,error,*)

    !Argument variables
    TYPE(VectorType), POINTER :: vector !<A pointer to the vector
    INTEGER(INTG), INTENT(IN) :: index !<The index of the vector to get
    REAL(DP), INTENT(OUT) :: value !<On return the value of the vector at the specified index
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Vector_ValuesGetDP",err,error,*999)

    CALL Vector_AssertIsFinished(vector,err,error,*999)
    CALL Vector_AssertIsDPData(vector,err,error,*999)
    IF(index<1.OR.index>vector%n) THEN
      localError="The index value of "//TRIM(NumberToVString(index,"*",err,error))// &
        & " is invalid. The index must be between 1 and "//TRIM(NumberToVString(vector%n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    value=vector%dataDP(index)

    EXITS("Vector_ValuesGetDP")
    RETURN
999 ERRORSEXITS("Vector_ValuesGetDP",err,error)
    RETURN 1
    
  END SUBROUTINE Vector_ValuesGetDP

  !
  !================================================================================================================================
  !

  !>Gets the values in a double precision real vector at the indices specified.
  SUBROUTINE Vector_ValuesGetDP1(vector,indices,values,err,error,*)

    !Argument variables
    TYPE(VectorType), POINTER :: vector !<A pointer to the vector
    INTEGER(INTG), INTENT(IN) :: indices(:) !<indices(i). The i'th index to get
    REAL(DP), INTENT(OUT) :: values(:) !<values(i). On return the i'th value to get
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,k
    TYPE(VARYING_STRING) :: localError

    ENTERS("Vector_ValuesGetDP1",err,error,*999)

    CALL Vector_AssertIsFinished(vector,err,error,*999)
    CALL Vector_AssertIsDPData(vector,err,error,*999)
    IF(SIZE(indices,1)/=SIZE(values,1)) THEN
      localError="The size of the indices array of "//TRIM(NumberToVString(SIZE(INDICES,1),"*",err,error))// &
        & " does not conform to the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    DO i=1,SIZE(indices,1)
      k=indices(i)
      IF(k<1.OR.k>vector%n) THEN
        localError="Index number "//TRIM(NumberToVString(i,"*",err,error))//" is invalid. The index is "// &
          & TRIM(NumberToVString(k,"*",err,error))//" and it must be between 1 and "// &
          & TRIM(NumberToVString(vector%n,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      values(i)=vector%dataDP(k)
    ENDDO !i

    EXITS("Vector_ValuesGetDP1")
    RETURN
999 ERRORSEXITS("Vector_ValuesGetDP1",err,error)
    RETURN 1
    
  END SUBROUTINE Vector_ValuesGetDP1

  !
  !================================================================================================================================
  !

  !>Gets a value in a logical vector at the location specified by the index
  SUBROUTINE Vector_ValuesGetL(vector,index,VALUE,err,error,*)

    !Argument variables
    TYPE(VectorType), POINTER :: vector !<A pointer to the vector
    INTEGER(INTG), INTENT(IN) :: index !<The index of the vector to get
    LOGICAL, INTENT(OUT) :: value !<On return the value of the vector at the specified index
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Vector_ValuesGetL",err,error,*999)

    CALL Vector_AssertIsFinished(vector,err,error,*999)
    CALL Vector_AssertIsLData(vector,err,error,*999)
    IF(index<1.OR.index>vector%n) THEN
      localError="The index value of "//TRIM(NumberToVString(index,"*",err,error))// &
        & " is invalid. The index must be between 1 and "//TRIM(NumberToVString(vector%n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    value=vector%dataL(index)

    EXITS("Vector_ValuesGetL")
    RETURN
999 ERRORSEXITS("Vector_ValuesGetL",err,error)
    RETURN 1
    
  END SUBROUTINE Vector_ValuesGetL

  !
  !================================================================================================================================
  !

  !>Gets the values in a logical vector at the indices specified.
  SUBROUTINE Vector_ValuesGetL1(vector,indices,values,err,error,*)

    !Argument variables
    TYPE(VectorType), POINTER :: vector !<A pointer to the vector
    INTEGER(INTG), INTENT(IN) :: indices(:) !<indices(i). The i'th index to get
    LOGICAL, INTENT(OUT) :: values(:) !<values(i). On return the i'th value to get
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,k
    TYPE(VARYING_STRING) :: localError

    ENTERS("Vector_ValuesGetL1",err,error,*999)

    CALL Vector_AssertIsFinished(vector,err,error,*999)
    CALL Vector_AssertIsLData(vector,err,error,*999)
    IF(SIZE(indices,1)/=SIZE(values,1)) THEN
      localError="The size of the indices array of "//TRIM(NumberToVString(SIZE(INDICES,1),"*",err,error))// &
        & " does not conform to the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    DO i=1,SIZE(indices,1)
      k=indices(i)
      IF(k<1.OR.k>vector%n) THEN
        localError="Index number "//TRIM(NumberToVString(i,"*",err,error))//" is invalid. The index is "// &
          & TRIM(NumberToVString(k,"*",err,error))//" and it must be between 1 and "// &
          & TRIM(NumberToVString(vector%n,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      values(i)=vector%dataL(k)
    ENDDO !i

    EXITS("Vector_ValuesGetL1")
    RETURN
999 ERRORSEXITS("Vector_ValuesGetL1",err,error)
    RETURN 1
    
  END SUBROUTINE Vector_ValuesGetL1

  !
  !================================================================================================================================
  !
  
  !>Sets a value in an integer vector at the specified index.
  SUBROUTINE Vector_ValuesSetIntg(vector,index,value,err,error,*)

    !Argument variables
    TYPE(VectorType), POINTER :: vector !<A pointer to the vector
    INTEGER(INTG), INTENT(IN) :: index !<The index to set
    INTEGER(INTG), INTENT(IN) :: value !<The value to set at the specified index
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Vector_ValuesSetIntg",err,error,*999)

    CALL Vector_AssertIsFinished(vector,err,error,*999)
    CALL Vector_AssertIsIntgData(vector,err,error,*999)
    IF(index<1.OR.index>vector%n) THEN
      localError="The index value of "//TRIM(NumberToVString(index,"*",err,error))// &
        & " is invalid. The index must be between 1 and  "//TRIM(NumberToVString(vector%n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    vector%dataIntg(index)=value
    
    EXITS("Vector_ValuesSetIntg")
    RETURN
999 ERRORSEXITS("Vector_ValuesSetIntg",err,error)
    RETURN 1
    
  END SUBROUTINE Vector_ValuesSetIntg

  !
  !================================================================================================================================
  !

  !>Sets the values in an integer vector at the specified indices.
  SUBROUTINE Vector_ValuesSetIntg1(vector,indices,values,err,error,*)

    !Argument variables
    TYPE(VectorType), POINTER :: vector !<A pointer to the vector
    INTEGER(INTG), INTENT(IN) :: indices(:) !<indices(i). The i'th index to set
    INTEGER(INTG), INTENT(IN) :: values(:) !<values(i). The i'th value to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,k
    TYPE(VARYING_STRING) :: localError

    ENTERS("Vector_ValuesSetIntg1",err,error,*999)

    CALL Vector_AssertIsFinished(vector,err,error,*999)
    CALL Vector_AssertIsIntgData(vector,err,error,*999)
    IF(SIZE(indices,1)/=SIZE(values,1)) THEN
      localError="The size of the indices array of "//TRIM(NumberToVString(SIZE(indices,1),"*",err,error))// &
        & " does not conform to the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    DO i=1,SIZE(indices,1)      
      k=indices(i)
      IF(k<1.OR.k>vector%n) THEN
        localError="Index number "//TRIM(NumberToVString(i,"*",err,error))//" is invalid. The index is "// &
          & TRIM(NumberToVString(k,"*",err,error))//" and it must be between 1 and "// &
          & TRIM(NumberToVString(vector%n,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ELSE
        vector%dataIntg(k)=values(i)
      ENDIF
    ENDDO !i

    EXITS("Vector_ValuesSetIntg1")
    RETURN
999 ERRORSEXITS("Vector_ValuesSetIntg1",err,error)
    RETURN 1
    
  END SUBROUTINE Vector_ValuesSetIntg1

  !
  !================================================================================================================================
  !
  
  !>Sets a value in a single precision real vector at the specified index.
  SUBROUTINE Vector_ValuesSetSP(vector,index,value,err,error,*)

    !Argument variables
    TYPE(VectorType), POINTER :: vector !<A pointer to the vector
    INTEGER(INTG), INTENT(IN) :: index !<The index to set
    REAL(SP), INTENT(IN) :: value !<The value to set at the specified index
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Vector_ValuesSetSP",err,error,*999)

    CALL Vector_AssertIsFinished(vector,err,error,*999)
    CALL Vector_AssertIsSPData(vector,err,error,*999)
    IF(index<1.OR.index>vector%n) THEN
      localError="The index value of "//TRIM(NumberToVString(index,"*",err,error))// &
        & " is invalid. The index must be between 1 and  "//TRIM(NumberToVString(vector%n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    vector%dataSP(index)=value
    
    EXITS("Vector_ValuesSetSP")
    RETURN
999 ERRORSEXITS("Vector_ValuesSetSP",err,error)
    RETURN 1
    
  END SUBROUTINE Vector_ValuesSetSP

  !
  !================================================================================================================================
  !

  !>Sets the values in a single precision real vector at the specified indices.
  SUBROUTINE Vector_ValuesSetSP1(vector,indices,values,err,error,*)

    !Argument variables
    TYPE(VectorType), POINTER :: vector !<A pointer to the vector
    INTEGER(INTG), INTENT(IN) :: indices(:) !<indices(i). The i'th index to set
    REAL(SP), INTENT(IN) :: values(:) !<values(i). The i'th value to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,k
    TYPE(VARYING_STRING) :: localError

    ENTERS("Vector_ValuesSetSP1",err,error,*999)

    CALL Vector_AssertIsFinished(vector,err,error,*999)
    CALL Vector_AssertIsSPData(vector,err,error,*999)
    IF(SIZE(indices,1)/=SIZE(values,1)) THEN
      localError="The size of the indices array of "//TRIM(NumberToVString(SIZE(indices,1),"*",err,error))// &
        & " does not conform to the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    DO i=1,SIZE(indices,1)      
      k=indices(i)
      IF(k<1.OR.k>vector%n) THEN
        localError="Index number "//TRIM(NumberToVString(i,"*",err,error))//" is invalid. The index is "// &
          & TRIM(NumberToVString(k,"*",err,error))//" and it must be between 1 and "// &
          & TRIM(NumberToVString(vector%n,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ELSE
        vector%dataSP(k)=values(i)
      ENDIF
    ENDDO !i

    EXITS("Vector_ValuesSetSP1")
    RETURN
999 ERRORSEXITS("Vector_ValuesSetSP1",err,error)
    RETURN 1
    
  END SUBROUTINE Vector_ValuesSetSP1

  !
  !================================================================================================================================
  !
  
  !>Sets a value in a double precision real vector at the specified index.
  SUBROUTINE Vector_ValuesSetDP(vector,index,value,err,error,*)

    !Argument variables
    TYPE(VectorType), POINTER :: vector !<A pointer to the vector
    INTEGER(INTG), INTENT(IN) :: index !<The index to set
    REAL(DP), INTENT(IN) :: value !<The value to set at the specified index
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Vector_ValuesSetDP",err,error,*999)

    CALL Vector_AssertIsFinished(vector,err,error,*999)
    CALL Vector_AssertIsDPData(vector,err,error,*999)
    IF(index<1.OR.index>vector%n) THEN
      localError="The index value of "//TRIM(NumberToVString(index,"*",err,error))// &
        & " is invalid. The index must be between 1 and  "//TRIM(NumberToVString(vector%n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    vector%dataDP(index)=value
    
    EXITS("Vector_ValuesSetDP")
    RETURN
999 ERRORSEXITS("Vector_ValuesSetDP",err,error)
    RETURN 1
    
  END SUBROUTINE Vector_ValuesSetDP

  !
  !================================================================================================================================
  !

  !>Sets the values in a double precision real vector at the specified indices.
  SUBROUTINE Vector_ValuesSetDP1(vector,indices,values,err,error,*)

    !Argument variables
    TYPE(VectorType), POINTER :: vector !<A pointer to the vector
    INTEGER(INTG), INTENT(IN) :: indices(:) !<indices(i). The i'th index to set
    REAL(DP), INTENT(IN) :: values(:) !<values(i). The i'th value to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,k
    TYPE(VARYING_STRING) :: localError

    ENTERS("Vector_ValuesSetDP1",err,error,*999)

    CALL Vector_AssertIsFinished(vector,err,error,*999)
    CALL Vector_AssertIsDPData(vector,err,error,*999)
    IF(SIZE(indices,1)/=SIZE(values,1)) THEN
      localError="The size of the indices array of "//TRIM(NumberToVString(SIZE(indices,1),"*",err,error))// &
        & " does not conform to the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    DO i=1,SIZE(indices,1)      
      k=indices(i)
      IF(k<1.OR.k>vector%n) THEN
        localError="Index number "//TRIM(NumberToVString(i,"*",err,error))//" is invalid. The index is "// &
          & TRIM(NumberToVString(k,"*",err,error))//" and it must be between 1 and "// &
          & TRIM(NumberToVString(vector%n,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ELSE
        vector%dataDP(k)=values(i)
      ENDIF
    ENDDO !i

    EXITS("Vector_ValuesSetDP1")
    RETURN
999 ERRORSEXITS("Vector_ValuesSetDP1",err,error)
    RETURN 1
    
  END SUBROUTINE Vector_ValuesSetDP1

  !
  !================================================================================================================================
  !
  
  !>Sets a value in a logical vector at the specified index.
  SUBROUTINE Vector_ValuesSetL(vector,index,value,err,error,*)

    !Argument variables
    TYPE(VectorType), POINTER :: vector !<A pointer to the vector
    INTEGER(INTG), INTENT(IN) :: index !<The index to set
    LOGICAL, INTENT(IN) :: value !<The value to set at the specified index
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Vector_ValuesSetL",err,error,*999)

    CALL Vector_AssertIsFinished(vector,err,error,*999)
    CALL Vector_AssertIsLData(vector,err,error,*999)
    IF(index<1.OR.index>vector%n) THEN
      localError="The index value of "//TRIM(NumberToVString(index,"*",err,error))// &
        & " is invalid. The index must be between 1 and  "//TRIM(NumberToVString(vector%n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    vector%dataL(index)=value
    
    EXITS("Vector_ValuesSetL")
    RETURN
999 ERRORSEXITS("Vector_ValuesSetL",err,error)
    RETURN 1
    
  END SUBROUTINE Vector_ValuesSetL

  !
  !================================================================================================================================
  !

  !>Sets the values in a logical vector at the specified indices.
  SUBROUTINE Vector_ValuesSetL1(vector,indices,values,err,error,*)

    !Argument variables
    TYPE(VectorType), POINTER :: vector !<A pointer to the vector
    INTEGER(INTG), INTENT(IN) :: indices(:) !<indices(i). The i'th index to set
    LOGICAL, INTENT(IN) :: values(:) !<values(i). The i'th value to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,k
    TYPE(VARYING_STRING) :: localError

    ENTERS("Vector_ValuesSetL1",err,error,*999)

    CALL Vector_AssertIsFinished(vector,err,error,*999)
    CALL Vector_AssertIsLData(vector,err,error,*999)
    IF(SIZE(indices,1)/=SIZE(values,1)) THEN
      localError="The size of the indices array of "//TRIM(NumberToVString(SIZE(indices,1),"*",err,error))// &
        & " does not conform to the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    DO i=1,SIZE(indices,1)      
      k=indices(i)
      IF(k<1.OR.k>vector%n) THEN
        localError="Index number "//TRIM(NumberToVString(i,"*",err,error))//" is invalid. The index is "// &
          & TRIM(NumberToVString(k,"*",err,error))//" and it must be between 1 and "// &
          & TRIM(NumberToVString(vector%n,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ELSE
        vector%dataL(k)=values(i)
      ENDIF
    ENDDO !i

    EXITS("Vector_ValuesSetL1")
    RETURN
999 ERRORSEXITS("Vector_ValuesSetL1",err,error)
    RETURN 1
    
  END SUBROUTINE Vector_ValuesSetL1

  !
  !================================================================================================================================
  !
  
END MODULE MatrixVector
