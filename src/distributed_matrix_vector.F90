!> \file
!> \author Chris Bradley
!> \brief This module handles all distributed matrix vector routines.
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

!> This module handles all distributed matrix vector routines.
MODULE DistributedMatrixVector

  USE BaseRoutines
  USE CmissMPI
  USE CmissPetsc
  USE ComputationEnvironment
  USE DistributedMatrixVectorAccessRoutines
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE ISO_C_BINDING
  USE Kinds
  USE Maths
  USE MatrixVector
#ifndef NOMPIMOD
  USE MPI
#endif
  USE Strings
  USE Types
  USE LINKEDLIST_ROUTINES

#include "macros.h"

  IMPLICIT NONE

  PRIVATE

#ifdef NOMPIMOD
#include "mpif.h"
#endif
#include "petscversion.h"
  
  !Module parameters

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
  
  !Module types

  !Module variables

  INTEGER(INTG), SAVE :: distributedDataId=100000000

  !Interfaces

  INTERFACE DISTRIBUTED_MATRIX_ALL_VALUES_SET
    MODULE PROCEDURE DistributedMatrix_AllValuesSetIntg
    MODULE PROCEDURE DistributedMatrix_AllValuesSetSP
    MODULE PROCEDURE DistributedMatrix_AllValuesSetDP
    MODULE PROCEDURE DistributedMatrix_AllValuesSetL
  END INTERFACE DISTRIBUTED_MATRIX_ALL_VALUES_SET

  INTERFACE DistributedMatrix_AllValuesSet
    MODULE PROCEDURE DistributedMatrix_AllValuesSetIntg
    MODULE PROCEDURE DistributedMatrix_AllValuesSetSP
    MODULE PROCEDURE DistributedMatrix_AllValuesSetDP
    MODULE PROCEDURE DistributedMatrix_AllValuesSetL
  END INTERFACE DistributedMatrix_AllValuesSet

  INTERFACE DISTRIBUTED_MATRIX_CREATE_FINISH
    MODULE PROCEDURE DistributedMatrix_CreateFinish
  END INTERFACE DISTRIBUTED_MATRIX_CREATE_FINISH

  INTERFACE DISTRIBUTED_MATRIX_CREATE_START
    MODULE PROCEDURE DistributedMatrix_CreateStart
  END INTERFACE DISTRIBUTED_MATRIX_CREATE_START

  INTERFACE DISTRIBUTED_MATRIX_DATA_GET
    MODULE PROCEDURE DistributedMatrix_DataGetIntg
    MODULE PROCEDURE DistributedMatrix_DataGetSP
    MODULE PROCEDURE DistributedMatrix_DataGetDP
    MODULE PROCEDURE DistributedMatrix_DataGetL
  END INTERFACE DISTRIBUTED_MATRIX_DATA_GET

  INTERFACE DistributedMatrix_DataGet
    MODULE PROCEDURE DistributedMatrix_DataGetIntg
    MODULE PROCEDURE DistributedMatrix_DataGetSP
    MODULE PROCEDURE DistributedMatrix_DataGetDP
    MODULE PROCEDURE DistributedMatrix_DataGetL
  END INTERFACE DistributedMatrix_DataGet

  INTERFACE DISTRIBUTED_MATRIX_DATA_RESTORE
    MODULE PROCEDURE DistributedMatrix_DataRestoreIntg
    MODULE PROCEDURE DistributedMatrix_DataRestoreSP
    MODULE PROCEDURE DistributedMatrix_DataRestoreDP
    MODULE PROCEDURE DistributedMatrix_DataRestoreL
  END INTERFACE DISTRIBUTED_MATRIX_DATA_RESTORE

  INTERFACE DistributedMatrix_DataRestore
    MODULE PROCEDURE DistributedMatrix_DataRestoreIntg
    MODULE PROCEDURE DistributedMatrix_DataRestoreSP
    MODULE PROCEDURE DistributedMatrix_DataRestoreDP
    MODULE PROCEDURE DistributedMatrix_DataRestoreL
  END INTERFACE DistributedMatrix_DataRestore

  INTERFACE DISTRIBUTED_MATRIX_DATA_TYPE_SET
    MODULE PROCEDURE DistributedMatrix_DataTypeSet
  END INTERFACE DISTRIBUTED_MATRIX_DATA_TYPE_SET

  INTERFACE DISTRIBUTED_MATRIX_DESTROY
    MODULE PROCEDURE DistributedMatrix_Destroy
  END INTERFACE DISTRIBUTED_MATRIX_DESTROY

  INTERFACE DISTRIBUTED_MATRIX_DUPLICATE
    MODULE PROCEDURE DistributedMatrix_Duplicate
  END INTERFACE DISTRIBUTED_MATRIX_DUPLICATE

  INTERFACE DISTRIBUTED_MATRIX_FORM
    MODULE PROCEDURE DistributedMatrix_Form
  END INTERFACE DISTRIBUTED_MATRIX_FORM

  INTERFACE DISTRIBUTED_MATRIX_GHOSTING_TYPE_SET
    MODULE PROCEDURE DistributedMatrix_GhostingTypeSet
  END INTERFACE DISTRIBUTED_MATRIX_GHOSTING_TYPE_SET

  INTERFACE DISTRIBUTED_MATRIX_LIBRARY_TYPE_SET
    MODULE PROCEDURE DistributedMatrix_LibraryTypeSet
  END INTERFACE DISTRIBUTED_MATRIX_LIBRARY_TYPE_SET

  INTERFACE DISTRIBUTED_MATRIX_MAX_COLUMNS_PER_ROW_GET
    MODULE PROCEDURE DistributedMatrix_MaxColumnsPerRowGet
  END INTERFACE DISTRIBUTED_MATRIX_MAX_COLUMNS_PER_ROW_GET

  INTERFACE DISTRIBUTED_MATRIX_NUMBER_NON_ZEROS_GET
    MODULE PROCEDURE DistributedMatrix_NumberOfNonZerosGet
  END INTERFACE DISTRIBUTED_MATRIX_NUMBER_NON_ZEROS_GET

  INTERFACE DISTRIBUTED_MATRIX_NUMBER_NON_ZEROS_SET
    MODULE PROCEDURE DistributedMatrix_NumberOfNonZerosSet
  END INTERFACE DISTRIBUTED_MATRIX_NUMBER_NON_ZEROS_SET

  INTERFACE DISTRIBUTED_MATRIX_OUTPUT
    MODULE PROCEDURE DistributedMatrix_Output
  END INTERFACE DISTRIBUTED_MATRIX_OUTPUT

  INTERFACE DISTRIBUTED_MATRIX_OVERRIDE_SET_ON
    MODULE PROCEDURE DistributedMatrix_OverrideSetOn
  END INTERFACE DISTRIBUTED_MATRIX_OVERRIDE_SET_ON

  INTERFACE DISTRIBUTED_MATRIX_OVERRIDE_SET_OFF
    MODULE PROCEDURE DistributedMatrix_OverrideSetOff
  END INTERFACE DISTRIBUTED_MATRIX_OVERRIDE_SET_OFF

  INTERFACE DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET
    MODULE PROCEDURE DistributedMatrix_StorageLocationsGet
  END INTERFACE DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET

  INTERFACE DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_SET
    MODULE PROCEDURE DistributedMatrix_StorageLocationsSet
  END INTERFACE DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_SET

  INTERFACE DISTRIBUTED_MATRIX_STORAGE_TYPE_GET
    MODULE PROCEDURE DistributedMatrix_StorageTypeGet
  END INTERFACE DISTRIBUTED_MATRIX_STORAGE_TYPE_GET

  INTERFACE DISTRIBUTED_MATRIX_STORAGE_TYPE_SET
    MODULE PROCEDURE DistributedMatrix_StorageTypeSet
  END INTERFACE DISTRIBUTED_MATRIX_STORAGE_TYPE_SET

  INTERFACE DISTRIBUTED_MATRIX_UPDATE_FINISH
    MODULE PROCEDURE DistributedMatrix_UpdateFinish
  END INTERFACE DISTRIBUTED_MATRIX_UPDATE_FINISH

  INTERFACE DISTRIBUTED_MATRIX_UPDATE_START
    MODULE PROCEDURE DistributedMatrix_UpdateStart
  END INTERFACE DISTRIBUTED_MATRIX_UPDATE_START

  INTERFACE DISTRIBUTED_MATRIX_UPDATE_ISFINISHED
    MODULE PROCEDURE DistributedMatrix_UpdateIsFinished
  END INTERFACE DISTRIBUTED_MATRIX_UPDATE_ISFINISHED

  INTERFACE DISTRIBUTED_MATRIX_UPDATE_WAITFINISHED
    MODULE PROCEDURE DistributedMatrix_UpdateWaitFinished
  END INTERFACE DISTRIBUTED_MATRIX_UPDATE_WAITFINISHED

  INTERFACE DISTRIBUTED_MATRIX_VALUES_ADD
    MODULE PROCEDURE DistributedMatrix_ValuesAddIntg
    MODULE PROCEDURE DistributedMatrix_ValuesAddIntg1    
    MODULE PROCEDURE DistributedMatrix_ValuesAddIntg2
    MODULE PROCEDURE DistributedMatrix_ValuesAddSP
    MODULE PROCEDURE DistributedMatrix_ValuesAddSP1    
    MODULE PROCEDURE DistributedMatrix_ValuesAddSP2
    MODULE PROCEDURE DistributedMatrix_ValuesAddDP
    MODULE PROCEDURE DistributedMatrix_ValuesAddDP1    
    MODULE PROCEDURE DistributedMatrix_ValuesAddDP2
    MODULE PROCEDURE DistributedMatrix_ValuesAddL
    MODULE PROCEDURE DistributedMatrix_ValuesAddL1    
    MODULE PROCEDURE DistributedMatrix_ValuesAddL2
  END INTERFACE DISTRIBUTED_MATRIX_VALUES_ADD

  INTERFACE DistributedMatrix_ValuesAdd
    MODULE PROCEDURE DistributedMatrix_ValuesAddIntg
    MODULE PROCEDURE DistributedMatrix_ValuesAddIntg1    
    MODULE PROCEDURE DistributedMatrix_ValuesAddIntg2
    MODULE PROCEDURE DistributedMatrix_ValuesAddSP
    MODULE PROCEDURE DistributedMatrix_ValuesAddSP1    
    MODULE PROCEDURE DistributedMatrix_ValuesAddSP2
    MODULE PROCEDURE DistributedMatrix_ValuesAddDP
    MODULE PROCEDURE DistributedMatrix_ValuesAddDP1    
    MODULE PROCEDURE DistributedMatrix_ValuesAddDP2
    MODULE PROCEDURE DistributedMatrix_ValuesAddL
    MODULE PROCEDURE DistributedMatrix_ValuesAddL1    
    MODULE PROCEDURE DistributedMatrix_ValuesAddL2
  END INTERFACE DistributedMatrix_ValuesAdd

  INTERFACE DISTRIBUTED_MATRIX_VALUES_GET
    MODULE PROCEDURE DistributedMatrix_ValuesGetIntg
    MODULE PROCEDURE DistributedMatrix_ValuesGetIntg1    
    MODULE PROCEDURE DistributedMatrix_ValuesGetIntg2
    MODULE PROCEDURE DistributedMatrix_ValuesGetSP
    MODULE PROCEDURE DistributedMatrix_ValuesGetSP1    
    MODULE PROCEDURE DistributedMatrix_ValuesGetSP2
    MODULE PROCEDURE DistributedMatrix_ValuesGetDP
    MODULE PROCEDURE DistributedMatrix_ValuesGetDP1    
    MODULE PROCEDURE DistributedMatrix_ValuesGetDP2
    MODULE PROCEDURE DistributedMatrix_ValuesGetL
    MODULE PROCEDURE DistributedMatrix_ValuesGetL1    
    MODULE PROCEDURE DistributedMatrix_ValuesGetL2
  END INTERFACE DISTRIBUTED_MATRIX_VALUES_GET

  INTERFACE DistributedMatrix_ValuesGet
    MODULE PROCEDURE DistributedMatrix_ValuesGetIntg
    MODULE PROCEDURE DistributedMatrix_ValuesGetIntg1    
    MODULE PROCEDURE DistributedMatrix_ValuesGetIntg2
    MODULE PROCEDURE DistributedMatrix_ValuesGetSP
    MODULE PROCEDURE DistributedMatrix_ValuesGetSP1    
    MODULE PROCEDURE DistributedMatrix_ValuesGetSP2
    MODULE PROCEDURE DistributedMatrix_ValuesGetDP
    MODULE PROCEDURE DistributedMatrix_ValuesGetDP1    
    MODULE PROCEDURE DistributedMatrix_ValuesGetDP2
    MODULE PROCEDURE DistributedMatrix_ValuesGetL
    MODULE PROCEDURE DistributedMatrix_ValuesGetL1    
    MODULE PROCEDURE DistributedMatrix_ValuesGetL2
  END INTERFACE DistributedMatrix_ValuesGet

  INTERFACE DISTRIBUTED_MATRIX_VALUES_SET
    MODULE PROCEDURE DistributedMatrix_ValuesSetIntg
    MODULE PROCEDURE DistributedMatrix_ValuesSetIntg1    
    MODULE PROCEDURE DistributedMatrix_ValuesSetIntg2
    MODULE PROCEDURE DistributedMatrix_ValuesSetSP
    MODULE PROCEDURE DistributedMatrix_ValuesSetSP1    
    MODULE PROCEDURE DistributedMatrix_ValuesSetSP2
    MODULE PROCEDURE DistributedMatrix_ValuesSetDP
    MODULE PROCEDURE DistributedMatrix_ValuesSetDP1    
    MODULE PROCEDURE DistributedMatrix_ValuesSetDP2
    MODULE PROCEDURE DistributedMatrix_ValuesSetL
    MODULE PROCEDURE DistributedMatrix_ValuesSetL1    
    MODULE PROCEDURE DistributedMatrix_ValuesSetL2
  END INTERFACE DISTRIBUTED_MATRIX_VALUES_SET

  INTERFACE DistributedMatrix_ValuesSet
    MODULE PROCEDURE DistributedMatrix_ValuesSetIntg
    MODULE PROCEDURE DistributedMatrix_ValuesSetIntg1    
    MODULE PROCEDURE DistributedMatrix_ValuesSetIntg2
    MODULE PROCEDURE DistributedMatrix_ValuesSetSP
    MODULE PROCEDURE DistributedMatrix_ValuesSetSP1    
    MODULE PROCEDURE DistributedMatrix_ValuesSetSP2
    MODULE PROCEDURE DistributedMatrix_ValuesSetDP
    MODULE PROCEDURE DistributedMatrix_ValuesSetDP1    
    MODULE PROCEDURE DistributedMatrix_ValuesSetDP2
    MODULE PROCEDURE DistributedMatrix_ValuesSetL
    MODULE PROCEDURE DistributedMatrix_ValuesSetL1    
    MODULE PROCEDURE DistributedMatrix_ValuesSetL2
  END INTERFACE DistributedMatrix_ValuesSet

  INTERFACE DISTRIBUTED_MATRIX_BY_VECTOR_ADD
    MODULE PROCEDURE DistributedMatrix_MatrixByVectorAdd
  END INTERFACE DISTRIBUTED_MATRIX_BY_VECTOR_ADD
  
  INTERFACE DISTRIBUTED_VECTOR_ALL_VALUES_SET
    MODULE PROCEDURE DistributedVector_AllValuesSetIntg
    MODULE PROCEDURE DistributedVector_AllValuesSetSP
    MODULE PROCEDURE DistributedVector_AllValuesSetDP
    MODULE PROCEDURE DistributedVector_AllValuesSetL
  END INTERFACE DISTRIBUTED_VECTOR_ALL_VALUES_SET

  INTERFACE DistributedVector_AllValuesSet
    MODULE PROCEDURE DistributedVector_AllValuesSetIntg
    MODULE PROCEDURE DistributedVector_AllValuesSetSP
    MODULE PROCEDURE DistributedVector_AllValuesSetDP
    MODULE PROCEDURE DistributedVector_AllValuesSetL
  END INTERFACE DistributedVector_AllValuesSet

  INTERFACE DISTRIBUTED_VECTOR_COPY
    MODULE PROCEDURE DistributedVector_CopyIntg
    MODULE PROCEDURE DistributedVector_CopySP
    MODULE PROCEDURE DistributedVector_CopyDP
    MODULE PROCEDURE DistributedVector_CopyL
  END INTERFACE DISTRIBUTED_VECTOR_COPY
  
  INTERFACE DistributedVector_Copy
    MODULE PROCEDURE DistributedVector_CopyIntg
    MODULE PROCEDURE DistributedVector_CopySP
    MODULE PROCEDURE DistributedVector_CopyDP
    MODULE PROCEDURE DistributedVector_CopyL
 END INTERFACE DistributedVector_Copy

  INTERFACE DISTRIBUTED_VECTOR_CREATE_FINISH
    MODULE PROCEDURE DistributedVector_CreateFinish
  END INTERFACE DISTRIBUTED_VECTOR_CREATE_FINISH

  INTERFACE DISTRIBUTED_VECTOR_CREATE_START
    MODULE PROCEDURE DistributedVector_CreateStart
  END INTERFACE DISTRIBUTED_VECTOR_CREATE_START
  
  INTERFACE DISTRIBUTED_VECTOR_DATA_GET
    MODULE PROCEDURE DistributedVector_DataGetIntg
    MODULE PROCEDURE DistributedVector_DataGetSP
    MODULE PROCEDURE DistributedVector_DataGetDP
    MODULE PROCEDURE DistributedVector_DataGetL
  END INTERFACE DISTRIBUTED_VECTOR_DATA_GET

  INTERFACE DistributedVector_DataGet
    MODULE PROCEDURE DistributedVector_DataGetIntg
    MODULE PROCEDURE DistributedVector_DataGetSP
    MODULE PROCEDURE DistributedVector_DataGetDP
    MODULE PROCEDURE DistributedVector_DataGetL
  END INTERFACE DistributedVector_DataGet

  INTERFACE DISTRIBUTED_VECTOR_DATA_RESTORE
    MODULE PROCEDURE DistributedVector_DataRestoreIntg
    MODULE PROCEDURE DistributedVector_DataRestoreSP
    MODULE PROCEDURE DistributedVector_DataRestoreDP
    MODULE PROCEDURE DistributedVector_DataRestoreL
  END INTERFACE DISTRIBUTED_VECTOR_DATA_RESTORE

  INTERFACE DistributedVector_DataRestore
    MODULE PROCEDURE DistributedVector_DataRestoreIntg
    MODULE PROCEDURE DistributedVector_DataRestoreSP
    MODULE PROCEDURE DistributedVector_DataRestoreDP
    MODULE PROCEDURE DistributedVector_DataRestoreL
  END INTERFACE DistributedVector_DataRestore

  INTERFACE DISTRIBUTED_VECTOR_DATA_TYPE_SET
    MODULE PROCEDURE DistributedVector_DataTypeSet
  END INTERFACE DISTRIBUTED_VECTOR_DATA_TYPE_SET

  INTERFACE DISTRIBUTED_VECTOR_DESTROY
    MODULE PROCEDURE DistributedVector_Destroy
  END INTERFACE DISTRIBUTED_VECTOR_DESTROY

  INTERFACE DISTRIBUTED_VECTOR_DUPLICATE
    MODULE PROCEDURE DistributedVector_Duplicate
  END INTERFACE DISTRIBUTED_VECTOR_DUPLICATE

  INTERFACE DISTRIBUTED_VECTOR_GHOSTING_TYPE_SET
    MODULE PROCEDURE DistributedVector_GhostingTypeSet
  END INTERFACE DISTRIBUTED_VECTOR_GHOSTING_TYPE_SET

  INTERFACE DISTRIBUTED_VECTOR_LIBRARY_TYPE_SET
    MODULE PROCEDURE DistributedVector_LibraryTypeSet
  END INTERFACE DISTRIBUTED_VECTOR_LIBRARY_TYPE_SET

  INTERFACE DISTRIBUTED_VECTOR_OUTPUT
    MODULE PROCEDURE DistributedVector_Output
  END INTERFACE DISTRIBUTED_VECTOR_OUTPUT

  INTERFACE DISTRIBUTED_VECTOR_OVERRIDE_SET_ON
    MODULE PROCEDURE DistributedVector_OverrideSetOn
  END INTERFACE DISTRIBUTED_VECTOR_OVERRIDE_SET_ON

  INTERFACE DISTRIBUTED_VECTOR_OVERRIDE_SET_OFF
    MODULE PROCEDURE DistributedVector_OverrideSetOff
  END INTERFACE DISTRIBUTED_VECTOR_OVERRIDE_SET_OFF

  INTERFACE DISTRIBUTED_VECTOR_UPDATE_FINISH
    MODULE PROCEDURE DistributedVector_UpdateFinish
  END INTERFACE DISTRIBUTED_VECTOR_UPDATE_FINISH

  INTERFACE DISTRIBUTED_VECTOR_UPDATE_START
    MODULE PROCEDURE DistributedVector_UpdateStart
  END INTERFACE DISTRIBUTED_VECTOR_UPDATE_START

  INTERFACE DISTRIBUTED_VECTOR_UPDATE_ISFINISHED
    MODULE PROCEDURE DistributedVector_UpdateIsFinished
  END INTERFACE DISTRIBUTED_VECTOR_UPDATE_ISFINISHED

  INTERFACE DISTRIBUTED_VECTOR_UPDATE_WAITFINISHED
    MODULE PROCEDURE DistributedVector_UpdateWaitFinished
  END INTERFACE DISTRIBUTED_VECTOR_UPDATE_WAITFINISHED

  INTERFACE DISTRIBUTED_VECTOR_VALUES_ADD
    MODULE PROCEDURE DistributedVector_ValuesAddIntg
    MODULE PROCEDURE DistributedVector_ValuesAddIntg1
    MODULE PROCEDURE DistributedVector_ValuesAddSP
    MODULE PROCEDURE DistributedVector_ValuesAddSP1
    MODULE PROCEDURE DistributedVector_ValuesAddDP
    MODULE PROCEDURE DistributedVector_ValuesAddDP1
    MODULE PROCEDURE DistributedVector_ValuesAddL
    MODULE PROCEDURE DistributedVector_ValuesAddL1
  END INTERFACE DISTRIBUTED_VECTOR_VALUES_ADD

  INTERFACE DistributedVector_ValuesAdd
    MODULE PROCEDURE DistributedVector_ValuesAddIntg
    MODULE PROCEDURE DistributedVector_ValuesAddIntg1
    MODULE PROCEDURE DistributedVector_ValuesAddSP
    MODULE PROCEDURE DistributedVector_ValuesAddSP1
    MODULE PROCEDURE DistributedVector_ValuesAddDP
    MODULE PROCEDURE DistributedVector_ValuesAddDP1
    MODULE PROCEDURE DistributedVector_ValuesAddL
    MODULE PROCEDURE DistributedVector_ValuesAddL1
  END INTERFACE DistributedVector_ValuesAdd

  INTERFACE DISTRIBUTED_VECTOR_VALUES_GET
    MODULE PROCEDURE DistributedVector_ValuesGetIntg
    MODULE PROCEDURE DistributedVector_ValuesGetIntg1
    MODULE PROCEDURE DistributedVector_ValuesGetSP
    MODULE PROCEDURE DistributedVector_ValuesGetSP1
    MODULE PROCEDURE DistributedVector_ValuesGetDP
    MODULE PROCEDURE DistributedVector_ValuesGetDP1
    MODULE PROCEDURE DistributedVector_ValuesGetL
    MODULE PROCEDURE DistributedVector_ValuesGetL1
  END INTERFACE DISTRIBUTED_VECTOR_VALUES_GET

  INTERFACE DistributedVector_ValuesGet
    MODULE PROCEDURE DistributedVector_ValuesGetIntg
    MODULE PROCEDURE DistributedVector_ValuesGetIntg1
    MODULE PROCEDURE DistributedVector_ValuesGetSP
    MODULE PROCEDURE DistributedVector_ValuesGetSP1
    MODULE PROCEDURE DistributedVector_ValuesGetDP
    MODULE PROCEDURE DistributedVector_ValuesGetDP1
    MODULE PROCEDURE DistributedVector_ValuesGetL
    MODULE PROCEDURE DistributedVector_ValuesGetL1
  END INTERFACE DistributedVector_ValuesGet

  INTERFACE DISTRIBUTED_VECTOR_VALUES_SET
    MODULE PROCEDURE DistributedVector_ValuesSetIntg
    MODULE PROCEDURE DistributedVector_ValuesSetIntg1
    MODULE PROCEDURE DistributedVector_ValuesSetSP
    MODULE PROCEDURE DistributedVector_ValuesSetSP1
    MODULE PROCEDURE DistributedVector_ValuesSetDP
    MODULE PROCEDURE DistributedVector_ValuesSetDP1
    MODULE PROCEDURE DistributedVector_ValuesSetL
    MODULE PROCEDURE DistributedVector_ValuesSetL1
  END INTERFACE DISTRIBUTED_VECTOR_VALUES_SET

  INTERFACE DistributedVector_ValuesSet
    MODULE PROCEDURE DistributedVector_ValuesSetIntg
    MODULE PROCEDURE DistributedVector_ValuesSetIntg1
    MODULE PROCEDURE DistributedVector_ValuesSetSP
    MODULE PROCEDURE DistributedVector_ValuesSetSP1
    MODULE PROCEDURE DistributedVector_ValuesSetDP
    MODULE PROCEDURE DistributedVector_ValuesSetDP1
    MODULE PROCEDURE DistributedVector_ValuesSetL
    MODULE PROCEDURE DistributedVector_ValuesSetL1
  END INTERFACE DistributedVector_ValuesSet

  INTERFACE DistributedVector_DotProduct
    MODULE PROCEDURE DistributedVector_DotProductIntg
    MODULE PROCEDURE DistributedVector_DotProductSp
    MODULE PROCEDURE DistributedVector_DotProductDp
  END INTERFACE DistributedVector_DotProduct

  INTERFACE DISTRIBUTED_MATRIX_LINKLIST_SET
    MODULE PROCEDURE DistributedMatrix_LinkListSet
  END INTERFACE DISTRIBUTED_MATRIX_LINKLIST_SET

  INTERFACE DISTRIBUTED_MATRIX_LINKLIST_GET
    MODULE PROCEDURE DistributedMatrix_LinkListGet
  END INTERFACE DISTRIBUTED_MATRIX_LINKLIST_GET

  PUBLIC DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE,DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE

  PUBLIC DISTRIBUTED_MATRIX_VECTOR_INTG_TYPE,DISTRIBUTED_MATRIX_VECTOR_SP_TYPE,DISTRIBUTED_MATRIX_VECTOR_DP_TYPE, &
    & DISTRIBUTED_MATRIX_VECTOR_L_TYPE,DISTRIBUTED_MATRIX_VECTOR_SPC_TYPE,DISTRIBUTED_MATRIX_VECTOR_DPC_TYPE

  PUBLIC DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE,DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE, &
    & DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE,DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE, &
    & DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE, &
    & DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE

  PUBLIC DISTRIBUTED_MATRIX_SYMMETRIC_TYPE,DISTRIBUTED_MATRIX_HERMITIAN_TYPE,DISTRIBUTED_MATRIX_SKEW_SYMMETRIC_TYPE, &
    & DISTRIBUTED_MATRIX_UNSYMMETRIC_TYPE

  PUBLIC DISTRIBUTED_MATRIX_VECTOR_INCLUDE_GHOSTS_TYPE,DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE
  
  PUBLIC DISTRIBUTED_MATRIX_ALL_VALUES_SET

  PUBLIC DistributedMatrix_AllValuesSet

  PUBLIC DISTRIBUTED_MATRIX_CREATE_FINISH,DISTRIBUTED_MATRIX_CREATE_START

  PUBLIC DistributedMatrix_CreateFinish,DistributedMatrix_CreateStart

  PUBLIC DISTRIBUTED_MATRIX_DATA_GET,DISTRIBUTED_MATRIX_DATA_RESTORE

  PUBLIC DistributedMatrix_DataGet,DistributedMatrix_DataRestore

  PUBLIC DISTRIBUTED_MATRIX_DATA_TYPE_SET

  PUBLIC DistributedMatrix_DataTypeGet,DistributedMatrix_DataTypeSet

  PUBLIC DistributedMatrix_DimensionsGet

  PUBLIC DISTRIBUTED_MATRIX_DESTROY

  PUBLIC DistributedMatrix_Destroy

  PUBLIC DISTRIBUTED_MATRIX_DUPLICATE

  PUBLIC DistributedMatrix_Duplicate

  PUBLIC DISTRIBUTED_MATRIX_FORM

  PUBLIC DistributedMatrix_Form

  PUBLIC DISTRIBUTED_MATRIX_GHOSTING_TYPE_SET

  PUBLIC DistributedMatrix_GhostingTypeSet

  PUBLIC DISTRIBUTED_MATRIX_LIBRARY_TYPE_SET

  PUBLIC DistributedMatrix_LibraryTypeSet

  PUBLIC DISTRIBUTED_MATRIX_LINKLIST_SET,DISTRIBUTED_MATRIX_LINKLIST_GET

  PUBLIC DistributedMatrix_LinkListSet,DistributedMatrix_LinkListGet

  PUBLIC DISTRIBUTED_MATRIX_MAX_COLUMNS_PER_ROW_GET

  PUBLIC DistributedMatrix_MaxColumnsPerRowGet

  PUBLIC DISTRIBUTED_MATRIX_NUMBER_NON_ZEROS_SET,DISTRIBUTED_MATRIX_NUMBER_NON_ZEROS_GET

  PUBLIC DistributedMatrix_NumberOfNonZerosGet,DistributedMatrix_NumberOfNonZerosSet

  PUBLIC DISTRIBUTED_MATRIX_OUTPUT

  PUBLIC DistributedMatrix_Output

  PUBLIC DISTRIBUTED_MATRIX_OVERRIDE_SET_ON,DISTRIBUTED_MATRIX_OVERRIDE_SET_OFF

  PUBLIC DistributedMatrix_OverrideSetOn,DistributedMatrix_OverrideSetOff

  PUBLIC DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET,DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_SET

  PUBLIC DistributedMatrix_StorageLocationsGet,DistributedMatrix_StorageLocationsSet

  PUBLIC DISTRIBUTED_MATRIX_STORAGE_TYPE_GET,DISTRIBUTED_MATRIX_STORAGE_TYPE_SET

  PUBLIC DistributedMatrix_StorageTypeGet,DistributedMatrix_StorageTypeSet

  PUBLIC DistributedMatrix_SymmetryTypeGet,DistributedMatrix_SymmetryTypeSet

  PUBLIC DISTRIBUTED_MATRIX_UPDATE_START,DISTRIBUTED_MATRIX_UPDATE_FINISH

  PUBLIC DistributedMatrix_UpdateFinish,DistributedMatrix_UpdateStart

  PUBLIC DISTRIBUTED_MATRIX_UPDATE_ISFINISHED,DISTRIBUTED_MATRIX_UPDATE_WAITFINISHED

  PUBLIC DistributedMatrix_UpdateIsFinished,DistributedMatrix_UpdateWaitFinished

  PUBLIC DISTRIBUTED_MATRIX_VALUES_ADD

  PUBLIC DistributedMatrix_ValuesAdd

  PUBLIC DISTRIBUTED_MATRIX_VALUES_GET,DISTRIBUTED_MATRIX_VALUES_SET

  PUBLIC DistributedMatrix_ValuesGet,DistributedMatrix_ValuesSet

  PUBLIC DISTRIBUTED_MATRIX_BY_VECTOR_ADD

  PUBLIC DistributedMatrix_MatrixByVectorAdd
  
  PUBLIC DISTRIBUTED_VECTOR_ALL_VALUES_SET

  PUBLIC DistributedVector_AllValuesSet

  PUBLIC DISTRIBUTED_VECTOR_COPY

  PUBLIC DistributedVector_Copy

  PUBLIC DISTRIBUTED_VECTOR_CREATE_FINISH,DISTRIBUTED_VECTOR_CREATE_START

  PUBLIC DistributedVector_CreateFinish,DistributedVector_CreateStart

  PUBLIC DISTRIBUTED_VECTOR_DATA_GET,DISTRIBUTED_VECTOR_DATA_RESTORE

  PUBLIC DistributedVector_DataGet,DistributedVector_DataRestore

  PUBLIC DISTRIBUTED_VECTOR_DATA_TYPE_SET

  PUBLIC DistributedVector_DataTypeGet,DistributedVector_DataTypeSet

  PUBLIC DISTRIBUTED_VECTOR_DESTROY

  PUBLIC DistributedVector_Destroy

  PUBLIC DistributedVector_DotProduct

  PUBLIC DISTRIBUTED_VECTOR_DUPLICATE

  PUBLIC DistributedVector_Duplicate

  PUBLIC DISTRIBUTED_VECTOR_GHOSTING_TYPE_SET

  PUBLIC DistributedVector_GhostingTypeSet

  PUBLIC DISTRIBUTED_VECTOR_LIBRARY_TYPE_SET

  PUBLIC DistributedVector_LibraryTypeSet

  PUBLIC DistributedVector_L2Norm

  PUBLIC DISTRIBUTED_VECTOR_OUTPUT

  PUBLIC DistributedVector_Output

  PUBLIC DISTRIBUTED_VECTOR_OVERRIDE_SET_ON,DISTRIBUTED_VECTOR_OVERRIDE_SET_OFF

  PUBLIC DistributedVector_OverrideSetOn,DistributedVector_OverrideSetOff

  PUBLIC DISTRIBUTED_VECTOR_UPDATE_START,DISTRIBUTED_VECTOR_UPDATE_FINISH

  PUBLIC DistributedVector_UpdateFinish,DistributedVector_UpdateStart

  PUBLIC DISTRIBUTED_VECTOR_UPDATE_ISFINISHED,DISTRIBUTED_VECTOR_UPDATE_WAITFINISHED

  PUBLIC DistributedVector_UpdateIsFinished,DistributedVector_UpdateWaitFinished
  
  PUBLIC DISTRIBUTED_VECTOR_VALUES_ADD

  PUBLIC DistributedVector_ValuesAdd

  PUBLIC DISTRIBUTED_VECTOR_VALUES_GET,DISTRIBUTED_VECTOR_VALUES_SET

  PUBLIC DistributedVector_ValuesGet,DistributedVector_ValuesSet

CONTAINS  
  
  !
  !================================================================================================================================
  !

  !>Sets all values in an integer distributed matrix to the specified value.
  SUBROUTINE DistributedMatrix_AllValuesSetIntg(distributedMatrix,value,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(IN) :: value !<The value to set 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedMatrix_AllValuesSetIntg",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_AllValuesSet(cmissMatrix%matrix,value,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot set all values for an integer PETSc distributed matrix.",err,error,*999)
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("DistributedMatrix_AllValuesSetIntg")
    RETURN
999 ERRORSEXITS("DistributedMatrix_AllValuesSetIntg",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_AllValuesSetIntg

  !
  !================================================================================================================================
  !

  !>Sets all values in a single precision distributed matrix to the specified value.
  SUBROUTINE DistributedMatrix_AllValuesSetSP(distributedMatrix,value,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    REAL(SP), INTENT(IN) :: value !<The value to set 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedMatrix_AllValuesSetSP",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_AllValuesSet(cmissMatrix%matrix,value,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot set all values for a single precision PETSc distributed matrix.",err,error,*999)
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("DistributedMatrix_AllValuesSetSP")
    RETURN
999 ERRORSEXITS("DistributedMatrix_AllValuesSetSP",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_AllValuesSetSP

  !
  !================================================================================================================================
  !

  !>Sets all values in a double precision distributed matrix to the specified value.
  SUBROUTINE DistributedMatrix_AllValuesSetDP(distributedMatrix,value,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    REAL(DP), INTENT(IN) :: value !<The value to set 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(DistributedMatrixPETScType), POINTER :: petscMatrix
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedMatrix_AllValuesSetDP",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_AllValuesSet(cmissMatrix%matrix,value,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      NULLIFY(petscMatrix)
      CALL DistributedMatrix_PETScMatrixGet(distributedMatrix,petscMatrix,err,error,*999)
      IF(ABS(VALUE)<=ZERO_TOLERANCE) THEN
        IF(petscMatrix%useOverrideMatrix) THEN
          CALL Petsc_MatZeroEntries(petscMatrix%overrideMatrix,err,error,*999)
        ELSE
          CALL Petsc_MatZeroEntries(petscMatrix%matrix,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Not implemented.",err,error,*999)
      ENDIF
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("DistributedMatrix_AllValuesSetDP")
    RETURN
999 ERRORSEXITS("DistributedMatrix_AllValuesSetDP",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_AllValuesSetDP

  !
  !================================================================================================================================
  !

  !>Sets all values in a logical distributed matrix to the specified value.
  SUBROUTINE DistributedMatrix_AllValuesSetL(distributedMatrix,value,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    LOGICAL, INTENT(IN) :: VALUE !<The value to set 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedMatrix_AllValuesSetL",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_AllValuesSet(cmissMatrix%matrix,value,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot set all values for a logical PETSc distributed matrix.",err,error,*999)
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("DistributedMatrix_AllValuesSetL")
    RETURN
999 ERRORSEXITS("DistributedMatrix_AllValuesSetL",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_AllValuesSetL

  !
  !================================================================================================================================
  !

  !>Finishes the creation of a CMISS distributed matrix.
  SUBROUTINE DistributedMatrix_CMISSCreateFinish(cmissMatrix,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix !<A pointer to the distributed CMISS matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: domainMapping
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("DistributedMatrix_CMISSCreateFinish",err,error,*998)

    IF(.NOT.ASSOCIATED(cmissMatrix)) CALL FlagError("Distributed matrix CMISS is not associated.",err,error,*998)
    
    cmissMatrix%baseTagNumber=distributedDataId
    NULLIFY(domainMapping)
    domainMapping=>cmissMatrix%distributedMatrix%rowDomainMapping
    IF(.NOT.ASSOCIATED(domainMapping)) CALL FlagError("Distributed matrix row domain mapping is not associated.",err,error,*998)
    IF(domainMapping%NUMBER_OF_DOMAINS==1) THEN
      distributedDataId=distributedDataId+1
    ELSE
      distributedDataId=distributedDataId+ &
        & domainMapping%ADJACENT_DOMAINS_PTR(domainMapping%NUMBER_OF_DOMAINS)
    END IF
    CALL Matrix_CreateFinish(cmissMatrix%matrix,err,error,*999)
   
    EXITS("DistributedMatrix_CMISSCreateFinish")
    RETURN
999 CALL DistributedMatrix_CMISSFinalise(cmissMatrix,dummyErr,dummyError,*998)
998 ERRORSEXITS("DistributedMatrix_CMISSCreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_CMISSCreateFinish

  !
  !================================================================================================================================
  !

  !>Finalise a CMISS distributed matrix.
  SUBROUTINE DistributedMatrix_CMISSFinalise(cmissMatrix,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix !<A pointer to the CMISS distributed matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("DistributedMatrix_CMISSFinalise",err,error,*999)

    IF(ASSOCIATED(cmissMatrix)) THEN
      CALL Matrix_Destroy(cmissMatrix%matrix,err,error,*999)
      DEALLOCATE(cmissMatrix)
    ENDIF
    
    EXITS("DistributedMatrix_CMISSFinalise")
    RETURN
999 ERRORSEXITS("DistributedMatrix_CMISSFinalise",err,error)
    RETURN 1
  END SUBROUTINE DistributedMatrix_CMISSFinalise

  !
  !================================================================================================================================
  !

  !>Intialises a CMISS distributed matrix.
  SUBROUTINE DistributedMatrix_CMISSInitialise(distributedMatrix,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: rowDomainMapping,columnDomainMapping
    TYPE(VARYING_STRING) :: dummyError,localError
    
    ENTERS("DistributedMatrix_CMISSInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*998)
    IF(ASSOCIATED(distributedMatrix%cmiss)) &
      & CALL FlagError("CMISS is already associated for this distributed matrix.",err,error,*998)
    
    rowDomainMapping=>distributedMatrix%rowDomainMapping    
    columnDomainMapping=>distributedMatrix%columnDomainMapping
    IF(.NOT.ASSOCIATED(rowDomainMapping)) &
      & CALL FlagError("Distributed matrix row domain mapping is not associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(columnDomainMapping)) &
      & CALL FlagError("Distributed matrix column domain mapping is not associated.",err,error,*998)
    
    ALLOCATE(distributedMatrix%cmiss,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate CMISS distributed matrix.",err,error,*999)
    distributedMatrix%cmiss%distributedMatrix=>distributedMatrix
    distributedMatrix%libraryType=DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE
    NULLIFY(distributedMatrix%cmiss%matrix)
    !Set the defaults
    CALL Matrix_CreateStart(distributedMatrix%cmiss%matrix,err,error,*999)
    CALL Matrix_DataTypeSet(distributedMatrix%cmiss%matrix,MATRIX_VECTOR_DP_TYPE,err,error,*999)
    CALL Matrix_StorageTypeSet(distributedMatrix%cmiss%matrix,MATRIX_BLOCK_STORAGE_TYPE,err,error,*999)
    SELECT CASE(distributedMatrix%ghostingType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_INCLUDE_GHOSTS_TYPE)
      CALL Matrix_SizeSet(distributedMatrix%cmiss%matrix,rowDomainMapping%TOTAL_NUMBER_OF_LOCAL, &
        & columnDomainMapping%NUMBER_OF_GLOBAL,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE)
      CALL Matrix_SizeSet(distributedMatrix%cmiss%matrix,rowDomainMapping%NUMBER_OF_LOCAL, &
        & columnDomainMapping%NUMBER_OF_GLOBAL,err,error,*999)
    CASE DEFAULT
      localError="The distributed matrix ghosting type of "// &
        & TRIM(NumberToVString(distributedMatrix%ghostingType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_CMISSInitialise")
    RETURN
999 IF(ASSOCIATED(distributedMatrix%cmiss)) &
      & CALL DistributedMatrix_CMISSFinalise(distributedMatrix%cmiss,dummyErr,dummyError,*999)
998 ERRORSEXITS("DistributedMatrix_CMISSInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_CMISSInitialise

  !
  !================================================================================================================================
  !

  !>Finishes the creation of a distributed matrix.
  SUBROUTINE DistributedMatrix_CreateFinish(distributedMatrix,err,error,*)

    !Argument variables 
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError,localError
    
    ENTERS("DistributedMatrix_CreateFinish",err,error,*998)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*998)
    IF(distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has been finished.",err,error,*998)
    
    SELECT CASE(distributedMatrix%libraryType)      
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      CALL DistributedMatrix_CMISSCreateFinish(distributedMatrix%cmiss,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL DistributedMatrix_PETScCreateFinish(distributedMatrix%petsc,err,error,*999)
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    distributedMatrix%matrixFinished=.TRUE.
    
    EXITS("DistributedMatrix_CreateFinish")
    RETURN
999 CALL DistributedMatrix_Finalise(distributedMatrix,dummyErr,dummyError,*998)
998 ERRORSEXITS("DistributedMatrix_CreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_CreateFinish

  !
  !================================================================================================================================
  !

  !>Starts the creation of a distributed matrix.
  SUBROUTINE DistributedMatrix_CreateStart(rowDomainMapping,columnDomainMapping,distributedMatrix,err,error,*)

    !Argument variables
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: rowDomainMapping !<A pointer to the row domain mapping to be used for the distribution
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: columnDomainMapping !<A pointer to the column domain mapping to be used for the distribution
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<On return a pointer to the distributed matrix being created
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError,localError
    
    ENTERS("DistributedMatrix_CreateStart",err,error,*998)

    IF(.NOT.ASSOCIATED(rowDomainMapping)) CALL FlagError("Row domain mapping is not associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(columnDomainMapping)) CALL FlagError("Column domain mapping is not associated.",err,error,*998)
    IF(ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is already associated.",err,error,*998)
    
    IF(rowDomainMapping%NUMBER_OF_DOMAINS==columnDomainMapping%NUMBER_OF_DOMAINS) THEN
      CALL DistributedMatrix_Initialise(rowDomainMapping,columnDomainMapping,distributedMatrix,err,error,*999)
      !Set the defaults
    ELSE
      localError="The number of domains in the row domain mapping ("// &
        & TRIM(NumberToVString(rowDomainMapping%NUMBER_OF_DOMAINS,"*",err,error))// &
        & ") does not match the number of domains in the column domain mapping ("// &
        & TRIM(NumberToVString(columnDomainMapping%NUMBER_OF_DOMAINS,"*",err,error))//")."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("DistributedMatrix_CreateStart")
    RETURN
999 CALL DistributedMatrix_Finalise(distributedMatrix,dummyErr,dummyError,*998)
998 ERRORSEXITS("DistributedMatrix_CreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_CreateStart

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the data of an integer distributed matrix. Note: the values can be used for read operations but a DistributedMatrix_ValuesSet call must be used to change any values. The pointer should not be deallocated and a DistributedMatrix_DataRestore call must be used after the data has finished being used.
  SUBROUTINE DistributedMatrix_DataGetIntg(distributedMatrix,data,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), POINTER :: data(:) !<On return a pointer to the distributed matrix data for this computational node. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedMatrix_DataGetIntg",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated",err,error,*999)
    IF(ASSOCIATED(data)) CALL FlagError("Data is already associated",err,error,*999)    
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_DataGet(distributedMatrix%cmiss%matrix,data,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot get data for an integer PETSc distributed matrix.",err,error,*999)
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_DataGetIntg")
    RETURN
999 ERRORSEXITS("DistributedMatrix_DataGetIntg",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_DataGetIntg

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the data of a single precision real distributed matrix. Note: the values can be used for read operations but a DistributedMatrix_ValuesSet call must be used to change any values. The pointer should not be deallocated and a DistributedMatrix_DataRestore call must be used after the data has finished being used.
  SUBROUTINE DistributedMatrix_DataGetSP(distributedMatrix,DATA,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    REAL(SP), POINTER :: data(:) !<On return a pointer to the distributed matrix data for this computational node. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedMatrix_DataGetSP",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated",err,error,*999)
    IF(ASSOCIATED(data)) CALL FlagError("Data is already associated",err,error,*999)    
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_DataGet(distributedMatrix%cmiss%matrix,data,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot get data for a single precision real PETSc distributed matrix.",err,error,*999)
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_DataGetSP")
    RETURN
999 ERRORSEXITS("DistributedMatrix_DataGetSP",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_DataGetSP

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the data of a double precision real distributed matrix. Note: the values can be used for read operations but a DistributedMatrix_ValuesSet call must be used to change any values. The pointer should not be deallocated and a DistributedMatrix_DataRestore call must be used after the data has finished being used.
  SUBROUTINE DistributedMatrix_DataGetDP(distributedMatrix,data,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    REAL(DP), POINTER :: data(:) !<On return a pointer to the distributed matrix data for this computational node Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    REAL(DP), POINTER :: petscData(:,:)
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(DistributedMatrixPETScType), POINTER :: petscMatrix
    TYPE(VARYING_STRING) :: localError
    TYPE(C_PTR) :: tempMe

    ENTERS("DistributedMatrix_DataGetDP",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated",err,error,*999)
    IF(ASSOCIATED(data)) CALL FlagError("Data is already associated",err,error,*999)    
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_DataGet(distributedMatrix%cmiss%matrix,DATA,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      NULLIFY(petscMatrix)
      CALL DistributedMatrix_PETScMatrixGet(distributedMatrix,petscMatrix,err,error,*999)
      IF(petscMatrix%useOverrideMatrix) THEN
        SELECT CASE(petscMatrix%storageType)
        CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
          CALL Petsc_MatDenseGetArrayF90(distributedMatrix%petsc%overrideMatrix,petscData,err,error,*999)
        CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
          CALL FlagError("Diagonal storage is not implemented for PETSc matrices.",err,error,*999)
        CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
          CALL FlagError("Column major storage is not implemented for PETSc matrices.",err,error,*999)
        CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
          CALL FlagError("Row major storage is not implemented for PETSc matrices.",err,error,*999)
        CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
          CALL Petsc_MatSeqAIJGetArrayF90(distributedMatrix%petsc%overrideMatrix,petscData,err,error,*999)
        CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
          CALL FlagError("Compressed column storage is not implemented for PETSc matrices.",err,error,*999)
        CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
          CALL FlagError("Row column storage is not implemented for PETSc matrices.",err,error,*999)
        CASE DEFAULT
          localError="The PETSc matrix storage type of "//TRIM(NumberToVString( &
            & petscMatrix%storageType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        SELECT CASE(petscMatrix%storageType)
        CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
          CALL Petsc_MatDenseGetArrayF90(petscMatrix%matrix,petscData,err,error,*999)
        CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
          CALL FlagError("Diagonal storage is not implemented for PETSc matrices.",err,error,*999)
        CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
          CALL FlagError("Column major storage is not implemented for PETSc matrices.",err,error,*999)
        CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
          CALL FlagError("Row major storage is not implemented for PETSc matrices.",err,error,*999)
        CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
          CALL Petsc_MatSeqAIJGetArrayF90(petscMatrix%matrix,petscData,err,error,*999)
        CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
          CALL FlagError("Compressed column storage is not implemented for PETSc matrices.",err,error,*999)
        CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
          CALL FlagError("Row column storage is not implemented for PETSc matrices.",err,error,*999)
        CASE DEFAULT
          localError="The PETSc matrix storage type of "//TRIM(NumberToVString( &
            & petscMatrix%storageType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDIF
      ! Convert 2D array from PETSc to 1D array
      ! Using C_F_POINTER(C_LOC(... is a bit ugly but transfer doesn't work with pointers
      SELECT CASE(petscMatrix%storageType)
      CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
        tempMe = C_LOC(petscData(1,1))
        CALL C_F_POINTER(tempMe,DATA,[petscMatrix%m*petscMatrix%n])
      CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
        CALL FlagError("Diagonal storage is not implemented for PETSc matrices.",err,error,*999)
      CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
        CALL FlagError("Column major storage is not implemented for PETSc matrices.",err,error,*999)
      CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
        CALL FlagError("Row major storage is not implemented for PETSc matrices.",err,error,*999)
      CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
        !PETSc returns an m * n matrix rather than number non-zeros by 1, so the returned
        !2D array actually contains junk data outside of the actual matrix.
        !This is a bug in PETSc but we can get the correct 1D data here
        tempMe = C_LOC(petscData(1,1))
        CALL C_F_POINTER(tempMe,DATA,[petscMatrix%numberOfNonZeros])
      CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
        CALL FlagError("Compressed column storage is not implemented for PETSc matrices.",err,error,*999)
      CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
        CALL FlagError("Row column storage is not implemented for PETSc matrices.",err,error,*999)
      CASE DEFAULT
        localError="The PETSc matrix storage type of "//TRIM(NumberToVString( &
          & petscMatrix%storageType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_DataGetDP")
    RETURN
999 ERRORSEXITS("DistributedMatrix_DataGetDP",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_DataGetDP

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the data of a logical distributed matrix. Note: the values can be used for read operations but a DistributedMatrix_ValuesSet call must be used to change any values. The pointer should not be deallocated and a DistributedMatrix_DataRestore call must be used after the data has finished being used.
  SUBROUTINE DistributedMatrix_DataGetL(distributedMatrix,DATA,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    LOGICAL, POINTER :: data(:) !<On return a pointer to the distributed matrix data for this computational node. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedMatrix_DataGetL",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated",err,error,*999)
    IF(ASSOCIATED(data)) CALL FlagError("Data is already associated",err,error,*999)    
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_DataGet(distributedMatrix%cmiss%matrix,data,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot get data for a logical PETSc distributed matrix.",err,error,*999)
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_DataGetL")
    RETURN
999 ERRORSEXITS("DistributedMatrix_DataGetL",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_DataGetL

  !
  !================================================================================================================================
  !

  !>Restores the integer data pointer returned from DistributedMatrix_DataGet once the data has finished being used.
  SUBROUTINE DistributedMatrix_DataRestoreIntg(distributedMatrix,data,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), POINTER :: data(:) !<The a pointer to the distributed matrix data for this computational node
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedMatrix_DataRestoreIntg",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(data)) CALL FlagError("Data is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(data)              
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot restore data for an integer PETSc distributed matrix.",err,error,*999)
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_DataRestoreIntg")
    RETURN
999 ERRORSEXITS("DistributedMatrix_DataRestoreIntg",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_DataRestoreIntg

  !
  !================================================================================================================================
  !

  !>Restores the single precision real data pointer returned from DistributedMatrix_DataGet once the data has finished being used.
  SUBROUTINE DistributedMatrix_DataRestoreSP(distributedMatrix,data,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    REAL(SP), POINTER :: data(:) !<The a pointer to the distributed matrix data for this computational node
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedMatrix_DataRestoreSP",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(data)) CALL FlagError("Data is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(data)              
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot restore data for a single precision real PETSc distributed matrix.",err,error,*999)
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_DataRestoreSP")
    RETURN
999 ERRORSEXITS("DistributedMatrix_DataRestoreSP",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_DataRestoreSP

  !
  !================================================================================================================================
  !

  !>Restores the double precision data pointer returned from DistributedMatrix_DataGet once the data has finished being used.
  SUBROUTINE DistributedMatrix_DataRestoreDP(distributedMatrix,DATA,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    REAL(DP), POINTER :: data(:) !<A pointer to the distributed matrix data for this computational node
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    REAL(DP), POINTER :: petscData(:,:)
    TYPE(DistributedMatrixPETScType), POINTER :: petscMatrix
    TYPE(VARYING_STRING) :: localError
    TYPE(C_PTR) :: tempMe
    
    ENTERS("DistributedMatrix_DataRestoreDP",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(data)) CALL FlagError("Data is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)      
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(data)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      NULLIFY(petscMatrix)
      CALL DistributedMatrix_PETScMatrixGet(distributedMatrix,petscMatrix,err,error,*999)
      SELECT CASE(petscMatrix%storageType)
      CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
        !Convert 1D array to 2D
        tempMe = C_LOC(data(1))
        CALL C_F_POINTER(tempMe,petscData,[petscMatrix%m,petscMatrix%n])
      CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
        CALL FlagError("Diagonal storage is not implemented for PETSc matrices.",err,error,*999)
      CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
        CALL FlagError("Column major storage is not implemented for PETSc matrices.",err,error,*999)
      CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
        CALL FlagError("Row major storage is not implemented for PETSc matrices.",err,error,*999)
      CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
        !PETSc expects an m * n 2D matrix rather than a 1D array with length equal to number of non-zeros
        !This is a bug in PETSc so we have to give it a 2D matrix with junk at the end
        tempMe = C_LOC(data(1))
        CALL C_F_POINTER(tempMe,petscData,[petscMatrix%m,petscMatrix%n])
      CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
        CALL FlagError("Compressed column storage is not implemented for PETSc matrices.",err,error,*999)
      CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
        CALL FlagError("Row column storage is not implemented for PETSc matrices.",err,error,*999)
      CASE DEFAULT
        localError="The PETSc matrix storage type of "//TRIM(NumberToVString(petscMatrix%storageType,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      IF(petscMatrix%useOverrideMatrix) THEN
        SELECT CASE(petscMatrix%storageType)
        CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
          CALL Petsc_MatDenseRestoreArrayF90(petscMatrix%overrideMatrix,petscData,err,error,*999)
        CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
          CALL FlagError("Diagonal storage is not implemented for PETSc matrices.",err,error,*999)
        CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
          CALL FlagError("Column major storage is not implemented for PETSc matrices.",err,error,*999)
        CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
          CALL FlagError("Row major storage is not implemented for PETSc matrices.",err,error,*999)
        CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
          CALL Petsc_MatSeqAIJRestoreArrayF90(petscMatrix%overrideMatrix,petscData,err,error,*999)
        CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
          CALL FlagError("Compressed column storage is not implemented for PETSc matrices.",err,error,*999)
        CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
          CALL FlagError("Row column storage is not implemented for PETSc matrices.",err,error,*999)
        CASE DEFAULT
          localError="The PETSc matrix storage type of "//TRIM(NumberToVString(petscMatrix%storageType,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        SELECT CASE(petscMatrix%storageType)
        CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
          CALL Petsc_MatDenseRestoreArrayF90(petscMatrix%matrix,petscData,err,error,*999)
        CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
          CALL FlagError("Diagonal storage is not implemented for PETSc matrices.",err,error,*999)
        CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
          CALL FlagError("Column major storage is not implemented for PETSc matrices.",err,error,*999)
        CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
          CALL FlagError("Row major storage is not implemented for PETSc matrices.",err,error,*999)
        CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
          CALL Petsc_MatSeqAIJRestoreArrayF90(petscMatrix%matrix,petscData,err,error,*999)
        CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
          CALL FlagError("Compressed column storage is not implemented for PETSc matrices.",err,error,*999)
        CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
          CALL FlagError("Row column storage is not implemented for PETSc matrices.",err,error,*999)
        CASE DEFAULT
          localError="The PETSc matrix storage type of "//TRIM(NumberToVString(petscMatrix%storageType,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDIF
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_DataRestoreDP")
    RETURN
999 ERRORSEXITS("DistributedMatrix_DataRestoreDP",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_DataRestoreDP

  !
  !================================================================================================================================
  !

  !>Restores the logical data pointer returned from DistributedMatrix_DataGet once the data has finished being used.
  SUBROUTINE DistributedMatrix_DataRestoreL(distributedMatrix,DATA,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    LOGICAL, POINTER :: data(:) !<The a pointer to the distributed matrix data for this computational node
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedMatrix_DataRestoreL",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(data)) CALL FlagError("Data is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(data)              
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot restore data for a logical PETSc distributed matrix.",err,error,*999)
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_DataRestoreL")
    RETURN
999 ERRORSEXITS("DistributedMatrix_DataRestoreL",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_DataRestoreL

  !
  !================================================================================================================================
  !

  !>Gets the data type of a distributed matrix.
  SUBROUTINE DistributedMatrix_DataTypeGet(distributedMatrix,dataType,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(OUT) :: dataType !<On return, the data type of the matrix. \see DISTRIBUTED_MATRIX_VECTOR_DataTypes,MATRIX_VECTOR
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("DistributedMatrix_DataTypeGet",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The matrix has not been finished.",err,error,*999)
    
    dataType=distributedMatrix%dataType
 
    EXITS("DistributedMatrix_DataTypeGet")
    RETURN
999 ERRORSEXITS("DistributedMatrix_DataTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_DataTypeGet

  !
  !================================================================================================================================
  !

  !>Sets/changes the data type of a distributed matrix.
  SUBROUTINE DistributedMatrix_DataTypeSet(distributedMatrix,dataType,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(IN) :: dataType !<The data type to set. \see MATRIX_VECTOR_DataTypes,MATRIX_VECTOR
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedMatrix_DataTypeSet",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has been finished.",err,error,*999)
     
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_DataTypeSet(cmissMatrix%matrix,dataType,err,error,*999)
      distributedMatrix%dataType=dataType
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      SELECT CASE(dataType)
      CASE(DISTRIBUTED_MATRIX_VECTOR_INTG_TYPE)
        CALL FlagError("An integer distributed PETSc matrix is not implemented.",err,error,*999)
      CASE(DISTRIBUTED_MATRIX_VECTOR_SP_TYPE)
        CALL FlagError("A single precision distributed PETSc matrix is not implemented.",err,error,*999)
      CASE(DISTRIBUTED_MATRIX_VECTOR_DP_TYPE)
        distributedMatrix%dataType=DISTRIBUTED_MATRIX_VECTOR_DP_TYPE
      CASE(DISTRIBUTED_MATRIX_VECTOR_L_TYPE)
        CALL FlagError("A logical distributed PETSc matrix is not implemented.",err,error,*999)
      CASE(DISTRIBUTED_MATRIX_VECTOR_SPC_TYPE)
        CALL FlagError("A single precision complex distributed PETSc matrix is not implemented.",err,error,*999)
      CASE(DISTRIBUTED_MATRIX_VECTOR_DPC_TYPE)
        CALL FlagError("A double precision complex distributed PETSc matrix is not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The specified data type of "//TRIM(NumberToVString(dataType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_DataTypeSet")
    RETURN
999 ERRORSEXITS("DistributedMatrix_DataTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_DataTypeSet

  !
  !================================================================================================================================
  !

  !>Gets the dimensions of a matrix on this computational node.
  SUBROUTINE DistributedMatrix_DimensionsGet(distributedMatrix,m,n,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix to get dimensions for
    INTEGER(INTG), INTENT(OUT) :: m !<On return, the number of rows in the matrix for this domain
    INTEGER(INTG), INTENT(OUT) :: n !<On return, the number of columns in the matrix for this domain
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(MATRIX_TYPE), POINTER :: matrix
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(DistributedMatrixPETScType), POINTER :: petscMatrix
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedMatrix_DimensionsGet",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
!!\TODO Move this to a matrix method
      matrix=>cmissMatrix%matrix
      IF(.NOT.ASSOCIATED(matrix)) CALL FlagError("Distributed matrix CMISS matrix is not associated.",err,error,*999)
      IF(.NOT.matrix%matrix_finished) CALL FlagError("The matrix has not been finished.",err,error,*999)
      m=matrix%m
      n=matrix%n
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

  !>Destroys a distributed matrix.
  SUBROUTINE DistributedMatrix_Destroy(distributedMatrix,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DistributedMatrix_Destroy",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    
    CALL DistributedMatrix_Finalise(distributedMatrix,err,error,*999)
   
    EXITS("DistributedMatrix_Destroy")
    RETURN
999 ERRORSEXITS("DistributedMatrix_Destroy",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_Destroy

  !
  !================================================================================================================================
  !

  !>Duplicates the structure of a distributed matrix and returns a pointer to the new matrix in newDistributedMatrix.
  SUBROUTINE DistributedMatrix_Duplicate(distributedMatrix,newDistributedMatrix,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix to duplicate
    TYPE(DistributedMatrixType), POINTER :: newDistributedMatrix !<On return, a pointer to the new duplicated distributed matrix. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix,newCmissMatrix
    TYPE(VARYING_STRING) :: dummyError,localError

    ENTERS("DistributedMatrix_Duplicate",err,error,*998)

    IF(ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*998)
    IF(ASSOCIATED(newDistributedMatrix)) CALL FlagError("New distributed matrix is already associated.",err,error,*998)
    
    SELECT CASE(distributedMatrix%libraryType)      
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*998)
      CALL DistributedMatrix_CreateStart(distributedMatrix%rowDomainMapping,distributedMatrix%columnDomainMapping, &
        & newDistributedMatrix,err,error,*999)
      NULLIFY(newCmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(newDistributedMatrix,newCmissMatrix,err,error,*999)
      CALL Matrix_Duplicate(cmissMatrix%matrix,newDistributedMatrix%cmiss%matrix,err,error,*999)
      CALL DistributedMatrix_CreateFinish(newDistributedMatrix,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL DistributedMatrix_CreateStart(distributedMatrix%rowDomainMapping,distributedMatrix%columnDomainMapping, &
        & newDistributedMatrix,err,error,*999)
      CALL DistributedMatrix_LibraryTypeSet(newDistributedMatrix,DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE,err,error,*999)
      CALL DistributedMatrix_CreateFinish(newDistributedMatrix,err,error,*999)
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_Duplicate")
    RETURN
999 CALL DistributedMatrix_Finalise(newDistributedMatrix,dummyErr,dummyError,*999)
998 ERRORSEXITS("DistributedMatrix_Duplicate",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_Duplicate

  !
  !================================================================================================================================
  !

  !>Finalises a distributed matrix and deallocates all memory.
  SUBROUTINE DistributedMatrix_Finalise(distributedMatrix,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DistributedMatrix_Finalise",err,error,*999)

    IF(ASSOCIATED(distributedMatrix)) THEN
      CALL DistributedMatrix_CMISSFinalise(distributedMatrix%cmiss,err,error,*999)
      CALL DistributedMatrix_PETScFinalise(distributedMatrix%petsc,err,error,*999)        
      DEALLOCATE(distributedMatrix)
    ENDIF
    
    EXITS("DistributedMatrix_Finalise")
    RETURN
999 ERRORSEXITS("DistributedMatrix_Finalise",err,error)
    RETURN 1
  END SUBROUTINE DistributedMatrix_Finalise

  !
  !================================================================================================================================
  !

  !>Forms a distributed matrix by initialising the structure of the matrix to zero
  SUBROUTINE DistributedMatrix_Form(distributedMatrix,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix to form.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: columnIdx,rowIdx
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(DistributedMatrixPETScType), POINTER :: petscMatrix
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedMatrix_Form",err,error,*999)

    IF(ASSOCIATED(distributedMatrix)) THEN
      IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
      SELECT CASE(distributedMatrix%libraryType)
      CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
        NULLIFY(cmissMatrix)
        CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
        NULLIFY(petscMatrix)
        CALL DistributedMatrix_PETScMatrixGet(distributedMatrix,petscMatrix,err,error,*999)
        SELECT CASE(petscMatrix%storageType)
        CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
          CALL Petsc_MatZeroEntries(petscMatrix%matrix,err,error,*999)
        CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
          CALL FlagError("Diagonal storage is not implemented for PETSc matrices.",err,error,*999)
        CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
          CALL FlagError("Column major storage is not implemented for PETSc matrices.",err,error,*999)
        CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
          CALL FlagError("Row major storage is not implemented for PETSc matrices.",err,error,*999)
        CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
          IF(petscMatrix%useOverrideMatrix) THEN
            DO rowIdx=1,petscMatrix%m
              DO columnIdx=petscMatrix%rowIndices(rowIdx),petscMatrix%rowIndices(rowIdx+1)-1
                CALL Petsc_MatSetValue(petscMatrix%overrideMatrix,petscMatrix%globalRowNumbers(rowIdx), &
                  & petscMatrix%columnIndices(columnIdx)-1,0.0_DP,PETSC_INSERT_VALUES, &
                  & err,error,*999) !PETSc uses 0 indicies
              ENDDO !columnIdx
            ENDDO !rowIdx
          ELSE
            DO rowIdx=1,petscMatrix%m
              DO columnIdx=petscMatrix%rowIndices(rowIdx),petscMatrix%rowIndices(rowIdx+1)-1
                CALL Petsc_MatSetValue(petscMatrix%matrix,petscMatrix%globalRowNumbers(rowIdx), &
                  & petscMatrix%columnIndices(columnIdx)-1,0.0_DP,PETSC_INSERT_VALUES, &
                  & err,error,*999) !PETSc uses 0 indicies
              ENDDO !columnIdx
            ENDDO !rowIdx
          ENDIF
        CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
          CALL FlagError("Compressed column storage is not implemented for PETSc matrices.",err,error,*999)
        CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
          CALL FlagError("Row column storage is not implemented for PETSc matrices.",err,error,*999)
        CASE DEFAULT
          localError="The PETSc matrix storage type of "//TRIM(NumberToVString(petscMatrix%storageType,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        IF(petscMatrix%useOverrideMatrix) THEN
          CALL Petsc_MatAssemblyBegin(petscMatrix%overrideMatrix,PETSC_MAT_FINAL_ASSEMBLY,err,error,*999)
          CALL Petsc_MatAssemblyEnd(petscMatrix%overrideMatrix,PETSC_MAT_FINAL_ASSEMBLY,err,error,*999)
        ELSE
          CALL Petsc_MatAssemblyBegin(petscMatrix%matrix,PETSC_MAT_FINAL_ASSEMBLY,err,error,*999)
          CALL Petsc_MatAssemblyEnd(petscMatrix%matrix,PETSC_MAT_FINAL_ASSEMBLY,err,error,*999)
        ENDIF
      CASE DEFAULT
        localError="The distributed matrix library type of "// &
          & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDIF
    
    EXITS("DistributedMatrix_Form")
    RETURN
999 ERRORSEXITS("DistributedMatrix_Form",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_Form

  !
  !================================================================================================================================
  !

  !>Sets/changes the ghosting type for a distributed matrix
  SUBROUTINE DistributedMatrix_GhostingTypeSet(distributedMatrix,ghostingType,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(IN) :: ghostingType !<The ghosting type \see DISTRIBUTED_MATRIX_VECTOR_GhostingTypes,DISTRIBUTED_MATRIX_VECTOR
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(DistributedMatrixPETScType), POINTER :: petscMatrix
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: rowDomainMapping,columnDomainMapping
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedMatrix_GhostingTypeSet",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has already been finished.",err,error,*999)

    NULLIFY(rowDomainMapping)
    CALL DistributedMatrix_RowMappingGet(distributedMatrix,rowDomainMapping,err,error,*999)
    NULLIFY(columnDomainMapping)
    CALL DistributedMatrix_ColumnMappingGet(distributedMatrix,columnDomainMapping,err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      SELECT CASE(ghostingType)
      CASE(DISTRIBUTED_MATRIX_VECTOR_INCLUDE_GHOSTS_TYPE)
        CALL Matrix_SizeSet(cmissMatrix%matrix,rowDomainMapping%TOTAL_NUMBER_OF_LOCAL, &
          & columnDomainMapping%NUMBER_OF_GLOBAL,err,error,*999)
      CASE(DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE)
        CALL Matrix_SizeSet(cmissMatrix%matrix,rowDomainMapping%NUMBER_OF_LOCAL, &
          & columnDomainMapping%NUMBER_OF_GLOBAL,err,error,*999)
      CASE DEFAULT
        localError="The given ghosting type of "//TRIM(NumberToVString(ghostingType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      NULLIFY(petscMatrix)
      CALL DistributedMatrix_PETScMatrixGet(distributedMatrix,petscMatrix,err,error,*999)
      SELECT CASE(ghostingType)
      CASE(DISTRIBUTED_MATRIX_VECTOR_INCLUDE_GHOSTS_TYPE)
        petscMatrix%n=rowDomainMapping%TOTAL_NUMBER_OF_LOCAL
      CASE(DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE)
        petscMatrix%n=rowDomainMapping%NUMBER_OF_LOCAL
      CASE DEFAULT
        localError="The given ghosting type of "//TRIM(NumberToVString(ghostingType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    distributedMatrix%ghostingType=ghostingType
    
    EXITS("DistributedMatrix_GhostingTypeSet")
    RETURN
999 ERRORSEXITS("DistributedMatrix_GhostingTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_GhostingTypeSet

  !
  !================================================================================================================================
  !

  !>Sets/changes the library type for a distributed matrix
  SUBROUTINE DistributedMatrix_LibraryTypeSet(distributedMatrix,libraryType,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix 
    INTEGER(INTG), INTENT(IN) :: libraryType !<The library type \see DISTRIBUTED_MATRIX_VECTOR_LibraryTypes,DISTRIBUTED_MATRIX_VECTOR
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,oldLibraryType
    TYPE(VARYING_STRING) :: dummyError,localError

    ENTERS("DistributedMatrix_LibraryTypeSet",err,error,*998)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*998)
    IF(distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has already been finished.",err,error,*998)
     
    oldLibraryType=distributedMatrix%libraryType
    IF(libraryType/=oldLibraryType) THEN
      !Initialise the new library type
      SELECT CASE(libraryType)
      CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
        CALL DistributedMatrix_CMISSInitialise(distributedMatrix,err,error,*999)
      CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
        CALL DistributedMatrix_PETScInitialise(distributedMatrix,err,error,*999)
      CASE DEFAULT
        localError="The library type of "//TRIM(NumberToVString(libraryType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Finalise the old library type
      SELECT CASE(oldLibraryType)
      CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
        CALL DistributedMatrix_CMISSFinalise(distributedMatrix%cmiss,err,error,*999)
      CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
        CALL DistributedMatrix_PETScFinalise(distributedMatrix%petsc,err,error,*999)
      CASE DEFAULT
        localError="The distributed matrix library type of "// &
          & TRIM(NumberToVString(oldLibraryType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      distributedMatrix%libraryType=libraryType
    ENDIF
    
    EXITS("DistributedMatrix_LibraryTypeSet")
    RETURN
999 SELECT CASE(libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      CALL DistributedMatrix_CMISSFinalise(distributedMatrix%cmiss,dummyErr,dummyError,*998)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL DistributedMatrix_PETScFinalise(distributedMatrix%petsc,dummyErr,dummyError,*998)
    END SELECT
998 ERRORSEXITS("DistributedMatrix_LibraryTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_LibraryTypeSet

  !
  !================================================================================================================================
  !

  !>Intialises a distributed matrix.
  SUBROUTINE DistributedMatrix_Initialise(rowDomainMapping,columnDomainMapping,distributedMatrix,err,error,*)

    !Argument variables
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: rowDomainMapping !<A pointer to the row domain mapping used to distribute this matrix
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: columnDomainMapping !<A pointer to the column domain mapping used to distribute this matrix
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("DistributedMatrix_Initialise",err,error,*998)

    IF(.NOT.ASSOCIATED(rowDomainMapping)) CALL FlagError("Row domain mapping is not associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(columnDomainMapping)) CALL FlagError("Column domain mapping is not associated.",err,error,*999)
    IF(ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is already associated.",err,error,*998)
       
    ALLOCATE(distributedMatrix,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocated the distributed matrix.",err,error,*999)
    distributedMatrix%matrixFinished=.FALSE.
    distributedMatrix%libraryType=0
    distributedMatrix%ghostingType=DISTRIBUTED_MATRIX_VECTOR_INCLUDE_GHOSTS_TYPE
    distributedMatrix%rowDomainMapping=>rowDomainMapping
    distributedMatrix%columnDomainMapping=>columnDomainMapping
    distributedMatrix%dataType=MATRIX_VECTOR_DP_TYPE
    NULLIFY(distributedMatrix%cmiss)
    NULLIFY(distributedMatrix%petsc)
    CALL DistributedMatrix_CMISSInitialise(distributedMatrix,err,error,*999)
    
    EXITS("DistributedMatrix_Initialise")
    RETURN
999 CALL DistributedMatrix_Finalise(distributedMatrix,dummyErr,dummyError,*998)
998 ERRORSEXITS("DistributedMatrix_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_Initialise

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

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed mtrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
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

  !>Sets/changes the number of non zeros for a distributed matrix.
  SUBROUTINE DistributedMatrix_NumberOfNonZerosSet(distributedMatrix,numberOfNonZeros,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(IN) :: numberOfNonZeros !<The number of non zeros in the matrix to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(DistributedMatrixPETScType), POINTER :: petscMatrix
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedMatrix_NumberOfNonZerosSet",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed mtrix is not associated.",err,error,*999)
    IF(distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has been finished.",err,error,*999)
    IF(numberOfNonZeros<=0) THEN
      localError="The specified number of non zeros of "//TRIM(NumberToVString(numberOfNonZeros,"*",err,error))// &
        & " is invalid. The number must be > 0."
      CALL FlagError(localError,err,error,*999)
    ENDIF
   
    SELECT CASE(distributedMatrix%libraryType)      
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_NumberOfNonZerosSet(distributedMatrix%cmiss%matrix,numberOfNonZeros,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      NULLIFY(petscMatrix)
      CALL DistributedMatrix_PETScMatrixGet(distributedMatrix,petscMatrix,err,error,*999)
      distributedMatrix%petsc%numberOfNonZeros=numberOfNonZeros
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("DistributedMatrix_NumberOfNonZerosSet")
    RETURN
999 ERRORSEXITS("DistributedMatrix_NumberOfNonZerosSet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_NumberOfNonZerosSet

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

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix is not finished.",err,error,*999)
    
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

  !>Sets/changes the LIST STRUCTURE for a distributed matrix.
  SUBROUTINE DistributedMatrix_LinkListSet(distributedMatrix,list,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    TYPE(LinkedList),POINTER :: list(:) 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(DistributedMatrixPETScType), POINTER :: petscMatrix
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedMatrix_LinkListSet",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed mtrix is not associated.",err,error,*999)
    IF(distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_LinkListSet(cmissMatrix%matrix,list,err,error,*999)
      !cmissMatrix%list=list 
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      NULLIFY(petscMatrix)
      CALL DistributedMatrix_PETScMatrixGet(distributedMatrix,petscMatrix,err,error,*999)
      petscMatrix%list=>list
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("DistributedMatrix_LinkListSet")
    RETURN
999 ERRORSEXITS("DistributedMatrix_LinkListSet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_LinkListSet

  !
  !================================================================================================================================
  !

  !>Gets the LINKLIST STURUCTURE for a distributed matrix.
  SUBROUTINE DistributedMatrix_LinkListGet(distributedMatrix,list,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    type(LinkedList),pointer :: list(:) 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(DistributedMatrixPETScType), POINTER :: petscMatrix
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedMatrix_LinkListGet",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix is not finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      !list=>cmissMatrix%list
      CALL Matrix_LinkListGet(cmissMatrix%matrix,list,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      NULLIFY(petscMatrix)
      CALL DistributedMatrix_PETScMatrixGet(distributedMatrix,petscMatrix,err,error,*999)
      list=>petscMatrix%list
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("DistributedMatrix_LinkListGet")
    RETURN
999 ERRORSEXITS("DistributedMatrix_LinkListGet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_LinkListGet

  !
  !================================================================================================================================
  !
  !>Outputs a distributed matrix.
  SUBROUTINE DistributedMatrix_Output(id,distributedMatrix,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: globalRow,numberOfColumns,rowIdx
    INTEGER(INTG), ALLOCATABLE :: columns(:)
    REAL(DP), ALLOCATABLE :: values(:)
    CHARACTER(LEN=9) :: rowString
    CHARACTER(LEN=39) :: initialString
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(DistributedMatrixPETScType), POINTER :: petscMatrix
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: rowMapping
    TYPE(VARYING_STRING) :: localError
     
    ENTERS("DistributedMatrix_Output",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("Distributed matrix has not been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_Output(id,cmissMatrix%matrix,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      NULLIFY(petscMatrix)
      CALL DistributedMatrix_PETScMatrixGet(distributedMatrix,petscMatrix,err,error,*999)
      NULLIFY(rowMapping)
      CALL DistributedMatrix_RowMappingGet(distributedMatrix,rowMapping,err,error,*999)
      ALLOCATE(columns(petscMatrix%maximumColumnIndicesPerRow),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate columns.",err,error,*999)
      ALLOCATE(values(petscMatrix%maximumColumnIndicesPerRow),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate values.",err,error,*999)
      DO rowIdx=1,petscMatrix%m
        globalRow=rowMapping%LOCAL_TO_GLOBAL_MAP(rowIdx)
        IF(petscMatrix%useOverrideMatrix) THEN
          CALL Petsc_MatGetRow(petscMatrix%overrideMatrix,globalRow-1,numberOfColumns,columns,values,err,error,*999)
        ELSE
          CALL Petsc_MatGetRow(petscMatrix%matrix,globalRow-1,numberOfColumns,columns,values,err,error,*999)
        ENDIF
        rowString=NumberToCharacter(rowIdx,"I9",err,error)
        initialString='("Matrix('//rowString//',:):",8(X,E13.6))'
        CALL WriteStringVector(id,1,1,numberOfColumns,8,8,values,initialString,'(20X,8(X,E13.6))',err,error,*999)
        IF(petscMatrix%useOverrideMatrix) THEN
          CALL Petsc_MatRestoreRow(petscMatrix%overrideMatrix,globalRow-1,numberOfColumns,columns,values,err,error,*999)
        ELSE
          CALL Petsc_MatRestoreRow(petscMatrix%matrix,globalRow-1,numberOfColumns,columns,values,err,error,*999)
        ENDIF
      ENDDO !rowIdx
      IF(ALLOCATED(values)) DEALLOCATE(values)
      IF(ALLOCATED(columns)) DEALLOCATE(columns)
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_Output")
    RETURN
999 IF(ALLOCATED(values)) DEALLOCATE(values)
    IF(ALLOCATED(columns)) DEALLOCATE(columns)
    ERRORSEXITS("DistributedMatrix_Output",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_Output

  !
  !================================================================================================================================
  !

  !>Sets the override matrix for a distributed matrix.
  SUBROUTINE DistributedMatrix_OverrideSetOn(distributedMatrix,overrideMatrix,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix to override
    TYPE(PetscMatType), INTENT(IN) :: overrideMatrix !<The override matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixPETScType), POINTER :: petscMatrix
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedMatrix_OverrideSetOn",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("Distributed matrix has not been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)          
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      NULLIFY(petscMatrix)
      CALL DistributedMatrix_PETScMatrixGet(distributedMatrix,petscMatrix,err,error,*999)
      IF(.NOT.petscMatrix%useOverrideMatrix) CALL FlagError("PETSc override matrix is already set.",err,error,*999)
      petscMatrix%useOverrideMatrix=.TRUE.
      petscMatrix%overrideMatrix=overrideMatrix
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
     
    EXITS("DistributedMatrix_OverrideSetOn")
    RETURN
999 ERRORSEXITS("DistributedMatrix_OverrideSetOn",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_OverrideSetOn

  !
  !================================================================================================================================
  !

  !>Turns off the override matrix for a distributed matrix.
  SUBROUTINE DistributedMatrix_OverrideSetOff(distributedMatrix,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix to override
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixPETScType), POINTER :: petscMatrix
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedMatrix_OverrideSetOff",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("Distributed matrix has not been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)          
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      NULLIFY(petscMatrix)
      CALL DistributedMatrix_PETScMatrixGet(distributedMatrix,petscMatrix,err,error,*999)
      IF(petscMatrix%useOverrideMatrix) CALL FlagError("PETSc override matrix is not set.",err,error,*999)
      distributedMatrix%petsc%useOverrideMatrix=.FALSE.
      CALL Petsc_MatInitialise(distributedMatrix%petsc%overrideMatrix,err,error,*999)
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_OverrideSetOff")
    RETURN
999 ERRORSEXITS("DistributedMatrix_OverrideSetOff",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_OverrideSetOff

  !
  !================================================================================================================================
  !

  !>Finishes the creation of a PETSc distributed matrix.
  SUBROUTINE DistributedMatrix_PETScCreateFinish(petscMatrix,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixPETScType), POINTER :: petscMatrix !<A pointer to the distributed PETSc matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,rowIdx
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: rowDomainMapping,columnDomainMapping
    TYPE(VARYING_STRING) :: dummyError,localError
    
    ENTERS("DistributedMatrix_PETScCreateFinish",err,error,*999)

    IF(.NOT.ASSOCIATED(petscMatrix)) CALL FlagError("Distributed matrix PETSc is not associated.",err,error,*998)

    NULLIFY(distributedMatrix)
    CALL DistributedMatrixPETSc_DistributedMatrixGet(petscMatrix,distributedMatrix,err,error,*999)
    NULLIFY(rowDomainMapping)
    CALL DistributedMatrix_RowMappingGet(distributedMatrix,rowDomainMapping,err,error,*999)
    NULLIFY(columnDomainMapping)
    CALL DistributedMatrix_ColumnMappingGet(distributedMatrix,columnDomainMapping,err,error,*999)
    
    SELECT CASE(petscMatrix%storageType)
    CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
      petscMatrix%numberOfNonZeros=petscMatrix%M*petscMatrix%globalN
      petscMatrix%maximumColumnIndicesPerRow=petscMatrix%globalN
      petscMatrix%dataSize=petscMatrix%numberOfNonZeros
      !Set up the Local to Global mappings
      ALLOCATE(petscMatrix%globalRowNumbers(petscMatrix%m),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate global row numbers for PETSc distributed matrix.",err,error,*999)
      DO rowIdx=1,petscMatrix%m
        petscMatrix%globalRowNumbers(rowIdx)=rowDomainMapping%LOCAL_TO_GLOBAL_MAP(rowIdx)-1 !PETSc uses 0 based indexing
      ENDDO !rowIdx
      !Set up the matrix
      ALLOCATE(petscMatrix%dataDP(petscMatrix%dataSize),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate PETSc matrix data.",err,error,*999)
      CALL Petsc_MatCreateDense(computationalEnvironment%mpiCommunicator,petscMatrix%m,petscMatrix%n, &
        & petscMatrix%globalM,petscMatrix%globalN,petscMatrix%dataDP,petscMatrix%matrix,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
      petscMatrix%numberOfNonZeros=petscMatrix%m
      petscMatrix%maximumColumnIndicesPerRow=1
      petscMatrix%dataSize=petscMatrix%numberOfNonZeros
      !Set up the Local to Global mappings
      ALLOCATE(petscMatrix%globalRowNumbers(petscMatrix%m),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate global row numbers for PETSc distributed matrix.",err,error,*999)
      DO rowIdx=1,petscMatrix%m
        petscMatrix%globalRowNumbers(rowIdx)=rowDomainMapping%LOCAL_TO_GLOBAL_MAP(rowIdx)-1 !PETSc uses 0 based indexing
      ENDDO !rowIdx
      !Set up the matrix
      ALLOCATE(petscMatrix%diagonalNumberOfNonZeros(petscMatrix%n),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate diagonal number of non zeros.",err,error,*999)
      ALLOCATE(petscMatrix%offdiagonalNumberOfNonZeros(petscMatrix%n),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate off diagonal number of non zeros.",err,error,*999)
      petscMatrix%diagonalNumberOfNonZeros=1
      petscMatrix%offdiagonalNumberOfNonZeros=0
      !Create the PETsc AIJ matrix
      CALL Petsc_MatCreateAIJ(computationalEnvironment%mpiCommunicator,petscMatrix%m,petscMatrix%n, &
        & petscMatrix%globalM,petscMatrix%globalN,PETSC_NULL_INTEGER,petscMatrix%diagonalNumberOfNonZeros, &
        & PETSC_NULL_INTEGER,petscMatrix%offdiagonalNumberOfNonZeros,petscMatrix%matrix,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
      CALL FlagError("Column major storage is not implemented for PETSc matrices.",err,error,*999)
    CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
      CALL FlagError("Row major storage is not implemented for PETSc matrices.",err,error,*999)
    CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
      IF(.NOT.ALLOCATED(petscMatrix%diagonalNumberOfNonZeros)) &
        & CALL FlagError("Matrix diagonal storage locations have not been set.",err,error,*999)
      IF(.NOT.ALLOCATED(petscMatrix%offdiagonalNumberOfNonZeros)) &
        & CALL FlagError("Matrix off diagonal storage locations have not been set.",err,error,*999)
      !Create the PETSc AIJ matrix
      CALL Petsc_MatCreateAIJ(computationalEnvironment%mpiCommunicator,petscMatrix%m,petscMatrix%n, &
        & petscMatrix%globalM,petscMatrix%globalN,PETSC_NULL_INTEGER,petscMatrix%diagonalNumberOfNonZeros, &
        & PETSC_NULL_INTEGER,petscMatrix%offdiagonalNumberOfNonZeros,petscMatrix%matrix,err,error,*999)
      !Set matrix options
      CALL Petsc_MatSetOption(petscMatrix%matrix,PETSC_MAT_NEW_NONZERO_LOCATION_ERR,.TRUE.,err,error,*999)
      CALL Petsc_MatSetOption(petscMatrix%matrix,PETSC_MAT_NEW_NONZERO_ALLOCATION_ERR,.TRUE.,err,error,*999)
      CALL Petsc_MatSetOption(petscMatrix%matrix,PETSC_MAT_UNUSED_NONZERO_LOCATION_ERR,.TRUE.,err,error,*999)
      !Set up the Local to Global mappings
      ALLOCATE(petscMatrix%globalRowNumbers(petscMatrix%m),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate global row numbers for PETSc distributed matrix.",err,error,*999)
      petscMatrix%maximumColumnIndicesPerRow=0
      DO rowIdx=1,petscMatrix%m
        petscMatrix%globalRowNumbers(rowIdx)=rowDomainMapping%LOCAL_TO_GLOBAL_MAP(rowIdx)-1 !PETSc uses 0 based indexing
        IF((petscMatrix%rowIndices(rowIdx+1)-petscMatrix%rowIndices(rowIdx))>petscMatrix%maximumColumnIndicesPerRow) &
          & petscMatrix%maximumColumnIndicesPerRow=petscMatrix%rowIndices(rowIdx+1)-petscMatrix%rowIndices(rowIdx)
      ENDDO !rowIdx
    CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
      CALL FlagError("Compressed column storage is not implemented for PETSc matrices.",err,error,*999)
    CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
      CALL FlagError("Row column storage is not implemented for PETSc matrices.",err,error,*999)
    CASE DEFAULT
      localError="The PETSc matrix storage type of "//TRIM(NumberToVString(petscMatrix%storageType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    SELECT CASE(petscMatrix%symmetryType)
    CASE(DISTRIBUTED_MATRIX_SYMMETRIC_TYPE)
      CALL Petsc_MatSetOption(petscMatrix%matrix,PETSC_MAT_SYMMETRIC,.TRUE.,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_HERMITIAN_TYPE)
      CALL Petsc_MatSetOption(petscMatrix%matrix,PETSC_MAT_HERMITIAN,.TRUE.,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_SKEW_SYMMETRIC_TYPE)
      CALL FlagError("Skew symmetric matrices are not implemented for PETSc matrices.",err,error,*999)
    CASE(DISTRIBUTED_MATRIX_UNSYMMETRIC_TYPE)
      CALL Petsc_MatSetOption(petscMatrix%matrix,PETSC_MAT_SYMMETRIC,.FALSE.,err,error,*999)
    CASE DEFAULT
      localError="The PETSc matrix symmetry type of "//TRIM(NumberToVString(petscMatrix%symmetryType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_PETScCreateFinish")
    RETURN
999 CALL DistributedMatrix_PETScFinalise(petscMatrix,dummyErr,dummyError,*998)
998 ERRORSEXITS("DistributedMatrix_PETScCreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_PETScCreateFinish

  !
  !================================================================================================================================
  !

  !>Finalise a PETSc distributed matrix.
  SUBROUTINE DistributedMatrix_PETScFinalise(petscMatrix,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixPETScType), POINTER :: petscMatrix !<A pointer to the PETSc distributed matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("DistributedMatrix_PETScFinalise",err,error,*999)

    IF(ASSOCIATED(petscMatrix)) THEN
      IF(ALLOCATED(petscMatrix%diagonalNumberOfNonZeros)) DEALLOCATE(petscMatrix%diagonalNumberOfNonZeros)
      IF(ALLOCATED(petscMatrix%offdiagonalNumberOfNonZeros)) DEALLOCATE(petscMatrix%offdiagonalNumberOfNonZeros)
      IF(ALLOCATED(petscMatrix%rowIndices)) DEALLOCATE(petscMatrix%rowIndices)
      IF(ALLOCATED(petscMatrix%columnIndices)) DEALLOCATE(petscMatrix%columnIndices)
      IF(ALLOCATED(petscMatrix%globalRowNumbers)) DEALLOCATE(petscMatrix%globalRowNumbers)
      CALL Petsc_MatFinalise(petscMatrix%matrix,err,error,*999)
      CALL Petsc_MatFinalise(petscMatrix%overrideMatrix,err,error,*999)
      DEALLOCATE(petscMatrix)
    ENDIF
    
    EXITS("DistributedMatrix_PETScFinalise")
    RETURN
999 ERRORSEXITS("DistributedMatrix_PETScFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_PETScFinalise

  !
  !================================================================================================================================
  !

  !>Intialises a PETSc distributed matrix.
  SUBROUTINE DistributedMatrix_PETScInitialise(distributedMatrix,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: rowDomainMapping,columnDomainMapping
    TYPE(VARYING_STRING) :: dummyError,localError
    
    ENTERS("DistributedMatrix_PETScInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*998)
    IF(ASSOCIATED(distributedMatrix%petsc))  &
      & CALL FlagError("PETSc is already associated for this distributed matrix",err,error,*998)

    NULLIFY(rowDomainMapping)
    CALL DistributedMatrix_RowMappingGet(distributedMatrix,rowDomainMapping,err,error,*999)
    NULLIFY(columnDomainMapping)
    CALL DistributedMatrix_ColumnMappingGet(distributedMatrix,columnDomainMapping,err,error,*999)
    
    ALLOCATE(distributedMatrix%petsc,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate PETSc distributed matrix.",err,error,*999)
    distributedMatrix%petsc%distributedMatrix=>distributedMatrix
    distributedMatrix%libraryType=DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE
    !Set the defaults          
    SELECT CASE(distributedMatrix%ghostingType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_INCLUDE_GHOSTS_TYPE)
      distributedMatrix%petsc%m=rowDomainMapping%TOTAL_NUMBER_OF_LOCAL            
    CASE(DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE)
      distributedMatrix%petsc%m=rowDomainMapping%NUMBER_OF_LOCAL
    CASE DEFAULT
      localError="The distributed matrix ghosting type of "// &
        & TRIM(NumberToVString(distributedMatrix%ghostingType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    distributedMatrix%petsc%n=columnDomainMapping%TOTAL_NUMBER_OF_LOCAL
    distributedMatrix%petsc%globalM=rowDomainMapping%NUMBER_OF_GLOBAL
    distributedMatrix%petsc%globalN=columnDomainMapping%NUMBER_OF_GLOBAL
    distributedMatrix%petsc%storageType=DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE
    distributedMatrix%petsc%symmetryType=DISTRIBUTED_MATRIX_SYMMETRIC_TYPE  !Should this be unsymmetric???
    distributedMatrix%petsc%dataSize=0
    distributedMatrix%petsc%maximumColumnIndicesPerRow=0
    distributedMatrix%petsc%useOverrideMatrix=.FALSE.
    CALL Petsc_MatInitialise(distributedMatrix%petsc%matrix,err,error,*999)
    CALL Petsc_MatInitialise(distributedMatrix%petsc%overrideMatrix,err,error,*999)
   
    EXITS("DistributedMatrix_PETScInitialise")
    RETURN
999 IF(ASSOCIATED(distributedMatrix%petsc)) &
      & CALL DistributedMatrix_PETScFinalise(distributedMatrix%petsc,dummyErr,dummyError,*998)
998 ERRORSEXITS("DistributedMatrix_PETScInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_PETScInitialise

  !
  !================================================================================================================================
  !

  !>Gets the storage locations (sparsity pattern) for a distributed matrix.
  SUBROUTINE DistributedMatrix_StorageLocationsGet(distributedMatrix,rowIndices,columnIndices,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), POINTER :: rowIndices(:) !<rowIndices(i). On return, the i'th row index of the matrix storage locations
    INTEGER(INTG), POINTER :: columnIndices(:) !<columnIndices(i). On return, the i'th column index of the matrix storage locations
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(DistributedMatrixPETScType), POINTER :: petscMatrix
    TYPE(VARYING_STRING) :: localError
     
    ENTERS("DistributedMatrix_StorageLocationsGet",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_StorageLocationsGet(cmissMatrix%matrix,rowIndices,columnIndices,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      NULLIFY(petscMatrix)
      CALL DistributedMatrix_PETScMatrixGet(distributedMatrix,petscMatrix,err,error,*999)
      SELECT CASE(petscMatrix%storageType)
      CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
        CALL FlagError("Cannot get matrix locations for a block storage matrix.",err,error,*999)
      CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
        CALL FlagError("Diagonal storage is not implemented for PETSc matrices.",err,error,*999)
      CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
        CALL FlagError("Column major storage is not implemented for PETSc matrices.",err,error,*999)
      CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
        CALL FlagError("Row major storage is not implemented for PETSc matrices.",err,error,*999)
      CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
        rowIndices=>petscMatrix%rowIndices
        columnIndices=>petscMatrix%columnIndices
      CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
        CALL FlagError("Compressed column storage is not implemented for PETSc matrices.",err,error,*999)
      CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
        CALL FlagError("Row column storage is not implemented for PETSc matrices.",err,error,*999)
      CASE DEFAULT
        localError="The matrix storage type of "// &
          & TRIM(NumberToVString(petscMatrix%storageType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_StorageLocationsGet")
    RETURN
999 ERRORSEXITS("DistributedMatrix_StorageLocationsGet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_StorageLocationsGet

  !
  !================================================================================================================================
  !

  !>Sets the storage locations (sparsity pattern) in a distributed matrix to that specified by the row and column indices.
  SUBROUTINE DistributedMatrix_StorageLocationsSet(distributedMatrix,rowIndices,columnIndices,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The i'th row index of the matrix storage locations
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(i). The i'th column index of the matrix storage locations
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: column,columnIdx,globalRowStart,globalRowFinish,rowIdx
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(DistributedMatrixPETScType), POINTER :: petscMatrix
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: rowDomainMapping,columnDomainMapping
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedMatrix_StorageLocationsSet",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has been finished.",err,error,*999)

    NULLIFY(rowDomainMapping)
    CALL DistributedMatrix_RowMappingGet(distributedMatrix,rowDomainMapping,err,error,*999)
    NULLIFY(columnDomainMapping)
    CALL DistributedMatrix_ColumnMappingGet(distributedMatrix,columnDomainMapping,err,error,*999)
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_StorageLocationsSet(cmissMatrix%matrix,rowIndices,columnIndices,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      NULLIFY(petscMatrix)
      CALL DistributedMatrix_PETScMatrixGet(distributedMatrix,petscMatrix,err,error,*999)
      SELECT CASE(petscMatrix%storageType)
      CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
        !Do nothing
      CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
        CALL FlagError("Diagonal storage is not implemented for PETSc matrices.",err,error,*999)
      CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
        CALL FlagError("Column major storage is not implemented for PETSc matrices.",err,error,*999)
      CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
        CALL FlagError("Row major storage is not implemented for PETSc matrices.",err,error,*999)
      CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
        IF(SIZE(rowIndices,1)/=petscMatrix%m+1) THEN
          localError="The supplied number of row indices of "// &
            & TRIM(NumberToVString(SIZE(rowIndices,1),"*",err,error))// &
            & " does not match the number of rows in the matrix + 1 of "// &
            & TRIM(NumberToVString(petscMatrix%m+1,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        IF(SIZE(columnIndices,1)/=petscMatrix%numberOfNonZeros) THEN
          localError="The supplied number of column indices of "// &
            & TRIM(NumberToVString(SIZE(columnIndices,1),"*",err,error))// &
            & " does not match the number of non-zeros in the matrix of "// &
            & TRIM(NumberToVString(petscMatrix%numberOfNonZeros,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        IF(rowIndices(1)/=1) THEN
          localError="Invalid row indices. The first row index of "// &
            & TRIM(NumberToVString(rowIndices(1),"*",err,error))//" does not equal 1."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        !Check the row indicies are correct
        IF(rowIndices(petscMatrix%m+1)/=petscMatrix%numberOfNonZeros+1) THEN
          localError="Invalid row indices. The last row index of "// &
            & TRIM(NumberToVString(rowIndices(petscMatrix%m+1),"*",err,error))// &
            & " does not equal the number of non-zeros + 1 of "// &
            & TRIM(NumberToVString(petscMatrix%numberOfNonZeros+1,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        DO rowIdx=2,petscMatrix%m+1
          IF(rowIndices(rowIdx)<rowIndices(rowIdx-1)) THEN
            localError="Invalid row indices. Row "//TRIM(NumberToVString(rowIdx,"*",err,error))// &
              & " index number "//TRIM(NumberToVString(rowIndices(rowIdx),"*",err,error))// &
              & " is less than row "//TRIM(NumberToVString(rowIdx-1,"*",err,error))//" index number "// &
              & TRIM(NumberToVString(rowIndices(rowIdx-1),"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ENDDO !rowIdx
        !Allocate the PETSc sparsity storage arrays
        ALLOCATE(petscMatrix%diagonalNumberOfNonZeros(petscMatrix%m),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate PETSc matrix diagonal number of non zeros.",err,error,*999)
        petscMatrix%diagonalNumberOfNonZeros=0
        ALLOCATE(petscMatrix%offdiagonalNumberOfNonZeros(petscMatrix%m),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate PETSc matrix off diagonal number of non zeros.", &
          & err,error,*999)
        petscMatrix%offdiagonalNumberOfNonZeros=0
        ALLOCATE(petscMatrix%rowIndices(petscMatrix%m+1),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate PETSc matrix row indices.",err,error,*999)
        petscMatrix%rowIndices(1:petscMatrix%m+1)=rowIndices(1:petscMatrix%m+1)
        ALLOCATE(petscMatrix%columnIndices(petscMatrix%numberOfNonZeros),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate PETSc matrix column indices.",err,error,*999)
        petscMatrix%columnIndices(1:petscMatrix%numberOfNonZeros)=columnIndices(1:petscMatrix%numberOfNonZeros)
        !Check the column indices are correct and calculate number of diagonal and off-diagonal columns
        globalRowStart=rowDomainMapping%LOCAL_TO_GLOBAL_MAP(1)
        globalRowFinish=rowDomainMapping%LOCAL_TO_GLOBAL_MAP(petscMatrix%M)
        DO rowIdx=1,petscMatrix%m
          DO columnIdx=rowIndices(rowIdx),rowIndices(rowIdx+1)-1
            column=columnIndices(columnIdx)
            IF(column<=0) THEN
              localError="Invalid column indices. Column index "//TRIM(NumberToVString(columnIdx,"*",err,error))// &
                & " gives column "//TRIM(NumberToVString(column,"*",err,error))//" which is less than zero."
              CALL FlagError(localError,err,error,*999)
            ENDIF
            IF(column>petscMatrix%globalN) THEN
              localError="Invalid column indices. Column index "//TRIM(NumberToVString(columnIdx,"*",err,error))// &
                & " gives column "//TRIM(NumberToVString(column,"*",err,error))// &
                & " which is greater than the number of columns of "// &
                & TRIM(NumberToVString(petscMatrix%globalN,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
            IF(column>=globalRowStart.AND.column<=globalRowFinish) THEN
              petscMatrix%diagonalNumberOfNonZeros(rowIdx)=petscMatrix%diagonalNumberOfNonZeros(rowIdx)+1
            ELSE
              petscMatrix%offdiagonalNumberOfNonZeros(rowIdx)=petscMatrix%offdiagonalNumberOfNonZeros(rowIdx)+1
            ENDIF
          ENDDO !columnIdx
          !Enforce a place for the diagonal entry.
          IF(petscMatrix%diagonalNumberOfNonZeros(rowIdx)==0) petscMatrix%diagonalNumberOfNonZeros(rowIdx)=1
        ENDDO !rowIdx
        IF(diagnostics3) THEN
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"PETSc distributed matrix sparsity:",err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Storage type = ",petscMatrix%storageType,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  m = ",petscMatrix%m,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  n = ",petscMatrix%n,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Global m = ",petscMatrix%globalM,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Global n = ",petscMatrix%globalN,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of non zeros = ",petscMatrix%numberOfNonZeros,err,error,*999)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,petscMatrix%m,8,8,petscMatrix%diagonalNumberOfNonZeros, &
            & '("  Diagonal number non zeros     :",8(X,I10))','(33X,8(X,I10))',err,error,*999)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,petscMatrix%m,8,8,petscMatrix%offdiagonalNumberOfNonZeros, &
            & '("  Off-diagonal number non zeros :",8(X,I10))','(33X,8(X,I10))',err,error,*999)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,petscMatrix%m+1,8,8,petscMatrix%rowIndices, &
            & '("  Row indices                   :",8(X,I10))','(33X,8(X,I10))',err,error,*999)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,petscMatrix%numberOfNonZeros,8,8,petscMatrix%columnIndices, &
            & '("  Column indices                :",8(X,I10))','(33X,8(X,I10))',err,error,*999)
        ENDIF
      CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
        CALL FlagError("Compressed column storage is not implemented for PETSc matrices.",err,error,*999)
      CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
        CALL FlagError("Row column storage is not implemented for PETSc matrices.",err,error,*999)
      CASE DEFAULT
        localError="The specified matrix storage type of "// &
          & TRIM(NumberToVString(petscMatrix%storageType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
   
    EXITS("DistributedMatrix_StorageLocationsSet")
    RETURN
999 IF(ASSOCIATED(petscMatrix)) THEN
      IF(ALLOCATED(petscMatrix%diagonalNumberOfNonZeros)) DEALLOCATE(petscMatrix%diagonalNumberOfNonZeros)
      IF(ALLOCATED(petscMatrix%offdiagonalNumberOfNonZeros)) DEALLOCATE(petscMatrix%offdiagonalNumberOfNonZeros)
    ENDIF
    ERRORSEXITS("DistributedMatrix_StorageLocationsSet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_StorageLocationsSet

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

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
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

  !>Sets/changes the storage type of a distributed matrix.
  SUBROUTINE DistributedMatrix_StorageTypeSet(distributedMatrix,storageType,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(IN) :: storageType !<The storage (sparsity) type to set. \see DISTRIBUTED_MATRIX_VECTOR_StorageTypes,DISTRIBUTED_MATRIX_VECTOR
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(DistributedMatrixPETScType), POINTER :: petscMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedMatrix_StorageTypeSet",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)  
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_StorageTypeSet(cmissMatrix%matrix,storageType,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      NULLIFY(petscMatrix)
      CALL DistributedMatrix_PETScMatrixGet(distributedMatrix,petscMatrix,err,error,*999)
      SELECT CASE(storageType)
      CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
        petscMatrix%storageType=DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE
      CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
        petscMatrix%storageType=DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE
      CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
        CALL FlagError("Column major storage is not implemented for PETSc matrices.",err,error,*999)
      CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
        CALL FlagError("Row major storage is not implemented for PETSc matrices.",err,error,*999)
      CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
        petscMatrix%storageType=DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE
      CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
        CALL FlagError("Compressed column storage is not implemented for PETSc matrices.",err,error,*999)
      CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
        CALL FlagError("Row column storage is not implemented for PETSc matrices.",err,error,*999)
      CASE DEFAULT
        localError="The specified matrix storage type of "//TRIM(NumberToVString(storageType,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_StorageTypeSet")
    RETURN
999 ERRORSEXITS("DistributedMatrix_StorageTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_StorageTypeSet

  !
  !================================================================================================================================
  !

  !>Gets the symetry type of a distributed matrix.
  SUBROUTINE DistributedMatrix_SymmetryTypeGet(distributedMatrix,symmetryType,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix to get the symmetry type for
    INTEGER(INTG), INTENT(OUT) :: symmetryType !<On return, the symmetry type of the distributed matrix. \see DISTRIBUTED_MATRIX_VECTOR_SymmetryTypes,DISTRIBUTED_MATRIX_VECTOR
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(DistributedMatrixPETScType), POINTER :: petscMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedMatrix_SymmetryTypeGet",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
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

  !>Sets/changes the symmetry type of a distributed matrix.
  SUBROUTINE DistributedMatrix_SymmetryTypeSet(distributedMatrix,symmetryType,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix to set the symmetry type for
    INTEGER(INTG), INTENT(IN) :: symmetryType !<The symmetry type to set. \see DISTRIBUTED_MATRIX_VECTOR_SymmetryTypes,DISTRIBUTED_MATRIX_VECTOR
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(DistributedMatrixPETScType), POINTER :: petscMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedMatrix_SymmetryTypeSet",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)      
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_SymmetryTypeSet(cmissMatrix%matrix,symmetryType,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      NULLIFY(petscMatrix)
      CALL DistributedMatrix_PETScMatrixGet(distributedMatrix,petscMatrix,err,error,*999)
      SELECT CASE(symmetryType)
      CASE(DISTRIBUTED_MATRIX_SYMMETRIC_TYPE)
        petscMatrix%symmetryType=DISTRIBUTED_MATRIX_SYMMETRIC_TYPE
      CASE(DISTRIBUTED_MATRIX_HERMITIAN_TYPE)
        IF(distributedMatrix%dataType==DISTRIBUTED_MATRIX_VECTOR_SPC_TYPE.OR. &
          & distributedMatrix%dataType==DISTRIBUTED_MATRIX_VECTOR_DPC_TYPE) THEN
          petscMatrix%symmetryType=DISTRIBUTED_MATRIX_HERMITIAN_TYPE
        ELSE
          CALL FlagError( &
            & "Cannot set the distributed matrix symmetry type to Hermitian as the matrix does not have a complex data type.", &
            & err,error,*999)
        ENDIF
      CASE(DISTRIBUTED_MATRIX_SKEW_SYMMETRIC_TYPE)
        CALL FlagError("Skew symmetric matrices are not implemented for PETSc matrices.",err,error,*999)
        petscMatrix%symmetryType=DISTRIBUTED_MATRIX_SKEW_SYMMETRIC_TYPE
      CASE(DISTRIBUTED_MATRIX_UNSYMMETRIC_TYPE)
        petscMatrix%symmetryType=DISTRIBUTED_MATRIX_UNSYMMETRIC_TYPE        
      CASE DEFAULT
        localError="The specified matrix symmetry type of "//TRIM(NumberToVString(symmetryType,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_SymmetryTypeSet")
    RETURN
999 ERRORSEXITS("DistributedMatrix_SymmetryTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_SymmetryTypeSet

  !
  !================================================================================================================================
  !

  !>Finishes the update procedure for a distributed matrix. This routine will wait until all transfers have completed!
  SUBROUTINE DistributedMatrix_UpdateFinish(distributedMatrix,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixPETScType), POINTER :: petscMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedMatrix_UpdateFinish",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      !Do nothing for now.
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      NULLIFY(petscMatrix)
      CALL DistributedMatrix_PETScMatrixGet(distributedMatrix,petscMatrix,err,error,*999)
      IF(petscMatrix%useOverrideMatrix) THEN
        CALL Petsc_MatAssemblyEnd(petscMatrix%overrideMatrix,PETSC_MAT_FINAL_ASSEMBLY,err,error,*999)
      ELSE
        CALL Petsc_MatAssemblyEnd(petscMatrix%matrix,PETSC_MAT_FINAL_ASSEMBLY,err,error,*999)
      ENDIF
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_UpdateFinish")
    RETURN
999 ERRORSEXITS("DistributedMatrix_UpdateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_UpdateFinish

  !
  !================================================================================================================================
  !

  !>Tests to see if a distributed matrix update has finised.
  SUBROUTINE DistributedMatrix_UpdateIsFinished(distributedMatrix,isFinished,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    LOGICAL, INTENT(OUT) :: isFinished !<On return, isFinished will be .TRUE. if the distributed matrix update has finished or .FALSE. if it has not
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("DistributedMatrix_UpdateIsFinished",err,error,*999)

    isFinished=.FALSE.
    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
    !Do nothing for now.
    isFinished=.TRUE.
    
    EXITS("DistributedMatrix_UpdateIsFinished")
    RETURN
999 ERRORSEXITS("DistributedMatrix_UpdateIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_UpdateIsFinished

  !
  !================================================================================================================================
  !

  !>Waits until a distributed matrix update has finised.
  SUBROUTINE DistributedMatrix_UpdateWaitFinished(distributedMatrix,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("DistributedMatrix_UpdateWaitFinished",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
    !Do nothing for now.
    
    EXITS("DistributedMatrix_UpdateWaitFinished")
    RETURN
999 ERRORSEXITS("DistributedMatrix_UpdateWaitFinished",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_UpdateWaitFinished

  !
  !================================================================================================================================
  !

  !>Starts the update procedure for a distributed matrix.
  SUBROUTINE DistributedMatrix_UpdateStart(distributedMatrix,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixPETScType), POINTER :: petscMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedMatrix_UpdateStart",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished)  CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      !Do nothing for now.
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      NULLIFY(petscMatrix)
      CALL DistributedMatrix_PETScMatrixGet(distributedMatrix,petscMatrix,err,error,*999)
      IF(petscMatrix%useOverrideMatrix) THEN
        CALL Petsc_MatAssemblyBegin(petscMatrix%overrideMatrix,PETSC_MAT_FINAL_ASSEMBLY,err,error,*999)
      ELSE
        CALL Petsc_MatAssemblyBegin(petscMatrix%matrix,PETSC_MAT_FINAL_ASSEMBLY,err,error,*999)
      ENDIF
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("DistributedMatrix_UpdateStart")
    RETURN
999 ERRORSEXITS("DistributedMatrix_UpdateStart",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_UpdateStart

  !
  !================================================================================================================================
  !

  !>Adds values to a distributed integer matrix.
  SUBROUTINE DistributedMatrix_ValuesAddIntg(distributedMatrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The i'th row index to add
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(i). The i'th column index to add
    INTEGER(INTG), INTENT(IN) :: values(:) !<values(i). The i'th value to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedMatrix_ValuesAddIntg",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_ValuesAdd(cmissMatrix%matrix,rowIndices,columnIndices,values,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Adding values to an integer PETSc distributed matrix is not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_ValuesAddIntg")
    RETURN
999 ERRORSEXITS("DistributedMatrix_ValuesAddIntg",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_ValuesAddIntg

  !
  !================================================================================================================================
  !

  !>Adds one value to a distributed integer matrix.
  SUBROUTINE DistributedMatrix_ValuesAddIntg1(distributedMatrix,rowIndex,columnIndex,value,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(IN) :: rowIndex !<The row index to add a value to 
    INTEGER(INTG), INTENT(IN) :: columnIndex !<The column index to add a value to
    INTEGER(INTG), INTENT(IN) :: value !<The value to add at the specified row and column
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(DistributedMatrixPETScType), POINTER :: petscMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedMatrix_ValuesAddIntg1",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_ValuesAdd(cmissMatrix%matrix,rowIndex,columnIndex,VALUE,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Adding values to an integer PETSc distributed matrix is not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_ValuesAddIntg1")
    RETURN
999 ERRORSEXITS("DistributedMatrix_ValuesAddIntg1",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_ValuesAddIntg1

  !
  !================================================================================================================================
  !

  !>Adds a matrix of values to a distributed integer matrix.
  SUBROUTINE DistributedMatrix_ValuesAddIntg2(distributedMatrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The ij'th row index to add
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(j). The ij'th column index to add
    INTEGER(INTG), INTENT(IN) :: values(:,:) !<values(i,j). The ij'th value to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedMatrix_ValuesAddIntg2",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_ValuesAdd(cmissMatrix%matrix,rowIndices,columnIndices,values,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Adding values to an integer PETSc distributed matrix is not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_ValuesAddIntg2")
    RETURN
999 ERRORSEXITS("DistributedMatrix_ValuesAddIntg2",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_ValuesAddIntg2

  !
  !================================================================================================================================
  !

  !>Adds values to a distributed single precision matrix.
  SUBROUTINE DistributedMatrix_ValuesAddSP(distributedMatrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The i'th row index to add
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(i). The i'th column index to add
    REAL(SP), INTENT(IN) :: values(:) !<values(i). The i'th value to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedMatrix_ValuesAddSP",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_ValuesAdd(cmissMatrix%matrix,rowIndices,columnIndices,values,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Adding values to a single precision PETSc distributed matrix is not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_ValuesAddSP")
    RETURN
999 ERRORSEXITS("DistributedMatrix_ValuesAddSP",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_ValuesAddSP

  !
  !================================================================================================================================
  !

  !>Adds one value to a distributed single precision matrix.
  SUBROUTINE DistributedMatrix_ValuesAddSP1(distributedMatrix,rowIndex,columnIndex,value,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(IN) :: rowIndex !<The row index to add a value to 
    INTEGER(INTG), INTENT(IN) :: columnIndex !<The column index to add a value to
    REAL(SP), INTENT(IN) :: value !<The value to add at the specified row and column
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedMatrix_ValuesAddSP1",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_ValuesAdd(cmissMatrix%matrix,rowIndex,columnIndex,value,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Adding values to a single precision PETSc distributed matrix is not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_ValuesAddSP1")
    RETURN
999 ERRORSEXITS("DistributedMatrix_ValuesAddSP1",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_ValuesAddSP1

  !
  !================================================================================================================================
  !

  !>Adds a matrix of values to a distributed single precision matrix.
  SUBROUTINE DistributedMatrix_ValuesAddSP2(distributedMatrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The ij'th row index to add
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(j). The ij'th column index to add
    REAL(SP), INTENT(IN) :: values(:,:) !<values(i,j). The ij'th value to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(VARYING_STRING) :: localError
   
    ENTERS("DistributedMatrix_ValuesAddSP2",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_ValuesAdd(cmissMatrix%matrix,rowIndices,columnIndices,values,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Adding values to a single precision PETSc distributed matrix is not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_ValuesAddSP2")
    RETURN
999 ERRORSEXITS("DistributedMatrix_ValuesAddSP2",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_ValuesAddSP2

  !
  !================================================================================================================================
  !

  !>Adds values to a distributed double precision matrix.
  SUBROUTINE DistributedMatrix_ValuesAddDP(distributedMatrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The i'th row index to add
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(i). The i'th column index to add
    REAL(DP), INTENT(IN) :: values(:) !<values(i). The i'th value to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: i
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(DistributedMatrixPETScType), POINTER :: petscMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedMatrix_ValuesAddDP",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_ValuesAdd(cmissMatrix%matrix,rowIndices,columnIndices,values,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      IF(SIZE(rowIndices,1)/=SIZE(values,1)) THEN
        localError="The size of the row indices array of "// &
          & TRIM(NumberToVString(SIZE(rowIndices,1),"*",err,error))// &
          & " does not conform to the size of the values array of "// &
          & TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      IF(SIZE(columnIndices,1)/=SIZE(values,1)) THEN
        localError="The size of the column indices array of "// &
          & TRIM(NumberToVString(SIZE(columnIndices,1),"*",err,error))// &
          & " does not conform to the size of the values array of "// &
          & TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      NULLIFY(petscMatrix)
      CALL DistributedMatrix_PETScMatrixGet(distributedMatrix,petscMatrix,err,error,*999)
      IF(petscMatrix%useOverrideMatrix) THEN
        DO i=1,SIZE(rowIndices,1)
          !Use global matrix row and column numbers
          CALL Petsc_MatSetValues(petscMatrix%overrideMatrix,1,petscMatrix%globalRowNumbers(rowIndices(i:i)), &
            & 1,columnIndices(i:i)-1,values(i:i),PETSC_ADD_VALUES,err,error,*999) !PETSc uses 0 indicies
        ENDDO !i
      ELSE
        DO i=1,SIZE(rowIndices,1)
          !Use global matrix row and column numbers
          CALL Petsc_MatSetValues(petscMatrix%matrix,1,petscMatrix%globalRowNumbers(rowIndices(i:i)), &
            & 1,columnIndices(i:i)-1,values(i:i),PETSC_ADD_VALUES,err,error,*999) !PETSc uses 0 indicies
        ENDDO !i
      ENDIF
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_ValuesAddDP")
    RETURN
999 ERRORSEXITS("DistributedMatrix_ValuesAddDP",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_ValuesAddDP

  !
  !================================================================================================================================
  !

  !>Adds one value to a distributed double precision matrix.
  SUBROUTINE DistributedMatrix_ValuesAddDP1(distributedMatrix,rowIndex,columnIndex,value,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(IN) :: rowIndex !<The row index to add a value to 
    INTEGER(INTG), INTENT(IN) :: columnIndex !<The column index to add a value to
    REAL(DP), INTENT(IN) :: value !<The value to add at the specified row and column
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: petscColIndex(1)
    REAL(DP) :: petscValue(1)
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(DistributedMatrixPETScType), POINTER :: petscMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedMatrix_ValuesAddDP1",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_ValuesAdd(cmissMatrix%matrix,rowIndex,columnIndex,VALUE,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      NULLIFY(petscMatrix)
      CALL DistributedMatrix_PETScMatrixGet(distributedMatrix,petscMatrix,err,error,*999)
      !Use global matrix row and column numbers
      petscColIndex(1)=columnIndex-1
      petscValue(1)=value
      IF(petscMatrix%useOverrideMatrix) THEN
        CALL Petsc_MatSetValue(petscMatrix%overrideMatrix,petscMatrix%globalRowNumbers(rowIndex),columnIndex-1, &
          & value,PETSC_ADD_VALUES,err,error,*999) !PETSc uses 0 based indices
      ELSE
        CALL Petsc_MatSetValue(petscMatrix%matrix,petscMatrix%globalRowNumbers(rowIndex),columnIndex-1, &
          & value,PETSC_ADD_VALUES,err,error,*999) !PETSc uses 0 based indices
      ENDIF
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_ValuesAddDP1")
    RETURN
999 ERRORSEXITS("DistributedMatrix_ValuesAddDP1",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_ValuesAddDP1

  !
  !================================================================================================================================
  !

  !>Adds a matrix of values to a distributed double precision matrix.
  SUBROUTINE DistributedMatrix_ValuesAddDP2(distributedMatrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The ij'th row index to add
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(j). The ij'th column index to add
    REAL(DP), INTENT(IN) :: values(:,:) !<values(i,j). The ij'th value to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: globalRowIndices(SIZE(rowIndices,1)),i
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(DistributedMatrixPETScType), POINTER :: petscMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedMatrix_ValuesAddDP2",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_ValuesAdd(cmissMatrix%matrix,rowIndices,columnIndices,values,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      IF(SIZE(rowIndices,1)/=SIZE(values,1)) THEN
        localError="The size of the row indices array of "// &
          & TRIM(NumberToVString(SIZE(rowIndices,1),"*",err,error))// &
          & " does not conform to the number of rows in the values array of "// &
          & TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      IF(SIZE(columnIndices,1)/=SIZE(values,2)) THEN                
        localError="The size of the column indices array of "// &
          & TRIM(NumberToVString(SIZE(columnIndices,1),"*",err,error))// &
          & " does not conform to the number of columns in the values array of "// &
          & TRIM(NumberToVString(SIZE(values,2),"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      NULLIFY(petscMatrix)
      CALL DistributedMatrix_PETScMatrixGet(distributedMatrix,petscMatrix,err,error,*999)
      DO i=1,SIZE(rowIndices,1)
        globalRowIndices(i)=petscMatrix%globalRowNumbers(rowIndices(i))
      ENDDO !i
      IF(petscMatrix%useOverrideMatrix) THEN
        CALL Petsc_MatSetValues(petscMatrix%overrideMatrix,SIZE(rowIndices,1),globalRowIndices, &
          & SIZE(columnIndices,1),columnIndices-1,values,PETSC_ADD_VALUES,err,error,*999) !PETSc uses 0 based indices
      ELSE
        CALL Petsc_MatSetValues(petscMatrix%matrix,SIZE(rowIndices,1),globalRowIndices, &
          & SIZE(columnIndices,1),columnIndices-1,values,PETSC_ADD_VALUES,err,error,*999) !PETSc uses 0 based indices
      ENDIF
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_ValuesAddDP2")
    RETURN
999 ERRORSEXITS("DistributedMatrix_ValuesAddDP2",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_ValuesAddDP2

  !
  !================================================================================================================================
  !

  !>Adds values to a distributed logical matrix.
  SUBROUTINE DistributedMatrix_ValuesAddL(distributedMatrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The i'th row index to add
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(i). The i'th column index to add
    LOGICAL, INTENT(IN) :: values(:) !<values(i). The i'th value to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedMatrix_ValuesAddL",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)      
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_ValuesAdd(cmissMatrix%matrix,rowIndices,columnIndices,values,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Adding values to a logical PETSc distributed matrix is not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_ValuesAddL")
    RETURN
999 ERRORSEXITS("DistributedMatrix_ValuesAddL",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_ValuesAddL

  !
  !================================================================================================================================
  !

  !>Adds one value to a distributed logical matrix.
  SUBROUTINE DistributedMatrix_ValuesAddL1(distributedMatrix,rowIndex,columnIndex,value,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(IN) :: rowIndex !<The row index to add a value to 
    INTEGER(INTG), INTENT(IN) :: columnIndex !<The column index to add a value to
    LOGICAL, INTENT(IN) :: value !<The value to add at the specified row and column
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedMatrix_ValuesAddL1",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_ValuesAdd(cmissMatrix%matrix,rowIndex,columnIndex,value,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Adding values to a logical PETSc distributed matrix is not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_ValuesAddL1")
    RETURN
999 ERRORSEXITS("DistributedMatrix_ValuesAddL1",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_ValuesAddL1

  !
  !================================================================================================================================
  !

  !>Adds a matrix of values to a distributed logical matrix.
  SUBROUTINE DistributedMatrix_ValuesAddL2(distributedMatrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The ij'th row index to add
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(j). The ij'th column index to add
    LOGICAL, INTENT(IN) :: values(:,:) !<values(i,j). The ij'th value to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedMatrix_ValuesAddL2",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_ValuesAdd(cmissMatrix%matrix,rowIndices,columnIndices,values,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Adding values to a logical PETSc distributed matrix is not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_ValuesAddL2")
    RETURN
999 ERRORSEXITS("DistributedMatrix_ValuesAddL2",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_ValuesAddL2

  !
  !================================================================================================================================
  !

  !>Gets values in a distributed integer matrix.
  SUBROUTINE DistributedMatrix_ValuesGetIntg(distributedMatrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The i'th row index to get
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(i). The i'th column index to get
    INTEGER(INTG), INTENT(OUT) :: values(:) !<values(i). On return the i'th value to get
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedMatrix_ValuesGetIntg",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)  
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_ValuesGet(cmissMatrix%matrix,rowIndices,columnIndices,values,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot get values for an integer PETSc distributed matrix.",err,error,*999)          
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
   
    EXITS("DistributedMatrix_ValuesGetIntg")
    RETURN
999 ERRORSEXITS("DistributedMatrix_ValuesGetIntg",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_ValuesGetIntg

  !
  !================================================================================================================================
  !

  !>Gets one value in a distributed integer matrix.
  SUBROUTINE DistributedMatrix_ValuesGetIntg1(distributedMatrix,rowIndex,columnIndex,value,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(IN) :: rowIndex !<The row index to get a value from
    INTEGER(INTG), INTENT(IN) :: columnIndex !<The column index to get a value from
    INTEGER(INTG), INTENT(OUT) :: value !<On return the value of the matrix at the specified row and column
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedMatrix_ValuesGetIntg1",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)      
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_ValuesGet(cmissMatrix%matrix,rowIndex,columnIndex,VALUE,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot get values for an integer PETSc distributed matrix.",err,error,*999)          
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_ValuesGetIntg1")
    RETURN
999 ERRORSEXITS("DistributedMatrix_ValuesGetIntg1",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_ValuesGetIntg1

  !
  !================================================================================================================================
  !

  !>Gets a matrix of values in a distributed integer matrix.
  SUBROUTINE DistributedMatrix_ValuesGetIntg2(distributedMatrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The ij'th row index to get
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(j). The ij'th column index to get
    INTEGER(INTG), INTENT(OUT) :: values(:,:) !<values(i,j). On return the ij'th value to get
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedMatrix_ValuesGetIntg2",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
     
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_ValuesGet(cmissMatrix%matrix,rowIndices,columnIndices,values,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot get values for an integer PETSc distributed matrix.",err,error,*999)          
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
   
    EXITS("DistributedMatrix_ValuesGetIntg2")
    RETURN
999 ERRORSEXITS("DistributedMatrix_ValuesGetIntg2",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_ValuesGetIntg2

  !
  !================================================================================================================================
  !

  !>Gets values in a distributed single precision matrix.
  SUBROUTINE DistributedMatrix_ValuesGetSP(distributedMatrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<The row index to get a value from
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<The column index to get a value from
    REAL(SP), INTENT(OUT) :: values(:) !<On return the value of the matrix at the specified row and column
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedMatrix_ValuesGetSP",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)      
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_ValuesGet(cmissMatrix%matrix,rowIndices,columnIndices,values,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot get values for a single precision PETSc distributed matrix.",err,error,*999)          
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_ValuesGetSP")
    RETURN
999 ERRORSEXITS("DistributedMatrix_ValuesGetSP",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_ValuesGetSP

  !
  !================================================================================================================================
  !

  !>Gets one value in a distributed single precision matrix.
  SUBROUTINE DistributedMatrix_ValuesGetSP1(distributedMatrix,rowIndex,columnIndex,value,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(IN) :: rowIndex !<The row index to get a value from
    INTEGER(INTG), INTENT(IN) :: columnIndex !<The column index to get a value from
    REAL(SP), INTENT(OUT) :: value !<On return the value of the matrix at the specified row and column
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedMatrix_ValuesGetSP1",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)      
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_ValuesGet(cmissMatrix%matrix,rowIndex,columnIndex,VALUE,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot get values for a single precision PETSc distributed matrix.",err,error,*999)          
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_ValuesGetSP1")
    RETURN
999 ERRORSEXITS("DistributedMatrix_ValuesGetSP1",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_ValuesGetSP1

  !
  !================================================================================================================================
  !

  !>Gets a matrix of values in a distributed single precision matrix.
  SUBROUTINE DistributedMatrix_ValuesGetSP2(distributedMatrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The ij'th row index to get
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(j). The ij'th column index to get
    REAL(SP), INTENT(OUT) :: values(:,:) !<values(i,j). On return the ij'th value to get
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedMatrix_ValuesGetSP2",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)  
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_ValuesGet(cmissMatrix%matrix,rowIndices,columnIndices,values,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot get values for a single precision PETSc distributed matrix.",err,error,*999)          
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_ValuesGetSP2")
    RETURN
999 ERRORSEXITS("DistributedMatrix_ValuesGetSP2",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_ValuesGetSP2

  !
  !================================================================================================================================
  !

  !>Gets values in a distributed double precision matrix.
  SUBROUTINE DistributedMatrix_ValuesGetDP(distributedMatrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The i'th row index to get
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(i). The i'th column index to get
    REAL(DP), INTENT(OUT) :: values(:) !<values(i). On return the i'th value to get
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: i
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(DistributedMatrixPETScType), POINTER :: petscMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedMatrix_ValuesGetDP",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
   
    SELECT CASE(distributedMatrix%libraryType)     
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_ValuesGet(cmissMatrix%matrix,rowIndices,columnIndices,values,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      IF(SIZE(rowIndices,1)/=SIZE(values,1)) THEN
        localError="The size of the row indices array of "// &
          & TRIM(NumberToVString(SIZE(rowIndices,1),"*",err,error))// &
          & " does not conform to the size of the values array of "// &
          & TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      IF(SIZE(columnIndices,1)/=SIZE(values,1)) THEN
        localError="The size of the column indices array of "// &
          & TRIM(NumberToVString(SIZE(columnIndices,1),"*",err,error))// &
          & " does not conform to the size of the values array of "// &
          & TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      NULLIFY(petscMatrix)
      CALL DistributedMatrix_PETScMatrixGet(distributedMatrix,petscMatrix,err,error,*999)
      IF(petscMatrix%useOverrideMatrix) THEN
        DO i=1,SIZE(rowIndices,1)
          CALL Petsc_MatGetValues(petscMatrix%overrideMatrix,1,petscMatrix%globalRowNumbers(rowIndices(i:i)), &
            & 1,columnIndices(i:i)-1,values(i:i),err,error,*999) !PETSc uses 0 based indices
        ENDDO !i
      ELSE
        DO i=1,SIZE(rowIndices,1)
          CALL Petsc_MatGetValues(petscMatrix%matrix,1,petscMatrix%globalRowNumbers(rowIndices(i:i)), &
            & 1,columnIndices(i:i)-1,values(i:i),err,error,*999) !PETSc uses 0 based indices
        ENDDO !i
      ENDIF
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_ValuesGetDP")
    RETURN
999 ERRORSEXITS("DistributedMatrix_ValuesGetDP",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_ValuesGetDP

  !
  !================================================================================================================================
  !

  !>Gets one value in a distributed double precision matrix.
  SUBROUTINE DistributedMatrix_ValuesGetDP1(distributedMatrix,rowIndex,columnIndex,value,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(IN) :: rowIndex !<The row index to get a value from
    INTEGER(INTG), INTENT(IN) :: columnIndex !<The column index to get a value from
    REAL(DP), INTENT(OUT) :: value !<On return the value of the matrix at the specified row and column
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: columnIndices(1)
    REAL(DP) :: values(1)
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(DistributedMatrixPETScType), POINTER :: petscMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedMatrix_ValuesGetDP1",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_ValuesGet(cmissMatrix%matrix,rowIndex,columnIndex,VALUE,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      NULLIFY(petscMatrix)
      CALL DistributedMatrix_PETScMatrixGet(distributedMatrix,petscMatrix,err,error,*999)
      columnIndices(1)=columnIndex-1
      IF(petscMatrix%useOverrideMatrix) THEN              
        CALL Petsc_MatGetValues(petscMatrix%overrideMatrix,1,petscMatrix%globalRowNumbers(rowIndex), &
          & 1,columnIndices,values,err,error,*999) !PETSc uses 0 based indices
      ELSE
        CALL Petsc_MatGetValues(petscMatrix%matrix,1,petscMatrix%globalRowNumbers(rowIndex), &
          & 1,columnIndices,values,err,error,*999) !PETSc uses 0 based indices
      ENDIF
      value=values(1)
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_ValuesGetDP1")
    RETURN
999 ERRORSEXITS("DistributedMatrix_ValuesGetDP1",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_ValuesGetDP1

  !
  !
  !================================================================================================================================
  !

  !>Gets a matrix of values in a distributed double precision matrix.
  SUBROUTINE DistributedMatrix_ValuesGetDP2(distributedMatrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The ij'th row index to get
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(j). The ij'th column index to get
    REAL(DP), INTENT(OUT) :: values(:,:) !<values(i,j). On return the ij'th value to get
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: globalRowIndices(SIZE(rowIndices,1)),i
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(DistributedMatrixPETScType), POINTER :: petscMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedMatrix_ValuesGetDP2",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_ValuesGet(cmissMatrix%matrix,rowIndices,columnIndices,values,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      IF(SIZE(rowIndices,1)/=SIZE(values,1)) THEN
        localError="The size of the row indices array of "// &
          & TRIM(NumberToVString(SIZE(rowIndices,1),"*",err,error))// &
          & " does not conform to the number of rows in the values array of "// &
          & TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      IF(SIZE(columnIndices,1)/=SIZE(values,2)) THEN
        localError="The size of the column indices array of "// &
          & TRIM(NumberToVString(SIZE(columnIndices,1),"*",err,error))// &
          & " does not conform to the number of columns in the values array of "// &
          & TRIM(NumberToVString(SIZE(values,2),"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      NULLIFY(petscMatrix)
      CALL DistributedMatrix_PETScMatrixGet(distributedMatrix,petscMatrix,err,error,*999)        
      DO i=1,SIZE(rowIndices,1)
        globalRowIndices(i)=petscMatrix%globalRowNumbers(rowIndices(i))
      ENDDO !i
      IF(petscMatrix%useOverrideMatrix) THEN
        CALL Petsc_MatGetValues(petscMatrix%overrideMatrix,SIZE(rowIndices,1),globalRowIndices, &
          & SIZE(columnIndices,1),columnIndices-1,values,err,error,*999) !PETSc uses 0 based row indices
      ELSE
        CALL Petsc_MatGetValues(petscMatrix%matrix,SIZE(rowIndices,1),globalRowIndices, &
          & SIZE(columnIndices,1),columnIndices-1,values,err,error,*999) !PETSc uses 0 based row indices
      ENDIF
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_ValuesGetDP2")
    RETURN
999 ERRORSEXITS("DistributedMatrix_ValuesGetDP2",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_ValuesGetDP2

  !
  !================================================================================================================================
  !

  !>Gets values in a distributed logical matrix.
  SUBROUTINE DistributedMatrix_ValuesGetL(distributedMatrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The i'th row index to get
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(i). The i'th column index to get
    LOGICAL, INTENT(OUT) :: values(:) !<values(i). On return the i'th value to get
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedMatrix_ValuesGetL",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_ValuesGet(cmissMatrix%matrix,rowIndices,columnIndices,values,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot get values for a logical PETSc distributed matrix.",err,error,*999)          
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_ValuesGetL")
    RETURN
999 ERRORSEXITS("DistributedMatrix_ValuesGetL",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_ValuesGetL

  !
  !================================================================================================================================
  !

  !>Gets one value in a distributed logical matrix.
  SUBROUTINE DistributedMatrix_ValuesGetL1(distributedMatrix,rowIndex,columnIndex,value,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(IN) :: rowIndex !<The row index to get a value from
    INTEGER(INTG), INTENT(IN) :: columnIndex !<The column index to get a value from
    LOGICAL, INTENT(OUT) :: value !<On return the value of the matrix at the specified row and column
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(VARYING_STRING) :: localError
   
    ENTERS("DistributedMatrix_ValuesGetL1",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_ValuesGet(cmissMatrix%matrix,rowIndex,columnIndex,VALUE,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot get values for a logical PETSc distributed matrix.",err,error,*999)          
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_ValuesGetL1")
    RETURN
999 ERRORSEXITS("DistributedMatrix_ValuesGetL1",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_ValuesGetL1

  !
  !================================================================================================================================
  !

  !>Gets a matrix of values in a distributed logical matrix.
  SUBROUTINE DistributedMatrix_ValuesGetL2(distributedMatrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The ij'th row index to get
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(j). The ij'th column index to get
    LOGICAL, INTENT(OUT) :: values(:,:) !<values(i,j). On return the ij'th value to get
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedMatrix_ValuesGetL2",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_ValuesGet(cmissMatrix%matrix,rowIndices,columnIndices,values,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot get values for a logical PETSc distributed matrix.",err,error,*999)                    
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_ValuesGetL2")
    RETURN
999 ERRORSEXITS("DistributedMatrix_ValuesGetL2",err,error)
    RETURN 1
  END SUBROUTINE DistributedMatrix_ValuesGetL2

  !
  !================================================================================================================================
  !

  !>Sets values in a distributed integer matrix.
  SUBROUTINE DistributedMatrix_ValuesSetIntg(distributedMatrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The i'th row index to set
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(i). The i'th column index to set
    INTEGER(INTG), INTENT(IN) :: values(:) !<values(i). The i'th value to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedMatrix_ValuesSetIntg",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_ValuesSet(cmissMatrix%matrix,rowIndices,columnIndices,values,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot get values for an integer PETSc distributed matrix.",err,error,*999)          
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_ValuesSetIntg")
    RETURN
999 ERRORSEXITS("DistributedMatrix_ValuesSetIntg",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_ValuesSetIntg

  !
  !================================================================================================================================
  !

  !>Sets one value in a distributed integer matrix.
  SUBROUTINE DistributedMatrix_ValuesSetIntg1(distributedMatrix,rowIndex,columnIndex,value,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(IN) :: rowIndex !<The row index to set a value to
    INTEGER(INTG), INTENT(IN) :: columnIndex !<The column index to set a value to
    INTEGER(INTG), INTENT(IN) :: value !<The value of the matrix to be set at the specified row and column
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedMatrix_ValuesSetIntg1",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_ValuesSet(cmissMatrix%matrix,rowIndex,columnIndex,value,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot get values for an integer PETSc distributed matrix.",err,error,*999)          
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_ValuesSetIntg1")
    RETURN
999 ERRORSEXITS("DistributedMatrix_ValuesSetIntg1",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_ValuesSetIntg1

  !
  !================================================================================================================================
  !

  !>Sets a matrix of values in a distributed integer matrix.
  SUBROUTINE DistributedMatrix_ValuesSetIntg2(distributedMatrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The ij'th row index to set
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(j). The ij'th column index to set
    INTEGER(INTG), INTENT(IN) :: values(:,:) !<values(i,j). The ij'th value to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedMatrix_ValuesSetIntg2",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_ValuesSet(cmissMatrix%matrix,rowIndices,columnIndices,values,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot get values for an integer PETSc distributed matrix.",err,error,*999)          
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_ValuesSetIntg2")
    RETURN
999 ERRORSEXITS("DistributedMatrix_ValuesSetIntg2",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_ValuesSetIntg2

  !
  !================================================================================================================================
  !

  !>Sets values in a distributed single precision matrix.
  SUBROUTINE DistributedMatrix_ValuesSetSP(distributedMatrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The i'th row index to set
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(i). The i'th column index to set
    REAL(SP), INTENT(IN) :: values(:) !<values(i). The i'th value to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedMatrix_ValuesSetSP",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_ValuesSet(cmissMatrix%matrix,rowIndices,columnIndices,values,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot get values for a single precision PETSc distributed matrix.",err,error,*999)          
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_ValuesSetSP")
    RETURN
999 ERRORSEXITS("DistributedMatrix_ValuesSetSP",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_ValuesSetSP

  !
  !================================================================================================================================
  !

  !>Sets one value in a distributed single precision matrix.
  SUBROUTINE DistributedMatrix_ValuesSetSP1(distributedMatrix,rowIndex,columnIndex,value,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(IN) :: rowIndex !<The row index to set a value to
    INTEGER(INTG), INTENT(IN) :: columnIndex !<The column index to set a value to
    REAL(SP), INTENT(IN) :: value !<The value of the matrix to be set at the specified row and column
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedMatrix_ValuesSetSP1",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_ValuesSet(cmissMatrix%matrix,rowIndex,columnIndex,value,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot get values for a single precision PETSc distributed matrix.",err,error,*999)          
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_ValuesSetSP1")
    RETURN
999 ERRORSEXITS("DistributedMatrix_ValuesSetSP1",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_ValuesSetSP1

  !
  !================================================================================================================================
  !

  !>Sets a matrix of values in a distributed single precision matrix.
  SUBROUTINE DistributedMatrix_ValuesSetSP2(distributedMatrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The ij'th row index to set
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(j). The ij'th column index to set
    REAL(SP), INTENT(IN) :: values(:,:) !<values(i,j). The ij'th value to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedMatrix_ValuesSetSP2",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
    
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_ValuesSet(cmissMatrix%matrix,rowIndices,columnIndices,values,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot get values for a single precision PETSc distributed matrix.",err,error,*999)          
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_ValuesSetSP2")
    RETURN
999 ERRORSEXITS("DistributedMatrix_ValuesSetSP2",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_ValuesSetSP2

  !
  !================================================================================================================================
  !

  !>Sets values in a distributed double precision matrix.
  SUBROUTINE DistributedMatrix_ValuesSetDP(distributedMatrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The i'th row index to set
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(i). The i'th column index to set
    REAL(DP), INTENT(IN) :: values(:) !<values(i). The i'th value to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: i
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(DistributedMatrixPETScType), POINTER :: petscMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedMatrix_ValuesSetDP",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
   
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_ValuesSet(cmissMatrix%matrix,rowIndices,columnIndices,values,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      IF(SIZE(rowIndices,1)/=SIZE(values,1)) THEN
        localError="The size of the row indices array of "// &
          & TRIM(NumberToVString(SIZE(rowIndices,1),"*",err,error))// &
          & " does not conform to the size of the values array of "// &
          & TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      IF(SIZE(columnIndices,1)/=SIZE(values,1)) THEN
        localError="The size of the column indices array of "// &
          & TRIM(NumberToVString(SIZE(columnIndices,1),"*",err,error))// &
          & " does not conform to the size of the values array of "// &
          & TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      NULLIFY(petscMatrix)
      CALL DistributedMatrix_PETScMatrixGet(distributedMatrix,petscMatrix,err,error,*999)
      IF(petscMatrix%useOverrideMatrix) THEN
        DO i=1,SIZE(rowIndices,1)
          CALL Petsc_MatSetValues(petscMatrix%overrideMatrix,1,petscMatrix%globalRowNumbers(rowIndices(i:i)), &
            & 1,columnIndices(i:i)-1,values(i:i),PETSC_INSERT_VALUES,err,error,*999) !0 based indices
        ENDDO !i
      ELSE
        DO i=1,SIZE(rowIndices,1)
          CALL Petsc_MatSetValues(petscMatrix%matrix,1,petscMatrix%globalRowNumbers(rowIndices(i:i)), &
            & 1,columnIndices(i:i)-1,values(i:i),PETSC_INSERT_VALUES,err,error,*999) !0 based indices
        ENDDO !i
      ENDIF
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_ValuesSetDP")
    RETURN
999 ERRORSEXITS("DistributedMatrix_ValuesSetDP",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_ValuesSetDP

  !
  !================================================================================================================================
  !

  !>Sets one value in a distributed double precision matrix.
  SUBROUTINE DistributedMatrix_ValuesSetDP1(distributedMatrix,rowIndex,columnIndex,value,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(IN) :: rowIndex !<The row index to set a value to
    INTEGER(INTG), INTENT(IN) :: columnIndex !<The column index to set a value to
    REAL(DP), INTENT(IN) :: value !<The value of the matrix to be set at the specified row and column
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(DistributedMatrixPETScType), POINTER :: petscMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedMatrix_ValuesSetDP1",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
   
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_ValuesSet(cmissMatrix%matrix,rowIndex,columnIndex,value,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      NULLIFY(petscMatrix)
      CALL DistributedMatrix_PETScMatrixGet(distributedMatrix,petscMatrix,err,error,*999)
      IF(petscMatrix%useOverrideMatrix) THEN
        CALL Petsc_MatSetValues(petscMatrix%overrideMatrix,1,petscMatrix%globalRowNumbers( &
          & rowIndex),1,[columnIndex-1],[value],PETSC_INSERT_VALUES,err,error,*999) !PETSc uses 0 based indices
      ELSE
        CALL Petsc_MatSetValues(petscMatrix%matrix,1,petscMatrix%globalRowNumbers(rowIndex), &
          & 1,[columnIndex-1],[value],PETSC_INSERT_VALUES,err,error,*999) !PETSc uses 0 based indices
      ENDIF
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_ValuesSetDP1")
    RETURN
999 ERRORSEXITS("DistributedMatrix_ValuesSetDP1",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_ValuesSetDP1

  !
  !================================================================================================================================
  !

  !>Sets a matrix of values in a distributed double precision matrix.
  SUBROUTINE DistributedMatrix_ValuesSetDP2(distributedMatrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The ij'th row index to set
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(j). The ij'th column index to set
    REAL(DP), INTENT(IN) :: values(:,:) !<values(i,j). The ij'th value to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: globalRowIndices(SIZE(rowIndices,1)),i
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(DistributedMatrixPETScType), POINTER :: petscMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedMatrix_ValuesSetDP2",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
   
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_ValuesSet(cmissMatrix%matrix,rowIndices,columnIndices,values,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      IF(SIZE(rowIndices,1)/=SIZE(values,1)) THEN
        localError="The size of the row indices array of "// &
          & TRIM(NumberToVString(SIZE(rowIndices,1),"*",err,error))// &
          & " does not conform to the number of rows in the values array of "// &
          & TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      IF(SIZE(columnIndices,1)/=SIZE(values,2)) THEN
        localError="The size of the column indices array of "// &
          & TRIM(NumberToVString(SIZE(columnIndices,1),"*",err,error))// &
          & " does not conform to the number of columns in the values array of "// &
          & TRIM(NumberToVString(SIZE(values,2),"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      NULLIFY(petscMatrix)
      CALL DistributedMatrix_PETScMatrixGet(distributedMatrix,petscMatrix,err,error,*999)
      DO i=1,SIZE(rowIndices,1)
        globalRowIndices(i)=petscMatrix%globalRowNumbers(rowIndices(i))
      ENDDO !i
      IF(petscMatrix%useOverrideMatrix) THEN
        CALL Petsc_MatSetValues(petscMatrix%overrideMatrix,SIZE(rowIndices,1),globalRowIndices, &
          & SIZE(columnIndices,1),columnIndices-1,values,PETSC_INSERT_VALUES,err,error,*999) !PETSc uses 0 based indices
      ELSE
        CALL Petsc_MatSetValues(petscMatrix%matrix,SIZE(rowIndices,1),globalRowIndices, &
          & SIZE(columnIndices,1),columnIndices-1,values,PETSC_INSERT_VALUES,err,error,*999) !PETSc uses 0 based indices
      ENDIF
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_ValuesSetDP2")
    RETURN
999 ERRORSEXITS("DistributedMatrix_ValuesSetDP2",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_ValuesSetDP2

  !
  !================================================================================================================================
  !

  !>Sets values in a distributed logical matrix.
  SUBROUTINE DistributedMatrix_ValuesSetL(distributedMatrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The i'th row index to set
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(i). The i'th column index to set
    LOGICAL, INTENT(IN) :: values(:) !<values(i). The i'th value to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedMatrix_ValuesSetL",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
   
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_ValuesSet(cmissMatrix%matrix,rowIndices,columnIndices,values,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot set values for a logical PETSc distributed matrix.",err,error,*999)          
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_ValuesSetL")
    RETURN
999 ERRORSEXITS("DistributedMatrix_ValuesSetL",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_ValuesSetL

  !
  !================================================================================================================================
  !

  !>Sets one value in a distributed logical matrix.
  SUBROUTINE DistributedMatrix_ValuesSetL1(distributedMatrix,rowIndex,columnIndex,value,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(IN) :: rowIndex !<The row index to set a value to
    INTEGER(INTG), INTENT(IN) :: columnIndex !<The column index to set a value to
    LOGICAL, INTENT(IN) :: value !<The value of the matrix to be set at the specified row and column
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedMatrix_ValuesSetL1",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
   
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_ValuesSet(cmissMatrix%matrix,rowIndex,columnIndex,VALUE,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot get values for a logical PETSc distributed matrix.",err,error,*999)          
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_ValuesSetL1")
    RETURN
999 ERRORSEXITS("DistributedMatrix_ValuesSetL1",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_ValuesSetL1

  !
  !================================================================================================================================
  !

  !>Sets a matrix of values in a distributed logical matrix.
  SUBROUTINE DistributedMatrix_ValuesSetL2(distributedMatrix,rowIndices,columnIndices,values,err,error,*)

    !Argument variables
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG), INTENT(IN) :: rowIndices(:) !<rowIndices(i). The ij'th row index to set
    INTEGER(INTG), INTENT(IN) :: columnIndices(:) !<columnIndices(j). The ij'th column index to set
    LOGICAL, INTENT(IN) :: values(:,:) !<values(i,j). The ij'th value to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedMatrix_ValuesSetL2",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("The distributed matrix has not been finished.",err,error,*999)
   
    SELECT CASE(distributedMatrix%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissMatrix)
      CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
      CALL Matrix_ValuesSet(cmissMatrix%matrix,rowIndices,columnIndices,values,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot get values for a logical PETSc distributed matrix.",err,error,*999)          
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedMatrix_ValuesSetL2")
    RETURN
999 ERRORSEXITS("DistributedMatrix_ValuesSetL2",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_ValuesSetL2

  !
  !================================================================================================================================
  !

  !>Calculates the matrix vector product of a distrubted matrix times a distributed vector and adds it to the distributed
  !>product vector. NOTE: This will only work for specific CMISS distributed matrices i.e., ones in which the columns of the
  !>matrix are distributed in the same way as the rows of the multiplied vector are distributed, and the rows of the matrix
  !>are distributed in the same way as the rows of the product vector.
  SUBROUTINE DistributedMatrix_MatrixByVectorAdd(rowSelectionType,alpha,distributedMatrix,distributedVector,distributedProduct, &
    & err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: rowSelectionType !<The row selection for the matrix-vector product \see DISTRIBUTED_MATRIX_VECTOR_GhostingTypes,DISTRIBUTED_MATRIX_VECTOR
    REAL(DP), INTENT(IN) :: alpha !<The multiplicative factor for the distributed matrix
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    TYPE(DistributedVectorType), POINTER :: distributedProduct !<On exit, the value of the matrix vector product added to the distributed product vector.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: columnIdx,localColumn,globalColumn,numberOfColumns,numberOfRows,row,rowIdx
    REAL(DP) :: sum
    TYPE(DistributedMatrixCMISSType), POINTER :: cmissMatrix
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector,cmissProduct
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: rowMapping,columnMapping,productMapping,vectorMapping
    TYPE(MATRIX_TYPE), POINTER :: matrix
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedMatrix_MatrixByVectorAdd",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated",err,error,*999)
    IF(.NOT.distributedMatrix%matrixFinished) CALL FlagError("Distributed matrix has not been finished.",err,error,*999)
    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distrubuted vector is not associated.",err,error,*999)
    IF(.NOT.distributedVector%vectorFinished) CALL FlagError("Distributed vector has not been finished.",err,error,*999)
    IF(.NOT.ASSOCIATED(distributedProduct)) CALL FlagError("The distributed product vector is not associated.",err,error,*999)
    IF(.NOT.distributedProduct%vectorFinished) &
      & CALL FlagError("The distributed product vector has not been finished.",err,error,*999)
    IF(distributedMatrix%libraryType/=distributedVector%libraryType) THEN
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))// &
        & " does not match the distributed matrix library type of "// &
        &  TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(distributedMatrix%libraryType/=distributedProduct%libraryType) THEN
      localError="The distributed product vector library type of "// &
        & TRIM(NumberToVString(distributedProduct%libraryType,"*",err,error))// &
        & " does not match the distributed matrix library type of "// &
        &  TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    NULLIFY(columnMapping)
    CALL DistributedMatrix_ColumnMappingGet(distributedMatrix,columnMapping,err,error,*999)
    NULLIFY(rowMapping)
    CALL DistributedMatrix_RowMappingGet(distributedMatrix,rowMapping,err,error,*999)
    NULLIFY(vectorMapping)
    CALL DistributedVector_RowMappingGet(distributedVector,vectorMapping,err,error,*999)
    NULLIFY(productMapping)
    CALL DistributedVector_RowMappingGet(distributedProduct,productMapping,err,error,*999)
    IF(.NOT.ASSOCIATED(columnMapping,distributedVector%domainMapping)) &
      & CALL FlagError("The distributed matrix and the distributed vector have different domain mappings.", &
      & err,error,*999)
    IF(.NOT.ASSOCIATED(rowMapping,distributedProduct%domainMapping)) &
      & CALL FlagError("The distributed matrix and the distributed product vector have different domain mappings.", &
      & err,error,*999)
      
    IF(ABS(alpha)>ZERO_TOLERANCE) THEN
      SELECT CASE(distributedMatrix%libraryType)
      CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
        NULLIFY(cmissMatrix)
        CALL DistributedMatrix_CMISSMatrixGet(distributedMatrix,cmissMatrix,err,error,*999)
        NULLIFY(cmissVector)
        CALL DistributedVector_CMISSVectorGet(distributedVector,cmissVector,err,error,*999)
        NULLIFY(cmissProduct)
        CALL DistributedVector_CMISSVectorGet(distributedProduct,cmissProduct,err,error,*999)
        matrix=>cmissMatrix%matrix
        IF(.NOT.ASSOCIATED(matrix)) CALL FlagError("CMISS matrix matrix is not associated.",err,error,*999)
        IF(matrix%data_Type/=distributedVector%dataType) THEN
          localError="The distributed vector data type of "// &
            & TRIM(NumberToVString(distributedVector%dataType,"*",err,error))// &
            & " does not match the distributed matrix data type of "// &
            & TRIM(NumberToVString(matrix%data_Type,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        IF(matrix%data_Type/=distributedProduct%dataType) THEN
          localError="The distributed product vector data type of "// &
            & TRIM(NumberToVString(distributedProduct%dataType,"*",err,error))// &
            & " does not match the distributed matrix data type of "// &
            & TRIM(NumberToVString(matrix%data_Type,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
          
        SELECT CASE(rowSelectionType)
        CASE(DISTRIBUTED_MATRIX_VECTOR_INCLUDE_GHOSTS_TYPE)
          numberOfRows=rowMapping%TOTAL_NUMBER_OF_LOCAL
        CASE(DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE)
          numberOfRows=rowMapping%NUMBER_OF_LOCAL
        CASE DEFAULT
          localError="The row selection type of "// &
            & TRIM(NumberToVString(rowSelectionType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        numberOfColumns=columnMapping%NUMBER_OF_GLOBAL
        SELECT CASE(matrix%data_Type)
        CASE(DISTRIBUTED_MATRIX_VECTOR_INTG_TYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(DISTRIBUTED_MATRIX_VECTOR_SP_TYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(DISTRIBUTED_MATRIX_VECTOR_DP_TYPE)
          SELECT CASE(matrix%storage_Type)
          CASE(MATRIX_BLOCK_STORAGE_TYPE)
            DO row=1,numberOfRows
              sum=0.0_DP
              DO localColumn=1,columnMapping%TOTAL_NUMBER_OF_LOCAL
                globalColumn=columnMapping%LOCAL_TO_GLOBAL_MAP(localColumn)
                sum=sum+matrix%data_DP(row+(globalColumn-1)*matrix%m)*cmissVector%dataDP(localColumn)
              ENDDO !localColumn
              cmissProduct%dataDP(row)=cmissProduct%dataDP(row)+alpha*sum
            ENDDO !row
          CASE(MATRIX_DIAGONAL_STORAGE_TYPE)
            DO row=1,numberOfRows
              sum=matrix%data_DP(row)*cmissVector%dataDP(row)
              cmissProduct%dataDP(row)=cmissProduct%dataDP(row)+alpha*sum
            ENDDO !row
          CASE(MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
            DO row=1,numberOfRows
              sum=0.0_DP
              DO localColumn=1,columnMapping%TOTAL_NUMBER_OF_LOCAL
                globalColumn=columnMapping%LOCAL_TO_GLOBAL_MAP(localColumn)
                sum=sum+matrix%data_DP(row+(globalColumn-1)*matrix%MAX_M)*cmissVector%dataDP(localColumn)
              ENDDO !localColumn
              cmissProduct%dataDP(row)=cmissProduct%dataDP(row)+alpha*sum
            ENDDO !row
          CASE(MATRIX_ROW_MAJOR_STORAGE_TYPE)
            DO row=1,numberOfRows
              sum=0.0_DP
              DO localColumn=1,columnMapping%TOTAL_NUMBER_OF_LOCAL
                globalColumn=columnMapping%LOCAL_TO_GLOBAL_MAP(localColumn)
                sum=sum+matrix%data_DP((row-1)*matrix%MAX_N+globalColumn)*cmissVector%dataDP(localColumn)
              ENDDO !localColumn
              cmissProduct%dataDP(row)=cmissProduct%dataDP(row)+alpha*sum
            ENDDO !row
          CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
            DO row=1,numberOfRows
              sum=0.0_DP
              DO columnIdx=matrix%row_Indices(row),matrix%row_Indices(row+1)-1
                globalColumn=matrix%column_Indices(columnIdx)
                !This ranks global to local mappings are stored in the first position
                localColumn=columnMapping%GLOBAL_TO_LOCAL_MAP(globalColumn)%LOCAL_NUMBER(1)
                sum=sum+matrix%data_DP(columnIdx)*cmissVector%dataDP(localColumn)
              ENDDO !localColumn
              cmissProduct%dataDP(row)=cmissProduct%dataDP(row)+alpha*sum
            ENDDO !row
          CASE(MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
            DO columnIdx=1,numberOfColumns
              DO rowIdx=matrix%column_Indices(columnIdx),matrix%column_Indices(columnIdx+1)-1
                row=matrix%row_Indices(rowIdx)
                localColumn=columnMapping%GLOBAL_TO_LOCAL_MAP(columnIdx)%LOCAL_NUMBER(1)
                sum=matrix%data_DP(rowIdx)*cmissVector%dataDP(localColumn)
                cmissProduct%dataDP(row)=cmissProduct%dataDP(row)+alpha*sum
              ENDDO !local_row
            ENDDO !columnIdx
          CASE(MATRIX_ROW_COLUMN_STORAGE_TYPE)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The matrix storage type of "// &
              & TRIM(NumberToVString(matrix%storage_Type,"*",err,error))//" is invalid."            
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(DISTRIBUTED_MATRIX_VECTOR_L_TYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(DISTRIBUTED_MATRIX_VECTOR_SPC_TYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(DISTRIBUTED_MATRIX_VECTOR_DPC_TYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE DEFAULT
          localError="The distributed matrix vector data type of "// &
            & TRIM(NumberToVString(matrix%data_Type,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
     CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
       CALL FlagError("Not implemented.",err,error,*999)
     CASE DEFAULT
       localError="The distributed matrix library type of "// &
         & TRIM(NumberToVString(distributedMatrix%libraryType,"*",err,error))//" is invalid"
       CALL FlagError(localError,err,error,*999)
     END SELECT
   ENDIF
    
    EXITS("DistributedMatrix_MatrixByVectorAdd")
    RETURN
999 ERRORSEXITS("DistributedMatrix_MatrixByVectorAdd",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_MatrixByVectorAdd
  
  !
  !================================================================================================================================
  !
  
  !>Sets all values in an integer distributed vector to the specified value.
  SUBROUTINE DistributedVector_AllValuesSetIntg(distributedVector,value,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    INTEGER(INTG), INTENT(IN) :: value !<The value to be set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedVector_AllValuesSetIntg",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(.NOT.distributedVector%vectorFinished) CALL FlagError("The distributed vector has not been finished.",err,error,*999)
    IF(distributedVector%dataType/=MATRIX_VECTOR_INTG_TYPE) THEN
      localError="The data type of "//TRIM(NumberToVString(distributedVector%dataType,"*",err,error))// &
        & " does not correspond to the integer data type of the given value."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(distributedVector%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissVector)
      CALL DistributedVector_CMISSVectorGet(distributedVector,cmissVector,err,error,*999)
      cmissVector%dataIntg=VALUE
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot get values for an integer PETSc distributed vector.",err,error,*999)          
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("DistributedVector_AllValuesSetIntg")
    RETURN
999 ERRORSEXITS("DistributedVector_AllValuesSetIntg",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_AllValuesSetIntg

  !
  !================================================================================================================================
  !

  !>Sets all values in a single precision distributed vector to the specified value.
  SUBROUTINE DistributedVector_AllValuesSetSP(distributedVector,value,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    REAL(SP), INTENT(IN) :: value !<The value to be set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedVector_AllValuesSetSP",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(.NOT.distributedVector%vectorFinished) CALL FlagError("The distributed vector has not been finished.",err,error,*999)
    IF(distributedVector%dataType/=MATRIX_VECTOR_SP_TYPE) THEN
      localError="The data type of "//TRIM(NumberToVString(distributedVector%dataType,"*",err,error))// &
        & " does not correspond to the single precision data type of the given value."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(distributedVector%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissVector)
      CALL DistributedVector_CMISSVectorGet(distributedVector,cmissVector,err,error,*999)
      cmissVector%dataSP=value
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot get values for a single precision PETSc distributed vector.",err,error,*999)          
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("DistributedVector_AllValuesSetSP")
    RETURN
999 ERRORSEXITS("DistributedVector_AllValuesSetSP",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_AllValuesSetSP

  !
  !================================================================================================================================
  !

  !>Sets all values in a double precision distributed vector to the specified value.
  SUBROUTINE DistributedVector_AllValuesSetDP(distributedVector,value,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    REAL(DP), INTENT(IN) :: value !<The value to be set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector
    TYPE(DistributedVectorPETScType), POINTER :: petscVector
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedVector_AllValuesSetDP",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(.NOT.distributedVector%vectorFinished) CALL FlagError("The distributed vector has not been finished.",err,error,*999)
    IF(distributedVector%dataType/=MATRIX_VECTOR_DP_TYPE) THEN
      localError="The data type of "//TRIM(NumberToVString(distributedVector%dataType,"*",err,error))// &
        & " does not correspond to the double precision data type of the given value."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(distributedVector%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissVector)
      CALL DistributedVector_CMISSVectorGet(distributedVector,cmissVector,err,error,*999)
      cmissVector%dataDP=VALUE
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      NULLIFY(petscVector)
      CALL DistributedVector_PETScVectorGet(distributedVector,petscVector,err,error,*999)
      IF(petscVector%useOverrideVector) THEN
        CALL Petsc_VecSet(petscVector%overrideVector,value,err,error,*999)
      ELSE
        CALL Petsc_VecSet(petscVector%vector,value,err,error,*999)
      ENDIF
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("DistributedVector_AllValuesSetDP")
    RETURN
999 ERRORSEXITS("DistributedVector_AllValuesSetDP",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_AllValuesSetDP

  !
  !================================================================================================================================
  !

  !>Sets all values in a logical distributed vector to the specified value.
  SUBROUTINE DistributedVector_AllValuesSetL(distributedVector,value,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    LOGICAL, INTENT(IN) :: value !<The value to be set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedVector_AllValuesSetL",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(.NOT.distributedVector%vectorFinished) CALL FlagError("The distributed vector has not been finished.",err,error,*999)
    IF(distributedVector%dataType/=MATRIX_VECTOR_L_TYPE) THEN
      localError="The data type of "//TRIM(NumberToVString(distributedVector%dataType,"*",err,error))// &
        & " does not correspond to the logical data type of the given value."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(distributedVector%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissVector)
      CALL DistributedVector_CMISSVectorGet(distributedVector,cmissVector,err,error,*999)
      cmissVector%dataL=value
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot get values for a logical PETSc distributed vector.",err,error,*999)          
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("DistributedVector_AllValuesSetL")
    RETURN
999 ERRORSEXITS("DistributedVector_AllValuesSetL",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_AllValuesSetL

  !
  !================================================================================================================================
  !

  !>Copies alpha times an integer distributed vector to another distributed vector.
  SUBROUTINE DistributedVector_CopyIntg(fromVector,toVector,alpha,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: fromVector !<A pointer to the distributed vector to copy from
    TYPE(DistributedVectorType), POINTER :: toVector !<A pointer to the distributed vector to copy to
    INTEGER(INTG), INTENT(IN) :: alpha !<The multiplicative factor for the copy.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedVectorCMISSType), POINTER :: fromCMISSVector,toCMISSVector
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: fromMapping,toMapping
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedVector_CopyIntg",err,error,*999)

    IF(.NOT.ASSOCIATED(fromVector)) CALL FlagError("From vector is not associated.",err,error,*999)
    IF(.NOT.fromVector%vectorFinished) CALL FlagError("From vector has not been finished.",err,error,*999)
    IF(.NOT.ASSOCIATED(toVector)) CALL FlagError("To vector is not associated.",err,error,*999)
    IF(.NOT.toVector%vectorFinished) CALL FlagError("To vector has not been finished.",err,error,*999)
    IF(fromVector%dataType/=toVector%dataType) THEN
      localError="The from vector data type of "// &
        & TRIM(NumberToVString(fromVector%dataType,"*",err,error))// &
        & " does not match the to vector data type of "// &
        & TRIM(NumberToVString(toVector%dataType,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(fromVector%dataType/=DISTRIBUTED_MATRIX_VECTOR_INTG_TYPE) THEN
      localError="The from vector data type of "//TRIM(NumberToVString(fromVector%dataType,"*",err,error))// &
        & " does not match the integer data type of the supplied alpha value."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(fromVector%libraryType/=toVector%libraryType) CALL FlagError("Not implemented.",err,error,*999) !Vectors are of from different library types
    NULLIFY(fromMapping)
    CALL DistributedVector_RowMappingGet(fromVector,fromMapping,err,error,*999)
    NULLIFY(toMapping)
    CALL DistributedVector_RowMappingGet(toVector,toMapping,err,error,*999)
    IF(.NOT.ASSOCIATED(fromVector%domainMapping,toVector%domainMapping)) &
      & CALL FlagError("The from vector does not have the same domain mapping as the to vector.",err,error,*999)
     
    SELECT CASE(fromVector%libraryType)      
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(fromCMISSVector)
      CALL DistributedVector_CMISSVectorGet(fromVector,fromCMISSVector,err,error,*999)
      NULLIFY(toCMISSVector)
      CALL DistributedVector_CMISSVectorGet(toVector,toCMISSVector,err,error,*999)
      toCMISSVector%dataIntg(1:toCMISSVector%n)=alpha*fromCMISSVector%dataIntg(1:fromCMISSVector%n)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot copy a vector fro an integer PETSc distributed vector.",err,error,*999)
    CASE DEFAULT
      localError="The from vector library type of "// &
        & TRIM(NumberToVString(fromVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
     
    EXITS("DistributedVector_CopyIntg")
    RETURN
999 ERRORSEXITS("DistributedVector_CopyIntg",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_CopyIntg

  !
  !================================================================================================================================
  !

  !>Copies alpha times a single precision distributed vector to another distributed vector.
  SUBROUTINE DistributedVector_CopySP(fromVector,toVector,alpha,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: fromVector !<A pointer to the distributed vector to copy from
    TYPE(DistributedVectorType), POINTER :: toVector !<A pointer to the distributed vector to copy to
    REAL(SP), INTENT(IN) :: alpha !<The multiplicative factor for the copy.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedVectorCMISSType), POINTER :: fromCMISSVector,toCMISSVector
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: fromMapping,toMapping
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedVector_CopySP",err,error,*999)

    IF(.NOT.ASSOCIATED(fromVector)) CALL FlagError("From vector is not associated.",err,error,*999)
    IF(.NOT.fromVector%vectorFinished) CALL FlagError("From vector has not been finished.",err,error,*999)
    IF(.NOT.ASSOCIATED(toVector)) CALL FlagError("To vector is not associated.",err,error,*999)
    IF(.NOT.toVector%vectorFinished) CALL FlagError("To vector has not been finished.",err,error,*999)
    IF(fromVector%dataType/=toVector%dataType) THEN
      localError="The from vector data type of "// &
        & TRIM(NumberToVString(fromVector%dataType,"*",err,error))// &
        & " does not match the to vector data type of "// &
        & TRIM(NumberToVString(toVector%dataType,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(fromVector%dataType/=DISTRIBUTED_MATRIX_VECTOR_SP_TYPE) THEN
      localError="The from vector data type of "//TRIM(NumberToVString(fromVector%dataType,"*",err,error))// &
        & " does not match the single precision real data type of the supplied alpha value."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(fromVector%libraryType/=toVector%libraryType) CALL FlagError("Not implemented.",err,error,*999) !Vectors are of from different library types
    NULLIFY(fromMapping)
    CALL DistributedVector_RowMappingGet(fromVector,fromMapping,err,error,*999)
    NULLIFY(toMapping)
    CALL DistributedVector_RowMappingGet(toVector,toMapping,err,error,*999)
    IF(.NOT.ASSOCIATED(fromVector%domainMapping,toVector%domainMapping)) &
      & CALL FlagError("The from vector does not have the same domain mapping as the to vector.",err,error,*999)

    
    SELECT CASE(fromVector%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(fromCMISSVector)
      CALL DistributedVector_CMISSVectorGet(fromVector,fromCMISSVector,err,error,*999)
      NULLIFY(toCMISSVector)
      CALL DistributedVector_CMISSVectorGet(toVector,toCMISSVector,err,error,*999)
      toCMISSVector%dataSP(1:toCMISSVector%n)=alpha*fromCMISSVector%dataSP(1:fromCMISSVector%n)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot copy a vector for a single precision PETSc distributed vector.",err,error,*999)
    CASE DEFAULT
      localError="The from vector library type of "// &
        & TRIM(NumberToVString(fromVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
     
    EXITS("DistributedVector_CopySP")
    RETURN
999 ERRORSEXITS("DistributedVector_CopySP",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_CopySP

  !
  !================================================================================================================================
  !

  !>Copies alpha times a double precision distributed vector to another distributed vector.
  SUBROUTINE DistributedVector_CopyDP(fromVector,toVector,alpha,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: fromVector !<A pointer to the distributed vector to copy from
    TYPE(DistributedVectorType), POINTER :: toVector !<A pointer to the distributed vector to copy to
    REAL(DP), INTENT(IN) :: alpha !<The multiplicative factor for the copy.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedVectorCMISSType), POINTER :: fromCMISSVector,toCMISSVector
    TYPE(DistributedVectorPETScType), POINTER :: fromPETScVector,toPETScVector
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: fromMapping,toMapping
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedVector_CopyDP",err,error,*999)

    IF(.NOT.ASSOCIATED(fromVector)) CALL FlagError("From vector is not associated.",err,error,*999)
    IF(.NOT.fromVector%vectorFinished) CALL FlagError("From vector has not been finished.",err,error,*999)
    IF(.NOT.ASSOCIATED(toVector)) CALL FlagError("To vector is not associated.",err,error,*999)
    IF(.NOT.toVector%vectorFinished) CALL FlagError("To vector has not been finished.",err,error,*999)
    IF(fromVector%dataType/=toVector%dataType) THEN
      localError="The from vector data type of "// &
        & TRIM(NumberToVString(fromVector%dataType,"*",err,error))// &
        & " does not match the to vector data type of "// &
        & TRIM(NumberToVString(toVector%dataType,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(fromVector%dataType/=DISTRIBUTED_MATRIX_VECTOR_DP_TYPE) THEN
      localError="The from vector data type of "//TRIM(NumberToVString(fromVector%dataType,"*",err,error))// &
        & " does not match the double precision real data type of the supplied alpha value."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(fromVector%libraryType/=toVector%libraryType) CALL FlagError("Not implemented.",err,error,*999) !Vectors are of from different library types
    NULLIFY(fromMapping)
    CALL DistributedVector_RowMappingGet(fromVector,fromMapping,err,error,*999)
    NULLIFY(toMapping)
    CALL DistributedVector_RowMappingGet(toVector,toMapping,err,error,*999)
    IF(.NOT.ASSOCIATED(fromVector%domainMapping,toVector%domainMapping)) &
      & CALL FlagError("The from vector does not have the same domain mapping as the to vector.",err,error,*999)
    
    SELECT CASE(fromVector%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(fromCMISSVector)
      CALL DistributedVector_CMISSVectorGet(fromVector,fromCMISSVector,err,error,*999)
      NULLIFY(toCMISSVector)
      CALL DistributedVector_CMISSVectorGet(toVector,toCMISSVector,err,error,*999)
      toCMISSVector%dataDP(1:toCMISSVector%n)=alpha*fromCMISSVector%dataDP(1:fromCMISSVector%n)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      NULLIFY(fromPETScVector)
      CALL DistributedVector_PETScVectorGet(fromVector,fromPETScVector,err,error,*999)
      NULLIFY(toPETScVector)
      CALL DistributedVector_PETScVectorGet(toVector,toPETScVector,err,error,*999)
      IF(fromPETScVector%useOverrideVector) THEN
        IF(toPETScVector%useOverrideVector) THEN
          CALL Petsc_VecCopy(fromPETScVector%overrideVector,toPETScVector%overrideVector,err,error,*999)
          CALL Petsc_VecScale(toPETScVector%overrideVector,alpha,err,error,*999)
        ELSE
          CALL Petsc_VecCopy(fromPETScVector%overrideVector,toPETScVector%vector,err,error,*999)
          CALL Petsc_VecScale(toPETScVector%vector,alpha,err,error,*999)
        ENDIF
      ELSE
        IF(toPETScVector%useOverrideVector) THEN
          CALL Petsc_VecCopy(fromPETScVector%vector,toPETScVector%overrideVector,err,error,*999)
          CALL Petsc_VecScale(toPETScVector%overrideVector,alpha,err,error,*999)
        ELSE
          CALL Petsc_VecCopy(fromPETScVector%vector,toPETScVector%vector,err,error,*999)
          CALL Petsc_VecScale(toPETScVector%vector,alpha,err,error,*999)
        ENDIF
      ENDIF
    CASE DEFAULT
      localError="The from vector library type of "// &
        & TRIM(NumberToVString(fromVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
     
    EXITS("DistributedVector_CopyDP")
    RETURN
999 ERRORSEXITS("DistributedVector_CopyDP",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_CopyDP

  !
  !================================================================================================================================
  !

  !>Copies alpha times a logical distributed vector to another distributed vector.
  SUBROUTINE DistributedVector_CopyL(fromVector,toVector,alpha,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: fromVector !<A pointer to the distributed vector to copy from
    TYPE(DistributedVectorType), POINTER :: toVector !<A pointer to the distributed vector to copy to
    LOGICAL, INTENT(IN) :: alpha !<The multiplicative factor for the copy.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedVectorCMISSType), POINTER :: fromCMISSVector,toCMISSVector
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: fromMapping,toMapping
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedVector_CopyL",err,error,*999)

    IF(.NOT.ASSOCIATED(fromVector)) CALL FlagError("From vector is not associated.",err,error,*999)
    IF(.NOT.fromVector%vectorFinished) CALL FlagError("From vector has not been finished.",err,error,*999)
    IF(.NOT.ASSOCIATED(toVector)) CALL FlagError("To vector is not associated.",err,error,*999)
    IF(.NOT.toVector%vectorFinished) CALL FlagError("To vector has not been finished.",err,error,*999)
    IF(fromVector%dataType/=toVector%dataType) THEN
      localError="The from vector data type of "// &
        & TRIM(NumberToVString(fromVector%dataType,"*",err,error))// &
        & " does not match the to vector data type of "// &
        & TRIM(NumberToVString(toVector%dataType,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(fromVector%dataType/=DISTRIBUTED_MATRIX_VECTOR_L_TYPE) THEN
      localError="The from vector data type of "//TRIM(NumberToVString(fromVector%dataType,"*",err,error))// &
        & " does not match the logical data type of the supplied alpha value."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(fromVector%libraryType/=toVector%libraryType) CALL FlagError("Not implemented.",err,error,*999) !Vectors are of from different library types
    NULLIFY(fromMapping)
    CALL DistributedVector_RowMappingGet(fromVector,fromMapping,err,error,*999)
    NULLIFY(toMapping)
    CALL DistributedVector_RowMappingGet(toVector,toMapping,err,error,*999)
    IF(.NOT.ASSOCIATED(fromVector%domainMapping,toVector%domainMapping)) &
      & CALL FlagError("The from vector does not have the same domain mapping as the to vector.",err,error,*999)
     
    SELECT CASE(fromVector%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(fromCMISSVector)
      CALL DistributedVector_CMISSVectorGet(fromVector,fromCMISSVector,err,error,*999)
      NULLIFY(toCMISSVector)
      CALL DistributedVector_CMISSVectorGet(toVector,toCMISSVector,err,error,*999)
      toCMISSVector%dataL(1:toCMISSVector%n)=alpha.AND.fromCMISSVector%dataL(1:fromCMISSVector%n)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot copy a vector for an integer PETSc distributed vector.",err,error,*999)
    CASE DEFAULT
      localError="The from vector library type of "// &
        & TRIM(NumberToVString(fromVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
     
    EXITS("DistributedVector_CopyL")
    RETURN
999 ERRORSEXITS("DistributedVector_CopyL",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_CopyL

  !
  !================================================================================================================================
  !

  !>Finalise a CMISS distributed vector.
  SUBROUTINE DistributedVector_CMISSFinalise(cmissVector,err,error,*)

    !Argument variables
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector !<A pointer to the CMISS distributed vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: domainIdx
    
    ENTERS("DistributedVector_CMISSFinalise",err,error,*999)

    IF(ASSOCIATED(cmissVector)) THEN
      IF(ALLOCATED(cmissVector%dataIntg)) DEALLOCATE(cmissVector%dataIntg)
      IF(ALLOCATED(cmissVector%dataSP)) DEALLOCATE(cmissVector%dataSP)
      IF(ALLOCATED(cmissVector%dataDP)) DEALLOCATE(cmissVector%dataDP)
      IF(ALLOCATED(cmissVector%dataL)) DEALLOCATE(cmissVector%dataL)
      IF(ALLOCATED(cmissVector%transfers)) THEN
        DO domainIdx=1,SIZE(cmissVector%transfers)
          CALL DistributedVector_CMISSTransferFinalise(cmissVector,domainIdx,err,error,*999)
        ENDDO !domain_idx
        DEALLOCATE(cmissVector%transfers)
      ENDIF
      DEALLOCATE(cmissVector)
    ENDIF
     
    EXITS("DistributedVector_CMISSFinalise")
    RETURN
999 ERRORSEXITS("DistributedVector_CMISSFinalise",err,error)
    RETURN 1
  END SUBROUTINE DistributedVector_CMISSFinalise

  !
  !================================================================================================================================
  !

  !>Intialises a CMISS distributed vector.
  SUBROUTINE DistributedVector_CMISSInitialise(distributedVector,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError,localError
    
    ENTERS("DistributedVector_CMISSInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*998)
    IF(ASSOCIATED(distributedVector%cmiss)) &
      & CALL FlagError("CMISS is already associated for this distributed vector.",err,error,*998)
    IF(.NOT.ASSOCIATED(distributedVector%domainMapping)) &
      & CALL FlagError("Distributed vector domain mapping is not associated.",err,error,*998)
    
    ALLOCATE(distributedVector%cmiss,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocated CMISS distributed vector.",err,error,*999)
    distributedVector%cmiss%distributedVector=>distributedVector
    distributedVector%libraryType=DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE
    !Set the defaults
    distributedVector%cmiss%baseTagNumber=0
    SELECT CASE(distributedVector%ghostingType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_INCLUDE_GHOSTS_TYPE)
      distributedVector%cmiss%n=distributedVector%domainMapping%TOTAL_NUMBER_OF_LOCAL
    CASE(DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE)
      distributedVector%cmiss%n=distributedVector%domainMapping%NUMBER_OF_LOCAL
    CASE DEFAULT
      localError="The distributed vector ghosting type of "// &
        & TRIM(NumberToVString(distributedVector%ghostingType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    distributedVector%cmiss%dataSize=0         
    
    EXITS("DistributedVector_CMISSInitialise")
    RETURN
999 CALL DistributedVector_CMISSFinalise(distributedVector%cmiss,dummyErr,dummyError,*998)
998 ERRORSEXITS("DistributedVector_CMISSInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_CMISSInitialise

  !
  !================================================================================================================================
  !

  !>Finishes the creation of a CMISS distributed vector
  SUBROUTINE DistributedVector_CMISSCreateFinish(cmissVector,err,error,*)

    !Argument variables
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector !<A pointer to the distributed CMISS vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: domainIdx,domainIdx2,domainNumber,dummyErr,myComputationalNodeNumber
    LOGICAL :: found
    TYPE(DistributedVectorType), POINTER :: distributedVector
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: domainMapping
    TYPE(VARYING_STRING) :: dummyError,localError
    
    ENTERS("DistributedVector_CMISSCreateFinish",err,error,*999)

    IF(.NOT.ASSOCIATED(cmissVector)) CALL FlagError("CMISS vector is not associated.",err,error,*999)
    NULLIFY(distributedVector)
    CALL DistributedVectorCMISS_DistributedVectorGet(cmissVector,distributedVector,err,error,*999)
    NULLIFY(domainMapping)
    CALL DistributedVector_RowMappingGet(distributedVector,domainMapping,err,error,*999)
    
    cmissVector%dataSize=cmissVector%n    
    SELECT CASE(distributedVector%dataType)
    CASE(MATRIX_VECTOR_INTG_TYPE)
      ALLOCATE(cmissVector%dataIntg(cmissVector%dataSize),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate CMISS distributed vector integer data.",err,error,*999)
    CASE(MATRIX_VECTOR_SP_TYPE)
      ALLOCATE(cmissVector%dataSP(cmissVector%dataSize),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate CMISS distributed vector single precsion data.",err,error,*999)
    CASE(MATRIX_VECTOR_DP_TYPE)
      ALLOCATE(cmissVector%dataDP(cmissVector%dataSize),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate CMISS distributed vector double precsion data.",err,error,*999)
    CASE(MATRIX_VECTOR_L_TYPE)
      ALLOCATE(cmissVector%dataL(cmissVector%dataSize),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate CMISS distributed vector logical data.",err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_SPC_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_DPC_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The distributed vector data type of "// &
        & TRIM(NumberToVString(distributedVector%dataType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    cmissVector%baseTagNumber=distributedDataId
    IF(domainMapping%NUMBER_OF_DOMAINS==1) THEN
      distributedDataId=distributedDataId+1
    ELSE
      distributedDataId=distributedDataId+domainMapping%ADJACENT_DOMAINS_PTR(domainMapping%NUMBER_OF_DOMAINS)
    END IF
    IF(domainMapping%NUMBER_OF_ADJACENT_DOMAINS>0) THEN
      myComputationalNodeNumber=ComputationalEnvironment_NodeNumberGet(err,error)
      IF(err/=0) GOTO 999
      IF(distributedVector%ghostingType==DISTRIBUTED_MATRIX_VECTOR_INCLUDE_GHOSTS_TYPE) THEN
        ALLOCATE(cmissVector%transfers(domainMapping%NUMBER_OF_ADJACENT_DOMAINS),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate CMISS distributed vector transfer buffers.",err,error,*999)
        DO domainIdx=1,domainMapping%NUMBER_OF_ADJACENT_DOMAINS
          CALL DistributedVector_CMISSTransferInitialise(cmissVector,domainIdx,err,error,*999)
          cmissVector%transfers(domainIdx)%sendBufferSize=domainMapping%ADJACENT_DOMAINS(domainIdx)%NUMBER_OF_SEND_GHOSTS
          cmissVector%transfers(domainIdx)%receiveBufferSize=domainMapping%ADJACENT_DOMAINS(domainIdx)%NUMBER_OF_RECEIVE_GHOSTS
          cmissVector%transfers(domainIdx)%dataType=distributedVector%dataType
          cmissVector%transfers(domainIdx)%sendTagNumber=cmissVector%baseTagNumber+ &
            & domainMapping%ADJACENT_DOMAINS_PTR(myComputationalNodeNumber)+domainIdx-1
          domainNumber=domainMapping%ADJACENT_DOMAINS(domainIdx)%DOMAIN_NUMBER
          found=.FALSE.
          DO domainIdx2=domainMapping%ADJACENT_DOMAINS_PTR(domainNumber),domainMapping%ADJACENT_DOMAINS_PTR(domainNumber+1)-1
            IF(domainMapping%ADJACENT_DOMAINS_LIST(domainIdx2)==myComputationalNodeNumber) THEN
              found=.TRUE.
              EXIT
            ENDIF
          ENDDO !domainIdx2
          IF(.NOT.found) CALL FlagError("Could not find domain to set the receive tag number.",err,error,*999)
            domainIdx2=domainIdx2-domainMapping%ADJACENT_DOMAINS_PTR(domainNumber)+1
            cmissVector%transfers(domainIdx)%receiveTagNumber=cmissVector%baseTagNumber+ &
              & domainMapping%ADJACENT_DOMAINS_PTR(domainNumber)+domainIdx2-1
          SELECT CASE(distributedVector%dataType)
          CASE(DISTRIBUTED_MATRIX_VECTOR_INTG_TYPE)
            ALLOCATE(cmissVector%transfers(domainIdx)%sendBufferIntg(cmissVector%transfers(domainIdx)%sendBufferSize),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate distributed vector send integer transfer buffer.",err,error,*999)
            ALLOCATE(cmissVector%transfers(domainIdx)%receiveBufferIntg(cmissVector%transfers(domainIdx)%receiveBufferSize), &
              & STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate distributed vector receive integer transfer buffer.",err,error,*999)
          CASE(DISTRIBUTED_MATRIX_VECTOR_SP_TYPE)
            ALLOCATE(cmissVector%transfers(domainIdx)%sendBufferSP(cmissVector%transfers(domainIdx)%sendBufferSize),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate distributed vector send single precision transfer buffer.", &
              & err,error,*999)
            ALLOCATE(cmissVector%transfers(domainIdx)%receiveBufferSP(cmissVector%transfers(domainIdx)%receiveBufferSize), &
              & STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate distributed vector receive single precision transfer buffer.", &
              & err,error,*999)
          CASE(DISTRIBUTED_MATRIX_VECTOR_DP_TYPE)
            ALLOCATE(cmissVector%transfers(domainIdx)%sendBufferDP(cmissVector%transfers(domainIdx)%sendBufferSize),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate distributed vector send double precision transfer buffer.", &
              & err,error,*999)
            ALLOCATE(cmissVector%transfers(domainIdx)%receiveBufferDP(cmissVector%transfers(domainIdx)%receiveBufferSize), &
              & STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate distributed vector receive double precision transfer buffer.", &
              & err,error,*999)
          CASE(DISTRIBUTED_MATRIX_VECTOR_L_TYPE)
            ALLOCATE(cmissVector%transfers(domainIdx)%sendBufferL(cmissVector%transfers(domainIdx)%sendBufferSize),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate distributed vector send logical transfer buffer.",err,error,*999)
            ALLOCATE(cmissVector%transfers(domainIdx)%receiveBufferL(cmissVector%transfers(domainIdx)%receiveBufferSize), &
              & STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate distributed vector receive logical transfer buffer.",err,error,*999)
          CASE(DISTRIBUTED_MATRIX_VECTOR_SPC_TYPE)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(DISTRIBUTED_MATRIX_VECTOR_DPC_TYPE)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The distributed vector data type of "// &
              & TRIM(NumberToVString(distributedVector%dataType,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ENDDO !domainIdx
      ENDIF
    ENDIF
    
    EXITS("DistributedVector_CMISSCreateFinish")
    RETURN
999 CALL DistributedVector_CMISSFinalise(cmissVector,dummyErr,dummyError,*998)      
998 ERRORSEXITS("DistributedVector_CMISSCreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_CMISSCreateFinish

  !
  !================================================================================================================================
  !

  !>Finishes the creation a distributed vector
  SUBROUTINE DistributedVector_CreateFinish(distributedVector,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector
    TYPE(DistributedVectorPETScType), POINTER :: petscVector
    TYPE(VARYING_STRING) :: dummyError,localError
    
    ENTERS("DistributedVector_CreateFinish",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(distributedVector%vectorFinished) CALL FlagError("The distributed vector has already been finished.",err,error,*999)
    
    SELECT CASE(distributedVector%libraryType)      
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissVector)
      CALL DistributedVector_CMISSVectorGet(distributedVector,cmissVector,err,error,*999)
      CALL DistributedVector_CMISSCreateFinish(cmissVector,err,error,*999)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      NULLIFY(petscVector)
      CALL DistributedVector_PETScVectorGet(distributedVector,petscVector,err,error,*999)
      CALL DistributedVector_PETScCreateFinish(petscVector,err,error,*999)
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    distributedVector%vectorFinished=.TRUE.
    
    EXITS("DistributedVector_CreateFinish")
    RETURN
999 CALL DistributedVector_Finalise(distributedVector,dummyErr,dummyError,*998)    
998 ERRORSEXITS("DistributedVector_CreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_CreateFinish

  !
  !================================================================================================================================
  !

  !>Starts the creation a distributed vector.
  SUBROUTINE DistributedVector_CreateStart(domainMapping,distributedVector,err,error,*)

    !Argument variables
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: domainMapping !<A pointer to the domain mapping used to distribute this vector
    TYPE(DistributedVectorType), POINTER :: distributedVector !<On return, a pointer to the created distributed vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("DistributedVector_CreateStart",err,error,*998)

    IF(.NOT.ASSOCIATED(domainMapping)) CALL FlagError("Domain mapping is not associated.",err,error,*998)
    IF(ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is already associated.",err,error,*998)
    
    CALL DistributedVector_Initialise(domainMapping,distributedVector,err,error,*999)
    !Set the default values
    
    EXITS("DistributedVector_CreateStart")
    RETURN
999 CALL DistributedVector_Finalise(distributedVector,dummyErr,dummyError,*998)
998 ERRORSEXITS("DistributedVector_CreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_CreateStart

  !
  !================================================================================================================================
  !

  !>Gets the data type of a distributed vector.
  SUBROUTINE DistributedVector_DataTypeGet(distributedVector,dataType,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    INTEGER(INTG), INTENT(OUT) :: dataType !<On return, the data type of the vector. \see DISTRIBUTED_MATRIX_VECTOR_DataTypes,MATRIX_VECTOR
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("DistributedVector_DataTypeGet",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(.NOT.distributedVector%vectorFinished) CALL FlagError("The distributed vector has not been finished.",err,error,*999)
    
    dataType=distributedVector%dataType

    EXITS("DistributedVector_DataTypeGet")
    RETURN
999 ERRORSEXITS("DistributedVector_DataTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_DataTypeGet

  !
  !================================================================================================================================
  !

  !>Sets/changes the data type of a distributed vector.
  SUBROUTINE DistributedVector_DataTypeSet(distributedVector,dataType,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    INTEGER(INTG), INTENT(IN) :: dataType !<The data type to be set \see DISTRIBUTED_MATRIX_VECTOR_DataTypes,DISTRIBUTED_MATRIX_VECTOR
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedVector_DataTypeSet",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(distributedVector%vectorFinished) CALL FlagError("The distributed vector has been finished.",err,error,*999)
    
    SELECT CASE(distributedVector%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      SELECT CASE(dataType)
      CASE(MATRIX_VECTOR_INTG_TYPE)
        distributedVector%dataType=MATRIX_VECTOR_INTG_TYPE
      CASE(MATRIX_VECTOR_SP_TYPE)
        distributedVector%dataType=MATRIX_VECTOR_SP_TYPE
      CASE(MATRIX_VECTOR_DP_TYPE)
        distributedVector%dataType=MATRIX_VECTOR_DP_TYPE
      CASE(MATRIX_VECTOR_L_TYPE)
        distributedVector%dataType=MATRIX_VECTOR_L_TYPE
      CASE(DISTRIBUTED_MATRIX_VECTOR_SPC_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(DISTRIBUTED_MATRIX_VECTOR_DPC_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The distributed data type of "//TRIM(NumberToVString(dataType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      SELECT CASE(dataType)
      CASE(MATRIX_VECTOR_INTG_TYPE)
        CALL FlagError("An integer distributed PETSc vector is not implemented.",err,error,*999)
      CASE(MATRIX_VECTOR_SP_TYPE)
        CALL FlagError("A single precision distributed PETSc vector is not implemented.",err,error,*999)
      CASE(MATRIX_VECTOR_DP_TYPE)
        distributedVector%dataType=MATRIX_VECTOR_DP_TYPE
      CASE(MATRIX_VECTOR_L_TYPE)
        CALL FlagError("A logical distributed PETSc vector is not implemented.",err,error,*999)
      CASE(DISTRIBUTED_MATRIX_VECTOR_SPC_TYPE)
        CALL FlagError("A single precision complex distributed PETSc vector is not implemented.",err,error,*999)
      CASE(DISTRIBUTED_MATRIX_VECTOR_DPC_TYPE)
        CALL FlagError("A double precision complex distributed PETSc vector is not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The distributed data type of "//TRIM(NumberToVString(dataType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedVector_DataTypeSet")
    RETURN
999 ERRORSEXITS("DistributedVector_DataTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_DataTypeSet

  !
  !================================================================================================================================
  !

  !>Destroys a distributed vector.
  SUBROUTINE DistributedVector_Destroy(distributedVector,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DistributedVector_Destroy",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    
    CALL DistributedVector_Finalise(distributedVector,err,error,*999)
    
    EXITS("DistributedVector_Destroy")
    RETURN
999 ERRORSEXITS("DistributedVector_Destroy",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_Destroy

  !
  !================================================================================================================================
  !

  !>Duplicates the structure of a distributed vector and returns a pointer to the new distributed vector in newDistributedVector.
  SUBROUTINE DistributedVector_Duplicate(distributedVector,newDistributedVector,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector to duplicate
    TYPE(DistributedVectorType), POINTER :: newDistributedVector !<On return a pointer to the new duplicated distributed vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("DistributedVector_Duplicate",err,error,*998)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*998)
    IF(ASSOCIATED(newDistributedVector)) CALL FlagError("New distributed vector is already associated.",err,error,*998)
    
    CALL DistributedVector_CreateStart(distributedVector%domainMapping,newDistributedVector,err,error,*999)    
    CALL DistributedVector_LibraryTypeSet(newDistributedVector,distributedVector%libraryType,err,error,*999)
    CALL DistributedVector_DataTypeSet(newDistributedVector,distributedVector%dataType,err,error,*999)
    CALL DistributedVector_CreateFinish(newDistributedVector,err,error,*999)
    
    EXITS("DistributedVector_Duplicate")
    RETURN
999 CALL DistributedVector_Finalise(newDistributedVector,dummyErr,dummyError,*998)
998 ERRORSEXITS("DistributedVector_Duplicate",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_Duplicate

  !
  !================================================================================================================================
  !

  !>Finalises a distributed vector and deallocates all memory.
  SUBROUTINE DistributedVector_Finalise(distributedVector,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DistributedVector_Finalise",err,error,*999)

    IF(ASSOCIATED(distributedVector)) THEN
      CALL DistributedVector_CMISSFinalise(distributedVector%cmiss,err,error,*999)
      CALL DistributedVector_PETScFinalise(distributedVector%petsc,err,error,*999)        
      DEALLOCATE(distributedVector)
    ENDIF
    
    EXITS("DistributedVector_Finalise")
    RETURN
999 ERRORSEXITS("DistributedVector_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises a distributed vector.
  SUBROUTINE DistributedVector_Initialise(domainMapping,distributedVector,err,error,*)

    !Argument variables
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: domainMapping !<A pointer to the domain mapping used to distribute this vector
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("DistributedVector_Initialise",err,error,*998)

    IF(.NOT.ASSOCIATED(domainMapping)) CALL FlagError("Domain mapping is not associated.",err,error,*998)
    IF(ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is already associated.",err,error,*998)
    
    ALLOCATE(distributedVector,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocated the distributed vector.",err,error,*999)
    distributedVector%vectorFinished=.FALSE.
    distributedVector%libraryType=0
    distributedVector%ghostingType=DISTRIBUTED_MATRIX_VECTOR_INCLUDE_GHOSTS_TYPE
    distributedVector%domainMapping=>domainMapping
    distributedVector%dataType=MATRIX_VECTOR_DP_TYPE
    NULLIFY(distributedVector%cmiss)
    NULLIFY(distributedVector%petsc)
    CALL DistributedVector_CMISSInitialise(distributedVector,err,error,*999)
    
    EXITS("DistributedVector_Initialise")
    RETURN
999 CALL DistributedVector_Finalise(distributedVector,dummyErr,dummyError,*998)
998 ERRORSEXITS("DistributedVector_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_Initialise

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the data of an integer distributed vector. Note: the values can be used for read operations but a DistributedVector_ValuesSet call must be used to change any values. The pointer should not be deallocated.
  SUBROUTINE DistributedVector_DataGetIntg(distributedVector,data,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    INTEGER(INTG), POINTER :: data(:) !<On return, a pointer to the data of the distributed vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedVector_DataGetIntg",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(ASSOCIATED(DATA)) CALL FlagError("Data is already associated.",err,error,*999)     
    IF(.NOT.distributedVector%vectorFinished) CALL FlagError("The distributed vector has not been finished.",err,error,*999)
    IF(distributedVector%dataType/=MATRIX_VECTOR_INTG_TYPE) THEN
      localError="The distributed data type of "// &
        & TRIM(NumberToVString(distributedVector%dataType,"*",err,error))// &
        & " does not correspond to the integer data type of the requested values."
      CALL FlagError(localError,err,error,*999)
    ENDIF
     
    SELECT CASE(distributedVector%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissVector)
      CALL DistributedVector_CMISSVectorGet(distributedVector,cmissVector,err,error,*999)
      data=>cmissVector%dataIntg
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot get data for an integer PETSc distributed vector.",err,error,*999)          
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedVector_DataGetIntg")
    RETURN
999 ERRORSEXITS("DistributedVector_DataGetIntg",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_DataGetIntg

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the data of a single precision distributed vector. Note: the values can be used for read operations but a DistributedVector_ValuesSet call must be used to change any values. The pointer should not be deallocated.
  SUBROUTINE DistributedVector_DataGetSP(distributedVector,data,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    REAL(SP), POINTER :: data(:) !<On return, a pointer to the data of the distributed vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedVector_DataGetSP",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(ASSOCIATED(DATA)) CALL FlagError("Data is already associated.",err,error,*999)     
    IF(.NOT.distributedVector%vectorFinished) CALL FlagError("The distributed vector has not been finished.",err,error,*999)
    IF(distributedVector%dataType/=MATRIX_VECTOR_SP_TYPE) THEN
      localError="The distributed data type of "// &
        & TRIM(NumberToVString(distributedVector%dataType,"*",err,error))// &
        & " does not correspond to the single precision real data type of the requested values."
      CALL FlagError(localError,err,error,*999)
    ENDIF
     
    SELECT CASE(distributedVector%libraryType)       
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissVector)
      CALL DistributedVector_CMISSVectorGet(distributedVector,cmissVector,err,error,*999)
      data=>cmissVector%dataSP
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot get values for a single precision PETSc distributed vector.",err,error,*999)          
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedVector_DataGetSP")
    RETURN
999 ERRORSEXITS("DistributedVector_DataGetSP",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_DataGetSP

  !
  !================================================================================================================================
  !

  !> Returns a pointer to the data of a double precision distributed vector. Note: the values can be used for read operations but a DistributedVector_ValuesSet call must be used to change any values. The pointer should not be deallocated.
  SUBROUTINE DistributedVector_DataGetDP(distributedVector,data,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    REAL(DP), POINTER :: data(:) !<On return, a pointer to the data of the distributed vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector
    TYPE(DistributedVectorPETScType), POINTER :: petscVector
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedVector_DataGetDP",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(ASSOCIATED(data)) CALL FlagError("Data is already associated.",err,error,*999)     
    IF(.NOT.distributedVector%vectorFinished) CALL FlagError("The distributed vector has not been finished.",err,error,*999)
    IF(distributedVector%dataType/=MATRIX_VECTOR_DP_TYPE) THEN
      localError="The distributed data type of "// &
        & TRIM(NumberToVString(distributedVector%dataType,"*",err,error))// &
        & " does not correspond to the double precision real data type of the requested values."
      CALL FlagError(localError,err,error,*999)
    ENDIF
     
    SELECT CASE(distributedVector%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissVector)
      CALL DistributedVector_CMISSVectorGet(distributedVector,cmissVector,err,error,*999)
      data=>cmissVector%dataDP
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      NULLIFY(petscVector)
      CALL DistributedVector_PETScVectorGet(distributedVector,petscVector,err,error,*999)
      IF(petscVector%useOverrideVector) THEN
        CALL Petsc_VecGetArrayReadF90(petscVector%overrideVector,data,err,error,*999)
      ELSE
        CALL Petsc_VecGetArrayReadF90(petscVector%vector,data,err,error,*999)
      ENDIF
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedVector_DataGetDP")
    RETURN
999 ERRORSEXITS("DistributedVector_DataGetDP",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_DataGetDP

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the data of a logical distributed vector. Note: the values can be used for read operations but a DistributedVector_ValuesSet call must be used to change any values. The pointer should not be deallocated.
  SUBROUTINE DistributedVector_DataGetL(distributedVector,data,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    LOGICAL, POINTER :: data(:) !<On return, a pointer to the data of the distributed vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedVector_DataGetL",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(ASSOCIATED(DATA)) CALL FlagError("Data is already associated.",err,error,*999)     
    IF(.NOT.distributedVector%vectorFinished) CALL FlagError("The distributed vector has not been finished.",err,error,*999)
    IF(distributedVector%dataType/=MATRIX_VECTOR_L_TYPE) THEN
      localError="The distributed data type of "// &
        & TRIM(NumberToVString(distributedVector%dataType,"*",err,error))// &
        & " does not correspond to the logical data type of the requested values."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(distributedVector%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissVector)
      CALL DistributedVector_CMISSVectorGet(distributedVector,cmissVector,err,error,*999)
      data=>cmissVector%dataL
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot get values for a logical PETSc distributed vector.",err,error,*999)          
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
   
    EXITS("DistributedVector_DataGetL")
    RETURN
999 ERRORSEXITS("DistributedVector_DataGetL",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_DataGetL

  !
  !================================================================================================================================
  !

  !>Restores the integer data pointer returned from DistributedVector_DataGet once the data has finished being used.
  SUBROUTINE DistributedVector_DataRestoreIntg(distributedVector,data,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    INTEGER(INTG), POINTER :: data(:) !<The a pointer to the distributed vector data for this computational node
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedVector_DataRestoreIntg",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(data)) CALL FlagError("Data is not associated.",err,error,*999)     
    IF(.NOT.distributedVector%vectorFinished) CALL FlagError("The distributed vector has not been finished.",err,error,*999)
    
    SELECT CASE(distributedVector%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(data)              
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot restore data for an integer PETSc distributed vector.",err,error,*999)
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedVector_DataRestoreIntg")
    RETURN
999 ERRORSEXITS("DistributedVector_DataRestoreIntg",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_DataRestoreIntg

  !
  !================================================================================================================================
  !

  !>Restores the single precision data pointer returned from DistributedVector_DataGet once the data has finished being used.
  SUBROUTINE DistributedVector_DataRestoreSP(distributedVector,data,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    REAL(SP), POINTER :: data(:) !<A pointer to the distributed vector data for this computational node
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedVector_DataRestoreSP",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(data)) CALL FlagError("Data is not associated.",err,error,*999)     
    IF(.NOT.distributedVector%vectorFinished) CALL FlagError("The distributed vector has not been finished.",err,error,*999)
    
    SELECT CASE(distributedVector%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(data)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot restore data for a single precision PETSc distributed vector.",err,error,*999)
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedVector_DataRestoreSP")
    RETURN
999 ERRORSEXITS("DistributedVector_DataRestoreSP",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_DataRestoreSP

  !
  !================================================================================================================================
  !

  !>Restores the double precision data pointer returned from DistributedVector_DataGet once the data has finished being used.
  SUBROUTINE DistributedVector_DataRestoreDP(distributedVector,data,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    REAL(DP), POINTER :: data(:) !<A pointer to the distributed vector data for this computational node
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedVectorPETScType), POINTER :: petscVector
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedVector_DataRestoreDP",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(data)) CALL FlagError("Data is not associated.",err,error,*999)     
    IF(.NOT.distributedVector%vectorFinished) CALL FlagError("The distributed vector has not been finished.",err,error,*999)
    
    SELECT CASE(distributedVector%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(data)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      NULLIFY(petscVector)
      CALL DistributedVector_PETScVectorGet(distributedVector,petscVector,err,error,*999)
      IF(distributedVector%petsc%useOverrideVector) THEN
        CALL Petsc_VecRestoreArrayReadF90(petscVector%overrideVector,data,err,error,*999)
      ELSE
        CALL Petsc_VecRestoreArrayReadF90(petscVector%vector,data,err,error,*999)
      ENDIF
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedVector_DataRestoreDP")
    RETURN
999 ERRORSEXITS("DistributedVector_DataRestoreDP",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_DataRestoreDP

  !
  !================================================================================================================================
  !

  !>Restores the logical data pointer returned from DistributedVector_DataGet once the data has finished being used.
  SUBROUTINE DistributedVector_DataRestoreL(distributedVector,data,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    LOGICAL, POINTER :: data(:) !<A pointer to the distributed vector data for this computational node
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("DistributedVector_DataRestoreL",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(data)) CALL FlagError("Data is not associated.",err,error,*999)     
    IF(.NOT.distributedVector%vectorFinished) CALL FlagError("The distributed vector has not been finished.",err,error,*999)
    
    SELECT CASE(distributedVector%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(data)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot restore data for a logical PETSc distributed vector.",err,error,*999)
    CASE DEFAULT
      localError="The distributed matrix library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedVector_DataRestoreL")
    RETURN
999 ERRORSEXITS("DistributedVector_DataRestoreL",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_DataRestoreL
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the ghosting type for a distributed vector
  SUBROUTINE DistributedVector_GhostingTypeSet(distributedVector,ghostingType,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector 
    INTEGER(INTG), INTENT(IN) :: ghostingType !<The ghosting type \see DISTRIBUTED_MATRIX_VECTOR_GhostingTypes,DISTRIBUTED_MATRIX_VECTOR
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector
    TYPE(DistributedVectorPETScType), POINTER :: petscVector
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: domainMapping
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedVector_GhostingTypeSet",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(distributedVector%vectorFinished) CALL FlagError("The distributed vector has already been finished.",err,error,*999)
    NULLIFY(domainMapping)
    CALL DistributedVector_RowMappingGet(distributedVector,domainMapping,err,error,*999)
    
    SELECT CASE(distributedVector%libraryType)      
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissVector)
      CALL DistributedVector_CMISSVectorGet(distributedVector,cmissVector,err,error,*999)
      SELECT CASE(ghostingType)
      CASE(DISTRIBUTED_MATRIX_VECTOR_INCLUDE_GHOSTS_TYPE)
        cmissVector%n=domainMapping%TOTAL_NUMBER_OF_LOCAL
      CASE(DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE)
        cmissVector%N=domainMapping%NUMBER_OF_LOCAL
      CASE DEFAULT
        localError="The given ghosting type of "//TRIM(NumberToVString(ghostingType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      NULLIFY(petscVector)
      CALL DistributedVector_PETScVectorGet(distributedVector,petscVector,err,error,*999)
      SELECT CASE(ghostingType)
      CASE(DISTRIBUTED_MATRIX_VECTOR_INCLUDE_GHOSTS_TYPE)
        petscVector%n=domainMapping%TOTAL_NUMBER_OF_LOCAL
      CASE(DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE)
        petscVector%n=domainMapping%NUMBER_OF_LOCAL
      CASE DEFAULT
        localError="The given ghosting type of "//TRIM(NumberToVString(ghostingType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    distributedVector%ghostingType=ghostingType
    
    EXITS("DistributedVector_GhostingTypeSet")
    RETURN
999 ERRORSEXITS("DistributedVector_GhostingTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_GhostingTypeSet

  !
  !================================================================================================================================
  !

  !>Sets/changes the library type for a distributed vector
  SUBROUTINE DistributedVector_LibraryTypeSet(distributedVector,libraryType,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector 
    INTEGER(INTG), INTENT(IN) :: libraryType !<The library type \see DISTRIBUTED_MATRIX_VECTOR_LibraryTypes,DISTRIBUTED_MATRIX_VECTOR
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,oldLibraryType
    TYPE(VARYING_STRING) :: dummyError,localError

    ENTERS("DistributedVector_LibraryTypeSet",err,error,*998)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*998)
    IF(distributedVector%vectorFinished) CALL FlagError("The distributed vector has already been finished.",err,error,*998)
    
    oldLibraryType=distributedVector%libraryType
    IF(libraryType/=oldLibraryType) THEN
      !Initialise the new library type
      SELECT CASE(libraryType)
      CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
        CALL DistributedVector_CMISSInitialise(distributedVector,err,error,*999)
      CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
        CALL DistributedVector_PETScInitialise(distributedVector,err,error,*999)
      CASE DEFAULT
        localError="The distributed vector library type of "//TRIM(NumberToVString(libraryType,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Finalise the old library type
      SELECT CASE(oldLibraryType)
      CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
        CALL DistributedVector_CMISSFinalise(distributedVector%cmiss,err,error,*999)
      CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
        CALL DistributedVector_PETScFinalise(distributedVector%petsc,err,error,*999)
      CASE DEFAULT
        localError="The distributed vector library type of "// &
          & TRIM(NumberToVString(oldLibraryType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      distributedVector%libraryType=libraryType
    ENDIF
    
    EXITS("DistributedVector_LibraryTypeSet")
    RETURN
999 SELECT CASE(libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      CALL DistributedVector_CMISSFinalise(distributedVector%cmiss,dummyErr,dummyError,*998)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL DistributedVector_PETScFinalise(distributedVector%petsc,dummyErr,dummyError,*998)
    END SELECT
998 ERRORSEXITS("DistributedVector_LibraryTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_LibraryTypeSet

  !
  !================================================================================================================================
  !

  !>Outputs a distributed vector to the specified output id.
  SUBROUTINE DistributedVector_Output(id,distributedVector,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The id of the output stream
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector to duplicate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    REAL(DP), POINTER :: vector(:)
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector
    TYPE(DistributedVectorPETScType), POINTER :: petscVector
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedVector_Output",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(.NOT.distributedVector%vectorFinished) CALL FlagError("Distributed vector has not been finished.",err,error,*999)
    
    SELECT CASE(distributedVector%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissVector)
      CALL DistributedVector_CMISSVectorGet(distributedVector,cmissVector,err,error,*999)
      SELECT CASE(distributedVector%dataType)
      CASE(DISTRIBUTED_MATRIX_VECTOR_INTG_TYPE)
        CALL WriteStringVector(id,1,1,cmissVector%n,8,8,cmissVector%dataIntg, &
          & '("Vector(:)          :",8(X,I13))','(20X,8(X,I13))',err,error,*999)
      CASE(DISTRIBUTED_MATRIX_VECTOR_SP_TYPE)
        CALL WriteStringVector(id,1,1,cmissVector%n,8,8,cmissVector%dataSP, &
          & '("Vector(:)          :",8(X,E13.6))','(20X,8(X,E13.6))',err,error,*999)
      CASE(DISTRIBUTED_MATRIX_VECTOR_DP_TYPE)
        CALL WriteStringVector(id,1,1,cmissVector%n,8,8,cmissVector%dataDP, &
          & '("Vector(:)          :",8(X,E13.6))','(20X,8(X,E13.6))',err,error,*999)
      CASE(DISTRIBUTED_MATRIX_VECTOR_L_TYPE)            
        CALL WriteStringVector(id,1,1,cmissVector%n,8,8,cmissVector%dataL, &
          & '("Vector(:)          :",8(X,L13))','(20X,8(X,L13))',err,error,*999)
      CASE(DISTRIBUTED_MATRIX_VECTOR_SPC_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(DISTRIBUTED_MATRIX_VECTOR_DPC_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The distributed vector data type of "// &
          & TRIM(NumberToVString(distributedVector%dataType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      NULLIFY(petscVector)
      CALL DistributedVector_PETScVectorGet(distributedVector,petscVector,err,error,*999)
      NULLIFY(vector)
      IF(petscVector%useOverrideVector) THEN
        CALL Petsc_VecGetArrayReadF90(petscVector%overrideVector,vector,err,error,*999)
      ELSE              
        CALL Petsc_VecGetArrayReadF90(petscVector%vector,vector,err,error,*999)
      ENDIF
      CALL WriteStringVector(id,1,1,petscVector%n,8,8,vector, &
        & '("Vector(:)          :",8(X,E13.6))','(20X,8(X,E13.6))',err,error,*999)
      IF(petscVector%useOverrideVector) THEN
        CALL Petsc_VecRestoreArrayReadF90(petscVector%overrideVector,vector,err,error,*999)
      ELSE              
        CALL Petsc_VecRestoreArrayReadF90(petscVector%vector,vector,err,error,*999)
      ENDIF
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
   
    EXITS("DistributedVector_Output")
    RETURN
999 ERRORSEXITS("DistributedVector_Output",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_Output

  !
  !================================================================================================================================
  !

  !>Sets the override vector for a distributed vector.
  SUBROUTINE DistributedVector_OverrideSetOn(distributedVector,overrideVector,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector to override
    TYPE(PetscVecType), INTENT(IN) :: overrideVector !<The override vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedVectorPETScType), POINTER :: petscVector
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedVector_OverrideSetOn",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(.NOT.distributedVector%vectorFinished) CALL FlagError("Distributed vector has not been finished.",err,error,*999)
    
    SELECT CASE(distributedVector%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)          
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      NULLIFY(petscVector)
      CALL DistributedVector_PETScVectorGet(distributedVector,petscVector,err,error,*999)
      IF(petscVector%useOverrideVector) CALL FlagError("The override vector is already set.",err,error,*999)      
      petscVector%useOverrideVector=.TRUE.
      petscVector%overrideVector=overrideVector
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedVector_OverrideSetOn")
    RETURN
999 ERRORSEXITS("DistributedVector_OverrideSetOn",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_OverrideSetOn

  !
  !================================================================================================================================
  !

  !>Turns off the override vector for a distributed vector.
  SUBROUTINE DistributedVector_OverrideSetOff(distributedVector,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector to override
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedVectorPETScType), POINTER :: petscVector
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedVector_OverrideSetOff",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(.NOT.distributedVector%vectorFinished) CALL FlagError("Distributed vector has not been finished.",err,error,*999)
    
    SELECT CASE(distributedVector%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)          
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      NULLIFY(petscVector)
      CALL DistributedVector_PETScVectorGet(distributedVector,petscVector,err,error,*999)
      IF(.NOT.petscVector%useOverrideVector) CALL FlagError("Distributed vector override is not set.",err,error,*999)      
      petscVector%useOverrideVector=.FALSE.
      CALL Petsc_VecInitialise(petscVector%overrideVector,err,error,*999)
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedVector_OverrideSetOff")
    RETURN
999 ERRORSEXITS("DistributedVector_OverrideSetOff",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_OverrideSetOff

  !
  !================================================================================================================================
  !

  !>Finishes the creation of a PETSc distributed vector.
  SUBROUTINE DistributedVector_PETScCreateFinish(petscVector,err,error,*)

    !Argument variables
    TYPE(DistributedVectorPETScType), POINTER :: petscVector !<A pointer to the distributed PETSc vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,i
    TYPE(DistributedVectorType), POINTER :: distributedVector
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: domainMapping
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("DistributedVector_PETScCreateFinish",err,error,*998)

    IF(.NOT.ASSOCIATED(petscVector)) CALL FlagError("PETSc vector is not associated.",err,error,*998)
    NULLIFY(distributedVector)
    CALL DistributedVectorPETSc_DistributedVectorGet(petscVector,distributedVector,err,error,*999)
    NULLIFY(domainMapping)
    CALL DistributedVector_RowMappingGet(distributedVector,domainMapping,err,error,*999)
    
    !Create the PETSc vector
    petscVector%dataSize=petscVector%n
    CALL Petsc_VecCreateMPI(computationalEnvironment%mpiCommunicator,petscVector%n,petscVector%globalN,petscVector%vector, &
      & err,error,*999)
    !Set up the Local to Global Mappings
    DO i=1,petscVector%n
      petscVector%globalNumbers(i)=domainMapping%LOCAL_TO_GLOBAL_MAP(i)-1
    ENDDO !i
   
    EXITS("DistributedVector_PETScCreateFinish")
    RETURN
999 CALL DistributedVector_PETScFinalise(petscVector,dummyErr,dummyError,*998)
998 ERRORSEXITS("DistributedVector_PETScCreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_PETScCreateFinish

  !
  !================================================================================================================================
  !

  !>Finalise a PETSc distributed vector.
  SUBROUTINE DistributedVector_PETScFinalise(petscVector,err,error,*)

    !Argument variables
    TYPE(DistributedVectorPETScType), POINTER :: petscVector !<A pointer to the PETSc distributed vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("DistributedVector_PETScFinalise",err,error,*999)

    IF(ASSOCIATED(petscVector)) THEN
      IF(ALLOCATED(petscVector%globalNumbers)) DEALLOCATE(petscVector%globalNumbers)
      CALL Petsc_VecFinalise(petscVector%vector,err,error,*999)
      CALL Petsc_VecFinalise(petscVector%overrideVector,err,error,*999)
      DEALLOCATE(petscVector)
    ENDIF
    
    EXITS("DistributedVector_PETScFinalise")
    RETURN
999 ERRORSEXITS("DistributedVector_PETScFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_PETScFinalise

  !
  !================================================================================================================================
  !

  !>Intialises a PETSc distributed vector.
  SUBROUTINE DistributedVector_PETScInitialise(distributedVector,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: domainMapping
    TYPE(VARYING_STRING) :: dummyError,localError
    
    ENTERS("DistributedVector_PETScInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated",err,error,*998)
    IF(ASSOCIATED(distributedVector%petsc)) &
      & CALL FlagError("PETSc is already associated for this distributed vector.",err,error,*998)
    NULLIFY(domainMapping)
    CALL DistributedVector_RowMappingGet(distributedVector,domainMapping,err,error,*999)
    
    ALLOCATE(distributedVector%petsc,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate PETSc distributed vector.",err,error,*999)
    distributedVector%petsc%distributedVector=>distributedVector
    distributedVector%libraryType=DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE
    !Set the defaults
    SELECT CASE(distributedVector%ghostingType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_INCLUDE_GHOSTS_TYPE)
      distributedVector%petsc%n=distributedVector%domainMapping%TOTAL_NUMBER_OF_LOCAL
    CASE(DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE)
      distributedVector%petsc%n=distributedVector%domainMapping%NUMBER_OF_LOCAL
    CASE DEFAULT
      localError="The distributed vector ghosting type of "// &
        & TRIM(NumberToVString(distributedVector%ghostingType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    distributedVector%petsc%globalN=distributedVector%domainMapping%NUMBER_OF_GLOBAL
    ALLOCATE(distributedVector%petsc%globalNumbers(distributedVector%petsc%n),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate PETSc distributed vector global numbers.",err,error,*999)
    distributedVector%petsc%useOverrideVector=.FALSE.
    CALL Petsc_VecInitialise(distributedVector%petsc%vector,err,error,*999)
    CALL Petsc_VecInitialise(distributedVector%petsc%overrideVector,err,error,*999)          
    
    EXITS("DistributedVector_PETScInitialise")
    RETURN
999 CALL DistributedVector_PETScFinalise(distributedVector%petsc,dummyErr,dummyError,*998)
998 ERRORSEXITS("DistributedVector_PETScInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_PETScInitialise

  !
  !================================================================================================================================
  !

  !>Finalises a CMISS distributed vector transfer information and deallocates all memory.
  SUBROUTINE DistributedVector_CMISSTransferFinalise(cmissVector,domainIdx,err,error,*)

    !Argument variables
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector !<A pointer to the CMISS distributed vector
    INTEGER(INTG), INTENT(IN) :: domainIdx !<The domain index of the distributed vector to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedVector_CMISSTransferFinalise",err,error,*999)

    IF(ASSOCIATED(cmissVector)) THEN
      IF(ALLOCATED(cmissVector%transfers)) THEN
        IF(domainIdx<=0.OR.domainIdx>SIZE(cmissVector%transfers,1)) THEN
          localError="The domain index of "//TRIM(NumberToVString(domainIdx,"*",err,error))// &
            & " is invalid. It must be between 1 and "//TRIM(NumberToVString(SIZE(cmissVector%transfers,1),"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        NULLIFY(cmissVector%transfers(domainIdx)%cmissVector)
        cmissVector%transfers(domainIdx)%dataType=0
        cmissVector%transfers(domainIdx)%receiveTagNumber=-1
        cmissVector%transfers(domainIdx)%sendTagNumber=-1
        cmissVector%transfers(domainIdx)%sendBufferSize=0
        cmissVector%transfers(domainIdx)%receiveBufferSize=0
        cmissVector%transfers(domainIdx)%mpiSendRequest=MPI_REQUEST_NULL
        cmissVector%transfers(domainIdx)%mpiReceiveRequest=MPI_REQUEST_NULL
        IF(ALLOCATED(cmissVector%transfers(domainIdx)%sendBufferIntg)) &
          & DEALLOCATE(cmissVector%transfers(domainIdx)%sendBufferIntg)
        IF(ALLOCATED(cmissVector%transfers(domainIdx)%sendBufferSP)) &
          & DEALLOCATE(cmissVector%transfers(domainIdx)%sendBufferSP)
        IF(ALLOCATED(cmissVector%transfers(domainIdx)%sendBufferDP)) &
          & DEALLOCATE(cmissVector%transfers(domainIdx)%sendBufferDP)
        IF(ALLOCATED(cmissVector%transfers(domainIdx)%sendBufferL)) &
          & DEALLOCATE(cmissVector%transfers(domainIdx)%sendBufferL)
        IF(ALLOCATED(cmissVector%transfers(domainIdx)%receiveBufferIntg)) &
          & DEALLOCATE(cmissVector%transfers(domainIdx)%receiveBufferIntg)
        IF(ALLOCATED(cmissVector%transfers(domainIdx)%receiveBufferSP)) &
          & DEALLOCATE(cmissVector%transfers(domainIdx)%receiveBufferSP)
        IF(ALLOCATED(cmissVector%transfers(domainIdx)%receiveBufferDP)) &
          & DEALLOCATE(cmissVector%transfers(domainIdx)%receiveBufferDP)
        IF(ALLOCATED(cmissVector%transfers(domainIdx)%receiveBufferL)) &
          & DEALLOCATE(cmissVector%transfers(domainIdx)%receiveBufferL)
      ENDIF
    ENDIF
    
    EXITS("DistributedVector_CMISSTransferFinalise")
    RETURN
999 ERRORSEXITS("DistributedVector_CMISSTransferFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_CMISSTransferFinalise

  !
  !================================================================================================================================
  !

  !>Initialises a CMISS distributed vector transfer information.
  SUBROUTINE DistributedVector_CMISSTransferInitialise(cmissVector,domainIdx,err,error,*)

    !Argument variables
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector !<A pointer to the CMISS distributed vector
    INTEGER(INTG), INTENT(IN) :: domainIdx !<The domain index to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedVector_CMISSTransferInitialise",err,error,*999)

    IF(.NOT.ASSOCIATED(cmissVector)) CALL FlagError("CMISS vector is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(cmissVector%transfers)) CALL FlagError("CMISS vector transfers is not allocated.",err,error,*999)
    IF(domainIdx<=0.OR.domainIdx>SIZE(cmissVector%transfers,1)) THEN
      localError="The domain index of "//TRIM(NumberToVString(domainIdx,"*",err,error))// &
        & " is invalid. It must be between 1 and "// &
        & TRIM(NumberToVString(SIZE(cmissVector%transfers,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    cmissVector%transfers(domainIdx)%cmissVector=>cmissVector
    cmissVector%transfers(domainIdx)%dataType=0
    cmissVector%transfers(domainIdx)%sendBufferSize=0
    cmissVector%transfers(domainIdx)%receiveBufferSize=0
    cmissVector%transfers(domainIdx)%sendTagNumber=-1
    cmissVector%transfers(domainIdx)%receiveTagNumber=-1
    cmissVector%transfers(domainIdx)%mpiSendRequest=MPI_REQUEST_NULL
    cmissVector%transfers(domainIdx)%mpiReceiveRequest=MPI_REQUEST_NULL
    
    EXITS("DistributedVector_CMISSTransferInitialise")
    RETURN
999 ERRORSEXITS("DistributedVector_CMISSTransferInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_CMISSTransferInitialise

  !
  !================================================================================================================================
  !

  !>Finishes the (ghost) update procedure for a distributed vector. This routine will wait until all transfers have completed!
  SUBROUTINE DistributedVector_UpdateFinish(distributedVector,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: domainIdx,i,numberOfComputationalNodes
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector
    TYPE(DistributedVectorPETScType), POINTER :: petscVector
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: domainMapping
    TYPE(VARYING_STRING) :: localError
   
    ENTERS("DistributedVector_UpdateFinish",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(.NOT.distributedVector%vectorFinished) CALL FlagError("The distributed vector has not been finished.",err,error,*999)
    NULLIFY(domainMapping)
    CALL DistributedVector_RowMappingGet(distributedVector,domainMapping,err,error,*999)
    
    SELECT CASE(distributedVector%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissVector)
      CALL DistributedVector_CMISSVectorGet(distributedVector,cmissVector,err,error,*999)
      numberOfComputationalNodes=ComputationalEnvironment_NumberOfNodesGet(err,error)
      IF(err/=0) GOTO 999
      IF(numberOfComputationalNodes>1) THEN
        CALL DistributedVector_UpdateWaitFinished(distributedVector,err,error,*999)
        !Copy the receive buffers back to the ghost positions in the data vector
        SELECT CASE(distributedVector%dataType)
        CASE(DISTRIBUTED_MATRIX_VECTOR_INTG_TYPE)
          DO domainIdx=1,domainMapping%NUMBER_OF_ADJACENT_DOMAINS
            DO i=1,domainMapping%ADJACENT_DOMAINS(domainIdx)%NUMBER_OF_RECEIVE_GHOSTS
              cmissVector%dataIntg(domainMapping%ADJACENT_DOMAINS(domainIdx)%LOCAL_GHOST_RECEIVE_INDICES(i))= &
                & cmissVector%transfers(domainIdx)%receiveBufferIntg(i)
            ENDDO !i
          ENDDO !domainIdx
        CASE(DISTRIBUTED_MATRIX_VECTOR_SP_TYPE)
          DO domainIdx=1,domainMapping%NUMBER_OF_ADJACENT_DOMAINS
            DO i=1,domainMapping%ADJACENT_DOMAINS(domainIdx)%NUMBER_OF_RECEIVE_GHOSTS
              cmissVector%dataSP(domainMapping%ADJACENT_DOMAINS(domainIdx)%LOCAL_GHOST_RECEIVE_INDICES(i))= &
                & cmissVector%transfers(domainIdx)%receiveBufferSP(i)
            ENDDO !i
          ENDDO !domainIdx
        CASE(DISTRIBUTED_MATRIX_VECTOR_DP_TYPE)
          DO domainIdx=1,domainMapping%NUMBER_OF_ADJACENT_DOMAINS
            DO i=1,domainMapping%ADJACENT_DOMAINS(domainIdx)%NUMBER_OF_RECEIVE_GHOSTS
              cmissVector%dataDP(domainMapping%ADJACENT_DOMAINS(domainIdx)%LOCAL_GHOST_RECEIVE_INDICES(i))= &
                & cmissVector%transfers(domainIdx)%receiveBufferDP(i)
            ENDDO !i
          ENDDO !domainIdx
        CASE(DISTRIBUTED_MATRIX_VECTOR_L_TYPE)
          DO domainIdx=1,domainMapping%NUMBER_OF_ADJACENT_DOMAINS
            DO i=1,domainMapping%ADJACENT_DOMAINS(domainIdx)%NUMBER_OF_RECEIVE_GHOSTS
              cmissVector%dataL(domainMapping%ADJACENT_DOMAINS(domainIdx)%LOCAL_GHOST_RECEIVE_INDICES(i))= &
                & cmissVector%transfers(domainIdx)%receiveBufferL(i)
            ENDDO !i
          ENDDO !domainIdx
        CASE(DISTRIBUTED_MATRIX_VECTOR_SPC_TYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(DISTRIBUTED_MATRIX_VECTOR_DPC_TYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE DEFAULT
          localError="The distributed vector data type of "// &
            & TRIM(NumberToVString(distributedVector%dataType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDIF
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      NULLIFY(petscVector)
      CALL DistributedVector_PETScVectorGet(distributedVector,petscVector,err,error,*999)
      IF(petscVector%useOverrideVector) THEN
        CALL Petsc_VecAssemblyEnd(petscVector%overrideVector,err,error,*999)
      ELSE
        CALL Petsc_VecAssemblyEnd(petscVector%vector,err,error,*999)
      ENDIF
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    IF(diagnostics1) THEN
      SELECT CASE(distributedVector%libraryType)
      CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Distributed vector :",err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Data type = ",distributedVector%dataType,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Base tag number = ",cmissVector%baseTagNumber,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of adjacent domains = ",domainMapping%NUMBER_OF_ADJACENT_DOMAINS, &
          & err,error,*999)
        DO domainIdx=1,domainMapping%NUMBER_OF_ADJACENT_DOMAINS
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Domain idx = ",domainIdx,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Domain number = ",domainMapping%ADJACENT_DOMAINS(domainIdx)% &
            & DOMAIN_NUMBER,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Receive tag number = ",cmissVector%transfers(domainIdx)% &
            & receiveTagNumber,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Send tag number = ",cmissVector%transfers(domainIdx)% &
            & sendTagNumber,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    MPI send request = ",cmissVector%transfers(domainIdx)% &
            & mpiSendRequest,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    MPI receive request = ",cmissVector%transfers(domainIdx)% &
            & mpiReceiveRequest,err,error,*999)
        ENDDO !domainIdx
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Data size = ",cmissVector%dataSize,err,error,*999)
        SELECT CASE(distributedVector%dataType)
        CASE(DISTRIBUTED_MATRIX_VECTOR_INTG_TYPE)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,cmissVector%dataSize,5,5,cmissVector%dataIntg, &
            & '("  Data :",5(X,I13))','(8X,5(X,I13))',err,error,*999)
        CASE(DISTRIBUTED_MATRIX_VECTOR_SP_TYPE)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,cmissVector%dataSize,5,5,cmissVector%dataSP, &
            & '("  Data :",5(X,E13.6))','(8X,5(X,E13.6))',err,error,*999)
        CASE(DISTRIBUTED_MATRIX_VECTOR_DP_TYPE)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,cmissVector%dataSize,5,5,cmissVector%dataDP, &
            & '("  Data :",5(X,E13.6))','(8X,5(X,E13.6))',err,error,*999)
        CASE(DISTRIBUTED_MATRIX_VECTOR_L_TYPE)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,cmissVector%dataSize,8,8,cmissVector%dataL, &
            & '("  Data :",8(X,L))','(8X,8(X,L))',err,error,*999)
        CASE(DISTRIBUTED_MATRIX_VECTOR_SPC_TYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(DISTRIBUTED_MATRIX_VECTOR_DPC_TYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE DEFAULT
          localError="The distributed vector data type of "// &
            & TRIM(NumberToVString(distributedVector%dataType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
        !Do nothing
      CASE DEFAULT
        localError="The distributed vector library type of "// &
          & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDIF
    
    EXITS("DistributedVector_UpdateFinish")
    RETURN
999 ERRORSEXITS("DistributedVector_UpdateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_UpdateFinish

  !
  !================================================================================================================================
  !

  !>Tests to see if a distributed vector update has finised! \todo USE MPI_TESTALL and store the request handles as big array.
  SUBROUTINE DistributedVector_UpdateIsFinished(distributedVector,isFinished,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    LOGICAL, INTENT(OUT) :: isFinished !<On return, is .TRUE. if all the transfer operations for the distributed vector have completed, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: domainIdx
    INTEGER(INTG) :: mpiIError,status(MPI_STATUS_SIZE)
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: domainMapping
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedVector_UpdateIsFinished",err,error,*999)

    isFinished=.FALSE.
    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(.NOT.distributedVector%vectorFinished) CALL FlagError("The distributed vector has not been finished.",err,error,*999)
    NULLIFY(domainMapping)
    CALL DistributedVector_RowMappingGet(distributedVector,domainMapping,err,error,*999)
    
    SELECT CASE(distributedVector%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissVector)
      CALL DistributedVector_CMISSVectorGet(distributedVector,cmissVector,err,error,*999)
!!TODO: USE MPI_TESTALL and store the request handles as big array.
      DO domainIdx=1,domainMapping%NUMBER_OF_ADJACENT_DOMAINS
        CALL MPI_TEST(cmissVector%transfers(domainIdx)%mpiReceiveRequest,isFinished,status,mpiIError)
        CALL MPI_ErrorCheck("MPI_TEST",mpiIError,err,error,*999)
        IF(.NOT.isFinished) EXIT
        !CALL MPI_TEST(cmissVector%transfers(domainIdx)%mpiSendRequest,isFinished,status,mpiIError)
        !CALL MPI_ErrorCheck("MPI_TEST",mpiIError,err,error,*999)
        !IF(.NOT.isFinished) EXIT
      ENDDO !domainIdx
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot test if update isfinished for a PETSc distributed vector.",err,error,*999)          
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedVector_UpdateIsFinished")
    RETURN
999 ERRORSEXITS("DistributedVector_UpdateIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_UpdateIsFinished

  !
  !================================================================================================================================
  !

  !>Waits until a distributed vector update has finised
  SUBROUTINE DistributedVector_UpdateWaitFinished(distributedVector,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: domainIdx
    INTEGER(INTG) :: mpiIError,status(MPI_STATUS_SIZE)
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: domainMapping
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedVector_UpdateWaitFinished",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(.NOT.distributedVector%vectorFinished) CALL FlagError("The distributed vector has not been finished.",err,error,*999)
    NULLIFY(domainMapping)
    CALL DistributedVector_RowMappingGet(distributedVector,domainMapping,err,error,*999)
    
    SELECT CASE(distributedVector%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissVector)
      CALL DistributedVector_CMISSVectorGet(distributedVector,cmissVector,err,error,*999)
!!TODO: USE MPI_WAITALL and store the request handles as big array.
      DO domainIdx=1,domainMapping%NUMBER_OF_ADJACENT_DOMAINS
        CALL MPI_WAIT(cmissVector%transfers(domainIdx)%mpiReceiveRequest,status,mpiIError)
        CALL MPI_ErrorCheck("MPI_WAIT",mpiIError,err,error,*999)
      ENDDO !domainIdx
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot wait for finished for a PETSc distributed vector.",err,error,*999)          
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
     
    EXITS("DistributedVector_UpdateWaitFinished")
    RETURN
999 ERRORSEXITS("DistributedVector_UpdateWaitFinished",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_UpdateWaitFinished

  !
  !================================================================================================================================
  !

  !>Starts the (ghost) update procedure for a distributed vector.
  SUBROUTINE DistributedVector_UpdateStart(distributedVector,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: domainIdx,i,mpiIError,numberOfComputationalNodes
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector
    TYPE(DistributedVectorPETScType), POINTER :: petscVector
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: domainMapping
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedVector_UpdateStart",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(.NOT.distributedVector%vectorFinished) CALL FlagError("The distributed vector has not been finished.",err,error,*999)
    NULLIFY(domainMapping)
    CALL DistributedVector_RowMappingGet(distributedVector,domainMapping,err,error,*999)
    
    SELECT CASE(distributedVector%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissVector)
      CALL DistributedVector_CMISSVectorGet(distributedVector,cmissVector,err,error,*999)
      numberOfComputationalNodes=ComputationalEnvironment_NumberOfNodesGet(err,error)
      IF(err/=0) GOTO 999
      IF(numberOfComputationalNodes>1) THEN
        IF(domainMapping%NUMBER_OF_ADJACENT_DOMAINS>0) THEN
          !Fill in the send buffers with the send ghost values
          DO domainIdx=1,domainMapping%NUMBER_OF_ADJACENT_DOMAINS
            SELECT CASE(distributedVector%dataType)
            CASE(DISTRIBUTED_MATRIX_VECTOR_INTG_TYPE)
              DO i=1,domainMapping%ADJACENT_DOMAINS(domainIdx)%NUMBER_OF_SEND_GHOSTS
                cmissVector%transfers(domainIdx)%sendBufferIntg(i)= &
                  & cmissVector%dataIntg(domainMapping%ADJACENT_DOMAINS(domainIdx)%LOCAL_GHOST_SEND_INDICES(i))
              ENDDO !i
            CASE(DISTRIBUTED_MATRIX_VECTOR_SP_TYPE)
              DO i=1,domainMapping%ADJACENT_DOMAINS(domainIdx)%NUMBER_OF_SEND_GHOSTS
                cmissVector%transfers(domainIdx)%sendBufferSP(i)= &
                  & cmissVector%dataSP(domainMapping%ADJACENT_DOMAINS(domainIdx)%LOCAL_GHOST_SEND_INDICES(i))
              ENDDO !i
            CASE(DISTRIBUTED_MATRIX_VECTOR_DP_TYPE)
              DO i=1,domainMapping%ADJACENT_DOMAINS(domainIdx)%NUMBER_OF_SEND_GHOSTS
                cmissVector%transfers(domainIdx)%sendBufferDP(i)= &
                  & cmissVector%dataDP(domainMapping%ADJACENT_DOMAINS(domainIdx)%LOCAL_GHOST_SEND_INDICES(i))
              ENDDO !i
            CASE(DISTRIBUTED_MATRIX_VECTOR_L_TYPE)
              DO i=1,domainMapping%ADJACENT_DOMAINS(domainIdx)%NUMBER_OF_SEND_GHOSTS
                cmissVector%transfers(domainIdx)%sendBufferL(i)= &
                  & cmissVector%dataL(domainMapping%ADJACENT_DOMAINS(domainIdx)%LOCAL_GHOST_SEND_INDICES(i))
              ENDDO !i
            CASE(DISTRIBUTED_MATRIX_VECTOR_SPC_TYPE)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(DISTRIBUTED_MATRIX_VECTOR_DPC_TYPE)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE DEFAULT
              localError="The distributed vector data type of "// &
                & TRIM(NumberToVString(distributedVector%dataType,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ENDDO !domainIdx
          !Post all the receive calls first and then the send calls.
          DO domainIdx=1,domainMapping%NUMBER_OF_ADJACENT_DOMAINS
            SELECT CASE(distributedVector%dataType)
            CASE(DISTRIBUTED_MATRIX_VECTOR_INTG_TYPE)
              CALL MPI_IRECV(cmissVector%transfers(domainIdx)%receiveBufferIntg, &
                & cmissVector%transfers(domainIdx)%receiveBufferSize,MPI_INTEGER, &
                & domainMapping%ADJACENT_DOMAINS(domainIdx)%DOMAIN_NUMBER, &
                & cmissVector%transfers(domainIdx)%receiveTagNumber, &
                & computationalEnvironment%mpiCommunicator, &
                & cmissVector%transfers(domainIdx)%mpiReceiveRequest,mpiIError)
              CALL MPI_ErrorCheck("MPI_IRECV",mpiIError,err,error,*999)
              IF(diagnostics5) THEN
                CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"MPI IRECV call posted:",err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Receive count = ",cmissVector%transfers(domainIdx)% &
                  & receiveBufferSize,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Receive datatype = ",MPI_INTEGER,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Receive source = ",domainMapping%ADJACENT_DOMAINS(domainIdx)% &
                  & DOMAIN_NUMBER,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Receive tag = ",cmissVector%transfers(domainIdx)% &
                  & receiveTagNumber,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Receive comm = ",computationalEnvironment%mpiCommunicator, &
                  & err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Receive request = ",cmissVector%transfers(domainIdx)% &
                  & mpiReceiveRequest,err,error,*999)
              ENDIF
            CASE(DISTRIBUTED_MATRIX_VECTOR_SP_TYPE)
              CALL MPI_IRECV(cmissVector%transfers(domainIdx)%receiveBufferSP, &
                & cmissVector%transfers(domainIdx)%receiveBufferSize,MPI_REAL, &
                & domainMapping%ADJACENT_DOMAINS(domainIdx)%DOMAIN_NUMBER, &
                & cmissVector%transfers(domainIdx)%receiveTagNumber, &
                & computationalEnvironment%mpiCommunicator, &
                & cmissVector%transfers(domainIdx)%mpiReceiveRequest,mpiIError)
              CALL MPI_ErrorCheck("MPI_IRECV",mpiIError,err,error,*999)
              IF(diagnostics5) THEN
                CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"MPI IRECV call posted:",err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Receive count = ",cmissVector%transfers(domainIdx)% &
                  & receiveBufferSize,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Receive datatype = ",MPI_REAL,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Receive source = ",domainMapping%ADJACENT_DOMAINS(domainIdx)% &
                  & DOMAIN_NUMBER,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Receive tag = ",cmissVector%transfers(domainIdx)% &
                  & receiveTagNumber,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Receive comm = ",computationalEnvironment%mpiCommunicator, &
                  & err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Receive request = ",cmissVector%transfers(domainIdx)% &
                  & mpiReceiveRequest,err,error,*999)
              ENDIF
            CASE(DISTRIBUTED_MATRIX_VECTOR_DP_TYPE)
              CALL MPI_IRECV(cmissVector%transfers(domainIdx)%receiveBufferDP, &
                & cmissVector%transfers(domainIdx)%receiveBufferSize,MPI_DOUBLE_PRECISION, &
                & domainMapping%ADJACENT_DOMAINS(domainIdx)%DOMAIN_NUMBER, &
                & cmissVector%transfers(domainIdx)%receiveTagNumber, &
                & computationalEnvironment%mpiCommunicator, &
                & cmissVector%transfers(domainIdx)%mpiReceiveRequest,mpiIError)
              CALL MPI_ErrorCheck("MPI_IRECV",mpiIError,err,error,*999)
              IF(diagnostics5) THEN
                CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"MPI IRECV call posted:",err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Receive count = ",cmissVector%transfers(domainIdx)% &
                  & receiveBufferSize,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Receive datatype = ",MPI_DOUBLE_PRECISION,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Receive source = ",domainMapping%ADJACENT_DOMAINS(domainIdx)% &
                  & DOMAIN_NUMBER,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Receive tag = ",cmissVector%transfers(domainIdx)%receiveTagNumber, &
                  & err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Receive comm = ",computationalEnvironment%mpiCommunicator, &
                  & err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Receive request = ",cmissVector%transfers(domainIdx)% &
                  & mpiReceiveRequest,err,error,*999)
              ENDIF
            CASE(DISTRIBUTED_MATRIX_VECTOR_L_TYPE)
              CALL MPI_IRECV(cmissVector%transfers(domainIdx)%receiveBufferL, &
                & cmissVector%transfers(domainIdx)%receiveBufferSize,MPI_LOGICAL, &
                & domainMapping%ADJACENT_DOMAINS(domainIdx)%DOMAIN_NUMBER, &
                & cmissVector%transfers(domainIdx)%receiveTagNumber, &
                & computationalEnvironment%mpiCommunicator, &
                & cmissVector%transfers(domainIdx)%mpiReceiveRequest,mpiIError)
              CALL MPI_ErrorCheck("MPI_IRECV",mpiIError,err,error,*999)
              IF(diagnostics5) THEN
                CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"MPI IRECV call posted:",err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Receive count = ",cmissVector%transfers(domainIdx)% &
                  & receiveBufferSize,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Receive datatype = ",MPI_LOGICAL,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Receive source = ",domainMapping%ADJACENT_DOMAINS(domainIdx)% &
                  & DOMAIN_NUMBER,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Receive tag = ",cmissVector%transfers(domainIdx)% &
                  & receiveTagNumber,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Receive comm = ",computationalEnvironment%mpiCommunicator, &
                  & err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Receive request = ",cmissVector%transfers(domainIdx)% &
                  & mpiReceiveRequest,err,error,*999)
              ENDIF
            CASE(DISTRIBUTED_MATRIX_VECTOR_SPC_TYPE)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(DISTRIBUTED_MATRIX_VECTOR_DPC_TYPE)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE DEFAULT
              localError="The distributed vector data type of "// &
                & TRIM(NumberToVString(distributedVector%dataType,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ENDDO !domainIdx
          !Post all the send calls.
          DO domainIdx=1,domainMapping%NUMBER_OF_ADJACENT_DOMAINS
            SELECT CASE(distributedVector%dataType)
            CASE(DISTRIBUTED_MATRIX_VECTOR_INTG_TYPE)
              CALL MPI_ISEND(cmissVector%transfers(domainIdx)%sendBufferIntg, &
                & cmissVector%transfers(domainIdx)%sendBufferSize,MPI_INTEGER, &
                & domainMapping%ADJACENT_DOMAINS(domainIdx)%DOMAIN_NUMBER, &
                & cmissVector%transfers(domainIdx)%sendTagNumber, &
                & computationalEnvironment%mpiCommunicator, &
                & cmissVector%transfers(domainIdx)%mpiSendRequest,mpiIError)
              CALL MPI_ErrorCheck("MPI_ISEND",mpiIError,err,error,*999)
              IF(diagnostics5) THEN
                CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"MPI ISEND call posted:",err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Send count = ",cmissVector%transfers(domainIdx)% &
                  & sendBufferSize,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Send datatype = ",MPI_INTEGER,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Send dest = ",domainMapping%ADJACENT_DOMAINS(domainIdx)% &
                  & DOMAIN_NUMBER,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Send tag = ",cmissVector%transfers(domainIdx)%sendTagNumber, &
                  & err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Send comm = ",computationalEnvironment%mpiCommunicator, &
                  & err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Send request = ",cmissVector%transfers(domainIdx)% &
                  & mpiSendRequest,err,error,*999)
              ENDIF
            CASE(DISTRIBUTED_MATRIX_VECTOR_SP_TYPE)
              CALL MPI_ISEND(cmissVector%transfers(domainIdx)%sendBufferSP, &
                & cmissVector%transfers(domainIdx)%sendBufferSize,MPI_REAL, &
                & domainMapping%ADJACENT_DOMAINS(domainIdx)%DOMAIN_NUMBER, &
                & cmissVector%transfers(domainIdx)%sendTagNumber, &
                & computationalEnvironment%mpiCommunicator, &
                & cmissVector%transfers(domainIdx)%mpiSendRequest,mpiIError)
              CALL MPI_ErrorCheck("MPI_ISEND",mpiIError,err,error,*999)
              IF(diagnostics5) THEN
                CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"MPI ISEND call posted:",err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Send count = ",cmissVector%transfers(domainIdx)% &
                  & sendBufferSize,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Send datatype = ",MPI_REAL,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Send dest = ",domainMapping%ADJACENT_DOMAINS(domainIdx)% &
                  & DOMAIN_NUMBER,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Send tag = ",cmissVector%transfers(domainIdx)%sendTagNumber, &
                  & err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Send comm = ",computationalEnvironment%mpiCommunicator, &
                  & err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Send request = ",cmissVector%transfers(domainIdx)% &
                  & mpiSendRequest,err,error,*999)
              ENDIF
            CASE(DISTRIBUTED_MATRIX_VECTOR_DP_TYPE)
              CALL MPI_ISEND(cmissVector%transfers(domainIdx)%sendBufferDP, &
                & cmissVector%transfers(domainIdx)%sendBufferSize,MPI_DOUBLE_PRECISION, &
                & domainMapping%ADJACENT_DOMAINS(domainIdx)%DOMAIN_NUMBER, &
                & cmissVector%transfers(domainIdx)%sendTagNumber, &
                & computationalEnvironment%mpiCommunicator, &
                & cmissVector%transfers(domainIdx)%mpiSendRequest,mpiIError)
              CALL MPI_ErrorCheck("MPI_ISEND",mpiIError,err,error,*999)
              IF(diagnostics5) THEN
                CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"MPI ISEND call posted:",err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Send count = ",cmissVector%transfers(domainIdx)% &
                  & sendBufferSize,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Send datatype = ",MPI_DOUBLE_PRECISION,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Send dest = ",domainMapping%ADJACENT_DOMAINS(domainIdx)% &
                  & DOMAIN_NUMBER,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Send tag = ",cmissVector%transfers(domainIdx)%sendTagNumber, &
                  & err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Send comm = ",computationalEnvironment%mpiCommunicator, &
                  & err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Send request = ",cmissVector%transfers(domainIdx)%mpiSendRequest, &
                  & err,error,*999)
              ENDIF
            CASE(DISTRIBUTED_MATRIX_VECTOR_L_TYPE)
              CALL MPI_ISEND(cmissVector%transfers(domainIdx)%sendBufferL, &
                & cmissVector%transfers(domainIdx)%sendBufferSize,MPI_LOGICAL, &
                & domainMapping%ADJACENT_DOMAINS(domainIdx)%DOMAIN_NUMBER, &
                & cmissVector%transfers(domainIdx)%sendTagNumber, &
                & computationalEnvironment%mpiCommunicator, &
                & cmissVector%transfers(domainIdx)%mpiSendRequest,mpiIError)
              CALL MPI_ErrorCheck("MPI_ISEND",mpiIError,err,error,*999)
              IF(diagnostics5) THEN
                CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"MPI ISEND call posted:",err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Send count = ",cmissVector%transfers(domainIdx)%sendBufferSize, &
                  & err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Send datatype = ",MPI_LOGICAL,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Send dest = ",domainMapping%ADJACENT_DOMAINS(domainIdx)% &
                  & DOMAIN_NUMBER,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Send tag = ",cmissVector%transfers(domainIdx)%sendTagNumber, &
                  & err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Send comm = ",computationalEnvironment%mpiCommunicator, &
                  & err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Send request = ",cmissVector%transfers(domainIdx)%mpiSendRequest, &
                  & err,error,*999)
              ENDIF
            CASE(DISTRIBUTED_MATRIX_VECTOR_SPC_TYPE)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(DISTRIBUTED_MATRIX_VECTOR_DPC_TYPE)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE DEFAULT
              localError="The distributed vector data type of "// &
                & TRIM(NumberToVString(distributedVector%dataType,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ENDDO !domainIdx
        ENDIF
      ENDIF
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)          
      NULLIFY(petscVector)
      CALL DistributedVector_PETScVectorGet(distributedVector,petscVector,err,error,*999)
      IF(petscVector%useOverrideVector) THEN
        CALL Petsc_VecAssemblyBegin(petscVector%overrideVector,err,error,*999)
      ELSE
        CALL Petsc_VecAssemblyBegin(petscVector%vector,err,error,*999)
      ENDIF
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    IF(diagnostics1) THEN
      SELECT CASE(distributedVector%libraryType)
      CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Distributed vector :",err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Data type = ",distributedVector%dataType,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Base tag number = ",cmissVector%baseTagNumber,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of adjacent domains = ",domainMapping%NUMBER_OF_ADJACENT_DOMAINS, &
          & err,error,*999)
        DO domainIdx=1,domainMapping%NUMBER_OF_ADJACENT_DOMAINS
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Domain idx = ",domainIdx,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Domain number = ",domainMapping%ADJACENT_DOMAINS(domainIdx)% &
            & DOMAIN_NUMBER,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Receive tag number = ",cmissVector%transfers(domainIdx)% &
            & receiveTagNumber,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Send tag number = ",cmissVector%transfers(domainIdx)%sendTagNumber, &
            & err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    MPI send request = ",cmissVector%transfers(domainIdx)%mpiSendRequest, &
            & err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    MPI receive request = ",cmissVector%transfers(domainIdx)% &
            & mpiReceiveRequest,err,error,*999)
        ENDDO !domainIdx
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Data size = ",cmissVector%dataSize,err,error,*999)
        SELECT CASE(distributedVector%dataType)
        CASE(DISTRIBUTED_MATRIX_VECTOR_INTG_TYPE)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,cmissVector%dataSize,5,5,cmissVector%dataIntg, &
            & '("  Data :",5(X,I13))','(8X,5(X,I13))',err,error,*999)
        CASE(DISTRIBUTED_MATRIX_VECTOR_SP_TYPE)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,cmissVector%dataSize,5,5,cmissVector%dataSP, &
            & '("  Data :",5(X,E13.6))','(8X,5(X,E13.6))',err,error,*999)
        CASE(DISTRIBUTED_MATRIX_VECTOR_DP_TYPE)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,cmissVector%dataSize,5,5,cmissVector%dataDP, &
            & '("  Data :",5(X,E13.6))','(8X,5(X,E13.6))',err,error,*999)
        CASE(DISTRIBUTED_MATRIX_VECTOR_L_TYPE)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,cmissVector%dataSize,8,8,cmissVector%dataL, &
            & '("  Data :",8(X,L))','(8X,8(X,L))',err,error,*999)
        CASE(DISTRIBUTED_MATRIX_VECTOR_SPC_TYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(DISTRIBUTED_MATRIX_VECTOR_DPC_TYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE DEFAULT
          localError="The distributed vector data type of "// &
            & TRIM(NumberToVString(distributedVector%dataType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
        !Do nothing
      CASE DEFAULT
        localError="The distributed vector library type of "// &
          & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDIF
    
    EXITS("DistributedVector_UpdateStart")
    RETURN
999 ERRORSEXITS("DistributedVector_UpdateStart",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_UpdateStart

  !
  !================================================================================================================================
  !

  !>Calculates the L2 norm of a distributed vector values on this computational node
  SUBROUTINE DistributedVector_L2Norm(distributedVector,norm,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), INTENT(IN), POINTER :: distributedVector !<A pointer to the distributed vector
    REAL(DP), INTENT(OUT) :: norm !<The L2 norm of values from this computational node
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedVector_L2Norm",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(.NOT.distributedVector%vectorFinished) CALL FlagError("The distributed vector has not been finished.",err,error,*999)
    
    SELECT CASE(distributedVector%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissVector)
      CALL DistributedVector_CMISSVectorGet(distributedVector,cmissVector,err,error,*999)
      SELECT CASE(distributedVector%dataType)
      CASE(DISTRIBUTED_MATRIX_VECTOR_DP_TYPE)
        CALL L2Norm(cmissVector%dataDP(1:cmissVector%dataSize),norm,err,error,*999)
      CASE(DISTRIBUTED_MATRIX_VECTOR_SP_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(DISTRIBUTED_MATRIX_VECTOR_INTG_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(DISTRIBUTED_MATRIX_VECTOR_L_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(DISTRIBUTED_MATRIX_VECTOR_SPC_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(DISTRIBUTED_MATRIX_VECTOR_DPC_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The distributed data type of "// &
          & TRIM(NumberToVString(distributedVector%dataType,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot calculate norm for a PETSc distributed vector.",err,error,*999)
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("DistributedVector_L2Norm")
    RETURN
999 ERRORSEXITS("DistributedVector_L2Norm",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_L2Norm
  
  !
  !================================================================================================================================
  !

  !>Calculates the dot product of 2 distributed integer vectors on this computational node
  SUBROUTINE DistributedVector_DotProductIntg(distributedVectorA,distributedVectorB,dotProduct,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), INTENT(IN), POINTER :: distributedVectorA !<A pointer to the distributed vector A
    TYPE(DistributedVectorType), INTENT(IN), POINTER :: distributedVectorB !<A pointer to the distributed vector B
    INTEGER(INTG), INTENT(OUT) :: dotProduct !<The dot product on this computational node
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: i
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVectorA,cmissVectorB
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedVector_DotProductIntg",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVectorA)) CALL FlagError("Distributed vector A is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(distributedVectorB)) CALL FlagError("Distributed vector B is not associated.",err,error,*999)    
    IF(.NOT.distributedVectorA%vectorFinished) CALL FlagError("Distributed vector A has not been finished.",err,error,*999)
    IF(.NOT.distributedVectorB%vectorFinished) CALL FlagError("Distributed vector B has not been finished.",err,error,*999)    
    IF(distributedVectorA%libraryType/=distributedVectorB%libraryType) THEN
      localError="The distributed vector A library type of "// &
        & TRIM(NumberToVString(distributedVectorA%libraryType,"*",err,error))// &
        & " does not match the distributed vector B library type of "// &
        & TRIM(NumberToVString(distributedVectorB%libraryType,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(distributedVectorA%dataType/=distributedVectorB%dataType) THEN
      localError="The distributed vector A data type of "// &
        & TRIM(NumberToVString(distributedVectorA%dataType,"*",err,error))// &
        & " does not match the distributed vector B data type of "// &
        & TRIM(NumberToVString(distributedVectorB%dataType,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(distributedVectorA%dataType/=MATRIX_VECTOR_INTG_TYPE) THEN
      localError="The distributed vectors data type of "// &
        & TRIM(NumberToVString(distributedVectorA%dataType,"*",err,error))// &
        & " does not match the integer data type of the supplied dot product value."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(distributedVectorA%libraryType)      
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissVectorA)
      CALL DistributedVector_CMISSVectorGet(distributedVectorA,cmissVectorA,err,error,*999)
      NULLIFY(cmissVectorB)
      CALL DistributedVector_CMISSVectorGet(distributedVectorB,cmissVectorB,err,error,*999)
      IF(cmissVectorA%dataSize/=cmissVectorB%dataSize) THEN
        localError="The distributed vector A data size of "// &
          & TRIM(NumberToVString(cmissVectorA%dataSize,"*",err,error))// &
          & " does not match the distributed vector B data size of "// &
          & TRIM(NumberToVString(cmissVectorB%dataSize,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      dotProduct=0
      DO i=1,cmissVectorA%dataSize
        dotProduct=dotProduct+cmissVectorA%dataIntg(i)*cmissVectorB%dataIntg(i)
      ENDDO !i
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot computer integer dot product for PETSc distributed vectors.",err,error,*999)
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVectorA%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("DistributedVector_DotProductIntg")
    RETURN
999 ERRORSEXITS("DistributedVector_DotProductIntg",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_DotProductIntg
  
  !
  !================================================================================================================================
  !

  !>Calculates the dot product of 2 distributed single-precision vectors on this computational node
  SUBROUTINE DistributedVector_DotProductSP(distributedVectorA,distributedVectorB,dotProduct,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), INTENT(IN), POINTER :: distributedVectorA !<A pointer to the distributed vector A
    TYPE(DistributedVectorType), INTENT(IN), POINTER :: distributedVectorB !<A pointer to the distributed vector B
    REAL(SP), INTENT(OUT) :: dotProduct !<The dot product on this computational node
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: i
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVectorA,cmissVectorB
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedVector_DotProductSP",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVectorA)) CALL FlagError("Distributed vector A is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(distributedVectorB)) CALL FlagError("Distributed vector B is not associated.",err,error,*999)    
    IF(.NOT.distributedVectorA%vectorFinished) CALL FlagError("Distributed vector A has not been finished.",err,error,*999)
    IF(.NOT.distributedVectorB%vectorFinished) CALL FlagError("Distributed vector B has not been finished.",err,error,*999)    
    IF(distributedVectorA%libraryType/=distributedVectorB%libraryType) THEN
      localError="The distributed vector A library type of "// &
        & TRIM(NumberToVString(distributedVectorA%libraryType,"*",err,error))// &
        & " does not match the distributed vector B library type of "// &
        & TRIM(NumberToVString(distributedVectorB%libraryType,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(distributedVectorA%dataType/=distributedVectorB%dataType) THEN
      localError="The distributed vector A data type of "// &
        & TRIM(NumberToVString(distributedVectorA%dataType,"*",err,error))// &
        & " does not match the distributed vector B data type of "// &
        & TRIM(NumberToVString(distributedVectorB%dataType,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(distributedVectorA%dataType/=MATRIX_VECTOR_SP_TYPE) THEN
      localError="The distributed vectors data type of "// &
        & TRIM(NumberToVString(distributedVectorA%dataType,"*",err,error))// &
        & " does not match the single precision real data type of the supplied dot product value."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(distributedVectorA%libraryType)      
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissVectorA)
      CALL DistributedVector_CMISSVectorGet(distributedVectorA,cmissVectorA,err,error,*999)
      NULLIFY(cmissVectorB)
      CALL DistributedVector_CMISSVectorGet(distributedVectorB,cmissVectorB,err,error,*999)
      IF(cmissVectorA%dataSize/=cmissVectorB%dataSize) THEN
        localError="The distributed vector A data size of "// &
          & TRIM(NumberToVString(cmissVectorA%dataSize,"*",err,error))// &
          & " does not match the distributed vector B data size of "// &
          & TRIM(NumberToVString(cmissVectorB%dataSize,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      dotProduct=0.0_SP
      DO i=1,cmissVectorA%dataSize
        dotProduct=dotProduct+cmissVectorA%dataSP(i)*cmissVectorB%dataSP(i)
      ENDDO !i
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot computer single precision real dot product for PETSc distributed vectors.",err,error,*999)
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVectorA%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("DistributedVector_DotProductSP")
    RETURN
999 ERRORSEXITS("DistributedVector_DotProductSP",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_DotProductSP

  !
  !================================================================================================================================
  !

  !>Calculates the dot product of 2 distributed double-precision vectors on this computational node
  SUBROUTINE DistributedVector_DotProductDp(distributedVectorA,distributedVectorB,dotProduct,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), INTENT(IN), POINTER :: distributedVectorA !<A pointer to the distributed vector A
    TYPE(DistributedVectorType), INTENT(IN), POINTER :: distributedVectorB !<A pointer to the distributed vector B
    REAL(DP), INTENT(OUT) :: dotProduct !<The dot product on this computational node
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: i
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVectorA,cmissVectorB
    TYPE(DistributedVectorPETScType), POINTER :: petscVectorA,petscVectorB
    TYPE(VARYING_STRING) :: localError

    ENTERS("DistributedVector_DotProductDp",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVectorA)) CALL FlagError("Distributed vector A is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(distributedVectorB)) CALL FlagError("Distributed vector B is not associated.",err,error,*999)    
    IF(.NOT.distributedVectorA%vectorFinished) CALL FlagError("Distributed vector A has not been finished.",err,error,*999)
    IF(.NOT.distributedVectorB%vectorFinished) CALL FlagError("Distributed vector B has not been finished.",err,error,*999)    
    IF(distributedVectorA%libraryType/=distributedVectorB%libraryType) THEN
      localError="The distributed vector A library type of "// &
        & TRIM(NumberToVString(distributedVectorA%libraryType,"*",err,error))// &
        & " does not match the distributed vector B library type of "// &
        & TRIM(NumberToVString(distributedVectorB%libraryType,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(distributedVectorA%dataType/=distributedVectorB%dataType) THEN
      localError="The distributed vector A data type of "// &
        & TRIM(NumberToVString(distributedVectorA%dataType,"*",err,error))// &
        & " does not match the distributed vector B data type of "// &
        & TRIM(NumberToVString(distributedVectorB%dataType,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(distributedVectorA%dataType/=MATRIX_VECTOR_DP_TYPE) THEN
      localError="The distributed vectors data type of "// &
        & TRIM(NumberToVString(distributedVectorA%dataType,"*",err,error))// &
        & " does not match the double precision real data type of the supplied dot product value."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(distributedVectorA%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissVectorA)
      CALL DistributedVector_CMISSVectorGet(distributedVectorA,cmissVectorA,err,error,*999)
      NULLIFY(cmissVectorB)
      CALL DistributedVector_CMISSVectorGet(distributedVectorB,cmissVectorB,err,error,*999)
      IF(cmissVectorA%dataSize/=cmissVectorB%dataSize) THEN
        localError="The distributed vector A data size of "// &
          & TRIM(NumberToVString(cmissVectorA%dataSize,"*",err,error))// &
          & " does not match the distributed vector B data size of "// &
          & TRIM(NumberToVString(cmissVectorB%dataSize,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      dotProduct=0.0_DP
      DO i=1,cmissVectorA%dataSize
        dotProduct=dotProduct+cmissVectorA%dataDP(i)*cmissVectorB%dataDP(i)
      ENDDO !i
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      NULLIFY(petscVectorA)
      CALL DistributedVector_PETScVectorGet(distributedVectorA,petscVectorA,err,error,*999)
      NULLIFY(petscVectorB)
      CALL DistributedVector_PETScVectorGet(distributedVectorB,petscVectorB,err,error,*999)
      CALL Petsc_VecDot(petscVectorA%vector,petscVectorB%vector,dotProduct,err,error,*999)
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVectorA%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("DistributedVector_DotProductDp")
    RETURN
999 ERRORSEXITS("DistributedVector_DotProductDp",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_DotProductDp
  
  !
  !================================================================================================================================
  !

  !>Adds values to a distributed integer vector.
  SUBROUTINE DistributedVector_ValuesAddIntg(distributedVector,indices,values,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    INTEGER(INTG), INTENT(IN) :: indices(:) !<indices(i). The i'th index to add
    INTEGER(INTG), INTENT(IN) :: values(:) !<values(i). The i'th value to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: i
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedVector_ValuesAddIntg",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(.NOT.distributedVector%vectorFinished) CALL FlagError("The distributed vector has not been finished.",err,error,*999)
    IF(SIZE(indices,1)/=SIZE(values,1)) THEN
      localError="The size of the indicies array of "//TRIM(NumberToVString(SIZE(indices,1),"*",err,error))// &
        & " does not conform to the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(distributedVector%dataType/=MATRIX_VECTOR_INTG_TYPE) THEN
      localError="The distributed data type of "// &
        & TRIM(NumberToVString(distributedVector%dataType,"*",err,error))// &
        & " does not correspond to the integer data type of the given values."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(distributedVector%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissVector)
      CALL DistributedVector_CMISSVectorGet(distributedVector,cmissVector,err,error,*999)
      DO i=1,SIZE(indices,1)
        !Allow all values added until dof mappings fixed. Ghost values that are added will not be propogated
        IF(indices(i)<=0.OR.indices(i)>cmissVector%dataSize) THEN
          localError="The index at position "//TRIM(NumberToVString(i,"*",err,error))//" of "// &
            & TRIM(NumberToVString(indices(i),"*",err,error))//" is invalid. The index must be between 1 and "// &
            & TRIM(NumberToVString(cmissVector%dataSize,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        cmissVector%dataIntg(indices(i))=cmissVector%dataIntg(indices(i))+values(i)
      ENDDO !i
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot add values for an integer PETSc distributed vector.",err,error,*999)                    
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedVector_ValuesAddIntg")
    RETURN
999 ERRORSEXITS("DistributedVector_ValuesAddIntg",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_ValuesAddIntg

  !
  !================================================================================================================================
  !

  !>Adds one value to a distributed integer vector.
  SUBROUTINE DistributedVector_ValuesAddIntg1(distributedVector,index,value,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    INTEGER(INTG), INTENT(IN) :: index !<The index to be added at
    INTEGER(INTG), INTENT(IN) :: value !<The value to be added
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedVector_ValuesAddIntg1",err,error,*999)
    
    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(.NOT.distributedVector%vectorFinished) CALL FlagError("The distributed vector has not been finished.",err,error,*999)
    IF(distributedVector%dataType/=MATRIX_VECTOR_SP_TYPE) THEN
      localError="The distributed data type of "// &
        & TRIM(NumberToVString(distributedVector%dataType,"*",err,error))// &
        & " does not correspond to the integer data type of the given values."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(distributedVector%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissVector)
      CALL DistributedVector_CMISSVectorGet(distributedVector,cmissVector,err,error,*999)
      !Allow all values to be added until dof mappings fixed. Ghost values that are added will not be propogated
      IF(index<=0.OR.index>cmissVector%dataSize) THEN
        localError="Index "//TRIM(NumberToVString(index,"*",err,error))// &
          & " is invalid. The index must be between 1 and "// &
          & TRIM(NumberToVString(cmissVector%dataSize,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      cmissVector%dataIntg(index)=cmissVector%dataIntg(index)+value
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot add values for an integer PETSc distributed vector.",err,error,*999)          
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedVector_ValuesAddIntg1")
    RETURN
999 ERRORSEXITS("DistributedVector_ValuesAddIntg1",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_ValuesAddIntg1

  !
  !================================================================================================================================
  !

  !>Adds values to a distributed single precision vector.
  SUBROUTINE DistributedVector_ValuesAddSP(distributedVector,indices,values,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    INTEGER(INTG), INTENT(IN) :: indices(:) !<indices(i). The i'th index to be added
    REAL(SP), INTENT(IN) :: values(:) !<values(i). The i'th value to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: i
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedVector_ValuesAddSP",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(.NOT.distributedVector%vectorFinished) CALL FlagError("The distributed vector has not been finished.",err,error,*999)
    IF(SIZE(indices,1)/=SIZE(values,1)) THEN
      localError="The size of the indicies array of "//TRIM(NumberToVString(SIZE(indices,1),"*",err,error))// &
        & " does not conform to the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(distributedVector%dataType/=MATRIX_VECTOR_SP_TYPE) THEN
      localError="The distributed data type of "// &
        & TRIM(NumberToVString(distributedVector%dataType,"*",err,error))// &
        & " does not correspond to the single precision real data type of the given values."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(distributedVector%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissVector)
      CALL DistributedVector_CMISSVectorGet(distributedVector,cmissVector,err,error,*999)
      DO i=1,SIZE(indices,1)
        !Allow all values to be added until dof mappings fixed. Ghost values that are added will not be propogated
        IF(indices(i)<=0.OR.indices(i)>cmissVector%dataSize) THEN
          localError="The index at position "//TRIM(NumberToVString(i,"*",err,error))//" of "// &
            & TRIM(NumberToVString(indices(i),"*",err,error))//" is invalid. The index must be between 1 and "// &
            & TRIM(NumberToVString(cmissVector%dataSize,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        cmissVector%dataSP(indices(i))=cmissVector%dataSP(indices(i))+values(i)
      ENDDO !i
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot add values for a single precision PETSc distributed vector.",err,error,*999)          
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedVector_ValuesAddSP")
    RETURN
999 ERRORSEXITS("DistributedVector_ValuesAddSP",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_ValuesAddSP

  !
  !================================================================================================================================
  !

  !>Adds one value to a distributed single precision vector.
  SUBROUTINE DistributedVector_ValuesAddSP1(distributedVector,index,value,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    INTEGER(INTG), INTENT(IN) :: index !<The index to be added
    REAL(SP), INTENT(IN) :: value !<The value to be added
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedVector_ValuesAddSP1",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(.NOT.distributedVector%vectorFinished) CALL FlagError("The distributed vector has not been finished.",err,error,*999)
    IF(distributedVector%dataType/=MATRIX_VECTOR_SP_TYPE) THEN
      localError="The distributed data type of "// &
        & TRIM(NumberToVString(distributedVector%dataType,"*",err,error))// &
        & " does not correspond to the single precision real data type of the given values."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(distributedVector%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissVector)
      CALL DistributedVector_CMISSVectorGet(distributedVector,cmissVector,err,error,*999)
      !Allow all values to be added until dof mappings fixed. Ghost values that are added will not be propogated
      IF(index<=0.OR.index>cmissVector%dataSize) THEN
        localError="Index "//TRIM(NumberToVString(index,"*",err,error))// &
          & " is invalid. The index must be between 1 and "// &
          & TRIM(NumberToVString(cmissVector%dataSize,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      cmissVector%dataSP(index)=cmissVector%dataSP(index)+VALUE
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot add values for a single precision PETSc distributed vector.",err,error,*999)          
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedVector_ValuesAddSP1")
    RETURN
999 ERRORSEXITS("DistributedVector_ValuesAddSP1",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_ValuesAddSP1

  !
  !================================================================================================================================
  !

  !>Adds values to a distributed double precision vector.
  SUBROUTINE DistributedVector_ValuesAddDP(distributedVector,indices,VALUES,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    INTEGER(INTG), INTENT(IN) :: indices(:) !<indices(i). The i'th index to be added
    REAL(DP), INTENT(IN) :: values(:) !<VALUES(i). The i'th value to added
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: i
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector
    TYPE(DistributedVectorPETscType), POINTER :: petscVector
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedVector_ValuesAddDP",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(.NOT.distributedVector%vectorFinished) CALL FlagError("The distributed vector has not been finished.",err,error,*999)
    IF(SIZE(indices,1)/=SIZE(values,1)) THEN
      localError="The size of the indicies array of "//TRIM(NumberToVString(SIZE(indices,1),"*",err,error))// &
        & " does not conform to the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(distributedVector%dataType/=MATRIX_VECTOR_DP_TYPE) THEN
      localError="The distributed data type of "// &
        & TRIM(NumberToVString(distributedVector%dataType,"*",err,error))// &
        & " does not correspond to the double precision real data type of the given values."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(distributedVector%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissVector)
      CALL DistributedVector_CMISSVectorGet(distributedVector,cmissVector,err,error,*999)
      DO i=1,SIZE(indices,1)
        !Allow all values to be added until dof mappings fixed. Ghost values that are added will not be propogated
        IF(indices(i)<=0.OR.indices(i)>cmissVector%dataSize) THEN
          localError="The index at position "//TRIM(NumberToVString(i,"*",err,error))//" of "// &
            & TRIM(NumberToVString(indices(i),"*",err,error))//" is invalid. The index must be between 1 and "// &
            & TRIM(NumberToVString(cmissVector%dataSize,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        cmissVector%dataDP(indices(i))=cmissVector%dataDP(indices(i))+VALUES(i)
      ENDDO !i
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      NULLIFY(petscVector)
      CALL DistributedVector_PETscVectorGet(distributedVector,petscVector,err,error,*999)
      IF(petscVector%useOverrideVector) THEN
        CALL Petsc_VecSetValues(petscVector%overrideVector,SIZE(indices,1),petscVector%globalNumbers(indices), &
          & values,PETSC_ADD_VALUES,err,error,*999)
      ELSE
        CALL Petsc_VecSetValues(petscVector%vector,SIZE(indices,1),petscVector%globalNumbers(indices), &
          & values,PETSC_ADD_VALUES,err,error,*999)
      ENDIF
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedVector_ValuesAddDP")
    RETURN
999 ERRORSEXITS("DistributedVector_ValuesAddDP",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_ValuesAddDP

  !
  !================================================================================================================================
  !

  !>Adds one value to a distributed double precision vector.
  SUBROUTINE DistributedVector_ValuesAddDP1(distributedVector,index,value,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    INTEGER(INTG), INTENT(IN) :: index !<The index to be added
    REAL(DP), INTENT(IN) :: value !<The value to be added
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    REAL(DP) :: petscValue(1)
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector
    TYPE(DistributedVectorPETscType), POINTER :: petscVector
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedVector_ValuesAddDP1",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(.NOT.distributedVector%vectorFinished) CALL FlagError("The distributed vector has not been finished.",err,error,*999)
    IF(distributedVector%dataType/=MATRIX_VECTOR_DP_TYPE) THEN
      localError="The distributed data type of "// &
        & TRIM(NumberToVString(distributedVector%dataType,"*",err,error))// &
        & " does not correspond to the double precision real data type of the given values."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(distributedVector%libraryType)  
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissVector)
      CALL DistributedVector_CMISSVectorGet(distributedVector,cmissVector,err,error,*999)
      !Allow all values to be added until dof mappings fixed. Ghost values that are added will not be propogated
      IF(index<=0.OR.index>cmissVector%dataSize) THEN
        localError="Index "//TRIM(NumberToVString(index,"*",err,error))// &
          & " is invalid. The index must be between 1 and "// &
          & TRIM(NumberToVString(cmissVector%dataSize,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      cmissVector%dataDP(index)=cmissVector%dataDP(index)+value
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      NULLIFY(petscVector)
      CALL DistributedVector_PETscVectorGet(distributedVector,petscVector,err,error,*999)
      petscValue(1)=value
      IF(petscVector%useOverrideVector) THEN
        CALL Petsc_VecSetValues(petscVector%overrideVector,1,petscVector%globalNumbers(index), &
          & petscValue,PETSC_ADD_VALUES,err,error,*999)
      ELSE
        CALL Petsc_VecSetValues(petscVector%vector,1,petscVector%globalNumbers(index), &
          & petscValue,PETSC_ADD_VALUES,err,error,*999)
      ENDIF
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedVector_ValuesAddDP1")
    RETURN
999 ERRORSEXITS("DistributedVector_ValuesAddDP1",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_ValuesAddDP1

  !
  !================================================================================================================================
  ! 

  !>Adds values to a distributed logical vector.
  SUBROUTINE DistributedVector_ValuesAddL(distributedVector,indices,VALUES,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    INTEGER(INTG), INTENT(IN) :: indices(:) !<indices(i). The i'th index to be added
    LOGICAL, INTENT(IN) :: values(:) !<VALUES(i). The i'th value to added
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: i
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedVector_ValuesAddL",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(.NOT.distributedVector%vectorFinished) CALL FlagError("The distributed vector has not been finished.",err,error,*999)
    IF(SIZE(indices,1)/=SIZE(values,1)) THEN
      localError="The size of the indicies array of "//TRIM(NumberToVString(SIZE(indices,1),"*",err,error))// &
        & " does not conform to the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(distributedVector%dataType/=MATRIX_VECTOR_L_TYPE) THEN
      localError="The distributed data type of "// &
        & TRIM(NumberToVString(distributedVector%dataType,"*",err,error))// &
        & " does not correspond to the logical data type of the given values."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(distributedVector%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissVector)
      CALL DistributedVector_CMISSVectorGet(distributedVector,cmissVector,err,error,*999)
      DO i=1,SIZE(indices,1)
        !Allow all values to be added until dof mappings fixed. Ghost values that are added will not be propogated
        IF(indices(i)<=0.OR.indices(i)>cmissVector%dataSize) THEN
          localError="The index at position "//TRIM(NumberToVString(i,"*",err,error))//" of "// &
            & TRIM(NumberToVString(indices(i),"*",err,error))//" is invalid. The index must be between 1 and "// &
            & TRIM(NumberToVString(cmissVector%dataSize,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        cmissVector%dataL(indices(i))=cmissVector%dataL(indices(i)).OR.values(i)
      ENDDO !i
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot add values for a logical PETSc distributed vector.",err,error,*999)          
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedVector_ValuesAddL")
    RETURN
999 ERRORSEXITS("DistributedVector_ValuesAddL",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_ValuesAddL

  !
  !================================================================================================================================
  !

  !>Adds one value to a distributed logical vector.
  SUBROUTINE DistributedVector_ValuesAddL1(distributedVector,index,value,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    INTEGER(INTG), INTENT(IN) :: index !<The index to be added
    LOGICAL, INTENT(IN) :: value !<The value to be added
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedVector_ValuesAddL1",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(.NOT.distributedVector%vectorFinished) CALL FlagError("The distributed vector has not been finished.",err,error,*999)
    IF(distributedVector%dataType/=MATRIX_VECTOR_L_TYPE) THEN
      localError="The distributed data type of "// &
        & TRIM(NumberToVString(distributedVector%dataType,"*",err,error))// &
        & " does not correspond to the logical data type of the given values."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(distributedVector%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissVector)
      CALL DistributedVector_CMISSVectorGet(distributedVector,cmissVector,err,error,*999)
      !Allow all values to be added until dof mappings fixed. Ghost values that are added will not be propogated
      IF(index<=0.OR.index>cmissVector%dataSize) THEN
        localError="Index "//TRIM(NumberToVString(index,"*",err,error))// &
          & " is invalid. The index must be between 1 and "// &
          & TRIM(NumberToVString(cmissVector%dataSize,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      cmissVector%dataL(index)=cmissVector%dataL(index).OR.value
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot add values for a logical PETSc distributed vector.",err,error,*999)          
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
   
    EXITS("DistributedVector_ValuesAddL1")
    RETURN
999 ERRORSEXITS("DistributedVector_ValuesAddL1",err,error)
    RETURN 1

  END SUBROUTINE DistributedVector_ValuesAddL1

  !
  !================================================================================================================================
  !

  !>Gets values in a distributed integer vector.
  SUBROUTINE DistributedVector_ValuesGetIntg(distributedVector,indices,values,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    INTEGER(INTG), INTENT(IN) :: indices(:) !<indices(i). The i'th index to be get
    INTEGER(INTG), INTENT(OUT) :: values(:) !<values(i). On return, the value at the i'th specified index
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: i
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedVector_ValuesGetIntg",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(.NOT.distributedVector%vectorFinished) CALL FlagError("The distributed vector has not been finished.",err,error,*999)
    IF(SIZE(indices,1)/=SIZE(values,1)) THEN
      localError="The size of the indicies array of "//TRIM(NumberToVString(SIZE(indices,1),"*",err,error))// &
        & " does not conform to the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(distributedVector%dataType/=MATRIX_VECTOR_INTG_TYPE) THEN
      localError="The distributed data type of "// &
        & TRIM(NumberToVString(distributedVector%dataType,"*",err,error))// &
        & " does not correspond to the integer data type of the given values."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(distributedVector%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissVector)
      CALL DistributedVector_CMISSVectorGet(distributedVector,cmissVector,err,error,*999)
      DO i=1,SIZE(indices,1)
        IF(indices(i)<=0.OR.indices(i)>cmissVector%dataSize) THEN
          localError="The index at position "//TRIM(NumberToVString(i,"*",err,error))//" of "// &
            & TRIM(NumberToVString(indices(i),"*",err,error))//" is invalid. The index must be between 1 and "// &
            & TRIM(NumberToVString(cmissVector%dataSize,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        values(i)=cmissVector%dataIntg(indices(i))
      ENDDO !i
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot set values for an integer PETSc distributed vector.",err,error,*999)          
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedVector_ValuesGetIntg")
    RETURN
999 ERRORSEXITS("DistributedVector_ValuesGetIntg",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_ValuesGetIntg

  !
  !================================================================================================================================
  !

  !>Gets one value in a distributed integer vector.
  SUBROUTINE DistributedVector_ValuesGetIntg1(distributedVector,index,value,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    INTEGER(INTG), INTENT(IN) :: index !<The index to be get
    INTEGER(INTG), INTENT(OUT) :: value !<On return, the value at the specified index
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedVector_ValuesGetIntg1",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(.NOT.distributedVector%vectorFinished) CALL FlagError("The distributed vector has not been finished.",err,error,*999)
    IF(distributedVector%dataType/=MATRIX_VECTOR_INTG_TYPE) THEN
      localError="The distributed data type of "// &
        & TRIM(NumberToVString(distributedVector%dataType,"*",err,error))// &
        & " does not correspond to the integer data type of the given values."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(distributedVector%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissVector)
      CALL DistributedVector_CMISSVectorGet(distributedVector,cmissVector,err,error,*999)
      IF(index<=0.OR.index>cmissVector%dataSize) THEN
        localError="Index "//TRIM(NumberToVString(index,"*",err,error))// &
          & " is invalid. The index must be between 1 and "// &
          & TRIM(NumberToVString(cmissVector%dataSize,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      value=cmissVector%dataIntg(index)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot set values for an integer PETSc distributed vector.",err,error,*999)          
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedVector_ValuesGetIntg1")
    RETURN
999 ERRORSEXITS("DistributedVector_ValuesGetIntg1",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_ValuesGetIntg1

  !
  !================================================================================================================================
  !

  !>Gets values in a distributed single precision vector.
  SUBROUTINE DistributedVector_ValuesGetSP(distributedVector,indices,values,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    INTEGER(INTG), INTENT(IN) :: indices(:) !<indices(i). The i'th index to be get
    REAL(SP), INTENT(OUT) :: values(:) !<values(i). On return, the value at the i'th specified index
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: i
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedVector_ValuesGetSP",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(.NOT.distributedVector%vectorFinished) CALL FlagError("The distributed vector has not been finished.",err,error,*999)
    IF(SIZE(indices,1)/=SIZE(values,1)) THEN
      localError="The size of the indicies array of "//TRIM(NumberToVString(SIZE(indices,1),"*",err,error))// &
        & " does not conform to the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(distributedVector%dataType/=MATRIX_VECTOR_SP_TYPE) THEN
      localError="The distributed data type of "// &
        & TRIM(NumberToVString(distributedVector%dataType,"*",err,error))// &
        & " does not correspond to the single precision real data type of the given values."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(distributedVector%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissVector)
      CALL DistributedVector_CMISSVectorGet(distributedVector,cmissVector,err,error,*999)
      DO i=1,SIZE(indices,1)
        IF(indices(i)<=0.OR.indices(i)>cmissVector%dataSize) THEN
          localError="The index at position "//TRIM(NumberToVString(i,"*",err,error))//" of "// &
            & TRIM(NumberToVString(indices(i),"*",err,error))//" is invalid. The index must be between 1 and "// &
            & TRIM(NumberToVString(cmissVector%dataSize,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        values(i)=cmissVector%dataSP(indices(i))
      ENDDO !i
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot get values for a single precision PETSc distributed vector.",err,error,*999)          
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedVector_ValuesGetSP")
    RETURN
999 ERRORSEXITS("DistributedVector_ValuesGetSP",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_ValuesGetSP

  !
  !================================================================================================================================
  !

  !>Gets one value in a distributed single precision vector.
  SUBROUTINE DistributedVector_ValuesGetSP1(distributedVector,index,value,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    INTEGER(INTG), INTENT(IN) :: index !<The index to be get
    REAL(SP), INTENT(OUT) :: value !<On return, the value at the specified index
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedVector_ValuesGetSP1",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(.NOT.distributedVector%vectorFinished) CALL FlagError("The distributed vector has not been finished.",err,error,*999)
    IF(distributedVector%dataType/=MATRIX_VECTOR_SP_TYPE) THEN
      localError="The distributed data type of "// &
        & TRIM(NumberToVString(distributedVector%dataType,"*",err,error))// &
        & " does not correspond to the single precision real data type of the given values."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(distributedVector%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissVector)
      CALL DistributedVector_CMISSVectorGet(distributedVector,cmissVector,err,error,*999)
      IF(index<=0.OR.index>cmissVector%dataSize) THEN
        localError="Index "//TRIM(NumberToVString(index,"*",err,error))// &
          & " is invalid. The index must be between 1 and "// &
          & TRIM(NumberToVString(cmissVector%dataSize,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      value=cmissVector%dataSP(index)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot set values for a single precision PETSc distributed vector.",err,error,*999)          
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedVector_ValuesGetSP1")
    RETURN
999 ERRORSEXITS("DistributedVector_ValuesGetSP1",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_ValuesGetSP1

  !
  !================================================================================================================================
  !

  !>Gets values in a distributed double precision vector.
  SUBROUTINE DistributedVector_ValuesGetDP(distributedVector,indices,values,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    INTEGER(INTG), INTENT(IN) :: indices(:) !<indices(i). The i'th index to be get
    REAL(DP), INTENT(OUT) :: values(:) !<values(i). On return, the value at the i'th specified index
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: i,petscIndices(SIZE(indices,1))
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector
    TYPE(DistributedVectorPETScType), POINTER :: petscVector
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: domainMapping
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedVector_ValuesGetDP",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(.NOT.distributedVector%vectorFinished) CALL FlagError("The distributed vector has not been finished.",err,error,*999)
    IF(SIZE(indices,1)/=SIZE(values,1)) THEN
      localError="The size of the indicies array of "//TRIM(NumberToVString(SIZE(indices,1),"*",err,error))// &
        & " does not conform to the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(distributedVector%dataType/=MATRIX_VECTOR_DP_TYPE) THEN
      localError="The distributed data type of "// &
        & TRIM(NumberToVString(distributedVector%dataType,"*",err,error))// &
        & " does not correspond to the double precision real data type of the given values."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(distributedVector%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissVector)
      CALL DistributedVector_CMISSVectorGet(distributedVector,cmissVector,err,error,*999)
      DO i=1,SIZE(indices,1)
        IF(indices(i)<=0.OR.indices(i)>cmissVector%dataSize) THEN
          localError="The index at position "//TRIM(NumberToVString(i,"*",err,error))//" of "// &
            & TRIM(NumberToVString(indices(i),"*",err,error))//" is invalid. The index must be between 1 and "// &
            & TRIM(NumberToVString(cmissVector%dataSize,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        values(i)=cmissVector%dataDP(indices(i))
      ENDDO !i
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      NULLIFY(petscVector)
      CALL DistributedVector_PETScVectorGet(distributedVector,petscVector,err,error,*999)
      NULLIFY(domainMapping)
      CALL DistributedVector_RowMappingGet(distributedVector,domainMapping,err,error,*999)
      DO i=1,SIZE(indices,1)
        petscIndices(i)=domainMapping%LOCAL_TO_GLOBAL_MAP(indices(i))-1 !PETSc uses global 0-based indices
      ENDDO !i
      IF(petscVector%useOverrideVector) THEN
        CALL Petsc_VecGetValues(petscVector%overrideVector,SIZE(indices,1),petscIndices,values,err,error,*999)
      ELSE
        CALL Petsc_VecGetValues(petscVector%vector,SIZE(indices,1),petscIndices,values,err,error,*999)
      ENDIF
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
   
    EXITS("DistributedVector_ValuesGetDP")
    RETURN
999 ERRORSEXITS("DistributedVector_ValuesGetDP",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_ValuesGetDP

  !
  !================================================================================================================================
  !

  !>Gets one value in a distributed double precision vector.
  SUBROUTINE DistributedVector_ValuesGetDP1(distributedVector,index,value,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    INTEGER(INTG), INTENT(IN) :: index !<The index to be get
    REAL(DP), INTENT(OUT) :: value !<On return, the value at the specified index
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: petscIndex(1)
    REAL(DP) :: petscValue(1)
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector
    TYPE(DistributedVectorPETScType), POINTER :: petscVector
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: domainMapping
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedVector_ValuesGetDP1",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(.NOT.distributedVector%vectorFinished) CALL FlagError("The distributed vector has not been finished.",err,error,*999)
    IF(distributedVector%dataType/=MATRIX_VECTOR_DP_TYPE) THEN
      localError="The distributed data type of "// &
        & TRIM(NumberToVString(distributedVector%dataType,"*",err,error))// &
        & " does not correspond to the double precision real data type of the given values."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(distributedVector%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissVector)
      CALL DistributedVector_CMISSVectorGet(distributedVector,cmissVector,err,error,*999)
      IF(index<=0.OR.index>cmissVector%dataSize) THEN
        localError="Index "//TRIM(NumberToVString(index,"*",err,error))// &
          & " is invalid. The index must be between 1 and "// &
          & TRIM(NumberToVString(cmissVector%dataSize,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      value=cmissVector%dataDP(index)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      NULLIFY(petscVector)
      CALL DistributedVector_PETScVectorGet(distributedVector,petscVector,err,error,*999)
      NULLIFY(domainMapping)
      CALL DistributedVector_RowMappingGet(distributedVector,domainMapping,err,error,*999)
      petscIndex=domainMapping%LOCAL_TO_GLOBAL_MAP(index)-1 !PETSc uses global 0-based indices
      IF(petscVector%useOverrideVector) THEN
        CALL Petsc_VecGetValues(petscVector%overrideVector,1,petscIndex,petscValue,err,error,*999)
      ELSE
        CALL Petsc_VecGetValues(petscVector%vector,1,petscIndex,petscValue,err,error,*999)
      ENDIF
      value=petscValue(1)
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedVector_ValuesGetDP1")
    RETURN
999 ERRORSEXITS("DistributedVector_ValuesGetDP1",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_ValuesGetDP1

  !
  !================================================================================================================================
  !

  !>Gets values in a distributed logical vector.
  SUBROUTINE DistributedVector_ValuesGetL(distributedVector,indices,values,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    INTEGER(INTG), INTENT(IN) :: indices(:) !<indices(i). The i'th index to be get
    LOGICAL, INTENT(OUT) :: values(:) !<values(i). On return, the value in the i'th specified index
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: i
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedVector_ValuesGetL",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(.NOT.distributedVector%vectorFinished) CALL FlagError("The distributed vector has not been finished.",err,error,*999)
    IF(SIZE(indices,1)/=SIZE(values,1)) THEN
      localError="The size of the indicies array of "//TRIM(NumberToVString(SIZE(indices,1),"*",err,error))// &
        & " does not conform to the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(distributedVector%dataType/=MATRIX_VECTOR_L_TYPE) THEN
      localError="The distributed data type of "// &
        & TRIM(NumberToVString(distributedVector%dataType,"*",err,error))// &
        & " does not correspond to the logical data type of the given values."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(distributedVector%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissVector)
      CALL DistributedVector_CMISSVectorGet(distributedVector,cmissVector,err,error,*999)
      DO i=1,SIZE(indices,1)
        IF(indices(i)<=0.OR.indices(i)>cmissVector%dataSize) THEN
          localError="The index at position "//TRIM(NumberToVString(i,"*",err,error))//" of "// &
            & TRIM(NumberToVString(indices(i),"*",err,error))//" is invalid. The index must be between 1 and "// &
            & TRIM(NumberToVString(cmissVector%dataSize,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        values(i)=cmissVector%dataL(indices(i))
      ENDDO !i
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot set values for a logical PETSc distributed vector.",err,error,*999)          
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedVector_ValuesGetL")
    RETURN
999 ERRORSEXITS("DistributedVector_ValuesGetL",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_ValuesGetL

  !
  !================================================================================================================================
  !

  !>Gets one value in a distributed logical vector.
  SUBROUTINE DistributedVector_ValuesGetL1(distributedVector,index,value,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    INTEGER(INTG), INTENT(IN) :: index !<The index to be get
    LOGICAL, INTENT(OUT) :: value !<On return, the value in the specified index
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedVector_ValuesGetL1",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(.NOT.distributedVector%vectorFinished) CALL FlagError("The distributed vector has not been finished.",err,error,*999)
    IF(distributedVector%dataType/=MATRIX_VECTOR_L_TYPE) THEN
      localError="The distributed data type of "// &
        & TRIM(NumberToVString(distributedVector%dataType,"*",err,error))// &
        & " does not correspond to the logical data type of the given values."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(distributedVector%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissVector)
      CALL DistributedVector_CMISSVectorGet(distributedVector,cmissVector,err,error,*999)
      IF(index<=0.OR.index>cmissVector%dataSize) THEN
        localError="Index "//TRIM(NumberToVString(index,"*",err,error))// &
          & " is invalid. The index must be between 1 and "// &
          & TRIM(NumberToVString(cmissVector%dataSize,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      value=cmissVector%dataL(index)
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot set values for a logical PETSc distributed vector.",err,error,*999)          
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedVector_ValuesGetL1")
    RETURN
999 ERRORSEXITS("DistributedVector_ValuesGetL1",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_ValuesGetL1

  !
  !================================================================================================================================
  !

  !>Sets values in a distributed integer vector.
  SUBROUTINE DistributedVector_ValuesSetIntg(distributedVector,indices,values,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    INTEGER(INTG), INTENT(IN) :: indices(:) !<indices(i). The i'th index to be set
    INTEGER(INTG), INTENT(IN) :: values(:) !<values(i). The i'th value to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: i
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedVector_ValuesSetIntg",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(.NOT.distributedVector%vectorFinished) CALL FlagError("The distributed vector has not been finished.",err,error,*999)
    IF(SIZE(indices,1)/=SIZE(values,1)) THEN
      localError="The size of the indicies array of "//TRIM(NumberToVString(SIZE(indices,1),"*",err,error))// &
        & " does not conform to the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(distributedVector%dataType/=MATRIX_VECTOR_INTG_TYPE) THEN
      localError="The distributed data type of "// &
        & TRIM(NumberToVString(distributedVector%dataType,"*",err,error))// &
        & " does not correspond to the integer data type of the given values."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(distributedVector%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissVector)
      CALL DistributedVector_CMISSVectorGet(distributedVector,cmissVector,err,error,*999)
      DO i=1,SIZE(indices,1)
        !Allow all values set until dof mappings fixed. Ghost values that are set will not be propogated
        IF(indices(i)<=0.OR.indices(i)>cmissVector%dataSize) THEN
          localError="The index at position "//TRIM(NumberToVString(i,"*",err,error))//" of "// &
            & TRIM(NumberToVString(indices(i),"*",err,error))//" is invalid. The index must be between 1 and "// &
            & TRIM(NumberToVString(cmissVector%dataSize,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        cmissVector%dataIntg(indices(i))=values(i)
      ENDDO !i
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot set values for an integer PETSc distributed vector.",err,error,*999)          
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedVector_ValuesSetIntg")
    RETURN
999 ERRORSEXITS("DistributedVector_ValuesSetIntg",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_ValuesSetIntg

  !
  !================================================================================================================================
  !

  !>Sets one value in a distributed integer vector.
  SUBROUTINE DistributedVector_ValuesSetIntg1(distributedVector,index,value,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    INTEGER(INTG), INTENT(IN) :: index !<The index to be set
    INTEGER(INTG), INTENT(IN) :: value !<The value to be set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedVector_ValuesSetIntg1",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(.NOT.distributedVector%vectorFinished) CALL FlagError("The distributed vector has not been finished.",err,error,*999)
    IF(distributedVector%dataType/=MATRIX_VECTOR_INTG_TYPE) THEN
      localError="The distributed data type of "// &
        & TRIM(NumberToVString(distributedVector%dataType,"*",err,error))// &
        & " does not correspond to the integer data type of the given values."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(distributedVector%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissVector)
      CALL DistributedVector_CMISSVectorGet(distributedVector,cmissVector,err,error,*999)
      !Allow all values set until dof mappings fixed. Ghost values that are set will not be propogated
      IF(index<=0.OR.index>cmissVector%dataSize) THEN
        localError="Index "//TRIM(NumberToVString(index,"*",err,error))// &
          & " is invalid. The index must be between 1 and "// &
          & TRIM(NumberToVString(cmissVector%dataSize,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      cmissVector%dataIntg(index)=value
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot set values for an integer PETSc distributed vector.",err,error,*999)          
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedVector_ValuesSetIntg1")
    RETURN
999 ERRORSEXITS("DistributedVector_ValuesSetIntg1",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_ValuesSetIntg1

  !
  !================================================================================================================================
  !

  !>Sets values in a distributed single precision vector.
  SUBROUTINE DistributedVector_ValuesSetSP(distributedVector,indices,values,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    INTEGER(INTG), INTENT(IN) :: indices(:) !<indices(i). The i'th index to be set
    REAL(SP), INTENT(IN) :: values(:) !<values(i). The i'th value to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: i
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedVector_ValuesSetSP",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(.NOT.distributedVector%vectorFinished) CALL FlagError("The distributed vector has not been finished.",err,error,*999)
    IF(SIZE(indices,1)/=SIZE(values,1)) THEN
      localError="The size of the indicies array of "//TRIM(NumberToVString(SIZE(indices,1),"*",err,error))// &
        & " does not conform to the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(distributedVector%dataType/=MATRIX_VECTOR_SP_TYPE) THEN
      localError="The distributed data type of "// &
        & TRIM(NumberToVString(distributedVector%dataType,"*",err,error))// &
        & " does not correspond to the single precision real data type of the given values."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(distributedVector%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissVector)
      CALL DistributedVector_CMISSVectorGet(distributedVector,cmissVector,err,error,*999)
      DO i=1,SIZE(indices,1)
        !Allow all values set until dof mappings fixed. Ghost values that are set will not be propogated
        IF(indices(i)<=0.OR.indices(i)>cmissVector%dataSize) THEN
          localError="The index at position "//TRIM(NumberToVString(i,"*",err,error))//" of "// &
            & TRIM(NumberToVString(indices(i),"*",err,error))//" is invalid. The index must be between 1 and "// &
            & TRIM(NumberToVString(cmissVector%dataSize,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        cmissVector%dataSP(indices(i))=values(i)
      ENDDO !i
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot get values for a single precision PETSc distributed vector.",err,error,*999)          
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedVector_ValuesSetSP")
    RETURN
999 ERRORSEXITS("DistributedVector_ValuesSetSP",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_ValuesSetSP

  !
  !================================================================================================================================
  !

  !>Sets one value in a distributed single precision vector.
  SUBROUTINE DistributedVector_ValuesSetSP1(distributedVector,index,value,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    INTEGER(INTG), INTENT(IN) :: index !<The index to be set
    REAL(SP), INTENT(IN) :: value !<The value to be set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedVector_ValuesSetSP1",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(.NOT.distributedVector%vectorFinished) CALL FlagError("The distributed vector has not been finished.",err,error,*999)
    IF(distributedVector%dataType/=MATRIX_VECTOR_SP_TYPE) THEN
      localError="The distributed data type of "// &
        & TRIM(NumberToVString(distributedVector%dataType,"*",err,error))// &
        & " does not correspond to the single precision real data type of the given values."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(distributedVector%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissVector)
      CALL DistributedVector_CMISSVectorGet(distributedVector,cmissVector,err,error,*999)
      !Allow all values set until dof mappings fixed. Ghost values that are set will not be propogated
      IF(index<=0.OR.index>cmissVector%dataSize) THEN
        localError="Index "//TRIM(NumberToVString(index,"*",err,error))// &
          & " is invalid. The index must be between 1 and "// &
          & TRIM(NumberToVString(cmissVector%dataSize,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      cmissVector%dataSP(index)=value
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot set values for a single precision PETSc distributed vector.",err,error,*999)          
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedVector_ValuesSetSP1")
    RETURN
999 ERRORSEXITS("DistributedVector_ValuesSetSP1",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_ValuesSetSP1

  !
  !================================================================================================================================
  !

  !>Sets values in a distributed double precision vector.
  SUBROUTINE DistributedVector_ValuesSetDP(distributedVector,indices,values,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    INTEGER(INTG), INTENT(IN) :: indices(:) !<indices(i). The i'th index to be set
    REAL(DP), INTENT(IN) :: values(:) !<values(i). The i'th value to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: i
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector
    TYPE(DistributedVectorPETScType), POINTER :: petscVector
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedVector_ValuesSetDP",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(.NOT.distributedVector%vectorFinished) CALL FlagError("The distributed vector has not been finished.",err,error,*999)
    IF(SIZE(indices,1)/=SIZE(values,1)) THEN
      localError="The size of the indicies array of "//TRIM(NumberToVString(SIZE(indices,1),"*",err,error))// &
        & " does not conform to the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(distributedVector%dataType/=MATRIX_VECTOR_DP_TYPE) THEN
      localError="The distributed data type of "// &
        & TRIM(NumberToVString(distributedVector%dataType,"*",err,error))// &
        & " does not correspond to the double precision real data type of the given values."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(distributedVector%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissVector)
      CALL DistributedVector_CMISSVectorGet(distributedVector,cmissVector,err,error,*999)
      DO i=1,SIZE(indices,1)
        !Allow all values set until dof mappings fixed. Ghost values that are set will not be propogated
        IF(indices(i)<=0.OR.indices(i)>cmissVector%dataSize) THEN
          localError="The index at position "//TRIM(NumberToVString(i,"*",err,error))//" of "// &
            & TRIM(NumberToVString(indices(i),"*",err,error))//" is invalid. The index must be between 1 and "// &
            & TRIM(NumberToVString(cmissVector%dataSize,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        cmissVector%dataDP(indices(i))=values(i)
      ENDDO !i
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      NULLIFY(petscVector)
      CALL DistributedVector_PETScVectorGet(distributedVector,petscVector,err,error,*999)
      IF(petscVector%useOverrideVector) THEN
        CALL Petsc_VecSetValues(petscVector%overrideVector,SIZE(indices,1),petscVector%globalNumbers(indices), &
          & values,PETSC_INSERT_VALUES,err,error,*999)
      ELSE
        CALL Petsc_VecSetValues(petscVector%vector,SIZE(indices,1),petscVector%globalNumbers(indices), &
          & values,PETSC_INSERT_VALUES,err,error,*999)
      ENDIF
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
   
    EXITS("DistributedVector_ValuesSetDP")
    RETURN
999 ERRORSEXITS("DistributedVector_ValuesSetDP",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_ValuesSetDP

  !
  !================================================================================================================================
  !

  !>Sets one value in a distributed double precision vector.
  SUBROUTINE DistributedVector_ValuesSetDP1(distributedVector,index,value,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    INTEGER(INTG), INTENT(IN) :: index !<The index to be set
    REAL(DP), INTENT(IN) :: value !<The value to be set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: petscIndex(1)
    REAL(DP) :: petscValue(1)
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector
    TYPE(DistributedVectorPETScType), POINTER :: petscVector
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedVector_ValuesSetDP1",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(.NOT.distributedVector%vectorFinished) CALL FlagError("The distributed vector has not been finished.",err,error,*999)
    IF(distributedVector%dataType/=MATRIX_VECTOR_DP_TYPE) THEN
      localError="The distributed data type of "// &
        & TRIM(NumberToVString(distributedVector%dataType,"*",err,error))// &
        & " does not correspond to the double precision real data type of the given values."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(distributedVector%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissVector)
      CALL DistributedVector_CMISSVectorGet(distributedVector,cmissVector,err,error,*999)
      !Allow all values set until dof mappings fixed. Ghost values that are set will not be propogated
      IF(index<=0.OR.index>cmissVector%dataSize) THEN
        localError="Index "//TRIM(NumberToVString(index,"*",err,error))// &
          & " is invalid. The index must be between 1 and "// &
          & TRIM(NumberToVString(cmissVector%dataSize,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      cmissVector%dataDP(index)=value
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      NULLIFY(petscVector)
      CALL DistributedVector_PETScVectorGet(distributedVector,petscVector,err,error,*999)
      petscIndex(1)=petscVector%globalNumbers(index)
      petscValue(1)=value
      IF(petscVector%useOverrideVector) THEN
        CALL Petsc_VecSetValues(petscVector%overrideVector,1,petscIndex,petscValue,PETSC_INSERT_VALUES,err,error,*999)
      ELSE
        CALL Petsc_VecSetValues(petscVector%vector,1,petscIndex,petscValue,PETSC_INSERT_VALUES,err,error,*999)
      ENDIF
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedVector_ValuesSetDP1")
    RETURN
999 ERRORSEXITS("DistributedVector_ValuesSetDP1",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_ValuesSetDP1

  !
  !================================================================================================================================
  !

  !>Sets values in a distributed logical vector.
  SUBROUTINE DistributedVector_ValuesSetL(distributedVector,indices,values,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    INTEGER(INTG), INTENT(IN) :: indices(:) !<indices(i). The i'th index to be set
    LOGICAL, INTENT(IN) :: values(:) !<values(i). The i'th value to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: i
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedVector_ValuesSetL",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(.NOT.distributedVector%vectorFinished) CALL FlagError("The distributed vector has not been finished.",err,error,*999)
    IF(SIZE(indices,1)/=SIZE(values,1)) THEN
      localError="The size of the indicies array of "//TRIM(NumberToVString(SIZE(indices,1),"*",err,error))// &
        & " does not conform to the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(distributedVector%dataType/=MATRIX_VECTOR_L_TYPE) THEN
      localError="The distributed data type of "// &
        & TRIM(NumberToVString(distributedVector%dataType,"*",err,error))// &
        & " does not correspond to the logical data type of the given values."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(distributedVector%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissVector)
      CALL DistributedVector_CMISSVectorGet(distributedVector,cmissVector,err,error,*999)
      DO i=1,SIZE(indices,1)
        !Allow all values set until dof mappings fixed. Ghost values that are set will not be propogated
        IF(indices(i)<=0.OR.indices(i)>cmissVector%dataSize) THEN
          localError="The index at position "//TRIM(NumberToVString(i,"*",err,error))//" of "// &
            & TRIM(NumberToVString(indices(i),"*",err,error))//" is invalid. The index must be between 1 and "// &
            & TRIM(NumberToVString(cmissVector%dataSize,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        cmissVector%dataL(indices(i))=values(i)
      ENDDO !i
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot set values for a logical PETSc distributed vector.",err,error,*999)          
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedVector_ValuesSetL")
    RETURN
999 ERRORSEXITS("DistributedVector_ValuesSetL",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_ValuesSetL

  !
  !================================================================================================================================
  !

  !>Sets one value in a distributed logical vector.
  SUBROUTINE DistributedVector_ValuesSetL1(distributedVector,index,value,err,error,*)

    !Argument variables
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    INTEGER(INTG), INTENT(IN) :: index !<The index to be set
    LOGICAL, INTENT(IN) :: value !<The value to be set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DistributedVector_ValuesSetL1",err,error,*999)

    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)
    IF(.NOT.distributedVector%vectorFinished) CALL FlagError("The distributed vector has not been finished.",err,error,*999)
    IF(distributedVector%dataType/=MATRIX_VECTOR_L_TYPE) THEN
      localError="The distributed data type of "// &
        & TRIM(NumberToVString(distributedVector%dataType,"*",err,error))// &
        & " does not correspond to the logical data type of the given values."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(distributedVector%libraryType)
    CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
      NULLIFY(cmissVector)
      CALL DistributedVector_CMISSVectorGet(distributedVector,cmissVector,err,error,*999)
      !Allow all values set until dof mappings fixed. Ghost values that are set will not be propogated
      IF(index<=0.OR.index>cmissVector%dataSize) THEN
        localError="Index "//TRIM(NumberToVString(index,"*",err,error))// &
          & " is invalid. The index must be between 1 and "// &
          & TRIM(NumberToVString(cmissVector%dataSize,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      cmissVector%dataL(index)=value
    CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
      CALL FlagError("Cannot set values for a logical PETSc distributed vector.",err,error,*999)          
    CASE DEFAULT
      localError="The distributed vector library type of "// &
        & TRIM(NumberToVString(distributedVector%libraryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DistributedVector_ValuesSetL1")
    RETURN
999 ERRORSEXITS("DistributedVector_ValuesSetL1",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_ValuesSetL1

  !
  !================================================================================================================================
  !
  
END MODULE DistributedMatrixVector
