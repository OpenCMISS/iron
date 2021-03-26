!> \file
!> \author Chris Bradley
!> \brief This module contains all solver matrices access method routines.
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

!> This module contains all solver matrices access method routines.
MODULE SolverMatricesAccessRoutines
  
  USE BaseRoutines
  USE ISO_VARYING_STRING
  USE Kinds
  USE Strings
  USE Types
  
#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup SolverMatricesRoutines_SelectMatricesTypes SolverMatricesRoutines::SelectMatricesTypes
  !> \brief The types of selection available for the solver matrices
  !> \see SolverMatricesRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_MATRICES_ALL=1 !<Select all the solver matrices and vectors \see SolverMatricesRoutines_SelectMatricesTypes,SolverMatricesRoutines
!  redundant when introducing dynamic nonlinear equations
!  INTEGER(INTG), PARAMETER :: SOLVER_MATRICES_DYNAMIC_ONLY=2 !<Select only the dynamic solver matrices and vectors \see SolverMatricesRoutines_SelectMatricesTypes,SolverMatricesRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_MATRICES_LINEAR_ONLY=3 !<Select only the linear solver matrices and vectors \see SolverMatricesRoutines_SelectMatricesTypes,SolverMatricesRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_MATRICES_NONLINEAR_ONLY=4 !<Select only the nonlinear solver matrices and vectors \see SolverMatricesRoutines_SelectMatricesTypes,SolverMatricesRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_MATRICES_JACOBIAN_ONLY=5 !<Select only the Jacobian solver matrix \see SolverMatricesRoutines_SelectMatricesTypes,SolverMatricesRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_MATRICES_RESIDUAL_ONLY=6 !<Select only the residual solver vector \see SolverMatricesRoutines_SelectMatricesTypes,SolverMatricesRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_MATRICES_RHS_ONLY=7 !<Select only the RHS solver vector \see SolverMatricesRoutines_SelectMatricesTypes,SolverMatricesRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_MATRICES_RHS_RESIDUAL_ONLY=8 !<Select only the residual and RHS solver vectors \see SolverMatricesRoutines_SelectMatricesTypes,SolverMatricesRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_MATRICES_LINEAR_RESIDUAL_ONLY=9 !<Select only the linear solver matrices and vectors plus the residual \see SolverMatricesRoutines_SelectMatricesTypes,SolverMatricesRoutines
  !>@}

  !Module types

  !Module variables

  !Interfaces

  INTERFACE SolverMatrices_StorageTypesGet
    MODULE PROCEDURE SolverMatrices_StorageTypesGet0
    MODULE PROCEDURE SolverMatrices_StorageTypesGet1
  END INTERFACE SolverMatrices_StorageTypesGet

  INTERFACE SolverMatrices_SymmetryTypesGet
    MODULE PROCEDURE SolverMatrices_SymmetryTypesGet0
    MODULE PROCEDURE SolverMatrices_SymmetryTypesGet1
  END INTERFACE SolverMatrices_SymmetryTypesGet

  PUBLIC SOLVER_MATRICES_ALL,SOLVER_MATRICES_LINEAR_ONLY,SOLVER_MATRICES_NONLINEAR_ONLY, &
    & SOLVER_MATRICES_JACOBIAN_ONLY,SOLVER_MATRICES_RESIDUAL_ONLY,SOLVER_MATRICES_RHS_ONLY, & 
    & SOLVER_MATRICES_RHS_RESIDUAL_ONLY,SOLVER_MATRICES_LINEAR_RESIDUAL_ONLY !,SOLVER_MATRICES_DYNAMIC_ONLY

  PUBLIC SolverMatrices_AssertIsFinished,SolverMatrices_AssertNotFinished
  
  PUBLIC SolverMatrices_LibraryTypeGet

  PUBLIC SolverMatrices_NumberOfSolverMatricesGet

  PUBLIC SolverMatrices_NumberOfRowsGet

  PUBLIC SolverMatrices_NumberOfGlobalRowsGet

  PUBLIC SolverMatrices_RHSDistributedVectorGet

  PUBLIC SolverMatrices_ResidualDistributedVectorGet

  PUBLIC SolverMatrices_SolverEquationsGet

  PUBLIC SolverMatrices_SolverMappingGet

  PUBLIC SolverMatrices_SolverMatrixGet

  PUBLIC SolverMatrices_StorageTypesGet

  PUBLIC SolverMatrices_SymmetryTypesGet

  PUBLIC SolverMatrices_UpdateResidualGet

  PUBLIC SolverMatrices_UpdateRHSGet

  PUBLIC SolverMatrix_MatrixNumberGet

  PUBLIC SolverMatrix_NumberOfColumnsGet

  PUBLIC SolverMatrix_SolverDistributedMatrixGet

  PUBLIC SolverMatrix_SolverDistributedVectorGet

  PUBLIC SolverMatrix_SolverMatricesGet

  PUBLIC SolverMatrix_StorageTypeGet

  PUBLIC SolverMatrix_SymmetryTypeGet

  PUBLIC SolverMatrix_UpdateMatrixGet
  
CONTAINS

  !
  !=================================================================================================================================
  !

  !>Assert that a solver matrices has been finished
  SUBROUTINE SolverMatrices_AssertIsFinished(solverMatrices,err,error,*)

    !Argument Variables
    TYPE(SolverMatricesType), POINTER, INTENT(IN) :: solverMatrices !<The solver matrices to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMatrices_AssertIsFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solverMatrices)) CALL FlagError("Solver matrices is not associated.",err,error,*999)
#endif    

    IF(.NOT.solverMatrices%solverMatricesFinished) CALL FlagError("Solver matrices has not been finished.",err,error,*999)
    
    EXITS("SolverMatrices_AssertIsFinished")
    RETURN
999 ERRORSEXITS("SolverMatrices_AssertIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMatrices_AssertIsFinished

  !
  !=================================================================================================================================
  !

  !>Assert that a solver matrices has not been finished
  SUBROUTINE SolverMatrices_AssertNotFinished(solverMatrices,err,error,*)

    !Argument Variables
    TYPE(SolverMatricesType), POINTER, INTENT(IN) :: solverMatrices !<The solver matrices to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMatrices_AssertNotFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solverMatrices)) CALL FlagError("Solver matrices is not associated.",err,error,*999)
#endif    

    IF(solverMatrices%solverMatricesFinished) CALL FlagError("Solver matrices has already been finished.",err,error,*999)
    
    EXITS("SolverMatrices_AssertNotFinished")
    RETURN
999 ERRORSEXITS("SolverMatrices_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMatrices_AssertNotFinished
  
  !
  !================================================================================================================================
  !
  
  !>Gets the library type for the solver matrices (and vectors)
  SUBROUTINE SolverMatrices_LibraryTypeGet(solverMatrices,libraryType,err,error,*)

    !Argument variables
    TYPE(SolverMatricesType), POINTER :: solverMatrices !<A pointer to the solver matrices.
    INTEGER(INTG), INTENT(OUT) :: libraryType !<On return, the library type of the specified solver matrices \see SolverRoutines_SolverLibraries
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverMatrices_LibraryTypeGet",err,error,*999)

    CALL SolverMatrices_AssertIsFinished(solverMatrices,err,error,*999)
    
    libraryType=solverMatrices%solverLibraryType
    
    EXITS("SolverMatrices_LibraryTypeGet")
    RETURN
999 ERRORSEXITS("SolverMatrices_LibraryTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMatrices_LibraryTypeGet
          
  !
  !================================================================================================================================
  !

  !>Get the number of solver matrices for the solver matrices
  SUBROUTINE SolverMatrices_NumberOfSolverMatricesGet(solverMatrices,numberOfMatrices,err,error,*)

    !Argument variables
    TYPE(SolverMatricesType), POINTER, INTENT(IN) :: solverMatrices !<The solver matrices to get the number of matrices for
    INTEGER(INTG), INTENT(OUT) :: numberOfMatrices !<On return, the number of matrices for the solver matrices
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables

    ENTERS("SolverMatrices_NumberOfSolverMatricesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solverMatrices)) CALL FlagError("Solver matrices are not associated.",err,error,*999)
#endif    
     
    numberOfMatrices=solverMatrices%numberOfMatrices

    EXITS("SolverMatrices_NumberOfSolverMatricesGet")
    RETURN
999 ERRORSEXITS("SolverMatrices_NumberOfSolverMatricesGet",err,error)
    RETURN 1

  END SUBROUTINE SolverMatrices_NumberOfSolverMatricesGet

  !
  !================================================================================================================================
  !

  !>Get the number of rows for the solver matrices
  SUBROUTINE SolverMatrices_NumberOfRowsGet(solverMatrices,numberOfRows,err,error,*)

    !Argument variables
    TYPE(SolverMatricesType), POINTER, INTENT(IN) :: solverMatrices !<The solver matrices to get the number of rows for
    INTEGER(INTG), INTENT(OUT) :: numberOfRows !<On return, the number of rows for the solver matrices
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables

    ENTERS("SolverMatrices_NumberOfRowsGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solverMatrices)) CALL FlagError("Solver matrices are not associated.",err,error,*999)
#endif    
     
    numberOfRows=solverMatrices%numberOfRows

    EXITS("SolverMatrices_NumberOfRowsGet")
    RETURN
999 ERRORSEXITS("SolverMatrices_NumberOfRowsGet",err,error)
    RETURN 1

  END SUBROUTINE SolverMatrices_NumberOfRowsGet

  !
  !================================================================================================================================
  !

  !>Get the number of global rows for the solver matrices
  SUBROUTINE SolverMatrices_NumberOfGlobalRowsGet(solverMatrices,numberOfGlobalRows,err,error,*)

    !Argument variables
    TYPE(SolverMatricesType), POINTER, INTENT(IN) :: solverMatrices !<The solver matrices to get the number of global rows for
    INTEGER(INTG), INTENT(OUT) :: numberOfGlobalRows !<On return, the number of global rows for the solver matrices
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables

    ENTERS("SolverMatrices_NumberOfGlobalRowsGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solverMatrices)) CALL FlagError("Solver matrices are not associated.",err,error,*999)
#endif    
     
    numberOfGlobalRows=solverMatrices%numberOfGlobalRows

    EXITS("SolverMatrices_NumberOfGlobalRowsGet")
    RETURN
999 ERRORSEXITS("SolverMatrices_NumberOfGlobalRowsGet",err,error)
    RETURN 1

  END SUBROUTINE SolverMatrices_NumberOfGlobalRowsGet

  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to specified solver RHS distributed vector for solver matrices.
  SUBROUTINE SolverMatrices_RHSDistributedVectorGet(solverMatrices,rhsDistributedVector,err,error,*)

    !Argument variables
    TYPE(SolverMatricesType), POINTER :: solverMatrices !<A pointer to the solver matrices to get the solver matrix for
    TYPE(DistributedVectorType), POINTER :: rhsDistributedVector !<On exit, a pointer to the solver RHS distributed vector. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMatrices_RHSDistributedVectorGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(rhsDistributedVector)) CALL FlagError("RHS distributed vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverMatrices)) CALL FlagError("Solver matrices is not associated.",err,error,*999)
#endif    
    
    rhsDistributedVector=>solverMatrices%rhsVector

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(rhsDistributedVector)) &
      & CALL FlagError("The RHS distributed vector is not associated for the solver matrices.",err,error,*999)
#endif    
      
    EXITS("SolverMatrices_RHSDistributedVectorGet")
    RETURN
999 NULLIFY(rhsDistributedVector)
998 ERRORSEXITS("SolverMatrices_RHSDistributedVectorGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMatrices_RHSDistributedVectorGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to specified solver residual distributed vector for solver matrices.
  SUBROUTINE SolverMatrices_ResidualDistributedVectorGet(solverMatrices,residualDistributedVector,err,error,*)

    !Argument variables
    TYPE(SolverMatricesType), POINTER :: solverMatrices !<A pointer to the solver matrices to get the residual distributed vector for
    TYPE(DistributedVectorType), POINTER :: residualDistributedVector !<On exit, a pointer to the solver residual distributed vector. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMatrices_ResidualDistributedVectorGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(residualDistributedVector)) CALL FlagError("Residual distributed vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverMatrices)) CALL FlagError("Solver matrices is not associated.",err,error,*999)
#endif    
    
    residualDistributedVector=>solverMatrices%residual

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(residualDistributedVector)) &
      & CALL FlagError("The residual distributed vector is not associated for the solver matrices.",err,error,*999)
#endif    
      
    EXITS("SolverMatrices_ResidualDistributedVectorGet")
    RETURN
999 NULLIFY(residualDistributedVector)
998 ERRORSEXITS("SolverMatrices_ResidualDistributedVectorGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMatrices_ResidualDistributedVectorGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the solver equations for solver matrices.
  SUBROUTINE SolverMatrices_SolverEquationsGet(solverMatrices,solverEquations,err,error,*)

    !Argument variables
    TYPE(SolverMatricesType), POINTER :: solverMatrices !<A pointer to the solver matrices to get the solver equations for
    TYPE(SolverEquationsType), POINTER :: solverEquations !<On exit, a pointer to the solver equations. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMatrices_SolverEquationsGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverMatrices)) CALL FlagError("Solver matrices is not associated.",err,error,*999)
#endif    

    solverEquations=>solverMatrices%solverEquations

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(solverEquations)) CALL FlagError("Solver matrices solver equations is not associated.",err,error,*999)
#endif    
      
    EXITS("SolverMatrices_SolverEquationsGet")
    RETURN
999 NULLIFY(solverEquations)
998 ERRORSEXITS("SolverMatrices_SolverEquationsGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMatrices_SolverEquationsGet

  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the solver mapping for solver matrices.
  SUBROUTINE SolverMatrices_SolverMappingGet(solverMatrices,solverMapping,err,error,*)

    !Argument variables
    TYPE(SolverMatricesType), POINTER :: solverMatrices !<A pointer to the solver matrices to get the solver mappings for
    TYPE(SolverMappingType), POINTER :: solverMapping !<On exit, a pointer to the solver mapping. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMatrices_SolverMappingGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(solverMapping)) CALL FlagError("Solver mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverMatrices)) CALL FlagError("Solver matrices is not associated.",err,error,*999)
#endif    

    solverMapping=>solverMatrices%solverMapping

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(solverMapping)) CALL FlagError("Solver matrices solver mapping is not associated.",err,error,*999)
#endif    
      
    EXITS("SolverMatrices_SolverMappingGet")
    RETURN
999 NULLIFY(solverMapping)
998 ERRORSEXITS("SolverMatrices_SolverMappingGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMatrices_SolverMappingGet

  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to specified solver matrix for solver matrices.
  SUBROUTINE SolverMatrices_SolverMatrixGet(solverMatrices,matrixIdx,solverMatrix,err,error,*)

    !Argument variables
    TYPE(SolverMatricesType), POINTER :: solverMatrices !<A pointer to the solver matrices to get the solver matrix for
    INTEGER(INTG), INTENT(IN) :: matrixIdx !<The matrix index of the solver matrix to get.
    TYPE(SolverMatrixType), POINTER :: solverMatrix !<On exit, a pointer to the solver matrix. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMatrices_SolverMatrixGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(solverMatrix)) CALL FlagError("Solver matrix is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverMatrices)) CALL FlagError("Solver matrices is not associated.",err,error,*999)
    IF(matrixIdx<1.OR.matrixIdx>solverMatrices%numberOfMatrices) THEN
      localError="The specified solver matrix index of "// &
        & TRIM(NumberToVString(matrixIdx,"*",err,error))//" is invalid. The matrix index needs to be >= 1 and <= "// &
        & TRIM(NumberToVString(solverMatrices%numberOfMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(solverMatrices%matrices)) CALL FlagError("Solver matrices matrices is not allocated.",err,error,*999)
#endif    
    
    solverMatrix=>solverMatrices%matrices(matrixIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(solverMatrix)) THEN
      localError="The solver matrix for index "//TRIM(NumberToVString(matrixIdx,"*",err,error))//" is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    EXITS("SolverMatrices_SolverMatrixGet")
    RETURN
999 NULLIFY(solverMatrix)
998 ERRORSEXITS("SolverMatrices_SolverMatrixGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMatrices_SolverMatrixGet
  
  !
  !================================================================================================================================
  !
  
  !>Gets the storage type (sparsity) of the solver matrices
  SUBROUTINE SolverMatrices_StorageTypesGet0(solverMatrices,storageType,err,error,*)

    !Argument variables
    TYPE(SolverMatricesType), POINTER :: solverMatrices !<A pointer to the solver matrices
    INTEGER(INTG), INTENT(OUT) :: storageType !<On return, the storage type for the solver matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx,storageTypes(1)
    TYPE(SolverMatrixType), POINTER :: solverMatrix
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("SolverMatrices_StorageTypesGet0",err,error,*999)

    CALL SolverMatrices_StorageTypesGet1(solverMatrices,storageTypes,err,error,*999)
    storageType=storageTypes(1)
    
    EXITS("SolverMatrices_StorageTypesGet0")
    RETURN
999 ERRORSEXITS("SolverMatrices_StorageTypesGet0",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMatrices_StorageTypesGet0

  !
  !================================================================================================================================
  !
  
  !>Gets the storage type (sparsity) of the solver matrices
  SUBROUTINE SolverMatrices_StorageTypesGet1(solverMatrices,storageTypes,err,error,*)

    !Argument variables
    TYPE(SolverMatricesType), POINTER :: solverMatrices !<A pointer to the solver matrices
    INTEGER(INTG), INTENT(OUT) :: storageTypes(:) !<storageTypes(matrixIdx). On return, the storage type for the matrixIdx'th solver matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(SolverMatrixType), POINTER :: solverMatrix
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("SolverMatrices_StorageTypesGet1",err,error,*999)

    CALL SolverMatrices_AssertIsFinished(solverMatrices,err,error,*999)
#ifdef WITH_PRECHECKS    
    IF(SIZE(storageTypes,1)<solverMatrices%numberOfMatrices) THEN
      localError="The size of storage types array is too small. The supplied size is "// &
        & TRIM(NumberToVString(SIZE(storageTypes,1),"*",err,error))//" and it needs to be >= "// &
        & TRIM(NumberToVString(solverMatrices%numberOfMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    DO matrixIdx=1,solverMatrices%numberOfMatrices
      NULLIFY(solverMatrix)
      CALL SolverMatrices_SolverMatrixGet(solverMatrices,matrixIdx,solverMatrix,err,error,*999)
      storageTypes(matrixIdx)=solverMatrix%storageType
    ENDDO !matrixIdx
    
    EXITS("SolverMatrices_StorageTypesGet1")
    RETURN
999 ERRORSEXITS("SolverMatrices_StorageTypesGet1",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMatrices_StorageTypesGet1

  !
  !================================================================================================================================
  !
  
  !>Gets the symmetry type of the solver matrices
  SUBROUTINE SolverMatrices_SymmetryTypesGet0(solverMatrices,symmetryType,err,error,*)

    !Argument variables
    TYPE(SolverMatricesType), POINTER :: solverMatrices !<A pointer to the solver matrices
    INTEGER(INTG), INTENT(OUT) :: symmetryType !<On return, the symmtry type for the solver matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx,symmetryTypes(1)
    TYPE(SolverMatrixType), POINTER :: solverMatrix
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("SolverMatrices_SymmetryTypesGet0",err,error,*999)

    CALL SolverMatrices_SymmetryTypesGet1(solverMatrices,symmetryTypes,err,error,*999)
    symmetryType=symmetryTypes(1)

    EXITS("SolverMatrices_SymmetryTypesGet0")
    RETURN
999 ERRORSEXITS("SolverMatrices_SymmetryTypesGet0",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMatrices_SymmetryTypesGet0

  !
  !================================================================================================================================
  !
  
  !>Gets the symmetry types of the solver matrices
  SUBROUTINE SolverMatrices_SymmetryTypesGet1(solverMatrices,symmetryTypes,err,error,*)

    !Argument variables
    TYPE(SolverMatricesType), POINTER :: solverMatrices !<A pointer to the solver matrices
    INTEGER(INTG), INTENT(OUT) :: symmetryTypes(:) !<symmetryTypes(matrixIdx). On return, the symmtry type for the matrixIdx'th solver matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(SolverMatrixType), POINTER :: solverMatrix
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("SolverMatrices_SymmetryTypesGet1",err,error,*999)

    CALL SolverMatrices_AssertIsFinished(solverMatrices,err,error,*999)
#ifdef WITH_PRECHECKS    
    IF(SIZE(symmetryTypes,1)<solverMatrices%numberOfMatrices) THEN
      localError="The size of symmetry types is too small. The supplied size is "// &
        & TRIM(NumberToVString(SIZE(symmetryTypes,1),"*",err,error))//" and it needs to be >= "// &
        & TRIM(NumberToVString(solverMatrices%numberOfMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    DO matrixIdx=1,solverMatrices%numberOfMatrices
      NULLIFY(solverMatrix)
      CALL SolverMatrices_SolverMatrixGet(solverMatrices,matrixIdx,solverMatrix,err,error,*999)
      symmetryTypes(matrixIdx)=solverMatrix%symmetryType
    ENDDO !matrixIdx
    
    EXITS("SolverMatrices_SymmetryTypesGet1")
    RETURN
999 ERRORSEXITS("SolverMatrices_SymmetryTypesGet1",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMatrices_SymmetryTypesGet1

  !
  !================================================================================================================================
  !
  
  !>Gets the update residual flag for solver matrices
  SUBROUTINE SolverMatrices_UpdateResidualGet(solverMatrices,updateResidual,err,error,*)

    !Argument variables
    TYPE(SolverMatricesType), POINTER :: solverMatrices !<A pointer to the solver matrices to get the update residual flag for
    LOGICAL, INTENT(OUT) :: updateResidual !<On return, the update flag for the solver residual vector.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("SolverMatrices_UpdateResidualGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solverMatrices)) CALL FlagError("Solver matrices is not associated.",err,error,*999)
#endif    

    updateResidual=solverMatrices%updateResidual
    
    EXITS("SolverMatrices_UpdateResidualGet")
    RETURN
999 ERRORSEXITS("SolverMatrices_UpdateResidualGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMatrices_UpdateResidualGet

  !
  !================================================================================================================================
  !
  
  !>Gets the update RHS flag for solver matrices
  SUBROUTINE SolverMatrices_UpdateRHSGet(solverMatrices,updateRHS,err,error,*)

    !Argument variables
    TYPE(SolverMatricesType), POINTER :: solverMatrices !<A pointer to the solver matrices to get the update RHS flag for
    LOGICAL, INTENT(OUT) :: updateRHS !<On return, the update flag for the solver RHS vector.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("SolverMatrices_UpdateRHSGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solverMatrices)) CALL FlagError("Solver matrices is not associated.",err,error,*999)
#endif    

    updateRHS=solverMatrices%updateRHSVector
    
    EXITS("SolverMatrices_UpdateRHSGet")
    RETURN
999 ERRORSEXITS("SolverMatrices_UpdateRHSGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMatrices_UpdateRHSGet

  !
  !================================================================================================================================
  !
  
  !>Gets the number of a solver matrix
  SUBROUTINE SolverMatrix_MatrixNumberGet(solverMatrix,matrixNumber,err,error,*)

    !Argument variables
    TYPE(SolverMatrixType), POINTER :: solverMatrix !<A pointer to the solver matrix to get the number for
    INTEGER(INTG), INTENT(OUT) :: matrixNumber !<On return, the matrix number of solver matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("SolverMatrix_MatrixNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solverMatrix)) CALL FlagError("Solver matrix is not associated.",err,error,*999)
#endif    

    matrixNumber=solverMatrix%matrixNumber
    
    EXITS("SolverMatrix_MatrixNumberGet")
    RETURN
999 ERRORSEXITS("SolverMatrix_MatrixNumberGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMatrix_MatrixNumberGet

  !
  !================================================================================================================================
  !
  
  !>Gets the number of columns in a solver matrix
  SUBROUTINE SolverMatrix_NumberOfColumnsGet(solverMatrix,numberOfColumns,err,error,*)

    !Argument variables
    TYPE(SolverMatrixType), POINTER :: solverMatrix !<A pointer to the solver matrix to get the number of columns for
    INTEGER(INTG), INTENT(OUT) :: numberOfColumns !<On return, the number of columns in the solver matrix.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
   
    ENTERS("SolverMatrix_NumberOfColumnsGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solverMatrix)) CALL FlagError("Solver matrix is not associated.",err,error,*999)
#endif    

    numberOfColumns=solverMatrix%numberOfColumns
    
    EXITS("SolverMatrix_NumberOfColumnsGet")
    RETURN
999 ERRORSEXITS("SolverMatrix_NumberOfColumnsGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMatrix_NumberOfColumnsGet

  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to solver distributed matrix for a solver matrix.
  SUBROUTINE SolverMatrix_SolverDistributedMatrixGet(solverMatrix,distributedMatrix,err,error,*)

    !Argument variables
    TYPE(SolverMatrixType), POINTER :: solverMatrix !<A pointer to the solver matrix to get the solver distributed matrix for
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<On exit, a pointer to the solver distributed matrix for the specified solver matrix. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMatrix_SolverDistributedMatrixGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverMatrix)) CALL FlagError("Solver matrix is not associated.",err,error,*999)
#endif    
    
    distributedMatrix=>solverMatrix%matrix

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(distributedMatrix)) &
      & CALL FlagError("The solver distributed matrix is not associated for the solver matrix.",err,error,*999)
#endif    
      
    EXITS("SolverMatrix_SolverDistributedMatrixGet")
    RETURN
999 NULLIFY(distributedMatrix)
998 ERRORSEXITS("SolverMatrix_SolverDistributedMatrixGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMatrix_SolverDistributedMatrixGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to solver distributed vector for a solver matrix.
  SUBROUTINE SolverMatrix_SolverDistributedVectorGet(solverMatrix,distributedVector,err,error,*)

    !Argument variables
    TYPE(SolverMatrixType), POINTER :: solverMatrix !<A pointer to the solver matrix to get the solver distributed vector for
    TYPE(DistributedVectorType), POINTER :: distributedVector !<On exit, a pointer to the solver distributed vector for the specified solver matrix. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMatrix_SolverDistributedVectorGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(distributedVector)) CALL FlagError("Solver distributed vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverMatrix)) CALL FlagError("Solver matrix is not associated.",err,error,*999)
#endif    
    
    distributedVector=>solverMatrix%solverVector

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(distributedVector)) &
      & CALL FlagError("The solver distributed vector is not associated for the solver matrix.",err,error,*999)
#endif    
      
    EXITS("SolverMatrix_SolverDistributedVectorGet")
    RETURN
999 NULLIFY(distributedVector)
998 ERRORSEXITS("SolverMatrix_SolverDistributedVectorGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMatrix_SolverDistributedVectorGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to solver matrices for a solver matrix.
  SUBROUTINE SolverMatrix_SolverMatricesGet(solverMatrix,solverMatrices,err,error,*)

    !Argument variables
    TYPE(SolverMatrixType), POINTER :: solverMatrix !<A pointer to the solver matrix to get the solver matrices for
    TYPE(SolverMatricesType), POINTER :: solverMatrices !<On exit, a pointer to the solver matrices for the specified solver matrix. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMatrix_SolverMatricesGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(solverMatrices)) CALL FlagError("Solver matrices is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverMatrix)) CALL FlagError("Solver matrix is not associated.",err,error,*999)
#endif    
    
    solverMatrices=>solverMatrix%solverMatrices

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(solverMatrices)) &
      & CALL FlagError("The solver matrices is not associated for the solver matrix.",err,error,*999)
#endif    
      
    EXITS("SolverMatrix_SolverMatricesGet")
    RETURN
999 NULLIFY(solverMatrices)
998 ERRORSEXITS("SolverMatrix_SolverMatricesGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMatrix_SolverMatricesGet
  
  !
  !================================================================================================================================
  !
  
  !>Gets the storage type for a solver matrix
  SUBROUTINE SolverMatrix_StorageTypeGet(solverMatrix,storageType,err,error,*)

    !Argument variables
    TYPE(SolverMatrixType), POINTER :: solverMatrix !<A pointer to the solver matrix to get the storage type for
    INTEGER(INTG), INTENT(OUT) :: storageType !<On return, the storage type the solver matrix.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("SolverMatrix_StorageTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solverMatrix)) CALL FlagError("Solver matrix is not associated.",err,error,*999)
#endif    

    storageType=solverMatrix%storageType
    
    EXITS("SolverMatrix_StorageTypeGet")
    RETURN
999 ERRORSEXITS("SolverMatrix_StorageTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMatrix_StorageTypeGet

  !
  !================================================================================================================================
  !
  
  !>Gets the symmetry type for a solver matrix
  SUBROUTINE SolverMatrix_SymmetryTypeGet(solverMatrix,symmetryType,err,error,*)

    !Argument variables
    TYPE(SolverMatrixType), POINTER :: solverMatrix !<A pointer to the solver matrix to get the symmetry type for
    INTEGER(INTG), INTENT(OUT) :: symmetryType !<On return, the symmetry type the solver matrix.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("SolverMatrix_SymmetryTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solverMatrix)) CALL FlagError("Solver matrix is not associated.",err,error,*999)
#endif    

    symmetryType=solverMatrix%symmetryType
    
    EXITS("SolverMatrix_SymmetryTypeGet")
    RETURN
999 ERRORSEXITS("SolverMatrix_SymmetryTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMatrix_SymmetryTypeGet

  !
  !================================================================================================================================
  !
  
  !>Gets the update flag for a solver matrix
  SUBROUTINE SolverMatrix_UpdateMatrixGet(solverMatrix,updateMatrix,err,error,*)

    !Argument variables
    TYPE(SolverMatrixType), POINTER :: solverMatrix !<A pointer to the solver matrix to get the update flag for
    LOGICAL, INTENT(OUT) :: updateMatrix !<On return, the update flag the solver matrix.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("SolverMatrix_UpdateMatrixGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solverMatrix)) CALL FlagError("Solver matrix is not associated.",err,error,*999)
#endif    

    updateMatrix=solverMatrix%updateMatrix
    
    EXITS("SolverMatrix_UpdateMatrixGet")
    RETURN
999 ERRORSEXITS("SolverMatrix_UpdateMatrixGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMatrix_UpdateMatrixGet

  !
  !================================================================================================================================
  !

END MODULE SolverMatricesAccessRoutines
