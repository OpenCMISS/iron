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

  PUBLIC SOLVER_MATRICES_ALL,SOLVER_MATRICES_LINEAR_ONLY,SOLVER_MATRICES_NONLINEAR_ONLY, &
    & SOLVER_MATRICES_JACOBIAN_ONLY,SOLVER_MATRICES_RESIDUAL_ONLY,SOLVER_MATRICES_RHS_ONLY, & 
    & SOLVER_MATRICES_RHS_RESIDUAL_ONLY,SOLVER_MATRICES_LINEAR_RESIDUAL_ONLY !,SOLVER_MATRICES_DYNAMIC_ONLY

  PUBLIC SolverMatrices_AssertIsFinished,SolverMatrices_AssertNotFinished
  
  PUBLIC SolverMatrices_LibraryTypeGet

  PUBLIC SolverMatrices_NumberOfMatricesGet

  PUBLIC SolverMatrices_RHSDistributedVectorGet

  PUBLIC SolverMatrices_ResidualDistributedVectorGet

  PUBLIC SolverMatrices_SolverMappingGet

  PUBLIC SolverMatrices_SolverMatrixGet

  PUBLIC SolverMatrices_StorageTypeGet

  PUBLIC SolverMatrices_SymmetryTypeGet

  PUBLIC SolverMatrix_SolverDistributedMatrixGet

  PUBLIC SolverMatrix_SolverDistributedVectorGet
  
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

    IF(.NOT.ASSOCIATED(solverMatrices)) CALL FlagError("Solver matrices is not associated.",err,error,*999)

    IF(.NOT.solverMatrices%solverMatricesFinished) &
      & CALL FlagError("Solver matrices has not been finished."
    
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

    IF(.NOT.ASSOCIATED(solverMatrices)) CALL FlagError("Solver matrices is not associated.",err,error,*999)

    IF(solverMatrices%solverMatricesFinished) &
      & CALL FlagError("Solver matrices has already been finished.",err,error,*999)
    
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
  SUBROUTINE SolverMatrices_NumberOfMatricesGet(solverMatrices,numberOfMatrices,err,error,*)

    !Argument variables
    TYPE(SolverMatricesType), POINTER, INTENT(IN) :: solverMatrices !<The solver matrices to get the number of matrices for
    INTEGER(INTG), INTENT(OUT) :: numberOfMatrices !<On return, the number of matrices for the solver matrices
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables

    ENTERS("SolverMatrices_NumberOfMatricesGet",err,error,*999)

    IF(.NOT.ASSOCIATED(solverMatrices)) CALL FlagError("Solver matrices are not associated.",err,error,*999)
     
    numberOfMatrices=solverMatrices%numberOfMatrices

    EXITS("SolverMatrices_NumberOfMatricesGet")
    RETURN
999 ERRORSEXITS("SolverMatrices_NumberOfMatricesGet",err,error)
    RETURN 1

  END SUBROUTINE SolverMatrices_NumberOfMatricesGet

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

    IF(ASSOCIATED(rhsDistributedVector)) CALL FlagError("RHS distributed vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverMatrices)) CALL FlagError("Solver matrices is not associated.",err,error,*999)
    
    rhsDistributedVector=>solverMatrices%rhsVector
    IF(.NOT.ASSOCIATED(rhsDistributedVector)) &
      & CALL FlagError("The RHS distributed vector is not associated for the solver matrices.",err,error,*999)
      
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

    IF(ASSOCIATED(residualDistributedVector)) CALL FlagError("Residual distributed vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverMatrices)) CALL FlagError("Solver matrices is not associated.",err,error,*999)
    
    residualDistributedVector=>solverMatrices%residual
    IF(.NOT.ASSOCIATED(residualDistributedVector)) &
      & CALL FlagError("The residual distributed vector is not associated for the solver matrices.",err,error,*999)
      
    EXITS("SolverMatrices_ResidualDistributedVectorGet")
    RETURN
999 NULLIFY(residualDistributedVector)
998 ERRORSEXITS("SolverMatrices_ResidualDistributedVectorGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMatrices_ResidualDistributedVectorGet
  
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

    IF(ASSOCIATED(solverMapping)) CALL FlagError("Solver mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverMatrices)) CALL FlagError("Solver matrices is not associated.",err,error,*999)

    solverMapping=>solverMatrices%solverMapping
    IF(.NOT.ASSOCIATED(solverMapping)) CALL FlagError("Solver matrices solver mapping is not associated.",err,error,*999)
      
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
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("SolverMatrices_SolverMatrixGet",err,error,*998)

    IF(ASSOCIATED(solverMatrix)) CALL FlagError("Solver matrix is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverMatrices)) CALL FlagError("Solver matrices is not associated.",err,error,*999)
    IF(matrixIdx<1.OR.matrixIdx>solverMatrices%numberOfMatrices) THEN
      localError="The specified solver matrix index of "// &
        & TRIM(NumberToVString(matrixIdx,"*",err,error))//" is invalid. The matrix index needs to be >= 1 and <= "// &
        & TRIM(NumberToVString(solverMatrices%numberOfMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(solverMatrices%matrices)) CALL FlagError("Solver matrices matrices is not allocated.",err,error,*999)
    
    solverMatrix=>solverMatrices%matrices(matrixIdx)%ptr
    IF(.NOT.ASSOCIATED(solverMatrix)) THEN
      localError="The solver matrix for index "//TRIM(NumberToVString(matrixIdx,"*",err,error))//" is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
      
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
  SUBROUTINE SolverMatrices_StorageTypeGet(solverMatrices,storageType,err,error,*)

    !Argument variables
    TYPE(SolverMatricesType), POINTER :: solverMatrices !<A pointer to the solver matrices
    INTEGER(INTG), INTENT(OUT) :: storageType(:) !<storageType(matrixIdx). On return, the storage type for the matrixIdx'th solver matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(SolverMatrixType), POINTER :: solverMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("SolverMatrices_StorageTypeGet",err,error,*999)

    CALL SolverMatrices_AssertIsFinished(solverMatrices,err,error,*999)
    IF(SIZE(storageType,1)<solverMatrices%numberOfMatrices) THEN
      localError="The size of storage type array is too small. The supplied size is "// &
        & TRIM(NumberToVString(SIZE(storageType,1),"*",err,error))//" and it needs to be >= "// &
        & TRIM(NumberToVString(solverMatrices%numberOfMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    DO matrixIdx=1,solverMatrices%numberOfMatrices
      NULLIFY(solverMatrix)
      CALL SolverMatrices_SolverMatrixGet(solverMatrices,matrixIdx,solverMatrix,err,error,*999)
      storageType(matrixIdx)=solverMatrix%storageType
    ENDDO !matrixIdx
    
    EXITS("SolverMatrices_StorageTypeGet")
    RETURN
999 ERRORSEXITS("SolverMatrices_StorageTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMatrices_StorageTypeGet

  !
  !================================================================================================================================
  !
  
  !>Gets the symmetry type of the solver matrices
  SUBROUTINE SolverMatrices_SymmetryTypeGet(solverMatrices,symmetryTypes,err,error,*)

    !Argument variables
    TYPE(SolverMatricesType), POINTER :: solverMatrices !<A pointer to the solver matrices
    INTEGER(INTG), INTENT(OUT) :: symmetryTypes(:) !<symmetryTypes(matrixIdx). On return, the symmtry type for the matrixIdx'th solver matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(SolverMatrixType), POINTER :: solverMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("SolverMatrices_SymmetryTypeGet",err,error,*999)

    CALL SolverMatrices_AssertIsFinished(solverMatrices,err,error,*999)
    IF(SIZE(symmetryTypes,1)<solverMatrices%numberOfMatrices) THEN
      localError="The size of symmetry types is too small. The supplied size is "// &
        & TRIM(NumberToVString(SIZE(symmetryTypes,1),"*",err,error))//" and it needs to be >= "// &
        & TRIM(NumberToVString(solverMatrices%numberOfMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    DO matrixIdx=1,solverMatrices%numberOfMatrices
      NULLIFY(solverMatrix)
      CALL SolverMatrices_SolverMatrixGet(solverMatrices,matrixIdx,solverMatrix,err,error,*999)
      symmetryTypes(matrixIdx)=solverMatrix%symmetryType
    ENDDO !matrixIdx
    
    EXITS("SolverMatrices_SymmetryTypeGet")
    RETURN
999 ERRORSEXITS("SolverMatrices_SymmetryTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMatrices_SymmetryTypeGet

  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to solver distributed matrix for a solver matrix.
  SUBROUTINE SolverMatrix_SolverDistributedMatrixGet(solverMatrix,solverDistributedMatrix,err,error,*)

    !Argument variables
    TYPE(SolverMatrixType), POINTER :: solverMatrix !<A pointer to the solver matrix to get the solver distributed matrix for
    TYPE(DistributedMatrixType), POINTER :: solverMatrix !<On exit, a pointer to the solver distributed matrix for the specified solver matrix. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMatrix_SolverDistributedMatrixGet",err,error,*998)

    IF(ASSOCIATED(solverDistributedMatrix)) CALL FlagError("Solver distributed matrix is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverMatrix)) CALL FlagError("Solver matrix is not associated.",err,error,*999)
    
    solverDistributedMatrix=>solverMatrix%matrix
    IF(.NOT.ASSOCIATED(solverDistributedMatrix)) &
      & CALL FlagError("The solver distributed matrix is not associated for the solver matrix.",err,error,*999)
      
    EXITS("SolverMatrix_SolverDistributedMatrixGet")
    RETURN
999 NULLIFY(solverDistributedMatrix)
998 ERRORSEXITS("SolverMatrix_SolverDistributedMatrixGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMatrix_SolverDistributedMatrixGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to solver distributed vector for a solver matrix.
  SUBROUTINE SolverMatrix_SolverDistributedVectorGet(solverMatrix,solverDistributedVector,err,error,*)

    !Argument variables
    TYPE(SolverMatrixType), POINTER :: solverMatrix !<A pointer to the solver matrix to get the solver distributed vector for
    TYPE(DistributedVectorType), POINTER :: solverVector !<On exit, a pointer to the solver distributed vector for the specified solver matrix. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMatrix_SolverDistributedVectorGet",err,error,*998)

    IF(ASSOCIATED(solverDistributedVector)) CALL FlagError("Solver distributed vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverMatrix)) CALL FlagError("Solver matrix is not associated.",err,error,*999)
    
    solverDistributedVector=>solverMatrix%solverVector
    IF(.NOT.ASSOCIATED(solverDistributedVector)) &
      & CALL FlagError("The solver distributed vector is not associated for the solver matrix.",err,error,*999)
      
    EXITS("SolverMatrix_SolverDistributedVectorGet")
    RETURN
999 NULLIFY(solverDistributedVector)
998 ERRORSEXITS("SolverMatrix_SolverDistributedVectorGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMatrix_SolverDistributedVectorGet
  
  !
  !================================================================================================================================
  !

END MODULE SolverMatricesAccessRoutines
