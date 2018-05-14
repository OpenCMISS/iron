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
  USE Kinds
  USE Strings
  USE Types
  
#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC SolverMatrices_LibraryTypeGet

  PUBLIC SolverMatrices_SolverMappingGet

  PUBLIC SolverMatrices_SolverMatrixGet

  
CONTAINS

  !
  !================================================================================================================================
  !
  
  !>Gets the library type for the solver matrices (and vectors)
  SUBROUTINE SolverMatrices_LibraryTypeGet(solverMatrices,libraryType,err,error,*)

    !Argument variables
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: solverMatrices !<A pointer to the solver matrices.
    INTEGER(INTG), INTENT(OUT) :: libraryType !<On return, the library type of the specified solver matrices \see SOLVER_ROUTINES_SolverLibraries
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverMatrices_LibraryTypeGet",err,error,*999)

    IF(.NOT.ASSOCIATED(solverMatrices)) CALL FlagError("Solver matrices is not associated.",err,error,*999)    
    IF(.NOT.solverMatrices%SOLVER_MATRICES_FINISHED) CALL FlagError("Solver matrices has not finished.",err,error,*999)
    
    libraryType=solverMatrices%solverLibraryType
    
    EXITS("SolverMatrices_LibraryTypeGet")
    RETURN
999 ERRORSEXITS("SolverMatrices_LibraryTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMatrices_LibraryTypeGet
          
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the solver mapping for solver matrices.
  SUBROUTINE SolverMatrices_SolverMappingGet(solverMatrices,solverMapping,err,error,*)

    !Argument variables
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: solverMatrices !<A pointer to the solver matrices to get the solver mappings for
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping !<On exit, a pointer to the solver mapping. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMatrices_SolverMappingGet",err,error,*998)

    IF(ASSOCIATED(solverMapping)) CALL FlagError("Solver mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverMatrices)) CALL FlagError("Solver matrices is not associated.",err,error,*999)

    solverMapping=>solverMatrices%SOLVER_MAPPING
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
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: solverMatrices !<A pointer to the solver matrices to get the solver matrix for
    INTEGER(INTG), INTENT(IN) :: matrixIdx !<The matrix index of the solver matrix to get.
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: solverMatrix !<On exit, a pointer to the solver matrix. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("SolverMatrices_SolverMatrixGet",err,error,*998)

    IF(ASSOCIATED(solverMatrix)) CALL FlagError("Solver matrix is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverMatrices)) CALL FlagError("Solver matrices is not associated.",err,error,*999)
    IF(matrixIdx<1.OR.matrixIdx>solverMatrices%NUMBER_OF_MATRICES) THEN
      localError="The specified solver matrix index of "// &
        & TRIM(NumberToVString(matrixIdx,"*",err,error))//" is invalid. The matrix index needs to be >= 1 and <= "// &
        & TRIM(NumberToVString(solverMatrices%NUMBER_OF_MATRICES,"*",err,error))//"."
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

END MODULE SolverMatricesAccessRoutines
