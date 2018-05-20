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

  PUBLIC DistributedMatrix_CMISSMatrixGet
  
  PUBLIC DistributedMatrix_ColumnMappingGet
  
  PUBLIC DistributedMatrix_PETScMatrixGet
    
  PUBLIC DistributedMatrix_RowMappingGet

  PUBLIC DistributedMatrixCMISS_DistributedMatrixGet

  PUBLIC DistributedMatrixPetsc_DistributedMatrixGet

  PUBLIC DistributedVector_CMISSVectorGet
  
  PUBLIC DistributedVector_PETScVectorGet
  
  PUBLIC DistributedVector_RowMappingGet

  PUBLIC DistributedVectorCMISS_DistributedVectorGet
  
  PUBLIC DistributedVectorPETSc_DistributedVectorGet  

CONTAINS

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

    IF(ASSOCIATED(cmissMatrix)) CALL FlagError("CMISS matrix is associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)    
    
    cmissMatrix=>distributedMatrix%cmiss
    IF(.NOT.ASSOCIATED(cmissMatrix))  CALL FlagError("CMISS matrix is not associated for the distributed matrix.",err,error,*999)
     
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
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: columnMapping !<On return, the column mapping for the distributed matrix. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DistributedMatrix_ColumnMappingGet",err,error,*998)

    IF(ASSOCIATED(columnMapping)) CALL FlagError("Column mapping is associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)    
    
    columnMapping=>distributedMatrix%columnDomainMapping
    IF(.NOT.ASSOCIATED(columnMapping)) CALL FlagError("Column mapping is not associated for the distributed matrix.",err,error,*999)
     
    EXITS("DistributedMatrix_ColumnMappingGet")
    RETURN
999 NULLIFY(columnMapping)
998 ERRORSEXITS("DistributedMatrix_ColumnMappingGet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_ColumnMappingGet
  
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

    IF(ASSOCIATED(petscMatrix)) CALL FlagError("PETSc matrix is associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)    
    
    petscMatrix=>distributedMatrix%petsc
    IF(.NOT.ASSOCIATED(petscMatrix))  CALL FlagError("PETSc matrix is not associated for the distributed matrix.",err,error,*999)
     
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
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: rowMapping !<On return, the row mapping for the distributed matrix. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DistributedMatrix_RowMappingGet",err,error,*998)

    IF(ASSOCIATED(rowMapping)) CALL FlagError("Row mapping is associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is not associated.",err,error,*999)    
    
    rowMapping=>distributedMatrix%rowDomainMapping
    IF(.NOT.ASSOCIATED(rowMapping))  CALL FlagError("Row mapping is not associated for the distributed matrix.",err,error,*999)
     
    EXITS("DistributedMatrix_RowMappingGet")
    RETURN
999 NULLIFY(rowMapping)
998 ERRORSEXITS("DistributedMatrix_RowMappingGet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrix_RowMappingGet
  
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

    IF(ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is associated.",err,error,*998)    
    IF(.NOT.ASSOCIATED(cmissMatrix)) CALL FlagError("CMISS matrix is not associated.",err,error,*999)
    
    distributedMatrix=>cmissMatrix%distributedMatrix
    IF(.NOT.ASSOCIATED(distributedMatrix)) &
      & CALL FlagError("Distributed matrix is not associated for the CMISS matrix.",err,error,*999)
     
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

    IF(ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is associated.",err,error,*998)    
    IF(.NOT.ASSOCIATED(petscMatrix)) CALL FlagError("PETSc matrix is not associated.",err,error,*999)
    
    distributedMatrix=>petscMatrix%distributedMatrix
    IF(.NOT.ASSOCIATED(distributedMatrix)) &
      & CALL FlagError("Distributed matrix is not associated for the PETSc matrix.",err,error,*999)
     
    EXITS("DistributedMatrixPETSc_DistributedMatrixGet")
    RETURN
999 NULLIFY(distributedMatrix)
998 ERRORSEXITS("DistributedMatrixPETSc_DistributedMatrixGet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedMatrixPETSc_DistributedMatrixGet
  
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

    IF(ASSOCIATED(cmissVector)) CALL FlagError("CMISS vector is associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)    
    
    cmissVector=>distributedVector%cmiss
    IF(.NOT.ASSOCIATED(cmissVector))  CALL FlagError("CMISS vector is not associated for the distributed vector.",err,error,*999)
     
    EXITS("DistributedVector_CMISSVectorGet")
    RETURN
999 NULLIFY(cmissVector)
998 ERRORSEXITS("DistributedVector_CMISSVectorGet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_CMISSVectorGet
  
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

    IF(ASSOCIATED(petscVector)) CALL FlagError("PETSc vector is associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)    
    
    petscVector=>distributedVector%petsc
    IF(.NOT.ASSOCIATED(petscVector))  CALL FlagError("PETSc vector is not associated for the distributed vector.",err,error,*999)
     
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
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: rowMapping !<On return, the row mapping for the distributed vector. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DistributedVector_RowMappingGet",err,error,*998)

    IF(ASSOCIATED(rowMapping)) CALL FlagError("Row mapping is associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is not associated.",err,error,*999)    
    
    rowMapping=>distributedVector%domainMapping
    IF(.NOT.ASSOCIATED(rowMapping)) CALL FlagError("Row mapping is not associated for the distributed vector.",err,error,*999)
     
    EXITS("DistributedVector_RowMappingGet")
    RETURN
999 NULLIFY(rowMapping)
998 ERRORSEXITS("DistributedVector_RowMappingGet",err,error)
    RETURN 1
    
  END SUBROUTINE DistributedVector_RowMappingGet
  
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

    IF(ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(cmissVector)) CALL FlagError("CMISS vector is not associated.",err,error,*999)    
    
    distributedVector=>cmissVector%distributedVector
    IF(.NOT.ASSOCIATED(distributedVector)) &
      & CALL FlagError("Distributed vector is not associated for the CMISS vector.",err,error,*999)
     
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

    IF(ASSOCIATED(distributedVector)) CALL FlagError("Distributed vector is associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(petscVector)) CALL FlagError("PETSc vector is not associated.",err,error,*999)    
    
    distributedVector=>petscVector%distributedVector
    IF(.NOT.ASSOCIATED(distributedVector)) &
      & CALL FlagError("Distributed vector is not associated for the PETSc vector.",err,error,*999)
     
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
