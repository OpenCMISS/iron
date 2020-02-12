!> \file
!> \author Chris Bradley
!> \brief This module contains all interface matrices access method routines.
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

!> This module contains all interface matrices access method routines.
MODULE InterfaceMatricesAccessRoutines
  
  USE BaseRoutines
  USE ISO_VARYING_STRING
  USE Kinds
  USE Strings
  USE Types

#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup InterfaceMatricesRoutines_InterfaceMatrixStructureTypes InterfaceMatricesRoutines::InterfaceMatrixStructureTypes
  !> \brief Interface matrices structure (sparsity) types
  !> \see InterfaceMatricesRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: INTERFACE_MATRIX_NO_STRUCTURE=1 !<No matrix structure - all elements can contain a value. \see InterfaceMatricesRoutines_InterfaceMatrixStructureTypes,InterfaceMatricesRoutines
  INTEGER(INTG), PARAMETER :: INTERFACE_MATRIX_FEM_STRUCTURE=2 !<Finite element matrix structure. \see InterfaceMatricesRoutines_InterfaceMatrixStructureTypes,InterfaceMatricesRoutines 
  !>@}

  !> \addtogroup InterfaceMatricesRoutines_InterfaceMatricesSparsityTypes InterfaceMatricesRoutines::InterfaceMatricesSparsityTypes
  !> \brief Interface matrices sparsity types
  !> \see InterfaceMatricesRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: INTERFACE_MATRICES_SPARSE_MATRICES=1 !<Use sparse interface matrices \see InterfaceMatricesRoutines_InterfaceMatricesSparsityTypes,InterfaceMatricesRoutines
  INTEGER(INTG), PARAMETER :: INTERFACE_MATRICES_FULL_MATRICES=2 !<Use fully populated interface matrices \see InterfaceMatricesRoutines_InterfaceMatricesSparsityTypes,InterfaceMatricesRoutines
  !>@}

  !> \addtogroup InterfaceMatricesRoutines_InterfaceMatricesTimeDependenceTypes InterfaceMatricesRoutines::InterfaceMatricesTimeDependenceTypes
  !> \brief Interface matrices time dependency types
  !> \see InterfaceMatricesRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: NUMBER_OF_INTERFACE_MATRIX_TYPES=4
  INTEGER(INTG), PARAMETER :: INTERFACE_MATRIX_STATIC=1 !<Interface matrix is of static type \see InterfaceMatricesRoutines_InterfaceMatricesTimeDependenceTypes,InterfaceMatricesRoutines
  INTEGER(INTG), PARAMETER :: INTERFACE_MATRIX_QUASI_STATIC=2 !<Interface matrix is of quasi-static type \see InterfaceMatricesRoutines_InterfaceMatricesTimeDependenceTypes,InterfaceMatricesRoutines
  INTEGER(INTG), PARAMETER :: INTERFACE_MATRIX_FIRST_ORDER_DYNAMIC=3 !<Interface matrix is of first order dynamic type \see InterfaceMatricesRoutines_InterfaceMatricesTimeDependenceTypes,InterfaceMatricesRoutines
  INTEGER(INTG), PARAMETER :: INTERFACE_MATRIX_SECOND_ORDER_DYNAMIC=4 !<Interface matrix is of second order dynamic type \see InterfaceMatricesRoutines_InterfaceMatricesTimeDependenceTypes,InterfaceMatricesRoutines
  !>@}

  !Module types

  !Module variables

  !Interfaces

  PUBLIC INTERFACE_MATRIX_NO_STRUCTURE,INTERFACE_MATRIX_FEM_STRUCTURE

  PUBLIC INTERFACE_MATRICES_SPARSE_MATRICES,INTERFACE_MATRICES_FULL_MATRICES

  PUBLIC INTERFACE_MATRIX_STATIC,INTERFACE_MATRIX_QUASI_STATIC,INTERFACE_MATRIX_FIRST_ORDER_DYNAMIC, &
    & INTERFACE_MATRIX_SECOND_ORDER_DYNAMIC,NUMBER_OF_INTERFACE_MATRIX_TYPES

  PUBLIC InterfaceMatrices_AssertIsFinished,InterfaceMatrices_AssertNotFinished

  PUBLIC InterfaceMatrices_InterfaceEquationsGet

  PUBLIC InterfaceMatrices_InterfaceMappingGet

  PUBLIC InterfaceMatrices_InterfaceMatrixGet
  
  PUBLIC InterfaceMatrices_RHSVectorExists

  PUBLIC InterfaceMatrices_RHSVectorGet

  PUBLIC InterfaceMatricesRHS_DistributedVectorGet

  PUBLIC InterfaceMatricesRHS_InterfaceMatricesGet

  PUBLIC InterfaceMatrix_DistributedMatrixGet
  
  PUBLIC InterfaceMatrix_InterfaceMatricesGet

  PUBLIC InterfaceMatrix_TempDistributedVectorGet
  
  PUBLIC InterfaceMatrix_TempTransposeDistributedVectorGet
  
  PUBLIC InterfaceMatrix_TransposeDistributedMatrixGet
  
CONTAINS
  
  !
  !================================================================================================================================
  !

  !>Assert that an interface matrices has been finished
  SUBROUTINE InterfaceMatrices_AssertIsFinished(interfaceMatrices,err,error,*)

    !Argument Variables
    TYPE(InterfaceMatricesType), POINTER, INTENT(IN) :: interfaceMatrices !<The interface matrices to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMatrices_AssertIsFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(interfaceMatrices)) CALL FlagError("Interface matrices is not associated.",err,error,*999)
#endif    

    IF(.NOT.interfaceMatrices%interfaceMatricesFinished) CALL FlagError("Interface matrices has not been finished.",err,error,*999)
    
    EXITS("InterfaceMatrices_AssertIsFinished")
    RETURN
999 ERRORSEXITS("InterfaceMatrices_AssertIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMatrices_AssertIsFinished

  !
  !================================================================================================================================
  !

  !>Assert that an interface matrices has not been finished
  SUBROUTINE InterfaceMatrices_AssertNotFinished(interfaceMatrices,err,error,*)

    !Argument Variables
    TYPE(InterfaceMatricesType), POINTER, INTENT(IN) :: interfaceMatrices !<The interface matrices to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMatrices_AssertNotFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(interfaceMatrices)) CALL FlagError("Interface matrices is not associated.",err,error,*999)
#endif    

    IF(interfaceMatrices%interfaceMatricesFinished) CALL FlagError("Interface matrices has already been finished.",err,error,*999)
    
    EXITS("InterfaceMatrices_AssertNotFinished")
    RETURN
999 ERRORSEXITS("InterfaceMatrices_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMatrices_AssertNotFinished

  !
  !================================================================================================================================
  !

  !>Gets the interface equations for an interface matrices.
  SUBROUTINE InterfaceMatrices_InterfaceEquationsGet(interfaceMatrices,interfaceEquations,err,error,*)

    !Argument variables
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices !<A pointer to the interface matrices to get the interface equations for
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations !<On exit, a pointer to the interface equations in the specified interface matrices. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMatrices_InterfaceEquationsGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(interfaceEquations)) CALL FlagError("Interface equations is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMatrices)) CALL FlagError("Interface matrices is not associated.",err,error,*999)
#endif    

    interfaceEquations=>interfaceMatrices%interfaceEquations

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(interfaceEquations)) &
      & CALL FlagError("Interface matrices interface equations is not associated.",err,error,*999)
#endif    
       
    EXITS("InterfaceMatrices_InterfaceEquationsGet")
    RETURN
999 NULLIFY(interfaceEquations)
998 ERRORSEXITS("InterfaceMatrices_InterfaceEquationsGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMatrices_InterfaceEquationsGet

  !
  !================================================================================================================================
  !

  !>Gets the interface mapping for an interface matrices.
  SUBROUTINE InterfaceMatrices_InterfaceMappingGet(interfaceMatrices,interfaceMapping,err,error,*)

    !Argument variables
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices !<A pointer to the interface matrices to get the interface mapping for
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping !<On exit, a pointer to the interface mapping in the specified interface matrices. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMatrices_InterfaceMappingGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(interfaceMapping)) CALL FlagError("Interface mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMatrices)) CALL FlagError("Interface matrices is not associated.",err,error,*999)
#endif    

    interfaceMapping=>interfaceMatrices%interfaceMapping

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(interfaceMapping)) &
      & CALL FlagError("Interface matrices interface mapping is not associated.",err,error,*999)
#endif    
       
    EXITS("InterfaceMatrices_InterfaceMappingGet")
    RETURN
999 NULLIFY(interfaceMapping)
998 ERRORSEXITS("InterfaceMatrices_InterfaceMappingGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMatrices_InterfaceMappingGet

  !
  !================================================================================================================================
  !

  !>Gets the specified interface matrix for interface matrices.
  SUBROUTINE InterfaceMatrices_InterfaceMatrixGet(interfaceMatrices,matrixIdx,interfaceMatrix,err,error,*)

    !Argument variables
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices !<A pointer to the interface matrices to get the interface matrix for
    INTEGER(INTG), INTENT(IN) :: matrixIdx !<The matrix index of the interface matrix to get
    TYPE(InterfaceMatrixType), POINTER :: interfaceMatrix !<On exit, a pointer to the interface matrix for the matrixIdx'th interface matrix. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("InterfaceMatrices_InterfaceMatrixGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(interfaceMatrix)) CALL FlagError("Interface matrix is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMatrices)) CALL FlagError("Interface matrices is not associated.",err,error,*999)
    IF(matrixIdx<1.OR.matrixIdx>interfaceMatrices%numberOfInterfaceMatrices) THEN
      localError="The specified matrix index of "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
        & " is invalid. The matrix index must be >= 1 and <= "// &
        & TRIM(NumberToVString(interfaceMatrices%numberOfInterfaceMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(interfaceMatrices%matrices)) CALL FlagError("Interface matrices matrices is not allocated.",err,error,*999)
#endif    
    
    interfaceMatrix=>interfaceMatrices%matrices(matrixIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(interfaceMatrix)) THEN
      localError="Interface matrix for interface matrix number "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
        & " is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    EXITS("InterfaceMatrices_InterfaceMatrixGet")
    RETURN
999 NULLIFY(interfaceMatrix)
998 ERRORSEXITS("InterfaceMatrices_InterfaceMatrixGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMatrices_InterfaceMatrixGet

  !
  !================================================================================================================================
  !

  !>Checks if the RHS vector exists for an interface matrices.
  SUBROUTINE InterfaceMatrices_RHSVectorExists(interfaceMatrices,rhsVector,err,error,*)

    !Argument variables
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices !<A pointer to the interface matrices to check the RHS vector for
    TYPE(InterfaceRHSType), POINTER :: rhsVector !<On exit, a pointer to the RHS vector in the specified interface matrices if it exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMatrices_RHSVectorExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(rhsVector)) CALL FlagError("RHS vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMatrices)) CALL FlagError("Interface matrices is not associated.",err,error,*999)
#endif    

    rhsVector=>interfaceMatrices%rhsVector
       
    EXITS("InterfaceMatrices_RHSVectorExists")
    RETURN
999 NULLIFY(rhsVector)
998 ERRORSEXITS("InterfaceMatrices_RHSVectorExists",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMatrices_RHSVectorExists

  !
  !================================================================================================================================
  !

  !>Gets the RHS vector for an interface matrices.
  SUBROUTINE InterfaceMatrices_RHSVectorGet(interfaceMatrices,rhsVector,err,error,*)

    !Argument variables
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices !<A pointer to the interface matrices to get the RHS vector for
    TYPE(InterfaceRHSType), POINTER :: rhsVector !<On exit, a pointer to the RHS vector in the specified interface matrices. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMatrices_RHSVectorGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(rhsVector)) CALL FlagError("RHS vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMatrices)) CALL FlagError("Interface matrices is not associated.",err,error,*999)
#endif    

    rhsVector=>interfaceMatrices%rhsVector

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(rhsVector)) CALL FlagError("Interface matrices RHS vector is not associated.",err,error,*999)
#endif    
       
    EXITS("InterfaceMatrices_RHSVectorGet")
    RETURN
999 NULLIFY(rhsVector)
998 ERRORSEXITS("InterfaceMatrices_RHSVectorGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMatrices_RHSVectorGet

  !
  !================================================================================================================================
  !

  !>Gets the RHS distributed vector for an interface matrices RHS.
  SUBROUTINE InterfaceMatricesRHS_DistributedVectorGet(interfaceMatricesRHS,distributedVector,err,error,*)

    !Argument variables
    TYPE(InterfaceRHSType), POINTER :: interfaceMatricesRHS !<A pointer to the interface matrices RHS to get the RHS distributed vector for
    TYPE(DistributedVectorType), POINTER :: distrbutedVector !<On exit, a pointer to the RHS distributed vector in the specified interface matrices RHS. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMatricesRHS_DistributedVectorGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(distributedVector)) CALL FlagError("Distrbuted vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMatricesRHS)) CALL FlagError("Interface matrices RHS is not associated.",err,error,*999)
#endif    

    distrbutedVector=>interfaceMatricesRHS%vector

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(distributedVector)) &
      & CALL FlagError("Interface matrices RHS distributed vector is not associated.",err,error,*999)
#endif    
       
    EXITS("InterfaceMatricesRHS_DistributedVectorGet")
    RETURN
999 NULLIFY(distributedVector)
998 ERRORSEXITS("InterfaceMatricesRHS_DistributedVectorGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMatricesRHS_DistributedVectorGet

  !
  !================================================================================================================================
  !

  !>Gets the interface matrices for an interface matrices RHS.
  SUBROUTINE InterfaceMatricesRHS_InterfaceMatricesGet(interfaceMatricesRHS,interfaceMatrices,err,error,*)

    !Argument variables
    TYPE(InterfaceRHSType), POINTER :: interfaceMatricesRHS !<A pointer to the interface matrices RHS to get the interface matrices for
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices !<On exit, a pointer to the interface matrices in the specified interface matrices RHS. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMatricesRHS_InterfaceMatricesGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(interfaceMatrices)) CALL FlagError("Interface matrices is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMatricesRHS)) CALL FlagError("Interface matrices RHS is not associated.",err,error,*999)
#endif    

    interfaceMatrices=>interfaceMatricesRHS%interfaceMatrices

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(interfaceMatrices)) &
      & CALL FlagError("Interface matrices RHS interface matrices is not associated.",err,error,*999)
#endif    
       
    EXITS("InterfaceMatricesRHS_InterfaceMatricesGet")
    RETURN
999 NULLIFY(interfaceMatrices)
998 ERRORSEXITS("InterfaceMatricesRHS_InterfaceMatricesGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMatricesRHS_InterfaceMatricesGet

  !
  !================================================================================================================================
  !

  !>Gets the distributed matrix for an interface matrix.
  SUBROUTINE InterfaceMatrix_DistributedMatrixGet(interfaceMatrix,distributedMatrix,err,error,*)

    !Argument variables
    TYPE(InterfaceMatrixType), POINTER :: interfaceMatrix !<A pointer to the interface matrix to get the distributed matrix for
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<On exit, a pointer to the distributed matrix in the specified interface matrix. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMatrix_DistributedMatrixGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMatrix)) CALL FlagError("Interface matrix is not associated.",err,error,*999)
#endif    

    distributedMatrix=>interfaceMatrix%matrix

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(distributedMatrix)) CALL FlagError("Interface matrix distributed matrix is not associated.",err,error,*999)
#endif    
       
    EXITS("InterfaceMatrix_DistributedMatrixGet")
    RETURN
999 NULLIFY(distributedMatrix)
998 ERRORSEXITS("InterfaceMatrices_DistributedMatrixGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMatrix_DistributedMatrixGet

  !
  !================================================================================================================================
  !

  !>Gets the interface matrices for an interface matrix.
  SUBROUTINE InterfaceMatrix_InterfaceMatricesGet(interfaceMatrix,interfaceMatrices,err,error,*)

    !Argument variables
    TYPE(InterfaceMatrixType), POINTER :: interfaceMatrix !<A pointer to the interface matrix to get the interfaces matrices for
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices !<On exit, a pointer to the interface matrices in the specified interface matrix. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMatrix_InterfaceMatricesGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(interfaceMatrices)) CALL FlagError("Interface matrices is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMatrix)) CALL FlagError("Interface matrix is not associated.",err,error,*999)
#endif    

    interfaceMatrices=>interfaceMatrix%interfaceMatrices

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(interfaceMatrices)) CALL FlagError("Interface matrix interface matrices is not associated.",err,error,*999)
#endif    
       
    EXITS("InterfaceMatrix_InterfaceMatricesGet")
    RETURN
999 NULLIFY(interfaceMatrices)
998 ERRORSEXITS("InterfaceMatrices_InterfaceMatricesGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMatrix_InterfaceMatricesGet

  !
  !================================================================================================================================
  !

  !>Gets the temporary distributed vector for an interface matrix.
  SUBROUTINE InterfaceMatrix_TempDistributedVectorGet(interfaceMatrix,temporaryDistributedVector,err,error,*)

    !Argument variables
    TYPE(InterfaceMatrixType), POINTER :: interfaceMatrix !<A pointer to the interface matrix to get the temporary distributed vector for
    TYPE(DistributedVectorType), POINTER :: temporaryDistributedVector !<On exit, a pointer to the temporary distributed vector in the specified interface matrix. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMatrix_TempDistributedVectorGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(temporaryDistributedVector)) CALL FlagError("Temporary distributed vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMatrix)) CALL FlagError("Interface matrix is not associated.",err,error,*999)
#endif    

    temporaryDistributedVector=>interfaceMatrix%tempVector

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(temporaryDistributedVector)) &
      & CALL FlagError("Interface matrix temporary distributed vector is not associated.",err,error,*999)
#endif    
       
    EXITS("InterfaceMatrix_TempDistributedVectorGet")
    RETURN
999 NULLIFY(temporaryDistributedVector)
998 ERRORSEXITS("InterfaceMatrices_TempDistributedVectorGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMatrix_TempDistributedVectorGet

  !
  !================================================================================================================================
  !

  !>Gets the temporary transpose distributed vector for an interface matrix.
  SUBROUTINE InterfaceMatrix_TempTransposeDistributedVectorGet(interfaceMatrix,tempTransDistributedVector,err,error,*)

    !Argument variables
    TYPE(InterfaceMatrixType), POINTER :: interfaceMatrix !<A pointer to the interface matrix to get the temporary transpose distributed vector for
    TYPE(DistributedVectorType), POINTER :: tempTransDistributedVector !<On exit, a pointer to the temporary transpose distributed vector in the specified interface matrix. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMatrix_TempTranposeDistributedVectorGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(tempTransDistributedVector)) &
      & CALL FlagError("Temporary transpose distributed vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMatrix)) CALL FlagError("Interface matrix is not associated.",err,error,*999)
#endif    

    tempTransDistributedVector=>interfaceMatrix%tempTranposeVector

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(tempTransDistributedVector)) &
      & CALL FlagError("Interface matrix temporary transpose distributed vector is not associated.",err,error,*999)
#endif    
       
    EXITS("InterfaceMatrix_TempTransposeDistributedVectorGet")
    RETURN
999 NULLIFY(tempTransDistributedVector)
998 ERRORSEXITS("InterfaceMatrices_TempTransposeDistributedVectorGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMatrix_TempTransposeDistributedVectorGet

  !
  !================================================================================================================================
  !

  !>Gets the transpose distributed matrix for an interface matrix.
  SUBROUTINE InterfaceMatrix_TransposeDistributedMatrixGet(interfaceMatrix,transposeDistributedMatrix,err,error,*)

    !Argument variables
    TYPE(InterfaceMatrixType), POINTER :: interfaceMatrix !<A pointer to the interface matrix to get the tranpsoe distributed matrix for
    TYPE(DistributedMatrixType), POINTER :: transposeDistributedMatrix !<On exit, a pointer to the transpose distributed matrix in the specified interface matrix. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMatrix_TransposeDistributedMatrixGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(transposeDistributedMatrix)) CALL FlagError("Transpose distributed matrix is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMatrix)) CALL FlagError("Interface matrix is not associated.",err,error,*999)
#endif    

    transposeDistributedMatrix=>interfaceMatrix%tranposeMatrix

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(transposeDistributedMatrix)) &
      & CALL FlagError("Interface matrix transpose distributed matrix is not associated.",err,error,*999)
#endif    
       
    EXITS("InterfaceMatrix_TransposeDistributedMatrixGet")
    RETURN
999 NULLIFY(tranposeDistributedMatrix)
998 ERRORSEXITS("InterfaceMatrices_TranposeDistributedMatrixGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMatrix_TransposeDistributedMatrixGet

  !
  !================================================================================================================================
  !

END MODULE InterfaceMatricesAccessRoutines
