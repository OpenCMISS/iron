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

  !Module types

  !Module variables

  !Interfaces

  PUBLIC InterfaceMatrices_InterfaceMatrixGet
  
CONTAINS
  
  !
  !================================================================================================================================
  !

  !>Gets the specified interface matrix for interface matrices.
  SUBROUTINE InterfaceMatrices_InterfaceMatrixGet(interfaceMatrices,matrixIdx,interfaceMatrix,err,error,*)

    !Argument variables
    TYPE(INTERFACE_MATRICES_TYPE), POINTER :: interfaceMatrices !<A pointer to the interface matrices to get the interface matrix for
    INTEGER(INTG), INTENT(IN) :: matrixIdx !<The matrix index of the interface matrix to get
    TYPE(INTERFACE_MATRIX_TYPE), POINTER :: interfaceMatrix !<On exit, a pointer to the interface matrix for the matrixIdx'th interface matrix. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("InterfaceMatrices_InterfaceMatrixGet",err,error,*998)

    IF(ASSOCIATED(interfaceMatrix)) CALL FlagError("Interface matrix is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMatrices)) CALL FlagError("Interface matrices is not associated.",err,error,*999)
    IF(matrixIdx<1.OR.matrixIdx>interfaceMatrices%NUMBER_OF_INTERFACE_MATRICES) THEN
      localError="The specified matrix index of "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
        & " is invalid. The matrix index must be >= 1 and <= "// &
        & TRIM(NumberToVString(interfaceMatrices%NUMBER_OF_INTERFACE_MATRICES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(interfaceMatrices%matrices)) CALL FlagError("Interface matrices matrices is not allocated.",err,error,*999)
    
    interfaceMatrix=>interfaceMatrices%matrices(matrixIdx)%ptr
    IF(.NOT.ASSOCIATED(interfaceMatrix)) THEN
      localError="Interface matrix for interface matrix number "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
        & " is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("InterfaceMatrices_InterfaceMatrixGet")
    RETURN
999 NULLIFY(interfaceMatrix)
998 ERRORSEXITS("InterfaceMatrices_InterfaceMatrixGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMatrices_InterfaceMatrixGet

  !
  !================================================================================================================================
  !

END MODULE InterfaceMatricesAccessRoutines
