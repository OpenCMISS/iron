!> \file
!> \author Chris Bradley
!> \brief This module contains all interface equations access method routines.
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

!> This module contains all interface equations access method routines.
MODULE InterfaceEquationsAccessRoutines
  
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

  PUBLIC InterfaceEquations_InterfaceMatricesGet

CONTAINS

  !
  !================================================================================================================================
  !

  !>Gets the interface matrices for an interface equations.
  SUBROUTINE InterfaceEquations_InterfaceMatricesGet(interfaceEquations,interfaceMatrices,err,error,*)

    !Argument variables
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: interfaceEquations !<A pointer to the interface equations to get the interface matrices for
    TYPE(INTERFACE_MATRICES_TYPE), POINTER :: interfaceMatrices !<On exit, a pointer to the interface matrices for the specified interface equations. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceEquations_InterfaceMatricesGet",err,error,*999)

    IF(.NOT.ASSOCIATED(interfaceEquations)) CALL FlagError("Interface equations is not associated.",err,error,*999)
    IF(ASSOCIATED(interfaceMatrices)) CALL FlagError("Interface matrices is already associated.",err,error,*999)

    interfaceMatrices=>interfaceEquations%INTERFACE_MATRICES
    IF(.NOT.ASSOCIATED(interfaceMatrices)) &
      & CALL FlagError("Interface equations interface matrices is not associated.",err,error,*999)
       
    EXITS("InterfaceEquations_InterfaceMatricesGet")
    RETURN
999 ERRORSEXITS("InterfaceEquations_InterfaceMatricesGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceEquations_InterfaceMatricesGet


END MODULE InterfaceEquationsAccessRoutines
