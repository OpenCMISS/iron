!> \file
!> \author Chris Bradley
!> \brief This module contains all equations set access method routines.
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

!> This module contains all equations set access method routines.
MODULE EquationsSetAccessRoutines
  
  USE BASE_ROUTINES
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

  INTERFACE EQUATIONS_SET_EQUATIONS_GET
    MODULE PROCEDURE EquationsSet_EquationsGet
  END INTERFACE EQUATIONS_SET_EQUATIONS_GET
  
  INTERFACE EQUATIONS_SET_USER_NUMBER_FIND
    MODULE PROCEDURE EquationsSet_UserNumberFind
  END INTERFACE EQUATIONS_SET_USER_NUMBER_FIND

  PUBLIC EquationsSet_EquationsGet

  PUBLIC EQUATIONS_SET_EQUATIONS_GET
  
  PUBLIC EquationsSet_UserNumberFind

  PUBLIC EQUATIONS_SET_USER_NUMBER_FIND

CONTAINS

  !
  !================================================================================================================================
  !

  !>Gets the equations for an equations set.
  SUBROUTINE EquationsSet_EquationsGet(equationsSet,equations,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to get the equations for
    TYPE(EQUATIONS_TYPE), POINTER :: equations !<On exit, a pointer to the equations in the specified equations set. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_EquationsGet",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(.NOT.equationsSet%EQUATIONS_SET_FINISHED) CALL FlagError("Equations set has not been finished.",err,error,*999)
    IF(ASSOCIATED(equations)) CALL FlagError("Equations is already associated.",err,error,*999)

    equations=>equationsSet%equations
    IF(.NOT.ASSOCIATED(equations)) THEN
      localError="Equations is not associated for equations set number "// &
      & TRIM(NumberToVString(equationsSet%USER_NUMBER,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("EquationsSet_EquationsGet")
    RETURN
999 ERRORSEXITS("EquationsSet_EquationsGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_EquationsGet

  !
  !================================================================================================================================
  !

  !>Finds and returns a pointer to the equations set identified by user number in the given region. If no equations set with that user number exists the equations set is left nullified.
  SUBROUTINE EquationsSet_UserNumberFind(userNumber,region,equationsSet,err,error,*)

    !Argument variables 
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the equation set to find.
    TYPE(REGION_TYPE), POINTER :: region !<The region to find the equations set in
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<On return, a pointer to the equations set if an equations set with the specified user number exists in the given region. If no equation set with the specified number exists a NULL pointer is returned. The pointer must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsSet_UserNumberFind",err,error,*999)

    IF(.NOT.ASSOCIATED(region)) CALL FlagError("Region is not associated.",err,error,*999)
    IF(ASSOCIATED(equationsSet)) CALL FlagError("Equations set is already associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(region%EQUATIONS_SETS)) THEN
      localError="The equations sets on region number "//TRIM(NumberToVString(region%USER_NUMBER,"*",err,error))// &
        & " are not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    NULLIFY(equationsSet)
    IF(ASSOCIATED(region%EQUATIONS_SETS%EQUATIONS_SETS)) THEN
      DO equationsSetIdx=1,region%EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS
        IF(ASSOCIATED(region%EQUATIONS_SETS%EQUATIONS_SETS(equationsSetIdx)%ptr)) THEN
          IF(region%EQUATIONS_SETS%EQUATIONS_SETS(equationsSetIdx)%ptr%USER_NUMBER==userNumber) THEN
            equationsSet=>region%EQUATIONS_SETS%EQUATIONS_SETS(equationsSetIdx)%ptr
            EXIT
          ENDIF
        ELSE
          localError="The equations set pointer in region equations sets is not associated for equations set index "// &
            & TRIM(NumberToVString(equationsSetIdx,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDDO !equationsSetIdx
    ENDIF
    
    EXITS("EquationsSet_UserNumberFind")
    RETURN
999 ERRORSEXITS("EquationsSet_UserNumberFind",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_UserNumberFind

  !
  !================================================================================================================================
  !

END MODULE EquationsSetAccessRoutines
