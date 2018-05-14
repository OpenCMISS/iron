!> \file
!> \author Chris Bradley
!> \brief This module contains all CellML access method routines.
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

!> This module contains all CellML access method routines.
MODULE CellMLAccessRoutines
  
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

  INTERFACE CELLML_USER_NUMBER_FIND
    MODULE PROCEDURE CellML_UserNumberFind
  END INTERFACE CELLML_USER_NUMBER_FIND

  PUBLIC CellML_UserNumberFind

  PUBLIC CELLML_USER_NUMBER_FIND

CONTAINS

  !
  !=================================================================================================================================
  !

  !>Finds and returns a pointer to the CellML environment identified by a user number on a region. If no CellML environment with that user number exists cellml is left nullified.
  SUBROUTINE CellML_UserNumberFind(userNumber,region,cellml,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number to find.
    TYPE(REGION_TYPE), POINTER :: region !<A pointer to the region to find the CellML user number.
    TYPE(CELLML_TYPE), POINTER :: cellml !<On return a pointer to the CellML environment with the given user number. If no CellML environment with that user number exists then the pointer is returned as NULL. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: cellmlIdx
    TYPE(VARYING_STRING) :: localError

    ENTERS("CellML_UserNumberFind",err,error,*999)

    IF(.NOT.ASSOCIATED(region)) CALL FlagError("Region is not associated.",err,error,*999)
    IF(ASSOCIATED(cellml)) CALL FlagError("CellML is already associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(region%CELLML_ENVIRONMENTS)) CALL FlagError("Region CellML environments is not associated.",err,error,*999)

    NULLIFY(cellml)
    IF(ALLOCATED(region%CELLML_ENVIRONMENTS%environments)) THEN
      DO cellmlIdx=1,region%CELLML_ENVIRONMENTS%NUMBER_OF_ENVIRONMENTS
        IF(ASSOCIATED(region%CELLML_ENVIRONMENTS%environments(cellmlIdx)%ptr)) THEN
          IF(region%CELLML_ENVIRONMENTS%environments(cellmlIdx)%ptr%USER_NUMBER==userNumber) THEN
            cellml=>region%CELLML_ENVIRONMENTS%environments(cellmlIdx)%ptr
            EXIT
          ENDIF
        ELSE
          localError="The cellml pointer in cellml environments is not associated for cellml index "// &
            & TRIM(NumberToVString(cellmlIdx,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDDO !cellmlIdx
    ENDIF

    EXITS("CellML_UserNumberFind")
    RETURN
999 ERRORSEXITS("CellML_UserNumberFind",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_UserNumberFind

  !
  !================================================================================================================================
  !

END MODULE CellMLAccessRoutines
