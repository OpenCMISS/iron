!> \file
!> \author Chris Bradley
!> \brief This module contains all coordinate system access method routines.
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

!> This module contains all coordinate system access method routines.
MODULE CoordinateSystemAccessRoutines
  
  USE BaseRoutines
  USE Kinds
  USE ISO_VARYING_STRING
  USE Strings
  USE Types
  
#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  INTERFACE COORDINATE_SYSTEM_DIMENSION_GET
    MODULE PROCEDURE CoordinateSystem_DimensionGet
  END INTERFACE COORDINATE_SYSTEM_DIMENSION_GET

  INTERFACE COORDINATE_SYSTEM_USER_NUMBER_FIND
    MODULE PROCEDURE CoordinateSystem_UserNumberFind
  END INTERFACE COORDINATE_SYSTEM_USER_NUMBER_FIND

  PUBLIC CoordinateSystem_CoordinateSystemsGet
  
  PUBLIC CoordinateSystem_DimensionGet

  PUBLIC COORDINATE_SYSTEM_DIMENSION_GET

  PUBLIC CoordinateSystem_Get

  PUBLIC CoordinateSystem_UserNumberFind

  PUBLIC COORDINATE_SYSTEM_USER_NUMBER_FIND


CONTAINS

  !
  !================================================================================================================================
  !

  !>Returns the coordinate systems for the coordinate system. 
  SUBROUTINE CoordinateSystem_CoordinateSystemsGet(coordinateSystem,coordinateSystems,err,error,*)

    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER, INTENT(IN) :: coordinateSystem !<A pointer to the coordinate system to get the coordinate systems for
    TYPE(CoordinateSystemsType), POINTER, INTENT(OUT) :: coordinateSystems !<On return, a pointer to the coordinate systems for the coordinate system. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("CoordinateSystem_CoordinateSystemsGet",err,error,*998)

    IF(ASSOCIATED(coordinateSystems)) CALL FlagError("Coordinate systems is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(coordinateSystem)) CALL FlagError("Coordinate system is not associated.",err,error,*999)

    coordinateSystems=>coordinateSystem%coordinateSystems
    IF(.NOT.ASSOCIATED(coordinateSystems)) THEN
      localError="Coordinate systems is not associated for coordinate system number "// &
        & TRIM(NumberToVString(coordinateSystem%USER_NUMBER,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("CoordinateSystem_CoordinateSystemsGet")
    RETURN
998 NULLIFY(coordinateSystems)
999 ERRORS("CoordinateSystem_CoordinateSystemsGet",err,error)
    EXITS("CoordinateSystem_CoordinateSystemsGet")
    RETURN 1

  END SUBROUTINE CoordinateSystem_CoordinateSystemsGet

  !
  !================================================================================================================================
  !

  !>Gets the coordinate system dimension. 
  SUBROUTINE CoordinateSystem_DimensionGet(coordinateSystem,numberOfDimensions,err,error,*)

    !Argument variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: coordinateSystem !<A pointer to the coordinate system to get the dimension for
    INTEGER(INTG), INTENT(OUT) :: numberOfDimensions !<On return, the number of dimensions in the coordinate system.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("CoordinateSystem_DimensionGet",err,error,*999)

    IF(.NOT.ASSOCIATED(coordinateSystem)) CALL FlagError("Coordinate system is not associated.",err,error,*999)
    IF(.NOT.coordinateSystem%COORDINATE_SYSTEM_FINISHED) CALL FlagError("Coordinate system has not been finished.",err,error,*999)

    numberOfDimensions=coordinateSystem%NUMBER_OF_DIMENSIONS
    
    EXITS("CoordinateSystem_DimensionGet")
    RETURN
999 ERRORSEXITS("CoordinateSystem_DimensionGet",err,error)
    RETURN 1

  END SUBROUTINE CoordinateSystem_DimensionGet

  !
  !================================================================================================================================
  !

  !>Finds and returns a pointer to the coordinate system with the given user number. 
  SUBROUTINE CoordinateSystem_Get(coordinateSystems,userNumber,coordinateSystem,err,error,*)

    !Argument variables
    TYPE(CoordinateSystemsType), POINTER :: coordinateSystems !<A pointer to the coordinate systems to get the user number for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the coordinate system to get
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: coordinateSystem !<On exit, a pointer to the coordinate system with the specified user number if it exists. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("CoordinateSystem_Get",err,error,*999)

    CALL CoordinateSystem_UserNumberFind(coordinateSystems,userNumber,coordinateSystem,err,error,*999)
    IF(.NOT.ASSOCIATED(coordinateSystem)) THEN
      localError="A coordinate system with an user number of "//TRIM(NumberToVString(userNumber,"*",err,error))//" does not exist."
      CALL FlagError(localError,err,error,*999)
    ENDIF
  
    EXITS("CoordinateSystem_Get")
    RETURN
999 ERRORSEXITS("CoordinateSystem_Get",err,error)
    RETURN 1
    
  END SUBROUTINE CoordinateSystem_Get

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the coordinate system identified by a user number. If a coordinate system with that number is not
  !>found then coordinate system is set to NULL.
  SUBROUTINE CoordinateSystem_UserNumberFind(coordinateSystems,userNumber,coordinateSystem,err,error,*)

    !Argument variables
    TYPE(CoordinateSystemsType), POINTER :: coordinateSystems !<The coordinate systems to find the user number for.
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the coordinate system to find.
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: coordinateSystem !<On exit, a pointer to the coordinate system with the specified user number if it exists. If no coordinate system has the specified user number the pointer is returned as NULL. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: coordinateSystemIdx
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("CoordinateSystem_UserNumberFind",ERR,ERROR,*999)

    IF(ASSOCIATED(coordinateSystem)) CALL FlagError("Coordinate system is already associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(coordinateSystems)) CALL FlagError("Coordinate systems is not associated.",err,error,*999)
   
    NULLIFY(coordinateSystem)
    IF(ASSOCIATED(coordinateSystems%coordinateSystems)) THEN
      DO coordinateSystemIdx=1,coordinateSystems%numberOfCoordinateSystems
        IF(ASSOCIATED(coordinateSystems%coordinateSystems(coordinateSystemIdx)%ptr)) THEN
          IF(coordinateSystems%coordinateSystems(coordinateSystemIdx)%ptr%USER_NUMBER==userNumber) THEN
            coordinateSystem=>coordinateSystems%coordinateSystems(coordinateSystemIdx)%ptr
            EXIT
          ENDIF
        ELSE
          localError="The coordinate system pointer in coordinate systems is not associated for coordinate system index "// &
            & TRIM(NumberToVString(coordinateSystemIdx,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDDO !coordinateSystemIdx
    ENDIF
    
    EXITS("CoordinateSystem_UserNumberFind")
    RETURN
999 ERRORS("CoordinateSystem_UserNumberFind",err,error)
    EXITS("CoordinateSystem_UserNumberFind")
    RETURN 1
    
  END SUBROUTINE CoordinateSystem_UserNumberFind

  !
  !================================================================================================================================
  !
  
END MODULE CoordinateSystemAccessRoutines
