!> \file
!> \author Chris Bradley
!> \brief This module contains all data point access method routines.
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

!> This module contains all data point access method routines.
MODULE DataPointAccessRoutines
  
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

  PUBLIC DataPoints_AssertIsFinished,DataPoints_AssertNotFinished

  PUBLIC DataPoints_CoordinateSystemGet

  PUBLIC DataPointSets_UserNumberFind

CONTAINS

  !
  !=================================================================================================================================
  !

  !>Assert that a data points has been finished
  SUBROUTINE DataPoints_AssertIsFinished(dataPoints,err,error,*)

    !Argument Variables
    TYPE(DataPointsType), POINTER, INTENT(IN) :: dataPoints !<The data points to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("DataPoints_AssertIsFinished",err,error,*999)

    IF(.NOT.ASSOCIATED(dataPoints)) CALL FlagError("Data points is not associated.",err,error,*999)

    IF(.NOT.dataPoints%dataPointsFinished) THEN
      localError="Data points number "//TRIM(NumberToVString(dataPoints%userNumber,"*",err,error))
      IF(ASSOCIATED(dataPoints%region)) localError=localError// &
        & " for region number "//TRIM(NumberToVString(dataPoints%region%userNumber,"*",err,error))
      IF(ASSOCIATED(dataPoints%interface)) localError=localError// &
        & " for interface number "//TRIM(NumberToVString(dataPoints%interface%userNumber,"*",err,error))
      localError=localError//" has not been finished."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("DataPoints_AssertIsFinished")
    RETURN
999 ERRORSEXITS("DataPoints_AssertIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE DataPoints_AssertIsFinished

  !
  !=================================================================================================================================
  !

  !>Assert that a data points has not been finished
  SUBROUTINE DataPoints_AssertNotFinished(dataPoints,err,error,*)

    !Argument Variables
    TYPE(DataPointsType), POINTER, INTENT(IN) :: dataPoints !<The data points to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("DataPoints_AssertNotFinished",err,error,*999)

    IF(.NOT.ASSOCIATED(dataPoints)) CALL FlagError("Data points is not associated.",err,error,*999)

    IF(dataPoints%dataPointsFinished) THEN
      localError="Data points number "//TRIM(NumberToVString(dataPoints%userNumber,"*",err,error))
      IF(ASSOCIATED(dataPoints%region)) localError=localError// &
        & " for region number "//TRIM(NumberToVString(dataPoints%region%userNumber,"*",err,error))
      IF(ASSOCIATED(dataPoints%interface)) localError=localError// &
        & " for interface number "//TRIM(NumberToVString(dataPoints%interface%userNumber,"*",err,error))
      localError=localError//" has already been finished."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("DataPoints_AssertNotFinished")
    RETURN
999 ERRORSEXITS("DataPoints_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE DataPoints_AssertNotFinished

  !
  !================================================================================================================================
  !

  !>Returns the coordinate system for data points accounting for regions and interfaces
  SUBROUTINE DataPoints_CoordinateSystemGet(dataPoints,coordinateSystem,err,error,*)

    !Argument variables
    TYPE(DataPointsType), POINTER :: dataPoints !<A pointer to the data points to get the coordinate system for
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem !<On exit, a pointer to the coordinate system for the data points. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DataPoints_CoordinateSystemGet",err,error,*998)

    IF(ASSOCIATED(coordinateSystem)) CALL FlagError("Coordinate system is already associated.",ERR,ERROR,*998)
    CALL DataPoints_AssertIsFinished(dataPoints,err,error,*999)

    IF(ASSOCIATED(dataPoints%region)) THEN
      IF(.NOT.dataPoints%region%regionFinished) CALL FlagError("Data points region has not been finished.",err,error,*999)
      coordinateSystem=>dataPoints%region%coordinateSystem
      IF(.NOT.ASSOCIATED(coordinateSystem)) THEN
        localError="The coordinate system for region number "// &
          & TRIM(NumberToVString(dataPoints%region%userNumber,"*",err,error))//" is not associated."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ELSE IF(ASSOCIATED(dataPoints%INTERFACE)) THEN
      IF(.NOT.dataPoints%interface%interfaceFinished) CALL FlagError("Data points interface has not been finished.",err,error,*999)
      coordinateSystem=>dataPoints%interface%coordinateSystem
      IF(.NOT.ASSOCIATED(coordinateSystem)) THEN
        localError="The coordinate system for interface number "// &
          & TRIM(NumberToVString(dataPoints%INTERFACE%userNumber,"*",err,error))//" is not associated."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Data points is not associated with a region or an interface.",err,error,*999)
    ENDIF
   
    EXITS("DataPoints_CoordinateSystemGet")
    RETURN
999 NULLIFY(coordinateSystem)
998 ERRORSEXITS("DataPoints_CoordinateSystemGet",err,error)
    RETURN 1
   
  END SUBROUTINE DataPoints_CoordinateSystemGet

  !
  !================================================================================================================================
  !
  
  !>Finds an returns a pointer to data points identified by a user number in the data point sets. 
  SUBROUTINE DataPointSets_UserNumberFind(dataPointSets,userNumber,dataPoints,err,error,*)

    !Argument variables
    TYPE(DataPointSetsType), POINTER :: dataPointSets !<A pointer to the data point sets to find the user number for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number to find
    TYPE(DataPointsType), POINTER :: dataPoints !<On exit, the pointer to the data points with the specified user number. If no data points with the specified user number exists the pointer is returned NULL. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: setIdx
    TYPE(DataPointsType), POINTER :: setDataPoints
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DataPointSets_UserNumberFind",err,error,*998)

    IF(ASSOCIATED(dataPointSets)) THEN
      IF(ASSOCIATED(dataPoints)) THEN
        CALL FlagError("Data points is already associated.",err,error,*998)
      ELSE
        NULLIFY(dataPoints)
        IF(ALLOCATED(dataPointSets%dataPointSets)) THEN
          setIdx=1
          DO WHILE(setIdx<=SIZE(dataPointSets%dataPointSets,1))
            setDataPoints=>dataPointSets%dataPointSets(setIdx)%ptr
            IF(ASSOCIATED(setDataPoints)) THEN
              IF(setDataPoints%userNumber==userNumber) THEN
                dataPoints=>dataPointSets%dataPointSets(setIdx)%ptr
                EXIT
              ENDIF
            ELSE
              localError="The data point sets data points is not associated for set index "// &
                & TRIM(NumberToVString(setIdx,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ENDDO
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Data point sets is not associated.",err,error,*999)
    ENDIF
    
    EXITS("DataPointSets_UserNumberFind")
    RETURN
999 NULLIFY(dataPoints)
998 ERRORSEXITS("DataPointSets_UserNumberFind",err,error)    
    RETURN 1
   
  END SUBROUTINE DataPointSets_UserNumberFind
  
  !
  !================================================================================================================================
  !
  
END MODULE DataPointAccessRoutines
