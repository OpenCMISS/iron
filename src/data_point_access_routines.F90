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
  USE DataProjectionAccessRoutines
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

  PUBLIC DataPoints_DataProjectionIndexGet

  PUBLIC DataPoints_DataProjectionUserGet

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

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(dataPoints)) CALL FlagError("Data points is not associated.",err,error,*999)
#endif    

    IF(.NOT.dataPoints%dataPointsFinished) THEN
      localError="Data points number "//TRIM(NumberToVString(dataPoints%userNumber,"*",err,error))
      IF(ASSOCIATED(dataPoints%region)) THEN
        localError=localError//" for region number "//TRIM(NumberToVString(dataPoints%region%userNumber,"*",err,error))
      ELSE IF(ASSOCIATED(dataPoints%INTERFACE)) THEN
        localError=localError//" for interface number "//TRIM(NumberToVString(dataPoints%INTERFACE%userNumber,"*",err,error))
        IF(ASSOCIATED(dataPoints%INTERFACE%parentRegion)) localError=localError// &
          & " of parent region number "//TRIM(NumberToVString(dataPoints%INTERFACE%parentRegion%userNumber,"*",err,error))        
      ENDIF
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

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(dataPoints)) CALL FlagError("Data points is not associated.",err,error,*999)
#endif    

    IF(dataPoints%dataPointsFinished) THEN
      localError="Data points number "//TRIM(NumberToVString(dataPoints%userNumber,"*",err,error))
      IF(ASSOCIATED(dataPoints%region)) THEN
        localError=localError//" for region number "//TRIM(NumberToVString(dataPoints%region%userNumber,"*",err,error))
      ELSE IF(ASSOCIATED(dataPoints%INTERFACE)) THEN
        localError=localError//" for interface number "//TRIM(NumberToVString(dataPoints%INTERFACE%userNumber,"*",err,error))
        IF(ASSOCIATED(dataPoints%INTERFACE%parentRegion)) localError=localError// &
          & " of parent region number "//TRIM(NumberToVString(dataPoints%INTERFACE%parentRegion%userNumber,"*",err,error))        
      ENDIF
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
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("DataPoints_CoordinateSystemGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(coordinateSystem)) CALL FlagError("Coordinate system is already associated.",err,error,*998)
#endif    
    CALL DataPoints_AssertIsFinished(dataPoints,err,error,*999)
    

    IF(ASSOCIATED(dataPoints%region)) THEN      
      IF(.NOT.dataPoints%region%regionFinished) CALL FlagError("Data points region has not been finished.",err,error,*999)
      coordinateSystem=>dataPoints%region%coordinateSystem
#ifdef WITH_POSTCHECKS      
      IF(.NOT.ASSOCIATED(coordinateSystem)) THEN
        localError="The coordinate system for region number "// &
          & TRIM(NumberToVString(dataPoints%region%userNumber,"*",err,error))//" is not associated."
        CALL FlagError(localError,err,error,*999)
      ENDIF
#endif      
    ELSE IF(ASSOCIATED(dataPoints%INTERFACE)) THEN
      IF(.NOT.dataPoints%interface%interfaceFinished) CALL FlagError("Data points interface has not been finished.",err,error,*999)
      coordinateSystem=>dataPoints%INTERFACE%coordinateSystem
#ifdef WITH_POSTCHECKS      
      IF(.NOT.ASSOCIATED(coordinateSystem)) THEN
        localError="The coordinate system for interface number "// &
          & TRIM(NumberToVString(dataPoints%INTERFACE%userNumber,"*",err,error))
        IF(ASSOCIATED(dataPoints%INTERFACE%parentRegion)) localError=localError// &
          & " of parent region number "//TRIM(NumberToVString(dataPoints%INTERFACE%parentRegion%userNumber,"*",err,error))
        localError=localError//" is not associated."
        CALL FlagError(localError,err,error,*999)
      ENDIF
#endif      
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

  !>Returns a data projection for data points
  SUBROUTINE DataPoints_DataProjectionIndexGet(dataPoints,dataProjectionIndex,dataProjection,err,error,*)

    !Argument variables
    TYPE(DataPointsType), POINTER :: dataPoints !<A pointer to the data points to get the coordinate system for
    INTEGER(INTG), INTENT(IN) :: dataProjectionIndex !<The index of the dat projection to get
    TYPE(DataProjectionType), POINTER :: dataProjection !<On exit, a pointer to the specified data projeciton for the data points. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("DataPoints_DataProjectionIndexGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(dataProjection)) CALL FlagError("Data projection is already associated.",err,error,*998)
    CALL DataPoints_AssertIsFinished(dataPoints,err,error,*999)
    IF(.NOT.ASSOCIATED(dataPoints%dataProjections)) THEN
      localError="The data projections is not associated for data points number "// &
        & TRIM(NumberToVString(dataPoints%userNumber,"*",err,error))
      IF(ASSOCIATED(dataPoints%region)) THEN
        localError=localError//" of region number "//TRIM(NumberToVString(dataPoints%region%userNumber,"*",err,error))
      ELSE IF(ASSOCIATED(dataPoints%INTERFACE)) THEN
        localError=localError//" of interface number "//TRIM(NumberToVString(dataPoints%INTERFACE%userNumber,"*",err,error))
        IF(ASSOCIATED(dataPoints%INTERFACE%parentRegion)) localError=localError// &
          & " of parent region number "//TRIM(NumberToVString(dataPoints%interface%parentRegion%userNumber,"*",err,error))        
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(dataProjectionIndex<=0.OR.dataProjectionIndex>dataPoints%dataProjections%numberOfDataProjections) THEN
      localError="The specified data projection index of "//TRIM(NumberToVString(dataProjectionIndex,"*",err,error))// &
        & " is invalid for data points number "//TRIM(NumberToVString(dataPoints%userNumber,"*",err,error))
      IF(ASSOCIATED(dataPoints%region)) THEN
        localError=localError//" of region number "//TRIM(NumberToVString(dataPoints%region%userNumber,"*",err,error))
      ELSE IF(ASSOCIATED(dataPoints%INTERFACE)) THEN
        localError=localError//" of interface number "//TRIM(NumberToVString(dataPoints%INTERFACE%userNumber,"*",err,error))
        IF(ASSOCIATED(dataPoints%INTERFACE%parentRegion)) localError=localError// &
          & " of parent region number "//TRIM(NumberToVString(dataPoints%interface%parentRegion%userNumber,"*",err,error))        
      ENDIF
      localError=localError//". The data projection index should be >=1 and <= "// &
        & TRIM(NumberToVString(dataPoints%dataProjections%numberOfDataProjections,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(dataPoints%dataProjections%dataProjections)) THEN
      localError="The data projections data projections have not been allocated for data points number "// &
        & TRIM(NumberToVString(dataPoints%userNumber,"*",err,error))
      IF(ASSOCIATED(dataPoints%region)) THEN
        localError=localError//" of region number "//TRIM(NumberToVString(dataPoints%region%userNumber,"*",err,error))
      ELSE IF(ASSOCIATED(dataPoints%INTERFACE)) THEN
        localError=localError//" of interface number "//TRIM(NumberToVString(dataPoints%INTERFACE%userNumber,"*",err,error))
        IF(ASSOCIATED(dataPoints%INTERFACE%parentRegion)) localError=localError// &
          & " of parent region number "//TRIM(NumberToVString(dataPoints%INTERFACE%parentRegion%userNumber,"*",err,error))        
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    dataProjection=>dataPoints%dataProjections%dataProjections(dataProjectionIndex)%ptr

#ifdef WITH_POSTCHECKS    
    IF(ASSOCIATED(dataProjection)) THEN
      localError="The data projection is not associated for data projection index "// &
        & TRIM(NumberToVString(dataProjectionIndex,"*",err,error))//" of data points number "// &
        & TRIM(NumberToVString(dataPoints%userNumber,"*",err,error))
      IF(ASSOCIATED(dataPoints%region)) THEN
        localError=localError//" of region number "//TRIM(NumberToVString(dataPoints%region%userNumber,"*",err,error))
      ELSE IF(ASSOCIATED(dataPoints%INTERFACE)) THEN
        localError=localError//" of interface number "//TRIM(NumberToVString(dataPoints%INTERFACE%userNumber,"*",err,error))
        IF(ASSOCIATED(dataPoints%INTERFACE%parentRegion)) localError=localError// &
          & " of parent region number "//TRIM(NumberToVString(dataPoints%INTERFACE%parentRegion%userNumber,"*",err,error))        
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
   
    EXITS("DataPoints_DataProjectionIndexGet")
    RETURN
999 NULLIFY(dataProjection)
998 ERRORSEXITS("DataPoints_DataProjectionIndexGet",err,error)
    RETURN 1
   
  END SUBROUTINE DataPoints_DataProjectionIndexGet

  !
  !================================================================================================================================
  !  

  !>Gets the data projection identified by a given user number. 
  SUBROUTINE DataPoints_DataProjectionUserGet(dataPoints,userNumber,dataProjection,err,error,*)

    !Argument variables
    TYPE(DataPointsType), POINTER :: dataPoints !<A pointer to the data points to get the data projection for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number to get the data projection for
    TYPE(DataProjectionType), POINTER :: dataProjection !<On exit, a pointer to the data projection for the data points. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DataPoints_DataProjectionUserGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(dataProjection)) CALL FlagError("Data projection is already associated.",err,error,*998)
#endif    
    CALL DataPoints_AssertIsFinished(dataPoints,err,error,*999)
    
    CALL DataProjection_UserNumberFind(dataPoints,userNumber,dataProjection,err,error,*999)
#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(dataProjection)) THEN
      localError="A data projection with a user number of "//TRIM(NumberToVString(userNumber,"*",err,error))// &
        & " does not exist for data points number "//TRIM(NumberToVString(dataPoints%userNumber,"*",err,error))
      IF(ASSOCIATED(dataPoints%region)) THEN
        localError=localError//" of region number "//TRIM(NumberToVString(dataPoints%region%userNumber,"*",err,error))
      ELSE IF(ASSOCIATED(dataPoints%INTERFACE)) THEN
        localError=localError//" of interface number "//TRIM(NumberToVString(dataPoints%INTERFACE%userNumber,"*",err,error))
        IF(ASSOCIATED(dataPoints%INTERFACE%parentRegion)) localError=localError// &
          & " of parent region number "//TRIM(NumberToVString(dataPoints%INTERFACE%parentRegion%userNumber,"*",err,error))        
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    EXITS("DataPoints_DataProjectionUserGet")
    RETURN
999 NULLIFY(dataProjection)   
998 ERRORSEXITS("DataPoints_DataProjectionUserGet",err,error)    
    RETURN 1
   
  END SUBROUTINE DataPoints_DataProjectionUserGet

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
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("DataPointSets_UserNumberFind",err,error,*998)

#ifdef WITH_PRECHECKS
    IF(ASSOCIATED(dataPoints)) CALL FlagError("Data points is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(dataPointSets)) CALL FlagError("Data point sets is not associated.",err,error,*999)
#endif    
    
    NULLIFY(dataPoints)
    IF(ALLOCATED(dataPointSets%dataPointSets)) THEN
      setIdx=1
      DO WHILE(setIdx<=SIZE(dataPointSets%dataPointSets,1))
#ifdef WITH_PRECHECKS        
        IF(.NOT.ASSOCIATED(dataPointSets%dataPointSets(setIdx)%ptr)) THEN
          localError="The data point sets data points is not associated for set index "// &
            & TRIM(NumberToVString(setIdx,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
#endif        
        IF(dataPointSets%dataPointSets(setIdx)%ptr%userNumber==userNumber) THEN
          dataPoints=>dataPointSets%dataPointSets(setIdx)%ptr
          EXIT
        ENDIF
      ENDDO
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
