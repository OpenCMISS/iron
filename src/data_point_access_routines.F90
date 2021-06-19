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

!>This module contains all data point access method routines.
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

  !>Gets the label for a data point identified by a given global number.
  INTERFACE DataPoints_DataLabelGet
    MODULE PROCEDURE DataPoints_DataLabelGetC
    MODULE PROCEDURE DataPoints_DataLabelGetVS
  END INTERFACE DataPoints_DataLabelGet

  !>Changes/sets the label for a data point identified by a given global number.
  INTERFACE DataPoints_DataLabelSet
    MODULE PROCEDURE DataPoints_DataLabelSetC
    MODULE PROCEDURE DataPoints_DataLabelSetVS
  END INTERFACE DataPoints_DataLabelSet

  PUBLIC DataPoints_AssertIsFinished,DataPoints_AssertNotFinished

  PUBLIC DataPoints_CoordinateSystemGet

  PUBLIC DataPoints_DataLabelGet,DataPoints_DataLabelSet

  PUBLIC DataPoints_DataPositionGet

  PUBLIC DataPoints_DataUserNumberGet

  PUBLIC DataPoints_DataWeightsGet

  PUBLIC DataPoints_DataProjectionIndexGet

  PUBLIC DataPoints_DataProjectionUserGet

  PUBLIC DataPoints_InterfaceGet

  PUBLIC DataPoints_IsInterfaceDataPoints

  PUBLIC DataPoints_IsRegionDataPoints

  PUBLIC DataPoints_NumberOfDataPointsGet

  PUBLIC DataPoints_NumberOfDimensionsGet

  PUBLIC DataPoints_RegionGet

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

  !>Gets the character label for a data point identified by a given global number. \see OpenCMISS::Iron::cmfe_DataPoints_LabelGet
  SUBROUTINE DataPoints_DataLabelGetC(dataPoints,globalNumber,label,err,error,*)

    !Argument variables
    TYPE(DataPointsType), POINTER :: dataPoints !<A pointer to the data points to get the label for
    INTEGER(INTG), INTENT(IN) :: globalNumber !<The global number to get the label for
    CHARACTER(LEN=*), INTENT(OUT) :: label !<On exit, the label of the specified global data point
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER :: cLength,vsLength
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("DataPoints_DataLabelGetC",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(dataPoints)) CALL FlagError("Data points is not associated.",err,error,*999)
    IF(globalNumber<1.OR.globalNumber>dataPoints%numberOfDataPoints) THEN
      localError="The specified data point global number of "//TRIM(NumberToVString(globalNumber,"*",err,error))// &
        & " is invalid for data points number "//TRIM(NumberToVString(dataPoints%userNumber,"*",err,error))
      IF(ASSOCIATED(dataPoints%region)) THEN
        localError=localError//" in region number "//TRIM(NumberToVString(dataPoints%region%userNumber,"*",err,error))
      ELSE IF(ASSOCIATED(dataPoints%INTERFACE)) THEN
        localError=localError//" in interface number "//TRIM(NumberToVString(dataPoints%interface%userNumber,"*",err,error))
        IF(ASSOCIATED(dataPoints%INTERFACE%parentRegion)) localError=localError// &
          & " of parent region number "//TRIM(NumberToVString(dataPoints%interface%parentRegion%userNumber,"*",err,error))        
      ENDIF
      localError=localError//". The data point global number should be >= 1 and <= "// &
        & TRIM(NumberToVString(dataPoints%numberOfDataPoints,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)      
    ENDIF
    IF(.NOT.ALLOCATED(dataPoints%dataPoints)) THEN
      localError="The data points array is not allocated for data points number "// &
        & TRIM(NumberToVString(dataPoints%userNumber,"*",err,error))
      IF(ASSOCIATED(dataPoints%region)) THEN
        localError=localError//" in region number "//TRIM(NumberToVString(dataPoints%region%userNumber,"*",err,error))
      ELSE IF(ASSOCIATED(dataPoints%INTERFACE)) THEN
        localError=localError//" in interface number "//TRIM(NumberToVString(dataPoints%interface%userNumber,"*",err,error))
        IF(ASSOCIATED(dataPoints%INTERFACE%parentRegion)) localError=localError// &
          & " of parent region number "//TRIM(NumberToVString(dataPoints%interface%parentRegion%userNumber,"*",err,error))        
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)      
    ENDIF
#endif
    
    cLength=LEN(label)
    vsLength=LEN_TRIM(dataPoints%dataPoints(globalNumber)%label)
    IF(cLength>vsLength) THEN
      label=CHAR(LEN_TRIM(dataPoints%dataPoints(globalNumber)%label))
    ELSE
      label=CHAR(dataPoints%dataPoints(globalNumber)%label,cLength)
    ENDIF
    
    EXITS("DataPoints_DataLabelGetC")
    RETURN
999 ERRORSEXITS("DataPoints_DataLabelGetC",err,error)    
    RETURN 1
   
  END SUBROUTINE DataPoints_DataLabelGetC
        
  !
  !================================================================================================================================
  !

  !>Gets the varying string label for a data point identified by a given global number. \see OpenCMISS::Iron::cmfe_DataPoints_LabelGet
  SUBROUTINE DataPoints_DataLabelGetVS(dataPoints,globalNumber,label,err,error,*)

    !Argument variables
    TYPE(DataPointsType), POINTER :: dataPoints !<A pointer to the data points to get the label for
    INTEGER(INTG), INTENT(IN) :: globalNumber !<The global number to get the label for
    TYPE(VARYING_STRING), INTENT(OUT) :: label !<On exit, the label of the specified global data point
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("DataPoints_DataLabelGetVS",err,error,*999)
    
#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(dataPoints)) CALL FlagError("Data points is not associated.",err,error,*999)
    IF(globalNumber<1.OR.globalNumber>dataPoints%numberOfDataPoints) THEN
      localError="The specified data point global number of "//TRIM(NumberToVString(globalNumber,"*",err,error))// &
        & " is invalid for data points number "//TRIM(NumberToVString(dataPoints%userNumber,"*",err,error))
      IF(ASSOCIATED(dataPoints%region)) THEN
        localError=localError//" in region number "//TRIM(NumberToVString(dataPoints%region%userNumber,"*",err,error))
      ELSE IF(ASSOCIATED(dataPoints%INTERFACE)) THEN
        localError=localError//" in interface number "//TRIM(NumberToVString(dataPoints%interface%userNumber,"*",err,error))
        IF(ASSOCIATED(dataPoints%INTERFACE%parentRegion)) localError=localError// &
          & " of parent region number "//TRIM(NumberToVString(dataPoints%interface%parentRegion%userNumber,"*",err,error))        
      ENDIF
      localError=localError//". The data point global number should be >= 1 and <= "// &
        & TRIM(NumberToVString(dataPoints%numberOfDataPoints,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)      
    ENDIF
    IF(.NOT.ALLOCATED(dataPoints%dataPoints)) THEN
      localError="The data points array is not allocated for data points number "// &
        & TRIM(NumberToVString(dataPoints%userNumber,"*",err,error))
      IF(ASSOCIATED(dataPoints%region)) THEN
        localError=localError//" in region number "//TRIM(NumberToVString(dataPoints%region%userNumber,"*",err,error))
      ELSE IF(ASSOCIATED(dataPoints%INTERFACE)) THEN
        localError=localError//" in interface number "//TRIM(NumberToVString(dataPoints%interface%userNumber,"*",err,error))
        IF(ASSOCIATED(dataPoints%INTERFACE%parentRegion)) localError=localError// &
          & " of parent region number "//TRIM(NumberToVString(dataPoints%interface%parentRegion%userNumber,"*",err,error))        
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)      
    ENDIF
#endif

    label=dataPoints%dataPoints(globalNumber)%label
     
    EXITS("DataPoints_DataLabelGetVS")
    RETURN
999 ERRORSEXITS("DataPoints_DataLabelGetVS",err,error)    
    RETURN 1
   
  END SUBROUTINE DataPoints_DataLabelGetVS

  !
  !================================================================================================================================
  !

  !>Changes/sets the character label for a data point identified by a given global number.
  SUBROUTINE DataPoints_DataLabelSetC(dataPoints,globalNumber,label,err,error,*)

    !Argument variables
    TYPE(DataPointsType), POINTER :: dataPoints !<A pointer to the data points to set the label for
    INTEGER(INTG), INTENT(IN) :: globalNumber !<The global number to set the label for
    CHARACTER(LEN=*), INTENT(IN) :: label !<The label to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
     
    ENTERS("DataPoints_DataLabelSetC",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(dataPoints)) CALL FlagError("Data points is not associated.",err,error,*999)
    IF(globalNumber<1.OR.globalNumber>dataPoints%numberOfDataPoints) THEN
      localError="The specified data point global number of "//TRIM(NumberToVString(globalNumber,"*",err,error))// &
        & " is invalid for data points number "//TRIM(NumberToVString(dataPoints%userNumber,"*",err,error))
      IF(ASSOCIATED(dataPoints%region)) THEN
        localError=localError//" in region number "//TRIM(NumberToVString(dataPoints%region%userNumber,"*",err,error))
      ELSE IF(ASSOCIATED(dataPoints%INTERFACE)) THEN
        localError=localError//" in interface number "//TRIM(NumberToVString(dataPoints%interface%userNumber,"*",err,error))
        IF(ASSOCIATED(dataPoints%INTERFACE%parentRegion)) localError=localError// &
          & " of parent region number "//TRIM(NumberToVString(dataPoints%interface%parentRegion%userNumber,"*",err,error))        
      ENDIF
      localError=localError//". The data point global number should be >= 1 and <= "// &
        & TRIM(NumberToVString(dataPoints%numberOfDataPoints,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)      
    ENDIF
    IF(.NOT.ALLOCATED(dataPoints%dataPoints)) THEN
      localError="The data points array is not allocated for data points number "// &
        & TRIM(NumberToVString(dataPoints%userNumber,"*",err,error))
      IF(ASSOCIATED(dataPoints%region)) THEN
        localError=localError//" in region number "//TRIM(NumberToVString(dataPoints%region%userNumber,"*",err,error))
      ELSE IF(ASSOCIATED(dataPoints%INTERFACE)) THEN
        localError=localError//" in interface number "//TRIM(NumberToVString(dataPoints%interface%userNumber,"*",err,error))
        IF(ASSOCIATED(dataPoints%INTERFACE%parentRegion)) localError=localError// &
          & " of parent region number "//TRIM(NumberToVString(dataPoints%interface%parentRegion%userNumber,"*",err,error))        
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)      
    ENDIF
#endif

    dataPoints%dataPoints(globalNumber)%label=label
     
    EXITS("DataPoints_DataLabelSetC")
    RETURN
999 ERRORSEXITS("DataPoints_DataLabelSetC",err,error)    
    RETURN 1
   
  END SUBROUTINE DataPoints_DataLabelSetC    
  
  !
  !================================================================================================================================
  !


  !>Changes/sets the varying string label for a data point identified by a given global number.
  SUBROUTINE DataPoints_DataLabelSetVS(dataPoints,globalNumber,label,err,error,*)

    !Argument variables
    TYPE(DataPointsType), POINTER :: dataPoints !<A pointer to the data points to set the label for
    INTEGER(INTG), INTENT(IN) :: globalNumber !<The global number to set the label for
    TYPE(VARYING_STRING), INTENT(IN) :: label !<The label to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
     
    ENTERS("DataPoints_DataLabelSetVS",err,error,*999)
    
#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(dataPoints)) CALL FlagError("Data points is not associated.",err,error,*999)
    IF(globalNumber<1.OR.globalNumber>dataPoints%numberOfDataPoints) THEN
      localError="The specified data point global number of "//TRIM(NumberToVString(globalNumber,"*",err,error))// &
        & " is invalid for data points number "//TRIM(NumberToVString(dataPoints%userNumber,"*",err,error))
      IF(ASSOCIATED(dataPoints%region)) THEN
        localError=localError//" in region number "//TRIM(NumberToVString(dataPoints%region%userNumber,"*",err,error))
      ELSE IF(ASSOCIATED(dataPoints%INTERFACE)) THEN
        localError=localError//" in interface number "//TRIM(NumberToVString(dataPoints%interface%userNumber,"*",err,error))
        IF(ASSOCIATED(dataPoints%INTERFACE%parentRegion)) localError=localError// &
          & " of parent region number "//TRIM(NumberToVString(dataPoints%interface%parentRegion%userNumber,"*",err,error))        
      ENDIF
      localError=localError//". The data point global number should be >= 1 and <= "// &
        & TRIM(NumberToVString(dataPoints%numberOfDataPoints,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)      
    ENDIF
    IF(.NOT.ALLOCATED(dataPoints%dataPoints)) THEN
      localError="The data points array is not allocated for data points number "// &
        & TRIM(NumberToVString(dataPoints%userNumber,"*",err,error))
      IF(ASSOCIATED(dataPoints%region)) THEN
        localError=localError//" in region number "//TRIM(NumberToVString(dataPoints%region%userNumber,"*",err,error))
      ELSE IF(ASSOCIATED(dataPoints%INTERFACE)) THEN
        localError=localError//" in interface number "//TRIM(NumberToVString(dataPoints%interface%userNumber,"*",err,error))
        IF(ASSOCIATED(dataPoints%INTERFACE%parentRegion)) localError=localError// &
          & " of parent region number "//TRIM(NumberToVString(dataPoints%interface%parentRegion%userNumber,"*",err,error))        
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)      
    ENDIF
#endif

    dataPoints%dataPoints(globalNumber)%label=label
    
    EXITS("DataPoints_DataLabelSetVS")
    RETURN
999 ERRORSEXITS("DataPoints_DataLabelSetVS",err,error)    
    RETURN 1
   
  END SUBROUTINE DataPoints_DataLabelSetVS
        
  !
  !================================================================================================================================
  !
  
  !>Gets the position for a data point identified by a given global number.
  SUBROUTINE DataPoints_DataPositionGet(dataPoints,globalNumber,position,err,error,*)

    !Argument variables
    TYPE(DataPointsType), POINTER :: dataPoints !<A pointer to the data points to set the number for
    INTEGER(INTG), INTENT(IN) :: globalNumber !<The global number to get the values for
    REAL(DP), INTENT(OUT) :: position(:) !<position(coordinateIdx). On exit, the position of the specified data point
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("DataPoints_DataPositionGet",err,error,*999)

    CALL DataPoints_AssertIsFinished(dataPoints,err,error,*999)
#ifdef WITH_PRECHECKS
    IF(globalNumber<1.OR.globalNumber>dataPoints%numberOfDataPoints) THEN
      localError="The specified data point global number of "//TRIM(NumberToVString(globalNumber,"*",err,error))// &
        & " is invalid for data points number "//TRIM(NumberToVString(dataPoints%userNumber,"*",err,error))
      IF(ASSOCIATED(dataPoints%region)) THEN
        localError=localError//" in region number "//TRIM(NumberToVString(dataPoints%region%userNumber,"*",err,error))
      ELSE IF(ASSOCIATED(dataPoints%INTERFACE)) THEN
        localError=localError//" in interface number "//TRIM(NumberToVString(dataPoints%interface%userNumber,"*",err,error))
        IF(ASSOCIATED(dataPoints%INTERFACE%parentRegion)) localError=localError// &
          & " of parent region number "//TRIM(NumberToVString(dataPoints%interface%parentRegion%userNumber,"*",err,error))        
      ENDIF
      localError=localError//". The data point global number should be >= 1 and <= "// &
        & TRIM(NumberToVString(dataPoints%numberOfDataPoints,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)      
    ENDIF
    IF(.NOT.ALLOCATED(dataPoints%dataPoints)) THEN
      localError="The data points array is not allocated for data points number "// &
        & TRIM(NumberToVString(dataPoints%userNumber,"*",err,error))
      IF(ASSOCIATED(dataPoints%region)) THEN
        localError=localError//" in region number "//TRIM(NumberToVString(dataPoints%region%userNumber,"*",err,error))
      ELSE IF(ASSOCIATED(dataPoints%INTERFACE)) THEN
        localError=localError//" in interface number "//TRIM(NumberToVString(dataPoints%interface%userNumber,"*",err,error))
        IF(ASSOCIATED(dataPoints%INTERFACE%parentRegion)) localError=localError// &
          & " of parent region number "//TRIM(NumberToVString(dataPoints%interface%parentRegion%userNumber,"*",err,error))        
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)      
    ENDIF
    IF(.NOT.ALLOCATED(dataPoints%dataPoints(globalNumber)%position)) THEN
      localError="The position array is not allocated for data point global number "// &
        & TRIM(NumberToVString(globalNumber,"*",err,error))//" of data points number "// &
        & TRIM(NumberToVString(dataPoints%userNumber,"*",err,error))
      IF(ASSOCIATED(dataPoints%region)) THEN
        localError=localError//" in region number "//TRIM(NumberToVString(dataPoints%region%userNumber,"*",err,error))
      ELSE IF(ASSOCIATED(dataPoints%INTERFACE)) THEN
        localError=localError//" in interface number "//TRIM(NumberToVString(dataPoints%interface%userNumber,"*",err,error))
        IF(ASSOCIATED(dataPoints%INTERFACE%parentRegion)) localError=localError// &
          & " of parent region number "//TRIM(NumberToVString(dataPoints%interface%parentRegion%userNumber,"*",err,error))        
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)      
    ENDIF
    IF(SIZE(position,1)<dataPoints%numberOfDimensions) THEN
      localError="The size of the specified position array of "//TRIM(NumberToVString(SIZE(position,1),"*",err,error))// &
        & " is too small. The array size needs to be >= "//TRIM(NumberToVString(dataPoints%numberOfDimensions,"*",err,error))// &
        & " for data points number "//TRIM(NumberToVString(dataPoints%userNumber,"*",err,error))
      IF(ASSOCIATED(dataPoints%region)) THEN
        localError=localError//" in region number "//TRIM(NumberToVString(dataPoints%region%userNumber,"*",err,error))
      ELSE IF(ASSOCIATED(dataPoints%INTERFACE)) THEN
        localError=localError//" in interface number "//TRIM(NumberToVString(dataPoints%interface%userNumber,"*",err,error))
        IF(ASSOCIATED(dataPoints%INTERFACE%parentRegion)) localError=localError// &
          & " of parent region number "//TRIM(NumberToVString(dataPoints%interface%parentRegion%userNumber,"*",err,error))        
      ENDIF
      localError=localError//"."         
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    position(1:dataPoints%numberOfDimensions)=dataPoints%dataPoints(globalNumber)%position(1:dataPoints%numberOfDimensions)
   
    EXITS("DataPoints_DataPositionGet")
    RETURN
999 ERRORSEXITS("DataPoints_DataPositionGet",err,error)    
    RETURN 1

  END SUBROUTINE DataPoints_DataPositionGet

  !
  !================================================================================================================================
  !  

  !>Gets the user number for a data point identified by a given global number. 
  SUBROUTINE DataPoints_DataUserNumberGet(dataPoints,globalNumber,userNumber,err,error,*)

    !Argument variables
    TYPE(DataPointsType), POINTER :: dataPoints !<A pointer to the data points to get the number for
    INTEGER(INTG), INTENT(IN) :: globalNumber !<The global number to get the user number for
    INTEGER(INTG), INTENT(OUT) :: userNumber !<On exit, the user number of the specified global data point
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("DataPoints_DataUserNumberGet",err,error,*999)

    CALL DataPoints_AssertIsFinished(dataPoints,err,error,*999)
#ifdef WITH_PRECHECKS
    IF(globalNumber<1.OR.globalNumber>dataPoints%numberOfDataPoints) THEN
      localError="The specified data point global number of "//TRIM(NumberToVString(globalNumber,"*",err,error))// &
        & " is invalid for data points number "//TRIM(NumberToVString(dataPoints%userNumber,"*",err,error))
      IF(ASSOCIATED(dataPoints%region)) THEN
        localError=localError//" in region number "//TRIM(NumberToVString(dataPoints%region%userNumber,"*",err,error))
      ELSE IF(ASSOCIATED(dataPoints%INTERFACE)) THEN
        localError=localError//" in interface number "//TRIM(NumberToVString(dataPoints%interface%userNumber,"*",err,error))
        IF(ASSOCIATED(dataPoints%INTERFACE%parentRegion)) localError=localError// &
          & " of parent region number "//TRIM(NumberToVString(dataPoints%interface%parentRegion%userNumber,"*",err,error))        
      ENDIF
      localError=localError//". The data point global number should be >= 1 and <= "// &
        & TRIM(NumberToVString(dataPoints%numberOfDataPoints,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)      
    ENDIF
    IF(.NOT.ALLOCATED(dataPoints%dataPoints)) THEN
      localError="The data points array is not allocated for data points number "// &
        & TRIM(NumberToVString(dataPoints%userNumber,"*",err,error))
      IF(ASSOCIATED(dataPoints%region)) THEN
        localError=localError//" in region number "//TRIM(NumberToVString(dataPoints%region%userNumber,"*",err,error))
      ELSE IF(ASSOCIATED(dataPoints%INTERFACE)) THEN
        localError=localError//" in interface number "//TRIM(NumberToVString(dataPoints%interface%userNumber,"*",err,error))
        IF(ASSOCIATED(dataPoints%INTERFACE%parentRegion)) localError=localError// &
          & " of parent region number "//TRIM(NumberToVString(dataPoints%interface%parentRegion%userNumber,"*",err,error))        
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)      
    ENDIF
#endif
    
    userNumber=dataPoints%dataPoints(globalNumber)%userNumber
    
    EXITS("DataPoints_DataUserNumberGet")
    RETURN
999 ERRORSEXITS("DataPoints_DataUserNumberGet",err,error)    
    RETURN 1
   
  END SUBROUTINE DataPoints_DataUserNumberGet
        
  !
  !================================================================================================================================
  !
  
  !>Gets the weights for a data point identified by a given global number.
  SUBROUTINE DataPoints_DataWeightsGet(dataPoints,globalNumber,weights,err,error,*)

    !Argument variables
    TYPE(DataPointsType), POINTER :: dataPoints !<A pointer to the data points to set the weights for
    INTEGER(INTG), INTENT(IN) :: globalNumber !<The global number to get the weights for
    REAL(DP), INTENT(OUT) :: weights(:) !<weights(coordinateIdx). On exit, the weights of the specified data point
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("DataPoints_DataWeightsGet",err,error,*999)

    CALL DataPoints_AssertIsFinished(dataPoints,err,error,*999)
#ifdef WITH_PRECHECKS
    IF(globalNumber<1.OR.globalNumber>dataPoints%numberOfDataPoints) THEN
      localError="The specified data point global number of "//TRIM(NumberToVString(globalNumber,"*",err,error))// &
        & " is invalid for data points number "//TRIM(NumberToVString(dataPoints%userNumber,"*",err,error))
      IF(ASSOCIATED(dataPoints%region)) THEN
        localError=localError//" in region number "//TRIM(NumberToVString(dataPoints%region%userNumber,"*",err,error))
      ELSE IF(ASSOCIATED(dataPoints%INTERFACE)) THEN
        localError=localError//" in interface number "//TRIM(NumberToVString(dataPoints%interface%userNumber,"*",err,error))
        IF(ASSOCIATED(dataPoints%INTERFACE%parentRegion)) localError=localError// &
          & " of parent region number "//TRIM(NumberToVString(dataPoints%interface%parentRegion%userNumber,"*",err,error))        
      ENDIF
      localError=localError//". The data point global number should be >= 1 and <= "// &
        & TRIM(NumberToVString(dataPoints%numberOfDataPoints,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)      
    ENDIF
    IF(.NOT.ALLOCATED(dataPoints%dataPoints)) THEN
      localError="The data points array is not allocated for data points number "// &
        & TRIM(NumberToVString(dataPoints%userNumber,"*",err,error))
      IF(ASSOCIATED(dataPoints%region)) THEN
        localError=localError//" in region number "//TRIM(NumberToVString(dataPoints%region%userNumber,"*",err,error))
      ELSE IF(ASSOCIATED(dataPoints%INTERFACE)) THEN
        localError=localError//" in interface number "//TRIM(NumberToVString(dataPoints%interface%userNumber,"*",err,error))
        IF(ASSOCIATED(dataPoints%INTERFACE%parentRegion)) localError=localError// &
          & " of parent region number "//TRIM(NumberToVString(dataPoints%interface%parentRegion%userNumber,"*",err,error))        
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)      
    ENDIF
    IF(.NOT.ALLOCATED(dataPoints%dataPoints(globalNumber)%weights)) THEN
      localError="The weights array is not allocated for data point global number "// &
        & TRIM(NumberToVString(globalNumber,"*",err,error))//" of data points number "// &
        & TRIM(NumberToVString(dataPoints%userNumber,"*",err,error))
      IF(ASSOCIATED(dataPoints%region)) THEN
        localError=localError//" in region number "//TRIM(NumberToVString(dataPoints%region%userNumber,"*",err,error))
      ELSE IF(ASSOCIATED(dataPoints%INTERFACE)) THEN
        localError=localError//" in interface number "//TRIM(NumberToVString(dataPoints%interface%userNumber,"*",err,error))
        IF(ASSOCIATED(dataPoints%INTERFACE%parentRegion)) localError=localError// &
          & " of parent region number "//TRIM(NumberToVString(dataPoints%interface%parentRegion%userNumber,"*",err,error))        
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)      
    ENDIF
    IF(SIZE(weights,1)<dataPoints%numberOfDimensions) THEN
      localError="The size of the specified weights array of "//TRIM(NumberToVString(SIZE(weights,1),"*",err,error))// &
        & " is too small. The array size needs to be >= "//TRIM(NumberToVString(dataPoints%numberOfDimensions,"*",err,error))// &
        & " for data points number "//TRIM(NumberToVString(dataPoints%userNumber,"*",err,error))
      IF(ASSOCIATED(dataPoints%region)) THEN
        localError=localError//" in region number "//TRIM(NumberToVString(dataPoints%region%userNumber,"*",err,error))
      ELSE IF(ASSOCIATED(dataPoints%INTERFACE)) THEN
        localError=localError//" in interface number "//TRIM(NumberToVString(dataPoints%interface%userNumber,"*",err,error))
        IF(ASSOCIATED(dataPoints%INTERFACE%parentRegion)) localError=localError// &
          & " of parent region number "//TRIM(NumberToVString(dataPoints%interface%parentRegion%userNumber,"*",err,error))        
      ENDIF
      localError=localError//"."         
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    weights(1:dataPoints%numberOfDimensions)=dataPoints%dataPoints(globalNumber)%weights(1:dataPoints%numberOfDimensions)
   
    EXITS("DataPoints_DataWeightsGet")
    RETURN
999 ERRORSEXITS("DataPoints_DataWeightsGet",err,error)    
    RETURN 1

  END SUBROUTINE DataPoints_DataWeightsGet

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
#endif    
    
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

  !>Returns the interface for data points
  SUBROUTINE DataPoints_InterfaceGet(dataPoints,interface,err,error,*)

    !Argument variables
    TYPE(DataPointsType), POINTER :: dataPoints !<A pointer to the data points to get the interface for
    TYPE(InterfaceType), POINTER :: interface !<On return, the data points interface. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DataPoints_InterfaceGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(ASSOCIATED(interface)) CALL FlagError("Interface is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(dataPoints)) CALL FlagError("Data points is not associated.",err,error,*999)
#endif    

    INTERFACE=>dataPoints%INTERFACE

#ifdef WITH_POSTCHECKS      
    IF(.NOT.ASSOCIATED(INTERFACE)) THEN
      localError="The interface for data points number "// &
        & TRIM(NumberToVString(dataPoints%userNumber,"*",err,error))//" is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
 
    EXITS("DataPoints_InterfaceGet")
    RETURN
999 NULLIFY(interface)
998 ERRORSEXITS("DataPoints_InterfaceGet",err,error)
    RETURN 1
    
  END SUBROUTINE DataPoints_InterfaceGet

  !
  !================================================================================================================================
  !

  !>Determines if the given data points are interface data points or not. 
  SUBROUTINE DataPoints_IsInterfaceDataPoints(dataPoints,interfaceDataPoints,err,error,*)

    !Argument variables
    TYPE(DataPointsType), POINTER :: dataPoints !<A pointer to the data points to determine if they are interface data points or not.
    LOGICAL :: interfaceDataPoints !<On exit, .TRUE. if the given data points are interface data points, .FALSE. if not. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("DataPoints_IsInterfaceDataPoints",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(dataPoints)) CALL FlagError("Data points is not associated.",err,error,*999)
#endif    

    interfaceDataPoints = ASSOCIATED(dataPoints%interface)
    
    EXITS("DataPoints_IsInterfaceDataPoints")
    RETURN
999 ERRORSEXITS("DataPoints_IsInterfaceDataPoints",err,error)
    RETURN 1
    
  END SUBROUTINE DataPoints_IsInterfaceDataPoints

  !
  !================================================================================================================================
  !

  !>Determines if the given data points are region data points or not. 
  SUBROUTINE DataPoints_IsRegionDataPoints(dataPoints,regionDataPoints,err,error,*)

    !Argument variables
    TYPE(DataPointsType), POINTER :: dataPoints !<A pointer to the data points to determine if they are region data points or not.
    LOGICAL :: regionDataPoints !<On exit, .TRUE. if the given data point are in a region, .FALSE. if not. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("DataPoints_IsRegionDataPoints",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(dataPoints)) CALL FlagError("Data points is not associated.",err,error,*999)
#endif    

    regionDataPoints = ASSOCIATED(dataPoints%region)
    
    EXITS("DataPoints_IsRegionDataPoints")
    RETURN
999 ERRORSEXITS("DataPoints_IsRegionDataPoints",err,error)
    RETURN 1
    
  END SUBROUTINE DataPoints_IsRegionDataPoints

  !
  !================================================================================================================================
  !

  !>Returns the number of data points. \see OpenCMISS::Iron::cmfe_DataPoints_NumberOfDataPointsGet
  SUBROUTINE DataPoints_NumberOfDataPointsGet(dataPoints,numberOfDataPoints,err,error,*)

    !Argument variables
    TYPE(DataPointsType), POINTER :: dataPoints !<A pointer to the data points to get the number of data points for
    INTEGER(INTG), INTENT(OUT) :: numberOfDataPoints !<On return, the number of data points
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("DataPoints_NumberOfDataPointsGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(dataPoints)) CALL FlagError("Data points is not associated.",err,error,*999)
#endif
    
    numberOfDataPoints=dataPoints%numberOfDataPoints
    
    EXITS("DataPoints_NumberOfDataPointsGet")
    RETURN
999 ERRORSEXITS("DataPoints_NumberOfDataPointsGet",err,error)    
    RETURN 1
   
  END SUBROUTINE DataPoints_NumberOfDataPointsGet

  !
  !================================================================================================================================
  !

  !>Returns the number of data point dimensions. 
  SUBROUTINE DataPoints_NumberOfDimensionsGet(dataPoints,numberOfDimensions,err,error,*)

    !Argument variables
    TYPE(DataPointsType), POINTER :: dataPoints !<A pointer to the data points to get the number of dimensions for
    INTEGER(INTG), INTENT(OUT) :: numberOfDimensions !<On return, the number of data point dimensions
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("DataPoints_NumberOfDimensionsGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(dataPoints)) CALL FlagError("Data points is not associated.",err,error,*999)
#endif
    
    numberOfDimensions=dataPoints%numberOfDimensions
    
    EXITS("DataPoints_NumberOfDimensionsGet")
    RETURN
999 ERRORSEXITS("DataPoints_NumberOfDimensionsGet",err,error)    
    RETURN 1
   
  END SUBROUTINE DataPoints_NumberOfDimensionsGet

  !
  !================================================================================================================================
  !

  !>Returns the region for a data points accounting for regions and interfaces
  SUBROUTINE DataPoints_RegionGet(dataPoints,region,err,error,*)

    !Argument variables
    TYPE(DataPointsType), POINTER :: dataPoints !<A pointer to the data points to get the region for
    TYPE(RegionType), POINTER :: region !<On return, the data points region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(InterfaceType), POINTER :: interface
    TYPE(VARYING_STRING) :: localError

    ENTERS("DataPoints_RegionGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(ASSOCIATED(region)) CALL FlagError("Region is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(dataPoints)) CALL FlagError("Field is not associated.",err,error,*999)
#endif    
        
    NULLIFY(region)
    NULLIFY(interface)
    region=>dataPoints%region
    IF(.NOT.ASSOCIATED(region)) THEN          
      INTERFACE=>dataPoints%INTERFACE
      IF(ASSOCIATED(INTERFACE)) THEN
        IF(ASSOCIATED(interface%parentRegion)) THEN
          region=>interface%parentRegion     
        ELSE
          localError="The parent region is not associated for data points number "// &
            & TRIM(NumberToVString(dataPoints%userNumber,"*",err,error))//" of interface number "// &
            & TRIM(NumberToVString(interface%userNumber,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        localError="A region or interface is not associated for data points number "// &
          & TRIM(NumberToVString(dataPoints%userNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDIF
    
    EXITS("DataPoints_RegionGet")
    RETURN
999 NULLIFY(region)    
998 ERRORSEXITS("DataPoints_RegionGet",err,error)
    RETURN 1
    
  END SUBROUTINE DataPoints_RegionGet

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
