!> \file
!> \author Chris Bradley
!> \brief This module contains all data projection access method routines.
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

!> This module contains all data projection access method routines.
MODULE DataProjectionAccessRoutines
  
  USE BaseRoutines
  USE Kinds
  USE ISO_VARYING_STRING
  USE Strings
  USE Types

#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup DataProjectionRoutines_DataProjectionExitTags DataProjectionRoutines::DataProjectionExitTags
  !> \brief Datapoint projection exit tags
  !> \see DataProjectionRoutines,OPENCMISS_DataProjectionExitTags
  !>@{ 
  INTEGER(INTG), PARAMETER :: DATA_PROJECTION_CANCELLED=0 !<Data projection has been cancelled. \see DataProjectionRoutines,OPENCMISS_DataProjectionExitTags
  INTEGER(INTG), PARAMETER :: DATA_PROJECTION_EXIT_TAG_CONVERGED=1 !<Data projection exited due to it being converged \see DataProjectionRoutines,OPENCMISS_DataProjectionExitTags
  INTEGER(INTG), PARAMETER :: DATA_PROJECTION_EXIT_TAG_BOUNDS=2 !<Data projection exited due to it hitting the bound and continue to travel out of the element. \see DataProjectionRoutines,OPENCMISS_DataProjectionExitTags
  INTEGER(INTG), PARAMETER :: DATA_PROJECTION_EXIT_TAG_MAX_ITERATION=3 !<Data projection exited due to it attaining maximum number of iteration specified by user. \see DataProjectionRoutines,OPENCMISS_DataProjectionExitTags
  INTEGER(INTG), PARAMETER :: DATA_PROJECTION_EXIT_TAG_NO_ELEMENT=4 !<Data projection exited due to no local element found, this happens when none of the candidate elements are within this computational node, and before MPI communication with other nodes. \see DataProjectionRoutines,OPENCMISS_DataProjectionExitTags
  !>@}

  !Module types

  !Module variables

  !Interfaces

  PUBLIC DATA_PROJECTION_CANCELLED,DATA_PROJECTION_EXIT_TAG_CONVERGED,DATA_PROJECTION_EXIT_TAG_BOUNDS, &
    & DATA_PROJECTION_EXIT_TAG_MAX_ITERATION,DATA_PROJECTION_EXIT_TAG_NO_ELEMENT

  PUBLIC DataProjection_AssertIsFinished,DataProjection_AssertNotFinished

  PUBLIC DataProjection_AssertIsProjected,DataProjection_AssertNotProjected

  PUBLIC DataProjection_DataPointsGet

  PUBLIC DataProjection_DecompositionGet

  PUBLIC DataProjection_ProjectionFieldGet

  PUBLIC DataProjection_ResultMaximumErrorGet

  PUBLIC DataProjection_ResultMinimumErrorGet

  PUBLIC DataProjection_ResultRMSErrorGet

  PUBLIC DataProjection_UserNumberFind
  
CONTAINS

  !
  !=================================================================================================================================
  !

  !>Assert that a data projection has been finished
  SUBROUTINE DataProjection_AssertIsFinished(dataProjection,err,error,*)

    !Argument Variables
    TYPE(DataProjectionType), POINTER, INTENT(IN) :: dataProjection !<The data projection to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("DataProjection_AssertIsFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
#endif    

    IF(.NOT.dataProjection%dataProjectionFinished) THEN
      localError="Data projection number "//TRIM(NumberToVString(dataProjection%userNumber,"*",err,error))
      IF(ASSOCIATED(dataProjection%projectionField)) localError=localError// &
        & " for field number "//TRIM(NumberToVString(dataProjection%projectionField%userNumber,"*",err,error))
      localError=localError//" has not been finished."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("DataProjection_AssertIsFinished")
    RETURN
999 ERRORSEXITS("DataProjection_AssertIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE DataProjection_AssertIsFinished

  !
  !=================================================================================================================================
  !

  !>Assert that a data projection has not been finished
  SUBROUTINE DataProjection_AssertNotFinished(dataProjection,err,error,*)

    !Argument Variables
    TYPE(DataProjectionType), POINTER, INTENT(IN) :: dataProjection !<The data projection to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("DataProjection_AssertNotFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
#endif    

    IF(dataProjection%dataProjectionFinished) THEN
      localError="Data projection number "//TRIM(NumberToVString(dataProjection%userNumber,"*",err,error))
      IF(ASSOCIATED(dataProjection%projectionField)) localError=localError// &
        & " for field number "//TRIM(NumberToVString(dataProjection%projectionField%userNumber,"*",err,error))
      localError=localError//" has already been finished."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("DataProjection_AssertNotFinished")
    RETURN
999 ERRORSEXITS("DataProjection_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE DataProjection_AssertNotFinished

  !
  !=================================================================================================================================
  !

  !>Assert that a data projection has been projected
  SUBROUTINE DataProjection_AssertIsProjected(dataProjection,err,error,*)

    !Argument Variables
    TYPE(DataProjectionType), POINTER, INTENT(IN) :: dataProjection !<The data projection to assert the projected status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("DataProjection_AssertIsProjected",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
#endif    

    IF(.NOT.dataProjection%dataProjectionProjected) THEN
      localError="Data projection number "//TRIM(NumberToVString(dataProjection%userNumber,"*",err,error))// &
        & " has not been projected."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("DataProjection_AssertIsProjected")
    RETURN
999 ERRORSEXITS("DataProjection_AssertIsProjected",err,error)
    RETURN 1
    
  END SUBROUTINE DataProjection_AssertIsProjected

  !
  !=================================================================================================================================
  !

  !>Assert that a data projection has not been projected
  SUBROUTINE DataProjection_AssertNotProjected(dataProjection,err,error,*)

    !Argument Variables
    TYPE(DataProjectionType), POINTER, INTENT(IN) :: dataProjection !<The data projection to assert the projected status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("DataProjection_AssertNotProjected",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
#endif    

    IF(dataProjection%dataProjectionProjected) THEN
      localError="Data projection number "//TRIM(NumberToVString(dataProjection%userNumber,"*",err,error))// &
        & " has already been projected."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("DataProjection_AssertNotProjected")
    RETURN
999 ERRORSEXITS("DataProjection_AssertNotProjected",err,error)
    RETURN 1
    
  END SUBROUTINE DataProjection_AssertNotProjected

  !
  !================================================================================================================================
  !

  !>Gets the data points for a data projection.
  SUBROUTINE DataProjection_DataPointsGet(dataProjection,dataPoints,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to get the data points for
    TYPE(DataPointsType), POINTER :: dataPoints !<On exit, the data points of the data projection. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
   
    ENTERS("DataProjection_DataPointsGet",err,error,*999)

#ifdef WITH_PRECHECK    
    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(ASSOCIATED(dataPoints)) CALL FlagError("Data points is already associated.",err,error,*999)
#endif    
    
    dataPoints=>dataProjection%dataPoints

#ifdef WITH_POSTCHECK    
    IF(.NOT.ASSOCIATED(dataPoints)) CALL FlagError("Data projection data points is not associated.",err,error,*999)
#endif    
 
    EXITS("DataProjection_DataPointsGet")
    RETURN
999 ERRORSEXITS("DataProjection_DataPointsGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_DataPointsGet

  !
  !================================================================================================================================
  !

  !>Gets the decompositon for a data projection.
  SUBROUTINE DataProjection_DecompositionGet(dataProjection,decomposition,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to get the decomposition for
    TYPE(DecompositionType), POINTER :: decomposition !<On exit, the decomposition of the data projection. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
   
    ENTERS("DataProjection_DecompositionGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(ASSOCIATED(decomposition)) CALL FlagError("Decomposition is already associated.",err,error,*999)
#endif    
    
    decomposition=>dataProjection%decomposition

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Data projection decomposition is not associated.",err,error,*999)
#endif    
 
    EXITS("DataProjection_DecompositionGet")
    RETURN
999 ERRORSEXITS("DataProjection_DecompositionGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_DecompositionGet

  !
  !================================================================================================================================
  !

  !>Gets the projection field for a data projection.
  SUBROUTINE DataProjection_ProjectionFieldGet(dataProjection,projectionField,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to get the projectionField for
    TYPE(FieldType), POINTER :: projectionField !<On exit, the projection field of the data projection. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
   
    ENTERS("DataProjection_ProjectionFieldGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(ASSOCIATED(projectionField)) CALL FlagError("Projection field is already associated.",err,error,*999)
#endif    
    
    projectionField=>dataProjection%projectionField

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(projectionField)) CALL FlagError("Data projection projection field is not associated.",err,error,*999)
#endif    
 
    EXITS("DataProjection_ProjectionFieldGet")
    RETURN
999 ERRORSEXITS("DataProjection_ProjectionFieldGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_ProjectionFieldGet

  !
  !================================================================================================================================
  !

  !>Gets the projection maximum error for a data projection.
  SUBROUTINE DataProjection_ResultMaximumErrorGet(dataProjection,maximumDataPoint,maximumError,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to get the maximum error for
    INTEGER(INTG), INTENT(OUT) :: maximumDataPoint !<On exit, the data point number of the maximum error.
    REAL(DP), INTENT(OUT) :: maximumError !<On exit, the maximum error for the data projection.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
   
    ENTERS("DataProjection_ResultMaximumErrorGet",err,error,*999)

    CALL DataProjection_AssertIsFinished(dataProjection,err,error,*999)
    CALL DataProjection_AssertIsProjected(dataProjection,err,error,*999)
    
    maximumDataPoint=dataProjection%maximumErrorDataPoint
    maximumError=dataProjection%maximumError
 
    EXITS("DataProjection_ResultMaximumErrorGet")
    RETURN
999 ERRORSEXITS("DataProjection_ResultMaximumErrorGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_ResultMaximumErrorGet

  !
  !================================================================================================================================
  !

  !>Gets the projection minimum error for a data projection.
  SUBROUTINE DataProjection_ResultMinimumErrorGet(dataProjection,minimumDataPoint,minimumError,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to get the minimum error for
    INTEGER(INTG), INTENT(OUT) :: minimumDataPoint !<On exit, the data point number of the minimum error.
    REAL(DP), INTENT(OUT) :: minimumError !<On exit, the maximum error for the data projection.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
   
    ENTERS("DataProjection_ResultMinimumErrorGet",err,error,*999)

    CALL DataProjection_AssertIsFinished(dataProjection,err,error,*999)
    CALL DataProjection_AssertIsProjected(dataProjection,err,error,*999)
    
    minimumDataPoint=dataProjection%minimumErrorDataPoint
    minimumError=dataProjection%minimumError
 
    EXITS("DataProjection_ResultMinimumErrorGet")
    RETURN
999 ERRORSEXITS("DataProjection_ResultMinimumErrorGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_ResultMinimumErrorGet

  !
  !================================================================================================================================
  !

  !>Gets the projection Root Mean Squared (RMS) error for a data projection.
  SUBROUTINE DataProjection_ResultRMSErrorGet(dataProjection,rmsError,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to get the RMS error for
    REAL(DP), INTENT(OUT) :: rmsError !<On exit, the RMS error for the data projection.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("DataProjection_ResultRMSErrorGet",err,error,*999)

    CALL DataProjection_AssertIsFinished(dataProjection,err,error,*999)
    CALL DataProjection_AssertIsProjected(dataProjection,err,error,*999)
    
    rmsError=dataProjection%rmsError
 
    EXITS("DataProjection_ResultRMSErrorGet")
    RETURN
999 ERRORSEXITS("DataProjection_ResultRMSErrorGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_ResultRMSErrorGet

  !
  !================================================================================================================================
  !

  !>Finds an returns a pointer to a data projection identified by a user number.
  SUBROUTINE DataProjection_UserNumberFind(dataPoints,userNumber,dataProjection,err,error,*)

    !Argument variables
    TYPE(DataPointsType), POINTER :: dataPoints !<A pointer to the data points to find the user number for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number to find
    TYPE(DataProjectionType), POINTER :: dataProjection !<On exit, the pointer to the data projection with the specified user number. If no data projection with the specified user number exists the pointer is returned NULL. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: projectionIdx
    TYPE(DataProjectionType), POINTER :: listDataProjection
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("DataProjection_UserNumberFind",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(dataPoints)) CALL FlagError("Data points is not associated.",err,error,*999)
    IF(ASSOCIATED(dataProjection)) CALL FlagError("Data projection is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(dataPoints%dataProjections)) &
      & CALL FlagError("Data points data projections is not associated.",err,error,*999)
#endif    
    
    NULLIFY(dataProjection)
    IF(ALLOCATED(dataPoints%dataProjections%dataProjections)) THEN
      projectionIdx=1
      DO WHILE(projectionIdx<=dataPoints%dataProjections%numberOfDataProjections)
        listDataProjection=>dataPoints%dataProjections%dataProjections(projectionIdx)%ptr
#ifdef WITH_PRECHECKS        
        IF(.NOT.ASSOCIATED(listDataProjection)) THEN
          localError="The data points data projections is not associated for projection index "// &
            & TRIM(NumberToVString(projectionIdx,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
#endif        
        IF(listDataProjection%userNumber==userNumber) THEN
          dataProjection=>dataPoints%dataProjections%dataProjections(projectionIdx)%ptr
          EXIT
        ENDIF
      ENDDO
    ENDIF
     
    EXITS("DataProjection_UserNumberFind")
    RETURN
999 NULLIFY(dataProjection)
998 ERRORSEXITS("DataProjection_UserNumberFind",err,error)    
    RETURN 1
   
  END SUBROUTINE DataProjection_UserNumberFind
  
  !
  !================================================================================================================================
  !

END MODULE DataProjectionAccessRoutines
