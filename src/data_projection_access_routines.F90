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

  PUBLIC DataProjection_DataPointsGet

  PUBLIC DataProjection_DecompositionGet

  PUBLIC DataProjection_ProjectionFieldGet

  PUBLIC DataProjection_ResultMaximumErrorGet

  PUBLIC DataProjection_ResultMinimumErrorGet

  PUBLIC DataProjection_ResultRMSErrorGet

CONTAINS

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

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(ASSOCIATED(dataPoints)) CALL FlagError("Data points is already associated.",err,error,*999)
    
    dataPoints=>dataProjection%dataPoints
    IF(.NOT.ASSOCIATED(dataPoints)) CALL FlagError("Data projection data points is not associated.",err,error,*999)
 
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
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition !<On exit, the decomposition of the data projection. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
   
    ENTERS("DataProjection_DecompositionGet",err,error,*999)

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(ASSOCIATED(decomposition)) CALL FlagError("Decomposition is already associated.",err,error,*999)
    
    decomposition=>dataProjection%decomposition
    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Data projection decomposition is not associated.",err,error,*999)
 
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
    TYPE(FIELD_TYPE), POINTER :: projectionField !<On exit, the projection field of the data projection. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
   
    ENTERS("DataProjection_ProjectionFieldGet",err,error,*999)

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(ASSOCIATED(projectionField)) CALL FlagError("Projection field is already associated.",err,error,*999)
    
    projectionField=>dataProjection%projectionField
    IF(.NOT.ASSOCIATED(projectionField)) CALL FlagError("Data projection projection field is not associated.",err,error,*999)
 
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

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionFinished) CALL FlagError("Data projection has not been finished.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionProjected) CALL FlagError("Data projection has not been projected.",err,error,*999)
    
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

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionFinished) CALL FlagError("Data projection has not been finished.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionProjected) CALL FlagError("Data projection has not been projected.",err,error,*999)
    
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

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionFinished) CALL FlagError("Data projection has not been finished.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionProjected) CALL FlagError("Data projection has not been projected.",err,error,*999)
    
    rmsError=dataProjection%rmsError
 
    EXITS("DataProjection_ResultRMSErrorGet")
    RETURN
999 ERRORSEXITS("DataProjection_ResultRMSErrorGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_ResultRMSErrorGet

  !
  !================================================================================================================================
  !

END MODULE DataProjectionAccessRoutines
