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

  !> \addtogroup DataProjectionRoutines_DataProjectionTypes DataProjectionRoutines::DataProjectionTypes
  !> \brief Datapoint projection definition type parameters
  !> \see DataProjectionRoutines,OPENCMISS_DataProjectionTypes
  !>@{ 
  INTEGER(INTG), PARAMETER :: DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE=1 !<The boundary line projection type for data projection, only projects to boundary lines of the mesh. \see DataProjectionRoutines,OPENCMISS_DataProjectionTypes
  INTEGER(INTG), PARAMETER :: DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE=2 !<The boundary face projection type for data projection, only projects to boundary faces of the mesh. \see DataProjectionRoutines,OPENCMISS_DataProjectionTypes
  INTEGER(INTG), PARAMETER :: DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE=3 !<The element projection type for data projection, projects to all elements in mesh. \see DataProjectionRoutines,OPENCMISS_DataProjectionTypes
  !>@}

  !> \addtogroup DataProjectionRoutines_DataProjectionDistanceRelations DataProjectionRoutines::DataProjectionDistanceRelations
 !> \brief Datapoint projection distance relations to select data points based on distance.
  !> \see DataProjectionRoutines,OPENCMISS_DataProjectionExitTags
  !>@{ 
  INTEGER(INTG), PARAMETER :: DATA_PROJECTION_DISTANCE_GREATER=1 !<Data projection distance relation is greater than \see DataProjectionRoutines,OPENCMISS_DataProjectionDistanceRelations
  INTEGER(INTG), PARAMETER :: DATA_PROJECTION_DISTANCE_GREATER_EQUAL=2 !<Data projection distance relation is greater than or equal \see DataProjectionRoutines,OPENCMISS_DataProjectionDistanceRelations
  INTEGER(INTG), PARAMETER :: DATA_PROJECTION_DISTANCE_LESS=3 !<Data projection distance relation is less than \see DataProjectionRoutines,OPENCMISS_DataProjectionDistanceRelations
  INTEGER(INTG), PARAMETER :: DATA_PROJECTION_DISTANCE_LESS_EQUAL=4 !<Data projection distance relation is less than or equal \see DataProjectionRoutines,OPENCMISS_DataProjectionDistanceRelations
  !>@}
  
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

  !>Gets the label for a data projection.
  INTERFACE DataProjection_LabelGet
    MODULE PROCEDURE DataProjection_LabelGetC
    MODULE PROCEDURE DataProjection_LabelGetVS
  END INTERFACE DataProjection_LabelGet
  
  PUBLIC DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE,DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE, &
    & DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE

  PUBLIC DATA_PROJECTION_DISTANCE_GREATER,DATA_PROJECTION_DISTANCE_GREATER_EQUAL,DATA_PROJECTION_DISTANCE_LESS, &
    & DATA_PROJECTION_DISTANCE_LESS_EQUAL
  
  PUBLIC DATA_PROJECTION_CANCELLED,DATA_PROJECTION_EXIT_TAG_CONVERGED,DATA_PROJECTION_EXIT_TAG_BOUNDS, &
    & DATA_PROJECTION_EXIT_TAG_MAX_ITERATION,DATA_PROJECTION_EXIT_TAG_NO_ELEMENT

  PUBLIC DataProjection_AbsoluteToleranceGet
  
  PUBLIC DataProjection_AssertIsFinished,DataProjection_AssertNotFinished

  PUBLIC DataProjection_AssertIsProjected,DataProjection_AssertNotProjected

  PUBLIC DataProjection_DataPointsGet

  PUBLIC DataProjection_DecompositionGet

  PUBLIC DataProjection_LabelGet

  PUBLIC DataProjection_MaximumInterationUpdateGet

  PUBLIC DataProjection_MaximumNumberOfIterationsGet

  PUBLIC DataProjection_NumberOfClosestElementsGet

  PUBLIC DataProjection_ProjectionFieldGet

  PUBLIC DataProjection_ProjectionTypeGet

  PUBLIC DataProjection_RelativeToleranceGet

  PUBLIC DataProjection_ResultDistanceGlobalGet
  
  PUBLIC DataProjection_ResultDistanceUserGet
  
  PUBLIC DataProjection_ResultElementNumberGlobalGet

  PUBLIC DataProjection_ResultElementNumberUserGet

  PUBLIC DataProjection_ResultElementFaceNumberGlobalGet

  PUBLIC DataProjection_ResultElementFaceNumberUserGet

  PUBLIC DataProjection_ResultElementLineNumberGlobalGet
  
  PUBLIC DataProjection_ResultElementLineNumberUserGet
  
  PUBLIC DataProjection_ResultElementXiGlobalGet

  PUBLIC DataProjection_ResultExitTagGlobalGet

  PUBLIC DataProjection_ResultExitTagUserGet

  PUBLIC DataProjection_ResultMaximumErrorGet

  PUBLIC DataProjection_ResultMinimumErrorGet

  PUBLIC DataProjection_ResultProjectionVectorGlobalGet
  
  PUBLIC DataProjection_ResultProjectionVectorUserGet
  
  PUBLIC DataProjection_ResultRMSErrorGet

  PUBLIC DataProjection_ResultXiGlobalGet

  PUBLIC DataProjection_ResultXiUserGet

  PUBLIC DataProjection_StartingXiGet

  PUBLIC DataProjection_UserNumberFind
  
CONTAINS

  !
  !================================================================================================================================
  !

  !>Gets the absolute tolerance for a data projection.
  SUBROUTINE DataProjection_AbsoluteToleranceGet(dataProjection,absoluteTolerance,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to get the absolute tolerance for
    REAL(DP), INTENT(OUT) :: absoluteTolerance !<On exit, the absolute tolerance of the specified data projection
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("DataProjection_AbsoluteToleranceGet",err,error,*999)

    CALL DataProjection_AssertIsFinished(dataProjection,err,error,*999)
    
    absoluteTolerance=dataProjection%absoluteTolerance       
     
    EXITS("DataProjection_AbsoluteToleranceGet")
    RETURN
999 ERRORSEXITS("DataProjection_AbsoluteToleranceGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_AbsoluteToleranceGet

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

  !>Gets the label for a data projection for varying string labels. \see OpenCMISS::Iron::cmfe_DataProjection_LabelGet
  SUBROUTINE DataProjection_LabelGetVS(dataProjection,label,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to get the label for
    TYPE(VARYING_STRING), INTENT(OUT) :: label !<the label to get
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("DataProjection_LabelGetVS",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
#endif    
    
    label=dataProjection%label
     
    EXITS("DataProjection_LabelGetVS")
    RETURN
999 ERRORSEXITS("DataProjection_LabelGetVS",err,error)
    RETURN 1

  END SUBROUTINE DataProjection_LabelGetVS

  !
  !================================================================================================================================
  !

  !>Gets the label for a data projection for character labels. \see OpenCMISS::Iron::cmfe_DataProjection_LabelGet
  SUBROUTINE DataProjection_LabelGetC(dataProjection,label,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to get the label for
    CHARACTER(LEN=*), INTENT(OUT) :: label !<the label to get
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: cLength,vsLength
    
    ENTERS("DataProjection_LabelGetC",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
#endif    
    
    cLength=LEN(label)
    vsLength=LEN_TRIM(dataProjection%label)
    IF(cLength>vsLength) THEN
      label=CHAR(LEN_TRIM(dataProjection%label))
    ELSE
      label=CHAR(dataProjection%label,cLength)
    ENDIF
    
    EXITS("DataProjection_LabelGetC")
    RETURN
999 ERRORSEXITS("DataProjection_LabelGetC",err,error)
    RETURN 1

  END SUBROUTINE DataProjection_LabelGetC

  !
  !================================================================================================================================
  !
  
  !>Gets the maximum iteration update for a data projection.
  SUBROUTINE DataProjection_MaximumInterationUpdateGet(dataProjection,maximumIterationUpdate,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to get the maximum iteration update for
    REAL(DP), INTENT(OUT) :: maximumIterationUpdate !<On exit, the maximum iteration update of the specified data projection
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("DataProjection_MaximumInterationUpdateGet",err,error,*999)

    CALL DataProjection_AssertIsFinished(dataProjection,err,error,*999)
    
    maximumIterationUpdate=dataProjection%maximumIterationUpdate       
     
    EXITS("DataProjection_MaximumInterationUpdateGet")
    RETURN
999 ERRORSEXITS("DataProjection_MaximumInterationUpdateGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_MaximumInterationUpdateGet

  !
  !================================================================================================================================
  !
  
  !>Gets the maximum number of iterations for a data projection.
  SUBROUTINE DataProjection_MaximumNumberOfIterationsGet(dataProjection,maximumNumberOfIterations,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to get the maximum number of iterations for
    INTEGER(INTG), INTENT(OUT) :: maximumNumberOfIterations !<On exit, the maximum number of iterations of the specified data projection
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("DataProjection_MaximumNumberOfIterationsGet",err,error,*999)

    CALL DataProjection_AssertIsFinished(dataProjection,err,error,*999)
    
    maximumNumberOfIterations=dataProjection%maximumNumberOfIterations       
    
    EXITS("DataProjection_MaximumNumberOfIterationsGet")
    RETURN
999 ERRORSEXITS("DataProjection_MaximumNumberOfIterationsGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_MaximumNumberOfIterationsGet

  !
  !================================================================================================================================
  !
  
  !>Gets the number of closest elements for a data projection.
  SUBROUTINE DataProjection_NumberOfClosestElementsGet(dataProjection,numberOfClosestElements,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to get the number of closest elements for
    INTEGER(INTG), INTENT(OUT) :: numberOfClosestElements !<On exit, the number of closest elements of the specified data projection
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("DataProjection_NumberOfClosestElementsGet",err,error,*999)

    CALL DataProjection_AssertIsFinished(dataProjection,err,error,*999)
    
    numberOfClosestElements=dataProjection%numberOfClosestElements       
    
    EXITS("DataProjection_NumberOfClosestElementsGet")
    RETURN
999 ERRORSEXITS("DataProjection_NumberOfClosestElementsGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_NumberOfClosestElementsGet

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

  !>Gets the projection type for a data projection.
  SUBROUTINE DataProjection_ProjectionTypeGet(dataProjection,projectionType,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to get the projection type for
    INTEGER(INTG), INTENT(OUT) :: projectionType !<On exit, the projection type of the specified data projection
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("DataProjection_ProjectionTypeGet",err,error,*999)

    CALL DataProjection_AssertIsFinished(dataProjection,err,error,*999)
    
    projectionType=dataProjection%projectionType       
    
    EXITS("DataProjection_ProjectionTypeGet")
    RETURN
999 ERRORSEXITS("DataProjection_ProjectionTypeGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_ProjectionTypeGet

  !
  !================================================================================================================================
  !

  !>Gets the relative tolerance for a data projection.
  SUBROUTINE DataProjection_RelativeToleranceGet(dataProjection,relativeTolerance,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to get the relative tolerance for
    REAL(DP), INTENT(OUT) :: relativeTolerance !<On exit, the relative tolerance of the specified data projection
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("DataProjection_RelativeToleranceGet",err,error,*999)

    CALL DataProjection_AssertIsFinished(dataProjection,err,error,*999)
    
    relativeTolerance=dataProjection%relativeTolerance       
   
    EXITS("DataProjection_RelativeToleranceGet")
    RETURN
999 ERRORSEXITS("DataProjection_RelativeToleranceGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_RelativeToleranceGet

  !
  !================================================================================================================================
  !

  !>Gets the projection distance for a data point identified by a given global number.
  SUBROUTINE DataProjection_ResultDistanceGlobalGet(dataProjection,dataPointGlobalNumber,projectionDistance,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection for which projection result is stored
    INTEGER(INTG), INTENT(IN) :: dataPointGlobalNumber !<The data projection global number to get the projection distance for
    REAL(DP), INTENT(OUT) :: projectionDistance !<On exit, the projection distance of the specified data point
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("DataProjection_ResultDistanceGlobalGet",err,error,*999)

    CALL DataProjection_AssertIsFinished(dataProjection,err,error,*999)
    CALL DataProjection_AssertIsProjected(dataProjection,err,error,*999)
#ifdef WITH_PRECHECKS    
    IF(.NOT.ALLOCATED(dataProjection%dataProjectionResults)) THEN
      localError="The data projection results array is not allocated for data projection number "// &
        & TRIM(NumberToVString(dataProjection%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(dataPointGlobalNumber<1.OR.dataPointGlobalNumber>SIZE(dataProjection%dataProjectionResults,1)) THEN
      localError="The specified data point global number of "//TRIM(NumberToVString(dataPointGlobalNumber,"*",err,error))// &
        & " is invalid. The global number should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(dataProjection%dataProjectionResults,1),"*",err,error))//"."
    ENDIF
#endif    
    
    projectionDistance=dataProjection%dataProjectionResults(dataPointGlobalNumber)%distance

    EXITS("DataProjection_ResultDistanceGlobalGet")
    RETURN
999 ERRORSEXITS("DataProjection_ResultDistanceGlobalGet",err,error)
    RETURN 1

  END SUBROUTINE DataProjection_ResultDistanceGlobalGet

  !
  !================================================================================================================================
  !

  !>Gets the projection distance for a data point identified by a given user number.
  SUBROUTINE DataProjection_ResultDistanceUserGet(dataProjection,dataPointUserNumber,projectionDistance,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection for which projection result is stored
    INTEGER(INTG), INTENT(IN) :: dataPointUserNumber !<The Data projection user number to get the projection distance for
    REAL(DP), INTENT(OUT) :: projectionDistance !<On exit, the projection distance of the specified data point
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointGlobalNumber
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("DataProjection_ResultDistanceUserGet",err,error,*999)

    CALL DataProjection_AssertIsFinished(dataProjection,err,error,*999)
    CALL DataProjection_AssertIsProjected(dataProjection,err,error,*999)
#ifdef WITH_PRECHECKS    
    IF(.NOT.ALLOCATED(dataProjection%dataProjectionResults)) THEN
      localError="The data projection results array is not allocated for data projection number "// &
        & TRIM(NumberToVString(dataProjection%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
     
    CALL DataProjection_DataPointGlobalNumberGet(dataProjection,dataPointUserNumber,dataPointGlobalNumber,err,error,*999)
    projectionDistance=dataProjection%dataProjectionResults(dataPointGlobalNumber)%distance

    EXITS("DataProjection_ResultDistanceUserGet")
    RETURN
999 ERRORSEXITS("DataProjection_ResultDistanceUserGet",err,error)
    RETURN 1

  END SUBROUTINE DataProjection_ResultDistanceUserGet

  !
  !================================================================================================================================
  !

  !>Gets the projection element number for a data point identified by a given global number.
  SUBROUTINE DataProjection_ResultElementNumberGlobalGet(dataProjection,dataPointGlobalNumber,projectionElementNumber,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection for which projection result is stored
    INTEGER(INTG), INTENT(IN) :: dataPointGlobalNumber !<The Data projection global number to get the projection element number for
    INTEGER(INTG), INTENT(OUT) :: projectionElementNumber !<On exit, the projection element number of the specified global data point
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
   
    ENTERS("DataProjection_ResultElementNumberGlobalGet",err,error,*999)

    CALL DataProjection_AssertIsFinished(dataProjection,err,error,*999)
    CALL DataProjection_AssertIsProjected(dataProjection,err,error,*999)
#ifdef WITH_PRECHECKS    
    IF(.NOT.ALLOCATED(dataProjection%dataProjectionResults)) THEN
      localError="The data projection results array is not allocated for data projection number "// &
        & TRIM(NumberToVString(dataProjection%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(dataPointGlobalNumber<1.OR.dataPointGlobalNumber>SIZE(dataProjection%dataProjectionResults,1)) THEN
      localError="The specified data point global number of "//TRIM(NumberToVString(dataPointGlobalNumber,"*",err,error))// &
        & " is invalid. The global number should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(dataProjection%dataProjectionResults,1),"*",err,error))//"."
    ENDIF
#endif    
    
    projectionElementNumber=dataProjection%dataProjectionResults(dataPointGlobalNumber)%elementNumber
 
    EXITS("DataProjection_ResultElementNumberGlobalGet")
    RETURN
999 ERRORSEXITS("DataProjection_ResultElementNumberGlobalGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_ResultElementNumberGlobalGet
  
  !
  !================================================================================================================================
  !

  !>Gets the projection element number for a data point identified by a given user number.
  SUBROUTINE DataProjection_ResultElementNumberUserGet(dataProjection,dataPointUserNumber,projectionElementNumber,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection for which projection result is stored
    INTEGER(INTG), INTENT(IN) :: dataPointUserNumber !<The Data projection user number to get the projection element number for
    INTEGER(INTG), INTENT(OUT) :: projectionElementNumber !<On exit, the projection element number of the specified user data point
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointGlobalNumber
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
   
    ENTERS("DataProjection_ResultElementNumberUserGet",err,error,*999)

    CALL DataProjection_AssertIsFinished(dataProjection,err,error,*999)
    CALL DataProjection_AssertIsProjected(dataProjection,err,error,*999)
#ifdef WITH_PRECHECKS    
    IF(.NOT.ALLOCATED(dataProjection%dataProjectionResults)) THEN
      localError="The data projection results array is not allocated for data projection number "// &
        & TRIM(NumberToVString(dataProjection%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    CALL DataProjection_DataPointGlobalNumberGet(dataProjection,dataPointUserNumber,dataPointGlobalNumber,err,error,*999)
    projectionElementNumber=dataProjection%dataProjectionResults(dataPointGlobalNumber)%elementNumber
 
    EXITS("DataProjection_ResultElementNumberUserGet")
    RETURN
999 ERRORSEXITS("DataProjection_ResultElementNumberUserGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_ResultElementNumberUserGet
  
  !
  !================================================================================================================================
  !

  !>Gets the projection element face number for a data point identified by a given global number.
  SUBROUTINE DataProjection_ResultElementFaceNumberGlobalGet(dataProjection,dataPointGlobalNumber,projectionElementFaceNumber &
    & ,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection for which projection result is stored
    INTEGER(INTG), INTENT(IN) :: dataPointGlobalNumber !<The data projection global number to get the projection element face number for
    INTEGER(INTG), INTENT(OUT) :: projectionElementFaceNumber !<On exit, the projection element face number of the specified global data point
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("DataProjection_ResultElementFaceNumberGlobalGet",err,error,*999)

    CALL DataProjection_AssertIsFinished(dataProjection,err,error,*999)
    CALL DataProjection_AssertIsProjected(dataProjection,err,error,*999)
#ifdef WITH_PRECHECKS    
    IF(.NOT.ALLOCATED(dataProjection%dataProjectionResults)) THEN
      localError="The data projection results array is not allocated for data projection number "// &
        & TRIM(NumberToVString(dataProjection%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(dataPointGlobalNumber<1.OR.dataPointGlobalNumber>SIZE(dataProjection%dataProjectionResults,1)) THEN
      localError="The specified data point global number of "//TRIM(NumberToVString(dataPointGlobalNumber,"*",err,error))// &
        & " is invalid. The global number should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(dataProjection%dataProjectionResults,1),"*",err,error))//"."
    ENDIF
    !Check if boundary faces projection type was set
    IF(dataProjection%projectionType/=DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE) THEN
      localError="The data projection type of "//TRIM(NumberToVString(dataProjection%projectionType,"*",err,error))// &
        & " for data projection number "//TRIM(NumberToVString(dataProjection%userNumber,"*",err,error))// &
        & " is not set to a boundary faces projection type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    projectionElementFaceNumber=dataProjection%dataProjectionResults(dataPointGlobalNumber)%elementLineFaceNumber
 
    EXITS("DataProjection_ResultElementFaceNumberGlobalGet")
    RETURN
999 ERRORSEXITS("DataProjection_ResultElementFaceNumberGlobalGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_ResultElementFaceNumberGlobalGet

  !
  !================================================================================================================================
  !

  !>Gets the projection element face number for a data point identified by a given user number.
  SUBROUTINE DataProjection_ResultElementFaceNumberUserGet(dataProjection,dataPointUserNumber,projectionElementFaceNumber &
    & ,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection for which projection result is stored
    INTEGER(INTG), INTENT(IN) :: dataPointUserNumber !<The Data projection user number to get the projection element face number for
    INTEGER(INTG), INTENT(OUT) :: projectionElementFaceNumber !<On exit, the projection element face number of the specified user data point
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointGlobalNumber
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("DataProjection_ResultElementFaceNumberUserGet",err,error,*999)

    CALL DataProjection_AssertIsFinished(dataProjection,err,error,*999)
    CALL DataProjection_AssertIsProjected(dataProjection,err,error,*999)
#ifdef WITH_PRECHECKS    
    IF(.NOT.ALLOCATED(dataProjection%dataProjectionResults)) THEN
      localError="The data projection results array is not allocated for data projection number "// &
        & TRIM(NumberToVString(dataProjection%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    !Check if boundary faces projection type was set
    IF(dataProjection%projectionType/=DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE) THEN
      localError="The data projection type of "//TRIM(NumberToVString(dataProjection%projectionType,"*",err,error))// &
        & " for data projection number "//TRIM(NumberToVString(dataProjection%userNumber,"*",err,error))// &
        & " is not set to a boundary faces projection type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    CALL DataProjection_DataPointGlobalNumberGet(dataProjection,dataPointUserNumber,dataPointGlobalNumber,err,error,*999)
    projectionElementFaceNumber=dataProjection%dataProjectionResults(dataPointGlobalNumber)%elementLineFaceNumber

    EXITS("DataProjection_ResultElementFaceNumberUserGet")
    RETURN
999 ERRORSEXITS("DataProjection_ResultElementFaceNumberUserGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_ResultElementFaceNumberUserGet

  !
  !================================================================================================================================
  !

  !>Gets the projection element line number for a data point identified by a given global number.
  SUBROUTINE DataProjection_ResultElementLineNumberGlobalGet(dataProjection,dataPointGlobalNumber,projectionElementLineNumber &
    & ,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection for which projection result is stored
    INTEGER(INTG), INTENT(IN) :: dataPointGlobalNumber !<The data projection user number to get the element line number distance for
    INTEGER(INTG), INTENT(OUT) :: projectionElementLineNumber !<On exit, the projection element line number of the specified global data point
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("DataProjection_ResultElementLineNumberGlobalGet",err,error,*999)

    CALL DataProjection_AssertIsFinished(dataProjection,err,error,*999)
    CALL DataProjection_AssertIsProjected(dataProjection,err,error,*999)
#ifdef WITH_PRECHECKS    
    IF(.NOT.ALLOCATED(dataProjection%dataProjectionResults)) THEN
      localError="The data projection results array is not allocated for data projection number "// &
        & TRIM(NumberToVString(dataProjection%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(dataPointGlobalNumber<1.OR.dataPointGlobalNumber>SIZE(dataProjection%dataProjectionResults,1)) THEN
      localError="The specified data point global number of "//TRIM(NumberToVString(dataPointGlobalNumber,"*",err,error))// &
        & " is invalid. The global number should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(dataProjection%dataProjectionResults,1),"*",err,error))//"."
    ENDIF
    !Check if boundary lines projection type was set
    IF(dataProjection%projectionType/=DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE) THEN
      localError="The data projection type of "//TRIM(NumberToVString(dataProjection%projectionType,"*",err,error))// &
        & " for data projection number "//TRIM(NumberToVString(dataProjection%userNumber,"*",err,error))// &
        & " is not set to a boundary lines projection type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    projectionElementLineNumber=dataProjection%dataProjectionResults(dataPointGlobalNumber)%elementLineFaceNumber

    EXITS("DataProjection_ResultElementLineNumberGlobalGet")
    RETURN
999 ERRORSEXITS("DataProjection_ResultElementLineNumberGlobalGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_ResultElementLineNumberGlobalGet

  !
  !================================================================================================================================
  !

  !>Gets the projection element line number for a data point identified by a given user number.
  SUBROUTINE DataProjection_ResultElementLineNumberUserGet(dataProjection,dataPointUserNumber,projectionElementLineNumber &
    & ,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection for which projection result is stored
    INTEGER(INTG), INTENT(IN) :: dataPointUserNumber !<The Data projection user number to get the element line number distance for
    INTEGER(INTG), INTENT(OUT) :: projectionElementLineNumber !<On exit, the projection element line number of the specified user data point
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointGlobalNumber
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("DataProjection_ResultElementLineNumberUserGet",err,error,*999)

    CALL DataProjection_AssertIsFinished(dataProjection,err,error,*999)
    CALL DataProjection_AssertIsProjected(dataProjection,err,error,*999)
#ifdef WITH_PRECHECKS    
    IF(.NOT.ALLOCATED(dataProjection%dataProjectionResults)) THEN
      localError="The data projection results array is not allocated for data projection number "// &
        & TRIM(NumberToVString(dataProjection%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    !Check if boundary lines projection type was set
    IF(dataProjection%projectionType/=DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE) THEN
      localError="The data projection type of "//TRIM(NumberToVString(dataProjection%projectionType,"*",err,error))// &
        & " for data projection number "//TRIM(NumberToVString(dataProjection%userNumber,"*",err,error))// &
        & " is not set to a boundary lines projection type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    CALL DataProjection_DataPointGlobalNumberGet(dataProjection,dataPointUserNumber,dataPointGlobalNumber,err,error,*999)
    projectionElementLineNumber=dataProjection%dataProjectionResults(dataPointGlobalNumber)%elementLineFaceNumber

    EXITS("DataProjection_ResultElementLineNumberUserGet")
    RETURN
999 ERRORSEXITS("DataProjection_ResultElementLineNumberUserGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_ResultElementLineNumberUserGet

  !
  !================================================================================================================================
  !

  !>Gets the projection element xi for a data point identified by a given global number.
  SUBROUTINE DataProjection_ResultElementXiGlobalGet(dataProjection,dataPointGlobalNumber,elementXi,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection for which projection result is stored
    INTEGER(INTG), INTENT(IN) :: dataPointGlobalNumber !<The data projection global number to get the projection element xi for
    REAL(DP), INTENT(OUT) :: elementXi(:) !<On exit, the projection element xi of the specified global data point
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: elementXiSize
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("DataProjection_ResultXiGlobalGet",err,error,*999)

    CALL DataProjection_AssertIsFinished(dataProjection,err,error,*999)
    CALL DataProjection_AssertIsProjected(dataProjection,err,error,*999)
#ifdef WITH_PRECHECKS    
    IF(.NOT.ALLOCATED(dataProjection%dataProjectionResults)) THEN
      localError="The data projection results array is not allocated for data projection number "// &
        & TRIM(NumberToVString(dataProjection%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(dataPointGlobalNumber<1.OR.dataPointGlobalNumber>SIZE(dataProjection%dataProjectionResults,1)) THEN
      localError="The specified data point global number of "//TRIM(NumberToVString(dataPointGlobalNumber,"*",err,error))// &
        & " is invalid. The global number should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(dataProjection%dataProjectionResults,1),"*",err,error))//"."
    ENDIF
    IF(.NOT.ALLOCATED(dataProjection%dataProjectionResults(dataPointGlobalNumber)%elementXi)) THEN
      localError="The element xi array is not allocated for global data point number "// &
        & TRIM(NumberToVString(dataPointGlobalNumber,"*",err,error))//" for data projection number "// &
        & TRIM(NumberToVString(dataProjection%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(elementXi,1)<SIZE(dataProjection%dataProjectionResults(dataPointGlobalNumber)%elementXi,1)) THEN
      localError="The specified element xi has size of "//TRIM(NumberToVString(SIZE(elementXi,1),"*",err,error))// &
        & " but it needs to have size of >= "// &
        & TRIM(NumberToVString(SIZE(dataProjection%dataProjectionResults(dataPointGlobalNumber)%elementXi,1),"*",err,error))//"." 
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    elementXiSize=SIZE(dataProjection%dataProjectionResults(dataPointGlobalNumber)%elementXi,1)
    elementXi(1:elementXiSize)=dataProjection%dataProjectionResults(dataPointGlobalNumber)%elementXi(1:elementXiSize)

    EXITS("DataProjection_ResultElementXiGlobalGet")
    RETURN
999 ERRORSEXITS("DataProjection_ResultElementXiGlobalGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_ResultElementXiGlobalGet

  !
  !================================================================================================================================
  !

  !>Gets the projection exit tag for a data point identified by a given global number.
  SUBROUTINE DataProjection_ResultExitTagGlobalGet(dataProjection,dataPointGlobalNumber,projectionExitTag,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection for which projection result is stored
    INTEGER(INTG), INTENT(IN) :: dataPointGlobalNumber !<The Data projection global number to get the projection exit tag for
    INTEGER(INTG), INTENT(OUT) :: projectionExitTag !<On exit, the projection exit tag of the specified global data point
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("DataProjection_ResultExitTagGlobalGet",err,error,*999)

    CALL DataProjection_AssertIsFinished(dataProjection,err,error,*999)
    CALL DataProjection_AssertIsProjected(dataProjection,err,error,*999)
#ifdef WITH_PRECHECKS    
    IF(.NOT.ALLOCATED(dataProjection%dataProjectionResults)) THEN
      localError="The data projection results array is not allocated for data projection number "// &
        & TRIM(NumberToVString(dataProjection%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(dataPointGlobalNumber<1.OR.dataPointGlobalNumber>SIZE(dataProjection%dataProjectionResults,1)) THEN
      localError="The specified data point global number of "//TRIM(NumberToVString(dataPointGlobalNumber,"*",err,error))// &
        & " is invalid. The global number should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(dataProjection%dataProjectionResults,1),"*",err,error))//"."
    ENDIF
#endif    
     
    projectionExitTag=dataProjection%dataProjectionResults(dataPointGlobalNumber)%exitTag

    EXITS("DataProjection_ResultExitTagGlobalGet")
    RETURN
999 ERRORSEXITS("DataProjection_ResultExitTagGlobalGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_ResultExitTagGlobalGet

  !
  !================================================================================================================================
  !

  !>Gets the projection exit tag for a data point identified by a given user number.
  SUBROUTINE DataProjection_ResultExitTagUserGet(dataProjection,dataPointUserNumber,projectionExitTag,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection for which projection result is stored
    INTEGER(INTG), INTENT(IN) :: dataPointUserNumber !<The Data projection user number to get the projection exit tag for
    INTEGER(INTG), INTENT(OUT) :: projectionExitTag !<On exit, the projection exit tag of the specified user data point
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointGlobalNumber
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("DataProjection_ResultExitTagUserGet",err,error,*999)

    CALL DataProjection_AssertIsFinished(dataProjection,err,error,*999)
    CALL DataProjection_AssertIsProjected(dataProjection,err,error,*999)
#ifdef WITH_PRECHECKS    
    IF(.NOT.ALLOCATED(dataProjection%dataProjectionResults)) THEN
      localError="The data projection results array is not allocated for data projection number "// &
        & TRIM(NumberToVString(dataProjection%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    CALL DataProjection_DataPointGlobalNumberGet(dataProjection,dataPointUserNumber,dataPointGlobalNumber,err,error,*999)
    projectionExitTag=dataProjection%dataProjectionResults(dataPointGlobalNumber)%exitTag

    EXITS("DataProjection_ResultExitTagUserGet")
    RETURN
999 ERRORSEXITS("DataProjection_ResultExitTagUserGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_ResultExitTagUserGet

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

  !>Gets the projection vector for a data point identified by a given global number.
  SUBROUTINE DataProjection_ResultProjectionVectorGlobalGet(dataProjection,dataPointGlobalNumber,projectionVector,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection for which projection result is stored
    INTEGER(INTG), INTENT(IN) :: dataPointGlobalNumber !<The data point global number to get the projection xi for
    REAL(DP), INTENT(OUT) :: projectionVector(:) !<On exit, the projection vector of the specified global data point
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("DataProjection_ResultProjectionVectorGlobalGet",err,error,*999)

    CALL DataProjection_AssertIsFinished(dataProjection,err,error,*999)
    CALL DataProjection_AssertIsProjected(dataProjection,err,error,*999)
#ifdef WITH_PRECHECKS    
    IF(.NOT.ALLOCATED(dataProjection%dataProjectionResults)) THEN
      localError="The data projection results array is not allocated for data projection number "// &
        & TRIM(NumberToVString(dataProjection%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(dataPointGlobalNumber<1.OR.dataPointGlobalNumber>SIZE(dataProjection%dataProjectionResults,1)) THEN
      localError="The specified data point global number of "//TRIM(NumberToVString(dataPointGlobalNumber,"*",err,error))// &
        & " is invalid. The global number should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(dataProjection%dataProjectionResults,1),"*",err,error))//"."
    ENDIF
    IF(SIZE(projectionVector,1)<dataProjection%numberOfCoordinates) THEN
      localError="The specified projection vector has a size of "// &
        & TRIM(NumberToVString(SIZE(projectionVector,1),"*",err,error))//" but it needs to have size of >= "// &
        & TRIM(NumberToVString(dataProjection%numberOfCoordinates,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
   
   projectionVector(1:dataProjection%numberOfCoordinates)=dataProjection%dataProjectionResults(dataPointGlobalNumber)% &
      & projectionVector(1:dataProjection%numberOfCoordinates)

    EXITS("DataProjection_ResultProjectionVectorGlobalGet")
    RETURN
999 ERRORSEXITS("DataProjection_ResultProjectionVectorGlobalGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_ResultProjectionVectorGlobalGet

  !
  !================================================================================================================================
  !

  !>Gets the projection vector for a data point identified by a given user number.
  SUBROUTINE DataProjection_ResultProjectionVectorUserGet(dataProjection,dataPointUserNumber,projectionVector,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection for which projection result is stored
    INTEGER(INTG), INTENT(IN) :: dataPointUserNumber !<The Data projection user number to get the projection xi for
    REAL(DP), INTENT(OUT) :: projectionVector(:) !<On exit, the projection vector of the specified user data point
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointGlobalNumber
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif
    
    ENTERS("DataProjection_ResultProjectionVectorUserGet",err,error,*999)

    CALL DataProjection_AssertIsFinished(dataProjection,err,error,*999)
    CALL DataProjection_AssertIsProjected(dataProjection,err,error,*999)
#ifdef WITH_PRECHECKS    
    IF(.NOT.ALLOCATED(dataProjection%dataProjectionResults)) THEN
      localError="The data projection results array is not allocated for data projection number "// &
        & TRIM(NumberToVString(dataProjection%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(projectionVector,1)<dataProjection%numberOfCoordinates) THEN
      localError="The specified projection vector has a size of "// &
        & TRIM(NumberToVString(SIZE(projectionVector,1),"*",err,error))//" but it needs to have size of >= "// &
        & TRIM(NumberToVString(dataProjection%numberOfCoordinates,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
   
    CALL DataProjection_DataPointGlobalNumberGet(dataProjection,dataPointUserNumber,dataPointGlobalNumber,err,error,*999)
    projectionVector(1:dataProjection%numberOfCoordinates)=dataProjection%dataProjectionResults(dataPointGlobalNumber)% &
      & projectionVector(1:dataProjection%numberOfCoordinates)

    EXITS("DataProjection_ResultProjectionVectorUserGet")
    RETURN
999 ERRORSEXITS("DataProjection_ResultProjectionVectorUserGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_ResultProjectionVectorUserGet

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

  !>Gets the projection xi for a data point identified by a given global number.
  SUBROUTINE DataProjection_ResultXiGlobalGet(dataProjection,dataPointGlobalNumber,projectionXi,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection for which projection result is stored
    INTEGER(INTG), INTENT(IN) :: dataPointGlobalNumber !<The data projection global number to get the projection xi for
    REAL(DP), INTENT(OUT) :: projectionXi(:) !<On exit, the projection xi of the specified global data point
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("DataProjection_ResultXiGlobalGet",err,error,*999)

    CALL DataProjection_AssertIsFinished(dataProjection,err,error,*999)
    CALL DataProjection_AssertIsProjected(dataProjection,err,error,*999)
#ifdef WITH_PRECHECKS    
    IF(.NOT.ALLOCATED(dataProjection%dataProjectionResults)) THEN
      localError="The data projection results array is not allocated for data projection number "// &
        & TRIM(NumberToVString(dataProjection%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(dataPointGlobalNumber<1.OR.dataPointGlobalNumber>SIZE(dataProjection%dataProjectionResults,1)) THEN
      localError="The specified data point global number of "//TRIM(NumberToVString(dataPointGlobalNumber,"*",err,error))// &
        & " is invalid. The global number should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(dataProjection%dataProjectionResults,1),"*",err,error))//"."
    ENDIF
    IF(.NOT.ALLOCATED(dataProjection%dataProjectionResults(dataPointGlobalNumber)%xi)) THEN
      localError="The xi array is not allocated for global data point number "// &
        & TRIM(NumberToVString(dataPointGlobalNumber,"*",err,error))//" for data projection number "// &
        & TRIM(NumberToVString(dataProjection%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(projectionXi,1)<dataProjection%numberOfXi) THEN
      localError="The specified projection xi has size of "//TRIM(NumberToVString(SIZE(projectionXi,1),"*",err,error))// &
        & " but it needs to have size of >= "//TRIM(NumberToVString(dataProjection%numberOfXi,"*",err,error))//"." 
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
        
    projectionXi(1:dataProjection%numberOfXi)=dataProjection%dataProjectionResults(dataPointGlobalNumber)% &
      & xi(1:dataProjection%numberOfXi)

    EXITS("DataProjection_ResultXiGlobalGet")
    RETURN
999 ERRORSEXITS("DataProjection_ResultXiGlobalGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_ResultXiGlobalGet

  !
  !================================================================================================================================
  !

  !>Gets the projection xi for a data point identified by a given user number.
  SUBROUTINE DataProjection_ResultXiUserGet(dataProjection,dataPointUserNumber,projectionXi,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection for which projection result is stored
    INTEGER(INTG), INTENT(IN) :: dataPointUserNumber !<The data point user number to get the projection xi for
    REAL(DP), INTENT(OUT) :: projectionXi(:) !<On exit, the projection xi of the specified user data point
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointGlobalNumber
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("DataProjection_ResultXiUserGet",err,error,*999)

    CALL DataProjection_AssertIsFinished(dataProjection,err,error,*999)
    CALL DataProjection_AssertIsProjected(dataProjection,err,error,*999)
    CALL DataProjection_DataPointGlobalNumberGet(dataProjection,dataPointUserNumber,dataPointGlobalNumber,err,error,*999)
#ifdef WITH_PRECHECKS    
    IF(.NOT.ALLOCATED(dataProjection%dataProjectionResults)) THEN
      localError="The data projection results array is not allocated for data projection number "// &
        & TRIM(NumberToVString(dataProjection%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(dataProjection%dataProjectionResults(dataPointGlobalNumber)%xi)) THEN
      localError="The xi array is not allocated for user data point number "// &
        & TRIM(NumberToVString(dataPointUserNumber,"*",err,error))//" for data projection number "// &
        & TRIM(NumberToVString(dataProjection%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(projectionXi,1)<dataProjection%numberOfXi) THEN
      localError="The specified projection xi has size of "//TRIM(NumberToVString(SIZE(projectionXi,1),"*",err,error))// &
        & " but it needs to have size of >= "//TRIM(NumberToVString(dataProjection%numberOfXi,"*",err,error))//"." 
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
        
    projectionXi(1:dataProjection%numberOfXi)=dataProjection%dataProjectionResults(dataPointGlobalNumber)% &
      & xi(1:dataProjection%numberOfXi)

    EXITS("DataProjection_ResultXiUserGet")
    RETURN
999 ERRORSEXITS("DataProjection_ResultXiUserGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_ResultXiUserGet

  !
  !================================================================================================================================
  !

  !>Gets the starting xi for a data projection.
  SUBROUTINE DataProjection_StartingXiGet(dataProjection,startingXi,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to get the starting xi for
    REAL(DP), INTENT(OUT) :: startingXi(:) !<On exit, the starting xi of the specified data projection
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DataProjection_StartingXiGet",err,error,*999)

    CALL DataProjection_AssertIsFinished(dataProjection,err,error,*999)
    IF(SIZE(startingXi,1)<SIZE(dataProjection%startingXi,1)) THEN
      localError="The size of the specified starting xi array of "//TRIM(NumberToVString(SIZE(startingXi,1),"*",err,error))// &
        & " is too small. The size must be >= "//TRIM(NumberToVString(SIZE(dataProjection%startingXi,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    startingXi(1:SIZE(dataProjection%startingXi,1))=dataProjection%startingXi(1:SIZE(dataProjection%startingXi,1))
    
    EXITS("DataProjection_StartingXiGet")
    RETURN
999 ERRORSEXITS("DataProjection_StartingXiGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_StartingXiGet

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
