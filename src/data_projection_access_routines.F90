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
  USE Trees
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

  PUBLIC DataProjection_DataPointGlobalNumberGet

  PUBLIC DataProjection_DataProjectionsGet

  PUBLIC DataProjection_DataResultGet

  PUBLIC DataProjection_DecompositionGet

  PUBLIC DataProjection_LabelGet

  PUBLIC DataProjection_MaximumInterationUpdateGet

  PUBLIC DataProjection_MaximumNumberOfIterationsGet

  PUBLIC DataProjection_NumberOfClosestElementsGet

  PUBLIC DataProjection_ProjectionFieldGet

  PUBLIC DataProjection_ProjectionFieldVariableGet

  PUBLIC DataProjection_ProjectionTypeGet

  PUBLIC DataProjection_RelativeToleranceGet

  PUBLIC DataProjection_ResultDistanceGet
  
  PUBLIC DataProjection_ResultElementNumberGet

  PUBLIC DataProjection_ResultElementFaceNumberGet

  PUBLIC DataProjection_ResultElementLineNumberGet
  
  PUBLIC DataProjection_ResultElementXiGet

  PUBLIC DataProjection_ResultExitTagGet

  PUBLIC DataProjection_ResultMaximumErrorGet

  PUBLIC DataProjection_ResultMinimumErrorGet

  PUBLIC DataProjection_ResultProjectionVectorGet
  
  PUBLIC DataProjection_ResultRMSErrorGet

  PUBLIC DataProjection_ResultXiGet

  PUBLIC DataProjection_StartingXiGet

  PUBLIC DataProjection_UserNumberFind
  
  PUBLIC DataProjectionResult_DistanceGet
  
  PUBLIC DataProjectionResult_ElementNumberGet

  PUBLIC DataProjectionResult_ElementFaceNumberGet

  PUBLIC DataProjectionResult_ElementLineNumberGet
  
  PUBLIC DataProjectionResult_ElementXiGet

  PUBLIC DataProjectionResult_ExitTagGet

  PUBLIC DataProjectionResult_ProjectionVectorGet
  
  PUBLIC DataProjectionResult_XiGet

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
#ifdef WITH_POSTCHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
  
    ENTERS("DataProjection_DataPointsGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(dataPoints)) CALL FlagError("Data points is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
#endif    
    
    dataPoints=>dataProjection%dataPoints

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(dataPoints)) THEN
      localError="The data points is not associated for data projection number "// &
        & TRIM(NumberToVString(dataProjection%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
 
    EXITS("DataProjection_DataPointsGet")
    RETURN
999 NULLIFY(dataPoints)
998 ERRORSEXITS("DataProjection_DataPointsGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_DataPointsGet

  !
  !================================================================================================================================
  !  

  !>Gets the user number for a data point identified by a given global number. 
  SUBROUTINE DataProjection_DataPointGlobalNumberGet(dataProjection,userNumber,globalNumber,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection whose data points to get the number for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the data point to get the global number for
    INTEGER(INTG), INTENT(OUT) :: globalNumber !<On exit, the global number of the specified user data point
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DataPointsType), POINTER  :: dataPoints
    TYPE(TreeNodeType), POINTER :: treeNode
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DataProjection_DataPointGlobalNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
#endif
    
    NULLIFY(dataPoints)
    CALL DataProjection_DataPointsGet(dataProjection,dataPoints,err,error,*999)    
    NULLIFY(treeNode)
    CALL Tree_Search(dataPoints%dataPointsTree,userNumber,treeNode,err,error,*999)
    IF(ASSOCIATED(treeNode)) THEN
      CALL Tree_NodeValueGet(dataPoints%dataPointsTree,treeNode,globalNumber,err,error,*999)
    ELSE
      localError="A data point with the user number of "//TRIM(NumberToVString(userNumber,"*",err,error))// &
        & " does not exist in data projection number "//TRIM(NumberToVString(dataProjection%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("DataProjection_DataPointGlobalNumberGet")
    RETURN
999 ERRORSEXITS("DataProjection_DataPointGlobalNumberGet",err,error)    
    RETURN 1
   
  END SUBROUTINE DataProjection_DataPointGlobalNumberGet
  
  !
  !================================================================================================================================
  !

  !>Gets the data projections for a data projection.
  SUBROUTINE DataProjection_DataProjectionsGet(dataProjection,dataProjections,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to get the data projections for
    TYPE(DataProjectionsType), POINTER :: dataProjections !<On exit, the data projections of the data projection. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
   
    ENTERS("DataProjection_DataProjectionsGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(dataProjections)) CALL FlagError("Data projectons is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
#endif    
    
    dataProjections=>dataProjection%dataProjections

#ifdef WITH_POSTCHECKS
    IF(.NOT.ASSOCIATED(dataProjections)) THEN
      localError="The data projections is not associated for data projection number "// &
        & TRIM(NumberToVString(dataProjection%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
 
    EXITS("DataProjection_DataProjectionsGet")
    RETURN
999 NULLIFY(dataProjections)
998 ERRORSEXITS("DataProjection_DataProjectionsGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_DataProjectionsGet

  !
  !================================================================================================================================
  !

  !>Gets a pointer the projection result for a data point identified by a given global number.
  SUBROUTINE DataProjection_DataResultGet(dataProjection,dataPointGlobalNumber,dataProjectionResult,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection for which projection result is stored
    INTEGER(INTG), INTENT(IN) :: dataPointGlobalNumber !<The data projection global number to get the projection result for
    TYPE(DataProjectionResultType), POINTER :: dataProjectionResult !<On exit, the a pointer to the projection result of the specified global data point. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("DataProjection_DataResultGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(dataProjectionResult)) CALL FlagError("Data projection result is already associated.",err,error,*998)
#endif    
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
        & " is invalid for data projection number "//TRIM(NumberToVString(dataProjection%userNumber,"*",err,error))// &
        & ". The global number should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(dataProjection%dataProjectionResults,1),"*",err,error))//"."
    ENDIF
#endif    
     
    dataProjectionResult=>dataProjection%dataProjectionResults(dataPointGlobalNumber)

#ifdef WITH_POSTCHECKS
    IF(.NOT.ASSOCIATED(dataProjectionResult)) THEN
      localError="The data projection result is not associated for the data projection global number "// &
        & TRIM(NumberToVString(dataPointGlobalNumber,"*",err,error))// &
        & " for data projection number "//TRIM(NumberToVString(dataProjection%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    EXITS("DataProjection_DataResultGet")
    RETURN
999 NULLIFY(dataProjectionResult)
998 ERRORSEXITS("DataProjection_DataResultGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_DataResultGet

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
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
   
    ENTERS("DataProjection_DecompositionGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(decomposition)) CALL FlagError("Decomposition is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
#endif    
    
    decomposition=>dataProjection%decomposition

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(decomposition)) THEN
      localError="The decomposition is not associated for data projection number "// &
        & TRIM(NumberToVString(dataProjection%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
 
    EXITS("DataProjection_DecompositionGet")
    RETURN
999 NULLIFY(decomposition)
998 ERRORSEXITS("DataProjection_DecompositionGet",err,error)    
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
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
   
    ENTERS("DataProjection_ProjectionFieldGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(projectionField)) CALL FlagError("Projection field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
#endif    
    
    projectionField=>dataProjection%projectionField

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(projectionField)) THEN
      localError="The projection field is not associated for data projection number "// &
        & TRIM(NumberToVString(dataProjection%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
 
    EXITS("DataProjection_ProjectionFieldGet")
    RETURN
999 NULLIFY(projectionField)
998 ERRORSEXITS("DataProjection_ProjectionFieldGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_ProjectionFieldGet

  !
  !================================================================================================================================
  !

  !>Gets the projection field variable for a data projection.
  SUBROUTINE DataProjection_ProjectionFieldVariableGet(dataProjection,projectionFieldVariable,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to get the projection field variable for
    TYPE(FieldVariableType), POINTER :: projectionFieldVariable !<On exit, the projection field variable of the data projection. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
   
    ENTERS("DataProjection_ProjectionFieldVariableGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(projectionFieldVariable)) CALL FlagError("Projection field variable is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
#endif    
    
    projectionFieldVariable=>dataProjection%projectionVariable

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(projectionFieldVariable)) THEN
      localError="The projection field variable is not associated for data projection number "// &
        & TRIM(NumberToVString(dataProjection%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
 
    EXITS("DataProjection_ProjectionFieldVariableGet")
    RETURN
999 NULLIFY(projectionFieldVariable)
998 ERRORSEXITS("DataProjection_ProjectionFieldVariableGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_ProjectionFieldVariableGet

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
  SUBROUTINE DataProjection_ResultDistanceGet(dataProjection,dataPointGlobalNumber,projectionDistance,err,error,*)

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
    
    ENTERS("DataProjection_ResultDistanceGet",err,error,*999)

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

    EXITS("DataProjection_ResultDistanceGet")
    RETURN
999 ERRORSEXITS("DataProjection_ResultDistanceGet",err,error)
    RETURN 1

  END SUBROUTINE DataProjection_ResultDistanceGet

  !
  !================================================================================================================================
  !

  !>Gets the projection element number for a data point identified by a given global number.
  SUBROUTINE DataProjection_ResultElementNumberGet(dataProjection,dataPointGlobalNumber,projectionElementNumber,err,error,*)

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
   
    ENTERS("DataProjection_ResultElementNumberGet",err,error,*999)

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
 
    EXITS("DataProjection_ResultElementNumberGet")
    RETURN
999 ERRORSEXITS("DataProjection_ResultElementNumberGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_ResultElementNumberGet
  
  !
  !================================================================================================================================
  !

  !>Gets the projection element face number for a data point identified by a given global number.
  SUBROUTINE DataProjection_ResultElementFaceNumberGet(dataProjection,dataPointGlobalNumber,projectionElementFaceNumber &
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
    
    ENTERS("DataProjection_ResultElementFaceNumberGet",err,error,*999)

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
 
    EXITS("DataProjection_ResultElementFaceNumberGet")
    RETURN
999 ERRORSEXITS("DataProjection_ResultElementFaceNumberGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_ResultElementFaceNumberGet

  !
  !================================================================================================================================
  !

  !>Gets the projection element line number for a data point identified by a given global number.
  SUBROUTINE DataProjection_ResultElementLineNumberGet(dataProjection,dataPointGlobalNumber,projectionElementLineNumber &
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
    
    ENTERS("DataProjection_ResultElementLineNumberGet",err,error,*999)

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

    EXITS("DataProjection_ResultElementLineNumberGet")
    RETURN
999 ERRORSEXITS("DataProjection_ResultElementLineNumberGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_ResultElementLineNumberGet

  !
  !================================================================================================================================
  !

  !>Gets the projection element xi for a data point identified by a given global number.
  SUBROUTINE DataProjection_ResultElementXiGet(dataProjection,dataPointGlobalNumber,numberOfElementXi,elementXi,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection for which projection result is stored
    INTEGER(INTG), INTENT(IN) :: dataPointGlobalNumber !<The data projection global number to get the projection element xi for
    INTEGER(INTG), INTENT(OUT) :: numberOfElementXi !<On exit, the number of element xi
    REAL(DP), INTENT(OUT) :: elementXi(:) !<On exit, the projection element xi of the specified global data point
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("DataProjection_ResultXiGet",err,error,*999)

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

    numberOfElementXi=SIZE(dataProjection%dataProjectionResults(dataPointGlobalNumber)%elementXi,1)
    elementXi(1:numberOfElementXi)=dataProjection%dataProjectionResults(dataPointGlobalNumber)%elementXi(1:numberOfElementXi)

    EXITS("DataProjection_ResultElementXiGet")
    RETURN
999 ERRORSEXITS("DataProjection_ResultElementXiGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_ResultElementXiGet

  !
  !================================================================================================================================
  !

  !>Gets the projection exit tag for a data point identified by a given global number.
  SUBROUTINE DataProjection_ResultExitTagGet(dataProjection,dataPointGlobalNumber,projectionExitTag,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection for which projection result is stored
    INTEGER(INTG), INTENT(IN) :: dataPointGlobalNumber !<The data projection global number to get the projection exit tag for
    INTEGER(INTG), INTENT(OUT) :: projectionExitTag !<On exit, the projection exit tag of the specified global data point
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("DataProjection_ResultExitTagGet",err,error,*999)

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

    EXITS("DataProjection_ResultExitTagGet")
    RETURN
999 ERRORSEXITS("DataProjection_ResultExitTagGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_ResultExitTagGet

  !
  !================================================================================================================================
  !

  !>Gets the projection maximum error for a data projection.
  SUBROUTINE DataProjection_ResultMaximumErrorGet(dataProjection,maximumDataPoint,maximumError,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to get the maximum error for
    INTEGER(INTG), INTENT(OUT) :: maximumDataPoint !<On exit, the data point global number of the maximum error.
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
    INTEGER(INTG), INTENT(OUT) :: minimumDataPoint !<On exit, the data point global number of the minimum error.
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
  SUBROUTINE DataProjection_ResultProjectionVectorGet(dataProjection,dataPointGlobalNumber,projectionVector,err,error,*)

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
    
    ENTERS("DataProjection_ResultProjectionVectorGet",err,error,*999)

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

    EXITS("DataProjection_ResultProjectionVectorGet")
    RETURN
999 ERRORSEXITS("DataProjection_ResultProjectionVectorGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_ResultProjectionVectorGet

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
  SUBROUTINE DataProjection_ResultXiGet(dataProjection,dataPointGlobalNumber,numberOfXi,projectionXi,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection for which projection result is stored
    INTEGER(INTG), INTENT(IN) :: dataPointGlobalNumber !<The data projection global number to get the projection xi for
    INTEGER(INTG), INTENT(OUT) :: numberOfXi !<On exit, the number of projection xi
    REAL(DP), INTENT(OUT) :: projectionXi(:) !<On exit, the projection xi of the specified global data point
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("DataProjection_ResultXiGet",err,error,*999)

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

    numberOfXi=dataProjection%numberOfXi
    projectionXi(1:numberOfXi)=dataProjection%dataProjectionResults(dataPointGlobalNumber)%xi(1:numberOfXi)

    EXITS("DataProjection_ResultXiGet")
    RETURN
999 ERRORSEXITS("DataProjection_ResultXiGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_ResultXiGet

  !
  !================================================================================================================================
  !

  !>Gets the starting xi for a data projection.
  SUBROUTINE DataProjection_StartingXiGet(dataProjection,numberOfXi,startingXi,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to get the starting xi for
    INTEGER(INTG), INTENT(OUT) :: numberOfXi !<On exit, the number of xi for the starting xi
    REAL(DP), INTENT(OUT) :: startingXi(:) !<On exit, the starting xi of the specified data projection
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DataProjection_StartingXiGet",err,error,*999)

    CALL DataProjection_AssertIsFinished(dataProjection,err,error,*999)
#ifdef WITH_PRECHECKS    
    IF(SIZE(startingXi,1)<SIZE(dataProjection%startingXi,1)) THEN
      localError="The size of the specified starting xi array of "//TRIM(NumberToVString(SIZE(startingXi,1),"*",err,error))// &
        & " is too small. The size must be >= "//TRIM(NumberToVString(SIZE(dataProjection%startingXi,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    numberOfXi=SIZE(dataProjection%startingXi,1)
    startingXi(1:numberOfXi)=dataProjection%startingXi(1:numberOfXi)
    
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
    TYPE(DataProjectionsType), POINTER :: dataProjections
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("DataProjection_UserNumberFind",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(dataProjection)) CALL FlagError("Data projection is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(dataPoints)) CALL FlagError("Data points is not associated.",err,error,*999)
#endif    

    NULLIFY(dataProjection)
    dataProjections=>dataPoints%dataProjections
    IF(ASSOCIATED(dataProjections)) THEN
      IF(ALLOCATED(dataProjections%dataProjections)) THEN
        projectionIdx=1
        DO WHILE(projectionIdx<=dataProjections%numberOfDataProjections)
          listDataProjection=>dataProjections%dataProjections(projectionIdx)%ptr
#ifdef WITH_PRECHECKS        
          IF(.NOT.ASSOCIATED(listDataProjection)) THEN
            localError="The data points data projections is not associated for projection index "// &
              & TRIM(NumberToVString(projectionIdx,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
#endif        
          IF(listDataProjection%userNumber==userNumber) THEN
            dataProjection=>dataProjections%dataProjections(projectionIdx)%ptr
            EXIT
          ENDIF
        ENDDO
      ENDIF
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

  !>Gets the projection distance for a data projection data point result.
  SUBROUTINE DataProjectionResult_DistanceGet(dataProjectionResult,projectionDistance,err,error,*)

    !Argument variables
    TYPE(DataProjectionResultType), POINTER :: dataProjectionResult !<A pointer to the data projection result to get the distance for
    REAL(DP), INTENT(OUT) :: projectionDistance !<On exit, the projection distance of the specified data point result
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
   
    ENTERS("DataProjectionResult_DistanceGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(dataProjectionResult)) CALL FlagError("Data projection result is not associated.",err,error,*999)
#endif    
    
    projectionDistance=dataProjectionResult%distance

    EXITS("DataProjectionResult_DistanceGet")
    RETURN
999 ERRORSEXITS("DataProjectionResult_DistanceGet",err,error)
    RETURN 1

  END SUBROUTINE DataProjectionResult_DistanceGet

  !
  !================================================================================================================================
  !

  !>Gets the projection element number for a data point result.
  SUBROUTINE DataProjectionResult_ElementNumberGet(dataProjectionResult,projectionElementNumber,err,error,*)

    !Argument variables
    TYPE(DataProjectionResultType), POINTER :: dataProjectionResult !<A pointer to the data projection result to get the projection element number for
    INTEGER(INTG), INTENT(OUT) :: projectionElementNumber !<On exit, the projection local element number of the specified global data point
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
   
    ENTERS("DataProjection_ResultElementNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(dataProjectionResult)) CALL FlagError("Data projection result is not associated.",err,error,*999)
#endif    
    
    projectionElementNumber=dataProjectionResult%elementNumber
 
    EXITS("DataProjectionResult_ElementNumberGet")
    RETURN
999 ERRORSEXITS("DataProjectionResult_ElementNumberGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjectionResult_ElementNumberGet
  
  !
  !================================================================================================================================
  !

  !>Gets the projection element face number for a data point result.
  SUBROUTINE DataProjectionResult_ElementFaceNumberGet(dataProjectionResult,projectionElementFaceNumber,err,error,*)

    !Argument variables
    TYPE(DataProjectionResultType), POINTER :: dataProjectionResult !<A pointer to the data projection result to get the element face number for
    INTEGER(INTG), INTENT(OUT) :: projectionElementFaceNumber !<On exit, the projection element face number of the specified data point result
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("DataProjection_ResultElementFaceNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(dataProjectionResult)) CALL FlagError("Data projection result is not associated.",err,error,*999)
#endif    
    
    projectionElementFaceNumber=dataProjectionResult%elementLineFaceNumber
 
    EXITS("DataProjectionResult_ElementFaceNumberGet")
    RETURN
999 ERRORSEXITS("DataProjectionResult_ElementFaceNumberGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjectionResult_ElementFaceNumberGet

  !
  !================================================================================================================================
  !

  !>Gets the projection element line number for a data point result
  SUBROUTINE DataProjectionResult_ElementLineNumberGet(dataProjectionResult,projectionElementLineNumber,err,error,*)

    !Argument variables
    TYPE(DataProjectionResultType), POINTER :: dataProjectionResult !<A pointer to the data projection result to get the element line number for
    INTEGER(INTG), INTENT(OUT) :: projectionElementLineNumber !<On exit, the projection element line number of the specified data point result
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("DataProjection_ResultElementLineNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(dataProjectionResult)) CALL FlagError("Data projection result is not associated.",err,error,*999)
#endif    
    
    projectionElementLineNumber=dataProjectionResult%elementLineFaceNumber

    EXITS("DataProjectionResult_ElementLineNumberGet")
    RETURN
999 ERRORSEXITS("DataProjectionResult_ElementLineNumberGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjectionResult_ElementLineNumberGet

  !
  !================================================================================================================================
  !

  !>Gets the projection element xi for a data point result.
  SUBROUTINE DataProjectionResult_ElementXiGet(dataProjectionResult,numberOfElementXi,elementXi,err,error,*)

    !Argument variables
    TYPE(DataProjectionResultType), POINTER :: dataProjectionResult !<A pointer to the data projection result to get the element xi for
    INTEGER(INTG), INTENT(OUT) :: numberOfElementXi !<One exit, the number of element xi.
    REAL(DP), INTENT(OUT) :: elementXi(:) !<On exit, the projection element xi of the specified data point result
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("DataProjectionResult_ElementXiGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(dataProjectionResult)) CALL FlagError("Data projection result is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(dataProjectionResult%elementXi)) &
      & CALL FlagError("The element xi array is not allocated for the data projection result.",err,error,*999)    
    IF(SIZE(elementXi,1)<SIZE(dataProjectionResult%elementXi,1)) THEN
      localError="The specified element xi has size of "//TRIM(NumberToVString(SIZE(elementXi,1),"*",err,error))// &
        & " but it needs to have size of >= "//TRIM(NumberToVString(SIZE(dataProjectionResult%elementXi,1),"*",err,error))//"." 
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif

    numberOfElementXi=SIZE(dataProjectionResult%elementXi,1)
    elementXi(1:numberOfElementXi)=dataProjectionResult%elementXi(1:numberOfElementXi)

    EXITS("DataProjectionResult_ElementXiGet")
    RETURN
999 ERRORSEXITS("DataProjectionResult_ElementXiGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjectionResult_ElementXiGet

  !
  !================================================================================================================================
  !

  !>Gets the projection exit tag for a data point result
  SUBROUTINE DataProjectionResult_ExitTagGet(dataProjectionResult,projectionExitTag,err,error,*)

    !Argument variables
    TYPE(DataProjectionResultType), POINTER :: dataProjectionResult !<A pointer to the data projection result to get the exit tag for.
    INTEGER(INTG), INTENT(OUT) :: projectionExitTag !<On exit, the projection exit tag of the specified data point result
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("DataProjectionResult_ExitTagGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(dataProjectionResult)) CALL FlagError("Data projection result is not associated.",err,error,*999)
#endif
         
    projectionExitTag=dataProjectionResult%exitTag

    EXITS("DataProjectionResult_ExitTagGet")
    RETURN
999 ERRORSEXITS("DataProjectionResult_ExitTagGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjectionResult_ExitTagGet

  !
  !================================================================================================================================
  !

  !>Gets the projection vector for a data point result.
  SUBROUTINE DataProjectionResult_ProjectionVectorGet(dataProjectionResult,projectionVector,err,error,*)

    !Argument variables
    TYPE(DataProjectionResultType), POINTER :: dataProjectionResult !<A pointer to the data projection result to get the projection  vector for
    REAL(DP), INTENT(OUT) :: projectionVector(:) !<On exit, the projection vector of the specified data point result
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfCoordinates
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("DataProjectionResult_ProjectionVectorGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(dataProjectionResult)) CALL FlagError("Data projection result is not associated.",err,error,*999)
    IF(SIZE(projectionVector,1)<SIZE(dataProjectionResult%projectionVector,1)) THEN
      localError="The specified projection vector has a size of "// &
        & TRIM(NumberToVString(SIZE(projectionVector,1),"*",err,error))//" but it needs to have size of >= "// &
        & TRIM(NumberToVString(SIZE(dataProjectionResult%projectionVector,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    numberOfCoordinates=SIZE(dataProjectionResult%projectionVector,1)
    projectionVector(1:numberOfCoordinates)=dataProjectionResult%projectionVector(1:numberOfCoordinates)

    EXITS("DataProjectionResult_ProjectionVectorGet")
    RETURN
999 ERRORSEXITS("DataProjectionResult_ProjectionVectorGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjectionResult_ProjectionVectorGet

  !
  !================================================================================================================================
  !

  !>Gets the projection xi for a data point result.
  SUBROUTINE DataProjectionResult_XiGet(dataProjectionResult,numberOfXi,projectionXi,err,error,*)

    !Argument variables
    TYPE(DataProjectionResultType), POINTER :: dataProjectionResult !<A pointer to the data projection result to get the xi for
    INTEGER(INTG), INTENT(OUT) :: numberOfXi !<On exit, the number of projection xi
    REAL(DP), INTENT(OUT) :: projectionXi(:) !<On exit, the projection xi of the specified data point result
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("DataProjectionResult_XiGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(dataProjectionResult)) CALL FlagError("Data projection result is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(dataProjectionResult%xi)) &
      & CALL FlagError("The xi array is not allocated for the data projection result.",err,error,*999)
    IF(SIZE(projectionXi,1)<SIZE(dataProjectionResult%xi,1)) THEN
      localError="The specified projection xi array has size of "//TRIM(NumberToVString(SIZE(projectionXi,1),"*",err,error))// &
        & " but it needs to have size of >= "//TRIM(NumberToVString(SIZE(dataProjectionResult%xi,1),"*",err,error))//"." 
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    numberOfXi=SIZE(dataProjectionResult%xi,1)
    projectionXi(1:numberOfXi)=dataProjectionResult%xi(1:numberOfXi)

    EXITS("DataProjectionResult_XiGet")
    RETURN
999 ERRORSEXITS("DataProjectionResult_XiGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjectionResult_XiGet

  !
  !================================================================================================================================
  !

END MODULE DataProjectionAccessRoutines
