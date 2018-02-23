!> \file
!> \author Tim Wu
!> \brief This module handles all data projection routines
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
!> Contributor(s): Chris Bradley, Kumar Mithraratne, Xiani (Nancy) Yan, Prasad Babarenda Gamage
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

!> This module handles all data projection routines.
MODULE DataProjectionRoutines

  USE BaseRoutines
  USE BasisRoutines
  USE BasisAccessRoutines
  USE CmissMPI  
  USE ComputationEnvironment
  USE Constants
  USE CoordinateSystemAccessRoutines
  USE DataPointAccessRoutines
  USE DataProjectionAccessRoutines
  USE DOMAIN_MAPPINGS
  USE FIELD_ROUTINES
  USE FieldAccessRoutines
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE Kinds
  USE MESH_ROUTINES
  USE MeshAccessRoutines
#ifndef NOMPIMOD
  USE MPI
#endif
  USE Sorting
  USE Strings
  USE Trees
  USE Types

#include "macros.h"

  IMPLICIT NONE

  PRIVATE

#ifdef NOMPIMOD
#include "mpif.h"
#endif

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

  !Module types

  !Module variables

  !Interfaces

  !>Gets the label for a data projection.
  INTERFACE DataProjection_LabelGet
    MODULE PROCEDURE DataProjection_LabelGetC
    MODULE PROCEDURE DataProjection_LabelGetVS
  END INTERFACE DataProjection_LabelGet
  
  !>Sets/changes the label for a data projection.
  INTERFACE DataProjection_LabelSet
    MODULE PROCEDURE DataProjection_LabelSetC
    MODULE PROCEDURE DataProjection_LabelSetVS
  END INTERFACE DataProjection_LabelSet

  !>Cancels the data pojection for data points based on the data point user numbers
  INTERFACE DataProjection_ProjectionCancelByDataPoints
    MODULE PROCEDURE DataProjection_ProjectionCancelByDataPoints0
    MODULE PROCEDURE DataProjection_ProjectionCancelByDataPoints1
  END INTERFACE DataProjection_ProjectionCancelByDataPoints

  !>Cancels the data pojection for data points based on the projection exit tag
  INTERFACE DataProjection_ProjectionCancelByExitTags
    MODULE PROCEDURE DataProjection_ProjectionCancelByExitTags0
    MODULE PROCEDURE DataProjection_ProjectionCancelByExitTags1
  END INTERFACE DataProjection_ProjectionCancelByExitTags

  PUBLIC DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE,DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE, &
    & DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE

  PUBLIC DATA_PROJECTION_DISTANCE_GREATER,DATA_PROJECTION_DISTANCE_GREATER_EQUAL,DATA_PROJECTION_DISTANCE_LESS, &
    & DATA_PROJECTION_DISTANCE_LESS_EQUAL
  
  PUBLIC DataProjection_AbsoluteToleranceGet,DataProjection_AbsoluteToleranceSet

  PUBLIC DataProjection_CreateFinish,DataProjection_CreateStart
  
  PUBLIC DataProjection_Destroy
  
  PUBLIC DataProjection_DataPointsProjectionEvaluate

  PUBLIC DataProjection_DataPointsPositionEvaluate
  
  PUBLIC DataProjection_ElementSet
  
  PUBLIC DataProjection_LabelGet,DataProjection_LabelSet

  PUBLIC DataProjection_MaximumInterationUpdateGet,DataProjection_MaximumInterationUpdateSet

  PUBLIC DataProjection_MaximumNumberOfIterationsGet,DataProjection_MaximumNumberOfIterationsSet

  PUBLIC DataProjection_NumberOfClosestElementsGet,DataProjection_NumberOfClosestElementsSet
  
  PUBLIC DataProjection_ProjectionCancelByDataPoints
  
  PUBLIC DataProjection_ProjectionCancelByDistance
  
  PUBLIC DataProjection_ProjectionCancelByExitTags
  
  PUBLIC DataProjection_ProjectionCandidateElementsSet
  
  PUBLIC DataProjection_ProjectionDataCandidateElementsSet
  
  PUBLIC DataProjection_ProjectionCandidateFacesSet
  
  PUBLIC DataProjection_ProjectionDataCandidateFacesSet
  
  PUBLIC DataProjection_ProjectionCandidateLinesSet
  
  PUBLIC DataProjection_ProjectionDataCandidateLinesSet
  
  PUBLIC DataProjection_ProjectionTypeGet,DataProjection_ProjectionTypeSet
  
  PUBLIC DataProjection_RelativeToleranceGet,DataProjection_RelativeToleranceSet

  PUBLIC DataProjection_ResultAnalysisOutput
  
  PUBLIC DataProjection_ResultDistanceGet
  
  PUBLIC DataProjection_ResultElementNumberGet,DataProjection_ResultElementFaceNumberGet,DataProjection_ResultElementLineNumberGet
  
  PUBLIC DataProjection_ResultExitTagGet

  PUBLIC DataProjection_ResultProjectionVectorGet
  
  PUBLIC DataProjection_ResultXiGet,DataProjection_ResultXiSet

  PUBLIC DataProjection_StartingXiGet,DataProjection_StartingXiSet

  PUBLIC DataProjection_UserNumberFind
  
  PUBLIC DataProjections_Initialise,DataProjections_Finalise
  
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

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionFinished) CALL FlagError("Data projection has not been finished.",err,error,*999)
    
    absoluteTolerance=dataProjection%absoluteTolerance       
     
    EXITS("DataProjection_AbsoluteToleranceGet")
    RETURN
999 ERRORSEXITS("DataProjection_AbsoluteToleranceGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_AbsoluteToleranceGet

  !
  !================================================================================================================================
  !

  !>Sets the absolute tolerance for a data projection.
  SUBROUTINE DataProjection_AbsoluteToleranceSet(dataProjection,absoluteTolerance,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to set the absolute tolerance for
    REAL(DP), INTENT(IN) :: absoluteTolerance !<the absolute tolerance to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("DataProjection_AbsoluteToleranceSet",err,error,*999)

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(dataProjection%dataProjectionFinished) CALL FlagError("Data projection has been finished.",err,error,*999)
    IF(absoluteTolerance<0.0_DP) &
      & CALL FlagError("The specified absolute tolerance is invalid. The tolerance must be > 0.0.",err,error,*999)
      
    dataProjection%absoluteTolerance=absoluteTolerance
    
    EXITS("DataProjection_AbsoluteToleranceSet")
    RETURN
999 ERRORSEXITS("DataProjection_AbsoluteToleranceSet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_AbsoluteToleranceSet

  !
  !================================================================================================================================
  !

  !>Find the closest elements to a data point based on starting xi guess.
  SUBROUTINE DataProjection_ClosestElementsFind(dataProjection,interpolatedPoint,dataPointLocation,numberOfCandidates, &
    & candidateElements,closestElements,closestDistances,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<Data projection problem to evaluate
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: interpolatedPoint !<The interpolated point for the data point
    REAL(DP), INTENT(IN) :: dataPointLocation(:) !<dataPointLocation(componentIdx). The location of data point in the region.
    INTEGER(INTG), INTENT(IN) :: numberOfCandidates !<The number of candidate elements.
    INTEGER(INTG), INTENT(IN) :: candidateElements(:) !<candidateElememnts(candidateIdx). The list of candidate elements for the projection.
    INTEGER(INTG), INTENT(OUT) :: closestElements(:) !<On exit, the list of element numbers with the shortest distances
    REAL(DP), INTENT(OUT) :: closestDistances(:) !<On exit, the list of shortest distances (squared).
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfCoordinates !<Region coordinate dimension
    INTEGER(INTG) :: numberOfClosestCandidates !<Number of closest elements to record
    REAL(DP) :: distanceVector(3) !<distance vector between data point and interpolated point, maximum dimension is 3
    REAL(DP) :: distance2 !<Distance squared
    INTEGER(INTG) :: closestElementIdx
    INTEGER(INTG) :: elementNumber,insertIdx
      
    ENTERS("DataProjection_ClosestElementsFind",err,error,*999)
    
    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionFinished) CALL FlagError("Data projection has not been finished.",err,error,*999)
    
    numberOfCoordinates=dataProjection%numberOfCoordinates
    numberOfClosestCandidates=MIN(numberOfCandidates,SIZE(closestElements,1))
    !loop through the first few elements
    DO closestElementIdx=1,numberOfClosestCandidates
      elementNumber=candidateElements(closestElementIdx)
      CALL Field_InterpolationParametersElementGet(dataProjection%projectionSetType,elementNumber, &
        & interpolatedPoint%INTERPOLATION_PARAMETERS,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
      CALL Field_InterpolateXi(NO_PART_DERIV,dataProjection%startingXi,interpolatedPoint,err,error,*999, &
        & FIELD_GEOMETRIC_COMPONENTS_TYPE)
      distanceVector(1:numberOfCoordinates) = interpolatedPoint%values(:,1)-dataPointLocation
      distance2 = DOT_PRODUCT(distanceVector(1:numberOfCoordinates),distanceVector(1:numberOfCoordinates))
      DO insertIdx=closestElementIdx,1,-1
        IF(insertIdx>1) THEN
          IF(distance2<closestDistances(insertIdx-1)) CYCLE
        ENDIF
        !insert the element number into the correct index
        IF(insertIdx<closestElementIdx) THEN
          closestDistances(insertIdx+1:closestElementIdx)=closestDistances(insertIdx:closestElementIdx-1)
          closestElements(insertIdx+1:closestElementIdx)=closestElements(insertIdx:closestElementIdx-1) 
        ENDIF
        closestDistances(insertIdx)=distance2
        closestElements(insertIdx)=elementNumber  
        EXIT                        
      ENDDO
    ENDDO
    !Loop through the rest of the elements
    DO closestElementIdx=numberOfClosestCandidates+1,numberOfCandidates
      elementNumber=candidateElements(closestElementIdx)
      CALL Field_InterpolationParametersElementGet(dataProjection%projectionSetType,elementNumber, &
        & interpolatedPoint%INTERPOLATION_PARAMETERS,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
      CALL Field_InterpolateXi(NO_PART_DERIV,dataProjection%startingXi,interpolatedPoint,err,error,*999, &
        & FIELD_GEOMETRIC_COMPONENTS_TYPE) 
      distanceVector(1:numberOfCoordinates)=interpolatedPoint%values(:,1)-dataPointLocation
      distance2 = DOT_PRODUCT(distanceVector(1:numberOfCoordinates),distanceVector(1:numberOfCoordinates))
      IF(distance2<closestDistances(numberOfClosestCandidates))THEN
        DO insertIdx=numberOfClosestCandidates,1,-1
          IF(insertIdx>1) THEN
            IF(distance2<closestDistances(insertIdx-1)) CYCLE
          ENDIF
          !insert the element into the correct index
          IF(insertIdx<numberOfClosestCandidates) THEN
            closestDistances(insertIdx+1:numberOfClosestCandidates)=closestDistances(insertIdx:numberOfClosestCandidates-1)
            closestElements(insertIdx+1:numberOfClosestCandidates)=closestElements(insertIdx:numberOfClosestCandidates-1)
          ENDIF
          closestDistances(insertIdx)=distance2
          closestElements(insertIdx)=elementNumber  
          EXIT 
        ENDDO
      ENDIF
    ENDDO
    !closestDistances=SQRT(closestDistances) !return shortest distances
    
    EXITS("DataProjection_ClosestElementsFind")
    RETURN
999 ERRORSEXITS("DataProjection_ClosestElementsFind",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_ClosestElementsFind
  
  !
  !================================================================================================================================
  !

  !>Find the closest faces to a data point base on starting xi guess.
  SUBROUTINE DataProjection_ClosestFacesFind(dataProjection,interpolatedPoint,dataPointLocation,numberOfCandidates, &
    & candidateElements,candidateElementFaces,closestElements,closestElementFaces,closestDistances,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<Data projection problem to evaluate
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: interpolatedPoint !<The interpolated point for the data point
    REAL(DP), INTENT(IN) :: dataPointLocation(:) !<dataPointLocation(componentIdx). The location of data point in the region.
    INTEGER(INTG), INTENT(IN) :: numberOfCandidates !<The number of candidate elements.
    INTEGER(INTG), INTENT(IN) :: candidateElements(:) !<candidateElememnts(candidateIdx). The list of candidate elements for the projection.
    INTEGER(INTG), INTENT(IN) :: candidateElementFaces(:) !<candidateElementFaces(candidateIdx). The list of candidate faces for the projection.
    INTEGER(INTG), INTENT(OUT) :: closestElements(:) !<On exit, the list of element numbers with the shortest distances
    INTEGER(INTG), INTENT(OUT) :: closestElementFaces(:) !<On exit, the list of face numbers with the shortest distances
    REAL(DP), INTENT(OUT) :: closestDistances(:) !<On exit, the list of shortest distances (squared).
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfCoordinates !<Region coordinate dimension
    INTEGER(INTG) :: numberOfClosestCandidates !<Number of closest elements to record
    REAL(DP) :: distanceVector(3) !<distance vector between data point and interpolated point, maximum dimension is 3
    REAL(DP) :: distance2 !<Distance squared
    INTEGER(INTG) :: closestElementIdx
    INTEGER(INTG) :: elementNumber,elementFaceNumber,faceNumber,insertIdx
      
    ENTERS("DataProjection_ClosestFacesFind",err,error,*999)
    
    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionFinished) CALL FlagError("Data projection has not been finished.",err,error,*999)
    
    numberOfCoordinates=dataProjection%numberOfCoordinates
    numberOfClosestCandidates=MIN(numberOfCandidates,SIZE(closestElements,1))
    !loop through the first few faces
    DO closestElementIdx=1,numberOfClosestCandidates
      elementNumber=candidateElements(closestElementIdx)
      elementFaceNumber=candidateElementFaces(closestElementIdx)
      faceNumber=interpolatedPoint%INTERPOLATION_PARAMETERS%field%decomposition%topology%elements%elements( &
        & elementNumber)%ELEMENT_FACES(elementFaceNumber)
      CALL Field_InterpolationParametersFaceGet(dataProjection%projectionSetType,faceNumber,interpolatedPoint% &
        & INTERPOLATION_PARAMETERS,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
      CALL Field_InterpolateXi(NO_PART_DERIV,dataProjection%startingXi,interpolatedPoint,err,error,*999, &
        & FIELD_GEOMETRIC_COMPONENTS_TYPE)
      distanceVector(1:numberOfCoordinates) = interpolatedPoint%values(:,1)-dataPointLocation
      distance2 = DOT_PRODUCT(distanceVector(1:numberOfCoordinates),distanceVector(1:numberOfCoordinates))
      DO insertIdx=closestElementIdx,1,-1
        IF(insertIdx>1) THEN
          IF(distance2<closestDistances(insertIdx-1)) CYCLE
        ENDIF
        !insert the element number into the correct index
        IF(insertIdx<closestElementIdx) THEN
          closestDistances(insertIdx+1:closestElementIdx)=closestDistances(insertIdx:closestElementIdx-1)
          closestElements(insertIdx+1:closestElementIdx)=closestElements(insertIdx:closestElementIdx-1)
          closestElementFaces(insertIdx+1:closestElementIdx)=closestElementFaces(insertIdx:closestElementIdx-1) 
        ENDIF
        closestDistances(insertIdx)=distance2
        closestElements(insertIdx)=elementNumber
        closestElementFaces(insertIdx)=elementFaceNumber
        EXIT                        
      ENDDO
    ENDDO
    !Loop through the rest of the faces
    DO closestElementIdx=numberOfClosestCandidates+1,numberOfCandidates
      elementNumber=candidateElements(closestElementIdx)
      elementFaceNumber=candidateElementFaces(closestElementIdx)
      faceNumber=interpolatedPoint%INTERPOLATION_PARAMETERS%field%decomposition%topology%elements%elements( &
        & elementNumber)%ELEMENT_FACES(elementFaceNumber)          
      CALL Field_InterpolationParametersFaceGet(dataProjection%projectionSetType,faceNumber,interpolatedPoint% &
        & INTERPOLATION_PARAMETERS,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
      CALL Field_InterpolateXi(NO_PART_DERIV,dataProjection%startingXi,interpolatedPoint,err,error,*999, &
        & FIELD_GEOMETRIC_COMPONENTS_TYPE) 
      distanceVector(1:numberOfCoordinates)=interpolatedPoint%values(:,1)-dataPointLocation
      distance2 = DOT_PRODUCT(distanceVector(1:numberOfCoordinates),distanceVector(1:numberOfCoordinates))
      IF(distance2<closestDistances(numberOfClosestCandidates))THEN
        DO insertIdx=numberOfClosestCandidates,1,-1
          IF(insertIdx>1) THEN
            IF(distance2<closestDistances(insertIdx-1)) CYCLE
          ENDIF
          !insert the element into the correct index
          IF(insertIdx<numberOfClosestCandidates) THEN
            closestDistances(insertIdx+1:numberOfClosestCandidates)=closestDistances(insertIdx:numberOfClosestCandidates-1)
            closestElements(insertIdx+1:numberOfClosestCandidates)=closestElements(insertIdx:numberOfClosestCandidates-1)
            closestElementFaces(insertIdx+1:numberOfClosestCandidates)=closestElementFaces(insertIdx: &
              & numberOfClosestCandidates-1)
          ENDIF
          closestDistances(insertIdx)=distance2
          closestElements(insertIdx)=elementNumber
          closestElementFaces(insertIdx)=elementFaceNumber
          EXIT 
        ENDDO
      ENDIF
    ENDDO
    !closestDistances=SQRT(closestDistances) !return shortest distances
    
    EXITS("DataProjection_ClosestFacesFind")
    RETURN
999 ERRORSEXITS("DataProjection_ClosestFacesFind",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_ClosestFacesFind
  
  !
  !================================================================================================================================
  !

  !>Find the closest lines to a data point base on starting xi guess.
  SUBROUTINE DataProjection_ClosestLinesFind(dataProjection,interpolatedPoint,dataPointLocation,numberOfCandidates, &
    & candidateElements,candidateElementLines,closestElements,closestElementLines,closestDistances,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<Data projection problem to evaluate
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: interpolatedPoint !<The interpolated point for the data point 
    REAL(DP), INTENT(IN) :: dataPointLocation(:) !<dataPointLocation(componentIdx). The location of data point in the region.
    INTEGER(INTG), INTENT(IN) :: numberOfCandidates !<The number of candidate elements.
    INTEGER(INTG), INTENT(IN) :: candidateElements(:) !<candidateElememnts(candidateIdx). The list of candidate elements for the projection.
    INTEGER(INTG), INTENT(IN) :: candidateElementLines(:) !<candidateElementLines(candidateIdx). The list of candidate lines for the projection.
    INTEGER(INTG), INTENT(OUT) :: closestElements(:) !<On exit, the list of element numbers with the shortest distances
    INTEGER(INTG), INTENT(OUT) :: closestElementLines(:) !<On exit, the list of lines numbers with the shortest distances
    REAL(DP), INTENT(OUT) :: closestDistances(:) !<On exit, the list of shortest distances (squared).
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfCoordinates !<Region coordinate dimension
    INTEGER(INTG) :: numberOfClosestCandidates !<Number of closest elements to record
    REAL(DP) :: distanceVector(3) !<distance vector between data point and interpolated point, maximum dimension is 3
    REAL(DP) :: distance2 !<Distance squared
    INTEGER(INTG) :: closestElementIdx
    INTEGER(INTG) :: elementNumber,elementLineNumber,lineNumber,insertIdx
      
    ENTERS("DataProjection_ClosestLinesFind",err,error,*999)
    
    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionFinished) CALL FlagError("Data projection has not been finished.",err,error,*999)
    
    numberOfCoordinates=dataProjection%numberOfCoordinates
    numberOfClosestCandidates=MIN(numberOfCandidates,SIZE(closestElements,1))
    !loop through the first few lines
    DO closestElementIdx=1,numberOfClosestCandidates
      elementNumber=candidateElements(closestElementIdx)
      elementLineNumber=candidateElementLines(closestElementIdx)
      lineNumber=interpolatedPoint%INTERPOLATION_PARAMETERS%field%decomposition%topology%elements%elements( &
        & elementNumber)%ELEMENT_LINES(elementLineNumber)
      CALL Field_InterpolationParametersLineGet(dataProjection%projectionSetType,lineNumber,interpolatedPoint% &
        & INTERPOLATION_PARAMETERS,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
      CALL Field_InterpolateXi(NO_PART_DERIV,dataProjection%startingXi,interpolatedPoint,err,error,*999, &
        & FIELD_GEOMETRIC_COMPONENTS_TYPE)
      distanceVector(1:numberOfCoordinates) = interpolatedPoint%values(:,1)-dataPointLocation
      distance2 = DOT_PRODUCT(distanceVector(1:numberOfCoordinates),distanceVector(1:numberOfCoordinates))
      DO insertIdx=closestElementIdx,1,-1
        IF(insertIdx>1) THEN
          IF(distance2<closestDistances(insertIdx-1)) CYCLE
        ENDIF
        !insert the element number into the correct index
        IF(insertIdx<closestElementIdx) THEN
          closestDistances(insertIdx+1:closestElementIdx)=closestDistances(insertIdx:closestElementIdx-1)
          closestElements(insertIdx+1:closestElementIdx)=closestElements(insertIdx:closestElementIdx-1)
          closestElementLines(insertIdx+1:closestElementIdx)=closestElementLines(insertIdx:closestElementIdx-1) 
        ENDIF
        closestDistances(insertIdx)=distance2
        closestElements(insertIdx)=elementNumber
        closestElementLines(insertIdx)=elementLineNumber
        EXIT                        
      ENDDO
    ENDDO
    !Loop through the rest of the lines
    DO closestElementIdx=numberOfClosestCandidates+1,numberOfCandidates
      elementNumber=candidateElements(closestElementIdx)
      elementLineNumber=candidateElementLines(closestElementIdx)
      lineNumber=interpolatedPoint%INTERPOLATION_PARAMETERS%field%decomposition%topology%elements%elements( &
        & elementNumber)%ELEMENT_LINES(elementLineNumber)          
      CALL Field_InterpolationParametersLineGet(dataProjection%projectionSetType,lineNumber,interpolatedPoint% &
        & INTERPOLATION_PARAMETERS,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
      CALL Field_InterpolateXi(NO_PART_DERIV,dataProjection%startingXi,interpolatedPoint,err,error,*999, &
        & FIELD_GEOMETRIC_COMPONENTS_TYPE) 
      distanceVector(1:numberOfCoordinates)=interpolatedPoint%values(:,1)-dataPointLocation 
      distance2 = DOT_PRODUCT(distanceVector(1:numberOfCoordinates),distanceVector(1:numberOfCoordinates))
      IF(distance2<closestDistances(numberOfClosestCandidates))THEN
        DO insertIdx=numberOfClosestCandidates,1,-1
          IF(insertIdx>1) THEN
            IF(distance2<closestDistances(insertIdx-1)) CYCLE
          ENDIF
          !insert the element into the correct index
          IF(insertIdx<numberOfClosestCandidates) THEN
            closestDistances(insertIdx+1:numberOfClosestCandidates)=closestDistances(insertIdx:numberOfClosestCandidates-1)
            closestElements(insertIdx+1:numberOfClosestCandidates)=closestElements(insertIdx:numberOfClosestCandidates-1)
            closestElementLines(insertIdx+1:numberOfClosestCandidates)=closestElementLines(insertIdx: &
              & numberOfClosestCandidates-1)
          ENDIF
          closestDistances(insertIdx)=distance2
          closestElements(insertIdx)=elementNumber
          closestElementLines(insertIdx)=elementLineNumber
          EXIT 
        ENDDO
      ENDIF
    ENDDO
    !closestDistances=SQRT(closestDistances) !return shortest distances
    
    EXITS("DataProjection_ClosestLinesFind")
    RETURN
999 ERRORSEXITS("DataProjection_ClosestLinesFind",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_ClosestLinesFind
  
  !
  !================================================================================================================================
  !

  !>Finishes the process of creating data projection.
  SUBROUTINE DataProjection_CreateFinish(dataProjection,err,error,*)
    
    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<On exit, a pointer to the created data projection
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DataPointsType), POINTER :: dataPoints
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition
    TYPE(VARYING_STRING) :: localError

    ENTERS("DataProjection_CreateFinish",err,error,*999)

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(dataProjection%dataProjectionFinished) CALL FlagError("Data projection has already been finished.",err,error,*999)
    NULLIFY(dataPoints)
    CALL DataProjection_DataPointsGet(dataProjection,dataPoints,err,error,*999)
    IF(.NOT.dataPoints%dataPointsFinished) &
      & CALL FlagError("Data projection data points have not been finished.",err,error,*999)
    NULLIFY(decomposition)
    CALL DataProjection_DecompositionGet(dataProjection,decomposition,err,error,*999)

    !Perform various checks
    SELECT CASE(dataProjection%projectionType)
    CASE(DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE)
      !Check we have lines calculated
      IF(.NOT.decomposition%CALCULATE_LINES) &
        & CALL FlagError("Need to set the decomposition calculate lines to true for a boundary lines projection type.", &
        & err,error,*999)
    CASE(DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE)
      !Check we have faces calculated
      IF(.NOT.decomposition%CALCULATE_FACES) &
        & CALL FlagError("Need to set the decomposition calculate faces to true for a boundary faces projection type.", &
        & err,error,*999)
    CASE(DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE)
      !Do nothing
    CASE DEFAULT
      localError="The data projection type of "//TRIM(NumberToVString(dataProjection%projectionType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    !Finish the data projection
    dataProjection%dataProjectionFinished=.TRUE.
    !Initialise data projection part in data points
    CALL DataProjection_DataProjectionResultsInitialise(dataProjection,err,error,*999) 
    
    EXITS("DataProjection_CreateFinish")
    RETURN
999 ERRORSEXITS("DataProjection_CreateFinish",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_CreateFinish
  
  !
  !================================================================================================================================
  !
  !>Starts the process of creating data projection.
  SUBROUTINE DataProjection_CreateStart(dataProjectionUserNumber,dataPoints,projectionField,projectionVariableType, &
    & dataProjection,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: dataProjectionUserNumber !<The user number of the data projection
    TYPE(DataPointsType), POINTER :: dataPoints !<A pointer to the data points in which to create data projection
    TYPE(FIELD_TYPE), POINTER :: projectionField !<A pointer to the field for the data projection
    INTEGER(INTG), INTENT(IN) :: projectionVariableType !<The field variable type of the projection field for the data projection \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    TYPE(DataProjectionType), POINTER :: dataProjection !<On exit, a pointer to the created data projection. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataProjectionIdx,dummyErr,insertStatus,numberOfCoordinates
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: dataPointsCoordinateSystem,fieldCoordinateSystem
    TYPE(DataProjectionPtrType), ALLOCATABLE :: newDataProjections(:)
    TYPE(DECOMPOSITION_TYPE), POINTER :: fieldDecomposition
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: projectionFieldVariable
    TYPE(VARYING_STRING) :: dummyError,localError
    
    ENTERS("DataProjection_CreateStart",err,error,*998)

    IF(.NOT.ASSOCIATED(dataPoints)) CALL FlagError("Data points is not associated.",err,error,*998)
    IF(.NOT.dataPoints%dataPointsFinished) CALL FlagError("Data points have not been finished.",err,error,*998)
    IF(.NOT.ASSOCIATED(projectionField)) CALL FlagError("Projection field is not associated.",err,error,*998)
    IF(ASSOCIATED(dataProjection)) CALL FlagError("Data projection is already associated.",err,error,*998)
    
    NULLIFY(projectionFieldVariable)
    CALL Field_VariableGet(projectionField,projectionVariableType,projectionFieldVariable,err,error,*999)          
    NULLIFY(dataPointsCoordinateSystem)
    CALL DataPoints_CoordinateSystemGet(dataPoints,dataPointsCoordinateSystem,err,error,*999)
    NULLIFY(fieldCoordinateSystem)
    CALL Field_CoordinateSystemGet(projectionField,fieldCoordinateSystem,err,error,*999)
    NULLIFY(fieldDecomposition)
    CALL Field_DecompositionGet(projectionField,fieldDecomposition,err,error,*999)
    !Check we are on the same coordinate system
    IF(ASSOCIATED(dataPointsCoordinateSystem,fieldCoordinateSystem)) THEN
      CALL CoordinateSystem_DimensionGet(dataPointsCoordinateSystem,numberOfCoordinates,err,error,*999)
      NULLIFY(dataProjection)
      CALL DataProjection_Initialise(dataProjection,err,error,*999)
      dataProjection%userNumber=dataProjectionUserNumber
      dataProjection%dataPoints=>dataPoints
      dataProjection%projectionField=>projectionField
      dataProjection%projectionVariableType=projectionVariableType
      dataProjection%decomposition=>fieldDecomposition
      dataProjection%numberOfCoordinates=numberOfCoordinates
      dataProjection%numberOfElementXi=fieldDecomposition%numberOfDimensions
      !Default always project to boundaries faces/lines when decomposition dimension is equal to region dimension.
      !If decomposition dimension is less, project to all elements            
      IF(fieldDecomposition%numberOfDimensions<numberOfCoordinates) THEN !decomposition dimension < data dimension
        dataProjection%numberOfXi=fieldDecomposition%numberOfDimensions
        dataProjection%projectionType=DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE
      ELSE
        SELECT CASE(fieldDecomposition%numberOfDimensions) !decomposition dimension = data dimension
        CASE(1)
          dataProjection%numberOfXi=1
          dataProjection%projectionType=DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE
        CASE(2) 
          dataProjection%numberOfXi=1
          dataProjection%projectionType=DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE
        CASE(3)
          dataProjection%numberOfXi=2
          dataProjection%projectionType=DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE
        CASE DEFAULT
          localError="The decomposition dimension of "// &
            & TRIM(NumberToVString(fieldDecomposition%numberOfDimensions,"*",err,error))// &
            & " is invalid. The dimension should be >=1 and <= 3."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDIF
      SELECT CASE(dataProjection%numberOfXi) !mesh dimension = data dimension
      CASE(1)
        dataProjection%numberOfClosestElements=2
      CASE(2)
        dataProjection%numberOfClosestElements=4  
      CASE(3)
        dataProjection%numberOfClosestElements=8
      CASE DEFAULT
        localError="The number of xi directions of "// &
          & TRIM(NumberToVString(fieldDecomposition%numberOfDimensions,"*",err,error))// &
          & " is invalid. The number of directions should be >=1 and <= 3."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      ALLOCATE(dataProjection%startingXi(dataProjection%numberOfXi),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate data points data projection starting xi.",err,error,*999)
      dataProjection%startingXi=0.5_DP !<initialised to 0.5 in each xi direction
      CALL DataProjection_DataProjectionCandidatesInitialise(dataProjection,err,error,*999)
      ALLOCATE(newDataProjections(dataPoints%dataProjections%numberOfDataProjections+1),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate new data projections.",err,error,*999)
      IF(dataPoints%dataProjections%numberOfDataProjections>0) THEN
        DO dataProjectionIdx=1,dataPoints%dataProjections%numberOfDataProjections
          newDataProjections(dataProjectionIdx)%ptr=>dataPoints%dataProjections%dataProjections(dataProjectionIdx)%ptr
        ENDDO !xiIdx
      ENDIF
      !Return the pointer
      newDataProjections(dataPoints%dataProjections%numberOfDataProjections+1)%ptr=>dataProjection
      dataPoints%dataProjections%numberOfDataProjections=dataPoints%dataProjections%numberOfDataProjections+1
      dataProjection%globalNumber=dataPoints%dataProjections%numberOfDataProjections
      CALL MOVE_ALLOC(newDataProjections,dataPoints%dataProjections%dataProjections)
    ELSE
      CALL FlagError("The data points and projection field have different coordinate systems.",err,error,*998)        
    ENDIF
     
    EXITS("DataProjection_CreateStart")
    RETURN
999 CALL DataProjection_Finalise(dataProjection,dummyErr,dummyError,*998)
998 IF(ALLOCATED(newDataProjections)) DEALLOCATE(newDataProjections)
    ERRORSEXITS("DataProjection_CreateStart",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_CreateStart 
  
  !
  !================================================================================================================================
  !

  !>Checks that a user data point number is defined for a specific data projection
  SUBROUTINE DataProjection_DataPointCheckExists(dataProjection,dataPointUserNumber,dataPointExists,dataPointGlobalNumber, &
    & err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection whose data points to check
    INTEGER(INTG) :: dataPointUserNumber !<The user data point number to check if it exists
    LOGICAL, INTENT(OUT) :: dataPointExists !<On exit, is .TRUE. if the data point user number exists in the region, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: dataPointGlobalNumber !<On exit, if the data point exists the global number corresponding to the user data point number. If the data point does not exist then global number will be 0.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DataPointsType), POINTER :: dataPoints
    TYPE(TREE_NODE_TYPE), POINTER :: treeNode
   
    ENTERS("DataProjection_DataPointCheckExists",err,ERROR,*999)
    
    dataPointExists=.FALSE.
    dataPointGlobalNumber=0
    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    NULLIFY(dataPoints)
    CALL DataProjection_DataPointsGet(dataProjection,dataPoints,err,error,*999)
    
    NULLIFY(treeNode)
    CALL Tree_Search(dataPoints%dataPointsTree,dataPointUserNumber,treeNode,err,error,*999)
    IF(ASSOCIATED(treeNode)) THEN
      CALL Tree_NodeValueGet(dataPoints%dataPointsTree,treeNode,dataPointGlobalNumber,err,error,*999)
      dataPointExists=.TRUE.
    ENDIF
 
    EXITS("DataProjection_DataPointCheckExists")
    RETURN
999 ERRORSEXITS("DataProjection_DataPointCheckExists",err,error)
    RETURN 1
   
  END SUBROUTINE DataProjection_DataPointCheckExists

  !
  !================================================================================================================================
  !
  
  !>Finalises the data projection candidate and deallocates all memory
  SUBROUTINE DataProjection_DataProjectionCandidateFinalise(dataProjectionCandidate,err,error,*)

    !Argument variables
    TYPE(DataProjectionCandidateType), INTENT(INOUT) :: dataProjectionCandidate !<The data candidate result to finalise 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("DataProjection_DataProjectionCandidateFinalise",err,error,*999)

    IF(ALLOCATED(dataProjectionCandidate%candidateElementNumbers)) DEALLOCATE(dataProjectionCandidate%candidateElementNumbers)
    IF(ALLOCATED(dataProjectionCandidate%localFaceLineNumbers)) DEALLOCATE(dataProjectionCandidate%localFaceLineNumbers)
    
    EXITS("DataProjection_DataProjectionCandidateFinalise")
    RETURN
999 ERRORS("DataProjection_DataProjectionCandidateFinalise",err,error)
    EXITS("DataProjection_DataProjectionCandidateFinalise")
    RETURN 1

  END SUBROUTINE DataProjection_DataProjectionCandidateFinalise
  
  !
  !================================================================================================================================
  !
  
  !>Initialises a data projection candidate
  SUBROUTINE DataProjection_DataProjectionCandidateInitialise(dataProjectionCandidate,err,error,*)

    !Argument variables
    TYPE(DataProjectionCandidateType) :: dataProjectionCandidate !<The data projection candidate to initialise 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("DataProjection_DataProjectionCandidateInitialise",err,error,*999)

    IF(ALLOCATED(dataProjectionCandidate%candidateElementNumbers)) DEALLOCATE(dataProjectionCandidate%candidateElementNumbers)
    IF(ALLOCATED(dataProjectionCandidate%localFaceLineNumbers)) DEALLOCATE(dataProjectionCandidate%localFaceLineNumbers)
    
    EXITS("DataProjection_DataProjectionCandidateInitialise")
    RETURN
999 ERRORS("DataProjection_DataProjectionCandidateInitialise",err,error)
    EXITS("DataProjection_DataProjectionCandidateInitialise")
    RETURN 1

  END SUBROUTINE DataProjection_DataProjectionCandidateInitialise
  
  !
  !================================================================================================================================
  !

  !>Finalise data projection candidates.
  SUBROUTINE DataProjection_DataProjectionCandidatesFinalise(dataProjection,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to finalise the data projection candidates for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointIdx
       
    ENTERS("DataProjection_DataProjectionCandidatesFinalise",err,error,*999)

    IF(ASSOCIATED(dataProjection)) THEN
      IF(ALLOCATED(dataProjection%dataProjectionCandidates)) THEN
        DO dataPointIdx=0,SIZE(dataProjection%dataProjectionCandidates,1)-1
          CALL DataProjection_DataProjectionCandidateFinalise(dataProjection%dataProjectionCandidates(dataPointIdx),err,error,*999)
        ENDDO !dataPointIdx
        DEALLOCATE(dataProjection%dataProjectionCandidates)
      ENDIF
    ENDIF
    
    EXITS("DataProjection_DataProjectionCandidatesFinalise")
    RETURN
999 ERRORS("DataProjection_DataProjectionCandidatesFinalise",err,error)
    EXITS("DataProjection_DataProjectionCandidatesFinalise")
    RETURN 1
    
  END SUBROUTINE DataProjection_DataProjectionCandidatesFinalise
  
  !
  !================================================================================================================================
  !
  
  !>Initialises the data projection candidates.
  SUBROUTINE DataProjection_DataProjectionCandidatesInitialise(dataProjection,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to initialise the projection candidates for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointIdx,dummyErr,numberOfDataPoints
    TYPE(DataPointsType), POINTER :: dataPoints
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("DataProjection_DataProjectionCandidatesInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*998)
    NULLIFY(dataPoints)
    CALL DataProjection_DataPointsGet(dataProjection,dataPoints,err,error,*998)
    
    numberOfDataPoints = dataPoints%numberOfDataPoints
    ALLOCATE(dataProjection%dataProjectionCandidates(0:numberOfDataPoints),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate data projection data projection candidates.",err,error,*999)
    DO dataPointIdx=1,numberOfDataPoints
      CALL DataProjection_DataProjectionCandidateInitialise(dataProjection%dataProjectionCandidates(dataPointIdx),err,error,*999)
    ENDDO !dataPointIdx
    
    EXITS("DataProjection_DataProjectionCandidatesInitialise")
    RETURN
999 CALL DataProjection_DataProjectionCandidatesFinalise(dataProjection,dummyErr,dummyError,*998)
998 ERRORS("DataProjection_DataProjectionCandidatesInitialise",err,error)
    EXITS("DataProjection_DataProjectionCandidatesInitialise")
    RETURN 1

  END SUBROUTINE DataProjection_DataProjectionCandidatesInitialise
  
  !
  !================================================================================================================================
  !
  
  !>Finalises the data projection result and deallocates all memory
  SUBROUTINE DataProjection_DataProjectionResultFinalise(dataProjectionResult,err,error,*)

    !Argument variables
    TYPE(DataProjectionResultType), INTENT(INOUT) :: dataProjectionResult !<The data projection result to finalise 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("DataProjection_DataProjectionResultFinalise",err,error,*999)

    dataProjectionResult%userNumber=0
    dataProjectionResult%distance=0.0
    dataProjectionResult%elementNumber=0
    dataProjectionResult%elementLineFaceNumber=0
    dataProjectionResult%exitTag=0
    IF(ALLOCATED(dataProjectionResult%xi)) DEALLOCATE(dataProjectionResult%xi)
    IF(ALLOCATED(dataProjectionResult%elementXi)) DEALLOCATE(dataProjectionResult%elementXi)
    IF(ALLOCATED(dataProjectionResult%projectionVector)) DEALLOCATE(dataProjectionResult%projectionVector)
    
    EXITS("DataProjection_DataProjectionResultFinalise")
    RETURN
999 ERRORS("DataProjection_DataProjectionResultFinalise",err,error)
    EXITS("DataProjection_DataProjectionResultFinalise")
    RETURN 1

  END SUBROUTINE DataProjection_DataProjectionResultFinalise
  
  !
  !================================================================================================================================
  !
  
  !>Initialises a data projection result
  SUBROUTINE DataProjection_DataProjectionResultInitialise(dataProjectionResult,err,error,*)

    !Argument variables
    TYPE(DataProjectionResultType) :: dataProjectionResult !<The data projection result to initialise 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("DataProjection_DataProjectionResultInitialise",err,error,*999)

    dataProjectionResult%userNumber=0
    dataProjectionResult%distance=0.0_DP
    dataProjectionResult%elementNumber=0
    dataProjectionResult%elementLineFaceNumber=0
    dataProjectionResult%exitTag=DATA_PROJECTION_EXIT_TAG_NO_ELEMENT
    
    EXITS("DataProjection_DataProjectionResultInitialise")
    RETURN
999 ERRORS("DataProjection_DataProjectionResultInitialise",err,error)
    EXITS("DataProjection_DataProjectionResultInitialise")
    RETURN 1

  END SUBROUTINE DataProjection_DataProjectionResultInitialise 
  
  !
  !================================================================================================================================
  !

  !>Finalise a data projection.
  SUBROUTINE DataProjection_DataProjectionResultsFinalise(dataProjection,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to finalise the data projection results for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointIdx
       
    ENTERS("DataProjection_DataProjectionResultsFinalise",err,error,*999)

    IF(ASSOCIATED(dataProjection)) THEN
      IF(ALLOCATED(dataProjection%dataProjectionResults)) THEN
        DO dataPointIdx=1,SIZE(dataProjection%dataProjectionResults,1)
          CALL DataProjection_DataProjectionResultFinalise(dataProjection%dataProjectionResults(dataPointIdx),err,error,*999)
        ENDDO !dataPointIdx
        DEALLOCATE(dataProjection%dataProjectionResults)
      ENDIF
    ENDIF
    
    EXITS("DataProjection_DataProjectionResultsFinalise")
    RETURN
999 ERRORS("DataProjection_DataProjectionResultsFinalise",err,error)
    EXITS("DataProjection_DataProjectionResultsFinalise")
    RETURN 1
    
  END SUBROUTINE DataProjection_DataProjectionResultsFinalise
  
  !
  !================================================================================================================================
  !
  
  !>Initialises the data projection part in a given data points.
  SUBROUTINE DataProjection_DataProjectionResultsInitialise(dataProjection,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to initialise the data projection result for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointIdx,dummyErr,numberOfDataPoints
    TYPE(DataPointsType), POINTER :: dataPoints
    TYPE(VARYING_STRING) :: dummyError,localError
    
    ENTERS("DataProjection_DataProjectionResultsInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*998)
    IF(.NOT.dataProjection%dataProjectionFinished) CALL FlagError("Data projection has not been finished.",err,error,*998)
    NULLIFY(dataPoints)
    CALL DataProjection_DataPointsGet(dataProjection,dataPoints,err,error,*999)
    IF(.NOT.dataPoints%dataPointsFinished) CALL FlagError("Data projection data points have not been finished.",err,error,*998)
    IF(.NOT.ALLOCATED(dataPoints%dataPoints)) CALL FlagError("Data points data points have not been allocated.",err,error,*999)
    
    numberOfDataPoints = dataPoints%numberOfDataPoints
    ALLOCATE(dataProjection%dataProjectionResults(numberOfDataPoints),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate data projection data projection results.",err,error,*999)
    DO dataPointIdx=1,numberOfDataPoints
      CALL DataProjection_DataProjectionResultInitialise(dataProjection%dataProjectionResults(dataPointIdx),err,error,*999)
      dataProjection%dataProjectionResults(dataPointIdx)%userNumber=dataPoints%dataPoints(dataPointIdx)%userNumber
      ALLOCATE(dataProjection%dataProjectionResults(dataPointIdx)%xi(dataProjection%numberOfXi),STAT=err)
      IF(err/=0) THEN
        localError="Could not allocate data projection results xi for data point index "// &
          & TRIM(NumberToVString(dataPointIdx,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      ALLOCATE(dataProjection%dataProjectionResults(dataPointIdx)%elementXi(dataProjection%numberOfElementXi),STAT=err)
      IF(err/=0) THEN
        localError="Could not allocate data projection results element xi for data point index "// &
          & TRIM(NumberToVString(dataPointIdx,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      ALLOCATE(dataProjection%dataProjectionResults(dataPointIdx)%projectionVector(dataProjection%numberOfCoordinates), &
        & STAT=err)
      IF(err/=0) THEN
        localError="Could not allocate data projection results projection vector for data point index "// &
          & TRIM(NumberToVString(dataPointIdx,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      dataProjection%dataProjectionResults(dataPointIdx)%xi(1:dataProjection%numberOfXi)= &
        & dataProjection%startingXi(1:dataProjection%numberOfXi)
      dataProjection%dataProjectionResults(dataPointIdx)%elementXi=0.0_DP
      dataProjection%dataProjectionResults(dataPointIdx)%projectionVector(1:dataProjection%numberOfCoordinates)=0.0_DP
    ENDDO !dataPointIdx

    EXITS("DataProjection_DataProjectionResultsInitialise")
    RETURN
999 CALL DataProjection_DataProjectionResultsFinalise(dataProjection,dummyErr,dummyError,*998)
998 ERRORS("DataProjection_DataProjectionResultsInitialise",err,error)
    EXITS("DataProjection_DataProjectionResultsInitialise")
    RETURN 1

  END SUBROUTINE DataProjection_DataProjectionResultsInitialise 
  
  !
  !================================================================================================================================
  !

  !>Destroys a data projection.
  SUBROUTINE DataProjection_Destroy(dataProjection,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: position,projectionIdx
    LOGICAL :: found
    TYPE(DataProjectionType), POINTER :: listDataProjection
    TYPE(DataProjectionPtrType), ALLOCATABLE :: newDataProjections(:)
    TYPE(DataProjectionsType), POINTER :: dataProjections
    TYPE(VARYING_STRING) :: localError
       
    ENTERS("DataProjection_Destroy",err,error,*999)

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    dataProjections=>dataProjection%dataProjections
    IF(.NOT.ASSOCIATED(dataProjections)) CALL FlagError("Data projection data projections is not associated.",err,error,*999)
    
    !Find the data projection in the list of data projections
    found=.FALSE.
    position=0
    DO WHILE(position<dataProjections%numberOfDataProjections.AND..NOT.found)
      position=position+1
      listDataProjection=>dataProjections%dataProjections(position)%ptr
      IF(ASSOCIATED(listDataProjection)) THEN
        IF(dataProjection%userNumber==listDataProjection%userNumber) THEN
          found=.TRUE.
          EXIT
        ENDIF
      ELSE
        localError="The data projection is not associated for data projections position "// &
          & TRIM(NumberToVString(position,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDDO
    IF(found) THEN
      !The data projection to destroy has been found. Finalise this projection.
      CALL DataProjection_Finalise(dataProjection,err,error,*999)
      !Remove the data projection from the list of data projections
      IF(dataProjections%numberOfDataProjections>1) THEN
        ALLOCATE(newDataProjections(dataProjections%numberOfDataProjections-1),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate new data projections.",err,error,*999)
        DO projectionIdx=1,dataProjections%numberOfDataProjections
          IF(projectionIdx<position) THEN
            newDataProjections(projectionIdx)%ptr=>dataProjections%dataProjections(projectionIdx)%ptr
          ELSE IF(projectionIdx>position) THEN
            IF(ASSOCIATED(dataProjections%dataProjections(projectionIdx)%ptr)) THEN
              dataProjections%dataProjections(projectionIdx)%ptr%globalNumber= &
                & dataProjections%dataProjections(projectionIdx)%ptr%globalNumber-1
            ELSE
              localError="The data projection is not associated for data projection index "// &
                & TRIM(NumberToVString(projectionIdx,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
            newDataProjections(projectionIdx-1)%ptr=>dataProjections%dataProjections(projectionIdx)%ptr
          ENDIF
        ENDDO !projectionIdx
        CALL MOVE_ALLOC(newDataProjections,dataProjections%dataProjections)
        dataProjections%numberOfDataProjections=dataProjections%numberOfDataProjections-1
      ELSE
        DEALLOCATE(dataProjections%dataProjections)
        dataProjections%numberOfDataProjections=0
      ENDIF
    ELSE
      localError="The data projection with user number "//TRIM(NumberToVString(dataProjection%userNumber,"*",err,error))// &
        & " could not be found."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("DataProjection_Destroy")
    RETURN
999 IF(ALLOCATED(newDataProjections)) DEALLOCATE(newDataProjections)
    ERRORSEXITS("DataProjection_Destroy",err,error)
    RETURN 1
    
  END SUBROUTINE DataProjection_Destroy
  
  !
  !================================================================================================================================
  !

  !>Calculate the errors for a data projection.
  SUBROUTINE DataProjection_ErrorsCalculate(dataProjection,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to calculate the errors for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointIdx,exitTag,maxDataPoint,minDataPoint,numberOfDataPoints
    REAL(DP) :: distance,maxError,minError,rmsError
    TYPE(DataPointsType), POINTER :: dataPoints
       
    ENTERS("DataProjection_ErrorsCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(dataProjection%dataProjectionResults)) &
      & CALL FlagError("Data projection projection results is not allocated.",err,error,*999)

    NULLIFY(dataPoints)
    CALL DataProjection_DataPointsGet(dataProjection,dataPoints,err,error,*999)
    !Evaluate errors
    rmsError=0.0_DP
    maxError=0.0_DP
    minError=HUGE(1.0_DP)
    maxDataPoint=0
    minDataPoint=0
    numberOfDataPoints=0
    DO dataPointIdx=1,dataPoints%numberOfDataPoints
      distance=dataProjection%dataProjectionResults(dataPointIdx)%distance
      exitTag=dataProjection%dataProjectionResults(dataPointIdx)%exitTag
      IF(exitTag/=DATA_PROJECTION_CANCELLED) THEN
        numberOfDataPoints=numberOfDataPoints+1
        rmsError=rmsError+distance*distance
        IF(distance>maxError) THEN
          maxError=distance
          maxDataPoint=dataPointIdx
        ENDIF
        IF(distance<minError) THEN
          minError=distance
          minDataPoint=dataPointIdx
        ENDIF
      ENDIF
    ENDDO !dataPointIdx
    IF(numberOfDataPoints>0) THEN
      rmsError=SQRT(rmsError/REAL(numberOfDataPoints,DP))
    ELSE
      rmsError=SQRT(rmsError)
    ENDIF
    dataProjection%rmsError=rmsError
    dataProjection%maximumError=maxError
    dataProjection%maximumErrorDataPoint=maxDataPoint
    dataProjection%minimumError=minError
    dataProjection%minimumErrorDataPoint=minDataPoint
    
    EXITS("DataProjection_ErrorsCalculate")
    RETURN
999 ERRORSEXITS("DataProjection_ErrorsCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE DataProjection_ErrorsCalculate

  !
  !================================================================================================================================
  !

  !>Finalise a data projection.
  SUBROUTINE DataProjection_Finalise(dataProjection,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
       
    ENTERS("DataProjection_Finalise",err,error,*999)

    IF(ASSOCIATED(dataProjection)) THEN
      IF(ALLOCATED(dataProjection%startingXi)) DEALLOCATE(dataProjection%startingXi)
      CALL DataProjection_DataProjectionCandidatesFinalise(dataProjection,err,error,*999)
      CALL DataProjection_DataProjectionResultsFinalise(dataProjection,err,error,*999)
      DEALLOCATE(dataProjection)
    ENDIF
    
    EXITS("DataProjection_Finalise")
    RETURN
999 ERRORSEXITS("DataProjection_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE DataProjection_Finalise
  
  !
  !================================================================================================================================
  !

  !>Initialise a data projection.
  SUBROUTINE DataProjection_Initialise(dataProjection,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
       
    ENTERS("DataProjection_Initialise",err,error,*998)

    IF(ASSOCIATED(dataProjection)) CALL FlagError("Data projection is already associated.",err,error,*998)

    ALLOCATE(dataProjection,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate data projection.",err,error,*999)
    dataProjection%globalNumber=0
    dataProjection%userNumber=0
    dataProjection%label=""
    dataProjection%dataProjectionFinished=.FALSE.
    dataProjection%dataProjectionProjected=.FALSE.
    NULLIFY(dataProjection%dataPoints)
    NULLIFY(dataProjection%dataProjections)
    NULLIFY(dataProjection%projectionField)
    dataProjection%projectionVariableType=0
    NULLIFY(dataProjection%decomposition)
    dataProjection%numberOfCoordinates=0
    dataProjection%numberOfXi=0
    dataProjection%maxNumberOfCandidates=0
    dataProjection%projectionType=0
    dataProjection%maximumIterationUpdate=0.5_DP
    dataProjection%maximumNumberOfIterations=25
    dataProjection%numberOfClosestElements=0
    dataProjection%absoluteTolerance=1.0E-8_DP
    dataProjection%relativeTolerance=1.0E-6_DP
    dataProjection%rmsError=0.0_DP
    dataProjection%maximumError=0.0_DP
    dataProjection%maximumErrorDataPoint=0
    dataProjection%minimumError=0.0_DP
    dataProjection%minimumErrorDataPoint=0
    
    EXITS("DataProjection_Initialise")
    RETURN
999 CALL DataProjection_Finalise(dataProjection,dummyErr,dummyError,*998)
998 ERRORSEXITS("DataProjection_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE DataProjection_Initialise
  
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
    TYPE(TREE_NODE_TYPE), POINTER :: treeNode
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DataProjection_DataPointGlobalNumberGet",err,error,*999)

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    NULLIFY(dataPoints)
    CALL DataProjection_DataPointsGet(dataProjection,dataPoints,err,error,*999)
    
    NULLIFY(treeNode)
    CALL Tree_Search(dataPoints%dataPointsTree,userNumber,treeNode,err,error,*999)
    IF(ASSOCIATED(treeNode)) THEN
      CALL Tree_NodeValueGet(dataPoints%dataPointsTree,treeNode,globalNumber,err,error,*999)
    ELSE
      localError="A data point with the user number of "//TRIM(NumberToVString(userNumber,"*",err,error))//" does not exist."
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

  !>Evaluates a data projection.
  SUBROUTINE DataProjection_DataPointsProjectionEvaluate(dataProjection,projectionFieldSetType,err,error,*) 
    
    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to evaluate
    INTEGER(INTG) :: projectionFieldSetType !<The parameter set type to evaluate the data projection for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(BASIS_TYPE), POINTER :: basis
    TYPE(FIELD_INTERPOLATION_PARAMETERS_PTR_TYPE), POINTER :: interpolationParameters(:)
    TYPE(FIELD_INTERPOLATED_POINT_PTR_TYPE), POINTER :: interpolatedPoints(:)
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: interpolatedPoint
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: parameterSet
    TYPE(DataPointsType), POINTER :: dataPoints
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition
    TYPE(DECOMPOSITION_ELEMENTS_TYPE), POINTER :: decompositionElements
    TYPE(DECOMPOSITION_FACES_TYPE), POINTER :: decompositionFaces
    TYPE(DECOMPOSITION_LINES_TYPE), POINTER :: decompositionLines
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: decompositionTopology    
    TYPE(DOMAIN_TYPE), POINTER :: domain
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: domainElements
    TYPE(DOMAIN_FACES_TYPE), POINTER :: domainFaces
    TYPE(DOMAIN_LINES_TYPE), POINTER :: domainLines
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: domainMappingElements
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: domainMappings
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: domainTopology
    INTEGER(INTG) :: myComputationalNode,numberOfComputationalNodes !<computational node/rank of the current process    
    INTEGER(INTG) :: numberOfDataPoints
    INTEGER(INTG) :: numberOfElements,numberOfFaces,numberOfLines,numberOfCandidates,dataNumberOfCandidates
    INTEGER(INTG) :: numberOfClosestCandidates,totalNumberOfClosestCandidates,reducedNumberOfCLosestCandidates
    INTEGER(INTG), ALLOCATABLE :: globalToLocalNumberOfClosestCandidates(:) 
    INTEGER(INTG), ALLOCATABLE :: candidateElements(:),closestElements(:,:),candidateLinesFaces(:),closestLinesFaces(:,:)
    INTEGER(INTG), ALLOCATABLE :: globalNumberOfClosestCandidates(:)
    INTEGER(INTG), ALLOCATABLE :: globalMPIDisplacements(:),sortingIndices1(:),sortingIndices2(:)
    INTEGER(INTG), ALLOCATABLE :: globalNumberOfProjectedPoints(:)
    INTEGER(INTG) :: MPIClosestDistances,dataProjectionGlobalNumber
    INTEGER(INTG) :: MPIIError
    INTEGER(INTG), ALLOCATABLE :: projectedElement(:),projectedLineFace(:),projectionExitTag(:)
    REAL(DP), ALLOCATABLE :: closestDistances(:,:),globalClosestDistances(:,:)
    REAL(DP), ALLOCATABLE :: projectedDistance(:,:),projectedXi(:,:),projectionVectors(:,:)   
    INTEGER(INTG) :: elementIdx,elementNumber,lineFaceIdx,lineFaceNumber,computationalNodeIdx,xiIdx,localElementNumber, &
      & localLineFaceNumber
    INTEGER(INTG) :: startIdx,finishIdx,dataPointIdx
    LOGICAL :: boundaryProjection,elementExists,ghostElement
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DataProjection_DataPointsProjectionEvaluate",err,error,*999)
    
    NULLIFY(interpolationParameters)
    NULLIFY(interpolatedPoints)
    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionFinished) CALL FlagError("Data projection has not been finished.",err,error,*999)      
    NULLIFY(parameterSet)
    CALL Field_ParameterSetGet(dataProjection%projectionField,dataProjection%projectionVariableType, &
      & projectionFieldSetType,parameterSet,err,error,*999)
    
    dataProjection%projectionSetType=projectionFieldSetType
    dataProjectionGlobalNumber=dataProjection%globalNumber
    dataPoints=>dataProjection%dataPoints
    CALL Field_InterpolationParametersInitialise(dataProjection%projectionField,interpolationParameters, &
      & err,error,*998,FIELD_GEOMETRIC_COMPONENTS_TYPE)
    CALL Field_InterpolatedPointsInitialise(interpolationParameters,interpolatedPoints,err,error,*998, &
      & FIELD_GEOMETRIC_COMPONENTS_TYPE)
    interpolatedPoint=>interpolatedPoints(dataProjection%projectionVariableType)%ptr
    decomposition=>dataProjection%decomposition
    NULLIFY(decompositionTopology)
    CALL Decomposition_TopologyGet(decomposition,decompositionTopology,err,error,*999)
    NULLIFY(domain)
    CALL Decomposition_DomainGet(decomposition,0,domain,err,error,*999)
    NULLIFY(domainTopology)
    CALL Domain_TopologyGet(domain,domainTopology,err,error,*999)
    NULLIFY(domainElements)
    CALL DomainTopology_ElementsGet(domainTopology,domainElements,err,error,*999)
    NULLIFY(domainFaces)
    CALL DomainTopology_FacesGet(domainTopology,domainFaces,err,error,*999)
    NULLIFY(domainLines)
    CALL DomainTopology_LinesGet(domainTopology,domainLines,err,error,*999)
    NULLIFY(domainMappings)
    CALL Domain_MappingsGet(domain,domainMappings,err,error,*999)
    NULLIFY(domainMappingElements)
    CALL DomainMappings_ElementsGet(domainMappings,domainMappingElements,err,error,*999)
    IF(.NOT.ASSOCIATED(domainElements%elements)) CALL FlagError("Domain elements elements is not associated.",err,error,*999)
    
    numberOfDataPoints=dataPoints%numberOfDataPoints
    myComputationalNode=ComputationalEnvironment_NodeNumberGet(err,error)
    IF(err/=0) GOTO 999
    numberOfComputationalNodes=ComputationalEnvironment_NumberOfNodesGet(err,error)
    IF(err/=0) GOTO 999
    boundaryProjection=(dataProjection%projectionType==DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE).OR. &
      & (dataProjection%projectionType==DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE)          
    !#########################################################################################################
    !Find elements/faces/lines inside the current computational node, get the boundary faces/lines only if
    !asked the elements/faces/lines are required to perform projection of points in the current computational
    !node the are all pre-allocated to the maximum array length (i.e., numberOfElements), but only up to the
    !numberOfCandidates'th index are assigned
    numberOfElements=domainElements%NUMBER_OF_ELEMENTS
    numberOfFaces=domainFaces%NUMBER_OF_FACES
    numberOfLines=domainLines%NUMBER_OF_LINES        
    numberOfCandidates=0
    IF(ALLOCATED(dataProjection%dataProjectionCandidates(0)%candidateElementNumbers).AND. &
      & ALLOCATED(dataProjection%dataProjectionCandidates(0)%localFaceLineNumbers)) THEN
      !The global candidates have been set.
      SELECT CASE(dataProjection%projectionType)
      CASE (DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE)
        !Identify all non-ghost boundary lines
        ALLOCATE(candidateElements(numberOfLines),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate candidate elements.",err,error,*999)
        ALLOCATE(candidateLinesFaces(numberOfLines),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate candidate lines faces.",err,error,*999)
        !Loop through all candidate element defined by user number
        DO elementIdx=1,SIZE(dataProjection%dataProjectionCandidates(0)%candidateElementNumbers,1)
          !Check if element exists on current domain, get local number
          CALL DecompositionTopology_ElementCheckExists(decompositionTopology,dataProjection% &
            & dataProjectionCandidates(0)%candidateElementNumbers(elementIdx),elementExists, &
            & localElementNumber,ghostElement,err,error,*999) 
          IF((elementExists).AND.(.NOT.ghostElement)) THEN
            !Get non-ghost elements
            numberOfCandidates=numberOfCandidates+1
            candidateElements(numberOfCandidates)=localElementNumber
            !Store element line number for line projection type                      
            candidateLinesFaces(numberOfCandidates)=dataProjection%dataProjectionCandidates(0)% &
              & localFaceLineNumbers(elementIdx) 
          ENDIF
        ENDDO !elementIdx
      CASE (DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE)
        !identify all non-ghost boundary faces
        IF(decomposition%numberOfDimensions>=2) THEN
          ALLOCATE(candidateElements(numberOfFaces),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate candidate elements.",err,error,*999)
          ALLOCATE(candidateLinesFaces(numberOfFaces),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate candidate lines faces.",err,error,*999)
          !Loop through all candidate elements defined by user number
          DO elementIdx=1,SIZE(dataProjection%dataProjectionCandidates(0)%candidateElementNumbers,1)
            !Check if element exists on current domain, get local number                      
            CALL DecompositionTopology_ElementCheckExists(decompositionTopology,dataProjection% &
              & dataProjectionCandidates(0)%candidateElementNumbers(elementIdx),elementExists, &
              & localElementNumber,ghostElement, &
              & err,error,*999)
            IF((elementExists).AND.(.NOT.ghostElement)) THEN
              !Get non-ghost elements
              numberOfCandidates=numberOfCandidates+1
              candidateElements(numberOfCandidates)=localElementNumber
              !Store element face number for face projection type
              candidateLinesFaces(numberOfCandidates)=dataProjection%dataProjectionCandidates(0)% &
                & localFaceLineNumbers(elementIdx)
            ENDIF
          ENDDO !elementIdx
        ELSE
          CALL FlagError("Decomposition number of dimensions has to be >= 2 for a faces projection type.",err,error,*999)        
        ENDIF
      CASE (DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE)
        !Identify all non-ghost elements
        IF(dataProjection%numberOfXi==decomposition%numberOfDimensions) THEN
          ALLOCATE(candidateElements(numberOfElements),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate candidate elements.",err,error,*999)
          !Loop through all candidate elements defined by user number                    
          DO elementIdx=1,SIZE(dataProjection%dataProjectionCandidates(0)%candidateElementNumbers,1)
            !Check if element exists on current domain, get local number
            CALL DecompositionTopology_ElementCheckExists(decompositionTopology,dataProjection% &
              & dataProjectionCandidates(0)%candidateElementNumbers(elementIdx),elementExists, &
              & localElementNumber,ghostElement,err,error,*999)
            IF((elementExists).AND.(.NOT.ghostElement)) THEN
              !Get non-ghost elements
              numberOfCandidates=numberOfCandidates+1
              candidateElements(numberOfCandidates)=localElementNumber
            ENDIF
          ENDDO !elementIdx
        ELSE
          CALL FlagError("Data projection number of xi has to equal to decomposition number of dimensions",err,error,*999)
        ENDIF
      CASE DEFAULT
        localError="The data projection type of "//TRIM(NumberToVString(dataProjection%projectionType,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      !If user didn't define candidate element number
      NULLIFY(decompositionElements)
      CALL DecompositionTopology_ElementsGet(decompositionTopology,decompositionElements,err,error,*999)
      SELECT CASE(dataProjection%projectionType)
      CASE(DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE)
        NULLIFY(decompositionLines)
        CALL DecompositionTopology_LinesGet(decompositionTopology,decompositionLines,err,error,*999)
        !identify all non-ghost boundary lines
        ALLOCATE(candidateElements(numberOfLines),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate candidate elements.",err,error,*999)
        ALLOCATE(candidateLinesFaces(numberOfLines),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate candidate lines faces.",err,error,*999)
        DO elementIdx=1,domainMappingElements%NUMBER_OF_LOCAL
          IF(decompositionElements%elements(elementIdx)%BOUNDARY_ELEMENT) THEN
            DO lineFaceIdx=1,SIZE(decompositionElements%elements(elementIdx)%ELEMENT_LINES,1)
              lineFaceNumber=decompositionElements%elements(elementIdx)%ELEMENT_LINES(lineFaceIdx)
              IF(decompositionLines%lines(lineFaceNumber)%BOUNDARY_LINE) THEN
                numberOfCandidates=numberOfCandidates+1
                candidateLinesFaces(numberOfCandidates)=lineFaceIdx
                candidateElements(numberOfCandidates)=elementIdx
              END IF
            ENDDO !lineFaceIdx
          ENDIF
        ENDDO !elementIdx
      CASE(DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE)
        NULLIFY(decompositionFaces)
        CALL DecompositionTopology_FacesGet(decompositionTopology,decompositionFaces,err,error,*999)                    
        !Identify all non-ghost boundary faces
        IF(decomposition%numberOfDimensions>=2) THEN
          ALLOCATE(candidateElements(numberOfFaces),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate candidate elements.",err,error,*999)
          ALLOCATE(candidateLinesFaces(numberOfFaces),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate candidate lines faces.",err,error,*999)
          DO elementIdx=1,domainMappingElements%NUMBER_OF_LOCAL
            IF(decompositionElements%elements(elementIdx)%BOUNDARY_ELEMENT) THEN
              DO lineFaceIdx=1,SIZE(decompositionElements%elements(elementIdx)%ELEMENT_FACES,1)
                lineFaceNumber=decompositionElements%elements(elementIdx)%ELEMENT_FACES(lineFaceIdx)
                IF(decompositionFaces%faces(lineFaceNumber)%BOUNDARY_FACE) THEN
                  numberOfCandidates=numberOfCandidates+1
                  candidateLinesFaces(numberOfCandidates)=lineFaceIdx
                  candidateElements(numberOfCandidates)=elementIdx
                ENDIF
              ENDDO !lineFaceIdx
            ENDIF
          ENDDO !elementIdx
        ELSE
          CALL FlagError("Decomposition number of dimensions has to be >= 2 for a faces projection type.",err,error,*999)        
        ENDIF
      CASE(DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE)
        !Identify all non-ghost elements
        IF(dataProjection%numberOfXi==decomposition%numberOfDimensions) THEN
          ALLOCATE(candidateElements(numberOfElements),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate candidate elements.",err,error,*999)
          DO elementIdx=1,domainMappingElements%NUMBER_OF_LOCAL
            numberOfCandidates=numberOfCandidates+1
            candidateElements(numberOfCandidates)=elementIdx
          ENDDO !elementIdx
        ELSE
          CALL FlagError("Data projection number of xi has to equal to decomposition mesh number of dimensions",err,error,*999)
        ENDIF
      CASE DEFAULT
        localError="The data projection type of "//TRIM(NumberToVString(dataProjection%projectionType,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDIF
    !##############################################################################################################
    !find the clostest elements/faces/lines for each point in the current computational node base on starting xi
    !the clostest elements/faces/lines are required to shrink down on the list of possible projection candiates
    numberOfClosestCandidates=MIN(dataProjection%numberOfClosestElements,dataProjection%maxNumberOfCandidates)
    !Allocated and store he information for each data point. The information has to be stored in the corresponding
    !rows for them to be contiguous in memory for easy MPI access
    ALLOCATE(closestElements(numberOfDataPoints,numberOfClosestCandidates),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate closest elements.",err,error,*999)
    closestElements=0
    IF(boundaryProjection) THEN
      ALLOCATE(closestLinesFaces(numberOfDataPoints,numberOfClosestCandidates),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate closest lines faces.",err,error,*999)
      closestLinesFaces=0
    ENDIF
    !Allocate and store the information for each data point. The information has to be stored in the corresponding
    !rows for them to be contiguous in memory for easy MPI access
    ALLOCATE(closestDistances(numberOfDataPoints,numberOfClosestCandidates),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate closest distances.",err,error,*999)
    closestDistances=0.0_DP
    SELECT CASE(dataProjection%projectionType)
    CASE(DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE)
      !Find closest candidate lines
      DO dataPointIdx=1,numberOfDataPoints
        IF(ALLOCATED(dataProjection%dataProjectionCandidates(dataPointIdx)%candidateElementNumbers).AND. &
          & ALLOCATED(dataProjection%dataProjectionCandidates(dataPointIdx)%localFaceLineNumbers)) THEN
          dataNumberOfCandidates=SIZE(dataProjection%dataProjectionCandidates(dataPointIdx)%candidateElementNumbers,1)
          CALL DataProjection_ClosestLinesFind(dataProjection,interpolatedPoint, &
            & dataPoints%dataPoints(dataPointIdx)%position,dataNumberOfCandidates,dataProjection% &
            & dataProjectionCandidates(dataPointIdx)%candidateElementNumbers, &
            & dataProjection%dataProjectionCandidates(dataPointIdx)%localFaceLineNumbers, &
            & closestElements(dataPointIdx,:),closestLinesFaces(dataPointIdx,:), &
            & closestDistances(dataPointIdx,:),err,error,*999)
        ELSE
          CALL DataProjection_ClosestLinesFind(dataProjection,interpolatedPoint, &
            & dataPoints%dataPoints(dataPointIdx)%position,numberOfCandidates,candidateElements, &
            & candidateLinesFaces,closestElements(dataPointIdx,:),closestLinesFaces(dataPointIdx,:), &
            & closestDistances(dataPointIdx,:),err,error,*999)
        ENDIF
      ENDDO !dataPointIdx
    CASE(DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE)
      !Find closest candidate faces      
      DO dataPointIdx=1,numberOfDataPoints
        IF(ALLOCATED(dataProjection%dataProjectionCandidates(dataPointIdx)%candidateElementNumbers).AND. &
          & ALLOCATED(dataProjection%dataProjectionCandidates(dataPointIdx)%localFaceLineNumbers)) THEN
          dataNumberOfCandidates=SIZE(dataProjection%dataProjectionCandidates(dataPointIdx)%candidateElementNumbers,1)
          CALL DataProjection_ClosestFacesFind(dataProjection,interpolatedPoint, &
            & dataPoints%dataPoints(dataPointIdx)%position,dataNumberOfCandidates,dataProjection% &
            & dataProjectionCandidates(dataPointIdx)%candidateElementNumbers, &
            & dataProjection%dataProjectionCandidates(dataPointIdx)%localFaceLineNumbers, &
            & closestElements(dataPointIdx,:),closestLinesFaces(dataPointIdx,:), &
            & closestDistances(dataPointIdx,:),err,error,*999)
        ELSE
          CALL DataProjection_ClosestFacesFind(dataProjection,interpolatedPoint, &
            & dataPoints%dataPoints(dataPointIdx)%position,numberOfCandidates,candidateElements, &
            & candidateLinesFaces,closestElements(dataPointIdx,:),closestLinesFaces(dataPointIdx,:), &
            & closestDistances(dataPointIdx,:),err,error,*999)
        ENDIF
      ENDDO !dataPointIdx
    CASE(DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE)
      !Find closest candidate elements
      DO dataPointIdx=1,numberOfDataPoints
        IF(ALLOCATED(dataProjection%dataProjectionCandidates(dataPointIdx)%candidateElementNumbers)) THEN
          dataNumberOfCandidates=SIZE(dataProjection%dataProjectionCandidates(dataPointIdx)%candidateElementNumbers,1)
          CALL DataProjection_ClosestElementsFind(dataProjection,interpolatedPoint, &
            & dataPoints%dataPoints(dataPointIdx)%position,dataNumberOfCandidates, &
            & dataProjection%dataProjectionCandidates(dataPointIdx)%candidateElementNumbers, &
            & closestElements(dataPointIdx,:),closestDistances(dataPointIdx,:),err,error,*999)
        ELSE
          CALL DataProjection_ClosestElementsFind(dataProjection,interpolatedPoint, &
            & dataPoints%dataPoints(dataPointIdx)%position,numberOfCandidates,candidateElements, &
            & closestElements(dataPointIdx,:),closestDistances(dataPointIdx,:),err,error,*999)
        ENDIF
      ENDDO !dataPointIdx
    CASE DEFAULT
      localError="The data projection type of "//TRIM(NumberToVString(dataProjection%projectionType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    !###################################################################################################################
    !Newton project data point to the list of closest elements, faces or lines
    !project the data points to each of the closest elements,
    !use MPI if number of computational nodes is greater than 1
    IF(numberOfComputationalNodes>1) THEN
      !Use MPI
      !Allocate arrays for MPI communication
      ALLOCATE(globalToLocalNumberOfClosestCandidates(numberOfDataPoints),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate global to local number of closest elements.",err,error,*999)
      ALLOCATE(globalNumberOfClosestCandidates(numberOfComputationalNodes),STAT=err) 
      IF(err/=0) CALL FlagError("Could not allocate global number of closest candidates.",err,error,*999)
      ALLOCATE(globalMPIDisplacements(numberOfComputationalNodes),STAT=err) 
      IF(err/=0) CALL FlagError("Could not allocate global MPI displacements.",err,error,*999)
      ALLOCATE(globalNumberOfProjectedPoints(numberOfComputationalNodes),STAT=err) 
      IF(err/=0) CALL FlagError("Could not allocate all number of projected points.",err,error,*999)
      ALLOCATE(projectionExitTag(numberOfDataPoints),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate projected.",err,error,*999)
      ALLOCATE(projectedElement(numberOfDataPoints),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate projected element.",err,error,*999)
      IF(boundaryProjection) THEN
        ALLOCATE(projectedLineFace(numberOfDataPoints),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate projected sub element.",err,error,*999)
      ENDIF
      !projectedDistance(2,:) stores the compuational node number, the information for each data point has to be stored
      !in the corresponding column for MPI_ALLREDUCE with location return to work
      ALLOCATE(projectedDistance(2,numberOfDataPoints),STAT=err) 
      IF(err/=0) CALL FlagError("Could not allocate projected distance.",err,error,*999)
      !the information for each data point is stored in the corresponding column to be consistent with projectedDistance
      ALLOCATE(projectedXi(dataProjection%numberOfXi,numberOfDataPoints),STAT=err)
      !ZJW: The projection vector for each data point is stored in the corresponding colum to be consistent with
      !projectedDistance. 
      ALLOCATE(projectionVectors(dataProjection%numberOfCoordinates,numberOfDataPoints), STAT=err) 
      IF(err/=0) CALL FlagError("Could not allocate projected vectors.",err,error,*999)
      ALLOCATE(sortingIndices2(numberOfDataPoints),STAT=err) 
      IF(err/=0) CALL FlagError("Could not allocate sorting indices 2.",err,error,*999)          
      !gather and distribute the number of closest elements from all computational nodes
      CALL MPI_ALLGATHER(numberOfClosestCandidates,1,MPI_INTEGER,globalNumberOfClosestCandidates,1,MPI_INTEGER, &
        & computationalEnvironment%mpiCommunicator,MPIIError)
      CALL MPI_ERROR_CHECK("MPI_ALLGATHER",MPIIError,err,error,*999)
      !Sum all number of closest candidates from all computational nodes
      totalNumberOfClosestCandidates=SUM(globalNumberOfClosestCandidates,1) 
      !Allocate arrays to store information gathered from all computational node
      !The information for each data point is stored in the corresponding row so they are contiguous in memory for
      !easy MPI access
      ALLOCATE(globalClosestDistances(numberOfDataPoints,totalNumberOfClosestCandidates),STAT=err) 
      IF(err/=0) CALL FlagError("Could not allocate all closest distances.",err,error,*999)
      ALLOCATE(sortingIndices1(totalNumberOfClosestCandidates),STAT=err) 
      IF(err/=0) CALL FlagError("Could not allocate sorting indices 1.",err,error,*999)
      !MPI:create and commit MPI_TYPE_CONTIGUOUS      
      CALL MPI_TYPE_CONTIGUOUS(numberOfDataPoints,MPI_DOUBLE_PRECISION,MPIClosestDistances,MPIIError)
      CALL MPI_ERROR_CHECK("MPI_TYPE_CONTIGUOUS",MPIIError,err,error,*999)        
      CALL MPI_TYPE_COMMIT(MPIClosestDistances,MPIIError)
      CALL MPI_ERROR_CHECK("MPI_TYPE_COMMIT",MPIIError,err,error,*999)
      !Create displacement vectors for MPI_ALLGATHERV
      globalMPIDisplacements(1)=0
      DO computationalNodeIdx=1,(numberOfComputationalNodes-1)
        globalMPIDisplacements(computationalNodeIdx+1)=globalMPIDisplacements(computationalNodeIdx)+ &
          & globalNumberOfClosestCandidates(computationalNodeIdx)
      ENDDO !computationalNodeIdx
      !Share closest element distances between all domains
      CALL MPI_ALLGATHERV(closestDistances(1,1),numberOfClosestCandidates,MPIClosestDistances, &
        & globalClosestDistances,globalNumberOfClosestCandidates,globalMPIDisplacements, &
        & MPIClosestDistances,computationalEnvironment%mpiCommunicator,MPIIError)
      CALL MPI_ERROR_CHECK("MPI_ALLGATHERV",MPIIError,err,error,*999)
      reducedNumberOfCLosestCandidates=MIN(dataProjection%numberOfClosestElements,totalNumberOfClosestCandidates)
      projectedDistance(2,:)=myComputationalNode
      !Find the globally closest distance in the current domain
      DO dataPointIdx=1,numberOfDataPoints
        CALL Sorting_BubbleIndexSort(globalClosestDistances(dataPointIdx,:),sortingIndices1,err,error,*999)
        sortingIndices1(1:totalNumberOfClosestCandidates)=sortingIndices1(1:totalNumberOfClosestCandidates)- &
          & globalMPIDisplacements(myComputationalNode+1) !shift the index to current computational node
        globalToLocalNumberOfClosestCandidates(dataPointIdx)=0
        DO elementIdx=1,reducedNumberOfCLosestCandidates
          !Sorted index indicates it is in the current computational domain
          IF((sortingIndices1(elementIdx)>=1).AND. &
            & (sortingIndices1(elementIdx)<=globalNumberOfClosestCandidates(myComputationalNode+1))) &
            & globalToLocalNumberOfClosestCandidates(dataPointIdx)= &
            & globalToLocalNumberOfClosestCandidates(dataPointIdx)+1
        ENDDO !elementIdx
        !Assign initial distance to something large                           
        projectedDistance(1,dataPointIdx)=globalClosestDistances(dataPointIdx,totalNumberOfClosestCandidates) 
      ENDDO !dataPointIdx
      SELECT CASE(dataProjection%projectionType)
      CASE(DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE)
        !Newton project to closest lines, and find miminum projection
        DO dataPointIdx=1,numberOfDataPoints
          numberOfClosestCandidates=globalToLocalNumberOfClosestCandidates(dataPointIdx)
          IF(numberOfClosestCandidates>0) THEN 
            CALL DataProjection_NewtonLinesEvaluate(dataProjection,interpolatedPoint, &
              & dataPoints%dataPoints(dataPointIdx)%position,closestElements( &
              & dataPointIdx,1:numberOfClosestCandidates),closestLinesFaces(dataPointIdx,1: &
              & numberOfClosestCandidates),projectionExitTag(dataPointIdx),projectedElement(dataPointIdx),  &
              & projectedLineFace(dataPointIdx),projectedDistance(1,dataPointIdx),projectedXi(:,dataPointIdx), &
              & projectionVectors(:,dataPointIdx),err,error,*999)
            !Map the element number to global number
            projectedElement(dataPointIdx)=domainMappingElements%LOCAL_TO_GLOBAL_MAP(projectedElement(dataPointIdx))
          ENDIF
        ENDDO !dataPointIdx
      CASE(DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE)
        !Newton project to closest faces, and find miminum projection
        DO dataPointIdx=1,numberOfDataPoints
          numberOfClosestCandidates=globalToLocalNumberOfClosestCandidates(dataPointIdx)
          IF(numberOfClosestCandidates>0) THEN 
            CALL DataProjection_NewtonFacesEvaluate(dataProjection,interpolatedPoint, &
              & dataPoints%dataPoints(dataPointIdx)%position,closestElements( &
              & dataPointIdx,1:numberOfClosestCandidates),closestLinesFaces(dataPointIdx, &
              & 1:numberOfClosestCandidates),projectionExitTag(dataPointIdx),projectedElement(dataPointIdx), &
              & projectedLineFace(dataPointIdx),projectedDistance(1,dataPointIdx),projectedXi(:,dataPointIdx), &
              & projectionVectors(:,dataPointIdx),err,error,*999)
            !Map the element number to global number
            projectedElement(dataPointIdx)=domainMappingElements%LOCAL_TO_GLOBAL_MAP(projectedElement(dataPointIdx))
          ENDIF
        ENDDO !dataPointIdx
      CASE(DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE)
        !Newton project to closest elements, and find miminum projection
        SELECT CASE(dataProjection%numberOfXi)
        CASE(1) !1D element
          DO dataPointIdx=1,numberOfDataPoints
            numberOfClosestCandidates=globalToLocalNumberOfClosestCandidates(dataPointIdx)
            IF(numberOfClosestCandidates>0) THEN 
              CALL DataProjection_NewtonElementsEvaluate_1(dataProjection,interpolatedPoint,dataPoints% &
                & dataPoints(dataPointIdx)%position,closestElements(dataPointIdx, &
                & 1:numberOfClosestCandidates),projectionExitTag(dataPointIdx), &
                & projectedElement(dataPointIdx),projectedDistance(1,dataPointIdx), &
                & projectedXi(:,dataPointIdx),projectionVectors(:,dataPointIdx),err,error,*999)
              !Map the element number to global number
              projectedElement(dataPointIdx)=domainMappingElements%LOCAL_TO_GLOBAL_MAP(projectedElement(dataPointIdx))
            ENDIF
          ENDDO !dataPointIdx
        CASE(2) !2D element
          DO dataPointIdx=1,numberOfDataPoints
            numberOfClosestCandidates=globalToLocalNumberOfClosestCandidates(dataPointIdx)
            IF(numberOfClosestCandidates>0) THEN 
              CALL DataProjection_NewtonElementsEvaluate_2(dataProjection,interpolatedPoint,dataPoints% &
                & dataPoints(dataPointIdx)%position,closestElements(dataPointIdx, &
                & 1:numberOfClosestCandidates),projectionExitTag(dataPointIdx), &
                & projectedElement(dataPointIdx),projectedDistance(1,dataPointIdx), &
                & projectedXi(:,dataPointIdx),projectionVectors(:,dataPointIdx), &
                & err,error,*999)
              !Map the element number to global number                        
              projectedElement(dataPointIdx)=domainMappingElements%LOCAL_TO_GLOBAL_MAP(projectedElement(dataPointIdx))
            ENDIF
          ENDDO !dataPointIdx
        CASE(3) !3D element
          DO dataPointIdx=1,numberOfDataPoints
            numberOfClosestCandidates=globalToLocalNumberOfClosestCandidates(dataPointIdx)
            IF(numberOfClosestCandidates>0) THEN 
              CALL DataProjection_NewtonElementsEvaluate_3(dataProjection,interpolatedPoint,dataPoints% &
                & dataPoints(dataPointIdx)%position,closestElements(dataPointIdx, &
                & 1:numberOfClosestCandidates),projectionExitTag(dataPointIdx), &
                & projectedElement(dataPointIdx),projectedDistance(1,dataPointIdx), &
                & projectedXi(:,dataPointIdx),projectionVectors(:,dataPointIdx),err,error,*999)
              !Map the element number to global number
              projectedElement(dataPointIdx)=domainMappingElements%LOCAL_TO_GLOBAL_MAP(projectedElement(dataPointIdx))
            ENDIF
          ENDDO !dataPointIdx
        CASE DEFAULT
          localError="The data projection number of xi directions of "// &
            & TRIM(NumberToVString(dataProjection%numberOfXi,"*",err,error))// &
            & " is invalid. The number of directions should be >= 1 and <= 3."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The data projection type of "//TRIM(NumberToVString(dataProjection%projectionType,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Find the shortest projected distance in all domains
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,projectedDistance,numberOfDataPoints,MPI_2DOUBLE_PRECISION,MPI_MINLOC, &
        & computationalEnvironment%mpiCommunicator,MPIIError)
      CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",MPIIError,err,error,*999)
      !Sort the computational node/rank from 0 to number of computational node
      CALL Sorting_BubbleIndexSort(projectedDistance(2,:),sortingIndices2,err,error,*999)
      DO computationalNodeIdx=0,(numberOfComputationalNodes-1)
        globalNumberOfProjectedPoints(computationalNodeIdx+1)=COUNT(ABS(projectedDistance(2,:)- &
          & REAL(computationalNodeIdx))<ZERO_TOLERANCE)
      ENDDO !computationalNodeIdx
      startIdx=SUM(globalNumberOfProjectedPoints(1:myComputationalNode))+1
      finishIdx=startIdx+globalNumberOfProjectedPoints(myComputationalNode+1)-1
      !create displacement vectors for MPI_ALLGATHERV          
      DO computationalNodeIdx=1,(numberOfComputationalNodes-1)
        globalMPIDisplacements(computationalNodeIdx+1)=globalMPIDisplacements(computationalNodeIdx)+ &
          & globalNumberOfProjectedPoints(computationalNodeIdx)
      ENDDO !computationalNodeIdx  
      !Shares minimum projection information between all domains
      CALL MPI_ALLGATHERV(projectedElement(sortingIndices2(startIdx:finishIdx)),globalNumberOfProjectedPoints( &
        & myComputationalNode+1),MPI_INTEGER,projectedElement,globalNumberOfProjectedPoints, &
        & globalMPIDisplacements,MPI_INTEGER,computationalEnvironment%mpiCommunicator,MPIIError) !projectedElement
      CALL MPI_ERROR_CHECK("MPI_ALLGATHERV",MPIIError,err,error,*999)
      IF(boundaryProjection) THEN
        CALL MPI_ALLGATHERV(projectedLineFace(sortingIndices2(startIdx:finishIdx)),globalNumberOfProjectedPoints( &
          & myComputationalNode+1),MPI_INTEGER,projectedLineFace,globalNumberOfProjectedPoints, &
          & globalMPIDisplacements,MPI_INTEGER,computationalEnvironment%mpiCommunicator,MPIIError) !projectedLineFace
        CALL MPI_ERROR_CHECK("MPI_ALLGATHERV",MPIIError,err,error,*999) 
      ENDIF
      DO xiIdx=1,dataProjection%numberOfXi
        CALL MPI_ALLGATHERV(projectedXi(xiIdx,sortingIndices2(startIdx:finishIdx)),globalNumberOfProjectedPoints( &
          & myComputationalNode+1),MPI_DOUBLE_PRECISION,projectedXi(xiIdx,:),globalNumberOfProjectedPoints, &
          & globalMPIDisplacements,MPI_DOUBLE_PRECISION,computationalEnvironment%mpiCommunicator,MPIIError) !projectedXi
        CALL MPI_ERROR_CHECK("MPI_ALLGATHERV",MPIIError,err,error,*999)
      ENDDO !xiIdx
      CALL MPI_ALLGATHERV(projectionExitTag(sortingIndices2(startIdx:finishIdx)),globalNumberOfProjectedPoints( &
        & myComputationalNode+1),MPI_INTEGER,projectionExitTag,globalNumberOfProjectedPoints, &
        & globalMPIDisplacements,MPI_INTEGER,computationalEnvironment%mpiCommunicator,MPIIError) !projectionExitTag
      CALL MPI_ERROR_CHECK("MPI_ALLGATHERV",MPIIError,err,error,*999)
      DO xiIdx=1,dataProjection%numberOfCoordinates
        CALL MPI_ALLGATHERV(projectionVectors(xiIdx, sortingIndices2(startIdx:finishIdx)), &
          & globalNumberOfProjectedPoints(myComputationalNode+1),MPI_DOUBLE_PRECISION,projectionVectors(xiIdx,:), &
          & globalNumberOfProjectedPoints,globalMPIDisplacements,MPI_DOUBLE_PRECISION,computationalEnvironment% &
          & mpiCommunicator,MPIIError)  !projectionVectors
        CALL MPI_ERROR_CHECK("MPI_ALLGATHERV",MPIIError,err,error,*999)
      ENDDO
      !Assign projection information to projected points
      DO dataPointIdx=1,numberOfDataPoints
        dataProjection%dataProjectionResults(sortingIndices2(dataPointIdx))%exitTag=projectionExitTag(dataPointIdx)
        dataProjection%dataProjectionResults(sortingIndices2(dataPointIdx))%elementNumber=projectedElement(dataPointIdx)
        dataProjection%dataProjectionResults(dataPointIdx)%DISTANCE=projectedDistance(1,dataPointIdx)
        dataProjection%dataProjectionResults(sortingIndices2(dataPointIdx))%xi(1:dataProjection%numberOfXi)= &
          & projectedXi(1:dataProjection%numberOfXi,dataPointIdx)
        dataProjection%dataProjectionResults(sortingIndices2(dataPointIdx))%projectionVector( &
          & 1:dataProjection%numberOfCoordinates)=projectionVectors(1:dataProjection%numberOfCoordinates,dataPointIdx)
      ENDDO !dataPointIdx
      projectedXi(:,sortingIndices2)=projectedXi
      projectionVectors(:, sortingIndices2)=projectionVectors
      projectedElement(sortingIndices2)=projectedElement       
      IF(dataProjection%projectionType==DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE) THEN
        DO dataPointIdx=1,numberOfDataPoints          
          dataProjection%dataProjectionResults(sortingIndices2(dataPointIdx))%elementLineFaceNumber=projectedLineFace(dataPointIdx)
        ENDDO !dataPointIdx
      ELSE IF(dataProjection%projectionType==DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE) THEN
        DO dataPointIdx=1,numberOfDataPoints          
          dataProjection%dataProjectionResults(sortingIndices2(dataPointIdx))%elementLineFaceNumber=projectedLineFace(dataPointIdx)
        ENDDO !dataPointIdx
      ENDIF
    ELSE
      !No need to use mpi
      SELECT CASE(dataProjection%projectionType)
      CASE(DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE)
        !Newton project to closest lines, and find miminum projection
        DO dataPointIdx=1,numberOfDataPoints
          CALL DataProjection_NewtonLinesEvaluate(dataProjection,interpolatedPoint,dataPoints%dataPoints( &
            & dataPointIdx)%position,closestElements(dataPointIdx,:),closestLinesFaces(dataPointIdx,:), &
            & dataProjection%dataProjectionResults(dataPointIdx)%exitTag,dataProjection% &
            & dataProjectionResults(dataPointIdx)%elementNumber,dataProjection% &
            & dataProjectionResults(dataPointIdx)%elementLineFaceNumber,dataProjection%dataProjectionResults( &
            & dataPointIdx)%DISTANCE,dataProjection%dataProjectionResults(dataPointIdx)%xi, &
            & dataProjection%dataProjectionResults(dataPointIdx)%projectionVector,err,error,*999)
        ENDDO !dataPointIdx
      CASE(DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE) 
        !Newton project to closest faces, and find miminum projection
        DO dataPointIdx=1,numberOfDataPoints
          CALL DataProjection_NewtonFacesEvaluate(dataProjection,interpolatedPoint,dataPoints%dataPoints( &
            & dataPointIdx)%position,closestElements(dataPointIdx,:),closestLinesFaces(dataPointIdx,:), &
            & dataProjection%dataProjectionResults(dataPointIdx)%exitTag,dataProjection% &
            & dataProjectionResults(dataPointIdx)%elementNumber,dataProjection% &
            & dataProjectionResults(dataPointIdx)%elementLineFaceNumber,dataProjection%dataProjectionResults( &
            & dataPointIdx)%DISTANCE,dataProjection%dataProjectionResults(dataPointIdx)%xi, &
            & dataProjection%dataProjectionResults(dataPointIdx)%projectionVector,err,error,*999)
        ENDDO !dataPointIdx
      CASE(DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE)        
        !Newton project to closest elements, and find miminum projection
        SELECT CASE(dataProjection%numberOfXi)
        CASE(1) !1D mesh
          DO dataPointIdx=1,numberOfDataPoints
            CALL DataProjection_NewtonElementsEvaluate_1(dataProjection,interpolatedPoint,dataPoints% &
              & dataPoints(dataPointIdx)%position,closestElements(dataPointIdx,:),dataProjection% &
              & dataProjectionResults(dataPointIdx)%exitTag,dataProjection%dataProjectionResults( &
              & dataPointIdx)%elementNumber,dataProjection%dataProjectionResults(dataPointIdx)%distance, &
              & dataProjection%dataProjectionResults(dataPointIdx)%xi, &
              & dataProjection%dataProjectionResults(dataPointIdx)%projectionVector,&
              & err,error,*999)
          ENDDO !dataPointIdx
        CASE(2) !2D mesh
          DO dataPointIdx=1,numberOfDataPoints
            CALL DataProjection_NewtonElementsEvaluate_2(dataProjection,interpolatedPoint,dataPoints% &
              & dataPoints(dataPointIdx)%position,closestElements(dataPointIdx,:),dataProjection% &
              & dataProjectionResults(dataPointIdx)%exitTag,dataProjection%dataProjectionResults( &
              & dataPointIdx)%elementNumber,dataProjection%dataProjectionResults(dataPointIdx)%distance, &
              & dataProjection%dataProjectionResults(dataPointIdx)%xi, &
              & dataProjection%dataProjectionResults(dataPointIdx)%projectionVector, &
              & err,error,*999)
          ENDDO !dataPointIdx
        CASE(3) !3D mesh
          DO dataPointIdx=1,numberOfDataPoints
            CALL DataProjection_NewtonElementsEvaluate_3(dataProjection,interpolatedPoint,dataPoints% &
              & dataPoints(dataPointIdx)%position,closestElements(dataPointIdx,:),dataProjection% &
              & dataProjectionResults(dataPointIdx)%exitTag,dataProjection%dataProjectionResults( &
              & dataPointIdx)%elementNumber,dataProjection%dataProjectionResults(dataPointIdx)%distance, &
              & dataProjection%dataProjectionResults(dataPointIdx)%xi, &
              & dataProjection%dataProjectionResults(dataPointIdx)%projectionVector, &
              & err,error,*999)
          ENDDO !dataPointIdx
        CASE DEFAULT
          localError="The data projection number of xi of "//TRIM(NumberToVString(dataProjection%numberOfXi,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The data projection type of "//TRIM(NumberToVString(dataProjection%projectionType,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDIF !numberOfComputationalNodes>1
    !Compute full elemental xi
    IF(dataProjection%numberOfXi==dataProjection%numberOfElementXi) THEN
      DO dataPointIdx=1,numberOfDataPoints
        dataProjection%dataProjectionResults(dataPointIdx)%elementXi= &
          & dataProjection%dataProjectionResults(dataPointIdx)%xi
      ENDDO !dataPointIdx
    ELSE
      DO dataPointIdx=1,numberOfDataPoints
        elementNumber=dataProjection%dataProjectionResults(dataPointIdx)%elementNumber
        localLineFaceNumber=dataProjection%dataProjectionResults(dataPointIdx)%elementLineFaceNumber
        basis=>domainElements%elements(elementNumber)%basis
        CALL Basis_BoundaryXiToXi(basis,localLineFaceNumber,dataProjection% &
          & dataProjectionResults(dataPointIdx)%xi(1:dataProjection%numberOfXi),dataProjection% &
          & dataProjectionResults(dataPointIdx)%elementXi,err,error,*999)
      ENDDO !dataPointIdx
    ENDIF
    !Evaluate errors
    CALL DataProjection_ErrorsCalculate(dataProjection,err,error,*999)
    !Set projected flag
    dataProjection%dataProjectionProjected=.TRUE.

    !Deallocate arrays used within this routine
    IF(ALLOCATED(candidateElements)) DEALLOCATE(candidateElements)
    IF(ALLOCATED(candidateLinesFaces)) DEALLOCATE(candidateLinesFaces)
    IF(ALLOCATED(closestElements)) DEALLOCATE(closestElements)
    IF(ALLOCATED(closestLinesFaces)) DEALLOCATE(closestLinesFaces)
    IF(ALLOCATED(closestDistances)) DEALLOCATE(closestDistances)
    IF(ALLOCATED(globalToLocalNumberOfClosestCandidates)) DEALLOCATE(globalToLocalNumberOfClosestCandidates)
    IF(ALLOCATED(globalNumberOfClosestCandidates)) DEALLOCATE(globalNumberOfClosestCandidates)
    IF(ALLOCATED(globalMPIDisplacements)) DEALLOCATE(globalMPIDisplacements)
    IF(ALLOCATED(globalNumberOfProjectedPoints)) DEALLOCATE(globalNumberOfProjectedPoints)
    IF(ALLOCATED(projectionExitTag)) DEALLOCATE(projectionExitTag)
    IF(ALLOCATED(projectedElement)) DEALLOCATE(projectedElement)
    IF(ALLOCATED(projectedLineFace)) DEALLOCATE(projectedLineFace)
    IF(ALLOCATED(projectedDistance)) DEALLOCATE(projectedDistance)
    IF(ALLOCATED(projectedXi)) DEALLOCATE(projectedXi)
    IF(ALLOCATED(projectionVectors)) DEALLOCATE(projectionVectors)
    IF(ALLOCATED(globalClosestDistances)) DEALLOCATE(globalClosestDistances)
    IF(ALLOCATED(sortingIndices1)) DEALLOCATE(sortingIndices1)
    IF(ALLOCATED(sortingIndices2)) DEALLOCATE(sortingIndices2)
    
    EXITS("DataProjection_DataPointsProjectionEvaluate")
    RETURN
999 IF(ALLOCATED(candidateElements)) DEALLOCATE(candidateElements)
    IF(ALLOCATED(candidateLinesFaces)) DEALLOCATE(candidateLinesFaces)
    IF(ALLOCATED(closestElements)) DEALLOCATE(closestElements)
    IF(ALLOCATED(closestLinesFaces)) DEALLOCATE(closestLinesFaces)
    IF(ALLOCATED(closestDistances)) DEALLOCATE(closestDistances)
    IF(ALLOCATED(globalToLocalNumberOfClosestCandidates)) DEALLOCATE(globalToLocalNumberOfClosestCandidates)
    IF(ALLOCATED(globalNumberOfClosestCandidates)) DEALLOCATE(globalNumberOfClosestCandidates)
    IF(ALLOCATED(globalMPIDisplacements)) DEALLOCATE(globalMPIDisplacements)
    IF(ALLOCATED(globalNumberOfProjectedPoints)) DEALLOCATE(globalNumberOfProjectedPoints)
    IF(ALLOCATED(projectionExitTag)) DEALLOCATE(projectionExitTag)
    IF(ALLOCATED(projectedElement)) DEALLOCATE(projectedElement)
    IF(ALLOCATED(projectedLineFace)) DEALLOCATE(projectedLineFace)
    IF(ALLOCATED(projectedDistance)) DEALLOCATE(projectedDistance)
    IF(ALLOCATED(projectedXi)) DEALLOCATE(projectedXi)
    IF(ALLOCATED(projectionVectors)) DEALLOCATE(projectionVectors)
    IF(ALLOCATED(globalClosestDistances)) DEALLOCATE(globalClosestDistances)
    IF(ALLOCATED(sortingIndices1)) DEALLOCATE(sortingIndices1)
    IF(ALLOCATED(sortingIndices2)) DEALLOCATE(sortingIndices2)
998 ERRORSEXITS("DataProjection_DataPointsProjectionEvaluate",err,error)
    RETURN 1
    
  END SUBROUTINE DataProjection_DataPointsProjectionEvaluate
  
  !
  !================================================================================================================================
  !
  
  !>Evaluate the data points position in a field based on data projection
  SUBROUTINE DataProjection_DataPointsPositionEvaluate(dataProjection,field,fieldVariableType,fieldParameterSetType,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<Data projection to give the xi locations and element number for the data points
    TYPE(FIELD_TYPE), POINTER :: field !<A pointer to the field to be interpolated
    INTEGER(INTG), INTENT(IN) :: fieldVariableType !<The field variable type to be interpolated
    INTEGER(INTG), INTENT(IN) :: fieldParameterSetType !<The parameter set to be interpolated
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointIdx,elementNumber,coordinateIdx
    TYPE(DataPointsType), POINTER :: dataPoints
    TYPE(FIELD_INTERPOLATED_POINT_PTR_TYPE), POINTER :: interpolatedPoints(:)
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: interpolatedPoint
    TYPE(FIELD_INTERPOLATION_PARAMETERS_PTR_TYPE), POINTER :: interpolationParameters(:)
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: fieldParameterSet
    
    ENTERS("DataProjection_DataPointsPositionEvaluate",err,error,*999)
    
    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionFinished) CALL FlagError("Data projection has not been finished.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionProjected) CALL FlagError("Data projection has not been evaluated.",err,error,*999)
    IF(.NOT.ALLOCATED(dataProjection%dataProjectionResults)) &
       CALL FlagError("Data projection projection results is not allocated.",err,error,*999)
    NULLIFY(dataPoints)
    CALL DataProjection_DataPointsGet(dataProjection,dataPoints,err,error,*999)
    IF(.NOT.ASSOCIATED(field)) CALL FlagError("Field is not associated.",err,error,*999)
    IF(field%TYPE/=FIELD_GEOMETRIC_TYPE.AND.field%TYPE/=FIELD_GEOMETRIC_GENERAL_TYPE) &
      & CALL FlagError("Cannot evaluate data points position on field other than geometric or geometric general type.", &
      & err,error,*999)
    NULLIFY(fieldParameterSet)
    CALL Field_ParameterSetGet(field,fieldVariableType,fieldParameterSetType,fieldParameterSet,err,error,*999)
    
    NULLIFY(interpolatedPoints)
    NULLIFY(interpolationParameters)
    CALL Field_InterpolationParametersInitialise(field,interpolationParameters,err,error,*999, &
      & FIELD_GEOMETRIC_COMPONENTS_TYPE)
    CALL Field_InterpolatedPointsInitialise(interpolationParameters,interpolatedPoints,err,error,*999, &
      & FIELD_GEOMETRIC_COMPONENTS_TYPE)
    interpolatedPoint=>interpolatedPoints(fieldVariableType)%ptr
    !Loop through data points 
    DO dataPointIdx=1,dataPoints%numberOfDataPoints
      elementNumber=dataProjection%dataProjectionResults(dataPointIdx)%elementNumber
      CALL Field_InterpolationParametersElementGet(fieldParameterSetType,elementNumber, &
        & interpolationParameters(fieldVariableType)%ptr,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
      CALL Field_InterpolateXi(NO_PART_DERIV,dataProjection%dataProjectionResults(dataPointIdx)%xi, &
        & interpolatedPoint,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
      DO coordinateIdx=1,SIZE(dataPoints%dataPoints(dataPointIdx)%position,1)
        dataPoints%dataPoints(dataPointIdx)%position(coordinateIdx)= &
          & interpolatedPoint%values(coordinateIdx,NO_PART_DERIV)
      ENDDO !coordinateIdx     
    ENDDO !dataPointIdx
    
    EXITS("DataProjection_DataPointsPositionEvaluate")
    RETURN
999 ERRORSEXITS("DataProjection_DataPointsPositionEvaluate",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_DataPointsPositionEvaluate
  
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

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionFinished) CALL FlagError("Data projection has not been finished.",err,error,*999)
    
    maximumIterationUpdate=dataProjection%maximumIterationUpdate       
     
    EXITS("DataProjection_MaximumInterationUpdateGet")
    RETURN
999 ERRORSEXITS("DataProjection_MaximumInterationUpdateGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_MaximumInterationUpdateGet

  !
  !================================================================================================================================
  !

  !>Sets the maximum iteration update for a data projection.
  SUBROUTINE DataProjection_MaximumInterationUpdateSet(dataProjection,maximumIterationUpdate,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to set the maximum iteration update for
    REAL(DP), INTENT(IN) :: maximumIterationUpdate !<the maximum iteration update to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DataProjection_MaximumInterationUpdateSet",err,error,*999)

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(dataProjection%dataProjectionFinished) CALL FlagError("Data projection has been finished.",err,error,*999)
    IF((maximumIterationUpdate<0.1_DP).OR.(maximumIterationUpdate>1.0_DP)) THEN
      localError="The specified maximum iteration update of "//TRIM(NumberToVString(maximumIterationUpdate,"*",err,error))// &
        & " is invalid. The value must be >= 0.1 and <= 1.0."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    dataProjection%maximumIterationUpdate=maximumIterationUpdate
     
    EXITS("DataProjection_MaximumInterationUpdateSet")
    RETURN
999 ERRORSEXITS("DataProjection_MaximumInterationUpdateSet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_MaximumInterationUpdateSet
  

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

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionFinished) CALL FlagError("Data projection has not been finished.",err,error,*999)
    
    maximumNumberOfIterations=dataProjection%maximumNumberOfIterations       
    
    EXITS("DataProjection_MaximumNumberOfIterationsGet")
    RETURN
999 ERRORSEXITS("DataProjection_MaximumNumberOfIterationsGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_MaximumNumberOfIterationsGet

  !
  !================================================================================================================================
  !
  
  !>Sets the maximum number of iterations for a data projection.
  SUBROUTINE DataProjection_MaximumNumberOfIterationsSet(dataProjection,maximumNumberOfIterations,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to set the maximum number of iterations for
    INTEGER(INTG), INTENT(IN) :: maximumNumberOfIterations !<the maximum number of iterations to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DataProjection_MaximumNumberOfIterationsSet",err,error,*999)

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(dataProjection%dataProjectionFinished) CALL FlagError("Data projection has been finished.",err,error,*999)
    IF(maximumNumberOfIterations<1) THEN
      localError="The specified maximum number of iterations of "// &
        & TRIM(NumberToVString(maximumNumberOfIterations,"*",err,error))//" is invalid. The value must be >= 1."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    dataProjection%maximumNumberOfIterations=maximumNumberOfIterations
    
    EXITS("DataProjection_MaximumNumberOfIterationsSet")
    RETURN
999 ERRORSEXITS("DataProjection_MaximumNumberOfIterationsSet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_MaximumNumberOfIterationsSet
  
  !
  !================================================================================================================================
  !
  
  !>Find the projection of a data point onto 1D elements
  SUBROUTINE DataProjection_NewtonElementsEvaluate_1(dataProjection,interpolatedPoint,dataPointLocation,candidateElements, &
    & projectionExitTag,projectionElementNumber,projectionDistance,projectionXi,projectionVector,err,error,*)
    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<Data projection problem to evaluate
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: interpolatedPoint !<The interpolated point for the data point
    REAL(DP), INTENT(IN) :: dataPointLocation(:) !<dataPointLocation(componentIdx). The location of data point in the region.
    INTEGER(INTG), INTENT(IN) :: candidateElements(:) !<candidateElememnts(candidateIdx). The list of candidate elements for the projection.
    INTEGER(INTG), INTENT(OUT) :: projectionExitTag !<On exit, the exit tag status for the data point projection
    INTEGER(INTG), INTENT(OUT) :: projectionElementNumber !<On exit, the element number of the data point projection
    REAL(DP), INTENT(OUT) :: projectionDistance !<On exit, the projection distance
    REAL(DP), INTENT(OUT) :: projectionXi(1) !<projectionXi(xiIdx). On exit, the xi location of the data point projection
    REAL(DP), INTENT(OUT) :: projectionVector(3) !<projectionVector(componentIdx). On exit, the projection vector for the data point projection
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string   
    !Local Variables
    LOGICAL :: insideRegion,converged
    INTEGER(INTG) :: elementNumber !local element number in current computational domain
    INTEGER(INTG) :: numberOfCoordinates
    INTEGER(INTG) :: bound,exitTag
    REAL(DP) :: xi(1),newXi(1),updateXi(1),updateXiNorm !<xi
    REAL(DP) :: distanceVector(3),relativeTolerance,absoluteTolerance !<tolerances
    REAL(DP) :: functionValue,newFunctionValue
    REAL(DP) :: functionGradient,functionHessian
    REAL(DP) :: maximumDelta,minimumDelta,delta !<trust region size
    REAL(DP) :: predictedReduction,predictionAccuracy
    
    INTEGER(INTG) :: elementIdx,iterationIdx1,iterationIdx2
    
    ENTERS("DataProjection_NewtonElementsEvaluate_1",err,error,*999)
              
    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionFinished) CALL FlagError("Data projection has not been finished.",err,error,*999)
    
    projectionExitTag=DATA_PROJECTION_EXIT_TAG_NO_ELEMENT
    numberOfCoordinates=dataProjection%numberOfCoordinates
    relativeTolerance=dataProjection%relativeTolerance
    absoluteTolerance=dataProjection%absoluteTolerance
    maximumDelta=dataProjection%maximumIterationUpdate
    minimumDelta=0.025_DP*maximumDelta !need to set a minimum, in case if it gets too small      
    !Project on each candidate elements
    DO elementIdx=1,SIZE(candidateElements,1) 
      elementNumber=candidateElements(elementIdx)
      IF(elementNumber>0) THEN
        exitTag=DATA_PROJECTION_EXIT_TAG_NO_ELEMENT
        converged=.FALSE.
        !Start at half the maximumDelta as we do not know if quadratic model is a good approximation yet
        delta=0.5_DP*maximumDelta 
        CALL Field_InterpolationParametersElementGet(dataProjection%projectionSetType,elementNumber,interpolatedPoint% &
          & INTERPOLATION_PARAMETERS,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
        xi=dataProjection%startingXi
        CALL Field_InterpolateXi(SECOND_PART_DERIV,xi,interpolatedPoint,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
        distanceVector(1:numberOfCoordinates)=interpolatedPoint%values(:,NO_PART_DERIV)-dataPointLocation
        functionValue=DOT_PRODUCT(distanceVector(1:numberOfCoordinates),distanceVector(1:numberOfCoordinates))       
        main_loop: DO iterationIdx1=1,dataProjection%maximumNumberOfIterations !(outer loop)
          !Check for bounds [0,1]
          IF(ABS(xi(1))<ZERO_TOLERANCE) THEN
            bound=-1 !bound at negative direction             
          ELSEIF(ABS(xi(1)-1.0_DP)<ZERO_TOLERANCE) THEN
            bound=1 !bound at positive direction
          ELSE !inside the bounds
            bound=0
          ENDIF
          !functionGradient 
          functionGradient=2.0_DP* &
            & (DOT_PRODUCT(distanceVector(1:numberOfCoordinates),interpolatedPoint%values(:,FIRST_PART_DERIV)))
          !functionHessian 
          functionHessian=2.0_DP* &
            & (DOT_PRODUCT(distanceVector(1:numberOfCoordinates),interpolatedPoint%values(:,SECOND_PART_DERIV))- &
            & DOT_PRODUCT(interpolatedPoint%values(:,FIRST_PART_DERIV),interpolatedPoint%values(:,FIRST_PART_DERIV)))
          !A model trust region approach, directly taken from CMISS CLOS22: V = -(H + EIGEN_SHIFT*I)g
          !The calculation of EIGEN_SHIFT are only approximated as opposed to the common trust region approach               
          !(inner loop: adjust region size) usually EXIT at 1 or 2 iterations
          DO iterationIdx2=1,dataProjection%maximumNumberOfIterations 
            insideRegion=.FALSE.
            IF(functionHessian>absoluteTolerance) THEN !positive: minimum exists
              updateXi(1)=-functionGradient/functionHessian
              updateXiNorm=DABS(updateXi(1))
              insideRegion=updateXiNorm<=delta
            ENDIF !positive                 
            IF(.NOT.insideRegion) THEN !minimum not in the region
              updateXi(1)=-DSIGN(delta,functionGradient)
              updateXiNorm=delta
            ENDIF
            IF((bound/=0).AND.(bound>0.EQV.updateXi(1)>0.0_DP)) THEN !projection go out of element bound
              exitTag=DATA_PROJECTION_EXIT_TAG_BOUNDS
              EXIT main_loop
            ENDIF
            converged=updateXiNorm<absoluteTolerance !first half of the convergence test (before collision detection)
            newXi=xi+updateXi !update xi
            IF(newXi(1)<0.0_DP) THEN !boundary collision check
              newXi(1)=0.0_DP
            ELSEIF(newXi(1)>1.0_DP) THEN
              newXi(1)=1.0_DP  
            ENDIF
            CALL Field_InterpolateXi(SECOND_PART_DERIV,newXi,interpolatedPoint,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
            distanceVector=interpolatedPoint%values(:,NO_PART_DERIV)-dataPointLocation
            newFunctionValue=DOT_PRODUCT(distanceVector(1:numberOfCoordinates),distanceVector(1:numberOfCoordinates))
            !second half of the convergence test
            converged=converged.AND.(DABS(newFunctionValue-functionValue)/(1.0_DP+functionValue)<relativeTolerance) 
            IF(converged) EXIT !converged: exit inner loop first
            IF((newFunctionValue-functionValue)>absoluteTolerance) THEN !bad model: reduce step size
              IF(delta<=minimumDelta) THEN
                !Something went wrong, minimumDelta too large? not likely to happen if minimumDelta is small
                exitTag=DATA_PROJECTION_EXIT_TAG_MAX_ITERATION !it will get stucked!!
                EXIT main_loop
              ENDIF
              delta=DMAX1(minimumDelta,0.25_DP*delta)
            ELSE
              predictedReduction=updateXi(1)*(functionGradient+0.5_DP*functionHessian*updateXi(1))
              predictionAccuracy=(newFunctionValue-functionValue)/predictedReduction
              IF(predictionAccuracy<0.01_DP) THEN !bad model: reduce region size
                IF(delta<=minimumDelta) THEN
                  !Something went wrong, minimumDelta too large? not likely to happen if minimumDelta is small
                  exitTag=DATA_PROJECTION_EXIT_TAG_MAX_ITERATION !it will get stucked!!
                  EXIT main_loop
                ENDIF
                delta=DMAX1(minimumDelta,0.5_DP*delta)
              ELSEIF(predictionAccuracy>0.9_DP.AND.predictionAccuracy<1.1_DP) THEN
                !good model: increase region size
                delta=DMIN1(maximumDelta,2.0_DP*delta)
                EXIT
              ELSE
                !ok model: keep the current region size
                EXIT
              ENDIF
            ENDIF
          ENDDO !iterationIdx2 (inner loop: adjust region size)
          functionValue=newFunctionValue
          xi=newXi
          IF(converged) THEN
            exitTag=DATA_PROJECTION_EXIT_TAG_CONVERGED
            EXIT
          ENDIF
        ENDDO main_loop !iterationIdx1 (outer loop)
        IF(exitTag==DATA_PROJECTION_EXIT_TAG_NO_ELEMENT.AND.iterationIdx1>=dataProjection%maximumNumberOfIterations) &
          & exitTag=DATA_PROJECTION_EXIT_TAG_MAX_ITERATION
        IF((projectionExitTag==DATA_PROJECTION_EXIT_TAG_NO_ELEMENT).OR.(DSQRT(functionValue)<projectionDistance)) THEN
          projectionExitTag=exitTag
          projectionElementNumber=elementNumber
          projectionDistance=DSQRT(functionValue)
          projectionXi=xi
          projectionVector=distanceVector
        ENDIF
      ENDIF
    ENDDO !elementIdx
    
    EXITS("DataProjection_NewtonElementsEvaluate_1")
    RETURN
999 ERRORSEXITS("DataProjection_NewtonElementsEvaluate_1",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_NewtonElementsEvaluate_1
  
  !
  !================================================================================================================================
  !
  
  !>Find the projection of a data point onto 2D elements
  SUBROUTINE DataProjection_NewtonElementsEvaluate_2(dataProjection,interpolatedPoint,dataPointLocation,candidateElements, &
    & projectionExitTag,projectionElementNumber,projectionDistance,projectionXi,projectionVector,err,error,*)
    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<Data projection problem to evaluate
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: interpolatedPoint !<The interpolated point for the data point 
    REAL(DP), INTENT(IN) :: dataPointLocation(:) !<dataPointLocation(componentIdx). The location of data point in the region.
    INTEGER(INTG), INTENT(IN) :: candidateElements(:) !<candidateElememnts(candidateIdx). The list of candidate elements for the projection.
    INTEGER(INTG), INTENT(OUT) :: projectionExitTag !<On exit, the exit tag status for the data point projection
    INTEGER(INTG), INTENT(OUT) :: projectionElementNumber !<On exit, the element number of the data point projection
    REAL(DP), INTENT(OUT) :: projectionDistance !<On exit, the projection distance
    REAL(DP), INTENT(OUT) :: projectionXi(2) !<projectionXi(xiIdx). On exit, the xi location of the data point projection
    REAL(DP), INTENT(OUT) :: projectionVector(3) !<projectionVector(componentIdx). On exit, the projection vector for the data point projection
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    LOGICAL :: free,converged,insideRegion
    INTEGER(INTG) :: elementNumber
    INTEGER(INTG) :: meshComponentNumber,numberOfCoordinates
    INTEGER(INTG) :: bound(2),exitTag
    REAL(DP) :: xi(2),newXi(2),updateXi(2),updateXiNorm !<xi
    REAL(DP) :: distanceVector(3),relativeTolerance,absoluteTolerance !<tolerances
    REAL(DP) :: functionValue,newFunctionValue
    REAL(DP) :: functionGradient(2),functionGradientNorm
    REAL(DP) :: functionHessian(2,2),hessianDiagonal(2)
    REAL(DP) :: temp1,temp2,determinant,minEigen,maxEigen,eigenShift
    REAL(DP) :: maximumDelta,minimumDelta,delta !<trust region size
    REAL(DP) :: predictedReduction,predictionAccuracy
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: domainMapping    
    
    INTEGER(INTG) :: elementIdx,xiIdx,fixedXiIdx,iterationIdx1,iterationIdx2
    
    ENTERS("DataProjection_NewtonElementsEvaluate_2",err,error,*999)
              
    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionFinished) CALL FlagError("Data projection has not been finished.",err,error,*999)
    
    projectionExitTag=DATA_PROJECTION_EXIT_TAG_NO_ELEMENT
    meshComponentNumber=interpolatedPoint%INTERPOLATION_PARAMETERS%field%decomposition%MESH_COMPONENT_NUMBER
    domainMapping=>interpolatedPoint%INTERPOLATION_PARAMETERS%field%decomposition%domain(meshComponentNumber)%ptr% &
      & mappings%elements
    numberOfCoordinates=dataProjection%numberOfCoordinates
    relativeTolerance=dataProjection%relativeTolerance
    absoluteTolerance=dataProjection%absoluteTolerance
    maximumDelta=dataProjection%maximumIterationUpdate
    minimumDelta=0.025_DP*maximumDelta !need to set a minimum, in case if it gets too small      
    DO elementIdx=1,SIZE(candidateElements,1) !project on each candidate elements
      elementNumber=candidateElements(elementIdx)
      IF(elementNumber>0) THEN
        exitTag=DATA_PROJECTION_EXIT_TAG_NO_ELEMENT
        converged=.FALSE.
        !start at half the maximumDelta as we do not know if quadratic model is a good approximation yet                      
        delta=0.5_DP*maximumDelta 
        CALL Field_InterpolationParametersElementGet(dataProjection%projectionSetType,elementNumber, &
          & interpolatedPoint%INTERPOLATION_PARAMETERS,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
        xi=dataProjection%startingXi
        CALL Field_InterpolateXi(SECOND_PART_DERIV,xi,interpolatedPoint,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
        distanceVector(1:numberOfCoordinates)=interpolatedPoint%values(1:numberOfCoordinates,NO_PART_DERIV)- &
          & dataPointLocation(1:numberOfCoordinates)
        functionValue=DOT_PRODUCT(distanceVector(1:numberOfCoordinates),distanceVector(1:numberOfCoordinates))
        !Outer loop
        main_loop: DO iterationIdx1=1,dataProjection%maximumNumberOfIterations 
          !Check for bounds [0,1]
          DO xiIdx=1,2 
            IF(ABS(xi(xiIdx))<ZERO_TOLERANCE) THEN
              bound(xiIdx)=-1 !bound at negative direction             
            ELSEIF(ABS(xi(xiIdx)-1.0_DP)<ZERO_TOLERANCE) THEN
              bound(xiIdx)=1 !bound at positive direction
            ELSE !inside the bounds
              bound(xiIdx)=0
            ENDIF
          ENDDO !xiIdx              
          !functionGradient 
          functionGradient(1)=2.0_DP*(DOT_PRODUCT(distanceVector(1:numberOfCoordinates), &
            & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S1)))
          functionGradient(2)=2.0_DP*(DOT_PRODUCT(distanceVector(1:numberOfCoordinates), &
            & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S2)))
          !functionHessian 
          functionHessian(1,1)=2.0_DP*(DOT_PRODUCT(distanceVector(1:numberOfCoordinates), &
            & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S1_S1))- &
            & DOT_PRODUCT(interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S1), &
            & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S1))) 
          functionHessian(1,2)=2.0_DP*(DOT_PRODUCT(distanceVector(1:numberOfCoordinates), &
            & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S1_S2))- &         
            & DOT_PRODUCT(interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S1), &
            & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S2)))
          functionHessian(2,2)=2.0_DP*(DOT_PRODUCT(distanceVector(1:numberOfCoordinates), &
            & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S2_S2))- &
            & DOT_PRODUCT(interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S2), &
            & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S2)))
          !A model trust region approach, a Newton step is taken if the minimum lies inside the trust region (delta),
          !if not, shift the step towards the steepest descent
          temp1=0.5_DP*(functionHessian(1,1)+functionHessian(2,2))
          temp2=DSQRT((0.5_DP*(functionHessian(1,1)-functionHessian(2,2)))**2+functionHessian(1,2)**2)
          minEigen=temp1-temp2
          maxEigen=temp1+temp2
          functionGradientNorm=DSQRT(DOT_PRODUCT(functionGradient,functionGradient))
          !Inner loop: adjust region size, usually EXIT at 1 or 2 iterations
          DO iterationIdx2=1,dataProjection%maximumNumberOfIterations 
            temp1=functionGradientNorm/delta
            !Estimate if the solution is inside the trust region without calculating a Newton step, this also guarantees
            !the Hessian matrix is positive definite
            insideRegion=(minEigen>=temp1).AND.(minEigen>absoluteTolerance) 
            IF(insideRegion) THEN
              determinant=minEigen*maxEigen !det(H)
              hessianDiagonal(1)=functionHessian(1,1)
              hessianDiagonal(2)=functionHessian(2,2)
            ELSE
              eigenShift=MAX(temp1-minEigen,absoluteTolerance) !shift towards steepest decent
              determinant=temp1*(maxEigen+eigenShift) !det(H)
              hessianDiagonal(1)=functionHessian(1,1)+eigenShift
              hessianDiagonal(2)=functionHessian(2,2)+eigenShift
            ENDIF
            updateXi(1)=-(hessianDiagonal(2)*functionGradient(1)-functionHessian(1,2)*functionGradient(2))/determinant
            updateXi(2)=(functionHessian(1,2)*functionGradient(1)-hessianDiagonal(1)*functionGradient(2))/determinant
            updateXiNorm=DSQRT(DOT_PRODUCT(updateXi,updateXi))
            free=.TRUE.
            DO xiIdx=1,2
              IF((bound(xiIdx)/=0).AND.(bound(xiIdx)>0.EQV.updateXi(xiIdx)>0.0_DP)) THEN
                !Projection will go out of element bounds
                IF(.NOT.free) THEN !both xi are fixed
                  exitTag=DATA_PROJECTION_EXIT_TAG_BOUNDS
                  EXIT main_loop
                ENDIF
                free=.FALSE.
                fixedXiIdx=xiIdx
              ENDIF
            ENDDO !xiIdx
            IF(free) THEN
              !Both xi are free
              IF(.NOT.insideRegion) THEN
                IF(updateXiNorm>0.0_DP) THEN
                  updateXi=delta/updateXiNorm*updateXi !readjust updateXi to lie on the region bound                      
                ENDIF
              ENDIF
            ELSE
              !xi are not free
              updateXi(fixedXiIdx)=0.0_DP
              xiIdx=3-fixedXiIdx
              insideRegion=.FALSE.
              IF(functionHessian(xiIdx,xiIdx)>0.0_DP) THEN
                !Positive: minimum exists in the unbounded direction                
                updateXi(xiIdx)=-functionGradient(xiIdx)/functionHessian(xiIdx,xiIdx)
                updateXiNorm=DABS(updateXi(xiIdx))
                insideRegion=updateXiNorm<=delta
              ENDIF
              IF(.NOT.insideRegion) THEN
                !Minimum not in the region
                updateXi(xiIdx)=-DSIGN(delta,functionGradient(xiIdx))
                updateXiNorm=delta
              ENDIF
            ENDIF
            !First half of the convergence test
            converged=updateXiNorm<absoluteTolerance 
            newXi=xi+updateXi !update xi
            DO xiIdx=1,2
              IF(newXi(xiIdx)<0.0_DP) THEN !boundary collision check
                newXi(xiIdx)=0.0_DP
                newXi(3-xiIdx)=xi(3-xiIdx)-updateXi(3-xiIdx)*xi(xiIdx)/updateXi(xiIdx)
              ELSEIF(newXi(xiIdx)>1.0_DP) THEN
                newXi(xiIdx)=1.0_DP  
                newXi(3-xiIdx)=xi(3-xiIdx)+updateXi(3-xiIdx)*(1.0_DP-xi(xiIdx))/updateXi(xiIdx)
              ENDIF
            ENDDO
            CALL Field_InterpolateXi(SECOND_PART_DERIV,newXi,interpolatedPoint,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
            distanceVector(1:numberOfCoordinates)=interpolatedPoint%values(:,NO_PART_DERIV)-dataPointLocation
            newFunctionValue=DOT_PRODUCT(distanceVector(1:numberOfCoordinates),distanceVector(1:numberOfCoordinates))
            !Second half of the convergence test (before collision detection)
            converged=converged.AND.(DABS(newFunctionValue-functionValue)/(1.0_DP+functionValue)<relativeTolerance) 
            IF(converged) EXIT !converged: exit inner loop first
            IF((newFunctionValue-functionValue)>absoluteTolerance) THEN
              !Bad model: reduce step size
              IF(delta<=minimumDelta) THEN
                !Something went wrong, minimumDelta too large? not likely to happen if minimumDelta is small
                exitTag=DATA_PROJECTION_EXIT_TAG_MAX_ITERATION !it will get stucked!!
                EXIT main_loop
              ENDIF
              delta=DMAX1(minimumDelta,0.25_DP*delta)
            ELSE
              predictedReduction=DOT_PRODUCT(functionGradient,updateXi)+ &
                & 0.5_DP*(updateXi(1)*(updateXi(1)*functionHessian(1,1)+2.0_DP*updateXi(2)*functionHessian(1,2))+ &
                & updateXi(2)**2*functionHessian(2,2))
              predictionAccuracy=(newFunctionValue-functionValue)/predictedReduction
              IF(predictionAccuracy<0.01_DP) THEN
                !Bad model: reduce region size
                IF(delta<=minimumDelta) THEN
                  !Something went wrong, minimumDelta too large? not likely to happen if minimumDelta is small
                  exitTag=DATA_PROJECTION_EXIT_TAG_MAX_ITERATION !it will get stucked!!
                  EXIT main_loop
                ENDIF
                delta=DMAX1(minimumDelta,0.5_DP*delta)
              ELSEIF(predictionAccuracy>0.9_DP.AND.predictionAccuracy<1.1_DP) THEN
                !Good model: increase region size
                delta=DMIN1(maximumDelta,2.0_DP*delta)
                EXIT
              ELSE
                !OK model: keep the current region size
                EXIT
              ENDIF
            ENDIF
          ENDDO !iterationIdx2 (inner loop: adjust region size)
          functionValue=newFunctionValue
          xi=newXi
          IF(converged) THEN
            exitTag=DATA_PROJECTION_EXIT_TAG_CONVERGED
            EXIT
          ENDIF
        ENDDO main_loop !iterationIdx1 (outer loop)
        IF(exitTag==DATA_PROJECTION_EXIT_TAG_NO_ELEMENT.AND.iterationIdx1>=dataProjection%maximumNumberOfIterations) &
          & exitTag=DATA_PROJECTION_EXIT_TAG_MAX_ITERATION
        IF((projectionExitTag==DATA_PROJECTION_EXIT_TAG_NO_ELEMENT).OR.(DSQRT(functionValue)<projectionDistance)) THEN
          projectionExitTag=exitTag
          projectionElementNumber=elementNumber
          projectionDistance=DSQRT(functionValue)
          projectionXi=xi
          projectionVector=distanceVector
        ENDIF
      ENDIF
    ENDDO !elementIdx
    
    EXITS("DataProjection_NewtonElementsEvaluate_2")
    RETURN
999 ERRORSEXITS("DataProjection_NewtonElementsEvaluate_2",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_NewtonElementsEvaluate_2

  !
  !================================================================================================================================
  !
  
  !>Find the projection of a data point onto 3D elements
  SUBROUTINE DataProjection_NewtonElementsEvaluate_3(dataProjection,interpolatedPoint,dataPointLocation,candidateElements, &
    & projectionExitTag,projectionElementNumber,projectionDistance,projectionXi,projectionVector,err,error,*)
    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<Data projection problem to evaluate
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: interpolatedPoint !<The interpolated point for the data point 
    REAL(DP), INTENT(IN) :: dataPointLocation(:) !<dataPointLocation(componentIdx). The location of data point in the region.
    INTEGER(INTG), INTENT(IN) :: candidateElements(:) !<candidateElememnts(candidateIdx). The list of candidate elements for the projection.
    INTEGER(INTG), INTENT(OUT) :: projectionExitTag !<On exit, the exit tag status for the data point projection
    INTEGER(INTG), INTENT(OUT) :: projectionElementNumber !<On exit, the element number of the data point projection
    REAL(DP), INTENT(OUT) :: projectionDistance !<On exit, the projection distance
    REAL(DP), INTENT(OUT) :: projectionXi(3) !<projectionXi(xiIdx). On exit, the xi location of the data point projection
    REAL(DP), INTENT(OUT) :: projectionVector(3) !<projectionVector(componentIdx). On exit, the projection vector for the data point projection
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
   
    !Local Variables
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: domainMapping
    LOGICAL :: free,converged,insideRegion
    INTEGER(INTG) :: elementNumber
    INTEGER(INTG) :: meshComponentNumber,numberOfCoordinates
    INTEGER(INTG) :: nBound,bound(3),exitTag
    REAL(DP) :: xi(3),newXi(3),updateXi(3),updateXiNorm !<xi
    REAL(DP) :: distanceVector(3),relativeTolerance,absoluteTolerance !<tolerances
    REAL(DP) :: functionValue,newFunctionValue
    REAL(DP) :: functionGradient(3),functionGradientNorm,functionGradient2(2)
    REAL(DP) :: functionHessian(3,3),hessianDiagonal(3),functionHessian2(2,2),hessianDiagonal2(2)
    REAL(DP) :: temp1,temp2,temp3,temp4,determinant,trace,trace2,minEigen,maxEigen,eigenShift    
    REAL(DP) :: maximumDelta,minimumDelta,delta !<trust region size
    REAL(DP) :: predictedReduction,predictionAccuracy
    INTEGER(INTG) :: elementIdx,xiIdx,faceXiIdxs(2),boundXiIdx,fixedXiIdx,fixedXiIdx2(2),iterationIdx1,iterationIdx2, &
      & fixedBoundXiIdx,xiIdx2
    
    ENTERS("DataProjection_NewtonElementsEvaluate_3",err,error,*999)
              
    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionFinished) CALL FlagError("Data projection has not been finished.",err,error,*999)
    
    projectionExitTag=DATA_PROJECTION_EXIT_TAG_NO_ELEMENT
    meshComponentNumber=interpolatedPoint%INTERPOLATION_PARAMETERS%field%decomposition%MESH_COMPONENT_NUMBER
    domainMapping=>interpolatedPoint%INTERPOLATION_PARAMETERS%field%decomposition%domain(meshComponentNumber)%ptr% &
      & mappings%elements
    numberOfCoordinates=dataProjection%numberOfCoordinates
    relativeTolerance=dataProjection%relativeTolerance
    absoluteTolerance=dataProjection%absoluteTolerance
    maximumDelta=dataProjection%maximumIterationUpdate
    minimumDelta=0.025_DP*maximumDelta !need to set a minimum, in case if it gets too small
    !Project on each candidate elements
    DO elementIdx=1,SIZE(candidateElements,1)
      elementNumber=candidateElements(elementIdx)
      IF(elementNumber>0) THEN
        exitTag=DATA_PROJECTION_EXIT_TAG_NO_ELEMENT
        converged=.FALSE.
        !start at half the maximumDelta as we do not know if quadratic model is a good approximation yet            
        delta=0.5_DP*maximumDelta 
        CALL Field_InterpolationParametersElementGet(dataProjection%projectionSetType,elementNumber, &
          & interpolatedPoint%INTERPOLATION_PARAMETERS,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
        xi=dataProjection%startingXi
        CALL Field_InterpolateXi(SECOND_PART_DERIV,xi,interpolatedPoint,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
        distanceVector(1:numberOfCoordinates)=interpolatedPoint%values(1:numberOfCoordinates,NO_PART_DERIV)- &
          & dataPointLocation(1:numberOfCoordinates)
        functionValue=DOT_PRODUCT(distanceVector,distanceVector)
        !Outer loop
        main_loop: DO iterationIdx1=1,dataProjection%maximumNumberOfIterations 
          !Check for bounds [0,1]
          DO xiIdx=1,3
            IF(ABS(xi(xiIdx))<ZERO_TOLERANCE) THEN
              bound(xiIdx)=-1 !bound at negative direction             
            ELSEIF(ABS(xi(xiIdx)-1.0_DP)<ZERO_TOLERANCE) THEN
              bound(xiIdx)=1 !bound at positive direction
            ELSE 
              bound(xiIdx)=0 !inside the bounds
            ENDIF
          ENDDO !xiIdx              
          !functionGradient 
          functionGradient(1)=2.0_DP*(DOT_PRODUCT(distanceVector(1:numberOfCoordinates), &
            & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S1)))
          functionGradient(2)=2.0_DP*(DOT_PRODUCT(distanceVector(1:numberOfCoordinates), &
            & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S2)))
          functionGradient(3)=2.0_DP*(DOT_PRODUCT(distanceVector(1:numberOfCoordinates), &
            & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S3)))
          !functionHessian 
          functionHessian(1,1)=2.0_DP*(DOT_PRODUCT(distanceVector(1:numberOfCoordinates), &
            & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S1_S1))- &
            & DOT_PRODUCT(interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S1), &
            & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S1))) 
          functionHessian(1,2)=2.0_DP*(DOT_PRODUCT(distanceVector(1:numberOfCoordinates), &
            & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S1_S2))- &         
            & DOT_PRODUCT(interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S1), &
            & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S2)))
          functionHessian(1,3)=2.0_DP*(DOT_PRODUCT(distanceVector(1:numberOfCoordinates), &
            & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S1_S3))- &         
            & DOT_PRODUCT(interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S1), &
            & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S3)))
          functionHessian(2,2)=2.0_DP*(DOT_PRODUCT(distanceVector(1:numberOfCoordinates), &
            & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S2_S2))- &
            & DOT_PRODUCT(interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S2), &
            & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S2)))
          functionHessian(2,3)=2.0_DP*(DOT_PRODUCT(distanceVector(1:numberOfCoordinates), &
            & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S2_S3))- &
            & DOT_PRODUCT(interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S2), &
            & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S3)))
          functionHessian(3,3)=2.0_DP*(DOT_PRODUCT(distanceVector(1:numberOfCoordinates), &
            & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S3_S3))- &
            & DOT_PRODUCT(interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S3), &
            & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S3)))    
          !A model trust region approach, a Newton step is taken if the solution lies inside the trust region (delta),
          !if not, shift the step towards the steepest descent
          trace=functionHessian(1,1)+functionHessian(2,2)+functionHessian(3,3) !tr(H)
          trace2=functionHessian(1,1)*functionHessian(2,2)+functionHessian(1,1)*functionHessian(3,3)+ &
            & functionHessian(2,2)*functionHessian(3,3)-functionHessian(1,2)**2-functionHessian(1,3)**2- &
            & functionHessian(2,3)**2 !tr(H**2)-(tr(H))**2
          determinant=functionHessian(1,1)*functionHessian(2,2)*functionHessian(3,3)- &
            & functionHessian(1,1)*functionHessian(2,3)**2-functionHessian(2,2)*functionHessian(1,3)**2- &
            & functionHessian(3,3)*functionHessian(1,2)**2+ &
            & 2.0_DP*functionHessian(1,2)*functionHessian(1,3)*functionHessian(2,3) !det(H)                 
          temp1=-trace/3.0_DP
          temp2=trace2/3.0_DP
          temp3=temp2-temp1**2 !<=0
          IF(temp3>-1.0E-5_DP) THEN !include some negatives for numerical errors
            minEigen=-temp1 !all eigenvalues are the same                
          ELSE
            temp3=DSQRT(-temp3)
            temp4=(determinant+3.0_DP*(temp1*temp2)-2.0_DP*temp1**3)/(2.0_DP*temp3**3)
            minEigen=2.0_DP*temp3*DCOS((DACOS(temp4)+TWOPI)/3.0_DP)-temp1                
          ENDIF
          functionGradientNorm=DSQRT(DOT_PRODUCT(functionGradient,functionGradient))
          !Inner loop: adjust region size, usually EXIT at 1 or 2 iterations
          DO iterationIdx2=1,dataProjection%maximumNumberOfIterations 
            temp1=functionGradientNorm/delta
            !Estimate if the solution is inside the trust region without calculating a Newton step, this also guarantees
            !the Hessian matrix is positive definite
            insideRegion=(minEigen>=temp1).AND.(minEigen>absoluteTolerance) 
            IF(insideRegion) THEN
              hessianDiagonal(1)=functionHessian(1,1)
              hessianDiagonal(2)=functionHessian(2,2)
              hessianDiagonal(3)=functionHessian(3,3)
            ELSE
              eigenShift=MAX(temp1-minEigen,absoluteTolerance) !shift towards steepest decent
              determinant=determinant+eigenShift*(trace2+eigenShift*(trace+eigenShift)) !shift the determinant
              hessianDiagonal(1)=functionHessian(1,1)+eigenShift
              hessianDiagonal(2)=functionHessian(2,2)+eigenShift
              hessianDiagonal(3)=functionHessian(3,3)+eigenShift
            ENDIF
            temp2=functionHessian(1,3)*functionHessian(2,3)-functionHessian(1,2)*hessianDiagonal(3)
            temp3=functionHessian(1,2)*functionHessian(2,3)-functionHessian(1,3)*hessianDiagonal(2)
            temp4=functionHessian(1,2)*functionHessian(1,3)-functionHessian(2,3)*hessianDiagonal(1)               
            updateXi(1)=((functionHessian(2,3)**2-hessianDiagonal(2)*hessianDiagonal(3))*functionGradient(1)- &
              & temp2*functionGradient(2)-temp3*functionGradient(3))/determinant
            updateXi(2)=((functionHessian(1,3)**2-hessianDiagonal(1)*hessianDiagonal(3))*functionGradient(2)- & 
              & temp2*functionGradient(1)-temp4*functionGradient(3))/determinant
            updateXi(3)=((functionHessian(1,2)**2-hessianDiagonal(1)*hessianDiagonal(2))*functionGradient(3)- &
              & temp3*functionGradient(1)-temp4*functionGradient(2))/determinant
            updateXiNorm=DSQRT(DOT_PRODUCT(updateXi,updateXi))
            free=.TRUE.
            nBound=0
            DO xiIdx=1,3
              IF((bound(xiIdx)/=0).AND.(bound(xiIdx)>0.EQV.updateXi(xiIdx)>0.0_DP)) THEN
                !Projection will go out of element bounds
                nBound=nBound+1
                free=.FALSE.
                IF(nBound<=2) THEN
                  fixedXiIdx2(nBound)=xiIdx
                ELSE
                  !All xi are fixed
                  exitTag=DATA_PROJECTION_EXIT_TAG_BOUNDS
                  EXIT main_loop
                ENDIF
              ENDIF
            ENDDO !xiIdx
            IF(free) THEN
              !All xi are free
              IF(.NOT.insideRegion) THEN
                IF(updateXiNorm>0.0_DP) THEN
                  updateXi=delta/updateXiNorm*updateXi !readjust updateXi to lie on the region bound                      
                ENDIF
              ENDIF
            ELSE
              !At least one of the xi are not free
              !Try 2D projection
              free=.TRUE.
              fixedXiIdx=fixedXiIdx2(1)
              IF(nBound==2) THEN
                !only fix the direction that is most strongly suggesting leaving the element
                IF(updateXi(fixedXiIdx2(2))>updateXi(fixedXiIdx2(1))) fixedXiIdx=fixedXiIdx2(2) 
              ENDIF
              updateXi(fixedXiIdx)=0.0_DP
              faceXiIdxs(1)=1+MOD(fixedXiIdx,3)
              faceXiIdxs(2)=1+MOD(fixedXiIdx+1,3)
              !functionGradient2
              functionGradient2(1)=functionGradient(faceXiIdxs(1))
              functionGradient2(2)=functionGradient(faceXiIdxs(2))
              !functionHessian2
              functionHessian2(1,1)=functionHessian(faceXiIdxs(1),faceXiIdxs(1))
              functionHessian2(1,2)=functionHessian(faceXiIdxs(1),faceXiIdxs(2))
              functionHessian2(2,2)=functionHessian(faceXiIdxs(2),faceXiIdxs(2))
              !re-estimate the trust solution in 2D
              temp1=0.5_DP*(functionHessian2(1,1)+functionHessian2(2,2))
              temp2=DSQRT((0.5_DP*(functionHessian2(1,1)-functionHessian2(2,2)))**2+functionHessian2(1,2)**2)
              minEigen=temp1-temp2
              maxEigen=temp1+temp2
              temp3=DSQRT(DOT_PRODUCT(functionGradient2,functionGradient2))/delta
              !Estimate if the solution is inside the trust region without calculating a Newton step, this also guarantees
              !the Hessian matrix is positive definite
              insideRegion=(minEigen>=temp3).AND.(minEigen>absoluteTolerance) 
              IF(insideRegion) THEN
                determinant=minEigen*maxEigen !determinant of functionHessian
                hessianDiagonal2(1)=functionHessian2(1,1)
                hessianDiagonal2(2)=functionHessian2(2,2)
              ELSE
                eigenShift=MAX(temp3-minEigen,absoluteTolerance) !shift towards steepest decent
                determinant=temp3*(maxEigen+eigenShift) !determinant of shifted functionHessian
                hessianDiagonal2(1)=functionHessian2(1,1)+eigenShift
                hessianDiagonal2(2)=functionHessian2(2,2)+eigenShift
              ENDIF
              updateXi(faceXiIdxs(1))=-(hessianDiagonal2(2)*functionGradient2(1)- &
                & functionHessian2(1,2)*functionGradient2(2))/determinant
              updateXi(faceXiIdxs(2))=(functionHessian2(1,2)*functionGradient2(1)- &
                & hessianDiagonal2(1)*functionGradient2(2))/determinant
              updateXiNorm=DSQRT(DOT_PRODUCT(updateXi,updateXi))
              !Check again for bounds
              DO boundXiIdx=1,2
                IF((bound(faceXiIdxs(boundXiIdx))/=0).AND.(bound(faceXiIdxs(boundXiIdx))>0.EQV. &
                  & updateXi(faceXiIdxs(boundXiIdx))>0.0_DP)) THEN
                  !Projection will go out of element bounds
                  IF(.NOT.free) THEN
                    !Both xi are fixed
                    exitTag=DATA_PROJECTION_EXIT_TAG_BOUNDS
                    EXIT main_loop
                  ENDIF
                  free=.FALSE.
                  fixedBoundXiIdx=boundXiIdx
                ENDIF
              ENDDO !xiIdx
              IF(free) THEN
                !Both xi are free
                IF(.NOT.insideRegion) THEN
                  IF(updateXiNorm>0.0_DP) THEN
                    updateXi=delta/updateXiNorm*updateXi !readjust updateXi to lie on the region bound                      
                  ENDIF
                ENDIF
              ELSE
                !xi are not free
                updateXi(faceXiIdxs(fixedBoundXiIdx))=0.0_DP
                xiIdx=faceXiIdxs(3-fixedBoundXiIdx)
                insideRegion=.FALSE.
                IF(functionHessian(xiIdx,xiIdx)>0.0_DP) THEN
                  !positive: minimum exists in the unbounded direction                
                  updateXi(xiIdx)=-functionGradient(xiIdx)/functionHessian(xiIdx,xiIdx)
                  updateXiNorm=DABS(updateXi(xiIdx))
                  insideRegion=updateXiNorm<=delta
                ENDIF
                IF(.NOT.insideRegion) THEN
                  !minimum not in the region
                  updateXi(xiIdx)=-DSIGN(delta,functionGradient(xiIdx))
                  updateXiNorm=delta
                ENDIF
              ENDIF !if xi are free (2D)
            ENDIF !if xi are free (3D)
            !First half of the convergence test
            converged=updateXiNorm<absoluteTolerance 
            newXi=xi+updateXi !update xi
            DO xiIdx=1,3
              IF(ABS(updateXi(xiIdx))<ZERO_TOLERANCE) THEN
                !FPE Handling
                IF(newXi(xiIdx)<0.0_DP) THEN
                  newXi(xiIdx)=0.0_DP
                  DO xiIdx2 = 1,3
                    IF(xiIdx2 /= xiIdx) THEN
                      newXi(xiIdx2)=xi(xiIdx2)-updateXi(xiIdx2)
                    ENDIF
                  ENDDO
                ELSEIF(newXi(xiIdx)>1.0_DP) THEN
                  newXi(xiIdx)=1.0_DP
                  DO xiIdx2 = 1,3
                    IF(xiIdx2 /= xiIdx) THEN
                      newXi(xiIdx2)=xi(xiIdx2)+updateXi(xiIdx2)
                    ENDIF
                  ENDDO
                ENDIF
              ELSE IF(newXi(xiIdx)<0.0_DP) THEN
                !boundary collision check
                newXi(xiIdx)=0.0_DP
                newXi=xi-updateXi*xi(xiIdx)/updateXi(xiIdx)
              ELSE IF(newXi(xiIdx)>1.0_DP) THEN
                newXi(xiIdx)=1.0_DP  
                newXi=xi+updateXi*(1.0_DP-xi(xiIdx))/updateXi(xiIdx)
              ENDIF
            ENDDO !xiIdx
            CALL Field_InterpolateXi(SECOND_PART_DERIV,newXi,interpolatedPoint,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
            distanceVector=interpolatedPoint%values(1:numberOfCoordinates,NO_PART_DERIV)-dataPointLocation
            newFunctionValue=DOT_PRODUCT(distanceVector,distanceVector)
            !Second half of the convergence test (before collision detection)
            converged=converged.AND.(DABS(newFunctionValue-functionValue)/(1.0_DP+functionValue)<relativeTolerance) 
            IF(converged) EXIT !converged: exit inner loop first
            IF((newFunctionValue-functionValue)>absoluteTolerance) THEN
              !bad model: reduce step size
              IF(delta<=minimumDelta) THEN
                !something went wrong, minimumDelta too large? not likely to happen if minimumDelta is small
                exitTag=DATA_PROJECTION_EXIT_TAG_MAX_ITERATION ! \it will get stucked!!
                EXIT main_loop
              ENDIF
              delta=DMAX1(minimumDelta,0.25_DP*delta)
            ELSE
              predictedReduction=DOT_PRODUCT(functionGradient,updateXi)+ &
                & 0.5_DP*(updateXi(1)*(updateXi(1)*functionHessian(1,1)+2.0_DP*updateXi(2)*functionHessian(1,2)+ &
                & 2.0_DP*updateXi(3)*functionHessian(1,3))+updateXi(2)*(updateXi(2)*functionHessian(2,2)+ &
                & 2.0_DP*updateXi(3)*functionHessian(2,3))+updateXi(2)**2*functionHessian(2,2))
              IF (ABS(predictedReduction) < ZERO_TOLERANCE) THEN
                converged = .TRUE.  !prediction reduction converged
                EXIT
              ENDIF
              predictionAccuracy=(newFunctionValue-functionValue)/predictedReduction
              IF(predictionAccuracy<0.01_DP) THEN 
                !bad model: reduce region size
                IF(delta<=minimumDelta) THEN
                  !something went wrong, minimumDelta too large? not likely to happen if minimumDelta is small
                  exitTag=DATA_PROJECTION_EXIT_TAG_MAX_ITERATION !it will get stucked!!
                  EXIT main_loop
                ENDIF
                delta=DMAX1(minimumDelta,0.5_DP*delta)
              ELSEIF(predictionAccuracy>0.9_DP.AND.predictionAccuracy<1.1_DP) THEN
                !good model: increase region size
                delta=DMIN1(maximumDelta,2.0_DP*delta)
                EXIT
              ELSE
                !ok model: keep the current region size
                EXIT
              ENDIF
            ENDIF
          ENDDO !iterationIdx2 (inner loop: adjust region size)
          functionValue=newFunctionValue
          xi=newXi
          IF(converged) THEN
            exitTag=DATA_PROJECTION_EXIT_TAG_CONVERGED
            EXIT
          ENDIF
        ENDDO main_loop !iterationIdx1 (outer loop)
        IF(exitTag==DATA_PROJECTION_EXIT_TAG_NO_ELEMENT.AND.iterationIdx1>=dataProjection%maximumNumberOfIterations) &
          & exitTag=DATA_PROJECTION_EXIT_TAG_MAX_ITERATION
        IF((projectionExitTag==DATA_PROJECTION_EXIT_TAG_NO_ELEMENT).OR.(DSQRT(functionValue)<projectionDistance)) THEN
          projectionExitTag=exitTag
          projectionElementNumber=elementNumber
          projectionDistance=DSQRT(functionValue)
          projectionXi=xi
          projectionVector=distanceVector
        ENDIF
      ENDIF
    ENDDO !elementIdx
    
    EXITS("DataProjection_NewtonElementsEvaluate_3")
    RETURN
999 ERRORSEXITS("DataProjection_NewtonElementsEvaluate_3",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_NewtonElementsEvaluate_3
  
  !
  !================================================================================================================================
  !

  !>Find the projection of a data point onto element faces (slight difference to DataProjection_NewtonElementsEvaluate_2)
  SUBROUTINE DataProjection_NewtonFacesEvaluate(dataProjection,interpolatedPoint,dataPointLocation,candidateElements, &
    & candidateElementFaces,projectionExitTag,projectionElementNumber,projectionElementFaceNumber,projectionDistance, &
    & projectionXi,projectionVector,err,error,*)
    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<Data projection problem to evaluate
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: interpolatedPoint !<The interpolation for the data point
    REAL(DP), INTENT(IN) :: dataPointLocation(:) !<dataPointLocation(componentIdx). The location of data point in the region.
    INTEGER(INTG), INTENT(IN) :: candidateElements(:)!<candidateElememnts(candidateIdx). The list of candidate elements for the projection.
    INTEGER(INTG), INTENT(IN) :: candidateElementFaces(:) !<candidateElementFaces(candidateIdx). The list of candidate faces for the projection.
    INTEGER(INTG), INTENT(OUT) :: projectionExitTag !<On exit, the exit tag status for the data point projection
    INTEGER(INTG), INTENT(OUT) :: projectionElementNumber !<On exit, the element number of the data point projection
    INTEGER(INTG), INTENT(OUT) :: projectionElementFaceNumber !<On exit, the face number of the data point projection
    REAL(DP), INTENT(OUT) :: projectionDistance !<On exit, the projection distance
    REAL(DP), INTENT(OUT) :: projectionXi(2) !<projectionXi(xiIdx). On exit, the xi location of the data point projection
    REAL(DP), INTENT(OUT) :: projectionVector(3) !<projectionVector(componentIdx). On exit, the projection vector for the data point projection
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string   
    !Local Variables
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: domainMapping
    LOGICAL :: free,converged,insideRegion
    INTEGER(INTG) :: elementNumber,elementFaceNumber,faceNumber
    INTEGER(INTG) :: meshComponentNumber,numberOfCoordinates
    INTEGER(INTG) :: bound(2),exitTag
    REAL(DP) :: xi(2),newXi(2),updateXi(2),updateXiNorm !<xi
    REAL(DP) :: distanceVector(3),relativeTolerance,absoluteTolerance !<tolerances
    REAL(DP) :: functionValue,newFunctionValue
    REAL(DP) :: functionGradient(2),functionGradientNorm
    REAL(DP) :: functionHessian(2,2),hessianDiagonal(2)
    REAL(DP) :: temp1,temp2,determinant,minEigen,maxEigen,eigenShift
    REAL(DP) :: maximumDelta,minimumDelta,delta !<trust region size
    REAL(DP) :: predictedReduction,predictionAccuracy       
    INTEGER(INTG) :: elementIdx,xiIdx,fixedXiIdx,iterationIdx1,iterationIdx2
    
    ENTERS("DataProjection_NewtonFacesEvaluate",err,error,*999)
              
    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionFinished) CALL FlagError("Data projection has not been finished.",err,error,*999)
    
    projectionExitTag=DATA_PROJECTION_EXIT_TAG_NO_ELEMENT
    meshComponentNumber=interpolatedPoint%INTERPOLATION_PARAMETERS%FIELD%decomposition%MESH_COMPONENT_NUMBER
    domainMapping=>interpolatedPoint%INTERPOLATION_PARAMETERS%FIELD%decomposition%domain(meshComponentNumber)% &
      & PTR%MAPPINGS%ELEMENTS
    numberOfCoordinates=dataProjection%numberOfCoordinates
    relativeTolerance=dataProjection%relativeTolerance
    absoluteTolerance=dataProjection%absoluteTolerance
    maximumDelta=dataProjection%maximumIterationUpdate
    minimumDelta=0.025_DP*maximumDelta !need to set a minimum, in case if it gets too small
    !Project on each candidate elements
    DO elementIdx=1,SIZE(candidateElements,1) 
      elementNumber=candidateElements(elementIdx)
      IF(elementNumber>0) THEN
        elementFaceNumber=candidateElementFaces(elementIdx)
        faceNumber=interpolatedPoint%INTERPOLATION_PARAMETERS%field%decomposition%topology%elements% &
          & elements(elementNumber)%ELEMENT_FACES(elementFaceNumber)     
        exitTag=DATA_PROJECTION_EXIT_TAG_NO_ELEMENT
        converged=.FALSE.
        !start at half the maximumDelta as we do not know if quadratic model is a good approximation yet
        delta=0.5_DP*maximumDelta 
        CALL Field_InterpolationParametersFaceGet(dataProjection%projectionSetType,faceNumber,interpolatedPoint% &
          & INTERPOLATION_PARAMETERS,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
        xi=dataProjection%startingXi
        CALL Field_InterpolateXi(SECOND_PART_DERIV,xi,interpolatedPoint,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
        distanceVector(1:numberOfCoordinates)=interpolatedPoint%values(1:numberOfCoordinates,NO_PART_DERIV)- &
          & dataPointLocation(1:numberOfCoordinates)
        functionValue=DOT_PRODUCT(distanceVector(1:numberOfCoordinates),distanceVector(1:numberOfCoordinates))
        !Outer loop
        main_loop: DO iterationIdx1=1,dataProjection%maximumNumberOfIterations
          !Check for bounds [0,1]
          DO xiIdx=1,2 
            IF(ABS(xi(xiIdx))<ZERO_TOLERANCE) THEN
              bound(xiIdx)=-1 !bound at negative direction             
            ELSEIF(ABS(xi(xiIdx)-1.0_DP)<ZERO_TOLERANCE) THEN
              bound(xiIdx)=1 !bound at positive direction
            ELSE 
              bound(xiIdx)=0 !inside the bounds
            ENDIF
          ENDDO !xiIdx              
          !functionGradient 
          functionGradient(1)=2.0_DP*(DOT_PRODUCT(distanceVector(1:numberOfCoordinates), &
            & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S1)))
          functionGradient(2)=2.0_DP*(DOT_PRODUCT(distanceVector(1:numberOfCoordinates), &
            & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S2)))
          !functionHessian 
          functionHessian(1,1)=2.0_DP*(DOT_PRODUCT(distanceVector(1:numberOfCoordinates), &
            & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S1_S1))- &
            & DOT_PRODUCT(interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S1), &
            & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S1))) 
          functionHessian(1,2)=2.0_DP*(DOT_PRODUCT(distanceVector(1:numberOfCoordinates), &
            & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S1_S2))- &         
            & DOT_PRODUCT(interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S1), &
            & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S2)))
          functionHessian(2,2)=2.0_DP*(DOT_PRODUCT(distanceVector(1:numberOfCoordinates), &
            & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S2_S2))- &
            & DOT_PRODUCT(interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S2), &
            & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S2)))
          !A model trust region approach, a Newton step is taken if the minimum lies inside the trust region (delta),
          !if not, shift the step towards the steepest descent
          temp1=0.5_DP*(functionHessian(1,1)+functionHessian(2,2))
          temp2=DSQRT((0.5_DP*(functionHessian(1,1)-functionHessian(2,2)))**2+functionHessian(1,2)**2)
          minEigen=temp1-temp2
          maxEigen=temp1+temp2
          functionGradientNorm=DSQRT(DOT_PRODUCT(functionGradient,functionGradient))
          !Inner loop: adjust region size, usually EXIT at 1 or 2 iterations
          DO iterationIdx2=1,dataProjection%maximumNumberOfIterations 
            temp1=functionGradientNorm/delta
            !Estimate if the solution is inside the trust region without calculating a Newton step, this also guarantees
            !the Hessian matrix is positive definite
            insideRegion=(minEigen>=temp1).AND.(minEigen>absoluteTolerance) 
            IF(insideRegion) THEN
              determinant=minEigen*maxEigen !det(H)
              hessianDiagonal(1)=functionHessian(1,1)
              hessianDiagonal(2)=functionHessian(2,2)
            ELSE
              eigenShift=MAX(temp1-minEigen,absoluteTolerance) !shift towards steepest decent
              determinant=temp1*(maxEigen+eigenShift) !det(H)
              hessianDiagonal(1)=functionHessian(1,1)+eigenShift
              hessianDiagonal(2)=functionHessian(2,2)+eigenShift
            ENDIF
            updateXi(1)=-(hessianDiagonal(2)*functionGradient(1)-functionHessian(1,2)*functionGradient(2))/determinant
            updateXi(2)=(functionHessian(1,2)*functionGradient(1)-hessianDiagonal(1)*functionGradient(2))/determinant
            updateXiNorm=DSQRT(DOT_PRODUCT(updateXi,updateXi))
            free=.TRUE.
            DO xiIdx=1,2
              IF((bound(xiIdx)/=0).AND.(bound(xiIdx)>0.EQV.updateXi(xiIdx)>0.0_DP)) THEN
                !Projection will go out of element bounds
                IF(.NOT.free) THEN !both xi are fixed
                  exitTag=DATA_PROJECTION_EXIT_TAG_BOUNDS
                  EXIT main_loop
                ENDIF
                free=.FALSE.
                fixedXiIdx=xiIdx
              ENDIF
            ENDDO !xiIdx
            IF(free) THEN
              !Both xi are free
              IF(.NOT.insideRegion) THEN
                IF(updateXiNorm>0.0_DP) THEN
                  updateXi=delta/updateXiNorm*updateXi !readjust updateXi to lie on the region bound                      
                ENDIF
              ENDIF
            ELSE
              !xi are not free
              updateXi(fixedXiIdx)=0.0_DP
              xiIdx=3-fixedXiIdx
              insideRegion=.FALSE.
              IF(functionHessian(xiIdx,xiIdx)>0.0_DP) THEN
                !Positive: minimum exists in the unbounded direction                
                updateXi(xiIdx)=-functionGradient(xiIdx)/functionHessian(xiIdx,xiIdx)
                updateXiNorm=DABS(updateXi(xiIdx))
                insideRegion=updateXiNorm<=delta
              ENDIF
              IF(.NOT.insideRegion) THEN
                !Minimum not in the region
                updateXi(xiIdx)=-DSIGN(delta,functionGradient(xiIdx))
                updateXiNorm=delta
              ENDIF
            ENDIF !if xi is free
            !First half of the convergence test
            converged=updateXiNorm<absoluteTolerance 
            newXi=xi+updateXi !update xi
            DO xiIdx=1,2
              !Boundary collision check
              IF(newXi(xiIdx)<0.0_DP) THEN 
                newXi(xiIdx)=0.0_DP
                newXi(3-xiIdx)=xi(3-xiIdx)-updateXi(3-xiIdx)*xi(xiIdx)/updateXi(xiIdx)
              ELSEIF(newXi(xiIdx)>1.0_DP) THEN
                newXi(xiIdx)=1.0_DP  
                newXi(3-xiIdx)=xi(3-xiIdx)+updateXi(3-xiIdx)*(1.0_DP-xi(xiIdx))/updateXi(xiIdx)
              ENDIF
            ENDDO
            CALL Field_InterpolateXi(SECOND_PART_DERIV,newXi,interpolatedPoint,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
            distanceVector=interpolatedPoint%values(1:numberOfCoordinates,NO_PART_DERIV)-dataPointLocation
            newFunctionValue=DOT_PRODUCT(distanceVector(1:numberOfCoordinates),distanceVector(1:numberOfCoordinates))
            !Second half of the convergence test (before collision detection)
            converged=converged.AND.(DABS(newFunctionValue-functionValue)/(1.0_DP+functionValue)<relativeTolerance) 
            IF(converged) THEN
              functionValue=newFunctionValue
              EXIT !converged: exit inner loop first
            ENDIF
            IF((newFunctionValue-functionValue)>absoluteTolerance) THEN
              !Bad model: reduce step size
              IF(delta<=minimumDelta) THEN
                !Something went wrong, minimumDelta too large? not likely to happen if minimumDelta is small
                exitTag=DATA_PROJECTION_EXIT_TAG_MAX_ITERATION !it will get stucked!!
                functionValue=newFunctionValue
                EXIT main_loop
              ENDIF
              delta=DMAX1(minimumDelta,0.25_DP*delta)
            ELSE
              predictedReduction=DOT_PRODUCT(functionGradient,updateXi)+ &
                & 0.5_DP*(updateXi(1)*(updateXi(1)*functionHessian(1,1)+2.0_DP*updateXi(2)*functionHessian(1,2))+ &
                & updateXi(2)**2*functionHessian(2,2))
              predictionAccuracy=(newFunctionValue-functionValue)/predictedReduction
              IF(predictionAccuracy<0.01_DP) THEN
                !Bad model: reduce region size
                IF(delta<=minimumDelta) THEN
                  !Something went wrong, minimumDelta too large? not likely to happen if minimumDelta is small
                  exitTag=DATA_PROJECTION_EXIT_TAG_MAX_ITERATION !it will get stucked!!
                  functionValue=newFunctionValue
                  EXIT main_loop
                ENDIF
                delta=DMAX1(minimumDelta,0.5_DP*delta)
              ELSEIF(predictionAccuracy>0.9_DP.AND.predictionAccuracy<1.1_DP) THEN
                !Good model: increase region size
                delta=DMIN1(maximumDelta,2.0_DP*delta)
                functionValue=newFunctionValue
                EXIT
              ELSE
                !OK model: keep the current region size
                functionValue=newFunctionValue
                EXIT
              ENDIF
            ENDIF
          ENDDO !iterationIdx2 (inner loop: adjust region size)
          functionValue=newFunctionValue
          xi=newXi
          IF(converged) THEN
            exitTag=DATA_PROJECTION_EXIT_TAG_CONVERGED
            EXIT
          ENDIF
        ENDDO main_loop !iterationIdx1 (outer loop)
        IF(exitTag==DATA_PROJECTION_EXIT_TAG_NO_ELEMENT.AND.iterationIdx1>=dataProjection%maximumNumberOfIterations) &
          & exitTag=DATA_PROJECTION_EXIT_TAG_MAX_ITERATION
        IF((projectionExitTag==DATA_PROJECTION_EXIT_TAG_NO_ELEMENT).OR.(DSQRT(functionValue)<projectionDistance)) THEN
          projectionExitTag=exitTag
          projectionElementNumber=elementNumber
          projectionElementFaceNumber=elementFaceNumber            
          projectionDistance=DSQRT(functionValue)
          projectionXi=xi
          projectionVector=distanceVector
        ENDIF
      ENDIF
    ENDDO !elementIdx
    
    EXITS("DataProjection_NewtonFacesEvaluate")
    RETURN
999 ERRORSEXITS("DataProjection_NewtonFacesEvaluate",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_NewtonFacesEvaluate

  !
  !================================================================================================================================
  !
  
  !>Find the projection of a data point onto element lines (slight difference to DataProjection_NewtonElementsEvaluate_1)
  SUBROUTINE DataProjection_NewtonLinesEvaluate(dataProjection,interpolatedPoint,dataPointLocation,candidateElements, &
    & candidateElementLines,projectionExitTag,projectionElementNumber,projectionElementLineNumber,projectionDistance, &
    & projectionXi,projectionVector,err,error,*)
    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<Data projection problem to evaluate
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: interpolatedPoint !<The interpolated point for the data point
    REAL(DP), INTENT(IN) :: dataPointLocation(:) !<dataPointLocation(componentIdx). The location of data point in the region.
    INTEGER(INTG), INTENT(IN) :: candidateElements(:) !<candidateElememnts(candidateIdx). The list of candidate elements for the projection.
    INTEGER(INTG), INTENT(IN) :: candidateElementLines(:) !<candidateElementLines(candidateIdx). The list of candidate lines for the projection.
    INTEGER(INTG), INTENT(OUT) :: projectionExitTag !<On exit, the exit tag status for the data point projection
    INTEGER(INTG), INTENT(OUT) :: projectionElementNumber !<On exit, the element number of the data point projection
    INTEGER(INTG), INTENT(OUT) :: projectionElementLineNumber !<On exit, the line number of the data point projection
    REAL(DP), INTENT(OUT) :: projectionDistance !<On exit, the projection distance
    REAL(DP), INTENT(OUT) :: projectionXi(1) !<projectionXi(xiIdx). On exit, the xi location of the data point projection
    REAL(DP), INTENT(OUT) :: projectionVector(3) !<projectionVector(componentIdx). On exit, the projection vector for the data point projection
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string   
    !Local Variables
    LOGICAL :: insideRegion,converged
    INTEGER(INTG) :: elementNumber,elementLineNumber,lineNumber !local element number in current computational domain
    INTEGER(INTG) :: numberOfCoordinates
    INTEGER(INTG) :: bound,exitTag
    REAL(DP) :: xi(1),newXi(1),updateXi(1),updateXiNorm !<xi
    REAL(DP) :: distanceVector(3),relativeTolerance,absoluteTolerance !<tolerances
    REAL(DP) :: functionValue,newFunctionValue
    REAL(DP) :: functionGradient,functionHessian
    REAL(DP) :: maximumDelta,minimumDelta,delta !<trust region size
    REAL(DP) :: predictedReduction,predictionAccuracy    
    INTEGER(INTG) :: elementIdx,iterationIdx1,iterationIdx2
    
    ENTERS("DataProjection_NewtonLinesEvaluate",err,error,*999)
              
    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionFinished) CALL FlagError("Data projection has not been finished.",err,error,*999)
    
    projectionExitTag=DATA_PROJECTION_EXIT_TAG_NO_ELEMENT
    numberOfCoordinates=dataProjection%numberOfCoordinates
    relativeTolerance=dataProjection%relativeTolerance
    absoluteTolerance=dataProjection%absoluteTolerance
    maximumDelta=dataProjection%maximumIterationUpdate
    minimumDelta=0.025_DP*maximumDelta !need to set a minimum, in case if it gets too small
    !Project on each candidate elements
    DO elementIdx=1,SIZE(candidateElements,1) 
      elementNumber=candidateElements(elementIdx)
      IF(elementNumber>0) THEN
        elementLineNumber=candidateElementLines(elementIdx)
        lineNumber=interpolatedPoint%INTERPOLATION_PARAMETERS%field%decomposition%topology%elements% &
          & elements(elementNumber)%ELEMENT_LINES(elementLineNumber)
        exitTag=DATA_PROJECTION_EXIT_TAG_NO_ELEMENT
        converged=.FALSE.
        !start at half the maximumDelta as we do not know if quadratic model is a good approximation yet
        delta=0.5_DP*maximumDelta
        CALL Field_InterpolationParametersLineGet(dataProjection%projectionSetType,lineNumber,interpolatedPoint% &
          & INTERPOLATION_PARAMETERS,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
        xi=dataProjection%startingXi
        CALL Field_InterpolateXi(SECOND_PART_DERIV,xi,interpolatedPoint,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
        distanceVector(1:numberOfCoordinates)=interpolatedPoint%values(1:numberOfCoordinates,NO_PART_DERIV)- &
          & dataPointLocation(1:numberOfCoordinates)
        functionValue=DOT_PRODUCT(distanceVector(1:numberOfCoordinates),distanceVector(1:numberOfCoordinates))
        !Outer loop
        main_loop: DO iterationIdx1=1,dataProjection%maximumNumberOfIterations 
          !Check for bounds [0,1]
          IF(ABS(xi(1))<ZERO_TOLERANCE) THEN
            bound=-1 !bound at negative direction             
          ELSEIF(ABS(xi(1)-1.0_DP)<ZERO_TOLERANCE) THEN
            bound=1 !bound at positive direction
          ELSE 
            bound=0 !inside the bounds
          ENDIF
          !functionGradient 
          functionGradient=2.0_DP*(DOT_PRODUCT(distanceVector(1:numberOfCoordinates), &
            & interpolatedPoint%values(1:numberOfCoordinates,FIRST_PART_DERIV)))
          !functionHessian 
          functionHessian=2.0_DP*(DOT_PRODUCT(distanceVector(1:numberOfCoordinates), &
            & interpolatedPoint%values(1:numberOfCoordinates,SECOND_PART_DERIV))- &
            & DOT_PRODUCT(interpolatedPoint%values(1:numberOfCoordinates,FIRST_PART_DERIV), &
            & interpolatedPoint%values(1:numberOfCoordinates,FIRST_PART_DERIV)))
          !A model trust region approach, directly taken from CMISS CLOS22: V = -(H + eigenShift*I)g
          !The calculation of eigenShift are only approximated as oppose to the common trust region approach
          !Inner loop: adjust region size, usually EXIT at 1 or 2 iterations
          DO iterationIdx2=1,dataProjection%maximumNumberOfIterations 
            insideRegion=.FALSE.
            IF(functionHessian>absoluteTolerance) THEN
              !Positive: minimum exists
              updateXi(1)=-functionGradient/functionHessian
              updateXiNorm=DABS(updateXi(1))
              insideRegion=updateXiNorm<=delta
            ENDIF
            IF(.NOT.insideRegion) THEN
              !minimum not in the region
              updateXi(1)=-DSIGN(delta,functionGradient)
              updateXiNorm=delta
            ENDIF
            IF((bound/=0).AND.(bound>0.EQV.updateXi(1)>0.0_DP)) THEN
              !Projection will go out of element bounds
              exitTag=DATA_PROJECTION_EXIT_TAG_BOUNDS
              EXIT main_loop
            ENDIF
            !First half of the convergence test (before collision detection)
            converged=updateXiNorm<absoluteTolerance 
            newXi=xi+updateXi !update xi
            !Boundary collision check
            IF(newXi(1)<0.0_DP) THEN 
              newXi(1)=0.0_DP
            ELSEIF(newXi(1)>1.0_DP) THEN
              newXi(1)=1.0_DP  
            ENDIF
            CALL Field_InterpolateXi(SECOND_PART_DERIV,newXi,interpolatedPoint,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
            distanceVector(1:numberOfCoordinates)=interpolatedPoint%values(1:numberOfCoordinates,NO_PART_DERIV)- &
              & dataPointLocation(1:numberOfCoordinates)
            newFunctionValue=DOT_PRODUCT(distanceVector(1:numberOfCoordinates),distanceVector(1:numberOfCoordinates))
            !Second half of the convergence test
            converged=converged.AND.(DABS(newFunctionValue-functionValue)/(1.0_DP+functionValue)<relativeTolerance) 
            IF(converged) EXIT !converged: exit inner loop first
            IF((newFunctionValue-functionValue)>absoluteTolerance) THEN
              !Bad model: reduce step size
              IF(delta<=minimumDelta) THEN
                !Something went wrong, minimumDelta too large? not likely to happen if minimumDelta is small
                exitTag=DATA_PROJECTION_EXIT_TAG_MAX_ITERATION !it will get stucked!!
                EXIT main_loop
              ENDIF
              delta=DMAX1(minimumDelta,0.25_DP*delta)
            ELSE
              predictedReduction=updateXi(1)*(functionGradient+0.5_DP*functionHessian*updateXi(1))
              predictionAccuracy=(newFunctionValue-functionValue)/predictedReduction
              IF(predictionAccuracy<0.01_DP) THEN
                !Bad model: reduce region size
                IF(delta<=minimumDelta) THEN
                  !Something went wrong, minimumDelta too large? not likely to happen if minimumDelta is small
                  exitTag=DATA_PROJECTION_EXIT_TAG_MAX_ITERATION !it will get stucked!!
                  EXIT main_loop
                ENDIF
                delta=DMAX1(minimumDelta,0.5_DP*delta)
              ELSEIF(predictionAccuracy>0.9_DP.AND.predictionAccuracy<1.1_DP) THEN
                !Good model: increase region size
                delta=DMIN1(maximumDelta,2.0_DP*delta)
                EXIT
              ELSE
                !OK model: keep the current region size
                EXIT
              ENDIF
            ENDIF
          ENDDO !iterationIdx2 (inner loop: adjust region size)
          functionValue=newFunctionValue
          xi=newXi
          IF(converged) THEN
            exitTag=DATA_PROJECTION_EXIT_TAG_CONVERGED
            EXIT
          ENDIF
        ENDDO main_loop !iterationIdx1 (outer loop)
        IF(exitTag==DATA_PROJECTION_EXIT_TAG_NO_ELEMENT.AND.iterationIdx1>=dataProjection%maximumNumberOfIterations) &
          & exitTag=DATA_PROJECTION_EXIT_TAG_MAX_ITERATION
        IF((projectionExitTag==DATA_PROJECTION_EXIT_TAG_NO_ELEMENT).OR.(DSQRT(functionValue)<projectionDistance)) THEN
          projectionExitTag=exitTag
          projectionElementNumber=elementNumber
          projectionElementLineNumber=elementLineNumber
          projectionDistance=DSQRT(functionValue)
          projectionXi=xi
          projectionVector=distanceVector
        ENDIF
      ENDIF
    ENDDO !elementIdx
    
    EXITS("DataProjection_NewtonLinesEvaluate")
    RETURN
999 ERRORSEXITS("DataProjection_NewtonLinesEvaluate",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_NewtonLinesEvaluate
  
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

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionFinished) CALL FlagError("Data projection has not been finished.",err,error,*999)
    
    numberOfClosestElements=dataProjection%numberOfClosestElements       
    
    EXITS("DataProjection_NumberOfClosestElementsGet")
    RETURN
999 ERRORSEXITS("DataProjection_NumberOfClosestElementsGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_NumberOfClosestElementsGet

  !
  !================================================================================================================================
  !
  
  !>Sets the number of closest elements for a data projection.
  SUBROUTINE DataProjection_NumberOfClosestElementsSet(dataProjection,numberOfClosestElements,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to set the number of closest elements for
    INTEGER(INTG), INTENT(IN) :: numberOfClosestElements !<the number of closest elements to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("DataProjection_NumberOfClosestElementsSet",err,error,*999)

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(dataProjection%dataProjectionFinished) CALL FlagError("Data projection has been finished.",err,error,*999)
    IF(numberOfClosestElements<1) CALL FlagError("Data projection number of closest elements must be >= 1.",err,error,*999)
        
    dataProjection%numberOfClosestElements=numberOfClosestElements
    
    EXITS("DataProjection_NumberOfClosestElementsSet")
    RETURN
999 ERRORSEXITS("DataProjection_NumberOfClosestElementsSet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_NumberOfClosestElementsSet
  
  !
  !================================================================================================================================
  !
  
  !>Cancel the data projection for data points based on a data point user number.
  SUBROUTINE DataProjection_ProjectionCancelByDataPoints0(dataProjection,dataPointUserNumber,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to cancel projections for
    INTEGER(INTG), INTENT(IN) :: dataPointUserNumber !<The data point user number to use to cancel projections 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("DataProjection_ProjectionCancelByDataPoints0",err,error,*999)

    CALL DataProjection_ProjectionCancelByDataPoints1(dataProjection,[dataPointUserNumber],err,error,*999)
    
    EXITS("DataProjection_ProjectionCancelByDataPoints0")
    RETURN
999 ERRORS("DataProjection_ProjectionCancelByDataPoints0",err,error)    
    EXITS("DataProjection_ProjectionCancelByDataPoints0")    
    RETURN 1

  END SUBROUTINE DataProjection_ProjectionCancelByDataPoints0
  
  !
  !================================================================================================================================
  !
  
  !>Cancel the data projection for data points based on data point user numbers.
  SUBROUTINE DataProjection_ProjectionCancelByDataPoints1(dataProjection,dataPointUserNumbers,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to cancel projections for
    INTEGER(INTG), INTENT(IN) :: dataPointUserNumbers(:) !<dataPointUserNumbers(dataPointIdx). The data point user numbers to use to cancel projections.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointIdx,dataPointGlobalNumber
   
    ENTERS("DataProjection_ProjectionCancelByDataPoints1",err,error,*999)

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionFinished) CALL FlagError("Data projection has not been finished.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionProjected) CALL FlagError("Data projection has not been evaluated.",err,error,*999)
    IF(.NOT.ALLOCATED(dataProjection%dataProjectionResults)) &
       CALL FlagError("Data projection projection results is not allocated.",err,error,*999)

    DO dataPointIdx=1,SIZE(dataPointUserNumbers,1)
      CALL DataProjection_DataPointGlobalNumberGet(dataProjection,dataPointUserNumbers(dataPointIdx),dataPointGlobalNumber, &
        & err,error,*999)
      dataProjection%dataProjectionResults(dataPointGlobalNumber)%exitTag=DATA_PROJECTION_CANCELLED
      dataProjection%dataProjectionResults(dataPointIdx)%elementNumber=0
      dataProjection%dataProjectionResults(dataPointIdx)%elementLineFaceNumber=0
      dataProjection%dataProjectionResults(dataPointIdx)%xi=0.0_DP
      dataProjection%dataProjectionResults(dataPointIdx)%elementXi=0.0_DP
      dataProjection%dataProjectionResults(dataPointIdx)%distance=0.0_DP
      dataProjection%dataProjectionResults(dataPointIdx)%projectionVector=0.0_DP          
    ENDDO !dataPointIdx
    
    !Re-evaluate errors
    CALL DataProjection_ErrorsCalculate(dataProjection,err,error,*999)
    
    EXITS("DataProjection_ProjectionCancelByDataPoints1")
    RETURN
999 ERRORS("DataProjection_ProjectionCancelByDataPoints1",err,error)    
    EXITS("DataProjection_ProjectionCancelByDataPoints1")    
    RETURN 1

  END SUBROUTINE DataProjection_ProjectionCancelByDataPoints1
  
  !
  !================================================================================================================================
  !
  
  !>Cancel the data projection for data points based on the projection distance.
  SUBROUTINE DataProjection_ProjectionCancelByDistance(dataProjection,distanceRelation,distance,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to cancel projections for
    INTEGER(INTG), INTENT(IN) :: distanceRelation !<The distance relation to use to cancel projections \see OPENCMISS_DataProjectionDistanceRelations
    REAL(DP), INTENT(IN) :: distance !<The distance by which to select the data points to cancel.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointIdx
    TYPE(DataPointsType), POINTER :: dataPoints
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DataProjection_ProjectionCancelByDistance",err,error,*999)

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionFinished) CALL FlagError("Data projection has not been finished.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionProjected) CALL FlagError("Data projection has not been evaluated.",err,error,*999)
    IF(.NOT.ALLOCATED(dataProjection%dataProjectionResults)) &
       CALL FlagError("Data projection projection results is not allocated.",err,error,*999)

    NULLIFY(dataPoints)
    CALL DataProjection_DataPointsGet(dataProjection,dataPoints,err,error,*999)
    SELECT CASE(distanceRelation)
    CASE(DATA_PROJECTION_DISTANCE_GREATER)
      DO dataPointIdx=1,dataPoints%numberOfDataPoints
        IF(dataProjection%dataProjectionResults(dataPointIdx)%distance>distance) THEN
          dataProjection%dataProjectionResults(dataPointIdx)%exitTag=DATA_PROJECTION_CANCELLED
          dataProjection%dataProjectionResults(dataPointIdx)%elementNumber=0
          dataProjection%dataProjectionResults(dataPointIdx)%elementLineFaceNumber=0
          dataProjection%dataProjectionResults(dataPointIdx)%xi=0.0_DP
          dataProjection%dataProjectionResults(dataPointIdx)%elementXi=0.0_DP
          dataProjection%dataProjectionResults(dataPointIdx)%distance=0.0_DP
          dataProjection%dataProjectionResults(dataPointIdx)%projectionVector=0.0_DP          
        ENDIF
      ENDDO !dataPointIdx
    CASE(DATA_PROJECTION_DISTANCE_GREATER_EQUAL)
      DO dataPointIdx=1,dataPoints%numberOfDataPoints
        IF(dataProjection%dataProjectionResults(dataPointIdx)%distance>=distance) THEN
          dataProjection%dataProjectionResults(dataPointIdx)%exitTag=DATA_PROJECTION_CANCELLED
          dataProjection%dataProjectionResults(dataPointIdx)%elementNumber=0
          dataProjection%dataProjectionResults(dataPointIdx)%elementLineFaceNumber=0
          dataProjection%dataProjectionResults(dataPointIdx)%xi=0.0_DP
          dataProjection%dataProjectionResults(dataPointIdx)%elementXi=0.0_DP
          dataProjection%dataProjectionResults(dataPointIdx)%distance=0.0_DP
          dataProjection%dataProjectionResults(dataPointIdx)%projectionVector=0.0_DP          
        ENDIF
      ENDDO !dataPointIdx
    CASE(DATA_PROJECTION_DISTANCE_LESS)
      DO dataPointIdx=1,dataPoints%numberOfDataPoints
        IF(dataProjection%dataProjectionResults(dataPointIdx)%distance<distance) THEN
          dataProjection%dataProjectionResults(dataPointIdx)%exitTag=DATA_PROJECTION_CANCELLED
          dataProjection%dataProjectionResults(dataPointIdx)%elementNumber=0
          dataProjection%dataProjectionResults(dataPointIdx)%elementLineFaceNumber=0
          dataProjection%dataProjectionResults(dataPointIdx)%xi=0.0_DP
          dataProjection%dataProjectionResults(dataPointIdx)%elementXi=0.0_DP
          dataProjection%dataProjectionResults(dataPointIdx)%distance=0.0_DP
          dataProjection%dataProjectionResults(dataPointIdx)%projectionVector=0.0_DP          
        ENDIF
      ENDDO !dataPointIdx
    CASE(DATA_PROJECTION_DISTANCE_LESS_EQUAL)
      DO dataPointIdx=1,dataPoints%numberOfDataPoints
        IF(dataProjection%dataProjectionResults(dataPointIdx)%distance<=distance) THEN
          dataProjection%dataProjectionResults(dataPointIdx)%exitTag=DATA_PROJECTION_CANCELLED
          dataProjection%dataProjectionResults(dataPointIdx)%elementNumber=0
          dataProjection%dataProjectionResults(dataPointIdx)%elementLineFaceNumber=0
          dataProjection%dataProjectionResults(dataPointIdx)%xi=0.0_DP
          dataProjection%dataProjectionResults(dataPointIdx)%elementXi=0.0_DP
          dataProjection%dataProjectionResults(dataPointIdx)%distance=0.0_DP
          dataProjection%dataProjectionResults(dataPointIdx)%projectionVector=0.0_DP          
        ENDIF
      ENDDO !dataPointIdx
    CASE DEFAULT
      localError="The specified distance relation of "//TRIM(NumberToVString(distanceRelation,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    !Re-evaluate errors
    CALL DataProjection_ErrorsCalculate(dataProjection,err,error,*999)
    
    EXITS("DataProjection_ProjectionCancelByDistance")
    RETURN
999 ERRORSEXITS("DataProjection_ProjectionCancelByDistance",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_ProjectionCancelByDistance
  
  !
  !================================================================================================================================
  !
  
  !>Cancel the data projection for data points based on the projection exit tag.
  SUBROUTINE DataProjection_ProjectionCancelByExitTags0(dataProjection,exitTag,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to cancel projections for
    INTEGER(INTG), INTENT(IN) :: exitTag !<The exit tag to use to cancel projections \see DataProjectionExitTags
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("DataProjection_ProjectionCancelByExitTags0",err,error,*999)

    CALL DataProjection_ProjectionCancelByExitTags1(dataProjection,[exitTag],err,error,*999)
    
    EXITS("DataProjection_ProjectionCancelByExitTags0")
    RETURN
999 ERRORSEXITS("DataProjection_ProjectionCancelByExitTags0",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_ProjectionCancelByExitTags0
  
  !
  !================================================================================================================================
  !
  
  !>Cancel the data projection for data points based on the projection exit tags.
  SUBROUTINE DataProjection_ProjectionCancelByExitTags1(dataProjection,exitTags,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to cancel projections for
    INTEGER(INTG), INTENT(IN) :: exitTags(:) !<exitTags(tagIdx). The exit tags to use to cancel projections \see DataProjectionExitTags
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointIdx,tagIdx
    TYPE(DataPointsType), POINTER :: dataPoints
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DataProjection_ProjectionCancelByExitTags1",err,error,*999)

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionFinished) CALL FlagError("Data projection has not been finished.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionProjected) CALL FlagError("Data projection has not been evaluated.",err,error,*999)
    IF(.NOT.ALLOCATED(dataProjection%dataProjectionResults)) &
       CALL FlagError("Data projection projection results is not allocated.",err,error,*999)
    DO tagIdx=1,SIZE(exitTags,1)
      SELECT CASE(exitTags(tagIdx))
      CASE(DATA_PROJECTION_CANCELLED)
        !Do nothing or should this be an error?
      CASE(DATA_PROJECTION_EXIT_TAG_CONVERGED,DATA_PROJECTION_EXIT_TAG_BOUNDS,DATA_PROJECTION_EXIT_TAG_MAX_ITERATION, &
        & DATA_PROJECTION_EXIT_TAG_NO_ELEMENT)
        !Do nothing
      CASE DEFAULT
        localError="The specified exit tag of "//TRIM(NumberToVString(exitTags(tagIdx),"*",err,error))// &
          & " at position "//TRIM(NumberToVString(tagIdx,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDDO !tagIdx
    
    NULLIFY(dataPoints)
    CALL DataProjection_DataPointsGet(dataProjection,dataPoints,err,error,*999)
    DO dataPointIdx=1,dataPoints%numberOfDataPoints
      DO tagIdx=1,SIZE(exitTags,1)
        IF(dataProjection%dataProjectionResults(dataPointIdx)%exitTag==exitTags(tagIdx)) THEN
          dataProjection%dataProjectionResults(dataPointIdx)%exitTag=DATA_PROJECTION_CANCELLED
          dataProjection%dataProjectionResults(dataPointIdx)%elementNumber=0
          dataProjection%dataProjectionResults(dataPointIdx)%elementLineFaceNumber=0
          dataProjection%dataProjectionResults(dataPointIdx)%xi=0.0_DP
          dataProjection%dataProjectionResults(dataPointIdx)%elementXi=0.0_DP
          dataProjection%dataProjectionResults(dataPointIdx)%distance=0.0_DP
          dataProjection%dataProjectionResults(dataPointIdx)%projectionVector=0.0_DP          
        ENDIF
      ENDDO !tagIdx
    ENDDO !dataPointIdx
    
    !Re-evaluate errors
    CALL DataProjection_ErrorsCalculate(dataProjection,err,error,*999)
    
    EXITS("DataProjection_ProjectionCancelByExitTags1")
    RETURN
999 ERRORSEXITS("DataProjection_ProjectionCancelByExitTags1",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_ProjectionCancelByExitTags1
  
  !
  !================================================================================================================================
  !
  
  !>Sets the candidates element numbers numbers for a data projection.
  SUBROUTINE DataProjection_ProjectionCandidateElementsSet(dataProjection,elementUserNumbers,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to set the projection type for
    INTEGER(INTG), INTENT(IN) :: elementUserNumbers(:) !<elementUserNumbers(elementIdx). The projection candidate user element numbers
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: elementIdx,elementLocalNumber,numberOfCandidates
    LOGICAL :: elementExists,ghostElement
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: decompositionTopology
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DataProjection_ProjectionCandidateElementsSet",err,error,*998)

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*998)
    IF(.NOT.ALLOCATED(dataProjection%dataProjectionCandidates)) &
      & CALL FlagError("Data projection candidates is not allocated.",err,error,*998)
    
    SELECT CASE(dataProjection%projectionType)
    CASE(DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE)
      CALL FlagError("Use the candidate lines routine not the candidate elements routine for a boundary lines projection.", &
        & err,error,*998)
    CASE(DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE)
      CALL FlagError("Use the candidate faces routine not the candidate elements routine for a boundary faces projection.", &
        & err,error,*998)
    CASE(DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE)
      NULLIFY(decomposition)
      CALL DataProjection_DecompositionGet(dataProjection,decomposition,err,error,*998)
      NULLIFY(decompositionTopology)
      CALL Decomposition_TopologyGet(decomposition,decompositionTopology,err,error,*998)
      IF(ALLOCATED(dataProjection%dataProjectionCandidates(0)%candidateElementNumbers)) &
        & DEALLOCATE(dataProjection%dataProjectionCandidates(0)%candidateElementNumbers)
      ALLOCATE(dataProjection%dataProjectionCandidates(0)%candidateElementNumbers(SIZE(elementUserNumbers,1)),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate candidiate element numbers.",err,error,*999)
      numberOfCandidates=0
      DO elementIdx=1,SIZE(elementUserNumbers,1)
        CALL DecompositionTopology_ElementCheckExists(decompositionTopology,elementUserNumbers(elementIdx), &
          & elementExists,elementLocalNumber,ghostElement,err,error,*999)       
        IF(elementExists) THEN
          IF(.NOT.ghostElement) THEN
            numberOfCandidates=numberOfCandidates+1
            dataProjection%dataProjectionCandidates(0)%candidateElementNumbers(numberOfCandidates)=elementLocalNumber
          ENDIF
        ELSE
          localError="Element user number "//TRIM(NumberToVString(elementUserNumbers(elementIdx),"*",err,error))// &
            & " does not exist."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDDO !elementIdx
      IF(numberOfCandidates>dataProjection%maxNumberOfCandidates) dataProjection%maxNumberOfCandidates=numberOfCandidates
    CASE DEFAULT
      localError="The data projection type of "//TRIM(NumberToVString(dataProjection%projectionType,"",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*998)
    END SELECT
    
    EXITS("DataProjection_ProjectionCandidateElementsSet")
    RETURN
999 IF(ALLOCATED(dataProjection%dataProjectionCandidates(0)%candidateElementNumbers)) &
      & DEALLOCATE(dataProjection%dataProjectionCandidates(0)%candidateElementNumbers)
998 ERRORS("DataProjection_ProjectionCandidateElementsSet",err,error)
    EXITS("DataProjection_ProjectionCandidateElementsSet")
    RETURN 1
    
  END SUBROUTINE DataProjection_ProjectionCandidateElementsSet
  
  !
  !================================================================================================================================
  !
  
  !>Sets the candidates element numbers numbers for data points in a data projection.
  SUBROUTINE DataProjection_ProjectionDataCandidateElementsSet(dataProjection,dataPointUserNumbers,elementUserNumbers,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to set the projection type for
    INTEGER(INTG), INTENT(IN) :: dataPointUserNumbers(:) !<dataPointUserNumbers(dataPointIdx). The data points for the projection candidate user element numbers
    INTEGER(INTG), INTENT(IN) :: elementUserNumbers(:) !<elementUserNumbers(elementIdx). The projection candidate user element numbers
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointGlobalNumber,dataPointIdx,elementIdx,elementLocalNumber,numberOfCandidates
    LOGICAL :: elementExists,ghostElement
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: decompositionTopology
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DataProjection_ProjectionDataCandidateElementsSet",err,error,*999)

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(dataProjection%dataProjectionCandidates)) &
      & CALL FlagError("Data projection candidates is not allocated.",err,error,*999)
    
    SELECT CASE(dataProjection%projectionType)
    CASE(DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE)
      CALL FlagError("Use the data candidate lines routine not the candidate elements routine for a boundary lines projection.", &
        & err,error,*999)
    CASE(DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE)
      CALL FlagError("Use the data candidate faces routine not the candidate elements routine for a boundary faces projection.", &
        & err,error,*999)
    CASE(DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE)
      NULLIFY(decomposition)
      CALL DataProjection_DecompositionGet(dataProjection,decomposition,err,error,*999)
      NULLIFY(decompositionTopology)
      CALL Decomposition_TopologyGet(decomposition,decompositionTopology,err,error,*999)
      DO dataPointIdx=1,SIZE(dataPointUserNumbers,1)
        CALL DataProjection_DataPointGlobalNumberGet(dataProjection,dataPointUserNumbers(dataPointIdx),dataPointGlobalNumber, &
          & err,error,*999)
        IF(ALLOCATED(dataProjection%dataProjectionCandidates(dataPointGlobalNumber)%candidateElementNumbers)) &
          & DEALLOCATE(dataProjection%dataProjectionCandidates(dataPointGlobalNumber)%candidateElementNumbers)
        ALLOCATE(dataProjection%dataProjectionCandidates(dataPointGlobalNumber)%candidateElementNumbers( &
          & SIZE(elementUserNumbers,1)),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate candidiate element numbers.",err,error,*999)
        numberOfCandidates=0
        DO elementIdx=1,SIZE(elementUserNumbers,1)
          CALL DecompositionTopology_ElementCheckExists(decompositionTopology,elementUserNumbers(elementIdx), &
            & elementExists,elementLocalNumber,ghostElement,err,error,*999)       
          IF(elementExists) THEN
            IF(.NOT.ghostElement) THEN
              numberOfCandidates=numberOfCandidates+1
              dataProjection%dataProjectionCandidates(dataPointGlobalNumber)%candidateElementNumbers(numberOfCandidates)= &
                & elementLocalNumber
            ENDIF
          ELSE
            localError="Element user number "//TRIM(NumberToVString(elementUserNumbers(elementIdx),"*",err,error))// &
              & " does not exist."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ENDDO !elementIdx
        IF(numberOfCandidates>dataProjection%maxNumberOfCandidates) dataProjection%maxNumberOfCandidates=numberOfCandidates
      ENDDO !dataPointIdx
    CASE DEFAULT
      localError="The data projection type of "//TRIM(NumberToVString(dataProjection%projectionType,"",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("DataProjection_ProjectionDataCandidateElementsSet")
    RETURN
999 ERRORS("DataProjection_ProjectionDataCandidateElementsSet",err,error)
    EXITS("DataProjection_ProjectionDataCandidateElementsSet")
    RETURN 1
    
  END SUBROUTINE DataProjection_ProjectionDataCandidateElementsSet
  
  !
  !================================================================================================================================
  !
  
  !>Sets the candidates element numbers and local faces for a data projection.
  SUBROUTINE DataProjection_ProjectionCandidateFacesSet(dataProjection,elementUserNumbers,localFaceNormals,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to set the projection type for
    INTEGER(INTG), INTENT(IN) :: elementUserNumbers(:) !<elementUserNumbers(elementIdx). The projection candidate user element numbers
    INTEGER(INTG), INTENT(IN) :: localFaceNormals(:) !<localFaceNormals(elementIdx). The projection candidate element face normals
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: elementIdx,elementLocalNumber,localFaceNumber,numberOfCandidates
    LOGICAL :: elementExists,ghostElement
    TYPE(BASIS_TYPE), POINTER :: basis
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: decompositionTopology
    TYPE(DOMAIN_TYPE), POINTER :: domain
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: domainElements
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: domainTopology
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DataProjection_ProjectionCandidateFacesSet",err,error,*998)

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*998)
    IF(.NOT.ALLOCATED(dataProjection%dataProjectionCandidates)) &
      & CALL FlagError("Data projection candidates is not allocated.",err,error,*998)
    
    SELECT CASE(dataProjection%projectionType)
    CASE(DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE)
      CALL FlagError("Use the candidate lines routine not the candidate faces routine for a boundary lines projection.", &
        & err,error,*998)
    CASE(DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE)
      CALL FlagError("Use the candidate elements routine not the candidate faces routine for an all elements projection.", &
        & err,error,*998)
    CASE(DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE)
      IF(SIZE(elementUserNumbers,1)/=SIZE(localFaceNormals,1)) THEN
        localError="The size of the element user numbers array of "// &
          & TRIM(NumberToVString(SIZE(elementUserNumbers,1),"*",err,error))// &
          & " does not match the size of the face normals array of "// &
          & TRIM(NumberToVString(SIZE(localFaceNormals,1),"*",err,error))//"."
        CALL FlagError(localError,err,error,*998)
      ENDIF
      NULLIFY(decomposition)
      CALL DataProjection_DecompositionGet(dataProjection,decomposition,err,error,*998)
      NULLIFY(decompositionTopology)
      CALL Decomposition_TopologyGet(decomposition,decompositionTopology,err,error,*998)
      NULLIFY(domain)
      CALL Decomposition_DomainGet(decomposition,0,domain,err,error,*998)
      NULLIFY(domainTopology)
      CALL Domain_TopologyGet(domain,domainTopology,err,error,*998)
      NULLIFY(domainElements)
      CALL DomainTopology_ElementsGet(domainTopology,domainElements,err,error,*998)
      IF(ALLOCATED(dataProjection%dataProjectionCandidates(0)%candidateElementNumbers)) &
        & DEALLOCATE(dataProjection%dataProjectionCandidates(0)%candidateElementNumbers)
      IF(ALLOCATED(dataProjection%dataProjectionCandidates(0)%localFaceLineNumbers)) &
        & DEALLOCATE(dataProjection%dataProjectionCandidates(0)%localFaceLineNumbers)
      ALLOCATE(dataProjection%dataProjectionCandidates(0)%candidateElementNumbers(SIZE(elementUserNumbers,1)),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate candidiate element numbers.",err,error,*999)
      ALLOCATE(dataProjection%dataProjectionCandidates(0)%localFaceLineNumbers(SIZE(localFaceNormals,1)),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate candidiate local face/line numbers.",err,error,*999)
      numberOfCandidates=0
      DO elementIdx=1,SIZE(elementUserNumbers,1)
        CALL DecompositionTopology_ElementCheckExists(decompositionTopology,elementUserNumbers(elementIdx), &
          & elementExists,elementLocalNumber,ghostElement,err,error,*999)       
        IF(elementExists) THEN
          IF(.NOT.ghostElement) THEN
            numberOfCandidates=numberOfCandidates+1           
            NULLIFY(basis)
            CALL DomainElements_BasisGet(domainElements,elementLocalNumber,basis,err,error,*999)
            CALL Basis_LocalFaceNumberGet(basis,localFaceNormals(elementIdx),localFaceNumber,err,error,*999)
            dataProjection%dataProjectionCandidates(0)%candidateElementNumbers(numberOfCandidates)=elementLocalNumber
            dataProjection%dataProjectionCandidates(0)%localFaceLineNumbers(numberOfCandidates)=localFaceNumber
          ENDIF
        ELSE
          localError="Element user number "//TRIM(NumberToVString(elementUserNumbers(elementIdx),"*",err,error))// &
            & " does not exist."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDDO !elementIdx
      IF(numberOfCandidates>dataProjection%maxNumberOfCandidates) dataProjection%maxNumberOfCandidates=numberOfCandidates
    CASE DEFAULT
      localError="The data projection type of "//TRIM(NumberToVString(dataProjection%projectionType,"",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*998)
    END SELECT
     
    EXITS("DataProjection_ProjectionCandidateFacesSet")
    RETURN
999 IF(ALLOCATED(dataProjection%dataProjectionCandidates(0)%candidateElementNumbers)) &
      & DEALLOCATE(dataProjection%dataProjectionCandidates(0)%candidateElementNumbers)
    IF(ALLOCATED(dataProjection%dataProjectionCandidates(0)%localFaceLineNumbers)) &
      & DEALLOCATE(dataProjection%dataProjectionCandidates(0)%localFaceLineNumbers)
998 ERRORSEXITS("DataProjection_ProjectionCandidateFacesSet",err,error)
    RETURN 1

  END SUBROUTINE DataProjection_ProjectionCandidateFacesSet
  
  !
  !================================================================================================================================
  !
  
  !>Sets the candidates element numbers and local faces for data points in a data projection.
  SUBROUTINE DataProjection_ProjectionDataCandidateFacesSet(dataProjection,dataPointUserNumbers,elementUserNumbers, &
    & localFaceNormals,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to set the projection type for
    INTEGER(INTG), INTENT(IN) :: dataPointUserNumbers(:) !<dataPointUserNumbers(dataPointIdx). The data points for the projection candidate user element numbers
    INTEGER(INTG), INTENT(IN) :: elementUserNumbers(:) !<elementUserNumbers(elementIdx). The projection candidate user element numbers
    INTEGER(INTG), INTENT(IN) :: localFaceNormals(:) !<localFaceNormals(elementIdx). The projection candidate element face normals
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointGlobalNumber,dataPointIdx,elementIdx,elementLocalNumber,localFaceNumber,numberOfCandidates
    LOGICAL :: elementExists,ghostElement
    TYPE(BASIS_TYPE), POINTER :: basis
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: decompositionTopology
    TYPE(DOMAIN_TYPE), POINTER :: domain
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: domainElements
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: domainTopology
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DataProjection_ProjectionDataCandidateFacesSet",err,error,*999)

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(dataProjection%dataProjectionCandidates)) &
      & CALL FlagError("Data projection candidates is not allocated.",err,error,*999)
    
    SELECT CASE(dataProjection%projectionType)
    CASE(DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE)
      CALL FlagError("Use the data candidate lines routine not the candidate faces routine for a boundary lines projection.", &
        & err,error,*999)
    CASE(DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE)
      CALL FlagError("Use the data candidate elements routine not the candidate faces routine for an all elements projection.", &
        & err,error,*999)
    CASE(DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE)
      IF(SIZE(elementUserNumbers,1)/=SIZE(localFaceNormals,1)) THEN
        localError="The size of the element user numbers array of "// &
          & TRIM(NumberToVString(SIZE(elementUserNumbers,1),"*",err,error))// &
          & " does not match the size of the face normals array of "// &
          & TRIM(NumberToVString(SIZE(localFaceNormals,1),"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      NULLIFY(decomposition)
      CALL DataProjection_DecompositionGet(dataProjection,decomposition,err,error,*999)
      NULLIFY(decompositionTopology)
      CALL Decomposition_TopologyGet(decomposition,decompositionTopology,err,error,*999)
      NULLIFY(domain)
      CALL Decomposition_DomainGet(decomposition,0,domain,err,error,*999)
      NULLIFY(domainTopology)
      CALL Domain_TopologyGet(domain,domainTopology,err,error,*999)
      NULLIFY(domainElements)
      CALL DomainTopology_ElementsGet(domainTopology,domainElements,err,error,*999)
      DO dataPointIdx=1,SIZE(dataPointUserNumbers,1)
        CALL DataProjection_DataPointGlobalNumberGet(dataProjection,dataPointUserNumbers(dataPointIdx),dataPointGlobalNumber, &
          & err,error,*999)
        IF(ALLOCATED(dataProjection%dataProjectionCandidates(dataPointGlobalNumber)%candidateElementNumbers)) &
          & DEALLOCATE(dataProjection%dataProjectionCandidates(dataPointGlobalNumber)%candidateElementNumbers)
        IF(ALLOCATED(dataProjection%dataProjectionCandidates(dataPointGlobalNumber)%localFaceLineNumbers)) &
          & DEALLOCATE(dataProjection%dataProjectionCandidates(dataPointGlobalNumber)%localFaceLineNumbers)
        ALLOCATE(dataProjection%dataProjectionCandidates(dataPointGlobalNumber)%candidateElementNumbers( &
          & SIZE(elementUserNumbers,1)),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate candidiate element numbers.",err,error,*999)
        ALLOCATE(dataProjection%dataProjectionCandidates(dataPointGlobalNumber)%localFaceLineNumbers( &
          & SIZE(localFaceNormals,1)),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate candidiate local face/line numbers.",err,error,*999)
        numberOfCandidates=0
        DO elementIdx=1,SIZE(elementUserNumbers,1)
          CALL DecompositionTopology_ElementCheckExists(decompositionTopology,elementUserNumbers(elementIdx), &
            & elementExists,elementLocalNumber,ghostElement,err,error,*999)       
          IF(elementExists) THEN
            IF(.NOT.ghostElement) THEN
              numberOfCandidates=numberOfCandidates+1           
              NULLIFY(basis)
              CALL DomainElements_BasisGet(domainElements,elementLocalNumber,basis,err,error,*999)
              CALL Basis_LocalFaceNumberGet(basis,localFaceNormals(elementIdx),localFaceNumber,err,error,*999)
              dataProjection%dataProjectionCandidates(dataPointGlobalNumber)%candidateElementNumbers(numberOfCandidates)= &
                & elementLocalNumber
              dataProjection%dataProjectionCandidates(dataPointGlobalNumber)%localFaceLineNumbers(numberOfCandidates)= &
                & localFaceNumber
            ENDIF
          ELSE
            localError="Element user number "//TRIM(NumberToVString(elementUserNumbers(elementIdx),"*",err,error))// &
              & " does not exist."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ENDDO !elementIdx
        IF(numberOfCandidates>dataProjection%maxNumberOfCandidates) dataProjection%maxNumberOfCandidates=numberOfCandidates
      ENDDO !dataPointIdx
    CASE DEFAULT
      localError="The data projection type of "//TRIM(NumberToVString(dataProjection%projectionType,"",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
     
    EXITS("DataProjection_ProjectionDataCandidateFacesSet")
    RETURN
999 ERRORS("DataProjection_ProjectionDataCandidateFacesSet",err,error)
    EXITS("DataProjection_ProjectionDataCandidateFacesSet")
    RETURN 1

  END SUBROUTINE DataProjection_ProjectionDataCandidateFacesSet
  
  !
  !================================================================================================================================
  !
  
  !>Sets the candidates element numbers and local lines for a data projection.
  SUBROUTINE DataProjection_ProjectionCandidateLinesSet(dataProjection,elementUserNumbers,localLineNormals,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to set the projection type for
    INTEGER(INTG), INTENT(IN) :: elementUserNumbers(:) !<elementUserNumbers(elementIdx). The projection candidate user element numbers
    INTEGER(INTG), INTENT(IN) :: localLineNormals(:,:) !<localLineNormals(normalIdx,elementIdx). The projection candidate line xi normals
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: elementIdx,elementLocalNumber,localLineNumber,numberOfCandidates
    LOGICAL :: elementExists,ghostElement
    TYPE(BASIS_TYPE), POINTER :: basis
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: decompositionTopology
    TYPE(DOMAIN_TYPE), POINTER :: domain
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: domainElements
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: domainTopology
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DataProjection_ProjectionCandidateLinesSet",err,error,*998)

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*998)
    IF(.NOT.ALLOCATED(dataProjection%dataProjectionCandidates)) &
      & CALL FlagError("Data projection candidates is not allocated.",err,error,*998)
     
    SELECT CASE(dataProjection%projectionType)
    CASE(DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE)
      CALL FlagError("Use the candidate faces routine not the candidate lines routine for a boundary faces projection.", &
        & err,error,*998)
    CASE(DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE)
      CALL FlagError("Use the candidate elements routine not the candidate lines routine for an all elements projection.", &
        & err,error,*998)
    CASE(DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE)
      IF(SIZE(elementUserNumbers,1)/=SIZE(localLineNormals,2)) THEN
        localError="The size of the element user numbers array of "// &
          & TRIM(NumberToVString(SIZE(elementUserNumbers,1),"*",err,error))// &
          & " does not match the second size of the line normals array of "// &
          & TRIM(NumberToVString(SIZE(localLineNormals,2),"*",err,error))//"."
        CALL FlagError(localError,err,error,*998)
      ENDIF
      IF(SIZE(localLineNormals,1)/=2) THEN
        localError="The first size of the line normals array of "// &
          & TRIM(NumberToVString(SIZE(localLineNormals,1),"*",err,error))//" is invalid. The size should be 2."
        CALL FlagError(localError,err,error,*998)
      ENDIF
      NULLIFY(decomposition)
      CALL DataProjection_DecompositionGet(dataProjection,decomposition,err,error,*998)
      NULLIFY(decompositionTopology)
      CALL Decomposition_TopologyGet(decomposition,decompositionTopology,err,error,*998)
      NULLIFY(domain)
      CALL Decomposition_DomainGet(decomposition,0,domain,err,error,*998)
      NULLIFY(domainTopology)
      CALL Domain_TopologyGet(domain,domainTopology,err,error,*998)
      NULLIFY(domainElements)
      CALL DomainTopology_ElementsGet(domainTopology,domainElements,err,error,*998)
      IF(ALLOCATED(dataProjection%dataProjectionCandidates(0)%candidateElementNumbers)) &
        & DEALLOCATE(dataProjection%dataProjectionCandidates(0)%candidateElementNumbers)
      IF(ALLOCATED(dataProjection%dataProjectionCandidates(0)%localFaceLineNumbers)) &
        & DEALLOCATE(dataProjection%dataProjectionCandidates(0)%localFaceLineNumbers)
      ALLOCATE(dataProjection%dataProjectionCandidates(0)%candidateElementNumbers(SIZE(elementUserNumbers,1)),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate candidiate element numbers.",err,error,*999)
      ALLOCATE(dataProjection%dataProjectionCandidates(0)%localFaceLineNumbers(SIZE(localLineNormals,1)),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate candidiate local face/line numbers.",err,error,*999)
      numberOfCandidates=0
      DO elementIdx=1,SIZE(elementUserNumbers,1)
        CALL DecompositionTopology_ElementCheckExists(decompositionTopology,elementUserNumbers(elementIdx), &
          & elementExists,elementLocalNumber,ghostElement,err,error,*999)       
        IF(elementExists) THEN
          IF(.NOT.ghostElement) THEN
            numberOfCandidates=numberOfCandidates+1
            NULLIFY(basis)
            CALL DomainElements_BasisGet(domainElements,elementLocalNumber,basis,err,error,*999)
            CALL Basis_LocalLineNumberGet(basis,localLineNormals(:,elementIdx),localLineNumber,err,error,*999)
            dataProjection%dataProjectionCandidates(0)%candidateElementNumbers(elementIdx)=elementLocalNumber
            dataProjection%dataProjectionCandidates(0)%localFaceLineNumbers(elementIdx)=localLineNumber
          ENDIF
        ELSE
          localError="Element user number "//TRIM(NumberToVString(elementUserNumbers(elementIdx),"*",err,error))// &
            & " does not exist."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDDO !elementIdx
      IF(numberOfCandidates>dataProjection%maxNumberOfCandidates) dataProjection%maxNumberOfCandidates=numberOfCandidates
    CASE DEFAULT
      localError="The data projection type of "//TRIM(NumberToVString(dataProjection%projectionType,"",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*998)
    END SELECT
     
    EXITS("DataProjection_ProjectionCandidateLinesSet")
    RETURN
999 IF(ALLOCATED(dataProjection%dataProjectionCandidates(0)%candidateElementNumbers)) &
      & DEALLOCATE(dataProjection%dataProjectionCandidates(0)%candidateElementNumbers)
    IF(ALLOCATED(dataProjection%dataProjectionCandidates(0)%localFaceLineNumbers)) &
      & DEALLOCATE(dataProjection%dataProjectionCandidates(0)%localFaceLineNumbers)
998 ERRORSEXITS("DataProjection_ProjectionCandidateLinesSet",err,error)
    RETURN 1

  END SUBROUTINE DataProjection_ProjectionCandidateLinesSet
  
  !
  !================================================================================================================================
  !
  
  !>Sets the candidates element numbers and local lines for data points of a data projection.
  SUBROUTINE DataProjection_ProjectionDataCandidateLinesSet(dataProjection,dataPointUserNumbers,elementUserNumbers, &
    & localLineNormals,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to set the projection type for
    INTEGER(INTG), INTENT(IN) :: dataPointUserNumbers(:) !<dataPointUserNumbers(dataPointIdx). The data points for the projection candidate user element numbers
    INTEGER(INTG), INTENT(IN) :: elementUserNumbers(:) !<elementUserNumbers(elementIdx). The projection candidate user element numbers
    INTEGER(INTG), INTENT(IN) :: localLineNormals(:,:) !<localLineNormals(normalIdx,elementIdx). The projection candidate line xi normals
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointGlobalNumber,dataPointIdx,elementIdx,elementLocalNumber,localLineNumber,numberOfCandidates
    LOGICAL :: elementExists,ghostElement
    TYPE(BASIS_TYPE), POINTER :: basis
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: decompositionTopology
    TYPE(DOMAIN_TYPE), POINTER :: domain
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: domainElements
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: domainTopology
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DataProjection_ProjectionDataCandidateLinesSet",err,error,*999)

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(dataProjection%dataProjectionCandidates)) &
      & CALL FlagError("Data projection candidates is not allocated.",err,error,*999)
     
    SELECT CASE(dataProjection%projectionType)
    CASE(DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE)
      CALL FlagError("Use the data candidate faces routine not the candidate lines routine for a boundary faces projection.", &
        & err,error,*999)
    CASE(DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE)
      CALL FlagError("Use the data candidate elements routine not the candidate lines routine for an all elements projection.", &
        & err,error,*999)
    CASE(DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE)
      IF(SIZE(elementUserNumbers,1)/=SIZE(localLineNormals,2)) THEN
        localError="The size of the element user numbers array of "// &
          & TRIM(NumberToVString(SIZE(elementUserNumbers,1),"*",err,error))// &
          & " does not match the second size of the line normals array of "// &
          & TRIM(NumberToVString(SIZE(localLineNormals,2),"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      IF(SIZE(localLineNormals,1)/=2) THEN
        localError="The first size of the line normals array of "// &
          & TRIM(NumberToVString(SIZE(localLineNormals,1),"*",err,error))//" is invalid. The size should be 2."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      NULLIFY(decomposition)
      CALL DataProjection_DecompositionGet(dataProjection,decomposition,err,error,*999)
      NULLIFY(decompositionTopology)
      CALL Decomposition_TopologyGet(decomposition,decompositionTopology,err,error,*999)
      NULLIFY(domain)
      CALL Decomposition_DomainGet(decomposition,0,domain,err,error,*999)
      NULLIFY(domainTopology)
      CALL Domain_TopologyGet(domain,domainTopology,err,error,*999)
      NULLIFY(domainElements)
      CALL DomainTopology_ElementsGet(domainTopology,domainElements,err,error,*999)
      DO dataPointIdx=1,SIZE(dataPointUserNumbers,1)
        CALL DataProjection_DataPointGlobalNumberGet(dataProjection,dataPointUserNumbers(dataPointIdx),dataPointGlobalNumber, &
          & err,error,*999)
        IF(ALLOCATED(dataProjection%dataProjectionCandidates(dataPointGlobalNumber)%candidateElementNumbers)) &
          & DEALLOCATE(dataProjection%dataProjectionCandidates(dataPointGlobalNumber)%candidateElementNumbers)
        IF(ALLOCATED(dataProjection%dataProjectionCandidates(dataPointGlobalNumber)%localFaceLineNumbers)) &
          & DEALLOCATE(dataProjection%dataProjectionCandidates(dataPointGlobalNumber)%localFaceLineNumbers)
        ALLOCATE(dataProjection%dataProjectionCandidates(dataPointGlobalNumber)%candidateElementNumbers( &
          & SIZE(elementUserNumbers,1)),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate candidiate element numbers.",err,error,*999)
        ALLOCATE(dataProjection%dataProjectionCandidates(dataPointGlobalNumber)%localFaceLineNumbers( &
          & SIZE(localLineNormals,1)),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate candidiate local face/line numbers.",err,error,*999)
        numberOfCandidates=0
        DO elementIdx=1,SIZE(elementUserNumbers,1)
          CALL DecompositionTopology_ElementCheckExists(decompositionTopology,elementUserNumbers(elementIdx), &
            & elementExists,elementLocalNumber,ghostElement,err,error,*999)       
          IF(elementExists) THEN
            IF(.NOT.ghostElement) THEN
              numberOfCandidates=numberOfCandidates+1
              NULLIFY(basis)
              CALL DomainElements_BasisGet(domainElements,elementLocalNumber,basis,err,error,*999)
              CALL Basis_LocalLineNumberGet(basis,localLineNormals(:,elementIdx),localLineNumber,err,error,*999)
              dataProjection%dataProjectionCandidates(dataPointGlobalNumber)%candidateElementNumbers(elementIdx)=elementLocalNumber
              dataProjection%dataProjectionCandidates(dataPointGlobalNumber)%localFaceLineNumbers(elementIdx)=localLineNumber
            ENDIF
          ELSE
            localError="Element user number "//TRIM(NumberToVString(elementUserNumbers(elementIdx),"*",err,error))// &
              & " does not exist."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ENDDO !elementIdx
        IF(numberOfCandidates>dataProjection%maxNumberOfCandidates) dataProjection%maxNumberOfCandidates=numberOfCandidates
      ENDDO !dataPointIdx
    CASE DEFAULT
      localError="The data projection type of "//TRIM(NumberToVString(dataProjection%projectionType,"",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
     
    EXITS("DataProjection_ProjectionDataCandidateLinesSet")
    RETURN
999 ERRORS("DataProjection_ProjectionDataCandidateLinesSet",err,error)
    EXITS("DataProjection_ProjectionDataCandidateLinesSet")
    RETURN 1

  END SUBROUTINE DataProjection_ProjectionDataCandidateLinesSet
  
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

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionFinished) CALL FlagError("Data projection has not been finished.",err,error,*999)
    
    projectionType=dataProjection%projectionType       
    
    EXITS("DataProjection_ProjectionTypeGet")
    RETURN
999 ERRORSEXITS("DataProjection_ProjectionTypeGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_ProjectionTypeGet

  !
  !================================================================================================================================
  !

  !>Sets the projection type for a data projection.
  SUBROUTINE DataProjection_ProjectionTypeSet(dataProjection,projectionType,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to set the projection type for
    INTEGER(INTG), INTENT(IN) :: projectionType !<the projection type to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    REAL(DP), ALLOCATABLE :: startingXi(:)
    
    ENTERS("DataProjection_ProjectionTypeSet",err,error,*999)

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(dataProjection%dataProjectionFinished) CALL FlagError("Data projection has been finished.",err,error,*999)

    dataProjection%projectionType=projectionType
    SELECT CASE(projectionType)
    CASE (DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE)
      dataProjection%numberOfXi=1
    CASE (DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE)
      dataProjection%numberOfXi=2
    CASE (DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE)
      dataProjection%numberOfXi=dataProjection%decomposition%numberOfDimensions
    CASE DEFAULT
      CALL FlagError("Input projection type is undefined.",err,error,*999)
    END SELECT
    IF(dataProjection%numberOfXi/=SIZE(dataProjection%startingXi,1)) THEN
      ALLOCATE(startingXi(dataProjection%numberOfXi),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate starting xi.",err,error,*999)
      IF(dataProjection%numberOfXi>SIZE(dataProjection%startingXi,1)) THEN
        startingXi(1:SIZE(dataProjection%startingXi,1))=dataProjection%startingXi(1:SIZE(dataProjection%startingXi,1))
        startingXi(SIZE(dataProjection%startingXi,1):dataProjection%numberOfXi)=0.5_DP
      ELSE
        startingXi(1:SIZE(dataProjection%startingXi,1))=dataProjection%startingXi(1:dataProjection%numberOfXi)
      ENDIF
      IF(ALLOCATED(dataProjection%startingXi)) DEALLOCATE(dataProjection%startingXi)
      ALLOCATE(dataProjection%startingXi(dataProjection%numberOfXi),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate data projection starting xi.",err,error,*999)
      dataProjection%startingXi(1:dataProjection%numberOfXi)=startingXi(1:dataProjection%numberOfXi)
      DEALLOCATE(startingXi)
    ENDIF
    
    EXITS("DataProjection_ProjectionTypeSet")
    RETURN
999 IF(ALLOCATED(startingXi)) DEALLOCATE(startingXi)
    ERRORSEXITS("DataProjection_ProjectionTypeSet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_ProjectionTypeSet
    
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

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionFinished) CALL FlagError("Data projection has not been finished.",err,error,*999)
    
    relativeTolerance=dataProjection%relativeTolerance       
   
    EXITS("DataProjection_RelativeToleranceGet")
    RETURN
999 ERRORSEXITS("DataProjection_RelativeToleranceGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_RelativeToleranceGet

  !
  !================================================================================================================================
  !

  !>Sets the relative tolerance for a data projection.
  SUBROUTINE DataProjection_RelativeToleranceSet(dataProjection,relativeTolerance,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to set the relative tolerance for
    REAL(DP), INTENT(IN) :: relativeTolerance !<the relative tolerance to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables    
    ENTERS("DataProjection_RelativeToleranceSet",err,error,*999)

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(dataProjection%dataProjectionFinished) CALL FlagError("Data projection has been finished.",err,error,*999)
    IF(relativeTolerance<0.0_DP) CALL FlagError("Data projection relative tolerance must be a positive real number.",err,error,*999)
    
    dataProjection%relativeTolerance=relativeTolerance
    
    EXITS("DataProjection_RelativeToleranceSet")
    RETURN
999 ERRORSEXITS("DataProjection_RelativeToleranceSet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_RelativeToleranceSet
  

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

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionFinished) CALL FlagError("Data projection has not been finished.",err,error,*999)
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

  !>Sets the starting xi for a data projection.
  SUBROUTINE DataProjection_StartingXiSet(dataProjection,startingXi,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to set the starting xi for
    REAL(DP), INTENT(IN) :: startingXi(:) !<The starting xi to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: xiIdx
    LOGICAL :: validInput
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DataProjection_StartingXiSet",err,error,*999)

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(dataProjection%dataProjectionFinished) CALL FlagError("Data projection has been finished.",err,error,*999)
    IF(SIZE(startingXi,1)/=dataProjection%numberOfXi) THEN
      localError="The size of the specified xi array of "//TRIM(NumberToVString(SIZE(startingXi,1),"*",err,error))// &
        & " is invalid. The size should be "//TRIM(NumberToVString(dataProjection%numberOfXi,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    validInput=.TRUE.
    DO xiIdx=1,dataProjection%numberOfXi
      IF((startingXi(xiIdx)<0.0_DP).OR.(startingXi(xiIdx)>1.0_DP)) THEN
        validInput=.FALSE.
        EXIT !this do
      ENDIF
    ENDDO !xiIdx
    IF(.NOT.validInput) CALL FlagError("Data projection starting xi must be between 0.0 and 1.0.",err,error,*999)
    
    dataProjection%startingXi(1:SIZE(startingXi))=startingXi(1:SIZE(startingXi))
    
    EXITS("DataProjection_StartingXiSet")
    RETURN
999 ERRORSEXITS("DataProjection_StartingXiSet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_StartingXiSet
  
  !
  !================================================================================================================================
  !

  !>Sets the element for a data projection.
  SUBROUTINE DataProjection_ElementSet(dataProjection,dataPointUserNumber,elementUserNumber,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to set the element for
    INTEGER(INTG), INTENT(IN) :: dataPointUserNumber !<data point user number
    INTEGER(INTG), INTENT(IN) :: elementUserNumber !<the user element number to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointGlobalNumber,localElementNumber
    LOGICAL :: elementExists,ghostElement
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DataProjection_ElementSet",err,error,*999)

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    NULLIFY(decomposition)
    CALL DataProjection_DecompositionGet(dataProjection,decomposition,err,error,*999)
    
    CALL DataProjection_DataPointGlobalNumberGet(dataProjection,dataPointUserNumber,dataPointGlobalNumber,err,error,*999)
    CALL DecompositionTopology_ElementCheckExists(decomposition%topology,elementUserNumber,elementExists,localElementNumber, &
      & ghostElement,err,error,*999)
    IF(elementExists) THEN
      IF(ghostElement) THEN
        localError="The specificed user element number of "//TRIM(NumberToVString(elementUserNumber,"*",err,error))// &
          & " cannot be set as a projection element as it is a ghost element."
        CALL FlagError(localError,err,error,*999)
      ELSE
        dataProjection%dataProjectionResults(dataPointGlobalNumber)%elementNumber=elementUserNumber
      ENDIF
    ELSE
      localError="The specificed user element number of "//TRIM(NumberToVString(elementUserNumber,"*",err,error))// &
        & " does not exist."
      CALL FlagError(localError,err,error,*999)
    ENDIF
   
    EXITS("DataProjection_ElementSet")
    RETURN
999 ERRORSEXITS("DataProjection_ElementSet",err,error)
    RETURN 1

  END SUBROUTINE DataProjection_ElementSet

  !
  !================================================================================================================================
  !

  !>Gets the projection distance for a data point identified by a given global number.
  SUBROUTINE DataProjection_ResultDistanceGet(dataProjection,dataPointUserNumber,projectionDistance,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection for which projection result is stored
    INTEGER(INTG), INTENT(IN) :: dataPointUserNumber !<The Data projection user number to get the projection distance for
    REAL(DP), INTENT(OUT) :: projectionDistance !<On exit, the projection distance of the specified data point
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointGlobalNumber
    
    ENTERS("DataProjection_ResultDistanceGet",err,error,*999)

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionFinished) CALL FlagError("Data projection has not been finished.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionProjected) CALL FlagError("Data projection has not been projected.",err,error,*999)
    
    CALL DataProjection_DataPointGlobalNumberGet(dataProjection,dataPointUserNumber,dataPointGlobalNumber,err,error,*999)
    projectionDistance=dataProjection%dataProjectionResults(dataPointGlobalNumber)%distance

    EXITS("DataProjection_ResultDistanceGet")
    RETURN
999 ERRORSEXITS("DataProjection_ResultDistanceGet",err,error)
    RETURN 1

  END SUBROUTINE DataProjection_ResultDistanceGet

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

    IF(ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    
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

    IF(ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    
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
  
  !>Sets the label for a data projection for varying string labels. \see OpenCMISS::Iron::cmfe_DataProjection_LabelSet
  SUBROUTINE DataProjection_LabelSetC(dataProjection,label,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to set the label for
    CHARACTER(LEN=*), INTENT(IN) :: label !<the label to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("DataProjection_LabelSetC",err,error,*999)

    IF(ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    
    dataProjection%label=label

    EXITS("DataProjection_LabelSetC")
    RETURN
999 ERRORSEXITS("DataProjection_LabelSetC",err,error)
    RETURN 1

  END SUBROUTINE DataProjection_LabelSetC

  !
  !================================================================================================================================
  !

  !>Sets the label for a data projection for varying string labels. \see OpenCMISS::Iron::cmfe_DataProjection_LabelSet
  SUBROUTINE DataProjection_LabelSetVS(dataProjection,label,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to set the label for
    TYPE(VARYING_STRING), INTENT(IN) :: label !<the label to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("DataProjection_LabelSetVS",err,error,*999)

    IF(ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    
    dataProjection%label=label

    EXITS("DataProjection_LabelSetVS")
    RETURN
999 ERRORSEXITS("DataProjection_LabelSetVS",err,error)
    RETURN 1

  END SUBROUTINE DataProjection_LabelSetVS

  !
  !================================================================================================================================
  !

  !>Outputs the analysis of data projection results. \see OpenCMISS::Iron::cmfe_DataProjection_ResultAnalysisOutput
  SUBROUTINE DataProjection_ResultAnalysisOutput(dataProjection,filename,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to output the result analysis for.
    CHARACTER(LEN=*), INTENT(IN) :: filename !<If not empty, the filename to output the data projection result analysis to. If empty, the analysis will be output to the standard output.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointExitTag,dataPointIdx,dataPointUserNumber,elementUserNumber,filenameLength,localElementNumber, &
      & localLineFaceNumber,myComputationalNodeNumber,normalIdx1,normalIdx2,numberOfComputationalNodes,outputID
    REAL(DP) :: distance
    CHARACTER(LEN=MAXSTRLEN) :: analFilename,format1,format2,localString
    TYPE(BASIS_TYPE), POINTER :: basis
    TYPE(DataPointsType), POINTER :: dataPoints
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition
    TYPE(DECOMPOSITION_ELEMENTS_TYPE), POINTER :: decompositionElements
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: decompositionTopology
    TYPE(DOMAIN_TYPE), POINTER :: domain
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: domainElements
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: domainTopology
    TYPE(FIELD_TYPE), POINTER :: projectionField
    TYPE(VARYING_STRING) :: localError
        
    ENTERS("DataProjection_ResultAnalysisOutput",err,error,*999)

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionFinished) CALL FlagError("Data projection has not been finished.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionProjected) CALL FlagError("Data projection has not been projected.",err,error,*999)
    
    NULLIFY(dataPoints)
    CALL DataProjection_DataPointsGet(dataProjection,dataPoints,err,error,*999)
    NULLIFY(projectionField)
    CALL DataProjection_ProjectionFieldGet(dataProjection,projectionField,err,error,*999)
    NULLIFY(decomposition)
    CALL Field_DecompositionGet(projectionField,decomposition,err,error,*999)
    NULLIFY(decompositionTopology)
    CALL Decomposition_TopologyGet(decomposition,decompositionTopology,err,error,*999)
    NULLIFY(decompositionElements)
    CALL DecompositionTopology_ElementsGet(decompositionTopology,decompositionElements,err,error,*999)
    NULLIFY(domain)
    CALL Decomposition_DomainGet(decomposition,0,domain,err,error,*999)
    NULLIFY(domainTopology)
    CALL Domain_TopologyGet(domain,domainTopology,err,error,*999)
    NULLIFY(domainElements)
    CALL DomainTopology_ElementsGet(domainTopology,domainElements,err,error,*999)
    numberOfComputationalNodes=ComputationalEnvironment_NumberOfNodesGet(err,error)
    IF(err/=0) GOTO 999
    myComputationalNodeNumber=ComputationalEnvironment_NodeNumberGet(err,error)
    IF(err/=0) GOTO 999
    !Find the correct output ID and open a file if necessary
    filenameLength=LEN_TRIM(filename)
    IF(filenameLength>=1) THEN
      IF(numberOfComputationalNodes>1) THEN
        WRITE(analFilename,('(A,A,I0)')) filename(1:filenameLength),".opdataproj.",myComputationalNodeNumber              
      ELSE
        analFilename=filename(1:filenameLength)//".opdataproj"
      ENDIF
      outputID=IO1_FILE_UNIT
      OPEN(UNIT=outputId,FILE=analFileName(1:LEN_TRIM(analFilename)),STATUS="REPLACE",FORM="FORMATTED",IOSTAT=err)
      IF(ERR/=0) CALL FlagError("Error opening data projection analysis output file.",err,error,*999)            
    ELSE
      outputID=GENERAL_OUTPUT_TYPE
    ENDIF

    CALL WriteString(outputID,"Data projection result analysis:",err,error,*999)
    CALL WriteString(outputID,"",err,error,*999)
    CALL WriteStringValue(outputID,"  Data projection number = ",dataProjection%userNumber,err,error,*999)
    CALL WriteStringValue(outputID,"  Data projection label = ",dataProjection%label,err,error,*999)
    CALL WriteStringValue(outputID,"  Data points number = ",dataPoints%userNumber,err,error,*999)
    CALL WriteStringValue(outputID,"  Data projection type = ",dataProjection%projectionType,err,error,*999)
    CALL WriteStringValue(outputID,"  Projection field number = ",dataProjection%projectionField%USER_NUMBER,err,error,*999)
    CALL WriteString(outputID,"",err,error,*999)
    CALL WriteStringValue(outputID,"  Number of data points = ",dataPoints%numberOfDataPoints,err,error,*999)
    CALL WriteString(outputID,"",err,error,*999)
    IF(dataPoints%numberOfDataPoints>0) THEN
      SELECT CASE(dataProjection%projectionType)
      CASE(DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE)
        localString="  Data pt#     User#  Exit tag  Element#  Line normals            Xi        Vector      Distance"
      CASE(DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE)
        localString="  Data pt#     User#  Exit tag  Element#  Face normal            Xi        Vector      Distance"
      CASE(DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE)
        localString="  Data pt#     User#  Exit tag  Element#            Xi        Vector      Distance"
      CASE DEFAULT
        localError="The data projection type of "// &
          & TRIM(NumberToVString(dataProjection%projectionType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      CALL WriteString(outputID,localString,err,error,*999)
      DO dataPointIdx=1,dataPoints%numberOfDataPoints
        dataPointUserNumber=dataProjection%dataProjectionResults(dataPointIdx)%userNumber
        dataPointExitTag=dataProjection%dataProjectionResults(dataPointIdx)%exitTag
        IF(dataPointExitTag==DATA_PROJECTION_CANCELLED) THEN
          WRITE(localString,'(2X,I8,2X,I8,8X,I2)') dataPointIdx,dataPointUserNumber,dataPointExitTag
          CALL WriteString(outputID,localString,err,error,*999)
        ELSE
          localElementNumber=dataProjection%dataProjectionResults(dataPointIdx)%elementNumber
          elementUserNumber=decompositionElements%elements(localElementNumber)%USER_NUMBER
          localLineFaceNumber=dataProjection%dataProjectionResults(dataPointIdx)%elementLineFaceNumber
          distance=dataProjection%dataProjectionResults(dataPointIdx)%distance
          basis=>domainElements%elements(localElementNumber)%basis
          IF(ASSOCIATED(basis)) THEN
            SELECT CASE(dataProjection%projectionType)
            CASE(DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE)
              normalIdx1=basis%localLineXiNormals(1,localLineFaceNumber)
              normalIdx2=basis%localLineXiNormals(2,localLineFaceNumber)
              WRITE(localString,'(2X,I8,2X,I8,8X,I2,2X,I8,9X,I2,X,I2,2X,E12.5,2X,E12.5,2X,E12.5)') dataPointIdx, &
                & dataPointUserNumber,dataPointExitTag,elementUserNumber,normalIdx1,normalIdx2, &
                & dataProjection%dataProjectionResults(dataPointIdx)%elementXi(1), &
                & dataProjection%dataProjectionResults(dataPointIdx)%projectionVector(1),distance
              format1="(56X,E12.5,2X,E12.5)"
              format2="(70X,E12.5)"
            CASE(DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE)
              normalIdx1=basis%localFaceXiNormal(localLineFaceNumber)
              WRITE(localString,'(2X,I8,2X,I8,8X,I2,2X,I8,11X,I2,2X,E12.5,2X,E12.5,2X,E12.5)') dataPointIdx, &
                & dataPointUserNumber,dataPointExitTag,elementUserNumber,normalIdx1, &
                & dataProjection%dataProjectionResults(dataPointIdx)%elementXi(1), &
                & dataProjection%dataProjectionResults(dataPointIdx)%projectionVector(1),distance
              format1="(55X,E12.5,2X,E12.5)"
              format2="(69X,E12.5)"
            CASE(DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE)
              WRITE(localString,'(2X,I8,2X,I8,8X,I2,2X,I8,2X,E12.5,2X,E12.5,2X,E12.5)') dataPointIdx, &
                & dataPointUserNumber,dataPointExitTag,elementUserNumber, &
                & dataProjection%dataProjectionResults(dataPointIdx)%elementXi(1), &
                & dataProjection%dataProjectionResults(dataPointIdx)%projectionVector(1),distance
              format1="(42X,E12.5,2X,E12.5)"
              format2="(56X,E12.5)"
            CASE DEFAULT
              localError="The data projection type of "// &
                & TRIM(NumberToVString(dataProjection%projectionType,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
            CALL WriteString(outputID,localString,err,error,*999)
            IF(dataProjection%numberOfElementXi>1) THEN
              WRITE(localString,format1) &
                & dataProjection%dataProjectionResults(dataPointIdx)%elementXi(2), &
                & dataProjection%dataProjectionResults(dataPointIdx)%projectionVector(2)
              CALL WriteString(outputID,localString,err,error,*999)
            ENDIF
            IF(dataProjection%numberOfElementXi>2) THEN
              WRITE(localString,format1) &
                & dataProjection%dataProjectionResults(dataPointIdx)%elementXi(3), &
                & dataProjection%dataProjectionResults(dataPointIdx)%projectionVector(3)
              CALL WriteString(outputID,localString,err,error,*999)
            ENDIF
            IF(dataProjection%numberOfCoordinates>dataProjection%numberOfElementXi) THEN
              WRITE(localString,format2) &
                & dataProjection%dataProjectionResults(dataPointIdx)%projectionVector(3)
              CALL WriteString(outputID,localString,err,error,*999)
            ENDIF
          ENDIF
        ENDIF
      ENDDO !dataPointIdx
      CALL WriteString(outputID,"",err,error,*999)
      CALL WriteString(outputID,"  Errors:",err,error,*999)
      CALL WriteString(outputID,"",err,error,*999)
      CALL WriteString(outputID,"  Error type            Value  Data user#",err,error,*999)
      WRITE(localString,'("  RMS error    ",2X,E12.5)') dataProjection%rmsError
      CALL WriteString(outputID,localString,err,error,*999)
      WRITE(localString,'("  Maximum error",2X,E12.5,4X,I8)') dataProjection%maximumError, &
        & dataProjection%dataProjectionResults(dataProjection%maximumErrorDataPoint)%userNumber
      CALL WriteString(outputID,localString,err,error,*999)
      WRITE(localString,'("  Minimum error",2X,E12.5,4X,I8)') dataProjection%minimumError, &
        & dataProjection%dataProjectionResults(dataProjection%minimumErrorDataPoint)%userNumber
      CALL WriteString(outputID,localString,err,error,*999)
      CALL WriteString(outputID,"",err,error,*999)
    ENDIF
    IF(filenameLength>=1) THEN
      CLOSE(UNIT=outputID,IOSTAT=err)
      IF(ERR/=0) CALL FlagError("Error closing data projection analysis output file.",err,error,*999)
    ENDIF

    EXITS("DataProjection_ResultAnalysisOutput")
    RETURN
999 ERRORSEXITS("DataProjection_ResultAnalysisOutput",err,error)
    RETURN 1

  END SUBROUTINE DataProjection_ResultAnalysisOutput

  !
  !================================================================================================================================
  !

  !>Gets the projection element number for a data point identified by a given global number.
  SUBROUTINE DataProjection_ResultElementNumberGet(dataProjection,dataPointUserNumber,projectionElementNumber,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection for which projection result is stored
    INTEGER(INTG), INTENT(IN) :: dataPointUserNumber !<The Data projection user number to get the projection element number for
    INTEGER(INTG), INTENT(OUT) :: projectionElementNumber !<On exit, the projection element number of the specified global data point
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointGlobalNumber
   
    ENTERS("DataProjection_ResultElementNumberGet",err,error,*999)

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionFinished) CALL FlagError("Data projection has not been finished.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionProjected) CALL FlagError("Data projection has not been projected.",err,error,*999)
    
    CALL DataProjection_DataPointGlobalNumberGet(dataProjection,dataPointUserNumber,dataPointGlobalNumber,err,error,*999)
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
  SUBROUTINE DataProjection_ResultElementFaceNumberGet(dataProjection,dataPointUserNumber,projectionElementFaceNumber &
    & ,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection for which projection result is stored
    INTEGER(INTG), INTENT(IN) :: dataPointUserNumber !<The Data projection user number to get the projection element face number for
    INTEGER(INTG), INTENT(OUT) :: projectionElementFaceNumber !<On exit, the projection element face number of the specified global data point
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointGlobalNumber
    
    ENTERS("DataProjection_ResultElementFaceNumberGet",err,error,*999)

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionFinished) CALL FlagError("Data projection has not been finished.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionProjected) CALL FlagError("Data projection has not been projected.",err,error,*999)
    
    !Check if boundary faces projection type was set
    IF(dataProjection%projectionType==DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE) THEN
      CALL DataProjection_DataPointGlobalNumberGet(dataProjection,dataPointUserNumber,dataPointGlobalNumber,err,error,*999)
      projectionElementFaceNumber=dataProjection%dataProjectionResults(dataPointGlobalNumber)%elementLineFaceNumber
    ELSE
      CALL FlagError("Data projection type is not set to a boundary faces projection type.",err,error,*999)
    ENDIF

    EXITS("DataProjection_ResultElementFaceNumberGet")
    RETURN
999 ERRORSEXITS("DataProjection_ResultElementFaceNumberGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_ResultElementFaceNumberGet

  !
  !================================================================================================================================
  !

  !>Gets the projection element line number for a data point identified by a given global number.
  SUBROUTINE DataProjection_ResultElementLineNumberGet(dataProjection,dataPointUserNumber,projectionElementLineNumber &
    & ,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection for which projection result is stored
    INTEGER(INTG), INTENT(IN) :: dataPointUserNumber !<The Data projection user number to get the element line number distance for
    INTEGER(INTG), INTENT(OUT) :: projectionElementLineNumber !<On exit, the projection element line number of the specified global data point
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointGlobalNumber
    
    ENTERS("DataProjection_ResultElementLineNumberGet",err,error,*999)

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionFinished) CALL FlagError("Data projection has not been finished.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionProjected) CALL FlagError("Data projection has not been projected.",err,error,*999)
    
    !Check if boundary lines projection type was set
    IF(dataProjection%projectionType==DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE) THEN
      CALL DataProjection_DataPointGlobalNumberGet(dataProjection,dataPointUserNumber,dataPointGlobalNumber,err,error,*999)
      projectionElementLineNumber=dataProjection%dataProjectionResults(dataPointGlobalNumber)%elementLineFaceNumber
    ELSE
      CALL FlagError("Data projection type is not set to a boundary lines projection type.",err,error,*999)
    ENDIF

    EXITS("DataProjection_ResultElementLineNumberGet")
    RETURN
999 ERRORSEXITS("DataProjection_ResultElementLineNumberGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_ResultElementLineNumberGet

  !
  !================================================================================================================================
  !

  !>Gets the projection exit tag for a data point identified by a given global number.
  SUBROUTINE DataProjection_ResultExitTagGet(dataProjection,dataPointUserNumber,projectionExitTag,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection for which projection result is stored
    INTEGER(INTG), INTENT(IN) :: dataPointUserNumber !<The Data projection user number to get the projection exit tag for
    INTEGER(INTG), INTENT(OUT) :: projectionExitTag !<On exit, the projection exit tag of the specified global data point
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointGlobalNumber
    
    ENTERS("DataProjection_ResultExitTagGet",err,error,*999)

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionFinished) CALL FlagError("Data projection has not been finished.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionProjected) CALL FlagError("Data projection has not been projected.",err,error,*999)
    
    CALL DataProjection_DataPointGlobalNumberGet(dataProjection,dataPointUserNumber,dataPointGlobalNumber,err,error,*999)
    projectionExitTag=dataProjection%dataProjectionResults(dataPointGlobalNumber)%exitTag

    EXITS("DataProjection_ResultExitTagGet")
    RETURN
999 ERRORSEXITS("DataProjection_ResultExitTagGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_ResultExitTagGet

  !
  !================================================================================================================================
  !

  !>Gets the projection xi for a data point identified by a given global number.
  SUBROUTINE DataProjection_ResultXiGet(dataProjection,dataPointUserNumber,projectionXi,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection for which projection result is stored
    INTEGER(INTG), INTENT(IN) :: dataPointUserNumber !<The Data projection user number to get the projection xi for
    REAL(DP), INTENT(OUT) :: projectionXi(:) !<On exit, the projection xi of the specified global data point
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointGlobalNumber
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DataProjection_ResultXiGet",err,error,*999)

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionFinished) CALL FlagError("Data projection has not been finished.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionProjected) CALL FlagError("Data projection has not been projected.",err,error,*999)
    IF(SIZE(projectionXi,1)<dataProjection%numberOfXi) THEN
      localError="The specified projection xi has size of "//TRIM(NumberToVString(SIZE(projectionXi,1),"*",err,error))// &
        & " but it needs to have size of >= "//TRIM(NumberToVString(dataProjection%numberOfXi,"*",err,error))//"." 
      CALL FlagError(localError,err,error,*999)
    ENDIF
        
    CALL DataProjection_DataPointGlobalNumberGet(dataProjection,dataPointUserNumber,dataPointGlobalNumber,err,error,*999)
    projectionXi(1:dataProjection%numberOfXi)=dataProjection%dataProjectionResults(dataPointGlobalNumber)% &
      & xi(1:dataProjection%numberOfXi)

    EXITS("DataProjection_ResultXiGet")
    RETURN
999 ERRORSEXITS("DataProjection_ResultXiGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_ResultXiGet

  !
  !================================================================================================================================
  !

  !>Sets the projection xi for a data point identified by a given global number.
  SUBROUTINE DataProjection_ResultXiSet(dataProjection,dataPointUserNumber,projectionXi,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection for which projection result is stored
    INTEGER(INTG), INTENT(IN) :: dataPointUserNumber !<The Data projection user number to set the projection xi for
    REAL(DP), INTENT(IN) :: projectionXi(:) !<The projection xi of the specified global data point to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointGlobalNumber
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DataProjection_ResultXiSet",err,error,*999)

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionFinished) CALL FlagError("Data projection has not been finished.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionProjected) CALL FlagError("Data projection has not been projected.",err,error,*999)
    IF(SIZE(projectionXi,1)/=dataProjection%numberOfXi) THEN
      localError="The specified projection xi has size of "//TRIM(NumberToVString(SIZE(projectionXi,1),"*",err,error))// &
        & "but it needs to have size of "//TRIM(NumberToVString(dataProjection%numberOfXi,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(ANY(projectionXi<0.0_DP).OR.ANY(projectionXi>1.0_DP)) &
      & CALL FlagError("The specified xi location is invalid. The xi location must be >= 0.0 and <= 1.0",err,error,*999)

    CALL DataProjection_DataPointGlobalNumberGet(dataProjection,dataPointUserNumber,dataPointGlobalNumber,err,error,*999)
    dataProjection%dataProjectionResults(dataPointGlobalNumber)%xi(1:dataProjection%numberOfXi)= &
      & projectionXi(1:dataProjection%numberOfXi)

    EXITS("DataProjection_ResultXiSet")
    RETURN
999 ERRORSEXITS("DataProjection_ResultXiSet",err,error)
    RETURN 1

  END SUBROUTINE DataProjection_ResultXiSet

  !
  !================================================================================================================================
  !

  !>Gets the projection vector for a data point identified by a given global number.
  SUBROUTINE DataProjection_ResultProjectionVectorGet(dataProjection,dataPointUserNumber,projectionVector,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection for which projection result is stored
    INTEGER(INTG), INTENT(IN) :: dataPointUserNumber !<The Data projection user number to get the projection xi for
    REAL(DP), INTENT(OUT) :: projectionVector(:) !<On exit, the projection vector of the specified global data point
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointGlobalNumber
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DataProjection_ResultProjectionVectorGet",err,error,*999)

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionFinished) CALL FlagError("Data projection has not been finished.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionProjected) CALL FlagError("Data projection has not been projected.",err,error,*999)
    IF(SIZE(projectionVector,1)<dataProjection%numberOfCoordinates) THEN
      localError="The specified projection vector has a size of "// &
        & TRIM(NumberToVString(SIZE(projectionVector,1),"*",err,error))//" but it needs to have size of "// &
        & TRIM(NumberToVString(dataProjection%numberOfCoordinates,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
   
    CALL DataProjection_DataPointGlobalNumberGet(dataProjection,dataPointUserNumber,dataPointGlobalNumber,err,error,*999)
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
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DataProjection_UserNumberFind",err,error,*998)

    IF(.NOT.ASSOCIATED(dataPoints)) CALL FlagError("Data points is not associated.",err,error,*999)
    IF(ASSOCIATED(dataProjection)) CALL FlagError("Data projection is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(dataPoints%dataProjections)) &
      & CALL FlagError("Data points data projections is not associated.",err,error,*999)
    
    NULLIFY(dataProjection)
    IF(ALLOCATED(dataPoints%dataProjections%dataProjections)) THEN
      projectionIdx=1
      DO WHILE(projectionIdx<=dataPoints%dataProjections%numberOfDataProjections)
        listDataProjection=>dataPoints%dataProjections%dataProjections(projectionIdx)%ptr
        IF(ASSOCIATED(listDataProjection)) THEN
          IF(listDataProjection%userNumber==userNumber) THEN
            dataProjection=>dataPoints%dataProjections%dataProjections(projectionIdx)%ptr
            EXIT
          ENDIF
        ELSE
          localError="The data points data projections is not associated for projection index "// &
            & TRIM(NumberToVString(projectionIdx,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
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

  !>Finalise a data projections.
  SUBROUTINE DataProjections_Finalise(dataProjections,err,error,*)

    !Argument variables
    TYPE(DataProjectionsType), POINTER :: dataProjections !<A pointer to the data projections to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: projectionIdx
       
    ENTERS("DataProjections_Finalise",err,error,*999)

    IF(ASSOCIATED(dataProjections)) THEN
      IF(ALLOCATED(dataProjections%dataProjections)) THEN
        DO projectionIdx=1,SIZE(dataProjections%dataProjections,1)
          CALL DataProjection_Finalise(dataProjections%dataProjections(projectionIdx)%ptr,err,error,*999)
        ENDDO !projectionIdx
        DEALLOCATE(dataProjections%dataProjections)
      ENDIF
      DEALLOCATE(dataProjections)
    ENDIF
    
    EXITS("DataProjections_Finalise")
    RETURN
999 ERRORSEXITS("DataProjections_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE DataProjections_Finalise
  
  !
  !================================================================================================================================
  !

  !>Initialise a data projections.
  SUBROUTINE DataProjections_Initialise(dataProjections,err,error,*)

    !Argument variables
    TYPE(DataProjectionsType), POINTER :: dataProjections !<A pointer to the data projections to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
       
    ENTERS("DataProjections_Initialise",err,error,*998)

    IF(ASSOCIATED(dataProjections)) CALL FlagError("Data projections is already associated.",err,error,*998)

    ALLOCATE(dataProjections,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate data projections.",err,error,*999)
    NULLIFY(dataProjections%dataPoints)
    dataProjections%numberOfDataProjections=0
    
    EXITS("DataProjections_Initialise")
    RETURN
999 CALL DataProjections_Finalise(dataProjections,dummyErr,dummyError,*998)
998 ERRORSEXITS("DataProjections_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE DataProjections_Initialise
  
  !
  !================================================================================================================================
  !

END MODULE DataProjectionRoutines

