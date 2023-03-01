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
!> Contributor(s): Tim Wu, Chris Bradley, Kumar Mithraratne, Xiani (Nancy) Yan, Prasad Babarenda Gamage
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
  USE ComputationRoutines
  USE ComputationAccessRoutines
  USE Constants
  USE CoordinateSystemAccessRoutines
  USE DataPointAccessRoutines
  USE DataProjectionAccessRoutines
  USE DecompositionRoutines
  USE DecompositionAccessRoutines
  USE DomainMappings
  USE FieldRoutines
  USE FieldAccessRoutines
  USE InputOutput
  USE ISO_VARYING_STRING
  USE Kinds
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


  !Module types

  !Module variables

  !Interfaces

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

  PUBLIC DataProjection_AbsoluteToleranceSet

  PUBLIC DataProjection_CreateFinish,DataProjection_CreateStart
  
  PUBLIC DataProjection_Destroy
  
  PUBLIC DataProjection_DataPointFieldVariableEvaluate
  
  PUBLIC DataProjection_DataPointsProjectionEvaluate

  PUBLIC DataProjection_DataPointsPositionEvaluate
  
  PUBLIC DataProjection_ElementSet
  
  PUBLIC DataProjection_LabelSet

  PUBLIC DataProjection_MaximumInterationUpdateSet

  PUBLIC DataProjection_MaximumNumberOfIterationsSet

  PUBLIC DataProjection_NumberOfClosestElementsSet
  
  PUBLIC DataProjection_ProjectionCancelByDataPoints
  
  PUBLIC DataProjection_ProjectionCancelByDistance
  
  PUBLIC DataProjection_ProjectionCancelByExitTags
  
  PUBLIC DataProjection_ProjectionCandidateElementsSet
  
  PUBLIC DataProjection_ProjectionDataCandidateElementsSet
  
  PUBLIC DataProjection_ProjectionCandidateFacesSet
  
  PUBLIC DataProjection_ProjectionDataCandidateFacesSet
  
  PUBLIC DataProjection_ProjectionCandidateLinesSet
  
  PUBLIC DataProjection_ProjectionDataCandidateLinesSet
  
  PUBLIC DataProjection_ProjectionTypeSet
  
  PUBLIC DataProjection_RelativeToleranceSet

  PUBLIC DataProjection_ResultAnalysisOutput
  
  PUBLIC DataProjection_ResultElementNumberSet
  
  PUBLIC DataProjection_ResultElementFaceNumberSet
  
  PUBLIC DataProjection_ResultElementLineNumberSet
 
  PUBLIC DataProjection_ResultProjectionXiSet

  PUBLIC DataProjection_StartingXiSet

  PUBLIC DataProjections_Initialise,DataProjections_Finalise
  
CONTAINS

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
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DataProjection_AbsoluteToleranceSet",err,error,*999)

    CALL DataProjection_AssertNotFinished(dataProjection,err,error,*999)
    IF(absoluteTolerance<0.0_DP) THEN
      localError="The specified absolute tolerance of "//TRIM(NumberToVString(absoluteTolerance,"*",err,error))// &
        & " is invalid. The tolerance must be >= 0.0."
      CALL FlagError(localError,err,error,*999)
    ENDIF
      
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
    & candidateElements,startingXi,closestElements,closestDistances,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<Data projection problem to evaluate
    TYPE(FieldInterpolatedPointType), POINTER :: interpolatedPoint !<The interpolated point for the data point
    REAL(DP), INTENT(IN) :: dataPointLocation(:) !<dataPointLocation(componentIdx). The location of data point in the region.
    INTEGER(INTG), INTENT(IN) :: numberOfCandidates !<The number of candidate elements.
    INTEGER(INTG), INTENT(IN) :: candidateElements(:) !<candidateElememnts(candidateIdx). The list of candidate local elements for the projection.
    REAL(DP), INTENT(IN) :: startingXi(:) !<startingXi(xiIdx). The starting xi location to find the closest elements for.
    INTEGER(INTG), INTENT(OUT) :: closestElements(:) !<On exit, the list of element local numbers with the shortest distances
    REAL(DP), INTENT(OUT) :: closestDistances(:) !<On exit, the list of shortest distances (squared).
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: closestElementIdx,elementNumber,insertIdx
    INTEGER(INTG) :: numberOfCoordinates !<Region coordinate dimension
    INTEGER(INTG) :: numberOfClosestCandidates !<Number of closest elements to record
    REAL(DP) :: distanceVector(3) !<distance vector between data point and interpolated point, maximum dimension is 3
    REAL(DP) :: distance2 !<Distance squared
      
    ENTERS("DataProjection_ClosestElementsFind",err,error,*999)
    
    CALL DataProjection_AssertIsFinished(dataProjection,err,error,*999)
    
    numberOfCoordinates=dataProjection%numberOfCoordinates
    numberOfClosestCandidates=MIN(numberOfCandidates,SIZE(closestElements,1))
    !loop through the first few elements
    DO closestElementIdx=1,numberOfClosestCandidates
      elementNumber=candidateElements(closestElementIdx)
      CALL Field_InterpolationParametersElementGet(dataProjection%projectionSetType,elementNumber, &
        & interpolatedPoint%interpolationParameters,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
      CALL Field_InterpolateXi(NO_PART_DERIV,startingXi,interpolatedPoint,err,error,*999, &
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
      ENDDO !insertIdx
    ENDDO !closestElementIdx
    !Loop through the rest of the elements
    DO closestElementIdx=numberOfClosestCandidates+1,numberOfCandidates
      elementNumber=candidateElements(closestElementIdx)
      CALL Field_InterpolationParametersElementGet(dataProjection%projectionSetType,elementNumber, &
        & interpolatedPoint%interpolationParameters,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
      CALL Field_InterpolateXi(NO_PART_DERIV,startingXi,interpolatedPoint,err,error,*999, &
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
    & candidateElements,candidateElementFaces,startingXi,closestElements,closestElementFaces,closestDistances,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<Data projection problem to evaluate
    TYPE(FieldInterpolatedPointType), POINTER :: interpolatedPoint !<The interpolated point for the data point
    REAL(DP), INTENT(IN) :: dataPointLocation(:) !<dataPointLocation(componentIdx). The location of data point in the region.
    INTEGER(INTG), INTENT(IN) :: numberOfCandidates !<The number of candidate elements.
    INTEGER(INTG), INTENT(IN) :: candidateElements(:) !<candidateElememnts(candidateIdx). The list of candidate local elements for the projection.
    INTEGER(INTG), INTENT(IN) :: candidateElementFaces(:) !<candidateElementFaces(candidateIdx). The list of candidate faces for the projection.
    REAL(DP), INTENT(IN) :: startingXi(:) !<startingXi(xiIdx). The starting xi location to find the closest faces from. 
    INTEGER(INTG), INTENT(OUT) :: closestElements(:) !<On exit, the list of element local numbers with the shortest distances
    INTEGER(INTG), INTENT(OUT) :: closestElementFaces(:) !<On exit, the list of face numbers with the shortest distances
    REAL(DP), INTENT(OUT) :: closestDistances(:) !<On exit, the list of shortest distances (squared).
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: closestElementIdx,elementNumber,elementFaceNumber,faceNumber,insertIdx
    INTEGER(INTG) :: numberOfCoordinates !<Region coordinate dimension
    INTEGER(INTG) :: numberOfClosestCandidates !<Number of closest elements to record
    REAL(DP) :: distanceVector(3) !<distance vector between data point and interpolated point, maximum dimension is 3
    REAL(DP) :: distance2 !<Distance squared
      
    ENTERS("DataProjection_ClosestFacesFind",err,error,*999)
    
    CALL DataProjection_AssertIsFinished(dataProjection,err,error,*999)
   
    numberOfCoordinates=dataProjection%numberOfCoordinates
    numberOfClosestCandidates=MIN(numberOfCandidates,SIZE(closestElements,1))
    !loop through the first few faces
    DO closestElementIdx=1,numberOfClosestCandidates
      elementNumber=candidateElements(closestElementIdx)
      elementFaceNumber=candidateElementFaces(closestElementIdx)
      faceNumber=interpolatedPoint%interpolationParameters%field%decomposition%topology%elements%elements( &
        & elementNumber)%elementFaces(elementFaceNumber)
      CALL Field_InterpolationParametersFaceGet(dataProjection%projectionSetType,faceNumber,interpolatedPoint% &
        & interpolationParameters,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
      CALL Field_InterpolateXi(NO_PART_DERIV,startingXi,interpolatedPoint,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
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
      ENDDO !insertIdx
    ENDDO !closestElementIdx
    !Loop through the rest of the faces
    DO closestElementIdx=numberOfClosestCandidates+1,numberOfCandidates
      elementNumber=candidateElements(closestElementIdx)
      elementFaceNumber=candidateElementFaces(closestElementIdx)
      faceNumber=interpolatedPoint%interpolationParameters%field%decomposition%topology%elements%elements( &
        & elementNumber)%elementFaces(elementFaceNumber)          
      CALL Field_InterpolationParametersFaceGet(dataProjection%projectionSetType,faceNumber,interpolatedPoint% &
        & interpolationParameters,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
      CALL Field_InterpolateXi(NO_PART_DERIV,startingXi,interpolatedPoint,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE) 
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
        ENDDO !insertIdx
      ENDIF
    ENDDO !closestElementIdx
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
    & candidateElements,candidateElementLines,startingXi,closestElements,closestElementLines,closestDistances,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<Data projection problem to evaluate
    TYPE(FieldInterpolatedPointType), POINTER :: interpolatedPoint !<The interpolated point for the data point 
    REAL(DP), INTENT(IN) :: dataPointLocation(:) !<dataPointLocation(componentIdx). The location of data point in the region.
    INTEGER(INTG), INTENT(IN) :: numberOfCandidates !<The number of candidate elements.
    INTEGER(INTG), INTENT(IN) :: candidateElements(:) !<candidateElememnts(candidateIdx). The list of candidate local elements for the projection.
    INTEGER(INTG), INTENT(IN) :: candidateElementLines(:) !<candidateElementLines(candidateIdx). The list of candidate lines for the projection.
    REAL(DP), INTENT(IN) :: startingXi(:) !<startingXi(xiIdx). The starting xi location to find the closest lines from.
    INTEGER(INTG), INTENT(OUT) :: closestElements(:) !<On exit, the list of element local numbers with the shortest distances
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
    
    CALL DataProjection_AssertIsFinished(dataProjection,err,error,*999)
    
    numberOfCoordinates=dataProjection%numberOfCoordinates
    numberOfClosestCandidates=MIN(numberOfCandidates,SIZE(closestElements,1))
    !loop through the first few lines
    DO closestElementIdx=1,numberOfClosestCandidates
      elementNumber=candidateElements(closestElementIdx)
      elementLineNumber=candidateElementLines(closestElementIdx)
      lineNumber=interpolatedPoint%interpolationParameters%field%decomposition%topology%elements%elements( &
        & elementNumber)%elementLines(elementLineNumber)
      CALL Field_InterpolationParametersLineGet(dataProjection%projectionSetType,lineNumber,interpolatedPoint% &
        & interpolationParameters,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
      CALL Field_InterpolateXi(NO_PART_DERIV,startingXi,interpolatedPoint,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
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
      ENDDO !insertIdx
    ENDDO !closestElementIdx
    !Loop through the rest of the lines
    DO closestElementIdx=numberOfClosestCandidates+1,numberOfCandidates
      elementNumber=candidateElements(closestElementIdx)
      elementLineNumber=candidateElementLines(closestElementIdx)
      lineNumber=interpolatedPoint%interpolationParameters%field%decomposition%topology%elements%elements( &
        & elementNumber)%elementLines(elementLineNumber)          
      CALL Field_InterpolationParametersLineGet(dataProjection%projectionSetType,lineNumber,interpolatedPoint% &
        & interpolationParameters,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
      CALL Field_InterpolateXi(NO_PART_DERIV,startingXi,interpolatedPoint,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE) 
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
        ENDDO !insertIdx
      ENDIF
    ENDDO !closestElementIdx
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
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(VARYING_STRING) :: localError

    ENTERS("DataProjection_CreateFinish",err,error,*999)

    CALL DataProjection_AssertNotFinished(dataProjection,err,error,*999)
    NULLIFY(dataPoints)
    CALL DataProjection_DataPointsGet(dataProjection,dataPoints,err,error,*999)
    CALL DataPoints_AssertIsFinished(dataPoints,err,error,*999)
    NULLIFY(decomposition)
    CALL DataProjection_DecompositionGet(dataProjection,decomposition,err,error,*999)

    !Perform various checks
    SELECT CASE(dataProjection%projectionType)
    CASE(DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE)
      !Check we have lines calculated
      CALL Decomposition_AssertCalculateLines(decomposition,err,error,*999)
    CASE(DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE)
      !Check we have faces calculated
      CALL Decomposition_AssertCalculateFaces(decomposition,err,error,*999)
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
    TYPE(FieldType), POINTER :: projectionField !<A pointer to the field for the data projection
    INTEGER(INTG), INTENT(IN) :: projectionVariableType !<The field variable type of the projection field for the data projection \see FieldRoutines_VariableTypes,FieldRoutines
    TYPE(DataProjectionType), POINTER :: dataProjection !<On exit, a pointer to the created data projection. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataProjectionIdx,dummyErr,numberOfCoordinates,numberOfDimensions
    TYPE(CoordinateSystemType), POINTER :: dataPointsCoordinateSystem,fieldCoordinateSystem
    TYPE(DataProjectionPtrType), ALLOCATABLE :: newDataProjections(:)
    TYPE(DecompositionType), POINTER :: fieldDecomposition
    TYPE(FieldVariableType), POINTER :: projectionFieldVariable
    TYPE(VARYING_STRING) :: dummyError,localError
    
    ENTERS("DataProjection_CreateStart",err,error,*998)

    CALL DataPoints_AssertIsFinished(dataPoints,err,error,*999)
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
    CALL Decomposition_NumberOfDimensionsGet(fieldDecomposition,numberOfDimensions,err,error,*999)
    !Check we are on the same coordinate system
    IF(.NOT.ASSOCIATED(dataPointsCoordinateSystem,fieldCoordinateSystem)) &
      & CALL FlagError("The data points and projection field have different coordinate systems.",err,error,*998)  
    CALL CoordinateSystem_DimensionGet(dataPointsCoordinateSystem,numberOfCoordinates,err,error,*999)
    NULLIFY(dataProjection)
    CALL DataProjection_Initialise(dataProjection,err,error,*999)
    dataProjection%userNumber=dataProjectionUserNumber
    dataProjection%dataPoints=>dataPoints
    dataProjection%projectionField=>projectionField
    dataProjection%projectionVariable=>projectionFieldVariable
    dataProjection%projectionVariableType=projectionVariableType
    dataProjection%decomposition=>fieldDecomposition
    dataProjection%numberOfCoordinates=numberOfCoordinates
    dataProjection%numberOfElementXi=numberOfDimensions
    !Default always project to boundaries faces/lines when decomposition dimension is equal to region dimension.
    !If decomposition dimension is less, project to all elements            
    IF(numberOfDimensions<numberOfCoordinates) THEN !decomposition dimension < data dimension
      dataProjection%numberOfXi=numberOfDimensions
      dataProjection%projectionType=DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE
    ELSE
      SELECT CASE(numberOfDimensions) !decomposition dimension = data dimension
      CASE(1)
        dataProjection%numberOfXi=1
        dataProjection%projectionType=DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE
        dataProjection%numberOfClosestElements=2
      CASE(2) 
        dataProjection%numberOfXi=1
        dataProjection%projectionType=DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE
        dataProjection%numberOfClosestElements=4  
      CASE(3)
        dataProjection%numberOfXi=2
        dataProjection%projectionType=DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE
        dataProjection%numberOfClosestElements=8
      CASE DEFAULT
        localError="The decomposition dimension of "//TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
          & " is invalid. The dimension should be >= 1 and <= 3."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDIF
    ALLOCATE(dataProjection%startingXi(dataProjection%numberOfXi,dataProjection%dataPoints%numberOfDataPoints),STAT=err)
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
    TYPE(TreeNodeType), POINTER :: treeNode
   
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

    CALL DataPoints_NumberOfDataPointsget(dataPoints,numberOfDataPoints,err,error,*999)
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
  
  !>Cancels a data projection result
  SUBROUTINE DataProjection_DataProjectionResultCancel(dataProjectionResult,err,error,*)

    !Argument variables
    TYPE(DataProjectionResultType) :: dataProjectionResult !<The data projection result to cancel 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("DataProjection_DataProjectionResultCancel",err,error,*999)

    CALL DataProjection_DataProjectionResultInitialise(dataProjectionResult,err,error,*999)
    dataProjectionResult%exitTag=DATA_PROJECTION_CANCELLED
    IF(ALLOCATED(dataProjectionResult%xi)) dataProjectionResult%xi=0.0_DP
    IF(ALLOCATED(dataProjectionResult%elementXi)) dataProjectionResult%elementXi=0.0_DP
    IF(ALLOCATED(dataProjectionResult%projectionVector)) dataProjectionResult%projectionVector=0.0_DP  
    
    EXITS("DataProjection_DataProjectionResultCancel")
    RETURN
999 ERRORS("DataProjection_DataProjectionResultCancel",err,error)
    EXITS("DataProjection_DataProjectionResultCancel")
    RETURN 1

  END SUBROUTINE DataProjection_DataProjectionResultCancel
  
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

    dataProjectionResult%userDataPointNumber=0
    dataProjectionResult%globalDataPointNumber=0
    dataProjectionResult%distance=0.0
    dataProjectionResult%elementLocalNumber=0
    dataProjectionResult%elementGlobalNumber=0
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

    dataProjectionResult%userDataPointNumber=0
    dataProjectionResult%globalDataPointNumber=0
    dataProjectionResult%distance=0.0_DP
    dataProjectionResult%elementLocalNumber=0
    dataProjectionResult%elementGlobalNumber=0
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
    INTEGER(INTG) :: dataPointIdx,dummyErr,numberOfDataPoints,dataPointUserNumber
    TYPE(DataPointsType), POINTER :: dataPoints
    TYPE(VARYING_STRING) :: dummyError,localError
    
    ENTERS("DataProjection_DataProjectionResultsInitialise",err,error,*998)

    CALL DataProjection_AssertIsFinished(dataProjection,err,error,*999)
    NULLIFY(dataPoints)    
    CALL DataProjection_DataPointsGet(dataProjection,dataPoints,err,error,*999)
    CALL DataPoints_AssertIsFinished(dataPoints,err,error,*999)
   
    CALL DataPoints_NumberOfDataPointsGet(dataPoints,numberOfDataPoints,err,error,*999)
    ALLOCATE(dataProjection%dataProjectionResults(numberOfDataPoints),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate data projection data projection results.",err,error,*999)
    DO dataPointIdx=1,numberOfDataPoints
      CALL DataProjection_DataProjectionResultInitialise(dataProjection%dataProjectionResults(dataPointIdx),err,error,*999)
      dataProjection%dataProjectionResults(dataPointIdx)%globalDataPointNumber=dataPointIdx
      CALL DataPoints_DataUserNumberGet(dataPoints,dataPointIdx,dataPointUserNumber,err,error,*999)
      dataProjection%dataProjectionResults(dataPointIdx)%userDataPointNumber=dataPointUserNumber
      ALLOCATE(dataProjection%dataProjectionResults(dataPointIdx)%xi(dataProjection%numberOfXi),STAT=err)
      IF(err/=0) THEN
        localError="Could not allocate data projection results xi for data point user number "// &
          & TRIM(NumberToVString(dataPointUserNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      ALLOCATE(dataProjection%dataProjectionResults(dataPointIdx)%elementXi(dataProjection%numberOfElementXi),STAT=err)
      IF(err/=0) THEN
        localError="Could not allocate data projection results element xi for data point user number "// &
          & TRIM(NumberToVString(dataPointUserNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      ALLOCATE(dataProjection%dataProjectionResults(dataPointIdx)%projectionVector(dataProjection%numberOfCoordinates),STAT=err)
      IF(err/=0) THEN
        localError="Could not allocate data projection results projection vector for data point user number "// &
          & TRIM(NumberToVString(dataPointUserNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      dataProjection%dataProjectionResults(dataPointIdx)%xi(1:dataProjection%numberOfXi)= &
        & dataProjection%startingXi(1:dataProjection%numberOfXi,dataPointIdx)
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

    NULLIFY(dataProjections)
    CALL DataProjection_DataProjectionsGet(dataProjection,dataProjections,err,error,*999)
    
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
    INTEGER(INTG) :: dataPointIdx,exitTag,maxDataPoint,minDataPoint,numberOfDataPoints,numberOfPoints
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
    numberOfPoints=0
    CALL DataPoints_NumberOfDataPointsGet(dataPoints,numberOfDataPoints,err,error,*999)
    DO dataPointIdx=1,numberOfDataPoints
      distance=dataProjection%dataProjectionResults(dataPointIdx)%distance
      exitTag=dataProjection%dataProjectionResults(dataPointIdx)%exitTag
      IF(exitTag/=DATA_PROJECTION_CANCELLED) THEN
        numberOfPoints=numberOfPoints+1
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
    IF(numberOfPoints>0) THEN
      rmsError=SQRT(rmsError/REAL(numberOfPoints,DP))
    ELSE
      rmsError=SQRT(rmsError)
    ENDIF
    IF(numberOfPoints==1) THEN
      !Set maxDataPoint and maxError equal to minDataPoint and minError, respectively.
      IF(maxDataPoint==0) THEN
        maxDataPoint=minDataPoint
        maxError=minError
      ELSE
        minDataPoint=maxDataPoint
        minError=maxError
      ENDIF
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
  
  !>Evaluate the data point in a field based on data projection
  SUBROUTINE DataProjection_DataPointFieldVariableEvaluate(dataProjection,fieldVariable,fieldParameterSetType, &
    & dataPointUserNumber,fieldVariableResult,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<Data projection to give the xi locations and element number for the data points
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to be interpolated
    INTEGER(INTG), INTENT(IN) :: fieldParameterSetType !<The parameter set to be interpolated
    INTEGER(INTG), INTENT(IN) :: dataPointUserNumber !<The data point user number to evaluate
    REAL(DP), INTENT(OUT) :: fieldResult(:) !<On exit, the field variable interpolated at the data point.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointGlobalNumber,elementNumber,coordinateIdx
    TYPE(DataPointsType), POINTER :: dataPoints
    TYPE(FieldInterpolatedPointType), POINTER :: interpolatedPoint
    TYPE(FieldInterpolationParametersType), POINTER :: interpolationParameters
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DataProjection_DataPointFieldVariableEvaluate",err,error,*999)

    CALL DataProjection_AssertIsFinished(dataProjection,err,error,*999)
    CALL DataProjection_AssertIsProjected(dataProjection,err,error,*999)
    CALL FieldVariable_NumberOfComponentsGet(fieldVariable,numberOfComponents,err,error,*999)
    
#ifdef WITH_PRECHECKS
    
    IF(.NOT.ALLOCATED(dataProjection%dataProjectionResults)) &
      CALL FlagError("Data projection projection results is not allocated.",err,error,*999)
    IF(SIZE(fieldResult,1)<numberOfComponents) THEN
      localError="The size of the supplied field result array of "//TRIM(NumberToVString(SIZE(fieldResult,1),"*",err,error))// &
        & " is too small for all field variable components. The size needs to be >= "// &
        & TRIM(NumberToVString(numberOfComponents,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    NULLIFY(dataPoints)
    CALL DataProjection_DataPointsGet(dataProjection,dataPoints,err,error,*999)
    CALL DataProjection_DataPointGlobalNumberGet(dataProjection,dataPointUserNumber,dataPointGlobalNumber,err,error,*999)
    
    NULLIFY(interpolationParameters)
    CALL FieldVariable_InterpolationParameterInitialise(fieldVariable,interpolationParameters,err,error,*999)
    NULLIFY(interpolatedPoint)
    CALL Field_InterpolatedPointsInitialise(interpolationParameters,interpolatedPoint,err,error,*999)
    
    elementNumber=dataProjection%dataProjectionResults(dataPointGlobalNumber)%elementLocalNumber
    CALL Field_InterpolationParametersElementGet(fieldParameterSetType,elementNumber,interpolationParameters,err,error,*999)
    CALL Field_InterpolateXi(NO_PART_DERIV,dataProjection%dataProjectionResults(dataPointGlobalNumber)%elementXi, &
      & interpolatedPoint,err,error,*999)
    DO coordinateIdx=1,numberOfComponents
      fieldResult(coordinateIdx)=interpolatedPoint%values(coordinateIdx,NO_PART_DERIV)
    ENDDO !coordinateIdx     
    
    CALL Field_InterpolatedPointFinalise(interpolatedPoint,err,error,*999)
    CALL FieldVariable_InterpolationParametersFinalise(interpolationParameters,err,error,*999)
    
    EXITS("DataProjection_DataPointFieldVariableEvaluate")
    RETURN
999 ERRORSEXITS("DataProjection_DataPointFieldVariableEvaluate",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_DataPointFieldVariableEvaluate
  
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
    INTEGER(INTG) :: computationNodeIdx,dataNumberOfCandidates,dataPointIdx,dataProjectionGlobalNumber,elementIdx,elementNumber, &
      & finishIdx,groupCommunicator,lineFaceIdx,lineFaceNumber,localElementNumber,localLineFaceNumber,MPIClosestDistances, &
      & MPIIError,myComputationNode,myComputationalNode,numberOfCandidates,numberOfClosestCandidates,numberOfDataDimensions, &
      & numberOfDataPoints,numberOfDimensions,numberOfElementFaces,numberOfElementLines,numberOfElements,numberOfFaces, &
      & numberOfGroupComputationNodes,numberOfLines,numberOfLocal,reducedNumberOfCLosestCandidates,startIdx, &
      & totalNumberOfClosestCandidates,totalNumberOfElements,xiIdx
    INTEGER(INTG), ALLOCATABLE :: candidateElements(:),candidateLinesFaces(:),closestElements(:,:),closestLinesFaces(:,:), &
      & globalMPIDisplacements(:),globalNumberOfClosestCandidates(:),globalNumberOfProjectedPoints(:), &
      & globalToLocalNumberOfClosestCandidates(:),projectedGlobalElement(:),projectionExitTag(:),projectedLineFace(:), &
      & projectedLocalElement(:),sortingIndices1(:),sortingIndices2(:)
    REAL(DP) :: position(3)
    REAL(DP), ALLOCATABLE :: closestDistances(:,:),globalClosestDistances(:,:),projectedDistance(:,:),projectedXi(:,:), &
      & projectionVectors(:,:)   
    LOGICAL :: boundaryElement,boundaryFace,boundaryLine,boundaryProjection
    TYPE(BasisType), POINTER :: basis
    TYPE(DataPointsType), POINTER :: dataPoints
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionElementsType), POINTER :: decompositionElements
    TYPE(DecompositionFacesType), POINTER :: decompositionFaces
    TYPE(DecompositionLinesType), POINTER :: decompositionLines
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology    
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainElementsType), POINTER :: domainElements
    TYPE(DomainFacesType), POINTER :: domainFaces
    TYPE(DomainLinesType), POINTER :: domainLines
    TYPE(DomainMappingType), POINTER :: domainMappingElements
    TYPE(DomainMappingsType), POINTER :: domainMappings
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(FieldInterpolationParametersType), POINTER :: interpolationParameters
    TYPE(FieldInterpolatedPointType), POINTER :: interpolatedPoint
    TYPE(FieldParameterSetType), POINTER :: parameterSet
    TYPE(VARYING_STRING) :: localError
    TYPE(WorkGroupType), POINTER :: workGroup
    
    ENTERS("DataProjection_DataPointsProjectionEvaluate",err,error,*999)
    
    CALL DataProjection_AssertIsFinished(dataProjection,err,error,*999)
    NULLIFY(parameterSet)
    CALL Field_ParameterSetGet(dataProjection%projectionField,dataProjection%projectionVariableType, &
      & projectionFieldSetType,parameterSet,err,error,*999)
    
    dataProjection%projectionSetType=projectionFieldSetType
    dataProjectionGlobalNumber=dataProjection%globalNumber
    NULLIFY(dataPoints)
    CALL DataProjection_DataPointsGet(dataProjection,dataPoints,err,error,*999)
    NULLIFY(interpolationParameters)
    CALL FieldVariable_InterpolationParameterInitialise(dataProjection%projectionVariable,interpolationParameters, &
      & err,error,*998,FIELD_GEOMETRIC_COMPONENTS_TYPE)
    NULLIFY(interpolatedPoint)
    CALL Field_InterpolatedPointInitialise(interpolationParameters,interpolatedPoint,err,error,*998,FIELD_GEOMETRIC_COMPONENTS_TYPE)
    NULLIFY(decomposition)
    CALL DataProjection_DecompositionGet(dataProjection,decomposition,err,error,*999)
    CALL Decomposition_NumberOfDimensionsGet(decomposition,numberOfDimensions,err,error,*999)
    NULLIFY(decompositionTopology)
    CALL Decomposition_DecompositionTopologyGet(decomposition,decompositionTopology,err,error,*999)
    NULLIFY(decompositionElements)
    CALL DecompositionTopology_DecompositionElementsGet(decompositionTopology,decompositionElements,err,error,*999)
    CALL DecompositionElements_NumberOfElementsGet(decompositionElements,numberOfElements,err,error,*999)
    CALL DecompositionElements_TotalNumberOfElementsGet(decompositionElements,totalNumberOfElements,err,error,*999)
    NULLIFY(domain)
    CALL Decomposition_DomainGet(decomposition,0,domain,err,error,*999)
    NULLIFY(domainTopology)
    CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
    NULLIFY(domainElements)
    CALL DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*999)
    CALL DomainElements_NumberOfElementsGet(domainElements,numberOfElements,err,error,*999)
    NULLIFY(domainFaces)
    CALL DomainTopology_DomainFacesGet(domainTopology,domainFaces,err,error,*999)
    CALL DomainFaces_NumberOfFacesGet(domainFaces,numberOfFaces,err,error,*999)
    NULLIFY(domainLines)
    CALL DomainTopology_DomainLinesGet(domainTopology,domainLines,err,error,*999)
    CALL DomainLines_NumberOfLinesGet(domainLines,numberOfLines,err,error,*999)
    NULLIFY(domainMappings)
    CALL Domain_DomainMappingsGet(domain,domainMappings,err,error,*999)
    NULLIFY(domainMappingElements)
    CALL DomainMappings_ElementsMappingGet(domainMappings,domainMappingElements,err,error,*999)
    CALL DomainMapping_NumberOfLocalGet(domainMappingElements,numberOfLocal,err,error,*999)
    CALL DataPoints_NumberOfDataPointsGet(dataPoints,numberOfDataPoints,err,error,*999)
    CALL DataPoints_NumberOfDimensionsGet(dataPoints,numberOfDataDimensions,err,error,*999)
    NULLIFY(workGroup)
    CALL Decomposition_WorkGroupGet(decomposition,workGroup,err,error,*999)
    CALL WorkGroup_GroupCommunicatorGet(workGroup,groupCommunicator,err,error,*999)
    CALL WorkGroup_NumberOfGroupNodesGet(workGroup,numberOfGroupComputationNodes,err,error,*999)
    CALL WorkGroup_GroupNodeNumberGet(workGroup,myComputationNode,err,error,*999)
    boundaryProjection=(dataProjection%projectionType==DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE).OR. &
      & (dataProjection%projectionType==DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE)          
    !#########################################################################################################
    !Find elements/faces/lines inside the current computation node, get the boundary faces/lines only if
    !asked the elements/faces/lines are required to perform projection of points in the current computation
    !node the are all pre-allocated to the maximum array length (i.e., numberOfElements), but only up to the
    !numberOfCandidates'th index are assigned
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
          localElementNumber=dataProjection%dataProjectionCandidates(0)%candidateElementNumbers(elementIdx)
          IF(localElementNumber<=numberOfElements) THEN
            !Get non-ghost elements
            numberOfCandidates=numberOfCandidates+1
            candidateElements(numberOfCandidates)=localElementNumber
            !Store element line number for line projection type                      
            candidateLinesFaces(numberOfCandidates)=dataProjection%dataProjectionCandidates(0)%localFaceLineNumbers(elementIdx) 
          ENDIF
        ENDDO !elementIdx
      CASE (DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE)
        IF(numberOfDimensions<2) THEN
          localError="The decomposition number of dimensions of "//TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
            & " is invalid for a faces projection type. The number of dimensions should be >= 2."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        ALLOCATE(candidateElements(numberOfFaces),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate candidate elements.",err,error,*999)
        ALLOCATE(candidateLinesFaces(numberOfFaces),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate candidate lines faces.",err,error,*999)
        !Loop through all candidate elements defined by user number
        DO elementIdx=1,SIZE(dataProjection%dataProjectionCandidates(0)%candidateElementNumbers,1)
          !Check if element exists on current domain, get local number
          localElementNumber=dataProjection%dataProjectionCandidates(0)%candidateElementNumbers(elementIdx)
          IF(localElementNumber>=1.AND.localElementNumber<=numberOfElements) THEN
            !Get non-ghost elements
            numberOfCandidates=numberOfCandidates+1
            candidateElements(numberOfCandidates)=localElementNumber
            !Store element face number for face projection type
            candidateLinesFaces(numberOfCandidates)=dataProjection%dataProjectionCandidates(0)%localFaceLineNumbers(elementIdx)
          ENDIF
        ENDDO !elementIdx
      CASE(DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE)
        !Identify all non-ghost elements
        IF(dataProjection%numberOfXi==numberOfDimensions) THEN
          localError="The data projection number of xi of "//TRIM(NumberToVString(dataProjection%numberOfXi,"*",err,error))// &
            & " does not match the decomposition number of dimensions of "// &
            & TRIM(NumberToVString(numberOfDimensions,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        ALLOCATE(candidateElements(numberOfElements),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate candidate elements.",err,error,*999)
        !Loop through all candidate elements defined by user number                    
        DO elementIdx=1,SIZE(dataProjection%dataProjectionCandidates(0)%candidateElementNumbers,1)
          !Check if element exists on current domain, get local number
          localElementNumber=dataProjection%dataProjectionCandidates(0)%candidateElementNumbers(elementIdx)
          IF(localElementNumber>=1.AND.localElementNumber<=numberOfElements) THEN
            !Get non-ghost elements
            numberOfCandidates=numberOfCandidates+1
            candidateElements(numberOfCandidates)=localElementNumber
          ENDIF
        ENDDO !elementIdx
      CASE DEFAULT
        localError="The data projection type of "//TRIM(NumberToVString(dataProjection%projectionType,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      !If user didn't define candidate element number
      SELECT CASE(dataProjection%projectionType)
      CASE(DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE)
        NULLIFY(decompositionLines)
        CALL DecompositionTopology_DecompositionLinesGet(decompositionTopology,decompositionLines,err,error,*999)
        !identify all non-ghost boundary lines
        ALLOCATE(candidateElements(numberOfLines),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate candidate elements.",err,error,*999)
        ALLOCATE(candidateLinesFaces(numberOfLines),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate candidate lines faces.",err,error,*999)
        DO elementIdx=1,numberOfLocal          
          CALL DecompositionElements_ElementBoundaryElementGet(decompositionElements,elementIdx,boundaryElement,err,error,*999)
          IF(boundaryElement) THEN
            NULLIFY(basis)
            CALL DomainElements_ElementBasisGet(domainElements,elementIdx,basis,err,error,*999)
            CALL Basis_NumberOfLocalLinesGet(basis,numberOfElementLines,err,error,*999)
            DO lineFaceIdx=1,numberOfElementLines
              CALL DecompositionElements_ElementLineNumberGet(decompositionElements,lineFaceIdx,elementIdx,lineFaceNumber, &
                & err,error,*999)
              CALL DecompositionLines_LineBoundaryLineGet(decompositionLines,lineFaceNumber,boundaryLine,err,error,*999)
              IF(boundaryLine) THEN
                numberOfCandidates=numberOfCandidates+1
                candidateLinesFaces(numberOfCandidates)=lineFaceIdx
                candidateElements(numberOfCandidates)=elementIdx
              ENDIF
            ENDDO !lineFaceIdx
          ENDIF
        ENDDO !elementIdx
      CASE(DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE)
        IF(numberOfDimensions<2) THEN
          localError="The decomposition number of dimensions of "//TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
            & " is invalid for a faces projection type. The number of dimensions should be >= 2."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        NULLIFY(decompositionFaces)
        CALL DecompositionTopology_DecompositionFacesGet(decompositionTopology,decompositionFaces,err,error,*999)
        !Identify all non-ghost boundary faces
        ALLOCATE(candidateElements(numberOfFaces),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate candidate elements.",err,error,*999)
        ALLOCATE(candidateLinesFaces(numberOfFaces),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate candidate lines faces.",err,error,*999)
        DO elementIdx=1,numberOfLocal
          CALL DecompositionElements_ElementBoundaryElementGet(decompositionElements,elementIdx,boundaryElement,err,error,*999)
          IF(boundaryElement) THEN
            NULLIFY(basis)
            CALL DomainElements_ElementBasisGet(domainElements,elementIdx,basis,err,error,*999)
            CALL Basis_NumberOfLocalFacesGet(basis,numberOfElementFaces,err,error,*999)
            DO lineFaceIdx=1,numberOfElementFaces
              CALL DecompositionElements_ElementFaceNumberGet(decompositionElements,lineFaceIdx,elementIdx,lineFaceNumber, &
                & err,error,*999)
              CALL DecompositionFaces_FaceBoundaryFaceGet(decompositionFaces,lineFaceNumber,boundaryFace,err,error,*999)
              IF(boundaryFace) THEN
                numberOfCandidates=numberOfCandidates+1
                candidateLinesFaces(numberOfCandidates)=lineFaceIdx
                candidateElements(numberOfCandidates)=elementIdx
              ENDIF
            ENDDO !lineFaceIdx
          ENDIF
        ENDDO !elementIdx
      CASE(DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE)
        !Identify all non-ghost elements
        IF(dataProjection%numberOfXi==numberOfDimensions) THEN
          localError="The data projection number of xi of "//TRIM(NumberToVString(dataProjection%numberOfXi,"*",err,error))// &
            & " does not match the decomposition number of dimensions of "// &
            & TRIM(NumberToVString(numberOfDimensions,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        ALLOCATE(candidateElements(numberOfElements),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate candidate elements.",err,error,*999)
        DO elementIdx=1,numberOfLocal
          numberOfCandidates=numberOfCandidates+1
          candidateElements(numberOfCandidates)=elementIdx
        ENDDO !elementIdx
      CASE DEFAULT
        localError="The data projection type of "//TRIM(NumberToVString(dataProjection%projectionType,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDIF
    IF(numberOfCandidates>dataProjection%maxNumberOfCandidates) dataProjection%maxNumberOfCandidates=numberOfCandidates
    !##############################################################################################################
    !find the clostest elements/faces/lines for each point in the current computation node base on starting xi
    !the clostest elements/faces/lines are required to shrink down on the list of possible projection candiates
    !numberOfClosestCandidates=MIN(dataProjection%numberOfClosestElements,dataProjection%maxNumberOfCandidates)
    numberOfClosestCandidates=dataProjection%maxNumberOfCandidates
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
        CALL DataPoints_DataPositionGet(dataPoints,dataPointIdx,position(1:numberOfDataDimensions),err,error,*999)
        IF(ALLOCATED(dataProjection%dataProjectionCandidates(dataPointIdx)%candidateElementNumbers).AND. &
          & ALLOCATED(dataProjection%dataProjectionCandidates(dataPointIdx)%localFaceLineNumbers)) THEN
          dataNumberOfCandidates=SIZE(dataProjection%dataProjectionCandidates(dataPointIdx)%candidateElementNumbers,1)
          CALL DataProjection_ClosestLinesFind(dataProjection,interpolatedPoint,position(1:numberOfDataDimensions), &
            & dataNumberOfCandidates,dataProjection%dataProjectionCandidates(dataPointIdx)%candidateElementNumbers, &
            & dataProjection%dataProjectionCandidates(dataPointIdx)%localFaceLineNumbers, &
            & dataProjection%startingXi(:,dataPointIdx),closestElements(dataPointIdx,:), &
            & closestLinesFaces(dataPointIdx,:),closestDistances(dataPointIdx,:),err,error,*999)
        ELSE
          CALL DataProjection_ClosestLinesFind(dataProjection,interpolatedPoint,position(1:numberOfDataDimensions), &
            & numberOfCandidates,candidateElements,candidateLinesFaces,dataProjection%startingXi(:,dataPointIdx), &
            & closestElements(dataPointIdx,:),closestLinesFaces(dataPointIdx,:),closestDistances(dataPointIdx,:), &
            & err,error,*999)
        ENDIF
      ENDDO !dataPointIdx
    CASE(DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE)
      !Find closest candidate faces      
      DO dataPointIdx=1,numberOfDataPoints
        CALL DataPoints_DataPositionGet(dataPoints,dataPointIdx,position(1:numberOfDataDimensions),err,error,*999)
        IF(ALLOCATED(dataProjection%dataProjectionCandidates(dataPointIdx)%candidateElementNumbers).AND. &
          & ALLOCATED(dataProjection%dataProjectionCandidates(dataPointIdx)%localFaceLineNumbers)) THEN
          dataNumberOfCandidates=SIZE(dataProjection%dataProjectionCandidates(dataPointIdx)%candidateElementNumbers,1)
          CALL DataProjection_ClosestFacesFind(dataProjection,interpolatedPoint,position(1:numberOfDataDimensions), &
            & dataNumberOfCandidates,dataProjection%dataProjectionCandidates(dataPointIdx)%candidateElementNumbers, &
            & dataProjection%dataProjectionCandidates(dataPointIdx)%localFaceLineNumbers, &
            & dataProjection%startingXi(:,dataPointIdx),closestElements(dataPointIdx,:), &
            & closestLinesFaces(dataPointIdx,:),closestDistances(dataPointIdx,:),err,error,*999)
        ELSE
          CALL DataProjection_ClosestFacesFind(dataProjection,interpolatedPoint,position(1:numberOfDataDimensions), &
            & numberOfCandidates,candidateElements,candidateLinesFaces,dataProjection%startingXi(:,dataPointIdx), &
            & closestElements(dataPointIdx,:),closestLinesFaces(dataPointIdx,:),closestDistances(dataPointIdx,:), &
            & err,error,*999)
        ENDIF
      ENDDO !dataPointIdx
    CASE(DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE)
      !Find closest candidate elements
      DO dataPointIdx=1,numberOfDataPoints
        CALL DataPoints_DataPositionGet(dataPoints,dataPointIdx,position(1:numberOfDataDimensions),err,error,*999)
        IF(ALLOCATED(dataProjection%dataProjectionCandidates(dataPointIdx)%candidateElementNumbers)) THEN
          dataNumberOfCandidates=SIZE(dataProjection%dataProjectionCandidates(dataPointIdx)%candidateElementNumbers,1)
          CALL DataProjection_ClosestElementsFind(dataProjection,interpolatedPoint,position(1:numberOfDataDimensions), &
            & dataNumberOfCandidates,dataProjection%dataProjectionCandidates(dataPointIdx)%candidateElementNumbers, &
            & dataProjection%startingXi(:,dataPointIdx),closestElements(dataPointIdx,:),closestDistances(dataPointIdx,:), &
            & err,error,*999)
        ELSE
          CALL DataProjection_ClosestElementsFind(dataProjection,interpolatedPoint,position(1:numberOfDataDimensions), &
            & numberOfCandidates,candidateElements,dataProjection%startingXi(:,dataPointIdx),closestElements(dataPointIdx,:), &
            & closestDistances(dataPointIdx,:),err,error,*999)
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
    !use MPI if number of computation nodes is greater than 1
    IF(numberOfGroupComputationNodes>1) THEN
      !Use MPI
      !Allocate arrays for MPI communication
      ALLOCATE(globalToLocalNumberOfClosestCandidates(numberOfDataPoints),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate global to local number of closest elements.",err,error,*999)
      ALLOCATE(globalNumberOfClosestCandidates(numberOfGroupComputationNodes),STAT=err) 
      IF(err/=0) CALL FlagError("Could not allocate global number of closest candidates.",err,error,*999)
      ALLOCATE(globalMPIDisplacements(numberOfGroupComputationNodes),STAT=err) 
      IF(err/=0) CALL FlagError("Could not allocate global MPI displacements.",err,error,*999)
      ALLOCATE(globalNumberOfProjectedPoints(numberOfGroupComputationNodes),STAT=err) 
      IF(err/=0) CALL FlagError("Could not allocate all number of projected points.",err,error,*999)
      ALLOCATE(projectionExitTag(numberOfDataPoints),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate projected.",err,error,*999)
      ALLOCATE(projectedLocalElement(numberOfDataPoints),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate projected local element.",err,error,*999)
      ALLOCATE(projectedGlobalElement(numberOfDataPoints),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate projected global element.",err,error,*999)
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
      !gather and distribute the number of closest elements from all computation nodes
      CALL MPI_ALLGATHER(numberOfClosestCandidates,1,MPI_INTEGER,globalNumberOfClosestCandidates,1,MPI_INTEGER, &
        & groupCommunicator,MPIIError)
      CALL MPI_ErrorCheck("MPI_ALLGATHER",MPIIError,err,error,*999)
      !Sum all number of closest candidates from all computation nodes
      totalNumberOfClosestCandidates=SUM(globalNumberOfClosestCandidates,1) 
      !Allocate arrays to store information gathered from all computation node
      !The information for each data point is stored in the corresponding row so they are contiguous in memory for
      !easy MPI access
      ALLOCATE(globalClosestDistances(numberOfDataPoints,totalNumberOfClosestCandidates),STAT=err) 
      IF(err/=0) CALL FlagError("Could not allocate all closest distances.",err,error,*999)
      ALLOCATE(sortingIndices1(totalNumberOfClosestCandidates),STAT=err) 
      IF(err/=0) CALL FlagError("Could not allocate sorting indices 1.",err,error,*999)
      !MPI:create and commit MPI_TYPE_CONTIGUOUS      
      CALL MPI_TYPE_CONTIGUOUS(numberOfDataPoints,MPI_DOUBLE_PRECISION,MPIClosestDistances,MPIIError)
      CALL MPI_ErrorCheck("MPI_TYPE_CONTIGUOUS",MPIIError,err,error,*999)        
      CALL MPI_TYPE_COMMIT(MPIClosestDistances,MPIIError)
      CALL MPI_ErrorCheck("MPI_TYPE_COMMIT",MPIIError,err,error,*999)      
      !Create displacement vectors for MPI_ALLGATHERV
      globalMPIDisplacements(1)=0
      DO computationNodeIdx=1,(numberOfGroupComputationNodes-1)
        globalMPIDisplacements(computationNodeIdx+1)=globalMPIDisplacements(computationNodeIdx)+ &
          & globalNumberOfClosestCandidates(computationNodeIdx)
      ENDDO !computationNodeIdx
      !Share closest element distances between all domains
      CALL MPI_ALLGATHERV(closestDistances(1,1),numberOfClosestCandidates,MPIClosestDistances, &
        & globalClosestDistances,globalNumberOfClosestCandidates,globalMPIDisplacements, &
        & MPIClosestDistances,groupCommunicator,MPIIError)
      CALL MPI_ErrorCheck("MPI_ALLGATHERV",MPIIError,err,error,*999)
      reducedNumberOfCLosestCandidates=MIN(dataProjection%numberOfClosestElements,totalNumberOfClosestCandidates)
      projectedDistance(2,:)=myComputationNode
      !Find the globally closest distance in the current domain
      DO dataPointIdx=1,numberOfDataPoints
        projectionExitTag(dataPointIdx)=dataProjection%dataProjectionResults(dataPointIdx)%exitTag
        projectedLocalElement(dataPointIdx)=dataProjection%dataProjectionResults(dataPointIdx)%elementLocalNumber
        CALL Sorting_BubbleIndexSort(globalClosestDistances(dataPointIdx,:),sortingIndices1,err,error,*999)
        sortingIndices1(1:totalNumberOfClosestCandidates)=sortingIndices1(1:totalNumberOfClosestCandidates)- &
          & globalMPIDisplacements(myComputationNode+1) !shift the index to current computation node
        globalToLocalNumberOfClosestCandidates(dataPointIdx)=0
        DO elementIdx=1,reducedNumberOfCLosestCandidates
          !Sorted index indicates it is in the current computation domain
          IF((sortingIndices1(elementIdx)>=1).AND. &
            & (sortingIndices1(elementIdx)<=globalNumberOfClosestCandidates(myComputationNode+1))) &
            & globalToLocalNumberOfClosestCandidates(dataPointIdx)=globalToLocalNumberOfClosestCandidates(dataPointIdx)+1
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
            CALL DataPoints_DataPositionGet(dataPoints,dataPointIdx,position,err,error,*999)
            projectedXi(:,dataPointIdx) = dataProjection%startingXi(:,dataPointIdx)
            CALL DataProjection_NewtonLinesEvaluate(dataProjection,interpolatedPoint,position(1:numberOfDataDimensions), &
              & closestElements(dataPointIdx,1:numberOfClosestCandidates),closestLinesFaces(dataPointIdx,1: &
              & numberOfClosestCandidates),projectionExitTag(dataPointIdx),projectedLocalElement(dataPointIdx),  &
              & projectedLineFace(dataPointIdx),projectedDistance(1,dataPointIdx),projectedXi(:,dataPointIdx), &
              & projectionVectors(:,dataPointIdx),err,error,*999)
            !Map the element number to global number
            CALL DomainMapping_LocalToGlobalGet(domainMappingElements,projectedLocalElement(dataPointIdx), &
              & projectedGlobalElement(dataPointIdx),err,error,*999)
          ENDIF
        ENDDO !dataPointIdx
      CASE(DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE)
        !Newton project to closest faces, and find miminum projection
        DO dataPointIdx=1,numberOfDataPoints
          numberOfClosestCandidates=globalToLocalNumberOfClosestCandidates(dataPointIdx)
          IF(numberOfClosestCandidates>0) THEN 
            CALL DataPoints_DataPositionGet(dataPoints,dataPointIdx,position,err,error,*999)
            projectedXi(:,dataPointIdx) = dataProjection%startingXi(:,dataPointIdx)
            CALL DataProjection_NewtonFacesEvaluate(dataProjection,interpolatedPoint,position(1:numberOfDataDimensions), &
              & closestElements(dataPointIdx,1:numberOfClosestCandidates),closestLinesFaces(dataPointIdx, &
              & 1:numberOfClosestCandidates),projectionExitTag(dataPointIdx),projectedLocalElement(dataPointIdx), &
              & projectedLineFace(dataPointIdx),projectedDistance(1,dataPointIdx),projectedXi(:,dataPointIdx), &
              & projectionVectors(:,dataPointIdx),err,error,*999)
            !Map the element number to global number
            CALL DomainMapping_LocalToGlobalGet(domainMappingElements,projectedLocalElement(dataPointIdx), &
              & projectedGlobalElement(dataPointIdx),err,error,*999)
          ENDIF
        ENDDO !dataPointIdx
      CASE(DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE)
        !Newton project to closest elements, and find miminum projection
        SELECT CASE(dataProjection%numberOfXi)
        CASE(1) !1D element
          DO dataPointIdx=1,numberOfDataPoints
            numberOfClosestCandidates=globalToLocalNumberOfClosestCandidates(dataPointIdx)
            IF(numberOfClosestCandidates>0) THEN 
              CALL DataPoints_DataPositionGet(dataPoints,dataPointIdx,position,err,error,*999)
              projectedXi(:,dataPointIdx) = dataProjection%startingXi(:,dataPointIdx)
              CALL DataProjection_NewtonElementsEvaluate_1(dataProjection,interpolatedPoint,position(1:numberOfDataDimensions), &
                & closestElements(dataPointIdx,1:numberOfClosestCandidates),projectionExitTag(dataPointIdx), &
                & projectedLocalElement(dataPointIdx),projectedDistance(1,dataPointIdx),projectedXi(:,dataPointIdx), &
                & projectionVectors(:,dataPointIdx),err,error,*999)
              !Map the element number to global number
              CALL DomainMapping_LocalToGlobalGet(domainMappingElements,projectedLocalElement(dataPointIdx), &
                & projectedGlobalElement(dataPointIdx),err,error,*999)
            ENDIF
          ENDDO !dataPointIdx
        CASE(2) !2D element
          DO dataPointIdx=1,numberOfDataPoints
            numberOfClosestCandidates=globalToLocalNumberOfClosestCandidates(dataPointIdx)
            IF(numberOfClosestCandidates>0) THEN 
              CALL DataPoints_DataPositionGet(dataPoints,dataPointIdx,position,err,error,*999)
              projectedXi(:,dataPointIdx) = dataProjection%startingXi(:,dataPointIdx)
              CALL DataProjection_NewtonElementsEvaluate_2(dataProjection,interpolatedPoint,position(1:numberOfDataDimensions), &
                & closestElements(dataPointIdx,1:numberOfClosestCandidates),projectionExitTag(dataPointIdx), &
                & projectedLocalElement(dataPointIdx),projectedDistance(1,dataPointIdx),projectedXi(:,dataPointIdx), &
                & projectionVectors(:,dataPointIdx), err,error,*999)
              !Map the element number to global number                        
              CALL DomainMapping_LocalToGlobalGet(domainMappingElements,projectedLocalElement(dataPointIdx), &
                & projectedGlobalElement(dataPointIdx),err,error,*999)
            ENDIF
          ENDDO !dataPointIdx
        CASE(3) !3D element
          DO dataPointIdx=1,numberOfDataPoints
            numberOfClosestCandidates=globalToLocalNumberOfClosestCandidates(dataPointIdx)
            IF(numberOfClosestCandidates>0) THEN 
              CALL DataPoints_DataPositionGet(dataPoints,dataPointIdx,position,err,error,*999)
              projectedXi(:,dataPointIdx) = dataProjection%startingXi(:,dataPointIdx)
              CALL DataProjection_NewtonElementsEvaluate_3(dataProjection,interpolatedPoint,position(1:numberOfDataDimensions), &
                & closestElements(dataPointIdx,1:numberOfClosestCandidates),projectionExitTag(dataPointIdx), &
                & projectedLocalElement(dataPointIdx),projectedDistance(1,dataPointIdx),projectedXi(:,dataPointIdx), &
                & projectionVectors(:,dataPointIdx),err,error,*999)
              !Map the element number to global number
              CALL DomainMapping_LocalToGlobalGet(domainMappingElements,projectedLocalElement(dataPointIdx), &
                & projectedGlobalElement(dataPointIdx),err,error,*999)
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
        & groupCommunicator,MPIIError)
      CALL MPI_ErrorCheck("MPI_ALLREDUCE",MPIIError,err,error,*999)
      !Sort the computation node/rank from 0 to number of computation node
      CALL Sorting_BubbleIndexSort(projectedDistance(2,:),sortingIndices2,err,error,*999)
      DO computationNodeIdx=0,(numberOfGroupComputationNodes-1)
        globalNumberOfProjectedPoints(computationNodeIdx+1)= &
          & COUNT(ABS(projectedDistance(2,:)-REAL(computationNodeIdx))<ZERO_TOLERANCE)
      ENDDO !computationNodeIdx
      startIdx=SUM(globalNumberOfProjectedPoints(1:myComputationNode))+1
      finishIdx=startIdx+globalNumberOfProjectedPoints(myComputationNode+1)-1
      !create displacement vectors for MPI_ALLGATHERV          
      DO computationNodeIdx=1,(numberOfGroupComputationNodes-1)
        globalMPIDisplacements(computationNodeIdx+1)=globalMPIDisplacements(computationNodeIdx)+ &
          & globalNumberOfProjectedPoints(computationNodeIdx)
      ENDDO !computationNodeIdx  
      !Shares minimum projection information between all domains
      CALL MPI_ALLGATHERV(projectedLocalElement(sortingIndices2(startIdx:finishIdx)),globalNumberOfProjectedPoints( &
        & myComputationNode+1),MPI_INTEGER,projectedLocalElement,globalNumberOfProjectedPoints, &
        & globalMPIDisplacements,MPI_INTEGER,groupCommunicator,MPIIError) !projectedLocalElement
      CALL MPI_ErrorCheck("MPI_ALLGATHERV",MPIIError,err,error,*999)
      CALL MPI_ALLGATHERV(projectedGlobalElement(sortingIndices2(startIdx:finishIdx)),globalNumberOfProjectedPoints( &
        & myComputationNode+1),MPI_INTEGER,projectedGlobalElement,globalNumberOfProjectedPoints, &
        & globalMPIDisplacements,MPI_INTEGER,groupCommunicator,MPIIError) !projectedGlobalElement
      CALL MPI_ErrorCheck("MPI_ALLGATHERV",MPIIError,err,error,*999)
      IF(boundaryProjection) THEN
        CALL MPI_ALLGATHERV(projectedLineFace(sortingIndices2(startIdx:finishIdx)),globalNumberOfProjectedPoints( &
          & myComputationNode+1),MPI_INTEGER,projectedLineFace,globalNumberOfProjectedPoints, &
          & globalMPIDisplacements,MPI_INTEGER,groupCommunicator,MPIIError) !projectedLineFace
        CALL MPI_ErrorCheck("MPI_ALLGATHERV",MPIIError,err,error,*999) 
      ENDIF
      DO xiIdx=1,dataProjection%numberOfXi
        CALL MPI_ALLGATHERV(projectedXi(xiIdx,sortingIndices2(startIdx:finishIdx)),globalNumberOfProjectedPoints( &
          & myComputationNode+1),MPI_DOUBLE_PRECISION,projectedXi(xiIdx,:),globalNumberOfProjectedPoints, &
          & globalMPIDisplacements,MPI_DOUBLE_PRECISION,groupCommunicator,MPIIError) !projectedXi
        CALL MPI_ErrorCheck("MPI_ALLGATHERV",MPIIError,err,error,*999)
      ENDDO !xiIdx
      CALL MPI_ALLGATHERV(projectionExitTag(sortingIndices2(startIdx:finishIdx)),globalNumberOfProjectedPoints( &
        & myComputationNode+1),MPI_INTEGER,projectionExitTag,globalNumberOfProjectedPoints, &
        & globalMPIDisplacements,MPI_INTEGER,groupCommunicator,MPIIError) !projectionExitTag
      CALL MPI_ErrorCheck("MPI_ALLGATHERV",MPIIError,err,error,*999)
      DO xiIdx=1,dataProjection%numberOfCoordinates
        CALL MPI_ALLGATHERV(projectionVectors(xiIdx,sortingIndices2(startIdx:finishIdx)), &
          & globalNumberOfProjectedPoints(myComputationNode+1),MPI_DOUBLE_PRECISION,projectionVectors(xiIdx,:), &
          & globalNumberOfProjectedPoints,globalMPIDisplacements,MPI_DOUBLE_PRECISION,groupCommunicator, &
          & MPIIError)  !projectionVectors
        CALL MPI_ErrorCheck("MPI_ALLGATHERV",MPIIError,err,error,*999)
      ENDDO
      !Assign projection information to projected points
      DO dataPointIdx=1,numberOfDataPoints
        dataProjection%dataProjectionResults(sortingIndices2(dataPointIdx))%exitTag=projectionExitTag(dataPointIdx)
        dataProjection%dataProjectionResults(sortingIndices2(dataPointIdx))%localElementNumber=projectedLocalElement(dataPointIdx)
        dataProjection%dataProjectionResults(sortingIndices2(dataPointIdx))%globalElementNumber=projectedGlobalElement(dataPointIdx)
        dataProjection%dataProjectionResults(dataPointIdx)%distance=projectedDistance(1,dataPointIdx)
        dataProjection%dataProjectionResults(sortingIndices2(dataPointIdx))%xi(1:dataProjection%numberOfXi)= &
          & projectedXi(1:dataProjection%numberOfXi,dataPointIdx)
        dataProjection%dataProjectionResults(sortingIndices2(dataPointIdx))%projectionVector( &
          & 1:dataProjection%numberOfCoordinates)=projectionVectors(1:dataProjection%numberOfCoordinates,dataPointIdx)
      ENDDO !dataPointIdx
      projectedXi(:,sortingIndices2)=projectedXi
      projectionVectors(:,sortingIndices2)=projectionVectors
      projectedLocalElement(sortingIndices2)=projectedLoalElement       
      projectedGlobalElement(sortingIndices2)=projectedGlobalElement       
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
          CALL DataPoints_DataPositionGet(dataPoints,dataPointIdx,position(1:numberOfDataDimensions),err,error,*999)
          CALL DataProjection_NewtonLinesEvaluate(dataProjection,interpolatedPoint,position(1:numberOfDataDimensions), &
            & closestElements(dataPointIdx,:),closestLinesFaces(dataPointIdx,:), &
            & dataProjection%dataProjectionResults(dataPointIdx)%exitTag,dataProjection% &
            & dataProjectionResults(dataPointIdx)%elementLocalNumber,dataProjection% &
            & dataProjectionResults(dataPointIdx)%elementLineFaceNumber,dataProjection%dataProjectionResults( &
            & dataPointIdx)%distance,dataProjection%dataProjectionResults(dataPointIdx)%xi, &
            & dataProjection%dataProjectionResults(dataPointIdx)%projectionVector,err,error,*999)
          dataProjection%dataProjectionResults(dataPointIdx)%elementGlobalNumber= &
            & dataProjection%dataProjectionResults(dataPointIdx)%elementLocalNumber
        ENDDO !dataPointIdx
      CASE(DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE) 
        !Newton project to closest faces, and find miminum projection
        DO dataPointIdx=1,numberOfDataPoints
          CALL DataPoints_DataPositionGet(dataPoints,dataPointIdx,position(1:numberOfDataDimensions),err,error,*999)
          CALL DataProjection_NewtonFacesEvaluate(dataProjection,interpolatedPoint,position(1:numberOfDataDimensions), &
            & closestElements(dataPointIdx,:),closestLinesFaces(dataPointIdx,:), &
            & dataProjection%dataProjectionResults(dataPointIdx)%exitTag,dataProjection% &
            & dataProjectionResults(dataPointIdx)%elementLocalNumber,dataProjection% &
            & dataProjectionResults(dataPointIdx)%elementLineFaceNumber,dataProjection%dataProjectionResults( &
            & dataPointIdx)%distance,dataProjection%dataProjectionResults(dataPointIdx)%xi, &
            & dataProjection%dataProjectionResults(dataPointIdx)%projectionVector,err,error,*999)
          dataProjection%dataProjectionResults(dataPointIdx)%elementGlobalNumber= &
            & dataProjection%dataProjectionResults(dataPointIdx)%elementLocalNumber
        ENDDO !dataPointIdx
      CASE(DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE)        
        !Newton project to closest elements, and find miminum projection
        SELECT CASE(dataProjection%numberOfXi)
        CASE(1) !1D mesh
          DO dataPointIdx=1,numberOfDataPoints
            CALL DataPoints_DataPositionGet(dataPoints,dataPointIdx,position(1:numberOfDataDimensions),err,error,*999)
            CALL DataProjection_NewtonElementsEvaluate_1(dataProjection,interpolatedPoint,position(1:numberOfDataDimensions), &
              & closestElements(dataPointIdx,:),dataProjection%dataProjectionResults(dataPointIdx)%exitTag, &
              & dataProjection%dataProjectionResults(dataPointIdx)%elementLocalNumber, &
              & dataProjection%dataProjectionResults(dataPointIdx)%distance, &
              & dataProjection%dataProjectionResults(dataPointIdx)%xi, &
              & dataProjection%dataProjectionResults(dataPointIdx)%projectionVector,err,error,*999)
            dataProjection%dataProjectionResults(dataPointIdx)%elementGlobalNumber= &
              & dataProjection%dataProjectionResults(dataPointIdx)%elementLocalNumber
          ENDDO !dataPointIdx
        CASE(2) !2D mesh
          DO dataPointIdx=1,numberOfDataPoints
            CALL DataPoints_DataPositionGet(dataPoints,dataPointIdx,position(1:numberOfDataDimensions),err,error,*999)
            CALL DataProjection_NewtonElementsEvaluate_2(dataProjection,interpolatedPoint,position(1:numberOfDataDimensions), &
              & closestElements(dataPointIdx,:),dataProjection%dataProjectionResults(dataPointIdx)%exitTag, &
              & dataProjection%dataProjectionResults(dataPointIdx)%elementLocalNumber, &
              & dataProjection%dataProjectionResults(dataPointIdx)%distance, &
              & dataProjection%dataProjectionResults(dataPointIdx)%xi, &
              & dataProjection%dataProjectionResults(dataPointIdx)%projectionVector,err,error,*999)
            dataProjection%dataProjectionResults(dataPointIdx)%elementGlobalNumber= &
              & dataProjection%dataProjectionResults(dataPointIdx)%elementLocalNumber
          ENDDO !dataPointIdx
        CASE(3) !3D mesh
          DO dataPointIdx=1,numberOfDataPoints
            CALL DataPoints_DataPositionGet(dataPoints,dataPointIdx,position(1:numberOfDataDimensions),err,error,*999)
            projectedXi(:,dataPointIdx) = dataProjection%startingXi(:,dataPointIdx)
            CALL DataProjection_NewtonElementsEvaluate_3(dataProjection,interpolatedPoint,position(1:numberOfDataDimensions), &
              & closestElements(dataPointIdx,:),dataProjection%dataProjectionResults(dataPointIdx)%exitTag, &
              & dataProjection%dataProjectionResults(dataPointIdx)%elementLocalNumber, &
              & dataProjection%dataProjectionResults(dataPointIdx)%distance, &
              & dataProjection%dataProjectionResults(dataPointIdx)%xi, &
              & dataProjection%dataProjectionResults(dataPointIdx)%projectionVector,err,error,*999)
            dataProjection%dataProjectionResults(dataPointIdx)%elementGNumber= &
              & dataProjection%dataProjectionResults(dataPointIdx)%elementLocalNumber
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
    ENDIF !numberOfGroupComputationNodes>1
    !Compute full elemental xi
    IF(dataProjection%numberOfXi==dataProjection%numberOfElementXi) THEN
      DO dataPointIdx=1,numberOfDataPoints
        dataProjection%dataProjectionResults(dataPointIdx)%elementXi=dataProjection%dataProjectionResults(dataPointIdx)%xi
      ENDDO !dataPointIdx
    ELSE
      DO dataPointIdx=1,numberOfDataPoints
        elementNumber=dataProjection%dataProjectionResults(dataPointIdx)%elementLocalNumber
        localLineFaceNumber=dataProjection%dataProjectionResults(dataPointIdx)%elementLineFaceNumber
        NULLIFY(basis)
        CALL DomainElements_ElementBasisGet(domainElements,elementNumber,basis,err,error,*999)
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
    IF(ALLOCATED(projectedLocalElement)) DEALLOCATE(projectedLocalElement)
    IF(ALLOCATED(projectedGlobalElement)) DEALLOCATE(projectedGlobalElement)
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
    IF(ALLOCATED(projectedLocalElement)) DEALLOCATE(projectedLocalElement)
    IF(ALLOCATED(projectedGlobalElement)) DEALLOCATE(projectedGlobalElement)
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
  SUBROUTINE DataProjection_DataPointsPositionEvaluate(dataProjection,field,fieldVariableType_,parameterSetType,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<Data projection to give the xi locations and element number for the data points
    TYPE(FieldType), POINTER :: field !<A pointer to the field to be interpolated
    INTEGER(INTG), INTENT(IN) :: fieldVariableType_ !<The field variable type to be interpolated
    INTEGER(INTG), INTENT(IN) :: parameterSetType !<The parameter set to be interpolated
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointIdx,elementNumber,fieldType_,numberOfDataDimensions,numberOfDataPoints
    TYPE(DataPointsType), POINTER :: dataPoints
    TYPE(FieldInterpolatedPointType), POINTER :: interpolatedPoint
    TYPE(FieldInterpolationParametersType), POINTER :: interpolationParameters
    TYPE(FieldParameterSetType), POINTER :: fieldParameterSet
    TYPE(FieldVariableType), POINTER :: fieldVariable
    
    ENTERS("DataProjection_DataPointsPositionEvaluate",err,error,*999)
    
    CALL DataProjection_AssertIsFinished(dataProjection,err,error,*999)
    CALL DataProjection_AssertIsProjected(dataProjection,err,error,*999)
    IF(.NOT.ALLOCATED(dataProjection%dataProjectionResults)) &
       CALL FlagError("Data projection projection results is not allocated.",err,error,*999)
    NULLIFY(dataPoints)
    CALL DataProjection_DataPointsGet(dataProjection,dataPoints,err,error,*999)
    CALL DataPoints_NumberOfDataPointsGet(dataPoints,numberOfDataPoints,err,error,*999)
    CALL DataPoints_NumberOfDimensionsGet(dataPoints,numberOfDataDimensions,err,error,*999)
    CALL Field_TypeGet(field,fieldType_,err,error,*999)
    IF(.NOT.ASSOCIATED(field)) CALL FlagError("Field is not associated.",err,error,*999)
    IF(fieldType_/=FIELD_GEOMETRIC_TYPE.AND.fieldType_/=FIELD_GEOMETRIC_GENERAL_TYPE) &
      & CALL FlagError("Cannot evaluate data points position on field other than geometric or geometric general type.", &
      & err,error,*999)
    NULLIFY(fieldVariable)
    CALL Field_VariableGet(field,fieldVariableType_,fieldVariable,err,error,*999)
    NULLIFY(fieldParameterSet)
    CALL FieldVariable_ParameterSetGet(fieldVariable,parameterSetType,fieldParameterSet,err,error,*999)
    
    NULLIFY(interpolatedPoint)
    NULLIFY(interpolationParameters)
    CALL FieldVariable_InterpolationParameterInitialise(fieldVariable,interpolationParameters,err,error,*999, &
      & FIELD_GEOMETRIC_COMPONENTS_TYPE)
    CALL Field_InterpolatedPointInitialise(interpolationParameters,interpolatedPoint,err,error,*999, &
      & FIELD_GEOMETRIC_COMPONENTS_TYPE)
    !Loop through data points 
    DO dataPointIdx=1,numberOfDataPoints
      elementNumber=dataProjection%dataProjectionResults(dataPointIdx)%elementLocalNumber
      CALL Field_InterpolationParametersElementGet(parameterSetType,elementNumber,interpolationParameters,err,error,*999, &
        & FIELD_GEOMETRIC_COMPONENTS_TYPE)
      CALL Field_InterpolateXi(NO_PART_DERIV,dataProjection%dataProjectionResults(dataPointIdx)%xi, &
        & interpolatedPoint,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
      dataPoints%dataPoints(dataPointIdx)%position(1:numberOfDataDimensions)= &
        & interpolatedPoint%values(1:numberOfDataDimensions,NO_PART_DERIV)
    ENDDO !dataPointIdx
    CALL Field_InterpolatedPointFinalise(interpolatedPoint,err,error,*999)
    CALL FieldVariable_InterpolationParameterFinalise(interpolationParameters,err,error,*999)
    
    EXITS("DataProjection_DataPointsPositionEvaluate")
    RETURN
999 ERRORSEXITS("DataProjection_DataPointsPositionEvaluate",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_DataPointsPositionEvaluate
  
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

    CALL DataProjection_AssertNotFinished(dataProjection,err,error,*999)
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

    CALL DataProjection_AssertNotFinished(dataProjection,err,error,*999)
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
    TYPE(FieldInterpolatedPointType), POINTER :: interpolatedPoint !<The interpolated point for the data point
    REAL(DP), INTENT(IN) :: dataPointLocation(:) !<dataPointLocation(componentIdx). The location of data point in the region.
    INTEGER(INTG), INTENT(IN) :: candidateElements(:) !<candidateElememnts(candidateIdx). The list of candidate elements for the projection.
    INTEGER(INTG), INTENT(OUT) :: projectionExitTag !<On exit, the exit tag status for the data point projection
    INTEGER(INTG), INTENT(OUT) :: projectionElementNumber !<On exit, the element number of the data point projection
    REAL(DP), INTENT(OUT) :: projectionDistance !<On exit, the projection distance
    REAL(DP), INTENT(INOUT) :: projectionXi(1) !<projectionXi(xiIdx). On exit, the xi location of the data point projection
    REAL(DP), INTENT(OUT) :: projectionVector(3) !<projectionVector(componentIdx). On exit, the projection vector for the data point projection
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string   
    !Local Variables
    INTEGER(INTG) :: bound,elementIdx,elementNumber,exitTag,iterationIdx1,iterationIdx2,numberOfCoordinates
    REAL(DP) :: absoluteTolerance,delta,distanceVector(3),functionGradient,functionHessian,functionValue,maximumDelta, &
      & minimumDelta,newFunctionValue,newXi(1),predictionAccuracy,predictedReduction,relativeTolerance,updateXi(1), &
      & updateXiNorm,xi(1)
    LOGICAL :: insideRegion,converged
    TYPE(FieldInterpolationParametersType), POINTER :: interpolationParameters
   
    ENTERS("DataProjection_NewtonElementsEvaluate_1",err,error,*999)
              
    CALL DataProjection_AssertIsFinished(dataProjection,err,error,*999)
    
    NULLIFY(interpolationParameters)
    CALL FieldInterpolatedPoint_InterpolationParametersGet(interpolatedPoint,interpolationParameters,err,error,*999)
    
    numberOfCoordinates=dataProjection%numberOfCoordinates
    IF(projectionExitTag==DATA_PROJECTION_USER_SPECIFIED) THEN
      IF(projectionElementNumber==0) CALL FlagError("The projection element number has not been set.",err,error,*999)
      CALL Field_InterpolationParametersElementGet(dataProjection%projectionSetType,projectionElementNumber, &
        & interpolationParameters,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
      CALL Field_InterpolateXi(NO_PART_DERIV,projectionXi,interpolatedPoint,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
      projectionVector(1:numberOfCoordinates)=dataPointLocation(1:numberOfCoordinates)- &
        & interpolatedPoint%values(1:numberOfCoordinates,NO_PART_DERIV)
      projectionDistance=DSQRT(DOT_PRODUCT(projectionVector(1:numberOfCoordinates),projectionVector(1:numberOfCoordinates)))
    ELSE
      projectionExitTag=DATA_PROJECTION_EXIT_TAG_NO_ELEMENT
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
          CALL Field_InterpolationParametersElementGet(dataProjection%projectionSetType,elementNumber,interpolationParameters, &
            & err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
          xi=projectionXi
          CALL Field_InterpolateXi(SECOND_PART_DERIV,xi,interpolatedPoint,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
          distanceVector(1:numberOfCoordinates)=dataPointLocation-interpolatedPoint%values(:,NO_PART_DERIV)
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
            functionGradient=-2.0_DP*(DOT_PRODUCT(distanceVector(1:numberOfCoordinates), &
              & interpolatedPoint%values(1:numberOfCoordinates,FIRST_PART_DERIV)))
            !functionHessian 
            functionHessian=-2.0_DP*(DOT_PRODUCT(distanceVector(1:numberOfCoordinates), &
              & interpolatedPoint%values(1:numberOfCoordinates,SECOND_PART_DERIV))- &
              & DOT_PRODUCT(interpolatedPoint%values(1:numberOfCoordinates,FIRST_PART_DERIV), &
              & interpolatedPoint%values(1:numberOfCoordinates,FIRST_PART_DERIV)))
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
              distanceVector(1:numberOfCoordinates)=dataPointLocation(1:numberOfCoordinates)- &
                & interpolatedPoint%values(1:numberOfCoordinates,NO_PART_DERIV)
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
    ENDIF
    
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
    TYPE(FieldInterpolatedPointType), POINTER :: interpolatedPoint !<The interpolated point for the data point 
    REAL(DP), INTENT(IN) :: dataPointLocation(:) !<dataPointLocation(componentIdx). The location of data point in the region.
    INTEGER(INTG), INTENT(IN) :: candidateElements(:) !<candidateElememnts(candidateIdx). The list of candidate elements for the projection.
    INTEGER(INTG), INTENT(OUT) :: projectionExitTag !<On exit, the exit tag status for the data point projection
    INTEGER(INTG), INTENT(OUT) :: projectionElementNumber !<On exit, the element number of the data point projection
    REAL(DP), INTENT(OUT) :: projectionDistance !<On exit, the projection distance
    REAL(DP), INTENT(INOUT) :: projectionXi(2) !<projectionXi(xiIdx). On exit, the xi location of the data point projection
    REAL(DP), INTENT(OUT) :: projectionVector(3) !<projectionVector(componentIdx). On exit, the projection vector for the data point projection
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: bound(2),elementIdx,elementNumber,exitTag,fixedXiIdx,iterationIdx1,iterationIdx2,meshComponentNumber, &
      & numberOfCoordinates,xiIdx
    REAL(DP) :: absoluteTolerance,delta,determinant,distanceVector(3),functionGradient(2),functionGradientNorm, &
      & functionHessian(2,2),functionValue,hessianDiagonal(2),maximumDelta,maxEigen,minimumDelta,minEigen,eigenShift, &
      & newFunctionValue,newXi(2),predictionAccuracy,predictedReduction,relativeTolerance,temp1,temp2,updateXi(2), &
      & updateXiNorm,xi(2)
    LOGICAL :: converged,free,insideRegion
    TYPE(FieldInterpolationParametersType), POINTER :: interpolationParameters
    
    ENTERS("DataProjection_NewtonElementsEvaluate_2",err,error,*999)

    CALL DataProjection_AssertIsFinished(dataProjection,err,error,*999)
    NULLIFY(interpolationParameters)
    CALL FieldInterpolatedPoint_InterpolationParametersGet(interpolatedPoint,interpolationParameters,err,error,*999)
    
    numberOfCoordinates=dataProjection%numberOfCoordinates
    IF(projectionExitTag==DATA_PROJECTION_USER_SPECIFIED) THEN
      IF(projectionElementNumber==0) CALL FlagError("The projection element number has not been set.",err,error,*999)
      CALL Field_InterpolationParametersElementGet(dataProjection%projectionSetType,projectionElementNumber, &
        & interpolationParameters,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
      CALL Field_InterpolateXi(NO_PART_DERIV,projectionXi,interpolatedPoint,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
      projectionVector(1:numberOfCoordinates)=dataPointLocation(1:numberOfCoordinates)- &
        & interpolatedPoint%values(1:numberOfCoordinates,NO_PART_DERIV)
      projectionDistance=DSQRT(DOT_PRODUCT(projectionVector(1:numberOfCoordinates),projectionVector(1:numberOfCoordinates)))
    ELSE
      projectionExitTag=DATA_PROJECTION_EXIT_TAG_NO_ELEMENT
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
            & interpolationParameters,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
          xi=projectionXi
          CALL Field_InterpolateXi(SECOND_PART_DERIV,xi,interpolatedPoint,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
          distanceVector(1:numberOfCoordinates)=dataPointLocation(1:numberOfCoordinates)- &
            & interpolatedPoint%values(1:numberOfCoordinates,NO_PART_DERIV)
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
            functionGradient(1)=-2.0_DP*(DOT_PRODUCT(distanceVector(1:numberOfCoordinates), &
              & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S1)))
            functionGradient(2)=-2.0_DP*(DOT_PRODUCT(distanceVector(1:numberOfCoordinates), &
              & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S2)))
            !functionHessian 
            functionHessian(1,1)=-2.0_DP*(DOT_PRODUCT(distanceVector(1:numberOfCoordinates), &
              & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S1_S1))- &
              & DOT_PRODUCT(interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S1), &
              & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S1))) 
            functionHessian(1,2)=-2.0_DP*(DOT_PRODUCT(distanceVector(1:numberOfCoordinates), &
              & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S1_S2))- &         
              & DOT_PRODUCT(interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S1), &
              & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S2)))
            functionHessian(2,2)=-2.0_DP*(DOT_PRODUCT(distanceVector(1:numberOfCoordinates), &
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
              converged=(updateXiNorm<absoluteTolerance)
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
              distanceVector(1:numberOfCoordinates)=dataPointLocation(1:numberOfCoordinates)- &
                & interpolatedPoint%values(1:numberOfCoordinates,NO_PART_DERIV)
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
                ELSE IF(predictionAccuracy>0.9_DP.AND.predictionAccuracy<1.1_DP) THEN
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
    ENDIF
    
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
    TYPE(FieldInterpolatedPointType), POINTER :: interpolatedPoint !<The interpolated point for the data point 
    REAL(DP), INTENT(IN) :: dataPointLocation(:) !<dataPointLocation(componentIdx). The location of data point in the region.
    INTEGER(INTG), INTENT(IN) :: candidateElements(:) !<candidateElememnts(candidateIdx). The list of candidate elements for the projection.
    INTEGER(INTG), INTENT(OUT) :: projectionExitTag !<On exit, the exit tag status for the data point projection
    INTEGER(INTG), INTENT(OUT) :: projectionElementNumber !<On exit, the element number of the data point projection
    REAL(DP), INTENT(OUT) :: projectionDistance !<On exit, the projection distance
    REAL(DP), INTENT(INOUT) :: projectionXi(3) !<projectionXi(xiIdx). On exit, the xi location of the data point projection
    REAL(DP), INTENT(OUT) :: projectionVector(3) !<projectionVector(componentIdx). On exit, the projection vector for the data point projection
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string   
    !Local Variables
    INTEGER(INTG) :: bound(3),boundXiIdx,elementIdx,elementNumber,exitTag,faceXiIdxs(2),fixedBoundXiIdx,fixedXiIdx, &
      & fixedXiIdx2(2),iterationIdx1,iterationIdx2,nBound,numberOfCoordinates,xiIdx,xiIdx2
    REAL(DP) :: absoluteTolerance,delta,determinant,distanceVector(3),eigenShift,functionGradient(3),functionGradient2(2), &
      & functionGradientNorm,functionHessian2(2,2),functionValue,functionHessian(3,3),hessianDiagonal(3),hessianDiagonal2(2), &
      & maxEigen,maximumDelta,minEigen,minimumDelta,newFunctionValue,newXi(3),predictionAccuracy,predictedReduction, &
      & relativeTolerance,temp1,temp2,temp3,temp4,trace,trace2,updateXi(3),updateXiNorm,xi(3)
    LOGICAL :: free,converged,insideRegion
    TYPE(FieldInterpolationParametersType), POINTER :: interpolationParameters
    
    ENTERS("DataProjection_NewtonElementsEvaluate_3",err,error,*999)
              
    CALL DataProjection_AssertIsFinished(dataProjection,err,error,*999)

    NULLIFY(interpolationParameters)
    CALL FieldInterpolatedPoint_InterpolationParametersGet(interpolatedPoint,interpolationParameters,err,error,*999)
    
    numberOfCoordinates=dataProjection%numberOfCoordinates
    IF(projectionExitTag==DATA_PROJECTION_USER_SPECIFIED) THEN
      IF(projectionElementNumber==0) CALL FlagError("The projection element number has not been set.",err,error,*999)
      CALL Field_InterpolationParametersElementGet(dataProjection%projectionSetType,projectionElementNumber, &
        & interpolatedPoint%INTERPOLATION_PARAMETERS,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
      CALL Field_InterpolateXi(NO_PART_DERIV,projectionXi,interpolatedPoint,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
      projectionVector(1:numberOfCoordinates)=dataPointLocation(1:numberOfCoordinates)- &
        & interpolatedPoint%values(1:numberOfCoordinates,NO_PART_DERIV)
      projectionDistance=DSQRT(DOT_PRODUCT(projectionVector(1:numberOfCoordinates),projectionVector(1:numberOfCoordinates)))
    ELSE
      projectionExitTag=DATA_PROJECTION_EXIT_TAG_NO_ELEMENT
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
            & interpolationParameters,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
          xi=projectionXi
          CALL Field_InterpolateXi(SECOND_PART_DERIV,xi,interpolatedPoint,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
          distanceVector(1:numberOfCoordinates)=dataPointLocation(1:numberOfCoordinates)- &
            interpolatedPoint%values(1:numberOfCoordinates,NO_PART_DERIV)
          functionValue=DOT_PRODUCT(distanceVector(1:numberOfCoordinates),distanceVector(1:numberOfCoordinates))
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
            functionGradient(1)=-2.0_DP*(DOT_PRODUCT(distanceVector(1:numberOfCoordinates), &
              & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S1)))
            functionGradient(2)=-2.0_DP*(DOT_PRODUCT(distanceVector(1:numberOfCoordinates), &
              & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S2)))
            functionGradient(3)=-2.0_DP*(DOT_PRODUCT(distanceVector(1:numberOfCoordinates), &
              & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S3)))
            !functionHessian 
            functionHessian(1,1)=-2.0_DP*(DOT_PRODUCT(distanceVector(1:numberOfCoordinates), &
              & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S1_S1))- &
              & DOT_PRODUCT(interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S1), &
              & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S1))) 
            functionHessian(1,2)=-2.0_DP*(DOT_PRODUCT(distanceVector(1:numberOfCoordinates), &
              & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S1_S2))- &         
              & DOT_PRODUCT(interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S1), &
              & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S2)))
            functionHessian(1,3)=-2.0_DP*(DOT_PRODUCT(distanceVector(1:numberOfCoordinates), &
              & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S1_S3))- &         
              & DOT_PRODUCT(interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S1), &
              & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S3)))
            functionHessian(2,2)=-2.0_DP*(DOT_PRODUCT(distanceVector(1:numberOfCoordinates), &
              & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S2_S2))- &
              & DOT_PRODUCT(interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S2), &
              & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S2)))
            functionHessian(2,3)=-2.0_DP*(DOT_PRODUCT(distanceVector(1:numberOfCoordinates), &
              & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S2_S3))- &
              & DOT_PRODUCT(interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S2), &
              & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S3)))
            functionHessian(3,3)=-2.0_DP*(DOT_PRODUCT(distanceVector(1:numberOfCoordinates), &
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
              converged=(updateXiNorm<absoluteTolerance)
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
                    ENDDO !xiIdx2
                  ELSE IF(newXi(xiIdx)>1.0_DP) THEN
                    newXi(xiIdx)=1.0_DP
                    DO xiIdx2 = 1,3
                      IF(xiIdx2 /= xiIdx) THEN
                        newXi(xiIdx2)=xi(xiIdx2)+updateXi(xiIdx2)
                      ENDIF
                    ENDDO !xiIdx2
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
              distanceVector(1:numberOfCoordinates)=dataPointLocation(1:numberOfCoordinates)- &
                & interpolatedPoint%values(1:numberOfCoordinates,NO_PART_DERIV)
              newFunctionValue=DOT_PRODUCT(distanceVector(1:numberOfCoordinates),distanceVector(1:numberOfCoordinates))
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
                ELSE IF(predictionAccuracy>0.9_DP.AND.predictionAccuracy<1.1_DP) THEN
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
    ENDIF
    
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
    TYPE(FieldInterpolatedPointType), POINTER :: interpolatedPoint !<The interpolation for the data point
    REAL(DP), INTENT(IN) :: dataPointLocation(:) !<dataPointLocation(componentIdx). The location of data point in the region.
    INTEGER(INTG), INTENT(IN) :: candidateElements(:)!<candidateElememnts(candidateIdx). The list of candidate local elements for the projection.
    INTEGER(INTG), INTENT(IN) :: candidateElementFaces(:) !<candidateElementFaces(candidateIdx). The list of candidate faces for the projection.
    INTEGER(INTG), INTENT(INOUT) :: projectionExitTag !<On exit, the exit tag status for the data point projection
    INTEGER(INTG), INTENT(OUT) :: projectionElementNumber !<On exit, the element local number of the data point projection
    INTEGER(INTG), INTENT(OUT) :: projectionElementFaceNumber !<On exit, the face number of the data point projection
    REAL(DP), INTENT(OUT) :: projectionDistance !<On exit, the projection distance
    REAL(DP), INTENT(INOUT) :: projectionXi(2) !<projectionXi(xiIdx). On exit, the xi location of the data point projection
    REAL(DP), INTENT(OUT) :: projectionVector(3) !<projectionVector(componentIdx). On exit, the projection vector for the data point projection
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string   
    !Local Variables
    INTEGER(INTG) :: bound(2),elementNumber,elementFaceNumber,elementIdx,exitTag,faceNumber,fixedXiIdx,iterationIdx1, &
      & iterationIdx2,numberOfCoordinates,xiIdx
    REAL(DP) :: absoluteTolerance,delta,determinant,distanceVector(3),eigenShift,functionGradient(2),functionGradientNorm, &
      & functionHessian(2,2),functionValue,hessianDiagonal(2),maximumDelta,maxEigen,minimumDelta,minEigen,newFunctionValue, &
      & newXi(2),predictionAccuracy,predictedReduction,relativeTolerance,temp1,temp2,updateXi(2),updateXiNorm,xi(2)
    LOGICAL :: free,converged,insideRegion
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionElementsType), POINTER :: decompositionElements
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology
    TYPE(FieldType), POINTER :: field
    TYPE(FieldInterpolationParametersType), POINTER :: interpolationParameters
    
    ENTERS("DataProjection_NewtonFacesEvaluate",err,error,*999)

    CALL DataProjection_AssertIsFinished(dataProjection,err,error,*999)

    NULLIFY(interpolationParameters)
    CALL FieldInterpolatedPoint_InterpolationParametersGet(interpolatedPoint,interpolationParameters,err,error,*999)
    NULLIFY(field)
    CALL FieldInterpolationParameters_FieldGet(interpolationParameters,field,err,error,*999)
    NULLIFY(decomposition)
    CALL Field_DecompositionGet(field,decomposition,err,error,*999)
    NULLIFY(decompositionTopology)
    CALL Decomposition_DecompositionTopologyGet(decomposition,decompositionTopology,err,error,*999)
    NULLIFY(decompositionElements)
    CALL DecompositionTopology_DecompositionElementsGet(decompositionTopology,decompositionElements,err,error,*999)
    
    numberOfCoordinates=dataProjection%numberOfCoordinates
    IF(projectionExitTag==DATA_PROJECTION_USER_SPECIFIED) THEN
      IF(projectionElementNumber==0) CALL FlagError("The projection element number has not been set.",err,error,*999)
      IF(projectionElementFaceNumber==0) CALL FlagError("The projection element face number has not been set.",err,error,*999)
      faceNumber=interpolatedPoint%INTERPOLATION_PARAMETERS%field%decomposition%topology%elements% &
        & elements(projectionElementNumber)%ELEMENT_FACES(projectionElementFaceNumber)
      CALL Field_InterpolationParametersFaceGet(dataProjection%projectionSetType,faceNumber,interpolatedPoint% &
        & INTERPOLATION_PARAMETERS,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
      CALL Field_InterpolateXi(NO_PART_DERIV,projectionXi,interpolatedPoint,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
      projectionVector(1:numberOfCoordinates)=dataPointLocation(1:numberOfCoordinates)- &
        & interpolatedPoint%values(1:numberOfCoordinates,NO_PART_DERIV)
      projectionDistance=DSQRT(DOT_PRODUCT(projectionVector(1:numberOfCoordinates),projectionVector(1:numberOfCoordinates)))
    ELSE
      projectionExitTag=DATA_PROJECTION_EXIT_TAG_NO_ELEMENT
      relativeTolerance=dataProjection%relativeTolerance
      absoluteTolerance=dataProjection%absoluteTolerance
      maximumDelta=dataProjection%maximumIterationUpdate
      minimumDelta=0.025_DP*maximumDelta !need to set a minimum, in case if it gets too small
      !Project on each candidate elements
      DO elementIdx=1,SIZE(candidateElements,1) 
        elementNumber=candidateElements(elementIdx)
        IF(elementNumber>0) THEN
          elementFaceNumber=candidateElementFaces(elementIdx)
          CALL DecompositionElements_ElementFaceNumberGet(decompositionElements,elementFaceNumber,elementNumber,faceNumber, &
            & err,error,*999)
          exitTag=DATA_PROJECTION_EXIT_TAG_NO_ELEMENT
          converged=.FALSE.
          !start at half the maximumDelta as we do not know if quadratic model is a good approximation yet
          delta=0.5_DP*maximumDelta 
          CALL Field_InterpolationParametersFaceGet(dataProjection%projectionSetType,faceNumber,interpolatedPoint% &
            & interpolationParameters,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
          xi=projectionXi
          CALL Field_InterpolateXi(SECOND_PART_DERIV,xi,interpolatedPoint,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
          distanceVector(1:numberOfCoordinates)=dataPointLocation(1:numberOfCoordinates)- &
            & interpolatedPoint%values(1:numberOfCoordinates,NO_PART_DERIV)
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
            functionGradient(1)=-2.0_DP*(DOT_PRODUCT(distanceVector(1:numberOfCoordinates), &
              & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S1)))
            functionGradient(2)=-2.0_DP*(DOT_PRODUCT(distanceVector(1:numberOfCoordinates), &
              & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S2)))
            !functionHessian 
            functionHessian(1,1)=-2.0_DP*(DOT_PRODUCT(distanceVector(1:numberOfCoordinates), &
              & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S1_S1))- &
              & DOT_PRODUCT(interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S1), &
              & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S1))) 
            functionHessian(1,2)=-2.0_DP*(DOT_PRODUCT(distanceVector(1:numberOfCoordinates), &
              & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S1_S2))- &         
              & DOT_PRODUCT(interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S1), &
              & interpolatedPoint%values(1:numberOfCoordinates,PART_DERIV_S2)))
            functionHessian(2,2)=-2.0_DP*(DOT_PRODUCT(distanceVector(1:numberOfCoordinates), &
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
              converged=(updateXiNorm<absoluteTolerance)
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
              distanceVector(1:numberOfCoordinates)=dataPointLocation(1:numberOfCoordinates)- &
                & interpolatedPoint%values(1:numberOfCoordinates,NO_PART_DERIV)
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
                ELSE IF(predictionAccuracy>0.9_DP.AND.predictionAccuracy<1.1_DP) THEN
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
    ENDIF
    
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
    TYPE(FieldInterpolatedPointType), POINTER :: interpolatedPoint !<The interpolated point for the data point
    REAL(DP), INTENT(IN) :: dataPointLocation(:) !<dataPointLocation(componentIdx). The location of data point in the region.
    INTEGER(INTG), INTENT(IN) :: candidateElements(:) !<candidateElememnts(candidateIdx). The list of candidate local elements for the projection.
    INTEGER(INTG), INTENT(IN) :: candidateElementLines(:) !<candidateElementLines(candidateIdx). The list of candidate lines for the projection.
    INTEGER(INTG), INTENT(INOUT) :: projectionExitTag !<On exit, the exit tag status for the data point projection
    INTEGER(INTG), INTENT(OUT) :: projectionElementNumber !<On exit, the element local number of the data point projection
    INTEGER(INTG), INTENT(OUT) :: projectionElementLineNumber !<On exit, the line number of the data point projection
    REAL(DP), INTENT(OUT) :: projectionDistance !<On exit, the projection distance
    REAL(DP), INTENT(INOUT) :: projectionXi(1) !<projectionXi(xiIdx). On exit, the xi location of the data point projection
    REAL(DP), INTENT(OUT) :: projectionVector(3) !<projectionVector(componentIdx). On exit, the projection vector for the data point projection
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string   
    !Local Variables
    INTEGER(INTG) :: bound,elementIdx,elementLineNumber,elementNumber,exitTag,iterationIdx1,iterationIdx2,lineNumber, &
      & numberOfCoordinates
    REAL(DP) :: absoluteTolerance,delta,distanceVector(3),functionGradient,functionHessian,functionValue,maximumDelta, &
      & minimumDelta,newFunctionValue,newXi(1),predictionAccuracy,predictedReduction,relativeTolerance,updateXi(1), &
      & updateXiNorm,xi(1)
    LOGICAL :: insideRegion,converged
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionElementsType), POINTER :: decompositionElements
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology
    TYPE(FieldType), POINTER :: field
    TYPE(FieldInterpolationParametersType), POINTER :: interpolationParameters
    
    ENTERS("DataProjection_NewtonLinesEvaluate",err,error,*999)
              
    CALL DataProjection_AssertIsFinished(dataProjection,err,error,*999)

    NULLIFY(interpolationParameters)
    CALL FieldInterpolatedPoint_InterpolationParametersGet(interpolatedPoint,interpolationParameters,err,error,*999)
    NULLIFY(field)
    CALL FieldInterpolationParameters_FieldGet(interpolationParameters,field,err,error,*999)
    NULLIFY(decomposition)
    CALL Field_DecompositionGet(field,decomposition,err,error,*999)
    NULLIFY(decompositionTopology)
    CALL Decomposition_DecompositionTopologyGet(decomposition,decompositionTopology,err,error,*999)
    NULLIFY(decompositionElements)
    CALL DecompositionTopology_DecompositionElementsGet(decompositionTopology,decompositionElements,err,error,*999)
    
    numberOfCoordinates=dataProjection%numberOfCoordinates
    IF(projectionExitTag==DATA_PROJECTION_USER_SPECIFIED) THEN
      IF(projectionElementNumber==0) CALL FlagError("The projection element number has not been set.",err,error,*999)
      IF(projectionElementLineNumber==0) CALL FlagError("The projection element line number has not been set.",err,error,*999)
      lineNumber=interpolatedPoint%INTERPOLATION_PARAMETERS%field%decomposition%topology%elements% &
        & elements(projectionElementNumber)%ELEMENT_LINES(projectionElementLineNumber)
      CALL Field_InterpolationParametersLineGet(dataProjection%projectionSetType,lineNumber,interpolationParameters, &
        & err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
      CALL Field_InterpolateXi(NO_PART_DERIV,projectionXi,interpolatedPoint,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
      projectionVector(1:numberOfCoordinates)=dataPointLocation(1:numberOfCoordinates)- &
        & interpolatedPoint%values(1:numberOfCoordinates,NO_PART_DERIV)
      projectionDistance=DSQRT(DOT_PRODUCT(projectionVector(1:numberOfCoordinates),projectionVector(1:numberOfCoordinates)))
    ELSE
      projectionExitTag=DATA_PROJECTION_EXIT_TAG_NO_ELEMENT
      relativeTolerance=dataProjection%relativeTolerance
      absoluteTolerance=dataProjection%absoluteTolerance
      maximumDelta=dataProjection%maximumIterationUpdate
      minimumDelta=0.025_DP*maximumDelta !need to set a minimum, in case if it gets too small
      !Project on each candidate elements
      DO elementIdx=1,SIZE(candidateElements,1) 
        elementNumber=candidateElements(elementIdx)
        IF(elementNumber>0) THEN
          elementLineNumber=candidateElementLines(elementIdx)
          CALL DecompositionElements_ElementLineNumberGet(decompositionElements,elementLineNumber,elementNumber,lineNumber, &
            & err,error,*999)
          exitTag=DATA_PROJECTION_EXIT_TAG_NO_ELEMENT
          converged=.FALSE.
          !start at half the maximumDelta as we do not know if quadratic model is a good approximation yet
          delta=0.5_DP*maximumDelta
          CALL Field_InterpolationParametersLineGet(dataProjection%projectionSetType,lineNumber,interpolationParameters, &
            & err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
          xi=projectionXi
          CALL Field_InterpolateXi(SECOND_PART_DERIV,xi,interpolatedPoint,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
          distanceVector(1:numberOfCoordinates)=dataPointLocation(1:numberOfCoordinates)- &
            & interpolatedPoint%values(1:numberOfCoordinates,NO_PART_DERIV)
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
            functionGradient=-2.0_DP*(DOT_PRODUCT(distanceVector(1:numberOfCoordinates), &
              & interpolatedPoint%values(1:numberOfCoordinates,FIRST_PART_DERIV)))
            !functionHessian 
            functionHessian=-2.0_DP*(DOT_PRODUCT(distanceVector(1:numberOfCoordinates), &
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
              ELSE IF(newXi(1)>1.0_DP) THEN
                newXi(1)=1.0_DP  
              ENDIF
              CALL Field_InterpolateXi(SECOND_PART_DERIV,newXi,interpolatedPoint,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
              distanceVector(1:numberOfCoordinates)=dataPointLocation(1:numberOfCoordinates)- &
                & interpolatedPoint%values(1:numberOfCoordinates,NO_PART_DERIV)
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
                ELSE IF(predictionAccuracy>0.9_DP.AND.predictionAccuracy<1.1_DP) THEN
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
    ENDIF
    
    EXITS("DataProjection_NewtonLinesEvaluate")
    RETURN
999 ERRORSEXITS("DataProjection_NewtonLinesEvaluate",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_NewtonLinesEvaluate
  
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

    CALL DataProjection_AssertNotFinished(dataProjection,err,error,*999)
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
  SUBROUTINE DataProjection_ProjectionCancelByDataPoints0(dataProjection,dataPointGlobalNumber,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to cancel projections for
    INTEGER(INTG), INTENT(IN) :: dataPointGlobalNumber !<The data point global number to use to cancel projections 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("DataProjection_ProjectionCancelByDataPoints0",err,error,*999)

    CALL DataProjection_ProjectionCancelByDataPoints1(dataProjection,[dataPointGlobalNumber],err,error,*999)
    
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
  SUBROUTINE DataProjection_ProjectionCancelByDataPoints1(dataProjection,dataPointGlobalNumbers,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to cancel projections for
    INTEGER(INTG), INTENT(IN) :: dataPointGlobalNumbers(:) !<dataPointUserNumbers(dataPointIdx). The data point global numbers to use to cancel projections.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointIdx,dataPointGlobalNumber
    TYPE(VARYING_STRING) :: localError
   
    ENTERS("DataProjection_ProjectionCancelByDataPoints1",err,error,*999)

    CALL DataProjection_AssertIsFinished(dataProjection,err,error,*999)
    CALL DataProjection_AssertIsProjected(dataProjection,err,error,*999)
#ifdef WITH_PRECHECKS
     IF(.NOT.ALLOCATED(dataProjection%dataProjectionResults)) &
       CALL FlagError("Data projection projection results is not allocated.",err,error,*999)
#endif
    
    DO dataPointIdx=1,SIZE(dataPointGlobalNumbers,1)
      CALL DataProjection_DataPointGlobalNumberGet(dataProjection,dataPointUserNumbers(dataPointIdx),dataPointGlobalNumber, &
        & err,error,*999)
      CALL DataProjection_DataProjectionResultCancel(dataProjection%dataProjectionResults(dataPointGlobalNumber),err,error,*999)
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

    CALL DataProjection_AssertIsFinished(dataProjection,err,error,*999)
    CALL DataProjection_AssertIsProjected(dataProjection,err,error,*999)
#ifdef WITH_PRECHECKS
     IF(.NOT.ALLOCATED(dataProjection%dataProjectionResults)) &
       CALL FlagError("Data projection projection results is not allocated.",err,error,*999)
#endif
 
    NULLIFY(dataPoints)
    CALL DataProjection_DataPointsGet(dataProjection,dataPoints,err,error,*999)
    SELECT CASE(distanceRelation)
    CASE(DATA_PROJECTION_DISTANCE_GREATER)
      DO dataPointIdx=1,dataPoints%numberOfDataPoints
        IF(dataProjection%dataProjectionResults(dataPointIdx)%distance>distance) THEN
          CALL DataProjection_DataProjectionResultCancel(dataProjection%dataProjectionResults(dataPointIdx),err,error,*999)
        ENDIF
      ENDDO !dataPointIdx
    CASE(DATA_PROJECTION_DISTANCE_GREATER_EQUAL)
      DO dataPointIdx=1,dataPoints%numberOfDataPoints
        IF(dataProjection%dataProjectionResults(dataPointIdx)%distance>=distance) THEN
          CALL DataProjection_DataProjectionResultCancel(dataProjection%dataProjectionResults(dataPointIdx),err,error,*999)
        ENDIF
      ENDDO !dataPointIdx
    CASE(DATA_PROJECTION_DISTANCE_LESS)
      DO dataPointIdx=1,dataPoints%numberOfDataPoints
        IF(dataProjection%dataProjectionResults(dataPointIdx)%distance<distance) THEN
          CALL DataProjection_DataProjectionResultCancel(dataProjection%dataProjectionResults(dataPointIdx),err,error,*999)
        ENDIF
      ENDDO !dataPointIdx
    CASE(DATA_PROJECTION_DISTANCE_LESS_EQUAL)
      DO dataPointIdx=1,dataPoints%numberOfDataPoints
        IF(dataProjection%dataProjectionResults(dataPointIdx)%distance<=distance) THEN
          CALL DataProjection_DataProjectionResultCancel(dataProjection%dataProjectionResults(dataPointIdx),err,error,*999)
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

    CALL DataProjection_AssertIsFinished(dataProjection,err,error,*999)
    CALL DataProjection_AssertIsProjected(dataProjection,err,error,*999)
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
          dataProjection%dataProjectionResults(dataPointIdx)%globalElementNumber=0
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
  SUBROUTINE DataProjection_ProjectionCandidateElementsSet(dataProjection,elementLocalNumbers,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to set the projection type for
    INTEGER(INTG), INTENT(IN) :: elementLocalNumbers(:) !<elementLocalNumbers(elementIdx). The projection candidate local element numbers
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: elementIdx,elementLocalNumber,numberOfCandidates
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
      IF(ALLOCATED(dataProjection%dataProjectionCandidates(0)%candidateElementNumbers)) &
        & DEALLOCATE(dataProjection%dataProjectionCandidates(0)%candidateElementNumbers)
      ALLOCATE(dataProjection%dataProjectionCandidates(0)%candidateElementNumbers(SIZE(elementLocalNumbers,1)),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate candidiate element numbers.",err,error,*999)
      numberOfCandidates=0
      DO elementIdx=1,SIZE(elementLocalNumbers,1)
        elementLocalNumber=elementLocalNumbers(elementIdx)
        numberOfCandidates=numberOfCandidates+1
        dataProjection%dataProjectionCandidates(0)%candidateElementNumbers(numberOfCandidates)=elementLocalNumber
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
  SUBROUTINE DataProjection_ProjectionDataCandidateElementsSet(dataProjection,dataPointGlobalNumbers,elementLocalNumbers, &
    & err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to set the projection type for
    INTEGER(INTG), INTENT(IN) :: dataPointGlobalNumbers(:) !<dataPointGlobalNumbers(dataPointIdx). The data points for the projection candidate local element numbers
    INTEGER(INTG), INTENT(IN) :: elementLocalNumbers(:) !<elementLocalNumbers(elementIdx). The projection candidate local element numbers
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointGlobalNumber,dataPointIdx,elementIdx,elementLocalNumber,numberOfCandidates,numberOfDataPoints, &
      & numberOfElements,totalNumberOfElements
    TYPE(DataPointsType), POINTER :: dataPoints
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionElementsType), POINTER :: decompositionElements
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology
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
      NULLIFY(dataPoints)
      CALL DataProjection_DataPointsGet(dataProjection,dataPoints,err,error,*999)
      CALL DataPoints_NumberOfDataPointsGet(dataPoints,numberOfDataPoints,err,error,*999)
      NULLIFY(decomposition)
      CALL DataProjection_DecompositionGet(dataProjection,decomposition,err,error,*999)
      NULLIFY(decompositionTopology)
      CALL Decomposition_DecompositionTopologyGet(decomposition,decompositionTopology,err,error,*999)
      NULLIFY(decompositionElements)
      CALL DecompositionTopology_DecompositionElementsGet(decompositionTopology,decompositionElements,err,error,*999)
      CALL DecompositionElements_NumberOfElementsGet(decompositionElements,numberOfElements,err,error,*999)
      CALL DecompositionElements_TotalNumberOfElementsGet(decompositionElements,totalNumberOfElements,err,error,*999)
      DO dataPointIdx=1,SIZE(dataPointGlobalNumbers,1)
        dataPointGlobalNumber=dataPointGlobalNumbers(dataPointIdx)
        IF(dataPointGlobalNumber<1.OR.dataPointGlobalNumber>numberOfDataPoints) THEN
          localError="The data point global number of "//TRIM(NumberToVString(dataPointGlobalNumber,"*",err,error))// &
            & " at position "//TRIM(NumberToVString(dataPointIdx,"*",err,error))// &
            & " in the list of data points is invalid. The data point global number should be >=1 and <= "// &
            & TRIM(NumberToVString(numberOfDataPoints,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        IF(ALLOCATED(dataProjection%dataProjectionCandidates(dataPointGlobalNumber)%candidateElementNumbers)) &
          & DEALLOCATE(dataProjection%dataProjectionCandidates(dataPointGlobalNumber)%candidateElementNumbers)
        ALLOCATE(dataProjection%dataProjectionCandidates(dataPointGlobalNumber)%candidateElementNumbers( &
          & SIZE(elementLocalNumbers,1)),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate candidiate element numbers.",err,error,*999)
        numberOfCandidates=0
        DO elementIdx=1,SIZE(elementLocalNumbers,1)
          elementLocalNumber=elementLocalNumbers(elementIdx)
          IF(elementLocalNumber<1.OR.elementLocalNumber>numberOfElements) THEN
            localError="The specified local error number of "//TRIM(NumberToVString(elementLocalNumber,"*",err,error))// &
              & " at position "//TRIM(NumberToVString(elementIdx,"*",err,error))// &
              & " is invalid. The local element number should be >= 1 and <= "// &
              & TRIM(NumberToVString(numberOfElements,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
          numberOfCandidates=numberOfCandidates+1
          dataProjection%dataProjectionCandidates(dataPointGlobalNumber)%candidateElementNumbers(numberOfCandidates)= &
            & elementLocalNumber
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
  SUBROUTINE DataProjection_ProjectionCandidateFacesSet(dataProjection,elementLocalNumbers,localFaceNormals,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to set the projection type for
    INTEGER(INTG), INTENT(IN) :: elementLocalNumbers(:) !<elementLocalNumbers(elementIdx). The projection candidate local element numbers
    INTEGER(INTG), INTENT(IN) :: localFaceNormals(:) !<localFaceNormals(elementIdx). The projection candidate element face normals
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: elementIdx,elementLocalNumber,localFaceNumber,numberOfCandidates,numberOfElements,totalNumberOfElements
    TYPE(BasisType), POINTER :: basis
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionElementsType), POINTER :: decompositionElements
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainElementsType), POINTER :: domainElements
    TYPE(DomainTopologyType), POINTER :: domainTopology
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
      IF(SIZE(elementLocalNumbers,1)/=SIZE(localFaceNormals,1)) THEN
        localError="The size of the element local numbers array of "// &
          & TRIM(NumberToVString(SIZE(elementLocalNumbers,1),"*",err,error))// &
          & " does not match the size of the face normals array of "// &
          & TRIM(NumberToVString(SIZE(localFaceNormals,1),"*",err,error))//"."
        CALL FlagError(localError,err,error,*998)
      ENDIF
      NULLIFY(decomposition)
      CALL DataProjection_DecompositionGet(dataProjection,decomposition,err,error,*998)
      NULLIFY(decompositionTopology)
      CALL Decomposition_DecompositionTopologyGet(decomposition,decompositionTopology,err,error,*998)
      NULLIFY(decompositionElements)
      CALL DecompositionTopology_DecompositionElementsGet(decompositionTopology,decompositionElements,err,error,*998)
      CALL DecompositionElements_NumberOfElementsGet(decompositionElements,numberOfElements,err,error,*998)
      CALL DecompositionElements_TotalNumberOfElementsGet(decompositionElements,totalNumberOfElements,err,error,*998)
      NULLIFY(domain)
      CALL Decomposition_DomainGet(decomposition,0,domain,err,error,*998)
      NULLIFY(domainTopology)
      CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*998)
      NULLIFY(domainElements)
      CALL DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*998)
      IF(ALLOCATED(dataProjection%dataProjectionCandidates(0)%candidateElementNumbers)) &
        & DEALLOCATE(dataProjection%dataProjectionCandidates(0)%candidateElementNumbers)
      IF(ALLOCATED(dataProjection%dataProjectionCandidates(0)%localFaceLineNumbers)) &
        & DEALLOCATE(dataProjection%dataProjectionCandidates(0)%localFaceLineNumbers)
      ALLOCATE(dataProjection%dataProjectionCandidates(0)%candidateElementNumbers(SIZE(elementLocalNumbers,1)),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate candidiate element numbers.",err,error,*999)
      ALLOCATE(dataProjection%dataProjectionCandidates(0)%localFaceLineNumbers(SIZE(localFaceNormals,1)),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate candidiate local face/line numbers.",err,error,*999)
      numberOfCandidates=0
      DO elementIdx=1,SIZE(elementLocalNumbers,1)
        elementLocalNumber=elementLocalNumbers(elementIdx)
        IF(elementLocalNumber<1.OR.elementLocalNumber>numberOfElements) THEN
          localError="The specified local error number of "//TRIM(NumberToVString(elementLocalNumber,"*",err,error))// &
            & " at position "//TRIM(NumberToVString(elementIdx,"*",err,error))// &
            & " is invalid. The local element number should be >= 1 and <= "// &
            & TRIM(NumberToVString(numberOfElements,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        numberOfCandidates=numberOfCandidates+1           
        NULLIFY(basis)
        CALL DomainElements_ElementBasisGet(domainElements,elementLocalNumber,basis,err,error,*999)
        CALL Basis_LocalFaceNumberGet(basis,localFaceNormals(elementIdx),localFaceNumber,err,error,*999)
        dataProjection%dataProjectionCandidates(0)%candidateElementNumbers(numberOfCandidates)=elementLocalNumber
        dataProjection%dataProjectionCandidates(0)%localFaceLineNumbers(numberOfCandidates)=localFaceNumber
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
  SUBROUTINE DataProjection_ProjectionDataCandidateFacesSet(dataProjection,dataPointGlobalNumbers,elementLocalNumbers, &
    & localFaceNormals,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to set the projection type for
    INTEGER(INTG), INTENT(IN) :: dataPointGlobalNumbers(:) !<dataPointGlobalNumbers(dataPointIdx). The data points for the projection candidate local element numbers
    INTEGER(INTG), INTENT(IN) :: elementLocalNumbers(:) !<elementLocalNumbers(elementIdx). The projection candidate local element numbers
    INTEGER(INTG), INTENT(IN) :: localFaceNormals(:) !<localFaceNormals(elementIdx). The projection candidate element face normals
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointGlobalNumber,dataPointIdx,elementIdx,elementLocalNumber,localFaceNumber,numberOfCandidates, &
      & numberOfDataPoints,numberOfElements,totalNumberOfElements
    TYPE(BasisType), POINTER :: basis
    TYPE(DataPointsType), POINTER :: dataPoints
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionElementsType), POINTER :: decompositionElements
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainElementsType), POINTER :: domainElements
    TYPE(DomainTopologyType), POINTER :: domainTopology
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
      IF(SIZE(elementLocalNumbers,1)/=SIZE(localFaceNormals,1)) THEN
        localError="The size of the element local numbers array of "// &
          & TRIM(NumberToVString(SIZE(elementLocalNumbers,1),"*",err,error))// &
          & " does not match the size of the face normals array of "// &
          & TRIM(NumberToVString(SIZE(localFaceNormals,1),"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      NULLIFY(dataPoints)
      CALL DataProjection_DataPointsGet(dataProjection,dataPoints,err,error,*999)
      CALL DataPoints_NumberOfDataPointsGet(dataPoints,numberOfDataPoints,err,error,*999)
      NULLIFY(decomposition)
      CALL DataProjection_DecompositionGet(dataProjection,decomposition,err,error,*999)
      NULLIFY(decompositionTopology)
      CALL Decomposition_DecompositionTopologyGet(decomposition,decompositionTopology,err,error,*999)
      NULLIFY(decompositionElements)
      CALL DecompositionTopology_DecompositionElementsGet(decompositionTopology,decompositionElements,err,error,*999)
      CALL DecompositionElements_NumberOfElementsGet(decompositionElements,numberOfElements,err,error,*999)
      CALL DecompositionElements_TotalNumberOfElementsGet(decompositionElements,totalNumberOfElements,err,error,*999)
      NULLIFY(domain)
      CALL Decomposition_DomainGet(decomposition,0,domain,err,error,*999)
      NULLIFY(domainTopology)
      CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
      NULLIFY(domainElements)
      CALL DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*999)
      DO dataPointIdx=1,SIZE(dataPointGlobalNumbers,1)
        dataPointGlobalNumber=dataPointGlobalNumbers(dataPointIdx)
        IF(dataPointGlobalNumber<1.OR.dataPointGlobalNumber>numberOfDataPoints) THEN
          localError="The data point global number of "//TRIM(NumberToVString(dataPointGlobalNumber,"*",err,error))// &
            & " at position "//TRIM(NumberToVString(dataPointIdx,"*",err,error))// &
            & " in the list of data points is invalid. The data point global number should be >=1 and <= "// &
            & TRIM(NumberToVString(numberOfDataPoints,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        IF(ALLOCATED(dataProjection%dataProjectionCandidates(dataPointGlobalNumber)%candidateElementNumbers)) &
          & DEALLOCATE(dataProjection%dataProjectionCandidates(dataPointGlobalNumber)%candidateElementNumbers)
        IF(ALLOCATED(dataProjection%dataProjectionCandidates(dataPointGlobalNumber)%localFaceLineNumbers)) &
          & DEALLOCATE(dataProjection%dataProjectionCandidates(dataPointGlobalNumber)%localFaceLineNumbers)
        ALLOCATE(dataProjection%dataProjectionCandidates(dataPointGlobalNumber)%candidateElementNumbers( &
          & SIZE(elementLocalNumbers,1)),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate candidiate element numbers.",err,error,*999)
        ALLOCATE(dataProjection%dataProjectionCandidates(dataPointGlobalNumber)%localFaceLineNumbers( &
          & SIZE(localFaceNormals,1)),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate candidiate local face/line numbers.",err,error,*999)
        numberOfCandidates=0
        DO elementIdx=1,SIZE(elementLocalNumbers,1)
          elementLocalNumber=elementLocalNumbers(elementIdx)
          IF(elementLocalNumber<1.OR.elementLocalNumber>numberOfElements) THEN
            localError="The specified local error number of "//TRIM(NumberToVString(elementLocalNumber,"*",err,error))// &
              & " at position "//TRIM(NumberToVString(elementIdx,"*",err,error))// &
              & " is invalid. The local element number should be >= 1 and <= "// &
              & TRIM(NumberToVString(numberOfElements,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
          numberOfCandidates=numberOfCandidates+1           
          NULLIFY(basis)
          CALL DomainElements_ElementBasisGet(domainElements,elementLocalNumber,basis,err,error,*999)
          CALL Basis_LocalFaceNumberGet(basis,localFaceNormals(elementIdx),localFaceNumber,err,error,*999)
          dataProjection%dataProjectionCandidates(dataPointGlobalNumber)%candidateElementNumbers(numberOfCandidates)= &
            & elementLocalNumber
          dataProjection%dataProjectionCandidates(dataPointGlobalNumber)%localFaceLineNumbers(numberOfCandidates)= &
            & localFaceNumber
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
  SUBROUTINE DataProjection_ProjectionCandidateLinesSet(dataProjection,elementLocalNumbers,localLineNormals,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to set the projection type for
    INTEGER(INTG), INTENT(IN) :: elementLocalNumbers(:) !<elementLocalNumbers(elementIdx). The projection candidate local element numbers
    INTEGER(INTG), INTENT(IN) :: localLineNormals(:,:) !<localLineNormals(normalIdx,elementIdx). The projection candidate line xi normals
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: elementIdx,elementLocalNumber,localLineNumber,numberOfCandidates,numberOfDimensions,numberOfElements, &
      & numberOfXi,totalNumberOfElements
    TYPE(BasisType), POINTER :: basis
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionElementsType), POINTER :: decompositionElements
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainElementsType), POINTER :: domainElements
    TYPE(DomainTopologyType), POINTER :: domainTopology
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
      IF(SIZE(elementLocalNumbers,1)/=SIZE(localLineNormals,2)) THEN
        localError="The size of the element user numbers array of "// &
          & TRIM(NumberToVString(SIZE(elementLocalNumbers,1),"*",err,error))// &
          & " does not match the second size of the line normals array of "// &
          & TRIM(NumberToVString(SIZE(localLineNormals,2),"*",err,error))//"."
        CALL FlagError(localError,err,error,*998)
      ENDIF
      NULLIFY(decomposition)
      CALL DataProjection_DecompositionGet(dataProjection,decomposition,err,error,*998)
      CALL Decomposition_NumberOfDimensionsGet(decomposition,numberOfDimensions,err,error,*998)
      IF(SIZE(localLineNormals,1)<(numberOfDimensions-1)) THEN
        localError="The first size of the line normals array of "// &
          & TRIM(NumberToVString(SIZE(localLineNormals,1),"*",err,error))//" is invalid. The size should be >= "// &
          & TRIM(NumberToVString(numberOfDimensions-1,"*",err,error))//"."
        CALL FlagError(localError,err,error,*998)
      ENDIF
      NULLIFY(decompositionTopology)
      CALL Decomposition_DecompositionTopologyGet(decomposition,decompositionTopology,err,error,*999)
      NULLIFY(decompositionElements)
      CALL DecompositionTopology_DecompositionElementsGet(decompositionTopology,decompositionElements,err,error,*999)
      CALL DecompositionElements_NumberOfElementsGet(decompositionElements,numberOfElements,err,error,*999)
      CALL DecompositionElements_TotalNumberOfElementsGet(decompositionElements,totalNumberOfElements,err,error,*999)
      NULLIFY(domain)
      CALL Decomposition_DomainGet(decomposition,0,domain,err,error,*998)
      NULLIFY(domainTopology)
      CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*998)
      NULLIFY(domainElements)
      CALL DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*998)
      IF(ALLOCATED(dataProjection%dataProjectionCandidates(0)%candidateElementNumbers)) &
        & DEALLOCATE(dataProjection%dataProjectionCandidates(0)%candidateElementNumbers)
      IF(ALLOCATED(dataProjection%dataProjectionCandidates(0)%localFaceLineNumbers)) &
        & DEALLOCATE(dataProjection%dataProjectionCandidates(0)%localFaceLineNumbers)
      ALLOCATE(dataProjection%dataProjectionCandidates(0)%candidateElementNumbers(SIZE(elementLocalNumbers,1)),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate candidiate element numbers.",err,error,*999)
      ALLOCATE(dataProjection%dataProjectionCandidates(0)%localFaceLineNumbers(SIZE(localLineNormals,1)),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate candidiate local face/line numbers.",err,error,*999)
      numberOfCandidates=0
      DO elementIdx=1,SIZE(elementLocalNumbers,1)
        elementLocalNumber=elementLocalNumbers(elementIdx)
        IF(elementLocalNumber<1.OR.elementLocalNumber>numberOfElements) THEN
          localError="The specified local error number of "//TRIM(NumberToVString(elementLocalNumber,"*",err,error))// &
            & " at position "//TRIM(NumberToVString(elementIdx,"*",err,error))// &
            & " is invalid. The local element number should be >= 1 and <= "// &
            & TRIM(NumberToVString(numberOfElements,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        numberOfCandidates=numberOfCandidates+1
        NULLIFY(basis)
        CALL DomainElements_ElementBasisGet(domainElements,elementLocalNumber,basis,err,error,*999)
        CALL Basis_NumberOfXiGet(basis,numberOfXi,err,error,*999)
        CALL Basis_LocalLineNumberGet(basis,localLineNormals(1:(numberOfXi-1),elementIdx),localLineNumber,err,error,*999)
        dataProjection%dataProjectionCandidates(0)%candidateElementNumbers(elementIdx)=elementLocalNumber
        dataProjection%dataProjectionCandidates(0)%localFaceLineNumbers(elementIdx)=localLineNumber
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
  SUBROUTINE DataProjection_ProjectionDataCandidateLinesSet(dataProjection,dataPointGlobalNumbers,elementLocalNumbers, &
    & localLineNormals,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to set the projection type for
    INTEGER(INTG), INTENT(IN) :: dataPointGlobalNumbers(:) !<dataPointGlobalNumbers(dataPointIdx). The data points for the projection candidate local element numbers
    INTEGER(INTG), INTENT(IN) :: elementLocalNumbers(:) !<elementLocalNumbers(elementIdx). The projection candidate local element numbers
    INTEGER(INTG), INTENT(IN) :: localLineNormals(:,:) !<localLineNormals(normalIdx,elementIdx). The projection candidate line xi normals
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointGlobalNumber,dataPointIdx,elementIdx,elementLocalNumber,localLineNumber,numberOfCandidates, &
      & numberOfDataPoints,numberOfDimensions,numberOfElements,numberOfXi,totalNumberOfElements
    TYPE(BasisType), POINTER :: basis
    TYPE(DataPointsType), POINTER :: dataPoints
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionElementsType), POINTER :: decompositionElements
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainElementsType), POINTER :: domainElements
    TYPE(DomainTopologyType), POINTER :: domainTopology
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
      IF(SIZE(elementLocalNumbers,1)/=SIZE(localLineNormals,2)) THEN
        localError="The size of the element user numbers array of "// &
          & TRIM(NumberToVString(SIZE(elementLocalNumbers,1),"*",err,error))// &
          & " does not match the second size of the line normals array of "// &
          & TRIM(NumberToVString(SIZE(localLineNormals,2),"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      NULLIFY(dataPoints)
      CALL DataProjection_DataPointsGet(dataProjection,dataPoints,err,error,*999)
      CALL DataPoints_NumberOfDataPointsGet(dataPoints,numberOfDataPoints,err,error,*999)
      NULLIFY(decomposition)
      CALL DataProjection_DecompositionGet(dataProjection,decomposition,err,error,*999)
      CALL Decomposition_NumberOfDimensionsGet(decomposition,numberOfDimensions,err,error,*999)
      IF(SIZE(localLineNormals,1)<(numberOfDimensions-1)) THEN
        localError="The first size of the line normals array of "// &
          & TRIM(NumberToVString(SIZE(localLineNormals,1),"*",err,error))//" is invalid. The size should be >= "// &
          & TRIM(NumberToVString(numberOfDimensions-1,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      NULLIFY(decompositionTopology)
      CALL Decomposition_DecompositionTopologyGet(decomposition,decompositionTopology,err,error,*999)
      NULLIFY(decompositionElements)
      CALL DecompositionTopology_DecompositionElementsGet(decompositionTopology,decompositionElements,err,error,*999)
      CALL DecompositionElements_NumberOfElementsGet(decompositionElements,numberOfElements,err,error,*999)
      CALL DecompositionElements_TotalNumberOfElementsGet(decompositionElements,totalNumberOfElements,err,error,*999)
      NULLIFY(domain)
      CALL Decomposition_DomainGet(decomposition,0,domain,err,error,*999)
      NULLIFY(domainTopology)
      CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
      NULLIFY(domainElements)
      CALL DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*999)
      DO dataPointIdx=1,SIZE(dataPointGlobalNumbers,1)
        dataPointGlobalNumber=dataPointGlobalNumbers(dataPointIdx)
        IF(dataPointGlobalNumber<1.OR.dataPointGlobalNumber>numberOfDataPoints) THEN
          localError="The data point global number of "//TRIM(NumberToVString(dataPointGlobalNumber,"*",err,error))// &
            & " at position "//TRIM(NumberToVString(dataPointIdx,"*",err,error))// &
            & " in the list of data points is invalid. The data point global number should be >=1 and <= "// &
            & TRIM(NumberToVString(numberOfDataPoints,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        IF(ALLOCATED(dataProjection%dataProjectionCandidates(dataPointGlobalNumber)%candidateElementNumbers)) &
          & DEALLOCATE(dataProjection%dataProjectionCandidates(dataPointGlobalNumber)%candidateElementNumbers)
        IF(ALLOCATED(dataProjection%dataProjectionCandidates(dataPointGlobalNumber)%localFaceLineNumbers)) &
          & DEALLOCATE(dataProjection%dataProjectionCandidates(dataPointGlobalNumber)%localFaceLineNumbers)
        ALLOCATE(dataProjection%dataProjectionCandidates(dataPointGlobalNumber)%candidateElementNumbers( &
          & SIZE(elementLocalNumbers,1)),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate candidiate element numbers.",err,error,*999)
        ALLOCATE(dataProjection%dataProjectionCandidates(dataPointGlobalNumber)%localFaceLineNumbers( &
          & SIZE(localLineNormals,1)),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate candidiate local face/line numbers.",err,error,*999)
        numberOfCandidates=0
        DO elementIdx=1,SIZE(elementLocalNumbers,1)
          elementLocalNumber=elementLocalNumbers(elementIdx)
          IF(elementLocalNumber<1.OR.elementLocalNumber>numberOfElements) THEN
            localError="The specified local error number of "//TRIM(NumberToVString(elementLocalNumber,"*",err,error))// &
              & " at position "//TRIM(NumberToVString(elementIdx,"*",err,error))// &
              & " is invalid. The local element number should be >= 1 and <= "// &
              & TRIM(NumberToVString(numberOfElements,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
          numberOfCandidates=numberOfCandidates+1
          NULLIFY(basis)
          CALL DomainElements_ElementBasisGet(domainElements,elementLocalNumber,basis,err,error,*999)
          CALL Basis_NumberOfXiGet(basis,numberOfXi,err,error,*999)
          CALL Basis_LocalLineNumberGet(basis,localLineNormals(1:(numberOfXi-1),elementIdx),localLineNumber,err,error,*999)
          dataProjection%dataProjectionCandidates(dataPointGlobalNumber)%candidateElementNumbers(elementIdx)=elementLocalNumber
          dataProjection%dataProjectionCandidates(dataPointGlobalNumber)%localFaceLineNumbers(elementIdx)=localLineNumber
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

  !>Sets the projection type for a data projection.
  SUBROUTINE DataProjection_ProjectionTypeSet(dataProjection,projectionType,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to set the projection type for
    INTEGER(INTG), INTENT(IN) :: projectionType !<the projection type to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    REAL(DP), ALLOCATABLE :: startingXi(:,:)
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DataProjection_ProjectionTypeSet",err,error,*999)

    CALL DataProjection_AssertNotFinished(dataProjection,err,error,*999)
    
    dataProjection%projectionType=projectionType
    SELECT CASE(projectionType)
    CASE (DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE)
      dataProjection%numberOfXi=1
    CASE (DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE)
      dataProjection%numberOfXi=2
    CASE (DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE)
      NULLIFY(decomposition)
      CALL DataProjection_DecompositionGet(dataProjection,decomposition,err,error,*999)
      CALL Decomposition_NumberOfDimensionsGet(decomposition,dataProjection%numberOfXi,err,error,*999)
    CASE DEFAULT
      localError="The specified projection type of "//TRIM(NumberToVString(projectionType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    IF(dataProjection%numberOfXi/=SIZE(dataProjection%startingXi,1)) THEN
      ALLOCATE(startingXi(dataProjection%numberOfXi,dataProjection%datapoints%numberOfDatapoints),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate starting xi.",err,error,*999)
      IF(dataProjection%numberOfXi>SIZE(dataProjection%startingXi,2)) THEN
        startingXi(1:SIZE(dataProjection%startingXi,1),:)=dataProjection%startingXi(1:SIZE(dataProjection%startingXi,1),:)
        startingXi(SIZE(dataProjection%startingXi,1):dataProjection%numberOfXi,:)=0.5_DP
      ELSE
        startingXi(1:SIZE(dataProjection%startingXi,1),:)=dataProjection%startingXi(1:dataProjection%numberOfXi,:)
      ENDIF
      CALL MOVE_ALLOC(startingXi,dataProjection%startingXi)
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

  !>Sets the relative tolerance for a data projection.
  SUBROUTINE DataProjection_RelativeToleranceSet(dataProjection,relativeTolerance,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to set the relative tolerance for
    REAL(DP), INTENT(IN) :: relativeTolerance !<the relative tolerance to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DataProjection_RelativeToleranceSet",err,error,*999)

    CALL DataProjection_AssertNotFinished(dataProjection,err,error,*999)
    IF(relativeTolerance<0.0_DP) THEN
      localError="The specified relative tolerance of "//TRIM(NumberToVString(relativeTolerance,"*",err,error))// &
        & " is invalid. The tolerance must be >= 0.0."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    dataProjection%relativeTolerance=relativeTolerance
    
    EXITS("DataProjection_RelativeToleranceSet")
    RETURN
999 ERRORSEXITS("DataProjection_RelativeToleranceSet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_RelativeToleranceSet

  !
  !================================================================================================================================
  !

!  !>Sets the starting xi for a data projection for a data point specified by a global number
  SUBROUTINE DataProjection_StartingXiSet(dataProjection,dataPointGlobalNumber,startingXi,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to set the starting xi for
    INTEGER(INTG) :: dataPointGlobalNumber !<The data point global number to set the starting xi for
    REAL(DP), INTENT(IN) :: startingXi(:) !<The starting xi to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("DataProjection_StartingXiSet",err,error,*999)

    CALL DataProjection_AssertNotFinished(dataProjection,err,error,*999)
    
#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(dataProjection%dataPoints)) CALL FlagError("Data projection data points is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(dataProject%dataPoints)) CALL FlagError("Data projection data points is not associated.",err,error,*999)
    IF(dataPointGlobalNumber<1.OR.dataPointGlobalNumber>dataProjection%dataPoints%numberOfDataPoints) THEN
      localError="The specified data point global number of "//TRIM(NumberToVString(dataPointGlobalNumber,"*",err,error))// &
        & " is invalid. The global number should be >= 1 and <= "// &
        & TRIM(NumberToVString(dataProjection%dataPoints%numberOfDataPoints,"*",err,error))//"."
    ENDIF
    IF(SIZE(startingXi,1)<dataProjection%numberOfXi) THEN
      localError="The size of the the specified xi array of "// &
        & TRIM(NumberToVString(SIZE(startingXi,1),"*",err,error))//" is invalid. The size should be >= "// &
        & TRIM(NumberToVString(dataProjection%numberOfXi,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF((ANY(startingXi(1:dataProjection%numberOfXi)<0.0_DP)).OR.(ANY(startingXi(1:dataProjection%numberOfXi)>1.0_DP))) &
      & CALL FlagError("The specified starting xi location must be between [0.0,1.0] for all xi coordinates.",err,error,*999)
#endif
   
    dataProjection%startingXi(1:dataProjection%numberOfXi,dataPointGlobalNumber)=startingXi(1:dataProjection%numberOfXi)
    
    EXITS("DataProjection_StartingXiSet")
    RETURN
999 ERRORSEXITS("DataProjection_StartingXiSet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_StartingXiSet
  
  !
  !================================================================================================================================
  !

  !>Sets the element for a data projection.
  SUBROUTINE DataProjection_ElementSet(dataProjection,dataPointGlobalNumber,elementLocalNumber,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to set the element for
    INTEGER(INTG), INTENT(IN) :: dataPointGlobalNumber !<data point user number
    INTEGER(INTG), INTENT(IN) :: elementLocalNumber !<the user element number to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfDataPoints,numberOfElements,totalNumberOfElements
    TYPE(DataPointsType), POINTER :: dataPoints
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionElementsType), POINTER :: decompositionElements
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DataProjection_ElementSet",err,error,*999)

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    NULLIFY(dataPoints)
    CALL DataProjection_DataPointsGet(dataProjection,dataPoints,err,error,*999)
    CALL DataPoints_NumberOfDataPointsGet(dataPoints,numberOfDataPoints,err,error,*999)
    NULLIFY(decomposition)
    CALL DataProjection_DecompositionGet(dataProjection,decomposition,err,error,*999)
    NULLIFY(decompositionTopology)
    CALL Decomposition_DecompositionTopologyGet(decomposition,decompositionTopology,err,error,*999)
    NULLIFY(decompositionElements)
    CALL DecompositionTopology_DecompositionElementsGet(decompositionTopology,decompositionElements,err,error,*999)
    CALL DecompositionElements_NumberOfElementsGet(decompositionElements,numberOfElements,err,error,*999)
    CALL DecompositionElements_TotalNumberOfElementsGet(decompositionElements,totalNumberOfElements,err,error,*999)
    
    IF(dataPointGlobalNumber<1.OR.dataPointGlobalNumber>numberOfDataPoints) THEN
      localError="The data point global number of "//TRIM(NumberToVString(dataPointGlobalNumber,"*",err,error))// &
        & " is invalid. The data point global number should be >=1 and <= "// &
        & TRIM(NumberToVString(numberOfDataPoints,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(elementLocalNumber<1.OR.elementLocalNumber>numberOfElements) THEN
      localError="The specified local error number of "//TRIM(NumberToVString(elementLocalNumber,"*",err,error))// &
        & " is invalid. The local element number should be >= 1 and <= "// &
        & TRIM(NumberToVString(numberOfElements,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    dataProjection%dataProjectionResults(dataPointGlobalNumber)%elementLocalNumber=elementLocalNumber
   
    EXITS("DataProjection_ElementSet")
    RETURN
999 ERRORSEXITS("DataProjection_ElementSet",err,error)
    RETURN 1

  END SUBROUTINE DataProjection_ElementSet

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
      & localLineFaceNumber,localLineXiNormals(2),myGroupComputationNodeNumber,normalIdx1,normalIdx2, &
      & numberOfCancelledDataPoints,numberOfDataPoints,numberOfGroupComputationNodes,outputID
    REAL(DP) :: distance
    CHARACTER(LEN=MAXSTRLEN) :: analFilename,format1,format2,localString
    TYPE(BasisType), POINTER :: basis
    TYPE(DataPointsType), POINTER :: dataPoints
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionElementsType), POINTER :: decompositionElements
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainElementsType), POINTER :: domainElements
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(FieldType), POINTER :: projectionField
    TYPE(VARYING_STRING) :: localError
    TYPE(WorkGroupType), POINTER :: workGroup
        
    ENTERS("DataProjection_ResultAnalysisOutput",err,error,*999)

    CALL DataProjection_AssertIsFinished(dataProjection,err,error,*999)
    CALL DataProjection_AssertIsProjected(dataProjection,err,error,*999)
    
    NULLIFY(dataPoints)
    CALL DataProjection_DataPointsGet(dataProjection,dataPoints,err,error,*999)
    CALL DataPoints_NumberOfDataPointsGet(dataPoints,numberOfDataPoints,err,error,*999)
    NULLIFY(projectionField)
    CALL DataProjection_ProjectionFieldGet(dataProjection,projectionField,err,error,*999)
    NULLIFY(decomposition)
    CALL Field_DecompositionGet(projectionField,decomposition,err,error,*999)
    NULLIFY(decompositionTopology)
    CALL Decomposition_DecompositionTopologyGet(decomposition,decompositionTopology,err,error,*999)
    NULLIFY(decompositionElements)
    CALL DecompositionTopology_DecompositionElementsGet(decompositionTopology,decompositionElements,err,error,*999)
    NULLIFY(domain)
    CALL Decomposition_DomainGet(decomposition,0,domain,err,error,*999)
    NULLIFY(domainTopology)
    CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
    NULLIFY(domainElements)
    CALL DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*999)
    NULLIFY(workGroup)
    CALL Decomposition_WorkGroupGet(decomposition,workGroup,err,error,*999)
    CALL WorkGroup_NumberOfGroupNodesGet(workGroup,numberOfGroupComputationNodes,err,error,*999)
    CALL WorkGroup_GroupNodeNumberGet(workGroup,myGroupComputationNodeNumber,err,error,*999)
    !Find the correct output ID and open a file if necessary
    filenameLength=LEN_TRIM(filename)
    IF(filenameLength>=1) THEN
      IF(numberOfGroupComputationNodes>1) THEN
        WRITE(analFilename,('(A,A,I0)')) filename(1:filenameLength),".opdataproj.",myGroupComputationNodeNumber              
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
    CALL WriteStringValue(outputID,"  Projection field number = ",dataProjection%projectionField%userNumber,err,error,*999)
    CALL WriteString(outputID,"",err,error,*999)
    CALL WriteStringValue(outputID,"  Number of data points = ",numberOfDataPoints,err,error,*999)
    CALL WriteString(outputID,"",err,error,*999)
    IF(numberOfDataPoints>0) THEN
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
      numberOfCancelledDataPoints=0
      DO dataPointIdx=1,numberOfDataPoints
        CALL DataPoints_DataUserNumberGet(dataPoints,dataPointIdx,dataPointUserNumber,err,error,*999)
        dataPointExitTag=dataProjection%dataProjectionResults(dataPointIdx)%exitTag
        IF(dataPointExitTag==DATA_PROJECTION_CANCELLED) THEN
          numberOfCancelledDataPoints=numberOfCancelledDataPoints+1
          WRITE(localString,'(2X,I8,2X,I8,8X,I2)') dataPointIdx,dataPointUserNumber,dataPointExitTag
          CALL WriteString(outputID,localString,err,error,*999)
        ELSE
          localElementNumber=dataProjection%dataProjectionResults(dataPointIdx)%elementLocalNumber
          elementUserNumber=decompositionElements%elements(localElementNumber)%userNumber
          localLineFaceNumber=dataProjection%dataProjectionResults(dataPointIdx)%elementLineFaceNumber
          distance=dataProjection%dataProjectionResults(dataPointIdx)%distance
          CALL DecompositionElements_ElementUserNumberGet(decompositionElements,localElementNumber,elementUserNumber, &
            & err,error,*999)
          NULLIFY(basis)
          CALL DomainElements_ElementBasisGet(domainElements,localElementNumber,basis,err,error,*999)
          SELECT CASE(dataProjection%projectionType)
          CASE(DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE)
            CALL Basis_LineXiNormalsGet(basis,localLineFaceNumber,localLineXiNormals,err,error,*999)
            normalIdx1=localLineXiNormals(1)
            normalIdx2=localLineXiNormals(2)
            WRITE(localString,'(2X,I8,2X,I8,8X,I2,2X,I8,9X,I2,X,I2,2X,E12.5,2X,E12.5,2X,E12.5)') dataPointIdx, &
              & dataPointUserNumber,dataPointExitTag,elementUserNumber,normalIdx1,normalIdx2, &
              & dataProjection%dataProjectionResults(dataPointIdx)%elementXi(1), &
              & dataProjection%dataProjectionResults(dataPointIdx)%projectionVector(1),distance
            format1="(56X,E12.5,2X,E12.5)"
            format2="(70X,E12.5)"
          CASE(DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE)
            CALL Basis_FaceXiNormalGet(basis,localLineFaceNumber,normalIdx1,err,error,*999)
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
      ENDDO !dataPointIdx
      CALL WriteString(outputID,"",err,error,*999)
      CALL WriteString(outputID,"  Number of cancelled data points = ",numberOfCancelledDataPoints,err,error,*999)
      CALL WriteString(outputID,"",err,error,*999)
      CALL WriteString(outputID,"  Errors:",err,error,*999)
      CALL WriteString(outputID,"",err,error,*999)
      CALL WriteString(outputID,"  Error type            Value  Data user#",err,error,*999)
      WRITE(localString,'("  RMS error    ",2X,E12.5)') dataProjection%rmsError
      CALL WriteString(outputID,localString,err,error,*999)
      IF(dataProjection%maximumErrorDataPoint/=0) THEN
        WRITE(localString,'("  Maximum error",2X,E12.5,4X,I8)') dataProjection%maximumError, &
          & dataProjection%dataProjectionResults(dataProjection%maximumErrorDataPoint)%userDataPointNumber
        CALL WriteString(outputID,localString,err,error,*999)
      ENDIF
      IF(dataProjection%minimumErrorDataPoint/=0) THEN
        WRITE(localString,'("  Minimum error",2X,E12.5,4X,I8)') dataProjection%minimumError, &
          & dataProjection%dataProjectionResults(dataProjection%minimumErrorDataPoint)%userDataPointNumber
        CALL WriteString(outputID,localString,err,error,*999)
      ENDIF
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

  !>Sets the projection element number for a data point identified by a given global number.
  SUBROUTINE DataProjection_ResultElementNumberSet(dataProjection,dataPointGlobalNumber,projectionElementUserNumber,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection for which projection result is stored
    INTEGER(INTG), INTENT(IN) :: dataPointGlobalNumber !<The data projection global number to set the projection element number for
    INTEGER(INTG), INTENT(IN) :: projectionElementUserNumber !<The projection element user number of the specified global data point to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointGlobalNumber,elementLocalNumber
    LOGICAL :: elementExists,ghostElement
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionElementsType), POINTER :: decompositionElements
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology
    TYPE(VARYING_STRING) :: localError
   
    ENTERS("DataProjection_ResultElementNumberSet",err,error,*999)

    CALL DataProjection_AssertIsFinished(dataProjection,err,error,*999)
    
    SELECT CASE(dataProjection%projectionType)
    CASE(DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE)
      CALL FlagError("Use the result element faces routine not the element routine for a boundary faces projection.", &
        & err,error,*999)
    CASE(DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE)
      CALL FlagError("Use the result element lines routine not the element routine for a boundary lines projection.", &
        & err,error,*999)
    CASE(DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE)
      NULLIFY(decomposition)
      CALL DataProjection_DecompositionGet(dataProjection,decomposition,err,error,*999)
      NULLIFY(decompositionTopology)
      CALL Decomposition_TopologyGet(decomposition,decompositionTopology,err,error,*999)
      NULLIFY(decompositionElements)
      CALL DecompositionTopology_DecompositionElementsGet(decompositionTopology,decompositionElements,err,error,*999)      
      CALL DecompositionElements_ElementCheckExists(decompositionElements,projectionElementUserNumber,elementExists, &
        & elementLocalNumber,ghostElement,err,error,*999)       
      IF(elementExists) THEN
        IF(.NOT.ghostElement) THEN
          CALL DecompositionElements_ElementGlobalNumberGet(decompositionElements,elementLocalNumber,elementGlobalNumber, &
            & err,error,*999)
          dataProjection%dataProjectionResults(dataPointGlobalNumber)%elementLocalNumber=elementLocalNumber
          dataProjection%dataProjectionResults(dataPointGlobalNumber)%elementGlobalNumber=elementGlobalNumber
          dataProjection%dataProjectionResults(dataPointGlobalNumber)%exitTag=DATA_PROJECTION_USER_SPECIFIED
        ENDIF
      ELSE
        localError="Element user number "//TRIM(NumberToVString(projectionElementUserNumber,"*",err,error))// &
          & " does not exist."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE DEFAULT
      localError="The data projection type of "//TRIM(NumberToVString(dataProjection%projectionType,"",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
      
    EXITS("DataProjection_ResultElementNumberSet")
    RETURN
999 ERRORSEXITS("DataProjection_ResultElementNumberSet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_ResultElementNumberSet
  
  !
  !================================================================================================================================
  !

  !>Sets the projection element face number for a data point identified by a given global number.
  SUBROUTINE DataProjection_ResultElementFaceNumberSet(dataProjection,dataPointGlobalNumber,projectionElementUserNumber, &
    & localFaceNormal,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection for which projection result is stored
    INTEGER(INTG), INTENT(IN) :: dataPointGlobalNumber !<The data point global number to get the projection element face number for
    INTEGER(INTG), INTENT(IN) :: projectionElementUserNumber !<The projection candidate user element number to set the result
    INTEGER(INTG), INTENT(IN) :: localFaceNormal !<The projection candidate element face normal to set the result for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: elementGlobalNumber,elementLocalNumber,localFaceNumber
    LOGICAL :: elementExists,ghostElement
    TYPE(BasisType), POINTER :: basis
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionElementsType), POINTER :: decompositionElements
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainElementsType), POINTER :: domainElements
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DataProjection_ResultElementFaceNumberSet",err,error,*999)

    CALL DataProjection_AssertIsFinished(dataProjection,err,error,*999)

    SELECT CASE(dataProjection%projectionType)
    CASE(DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE)
      CALL FlagError("Use the result element lines routine not the element faces routine for a boundary faces projection.", &
        & err,error,*999)
    CASE(DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE)
      CALL FlagError("Use the result elements routine not the element faces routine for an all elements projection.", &
        & err,error,*999)
    CASE(DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE)
      NULLIFY(decomposition)
      CALL DataProjection_DecompositionGet(dataProjection,decomposition,err,error,*999)
      NULLIFY(decompositionTopology)
      CALL Decomposition_TopologyGet(decomposition,decompositionTopology,err,error,*999)
      NULLIFY(decompositionElements)
      CALL DecompositionTopology_DecompositionElementsGet(decompositionTopology,decompositionElements,err,error,*999)      
      CALL DecompositionElements_ElementCheckExists(decompositionElements,elementUserNumber,elementExists, &
        & elementLocalNumber,ghostElement,err,error,*999)       
      IF(elementExists) THEN
        IF(.NOT.ghostElement) THEN
          CALL DecompositionElements_ElementGlobalNumberGet(decompositionElements,elementLocalNumber,elementGlobalNumber, &
            & err,error,*999)
          NULLIFY(domain)
          CALL Decomposition_DomainGet(decomposition,0,domain,err,error,*999)
          NULLIFY(domainTopology)
          CALL Domain_TopologyGet(domain,domainTopology,err,error,*999)
          NULLIFY(domainElements)
          CALL DomainTopology_ElementsGet(domainTopology,domainElements,err,error,*999)
          NULLIFY(basis)
          CALL DomainElements_BasisGet(domainElements,elementLocalNumber,basis,err,error,*999)
          CALL Basis_LocalFaceNumberGet(basis,localFaceNormal,localFaceNumber,err,error,*999)
          dataProjection%dataProjectionResults(dataPointGlobalNumber)%elementLocalNumber=elementLocalNumber
          dataProjection%dataProjectionResults(dataPointGlobalNumber)%elementGlobalNumber=elementGlobalNumber
          dataProjection%dataProjectionResults(dataPointGlobalNumber)%elementLineFaceNumber=localFaceNumber
          dataProjection%dataProjectionResults(dataPointGlobalNumber)%exitTag=DATA_PROJECTION_USER_SPECIFIED
        ENDIF
      ELSE
        localError="Element user number "//TRIM(NumberToVString(elementUserNumber,"*",err,error))// &
          & " does not exist."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE DEFAULT
      localError="The data projection type of "//TRIM(NumberToVString(dataProjection%projectionType,"",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("DataProjection_ResultElementFaceNumberSet")
    RETURN
999 ERRORSEXITS("DataProjection_ResultElementFaceNumberSet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_ResultElementFaceNumberSet

  !
  !================================================================================================================================
  !

  !>Sets the projection element line number for a data point identified by a given global number.
  SUBROUTINE DataProjection_ResultElementLineNumberSet(dataProjection,dataPointGlobalNumber,projectionElementUserNumber, &
    & localLineNormals,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection for which projection result is stored
    INTEGER(INTG), INTENT(IN) :: dataPointGlobalNumber !<The data projection global number to set the element line number for
    INTEGER(INTG), INTENT(IN) :: projectedElementUserNumber !<The projection candidate user element number to set the result
    INTEGER(INTG), INTENT(IN) :: localLineNormals(:) !<localLineNormals(normalDirectionIdx). The projection candidate element line normals to set the result for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: elementGlobalNumber,elementLocalNumber,localLineNumber
    LOGICAL :: elementExists,ghostElement
    TYPE(BasisType), POINTER :: basis
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionElementsType), POINTER :: decompositionElements
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainElementsType), POINTER :: domainElements
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(VARYING_STRING) :: localError
     
    ENTERS("DataProjection_ResultElementLineNumberSet",err,error,*999)

    CALL DataProjection_AssertIsFinished(dataProjection,err,error,*999)
    
#ifdef WITH_PRECHECKS    
    IF(SIZE(localLineNormals,1)/=2) THEN
      localError="The size of the line normals array of "// &
        & TRIM(NumberToVString(SIZE(localLineNormals,1),"*",err,error))//" is invalid. The size should be 2."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    SELECT CASE(dataProjection%projectionType)
    CASE(DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE)
      CALL FlagError("Use the result element faces routine not the element lines routine for a boundary lines projection.", &
        & err,error,*999)
    CASE(DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE)
      CALL FlagError("Use the result elements routine not the element lines routine for an all elements projection.", &
        & err,error,*999)
    CASE(DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE)
      NULLIFY(decomposition)
      CALL DataProjection_DecompositionGet(dataProjection,decomposition,err,error,*999)
      NULLIFY(decompositionTopology)
      CALL Decomposition_TopologyGet(decomposition,decompositionTopology,err,error,*999)
      NULLIFY(decompositionElements)
      CALL DecompositionTopology_DecompositionElementsGet(decompositionTopology,decompositionElements,err,error,*999)
      CALL DecompositionElements_ElementCheckExists(decompositionElements,projectionElementUserNumber,elementExists, &
        & elementLocalNumber,ghostElement,err,error,*999)       
      IF(elementExists) THEN
        IF(.NOT.ghostElement) THEN
          CALL DecompositionElements_ElementGlobalNumberGet(decompositionElements,elementLocalNumber,elementGlobalNumber, &
            & err,error,*999)
          NULLIFY(domain)
          CALL Decomposition_DomainGet(decomposition,0,domain,err,error,*999)
          NULLIFY(domainTopology)
          CALL Domain_TopologyGet(domain,domainTopology,err,error,*999)
          NULLIFY(domainElements)
          CALL DomainTopology_ElementsGet(domainTopology,domainElements,err,error,*999)
          NULLIFY(basis)
          CALL DomainElements_BasisGet(domainElements,elementLocalNumber,basis,err,error,*999)
          CALL Basis_LocalLineNumberGet(basis,localLineNormals,localLineNumber,err,error,*999)
          dataProjection%dataProjectionResults(dataPointGlobalNumber)%elementLocalNumber=elementLocalNumber
          dataProjection%dataProjectionResults(dataPointGlobalNumber)%elementGlobalNumber=elementGlobalNumber
          dataProjection%dataProjectionResults(dataPointGlobalNumber)%elementLineFaceNumber=localLineNumber
          dataProjection%dataProjectionResults(dataPointGlobalNumber)%exitTag=DATA_PROJECTION_USER_SPECIFIED
        ENDIF
      ELSE
        localError="Element user number "//TRIM(NumberToVString(elementUserNumber,"*",err,error))// &
          & " does not exist."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE DEFAULT
      localError="The data projection type of "//TRIM(NumberToVString(dataProjection%projectionType,"",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("DataProjection_ResultElementLineNumberSet")
    RETURN
999 ERRORSEXITS("DataProjection_ResultElementLineNumberSet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_ResultElementLineNumberSet

  !
  !================================================================================================================================
  !

  !>Sets the projection xi for a data point identified by a given global number.
  SUBROUTINE DataProjection_ResultProjectionXiSet(dataProjection,dataPointGlobalNumber,projectionXi,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection for which projection result is stored
    INTEGER(INTG), INTENT(IN) :: dataPointGlobalNumber !<The data point global number to set the projection xi for
    REAL(DP), INTENT(IN) :: projectionXi(:) !<projectionXi(xiIdx). The projection xi of the specified global data point to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: elementLocalNumber, localLineFaceNumber
    TYPE(BasisType), POINTER :: basis
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainElementsType), POINTER :: domainElements
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DataProjection_ResultProjectionXiSet",err,error,*999)

    CALL DataProjection_AssertNotFinished(dataProjection,err,error,*999)
#ifdef WITH_PRECHECKS
    IF(.NOT.ALLOCATED(dataProjection%dataProjectionResults)) &
      & CALL FlagError("Data projection data projection results is not allocated.",err,error,*999)
    IF(SIZE(projectionXi,1)/=dataProjection%numberOfXi) THEN
      localError="The specified projection xi has size of "//TRIM(NumberToVString(SIZE(projectionXi,1),"*",err,error))// &
        & "but it needs to have size of "//TRIM(NumberToVString(dataProjection%numberOfXi,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(ANY(projectionXi<0.0_DP).OR.ANY(projectionXi>1.0_DP)) &
      & CALL FlagError("The specified xi location is invalid. The xi location must be >= 0.0 and <= 1.0",err,error,*999)
#endif    

    dataProjection%dataProjectionResults(dataPointGlobalNumber)%xi(1:dataProjection%numberOfXi)= &
      & projectionXi(1:dataProjection%numberOfXi)
    
    !Compute full elemental xi
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

    IF(dataProjection%numberOfXi==dataProjection%numberOfElementXi) THEN
      dataProjection%dataProjectionResults(dataPointGlobalNumber)%elementXi= &
        & dataProjection%dataProjectionResults(dataPointGlobalNumber)%xi
    ELSE
      elementLocalNumber=dataProjection%dataProjectionResults(dataPointGlobalNumber)%elementLocalNumber
      localLineFaceNumber=dataProjection%dataProjectionResults(dataPointGlobalNumber)%elementLineFaceNumber
      NULLIFY(basis)
      CALL DomainElements_ElementBasisGet(domainElements,elementLocalNumber,basis,err,error,*999)
      CALL Basis_BoundaryXiToXi(basis,localLineFaceNumber,dataProjection% &
        & dataProjectionResults(dataPointGlobalNumber)%xi(1:dataProjection%numberOfXi),dataProjection% &
        & dataProjectionResults(dataPointGlobalNumber)%elementXi,err,error,*999)
    ENDIF

    EXITS("DataProjection_ResultProjectionXiSet")
    RETURN
999 ERRORSEXITS("DataProjection_ResultProjectionXiSet",err,error)
    RETURN 1

  END SUBROUTINE DataProjection_ResultProjectionXiSet

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

