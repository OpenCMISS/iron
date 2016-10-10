!> \file
!> \author Tim Wu
!> \brief This module handles all data point routines.
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

!> This module handles all data point routines.

MODULE DataPointRoutines

  USE BASE_ROUTINES
  USE COMP_ENVIRONMENT
  USE COORDINATE_ROUTINES
  USE DataProjectionRoutines
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE Kinds
  USE Strings
  USE Trees
  USE Types

#include "macros.h"

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  !>Starts the process of creating data points for an interface or region
  INTERFACE DataPoints_CreateStart
    MODULE PROCEDURE DataPoints_CreateStartRegion
    MODULE PROCEDURE DataPoints_CreateStartInterface
  END INTERFACE DataPoints_CreateStart

  !>Initialises data points for an interface or region
  INTERFACE DataPoints_Initialise
    MODULE PROCEDURE DataPoints_InitialiseRegion
    MODULE PROCEDURE DataPoints_InitialiseInterface
  END INTERFACE DataPoints_Initialise

  !>Gets the label for a data point identified by a given global number.
  INTERFACE DataPoints_LabelGet
    MODULE PROCEDURE DataPoints_LabelGetC
    MODULE PROCEDURE DataPoints_LabelGetVS
  END INTERFACE DataPoints_LabelGet

  !>Changes/sets the label for a data point identified by a given global number.
  INTERFACE DataPoints_LabelSet
    MODULE PROCEDURE DataPoints_LabelSetC
    MODULE PROCEDURE DataPoints_LabelSetVS
  END INTERFACE DataPoints_LabelSet

  PUBLIC DataPoint_CheckExists

  PUBLIC DataPoints_CreateFinish,DataPoints_CreateStart,DataPoints_Destroy
  
  PUBLIC DataPoints_DataProjectionGet,DataPoints_DataProjectionGlobalNumberGet

  PUBLIC DataPoints_GlobalNumberGet

  PUBLIC DataPoints_LabelGet,DataPoints_LabelSet
  
  PUBLIC DataPoints_NumberOfDataPointsGet
  
  PUBLIC DataPoints_PositionGet,DataPoints_PositionSet

  PUBLIC DataPoints_UserNumberGet,DataPoints_UserNumberSet
  
  PUBLIC DataPoints_WeightsGet,DataPoints_WeightsSet

CONTAINS

  !
  !================================================================================================================================
  !

  !>Checks that a user data point number is defined on the specified region.
  SUBROUTINE DataPoint_CheckExists(dataPoints,userNumber,dataPointExists,globalNumber,err,error,*)

    !Argument variables
    TYPE(DataPointsType), POINTER :: dataPoints !<A pointer to the data points to check
    INTEGER(INTG) :: userNumber !<The user data point number to check if it exists
    LOGICAL, INTENT(OUT) :: dataPointExists !<On exit, is .TRUE. if the data point user number exists in the region, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: globalNumber !<On exit, if the data point exists the global number corresponding to the user data point number. If the data point does not exist then global number will be 0.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(TREE_NODE_TYPE), POINTER :: treeNode
   
    ENTERS("DataPoint_CheckExists",err,error,*999)

    dataPointExists=.FALSE.
    globalNumber=0
    IF(ASSOCIATED(dataPoints)) THEN
      NULLIFY(treeNode)
      CALL Tree_Search(dataPoints%dataPointsTree,userNumber,treeNode,err,error,*999)
      IF(ASSOCIATED(treeNode)) THEN
        CALL Tree_NodeValueGet(dataPoints%dataPointsTree,treeNode,globalNumber,err,error,*999)
        dataPointExists=.TRUE.
      ENDIF
    ELSE
      CALL FlagError("Data points is not associated.",err,error,*999)
    ENDIF

    EXITS("DataPoint_CheckExists")
    RETURN
999 ERRORSEXITS("DataPoint_CheckExists",err,error)
    RETURN 1
   
  END SUBROUTINE DataPoint_CheckExists

  !
  !================================================================================================================================
  !

  !>Finalises a data point and deallocates all memory
  SUBROUTINE DataPoint_Finalise(dataPoint,err,error,*)
    
    !Argument variables
    TYPE(DataPointType),INTENT(OUT) :: dataPoint !<The data point to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DataPoint_Finalise",err,error,*999)

    dataPoint%globalNumber=0
    dataPoint%userNumber=0
    dataPoint%label=""
    IF(ALLOCATED(dataPoint%position)) DEALLOCATE(dataPoint%position)
    IF(ALLOCATED(dataPoint%weights)) DEALLOCATE(dataPoint%weights)
    
    EXITS("DataPoint_Finalise")
    RETURN
999 ERRORSEXITS("DataPoint_Finalise",err,error)
    RETURN 1  
 
  END SUBROUTINE DataPoint_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises a data point
  SUBROUTINE DataPoint_Initialise(dataPoint,err,error,*)
    
    !Argument variables
    TYPE(DataPointType),INTENT(OUT) :: dataPoint !<The data point to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DataPoint_Initialise",err,error,*999)

    dataPoint%globalNumber=0
    dataPoint%userNumber=0
    dataPoint%label=""
    
    EXITS("DataPoint_Initialise")
    RETURN
999 ERRORSEXITS("DataPoint_Initialise",err,error)
    RETURN 1  
 
  END SUBROUTINE DataPoint_Initialise

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating data points in the region.
  SUBROUTINE DataPoints_CreateFinish(dataPoints,err,error,*)

    !Argument variables
    TYPE(DataPointsType), POINTER :: dataPoints !<A pointer to the data points to be finished
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointIdx
    
    ENTERS("DataPoints_CreateFinish",err,error,*999)

    IF(ASSOCIATED(dataPoints)) THEN
      IF(dataPoints%dataPointsFinished) THEN
        CALL FlagError("Data points have already been finished.",err,error,*999)
      ELSE
        dataPoints%dataPointsFinished=.TRUE.
      ENDIF
    ELSE
      CALL FlagError("Data points is not associated.",err,error,*999)
    ENDIF
    
    IF(diagnostics1) THEN 
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,  "Data points:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,  "Number of dimensions  = ",dataPoints%numberOfDimensions,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,  "Number of data points = ",dataPoints%numberOfDataPoints,err,error,*999)
      DO dataPointIdx=1,dataPoints%numberOfDataPoints
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Data Point : ",dataPointIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Global number = ",dataPoints%dataPoints(dataPointIdx)% &
          & globalNumber,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    User number   = ",dataPoints%dataPoints(dataPointIdx)% &
          & userNumber,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Label         = ",dataPoints%dataPoints(dataPointIdx)%label, &
          & err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,dataPoints%numberOfDimensions,3,3,dataPoints%dataPoints(dataPointIdx)% &
          & position,'("    Position      =",3(X,E13.6))','(19X,3(X,E13.6))',err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,dataPoints%numberOfDimensions,3,3,dataPoints%dataPoints(dataPointIdx)% &
          & weights,'("    Weights       =",3(X,E13.6))','(19X,3(X,E13.6))',err,error,*999)
      ENDDO !dataPointIdx
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Data points User->Global number tree:",err,error,*999)
      CALL Tree_Output(DIAGNOSTIC_OUTPUT_TYPE,dataPoints%dataPointsTree,err,error,*999)
    ENDIF

    EXITS("DataPoints_CreateFinish")
    RETURN
999 ERRORSEXITS("DataPoints_CreateFinish",err,error)
    RETURN 1
   
  END SUBROUTINE DataPoints_CreateFinish

  !
  !================================================================================================================================
  !

  !>Starts the process of creating generic data points
  SUBROUTINE DataPoints_CreateStartGeneric(numberOfDataPoints,numberOfDimensions,dataPoints,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: numberOfDataPoints !<The number of data points to create
    INTEGER(INTG), INTENT(IN) :: numberOfDimensions !<The number of dimensions for data points values
    TYPE(DataPointsType), POINTER :: dataPoints !<On exit, a pointer to the created data points
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: insertStatus,dataPointIdx,coord_idx
    TYPE(VARYING_STRING) :: localError

    ENTERS("DataPoints_CreateStartGeneric",err,error,*999)

    IF(ASSOCIATED(dataPoints)) THEN
      IF(numberOfDataPoints>0) THEN
        IF(numberOfDimensions>=1.AND.numberOfDimensions<=3) THEN          
          ALLOCATE(dataPoints%dataPoints(numberOfDataPoints),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate data points data points.",err,error,*999)
          dataPoints%numberOfDataPoints=numberOfDataPoints
          dataPoints%numberOfDimensions=numberOfDimensions
          dataPoints%numberOfDataProjections=0
          IF(ALLOCATED(dataPoints%dataProjections)) DEALLOCATE(dataPoints%dataProjections)
          CALL Tree_CreateStart(dataPoints%dataProjectionsTree,err,error,*999)
          CALL Tree_InsertTypeSet(dataPoints%dataProjectionsTree,TREE_NO_DUPLICATES_ALLOWED,err,error,*999)
          CALL Tree_CreateFinish(dataPoints%dataProjectionsTree,err,error,*999)
          CALL Tree_CreateStart(dataPoints%dataPointsTree,err,error,*999)
          CALL Tree_InsertTypeSet(dataPoints%dataPointsTree,TREE_NO_DUPLICATES_ALLOWED,err,error,*999)
          CALL Tree_CreateFinish(dataPoints%dataPointsTree,err,error,*999)
          !Set default data point numbers
          DO dataPointIdx=1,dataPoints%numberOfDataPoints
            CALL DataPoint_Initialise(dataPoints%dataPoints(dataPointIdx),err,error,*999)
            !Default the user number to the global number
            dataPoints%dataPoints(dataPointIdx)%globalNumber=dataPointIdx
            dataPoints%dataPoints(dataPointIdx)%userNumber=dataPointIdx
            ALLOCATE(dataPoints%dataPoints(dataPointIdx)%position(numberOfDimensions),STAT=err)
            IF(err/=0) THEN
              localError="Could not allocate data points data position for data point number "// &
                & TRIM(NumberToVString(dataPointIdx,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
            ALLOCATE(dataPoints%dataPoints(dataPointIdx)%weights(numberOfDimensions),STAT=err)
            IF(err/=0) THEN
              localError="Could not allocate data points data weights for data point number "// &
                & TRIM(NumberToVString(dataPointIdx,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
            !Initialise data points position to 0.0 and weights to 1.0
            dataPoints%dataPoints(dataPointIdx)%position=0.0_DP
            dataPoints%dataPoints(dataPointIdx)%weights=1.0_DP
            CALL Tree_ItemInsert(dataPoints%dataPointsTree,dataPointIdx,dataPointIdx,insertStatus,err,error,*999)
          ENDDO !dataPointIdx
        ELSE
          localError="The specified number of dimensions of "//TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
            & " is invalid. The number of dimensions must be >= 1 and <= 3."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        localError="The specified number of data points of "//TRIM(NumberToVString(numberOfDataPoints,"*",err,error))// &
          & " is invalid. The number of data points must be > 0."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Data points is not associated.",err,error,*999)
    ENDIF

    EXITS("DataPoints_CreateStartGeneric")
    RETURN  
999 ERRORSEXITS("DataPoints_CreateStartGeneric",err,error)
    RETURN 1
   
  END SUBROUTINE DataPoints_CreateStartGeneric

  !
  !================================================================================================================================
  !

  !>Starts the process of creating data points in an interface.
  SUBROUTINE DataPoints_CreateStartInterface(interface,numberOfDataPoints,dataPoints,err,error,*)

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: interface !<A pointer to the interface in which to create the data points
    INTEGER(INTG), INTENT(IN) :: numberOfDataPoints !<The number of data points to create
    TYPE(DataPointsType), POINTER :: dataPoints !<On exit, a pointer to the created data points. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("DataPoints_CreateStartInterface",err,error,*998)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(ASSOCIATED(INTERFACE%COORDINATE_SYSTEM)) THEN
        IF(ASSOCIATED(dataPoints)) THEN
          CALL FlagError("Data points is already associated.",err,error,*999)
        ELSE
          IF(ASSOCIATED(INTERFACE%dataPoints)) THEN
            CALL FlagError("Interface already has data points associated.",err,error,*998)
          ELSE
            !Initialise the data points for the interface
            CALL DataPoints_Initialise(INTERFACE,err,error,*999)
            !Create the data points 
            CALL DataPoints_CreateStartGeneric(numberOfDataPoints,INTERFACE%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS, &
              & INTERFACE%dataPoints,err,error,*999)
            !Return the pointer        
            dataPoints=>interface%dataPoints
          ENDIF
        ENDIF
      ELSE
        CALL FlagError("Interface coordinate system is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Interface is not associated.",err,error,*998)
    ENDIF
    
    EXITS("DataPoints_CreateStartInterface")
    RETURN
999 CALL DataPoints_Finalise(interface%dataPoints,dummyErr,dummyError,*998)    
998 ERRORSEXITS("DataPoints_CreateStartInterface",err,error)
    RETURN 1
   
  END SUBROUTINE DataPoints_CreateStartInterface

  !
  !================================================================================================================================
  !

  !>Starts the process of creating data points in an region.
  SUBROUTINE DataPoints_CreateStartRegion(region,numberOfDataPoints,dataPoints,err,error,*)

    !Argument variables
    TYPE(region_TYPE), POINTER :: region !<A pointer to the region in which to create the data points
    INTEGER(INTG), INTENT(IN) :: numberOfDataPoints !<The number of data points to create
    TYPE(DataPointsType), POINTER :: dataPoints !<On exit, a pointer to the created data points. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("DataPoints_CreateStartRegion",err,error,*998)

    IF(ASSOCIATED(region)) THEN
      IF(ASSOCIATED(region%COORDINATE_SYSTEM)) THEN
        IF(ASSOCIATED(dataPoints)) THEN
          CALL FlagError("Data points is already associated.",err,error,*999)
        ELSE
          IF(ASSOCIATED(region%dataPoints)) THEN
            CALL FlagError("Region already has data points associated.",err,error,*998)
          ELSE
            !Initialise the data points for the region
            CALL DataPoints_Initialise(region,err,error,*999)
            !Create the data points 
            CALL DataPoints_CreateStartGeneric(numberOfDataPoints,region%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS, &
              & region%dataPoints,err,error,*999)
            !Return the pointer        
            dataPoints=>region%dataPoints
          ENDIF
        ENDIF
      ELSE
        CALL FlagError("Region coordinate system is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Region is not associated.",err,error,*998)
    ENDIF
    
    EXITS("DataPoints_CreateStartRegion")
    RETURN
999 CALL DataPoints_Finalise(region%dataPoints,dummyErr,dummyError,*998)    
998 ERRORSEXITS("DataPoints_CreateStartRegion",err,error)
    RETURN 1
   
  END SUBROUTINE DataPoints_CreateStartRegion
     
  !
  !================================================================================================================================
  !

  !>Destroys data points. \see OPENCMISS::Iron::cmfe_DataPointsDestroy
  SUBROUTINE DataPoints_Destroy(dataPoints,err,error,*)

    !Argument variables
    TYPE(DataPointsType), POINTER :: dataPoints !<A pointer to the data points to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("DataPoints_Destroy",err,error,*999)

    IF(ASSOCIATED(dataPoints)) THEN
      IF(ASSOCIATED(dataPoints%region)) THEN
        NULLIFY(dataPoints%region%dataPoints)
      ELSE
        IF(ASSOCIATED(dataPoints%INTERFACE)) THEN
          NULLIFY(dataPoints%interface%dataPoints)
        ELSE
          CALL FLAG_ERROR("Data_points region or interface is not associated.",err,error,*999)
        ENDIF
      ENDIF
      CALL DataPoints_Finalise(dataPoints,err,error,*999)
    ENDIF
   
    EXITS("DataPoints_Destroy")
    RETURN
999 ERRORSEXITS("DataPoints_Destroy",err,error)
    RETURN 1
   
  END SUBROUTINE DataPoints_Destroy

  !
  !===============================================================================================================================
  !

  !>Finalises the data points and deallocates any memory. 
  SUBROUTINE DataPoints_Finalise(dataPoints,err,error,*)

    !Argument variables
    TYPE(DataPointsType), POINTER :: dataPoints !<A pointer to the data points to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointIdx,dataProjectionIdx

    ENTERS("DataPoints_Finalise",err,error,*999)

    IF(ASSOCIATED(dataPoints)) THEN
      IF(ALLOCATED(dataPoints%dataPoints)) THEN
        DO dataPointIdx=1,SIZE(dataPoints%dataPoints,1)
          CALL DataPoint_Finalise(dataPoints%dataPoints(dataPointIdx),err,error,*999)
        ENDDO !dataPointIdx
        DEALLOCATE(dataPoints%dataPoints)
      ENDIF
      IF(ASSOCIATED(dataPoints%dataPointsTree)) CALL Tree_Destroy(dataPoints%dataPointsTree,err,error,*999)
      dataPoints%numberOfDataProjections=0
      IF(ALLOCATED(dataPoints%dataProjections)) THEN
        DO dataProjectionIdx=1,SIZE(dataPoints%dataProjections,1)
          CALL DataProjection_Destroy(dataPoints%dataProjections(dataProjectionIdx)%PTR,err,error,*999)
        ENDDO !dataProjectionIdx
        DEALLOCATE(dataPoints%dataProjections)
      ENDIF
      IF(ASSOCIATED(dataPoints%dataProjectionsTree)) CALL Tree_Destroy(dataPoints%dataProjectionsTree,err,error,*999)
      DEALLOCATE(dataPoints)
    ENDIF
    
    EXITS("DataPoints_Finalise")
    RETURN
999 ERRORSEXITS("DataPoints_Finalise",err,error)
    RETURN 1

  END SUBROUTINE DataPoints_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises the data points.
  SUBROUTINE DataPoints_InitialiseGeneric(dataPoints,err,error,*)

    !Argument variables
    TYPE(DataPointsType), POINTER :: dataPoints !<A pointer to the data points to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DataPoints_InitialiseGeneric",err,error,*999)

    IF(ASSOCIATED(dataPoints)) THEN
      NULLIFY(dataPoints%region)
      NULLIFY(dataPoints%interface)
      dataPoints%dataPointsFinished=.FALSE.
      dataPoints%numberOfDimensions=0
      dataPoints%numberOfDataPoints=0
      NULLIFY(dataPoints%dataPointsTree)
      dataPoints%numberOfDataProjections=0
      NULLIFY(dataPoints%dataProjectionsTree)
    ELSE
      CALL FlagError("Data points is not associated.",err,error,*999)
    ENDIF
    
    EXITS("DataPoints_InitialiseGeneric")
    RETURN
999 ERRORSEXITS("DataPoints_InitialiseGeneric",err,error)
    RETURN 1
  END SUBROUTINE DataPoints_InitialiseGeneric

  !
  !================================================================================================================================
  !

  !>Initialises the data points in a given interface.
  SUBROUTINE DataPoints_InitialiseInterface(INTERFACE,err,error,*)

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: interface !<A pointer to the interface to initialise the data points for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DataPoints_InitialiseInterface",err,error,*999)

    IF(ASSOCIATED(interface)) THEN
      IF(ASSOCIATED(interface%dataPoints)) THEN
        CALL FlagError("Interface already has associated data points.",err,error,*999)
      ELSE
        ALLOCATE(INTERFACE%dataPoints,STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate interface data points.",err,error,*999)
        CALL DataPoints_InitialiseGeneric(interface%dataPoints,err,error,*999)
        interface%dataPoints%interface=>interface
      ENDIF
    ELSE
      CALL FlagError("Interface is not associated.",err,error,*999)
    ENDIF
    
    EXITS("DataPoints_InitialiseInterface")
    RETURN
999 ERRORSEXITS("DataPoints_InitialiseInterface",err,error)
    RETURN 1
    
  END SUBROUTINE DataPoints_InitialiseInterface

  !
  !================================================================================================================================
  !

  !>Initialises the data points in a given region.
  SUBROUTINE DataPoints_InitialiseRegion(region,err,error,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: region !<A pointer to the region to initialise the data points for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DataPoints_InitialiseRegion",err,error,*999)

    IF(ASSOCIATED(region)) THEN
      IF(ASSOCIATED(region%dataPoints)) THEN
        CALL FlagError("Region has associated data points.",err,error,*999)
      ELSE
        ALLOCATE(region%dataPoints,STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate region data points.",err,error,*999)
        CALL DataPoints_InitialiseGeneric(region%dataPoints,err,error,*999)
        region%dataPoints%region=>region
      ENDIF
    ELSE
      CALL FlagError("Region is not associated.",err,error,*999)
    ENDIF
    
    EXITS("DataPoints_InitialiseRegion")
    RETURN
999 ERRORSEXITS("DataPoints_InitialiseRegion",err,error)
    RETURN 1
    
  END SUBROUTINE DataPoints_InitialiseRegion

  !
  !================================================================================================================================
  !  

  !>Gets the global number for a data point identified by a given user number. 
  SUBROUTINE DataPoints_GlobalNumberGet(dataPoints,userNumber,globalNumber,err,error,*)

    !Argument variables
    TYPE(DataPointsType), POINTER :: dataPoints !<A pointer to the data points to get the number for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number to get the global number for
    INTEGER(INTG), INTENT(OUT) :: globalNumber !<On exit, the global number of the specified user data point
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    LOGICAL :: dataPointExists
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DataPoints_GlobalNumberGet",err,error,*999)

    IF(ASSOCIATED(dataPoints)) THEN
      CALL DataPoint_CheckExists(dataPoints,userNumber,dataPointExists,globalNumber,err,error,*999)
      IF(.NOT.dataPointExists) THEN
        localError="The data point with user number "//TRIM(NumberToVString(userNumber,"*",err,error))//" does not exist."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Data points is not associated.",err,error,*999)
    ENDIF
    
    EXITS("DataPoints_GlobalNumberGet")
    RETURN
999 ERRORSEXITS("DataPoints_GlobalNumberGet",err,error)    
    RETURN 1
   
  END SUBROUTINE DataPoints_GlobalNumberGet
       
  !
  !================================================================================================================================
  !

  !>Gets the character label for a data point identified by a given user number.
  SUBROUTINE DataPoints_LabelGetC(dataPoints,userNumber,label,err,error,*)

    !Argument variables
    TYPE(DataPointsType), POINTER :: dataPoints !<A pointer to the data points to get the label for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number to get the label for
    CHARACTER(LEN=*), INTENT(OUT) :: label !<On exit, the label of the specified user data point
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER :: cLength,globalNumber,vsLength
    
    ENTERS("DataPoints_LabelGetC",err,error,*999)

    IF(ASSOCIATED(dataPoints)) THEN
      IF(dataPoints%dataPointsFinished) THEN
        CALL DataPoints_GlobalNumberGet(dataPoints,userNumber,globalNumber,err,error,*999)
        cLength=LEN(label)
        vsLength=LEN_TRIM(dataPoints%dataPoints(globalNumber)%label)
        IF(cLength>vsLength) THEN
          label=CHAR(LEN_TRIM(dataPoints%dataPoints(globalNumber)%label))
        ELSE
          label=CHAR(dataPoints%dataPoints(globalNumber)%label,cLength)
        ENDIF
      ELSE
        CALL FlagError("Data points have not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Data points is not associated.",err,error,*999)
    ENDIF
    
    EXITS("DataPoints_LabelGetC")
    RETURN
999 ERRORSEXITS("DataPoints_LabelGetC",err,error)    
    RETURN 1
   
  END SUBROUTINE DataPoints_LabelGetC
        
  !
  !================================================================================================================================
  !

  !>Gets the varying string label for a data point identified by a given user number. 
  SUBROUTINE DataPoints_LabelGetVS(dataPoints,userNumber,label,err,error,*)

    !Argument variables
    TYPE(DataPointsType), POINTER :: dataPoints !<A pointer to the data points to get the label for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number to get the label for
    TYPE(VARYING_STRING), INTENT(OUT) :: label !<On exit, the label of the specified global data point
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: globalNumber
    
    ENTERS("DataPoints_LabelGetVS",err,error,*999)

    IF(ASSOCIATED(dataPoints)) THEN
      IF(dataPoints%dataPointsFinished) THEN
        CALL DataPoints_GlobalNumberGet(dataPoints,userNumber,globalNumber,err,error,*999)
        label=dataPoints%dataPoints(globalNumber)%label
      ELSE
        CALL FlagError("Data points have not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Data points is not associated.",err,error,*999)
    ENDIF
    
    EXITS("DataPoints_LabelGetVS")
    RETURN
999 ERRORSEXITS("DataPoints_LabelGetVS",err,error)    
    RETURN 1
   
  END SUBROUTINE DataPoints_LabelGetVS

  !
  !================================================================================================================================
  !

  !>Changes/sets the character label for a data point identified by a given user number.
  SUBROUTINE DataPoints_LabelSetC(dataPoints,userNumber,label,err,error,*)

    !Argument variables
    TYPE(DataPointsType), POINTER :: dataPoints !<A pointer to the data points to set the label for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number to set the label for
    CHARACTER(LEN=*), INTENT(IN) :: label !<The label to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: globalNumber
    
    ENTERS("DataPoints_LabelSetC",err,error,*999)

    IF(ASSOCIATED(dataPoints)) THEN
      IF(dataPoints%dataPointsFinished) THEN
        CALL FlagError("Data points have been finished.",err,error,*999)
      ELSE
        CALL DataPoints_GlobalNumberGet(dataPoints,userNumber,globalNumber,err,error,*999)
        dataPoints%dataPoints(globalNumber)%label=label
      ENDIF
    ELSE
      CALL FlagError("Data points is not associated.",err,error,*999)
    ENDIF
    
    EXITS("DataPoints_LabelSetC")
    RETURN
999 ERRORSEXITS("DataPoints_LabelSetC",err,error)    
    RETURN 1
   
  END SUBROUTINE DataPoints_LabelSetC    
  
  !
  !================================================================================================================================
  !


  !>Changes/sets the varying string label for a data point identified by a given user number.
  SUBROUTINE DataPoints_LabelSetVS(dataPoints,userNumber,label,err,error,*)

    !Argument variables
    TYPE(DataPointsType), POINTER :: dataPoints !<A pointer to the data points to set the label for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number to set the label for
    TYPE(VARYING_STRING), INTENT(IN) :: label !<The label to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: globalNumber
    
    ENTERS("DataPoints_LabelSetVS",err,error,*999)

    IF(ASSOCIATED(dataPoints)) THEN
      IF(dataPoints%dataPointsFinished) THEN
        CALL FlagError("Data points have been finished.",err,error,*999)
      ELSE
        CALL DataPoints_GlobalNumberGet(dataPoints,userNumber,globalNumber,err,error,*999)
        dataPoints%dataPoints(globalNumber)%label=label
      ENDIF
    ELSE
      CALL FlagError("Data points is not associated.",err,error,*999)
    ENDIF
    
    EXITS("DataPoints_LabelSetVS")
    RETURN
999 ERRORSEXITS("DataPoints_LabelSetVS",err,error)    
    RETURN 1
   
  END SUBROUTINE DataPoints_LabelSetVS
        
  !
  !================================================================================================================================
  !
  
  !>Gets the position for a data point identified by a given user number.
  SUBROUTINE DataPoints_PositionGet(dataPoints,userNumber,position,err,error,*)

    !Argument variables
    TYPE(DataPointsType), POINTER :: dataPoints !<A pointer to the data points to set the number for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number to get the values for
    REAL(DP), INTENT(OUT) :: position(:) !<position(coordinateIdx). On exit, the position of the specified data point
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: globalNumber
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DataPoints_PositionGet",err,error,*999)

    IF(ASSOCIATED(dataPoints)) THEN
      IF(dataPoints%dataPointsFinished) THEN
        IF(SIZE(position,1)>=dataPoints%numberOfDimensions) THEN
          CALL DataPoints_GlobalNumberGet(dataPoints,userNumber,globalNumber,err,error,*999)
          position(1:dataPoints%numberOfDimensions)=dataPoints%dataPoints(globalNumber)%position(1:dataPoints%numberOfDimensions)
        ELSE
          localError="The size of the specified position array of "//TRIM(NumberToVString(SIZE(position,1),"*",err,error))// &
            & " is too small. The array size needs to be >= "// &
            & TRIM(NumberToVString(dataPoints%numberOfDimensions,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Data points have not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Data points is not associated.",err,error,*999)
    ENDIF
    
    EXITS("DataPoints_PositionGet")
    RETURN
999 ERRORSEXITS("DataPoints_PositionGet",err,error)    
    RETURN 1

  END SUBROUTINE DataPoints_PositionGet

  !
  !================================================================================================================================
  !

  !>Changes/sets the position for a data point identified by a given user number.
  SUBROUTINE DataPoints_PositionSet(dataPoints,userNumber,position,err,error,*)

    !Argument variables
    TYPE(DataPointsType), POINTER :: dataPoints !<A pointer to the data points to set the number for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number to set the values for
    REAL(DP), INTENT(IN) :: position(:) !<position(coordinateIdx). The data point position to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: globalNumber
    TYPE(VARYING_STRING) :: localError
   
    ENTERS("DataPoints_PositionSet",err,error,*999)

    IF(ASSOCIATED(dataPoints)) THEN   
      IF(dataPoints%dataPointsFinished) THEN
        CALL FlagError("Data points have been finished.",err,error,*999)
      ELSE
        IF(SIZE(position,1)>=dataPoints%numberOfDimensions) THEN
          CALL DataPoints_GlobalNumberGet(dataPoints,userNumber,globalNumber,err,error,*999)
          dataPoints%dataPoints(globalNumber)%position(1:dataPoints%numberOfDimensions)=position(1:dataPoints%numberOfDimensions)
        ELSE
          localError="The size of the specified position array of "//TRIM(NumberToVString(SIZE(position,1),"*",err,error))// &
            & " is too small. The array size needs to be >= "// &
            & TRIM(NumberToVString(dataPoints%numberOfDimensions,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Data points is not associated.",err,error,*999)
    ENDIF
    
    EXITS("DataPoints_PositionSet")
    RETURN
999 ERRORSEXITS("DataPoints_PositionSet",err,error)    
    RETURN 1

  END SUBROUTINE DataPoints_PositionSet

  !
  !================================================================================================================================
  !

  !>Returns the number of data points. \see OPENCMISS::Iron::cmfe_DataPointsNumberOfDataPointsGet
  SUBROUTINE DataPoints_NumberOfDataPointsGet(dataPoints,numberOfDataPoints,err,error,*)

    !Argument variables
    TYPE(DataPointsType), POINTER :: dataPoints !<A pointer to the data points to get the number of data points for
    INTEGER(INTG), INTENT(OUT) :: numberOfDataPoints !<On return, the number of data points
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("DataPoints_NumberOfDataPointsGet",err,error,*999)

    IF(ASSOCIATED(dataPoints)) THEN
      IF(dataPoints%dataPointsFinished) THEN
        numberOfDataPoints=dataPoints%numberOfDataPoints
      ELSE
        CALL FlagError("Data points have not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Data points is not associated.",err,error,*999)
    ENDIF
    
    EXITS("DataPoints_NumberOfDataPointsGet")
    RETURN
999 ERRORSEXITS("DataPoints_NumberOfDataPointsGet",err,error)    
    RETURN 1
   
  END SUBROUTINE DataPoints_NumberOfDataPointsGet

  !
  !================================================================================================================================
  !  

  !>Gets the user number for a data point identified by a given global number. 
  SUBROUTINE DataPoints_UserNumberGet(dataPoints,globalNumber,userNumber,err,error,*)

    !Argument variables
    TYPE(DataPointsType), POINTER :: dataPoints !<A pointer to the data points to get the number for
    INTEGER(INTG), INTENT(IN) :: globalNumber !<The global number to get the user number for
    INTEGER(INTG), INTENT(OUT) :: userNumber !<On exit, the user number of the specified global data point
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DataPoints_UserNumberGet",err,error,*999)

    IF(ASSOCIATED(dataPoints)) THEN
      IF(dataPoints%dataPointsFinished) THEN
        IF(globalNumber>=1.AND.globalNumber<=dataPoints%numberOfDataPoints) THEN
          userNumber=dataPoints%dataPoints(globalNumber)%userNumber
        ELSE
          localError="The specified global data point number of "//TRIM(NumberToVString(globalNumber,"*",err,error))// &
            & " is invalid. The global data point number should be between 1 and "// &
            & TRIM(NumberToVString(dataPoints%numberOfDataPoints,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Data points have not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Data points is not associated.",err,error,*999)
    ENDIF
    
    EXITS("DataPoints_UserNumberGet")
    RETURN
999 ERRORSEXITS("DataPoints_UserNumberGet",err,error)    
    RETURN 1
   
  END SUBROUTINE DataPoints_UserNumberGet
        
  !
  !================================================================================================================================
  !

  !>Changes/sets the user number for a data point identified by a given global number.
  SUBROUTINE DataPoints_UserNumberSet(dataPoints,globalNumber,userNumber,err,error,*)

    !Argument variables
    TYPE(DataPointsType), POINTER :: dataPoints !<A pointer to the data points to set the number for
    INTEGER(INTG), INTENT(IN) :: globalNumber !<The global number to set the user number for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: insertStatus,otherGlobalNumber
    LOGICAL :: dataPointExists
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DataPoints_UserNumberSet",err,error,*999)

    IF(ASSOCIATED(dataPoints)) THEN
      IF(dataPoints%dataPointsFinished) THEN
        CALL FlagError("Data points have been finished.",err,error,*999)
      ELSE
        IF(globalNumber>=1.AND.globalNumber<=dataPoints%numberOfDataPoints) THEN
          !Check the data point user number is not already used
          CALL DataPoint_CheckExists(dataPoints,userNumber,dataPointExists,otherGlobalNumber,err,error,*999)
          IF(dataPointExists) THEN
            IF(otherGlobalNumber/=globalNumber) THEN
              localError="The specified data point user number of "//TRIM(NumberToVString(userNumber,"*",err,error))// &
                & " is already used by global data point number "//TRIM(NumberToVString(otherGlobalNumber,"*",err,error))// &
                & ". User data point numbers must be unique."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ELSE
            CALL Tree_ItemDelete(dataPoints%dataPointsTree,dataPoints%dataPoints(globalNumber)%userNumber,err,error,*999)
            CALL Tree_ItemInsert(dataPoints%dataPointsTree,userNumber,globalNumber,insertStatus,err,error,*999)
            IF(insertStatus/=TREE_NODE_INSERT_SUCESSFUL) CALL FlagError("Unsucessful data points tree insert.",err,error,*999)
            dataPoints%dataPoints(globalNumber)%userNumber=userNumber
          ENDIF
        ELSE
          localError="The specified global data point number of "//TRIM(NumberToVString(globalNumber,"*",err,error))// &
            & " is invalid. The global data point number should be between 1 and "// &
            & TRIM(NumberToVString(dataPoints%numberOfDataPoints,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Data points is not associated.",err,error,*999)
    ENDIF
    
    EXITS("DataPoints_UserNumberSet")
    RETURN
999 ERRORSEXITS("DataPoints_UserNumberSet",err,error)    
    RETURN 1
   
  END SUBROUTINE DataPoints_UserNumberSet
  
  !
  !================================================================================================================================
  !
  
  !>Gets the weights for a data point identified by a given user number.
  SUBROUTINE DataPoints_WeightsGet(dataPoints,userNumber,weights,err,error,*)

    !Argument variables
    TYPE(DataPointsType), POINTER :: dataPoints !<A pointer to the data points to set the number for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number to get the weights for
    REAL(DP), INTENT(OUT) :: weights(:) !<weights(coordinateIdx). On exit, the weights of the specified data point
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: globalNumber
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DataPoints_WeightsGet",err,error,*999)

    IF(ASSOCIATED(dataPoints)) THEN
      IF(dataPoints%dataPointsFinished) THEN
        IF(SIZE(weights,1)>=dataPoints%numberOfDimensions) THEN
          CALL DataPoints_GlobalNumberGet(dataPoints,userNumber,globalNumber,err,error,*999)
          weights(1:dataPoints%numberOfDimensions)=dataPoints%dataPoints(globalNumber)%weights(1:dataPoints%numberOfDimensions)
        ELSE
          localError="The size of the specified weights array of "//TRIM(NumberToVString(SIZE(weights,1),"*",err,error))// &
            & " is too small. The array size needs to be >= "// &
            & TRIM(NumberToVString(dataPoints%numberOfDimensions,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Data points have not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Data points is not associated.",err,error,*999)
    ENDIF
    
    EXITS("DataPoints_WeightsGet")
    RETURN
999 ERRORSEXITS("DataPoints_WeightsGet",err,error)    
    RETURN 1

  END SUBROUTINE DataPoints_WeightsGet

  !
  !================================================================================================================================
  !

  !>Changes/sets the weights for a data point identified by a given user number.
  SUBROUTINE DataPoints_WeightsSet(dataPoints,userNumber,weights,err,error,*)

    !Argument variables
    TYPE(DataPointsType), POINTER :: dataPoints !<A pointer to the data points to set the number for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number to set the weights for
    REAL(DP), INTENT(IN) :: weights(:) !<weights(coordinateIdx). The data point weights to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: globalNumber
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DataPoints_WeightsSet",err,error,*999)

    IF(ASSOCIATED(dataPoints)) THEN
      IF(dataPoints%dataPointsFinished) THEN
        CALL FlagError("Data points have been finished.",err,error,*999)
      ELSE
        IF(SIZE(weights,1)>=dataPoints%numberOfDimensions) THEN
          CALL DataPoints_GlobalNumberGet(dataPoints,userNumber,globalNumber,err,error,*999)
          dataPoints%dataPoints(globalNumber)%weights(1:dataPoints%numberOfDimensions)=weights(1:dataPoints%numberOfDimensions)
        ELSE
          localError="The size of the specified weights array of "//TRIM(NumberToVString(SIZE(weights,1),"*",err,error))// &
            & " is too small. The array size needs to be >= "// &
            & TRIM(NumberToVString(dataPoints%numberOfDimensions,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Data points is not associated.",err,error,*999)
    ENDIF
    
    EXITS("DataPoints_WeightsSet")
    RETURN
999 ERRORSEXITS("DataPoints_WeightsSet",err,error)    
    RETURN 1

  END SUBROUTINE DataPoints_WeightsSet

  !
  !================================================================================================================================
  !  

  !>Gets the data projection identified by a given user number. 
  SUBROUTINE DataPoints_DataProjectionGet(dataPoints,userNumber,dataProjection,err,error,*)

    !Argument variables
    TYPE(DataPointsType), POINTER :: dataPoints !<A pointer to the data points to get the data projection for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number to get the data projection for
    TYPE(DataProjectionType), POINTER :: dataProjection !<On exit, a pointer to the data projection for the data points. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: globalNumber
    TYPE(TREE_NODE_TYPE), POINTER :: treeNode
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DataPoints_DataProjectionGet",err,error,*999)
    
    IF(ASSOCIATED(dataPoints)) THEN
      IF(dataPoints%dataPointsFinished) THEN 
        IF(ASSOCIATED(dataProjection)) THEN
          CALL FlagError("Data projection is already associated.",err,error,*999)
        ELSE
          NULLIFY(treeNode)
          CALL Tree_Search(dataPoints%dataProjectionsTree,userNumber,treeNode,err,error,*999)
          IF(ASSOCIATED(treeNode)) THEN
            CALL Tree_NodeValueGet(dataPoints%dataProjectionsTree,treeNode,globalNumber,err,error,*999)
          ELSE
            localError="The data point projection with user number "//TRIM(NumberToVString(userNumber,"*", &
              & err,error))//" does not exist."
            CALL FlagError(localError,err,error,*999)
          ENDIF
          dataProjection=>dataPoints%dataProjections(globalNumber)%ptr
          IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data points data projection is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Data points has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Data points is not associated.",err,error,*999)
    ENDIF    
    
    EXITS("DataPoints_DataProjectionGet")
    RETURN
999 ERRORSEXITS("DataPoints_DataProjectionGet",err,error)    
    RETURN 1
   
  END SUBROUTINE DataPoints_DataProjectionGet

  !
  !================================================================================================================================
  !  

  !>Gets the user number for a data point identified by a given global number. \todo Is this routine necessary?
  SUBROUTINE DataPoints_DataProjectionGlobalNumberGet(dataPoints,userNumber,globalNumber,err,error,*)

    !Argument variables
    TYPE(DataPointsType), POINTER :: dataPoints !<A pointer to the data points to get the number for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number to get the global number foror
    INTEGER(INTG), INTENT(OUT) :: globalNumber !<On exit, the global number of the specified user data projection
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    TYPE(TREE_NODE_TYPE), POINTER :: treeNode
    
    ENTERS("DataPoints_DataProjectionGlobalNumberGet",err,error,*999)

    IF(ASSOCIATED(dataPoints)) THEN
      IF(dataPoints%dataPointsFinished) THEN
        NULLIFY(treeNode)
        CALL Tree_Search(dataPoints%dataProjectionsTree,userNumber,treeNode,err,error,*999)
        IF(ASSOCIATED(treeNode)) THEN
          CALL Tree_NodeValueGet(dataPoints%dataProjectionsTree,treeNode,globalNumber,err,error,*999)
        ELSE
          localError="The data point projection with user number "//TRIM(NumberToVString(userNumber,"*", &
            & err,error))//" does not exist."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Data points have not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Data points is not associated.",err,error,*999)
    ENDIF
    
    EXITS("DataPoints_DataProjectionGlobalNumberGet")
    RETURN
999 ERRORSEXITS("DataPoints_DataProjectionGlobalNumberGet",err,error)    
    RETURN 1
   
  END SUBROUTINE DataPoints_DataProjectionGlobalNumberGet

  !
  !================================================================================================================================
  !

END MODULE DataPointRoutines



