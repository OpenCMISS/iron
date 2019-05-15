!> \file
!> \author Chris Bradley
!> \brief This module contains all control loop access method routines.
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

!> This module contains all control loop access method routines.
MODULE ControlLoopAccessRoutines
  
  USE BaseRoutines
  USE ISO_VARYING_STRING
  USE Kinds
  USE PROBLEM_CONSTANTS
  USE Strings
  USE Types
  
#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup OpenCMISS_ControlLoopConstants OpenCMISS::Iron::ControlLoop::Constants
  !> \brief Control loop constants.
  !>@{
  !> \addtogroup ControlLoopRoutines_ControlLoopIdentifiers ControlLoopRoutines::ControlLoopIdentifiers
  !> \brief The control loop identification parameters
  !>@{
  INTEGER(INTG), PARAMETER :: CONTROL_LOOP_NODE=0 !<The identifier for a each "leaf" node in a control loop. \see ControlLoopRoutines_ControlLoopIdentifiers,ControlLoopRoutines
  !>@}
  !> \addtogroup ControlLoopRoutines_ControlLoopTypes ControlLoopRoutines::ControlLoopTypes
  !> \brief Control loop type parameters
  !> \see ControlLoopRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: CONTROL_SIMPLE_TYPE=1 !<Simple, one iteration control loop. \see ControlLoopRoutines_ControlLoopTypes,ControlLoopRoutines
  INTEGER(INTG), PARAMETER :: CONTROL_FIXED_LOOP_TYPE=2 !<Fixed iteration control loop. \see ControlLoopRoutines_ControlLoopTypes,ControlLoopRoutines
  INTEGER(INTG), PARAMETER :: CONTROL_TIME_LOOP_TYPE=3 !<Time control loop. \see ControlLoopRoutines_ControlLoopTypes,ControlLoopRoutines
  INTEGER(INTG), PARAMETER :: CONTROL_WHILE_LOOP_TYPE=4 !<While control loop. \see ControlLoopRoutines_ControlLoopTypes,ControlLoopRoutines
  INTEGER(INTG), PARAMETER :: CONTROL_LOAD_INCREMENT_LOOP_TYPE=5 !<Load increment control loop. \see ControlLoopRoutines_ControlLoopTypes,ControlLoopRoutines
  !>@}
  !> \addtogroup ControlLoop_OutputTypes OpenCMISS::Iron::ControlLoop::OutputTypes
  !> \brief The types of output for a control loop.
  !> \see ControlLoop
  !>@{
  INTEGER(INTG), PARAMETER :: CONTROL_LOOP_NO_OUTPUT=0 !<No output from the control loop \see ControlLoop_OutputTypes,ControlLoop
  INTEGER(INTG), PARAMETER :: CONTROL_LOOP_PROGRESS_OUTPUT=1 !<Progress output from control loop \see ControlLoop_OutputTypes,ControlLoop
  INTEGER(INTG), PARAMETER :: CONTROL_LOOP_TIMING_OUTPUT=2 !<Timing output from the control loop \see ControlLoop_OutputTypes,ControlLoop
  !>@}

  !> \addtogroup ControlLoop_FieldVariableLinearityTypes OpenCMISS::Iron::ControlLoop::FieldVariableLinearityTypes
  !> \brief The linearity type of control loop field variables
  !> \see ControlLoop
  !>@{
  INTEGER(INTG), PARAMETER :: CONTROL_LOOP_FIELD_VARIABLE_LINEAR=1 !<The control loop field variable is linear \see ControlLoop_FieldVariableLinearityTypes,ControlLoop
  INTEGER(INTG), PARAMETER :: CONTROL_LOOP_FIELD_VARIABLE_NONLINEAR=2 !<The control loop field variable is nonlinear \see ControlLoop_FieldVariableLinearityTypes,ControlLoop
  !>@}

  !> \addtogroup ControlLoop_FieldVariableTimeDependenceTypes OpenCMISS::Iron::ControlLoop::FieldVariableTimeDependenceTypes
  !> \brief The time dependence type of control loop field variables
  !> \see ControlLoop
  !>@{
  INTEGER(INTG), PARAMETER :: CONTROL_LOOP_FIELD_VARIABLE_STATIC=1 !<The control loop field variable is static \see ControlLoop_FieldVariableTimeDependenceTypes,ControlLoop
  INTEGER(INTG), PARAMETER :: CONTROL_LOOP_FIELD_VARIABLE_QUASISTATIC=2 !<The control loop field variable is quasistatic \see ControlLoop_FieldVariableTimeDependenceTypes,ControlLoop
  INTEGER(INTG), PARAMETER :: CONTROL_LOOP_FIELD_VARIABLE_FIRST_DEGREE_DYNAMIC=3 !<The control loop field variable is first degree dynamic i.e., we use first time derivatives \see ControlLoop_FieldVariableTimeDependenceTypes,ControlLoop
  INTEGER(INTG), PARAMETER :: CONTROL_LOOP_FIELD_VARIABLE_SECOND_DEGREE_DYNAMIC=4 !<The control loop field variable is second degree dynamic i.e., we use second time derivatives \see ControlLoop_FieldVariableTimeDependenceTypes,ControlLoop
  !>@}
  !>@}
  
  !Module types

  !Module variables
  
  !Interfaces

  !>Returns the specified control loop as indexed by the control loop identifier from the control loop root. \see OpenCMISS::Iron::cmge_ControlLoop_Get
  INTERFACE CONTROL_LOOP_GET
    MODULE PROCEDURE ControlLoop_Get0
    MODULE PROCEDURE ControlLoop_Get1
  END INTERFACE CONTROL_LOOP_GET

  !>Returns the specified control loop as indexed by the control loop identifier from the control loop root. \see OpenCMISS::Iron::cmfe_ControlLoop_Get
  INTERFACE ControlLoop_Get
    MODULE PROCEDURE ControlLoop_Get0
    MODULE PROCEDURE ControlLoop_Get1
  END INTERFACE ControlLoop_Get
  
  INTERFACE CONTROL_LOOP_CURRENT_TIMES_GET
    MODULE PROCEDURE ControlLoop_CurrentTimesGet
  END INTERFACE CONTROL_LOOP_CURRENT_TIMES_GET

  INTERFACE CONTROL_LOOP_SOLVERS_GET
    MODULE PROCEDURE ControlLoop_SolversGet
  END INTERFACE CONTROL_LOOP_SOLVERS_GET

  PUBLIC CONTROL_LOOP_NODE

  PUBLIC CONTROL_SIMPLE_TYPE,CONTROL_FIXED_LOOP_TYPE,CONTROL_TIME_LOOP_TYPE,CONTROL_WHILE_LOOP_TYPE, &
    & CONTROL_LOAD_INCREMENT_LOOP_TYPE

  PUBLIC CONTROL_LOOP_NO_OUTPUT,CONTROL_LOOP_PROGRESS_OUTPUT,CONTROL_LOOP_TIMING_OUTPUT
  
  PUBLIC CONTROL_LOOP_FIELD_VARIABLE_LINEAR,CONTROL_LOOP_FIELD_VARIABLE_NONLINEAR

  PUBLIC CONTROL_LOOP_FIELD_VARIABLE_STATIC,CONTROL_LOOP_FIELD_VARIABLE_QUASISTATIC, &
    & CONTROL_LOOP_FIELD_VARIABLE_FIRST_DEGREE_DYNAMIC,CONTROL_LOOP_FIELD_VARIABLE_SECOND_DEGREE_DYNAMIC

  PUBLIC ControlLoop_AssertIsFinished,ControlLoop_AssertNotFinished

  PUBLIC ControlLoop_AssertIsFixedLoop

  PUBLIC ControlLoop_AssertIsLoadIncrementLoop

  PUBLIC ControlLoop_AssertIsSimpleLoop

  PUBLIC ControlLoop_AssertIsTimeLoop

  PUBLIC ControlLoop_AssertIsWhileLoop

  PUBLIC ControlLoop_Get

  PUBLIC CONTROL_LOOP_GET

  PUBLIC CONTROL_LOOP_CURRENT_TIMES_GET

  PUBLIC ControlLoop_CurrentTimesGet

  PUBLIC ControlLoop_CurrentTimeInformationGet

  PUBLIC ControlLoop_FixedLoopGet

  PUBLIC ControlLoop_LoadIncrementLoopGet

  PUBLIC ControlLoop_NumberOfIterationsGet

  PUBLIC ControlLoop_NumberOfSubLoopsGet

  PUBLIC ControlLoop_OutputTypeGet

  PUBLIC ControlLoop_ProblemGet

  PUBLIC ControlLoop_SimpleLoopGet

  PUBLIC ControlLoop_SolversGet

  PUBLIC CONTROL_LOOP_SOLVERS_GET

  PUBLIC ControlLoop_SubLoopGet

  PUBLIC ControlLoop_TimeLoopGet

  PUBLIC ControlLoop_WhileLoopGet

CONTAINS

  !
  !================================================================================================================================
  !

  !>Assert that a control loop has been finished
  SUBROUTINE ControlLoop_AssertIsFinished(controlLoop,err,error,*)

    !Argument Variables
    TYPE(ControlLoopType), POINTER, INTENT(IN) :: controlLoop !<The control loop to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("ControlLoop_AssertIsFinished",err,error,*999)

    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("Control loop is not associated.",err,error,*999)

    IF(.NOT.controlLoop%controlLoopFinished) CALL FlagError("Control loop has not been finished.",err,error,*999)
    
    EXITS("ControlLoop_AssertIsFinished")
    RETURN
999 ERRORSEXITS("ControlLoop_AssertIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_AssertIsFinished

  !
  !================================================================================================================================
  !

  !>Assert that a control loop has not been finished
  SUBROUTINE ControlLoop_AssertNotFinished(controlLoop,err,error,*)

    !Argument Variables
    TYPE(ControlLoopType), POINTER, INTENT(IN) :: controlLoop !<The control loop to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("ControlLoop_AssertNotFinished",err,error,*999)

    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("Control loop is not associated.",err,error,*999)

    IF(controlLoop%controlLoopFinished) CALL FlagError("Control loop has already been finished.",err,error,*999)
    
    EXITS("ControlLoop_AssertNotFinished")
    RETURN
999 ERRORSEXITS("ControlLoop_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_AssertNotFinished

  !
  !================================================================================================================================
  !

  !>Assert that a control loop is a fixed loop
  SUBROUTINE ControlLoop_AssertIsFixedLoop(controlLoop,err,error,*)

    !Argument Variables
    TYPE(ControlLoopType), POINTER, INTENT(IN) :: controlLoop !<The control loop to assert the fixed loop for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("ControlLoop_AssertIsFixedLoop",err,error,*999)

    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("Control loop is not associated.",err,error,*999)

    IF(controlLoop%loopType/=CONTROL_FIXED_LOOP_TYPE) &
      & CALL FlagError("The specified control loop is not a fixed control loop.",err,error,*999)
          
    EXITS("ControlLoop_AssertIsFixedLoop")
    RETURN
999 ERRORSEXITS("ControlLoop_AssertIsFixedLoop",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_AssertIsFixedLoop

  !
  !================================================================================================================================
  !

  !>Assert that a control loop is a load increment loop
  SUBROUTINE ControlLoop_AssertIsLoadIncrementLoop(controlLoop,err,error,*)

    !Argument Variables
    TYPE(ControlLoopType), POINTER, INTENT(IN) :: controlLoop !<The control loop to assert the load increment loop for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("ControlLoop_AssertIsLoadIncrementLoop",err,error,*999)

    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("Control loop is not associated.",err,error,*999)

    IF(controlLoop%loopType/=CONTROL_LOAD_INCREMENT_LOOP_TYPE) &
      & CALL FlagError("The specified control loop is not a load increment control loop.",err,error,*999)
          
    EXITS("ControlLoop_AssertIsLoadIncrementLoop")
    RETURN
999 ERRORSEXITS("ControlLoop_AssertIsLoadIncrementLoop",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_AssertIsLoadIncrementLoop

  !
  !================================================================================================================================
  !

  !>Assert that a control loop is a simple loop
  SUBROUTINE ControlLoop_AssertIsSimpleLoop(controlLoop,err,error,*)

    !Argument Variables
    TYPE(ControlLoopType), POINTER, INTENT(IN) :: controlLoop !<The control loop to assert the simple loop for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("ControlLoop_AssertIsSimpleLoop",err,error,*999)

    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("Control loop is not associated.",err,error,*999)

    IF(controlLoop%loopType/=CONTROL_SIMPLE_TYPE) &
      & CALL FlagError("The specified control loop is not a simple control loop.",err,error,*999)
          
    EXITS("ControlLoop_AssertIsSimpleLoop")
    RETURN
999 ERRORSEXITS("ControlLoop_AssertIsSimpleLoop",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_AssertIsSimpleLoop

  !
  !================================================================================================================================
  !

  !>Assert that a control loop is a time loop
  SUBROUTINE ControlLoop_AssertIsTimeLoop(controlLoop,err,error,*)

    !Argument Variables
    TYPE(ControlLoopType), POINTER, INTENT(IN) :: controlLoop !<The control loop to assert the time loop for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("ControlLoop_AssertIsTimeLoop",err,error,*999)

    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("Control loop is not associated.",err,error,*999)

    IF(controlLoop%loopType/=CONTROL_TIME_LOOP_TYPE) &
      & CALL FlagError("The specified control loop is not a time control loop.",err,error,*999)
          
    EXITS("ControlLoop_AssertIsTimeLoop")
    RETURN
999 ERRORSEXITS("ControlLoop_AssertIsTimeLoop",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_AssertIsTimeLoop

  !
  !================================================================================================================================
  !

  !>Assert that a control loop is a while loop
  SUBROUTINE ControlLoop_AssertIsWhileLoop(controlLoop,err,error,*)

    !Argument Variables
    TYPE(ControlLoopType), POINTER, INTENT(IN) :: controlLoop !<The control loop to assert the while loop for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("ControlLoop_AssertIsWhileLoop",err,error,*999)

    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("Control loop is not associated.",err,error,*999)

    IF(controlLoop%loopType/=CONTROL_WHILE_LOOP_TYPE) &
      & CALL FlagError("The specified control loop is not a while control loop.",err,error,*999)
          
    EXITS("ControlLoop_AssertIsWhileLoop")
    RETURN
999 ERRORSEXITS("ControlLoop_AssertIsWhileLoop",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_AssertIsWhileLoop

  !
  !================================================================================================================================
  !

  !>Returns the specified control loop as indexed by the control loop identifier from the control loop root.
  SUBROUTINE ControlLoop_Get0(controlLoopRoot,controlLoopIdentifier,controlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER, INTENT(IN) :: controlLoopRoot!<A pointer to the control loop to root
    INTEGER(INTG), INTENT(IN) :: controlLoopIdentifier !<The control loop identifier
    TYPE(ControlLoopType), POINTER :: controlLoop !<On exit, the specified control loop
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("ControlLoop_Get0",err,error,*999)

    CALL ControlLoop_Get1(controlLoopRoot,[controlLoopIdentifier],controlLoop,err,error,*999)
       
    EXITS("ControlLoop_Get0")
    RETURN
999 ERRORSEXITS("ControlLoop_Get0",err,error)
    RETURN 1
  END SUBROUTINE ControlLoop_Get0

  !
  !================================================================================================================================
  !

  !>Returns the specified control loop as indexed by the control loop identifier from the control loop root.
  SUBROUTINE ControlLoop_Get1(controlLoopRoot,controlLoopIdentifiers,controlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER, INTENT(IN) :: controlLoopRoot !<A pointer to the control loop to root
    INTEGER(INTG), INTENT(IN) :: controlLoopIdentifiers(:) !<controlLoopIdentifiers(identifierIdx). The control loop identifiers
    TYPE(ControlLoopType), POINTER :: controlLoop !<On exit, the specified control loop. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: controlLoopIdx
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("ControlLoop_Get1",err,error,*998)

    IF(.NOT.ASSOCIATED(controlLoopRoot)) CALL FlagError("Control loop root is not associated.",err,error,*998)
    IF(ASSOCIATED(controlLoop)) CALL FlagError("Control loop is already associated.",err,error,*998)      
    IF(.NOT.COUNT(controlLoopIdentifiers==CONTROL_LOOP_NODE)==1) THEN
      localError="Invalid control loop identifier. The control loop identifier has "// &
        & TRIM(NumberToVString(COUNT(controlLoopIdentifiers==CONTROL_LOOP_NODE),"*",err,error))// &
        & " control loop node identifiers and it should only have 1."
      CALL FlagError(localError,err,error,*998)
    ENDIF
    IF(.NOT.controlLoopIdentifiers(SIZE(controlLoopIdentifiers,1))==CONTROL_LOOP_NODE) THEN
      localError="Invalid control loop identifier. The last value in the identifier vector is "// &
        & TRIM(NumberToVString(controlLoopIdentifiers(SIZE(controlLoopIdentifiers,1)),"*",err,error))// &
        & " and it should be "//TRIM(NumberToVString(CONTROL_LOOP_NODE,"*",err,error))//"."
      CALL FlagError(localError,err,error,*998)
    ENDIF
    
    controlLoop=>controlLoopRoot
    DO controlLoopIdx=1,SIZE(controlLoopIdentifiers,1)
      IF(controlLoopIdentifiers(controlLoopIdx)==CONTROL_LOOP_NODE) THEN
        EXIT
      ELSE
        IF(controlLoopIdentifiers(controlLoopIdx)>0.AND. &
          & controlLoopIdentifiers(controlLoopIdx)<=controlLoop%numberOfSubLoops) THEN
          controlLoop=>controlLoop%subLoops(controlLoopIdentifiers(controlLoopIdx))%ptr
          IF(.NOT.ASSOCIATED(controlLoop)) THEN
            localError="Control sub loop number "//TRIM(NumberToVString(controlLoopIdentifiers(controlLoopIdx),"*",err,error))// &
              & " at identifier index "//TRIM(NumberToVString(controlLoopIdx,"*",err,error))//" is not associated."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ELSE
          localError="Invalid control loop identifier. The identifier at index "// &
            & TRIM(NumberToVString(controlLoopIdx,"*",err,error))//" is "// &
            & TRIM(NumberToVString(controlLoopIdentifiers(controlLoopIdx),"*",err,error))// &
            & ". The identifier must be between 1 and "//TRIM(NumberToVString(controlLoop%numberOfSubLoops,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDIF
    ENDDO !controlLoopIdx
        
    EXITS("ControlLoop_Get1")
    RETURN
999 NULLIFY(controlLoop)
998 ERRORSEXITS("ControlLoop_Get1",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_Get1

  !
  !================================================================================================================================
  !

  !>Gets the current time parameters for a time control loop. If the specified loop is not a time loop the next time loop up the chain will be used. \see OpenCMISS::cmfe_ControlLoop_CurrentTimesGet
  SUBROUTINE ControlLoop_CurrentTimesGet(controlLoop,currentTime,timeIncrement,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER, INTENT(IN) :: controlLoop !<The control loop to get the current times for
    REAL(DP), INTENT(OUT) :: currentTime !<On exit, the current time.
    REAL(DP), INTENT(OUT) :: timeIncrement !<On exit, the current time increment.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: currentIteration,outputIteration
    REAL(DP) :: startTime,stopTime

    ENTERS("ControlLoop_CurrentTimesGet",err,error,*999)

    CALL ControlLoop_CurrentTimeInformationGet(controlLoop,currentTime,timeIncrement,startTime,stopTime,currentIteration, &
      & outputIteration,err,error,*999)
       
    EXITS("ControlLoop_CurrentTimesGet")
    RETURN
999 ERRORSEXITS("ControlLoop_CurrentTimesGet",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_CurrentTimesGet
  
  !
  !================================================================================================================================
  !

  !>Gets the current loop information for a time control loop. If the specified loop is not a time loop the next time loop up the chain will be used.
  SUBROUTINE ControlLoop_CurrentTimeInformationGet(controlLoop,currentTime,timeIncrement,startTime,stopTime,currentIteration, &
    & outputIteration,err,error,*)
    
    !Argument variables
    TYPE(ControlLoopType), POINTER, INTENT(IN) :: controlLoop !<The control loop to get the time information for
    REAL(DP), INTENT(OUT) :: currentTime !<On exit, the current time.
    REAL(DP), INTENT(OUT) :: timeIncrement !<On exit, the current time increment.
    REAL(DP), INTENT(OUT) :: startTime !<On exit, the start time for the loop
    REAL(DP), INTENT(OUT) :: stopTime !<On exit, the stop time for the loop
    INTEGER(INTG), INTENT(OUT) :: currentIteration !<On exit, the current iteration number for the loop
    INTEGER(INTG), INTENT(OUT) :: outputIteration !<On exit, the output iteration number for the loop
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables    
    INTEGER(INTG) :: controlLoopLevel,levelIdx
    TYPE(ControlLoopType), POINTER :: parentLoop
    TYPE(ControlLoopTimeType), POINTER :: timeLoop

    ENTERS("ControlLoop_CurrentTimeInformationGet",err,error,*999)

    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("Control loop is not associated.",err,error,*999)
    IF(.NOT.controlLoop%controlLoopFinished) CALL FlagError("Control loop has not been finished.",err,error,*999)

    !Find a time loop from either the specified control loop or the next time loop up the chain.
    controlLoopLevel=controlLoop%controlLoopLevel
    parentLoop=>controlLoop
    DO levelIdx=controlLoopLevel,1,-1
      IF(controlLoopLevel==0) THEN
        CALL FlagError("Could not find a time loop for the specified control loop.",err,error,*999)
      ELSE
        IF(parentLoop%loopType==CONTROL_TIME_LOOP_TYPE) THEN
          timeLoop=>parentLoop%timeLoop
          IF(ASSOCIATED(timeLoop)) THEN
            currentTime=timeLoop%currentTime
            timeIncrement=timeLoop%timeIncrement
            startTime=timeLoop%startTime
            stopTime=timeLoop%stopTime
            currentIteration=timeLoop%iterationNumber
            outputIteration=timeLoop%outputNumber
          ELSE
            CALL FlagError("Control loop time loop is not associated.",err,error,*999)
          ENDIF
          EXIT
        ELSE
          parentLoop=>parentLoop%parentLoop
        ENDIF
      ENDIF
    ENDDO !levelIdx
       
    EXITS("ControlLoop_CurrentTimeInformationGet")
    RETURN
999 ERRORSEXITS("ControlLoop_CurrentTimeInformationGet",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_CurrentTimeInformationGet
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the fixed loop for a control loop.
  SUBROUTINE ControlLoop_FixedLoopGet(controlLoop,fixedLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER, INTENT(IN) :: controlLoop !<A pointer to control loop to get the problem for.
    TYPE(ControlLoopFixedType), POINTER :: fixedLoop !<On exit, a pointer to the control loop fixed loop. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("ControlLoop_FixedLoopGet",err,error,*998)

    IF(ASSOCIATED(fixedLoop)) CALL FlagError("Fixed loop is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("Control loop is not associated.",err,error,*999)

    fixedLoop=>controlLoop%fixedLoop
    IF(.NOT.ASSOCIATED(fixedLoop)) CALL FlagError("Control loop fixed loop is not associated.",err,error,*999)
       
    EXITS("ControlLoop_FixedLoopGet")
    RETURN
999 NULLIFY(fixedLoop)
998 ERRORSEXITS("ControlLoop_FixedLoopGet",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_FixedLoopGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the load increment loop for a control loop.
  SUBROUTINE ControlLoop_LoadIncrementLoopGet(controlLoop,loadIncrementLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER, INTENT(IN) :: controlLoop !<A pointer to control loop to get the problem for.
    TYPE(ControlLoopLoadIncrementType), POINTER :: loadIncrementLoop !<On exit, a pointer to the control loop load increment loop. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("ControlLoop_LoadIncrementLoopGet",err,error,*998)

    IF(ASSOCIATED(loadIncrementLoop)) CALL FlagError("Load increment loop is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("Control loop is not associated.",err,error,*999)

    loadIncrementLoop=>controlLoop%loadIncrementLoop
    IF(.NOT.ASSOCIATED(loadIncrementLoop)) CALL FlagError("Control loop load increment loop is not associated.",err,error,*999)
       
    EXITS("ControlLoop_LoadIncrementLoopGet")
    RETURN
999 NULLIFY(loadIncrementLoop)
998 ERRORSEXITS("ControlLoop_LoadIncrementLoopGet",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_LoadIncrementLoopGet

  !
  !================================================================================================================================
  !

  !>Gets the number of iterations for a time type control loop. If the value is not set to something /=0, it will be computed the first time the loop is executed. If it is retrieved earlier and 0 is returned, this means the value was not yet computed. \see OpenCMISS_ControlLoop_NumberOfIterationsGet
  SUBROUTINE ControlLoop_NumberOfIterationsGet(controlLoop,numberOfIterations,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER, INTENT(IN) :: controlLoop !<A pointer to time control loop to set the number of iterations for
    INTEGER(INTG), INTENT(OUT) :: numberOfIterations !<The number of iterations for the time control loop.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ControlLoopTimeType), POINTER :: timeLoop
 
    ENTERS("ControlLoop_NumberOfIterationsGet",err,error,*999)

    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("Control loop is not associated.",err,error,*999)
    CALL ControlLoop_AssertIsTimeLoop(controlLoop,err,error,*999)
    NULLIFY(timeLoop)
    CALL ControlLoop_TimeLoopGet(controlLoop,timeLoop,err,error,*999)
    
    numberOfIterations=timeLoop%numberOfIterations
       
    EXITS("ControlLoop_NumberOfIterationsGet")
    RETURN
999 ERRORSEXITS("ControlLoop_NumberOfIterationsGet",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_NumberOfIterationsGet
  
  !
  !================================================================================================================================
  !

  !>Gets the number of sub loops for a control loop. \see OpenCMISS::Iron::cmfe_ControlLoop_NumberOfSubLoopsGet
  SUBROUTINE ControlLoop_NumberOfSubLoopsGet(controlLoop,numberOfSubLoops,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER, INTENT(IN) :: controlLoop !<A pointer to control loop to get the number of sub loops for
    INTEGER(INTG), INTENT(OUT) :: numberOfSubLoops !<On return, the number of sub loops for the specified control loop
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("ControlLoop_NumberOfSubLoopsGet",err,error,*999)

    CALL ControlLoop_AssertIsFinished(controlLoop,err,error,*999)

    numberOfSubLoops=controlLoop%numberOfSubLoops
        
    EXITS("ControlLoop_NumberOfSubLoopsGet")
    RETURN
999 ERRORSEXITS("ControlLoop_NumberOfSubLoopsGet",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_NumberOfSubLoopsGet

  !
  !================================================================================================================================
  !

  !>Gets the output type for a control loop. \see OpenCMISS::Iron::cmfe_ControlLoop_OutputTypeGet
  SUBROUTINE ControlLoop_OutputTypeGet(controlLoop,outputType,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER, INTENT(IN) :: controlLoop !<A pointer to the control loop to get the output type for.
    INTEGER(INTG), INTENT(OUT) :: outputType !<On exit, the output type of the control loop \see ControlLoopRoutines_OutputTypes,ControlLoopRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("ControlLoop_OutputTypeGet",err,error,*999)

    CALL ControlLoop_AssertIsFinished(controlLoop,err,error,*999)

    outputType=controlLoop%outputType
       
    EXITS("ControlLoop_OutputTypeGet")
    RETURN
999 ERRORSEXITS("ControlLoop_OutputTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_OutputTypeGet
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the problem for a control loop.
  SUBROUTINE ControlLoop_ProblemGet(controlLoop,problem,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER, INTENT(IN) :: controlLoop !<A pointer to control loop to get the problem for.
    TYPE(ProblemType), POINTER :: problem !<On exit, a pointer to the control loop problem. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("ControlLoop_ProblemGet",err,error,*998)

    IF(ASSOCIATED(problem)) CALL FlagError("Problem is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("Control loop is not associated.",err,error,*999)

    problem=>controlLoop%problem
    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Control loop problem is not associated.",err,error,*999)
       
    EXITS("ControlLoop_ProblemGet")
    RETURN
998 NULLIFY(problem)
999 ERRORSEXITS("ControlLoop_ProblemGet",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_ProblemGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the simple loop for a control loop.
  SUBROUTINE ControlLoop_SimpleLoopGet(controlLoop,simpleLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER, INTENT(IN) :: controlLoop !<A pointer to control loop to get the problem for.
    TYPE(ControlLoopSimpleType), POINTER :: simpleLoop !<On exit, a pointer to the control loop simple loop. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("ControlLoop_SimpleLoopGet",err,error,*998)

    IF(ASSOCIATED(simpleLoop)) CALL FlagError("Simple loop is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("Control loop is not associated.",err,error,*999)

    simpleLoop=>controlLoop%simpleLoop
    IF(.NOT.ASSOCIATED(simpleLoop)) CALL FlagError("Control loop simple loop is not associated.",err,error,*999)
       
    EXITS("ControlLoop_SimpleLoopGet")
    RETURN
999 NULLIFY(simpleLoop)
998 ERRORSEXITS("ControlLoop_SimpleLoopGet",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_SimpleLoopGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the solvers for a control loop.
  SUBROUTINE ControlLoop_SolversGet(controlLoop,solvers,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER, INTENT(IN) :: controlLoop !<A pointer to control loop to get the solvers for.
    TYPE(SOLVERS_TYPE), POINTER :: solvers !<On exit, a pointer to the control loop solvers. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("ControlLoop_SolversGet",err,error,*998)

    IF(ASSOCIATED(solvers)) CALL FlagError("Solvers is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("Control loop is not associated.",err,error,*999)
 
    solvers=>controlLoop%solvers
    IF(.NOT.ASSOCIATED(solvers)) CALL FlagError("Control loop solvers is not associated.",err,error,*999)
       
    EXITS("ControlLoop_SolversGet")
    RETURN
998 NULLIFY(solvers)
999 ERRORSEXITS("ControlLoop_SolversGet",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_SolversGet

  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the specified sub loop for a control loop.
  SUBROUTINE ControlLoop_SubLoopGet(controlLoop,subLoopIdx,subLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the control loop to get the sub loop for
    INTEGER(INTG), INTENT(IN) :: subLoopIdx !<The sub loop index in the control loop to get the sub loop for
    TYPE(ControlLoopType), POINTER :: subLoop !<On exit, a pointer to the specified sub loop. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("ControlLoop_SubLoopGet",err,error,*998)

    IF(ASSOCIATED(subLoop)) CALL FlagError("Sub loop is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("Control loop is not associated.",err,error,*999)
    IF(subLoopIdx<=0.OR.subLoopIdx>controlLoop%numberOfSubLoops) THEN
      localError="The specified sub loop index of "//TRIM(NumberToVString(subLoopIdx,"*",err,error))// &
        & " is invalid. The index must be > 0 and <= "// TRIM(NumberToVString(controlLoop%numberOfSubLoops,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(controlLoop%subLoops)) &
      & CALL FlagError("Control loop sub loops is not allocated.",err,error,*999)

    subLoop=>controlLoop%subLoops(subLoopIdx)%ptr
    IF(.NOT.ASSOCIATED(subLoop)) THEN
      localError="The sub loop for the specified sub loop index of "// &
        & TRIM(NumberToVString(subLoopIdx,"*",err,error))//" is not associated."      
      CALL FlagError(localError,err,error,*999)
    ENDIF
      
    EXITS("ControlLoop_SubLoopGet")
    RETURN
999 NULLIFY(subLoop)
998 ERRORSEXITS("ControlLoop_SubLoopGet",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_SubLoopGet
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the time loop for a control loop.
  SUBROUTINE ControlLoop_TimeLoopGet(controlLoop,timeLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER, INTENT(IN) :: controlLoop !<A pointer to control loop to get the problem for.
    TYPE(ControlLoopTimeType), POINTER :: timeLoop !<On exit, a pointer to the control loop time loop. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("ControlLoop_TimeLoopGet",err,error,*998)

    IF(ASSOCIATED(timeLoop)) CALL FlagError("Time loop is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("Control loop is not associated.",err,error,*999)

    timeLoop=>controlLoop%timeLoop
    IF(.NOT.ASSOCIATED(timeLoop)) CALL FlagError("Control loop time loop is not associated.",err,error,*999)
       
    EXITS("ControlLoop_TimeLoopGet")
    RETURN
999 NULLIFY(timeLoop)
998 ERRORSEXITS("ControlLoop_TimeLoopGet",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_TimeLoopGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the while loop for a control loop.
  SUBROUTINE ControlLoop_WhileLoopGet(controlLoop,whileLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER, INTENT(IN) :: controlLoop !<A pointer to control loop to get the problem for.
    TYPE(ControlLoopWhileType), POINTER :: whileLoop !<On exit, a pointer to the control loop while loop. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("ControlLoop_WhileLoopGet",err,error,*998)

    IF(ASSOCIATED(whileLoop)) CALL FlagError("While loop is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("Control loop is not associated.",err,error,*999)

    whileLoop=>controlLoop%whileLoop
    IF(.NOT.ASSOCIATED(whileLoop)) CALL FlagError("Control loop while loop is not associated.",err,error,*999)
       
    EXITS("ControlLoop_WhileLoopGet")
    RETURN
999 NULLIFY(whileLoop)
998 ERRORSEXITS("ControlLoop_WhileLoopGet",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_WhileLoopGet

  !
  !================================================================================================================================
  !

END MODULE ControlLoopAccessRoutines
